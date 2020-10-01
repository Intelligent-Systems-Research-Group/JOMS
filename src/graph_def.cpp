
#include "graph_def.h"

extern "C" {
#include "Opt.h"
}

#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "types.h"
#include "MeshTopology.h"
#include "corpora.h"
#include <random>
#include "skeleton.h"
#include "VertexGroupLoader.h"
#include "rotations.h"
#include "model.h"
#include "exporter.h"
#include <time.h>
#include <fstream>
#include <memory>
#include "SubdivEvaluator.h"
#include "correspondences.h"
#include "calculus.h"
#include <sys/stat.h>
#include <sys/types.h>
//#include "InputData.h"
#include <chrono>
//#include <gperftools/profiler.h>
#include "Camera.h"
#include "Pose2d.h"
#include "Icp.h"
#include "InputData.h"
#include "TermBuilder.h"
#include "cost.h"


Vector3 s2v(Scalar3 u) {
    Vector3 v;
    v[0] = u.x;
    v[1] = u.y;
    v[2] = u.z;
    return v;
}

const  Eigen::IOFormat XYZFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n");



void dumpCloudA(std::string name, int iter, const Matrix3X& cloud, const Matrix3X* normals = NULL) {
    std::ofstream file(name + std::string(".xyz"));
    if (normals == NULL) {
        file << cloud.transpose().format(XYZFormat);
    }
    else {
        MatrixX M(6, cloud.cols());
        M.block(0, 0, 3, cloud.cols()) = cloud;
        M.block(3, 0, 3, cloud.cols()) = *normals;
        file << M.transpose().format(XYZFormat);
    }
}

void dumpSubivA(std::string name, int iter, std::unique_ptr<SubdivEvaluator>& eval, const Matrix3X& vertices) {
    MeshTopology resultTopology;
    Matrix3X resultVertices;
    eval->generate_refined_mesh(vertices, 2, &resultTopology, &resultVertices);
    auto path = name + std::string("_")  + std::to_string(iter) + std::string(".obj");
    saveObj(path, &resultTopology, &resultVertices);
}


void unravel(const Matrix3X& joints, std::vector<Vector3>& locs) {

    locs.resize(joints.cols());
    for(int i = 0; i < joints.cols(); i++) {
        locs[i] = joints.col(i);
    }
}


Problem::Problem(const InputParams& ip, const Corpora& _corpora) : inputParams(ip),
                                                              vmarkers(ip.nscans),
                                                              vreal(ip.nscans),
                                                              cloud(ip.nscans),
                                                              model2Data(ip.nscans),
                                                              cloudNorm(ip.nscans),
                                                              topologies(ip.nscans),
                                                              corpora(_corpora),
															  skeleton(ip.template_skeleton_path),
                                                              REST_SIZE( (1+ip.ncomp+(ip.ncomp+1)+ip.ndeformjoints*(ip.ndeformshapes+2 + 1) + 1) ),
                                                              MODEL_SIZE (23),
                                                              PJOINT_SIZE (1 +1+2*ip.ncomp),
                                                              SCALE(1.0f),
                                                              LAMBDA(100)

    {
        std::cerr << "Problem Created" << std::endl;
        state = NULL;
        problem = NULL;
        plan = NULL;

        Opt_InitializationParameters param = {};
#ifndef USE_DOUBLE
        param.doublePrecision = 0;
#else
        param.doublePrecision = 1;
#endif
        param.verbosityLevel = 0;
        param.collectPerKernelTimingInfo = 0;
        //param.threadsPerBlock = 512;

        std::cout << "doublePrecision " << param.doublePrecision << std::endl;
        std::cout << "verbosityLevel " << param.verbosityLevel << std::endl;
        std::cout << "collectPerKernelTimingInfo " << param.collectPerKernelTimingInfo << std::endl;
        std::cout << "threadsPerBlock " << param.threadsPerBlock << std::endl;

        std::cout << "Start NewState" << std::endl;
        state = Opt_NewState(param);
        problem = Opt_ProblemDefine(state, ip.optscript_path.c_str(), "LMGPU"); //"gaussNewtonGPU"


#ifdef PRECOMPUTE_SIZE
        initOpt();
#endif

        //exit(0);

        //corpora.scanCount.clear();
        //for(int i = 0; i < ip.npersons; i++) {
        //    corpora.scanCount.push_back(10);
        //}
        w_fitSqrt = 0; //.5;
        w_surface = 1;
        termIdx = 0;
        reverse_offset = -1;
    }

    Matrix3X Problem::calcJointInit() {
        auto jointMap = calculateJoints(border_groups, &control_vertices);
        Matrix3X jointPriorCloud(3,inputParams.njoints);
        jointPriorCloud.setConstant(0);
        for (size_t i = 0; i < inputParams.njoints-1; i++) {
            auto key = skeleton.getJoint(i + 1);
            auto secondKey = key;
            if (border_groups.find(key) == border_groups.end()) {
                secondKey.first = key.second;
                secondKey.second = key.first;
            }
            Vector3 lu = jointMap.at(secondKey);
            jointPriorCloud.col(i+1) = lu;
        }
        return jointPriorCloud;
    }

    MatrixX Problem::calcSkinningWeights(int njoints, int nverts) {
        MatrixX ws; ws.resize(njoints,nverts);
        for(int i = 0; i < nverts; i++) {
            auto w = skeleton.extractVertexWeights(weights,i);
            for(int j = 0; j < njoints; j++) {
                ws(j,i) = w[j];
            }
            std::cout << std::endl;
        }
        return ws;
    }

    void Problem::createModelWeights(const Model& initModel, int isFree)  {
        bool freeVars = !inputParams.freeze_weights;
        if(isFree != freeVars) {
            return;
        }


        if(freeVars) {
            barriers.shape = barriers.pose;
        } else {
            barriers.shape = barriers.endD3;
        }

        for(int i = 0; i < inputParams.ncomp+1; i++) {
            for(int j = 0; j < inputParams.nverts; j++) {
                Scalar3 vert;
                if(i == 0) {
                    vert.x = initModel.mean(0,j);
                    vert.y = initModel.mean(1,j);
                    vert.z = initModel.mean(2,j);
                } else {
                    vert.x = initModel.pcs[i-1](0,j); //distribution(generator);
                    vert.y = initModel.pcs[i-1](1,j);
                    vert.z = initModel.pcs[i-1](2,j);
                }
                inputData.addDynamic(vert, freeVars);
            }
        }

        barriers.joints = barriers.shape + (inputParams.ncomp+1)*inputParams.nverts;

        for (size_t i = 0; i < inputParams.njoints; i++) {
            Vector3 lu = initModel.joints.col(i);
            Scalar3 pt;
            pt.x = lu(0);
            pt.y = lu(1);
            pt.z = lu(2);
            inputData.addDynamic(pt, freeVars);
            //std::cout << lu.transpose() << std::endl;
        }

        for(int i = 0; i < inputParams.ncomp; i++) {
            for(int j = 0; j < inputParams.njoints; j++) {
                Vector3 pt = initModel.jointpcs[i].col(j);
                Scalar3 vert;
                vert.x = pt[0];//distribution(generator);
                vert.y = pt[1];//distribution(generator);
                vert.z = pt[2];//distribution(generator);
                inputData.addDynamic(vert, freeVars);
            }
        }

        barriers.deform = barriers.joints+(inputParams.ncomp+1)*inputParams.njoints;

        int defn = inputParams.ndeformjoints*
                   inputParams.ndeformshapes;
        for(int k = 0; k < defn; k++) {
            for(int i = 0; i < inputParams.nverts; i++) {
                Vector3 pt = initModel.deform[k].col(i);
                Scalar3 vert;
                vert.x = pt[0];//distribution(generator)
                vert.y = pt[1];//distribution(generator);
                vert.z = pt[2];//distribution(generator);
                inputData.addDynamic(vert, freeVars);
            }
        }


        barriers.poseMean = barriers.deform+inputParams.nverts*inputParams.ndeformshapes*inputParams.ndeformjoints;
        for(int j = 1; j < inputParams.njoints; j++) {
            Scalar3 vert;
            vert.x = 0;
            vert.y = 0;
            vert.z = 0;
            inputData.addDynamic(vert, freeVars);
        }

        if(freeVars) {
            barriers.pose = barriers.poseMean + (inputParams.njoints - 1);
        } else {
            barriers.endD3 = barriers.poseMean + (inputParams.njoints - 1);
        }

        if(freeVars) {
            barriers.weights = barriers.end1d;
        } else {
            barriers.weights = barriers.endD1;
        }

        for(int i = 0; i < inputParams.nverts; i++) {
            for(int j = 0; j < inputParams.njoints; j++) {
                auto wi = initModel.weights(j,i);
                if(wi < WMIN) wi = WMIN;
                if(wi > WMAX) wi = WMAX;
                auto y = (2*wi-1)*sqrt(wi*(1-wi))/(2*wi-2*wi*wi);
                //std::cout << LAMBDA*y << "\t";
                //std::cout << LAMBDA*y << "//" << wi << "\t";
                std::cout << wi << "\t";
                inputData.addDynamic(LAMBDA*y, freeVars);
            }
            std::cout << std::endl;
        }
        if(freeVars) {
            barriers.end1d = barriers.weights + inputParams.njoints*inputParams.nverts;
        } else {
            barriers.endD1 = barriers.weights + inputParams.njoints*inputParams.nverts;
        }



    }

    bool Problem::loadData(MeshTopology* top) {
        InputParams& ip = inputParams;

        eval.reset(new SubdivEvaluator(*top));
        Matrix3X samplingDebug;
        auto jointMap = calculateJoints(border_groups, &control_vertices);
        Matrix3X myJointPrior = calcJointInit();
        MatrixX skinningWeights =calcSkinningWeights(ip.njoints, ip.nverts);
        Model initModel;
        assert(ip.ndeformjoints*ip.ndeformshapes == 18 || ip.ndeformjoints*ip.ndeformshapes == 0);
        initModel.init(control_vertices,
                       myJointPrior, skinningWeights,
                       ip.ncomp,
                       ip.ndeformjoints*ip.ndeformshapes,
						ip.npersons, ip.nscans);
        if(ip.model_weights_path != "") {
            initModel.read(ip.model_weights_path);
        }
        eval->init(initModel.mean, &samplingDebug);

        bspline_control = eval->getBase();
        std::cout << "bspline control size " << bspline_control.size() << std::endl;
        dumpCloudA(ip.out_path  + std::string("/samplingDebug"), 0, samplingDebug);

        //TODO to config file
        const std::string base = ip.base_dataset_folder + "/reg/";  // /50026/";
        const std::string scanbase = ip.base_dataset_folder + "/scans/";
        const std::string maskbase = ip.base_dataset_folder + "/masks/";
        const std::string pose2dbase = ip.base_dataset_folder + "/pose2d/";

        const int m = inputParams.nscans;//corpora.personCount;

        //std::string path = base+corpora.reg_names[0];
        //std::cout << path << std::endl;
        std::string ext = ".xyz";//".obj";
        std::string ext2 = ".xyz";

        bmarkers.resize(ip.nscans*ip.nlabels);

        std::cout << marker_vertex_groups.size() << std::endl;
        std::cout << ip.nlabels << std::endl;
        assert(marker_vertex_groups.size() == ip.nlabels);
        //markerIdxs.resize(ip.nscans);
        static int iter = 0;
        //loadObj("dummy.obj", top, &v[0],NULL,SCALE);
        // TODO refactor
        for(int i = 0; i < m; i++) {

            int current_person = corpora.scan_to_person_id[i];
            model2Data[i].resize(ip.nverts);
            for(int j = 0; j < ip.nverts; j++) {
                model2Data[i][j] = 0;
            }

            std::cerr << "LOADING: " << i << " " << iter << std::endl;
            vmarkers[i].resize(3, marker_vertex_groups.size());
            std::cout << "DEBUG: " << i << " " << iter << std::endl;
            iter++;
            std::string loadpath = base+ corpora.folder[corpora.scan_to_person_id[i]] + "/"
                                   +corpora.reg_names[i] + ext;
			std::string labelpath;

			if(corpora.person_fist[current_person]) {
				labelpath = base + corpora.folder[corpora.scan_to_person_id[i]] + "/"
                                   +corpora.folder[corpora.scan_to_person_id[i]]+".txt";
			} else {
				labelpath = ip.template_data_marker_path;
			}
            std::string maskpath;

            std::vector<bool> scanMask;
            if(corpora.person_has_mask[current_person]) {
                maskpath = maskbase + corpora.folder[corpora.scan_to_person_id[i]] + "/"
                           +corpora.reg_names[i] + ".txt";
                readVector(maskpath, &scanMask);
            } else {
                //maskpath = maskbase + "csr4000a.txt";
            }

            std::string scanpath = scanbase + corpora.folder[corpora.scan_to_person_id[i]] + "/"
                                   +corpora.scan_names[i] + ext;
            std::cout << i << " try start loading " <<  loadpath << std::endl;

            Matrix3X markers, gtverts, points, normals;
            MeshTopology tempTop;
            //loadObj(loadpath, &topologies[i], &markers,NULL,SCALE); //,1, false

            std::cout << i << " try start loading " <<  scanpath  << " flipped: " << corpora.flipped[i] << std::endl;

            if(corpora.person_has_labels_3d[current_person]) {
                readCloud(loadpath, &gtverts, NULL, 1, corpora.flipped[i]);
                assert(gtverts.cols() ==  6890);
                //std::cout << gtverts.transpose() << std::endl;
                //exit(0);
            }
            readCloud(scanpath, &points, &normals, 1, corpora.flipped[i]); //,1, false

            if(scanMask.size() == 0) {
                for(int k = 0; k < points.cols(); k++) {
                    scanMask.push_back(true);
                }
            }
            assert(points.cols() == scanMask.size());

#ifdef USE_MARKERS
            //filter(&points,&normals,gtverts);
            std::cout << "normals: " << normals.cols() << "\t Masks" << scanMask.size() << std::endl;
            //assert(normals.cols() == scanMask.size());
            if(false) { //ip.nscans-2
                //filter(&points,&normals,scanMask);
            }

            //normals = points;
            //normals.colwise() -= ori;


            filter(&points,&normals,scanMask);
			//if(!corpora.person_fist[current_person]) {
				filterGroundPlane(&points,&normals,
					corpora.person_ground_axes[current_person],
					corpora.person_ground_planes[current_person]); //-.583
			//}
            fill(&points,&normals,ip.npoints);

#endif
            //corruptNonUniform(&points,&normals);
            std::cout << "Scanstats " << points.cols() << " " << normals.cols() << std::endl;

            //assert(points.cols() >= ip.npoints);
            assert(points.cols() == normals.cols());

            std::map<std::string, std::vector<size_t>> mynewgroup;

            //std::string labelpath = "template/markers_data.txt"; //"labels.txt.old"
#ifdef USE_MARKERS
			if(corpora.person_has_labels_3d[current_person]) {
            	loadVertexGroups(labelpath, &mynewgroup, NULL, NULL);
				int testSize = mynewgroup.size();
			}
#endif


            vmarkers[i].resize(3, ip.nlabels);

            std::cout << marker_vertex_groups.size() << std::endl;
            int idx = 0;
            {
                vmarkers[i].setZero();
                for(auto it = marker_vertex_groups.begin(); it != marker_vertex_groups.end(); ++it) {
                    bool markerExists = !(mynewgroup.find(it->first) == mynewgroup.end());
                    bmarkers[i*ip.nlabels+idx] = markerExists;
                    if(markerExists) {
                        //std::cout << "Marker Exists!" << std::endl;
                        if(corpora.person_has_labels_3d[current_person]) {
							int midx = mynewgroup[it->first][0];
                            vmarkers[i].col(idx) = gtverts.col(midx);
                        }
                        else {
                            vmarkers[i].col(idx).setConstant(0);
                        }
                    }
                    idx++;
                }
            }
            /*
            std::vector<int> idxlist;
            for(int j = 0; j < ip.nverts; j++) {
                idxlist.push_back(j);
            }
            std::random_shuffle(idxlist.begin(), idxlist.end());


            v[i].resize(3,ip.nlabels);
            for(int j = 0; j < ip.nlabels; j++) {
                v[i].col(j) = markers.col(idxlist[j]);
            }
            markerIdxs[i] = idxlist;
            */

            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(points.cols());
            perm.setIdentity();
            std::random_shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size());

            cloud[i] = (points * perm).leftCols(ip.npoints);
            cloudNorm[i] = (normals * perm).leftCols(ip.npoints);

            //greedyFarthestSampling(points, normals, ip.npoints, &cloud[i], &cloudNorm[i]);


            //saveObj(std::string("out/gtmeth")+ std::to_string(i) +std::string(".obj"),  top,  &v[i]);
            std::cout << "load success " << corpora.rot[i].transpose() << std::endl;
        }

        OpenSubdiv::Far::StencilTable const * localStencilTable = eval->getLocalStencilTable();
        localPointCount = localStencilTable->GetNumStencils();

        barriers.localStencil = BarrierIndex(DW);
        std::cout << "LOCAL STENCIL SIZE: " << localPointCount << std::endl;

        localIdx.resize(localPointCount*LOCAL_SIZE);
        for(int i = 0; i < localPointCount; i++) {
            auto st = localStencilTable->GetStencil(i);
            const int *ind = st.GetVertexIndices();
            float const *wei = st.GetWeights();

            std::cout << "LOCAL SIZE: " << st.GetSize() << std::endl;
            assert(st.GetSize() <= LOCAL_SIZE);

            std::array<Scalar,LOCAL_SIZE> local;
            for(int j = 0; j < LOCAL_SIZE; j++) {
                local[j] = 0; //set to arbitrary local idx with weight 0
                localIdx[i*LOCAL_SIZE+j] = ind[0];
            }

            for(int j = 0; j < st.GetSize(); j++) {
                local[j] = wei[j];
                localIdx[i*LOCAL_SIZE+j] = ind[j];
                assert(ind[j] >= 0 && ind[j] < ip.nverts);
            }
            inputData.addPoint(local);
        }

        //localBarrier...

        //barriers.shapePrior = 0;
        //TODO strange!
        barriers.labelBase = BarrierIndex(D3);
        barriers.pointBase = BarrierIndex(D3);


        for(int i = 0; i < m; i++) {
            std::cout << "Markers " << i << " : " <<  vmarkers[i].cols() << std::endl;
            for(int j = 0; j < ip.nlabels/*v[i].cols()*/; j++) {
                Scalar3 pt;
                pt.x = vmarkers[i](0,j);
                pt.y = vmarkers[i](1,j);
                pt.z = vmarkers[i](2,j);
                inputData.addPoint(pt);

                //barriers.shapePrior++;
	            barriers.pointBase++;
            }
        }

        int maxSize = -1;
        srand(time(NULL));

        std::shared_ptr<RayTracer> rayTracer;
        rayTracer.reset(new RayTracer());
        barriers.surface = BarrierIndex(D2);
        for(int i = 0; i < m; i++) {

            std::ostringstream ostr;
            //ostr << std::internal << std::setfill('0') << std::setw(5) << i;
            ostr << corpora.folder[corpora.scan_to_person_id[i]] << "."
                                   << corpora.scan_names[i];
			std::string folder = ip.out_path + std::string("/") + ostr.str();
            mkdir(folder.c_str(),0777);
            folder = folder + "/";

            //cloud[i].resize(3,ip.npoints);
            //cloudNorm[i].resize(3,ip.npoints); cloudNorm[i].setZero();
            Matrix3X subdivCloud; subdivCloud.resize(3,ip.npoints);

            for(int p = 0; p < ip.npoints; p++) {
                int quadId = rand() % (eval->ptexnum);
                face_idx.push_back(quadId);
                patch_idx.push_back(0);
                tangBasis.push_back(Matrix22());
                Scalar2 u;
                u.x = .5;
                u.y = .5;
                inputData.addPoint(u);

                Scalar3 vert;
                vert.x = cloud[i](0,p);
                vert.y = cloud[i](1,p);
                vert.z = cloud[i](2,p);

                inputData.addPoint(vert);

                Scalar3 norm;
                norm.x = cloudNorm[i](0,p);
                norm.y = cloudNorm[i](1,p);
                norm.z = cloudNorm[i](2,p);

                inputData.addPoint(norm);

            }

            dumpCloudA(folder + std::string("labels")+std::to_string(i), i, vmarkers[i]);

            std::cout << "scan: " << i << std::endl;
            dumpCloudA(folder + std::string("cloud")+std::to_string(i), i, cloud[i],&cloudNorm[i]);

            surfaceState.addData(&(*eval), rayTracer,  top, &cloud[i], &cloudNorm[i]);
            std::cout << "Cloud Data configured: " << i << std::endl;
        }

		barriers.surfaceScale =  barriers.surface + ip.npoints*ip.nscans;
		for(int i = 0; i < m; i++) {
            for(int p = 0; p < ip.npoints; p++) {
                Scalar2 u;
                u.x = 1.0;
                u.y = 1.0;
                inputData.addPoint(u); //step scale
            }
        }


        barriers.surfaceLabel = barriers.surfaceScale + ip.npoints*ip.nscans;
        assert(marker_vertex_groups.size() == ip.nlabels);
        for(auto markerIt = marker_vertex_groups.begin(); markerIt != marker_vertex_groups.end(); ++markerIt) {
            double u,v;
            int fidx, pidx;
            Vector2 pt;
            std::cout << "Before coord " << std::endl;
            top->vertexToSurfaceCoordinate(markerIt->second[0], fidx, u, v);
            pt[0] = u; pt[1] = v;
            pt = eval->normalize(fidx, pt, &pidx);
            //TODO convert base!
            Scalar2 vert;
            vert.x = pt[0];
            vert.y = pt[1];
            std::cout << vert.x << "," << vert.y << std::endl;
            inputData.addPoint(vert);
        }
        //vertexLabel has to come after surfaceLabel
        barriers.vertexLabel = barriers.surfaceLabel + inputParams.nlabels;
        for(int i = 0; i < ip.nverts; i++) {

            int fidx;
            double u,v;
            top->vertexToSurfaceCoordinate(i,fidx,u,v);
            Scalar2 vert;
            vert.x = u;
            vert.y = v;
            inputData.addPoint(vert);
        }

        barriers.label2d = barriers.vertexLabel + inputParams.nverts;
        for(int v = 0; v < ip.nviews; v++) {
            for(int i = 0; i < m; i++) {
                int p = corpora.scan_to_person_id[i];
                std::string pose2d_path = pose2dbase + corpora.folder[p] + "/"
                                          + corpora.scan_names[i] + std::string("_") +
                                          std::to_string(v) + "_keypoints.json";
                Pose2d pose2d(ip.template_surface_map_path);
                int current_person = corpora.scan_to_person_id[i];
                if(corpora.person_has_labels_2d[current_person])
                    pose2d.read(pose2d_path);
                for(int j = 0; j <  ip.njointlabels; j++) {
					if(!corpora.person_has_labels_2d[current_person]) {
						Scalar2 vert;
                    	vert.x = 0;
                    	vert.y = 0;
                    	inputData.addPoint(vert);
						continue;
					}
                    std::string joint_name = skeleton.getJoint(j).first;
                    Vector3 joint = pose2d.getJoint(joint_name);
                    Scalar2 vert;
                    vert.x = joint[0];
                    vert.y = joint[1];
                    inputData.addPoint(vert);
                }
            }
        }
        barriers.labelSurface2d = barriers.label2d +
                                  ip.nviews*ip.nscans*ip.njointlabels;
        for(int v = 0; v < ip.nviews; v++) {
            for(int i = 0; i < m; i++) {
                std::string pose2d_path = pose2dbase + corpora.folder[corpora.scan_to_person_id[i]] + "/"
                                          + corpora.scan_names[i] + std::string("_") +
                                          std::to_string(v) + "_keypoints.json";
                Pose2d pose2d(ip.template_surface_map_path);
                int current_person = corpora.scan_to_person_id[i];
                if(corpora.person_has_labels_2d[current_person])
                    pose2d.read(pose2d_path);
                for(int j = 0; j <  ip.nverts; j++) {
					if(!corpora.person_has_labels_2d[current_person]) {
						Scalar2 vert;
                    	vert.x = 0;
                    	vert.y = 0;
                    	inputData.addPoint(vert);
						continue;
					}
                    Vector3 surf = pose2d.getSurface(j);
                    Scalar2 vert;
                    vert.x = surf[0];
                    vert.y = surf[1];
                    inputData.addPoint(vert);
                }
            }
        }



        //exit(1);
        std::cout << "Init Subdiv weights " << std::endl;

        barriers.subdivLabelWeights = BarrierIndex(D1);

        barriers.labelMask = barriers.subdivLabelWeights;
        for(int i = 0; i < ip.nscans; i++) {
            for(int j = 0; j < ip.nlabels; j++) {
                inputData.addPoint(bmarkers[i*ip.nlabels+j] ? 1 : 0);
            }
        }
        barriers.singleZero = barriers.labelMask + ip.nlabels*ip.nscans;
        barriers.singleOne = barriers.singleZero + 1;
		barriers.singleTwo = barriers.singleOne + 1;
        inputData.addPoint(0);
        inputData.addPoint(1);
		inputData.addPoint(2);
		barriers.groundPlanes = barriers.singleTwo + 1;
		for(int i = 0; i < ip.npersons; i++) {
			inputData.addPoint(corpora.person_ground_planes[i]);
		}

        barriers.personCount = barriers.groundPlanes + ip.npersons;
        for(int i = 0; i < ip.npersons; i++) {
            inputData.addPoint(corpora.scanCount[i]);
        }

        barriers.labelWeights2d = barriers.personCount + ip.npersons;
        for(int v = 0; v < ip.nviews; v++) {
            for(int i = 0; i < m; i++) {
                std::string pose2d_path = pose2dbase + corpora.folder[corpora.scan_to_person_id[i]] + "/"
                                          + corpora.scan_names[i] + std::string("_")
                                          + std::to_string(v) + "_keypoints.json";
                Pose2d pose2d(ip.template_surface_map_path);
                int current_person = corpora.scan_to_person_id[i];
                if(corpora.person_has_labels_2d[current_person]) {
                    pose2d.read(pose2d_path);
                }

                for(int j = 0; j <  ip.njointlabels; j++) {
					if(!corpora.person_has_labels_2d[current_person]) {
						inputData.addPoint(0);
						continue;
					}
                    std::string joint_name = skeleton.getJoint(j).first;
                    Vector3 joint = pose2d.getJoint(joint_name);
					Scalar score = joint[2];
					//score = score > .9 ? score : 0;
                    inputData.addPoint(score);
                }
            }
        }
        barriers.labelSurfaceWeights2d = barriers.labelWeights2d +
                                         ip.nviews*ip.njointlabels*ip.nscans;
        for(int v = 0; v < ip.nviews; v++) {
            for(int i = 0; i < m; i++) {
                std::string pose2d_path = pose2dbase + corpora.folder[corpora.scan_to_person_id[i]] + "/"
                                          + corpora.scan_names[i] + std::string("_")
                                          + std::to_string(v) + "_keypoints.json";
                Pose2d pose2d(ip.template_surface_map_path);
                int current_person = corpora.scan_to_person_id[i];
                if(corpora.person_has_labels_2d[current_person]) {
                    pose2d.read(pose2d_path);
                }

                for(int j = 0; j <  ip.nverts; j++) {
					if(!corpora.person_has_labels_2d[current_person]) {
						Vector3 surf; surf.setConstant(0);
						inputData.addPoint(surf[2]);
						continue;
					}
                    Vector3 surf = pose2d.getSurface(j);
                    inputData.addPoint(surf[2]);
                }
            }
        }

        barriers.metric = barriers.labelSurfaceWeights2d
                          + ip.nviews*ip.nverts*ip.nscans;
        //std::cout << "Check boundary sanity " << barriers.metric << " : " << bsplines_idx.size() << std::endl;

        VectorX M1 = 2 * loadVector(ip.template_laplace_path); //.1 2
        std::cout << M1 << std::endl;
        assert(M1.size() == top->quads.cols()*16);
        int pointer = 0;
        for(int i = 0; i < top->quads.cols(); i++) {
            for(int j = 0; j < 16; j++) {
                //std::cout << M1[pointer] << std::endl;
                //Scalar w = (((j%4) == (j/4)) ? 2.5 : 0); //0.001
                //inputData.addPoint(w); //M1[pointer++]
                inputData.addPoint(M1[pointer++]);
                //M1[pointer++] = w;
                //inputData.addPoint(M1[pointer++]);
            }
        }

        std::cout << "end: " << (barriers.metric+ M1.size()) << std::endl;


        //assert(pointer == M1.size()-1);

        barriers.weightInit = barriers.metric + top->quads.cols()*16;

        for(int i = 0; i < inputParams.nverts; i++) {
            for(int j = 0; j < inputParams.njoints; j++) {
                auto wi = initModel.weights(j,i);
                if(wi < WMIN) wi = WMIN;
                if(wi > WMAX) wi = WMAX;
                //Scalar y = (2*wi-1)*sqrt(wi*(1-wi))/(2*wi-2*wi*wi);

                //inputData.addPoint(LAMBDA*y);
                inputData.addPoint(wi);
                std::cout << wi << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "MaxSize: "<< maxSize << std::endl;

        barriers.surfaceBase = barriers.pointBase + 2*m*ip.npoints;

        /*
        for(int i = 0; i < m; i++) {
            std::vector<SurfacePoint> surfacePoints;
            for(int p = 0; p < ip.npoints; p++) {
                Scalar3 vert;
                vert.x = 0;
                vert.y = 0;
                vert.z = 0;             scan_names
                inputData.addPoint(vert);
                inputData.addPoint(vert);
                inputData.addPoint(vert);
                inputData.addPoint(vert);
                inputData.addPoint(vert);
            }
        }*/
        barriers.poseBase = barriers.surfaceBase;// + 5*m*ip.npoints;

        for(int i = 0; i < 2*inputParams.nscans; i++) {
            int current_person = corpora.scan_to_person_id[i % inputParams.nscans];
            for(int j = 0; j < inputParams.njoints; j++) { //+trans

                Scalar3 vert;
                Vector3 rot;
				rot = initModel.scans[i%inputParams.nscans].posebase.col(j);
                //rot.setZero();
				assert(initModel.scans[i%inputParams.nscans].posebase.cols() == inputParams.njoints);
				std::cout << rot << std::endl;
                vert.x = rot[0];
                vert.y = rot[1];
                vert.z = rot[2];
                inputData.addPoint(vert);
            }
			std::cout << std::endl;
        }
		//exit(1);
        barriers.absolutePoseBase = barriers.poseBase + ip.nscans*ip.njoints;
        barriers.poseBaseMean =  barriers.absolutePoseBase + ip.nscans*ip.njoints;
        {
            for(int j = 1; j < inputParams.njoints; j++) { //+trans
                Scalar3 vert;
                vert.x = initModel.meanPose(0,j-1);
                vert.y = initModel.meanPose(1,j-1);
                vert.z = initModel.meanPose(2,j-1);
                inputData.addPoint(vert);
            }
        }

        barriers.zero = barriers.poseBaseMean + (inputParams.njoints - 1);
        {
            Scalar3 vert;
            vert.x = 0;
            vert.y = 0;
            vert.z = 0;

            inputData.addPoint(vert);
        }

        barriers.one = barriers.zero + 1;

        {
            Scalar3 vert;
            vert.x = 1;
            vert.y = 1;
            vert.z = 1;

            inputData.addPoint(vert);
        }

        barriers.two = barriers.one + 1;
        {
            Scalar3 vert;
            vert.x = 2;
            vert.y = 2;
            vert.z = 2;

            inputData.addPoint(vert);
        }
        barriers.three = barriers.two + 1;
        {
            Scalar3 vert;
            vert.x = 3;
            vert.y = 3;
            vert.z = 3;

            inputData.addPoint(vert);
        }

        barriers.cameraViews = barriers.three + 1;
		
		for(int camId = 0; camId < ip.camPaths.size(); camId++) {
        for(int v = 0; v < ip.nviews; v++) {
			//std::string camPath = corpora.person_fist[current_person] ? 
			//	"templates/caesar_cameras.json" : "templates/cameras.json"; 
            Camera cam(ip.camPaths[camId], v);
            Vector3 pos = cam.getCameraPosition();
            Vector3 ori  = cam.getCameraOrientation();
            Scalar3 p = toScalar3(pos);
            Scalar3 o = toScalar3(ori);
            inputData.addPoint(p);
            inputData.addPoint(o);
        }
		}

        barriers.meanShapePrior = barriers.cameraViews + 2*ip.nviews*ip.camPaths.size();
        for(int i = 0; i < inputParams.nverts; i++) {
            Scalar3 pt;
            pt.x = initModel.mean(0,i);
            pt.y = initModel.mean(1,i);
            pt.z = initModel.mean(2,i);
            inputData.addPoint(pt);
        }
        saveObj(ip.out_path +std::string("/meanShapePrior.obj"), top, &initModel.mean);

        std::cout << "End Data Insertion: " << (barriers.meanShapePrior + inputParams.nverts) << std::endl;


        barriers.jointPrior = barriers.meanShapePrior + inputParams.nverts;
        //auto jointMap = calculateJoints(border_groups, &control_vertices);

        for (size_t i = 0; i < inputParams.njoints; i++) {
            Vector3 lu = initModel.joints.col(i);
            Scalar3 pt;
            pt.x = lu(0);
            pt.y = lu(1);
            pt.z = lu(2);
            inputData.addPoint(pt);
        }
        dumpCloudA(ip.out_path + std::string("/jointPrior"), 0, initModel.joints);

        barriers.endD1 = barriers.weightInit + inputParams.nverts*inputParams.njoints;
        barriers.endD2 = barriers.labelSurface2d + ip.nviews*ip.nscans*ip.nverts;
        barriers.endD3 = barriers.jointPrior + (inputParams.njoints);
        barriers.endDW = barriers.localStencil + localPointCount;

        createModelWeights(initModel, false); //DO NOT MOVE THIS CALL

        inputData.verifyDataSize(
                barriers.endD1,
                barriers.endD2,
                barriers.endD3,
                barriers.endDW
        );
        return true;


        //barriers.poseBase = barriers.jointPrior+inputParams.njoints-1;
    }

    void Problem::createUnknowns(MeshTopology* top) {
        InputParams& ip = inputParams;
        auto jointMap = calculateJoints(border_groups, &control_vertices);
        MatrixX skinningWeights =calcSkinningWeights(ip.njoints, ip.nverts);
        Matrix3X myJointPrior = calcJointInit();
        Model initModel;
        initModel.init(control_vertices, myJointPrior,
                       skinningWeights, ip.ncomp,
                       ip.ndeformjoints*ip.ndeformshapes,
					   ip.npersons,ip.nscans);
        if(ip.model_weights_path != "") {
            initModel.read(ip.model_weights_path);
        }

        //std::array<Scalar,PATCH_SIZE> patch;
        //for(int i = 0; i < PATCH_SIZE; i++) {
        //    patch[i] = 0;
        //}
        //inputData.addUnknown(patch);

        /*{
            Scalar2 vert;
            vert.x = 0;
            vert.y = 0;

            inputData.addUnknown(vert);
        }*/
        barriers.surfaceUnknown = BarrierIndex(M2);
        barriers.robust = BarrierIndex(M1);


        for(int i = 0; i < inputParams.nscans; i++) {
            for(int j = 0; j < inputParams.npoints; j++) {
                Scalar2 vert;
                vert.x = 0;
                vert.y = 0;

                inputData.addUnknown(vert);

                Scalar rob = 1;
                inputData.addUnknown(rob);
            }
        }
        barriers.end2d = barriers.surfaceUnknown +
                         inputParams.nscans * inputParams.npoints;

        barriers.registration = BarrierIndex(M3);

        for(int i = 0; i < inputParams.nscans; i++) {

            std::vector<Matrix3> localmat(1);

            Matrix3X angle(3,1);
            Matrix3X euler(3,1);
            angle.col(0) = corpora.rot[i];
            angle2euler(angle,euler);

            rod2Mat(euler,localmat);
            Matrix3X rotated = localmat[0]*initModel.mean;

            Vector3 rmean = cloud[i].rowwise().mean() - rotated.rowwise().mean();
            Scalar mx, my, mz;
            mx = rmean[0];
            my = rmean[1];
            mz = rmean[2];

            for(int j = 0; j < inputParams.nverts; j++) {
                rotated(0,j) += mx;
                rotated(1,j) += my;
                rotated(2,j) += mz;
            }

            for(int j = 0; j < inputParams.nverts; j++) {
                Scalar3 vert;
                vert.x = rotated(0,j);
                vert.y = rotated(1,j);
                vert.z = rotated(2,j);
                inputData.addUnknown(vert);
            }
            //surfaceState.setCurrentModel(i,rotated, false);

            //saveObj(std::string("out/Ta")+ std::to_string(iter) + std::string("_") + std::to_string(i)+std::string(".obj"), &top, &rotated);

            /*
            for(int j = 0; j < inputParams.nverts; j++) {
                Scalar3 vert;
                vert.x = vreal[i](0,j);
                vert.y = vreal[i](1,j);
                vert.z = vreal[i](2,j);
                inputData.addUnknown(vert);
            }
            surfaceState.setCurrentModel(i,vreal[i], false);
            */
        }

        barriers.localRegistration = barriers.registration + inputParams.nscans * inputParams.nverts;

        OpenSubdiv::Far::StencilTable const * localStencilTable = eval->getLocalStencilTable();
        std::cout << "LOCAL STENCIL SIZE: " << localPointCount << std::endl;

        for(int k = 0; k < inputParams.nscans; k++) {
            std::cout << "Compute local points for scan: " << k << std::endl;
            Matrix3X T; T.resize(3,inputParams.nverts);
            inputData.read(T,barriers.registration + k*inputParams.nverts,inputParams.nverts);
            Matrix3X localPts(3,localPointCount);
            for(int i = 0; i < localPointCount; i++) {
                auto st = localStencilTable->GetStencil(i);
                //std::cout <<  st.GetSize() << std::endl;
                auto *ind = st.GetVertexIndices();
                auto *wei = st.GetWeights();
                Vector3 localPt; localPt.setZero();
                for(int j = 0; j < st.GetSize(); j++) {
                    //std::cout <<  wei[j] << " " << ind[j] << std::endl;
                    localPt += wei[j]*T.col(ind[j]);
                }
                Scalar3 val;
                val.x = localPt[0];
                val.y = localPt[1];
                val.z = localPt[2];
                inputData.addUnknown(val);
                localPts.col(i) = localPt;
            }
            surfaceState.setCurrentModel(k, T, localPts, true); //false

            //saveObj(std::string("out/gtmeth")+ std::to_string(k) +std::string(".obj"),  top,  &T);
            //dumpCloudA(std::string("out/p"), k, localPts);
            //Matrix3X direct = eval->localPoints(T);
            //dumpCloudA(std::string("out/direct"), k, direct);
            //Matrix3X custom = eval->evalPointsCustom(T);
            //dumpCloudA(std::string("out/customF"), k, custom);
            //dumpCloudA(std::string("out/customT"), k, T);
        }

        barriers.debugRegistration = barriers.localRegistration + inputParams.nscans*localPointCount;

        for(int i = 0; i < inputParams.nlabels*inputParams.nscans; i++) {
            Scalar3 val;
            val.x = 0;
            val.y = 0;
            val.z = 0;
            inputData.addUnknown(val);
        }

        barriers.coeff = barriers.robust + inputParams.nscans*inputParams.npoints;
        for(int i = 0; i < inputParams.npersons; i++) {
            for(int j = 0; j < inputParams.ncomp; j++) {
                Scalar coeff = initModel.persons[i].coeff[j];
				std::cout << coeff << ",";
                inputData.addUnknown(coeff);
            }
			std::cout << std::endl;
        }
		//exit(1);

        barriers.pose = barriers.debugRegistration + inputParams.nlabels*inputParams.nscans;
        barriers.end1d = barriers.coeff + inputParams.npersons*inputParams.ncomp;

        createModelWeights(initModel, true);

        for(int i = 0; i < 2*inputParams.nscans; i++) {
            std::vector<Matrix3> localmat(1);
            int scan_id = i % inputParams.nscans;
			int current_person = corpora.scan_to_person_id[scan_id];


            Matrix3X angle(3,1);
            Matrix3X euler(3,1);
            angle.col(0) = corpora.rot[scan_id];
            angle2euler(angle,euler);


            rod2Mat(euler,localmat);
            Matrix3X rotated = localmat[0]*control_vertices;

            auto mean = cloud[scan_id].rowwise().mean() - rotated.rowwise().mean();
            //std::cout << mean.rows() << " : " << mean.cols() << std::endl;
            Scalar3 vert;
		
			
			if(corpora.person_fist[current_person]) {
				//vert.x = 0;
                //vert.y = 0;
                //vert.z = -1;
				Vector3 tr = initModel.scans[scan_id].trans.col(0);
                vert.x = tr[0];
                vert.y = tr[1];
                vert.z = tr[2];

            } else {
            	//vert.x = mean(0,0);//distribution(generator);
	            //vert.y = mean(1,0);//distribution(generator);
    	        //vert.z = mean(2,0);//distribution(generator);
				Vector3 tr = initModel.scans[scan_id].trans.col(0);
				vert.x = tr[0];
				vert.y = tr[1];
				vert.z = tr[2];
            }

            inputData.addUnknown(vert);
            for(int j = 1; j < inputParams.njoints+1; j++) { //+trans
                Scalar3 vert;
                vert.x = 0;//distribution(generator);
                vert.y = 0;//distribution(generator);
                vert.z = 0;//distribution(generator);
                inputData.addUnknown(vert);
            }
        }

        barriers.absolutePose = barriers.pose + (inputParams.njoints+1)*inputParams.nscans;


        barriers.pjoint  = barriers.absolutePose + (inputParams.njoints+1)*inputParams.nscans;

        for(int j = 0; j < inputParams.npersons + inputParams.nscans; j++) {
            for (size_t i = 0; i < inputParams.njoints; i++) {
                Vector3 lu = initModel.joints.col(i);
                Scalar3 pt;
                pt.x = lu(0);
                pt.y = lu(1);
                pt.z = lu(2);
                inputData.addUnknown(pt);
                //std::cout << lu.transpose() << std::endl;
            }
        }

        barriers.sjoint = barriers.pjoint + inputParams.njoints*inputParams.npersons;
        barriers.pshape = barriers.sjoint + inputParams.njoints*inputParams.nscans;
        for(int i = 0; i < inputParams.nscans; i++) {
            for(int j = 0; j < inputParams.nverts; j++) {
                Scalar3 vert;
                vert.x = initModel.mean(0,j);
                vert.y = initModel.mean(1,j);
                vert.z = initModel.mean(2,j);
                inputData.addUnknown(vert);
            }
        }
        barriers.end3d = barriers.pshape + (inputParams.nverts)*inputParams.nscans;
        //barriers.end3d = barriers.pjoint + (inputParams.njoints)*inputParams.npersons; //inputParams.njoints+1
        //std::cout << barriers.end3d << ":" << std::endl;
        //exit(-1);

        inputData.verifySize(
                barriers.end3d,
                barriers.end2d,
                barriers.end1d
        );

    }


    void Problem::createReducedSubdivTerm(MeshTopology* top) {
        InputParams& ip = inputParams;
        int edgeSize = 6+16; //4+16;
        std::cout << "Reduced Term Edge Size: " << edgeSize << std::endl;
        std::vector<BarrierIndex> edge(edgeSize);
        int tid = termIdx++;
        for(int i = 0; i < ip.nscans; i++) {
            for(int p = 0; p < ip.npoints; p++) {
                int offset3d = 0;
                edge[offset3d++] = barriers.pointBase + 2*i*ip.npoints + 2*p+0;
                edge[offset3d++] = barriers.pointBase + 2*i*ip.npoints + 2*p+1;

                edge[offset3d++] = barriers.surfaceUnknown + i*ip.npoints + p;
                edge[offset3d++] = barriers.surface + i*ip.npoints + p;
				edge[offset3d++] = barriers.surfaceScale + i*ip.npoints + p;
                edge[offset3d++] = barriers.robust + i*ip.npoints + p;

                int quadId = face_idx[i*ip.npoints+p];
                auto quad = top->quads.col(quadId);
                for(int q = 0; q < 16; q++) {
                    //q is just a dummy
                    edge[offset3d++] = barriers.localRegistration + i*localPointCount + q;
                }

                assert(edge.size() == edgeSize);
                terms[tid].add(edge);
            }
        }

    }

    void Problem::createReverseDataEdges(MeshTopology* top, int tid) {
        InputParams& ip = inputParams;
        int edgeSize = 4+16;
        std::vector<BarrierIndex> edge(edgeSize);
        for(int i = 0; i < ip.nscans; i++) {
            for(int v = 0; v < ip.nverts; v++) {
                int offset3d = 0;
                double u1,u2;
                int fidx, pidx;
                fidx = pidx = -1;
                top->vertexToSurfaceCoordinate(v, fidx, u1, u2);
                Vector2 pt;
                pt[0] = u1;
                pt[1] = u2;
                eval->normalize(fidx, pt, &pidx); //is that  required?
                //set by reversed ICP
                edge[offset3d++] = barriers.pointBase + 2*i*ip.npoints;
                edge[offset3d++] = barriers.singleOne;
                edge[offset3d++] = barriers.singleOne; //switch to indicate reverse
                edge[offset3d++] = barriers.vertexLabel + v;

                std::cout << top->num_faces() << "," << i << ","
                          << v << "," << fidx << "," << pidx << std::endl;

                for(int q = 0; q < 16; q++) {
                    int offset = bspline_control.at(16*pidx + q);
                    assert(offset < localPointCount);
                    edge[offset3d++] = barriers.localRegistration + i*localPointCount + offset;
                }
                terms[tid].add(edge);
            }
        }
    }

    void Problem::createLabelTerm(MeshTopology* top) {
        InputParams& ip = inputParams;
        int edgeSize = 4+16;
        std::cout << "Bilinear Term Edge Size: " << edgeSize << std::endl;
        std::vector<BarrierIndex> edge(edgeSize);
        int tid = termIdx++;

        int seqId = 0;
		int scanLabelCount = 0;
        //assert(corpora.markerFrames <= ip.nscans);
        //int size = corpora.markerFrames < ip.nscans ? corpora.markerFrames : ip.nscans;
        for(int i = 0; i < ip.nscans; i++) {
            int current_person = corpora.scan_to_person_id[i];
            if(!corpora.person_has_labels_3d[current_person])
                continue;
            auto markerIt = marker_vertex_groups.begin();
            std::cout << "Label prior for t = " << i << std::endl;
            for(int p = 0; p < ip.nlabels; p++) {

                double u,v;
                int fidx;
                int pidx;
                top->vertexToSurfaceCoordinate(markerIt->second[0], fidx, u, v);
                Vector2 pt;
                pt[0] = u;
                pt[1] = v;
                eval->normalize(fidx, pt, &pidx); //is that  required?
                int offset3d = 0;
                //edge[offset3d++] = barriers.debugRegistration + i*ip.nlabels + p;
                edge[offset3d++] = barriers.labelBase + i*ip.nlabels + p;
                edge[offset3d++] = barriers.labelMask + i*ip.nlabels + p;
                edge[offset3d++] = barriers.singleZero;
                edge[offset3d++] = barriers.surfaceLabel + p;
                //int quadId = face_idx[i*ip.npoints+p];
                //auto quad = top->quads.col(quadId);
                for(int q = 0; q < 16; q++) {
                    int offset = bspline_control.at(16*pidx + q);
                    assert(offset < localPointCount);
                    edge[offset3d++] = barriers.localRegistration + i*localPointCount + offset;
                }
                terms[tid].add(edge);
                assert(edge.size() == edgeSize);
                ++markerIt;
				scanLabelCount++;
            }
        }
        reverse_offset = edgeSize*scanLabelCount*ip.nlabels;//terms[1].count();
        if(ip.useModel2Data)
            createReverseDataEdges(top, tid);
    }

    void Problem::createLabel2dTerm() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = 7;
        std::vector<BarrierIndex> edge(edgeSize);
        for(int v = 0; v < ip.nviews; v++) {
            for(int i = 0; i < ip.nscans; i++) {
                int current_person = corpora.scan_to_person_id[i];
                if(!corpora.person_has_labels_2d[current_person])
                    continue;
                for(int j = 0; j < ip.njointlabels; j++) {
                    int offset3d = 0;
                    int joint = j;
					int camId = corpora.person_fist[current_person] ? 1 : 0; 
					std::cout << camId << std::endl;
                    edge[offset3d++] = barriers.label2d + v*ip.nscans*ip.njointlabels + (ip.njointlabels*i) + j;
                    edge[offset3d++] = barriers.labelWeights2d + v*ip.nscans*ip.njointlabels + (ip.njointlabels*i) + j;
                    edge[offset3d++] = barriers.sjoint + (ip.njoints)*i + joint;
                    edge[offset3d++] = barriers.absolutePose + (ip.njoints+1)*i;
                    edge[offset3d++] = barriers.singleZero;
                    edge[offset3d++] = barriers.cameraViews+2*ip.nviews*camId +2*v;
                    edge[offset3d++] = barriers.cameraViews+ 2*ip.nviews*camId   +2*v+1;
                    terms[tid].add(edge);
                }

                for(int j = 0; j < ip.nverts; j++) {
                    int offset3d = 0;
					int camId = corpora.person_fist[current_person] ? 1 : 0;
                    edge[offset3d++] = barriers.labelSurface2d + v*ip.nscans*ip.nverts + (ip.nverts*i) + j;
                    edge[offset3d++] = barriers.labelSurfaceWeights2d + v*ip.nscans*ip.nverts + (ip.nverts*i) + j;
                    edge[offset3d++] = barriers.registration + i*ip.nverts + j;
                    edge[offset3d++] = barriers.absolutePose + (ip.njoints+1)*i;
                    edge[offset3d++] = barriers.singleOne;
                    edge[offset3d++] = barriers.cameraViews+2*ip.nviews*camId +2*v;
                    edge[offset3d++] = barriers.cameraViews+2*ip.nviews*camId +2*v+1;
                    terms[tid].add(edge);
                }

            }
        }
    }
    void Problem::createRestModelTerm() {
		modelTermIndex = termIdx;
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = REST_SIZE;
        std::cout << "Model Term Edge Size: " << edgeSize << std::endl;
        std::vector<BarrierIndex> edge(edgeSize);
        int offset3d = 0;
        int person = 0;
        int relativeScan = 0;
        for(int i = 0; i < ip.nscans; i++) {
            std::cout <<  "Scan Id: " << i << "\tPerson Id: " << person << std::endl;
            for(int j = 0; j < ip.nverts; j++) {
                offset3d = 0;
                //edge[offset3d++] = i*ip.nverts + j;
                edge[offset3d++] = barriers.pshape + i*ip.nverts + j;
                edge[offset3d++] = barriers.shape + j;
                for(int k = 0; k < ip.ncomp; k++) {
                    edge[offset3d++] = barriers.shape + (k+1)*ip.nverts + j;
                }
                for(int k = 0; k < ip.ndeformjoints*ip.ndeformshapes; k++) {
                    edge[offset3d++] = barriers.deform + k*ip.nverts + j;
                }
                std::cout << i << ": ";
                for(int k = 0; k < ip.ndeformjoints; k++) {
                    int joint = closestDeformJoints(k,j); //k+1
                    assert(joint != 0);
                    assert(ip.ndeformjoints < 4);
                    std::cout << joint << ",";
                    edge[offset3d++] = barriers.poseBase + (ip.njoints)*i + joint;
                    edge[offset3d++] = barriers.pose + (ip.njoints+1)*i + joint + 1;
					bool fist = corpora.person_fist[person];
                    edge[offset3d++] = !skeleton.isRigid(joint, fist) ? barriers.one : barriers.zero;
                }
                std::cout << std::endl;
                for(int k = 0; k < ip.ncomp; k++) {
                    edge[offset3d++] = barriers.coeff + ip.ncomp*person+k;
                }
                edge[offset3d++] = i < corpora.markerFrames ?  barriers.one : barriers.zero;
                //std::cout << "Size Rest Model Debug: " << offset3d << "/" << edge.size() << std::endl;
                assert(edge.size() == edgeSize && edgeSize == offset3d);
                terms[tid].add(edge);
            }
            relativeScan++;
            if(relativeScan == corpora.scanCount[person]) {
                relativeScan = 0;
                person++;
            }
        }
    }

    bool Problem::jointInInfluence(int vid, int joint) {
        SVector4 influence = closestJoints.col(vid).cast<int>();
        int* regIdx_ptr = std::find(&influence[0], &influence[0]+4, joint);
        int offset = regIdx_ptr - &influence[0];
        return regIdx_ptr != &influence[0]+4;
    }
    void Problem::createModelTerm() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = MODEL_SIZE;
        std::cout << "Model Term Edge Size: " << edgeSize << std::endl;
        std::vector<BarrierIndex> edge(edgeSize);
        int offset3d = 0;
        int person = 0;
        int relativeScan = 0;
        for(int i = 0; i < ip.nscans; i++) {
            std::cout <<  "Scan Id: " << i << "\tPerson Id: " << person << std::endl;
            for(int j = 0; j < ip.nverts; j++) {
                offset3d = 0;
                edge[offset3d++] = barriers.registration + i*ip.nverts + j;
                edge[offset3d++] = barriers.pshape + i*ip.nverts + j;
                int weightDebug = 0;
                for(int joint = 0; joint < ip.njoints; joint++) {
                    if(jointInInfluence(j,joint)) {
                        edge[offset3d++] = barriers.weights + j*ip.njoints + joint;
                        weightDebug++;
                    }
                }
				if(weightDebug != 4) {
					std::cout << j << " : " << weightDebug << std::endl;
				}
                assert(weightDebug == 4);
                for(int joint = 0; joint < ip.njoints; joint++) {
                    if(jointInInfluence(j,joint)) {
                        edge[offset3d++] = barriers.pjoint + (ip.njoints)*person + joint;
                    }
                }
                for(int joint = 0; joint < ip.njoints; joint++) {
                    if(jointInInfluence(j,joint)) {
                        edge[offset3d++] = barriers.sjoint + (ip.njoints)*i + joint;
                    }
                }
                edge[offset3d++] = barriers.absolutePose + (ip.njoints+1)*i;
                for(int joint = 0; joint < ip.njoints; joint++) {
                    if(jointInInfluence(j,joint)) {
                        edge[offset3d++] = barriers.absolutePose + (ip.njoints+1)*i + joint + 1;
                    }
                }
                for(int joint = 0; joint < ip.njoints; joint++) {
                    if(jointInInfluence(j,joint)) {
                        edge[offset3d++] = barriers.absolutePoseBase + (ip.njoints)*i + joint;
                    }
                }
                assert(edge.size() == edgeSize && edgeSize == offset3d);
                terms[tid].add(edge);
            }
            relativeScan++;
            if(relativeScan == corpora.scanCount[person]) {
                relativeScan = 0;
                person++;
            }
        }
    }

    void Problem::createAdvancedPoseTerm() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = ip.freeze_weights ? 3 : 4;
        std::vector<BarrierIndex> edge(edgeSize);
        for(int i = 0; i < ip.nscans; i++) {
            for(int j = 1; j < ip.njoints; j++) {
                int offset3d = 0;
                edge[offset3d++] = barriers.poseBaseMean + (j-1);
                if(!ip.freeze_weights) {
                    edge[offset3d++] = barriers.poseMean + (j-1);
                }
                edge[offset3d++] = barriers.poseBase + i*ip.njoints + j;
                edge[offset3d++] = barriers.pose + (ip.njoints+1)*i + j + 1;
                terms[tid].add(edge);
            }
        }
    }

    void Problem::createPJointTerm() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = PJOINT_SIZE;
        std::cout << "PJoint Term Edge Size: " << edgeSize << std::endl;
        std::vector<BarrierIndex> edge(edgeSize);
        int offset3d = 0;
        for(int i = 0; i < ip.npersons; i++) {
            for(int joint = 0; joint < ip.njoints; joint++) {
                offset3d = 0;
                edge[offset3d++] = barriers.pjoint + (ip.njoints)*i + joint;
                edge[offset3d++] = barriers.joints + joint;
                for(int k = 0; k < ip.ncomp; k++) {
                    edge[offset3d++] = barriers.joints + (1+k)*ip.njoints + joint;
                    edge[offset3d++] = barriers.coeff + ip.ncomp*i+k;
                }
                terms[tid].add(edge);
            }
        }
        std::cout << "End PJoint Term Edge Size: " << edgeSize << std::endl;
    }
    void Problem::createWeightPriorTerm() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = 2;
        std::vector<BarrierIndex> edge(edgeSize);
        for(int i = 0; i < ip.nverts; i++) {
            SVector4 influence = closestJoints.col(i).cast<int>();
            for(int k = 0; k < ip.njoints; k++) {
                bool dead = true;
                for(int j = 0; j < 4; j++) {
                    int joint = influence[j];
                    if(joint == k) dead = false;
                }
                if(dead) {
                    edge[0] = barriers.weights + i*ip.njoints + k;
                    edge[1] = barriers.weightInit + i*ip.njoints + k;
                    terms[tid].add(edge);
                }
            }
        }
    }

    void Problem::createWeightNormTerm() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = 4;
        std::vector<BarrierIndex> edge(edgeSize);
        for(int i = 0; i < ip.nverts; i++) {
            SVector4 influence = closestJoints.col(i).cast<int>();
            for(int j = 0; j < 4; j++) {
                int joint = influence[j];
                edge[j] = barriers.weights + i*ip.njoints + joint;
            }
            terms[tid].add(edge);
        }
    }
    void Problem::createMeanShapeTerm(MeshTopology* top) {
        int tid = termIdx++;
        std::vector<BarrierIndex> edge(4+2*4+4+1);
        InputParams& ip = inputParams;
        for(int i = 0; i < top->quads.cols(); i++) {
            auto quad = top->quads.col(i);
            for(int m = 0; m < 4; m++) {
                bool head = true;
                for(int j = 0; j < 4; j++) {
                    int u = quad[j];
                    head = head && rigidVertices[u];
                    edge[0+j] = barriers.shape + u;
                    edge[4+j] = barriers.meanShapePrior + u;
                    edge[8+j] = barriers.metric + 16*i + 4*m + j;
                    edge[12+j]= barriers.one;
                }
                edge[16] = head ? barriers.one : barriers.zero;
                terms[tid].add(edge);
            }

            for(int c = 0; c < ip.ncomp; c++) {
                for(int m = 0; m < 4; m++) {
                    bool head = true;
                    for(int j = 0; j < 4; j++) {
                        int u = quad[j];
                        head = head && rigidVertices[u]; //remove headVertices
                        edge[0+j] = barriers.shape + ip.nverts*(c+1) + u;
                        edge[4+j] = barriers.zero;
                        edge[8+j] = barriers.metric + 16*i + 4*m + j;
                        edge[12+j]= barriers.one;
                    }
                    edge[16] = head ? barriers.one : barriers.zero;
                    terms[tid].add(edge);
                }
            }
        }

        if(ip.ndeformjoints == 0) return;
        //for(int k = 0; k < ip.ndeformjoints*ip.ndeformshapes; k++) {
        //    int jointPos = k / 9;
        for(int i = 0; i < top->quads.cols(); i++) {
            auto quad = top->quads.col(i);
            std::vector<int> joints;
            for(int joint = 0; joint < ip.ndeformjoints; joint++) {
                for(int j = 0; j < 4; j++) {
                    int u = quad[j];
                    int jdiff = closestDeformJoints(joint,u);
                    if(std::find(joints.begin(), joints.end(), jdiff) == joints.end()) {
                        joints.push_back(jdiff);
                    }
                }
            }
            for(int k=0;k<9;k++) {
                for(auto joint : joints) {
                    for(int m = 0; m < 4; m++) {
                        bool head = true;
                        for(int j = 0; j < 4; j++) {
                            int u = quad[j];
                            int v = quad[(j+1) % 4];
                            //head = head && headVertices[u];
                            head = head && rigidVertices[u];

                            assert(ip.ndeformjoints == 2);
                            assert(ip.ndeformshapes == 9);

                            int joint_u1 = closestDeformJoints(0,u);
                            int joint_u2 = closestDeformJoints(1,u);
                            int kk = -1;
                            if(joint_u1 == joint) kk = 0;
                            if(joint_u2 == joint) kk = 1;
                            if(kk == -1) { /*kk == -1*/
                                edge[0+j] = barriers.deform + (9*0 + k)*ip.nverts + u; //fake
                                edge[4+j] = barriers.zero;
                                edge[8+j] = barriers.metric + 16*i + 4*m + j;
                                edge[12+j]= barriers.zero;
                            } else {
                                edge[0+j] = barriers.deform + (9*kk + k)*ip.nverts + u;
                                edge[4+j] = barriers.zero;
                                edge[8+j] = barriers.metric + 16*i + 4*m + j;
                                edge[12+j]= barriers.one;
                            }
                        }
                        edge[16] = head ? barriers.three : barriers.two;
                        terms[tid].add(edge);
                    }
                }
            }
        }
    }

    void Problem::evalMeanShapeTerm(const Model& model) {
        InputParams& ip = inputParams;
        double sum = 0;
        for(int i = 0; i < model.top.quads.cols(); i++) {
            auto quad = model.top.quads.col(i);
            for(int m = 0; m < 4; m++) {
                Matrix3X U; U.resize(3,4);
                Matrix3X V; V.resize(3,4);
                for(int j = 0; j < 4; j++) {
                    int u = quad[j];
                    int u2 = quad[(j+1)%4];
                    U.col(j) = -model.mean.col(u) + model.mean.col(u2);
                    V.col(j) = -control_vertices.col(u) + control_vertices.col(u2);
                    //Matrix4 M;
                    //Scalar m1 =
                    //M <<  = barriers.metric + 16*i + 4*m + j;
                }
                sum += .1*1.*(U.col(m).transpose()*V.col(m)).squaredNorm();
            }
/*
            for(int c = 0; c <= ip.ncomp; c++) {
                for(int m = 0; m < 4; m++) {
                    for(int j = 0; j < 4; j++) {
                        int u = quad[j];
                        edge[0+j] = barriers.shape + ip.nverts*(c+1) + u;
                        edge[4+j] = barriers.zero;
                        edge[8+j] = barriers.metric + 16*i + 4*m + j;
                    }
                }
            }
*/
        }
        std::cout << "Shape Regularization Error: " << sum << std::endl;
    }

    void Problem::createCoeffRegTerm() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = 2;
        std::vector<BarrierIndex> edge(edgeSize);
        for(int i = 0; i < ip.npersons; i++) {
            for(int k = 0; k < ip.ncomp; k++) {
                edge[0] = barriers.coeff + ip.ncomp*i+k;
                edge[1] = barriers.personCount + i;
				std::cout << "Person Count: " <<
				inputData.readPoint1d(barriers.personCount + i)
				<< std::endl;
                terms[tid].add(edge);
            }
        }
		//exit(1);
    }

	void Problem::createGroundRepelTerm(MeshTopology* top) {
		int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = 3+16;
        std::vector<BarrierIndex> edge(edgeSize);
        for(int i = 0; i < ip.nscans; i++) {
			int person = corpora.scan_to_person_id[i];
			//if (!corpora.person_fist[person])
			//	continue;
            for(int v = 0; v < ip.nverts; v++) {
                int offset3d = 0;
                double u1,u2;
                int fidx, pidx;
                fidx = pidx = -1;
                top->vertexToSurfaceCoordinate(v, fidx, u1, u2);
                Vector2 pt;
                pt[0] = u1;
                pt[1] = u2;
                eval->normalize(fidx, pt, &pidx); //is that  required?
                //set by reversed ICP
                //edge[offset3d++] = barriers.pointBase + 2*i*ip.npoints;
				edge[offset3d++] = barriers.vertexLabel + v;
				edge[offset3d++] = barriers.singleZero + corpora.person_ground_axes[person];
				edge[offset3d++] = barriers.groundPlanes + person;

                std::cout << top->num_faces() << "," << i << ","
                          << v << "," << fidx << "," << pidx << std::endl;

                for(int q = 0; q < 16; q++) {
                    int offset = bspline_control.at(16*pidx + q);
                    assert(offset < localPointCount);
                    edge[offset3d++] = barriers.localRegistration + i*localPointCount + offset;
                }
                terms[tid].add(edge);
            }
        }
    }


    void Problem::createJointRingReg() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = 1+2*ip.max_vertex_ring;
        std::vector<BarrierIndex> edge(edgeSize);
        for (int c = 0; c < ip.ncomp+1; c++) {
            for (size_t i = 0; i < ip.njoints-1; i++) {
                auto key = skeleton.getJoint(i + 1);
                auto secondKey = key;
                if (border_groups.find(key) == border_groups.end()) {
                    secondKey.first = key.second;
                    secondKey.second = key.first;
                }
                auto vertex_ring = border_groups.at(secondKey);
                std::cout << key.first << " " << vertex_ring.size()  << std::endl;
                edge[0] = barriers.joints + ip.njoints*c + i + 1;
                for (size_t j = 0; j < ip.max_vertex_ring; j++)  {
                    assert(vertex_ring.size() <= ip.max_vertex_ring);
                    if (j >= vertex_ring.size()) {
                        edge[1+2*j] = barriers.shape;
                        edge[2+2*j] = barriers.singleZero;
                    } else {
                        size_t vertex_index = vertex_ring[j];
                        edge[1+2*j] = barriers.shape + c*ip.nverts + vertex_index;
                        edge[2+2*j] = barriers.singleOne;
                    }
                }
                terms[tid].add(edge);

            }
        }
        //exit(0);

    }

    void Problem::createJointReg() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = 4;
        std::vector<BarrierIndex> edge(edgeSize);
        for (size_t i = 0; i < ip.njoints-1; i++) {
            auto key = skeleton.getJoint(i + 1);
            auto secondKey = key;
            if (border_groups.find(key) == border_groups.end()) {
                secondKey.first = key.second;
                secondKey.second = key.first;
            }
            auto vertex_ring = border_groups.at(secondKey);
            //std::cout << secondKey.first << ":" << secondKey.second << ":" << vertex_ring.size() <<std::endl;
            //auto lu = jointMap.at(secondKey);
            edge[0] = barriers.jointPrior + 1 + i;
            edge[1] = barriers.joints + i + 1;
            for (size_t j = 0; j < vertex_ring.size(); j++)  {
                size_t vertex_index = vertex_ring[j];
                std::cout << (i+1) << " -> " << vertex_index << std::endl;
                edge[2] = barriers.meanShapePrior + vertex_index;
                edge[3] = barriers.shape + vertex_index;
                terms[tid].add(edge);
            }

            for(int c = 0; c < ip.ncomp; c++) {
                edge[0] = barriers.zero;
                edge[1] = barriers.joints + ip.njoints*(c+1) + i + 1;
                for (size_t j = 0; j < vertex_ring.size(); j++) {
                    size_t vertex_index = vertex_ring[j];
                    edge[2] = barriers.zero;
                    edge[3] = barriers.shape + (c+1)*ip.nverts + vertex_index;
                    terms[tid].add(edge);
                }
            }

        }
    }

    void Problem::createPoseReg() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        int edgeSize = 2;
        std::vector<BarrierIndex> edge(edgeSize);
        for(int i = 0; i < ip.nscans; i++) {

            for(int joint = 1; joint < ip.njoints; joint++) {
                edge[0] = barriers.pose + (ip.njoints+1)*i + 1 +joint;
                edge[1] = barriers.poseBase + i*inputParams.njoints + joint;
                terms[tid].add(edge);
            }

            //hands!
            /*
            edge[0] = barriers.pose + (ip.njoints+1)*i + 1 +4;
            edge[1] = barriers.poseBase + i*inputParams.njoints + 4;
            terms[tid].add(edge);
            edge[0] = barriers.pose + (ip.njoints+1)*i + 1 +7;
            edge[1] = barriers.poseBase + i*inputParams.njoints + 7;
            terms[tid].add(edge);
            */
        }
    }

    void Problem::centerReg() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        std::vector<BarrierIndex> edge(1);
        //for(int i = 0; i < ip.nverts; i++) {
        //edge[0] = barriers.shape + 613;
        //terms[tid].add(edge);

        //edge[0] = barriers.joints;
        //terms[tid].add(edge);

        //terms[tid].add(edge);
        //}
        for(int c = 0; c < ip.ncomp+1; c++) {
            edge[0] = barriers.joints + c*ip.njoints;
            terms[tid].add(edge);

			edge[0] = barriers.shape + c*ip.nverts + ip.center_vertex_id; //613 259
            terms[tid].add(edge);
        }
    }

    void Problem::createShapeReg() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        std::vector<BarrierIndex> edge(1);
        for(int i = 0; i < ip.nverts; i++) {
            for(int k = 0; k < ip.ncomp; k++) {
                edge[0] = barriers.shape + ip.ncomp*(k+1)+i;
                terms[tid].add(edge);
            }
        }
    }

    void Problem::createLocalTerm() {
        int edgeSize = 2+LOCAL_SIZE;
        std::vector<BarrierIndex> edge(edgeSize);
        InputParams& ip = inputParams;
        int tid = termIdx++;
        for(int i = 0; i < ip.nscans; i++) {
            for(int j = 0; j < localPointCount; j++) {
                int offset3d = 0;
                edge[offset3d++] = barriers.localRegistration + i*localPointCount + j;
                edge[offset3d++] = barriers.localStencil + j;
                for(int k = 0; k < LOCAL_SIZE; k++) {
                    assert(localIdx[j*LOCAL_SIZE+k] >= 0 && localIdx[j*LOCAL_SIZE+k] < ip.nverts);
                    edge[offset3d++] = barriers.registration + i*ip.nverts + localIdx[j*LOCAL_SIZE+k];
                }
                terms[tid].add(edge);
            }
        }
    }

    void Problem::createTemporalPrior() {
        InputParams& ip = inputParams;
        int tid = termIdx++;
        //int edgeSize = 2*(ip.njoints+1+1);
        int edgeSize = 4;
        std::vector<BarrierIndex> edge(edgeSize);
        int seqId = 0;
        for(int i = 0; i < ip.nscans-1; i++) {
			/*
            if(corpora.sequences[seqId] == i) {
                seqId++;
                continue;
            } else {

            }*/

            std::cout << "Temporal prior for t = " << i << std::endl;
            int j = i + 1;

            for(int joint = 0; joint < ip.njoints; joint++) {
                //if(skeleton.isRigid(joint)) continue;
                int offset3d = 0;
                //edge[offset3d++] = barriers.pose + (ip.njoints+1)*i;
                edge[offset3d++] = barriers.pose + (ip.njoints+1)*i + 1 + joint;
                edge[offset3d++] = barriers.poseBase + (ip.njoints)*i + joint;
                //}
                //edge[offset3d++] = i;

                //for(int joint = 0; joint < ip.njoints+1; joint++) {
                //edge[offset3d++] = barriers.pose + (ip.njoints+1)*j;
                edge[offset3d++] = barriers.pose + (ip.njoints+1)*j + 1 + joint;
                edge[offset3d++] = barriers.poseBase + (ip.njoints)*j + joint;
                //}
                //edge[offset3d++] = j;
                terms[tid].add(edge);
            }
/*            
			int offset3d = 0;
            edge[offset3d++] = barriers.pose + (ip.njoints+1)*i;
            edge[offset3d++] = barriers.zero;
            edge[offset3d++] = barriers.pose + (ip.njoints+1)*j;
            edge[offset3d++] = barriers.zero;
            terms[tid].add(edge);
*/
        }
    }

    void Problem::createSymmetryTerm() {

        InputParams& ip = inputParams;
        int tid = termIdx++;
        int edgeSize = 2;
        std::vector<BarrierIndex> edge(edgeSize);

        Matrix3X dummyVerts;
        loadObj(ip.template_mesh_path, NULL, &dummyVerts,NULL,SCALE);

        std::vector<int> sym; sym.resize(dummyVerts.cols());
        findSymmetric(dummyVerts, &sym);
		for(int c = 0; c < ip.ncomp+1; c++) {
        for(int i = 0; i < sym.size(); i++) {
            if(sym[i] >= 0) {
                //std::cout << i << " <-> " << sym[i] << std::endl;
                edge[0] = barriers.shape + c*ip.nverts + i;
                edge[1] = barriers.shape + c*ip.nverts + sym[i];
                //terms[tid].add(edge);
            } else {
                edge[0] = barriers.shape + c*ip.nverts + i;
                edge[1] = barriers.shape + c*ip.nverts + i;
			}
			terms[tid].add(edge);
        }
		for(int i = 0; i < ip.njoints; i++) {
            edge[0] = barriers.joints + c*ip.njoints + i;
			int j = skeleton.getJointMirror(i);
            edge[1] = barriers.joints + c*ip.njoints + j;

			auto jI = skeleton.getJoint(i);
			auto jJ = skeleton.getJoint(j);
			std::cout << "(" << jI.first << "," << jI.second << ") <-> (" << 
				jJ.first << "," << jJ.second << ")" << std::endl;
            terms[tid].add(edge);
        }

		}
    }

    void Problem::relativeToAbsoluteTerm() {
        int tid = termIdx++;
        jointSPTermIdx = tid;
        InputParams& ip = inputParams;
        std::vector<BarrierIndex> edge(12);
        int relativeScan = 0;
        int person = 0;
        for(int i = 0; i < ip.nscans; i++) {
            std::cout << person <<  " -> " << i << std::endl;
            for(int j = 0; j < ip.njoints; j++) {
                bool root = j == 0;
                int k = root ? 0 : skeleton.getParentJoint(j);
				bool fist = corpora.person_fist[person];
                std::cout << person <<  " -> " << i
                	<< " -> " << k << " -> " << j 
					<< " " << skeleton.isRigid(j,fist) << std::endl;
                edge[0] = barriers.absolutePose + i*(ip.njoints+1) + k + 1;
                edge[1] = barriers.absolutePoseBase + i*(ip.njoints) + k;
                edge[2] = barriers.pose + i*(ip.njoints+1) + j + 1;
                edge[3] = barriers.poseBase + i*(ip.njoints) + j;
                edge[4] = barriers.absolutePose + i*(ip.njoints+1) + j + 1;
                edge[5] = barriers.absolutePoseBase + i*(ip.njoints) + j;

                edge[6] = barriers.sjoint + i*ip.njoints + k;
                edge[7] = barriers.sjoint + i*ip.njoints + j;

                edge[8] = barriers.pjoint + person*ip.njoints + j;
                edge[9] = barriers.pjoint + person*ip.njoints + k;

                edge[10] = root ? barriers.one : barriers.zero;
                edge[11] = skeleton.isRigid(j,fist) ? barriers.one : barriers.zero;
                terms[tid].add(edge);
            }
            relativeScan++;
            if(relativeScan == corpora.scanCount[person]) {
                relativeScan = 0;
                person++;
            }
        }
		//exit(1);
    }
    void Problem::linkVariableTerm() {
        int tid = termIdx++;
        InputParams& ip = inputParams;
        std::vector<BarrierIndex> edge(2);
        for(int i = 0; i < ip.nscans; i++) {
            edge[0] = barriers.absolutePose + i*(ip.njoints+1);
            edge[1] = barriers.pose + i*(ip.njoints+1);
            terms[tid].add(edge);
        }
    }


    void Problem::initOpt() {
        InputParams& ip = inputParams;

        int   nLinearIterations = 200; //10000
        int   nNonLinearIterations = 10000;
        float funtol = 1E-6;
        float tol = 0.001;
        unsigned int dims[5];
        inputData.getDims(dims);
        plan = Opt_ProblemPlan(state, problem, dims); //&dim
        std::cout << "Problem plan exited" << std::endl;
        Opt_SetSolverParameter(state, plan, "nIterations", (void*)&nNonLinearIterations);
        Opt_SetSolverParameter(state, plan, "lIterations", (void*)&nLinearIterations);
        Opt_SetSolverParameter(state, plan, "q_tolerance", (void*)&tol);
		for(int i = 0; i < 5; i++) std::cout << dims[i] << std::endl;
		//exit(0);
    }

    void Problem::prepare() {
        inputData.allocate();
        for(int i = 0; i < terms.size(); i++) {
            terms[i].allocate();
        }

        copy();
    }

	void Problem::copy() {
        inputData.copy();
        for(int i = 0; i < terms.size(); i++) {
            terms[i].copy();
        }
    }

    void Problem::createOptParams() {
        const int inputCount = 8+1;

        int totalSize = inputCount;
        for(int i = 0; i < terms.size(); i++) {
            totalSize += terms[i].countIdxArrays() + 1;
            std::cout << "Term " << i <<
                      " with " << (terms[i].countIdxArrays() + 1)
                      << " params" << std::endl;
        }
        optinput.resize(totalSize);

        inputData.assign(&optinput[0]);
        optinput[7] = &w_fitSqrt;
        optinput[8] = &w_surface;
        int offset = inputCount;
        for(int i = 0; i < terms.size(); i++) {
            terms[i].assign(&optinput[offset]);
            offset += terms[i].countIdxArrays() + 1;
        }
        std::cout << "Total Size is: " << totalSize << std::endl;
    }

    void Problem::retrieve(Model& model) {
        InputParams& ip = inputParams;
        inputData.retrieve();
		std::cout << "Start Retrieve Model" << std::endl;
        model.closestJoints = closestJoints;
        model.closestDeformJoints = closestDeformJoints;
        model.skeleton = &skeleton;

        model.mean.resize(3,ip.nverts);
        inputData.readDynamic(model.mean,barriers.shape, ip.nverts);

        model.pcs.resize(inputParams.ncomp);
        model.jointpcs.resize(inputParams.ncomp);

        model.meanPose.resize(3,inputParams.njoints-1);

        for(int j = 0; j < inputParams.njoints-1; j++) {
            Scalar3 meanBase =
                    inputData.readPoint3d(barriers.poseBaseMean+j);
            Matrix3X X; X.resize(3,1);
            X << meanBase.x, meanBase.y,meanBase.z;

            Matrix3 U;
            Matrix3X u; u.resize(3,1); u.setConstant(0);
            if(!ip.freeze_weights) {
                inputData.read(u,barriers.poseMean+j,1);
            }

            angle_axis_to_rotation_matrix(u.col(0),&U);
            std::vector<Matrix3> bases;
            rod2Mat(X,bases);
            Matrix3 result = bases[0]*U;
            Vector3 euler = result.eulerAngles(2, 1, 0);
            //std::cout << j << " " << euler.transpose() << std::endl;
            model.meanPose.col(j) = euler;
        }

        for(int i = 0; i < ip.ncomp; i++) {
            model.pcs[i].resize(3,ip.nverts);
            inputData.readDynamic(model.pcs[i],barriers.shape + (1+i)*ip.nverts,ip.nverts);
        }

        for(int i = 0; i < ip.npersons; i++) {
            Person person;
			person.coeff.resize(ip.ncomp);
            for(int j = 0; j < ip.ncomp; j++) {
                Scalar coeff = inputData.read(barriers.coeff + i*ip.ncomp+j);
                person.coeff[j] = coeff;
            }
            person.joints.resize(3,ip.njoints);
            inputData.read(person.joints,barriers.pjoint + i*inputParams.njoints,inputParams.njoints);
            model.persons.push_back(person);
        }

        for(int i = 0; i < ip.ncomp+1; i++) {

            if(i == 0) {
                model.joints.resize(3,ip.njoints);
                inputData.readDynamic(model.joints, barriers.joints, inputParams.njoints);
            } else {
                model.jointpcs[i-1].resize(3,ip.njoints);
                inputData.readDynamic(model.jointpcs[i-1],barriers.joints + i*inputParams.njoints,inputParams.njoints);
            }

        }

        //barriers.pose = barriers.joints+(inputParams.ncomp+1)*inputParams.njoints;
        for(int i = 0; i < inputParams.nscans; i++) {
            Scan scan;
            scan.trans.resize(3,1);
            scan.pose.resize(3,inputParams.njoints);
            scan.posebase.resize(3,inputParams.njoints);
            scan.absoluteTrans.resize(3,1);
            scan.absolutePose.resize(3,inputParams.njoints);
            scan.absolutePosebase.resize(3,inputParams.njoints);
            scan.joints.resize(3,inputParams.njoints);

            inputData.read(scan.trans,barriers.pose+(inputParams.njoints+1)*i,1);
            inputData.read(scan.absoluteTrans,barriers.absolutePose+(inputParams.njoints+1)*i,1);
            inputData.read(scan.joints,barriers.sjoint+i*inputParams.njoints,inputParams.njoints);

			scan.T.resize(3,inputParams.nverts);
            inputData.read(scan.T,barriers.registration+(inputParams.nverts)*i,inputParams.nverts);

            scan.label2dRays.resize(ip.nviews);

			int current_person = corpora.scan_to_person_id[i];

            for(int v = 0; v < ip.nviews; v++) {
                scan.label2dRays[v].resize(3,ip.njointlabels+ip.nverts+1);
				if(!corpora.person_has_labels_2d[current_person]) {
					scan.label2dRays[v].setConstant(0);
					continue;
				}
				int camId = corpora.person_fist[current_person] ? 1 : 0; 
				std::string camPath = ip.camPaths[camId];
                Camera camera(camPath, v);
                {
                    scan.label2dRays[v].col(ip.njointlabels+ip.nverts) = camera.getCameraPosition();
                }
                for(int j = 0; j < ip.njointlabels; j++) {
                    Scalar2 a = inputData.readPoint2d(barriers.label2d
                                                      +v*ip.nscans*ip.njointlabels+ip.njointlabels*i + j);
                    Vector2 pt;
                    pt[0] = a.x;
                    pt[1] = a.y;
                    scan.label2dRays[v].col(j) = camera.computeRay(pt);
                }
				for(int j = 0; j < ip.nverts; j++) {
                    Scalar2 a = inputData.readPoint2d(barriers.labelSurface2d
                                                      +v*ip.nscans*ip.nverts+ip.nverts*i + j);
                    Vector2 pt;
                    pt[0] = a.x;
                    pt[1] = a.y;
                    scan.label2dRays[v].col(ip.njointlabels+j) = camera.computeRay(pt);
                }

            }

            //inputData.read(scan.pose,barriers.pose+(inputParams.njoints+1)*i+1,inputParams.njoints);


            for(int j = 0; j < inputParams.njoints; j++) {
                Scalar3 thetabase  = inputData.readPoint3d(barriers.poseBase + i*inputParams.njoints + j);
                scan.posebase(0,j) = thetabase.x;
                scan.posebase(1,j) = thetabase.y;
                scan.posebase(2,j) = thetabase.z;

				Matrix3X theta; theta.resize(3,1);
                inputData.read(theta,barriers.pose + i*(inputParams.njoints+1) + 1 + j,1);
                scan.pose(0,j) = theta(0,0);
                scan.pose(1,j) = theta(1,0);
                scan.pose(2,j) = theta(2,0);

            }
            for(int j = 0; j < inputParams.njoints; j++) {
				Scalar3 thetabase  = inputData.readPoint3d(barriers.absolutePoseBase + i*inputParams.njoints + j);
                scan.absolutePosebase(0,j) = thetabase.x;
                scan.absolutePosebase(1,j) = thetabase.y;
                scan.absolutePosebase(2,j) = thetabase.z;

                Matrix3X theta; theta.resize(3,1);
                inputData.read(theta,barriers.absolutePose + i*(inputParams.njoints+1) + 1 + j,1);
                scan.absolutePose(0,j) = theta(0,0);
                scan.absolutePose(1,j) = theta(1,0);
                scan.absolutePose(2,j) = theta(2,0);
            }
            //std::cout << scan.pose.transpose() << std::endl;
            scan.Trest.resize(3,inputParams.nverts);
            inputData.read(scan.Trest,barriers.pshape+(inputParams.nverts)*i,inputParams.nverts);

            //scan.T.resize(3,inputParams.nverts);
            //inputData.read(scan.T,barriers.registration+(inputParams.nverts)*i,inputParams.nverts);

            scan.Tlocal.resize(3,localPointCount);
            inputData.read(scan.Tlocal,barriers.localRegistration+localPointCount*i,localPointCount);

            scan.points.resize(3, inputParams.npoints);
			scan.patchPoints.resize(3, inputParams.npoints);

            std::vector<SurfacePoint> surfacePoints;
            for(int p = 0; p < ip.npoints; p++) {
                //transform to unnormalized space
                Scalar2 ur = inputData.readPoint2d(
                        barriers.surface + i*ip.npoints+p);
                Vector2 us; us[0] = ur.x; us[1] = ur.y;
                const int patch = patch_idx[i*ip.npoints+p];
                us = eval->unnormalize(patch, us);

                SurfacePoint surfacePoint;
                surfacePoint.u = us;
                surfacePoint.face = face_idx[i*ip.npoints+p];
                surfacePoints.push_back(surfacePoint);
            }
	
            //u = u * 100;

			Matrix3X Xdu; Xdu.resize(3,ip.npoints);
			Matrix3X Xdv; Xdv.resize(3,ip.npoints);

            //eval->evaluateSubdivSurface(scan.T,
            //    surfacePoints, &scan.points,NULL,NULL,NULL,
			//	&Xdu, &Xdv);
			eval->evalPointsPartial(scan.Tlocal, surfacePoints,&scan.patchPoints,&Xdu,&Xdv);
			Matrix3X tangentPoint; tangentPoint.resize(3,ip.npoints);

			Matrix2X deltaU; deltaU.resize(2,ip.npoints);
            inputData.read(deltaU, barriers.surfaceUnknown + i*ip.npoints, ip.npoints);

			for(int p = 0; p < ip.npoints; p++) {
				Vector2 us;
                {
                    Scalar2 uS = inputData.readPoint2d(barriers.surfaceScale + i*ip.npoints+p);
                    us[0] = uS.x;
                    us[1] = uS.y;
                }

				tangentPoint.col(p) = scan.patchPoints.col(p) + 
					Xdu.col(p)*deltaU(0,p)*ip.continuous_step_size*us[0] + 
					Xdv.col(p)*deltaU(1,p)*ip.continuous_step_size*us[1];
            }
			scan.tangentStep = tangentPoint;

			scan.labels.resize(3,ip.nlabels); scan.labels.setZero();
            std::vector<SurfacePoint> labelPoints;

            auto markerIt = marker_vertex_groups.begin();
            for(int p = 0; p < ip.nlabels; p++) {

                double u,v;
                int fidx;
                int pidx;
                model.top.vertexToSurfaceCoordinate(markerIt->second[0], fidx, u, v);
                Vector2 pt;
                pt[0] = u;
                pt[1] = v;

                SurfacePoint surfacePoint;
                surfacePoint.u = pt;
                surfacePoint.face = fidx;
                labelPoints.push_back(surfacePoint);

                //std::cout << std::endl;
                ++markerIt;
            }
            eval->evaluateSubdivSurface(scan.T,
                                        labelPoints, &scan.labels);

            scan.labels_initial.resize(3,ip.nlabels);
            inputData.read(scan.labels_initial,barriers.debugRegistration+ip.nlabels*i,ip.nlabels);

            scan.labels_partial.resize(3,ip.nlabels);
            //scan.labels_partial = eval->evalPointsPartial(scan.Tlocal, labelPoints);
            //scan.labels_partial = eval->evalPointsCustom(scan.T, labelPoints);
			
			eval->evalPointsPartial(scan.Tlocal,labelPoints,&scan.labels_partial,NULL,NULL);

            scan.robust.resize(ip.npoints);
            inputData.read(scan.robust, barriers.robust+ip.npoints*i,ip.npoints);

            model.scans.push_back(scan);
        }

        static int bla = 0;
        model.weights.resize(ip.njoints,ip.nverts); model.weights.setZero();
        for(int i = 0; i < ip.nverts; i++) {
            SVector4 influence = closestJoints.col(i).cast<int>();

            Scalar denom = 0;
            for(int j = 0; j < 4; j++) {
                int joint = influence[j];
                Scalar w = inputData.readDynamic1d(barriers.weights+inputParams.njoints*i + joint) / LAMBDA;
                //model.weights(j,i) = model.weights(j,i)*model.weights(j,i);
                //denom += exp(w);
                denom += .5 + .5*w / sqrt(1 + w*w);
            }
            if(bla % 100 == 0) {
                //std::cout << i << "\t" << denom << std::endl;
            }
            for(int j = 0; j < 4; j++) {
                int joint = influence[j];
                Scalar w = inputData.readDynamic1d(barriers.weights+inputParams.njoints*i + joint) / LAMBDA;
                //Scalar wp = inputData.read(barriers.weights+inputParams.njoints*i + joint) / LAMBDA;
                //model.weights(j,i) = model.weights(j,i)*model.weights(j,i);
                //model.weights(joint,i) = 1.0 / (1 + exp(-LAMBDA*w));
                model.weights(joint,i) = .5 + .5*w / sqrt(1 + w*w);
                //model.weights(joint,i) = model.weights(joint,i) / denom;
            }
        }


        MatrixX raw; raw.resize(ip.njoints,ip.nverts);
        for(int i = 0; i < ip.nverts; i++) {
            for(int j = 0; j < ip.njoints; j++) {
                raw(j,i) = inputData.readDynamic1d(barriers.weights+inputParams.njoints*i + j);
            }
        }
        if(bla % 50 == 0) {
            //std::cout << raw.transpose() << std::endl;
            //std::cout << model.weights.transpose() << std::endl;
        }
        bla++;

        model.deform.resize(ip.ndeformjoints*ip.ndeformshapes);
        for(int i = 0; i < ip.ndeformjoints*ip.ndeformshapes; i++) {
            model.deform[i].resize(3,ip.nverts);
            inputData.readDynamic(model.deform[i],barriers.deform+ip.nverts*i,ip.nverts);
        }



        int personId = 0;
        int relativeScan = 0;

        Matrix3X mymean = model.mean;
        for(int i = 0; i < ip.nscans; i++) {
            Matrix3X C1 = model.mean.replicate(1,1); //C1.setConstant(0);
            Matrix3X skinned(3,ip.nverts);
            const Person& person = model.persons[personId];
            Matrix3X J1;
            J1 = model.joints;
            for(int j = 0; j < ip.ncomp; j++) {
                Scalar coeff = person.coeff[j];
                const Matrix3X& pc = model.pcs[j];
                C1 = C1 + coeff*pc;
                const Matrix3X& jpc = model.jointpcs[j];
                //std::cout << jpc.transpose() << std::endl;
                J1 = J1 + coeff*jpc;
            }
            std::vector<Vector3> skinjoints;
            std::vector<Matrix3> localrot;
            rod2Mat(model.scans[i].posebase,localrot);
            for(int j = 0; j < ip.njoints; j++) {
                Matrix3 so3; so3.setIdentity();
                Vector3 r = model.scans[i].pose.col(j);
                so3(0,1) = -r[2];
                so3(1,0) = r[2];
                so3(0,2) = r[1];
                so3(2,0) = -r[1];
                so3(1,2) = -r[0];
                so3(2,1) = r[0];
                localrot[j] = localrot[j]*so3;
            }
            unravel(J1,skinjoints);

            for(int j = 0; j < ip.ndeformjoints; j++) {
                for(int k = 0; k < ip.ndeformshapes; k++) {
                    for(int v = 0; v < ip.nverts; v++) {
                        int joint = closestDeformJoints(j,v);
                        Scalar melem = localrot[joint](k/3,k%3);
                        if(k % 4 == 0) melem -= 1;
                        C1.col(v) = C1.col(v) + melem*model.deform[9*j+k].col(v);
                    }
                }
            }

            std::vector<Matrix4> absolutes(ip.njoints);

            skeleton.calculate_kinematic_chain(skinjoints,
                                               localrot, absolutes);

            for(int v = 0; v < ip.nverts; v++) {
                VectorX wv = model.weights.col(v); wv.setZero();
                SVector4 influence = closestJoints.col(v).cast<int>();
                for(int infu =   0; infu < 4; infu++) {
                    int joint = influence[infu];
                    wv[joint] = model.weights(joint,v);
                }
                Vector3 vert = skeleton.skin_vertex(absolutes,wv,C1.col(v));
                skinned.col(v) = vert + model.scans[i].trans;
            }
            model.scans[i].T_estimate = skinned;

            model.scans[i].Tlocal_estimate.resize(3,localPointCount);
            OpenSubdiv::Far::StencilTable const * localStencilTable = eval->getLocalStencilTable();
            for(int k = 0; k < localPointCount; k++) {
                auto st = localStencilTable->GetStencil(k);
                const OpenSubdiv::Vtr::Index *ind = st.GetVertexIndices();
                float const *wei = st.GetWeights();
                Vector3 localPt; localPt.setZero();
                for(int j = 0; j < st.GetSize(); j++) {
                    localPt += wei[j]* model.scans[i].T_estimate.col(ind[j]);
                    //std::cout << ind[j] << ": " << wei[j] << ", ";
                }
                //std::cout << std::endl;
                model.scans[i].Tlocal_estimate.col(k) = localPt;
            }

            Matrix3X customEval = eval->evalPointsCustom(model.scans[i].T);
            Matrix3X partialEval = eval->evalPointsPartial(model.scans[i].Tlocal);

            relativeScan++;
            if(relativeScan == corpora.scanCount[personId]) {
                relativeScan = 0;
                personId++;
            }
        }

		std::cout << "End Retrieve Model" << std::endl;
    }

    void Problem::release() {
        inputData.release();
        for(int i = 0; i < terms.size(); i++) {
            terms[i].release();
        }
    }

    void Problem::manifold_update() {
        InputParams& ip = inputParams;
        inputData.retrieve();
		std::cout << "Start Manifold Update" << std::endl;

        for(int i = 0; i < 2*inputParams.nscans; i++) {
            //std::cout << "Result" << std::endl;

            //inputData.clamp2d(i*ip.npoints,ip.npoints);

            for(int j = 0; j < inputParams.njoints; j++) { //+trans
                BarrierIndex varOffset = barriers.pose + (ip.njoints+1)*i + j+1;
                //std::cout << "Offset is: " << varOffset << std::endl;
                Matrix3X u; u.resize(3,1);
                Matrix3 U;
                inputData.read(u,varOffset,1);
                //std::cout << j << " " << u.transpose() << std::endl;

                //std::cout << i << ", " << j << " -> " << X.transpose() <<" * " << u.transpose()  << std::endl;
                angle_axis_to_rotation_matrix(u.col(0),&U);
                //std::cout << "\n" << U << "\n" << std::endl;

				BarrierIndex baseOffset = barriers.poseBase + (ip.njoints)*i + j;
				inputData.readPoints(u,baseOffset,1);
                std::vector<Matrix3> bases;
                rod2Mat(u,bases);
                //angle_axis_to_rotation_matrix
                //std::cout << j << " -> " <<  bases[0] << std::endl;
                Matrix3 result = bases[0]*U;
                Vector3 euler = result.eulerAngles(2, 1, 0);
                //Scalar3 toset;
                //std::cout << j << " " << euler.transpose() << std::endl;
                Scalar3 out;
                out.x = euler(0);
                out.y = euler(1);
                out.z = euler(2);
                inputData.setPoint(barriers.poseBase + i*ip.njoints + j,out);
                inputData.resetUnknown3d(varOffset);
            }
        }
        for(int j = 0; j < inputParams.njoints-1; j++) {
            Scalar3 meanBase =
                    inputData.readPoint3d(barriers.poseBaseMean+j);
            Matrix3X X; X.resize(3,1);
            X << meanBase.x, meanBase.y,meanBase.z;

            Matrix3 U;
            Matrix3X u; u.resize(3,1); u.setConstant(0);
            if(!ip.freeze_weights) {
                inputData.read(u,barriers.poseMean+j,1);
            }

            angle_axis_to_rotation_matrix(u.col(0),&U);
            std::vector<Matrix3> bases;
            rod2Mat(X,bases);
            Matrix3 result = bases[0]*U;
            Vector3 euler = result.eulerAngles(2, 1, 0);
            //std::cout << j << " " << euler.transpose() << std::endl;
            Scalar3 out;
            out.x = euler(0);
            out.y = euler(1);
            out.z = euler(2);
            inputData.setPoint(barriers.poseBaseMean+j,out);
            if(!ip.freeze_weights) {
                inputData.resetUnknown3d(barriers.poseMean+j);
            }
        }

        //std::cout << "End" << std::endl;
        //}
        copy();
		std::cout << "End Manifold Update" << std::endl;
    }

#define TIC t1 = std::chrono::high_resolution_clock::now();
#define TOC(NAME) t2 = std::chrono::high_resolution_clock::now(); std::cout << (NAME) << " " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;

//#define TIC ;
//#define TOC(NAME) ;
    void Problem::removeRigidJoints() {
        //TODO retrieve and save
        inputData.retrieve();
        InputParams& ip = inputParams;
        //std::vector<std::string> jointsToRemove =
        //        {"m_rhand","m_lhand","m_rfoot","m_lfoot"};
		std::vector<std::string> jointsToRemove =
        {
	            "m_neck",
                "m_head",
                "m_rlowerarm",
                "m_rhand", 
/*
                "m_rthumb",
                "m_rthumb-proximal",
                "m_rthumb-distal",
                "m_rfist",
                "m_rindex-middle",
                "m_rmiddle-middle",
                "m_rring-middle",
                "m_rpinky-middle",
*/                
                "m_llowerarm",
                "m_lhand", 
/*
               "m_lthumb",
               "m_lthumb-proximal",
               "m_lthumb-distal",
               "m_lfist",
               "m_lindex-middle",
               "m_lmiddle-middle",
               "m_lring-middle",
               "m_lpinky-middle",
*/
                "m_hip",
                "m_pelvis",
                "m_rupperleg",
                "m_rlowerleg",
                "m_rfoot",
                "m_lupperleg",
                "m_llowerleg",
                "m_lfoot"

		};
        int restJointColBase3 = 2+ip.ncomp+ip.ndeformjoints*ip.ndeformshapes+2;
        for(int i = 0; i < ip.nscans; i++) {
            for(int j = 0; j < ip.nverts; j++) {
                for(int k = 0; k < ip.ndeformjoints; k++) {
                    int jointIdx = closestDeformJoints(k,j);
                    int restJointCol3 = restJointColBase3 + 3*k;
                    int restJointRow3= i*ip.nverts+j;
                    std::string joint = skeleton.getJoint(jointIdx).first;
                    std::vector<std::string>::iterator it = std::find(jointsToRemove.begin(), jointsToRemove.end(), joint);
                    if(it == jointsToRemove.end())
                        continue;
                    //assert(jointIdx == *it);
					
                    terms[modelTermIndex].set(restJointCol3, restJointRow3,
                                 barriers.one);

                }
            }
            for(int k = 0; k < ip.njoints; k++) {
                terms[jointSPTermIdx].set(11,i*ip.njoints+k,barriers.zero);
            }
        }
        copy();
    }

    void Problem::reverse_icp_update(MeshTopology* top) {
        InputParams& ip = inputParams;
        for(int i = 0; i < ip.nscans; i++) {
            Matrix3X T; T.resize(3,ip.nverts);
            inputData.read(T,barriers.registration+(inputParams.nverts)*i,
                           inputParams.nverts);
            Matrix3X modelpts; modelpts.resize(3, 1);
            Icp icp(cloud[i]);

            for(int v = 0; v < ip.nverts; v++) {
                double u1,u2;
                int fidx, pidx;
                std::vector<SurfacePoint> surfacePoints;
                SurfacePoint sp;

                top->vertexToSurfaceCoordinate(v, fidx, u1, u2);
                Vector2 pt;
                pt[0] = u1;
                pt[1] = u2;
                sp.u = pt;
                sp.face = fidx;
                surfacePoints.push_back(sp);
                //eval->normalize(fidx, pt, &pidx); Normalize Required?
                eval->evaluateSubdivSurface(T,surfacePoints,&modelpts);
                int idx = icp.findClosestPoint(modelpts.col(0));
                model2Data[i][v] = idx;
                terms[1].set(0, reverse_offset + i*ip.nverts+v,
                             barriers.pointBase + 2*i*ip.npoints + 2*idx+0);
            }
        }
    }


    void Problem::icp_update(MeshTopology* top, bool alternating) {
        //return;
        InputParams& ip = inputParams;
        std::chrono::high_resolution_clock::time_point t1, t2;
        //TIC
        inputData.retrieve();
        //TOC("Retrieve memory")

        std::vector<SurfacePoint> surfacePoints; surfacePoints.resize(ip.npoints);
        SubdivEvaluator::triplets_t dSdX, dSudX, dSvdX;
        Matrix3X T; T.resize(3,ip.nverts);
        Matrix3X localPts(3,localPointCount);

        Matrix3X surface3d; surface3d.resize(3,ip.npoints);
        Matrix3X dSdu; dSdu.resize(3,ip.npoints);
        Matrix3X dSdv; dSdv.resize(3,ip.npoints);
        Matrix3X dSduu; dSduu.resize(3,ip.npoints);
        Matrix3X dSduv; dSduv.resize(3,ip.npoints);
        Matrix3X dSdvv; dSdvv.resize(3,ip.npoints);

        Matrix3X c1; c1.resize(3, ip.npoints);
        Matrix3X c2; c2.resize(3, ip.npoints);

        for(int i = 0; i < ip.nscans; i++) {
            //TIC
            //std::cout << "Start icp update" << std::endl;
            BarrierIndex varOffset = barriers.registration + i*ip.nverts;
            inputData.read(T,varOffset,ip.nverts);
            inputData.read(localPts,barriers.localRegistration + i*localPointCount,localPointCount);
            surfaceState.setCurrentModel(i,T, localPts, true);
            //TOC(std::string("DISCRETE UPDATE SET MODEL ")+std::to_string(i))
            //TIC
            //std::cout << "current model set" << std::endl;
            int dicreteCount = 0;

            for(int p = 0; p < ip.npoints; p++) {
                Vector2 uold;
                {
                    Scalar2 ur = inputData.readPoint2d(
                            barriers.surface + i*ip.npoints+p);
                    uold[0] = ur.x;// = 0;
                    uold[1] = ur.y;// = 0;
                }
                int fold = face_idx[i*ip.npoints+p];
                int patch = patch_idx[i*ip.npoints+p];

                uold = eval->unnormalize(patch, uold);
                //f = 0;
                //std::cout << f << "\t" << u.transpose() << "\t";
                Vector2 u;
                int f;
                bool res;
                u = uold;
                f = fold;
                Vector2 d; d.setConstant(0);
                Vector2 unew; int fnew;

                res = surfaceState.discreteUpdate(i,p, u, f);

                if(0) {
                    surfaceState.alternatingUpdate(i, p, u, f, d, unew,fnew);
                    u = unew;
                    f = fnew;
                }


                if(res)  {
                    dicreteCount++;
                }

                //std::cout << f << "\t" << u.transpose() << std::endl;
                auto quad = top->quads.col(f);
                face_idx[i*ip.npoints+p] = f;
                {
                    Scalar2 ur;
                    u = eval->normalize(f, u, &patch_idx[i*ip.npoints+p]);
                    ur.x = u[0];
                    ur.y = u[1];
                    inputData.setPoint(
                            barriers.surface + i*ip.npoints+p, ur);
                }

                surfacePoints[p].u[0] = u[0];
                surfacePoints[p].u[1] = u[1];
                surfacePoints[p].face = face_idx[i*ip.npoints+p];
            }
            //TOC(std::string("DISCRETE UPDATE SUBDIV ")+std::to_string(i))
            //std::cout << dicreteCount << " : " << ip.npoints << std::endl;



            //TIC
            //eval->evaluateSubdivSurface(T,
            //    surfacePoints, &surface3d, &dSdX, &dSudX, &dSvdX,
            //    &dSdu, &dSdv, &dSduu, &dSduv, &dSdvv);
			Matrix3X X; X.resize(3, ip.npoints);
            Matrix3X Xdu; Xdu.resize(3, ip.npoints);
            Matrix3X Xdv; Xdv.resize(3, ip.npoints);
            eval->evalPointsPartial(localPts, surfacePoints,&X,&Xdu,&Xdv);

            int dataId = 0;
            for(int p = 0; p < ip.npoints; p++) {
				Scalar2 us;
                us.x = 1.0; // / Xdu.col(p).stableNorm();
                us.y = 1.0; // / Xdv.col(p).stableNorm();
                inputData.setPoint(barriers.surfaceScale + i*ip.npoints+p, us);


                int pidx = patch_idx[i*ip.npoints+p];
                //std::cout << p << "\t" << quadId << "\t" << dSdX.size() << std::endl;
                for(int j = 0; j < 16; j++) {
                    //Matrix3X T; T.resize(3,ip.nverts);
                    //Matrix3X T_old; T_old.resize(3,ip.nverts);
                    //inputData.read(T,varOffset,ip.nverts);
                    //inputData.read(T_old,varOffset,ip.nverts,true);

                    int offset = bspline_control.at(16*pidx + j);
                    assert(offset < localPointCount);
                    BarrierIndex value = barriers.localRegistration + i*localPointCount + offset;
                    //terms[0].set(8+j, i*ip.npoints+p, value);
                    terms[0].set(6+j, i*ip.npoints+p, value); //4
                }

            }
            //std::cout << i << "," << std::flush;
            //TOC(std::string("DISCRETE UPDATE SET VALUES ")+std::to_string(i))
        }
        if(ip.useModel2Data) {
            reverse_icp_update(top);
        }
		//std::cout << std::endl;
        //TIC
        copy();
        //TOC("COPY MEMORY")
    }


    void Problem::continuous_update(int iteration, int failCount, MeshTopology* top, bool fake) {
        InputParams& ip = inputParams;
        inputData.retrieve();
        double error = 0;
        double modelError = 0;
        double oldError = 0;
        for(int i = 0; i < ip.nscans; i++) {
            //std::cout << "Start continuous update" << std::endl;
            std::vector<SurfacePoint> surfacePoints;
            std::vector<SurfacePoint> oldPoints;

            SubdivEvaluator::triplets_t dSdX, dSudX, dSvdX;
            BarrierIndex varOffset = barriers.registration + i*ip.nverts;
            Matrix3X T; T.resize(3,ip.nverts);
            //Matrix3X T_old; T_old.resize(3,ip.nverts);
            Matrix3X localPts(3,localPointCount);
            inputData.read(T,varOffset,ip.nverts, false); //true
            //inputData.read(T_old,varOffset,ip.nverts,false);
            inputData.read(localPts,barriers.localRegistration + i*localPointCount,localPointCount, false); //true
            //surfaceState.setCurrentModel(i,T_old, false);
            surfaceState.setCurrentModel(i,T,localPts, false); //false
            //std::cout << "current model set" << std::endl;

            //Matrix3X tangent; tangent.resize(3,ip.npoints);

            Matrix2X deltaU; deltaU.resize(2,ip.npoints);
            //deltaU.setZero(); //REMOVE!!!
            inputData.read(deltaU, barriers.surfaceUnknown + i*ip.npoints, ip.npoints);
            deltaU = ip.continuous_step_size*deltaU;
            //deltaU.setConstant(1);

            //if(fake) deltaU = 0*deltaU;
            //std::cout << "delU " << deltaU.transpose() << std::endl;
            Scalar maxMetric = 0;
            Scalar respLocal = 0;
            Scalar maxLocal = 0;
            Scalar repMetric = 0;
            Scalar meanMetric = 0;
            Scalar meanLocal = 0;
            int faceJump = 0;

            double planarUpdate = 0;


            for(int p = 0; p < ip.npoints; p++) {
                inputData.resetUnknown2d(
                        barriers.surfaceUnknown + i*ip.npoints+p);
                Vector2 u;
                {
                    Scalar2 ur = inputData.readPoint2d(barriers.surface + i*ip.npoints+p);
                    u[0] = ur.x;
                    u[1] = ur.y;
                    SurfacePoint surfacePoint;
                    surfacePoint.u[0] = u[0];
                    surfacePoint.u[1] = u[1];
                    surfacePoint.face = face_idx[i*ip.npoints+p];
                    oldPoints.push_back(surfacePoint);
                }
                int f = face_idx[i*ip.npoints+p];
                int fnew;
                Vector2 unew;
                Scalar singleDistance;

                maxLocal = fmax(deltaU.col(p).norm(), maxLocal);
				
				Vector2 us;
				{
					Scalar2 uS = inputData.readPoint2d(barriers.surfaceScale + i*ip.npoints+p);
					us[0] = uS.x;
					us[1] = uS.y;
				}

				Vector2 step = deltaU.col(p);
				step[0] *= us[0];
				step[1] *= us[1];
                int jumpCount = surfaceState.move(i, u, f, step, unew, fnew, &singleDistance);
                if(f != fnew) faceJump++;
                f = fnew;
                u = unew;

                auto quad = top->quads.col(f);
                //face_idx[i*ip.npoints+p] = f;
                /*{
                    Scalar2 ur;
                    ur.x = u[0];
                    ur.y = u[1];
                    inputData.setPoint(i*ip.npoints+p, ur);
                }*/

                SurfacePoint surfacePoint;
                surfacePoint.u[0] = u[0];
                surfacePoint.u[1] = u[1];
                surfacePoint.face = f;//face_idx[i*ip.npoints+p];
                surfacePoints.push_back(surfacePoint);
            }

            Matrix3X X; X.resize(3, ip.npoints);
            Matrix3X Y; Y.resize(3, ip.npoints);
			Matrix3X Xdu; Xdu.resize(3, ip.npoints);
			Matrix3X Xdv; Xdv.resize(3, ip.npoints);

            eval->evalPointsPartial(localPts, surfacePoints,&X,&Xdu,&Xdv);
            eval->evalPointsPartial(localPts, oldPoints,&Y,NULL,NULL);
			/*
			std::ostringstream ostr;
            ostr << std::internal << std::setfill('0') << std::setw(5) << i;
			std::string folder = ip.out_path + std::string("/") + ostr.str();
            mkdir(folder.c_str(),0777);
            folder = folder + "/";
			std::string path = folder + std::string("cont_")+ std::to_string(iteration) + std::string("_") + 
				std::to_string(i)+std::string("fail") + std::to_string(failCount) + std::string(".obj");
			quiver(path,X,Y);

            for(int p = 0; p < ip.npoints; p++) {
                Scalar d = (Y.col(p) - X.col(p)).norm();
                maxMetric = fmax(d, maxMetric);
            }
			*/
            //std::cout << "Max Update Length: " << maxLocal << std::endl;
            //std::cout << "Max R^3 Length: " << maxMetric << std::endl;
            //std::cout << "Face Jumps: " << faceJump << std::endl;

            surfaceState.setCurrentModel(i,T, localPts, false);
			double average = 0.0;
            for(int p = 0; p < ip.npoints; p++) {
                face_idx[i*ip.npoints+p] = surfacePoints[p].face;
                patch_idx[i*ip.npoints+p] = surfacePoints[p].face;
                //tangBasis[i*ip.npoints+p] = basis;
                Scalar2 ur;
                ur.x = surfacePoints[p].u[0];
                ur.y = surfacePoints[p].u[1];
                inputData.setPoint(barriers.surface + i*ip.npoints+p, ur);

				Scalar2 us;
                us.x = 1.0; // / Xdu.col(p).stableNorm();
                us.y = 1.0; // / Xdv.col(p).stableNorm();
				average += Xdu.col(p).stableNorm()/ip.npoints;
                inputData.setPoint(barriers.surfaceScale + i*ip.npoints+p, us);
            }
			//std::cout << "Average Xdu length: " << average << std::endl;
            int dataId = 0;
            for(int p = 0; p < ip.npoints; p++) {

                int pidx = patch_idx[i*ip.npoints+p];
                for(int j = 0; j < 16; j++) {
                    int offset = bspline_control.at(16*pidx + j);
                    assert(offset < localPointCount);
                    BarrierIndex value = barriers.localRegistration + i*localPointCount + offset;
                    terms[0].set(6+j, i*ip.npoints+p, value); //3
                }

            }
            ///std::cout << i << "," << std::flush;
        }
        //std::cout << std::endl;
        copy();
    }



    //TODO Broken

    void Problem::implicit_update(int iteration, int failCount, MeshTopology* top, bool fake) {
        InputParams& ip = inputParams;
        inputData.retrieve();
        double error = 0;
        double modelError = 0;
        double oldError = 0;
        for(int i = 0; i < ip.nscans; i++) {
            //std::cout << "Start continuous update" << std::endl;
            std::vector<SurfacePoint> surfacePoints;
            std::vector<SurfacePoint> oldPoints;

            SubdivEvaluator::triplets_t dSdX, dSudX, dSvdX;
            BarrierIndex varOffset = barriers.registration + i*ip.nverts;
            Matrix3X T; T.resize(3,ip.nverts);
            //Matrix3X T_old; T_old.resize(3,ip.nverts);
            Matrix3X localPts(3,localPointCount);
            inputData.read(T,varOffset,ip.nverts, true);
            //inputData.read(T_old,varOffset,ip.nverts,false);
            inputData.read(localPts,barriers.localRegistration + i*localPointCount,localPointCount, true);
            //surfaceState.setCurrentModel(i,T_old, false);
            surfaceState.setCurrentModel(i,T,localPts, true); //false
            //std::cout << "current model set" << std::endl;

            //Matrix3X tangent; tangent.resize(3,ip.npoints);

            Matrix2X deltaU; deltaU.resize(2,ip.npoints);
            deltaU.setZero();

			//Matrix3X aS; aS.resize(3, ip.npoints);
			//Matrix3X aSdu; aSdu.resize(3, ip.npoints);
			//Matrix3X aSdv; aSdv.resize(3, ip.npoints);
			//eval->evalPointsPartial(localPts, surfacePoints,&aS,&aSdu,&aSdv);

            Scalar maxMetric = 0;
            Scalar respLocal = 0;
            Scalar maxLocal = 0;
            Scalar repMetric = 0;
            Scalar meanMetric = 0;
            Scalar meanLocal = 0;
            int faceJump = 0;

            double planarUpdate = 0;


            for(int p = 0; p < ip.npoints; p++) {
                inputData.resetUnknown2d(
                        barriers.surfaceUnknown + i*ip.npoints+p);
                Vector2 u;
                {
                    Scalar2 ur = inputData.readPoint2d(barriers.surface + i*ip.npoints+p);
                    u[0] = ur.x;
                    u[1] = ur.y;
                    SurfacePoint surfacePoint;
                    surfacePoint.u[0] = u[0];
                    surfacePoint.u[1] = u[1];
                    surfacePoint.face = face_idx[i*ip.npoints+p];
                    oldPoints.push_back(surfacePoint);
                }
                int f = face_idx[i*ip.npoints+p];
                int fnew;
                Vector2 unew;
                Scalar singleDistance;

                maxLocal = fmax(deltaU.col(p).norm(), maxLocal);

                int jumpCount = surfaceState.move(i, u, f, deltaU.col(p), unew, fnew, &singleDistance);
                if(f != fnew) faceJump++;
                f = fnew;
                u = unew;

                auto quad = top->quads.col(f);
                //face_idx[i*ip.npoints+p] = f;
                /*{
                    Scalar2 ur;
                    ur.x = u[0];
                    ur.y = u[1];
                    inputData.setPoint(i*ip.npoints+p, ur);
                }*/

                SurfacePoint surfacePoint;
                surfacePoint.u[0] = u[0];
                surfacePoint.u[1] = u[1];
                surfacePoint.face = f;//face_idx[i*ip.npoints+p];
                surfacePoints.push_back(surfacePoint);
            }

            Matrix3X X; X.resize(3, ip.npoints);
            Matrix3X Y; Y.resize(3, ip.npoints);

            eval->evalPointsPartial(localPts, surfacePoints,&X,NULL,NULL);
            eval->evalPointsPartial(localPts, oldPoints,&Y,NULL,NULL);


            for(int p = 0; p < ip.npoints; p++) {
                Scalar d = (Y.col(p) - X.col(p)).norm();
                maxMetric = fmax(d, maxMetric);
            }

            //std::cout << "Max Update Length: " << maxLocal << std::endl;
            //std::cout << "Max R^3 Length: " << maxMetric << std::endl;
            //std::cout << "Face Jumps: " << faceJump << std::endl;

            surfaceState.setCurrentModel(i,T, localPts, false);

            for(int p = 0; p < ip.npoints; p++) {
                face_idx[i*ip.npoints+p] = surfacePoints[p].face;
                patch_idx[i*ip.npoints+p] = surfacePoints[p].face;
                //tangBasis[i*ip.npoints+p] = basis;
                Scalar2 ur;
                ur.x = surfacePoints[p].u[0];
                ur.y = surfacePoints[p].u[1];
                inputData.setPoint(barriers.surface + i*ip.npoints+p, ur);
            }
            int dataId = 0;
            for(int p = 0; p < ip.npoints; p++) {

                int pidx = patch_idx[i*ip.npoints+p];
                for(int j = 0; j < 16; j++) {
                    int offset = bspline_control.at(16*pidx + j);
                    assert(offset < localPointCount);
                    BarrierIndex value = barriers.localRegistration + i*localPointCount + offset;
                    terms[0].set(6+j, i*ip.npoints+p, value); //3
                }

            }
            ///std::cout << i << "," << std::flush;
        }
        //std::cout << std::endl;
        copy();
    }

    void Problem::renormalize() {
        InputParams& ip = inputParams;
        inputData.retrieve();
        for(int v = 0; v < ip.nverts; v++) {
            SVector4 influence = closestJoints.col(v).cast<int>();
            /*
            float maximum = -1000000;
            for(int infu =   0; infu < 4; infu++) {
                int joint = influence[infu];
                float w = inputData.read(barriers.weights+inputParams.njoints*v + joint);
                maximum = fmax(maximum,w);
            }
            for(int infu =   0; infu < 4; infu++) {
                int joint = influence[infu];
                int offset = barriers.weights+inputParams.njoints*v + joint;
                float w = inputData.read(offset);
                w = w - maximum;
                inputData.write(w, offset);
            }
            */
            Scalar sum = 0;
            for(int infu =   0; infu < 4; infu++) {
                int joint = influence[infu];
                Scalar w = inputData.read(barriers.weights+inputParams.njoints*v + joint);
                sum += exp(w);
            }
            for(int infu =   0; infu < 4; infu++) {
                int joint = influence[infu];
                BarrierIndex offset = barriers.weights+inputParams.njoints*v + joint;
                Scalar w = inputData.read(offset);
                Scalar out = log(exp(w) / sum);
                inputData.write(out, offset);
            }
        }
        copy();

    }

    void Problem::getPersonStart(int query, int* start, int* end) {
        InputParams& ip = inputParams;
        int relativeScan = 0;
        int personId = 0;
        *start = -1;
        *end = -1;

        for(int i = 0; i < ip.nscans; i++) {
            if(personId == query && relativeScan == 0) {
                *start = i;
            }
            if(personId == (query+1) && relativeScan == 0) {
                *end = i;
                return;
            }
            relativeScan++;
            if(relativeScan == corpora.scanCount[personId]) {
                relativeScan = 0;
                personId++;
            }
        }
        *end = ip.nscans;
    }

    void Problem::propagateModel(MeshTopology* top) {
        InputParams& ip = inputParams;
        Model model;
        model.top = *top;
        retrieve(model);

        for(int i = 0; i < inputParams.nscans; i++) {
            Matrix3X T = model.scans[i].T_estimate;
            for(int j = 0; j < inputParams.nverts; j++) {
                Scalar3 vert;
                vert.x = T(0,j);
                vert.y = T(1,j);
                vert.z = T(2,j);
                //inputData.addUnknown(vert);
                BarrierIndex idx = barriers.registration + i*ip.nverts + j;
                inputData.setUnknown3d(idx, vert);
            }
        }

        for(int k = 0; k < inputParams.nscans; k++) {
            Matrix3X Tlocal_estimate = model.scans[k].Tlocal_estimate;
            //std::cout << "Compute local points for scan: " << k << std::endl;
            for(int i = 0; i < localPointCount; i++) {
                Scalar3 val;
                val.x = Tlocal_estimate(0,i);
                val.y = Tlocal_estimate(1,i);
                val.z = Tlocal_estimate(2,i);
                BarrierIndex index =
                        barriers.localRegistration + k*localPointCount + i;
                inputData.setUnknown3d(index, val);
            }

        }
        copy();
    }

    void Problem::dumpGraph() {

        {
            std::filebuf fb;
            fb.open("metadata.txt", std::ios::out);
            std::ostream outfile(&fb);
            for(int i = 0; i < terms.size(); i++) {
                terms[i].dumpMetatdata(outfile);
            }
        }
        {
            for(int i = 0; i < terms.size(); i++) {
                std::filebuf fb;
                std::string filename = "graphs/nodes";
                fb.open(filename+terms[i].getName()+".csv", std::ios::out);
                std::ostream outfile(&fb);
                outfile << "Id;Label\n";
                terms[i].dumpNodes(outfile);
            }

        }

        {

            for(int i = 0; i < terms.size(); i++) {
                std::filebuf fb;
                std::string filename = "graphs/edges";
                fb.open(filename+terms[i].getName()+".csv", std::ios::out);
                std::ostream outfile(&fb);
                outfile << "Source;Target;Label\n";
                terms[i].dumpEdges(outfile);
            }
        }

    }
    void Problem::dumpModel(int i, MeshTopology* top) {
        InputParams& ip = inputParams;
        Model model;
        model.top = *top;
        retrieve(model);

#ifdef FULL_DUMP
        std::string weightDir = ip.out_path + "/weights";
        mkdir(weightDir.c_str(),0777);
        model.write(weightDir + "/weights" + std::to_string(i) );
#endif
#ifdef VERBOSE_DUMP
        //evalMeanShapeTerm(model)
		//for(int q = 0; q <= ip.ndeformshapes; q+=9) {
		//	dump(model, *top, true, i, ip.ncomp,q);
		//}
		for(int c = 0; c <= ip.ncomp; c++) {
        	dump(model, *top, true, i, c,ip.ndeformshapes);
		}
#else
		dump(model, *top, true, i, ip.ncomp,ip.ndeformshapes);
#endif

#ifdef FULL_DUMP
        ////dump(model, top, false, 2*i+1);
        int mode[3] = {2,1,0};
        for(int p = 0; p < ip.npersons; p++) {
            int scanStart;
            int scanEnd;
            getPersonStart(p, &scanStart, &scanEnd);
            //std::cout << dumppath << " at " << scanStart << ":" << scanEnd << std::endl;

            std::ostringstream ostr;
            //ostr << std::internal << std::setfill('0') << std::setw(5) << p;
			ostr << corpora.folder[p];
            std::string folder = ip.out_path +std::string("/") + ostr.str();
            mkdir(folder.c_str(),0777);


            std::string dumppath =
                    std::string(folder + std::string("/p_")+std::to_string(i));

            std::string dumppath2 =
                    std::string(folder + std::string("/w_")+std::to_string(i));
            std::cout << "Writing to " << dumppath << std::endl;
            writePersonFbx(dumppath, mode,  &model, p, scanStart, scanEnd);
            //writePersonFbx(dumppath2, mode,  &model, p, 0, ip.nscans);
        }
#endif

    }

    void Problem::testUpdate(MeshTopology* top) {
        InputParams& ip = inputParams;
        Matrix3X localPts(3,localPointCount);
        inputData.read(localPts,barriers.localRegistration + 0*localPointCount,localPointCount);

        SurfacePoint p;
        p.face = rand() % top->num_faces();
        p.u.setConstant(0);

        Matrix3X target; target.resize(3,1);
        Matrix3X surf; surf.resize(3,1);
        target.col(0) = cloud[0].col(0);

        SurfacePoint q;
        Matrix3X pathCloud;
        eval->minimizeSurfacePoint(*top, p, localPts, target, q, &pathCloud);

        std::vector<SurfacePoint> pts;
        pts.push_back(q);
        eval->evalPointsPartial(localPts, pts,&surf);

        dumpCloudA(ip.out_path+std::string("/target"), 0, target);

        dumpCloudA(ip.out_path+std::string("/surf"), 0, surf);

        dumpCloudA(ip.out_path +std::string("/path"), 0, pathCloud);

    }

    void Problem::run() {
        InputParams& ip = inputParams;
        bool optimize = true;
        MeshTopology top;
        //loadObj("dummy.obj", &top, &control_vertices,NULL,SCALE);
        loadObj(ip.template_mesh_path, &top, &control_vertices,NULL,SCALE);
        //control_vertices = 0.1 * control_vertices;
        std::map<std::string, std::vector<size_t>> vertex_groups;
        loadVertexGroups(ip.template_vgroups_path, &vertex_groups, &border_groups, &weights); //vgroups.txt _fixmid.txt

        SparseMatrix sparseWeights = skeleton.extractWeightMatrix(weights, ip.nverts);

        headVertices.resize(ip.nverts);
        for(int i = 0; i < ip.nverts; i++) {
            headVertices[i] = false;
        }
        rigidVertices =  headVertices;


        //for(size_t i : vertex_groups["m_head"]) {
        //    headVertices[i] = true;
        //}
        skeleton.add(vertex_groups,&rigidVertices);



        //for (auto i : border_groups) {
        //    for (auto j : i.second) {
        //        headVertices[j] = false;
        //        rigidVertices[j] = false;
        //    }
        //}
        int headCount = 0;
        for(int i = 0; i < ip.nverts; i++) {
            if(headVertices[i]) headCount++;
        }
        std::cout << "Head Vertex Count: " << headCount << std::endl;
        //exit(0);


#ifdef USE_MARKERS
        loadVertexGroups(ip.template_marker_path, &marker_vertex_groups, NULL, NULL); //markers.txt.old
#endif

        std::cout << marker_vertex_groups.size() << std::endl;

        {
			std::cout << "Border_groups_size " << border_groups.size() << std::endl;
            std::vector<std::vector<int>> borderIndices(border_groups.size());
            for (auto i : border_groups) {
				std::cout << i.first.first << " : " << i.first.second << std::endl;
                int j = skeleton.getJointIndex(i.first.first, i.first.second);
				std::cout << j << std::endl;
                borderIndices[j].insert(borderIndices[j].end(), i.second.begin(), i.second.end());
            }

            Matrix4X jointDistance;
            closestJoints = top.findClosest(control_vertices, borderIndices, sparseWeights, &jointDistance);
        }
        {
            std::vector<std::vector<int>> borderIndices(border_groups.size());
            for (auto i : border_groups) {
                if(i.first.first == "m_root" || i.first.second == "m_root") continue;
                int j = skeleton.getJointIndex(i.first.first, i.first.second);
                borderIndices[j].insert(borderIndices[j].end(), i.second.begin(), i.second.end());
            }

            Matrix4X jointDistance;
            closestDeformJoints = top.findClosest(control_vertices, borderIndices, sparseWeights, &jointDistance);
        }

        Model model2, model3;



        if(!loadData(&top)) return;
        /* FIRST TWO TERMS HAVE TO STAY FIRST */
        terms.push_back(TermBuilder("A",6+16,&inputData)); //4+16
        terms.push_back(TermBuilder("B",4+16,&inputData));


        terms.push_back(TermBuilder("B2",7,&inputData));
        
		terms.push_back(TermBuilder("C2",REST_SIZE,&inputData));
        terms.push_back(TermBuilder("C",MODEL_SIZE,&inputData));

        terms.push_back(TermBuilder("PJ",PJOINT_SIZE,&inputData));


        if(!ip.freeze_weights) {
            terms.push_back(TermBuilder("D",2,&inputData));
            terms.push_back(TermBuilder("E",4,&inputData));
            terms.push_back(TermBuilder("F",16+1,&inputData)); //shape reg term
        }
        ////terms.push_back(TermBuilder("D",2,&inputData));//deform pseudo smooth term

        terms.push_back(TermBuilder("SR",2,&inputData)); //coeff reg term
		terms.push_back(TermBuilder("G2",16+3,&inputData)); //Repel term

        if(!ip.freeze_weights) {
#ifdef USE_RINGREG
            terms.push_back(TermBuilder("J",1+2*ip.max_vertex_ring,&inputData));
#else
            terms.push_back(TermBuilder("J",4,&inputData)); //JointRegTerm
#endif
            terms.push_back(TermBuilder("G",1,&inputData));
        }

        terms.push_back(TermBuilder("H",2+LOCAL_SIZE,&inputData)); //Local Term

		if(ip.useTemporal) {
        	terms.push_back(TermBuilder("T",4));
		}
		if(!ip.freeze_weights && ip.use_symmetry) {
        	terms.push_back(TermBuilder("S",2,&inputData));
        }
		 ////terms.push_back(TermBuilder("P",2,&inputData));
        terms.push_back(TermBuilder("PA",ip.freeze_weights ? 3 : 4,&inputData));


        //terms.push_back(TermBuilder(1+16*16));
        terms.push_back(TermBuilder("RA",12,&inputData));
        terms.push_back(TermBuilder("LV",2,&inputData));
        createUnknowns(&top);

        /* FIRST TWO TERMS HAVE TO STAY FIRST */
        createReducedSubdivTerm(&top);
        createLabelTerm(&top);

        createLabel2dTerm();
        createRestModelTerm();
        createModelTerm();

        createPJointTerm();
        if(!ip.freeze_weights) {
            createWeightPriorTerm();
            createWeightNormTerm();
            createMeanShapeTerm(&top);
        }
        createCoeffRegTerm();

		createGroundRepelTerm(&top);

        if(!ip.freeze_weights) {
#ifdef USE_RINGREG
            createJointRingReg();
#else
            createJointReg();
#endif
            centerReg();
        }

        createLocalTerm();
		
		if(ip.useTemporal) {
        	createTemporalPrior();
		}
		if(!ip.freeze_weights && ip.use_symmetry) {
        	createSymmetryTerm();
		}

        createAdvancedPoseTerm();
        relativeToAbsoluteTerm();
        linkVariableTerm();

        barriers.print();

        prepare();
        createOptParams();
#ifndef PRECOMPUTE_SIZE
        initOpt();
#else
        inputData.assertDims(precomputed_dims);
#endif

        w_fitSqrt = 0;
        w_surface = 1;
        const int THRESH = 30;
        const int THRESH_MUL = 6; //6

        if(optimize) {
            std::cout << "nparams: " << optinput.size() << std::endl;
            //Opt_ProblemInit(state, plan, &optinput[0]);
            inputData.checkpoint();
            manifold_update();
            propagateModel(&top);
            icp_update(&top, false);
            Opt_ProblemInit(state, plan, &optinput[0]);
            dumpModel(0,&top);

            int failCount = 0;
            for(int i = 0; i < inputParams.iteration_end; i++) { //150

                w_fitSqrt = i;
                w_surface = i > inputParams.jointly - 2
					? (i % inputParams.alternating == 0 ? 0 : 1)
					: 0;
                std::cout << "Iteration: " << i << std::endl;
				std::cout << "Surface Opt:  " << w_surface << std::endl;
                int result = 1;

                if(i == inputParams.raycasting) { //THRESH_MUL
                    for(int j = 0; j < ip.nscans; j++) {
                        surfaceState.disableRayTracing(j);
                    }
                }
				//if(i == 15) //60 15
				//	removeRigidJoints();
                //if(w_surface > 0 && i % inputParams.alternating == 0) { //%5
                //    icp_update(&top, false);
                //}


                auto oldFaceIdx = face_idx;
                auto oldPatchIdx = patch_idx;
                auto oldTangBasis = tangBasis;
                auto oldModel2Data = model2Data;
                terms[0].checkpoint();
                if(ip.useModel2Data) {
                    terms[1].checkpoint();
                }
                inputData.checkpoint();

				Cost cpu_cost(&inputData, &terms);
                while(true) {
                    inputData.dim();
                    double cost = Opt_ProblemCurrentCost(state, plan, &optinput[0]);
					double dataBefore = cpu_cost.compute(w_fitSqrt+1, w_surface);
					Opt_ProblemHalfstep(state, plan, &optinput[0]);
					double dataAfter = cpu_cost.compute(w_fitSqrt+1, w_surface);
                    double euclidCost = Opt_ProblemCurrentCost(state, plan, &optinput[0]);
                    
					//if(w_surface == 0)
					manifold_update();
					double dataManifold = cpu_cost.compute(w_fitSqrt+1, w_surface);

					double afterManifold = Opt_ProblemCurrentCost(state, plan, &optinput[0]);

                    double afterContinuous = afterManifold;
                    if(w_surface) { //8
						if(/*result == 0 ||*/ (i % inputParams.dump_progress == 0 && i > 0))
							dumpModel(1000*i+failCount,&top);
                        continuous_update(i, failCount, &top, false);
						if(/*result == 0 ||*/ (i % inputParams.dump_progress == 0 && i > 0))
							dumpModel(100000*i+failCount,&top);
                        afterContinuous = Opt_ProblemCurrentCost(state, plan, &optinput[0]);
                    } else {
                        afterContinuous = Opt_ProblemCurrentCost(state, plan, &optinput[0]);
                    }
					double dataCont = cpu_cost.compute(w_fitSqrt+1, w_surface);

                    double newcost = Opt_ProblemCurrentCost(state, plan, &optinput[0]);
                    std::cerr << i << " : " << cost << " -> " << euclidCost << " -> " << afterManifold << " -> " <<
                              afterContinuous << " -> " << newcost << std::endl;
					std::cerr << "Data " << i << " : " << dataBefore << " -> " 
							  << dataAfter << " -> " <<
                              dataManifold << " -> " << dataCont << std::endl;

                    int accepted = Opt_ProblemDecreaseAccepted(state,plan,&optinput[0]);
                    if(accepted == 1 || result == 0) {
                        failCount = 0;
                        std::cout << "TR Step going to be  accepted" << std::endl;
                        result = Opt_ProblemUpdate(state, plan, &optinput[0]);
                        //        	if(i == inputParams.raycasting) { //THRESH_MUL
                        //	for(int j = 0; j < ip.nscans; j++) {
                        //	surfaceState.disableRayTracing(j);
                        //	}
                        //}
                        if(w_surface == 0 && i >= ip.iteration_data - 1) {
                            icp_update(&top, i > 4*THRESH);
                        }
                        //if(w_surface > 0 && i % inputParams.alternating == 0) { //%5
                        //	icp_update(&top, false);
                        //}
                        break;
                    }
                    std::cout << "TR Step going to be  rejected" << std::endl;
                    result = Opt_ProblemUpdate(state, plan, &optinput[0]);

                    terms[0].restore();
                    if(ip.useModel2Data) {
                        terms[1].restore();
                    }
                    inputData.restore();
                    face_idx = oldFaceIdx;
                    patch_idx = oldPatchIdx;
                    tangBasis = oldTangBasis;
                    model2Data = oldModel2Data;
                    failCount++;

                    if(result == 0) {
                        if( i < 2*THRESH-1) {
                            i = 2*THRESH-1;
                        }
                        break;
                    }

                }
                if(/*result == 0 ||*/ (i % inputParams.dump_progress == 0 && i > 0)) {
                    dumpModel(i+1,&top);
                }
                if(result == 0) {
                    break;
                }
            }

            std::cout << "Try Free" << std::endl;
            Opt_PlanFree(state, plan);
            std::cout << "Try Delete" << std::endl;
            Opt_ProblemDelete(state, problem);
        }
        std::cout << "Try Release" << std::endl;
        release();
        std::cout << "End Run" << std::endl;

    }
    void Problem::dump(const Model& model, const MeshTopology& top, bool real, int iter, int pc_limit, int qc_limit) {
        InputParams& ip = inputParams;


        int personId = 0;
        int relativeScan = 0;


        Matrix3X mymean = model.mean;
        Matrix3X mymeanjoint = model.joints;
        saveObj(ip.out_path +std::string("/c")+ std::to_string(iter) + std::string("mean.obj"), &top, &mymean);
        dumpCloud(ip.out_path +std::string("/cJ_")+ std::string("mean"), iter, mymeanjoint);
        for(int j = 0; j < pc_limit; j++) {
            Matrix3X pc = model.pcs[j];
            Matrix3X temp = mymean + 3*pc;
            saveObj(ip.out_path +std::string("/c")+ std::to_string(iter) + std::string("pc") + std::to_string(j) + std::string(".obj"), &top, &temp);
            const Matrix3X& jpc = model.jointpcs[j];
            Matrix3X jtemp = mymeanjoint + 3*jpc;
            dumpCloud(ip.out_path +std::string("/cJ") + std::string("pc") + std::to_string(j) + "_", iter, jtemp);
        }

        for(int i = 0; i < ip.nscans; i++) {
            std::ostringstream ostr;
            //ostr << std::internal << std::setfill('0') << std::setw(5) << i;
            ostr << corpora.folder[corpora.scan_to_person_id[i]] << "."
                                   << corpora.scan_names[i];

			std::string folder = ip.out_path + std::string("/") + ostr.str();
            mkdir(folder.c_str(),0777);
            folder = folder + "/";

            Matrix3X C1 = model.mean.replicate(1,1); //C1.setConstant(0);
            Matrix3X skinned(3,ip.nverts);
            Matrix3X unlinked_mesh(3, ip.nverts);

            const Person& person = model.persons[personId];

            if(relativeScan == 0) {
                std::cout << personId << " : " << std::endl;
                for(int j = 0; j < pc_limit; j++) {
                    std::cerr << person.coeff[j] << ", ";
                }
                std::cerr << std::endl;
            }

            Matrix3X J1;
            if(real) {
                J1 = model.joints;
            } else {
                J1 = person.joints;
            }


            for(int j = 0; j < pc_limit; j++) {
                Scalar coeff = person.coeff[j];
                const Matrix3X& pc = model.pcs[j];
                C1 = C1 + coeff*pc;
                const Matrix3X& jpc = model.jointpcs[j];
                if(real) {
                    J1 = J1 + coeff*jpc;
                }
            }
#ifdef FULL_DUMP
			dumpCloudA(folder + std::string("jrest")+ std::to_string(iter) + std::string("_") + std::to_string(i), i, J1);


            saveObj(folder + std::string("Prest")+ std::to_string(iter) + std::string("_") + std::to_string(i)+std::string(".obj"), &top, &C1);
#endif
            std::vector<Vector3> skinjoints;
            std::vector<Matrix3> localrot;
            rod2Mat(model.scans[i].posebase,localrot);

            for(int j = 0; j < ip.njoints; j++) {
                Matrix3 so3; so3.setIdentity();
                Vector3 r = model.scans[i].pose.col(j);
                so3(0,1) = -r[2];
                so3(1,0) = r[2];
                so3(0,2) = r[1];
                so3(2,0) = -r[1];
                so3(1,2) = -r[0];
                so3(2,1) = r[0];
                localrot[j] = localrot[j]*so3;
            }

            std::vector<Matrix3> absrot;
            rod2Mat(model.scans[i].absolutePosebase,absrot);
            for(int j = 0; j < ip.njoints; j++) {
                Matrix3 so3; so3.setIdentity();
                Vector3 r = model.scans[i].absolutePose.col(j);
                so3(0,1) = -r[2];
                so3(1,0) = r[2];
                so3(0,2) = r[1];
                so3(2,0) = -r[1];
                so3(1,2) = -r[0];
                so3(2,1) = r[0];
                absrot[j] = absrot[j]*so3;
            }


            unravel(J1,skinjoints);
            Matrix3X shapeForm(C1);

            for(int j = 0; j < ip.ndeformjoints; j++) {
                for(int k = 0; k < qc_limit; k++) {
                    for(int v = 0; v < ip.nverts; v++) {
                        int joint = closestDeformJoints(j,v);
                        Scalar melem = localrot[joint](k/3,k%3);
                        if(k % 4 == 0) melem -= 1;
                        C1.col(v) = C1.col(v) + melem*model.deform[9*j+k].col(v);
                    }
                }
            }

            std::vector<Matrix4> absolutes(ip.njoints);

            skeleton.calculate_kinematic_chain(skinjoints,
                                               localrot, absolutes);

            std::vector<Matrix4> unlinked(ip.njoints);
            skeleton.calculate_kinematic_chain(skinjoints,
                                               absrot, unlinked, model.scans[i].joints);

            for(int v = 0; v < ip.nverts; v++) {
                VectorX wv = model.weights.col(v); wv.setZero();
                SVector4 influence = closestJoints.col(v).cast<int>();
                for(int infu =   0; infu < 4; infu++) {
                    int joint = influence[infu];
                    wv[joint] = model.weights(joint,v);
                }
                Vector3 vert = skeleton.skin_vertex(absolutes,wv,C1.col(v));
                skinned.col(v) = vert + model.scans[i].trans;

                Vector3 unlinked_vert = skeleton.skin_vertex(unlinked,wv,C1.col(v));
                unlinked_mesh.col(v) = unlinked_vert + model.scans[i].absoluteTrans;


            }


            Matrix3X jointsOut; jointsOut.resize(3,absolutes.size());
            for(int j = 0; j < absolutes.size(); j++) {
                Vector4 joint;
                joint[0] = skinjoints[j][0];
                joint[1] = skinjoints[j][1];
                joint[2] = skinjoints[j][2];
                joint[3] = 1;

                Vector4 jointRes = absolutes[j]*joint;
                jointsOut(0,j) = jointRes[0] + model.scans[i].trans(0,0);
                jointsOut(1,j) = jointRes[1] + model.scans[i].trans(1,0);
                jointsOut(2,j) = jointRes[2] + model.scans[i].trans(2,0);
            }
#ifdef FULL_DUMP
            dumpCloudA(folder + std::string("pj")+ std::to_string(iter) + std::string("_") + std::to_string(i), i, person.joints);
            dumpCloudA(folder + std::string("j")+ std::to_string(iter) + std::string("_") + std::to_string(i), i, jointsOut);
			//dumpCloudA(folder + std::string("sj")+ std::to_string(iter) + std::string("_") + std::to_string(i), i, jointsOut);
			if(corpora.person_has_labels_2d[personId]) {
            for(int v = 0; v < ip.nviews; v++){
                Matrix3X sj = model.scans[i].joints;
                sj.colwise() += model.scans[i].trans.col(0);

                Matrix3X footpoints; footpoints.resize(3,ip.njointlabels);
                for(int k = 0; k < ip.njointlabels; k++) {
	                int camId = corpora.person_fist[personId] ? 1 : 0; 
	                std::string camPath = ip.camPaths[camId];

                    Camera cam(camPath, v);
                    footpoints.col(k) = cam.projectJointOnRay(model.scans[i].label2dRays[v].col(k), sj.col(k));
                }

                Matrix3X sjn = -footpoints + sj;

                dumpCloudA(folder + std::string("latent joint")+
                           std::to_string(iter) + std::string("_") +
                           std::to_string(v)+"_", i, sj, &sjn);
                //dumpCloudA(folder + std::string("sj2d")+
                //           std::to_string(iter) + std::string("_") +
                //           std::to_string(v) + "_",
                //           i, footpoints);
				dumpCloudA(folder + std::string("jointProjection")+
                           std::to_string(iter) + std::string("_") +
                           std::to_string(v)+"_", i, footpoints, &sjn);
            }
			}
#endif
            if(ip.useModel2Data) {
                Matrix3X m2d; m2d.resize(3,ip.nverts);
                for(int k = 0; k < ip.nverts; k++) {
                    m2d.col(k) = cloud[i].col(model2Data[i][k]);
                }
                Matrix3X m2dn = model.scans[i].T - m2d;
                dumpCloudA(folder + std::string("m2d")+
                           std::to_string(iter) + std::string("_") + std::to_string(i),
                           i, m2d,&m2dn);
            }
            saveObj(folder + std::string("c")+ std::to_string(iter) + std::string("_pc") + 
				std::to_string(pc_limit) + std::string("_") + std::string("_qc") + 
                std::to_string(qc_limit) + std::string("_") +
				std::to_string(i)+std::string(".obj"), &top, &skinned);
			//dumpSubivA(folder + std::string("csub")+ std::to_string(iter) + std::string("_") + std::to_string(i)+std::string(".obj"), 
            //    iter, eval, skinned);
            saveObj(folder + std::string("u")+ std::to_string(iter) + std::string("_") + std::to_string(i)+std::string(".obj"), &top, &unlinked_mesh);

            saveObj(folder + std::string("T")+ std::to_string(iter) + std::string("_") + std::to_string(i)+std::string(".obj"), &top, &model.scans[i].T);
			//dumpSubivA(folder + std::string("Tsub")+ std::to_string(iter) + std::string("_") + std::to_string(i)+std::string(".obj"), 
			//	iter, eval, model.scans[i].T);
#ifdef FULL_DUMP
            saveObj(folder + std::string("Trest")+ std::to_string(iter) + std::string("_") + std::to_string(i)+std::string(".obj"), &top, &model.scans[i].Trest);
#endif
            Matrix3X normals = cloud[i] - model.scans[i].patchPoints;
            dumpCloudA(folder + std::string("S")+ std::to_string(iter) + std::string("_") + std::to_string(i), i, model.scans[i].patchPoints,&normals);
            Matrix3X normLabel = vmarkers[i] - model.scans[i].labels;


			MeshTopology patchMesh;
			Matrix3X patchVertices;
			eval->createPatches(model.scans[i].Tlocal,
    			&patchMesh, &patchVertices);
			saveObj(folder + std::string("Patch")+ std::to_string(iter) + std::string("_") + std::to_string(i)+std::string(".obj"), &patchMesh, &patchVertices);


			Matrix3X tangentNormals = cloud[i] - model.scans[i].tangentStep;
			dumpCloudA(folder + std::string("TanS")+ std::to_string(iter) + std::string("_") + std::to_string(i),i,
				model.scans[i].tangentStep,&tangentNormals);
#ifdef FULL_DUMP

            dumpCloudA(folder + std::string("Tlabels")+ std::to_string(iter) + std::string("_") + std::to_string(i), i, model.scans[i].labels, &normLabel);
            dumpCloudA(folder + std::string("Tlabels_init")+ std::to_string(iter) + std::string("_") + std::to_string(i), i, model.scans[i].labels_initial);
            dumpCloudA(folder +std::string("Tlabels_partial")+ std::to_string(iter) + std::string("_") + std::to_string(i), i, model.scans[i].labels_partial);
			if(corpora.person_has_labels_2d[personId]) {
            for(int v = 0; v < ip.nviews; v++){
                Matrix3X ray = model.scans[i].label2dRays[v];
                ray.colwise() -= model.scans[i].label2dRays[v].col(ip.njointlabels+ip.nverts);
                ray *= 10;

                dumpCloud(folder + std::string("Rays")+ std::to_string(iter) + std::string("_") + std::to_string(v) + "_",
                          i,  model.scans[i].label2dRays[v], &ray);
            }
			}
            Matrix3X customEval = eval->evalPointsCustom(model.scans[i].T);
            Matrix3X partialEval = eval->evalPointsPartial(model.scans[i].Tlocal);
#endif
            relativeScan++;
            if(relativeScan == corpora.scanCount[personId]) {
                relativeScan = 0;
                personId++;
            }
        }
    }

    void Problem::evalJointReg(const Model& model) {
        InputParams& ip = inputParams;

        double jointReg = 0;
        for (size_t i = 0; i < ip.njoints-1; i++) {
            auto key = skeleton.getJoint(i + 1);
            auto secondKey = key;
            if (border_groups.find(key) == border_groups.end()) {
                secondKey.first = key.second;
                secondKey.second = key.first;
            }
            auto vertex_ring = border_groups.at(secondKey);

            {
                Vector3 jp = s2v(inputData.readPoint3d(barriers.jointPrior + i));
                Matrix3X jv; jv.resize(3,1);
                inputData.read(jv, barriers.joints + i + 1, 1);
                for (size_t j = 0; j < vertex_ring.size(); j++) {
                    size_t vertex_index = vertex_ring[j];
                    Vector3 vp = s2v(inputData.readPoint3d(barriers.meanShapePrior + vertex_index));
                    Matrix3X vv; vv.resize(3,1);
                    inputData.read(vv, barriers.shape + vertex_index, 1);

                    for(int ii = 0; ii < 3; ii++) {
                        Scalar delta = (jp(ii,0)-jv(ii,0)) - (vp(ii,0) - vv(ii,0));
                        jointReg += delta*delta;
                    }
                }
            }
            /*~/projects/Optlang/Opt/tests/shapefit/graph.csv
            for(int c = 0; c < ip.ncomp; c++) {
                edge[0] = barriers.zero;
		        //auto lv = global.jointRegressor.block<3, 1>(3 * i, 0);
                edge[1] = barriers.joints + ip.njoints*(c+1) + i + 1;
		        //auto& R = global.RSkel[i];
		        for (size_t j = 0; j < vertex_ring.size(); j++) {
			        size_t vertex_index = vertex_ring[j];
                    edge[2] = barriers.zero;
			        //auto& u = control_vertices.col(vertex_index); //.cast<Diffable>()
                    edge[3] = barriers.shape + (c+1)*ip.nverts + vertex_index;
                    //auto& v = global.meanShape.col(vertex_index);
			        //residual.col(r_idx) = factor / sqrt(vertex_ring.size() + global.joint_count)*((v - lv) - R*(u - lu)); //1*hyperparam
			        //r_idx++;
                    terms[tid].add(edge);
		        }
            }
            */
        }

        jointReg *= 1;
        std::cout << "Joint Reg Error: " << jointReg << std::endl;
    }
