#pragma once

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


class Problem {
private:
    const size_t REST_SIZE;
    const size_t MODEL_SIZE;
    const size_t PJOINT_SIZE;
    const float SCALE;
    const int LAMBDA;

    InputParams inputParams;
    std::vector<TermBuilder> terms;
    Corpora corpora;

    //TermBuilder dataTerm;
    //TermBuilder weightNormTerm;
    //TermBuilder coeffReg;
    //TermBuilder pcReg;
    Barriers barriers;

    std::vector<void*> optinput;

    InputData inputData;
    Opt_State* state;
    Opt_Problem* problem;
    Opt_Plan* plan;

    
    std::map<std::pair<std::string, std::string>, std::vector<size_t>> border_groups;
    hyp::Skeleton skeleton;
    Weights weights;
    Matrix3X control_vertices;
    IMatrix4X closestJoints;
    IMatrix4X closestDeformJoints;
    std::vector<bool> headVertices;
    std::vector<bool> rigidVertices;
    
    std::unique_ptr<SubdivEvaluator> eval;

    int localPointCount;
    std::vector<int> face_idx;
    std::vector<int> patch_idx;
    std::vector<int> bspline_control;

    std::vector<Matrix22> tangBasis;
    std::vector<int> bsplines_idx;
    std::vector<int> localIdx;
    std::vector<int> labelSubdiv_idx;
	std::vector<std::vector<int>> model2Data;
    std::map<std::string, std::vector<size_t>> marker_vertex_groups;

    
    std::vector<MeshTopology> topologies;
    std::vector<Matrix3X> vmarkers;
    std::vector<bool> bmarkers;
    std::vector<Matrix3X> vreal;
    std::vector<Matrix3X> cloud;
    std::vector<Matrix3X> cloudNorm;
        
    //std::vector<std::vector<int>> markerIdxs;

    SurfaceState surfaceState;

    Scalar w_fitSqrt;
    Scalar w_surface;
    int termIdx;
	int reverse_offset;
	int jointSPTermIdx;
	int modelTermIndex;
    //Scalar trustRadius;

    Matrix3X calcJointInit();

	MatrixX calcSkinningWeights(int njoints, int nverts);

    void createModelWeights(const Model& initModel, int isFree);
	bool loadData(MeshTopology* top);
        
    void createUnknowns(MeshTopology* top);


    void createReducedSubdivTerm(MeshTopology* top);

	void createReverseDataEdges(MeshTopology* top, int tid);
    
	void createLabelTerm(MeshTopology* top);

	void createLabel2dTerm();
    void createRestModelTerm();

    bool jointInInfluence(int vid, int joint);
    void createModelTerm();
    
    void createAdvancedPoseTerm();
	void createTemporalPrior();    
    void createPJointTerm();
    void createWeightPriorTerm();
    
    void createWeightNormTerm();
    void createMeanShapeTerm(MeshTopology* top);
    
    void evalMeanShapeTerm(const Model& model);
    
    void createCoeffRegTerm();

	void createGroundRepelTerm(MeshTopology* top);

	void createJointRingReg();

    void createJointReg();

    void createPoseReg();

    void centerReg();

    void createShapeReg();
    void createLocalTerm();
    
	
	void createSymmetryTerm();

    void relativeToAbsoluteTerm();
    void linkVariableTerm();


    void initOpt();

    void prepare();

    void copy();

    void createOptParams();

    void retrieve(Model& model);
    
    void release();

    void manifold_update() ;

	void removeRigidJoints() ;
	void reverse_icp_update(MeshTopology* top);


    void icp_update(MeshTopology* top, bool alternating);

	void implicit_update(int iteration, int failCount, MeshTopology* top, bool fake);

    void continuous_update(int iteration, int failCount, MeshTopology* top, bool fake);
    void renormalize() ;

    void getPersonStart(int query, int* start, int* end);

    void propagateModel(MeshTopology* top);

    void dumpGraph();
    void dumpModel(int i, MeshTopology* top);
    void testUpdate(MeshTopology* top);

public:
    Problem(const InputParams& ip, const Corpora& _corpora);
    void run();
    void dump(const Model& model, const MeshTopology& top, bool real, int iter, int pcs, int qcs);
    
    void evalJointReg(const Model& model);
};
