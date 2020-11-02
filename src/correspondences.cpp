#include "correspondences.h"
#include "MeshTopology.h"

SurfaceState::SingleState::SingleState(int _index, SubdivEvaluator* eval, std::shared_ptr<RayTracer> _rayTracer, MeshTopology* mesh, Matrix3X* scan, Matrix3X* scan_normals) {
	this->eval = eval;
	this->mesh = mesh;
	this->scan = scan;
	this->scan_normals = scan_normals;
	this->npts = eval->controlVertexPtex.size();
    this->index = _index;
	this->modelFaceMask.resize(this->mesh->num_faces());
#ifndef USE_NORMAL
	data.resize(npts, 3);
#else
    data.resize(npts, 6);
#endif
	rayTracer = _rayTracer;
	shooting = true;
}

int SurfaceState::SingleState::move(Vector2 u, int f, Vector2 d, Vector2& out, int& fout, 
        Scalar* distance,std::vector<Vector3>* path) {

	Vector2 dir;
	//std::cout << u.transpose() << "->" << d.transpose() << std::endl;
	int jumps = eval->increment_u_crossing_edges(*mesh, localPts, f, u, d, &fout, &out, &dir, distance, path);
	//int jumps = eval->increment_u_crossing_edges(*mesh, pts, f, u, d, &fout, &out, &dir, distance, path);	

	//int jumps = eval->increment_accurate(*mesh, pts, f, u, d, &fout, &out, &dir, distance, path);
	//int jumps = eval->increment_u_crossing_edges(*mesh, pts, f, u, d, &fout, &out, &dir, distance, path);
    return jumps;
}

Vector3 SurfaceState::SingleState::moveTangent(Vector2 u, int f, Vector2 d) {
    SurfacePoint sp;
	sp.face = f;
	sp.u = u;
	std::vector<SurfacePoint> sps;
	sps.push_back(sp);

    Matrix3X outS(3, 1);
	Matrix3X out_Su(3, 1);
	Matrix3X out_Sv(3, 1);

    eval->evaluateSubdivSurface(pts,
                sps, &outS, NULL, NULL, NULL,
                &out_Su, &out_Sv);

    Matrix3X result = outS + d[0]*out_Su + d[1]*out_Sv;
    return result.col(0);
}

void SurfaceState::SingleState::readModelPoints(Matrix3X& outS, Matrix3X& n) {
	outS.resize(3, npts);
#ifndef USE_NORMAL
    outS = eval->evalPointsPartial(localPts);
    /*eval->evaluatePtex(pts,
		&outS,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
        true);*/
#else
	Matrix3X out_Su(3, npts);
	Matrix3X out_Sv(3, npts);
	/*
	eval->evaluatePtex(pts,
		&outS,
		NULL,
		NULL,
		NULL,
		&out_Su,
		&out_Sv,
        true);
    */
    eval->evalPointsPartial(localPts,eval->controlVertexPtex,&outS,&out_Su,&out_Sv);
    
	n.resize(3, npts);
	n.setZero();
	
	if(scan_normals != NULL) {
		for (int i = 0; i < npts; i++) {
			n.col(i) = out_Su.col(i).cross(out_Sv.col(i));
			n.col(i).normalize();
		}
	}
#endif
}

Scalar SurfaceState::SingleState::evaluateCurrent(int dataIndex, Vector2 u, int face) {
#ifndef USE_NORMAL
    Vector3 p;
	Matrix3X outS(3, 1);
    SurfacePoint sp;
	sp.face = face;
	sp.u = u;
	std::vector<SurfacePoint> sps;
	sps.push_back(sp);
	eval->evaluateSubdivSurface(pts,
		sps,
		&outS);
    return (w*(p - scan->col(dataIndex))).squaredNorm();
#else
	Vector3 p, n;
	Matrix3X outS(3, 1);
	Matrix3X out_Su(3, 1);
	Matrix3X out_Sv(3, 1);
	SurfacePoint sp;
	sp.face = face;
	sp.u = u;
	std::vector<SurfacePoint> sps;
	sps.push_back(sp);
	/*
	eval->evaluateSubdivSurface(pts,
		sps,
		&outS,
		NULL,
		NULL,
		NULL,
		&out_Su,
		&out_Sv);
	*/
	
	eval->evalPointsPartial(localPts,sps,&outS,&out_Su,&out_Sv);
	Vector3 scanNormal; scanNormal.setZero();
	n.setZero();
	if(scan_normals != NULL) {
		n = out_Su.col(0).cross(out_Sv.col(0));
		n.normalize();
		scanNormal = scan_normals->col(dataIndex);
		out_Su.normalize();
		out_Sv.normalize();
	}
	p = outS.col(0);
	//Scalar o1 = dn(0)*out_Su(0,0)+dn(1)*out_Su(1,0)+dn(2)*out_Su(2,0);
	//Scalar o2 = dn(0)*out_Sv(0,0)+dn(1)*out_Sv(1,0)+dn(2)*out_Sv(2,0);
	//Scalar resN = nw*o1*nw*o1 + nw*o2+nw*o2;
	//return resN + w*(p - scan->col(dataIndex)).squaredNorm();
	if(shooting) {
		return w*(p - scan->col(dataIndex)).squaredNorm();
	} else {
		return (nw*(n - scanNormal)).squaredNorm() + (w*(p - scan->col(dataIndex))).squaredNorm();
	}
	
	//Scalar en = nw*(n(0)*dn(0)+n(1)*dn(1)+n(2)*dn(2)-1);
	//en = en*en;
	//return en + w*(p - scan->col(dataIndex)).squaredNorm();
        
#endif
}

void SurfaceState::SingleState::setCurrentModel(Matrix3X& pts, Matrix3X& localPts, const std::vector<int>& excludedVertices, bool rebuild_tree) {

 	this->pts = pts;
    this->localPts = localPts;
	for(int i = 0; i < this->mesh->num_faces();i++) {
        this->modelFaceMask[i] = true;
    }
	mesh->computeMaskFaces(excludedVertices, &modelFaceMask);


	Matrix3X outS, n;
    //std::cout << "Start Point Sampling" << std::endl;
	readModelPoints(outS, n);
    //std::cout << "End Point Sampling" << std::endl;
	data.block(0, 0, npts, 3) = w*outS.transpose();

	//HACK: move masked sample points far away. This way, they will never be nearest neighbors
	for(int i = 0; i < npts; i++) {
		int face = eval->controlVertexPtex[i].face;
		if(!this->modelFaceMask[face]) {
			data(i,0) = 1000000;
			data(i,1) = 1000000;
			data(i,2) = 1000000;
		}
	}
#ifdef USE_NORMAL
	data.block(0, 3, npts, 3) = nw*n.transpose();
#endif

	/*for (int i = 0; i < 3; i++) {
	for (int j = 0; j < npts; j++) {
	data(j, i) = w*outS(i, j);
	data(j, i+3) = nw*n(i, j);
	}
	}*/

	//mat_index.reset(new my_kd_tree_t(6, data, 10));

    //std::cout << "Start Tree Building" << std::endl;
    if(rebuild_tree) {
        mat_index.reset(new my_kd_tree_t(6, data, 10));
	    mat_index->index->buildIndex();
    }
    //std::cout << "End Tree Building" << std::endl;
	/*Scalar out_dists_sqr;
	size_t   ret_index;
	nanoflann::KNNResultSet<Scalar> resultSet(1);
	resultSet.init(&ret_index, &out_dists_sqr);
	double query_fuse[6];
	for (int i = 0; i < 3; i++) {
	query_fuse[i] = w*outS(i, 100);
	query_fuse[3 + i] = nw*n(i, 100);
	}
	for (int i = 0; i < 6; i++) {
	std::cout << query_fuse[i] << std::endl;
	}
	std::cout << data.row(100) << std::endl;

	mat_index->index->findNeighbors(resultSet, &query_fuse[0], nanoflann::SearchParams(10));
	std::cout << ret_index << 
	<< out_dists_sqr << std::endl; */
	//dumpCloud("out/test", this->index, outS, &n); //&n
	//exit(0);
	
	
	//eval->generate_refined_mesh(pts, 3, &resultTopology, &subdVertices);
	char cfg[50];
	strcpy(cfg,"");
	rayTracer->device_cleanup();
	//rayTracer.device_init(cfg, subdVertices, resultTopology);
	rayTracer->device_init(cfg, pts, *mesh);
}

void SurfaceState::SingleState::alternatingUpdate(int dataIndex, Vector2 u, int f, Vector2 d, Vector2& out, int& fout) {
	Vector3 query_pt = scan->col(dataIndex);
	SurfacePoint p, q;
	p.u = d;
	p.face = f;
	
	Matrix3X pt; pt.resize(3,1);
	pt.col(0) = query_pt;
	
	eval->minimizeSurfacePoint(*this->mesh, p, localPts, pt, q, NULL);
	out = q.u;
	fout = q.face;
}

bool SurfaceState::SingleState::discreteUpdate(int dataIndex, Vector2& out, int& fout) {
	//shooting = false;
	
	Vector3 query_pt = scan->col(dataIndex);
    #ifdef USE_NORMAL
	Vector3 query_normal; query_normal.setZero();
	
	if(scan_normals != NULL)
		query_normal = scan_normals->col(dataIndex);
	
	//if(index == 0) {
	//std::cout << "UPDATE " << dataIndex << std::endl;
	
	Vector3 rayHitPos;
	if(shooting) {
		Vector2 utemp; utemp.setConstant(-1);
		int ftemp = -1;
		rayHitPos = rayTracer->shootNormals(query_pt, query_normal, &utemp, &ftemp);
		if(ftemp >= 0 && this->modelFaceMask[ftemp]) {
			assert(out(0,0) >= 0.0 && out(0,0) <= 1.0 && out(1,0) >= 0.0 && out(1,0) <= 1.0);
			assert(utemp(0,0) >= 0.0 && utemp(0,0) <= 1.0 && utemp(1,0) >= 0.0 && utemp(1,0) <= 1.0);
			Scalar oldDist = evaluateCurrent(dataIndex, out, fout);
			Scalar tempDist = evaluateCurrent(dataIndex, utemp, ftemp);
			if(tempDist < oldDist) {
				out = utemp;
				fout = ftemp;
				return true;
			}
		}
	} else {
		rayHitPos = query_pt;
	}
	//} 
	
	Scalar query_fuse[6];
    
	for (int i = 0; i < 3; i++) {
		query_fuse[i] = w*rayHitPos[i];
		query_fuse[3 + i] = nw * query_normal[i];
	}

    #else
	Scalar query_fuse[3];
	for (int i = 0; i < 3; i++) {
		query_fuse[i] = w*rayHitPos[i];
	}
    #endif

	Scalar out_dists_sqr = -1.0f;
	size_t   ret_index = 0;
	nanoflann::KNNResultSet<Scalar> resultSet(1);
	resultSet.init(&ret_index, &out_dists_sqr);
	mat_index->index->findNeighbors(resultSet, &query_fuse[0], 
		nanoflann::SearchParams(10)); 
	
	Scalar oldDist = evaluateCurrent(dataIndex, out, fout);
	
	Vector2 outnew = eval->controlVertexPtex[ret_index].u;
	int fnew = eval->controlVertexPtex[ret_index].face;
	Scalar realDist = evaluateCurrent(dataIndex, outnew, fnew);
	
	assert(out_dists_sqr != -1.0f);
		
	//if(fabs(out_dists_sqr - realDist) > 0.001) {
	//	std::cout << "knn distance deviates from computed distance!" << std::endl;
	//}

	if (oldDist > out_dists_sqr) {
		out = outnew;
		fout = fnew;
		//std::cout << fout << " : " << out.transpose() << std::endl;
        return true;
	}
	else {
        return false;
	}
	assert(0.0f <= out[0] && out[0] <= 1.0f);
	assert(0.0f <= out[1] && out[1] <= 1.0f);
	assert(0 <= fout && fout <= eval->ptexnum);
}

	void SurfaceState::SingleState::disableRayTracing() {
		shooting = false;
	}

	SurfaceState::SurfaceState() {
			
	}
	void SurfaceState::addData(SubdivEvaluator* eval, std::shared_ptr<RayTracer> _rayTracer, MeshTopology* mesh, Matrix3X* scan, Matrix3X* scan_normals) {
		states.push_back(SingleState(states.size()-1,eval, _rayTracer, mesh, scan, scan_normals));
	}
	int SurfaceState::move(int stateId, Vector2 u, int f, Vector2 d, Vector2& out, int& fout, Scalar* distance,
            std::vector<Vector3>* path) {
		return states[stateId].move(u, f, d, out, fout, distance,path);
           
	}
	
	void SurfaceState::setCurrentModel(int stateId, Matrix3X& pts, Matrix3X& localPts, 
		const std::vector<int>& excludedVertices, bool rebuild_tree) {
        //std::cout << "Set Model " << stateId << std::endl;
		states[stateId].setCurrentModel(pts, localPts, excludedVertices, rebuild_tree);
	}
	
	void SurfaceState::alternatingUpdate(int stateId, int dataIndex, Vector2 u, int f, Vector2 d, Vector2& out, int& fout) {
		states[stateId].alternatingUpdate(dataIndex, u, f, d, out, fout);
	}

	bool SurfaceState::discreteUpdate(int stateId, int dataIndex, Vector2& out, int& fout) {
		return states[stateId].discreteUpdate(dataIndex, out, fout);
	}

    Vector3 SurfaceState::moveTangent(int stateId, Vector2 u, int f, Vector2 d) {
        return states[stateId].moveTangent(u,f,d);
    }
    
    void SurfaceState::disableRayTracing(int stateId) {
		states[stateId].disableRayTracing();
	}

	MeshTopology* SurfaceState::getTopology(int stateId) {
		return states[stateId].mesh;
	}
