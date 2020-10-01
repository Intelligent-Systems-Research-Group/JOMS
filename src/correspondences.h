#pragma once

#include <nanoflann.hpp>
#include "types.h"
#include "MeshTopology.h"
#include "SubdivEvaluator.h"
#include <memory>
#include "raytracer.h"

class ISurfaceState {
public:
	virtual int move(int stateId, Vector2 u, int f, Vector2 d, Vector2& out, int& fout, Scalar* dist,
        std::vector<Vector3>* path = NULL) = 0;
    virtual Vector3 moveTangent(int stateId, Vector2 u, int f, Vector2 d) = 0;
    virtual void alternatingUpdate(int stateId, int dataIndex,Vector2 u, int f, Vector2 d, Vector2& out, int& fout) = 0;
	virtual bool discreteUpdate(int stateId, int dataIndex, Vector2& out, int& fout) = 0;
	virtual void setCurrentModel(int stateId, Matrix3X& pts, Matrix3X& localPts, bool rebuild_tree) = 0;
	virtual MeshTopology* getTopology(int stateId) = 0;
	virtual void disableRayTracing(int stateId) = 0;
};


class SurfaceState : public ISurfaceState {
	class SingleState {
		typedef nanoflann::KDTreeEigenMatrixAdaptor<MatrixX>  my_kd_tree_t;
	private:
		Matrix3X pts;
        Matrix3X localPts;
        MeshTopology resultTopology;
		Matrix3X subdVertices;
		
		MatrixX data;
		Matrix3X* scan;
		Matrix3X* scan_normals;
		SubdivEvaluator* eval;
		std::unique_ptr<my_kd_tree_t> mat_index;
		int npts;
		const Scalar nw = 0.00; //001
		const Scalar w = 1;
        int index;
		std::shared_ptr<RayTracer> rayTracer;
		bool shooting;
	public:
		MeshTopology* mesh;
		SingleState(int index, SubdivEvaluator* eval, std::shared_ptr<RayTracer> _rayTracer, MeshTopology* mesh, Matrix3X* scan, Matrix3X* scan_normals);
		virtual int move(Vector2 u, int f, Vector2 d, Vector2& out, int& fout, Scalar* distance,
            std::vector<Vector3>* path = NULL);
        virtual Vector3 moveTangent(Vector2 u, int f, Vector2 d);
		virtual void readModelPoints(Matrix3X& outS, Matrix3X& n);
		virtual Scalar evaluateCurrent(int dataIndex, Vector2 u, int face);
		virtual void setCurrentModel(Matrix3X& pts, Matrix3X& localPts, bool rebuild_tree);
		virtual void alternatingUpdate(int dataIndex, Vector2 u, int f, Vector2 d, Vector2& out, int& fout);
		virtual bool discreteUpdate(int dataIndex, Vector2& out, int& fout);
		virtual void disableRayTracing();
    };
	
    std::vector<SingleState> states;
public:
	SurfaceState();
	void addData(SubdivEvaluator* eval, std::shared_ptr<RayTracer> _rayTracer, MeshTopology* mesh, Matrix3X* scan, Matrix3X* scan_normals);
	virtual int move(int stateId, Vector2 u, int f, Vector2 d, Vector2& out, int& fout, Scalar* distance,
        std::vector<Vector3>* path = NULL);
    virtual Vector3 moveTangent(int stateId, Vector2 u, int f, Vector2 d);
	virtual void setCurrentModel(int stateId, Matrix3X& pts,Matrix3X& localPts,  bool rebuild_tree);
	virtual void alternatingUpdate(int stateId, int dataIndex,Vector2 u, int f, Vector2 d, Vector2& out, int& fout);
	virtual bool discreteUpdate(int stateId, int dataIndex, Vector2& out, int& fout);
	virtual void disableRayTracing(int stateId);
	virtual MeshTopology* getTopology(int stateId);
};

