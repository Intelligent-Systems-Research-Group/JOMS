#pragma once
#include "types.h"
#include <Eigen/Eigen>

//#include <Eigen\src\SparseCore\SparseUtil.h>

#include <iso646.h> //To define the words and, not, etc. as operators in Windows
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/stencilTableFactory.h>
#include <opensubdiv/osd/cpuEvaluator.h>
#include <opensubdiv/osd/cpuVertexBuffer.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/stencilTable.h>

#include "eigen_extras.h"
#include "MeshTopology.h"
#include "particles.h"
//#include "MarkerTrajectories.h"


using namespace OpenSubdiv;

#define MAX_NUM_W  28		//If using ENDCAP_BSPLINE_BASIS

// Vertex container implementation for OSD
struct OSD_Vertex {
  OSD_Vertex() { }

  void Clear(void * = 0) {
    point.setZero();
  }

  void AddWithWeight(OSD_Vertex const & src, float weight) {
    point += weight * src.point;
  }

  void SetPosition(float x, float y, float z) {
    point << x, y, z;
  }

  Vector3 point;
};

// A (u,v) point on a particular face of a facewise parametric surface
struct SurfacePoint {
  int face;
  Vector2 u;

  //  SurfacePoint() {}
  //  SurfacePoint(int face, Vector2 const& u) :face(face), u(u) {}
};

struct SubdivEvaluator {
  typedef Eigen::Triplet<Scalar> triplet_t;
  typedef std::vector<Eigen::Triplet<Scalar>> triplets_t;

  OpenSubdiv::Far::PatchTable *patchTable;
  OpenSubdiv::Far::PatchMap *patchMap;
  OpenSubdiv::Far::StencilTable const *cvstencils;

  size_t  nVertices;
  size_t  num_faces;
  size_t  nRefinerVertices;
  size_t  nLocalPoints;
  int ptexnum;

  std::vector<Vtr::Index> compressed_faces;

  Eigen::VectorXi vertsperface;

  mutable std::vector<OSD_Vertex> evaluation_verts_buffer;
  static const int maxlevel = 3;
  OpenSubdiv::Far::TopologyRefiner * refiner;
  Far::TopologyRefiner * refiner2;
  std::vector<SurfacePoint> controlVertexPtex;

  std::vector<STParticles::FaceInfo> adjacency;

  void generate_refined_mesh(Matrix3X const& vert_coords, int levels, MeshTopology* mesh_out, Matrix3X* verts_out);
  void init(Matrix3X const& vert_coords, Matrix3X* outS);
  //void move(const MeshTopology& mesh, const Matrix3X& vert_coords, const Matrix2X& du, Matrix2X* dir);
  int increment_u_crossing_edges(const MeshTopology& mesh, Matrix3X const& X, int face, const Vector2& u, const Vector2& du, int* new_face_out, Vector2* new_u_out, Vector2* dir, Scalar* distance, std::vector<Vector3>* path = NULL);
  int increment_u_crossing_edges_simple(const MeshTopology& mesh, Matrix3X const& X, int face, const Vector2& u, const Vector2& du, int* new_face_out, Vector2* new_u_out, Vector2* dir, Scalar* distance, std::vector<Vector3>* path = NULL);

  int increment_accurate(const MeshTopology& mesh, Matrix3X const& X, int face, const Vector2& u, const Vector2& du, int* new_face_out, Vector2* new_u_out, Vector2* dir, Scalar* distance, std::vector<Vector3>* path = NULL);
  int minimizeSurfacePoint(const MeshTopology& mesh, SurfacePoint p, const Matrix3X& X, const Matrix3X& target, SurfacePoint& q, Matrix3X* pathCloud);
  SubdivEvaluator(MeshTopology const& mesh);
  void evaluateSubdivSurface(Matrix3X const& vert_coords,
    std::vector<SurfacePoint> const& uv,
    Matrix3X* out_S,
    triplets_t* out_dSdX = 0,
    triplets_t* out_dSudX = 0,
    triplets_t* out_dSvdX = 0,
    Matrix3X* out_Su = 0,
    Matrix3X* out_Sv = 0,
    Matrix3X* out_Suu = 0,
    Matrix3X* out_Suv = 0,
    Matrix3X* out_Svv = 0,
    Matrix3X* out_N = 0,
    Matrix3X* out_Nu = 0,
    Matrix3X* out_Nv = 0,
    bool printDebug = false) const;

  void evaluatePtex(const Matrix3X& vert_coords,
	  Matrix3X* out_S,
	  triplets_t* out_dSdX=NULL,
	  triplets_t* out_dSudX=NULL,
	  triplets_t* out_dSvdX = NULL,
	  Matrix3X* out_Su = NULL,
	  Matrix3X* out_Sv = NULL,
      bool printDebug = false
  ) const; 


	//
	//triplets_t* out_dSudX,
	//	  triplets_t* out_dSvdX,
	// Matrix3X* out_Su,
	//  Matrix3X* out_Sv,
	//  Matrix3X* out_Suu,
	//  Matrix3X* out_Suv,
	//  Matrix3X* out_Svv,
	// Matrix3X* out_N,  Matrix3X* out_Nu,  Matrix3X* out_Nv


  const OpenSubdiv::Far::StencilTable*  getLocalStencilTable();
  OpenSubdiv::Far::StencilTable const * GetLocalPointStencilTable();
  Matrix3X localPoints(const Matrix3X& base);
  Matrix3X evalPointsCustom(const Matrix3X& base);
  Matrix3X evalPointsCustom(const Matrix3X& base, const std::vector<SurfacePoint>& surfacePoints);
  Matrix3X evalPointsPartial(const Matrix3X& base);  
  Matrix3X evalPointsPartial(const Matrix3X& base, const std::vector<SurfacePoint>& surfacePoints);
  void evalPointsPartial(const Matrix3X& localPts, const std::vector<SurfacePoint>& surfacePoints,
	Matrix3X* out_S, Matrix3X* out_Su=NULL, Matrix3X* out_Sv=NULL);
  
  Vector2 unnormalize(int patch_idx, Vector2 u);
  Vector2   normalize(int face_idx, Vector2 u, int* patch_idx);
  std::vector<int> getBase();

  SubdivEvaluator(SubdivEvaluator const& that);
  SubdivEvaluator& operator=(SubdivEvaluator const& that);

  ~SubdivEvaluator();

 void createPatches(const Matrix3X& localPts,
    MeshTopology* mesh, Matrix3X* verts);
};
