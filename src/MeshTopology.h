#pragma once

#include "types.h"
#include <set>

struct texcoord {
  int faceId;
  double s;
  double t;
};
struct MeshTopology {
  Eigen::Array<int, 5, Eigen::Dynamic> quads;
  IMatrix4X face_adj;

  std::vector<int> compressed_faces;
  std::vector<size_t> nverts_per_face;
  std::vector<texcoord> texcoords;

  size_t num_vertices;
  size_t num_faces() const { return quads.cols(); }

  void update_adjacencies();
  void createEdgeList(IMatrix2X &edges, IMatrix4X *faceAdj = NULL);
  void vertexAdjacency(std::map<size_t, std::set<size_t>> &adjacency,
                       IMatrix4X *faceAdj = NULL);
  void vertexToSurfaceCoordinate(int vidx, int &fidx, double &u, double &v);
  void initVertexToSurfaceCoordinate(int vidx, int &fidx, double &u, double &v);
  std::map<std::pair<size_t, size_t>, double> calculateLaplaceWeights();
  IMatrix4X findClosest(const Matrix3X &verts,
                        const std::vector<std::vector<int>> &vgroups,
                        const SparseMatrix &, Matrix4X *distance = NULL);
  void computeMaskFaces(const std::vector<int> &invMaskedVertices,
                        std::vector<bool> *outMaskedFaces);
};

void makeCube(MeshTopology *mesh, Matrix3X *verts);
void makeCube2(MeshTopology *mesh, Matrix3X *verts);
void makeTriangle(MeshTopology *mesh, Matrix3X *verts);
void makePath(MeshTopology *mesh, Matrix3X *verts);
bool loadObj(std::string path, MeshTopology *mesh, Matrix3X *verts,
             Matrix3X *normals = NULL, double scale = 1.0, bool mirror = false,
             bool rot = false);
bool saveObj(std::string path, const MeshTopology *mesh, const Matrix3X *verts);
void dumpCloud(std::string name, int iter, const Matrix3X &cloud,
               const Matrix3X *normals = NULL);
bool readCloud(std::string name, Matrix3X *cloud, Matrix3X *normals,
               double scale, bool mirror, bool rot = false);
bool readMatrix(std::string name, MatrixX *matrix);
bool readVector(std::string name, std::vector<bool> *result);
void filter(Matrix3X *data, Matrix3X *normal, std::vector<bool> mask);
void filter(Matrix3X *data, Matrix3X *normal);
void filterGroundPlane(Matrix3X *data, Matrix3X *normal, int axis,
                       Scalar groundPlane);
void fill(Matrix3X *data, Matrix3X *normal, int n);
void findSymmetric(const Matrix3X &verts, std::vector<int> *sym);
void filter(Matrix3X *data, Matrix3X *normal, const Matrix3X &model);
void corruptNonUniform(Matrix3X *data, Matrix3X *normal);
std::vector<int> greedyFarthestSampling(const Matrix3X &data,
                                        const Matrix3X &normal, int k,
                                        Matrix3X *out, Matrix3X *nout);
void quiver(std::string path, const Matrix3X &X, const Matrix3X &Y);
