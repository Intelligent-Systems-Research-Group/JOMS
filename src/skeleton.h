#pragma once

#include "types.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <map>
#include <string>
#include <vector>

Vector3 vector4to3(const Vector4 &vec);
Vector4 vector3to4(const Vector3 &vec);

void angle_axis_to_rotation_matrix(const Vector3 &angle_axis, Matrix3 *R);

Matrix4 build_rigid_transform(const Matrix3 &R, const Vector3 &t);
Matrix4 rigid_inverse(const Matrix4 &T);

void rod2Mat(const Matrix3X &rod, std::vector<Matrix3> &mats);
void angle2euler(const Matrix3X &angle, Matrix3X &euler);

template <typename U, typename V, typename W>
void mul4x4(W &w, const U &u, const V &v) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      W(i, j) = U(i, 0) * V(0, j) + U(i, 1) * V(1, j) + U(i, 2) * V(2, j) +
                U(i, 3) * V(3, j);
    }
  }
}

namespace hyp {

class Skeleton {
private:
  std::vector<std::string> part_names;
  std::vector<std::string> rigid;
  std::vector<std::string> halfrigid;
  std::map<std::string, std::string> ancestors;
  std::map<std::string, std::vector<std::string>> influence;
  std::map<std::string, std::string> symmetry;

  void set_joint_location(std::string from, std::string to, Vector3 position,
                          std::vector<Vector3> &joint_positions);

public:
  Skeleton(std::string path);
  void calculate_kinematic_chain(std::vector<Vector3> &joint_positions,
                                 std::vector<Matrix3> &rotations,
                                 std::vector<Matrix4> &result,
                                 Matrix3X resJoints = Matrix3X());
  void add(std::map<std::string, std::vector<size_t>> &vertex_groups,
           std::vector<bool> *output);
  bool isRigid(size_t index, bool fullRigid);
  size_t getParentJoint(size_t index);
  std::pair<std::string, std::string> getJoint(size_t index);
  int getJointIndex(std::string from, std::string to);
  int getJointMirror(size_t index);
  size_t getJointCount();
  std::vector<Scalar> extractVertexWeights(const Weights &weights,
                                           size_t index);
  SparseMatrix extractWeightMatrix(const Weights &weights, size_t n);
  void setPose(
      std::vector<Matrix3> &vec_theta,
      std::map<std::pair<std::string, std::string>, Vector3> &joint_locations,
      std::vector<Matrix4> &absolutes);
  Vector3 skin_vertex(std::vector<Matrix4> &absolutes, const VectorX &weights,
                      const Vector3 &position);
  Vector3 rigidOrientation(int j);
  Vector3 makeFist(int jointIdx);
};
} // namespace hyp
