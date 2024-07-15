#include "skeleton.h"

#include "rotations.h"
#include "types.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

Vector3 vector4to3(const Vector4 &vec) {
  Vector3 out;
  out[0] = vec[0];
  out[1] = vec[1];
  out[2] = vec[2];
  return out;
}
Vector4 vector3to4(const Vector3 &vec) {
  Vector4 out;
  out[0] = vec[0];
  out[1] = vec[1];
  out[2] = vec[2];
  out[3] = 1;
  return out;
}
/*
void angle_axis_to_rotation_matrix(
        const DVector3& angle_axis,
        DMatrix3 *R)
{
        Diffable norm = angle_axis.norm();

        //if (norm < .0001)
        //{
        //R->setIdentity();
        //return;
        //}


        Diffable x = angle_axis[0] / norm;
        Diffable y = angle_axis[1] / norm;
        Diffable z = angle_axis[2] / norm;

        Diffable s = sin(norm);
        Diffable c = cos(norm);

        *R << x*x + (1 - x*x)*c, x*y*(1 - c) - z*s, x*z*(1 - c) + y*s,
                x*y*(1 - c) + z*s, y*y + (1 - y*y)*c, y*z*(1 - c) - x*s,
                x*z*(1 - c) - y*s, z*y*(1 - c) + x*s, z*z + (1 - z*z)*c;
}
*/

Matrix4 build_rigid_transform(const Matrix3 &R, const Vector3 &t) {
  Matrix4 result;
  result.setIdentity();
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      result(i, j) = R(i, j);
    }
  }
  for (size_t i = 0; i < 3; i++) {
    result(i, 3) = t[i];
  }
  return result;
}

Matrix4 rigid_inverse(const Matrix4 &T) {
  Matrix3 R = T.block<3, 3>(0, 0);
  Vector3 t = T.block<3, 1>(0, 3);
  Matrix4 result;
  result.setIdentity();
  result.block<3, 3>(0, 0) = R.transpose();
  result.block<3, 1>(0, 3) = -R.transpose() * t;
  return result;
}

namespace hyp {

void Skeleton::calculate_kinematic_chain(std::vector<Vector3> &joint_positions,
                                         std::vector<Matrix3> &rotations,
                                         std::vector<Matrix4> &result,
                                         Matrix3X resJoints) {
  // std::cout << "calc kinematic chain" << std::endl;
  bool isAbsolute = resJoints.cols() != 0;
  int jointCount = joint_positions.size();

  std::vector<Matrix4> relatives(jointCount);
  std::vector<Matrix4> fixedJointTransforms(jointCount);

  for (size_t i = 0; i < jointCount; i++) {
    Matrix3 &rotation = rotations[i];
    Matrix4 relative;
    Vector3 joint_rel;
    if (i == 0)
      joint_rel = joint_positions[0];
    else {
      std::string ancestor = this->ancestors[part_names[i]];
      int j = std::find(part_names.begin(), part_names.end(), ancestor) -
              part_names.begin();
      joint_rel = joint_positions[i] - joint_positions[j];
    }
    if (isAbsolute) {
      // Vector3 base = resJoints.block<3,1>(0,0);
      Vector3 current = resJoints.block<3, 1>(0, i);
      // auto rt = (i == 0) ? base : current - base;
      relative = build_rigid_transform(rotation, current);
    } else {
      relative = build_rigid_transform(rotation, joint_rel);
    }
    relatives[i] = relative;
    Matrix4 fixedJointTransform;
    fixedJointTransform.setIdentity();
    for (size_t k = 0; k < 3; k++) {
      // fixedJointTransform(k, 3) = joint_positions[i][k];
      fixedJointTransform(k, 3) = joint_rel[k];
      // fixedJointTransform(k, 3) =
      //	joint_positions[i][k];
      //	joint_positions[0][k]
    }
    fixedJointTransforms[i] = fixedJointTransform;
  }

  std::vector<Matrix4> absolutes(jointCount);
  std::vector<Matrix4> absolutesFixed(jointCount);

  absolutes[0] = relatives[0];
  absolutesFixed[0] = fixedJointTransforms[0];
  // std::cout << "survived relative" << std::endl;
  for (size_t i = 1; i < jointCount; i++) {
    // std::cout << "relid " << i << std::endl;
    std::string ancestor = this->ancestors[part_names[i]];
    int j = std::find(part_names.begin(), part_names.end(), ancestor) -
            part_names.begin();
    if (isAbsolute) {
      // absolutes[i].block<3,3>(0,0) = relatives[i].block<3,3>(0,0);
      // absolutes[i].block<3,1>(0,3) = absolutes[j].block<3,3>(0,0) *
      // relatives[i].block<3,1>(0,3)
      //	+ absolutes[j].block<3,1>(0,3);
      // absolutes[i](3,3) = 1;
      absolutes[i] = relatives[i];
    } else {
      absolutes[i] = absolutes[j] * relatives[i];
    }
    absolutesFixed[i] = absolutesFixed[j] * fixedJointTransforms[i];
  }

  for (size_t i = 0; i < jointCount; i++) {
    // std::cout << absolutes[i].col(3).transpose() << std::endl;
    result[i] = absolutes[i] * rigid_inverse(absolutesFixed[i]);
  }
}

void Skeleton::set_joint_location(std::string from, std::string to,
                                  Vector3 position,
                                  std::vector<Vector3> &joint_positions) {
  assert(ancestors[to] == from || ancestors[from] == to);
  std::string joint_name = ancestors[to] == from ? to : from;
  int index = std::find(part_names.begin(), part_names.end(), joint_name) -
              part_names.begin();
  joint_positions[index] = position;
  // std::cout << index <<  ":" << joint_name << " -> " <<
  // joint_positions[index].transpose() << std::endl;
}
Skeleton::Skeleton(std::string path) {
  std::ifstream file(path);
  json f;
  file >> f;
  f["part_names"].get_to(part_names);
  f["ancestors"].get_to(ancestors);
  f["influence"].get_to(influence);
  f["rigid"].get_to(rigid);
  f["halfrigid"].get_to(halfrigid);
  f["symmetry"].get_to(symmetry);
  // std::cout << part_names << std::endl;
  // std::cout << ancestors << std::endl;
  // std::cout << influence << std::endl;
  // std::cout << rigid << std::endl;
  // exit(1);
}

void Skeleton::add(std::map<std::string, std::vector<size_t>> &vertex_groups,
                   std::vector<bool> *output) {
  for (std::string name : rigid)
    for (size_t i : vertex_groups[name]) {
      (*output)[i] = true;
    }
}

bool Skeleton::isRigid(size_t index, bool fullRigid) {
  std::string name = part_names[index];
  std::vector<std::string> &temp = fullRigid ? rigid : halfrigid;
  auto it = std::find(temp.begin(), temp.end(), name);
  return it != temp.end();
}

int Skeleton::getJointMirror(size_t index) {
  std::string name = part_names[index];
  std::string mname = symmetry[name];
  auto it = std::find(part_names.begin(), part_names.end(), mname);
  assert(it != part_names.end());
  int ret_idx = it - part_names.begin();
  return ret_idx;
}

std::pair<std::string, std::string> Skeleton::getJoint(size_t index) {
  return std::pair<std::string, std::string>(part_names[index],
                                             ancestors[part_names[index]]);
}

size_t Skeleton::getParentJoint(size_t index) {
  auto parentName = ancestors[part_names[index]];
  return std::find(part_names.begin(), part_names.end(), parentName) -
         part_names.begin();
}

int Skeleton::getJointIndex(std::string from, std::string to) {
  std::cout << "Ancestor size: " << ancestors.size() << std::endl;
  std::map<std::string, std::string>::iterator itera = ancestors.find(from);
  std::map<std::string, std::string>::iterator iterb = ancestors.find(to);
  std::map<std::string, std::string>::iterator begin = ancestors.begin();
  if (itera != ancestors.end() && itera->second == to) {
    for (size_t i = 0; i < part_names.size(); i++) {
      if (itera->first == part_names[i])
        return i;
    }
  } else if (iterb != ancestors.end() && iterb->second == from) {
    for (size_t i = 0; i < part_names.size(); i++) {
      if (iterb->first == part_names[i])
        return i;
    }
  }
  std::cout << "Test:" << std::endl;
  std::cout << iterb->first << std::endl;
  std::cout << (itera != ancestors.end() && itera->second == to) << std::endl;
  std::cout << (iterb != ancestors.end() && iterb->second == from) << std::endl;
  std::cout << from << " -> " << to << std::endl;
  std::cout << "Error!" << std::endl;
  return -1;
}

size_t Skeleton::getJointCount() { return part_names.size(); }
/*
void print() {
for (auto i = joint_positions.begin(); i != joint_positions.end(); i++) {
std::cout << i->transpose() << std::endl;
}
}
*/

std::vector<Scalar> Skeleton::extractVertexWeights(const Weights &weights,
                                                   size_t index) {
  std::vector<Scalar> result(part_names.size());
  int hits = 0;
  for (size_t i = 0; i < part_names.size(); i++) {
    JointVertexPair pair(part_names[i], index);
    Weights::const_iterator iter = weights.find(pair);
    if (iter != weights.end()) {
      // std::cout << pair.first << " " << pair.second << " " << iter->second <<
      // std::endl;
      result[i] = iter->second;
      hits++;
    }
  }
  std::cout << index << " " << hits << std::endl;
  assert(hits > 0);
  return result;
}

SparseMatrix Skeleton::extractWeightMatrix(const Weights &weights, size_t n) {
  SparseMatrix result;
  Weights::const_iterator it;
  for (it = weights.begin(); it != weights.end(); ++it) {
    size_t index = it->first.second;
    std::string name = it->first.first;
    auto iter = std::find(part_names.begin(), part_names.end(), name);
    assert(iter != part_names.end());
    size_t jointIndex = iter - part_names.begin();
    result[IndexPair(index, jointIndex)] = it->second;
  }
  return result;
}

void Skeleton::setPose(
    std::vector<Matrix3> &vec_theta,
    std::map<std::pair<std::string, std::string>, Vector3> &joint_locations,
    std::vector<Matrix4> &absolutes) {
  // std::vector<DMatrix4> absolutes(part_names.size());
  std::vector<Vector3> joint_positions(part_names.size());
  for (auto it_joint = joint_locations.begin();
       it_joint != joint_locations.end(); ++it_joint) {
    // std::cout << it_joint->first.first << " <-> " << it_joint->first.second
    // << " : " << it_joint->second.transpose() << endl;
    set_joint_location(it_joint->first.first, it_joint->first.second,
                       it_joint->second, joint_positions);
  }
  calculate_kinematic_chain(joint_positions, vec_theta, absolutes);
}
Vector3 Skeleton::skin_vertex(std::vector<Matrix4> &absolutes,
                              const VectorX &weights, const Vector3 &position) {
  // assert(vec_theta.size == part_names.size());
  // assert(weights.size == part_names.size());
  Vector3 result;

  result.setZero();
  Vector4 homo = vector3to4(position);
  for (size_t i = 0; i < absolutes.size(); i++) {
    // if(weights[i] > 1E-15) {
    Vector4 transformed = absolutes[i] * homo;
    // result += weights[i] * transformed;
    result[0] += weights[i] * transformed[0];
    result[1] += weights[i] * transformed[1];
    result[2] += weights[i] * transformed[2];
    //}
  }
  return result;
}

Vector3 Skeleton::rigidOrientation(int j) {
  Vector3 rot;
  rot.setZero();
  Matrix3X angle(3, 1);
  Matrix3X euler(3, 1);
  euler.setZero();
  return euler;
  Matrix3 R, R2, Q;
  Matrix3X angle2(3, 1);

  switch (j) {
  case 0:
    angle.col(0) << 90 * 3.141 / 180.0, 0, 0;
    angle_axis_to_rotation_matrix(angle, &R);
    angle2.col(0) << 0, 0, 45 * 3.141 / 180.0;
    angle_axis_to_rotation_matrix(angle2, &R2);
    Q = R2 * R;
    euler.col(0) = Q.eulerAngles(2, 1, 0);
    // angle2euler(angle,euler);
    break;
  case 3:
    angle.col(0) << 0, 0, -70 * 3.141 / 180.0;
    angle2euler(angle, euler);
    break;
  case 4:
    angle.col(0) << -45 * 3.141 / 180.0, 0, 0;
    angle2euler(angle, euler);
    break;
  case 5:
    angle.col(0) << 45 * 3.141 / 180.0, 0, 0;
    angle2euler(angle, euler);
    break;
  case 14:
    angle.col(0) << 0, 0, 70 * 3.141 / 180.0;
    angle2euler(angle, euler);
    break;
  case 15:
    angle.col(0) << -45 * 3.141 / 180.0, 0, 0;
    angle2euler(angle, euler);
    break;
  case 16:
    angle.col(0) << 45 * 3.141 / 180.0, 0, 0;
    angle2euler(angle, euler);
    break;
  }
  rot[0] = euler(0, 0);
  rot[1] = euler(1, 0);
  rot[2] = euler(2, 0);
  // rot.setZero();
  return rot;
  // if(j != 0) {
  //}
  // return Vector3(0,0,3.141f/2);
}
Vector3 Skeleton::makeFist(int j) {

  Vector3 rot;
  rot.setZero();
  Matrix3X angle(3, 1);
  Matrix3X euler(3, 1);
  euler.setZero();

  switch (j) {
  case 3:
    angle.col(0) << 0, 0, -80 * 3.141 / 180.0;
    angle2euler(angle, euler);
    break;
  case 14:
    angle.col(0) << 0, 0, 80 * 3.141 / 180.0;
    angle2euler(angle, euler);
    break;

  case 6:
    angle.col(0) << 60 * 3.141 / 180.0, 0, 0;
    angle2euler(angle, euler);
    break;
  case 7:
    angle.col(0) << 120 * 3.141 / 180.0, 0, 0;
    angle2euler(angle, euler);
    break;
  // case 8:
  case 9:
  case 10:
  case 11:
  case 12:
  case 13:
    angle.col(0) << 0, 0, -90 * 3.141 / 180.0;
    angle2euler(angle, euler);
    break;

  case 17:
    angle.col(0) << 60 * 3.141 / 180.0, 0, 0;
    angle2euler(angle, euler);
    break;
  case 18:
    angle.col(0) << 120 * 3.141 / 180.0, 0, 0;
    angle2euler(angle, euler);
    break;
  // case 19:
  case 20:
  case 21:
  case 22:
  case 23:
  case 24:
    angle.col(0) << 0, 0, 90 * 3.141 / 180.0;
    angle2euler(angle, euler);
    break;
  default:
    break;
  }

  rot[0] = euler(0, 0);
  rot[1] = euler(1, 0);
  rot[2] = euler(2, 0);
  return rot;
}
} // namespace hyp

void rod2Mat(const Matrix3X &rod, std::vector<Matrix3> &mats) {
  mats.resize(rod.cols());
  for (int i = 0; i < rod.cols(); i++) {
    VectorX axis = rod.col(i);
    euler_to_rotation_matrix(axis, &mats[i]);
  }
}

void angle2euler(const Matrix3X &angle, Matrix3X &euler) {
  for (int i = 0; i < angle.cols(); i++) {
    Matrix3 R;
    angle_axis_to_rotation_matrix(angle.col(i), &R);
    Vector3 e = R.eulerAngles(2, 1, 0);

    euler.col(i) = e;
  }
}
