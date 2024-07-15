#include "VertexGroupLoader.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

void loadVertexGroups(std::string path,
                      std::map<std::string, std::vector<size_t>> *vertex_groups,
                      std::map<std::pair<std::string, std::string>,
                               std::vector<size_t>> *border_groups,
                      Weights *weights) {

  std::ifstream infile(path);
  assert(!infile.fail());
  std::string line;
  int num_verts = 0;
  int num_faces = 0;
  int vert_idx = 0;
  int face_idx = 0;

  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    size_t index;
    if (!(iss >> index)) {
      std::cout << "error" << std::endl;
      assert(0);
    }
    std::string group_name;
    std::string border_names[2];
    size_t i = 0;

    while ((iss >> group_name)) {
      assert(i < 2);
      // std::cout << group_name << std::endl;
      (*vertex_groups)[group_name].push_back(index);
      border_names[i] = group_name;
      i++;
    }
    if (i == 1) {
      if (weights)
        (*weights)[JointVertexPair(border_names[0], index)] = 1.0;
    }
    if (i == 2) {
      if (weights) {
        if (border_names[0] != "m_root" && border_names[1] != "m_root") {
          (*weights)[JointVertexPair(border_names[0], index)] = 0.5;
          (*weights)[JointVertexPair(border_names[1], index)] = 0.5;
        } else if (border_names[0] == "m_root" && border_names[1] != "m_root") {
          (*weights)[JointVertexPair(border_names[1], index)] = 1;
        } else if (border_names[0] != "m_root" && border_names[1] == "m_root") {
          (*weights)[JointVertexPair(border_names[0], index)] = 1;
        } else {
          std::cout << "ERROR" << std::endl;
        }
      }
      string first =
          border_names[0] < border_names[1] ? border_names[0] : border_names[1];
      string second =
          border_names[0] < border_names[1] ? border_names[1] : border_names[0];
      if (border_groups)
        (*border_groups)[std::pair<std::string, std::string>(first, second)]
            .push_back(index);
    }
  }
}

std::map<std::pair<std::string, std::string>, Vector3> calculateJoints(
    const std::map<std::pair<std::string, std::string>, std::vector<size_t>>
        &vertex_groups,
    Matrix3X *vertices) {
  std::map<std::pair<std::string, std::string>, Vector3> joints;
  std::map<std::pair<std::string, std::string>,
           std::vector<size_t>>::const_iterator it;
  for (it = vertex_groups.begin(); it != vertex_groups.end(); ++it) {
    const std::vector<size_t> &indices = it->second;
    Vector3 joint_pos;
    joint_pos.setZero();
    for (size_t i = 0; i < indices.size(); i++) {
      Vector3 vec = vertices->col(indices[i]);
      joint_pos += vec;
    }
    joint_pos /= indices.size();
    joints[it->first] = joint_pos;
  }
  return joints;
}

std::map<std::pair<std::string, std::string>, VectorX>
calculateInitialRegressor(const std::map<std::pair<std::string, std::string>,
                                         std::vector<size_t>> &vertex_groups,
                          Matrix3X *vertices) {
  std::map<std::pair<std::string, std::string>, VectorX> jointsRegressor;
  std::map<std::pair<std::string, std::string>,
           std::vector<size_t>>::const_iterator it;
  for (it = vertex_groups.begin(); it != vertex_groups.end(); ++it) {
    const std::vector<size_t> &indices = it->second;
    VectorX jointRegressorPart(vertices->size());
    jointRegressorPart.setZero();
    double n = indices.size();
    for (size_t i = 0; i < indices.size(); i++) {
      int idx = indices[i];
      jointRegressorPart[i] = 1.0 / n;
      jointRegressorPart[i + 1 * vertices->cols()] = 1.0 / n;
      jointRegressorPart[i + 2 * vertices->cols()] = 1.0 / n;
    }
    jointsRegressor[it->first] = jointRegressorPart;
  }
  return jointsRegressor;
}

VectorX loadVector(std::string path) {

  int size = 0;
  std::string line;
  float x;
  {
    std::ifstream infile(path);
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (!(iss >> x)) {
        std::cout << "error" << std::endl;
        assert(0);
      }
      size++;
    }
  }

  VectorX y(size);
  y.setConstant(-1);

  {
    std::ifstream infile(path);
    int i = 0;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (i >= size || !(iss >> x)) {
        std::cout << "error" << std::endl;
        assert(0);
      }
      y[i++] = x;
    }
  }
  return y;
}
