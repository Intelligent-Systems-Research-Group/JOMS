#pragma once
#include "types.h"
#include <string>

class Pose2d {
private:
  Matrix3X joints;
  std::string convertJointName(std::string joint);
  std::map<std::string, std::string> skel_to_openpose;
  std::map<std::string, unsigned int> openpose_to_index;
  void init(bool base25);

public:
  Pose2d(std::string surfaceMapPath);
  void read(std::string path);
  Vector3 getJoint(std::string joint);
  std::vector<size_t> getSurfaceIndices();
  Vector3 getSurface(size_t vertexIdx);
};
