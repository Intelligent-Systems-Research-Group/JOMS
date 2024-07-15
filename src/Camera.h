#pragma once
#include "types.h"

class Camera {
private:
  std::vector<Matrix3> Rs;
  std::vector<Vector3> ts;

  Scalar fx;
  Scalar fy;
  Scalar cx;
  Scalar cy;

  Matrix3 R;
  Vector3 t;

public:
  Camera(std::string configPath, int camId = 0);
  Vector3 getCameraPosition();
  // Return the Camera Orientation in euler angles
  Vector3 getCameraOrientation();
  //
  Vector3 computeRay(Vector2 pt);
  // imagePt in world space
  Vector3 projectJointOnRay(Vector3 imagePt, Vector3 joint);
};
