#include "Camera.h"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <stdlib.h>
using json = nlohmann::json;

/*
Camera::Camera() {
        fx = 1498.2242623f;
    fy = 1498.2242623f;
    cx = 790.263706f;
    cy = 578.90334f;

        R(0,0) = 0.99924469;
        R(0,1) = -0.01071995;
        R(0,2) = 0.03735144;
        R(1,0) = -0.00488005;
        R(1,1) = -0.98820433;
        R(1,2) = -0.15306333;
        R(2,0) = 0.03855168;
        R(2,1) = 0.15276545;
        R(2,2) = -0.98751025;

        R.transposeInPlace();
        t[0] = -0.03609917;
        t[1] = 0.43416458;
        t[2] = 2.37101226;
}
*/
Camera::Camera(std::string configPath, int camId) {
  // std::cout << "Create Cam: " << configPath
  //	<< "View : " << camId << std::endl;
  if (Rs.size() == 0) {
    std::ifstream file(configPath);
    json j;
    file >> j;
    for (int i = 0; i < j["setup"].size(); i++) {
      {
        auto raw_joints = j["setup"][i]["R"];
        std::vector<Scalar> vec_joints;
        raw_joints.get_to(vec_joints);
        Rs.push_back(Matrix3(vec_joints.data()));
      }
      {
        auto raw_joints = j["setup"][i]["t"];
        std::vector<Scalar> vec_joints;
        raw_joints.get_to(vec_joints);
        ts.push_back(Vector3(vec_joints.data()));
      }
    }
  }
  fx = 2 * 800.0f;
  fy = 2 * 800.0f;
  cx = 2 * 395.131f;
  cy = 2 * 310.548f;

  R = Rs[camId];
  t = ts[camId];
  R.transposeInPlace();

  /*
R(0,0) = 0.99924469;
R(0,1) = -0.01071995;
R(0,2) = 0.03735144;
R(1,0) = -0.00488005;
R(1,1) = -0.98820433;
R(1,2) = -0.15306333;
R(2,0) = 0.03855168;
R(2,1) = 0.15276545;
R(2,2) = -0.98751025;
  t << -0.04783459,0.79178219,2.27646524+1.5;
  */
}

Vector3 Camera::getCameraOrientation() { return R.eulerAngles(2, 1, 0); }

Vector3 Camera::getCameraPosition() {
  // return R.transpose()*(-t);
  return t;
}

Vector3 Camera::computeRay(Vector2 pt) {
  std::cout << pt.transpose() << std::endl;
  Matrix3 C;
  C.setConstant(0);
  C(0, 0) = fx;
  C(0, 2) = cx;
  C(1, 1) = fy;
  C(1, 2) = cy;
  C(2, 2) = 1;

  Vector3 ph;
  ph[0] = pt[0];
  ph[1] = pt[1];
  ph[2] = 1;

  ph = C.inverse() * ph;

  // ph = ph - t;

  // ph = R.transpose() * ph;

  return R * ph + t;
}

Vector3 Camera::projectJointOnRay(Vector3 imagePt, Vector3 joint) {

  Vector3 x1 = getCameraPosition();
  Vector3 x2 = imagePt;
  Vector3 x3 = joint;

  Vector3 u = x2 - x1;
  u.normalize();
  Scalar d = -u.transpose() * x3;
  Scalar k = (-u.transpose() * x1 - d);
  Vector3 result = k * u + x1;
  return result;
}
