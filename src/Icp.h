#pragma once

#include "types.h"
#include <memory>
#include <nanoflann.hpp>

class Icp {
public:
  Icp(Matrix3X cloud);
  int findClosestPoint(Vector3 query);

private:
  typedef nanoflann::KDTreeEigenMatrixAdaptor<MatrixX> my_kd_tree_t;
  std::unique_ptr<my_kd_tree_t> mat_index;
  MatrixX data;
};
