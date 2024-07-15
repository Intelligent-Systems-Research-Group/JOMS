#include "model.h"
#include <Eigen/Dense>
#include <H5Cpp.h>

#include "eigen3-hdf5.hpp"
#include <random>

void Model::write(std::string path) {
  H5::H5File file(path, H5F_ACC_TRUNC);
  EigenHDF5::save(file, "mean", mean);
  for (int i = 0; i < pcs.size(); i++) {
    EigenHDF5::save(file, "pc" + std::to_string(i), pcs[i]);
  }
  EigenHDF5::save(file, "joints", joints);
  for (int i = 0; i < jointpcs.size(); i++) {
    EigenHDF5::save(file, "jointpc" + std::to_string(i), jointpcs[i]);
  }
  EigenHDF5::save(file, "weights", weights);
  for (int i = 0; i < deform.size(); i++) {
    EigenHDF5::save(file, "deform" + std::to_string(i), deform[i]);
  }
  EigenHDF5::save(file, "meanPose", meanPose);

  int npersons = persons.size();
  int npcs = pcs.size();
  MatrixX betas;
  betas.resize(npersons, npcs);
  for (int i = 0; i < npersons; i++) {
    for (int c = 0; c < npcs; c++) {
      betas(i, c) = persons[i].coeff[c];
    }
  }
  EigenHDF5::save(file, "betas", betas);

  int njoints = scans[0].posebase.cols();
  Matrix3X ts;
  ts.resize(3, scans.size());
  MatrixX Rs;
  Rs.resize(3 * njoints, scans.size());
  for (int i = 0; i < scans.size(); i++) {
    ts.col(i) = scans[i].trans;
    for (int j = 0; j < njoints; j++) {
      Rs(3 * j + 0, i) = scans[i].posebase(0, j);
      Rs(3 * j + 1, i) = scans[i].posebase(1, j);
      Rs(3 * j + 2, i) = scans[i].posebase(2, j);
    }
  }
  EigenHDF5::save(file, "ts", ts);
  EigenHDF5::save(file, "Rs", Rs);
}

void Model::read(std::string path) {
  H5::H5File file(path, H5F_ACC_RDONLY);
  EigenHDF5::load(file, "mean", mean);
  for (int i = 0; i < pcs.size(); i++) {
    EigenHDF5::load(file, "pc" + std::to_string(i), pcs[i]);
  }
  EigenHDF5::load(file, "joints", joints);
  for (int i = 0; i < jointpcs.size(); i++) {
    EigenHDF5::load(file, "jointpc" + std::to_string(i), jointpcs[i]);
  }
  EigenHDF5::load(file, "weights", weights);
  for (int i = 0; i < deform.size(); i++) {
    EigenHDF5::load(file, "deform" + std::to_string(i), deform[i]);
  }
  EigenHDF5::load(file, "meanPose", meanPose);
#ifndef VERBOSE_DUMP
  return;
#endif
  MatrixX betas;
  Matrix3X ts;
  MatrixX Rs;
  std::cout << "Before Read betas" << std::endl;
  EigenHDF5::load(file, "betas", betas);
  std::cout << "Before Read ts" << std::endl;
  EigenHDF5::load(file, "ts", ts);
  std::cout << "Before Read Rs" << std::endl;
  EigenHDF5::load(file, "Rs", Rs);
  assert(persons.size() == betas.rows() && "person count != beta count");
  assert(scans.size() == ts.cols() && "scan count != ts count");
  assert(scans.size() == Rs.cols() && "scan count != Rs count");
  assert(scans[0].posebase.cols() == Rs.rows() / 3 &&
         "pose joints != Rs.rows()/3");
  assert(betas.cols() >= persons[0].coeff.size());
  int betaDims = persons[0].coeff.size();

  for (int i = 0; i < persons.size(); i++) {
    persons[i].coeff = betas.row(i).segment(0, betaDims);
  }
  for (int i = 0; i < scans.size(); i++) {
    scans[i].trans = ts.col(i);
    for (int j = 0; j < Rs.rows() / 3; j++) {
      scans[i].posebase(0, j) = Rs(3 * j + 0, i);
      scans[i].posebase(1, j) = Rs(3 * j + 1, i);
      scans[i].posebase(2, j) = Rs(3 * j + 2, i);
    }
  }
}

void Model::init(const Matrix3X &_mean, const Matrix3X &_joints,
                 const MatrixX &_weights, int ncomp, int ndefcomp, int npersons,
                 int nscans) {

  mean = _mean;
  joints = _joints;
  weights = _weights;

  pcs.resize(ncomp);
  jointpcs.resize(ncomp);
  deform.resize(ndefcomp);

  meanPose.resize(3, joints.cols() - 1);
  meanPose.setConstant(0);

  for (int i = 0; i < pcs.size(); i++) {
    pcs[i].resize(3, mean.cols());
    pcs[i].setConstant(0);
  }
  for (int i = 0; i < jointpcs.size(); i++) {
    jointpcs[i].resize(3, joints.cols());
    jointpcs[i].setConstant(0);
  }
  for (int i = 0; i < deform.size(); i++) {
    deform[i].resize(3, mean.cols());
    deform[i].setConstant(0);
  }

  persons.resize(npersons);
  scans.resize(nscans);
  std::normal_distribution<double> distribution(0, 0.001);
  std::default_random_engine generator;

  for (int i = 0; i < persons.size(); i++) {
    persons[i].coeff.resize(ncomp);
    for (int j = 0; j < ncomp; j++) {
      persons[i].coeff[j] = distribution(generator);
    }
  }
  for (int i = 0; i < scans.size(); i++) {
    scans[i].trans.resize(3, 1);
    scans[i].trans.setConstant(0);
    scans[i].posebase.resize(3, joints.cols());
    scans[i].posebase.setConstant(0);
  }
}

std::ostream &operator<<(std::ostream &os, enum BarrierType bi) {
  const char *names[] = {"NONE", "D1", "D2", "D3", "DW",
                         "DP",   "M1", "M2", "M3"};
  os << names[bi];
  return os;
}

std::ostream &operator<<(std::ostream &os, const BarrierIndex &bi) {
  os << bi.type << " " << bi.index;
  return os;
}

void Barriers::print() {

  std::cout << "labelBase: " << labelBase << std::endl;
  std::cout << "pointBase: " << pointBase << std::endl;
  std::cout << "surfaceBase : " << surfaceBase << std::endl;
  // std::cout << "shapePrior : " << shapePrior << std::endl;
  std::cout << "jointPrior: " << jointPrior << std::endl;
  std::cout << "poseBase: " << poseBase << std::endl;
  std::cout << "absolutePoseBase " << std::endl;
  std::cout << "zero: " << zero << std::endl;
  std::cout << "one: " << one << std::endl;
  std::cout << "two: " << two << std::endl;
  std::cout << "three: " << three << std::endl;
  std::cout << "meanShapePrior: " << meanShapePrior << std::endl;

  // std::cout << "metric: " << metric << std::endl;
  // std::cout << "weightInit: " << weightInit << std::endl;

  std::cout << "registration_barrier: " << registration << std::endl;
  std::cout << "registration_local_barrier: " << localRegistration << std::endl;
  std::cout << "shape_barrier: " << shape << std::endl;
  std::cout << "joints barrier : " << joints << std::endl;
  std::cout << "deform barrier : " << deform << std::endl;
  std::cout << "pose barrier: " << pose << std::endl;
  // std::cout << "pjoint barrier: " << pjoint << std::endl;
  std::cout << "end3d barrier: " << end3d << std::endl;

  std::cout << "coeff barrier: " << coeff << std::endl;
  std::cout << "weights barrier: " << weights << std::endl;
  std::cout << "end1d barrier: " << end1d << std::endl;
}
