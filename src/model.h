#pragma once

#include "MeshTopology.h"
#include "skeleton.h"
#include <iostream>
#include <map>
#include <string>

//@TODO: save and load
struct Scan {
  Matrix3X trans;
  Matrix3X absoluteTrans;
  Matrix3X pose;
  Matrix3X posebase;
  Matrix3X absolutePose;
  Matrix3X absolutePosebase;
  Matrix3X Trest;
  Matrix3X T;
  Matrix3X Tlocal;
  Matrix3X labels;
  Matrix3X labels_initial;
  Matrix3X labels_partial;

  Matrix3X T_estimate;
  Matrix3X Tlocal_estimate;
  Matrix3X points;
  Matrix3X patchPoints; // points computed using lifting
  Matrix3X tangentStep;
  VectorX robust;

  Matrix3X joints;
  std::vector<Matrix3X> label2dRays; // jointabels projected in world space
};

struct Person {
  // std::vector<float> coeff;
  VectorX coeff;
  // Matrix3X restShape;
  Matrix3X joints;
};

struct Model {
  IMatrix4X closestJoints;
  IMatrix4X closestDeformJoints;
  hyp::Skeleton *skeleton;
  MeshTopology top;

  /*
      Model Parameters
  */
  Matrix3X mean;
  std::vector<Matrix3X> pcs;
  Matrix3X joints;
  std::vector<Matrix3X> jointpcs;
  MatrixX weights;
  std::vector<Matrix3X> deform;
  Matrix3X meanPose;
  /*
      End Model Parameters
  */

  std::vector<Person> persons;
  std::vector<Scan> scans;

  void write(std::string path);
  void read(std::string path);
  void init(const Matrix3X &_mean, const Matrix3X &_joints,
            const MatrixX &_weights, int ncomp, int ndefcomp, int nperson,
            int nscans);
};

struct InputParams {
  int npoints;
  int nlabels;
  int njointlabels;
  int nviews;
  int nverts;
  int njoints;
  int nvarjoints;
  int ncomp;
  int nscans;
  int ndeformjoints;
  int ndeformshapes;
  int nbasis;
  float continuous_step_size;
  int iteration_end;
  int npersons;
  std::string optscript_path;
  std::string out_path;
  std::string model_weights_path;
  std::string base_dataset_folder;
  bool freeze_weights;
  bool use_symmetry;
  int dump_progress;
  int raycasting;
  int jointly;
  int alternating;
  bool useModel2Data;
  bool useTemporal;
  int iteration_data;
  int max_vertex_ring;
  int center_vertex_id;
  int vertex_masking_end;
  std::vector<std::string> exclude_icp_vertex_groups;
  std::vector<std::string> camPaths;

  std::string template_path;
  std::string template_data_marker_path;
  std::string template_marker_path;
  std::string template_vgroups_path;
  std::string template_mesh_path;
  std::string template_laplace_path;
  std::string template_skeleton_path;
  std::string template_surface_map_path;
};

enum BarrierType { NONE, D1, D2, D3, DW, M1, M2, M3 };

std::ostream &operator<<(std::ostream &os, enum BarrierType bi);

class BarrierIndex {
public:
  int index;
  BarrierType type;

public:
  BarrierIndex(BarrierType _type = NONE, int _index = 0)
      : type(_type), index(_index) {}
  BarrierIndex operator+(int offset) const {
    return BarrierIndex(this->type, this->index + offset);
  }

  BarrierIndex &operator+=(int offset) {
    this->index += offset;
    return *this;
  }

  BarrierIndex operator++(int) {
    BarrierIndex temp(*this);
    this->index++;
    return temp;
  }
  int operator()() const { return this->index; }
  bool hasType(BarrierType _type) const { return _type == type; }

  bool sharesType(const BarrierIndex &bi) const { return type == bi.type; }

  BarrierType getType() const { return type; }
  // operator size_t() const { return this->index; }

  bool isFree() const {
    const bool var[] = {false, false, false, false, false,
                        false, true,  true,  true};
    return var[type];
  }
};

std::ostream &operator<<(std::ostream &os, const BarrierIndex &bi);

struct Barriers {
  BarrierIndex surface;
  BarrierIndex surfaceScale;
  BarrierIndex surfaceLabel; // DO NOT CHANGE ORDER IN MEMORY
  BarrierIndex vertexLabel;  // DO NOT CHANGE ORDER IN MEMORY
  BarrierIndex labelBase;
  BarrierIndex pointBase;
  BarrierIndex surfaceBase;

  // 2D pose estimation labels for joint keypointso
  BarrierIndex label2d;

  // 2D pose estimation labels for surface keypoints
  BarrierIndex labelSurface2d;
  // int shapePrior;
  BarrierIndex jointPrior;
  BarrierIndex poseBase;
  BarrierIndex absolutePoseBase;
  BarrierIndex poseBaseMean;
  //(0,0,0)
  BarrierIndex zero;
  //(1,1,1)
  BarrierIndex one;
  //(2,2,2)
  BarrierIndex two;
  //(3,3,3)
  BarrierIndex three;
  // nviews*2 (R0,t0,R1,t1,...,Rn,tn) openpose camera extrinsics
  BarrierIndex cameraViews;
  BarrierIndex meanShapePrior;

  BarrierIndex localStencil;

  BarrierIndex subdivLabelWeights;
  BarrierIndex labelMask;
  // always set to 2
  BarrierIndex singleTwo;
  // always set to 1
  BarrierIndex singleOne;
  // always set to 0
  BarrierIndex singleZero;
  BarrierIndex groundPlanes;
  BarrierIndex personCount;

  // 2D pose estimation weights for joint keypoints
  BarrierIndex labelWeights2d;
  // 2D pose estimation weights for surface keypoints
  BarrierIndex labelSurfaceWeights2d;

  // Laplace-Beltrami coefficients for shape regularization
  BarrierIndex metric;
  BarrierIndex weightInit;

  BarrierIndex endD1;
  BarrierIndex endD2;
  BarrierIndex endD3;
  BarrierIndex endDW;

  BarrierIndex registration;
  BarrierIndex localRegistration;
  BarrierIndex debugRegistration;
  BarrierIndex shape;
  BarrierIndex joints;
  BarrierIndex deform;
  BarrierIndex pose;
  BarrierIndex absolutePose;
  BarrierIndex poseMean;
  BarrierIndex pjoint;
  BarrierIndex sjoint;
  BarrierIndex pshape;
  BarrierIndex end3d;

  BarrierIndex surfaceUnknown;
  BarrierIndex end2d;

  BarrierIndex robust;
  BarrierIndex coeff;
  BarrierIndex weights;
  BarrierIndex end1d;

  void print();
};
