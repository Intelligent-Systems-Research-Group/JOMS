#pragma once
#include "types.h"
#include <string>
#include <vector>

//@TODO load with config
struct Corpora {

  std::string person_config_folder;
  std::vector<std::string> person_ids;
  std::vector<std::string> person_config_suffixes;
  std::vector<bool> person_has_mask;
  std::vector<bool> person_has_labels_3d;
  std::vector<bool> person_has_labels_2d;
  std::vector<bool> person_fist;
  std::vector<Scalar> person_ground_planes;
  std::vector<int> person_ground_axes;

  int personCount;
  int markerFrames;

  // scans per person
  std::vector<int> scanCount;

  // file names for markers
  std::vector<std::string> reg_names;

  // file names for scans
  std::vector<std::string> scan_names;

  std::vector<int> scan_to_person_id;
  std::vector<std::string> folder;
  std::vector<bool> flipped;
  std::vector<Vector3> rot;
  std::vector<int> sequences;
};

Corpora constructFromConfig(std::string path);
void printCorpora(const Corpora &corpora);
