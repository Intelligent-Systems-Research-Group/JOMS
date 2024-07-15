#include "corpora.h"

#include "rotations.h"
#include <Eigen/Geometry>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>

#include "INIReader.h"

int loadCorpora(std::string path, bool isFlipped,
                std::vector<std::string> *result,
                std::vector<int> *personScanCount, std::vector<bool> *flipped,
                std::vector<Vector3> *pos, std::vector<int> *person_ids,
                int person_id) { //, std::vector<int>* seq
  std::string line;
  std::ifstream infile(path);
  int todo = 0;
  std::getline(infile, line);
  int scans = 1;
  std::stoi(line);

  int computePerson = -1;
  int cap = 0;
  int sequence = 0;
  // int maxcap = 541;//473; 524 541
  while (std::getline(infile, line) /* && cap < maxcap*/) {
    std::cout << line << std::endl;
    std::cout << "Line: " << line << std::endl;
    if (line == "")
      break;
    else if (todo == 0) {
      todo = std::stoi(line);
      std::cout << "Scans of the Person chars: " << todo << std::endl;
      // todo = maxcap;
      computePerson++;
      if (computePerson == 1) {
        break;
      }
      if (personScanCount) {
        personScanCount->push_back(todo);
      }
      // if (seq) {
      //}
    } else {
      if (person_ids)
        person_ids->push_back(person_id);
      result->push_back(line);
      if (pos != NULL)
        pos->push_back(Vector3(0, 0, 0));
      cap++;
      if (flipped)
        flipped->push_back(false);
      // sequence++;

      // if(seq)
      todo--;
    }
  }
  if (personScanCount) {
    std::cout << "finished " << personScanCount->size() << std::endl;
    std::cout << "Persons: " << personScanCount->size()
              << " Scans: " << (*personScanCount)[0]
              << "Real Scans: " << result->size() << std::endl;
    // for(int i = 0; i < result->size(); i++) {
    //	std::cout << i << " : " << (*result)[i] << std::endl;
    //}
  }
  return scans;
}

std::vector<std::string> loadFilePaths(std::string path) {
  std::ifstream infile(path);

  if (!infile.is_open()) {
    std::cout << "file not found" << std::endl;
  }

  int i = 0;

  std::vector<std::string> result;
  std::string line;
  while (std::getline(infile, line)) {
    result.push_back(line);
    i++;
    // if(i > 100) break;
  }
  return result;
}

static void split(const std::string &s, char c, std::vector<std::string> &v) {
  using namespace std;
  string::size_type i = 0;
  string::size_type j = s.find(c);
  if (j == string::npos) {
    v.push_back(s);
    return;
  }
  while (j != string::npos) {
    v.push_back(s.substr(i, j - i));
    i = ++j;
    j = s.find(c, j);

    if (j == string::npos)
      v.push_back(s.substr(i, s.length()));
  }
}

Corpora constructFromConfig(std::string path) {
  Corpora corpora;
  corpora.personCount = 0;

  std::vector<std::string> paths;
  split(path, ':', paths);
  assert(paths.size() > 0);
  for (auto p : paths) {
    std::cout << "Reading Paths " << p << std::endl;
    INIReader reader(p);
    auto person_config_folder =
        reader.GetString("dataset", "person_config_path", "NaN");
    auto person_ids = reader.GetStringVector("dataset", "person_ids");
    auto person_config_suffixes =
        reader.GetStringVector("dataset", "person_config_suffixes");
    auto person_has_mask = reader.GetBoolVector("dataset", "person_has_mask");
    auto person_has_labels_3d =
        reader.GetBoolVector("dataset", "person_has_labels_3d");
    auto person_has_labels_2d =
        reader.GetBoolVector("dataset", "person_has_labels_2d");
    auto person_fist = reader.GetBoolVector("dataset", "person_fist");
    auto person_ground_planes =
        reader.GetFloatVector("dataset", "person_ground_planes");
    auto person_ground_axes =
        reader.GetIntVector("dataset", "person_ground_axes");

    assert(corpora.person_config_folder == "" ||
           corpora.person_config_folder == person_config_folder);
    corpora.person_config_folder = person_config_folder;
    corpora.person_ids.insert(corpora.person_ids.end(), person_ids.begin(),
                              person_ids.end());
    corpora.person_config_suffixes.insert(corpora.person_config_suffixes.end(),
                                          person_config_suffixes.begin(),
                                          person_config_suffixes.end());
    corpora.person_has_mask.insert(corpora.person_has_mask.end(),
                                   person_has_mask.begin(),
                                   person_has_mask.end());
    corpora.person_has_labels_3d.insert(corpora.person_has_labels_3d.end(),
                                        person_has_labels_3d.begin(),
                                        person_has_labels_3d.end());
    corpora.person_has_labels_2d.insert(corpora.person_has_labels_2d.end(),
                                        person_has_labels_2d.begin(),
                                        person_has_labels_2d.end());
    corpora.person_fist.insert(corpora.person_fist.end(), person_fist.begin(),
                               person_fist.end());
    corpora.person_ground_planes.insert(corpora.person_ground_planes.end(),
                                        person_ground_planes.begin(),
                                        person_ground_planes.end());
    corpora.person_ground_axes.insert(corpora.person_ground_axes.end(),
                                      person_ground_axes.begin(),
                                      person_ground_axes.end());
  }

  int npersons = corpora.person_ids.size();
  for (int i = 0; i < npersons; i++) {
    std::string person_config = corpora.person_config_folder +
                                corpora.person_ids[i] + "_" +
                                corpora.person_config_suffixes[i] + ".txt";
    std::cout << "CONFIG: " << person_config << "\n";
    bool flipped = corpora.person_fist[i];
    loadCorpora(person_config, flipped, &corpora.scan_names, &corpora.scanCount,
                &corpora.flipped, &corpora.rot, &corpora.scan_to_person_id, i);
    loadCorpora(person_config, flipped, &corpora.reg_names, NULL, NULL, NULL,
                NULL, i);
    corpora.folder.push_back(corpora.person_ids[i]);
  }
  std::cout << "DFaust stats, nscans: " << corpora.scan_names.size()
            << std::endl;
  // exit(0);
  corpora.personCount = npersons;
  corpora.markerFrames = corpora.scan_names.size();

  // corpora.personCount += newPersons;
  corpora.sequences.resize(14);

  const int seq[14] = {25, 77, 27, 52, 54, 27, 23, 61, 28, 35, 36, 26, 27, 28};
  int nextStart = 0;
  for (int i = 0; i < 14; i++) {
    corpora.sequences[i] = nextStart;
    nextStart += seq[i];
  }
  return corpora;
}

void printCorpora(const Corpora &corpora) {
  std::cout << "Person Count " << corpora.personCount << std::endl;
  std::cout << "scan_names"
            << " " << corpora.scan_names.size() << std::endl;
  std::cout << "scanCount " << corpora.scanCount.size() << std::endl;
  for (int i = 0; i < corpora.personCount; i++) {
    std::cout << "Person " << i << "Scans" << corpora.scanCount[i];
    for (int j = 0; j < corpora.scanCount[i]; j++) {
      std::cout << j << " " << corpora.scan_names[j] << " "
                << corpora.reg_names[j] << " " << corpora.scan_names[j] << " "
                << std::endl;
    }
  }
}
