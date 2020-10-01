#pragma once
#include <string>
#include "model.h"


//bool writeDebug(std::string path, Global* global);
//bool writeFbx(std::string path, int* mode, Global* global);
//bool readFbx(std::string path, int* mode, Global* global);
bool writePersonFbx(std::string path, int* mode, Model* model, int personId, int poseId, int endId);
//bool convertFbx(std::string path, std::string outpath);
bool reparameterize(const Model& source, Model* target);
