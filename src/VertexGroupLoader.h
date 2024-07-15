#pragma once
#include "types.h"
#include <map>
#include <string>

void loadVertexGroups(std::string path,
                      std::map<std::string, std::vector<size_t>> *,
                      std::map<std::pair<std::string, std::string>,
                               std::vector<size_t>> *border_groups,
                      Weights *weights);
std::map<std::pair<std::string, std::string>, Vector3> calculateJoints(
    const std::map<std::pair<std::string, std::string>, std::vector<size_t>>
        &vertex_groups,
    Matrix3X *vertices);
std::map<std::pair<std::string, std::string>, VectorX>
calculateInitialRegressor(const std::map<std::pair<std::string, std::string>,
                                         std::vector<size_t>> &vertex_groups,
                          Matrix3X *vertices);

VectorX loadVector(std::string path);
