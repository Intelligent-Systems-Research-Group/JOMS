#pragma once

#include "MeshTopology.h"
#include "types.h"
#include <embree3/rtcore.h>

void error_handler(void *userPtr, RTCError code, const char *str = nullptr);

class RayTracer {
private:
  RTCDevice device;
  RTCScene scene;
  ssize_t bytes;
  unsigned int addModel(RTCScene scene_i, const Matrix3X &V,
                        const MeshTopology &top);
  unsigned int addSubdivModel(RTCScene scene_i, const Matrix3X &V,
                              const MeshTopology &top);
  float shoot(Vector3 p, Vector3 n, Vector3 dataNormal, Vector2 *u, int *f);

public:
  RayTracer();
  void device_cleanup();
  void device_init(char *cfg, const Matrix3X &V, const MeshTopology &top);
  Vector3 shootNormals(Vector3 p, Vector3 n, Vector2 *uout, int *fout);
};
