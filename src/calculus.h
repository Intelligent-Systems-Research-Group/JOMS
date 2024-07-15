#pragma once

#include "types.h"

Matrix32 principleCurvatureAmbient(Vector3 p, Vector3 dSdu, Vector3 dSdv,
                                   Vector3 dSduu, Vector3 dSduv, Vector3 dSdvv);

Matrix22 principleCurvatureTangent(Vector3 p, Vector3 dSdu, Vector3 dSdv,
                                   Vector3 dSduu, Vector3 dSduv, Vector3 dSdvv);
