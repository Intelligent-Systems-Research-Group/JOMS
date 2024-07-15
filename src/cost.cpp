#include "cost.h"
#include <iostream>

Cost::Cost(InputData *_inputData, std::vector<TermBuilder> *_terms) {
  inputData = _inputData;
  terms = _terms;
}

Scalar Cost::compute(Scalar _w_fitSqrt, Scalar _w_surface) {
  w_fitSqrt = _w_fitSqrt;
  w_surface = _w_surface;
  inputData->retrieve();

  double cost = 0;
  for (int i = 0; i < 1; i++) {
    TermBuilder &builder = (*terms)[i];
    double termError = 0;
    for (int j = 0; j < builder.count(); j++) {
      std::vector<BarrierIndex> args;
      for (int k = 0; k < builder.countIdxArrays(); k++) {
        args.push_back(builder.getTyped(k, j));
      }
      termError += catmulSimple(args);
    }
    std::cout << "Term Index " << i << " , Error: " << termError << std::endl;
    cost += termError;
  }
  return cost;
}
Scalar Cost::catmulSimple(const std::vector<BarrierIndex> &args) {
  Vector3 d;
  Vector3 dn;
  Vector2 du;
  Vector2 u;
  Vector2 su;
  Scalar w;
  Vector3 Ts[16];
  int offset = 0;

  Scalar LABEL_ITER_MUL = 85.0 / 30.0;
  Scalar THRESH = 30;
  Scalar ENABLE_DATA_TERM = 0; // 20;
  Scalar continuous_step_size = 1.0;

  d = inputData->readPoint3dE(args[offset++]);
  dn = inputData->readPoint3dE(args[offset++]);
  du = inputData->readVector2(args[offset++]);
  u = inputData->readPoint2dE(args[offset++]);
  su = inputData->readPoint2dE(args[offset++]);
  w = inputData->read(args[offset++]);
  for (int i = 0; i < 16; i++) {
    Ts[i] = inputData->readVector3(args[offset++]);
  }
  assert(offset == args.size());
  // fill values

  Vector2 step = w_surface * continuous_step_size * du.cwiseProduct(su);
  // if (w_surface > 0.5)
  //	std::cout << w_surface << " " << step.transpose() << std::endl;
  u = u + step;
  Scalar us[4];
  Scalar vs[4];
  Scalar t, t2, t3;

  t = u[0];
  t2 = t * t;
  t3 = t2 * t;
  us[0] = 1.0 / 6.0 * (1.0 - 3.0 * (t - t2) - t3);
  us[1] = 1.0 / 6.0 * (4.0 - 6.0 * t2 + 3.0 * t3);
  us[2] = 1.0 / 6.0 * (1.0 + 3.0 * (t + t2 - t3));
  us[3] = 1.0 / 6.0 * (t3);

  t = u[1];
  t2 = t * t;
  t3 = t2 * t;
  vs[0] = 1.0 / 6.0 * (1.0 - 3.0 * (t - t2) - t3);
  vs[1] = 1.0 / 6.0 * (4.0 - 6.0 * t2 + 3.0 * t3);
  vs[2] = 1.0 / 6.0 * (1.0 + 3.0 * (t + t2 - t3));
  vs[3] = 1.0 / 6.0 * (t3);

  Vector3 qpbase;
  qpbase.setConstant(0);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      qpbase += us[j] * vs[i] * Ts[4 * i + j];
    }
  }
  Vector3 dp = d - qpbase;
  Vector3 r = dp;

  Scalar frac = 1;
  Scalar cw = w_fitSqrt <= ENABLE_DATA_TERM ? 0 : sqrtf(7.0 / 8.0) * frac;
  cw = .1 * (w_fitSqrt <= LABEL_ITER_MUL * THRESH ? cw : frac);

  Scalar tau = w_fitSqrt <= 70 ? 100 : 50;
  tau = w_fitSqrt <= 80 ? tau : 25;
  tau = w_fitSqrt <= 90 ? tau : 10;
  tau = w_fitSqrt <= 100 ? tau : 2;
  Scalar wreg = cw * (tau / sqrtf(2)) * (w * w - 1);

  Vector3 re = sqrtf(20 * cw) * r * w;
  // return .5*(re.dot(re) + wreg*wreg);
  return r.dot(r);
}
