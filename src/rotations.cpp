#include "rotations.h"

/*
void rotation_matrix_to_angle_axis(
        const Matrix3& R,
        Vector3 *angle_axis)
{
        double angle = acos((R.trace() - 1) / 2);

        if (angle < .0000001)
        {
                angle_axis->setZero();
                return;
        }

        double x = R(2, 1) - R(1, 2);
        double y = R(0, 2) - R(2, 0);
        double z = R(1, 0) - R(0, 1);

        double s = 1 / (2 * sin(angle));

        (*angle_axis)(0, 0) = s*x;
        (*angle_axis)(1, 0) = s*y;
        (*angle_axis)(2, 0) = s*z;
}
*/

void angle_axis_to_rotation_matrix(const Vector3 &angle_axis, Matrix3 *R) {
  Scalar norm = angle_axis.norm();

  if (norm < .0000001) {
    R->setIdentity();
    return;
  }

  Scalar x = angle_axis[0] / norm;
  Scalar y = angle_axis[1] / norm;
  Scalar z = angle_axis[2] / norm;

  Scalar s = sin(norm);
  Scalar c = cos(norm);

  *R << x * x + (1 - x * x) * c, x * y * (1 - c) - z * s,
      x * z * (1 - c) + y * s, x * y * (1 - c) + z * s, y * y + (1 - y * y) * c,
      y * z * (1 - c) - x * s, x * z * (1 - c) - y * s, z * y * (1 - c) + x * s,
      z * z + (1 - z * z) * c;
}

void euler_to_rotation_matrix(const Vector3 &angle_axis, Matrix3 *R) {

  Scalar alpha = angle_axis[0];
  Scalar beta = angle_axis[1];
  Scalar gamma = angle_axis[2];
  Scalar c1 = cos(alpha);
  Scalar c2 = cos(beta);
  Scalar c3 = cos(gamma);
  Scalar s1 = sin(alpha);
  Scalar s2 = sin(beta);
  Scalar s3 = sin(gamma);

  *R << c1 * c2, c1 * s2 * s3 - c3 * s1, s1 * s3 + c1 * c3 * s2, c2 * s1,
      c1 * c3 + s1 * s2 * s3, c3 * s1 * s2 - c1 * s3, -s2, c2 * s3, c2 * c3;
}

/*
void angle_axis_to_rotation_matrix(
        const Vector3& angle_axis,
        Matrix3 *R)
{

    Scalar alpha = angle_axis[0];
    Scalar beta = angle_axis[1];
    Scalar gamma= angle_axis[2];
    Scalar CosAlpha = cos(alpha);
    Scalar CosBeta = cos(beta);
    Scalar CosGamma = cos(gamma);
    Scalar SinAlpha = sin(alpha);
    Scalar SinBeta = sin(beta);
    Scalar SinGamma = sin(gamma);

        *R << CosGamma*CosBeta,
            -SinGamma*CosAlpha + CosGamma*SinBeta*SinAlpha,
            SinGamma*SinAlpha + CosGamma*SinBeta*CosAlpha,
            SinGamma*CosBeta,
            CosGamma*CosAlpha + SinGamma*SinBeta*SinAlpha,
            -CosGamma*SinAlpha + SinGamma*SinBeta*CosAlpha,
            -SinBeta,
            CosBeta*SinAlpha,
            CosBeta*CosAlpha;
}
*/

void rotation_matrix_to_angle_axis(const Matrix3 &R, MatrixX *angle_axis) {
  Scalar e = 0.000001;
  Scalar t = R.trace();

  Scalar rx = R(2, 1) - R(1, 2);
  Scalar ry = R(0, 2) - R(2, 0);
  Scalar rz = R(1, 0) - R(0, 1);

  // std::cout << "t " << t << std::endl;

  if (t >= (3 - e)) {
    (*angle_axis)(0, 0) = (0.5 - (t - 3) / 12) * rx;
    (*angle_axis)(1, 0) = (0.5 - (t - 3) / 12) * ry;
    (*angle_axis)(2, 0) = (0.5 - (t - 3) / 12) * rz;
    return;
  }
  if ((3 - e) > t && t > (-1 + e)) {
    Scalar theta = acos((t - 1.0) / 2.0);
    Scalar f = theta / (2 * sin(theta));
    (*angle_axis)(0, 0) = f * rx;
    (*angle_axis)(1, 0) = f * ry;
    (*angle_axis)(2, 0) = f * rz;
    return;
  }
  int a = 0;
  if (R(0, 0) >= R(1, 1) && R(0, 0) >= R(2, 2)) {
    a = 0;
  } else if (R(1, 1) >= R(0, 0) && R(1, 1) >= R(2, 2)) {
    a = 1;
  } else {
    a = 2;
  }
  int b = (a + 1) % 3;
  int c = (a + 2) % 3;

  Scalar s = sqrt(R(a, a) - R(b, b) - R(c, c) + 1);
  Scalar va = s / 2;
  Scalar vb = 1 / (2 * s) * (R(b, a) - R(a, b));
  Scalar vc = 1 / (2 * s) * (R(c, a) - R(a, c));
  Scalar vnorm = sqrt(va * va + vb * vb + vc * vc);

  (*angle_axis)(0, 0) = (M_PI / vnorm) * va;
  (*angle_axis)(1, 0) = (M_PI / vnorm) * vb;
  (*angle_axis)(2, 0) = (M_PI / vnorm) * vc;
}
