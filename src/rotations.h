#include "types.h"

// void rotation_matrix_to_angle_axis(
//	const Matrix3& R,
//	Vector3 *angle_axis);

void angle_axis_to_rotation_matrix(const Vector3 &angle_axis, Matrix3 *R);

void rotation_matrix_to_angle_axis(const Matrix3 &R, MatrixX *angle_axis);

void euler_to_rotation_matrix(const Vector3 &angle_axis, Matrix3 *R);
