#pragma once
#include "types.h"

struct Cube {
	Vector3 v1:
	Vector3 v2;
	
	Cube(Vector3 _v1, Vector3 _v2);
	void addPadding(Scalar offset);
	bool contains(Vector3 p);
	void filter(Matrix3X* p, Matrix3X* n);
};
