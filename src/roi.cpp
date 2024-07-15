#include "roi.h"

Cube::Cube(Vector3 _v1, Vector3 _v2) : v1(_v1), v2(_v2) {}
void Cube::addPadding(Scalar offset) {
  v1 -= offset;
  v2 -= offset;
}
bool contains(Vector3 p) {
  return v1[0] < p[0] && v1[1] < p[1] && v1[2] < p[2] && v2[0] > p[0] &&
         v2[1] > p[1] && v2[2] > p[2];
}
void filter(Matrix3X *p, Matrix3X *n) {
  Matrix3X p2 = *p;
  Matrix3X n2 = *n;

  int j = 0;
  for (int i = 0; i < p->size(); i++) {
    Vector3 q = p->col(i);
    Vector3 qn = n->col(i);
    if (contains(q)) {
      p2.col(j) = q;
      n2.col(j) = qn;
      j++;
    }
  }
  p2.conservativeResize(3, j);
  n2.conservativeResize(3, j);

  *p = p2;
  *n = n2;
}
