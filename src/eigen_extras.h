#pragma once

#include <Eigen/Eigen>
/*
typedef double Scalar;
// Like the eigen typedefs, but using the Scalar template parameter
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::DontAlign>
MatrixXX; typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic, Eigen::DontAlign>
Matrix3X; typedef Eigen::Matrix<Scalar, 2, Eigen::Dynamic, Eigen::DontAlign>
Matrix2X; typedef Eigen::Matrix<Scalar, -1, 3, Eigen::Dynamic, Eigen::DontAlign>
MatrixX3; typedef Eigen::Matrix<Scalar, -1, 2, Eigen::Dynamic, Eigen::DontAlign>
MatrixX2; typedef Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign> Matrix22;
typedef Eigen::Matrix<Scalar, 3, 2, Eigen::DontAlign> Matrix32;
typedef Eigen::Matrix<Scalar, 2, 3, Eigen::DontAlign> Matrix23;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Eigen::DontAlign> VectorX;
typedef Eigen::Matrix<Scalar, 2, 1, Eigen::DontAlign> Vector2;
typedef Eigen::Matrix<Scalar, 3, 1, Eigen::DontAlign> Vector3;


template <typename T, int _Options, typename _Index>
void write(Eigen::SparseMatrix<T, _Options, _Index> const& J, char const*
filename)
{
  std::ofstream f{ filename };
  if (!f.good()) {
    std::cerr << "Failed to open [" << filename << "] for writing\n";
    return;
  }
  std::cout << "Writing " << J.rows() << "x" << J.cols() << " sparse to [" <<
filename << "]\n"; for (int k = 0; k < J.outerSize(); ++k) for
(Eigen::SparseMatrix<T, _Options, _Index>::InnerIterator it(J, k); it; ++it) f
<< it.row() << "\t" << it.col() << "\t" << it.value() << std::endl;
}

template <typename Derived>
void write(Eigen::MatrixBase<Derived> const& J, char const* filename)
{
  std::ofstream f{ filename };
  if (!f.good()) {
    std::cerr << "Failed to open [" << filename << "] for writing\n";
    return;
  }
  std::cout << "Writing " << J.rows() << "x" << J.cols() << " dense to [" <<
filename << "]\n"; f << J;
}
*/

/*
template <typename M1, typename M2>
auto hcat(M1 const& A, M2 const& B)
{
  M1 out(A.rows(), A.cols() + B.cols());
  out << A, B;
  return out;
}
*/
