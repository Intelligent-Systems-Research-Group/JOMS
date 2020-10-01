#pragma once

//#define NDEBUG
#include <float.h>


#include <stdexcept>
//#undef eigen_assert
//#define eigen_assert(x) \
//  if (!(x)) { throw (std::runtime_error("Put your message here")); }

//#include <gperftools/heap_profiler.h>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>
#include <map>
#include <cuda_runtime.h>
//#include <adolc/adolc.h>
//#include <adolc/adolc_sparse.h>


//#define VERBOSE_DUMP
#define FULL_DUMP
#define M_PI 3.14159265358979323846
#define SMALL_CONFIG

#define USE_MARKERS

#define SPLIT_JOINTS
#define USE_RINGREG
//#define PRECOMPUTE_SIZE
#define USE_NORMAL
//#define USE_DOUBLE

#define MODE (0)
#define MODEL_FREE (MODE == 0)
#define PREREGISTERED (MODE == 1)
#define FULL_TRAINING (MODE == 2)

#define DATASET (1)
#define DATASET_SINGLE (DATASET == 0)
#define DATASET_FULL (DATASET == 1)
#define DATASET_PARTIAL (DATASET == 2)

#define PARTIAL_FREE
//#define POSE_DEFORM
#define REG_SMOOTH
//6890
//#define N_DATA_POINTS 10000
#define NTHREADS (4)
#define SAMPLING (0.005)
//0.005
//.1
#define RAY_MAX (.081)


//#define NTHREADS (omp_get_max_threads())

#define MyAlignment Eigen::DontAlign
//typedef float Scalar;
#define PATCH_SIZE (16*3)
//#define LOCAL_SIZE (28)
#define LOCAL_SIZE (28)
//old pose size is 16*3 new was 33*3
#define POSE_SIZE (NJOINTS*3) 
//#define MAX_VERTEX_RING (16)


#ifdef SMALL_CONFIG
#define NO_LABEL_TAIL (0)
#else
#define NO_LABEL_TAIL (2)
#endif

#define WMIN (0.001)
#define WMAX (1-WMIN)
//(i >= 6 and i < 15) or (i >=17 and i < 26)
//#define RIGID(i) (((i) >= 5 && (i) < 14) || ((i) >=16 && (i) < 25))
//#define RIGID(i) (false)

#ifndef USE_DOUBLE
typedef float  Scalar;
typedef float2 Scalar2;
typedef float3 Scalar3;
#else
typedef double  Scalar;
typedef double2 Scalar2;
typedef double3 Scalar3;
#endif

#define logging (std::cout)

typedef Eigen::Matrix<Scalar, 2, 1, MyAlignment> Vector2;
typedef Eigen::Matrix<Scalar, 3, 1, MyAlignment> Vector3;
typedef Eigen::Matrix<Scalar, 4, 1, MyAlignment> Vector4;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, MyAlignment> VectorX;
typedef Eigen::Matrix<Scalar, 2, 2, MyAlignment> Matrix2;
typedef Eigen::Matrix<Scalar, 3, 3, MyAlignment> Matrix3;
typedef Eigen::Matrix<Scalar, 4, 4, MyAlignment> Matrix4;

typedef Eigen::Matrix<Scalar, 2, 2, MyAlignment> Matrix22;
typedef Eigen::Matrix<Scalar, 3, 2, MyAlignment> Matrix32;
typedef Eigen::Matrix<Scalar, 3, 4, MyAlignment> Matrix34;

typedef Eigen::Matrix<Scalar, 2, Eigen::Dynamic, MyAlignment> Matrix2X;
typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic, MyAlignment> Matrix3X;
typedef Eigen::Matrix<Scalar, 4, Eigen::Dynamic, MyAlignment> Matrix4X;
typedef Eigen::Matrix<Scalar, PATCH_SIZE, Eigen::Dynamic, MyAlignment> MatrixPX;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, MyAlignment> MatrixX;


//typedef pdouble Modifyable;
//typedef Eigen::Matrix<Modifyable, 3, Eigen::Dynamic, MyAlignment> PMatrix3X;

typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, MyAlignment> IMatrixX;
typedef Eigen::Matrix<size_t, 2, Eigen::Dynamic, MyAlignment> IMatrix2X;
typedef Eigen::Matrix<size_t, 4, Eigen::Dynamic, MyAlignment> IMatrix4X;

typedef Eigen::Matrix<int, 4, 1, MyAlignment> SVector4;
typedef Eigen::Matrix<int, -1, 1, MyAlignment> IVectorX;

inline Scalar3 toScalar3(Vector3 x) {
	Scalar3 vert;
    vert.x = x[0];
    vert.y = x[1];
    vert.z = x[2];
	return vert;
}


template<typename Derived>
inline bool is_finite(const Eigen::MatrixBase<Derived>& x)
{
	return ((x - x).array() == (x - x).array()).all();
}

template<typename Derived>
inline bool is_nan(const Eigen::MatrixBase<Derived>& x)
{
	return !((x.array() == x.array())).all();
}

template<typename Derived>
inline bool is_sane(const Eigen::MatrixBase<Derived>& x)
{
        return is_finite(x) && !is_nan(x);
}

template<typename U, typename V>
void setConstant(U& u, V& v)
{
	assert(u.rows() == v.rows() && u.cols() == v.cols());
	for (size_t i = 0; i < u.rows(); i++) {
		for (size_t j = 0; j < v.cols(); j++) {
			u(i, j) = v(i, j);
		}
	}
}

template<typename U, typename V>
void assign(U& u, const V& v)
{
	assert(u.rows() == v.rows() && u.cols() == v.cols());
	for (size_t i = 0; i < u.rows(); i++) {
		for (size_t j = 0; j < v.cols(); j++) {
			u(i, j) <<= v(i, j);
		}
	}
}

template<typename U, typename V>
void extract(U& u, V& v){
	assert(u.rows() == v.rows() && u.cols() == v.cols());
	for (size_t i = 0; i < u.rows(); i++) {
		for (size_t j = 0; j < v.cols(); j++) {
			v(i, j) >>= u(i, j);
		}
	}
}

template<typename U, typename V>
void read(U& u, V& v) {
	assert(u.rows() == v.rows() && u.cols() == v.cols());
	for (size_t i = 0; i < u.rows(); i++) {
		for (size_t j = 0; j < v.cols(); j++) {
			u(i, j) = v(i, j).getValue();
		}
	}
}

typedef std::pair<std::string, std::string> JointName;
typedef std::pair<std::string, size_t> JointVertexPair;
typedef std::map<JointVertexPair, Scalar> Weights;
typedef std::pair<size_t, size_t> IndexPair;
typedef std::map<IndexPair, Scalar> SparseMatrix;

