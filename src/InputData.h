#pragma once
#include "types.h"
#include "model.h"

#include <vector>

class InputData {

    //std::vector<std::array<Scalar,PATCH_SIZE>> fpd;
    std::vector <Scalar3> f3d;
    std::vector <Scalar2> f2d;
    std::vector <Scalar> f1d;

    std::vector <std::array<Scalar, LOCAL_SIZE>> dwd;
    std::vector <Scalar> d1d;
    std::vector <Scalar2> d2d;
    std::vector <Scalar3> d3d;


    //std::vector<std::array<Scalar,PATCH_SIZE>> backup_fpd;
    std::vector <Scalar> backup_f1d;
    std::vector <Scalar2> backup_f2d;
    std::vector <Scalar3> backup_f3d;

    std::vector <std::array<Scalar, LOCAL_SIZE>> backup_dwd;
    std::vector <Scalar> backup_d1d;
    std::vector <Scalar2> backup_d2d;
    std::vector <Scalar3> backup_d3d;


    //Scalar* cfpd;
    Scalar3 *cf3d;
    Scalar2 *cf2d;
    Scalar *cf1d;

    Scalar *cdwd;
    Scalar3 *cd3d;
    Scalar2 *cd2d;
    Scalar *cd1d;
public:
    void verifySize(BarrierIndex f3, BarrierIndex f2, BarrierIndex f1);

    void verifyDataSize(BarrierIndex d1, BarrierIndex d2,
                        BarrierIndex d3, BarrierIndex dw);

    bool verifyIndex(BarrierIndex i);

    void checkpoint();

    void restore();

    void addUnknown(Scalar3 _f3d);

    void addUnknown(Scalar2 _f2d);

    void addUnknown(Scalar _f1d);

    void addPoint(std::array <Scalar, LOCAL_SIZE> _dwd);


    void addPoint(Scalar3 _d3d);

    void addPoint(Scalar2 _d2d);

    void addPoint(Scalar _d1d);

    void addDynamic(Scalar x, bool isUnknown);

    void addDynamic(Scalar2 x, bool isUnknown);

    void addDynamic(Scalar3 x, bool isUnknown);


    void setPoint(BarrierIndex index, Scalar3 _d3d);

    void setPoint(BarrierIndex index, Scalar2 _d2d);

    void setPoint(BarrierIndex index, Scalar _d1d);

    std::array <Scalar, LOCAL_SIZE> readPoint(BarrierIndex index);


    Scalar3 readPoint3d(BarrierIndex index);

    Scalar2 readPoint2d(BarrierIndex index);

    Scalar readPoint1d(BarrierIndex index);
	
	Vector3 readPoint3dE(BarrierIndex index);

    Vector2 readPoint2dE(BarrierIndex index);

    void resetUnknown3d(BarrierIndex index);

    void setUnknown3d(BarrierIndex index, Scalar3 _f3d);

    void resetUnknown2d(BarrierIndex index);

    void allocate();

    void release();

    void copy();

    void retrieve();

    void getDims(unsigned int *dims);

    void assertDims(unsigned int *dims);

    unsigned int dim() const;

    void assign(void **result);

    void read(Matrix3X &A, BarrierIndex start, int size, bool backup = false);

    void read(Matrix2X &A, BarrierIndex start, int size);

	Vector3 readVector3(BarrierIndex start, bool backup = false);
    Vector2 readVector2(BarrierIndex start);


    void readPoints(Matrix3X &A, BarrierIndex start, int size);
	void readPoints(Matrix2X &A, BarrierIndex start, int size);


    void clamp2d(BarrierIndex start, int size);

    void write(const VectorX &A, BarrierIndex start, int size);

    void readPoints(VectorX &A, BarrierIndex start, int size);

    void read(VectorX &A, BarrierIndex start, int size);

    Scalar read(BarrierIndex i);

    Scalar readDynamic1d(BarrierIndex i);


    void readDynamic(Matrix3X &A, BarrierIndex start, int size);

    void readDynamic(VectorX &A, BarrierIndex start, int size);

    void write(Scalar x, BarrierIndex i);
};
