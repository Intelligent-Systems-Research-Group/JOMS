#include "InputData.h"
#include "types.h"
#include <cuda_runtime.h>
//extern "C" {
//#include "Opt.h"
//}


inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

void InputData::verifySize(BarrierIndex f3, BarrierIndex f2, BarrierIndex f1) {
    dim();
    std::cout << "F3: " << f3() << std::endl;
    std::cout << "F2: " << f2() << std::endl;
    std::cout << "F1: " << f1() << std::endl;

    assert(f3d.size() == f3());
    assert(f2d.size() == f2());
    assert(f1d.size() == f1());
}
void InputData::verifyDataSize(BarrierIndex d1, BarrierIndex d2,
                               BarrierIndex d3, BarrierIndex dw) {

    assert(d1d.size() == d1());
    assert(d2d.size() == d2());
    assert(d3d.size() == d3());
    assert(dwd.size() == dw());
}


bool InputData::verifyIndex(BarrierIndex i) {
    int idx = i();
    assert(idx >= 0);
    switch(i.type) {
        case D1:
            return idx < d1d.size();
            break;
        case D2:
            return idx < d2d.size();
            break;
        case D3:
            return idx < d3d.size();
            break;
        case DW:
            return idx < dwd.size();
            break;
        case M1:
            return idx < f1d.size();
            break;
        case M2:
            return idx < f2d.size();
            break;
        case M3:
            return idx < f3d.size();
            break;
        default:
            assert(false);
    }
}

void InputData::checkpoint() {

    backup_dwd = dwd;
    backup_d1d = d1d;
    backup_d2d = d2d;
    backup_d3d = d3d;
    backup_f1d = f1d;
    backup_f2d = f2d;
    backup_f3d = f3d;
}
void InputData::restore() {
    dwd = backup_dwd;
    d1d = backup_d1d;
    d2d = backup_d2d;
    d3d = backup_d3d;
    f1d = backup_f1d;
    f2d = backup_f2d;
    f3d = backup_f3d;
    copy();
}

void InputData::addUnknown(Scalar3 _f3d) {
    f3d.push_back(_f3d);
}

void InputData::addUnknown(Scalar2 _f2d) {
    f2d.push_back(_f2d);
}

void InputData::addUnknown(Scalar _f1d) {
    f1d.push_back(_f1d);
}

void InputData::addPoint(std::array<Scalar,LOCAL_SIZE> _dwd) {
    dwd.push_back(_dwd);
}

void InputData::addPoint(Scalar3 _d3d) {
    d3d.push_back(_d3d);
}

void InputData::addPoint(Scalar2 _d2d) {
    d2d.push_back(_d2d);
}

void InputData::addPoint(Scalar _d1d) {
    d1d.push_back(_d1d);
}

void InputData::addDynamic(Scalar x, bool isUnknown) {
    isUnknown ? addUnknown(x) : addPoint(x);
}
void InputData::addDynamic(Scalar2 x, bool isUnknown) {
    isUnknown ? addUnknown(x) : addPoint(x);
}
void InputData::addDynamic(Scalar3 x, bool isUnknown) {
    isUnknown ? addUnknown(x) : addPoint(x);
}


void InputData::setPoint(BarrierIndex index, Scalar3 _d3d) {
    assert(index.hasType(D3));
    assert(index() < d3d.size());
    d3d[index()] = _d3d;
}

void InputData::setPoint(BarrierIndex index, Scalar2 _d2d) {
    assert(index.hasType(D2));
    assert(index() < d2d.size());
    d2d[index()] = _d2d;
}

void InputData::setPoint(BarrierIndex index, Scalar _d1d) {
    assert(index.hasType(D1));
    assert(index() < d1d.size());
    d1d[index()] = _d1d;
}

std::array<Scalar,LOCAL_SIZE> InputData::readPoint(BarrierIndex index) {
    assert(index.hasType(DW));
    return dwd[index()];
}


Scalar3 InputData::readPoint3d(BarrierIndex index) {
    assert(index.hasType(D3));
    assert(index() < d3d.size());
    return d3d[index()];
}

Scalar2 InputData::readPoint2d(BarrierIndex index) {
    assert(index.hasType(D2));
    assert(index() < d2d.size());
    return d2d[index()];
}

Scalar InputData::readPoint1d(BarrierIndex index) {
    assert(index.hasType(D1));
    return d1d[index()];
}

Vector3 InputData::readPoint3dE(BarrierIndex index) {
    assert(index.hasType(D3));
    assert(index() < d3d.size());
    Scalar3 res = d3d[index()];
	Vector3 x;
	x[0] = res.x;
	x[1] = res.y;
	x[2] = res.z;
	return x;
}

Vector2 InputData::readPoint2dE(BarrierIndex index) {
    assert(index.hasType(D2));
    assert(index() < d2d.size());
    Scalar2 res = d2d[index()];
	Vector2 x;
    x[0] = res.x;
    x[1] = res.y;
    return x;
}


void InputData::resetUnknown3d(BarrierIndex index) {
    assert(index.hasType(M3));
    assert(index() < f3d.size());
    f3d[index()].x = 0;
    f3d[index()].y = 0;
    f3d[index()].z = 0;
}

void InputData::setUnknown3d(BarrierIndex index, Scalar3 _f3d) {
    assert(index.hasType(M3));
    assert(index() < f3d.size());
    f3d[index()].x = _f3d.x;
    f3d[index()].y = _f3d.y;
    f3d[index()].z = _f3d.z;
}

void InputData::resetUnknown2d(BarrierIndex index) {
    assert(index.hasType(M2));
    assert(index() < f2d.size());
    f2d[index()].x = 0;
    f2d[index()].y = 0;
}

void InputData::allocate() {
    int size = dim();
    gpuErrchk(cudaMalloc(&cf3d, size*sizeof(Scalar3)));
    gpuErrchk(cudaMalloc(&cf2d, size*sizeof(Scalar2)));
    gpuErrchk(cudaMalloc(&cf1d, size*sizeof(Scalar)));

    gpuErrchk(cudaMalloc(&cdwd, dwd.size()*sizeof(Scalar)*LOCAL_SIZE));
    gpuErrchk(cudaMalloc(&cd3d, d3d.size()*sizeof(Scalar3)));
    gpuErrchk(cudaMalloc(&cd2d, d2d.size()*sizeof(Scalar2)));
    gpuErrchk(cudaMalloc(&cd1d, d1d.size()*sizeof(Scalar)));

    cudaMemset(cf3d, 0, size*sizeof(Scalar3));
    cudaMemset(cf2d, 0, size*sizeof(Scalar2));
    cudaMemset(cf1d, 0, size*sizeof(Scalar));

    cudaMemset(cdwd, 0, dwd.size()*sizeof(Scalar)*LOCAL_SIZE);
    cudaMemset(cd3d, 0, d3d.size()*sizeof(Scalar3));
    cudaMemset(cd2d, 0, d2d.size()*sizeof(Scalar2));
    cudaMemset(cd1d, 0, d1d.size()*sizeof(Scalar));
}

void InputData::release() {
    cudaFree(cf3d);
    cudaFree(cf2d);
    cudaFree(cf1d);
    cudaFree(cdwd);
    cudaFree(cd3d);
    cudaFree(cd2d);
    cudaFree(cd1d);
}

void InputData::copy() {
    std::cout << "Start Copy Data" << std::endl;
    cudaMemcpy(cf3d, f3d.data(), f3d.size()*sizeof(Scalar3), cudaMemcpyHostToDevice);
    cudaMemcpy(cf2d, f2d.data(), f2d.size()*sizeof(Scalar2), cudaMemcpyHostToDevice);
    cudaMemcpy(cf1d, f1d.data(), f1d.size()*sizeof(Scalar), cudaMemcpyHostToDevice);

    assert(dwd[1][0] == dwd[0][LOCAL_SIZE]);

    cudaMemcpy(cdwd, dwd.data(), dwd.size()*sizeof(Scalar)*LOCAL_SIZE, cudaMemcpyHostToDevice);
    cudaMemcpy(cd3d, d3d.data(), d3d.size()*sizeof(Scalar3), cudaMemcpyHostToDevice);
    cudaMemcpy(cd2d, d2d.data(), d2d.size()*sizeof(Scalar2), cudaMemcpyHostToDevice);
    cudaMemcpy(cd1d, d1d.data(), d1d.size()*sizeof(Scalar), cudaMemcpyHostToDevice);
    std::cout << "End Copy Data" << std::endl;
}

void InputData::retrieve() {
    std::cout << "Start Retrieving Data" << std::endl;
    cudaMemcpy(f3d.data(), cf3d, f3d.size()*sizeof(Scalar3), cudaMemcpyDeviceToHost);
    cudaMemcpy(f2d.data(), cf2d, f2d.size()*sizeof(Scalar2), cudaMemcpyDeviceToHost);
    cudaMemcpy(f1d.data(), cf1d, f1d.size()*sizeof(Scalar), cudaMemcpyDeviceToHost);

    cudaMemcpy(dwd.data(), cdwd, dwd.size()*sizeof(Scalar)*LOCAL_SIZE, cudaMemcpyDeviceToHost);
    cudaMemcpy(d1d.data(), cd1d, d1d.size()*sizeof(Scalar), cudaMemcpyDeviceToHost);
    cudaMemcpy(d2d.data(), cd2d, d2d.size()*sizeof(Scalar2), cudaMemcpyDeviceToHost);
    cudaMemcpy(d3d.data(), cd3d, d3d.size()*sizeof(Scalar3), cudaMemcpyDeviceToHost);
    std::cout << "End Retrieving Data" << std::endl;
}
void InputData::getDims(unsigned int* dims) {
    dims[0] = dim();
    dims[1] = dwd.size();
    dims[2] = d3d.size();
    dims[3] = d2d.size();
    dims[4] = d1d.size();

}

void InputData::assertDims(unsigned int* dims) {
    assert(dims[0] == dim());
    assert(dims[1] == dwd.size());
    assert(dims[2] == d3d.size());
    assert(dims[3] == d2d.size());
    assert(dims[4] == d1d.size());
}
unsigned int InputData::dim() const {
    std::cout << "Size: " <<
              f3d.size() << " " <<
              f2d.size() << " " <<
              f1d.size() << " " <<
              dwd.size() << " " <<
              d3d.size() << " " <<
              d2d.size() << " " <<
              d1d.size() << " " << std::endl;
    //unsigned long a = std::max(fpd.size(),f3d.size());
    unsigned long a = f3d.size();
    a = std::max(a,f2d.size());
    a = std::max(a,f1d.size());
    if((unsigned int) a != a) {
        assert(false);
    } else {
        std::cout << "Dimsize: " << a << std::endl;
    }
    return a;
}
void InputData::assign(void** result) {
    int idx = 0;
    result[idx++] = cf3d;
    result[idx++] = cf2d;
    result[idx++] = cf1d;
    result[idx++] = cdwd;
    result[idx++] = cd3d;
    result[idx++] = cd2d;
    result[idx++] = cd1d;
}
void InputData::read(Matrix3X& A, BarrierIndex start, int size, bool backup) {
    assert(start.hasType(M3));
    assert(start() + size <= f3d.size());
    if(!backup) {
        for(int i = 0; i < size; i++) {
            A(0,i) = f3d[start()+i].x;
            A(1,i) = f3d[start()+i].y;
            A(2,i) = f3d[start()+i].z;
        }
    } else {
        //std::cout << start << " " << size << " " << backup_f3d.size() << std::endl;
        for(int i = 0; i < size; i++) {
            A(0,i) = backup_f3d[start()+i].x;
            A(1,i) = backup_f3d[start()+i].y;
            A(2,i) = backup_f3d[start()+i].z;
        }
    }
}

void InputData::read(Matrix2X& A, BarrierIndex start, int size) {
    assert(start.hasType(M2));
    assert(start() + size <= f2d.size());
    for(int i = 0; i < size; i++) {
        A(0,i) = f2d[start()+i].x;
        A(1,i) = f2d[start()+i].y;
    }
}

Vector3 InputData::readVector3(BarrierIndex start, bool backup) {
	assert(start.hasType(M3));
	Vector3 x;
	if(!backup) {
		x[0] = f3d[start()].x;
		x[1] = f3d[start()].y;
		x[2] = f3d[start()].z;
	} else {
		x[0] = backup_f3d[start()].x;
        x[1] = backup_f3d[start()].y;
        x[2] = backup_f3d[start()].z;
	}
	return x;
}

Vector2 InputData::readVector2(BarrierIndex start) {
	assert(start.hasType(M2));
	Vector2 x;
	x[0] = f2d[start()].x;
	x[1] = f2d[start()].y;
	return x;
}

void InputData::readPoints(Matrix3X& A, BarrierIndex start, int size) {
    assert(start.hasType(D3));
    assert(start() + size <= d3d.size());
    for(int i = 0; i < size; i++) {
        A(0,i) = d3d[start()+i].x;
        A(1,i) = d3d[start()+i].y;
        A(2,i) = d3d[start()+i].z;
    }
}

void InputData::readPoints(Matrix2X& A, BarrierIndex start, int size) {
    assert(start.hasType(D2));
    assert(start() + size <= d2d.size());
    for(int i = 0; i < size; i++) {
        A(0,i) = d2d[start()+i].x;
        A(1,i) = d2d[start()+i].y;
    }
}


void InputData::clamp2d(BarrierIndex start, int size) {
    assert(start.hasType(M2));
    for(int i = 0; i < size; i++) {
        f2d[start()+i].x = fmin(fmax(f2d[start()+i].x,0),1);
        f2d[start()+i].y = fmin(fmax(f2d[start()+i].y,0),1);
    }
}

void InputData::write(const VectorX& A, BarrierIndex start, int size) {
    assert(start.hasType(M1));
    assert(start() + size <= f1d.size());
    for(int i = 0; i < size; i++) {
        f1d[start()+i] = A(i);
    }
}

void InputData::readPoints(VectorX& A, BarrierIndex start, int size) {
    assert(start.hasType(D1));
    assert(start() + size <= d1d.size());
    for(int i = 0; i < size; i++) {
        A(i) = d1d[start()+i];
    }
}


void InputData::read(VectorX& A, BarrierIndex start, int size) {
    assert(start.hasType(M1));
    assert(start() + size <= f1d.size());
    for(int i = 0; i < size; i++) {
        A(i) = f1d[start()+i];
    }
}

Scalar InputData::read(BarrierIndex i) {
    assert(i.hasType(M1));
    assert(i() < f1d.size());
    return f1d[i()];
}

Scalar InputData::readDynamic1d(BarrierIndex i) {
    assert(i.hasType(M1) || i.hasType(D1));
    return i.hasType(M1) ? read(i) : readPoint1d(i);
}



void InputData::readDynamic(Matrix3X& A, BarrierIndex start, int size) {
    assert(start.hasType(D3) || start.hasType(M3));
    if(start.hasType(D3)) {
        readPoints(A, start, size);
    } else {
        read(A, start, size);
    }
}

void InputData::readDynamic(VectorX& A, BarrierIndex start, int size) {
    assert(start.hasType(D1) || start.hasType(M1));
    if(start.hasType(D1)) {
        readPoints(A, start, size);
    } else {
        read(A, start, size);
    }
}

void InputData::write(Scalar x, BarrierIndex i) {
    assert(i.hasType(M1));
    assert(i() < f1d.size());
    f1d[i()] = x;
}

#undef gpuErrchk
