#include "calculus.h"
#include <iostream>

Matrix32 principleCurvatureAmbient(Vector3 p, Vector3 dSdu, Vector3 dSdv, Vector3 dSduu, Vector3 dSduv, Vector3 dSdvv) {
    Matrix32 dSdx;
    dSdx.col(0) = dSdu;
    dSdx.col(1) = dSdv;
    //Matrix32 dSdxN = dSdx;
    //dSdxN.col(0) /= dSdxN.col(0).norm();
    //dSdxN.col(1) /= dSdxN.col(1).norm();

    Matrix22 I = dSdx.transpose()*dSdx;
    Vector3 n = dSdu.cross(dSdv);
    n = n / n.norm();
    Matrix22 II;
    II(0,0) = n.dot(dSduu);
    II(0,1) = II(1,0) = n.dot(dSduv);
    II(1,1) = n.dot(dSdvv);
    Matrix22 S = -I.inverse()*II;
    Eigen::EigenSolver<Matrix22> es;
    es.compute(S,true);
    Matrix32 principalDirections;
    principalDirections.col(0) = dSdx*es.eigenvectors().real().col(0); //* sqrt(es.eigenvalues().real()[0]
    principalDirections.col(1) = dSdx*es.eigenvectors().real().col(1); // * sqrt(es.eigenvalues().real()[1]);
    principalDirections.col(0) /= principalDirections.col(0).norm();
    principalDirections.col(1) /= principalDirections.col(1).norm();
    principalDirections.col(0) *= fabs(es.eigenvalues().real()[0]);
    principalDirections.col(1) *= fabs(es.eigenvalues().real()[1]);
    Scalar k1 = fabs(es.eigenvalues().real()[0]);
    Scalar k2 = fabs(es.eigenvalues().real()[1]);
    if(1E-6 > k1  || 1E-6 > k2) {
        std::cout << "extremely low curvature" << std::endl;
    }
    if(k1 > 1E6 || k2 > 1E6) {
        std::cout << "extremely high curvature" << std::endl;
    }
    return principalDirections;
}

Matrix22 principleCurvatureTangent(Vector3 p, Vector3 dSdu, Vector3 dSdv, Vector3 dSduu, Vector3 dSduv, Vector3 dSdvv) {
    Matrix32 dSdx;
    dSdx.col(0) = dSdu;
    dSdx.col(1) = dSdv;
    //Matrix32 dSdxN = dSdx;
    //dSdxN.col(0) /= dSdxN.col(0).norm();
    //dSdxN.col(1) /= dSdxN.col(1).norm();

    Matrix22 I = dSdx.transpose()*dSdx;
    Vector3 n = dSdu.cross(dSdv);
    n = n / n.norm();
    Matrix22 II;
    II(0,0) = n.dot(dSduu);
    II(0,1) = II(1,0) = n.dot(dSduv);
    II(1,1) = n.dot(dSdvv);
    Matrix22 S = -I.inverse()*II;
    Eigen::EigenSolver<Matrix22> es;
    es.compute(S,true);

    Matrix32 pr;
    pr.col(0) = dSdx*es.eigenvectors().real().col(0); //* sqrt(es.eigenvalues().real()[0]
    pr.col(1) = dSdx*es.eigenvectors().real().col(1); // * sqrt(es.eigenvalues().real()[1]);

    Matrix22 principalDirections;
    principalDirections.col(0) = es.eigenvectors().real().col(0); //* sqrt(es.eigenvalues().real()[0]
    principalDirections.col(1) = es.eigenvectors().real().col(1); // * sqrt(es.eigenvalues().real()[1]);
    principalDirections.col(0) /= pr.col(0).norm();
    principalDirections.col(1) /= pr.col(1).norm();
    //principalDirections = principalDirections.inverse().eval();
    return principalDirections;
}
