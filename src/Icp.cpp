#include "Icp.h"

Icp::Icp(Matrix3X cloud) {
	data.resize(cloud.cols(), 3);

	for (int j = 0; j < cloud.cols(); j++) {
		data(j,0) = cloud(0,j);
		data(j,1) = cloud(1,j);
		data(j,2) = cloud(2,j);
	}
	mat_index.reset(new my_kd_tree_t(data, 10));
	mat_index->index->buildIndex();	
}

int Icp::findClosestPoint(Vector3 query) {
	Scalar query_fuse[3];
	query_fuse[0] = query[0];
	query_fuse[1] = query[1];
	query_fuse[2] = query[2];
	
	Scalar out_dists_sqr = -1.0f;
	size_t ret_index = 0;

	nanoflann::KNNResultSet<Scalar> resultSet(1);
	resultSet.init(&ret_index, &out_dists_sqr);
	mat_index->index->findNeighbors(resultSet, &query_fuse[0],
		nanoflann::SearchParams(10));

	return ret_index;
}

