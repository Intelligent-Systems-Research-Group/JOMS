#include "MeshTopology.h"
#include "types.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>
#include <nanoflann.hpp>
#include <memory>
#include <stdlib.h>

void makeCube(MeshTopology* mesh, Matrix3X* verts)
{
  //Initial mesh - A cube/parallelepiped
  size_t num_verts = 8;
  size_t num_faces = 6;
  mesh->num_vertices = num_verts;

  // Init vertices
  verts->resize(3, num_verts);
  verts->col(0) << -1.0, -1.0, +1.0;
  verts->col(1) << +1.0, -1.0, +1.0;
  verts->col(2) << -1.0, +1.0, +1.0;
  verts->col(3) << +1.0, +1.0, +1.0;
  verts->col(4) << -1.0, +1.0, -1.0;
  verts->col(5) << +1.0, +1.0, -1.0;
  verts->col(6) << -1.0, -1.0, -1.0;
  verts->col(7) << +1.0, -1.0, -1.0;

  //Fill the vertices per face
  mesh->quads.resize(5, num_faces);
  mesh->quads.col(0) << 0, 1, 3, 2, -1;
  mesh->quads.col(1) << 2, 3, 5, 4, -1;
  mesh->quads.col(2) << 4, 5, 7, 6, -1;
  mesh->quads.col(3) << 6, 7, 1, 0, -1;
  mesh->quads.col(4) << 1, 7, 5, 3, -1;
  mesh->quads.col(5) << 6, 0, 2, 4, -1;

  mesh->update_adjacencies();
}

void makeCube2(MeshTopology* mesh, Matrix3X* verts) {
	//Initial mesh - A cube/parallelepiped
	size_t num_verts = 8;
	size_t num_faces = 7;
	mesh->num_vertices = num_verts;

	// Init vertices
	verts->resize(3, num_verts);
	verts->col(0) << -1.0, -1.0, +1.0;
	verts->col(1) << +1.0, -1.0, +1.0;
	verts->col(2) << -1.0, +1.0, +1.0;
	verts->col(3) << +1.0, +1.0, +1.0;
	verts->col(4) << -1.0, +1.0, -1.0;
	verts->col(5) << +1.0, +1.0, -1.0;
	verts->col(6) << -1.0, -1.0, -1.0;
	verts->col(7) << +1.0, -1.0, -1.0;

	//Fill the vertices per face
	mesh->quads.resize(4, num_faces);
	mesh->quads.col(0) << 0, 1, 3, 2;
	mesh->quads.col(1) << 2, 3, 5, 4;
	mesh->quads.col(2) << 4, 5, 7, 6;
	mesh->quads.col(3) << 6, 7, 1, 0;
	mesh->quads.col(4) << 1, 7, 5, 3;
	mesh->quads.col(5) << 6, 0, 2, -1;
	mesh->quads.col(6) << 6, 2, 4, -1;

	mesh->update_adjacencies();
}

void makeTriangle(MeshTopology* mesh, Matrix3X* verts) {
	//Initial mesh - A cube/parallelepiped
	size_t num_verts = 8;
	size_t num_faces = 12;
	mesh->num_vertices = num_verts;

	// Init vertices
	verts->resize(3, num_verts);
	verts->col(0) << -1.0, -1.0, +1.0;
	verts->col(1) << +1.0, -1.0, +1.0;
	verts->col(2) << -1.0, +1.0, +1.0;
	verts->col(3) << +1.0, +1.0, +1.0;
	verts->col(4) << -1.0, +1.0, -1.0;
	verts->col(5) << +1.0, +1.0, -1.0;
	verts->col(6) << -1.0, -1.0, -1.0;
	verts->col(7) << +1.0, -1.0, -1.0;

	//Fill the vertices per face
	mesh->quads.resize(4, num_faces);
	mesh->quads.col(0) << 0, 1, 3, -1;
	mesh->quads.col(1) << 0, 3, 2, -1;
	mesh->quads.col(2) << 2, 3, 5, -1;
	mesh->quads.col(3) << 2, 5, 4, -1;
	mesh->quads.col(4) << 4, 5, 7, -1;
	mesh->quads.col(5) << 4, 7, 6, -1;
	mesh->quads.col(6) << 6, 7, 1, -1;
	mesh->quads.col(7) << 6, 1, 0, -1;
	mesh->quads.col(8) << 1, 7, 5, -1;
	mesh->quads.col(9) << 1, 5, 3, -1;
	mesh->quads.col(10) << 6, 0, 2, -1;
	mesh->quads.col(11) << 6, 2, 4, -1;

	mesh->update_adjacencies();
}

void makePath(MeshTopology* mesh, Matrix3X* verts) {
	mesh->num_vertices = verts->cols();
	mesh->quads.resize(5, mesh->num_vertices - 2);
	mesh->quads.setConstant(-1);
	mesh->nverts_per_face.resize(mesh->num_vertices - 2);
	for (int j = 0; j < mesh->num_vertices - 2; j++) {
		mesh->quads(0, j) = j;
		mesh->quads(1, j) = j+1;
		mesh->quads(2, j) = j+2;
		mesh->nverts_per_face[j] = 3;
	}
}

bool saveObj(std::string path, const MeshTopology* mesh, const Matrix3X* verts) {
	std::filebuf fb;
	fb.open(path, std::ios::out);
	std::ostream outfile(&fb);
	int nVerts = mesh->num_vertices;
	int nFaces = mesh->num_faces();
	for (int i = 0; i < nVerts; i++) {
		outfile << "v ";
		for(int j = 0; j < 3; j++) {
			outfile << (*verts)(j,i) << " ";
		}
		outfile << "\n";
	}
	for (int i = 0; i < nFaces; i++) {
		outfile << "f ";
		//std::cout << mesh->nverts_per_face[i] << std::endl;
		for (int j = 0; j < mesh->nverts_per_face[i]; j++) {
			outfile << (mesh->quads(j,i) + 1) << " ";
		}
		outfile << "\n";
	}
	return true;
}
bool loadObj(std::string path, MeshTopology* mesh, Matrix3X* verts, Matrix3X* normals, double scale, bool mirror, bool rot) {
	std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cout << "reading failed for: " << path << std::endl;
        return false;
    } else {
        //std::cout << "start reading: " << path << std::endl;
    }
	
	std::string line;
	int num_verts = 0;
    int num_normals = 0;
	int num_faces = 0;
	int vert_idx = 0;
	int normal_idx = 0;
	int face_idx = 0;

	while (std::getline(infile, line))
	{
        //std::cout << line << std::endl;
		if (line == "" || line == "\r") continue;
		std::istringstream iss(line);
		std::string prefix;
		if (!(iss >> prefix)) {
            for(int i = 0; i < line.size(); i++) {
                std::cout << ((int)line[i]) << std::endl;
            }
			std::cout << line << " error while loading file" << std::endl;
            assert(false);
			return false;
            continue;
		}
		if (prefix == "v") num_verts++;
        else if(prefix == "vn") num_normals++;
		else if (prefix == "f") num_faces++;
		else continue;
	}
    std::cout << "Verts: " << num_verts << " Normals: " << num_normals << std::endl;
    assert(normals == NULL || num_verts == num_normals);

	verts->resize(3, num_verts);
	if(normals != NULL)
		normals->resize(3, num_verts);
    if(mesh != NULL) {
	    mesh->num_vertices = num_verts;
	    mesh->quads.resize(5, num_faces);
    }

	infile = std::ifstream(path);
	while (std::getline(infile, line))
	{
		if (line == "" || line == "\r") continue;
		std::istringstream iss(line);
		std::string prefix;
		if (!(iss >> prefix)) {
			std::cout << "error while loading file for real" << std::endl;
            return false;
		}
		if (prefix == "v") {
			Scalar x, y, z;
			if (!(iss >> x >> y >> z)) {
				std::cout << "error while parsing vertex" << std::endl;
				return false;
			}
			if(!mirror) {
				if (rot) {
					verts->col(vert_idx) << scale*x, scale*z, -scale*y;
				}
				else {

					verts->col(vert_idx) << scale*x, scale*y, scale*z;
				}
			}
			else {
				if (rot) {
					verts->col(vert_idx) << -scale*x, scale*z, -scale*y;
				}
				else {

					verts->col(vert_idx) << -scale*x, scale*y, scale*z;
				}
			}
			vert_idx++;
		} else if(prefix == "vn") {
			if(normals != NULL) {
				Scalar x, y, z;
				if (!(iss >> x >> y >> z)) {
					std::cout << "error while parsing normal" << std::endl;
					return false;
				}
				if(!mirror) {
					if (rot) {
						normals->col(normal_idx) << x, z, -y;
					}
					else {
						normals->col(normal_idx) << x, y, z;
					}
				} else {
					if (rot) {
						normals->col(normal_idx) << -x, z, -y;
					}
					else {
						normals->col(normal_idx) << -x, y, z;
					}
				}
				normal_idx++;
			}
		} else if (prefix == "f" && mesh != NULL) {
			std::string token;
			std::string tok;
			//std::istringstream linestream(std::string(line.begin()+2,line.end()));
			int vert_idxs[5] = {-1, -1, -1, -1, -1};
			for (int i = 0; i < 5; i++) {
				if (!(iss >> tok)) {
					if (i < 3) {
						std::cout << "error while parsing face" << std::endl;
						return false;
					}
					break;
				}

				std::istringstream linestream(tok);
				while (std::getline(linestream, token, '/')) {
					std::istringstream tokenstream(token);
					tokenstream >> vert_idxs[i];
					vert_idxs[i]--;

                    
					break;
				}
			}

			mesh->quads.col(face_idx) << vert_idxs[0], vert_idxs[1], vert_idxs[2], vert_idxs[3], vert_idxs[4];
			face_idx++;
		}
		else continue;
	}
    if(mesh != NULL)
	    mesh->update_adjacencies();
	return true;
}

void MeshTopology::createEdgeList(IMatrix2X& edges, IMatrix4X* faceAdj) {
	std::set<std::set<size_t>> uedges;
	std::map<std::set<size_t>, std::set<size_t>> dual;
	if(faceAdj != NULL)
		faceAdj->resize(4, quads.cols());


	for (int i = 0; i < quads.cols(); i++) {
		size_t verts_per_face = nverts_per_face[i];
		for (int j = 0; j < verts_per_face; j++) {
			//for (int k = j+1; k < nverts_per_face[i] - 1; k++) {
			std::set<size_t> pair; 
			int k = (j + 1) % (verts_per_face);
			size_t from = quads(j, i);
			size_t to = quads(k, i);
			//if (from == 0 || to == 0) {
			//	std::cout << from << ":" << to << " " << i << ":" << verts_per_face << std::endl;
			//}
			pair.insert(from);
			pair.insert(to);
			uedges.insert(pair);
			if(faceAdj != NULL)
				dual[pair].insert(i);
		}
	}

	if (faceAdj != NULL) {
		for (int i = 0; i < quads.cols(); i++) {
			size_t verts_per_face = nverts_per_face[i];
			for (int j = 0; j < verts_per_face; j++) {
				std::set<size_t> pair;
				int k = (j + 1) % (verts_per_face);
				size_t from = quads(j, i);
				size_t to = quads(k, i);

				pair.insert(from);
				pair.insert(to);

				auto& faceEdges = dual[pair];
				for (auto it = faceEdges.begin(); it != faceEdges.end(); it++) {
					if (*it == i) continue;
					(*faceAdj)(j, i) = *it;
				}
			}
		}
	}

	edges.resize(2, uedges.size());
	size_t i = 0;
	for (auto it = uedges.begin(); it != uedges.end(); it++) {
		auto pt = it->begin();
		edges(0, i) = *pt;
		pt++;
		edges(1, i) = *pt;
		i++;
	}
}
void MeshTopology::vertexAdjacency(std::map<size_t, std::set<size_t>>& adjacency, IMatrix4X* faceAdj) {
	IMatrix2X edges;
	createEdgeList(edges, faceAdj);
	for (size_t i = 0; i < edges.cols(); i++) {
		size_t from = edges(0, i);
		size_t to = edges(1, i);
		/*
		if (from == 0) {
			std::cout << from << ":" << to << std::endl;
		}
		if (to == 0) {
			std::cout << from << ":" << to << std::endl;
		}
		*/
		adjacency[from].insert(to);
		adjacency[to].insert(from);
	}
}

void MeshTopology::initVertexToSurfaceCoordinate(int vidx, int& fidx, double &u, double& v) {
	double e = 0.001;
	for (int i = 0; i < quads.cols(); i++) {
		for (int j = 0; j < 4; j++) {
			if (quads(j, i) == vidx) {
				fidx = i;
				switch (j) {
				case 0:
					u = 0 + e;
					v = 0 + e;
					break;
				case 1:
					u = 1 - e;
					v = 0 + e;
					break;
				case 2:
					u = 1 - e;
					v = 1 - e;
					break;
				case 3:
					u = 0 + e;
					v = 1 - e;
					break;
				}
				return;
			}
		}
	}
	return;
}
void MeshTopology::vertexToSurfaceCoordinate(int vidx, int& fidx, double &u, double& v) {
	fidx = texcoords[vidx].faceId;
	u = texcoords[vidx].s;
	v = texcoords[vidx].t;
}

std::map<std::pair<size_t, size_t>, double> MeshTopology::calculateLaplaceWeights() {
	std::map<std::pair<size_t, size_t>, double> w;
	std::map<size_t, std::set<size_t>> adjacency;
	vertexAdjacency(adjacency);

	for (size_t i = 0; i < quads.cols(); i++) {
		auto& neighbours = adjacency.at(i);
		for (auto it = neighbours.begin(); it != neighbours.end(); it++) {
			size_t vertex_index = *it;
			std::pair<size_t, size_t> key = std::make_pair(i, vertex_index);
		}
	}

	return w;
}


void MeshTopology::update_adjacencies()
{
	//std::cout << quads.transpose() << std::endl;
	compressed_faces.clear();
	nverts_per_face.clear();
	for (size_t f = 0; f < num_faces(); f++)
	{
		size_t nverts = 0;
		for (size_t k = 0; k < 5; k++)
		{
			size_t idx = quads(k, f);
			if (idx != -1) {
				compressed_faces.push_back(idx);
				nverts++;
			}
		}
		nverts_per_face.push_back(nverts);
	}
  //Find the adjacent faces to every face
  face_adj.resize(4, num_faces());
  face_adj.fill(-1);

  for (size_t f = 0; f < num_faces(); f++) {
	  if (nverts_per_face[f] != 4) return;
  }

  IMatrix2X edges;
  createEdgeList(edges, &face_adj);

  for (int i = 0; i < num_vertices/*nverts_per_face.size()*/; i++) {
	  texcoord coord;
	  initVertexToSurfaceCoordinate(i, coord.faceId, coord.s, coord.t);
	  texcoords.push_back(coord);
  }
  
  /*
  for (size_t f = 0; f < num_faces(); f++)
    for (size_t k = 0; k < (int) (nverts_per_face[f]); k++)
    {
      // Find kth edge 
      unsigned int kinc = (int(k) + 1) % ((int)(nverts_per_face[f]));
      int edge[2] = { quads(k,f), quads(kinc,f) };

      // And find the face that shares its reverse
      int found = 0;
      int other = -1;

      for (size_t fa = 0; fa < num_faces(); fa++)
      {
		if (f == fa) continue;
        for (size_t l = 0; l < (int)(nverts_per_face[fa]); l++)
          if ((quads(l, fa) == edge[1]) && (quads((l + 1) % (int)(nverts_per_face[fa]), fa) == edge[0])) {
            other = (int) fa;
            found++;
          }
      }
      assert(found == 1);

      face_adj(k, f) = other;
    }
	*/
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	for (int i = 0; i < v.size(); i++) idx[i] = i;

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}

struct edge { int to; double length; };

std::vector<double> dijkstra(const std::vector< std::vector<edge> > &graph, int source) {
	std::vector<double> min_distance(graph.size(), INT_MAX);
	min_distance[source] = 0;
	std::set< std::pair<int, int> > active_vertices;
	active_vertices.insert({ 0,source });

	while (!active_vertices.empty()) {
		int where = active_vertices.begin()->second;
		active_vertices.erase(active_vertices.begin());
		for (auto ed : graph[where])
			if (min_distance[ed.to] > min_distance[where] + ed.length) {
				active_vertices.erase({ min_distance[ed.to], ed.to });
				min_distance[ed.to] = min_distance[where] + ed.length;
				active_vertices.insert({ min_distance[ed.to], ed.to });
			}
	}
	return min_distance;
}

IMatrix4X MeshTopology::findClosest(const Matrix3X& verts, const std::vector<std::vector<int>>& vgroups, const SparseMatrix& sparseWeights,  Matrix4X* dist) {
	IMatrix4X output(4, verts.cols());
	if (dist != NULL) {
		dist->resize(4, verts.cols());
	}
	std::vector<std::vector<double>> distances(verts.cols());
	std::vector< std::vector<edge> > graph(verts.cols());
	IMatrix2X edges;
	createEdgeList(edges, NULL);
	for (int i = 0; i < edges.cols(); i++) {
		edge e1, e2;
		int from = edges(0, i);
		int to = edges(1, i);
		e1.to = to;
		e1.length = (verts.col(to) - verts.col(from)).norm();
		e2.to = from;
		e2.length = e1.length;
		graph[from].push_back(e1);
		graph[to].push_back(e2);
	}

	for (int i = 0; i < vgroups.size(); i++) {
		graph.push_back(std::vector<edge>());
		int helperIdx = graph.size() - 1;
		for (int j = 0; j < vgroups[i].size(); j++) {
			edge e;
			e.to = vgroups[i][j];
			e.length = 0;
			graph[helperIdx].push_back(e);
		}
	}

	for (int i = 0; i < vgroups.size(); i++) {
		auto dist = dijkstra(graph, verts.cols() + i);
		for(int j = 0; j < verts.cols(); j++) {
            auto distance = dist[j];
            if(sparseWeights.find(IndexPair(j,i)) != sparseWeights.end()) {
                distance -= 1000000;
            }
			distances[j].push_back(distance);
		}
	}
	for (int i = 0; i < verts.cols(); i++) {
		std::vector<size_t> sorted = sort_indexes(distances[i]);
		for (int j = 0; j < 4; j++) {
			output(j, i) = sorted[j];
			if (dist != NULL) {
				(*dist)(j, i) = distances[i][output(j, i)];
			}
		}
	}

	return output;
}
const static Eigen::IOFormat XYZFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n");
void dumpCloud(std::string name, int iter, const Matrix3X& cloud, const Matrix3X* normals) {
	std::ofstream file(name + std::to_string(iter) + std::string(".xyz"));
	if(normals == NULL) {
		file << cloud.transpose().format(XYZFormat);
	}
	else {
		MatrixX M(6, cloud.cols());
		M.block(0, 0, 3, cloud.cols()) = cloud;
		M.block(3, 0, 3, cloud.cols()) = *normals;
		file << M.transpose().format(XYZFormat);
	}
}


bool readCloud(std::string path, Matrix3X* cloud, Matrix3X* normals,double scale, bool mirror, bool rot) {
    std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cout << "reading failed for: " << path << std::endl;
    } else {
        //std::cout << "start reading: " << path << std::endl;
    }
	
	std::string line;
	int num_verts = 0;

	while (std::getline(infile, line))
	{
        if (line.find("nan") == std::string::npos) {
			
            num_verts++;
        } else {
			std::cout << "NAN found" << std::endl;
			exit(1);
		}
	}
    cloud->resize(3,num_verts);
    if(normals != NULL)
		normals->resize(3,num_verts);

    int idx = 0;
    
    std::ifstream infile2(path);
    if (!infile2.is_open()) {
        std::cout << "reading failed for: " << path << std::endl;
    } else {
        //std::cout << "start reading: " << path << std::endl;
    }
    while (std::getline(infile2, line)) {
       //std::cout << line << std::endl;
		if(line == "") break;
        if(line.find("nan") != std::string::npos) {
            std::cout << "NAN found" << std::endl;
            exit(1);
            continue;
        }
		std::istringstream iss(line);

        Scalar x, y, z, nx,ny,nz;
	    if (!(iss >> x >> y >> z >> nx >> ny >> nz)) {            
		    std::cout << line << " error while loading file" << std::endl;
            assert(false);
		    return false;
	    }
	    /*
        if(nx == 0 && ny == 0 && nz == 0) {
            num_verts--;
            continue;
        }*/

		if(!mirror) {
			if (rot) {
				cloud->col(idx) << scale*x, scale*z, -scale*y;
			}
			else {

				cloud->col(idx) << scale*x, scale*y, scale*z;
			}
		}
		else {
			if (rot) {
				cloud->col(idx) << -scale*x, scale*z, -scale*y;
			}
			else {

				cloud->col(idx) << -scale*x, scale*y, scale*z;
			}
		}
		if(normals != NULL) {
			Scalar len = sqrt(nx*nx+ny*ny+nz*nz);
			if(len > 0) {
				nx /= len;
				ny /= len;
				nz /= len;
			}
			
			if(!mirror) {
				if (rot) {
					normals->col(idx) << nx, nz, -ny;
				}
				else {
					normals->col(idx) << nx, ny, nz;
				}
			} else {
				if (rot) {
					normals->col(idx) << -nx, nz, -ny;
				}
				else {
					normals->col(idx) << -nx, ny, nz;
				}
			}
		}
        idx++;
		
	}
    cloud->conservativeResize(3,num_verts);
    if(normals != NULL) {
		normals->conservativeResize(3,num_verts);
	}
    assert(is_sane(*cloud));
    assert(normals == NULL || is_sane(*normals));
    return true;
}

bool readMatrix(std::string path, MatrixX* matrix) {
	std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cout << "reading failed for: " << path << std::endl;
    } else {
        //std::cout << "start reading: " << path << std::endl;
    }
    
    std::string line;
	int rows = 0;
	int cols = 3;

	while (std::getline(infile, line))
	{
        if (line.find("nan") == std::string::npos) {
            rows++;
        }
	}
	
	matrix->resize(rows,cols);
	
	std::ifstream infile2(path);
    if (!infile2.is_open()) {
        std::cout << "reading failed for: " << path << std::endl;
    } else {
        //std::cout << "start reading: " << path << std::endl;
    }
    int idx = 0;
	while (std::getline(infile2, line)) {
       //std::cout << line << std::endl;
		if(line == "") break;
        if(line.find("nan") != std::string::npos) {
            std::cout << "NAN found" << std::endl;
            continue;
        }
		std::istringstream iss(line);

        Scalar x, y, z, nx,ny,nz, c;
	    if (!(iss >> x >> y >> z)) {            
		    std::cout << line << " error while loading file" << std::endl;
            assert(false);
		    return false;
	    }
        matrix->row(idx) << x, y, z, nx , ny, nz, c;
        idx++;
	}
	return true;
}

bool readVector(std::string path, std::vector<bool>* result) {
	std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cout << "reading failed for: " << path << std::endl;
        exit(1);
        return false;
    }
    std::string line;
    int sum = 0;
	while (std::getline(infile, line))
	{
        switch(line[0]) {
		case '0':
			result->push_back(false);
			break;
		case '1':
			result->push_back(true);
                        sum++;
			break;
		default:
			exit(1);        
			return false;
		}
	}
        std::cout << "Ground Truth point count: " << sum << std::endl;
	return true;

}

void filter(Matrix3X* data, Matrix3X* normal, std::vector<bool> mask) {
	int j = 0;
	for(int i = 0; i < data->cols(); i++) {
		bool validNormal = normal->col(i).norm() > 0;
		if(mask[i] && validNormal) {
			data->col(j) = data->col(i);
			normal->col(j) = normal->col(i);
			j++;
		}
	}
	data->conservativeResize(3, j);
	normal->conservativeResize(3,j);

	assert(is_sane(*data));
        assert(is_sane(*normal));
}

void filter(Matrix3X* data, Matrix3X* normal) {
	int j = 0;
	for(int i = 0; i < data->cols(); i++) {
		bool validNormal = normal->col(i).norm() > 0;
		if(validNormal) {
			data->col(j) = data->col(i);
			normal->col(j) = normal->col(i);
			j++;
		}
	}
	data->conservativeResize(3, j);
	normal->conservativeResize(3,j);

	assert(is_sane(*data));
    assert(is_sane(*normal));
}

void filterGroundPlane(Matrix3X* data, Matrix3X* normal, int axis, Scalar groundPlane) {
	int j = 0;
    for(int i = 0; i < data->cols(); i++) {
        bool overGround = (*data)(axis,i) > groundPlane;
        if(overGround) {
            data->col(j) = data->col(i);
            normal->col(j) = normal->col(i);
            j++;
        }
    }
    data->conservativeResize(3, j);
    normal->conservativeResize(3,j);

    assert(is_sane(*data));
    assert(is_sane(*normal));

}

void fill(Matrix3X* data, Matrix3X* normal, int n) {
	int m = data->cols();
	if(m >= n) return;
	
	data->conservativeResize(3, n);
	normal->conservativeResize(3,n);
	//std::cout << "m" << m <<"n"<<n <<std::endl;
	for(int i = m; i < n; i++) {
		int j = rand() % m;
		//std::cout << "i" << i <<"j"<<j<<"n"<<n <<std::endl;
		data->col(i) = data->col(j);
		normal->col(i) = normal->col(j);
	}
}

void findSymmetric(const Matrix3X& verts, std::vector<int>* sym) {
	float maxError = 0;
	for(int i = 0; i < verts.cols(); i++) {
		Vector3 query = verts.col(i);
		query[0] *= -1;
		int result = -1;
		Scalar dist = 1000000000;
		for(int j = 0; j < verts.cols(); j++) {
			Vector3 pt =  verts.col(j);
			if((pt-query).stableNorm() < dist) {
				result = j;
				dist = (pt-query).stableNorm();
			}
		}
		//std::cout << i << " <-> " << result << " | " << dist << std::endl;
		assert(result >= 0);
		assert(dist <= 0.01);
		maxError = fmax(maxError,dist);
		(*sym)[i] = result;
	}
	//std::cout << maxError << std::endl;
	//exit(1);
}

void filter(Matrix3X* data, Matrix3X* normal, const Matrix3X& model) {
	assert(data->cols() == normal->cols());
	typedef nanoflann::KDTreeEigenMatrixAdaptor<MatrixX>  my_kd_tree_t;
	std::unique_ptr<my_kd_tree_t> mat_index;
	MatrixX input = model.transpose();
	mat_index.reset(new my_kd_tree_t(input, 10));
	mat_index->index->buildIndex();
	
	Scalar out_dists_sqr = -1.0f;
	size_t   ret_index = 0;
	int j = 0;
	for(int i = 0; i < data->cols(); i++) {
		nanoflann::KNNResultSet<Scalar> resultSet(1);
		resultSet.init(&ret_index, &out_dists_sqr);
		Scalar query_fuse[3];
		query_fuse[0] = (*data)(0,i); 
		query_fuse[1] = (*data)(1,i);
		query_fuse[2] = (*data)(2,i);
		mat_index->index->findNeighbors(resultSet, &query_fuse[0], nanoflann::SearchParams(10));
		if(out_dists_sqr < 0.01) {
			data->col(j) = data->col(i);
			normal->col(j) = normal->col(i);
			j++;
		}
	}
	data->conservativeResize(3, j);
	normal->conservativeResize(3,j);
}

void corruptNonUniform(Matrix3X* data, Matrix3X* normal) {
	int originIdx = rand() % data->cols();
	Vector3 origin = data->col(originIdx);
	const Scalar radius = .5;
	Matrix3X result; result.resize(3, data->cols());
	Matrix3X resultNormal; resultNormal.resize(3, normal->cols());
	int j = 0;
	for(int i = 0; i < data->cols(); i++) {
		if((data->col(i) - origin).norm() > radius) {
			result.col(j) = data->col(i);
			resultNormal.col(j) = normal->col(i);
			j++;
		}
	}
	result.conservativeResize(3,j);
	resultNormal.conservativeResize(3,j);
	*data = result;
	*normal = resultNormal;
}
std::vector<int> greedyFarthestSampling(const Matrix3X& data, const Matrix3X& normal, int k, Matrix3X* out, Matrix3X* nout) {
	int n = data.cols();
	std::vector<int> r;
	r.resize(k);
	for(int j = 0; j < k; j++) {
		r[j] = rand() % n;
	}
	
	for(int i = 0; i< 1000; i++) {
		int jidx = rand() % k;
		int oidx = r[jidx];
		int nidx = rand() % n;
		double cold = 0;
		double cnew = 0;
		for(int j = 0; j < k; j++) {
			if(jidx == j) continue;
			cold += (data.col(oidx) - data.col(j)).norm();
			cnew += (data.col(nidx) - data.col(j)).norm();
		}
		if(cnew > cold) r[jidx] = nidx;
	}
	out->resize(3,k);
	nout->resize(3,k);
	for(int j = 0; j < k; j++) {
		out->col(j) = data.col(r[j]);
		nout->col(j) = normal.col(r[j]);
	}
	return r;
}

void quiver(std::string path, const Matrix3X& X, const Matrix3X& Y) {
	Matrix3X X1 = X;
	X1.row(0) = X1.row(0).array() + 0.001f;
	Matrix3X X2 = X;
    X2.row(0) = X2.row(0).array() - 0.001f;
	Matrix3X Y1 = Y;
    Y1.row(0) = Y1.row(0).array() + 0.0005f;
    Matrix3X Y2 = Y;
    Y2.row(0) = Y2.row(0).array() - 0.0005f;
	Matrix3X V; V.resize(3,4*X.cols());

	MeshTopology mesh;
	mesh.quads.resize(5, X.cols());
	mesh.quads.setConstant(-1);
	for(int i = 0; i < X.cols(); i++) {
		mesh.quads(0,i) = 4*i+0;
		mesh.quads(1,i) = 4*i+1;
		mesh.quads(2,i) = 4*i+2;
		mesh.quads(3,i) = 4*i+3;
		V.col(4*i+0) = X1.col(i);
		V.col(4*i+1) = X2.col(i);
		V.col(4*i+2) = Y2.col(i);
		V.col(4*i+3) = Y1.col(i);
	}
	mesh.num_vertices = V.cols();
	mesh.update_adjacencies();
	saveObj(path,&mesh,&V);
}
