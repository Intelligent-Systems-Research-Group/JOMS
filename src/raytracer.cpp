#include "raytracer.h"
#include <iostream>

void error_handler(void* userPtr, RTCError code, const char* str) {
	if(code == RTC_ERROR_NONE) return;
	std::cerr << "Embree3 ErrorCode: " << code << std::endl;
	if(str != nullptr) {
		std::cerr << str << std::endl;
	}
}

bool memory_handler(void* userPtr, ssize_t bytes, bool post) {
	ssize_t* pcount = (ssize_t*) userPtr;
	if(post)
		*pcount += bytes;
	return true;
}

struct Triangle { int v0, v1, v2; };
struct Vertex   { float x,y,z,r;  };

unsigned int RayTracer::addModel (RTCScene scene_i, const Matrix3X& V, const MeshTopology& top) {
	RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
	
	//face_colors = (Vec3fa*) alignedMalloc(12*sizeof(Vec3fa));
	//vertex_colors = (Vec3fa*) alignedMalloc(8*sizeof(Vec3fa));
	int nfaces = 2*top.quads.cols();
	Vertex* vertices = (Vertex*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,sizeof(Vertex),V.cols()); 
	for(int i = 0; i < V.cols(); i++) {
		vertices[i].x = V(0,i);
		vertices[i].y = V(1,i);
		vertices[i].z = V(2,i);
	}
	Triangle* triangles = (Triangle*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,sizeof(Triangle),nfaces);
	
	for(int i = 0; i < top.quads.cols(); i++) {
		triangles[2*i].v0 = top.quads(0,i);
		triangles[2*i].v1 = top.quads(1,i);
		triangles[2*i].v2 = top.quads(2,i);
		
		triangles[2*i+1].v0 = top.quads(0,i);
		triangles[2*i+1].v1 = top.quads(2,i);
		triangles[2*i+1].v2 = top.quads(3,i);
	}
	
	rtcCommitGeometry(mesh);
	unsigned int geomID = rtcAttachGeometry(scene_i,mesh);
	rtcReleaseGeometry(mesh);
	return geomID;
}

unsigned int RayTracer::addSubdivModel (RTCScene scene_i, const Matrix3X& V, const MeshTopology& top) {
	RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SUBDIVISION); //SUBDIVISION
	rtcSetGeometryTessellationRate(mesh,8.0f);
	
	//face_colors = (Vec3fa*) alignedMalloc(12*sizeof(Vec3fa));
	//vertex_colors = (Vec3fa*) alignedMalloc(8*sizeof(Vec3fa));
	int nfaces = top.quads.cols();
	Vertex* vertices = (Vertex*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,sizeof(Vertex),V.cols()); 
	for(int i = 0; i < V.cols(); i++) {
		vertices[i].x = V(0,i);
		vertices[i].y = V(1,i);
		vertices[i].z = V(2,i);
	}
	
	unsigned int* edgeCount = (unsigned*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_FACE,0,RTC_FORMAT_UINT,sizeof(unsigned int),nfaces); 
	for(int i = 0; i < nfaces; i++) {
		edgeCount[i] = 4;
	}
	
	unsigned int* indices = (unsigned int*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT,sizeof(unsigned int),4*nfaces);
	
	for(int i = 0; i < top.quads.cols(); i++) {
		for(int j = 0; j < 4; j++) {
			indices[4*i + j] = top.quads(j,i);
		}
	}
	
	rtcCommitGeometry(mesh);
	unsigned int geomID = rtcAttachGeometry(scene_i,mesh);
	rtcReleaseGeometry(mesh);
	return geomID;
}

void RayTracer::device_init (char* cfg, const Matrix3X& V, const MeshTopology& top) {
	//std::cout << "Embree3 Current Memory " << bytes << std::endl;
	
	
	assert(device == nullptr);
	device = rtcNewDevice(cfg);
	error_handler(nullptr,rtcGetDeviceError(device));
	
	rtcSetDeviceErrorFunction(device,error_handler,nullptr);
	//rtcSetDeviceMemoryMonitorFunction(device,memory_handler,(void*) &bytes);


	assert(rtcGetDeviceProperty(device,RTC_DEVICE_PROPERTY_SUBDIVISION_GEOMETRY_SUPPORTED) != 0);

	/* create scene */
	assert(scene == nullptr);
	scene = rtcNewScene(device);

	/* add cube */
	addSubdivModel(scene, V, top);

	/* commit changes to scene */
	rtcCommitScene (scene);
}

void RayTracer::device_cleanup ()
{
  if(scene != nullptr)
	rtcReleaseScene(scene); 
  scene = nullptr;
  
  
  if(device != nullptr) {
	  //rtcSetDeviceMemoryMonitorFunction(device,nullptr,(void*) nullptr);
	  rtcSetDeviceErrorFunction(device,nullptr,nullptr);
	  rtcReleaseDevice(device); 
      device = nullptr;
  }
}

RayTracer::RayTracer() {
	device = nullptr;
	scene = nullptr;
	bytes = 0;
}

float RayTracer::shoot(Vector3 p, Vector3 n, Vector3 dn, Vector2* u, int* f) {
	const float range = RAY_MAX+.019;
	assert(n.norm() > 0);
	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	
	RTCRayHit rayHit;
	RTCRay ray;
	ray.org_x = p[0];
	ray.org_y = p[1];
	ray.org_z = p[2];
	ray.tnear = 0.0;
	
	ray.dir_x = n[0];
	ray.dir_y = n[1];
	ray.dir_z = n[2];
	ray.time = 0.0;
	
	ray.tfar = range;
	ray.mask = 0;
	ray.id = 0;
	ray.flags = 0;

	rayHit.ray = ray;

	rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
	rayHit.hit.primID = RTC_INVALID_GEOMETRY_ID;
	for(int i = 0; i < RTC_MAX_INSTANCE_LEVEL_COUNT; i++) {
		rayHit.hit.instID[i] = RTC_INVALID_GEOMETRY_ID;
	}
	
	rtcIntersect1(scene,&context,&rayHit);
	float dist = 2*range;
	//std::cout << "Ray: " << p.transpose() << " -> " << (-n).transpose() << std::endl;
	if(rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
		if(rayHit.hit.Ng_x*dn[0] + rayHit.hit.Ng_y*dn[1] + rayHit.hit.Ng_z*dn[2] > 0) {
			dist = rayHit.ray.tfar;
			*f = rayHit.hit.primID;
			(*u)(0,0) = rayHit.hit.u;
			(*u)(1,0) = rayHit.hit.v;
			if((*u)(0,0) >= 0.0 && (*u)(0,0) <= 1.0 && (*u)(1,0) >= 0.0 && (*u)(1,0) <= 1.0) {
			} else {
		 		std::cerr << (*f) << ":" <<  u->transpose() << std::endl;
				(*u)(0,0) = fmax((*u)(0,0),0);
				(*u)(1,0) = fmax((*u)(1,0),0);
				(*u)(0,0) = fmin((*u)(0,0),1);
                (*u)(1,0) = fmin((*u)(1,0),1);
			}
			assert((*u)(0,0) >= 0.0 && (*u)(0,0) <= 1.0);
			assert((*u)(1,0) >= 0.0 && (*u)(1,0) <= 1.0);
		}
		//std::cout << "primID: " << rayHit.hit.primID << " -> " << rayHit.ray.tfar << std::endl;
	} else {
		//std::cout << "ray missed" << std::endl;
	}
	return dist;
}

Vector3 RayTracer::shootNormals(Vector3 p, Vector3 n, Vector2* uout, int* fout) {
	const float range = RAY_MAX; //.21  .081
	Vector2 u1, u2;
	int f1, f2;
	
	float d1 = shoot(p,n, n, &u1, &f1);
	float d2 = -shoot(p,-n, n, &u2, &f2);
	float d = fabs(d1) < fabs(d2) ? d1 : d2;
	Vector2 u = fabs(d1) < fabs(d2) ? u1 : u2;
	int f = fabs(d1) < fabs(d2) ? f1 : f2;
	
	if(fabs(d) > 1.5*range) {
		d = 0;
		*fout = -1;
	} else {
		//std::cout << f << " : " << u.transpose() << std::endl;
		*fout = f;
		*uout = u;
	}
	Vector3 pnew = p + d*n;
		
	return pnew;
}
