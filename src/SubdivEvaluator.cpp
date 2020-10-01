#include "SubdivEvaluator.h"


SubdivEvaluator::SubdivEvaluator(SubdivEvaluator const& that) {
	*this = that;
}

SubdivEvaluator& SubdivEvaluator::operator=(SubdivEvaluator const& that) {
	this->nVertices = that.nVertices;
	this->nRefinerVertices = that.nRefinerVertices;
	this->nLocalPoints = that.nLocalPoints;
	this->compressed_faces = that.compressed_faces;
	this->num_faces = that.num_faces;
	this->vertsperface = that.vertsperface;
	//this->init();
	//this->patchTable = new OpenSubdiv::Far::PatchTable(*that.patchTable);
	//this->refiner = new OpenSubdiv::Far::TopologyRefiner(*that.refiner)
	//this->evaluation_verts_buffer = that.evaluation_verts_buffer;
	return *this;
}

SubdivEvaluator::~SubdivEvaluator() {
	delete patchMap;
	delete patchTable;
	delete refiner;
	// xxawf delete refiner2
}

SubdivEvaluator::SubdivEvaluator(MeshTopology const& mesh)
{
	nVertices = mesh.num_vertices;

	num_faces = mesh.num_faces();

	//Fill the topology of the mesh

	vertsperface.resize((int)num_faces);
	//vertsperface.setConstant((int) mesh.quads.rows());
	for (int i = 0; i < (int)num_faces; i++) {
		vertsperface[i] = (int)mesh.nverts_per_face[i];
	}
	compressed_faces = mesh.compressed_faces;
	//init();
}

const OpenSubdiv::Far::StencilTable*  SubdivEvaluator::getLocalStencilTable() {
    return cvstencils;
}
OpenSubdiv::Far::StencilTable const * SubdivEvaluator::GetLocalPointStencilTable() {
    return patchTable->GetLocalPointStencilTable();
}

Matrix3X SubdivEvaluator::localPoints(const Matrix3X& base) {
    //std::cout << "Local Point Computation: " << nRefinerVertices << ", " << nLocalPoints 
    //    << ", " << patchTable->GetLocalPointStencilTable()->GetNumStencils() << std::endl;
    evaluation_verts_buffer.resize(nRefinerVertices + nLocalPoints);

    for(int i = 0; i < nRefinerVertices; i++)
        evaluation_verts_buffer[i].point = base.col(i);
    /*
    int nRefinedLevels = refiner->GetNumLevels();

    // Interpolate vertex primvar data : they are the control vertices
    // of the limit patches (see far_tutorial_0 for details)
    OSD_Vertex * src = &evaluation_verts_buffer[0];
    for (int level = 1; level < nRefinedLevels; ++level) {
        OSD_Vertex * dst = src + refiner->GetLevel(level-1).GetNumVertices();
        Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
        src = dst;
    }
    */
	// Evaluate local points from interpolated vertex primvars.
	int numStencils = patchTable->GetLocalPointStencilTable()->GetNumStencils();
    assert(nRefinerVertices + nLocalPoints == numStencils);

	patchTable->ComputeLocalPointValues(&evaluation_verts_buffer[0], &evaluation_verts_buffer[nRefinerVertices]);
    Matrix3X output(3,nRefinerVertices + nLocalPoints);
    for(int i = 0; i < nRefinerVertices + nLocalPoints; i++) {
        output.col(i) = evaluation_verts_buffer[i].point;
    }
    return output;
}

void SubdivEvaluator::init(Matrix3X const& vert_coords, Matrix3X* samples) {
	Far::TopologyDescriptor desc;
	desc.numVertices = (int)nVertices;
	desc.numFaces = (int)num_faces;

	desc.numVertsPerFace = vertsperface.data();
	//desc.vertIndicesPerFace = mesh.quads.data();
	desc.vertIndicesPerFace = compressed_faces.data();

	//Instantiate a FarTopologyRefiner from the descriptor.
	Sdc::SchemeType type = OpenSubdiv::Sdc::SCHEME_CATMARK;
	// Adpative refinement is only supported for CATMARK
	// Scheme LOOP is only supported if the mesh is purely composed of triangles

	Sdc::Options options;
	options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_NONE);
	typedef Far::TopologyRefinerFactory<Far::TopologyDescriptor> Refinery;
	//OpenSubdiv::Far::TopologyRefiner *refiner = Refinery::Create(desc, Refinery::Options(type, options));
	this->refiner = Refinery::Create(desc, Refinery::Options(type, options));

	const int maxIsolation = 0; //Don't change it!
	refiner->RefineAdaptive(Far::TopologyRefiner::AdaptiveOptions(maxIsolation));
	Far::TopologyLevel toplevel = refiner->GetLevel(0);
	// Generate a set of Far::PatchTable that we will use to evaluate the surface limit
	Far::PatchTableFactory::Options patchOptions;
	patchOptions.useInfSharpPatch = false;
	patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_BSPLINE_BASIS;
    //patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
	//patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_LEGACY_GREGORY;

	patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
	//std::cout << patchTable->GetNumPtexFaces() << std::endl;
	//std::cout << refiner->GetNumVerticesTotal() << std::endl;

	// Compute the total number of points we need to evaluate patchtable.
	// we use local points around extraordinary features.
	nRefinerVertices = refiner->GetNumVerticesTotal();
	nLocalPoints = patchTable->GetNumLocalPoints();


	/*OpenSubdiv::Far::StencilTableFactory::Options options;
	options.S
	options.generateIntermediateLevels = false;
	options.generateControlVerts = true;
	options.generateOffsets = true;*/

	//cvstencils = patchTable->GetLocalPointStencilTable();
	OpenSubdiv::Far::StencilTableFactory::Options opt;
	opt.generateIntermediateLevels = false;
	opt.generateControlVerts = true;
	opt.generateOffsets = true;
	cvstencils = OpenSubdiv::Far::StencilTableFactory::Create(*refiner, opt);
	if (OpenSubdiv::Far::StencilTable const *localPointStencilTable =
		patchTable->GetLocalPointStencilTable()) {
		OpenSubdiv::Far::StencilTable const *table =
			OpenSubdiv::Far::StencilTableFactory::AppendLocalPointStencilTable(
				*refiner, cvstencils, localPointStencilTable);
        delete cvstencils;
		cvstencils = table;
	}

    for(int i = 0; i < cvstencils->GetNumStencils(); i++) {
        auto st = cvstencils->GetStencil(i);
        auto *ind = st.GetVertexIndices();
        float const *wei = st.GetWeights();
    }


	// Create a buffer to hold the position of the refined verts and
	// local points.
	evaluation_verts_buffer.resize(nRefinerVertices + nLocalPoints);


	// xxaqwf delete refiner here?

	// This refiner is to generate subdivided meshes
	// Instantiate a FarTopologyRefiner from the descriptor
	this->refiner2 = Refinery::Create(desc, Refinery::Options(type, options));

	// Uniformly refine the topolgy up to 'maxlevel'

	refiner2->RefineUniform(Far::TopologyRefiner::UniformOptions(maxlevel));
    ptexnum = patchTable->GetNumPtexFaces();

	const int ptexnum = patchTable->GetNumPtexFaces();
	const int patchnum = patchTable->GetNumPatchesTotal();

    std::vector<SurfacePoint> boundaries;
    Matrix3X out_S(3,3*ptexnum);
    

    for (int i = 0; i < ptexnum; i+=1) {
        SurfacePoint sp;
        sp.face = i;
        sp.u[0] = 0;
        sp.u[1] = 0;
        boundaries.push_back(sp);
        sp.u[0] = 1;
        sp.u[1] = 0;
        boundaries.push_back(sp);
        sp.u[0] = 0;
        sp.u[1] = 1;
        boundaries.push_back(sp);
    }
    evaluateSubdivSurface(vert_coords,
	boundaries,	&out_S);

    const Scalar sampling = SAMPLING; //0.01 0.005
    const int minpts = 1;
    //int perFace = 33;
	for (int i = 0; i < ptexnum; i+=1) { //ptexnum
        Scalar lu = (out_S.col(3*i)-out_S.col(3*i+1)).norm();
        Scalar lv = (out_S.col(3*i)-out_S.col(3*i+2)).norm();
        int su = (int)(lu / sampling);
        int sv = (int)(lv / sampling);
        
        su = su <= minpts ? minpts : su;
        sv = sv <= minpts ? minpts : sv;
		for (int u = 0; u <  su; u++) {
			for (int v = 0; v < sv; v++) {
				SurfacePoint sp;
				sp.face = i; //rand() % ptexnum;
				sp.u[0] = .5*(1.0/su) + (1.0/su)*u;
				sp.u[1] = .5*(1.0/sv) + (1.0/sv)*v;
				controlVertexPtex.push_back(sp);
			}
		}
		//SurfacePoint sp;
		//sp.face = i;
		//sp.u[0] = .5;
		//sp.u[1] = .5;
		//controlVertexPtex.push_back(sp);
	}
    if(samples != NULL) samples->resize(3,controlVertexPtex.size());
    evaluateSubdivSurface(vert_coords,
	controlVertexPtex,	samples);

	this->patchMap = new Far::PatchMap(*patchTable);
	Far::PatchMap patchmap(*patchTable);
	const int numParticles = ptexnum; //* 5 * 5;
	STParticles particles(*refiner, patchTable, numParticles);
	adjacency = particles.GetAdjacency();
	auto pos = particles.GetPositions();
	/*
	for (int i = 0; i < pos.size(); i++) {
	OpenSubdiv::Far::PatchTable::PatchHandle const *handle =
	patchmap.FindPatch(pos[i].ptexIndex, pos[i].s, pos[i].t);
	if (handle) {
	float pWeights[MAX_NUM_W];
	auto coord = OpenSubdiv::Osd::PatchCoord(*handle, pos[i].s, pos[i].t);
	patchTable->EvaluateBasis(handle, pos[i].s, pos[i].t, pWeights);
	}
	}
	*/
	/*
	for (int i = 0; i < numParticles; i++) {
	controlVertexPtex[i].face = pos[i].ptexIndex;
	controlVertexPtex[i].u[0] = pos[i].s;
	controlVertexPtex[i].u[1] = pos[i].t;
	}
	*/

}

void SubdivEvaluator::generate_refined_mesh(Matrix3X const& vert_coords, int levels, MeshTopology* mesh_out, Matrix3X* verts_out)
{
	if (levels > maxlevel) {
		std::cerr << "SubdivEvaluator::generate_refined_mesh: level too high\n";
		levels = maxlevel;
	}

	// Allocate a buffer for vertex primvar data. The buffer length is set to
	// be the sum of all children vertices up to the highest level of refinement.
	std::vector<OSD_Vertex> vbuffer(refiner2->GetNumVerticesTotal());
	OSD_Vertex * verts = &vbuffer[0];


	// Initialize coarse mesh positions
	int nCoarseVerts = (int)vert_coords.cols();
	for (int i = 0; i<nCoarseVerts; ++i)
		verts[i].point = vert_coords.col(i);

	// Interpolate vertex primvar data
	Far::PrimvarRefiner primvarRefiner(*refiner2);

	OSD_Vertex * src = verts;
	for (int level = 1; level <= levels; ++level) {
		OSD_Vertex * dst = src + refiner2->GetLevel(level - 1).GetNumVertices();
		primvarRefiner.Interpolate(level, src, dst);
		src = dst;
	}

	Far::TopologyLevel const & refLastLevel = refiner2->GetLevel(levels);

	int nverts = refLastLevel.GetNumVertices();
	int nfaces = refLastLevel.GetNumFaces();

	// Print vertex positions
	//int firstOfLastVerts = refiner2->GetNumVerticesTotal() - nverts;
	// src -= nverts;

	verts_out->resize(3, nverts);
	for (int vert = 0; vert < nverts; ++vert)
		verts_out->col(vert) = src[vert].point;

	// Print faces
	mesh_out->num_vertices = nverts;
	mesh_out->quads.resize(5, nfaces);
	mesh_out->quads.setConstant(-1);
	for (int face = 0; face < nfaces; ++face) {

		Far::ConstIndexArray fverts = refLastLevel.GetFaceVertices(face);

		// all refined Catmark faces should be quads
		//assert(fverts.size() == 4);

		for (int vert = 0; vert < fverts.size(); ++vert)
			mesh_out->quads(vert, face) = fverts[vert];
	}
	mesh_out->update_adjacencies();
}

void SubdivEvaluator::evaluatePtex(const Matrix3X& vert_coords,
	Matrix3X* out_S,
	triplets_t* out_dSdX,
	triplets_t* out_dSudX,
	triplets_t* out_dSvdX,
	Matrix3X* out_Su,
	Matrix3X* out_Sv,
    bool printDebug
) const {
	evaluateSubdivSurface(vert_coords, controlVertexPtex, out_S,
		out_dSdX, out_dSudX, out_dSvdX, out_Su, out_Sv, NULL, NULL, NULL, NULL, NULL, NULL,printDebug);
}
/*
void SubdivEvaluator::move(const MeshTopology& mesh, const Matrix3X& vert_coords, const Matrix2X& du, Matrix2X* dir) {
	Matrix3X out_S(3, controlVertexPtex.size());
	int nhops = 0;
	int loopers = 0;
	int totalhops = 0;
	evaluateSubdivSurface(vert_coords, controlVertexPtex, &out_S);
	for (int i = 0; i < controlVertexPtex.size(); i++) { //controlVertexPtex.size()
		int face = controlVertexPtex[i].face;
		int newFace;
		Vector2 new_u;
		Vector2 idir; idir.setConstant(-100);
		int nhops = increment_u_crossing_edges(mesh, vert_coords, face, controlVertexPtex[i].u, du.col(i), &newFace, &new_u, &idir);
		dir->col(i) = idir;
		if (nhops < 0)
			++loopers;
		totalhops += std::abs(nhops);
		controlVertexPtex[i].face = newFace;
		controlVertexPtex[i].u = new_u;
	}
	if (loopers > 0)
		std::cerr << "[" << totalhops / Scalar(controlVertexPtex.size()) << " hops, " << loopers << " points looped]";
	else if (totalhops > 0)
		std::cerr << "[" << totalhops << "/" << Scalar(controlVertexPtex.size()) << " hops]";

}
*/
int SubdivEvaluator::increment_u_crossing_edges(const MeshTopology& mesh, Matrix3X const& X, int face, const Vector2& u, const Vector2& du, int* new_face_out, Vector2* new_u_out, Vector2* dir, Scalar* distance,std::vector<Vector3>* path)
{
	const int MAX_HOPS = 1000;
	const int MAX_FIX_HOPS = 3; //100

    assert((u)[0] >= 0 && (u)[1] <= 1.0 && (u)[1] >= 0 && (u)[1] <= 1.0);

	Scalar u1_old = u[0];
	Scalar u2_old = u[1];
	Scalar du1 = du[0];
	Scalar du2 = du[1];
	Scalar u1_new = u1_old + du1;
	Scalar u2_new = u2_old + du2;
    Vector3 pcurrent;
    /*
    {
        std::vector<SurfacePoint> pts;
	    pts.push_back({ face,{ u[0], u[1] } });
	    Matrix3X S(3, 1);
	    this->evaluateSubdivSurface(X, pts, &S);
        pcurrent = S.col(0);
    }
    */
    if(distance != NULL)
        *distance = 0;

	dir->setZero();

	for (int count = 0; ; ++count) {
		bool crossing = (u1_new < 0.f) || (u1_new > 1.f) || (u2_new < 0.f) || (u2_new > 1.f);

		if (!crossing) {
			*new_face_out = face;
			*new_u_out << u1_new, u2_new;
			(*dir)[0] = du1;
			(*dir)[1] = du2;
			dir->normalize();
            
            if(0){
                std::vector<SurfacePoint> pts;
	            pts.push_back({ face,{ u1_new, u2_new } });
	            Matrix3X S(3, 1);
	            this->evaluateSubdivSurface(X, pts, &S);
                *distance += (S.col(0)-pcurrent).norm();
            }
            //std::cout << (*new_face_out) << " -> " << new_u_out->transpose()  << std::endl;

            assert((*new_u_out)[0] >= 0 && (*new_u_out)[1] <= 1.0 && (*new_u_out)[1] >= 0 && (*new_u_out)[1] <= 1.0);
			return count;
		}

		//Find the new face	and the coordinates of the crossing point within the old face and the new face
		int face_new;

		bool face_found = false;

		Scalar dif, aux, u1_cross, u2_cross;

		if (u1_new < 0.f)
		{
			dif = u1_old;
			const Scalar u2t = u2_old - du2*dif / du1;
			if ((u2t >= 0.f) && (u2t <= 1.f))
			{
				//face_new = adjacency[face].adjface(3); aux = u2t; face_found = true;
				face_new = mesh.face_adj(3, face); aux = u2t; face_found = true;
				u1_cross = 0.f; u2_cross = u2t;
			}
		}
		if ((u1_new > 1.f) && (!face_found))
		{
			dif = 1.f - u1_old;
			const Scalar u2t = u2_old + du2*dif / du1;
			if ((u2t >= 0.f) && (u2t <= 1.f))
			{
				//face_new = adjacency[face].adjface(1); aux = 1.f - u2t; face_found = true;
				face_new = mesh.face_adj(1, face); aux = 1.f - u2t; face_found = true;
				u1_cross = 1.f; u2_cross = u2t;
			}
		}
		if ((u2_new < 0.f) && (!face_found))
		{
			dif = u2_old;
			const Scalar u1t = u1_old - du1*dif / du2;
			if ((u1t >= 0.f) && (u1t <= 1.f))
			{
				//face_new = adjacency[face].adjface(0); aux = 1.f - u1t; face_found = true;
				face_new = mesh.face_adj(0, face); aux = 1.f - u1t; face_found = true;
				u1_cross = u1t; u2_cross = 0.f;
			}
		}
		if ((u2_new > 1.f) && (!face_found))
		{
			dif = 1.f - u2_old;
			const Scalar u1t = u1_old + du1*dif / du2;
			if ((u1t >= 0.f) && (u1t <= 1.f))
			{
				//face_new = adjacency[face].adjface(2); aux = u1t; face_found = true;
				face_new = mesh.face_adj(2, face); aux = u1t; face_found = true;
				u1_cross = u1t; u2_cross = 1.f;
			}
		}
		assert(face_found);
		if (!face_found) {
			std::cout << "error" << std::endl;
		}

		// Find the coordinates of the crossing point as part of the new face, and update u_old (as that will be new u in next iter).
		unsigned int conf = 1000000;
		for (unsigned int f = 0; f < 4; f++)
			if (mesh.face_adj(f, face_new) == face) { conf = f; }
		//if (adjacency[face_new].adjface(face) == face) { conf = f; }

		switch (conf)
		{
		case 0: u1_old = aux; u2_old = 0.f; break;
		case 1: u1_old = 1.f; u2_old = aux; break;
		case 2:	u1_old = 1.f - aux; u2_old = 1.f; break;
		case 3:	u1_old = 0.f; u2_old = 1.f - aux; break;
        default: std::cout << "ERROR! " << std::endl;
		}

		// Evaluate the subdivision surface at the edge (with respect to the original face)
		std::vector<SurfacePoint> pts;
		pts.push_back({ face,{ u1_cross, u2_cross } });
		pts.push_back({ face_new,{ u1_old, u2_old } });
		Matrix3X S(3, 2);
		Matrix3X Su(3, 2);
		Matrix3X Sv(3, 2);
		//this->evaluateSubdivSurface(X, pts, &S, NULL, NULL, NULL, &Su, &Sv);
#if 1
		evalPointsPartial(X, pts,&S,&Su, &Sv);
#else
		evaluateSubdivSurface(X, pts, &S, NULL, NULL, NULL, &Su, &Sv);
#endif

        if(path != NULL) {
            path->push_back(S.col(0));
            path->push_back(S.col(1));
        }
    
		Matrix32 J_Sa;
		J_Sa.col(0) = Su.col(0);
		J_Sa.col(1) = Sv.col(0);
        //Matrix22 F1 = J_SA.transpose() * J_SA;

		Matrix32 J_Sb;
		J_Sb.col(0) = Su.col(1);
		J_Sb.col(1) = Sv.col(1);

        // *distance += (pcurrent-S.col Matrix3X localPoints(const Matrix3X& base);(0)).norm();
        pcurrent = S.col(0);

		//Compute the new u increments
		Vector2 du_remaining;
		du_remaining << u1_new - u1_cross, u2_new - u2_cross;
		Vector3 prod = J_Sa*du_remaining;
		Matrix22 AtA = J_Sb.transpose()*J_Sb;
		Vector2 AtB = J_Sb.transpose()*prod;

		//Vector2 du_new = AtA.ldlt().solve(AtB);
		Vector2  u_incr = AtA.inverse()*AtB;

		du1 = u_incr[0];
		du2 = u_incr[1];
		if (count == MAX_FIX_HOPS) {
			auto dmax = std::max(du1, du2);
			//std::cout << "nudge a lot" << std::endl;
			Scalar scale = Scalar(0.5 / dmax);
			*new_face_out = face;
			//(*new_u_out) << (u1_old + du1 * scale), (u2_old + du2 * scale);
			(*new_u_out) << 0.5, 0.5;
			//(*new_u_out)[0] = .5;//u1_old;
			//(*new_u_out)[1] = .5;//u2_old;
            //std::cerr << "Problem!!! Many jumps between the mesh faces for the update of one correspondence. I remove the remaining u_increment!\n";
			return -count;
		}
		else if (count == MAX_HOPS) {
			//std::cerr << "Problem!!! Many jumps between the mesh faces for the update of one correspondence. I remove the remaining u_increment!\n";
			auto dmax = std::max(du1, du2);
			Scalar scale = Scalar(0.5 / dmax);
			*new_face_out = face;
			// *new_u_out << u1_old + du1 * scale, u2_old + du2 * scale;
			*new_u_out << 0.5, 0.5;

			assert((*new_u_out)[0] >= 0 && (*new_u_out)[1] <= 1.0 && (*new_u_out)[1] >= 0 && (*new_u_out)[1] <= 1.0);
			return -count;
		}
		else {

			u1_new = u1_old + du1;
			u2_new = u2_old + du2;
			face = face_new;
		}
	}
}

int SubdivEvaluator::increment_u_crossing_edges_simple(const MeshTopology& mesh, Matrix3X const& X, int face, const Vector2& u, const Vector2& du, int* new_face_out, Vector2* new_u_out, Vector2* dir, Scalar* distance,std::vector<Vector3>* path)
{
    const int MAX_HOPS = 1000;
    const int MAX_FIX_HOPS = 3; //100

    assert((u)[0] >= 0 && (u)[1] <= 1.0 && (u)[1] >= 0 && (u)[1] <= 1.0);

    Scalar u1_old = u[0];
    Scalar u2_old = u[1];
    Scalar du1 = du[0];
    Scalar du2 = du[1];
    Scalar u1_new = u1_old + du1;
    Scalar u2_new = u2_old + du2;
    Vector3 pcurrent;
    /*
    {
        std::vector<SurfacePoint> pts;
        pts.push_back({ face,{ u[0], u[1] } });
        Matrix3X S(3, 1);
        this->evaluateSubdivSurface(X, pts, &S);
        pcurrent = S.col(0);
    }
    */
    if(distance != NULL)
        *distance = 0;

    dir->setZero();

    for (int count = 0; ; ++count) {
        bool crossing = (u1_new < 0.f) || (u1_new > 1.f) || (u2_new < 0.f) || (u2_new > 1.f);

        if (!crossing) {
            *new_face_out = face;
            *new_u_out << u1_new, u2_new;
            (*dir)[0] = du1;
            (*dir)[1] = du2;
            dir->normalize();
            
            if(0){
                std::vector<SurfacePoint> pts;
                pts.push_back({ face,{ u1_new, u2_new } });
                Matrix3X S(3, 1);
                this->evaluateSubdivSurface(X, pts, &S);
                *distance += (S.col(0)-pcurrent).norm();
            }
            //std::cout << (*new_face_out) << " -> " << new_u_out->transpose()  << std::endl;

            assert((*new_u_out)[0] >= 0 && (*new_u_out)[1] <= 1.0 && (*new_u_out)[1] >= 0 && (*new_u_out)[1] <= 1.0);
            return count;
        }

        //Find the new face and the coordinates of the crossing point within the old face and the new face
        int face_new;

        bool face_found = false;

        Scalar dif, aux, u1_cross, u2_cross;

        if (u1_new < 0.f)
        {
            dif = u1_old;
            const Scalar u2t = u2_old - du2*dif / du1;
            if ((u2t >= 0.f) && (u2t <= 1.f))
            {
                //face_new = adjacency[face].adjface(3); aux = u2t; face_found = true;
                face_new = mesh.face_adj(3, face); aux = u2t; face_found = true;
                u1_cross = 0.f; u2_cross = u2t;
            }
        }
        if ((u1_new > 1.f) && (!face_found))
        {
            dif = 1.f - u1_old;
            const Scalar u2t = u2_old + du2*dif / du1;
            if ((u2t >= 0.f) && (u2t <= 1.f))
            {
                //face_new = adjacency[face].adjface(1); aux = 1.f - u2t; face_found = true;
                face_new = mesh.face_adj(1, face); aux = 1.f - u2t; face_found = true;
                u1_cross = 1.f; u2_cross = u2t;
            }
        }
        if ((u2_new < 0.f) && (!face_found))
        {
            dif = u2_old;
            const Scalar u1t = u1_old - du1*dif / du2;
            if ((u1t >= 0.f) && (u1t <= 1.f))
            {
                //face_new = adjacency[face].adjface(0); aux = 1.f - u1t; face_found = true;
                face_new = mesh.face_adj(0, face); aux = 1.f - u1t; face_found = true;
                u1_cross = u1t; u2_cross = 0.f;
            }
        }
        if ((u2_new > 1.f) && (!face_found))
        {
            dif = 1.f - u2_old;
            const Scalar u1t = u1_old + du1*dif / du2;
            if ((u1t >= 0.f) && (u1t <= 1.f))
            {
                //face_new = adjacency[face].adjface(2); aux = u1t; face_found = true;
                face_new = mesh.face_adj(2, face); aux = u1t; face_found = true;
                u1_cross = u1t; u2_cross = 1.f;
            }
        }
        assert(face_found);
        if (!face_found) {
            std::cout << "error" << std::endl;
        }

        // Find the coordinates of the crossing point as part of the new face, and update u_old (as that will be new u in next iter).
        unsigned int conf = 1000000;
        for (unsigned int f = 0; f < 4; f++)
            if (mesh.face_adj(f, face_new) == face) { conf = f; }
        //if (adjacency[face_new].adjface(face) == face) { conf = f; }

        switch (conf)
        {
        case 0: u1_old = aux; u2_old = 0.f; break;
        case 1: u1_old = 1.f; u2_old = aux; break;
        case 2: u1_old = 1.f - aux; u2_old = 1.f; break;
        case 3: u1_old = 0.f; u2_old = 1.f - aux; break;
        default: std::cout << "ERROR! " << std::endl;
        }

        // Evaluate the subdivision surface at the edge (with respect to the original face)
        std::vector<SurfacePoint> pts;
        //pts.push_back({ face,{ u1_cross, u2_cross } });
        //pts.push_back({ face_new,{ u1_old, u2_old } });

        du1 = 0.5;
        du2 = 0.5;
        if (count == MAX_FIX_HOPS) {
            auto dmax = std::max(du1, du2);
            //std::cout << "nudge a lot" << std::endl;
            Scalar scale = Scalar(0.5 / dmax);
            *new_face_out = face;
            //(*new_u_out) << (u1_old + du1 * scale), (u2_old + du2 * scale);
            (*new_u_out) << 0.5, 0.5;
            //(*new_u_out)[0] = .5;//u1_old;
            //(*new_u_out)[1] = .5;//u2_old;
            //std::cerr << "Problem!!! Many jumps between the mesh faces for the update of one correspondence. I remove the remaining u_increment!\n";
            return -count;
        }
        else if (count == MAX_HOPS) {
            //std::cerr << "Problem!!! Many jumps between the mesh faces for the update of one correspondence. I remove the remaining u_increment!\n";
            auto dmax = std::max(du1, du2);
            Scalar scale = Scalar(0.5 / dmax);
            *new_face_out = face;
            // *new_u_out << u1_old + du1 * scale, u2_old + du2 * scale;
            *new_u_out << 0.5, 0.5;

            assert((*new_u_out)[0] >= 0 && (*new_u_out)[1] <= 1.0 && (*new_u_out)[1] >= 0 && (*new_u_out)[1] <= 1.0);
            return -count;
        }
        else {

            u1_new = u1_old + du1;
            u2_new = u2_old + du2;
            face = face_new;
        }
    }
}



int SubdivEvaluator::increment_accurate(const MeshTopology& mesh, Matrix3X const& X, int face, const Vector2& u, const Vector2& du, int* new_face_out, Vector2* new_u_out, Vector2* dir, Scalar* distance,std::vector<Vector3>* path)
{
	const int MAX_HOPS = 1000;
	const int MAX_FIX_HOPS = 3; //100

    assert((u)[0] >= 0 && (u)[1] <= 1.0 && (u)[1] >= 0 && (u)[1] <= 1.0);

	Scalar u1_old = u[0];
	Scalar u2_old = u[1];
    Scalar du1_total = du[0];
	Scalar du2_total = du[1];
	Scalar du1 = du[0];
	Scalar du2 = du[1];
	Scalar u1_new = u1_old + du1;
	Scalar u2_new = u2_old + du2;
    Vector3 pcurrent;
    /*
    {
        std::vector<SurfacePoint> pts;
	    pts.push_back({ face,{ u[0], u[1] } });
	    Matrix3X S(3, 1);
	    this->evaluateSubdivSurface(X, pts, &S);
        pcurrent = S.col(0);
    }
    */
    if(distance != NULL)
        *distance = 0;

	dir->setZero();

    while(du1_total != 0 && du2_total != 0) {

        du1 = fmax(fmin(du1_total, .01f),-.01);
        du2 = fmax(fmin(du2_total, .01f),-.01);
        du1_total -= du1;
        du2_total -= du2;

	    for (int count = 0; ; ++count) {
		    bool crossing = (u1_new < 0.f) || (u1_new > 1.f) || (u2_new < 0.f) || (u2_new > 1.f);

		    if (!crossing) {
			    *new_face_out = face;
			    *new_u_out << u1_new, u2_new;
			    //(*dir)[0] = du1;
			    //(*dir)[1] = du2;
			    //dir->normalize();
                
                if(0){
                    std::vector<SurfacePoint> pts;
	                pts.push_back({ face,{ u1_new, u2_new } });
	                Matrix3X S(3, 1);
	                this->evaluateSubdivSurface(X, pts, &S);
                    path->push_back(S.col(0));
                    // *distance += (S.col(0)-pcurrent).norm();
                }
                //std::cout << (*new_face_out) << " -> " << new_u_out->transpose()  << std::endl;
                assert((*new_u_out)[0] >= 0 && (*new_u_out)[1] <= 1.0 && (*new_u_out)[1] >= 0 && (*new_u_out)[1] <= 1.0);
			    break;
		    }

		    //Find the new face	and the coordinates of the crossing point within the old face and the new face
		    int face_new;

		    bool face_found = false;

		    Scalar dif, aux, u1_cross, u2_cross;

		    if (u1_new < 0.f)
		    {
			    dif = u1_old;
			    const Scalar u2t = u2_old - du2*dif / du1;
			    if ((u2t >= 0.f) && (u2t <= 1.f))
			    {
				    //face_new = adjacency[face].adjface(3); aux = u2t; face_found = true;
				    face_new = mesh.face_adj(3, face); aux = u2t; face_found = true;
				    u1_cross = 0.f; u2_cross = u2t;
			    }
		    }
		    if ((u1_new > 1.f) && (!face_found))
		    {
			    dif = 1.f - u1_old;
			    const Scalar u2t = u2_old + du2*dif / du1;
			    if ((u2t >= 0.f) && (u2t <= 1.f))
			    {
				    //face_new = adjacency[face].adjface(1); aux = 1.f - u2t; face_found = true;
				    face_new = mesh.face_adj(1, face); aux = 1.f - u2t; face_found = true;
				    u1_cross = 1.f; u2_cross = u2t;
			    }
		    }
		    if ((u2_new < 0.f) && (!face_found))
		    {
			    dif = u2_old;
			    const Scalar u1t = u1_old - du1*dif / du2;
			    if ((u1t >= 0.f) && (u1t <= 1.f))
			    {
				    //face_new = adjacency[face].adjface(0); aux = 1.f - u1t; face_found = true;
				    face_new = mesh.face_adj(0, face); aux = 1.f - u1t; face_found = true;
				    u1_cross = u1t; u2_cross = 0.f;
			    }
		    }
		    if ((u2_new > 1.f) && (!face_found))
		    {
			    dif = 1.f - u2_old;
			    const Scalar u1t = u1_old + du1*dif / du2;
			    if ((u1t >= 0.f) && (u1t <= 1.f))
			    {
				    //face_new = adjacency[face].adjface(2); aux = u1t; face_found = true;
				    face_new = mesh.face_adj(2, face); aux = u1t; face_found = true;
				    u1_cross = u1t; u2_cross = 1.f;
			    }
		    }
		    assert(face_found);
		    if (!face_found) {
			    std::cout << "error" << std::endl;
				std::cout << u1_new << " , " << u2_new << std::endl;
				assert(face_found);
		    }

		    // Find the coordinates of the crossing point as part of the new face, and update u_old (as that will be new u in next iter).
		    unsigned int conf = 1000000;
		    for (unsigned int f = 0; f < 4; f++)
			    if (mesh.face_adj(f, face_new) == face) { conf = f; }
		    //if (adjacency[face_new].adjface(face) == face) { conf = f; }

		    switch (conf)
		    {
		    case 0: u1_old = aux; u2_old = 0.f; break;
		    case 1: u1_old = 1.f; u2_old = aux; break;
		    case 2:	u1_old = 1.f - aux; u2_old = 1.f; break;
		    case 3:	u1_old = 0.f; u2_old = 1.f - aux; break;
            default: std::cout << "ERROR! " << std::endl;
		    }

		    // Evaluate the subdivision surface at the edge (with respect to the original face)
		    std::vector<SurfacePoint> pts;
		    pts.push_back({ face,{ u1_cross, u2_cross } });
		    pts.push_back({ face_new,{ u1_old, u2_old } });
		    Matrix3X S(3, 2);
		    Matrix3X Su(3, 2);
		    Matrix3X Sv(3, 2);
		    this->evaluateSubdivSurface(X, pts, &S, NULL, NULL, NULL, &Su, &Sv);
            if(path != NULL) {
                path->push_back(S.col(0));
                path->push_back(S.col(1));
            }
        
		    Matrix32 J_Sa;
		    J_Sa.col(0) = Su.col(0);
		    J_Sa.col(1) = Sv.col(0);

		    Matrix32 J_Sb;
		    J_Sb.col(0) = Su.col(1);
		    J_Sb.col(1) = Sv.col(1);

            // *distance += (pcurrent-S.col(0)).norm();
            pcurrent = S.col(0);

		    //Compute the new u increments
		    Vector2 du_remaining;
		    du_remaining << u1_new - u1_cross, u2_new - u2_cross;
		    Vector3 prod = J_Sa*du_remaining;
		    Matrix22 AtA = J_Sb.transpose()*J_Sb;
		    Vector2 AtB = J_Sb.transpose()*prod;

		    //Vector2 du_new = AtA.ldlt().solve(AtB);
		    Vector2  u_incr = AtA.inverse()*AtB;

		    du1 = u_incr[0];
		    du2 = u_incr[1];
		    if (count == MAX_FIX_HOPS) {
			    auto dmax = std::max(du1, du2);
			    //std::cout << "nudge a lot" << std::endl;
			    Scalar scale = Scalar(0.5 / dmax);
			    *new_face_out = face;
			    //(*new_u_out) << (u1_old + du1 * scale), (u2_old + du2 * scale);
			    (*new_u_out) << 0.5, 0.5;
			    //(*new_u_out)[0] = .5;//u1_old;
			    //(*new_u_out)[1] = .5;//u2_old;
                //std::cerr << "Problem!!! Many jumps between the mesh faces for the update of one correspondence. I remove the remaining u_increment!\n";
			    break;
		    }
		    else if (count == MAX_HOPS) {
			    //std::cerr << "Problem!!! Many jumps between the mesh faces for the update of one correspondence. I remove the remaining u_increment!\n";
			    auto dmax = std::max(du1, du2);
			    Scalar scale = Scalar(0.5 / dmax);
			    *new_face_out = face;
			    // *new_u_out << u1_old + du1 * scale, u2_old + du2 * scale;
			    *new_u_out << 0.5, 0.5;

			    assert((*new_u_out)[0] >= 0 && (*new_u_out)[1] <= 1.0 && (*new_u_out)[1] >= 0 && (*new_u_out)[1] <= 1.0);
			    break;
		    }
		    else {

			    u1_new = u1_old + du1;
			    u2_new = u2_old + du2;
			    face = face_new;
		    }
	    }
    }
    return 0;
}

Matrix3X SubdivEvaluator::evalPointsCustom(const Matrix3X& base) {
    OpenSubdiv::Far::StencilTable const * localStencilTable = getLocalStencilTable();
    int localPointCount = localStencilTable->GetNumStencils();
    
    //std::cout << "LOCAL STENCIL SIZE: " << localPointCount << std::endl;
    int maxStencilBase = -1;
      
    Matrix3X localPts(3,localPointCount);
    for(int i = 0; i < localPointCount; i++) {
        auto st = localStencilTable->GetStencil(i);
        //std::cout <<  st.GetSize() << std::endl;
        auto *ind = st.GetVertexIndices();
        auto *wei = st.GetWeights();
        Vector3 localPt; localPt.setZero();
        //assert(st.GetSize() < 20);
        if(maxStencilBase < st.GetSize()) maxStencilBase = st.GetSize();
        for(int j = 0; j < st.GetSize(); j++) {
            assert(ind[j] < base.cols());
            //std::cout <<  wei[j] << " " << ind[j] << std::endl;
            localPt += wei[j]*base.col(ind[j]);
        }
        localPts.col(i) = localPt;
    }
    //std::cout << "Max Stencil Base: " << maxStencilBase << std::endl;
    //Matrix3X localPts = localPoints(base);

    Far::PatchMap patchmap(*patchTable);
    Matrix3X result(3,controlVertexPtex.size()); result.setZero();
     
    float pWeights[16];
    for(int k = 0; k < controlVertexPtex.size(); k++) {
        auto sp = controlVertexPtex[k];
        float u = sp.u[0];
        float v = sp.u[1];
        Far::PatchTable::PatchHandle const * handle = patchmap.FindPatch(sp.face, sp.u[0], sp.u[1]);
        assert(handle);
        Far::PatchParam param = patchTable->GetPatchParam(*handle);
        param.Normalize(u,v);

        
        patchTable->EvaluateBasis(*handle, sp.u[0], sp.u[1], pWeights);
        Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);

        float us[4];
        float vs[4];
        float u2 = u*u;
        float u3 = u2*u;
        us[0] = 1.0f/6.0f * (1.0f - 3.0f*(u - u2) - u3);
        us[1] = 1.0f/6.0f * (4.0f - 6.0f*u2 + 3.0f*u3);
        us[2] = 1.0f/6.0f * (1.0f + 3.0f*(u + u2 - u3));
        us[3] = 1.0f/6.0f * (u3);

        float v2 = v*v;
        float v3 = v2*v;
        vs[0] = 1.0f/6.0f * (1.0f - 3.0f*(v - v2) - v3);
        vs[1] = 1.0f/6.0f * (4.0f - 6.0f*v2 + 3.0f*v3);
        vs[2] = 1.0f/6.0f * (1.0f + 3.0f*(v + v2 - v3));
        vs[3] = 1.0f/6.0f * (v3);


        for(int i = 0; i < cvs.size(); i++) {
            float w = vs[i / 4] * us[i % 4];
            result.col(k) += /*pWeights[i]*/ w * localPts.col(cvs[i]);
        }
    }
    return result;
}

int SubdivEvaluator::minimizeSurfacePoint(const MeshTopology& mesh, SurfacePoint p, const Matrix3X& X, const Matrix3X& target, SurfacePoint& q, Matrix3X* pathCloud) {
	
    std::vector<Vector3> path;
	//q = p;
	SurfacePoint x = p;
	
	Scalar dOld, dNew;
	Scalar alpha = 1;
	const int cap = 6; 
	for(int iter = 0; iter < cap; iter++) {
		std::vector<SurfacePoint> inPoints; inPoints.push_back(x);
		Matrix3X S; S.resize(3,1);
		Matrix3X Su; Su.resize(3,1);
		Matrix3X Sv; Sv.resize(3,1);
		
		evalPointsPartial(X, inPoints, &S, &Su, &Sv);
		path.push_back(S.col(0));
		Matrix3X r = S - target;
		
		if(iter == 0) dOld = r.norm();
		if(iter == cap-1) {
			dNew = r.norm();
			break;
		}
		
		VectorX f; f.resize(3*r.cols());
		for(int i = 0; i < r.cols(); i++) {
			f[3*i+0] = r(0,i);
			f[3*i+1] = r(1,i);
			f[3*i+2] = r(2,i);
		}
		
		//std::cout << "F " << iter << ": " <<  f.stableNorm() << std::endl;
		
		Matrix2X JT; JT.resize(2,3*r.cols());
		for(int i = 0; i < r.cols(); i++) {
			JT(0,3*i+0) = Su(0,i);
			JT(0,3*i+1) = Su(1,i);
			JT(0,3*i+2) = Su(2,i);
			JT(1,3*i+0) = Sv(0,i);
			JT(1,3*i+1) = Sv(1,i);
			JT(1,3*i+2) = Sv(2,i);
		}
		auto J = JT.transpose();
		
		auto JTJ = JT*J;
		
		Vector2 g = JT*f;
		//Vector2 du = -alpha*g;
		
		Vector2  du = JTJ.inverse()*((-JT)*f);
		
		
        Scalar distance;

		Vector2 dir;
		int fout;
		Vector2 uout;
		int jumps = increment_u_crossing_edges(mesh, X, x.face, x.u, du, &fout, &uout, &dir, &distance, &path);
		x.face = fout;
		x.u = uout;
	}
	if(dOld < dNew) q = p;
	else q = x;
	/*
	if(pathCloud != NULL) {
		pathCloud->resize(3,path.size());
		for(int i = 0; i < path.size(); i++) {
			pathCloud->col(i) = path[i]; 
		}
	}
	*/ 
	return 0;
}


Matrix3X SubdivEvaluator::evalPointsCustom(const Matrix3X& base, const std::vector<SurfacePoint>& surfacePoints) {
    OpenSubdiv::Far::StencilTable const * localStencilTable = getLocalStencilTable();
    int localPointCount = localStencilTable->GetNumStencils();
    
    //std::cout << "LOCAL STENCIL SIZE: " << localPointCount << std::endl;
    int maxStencilBase = -1;
      
    Matrix3X localPts(3,localPointCount);
    for(int i = 0; i < localPointCount; i++) {
        auto st = localStencilTable->GetStencil(i);
        //std::cout <<  st.GetSize() << std::endl;
        auto *ind = st.GetVertexIndices();
        auto *wei = st.GetWeights();
        Vector3 localPt; localPt.setZero();
        //assert(st.GetSize() < 20);
        if(maxStencilBase < st.GetSize()) maxStencilBase = st.GetSize();
        for(int j = 0; j < st.GetSize(); j++) {
            assert(ind[j] < base.cols());
            //std::cout <<  wei[j] << " " << ind[j] << std::endl;
            localPt += wei[j]*base.col(ind[j]);
        }
        localPts.col(i) = localPt;
    }
    //std::cout << "Max Stencil Base: " << maxStencilBase << std::endl;
    //Matrix3X localPts = localPoints(base);

    Far::PatchMap patchmap(*patchTable);
    Matrix3X result(3,surfacePoints.size()); result.setZero();
     
    float pWeights[16];
    for(int k = 0; k < surfacePoints.size(); k++) {
        auto sp = surfacePoints[k];
        float u = sp.u[0];
        float v = sp.u[1];
	if(!((u>=0.0f) && (u<=1.0f) && (v>=0.0f) && (v<=1.0f) )) {
	    std::cerr << k << ": " << u  << "," << v << std::endl;
	}
        Far::PatchTable::PatchHandle const * handle = patchmap.FindPatch(sp.face, u, v);
        assert(handle);
        Far::PatchParam param = patchTable->GetPatchParam(*handle);
        param.Normalize(u,v);

        patchTable->EvaluateBasis(*handle, sp.u[0], sp.u[1], pWeights);
        Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);

        float us[4];
        float vs[4];
        float u2 = u*u;
        float u3 = u2*u;
        us[0] = 1.0f/6.0f * (1.0f - 3.0f*(u - u2) - u3);
        us[1] = 1.0f/6.0f * (4.0f - 6.0f*u2 + 3.0f*u3);
        us[2] = 1.0f/6.0f * (1.0f + 3.0f*(u + u2 - u3));
        us[3] = 1.0f/6.0f * (u3);

        float v2 = v*v;
        float v3 = v2*v;
        vs[0] = 1.0f/6.0f * (1.0f - 3.0f*(v - v2) - v3);
        vs[1] = 1.0f/6.0f * (4.0f - 6.0f*v2 + 3.0f*v3);
        vs[2] = 1.0f/6.0f * (1.0f + 3.0f*(v + v2 - v3));
        vs[3] = 1.0f/6.0f * (v3);


        for(int i = 0; i < cvs.size(); i++) {
            float w = vs[i / 4] * us[i % 4];
            result.col(k) += /*pWeights[i]*/ w * localPts.col(cvs[i]);
        }
    }
    return result;
}

Matrix3X SubdivEvaluator::evalPointsPartial(const Matrix3X& localPts) {
    OpenSubdiv::Far::StencilTable const * localStencilTable = getLocalStencilTable();
    int localPointCount = localStencilTable->GetNumStencils();
    
    //std::cout << "LOCAL STENCIL SIZE: " << localPointCount << std::endl;

    //Matrix3X localPts = localPoints(base);

    
    Matrix3X result(3,controlVertexPtex.size()); result.setZero();
     
    float pWeights[16];
    for(int k = 0; k < controlVertexPtex.size(); k++) {
        auto sp = controlVertexPtex[k];
        float u = sp.u[0];
        float v = sp.u[1];
	if(!((u>=0.0f) && (u<=1.0f) && (v>=0.0f) && (v<=1.0f) )) {
            std::cerr << k << ": "<< u  << "," << v << std::endl;
        }
        Far::PatchTable::PatchHandle const * handle = this->patchMap->FindPatch(sp.face, u, v);
        assert(handle);
        Far::PatchParam param = patchTable->GetPatchParam(*handle);
        param.Normalize(u,v);

        
        patchTable->EvaluateBasis(*handle, sp.u[0], sp.u[1], pWeights);
        Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);

        float us[4];
        float vs[4];
        float u2 = u*u;
        float u3 = u2*u;
        us[0] = 1.0f/6.0f * (1.0f - 3.0f*(u - u2) - u3);
        us[1] = 1.0f/6.0f * (4.0f - 6.0f*u2 + 3.0f*u3);
        us[2] = 1.0f/6.0f * (1.0f + 3.0f*(u + u2 - u3));
        us[3] = 1.0f/6.0f * (u3);

        float v2 = v*v;
        float v3 = v2*v;
        vs[0] = 1.0f/6.0f * (1.0f - 3.0f*(v - v2) - v3);
        vs[1] = 1.0f/6.0f * (4.0f - 6.0f*v2 + 3.0f*v3);
        vs[2] = 1.0f/6.0f * (1.0f + 3.0f*(v + v2 - v3));
        vs[3] = 1.0f/6.0f * (v3);


        for(int i = 0; i < cvs.size(); i++) {
            float w = vs[i / 4] * us[i % 4];
            result.col(k) += /*pWeights[i]*/ w * localPts.col(cvs[i]);
        }
    }
    return result;
}

void SubdivEvaluator::evalPointsPartial(const Matrix3X& localPts, const std::vector<SurfacePoint>& surfacePoints,
	Matrix3X* out_S, Matrix3X* out_Su, Matrix3X* out_Sv) {
    OpenSubdiv::Far::StencilTable const * localStencilTable = getLocalStencilTable();
    int localPointCount = localStencilTable->GetNumStencils();
    
    //std::cout << "LOCAL STENCIL SIZE: " << localPointCount << std::endl;

    //Matrix3X localPts = localPoints(base);
    
    assert(surfacePoints.size() == out_S->cols());
	assert(!out_Su || (surfacePoints.size() == out_Su->cols()));
	assert(!out_Sv || (surfacePoints.size() == out_Sv->cols()));

    Matrix3X result(3,surfacePoints.size()); result.setZero();
    if(out_S != NULL) out_S->setZero();
    if(out_Su != NULL) out_Su->setZero();
    if(out_Sv != NULL) out_Sv->setZero();
    
    
     
    float pWeights[16];
    for(int k = 0; k < surfacePoints.size(); k++) {
        auto sp = surfacePoints[k];
        float u = sp.u[0];
        float v = sp.u[1];
	if(!((u>=0.0f) && (u<=1.0f) && (v>=0.0f) && (v<=1.0f) )) {
            std::cerr << k << ": " << u  << "," << v << std::endl;
        }

        Far::PatchTable::PatchHandle const * handle = patchMap->FindPatch(sp.face, u, v);
        assert(handle);
        Far::PatchParam param = patchTable->GetPatchParam(*handle);
        param.Normalize(u,v);

        
        patchTable->EvaluateBasis(*handle, sp.u[0], sp.u[1], pWeights);
        Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);

        float us[4];
        float vs[4];
        float u2 = u*u;
        float u3 = u2*u;
        us[0] = 1.0f/6.0f * (1.0f - 3.0f*(u - u2) - u3);
        us[1] = 1.0f/6.0f * (4.0f - 6.0f*u2 + 3.0f*u3);
        us[2] = 1.0f/6.0f * (1.0f + 3.0f*(u + u2 - u3));
        us[3] = 1.0f/6.0f * (u3);
        
        
        float dus[4];
		dus[0] = -0.5*u2 +     u    -.5;
		dus[1] =  1.5*u2 - 2.0*u;
		dus[2] = -1.5*u2 +     u    +.5;
		dus[3] =  0.5*u2;

        float v2 = v*v;
        float v3 = v2*v;
        vs[0] = 1.0f/6.0f * (1.0f - 3.0f*(v - v2) - v3);
        vs[1] = 1.0f/6.0f * (4.0f - 6.0f*v2 + 3.0f*v3);
        vs[2] = 1.0f/6.0f * (1.0f + 3.0f*(v + v2 - v3));
        vs[3] = 1.0f/6.0f * (v3);
        
       
		float dvs[4];
		dvs[0] = -0.5*v2 +     v    -.5;
		dvs[1] =  1.5*v2 - 2.0*v;
		dvs[2] = -1.5*v2 +     v    +.5;
		dvs[3] =  0.5*v2;


        for(int i = 0; i < cvs.size(); i++) {
			if(out_S != NULL)
				 out_S->col(k) += vs[i / 4] * us[i % 4] * localPts.col(cvs[i]);
			if(out_Su != NULL)
				 out_Su->col(k) += vs[i / 4] * dus[i % 4] * localPts.col(cvs[i]);
			if(out_Sv != NULL)
				 out_Sv->col(k) += dvs[i / 4] * us[i % 4] * localPts.col(cvs[i]);
           
        }
    }
}

Vector2 SubdivEvaluator::unnormalize(int patch_idx, Vector2 u) {
	/*
    Far::PatchMap patchmap(*patchTable);
    std::vector<Far::PatchMap::Handle>& _handles = patchmap._handles;
    Far::PatchParam param = patchTable->GetPatchParam(_handles[patch_idx]);
    Vector2 result = u;
    float u1 = u[0];
    float u2 = u[1];
    param.Unnormalize(u1,u2);
    result[0] = u1;
    result[1] = u2;
	assert(u == result);
    return result;
    */
	return u;
}

Vector2 SubdivEvaluator::normalize(int face_idx, Vector2 u, int* patch_idx) {
	/*
	//std::cout << patchTable->GetNumPatchesTotal() << ":" << patchTable->GetNumPtexFaces() << std::endl;
    Far::PatchMap patchmap(*patchTable);
    Far::PatchTable::PatchHandle const * handle = patchmap.FindPatch(face_idx, u[0], u[1]);
    Far::PatchParam param = patchTable->GetPatchParam(*handle);
    float u1 = u[0];
    float u2 = u[1];
    *patch_idx = handle->patchIndex;
    param.Normalize(u1, u2);
    Vector2 result;
    result[0] = u1;
    result[1] = u2;
    assert(u == result);
	assert(*patch_idx == face_idx);
    return result;
    */
    *patch_idx = face_idx;
    return u;
}

std::vector<int> SubdivEvaluator::getBase() {
    Far::PatchMap patchmap(*patchTable);
    std::vector<Far::PatchMap::Handle>& _handles = patchmap._handles;
    std::vector<int> patchBase(16*_handles.size());
    for(int i = 0; i < _handles.size(); i++) {
		//std::cout << _handles[i].patchIndex << std::endl;
        Far::ConstIndexArray cvs = patchTable->GetPatchVertices(_handles[i]);
        for(int j = 0; j < 16; j++) {
            patchBase[i*16+j] = cvs[j];
        }
    }
    return patchBase;
}

void SubdivEvaluator::evaluateSubdivSurface(Matrix3X const& vert_coords,
	std::vector<SurfacePoint> const& uv,
	Matrix3X* out_S,
	triplets_t* out_dSdX,
	triplets_t* out_dSudX,
	triplets_t* out_dSvdX,
	Matrix3X* out_Su,
	Matrix3X* out_Sv,
	Matrix3X* out_Suu,
	Matrix3X* out_Suv,
	Matrix3X* out_Svv,
	Matrix3X* out_N,
	Matrix3X* out_Nu,
	Matrix3X* out_Nv,
    bool printDebug) const
{
	// Check it's the same size vertex array
	assert(vert_coords.cols() == nVertices);
	// Check output size matches input
	assert(uv.size() == out_S->cols());
	assert(!out_Su || (uv.size() == out_Su->cols()));
	assert(!out_Sv || (uv.size() == out_Sv->cols()));
	assert(!out_Suu || (uv.size() == out_Suu->cols()));
	assert(!out_Suv || (uv.size() == out_Suv->cols()));
	assert(!out_Svv || (uv.size() == out_Svv->cols()));
	// Check that we can use raw pointers as iterators over the evaluation_verts
	//assert((uint8_t*) &evaluation_verts_buffer[1] - (uint8_t*) &evaluation_verts_buffer[0] == 12);

	if (0) {
		for (int i = 0; i < uv.size(); ++i) {
			out_S->col(i)[0] = uv[i].u[0];
			out_S->col(i)[1] = uv[i].u[1];
			out_S->col(i)[2] = 0.0;

			if (!out_Su) continue;

			out_Su->col(i)[0] = 1;
			out_Su->col(i)[1] = 0;
			out_Su->col(i)[2] = 0;

			out_Sv->col(i)[0] = 0;
			out_Sv->col(i)[1] = 1;
			out_Sv->col(i)[2] = 0;

		}
		return;
	}
	
	Far::LimitStencilTableFactory::LocationArrayVec locationArray;
	std::vector<float> myus;
	std::vector<float> myvs;
	for (size_t i = 0; i < uv.size(); i++) {
		myus.push_back(uv[i].u.x());
		myvs.push_back(uv[i].u.y());
	}
	for (size_t i = 0; i < uv.size(); i++) {
		Far::LimitStencilTableFactory::LocationArray locar;
		locar.ptexIdx = uv[i].face;
		locar.s = &myus[i];
		locar.t = &myvs[i];
		locar.numLocations = 1;
		locationArray.push_back(locar);
	}

	float
		pWeights[MAX_NUM_W],
		dsWeights[MAX_NUM_W],
		dtWeights[MAX_NUM_W],
		dssWeights[MAX_NUM_W],
		dttWeights[MAX_NUM_W],
		dstWeights[MAX_NUM_W];

	// Zero the output arrays
	if (out_S) out_S->setZero();
	if (out_Su) out_Su->setZero();
	if (out_Sv) out_Sv->setZero();
	if (out_Suu) out_Suu->setZero();
	if (out_Suv) out_Suv->setZero();
	if (out_Svv) out_Svv->setZero();
	if (out_N) out_N->setZero();
	if (out_Nu) out_Nu->setZero();
	if (out_Nv) out_Nv->setZero();

	if (uv.size() == 0) {
		return;
	}

	// Preallocate triplet vectors to max feasibly needed
#define CLEAR(VAR)\
  if (VAR) {\
    VAR->reserve(MAX_NUM_W*uv.size());\
    VAR->resize(0);\
  }
	CLEAR(out_dSdX);
	CLEAR(out_dSudX);
	CLEAR(out_dSvdX);
#undef CLEAR


	Far::LimitStencilTableFactory::Options opts;
	bool firstDeriv = out_dSudX || out_dSvdX || out_Su || out_Sv;
	bool secondDeriv = out_Suu || out_Svv || out_Suv;
	opts.generate1stDerivatives = firstDeriv;
	opts.generate2ndDerivatives = secondDeriv;
	Far::LimitStencilTable const * limitTable = Far::LimitStencilTableFactory::Create(
		*refiner, locationArray, cvstencils, patchTable, opts);

	const std::vector<OpenSubdiv::Vtr::Index> idx_list = limitTable->GetControlIndices();
	const std::vector<OpenSubdiv::Vtr::Index> idx_offsets = limitTable->GetOffsets();

	const size_t limit_stencil_count = limitTable->GetNumStencils();
    int maxSize = 0;
	for (size_t i = 0; i < limit_stencil_count; i++) {
		Far::LimitStencil limitStencil = limitTable->GetLimitStencil(i);
		const OpenSubdiv::Vtr::Index* myidx = limitStencil.GetVertexIndices();
		if (limitStencil.GetSize() > MAX_NUM_W) {
			std::cerr << limitStencil.GetSize() << std::endl;
		}
		const int control_range = limitStencil.GetSize();
		const float* w = limitStencil.GetWeights();
		const float* wu = limitStencil.GetDuWeights();
		const float* wv = limitStencil.GetDvWeights();
		const float* wuu = limitStencil.GetDuuWeights();
		const float* wvv = limitStencil.GetDvvWeights();
		const float* wuv = limitStencil.GetDuvWeights();
        if(printDebug) {
            if( control_range > maxSize) {
                maxSize = control_range;
            }
        }
        assert(control_range <= MAX_NUM_W);

		for (int j = 0; j < control_range; j++) {
            if(j >= MAX_NUM_W) {
                std::cerr << "OVERFLOW!!!" << std::endl;
            }
			const OpenSubdiv::Vtr::Index idx = myidx[j];

			const Matrix3X& coord = vert_coords.col(idx);
			if (out_S) out_S->col(i) += w[j] * coord;
			if (out_Su) out_Su->col(i) += wu[j] * coord;
			if (out_Sv) out_Sv->col(i) += wv[j] * coord;
			if (out_Suu) out_Suu->col(i) += wuu[j] * coord;
			if (out_Svv) out_Svv->col(i) += wvv[j] * coord;
			if (out_Suv) out_Suv->col(i) += wuv[j] * coord;
			//if (out_dSdX && w[j] != 0) {
            if (out_dSdX) {
                out_dSdX->push_back(triplet_t(i, idx, w[j]));
            }
			//if (out_dSudX && wu[j] != 0) {
            if (out_dSudX) {
                out_dSudX->push_back(triplet_t(i, idx, wu[j]));
            }
            //if (out_dSvdX && wv[j] != 0) {
			if (out_dSvdX) {
                out_dSvdX->push_back(triplet_t(i, idx, wv[j]));
            }
			//if (out_dSdX) out_dSdX->add(idx, i, w[j]);
		}
	}
    if(printDebug) {
        //std::cout << "Biggest Influence Size: " << maxSize << std::endl;
    }

	delete limitTable;

	return;
}

void SubdivEvaluator::createPatches(const Matrix3X& localPts,
	MeshTopology* mesh, Matrix3X* verts) {

	const int WIDTH = 3;
    const int HEIGHT = 3;
    const int VERTEX_COUNT = WIDTH*HEIGHT*num_faces;
	int k = 0;
	mesh->num_vertices = VERTEX_COUNT;
    mesh->quads.resize(5, (WIDTH-1)*(HEIGHT-1)*num_faces);
    mesh->quads.setConstant(-1);
    for (int f = 0; f < num_faces; ++f) {
		for(int x = 0; x < WIDTH-1; x++) {
        for(int y = 0; y < HEIGHT-1; y++) {
			int j = f*(WIDTH)*(HEIGHT);
    		mesh->quads(0, k) = j+(HEIGHT)*x + y;
			mesh->quads(1, k) = j+(HEIGHT)*(x+1) + y;
			mesh->quads(2, k) = j+(HEIGHT)*(x+1) + y + 1;
			mesh->quads(3, k) = j+(HEIGHT)*x + y + 1;
			k++;
		}
		}
    }
    mesh->update_adjacencies();


	OpenSubdiv::Far::StencilTable const * localStencilTable = getLocalStencilTable();
	float pWeights[16];

	k = 0;
	verts->resize(3,VERTEX_COUNT); verts->setConstant(0);
	for(int f = 0; f < num_faces; f++) {
		for(int x = 0; x < WIDTH; x++) {
		for(int y = 0; y < HEIGHT; y++) {
			float u = x / (WIDTH-1.0);
			float v = y / (HEIGHT-1.0);
			Far::PatchTable::PatchHandle const * handle = patchMap->FindPatch(f, u, v);
        	assert(handle);
        	Far::PatchParam param = patchTable->GetPatchParam(*handle);
        	param.Normalize(u,v);

        
        	patchTable->EvaluateBasis(*handle, u, v, pWeights);
        	Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);

			float us[4];
        	float vs[4];
        	float u2 = u*u;
        	float u3 = u2*u;
        	us[0] = 1.0f/6.0f * (1.0f - 3.0f*(u - u2) - u3);
        	us[1] = 1.0f/6.0f * (4.0f - 6.0f*u2 + 3.0f*u3);
        	us[2] = 1.0f/6.0f * (1.0f + 3.0f*(u + u2 - u3));
        	us[3] = 1.0f/6.0f * (u3);
        
        
        	float dus[4];
        	dus[0] = -0.5*u2 +     u    -.5;
        	dus[1] =  1.5*u2 - 2.0*u;
        	dus[2] = -1.5*u2 +     u    +.5;
        	dus[3] =  0.5*u2;

        	float v2 = v*v;
        	float v3 = v2*v;
        	vs[0] = 1.0f/6.0f * (1.0f - 3.0f*(v - v2) - v3);
        	vs[1] = 1.0f/6.0f * (4.0f - 6.0f*v2 + 3.0f*v3);
        	vs[2] = 1.0f/6.0f * (1.0f + 3.0f*(v + v2 - v3));
        	vs[3] = 1.0f/6.0f * (v3);
         
       
        	float dvs[4];
        	dvs[0] = -0.5*v2 +     v    -.5;
        	dvs[1] =  1.5*v2 - 2.0*v;
        	dvs[2] = -1.5*v2 +     v    +.5;
        	dvs[3] =  0.5*v2;

			for(int i = 0; i < cvs.size(); i++) {
            	verts->col(k) += vs[i / 4] * us[i % 4] * localPts.col(cvs[i]);
        	}
			k++;

			}
			}
	}
}
