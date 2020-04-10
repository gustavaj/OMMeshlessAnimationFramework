#include "Lattice.h"

#include <random>

namespace OML
{
	Lattice::Lattice()
		: Lattice("Lattice")
	{
	}

	Lattice::Lattice(std::string name)
		: m_name(name)
	{
		addLatticeProperties();
	}

	Lattice::~Lattice()
	{
	}

	void Lattice::addPatch(Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight)
	{
		OpenMesh::VertexHandle tl, tr, bl, br;
		for (auto it = vertices_begin(); it != vertices_end(); it++)
		{
			auto p = point(*it);
			if ((p - topLeft).length() < 1e-5) tl = *it;
			if ((p - topRight).length() < 1e-5) tr = *it;
			if ((p - bottomLeft).length() < 1e-5) bl = *it;
			if ((p - bottomRight).length() < 1e-5) br = *it;
		}

		if (!tl.is_valid()) tl = add_vertex(topLeft);
		if (!tr.is_valid()) tr = add_vertex(topRight);
		if (!bl.is_valid()) bl = add_vertex(bottomLeft);
		if (!br.is_valid()) br = add_vertex(bottomRight);

		add_face(tl, bl, br, tr);
	}

	void Lattice::addPatch(Vec2f topLeft, Vec2f topRight, Vec2f bottomLeft, Vec2f bottomRight)
	{
		addPatch(Vec3f(topLeft[0], topLeft[1], 0.0f), Vec3f(topRight[0], topRight[1], 0.0f),
			Vec3f(bottomLeft[0], bottomLeft[1], 0.0f), Vec3f(bottomRight[0], bottomRight[1], 0.0f));
	}

	void Lattice::addPatch(Vec2f topLeft, float width, float height)
	{
		addPatch(Vec3f(topLeft[0], topLeft[1], 0.0f), width, height, 0.0f);
	}

	void Lattice::addPatch(Vec3f topLeft, float width, float height, float depth)
	{
		Vec3f right(width, 0.0f, 0.0f);
		Vec3f down(0.0f, height, 0.0f);
		Vec3f in(0.0f, 0.0f, depth);
		addPatch(topLeft + in, topLeft + right + in, topLeft + down, topLeft + right + down);
	}

	void Lattice::addGrid(Vec2f topLeft, float width, float height, int rows, int cols)
	{
		float w = width / (float)cols;
		float h = height / (float)rows;
		Vec2f dw(w, 0.0f);
		Vec2f dh(0.0f, h);

		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				addPatch(topLeft + j * dw, w, h);
			}
			topLeft += dh;
		}
	}

	void Lattice::addDonut(Vec2f topLeft, float outerRadius, float innerRadius)
	{
		float d = (outerRadius - innerRadius) / 2;
		float r = innerRadius;
		Vec2f dw(d, 0.0f);
		Vec2f dh(0.0f, d);
		Vec2f rw(r, 0.0f);
		Vec2f rh(0.0f, r);

		addPatch(topLeft, d, d);
		addPatch(topLeft + dw, r, d);
		addPatch(topLeft + dw + rw, d, d);

		topLeft += dh;
		addPatch(topLeft, d, r);
		addPatch(topLeft + dw + rw, d, r);

		topLeft += rh;
		addPatch(topLeft, d, d);
		addPatch(topLeft + dw, r, d);
		addPatch(topLeft + dw + rw, d, d);
	}

	void Lattice::addCylinder(Vec3f center, float radius, float height, int rows, int cols)
	{
		float dh = height / (float)rows;
		float d_angle = 2 * M_PI / cols;

		std::vector<std::vector<Vec3f>> coords(rows + 1);

		for (size_t i = 0; i < rows + 1; i++) {
			coords[i] = std::vector<Vec3f>(cols + 1);
			for (size_t j = 0; j < cols; j++) {
				float theta = j * d_angle;
				coords[i][j] = Vec3f(radius * std::cos(theta), dh * i, radius * std::sin(theta)) + center;
			}
			coords[i][cols] = coords[i][0];
		}

		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				addPatch(coords[i + 1][j], coords[i + 1][j + 1], coords[i][j], coords[i][j + 1]);
			}
		}
	}

	void Lattice::addSphere(Vec3f center, float radius, int segments, int slices)
	{
		std::vector<std::vector<Vec3f>> points(slices + 1);

		// https://stackoverflow.com/questions/43412525/algorithm-to-draw-a-sphere-using-quadrilaterals
		for (int slice = 0; slice <= slices; slice++)
		{
			points[slice] = std::vector<Vec3f>(segments + 1);

			float lat0 = M_PI * (((float)slice / slices) - 0.5f);
			if (slice == 0) lat0 = M_PI * (0.001f - 0.5f);
			else if (slice == slices) lat0 = M_PI * (0.999f - 0.5f);
			float y0 = std::sin(lat0) * radius;
			if (slice == 0) y0 = -radius * ((float)(slices - 0.5f) / slices);
			else if (slice == slices) y0 = radius * ((float)(slices - 0.5f) / slices);
			float yr0 = std::cos(lat0);

			for (int segment = 0; segment < segments; segment++)
			{
				float long0 = 2.0f * M_PI * ((float)segment / segments);
				float x0 = std::cos(long0) * radius;
				float z0 = std::sin(long0) * radius;

				points[slice][segment] = Vec3f(x0 * yr0, y0, z0 * yr0) + center;
			}

			points[slice][segments] = points[slice][0];
		}

		for (size_t i = 0; i < slices; i++) {
			for (size_t j = 0; j < segments; j++) {
				addPatch(points[i + 1][j], points[i + 1][j + 1], points[i][j], points[i][j + 1]);
			}
		}
	}

	void Lattice::removeNormalSinSimulation()
	{
		m_simulators.resize(0);
		resetMatrices();
	}

	void Lattice::addNormalSinSimulation()
	{
		resetMatrices();

		m_simulators.resize(0);
		for (auto& locus : m_loci)
		{
			double t = static_cast<double>(rand() % 100);
			double amp, speed;
			if (m_maxAmp > m_minAmp) {
				amp = m_minAmp + static_cast<double>(rand() % (int)(m_maxAmp - m_minAmp));
			}
			else {
				amp = static_cast<double>(rand() % (int)m_maxAmp);
			}
			if (m_maxSpeed > m_minSpeed) {
				speed = m_minSpeed + static_cast<double>(rand() % (int)(m_maxSpeed - m_minSpeed));
			}
			else {
				speed = static_cast<double>(rand() % (int)m_maxSpeed);
			}
			m_simulators.push_back(NormalSinSimulator(0.0, amp, speed, locus.normal));
		}
	}

	void Lattice::update(double dt)
	{
		if (m_animate) 
		{
			for (size_t i = 0; i < m_simulators.size(); i++)
			{
				m_simulators[i].simulate(dt, m_matrices[i]);
			}
		}

		localUpdate(dt);
	}

	void Lattice::induceLattice()
	{
		// Set the edge color of gridlines depending on if the edge is on the boundary or not
		setupEdgeColor();

		// Setup the valence property of the loci, and set the color of gridpoints based on the point's valence
		setupLociValenceAndPointColor();

		// Find and resolve T-loci, and add local surface for them
		handleTLoci();

		// Add regular local surfaces for each loci, and setup patches
		setupLocalSurfacesAndPatches();
	}

	void Lattice::resetMatrices()
	{
		m_matrices.resize(0);
		for (auto mat : m_initialMatrices)
		{
			m_matrices.push_back(mat);
		}
	}

	void Lattice::setupEdgeColor()
	{
		// Set the color of the edges based on if they are on the boundary or not
		for (auto e_itr = edges_begin(); e_itr != edges_end(); e_itr++)
		{
			set_color(*e_itr, is_boundary(*e_itr) ? BOUNDARY_EDGE_COLOR : INNER_EDGE_COLOR);
		}
	}

	void Lattice::setupLociValenceAndPointColor()
	{
		// Loop over vertices, set valence and color based on valence.
		for (auto v_itr = vertices_begin(); v_itr != vertices_end(); v_itr++)
		{
			auto vh = (*v_itr);
			auto valenceProp = property(LatticeProperties::VertexValence, vh);

			valenceProp = valence(vh);

			property(LatticeProperties::VertexValence, vh) = valenceProp;

			LocusType type = LocusType::Unresolved;
			switch (valenceProp)
			{
			case 2: type = LocusType::Corner; break;
			case 3: type = LocusType::Boundary; break;
			case 4: type = LocusType::Inner; break;
			}
			property(LatticeProperties::Type, vh) = type;

			set_color(vh, LOCUS_VALENCE_COLOR[valenceProp]);

			// Initialize local surface idx to -1
			property(LatticeProperties::LocusIndex, vh) = -1;
		}
	}

	void Lattice::handleTLoci()
	{
		// Look for T-loci, and handle it!
		//for (auto v_itr = vertices_begin(); v_itr != vertices_end(); v_itr++)
		//{
		//	auto vh = *v_itr;
		//	if (property(LatticeProperties::LocusValence, vh) == 5
		//		|| (property(LatticeProperties::LocusValence, vh) == 4 && is_boundary(vh)))
		//	{
		//		// Find the opposite extreme point vertex
		//		OpenMesh::VertexHandle vh2, vht;
		//		for (auto vv_itr = vv_iter(vh); vv_itr.is_valid(); vv_itr++)
		//		{
		//			if (property(LatticeProperties::LocusValence, *vv_itr) == 5 ||
		//				(property(LatticeProperties::LocusValence, *vv_itr) == 4 && is_boundary(*vv_itr)))
		//			{
		//				vh2 = *vv_itr;
		//				break;
		//			}
		//		}

		//		auto heh12 = find_halfedge(vh, vh2);

		//		OpenMesh::HalfedgeHandle heht1, heht2;

		//		// Find the other t-point
		//		for (auto vv_itr = vv_iter(vh); vv_itr.is_valid(); vv_itr++)
		//		{
		//			heht1 = find_halfedge(*vv_itr, vh);
		//			heht2 = find_halfedge(*vv_itr, vh2);
		//			if (heht1.is_valid() && heht2.is_valid())
		//			{
		//				vht = *vv_itr;
		//				break;
		//			}
		//		}

		//		if (vht.is_valid()) {
		//			// Change color of vertices
		//			// Yellow for extreme points, purple for t-poitns
		//			set_color(vh, Col3(255, 255, 0));
		//			set_color(vh2, Col3(255, 255, 0));
		//			set_color(vht, Col3(255, 0, 255));

		//			// heh12
		//			auto heh21 = opposite_halfedge_handle(heh12);

		//			auto vh3 = from_vertex_handle(prev_halfedge_handle(heh21));
		//			auto vh4 = to_vertex_handle(next_halfedge_handle(heh21));

		//			delete_edge(edge_handle(heh12));

		//			garbage_collection(true, true, true);

		//			std::vector<OpenMesh::VertexHandle> vhs{ vh4, vh3, vh2, vht, vh };
		//			add_face(vhs);

		//			// Update valence
		//			property(LatticeProperties::LocusValence, vh) -= 1;
		//			property(LatticeProperties::LocusValence, vh2) -= 1;


		//			auto heh_1_to_t = find_halfedge(vh, vht);
		//			auto heh_t_to_2 = find_halfedge(vht, vh2);

		//			// Create local surface and update indices for all affected parts.
		//			auto vh_before = from_vertex_handle(prev_halfedge_handle(opposite_halfedge_handle(prev_halfedge_handle(heh_1_to_t))));
		//			auto vh_after = to_vertex_handle(next_halfedge_handle(opposite_halfedge_handle(next_halfedge_handle(heh_t_to_2))));

		//			auto v1 = calc_edge_vector(opposite_halfedge_handle(heh_1_to_t)) + calc_edge_vector(find_halfedge(vh, vh_before));
		//			auto v2 = calc_edge_vector(heh_t_to_2) + calc_edge_vector(find_halfedge(vh2, vh_after));

		//			auto u1 = calc_edge_vector(next_halfedge_handle(opposite_halfedge_handle(heh_1_to_t)));
		//			auto u2 = calc_edge_vector(opposite_halfedge_handle(prev_halfedge_handle(heh_1_to_t)));
		//			auto u3 = calc_edge_vector(opposite_halfedge_handle(prev_halfedge_handle(opposite_halfedge_handle(heh_t_to_2))));
		//			auto u4 = calc_edge_vector(next_halfedge_handle(heh_t_to_2));

		//			addLocalSurface(v1 + u1, v1 + u2, v2 + u3, v2 + u4, point(vht));
		//			int ls_idx = m_localSurfaces.size() - 1;
		//			property(LatticeProperties::LocusLocalSurfaceIdx, vh) = ls_idx;
		//			property(LatticeProperties::LocusLocalSurfaceIdx, vht) = ls_idx;
		//			property(LatticeProperties::LocusLocalSurfaceIdx, vh2) = ls_idx;
		//		}
		//	}
		//}
	}

	void Lattice::setupLocalSurfacesAndPatches()
	{
		// init values
		m_numControlPoints = 0;
		m_numLoci = 0;
		m_numPatches = 0;

		// Add local surfaces
		for (auto f_itr = faces_begin(); f_itr != faces_end(); f_itr++)
		{
			// Face
			auto fh = *f_itr;
			auto f = face(fh);

			// Halfedge/edge handles
			auto topheh = halfedge_handle(fh);
			auto topeh = edge_handle(topheh);

			auto leftheh = next_halfedge_handle(topheh);
			auto lefteh = edge_handle(leftheh);

			auto bottomheh = next_halfedge_handle(leftheh);
			auto boteh = edge_handle(bottomheh);

			auto rightheh = next_halfedge_handle(bottomheh);
			auto righteh = edge_handle(rightheh);

			// Add local surfaces on vertices
			auto vh1 = to_vertex_handle(topheh); // Top Left
			auto vh2 = to_vertex_handle(rightheh); // Top Right
			auto vh3 = to_vertex_handle(leftheh); // Bottom left
			auto vh4 = to_vertex_handle(bottomheh); // Bottom right

			// Local surfaces are handled differently based on where it is in the patch
			if (property(LatticeProperties::LocusIndex, vh1) == -1)
				addLocusOnVertex(vh1, vh3, vh2, 1);
			if (property(LatticeProperties::LocusIndex, vh2) == -1)
				addLocusOnVertex(vh2, vh1, vh4, 2);
			if (property(LatticeProperties::LocusIndex, vh3) == -1)
				addLocusOnVertex(vh3, vh4, vh1, 3);
			if (property(LatticeProperties::LocusIndex, vh4) == -1)
				addLocusOnVertex(vh4, vh2, vh3, 4);

			// Setup patch vertices for face.
			Patch patch;
			patch.lociIndices = {
				property(LatticeProperties::LocusIndex, vh1),
				property(LatticeProperties::LocusIndex, vh2),
				property(LatticeProperties::LocusIndex, vh3),
				property(LatticeProperties::LocusIndex, vh4)
			};
			patch.fh = fh;
			m_patches.push_back(std::move(patch));

			m_boundaries.push_back(m_loci[patch.lociIndices[0]].faceMappings[fh]);
			m_boundaries.push_back(m_loci[patch.lociIndices[1]].faceMappings[fh]);
			m_boundaries.push_back(m_loci[patch.lociIndices[2]].faceMappings[fh]);
			m_boundaries.push_back(m_loci[patch.lociIndices[3]].faceMappings[fh]);

			setupPatchVertices(patch);

			m_numPatches++;
		}
	}

	void Lattice::addLatticeProperties()
	{
		add_property(LatticeProperties::VertexValence);
		add_property(LatticeProperties::LocusIndex);
		add_property(LatticeProperties::Type);
	}

	void Lattice::addLocus(
		OpenMesh::VertexHandle vertex, std::vector<Vec3f>& controlPoints,
		std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo>& faceMappings, Vec3f offset)
	{
		Locus locus;
		locus.controlPointIndex = m_numControlPoints;
		locus.controlPointCount = controlPoints.size();
		locus.vh = vertex;
		locus.faceMappings = faceMappings;
		locus.matrixIndex = m_numLoci;

		// Find the average normal vector over all the patches for the locus
		Vec3f p = point(vertex);
		Vec3f normal = Vec3f(0, 0, 0);
		auto first_adj_vertex = vv_iter(vertex);
		auto type = property(LatticeProperties::Type, vertex);
		bool loopAround = !(type == LocusType::Corner || type == LocusType::Boundary);
		//Vec3f first_adj_point = point(*first_adj_vertex);
		for (auto adj_vertex = first_adj_vertex; adj_vertex.is_valid(); ) {
			Vec3f v1 = (point(*adj_vertex) - p);
			adj_vertex++;
			if (!adj_vertex.is_valid() && !loopAround)
				break;
			Vec3f next_adj_vertex = adj_vertex.is_valid() ? point(*adj_vertex) : point(*first_adj_vertex);
			Vec3f v2 = (next_adj_vertex - p);
			normal += v1 % v2;
		}
		locus.normal = normal.normalize();

		property(LatticeProperties::LocusIndex, vertex) = m_numLoci;

		m_loci.push_back(locus);
		m_numLoci++;

		for (auto& p : controlPoints) {
			m_controlPoints.push_back(glm::vec4(p[0], p[1], p[2], 1.0f));
		}
		//m_patchUniforms.controlPoints.insert(m_patchUniforms.controlPoints.end(), controlPoints.begin(), controlPoints.end());
		m_numControlPoints += controlPoints.size();

		m_matrices.push_back(
			glm::translate(glm::mat4(1.0f), glm::vec3(offset[0], offset[1], offset[2])));
		m_initialMatrices.push_back(
			glm::translate(glm::mat4(1.0f), glm::vec3(offset[0], offset[1], offset[2])));

		setupLocalSurfaceVertex(locus);
	}

	void Lattice::addLocusOnVertex(
		OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
		OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace)
	{
		auto valence = property(LatticeProperties::VertexValence, vertex);
		switch (valence) {
			case 2: { addLocusOnCornerVertex(vertex, next_vertex, prev_vertex, vertexIndexOnFace); break; }
			case 3: { addLocusOnBoundaryVertex(vertex, next_vertex, prev_vertex, vertexIndexOnFace); break; }
			case 4: { addLocusOnInnerVertex(vertex, next_vertex, prev_vertex, vertexIndexOnFace); break; }
		}
	}

	void Lattice::addLocusOnCornerVertex(
		OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
		OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace)
	{
		auto heh_locus_to_next = find_halfedge(vertex, next_vertex);
		Vec3f offset = point(vertex);
		Vec3f zero(0.0f, 0.0f, 0.0f);
		Vec3f o_to_n = calc_edge_vector(heh_locus_to_next);
		Vec3f o_to_n_next = calc_edge_vector(next_halfedge_handle(heh_locus_to_next));
		Vec3f o_to_p = -calc_edge_vector(prev_halfedge_handle(heh_locus_to_next));
		Vec3f v_inner = o_to_n + o_to_n_next;

		std::vector<Vec3f> controlPoints;
		switch (vertexIndexOnFace) {
			case 1: { controlPoints = createLocalSurfaceControlPoints(zero, o_to_p, o_to_n, v_inner); break; }
			case 2: { controlPoints = createLocalSurfaceControlPoints(o_to_n, zero, v_inner, o_to_p); break; }
			case 3: { controlPoints = createLocalSurfaceControlPoints(o_to_p, v_inner, zero, o_to_n); break; }
			case 4: { controlPoints = createLocalSurfaceControlPoints(v_inner, o_to_n, o_to_p, zero); break; }
		}

		std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo> faceMappings;
		faceMappings.insert({ face_handle(heh_locus_to_next), BoundaryInfo(0.0f, 1.0f, 0.0f, 1.0f) });
		addLocus(vertex, controlPoints, faceMappings, offset);
	}

	void Lattice::addLocusOnBoundaryVertex(
		OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
		OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace)
	{
		Vec3f offset = point(vertex);
		Vec3f zero(0.0f, 0.0f, 0.0f);

		auto heh_next = find_halfedge(vertex, next_vertex);

		BoundaryInfo boundaryFace(0.0f, 1.0f, 0.0f, 1.0f);
		BoundaryInfo boundaryAdjFace(0.0f, 1.0f, 0.0f, 1.0f);

		std::vector<Vec3f> controlPoints;
		std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo> faceMappings;

		if (is_boundary(edge_handle(heh_next))) {
			auto heh_prev = find_halfedge(vertex, prev_vertex);
			Vec3f n = calc_edge_vector(heh_next); // vector from locus to next_locus
			Vec3f no = -calc_edge_vector(prev_halfedge_handle(heh_prev)); // vector from locus to the locus on the adjacent face
			Vec3f nop = -calc_edge_vector(next_halfedge_handle(next_halfedge_handle(heh_prev)));
			Vec3f p = calc_edge_vector(heh_prev);
			Vec3f np = calc_edge_vector(next_halfedge_handle(heh_next));

			Vec3f p_2 = p * 0.5f;
			Vec3f pnop = nop;
			Vec3f pnop_2 = pnop * 0.5f;
			Vec3f pnp = np;
			Vec3f pnp_2 = pnp * 0.5f;

			//float halfwayPoint = 0.5f + (n.length() / (n - no).length() - 0.5f) / 2;
			float halfwayPoint = 0.5f;
			switch (vertexIndexOnFace) {
				case 1: { 
					controlPoints = createLocalSurfaceControlPoints(no, no + pnop_2, no + pnop, zero, p_2, p, n, n + pnp_2, n + pnp);
					boundaryFace.vs = halfwayPoint;
					boundaryAdjFace.ve = halfwayPoint;
					break; 
				}
				case 2: { 
					controlPoints = createLocalSurfaceControlPoints(n, zero, no, n + pnp_2, p_2, no + pnop_2, n + pnp, p, no + pnop);
					boundaryFace.ue = halfwayPoint;
					boundaryAdjFace.us = halfwayPoint;
					break; 
				}
				case 3: { 
					controlPoints = createLocalSurfaceControlPoints(no + pnop, p, n + pnp, no + pnop_2, p_2, n + pnp_2, no, zero, n);
					boundaryFace.us = halfwayPoint;
					boundaryAdjFace.ue = halfwayPoint;
					break; 
				}
				case 4: { 
					controlPoints = createLocalSurfaceControlPoints(n + pnp, n + pnp_2, n, p, p_2, zero, no + pnop, no + pnop_2, no); 
					boundaryFace.ve = halfwayPoint;
					boundaryAdjFace.vs = halfwayPoint;
					break; 
				}
			}

			faceMappings.insert({ face_handle(heh_next), boundaryFace });
			faceMappings.insert({ face_handle(heh_prev), boundaryAdjFace });
		}
		else {
			auto heh_prev = find_halfedge(prev_vertex, vertex);
			Vec3f p = -calc_edge_vector(heh_prev);
			Vec3f po = calc_edge_vector(next_halfedge_handle(opposite_halfedge_handle(heh_next)));
			Vec3f pon = calc_edge_vector(next_halfedge_handle(next_halfedge_handle(opposite_halfedge_handle(heh_next))));
			Vec3f n = calc_edge_vector(heh_next);
			Vec3f pn = -calc_edge_vector(prev_halfedge_handle(heh_prev));

			Vec3f n_2 = n * 0.5f;
			Vec3f npn = pn;
			Vec3f npn_2 = npn * 0.5f;
			Vec3f npon = pon;
			Vec3f npon_2 = npon * 0.5f;

			//float halfwayPoint = 0.5f + (p.length() / (p - po).length() - 0.5f) / 2;
			float halfwayPoint = 0.5f;

			switch (vertexIndexOnFace) {
				case 1: { 
					controlPoints = createLocalSurfaceControlPoints(po, zero, p, po + npon_2, n_2, p + npn_2, po + npon, n, p + npn);
					boundaryFace.us = halfwayPoint;
					boundaryAdjFace.ue = halfwayPoint;
					break; 
				}
				case 2: { 
					controlPoints = createLocalSurfaceControlPoints(po + npon, po + npon_2, po, n, n_2, zero, p + npn, p + npn_2, p);
					boundaryFace.vs = halfwayPoint;
					boundaryAdjFace.ve = halfwayPoint;
					break; 
				}
				case 3: { 
					controlPoints = createLocalSurfaceControlPoints(p, p + npn_2, p + npn, zero, n_2, n, po, po + npon_2, po + npon);
					boundaryFace.ve = halfwayPoint;
					boundaryAdjFace.vs = halfwayPoint;
					break; 
				}
				case 4: { 
					controlPoints = createLocalSurfaceControlPoints(p + npn, n, po + npon, p + npn_2, n_2, po + npon_2, p, zero, po);
					boundaryFace.ue = halfwayPoint;
					boundaryAdjFace.us = halfwayPoint;
					break; 
				}
			}

			faceMappings.insert({ face_handle(heh_next), boundaryFace });
			faceMappings.insert({ face_handle(opposite_halfedge_handle(heh_next)), boundaryAdjFace });
		}

		addLocus(vertex, controlPoints, faceMappings, offset);
	}

	void Lattice::addLocusOnInnerVertex(
		OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
		OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace)
	{
		Vec3f offset = point(vertex);
		Vec3f zero(0.0f, 0.0f, 0.0f);

		auto hehf1 = find_halfedge(vertex, next_vertex);
		auto hehf2 = opposite_halfedge_handle(hehf1);
		auto hehf3 = prev_halfedge_handle(opposite_halfedge_handle(prev_halfedge_handle(hehf1)));
		auto hehf4 = next_halfedge_handle(opposite_halfedge_handle(next_halfedge_handle(hehf2)));

		BoundaryInfo boundaryF1, boundaryF2, boundaryF3, boundaryF4;

		auto u1 = calc_edge_vector(next_halfedge_handle(hehf4));
		auto u2 = -calc_edge_vector(prev_halfedge_handle(hehf3));
		auto u3 = calc_edge_vector(next_halfedge_handle(hehf2));
		auto u4 = -calc_edge_vector(prev_halfedge_handle(hehf1));
		auto u5 = -calc_edge_vector(prev_halfedge_handle(hehf2));
		auto u6 = calc_edge_vector(next_halfedge_handle(hehf1));
		auto v1 = calc_edge_vector(hehf4);
		auto v2 = calc_edge_vector(hehf1);

		auto u42 = u2;
		auto u31 = u1;
		auto u46 = u6;
		auto u35 = u5;

		//float halfwayPointU = 0.5f + (u4.length() / (u3 - u4).length() - 0.5f) / 2;
		//float halfwayPointV = 0.5f + (v2.length() / (v1 - v2).length() - 0.5f) / 2;
		float halfwayPointU = 0.5f;
		float halfwayPointV = 0.5f;
		std::vector<Vec3f> controlPoints;

		switch (vertexIndexOnFace) {
			case 1: { 
				controlPoints = createLocalSurfaceControlPoints(u31 + v1, v1, u42 + v1, u3, zero, u4, u35 + v2, v2, u46 + v2);
				boundaryF1.us = halfwayPointU;
				boundaryF1.vs = halfwayPointV;
				boundaryF2.ue = halfwayPointU;
				boundaryF2.vs = halfwayPointV;
				boundaryF3.us = halfwayPointU;
				boundaryF3.ve = halfwayPointV;
				boundaryF4.ue = halfwayPointU;
				boundaryF4.ve = halfwayPointV;
				break; 
			}
			case 2: { 
				controlPoints = createLocalSurfaceControlPoints(u35 + v2, u3, u31 + v1, v2, zero, v1, u46 + v2, u4, u42 + v1);
				boundaryF1.ue = halfwayPointU;
				boundaryF1.vs = halfwayPointV;
				boundaryF2.ue = halfwayPointU;
				boundaryF2.ve = halfwayPointV;
				boundaryF3.us = halfwayPointU;
				boundaryF3.vs = halfwayPointV;
				boundaryF4.us = halfwayPointU;
				boundaryF4.ve = halfwayPointV;
				break;
			}
			case 3: {
				controlPoints = createLocalSurfaceControlPoints(u42 + v1, u4, u46 + v2, v1, zero, v2, u31 + v1, u3, u35 + v2);
				boundaryF1.us = halfwayPointU;
				boundaryF1.ve = halfwayPointV;
				boundaryF2.us = halfwayPointU;
				boundaryF2.vs = halfwayPointV;
				boundaryF3.ue = halfwayPointU;
				boundaryF3.ve = halfwayPointV;
				boundaryF4.ue = halfwayPointU;
				boundaryF4.vs = halfwayPointV;
				break;
			}
			case 4: {
				controlPoints = createLocalSurfaceControlPoints(u46 + v2, v2, u35 + v2, u4, zero, u3, u42 + v1, v1, u31 + v1);
				boundaryF1.ue = halfwayPointU;
				boundaryF1.ve = halfwayPointV;
				boundaryF2.us = halfwayPointU;
				boundaryF2.ve = halfwayPointV;
				boundaryF3.ue = halfwayPointU;
				boundaryF3.vs = halfwayPointV;
				boundaryF4.us = halfwayPointU;
				boundaryF4.vs = halfwayPointV;
				break;
			}
		}

		std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo> faceMappings;
		faceMappings.insert({ face_handle(hehf1), boundaryF1 });
		faceMappings.insert({ face_handle(hehf2), boundaryF2 });
		faceMappings.insert({ face_handle(hehf3), boundaryF3 });
		faceMappings.insert({ face_handle(hehf4), boundaryF4 });

		addLocus(vertex, controlPoints, faceMappings, offset);
	}

	std::vector<Vec3f> Lattice::createLocalSurfaceControlPoints(
		Vec3f topLeft, Vec3f topRight,
		Vec3f bottomLeft, Vec3f bottomRight)
	{
		Vec3f topMiddle = (topLeft + topRight) / 2;
		Vec3f middleLeft = (topLeft + bottomLeft) / 2;
		Vec3f middleRight = (topRight + bottomRight) / 2;
		Vec3f bottomMiddle = (bottomLeft + bottomRight) / 2;
		Vec3f middle = ((topMiddle + bottomMiddle) / 2 + (middleLeft + middleRight) / 2) / 2;
		return {
			topLeft, topMiddle, topRight,
			middleLeft, middle, middleRight,
			bottomLeft, bottomMiddle, bottomRight
		};
	}

	std::vector<Vec3f> Lattice::createLocalSurfaceControlPoints(
		Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
		Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
		Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight)
	{
		return {
			topLeft, topMiddle, topRight,
			middleLeft, middle, middleRight,
			bottomLeft, bottomMiddle, bottomRight
		};
	}
}