#include "Lattice.h"

#include <chrono>
#include <iostream>

namespace OML
{
	Lattice::Lattice()
		: Lattice("Lattice")
	{
	}

	Lattice::Lattice(std::string name)
		: m_name(name), m_matrix(1.0f), m_color(0.8f, 0.2f, 0.4f)
	{
		addLatticeProperties();
	}

	Lattice::~Lattice()
	{
	}

	void Lattice::addPatch(Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight)
	{
		// TODO: Should probably try to improve this, or not.
		glm::vec3 TL = glm::vec3(topLeft[0], topLeft[1], topLeft[2]);
		glm::vec3 TR = glm::vec3(topRight[0], topRight[1], topRight[2]);
		glm::vec3 BL = glm::vec3(bottomLeft[0], bottomLeft[1], bottomLeft[2]);
		glm::vec3 BR = glm::vec3(bottomRight[0], bottomRight[1], bottomRight[2]);

		OpenMesh::VertexHandle tlvh, trvh, blvh, brvh;

		auto TLItr = m_uniquePointsIndexMap.insert({ TL, m_numUniquePointsAdded });
		if (TLItr.second) {
			tlvh = add_vertex(topLeft);
			m_numUniquePointsAdded++;
		}
		else tlvh = vertex_handle(TLItr.first->second);

		auto TRItr = m_uniquePointsIndexMap.insert({ TR, m_numUniquePointsAdded });
		if (TRItr.second) {
			trvh = add_vertex(topRight);
			m_numUniquePointsAdded++;
		}
		else trvh = vertex_handle(TRItr.first->second);

		auto BLItr = m_uniquePointsIndexMap.insert({ BL, m_numUniquePointsAdded });
		if (BLItr.second) {
			blvh = add_vertex(bottomLeft);
			m_numUniquePointsAdded++;
		}
		else blvh = vertex_handle(BLItr.first->second);

		auto BRItr = m_uniquePointsIndexMap.insert({ BR, m_numUniquePointsAdded });
		if (BRItr.second) {
			brvh = add_vertex(bottomRight);
			m_numUniquePointsAdded++;
		}
		else brvh = vertex_handle(BRItr.first->second);

		add_face(tlvh, blvh, brvh, trvh);

		/*OpenMesh::VertexHandle tl, tr, bl, br;
		for (auto it = vertices_begin(); it != vertices_end(); it++)
		{
			auto p = point(*it);
			if (equal(p, topLeft)) tl = *it;
			if (equal(p, topRight)) tr = *it;
			if (equal(p, bottomLeft)) bl = *it;
			if (equal(p, bottomRight)) br = *it;
		}

		if (!tl.is_valid()) tl = add_vertex(topLeft);
		if (!tr.is_valid()) tr = add_vertex(topRight);
		if (!bl.is_valid()) bl = add_vertex(bottomLeft);
		if (!br.is_valid()) br = add_vertex(bottomRight);

		add_face(tl, bl, br, tr);*/
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
		/*
			Improvements of this for a 31x31 grid:
			-time start: ~3000ms
			-improve point equality check in addPatch(), dont use .length(): ~800ms
			-dont compare points, instead keep a map of unique points: ~60ms
		*/

		auto start = std::chrono::high_resolution_clock::now();

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

		auto end = std::chrono::high_resolution_clock::now();
		auto time = std::chrono::duration<double, std::milli>(end - start).count();

		std::cout << "Lattice::addGrid(rows: " << rows << ", cols: " << cols << ") time: " << time << "ms" << std::endl;
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

	void Lattice::removeSimulator(SimulatorTypes simulatorType)
	{
		auto it = m_simulators.find(simulatorType);
		if (it != m_simulators.end()) {
			for (auto& sim : m_simulators[simulatorType]) {
				sim.second->undoTransformation(m_matrices[sim.first]);
			}
			m_simulators.erase(it);
		}
	}

	void Lattice::addNormalSinSimulation()
	{
		removeSimulator(SimulatorTypes::NormalSin);

		std::unordered_map<uint32_t, std::shared_ptr<Simulator>> sims;
		for (auto& locus : m_loci)
		{
			sims.insert({ locus.matrixIndex, 
				std::make_shared<NormalSinSimulator>(glm::vec3(locus.normal[0], locus.normal[1], locus.normal[2])) });
		}
		m_simulators.insert({ SimulatorTypes::NormalSin, sims });
	}

	void Lattice::addRandomSphereSimulation()
	{
		removeSimulator(SimulatorTypes::RandomSphere);

		std::unordered_map<uint32_t, std::shared_ptr<Simulator>> sims;
		for (auto& locus : m_loci)
		{
			sims.insert({ locus.matrixIndex,
				std::make_shared<RandomSphereSimulator>() });
		}
		m_simulators.insert({ SimulatorTypes::RandomSphere, sims });
	}

	void Lattice::addNormalRotationSimulation()
	{
		removeSimulator(SimulatorTypes::Rotation);

		std::unordered_map<uint32_t, std::shared_ptr<Simulator>> sims;
		for (auto& locus : m_loci)
		{
			sims.insert({ locus.matrixIndex,
				std::make_shared<RangeRotationSimulator>(glm::vec3(locus.normal[0], locus.normal[1], locus.normal[2])) });
		}
		m_simulators.insert({ SimulatorTypes::Rotation, sims });
	}

	void Lattice::addXYScalingSimulation()
	{
		removeSimulator(SimulatorTypes::XYScale);

		std::unordered_map<uint32_t, std::shared_ptr<Simulator>> sims;
		for (auto& locus : m_loci)
		{
			sims.insert({ locus.matrixIndex, std::make_shared<XYScalingSimulator>() });
		}
		m_simulators.insert({ SimulatorTypes::XYScale, sims });
	}

	void Lattice::update(double dt)
	{
		if (m_simulate)
		{
			for (auto& simT : m_simulators)
			{
				for (auto& sim : simT.second) {
					sim.second->simulate(dt, m_matrices[sim.first]);
				}
			}
		}

		localUpdate(dt);
	}

	void Lattice::induceLattice()
	{
		auto start = std::chrono::high_resolution_clock::now();

		std::cout << "Number of unique points added: " << m_numUniquePointsAdded << std::endl;
		// Clear out the pointIndexMap
		m_uniquePointsIndexMap.clear();

		// Setup the valence property of the loci, and set the color of gridpoints based on the point's valence
		setupLociValenceAndPointColor();

		// Find and resolve T-loci, and add local surface for them
		//handleTLoci();

		// Set the edge color of gridlines depending on if the edge is on the boundary or not
		setupEdgeColor();

		// Add regular local surfaces for each loci, and setup patches
		setupLocalSurfacesAndPatches();

		// Clean up boundary map
		m_boundaryMap.clear();
		std::cout << "Unique boundary infos added: " << m_numUniqueBoundaries << std::endl;

		auto end = std::chrono::high_resolution_clock::now();
		auto time = std::chrono::duration<double, std::milli>(end - start).count();

		std::cout << "Lattice::induceLattice() time: " << time << "ms" << std::endl;
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
		for (auto v_itr = vertices_begin(); v_itr != vertices_end(); v_itr++)
		{
			auto vh_terminal = *v_itr;
			// Look for a terminal point for a T-vertex
			if (property(LatticeProperties::VertexValence, vh_terminal) == 5
				|| (property(LatticeProperties::VertexValence, vh_terminal) == 4 && is_boundary(vh_terminal)))
			{
				// Find the other terminal point vertex
				OpenMesh::VertexHandle vh_terminal2, vh_T;
				for (auto vv_itr = vv_iter(vh_terminal); vv_itr.is_valid(); vv_itr++)
				{
					if (property(LatticeProperties::VertexValence, *vv_itr) == 5 ||
						(property(LatticeProperties::VertexValence, *vv_itr) == 4 && is_boundary(*vv_itr)))
					{
						vh_terminal2 = *vv_itr;
						break;
					}
				}

				auto heh12 = find_halfedge(vh_terminal, vh_terminal2);

				OpenMesh::HalfedgeHandle heht1, heht2;

				// Find the T-point
				for (auto vv_itr = vv_iter(vh_terminal); vv_itr.is_valid(); vv_itr++)
				{
					heht1 = find_halfedge(*vv_itr, vh_terminal);
					heht2 = find_halfedge(*vv_itr, vh_terminal2);
					if (heht1.is_valid() && heht2.is_valid())
					{
						vh_T = *vv_itr;
						break;
					}
				}

				if (vh_T.is_valid()) {
					// Change color of vertices
					// Yellow for terminal points, purple for t-poitns
					set_color(vh_terminal, Col3(255, 255, 0));
					set_color(vh_terminal2, Col3(255, 255, 0));
					set_color(vh_T, Col3(255, 0, 255));
					property(LatticeProperties::Type, vh_terminal) = LocusType::T_Terminal;
					property(LatticeProperties::Type, vh_terminal2) = LocusType::T_Terminal;
					property(LatticeProperties::Type, vh_T) = LocusType::T;

					// heh12
					auto heh21 = opposite_halfedge_handle(heh12);

					auto vh3 = from_vertex_handle(prev_halfedge_handle(heh21));
					auto vh4 = to_vertex_handle(next_halfedge_handle(heh21));

					delete_edge(edge_handle(heh12));

					garbage_collection(false, true, true);

					std::vector<OpenMesh::VertexHandle> vhs{ vh4, vh3, vh_terminal2, vh_T, vh_terminal };
					add_face(vhs);

					// Update valence
					property(LatticeProperties::VertexValence, vh_terminal) -= 1;
					property(LatticeProperties::VertexValence, vh_terminal2) -= 1;

					bool terminal1IsBoundary = is_boundary(vh_terminal);
					bool terminal2IsBoundary = is_boundary(vh_terminal2);

					auto heh_1_to_T = find_halfedge(vh_terminal, vh_T);
					auto heh_T_to_2 = find_halfedge(vh_T, vh_terminal2);
					OpenMesh::HalfedgeHandle heh_over_1_to_T, heh_below_T_to_2;
					Vec3f v1, v1_1, v1_2, v2, v2_1, v2_2;
					if (terminal1IsBoundary) {
						v1_1 = Vec3f(0.0f);
						v1 = -calc_edge_vector(heh_1_to_T);
					}
					else {
						heh_over_1_to_T = prev_halfedge_handle(opposite_halfedge_handle(prev_halfedge_handle(heh_1_to_T)));
						v1_1 = -calc_edge_vector(heh_1_to_T);
						v1_2 = -calc_edge_vector(heh_over_1_to_T);
						v1 = v1_1 + v1_2;
					}
					if (terminal2IsBoundary) {
						v2_2 = Vec3f(0.0f);
						v2 = calc_edge_vector(heh_T_to_2);
					}
					else {
						heh_below_T_to_2 = next_halfedge_handle(opposite_halfedge_handle(next_halfedge_handle(heh_T_to_2)));
						v2_1 = calc_edge_vector(heh_T_to_2);
						v2_2 = calc_edge_vector(heh_below_T_to_2);
						v2 = v2_1 + v2_2;
					}
					auto u1 = calc_edge_vector(next_halfedge_handle(opposite_halfedge_handle(heh_1_to_T)));
					auto u2 = calc_edge_vector(next_halfedge_handle(heh_1_to_T));

					// Swap around if needed
					/*if (smaller(u2, u1)) {
						std::swap(u1, u2);
					}
					if (smaller(v2, v1)) {
						std::swap(v1, v2);
					}
					if (smaller(u1, v1)) {
						std::swap(u1, v1); 
						std::swap(u2, v2);
					}*/


					float halfwayU = u1.length() / (u1.length() + u2.length());
					float totalVLength = (v1.length() + v2.length());
					float q0V = 0.0f;
					float q1V = v1_1.length() / totalVLength;
					float q2V = v1.length() / totalVLength;
					float q3V = 1.0f - (v2_2.length() / totalVLength);
					float q4V = 1.0f;

					auto offset = point(vh_T);
					uint32_t controlPointIndex = createLocalSurfaceControlPoints(
						u1 + v1, u2 + v1, u1 + v2, u2 + v2);
					uint32_t controlPointCount = 9;

					std::unordered_map<OpenMesh::FaceHandle, uint32_t> faceMappings_terminal1;
					if (!terminal1IsBoundary) {
						faceMappings_terminal1.insert({ face_handle(opposite_halfedge_handle(heh_over_1_to_T)), 
							addBoundaryInfo(BoundaryInfo(0.0f, halfwayU, q0V, q1V)) });
						faceMappings_terminal1.insert({ face_handle(heh_over_1_to_T), 
							addBoundaryInfo(BoundaryInfo(halfwayU, 1.0f, q0V, q1V)) });
					}
					faceMappings_terminal1.insert({ face_handle(opposite_halfedge_handle(heh_1_to_T)), 
						addBoundaryInfo(BoundaryInfo(0.0f, halfwayU, q1V, q3V)) });
					faceMappings_terminal1.insert({ face_handle(heh_1_to_T), 
						addBoundaryInfo(BoundaryInfo(halfwayU, 1.0f, q1V, q2V)) });

					std::unordered_map<OpenMesh::FaceHandle, uint32_t> faceMappings_T;
					faceMappings_T.insert({ face_handle(heh_1_to_T), 
						addBoundaryInfo(BoundaryInfo(halfwayU, 1.0f, q1V, q2V)) });
					faceMappings_T.insert({ face_handle(heh_T_to_2), 
						addBoundaryInfo(BoundaryInfo(halfwayU, 1.0f, q2V, q3V)) });

					std::unordered_map<OpenMesh::FaceHandle, uint32_t> faceMappings_terminal2;
					faceMappings_terminal2.insert({ face_handle(opposite_halfedge_handle(heh_T_to_2)), 
						addBoundaryInfo(BoundaryInfo(0.0f, halfwayU, q1V, q3V)) });
					faceMappings_terminal2.insert({ face_handle(heh_T_to_2), 
						addBoundaryInfo(BoundaryInfo(halfwayU, 1.0f, q2V, q3V)) });
					if (!terminal2IsBoundary) {
						faceMappings_terminal2.insert({ face_handle(opposite_halfedge_handle(heh_below_T_to_2)), 
							addBoundaryInfo(BoundaryInfo(0.0f, halfwayU, q3V, q4V)) });
						faceMappings_terminal2.insert({ face_handle(heh_below_T_to_2), 
							addBoundaryInfo(BoundaryInfo(halfwayU, 1.0f, q3V, q4V)) });
					}

					addLocus(vh_terminal, controlPointIndex, controlPointCount, faceMappings_terminal1, offset, true);
					addLocus(vh_terminal2, controlPointIndex, controlPointCount, faceMappings_terminal2, offset, false);
					addLocus(vh_T, controlPointIndex, controlPointCount, faceMappings_T, offset, false);
				}
			}
		}
	}

	void Lattice::setupLocalSurfacesAndPatches()
	{
		// Add local surfaces
		for (auto f_itr = faces_begin(); f_itr != faces_end(); f_itr++)
		{
			// Face
			auto fh = *f_itr;
			auto f = face(fh);

			// Halfedge/edge handles
			auto topheh = halfedge_handle(fh);

			auto leftheh = next_halfedge_handle(topheh);

			auto bottomheh = next_halfedge_handle(leftheh);

			auto rightheh = prev_halfedge_handle(topheh);

			// Add local surfaces on vertices
			auto vh1 = to_vertex_handle(topheh); // Top Left
			auto vh2 = to_vertex_handle(rightheh); // Top Right
			auto vh3 = to_vertex_handle(leftheh); // Bottom left
			auto vh4 = to_vertex_handle(bottomheh); // Bottom right

			Vec3f L2RDir = -calc_edge_vector(topheh).normalize();
			Vec3f T2BDir = calc_edge_vector(leftheh).normalize();

			// Local surfaces are handled differently based on where it is in the patch
			if (property(LatticeProperties::LocusIndex, vh1) == -1)
				addLocusOnVertex(vh1, vh3, vh2, 1, L2RDir, T2BDir);
			if (property(LatticeProperties::LocusIndex, vh2) == -1)
				addLocusOnVertex(vh2, vh1, vh4, 2, L2RDir, T2BDir);
			if (property(LatticeProperties::LocusIndex, vh3) == -1)
				addLocusOnVertex(vh3, vh4, vh1, 3, L2RDir, T2BDir);
			if (property(LatticeProperties::LocusIndex, vh4) == -1)
				addLocusOnVertex(vh4, vh2, vh3, 4, L2RDir, T2BDir);

			// Setup patch vertices for face.
			Patch patch;
			patch.lociIndices = {
				property(LatticeProperties::LocusIndex, vh1),
				property(LatticeProperties::LocusIndex, vh2),
				property(LatticeProperties::LocusIndex, vh3),
				property(LatticeProperties::LocusIndex, vh4)
			};
			patch.fh = fh;
			if (m_useRandomPatchColors) {
				patch.color = glm::vec3(m_rng.random(0.0f, 1.0f), m_rng.random(0.0f, 1.0f), m_rng.random(0.0f, 1.0f));
				/*float fac = (float)m_numPatches / (float)n_faces();
				if (fac < 0.5f) {
					fac *= 2.0f;
					patch.color = glm::vec3(1.0f - fac, fac, 0.0f);
				}
				else {
					fac = (fac - 0.5f) * 2.0f;
					patch.color = glm::vec3(1.0f - fac, 1.0f - fac, fac);
				}*/
			}
			else {
				patch.color = m_color;
			}
			m_patches.push_back(std::move(patch));
		}
	}

	void Lattice::addLatticeProperties()
	{
		add_property(LatticeProperties::VertexValence);
		add_property(LatticeProperties::LocusIndex);
		add_property(LatticeProperties::Type);
	}

	void Lattice::addLocus(
		OpenMesh::VertexHandle vertex, uint32_t controlPointIndex, uint32_t controlPointCount,
		std::unordered_map<OpenMesh::FaceHandle, uint32_t>& boundaryIndices, Vec3f offset, bool addMatrix)
	{
		if (addMatrix) {
			m_matrices.push_back(
				glm::translate(glm::mat4(1.0f), glm::vec3(offset[0], offset[1], offset[2])));
			m_initialMatrices.push_back(
				glm::translate(glm::mat4(1.0f), glm::vec3(offset[0], offset[1], offset[2])));
		}

		Locus locus;
		locus.controlPointIndex = controlPointIndex;
		locus.controlPointCount = controlPointCount;
		locus.vh = vertex;
		locus.boundaryIndices = boundaryIndices;
		locus.matrixIndex = m_matrices.size() - 1;

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

		property(LatticeProperties::LocusIndex, vertex) = m_loci.size();

		m_loci.push_back(locus);
	}

	size_t Lattice::addBoundaryInfo(BoundaryInfo boundary)
	{
		auto res = m_boundaryMap.insert({ boundary, m_numUniqueBoundaries });
		if (res.second) {
			m_boundaries.push_back(boundary);
			m_numUniqueBoundaries++;
		}

		return res.first->second;
	}

	void Lattice::addLocusOnVertex(
		OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
		OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace,
		Vec3f L2RDir, Vec3f T2BDir)
	{
		auto valence = property(LatticeProperties::VertexValence, vertex);
		switch (valence) {
			case 2: { addLocusOnCornerVertex(vertex, next_vertex, prev_vertex, vertexIndexOnFace, L2RDir, T2BDir); break; }
			case 3: { addLocusOnBoundaryVertex(vertex, next_vertex, prev_vertex, vertexIndexOnFace, L2RDir, T2BDir); break; }
			case 4: { addLocusOnInnerVertex(vertex, next_vertex, prev_vertex, vertexIndexOnFace, L2RDir, T2BDir); break; }
		}
	}

	void Lattice::addLocusOnCornerVertex(
		OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
		OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace,
		Vec3f L2RDir, Vec3f T2BDir)
	{
		addLocalSurfaceOnLoci(vertex, L2RDir, T2BDir);

		//auto heh_locus_to_next = find_halfedge(vertex, next_vertex);
		//Vec3f offset = point(vertex);
		//Vec3f zero(0.0f, 0.0f, 0.0f);
		//Vec3f o_to_n = calc_edge_vector(heh_locus_to_next);
		//Vec3f o_to_n_next = calc_edge_vector(next_halfedge_handle(heh_locus_to_next));
		//Vec3f o_to_p = -calc_edge_vector(prev_halfedge_handle(heh_locus_to_next));
		//Vec3f v_inner = o_to_n + o_to_n_next;

		//auto p1 = zero;
		//auto p2 = o_to_p;
		//auto p3 = o_to_n;
		//auto p4 = v_inner;

		//if (smaller(p2, p1, L2RDir)) std::swap(p1, p2);
		//if (smaller(p4, p3, L2RDir)) std::swap(p3, p4);
		//if (smaller(p3, p1, T2BDir)) {
		//	std::swap(p1, p3);
		//	std::swap(p2, p4);
		//}

		////uint32_t controlPointIndex = 0;
		//uint32_t controlPointCount = 9;
		//uint32_t controlPointIndex = createLocalSurfaceControlPoints(zero, o_to_p, o_to_n, v_inner);
		///*switch (vertexIndexOnFace) {
		//	case 1: { controlPointIndex = createLocalSurfaceControlPoints(zero, o_to_p, o_to_n, v_inner); break; }
		//	case 2: { controlPointIndex = createLocalSurfaceControlPoints(o_to_n, zero, v_inner, o_to_p); break; }
		//	case 3: { controlPointIndex = createLocalSurfaceControlPoints(o_to_p, v_inner, zero, o_to_n); break; }
		//	case 4: { controlPointIndex = createLocalSurfaceControlPoints(v_inner, o_to_n, o_to_p, zero); break; }
		//}*/

		//std::unordered_map<OpenMesh::FaceHandle, uint32_t> boundaryIndices;
		//m_boundaries.push_back(BoundaryInfo( 0.0f, 1.0f, 0.0f, 1.0f ));
		//boundaryIndices.insert({ face_handle(heh_locus_to_next), m_boundaries.size() - 1 });
		//addLocus(vertex, controlPointIndex, controlPointCount, boundaryIndices, offset);
	}

	void Lattice::addLocusOnBoundaryVertex(
		OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
		OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace,
		Vec3f L2RDir, Vec3f T2BDir)
	{
		/*std::vector<OpenMesh::FaceHandle> faces;
		auto heh_next = find_halfedge(vertex, next_vertex);
		OpenMesh::HalfedgeHandle heh_prev;
		if (is_boundary(edge_handle(heh_next))) {
			heh_prev = find_halfedge(vertex, prev_vertex);
		}
		else {
			heh_prev = find_halfedge(prev_vertex, vertex);
		}
		faces.push_back(face_handle(heh_next));
		faces.push_back(face_handle(heh_prev));
		addLocalSurfaceOnBoundaryLocus(vertex, faces, L2RDir, T2BDir);*/

		Vec3f offset = point(vertex);
		Vec3f zero(0.0f, 0.0f, 0.0f);

		auto heh_next = find_halfedge(vertex, next_vertex);

		BoundaryInfo boundaryFace(0.0f, 1.0f, 0.0f, 1.0f);
		BoundaryInfo boundaryAdjFace(0.0f, 1.0f, 0.0f, 1.0f);

		uint32_t controlPointIndex = 0;
		uint32_t controlPointCount = 9;
		std::unordered_map<OpenMesh::FaceHandle, uint32_t> boundaryIndices;

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

			float halfwayPoint = 0.5f + (n.length() / (n - no).length() - 0.5f) / 2;
			//float halfwayPoint = 0.5f;
			switch (vertexIndexOnFace) {
				case 1: { 
					controlPointIndex = createLocalSurfaceControlPoints(no, no + pnop_2, no + pnop, zero, p_2, p, n, n + pnp_2, n + pnp);
					boundaryFace.vs = halfwayPoint;
					boundaryAdjFace.ve = halfwayPoint;
					break; 
				}
				case 2: { 
					controlPointIndex = createLocalSurfaceControlPoints(n, zero, no, n + pnp_2, p_2, no + pnop_2, n + pnp, p, no + pnop);
					boundaryFace.ue = halfwayPoint;
					boundaryAdjFace.us = halfwayPoint;
					break; 
				}
				case 3: { 
					controlPointIndex = createLocalSurfaceControlPoints(no + pnop, p, n + pnp, no + pnop_2, p_2, n + pnp_2, no, zero, n);
					boundaryFace.us = halfwayPoint;
					boundaryAdjFace.ue = halfwayPoint;
					break; 
				}
				case 4: { 
					controlPointIndex = createLocalSurfaceControlPoints(n + pnp, n + pnp_2, n, p, p_2, zero, no + pnop, no + pnop_2, no);
					boundaryFace.ve = halfwayPoint;
					boundaryAdjFace.vs = halfwayPoint;
					break; 
				}
			}

			boundaryIndices.insert({ face_handle(heh_next), addBoundaryInfo(boundaryFace) });
			boundaryIndices.insert({ face_handle(heh_prev), addBoundaryInfo(boundaryAdjFace) });
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

			float halfwayPoint = 0.5f + (p.length() / (p - po).length() - 0.5f) / 2;
			//float halfwayPoint = 0.5f;

			switch (vertexIndexOnFace) {
				case 1: { 
					controlPointIndex = createLocalSurfaceControlPoints(po, zero, p, po + npon_2, n_2, p + npn_2, po + npon, n, p + npn);
					boundaryFace.us = halfwayPoint;
					boundaryAdjFace.ue = halfwayPoint;
					break; 
				}
				case 2: { 
					controlPointIndex = createLocalSurfaceControlPoints(po + npon, po + npon_2, po, n, n_2, zero, p + npn, p + npn_2, p);
					boundaryFace.vs = halfwayPoint;
					boundaryAdjFace.ve = halfwayPoint;
					break; 
				}
				case 3: { 
					controlPointIndex = createLocalSurfaceControlPoints(p, p + npn_2, p + npn, zero, n_2, n, po, po + npon_2, po + npon);
					boundaryFace.ve = halfwayPoint;
					boundaryAdjFace.vs = halfwayPoint;
					break; 
				}
				case 4: { 
					controlPointIndex = createLocalSurfaceControlPoints(p + npn, n, po + npon, p + npn_2, n_2, po + npon_2, p, zero, po);
					boundaryFace.ue = halfwayPoint;
					boundaryAdjFace.us = halfwayPoint;
					break; 
				}
			}

			boundaryIndices.insert({ face_handle(heh_next), addBoundaryInfo(boundaryFace) });
			boundaryIndices.insert({ face_handle(opposite_halfedge_handle(heh_next)), addBoundaryInfo(boundaryAdjFace) });
		}

		addLocus(vertex, controlPointIndex, controlPointCount, boundaryIndices, offset);
	}

	void Lattice::addLocusOnInnerVertex(
		OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
		OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace,
		Vec3f L2RDir, Vec3f T2BDir)
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

		float halfwayPointU = 0.5f + (u4.length() / (u3 - u4).length() - 0.5f) / 2;
		float halfwayPointV = 0.5f + (v2.length() / (v1 - v2).length() - 0.5f) / 2;
		//float halfwayPointU = 0.5f;
		//float halfwayPointV = 0.5f;
		uint32_t controlPointIndex = 0;
		uint32_t controlPointCount = 9;

		switch (vertexIndexOnFace) {
			case 1: { 
				controlPointIndex = createLocalSurfaceControlPoints(u31 + v1, v1, u42 + v1, u3, zero, u4, u35 + v2, v2, u46 + v2);
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
				controlPointIndex = createLocalSurfaceControlPoints(u35 + v2, u3, u31 + v1, v2, zero, v1, u46 + v2, u4, u42 + v1);
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
				controlPointIndex = createLocalSurfaceControlPoints(u42 + v1, u4, u46 + v2, v1, zero, v2, u31 + v1, u3, u35 + v2);
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
				controlPointIndex = createLocalSurfaceControlPoints(u46 + v2, v2, u35 + v2, u4, zero, u3, u42 + v1, v1, u31 + v1);
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

		std::unordered_map<OpenMesh::FaceHandle, uint32_t> boundaryIndices;
		boundaryIndices.insert({ face_handle(hehf1), addBoundaryInfo(boundaryF1) });
		boundaryIndices.insert({ face_handle(hehf2), addBoundaryInfo(boundaryF2) });
		boundaryIndices.insert({ face_handle(hehf3), addBoundaryInfo(boundaryF3) });
		boundaryIndices.insert({ face_handle(hehf4), addBoundaryInfo(boundaryF4) });

		addLocus(vertex, controlPointIndex, controlPointCount, boundaryIndices, offset);
	}

	uint32_t Lattice::createLocalSurfaceControlPoints(
		Vec3f topLeft, Vec3f topRight,
		Vec3f bottomLeft, Vec3f bottomRight)
	{
		Vec3f topMiddle = (topLeft + topRight) / 2;
		Vec3f middleLeft = (topLeft + bottomLeft) / 2;
		Vec3f middleRight = (topRight + bottomRight) / 2;
		Vec3f bottomMiddle = (bottomLeft + bottomRight) / 2;
		Vec3f middle = ((topMiddle + bottomMiddle) / 2 + (middleLeft + middleRight) / 2) / 2;

		uint32_t idx = m_controlPoints.size();
		m_controlPoints.push_back(glm::vec4(topLeft[0], topLeft[1], topLeft[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(topMiddle[0], topMiddle[1], topMiddle[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(topRight[0], topRight[1], topRight[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(middleLeft[0], middleLeft[1], middleLeft[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(middle[0], middle[1], middle[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(middleRight[0], middleRight[1], middleRight[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(bottomLeft[0], bottomLeft[1], bottomLeft[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(bottomMiddle[0], bottomMiddle[1], bottomMiddle[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(bottomRight[0], bottomRight[1], bottomRight[2], 1.0f));
		return idx;
	}

	uint32_t Lattice::createLocalSurfaceControlPoints(
		Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
		Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
		Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight)
	{
		uint32_t idx = m_controlPoints.size();
		m_controlPoints.push_back(glm::vec4(topLeft[0], topLeft[1], topLeft[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(topMiddle[0], topMiddle[1], topMiddle[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(topRight[0], topRight[1], topRight[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(middleLeft[0], middleLeft[1], middleLeft[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(middle[0], middle[1], middle[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(middleRight[0], middleRight[1], middleRight[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(bottomLeft[0], bottomLeft[1], bottomLeft[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(bottomMiddle[0], bottomMiddle[1], bottomMiddle[2], 1.0f));
		m_controlPoints.push_back(glm::vec4(bottomRight[0], bottomRight[1], bottomRight[2], 1.0f));
		return idx;
	}

	void Lattice::addLocalSurfaceOnLoci(OpenMesh::VertexHandle vh, Vec3f L2RDir, Vec3f T2BDir)
	{
		std::vector<OpenMesh::FaceHandle> faces;
		for (auto adj_face_itr = vf_iter(vh); adj_face_itr.is_valid(); adj_face_itr++) {
			faces.push_back(*adj_face_itr);
		}

		switch (faces.size()) {
		case 1: { addLocalSurfaceOnCornerLocus(vh, faces[0]); break; }
		case 2: { addLocalSurfaceOnBoundaryLocus(vh, faces, L2RDir, T2BDir); break; }
		case 4: { addLocalSurfaceOnInnerLocus(vh, faces, L2RDir, T2BDir); break; }
		}
	}

	void Lattice::getCornerPointsOfFaceL2RT2B(OpenMesh::FaceHandle& fh, std::vector<Vec3f>& points)
	{
		// Halfedge/edge handles
		OpenMesh::HalfedgeHandle topheh = halfedge_handle(fh);
		OpenMesh::HalfedgeHandle leftheh = next_halfedge_handle(topheh);
		OpenMesh::HalfedgeHandle bottomheh = next_halfedge_handle(leftheh);
		OpenMesh::HalfedgeHandle rightheh = prev_halfedge_handle(topheh);

		if (points.capacity() != 4) points.resize(4);

		points[0] = point(to_vertex_handle(topheh));	// Top Left
		points[1] = point(to_vertex_handle(rightheh));	// Top Right
		points[2] = point(to_vertex_handle(leftheh));	// Bottom left
		points[3] = point(to_vertex_handle(bottomheh));	// Bottom right
	}

	void Lattice::addLocalSurfaceOnCornerLocus(OpenMesh::VertexHandle vh, OpenMesh::FaceHandle fh)
	{
		std::vector<Vec3f> points(4);
		getCornerPointsOfFaceL2RT2B(fh, points);
		Vec3f offset = point(vh);

		//uint32_t controlPointIndex = 0;
		uint32_t controlPointCount = 9;
		uint32_t controlPointIndex = createLocalSurfaceControlPoints(points[0] - offset, 
			points[1] - offset, points[2] - offset, points[3] - offset);

		std::unordered_map<OpenMesh::FaceHandle, uint32_t> boundaryIndices;
		boundaryIndices.insert({ fh, addBoundaryInfo(BoundaryInfo(0.0f, 1.0f, 0.0f, 1.0f)) });
		addLocus(vh, controlPointIndex, controlPointCount, boundaryIndices, offset);
	}
	void Lattice::addLocalSurfaceOnBoundaryLocus(
		OpenMesh::VertexHandle& vh, std::vector<OpenMesh::FaceHandle>& faces, Vec3f L2RDir, Vec3f T2BDir)
	{
		Vec3f f1TopLeft = point(to_vertex_handle(halfedge_handle(faces[0])));
		Vec3f f2TopLeft = point(to_vertex_handle(halfedge_handle(faces[1])));
		float l2rDiff = ((f1TopLeft - f2TopLeft) * L2RDir).sqrnorm();
		float t2bDiff = ((f1TopLeft - f2TopLeft) * T2BDir).sqrnorm();
		Vec3f offset = point(vh);

		std::vector<Vec3f> f1_points(4), f2_points(4);

		BoundaryInfo boundaryF1(0.0f, 1.0f, 0.0f, 1.0f);
		BoundaryInfo boundaryF2(0.0f, 1.0f, 0.0f, 1.0f);

		std::unordered_map<OpenMesh::FaceHandle, uint32_t> boundaryIndices;

		Vec3f p00, p10, p20, p01, p11, p21, p02, p12, p22;
		if (l2rDiff > t2bDiff) {
			// Faces are horizontal
			if (smaller(faces[1], faces[0], L2RDir)) std::swap(faces[0], faces[1]);
			getCornerPointsOfFaceL2RT2B(faces[0], f1_points);
			getCornerPointsOfFaceL2RT2B(faces[1], f2_points);
			p00 = f1_points[0]; p10 = f1_points[1]; p20 = f2_points[1]; // Top row
			p02 = f1_points[2]; p12 = f1_points[3]; p22 = f2_points[3]; // Bottom row
			p01 = p00 + (p02 - p00) * 0.5f; p11 = p10 + (p12 - p10) * 0.5f; p21 = p20 + (p22 - p20) * 0.5f;
			float halfwayU = p10.length() / (p20 - p00).length();
			boundaryF1.ue = halfwayU; boundaryF2.us = halfwayU;
		}
		else {
			// Faces are vertical
			if (smaller(faces[1], faces[0], T2BDir)) std::swap(faces[0], faces[1]);
			getCornerPointsOfFaceL2RT2B(faces[0], f1_points);
			getCornerPointsOfFaceL2RT2B(faces[1], f2_points);
			p00 = f1_points[0]; p01 = f1_points[2]; p02 = f2_points[2]; // Left col
			p20 = f1_points[1]; p21 = f1_points[3]; p22 = f2_points[3]; // Right col
			p10 = p00 + (p20 - p00) * 0.5f; p11 = p01 + (p21 - p01) * 0.5f; p12 = p02 + (p22 - p02) * 0.5f;
			float halfwayV = p11.length() / (p21 - p01).length();
			boundaryF1.ve = halfwayV; boundaryF2.vs = halfwayV;
		}

		boundaryIndices.insert({ faces[0], addBoundaryInfo(boundaryF1) });
		boundaryIndices.insert({ faces[1], addBoundaryInfo(boundaryF2) });

		uint32_t controlPointCount = 9;
		uint32_t controlPointIndex = createLocalSurfaceControlPoints(
			p00 - offset, p10 - offset, p20 - offset,
			p01 - offset, p11 - offset, p21 - offset,
			p02 - offset, p12 - offset, p22 - offset);

		addLocus(vh, controlPointIndex, controlPointCount, boundaryIndices, offset);
	}
}