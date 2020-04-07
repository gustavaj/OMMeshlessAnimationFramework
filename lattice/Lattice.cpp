#include "Lattice.h"

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
			if (p == topLeft) tl = *it;
			if (p == topRight) tr = *it;
			if (p == bottomLeft) bl = *it;
			if (p == bottomRight) br = *it;
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
			coords[i] = std::vector<Vec3f>(cols);
			for (size_t j = 0; j < cols; j++) {
				float theta = j * d_angle;
				coords[i][j] = Vec3f(radius * std::cos(theta), dh * i, radius * std::sin(theta)) + center;
			}
		}

		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols - 1; j++) {
				addPatch(coords[i + 1][j], coords[i + 1][j + 1], coords[i][j], coords[i][j + 1]);
				if (j == (cols - 2))
					addPatch(coords[i + 1][j + 1], coords[i + 1][0], coords[i][j + 1], coords[i][0]);
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

	void Lattice::initVulkanStuff(
		VkDevice* device, vks::VulkanDevice* vulkanDevice,
		VkQueue* queue, VkCommandPool* commandPool,
		VkDescriptorPool* descriptorPool, VkRenderPass* renderPass,
		VkAllocationCallbacks* allocator)
	{
		m_device = device;
		m_vulkanDevice = vulkanDevice;
		m_queue = queue;
		m_commandPool = commandPool;
		m_descriptorPool = descriptorPool;
		m_renderPass = renderPass;
		m_allocator = allocator;

		createBuffers();
		prepareUniformBuffers();
		setupDescriptorSetLayouts();
		preparePipelines();
		setupDescriptorPool();
		setupDescriptorSets();
	}

	void Lattice::destroyVulkanStuff()
	{
		vkDestroyPipeline(*m_device, m_pointsPipeline, m_allocator);
		vkDestroyPipeline(*m_device, m_linesPipeline, m_allocator);

		vkDestroyPipelineLayout(*m_device, m_pipelineLayout, m_allocator);
		vkDestroyDescriptorSetLayout(*m_device, m_descriptorSetLayout, m_allocator);
	}

	void Lattice::addToCommandbuffer(VkCommandBuffer& commandBuffer)
	{
		if (!m_draw) return;

		VkDeviceSize offsets[1] = { 0 };

		vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &m_descriptorSet, 0, NULL);

		// Draw local surfaces
		if (m_drawLocalSurfaces)
		{
			if (m_wireframe)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_localSurfaceWireframePipeline);
			else
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_localSurfacePipeline);
			vkCmdBindVertexBuffers(commandBuffer, 0, 1, &m_localSurfaceBuffer.buffer, offsets);
			vkCmdDraw(commandBuffer, m_localSurfaceBuffer.count, 1, 0, 0);
		}

		// Draw Lattice Grid
		if (m_drawLatticeGrid)
		{
			// Draw gridlines
			vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_linesPipeline);
			vkCmdBindVertexBuffers(commandBuffer, 0, 1, &m_linesBuffer.buffer, offsets);
			vkCmdDraw(commandBuffer, m_linesBuffer.count, 1, 0, 0);

			// Draw gridpoints
			vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pointsPipeline);
			vkCmdBindVertexBuffers(commandBuffer, 0, 1, &m_pointsBuffer.buffer, offsets);
			vkCmdDraw(commandBuffer, m_pointsBuffer.count, 1, 0, 0);
		}
	}

	void Lattice::onViewChanged(glm::mat4 projection, glm::mat4 view)
	{
		m_ubo.projection = projection;
		m_ubo.modelview = view; // Model matrix for this?
		updateUniformBuffer();
	}

	bool Lattice::onUpdateUIOverlay(vks::UIOverlay* overlay)
	{
		bool rebuildCmd = false;

		if (overlay->header(m_name.c_str(), false))
		{
			if (overlay->checkBox("Render", &m_draw)) rebuildCmd = true;
			if (overlay->header("Settings", false))
			{
				if (overlay->checkBox("Draw Surface", &m_drawSurface)) rebuildCmd = true;
				if (overlay->checkBox("Draw Local", &m_drawLocalSurfaces)) rebuildCmd = true;
				if (overlay->checkBox("Draw Grid", &m_drawLatticeGrid)) rebuildCmd = true;
				if (overlay->checkBox("Draw Wireframe", &m_wireframe)) rebuildCmd = true;
				if (overlay->checkBox("Draw Normals", &m_drawNormals)) rebuildCmd = true;
				if (overlay->checkBox("Pixel-Accurate", &m_drawPixelAccurate)) rebuildCmd = true;
				if (!m_drawPixelAccurate)
				{
					if (overlay->sliderInt("TessInner", &m_ubo.tessInner, 0, 64)) updateUniformBuffer();
					if (overlay->sliderInt("TessOuter", &m_ubo.tessOuter, 0, 64)) updateUniformBuffer();
				}
				if (overlay->comboBox("B-Function", &m_ubo.bFunctionIndex, BFunctionNames)) updateUniformBuffer();
			}
		}

		// TODO: Local Surface Editing

		return rebuildCmd;
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
			auto valenceProp = property(LatticeProperties::LocusValence, vh);

			valenceProp = valence(vh);

			property(LatticeProperties::LocusValence, vh) = valenceProp;
			set_color(vh, LOCUS_VALENCE_COLOR[valenceProp]);

			// Initialize local surface idx to -1
			property(LatticeProperties::LocusLocalSurfaceIdx, vh) = -1;
		}
	}

	void Lattice::handleTLoci()
	{
		// Look for T-loci, and handle it!
		for (auto v_itr = vertices_begin(); v_itr != vertices_end(); v_itr++)
		{
			auto vh = *v_itr;
			if (property(LatticeProperties::LocusValence, vh) == 5
				|| (property(LatticeProperties::LocusValence, vh) == 4 && is_boundary(vh)))
			{
				// Find the opposite extreme point vertex
				OpenMesh::VertexHandle vh2, vht;
				for (auto vv_itr = vv_iter(vh); vv_itr.is_valid(); vv_itr++)
				{
					if (property(LatticeProperties::LocusValence, *vv_itr) == 5 ||
						(property(LatticeProperties::LocusValence, *vv_itr) == 4 && is_boundary(*vv_itr)))
					{
						vh2 = *vv_itr;
						break;
					}
				}

				auto heh12 = find_halfedge(vh, vh2);

				OpenMesh::HalfedgeHandle heht1, heht2;

				// Find the other t-point
				for (auto vv_itr = vv_iter(vh); vv_itr.is_valid(); vv_itr++)
				{
					heht1 = find_halfedge(*vv_itr, vh);
					heht2 = find_halfedge(*vv_itr, vh2);
					if (heht1.is_valid() && heht2.is_valid())
					{
						vht = *vv_itr;
						break;
					}
				}

				if (vht.is_valid()) {
					// Change color of vertices
					// Yellow for extreme points, purple for t-poitns
					set_color(vh, Col3(255, 255, 0));
					set_color(vh2, Col3(255, 255, 0));
					set_color(vht, Col3(255, 0, 255));

					// heh12
					auto heh21 = opposite_halfedge_handle(heh12);

					auto vh3 = from_vertex_handle(prev_halfedge_handle(heh21));
					auto vh4 = to_vertex_handle(next_halfedge_handle(heh21));

					delete_edge(edge_handle(heh12));

					garbage_collection(true, true, true);

					std::vector<OpenMesh::VertexHandle> vhs{ vh4, vh3, vh2, vht, vh };
					add_face(vhs);

					// Update valence
					property(LatticeProperties::LocusValence, vh) -= 1;
					property(LatticeProperties::LocusValence, vh2) -= 1;


					auto heh_1_to_t = find_halfedge(vh, vht);
					auto heh_t_to_2 = find_halfedge(vht, vh2);

					// Create local surface and update indices for all affected parts.
					auto vh_before = from_vertex_handle(prev_halfedge_handle(opposite_halfedge_handle(prev_halfedge_handle(heh_1_to_t))));
					auto vh_after = to_vertex_handle(next_halfedge_handle(opposite_halfedge_handle(next_halfedge_handle(heh_t_to_2))));

					auto v1 = calc_edge_vector(opposite_halfedge_handle(heh_1_to_t)) + calc_edge_vector(find_halfedge(vh, vh_before));
					auto v2 = calc_edge_vector(heh_t_to_2) + calc_edge_vector(find_halfedge(vh2, vh_after));

					auto u1 = calc_edge_vector(next_halfedge_handle(opposite_halfedge_handle(heh_1_to_t)));
					auto u2 = calc_edge_vector(opposite_halfedge_handle(prev_halfedge_handle(heh_1_to_t)));
					auto u3 = calc_edge_vector(opposite_halfedge_handle(prev_halfedge_handle(opposite_halfedge_handle(heh_t_to_2))));
					auto u4 = calc_edge_vector(next_halfedge_handle(heh_t_to_2));

					addLocalSurface(v1 + u1, v1 + u2, v2 + u3, v2 + u4, point(vht));
					int ls_idx = m_localSurfaces.size() - 1;
					property(LatticeProperties::LocusLocalSurfaceIdx, vh) = ls_idx;
					property(LatticeProperties::LocusLocalSurfaceIdx, vht) = ls_idx;
					property(LatticeProperties::LocusLocalSurfaceIdx, vh2) = ls_idx;
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
			if (property(LatticeProperties::LocusLocalSurfaceIdx, vh1) == -1)
				addLocalSurfaceOnLocus(vh1, vh3, vh2, 1);
			if (property(LatticeProperties::LocusLocalSurfaceIdx, vh2) == -1)
				addLocalSurfaceOnLocus(vh2, vh1, vh4, 2);
			if (property(LatticeProperties::LocusLocalSurfaceIdx, vh3) == -1)
				addLocalSurfaceOnLocus(vh3, vh4, vh1, 3);
			if (property(LatticeProperties::LocusLocalSurfaceIdx, vh4) == -1)
				addLocalSurfaceOnLocus(vh4, vh2, vh3, 4);

			// Add BlendingSpline with the local surfaces of the vertices, edges and faces
			// Create vector of indices for the blending spline surface
			//std::vector<uint32_t> indices;

			//indices.push_back(property(LatticeProperties::LocusLocalSurfaceIdx, vh1));
			//indices.push_back(property(LatticeProperties::LocusLocalSurfaceIdx, vh2));
			//indices.push_back(property(LatticeProperties::LocusLocalSurfaceIdx, vh3));
			//indices.push_back(property(LatticeProperties::LocusLocalSurfaceIdx, vh4));

			//bool isBoundaryTop = is_boundary(topeh);
			//bool isBoundaryLeft = is_boundary(lefteh);
			//bool isBoundaryBot = is_boundary(boteh);
			//bool isBoundaryRight = is_boundary(righteh);

			//// If both vertices of a boundary edge has a valence of four, the edge should not be treated as a boundary edge.
			//if (property(LatticeProperties::LocusValence, from_vertex_handle(halfedge_handle(topeh, 0))) == 4
			//	&& property(LatticeProperties::LocusValence, to_vertex_handle(halfedge_handle(topeh, 0))) == 4) {
			//	isBoundaryTop = false;
			//}
			//if (property(LatticeProperties::LocusValence, from_vertex_handle(halfedge_handle(lefteh, 0))) == 4
			//	&& property(LatticeProperties::LocusValence, to_vertex_handle(halfedge_handle(lefteh, 0))) == 4) {
			//	isBoundaryLeft = false;
			//}
			//if (property(LatticeProperties::LocusValence, from_vertex_handle(halfedge_handle(boteh, 0))) == 4
			//	&& property(LatticeProperties::LocusValence, to_vertex_handle(halfedge_handle(boteh, 0))) == 4) {
			//	isBoundaryBot = false;
			//}
			//if (property(LatticeProperties::LocusValence, from_vertex_handle(halfedge_handle(righteh, 0))) == 4
			//	&& property(LatticeProperties::LocusValence, to_vertex_handle(halfedge_handle(righteh, 0))) == 4) {
			//	isBoundaryRight = false;
			//}

			/*UVBoundary boundary(
				isBoundaryTop,
				isBoundaryLeft,
				isBoundaryBot,
				isBoundaryRight,
				idx);

			m_blendingSurfaces.push_back(BlendingSurface(indices, boundary));
			idx += 1.0f;*/
		}
	}

	void Lattice::addLatticeProperties()
	{
		add_property(LatticeProperties::LocusValence);
		add_property(LatticeProperties::LocusLocalSurfaceIdx);
		add_property(LatticeProperties::FaceMappings);
	}

	void Lattice::addLocalSurfaceOnLocus(
		OpenMesh::VertexHandle locus, OpenMesh::VertexHandle next_locus,
		OpenMesh::VertexHandle prev_locus, int locusIndexOnFace)
	{
		auto valence = property(LatticeProperties::LocusValence, locus);
		switch (valence) {
			case 2: { addLocalSurfaceOnCornerLocus(locus, next_locus, prev_locus, locusIndexOnFace); break; }
			case 3: { addLocalSurfaceOnEdgeLocus(locus, next_locus, prev_locus, locusIndexOnFace); break; }
			case 4: { addLocalSurfaceOnInnerLocus(locus, next_locus, prev_locus, locusIndexOnFace); break; }
		}
	}

	void Lattice::addLocalSurfaceOnCornerLocus(
		OpenMesh::VertexHandle locus, OpenMesh::VertexHandle next_locus,
		OpenMesh::VertexHandle prev_locus, int locusIndexOnFace)
	{
		auto heh_locus_to_next = find_halfedge(locus, next_locus);
		Vec3f offset = point(locus);
		Vec3f zero(0.0f, 0.0f, 0.0f);
		Vec3f o_to_n = calc_edge_vector(heh_locus_to_next);
		Vec3f o_to_p = -calc_edge_vector(prev_halfedge_handle(heh_locus_to_next));
		Vec3f v_inner = o_to_n + o_to_p;

		switch (locusIndexOnFace) {
			case 1: { addLocalSurface(zero, o_to_p, o_to_n, v_inner, offset); break; }
			case 2: { addLocalSurface(o_to_n, zero, v_inner, o_to_p, offset); break; }
			case 3: { addLocalSurface(o_to_p, v_inner, zero, o_to_n, offset); break; }
			case 4: { addLocalSurface(v_inner, o_to_n, o_to_p, zero, offset); break; }
		}

		property(LatticeProperties::LocusLocalSurfaceIdx, locus) = m_localSurfaces.size() - 1;
		property(LatticeProperties::FaceMappings, locus).insert({face_handle(heh_locus_to_next), BoundaryInfo(0.0f, 1.0f, 0.0f, 1.0f)});
	}

	void Lattice::addLocalSurfaceOnEdgeLocus(
		OpenMesh::VertexHandle locus, OpenMesh::VertexHandle next_locus,
		OpenMesh::VertexHandle prev_locus, int locusIndexOnFace)
	{
		Vec3f center = point(locus);
		Vec3f zero(0.0f, 0.0f, 0.0f);

		auto heh_next = find_halfedge(locus, next_locus);

		BoundaryInfo boundaryFace(0.0f, 1.0f, 0.0f, 1.0f);
		BoundaryInfo boundaryAdjFace(0.0f, 1.0f, 0.0f, 1.0f);

		if (is_boundary(edge_handle(heh_next))) {
			auto heh_prev = find_halfedge(locus, prev_locus);
			Vec3f n = calc_edge_vector(heh_next); // vector from locus to next_locus
			Vec3f no = -calc_edge_vector(prev_halfedge_handle(heh_prev)); // vector from locus to the locus on the adjacent face
			Vec3f nop = -calc_edge_vector(next_halfedge_handle(next_halfedge_handle(heh_prev)));
			Vec3f p = calc_edge_vector(heh_prev);
			Vec3f np = calc_edge_vector(next_halfedge_handle(heh_next));

			Vec3f p_2 = p * 0.5f;
			Vec3f pnop = lerp(p, nop, 1.0f);
			Vec3f pnop_2 = pnop * 0.5f;
			Vec3f pnp = lerp(p, np, 1.0f);
			Vec3f pnp_2 = pnp * 0.5f;

			float halfwayPoint = (n - no).length() / no.length();

			switch (locusIndexOnFace) {
				case 1: { 
					addLocalSurface(no, no + pnop_2, no + pnop, zero, p_2, p, n, n + pnp_2, n + pnp, center); 
					boundaryFace.vs = halfwayPoint;
					boundaryAdjFace.ve = halfwayPoint;
					break; 
				}
				case 2: { 
					addLocalSurface(n, zero, no, n + pnp_2, p_2, no + pnop_2, n + pnp, p, no + pnop, center); 
					boundaryFace.ue = halfwayPoint;
					boundaryAdjFace.us = halfwayPoint;
					break; 
				}
				case 3: { 
					addLocalSurface(no + pnop, p, n + pnp, no + pnop_2, p_2, n + pnp_2, no, zero, n, center); 
					boundaryFace.us = halfwayPoint;
					boundaryAdjFace.ue = halfwayPoint;
					break; 
				}
				case 4: { 
					addLocalSurface(n + pnp, n + pnp_2, n, p, p_2, zero, no + pnop, no + pnop_2, no, center); 
					boundaryFace.ve = halfwayPoint;
					boundaryAdjFace.vs = halfwayPoint;
					break; 
				}
			}

			property(LatticeProperties::FaceMappings, locus).insert({ face_handle(heh_next),
				boundaryFace });
			property(LatticeProperties::FaceMappings, locus).insert({ face_handle(heh_prev),
				boundaryAdjFace });
		}
		else {
			auto heh_prev = find_halfedge(prev_locus, locus);
			Vec3f p = -calc_edge_vector(heh_prev);
			Vec3f po = calc_edge_vector(next_halfedge_handle(opposite_halfedge_handle(heh_next)));
			Vec3f pon = calc_edge_vector(next_halfedge_handle(next_halfedge_handle(opposite_halfedge_handle(heh_next))));
			Vec3f n = calc_edge_vector(heh_next);
			Vec3f pn = -calc_edge_vector(prev_halfedge_handle(heh_prev));

			Vec3f n_2 = n * 0.5f;
			Vec3f npn = lerp(n, pn, 1.0f);
			Vec3f npn_2 = npn * 0.5f;
			Vec3f npon = lerp(n, pon, 1.0f);
			Vec3f npon_2 = npon * 0.5f;

			float halfwayPoint = (p - po).length() / p.length();

			switch (locusIndexOnFace) {
				case 1: { 
					addLocalSurface(po, zero, p, po + npon_2, n_2, p + npn_2, po + npon, n, p + npn, center);
					boundaryFace.us = halfwayPoint;
					boundaryAdjFace.ue = halfwayPoint;
					break; 
				}
				case 2: { 
					addLocalSurface(po + npon, po + npon_2, po, n, n_2, zero, p + npn, p + npn_2, p, center); 
					boundaryFace.vs = halfwayPoint;
					boundaryAdjFace.ve = halfwayPoint;
					break; 
				}
				case 3: { 
					addLocalSurface(p, p + npn_2, p + npn, zero, n_2, n, po, po + npon_2, po + npon, center); 
					boundaryFace.ve = halfwayPoint;
					boundaryAdjFace.vs = halfwayPoint;
					break; 
				}
				case 4: { 
					addLocalSurface(p + npn, n, po + npon, p + npn_2, n_2, po + npon_2, p, zero, po, center); 
					boundaryFace.ue = halfwayPoint;
					boundaryAdjFace.us = halfwayPoint;
					break; 
				}
			}

			property(LatticeProperties::FaceMappings, locus).insert({ face_handle(heh_next),
				boundaryFace });
			property(LatticeProperties::FaceMappings, locus).insert({ face_handle(opposite_halfedge_handle(heh_next)),
				boundaryAdjFace });
		}

		property(LatticeProperties::LocusLocalSurfaceIdx, locus) = m_localSurfaces.size() - 1;
	}

	void Lattice::addLocalSurfaceOnInnerLocus(
		OpenMesh::VertexHandle locus, OpenMesh::VertexHandle next_locus,
		OpenMesh::VertexHandle prev_locus, int locusIndexOnFace)
	{
		Vec3f po = point(locus);
		Vec3f zero(0.0f, 0.0f, 0.0f);

		auto hehf1 = find_halfedge(locus, next_locus);
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

		auto u42 = lerp(u4, u2, 1.0f);
		auto u31 = lerp(u3, u1, 1.0f);
		auto u46 = lerp(u4, u6, 1.0f);
		auto u35 = lerp(u3, u5, 1.0f);

		float halfwayPointU = (u3 - u4).length() / u3.length();
		float halfwayPointV = (v1 - v2).length() / v2.length();

		switch (locusIndexOnFace) {
			case 1: { 
				addLocalSurface(u31 + v1, v1, u42 + v1, u3, zero, u4, u35 + v2, v2, u46 + v2, po);
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
				addLocalSurface(u35 + v2, u3, u31 + v1, v2, zero, v1, u46 + v2, u4, u42 + v1, po);
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
				addLocalSurface(u42 + v1, u4, u46 + v2, v1, zero, v2, u31 + v1, u3, u35 + v2, po);
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
				addLocalSurface(u46 + v2, v2, u35 + v2, u4, zero, u3, u42 + v1, v1, u31 + v1, po);
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

		property(LatticeProperties::FaceMappings, locus).insert({ face_handle(hehf1), boundaryF1 });
		property(LatticeProperties::FaceMappings, locus).insert({ face_handle(hehf2), boundaryF2 });
		property(LatticeProperties::FaceMappings, locus).insert({ face_handle(hehf3), boundaryF3 });
		property(LatticeProperties::FaceMappings, locus).insert({ face_handle(hehf4), boundaryF4 });

		property(LatticeProperties::LocusLocalSurfaceIdx, locus) = m_localSurfaces.size() - 1;
	}

	void Lattice::addLocalSurface(
		Vec3f topLeft, Vec3f topRight,
		Vec3f bottomLeft, Vec3f bottomRight, Vec3f offset)
	{
		m_localSurfaces.push_back(
			LocalSurface(topLeft, topRight, bottomLeft, bottomRight, offset, Vec2f(0.5, 0.5))
		);
	}

	void Lattice::addLocalSurface(
		Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
		Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
		Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight, Vec3f offset)
	{
		m_localSurfaces.push_back(
			LocalSurface(
				topLeft, topMiddle, topRight,
				middleLeft, middle, middleRight,
				bottomLeft, bottomMiddle, bottomRight, offset,
				Vec2f((topLeft - topMiddle).length() / (topLeft - topRight).length(), (topLeft - middleLeft).length() / (topLeft - bottomLeft).length()))
		);
	}

	void Lattice::createDeviceLocalBuffer(
		VkBuffer& buffer, VkDeviceMemory& memory, void* data,
		uint32_t bufferSize, VkBufferUsageFlagBits usage)
	{
		if (buffer != VK_NULL_HANDLE) {
			vkDestroyBuffer(*m_device, buffer, m_allocator);
		}
		if (memory != VK_NULL_HANDLE) {
			vkFreeMemory(*m_device, memory, m_allocator);
		}

		struct {
			VkBuffer buffer;
			VkDeviceMemory memory;
		} stagingBuffer;

		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			bufferSize,
			&stagingBuffer.buffer,
			&stagingBuffer.memory,
			data));

		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			usage | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
			VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
			bufferSize,
			&buffer,
			&memory));

		VkCommandBuffer cmdBuffer;

		VkCommandBufferAllocateInfo cmdBufAllocateInfo =
			vks::initializers::commandBufferAllocateInfo(
				*m_commandPool,
				VK_COMMAND_BUFFER_LEVEL_PRIMARY,
				1);

		VK_CHECK_RESULT(vkAllocateCommandBuffers(*m_device, &cmdBufAllocateInfo, &cmdBuffer));

		VkCommandBufferBeginInfo cmdBufInfo = vks::initializers::commandBufferBeginInfo();
		VK_CHECK_RESULT(vkBeginCommandBuffer(cmdBuffer, &cmdBufInfo));

		VkBufferCopy copyRegion = {};

		copyRegion.size = bufferSize;
		vkCmdCopyBuffer(
			cmdBuffer,
			stagingBuffer.buffer,
			buffer,
			1,
			&copyRegion);

		VK_CHECK_RESULT(vkEndCommandBuffer(cmdBuffer));

		VkSubmitInfo submitInfo = {};
		submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
		submitInfo.commandBufferCount = 1;
		submitInfo.pCommandBuffers = &cmdBuffer;

		VK_CHECK_RESULT(vkQueueSubmit(*m_queue, 1, &submitInfo, VK_NULL_HANDLE));
		VK_CHECK_RESULT(vkQueueWaitIdle(*m_queue));

		vkFreeCommandBuffers(*m_device, *m_commandPool, 1, &cmdBuffer);

		vkDestroyBuffer(*m_device, stagingBuffer.buffer, nullptr);
		vkFreeMemory(*m_device, stagingBuffer.memory, nullptr);
	}

	void Lattice::createBuffers()
	{
		// Lattice grid points
		std::vector<GridVertex> gridPoints;

		for (auto v_itr = vertices_begin(); v_itr != vertices_end(); v_itr++)
		{
			GridVertex v;
			v.pos = point(*v_itr);
			v.col = color(*v_itr);
			gridPoints.push_back(v);
		}

		uint32_t gridPointVertexBufferSize = gridPoints.size() * sizeof(OML::GridVertex);
		createDeviceLocalBuffer(m_pointsBuffer.buffer, m_pointsBuffer.memory,
			static_cast<void*>(gridPoints.data()), gridPointVertexBufferSize, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
		m_pointsBuffer.count = gridPoints.size();

		// Lattice grid lines
		std::vector<GridVertex> gridLines;

		for (auto e_itr = edges_begin(); e_itr != edges_end(); e_itr++)
		{
			GridVertex v1;
			v1.pos = point(from_vertex_handle(halfedge_handle(*e_itr, 0)));
			v1.col = color(*e_itr);
			gridLines.push_back(v1);

			GridVertex v2;
			v2.pos = point(to_vertex_handle(halfedge_handle(*e_itr, 0)));
			v2.col = color(*e_itr);
			gridLines.push_back(v2);
		}

		uint32_t gridLineVertexBufferSize = gridLines.size() * sizeof(OML::GridVertex);
		createDeviceLocalBuffer(m_linesBuffer.buffer, m_linesBuffer.memory,
			static_cast<void*>(gridLines.data()), gridLineVertexBufferSize, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
		m_linesBuffer.count = gridLines.size();

		// Local surfaces
		m_localSurfaceBuffer.count = m_localSurfaces.size();
		std::vector<Bezier3x3Vertex> localSurfaceVertices(m_localSurfaceBuffer.count);
		for (size_t i = 0; i < m_localSurfaceBuffer.count; i++)
		{
			localSurfaceVertices[i] = m_localSurfaces[i].vertex();
		}
		uint32_t localSurfaceVertexBufferSize = m_localSurfaceBuffer.count * sizeof(Bezier3x3Vertex);
		createDeviceLocalBuffer(m_localSurfaceBuffer.buffer, m_localSurfaceBuffer.memory,
			static_cast<void*>(localSurfaceVertices.data()), localSurfaceVertexBufferSize, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
	}

	void Lattice::prepareUniformBuffers()
	{
		// Shared tessellation shader stages uniform buffer
		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&m_uniformBuffer,
			sizeof(m_ubo)));

		// Map persistent
		VK_CHECK_RESULT(m_uniformBuffer.map());

		updateUniformBuffer();
	}

	void Lattice::setupDescriptorSetLayouts()
	{
		VkDescriptorSetLayoutCreateInfo descriptorLayout;
		VkPipelineLayoutCreateInfo pipelineLayoutCreateInfo;
		std::vector<VkDescriptorSetLayoutBinding> setLayoutBindings;

		// Terrain
		setLayoutBindings =
		{
			// Binding 0 : Shared shader ubo
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				VK_SHADER_STAGE_ALL_GRAPHICS,
				0)
		};

		descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(setLayoutBindings.data(),
			static_cast<uint32_t>(setLayoutBindings.size()));
		VK_CHECK_RESULT(vkCreateDescriptorSetLayout(*m_device, &descriptorLayout,
			nullptr, &m_descriptorSetLayout));
		pipelineLayoutCreateInfo =
			vks::initializers::pipelineLayoutCreateInfo(&m_descriptorSetLayout, 1);

		// Push constants
		/*VkPushConstantRange pushConstantRange = vks::initializers::pushConstantRange(
			VK_SHADER_STAGE_ALL_GRAPHICS, sizeof(OML::UVBoundary), 0);
		pipelineLayoutCreateInfo.pushConstantRangeCount = 1;
		pipelineLayoutCreateInfo.pPushConstantRanges = &pushConstantRange;*/

		VK_CHECK_RESULT(vkCreatePipelineLayout(*m_device, &pipelineLayoutCreateInfo, m_allocator, &m_pipelineLayout));
	}

	VkPipelineShaderStageCreateInfo Lattice::loadShader(std::string fileName, VkShaderStageFlagBits stage)
	{
		VkPipelineShaderStageCreateInfo shaderStage = {};
		shaderStage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
		shaderStage.stage = stage; shaderStage.module = vks::tools::loadShader(fileName.c_str(), *m_device);
		shaderStage.pName = "main";
		assert(shaderStage.module != VK_NULL_HANDLE);
		return shaderStage;
	}

	void Lattice::preparePipelines()
	{
		// Input Assembly States
		VkPipelineInputAssemblyStateCreateInfo lineInputAssemblyState =
			vks::initializers::pipelineInputAssemblyStateCreateInfo(
				VK_PRIMITIVE_TOPOLOGY_LINE_LIST,
				0,
				VK_FALSE);
		VkPipelineInputAssemblyStateCreateInfo pointInputAssemblyState =
			vks::initializers::pipelineInputAssemblyStateCreateInfo(
				VK_PRIMITIVE_TOPOLOGY_POINT_LIST,
				0,
				VK_FALSE);
		// Local surface pipeline
		VkPipelineInputAssemblyStateCreateInfo patchInputAssemblyState =
			vks::initializers::pipelineInputAssemblyStateCreateInfo(
				VK_PRIMITIVE_TOPOLOGY_PATCH_LIST,
				0,
				VK_FALSE);

		VkPipelineRasterizationStateCreateInfo rasterizationState =
			vks::initializers::pipelineRasterizationStateCreateInfo(
				VK_POLYGON_MODE_FILL,
				VK_CULL_MODE_NONE,
				VK_FRONT_FACE_COUNTER_CLOCKWISE,
				0);

		VkPipelineColorBlendAttachmentState blendAttachmentState =
			vks::initializers::pipelineColorBlendAttachmentState(
				0xf,
				VK_FALSE);

		VkPipelineColorBlendStateCreateInfo colorBlendState =
			vks::initializers::pipelineColorBlendStateCreateInfo(
				1,
				&blendAttachmentState);

		VkPipelineDepthStencilStateCreateInfo depthStencilState =
			vks::initializers::pipelineDepthStencilStateCreateInfo(
				VK_TRUE,
				VK_TRUE,
				VK_COMPARE_OP_LESS_OR_EQUAL);

		VkPipelineViewportStateCreateInfo viewportState =
			vks::initializers::pipelineViewportStateCreateInfo(1, 1, 0);

		VkPipelineMultisampleStateCreateInfo multisampleState =
			vks::initializers::pipelineMultisampleStateCreateInfo(
				VK_SAMPLE_COUNT_1_BIT,
				0);

		std::vector<VkDynamicState> dynamicStateEnables = {
			VK_DYNAMIC_STATE_VIEWPORT,
			VK_DYNAMIC_STATE_SCISSOR,
			VK_DYNAMIC_STATE_LINE_WIDTH
		};

		VkPipelineDynamicStateCreateInfo dynamicState =
			vks::initializers::pipelineDynamicStateCreateInfo(
				dynamicStateEnables.data(),
				static_cast<uint32_t>(dynamicStateEnables.size()),
				0);

		// Vertex bindings and attributes
		// Binding description
		std::vector<VkVertexInputBindingDescription> vertexInputBindings = {
			OML::GridVertex::GetBindingDescription()
		};

		// Attribute descriptions
		std::vector<VkVertexInputAttributeDescription> vertexInputAttributes =
			OML::GridVertex::GetAttributeDesctiptions();

		VkPipelineVertexInputStateCreateInfo vertexInputState = vks::initializers::pipelineVertexInputStateCreateInfo();
		vertexInputState.vertexBindingDescriptionCount = static_cast<uint32_t>(vertexInputBindings.size());
		vertexInputState.pVertexBindingDescriptions = vertexInputBindings.data();
		vertexInputState.vertexAttributeDescriptionCount = static_cast<uint32_t>(vertexInputAttributes.size());
		vertexInputState.pVertexAttributeDescriptions = vertexInputAttributes.data();

		std::array<VkPipelineShaderStageCreateInfo, 2> shaderStages;
		shaderStages[0] = loadShader("./shaders/LatticeHelpers/poscolorpass/poscolorpass.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		shaderStages[1] = loadShader("./shaders/LatticeHelpers/poscolorpass/poscolorpass.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);

		VkGraphicsPipelineCreateInfo pipelineCreateInfo =
			vks::initializers::pipelineCreateInfo(m_pipelineLayout, *m_renderPass, 0);

		pipelineCreateInfo.pVertexInputState = &vertexInputState;
		pipelineCreateInfo.pRasterizationState = &rasterizationState;
		pipelineCreateInfo.pColorBlendState = &colorBlendState;
		pipelineCreateInfo.pMultisampleState = &multisampleState;
		pipelineCreateInfo.pViewportState = &viewportState;
		pipelineCreateInfo.pDepthStencilState = &depthStencilState;
		pipelineCreateInfo.pDynamicState = &dynamicState;
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(shaderStages.size());
		pipelineCreateInfo.pStages = shaderStages.data();
		pipelineCreateInfo.renderPass = *m_renderPass;

		// Grid lines pipeline
		pipelineCreateInfo.pInputAssemblyState = &lineInputAssemblyState;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_linesPipeline));

		// Grid point pipeline
		pipelineCreateInfo.pInputAssemblyState = &pointInputAssemblyState;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_pointsPipeline));

		
		pipelineCreateInfo.pInputAssemblyState = &patchInputAssemblyState;

		// Local Surface Pipelines
		std::vector<VkVertexInputBindingDescription> localVertexBindings = Bezier3x3Vertex::GetBindingDescription();
		std::vector<VkVertexInputAttributeDescription> localInputAttributes = Bezier3x3Vertex::GetAttributeDescriptions();
		VkPipelineVertexInputStateCreateInfo localInputState = vks::initializers::pipelineVertexInputStateCreateInfo();
		localInputState.vertexBindingDescriptionCount = static_cast<uint32_t>(localVertexBindings.size());
		localInputState.pVertexBindingDescriptions = localVertexBindings.data();
		localInputState.vertexAttributeDescriptionCount = static_cast<uint32_t>(localInputAttributes.size());
		localInputState.pVertexAttributeDescriptions = localInputAttributes.data();
		pipelineCreateInfo.pVertexInputState = &localInputState;

		VkPipelineTessellationStateCreateInfo tessellationState =
			vks::initializers::pipelineTessellationStateCreateInfo(1);
		pipelineCreateInfo.pTessellationState = &tessellationState;

		std::array<VkPipelineShaderStageCreateInfo, 4> localShaderStages;
		struct SpecializationData {
			int numLocalSurfaces;
		} specializationData;
		specializationData.numLocalSurfaces = m_localSurfaces.size();
		std::array<VkSpecializationMapEntry, 1> specializationMapEntries;
		specializationMapEntries[0].constantID = 0;
		specializationMapEntries[0].size = sizeof(specializationData.numLocalSurfaces);
		specializationMapEntries[0].offset = offsetof(SpecializationData, numLocalSurfaces);

		VkSpecializationInfo specializationInfo = {};
		specializationInfo.dataSize = sizeof(specializationData);
		specializationInfo.mapEntryCount = static_cast<uint32_t>(specializationMapEntries.size());
		specializationInfo.pMapEntries = specializationMapEntries.data();
		specializationInfo.pData = &specializationData;

		localShaderStages[0] = loadShader("./shaders/LatticeHelpers/localsurfaces/bezier3x3.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		localShaderStages[1] = loadShader("./shaders/LatticeHelpers/localsurfaces/bezier3x3.tesc.spv", VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		localShaderStages[2] = loadShader("./shaders/LatticeHelpers/localsurfaces/bezier3x3.tese.spv", VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		localShaderStages[2].pSpecializationInfo = &specializationInfo;
		localShaderStages[3] = loadShader("./shaders/LatticeHelpers/localsurfaces/bezier3x3.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(localShaderStages.size());
		pipelineCreateInfo.pStages = localShaderStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_localSurfacePipeline));

		rasterizationState.polygonMode = VK_POLYGON_MODE_LINE;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_localSurfaceWireframePipeline));

		/*
		// Blending surfaces
		VkPipelineTessellationStateCreateInfo bsTessState =
			vks::initializers::pipelineTessellationStateCreateInfo(numLocalSurfacesPerDir * numLocalSurfacesPerDir);
		pipelineCreateInfo.pTessellationState = &bsTessState;

		std::array<VkPipelineShaderStageCreateInfo, 4> bsShaderStages;
		struct BSSpecData {
			int numLocalSurfacesPerDir;
			int totalPatches;
		} bsSpecData;
		bsSpecData.numLocalSurfacesPerDir = numLocalSurfacesPerDir;
		bsSpecData.totalPatches = lattice.numPatches();
		std::array<VkSpecializationMapEntry, 2> bsSpecMapEntries;
		bsSpecMapEntries[0].constantID = 0;
		bsSpecMapEntries[0].size = sizeof(bsSpecData.numLocalSurfacesPerDir);
		bsSpecMapEntries[0].offset = offsetof(BSSpecData, numLocalSurfacesPerDir);
		bsSpecMapEntries[1].constantID = 1;
		bsSpecMapEntries[1].size = sizeof(bsSpecData.totalPatches);
		bsSpecMapEntries[1].offset = offsetof(BSSpecData, totalPatches);

		VkSpecializationInfo bsSpecInfo = {};
		bsSpecInfo.dataSize = sizeof(bsSpecData);
		bsSpecInfo.mapEntryCount = static_cast<uint32_t>(bsSpecMapEntries.size());
		bsSpecInfo.pMapEntries = bsSpecMapEntries.data();
		bsSpecInfo.pData = &bsSpecData;

		bsShaderStages[0] = loadShader("./shaders/OMLattice/lattice.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		bsShaderStages[1] = loadShader("./shaders/OMLattice/lattice.tesc.spv", VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		bsShaderStages[2] = loadShader("./shaders/OMLattice/lattice.tese.spv", VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		bsShaderStages[2].pSpecializationInfo = &bsSpecInfo;
		bsShaderStages[3] = loadShader("./shaders/OMLattice/lattice.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(bsShaderStages.size());
		pipelineCreateInfo.pStages = bsShaderStages.data();

		rasterizationState.polygonMode = VK_POLYGON_MODE_LINE;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(device, pipelineCache, 1, &pipelineCreateInfo, nullptr, &blendingSurfaceWireframePipeline));

		rasterizationState.polygonMode = VK_POLYGON_MODE_FILL;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(device, pipelineCache, 1, &pipelineCreateInfo, nullptr, &blendingSurfacePipeline));

		std::array<VkPipelineShaderStageCreateInfo, 5> normalShaderStages;
		normalShaderStages[0] = loadShader("./shaders/OMLattice/lattice.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		normalShaderStages[1] = loadShader("./shaders/OMLattice/lattice.tesc.spv", VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		normalShaderStages[2] = loadShader("./shaders/OMLattice/lattice.tese.spv", VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		normalShaderStages[2].pSpecializationInfo = &bsSpecInfo;
		normalShaderStages[3] = loadShader("./shaders/OMLattice/normals.geom.spv", VK_SHADER_STAGE_GEOMETRY_BIT);
		normalShaderStages[4] = loadShader("./shaders/OMLattice/normals.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(normalShaderStages.size());
		pipelineCreateInfo.pStages = normalShaderStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(device, pipelineCache, 1, &pipelineCreateInfo, nullptr, &normalPipeline));
		*/
	}

	void Lattice::setupDescriptorPool()
	{
		std::vector<VkDescriptorPoolSize> poolSizes =
		{
			vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 3),
			vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 3)
		};

		VkDescriptorPoolCreateInfo descriptorPoolInfo =
			vks::initializers::descriptorPoolCreateInfo(
				static_cast<uint32_t>(poolSizes.size()),
				poolSizes.data(),
				2);

		VK_CHECK_RESULT(vkCreateDescriptorPool(*m_device, &descriptorPoolInfo, nullptr, m_descriptorPool));
	}

	void Lattice::setupDescriptorSets()
	{
		VkDescriptorSetAllocateInfo allocInfo;
		std::vector<VkWriteDescriptorSet> writeDescriptorSets;

		// Terrain
		allocInfo = vks::initializers::descriptorSetAllocateInfo(*m_descriptorPool, &m_descriptorSetLayout, 1);
		VK_CHECK_RESULT(vkAllocateDescriptorSets(*m_device, &allocInfo, &m_descriptorSet));

		writeDescriptorSets =
		{
			// Binding 0 : Shared shader ubo
			vks::initializers::writeDescriptorSet(
				m_descriptorSet,
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				0,
				&m_uniformBuffer.descriptor)
		};
		vkUpdateDescriptorSets(*m_device, static_cast<uint32_t>(writeDescriptorSets.size()), writeDescriptorSets.data(), 0, NULL);
	}

	void Lattice::updateUniformBuffer()
	{
		memcpy(m_uniformBuffer.mapped, &m_ubo, sizeof(m_ubo));
	}
}