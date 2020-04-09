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
		vkDestroyPipeline(*m_device, m_localSurfacePipeline, m_allocator);
		vkDestroyPipeline(*m_device, m_localSurfaceWireframePipeline, m_allocator);
		vkDestroyPipeline(*m_device, m_patchPipeline, m_allocator);
		vkDestroyPipeline(*m_device, m_patchWireframePipeline, m_allocator);
		vkDestroyPipeline(*m_device, m_normalPipeline, m_allocator);

		vkDestroyPipelineLayout(*m_device, m_pipelineLayout, m_allocator);
		vkDestroyDescriptorSetLayout(*m_device, m_descriptorSetLayout, m_allocator);

		if (m_pointsBuffer.buffer != VK_NULL_HANDLE) {
			vkDestroyBuffer(*m_device, m_pointsBuffer.buffer, m_allocator);
			vkFreeMemory(*m_device, m_pointsBuffer.memory, m_allocator);
		}
		if (m_linesBuffer.buffer != VK_NULL_HANDLE) {
			vkDestroyBuffer(*m_device, m_linesBuffer.buffer, m_allocator);
			vkFreeMemory(*m_device, m_linesBuffer.memory, m_allocator);
		}
		if (m_localSurfaceVertexBuffer.buffer != VK_NULL_HANDLE) {
			vkDestroyBuffer(*m_device, m_localSurfaceVertexBuffer.buffer, m_allocator);
			vkFreeMemory(*m_device, m_localSurfaceVertexBuffer.memory, m_allocator);
		}
		if (m_patchVertexBuffer.buffer != VK_NULL_HANDLE) {
			vkDestroyBuffer(*m_device, m_patchVertexBuffer.buffer, m_allocator);
			vkFreeMemory(*m_device, m_patchVertexBuffer.memory, m_allocator);
		}
	}

	void Lattice::addToCommandbuffer(VkCommandBuffer& commandBuffer)
	{
		if (!m_draw) return;

		VkDeviceSize offsets[1] = { 0 };

		vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &m_descriptorSet, 0, NULL);

		if (m_drawSurface)
		{
			if (m_wireframe)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_patchWireframePipeline);
			else
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_patchPipeline);
			vkCmdBindVertexBuffers(commandBuffer, 0, 1, &m_patchVertexBuffer.buffer, offsets);
			vkCmdDraw(commandBuffer, m_patchVertexBuffer.count, 1, 0, 0);
			
			if (m_drawNormals)
			{
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_normalPipeline);
				vkCmdDraw(commandBuffer, m_patchVertexBuffer.count, 1, 0, 0);
			}
		}

		// Draw local surfaces
		if (m_drawLocalSurfaces)
		{
			if (m_wireframe)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_localSurfaceWireframePipeline);
			else
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_localSurfacePipeline);
			vkCmdBindVertexBuffers(commandBuffer, 0, 1, &m_localSurfaceVertexBuffer.buffer, offsets);
			vkCmdDraw(commandBuffer, m_localSurfaceVertexBuffer.count, 1, 0, 0);
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
		m_latticeUniforms.projection = projection;
		m_latticeUniforms.modelview = view; // Model matrix for this?
		updateLatticeUniformBuffer();
	}

	void Lattice::update(double dt)
	{
		if (m_animate) {
			for (auto& mat : m_matrixUniforms.matrices)
				mat = glm::rotate(mat, glm::radians(60.0f) * (float)dt, glm::vec3(0.0f, 0.0f, 1.0f));
			updateMatrixUniformBuffer();
		}
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
				overlay->checkBox("Animate", &m_animate);
				if (!m_drawPixelAccurate)
				{
					if (overlay->sliderInt("TessInner", &m_latticeUniforms.tessInner, 0, 64)) updateLatticeUniformBuffer();
					if (overlay->sliderInt("TessOuter", &m_latticeUniforms.tessOuter, 0, 64)) updateLatticeUniformBuffer();
				}
				if (overlay->comboBox("B-Function", &m_latticeUniforms.bFunctionIndex, BFunctionNames)) updateLatticeUniformBuffer();
			}

			if (overlay->header("Surfaces", false))
			{
				for (size_t i = 0; i < m_patches.size(); i++)
				{
					auto& patch = m_patches[i];
					std::string patchTitle = "Patch " + std::to_string(i);
					if (overlay->header(patchTitle.c_str(), false))
					{
						for (size_t j = 0; j < 4; j++)
						{
							std::string lociTitle = "Local " + std::to_string(j);
							if (overlay->header(lociTitle.c_str(), false))
							{
								auto idx = m_loci[patch.lociIndices[j]].matrixIndex;
								auto& mat = m_matrixUniforms.matrices[idx];
								float* x = &mat[3][0];
								//auto* trans = &mat[3];
								//auto* x = &trans[0];
								if (overlay->inputVec3("Trans", x, 1.0f, 1))
								{
									updateMatrixUniformBuffer();
								}
							}
						}
					}
				}
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
			auto valenceProp = property(LatticeProperties::VertexValence, vh);

			valenceProp = valence(vh);

			property(LatticeProperties::VertexValence, vh) = valenceProp;
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
		// Set up all values
		m_matrixUniforms.matrices.resize(0);
		m_patchUniforms.controlPoints.resize(0);
		m_patchUniforms.boundaries.resize(0);

		m_numControlPoints = 0;
		m_numLoci = 0;
		m_numPatches = 0;

		m_loci.resize(0);
		m_patches.resize(0);

		m_localSurfaceVertices.resize(0);
		m_patchVertices.resize(0);

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

			auto& loci00 = m_loci[property(LatticeProperties::LocusIndex, vh1)];
			auto& loci10 = m_loci[property(LatticeProperties::LocusIndex, vh2)];
			auto& loci01 = m_loci[property(LatticeProperties::LocusIndex, vh3)];
			auto& loci11 = m_loci[property(LatticeProperties::LocusIndex, vh4)];

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

			LocalSurfaceVertex localVert00;
			localVert00.controlPointIndex = loci00.controlPointIndex;
			localVert00.controlPointCount = loci00.controlPointCount;
			localVert00.matrixIndex = loci00.matrixIndex;
			localVert00.boundaryIndex = m_patchUniforms.boundaries.size();
			m_patchUniforms.boundaries.push_back(loci00.faceMappings[fh]);
			m_patchVertices.push_back(localVert00);

			LocalSurfaceVertex localVert10;
			localVert10.controlPointIndex = loci10.controlPointIndex;
			localVert10.controlPointCount = loci10.controlPointCount;
			localVert10.matrixIndex = loci10.matrixIndex;
			localVert10.boundaryIndex = m_patchUniforms.boundaries.size();
			m_patchUniforms.boundaries.push_back(loci10.faceMappings[fh]);
			m_patchVertices.push_back(localVert10);

			LocalSurfaceVertex localVert01;
			localVert01.controlPointIndex = loci01.controlPointIndex;
			localVert01.controlPointCount = loci01.controlPointCount;
			localVert01.matrixIndex = loci01.matrixIndex;
			localVert01.boundaryIndex = m_patchUniforms.boundaries.size();
			m_patchUniforms.boundaries.push_back(loci01.faceMappings[fh]);
			m_patchVertices.push_back(localVert01);

			LocalSurfaceVertex localVert11;
			localVert11.controlPointIndex = loci11.controlPointIndex;
			localVert11.controlPointCount = loci11.controlPointCount;
			localVert11.matrixIndex = loci11.matrixIndex;
			localVert11.boundaryIndex = m_patchUniforms.boundaries.size();
			m_patchUniforms.boundaries.push_back(loci11.faceMappings[fh]);
			m_patchVertices.push_back(localVert11);

			m_numPatches++;
		}

		std::cout << "Done" << std::endl;
	}

	void Lattice::addLatticeProperties()
	{
		add_property(LatticeProperties::VertexValence);
		add_property(LatticeProperties::LocusIndex);
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

		property(LatticeProperties::LocusIndex, vertex) = m_numLoci;

		m_loci.push_back(locus);
		m_numLoci++;

		for (auto& p : controlPoints) {
			m_patchUniforms.controlPoints.push_back(glm::vec4(p[0], p[1], p[2], 0.0f));
		}
		//m_patchUniforms.controlPoints.insert(m_patchUniforms.controlPoints.end(), controlPoints.begin(), controlPoints.end());
		m_numControlPoints += controlPoints.size();

		m_matrixUniforms.matrices.push_back(
			glm::translate(glm::mat4(1.0f), glm::vec3(offset[0], offset[1], offset[2])));

		m_localSurfaceVertices.push_back(LocalSurfaceVertex(
			locus.controlPointIndex, locus.controlPointCount, locus.matrixIndex, 0));
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
		Vec3f o_to_p = -calc_edge_vector(prev_halfedge_handle(heh_locus_to_next));
		Vec3f v_inner = o_to_n + o_to_p;

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

			float halfwayPoint = no.length() / (n - no).length();
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

			float halfwayPoint = p.length() / (p - po).length();

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

		float halfwayPointU = u3.length() / (u3 - u4).length();
		float halfwayPointV = v2.length() / (v1 - v2).length();

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

		uint32_t gridPointVertexBufferSize = gridPoints.size() * sizeof(GridVertex);
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

		uint32_t gridLineVertexBufferSize = gridLines.size() * sizeof(GridVertex);
		createDeviceLocalBuffer(m_linesBuffer.buffer, m_linesBuffer.memory,
			static_cast<void*>(gridLines.data()), gridLineVertexBufferSize, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
		m_linesBuffer.count = gridLines.size();

		// Local surface vertices
		uint32_t localSurfaceVertexBufferSize = m_localSurfaceVertices.size() * sizeof(LocalSurfaceVertex);
		createDeviceLocalBuffer(m_localSurfaceVertexBuffer.buffer, m_localSurfaceVertexBuffer.memory,
			static_cast<void*>(m_localSurfaceVertices.data()), localSurfaceVertexBufferSize, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
		m_localSurfaceVertexBuffer.count = m_localSurfaceVertices.size();

		// Patch vertices
		uint32_t patchVertexBufferSize = m_patchVertices.size() * sizeof(LocalSurfaceVertex);
		createDeviceLocalBuffer(m_patchVertexBuffer.buffer, m_patchVertexBuffer.memory,
			static_cast<void*>(m_patchVertices.data()), patchVertexBufferSize, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
		m_patchVertexBuffer.count = m_patchVertices.size();
	}

	void Lattice::prepareUniformBuffers()
	{
		// Shared tessellation shader stages uniform buffer
		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&m_latticeUniformBuffer,
			sizeof(m_latticeUniforms)));

		// Map persistent
		VK_CHECK_RESULT(m_latticeUniformBuffer.map());

		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&m_matrixUniformBuffer,
			sizeof(glm::mat4) * m_numLoci
		));

		VK_CHECK_RESULT(m_matrixUniformBuffer.map());

		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&m_patchUniformBuffer,
			sizeof(glm::vec4) * m_numControlPoints + sizeof(BoundaryInfo) * m_numPatches * 4
		));

		VK_CHECK_RESULT(m_patchUniformBuffer.map());

		updateLatticeUniformBuffer();
		updateMatrixUniformBuffer();
		updatePatchUniformBuffer();
	}

	void Lattice::setupDescriptorSetLayouts()
	{
		std::vector<VkDescriptorSetLayout> descriptorLayouts;

		// LAttice
		std::vector<VkDescriptorSetLayoutBinding> setLayoutBindings = {
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT | 
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_GEOMETRY_BIT,
				0), 
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
				1), 
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
				2)
		};

		VkDescriptorSetLayoutCreateInfo descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(
			setLayoutBindings.data(), static_cast<uint32_t>(setLayoutBindings.size()));
		VK_CHECK_RESULT(vkCreateDescriptorSetLayout(*m_device, &descriptorLayout,
			nullptr, &m_descriptorSetLayout));
		descriptorLayouts.push_back(m_descriptorSetLayout);

		// Pipeline Layout
		VkPipelineLayoutCreateInfo pipelineLayoutCreateInfo =
			vks::initializers::pipelineLayoutCreateInfo(descriptorLayouts.data(), descriptorLayouts.size());

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
				VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT,
				VK_TRUE);
		blendAttachmentState.srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
		blendAttachmentState.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
		blendAttachmentState.colorBlendOp = VK_BLEND_OP_ADD;
		blendAttachmentState.srcAlphaBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
		blendAttachmentState.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
		blendAttachmentState.alphaBlendOp = VK_BLEND_OP_SUBTRACT;

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
		shaderStages[0] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/LatticeHelpers/poscolorpass/poscolorpass.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		shaderStages[1] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/LatticeHelpers/poscolorpass/poscolorpass.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);

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
		std::vector<VkVertexInputBindingDescription> localVertexBindings = LocalSurfaceVertex::GetBindingDescriptions();
		std::vector<VkVertexInputAttributeDescription> localInputAttributes = LocalSurfaceVertex::GetAttributeDescriptions();
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
			int numLocalSurfacesControlPoints;
			int numLocalSurfaces;
			int numPatches;
		} specializationData;
		specializationData.numLocalSurfacesControlPoints = m_numControlPoints;
		specializationData.numLocalSurfaces = m_numLoci;
		specializationData.numPatches = m_numPatches;
		std::array<VkSpecializationMapEntry, 3> specializationMapEntries;
		specializationMapEntries[0].constantID = 0;
		specializationMapEntries[0].size = sizeof(specializationData.numLocalSurfacesControlPoints);
		specializationMapEntries[0].offset = offsetof(SpecializationData, numLocalSurfacesControlPoints);
		specializationMapEntries[1].constantID = 1;
		specializationMapEntries[1].size = sizeof(specializationData.numLocalSurfaces);
		specializationMapEntries[1].offset = offsetof(SpecializationData, numLocalSurfaces);
		specializationMapEntries[2].constantID = 2;
		specializationMapEntries[2].size = sizeof(specializationData.numPatches);
		specializationMapEntries[2].offset = offsetof(SpecializationData, numPatches);

		VkSpecializationInfo specializationInfo = {};
		specializationInfo.dataSize = sizeof(specializationData);
		specializationInfo.mapEntryCount = static_cast<uint32_t>(specializationMapEntries.size());
		specializationInfo.pMapEntries = specializationMapEntries.data();
		specializationInfo.pData = &specializationData;

		localShaderStages[0] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/LatticeHelpers/localsurfaces/bezier3x3.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		localShaderStages[1] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/LatticeHelpers/localsurfaces/bezier3x3.tesc.spv", VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		localShaderStages[2] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/LatticeHelpers/localsurfaces/bezier3x3.tese.spv", VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		localShaderStages[2].pSpecializationInfo = &specializationInfo;
		localShaderStages[3] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/LatticeHelpers/localsurfaces/bezier3x3.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(localShaderStages.size());
		pipelineCreateInfo.pStages = localShaderStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_localSurfacePipeline));

		rasterizationState.polygonMode = VK_POLYGON_MODE_LINE;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_localSurfaceWireframePipeline));

		// Patch pipelines
		VkPipelineTessellationStateCreateInfo bsTessState =
			vks::initializers::pipelineTessellationStateCreateInfo(4);
		pipelineCreateInfo.pTessellationState = &bsTessState;

		std::array<VkPipelineShaderStageCreateInfo, 4> bsShaderStages;
		bsShaderStages[0] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/lattice.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		bsShaderStages[1] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/lattice.tesc.spv", VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		bsShaderStages[2] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/lattice.tese.spv", VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		bsShaderStages[2].pSpecializationInfo = &specializationInfo;
		bsShaderStages[3] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/lattice.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(bsShaderStages.size());
		pipelineCreateInfo.pStages = bsShaderStages.data();

		rasterizationState.polygonMode = VK_POLYGON_MODE_LINE;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, nullptr, &m_patchWireframePipeline));

		rasterizationState.polygonMode = VK_POLYGON_MODE_FILL;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, nullptr, &m_patchPipeline));

		std::array<VkPipelineShaderStageCreateInfo, 5> normalShaderStages;
		normalShaderStages[0] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/lattice.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		normalShaderStages[1] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/lattice.tesc.spv", VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		normalShaderStages[2] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/lattice.tese.spv", VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		normalShaderStages[2].pSpecializationInfo = &specializationInfo;
		normalShaderStages[3] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/normals.geom.spv", VK_SHADER_STAGE_GEOMETRY_BIT);
		normalShaderStages[4] = loadShader("P:/Projects/Visual Studio/OMMeshlessAnimationFramework/shaders/Lattice/normals.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(normalShaderStages.size());
		pipelineCreateInfo.pStages = normalShaderStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_normalPipeline));
	}

	void Lattice::setupDescriptorPool()
	{
		std::vector<VkDescriptorPoolSize> poolSizes =
		{
			vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 3)
		};

		VkDescriptorPoolCreateInfo descriptorPoolInfo =
			vks::initializers::descriptorPoolCreateInfo(
				static_cast<uint32_t>(poolSizes.size()),
				poolSizes.data(),
				3);

		VK_CHECK_RESULT(vkCreateDescriptorPool(*m_device, &descriptorPoolInfo, nullptr, m_descriptorPool));
	}

	void Lattice::setupDescriptorSets()
	{
		VkDescriptorSetAllocateInfo allocInfo;
		std::vector<VkWriteDescriptorSet> writeDescriptorSets;

		allocInfo = vks::initializers::descriptorSetAllocateInfo(*m_descriptorPool, &m_descriptorSetLayout, 1);
		VK_CHECK_RESULT(vkAllocateDescriptorSets(*m_device, &allocInfo, &m_descriptorSet));

		writeDescriptorSets =
		{
			vks::initializers::writeDescriptorSet(
				m_descriptorSet,
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				0,
				&m_latticeUniformBuffer.descriptor),
			vks::initializers::writeDescriptorSet(
				m_descriptorSet,
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				1,
				&m_matrixUniformBuffer.descriptor),
			vks::initializers::writeDescriptorSet(
				m_descriptorSet,
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				2,
				&m_patchUniformBuffer.descriptor),
		};

		vkUpdateDescriptorSets(*m_device, static_cast<uint32_t>(writeDescriptorSets.size()), writeDescriptorSets.data(), 0, NULL);
	}

	void Lattice::updateLatticeUniformBuffer()
	{
		memcpy(m_latticeUniformBuffer.mapped, &m_latticeUniforms, sizeof(m_latticeUniforms));
	}

	void Lattice::updateMatrixUniformBuffer()
	{
		memcpy(m_matrixUniformBuffer.mapped, &m_matrixUniforms.matrices[0][0], sizeof(glm::mat4) * m_numLoci);
	}

	void Lattice::updatePatchUniformBuffer()
	{
		std::vector<glm::vec4> patchUniforms;
		patchUniforms.insert(patchUniforms.end(), m_patchUniforms.controlPoints.begin(), m_patchUniforms.controlPoints.end());
		for (auto& boundary : m_patchUniforms.boundaries)
			patchUniforms.push_back(glm::vec4(boundary.us, boundary.ue, boundary.vs, boundary.ve));
		memcpy(m_patchUniformBuffer.mapped, &patchUniforms[0],
			sizeof(glm::vec4) * m_numControlPoints + sizeof(BoundaryInfo) * m_numPatches * 4);
	}
}