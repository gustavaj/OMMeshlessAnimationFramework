#include "SWVulkanLattice.h"

namespace SWVL
{
	SWVulkanLattice::SWVulkanLattice()
		: SWVulkanLattice("")
	{
	}

	SWVulkanLattice::SWVulkanLattice(std::string name)
		: OML::Lattice(name), m_pointsPipeline(VK_NULL_HANDLE), m_linesPipeline(VK_NULL_HANDLE), m_localSurfacePipeline(VK_NULL_HANDLE),
		  m_localSurfaceWireframePipeline(VK_NULL_HANDLE), m_patchPipeline(VK_NULL_HANDLE), m_normalPipeline(VK_NULL_HANDLE),
		  m_pipelineLayout(VK_NULL_HANDLE), m_descriptorSetLayout(VK_NULL_HANDLE), m_descriptorSet(VK_NULL_HANDLE),
		  m_device(nullptr), m_vulkanDevice(nullptr), m_descriptorPool(nullptr), m_renderPass(nullptr),
		  m_queue(nullptr), m_commandPool(nullptr), m_allocator(nullptr), m_selectedSurface(0)
	{
	}

	SWVulkanLattice::~SWVulkanLattice()
	{
	}

	void SWVulkanLattice::initVulkanStuff(VkDevice* device, vks::VulkanDevice* vulkanDevice, VkQueue* queue, VkCommandPool* commandPool, VkDescriptorPool* descriptorPool, VkRenderPass* renderPass, VkAllocationCallbacks* allocator)
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

		for (size_t i = 0; i < m_numPatches; i++)
		{
			m_listItems.push_back("Patch " + std::to_string(i) + " - p00");
			m_listItems.push_back("Patch " + std::to_string(i) + " - p10");
			m_listItems.push_back("Patch " + std::to_string(i) + " - p01");
			m_listItems.push_back("Patch " + std::to_string(i) + " - p11");
		}
	}

	void SWVulkanLattice::destroyVulkanStuff()
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

	void SWVulkanLattice::addToCommandbuffer(VkCommandBuffer& commandBuffer)
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

	void SWVulkanLattice::onViewChanged(glm::mat4 projection, glm::mat4 view)
	{
		m_uniforms.projection = projection;
		m_uniforms.modelview = view;
		updateLatticeUniformBuffer();
	}

	void SWVulkanLattice::update(double dt)
	{
		if (m_animate) {
			for (auto& mat : m_matrices)
				mat = glm::rotate(mat, glm::radians(60.0f) * (float)dt, glm::vec3(0.0f, 0.0f, 1.0f));
			updateMatrixUniformBuffer();
		}
	}

	bool SWVulkanLattice::onUpdateUIOverlay(vks::UIOverlay* overlay)
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
					if (overlay->sliderInt("TessInner", &m_uniforms.tessInner, 0, 64)) updateLatticeUniformBuffer();
					if (overlay->sliderInt("TessOuter", &m_uniforms.tessOuter, 0, 64)) updateLatticeUniformBuffer();
				}
				if (overlay->comboBox("B-Function", &m_uniforms.bFunctionIndex, BFunctionNames)) updateLatticeUniformBuffer();
			}

			// TODO: Make this nice
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
								auto& mat = m_matrices[idx];
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

			/*if (overlay->header("Surfaces2", true))
			{
				overlay->listSurfaces(m_listItems, &m_selectedSurface);
			}*/
		}

		return rebuildCmd;
	}

	void SWVulkanLattice::setupLocalSurfaceVertex(OML::Locus& locus)
	{
		m_localSurfaceVertices.push_back(LocalSurfaceVertex(
			locus.controlPointIndex, locus.controlPointCount, locus.matrixIndex, 0));
	}

	void SWVulkanLattice::setupPatchVertices(OML::Patch& patch)
	{
		OML::Locus& locus00 = m_loci[patch.lociIndices[0]];
		OML::Locus& locus10 = m_loci[patch.lociIndices[1]];
		OML::Locus& locus01 = m_loci[patch.lociIndices[2]];
		OML::Locus& locus11 = m_loci[patch.lociIndices[3]];

		LocalSurfaceVertex localVert00;
		localVert00.controlPointIndex = locus00.controlPointIndex;
		localVert00.controlPointCount = locus00.controlPointCount;
		localVert00.matrixIndex = locus00.matrixIndex;
		localVert00.boundaryIndex = m_boundaries.size() - 4;
		m_patchVertices.push_back(localVert00);

		LocalSurfaceVertex localVert10;
		localVert10.controlPointIndex = locus10.controlPointIndex;
		localVert10.controlPointCount = locus10.controlPointCount;
		localVert10.matrixIndex = locus10.matrixIndex;
		localVert10.boundaryIndex = m_boundaries.size() - 3;
		m_patchVertices.push_back(localVert10);

		LocalSurfaceVertex localVert01;
		localVert01.controlPointIndex = locus01.controlPointIndex;
		localVert01.controlPointCount = locus01.controlPointCount;
		localVert01.matrixIndex = locus01.matrixIndex;
		localVert01.boundaryIndex = m_boundaries.size() - 2;
		m_patchVertices.push_back(localVert01);

		LocalSurfaceVertex localVert11;
		localVert11.controlPointIndex = locus11.controlPointIndex;
		localVert11.controlPointCount = locus11.controlPointCount;
		localVert11.matrixIndex = locus11.matrixIndex;
		localVert11.boundaryIndex = m_boundaries.size() - 1;
		m_patchVertices.push_back(localVert11);
	}

	void SWVulkanLattice::createDeviceLocalBuffer(VkBuffer& buffer, VkDeviceMemory& memory, void* data, uint32_t bufferSize, VkBufferUsageFlagBits usage)
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

	void SWVulkanLattice::createBuffers()
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

	void SWVulkanLattice::prepareUniformBuffers()
	{
		// Shared tessellation shader stages uniform buffer
		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&m_latticeUniformBuffer,
			sizeof(m_uniforms)));

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
			sizeof(glm::vec4) * m_numControlPoints + sizeof(OML::BoundaryInfo) * m_numPatches * 4
		));

		VK_CHECK_RESULT(m_patchUniformBuffer.map());

		updateLatticeUniformBuffer();
		updateMatrixUniformBuffer();
		updatePatchUniformBuffer();
	}

	void SWVulkanLattice::setupDescriptorSetLayouts()
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

	VkPipelineShaderStageCreateInfo SWVulkanLattice::loadShader(std::string fileName, VkShaderStageFlagBits stage)
	{
		VkPipelineShaderStageCreateInfo shaderStage = {};
		shaderStage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
		shaderStage.stage = stage; shaderStage.module = vks::tools::loadShader(fileName.c_str(), *m_device);
		shaderStage.pName = "main";
		assert(shaderStage.module != VK_NULL_HANDLE);
		return shaderStage;
	}

	void SWVulkanLattice::preparePipelines()
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
			GridVertex::GetBindingDescription()
		};

		// Attribute descriptions
		std::vector<VkVertexInputAttributeDescription> vertexInputAttributes =
			GridVertex::GetAttributeDesctiptions();

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

	void SWVulkanLattice::setupDescriptorPool()
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

	void SWVulkanLattice::setupDescriptorSets()
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

	void SWVulkanLattice::updateLatticeUniformBuffer()
	{
		memcpy(m_latticeUniformBuffer.mapped, &m_uniforms, sizeof(m_uniforms));
	}

	void SWVulkanLattice::updateMatrixUniformBuffer()
	{
		memcpy(m_matrixUniformBuffer.mapped, &m_matrices[0][0], sizeof(glm::mat4) * m_numLoci);
	}

	void SWVulkanLattice::updatePatchUniformBuffer()
	{
		std::vector<glm::vec4> patchUniforms;
		patchUniforms.insert(patchUniforms.end(), m_controlPoints.begin(), m_controlPoints.end());
		for (auto& boundary : m_boundaries)
			patchUniforms.push_back(glm::vec4(boundary.us, boundary.ue, boundary.vs, boundary.ve));
		memcpy(m_patchUniformBuffer.mapped, &patchUniforms[0],
			sizeof(glm::vec4) * m_numControlPoints + sizeof(OML::BoundaryInfo) * m_numPatches * 4);
	}
}