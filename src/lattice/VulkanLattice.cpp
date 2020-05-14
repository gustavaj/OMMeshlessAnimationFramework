#include "VulkanLattice.h"

#include "Shaders.h"
#include <chrono>

namespace OML
{

	std::vector<VkVertexInputBindingDescription> LocalSurfaceVertex::GetBindingDescriptions()
	{
		std::vector<VkVertexInputBindingDescription> bindings;
		bindings.resize(1);
		bindings[0].binding = 0;
		bindings[0].stride = sizeof(LocalSurfaceVertex);
		bindings[0].inputRate = VK_VERTEX_INPUT_RATE_VERTEX;
		return bindings;
	}

	std::vector<VkVertexInputAttributeDescription> LocalSurfaceVertex::GetAttributeDescriptions()
	{
		std::vector<VkVertexInputAttributeDescription> attributeDescriptions = {};
		attributeDescriptions.resize(2);
		attributeDescriptions[0].binding = 0;
		attributeDescriptions[0].location = 0;
		attributeDescriptions[0].format = VK_FORMAT_R32G32B32A32_UINT;
		attributeDescriptions[0].offset = offsetof(LocalSurfaceVertex, controlPointIndex);
		attributeDescriptions[1].binding = 0;
		attributeDescriptions[1].location = 1;
		attributeDescriptions[1].format = VK_FORMAT_R32G32B32_SFLOAT;
		attributeDescriptions[1].offset = offsetof(LocalSurfaceVertex, color);

		return attributeDescriptions;
	}

	VkVertexInputBindingDescription GridVertex::GetBindingDescription() {
		VkVertexInputBindingDescription bindingDescription = {};
		bindingDescription.binding = 0;
		bindingDescription.stride = sizeof(GridVertex);
		bindingDescription.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;
		return bindingDescription;
	}

	std::vector<VkVertexInputAttributeDescription> GridVertex::GetAttributeDesctiptions() {
		std::vector<VkVertexInputAttributeDescription> attributeDescriptions = {};
		attributeDescriptions.resize(2);

		attributeDescriptions[0].binding = 0;
		attributeDescriptions[0].location = 0;
		attributeDescriptions[0].format = VK_FORMAT_R32G32B32_SFLOAT;
		attributeDescriptions[0].offset = offsetof(GridVertex, pos);

		attributeDescriptions[1].binding = 0;
		attributeDescriptions[1].location = 1;
		attributeDescriptions[1].format = VK_FORMAT_R8G8B8_UINT;
		attributeDescriptions[1].offset = offsetof(GridVertex, col);

		return attributeDescriptions;
	}

	VulkanLattice::VulkanLattice()
		: VulkanLattice("")
	{
	}

	VulkanLattice::VulkanLattice(std::string name, LocalSurfaceType lsType, EvaluationMethod evalMethod)
		: Lattice(name, lsType), m_evalMethod(evalMethod), m_pointsPipeline(VK_NULL_HANDLE), m_linesPipeline(VK_NULL_HANDLE), m_localSurfacePipeline(VK_NULL_HANDLE),
		m_localSurfaceWireframePipeline(VK_NULL_HANDLE), m_patchPipeline(VK_NULL_HANDLE), m_patchWireframePipeline(VK_NULL_HANDLE),
		m_normalPipeline(VK_NULL_HANDLE), m_pipelineLayout(VK_NULL_HANDLE), m_descriptorSetLayout(VK_NULL_HANDLE), m_descriptorSet(VK_NULL_HANDLE),
		m_device(nullptr), m_vulkanDevice(nullptr), m_descriptorPool(nullptr), m_renderPass(nullptr),
		m_queue(nullptr), m_commandPool(nullptr), m_allocator(nullptr), m_selectedSurface(0)
	{
		m_menuSuffix = "##" + std::to_string(Lattice::Index++);
	}

	VulkanLattice::~VulkanLattice()
	{
	}

	void VulkanLattice::initVulkan(
		VkDevice* device, vks::VulkanDevice* vulkanDevice, VkQueue* queue, 
		VkCommandPool* commandPool, VkDescriptorPool* descriptorPool, 
		VkRenderPass* renderPass, VkAllocationCallbacks* allocator)
	{
		if (m_vulkanInitiated) return;

		m_device = device;
		m_vulkanDevice = vulkanDevice;
		m_queue = queue;
		m_commandPool = commandPool;
		m_descriptorPool = descriptorPool;
		m_renderPass = renderPass;
		m_allocator = allocator;

		setupVertices();
		setupQueryResultBuffer();
		createBuffers();
		prepareUniformBuffers();
		setupDescriptorSetLayouts();
		preparePipelines();
		setupDescriptorPool();
		setupDescriptorSets();

		m_listItems.resize(0);
		m_listItems.push_back(m_name);
		for (size_t i = 0; i < m_patches.size(); i++)
		{
			auto idx = std::to_string(i);
			m_listItems.push_back("Patch " + idx + " - p00");
			m_listItems.push_back("Patch " + idx + " - p10");
			m_listItems.push_back("Patch " + idx + " - p01");
			m_listItems.push_back("Patch " + idx + " - p11");
		}

		m_vulkanInitiated = true;
		m_destroyed = false;
	}

	void VulkanLattice::destroyVulkan()
	{
		if (m_destroyed) return;

		for (auto& lst : m_localSurfaceTextures) {
			lst.second.destroy();
		}

		for (auto& lst : m_batchTextures) {
			lst.destroy();
		}

		for (auto& shader : m_shaderModules) {
			vkDestroyShaderModule(*m_device, shader.second, m_allocator);
		}

		if (m_pointsPipeline) {
			vkDestroyPipeline(*m_device, m_pointsPipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_linesPipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_localSurfacePipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_localSurfaceWireframePipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_patchPipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_patchWireframePipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_normalPipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_displayPixelAccuracyPipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_displayTriSizePipeline, m_allocator);
			vkDestroyPipeline(*m_device, m_displaySurfaceAccuracyPipeline, m_allocator);

			vkDestroyPipelineLayout(*m_device, m_pipelineLayout, m_allocator);
			vkDestroyDescriptorSetLayout(*m_device, m_descriptorSetLayout, m_allocator);
			if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image) {
				vkDestroyDescriptorSetLayout(*m_device, m_localSamplerDescriptor.setLayout, m_allocator);
				vkDestroyDescriptorSetLayout(*m_device, m_samplerDescriptor.setLayout, m_allocator);
			}
			if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched) {
				vkDestroyDescriptorSetLayout(*m_device, m_batchDescriptor.setLayout, m_allocator);
			}
		}

		m_latticeUniformBuffer.destroy();
		m_matrixUniformBuffer.destroy();
		m_controlPointBuffer.destroy();
		m_boundariesBuffer.destroy();
		m_localSurfaceDataBuffer.destroy();

		if (m_pointsBuffer.buffer) {
			vkDestroyBuffer(*m_device, m_pointsBuffer.buffer, m_allocator);
			vkFreeMemory(*m_device, m_pointsBuffer.memory, m_allocator);
		}
		if (m_linesBuffer.buffer) {
			vkDestroyBuffer(*m_device, m_linesBuffer.buffer, m_allocator);
			vkFreeMemory(*m_device, m_linesBuffer.memory, m_allocator);
		}
		if (m_localSurfaceVertexBuffer.buffer) {
			vkDestroyBuffer(*m_device, m_localSurfaceVertexBuffer.buffer, m_allocator);
			vkFreeMemory(*m_device, m_localSurfaceVertexBuffer.memory, m_allocator);
		}
		if (m_patchVertexBuffer.buffer) {
			vkDestroyBuffer(*m_device, m_patchVertexBuffer.buffer, m_allocator);
			vkFreeMemory(*m_device, m_patchVertexBuffer.memory, m_allocator);
		}

		if (m_queryPool) {
			vkDestroyQueryPool(*m_device, m_queryPool, m_allocator);
			vkDestroyBuffer(*m_device, m_queryResult.buffer, m_allocator);
			vkFreeMemory(*m_device, m_queryResult.memory, m_allocator);
		}
		if (m_timingPool) {
			vkDestroyQueryPool(*m_device, m_timingPool, m_allocator);
			vkDestroyBuffer(*m_device, m_timingResult.buffer, m_allocator);
			vkFreeMemory(*m_device, m_timingResult.memory, m_allocator);
		}

		m_destroyed = true;
		m_vulkanInitiated = false;
	}

	void VulkanLattice::addToCommandbufferPreRenderpass(VkCommandBuffer& commandBuffer)
	{
		if (!m_draw) return;

		if (m_doPipelineQueries) {
			vkCmdResetQueryPool(commandBuffer, m_queryPool, 0, m_queryResult.count);
		}
		if (m_doPipelineTimings) {
			vkCmdResetQueryPool(commandBuffer, m_timingPool, 0, m_timingResult.count);
		}
	}

	void VulkanLattice::addToCommandbuffer(VkCommandBuffer& commandBuffer)
	{
		if (!m_draw) return;

		VkDeviceSize offsets[1] = { 0 };

		// Pipeline queries and timings
		if (m_doPipelineQueries) {
			vkCmdBeginQuery(commandBuffer, m_queryPool, 0, 0);
		}
		if (m_doPipelineTimings) {
			vkCmdWriteTimestamp(commandBuffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, m_timingPool, 0);
		}

		// Bind the descriptor set containing the uniform buffer and the common storage buffers
		vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, 
			m_pipelineLayout, 0, 1, &m_descriptorSet, 0, NULL);

		// Draw surface
		if (m_drawSurface)
		{
			// Bind the correct pipeline
			if (m_evalMethod != EvaluationMethod::Direct && m_surfaceColor == SurfaceColor::SurfAccuracy)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_displaySurfaceAccuracyPipeline);
			else if (m_surfaceColor == SurfaceColor::PixelAccuracy)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_displayPixelAccuracyPipeline);
			else if (m_surfaceColor == SurfaceColor::TriSize)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_displayTriSizePipeline);
			else if (m_wireframe)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_patchWireframePipeline);
			else
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_patchPipeline);
			vkCmdBindVertexBuffers(commandBuffer, 0, 1, &m_patchVertexBuffer.buffer, offsets);
			// EvalMethod dependent
			switch (m_evalMethod)
			{
			case EvaluationMethod::Direct:
			case EvaluationMethod::Pre_Sampled_Buffer: {
				vkCmdDraw(commandBuffer, m_patchVertexBuffer.count, 1, 0, 0);
				break;
			}
			case EvaluationMethod::Pre_Sampled_Image: {
				for (size_t i = 0; i < m_samplerDescriptor.sets.size(); i++)
				{
					vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, 
						m_pipelineLayout, 1, 1, &m_samplerDescriptor.sets[i], 0, NULL);
					vkCmdDraw(commandBuffer, 4, 1, i * 4, 0);
				}
				break;
			}
			case EvaluationMethod::Pre_Sampled_Image_Batched: {
				for (size_t i = 0; i < m_batchDescriptor.sets.size(); i++)
				{
					vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, 
						m_pipelineLayout, 1, 1, &m_batchDescriptor.sets[i], 0, NULL);
					vkCmdDraw(commandBuffer, m_batchInfo[i].count, 1, m_batchInfo[i].index, 0);
				}
				break;
			}
			}


			if (m_drawNormals)
			{
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_normalPipeline);
				// EvalMethod dependent
				switch (m_evalMethod)
				{
				case EvaluationMethod::Direct:
				case EvaluationMethod::Pre_Sampled_Buffer: {
					vkCmdDraw(commandBuffer, m_patchVertexBuffer.count, 1, 0, 0);
					break;
				}
				case EvaluationMethod::Pre_Sampled_Image: {
					for (size_t i = 0; i < m_samplerDescriptor.sets.size(); i++)
					{
						vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS,
							m_pipelineLayout, 1, 1, &m_samplerDescriptor.sets[i], 0, NULL);
						vkCmdDraw(commandBuffer, 4, 1, i * 4, 0);
					}
					break;
				}
				case EvaluationMethod::Pre_Sampled_Image_Batched: {
					for (size_t i = 0; i < m_batchDescriptor.sets.size(); i++)
					{
						vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS,
							m_pipelineLayout, 1, 1, &m_batchDescriptor.sets[i], 0, NULL);
						vkCmdDraw(commandBuffer, m_batchInfo[i].count, 1, m_batchInfo[i].index, 0);
					}
					break;
				}
				}
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
			// EvalMethod dependent
			if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image)
			{
				for (size_t i = 0; i < m_localSamplerDescriptor.sets.size(); i++)
				{
					vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, 
						m_pipelineLayout, 2, 1, &m_localSamplerDescriptor.sets[i], 0, NULL);
					vkCmdDraw(commandBuffer, 1, 1, i, 0);
				}
			}
			else
			{
				vkCmdDraw(commandBuffer, m_localSurfaceVertexBuffer.count, 1, 0, 0);
			}
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

		if (m_doPipelineQueries) {
			vkCmdEndQuery(commandBuffer, m_queryPool, 0);
		}
		if (m_doPipelineTimings) {
			vkCmdWriteTimestamp(commandBuffer, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, m_timingPool, 1);
		}
	}

	void VulkanLattice::onViewChanged(glm::mat4 projection, glm::mat4 view)
	{
		m_uniforms.projection = projection;
		m_view = view;
		m_frustum.update(projection * view);
		m_uniforms.frustumPlanes = m_frustum.planes;
		updateLatticeUniformBuffer();
	}

	void VulkanLattice::onWindowResized(float width, float height)
	{
		m_uniforms.windowSize = glm::vec2(width, height);
	}

	bool VulkanLattice::onUpdateUIOverlay(vks::UIOverlay* overlay)
	{
		bool rebuildCmd = false;

		ImGui::Separator();

		ImGui::PushID(m_menuSuffix.c_str());

		if (overlay->header(m_name.c_str(), false))
		{
			// Settings
			ImGui::Spacing(); ImGui::Spacing(); ImGui::SameLine();
			ImGui::BeginGroup();
			if (overlay->header("Settings", false))
			{
				if (overlay->checkBox("Render", &m_draw)) rebuildCmd = true;
				ImGui::Separator();
				if (overlay->checkBox("Draw Surface", &m_drawSurface)) rebuildCmd = true;
				if (overlay->checkBox("Draw Local", &m_drawLocalSurfaces)) rebuildCmd = true;
				if (overlay->checkBox("Draw Grid", &m_drawLatticeGrid)) rebuildCmd = true;
				if (overlay->checkBox("Draw Wireframe", &m_wireframe)) rebuildCmd = true;
				if (overlay->checkBox("Draw Normals", &m_drawNormals)) rebuildCmd = true;
				if (m_drawNormals)
				{
					if (overlay->sliderFloat("Length", &m_uniforms.normalLength, 1.0f, 10.0f)) updateLatticeUniformBuffer();
				}
				
				ImGui::Separator();
				
				if (overlay->comboBox("Surface Color", &m_surfaceColor, SurfaceColorNames)) {
					if (m_surfaceColor == SurfaceColor::SurfAccuracy && m_evalMethod == EvaluationMethod::Direct) {
						m_surfaceColor = SurfaceColor::Default;
					}
					rebuildCmd = true;
				}
				if (m_surfaceColor == SurfaceColor::SurfAccuracy)
				{
					if (overlay->sliderFloat("Error(e)", &m_uniforms.maxError, 0.1f, 10.0f)) updateLatticeUniformBuffer();
					overlay->text("w<=0.1e, g<=0.2e, b<=0.5e\ny<=e, r>e");
				}
				else if (m_surfaceColor == SurfaceColor::PixelAccuracy)
				{
					overlay->text("w<=0.5, g<=1.0, b<=2.0\ny<=5.0, r>5.0");
				}
				else if (m_surfaceColor == SurfaceColor::TriSize)
				{
					overlay->text("w<=1.0, g<=5.0, b<=10.0\ny<=20, r>20.0");
				}

				ImGui::Separator();
				if (overlay->checkBox("Frustum Culling", &m_uniforms.doCulling)) updateLatticeUniformBuffer();
				if (overlay->comboBox("TessFactorMethod", &m_uniforms.tessFactorMethod, TessFactorMethodNames)) updateLatticeUniformBuffer();
				switch (m_uniforms.tessFactorMethod)
				{
				case TessFactorMethod::Static: {
					if (overlay->sliderInt("TessInner", &m_uniforms.tessInner, 0, 64)) updateLatticeUniformBuffer();
					if (overlay->sliderInt("TessOuter", &m_uniforms.tessOuter, 0, 64)) updateLatticeUniformBuffer();
					break;
				}
				case TessFactorMethod::Dynamic: {
					if (overlay->sliderInt("PixelsPerEdge", &m_uniforms.pixelsPerEdge, 1, 50)) updateLatticeUniformBuffer();
					break;
				}
				}
				ImGui::Separator();
				if (overlay->comboBox("B-Function", &m_uniforms.bFunctionIndex, OML::BFunctionNames)) updateLatticeUniformBuffer();
				ImGui::Separator();
			}

			// Simulation
			if (overlay->header("Simulation", false))
			{
				overlay->checkBox("Simulate", &m_simulate);
				overlay->comboBox("Simulator", &m_simulatorIndex, OML::SimulatorNames);
				if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::NormalSin)) {
					overlay->sliderFloat("Min amp", &OML::Simulator::MinAmp, OML::Simulator::AmpRange.x, OML::Simulator::MaxAmp);
					overlay->sliderFloat("Max amp", &OML::Simulator::MaxAmp, OML::Simulator::MinAmp, OML::Simulator::AmpRange.y);
				}
				else if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::NormalRotation)) {
					overlay->sliderFloat("Min angle", &OML::Simulator::MinAngle, OML::Simulator::AngleRange.x, OML::Simulator::MaxAngle);
					overlay->sliderFloat("Max angle", &OML::Simulator::MaxAngle, OML::Simulator::MinAngle, OML::Simulator::AngleRange.y);
				}
				else if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::XYScale)) {
					overlay->sliderFloat("Min scale - 1", &OML::Simulator::MinScale, OML::Simulator::ScaleRange.x, OML::Simulator::MaxScale);
					overlay->sliderFloat("Max scale - 1", &OML::Simulator::MaxScale, OML::Simulator::MinScale, OML::Simulator::ScaleRange.y);
				}
				overlay->sliderFloat("Min speed", &OML::Simulator::MinSpeed, OML::Simulator::SpeedRange.x, OML::Simulator::MaxSpeed);
				overlay->sliderFloat("Max speed", &OML::Simulator::MaxSpeed, OML::Simulator::MinSpeed, OML::Simulator::SpeedRange.y);
				bool simulatorActive = m_simulators.find(static_cast<OML::SimulatorTypes>(m_simulatorIndex)) != m_simulators.end();
				std::string addBtnTitle = simulatorActive ? "Update" : "Add";
				if (overlay->button(addBtnTitle.c_str())) {
					if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::NormalSin)) 
						addNormalSinSimulation();
					else if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::NormalRotation)) 
						addNormalRotationSimulation();
					else if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::XYScale)) 
						addXYScalingSimulation();
				}
				if (simulatorActive) {
					ImGui::SameLine();
					if (overlay->button("Remove")) {
						if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::NormalSin)) 
							removeSimulator(OML::SimulatorTypes::NormalSin);
						else if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::NormalRotation))
							removeSimulator(OML::SimulatorTypes::NormalRotation);
						else if (m_simulatorIndex == static_cast<int>(OML::SimulatorTypes::XYScale)) 
							removeSimulator(OML::SimulatorTypes::XYScale);
					}
				}
			}

			// Pipeline stats
			if (overlay->header("Statistics", false))
			{
				if (overlay->checkBox("Stats", &m_doPipelineQueries)) rebuildCmd = true;
				if (m_doPipelineQueries) {
					overlay->text("Input Assembler Vertices: %d", m_pipelineStats[0]);
					overlay->text("Input Assembler Primitives: %d", m_pipelineStats[1]);
					overlay->text("Vert invocations: %d", m_pipelineStats[2]);
					overlay->text("TCS patches: %d", m_pipelineStats[8]);
					overlay->text("TES invocations: %d", m_pipelineStats[9]);
					overlay->text("Geom invocations: %d", m_pipelineStats[3]);
					overlay->text("Geom primitives: %d", m_pipelineStats[4]);
					overlay->text("Clipping invocations: %d", m_pipelineStats[5]);
					overlay->text("Clipping primitives: %d", m_pipelineStats[6]);
					overlay->text("Frag invocations: %d", m_pipelineStats[7]);

					ImGui::Separator();
				}

				if (overlay->checkBox("Timings", &m_doPipelineTimings)) rebuildCmd = true;
				if (m_doPipelineTimings) {
					const uint64_t diff = m_pipelineTimings[1] - m_pipelineTimings[0];
					double micros = diff * m_timestampPeriod / 1000000;
					if (micros > 10000000.0f) micros = 0.0f;
					overlay->text("Rendering time: %.4f ms", micros);
				}

				overlay->text("GPU memory usage: \n%d bytes", m_deviceMemoryUsage);
			}

			// For editing local surfaces
			if (overlay->header("Surfaces", false))
			{
				std::vector<const char*> items(m_listItems.size());
				for (size_t i = 0; i < m_listItems.size(); i++) items[i] = m_listItems[i].c_str();
				ImGui::ListBox("", &m_selectedSurface, items.data(), items.size(), 9);
				glm::mat4* mat;
				if (m_selectedSurface == 0) {
					mat = &m_matrix;
				}
				else {
					int matIdx = m_loci[m_patches[(m_selectedSurface - 1) / 4].lociIndices[(m_selectedSurface - 1) % 4]].matrixIndex;
					mat = &m_matrices[matIdx];
				}

				bool update = false;

				ImGui::SameLine();
				ImGui::BeginGroup();

				ImGui::Text(m_listItems[m_selectedSurface].c_str());
				ImGui::Separator();

				ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(5.0f, 1.0f));

				// Translation
				ImGui::BeginGroup();

				ImGui::Text("Translate");
				if (ImGui::InputFloat("x", &(*mat)[3][0], 0.5f, 1.0f, "%.1f")) update = true;
				if (ImGui::InputFloat("y", &(*mat)[3][1], 0.5f, 1.0f, "%.1f")) update = true;
				if (ImGui::InputFloat("z", &(*mat)[3][2], 0.5f, 1.0f, "%.1f")) update = true;

				ImGui::EndGroup();

				const float button_size = ImGui::GetFrameHeight();

				// Rotation
				ImGui::BeginGroup();

				float itemInnerSpacingX = ImGui::GetStyle().ItemInnerSpacing.x;

				ImGui::Text("Rot");
				ImGui::Text("x");
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("-##rx", ImVec2(button_size, button_size))) {
					(*mat) = glm::rotate((*mat), glm::radians(-10.0f), glm::vec3(1, 0, 0));
					update = true;
				}
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("+##rx", ImVec2(button_size, button_size))) {
					(*mat) = glm::rotate((*mat), glm::radians(10.0f), glm::vec3(1, 0, 0));
					update = true;
				}
				ImGui::Text("y");
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("-##ry", ImVec2(button_size, button_size))) {
					(*mat) = glm::rotate((*mat), glm::radians(-10.0f), glm::vec3(0, 1, 0));
					update = true;
				}
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("+##ry", ImVec2(button_size, button_size))) {
					(*mat) = glm::rotate((*mat), glm::radians(10.0f), glm::vec3(0, 1, 0));
					update = true;
				}
				ImGui::Text("z");
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("-##rz", ImVec2(button_size, button_size))) {
					(*mat) = glm::rotate((*mat), glm::radians(-10.0f), glm::vec3(0, 0, 1));
					update = true;
				}
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("+##rz", ImVec2(button_size, button_size))) {
					(*mat) = glm::rotate((*mat), glm::radians(10.0f), glm::vec3(0, 0, 1));
					update = true;
				}

				ImGui::EndGroup();

				ImGui::SameLine();
				ImGui::Spacing();
				ImGui::SameLine();

				// Scaling
				ImGui::BeginGroup();

				ImGui::Text("Scale");
				ImGui::Text("x");
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("-##sx", ImVec2(button_size, button_size))) {
					(*mat) = glm::scale((*mat), glm::vec3(0.9, 1.0, 1.0));
					update = true;
				}
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("+##sx", ImVec2(button_size, button_size))) {
					(*mat) = glm::scale((*mat), glm::vec3(1.1, 1.0, 1.0));
					update = true;
				}
				ImGui::Text("y");
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("-##sy", ImVec2(button_size, button_size))) {
					(*mat) = glm::scale((*mat), glm::vec3(1.0, 0.9, 1.0));
					update = true;
				}
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("+##sy", ImVec2(button_size, button_size))) {
					(*mat) = glm::scale((*mat), glm::vec3(1.0, 1.1, 1.0));
					update = true;
				}
				ImGui::Text("z");
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("-##sz", ImVec2(button_size, button_size))) {
					(*mat) = glm::scale((*mat), glm::vec3(1.0, 1.0, 0.9));
					update = true;
				}
				ImGui::SameLine(0, itemInnerSpacingX);
				if (ImGui::Button("+##sz", ImVec2(button_size, button_size))) {
					(*mat) = glm::scale((*mat), glm::vec3(1.0, 1.0, 1.1));
					update = true;
				}

				ImGui::EndGroup();

				ImGui::PopStyleVar(1);

				ImGui::EndGroup();

				if (update) {
					if (m_selectedSurface == 0) updateLatticeUniformBuffer();
					else updateMatrixUniformBuffer();
				}
			}

			ImGui::EndGroup();
		}

		//ImGui::ShowStyleEditor();
		ImGui::PopID();

		return rebuildCmd;
	}

	void VulkanLattice::CheckAndSetupRequiredPhysicalDeviceFeatures(VkPhysicalDeviceFeatures& deviceFeatures, VkPhysicalDeviceFeatures& enabledFeatures)
	{
		assert(deviceFeatures.tessellationShader && deviceFeatures.geometryShader &&
			deviceFeatures.fillModeNonSolid && deviceFeatures.pipelineStatisticsQuery &&
			deviceFeatures.vertexPipelineStoresAndAtomics && deviceFeatures.fragmentStoresAndAtomics);

		enabledFeatures.tessellationShader = VK_TRUE;				// For tessellation
		enabledFeatures.geometryShader = VK_TRUE;					// For drawing normals
		enabledFeatures.fillModeNonSolid = VK_TRUE;					// For drawing wireframe
		enabledFeatures.pipelineStatisticsQuery = VK_TRUE;			// For getting stats from the rendering
		enabledFeatures.vertexPipelineStoresAndAtomics = VK_TRUE;	// For using shader storage buffers
		enabledFeatures.fragmentStoresAndAtomics = VK_TRUE;			// For using shader storage buffers in fragment shader
	}

	void VulkanLattice::localUpdate(double dt)
	{
		if (m_simulate) {
			updateMatrixUniformBuffer();
		}

		getQueryResults();
	}

	void VulkanLattice::createDeviceLocalBuffer(VkBuffer& buffer, VkDeviceMemory& memory, void* data, uint32_t bufferSize, VkBufferUsageFlagBits usage)
	{
		if (buffer != VK_NULL_HANDLE) {
			vkDestroyBuffer(*m_device, buffer, m_allocator);
		}
		if (memory != VK_NULL_HANDLE) {
			vkFreeMemory(*m_device, memory, m_allocator);
		}

		m_deviceMemoryUsage += bufferSize;

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

		vkDestroyBuffer(*m_device, stagingBuffer.buffer, m_allocator);
		vkFreeMemory(*m_device, stagingBuffer.memory, m_allocator);
	}

	void VulkanLattice::setupVertices()
	{
		Timer::Start("setupVertices", "VulkanLattice::setupVertices()");

		// Eval method dependent

		// Set up local surface vertices
		m_localSurfaceVertices.resize(m_loci.size());
		for (size_t i = 0; i < m_loci.size(); i++)
		{
			if (m_evalMethod == EvaluationMethod::Pre_Sampled_Buffer) {
				// Evaluate local surfaces and load them as textures.
				auto res = m_LSIdxToLSBufferMap.find(m_loci[i].controlPointIndex);
				if (res == m_LSIdxToLSBufferMap.end()) {
					std::vector<glm::vec3> controlPoints(m_loci[i].controlPointCount);
					for (size_t j = 0; j < controlPoints.size(); j++)
					{
						controlPoints[j] = glm::vec3(m_controlPoints[m_loci[i].controlPointIndex + j]);
					}
					auto it = m_LSIdxToLSBufferMap.insert({ m_loci[i].controlPointIndex,
						m_LSBuffer.addLocalSurface(controlPoints, m_lsType) });
				}

				m_localSurfaceVertices[i] = LocalSurfaceVertex(m_LSIdxToLSBufferMap[m_loci[i].controlPointIndex],
					m_loci[i].controlPointCount, m_loci[i].matrixIndex, 0, m_loci[i].color);
			}
			else
			{
				m_localSurfaceVertices[i] = LocalSurfaceVertex(m_loci[i].controlPointIndex,
					m_loci[i].controlPointCount, m_loci[i].matrixIndex, 0, m_loci[i].color);
			}

			if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image)
			{
				// Evaluate local surfaces and load them as textures.
				auto res = m_localSurfaceTextures.find(m_loci[i].controlPointIndex);
				if (res == m_localSurfaceTextures.end()) {
					auto it = m_localSurfaceTextures.insert({ m_loci[i].controlPointIndex,
						LocalSurfaceTexture(m_device, m_vulkanDevice, m_commandPool, m_queue, m_allocator) }).first;
					std::vector<glm::vec3> controlPoints(m_loci[i].controlPointCount);
					for (size_t j = 0; j < controlPoints.size(); j++)
					{
						controlPoints[j] = glm::vec3(m_controlPoints[m_loci[i].controlPointIndex + j]);
					}
					it->second.loadLocalSurface(controlPoints, NUM_SAMPLES_U, NUM_SAMPLES_V, m_lsType);
					m_deviceMemoryUsage += it->second.memoryUsage();
				}
			}
		}

		// Set up patch vertices
		m_patchVertices.resize(0);
		int batch = -1; // Gets ++'ed when i == 0
		for (size_t i = 0; i < m_patches.size(); i++)
		{
			OML::Locus& locus00 = m_loci[m_patches[i].lociIndices[0]];
			OML::Locus& locus10 = m_loci[m_patches[i].lociIndices[1]];
			OML::Locus& locus01 = m_loci[m_patches[i].lociIndices[2]];
			OML::Locus& locus11 = m_loci[m_patches[i].lociIndices[3]];

			LocalSurfaceVertex localVert00;
			localVert00.controlPointIndex = locus00.controlPointIndex;
			localVert00.controlPointCount = locus00.controlPointCount;
			localVert00.matrixIndex = locus00.matrixIndex;
			localVert00.boundaryIndex = locus00.boundaryIndices[m_patches[i].faceIdx];
			localVert00.color = m_patches[i].color;

			LocalSurfaceVertex localVert10;
			localVert10.controlPointIndex = locus10.controlPointIndex;
			localVert10.controlPointCount = locus10.controlPointCount;
			localVert10.matrixIndex = locus10.matrixIndex;
			localVert10.boundaryIndex = locus10.boundaryIndices[m_patches[i].faceIdx];
			localVert10.color = m_patches[i].color;

			LocalSurfaceVertex localVert01;
			localVert01.controlPointIndex = locus01.controlPointIndex;
			localVert01.controlPointCount = locus01.controlPointCount;
			localVert01.matrixIndex = locus01.matrixIndex;
			localVert01.boundaryIndex = locus01.boundaryIndices[m_patches[i].faceIdx];
			localVert01.color = m_patches[i].color;

			LocalSurfaceVertex localVert11;
			localVert11.controlPointIndex = locus11.controlPointIndex;
			localVert11.controlPointCount = locus11.controlPointCount;
			localVert11.matrixIndex = locus11.matrixIndex;
			localVert11.boundaryIndex = locus11.boundaryIndices[m_patches[i].faceIdx];
			localVert11.color = m_patches[i].color;

			if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
			{
				if (i % NUM_PATCHES_PER_BATCH == 0) {
					batch++;
					// TODO: If the last batch can not fill up all the spaces, make it smaller!
					m_batchTextures.push_back(LocalSurfaceTextureBatch(
						BATCH_ROWS, BATCH_COLS, NUM_SAMPLES_U, NUM_SAMPLES_V,
						m_device, m_vulkanDevice, m_commandPool, m_queue, m_allocator
					));
					if (m_patches.size() - i < NUM_PATCHES_PER_BATCH) {
						m_batchInfo.push_back({ static_cast<uint32_t>(i * 4), static_cast<uint32_t>(m_patches.size() - i) * 4 });
					}
					else {
						m_batchInfo.push_back({ static_cast<uint32_t>(i * 4), NUM_PATCHES_PER_BATCH * 4 });
					}
				}

				auto coords = m_batchTextures[batch].addPatch(
					{ m_controlPoints.begin() + locus00.controlPointIndex,
					m_controlPoints.begin() + locus00.controlPointIndex + locus00.controlPointCount },
					m_boundaries[locus00.boundaryIndices[m_patches[i].faceIdx]],
					{ m_controlPoints.begin() + locus10.controlPointIndex,
					m_controlPoints.begin() + locus10.controlPointIndex + locus10.controlPointCount },
					m_boundaries[locus10.boundaryIndices[m_patches[i].faceIdx]],
					{ m_controlPoints.begin() + locus01.controlPointIndex,
					m_controlPoints.begin() + locus01.controlPointIndex + locus01.controlPointCount },
					m_boundaries[locus01.boundaryIndices[m_patches[i].faceIdx]],
					{ m_controlPoints.begin() + locus11.controlPointIndex,
					m_controlPoints.begin() + locus11.controlPointIndex + locus11.controlPointCount },
					m_boundaries[locus11.boundaryIndices[m_patches[i].faceIdx]],
					m_lsType
				);

				localVert00.controlPointIndex = coords.first;
				localVert00.controlPointCount = coords.second;
				localVert00.boundaryIndex = NUM_SAMPLES_U;
				localVert10.controlPointIndex = coords.first;
				localVert10.controlPointCount = coords.second;
				localVert10.boundaryIndex = NUM_SAMPLES_U;
				localVert01.controlPointIndex = coords.first;
				localVert01.controlPointCount = coords.second;
				localVert01.boundaryIndex = NUM_SAMPLES_U;
				localVert11.controlPointIndex = coords.first;
				localVert11.controlPointCount = coords.second;
				localVert11.boundaryIndex = NUM_SAMPLES_U;
			}
			else if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image)
			{
				PatchSamplerInfo sampler;
				sampler.p00Sampler = &m_localSurfaceTextures[locus00.controlPointIndex];
				sampler.p10Sampler = &m_localSurfaceTextures[locus10.controlPointIndex];
				sampler.p01Sampler = &m_localSurfaceTextures[locus01.controlPointIndex];
				sampler.p11Sampler = &m_localSurfaceTextures[locus11.controlPointIndex];
				m_patchSamplers.push_back(sampler);
			}
			else if (m_evalMethod == EvaluationMethod::Pre_Sampled_Buffer)
			{
				localVert00.controlPointIndex = m_LSIdxToLSBufferMap[locus00.controlPointIndex];
				localVert10.controlPointIndex = m_LSIdxToLSBufferMap[locus10.controlPointIndex];
				localVert01.controlPointIndex = m_LSIdxToLSBufferMap[locus01.controlPointIndex];
				localVert11.controlPointIndex = m_LSIdxToLSBufferMap[locus11.controlPointIndex];
			}


			m_patchVertices.push_back(localVert00);
			m_patchVertices.push_back(localVert10);
			m_patchVertices.push_back(localVert01);
			m_patchVertices.push_back(localVert11);
		}

		if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
		{
			for (auto& batch : m_batchTextures)
			{
				batch.allocateMemory();
				m_deviceMemoryUsage += batch.memoryUsage();
			}
			m_batchDescriptor.sets.resize(m_batchInfo.size());
		}
		else if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image)
		{
			m_samplerDescriptor.sets.resize(m_patchSamplers.size());
			m_localSamplerDescriptor.sets.resize(m_localSurfaceVertices.size());
		}

		Timer::Stop("setupVertices", "VulkanLattice::setupVertices()");
	}

	void VulkanLattice::createBuffers()
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

	void VulkanLattice::prepareUniformBuffers()
	{
		// Shared tessellation shader stages uniform buffer
		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&m_latticeUniformBuffer,
			sizeof(m_uniforms)));

		// Map persistent
		VK_CHECK_RESULT(m_latticeUniformBuffer.map());
		m_deviceMemoryUsage += m_latticeUniformBuffer.size;

		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&m_matrixUniformBuffer,
			sizeof(glm::mat4) * m_matrices.size()
		));

		VK_CHECK_RESULT(m_matrixUniformBuffer.map());
		m_deviceMemoryUsage += m_matrixUniformBuffer.size;

		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
			VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
			&m_controlPointBuffer,
			sizeof(glm::vec4) * m_controlPoints.size()
		));

		VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
			VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
			&m_boundariesBuffer,
			sizeof(OML::BoundaryInfo) * m_boundaries.size()
		));

		if (m_evalMethod == EvaluationMethod::Pre_Sampled_Buffer)
		{
			VK_CHECK_RESULT(m_vulkanDevice->createBuffer(
				VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
				VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
				&m_localSurfaceDataBuffer,
				m_LSBuffer.bufferSize()
			));
		}

		updateLatticeUniformBuffer();
		updateMatrixUniformBuffer();
		uploadStorageBuffers();
	}

	void VulkanLattice::uploadStorageBuffers()
	{
		createDeviceLocalBuffer(m_controlPointBuffer.buffer, m_controlPointBuffer.memory,
			m_controlPoints.data(), m_controlPointBuffer.size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
		m_controlPointBuffer.descriptor.buffer = m_controlPointBuffer.buffer;
		createDeviceLocalBuffer(m_boundariesBuffer.buffer, m_boundariesBuffer.memory,
			m_boundaries.data(), m_boundariesBuffer.size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
		m_boundariesBuffer.descriptor.buffer = m_boundariesBuffer.buffer;

		if (m_evalMethod == EvaluationMethod::Pre_Sampled_Buffer)
		{
			auto lsData = m_LSBuffer.getData();
			createDeviceLocalBuffer(m_localSurfaceDataBuffer.buffer, m_localSurfaceDataBuffer.memory,
				lsData.data(), m_localSurfaceDataBuffer.size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
			m_localSurfaceDataBuffer.descriptor.buffer = m_localSurfaceDataBuffer.buffer;
		}
	}

	void VulkanLattice::setupDescriptorSetLayouts()
	{
		std::vector<VkDescriptorSetLayout> descriptorLayouts;

		// LAttice
		std::vector<VkDescriptorSetLayoutBinding> setLayoutBindings = {
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT |
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_GEOMETRY_BIT |
				VK_SHADER_STAGE_FRAGMENT_BIT,
				0),
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT | VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_FRAGMENT_BIT,
				1),
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT | VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_FRAGMENT_BIT,
				2),
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT | VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_FRAGMENT_BIT,
				3)
		};
		if (m_evalMethod == EvaluationMethod::Pre_Sampled_Buffer)
		{
			setLayoutBindings.push_back(
				vks::initializers::descriptorSetLayoutBinding(
					VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
					VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_FRAGMENT_BIT,
					4)
			);
		}

		VkDescriptorSetLayoutCreateInfo descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(
			setLayoutBindings.data(), static_cast<uint32_t>(setLayoutBindings.size()));
		VK_CHECK_RESULT(vkCreateDescriptorSetLayout(*m_device, &descriptorLayout,
			m_allocator, &m_descriptorSetLayout));
		descriptorLayouts.push_back(m_descriptorSetLayout);

		if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image)
		{
			setLayoutBindings = {
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
				0),
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
				1),
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
				2),
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
				3)
			};
			descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(
				setLayoutBindings.data(), static_cast<uint32_t>(setLayoutBindings.size()));
			VK_CHECK_RESULT(vkCreateDescriptorSetLayout(*m_device, &descriptorLayout,
				m_allocator, &m_samplerDescriptor.setLayout));
			descriptorLayouts.push_back(m_samplerDescriptor.setLayout);

			setLayoutBindings = {
				vks::initializers::descriptorSetLayoutBinding(
					VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
					VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
					0)
			};
			descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(
				setLayoutBindings.data(), static_cast<uint32_t>(setLayoutBindings.size()));
			VK_CHECK_RESULT(vkCreateDescriptorSetLayout(*m_device, &descriptorLayout,
				m_allocator, &m_localSamplerDescriptor.setLayout));
			descriptorLayouts.push_back(m_localSamplerDescriptor.setLayout);
		}
		else if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
		{
			setLayoutBindings = {
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
				0)
			};
			descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(
				setLayoutBindings.data(), static_cast<uint32_t>(setLayoutBindings.size()));
			VK_CHECK_RESULT(vkCreateDescriptorSetLayout(*m_device, &descriptorLayout,
				m_allocator, &m_batchDescriptor.setLayout));
			descriptorLayouts.push_back(m_batchDescriptor.setLayout);
		}

		// Pipeline Layout
		VkPipelineLayoutCreateInfo pipelineLayoutCreateInfo =
			vks::initializers::pipelineLayoutCreateInfo(descriptorLayouts.data(), descriptorLayouts.size());

		VK_CHECK_RESULT(vkCreatePipelineLayout(*m_device, &pipelineLayoutCreateInfo, m_allocator, &m_pipelineLayout));
	}

	VkPipelineShaderStageCreateInfo VulkanLattice::loadShader(std::pair<std::string, std::vector<uint32_t>&> src, VkShaderStageFlagBits stage)
	{
		VkPipelineShaderStageCreateInfo shaderStage = {};
		shaderStage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
		shaderStage.stage = stage;
		if (m_shaderModules.find(src.first) == m_shaderModules.end()) {
			VkShaderModule shaderModule;
			VkShaderModuleCreateInfo moduleCreateInfo{};
			moduleCreateInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
			moduleCreateInfo.codeSize = src.second.size() * sizeof(uint32_t);
			moduleCreateInfo.pCode = src.second.data();

			VK_CHECK_RESULT(vkCreateShaderModule(*m_device, &moduleCreateInfo, NULL, &shaderModule));
			m_shaderModules.insert({ src.first, shaderModule });
		}
		shaderStage.module = m_shaderModules[src.first];
		shaderStage.pName = "main";
		assert(shaderStage.module != VK_NULL_HANDLE);
		return shaderStage;
	}

	void VulkanLattice::preparePipelines()
	{
		Timer::Start("preparePipelines", "VulkanLattice::preparePipelines()");

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

		// TODO: Probably remove this dynamic state as it is never changed, aside from maybe viewport?
		// Probably has some performance penalty.
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
		shaderStages[0] = loadShader(OML::Shaders::GetPosColorPassVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		shaderStages[1] = loadShader(OML::Shaders::GetUintFlatColorFragShader(), VK_SHADER_STAGE_FRAGMENT_BIT);

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

		OML::ShaderOptions options = {};
		options.numControl = m_controlPoints.size();
		options.numLocal = m_matrices.size();
		options.numPatches = m_patches.size();
		options.numSamplesU = NUM_SAMPLES_U;
		options.numSamplesV = NUM_SAMPLES_V;

		// Local surface pipeline
		std::array<VkPipelineShaderStageCreateInfo, 4> localShaderStages;
		localShaderStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		localShaderStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(1, m_lsType, options), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		localShaderStages[2] = loadShader(OML::Shaders::GetTeseShader(m_lsType, OML::TeseShaderType::Local, 
			(m_evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched ? OML::EvaluationMethod::Direct : m_evalMethod), options),
			VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		localShaderStages[3] = loadShader(OML::Shaders::GetFlatColorFragShader(), VK_SHADER_STAGE_FRAGMENT_BIT);
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
		bsShaderStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		bsShaderStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(4, m_lsType, options), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		bsShaderStages[2] = loadShader(OML::Shaders::GetTeseShader(m_lsType, OML::TeseShaderType::Lattice, m_evalMethod, options),
			VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		bsShaderStages[3] = loadShader(OML::Shaders::GetShadedColorFragShader(), VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(bsShaderStages.size());
		pipelineCreateInfo.pStages = bsShaderStages.data();

		rasterizationState.polygonMode = VK_POLYGON_MODE_LINE;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_patchWireframePipeline));

		rasterizationState.polygonMode = VK_POLYGON_MODE_FILL;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_patchPipeline));

		// Surface accuracy pipeline
		if (m_evalMethod != EvaluationMethod::Direct)
		{
			std::array<VkPipelineShaderStageCreateInfo, 4> surfaceAccuracyStages;
			surfaceAccuracyStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
			surfaceAccuracyStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(4, m_lsType, options), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
			surfaceAccuracyStages[2] = loadShader(OML::Shaders::GetTeseShader(m_lsType, OML::TeseShaderType::Surf_Accuracy, m_evalMethod, options),
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
			surfaceAccuracyStages[3] = loadShader(OML::Shaders::GetFlatColorFragShader(), VK_SHADER_STAGE_FRAGMENT_BIT);
			pipelineCreateInfo.stageCount = static_cast<uint32_t>(surfaceAccuracyStages.size());
			pipelineCreateInfo.pStages = surfaceAccuracyStages.data();

			VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_displaySurfaceAccuracyPipeline));
		}

		// Pixel accuracy display pipeline
		std::array<VkPipelineShaderStageCreateInfo, 4> pixelAccuracyStages;
		pixelAccuracyStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		pixelAccuracyStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(4, m_lsType, options), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		pixelAccuracyStages[2] = loadShader(OML::Shaders::GetTeseShader(m_lsType, OML::TeseShaderType::Pixel_Accuracy, m_evalMethod, options),
			VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		pixelAccuracyStages[3] = loadShader(OML::Shaders::GetBiQuadLatticePixelAccuracyFragShader(
			m_lsType, options.numControl, options.numLocal, options.numPatches), VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(pixelAccuracyStages.size());
		pipelineCreateInfo.pStages = pixelAccuracyStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_displayPixelAccuracyPipeline));

		// Triangle Size display pipeline
		std::array<VkPipelineShaderStageCreateInfo, 5> triSizeStages;
		triSizeStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		triSizeStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(4, m_lsType, options), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		triSizeStages[2] = loadShader(OML::Shaders::GetTeseShader(m_lsType, OML::TeseShaderType::Lattice, m_evalMethod, options), VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		triSizeStages[3] = loadShader(OML::Shaders::GetTriSizeGeomShader(), VK_SHADER_STAGE_GEOMETRY_BIT);
		triSizeStages[4] = loadShader(OML::Shaders::GetShadedColorFragShader(), VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(triSizeStages.size());
		pipelineCreateInfo.pStages = triSizeStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_displayTriSizePipeline));

		// Normals pipeline
		std::array<VkPipelineShaderStageCreateInfo, 5> normalShaderStages;
		normalShaderStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		normalShaderStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(4, m_lsType, options), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		normalShaderStages[2] = loadShader(OML::Shaders::GetTeseShader(m_lsType, OML::TeseShaderType::Normals, m_evalMethod, options),
			VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		normalShaderStages[3] = loadShader(OML::Shaders::GetLatticeNormalsGeomShader(), VK_SHADER_STAGE_GEOMETRY_BIT);
		normalShaderStages[4] = loadShader(OML::Shaders::GetFlatColorFragShader(), VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(normalShaderStages.size());
		pipelineCreateInfo.pStages = normalShaderStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_normalPipeline));

		Shaders::SortShaderNames();

		Timer::Stop("preparePipelines", "VulkanLattice::preparePipelines()");
	}

	void VulkanLattice::setupDescriptorPool()
	{
		std::vector<VkDescriptorPoolSize> poolSizes =
		{
			vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1),
			vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 
				m_evalMethod == EvaluationMethod::Pre_Sampled_Buffer ? 4 : 3)
		};

		uint32_t maxSets = 1;

		if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image) 
		{
			maxSets = 4 * m_samplerDescriptor.sets.size() + m_localSamplerDescriptor.sets.size() + 1;
			poolSizes.push_back(
				vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, maxSets - 1));
		}
		else if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
		{
			maxSets = MAX_BATCHES + 1;
			poolSizes.push_back(
				vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, maxSets - 1));
		}

		VkDescriptorPoolCreateInfo descriptorPoolInfo =
			vks::initializers::descriptorPoolCreateInfo(
				static_cast<uint32_t>(poolSizes.size()),
				poolSizes.data(),
				maxSets);

		VK_CHECK_RESULT(vkCreateDescriptorPool(*m_device, &descriptorPoolInfo, m_allocator, m_descriptorPool));
	}

	void VulkanLattice::setupDescriptorSets()
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
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				1,
				&m_matrixUniformBuffer.descriptor),
			vks::initializers::writeDescriptorSet(
				m_descriptorSet,
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				2,
				&m_controlPointBuffer.descriptor),
			vks::initializers::writeDescriptorSet(
				m_descriptorSet,
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				3,
				&m_boundariesBuffer.descriptor)
		};
		if (m_evalMethod == EvaluationMethod::Pre_Sampled_Buffer)
		{
			writeDescriptorSets.push_back(
				vks::initializers::writeDescriptorSet(
					m_descriptorSet,
					VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
					4,
					&m_localSurfaceDataBuffer.descriptor)
			);
		}

		vkUpdateDescriptorSets(*m_device, static_cast<uint32_t>(writeDescriptorSets.size()), 
			writeDescriptorSets.data(), 0, NULL);

		if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image)
		{
			for (size_t i = 0; i < m_samplerDescriptor.sets.size(); i++)
			{
				allocInfo = vks::initializers::descriptorSetAllocateInfo(
					*m_descriptorPool, &m_samplerDescriptor.setLayout, 1);
				VK_CHECK_RESULT(vkAllocateDescriptorSets(*m_device, &allocInfo, &m_samplerDescriptor.sets[i]));

				writeDescriptorSets = {
					vks::initializers::writeDescriptorSet(
						m_samplerDescriptor.sets[i],
						VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
						0,
						m_patchSamplers[i].p00Sampler->descriptor()),
					vks::initializers::writeDescriptorSet(
						m_samplerDescriptor.sets[i],
						VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
						1,
						m_patchSamplers[i].p10Sampler->descriptor()),
					vks::initializers::writeDescriptorSet(
						m_samplerDescriptor.sets[i],
						VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
						2,
						m_patchSamplers[i].p01Sampler->descriptor()),
					vks::initializers::writeDescriptorSet(
						m_samplerDescriptor.sets[i],
						VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
						3,
						m_patchSamplers[i].p11Sampler->descriptor())
				};

				vkUpdateDescriptorSets(*m_device, static_cast<uint32_t>(writeDescriptorSets.size()), 
					writeDescriptorSets.data(), 0, NULL);
			}

			for (size_t i = 0; i < m_localSamplerDescriptor.sets.size(); i++)
			{
				allocInfo = vks::initializers::descriptorSetAllocateInfo(
					*m_descriptorPool, &m_localSamplerDescriptor.setLayout, 1);
				VK_CHECK_RESULT(vkAllocateDescriptorSets(*m_device, &allocInfo, &m_localSamplerDescriptor.sets[i]));

				writeDescriptorSets = {
					vks::initializers::writeDescriptorSet(
						m_localSamplerDescriptor.sets[i],
						VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
						0,
						m_localSurfaceTextures[m_localSurfaceVertices[i].controlPointIndex].descriptor())
				};

				vkUpdateDescriptorSets(*m_device, static_cast<uint32_t>(writeDescriptorSets.size()), 
					writeDescriptorSets.data(), 0, NULL);
			}
		}
		else if (m_evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
		{
			for (size_t i = 0; i < m_batchDescriptor.sets.size(); i++)
			{
				allocInfo = vks::initializers::descriptorSetAllocateInfo(
					*m_descriptorPool, &m_batchDescriptor.setLayout, 1);
				VK_CHECK_RESULT(vkAllocateDescriptorSets(*m_device, &allocInfo, &m_batchDescriptor.sets[i]));

				writeDescriptorSets = {
					vks::initializers::writeDescriptorSet(
						m_batchDescriptor.sets[i],
						VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
						0,
						m_batchTextures[i].descriptor())
				};

				vkUpdateDescriptorSets(*m_device, static_cast<uint32_t>(writeDescriptorSets.size()),
					writeDescriptorSets.data(), 0, NULL);
			}
		}
	}

	void VulkanLattice::updateLatticeUniformBuffer()
	{
		vkDeviceWaitIdle(*m_device); // Just do it
		m_uniforms.modelview = m_view * m_matrix;
		m_uniforms.normal = glm::transpose(glm::inverse(m_matrix));
		memcpy(m_latticeUniformBuffer.mapped, &m_uniforms, sizeof(m_uniforms));
	}

	void VulkanLattice::updateMatrixUniformBuffer()
	{
		memcpy(m_matrixUniformBuffer.mapped, &m_matrices[0][0], sizeof(glm::mat4) * m_matrices.size());
	}

	void VulkanLattice::setupQueryResultBuffer()
	{
		m_queryResult.count = 10;
		m_timingResult.count = 2;

		uint32_t bufSize = m_queryResult.count * sizeof(uint64_t);

		VkMemoryRequirements memReqs;
		VkMemoryAllocateInfo memAlloc = vks::initializers::memoryAllocateInfo();
		VkBufferCreateInfo bufferCreateInfo =
			vks::initializers::bufferCreateInfo(
				VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
				bufSize);

		// Results are saved in a host visible buffer for easy access by the application
		VK_CHECK_RESULT(vkCreateBuffer(*m_device, &bufferCreateInfo, m_allocator, &m_queryResult.buffer));
		vkGetBufferMemoryRequirements(*m_device, m_queryResult.buffer, &memReqs);
		memAlloc.allocationSize = memReqs.size;
		memAlloc.memoryTypeIndex = m_vulkanDevice->getMemoryType(memReqs.memoryTypeBits,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
		VK_CHECK_RESULT(vkAllocateMemory(*m_device, &memAlloc, m_allocator, &m_queryResult.memory));
		VK_CHECK_RESULT(vkBindBufferMemory(*m_device, m_queryResult.buffer, m_queryResult.memory, 0));

		uint32_t timeBufSize = m_timingResult.count * sizeof(uint64_t);
		VkMemoryRequirements timeMemReqs;
		VkMemoryAllocateInfo timeMemAlloc = vks::initializers::memoryAllocateInfo();
		VkBufferCreateInfo timeBufferCreateInfo =
			vks::initializers::bufferCreateInfo(
				VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
				timeBufSize
			);
		VK_CHECK_RESULT(vkCreateBuffer(*m_device, &timeBufferCreateInfo, m_allocator, &m_timingResult.buffer));
		vkGetBufferMemoryRequirements(*m_device, m_timingResult.buffer, &timeMemReqs);
		timeMemAlloc.allocationSize = timeMemReqs.size;
		timeMemAlloc.memoryTypeIndex = m_vulkanDevice->getMemoryType(timeMemReqs.memoryTypeBits,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
		VK_CHECK_RESULT(vkAllocateMemory(*m_device, &timeMemAlloc, m_allocator, &m_timingResult.memory));
		VK_CHECK_RESULT(vkBindBufferMemory(*m_device, m_timingResult.buffer, m_timingResult.memory, 0));
		m_timestampPeriod = m_vulkanDevice->properties.limits.timestampPeriod;

		VkQueryPoolCreateInfo queryPoolInfo = {};
		queryPoolInfo.sType = VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO;
		queryPoolInfo.queryType = VK_QUERY_TYPE_PIPELINE_STATISTICS;
		queryPoolInfo.pipelineStatistics =
			VK_QUERY_PIPELINE_STATISTIC_INPUT_ASSEMBLY_VERTICES_BIT |
			VK_QUERY_PIPELINE_STATISTIC_INPUT_ASSEMBLY_PRIMITIVES_BIT |
			VK_QUERY_PIPELINE_STATISTIC_VERTEX_SHADER_INVOCATIONS_BIT |
			VK_QUERY_PIPELINE_STATISTIC_GEOMETRY_SHADER_INVOCATIONS_BIT |
			VK_QUERY_PIPELINE_STATISTIC_GEOMETRY_SHADER_PRIMITIVES_BIT |
			VK_QUERY_PIPELINE_STATISTIC_CLIPPING_INVOCATIONS_BIT |
			VK_QUERY_PIPELINE_STATISTIC_CLIPPING_PRIMITIVES_BIT |
			VK_QUERY_PIPELINE_STATISTIC_FRAGMENT_SHADER_INVOCATIONS_BIT |
			VK_QUERY_PIPELINE_STATISTIC_TESSELLATION_CONTROL_SHADER_PATCHES_BIT |
			VK_QUERY_PIPELINE_STATISTIC_TESSELLATION_EVALUATION_SHADER_INVOCATIONS_BIT;

		queryPoolInfo.queryCount = m_queryResult.count;
		VK_CHECK_RESULT(vkCreateQueryPool(*m_device, &queryPoolInfo, NULL, &m_queryPool));

		VkQueryPoolCreateInfo timingPoolInfo = {};
		timingPoolInfo.sType = VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO;
		timingPoolInfo.queryType = VK_QUERY_TYPE_TIMESTAMP;
		timingPoolInfo.queryCount = m_timingResult.count;
		VK_CHECK_RESULT(vkCreateQueryPool(*m_device, &timingPoolInfo, NULL, &m_timingPool));

		m_deviceMemoryUsage += bufSize + timeBufSize;
	}

	void VulkanLattice::getQueryResults()
	{
		if (m_doPipelineQueries) {
			vkGetQueryPoolResults(
				*m_device,
				m_queryPool,
				0,
				1,
				sizeof(m_pipelineStats),
				m_pipelineStats,
				sizeof(uint64_t),
				VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT);
		}

		if (m_doPipelineTimings) {
			vkGetQueryPoolResults(
				*m_device,
				m_timingPool,
				0,
				m_timingResult.count,
				sizeof(m_pipelineTimings),
				m_pipelineTimings,
				sizeof(uint64_t),
				VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT
			);
		}
	}

}