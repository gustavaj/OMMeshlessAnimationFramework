#include "SWVulkanLattice.h"

#include "Shaders.h"
#include <chrono>

namespace SWVL
{

	SWVulkanLattice::SWVulkanLattice()
		: SWVulkanLattice("")
	{
	}

	SWVulkanLattice::SWVulkanLattice(std::string name)
		: OML::Lattice(name), m_pointsPipeline(VK_NULL_HANDLE), m_linesPipeline(VK_NULL_HANDLE), m_localSurfacePipeline(VK_NULL_HANDLE),
		  m_localSurfaceWireframePipeline(VK_NULL_HANDLE), m_patchPipeline(VK_NULL_HANDLE), m_patchWireframePipeline(VK_NULL_HANDLE),
		  m_normalPipeline(VK_NULL_HANDLE), m_pipelineLayout(VK_NULL_HANDLE), m_descriptorSetLayout(VK_NULL_HANDLE), m_descriptorSet(VK_NULL_HANDLE),
		  m_device(nullptr), m_vulkanDevice(nullptr), m_descriptorPool(nullptr), m_renderPass(nullptr),
		  m_queue(nullptr), m_commandPool(nullptr), m_allocator(nullptr), m_selectedSurface(0)
	{
		m_menuSuffix = "##" + std::to_string(OML::Lattice::Index++);
	}

	SWVulkanLattice::~SWVulkanLattice()
	{
	}

	void SWVulkanLattice::initVulkanStuff(VkDevice* device, vks::VulkanDevice* vulkanDevice, VkQueue* queue, VkCommandPool* commandPool, VkDescriptorPool* descriptorPool, VkRenderPass* renderPass, VkAllocationCallbacks* allocator)
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

	void SWVulkanLattice::destroyVulkanStuff()
	{
		if (m_destroyed) return;

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

			vkDestroyPipelineLayout(*m_device, m_pipelineLayout, m_allocator);
			vkDestroyDescriptorSetLayout(*m_device, m_descriptorSetLayout, m_allocator);
		}

		m_latticeUniformBuffer.destroy();
		m_matrixUniformBuffer.destroy();
		m_controlPointBuffer.destroy();
		m_boundariesBuffer.destroy();

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

	void SWVulkanLattice::addToCommandbufferPreRenderpass(VkCommandBuffer& commandBuffer)
	{
		if (m_doPipelineQueries) {
			vkCmdResetQueryPool(commandBuffer, m_queryPool, 0, m_queryResult.count);
		}
		if (m_doPipelineTimings) {
			vkCmdResetQueryPool(commandBuffer, m_timingPool, 0, m_timingResult.count);
		}
	}

	void SWVulkanLattice::addToCommandbuffer(VkCommandBuffer& commandBuffer)
	{
		if (!m_draw) return;

		VkDeviceSize offsets[1] = { 0 };

		if (m_doPipelineQueries) {
			vkCmdBeginQuery(commandBuffer, m_queryPool, 0, 0);
		}
		if (m_doPipelineTimings) {
			vkCmdWriteTimestamp(commandBuffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, m_timingPool, 0);
		}

		vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &m_descriptorSet, 0, NULL);

		if (m_drawSurface)
		{
			if (m_displayPixelAccuracy)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_displayPixelAccuracyPipeline);
			else if (m_wireframe)
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_patchWireframePipeline);
			else
				vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_patchPipeline);
			vkCmdBindVertexBuffers(commandBuffer, 0, 1, &m_patchVertexBuffer.buffer, offsets);
			vkCmdDraw(commandBuffer, m_patchVertexBuffer.count, 1, 0, 0);
			/*for (size_t i = 0; i < m_patchVertices.size() / 4; i++)
			{
				vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &m_descriptorSet, 0, NULL);
				vkCmdDraw(commandBuffer, 4, 1, i * 4, 0);
			}*/

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

		if (m_doPipelineQueries) {
			vkCmdEndQuery(commandBuffer, m_queryPool, 0);
		}
		if (m_doPipelineTimings) {
			vkCmdWriteTimestamp(commandBuffer, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, m_timingPool, 1);
		}
	}

	void SWVulkanLattice::onViewChanged(glm::mat4 projection, glm::mat4 view)
	{
		m_uniforms.projection = projection;
		m_view = view;
		updateLatticeUniformBuffer();
	}

	void SWVulkanLattice::onWindowResized(float width, float height)
	{
		m_uniforms.windowSize = glm::vec2(width, height);
	}

	void SWVulkanLattice::localUpdate(double dt)
	{
		if (m_simulate) {
			updateMatrixUniformBuffer();
		}

		getQueryResults();
	}

	bool SWVulkanLattice::onUpdateUIOverlay(vks::UIOverlay* overlay)
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
				ImGui::Separator();
				if (overlay->checkBox("Pixel-Accurate", &m_drawPixelAccurate)) rebuildCmd = true;
				if (overlay->checkBox("Color by Surf Acc", &m_displaySurfaceAccuracy)) rebuildCmd = true;
				if (overlay->checkBox("Color by Pixel Acc", &m_displayPixelAccuracy)) rebuildCmd = true;
				if (!m_drawPixelAccurate)
				{
					if (overlay->sliderInt("TessInner", &m_uniforms.tessInner, 0, 64)) updateLatticeUniformBuffer();
					if (overlay->sliderInt("TessOuter", &m_uniforms.tessOuter, 0, 64)) updateLatticeUniformBuffer();
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
				if (m_simulatorIndex == 0 || m_simulatorIndex == 1) {
					overlay->sliderFloat("Min amp", &OML::Simulator::MinAmp, OML::Simulator::AmpRange.x, OML::Simulator::MaxAmp);
					overlay->sliderFloat("Max amp", &OML::Simulator::MaxAmp, OML::Simulator::MinAmp, OML::Simulator::AmpRange.y);
				}
				else if (m_simulatorIndex == 2) {
					overlay->sliderFloat("Min angle", &OML::Simulator::MinAngle, OML::Simulator::AngleRange.x, OML::Simulator::MaxAngle);
					overlay->sliderFloat("Max angle", &OML::Simulator::MaxAngle, OML::Simulator::MinAngle, OML::Simulator::AngleRange.y);
				}
				else if (m_simulatorIndex == 3) {
					overlay->sliderFloat("Min scale - 1", &OML::Simulator::MinScale, OML::Simulator::ScaleRange.x, OML::Simulator::MaxScale);
					overlay->sliderFloat("Max scale - 1", &OML::Simulator::MaxScale, OML::Simulator::MinScale, OML::Simulator::ScaleRange.y);
				}
				overlay->sliderFloat("Min speed", &OML::Simulator::MinSpeed, OML::Simulator::SpeedRange.x, OML::Simulator::MaxSpeed);
				overlay->sliderFloat("Max speed", &OML::Simulator::MaxSpeed, OML::Simulator::MinSpeed, OML::Simulator::SpeedRange.y);
				bool simulatorActive = m_simulators.find(static_cast<OML::SimulatorTypes>(m_simulatorIndex)) != m_simulators.end();
				std::string addBtnTitle = simulatorActive ? "Update" : "Add";
				if (overlay->button(addBtnTitle.c_str())) {
					if (m_simulatorIndex == 0) addNormalSinSimulation();
					else if (m_simulatorIndex == 1) addRandomSphereSimulation();
					else if (m_simulatorIndex == 2) addNormalRotationSimulation();
					else if (m_simulatorIndex == 3) addXYScalingSimulation();
				}
				if (simulatorActive) {
					ImGui::SameLine();
					if (overlay->button("Remove")) {
						if (m_simulatorIndex == 0) removeSimulator(OML::SimulatorTypes::NormalSin);
						else if (m_simulatorIndex == 1) removeSimulator(OML::SimulatorTypes::RandomSphere);
						else if (m_simulatorIndex == 2) removeSimulator(OML::SimulatorTypes::Rotation);
						else if (m_simulatorIndex == 3) removeSimulator(OML::SimulatorTypes::XYScale);
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

	void SWVulkanLattice::CheckAndSetupRequiredPhysicalDeviceFeatures(
		VkPhysicalDeviceFeatures& deviceFeatures, VkPhysicalDeviceFeatures& enabledFeatures)
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

		vkDestroyBuffer(*m_device, stagingBuffer.buffer, m_allocator);
		vkFreeMemory(*m_device, stagingBuffer.memory, m_allocator);
	}

	void SWVulkanLattice::setupVertices()
	{
		// Set up local surface vertices
		m_localSurfaceVertices.resize(0);
		for (size_t i = 0; i < m_loci.size(); i++)
		{
			m_localSurfaceVertices.push_back(LocalSurfaceVertex(
				m_loci[i].controlPointIndex, m_loci[i].controlPointCount, m_loci[i].matrixIndex, 0));
			m_localSurfaceVertices.back().color = m_loci[i].color;
		}

		// Set up patch vertices
		m_patchVertices.resize(0);
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
			m_patchVertices.push_back(localVert00);

			LocalSurfaceVertex localVert10;
			localVert10.controlPointIndex = locus10.controlPointIndex;
			localVert10.controlPointCount = locus10.controlPointCount;
			localVert10.matrixIndex = locus10.matrixIndex;
			localVert10.boundaryIndex = locus10.boundaryIndices[m_patches[i].faceIdx];
			localVert10.color = m_patches[i].color;
			m_patchVertices.push_back(localVert10);

			LocalSurfaceVertex localVert01;
			localVert01.controlPointIndex = locus01.controlPointIndex;
			localVert01.controlPointCount = locus01.controlPointCount;
			localVert01.matrixIndex = locus01.matrixIndex;
			localVert01.boundaryIndex = locus01.boundaryIndices[m_patches[i].faceIdx];
			localVert01.color = m_patches[i].color;
			m_patchVertices.push_back(localVert01);

			LocalSurfaceVertex localVert11;
			localVert11.controlPointIndex = locus11.controlPointIndex;
			localVert11.controlPointCount = locus11.controlPointCount;
			localVert11.matrixIndex = locus11.matrixIndex;
			localVert11.boundaryIndex = locus11.boundaryIndices[m_patches[i].faceIdx];
			localVert11.color = m_patches[i].color;
			m_patchVertices.push_back(localVert11);
		}
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
			VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&m_matrixUniformBuffer,
			sizeof(glm::mat4) * m_matrices.size()
		));

		VK_CHECK_RESULT(m_matrixUniformBuffer.map());

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

		updateLatticeUniformBuffer();
		updateMatrixUniformBuffer();
		uploadStorageBuffers();
	}

	void SWVulkanLattice::setupDescriptorSetLayouts()
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
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_FRAGMENT_BIT,
				1),
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_FRAGMENT_BIT,
				2),
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT | VK_SHADER_STAGE_FRAGMENT_BIT,
				3)
		};

		VkDescriptorSetLayoutCreateInfo descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(
			setLayoutBindings.data(), static_cast<uint32_t>(setLayoutBindings.size()));
		VK_CHECK_RESULT(vkCreateDescriptorSetLayout(*m_device, &descriptorLayout,
			m_allocator, &m_descriptorSetLayout));
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
		shaderStage.stage = stage;
		if (m_shaderModules.find(fileName) == m_shaderModules.end()) {
			m_shaderModules.insert({ fileName, vks::tools::loadShader(fileName.c_str(), *m_device) });
		}
		shaderStage.module = m_shaderModules[fileName];
		shaderStage.pName = "main";
		assert(shaderStage.module != VK_NULL_HANDLE);
		return shaderStage;
	}

	VkPipelineShaderStageCreateInfo SWVulkanLattice::loadShader(std::pair<std::string, std::vector<uint32_t>&> src, VkShaderStageFlagBits stage)
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

	void SWVulkanLattice::preparePipelines()
	{
		auto start = std::chrono::high_resolution_clock::now();
		std::cout << "Preparing pipelines" << std::endl;

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
		options.numSamplesU = 0;
		options.numSamplesV = 0;
		options.maxError = 1.0f;
		options.normalLength = 10.0f;

		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier;
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct;

		std::array<VkPipelineShaderStageCreateInfo, 4> localShaderStages;
		localShaderStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		localShaderStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(1), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		localShaderStages[2] = loadShader(OML::Shaders::GetTeseShader(lsType, OML::TeseShaderType::Local, evalMethod, options), 
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
		bsShaderStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(4), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		bsShaderStages[2] = loadShader(OML::Shaders::GetTeseShader(lsType, OML::TeseShaderType::Lattice, evalMethod, options), 
			VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		bsShaderStages[3] = loadShader(OML::Shaders::GetShadedColorFragShader(), VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(bsShaderStages.size());
		pipelineCreateInfo.pStages = bsShaderStages.data();

		rasterizationState.polygonMode = VK_POLYGON_MODE_LINE;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_patchWireframePipeline));

		rasterizationState.polygonMode = VK_POLYGON_MODE_FILL;
		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_patchPipeline));
		
		std::array<VkPipelineShaderStageCreateInfo, 4> pixelAccuracyStages;
		pixelAccuracyStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		pixelAccuracyStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(4), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		pixelAccuracyStages[2] = loadShader(OML::Shaders::GetTeseShader(lsType, OML::TeseShaderType::Pixel_Accuracy, evalMethod, options), 
			VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		pixelAccuracyStages[3] = loadShader(OML::Shaders::GetBiQuadLatticePixelAccuracyFragShader(
			options.numControl, options.numLocal, options.numPatches), VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(pixelAccuracyStages.size());
		pipelineCreateInfo.pStages = pixelAccuracyStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_displayPixelAccuracyPipeline));

		std::array<VkPipelineShaderStageCreateInfo, 5> normalShaderStages;
		normalShaderStages[0] = loadShader(OML::Shaders::GetLocalSurfaceInfoVertShader(), VK_SHADER_STAGE_VERTEX_BIT);
		normalShaderStages[1] = loadShader(OML::Shaders::GetLocalSurfaceInfoTescShader(4), VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT);
		normalShaderStages[2] = loadShader(OML::Shaders::GetTeseShader(lsType, OML::TeseShaderType::Normals, evalMethod, options), 
			VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT);
		normalShaderStages[3] = loadShader(OML::Shaders::GetLatticeNormalsGeomShader(8.0f), VK_SHADER_STAGE_GEOMETRY_BIT);
		normalShaderStages[4] = loadShader(OML::Shaders::GetFlatColorFragShader(), VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(normalShaderStages.size());
		pipelineCreateInfo.pStages = normalShaderStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(*m_device, /*pipelineCache*/nullptr, 1, &pipelineCreateInfo, m_allocator, &m_normalPipeline));

		auto end = std::chrono::high_resolution_clock::now();
		auto time = std::chrono::duration<double, std::milli>(end - start).count();

		std::cout << "Preparing pipelines time: " << time << "ms" << std::endl;
	}

	void SWVulkanLattice::setupDescriptorPool()
	{
		std::vector<VkDescriptorPoolSize> poolSizes =
		{
			vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1),
			vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 3)
		};

		VkDescriptorPoolCreateInfo descriptorPoolInfo =
			vks::initializers::descriptorPoolCreateInfo(
				static_cast<uint32_t>(poolSizes.size()),
				poolSizes.data(),
				1);

		VK_CHECK_RESULT(vkCreateDescriptorPool(*m_device, &descriptorPoolInfo, m_allocator, m_descriptorPool));
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

		vkUpdateDescriptorSets(*m_device, static_cast<uint32_t>(writeDescriptorSets.size()), writeDescriptorSets.data(), 0, NULL);
	}

	void SWVulkanLattice::updateLatticeUniformBuffer()
	{
		vkDeviceWaitIdle(*m_device); // Just do it
		m_uniforms.modelview = m_view * m_matrix;
		m_uniforms.normal = glm::transpose(glm::inverse(m_matrix));
		memcpy(m_latticeUniformBuffer.mapped, &m_uniforms, sizeof(m_uniforms));
	}

	void SWVulkanLattice::updateMatrixUniformBuffer()
	{
		memcpy(m_matrixUniformBuffer.mapped, &m_matrices[0][0], sizeof(glm::mat4) * m_matrices.size());
	}

	void SWVulkanLattice::uploadStorageBuffers()
	{
		createDeviceLocalBuffer(m_controlPointBuffer.buffer, m_controlPointBuffer.memory,
			m_controlPoints.data(), m_controlPointBuffer.size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
		m_controlPointBuffer.descriptor.buffer = m_controlPointBuffer.buffer;
		createDeviceLocalBuffer(m_boundariesBuffer.buffer, m_boundariesBuffer.memory,
			m_boundaries.data(), m_boundariesBuffer.size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
		m_boundariesBuffer.descriptor.buffer = m_boundariesBuffer.buffer;
	}

	void SWVulkanLattice::setupQueryResultBuffer()
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
	}

	void SWVulkanLattice::getQueryResults()
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