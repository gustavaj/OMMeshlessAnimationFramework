#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <algorithm>

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vulkan/vulkan.h>

#include "vulkan/VulkanBaseExample.h"
#include "vulkan/VulkanBuffer.hpp"
#include "vulkan/VulkanModel.hpp"

#include "lattice/Lattice.h"

class LatticeExample : public VulkanExampleBase
{
public:

	std::vector<OML::Lattice> lattices;

	struct {
		VkBuffer buffer;
		VkDeviceMemory memory;
	} queryResult, timingResult;

	VkQueryPool queryPool = VK_NULL_HANDLE;
	VkQueryPool timingPool = VK_NULL_HANDLE;
	int pipelineStatsCount = 10;
	uint64_t pipelineStats[10] = { 0 };
	int timepointStatsCount = 2;
	uint64_t timepointStats[2] = { 0 };
	// VkPhysicalDeviceLimits::timestampPeriod (the number of nanoseconds it takes for a timestamp value to be incremented by 1)
	double timestampPeriod;


	LatticeExample(bool enableValidation)
		: VulkanExampleBase(enableValidation)
	{
		zoom = -100.0f;
		zoomSpeed = 20.0f;
		title = "OpenMesh Lattice Example";
		camera.setPerspective(60.0f, (float)width / (float)height, 1.0f, 10000.0f);
		//camera.position = glm::vec3(0.0f, 0.0f, -30.0f);
		/*camera.type = Camera::CameraType::firstperson;*/
		camera.setRotation(glm::vec3(0.0f, 0.0f, 0.0f));
		camera.setTranslation(glm::vec3(0.0f, 0.0f, -50.0f));
		/*camera.movementSpeed = 0.025f;*/
		settings.overlay = true;
	}

	~LatticeExample()
	{
		for (auto& lat : lattices)
			lat.destroyVulkanStuff();

		if (queryPool != VK_NULL_HANDLE) {
			vkDestroyQueryPool(device, queryPool, nullptr);
			vkDestroyBuffer(device, queryResult.buffer, nullptr);
			vkFreeMemory(device, queryResult.memory, nullptr);
		}
		if (timingPool != VK_NULL_HANDLE) {
			vkDestroyQueryPool(device, timingPool, nullptr);
			vkDestroyBuffer(device, timingResult.buffer, nullptr);
			vkFreeMemory(device, timingResult.memory, nullptr);
		}
	}

	virtual void getEnabledFeatures() {
		// Tessellation shader support is required for this example
		if (deviceFeatures.tessellationShader && deviceFeatures.geometryShader &&
			deviceFeatures.fillModeNonSolid && deviceFeatures.pipelineStatisticsQuery) {
			enabledFeatures.tessellationShader = VK_TRUE;
			enabledFeatures.geometryShader = VK_TRUE;
			enabledFeatures.fillModeNonSolid = VK_TRUE;
			enabledFeatures.pipelineStatisticsQuery = VK_TRUE;
		}
		else {
			vks::tools::exitFatal("Selected GPU does not support tessellation shaders!", VK_ERROR_FEATURE_NOT_PRESENT);
		}
	}

	void setupQueryResultBuffer() {
		uint32_t bufSize = pipelineStatsCount * sizeof(uint64_t);

		VkMemoryRequirements memReqs;
		VkMemoryAllocateInfo memAlloc = vks::initializers::memoryAllocateInfo();
		VkBufferCreateInfo bufferCreateInfo =
			vks::initializers::bufferCreateInfo(
				VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
				bufSize);

		// Results are saved in a host visible buffer for easy access by the application
		VK_CHECK_RESULT(vkCreateBuffer(device, &bufferCreateInfo, nullptr, &queryResult.buffer));
		vkGetBufferMemoryRequirements(device, queryResult.buffer, &memReqs);
		memAlloc.allocationSize = memReqs.size;
		memAlloc.memoryTypeIndex = vulkanDevice->getMemoryType(memReqs.memoryTypeBits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
		VK_CHECK_RESULT(vkAllocateMemory(device, &memAlloc, nullptr, &queryResult.memory));
		VK_CHECK_RESULT(vkBindBufferMemory(device, queryResult.buffer, queryResult.memory, 0));

		uint32_t timeBufSize = timepointStatsCount * sizeof(uint64_t);
		VkMemoryRequirements timeMemReqs;
		VkMemoryAllocateInfo timeMemAlloc = vks::initializers::memoryAllocateInfo();
		VkBufferCreateInfo timeBufferCreateInfo =
			vks::initializers::bufferCreateInfo(
				VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
				timeBufSize
			);
		VK_CHECK_RESULT(vkCreateBuffer(device, &timeBufferCreateInfo, nullptr, &timingResult.buffer));
		vkGetBufferMemoryRequirements(device, timingResult.buffer, &timeMemReqs);
		timeMemAlloc.allocationSize = timeMemReqs.size;
		timeMemAlloc.memoryTypeIndex = vulkanDevice->getMemoryType(timeMemReqs.memoryTypeBits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
		VK_CHECK_RESULT(vkAllocateMemory(device, &timeMemAlloc, nullptr, &timingResult.memory));
		VK_CHECK_RESULT(vkBindBufferMemory(device, timingResult.buffer, timingResult.memory, 0));
		const auto timeperiod = vulkanDevice->properties.limits.timestampPeriod;
		timestampPeriod = timeperiod;

		// Create query pool
		if (deviceFeatures.pipelineStatisticsQuery) {
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

			queryPoolInfo.queryCount = pipelineStatsCount;
			VK_CHECK_RESULT(vkCreateQueryPool(device, &queryPoolInfo, NULL, &queryPool));

			VkQueryPoolCreateInfo timingPoolInfo = {};
			timingPoolInfo.sType = VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO;
			timingPoolInfo.queryType = VK_QUERY_TYPE_TIMESTAMP;
			timingPoolInfo.queryCount = timepointStatsCount;
			VK_CHECK_RESULT(vkCreateQueryPool(device, &timingPoolInfo, NULL, &timingPool));
		}
	}

	void getQueryResults() {
		// We use vkGetQueryResults to copy the results into a host visible buffer
		vkGetQueryPoolResults(
			device,
			queryPool,
			0,
			1,
			sizeof(pipelineStats),
			pipelineStats,
			sizeof(uint64_t),
			VK_QUERY_RESULT_64_BIT);
		vkGetQueryPoolResults(
			device,
			timingPool,
			0,
			timepointStatsCount,
			sizeof(timepointStats),
			timepointStats,
			sizeof(uint64_t),
			VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT
		);
	}

	void loadAssets() {
		createLatticeGeometry();

		for (auto& lat : lattices)
		{
			lat.induceLattice();
			lat.initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
		}
	}

	virtual void createLatticeGeometry()
	{
		lattices.push_back(OML::Lattice("Test Lattice"));
		lattices.back().addPatch(OML::Vec2f(-5.0f, -5.0f), 10.0f, 10.0f);
	}

	void reBuildCommangBuffers() {
		if (!checkCommandBuffers())
		{
			destroyCommandBuffers();
			createCommandBuffers();
		}
		buildCommandBuffers();
	}

	void buildCommandBuffers() {
		VkCommandBufferBeginInfo cmdBufInfo = vks::initializers::commandBufferBeginInfo();

		VkClearValue clearValues[2];
		clearValues[0].color = { 0.2f, 0.67f, 0.92f, 1.0f };
		clearValues[1].depthStencil = { 1.0f, 0 };

		VkRenderPassBeginInfo renderPassBeginInfo = vks::initializers::renderPassBeginInfo();
		renderPassBeginInfo.renderPass = renderPass;
		renderPassBeginInfo.renderArea.offset.x = 0;
		renderPassBeginInfo.renderArea.offset.y = 0;
		renderPassBeginInfo.renderArea.extent.width = width;
		renderPassBeginInfo.renderArea.extent.height = height;
		renderPassBeginInfo.clearValueCount = 2;
		renderPassBeginInfo.pClearValues = clearValues;

		for (int32_t i = 0; i < drawCmdBuffers.size(); ++i)
		{
			renderPassBeginInfo.framebuffer = frameBuffers[i];

			VK_CHECK_RESULT(vkBeginCommandBuffer(drawCmdBuffers[i], &cmdBufInfo));
			
			vkCmdResetQueryPool(drawCmdBuffers[i], queryPool, 0, pipelineStatsCount);
			vkCmdResetQueryPool(drawCmdBuffers[i], timingPool, 0, timepointStatsCount);
			vkCmdWriteTimestamp(drawCmdBuffers[i], VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, timingPool, 0);

			vkCmdBeginRenderPass(drawCmdBuffers[i], &renderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

			VkViewport viewport = vks::initializers::viewport((float)width, (float)height, 0.0f, 1.0f);
			vkCmdSetViewport(drawCmdBuffers[i], 0, 1, &viewport);

			VkRect2D scissor = vks::initializers::rect2D(width, height, 0, 0);
			vkCmdSetScissor(drawCmdBuffers[i], 0, 1, &scissor);

			vkCmdSetLineWidth(drawCmdBuffers[i], 1.0f);

			vkCmdBeginQuery(drawCmdBuffers[i], queryPool, 0, 0);

			for (auto& lat : lattices)
			{
				lat.addToCommandbuffer(drawCmdBuffers[i]);
			}

			vkCmdEndQuery(drawCmdBuffers[i], queryPool, 0);

			drawUI(drawCmdBuffers[i]);

			vkCmdEndRenderPass(drawCmdBuffers[i]);

			vkCmdWriteTimestamp(drawCmdBuffers[i], VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, timingPool, 1);

			VK_CHECK_RESULT(vkEndCommandBuffer(drawCmdBuffers[i]));
		}
	}

	void draw() {
		VulkanExampleBase::prepareFrame();

		// Command buffer to be sumitted to the queue
		submitInfo.commandBufferCount = 1;
		submitInfo.pCommandBuffers = &drawCmdBuffers[currentBuffer];

		// Submit to queue
		VK_CHECK_RESULT(vkQueueSubmit(queue, 1, &submitInfo, VK_NULL_HANDLE));

		if (deviceFeatures.pipelineStatisticsQuery) {
			// Read query results for displaying in next frame
			getQueryResults();
		}

		VulkanExampleBase::submitFrame();
	}

	void prepare() {
		VulkanExampleBase::prepare();
		loadAssets();
		if (deviceFeatures.pipelineStatisticsQuery) {
			setupQueryResultBuffer();
		}
		buildCommandBuffers();
		viewChanged();
		prepared = true;
	}

	virtual void render() {
		if (!prepared)
			return;
		for (auto& lat : lattices)
			lat.update(frameTimer);
		draw();
	}

	virtual void viewChanged() {
		for (auto& lat : lattices)
		{
			lat.onViewChanged(camera.matrices.perspective, camera.matrices.view);
		}
	}

	virtual void OnUpdateUIOverlay(vks::UIOverlay* overlay) {

		if (overlay->header("Pipeline Stats", false)) {
			overlay->text("Input Assembler Vertices: %d", pipelineStats[0]);
			overlay->text("Input Assembler Primitives: %d", pipelineStats[1]);
			overlay->text("Vert invocations: %d", pipelineStats[2]);
			overlay->text("TCS patches: %d", pipelineStats[8]);
			overlay->text("TES invocations: %d", pipelineStats[9]);
			overlay->text("Geom invocations: %d", pipelineStats[3]);
			overlay->text("Geom primitives: %d", pipelineStats[4]);
			overlay->text("Clipping invocations: %d", pipelineStats[5]);
			overlay->text("Clipping primitives: %d", pipelineStats[6]);
			overlay->text("Frag invocations: %d", pipelineStats[7]);
		}

		if (overlay->header("Timings", false)) {
			const uint64_t diff = timepointStats[1] - timepointStats[0];
			double nanos1 = diff * timestampPeriod / 1000000;
			if (nanos1 > 100000.0f) nanos1 = 0.0f;
			overlay->text("Renderpass time: %.4f ms", nanos1);
		}

		for (auto& lat : lattices)
		{
			if (lat.onUpdateUIOverlay(overlay)) {
				buildCommandBuffers();
			}
		}		
	}

};

class GridLatticeExample : public LatticeExample
{
public:
	int m_rows, m_cols;
	float m_width, m_height;

	GridLatticeExample(bool enableValidation, float width, float height, int rows, int cols)
		: LatticeExample(enableValidation), m_width(width), m_height(height), m_rows(rows), m_cols(cols)
	{
	}

	virtual void createLatticeGeometry() override
	{
		lattices.push_back(OML::Lattice("Grid Lattice"));
		lattices.back().addGrid(OML::Vec2f(-m_width / 2, -m_height / 2), m_width, m_height, m_rows, m_cols);
	}
};

class CylinderLatticeExample : public LatticeExample
{
public:
	int m_rows, m_cols;
	float m_radius, m_height;

	CylinderLatticeExample(bool enableValidation, float radius, float height, int rows, int cols)
		: LatticeExample(enableValidation), m_radius(radius), m_height(height), m_rows(rows), m_cols(cols)
	{
	}

	virtual void createLatticeGeometry() override
	{
		lattices.push_back(OML::Lattice("Cylinder Lattice"));
		lattices.back().addCylinder(OML::Vec3f(0.0f, 0.0f, 0.0f), m_radius, m_height, m_rows, m_cols);
	}
};

class SphereLatticeExample : public LatticeExample
{
public:
	int m_segments, m_slices;
	float m_radius;

	SphereLatticeExample(bool enableValidation, float radius, int segments, int slices)
		: LatticeExample(enableValidation), m_radius(radius), m_segments(segments), m_slices(slices)
	{
	}

	virtual void createLatticeGeometry() override
	{
		lattices.push_back(OML::Lattice("Sphere Lattice"));
		lattices.back().addSphere(OML::Vec3f(0.0f, 0.0f, 0.0f), m_radius, m_segments, m_slices);
	}
};

class NonUniformGridExample : public LatticeExample
{
public:
	NonUniformGridExample(bool enableValidation)
		: LatticeExample(enableValidation)
	{
	}

	virtual void createLatticeGeometry() override
	{
		lattices.push_back(OML::Lattice("Non uniform Lattice"));
		OML::Vec3f p00(-30, -20, 0), p10(0, -20, 0), p20(10, -20, 0), p30(20, -20, 0);
		OML::Vec3f p01(-30, 0, 0),	p11(0, 0, 0),	p21(10, 0, 0), p31(20, 0, 0);
		OML::Vec3f p02(-30, 10, 0),	p12(0, 10, 0),	p22(10, 10, 0), p32(20, 10, 0);
		OML::Vec3f p03(-30, 20, 0),	p13(0, 20, 0),	p23(10, 20, 0), p33(20, 20, 0);
		lattices[0].addPatch(p00, p10, p01, p11);
		lattices[0].addPatch(p10, p20, p11, p21);
		lattices[0].addPatch(p20, p30, p21, p31);
		lattices[0].addPatch(p01, p11, p02, p12);
		lattices[0].addPatch(p11, p21, p12, p22);
		lattices[0].addPatch(p21, p31, p22, p32);
		lattices[0].addPatch(p02, p12, p03, p13);
		lattices[0].addPatch(p12, p22, p13, p23);
		lattices[0].addPatch(p22, p32, p23, p33);
	}
};

class NonRectangularExample : public LatticeExample
{
public:
	NonRectangularExample(bool enableValidation)
		: LatticeExample(enableValidation)
	{
	}

	virtual void createLatticeGeometry() override
	{
		lattices.push_back(OML::Lattice("Non-rectangular Lattice"));
		// Grid where quads have angles != 90
		lattices[0].addPatch(OML::Vec2f(-6,-15), OML::Vec2f(0,-10), OML::Vec2f(-15,-5), OML::Vec2f(0,0));
		lattices[0].addPatch(OML::Vec2f(0,-10), OML::Vec2f(5,-6), OML::Vec2f(0,0), OML::Vec2f(10,0));
		lattices[0].addPatch(OML::Vec2f(5, -6), OML::Vec2f(17,-5), OML::Vec2f(10,0), OML::Vec2f(14,-1));
		lattices[0].addPatch(OML::Vec2f(-15,-5), OML::Vec2f(0,0), OML::Vec2f(-5,13), OML::Vec2f(2,10));
		lattices[0].addPatch(OML::Vec2f(0,0), OML::Vec2f(10,0), OML::Vec2f(2,10), OML::Vec2f(10,10));
		lattices[0].addPatch(OML::Vec2f(10, 0), OML::Vec2f(14, -1), OML::Vec2f(10, 10), OML::Vec2f(20, 10));
		lattices[0].addPatch(OML::Vec2f(-5, 13), OML::Vec2f(2, 10), OML::Vec2f(-8, 20), OML::Vec2f(2, 20));
		lattices[0].addPatch(OML::Vec2f(2, 10), OML::Vec2f(10, 10), OML::Vec2f(2, 20), OML::Vec2f(10, 15));
		lattices[0].addPatch(OML::Vec2f(10, 10), OML::Vec2f(20, 10), OML::Vec2f(10, 15), OML::Vec2f(15, 18));
	}
};

class TLocusExample : public LatticeExample
{
public:
	TLocusExample(bool enableValidation)
		: LatticeExample(enableValidation)
	{
	}

	virtual void createLatticeGeometry() override
	{
		lattices.push_back(OML::Lattice("T-locus Lattice"));

		lattices[0].addPatch(OML::Vec2f(-20.0f, -20.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(-10.0f, -20.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(0.0f, -20.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(10.0f, -20.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(-20.0f, -10.0f), 10.0f, 20.0f);

		OML::Vec2f p1(-10.0f, -10.0f);
		OML::Vec2f p2(0.0f, -10.0f);
		OML::Vec2f p3(0.0f, 0.0f);
		lattices[0].addPatch(p1, 10.0f, 20.0f);
		lattices[0].addPatch(p2, 10.0f, 10.0f);
		lattices[0].addPatch(p3, 10.0f, 10.0f);

		lattices[0].addPatch(OML::Vec2f(10.0f, -10.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(10.0f, 0.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(-20.0f, 10.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(-10.0f, 10.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(0.0f, 10.0f), 10.0f, 10.0f);
		lattices[0].addPatch(OML::Vec2f(10.0f, 10.0f), 10.0f, 10.0f);
	}
};




class MultiLatticeExample : public LatticeExample
{
public:
	MultiLatticeExample(bool enableValidation) : LatticeExample(enableValidation) {}

	virtual void createLatticeGeometry() override
	{

		lattices.push_back(OML::Lattice("Grid Lattice"));
		lattices[0].addGrid(OML::Vec2f(-25.0f, -25.0f), 20.0f, 20.0f, 3, 3);
		lattices[0].induceLattice();
		lattices[0].initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);

		lattices.push_back(OML::Lattice("Cylinder Lattice"));
		lattices[1].addCylinder(OML::Vec3f(20.0f, 0.0f, -30.0f), 10.0f, 30.0f, 4, 8);
		lattices[1].induceLattice();
		lattices[1].initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);

		lattices.push_back(OML::Lattice("Sphere Lattice"));
		lattices[2].addSphere(OML::Vec3f(0.0f, -20.0f, 20.0f), 8.0f, 6, 6);
		lattices[2].induceLattice();
		lattices[2].initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
	}
};