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

#include "lattice/SWVulkanLattice.h"

class LatticeExample : public VulkanExampleBase
{
public:

	std::vector<SWVL::SWVulkanLattice> lattices;

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
	}

	virtual void getEnabledFeatures() {
		// Tessellation shader support is required for this example
		if (deviceFeatures.tessellationShader && deviceFeatures.geometryShader &&
			deviceFeatures.fillModeNonSolid && deviceFeatures.pipelineStatisticsQuery &&
			deviceFeatures.vertexPipelineStoresAndAtomics) {
			enabledFeatures.tessellationShader = VK_TRUE;
			enabledFeatures.geometryShader = VK_TRUE;
			enabledFeatures.fillModeNonSolid = VK_TRUE;
			enabledFeatures.pipelineStatisticsQuery = VK_TRUE;
			enabledFeatures.vertexPipelineStoresAndAtomics = VK_TRUE;
		}
		else {
			vks::tools::exitFatal("Selected GPU does not support tessellation shaders!", VK_ERROR_FEATURE_NOT_PRESENT);
		}
	}

	void setupQueryResultBuffer() {
		
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
		lattices.push_back(SWVL::SWVulkanLattice("Test Lattice"));
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

			for (auto& lat : lattices)
			{
				lat.addToCommandbufferPreRenderpass(drawCmdBuffers[i]);
			}

			vkCmdBeginRenderPass(drawCmdBuffers[i], &renderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

			VkViewport viewport = vks::initializers::viewport((float)width, (float)height, 0.0f, 1.0f);
			vkCmdSetViewport(drawCmdBuffers[i], 0, 1, &viewport);

			VkRect2D scissor = vks::initializers::rect2D(width, height, 0, 0);
			vkCmdSetScissor(drawCmdBuffers[i], 0, 1, &scissor);

			vkCmdSetLineWidth(drawCmdBuffers[i], 1.0f);

			for (auto& lat : lattices)
			{
				lat.addToCommandbuffer(drawCmdBuffers[i]);
			}

			drawUI(drawCmdBuffers[i]);

			vkCmdEndRenderPass(drawCmdBuffers[i]);

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
		lattices.push_back(SWVL::SWVulkanLattice("Grid Lattice"));
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
		lattices.push_back(SWVL::SWVulkanLattice("Cylinder Lattice"));
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
		lattices.push_back(SWVL::SWVulkanLattice("Sphere Lattice"));
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
		lattices.push_back(SWVL::SWVulkanLattice("Non uniform Lattice"));
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
		lattices.push_back(SWVL::SWVulkanLattice("Non-rectangular Lattice"));
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
		lattices.push_back(SWVL::SWVulkanLattice("T-locus Lattice"));

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

		lattices.push_back(SWVL::SWVulkanLattice("Grid Lattice"));
		lattices[0].addGrid(OML::Vec2f(-25.0f, -25.0f), 20.0f, 20.0f, 3, 3);
		lattices[0].induceLattice();
		lattices[0].initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);

		lattices.push_back(SWVL::SWVulkanLattice("Cylinder Lattice"));
		lattices[1].addCylinder(OML::Vec3f(20.0f, 0.0f, -30.0f), 10.0f, 30.0f, 4, 8);
		lattices[1].induceLattice();
		lattices[1].initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);

		lattices.push_back(SWVL::SWVulkanLattice("Sphere Lattice"));
		lattices[2].addSphere(OML::Vec3f(0.0f, -20.0f, 20.0f), 8.0f, 6, 6);
		lattices[2].induceLattice();
		lattices[2].initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
	}
};