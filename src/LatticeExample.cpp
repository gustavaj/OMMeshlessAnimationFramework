#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <sstream>

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

const std::vector<std::string> LATTICE_TYPES = { "Grid", "Cylinder", "Sphere", "Patches" };

using Lattice = SWVL::SWVulkanLattice;

class LatticeExample : public VulkanExampleBase
{
public:

	std::vector<Lattice> lattices;

	// Variables used for lattice creation
	int latticeCreateTypeIndex = 0;
	int c_rows = 2, c_cols = 2;
	float c_width = 10.0f, c_height = 10.0f;
	glm::vec3 c_pos = glm::vec3(0.0f);
	glm::vec3 c_rot = glm::vec3(0.0f);
	glm::vec3 c_col = glm::vec3(0.8f, 0.2f, 0.5f);
	bool perPatchColors = false;
	int latticeIdx = 0;
	int latticeDeleteIdx = 0;
	char latticeName[64] = "Lattice";
	std::vector<std::string> latticeNames;

	float patchTopLeftX = 0.0f, patchTopLeftY = 0.0f;
	float patchWidth = 10.0f, patchHeight = 10.0f;
	int patchInputIdx = 0;
	std::vector<OML::Vec2f> patchInputs;
	std::vector<std::string> patchInputStrings;
	std::vector<const char*> patchInputChars;

	LatticeExample(bool enableValidation)
		: VulkanExampleBase(enableValidation)
	{
		title = "Lattice Example";

		// Camera
		zoom = -100.0f;
		zoomSpeed = 20.0f;
		camera.setPerspective(60.0f, (float)width / (float)height, 1.0f, 10000.0f);
		/*camera.type = Camera::CameraType::firstperson;*/
		/*camera.movementSpeed = 0.025f;*/
		camera.setRotation(glm::vec3(0.0f, 0.0f, 0.0f));
		camera.setTranslation(glm::vec3(0.0f, 0.0f, -50.0f));

		settings.overlay = true; // ImGui overlay
	}

	~LatticeExample()
	{
		for (auto& lat : lattices)
			lat.destroyVulkanStuff();
	}

	void setupAnimation()
	{
		for (auto& lat : lattices)
		{
			lat.addNormalSinSimulation();
			lat.setAnimate(true);
			lat.setDrawLatticeGrid(false);
		}
	}

	virtual void getEnabledFeatures() {
		Lattice::CheckAndSetupRequiredPhysicalDeviceFeatures(deviceFeatures, enabledFeatures);
	}

	void loadAssets() {
		createLatticeGeometry();
	}

	virtual void createLatticeGeometry()
	{
	}

	void buildCommandBuffers() {
		vkDeviceWaitIdle(device);

		if (!checkCommandBuffers())
		{
			destroyCommandBuffers();
			createCommandBuffers();
		}

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

		// Submit to queue
		VK_CHECK_RESULT(vkQueueSubmit(queue, 1, &submitInfo, inFlightFences[currentFrame]));

		VulkanExampleBase::submitFrame();
	}

	void prepare() {
		VulkanExampleBase::prepare();
		loadAssets();
		if(!checkCommandBuffers())
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

	void addLattice(Lattice& lattice) 
	{
		lattice.setName(lattice.name() + "##" + std::to_string(latticeIdx));
		lattice.induceLattice();
		lattice.initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
		lattice.onViewChanged(camera.matrices.perspective, camera.matrices.view);
		lattices.push_back(std::move(lattice));
		latticeNames.push_back(lattice.name());
		latticeIdx++;
		buildCommandBuffers();
	}

	void removeLattice(std::string name)
	{
		vkQueueWaitIdle(queue);
		for (auto lat_it = lattices.begin(); lat_it != lattices.end(); lat_it++)
		{
			if (lat_it->name() == name)
			{
				lat_it->destroyVulkanStuff();
				lattices.erase(lat_it);
				break;
			}
		}
		for (auto lat_name_it = latticeNames.begin(); lat_name_it != latticeNames.end(); lat_name_it++)
		{
			if ((*lat_name_it) == name)
			{
				latticeNames.erase(lat_name_it);
				break;
			}
		}
	}

	virtual void OnUpdateUIOverlay(vks::UIOverlay* overlay) {
		// Lattice creation
		if (overlay->header("Add/Delete", false))
		{
			if (latticeCreateTypeIndex == 3)
			{
				ImGui::PushItemWidth(260.0f * UIOverlay.scale);
				ImGui::ListBox("", &patchInputIdx, patchInputChars.data(), patchInputChars.size(), 4);
				ImGui::PushItemWidth(110.0f * UIOverlay.scale);
			}

			ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(5.0f, 1.0f));

			ImGui::BeginGroup();

			overlay->comboBox("Type", &latticeCreateTypeIndex, LATTICE_TYPES);
			ImGui::InputText("Name", latticeName, IM_ARRAYSIZE(latticeName));

			if (latticeCreateTypeIndex == 0)
			{
				overlay->sliderFloat("width", &c_width, 10.0f, 200.0f);
				overlay->sliderFloat("height", &c_height, 10.0f, 200.0f);
				overlay->sliderInt("rows", &c_rows, 1, 100);
				overlay->sliderInt("cols", &c_cols, 1, 100);
			}
			else if (latticeCreateTypeIndex == 1)
			{
				overlay->sliderFloat("radius", &c_width, 10.0f, 100.0f);
				overlay->sliderFloat("height", &c_height, 10.0f, 100.0f);
				overlay->sliderInt("slices", &c_rows, 4, 30);
				overlay->sliderInt("segments", &c_cols, 4, 30);
			}
			else if (latticeCreateTypeIndex == 2)
			{
				overlay->sliderFloat("radius", &c_width, 10.0f, 100.0f);
				overlay->sliderInt("slices", &c_rows, 4, 30);
				overlay->sliderInt("segments", &c_cols, 4, 30);
				ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
			}
			else if (latticeCreateTypeIndex == 3)
			{
				overlay->inputFloat("x", &patchTopLeftX, 1.0f, 1);
				overlay->inputFloat("y", &patchTopLeftY, 1.0f, 1);
				overlay->inputFloat("Width", &patchWidth, 1.0f, 1);
				overlay->inputFloat("Height", &patchHeight, 1.0f, 1);
				if (overlay->button("Add patch")) {
					patchInputs.push_back(OML::Vec2f(patchTopLeftX, patchTopLeftY));
					patchInputs.push_back(OML::Vec2f(patchWidth, patchHeight));
					std::stringstream out; out << std::setprecision(1) << std::fixed;
					out << patchTopLeftX; std::string topLeftXStr = out.str(); out.str(""); out.clear();
					out << patchTopLeftY; std::string topLeftYStr = out.str(); out.str(""); out.clear();
					out << patchWidth; std::string widthStr = out.str(); out.str(""); out.clear();
					out << patchHeight; std::string heightStr = out.str(); out.str(""); out.clear();
					patchInputStrings.push_back("addPatch((" + topLeftXStr + ", " + topLeftYStr + "), " + 
						widthStr + ", " + heightStr + ");");
					patchInputChars.push_back(patchInputStrings.back().c_str());
				}
				ImGui::SameLine();
				if (overlay->button("Delete patch")) {
					if (!patchInputs.empty()) {
						patchInputs.erase(patchInputs.begin() + patchInputIdx);
						patchInputStrings.erase(patchInputStrings.begin() + patchInputIdx);
						patchInputChars.erase(patchInputChars.begin() + patchInputIdx);
					}
				}
			}

			overlay->checkBox("Random colors", &perPatchColors);
			if (!perPatchColors) {
				ImGui::ColorEdit3("Color", &c_col[0]);
			}

			ImGui::EndGroup();
			ImGui::SameLine();
			ImGui::BeginGroup();

			ImGui::Text("Translate");
			ImGui::InputFloat("x##t", &c_pos[0], 0.5f, 1.0f, "%.1f");
			ImGui::InputFloat("y##t", &c_pos[1], 0.5f, 1.0f, "%.1f");
			ImGui::InputFloat("z##t", &c_pos[2], 0.5f, 1.0f, "%.1f");

			ImGui::Text("Rotation");
			ImGui::InputFloat("x##r", &c_rot[0], 1.0f, 1.0f, "%.1f");
			ImGui::InputFloat("y##r", &c_rot[1], 1.0f, 1.0f, "%.1f");
			ImGui::InputFloat("z##r", &c_rot[2], 1.0f, 1.0f, "%.1f");

			ImGui::EndGroup();

			if (overlay->button("Add"))
			{
				Lattice lat;
				lat.setName(std::string(latticeName));
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), c_pos);
				if (c_rot[2] != 0.0f) mat = glm::rotate(mat, c_rot[2], glm::vec3(0, 0, 1));
				if (c_rot[1] != 0.0f) mat = glm::rotate(mat, c_rot[1], glm::vec3(0, 1, 0));
				if (c_rot[0] != 0.0f) mat = glm::rotate(mat, c_rot[0], glm::vec3(1, 0, 0));
				lat.setMatrix(mat);
				lat.setUseRandomPatchColors(perPatchColors);
				lat.setPatchColor(c_col);
				if (latticeCreateTypeIndex == 0) {
					lat.addGrid(OML::Vec2f(-c_width / 2, -c_height / 2), c_width, c_height, c_rows, c_cols);
				}
				else if (latticeCreateTypeIndex == 1) {
					lat.addCylinder(OML::Vec3f(0.0f), c_width, c_height, c_rows, c_cols);
				}
				else if (latticeCreateTypeIndex == 2) {
					lat.addSphere(OML::Vec3f(0.0f), c_width, c_rows, c_cols);
				}
				else if (latticeCreateTypeIndex == 3) {
					for (size_t i = 0; i < patchInputs.size(); i += 2) {
						lat.addPatch(patchInputs[i], patchInputs[i + 1][0], patchInputs[i + 1][1]);
					}
					patchInputs.resize(0); patchInputStrings.resize(0); patchInputChars.resize(0);
				}
				addLattice(lat);
			}

			if (!latticeNames.empty())
			{
				ImGui::Separator();
				overlay->comboBox("", &latticeDeleteIdx, latticeNames);
				ImGui::SameLine();
				if (overlay->button("Delete"))
				{
					removeLattice(latticeNames[latticeDeleteIdx]);
				}
			}

			ImGui::PopStyleVar(1);
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
		Lattice lat("Grid Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addGrid(OML::Vec2f(-m_width / 2, -m_height / 2), m_width, m_height, m_rows, m_cols);
		addLattice(lat);
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
		Lattice lat("Cylinder Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addCylinder(OML::Vec3f(0.0f, 0.0f, 0.0f), m_radius, m_height, m_rows, m_cols);
		addLattice(lat);
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
		Lattice lat("Sphere Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addSphere(OML::Vec3f(0.0f, 0.0f, 0.0f), m_radius, m_segments, m_slices);
		addLattice(lat);
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
		Lattice lat("Non uniform Lattice");
		lat.setUseRandomPatchColors(true);
		OML::Vec3f p00(-30, -20, 0), p10(0, -20, 0), p20(10, -20, 0), p30(20, -20, 0);
		OML::Vec3f p01(-30, 0, 0),	p11(0, 0, 0),	p21(10, 0, 0), p31(20, 0, 0);
		OML::Vec3f p02(-30, 10, 0),	p12(0, 10, 0),	p22(10, 10, 0), p32(20, 10, 0);
		OML::Vec3f p03(-30, 20, 0),	p13(0, 20, 0),	p23(10, 20, 0), p33(20, 20, 0);
		lat.addPatch(p00, p10, p01, p11);
		lat.addPatch(p10, p20, p11, p21);
		lat.addPatch(p20, p30, p21, p31);
		lat.addPatch(p01, p11, p02, p12);
		lat.addPatch(p11, p21, p12, p22);
		lat.addPatch(p21, p31, p22, p32);
		lat.addPatch(p02, p12, p03, p13);
		lat.addPatch(p12, p22, p13, p23);
		lat.addPatch(p22, p32, p23, p33);
		addLattice(lat);
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
		Lattice lat("Non-rectangular Lattice");
		lat.setUseRandomPatchColors(true);
		// Grid where quads have angles != 90
		lat.addPatch(OML::Vec2f(-6,-15), OML::Vec2f(0,-10), OML::Vec2f(-15,-5), OML::Vec2f(0,0));
		lat.addPatch(OML::Vec2f(0,-10), OML::Vec2f(5,-6), OML::Vec2f(0,0), OML::Vec2f(10,0));
		lat.addPatch(OML::Vec2f(5, -6), OML::Vec2f(17,-5), OML::Vec2f(10,0), OML::Vec2f(14,-1));
		lat.addPatch(OML::Vec2f(-15,-5), OML::Vec2f(0,0), OML::Vec2f(-5,13), OML::Vec2f(2,10));
		lat.addPatch(OML::Vec2f(0,0), OML::Vec2f(10,0), OML::Vec2f(2,10), OML::Vec2f(10,10));
		lat.addPatch(OML::Vec2f(10, 0), OML::Vec2f(14, -1), OML::Vec2f(10, 10), OML::Vec2f(20, 10));
		lat.addPatch(OML::Vec2f(-5, 13), OML::Vec2f(2, 10), OML::Vec2f(-8, 20), OML::Vec2f(2, 20));
		lat.addPatch(OML::Vec2f(2, 10), OML::Vec2f(10, 10), OML::Vec2f(2, 20), OML::Vec2f(10, 15));
		lat.addPatch(OML::Vec2f(10, 10), OML::Vec2f(20, 10), OML::Vec2f(10, 15), OML::Vec2f(15, 18));
		addLattice(lat);
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
		Lattice lat("T-locus Lattice");
		lat.setUseRandomPatchColors(true);

		/*lat.addPatch(OML::Vec2f(-20.0f, -20.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(-10.0f, -20.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(0.0f, -20.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(10.0f, -20.0f), 10.0f, 10.0f);*/
		lat.addPatch(OML::Vec2f(-20.0f, -10.0f), 10.0f, 20.0f);

		OML::Vec2f p1(-10.0f, -10.0f);
		OML::Vec2f p2(0.0f, -10.0f);
		OML::Vec2f p3(0.0f, 0.0f);
		lat.addPatch(p1, 10.0f, 20.0f);
		lat.addPatch(p2, 10.0f, 10.0f);
		lat.addPatch(p3, 10.0f, 10.0f);

		lat.addPatch(OML::Vec2f(10.0f, -10.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(10.0f, 0.0f), 10.0f, 10.0f);
		/*lat.addPatch(OML::Vec2f(-20.0f, 10.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(-10.0f, 10.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(0.0f, 10.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(10.0f, 10.0f), 10.0f, 10.0f);*/

		addLattice(lat);
	}
};




class MultiLatticeExample : public LatticeExample
{
public:
	MultiLatticeExample(bool enableValidation) : LatticeExample(enableValidation) {}

	virtual void createLatticeGeometry() override
	{

		Lattice lat("Grid Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addGrid(OML::Vec2f(-25.0f, -25.0f), 20.0f, 20.0f, 3, 3);
		lat.induceLattice();
		lat.initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
		addLattice(lat);

		Lattice lat2("Cylinder Lattice");
		lat2.setUseRandomPatchColors(true);
		lat2.addCylinder(OML::Vec3f(20.0f, 0.0f, -30.0f), 10.0f, 30.0f, 4, 8);
		lat2.induceLattice();
		lat2.initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
		addLattice(lat2);

		Lattice lat3("Sphere Lattice");
		lat3.setUseRandomPatchColors(true);
		lat3.addSphere(OML::Vec3f(0.0f, -20.0f, 20.0f), 8.0f, 6, 6);
		lat3.induceLattice();
		lat3.initVulkanStuff(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
		addLattice(lat3);
	}
};