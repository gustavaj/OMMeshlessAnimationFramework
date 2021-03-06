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

#include "lattice/VulkanLattice.h"

#include "lattice/Shaders.h"

enum LATTICE_CREATE_TYPE { Grid = 0, Cylinder, Sphere, Torus, Patches };
const std::vector<std::string> LATTICE_TYPES = { "Grid", "Cylinder", "Sphere", "Torus", "Patches" };
const std::vector<std::string> EVALUATION_METHODS = { "Direct", "PreSampledImg", "PreSampledImgBatched", "PreSampledBuffer" };
const std::vector<std::string> LOCAL_SURFACE_TYPES = { "Biquadratic Bezier", "Bicubic Bezier", "Plane" };

class LatticeExample : public VulkanExampleBase
{
public:

	std::vector<OML::VulkanLattice> lattices;

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

	int evalMethodIdx = 0;
	int lsTypeIdx = 0;

	float patchTopLeftX = 0.0f, patchTopLeftY = 0.0f;
	float patchWidth = 10.0f, patchHeight = 10.0f;
	int patchInputIdx = 0;
	std::vector<OML::Vec2f> patchInputs;
	std::vector<std::string> patchInputStrings;
	std::vector<const char*> patchInputChars;

	bool ShaderWindowOpen = false;
	int ShaderItemIndex = 0;

	std::vector<ImVec4> ShaderWordColors = {
		ImVec4(0.9, 0.9, 0.9, 1.0), // Normal text
		ImVec4(0.6, 0.6, 0.6, 1.0), // line numbers
		ImVec4(1.0, 0.3, 0.2, 1.0), // Data types
		ImVec4(1.0, 0.6, 0.2, 1.0), // Built-in functions
		ImVec4(0.7, 0.0, 1.0, 1.0), // Macros
		ImVec4(0.0, 0.7, 0.7, 1.0), // Something
		ImVec4(0.0, 0.7, 0.0, 1.0), // Comment
		ImVec4(0.9, 0.8, 0.3, 1.0), // gl_ stuff
		ImVec4(0.6, 0.7, 0.9, 1.0), // User defined Buffers
		ImVec4(0.7, 0.2, 0.3, 1.0), // User defined types
		ImVec4(0.8, 0.4, 0.8, 1.0), // User defined constants
	};
	std::unordered_map<std::string, size_t> GLSLWordMap = {
		// Data types
		{"bool", 2}, {"uint", 2}, {"uvec2", 2}, {"uvec3", 2}, {"uvec4", 2}, {"int", 2}, {"void", 2}, {"vec2", 2}, {"vec3", 2},
		{"vec4", 2}, {"sampler2D", 2}, {"sampler2DArray", 2}, {"mat2", 2}, {"mat3", 2}, {"mat4", 2}, {"float", 2}, {"double", 2},

		// built-in functions
		{"radians", 3}, {"degrees", 3}, {"sin", 3}, {"cos", 3}, {"tan", 3}, {"pow", 3}, {"exp", 3}, {"sqrt", 3}, {"abs", 3},
		{"sign", 3}, {"floor", 3}, {"round", 3}, {"ceil", 3}, {"fract", 3}, {"mod", 3}, {"min", 3}, {"max", 3}, {"clamp", 3},
		{"mix", 3}, {"step", 3}, {"length", 3}, {"distance", 3}, {"dot", 3}, {"cross", 3}, {"normalize", 3}, {"reflect", 3},
		{"transpose", 3}, {"determinant", 3}, {"inverse", 3}, {"texture", 3}, {"textureSize", 3}, {"texture2D", 3}, {"EmitVertex", 3},
		{"EndPrimitive", 3},

		// Macros
		{"#", 4}, {"#define", 4}, {"#if", 4}, {"#ifdef", 4}, {"#ifndef", 4}, {"#else", 4}, {"#elif", 4}, {"#endif", 4}, {"#pragma", 4}, {"#version", 4},

		// Something
		{"attribute", 5}, {"const", 5}, {"uniform", 5}, {"varying", 5}, {"buffer", 5}, {"readonly", 5}, {"writeonly", 5}, {"struct", 5},
		{"layout", 5}, {"flat", 5}, {"smooth", 5}, {"patch", 5}, {"sample", 5}, {"break", 5}, {"continue", 5}, {"do", 5}, {"for", 5},
		{"while", 5}, {"switch", 5}, {"case", 5}, {"default", 5}, {"if", 5}, {"else", 5}, {"in", 5}, {"out", 5}, {"inout", 5}, {"true", 5},
		{"false", 5}, {"invariant", 5}, {"discard", 5}, {"return", 5},

		// Comments
		{"//", 6}, {"/*", 6}, {"*/", 6},

		// gl_ stuff
		{"gl_Position", 7}, {"gl_PointSize", 7}, {"gl_InvocationID", 7}, {"gl_TessLevelInner", 7}, {"gl_TessLevelOuter", 7}, {"gl_TessCoord", 7}, {"gl_FragCoord", 7},

		// User defined Buffers
		{"ubo", 8}, {"lsDataBuffer", 8}, {"matrixBuffer", 8}, {"boundaryBuffer", 8}, {"controlPointBuffer", 8},
		
		// User defined types
		{"LocalSurfaceInfo", 9}, {"Sampler", 9}, {"BoundaryInfo", 9},

		// User defined constants
		{"BFunction_B1Poly", 10}, {"BFunction_B2Poly", 10}, {"BFunction_Lerbs", 10}, 
		{"TessFactorMethod_Static", 10}, {"TessFactorMethod_Dynamic", 10}, {"TessFactorMethod_PixelAccurate", 10},
	};

	LatticeExample(bool enableValidation, OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier, 
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct, bool frontCamera = false)
		: VulkanExampleBase(enableValidation), evalMethodIdx(static_cast<uint32_t>(evalMethod)), 
		lsTypeIdx(static_cast<uint32_t>(lsType))
	{
		title = "Lattice Example";

		// Camera
		zoom = 0.0f;
		zoomSpeed = 100.0f;
		camera.setPerspective(60.0f, (float)width / (float)height, 1.0f, 10000.0f);
		//camera.type = Camera::CameraType::firstperson;
		//camera.movementSpeed = 100.0f;
		//camera.setRotation(glm::vec3(0.0f, 0.0f, 0.0f));
		//camera.setTranslation(glm::vec3(0.0f, 0.0f, -50.0f));
		if (frontCamera)
		{
			// front
			camera.setRotation(glm::vec3(70, 0.0f, 0.0f));
			camera.setTranslation(glm::vec3(0.0f, 0.0f, -250.0f));
		}
		else
		{
			// top
			camera.setRotation(glm::vec3(0.0f, 0.0f, 0.0f));
			camera.setTranslation(glm::vec3(0.0f, 0.0f, -250.0f));
		}

		settings.overlay = true; // ImGui overlay
	}

	~LatticeExample()
	{
		for (auto& lat : lattices) {
			lat.destroyVulkan();
		}
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

	virtual void getEnabledFeatures() override {
		OML::VulkanLattice::CheckAndSetupRequiredPhysicalDeviceFeatures(deviceFeatures, enabledFeatures);
	}

	void loadAssets() {
		createLatticeGeometry();

		windowResized();
	}

	virtual void createLatticeGeometry()
	{
	}

	void buildCommandBuffers() {
		//vkDeviceWaitIdle(device);
		vkWaitForFences(device, MAX_FRAMES_IN_FLIGHT, &inFlightFences[0], VK_TRUE, UINT64_MAX);

		if (!checkCommandBuffers())
		{
			destroyCommandBuffers();
			createCommandBuffers();
		}

		VkCommandBufferBeginInfo cmdBufInfo = vks::initializers::commandBufferBeginInfo();

		VkClearValue clearValues[2];
		//clearValues[0].color = { 1.0f, 1.0f, 1.0f, 1.0f };
		clearValues[0].color = { 0.925f, 0.925f, 0.925f, 1.0f };
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
		// Note: The synchronization code with semaphores and fences that is inside VulkanBaseExample
		// is taken from https://vulkan-tutorial.com/Drawing_a_triangle/Drawing/Rendering_and_presentation#page_Synchronization
		VK_CHECK_RESULT(vkQueueSubmit(queue, 1, &submitInfo, inFlightFences[currentFrame]));

		VulkanExampleBase::submitFrame();
	}

	void prepare() {
		VulkanExampleBase::prepare();
		loadAssets();
		if (!checkCommandBuffers())
			rebuildCmdBuffers = true;
		viewChanged();
		prepared = true;
	}

	virtual void render() override {
		if (!prepared)
			return;
		if (rebuildCmdBuffers)
		{
			buildCommandBuffers();
		}
		for (auto& lat : lattices)
			lat.update(frameTimer);
		draw();
	}

	virtual void viewChanged() override {
		for (auto& lat : lattices)
		{
			lat.onViewChanged(camera.matrices.perspective, camera.matrices.view);
		}
	}

	virtual void windowResized() override {
		for (auto& lat : lattices)
		{
			lat.onWindowResized(width, height);
		}
	}

	void addLattice(OML::VulkanLattice& lattice)
	{
		lattice.setName(lattice.name() + "##" + std::to_string(latticeIdx));
		lattice.setLocalSurfaceType(static_cast<OML::LocalSurfaceType>(lsTypeIdx));
		lattice.setEvaluationMethod(static_cast<OML::EvaluationMethod>(evalMethodIdx));
		lattice.induceLattice();
		lattice.initVulkan(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
		lattice.onViewChanged(camera.matrices.perspective, camera.matrices.view);
		lattices.push_back(std::move(lattice));
		latticeNames.push_back(lattice.name());
		latticeIdx++;
		rebuildCmdBuffers = true;
	}

	void removeLattice(std::string name)
	{
		vkQueueWaitIdle(queue);
		for (auto lat_it = lattices.begin(); lat_it != lattices.end(); lat_it++)
		{
			if (lat_it->name() == name)
			{
				lat_it->destroyVulkan();
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

	virtual void OnUpdateUIOverlay(vks::UIOverlay* overlay) override {
		// Lattice creation
		if (overlay->header("Add/Delete", false))
		{
			if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Patches)
			{
				ImGui::PushItemWidth(260.0f * UIOverlay.scale);
				ImGui::ListBox("", &patchInputIdx, patchInputChars.data(), patchInputChars.size(), 4);
				ImGui::PushItemWidth(110.0f * UIOverlay.scale);
			}

			ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(5.0f, 1.0f));

			ImGui::BeginGroup();

			overlay->comboBox("lsType", &lsTypeIdx, LOCAL_SURFACE_TYPES);

			overlay->comboBox("Type", &latticeCreateTypeIndex, LATTICE_TYPES);
			ImGui::InputText("Name", latticeName, IM_ARRAYSIZE(latticeName));

			if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Grid)
			{
				overlay->sliderFloat("width", &c_width, 10.0f, 200.0f);
				overlay->sliderFloat("height", &c_height, 10.0f, 200.0f);
				overlay->sliderInt("rows", &c_rows, 1, 100);
				overlay->sliderInt("cols", &c_cols, 1, 100);
			}
			else if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Cylinder)
			{
				overlay->sliderFloat("radius", &c_width, 10.0f, 100.0f);
				overlay->sliderFloat("height", &c_height, 10.0f, 100.0f);
				overlay->sliderInt("slices", &c_rows, 4, 30);
				overlay->sliderInt("segments", &c_cols, 4, 30);
			}
			else if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Sphere)
			{
				overlay->sliderFloat("radius", &c_width, 10.0f, 100.0f);
				overlay->sliderInt("slices", &c_rows, 4, 30);
				overlay->sliderInt("segments", &c_cols, 4, 30);
				ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
			}
			else if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Torus)
			{
				overlay->inputFloat("radius", &c_width, 5.0f, 0);
				overlay->inputFloat("wheelRadius", &c_height, 5.0f, 0);
				ImGui::InputInt("segments", &c_rows);
				ImGui::InputInt("slices", &c_cols);
			}
			else if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Patches)
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
			overlay->comboBox("evalMethod", &evalMethodIdx, EVALUATION_METHODS);

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
				OML::VulkanLattice lat;
				lat.setName(std::string(latticeName));
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), c_pos);
				if (c_rot[2] != 0.0f) mat = glm::rotate(mat, c_rot[2], glm::vec3(0, 0, 1));
				if (c_rot[1] != 0.0f) mat = glm::rotate(mat, c_rot[1], glm::vec3(0, 1, 0));
				if (c_rot[0] != 0.0f) mat = glm::rotate(mat, c_rot[0], glm::vec3(1, 0, 0));
				lat.setMatrix(mat);
				lat.setUseRandomPatchColors(perPatchColors);
				lat.setPatchColor(c_col);
				if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Grid) {
					lat.addGrid(OML::Vec2f(-c_width / 2, -c_height / 2), c_width, c_height, c_rows, c_cols);
				}
				else if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Cylinder) {
					lat.addCylinder(OML::Vec3f(0.0f), c_width, c_height, c_rows, c_cols);
				}
				else if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Sphere) {
					lat.addSphere(OML::Vec3f(0.0f), c_width, c_rows, c_cols);
				}
				else if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Torus)
				{
					lat.addTorus(OML::Vec3f(0.0f), c_width, c_height, c_rows, c_cols);
				}
				else if (latticeCreateTypeIndex == LATTICE_CREATE_TYPE::Patches) {
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
				rebuildCmdBuffers = true;
			}
		}

		ImGui::Separator();

		if (ImGui::Button("Shaders"))
		{
			ShaderWindowOpen = !ShaderWindowOpen;
		}



		if(ShaderWindowOpen) {
			ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.1, 0.1, 0.1, 1.0));

			ImGui::SetNextWindowSize(ImVec2(800, 800));

			ImGui::Begin("Shaders Viewer");

			overlay->comboBox("Shaders##combo", &ShaderItemIndex, OML::Shaders::ShaderNames);
			std::string code = OML::Shaders::ShaderSources[OML::Shaders::ShaderNames[ShaderItemIndex]];

			int lineNr = 0;
			std::string lineNrStr;
			std::string line;
			std::string word;

			bool isComment = false;
			size_t commentColIdx = 6;

			std::string nastyChars = "(){}[].;";

			std::stringstream liness(code);

			OML::Random rnd;

			while (std::getline(liness, line, '\n'))
			{
				lineNrStr = (lineNr < 10 ? "00" : (lineNr < 100 ? "0" : "")) + std::to_string(++lineNr) + ": ";

				ImGui::TextColored(ShaderWordColors[1], lineNrStr.c_str());

				std::stringstream wordss(line);
				while (std::getline(wordss, word, ' '))
				{
					if (isComment || word.find("//") != std::string::npos) {
						isComment = true;
						ImGui::SameLine();
						ImGui::TextColored(ShaderWordColors[commentColIdx], word.c_str());
						continue;
					}

					int colIdx = 0;
					size_t pos = 0;
					size_t idx = word.find_first_of(nastyChars, pos);
					if (idx == std::string::npos || (idx == pos && word.size() == 1)) {
						auto res = GLSLWordMap.find(word);
						if (res != GLSLWordMap.end()) colIdx = res->second;
						ImGui::SameLine();
						ImGui::TextColored(ShaderWordColors[colIdx], word.c_str());
					}
					else {
						bool first = true;
						while (idx != std::string::npos)
						{
							std::string subword = word.substr(pos, idx - pos);
							pos = idx + 1;
							auto res = GLSLWordMap.find(subword);
							if (res != GLSLWordMap.end()) colIdx = res->second;
							if (first) {
								ImGui::SameLine();
								first = false;
							}
							else {
								ImGui::SameLine(0.0f, 0.0f);
							}
							ImGui::TextColored(ShaderWordColors[colIdx], subword.c_str());
							colIdx = 0;
							ImGui::SameLine(0.0f, 0.0f);
							//ImGui::SetCursorPosX(ImGui::GetCursorPosX() - 5.0f);
							std::string specialChar = word.substr(idx, 1);
							ImGui::TextColored(ShaderWordColors[colIdx], specialChar.c_str());
							idx = word.find_first_of(nastyChars, pos);
							if (idx == word.size() - 1) break;
						}
						if (pos != word.size())
						{
							std::string lastsubword = word.substr(pos);
							auto res = GLSLWordMap.find(lastsubword);
							if (res != GLSLWordMap.end()) colIdx = res->second;
							ImGui::SameLine(0.0f, 0.0f);
							ImGui::TextColored(ShaderWordColors[colIdx], lastsubword.c_str());
						}
					}
				}

				isComment = false;
			}

			if (ImGui::Button("Close##Shaders")) ShaderWindowOpen = false;

			ImGui::End();

			ImGui::PopStyleColor();
		}
	}

};

class GridLatticeExample : public LatticeExample
{
public:
	int m_rows, m_cols;
	float m_width, m_height;

	GridLatticeExample(bool enableValidation, float width, float height, int rows, int cols,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod, false), 
		m_width(width), m_height(height), m_rows(rows), m_cols(cols)
	{
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("Grid Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addGrid(OML::Vec2f(-m_width / 2, -m_height / 2), m_width, m_height, m_rows, m_cols);
		addLattice(lat);
	}
};

class RandomGridLatticeExample : public LatticeExample
{
public:
	int m_rows, m_cols;
	float m_width, m_height;

	RandomGridLatticeExample(bool enableValidation, float width, float height, int rows, int cols,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod), 
		m_width(width), m_height(height), m_rows(rows), m_cols(cols)
	{
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("Random Grid Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addGridRandom(OML::Vec2f(-m_width / 2, -m_height / 2), m_width, m_height, m_rows, m_cols);
		addLattice(lat);
	}
};

class CylinderLatticeExample : public LatticeExample
{
public:
	int m_rows, m_cols;
	float m_radius, m_height;

	CylinderLatticeExample(bool enableValidation, float radius, float height, int rows, int cols,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod), 
		m_radius(radius), m_height(height), m_rows(rows), m_cols(cols)
	{
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("Cylinder Lattice");
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

	SphereLatticeExample(bool enableValidation, float radius, int segments, int slices,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod), 
		m_radius(radius), m_segments(segments), m_slices(slices)
	{
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("Sphere Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addSphere(OML::Vec3f(0.0f, 0.0f, 0.0f), m_radius, m_segments, m_slices);
		addLattice(lat);
	}
};

class TorusLatticeExample : public LatticeExample
{
public:
	int m_segments, m_slices;
	float m_radius, m_wheelRadius;

	TorusLatticeExample(bool enableValidation, float radius, float wheelRadius, int segments, int slices,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod),
		m_radius(radius), m_wheelRadius(wheelRadius), m_segments(segments), m_slices(slices)
	{
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("Torus Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addTorus(OML::Vec3f(0.0f, 0.0f, 0.0f), m_radius, m_wheelRadius, m_segments, m_slices);
		addLattice(lat);
	}
};

class NonUniformGridExample : public LatticeExample
{
public:
	NonUniformGridExample(bool enableValidation,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod)
	{
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("Non uniform Lattice");
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
	NonRectangularExample(bool enableValidation,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod)
	{
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("Non-rectangular Lattice");
		lat.setUseRandomPatchColors(true);
		// Grid where quads have angles != 90
		lat.addPatch(OML::Vec2f(-60,-150), OML::Vec2f(0,-100), OML::Vec2f(-150,-50), OML::Vec2f(0,0));
		lat.addPatch(OML::Vec2f(0,-100), OML::Vec2f(50,-60), OML::Vec2f(0,0), OML::Vec2f(100,0));
		lat.addPatch(OML::Vec2f(50, -60), OML::Vec2f(170,-50), OML::Vec2f(100,0), OML::Vec2f(140,-10));
		lat.addPatch(OML::Vec2f(-150,-50), OML::Vec2f(0,0), OML::Vec2f(-50,130), OML::Vec2f(20,100));
		lat.addPatch(OML::Vec2f(0,0), OML::Vec2f(100,0), OML::Vec2f(20,100), OML::Vec2f(100,100));
		lat.addPatch(OML::Vec2f(100, 0), OML::Vec2f(140, -10), OML::Vec2f(100, 100), OML::Vec2f(200, 100));
		lat.addPatch(OML::Vec2f(-50, 130), OML::Vec2f(20, 100), OML::Vec2f(-80, 200), OML::Vec2f(20, 200));
		lat.addPatch(OML::Vec2f(20, 100), OML::Vec2f(100, 100), OML::Vec2f(20, 200), OML::Vec2f(100, 150));
		lat.addPatch(OML::Vec2f(100, 100), OML::Vec2f(200, 100), OML::Vec2f(100, 150), OML::Vec2f(150, 180));
		addLattice(lat);
	}
};

class TLocusExample : public LatticeExample
{
public:
	TLocusExample(bool enableValidation,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod)
	{
	}

	virtual void createLatticeGeometry() override
	{
		//OML::VulkanLattice bkg("background");
		//bkg.setUseRandomPatchColors(false);
		//bkg.addPatch(OML::Vec2f(-20, -10), 30, 20);
		//addLattice(bkg);

		OML::VulkanLattice lat("T-locus Lattice");
		lat.setUseRandomPatchColors(false);

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

		//lat.addPatch(OML::Vec2f(10.0f, -10.0f), 10.0f, 10.0f);
		//lat.addPatch(OML::Vec2f(10.0f, 0.0f), 10.0f, 10.0f);
		/*lat.addPatch(OML::Vec2f(-20.0f, 10.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(-10.0f, 10.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(0.0f, 10.0f), 10.0f, 10.0f);
		lat.addPatch(OML::Vec2f(10.0f, 10.0f), 10.0f, 10.0f);*/

		addLattice(lat);
	}
};

class TLocusx4Example : public LatticeExample
{
public:
	TLocusx4Example(bool enableValidation,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod)
	{
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("T-locus x4 Lattice");
		lat.setUseRandomPatchColors(true);

		lat.addPatch(OML::Vec2f(-30, -30), 10, 10);
		lat.addPatch(OML::Vec2f(-20, -30), 10, 10);
		lat.addPatch(OML::Vec2f(-10, -30), 10, 10);
		lat.addPatch(OML::Vec2f(0, -30), 10, 10);
		lat.addPatch(OML::Vec2f(10, -30), 10, 10);
		lat.addPatch(OML::Vec2f(20, -30), 10, 10);

		lat.addPatch(OML::Vec2f(-30, -20), 10, 10);
		lat.addPatch(OML::Vec2f(-20, -20), 10, 10);
		lat.addPatch(OML::Vec2f(-10, -20), 20, 10);
		lat.addPatch(OML::Vec2f(10, -20), 10, 10);
		lat.addPatch(OML::Vec2f(20, -20), 10, 10);

		lat.addPatch(OML::Vec2f(-30, -10), 10, 10);
		lat.addPatch(OML::Vec2f(-20, -10), 10, 20);
		lat.addPatch(OML::Vec2f(-10, -10), 20, 20);
		lat.addPatch(OML::Vec2f(10, -10), 10, 20);
		lat.addPatch(OML::Vec2f(20, -10), 10, 10);

		lat.addPatch(OML::Vec2f(-30, 0), 10, 10);
		lat.addPatch(OML::Vec2f(20, 0), 10, 10);

		lat.addPatch(OML::Vec2f(-30, 10), 10, 10);
		lat.addPatch(OML::Vec2f(-20, 10), 10, 10);
		lat.addPatch(OML::Vec2f(-10, 10), 20, 10);
		lat.addPatch(OML::Vec2f(10, 10), 10, 10);
		lat.addPatch(OML::Vec2f(20, 10), 10, 10);

		lat.addPatch(OML::Vec2f(-30, 20), 10, 10);
		lat.addPatch(OML::Vec2f(-20, 20), 10, 10);
		lat.addPatch(OML::Vec2f(-10, 20), 10, 10);
		lat.addPatch(OML::Vec2f(0, 20), 10, 10);
		lat.addPatch(OML::Vec2f(10, 20), 10, 10);
		lat.addPatch(OML::Vec2f(20, 20), 10, 10);

		// 29 patches, 44 loci, 36 local

		addLattice(lat);
	}
};




class MultiLatticeExample : public LatticeExample
{
public:
	MultiLatticeExample(bool enableValidation,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct) 
		: LatticeExample(enableValidation, lsType, evalMethod) {}

	virtual void createLatticeGeometry() override
	{

		OML::VulkanLattice lat("Grid Lattice");
		lat.setUseRandomPatchColors(true);
		lat.addGrid(OML::Vec2f(-25.0f, -25.0f), 20.0f, 20.0f, 3, 3);
		addLattice(lat);

		OML::VulkanLattice lat2("Cylinder Lattice");
		lat2.setUseRandomPatchColors(true);
		lat2.addCylinder(OML::Vec3f(20.0f, 0.0f, -30.0f), 10.0f, 30.0f, 4, 8);
		addLattice(lat2);

		OML::VulkanLattice lat3("Sphere Lattice");
		lat3.setUseRandomPatchColors(true);
		lat3.addSphere(OML::Vec3f(0.0f, -20.0f, 20.0f), 8.0f, 6, 6);
		addLattice(lat3);
	}
};



class AllMethodsGridExample : public LatticeExample
{
public:
	AllMethodsGridExample(bool enableValidation, bool moveLocals,
		OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Quadratic_Bezier,
		OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct)
		: LatticeExample(enableValidation, lsType, evalMethod, false), m_moveLocals(moveLocals) {}

	bool m_moveLocals = false;

	virtual void createLatticeGeometry() override
	{
		std::vector<std::string> names = { "Direct", "Pre_Img", "Pre_Batch", "Pre_Buf", "Pre_Buf_NoInterp" };
		std::vector<OML::EvaluationMethod> evals = { OML::EvaluationMethod::Direct, OML::EvaluationMethod::Pre_Sampled_Image,
			OML::EvaluationMethod::Pre_Sampled_Image_Batched, OML::EvaluationMethod::Pre_Sampled_Buffer,
			OML::EvaluationMethod::Pre_Sampled_Buffer_No_Interpolation };
		std::vector<glm::vec3> colors = {glm::vec3(0,0,0), glm::vec3(0,1,0), glm::vec3(0,0,1), glm::vec3(1,1,0), glm::vec3(1,0,1) };
		for (size_t i = 0; i < 5; i++)
		{
			OML::VulkanLattice lat(names[i]);
			lat.setUseRandomPatchColors(true);
			lat.setPatchColor(colors[i]);
			lat.setLocalSurfaceType(static_cast<OML::LocalSurfaceType>(lsTypeIdx));
			lat.setDraw(false);
			lat.setEvaluationMethod(evals[i]);
			lat.addGrid(OML::Vec2f(-100.0f, -100.0f), 200.0f, 200.0f, 3, 3);
			lat.induceLattice();
			lat.initVulkan(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
			lat.onViewChanged(camera.matrices.perspective, camera.matrices.view);

			if (m_moveLocals)
			{
				lat.translateLocalSurface(0, glm::vec3(0, 0, 50));
				lat.translateLocalSurface(1, glm::vec3(0, 0, -50));
				lat.translateLocalSurface(2, glm::vec3(0, 0, 25));
				lat.translateLocalSurface(3, glm::vec3(0, 0, 0));
				lat.translateLocalSurface(4, glm::vec3(0, 0, -25));
				lat.translateLocalSurface(5, glm::vec3(0, 0, 10));
				lat.translateLocalSurface(6, glm::vec3(0, 0, -12));
				lat.translateLocalSurface(7, glm::vec3(0, 0, 10));
				lat.translateLocalSurface(8, glm::vec3(0, 0, 10));
				lat.translateLocalSurface(9, glm::vec3(0, 0, 25));
				lat.translateLocalSurface(10, glm::vec3(0, 0, -12));
				lat.translateLocalSurface(11, glm::vec3(0, 0, 0));
				lat.translateLocalSurface(12, glm::vec3(0, 0, 30));
				lat.translateLocalSurface(13, glm::vec3(0, 0, -12));
				lat.translateLocalSurface(14, glm::vec3(0, 0, 2));
				lat.translateLocalSurface(15, glm::vec3(0, 0, 21));
				lat.updateMatrixBuffer();
			}


			lattices.push_back(std::move(lat));
			latticeNames.push_back(names[i]);
			latticeIdx++;
		}
		rebuildCmdBuffers = true;
	}
};

class PixelAccTest : public LatticeExample
{
public:
	PixelAccTest(bool enableValidation)
		: LatticeExample(enableValidation, OML::LocalSurfaceType::Plane, OML::EvaluationMethod::Direct){
		camera.setPerspective(60.0f, (float)width / (float)height, 1.0f, 10000.0f);
		//camera.setTranslation(glm::vec3(0, 100.0f, -750.0f));
	}

	virtual void createLatticeGeometry() override
	{
		camera.setRotation(glm::vec3(70, 0.0f, 0.0f));
		camera.setTranslation(glm::vec3(0.0f, 0.0f, -130.0f));
		zoomSpeed = 25.0f;

		OML::VulkanLattice lat("PixelAccTestLattice");
		lat.setUseRandomPatchColors(true);
		lat.setLocalSurfaceType(OML::LocalSurfaceType::Cubic_Bezier);
		lat.setEvaluationMethod(OML::EvaluationMethod::Direct);
		int rows = 20;
		int cols = 20;
		lat.addGrid(OML::Vec2f(-750.0f, -750.0f), 1000.0f, 1000.0f, rows, cols);
		lat.induceLattice();
		lat.initVulkan(&device, vulkanDevice, &queue, &cmdPool, &descriptorPool, &renderPass, nullptr);
		lat.onViewChanged(camera.matrices.perspective, camera.matrices.view);

		for (size_t j = 0; j < (rows + 1); j++)
		{
			for (size_t i = 0; i < (cols + 1); i++)
			{
				float ti = ((float)i / (float)(rows)) * 2 * M_PI;
				float tj = ((float)j / (float)(cols)) * 2 * M_PI;
				lat.translateLocalSurface(j * cols + i, glm::vec3(
					0, 0, std::sin(4*ti) * std::sin(5*tj) * 30.0f));
			}
		}
		/*for (size_t i = 0; i < 2601; i++)
		{
			float t = ((float)i / 2601.0f) * 2 * M_PI;
			lat.translateLocalSurface(i, glm::vec3(0, 0, std::sin(t) * 1000.0f));
		}*/
		lat.updateMatrixBuffer();

		lattices.push_back(std::move(lat));
		latticeNames.push_back("PixelAccTestLattice");
		latticeIdx++;
		rebuildCmdBuffers = true;
	}
};



/*
	Benchmark setup:
		- fullscreen
		- no validation layers
		- 64-bit release build
		- 2 sec warmup
		- 5 sec duration
		- B2Poly B-function
		- no gui
		- no pipeline queries/timings

*/
class BenchmarkGrid : public LatticeExample
{
public:
	int tess, rowcol;
	std::string resultDir;
	double warmup = 2.0;
	double duration = 5.0;

	BenchmarkGrid(int rowcol, int tess, std::string resultDir, OML::LocalSurfaceType lsType, OML::EvaluationMethod evalMethod)
		: LatticeExample(false, lsType, evalMethod), rowcol(rowcol), tess(tess), resultDir(resultDir)
	{
		settings.overlay = false;
		settings.fullscreen = true;
	}

	virtual void createLatticeGeometry() override
	{
		OML::VulkanLattice lat("Grid");
		lat.addGrid(OML::Vec2f(-250.0f, -250.0f), 500.0f, 500.0f, rowcol, rowcol);
		lat.setUseRandomPatchColors(true);
		lat.setTessellationFactors(tess, tess);

		std::string benchmarkFilename = resultDir;
		benchmarkFilename += std::to_string(rowcol) + "x" + std::to_string(rowcol) + "grid_";
		benchmarkFilename += std::to_string(tess) + "x" + std::to_string(tess) + "tess_";
		benchmarkFilename += "lsType" + std::to_string(lsTypeIdx) + "_evalMethod"
			+ std::to_string(evalMethodIdx) + "_" + std::to_string(static_cast<int>(warmup)) 
			+ "w_" + std::to_string(static_cast<int>(duration)) + "d.txt";

		benchmark.active = true;
		benchmark.warmup = warmup;
		benchmark.duration = duration;
		benchmark.filename = benchmarkFilename;

		//camera.setPerspective(60.0f, (float)width / (float)height, 1.0f, 10000.0f);
		camera.setRotation(glm::vec3(0, 0, 0));
		camera.setTranslation(glm::vec3(0, 0, -500));

		addLattice(lat);
	}
};