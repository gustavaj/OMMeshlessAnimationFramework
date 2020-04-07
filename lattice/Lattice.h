#pragma once

#include <vulkan/vulkan.h>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <unordered_map>

// Sascha willems
#include "../vulkan/VulkanDevice.hpp"
#include "../vulkan/VulkanBuffer.hpp"
#include "../vulkan/VulkanUIOverlay.h"

#include "LocalSurface.h"


namespace OML {

	using Vec2f = OpenMesh::Vec2f;
	using Vec3f = OpenMesh::Vec3f;
	using Col3 = OpenMesh::Vec3uc;

	// TODO: Move maybe?
	struct BoundaryInfo
	{
		BoundaryInfo()
			: us(0.0f), ue(1.0f), vs(0.0f), ve(1.0f) {}
		BoundaryInfo(float us, float ue, float vs, float ve)
			: us(us), ue(ue), vs(vs), ve(ve) {}

		float us, ue, vs, ve;
	};

	struct Locus
	{
		LocalSurface* localSurface;
		OpenMesh::VertexHandle vh;
		std::map<OpenMesh::FaceHandle, BoundaryInfo> faceMappings;
	};

	struct Patch
	{
		std::array<size_t, 4> lociIndices;
		OpenMesh::FaceHandle fh;
	};

	// End TODO




	// Custom traits
	struct LatticeTraits : public OpenMesh::DefaultTraits
	{
		typedef Vec3f Point;

		VertexAttributes(OpenMesh::Attributes::Color | OpenMesh::Attributes::Status); // Use color for vailence

		EdgeAttributes(OpenMesh::Attributes::Color | OpenMesh::Attributes::Status); // Use color for on/inside boundary

		HalfedgeAttributes(OpenMesh::Attributes::Status); // Use for deleting half edges

		FaceAttributes(OpenMesh::Attributes::Status); // Use for deleting faces
	};

	namespace LatticeProperties {
		static OpenMesh::VPropHandleT<size_t> LocusValence;
		static OpenMesh::VPropHandleT<uint32_t> LocusLocalSurfaceIdx;
		static OpenMesh::VPropHandleT<std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo>> FaceMappings;
	};

	const Col3 BOUNDARY_EDGE_COLOR = Col3(0, 0, 0);
	const Col3 INNER_EDGE_COLOR = Col3(200, 200, 200);
	const std::vector<Col3> LOCUS_VALENCE_COLOR = {
		Col3(0,255,255), Col3(255,255,255), Col3(255,0,0), Col3(0,255,0), Col3(0,0,255), Col3(0,255,255)
	};

	const std::vector<std::string> BFunctionNames = {
		"3x^2-2x^3", "6x^5-15x^4+10x^3", "1/(1+e^(1/x-1/(1-x)))"
	};

	// Struct containing info for a vertex used for points and lines in the lattice grid.
	struct GridVertex
	{
		Vec3f pos;
		Col3 col;

		static VkVertexInputBindingDescription GetBindingDescription() {
			VkVertexInputBindingDescription bindingDescription = {};
			bindingDescription.binding = 0;
			bindingDescription.stride = sizeof(GridVertex);
			bindingDescription.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;
			return bindingDescription;
		}

		static std::vector<VkVertexInputAttributeDescription> GetAttributeDesctiptions() {
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
	};

	class Lattice : public OpenMesh::PolyMesh_ArrayKernelT<LatticeTraits>
	{
	public:
		Lattice();
		Lattice(std::string name);
		~Lattice();

		// Framework
		// Functions for creating patches
		void addPatch(Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight);
		void addPatch(Vec2f topLeft, Vec2f topRight, Vec2f bottomLeft, Vec2f bottomRight);
		void addPatch(Vec2f topLeft, float width, float height);
		void addPatch(Vec3f topLeft, float width, float height, float depth);
		void addGrid(Vec2f topLeft, float width, float height, int rows, int cols);
		void addDonut(Vec2f topLeft, float outerRadius, float innerRadius);
		void addCylinder(Vec3f center, float radius, float height, int rows, int cols);
		void addSphere(Vec3f center, float radius, int segments, int slices);

		void induceLattice();

		void initVulkanStuff(
			VkDevice* device, vks::VulkanDevice* vulkanDevice, 
			VkQueue* queue, VkCommandPool* commandPool,
			VkDescriptorPool* descriptorPool, VkRenderPass* renderPass,
			VkAllocationCallbacks* allocator);
		void destroyVulkanStuff();
		void addToCommandbuffer(VkCommandBuffer& commandBuffer);
		void onViewChanged(glm::mat4 projection, glm::mat4 view);

		// ImGui
		bool onUpdateUIOverlay(vks::UIOverlay* overlay);

		// Setters
		void setName(std::string name) { m_name = name; }
		void setDraw(bool draw) { m_draw = draw; }
		void setDrawLatticeGrid(bool drawLatticeGrid) { m_drawLatticeGrid = drawLatticeGrid; }
		void setDrawLocalSurfaces(bool drawLocalSurfaces) { m_drawLocalSurfaces = drawLocalSurfaces; }
		void setDrawSurface(bool drawSurface) { m_drawSurface = drawSurface; }
		bool setDrawNormals(bool drawNormals) { m_drawNormals = drawNormals; }
		void setWireframe(bool wireframe) { m_wireframe = wireframe; }
		void setDrawPixelAccurate(bool drawPixelAccurate) { m_drawPixelAccurate = drawPixelAccurate; }

		int numPatches() { return n_faces(); }

	private:
		void setupEdgeColor();
		void setupLociValenceAndPointColor();
		void handleTLoci();
		void setupLocalSurfacesAndPatches();

		void addLatticeProperties();

		void addLocalSurfaceOnLocus(
			OpenMesh::VertexHandle locus, OpenMesh::VertexHandle next_locus,
			OpenMesh::VertexHandle prev_locus, int locusIndexOnFace);
		void addLocalSurfaceOnCornerLocus(
			OpenMesh::VertexHandle locus, OpenMesh::VertexHandle next_locus,
			OpenMesh::VertexHandle prev_locus, int locusIndexOnFace);
		void addLocalSurfaceOnEdgeLocus(
			OpenMesh::VertexHandle locus, OpenMesh::VertexHandle next_locus,
			OpenMesh::VertexHandle prev_locus, int locusIndexOnFace);
		void addLocalSurfaceOnInnerLocus(
			OpenMesh::VertexHandle locus, OpenMesh::VertexHandle next_locus,
			OpenMesh::VertexHandle prev_locus, int locusIndexOnFace);

		void addLocalSurface(
			Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight, Vec3f offset);
		void addLocalSurface(
			Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
			Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
			Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight, Vec3f offset);

		inline Vec3f lerp(Vec3f a, Vec3f b, float t) {
			if (t < 1e-5) return a;
			else if ((1.0f - t) < 1e-5) return b;
			else return a + t * (b - a);
		}

		// TODO: Change this into its own class?
		struct {
			VkBuffer buffer = VK_NULL_HANDLE;
			VkDeviceMemory memory = VK_NULL_HANDLE;
			int count;
		} m_pointsBuffer, m_linesBuffer, m_localSurfaceBuffer;

		// Vulkan functions
		void createDeviceLocalBuffer(
			VkBuffer& buffer, VkDeviceMemory& memory, void* data, 
			uint32_t bufferSize, VkBufferUsageFlagBits usage);
		void createBuffers();
		void prepareUniformBuffers();
		void setupDescriptorSetLayouts();
		VkPipelineShaderStageCreateInfo loadShader(std::string fileName, VkShaderStageFlagBits stage);
		void preparePipelines();
		void setupDescriptorPool();
		void setupDescriptorSets();
		void updateUniformBuffer();


		struct {
			int tessInner = 5;
			int tessOuter = 5;
			int bFunctionIndex = 0;
			alignas(16) glm::mat4 projection = glm::mat4(1.0f);
			alignas(16) glm::mat4 modelview = glm::mat4(1.0f);
		} m_ubo;
		vks::Buffer m_uniformBuffer;

		VkPipeline m_pointsPipeline;
		VkPipeline m_linesPipeline;
		VkPipeline m_localSurfacePipeline;
		VkPipeline m_localSurfaceWireframePipeline;

		VkPipelineLayout m_pipelineLayout;
		VkDescriptorSetLayout m_descriptorSetLayout;
		VkDescriptorSet m_descriptorSet;


		std::string m_name;
		bool m_draw = true;
		bool m_drawLatticeGrid = true;
		bool m_drawLocalSurfaces = false;
		bool m_drawSurface = true;
		bool m_drawNormals = false;
		bool m_wireframe = false;
		bool m_drawPixelAccurate = false;


		std::vector<LocalSurface> m_localSurfaces;
		std::vector<Locus> m_loci;
		std::vector<Patch> m_patches;



		// Vulkan stuff
		VkDevice* m_device;
		vks::VulkanDevice* m_vulkanDevice;
		VkCommandPool* m_commandPool;
		VkDescriptorPool* m_descriptorPool;
		VkRenderPass* m_renderPass;
		VkQueue* m_queue;
		VkAllocationCallbacks* m_allocator;
	};

}