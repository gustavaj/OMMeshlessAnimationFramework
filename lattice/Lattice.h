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

//#include "LocalSurface.h"


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

	enum class LocalSurfaceType
	{
		Bezier3x3 = 0
	};

	struct Locus
	{
		Locus() : controlPointIndex(0), controlPointCount(0), matrixIndex(0) {}

		uint32_t controlPointIndex;
		uint32_t controlPointCount;
		uint32_t matrixIndex;
		OpenMesh::VertexHandle vh;
		std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo> faceMappings;
	};

	struct Patch
	{
		std::array<Locus*, 4> loci;
		OpenMesh::FaceHandle fh;
	};

	// Holds the information for one local surface that is sent to the vertex shader
	struct LocalSurfaceVertex
	{
		LocalSurfaceVertex() : LocalSurfaceVertex(0, 0, 0, 0) {}
		LocalSurfaceVertex(
			uint32_t controlPointIndex, uint32_t controlPointCount,
			uint32_t matrixIndex, uint32_t boundaryIndex)
			: controlPointIndex(controlPointIndex), controlPointCount(controlPointCount),
			matrixIndex(matrixIndex), boundaryIndex(boundaryIndex) {}
		uint32_t controlPointIndex, controlPointCount, matrixIndex, boundaryIndex;

		static std::vector<VkVertexInputBindingDescription> GetBindingDescriptions()
		{
			std::vector<VkVertexInputBindingDescription> bindings;
			bindings.resize(1);
			bindings[0].binding = 0;
			bindings[0].stride = sizeof(LocalSurfaceVertex);
			bindings[0].inputRate = VK_VERTEX_INPUT_RATE_VERTEX;
			return bindings;
		}

		static std::vector<VkVertexInputAttributeDescription> GetAttributeDescriptions()
		{
			std::vector<VkVertexInputAttributeDescription> attributeDescriptions = {};
			attributeDescriptions.resize(1);
			attributeDescriptions[0].binding = 0;
			attributeDescriptions[0].location = 0;
			attributeDescriptions[0].format = VK_FORMAT_R32G32B32A32_UINT;
			attributeDescriptions[0].offset = offsetof(LocalSurfaceVertex, controlPointIndex);

			return attributeDescriptions;
		}
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
		static OpenMesh::VPropHandleT<size_t> VertexValence;
		static OpenMesh::VPropHandleT<uint32_t> LocusIndex;
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
		// Framework
		Lattice();
		Lattice(std::string name);
		~Lattice();

		// Functions for creating patches
		/* Add patch in space */
		void addPatch(Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight);
		/* Add patch in xy-plane */
		void addPatch(Vec2f topLeft, Vec2f topRight, Vec2f bottomLeft, Vec2f bottomRight);
		/* Add patch in xy-plane */
		void addPatch(Vec2f topLeft, float width, float height);
		/* Add patch in space */
		void addPatch(Vec3f topLeft, float width, float height, float depth);
		/* Add grid in xy-plane */
		void addGrid(Vec2f topLeft, float width, float height, int rows, int cols);
		/* Add donut in xy-plane */
		void addDonut(Vec2f topLeft, float outerRadius, float innerRadius);
		/* Add cylinder in space */
		void addCylinder(Vec3f center, float radius, float height, int rows, int cols);
		/* Add sphere in space */
		void addSphere(Vec3f center, float radius, int segments, int slices);

		/* Finalize lattice creation */
		void induceLattice();

		/* Create buffers and pipelines and stuff for Vulkan */
		void initVulkanStuff(
			VkDevice* device, vks::VulkanDevice* vulkanDevice, 
			VkQueue* queue, VkCommandPool* commandPool,
			VkDescriptorPool* descriptorPool, VkRenderPass* renderPass,
			VkAllocationCallbacks* allocator);
		/* Clean up vulkan stuff */
		void destroyVulkanStuff();
		/* Add drawing commands to a command buffer */
		void addToCommandbuffer(VkCommandBuffer& commandBuffer);
		/* Update projection and view matrices */
		void onViewChanged(glm::mat4 projection, glm::mat4 view);

		/* Update */
		void update(double dt);

		/* Update UIOverlay with information about the lattice */
		bool onUpdateUIOverlay(vks::UIOverlay* overlay);

		/* Set the name of the lattice */
		void setName(std::string name) { m_name = name; }
		/* Sets wheter the lattice is displayed */
		void setDraw(bool draw) { m_draw = draw; }
		/* Sets wheter the lattice should be animated */
		void setAnimate(bool animate) { m_animate = animate; }
		/* Sets wheter the grid is displayed */
		void setDrawLatticeGrid(bool drawLatticeGrid) { m_drawLatticeGrid = drawLatticeGrid; }
		/* Sets wheter local surfaces are displayed */
		void setDrawLocalSurfaces(bool drawLocalSurfaces) { m_drawLocalSurfaces = drawLocalSurfaces; }
		/* Sets wheter the Blending surface is displayed */
		void setDrawSurface(bool drawSurface) { m_drawSurface = drawSurface; }
		/* Sets wheter normals are displayed */
		bool setDrawNormals(bool drawNormals) { m_drawNormals = drawNormals; }
		/* Sets wireframe mode */
		void setWireframe(bool wireframe) { m_wireframe = wireframe; }
		/* Sets wheter the surface is drawn pixel-accurate */
		void setDrawPixelAccurate(bool drawPixelAccurate) { m_drawPixelAccurate = drawPixelAccurate; }

	private:
		// Helper functions for setting up the lattice stuff
		void setupEdgeColor();
		void setupLociValenceAndPointColor();
		void handleTLoci();
		void setupLocalSurfacesAndPatches();

		void addLatticeProperties();

		void addLocusOnVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace);
		void addLocusOnCornerVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace);
		void addLocusOnBoundaryVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace);
		void addLocusOnInnerVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace);

		void addLocus(OpenMesh::VertexHandle vertex, std::vector<Vec3f>& controlPoints,
			std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo>& faceMappings, Vec3f offset);

		std::vector<Vec3f> createLocalSurfaceControlPoints(
			Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight);
		std::vector<Vec3f> createLocalSurfaceControlPoints(
			Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
			Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
			Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight);

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
		void updateLatticeUniformBuffer();
		void updateMatrixUniformBuffer();
		void updatePatchUniformBuffer();



		// Buffers
		struct {
			VkBuffer buffer = VK_NULL_HANDLE;
			VkDeviceMemory memory = VK_NULL_HANDLE;
			int count;
		} m_pointsBuffer, m_linesBuffer, m_localSurfaceVertexBuffer, m_patchVertexBuffer;

		// Uniform buffer stuff
		struct {
			int tessInner = 5;
			int tessOuter = 5;
			int bFunctionIndex = 0;
			alignas(16) glm::mat4 projection = glm::mat4(1.0f);
			alignas(16) glm::mat4 modelview = glm::mat4(1.0f);
		} m_latticeUniforms;
		vks::Buffer m_latticeUniformBuffer;

		struct {
			alignas(16) std::vector<glm::mat4> matrices;
		} m_matrixUniforms;
		vks::Buffer m_matrixUniformBuffer;
		
		struct {
			std::vector<glm::vec4> controlPoints;
			std::vector<BoundaryInfo> boundaries;
		} m_patchUniforms;
		vks::Buffer m_patchUniformBuffer;

		// Pipelines and descriptor set stuff
		VkPipeline m_pointsPipeline;
		VkPipeline m_linesPipeline;
		VkPipeline m_localSurfacePipeline;
		VkPipeline m_localSurfaceWireframePipeline;

		VkPipelineLayout m_pipelineLayout;
		VkDescriptorSetLayout m_descriptorSetLayout;
		VkDescriptorSet m_descriptorSet;

		// Stuff passed from class creating the lattice, used for creating vulkan stuff
		VkDevice* m_device;
		vks::VulkanDevice* m_vulkanDevice;
		VkCommandPool* m_commandPool;
		VkDescriptorPool* m_descriptorPool;
		VkRenderPass* m_renderPass;
		VkQueue* m_queue;
		VkAllocationCallbacks* m_allocator;





		// Lattice stuff
		std::string m_name;

		bool m_draw = true;
		bool m_animate = false;
		bool m_drawLatticeGrid = true;
		bool m_drawLocalSurfaces = false;
		bool m_drawSurface = true;
		bool m_drawNormals = false;
		bool m_wireframe = false;
		bool m_drawPixelAccurate = false;

		uint32_t m_numControlPoints;
		uint32_t m_numLoci;
		uint32_t m_numPatches;

		std::vector<Locus> m_loci;
		std::vector<Patch> m_patches;

		std::vector<LocalSurfaceVertex> m_localSurfaceVertices;
		std::vector<LocalSurfaceVertex> m_patchVertices;
	};

}