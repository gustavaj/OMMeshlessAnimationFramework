#include <vulkan/vulkan.h>

// Sascha willems
#include "../vulkan/VulkanDevice.hpp"
#include "../vulkan/VulkanBuffer.hpp"
#include "../vulkan/VulkanUIOverlay.h"

#include "Lattice.h"

namespace SWVL
{
	struct Buffer {
		Buffer() : buffer(VK_NULL_HANDLE), memory(VK_NULL_HANDLE), count(0) {}

		VkBuffer buffer;
		VkDeviceMemory memory;
		int count = 0;
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

	// Struct containing info for a vertex used for points and lines in the lattice grid.
	struct GridVertex
	{
		OML::Vec3f pos;
		OML::Col3 col;

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

	const std::vector<std::string> BFunctionNames = {
		"3x^2-2x^3", "6x^5-15x^4+10x^3", "1/(1+e^(1/x-1/(1-x)))"
	};

	class SWVulkanLattice : public OML::Lattice
	{
	public:
		SWVulkanLattice();
		SWVulkanLattice(std::string name);
		~SWVulkanLattice();

		// Vulkan framework
		/* Create buffers and pipelines and stuff for Vulkan */
		void initVulkanStuff(
			VkDevice* device, vks::VulkanDevice* vulkanDevice,
			VkQueue* queue, VkCommandPool* commandPool,
			VkDescriptorPool* descriptorPool, VkRenderPass* renderPass,
			VkAllocationCallbacks* allocator);
		/* Clean up vulkan stuff */
		void destroyVulkanStuff();
		/* Add commands before renderpass has started, like e.g. resetting query pools */
		void addToCommandbufferPreRenderpass(VkCommandBuffer& commandBuffer);
		/* Add drawing commands to a command buffer */
		void addToCommandbuffer(VkCommandBuffer& commandBuffer);
		/* Update projection and view matrices */
		void onViewChanged(glm::mat4 projection, glm::mat4 view);

		/* Update UIOverlay with information about the lattice */
		bool onUpdateUIOverlay(vks::UIOverlay* overlay);

	protected:
		virtual void setupLocalSurfaceVertex(OML::Locus& locus) override;
		virtual void setupPatchVertices(OML::Patch& patch) override;
		virtual void localUpdate(double dt) override;

	private:
		// Vulkan functions
		void createDeviceLocalBuffer(
			VkBuffer& buffer, VkDeviceMemory& memory, void* data,
			uint32_t bufferSize, VkBufferUsageFlagBits usage);
		void createBuffers();
		void prepareUniformBuffers();
		void uploadStorageBuffers();
		void setupDescriptorSetLayouts();
		VkPipelineShaderStageCreateInfo loadShader(std::string fileName, VkShaderStageFlagBits stage);
		void preparePipelines();
		void setupDescriptorPool();
		void setupDescriptorSets();
		void updateLatticeUniformBuffer();
		void updateMatrixUniformBuffer();
		void setupQueryResultBuffer();
		void getQueryResults();

		// Buffers
		Buffer m_pointsBuffer;
		Buffer m_linesBuffer;
		Buffer m_localSurfaceVertexBuffer;
		Buffer m_patchVertexBuffer;
		Buffer m_queryResult;
		Buffer m_timingResult;

		VkQueryPool m_queryPool = VK_NULL_HANDLE;
		VkQueryPool m_timingPool = VK_NULL_HANDLE;
		bool m_doPipelineTimings = false;
		bool m_doPipelineQueries = false;
		uint64_t m_pipelineStats[10] = { 0 };
		uint64_t m_pipelineTimings[2] = { 0 };
		double m_timestampPeriod = 1.0;

		// Uniform buffers
		vks::Buffer m_latticeUniformBuffer;
		vks::Buffer m_matrixUniformBuffer;
		vks::Buffer m_controlPointBuffer;
		vks::Buffer m_boundariesBuffer;

		// Pipelines and descriptor set stuff
		VkPipeline m_pointsPipeline;
		VkPipeline m_linesPipeline;
		VkPipeline m_localSurfacePipeline;
		VkPipeline m_localSurfaceWireframePipeline;
		VkPipeline m_patchPipeline;
		VkPipeline m_patchWireframePipeline;
		VkPipeline m_normalPipeline;

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

		std::vector<LocalSurfaceVertex> m_localSurfaceVertices;
		std::vector<LocalSurfaceVertex> m_patchVertices;

		std::vector<std::string> m_listItems;
		int m_selectedSurface;

		bool m_destroyed = false;
	};
}