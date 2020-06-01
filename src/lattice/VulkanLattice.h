#pragma once

#include <vulkan/vulkan.h>

// TODO: Remove all dependencies on Sascha Willems' framework. Then rename.. Maybe not.
// Sascha willems
#include "../vulkan/VulkanDevice.hpp"
#include "../vulkan/VulkanBuffer.hpp"
#include "../vulkan/VulkanUIOverlay.h"

#include "Lattice.h"
#include "LocalSurfaceTexture.h"
#include "LocalSurfaceTextureBatch.h"
#include "LocalSurfaceBuffer.h"

//#define DIRECT_MULTIPLE_DRAW_CALLS
//#define ADD_DUMMY_DATA_TO_CONTROL_POINT_BUFFER
//const size_t DUMMY_DATA_SIZE_BYTES = 500000000;

namespace OML {

	// Parameters used for pre-evaluation
	const int NUM_SAMPLES_U = 16;
	const int NUM_SAMPLES_V = 16;
	const int NUM_PATCHES_PER_BATCH = 1024;
	const int BATCH_ROWS = 32;
	const int BATCH_COLS = 32;
	const int MAX_BATCHES = 1;

	// Various structs used by the vulkan implementation

	// Simple class holding the data used for a device local vulkan buffer.
	struct Buffer {
		Buffer() : buffer(VK_NULL_HANDLE), memory(VK_NULL_HANDLE), count(0) {}

		VkBuffer buffer;
		VkDeviceMemory memory;
		int count = 0;
	};

	/*
	Holds the information for one local surface that is sent to the vertex shader
	A local surface is defined by a vec4 in the shader, by the index and count of their
	control points from the m_controlPoints vector and the indices of the transformation matrix
	and in the m_matrices vector and the BoundaryInfo in the m_boundaries vector. Also the color of the patch.

	For rendering local surfaces the TCS takes one vertex per patch, and for a patch(Block) it takes 4 vertices per patch
	*/
	struct LocalSurfaceVertex
	{
		LocalSurfaceVertex() : LocalSurfaceVertex(0, 0, 0, 0, glm::vec3(1.0f)) {}
		LocalSurfaceVertex(uint32_t controlPointIndex, uint32_t controlPointCount,
						   uint32_t matrixIndex, uint32_t boundaryIndex, glm::vec3 color,
						   uint32_t s = 0, uint32_t t = 0, uint32_t numSamples = 0, uint32_t dataIndex = 0)
			: controlPointIndex(controlPointIndex), controlPointCount(controlPointCount),
			  matrixIndex(matrixIndex), boundaryIndex(boundaryIndex), color(color),
			  s(s), t(t), numSamples(0), dataIndex(dataIndex) {}

		uint32_t controlPointIndex;
		uint32_t controlPointCount; 
		uint32_t matrixIndex;
		uint32_t boundaryIndex;
		glm::vec3 color;
		uint32_t s;
		uint32_t t;
		uint32_t numSamples;
		uint32_t dataIndex;

		static std::vector<VkVertexInputBindingDescription> GetBindingDescriptions();
		static std::vector<VkVertexInputAttributeDescription> GetAttributeDescriptions();
	};

	// Struct containing info for a vertex used for points and lines in the lattice grid.
	struct GridVertex
	{
		OML::Vec3f pos;
		OML::Col3 col;

		static VkVertexInputBindingDescription GetBindingDescription();
		static std::vector<VkVertexInputAttributeDescription> GetAttributeDesctiptions();
	};

	// Pointers to the four samplers used for the local surfaces of a given patch
	struct PatchSamplerInfo
	{
		LocalSurfaceTexture* p00Sampler;
		LocalSurfaceTexture* p10Sampler;
		LocalSurfaceTexture* p01Sampler;
		LocalSurfaceTexture* p11Sampler;
	};


	struct BatchInfo {
		uint32_t index;
		uint32_t count;
	};




	class VulkanLattice : public OML::Lattice
	{
	public:
		VulkanLattice();
		VulkanLattice(std::string name, 
			LocalSurfaceType lsType = LocalSurfaceType::Quadratic_Bezier,
			EvaluationMethod evalMethod = EvaluationMethod::Direct);
		~VulkanLattice();

		// Vulkan framework
		void initVulkan(
			VkDevice* device, vks::VulkanDevice* vulkanDevice,
			VkQueue* queue, VkCommandPool* commandPool,
			VkDescriptorPool* descriptorPool, VkRenderPass* renderPass,
			VkAllocationCallbacks* allocator);
		/* Clean up vulkan stuff */
		void destroyVulkan();
		/* Add commands before renderpass has started, like e.g. resetting query pools */
		void addToCommandbufferPreRenderpass(VkCommandBuffer& commandBuffer);
		/* Add drawing commands to a command buffer */
		void addToCommandbuffer(VkCommandBuffer& commandBuffer);
		/* Update projection and view matrices */
		void onViewChanged(glm::mat4 projection, glm::mat4 view);
		/* Update windowSize uniform */
		void onWindowResized(float width, float height);

		/* Update UIOverlay with information about the lattice */
		bool onUpdateUIOverlay(vks::UIOverlay* overlay);

		// Sets the evaluation method to be used, must be called before initVulkan() to have any effect
		void setEvaluationMethod(EvaluationMethod evalMethod) { m_evalMethod = evalMethod; }

		/*
		Takes a VkPhysicalDeviceFeatures structure of the devices' capabilities, and checks if all the
		required features are supported. Writes VK_TRUE to all the required features in the enabledFeatures struct.
		*/
		static void CheckAndSetupRequiredPhysicalDeviceFeatures(
			VkPhysicalDeviceFeatures& deviceFeatures, VkPhysicalDeviceFeatures& enabledFeatures);


		/* Set tessellation factors */
		void setTessellationFactors(int inner, int outer) {
			m_uniforms.tessInner = inner; 
			m_uniforms.tessOuter = outer;
		}
		void updateUniformBuffer() {
			updateLatticeUniformBuffer();
		}
		void updateMatrixBuffer() {
			updateMatrixUniformBuffer();
		}

	protected:
		// Called from the update function in the the Lattice class.
		virtual void localUpdate(double dt) override;

	private:
		EvaluationMethod m_evalMethod;
		// Boolean to prevent vulkan stuff to be initiated more than once.
		bool m_vulkanInitiated = false;
		// Boolean to prevent vulkan stuff from being destroyed if it has not been initiated first.
		bool m_destroyed = true;

		// Note: The size of the memory usage is equal to the size of the data stored in GPU memory
		// The size of the memory allocations may differ.
		uint64_t m_deviceMemoryUsage = 0;

		// Menu variables collected here--
		// Suffix to use on all menu items to make their id unique
		// Used so menu controls are responsive even when multiple Lattices are 
		// Created
		std::string m_menuSuffix;
		// Index of curretn simulator in menu
		int m_simulatorIndex = 0;
		// Holds the names used in the ImGui ui for the local surfaces.
		std::vector<std::string> m_listItems;
		// Index of current surface in menu
		int m_selectedSurface = 0;


		// Create a buffer that is local on the GPU
		void createDeviceLocalBuffer(
			VkBuffer& buffer, VkDeviceMemory& memory, void* data,
			uint32_t bufferSize, VkBufferUsageFlagBits usage);
		// Setup the vertices used for rendering. Grid, local surfaces, surfaces.
		void setupVertices();
		void createBuffers();
		void prepareUniformBuffers();
		void uploadStorageBuffers();
		void setupDescriptorSetLayouts();
		// Helper function to load a shader.
		VkPipelineShaderStageCreateInfo loadShader(std::pair<std::string, std::vector<uint32_t>&> src, VkShaderStageFlagBits stage);
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
		// Enables/disables doing timings of the render pass
		bool m_doPipelineTimings = false;
		// Enables/disables doing queries for the render pass.
		bool m_doPipelineQueries = false;
		uint64_t m_pipelineStats[10] = { 0 };
		uint64_t m_pipelineTimings[2] = { 0 };
		// The time it takes for the gpu to increase the timestamp in nanoseconds.
		double m_timestampPeriod = 1.0;

		// Uniform/Storage buffers
		vks::Buffer m_latticeUniformBuffer;
		vks::Buffer m_matrixUniformBuffer;
		vks::Buffer m_controlPointBuffer;
		vks::Buffer m_boundariesBuffer;
		vks::Buffer m_localSurfaceDataBuffer;

		// Pipelines and descriptor set stuff
		VkPipeline m_pointsPipeline;
		VkPipeline m_linesPipeline;
		VkPipeline m_localSurfacePipeline;
		VkPipeline m_localSurfaceWireframePipeline;
		VkPipeline m_patchPipeline;
		VkPipeline m_patchWireframePipeline;
		VkPipeline m_normalPipeline;
		VkPipeline m_displaySurfaceAccuracyPipeline;
		VkPipeline m_displayPixelAccuracyPipeline;
		VkPipeline m_displayTriSizePipeline;

		VkPipelineLayout m_pipelineLayout;
		VkDescriptorSetLayout m_descriptorSetLayout;
		VkDescriptorSet m_descriptorSet;

		// Stuff that is specific to a certain evaluation method
		struct {
			VkDescriptorSetLayout setLayout;
			std::vector<VkDescriptorSet> sets;
		} m_localSamplerDescriptor, m_samplerDescriptor, m_batchDescriptor;

		std::unordered_map<uint32_t, LocalSurfaceTexture> m_localSurfaceTextures;
		std::vector<PatchSamplerInfo> m_patchSamplers;

		LocalSurfaceBuffer m_LSBuffer;
		std::unordered_map<uint32_t, uint32_t> m_LSIdxToLSBufferMap;

		// Map for holding the created shader modules.
		std::unordered_map<std::string, VkShaderModule> m_shaderModules;

		std::vector<LocalSurfaceTextureBatch> m_batchTextures;
		std::vector<BatchInfo> m_batchInfo;

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
	};

}