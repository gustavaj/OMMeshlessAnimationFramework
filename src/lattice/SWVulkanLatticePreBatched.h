#pragma once

#include "SWVulkanLattice.h"


#include "LocalSurfaceTextureBatch.h"

/*
	If running with debug layers turned on there will be a large CPU overhead from
	the layer checking if the images are in the correct layout.


	NEW PLAN:
	-Bigger textures!
	-Do something like 100(1000?) patches per texture (3 layers), batch draw calls
	-For every patch, evaluate the local surfaces for that patch, i.e. ls1 from [0-0.5]^2
	-Only buffer needed is matrices, maybe add normal matrices too?
*/

namespace SWVL
{
	// !! Should probably not try different samling rates.
	const int NUM_SAMPLES_U_BATCH = 8;
	const int NUM_SAMPLES_V_BATCH = 8;
	const int NUM_PATCHES_PER_BATCH = 1024;
	const int BATCH_ROWS = 32;
	const int BATCH_COLS = 32;
	const int MAX_BATCHES = 4;

	class SWVulkanLatticePreBatched : public OML::Lattice
	{
	public:
		SWVulkanLatticePreBatched();
		SWVulkanLatticePreBatched(std::string name);
		~SWVulkanLatticePreBatched();

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
		/* Update windowSize uniform */
		void onWindowResized(float width, float height);

		/* Update UIOverlay with information about the lattice */
		bool onUpdateUIOverlay(vks::UIOverlay* overlay);

		/*
		Takes a VkPhysicalDeviceFeatures structure of the devices' capabilities, and checks if all the
		required features are supported. Writes VK_TRUE to all the required features in the enabledFeatures struct.
		*/
		static void CheckAndSetupRequiredPhysicalDeviceFeatures(
			VkPhysicalDeviceFeatures& deviceFeatures, VkPhysicalDeviceFeatures& enabledFeatures);

	protected:
		// Called from the update function in the the Lattice class.
		virtual void localUpdate(double dt) override;

	private:
		// Boolean to prevent vulkan stuff to be initiated more than once.
		bool m_vulkanInitiated = false;
		// Boolean to prevent vulkan stuff from being destroyed if it has not been initiated first.
		bool m_destroyed = true;

		// Menu variables collected here--
		// Suffix to use on all menu items to make their id unique
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

		VkPipelineLayout m_pipelineLayout;
		VkDescriptorSetLayout m_descriptorSetLayout;
		VkDescriptorSet m_descriptorSet;

		VkDescriptorSetLayout m_batchDescSetLayout;
		std::vector<VkDescriptorSet> m_batchDescSets;

		std::vector<LocalSurfaceTextureBatch> m_batchTextures;

		struct BatchInfo {
			uint32_t index;
			uint32_t count;
		};
		std::vector<BatchInfo> m_batchInfo;


		// Map for holding the created shader modules.
		std::unordered_map<std::string, VkShaderModule> m_shaderModules;

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