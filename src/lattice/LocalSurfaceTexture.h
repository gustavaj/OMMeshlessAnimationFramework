#pragma once

#include <vulkan/vulkan.h>

#include <glm/glm.hpp>

#include <vector>

#include "../vulkan/VulkanDevice.hpp"

// https://github.com/SaschaWillems/Vulkan/blob/master/examples/texturearray/texturearray.cpp

namespace SWVL {

	class LocalSurfaceTexture
	{
	public:
		LocalSurfaceTexture();
		LocalSurfaceTexture(VkDevice* device, vks::VulkanDevice* vulkanDevice, 
			VkCommandPool* commandPool, VkQueue* queue, VkAllocationCallbacks* allocator);
		~LocalSurfaceTexture();

		void loadBezier3x3(std::vector<glm::vec3> controlPoints, uint32_t numSamplesU, uint32_t numSamplesV);

		const VkImage& image() { return m_image; }
		const VkSampler& sampler() { return m_sampler; }
		const VkImageView& imageView() { return m_imageView; }

	private:
		glm::vec3 bezBasis(float t);
		glm::vec3 bezBasisDer(float t);

		void createImage(std::vector<glm::vec3> positions, std::vector<glm::vec3> normals);

		uint32_t m_index;
		uint32_t m_width;
		uint32_t m_height;
		uint32_t m_layers;

		VkImage m_image;
		VkDeviceMemory m_memory;
		VkSampler m_sampler;
		VkImageView m_imageView;

		// Vulkan pointers
		VkDevice* m_device;
		vks::VulkanDevice* m_vulkanDevice;
		VkCommandPool* m_commandPool;
		VkQueue* m_queue;
		VkAllocationCallbacks* m_allocator;
	};

}