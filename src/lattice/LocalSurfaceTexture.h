#pragma once

#include <vulkan/vulkan.h>

#include <glm/glm.hpp>

#include <vector>

#include "../vulkan/VulkanDevice.hpp"
#include "../vulkan/VulkanTexture.hpp"

#include "LatticeUtility.h"

// https://github.com/SaschaWillems/Vulkan/blob/master/examples/texturearray/texturearray.cpp

namespace OML {

	class LocalSurfaceTexture
	{
	public:
		LocalSurfaceTexture();
		LocalSurfaceTexture(VkDevice* device, vks::VulkanDevice* vulkanDevice,
			VkCommandPool* commandPool, VkQueue* queue, VkAllocationCallbacks* allocator);
		~LocalSurfaceTexture();

		void loadLocalSurface(std::vector<glm::vec3>& controlPoints,
			uint32_t numSamplesU, uint32_t numSamplesV, LocalSurfaceType lsType);
		void loadBezier3x3(std::vector<glm::vec3>& controlPoints, uint32_t numSamplesU, uint32_t numSamplesV);
		void loadBezier4x4(std::vector<glm::vec3>& controlPoints, uint32_t numSamplesU, uint32_t numSamplesV);
		void loadPlane(std::vector<glm::vec3>& controlPoints, uint32_t numSamplesU, uint32_t numSamplesV);

		void destroy();

		uint32_t memoryUsage() const { return m_memoryUsage; }

		const VkImage& image() { return m_texture.image; }
		const VkSampler& sampler() { return m_texture.sampler; }
		const VkImageView& imageView() { return m_texture.view; }
		VkDescriptorImageInfo* descriptor() { return &m_texture.descriptor; }

	private:
		glm::vec3 bezBasis3(float t);
		glm::vec3 bezBasisDer3(float t);
		glm::vec4 bezBasis4(float t);
		glm::vec4 bezBasisDer4(float t); 
		inline glm::vec3 mix(glm::vec3 a, glm::vec3 b, float t) {
			return (1.0f - t) * a + (t * b);
		}

		void createImage(std::vector<glm::vec4> samples);

		uint32_t m_index;
		uint32_t m_width;
		uint32_t m_height;
		uint32_t m_layers;

		uint32_t m_memoryUsage = 0;

		VkImage m_image;
		VkDeviceMemory m_memory;
		VkSampler m_sampler;
		VkImageView m_imageView;
		VkDescriptorImageInfo m_descriptor;

		vks::Texture2DArray m_texture;

		// Vulkan pointers
		VkDevice* m_device;
		vks::VulkanDevice* m_vulkanDevice;
		VkCommandPool* m_commandPool;
		VkQueue* m_queue;
		VkAllocationCallbacks* m_allocator;

		std::vector<glm::vec4> m_mapped;
	};

}