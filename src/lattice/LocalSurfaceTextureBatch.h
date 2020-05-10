#pragma once

#include <vulkan/vulkan.h>

#include <glm/glm.hpp>

#include <unordered_map>
#include <vector>

#include "../vulkan/VulkanDevice.hpp"
#include "../vulkan/VulkanTexture.hpp"

#include "LatticeUtility.h"

// https://github.com/SaschaWillems/Vulkan/blob/master/examples/texturearray/texturearray.cpp

namespace OML {

	class LocalSurfaceTextureBatch
	{
	public:
		LocalSurfaceTextureBatch();
		LocalSurfaceTextureBatch(uint32_t rows, uint32_t cols, uint32_t numSamplesU, uint32_t numSamplesV,
			VkDevice* device, vks::VulkanDevice* vulkanDevice, 
			VkCommandPool* commandPool, VkQueue* queue, VkAllocationCallbacks* allocator);
		~LocalSurfaceTextureBatch();

		std::pair<uint32_t, uint32_t> addPatch(
			std::vector<glm::vec3> p00Points, OML::BoundaryInfo& p00Boundary,
			std::vector<glm::vec3> p10Points, OML::BoundaryInfo& p10Boundary,
			std::vector<glm::vec3> p01Points, OML::BoundaryInfo& p01Boundary,
			std::vector<glm::vec3> p11Points, OML::BoundaryInfo& p11Boundary,
			LocalSurfaceType type
			);

		void allocateMemory();

		void destroy();

		uint32_t memoryUsage() const { return m_memoryUsage; }

		// Returns true if the given sample size can fit into this image
		bool isFull() { return m_rows == m_curX && m_cols == m_curY; }

		const VkImage& image() { return m_texture.image; }
		const VkSampler& sampler() { return m_texture.sampler; }
		const VkImageView& imageView() { return m_texture.view; }
		VkDescriptorImageInfo* descriptor() { return &m_texture.descriptor; }

	private:
		std::vector<glm::vec4> m_data;

		inline float mix(float a, float b, float t) {
			return (1.0f - t) * a + (t * b);
		}
		inline glm::vec3 mix(glm::vec3 a, glm::vec3 b, float t) {
			return (1.0f - t) * a + (t * b);
		}

		glm::vec3 bezBasis3(float t);
		glm::vec3 bezBasisDer3(float t);
		glm::vec4 bezBasis4(float t);
		glm::vec4 bezBasisDer4(float t);

		void createImage();

		void loadBezier3x3(std::vector<glm::vec3>& controlPoints, 
			OML::BoundaryInfo& boundary, uint32_t baseLayer);
		void loadBezier4x4(std::vector<glm::vec3>& controlPoints,
			BoundaryInfo& boundary, uint32_t baseLayer);
		void loadPlane(std::vector<glm::vec3>& controlPoints,
			BoundaryInfo& boundary, uint32_t baseLayer);

		uint32_t m_numSamplesU;
		uint32_t m_numSamplesV;
		uint32_t m_layers = 12;

		uint32_t m_memoryUsage = 0;
		
		uint32_t m_rows;
		uint32_t m_cols;
		uint32_t m_curX = 0;
		uint32_t m_curY = 0;

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