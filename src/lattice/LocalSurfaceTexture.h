#pragma once

#include <vulkan/vulkan.h>

#include <glm/glm.hpp>

#include <unordered_map>
#include <vector>

#include "../vulkan/VulkanDevice.hpp"
#include "../vulkan/VulkanTexture.hpp"

#include "Lattice.h"

// https://github.com/SaschaWillems/Vulkan/blob/master/examples/texturearray/texturearray.cpp

/*
	TODO:
	-Make texture much bigger, like 4096x4096 or something
	-Add several local surfaces to the same image
*/

namespace SWVL {

	// indices in x and y direction plus width and height of the part
	struct ImageIndex {
		float x;
		float y;
		float w;
		float h;
	};

	class LocalSurfaceTexture
	{
	public:
		LocalSurfaceTexture();
		LocalSurfaceTexture(/*uint32_t width, uint32_t height, uint32_t numImgU, uint32_t numImgV,*/
			VkDevice* device, vks::VulkanDevice* vulkanDevice, 
			VkCommandPool* commandPool, VkQueue* queue, VkAllocationCallbacks* allocator);
		~LocalSurfaceTexture();

		void addPatch(std::vector<)

		void loadBezier3x3(std::vector<glm::vec3> controlPoints, uint32_t numSamplesU, uint32_t numSamplesV);

		void destroy();

		// Returns true if the given sample size can fit into this image
		bool isFull() { return m_rows == m_curX && m_cols == m_curY; }

		const VkImage& image() { return m_texture.image; }
		const VkSampler& sampler() { return m_texture.sampler; }
		const VkImageView& imageView() { return m_texture.view; }
		VkDescriptorImageInfo* descriptor() { return &m_texture.descriptor; }

	private:
		glm::vec3 bezBasis(float t);
		glm::vec3 bezBasisDer(float t);

		void createImage(std::vector<glm::vec4> samples);

		// Map from a pair of control point index and face index to a struct containing 
		// the index for the data in the image.
		std::unordered_map<uint32_t, ImageIndex> m_map;

		uint32_t m_index;

		uint32_t m_width;
		uint32_t m_height;
		uint32_t m_layers = 3;
		
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