#include "LocalSurfaceTexture.h"

namespace SWVL {

	LocalSurfaceTexture::LocalSurfaceTexture()
		: LocalSurfaceTexture(VK_NULL_HANDLE, nullptr, VK_NULL_HANDLE, VK_NULL_HANDLE, nullptr)
	{
	}

	LocalSurfaceTexture::LocalSurfaceTexture(VkDevice* device, vks::VulkanDevice* vulkanDevice, 
		VkCommandPool* commandPool, VkQueue* queue, VkAllocationCallbacks* allocator)
		: m_device(device), m_vulkanDevice(vulkanDevice), m_commandPool(commandPool), 
		m_queue(queue), m_allocator(allocator)
	{
	}

	LocalSurfaceTexture::~LocalSurfaceTexture()
	{
		if (m_device)
		{
			vkDestroyImageView(*m_device, m_imageView, m_allocator);
			vkDestroyImage(*m_device, m_image, m_allocator);
			vkDestroySampler(*m_device, m_sampler, m_allocator);
			vkFreeMemory(*m_device, m_memory, m_allocator);
		}
	}

	void LocalSurfaceTexture::loadBezier3x3(
		std::vector<glm::vec3> controlPoints, uint32_t numSamplesU, uint32_t numSamplesV)
	{
		m_width = numSamplesU;
		m_height = numSamplesV;
		m_layers = 2;

		glm::vec3& p00 = controlPoints[0], p10 = controlPoints[1], p20 = controlPoints[2];
		glm::vec3& p01 = controlPoints[3], p11 = controlPoints[4], p21 = controlPoints[5];
		glm::vec3& p02 = controlPoints[6], p12 = controlPoints[7], p22 = controlPoints[8];

		std::vector<glm::vec3> positions(numSamplesU * numSamplesV), normals(numSamplesU * numSamplesV);

		float du = 1.0f / (float)numSamplesU;
		float dv = 1.0f / (float)numSamplesV;

		for (size_t i = 0; i < numSamplesU; i++)
		{
			float u = du * i;
			glm::vec3 bu = bezBasis(u);
			glm::vec3 bud = bezBasis(u);
			for (size_t j = 0; j < numSamplesV; j++)
			{
				float v = dv * j;
				glm::vec3 bv = bezBasis(v);
				glm::vec3 bvd = bezBasisDer(v);

				glm::vec3 pos =
					p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] +
					p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] +
					p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2];

				glm::vec3 dpdu =
					p00 * bud[0] * bv[0] + p01 * bud[0] * bv[1] + p02 * bud[0] * bv[2] +
					p10 * bud[1] * bv[0] + p11 * bud[1] * bv[1] + p12 * bud[1] * bv[2] +
					p20 * bud[2] * bv[0] + p21 * bud[2] * bv[1] + p22 * bud[2] * bv[2];

				glm::vec3 dpdv =
					p00 * bu[0] * bvd[0] + p01 * bu[0] * bvd[1] + p02 * bu[0] * bvd[2] +
					p10 * bu[1] * bvd[0] + p11 * bu[1] * bvd[1] + p12 * bu[1] * bvd[2] +
					p20 * bu[2] * bvd[0] + p21 * bu[2] * bvd[1] + p22 * bu[2] * bvd[2];

				positions[i * numSamplesU + j] = pos;
				normals[i * numSamplesU + j] = glm::normalize(glm::cross(dpdu, dpdv));
			}
		}

		createImage(positions, normals);
	}

	glm::vec3 LocalSurfaceTexture::bezBasis(float t)
	{
		return glm::vec3(std::pow(1 - t, 2), 2 * t * (1 - t), std::pow(t, 2));
	}

	glm::vec3 LocalSurfaceTexture::bezBasisDer(float t)
	{
		return glm::vec3(2 * t - 2, 2 - 4 * t, 2 * t);
	}

	void LocalSurfaceTexture::createImage(std::vector<glm::vec3> positions, std::vector<glm::vec3> normals)
	{
		VkDeviceSize imageSize = m_width * m_height * sizeof(glm::vec3) * 2;

		positions.insert(positions.end(), normals.begin(), normals.end());

		VkMemoryAllocateInfo memAllocInfo = vks::initializers::memoryAllocateInfo();
		VkMemoryRequirements memReqs;

		// Create staging buffer
		VkBuffer stagingBuffer;
		VkDeviceMemory stagingMemory;

		VkBufferCreateInfo bufferCreateInfo = {};
		bufferCreateInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
		bufferCreateInfo.size = imageSize;
		bufferCreateInfo.usage = VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
		bufferCreateInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

		VK_CHECK_RESULT(vkCreateBuffer(*m_device, &bufferCreateInfo, m_allocator, &stagingBuffer));

		vkGetBufferMemoryRequirements(*m_device, stagingBuffer, &memReqs);

		memAllocInfo.allocationSize = memReqs.size;
		memAllocInfo.memoryTypeIndex = m_vulkanDevice->getMemoryType(memReqs.memoryTypeBits, 
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);

		VK_CHECK_RESULT(vkAllocateMemory(*m_device, &memAllocInfo, m_allocator, &stagingMemory));
		VK_CHECK_RESULT(vkBindBufferMemory(*m_device, stagingBuffer, stagingMemory, 0));


		// Transfer data to buffer
		float* data = nullptr;
		VK_CHECK_RESULT(vkMapMemory(*m_device, stagingMemory, 0, memReqs.size, 0, (void**)data));
		memcpy(data, positions.data(), imageSize);
		vkUnmapMemory(*m_device, stagingMemory);

		// Setup buffer copy regions for array layers
		std::vector<VkBufferImageCopy> bufferCopyRegions;

		for (uint32_t layer = 0; layer < m_layers; layer++)
		{
			// Setup a buffer image copy structure for the two layers
			uint32_t offset = m_width * m_height * sizeof(glm::vec3);
			VkBufferImageCopy bufferCopyRegion = {};
			bufferCopyRegion.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
			bufferCopyRegion.imageSubresource.mipLevel = 0;
			bufferCopyRegion.imageSubresource.baseArrayLayer = layer;
			bufferCopyRegion.imageSubresource.layerCount = 1;
			bufferCopyRegion.imageExtent = { m_width, m_height, 1 };
			bufferCopyRegion.bufferOffset = offset;
			bufferCopyRegions.push_back(bufferCopyRegion);
		}


		// Create image texture
		VkImageCreateInfo imageInfo = {};
		imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
		imageInfo.imageType = VK_IMAGE_TYPE_2D;
		imageInfo.format = VK_FORMAT_R32G32B32_SFLOAT;
		imageInfo.mipLevels = 1;
		imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
		imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
		imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
		imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
		imageInfo.extent = { m_width, m_height, 1 };
		imageInfo.usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
		imageInfo.arrayLayers = m_layers;

		VK_CHECK_RESULT(vkCreateImage(*m_device, &imageInfo, m_allocator, &m_image));

		vkGetImageMemoryRequirements(*m_device, m_image, &memReqs);
		
		// Allocate memory, use vulkanDevice*?
		VkMemoryAllocateInfo allocInfo = {};
		allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
		allocInfo.allocationSize = memReqs.size;
		allocInfo.memoryTypeIndex = m_vulkanDevice->getMemoryType(
			memReqs.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);

		VK_CHECK_RESULT(vkAllocateMemory(*m_device, &memAllocInfo, m_allocator, &m_memory));
		VK_CHECK_RESULT(vkBindImageMemory(*m_device, m_image, m_memory, 0));


		// Create a copy command
		VkCommandBuffer copyCmd;
		VkCommandBufferAllocateInfo cmdBufAllocateInfo =
			vks::initializers::commandBufferAllocateInfo(
				*m_commandPool,
				VK_COMMAND_BUFFER_LEVEL_PRIMARY,
				1);

		VK_CHECK_RESULT(vkAllocateCommandBuffers(*m_device, &cmdBufAllocateInfo, &copyCmd));

		VkCommandBufferBeginInfo cmdBufInfo = vks::initializers::commandBufferBeginInfo();
		VK_CHECK_RESULT(vkBeginCommandBuffer(copyCmd, &cmdBufInfo));

		// Image barrier for optimal image (target)
		// Set initial layout for all array layers (faces) of the optimal (target) tiled texture
		VkImageSubresourceRange subresourceRange = {};
		subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
		subresourceRange.baseMipLevel = 0;
		subresourceRange.levelCount = 1;
		subresourceRange.layerCount = m_layers;

		vks::tools::setImageLayout(copyCmd, m_image, VK_IMAGE_LAYOUT_UNDEFINED,
			VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, subresourceRange);

		vkCmdCopyBufferToImage(copyCmd, stagingBuffer, m_image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
			bufferCopyRegions.size(), bufferCopyRegions.data());

		// Change image layout to shader read after all layers have been copied
		vks::tools::setImageLayout(copyCmd, m_image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
			VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, subresourceRange);

		VK_CHECK_RESULT(vkEndCommandBuffer(copyCmd));

		VkSubmitInfo submitInfo = {};
		submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
		submitInfo.commandBufferCount = 1;
		submitInfo.pCommandBuffers = &copyCmd;

		VK_CHECK_RESULT(vkQueueSubmit(*m_queue, 1, &submitInfo, VK_NULL_HANDLE));
		VK_CHECK_RESULT(vkQueueWaitIdle(*m_queue));

		vkFreeCommandBuffers(*m_device, *m_commandPool, 1, &copyCmd);

		// Create Sampler
		VkSamplerCreateInfo samplerCreateInfo = {};
		samplerCreateInfo.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
		samplerCreateInfo.magFilter = VK_FILTER_LINEAR;
		samplerCreateInfo.minFilter = VK_FILTER_LINEAR;
		samplerCreateInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_LINEAR;
		samplerCreateInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
		samplerCreateInfo.addressModeV = samplerCreateInfo.addressModeU;
		samplerCreateInfo.addressModeW = samplerCreateInfo.addressModeU;
		samplerCreateInfo.mipLodBias = 0.0f;
		samplerCreateInfo.maxAnisotropy = 8;
		samplerCreateInfo.compareOp = VK_COMPARE_OP_NEVER;
		samplerCreateInfo.minLod = 0.0f;
		samplerCreateInfo.maxLod = 0.0f;
		samplerCreateInfo.borderColor = VK_BORDER_COLOR_FLOAT_OPAQUE_WHITE;
		VK_CHECK_RESULT(vkCreateSampler(*m_device, &samplerCreateInfo, m_allocator, &m_sampler));

		// Create image view
		VkImageViewCreateInfo viewCreateInfo = {};
		viewCreateInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
		viewCreateInfo.format = VK_FORMAT_R32G32B32_SFLOAT;
		viewCreateInfo.components = { VK_COMPONENT_SWIZZLE_R, VK_COMPONENT_SWIZZLE_G, VK_COMPONENT_SWIZZLE_B };
		viewCreateInfo.subresourceRange = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1 };
		viewCreateInfo.subresourceRange.layerCount = m_layers;
		viewCreateInfo.subresourceRange.levelCount = 1;
		viewCreateInfo.image = m_image;
		VK_CHECK_RESULT(vkCreateImageView(*m_device, &viewCreateInfo, m_allocator, &m_imageView));

		// Clean up staging resources
		vkFreeMemory(*m_device, stagingMemory, m_allocator);
		vkDestroyBuffer(*m_device, stagingBuffer, m_allocator);
	}

}