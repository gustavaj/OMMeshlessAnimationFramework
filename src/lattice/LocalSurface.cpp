#include "LocalSurface.h"


namespace OML
{
	std::vector<VkVertexInputBindingDescription> Bezier3x3Vertex::GetBindingDescription()
	{
		std::vector<VkVertexInputBindingDescription> bindings;

		bindings.resize(1);
		bindings[0].binding = 0;
		bindings[0].stride = sizeof(Bezier3x3Vertex);
		bindings[0].inputRate = VK_VERTEX_INPUT_RATE_VERTEX;

		return bindings;
	}

	std::vector<VkVertexInputAttributeDescription> Bezier3x3Vertex::GetAttributeDescriptions()
	{
		std::vector<VkVertexInputAttributeDescription> attributeDescriptions = {};
		attributeDescriptions.resize(13);
		for (size_t i = 0; i < 13; i++) {
			attributeDescriptions[i].binding = 0;
			attributeDescriptions[i].location = i;
			attributeDescriptions[i].format = VK_FORMAT_R32G32B32_SFLOAT;
		}
		attributeDescriptions[0].offset = offsetof(Bezier3x3Vertex, p00);
		attributeDescriptions[1].offset = offsetof(Bezier3x3Vertex, p10);
		attributeDescriptions[2].offset = offsetof(Bezier3x3Vertex, p20);
		attributeDescriptions[3].offset = offsetof(Bezier3x3Vertex, p01);
		attributeDescriptions[4].offset = offsetof(Bezier3x3Vertex, p11);
		attributeDescriptions[5].offset = offsetof(Bezier3x3Vertex, p21);
		attributeDescriptions[6].offset = offsetof(Bezier3x3Vertex, p02);
		attributeDescriptions[7].offset = offsetof(Bezier3x3Vertex, p12);
		attributeDescriptions[8].offset = offsetof(Bezier3x3Vertex, p22);
		attributeDescriptions[9].offset = offsetof(Bezier3x3Vertex, trans);
		attributeDescriptions[10].offset = offsetof(Bezier3x3Vertex, rot);
		attributeDescriptions[11].offset = offsetof(Bezier3x3Vertex, scale);
		attributeDescriptions[12].format = VK_FORMAT_R32G32_SFLOAT;
		attributeDescriptions[12].offset = offsetof(Bezier3x3Vertex, boundaryUV);


		return attributeDescriptions;
	}

	Bezier3x3Vertex Bezier3x3Vertex::CreateFlatSurface(Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight, Vec3f offset, Vec2f boundaryUV)
	{
		Vec3f topLeftToRight = topRight - topLeft;
		Vec3f topToBotLeft = bottomLeft - topLeft;
		Vec3f topToBotRight = bottomRight - topRight;
		Vec3f botLeftToRight = bottomRight - bottomLeft;
		return Bezier3x3Vertex(
			topLeft,												// p00
			topLeft + topLeftToRight * 0.5,							// p10
			topRight,												// p20
			topLeft + topToBotLeft * 0.5,							// p01
			topLeft + topLeftToRight * 0.5 + topToBotLeft * 0.5,	// p11
			topRight + topToBotRight * 0.5,							// p21
			bottomLeft,												// p02
			bottomLeft + botLeftToRight * 0.5,						// p12
			bottomRight,											// p22
			offset,													// trans
			Vec3f(0.0f, 0.0f, 0.0f),								// rot
			Vec3f(1.0f, 1.0f, 1.0f),								// scale
			boundaryUV
		);
	}

	Bezier3x3Vertex Bezier3x3Vertex::CreateFlatSurface(Vec3f topLeft, Vec3f topMiddle, Vec3f topRight, Vec3f middleLeft, Vec3f middle, Vec3f middleRight, Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight, Vec3f offset, Vec2f boundaryUV)
	{
		return Bezier3x3Vertex(
			topLeft,					// p00
			topMiddle,					// p10
			topRight,					// p20
			middleLeft,					// p01
			middle,						// p11
			middleRight,				// p21
			bottomLeft,					// p02
			bottomMiddle,				// p12
			bottomRight,				// p22
			offset,						// trans
			Vec3f(0.0f, 0.0f, 0.0f),	// rot
			Vec3f(1.0f, 1.0f, 1.0f),	// scale
			boundaryUV
		);
	}
}

