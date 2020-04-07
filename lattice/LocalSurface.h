#pragma once

#include <vulkan/vulkan.h>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <vector>

namespace OML
{
	using Vec3f = OpenMesh::Vec3f;
	using Vec2f = OpenMesh::Vec2f;

	struct Bezier3x3Vertex
	{
		Bezier3x3Vertex() {}
		Bezier3x3Vertex(Vec3f p00, Vec3f p10, Vec3f p20,
			Vec3f p01, Vec3f p11, Vec3f p21,
			Vec3f p02, Vec3f p12, Vec3f p22,
			Vec3f trans, Vec3f rot, Vec3f scale, Vec2f boundaryUV)
			: p00(p00), p10(p10), p20(p20), 
			  p01(p01), p11(p11), p21(p21),
			  p02(p02), p12(p12), p22(p22), 
			  trans(trans), rot(rot), scale(scale), boundaryUV(boundaryUV) {}

		Vec3f p00, p10, p20;
		Vec3f p01, p11, p21;
		Vec3f p02, p12, p22;
		Vec3f trans, rot, scale;
		Vec2f boundaryUV;

		static std::vector<VkVertexInputBindingDescription> GetBindingDescription();
		static std::vector<VkVertexInputAttributeDescription> GetAttributeDescriptions();

		static Bezier3x3Vertex CreateFlatSurface(Vec3f topLeft, Vec3f topRight,
			Vec3f bottomLeft, Vec3f bottomRight,
			Vec3f offset, Vec2f boundaryUV);
		static Bezier3x3Vertex CreateFlatSurface(
			Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
			Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
			Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight,
			Vec3f offset, Vec2f boundaryUV);
	};

	class LocalSurface
	{
	public:

		LocalSurface() {}
		LocalSurface(Vec3f topLeft, Vec3f topRight,
			Vec3f bottomLeft, Vec3f bottomRight, Vec3f offset, Vec2f boundaryUV)
		{
			m_vertex = Bezier3x3Vertex::CreateFlatSurface(
				topLeft, topRight,
				bottomLeft, bottomRight, offset, boundaryUV);
		}
		LocalSurface(
			Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
			Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
			Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight, Vec3f offset, Vec2f boundaryUV)
		{
			m_vertex = Bezier3x3Vertex::CreateFlatSurface(
				topLeft, topMiddle, topRight,
				middleLeft, middle, middleRight,
				bottomLeft, bottomMiddle, bottomRight, offset, boundaryUV);
		}
		~LocalSurface() {}

		Bezier3x3Vertex& vertex() { return m_vertex; }

		float getStartU() { return (m_vertex.p01 - m_vertex.p11).length() / (m_vertex.p21 - m_vertex.p01).length(); }
		float getEndU() { return (m_vertex.p01 - m_vertex.p11).length() / (m_vertex.p21 - m_vertex.p01).length(); }
		float getStartV() { return (m_vertex.p10 - m_vertex.p11).length() / (m_vertex.p12 - m_vertex.p10).length(); }
		float getEndV() { return (m_vertex.p10 - m_vertex.p11).length() / (m_vertex.p12 - m_vertex.p10).length(); }

		uint32_t vertexSize() { return sizeof(Bezier3x3Vertex); }

	protected:
		Bezier3x3Vertex m_vertex;
	};
}