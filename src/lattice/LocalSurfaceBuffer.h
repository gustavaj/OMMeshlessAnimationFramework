#pragma once

#include <glm/glm.hpp>

#include <vector>

namespace OML {

	struct LocalSurfaceData {
		std::vector<glm::vec4> positions;
		std::vector<glm::vec4> partialU;
		std::vector<glm::vec4> partialV;
	};

	class LocalSurfaceBuffer {
	public:
		LocalSurfaceBuffer() : LocalSurfaceBuffer(16, 16) {}
		LocalSurfaceBuffer(uint32_t numSamplesU, uint32_t numSamplesV)
			: m_numSamplesU(numSamplesU), m_numSamplesV(numSamplesV) {}

		uint32_t addBezier3x3(std::vector<glm::vec3> controlPoints);

		uint32_t bufferSize() 
		{
			return m_data.size() * m_numSamplesU * m_numSamplesV * 3 * sizeof(glm::vec4);
		}

		std::vector<glm::vec4> getData();
		
	private:
		uint32_t m_numSamplesU;
		uint32_t m_numSamplesV;

		std::vector<LocalSurfaceData> m_data;

		glm::vec3 bezBasis(float t);
		glm::vec3 bezBasisDer(float t);
	};

}