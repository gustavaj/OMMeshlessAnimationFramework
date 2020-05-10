#pragma once

#include <glm/glm.hpp>

#include <vector>

#include "LatticeUtility.h"

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

		uint32_t addLocalSurface(std::vector<glm::vec3>& controlPoints, LocalSurfaceType lsType);
		uint32_t addBezier3x3(std::vector<glm::vec3>& controlPoints);
		uint32_t addBezier4x4(std::vector<glm::vec3>& controlPoints);
		uint32_t addPlane(std::vector<glm::vec3>& controlPoints);

		uint32_t bufferSize() 
		{
			return m_data.size() * m_numSamplesU * m_numSamplesV * 3 * sizeof(glm::vec4);
		}

		std::vector<glm::vec4> getData();
		
	private:
		uint32_t m_numSamplesU;
		uint32_t m_numSamplesV;

		std::vector<LocalSurfaceData> m_data;

		glm::vec3 bezBasis3(float t);
		glm::vec3 bezBasisDer3(float t);
		glm::vec4 bezBasis4(float t);
		glm::vec4 bezBasisDer4(float t);

		inline glm::vec3 mix(glm::vec3 a, glm::vec3 b, float t) {
			return (1.0f - t) * a + (t * b);
		}
	};

}