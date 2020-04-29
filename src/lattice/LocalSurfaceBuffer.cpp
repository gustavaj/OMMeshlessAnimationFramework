#include "LocalSurfaceBuffer.h"

namespace OML {

	uint32_t LocalSurfaceBuffer::addBezier3x3(std::vector<glm::vec3> controlPoints)
	{
		int numSamples = m_numSamplesU * m_numSamplesV;

		LocalSurfaceData lsData;
		lsData.positions.resize(numSamples);
		lsData.partialU.resize(numSamples);
		lsData.partialV.resize(numSamples);

		glm::vec3& p00 = controlPoints[0], p10 = controlPoints[1], p20 = controlPoints[2];
		glm::vec3& p01 = controlPoints[3], p11 = controlPoints[4], p21 = controlPoints[5];
		glm::vec3& p02 = controlPoints[6], p12 = controlPoints[7], p22 = controlPoints[8];

		float du = 1.0f / (float)(m_numSamplesU - 1);
		float dv = 1.0f / (float)(m_numSamplesV - 1);

		for (size_t j = 0; j < m_numSamplesV; j++)
		{
			float v = dv * j;
			glm::vec3 bv = bezBasis(v);
			glm::vec3 bvd = bezBasisDer(v);
			for (size_t i = 0; i < m_numSamplesU; i++)
			{
				float u = du * i;
				glm::vec3 bu = bezBasis(u);
				glm::vec3 bud = bezBasisDer(u);

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

				lsData.positions[j * m_numSamplesV + i] = glm::vec4(pos.x, pos.y, pos.z, 1.0f);
				lsData.partialU[j * m_numSamplesV + i] = glm::vec4(dpdu.x, dpdu.y, dpdu.z, 0.0f);
				lsData.partialV[j * m_numSamplesV + i] = glm::vec4(dpdv.x, dpdv.y, dpdv.z, 0.0f);
			}
		}

		m_data.push_back(std::move(lsData));

		return m_data.size() - 1;
	}

	std::vector<glm::vec4> LocalSurfaceBuffer::getData()
	{
		std::vector<glm::vec4> data;

		for (auto& localSurfaceData : m_data)
		{
			data.insert(data.end(), localSurfaceData.positions.begin(), localSurfaceData.positions.end());
			data.insert(data.end(), localSurfaceData.partialU.begin(), localSurfaceData.partialU.end());
			data.insert(data.end(), localSurfaceData.partialV.begin(), localSurfaceData.partialV.end());
		}

		return data;
	}
	
	glm::vec3 LocalSurfaceBuffer::bezBasis(float t)
	{
		return glm::vec3(std::pow(1 - t, 2), 2 * t * (1 - t), std::pow(t, 2));
	}

	glm::vec3 LocalSurfaceBuffer::bezBasisDer(float t)
	{
		return glm::vec3(2 * t - 2, 2 - 4 * t, 2 * t);
	}

}