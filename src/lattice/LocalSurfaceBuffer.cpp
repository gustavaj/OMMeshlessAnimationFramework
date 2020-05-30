#include "LocalSurfaceBuffer.h"

namespace OML {
	uint32_t LocalSurfaceBuffer::addLocalSurface(std::vector<glm::vec3>& controlPoints, LocalSurfaceType lsType, int numSamplesU, int numSamplesV)
	{
		m_numSamplesU = numSamplesU;
		m_numSamplesV = numSamplesV;
		switch (lsType)
		{
		case LocalSurfaceType::Quadratic_Bezier: return addBezier3x3(controlPoints);
		case LocalSurfaceType::Cubic_Bezier: return addBezier4x4(controlPoints);
		case LocalSurfaceType::Plane: return addPlane(controlPoints);
		}
	}

	uint32_t LocalSurfaceBuffer::addBezier3x3(std::vector<glm::vec3>& controlPoints)
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
			glm::vec3 bv = bezBasis3(v);
			glm::vec3 bvd = bezBasisDer3(v);
			for (size_t i = 0; i < m_numSamplesU; i++)
			{
				float u = du * i;
				glm::vec3 bu = bezBasis3(u);
				glm::vec3 bud = bezBasisDer3(u);

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

	uint32_t LocalSurfaceBuffer::addBezier4x4(std::vector<glm::vec3>& controlPoints)
	{
		int numSamples = m_numSamplesU * m_numSamplesV;

		LocalSurfaceData lsData;
		lsData.positions.resize(numSamples);
		lsData.partialU.resize(numSamples);
		lsData.partialV.resize(numSamples);

		glm::vec3& p00 = controlPoints[0], p10 = controlPoints[1], p20 = controlPoints[2], p30 = controlPoints[3];
		glm::vec3& p01 = controlPoints[4], p11 = controlPoints[5], p21 = controlPoints[6], p31 = controlPoints[7];
		glm::vec3& p02 = controlPoints[8], p12 = controlPoints[9], p22 = controlPoints[10], p32 = controlPoints[11];
		glm::vec3& p03 = controlPoints[12], p13 = controlPoints[13], p23 = controlPoints[14], p33 = controlPoints[15];

		float du = 1.0f / (float)(m_numSamplesU - 1);
		float dv = 1.0f / (float)(m_numSamplesV - 1);

		for (size_t j = 0; j < m_numSamplesV; j++)
		{
			float v = dv * j;
			glm::vec4 bv = bezBasis4(v);
			glm::vec4 bvd = bezBasisDer4(v);
			for (size_t i = 0; i < m_numSamplesU; i++)
			{
				float u = du * i;
				glm::vec4 bu = bezBasis4(u);
				glm::vec4 bud = bezBasisDer4(u);

				glm::vec3 pos =
					p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] + p03 * bu[0] * bv[3] +
					p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] + p13 * bu[1] * bv[3] +
					p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2] + p23 * bu[2] * bv[3] +
					p30 * bu[3] * bv[0] + p31 * bu[3] * bv[1] + p32 * bu[3] * bv[2] + p33 * bu[3] * bv[3];

				glm::vec3 dpdu =
					p00 * bud[0] * bv[0] + p01 * bud[0] * bv[1] + p02 * bud[0] * bv[2] + p03 * bud[0] * bv[3] +
					p10 * bud[1] * bv[0] + p11 * bud[1] * bv[1] + p12 * bud[1] * bv[2] + p13 * bud[1] * bv[3] +
					p20 * bud[2] * bv[0] + p21 * bud[2] * bv[1] + p22 * bud[2] * bv[2] + p23 * bud[2] * bv[3] +
					p30 * bud[3] * bv[0] + p31 * bud[3] * bv[1] + p32 * bud[3] * bv[2] + p33 * bud[3] * bv[3];

				glm::vec3 dpdv =
					p00 * bu[0] * bvd[0] + p01 * bu[0] * bvd[1] + p02 * bu[0] * bvd[2] + p03 * bu[0] * bvd[3] +
					p10 * bu[1] * bvd[0] + p11 * bu[1] * bvd[1] + p12 * bu[1] * bvd[2] + p13 * bu[1] * bvd[3] +
					p20 * bu[2] * bvd[0] + p21 * bu[2] * bvd[1] + p22 * bu[2] * bvd[2] + p23 * bu[2] * bvd[3] +
					p30 * bu[3] * bvd[0] + p31 * bu[3] * bvd[1] + p32 * bu[3] * bvd[2] + p33 * bu[3] * bvd[3];

				lsData.positions[j * m_numSamplesU + i] = glm::vec4(pos.x, pos.y, pos.z, 1.0f);
				lsData.partialU[j * m_numSamplesU + i] = glm::vec4(dpdu.x, dpdu.y, dpdu.z, 0.0f);
				lsData.partialV[j * m_numSamplesU + i] = glm::vec4(dpdv.x, dpdv.y, dpdv.z, 0.0f);
			}
		}

		m_data.push_back(std::move(lsData));

		return m_data.size() - 1;
	}

	uint32_t LocalSurfaceBuffer::addPlane(std::vector<glm::vec3>& controlPoints)
	{
		int numSamples = m_numSamplesU * m_numSamplesV;

		LocalSurfaceData lsData;
		lsData.positions.resize(numSamples);
		lsData.partialU.resize(numSamples);
		lsData.partialV.resize(numSamples);

		glm::vec3& p00 = controlPoints[0], p10 = controlPoints[1];
		glm::vec3& p01 = controlPoints[2], p11 = controlPoints[3];

		float du = 1.0f / (float)(m_numSamplesU - 1);
		float dv = 1.0f / (float)(m_numSamplesV - 1);

		for (size_t j = 0; j < m_numSamplesV; j++)
		{
			float v = dv * j;
			for (size_t i = 0; i < m_numSamplesU; i++)
			{
				float u = du * i;

				glm::vec3 pos = mix(mix(p00, p10, u), mix(p01, p11, u), v);
				glm::vec3 dpdu = p10 - p00;
				glm::vec3 dpdv = p01 - p00;

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

	glm::vec3 LocalSurfaceBuffer::bezBasis3(float t)
	{
		return glm::vec3(std::pow(1 - t, 2), 2 * t * (1 - t), std::pow(t, 2));
	}

	glm::vec3 LocalSurfaceBuffer::bezBasisDer3(float t)
	{
		return glm::vec3(2 * t - 2, 2 - 4 * t, 2 * t);
	}

	glm::vec4 LocalSurfaceBuffer::bezBasis4(float t)
	{
		return glm::vec4(std::pow(1 - t, 3), 3 * t * std::pow(1 - t, 2), 3 * t * t * (1 - t), std::pow(t, 3));
	}

	glm::vec4 LocalSurfaceBuffer::bezBasisDer4(float t)
	{
		return glm::vec4(-3 * std::pow(1 - t, 2), 3 * (1 - t) * (1 - 3 * t), 3 * t * (2 - 3 * t), 3 * t * t);
	}

}