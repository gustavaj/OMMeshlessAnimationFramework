#include "Simulators.h"

#include <glm/gtc/matrix_transform.hpp>

namespace OML {

	// Parameters used for constructing simulators.
	Random Simulator::RNG{};

	float Simulator::MinSpeed = 0.0f;
	float Simulator::MaxSpeed = 5.0f;
	glm::vec2 Simulator::SpeedRange = glm::vec2(0.0f, 20.0f);

	float Simulator::MinAmp = 0.0f;
	float Simulator::MaxAmp = 5.0f;
	glm::vec2 Simulator::AmpRange = glm::vec2(0.0f, 30.0f);

	float Simulator::MinAngle = 10.0f;
	float Simulator::MaxAngle = 60.0f;
	glm::vec2 Simulator::AngleRange = glm::vec2(0.0f, 90.0f);

	float Simulator::MinScale = 0.1f;
	float Simulator::MaxScale = 0.5f;
	glm::vec2 Simulator::ScaleRange = glm::vec2(0.0f, 1.0f);

	void NormalSinSimulator::simulate(double dt, glm::mat4& matrix)
	{
		m_t += dt;
		float f = std::sin(m_t * m_speed) * m_amp;
		glm::vec3 offset = (m_normal * f);
		glm::vec4 trans(offset[0] - m_lastOffset[0], offset[1] - m_lastOffset[1], offset[2] - m_lastOffset[2], 0.0f);
		matrix[3] += trans;
		m_lastOffset = offset;
	}

	void NormalSinSimulator::undoTransformation(glm::mat4& matrix)
	{
		glm::vec4 trans(-m_lastOffset[0], -m_lastOffset[1], -m_lastOffset[2], 0.0f);
		matrix[3] += trans;
		m_lastOffset = glm::vec3(0.0f);
	}




	void RandomSphereSimulator::simulate(double dt, glm::mat4& matrix)
	{
		m_t += dt;

		glm::vec3 offset = m_direction * (float)(dt * m_speed);
		glm::vec4 trans(offset[0], offset[1], offset[2], 0.0f);

		matrix[3] += trans;
		m_lastOffset += offset;
		if (m_lastOffset.length() > m_amp) {
			m_direction = -m_direction;// glm::normalize((-m_direction +
				//glm::vec3(RNG.random(0.0, 1.0), RNG.random(0.0, 1.0), RNG.random(-1.0, 1.0))));
		}
	}

	void RandomSphereSimulator::undoTransformation(glm::mat4& matrix)
	{
		glm::vec4 trans(-m_lastOffset[0], -m_lastOffset[1], -m_lastOffset[2], 0.0f);
		matrix[3] += trans;
		m_lastOffset = glm::vec3(0.0f);
	}




	void RangeRotationSimulator::simulate(double dt, glm::mat4& matrix)
	{
		m_t += dt;

		float angle = std::sin(m_t * m_speed) * m_maxAngle;
		matrix = glm::rotate(matrix, glm::radians(angle - m_lastAngle), m_rotAxis);
		m_lastAngle = angle;
	}

	void RangeRotationSimulator::undoTransformation(glm::mat4& matrix)
	{
		matrix = glm::rotate(matrix, glm::radians(-m_lastAngle), m_rotAxis);
		m_lastAngle = 0.0f;
	}




	void XYScalingSimulator::simulate(double dt, glm::mat4& matrix)
	{
		m_t += dt;

		float newScale = 1.0f + std::sin(m_t * m_speed) * m_maxScale;
		float s = newScale / m_lastScale;
		glm::vec3 scale = glm::vec3(s, s, 1.0f);
		matrix = glm::scale(matrix, scale);
		m_lastScale = newScale;
	}

	void XYScalingSimulator::undoTransformation(glm::mat4& matrix)
	{
		float s = 1.0f / m_lastScale;
		glm::vec3 scale = glm::vec3(s, s, 1.0f);
		matrix = glm::scale(matrix, scale);
		m_lastScale = 1.0f;
	}

}