#pragma once

#include <glm/glm.hpp>

#include <random>

namespace OML {

	struct Random {
		Random() {
			std::random_device rd;
			m_mt = std::mt19937(rd());
			m_randomDist = std::uniform_real_distribution<double>(0.0, 1.0);
		}

		std::mt19937 m_mt;
		std::uniform_real_distribution<double> m_randomDist;

		double random(double a, double b)
		{
			return m_randomDist(m_mt) * std::max(1.0, (b - a)) + a;
		}

		float random(float a, float b)
		{
			return m_randomDist(m_mt) * std::max(1.0f, (b - a)) + a;
		}
	};


	enum class SimulatorTypes {
		None = 0, NormalSin, RandomSphere
	};

	const std::vector<std::string> SimulatorNames = {
		"None", "NormalSin", "RandomSphere"
	};

	class Simulator {
	public:
		Simulator() : Simulator(0.0, 1.0) {}
		Simulator(double t, double speed)
			: m_t(t), m_speed(speed), m_lastOffset(glm::vec3(0.0f)) {}
		~Simulator() {}

		virtual void simulate(double dt, glm::mat4& matrix) = 0;
	protected:
		double m_t;
		double m_speed;
		glm::vec3 m_lastOffset;
	};

	class NormalSinSimulator : public Simulator {
	public:
		NormalSinSimulator() : NormalSinSimulator(0.0, 0.0, 0.0, glm::vec3(0.0f, 0.0f, 0.0f)) {}
		NormalSinSimulator(double t, double speed, double max, glm::vec3 normal)
			: Simulator(t, speed), m_max(max), m_normal(normal[0], normal[1], normal[2]) {}
		~NormalSinSimulator() {}

		virtual void simulate(double dt, glm::mat4& matrix) override {
			m_t += dt;
			float f = std::sin(m_t * m_speed) * m_max;
			glm::vec3 offset = (m_normal * f);
			glm::vec4 trans(offset[0] - m_lastOffset[0], offset[1] - m_lastOffset[1], offset[2] - m_lastOffset[2], 0.0f);
			matrix[3] -= trans;
			m_lastOffset = offset;
		}

	private:
		double m_max;
		glm::vec3 m_normal;
	};

	class RandomSphereSimulator : public Simulator {
	public:
		RandomSphereSimulator() : RandomSphereSimulator(0.0, 1.0, 1.0f, glm::vec3(1.0f)) {}
		RandomSphereSimulator(double t, double speed, double max, glm::vec3 dir)
			: Simulator(t, speed), m_max(max), m_direction(glm::normalize(dir)) {}
		~RandomSphereSimulator() {}

		virtual void simulate(double dt, glm::mat4& matrix) override {
			m_t += dt;

			glm::vec3 offset = m_direction * (float)(dt * m_speed);
			glm::vec4 trans(offset[0], offset[1], offset[2], 0.0f);

			matrix[3] -= trans;
			m_lastOffset += offset;
			if (m_lastOffset.length() > m_max) {
				m_direction = -m_direction;// glm::normalize((-m_direction +
					//glm::vec3(m_rng.random(0.0, 1.0), m_rng.random(0.0, 1.0), m_rng.random(-1.0, 1.0))));
			}
		}

	private:
		double m_max;
		glm::vec3 m_direction;

		Random m_rng;
	};
}