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
			return m_randomDist(m_mt) * std::max(0.0, (b - a)) + a;
		}

		float random(float a, float b)
		{
			return m_randomDist(m_mt) * std::max(0.0f, (b - a)) + a;
		}
	};


	enum class SimulatorTypes {
		NormalSin = 0, RandomSphere, Rotation, XYScale
	};

	const std::vector<std::string> SimulatorNames = {
		"NormalSin", "RandomSphere", "NormalRotation", "XYScale"
	};

	class Simulator {
	public:
		Simulator()
			: m_t(0.0), m_speed(RNG.random(Simulator::MinSpeed, Simulator::MaxSpeed)) {}
		~Simulator() {}

		/* Performs simulation on the given matrix */
		virtual void simulate(double dt, glm::mat4& matrix) = 0;
		/* Undoes the transformation caused by this simulator on the given matrix */
		virtual void undoTransformation(glm::mat4& matrix) = 0;
		
		// Static member variables used to construct random simulators
		static glm::vec2 SpeedRange;
		static float MinSpeed, MaxSpeed;
		static glm::vec2 AmpRange;
		static float MinAmp, MaxAmp;
		static glm::vec2 AngleRange;
		static float MinAngle, MaxAngle;
		static glm::vec2 ScaleRange;
		static float MinScale, MaxScale;

	protected:
		double m_t;
		float m_speed;

		static Random RNG;
	};

	class NormalSinSimulator : public Simulator {
	public:
		NormalSinSimulator() : NormalSinSimulator(glm::vec3(0.0f, 0.0f, 1.0f)) {}
		NormalSinSimulator(glm::vec3 normal) : 
			Simulator(), m_amp(RNG.random(Simulator::MinAmp, Simulator::MaxAmp)), 
			m_normal(normal[0], normal[1], normal[2]), m_lastOffset(glm::vec3(0.0f)) {}
		~NormalSinSimulator() {}

		virtual void simulate(double dt, glm::mat4& matrix) override;
		virtual void undoTransformation(glm::mat4& matrix) override;

	private:
		double m_amp;
		glm::vec3 m_normal;
		glm::vec3 m_lastOffset;
	};

	class RandomSphereSimulator : public Simulator {
	public:
		RandomSphereSimulator() : 
			Simulator(), m_amp(RNG.random(Simulator::MinAmp, Simulator::MaxAmp)), 
			m_direction(glm::normalize(glm::vec3(RNG.random(0.0f, 1.0f), 
				RNG.random(0.0f, 1.0f), RNG.random(0.0f, 1.0f)))),
			m_lastOffset(glm::vec3(0.0f)) {}
		~RandomSphereSimulator() {}

		virtual void simulate(double dt, glm::mat4& matrix) override;
		virtual void undoTransformation(glm::mat4& matrix) override;

	private:
		double m_amp;
		glm::vec3 m_direction;
		glm::vec3 m_lastOffset;
	};

	class RangeRotationSimulator : public Simulator {
	public:
		RangeRotationSimulator() : RangeRotationSimulator(glm::vec3(0.0f, 0.0f, 1.0f)) {}
		RangeRotationSimulator(glm::vec3 rotAxis = glm::vec3(0.0f, 0.0f, 1.0f)) : 
			Simulator(), m_maxAngle(RNG.random(Simulator::MinAngle, Simulator::MaxAngle)), 
			m_rotAxis(rotAxis), m_lastAngle(0.0) {}
		~RangeRotationSimulator() {}

		virtual void simulate(double dt, glm::mat4& matrix) override;
		virtual void undoTransformation(glm::mat4& matrix) override;

	private:
		double m_maxAngle;
		glm::vec3 m_rotAxis;
		float m_lastAngle;
	};

	class XYScalingSimulator : public Simulator {
	public:
		XYScalingSimulator() : 
			Simulator(), m_maxScale(RNG.random(Simulator::MinScale, Simulator::MaxScale)), 
			m_lastScale(1.0f) {}
		~XYScalingSimulator() {}

		virtual void simulate(double dt, glm::mat4& matrix) override;
		virtual void undoTransformation(glm::mat4& matrix) override;

	private:
		double m_maxScale;
		float m_lastScale;
	};
}