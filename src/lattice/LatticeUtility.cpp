#include "LatticeUtility.h"

#include <iostream>
#include <string>

namespace OML {

	std::unordered_map<std::string, std::chrono::steady_clock::time_point> Timer::TimePoints;

	void Timer::Start(std::string name, std::string out)
	{
		Timer::TimePoints.insert({ name, std::chrono::high_resolution_clock::now() });
		std::cout << out << " start" << std::endl;
	}

	void Timer::Stop(std::string name, std::string out)
	{
		auto it = Timer::TimePoints.find(name);
		if (it != Timer::TimePoints.end())
		{
			auto end = std::chrono::high_resolution_clock::now();
			auto time = std::chrono::duration<double, std::milli>(end - it->second).count();

            std::cout << out << " time: " << time << "ms" << std::endl;

			Timer::TimePoints.erase(it);
		}
	}

}
