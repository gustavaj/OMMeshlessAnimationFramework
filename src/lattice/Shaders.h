#include <shaderc/shaderc.hpp>

#include <unordered_map>
#include <vector>
#include <string>

namespace OML {

	class Shaders {
	public:
		/* Names */
		static std::string PosColorPassVertName;
		static std::string PosColorPassFragName;

		/* Shaders */
		// Vertex shaders
		static std::vector<uint32_t>& GetPosColorPassVert();

		// Tessellation Control Shaders

		// Tessellation Evaluation Shaders

		// Fragment Shaders
		static std::vector<uint32_t>& GetPosColorPassFrag();


		static void LoadSpirv(std::string name, std::string source, shaderc_shader_kind type);

		static std::unordered_map<std::string, std::vector<uint32_t>> SpirvMap;
		static std::unordered_map<std::string, std::string> ShaderSources;
		/*static shaderc::Compiler Compiler;
		static shaderc::CompileOptions Options;*/

		// Common shader code
		static std::string CommonUniformBufferString;
	};

}