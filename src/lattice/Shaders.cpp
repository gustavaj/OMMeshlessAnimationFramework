#include "Shaders.h"

#include <iostream>
#include <cassert>

namespace OML {

	std::unordered_map<std::string, std::vector<uint32_t>> Shaders::SpirvMap;
	std::unordered_map<std::string, std::string> Shaders::ShaderSources;
	/*shaderc::Compiler Shaders::Compiler;
	shaderc::CompileOptions Shaders::Options;*/

	// Shader names
	std::string Shaders::PosColorPassVertName = "PosColorPassVert";
	std::string Shaders::PosColorPassFragName = "PosColorPassFrag";

	std::vector<uint32_t>& Shaders::GetPosColorPassVert()
	{
		if (Shaders::SpirvMap.find(Shaders::PosColorPassVertName) == Shaders::SpirvMap.end())
		{
			std::string src =
				"#version 450\n"
				"\n"
				"layout(location = 0) in vec3 inPos;\n"
				"layout(location = 1) in uvec3 inCol;\n"
				"\n"
				"layout(location = 0) out uvec3 vCol;\n"
				"\n"
				+ Shaders::CommonUniformBufferString +
				"\n"
				"void main() {\n"
				"  gl_Position = ubo.projection * ubo.modelview * vec4(inPos, 1.0f);\n"
				"  gl_PointSize = 10.0f;\n"
				"  vCol = inCol;\n"
				"}\n";
			Shaders::LoadSpirv(Shaders::PosColorPassVertName, src, shaderc_shader_kind::shaderc_vertex_shader);
		}

		return Shaders::SpirvMap[Shaders::PosColorPassVertName];
	}

	std::vector<uint32_t>& Shaders::GetPosColorPassFrag()
	{
		if (Shaders::SpirvMap.find(Shaders::PosColorPassFragName) == Shaders::SpirvMap.end())
		{
			std::string src =
				"#version 450\n"
				"\n"
				"layout(location = 0) flat in uvec3 vCol;\n"
				"\n"
				"layout(location = 0) out vec4 FragColor;\n"
				"\n"
				"void main() {\n"
				"  if(vCol.r == 255 && vCol.g == 255 && vCol.b == 255) {\n"
				"    discard;\n"
				"  }\n"
				"  FragColor = vec4(vCol / 255.0f, 1.0f);\n"
				"}\n";
			Shaders::LoadSpirv(Shaders::PosColorPassFragName, src, shaderc_shader_kind::shaderc_fragment_shader);
		}

		return Shaders::SpirvMap[Shaders::PosColorPassFragName];
	}

	void Shaders::LoadSpirv(std::string name, std::string source, shaderc_shader_kind type)
	{
		shaderc::Compiler compiler;
		shaderc::CompileOptions options;
		options.SetSourceLanguage(shaderc_source_language_glsl);
		shaderc::SpvCompilationResult res = compiler.CompileGlslToSpv(source, type, name.c_str(), options);
		if (res.GetCompilationStatus() != shaderc_compilation_status_success)
		{
			std::cout << "Error when compiling shader " << name << "\n"
				<< res.GetNumWarnings() << " Warnings and "
				<< res.GetNumErrors() << " Errors found\nError message:\n"
				<< res.GetErrorMessage() << std::endl;
			assert(res.GetCompilationStatus() != shaderc_compilation_status_success);
		}

		std::vector<uint32_t> spirv;
		spirv.assign(res.cbegin(), res.cend());

		Shaders::SpirvMap.insert({ name, spirv });
		Shaders::ShaderSources.insert({ name, source });
	}




	// Common shader code
	std::string Shaders::CommonUniformBufferString =
		"layout(set = 0, binding = 0) uniform UniformBufferObject {\n"
		"  int tessInner;\n"
		"  int tessOuter;\n"
		"  int bFunctionIndex;\n"
		"  mat4 projection;\n"
		"  mat4 modelview;\n"
		"  mat4 normal;\n"
		"  vec2 windowSize;\n"
		"} ubo;\n";

}