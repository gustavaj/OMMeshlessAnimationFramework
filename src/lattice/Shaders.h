#pragma once

#include <shaderc/shaderc.hpp>

#include <unordered_map>
#include <vector>
#include <string>

#include "LatticeUtility.h"

namespace OML {

	using NameSpirvPair = std::pair<std::string, std::vector<uint32_t>&>;

	struct ShaderOptions {
		uint32_t numControl = 0;
		uint32_t numLocal = 0;
		uint32_t numPatches = 0;
		uint32_t numBoundaries = 0;
		uint32_t numSamplesU = 0;
		uint32_t numSamplesV = 0;

		float maxError = 1.0f;
		float normalLength = 8.0f;
	};

	enum class TeseShaderType {
		Lattice = 0,
		Local,
		Pixel_Accuracy,
		Surf_Accuracy,
		Normals
	};

	class Shaders {
	public:
		// Names
		static std::string PosColorPassVertName;
		static std::string UintFlatColorFragName;
		static std::string ShadedColorFragName;
		static std::string ShadedColorNoInterpFragName;
		static std::string LocalSurfaceInfoVertName;
		static std::string FlatColorFragName;
		static std::string FlatColorNoInterpFragName;
		static std::string LocalSurfaceInfoTescName;
		static std::string BiQuadLatticeNormalsGeomName;
		static std::string TriSizeGeomName;

		static void PrintShader(std::string name);
		static void PrintShader(std::string name, std::string& source);

		// Shaders

		// Vertex shaders
		static NameSpirvPair GetPosColorPassVertShader();
		static NameSpirvPair GetLocalSurfaceInfoVertShader();

		// Tessellation Control Shaders
		static NameSpirvPair GetLocalSurfaceInfoTescShader(
			uint32_t patchVertices, LocalSurfaceType lsType, ShaderOptions& options);

		// Tessellation Evaluation Shaders
		static std::string GetTeseShaderName(LocalSurfaceType lsType, TeseShaderType teseType,
			EvaluationMethod evalMethod, ShaderOptions& options);
		static std::string GetTeseDescriptorSets(LocalSurfaceType lsType, TeseShaderType teseType,
			EvaluationMethod evalMethod, ShaderOptions& options);
		static std::string GetTeseLocalSurfaceEvaluator(LocalSurfaceType lsType, TeseShaderType teseType,
			EvaluationMethod evalMethod, ShaderOptions& options);
		static std::string GetTeseMainBody(LocalSurfaceType lsType, TeseShaderType teseType,
			EvaluationMethod evalMethod, ShaderOptions& options);
		static NameSpirvPair GetTeseShader(LocalSurfaceType lsType, TeseShaderType teseType,
			EvaluationMethod evalMethod, ShaderOptions& options);

		// Geometry Shaders
		static NameSpirvPair GetLatticeNormalsGeomShader();
		static NameSpirvPair GetTriSizeGeomShader();

		// Fragment Shaders
		static NameSpirvPair GetUintFlatColorFragShader();
		static NameSpirvPair GetFlatColorFragShader();
		static NameSpirvPair GetFlatColorNoInterpFragShader();
		static NameSpirvPair GetShadedColorFragShader();
		static NameSpirvPair GetShadedColorNoInterpFragShader();
		static NameSpirvPair GetBiQuadLatticePixelAccuracyFragShader(LocalSurfaceType lsType,
			uint32_t numLocalsurfaceControlPoints, uint32_t numLocalSurfaces, uint32_t numPatches);


		static inline bool NotInMap(std::string name)
		{
			return Shaders::SpirvMap.find(name) == Shaders::SpirvMap.end();
		}

		static void LoadSpirv(std::string name, std::string& source, shaderc_shader_kind type);

		static std::unordered_map<std::string, std::vector<uint32_t>> SpirvMap;
		static std::unordered_map<std::string, std::string> ShaderSources;
		static std::vector<std::string> ShaderNames;
		/*static shaderc::Compiler Compiler;
		static shaderc::CompileOptions Options;*/

		static void SortShaderNames();

		// Common shader code
		static std::string ShaderHeader();
		static std::string TeseHeader();
		
		static std::string CommonUniformBufferString();
		
		static std::string LocalSurfaceInfoStruct();
		static std::string SamplerStruct();
		static std::string BoundaryInfoStruct();
		static std::string QuadUniformSampler2DArrayLayout(std::string set, uint32_t firstBinding);

		static std::string MatrixStorageBufferString(
			std::string numLocalSurfaces, std::string set, std::string binding);
		static std::string ControlPointStorageBufferString(
			std::string numControlPoints, std::string set, std::string binding);
		static std::string BoundaryBufferString(
			std::string numBoundaries, std::string set, std::string binding);
		static std::string LocalSurfaceDataBufferString(
			std::string numDataPoints, std::string set, std::string binding);

		static std::string BlendingFunctionString();

		static std::string BiQuadEvaluatorString(std::string prefix);
		static std::string BiQuadEvaluatorOnlyPosString(std::string prefix);
		static std::string BiQuadEvaluatorWithDer2();
		static std::string BiCubicEvaluatorString(std::string prefix);
		static std::string BiCubicEvaluatorOnlyPosString(std::string prefix);
		static std::string BiCubicEvaluatorWithDer2();
		static std::string PlaneEvaluatorString(std::string prefix);
		static std::string PlaneEvaluatorWithDer2();

		static std::string SampleLocalSurfaceEvaluator(std::string prefix);
		static std::string BatchedSampleLocalSurfaceEvaluator(std::string prefix);
		static std::string SampleBufferEvaluator(std::string prefix, std::string nSampU, std::string nSampV);
		static std::string SampleBufferNoInterpEvaluator(std::string prefix, std::string nSampU, std::string nSampV);

		static std::string EvaluateLocalSurfacesString(std::string prefix);
		static std::string SampleLocalSurfacesString(std::string prefix);
		static std::string SampleBatchedLocalSurfacesString(std::string prefix);

		static std::string BlendLocalSurfacesString(std::string prefix);
		static std::string GetPixelAccGlobalSurfEvaluatorString();

		static std::string EvaluateLocalSurfacesOnlyPosString(std::string prefix);
		static std::string BlendLocalSurfacesOnlyPos(std::string prefix);

		static std::string ColorByMaxErrorString();

		static std::string ClipToWindowFunctionString();
	};

}