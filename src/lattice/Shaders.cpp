#include "Shaders.h"

#include <iostream>
#include <cassert>

namespace OML {

	std::unordered_map<std::string, std::vector<uint32_t>> Shaders::SpirvMap;
	std::unordered_map<std::string, std::string> Shaders::ShaderSources;
	std::vector<std::string> Shaders::ShaderNames;

	// Shader names
	std::string Shaders::PosColorPassVertName = "PosColorPassVert";
	std::string Shaders::UintFlatColorFragName = "PosColorPassFrag";
	std::string Shaders::LocalSurfaceInfoVertName = "LocalSurfaceInfoVert";
	std::string Shaders::FlatColorFragName = "FlatColorFrag";
	std::string Shaders::LocalSurfaceInfoTescName = "LocalSurfaceInfoTesc";
	std::string Shaders::ShadedColorFragName = "ShadedColorFrag";
	std::string Shaders::BiQuadLatticeNormalsGeomName = "BiQuadLatticeNormalsGeom";
	std::string Shaders::BiQuadLatticePixelAccuracyFragName = "BiQuadLatticePixelAccuracyFrag";

	void Shaders::PrintShader(std::string name)
	{
		if (!NotInMap(name))
		{
			Shaders::PrintShader(name, Shaders::ShaderSources[name]);
		}
	}

	void Shaders::PrintShader(std::string name, std::string& source)
	{
		std::cout << "\nShader " << name << " Source:\n";
		int line = 0;
		std::cout << ++line << ": ";
		for (const char& c : source)
		{
			if (c == '\n') {
				std::cout << "\n" << ++line << ": ";
			}
			else {
				std::cout << c;
			}
		}
		std::cout << "\n" << std::endl;
	}

	NameSpirvPair Shaders::GetPosColorPassVertShader()
	{
		if (Shaders::NotInMap(Shaders::PosColorPassVertName))
		{
			std::string src =
				Shaders::ShaderHeader() +
				" \n"
				"layout ( location = 0 ) in vec3 inPos;\n"
				"layout ( location = 1 ) in uvec3 inCol;\n"
				" \n"
				"layout ( location = 0 ) out uvec3 vCol;\n"
				" \n"
				+ Shaders::CommonUniformBufferString() +
				" \n"
				"void main() {\n"
				"    gl_Position = ubo.projection * ubo.modelview * vec4 ( inPos, 1.0f );\n"
				"    gl_PointSize = 10.0f;\n"
				"    vCol = inCol;\n"
				"}";
			Shaders::LoadSpirv(Shaders::PosColorPassVertName, src, shaderc_shader_kind::shaderc_vertex_shader);
		}
		return { Shaders::PosColorPassVertName, Shaders::SpirvMap[Shaders::PosColorPassVertName] };
	}

	NameSpirvPair Shaders::GetLocalSurfaceInfoVertShader()
	{
		if (Shaders::NotInMap(Shaders::LocalSurfaceInfoVertName))
		{
			std::string src =
				Shaders::ShaderHeader() +
				" \n"
				+ Shaders::LocalSurfaceInfoStruct() +
				" \n"
				"layout ( location = 0 ) in uvec4 inLSInfo;\n"
				"layout ( location = 1 ) in vec3 inColor;\n"
				" \n"
				"layout ( location = 0 ) out LocalSurfaceInfo vLSInfo;\n"
				"layout ( location = 4 ) out vec3 vColor;\n"
				" \n"
				"void main() {\n"
				"    LocalSurfaceInfo lsInfo;\n"
				"    lsInfo.controlPointIndex = inLSInfo.x;\n"
				"    lsInfo.controlPointCount = inLSInfo.y;\n"
				"    lsInfo.matrixIndex = inLSInfo.z;\n"
				"    lsInfo.boundaryIndex = inLSInfo.w;\n"
				" \n"
				"    vLSInfo = lsInfo;\n"
				"    vColor = inColor;\n"
				"}";

			Shaders::LoadSpirv(Shaders::LocalSurfaceInfoVertName, src, shaderc_shader_kind::shaderc_vertex_shader);
		}
		return { Shaders::LocalSurfaceInfoVertName, Shaders::SpirvMap[Shaders::LocalSurfaceInfoVertName] };
	}

	NameSpirvPair Shaders::GetLocalSurfaceInfoTescShader(uint32_t patchVertices)
	{
		std::string verts = std::to_string(patchVertices);
		std::string name = Shaders::LocalSurfaceInfoTescName + verts;
		if (Shaders::NotInMap(name))
		{
			std::string src =
				Shaders::ShaderHeader() +
				" \n"
				+ Shaders::LocalSurfaceInfoStruct() +
				" \n"
				"layout ( vertices = " + verts + " ) out;\n"
				" \n"
				"layout ( location = 0 ) in LocalSurfaceInfo vLSInfo[];\n"
				"layout ( location = 4 ) in vec3 vColor[];\n"
				" \n"
				"layout ( location = 0 ) out LocalSurfaceInfo tcLSInfo[];\n"
				"layout ( location = 4 ) out vec3 tcColor[];\n"
				" \n"
				+ Shaders::CommonUniformBufferString() +
				" \n"
				"void main() {\n"
				"    tcLSInfo[ gl_InvocationID ] = vLSInfo[ gl_InvocationID ];\n"
				"    tcColor[ gl_InvocationID ] = vColor[ gl_InvocationID ];\n"
				" \n"
				"    if( gl_InvocationID == 0 ) {\n"
				"        gl_TessLevelInner [0] = ubo.tessInner;\n"
				"        gl_TessLevelInner [1] = ubo.tessInner;\n"
				"        gl_TessLevelOuter [0] = ubo.tessOuter;\n"
				"        gl_TessLevelOuter [1] = ubo.tessOuter;\n"
				"        gl_TessLevelOuter [2] = ubo.tessOuter;\n"
				"        gl_TessLevelOuter [3] = ubo.tessOuter;\n"
			    "    }\n"
				"}";
			Shaders::LoadSpirv(name, src, shaderc_shader_kind::shaderc_tess_control_shader);
		}
		return { name, Shaders::SpirvMap[name] };
	}

	std::string Shaders::GetTeseShaderName(
		LocalSurfaceType lsType, TeseShaderType teseType, EvaluationMethod evalMethod, ShaderOptions& options)
	{
		std::string name;

		if (evalMethod == EvaluationMethod::Direct)
		{
			name += "Direct";

			switch (teseType)
			{
			case TeseShaderType::Lattice:			name += "Lattice"; break;
			case TeseShaderType::Local:				name += "Local"; break;
			case TeseShaderType::Pixel_Accuracy:	name += "LatPixAcc"; break;
			case TeseShaderType::Surf_Accuracy:		name += "LatSurfAcc"; break;
			case TeseShaderType::Normals:			name += "LateNormals"; break;
			}

			switch (lsType)
			{
			case LocalSurfaceType::Quadratic_Bezier:  name += "QuadBezier"; break;
			case LocalSurfaceType::Cubic_Bezier: name += "CubicBezier"; break;
			}
		} 
		else 
		{
			switch (evalMethod)
			{
			case EvaluationMethod::Pre_Sampled_Image:			name += "PreImg"; break;
			case EvaluationMethod::Pre_Sampled_Image_Batched:	name += "PreImgBatch"; break;
			case EvaluationMethod::Pre_Sampled_Buffer:			name += "PreBuf"; break;
			}

			switch (teseType)
			{
			case TeseShaderType::Lattice:			name += "Lattice"; break;
			case TeseShaderType::Local:				name += "Local"; break;
			case TeseShaderType::Pixel_Accuracy:	name += "LatPixAcc"; break;
			case TeseShaderType::Surf_Accuracy:		name += "LatSurfAcc"; break;
			case TeseShaderType::Normals:			name += "LatNormals"; break;
			}

			name += "_" + std::to_string(options.numSamplesU) + "x" + std::to_string(options.numSamplesV) + "S";
		}

		name += "_" + std::to_string(options.numControl) + "C_" + std::to_string(options.numLocal) + "L_"
			+ std::to_string(options.numPatches) + "P_Tese";

		return name;
	}

	std::string Shaders::GetTeseDescriptorSets(
		LocalSurfaceType lsType, TeseShaderType teseType, EvaluationMethod evalMethod, ShaderOptions& options)
	{
		std::string nControl = std::to_string(options.numControl);
		std::string nLocal = std::to_string(options.numLocal);
		std::string nPatches = std::to_string(options.numPatches);
		std::string nBoundaries = std::to_string(options.numPatches * 4);
		std::string desc =
			Shaders::CommonUniformBufferString() + " \n" +
			Shaders::MatrixStorageBufferString(nLocal, "0", "1") + " \n" +
			Shaders::ControlPointStorageBufferString(nControl, "0", "2") + " \n" +
			Shaders::BoundaryBufferString(nBoundaries, "0", "3") + "\n\n";

		if (evalMethod == EvaluationMethod::Pre_Sampled_Image)
		{
			if (teseType == TeseShaderType::Local)
			{
				desc += "layout ( set = 2, binding = 0 ) uniform sampler2DArray localSampler;\n\n";
			}
			else {
				desc += Shaders::QuadUniformSampler2DArrayLayout("1", 0) + "\n\n";
			}
		}
		else if (evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
		{
			desc += "layout ( set = 1, binding = 0 ) uniform sampler2DArray surfSampler;\n\n";
		}
		else if (evalMethod == EvaluationMethod::Pre_Sampled_Buffer)
		{
			std::string nLocalSurfaceData = std::to_string(options.numLocal * options.numSamplesU * options.numSamplesV * 3);
			desc += Shaders::LocalSurfaceDataBufferString(nLocalSurfaceData, "0", "4") + "\n\n";
		}

		switch (teseType)
		{
		case TeseShaderType::Normals:
		case TeseShaderType::Lattice: {
			desc +=
				"layout ( location = 0 ) in LocalSurfaceInfo tcLSInfo[];\n"
				"layout ( location = 4 ) in vec3 tcColor[];\n"
				" \n"
				"layout ( location = 0 ) out vec3 teColor;\n"
				"layout ( location = 1 ) out vec3 teNormal;\n"
				"layout ( location = 2 ) out vec3 tePosition;\n";
			break;
		}
		case TeseShaderType::Surf_Accuracy:
		case TeseShaderType::Local: {
			desc +=
				"layout ( location = 0 ) in LocalSurfaceInfo tcLSInfo[];\n"
				"layout ( location = 4 ) in vec3 tcColor[];\n"
				" \n"
				"layout ( location = 0 ) out vec3 teColor;\n";
			break;
		}
		case TeseShaderType::Pixel_Accuracy: {
			desc +=
				"layout ( location = 0 ) in LocalSurfaceInfo tcLSInfo[];\n"
				"layout ( location = 4 ) in vec3 tcColor[];\n"
				" \n"
				"layout ( location = 0 ) out vec3 teColor;\n"
				"layout ( location = 1 ) out vec3 teNormal;\n"
				"layout ( location = 2 ) out vec3 tePosition;\n"
				"layout ( location = 3 ) out vec2 teUVCoords;\n"
				"layout ( location = 4 ) out LocalSurfaceInfo ls00Info;\n"
				"layout ( location = 8 ) out LocalSurfaceInfo ls10Info;\n"
				"layout ( location = 12 ) out LocalSurfaceInfo ls01Info;\n"
				"layout ( location = 16 ) out LocalSurfaceInfo ls11Info;\n";
			break;
		}
		}

		return desc;
	}

	std::string Shaders::GetTeseLocalSurfaceEvaluator(
		LocalSurfaceType lsType, TeseShaderType teseType, EvaluationMethod evalMethod, ShaderOptions& options)
	{
		std::string eval;

		if (teseType != TeseShaderType::Local) {
			eval +=
				Shaders::BlendingFunctionTypes() + "\n\n" +
				Shaders::BlendingFunctionString() + "\n\n";
		}

		std::string nSampU = std::to_string(options.numSamplesU);
		std::string nSampV = std::to_string(options.numSamplesV);

		if (teseType == TeseShaderType::Surf_Accuracy)
		{
			switch (evalMethod)
			{
			case EvaluationMethod::Pre_Sampled_Image: {
				eval += Shaders::BiQuadEvaluatorOnlyPosString("ref_") + " \n";
				eval += Shaders::SampleLocalSurfaceEvaluator("");
				break;
			}
			case EvaluationMethod::Pre_Sampled_Image_Batched: {
				eval += Shaders::BiQuadEvaluatorOnlyPosString("ref_") + " \n";
				eval += Shaders::BatchedSampleLocalSurfaceEvaluator("");
				break;
			}
			case EvaluationMethod::Pre_Sampled_Buffer: {
				eval += Shaders::BiQuadEvaluatorOnlyPosString("ref_") + " \n";
				eval += Shaders::SampleBufferEvaluator("", nSampU, nSampV);
				break;
			}
			}

		}
		else
		{
			switch (evalMethod)
			{
			case EvaluationMethod::Direct: eval += Shaders::BiQuadEvaluatorString(""); break;
			case EvaluationMethod::Pre_Sampled_Image: eval += Shaders::SampleLocalSurfaceEvaluator(""); break;
			case EvaluationMethod::Pre_Sampled_Image_Batched: eval += Shaders::BatchedSampleLocalSurfaceEvaluator(""); break;
			case EvaluationMethod::Pre_Sampled_Buffer: eval += Shaders::SampleBufferEvaluator("", nSampU, nSampV); break;
			}
		}

		return eval;
	}

	std::string Shaders::GetTeseMainBody(
		LocalSurfaceType lsType, TeseShaderType teseType, EvaluationMethod evalMethod, ShaderOptions& options)
	{
		std::string main =
			"\nvoid main() {\n"
			"    float u = gl_TessCoord.x;\n"
			"    float v = gl_TessCoord.y;\n\n";

		if (teseType == TeseShaderType::Local)
		{
			if (evalMethod == EvaluationMethod::Pre_Sampled_Image)
			{
				main += "    Sampler local = evaluateLocal(tcLSInfo[0], localSampler, u, v);\n";
			}
			else
			{
				main += "    Sampler local = evaluateLocal(tcLSInfo[0], u, v);\n";
			}
			main +=
				"    gl_Position = ubo.projection * ubo.modelview * vec4 ( local.p, 1.0f );\n"
				"    teColor = tcColor[0];\n"
				"}";
		}
		else if (teseType == TeseShaderType::Surf_Accuracy)
		{
			if (evalMethod == EvaluationMethod::Pre_Sampled_Image)
			{
				main += Shaders::SampleLocalSurfacesString("") + " \n";
			}
			else if (evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
			{
				main += Shaders::SampleBatchedLocalSurfacesString("") + " \n";
			}
			else
			{
				main += Shaders::EvaluateLocalSurfacesOnlyPosString("") + " \n";
			}

			main +=
				Shaders::BlendLocalSurfacesOnlyPos("") + " \n" +
				"    gl_Position = ubo.projection * ubo.modelview * vec4( pos, 1.0f );\n" +
				Shaders::EvaluateLocalSurfacesOnlyPosString("ref_") + " \n" +
				Shaders::BlendLocalSurfacesOnlyPos("ref_") + " \n" +
				"    float diff = distance ( pos, ref_pos );\n"
				"    float maxError = " + std::to_string(options.maxError) + ";\n" +
				Shaders::ColorByMaxErrorString() + " \n" +
				"}";
		}
		else
		{
			if (evalMethod == EvaluationMethod::Pre_Sampled_Image)
			{
				main += Shaders::SampleLocalSurfacesString("") + " \n";
			}
			else if (evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
			{
				main += Shaders::SampleBatchedLocalSurfacesString("") + " \n";
			}
			else
			{
				main += Shaders::EvaluateLocalSurfacesString("") + " \n";
			}

			main +=
				Shaders::BlendLocalSurfacesString("") + " \n" +
				"    gl_Position = ubo.projection * ubo.modelview * vec4 (pos, 1.0f);\n"
				"    tePosition = vec3 (ubo.modelview * vec4 (pos, 1.0f));\n";
			if (teseType == TeseShaderType::Normals)
			{
				main += "    teNormal = normalize ( cross (dpdu, dpdv) );\n";
			}
			else
			{
				main += "    teNormal = vec3 (ubo.normal * vec4 ( normalize ( cross (dpdu, dpdv) ), 0.0f) );\n";
			}
			main += "    teColor = tcColor[0];\n";

			if (teseType == TeseShaderType::Pixel_Accuracy)
			{
				main +=
					"    teUVCoords = vec2 (u, v);\n"
					"    ls00Info = tcLSInfo[0];\n"
					"    ls10Info = tcLSInfo[1];\n"
					"    ls01Info = tcLSInfo[2];\n"
					"    ls11Info = tcLSInfo[3];\n";
			}

			main += "}";
		}

		return main;
	}

	NameSpirvPair Shaders::GetTeseShader(
		LocalSurfaceType lsType, TeseShaderType teseType, EvaluationMethod evalMethod, ShaderOptions& options)
	{
		std::string name = Shaders::GetTeseShaderName(lsType, teseType, evalMethod, options);
		if (Shaders::NotInMap(name))
		{
			std::string src =
				Shaders::ShaderHeader() + " \n" +
				Shaders::TeseHeader() + " \n" +
				Shaders::LocalSurfaceInfoStruct() + " \n" +
				Shaders::SamplerStruct() + " \n" +
				Shaders::BoundaryInfoStruct() + " \n" +
				Shaders::GetTeseDescriptorSets(lsType, teseType, evalMethod, options) + " \n" +
				Shaders::GetTeseLocalSurfaceEvaluator(lsType, teseType, evalMethod, options) + " \n" +
				Shaders::GetTeseMainBody(lsType, teseType, evalMethod, options);

			Shaders::LoadSpirv(name, src, shaderc_shader_kind::shaderc_tess_evaluation_shader);
		}
		return { name, Shaders::SpirvMap[name] };
	}

	NameSpirvPair Shaders::GetLatticeNormalsGeomShader(float length)
	{
		if (Shaders::NotInMap(Shaders::BiQuadLatticeNormalsGeomName))
		{
			std::string l = std::to_string(length);
			std::string src =
				Shaders::ShaderHeader() + " \n" +
				"layout ( triangles ) in;\n"
				" \n"
				"layout ( location = 1 ) in vec3 teNormal[];\n"
				"layout ( location = 2 ) in vec3 tePosition[];\n"
				" \n" +
				Shaders::CommonUniformBufferString() + " \n" +
				" \n"
				"layout ( line_strip , max_vertices = 18 ) out;\n"
				" \n"
				"layout ( location = 0 ) out vec3 gColor;\n"
				" \n"
				"void main() {\n"
				"    for ( int vertex = 0; vertex < 3; ++vertex ) {\n"
				"        gl_Position = ubo.projection * vec4 (tePosition[vertex], 1);\n"
				"        gColor = vec3 (0.2);\n"
				"        EmitVertex ();\n"
				" \n"
				"        gl_Position = ubo.projection * ( vec4 (tePosition[vertex], 1)\n"
				"            + ( transpose ( inverse (ubo.modelview)) * vec4 (teNormal[vertex] * " + l + ", 0)));\n"
				"        gColor = vec3 (0.6);\n"
				"        EmitVertex ();\n"
				" \n"
				"        EndPrimitive ();\n"
				"    }\n"
				"}";
			Shaders::LoadSpirv(Shaders::BiQuadLatticeNormalsGeomName, src, shaderc_shader_kind::shaderc_geometry_shader);
		}
		return { Shaders::BiQuadLatticeNormalsGeomName, Shaders::SpirvMap[Shaders::BiQuadLatticeNormalsGeomName] };
	}

	NameSpirvPair Shaders::GetUintFlatColorFragShader()
	{
		if (Shaders::NotInMap(Shaders::UintFlatColorFragName))
		{
			std::string src =
				Shaders::ShaderHeader() +
				" \n"
				"layout ( location = 0) flat in uvec3 vCol;\n"
				" \n"
				"layout ( location = 0) out vec4 FragColor;\n"
				" \n"
				"void main() {\n"
				"    if (vCol.r == 255 && vCol.g == 255 && vCol.b == 255) {\n"
				"        discard;\n"
				"    }\n"
				"    FragColor = vec4 (vCol / 255.0f, 1.0f);\n"
				"}";
			Shaders::LoadSpirv(Shaders::UintFlatColorFragName, src, shaderc_shader_kind::shaderc_fragment_shader);
		}
		return { Shaders::UintFlatColorFragName, Shaders::SpirvMap[Shaders::UintFlatColorFragName] };
	}

	NameSpirvPair Shaders::GetFlatColorFragShader()
	{
		if (Shaders::NotInMap(Shaders::FlatColorFragName))
		{
			std::string src =
				Shaders::ShaderHeader() + 
				" \n"
				"layout ( location = 0 ) in vec3 teColor;\n"
				" \n"
				"layout ( location = 0 ) out vec4 FragColor;\n"
				" \n"
				"void main() {\n"
				"    FragColor = vec4 (teColor, 1.0);\n"
				"}";
			Shaders::LoadSpirv(Shaders::FlatColorFragName, src, shaderc_shader_kind::shaderc_fragment_shader);
		}
		return { Shaders::FlatColorFragName, Shaders::SpirvMap[Shaders::FlatColorFragName] };
	}

	NameSpirvPair Shaders::GetShadedColorFragShader()
	{
		if (Shaders::NotInMap(Shaders::ShadedColorFragName))
		{
			std::string src = 
				Shaders::ShaderHeader() +
				" \n"
				"layout ( location = 0 ) in vec3 teColor;\n"
				"layout ( location = 1 ) in vec3 teNormal;\n"
				" \n"
				"layout ( location = 0 ) out vec4 FragColor;\n"
				" \n"
				"void main() {\n"
				"    vec3 light = normalize ( vec3 (0.5f, 0.5f, 1.0f));\n"
				"    vec3 ambient = vec3 (0.3f);\n"
				"    vec3 diffuse = max ( dot (light, teNormal), 0.0f) * vec3 (0.8f);\n"
				" \n"
				"    FragColor = vec4 ((ambient + diffuse) * teColor, 1.0);\n"
				"}";
			Shaders::LoadSpirv(Shaders::ShadedColorFragName, src, shaderc_shader_kind::shaderc_fragment_shader);
		}
		return { Shaders::ShadedColorFragName, Shaders::SpirvMap[Shaders::ShadedColorFragName] };
	}

	NameSpirvPair Shaders::GetBiQuadLatticePixelAccuracyFragShader(
		uint32_t numLocalsurfaceControlPoints, uint32_t numLocalSurfaces, uint32_t numPatches)
	{
		std::string nControl = std::to_string(numLocalsurfaceControlPoints);
		std::string nLocal = std::to_string(numLocalSurfaces);
		std::string nPatches = std::to_string(numPatches);
		std::string name = Shaders::BiQuadLatticePixelAccuracyFragName
			+ "_" + nControl + "C_" + nLocal + "L_" + nPatches + "P";
		if (Shaders::NotInMap(name))
		{
			std::string nBoundaries = std::to_string(numPatches * 4);
			std::string src =
				Shaders::ShaderHeader() + " \n" +
				Shaders::LocalSurfaceInfoStruct() + " \n" +
				Shaders::SamplerStruct() + " \n" +
				Shaders::BoundaryInfoStruct() + " \n" +
				Shaders::BlendingFunctionTypes() + " \n" +
				" \n"
				"layout ( location = 0 ) in vec3 teColor;\n"
				"layout ( location = 1 ) in vec3 teNormal;\n"
				"layout ( location = 2 ) in vec3 tePosition;\n"
				"layout ( location = 3 ) in vec2 teUVCoords;\n"
				"layout ( location = 4 ) flat in LocalSurfaceInfo ls00Info;\n"
				"layout ( location = 8 ) flat in LocalSurfaceInfo ls10Info;\n"
				"layout ( location = 12 ) flat in LocalSurfaceInfo ls01Info;\n"
				"layout ( location = 16 ) flat in LocalSurfaceInfo ls11Info;\n"
				" \n"
				"layout ( location = 0 ) out vec4 FragColor;\n"
				" \n" +
				Shaders::CommonUniformBufferString() + " \n" +
				Shaders::MatrixStorageBufferString(nLocal, "0", "1") + " \n" +
				Shaders::ControlPointStorageBufferString(nControl, "0", "2") + " \n" +
				Shaders::BoundaryBufferString(nBoundaries, "0", "3") + " \n" +
				Shaders::BlendingFunctionString() + " \n" +
				Shaders::BiQuadEvaluatorOnlyPosString("") + " \n" +
				"void main() {\n"
				"    float u = teUVCoords[0];\n"
				"    float v = teUVCoords[1];\n"
				" \n" +
				"    Sampler s00 = evaluateLocal(ls00Info, u, v);\n"
				"    Sampler s10 = evaluateLocal(ls10Info, u, v);\n"
				"    Sampler s01 = evaluateLocal(ls01Info, u, v);\n"
				"    Sampler s11 = evaluateLocal(ls11Info, u, v);\n"
				" \n"
				"    float Bu = 1.0 - bFunction(u);\n"
				"    float Bv = 1.0 - bFunction(v);\n"
				" \n"
				"    vec3 Su0 = s10.p + (s00.p - s10.p) * Bu;\n"
				"    vec3 Su1 = s11.p + (s01.p - s11.p) * Bu;\n"
				"    vec3 p = Su1 + (Su0 - Su1) * Bv;\n"
				" \n"
				"    vec4 pos = ubo.projection * ubo.modelview * vec4 (p, 1.0f);\n"
				"    vec4 ndc = pos / pos.w;\n"
				"    float windowX = (ubo.windowSize.x / 2) * (ndc.x + 1);\n"
				"    float windowY = (ubo.windowSize.y / 2) * (ndc.y + 1);\n"
				"    float windowZ = 0.5 * ndc.z + 0.5;\n"
				"    vec4 windowPos = vec4 (windowX, windowY, windowZ, 1.0f);\n"
				" \n"
				"    float error = distance (gl_FragCoord.xy, windowPos.xy);\n"
				"    vec3 color;\n"
				" \n"
				"    if (error <= 0.5) {\n"
				"        color = vec3 (0.7, 0.7, 0.7);\n"
				"    }\n"
				"    else if (error <= 1.0) {\n"
				"        color = vec3 (0.0, 1.0, 0.0);\n"
				"    }\n"
				"    else if (error <= 2.0) {\n"
				"        color = vec3 (0.0, 0.0, 1.0);\n"
				"    }\n"
				"    else if (error <= 5.0) {\n"
				"        color = vec3 (1.0, 1.0, 0.0);\n"
				"    }\n"
				"    else {\n"
				"        color = vec3 (1.0, 0.0, 0.0);\n"
				"    }\n"
				" \n"
				"    FragColor = vec4 (color, 1.0);\n"
				"}";
			Shaders::LoadSpirv(name, src, shaderc_shader_kind::shaderc_fragment_shader);
		}
		return { name, Shaders::SpirvMap[name] };
	}

	void Shaders::LoadSpirv(std::string name, std::string& source, shaderc_shader_kind type)
	{
		shaderc::Compiler compiler;
		//shaderc::CompileOptions options;
		// This optimization thing adds a lot of cpu cost, like 10x, should try to test if it gives any benefit.
		//options.SetOptimizationLevel(shaderc_optimization_level::shaderc_optimization_level_performance);
		shaderc::SpvCompilationResult res = compiler.CompileGlslToSpv(source, type, name.c_str()/*, options*/);
		if (res.GetCompilationStatus() != shaderc_compilation_status_success)
		{
			std::cout << "Error when compiling shader " << name << "\n"
				<< res.GetNumWarnings() << " Warnings and "
				<< res.GetNumErrors() << " Errors found\nError message:\n"
				<< res.GetErrorMessage() << std::endl;
			Shaders::PrintShader(name, source);
			assert(false);
		}

		std::vector<uint32_t> spirv;
		spirv.assign(res.cbegin(), res.cend());

		Shaders::SpirvMap.insert({ name, spirv });
		Shaders::ShaderSources.insert({ name, source });
		Shaders::ShaderNames.push_back(name);
	}




	// Common shader code
	inline std::string Shaders::ShaderHeader()
	{ 
		return "#version 450\n";
	}

	inline std::string Shaders::TeseHeader()
	{
		return "layout ( quads, equal_spacing, ccw ) in;\n";
	}

	inline std::string Shaders::CommonUniformBufferString()
	{
		return
			"layout ( set = 0, binding = 0 ) uniform UniformBufferObject {\n"
			"    int tessInner;\n"
			"    int tessOuter;\n"
			"    int bFunctionIndex;\n"
			"    mat4 projection;\n"
			"    mat4 modelview;\n"
			"    mat4 normal;\n"
			"    vec2 windowSize;\n"
			"} ubo;\n";
	}

	inline std::string Shaders::LocalSurfaceInfoStruct() 
	{
		return 
			"struct LocalSurfaceInfo {\n"
			"    uint controlPointIndex;\n"
			"    uint controlPointCount;\n"
			"    uint matrixIndex;\n"
			"    uint boundaryIndex;\n"
			"};\n";
	}

	inline std::string Shaders::SamplerStruct() 
	{
		return 
			"struct Sampler {\n"
			"    vec3 p, u, v;\n"
			"};\n";
	}

	inline std::string Shaders::BoundaryInfoStruct() 
	{
		return
			"struct BoundaryInfo {\n"
			"    float us, ue, vs, ve;\n"
			"};\n";
	}

	inline std::string Shaders::QuadUniformSampler2DArrayLayout(std::string set, uint32_t firstBinding)
	{
		return
			"layout ( set = " + set + ", binding = " + std::to_string(firstBinding) + " ) uniform sampler2DArray p00Sampler;\n"
			"layout ( set = " + set + ", binding = " + std::to_string(firstBinding+1) + " ) uniform sampler2DArray p10Sampler;\n"
			"layout ( set = " + set + ", binding = " + std::to_string(firstBinding+2) + " ) uniform sampler2DArray p01Sampler;\n"
			"layout ( set = " + set + ", binding = " + std::to_string(firstBinding+3) + " ) uniform sampler2DArray p11Sampler;\n";
	}

	inline std::string Shaders::MatrixStorageBufferString(
		std::string numLocalSurfaces, std::string set, std::string binding) 
	{
		return
			"layout ( set = " + set + ", binding = 1 ) buffer MatrixBuffer {\n"
			"    mat4 matrices[" + numLocalSurfaces + "];\n"
			"} matrixBuffer;\n";
	}

	inline std::string Shaders::ControlPointStorageBufferString(
		std::string numControlPoints, std::string set, std::string binding)
	{
		return
			"layout ( set = " + set + ", binding = " + binding + " ) buffer ControlPointBuffer {\n"
			"    vec4 controlPoints[" + numControlPoints + "];\n"
			"} controlPointBuffer;\n";
	}

	inline std::string Shaders::BoundaryBufferString(
		std::string numBoundaries, std::string set, std::string binding)
	{
		return
			"layout ( set = " + set + ", binding = " + binding + " ) buffer BoundaryBuffer {\n"
			"    BoundaryInfo boundaries[" + numBoundaries + "];\n"
			"} boundaryBuffer;\n";
	}

	inline std::string Shaders::LocalSurfaceDataBufferString(
		std::string numDataPoints, std::string set, std::string binding)
	{
		return
			"layout ( set = " + set + ", binding = " + binding + " ) buffer LocalSurfaceDataBuffer {\n"
			"    vec4 localSurfaceData[" + numDataPoints + "];\n"
			"} lsDataBuffer;\n";
	}

	inline std::string Shaders::BlendingFunctionString()
	{
		return
			"float bFunction( float w ) {\n"
			"    if (w < 1e-5) return 0.0f;\n"
			"    if (1 - w < 1e-5) return 1.0f;\n"
			"    switch (ubo.bFunctionIndex) {\n"
			"        case B1Poly: { return (3 * pow (w, 2) - 2 * pow (w, 3)); }\n"
			"        case B2Poly: { return 6 * pow (w, 5) - 15 * pow (w, 4) + 10 * pow (w, 3); }\n"
			"        case Lerbs: { return 1.0 / (1.0 + exp (1.0 / w - 1.0 / (1.0 - w))); }\n"
			"    }\n"
			"}\n"
			" \n"
			"float bDerivative( float w ) {\n"
			"    if (w < 1e-5) return 0.0f;\n"
			"    if (1 - w < 1e-5) return 0.0f;\n"
			"    switch (ubo.bFunctionIndex) {\n"
			"    case B1Poly: { return 6 * w - 6 * pow (w, 2); }\n"
			"    case B2Poly: { return 30 * pow (w, 4) - 60 * pow (w, 3) + 30 * pow (w, 2); }\n"
			"    case Lerbs: {\n"
			"        return (2 * exp (-((1.0) / ((w - 1.0) * w))) *\n"
			"            (pow (w, 2.0) - w + 0.5)) /\n"
			"            (pow (exp (-(1.0 / (w - 1.0))) + exp (1.0 / w), 2.0) *\n"
			"             pow (w - 1.0, 2.0) * pow (w, 2.0)); }\n"
			"    }\n"
			"}\n";
	}

	inline std::string Shaders::BlendingFunctionTypes()
	{
		return
			"const int B1Poly = 0;\n"
			"const int B2Poly = 1;\n"
			"const int Lerbs = 2;\n";
	}

	inline std::string Shaders::BiQuadEvaluatorString(std::string prefix)
	{
		return
			"vec3 bez3Basis( float t ) {\n"
			"    return vec3 ( pow(1-t, 2), 2*t*(1-t), pow (t, 2));\n"
			"}\n"
			" \n"
			"vec3 bez3BasisDer( float t ) {\n"
			"    return vec3 (2*t-2, 2-4*t, 2*t);\n"
			"}\n"
			" \n"
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"    BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"    float local_u = mix (bi.us, bi.ue, u);\n"
			"    float local_v = mix (bi.vs, bi.ve, v);\n"
			" \n"
			"    vec3 bu = bez3Basis(local_u);\n"
			"    vec3 bud = bez3BasisDer(local_u);\n"
			"    vec3 bv = bez3Basis(local_v);\n"
			"    vec3 bvd = bez3BasisDer(local_v);\n"
			" \n"
			"    vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"    vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"    vec3 p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"    vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			"    vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz;\n"
			"    vec3 p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;\n"
			"    vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz;\n"
			"    vec3 p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz;\n"
			"    vec3 p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;\n"
			" \n"
			"    vec3 pos =\n"
			"        p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] +\n"
			"        p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] +\n"
			"        p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2];\n"
			" \n"
			"    vec3 dpdu =\n"
			"        p00 * bud[0] * bv[0] + p01 * bud[0] * bv[1] + p02 * bud[0] * bv[2] +\n"
			"        p10 * bud[1] * bv[0] + p11 * bud[1] * bv[1] + p12 * bud[1] * bv[2] +\n"
			"        p20 * bud[2] * bv[0] + p21 * bud[2] * bv[1] + p22 * bud[2] * bv[2];\n"
			" \n"
			"    vec3 dpdv =\n"
			"		p00 * bu[0] * bvd[0] + p01 * bu[0] * bvd[1] + p02 * bu[0] * bvd[2] +\n"
			"		p10 * bu[1] * bvd[0] + p11 * bu[1] * bvd[1] + p12 * bu[1] * bvd[2] +\n"
			"		p20 * bu[2] * bvd[0] + p21 * bu[2] * bvd[1] + p22 * bu[2] * bvd[2];\n"
			" \n"
			"    mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"	 pos = vec3 (matrix * vec4 (pos, 1.0f));\n"
			"    // mat4 normalMat = transpose ( inverse (matrix));\n"
			"    dpdu = vec3 (matrix * vec4 (dpdu, 0.0f));\n"
			"    dpdv = vec3 (matrix * vec4 (dpdv, 0.0f));\n"
			" \n"
			"    return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);\n"
			"}\n";
	}

	inline std::string Shaders::BiQuadEvaluatorOnlyPosString(std::string prefix)
	{
		return
			"vec3 bez3Basis( float t ) {\n"
			"    return vec3 ( pow (1-t, 2), 2*t*(1-t), pow (t, 2));\n"
			"}\n"
			" \n"
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"    BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"    float local_u = mix (bi.us, bi.ue, u);\n"
			"    float local_v = mix (bi.vs, bi.ve, v);\n"
			" \n"
			"    vec3 bu = bez3Basis(local_u);\n"
			"    vec3 bv = bez3Basis(local_v);\n"
			" \n"
			"    vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"    vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"    vec3 p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"    vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			"    vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz;\n"
			"    vec3 p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;\n"
			"    vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz;\n"
			"    vec3 p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz;\n"
			"    vec3 p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;\n"
			" \n"
			"    vec3 pos =\n"
			"		p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] +\n"
			"		p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] +\n"
			"		p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2];\n"
			" \n"
			"    mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"    vec4 p = matrix * vec4 (pos, 1.0f);\n"
			" \n"
			"    return Sampler(p.xyz, vec3(0.0), vec3(0.0));\n"
			"}\n";
	}

	inline std::string Shaders::SampleLocalSurfaceEvaluator(std::string prefix)
	{
		return
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, sampler2DArray surfSampler, float u, float v) {\n"
			"    BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			" \n"
			"    vec2 coords = vec2 ( mix (bi.us, bi.ue, u), mix (bi.vs, bi.ve, v));\n"
			" \n"
			"    vec3 pos = texture (surfSampler, vec3 (coords, 0)).xyz;\n"
			"    vec3 dpdu = texture (surfSampler, vec3 (coords, 1)).xyz;\n"
			"    vec3 dpdv = texture (surfSampler, vec3 (coords, 2)).xyz;\n"
			" \n"
			"    mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"    pos = vec3 (matrix * vec4 (pos, 1.0f));\n"
			"    // mat4 normalMat = transpose ( inverse (matrix));\n"
			"    dpdu = vec3 (matrix * vec4 (dpdu, 0.0f));\n"
			"    dpdv = vec3 (matrix * vec4 (dpdv, 0.0f));\n"
			" \n"
			"    return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);\n"
			"}\n";
	}

	inline std::string Shaders::BatchedSampleLocalSurfaceEvaluator(std::string prefix)
	{
		return
			"Sampler evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v, int layer) {\n"
			"    ivec3 dim = textureSize (surfSampler, 0);\n"
			"    float ss = float (lsInfo.controlPointIndex + 0.5);\n"
			"    float se = float (lsInfo.controlPointIndex + lsInfo.boundaryIndex - 0.5);\n"
			"    float ts = float (lsInfo.controlPointCount + 0.5);\n"
			"    float te = float (lsInfo.controlPointCount + lsInfo.boundaryIndex - 0.5);\n"
			"    float s = mix (ss, se, u) / float (dim.x);\n"
			"    float t = mix (ts, te, v) / float (dim.y);\n"
			" \n"
			"    vec3 pos = texture (surfSampler, vec3 (s, t, layer)).xyz;\n"
			"    vec3 dpdu = texture (surfSampler, vec3 (s, t, layer + 1)).xyz;\n"
			"    vec3 dpdv = texture (surfSampler, vec3 (s, t, layer + 2)).xyz;\n"
			" \n"
			"    mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"    pos = vec3 (matrix * vec4 (pos, 1.0f));\n"
			"    // mat4 normalMat = transpose ( inverse (matrix));\n"
			"    dpdu = vec3 (matrix * vec4 (dpdu, 0.0f));\n"
			"    dpdv = vec3 (matrix * vec4 (dpdv, 0.0f));\n"
			" \n"
			"    return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);\n"
			"}\n";
	}

	inline std::string Shaders::SampleBufferEvaluator(std::string prefix, std::string nSampU, std::string nSampV)
	{
		return
			"// Takes in a u and v coordinate and returns a Sampler object (p,u,v)\n"
			"// Uses linear interpolation to mix the neighbouring values when the u and\n"
			"// v coordinates are between sample values\n"
			"Sampler sampleBuffer( int idx, float u, float v ) {\n"
			"    Sampler res;\n"
			"    // Step one, transform u,v coordinates [0,1]^2 to buffer coordinates [0,numSamplesU]x[0,numSamplesV]\n"
			"    float fi = u * (" + nSampU + " - 1);\n"
			"    float fj = v * (" + nSampV + " - 1);\n"
			"    int i = int ( floor(fi));\n"
			"    int j = int ( floor(fj));\n"
			"    float fu = fract (fi);\n"
			"    float fv = fract (fj);\n"
			" \n"
			"    int c00 = j * " + nSampU + " + i;\n"
			"    int c10 = j * " + nSampU + " + i + 1;\n"
			"    int c01 = (j + 1) * " + nSampU + " + i;\n"
			"    int c11 = (j + 1) * " + nSampU + " + i + 1;\n"
			" \n"
			"    // Step two, linearly interpolate between neighbouring values\n"
			"    vec3 p00 = lsDataBuffer.localSurfaceData[idx + c00].xyz;\n"
			"    vec3 p10 = lsDataBuffer.localSurfaceData[idx + c10].xyz;\n"
			"    vec3 p01 = lsDataBuffer.localSurfaceData[idx + c01].xyz;\n"
			"    vec3 p11 = lsDataBuffer.localSurfaceData[idx + c11].xyz;\n"
			"    res.p = (((1 - fu) * p00 + fu * p10) * (1 - fv) + ((1 - fu) * p01 + fu * p11) * fv);\n"
			" \n"
			"    int offsetU = " + nSampU + " * " + nSampV + ";\n"
			"    vec3 u00 = lsDataBuffer.localSurfaceData[idx + c00 + offsetU].xyz;\n"
			"    vec3 u10 = lsDataBuffer.localSurfaceData[idx + c10 + offsetU].xyz;\n"
			"    vec3 u01 = lsDataBuffer.localSurfaceData[idx + c01 + offsetU].xyz;\n"
			"    vec3 u11 = lsDataBuffer.localSurfaceData[idx + c11 + offsetU].xyz;\n"
			"    res.u = (((1 - fu) * u00 + fu * u10) * (1 - fv) + ((1 - fu) * u01 + fu * u11) * fv);\n"
			" \n"
			"    int offsetV = offsetU * 2;\n"
			"    vec3 v00 = lsDataBuffer.localSurfaceData[idx + c00 + offsetV].xyz;\n"
			"    vec3 v10 = lsDataBuffer.localSurfaceData[idx + c10 + offsetV].xyz;\n"
			"    vec3 v01 = lsDataBuffer.localSurfaceData[idx + c01 + offsetV].xyz;\n"
			"    vec3 v11 = lsDataBuffer.localSurfaceData[idx + c11 + offsetV].xyz;\n"
			"    res.v = (((1 - fu) * v00 + fu * v10) * (1 - fv) + ((1 - fu) * v01 + fu * v11) * fv);\n"
			" \n"
			"    // Return the sampled values\n"
			"    return res;\n"
			"}\n"
			" \n"
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"    BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			" \n"
			"    Sampler s = sampleBuffer(int(lsInfo.controlPointIndex * " + nSampU + " * " + nSampV + " * 3), \n"
			"        mix (bi.us, bi.ue, u), mix(bi.vs, bi.ve, v));\n"
			" \n"
			"    mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"    s.p = vec3 (matrix * vec4 (s.p, 1.0f));\n"
			"    // Should probbly do this! Or make another matrix buffer for normal matrices\n"
			"    // mat4 normalMat = transpose(inverse(matrix));\n"
			"    s.u = vec3 (matrix * vec4 (s.u, 0.0f));\n"
			"    s.v = vec3 (matrix * vec4 (s.v, 0.0f));\n"
			" \n"
			"    return s;\n"
			"}\n";
	}

	inline std::string Shaders::EvaluateLocalSurfacesString(std::string prefix)
	{
		return
			"    // Evaluate Local Surfaces\n"
			"    Sampler " + prefix + "s00 = " + prefix + "evaluateLocal(tcLSInfo[0], u, v);\n"
			"    Sampler " + prefix + "s10 = " + prefix + "evaluateLocal(tcLSInfo[1], u, v);\n"
			"    Sampler " + prefix + "s01 = " + prefix + "evaluateLocal(tcLSInfo[2], u, v);\n"
			"    Sampler " + prefix + "s11 = " + prefix + "evaluateLocal(tcLSInfo[3], u, v);\n";
	}

	inline std::string Shaders::SampleLocalSurfacesString(std::string prefix)
	{
		return
			"    // Evaluate Local Surfaces\n"
			"    Sampler " + prefix + "s00 = " + prefix + "evaluateLocal(tcLSInfo[0], p00Sampler, u, v);\n"
			"    Sampler " + prefix + "s10 = " + prefix + "evaluateLocal(tcLSInfo[1], p10Sampler, u, v);\n"
			"    Sampler " + prefix + "s01 = " + prefix + "evaluateLocal(tcLSInfo[2], p01Sampler, u, v);\n"
			"    Sampler " + prefix + "s11 = " + prefix + "evaluateLocal(tcLSInfo[3], p11Sampler, u, v);\n";
	}

	std::string Shaders::SampleBatchedLocalSurfacesString(std::string prefix)
	{
		return
			"    // Evaluate Local Surfaces\n"
			"    Sampler " + prefix + "s00 = " + prefix + "evaluateLocal(tcLSInfo[0], u, v, 0);\n"
			"    Sampler " + prefix + "s10 = " + prefix + "evaluateLocal(tcLSInfo[1], u, v, 3);\n"
			"    Sampler " + prefix + "s01 = " + prefix + "evaluateLocal(tcLSInfo[2], u, v, 6);\n"
			"    Sampler " + prefix + "s11 = " + prefix + "evaluateLocal(tcLSInfo[3], u, v, 9);\n";
	}

	inline std::string Shaders::BlendLocalSurfacesString(std::string prefix)
	{
		return
			"    // Evaluate tensor product blending spline surface\n"
			"    float " + prefix + "Bu = 1.0 - bFunction(u);\n"
			"    float " + prefix + "Bv = 1.0 - bFunction(v);\n"
			"    float " + prefix + "Bu1 = -bDerivative(u);\n"
			"    float " + prefix + "Bv1 = -bDerivative(v);\n"
			" \n"
			"    vec3 " + prefix + "Su0 = " + prefix + "s10.p + (" + prefix + "s00.p - " + prefix + "s10.p) * " + prefix + "Bu;\n"
			"    vec3 " + prefix + "Su1 = " + prefix + "s11.p + (" + prefix + "s01.p - " + prefix + "s11.p) * " + prefix + "Bu;\n"
			" \n"
			"    vec3 " + prefix + "Uu0 = " + prefix + "s10.u + (" + prefix + "s00.u - " + prefix + "s10.u) * " + prefix + "Bu\n"
			"        + (" + prefix + "s00.p - " + prefix + "s10.p) * " + prefix + "Bu1;\n"
			"    vec3 " + prefix + "Uu1 = " + prefix + "s11.u + (" + prefix + "s01.u - " + prefix + "s11.u) * " + prefix + "Bu\n"
			"        + (" + prefix + "s01.p - " + prefix + "s11.p) * " + prefix + "Bu1;\n"
			" \n"
			"    vec3 " + prefix + "Vu0 = " + prefix + "s10.v + (" + prefix + "s00.v - " + prefix + "s10.v) * " + prefix + "Bu;\n"
			"    vec3 " + prefix + "Vu1 = " + prefix + "s11.v + (" + prefix + "s01.v - " + prefix + "s11.v) * " + prefix + "Bu;\n"
			" \n"
			"    vec3 " + prefix + "pos = " + prefix + "Su1 + (" + prefix + "Su0 - " + prefix + "Su1) * " + prefix + "Bv;\n"
			"    vec3 " + prefix + "dpdu = " + prefix + "Uu1 + (" + prefix + "Uu0 - " + prefix + "Uu1) * " + prefix + "Bv;\n"
			"    vec3 " + prefix + "dpdv = " + prefix + "Vu1 + (" + prefix + "Vu0 - " + prefix + "Vu1) * " + prefix + "Bv\n"
			"        + (" + prefix + "Su0 - " + prefix + "Su1) * " + prefix + "Bv1;\n";
	}

	inline std::string Shaders::EvaluateLocalSurfacesOnlyPosString(std::string prefix)
	{
		return
			"    // Evaluate Local Surfaces\n"
			"    Sampler " + prefix + "s00 = " + prefix + "evaluateLocal(tcLSInfo[0], u, v);\n"
			"    Sampler " + prefix + "s10 = " + prefix + "evaluateLocal(tcLSInfo[1], u, v);\n"
			"    Sampler " + prefix + "s01 = " + prefix + "evaluateLocal(tcLSInfo[2], u, v);\n"
			"    Sampler " + prefix + "s11 = " + prefix + "evaluateLocal(tcLSInfo[3], u, v);\n";
	}

	inline std::string Shaders::BlendLocalSurfacesOnlyPos(std::string prefix)
	{
		return
			"    // Evaluate tensor product blending spline surface\n"
			"    float " + prefix + "Bu = 1.0 - bFunction(u);\n"
			"    float " + prefix + "Bv = 1.0 - bFunction(v);\n"
			" \n"
			"    vec3 " + prefix + "Su0 = " + prefix + "s10.p + (" + prefix + "s00.p - " + prefix + "s10.p) * " + prefix + "Bu;\n"
			"    vec3 " + prefix + "Su1 = " + prefix + "s11.p + (" + prefix + "s01.p - " + prefix + "s11.p) * " + prefix + "Bu;\n"
			" \n"
			"    vec3 " + prefix + "pos = " + prefix + "Su1 + (" + prefix + "Su0 - " + prefix + "Su1) * " + prefix + "Bv;\n";
	}

	inline std::string Shaders::ColorByMaxErrorString()
	{
		return
			"    if(diff <= maxError * 0.1) {\n"
			"        teColor = vec3 (0.5, 0.5, 0.5);\n"
			"    }\n"
			"    else if (diff <= maxError * 0.2) {\n"
			"        teColor = vec3 (0.0, 1.0, 0.0);\n"
			"    }\n"
			"    else if (diff <= maxError * 0.5) {\n"
			"        teColor = vec3 (0.0, 0.0, 1.0);\n"
			"    }\n"
			"    else if (diff <= maxError) {\n"
			"        teColor = vec3 (1.0, 1.0, 0.0);\n"
			"    }\n"
			"    else {\n"
			"        teColor = vec3 (1.0, 0.0, 0.0);\n"
			"    }\n";
	}

}