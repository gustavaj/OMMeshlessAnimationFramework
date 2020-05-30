#include "Shaders.h"

#include <iostream>
#include <cassert>
#include <algorithm>

namespace OML {

	std::unordered_map<std::string, std::vector<uint32_t>> Shaders::SpirvMap;
	std::unordered_map<std::string, std::string> Shaders::ShaderSources;
	std::vector<std::string> Shaders::ShaderNames;

	// Shader names
	std::string Shaders::PosColorPassVertName = "Vert_PosColorPass";
	std::string Shaders::UintFlatColorFragName = "Frag_PosColorPass";
	std::string Shaders::LocalSurfaceInfoVertName = "Vert_LocalSurfaceInfo";
	std::string Shaders::FlatColorFragName = "Frag_FlatColor";
	std::string Shaders::FlatColorNoInterpFragName = "Frag_FlatColorNoInterp";
	std::string Shaders::LocalSurfaceInfoTescName = "LocalSurfaceInfo";
	std::string Shaders::ShadedColorFragName = "Frag_ShadedColor";
	std::string Shaders::ShadedColorNoInterpFragName = "Frag_ShadedColorNoInterp";
	std::string Shaders::BiQuadLatticeNormalsGeomName = "Geom_LatticeNormals";
	std::string Shaders::TriSizeGeomName = "Geom_TriSize";

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
				"layout(location = 0) in vec3 inPos;\n"
				"layout(location = 1) in uvec3 inCol;\n"
				" \n"
				"layout(location = 0) out uvec3 vCol;\n"
				" \n"
				+ Shaders::CommonUniformBufferString() +
				" \n"
				"void main() {\n"
				"  gl_Position = ubo.projection * ubo.modelview * vec4(inPos, 1.0f);\n"
				"  gl_PointSize = 10.0f;\n"
				"  vCol = inCol;\n"
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
				"layout(location = 0) in uvec4 inLSInfo;\n"
				"layout(location = 1) in vec3 inColor;\n"
				"layout(location = 2) in uvec4 inLSInfo2;\n"
				" \n"
				"layout(location = 0) out LocalSurfaceInfo vLSInfo;\n"
				"layout(location = 5) out vec3 vColor;\n"
				" \n"
				"void main() {\n"
				"  LocalSurfaceInfo lsInfo;\n"
				"  lsInfo.controlPointIndex = inLSInfo.x;\n"
				"  lsInfo.matrixIndex = inLSInfo.z;\n"
				"  lsInfo.boundaryIndex = inLSInfo.w;\n"
				"  lsInfo.batchInfo = uvec3(inLSInfo2.x, inLSInfo2.y, inLSInfo2.z);\n"
				"  lsInfo.dataIndex = inLSInfo2.w;\n"
				" \n"
				"  vLSInfo = lsInfo;\n"
				"  vColor = inColor;\n"
				"}";

			Shaders::LoadSpirv(Shaders::LocalSurfaceInfoVertName, src, shaderc_shader_kind::shaderc_vertex_shader);
		}
		return { Shaders::LocalSurfaceInfoVertName, Shaders::SpirvMap[Shaders::LocalSurfaceInfoVertName] };
	}

	NameSpirvPair Shaders::GetLocalSurfaceInfoTescShader(uint32_t patchVertices, LocalSurfaceType lsType, ShaderOptions& options)
	{
		std::string verts = std::to_string(patchVertices);
		std::string nControl = std::to_string(options.numControl);
		std::string nLocal = std::to_string(options.numLocal);
		std::string nPatches = std::to_string(options.numPatches);
		std::string nBoundaries = std::to_string(options.numPatches * 4);
		std::string evalString;
		std::string name = "Tesc_";
		/*switch (lsType)
		{
		case LocalSurfaceType::Quadratic_Bezier: {
			name += "QuadBezier_";
			evalString = Shaders::BiQuadEvaluatorWithDer2();
			break;
		}
		case LocalSurfaceType::Cubic_Bezier: {
			name += "CubicBezier_";
			evalString = Shaders::BiCubicEvaluatorWithDer2();
			break;
		}
		case LocalSurfaceType::Plane: {
			name += "Plane_";
			evalString = Shaders::PlaneEvaluatorWithDer2();
			break;
		}
		}*/
		name += "CubicBezier_";
		evalString = Shaders::BiCubicEvaluatorWithDer2();
		name += Shaders::LocalSurfaceInfoTescName + "_" + verts + "V_" + nControl + "C_" + nLocal + "L_" + nPatches + "P";
		if (Shaders::NotInMap(name))
		{
			std::string src =
				Shaders::ShaderHeader() +
				" \n" +
				Shaders::LocalSurfaceInfoStruct() + " \n" +
				Shaders::SamplerStruct() + " \n" +
				Shaders::BoundaryInfoStruct() + " \n" +
				" \n"
				"layout(vertices = " + verts + ") out;\n"
				" \n"
				"layout(location = 0) in LocalSurfaceInfo vLSInfo[];\n"
				"layout(location = 5) in vec3 vColor[];\n"
				" \n"
				"layout(location = 0) out LocalSurfaceInfo tcLSInfo[];\n"
				"layout(location = 5) out vec3 tcColor[];\n"
				" \n"
				"const int TessFactorMethod_Static = 0;\n"
				"const int TessFactorMethod_Dynamic = 1;\n"
				"const int TessFactorMethod_PixelAccurate = 2;\n"
				" \n" +
				Shaders::CommonUniformBufferString() + " \n" +
				Shaders::MatrixStorageBufferString(nLocal, "0", "1") + " \n" +
				Shaders::ControlPointStorageBufferString(nControl, "0", "2") + " \n" +
				Shaders::BoundaryBufferString(nBoundaries, "0", "3") + " \n" +
				Shaders::BlendingFunctionString() + " \n" +
				Shaders::GetPixelAccGlobalSurfEvaluatorString() + " \n" +
				evalString + " \n" +
				" \n"
				"float getPostProjectionSphereExtent(vec3 origin, float diameter) {\n"
				"  vec4 clipPos = ubo.projection * ubo.modelview * vec4(origin, 1.0f);\n"
				"  return abs(diameter * ubo.projection[1][1] / clipPos.w);\n"
				"}\n"
				" \n"
				"float calculateTessellationFactor(vec3 control0, vec3 control1) {\n"
				"  float e0 = distance(control0, control1);\n"
				"  vec3 m0 = (control0 + control1) / 2;\n"
				"  return max(1.0f, (ubo.windowSize.y / ubo.pixelsPerEdge) * getPostProjectionSphereExtent(m0, e0));\n"
				"}\n"
				" \n"
				"void main() {\n"
				"  tcLSInfo[gl_InvocationID] = vLSInfo[gl_InvocationID];\n"
				"  tcColor[gl_InvocationID] = vColor[gl_InvocationID];\n"
				" \n"
				"  if(gl_InvocationID == 0) {\n";
			if (patchVertices == 1)
			{
				src +=
					"    gl_TessLevelInner[0] = ubo.tessInner;\n"
					"    gl_TessLevelInner[1] = ubo.tessInner;\n"
					"    gl_TessLevelOuter[0] = ubo.tessOuter;\n"
					"    gl_TessLevelOuter[1] = ubo.tessOuter;\n"
					"    gl_TessLevelOuter[2] = ubo.tessOuter;\n"
					"    gl_TessLevelOuter[3] = ubo.tessOuter;\n"
					"  }\n"
					"}\n";
			}
			else 
			{
				src +=
					"    if(ubo.tessFactorMethod == TessFactorMethod_Static) {\n"
					"      gl_TessLevelInner[0] = ubo.tessInner;\n"
					"      gl_TessLevelInner[1] = ubo.tessInner;\n"
					"      gl_TessLevelOuter[0] = ubo.tessOuter;\n"
					"      gl_TessLevelOuter[1] = ubo.tessOuter;\n"
					"      gl_TessLevelOuter[2] = ubo.tessOuter;\n"
					"      gl_TessLevelOuter[3] = ubo.tessOuter;\n"
					"    }\n"
					"    else if(ubo.tessFactorMethod == TessFactorMethod_Dynamic) {\n"
					"      // https://developer.nvidia.com/content/dynamic-hardware-tessellation-basics \n"
					"      // Calculate the four corner points\n"
					"      vec3 topLeft  = evaluateLocal(tcLSInfo[0], 0.0, 0.0).p;\n"
					"      vec3 topRight = evaluateLocal(tcLSInfo[1], 1.0, 0.0).p;\n"
					"      vec3 botLeft  = evaluateLocal(tcLSInfo[2], 0.0, 1.0).p;\n"
					"      vec3 botRight = evaluateLocal(tcLSInfo[3], 1.0, 1.0).p;\n"
					" \n"
					"      // Use corner points to calculate outer edge tess factors\n"
					"      float tf0 = calculateTessellationFactor(topLeft, botLeft);\n"
					"      float tf1 = calculateTessellationFactor(topLeft, topRight);\n"
					"      float tf2 = calculateTessellationFactor(topRight, botRight);\n"
					"      float tf3 = calculateTessellationFactor(botLeft, botRight);\n"
					"      gl_TessLevelOuter[0] = tf0;\n"
					"      gl_TessLevelOuter[1] = tf1;\n"
					"      gl_TessLevelOuter[2] = tf2;\n"
					"      gl_TessLevelOuter[3] = tf3;\n"
					" \n"
					"      // Set inner tess factors based on the max factor in the respective directions\n"
					"      // Should probably be no problems with cracks since the edges are the same\n"
					"      gl_TessLevelInner[0] = max(tf1, tf3);\n"
					"      gl_TessLevelInner[1] = max(tf0, tf2);\n"
					"    }\n"
					"    else if(ubo.tessFactorMethod == TessFactorMethod_PixelAccurate) {\n"
					"      // Find the bounds on the second derivative and min/max values for surface\n"
					"      // By evaluating the surface ubo.numPASamples^2 times\n"
					"      float Zuu_max = -1.0f;\n"
					"      float Zuv_max = -1.0f;\n"
					"      float Zvv_max = -1.0f;\n"
					"      float Xuu_max = -1.0f;\n"
					"      float Xuv_max = -1.0f;\n"
					"      float Xvv_max = -1.0f;\n"
					"      float X_max = -1.0f;\n"
					"      float Z_min = 100000000.0f;\n"
					" \n"
					"      for(int i = 0; i < ubo.numPASamples; i++) {\n"
					"        float u = float(i) / float(ubo.numPASamples - 1);\n"
					"        for(int j = 0; j < ubo.numPASamples; j++) {\n"
					"          float v = float(j) / float(ubo.numPASamples - 1);\n"
					"          Sampler2 res = evaluateGlobal(\n"
					"            evaluateLocal(tcLSInfo[0], u, v),\n"
					"            evaluateLocal(tcLSInfo[1], u, v),\n"
					"            evaluateLocal(tcLSInfo[2], u, v),\n"
					"            evaluateLocal(tcLSInfo[3], u, v), u, v);\n"
					"          res.uu = vec3(ubo.modelview * vec4(res.uu, 0.0f));\n"
					"          res.uv = vec3(ubo.modelview * vec4(res.uv, 0.0f));\n"
					"          res.vv = vec3(ubo.modelview * vec4(res.vv, 0.0f));\n"
					"          res.p = vec3(ubo.modelview * vec4(res.p, 1.0f));\n"
					"          Zuu_max = max(Zuu_max, abs(res.uu.z));\n"
					"          Zuv_max = max(Zuv_max, abs(res.uv.z));\n"
					"          Zvv_max = max(Zvv_max, abs(res.vv.z));\n"
					"          Xuu_max = max(Xuu_max, abs(res.uu.x));\n"
					"          Xuv_max = max(Xuv_max, abs(res.uv.x));\n"
					"          Xvv_max = max(Xvv_max, abs(res.vv.x));\n"
					"          X_max = max(X_max, abs(res.p.x));\n"
					"          Z_min = min(Z_min, abs(res.p.z));\n"
					"        }\n"
					"      }\n"
					" \n"
					"      // Find the tess factor using Hjelmervik's method\n"
					"      float eps_x = 0.5f / ubo.windowSize.x;\n"
					"      float A = ubo.projection[0][0];\n"
					"      float target_error = (8 * eps_x) / A;\n"
					"      int tess_u = 64;\n"
					"      int tess_v = 64;\n"
					"      float du = 1.0f / float(tess_u);\n"
					"      float dv = 1.0f / float(tess_v);\n"
					"      float a = ((Zuu_max * X_max) / pow(Z_min, 2)) + (Xuu_max / Z_min);\n"
					"      float b = ((Zvv_max * X_max) / pow(Z_min, 2)) + (Xvv_max / Z_min);\n"
					"      float c = ((Zuv_max * X_max) / pow(Z_min, 2)) + (Xuv_max / Z_min);\n"
					"      float error = pow(du, 2) * a + pow(dv, 2) * b + 2 * du * dv * c;\n"
					" \n"
					"      while(target_error >= error) {\n"
					"        int new_tess_u = max(1, tess_u - 1);\n"
					"        int new_tess_v = max(1, tess_v - 1);\n"
					"        float du_ = 1.0f / float(new_tess_u);\n"
					"        float dv_ = 1.0f / float(new_tess_v);\n"
					"        //error = pow(du_, 2) * a + pow(dv_, 2) * b + 2 * du_ * dv_ * c;\n"
					"        float error_u = pow(du_, 2) * a + pow(dv, 2) * b + 2 * du_ * dv * c;\n"
					"        float error_v = pow(du, 2) * a + pow(dv_, 2) * b + 2 * du * dv_ * c;\n"
					"        if(error_u < error_v && error_u < target_error && tess_u != 1) {\n"
					"          tess_u = new_tess_u;\n"
					"          du = du_;\n"
					"          error = error_u;\n"
					"        }\n"
					"        else if(error_v < target_error && tess_v != 1) {\n"
					"          tess_v = new_tess_v;\n"
					"          dv = dv_;\n"
					"          error = error_v;\n"
					"        }\n"
					"        else {\n"
					"          break;\n"
					"        }\n"
					"      }\n"
					" \n"
					"      // Set tess levels, needs to cooperate with neighbours to set\n"
					"      // outer levels, for now, just set it as the same as the inner.\n"
					"      gl_TessLevelInner[0] = tess_u;\n"
					"      gl_TessLevelOuter[1] = tess_u;\n"
					"      gl_TessLevelOuter[3] = tess_u;\n"
					"      gl_TessLevelInner[1] = tess_v;\n"
					"      gl_TessLevelOuter[0] = tess_v;\n"
					"      gl_TessLevelOuter[2] = tess_v;\n"
					"    }\n"
					"  }\n"
					"}";
			}
			Shaders::LoadSpirv(name, src, shaderc_shader_kind::shaderc_tess_control_shader);
		}
		return { name, Shaders::SpirvMap[name] };
	}

	std::string Shaders::GetTeseShaderName(
		LocalSurfaceType lsType, TeseShaderType teseType, EvaluationMethod evalMethod, ShaderOptions& options)
	{
		std::string name = "Tese_";

		if (evalMethod == EvaluationMethod::Direct)
		{
			name += "Direct_";

			switch (teseType)
			{
			case TeseShaderType::Lattice:			name += "Lattice_"; break;
			case TeseShaderType::Local:				name += "Local_"; break;
			case TeseShaderType::Pixel_Accuracy:	name += "LatPixAcc_"; break;
			case TeseShaderType::Surf_Accuracy:		name += "LatSurfAcc_"; break;
			case TeseShaderType::Normals:			name += "LateNormals_"; break;
			}

			switch (lsType)
			{
			case LocalSurfaceType::Quadratic_Bezier:  name += "QuadBezier_"; break;
			case LocalSurfaceType::Cubic_Bezier: name += "CubicBezier_"; break;
			case LocalSurfaceType::Plane: name += "Plane_"; break;
			}
		} 
		else 
		{
			switch (evalMethod)
			{
			case EvaluationMethod::Pre_Sampled_Image:			name += "PreImg_"; break;
			case EvaluationMethod::Pre_Sampled_Image_Batched:	name += "PreImgBatch_"; break;
			case EvaluationMethod::Pre_Sampled_Buffer:			name += "PreBuf_"; break;
			}

			switch (teseType)
			{
			case TeseShaderType::Lattice:			name += "Lattice_"; break;
			case TeseShaderType::Local:				name += "Local_"; break;
			case TeseShaderType::Pixel_Accuracy:	name += "LatPixAcc_"; break;
			case TeseShaderType::Surf_Accuracy:		name += "LatSurfAcc_"; break;
			case TeseShaderType::Normals:			name += "LatNormals_"; break;
			}

			name += std::to_string(options.numSamplesU) + "x" + std::to_string(options.numSamplesV) + "S";
		}

		name += std::to_string(options.numControl) + "C_" + std::to_string(options.numLocal) + "L_"
			+ std::to_string(options.numPatches) + "P";

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
				desc += "layout(set = 2, binding = 0) uniform sampler2DArray localSampler;\n\n";
			}
			else {
				desc += Shaders::QuadUniformSampler2DArrayLayout("1", 0) + "\n\n";
			}
		}
		else if (evalMethod == EvaluationMethod::Pre_Sampled_Image_Batched)
		{
			desc += "layout(set = 1, binding = 0) uniform sampler2DArray surfSampler;\n\n";
		}
		else if (evalMethod == EvaluationMethod::Pre_Sampled_Buffer)
		{
			std::string nLocalSurfaceData = std::to_string(options.numLocal * options.numSamplesU * options.numSamplesV * 3);
			desc += Shaders::LocalSurfaceDataBufferString(nLocalSurfaceData, "0", "4") + "\n\n";
		}

		switch (teseType)
		{
		case TeseShaderType::Normals:
		case TeseShaderType::Local:
		case TeseShaderType::Lattice: {
			desc +=
				"layout(location = 0) in LocalSurfaceInfo tcLSInfo[];\n"
				"layout(location = 5) in vec3 tcColor[];\n"
				" \n"
				"layout(location = 0) out vec3 teColor;\n"
				"layout(location = 1) out vec3 teNormal;\n"
				"layout(location = 2) out vec3 tePosition;\n";
			break;
		}
		case TeseShaderType::Surf_Accuracy: {
			desc +=
				"layout(location = 0) in LocalSurfaceInfo tcLSInfo[];\n"
				"layout(location = 5) in vec3 tcColor[];\n"
				" \n"
				"layout(location = 0) out vec3 teColor;\n"
				"layout(location = 1) out vec3 teNormal;\n";
			break;
		}
		case TeseShaderType::Pixel_Accuracy: {
			desc +=
				"layout(location = 0) in LocalSurfaceInfo tcLSInfo[];\n"
				"layout(location = 5) in vec3 tcColor[];\n"
				" \n"
				"layout(location = 0) out vec3 teNormal;\n"
				"layout(location = 1) out vec2 teUVCoords;\n"
				"layout(location = 6) out LocalSurfaceInfo ls00Info;\n"
				"layout(location = 11) out LocalSurfaceInfo ls10Info;\n"
				"layout(location = 16) out LocalSurfaceInfo ls01Info;\n"
				"layout(location = 21) out LocalSurfaceInfo ls11Info;\n";
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
				Shaders::BlendingFunctionString() + "\n\n";
		}

		std::string nSampU = std::to_string(options.numSamplesU);
		std::string nSampV = std::to_string(options.numSamplesV);

		if (teseType == TeseShaderType::Surf_Accuracy)
		{
			std::string refEval;
			//refEval += Shaders::ClipToWindowFunctionString() + " \n";

			switch (lsType)
			{
			case LocalSurfaceType::Quadratic_Bezier: refEval += Shaders::BiQuadEvaluatorString("ref_"); break;
			case LocalSurfaceType::Cubic_Bezier: refEval += Shaders::BiCubicEvaluatorString("ref_"); break;
			case LocalSurfaceType::Plane: refEval += Shaders::PlaneEvaluatorString("ref_"); break;
			}

			switch (evalMethod)
			{
			case EvaluationMethod::Pre_Sampled_Image: {
				eval += refEval + " \n";
				eval += Shaders::SampleLocalSurfaceEvaluator("");
				break;
			}
			case EvaluationMethod::Pre_Sampled_Image_Batched: {
				eval += refEval + " \n";
				eval += Shaders::BatchedSampleLocalSurfaceEvaluator("");
				break;
			}
			case EvaluationMethod::Pre_Sampled_Buffer: {
				eval += refEval + " \n";
				eval += Shaders::SampleBufferEvaluator("", nSampU, nSampV);
				break;
			}
			}

		}
		else
		{
			switch (evalMethod)
			{
			case EvaluationMethod::Direct: {
				switch (lsType)
				{
				case LocalSurfaceType::Quadratic_Bezier: eval += Shaders::BiQuadEvaluatorString(""); break;
				case LocalSurfaceType::Cubic_Bezier: eval += Shaders::BiCubicEvaluatorString(""); break;
				case LocalSurfaceType::Plane: eval += Shaders::PlaneEvaluatorString(""); break;
				}
				break;
			}
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
			"  float u = gl_TessCoord.x;\n"
			"  float v = gl_TessCoord.y;\n \n";

		if (teseType == TeseShaderType::Local)
		{
			if (evalMethod == EvaluationMethod::Pre_Sampled_Image)
			{
				main += "  Sampler local = evaluateLocal(tcLSInfo[0], localSampler, u, v);\n";
			}
			else
			{
				main += "  Sampler local = evaluateLocal(tcLSInfo[0], u, v);\n";
			}
			main +=
				"  gl_Position = ubo.projection * ubo.modelview * vec4(local.p, 1.0f);\n"
				"  tePosition = vec3(ubo.modelview * vec4(local.p, 1.0f));\n"
				"  teNormal = vec3(ubo.normal * vec4(normalize(cross(local.u, local.v)), 0.0f));\n"
				"  teColor = tcColor[0];\n"
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
				Shaders::BlendLocalSurfacesString("") + " \n" +
				"  gl_Position = ubo.projection * ubo.modelview * vec4(pos, 1.0f);\n"
				"  teNormal = vec3(ubo.normal * vec4(normalize(cross(dpdu, dpdv)), 0.0f));\n \n"
				"  float ref_u = gl_TessCoord.x;\n"
				"  float ref_v = gl_TessCoord.y;\n \n"
				"  // Evaluate Local Surfaces\n"
				"  Sampler ref_s00 = ref_evaluateLocal(tcLSInfo[0], ref_u, ref_v);\n"
				"  Sampler ref_s10 = ref_evaluateLocal(tcLSInfo[1], ref_u, ref_v);\n"
				"  Sampler ref_s01 = ref_evaluateLocal(tcLSInfo[2], ref_u, ref_v);\n"
				"  Sampler ref_s11 = ref_evaluateLocal(tcLSInfo[3], ref_u, ref_v);\n \n"
				"  // Evaluate tensor product blending spline surface\n"
				"  float ref_Bu = 1.0 - bFunction(ref_u);\n"
				"  float ref_Bv = 1.0 - bFunction(ref_v);\n"
				" \n"
				"  vec3 ref_Su0 = ref_s10.p + (ref_s00.p - ref_s10.p) * ref_Bu;\n"
				"  vec3 ref_Su1 = ref_s11.p + (ref_s01.p - ref_s11.p) * ref_Bu;\n"
				" \n"
				"  vec3 ref_pos = ref_Su1 + (ref_Su0 - ref_Su1) * ref_Bv;\n \n"
				//"  vec4 clip = ubo.projection * ubo.modelview * vec4(pos, 1.0f);\n"
				//"  vec4 ref_clip = ubo.projection * ubo.modelview * vec4(ref_pos, 1.0f);\n"
				//"  vec3 pos_v = vec3(ubo.modelview * vec4(pos, 1.0f));\n"
				//"  vec3 ref_pos_v = vec3(ubo.modelview * vec4(ref_pos, 1.0f));\n"
				//"  float diff = distance(clipToWindow(clip), clipToWindow(ref_clip));\n" + " \n" +
				"  float diff = distance(pos, ref_pos);\n \n" +
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
				"  gl_Position = ubo.projection * ubo.modelview * vec4(pos, 1.0f);\n";

			if (teseType == TeseShaderType::Pixel_Accuracy)
			{
				main +=
					"  teNormal = vec3(ubo.normal * vec4(normalize(cross(dpdu, dpdv)), 0.0f));\n"
					"  teUVCoords = vec2(u, v);\n"
					"  ls00Info = tcLSInfo[0];\n"
					"  ls10Info = tcLSInfo[1];\n"
					"  ls01Info = tcLSInfo[2];\n"
					"  ls11Info = tcLSInfo[3];\n";
			}
			else
			{
				main +=
					"  tePosition = vec3(ubo.modelview * vec4(pos, 1.0f));\n"
					"  teColor = tcColor[0];\n";
				if (teseType == TeseShaderType::Normals)
				{
					main += "  teNormal = normalize(cross(dpdu, dpdv));\n";
				}
				else
				{
					main += "  teNormal = vec3(ubo.normal * vec4(normalize(cross(dpdu, dpdv)), 0.0f));\n";
				}
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

	NameSpirvPair Shaders::GetLatticeNormalsGeomShader()
	{
		if (Shaders::NotInMap(Shaders::BiQuadLatticeNormalsGeomName))
		{
			std::string src =
				Shaders::ShaderHeader() + " \n" +
				"layout(triangles) in;\n"
				" \n"
				"layout(location = 1) in vec3 teNormal[];\n"
				"layout(location = 2) in vec3 tePosition[];\n"
				" \n" +
				Shaders::CommonUniformBufferString() + " \n" +
				" \n"
				"layout(line_strip, max_vertices = 6) out;\n"
				" \n"
				"layout(location = 0) out vec3 gColor;\n"
				" \n"
				"void main() {\n"
				"  for(int vertex = 0; vertex < 3; ++vertex) {\n"
				"    gl_Position = ubo.projection * vec4(tePosition[vertex], 1);\n"
				"    gColor = vec3(0.2);\n"
				"    EmitVertex();\n"
				" \n"
				"    gl_Position = ubo.projection * (vec4(tePosition[vertex], 1)\n"
				"      + (transpose(inverse(ubo.modelview)) * vec4(teNormal[vertex] * ubo.normalLength, 0)));\n"
				"    gColor = vec3(0.6);\n"
				"    EmitVertex();\n"
				" \n"
				"    EndPrimitive();\n"
				"  }\n"
				"}";
			Shaders::LoadSpirv(Shaders::BiQuadLatticeNormalsGeomName, src, shaderc_shader_kind::shaderc_geometry_shader);
		}
		return { Shaders::BiQuadLatticeNormalsGeomName, Shaders::SpirvMap[Shaders::BiQuadLatticeNormalsGeomName] };
	}

	NameSpirvPair Shaders::GetTriSizeGeomShader()
	{
		if (Shaders::NotInMap(Shaders::TriSizeGeomName))
		{
			std::string src =
				Shaders::ShaderHeader() + " \n" +
				"layout(triangles) in;\n"
				" \n"
				"layout(location = 1) in vec3 teNormal[];\n"
				" \n" +
				Shaders::CommonUniformBufferString() + " \n" +
				" \n"
				"layout(triangle_strip, max_vertices = 3) out;\n"
				" \n"
				"layout(location = 0) out vec3 gColor;\n"
				"layout(location = 1) out vec3 gNormal;\n"
				" \n" +
				Shaders::ClipToWindowFunctionString() + "\n" + 
				"void main() {\n"
				"  vec4 p1 = gl_in[0].gl_Position;\n"
				"  vec4 p2 = gl_in[1].gl_Position;\n"
				"  vec4 p3 = gl_in[2].gl_Position;\n"
				" \n"
				"  // Set the color based on the longest edge of the triangle,\n"
				"  // when projected to screen space.\n"
				"  vec4 p1_w = clipToWindow(p1);\n"
				"  vec4 p2_w = clipToWindow(p2);\n"
				"  vec4 p3_w = clipToWindow(p3);\n"
				"  float maxEdgeLength = max(distance(p1_w, p2_w), \n"
				"    max(distance(p2_w, p3_w), distance(p1_w, p3_w)));\n"
				"  vec3 color = vec3(1, 0, 0);\n"
				"  if(maxEdgeLength < 1.0f) {\n"
				"    color = vec3(0.6, 0.6, 0.6);\n"
				"  }\n"
				"  else if(maxEdgeLength < 5.0f) {\n"
				"    color = vec3(0.0, 1.0, 0.0);\n"
				"  }\n"
				"  else if(maxEdgeLength < 10.0f) {\n"
				"    color = vec3(0.0, 0.0, 1.0);\n"
				"  }\n"
				"  else if(maxEdgeLength < 20.0f) {\n"
				"    color = vec3(1.0, 1.0, 0.0);\n"
				"  }\n"
				" \n"
				"  // Just pass the positions through the shader\n"
				"  gl_Position = p1;\n"
				"  gColor = color;\n"
				"  gNormal = teNormal[0];\n"
				"  EmitVertex();\n"
				"  gl_Position = p2;\n"
				"  gColor = color;\n"
				"  gNormal = teNormal[1];\n"
				"  EmitVertex();\n"
				"  gl_Position = p3;\n"
				"  gColor = color;\n"
				"  gNormal = teNormal[2];\n"
				"  EmitVertex();\n"
				"  EndPrimitive();\n"
				"}";
			Shaders::LoadSpirv(Shaders::TriSizeGeomName, src, shaderc_shader_kind::shaderc_geometry_shader);
		}
		return { Shaders::TriSizeGeomName, Shaders::SpirvMap[Shaders::TriSizeGeomName] };
	}

	NameSpirvPair Shaders::GetUintFlatColorFragShader()
	{
		if (Shaders::NotInMap(Shaders::UintFlatColorFragName))
		{
			std::string src =
				Shaders::ShaderHeader() +
				" \n"
				"layout(location = 0) flat in uvec3 vCol;\n"
				" \n"
				"layout(location = 0) out vec4 FragColor;\n"
				" \n"
				"void main() {\n"
				"  if(vCol.r == 255 && vCol.g == 255 && vCol.b == 255) {\n"
				"    discard;\n"
				"  }\n"
				"  FragColor = vec4(vCol / 255.0f, 1.0f);\n"
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
				"layout(location = 0) in vec3 teColor;\n"
				" \n"
				"layout(location = 0) out vec4 FragColor;\n"
				" \n"
				"void main() {\n"
				"  FragColor = vec4(teColor, 1.0);\n"
				"}";
			Shaders::LoadSpirv(Shaders::FlatColorFragName, src, shaderc_shader_kind::shaderc_fragment_shader);
		}
		return { Shaders::FlatColorFragName, Shaders::SpirvMap[Shaders::FlatColorFragName] };
	}

	NameSpirvPair Shaders::GetFlatColorNoInterpFragShader()
	{
		if (Shaders::NotInMap(Shaders::FlatColorNoInterpFragName))
		{
			std::string src =
				Shaders::ShaderHeader() +
				" \n"
				"layout(location = 0) flat in vec3 teColor;\n"
				" \n"
				"layout(location = 0) out vec4 FragColor;\n"
				" \n"
				"void main() {\n"
				"  FragColor = vec4(teColor, 1.0);\n"
				"}";
			Shaders::LoadSpirv(Shaders::FlatColorNoInterpFragName, src, shaderc_shader_kind::shaderc_fragment_shader);
		}
		return { Shaders::FlatColorNoInterpFragName, Shaders::SpirvMap[Shaders::FlatColorNoInterpFragName] };
	}

	NameSpirvPair Shaders::GetShadedColorFragShader()
	{
		if (Shaders::NotInMap(Shaders::ShadedColorFragName))
		{
			std::string src = 
				Shaders::ShaderHeader() +
				" \n"
				"layout(location = 0) in vec3 teColor;\n"
				"layout(location = 1) in vec3 teNormal;\n"
				" \n"
				"layout(location = 0) out vec4 FragColor;\n"
				" \n"
				"void main() {\n"
				"  vec3 light = normalize(vec3 (0.5f, 0.5f, 1.0f));\n"
				"  vec3 ambient = vec3(0.3f);\n"
				"  vec3 diffuse = max(dot(light, teNormal), 0.0f) * vec3(0.8f);\n"
				" \n"
				"  FragColor = vec4((ambient + diffuse) * teColor, 1.0);\n"
				"}";
			Shaders::LoadSpirv(Shaders::ShadedColorFragName, src, shaderc_shader_kind::shaderc_fragment_shader);
		}
		return { Shaders::ShadedColorFragName, Shaders::SpirvMap[Shaders::ShadedColorFragName] };
	}

	NameSpirvPair Shaders::GetShadedColorNoInterpFragShader()
	{
		if (Shaders::NotInMap(Shaders::ShadedColorNoInterpFragName))
		{
			std::string src =
				Shaders::ShaderHeader() +
				" \n"
				"layout(location = 0) flat in vec3 teColor;\n"
				"layout(location = 1) in vec3 teNormal;\n"
				" \n"
				"layout(location = 0) out vec4 FragColor;\n"
				" \n"
				"void main() {\n"
				"  vec3 light = normalize(vec3 (0.5f, 0.5f, 1.0f));\n"
				"  vec3 ambient = vec3(0.3f);\n"
				"  vec3 diffuse = max(dot(light, teNormal), 0.0f) * vec3(0.8f);\n"
				" \n"
				"  FragColor = vec4((ambient + diffuse) * teColor, 1.0);\n"
				"}";
			Shaders::LoadSpirv(Shaders::ShadedColorNoInterpFragName, src, shaderc_shader_kind::shaderc_fragment_shader);
		}
		return { Shaders::ShadedColorNoInterpFragName, Shaders::SpirvMap[Shaders::ShadedColorNoInterpFragName] };
	}

	NameSpirvPair Shaders::GetBiQuadLatticePixelAccuracyFragShader(LocalSurfaceType lsType,
		uint32_t numLocalsurfaceControlPoints, uint32_t numLocalSurfaces, uint32_t numPatches)
	{
		std::string nControl = std::to_string(numLocalsurfaceControlPoints);
		std::string nLocal = std::to_string(numLocalSurfaces);
		std::string nPatches = std::to_string(numPatches);
		std::string name;
		switch (lsType)
		{
		case LocalSurfaceType::Quadratic_Bezier:  name = "Frag_QuadBezier_"; break;
		case LocalSurfaceType::Cubic_Bezier:	  name = "Frag_CubicBezier_"; break;
		case LocalSurfaceType::Plane:			  name = "Frag_Plane_"; break;
		}
		name += "PixelAccuracy_" + nControl + "C_" + nLocal + "L_" + nPatches + "P";
		std::string evalString;
		switch (lsType)
		{
		case LocalSurfaceType::Quadratic_Bezier: evalString = Shaders::BiQuadEvaluatorString(""); break;
		case LocalSurfaceType::Cubic_Bezier: evalString = Shaders::BiCubicEvaluatorString(""); break;
		case LocalSurfaceType::Plane: evalString = Shaders::PlaneEvaluatorString(""); break;
		}
		if (Shaders::NotInMap(name))
		{
			std::string nBoundaries = std::to_string(numPatches * 4);
			std::string src =
				Shaders::ShaderHeader() + " \n" +
				Shaders::LocalSurfaceInfoStruct() + " \n" +
				Shaders::SamplerStruct() + " \n" +
				Shaders::BoundaryInfoStruct() + " \n" +
				" \n"
				"layout(location = 0) in vec3 teNormal;\n"
				"layout(location = 1) in vec2 teUVCoords;\n"
				"layout(location = 6) flat in LocalSurfaceInfo ls00Info;\n"
				"layout(location = 11) flat in LocalSurfaceInfo ls10Info;\n"
				"layout(location = 16) flat in LocalSurfaceInfo ls01Info;\n"
				"layout(location = 21) flat in LocalSurfaceInfo ls11Info;\n"
				" \n"
				"layout(location = 0) out vec4 FragColor;\n"
				" \n" +
				Shaders::CommonUniformBufferString() + " \n" +
				Shaders::MatrixStorageBufferString(nLocal, "0", "1") + " \n" +
				Shaders::ControlPointStorageBufferString(nControl, "0", "2") + " \n" +
				Shaders::BoundaryBufferString(nBoundaries, "0", "3") + " \n" +
				Shaders::BlendingFunctionString() + " \n" +
				evalString + " \n" +
				Shaders::ClipToWindowFunctionString() + "\n" +
				"void main() {\n"
				"  float u = teUVCoords[0];\n"
				"  float v = teUVCoords[1];\n"
				" \n" +
				"  Sampler s00 = evaluateLocal(ls00Info, u, v);\n"
				"  Sampler s10 = evaluateLocal(ls10Info, u, v);\n"
				"  Sampler s01 = evaluateLocal(ls01Info, u, v);\n"
				"  Sampler s11 = evaluateLocal(ls11Info, u, v);\n"
				" \n"
				"  float Bu = 1.0 - bFunction(u);\n"
				"  float Bv = 1.0 - bFunction(v);\n"
				" \n"
				"  vec3 Su0 = s10.p + (s00.p - s10.p) * Bu;\n"
				"  vec3 Su1 = s11.p + (s01.p - s11.p) * Bu;\n"
				"  vec3 p = Su1 + (Su0 - Su1) * Bv;\n"
				" \n"
				"  vec4 pos = ubo.projection * ubo.modelview * vec4(p, 1.0f);\n"
				"  vec4 windowPos = clipToWindow(pos);\n"
				" \n"
				"  float error = distance(gl_FragCoord.xy, windowPos.xy);\n"
				"  vec3 color;\n"
				" \n"
				"  if (error <= 0.1) {\n"
				"    color = vec3(0.7, 0.7, 0.7);\n"
				"  }\n"
				"  else if (error <= 0.5) {\n"
				"    color = vec3(0.0, 1.0, 0.0);\n"
				"  }\n"
				"  else if (error <= 2.0) {\n"
				"    color = vec3(0.0, 0.0, 1.0);\n"
				"  }\n"
				"  else if (error <= 5.0) {\n"
				"    color = vec3(1.0, 1.0, 0.0);\n"
				"  }\n"
				"  else {\n"
				"    color = vec3(1.0, 0.0, 0.0);\n"
				"  }\n"
				" \n"
				"  FragColor = vec4(color, 1.0);\n"
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



	// Just sorts the vector of shaderNames alphabetically
	// For nicer display in the shader viewer
	void Shaders::SortShaderNames()
	{
		std::sort(Shaders::ShaderNames.begin(), Shaders::ShaderNames.end(),
			[](std::string a, std::string b) { return a < b; });
	}

	// Common shader code
	inline std::string Shaders::ShaderHeader()
	{ 
		return "#version 450\n";
	}

	inline std::string Shaders::TeseHeader()
	{
		return "layout(quads, fractional_even_spacing, ccw) in;\n";
	}

	inline std::string Shaders::CommonUniformBufferString()
	{
		return
			"layout(set = 0, binding = 0) uniform UniformBufferObject {\n"
			"  int tessInner;\n"
			"  int tessOuter;\n"
			"  int tessFactorMethod;\n"
			"  int pixelsPerEdge;\n"
			"  int numPASamples;\n"
			"  float maxError;\n"
			"  float normalLength;\n"
			"  mat4 projection;\n"
			"  mat4 modelview;\n"
			"  mat4 normal;\n"
			"  vec2 windowSize;\n"
			"} ubo;\n";
	}

	inline std::string Shaders::LocalSurfaceInfoStruct() 
	{
		return 
			"struct LocalSurfaceInfo {\n"
			"  uint controlPointIndex;\n"
			"  uint matrixIndex;\n"
			"  uint boundaryIndex;\n"
			"  uvec3 batchInfo;\n"
			"  uint dataIndex;\n"
			"};\n";
	}

	inline std::string Shaders::SamplerStruct() 
	{
		return 
			"struct Sampler {\n"
			"  vec3 p;\n"
			"  vec3 u;\n"
			"  vec3 v;\n"
			"};\n";
	}

	inline std::string Shaders::BoundaryInfoStruct() 
	{
		return
			"struct BoundaryInfo {\n"
			"  float us;\n"
			"  float ue;\n"
			"  float vs;\n"
			"  float ve;\n"
			"};\n";
	}

	inline std::string Shaders::QuadUniformSampler2DArrayLayout(std::string set, uint32_t firstBinding)
	{
		return
			"layout(set = " + set + ", binding = " + std::to_string(firstBinding) + ") uniform sampler2DArray p00Sampler;\n"
			"layout(set = " + set + ", binding = " + std::to_string(firstBinding+1) + ") uniform sampler2DArray p10Sampler;\n"
			"layout(set = " + set + ", binding = " + std::to_string(firstBinding+2) + ") uniform sampler2DArray p01Sampler;\n"
			"layout(set = " + set + ", binding = " + std::to_string(firstBinding+3) + ") uniform sampler2DArray p11Sampler;\n";
	}

	inline std::string Shaders::MatrixStorageBufferString(
		std::string numLocalSurfaces, std::string set, std::string binding) 
	{
		return
			"layout(set = " + set + ", binding = 1) buffer MatrixBuffer {\n"
			"  mat4 matrices[" + numLocalSurfaces + "];\n"
			"} matrixBuffer;\n";
	}

	inline std::string Shaders::ControlPointStorageBufferString(
		std::string numControlPoints, std::string set, std::string binding)
	{
		return
			"layout(set = " + set + ", binding = " + binding + ") buffer ControlPointBuffer {\n"
			"  vec4 controlPoints[" + numControlPoints + "];\n"
			"} controlPointBuffer;\n";
	}

	inline std::string Shaders::BoundaryBufferString(
		std::string numBoundaries, std::string set, std::string binding)
	{
		return
			"layout(set = " + set + ", binding = " + binding + ") buffer BoundaryBuffer {\n"
			"  BoundaryInfo boundaries[" + numBoundaries + "];\n"
			"} boundaryBuffer;\n";
	}

	inline std::string Shaders::LocalSurfaceDataBufferString(
		std::string numDataPoints, std::string set, std::string binding)
	{
		return
			"layout(set = " + set + ", binding = " + binding + ") buffer LocalSurfaceDataBuffer {\n"
			"  vec4 localSurfaceData[" + numDataPoints + "];\n"
			"} lsDataBuffer;\n";
	}

	inline std::string Shaders::BlendingFunctionString()
	{
		return
			"float bFunction(float w) {\n"
			"  return 6 * pow(w, 5) - 15 * pow(w, 4) + 10 * pow(w, 3);\n"
			"}\n"
			" \n"
			"float bDerivative(float w) {\n"
			"  return 30 * pow(w, 4) - 60 * pow(w, 3) + 30 * pow(w, 2);\n"
			"}\n";
	}

	inline std::string Shaders::BiQuadEvaluatorString(std::string prefix)
	{
		return
			"vec3 " + prefix + "bez3Basis(float t) {\n"
			"  return vec3(pow(1-t, 2), 2*t*(1-t), pow(t, 2));\n"
			"}\n"
			" \n"
			"vec3 " + prefix + "bez3BasisDer(float t) {\n"
			"  return vec3(2*t-2, 2-4*t, 2*t);\n"
			"}\n"
			" \n"
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"  float local_u = mix(bi.us, bi.ue, u);\n"
			"  float local_v = mix(bi.vs, bi.ve, v);\n"
			" \n"
			"  vec3 bu = " + prefix + "bez3Basis(local_u);\n"
			"  vec3 bud = " + prefix + "bez3BasisDer(local_u);\n"
			"  vec3 bv = " + prefix + "bez3Basis(local_v);\n"
			"  vec3 bvd = " + prefix + "bez3BasisDer(local_v);\n"
			" \n"
			"  vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"  vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"  vec3 p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"  vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			"  vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz;\n"
			"  vec3 p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;\n"
			"  vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz;\n"
			"  vec3 p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz;\n"
			"  vec3 p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;\n"
			" \n"
			"  vec3 pos =\n"
			"    p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] +\n"
			"    p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] +\n"
			"    p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2];\n"
			" \n"
			"  vec3 dpdu =\n"
			"    p00 * bud[0] * bv[0] + p01 * bud[0] * bv[1] + p02 * bud[0] * bv[2] +\n"
			"    p10 * bud[1] * bv[0] + p11 * bud[1] * bv[1] + p12 * bud[1] * bv[2] +\n"
			"    p20 * bud[2] * bv[0] + p21 * bud[2] * bv[1] + p22 * bud[2] * bv[2];\n"
			" \n"
			"  vec3 dpdv =\n"
			"    p00 * bu[0] * bvd[0] + p01 * bu[0] * bvd[1] + p02 * bu[0] * bvd[2] +\n"
			"    p10 * bu[1] * bvd[0] + p11 * bu[1] * bvd[1] + p12 * bu[1] * bvd[2] +\n"
			"    p20 * bu[2] * bvd[0] + p21 * bu[2] * bvd[1] + p22 * bu[2] * bvd[2];\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = vec3(matrix * vec4(pos, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  dpdu = vec3(matrix * vec4(dpdu, 0.0f));\n"
			"  dpdv = vec3(matrix * vec4(dpdv, 0.0f));\n"
			" \n"
			"  return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);\n"
			"}\n";
	}

	inline std::string Shaders::BiQuadEvaluatorOnlyPosString(std::string prefix)
	{
		return
			"vec3 " + prefix + "bez3Basis(float t) {\n"
			"  return vec3(pow (1-t, 2), 2*t*(1-t), pow(t, 2));\n"
			"}\n"
			" \n"
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"  float local_u = mix(bi.us, bi.ue, u);\n"
			"  float local_v = mix(bi.vs, bi.ve, v);\n"
			" \n"
			"  vec3 bu = " + prefix + "bez3Basis(local_u);\n"
			"  vec3 bv = " + prefix + "bez3Basis(local_v);\n"
			" \n"
			"  vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"  vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"  vec3 p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"  vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			"  vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz;\n"
			"  vec3 p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;\n"
			"  vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz;\n"
			"  vec3 p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz;\n"
			"  vec3 p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;\n"
			" \n"
			"  vec3 pos =\n"
			"    p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] +\n"
			"    p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] +\n"
			"    p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2];\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  vec4 p = matrix * vec4(pos, 1.0f);\n"
			" \n"
			"  return Sampler(p.xyz, vec3(0.0), vec3(0.0));\n"
			"}\n";
	}

	std::string Shaders::BiQuadEvaluatorWithDer2()
	{
		return
			"vec3 bez3Basis(float t) {\n"
			"  return vec3(pow(1-t, 2), 2*t*(1-t), pow(t, 2));\n"
			"}\n"
			" \n"
			"vec3 bez3BasisDer(float t) {\n"
			"  return vec3(2*t-2, 2-4*t, 2*t);\n"
			"}\n"
			" \n"
			"vec3 bez3BasisDer2(float t) {\n"
			"  return vec3(4, -8, 4);\n"
			"}\n"
			" \n"
			"Sampler2 evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"  float local_u = mix(bi.us, bi.ue, u);\n"
			"  float local_v = mix(bi.vs, bi.ve, v);\n"
			" \n"
			"  vec3 bu = bez3Basis(local_u);\n"
			"  vec3 bud = bez3BasisDer(local_u);\n"
			"  vec3 bud2 = bez3BasisDer2(local_u);\n"
			"  vec3 bv = bez3Basis(local_v);\n"
			"  vec3 bvd = bez3BasisDer(local_v);\n"
			"  vec3 bvd2 = bez3BasisDer2(local_v);\n"
			" \n"
			"  vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"  vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"  vec3 p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"  vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			"  vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz;\n"
			"  vec3 p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;\n"
			"  vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz;\n"
			"  vec3 p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz;\n"
			"  vec3 p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;\n"
			" \n"
			"  vec3 pos =\n"
			"    p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] +\n"
			"    p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] +\n"
			"    p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2];\n"
			" \n"
			"  vec3 dpdu =\n"
			"    p00 * bud[0] * bv[0] + p01 * bud[0] * bv[1] + p02 * bud[0] * bv[2] +\n"
			"    p10 * bud[1] * bv[0] + p11 * bud[1] * bv[1] + p12 * bud[1] * bv[2] +\n"
			"    p20 * bud[2] * bv[0] + p21 * bud[2] * bv[1] + p22 * bud[2] * bv[2];\n"
			" \n"
			"  vec3 dpdv =\n"
			"    p00 * bu[0] * bvd[0] + p01 * bu[0] * bvd[1] + p02 * bu[0] * bvd[2] +\n"
			"    p10 * bu[1] * bvd[0] + p11 * bu[1] * bvd[1] + p12 * bu[1] * bvd[2] +\n"
			"    p20 * bu[2] * bvd[0] + p21 * bu[2] * bvd[1] + p22 * bu[2] * bvd[2];\n"
			" \n"
			"  vec3 dpdudv =\n"
			"    p00 * bud[0] * bvd[0] + p01 * bud[0] * bvd[1] + p02 * bud[0] * bvd[2] +\n"
			"    p10 * bud[1] * bvd[0] + p11 * bud[1] * bvd[1] + p12 * bud[1] * bvd[2] +\n"
			"    p20 * bud[2] * bvd[0] + p21 * bud[2] * bvd[1] + p22 * bud[2] * bvd[2];\n"
			" \n"
			"  vec3 dpdudu =\n"
			"    p00 * bud2[0] * bv[0] + p01 * bud2[0] * bv[1] + p02 * bud2[0] * bv[2] +\n"
			"    p10 * bud2[1] * bv[0] + p11 * bud2[1] * bv[1] + p12 * bud2[1] * bv[2] +\n"
			"    p20 * bud2[2] * bv[0] + p21 * bud2[2] * bv[1] + p22 * bud2[2] * bv[2];\n"
			" \n"
			"  vec3 dpdvdv =\n"
			"    p00 * bu[0] * bvd2[0] + p01 * bu[0] * bvd2[1] + p02 * bu[0] * bvd2[2] +\n"
			"    p10 * bu[1] * bvd2[0] + p11 * bu[1] * bvd2[1] + p12 * bu[1] * bvd2[2] +\n"
			"    p20 * bu[2] * bvd2[0] + p21 * bu[2] * bvd2[1] + p22 * bu[2] * bvd2[2];\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = vec3(matrix * vec4(pos, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  dpdu = vec3(matrix * vec4(dpdu, 0.0f));\n"
			"  dpdv = vec3(matrix * vec4(dpdv, 0.0f));\n"
			"  dpdudv = vec3(matrix * vec4(dpdudv, 0.0f));\n"
			"  dpdudu = vec3(matrix * vec4(dpdudu, 0.0f));\n"
			"  dpdvdv = vec3(matrix * vec4(dpdvdv, 0.0f));\n"
			" \n"
			"  return Sampler2(pos, dpdu, dpdv, dpdudv, dpdudu, dpdvdv);\n"
			"}\n";
	}

	std::string Shaders::BiCubicEvaluatorString(std::string prefix)
	{
		return
			"vec4 " + prefix + "bez4Basis(float t) {\n"
			"  return vec4(pow(1-t, 3), 3*t*pow((1-t),2), 3*t*t*(1-t), pow(t, 3));\n"
			"}\n"
			" \n"
			"vec4 " + prefix + "bez4BasisDer(float t) {\n"
			"  return vec4(-3*pow(1-t, 2), 3*(1-t)*(1-3*t), 3*t*(2-3*t), 3*t*t);\n"
			"}\n"
			" \n"
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"  float local_u = mix(bi.us, bi.ue, u);\n"
			"  float local_v = mix(bi.vs, bi.ve, v);\n"
			" \n"
			"  vec4 bu = " + prefix + "bez4Basis(local_u);\n"
			"  vec4 bud = " + prefix + "bez4BasisDer(local_u);\n"
			"  vec4 bv = " + prefix + "bez4Basis(local_v);\n"
			"  vec4 bvd = " + prefix + "bez4BasisDer(local_v);\n"
			" \n"
			"  vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"  vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"  vec3 p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"  vec3 p30 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			"  vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz;\n"
			"  vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;\n"
			"  vec3 p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz;\n"
			"  vec3 p31 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz;\n"
			"  vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;\n"
			"  vec3 p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 9].xyz;\n"
			"  vec3 p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 10].xyz;\n"
			"  vec3 p32 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 11].xyz;\n"
			"  vec3 p03 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 12].xyz;\n"
			"  vec3 p13 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 13].xyz;\n"
			"  vec3 p23 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 14].xyz;\n"
			"  vec3 p33 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 15].xyz;\n"
			" \n"
			"  vec3 pos =\n"
			"    p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] + p03 * bu[0] * bv[3] +\n"
			"    p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] + p13 * bu[1] * bv[3] +\n"
			"    p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2] + p23 * bu[2] * bv[3] +\n"
			"    p30 * bu[3] * bv[0] + p31 * bu[3] * bv[1] + p32 * bu[3] * bv[2] + p33 * bu[3] * bv[3];\n"
			" \n"
			"  vec3 dpdu =\n"
			"    p00 * bud[0] * bv[0] + p01 * bud[0] * bv[1] + p02 * bud[0] * bv[2] + p03 * bud[0] * bv[3] +\n"
			"    p10 * bud[1] * bv[0] + p11 * bud[1] * bv[1] + p12 * bud[1] * bv[2] + p13 * bud[1] * bv[3] +\n"
			"    p20 * bud[2] * bv[0] + p21 * bud[2] * bv[1] + p22 * bud[2] * bv[2] + p23 * bud[2] * bv[3] +\n"
			"    p30 * bud[3] * bv[0] + p31 * bud[3] * bv[1] + p32 * bud[3] * bv[2] + p33 * bud[3] * bv[3];\n"
			" \n"
			"  vec3 dpdv =\n"
			"    p00 * bu[0] * bvd[0] + p01 * bu[0] * bvd[1] + p02 * bu[0] * bvd[2] + p03 * bu[0] * bvd[3] +\n"
			"    p10 * bu[1] * bvd[0] + p11 * bu[1] * bvd[1] + p12 * bu[1] * bvd[2] + p13 * bu[1] * bvd[3] +\n"
			"    p20 * bu[2] * bvd[0] + p21 * bu[2] * bvd[1] + p22 * bu[2] * bvd[2] + p23 * bu[2] * bvd[3] +\n"
			"    p30 * bu[3] * bvd[0] + p31 * bu[3] * bvd[1] + p32 * bu[3] * bvd[2] + p33 * bu[3] * bvd[3];\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = vec3(matrix * vec4(pos, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  dpdu = vec3(matrix * vec4(dpdu, 0.0f));\n"
			"  dpdv = vec3(matrix * vec4(dpdv, 0.0f));\n"
			" \n"
			"  return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);\n"
			"}\n";
	}

	std::string Shaders::BiCubicEvaluatorOnlyPosString(std::string prefix)
	{
		return
			"vec4 " + prefix + "bez4Basis(float t) {\n"
			"  return vec4(pow(1-t, 3), 3*t* pow((1-t),2), 3*t*t*(1-t), pow(t, 3));\n"
			"}\n"
			" \n"
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"  float local_u = mix(bi.us, bi.ue, u);\n"
			"  float local_v = mix(bi.vs, bi.ve, v);\n"
			" \n"
			"  vec4 bu = " + prefix + "bez4Basis(local_u);\n"
			"  vec4 bv = " + prefix + "bez4Basis(local_v);\n"
			" \n"
			"  vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"  vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"  vec3 p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"  vec3 p30 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			"  vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz;\n"
			"  vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;\n"
			"  vec3 p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz;\n"
			"  vec3 p31 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz;\n"
			"  vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;\n"
			"  vec3 p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 9].xyz;\n"
			"  vec3 p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 10].xyz;\n"
			"  vec3 p32 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 11].xyz;\n"
			"  vec3 p03 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 12].xyz;\n"
			"  vec3 p13 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 13].xyz;\n"
			"  vec3 p23 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 14].xyz;\n"
			"  vec3 p33 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 15].xyz;\n"
			" \n"
			"  vec3 pos =\n"
			"    p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] + p03 * bu[0] * bv[3] +\n"
			"    p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] + p13 * bu[1] * bv[3] +\n"
			"    p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2] + p23 * bu[2] * bv[3] +\n"
			"    p30 * bu[3] * bv[0] + p31 * bu[3] * bv[1] + p32 * bu[3] * bv[2] + p33 * bu[3] * bv[3];\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = matrix * vec4(pos, 1.0f);\n"
			" \n"
			"  return Sampler(pos.xyz, vec3(0.0), vec3(0.0));\n"
			"}\n";
	}

	std::string Shaders::BiCubicEvaluatorWithDer2()
	{
		return
			"vec4 bez4Basis(float t) {\n"
			"  return vec4(pow(1-t, 3), 3*t*pow((1-t),2), 3*t*t*(1-t), pow(t, 3));\n"
			"}\n"
			" \n"
			"vec4 bez4BasisDer(float t) {\n"
			"  return vec4(-3*pow(1-t, 2), 3*(1-t)*(1-3*t), 3*t*(2-3*t), 3*t*t);\n"
			"}\n"
			" \n"
			"vec4 bez4BasisDer2(float t) {\n"
			"  return vec4(6-6*t, 18*t - 12, 6 - 18 * t, 6 * t);\n"
			"}\n"
			" \n"
			"Sampler2 evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"  float local_u = mix(bi.us, bi.ue, u);\n"
			"  float local_v = mix(bi.vs, bi.ve, v);\n"
			" \n"
			"  vec4 bu = bez4Basis(local_u);\n"
			"  vec4 bud = bez4BasisDer(local_u);\n"
			"  vec4 bud2 = bez4BasisDer2(local_u);\n"
			"  vec4 bv = bez4Basis(local_v);\n"
			"  vec4 bvd = bez4BasisDer(local_v);\n"
			"  vec4 bvd2 = bez4BasisDer2(local_v);\n"
			" \n"
			"  vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"  vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"  vec3 p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"  vec3 p30 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			"  vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz;\n"
			"  vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;\n"
			"  vec3 p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz;\n"
			"  vec3 p31 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz;\n"
			"  vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;\n"
			"  vec3 p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 9].xyz;\n"
			"  vec3 p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 10].xyz;\n"
			"  vec3 p32 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 11].xyz;\n"
			"  vec3 p03 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 12].xyz;\n"
			"  vec3 p13 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 13].xyz;\n"
			"  vec3 p23 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 14].xyz;\n"
			"  vec3 p33 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 15].xyz;\n"
			" \n"
			"  vec3 pos =\n"
			"    p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] + p03 * bu[0] * bv[3] +\n"
			"    p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] + p13 * bu[1] * bv[3] +\n"
			"    p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2] + p23 * bu[2] * bv[3] +\n"
			"    p30 * bu[3] * bv[0] + p31 * bu[3] * bv[1] + p32 * bu[3] * bv[2] + p33 * bu[3] * bv[3];\n"
			" \n"
			"  vec3 dpdu =\n"
			"    p00 * bud[0] * bv[0] + p01 * bud[0] * bv[1] + p02 * bud[0] * bv[2] + p03 * bud[0] * bv[3] +\n"
			"    p10 * bud[1] * bv[0] + p11 * bud[1] * bv[1] + p12 * bud[1] * bv[2] + p13 * bud[1] * bv[3] +\n"
			"    p20 * bud[2] * bv[0] + p21 * bud[2] * bv[1] + p22 * bud[2] * bv[2] + p23 * bud[2] * bv[3] +\n"
			"    p30 * bud[3] * bv[0] + p31 * bud[3] * bv[1] + p32 * bud[3] * bv[2] + p33 * bud[3] * bv[3];\n"
			" \n"
			"  vec3 dpdv =\n"
			"    p00 * bu[0] * bvd[0] + p01 * bu[0] * bvd[1] + p02 * bu[0] * bvd[2] + p03 * bu[0] * bvd[3] +\n"
			"    p10 * bu[1] * bvd[0] + p11 * bu[1] * bvd[1] + p12 * bu[1] * bvd[2] + p13 * bu[1] * bvd[3] +\n"
			"    p20 * bu[2] * bvd[0] + p21 * bu[2] * bvd[1] + p22 * bu[2] * bvd[2] + p23 * bu[2] * bvd[3] +\n"
			"    p30 * bu[3] * bvd[0] + p31 * bu[3] * bvd[1] + p32 * bu[3] * bvd[2] + p33 * bu[3] * bvd[3];\n"
			" \n"
			"  vec3 dpdudv =\n"
			"    p00 * bud[0] * bvd[0] + p01 * bud[0] * bvd[1] + p02 * bud[0] * bvd[2] + p03 * bud[0] * bvd[3] +\n"
			"    p10 * bud[1] * bvd[0] + p11 * bud[1] * bvd[1] + p12 * bud[1] * bvd[2] + p13 * bud[1] * bvd[3] +\n"
			"    p20 * bud[2] * bvd[0] + p21 * bud[2] * bvd[1] + p22 * bud[2] * bvd[2] + p23 * bud[2] * bvd[3] +\n"
			"    p30 * bud[3] * bvd[0] + p31 * bud[3] * bvd[1] + p32 * bud[3] * bvd[2] + p33 * bud[3] * bvd[3];\n"
			" \n"
			"  vec3 dpdudu =\n"
			"    p00 * bud2[0] * bv[0] + p01 * bud2[0] * bv[1] + p02 * bud2[0] * bv[2] + p03 * bud2[0] * bv[3] +\n"
			"    p10 * bud2[1] * bv[0] + p11 * bud2[1] * bv[1] + p12 * bud2[1] * bv[2] + p13 * bud2[1] * bv[3] +\n"
			"    p20 * bud2[2] * bv[0] + p21 * bud2[2] * bv[1] + p22 * bud2[2] * bv[2] + p23 * bud2[2] * bv[3] +\n"
			"    p30 * bud2[3] * bv[0] + p31 * bud2[3] * bv[1] + p32 * bud2[3] * bv[2] + p33 * bud2[3] * bv[3];\n"
			" \n"
			"  vec3 dpdvdv =\n"
			"    p00 * bu[0] * bvd2[0] + p01 * bu[0] * bvd2[1] + p02 * bu[0] * bvd2[2] + p03 * bu[0] * bvd2[3] +\n"
			"    p10 * bu[1] * bvd2[0] + p11 * bu[1] * bvd2[1] + p12 * bu[1] * bvd2[2] + p13 * bu[1] * bvd2[3] +\n"
			"    p20 * bu[2] * bvd2[0] + p21 * bu[2] * bvd2[1] + p22 * bu[2] * bvd2[2] + p23 * bu[2] * bvd2[3] +\n"
			"    p30 * bu[3] * bvd2[0] + p31 * bu[3] * bvd2[1] + p32 * bu[3] * bvd2[2] + p33 * bu[3] * bvd2[3];\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = vec3(matrix * vec4(pos, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  dpdu = vec3(matrix * vec4(dpdu, 0.0f));\n"
			"  dpdv = vec3(matrix * vec4(dpdv, 0.0f));\n"
			"  dpdudv = vec3(matrix * vec4(dpdudv, 0.0f));\n"
			"  dpdudu = vec3(matrix * vec4(dpdudu, 0.0f));\n"
			"  dpdvdv = vec3(matrix * vec4(dpdvdv, 0.0f));\n"
			" \n"
			"  return Sampler2(pos, dpdu, dpdv, dpdudv, dpdudu, dpdvdv);\n"
			"}\n";
	}

	std::string Shaders::PlaneEvaluatorString(std::string prefix)
	{
		return
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"  float local_u = mix(bi.us, bi.ue, u);\n"
			"  float local_v = mix(bi.vs, bi.ve, v);\n"
			" \n"
			"  vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"  vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"  vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"  vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			" \n"
			"  vec3 pos = mix(mix(p00, p10, local_u), mix(p01, p11, local_u), local_v);\n"
			" \n"
			"  vec3 dpdu = p10 - p00;\n"
			" \n"
			"  vec3 dpdv = p01 - p00;\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = vec3(matrix * vec4(pos, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  dpdu = vec3(matrix * vec4(dpdu, 0.0f));\n"
			"  dpdv = vec3(matrix * vec4(dpdv, 0.0f));\n"
			" \n"
			"  return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);\n"
			"}\n";
	}

	std::string Shaders::PlaneEvaluatorWithDer2()
	{
		return
			"Sampler2 evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			"  float local_u = mix(bi.us, bi.ue, u);\n"
			"  float local_v = mix(bi.vs, bi.ve, v);\n"
			" \n"
			"  vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz;\n"
			"  vec3 p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz;\n"
			"  vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;\n"
			"  vec3 p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz;\n"
			" \n"
			"  vec3 pos = mix(mix(p00, p10, local_u), mix(p01, p11, local_u), local_v);\n"
			" \n"
			"  vec3 dpdu = p10 - p00;\n"
			" \n"
			"  vec3 dpdv = p01 - p00;\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = vec3(matrix * vec4(pos, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  dpdu = vec3(matrix * vec4(dpdu, 0.0f));\n"
			"  dpdv = vec3(matrix * vec4(dpdv, 0.0f));\n"
			" \n"
			"  return Sampler2(pos.xyz, dpdu.xyz, dpdv.xyz, vec3(0.0f), vec3(0.0f), vec3(0.0f));\n"
			"}\n";
	}

	inline std::string Shaders::SampleLocalSurfaceEvaluator(std::string prefix)
	{
		return
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, sampler2DArray surfSampler, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			" \n"
			"  vec2 coords = vec2(mix(bi.us, bi.ue, u), mix(bi.vs, bi.ve, v));\n"
			" \n"
			"  vec3 pos = texture(surfSampler, vec3(coords, 0)).xyz;\n"
			"  vec3 dpdu = texture(surfSampler, vec3(coords, 1)).xyz;\n"
			"  vec3 dpdv = texture(surfSampler, vec3(coords, 2)).xyz;\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = vec3 (matrix * vec4(pos, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  dpdu = vec3(matrix * vec4(dpdu, 0.0f));\n"
			"  dpdv = vec3(matrix * vec4(dpdv, 0.0f));\n"
			" \n"
			"  return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);\n"
			"}\n";
	}

	inline std::string Shaders::BatchedSampleLocalSurfaceEvaluator(std::string prefix)
	{
		return
			"Sampler evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v, int layer) {\n"
			"  ivec3 dim = textureSize(surfSampler, 0);\n"
			"  float ss = float(lsInfo.batchInfo[0] + 0.5);\n"
			"  float se = float(lsInfo.batchInfo[0] + lsInfo.batchInfo[2] - 0.5);\n"
			"  float ts = float(lsInfo.batchInfo[1] + 0.5);\n"
			"  float te = float(lsInfo.batchInfo[1] + lsInfo.batchInfo[2] - 0.5);\n"
			"  float s = mix(ss, se, u) / float(dim.x);\n"
			"  float t = mix(ts, te, v) / float(dim.y);\n"
			" \n"
			"  vec3 pos = texture(surfSampler, vec3(s, t, layer)).xyz;\n"
			"  vec3 dpdu = texture(surfSampler, vec3(s, t, layer + 1)).xyz;\n"
			"  vec3 dpdv = texture(surfSampler, vec3(s, t, layer + 2)).xyz;\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  pos = vec3(matrix * vec4(pos, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  dpdu = vec3(matrix * vec4(dpdu, 0.0f));\n"
			"  dpdv = vec3(matrix * vec4(dpdv, 0.0f));\n"
			" \n"
			"  return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);\n"
			"}\n";
	}

	inline std::string Shaders::SampleBufferEvaluator(std::string prefix, std::string nSampU, std::string nSampV)
	{
		return
			"// Takes in a u and v coordinate and returns a Sampler object (p,u,v)\n"
			"// Uses linear interpolation to mix the neighbouring values when the u and\n"
			"// v coordinates are between sample values\n"
			"Sampler sampleBuffer(int idx, float u, float v) {\n"
			"  Sampler res;\n"
			"  // Step one, transform u,v coordinates [0,1]^2 to buffer coordinates [0,numSamplesU]x[0,numSamplesV]\n"
			"  float fi = u * (" + nSampU + " - 1);\n"
			"  float fj = v * (" + nSampV + " - 1);\n"
			"  int i = int(floor(fi));\n"
			"  int j = int(floor(fj));\n"
			"  float fu = fract(fi);\n"
			"  float fv = fract(fj);\n"
			" \n"
			"  int c00 = j * " + nSampU + " + i;\n"
			"  int c10 = j * " + nSampU + " + i + 1;\n"
			"  int c01 = (j + 1) * " + nSampU + " + i;\n"
			"  int c11 = (j + 1) * " + nSampU + " + i + 1;\n"
			" \n"
			"  // Step two, linearly interpolate between neighbouring values\n"
			"  vec3 p00 = lsDataBuffer.localSurfaceData[idx + c00].xyz;\n"
			"  vec3 p10 = lsDataBuffer.localSurfaceData[idx + c10].xyz;\n"
			"  vec3 p01 = lsDataBuffer.localSurfaceData[idx + c01].xyz;\n"
			"  vec3 p11 = lsDataBuffer.localSurfaceData[idx + c11].xyz;\n"
			"  res.p = (((1 - fu) * p00 + fu * p10) * (1 - fv) + ((1 - fu) * p01 + fu * p11) * fv);\n"
			" \n"
			"  int offsetU = " + nSampU + " * " + nSampV + ";\n"
			"  vec3 u00 = lsDataBuffer.localSurfaceData[idx + c00 + offsetU].xyz;\n"
			"  vec3 u10 = lsDataBuffer.localSurfaceData[idx + c10 + offsetU].xyz;\n"
			"  vec3 u01 = lsDataBuffer.localSurfaceData[idx + c01 + offsetU].xyz;\n"
			"  vec3 u11 = lsDataBuffer.localSurfaceData[idx + c11 + offsetU].xyz;\n"
			"  res.u = (((1 - fu) * u00 + fu * u10) * (1 - fv) + ((1 - fu) * u01 + fu * u11) * fv);\n"
			" \n"
			"  int offsetV = offsetU * 2;\n"
			"  vec3 v00 = lsDataBuffer.localSurfaceData[idx + c00 + offsetV].xyz;\n"
			"  vec3 v10 = lsDataBuffer.localSurfaceData[idx + c10 + offsetV].xyz;\n"
			"  vec3 v01 = lsDataBuffer.localSurfaceData[idx + c01 + offsetV].xyz;\n"
			"  vec3 v11 = lsDataBuffer.localSurfaceData[idx + c11 + offsetV].xyz;\n"
			"  res.v = (((1 - fu) * v00 + fu * v10) * (1 - fv) + ((1 - fu) * v01 + fu * v11) * fv);\n"
			" \n"
			"  // Return the sampled values\n"
			"  return res;\n"
			"}\n"
			" \n"
			"Sampler " + prefix + "evaluateLocal(LocalSurfaceInfo lsInfo, float u, float v) {\n"
			"  BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];\n"
			" \n"
			"  Sampler s = sampleBuffer(int(lsInfo.dataIndex * " + nSampU + " * " + nSampV + " * 3), \n"
			"    mix(bi.us, bi.ue, u), mix(bi.vs, bi.ve, v));\n"
			" \n"
			"  mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];\n"
			" \n"
			"  s.p = vec3(matrix * vec4(s.p, 1.0f));\n"
			"  // mat4 normalMat = transpose(inverse(matrix));\n"
			"  s.u = vec3(matrix * vec4(s.u, 0.0f));\n"
			"  s.v = vec3(matrix * vec4(s.v, 0.0f));\n"
			" \n"
			"  return s;\n"
			"}\n";
	}

	inline std::string Shaders::EvaluateLocalSurfacesString(std::string prefix)
	{
		return
			"  // Evaluate Local Surfaces\n"
			"  Sampler " + prefix + "s00 = " + prefix + "evaluateLocal(tcLSInfo[0], u, v);\n"
			"  Sampler " + prefix + "s10 = " + prefix + "evaluateLocal(tcLSInfo[1], u, v);\n"
			"  Sampler " + prefix + "s01 = " + prefix + "evaluateLocal(tcLSInfo[2], u, v);\n"
			"  Sampler " + prefix + "s11 = " + prefix + "evaluateLocal(tcLSInfo[3], u, v);\n";
	}

	inline std::string Shaders::SampleLocalSurfacesString(std::string prefix)
	{
		return
			"  // Evaluate Local Surfaces\n"
			"  Sampler " + prefix + "s00 = " + prefix + "evaluateLocal(tcLSInfo[0], p00Sampler, u, v);\n"
			"  Sampler " + prefix + "s10 = " + prefix + "evaluateLocal(tcLSInfo[1], p10Sampler, u, v);\n"
			"  Sampler " + prefix + "s01 = " + prefix + "evaluateLocal(tcLSInfo[2], p01Sampler, u, v);\n"
			"  Sampler " + prefix + "s11 = " + prefix + "evaluateLocal(tcLSInfo[3], p11Sampler, u, v);\n";
	}

	std::string Shaders::SampleBatchedLocalSurfacesString(std::string prefix)
	{
		return
			"  // Evaluate Local Surfaces\n"
			"  Sampler " + prefix + "s00 = " + prefix + "evaluateLocal(tcLSInfo[0], u, v, 0);\n"
			"  Sampler " + prefix + "s10 = " + prefix + "evaluateLocal(tcLSInfo[1], u, v, 3);\n"
			"  Sampler " + prefix + "s01 = " + prefix + "evaluateLocal(tcLSInfo[2], u, v, 6);\n"
			"  Sampler " + prefix + "s11 = " + prefix + "evaluateLocal(tcLSInfo[3], u, v, 9);\n";
	}

	inline std::string Shaders::BlendLocalSurfacesString(std::string prefix)
	{
		return
			"  // Evaluate tensor product blending spline surface\n"
			"  float " + prefix + "Bu = 1.0 - bFunction(u);\n"
			"  float " + prefix + "Bv = 1.0 - bFunction(v);\n"
			"  float " + prefix + "Bu1 = -bDerivative(u);\n"
			"  float " + prefix + "Bv1 = -bDerivative(v);\n"
			" \n"
			"  vec3 " + prefix + "Su0 = " + prefix + "s10.p + (" + prefix + "s00.p - " + prefix + "s10.p) * " + prefix + "Bu;\n"
			"  vec3 " + prefix + "Su1 = " + prefix + "s11.p + (" + prefix + "s01.p - " + prefix + "s11.p) * " + prefix + "Bu;\n"
			" \n"
			"  vec3 " + prefix + "Uu0 = " + prefix + "s10.u + (" + prefix + "s00.u - " + prefix + "s10.u) * " + prefix + "Bu\n"
			"    + (" + prefix + "s00.p - " + prefix + "s10.p) * " + prefix + "Bu1;\n"
			"  vec3 " + prefix + "Uu1 = " + prefix + "s11.u + (" + prefix + "s01.u - " + prefix + "s11.u) * " + prefix + "Bu\n"
			"    + (" + prefix + "s01.p - " + prefix + "s11.p) * " + prefix + "Bu1;\n"
			" \n"
			"  vec3 " + prefix + "Vu0 = " + prefix + "s10.v + (" + prefix + "s00.v - " + prefix + "s10.v) * " + prefix + "Bu;\n"
			"  vec3 " + prefix + "Vu1 = " + prefix + "s11.v + (" + prefix + "s01.v - " + prefix + "s11.v) * " + prefix + "Bu;\n"
			" \n"
			"  vec3 " + prefix + "pos = " + prefix + "Su1 + (" + prefix + "Su0 - " + prefix + "Su1) * " + prefix + "Bv;\n"
			"  vec3 " + prefix + "dpdu = " + prefix + "Uu1 + (" + prefix + "Uu0 - " + prefix + "Uu1) * " + prefix + "Bv;\n"
			"  vec3 " + prefix + "dpdv = " + prefix + "Vu1 + (" + prefix + "Vu0 - " + prefix + "Vu1) * " + prefix + "Bv\n"
			"    + (" + prefix + "Su0 - " + prefix + "Su1) * " + prefix + "Bv1;\n";
	}

	std::string Shaders::GetPixelAccGlobalSurfEvaluatorString()
	{
		return
			"float bDerivative2(float w) {\n"
			"  return 120 * pow(w, 3) - 180 * pow(w, 2) + 60 * w;\n"
			"}\n"
			" \n"
			"struct Sampler2 {\n"
			"  vec3 p, u, v, uv, uu, vv;\n"
			"};\n"
			" \n"
			"Sampler2 evaluateGlobal(Sampler2 s00, Sampler2 s10,\n"
			"  Sampler2 s01, Sampler2 s11, float u, float v) {\n"
			"  Sampler2 res;\n"
			"  float Bu = 1.0 - bFunction(u);\n"
			"  float Bv = 1.0 - bFunction(v);\n"
			"  float Bu1 = -bDerivative(u);\n"
			"  float Bv1 = -bDerivative(v);\n"
			"  float Bu2 = -bDerivative2(u);\n"
			"  float Bv2 = -bDerivative2(v);\n"
			" \n"
			"  vec3 Su0 = s10.p + (s00.p - s10.p) * Bu;\n"
			"  vec3 Su1 = s11.p + (s01.p - s11.p) * Bu;\n"
			" \n"
			"  vec3 Uu0 = s10.u + (s00.u - s10.u) * Bu\n"
			"    + (s00.p - s10.p) * Bu1;\n"
			"  vec3 Uu1 = s11.u + (s01.u - s11.u) * Bu\n"
			"    + (s01.p - s11.p) * Bu1;\n"
			" \n"
			"  vec3 Vu0 = s10.v + (s00.v - s10.v) * Bu;\n"
			"  vec3 Vu1 = s11.v + (s01.v - s11.v) * Bu;\n"
			" \n"
			"  vec3 UUu0 = s10.uu + (s00.uu - s10.uu) * Bu\n"
			"    + 2 * (s00.u - s10.u) * Bu1 + (s00.p - s10.p) * Bu2;\n"
			"  vec3 UUu1 = s11.uu + (s01.uu - s11.uu) * Bu\n"
			"    + 2 * (s01.u - s11.u) * Bu1 + (s01.p - s11.p) * Bu2;\n"
			" \n"
			"  vec3 VVu0 = s10.vv + (s00.vv - s10.vv) * Bu;\n"
			"  vec3 VVu1 = s11.vv + (s01.vv - s11.vv) * Bu;\n"
			" \n"
			"  vec3 UVu0 = s10.uv + (s00.uv - s10.uv) * Bu\n"
			"    + (s00.p - s10.p) * Bu1;\n"
			"  vec3 UVu1 = s11.uv + (s01.uv - s11.uv) * Bu\n"
			"    + (s01.p - s11.p) * Bu1;\n"
			" \n"
			"  res.p = Su1 + (Su0 - Su1) * Bv;\n"
			"  res.u = Uu1 + (Uu0 - Uu1) * Bv;\n"
			"  res.v = Vu1 + (Vu0 - Vu1) * Bv + (Su0 - Su1) * Bv1;\n"
			"  res.uv = UVu1 + (UVu0 - UVu1) * Bv + (Su0 - Su1) * Bv1;\n"
			"  res.uu = UUu1 + (UUu0 - UUu1) * Bv;\n"
			"  res.vv = VVu1 + (VVu0 - VVu1) * Bv + 2 * ((Vu0 - Vu1) * Bv1) + (Su0 - Su1) * Bv2;\n"
			" \n"
			"  return res;\n"
			"}\n";
	}

	inline std::string Shaders::EvaluateLocalSurfacesOnlyPosString(std::string prefix)
	{
		return
			"  // Evaluate Local Surfaces\n"
			"  Sampler " + prefix + "s00 = " + prefix + "evaluateLocal(tcLSInfo[0], u, v);\n"
			"  Sampler " + prefix + "s10 = " + prefix + "evaluateLocal(tcLSInfo[1], u, v);\n"
			"  Sampler " + prefix + "s01 = " + prefix + "evaluateLocal(tcLSInfo[2], u, v);\n"
			"  Sampler " + prefix + "s11 = " + prefix + "evaluateLocal(tcLSInfo[3], u, v);\n";
	}

	inline std::string Shaders::BlendLocalSurfacesOnlyPos(std::string prefix)
	{
		return
			"  // Evaluate tensor product blending spline surface\n"
			"  float " + prefix + "Bu = 1.0 - bFunction(u);\n"
			"  float " + prefix + "Bv = 1.0 - bFunction(v);\n"
			" \n"
			"  vec3 " + prefix + "Su0 = " + prefix + "s10.p + (" + prefix + "s00.p - " + prefix + "s10.p) * " + prefix + "Bu;\n"
			"  vec3 " + prefix + "Su1 = " + prefix + "s11.p + (" + prefix + "s01.p - " + prefix + "s11.p) * " + prefix + "Bu;\n"
			" \n"
			"  vec3 " + prefix + "pos = " + prefix + "Su1 + (" + prefix + "Su0 - " + prefix + "Su1) * " + prefix + "Bv;\n";
	}

	inline std::string Shaders::ColorByMaxErrorString()
	{
		return
			"  if(diff <= (ubo.maxError * 0.1f)) {\n"
			"    teColor = vec3(0.5, 0.5, 0.5);\n"
			"  }\n"
			"  else if (diff <= (ubo.maxError * 0.2f)) {\n"
			"    teColor = vec3(0.0, 1.0, 0.0);\n"
			"  }\n"
			"  else if (diff <= (ubo.maxError * 0.5f)) {\n"
			"    teColor = vec3(0.0, 0.0, 1.0);\n"
			"  }\n"
			"  else if (diff <= ubo.maxError) {\n"
			"    teColor = vec3(1.0, 1.0, 0.0);\n"
			"  }\n"
			"  else {\n"
			"    teColor = vec3(1.0, 0.0, 0.0);\n"
			"  }\n";
	}

	std::string Shaders::ClipToWindowFunctionString()
	{
		return
			"vec4 clipToWindow(vec4 clip) {\n"
			"  vec4 ndc = clip / clip.w;\n"
			"  float windowX = (ubo.windowSize.x / 2.0f) * (ndc.x + 1.0f);\n"
			"  float windowY = (ubo.windowSize.y / 2.0f) * (ndc.y + 1.0f);\n"
			"  float windowZ = 0.5f * ndc.z + 0.5f;\n"
			"  return vec4(windowX, windowY, windowZ, 1.0f);\n"
			"}\n";
	}

}