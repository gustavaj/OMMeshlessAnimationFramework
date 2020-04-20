#version 450

layout(quads, equal_spacing, ccw) in;

struct LocalSurfaceInfo
{
	uint controlPointIndex;
	uint controlPointCount;
	uint matrixIndex;
	uint boundaryIndex;
};

struct Sampler
{
	vec3 p, u, v;
};

struct BoundaryInfo
{
	float us, ue, vs, ve;
};

layout(location = 0) in LocalSurfaceInfo tcLSInfo[];

layout(location = 0) out vec3 teColor;

layout(set = 0, binding = 0) uniform LatticeUBO
{
	int tessInnerLevel;
	int tessOuterLevel;
	int bFunctionIndex;
	mat4 projection;
	mat4 modelview;
	mat4 normal;
	vec2 windowSize;
} latticeUbo;

layout(constant_id = 0) const int numLocalSurfaceControlPoints = 36000;
layout(constant_id = 1) const int numLocalSurfaces = 1000;
layout(constant_id = 2) const int numPatches = 1000;

layout(set = 0, binding = 1) buffer MatrixBuffer
{
	mat4 matrices[numLocalSurfaces];
} matrixBuffer;

layout(set = 0, binding = 2) buffer ControlPointBuffer
{
	vec4 controlPoints[numLocalSurfaceControlPoints];
} controlPointBuffer;

layout(set = 0, binding = 3) buffer BoundaryBuffer
{
	BoundaryInfo boundaries[numPatches * 4];
} boundaryBuffer;

layout(set = 2, binding = 0) uniform sampler2DArray localSampler;

Sampler evaluateBiquadraticBezier(LocalSurfaceInfo lsInfo, float u, float v)
{	
	BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];
	
	vec2 coords = vec2(mix(bi.us, bi.ue, u), mix(bi.vs, bi.ve, v));
	
	vec3 pos = texture(localSampler, vec3(coords, 0)).xyz;
	vec3 dpdu = texture(localSampler, vec3(coords, 1)).xyz;
	vec3 dpdv = texture(localSampler, vec3(coords, 2)).xyz;
			 
	mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];
	
	pos  = vec3(matrix * vec4(pos, 1.0f));
	//mat4 normalMat = transpose(inverse(matrix));
	dpdu = vec3(matrix * vec4(dpdu, 0.0f));
	dpdv = vec3(matrix * vec4(dpdv, 0.0f));
	
	return Sampler(pos.xyz, dpdu.xyz, dpdv.xyz);
}

void main()
{
	float u = gl_TessCoord.x;
	float v = gl_TessCoord.y;
	
	Sampler local = evaluateBiquadraticBezier(tcLSInfo[0], u, v);
		
	gl_Position = latticeUbo.projection * latticeUbo.modelview * vec4(local.p, 1.0f);
	
	// if(fac <= 0.5) {
		// fac = fac * 2;
		// teColor = vec3(1 - fac, fac, 0.0f);
	// }
	// else {
		// fac = (fac - 0.5) * 2;
		// teColor = vec3(1 - fac, 1 - fac, fac);
	// }
	teColor = vec3(u, 0, v);
}