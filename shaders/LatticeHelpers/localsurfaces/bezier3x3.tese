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

Sampler evaluateBiquadraticBezier(LocalSurfaceInfo lsInfo, float u, float v)
{	
	float bu0 = (1-u) * (1-u);
	float bu1 = 2 * u * (1-u);
	float bu2 = u * u;
	
	float dbu0 = -2 * (1-u);
	float dbu1 = 2 - 4 * u;
	float dbu2 = 2 * u;
	
	float bv0 = (1-v) * (1-v);
	float bv1 = 2 * v * (1-v);
	float bv2 = v * v;
	
	float dbv0 = -2 * (1-v);
	float dbv1 = 2 - 4 * v;
	float dbv2 = 2 * v;
	
	vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz, p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz, p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;
	vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz, p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz, p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;
	vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz, p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz, p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;
	
	vec3 pos = bv0 * (bu0 * p00 + bu1 * p10 + bu2 * p20)
			 + bv1 * (bu0 * p01 + bu1 * p11 + bu2 * p21)
			 + bv2 * (bu0 * p02 + bu1 * p12 + bu2 * p22);
			 
	vec3 dpdu = bv0 * (dbu0 * p00 + dbu1 * p10 + dbu2 * p20)
			  + bv1 * (dbu0 * p01 + dbu1 * p11 + dbu2 * p21)
			  + bv2 * (dbu0 * p02 + dbu1 * p12 + dbu2 * p22);
	
	vec3 dpdv = dbv0 * (bu0 * p00 + bu1 * p10 + bu2 * p20)
			  + dbv1 * (bu0 * p01 + bu1 * p11 + bu2 * p21)
			  + dbv2 * (bu0 * p02 + bu1 * p12 + bu2 * p22);
			 
	mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];
	
	pos  = vec3(matrix * vec4(pos, 1.0f));
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