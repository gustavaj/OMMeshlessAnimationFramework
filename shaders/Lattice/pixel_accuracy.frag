#version 450

struct LocalSurfaceInfo
{
	uint controlPointIndex;
	uint controlPointCount;
	uint matrixIndex;
	uint boundaryIndex;
};

layout(location = 0) in vec3 teColor;
layout(location = 1) in vec3 teNormal;
layout(location = 2) in vec4 tePosition;
layout(location = 3) in vec2 teUVCoords;
layout(location = 4)  flat in LocalSurfaceInfo ls00Info;
layout(location = 8)  flat in LocalSurfaceInfo ls10Info;
layout(location = 12) flat in LocalSurfaceInfo ls01Info;
layout(location = 16) flat in LocalSurfaceInfo ls11Info;

layout(location = 0) out vec4 FragColor;

struct Sampler
{
	vec3 p, u, v;
};

struct BoundaryInfo
{
	float us, ue, vs, ve;
};

const int B1Poly = 0;
const int B2Poly = 1;
const int Lerbs = 2;

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

float bFunction(float w)
{
	if(w < 1e-5) return 0.0f;
	if(1 - w < 1e-5) return 1.0f;
	switch(latticeUbo.bFunctionIndex) {
		case B1Poly: { return (3*pow(w, 2)-2*pow(w,3)); }
		case B2Poly: { return 6 * pow(w, 5) - 15 * pow(w, 4) + 10 * pow(w, 3); }
		case Lerbs: { return 1.0 / (1.0 + exp(1.0 / w - 1.0 / (1.0 - w))); }
	}	
}

vec3 bezBasis(float t) {
	return vec3(pow(1-t, 2), 2*t*(1-t), pow(t, 2));
}

vec3 evaluateBiquadraticBezier(LocalSurfaceInfo lsInfo, float u, float v)
{	
	BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];
	float local_u = mix(bi.us, bi.ue, u);
	float local_v = mix(bi.vs, bi.ve, v);

	vec3 bu = bezBasis(local_u);
	vec3 bv = bezBasis(local_v);
	
	vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz, p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz, p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;
	vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz, p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz, p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;
	vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz, p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz, p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;
	
	vec3 pos = 
		p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] +
		p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] +
		p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2];
			 
	mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];
	
	return vec3(matrix * vec4(pos, 1.0f));
}

void main()
{
	float u = teUVCoords[0], v = teUVCoords[1];

	// Evaluate local surfaces
	vec3 s00 = evaluateBiquadraticBezier(ls00Info, u, v);
	vec3 s10 = evaluateBiquadraticBezier(ls10Info, u, v);
	vec3 s01 = evaluateBiquadraticBezier(ls01Info, u, v);
	vec3 s11 = evaluateBiquadraticBezier(ls11Info, u, v);
	
	// Evaluate tensor product blending spline surface
	float Bu = 1.0 - bFunction(u);
	float Bv = 1.0 - bFunction(v);
	
	vec3 Su0 = s10 + (s00 - s10) * Bu;
	vec3 Su1 = s11 + (s01 - s11) * Bu;
	vec3 p  = Su1 + (Su0 - Su1) * Bv;
			  
	vec4 pos = latticeUbo.projection * latticeUbo.modelview * vec4(p, 1.0f);
	vec4 ndc = pos / pos.w;
	float windowX = (latticeUbo.windowSize.x / 2) * (ndc.x + 1);
	float windowY = (latticeUbo.windowSize.y / 2) * (ndc.y + 1);
	float windowZ = 0.5 * ndc.z + 0.5;
	vec4 windowPos = vec4(windowX, windowY, windowZ, 1.0f);
	
	float error = distance(gl_FragCoord.xy, windowPos.xy);
	vec3 color;
	
	if(error <= 0.5) {
		color = vec3(0.7, 0.7, 0.7);
	}
	else if(error <= 1.0) {
		color = vec3(0.0, 1.0, 0.0);
	}
	else if(error <= 2.0) {
		color = vec3(0.0, 0.0, 1.0);
	}
	else if(error <= 5.0) {
		color = vec3(1.0, 1.0, 0.0);
	}
	else {
		color = vec3(1.0, 0.0, 0.0);
	}
	
    FragColor = vec4(color, 1.0);
}