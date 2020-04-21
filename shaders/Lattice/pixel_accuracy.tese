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

const int B1Poly = 0;
const int B2Poly = 1;
const int Lerbs = 2;

layout(location = 0) in LocalSurfaceInfo tcLSInfo[];
layout(location = 4) in vec3[] tcColor;

layout(location = 0) out vec3 teColor;
layout(location = 1) out vec3 teNormal;
layout(location = 2) out vec4 tePosition;
layout(location = 3) out vec2 teUVCoords;
layout(location = 4) out LocalSurfaceInfo ls00Info;
layout(location = 8) out LocalSurfaceInfo ls10Info;
layout(location = 12) out LocalSurfaceInfo ls01Info;
layout(location = 16) out LocalSurfaceInfo ls11Info;

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

float bDerivative(float w)
{
	if(w < 1e-5) return 0.0f;
	if(1 - w < 1e-5) return 0.0f;
	switch(latticeUbo.bFunctionIndex) {
		case B1Poly: { return 6*w-6*pow(w,2); }
		case B2Poly: { return 30 * pow(w, 4) - 60 * pow(w, 3) + 30 * pow(w, 2); }
		case Lerbs: {
			return (2 * exp( - ( (1.0) / ((w-1.0)*w) ) ) *
				(pow(w,2.0) - w + 0.5)) /
				(pow( exp( -(1.0/(w-1.0)) ) + exp(1.0/w), 2.0) *
				pow(w-1.0,2.0) * pow(w,2.0)); }
	}	
}

vec3 bezBasis(float t) {
	return vec3(pow(1-t, 2), 2*t*(1-t), pow(t, 2));
}

vec3 bezBasisDer(float t) {
	return vec3(2*t-2, 2-4*t, 2*t);
}

Sampler evaluateBiquadraticBezier(LocalSurfaceInfo lsInfo, float u, float v)
{	
	BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];
	float local_u = mix(bi.us, bi.ue, u);
	float local_v = mix(bi.vs, bi.ve, v);

	vec3 bu = bezBasis(local_u);
	vec3 bud = bezBasisDer(local_u);
	vec3 bv = bezBasis(local_v);
	vec3 bvd = bezBasisDer(local_v);
	
	vec3 p00 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 0].xyz, p10 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 1].xyz, p20 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 2].xyz;
	vec3 p01 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 3].xyz, p11 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 4].xyz, p21 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 5].xyz;
	vec3 p02 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 6].xyz, p12 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 7].xyz, p22 = controlPointBuffer.controlPoints[lsInfo.controlPointIndex + 8].xyz;
	
	vec3 pos = 
		p00 * bu[0] * bv[0] + p01 * bu[0] * bv[1] + p02 * bu[0] * bv[2] +
		p10 * bu[1] * bv[0] + p11 * bu[1] * bv[1] + p12 * bu[1] * bv[2] +
		p20 * bu[2] * bv[0] + p21 * bu[2] * bv[1] + p22 * bu[2] * bv[2];
			  
	vec3 dpdu = 
		p00 * bud[0] * bv[0] + p01 * bud[0] * bv[1] + p02 * bud[0] * bv[2] +
		p10 * bud[1] * bv[0] + p11 * bud[1] * bv[1] + p12 * bud[1] * bv[2] +
		p20 * bud[2] * bv[0] + p21 * bud[2] * bv[1] + p22 * bud[2] * bv[2];
			  
	vec3 dpdv = 
		p00 * bu[0] * bvd[0] + p01 * bu[0] * bvd[1] + p02 * bu[0] * bvd[2] +
		p10 * bu[1] * bvd[0] + p11 * bu[1] * bvd[1] + p12 * bu[1] * bvd[2] +
		p20 * bu[2] * bvd[0] + p21 * bu[2] * bvd[1] + p22 * bu[2] * bvd[2];
			 
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
	
	// Evaluate local surfaces
	Sampler s00 = evaluateBiquadraticBezier(tcLSInfo[0], u, v);
	Sampler s10 = evaluateBiquadraticBezier(tcLSInfo[1], u, v);
	Sampler s01 = evaluateBiquadraticBezier(tcLSInfo[2], u, v);
	Sampler s11 = evaluateBiquadraticBezier(tcLSInfo[3], u, v);
	
	// Evaluate tensor product blending spline surface
	float Bu = 1.0 - bFunction(u);
	float Bv = 1.0 - bFunction(v);
	float Bu1 = -bDerivative(u);
	float Bv1 = -bDerivative(v);
	
	vec3 Su0 = s10.p + (s00.p - s10.p) * Bu;
	vec3 Su1 = s11.p + (s01.p - s11.p) * Bu;

	vec3 Uu0 = s10.u + (s00.u - s10.u) * Bu + (s00.p - s10.p) * Bu1;
	vec3 Uu1 = s11.u + (s01.u - s11.u) * Bu + (s01.p - s11.p) * Bu1;

	vec3 Vu0 = s10.v + (s00.v - s10.v) * Bu;
	vec3 Vu1 = s11.v + (s01.v - s11.v) * Bu;
	
	vec3 pos  = Su1 + (Su0 - Su1) * Bv;
	vec3 dpdu = Uu1 + (Uu0 - Uu1) * Bv;
	vec3 dpdv = Vu1 + (Vu0 - Vu1) * Bv + (Su0 - Su1) * Bv1;
			  
	gl_Position = latticeUbo.projection * latticeUbo.modelview * vec4(pos, 1.0f);
	teNormal = vec3(latticeUbo.normal * vec4(normalize(cross(dpdu, dpdv)), 0.0f));
	tePosition = gl_Position;
	teUVCoords = vec2(u, v);
	
	ls00Info = tcLSInfo[0];
	ls10Info = tcLSInfo[1];
	ls01Info = tcLSInfo[2];
	ls11Info = tcLSInfo[3];
	
	// float fac = gl_PrimitiveID / float(numPatches);
	// if(fac < 0.5f) {
		// fac = fac * 2.0f;
		// teColor = vec3(1.0f - fac, fac, 0.0f);
	// }
	// else if (fac < 1.0f){
		// fac = (fac - 0.5f) * 2.0f;
		// teColor = vec3(1.0f - fac, 1.0f - fac, fac);
	// }
	
	teColor = tcColor[0];
}