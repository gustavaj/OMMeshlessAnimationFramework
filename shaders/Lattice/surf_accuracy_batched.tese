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

layout(set = 1, binding = 0) uniform sampler2DArray p00Sampler;
layout(set = 1, binding = 1) uniform sampler2DArray p10Sampler;
layout(set = 1, binding = 2) uniform sampler2DArray p01Sampler;
layout(set = 1, binding = 3) uniform sampler2DArray p11Sampler;

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

Sampler evaluateBiquadraticBezierDirect(LocalSurfaceInfo lsInfo, float u, float v)
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

Sampler evaluateBiquadraticBezierSampled(LocalSurfaceInfo lsInfo, sampler2DArray surfSampler, float u, float v)
{	
	BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];
	
	vec2 coords = vec2(mix(bi.us, bi.ue, u), mix(bi.vs, bi.ve, v));
	
	vec3 pos = texture(surfSampler, vec3(coords, 0)).xyz;
	vec3 dpdu = texture(surfSampler, vec3(coords, 1)).xyz;
	vec3 dpdv = texture(surfSampler, vec3(coords, 2)).xyz;
			 
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
	
	// Sample local surfaces
	Sampler s00_samp = evaluateBiquadraticBezierSampled(tcLSInfo[0], p00Sampler, u, v);
	Sampler s10_samp = evaluateBiquadraticBezierSampled(tcLSInfo[1], p10Sampler, u, v);
	Sampler s01_samp = evaluateBiquadraticBezierSampled(tcLSInfo[2], p01Sampler, u, v);
	Sampler s11_samp = evaluateBiquadraticBezierSampled(tcLSInfo[3], p11Sampler, u, v);
	
	// Evaluate tensor product blending spline surface
	float Bu = 1.0 - bFunction(u);
	float Bv = 1.0 - bFunction(v);
	float Bu1 = -bDerivative(u);
	float Bv1 = -bDerivative(v);
	
	// Evaluate from the sampled data
	vec3 Su0_samp = s10_samp.p + (s00_samp.p - s10_samp.p) * Bu;
	vec3 Su1_samp = s11_samp.p + (s01_samp.p - s11_samp.p) * Bu;

	vec3 Uu0_samp = s10_samp.u + (s00_samp.u - s10_samp.u) * Bu + (s00_samp.p - s10_samp.p) * Bu1;
	vec3 Uu1_samp = s11_samp.u + (s01_samp.u - s11_samp.u) * Bu + (s01_samp.p - s11_samp.p) * Bu1;

	vec3 Vu0_samp = s10_samp.v + (s00_samp.v - s10_samp.v) * Bu;
	vec3 Vu1_samp = s11_samp.v + (s01_samp.v - s11_samp.v) * Bu;
	
	vec3 pos_samp  = Su1_samp + (Su0_samp - Su1_samp) * Bv;
	vec3 dpdu = Uu1_samp + (Uu0_samp - Uu1_samp) * Bv;
	vec3 dpdv = Vu1_samp + (Vu0_samp - Vu1_samp) * Bv + (Su0_samp - Su1_samp) * Bv1;
			  
	gl_Position = latticeUbo.projection * latticeUbo.modelview * vec4(pos_samp, 1.0f);
	
	// Evaluate local surfaces direct
	Sampler s00_dir = evaluateBiquadraticBezierDirect(tcLSInfo[0], u, v);
	Sampler s10_dir = evaluateBiquadraticBezierDirect(tcLSInfo[1], u, v);
	Sampler s01_dir = evaluateBiquadraticBezierDirect(tcLSInfo[2], u, v);
	Sampler s11_dir = evaluateBiquadraticBezierDirect(tcLSInfo[3], u, v);
	
	// Evaluate from the direct data
	vec3 Su0_dir = s10_dir.p + (s00_dir.p - s10_dir.p) * Bu;
	vec3 Su1_dir = s11_dir.p + (s01_dir.p - s11_dir.p) * Bu;
	
	vec3 pos_dir  = Su1_dir + (Su0_dir - Su1_dir) * Bv;
	
	// Longest side
	// Sampler topLeft = evaluateBiquadraticBezierDirect(tcLSInfo[0], 0.0, 0.0);
	// Sampler topRight = evaluateBiquadraticBezierDirect(tcLSInfo[1], 1.0, 0.0);
	// Sampler bottomLeft = evaluateBiquadraticBezierDirect(tcLSInfo[2], 0.0, 1.0);
	// float longestSide = max(distance(topLeft.p, topRight.p), distance(topLeft.p, bottomLeft.p));
	
	// Difference between the two
	float diff = distance(pos_samp, pos_dir);
	float maxError = 1.0f;
	
	if(diff <= maxError * 0.1) {
		teColor = vec3(0.5, 0.5, 0.5);
	}
	else if(diff <= maxError * 0.2) {
		teColor = vec3(0.0, 1.0, 0.0);
	}
	else if(diff <= maxError * 0.5) {
		teColor = vec3(0.0, 0.0, 1.0);
	}
	else if(diff <= maxError) {
		teColor = vec3(1.0, 1.0, 0.0);
	}
	else {
		teColor = vec3(1.0, 0.0, 0.0);
	}
}