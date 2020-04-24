#version 450

layout(quads, equal_spacing, ccw) in;

struct LocalSurfaceInfo
{
	uint s;
	uint t;
	uint matrixIndex;
	uint boundaryIndex;
};

struct Sampler
{
	vec3 p, u, v;
};

const int B1Poly = 0;
const int B2Poly = 1;
const int Lerbs = 2;

layout(location = 0) in LocalSurfaceInfo tcLSInfo[];
layout(location = 4) in vec3[] tcColor;

layout(location = 0) out vec3 teColor;
layout(location = 1) out vec3 teNormal;
layout(location = 2) out vec3 tePosition;

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

layout(constant_id = 1) const int numLocalSurfaces = 1000;
layout(constant_id = 3) const int numSamplesU = 64;
layout(constant_id = 4) const int numSamplesV = 64;

layout(set = 0, binding = 1) buffer MatrixBuffer
{
	mat4 matrices[numLocalSurfaces];
} matrixBuffer;

layout(set = 1, binding = 0) uniform sampler2DArray surfSampler;

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

Sampler evaluateBiquadraticBezier(LocalSurfaceInfo lsInfo, float u, float v, int layer)
{
	vec3 pos  = texture(surfSampler, vec3(u, v, layer)).xyz;
	vec3 dpdu = texture(surfSampler, vec3(u, v, layer + 1)).xyz;
	vec3 dpdv = texture(surfSampler, vec3(u, v, layer + 2)).xyz;
			 
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
	Sampler s00 = evaluateBiquadraticBezier(tcLSInfo[0], u, v, 0);
	Sampler s10 = evaluateBiquadraticBezier(tcLSInfo[1], u, v, 3);
	Sampler s01 = evaluateBiquadraticBezier(tcLSInfo[2], u, v, 6);
	Sampler s11 = evaluateBiquadraticBezier(tcLSInfo[3], u, v, 9);
	
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
	tePosition = vec3(latticeUbo.modelview * vec4(pos, 1.0f));
	teNormal = vec3(latticeUbo.normal * vec4(normalize(cross(dpdu, dpdv)), 0.0f));
	
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