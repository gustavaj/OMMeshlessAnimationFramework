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

layout(constant_id = 0) const int numLocalSurfaceControlPoints = 36000;
layout(constant_id = 1) const int numLocalSurfaces = 1000;
layout(constant_id = 2) const int numPatches = 1000;
layout(constant_id = 3) const int numSamplesU = 64;
layout(constant_id = 4) const int numSamplesV = 64;

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

layout(set = 0, binding = 4) buffer LocalSurfaceDataBuffer {
	vec4 localSurfaceData[numLocalSurfaces * numSamplesU * numSamplesV * 3];
} lsDataBuffer;

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

// Takes in a u and v coordinate and returns a Sampler object (p,u,v)
// Uses linear interpolation to mix the neighbouring values when the u and
// v coordinates are between sample values
Sampler sampleBuffer(int idx, float u, float v) {
	Sampler res;
	// Step one, transform u,v coordinates [0,1]^2 to buffer coordinates [0,numSamplesU]x[0,numSamplesV]
	float fi = u * (numSamplesU - 1);
	float fj = v * (numSamplesV - 1);
	int i = int(floor(fi));
	int j = int(floor(fj));	
	float fu = fract(fi);
	float fv = fract(fj);
	
	int c00 = j * numSamplesU + i;
	int c10 = j * numSamplesU + i + 1;
	int c01 = (j + 1) * numSamplesU + i;
	int c11 = (j + 1) * numSamplesU + i + 1;
	
	// Step two, linearly interpolate between neighbouring values	
	vec3 p00 = lsDataBuffer.localSurfaceData[idx + c00].xyz;
	vec3 p10 = lsDataBuffer.localSurfaceData[idx + c10].xyz;
	vec3 p01 = lsDataBuffer.localSurfaceData[idx + c01].xyz;
	vec3 p11 = lsDataBuffer.localSurfaceData[idx + c11].xyz;
	res.p = (((1-fu)*p00+fu*p10)*(1-fv)+((1-fu)*p01+fu*p11)*fv);
	
	int offsetU = numSamplesU * numSamplesV;	
	vec3 u00 = lsDataBuffer.localSurfaceData[idx + c00 + offsetU].xyz;
	vec3 u10 = lsDataBuffer.localSurfaceData[idx + c10 + offsetU].xyz;
	vec3 u01 = lsDataBuffer.localSurfaceData[idx + c01 + offsetU].xyz;
	vec3 u11 = lsDataBuffer.localSurfaceData[idx + c11 + offsetU].xyz;
	res.u = (((1-fu)*u00+fu*u10)*(1-fv)+((1-fu)*u01+fu*u11)*fv);
	
	int offsetV = offsetU * 2;	
	vec3 v00 = lsDataBuffer.localSurfaceData[idx + c00 + offsetV].xyz;
	vec3 v10 = lsDataBuffer.localSurfaceData[idx + c10 + offsetV].xyz;
	vec3 v01 = lsDataBuffer.localSurfaceData[idx + c01 + offsetV].xyz;
	vec3 v11 = lsDataBuffer.localSurfaceData[idx + c11 + offsetV].xyz;
	res.v = (((1-fu)*v00+fu*v10)*(1-fv)+((1-fu)*v01+fu*v11)*fv);
	
	// Return the sampled values
	return res;
}

Sampler evaluateBiquadraticBezier(LocalSurfaceInfo lsInfo, float u, float v)
{	
	BoundaryInfo bi = boundaryBuffer.boundaries[lsInfo.boundaryIndex];
	
	Sampler s = sampleBuffer(int(lsInfo.controlPointIndex * numSamplesU * numSamplesV * 3), mix(bi.us, bi.ue, u), mix(bi.vs, bi.ve, v));
			 
	mat4 matrix = matrixBuffer.matrices[lsInfo.matrixIndex];
	
	s.p  = vec3(matrix * vec4(s.p, 1.0f));
	// Should probbly do this! Or make another matrix buffer for normal matrices
	//mat4 normalMat = transpose(inverse(matrix));
	s.u = vec3(matrix * vec4(s.u, 0.0f));
	s.v = vec3(matrix * vec4(s.v, 0.0f));
	
	return s;
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