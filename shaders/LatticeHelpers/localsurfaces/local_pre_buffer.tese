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

// Struct to make it easier to work with pre evaluated local surfaces
// Data is layed out in a row by row manner
// struct LocalSurfaceData {
	// vec4 positions[numSamplesU * numSamplesV];
	// vec4 partialU[numSamplesU * numSamplesV];
	// vec4 partialV[numSamplesU * numSamplesV];
// };

layout(set = 0, binding = 4) buffer LocalSurfaceDataBuffer {
	vec4 localSurfaceData[numLocalSurfaces * numSamplesU * numSamplesV * 3];
} lsDataBuffer;

// Takes in a u and v coordinate and returns a Sampler object (p,u,v)
// Uses linear interpolation to mix the neighbouring values when the u and
// v coordinates are between sample values
Sampler sampleBuffer(int idx, float u, float v) {
	Sampler res;
	// Step one, transform u,v coordinates [0,1]^2 to buffer coordinates [0,numSamplesU]x[0,numSamplesV]
	float fi = u * (numSamplesU - 1);
	float fj = v * (numSamplesV - 1);
	int i = int(floor(fi));
	int i_1 = int(ceil(fi));
	int j = int(floor(fj));
	int j_1 = int(ceil(fj));
	
	// Dont interpolate if the values are close to being integers?
	
	float factorU = fi - i;
	float factorV = fj - j;
	
	int c00 = j * numSamplesU + i;
	int c10 = j * numSamplesU + i_1;
	int c01 = j_1 * numSamplesU + i;
	int c11 = j_1 * numSamplesU + i_1;
	
	// Step two, linearly interpolate between neighbouring values
	vec4 p0 = mix(lsDataBuffer.localSurfaceData[idx + c00], lsDataBuffer.localSurfaceData[idx + c10], factorU);
	vec4 p1 = mix(lsDataBuffer.localSurfaceData[idx + c01], lsDataBuffer.localSurfaceData[idx + c11], factorU);
	res.p = mix(p0, p1, factorV).xyz;
	
	int offsetU = numSamplesU * numSamplesV;
	vec4 u0 = mix(lsDataBuffer.localSurfaceData[idx + c00 + offsetU], lsDataBuffer.localSurfaceData[idx + c10 + offsetU], factorU);
	vec4 u1 = mix(lsDataBuffer.localSurfaceData[idx + c01 + offsetU], lsDataBuffer.localSurfaceData[idx + c11 + offsetU], factorU);
	res.u = mix(u0, u1, factorV).xyz;
	
	int offsetV = offsetU * 2;
	vec4 v0 = mix(lsDataBuffer.localSurfaceData[idx + c00 + offsetV], lsDataBuffer.localSurfaceData[idx + c10 + offsetV], factorU);
	vec4 v1 = mix(lsDataBuffer.localSurfaceData[idx + c01 + offsetV], lsDataBuffer.localSurfaceData[idx + c11 + offsetV], factorU);
	res.v = mix(v0, v1, factorV).xyz;
	
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