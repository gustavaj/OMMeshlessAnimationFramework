#version 450

struct LocalSurfaceInfo
{
	uint controlPointIndex;
	uint controlPointCount;
	uint matrixIndex;
	uint boundaryIndex;
};

layout(vertices = 1) out;

layout(location = 0) in LocalSurfaceInfo[] vLSInfo;

layout(location = 0) out LocalSurfaceInfo[] tcLSInfo;

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

void main()
{
    tcLSInfo[gl_InvocationID] = vLSInfo[gl_InvocationID];

	if(gl_InvocationID == 0) {
		gl_TessLevelInner[0] = latticeUbo.tessInnerLevel;
		gl_TessLevelInner[1] = latticeUbo.tessInnerLevel;
		gl_TessLevelOuter[0] = latticeUbo.tessOuterLevel;
		gl_TessLevelOuter[1] = latticeUbo.tessOuterLevel;
		gl_TessLevelOuter[2] = latticeUbo.tessOuterLevel;
		gl_TessLevelOuter[3] = latticeUbo.tessOuterLevel;
	}
}