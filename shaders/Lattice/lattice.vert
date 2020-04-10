#version 450

struct LocalSurfaceInfo
{
	uint controlPointIndex;
	uint controlPointCount;
	uint matrixIndex;
	uint boundaryIndex;
};

layout(location = 0) in uvec4 inLSInfo;
layout(location = 1) in vec3 inColor;

layout(location = 0) out LocalSurfaceInfo vLSInfo;
layout(location = 4) out vec3 vColor;

void main()
{
	LocalSurfaceInfo lsInfo;
	lsInfo.controlPointIndex = inLSInfo.x;
	lsInfo.controlPointCount = inLSInfo.y;
	lsInfo.matrixIndex = inLSInfo.z;
	lsInfo.boundaryIndex = inLSInfo.w;
	
    vLSInfo = lsInfo;
	vColor = inColor;
}