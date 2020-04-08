#version 450

struct LocalSurfaceInfo
{
	uint controlPointIndex;
	uint controlPointCount;
	uint matrixIndex;
	uint boundaryIndex;
};

layout(location = 0) in uvec4 inLSInfo;

layout(location = 0) out LocalSurfaceInfo vLSInfo;

void main()
{
	LocalSurfaceInfo lsInfo;
	lsInfo.controlPointIndex = inLSInfo.x;
	lsInfo.controlPointCount = inLSInfo.y;
	lsInfo.matrixIndex = inLSInfo.z;
	lsInfo.boundaryIndex = inLSInfo.w;
	
    vLSInfo = lsInfo;
}