#version 450

struct LocalSurface {
	vec3 p00;
	vec3 p10;
	vec3 p20;
	vec3 p01;
	vec3 p11;
	vec3 p21;
	vec3 p02;
	vec3 p12;
	vec3 p22;
	vec3 trans;
	vec3 rot;
	vec3 scale;
	vec2 boundaryUV;
};

layout(vertices = 25) out;

layout(location = 0) in LocalSurface vSurf[];

layout(location = 0) out LocalSurface tcSurfs[];

layout(set = 0, binding = 0) uniform UBO
{
	int tessInnerLevel;
	int tessOuterLevel;
	int bFunctionIndex;
	mat4 projection;
	mat4 modelview;
	mat4[4] localSurfaceMatrices;
} ubo;

void main()
{
    tcSurfs[gl_InvocationID] = vSurf[gl_InvocationID];

	if(gl_InvocationID == 0) {
		gl_TessLevelInner[0] = ubo.tessInnerLevel;
		gl_TessLevelInner[1] = ubo.tessInnerLevel;
		gl_TessLevelOuter[0] = ubo.tessOuterLevel;
		gl_TessLevelOuter[1] = ubo.tessOuterLevel;
		gl_TessLevelOuter[2] = ubo.tessOuterLevel;
		gl_TessLevelOuter[3] = ubo.tessOuterLevel;
	}
}