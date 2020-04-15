#version 450

layout(location = 0) in vec3 inPos;
layout(location = 1) in uvec3 inCol;

layout(location = 0) out uvec3 vCol;

layout(set = 0, binding = 0) uniform UBO
{
	int tessInner;
	int tessOuter;
	int bFunctionIndex;
	mat4 projection;
	mat4 modelview;
	mat4 normal;
} ubo;

void main()
{
    gl_Position = ubo.projection * ubo.modelview * vec4(inPos, 1.0f);
	gl_PointSize = 10.0f;
	vCol = inCol;
}