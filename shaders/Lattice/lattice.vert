#version 450

layout(location = 0) in vec3 p00;
layout(location = 1) in vec3 p10;
layout(location = 2) in vec3 p20;
layout(location = 3) in vec3 p01;
layout(location = 4) in vec3 p11;
layout(location = 5) in vec3 p21;
layout(location = 6) in vec3 p02;
layout(location = 7) in vec3 p12;
layout(location = 8) in vec3 p22;
layout(location = 9) in vec3 trans;
layout(location = 10) in vec3 rot;
layout(location = 11) in vec3 scale;
layout(location = 12) in vec2 boundaryUV;

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

layout(location = 0) out LocalSurface vSurf;

void main()
{
    LocalSurface surf;
	
	surf.p00 = p00;
	surf.p10 = p10;
	surf.p20 = p20;
	surf.p01 = p01;
	surf.p11 = p11;
	surf.p21 = p21;
	surf.p02 = p02;
	surf.p12 = p12;
	surf.p22 = p22;
	surf.trans = trans;
	surf.rot = rot;
	surf.scale = scale;
	surf.boundaryUV = boundaryUV;
	
	vSurf = surf;
}