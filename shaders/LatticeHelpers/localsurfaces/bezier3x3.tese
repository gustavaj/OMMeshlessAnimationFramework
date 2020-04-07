#version 450

layout(quads, equal_spacing, ccw) in;

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

layout(location = 0) in LocalSurface tcSurfs[];

layout(location = 0) out vec3 teColor;

layout(constant_id = 0) const int numLocalSurfaces = 9;

layout(set = 0, binding = 0) uniform UBO
{
	int tessInnerLevel;
	int tessOuterLevel;
	int bFunctionIndex;
	mat4 projection;
	mat4 modelview;
} ubo;

mat4 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

mat3 evaluateBiquadraticBezier(LocalSurface surf, float u, float v)
{	
	float bu0 = (1-u) * (1-u);
	float bu1 = 2 * u * (1-u);
	float bu2 = u * u;
	
	float dbu0 = -2 * (1-u);
	float dbu1 = 2 - 4 * u;
	float dbu2 = 2 * u;
	
	float bv0 = (1-v) * (1-v);
	float bv1 = 2 * v * (1-v);
	float bv2 = v * v;
	
	float dbv0 = -2 * (1-v);
	float dbv1 = 2 - 4 * v;
	float dbv2 = 2 * v;
	
	vec3 pos = bv0 * (bu0 * surf.p00 + bu1 * surf.p10 + bu2 * surf.p20)
			 + bv1 * (bu0 * surf.p01 + bu1 * surf.p11 + bu2 * surf.p21)
			 + bv2 * (bu0 * surf.p02 + bu1 * surf.p12 + bu2 * surf.p22);
			 
	vec3 dpdu = bv0 * (dbu0 * surf.p00 + dbu1 * surf.p10 + dbu2 * surf.p20)
			  + bv1 * (dbu0 * surf.p01 + dbu1 * surf.p11 + dbu2 * surf.p21)
			  + bv2 * (dbu0 * surf.p02 + dbu1 * surf.p12 + dbu2 * surf.p22);
	
	vec3 dpdv = dbv0 * (bu0 * surf.p00 + bu1 * surf.p10 + bu2 * surf.p20)
			  + dbv1 * (bu0 * surf.p01 + bu1 * surf.p11 + bu2 * surf.p21)
			  + dbv2 * (bu0 * surf.p02 + bu1 * surf.p12 + bu2 * surf.p22);
			 
	mat4 trans = mat4(1.0f, 		0.0f,         0.0f,    		0.0f,
					  0.0f,	   		1.0f, 		  0.0f,    		0.0f,
					  0.0f,    		0.0f,    	  1.0f, 		0.0f,
					  surf.trans.x, surf.trans.y, surf.trans.z, 1.0f);			 
	
	mat4 rot = mat4(1.0f);
	if(surf.rot.z != 0) rot = rot * rotationMatrix(vec3(0, 0, 1), radians(surf.rot.z));
	if(surf.rot.y != 0) rot = rot * rotationMatrix(vec3(0, 1, 0), radians(surf.rot.y));
	if(surf.rot.x != 0) rot = rot * rotationMatrix(vec3(1, 0, 0), radians(surf.rot.x));
			 
	mat4 scale = mat4(surf.scale.x, 0.0f,         0.0f,    		0.0f,
					  0.0f,	   		surf.scale.y, 0.0f,    		0.0f,
					  0.0f,    		0.0f,    	  surf.scale.z, 0.0f,
					  0.0f,         0.0f,         0.0f, 	    1.0f);
					  
	mat4 normalScale = mat4(1/surf.scale.x, 0.0f,         0.0f,    		0.0f,
							0.0f,	   		1/surf.scale.y, 0.0f,    		0.0f,
							0.0f,    		0.0f,    	  1/surf.scale.z, 0.0f,
							0.0f,         0.0f,         0.0f, 	    1.0f);
					  
	mat4 model = trans * rot * scale;
	
	pos = vec3(model * vec4(pos, 1.0f));
	dpdu = vec3(model * vec4(dpdu, 0.0f));
	dpdv = vec3(model * vec4(dpdv, 0.0f));
	
	return mat3(pos.xyz, dpdu.xyz, dpdv.xyz);
}

void main()
{
	float u = gl_TessCoord.x;
	float v = gl_TessCoord.y;
	
	mat3 p = evaluateBiquadraticBezier(tcSurfs[0], u, v);
		
	gl_Position = ubo.projection * ubo.modelview * vec4(p[0], 1.0f);
	
	float fac = gl_PrimitiveID / float(numLocalSurfaces);
	
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