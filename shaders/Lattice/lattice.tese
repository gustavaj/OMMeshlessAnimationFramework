#version 450

layout(quads, equal_spacing, ccw) in;

const int B1Poly = 0;
const int B2Poly = 1;
const int Lerbs = 2;

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
layout(location = 1) out vec3 teNormal;
layout(location = 2) out vec3 tePosition;

layout(set = 0, binding = 0) uniform UBO
{
	int tessInnerLevel;
	int tessOuterLevel;
	int bFunctionIndex;
	mat4 projection;
	mat4 modelview;
	mat4 localSurfaceMatrices[4];
} ubo;

layout (constant_id = 0) const int degree = 2;
layout (constant_id = 1) const int totalPatches = 9;

layout( push_constant ) uniform PushConstants 
{
	int isLeftBoundary;
	int isRightBoundary;
	int isTopBoundary;
	int isBottomBoundary;
	float idx;
} push;

int indexI(float u) 
{
	return int(max(1, ceil(u * (degree - 1))));
}

int indexJ(float v) 
{
	return int(max(1, ceil(v * (degree - 1))));
}

float bFunction(float w)
{
	if(w < 1e-5) return 0.0f;
	if(1 - w < 1e-5) return 1.0f;
	switch(ubo.bFunctionIndex) {
		case B1Poly: { return (3*pow(w, 2)-2*pow(w,3)); }
		case B2Poly: { return 6 * pow(w, 5) - 15 * pow(w, 4) + 10 * pow(w, 3); }
		case Lerbs: { return 1.0 / (1.0 + exp(1.0 / w - 1.0 / (1.0 - w))); }
	}	
}

float bDerivative(float w)
{
	if(w < 1e-5) return 0.0f;
	if(1 - w < 1e-5) return 1.0f;
	switch(ubo.bFunctionIndex) {
		case B1Poly: { return 6*w-6*pow(w,2); }
		case B2Poly: { return 30 * pow(w, 4) - 60 * pow(w, 3) + 30 * pow(w, 2); }
		case Lerbs: {
			return (2 * exp( - ( (1.0) / ((w-1.0)*w) ) ) *
				(pow(w,2.0) - w + 0.5)) /
				(pow( exp( -(1.0/(w-1.0)) ) + exp(1.0/w), 2.0) *
				pow(w-1.0,2.0) * pow(w,2.0)); }
	}	
}

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

float getLocalU(LocalSurface surf, int i, float u, bool isLeft)
{
	if(isLeft) {
		if(i == 1) {
			if(push.isLeftBoundary == 1) { return u; }
			else { return mix(0.5f + (surf.boundaryUV[0] - 0.5f) / 2, 1.0, u); }
		}
		else { return mix(0.5, 1.0, u); }
	}
	else {
		if(i == (degree - 1)) {
			if(push.isRightBoundary == 1) { return u; }
			else { return mix(0.0, 0.5f + (surf.boundaryUV[0] - 0.5f) / 2, u); }
		}
		else { return mix(0.0, 0.5, u); }
	}
}

float getLocalV(LocalSurface surf, int j, float v, bool isTop)
{
	if(isTop) {
		if(j == 1) {
			if(push.isTopBoundary == 1) { return v; }
			else { return mix(0.5f + (surf.boundaryUV[1] - 0.5f) / 2, 1.0, v); }
		}
		else { return mix(0.5, 1.0, v); }
	}
	else {
		if(j == (degree - 1)) {
			if(push.isBottomBoundary == 1) { return v; }
			else { return mix(0.0, 0.5f + (surf.boundaryUV[1] - 0.5f) / 2, v); }
		}
		else { return mix(0.0, 0.5, v); }
	}
}

vec3 bezBasis(float t) {
	return vec3(pow(1-t, 2), 2*t*(1-t), pow(t, 2));
}

vec3 bezBasisDer(float t) {
	return vec3(2*t-2, 2-4*t, 2*t);
}

mat3 evaluateBiquadraticBezier(LocalSurface surf, float u, float v, int index)
{
	vec3 bu = bezBasis(u);
	vec3 bud = bezBasisDer(u);
	vec3 bv = bezBasis(v);
	vec3 bvd = bezBasisDer(v);
			  
	vec3 pos = 
		surf.p00 * bu[0] * bv[0] + surf.p01 * bu[0] * bv[1] + surf.p02 * bu[0] * bv[2] +
		surf.p10 * bu[1] * bv[0] + surf.p11 * bu[1] * bv[1] + surf.p12 * bu[1] * bv[2] +
		surf.p20 * bu[2] * bv[0] + surf.p21 * bu[2] * bv[1] + surf.p22 * bu[2] * bv[2];
			  
	vec3 dpdu = 
		surf.p00 * bud[0] * bv[0] + surf.p01 * bud[0] * bv[1] + surf.p02 * bud[0] * bv[2] +
		surf.p10 * bud[1] * bv[0] + surf.p11 * bud[1] * bv[1] + surf.p12 * bud[1] * bv[2] +
		surf.p20 * bud[2] * bv[0] + surf.p21 * bud[2] * bv[1] + surf.p22 * bud[2] * bv[2];
			  
	vec3 dpdv = 
		surf.p00 * bu[0] * bvd[0] + surf.p01 * bu[0] * bvd[1] + surf.p02 * bu[0] * bvd[2] +
		surf.p10 * bu[1] * bvd[0] + surf.p11 * bu[1] * bvd[1] + surf.p12 * bu[1] * bvd[2] +
		surf.p20 * bu[2] * bvd[0] + surf.p21 * bu[2] * bvd[1] + surf.p22 * bu[2] * bvd[2];
	
	pos  = vec3(ubo.localSurfaceMatrices[index] * vec4(pos, 1.0f));
	dpdu = vec3(ubo.localSurfaceMatrices[index] * vec4(dpdu, 0.0f));
	dpdv = vec3(ubo.localSurfaceMatrices[index] * vec4(dpdv, 0.0f));
	
	return mat3(pos.xyz, dpdu.xyz, dpdv.xyz);
}

void main()
{
	float u = gl_TessCoord.x;
	float v = gl_TessCoord.y;
	
	int i = indexI(u);
	int j = indexJ(v);	
	
	// Scale U/V, "w()"
	float stepUV = 1.0f / float(degree - 1);
	float scaledU = (u - (i - 1) * stepUV) / stepUV;
	float scaledV = (v - (j - 1) * stepUV) / stepUV;
	
	LocalSurface local00 = tcSurfs[(i - 1) + (j - 1) * degree];
	LocalSurface local10 = tcSurfs[i + (j - 1) * degree];
	LocalSurface local01 = tcSurfs[(i - 1) + j * degree];
	LocalSurface local11 = tcSurfs[i + j * degree];
	
	// Evaluate local surfaces
	mat3 s00 = evaluateBiquadraticBezier(local00, getLocalU(local00, i, scaledU, true),  getLocalV(local00, j, scaledV, true),  (i - 1) + (j - 1) * degree);
	mat3 s10 = evaluateBiquadraticBezier(local10, getLocalU(local10, i, scaledU, false), getLocalV(local10, j, scaledV, true),   i + (j - 1) * degree);
	mat3 s01 = evaluateBiquadraticBezier(local01, getLocalU(local01, i, scaledU, true),  getLocalV(local01, j, scaledV, false), (i - 1) + j * degree);
	mat3 s11 = evaluateBiquadraticBezier(local11, getLocalU(local11, i, scaledU, false), getLocalV(local11, j, scaledV, false),  i + j * degree);
	
	float bu = bFunction(scaledU);
	float bv = bFunction(scaledV);
	float bud = bDerivative(scaledU);
	float bvd = bDerivative(scaledV);
	
	// If the u,v value is on the boundary, interpolate only those points to ensure smoothness over the edges.
	vec3 pos, dpdu, dpdv;
	if(scaledU < 1e-5) {
		if(scaledV < 1e-5) {
			pos = s00[0]; dpdu = s00[1]; dpdv = s00[2]; // Interpolate top left corner
		}
		else if(1.0f - scaledV < 1e-5) {
			pos = s01[0]; dpdu = s01[1]; dpdv = s01[2]; // Interpolate bottom left corner
		}
		else {
			// Interpolate left edge
			pos  = s00[0] + bv * (s01[0] - s00[0]);
			dpdu = s00[1] + bv * (s01[1] - s00[1]);
			dpdv = s00[2] + bv * (s01[2] - s00[2]) + bvd * (s01[0] - s00[0]);			
		}
	}
	else if(1.0f - scaledU < 1e-5) {
		if(scaledV < 1e-5) {
			pos = s10[0]; dpdu = s10[1]; dpdv = s10[2]; // Interpolate top right corner
		}
		else if(1.0f - scaledV < 1e-5) {
			pos = s11[0]; dpdu = s11[1]; dpdv = s11[2]; // Interpolate bottom right corner
		}
		else {
			// Interpolate right edge
			pos  = s10[0] + bv * (s11[0] - s10[0]);
			dpdu = s10[1] + bv * (s11[1] - s10[1]);
			dpdv = s10[2] + bv * (s11[2] - s10[2]) + bvd * (s11[0] - s10[0]);
		}
	}
	else if(scaledV < 1e-5) {
		// Interpolate top edge
		pos  = s00[0] + bu * (s10[0] - s00[0]);
		dpdu = s00[1] + bu * (s10[1] - s00[1]) + bud * (s10[0] - s00[0]);
		dpdv = s00[2] + bu * (s10[2] - s00[2]);
	}
	else if(1.0f - scaledV < 1e-5) {
		// Interpolate bottom edge
		pos  = s01[0] + bu * (s11[0] - s01[0]);
		dpdu = s01[1] + bu * (s11[1] - s01[1]) + bud * (s11[0] - s01[0]);
		dpdv = s01[2] + bu * (s11[2] - s01[2]);
	}
	else {
		// Blend for final position
		pos = (((1 - bu) * s00[0] + bu * s10[0]) * (1 - bv)
			+  ((1 - bu) * s01[0] + bu * s11[0]) *      bv);

		dpdu = (((-bud)   * s00[0] + bud * s10[0]) * (1 - bv)
		     +  ((-bud)   * s01[0] + bud * s11[0]) *      bv)
			 + (((1 - bu) * s00[1] + bu  * s10[1]) * (1 - bv)
			 +  ((1 - bu) * s01[1] + bu  * s11[1]) *      bv);
			 
		dpdv = (((1 - bu) * s00[0] + bu * s10[0]) *   (-bvd)
			 +  ((1 - bu) * s01[0] + bu * s11[0]) *     bvd)
			 + (((1 - bu) * s00[2] + bu * s10[2]) * (1 - bv)
			 +  ((1 - bu) * s01[2] + bu * s11[2]) *      bv);
	}
				 
	
			  
	gl_Position = ubo.projection * ubo.modelview * vec4(pos, 1.0f);
	tePosition = vec3(ubo.modelview * vec4(pos, 1.0f));
	teNormal = normalize(cross(dpdu, dpdv));
	
	float fac = push.idx / float(totalPatches);
	if(fac < 0.5f) {
		fac = fac * 2.0f;
		teColor = vec3(1.0f - fac, fac, 0.0f);
	}
	else if (fac < 1.0f){
		fac = (fac - 0.5f) * 2.0f;
		teColor = vec3(1.0f - fac, 1.0f - fac, fac);
	}
}