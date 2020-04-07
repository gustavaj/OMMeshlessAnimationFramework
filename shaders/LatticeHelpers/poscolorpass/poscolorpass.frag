#version 450

layout(location = 0) flat in uvec3 vColor;

layout(location = 0) out vec4 FragColor;

void main()
{
	if(vColor.r == 255 && vColor.g == 255 && vColor.b == 255) {
		discard;
	}
	FragColor = vec4(vColor.r / 255.0f, vColor.g / 255.0f, vColor.b / 255.0f, 1.0);
}