#version 450

layout(location = 0) in vec3 teColor;
layout(location = 1) in vec3 teNormal;

layout(location = 0) out vec4 FragColor;

void main()
{
	vec3 light = normalize(vec3(0.5f, 0.5f, 1.0f));
	vec3 ambient = vec3(0.3f);
	vec3 diffuse = max(dot(light, teNormal), 0.0f) * vec3(0.8f);
	
    FragColor = vec4((ambient + diffuse) * teColor, 1.0);
}