#version 450

layout(location = 0) in vec3 teColor;

layout(location = 0) out vec4 FragColor;


void main()
{
    FragColor = vec4(teColor, 0.4);
}