#version 450 

layout( triangles ) in; 

layout( location = 1 ) in vec3 vert_normal[]; 
layout( location = 2 ) in vec3 vert_position[];

layout(set = 0, binding = 0) uniform UBO
{
	int tessInnerLevel;
	int tessOuterLevel;
	int bFunctionIndex;
	mat4 projection;
	mat4 modelview;
	mat4 normal;
	vec2 windowSize;
} ubo;

layout( line_strip, max_vertices = 6 ) out; 

layout( location = 0 ) out vec4 geom_color; 

void main() { 
  for( int vertex = 0; vertex < 3; ++vertex ) { 
    gl_Position = ubo.projection * vec4(vert_position[vertex], 1);
    geom_color = vec4( 0.2 ); 
    EmitVertex(); 

    gl_Position = ubo.projection * (vec4(vert_position[vertex], 1) + (transpose(inverse(ubo.modelview)) * vec4(vert_normal[vertex] * 1.5f, 0)));
    geom_color = vec4( 0.6 ); 
    EmitVertex(); 

    EndPrimitive(); 
  } 
}