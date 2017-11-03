#version 150 core

in vec3 v_pos;

layout (std140)
uniform Locals {
        vec4 u_color;
	mat4 u_model;
	mat4 u_view;
	mat4 u_proj;
};

void main() {
        gl_Position = u_proj * u_view * u_model * vec4(v_pos, 1.0);
}
