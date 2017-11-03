#version 150 core

out vec4 Target0;

layout (std140)
uniform Locals {
        vec4 u_color;
	mat4 u_model;
	mat4 u_view;
	mat4 u_proj;
};

void main() {
        Target0 = u_color;
}
