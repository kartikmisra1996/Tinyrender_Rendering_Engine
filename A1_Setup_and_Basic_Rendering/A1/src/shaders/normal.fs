/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core

in vec3 vNormal;
out vec3 color;

void main() {
	color = vNormal;

	//gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}
