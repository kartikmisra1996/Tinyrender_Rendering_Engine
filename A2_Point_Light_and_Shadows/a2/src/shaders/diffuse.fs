/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core

#define pi 3.14159265358979323846

uniform vec3 camPos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform vec3 albedo;


in vec3 vNormal;
in vec3 vPos;
out vec3 color;

void main() {

    vec3 lightDirection = (lightPos - vPos);


    float R = length(lightDirection);
	float result = dot(normalize(lightDirection), normalize(vNormal));
    vec3 brdf = albedo / (pi) * dot(normalize(lightDirection), normalize(vNormal));
	vec3 shading = (lightIntensity/ pow(R,2));
    color = brdf * shading;

}