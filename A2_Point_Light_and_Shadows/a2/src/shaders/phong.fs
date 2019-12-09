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
uniform float exponent;
uniform vec3 albedo;
uniform vec3 rho_d;
uniform vec3 rho_s;
float invPi = 1/pi;

in vec3 vNormal;
in vec3 vPos;
out vec3 color;

void main() {
  vec3 lightDirection = (lightPos - vPos);
  vec3 camDirection = (camPos - vPos);
  vec3 reflection = reflect((lightDirection), (vNormal));

  float specularAngle = dot(normalize(reflection), normalize(camDirection));
  vec3 phongBSDFEval = (rho_d * invPi) + rho_s * (exponent + 2) / (2 * pi) * pow(cos(specularAngle), exponent);
  vec3 colorDiffuse = (phongBSDFEval) * dot(normalize(lightDirection), normalize(vNormal));
  
  color = colorDiffuse * (lightIntensity / pow(length(lightDirection), 2));
}
