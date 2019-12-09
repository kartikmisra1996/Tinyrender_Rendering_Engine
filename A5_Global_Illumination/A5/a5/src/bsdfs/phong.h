/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

	v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);
		v3f perfectReflection = PhongBSDF::reflect(i.wo);
		float specularAngle = glm::dot(glm::normalize(i.wi), glm::normalize(perfectReflection));

		if (specularAngle < 0) {
			specularAngle = 0;
		}


		if ((Frame::cosTheta(glm::normalize(i.wo)) > 0.f && Frame::cosTheta(glm::normalize(i.wi))) > 0.f) {

			v3f phongReflectance = (diffuseReflectance->eval(worldData, i) * INV_PI)  //diffuse part
				+ (specularReflectance->eval(worldData, i)                     // >
					* ((exponent->eval(worldData, i) + 2) / (2 * M_PI))              // >> specular part
					* std::pow(specularAngle, exponent->eval(worldData, i)));     // >

			val = scale * phongReflectance * (glm::dot(i.wi, glm::normalize(i.frameNs.toLocal(i.frameNs.n)))); //multiply by cosine factor
		}
		else {
			val = v3f(0.f);

		}

		return val;
	}

	float pdf(const SurfaceInteraction& i) const override {
		float pdf = 0.f;

		v3f refl = reflect(i.wo);
		v3f v = glm::toMat4(glm::quat(refl, v3f(0, 0, 1))) * v4f(i.wi, 1);

		float expVal = exponent->eval(worldData, i);
		pdf = Warp::squareToPhongLobePdf(v, expVal);
		return pdf;
	}

	v3f sample(SurfaceInteraction& i, const v2f& _sample, float* pdf1) const override {
		v3f val(0.f);

		float expVal = exponent->eval(worldData, i);
		v3f wi = Warp::squareToPhongLobe(_sample, expVal);

		glm::quat temp = glm::quat(v3f(0, 0, 1), reflect(i.wo));
		wi = glm::toMat4(temp) * v4f(wi, 1);
		i.wi = wi;

		*pdf1 = pdf(i);


		if (pdf(i) == 0.f) {
			return v3f(0);
		}

		else {
			val = eval(i) / pdf(i);
			return val;
		}
	}

    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END