/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model + Diffuse
 */
struct MixtureBSDF : BSDF {
    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    MixtureBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
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
		v3f result(0.f);
		v3f wi = i.wi;
		v3f wo = i.wo;

		if (Frame::cosTheta(wi) > 0 && Frame::cosTheta(wo) > 0) {

			v3f diffuse = diffuseReflectance->eval(worldData, i);
			v3f specular = specularReflectance->eval(worldData, i);

			float n = exponent->eval(worldData, i);
			float max = fmax(glm::dot(glm::normalize(reflect(wo)), wi), 0.0);

			max = fmin(1, max);
			v3f vTemp = diffuse * INV_PI + specular ;
			vTemp *= INV_TWOPI * (n + 2)*pow(max, n);
			result = vTemp * Frame::cosTheta(wi) * scale;

		}
		else {
			result = v3f(0.f);
		}


		return result;
	}

	float pdf(const SurfaceInteraction& i) const override {
		float pdf = 0.f;
		v3f v = glm::toMat4(glm::quat(reflect(i.wo), v3f(0, 0, 1))) * v4f(i.wi, 1);

		float spec = fmax(Warp::squareToPhongLobePdf(v, exponent->eval(worldData, i)), 0);
		float diffusePdf = fmax(Warp::squareToCosineHemispherePdf(i.wi), 0);

		pdf = spec * specularSamplingWeight + diffusePdf * (1 - specularSamplingWeight);
		return pdf;
	}

	v3f sample(SurfaceInteraction& i, const v2f& _sample, float* pdfv) const override {
		v3f val(0.f);
		if (_sample.x < specularSamplingWeight) {

			float newSample = _sample.x / specularSamplingWeight;
			v2f sample = v2f(newSample, _sample.y);

			float exp = exponent->eval(worldData, i);
			v3f wi = Warp::squareToPhongLobe(sample, exp);

			v3f ref = reflect(i.wo);
			wi = glm::toMat4(glm::quat(v3f(0, 0, 1), ref)) * v4f(wi, 1);
			i.wi = wi;
		}
		else {

			float x1 = _sample.x - specularSamplingWeight;
			float x2 = 1 - specularSamplingWeight;
			v2f sampleModifed = v2f(x1 / x2, _sample.y);

			i.wi = Warp::squareToCosineHemisphere(sampleModifed);
		}

		*pdfv = pdf(i);
		if (pdf(i) == 0.f) {
			return v3f(0.f);
		}
		val = eval(i) / pdf(i);
		return val;
	}

    std::string toString() const override { return "Mixture"; }
};

TR_NAMESPACE_END