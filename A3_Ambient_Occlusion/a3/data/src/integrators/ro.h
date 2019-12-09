/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include <random>

TR_NAMESPACE_BEGIN

/**
 * Reflection occlusion integrator
 */
struct ROIntegrator : Integrator {

    float m_exponent;

    explicit ROIntegrator(const Scene& scene) : Integrator(scene) {
        m_exponent = scene.config.integratorSettings.ro.exponent;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f Li(0.f);

		SurfaceInteraction surfaceInfo;

		bool intersect1 = scene.bvh->intersect(ray, surfaceInfo);
		//check if the ray intersects
		if (intersect1) {

			// No for loop reaquired as A1 bonus samples each pixel 16 times already
			SurfaceInteraction secondInteraction;
			p2f p(sampler.next2D());

			v3f wi = Warp::squareToPhongLobe(p, m_exponent);
			v3f wr = reflect(surfaceInfo.wo);
			v3f wrW = glm::normalize(surfaceInfo.frameNs.toWorld(wr));

			Frame rFrame = Frame(wrW);

			wi = glm::toMat4(glm::quat(v3f(0, 0, 1), reflect(surfaceInfo.wo))) * v4f(wi, 1);
			wi = glm::normalize(wi);

			v3f wiW = glm::normalize(rFrame.toWorld(wi));

			// casting secondary ray from intersect point 1 
			Ray ray = Ray(surfaceInfo.p + Epsilon, wiW, Epsilon);

			bool intersect2 = scene.bvh->intersect(ray, secondInteraction);
			// check if the bouncing ray intersects
			if (intersect2) {
				Li = Li + v3f(0.f);
			}
			else {
				float cosT = glm::dot(surfaceInfo.frameNs.n, surfaceInfo.frameNs.toWorld(wi));
				float cosA = fmax(wi.z, 0.f);

				//compute Li unless cosT is negative
				if (cosT > 0) {
					float inter = Warp::squareToPhongLobePdf(wi, m_exponent);
					Li = Li + cosT * 1 * (m_exponent + 2) * INV_TWOPI * max(0.0f, pow(cosA, m_exponent)) / inter;
				}
			}
		}
		return Li;
    }
};

TR_NAMESPACE_END