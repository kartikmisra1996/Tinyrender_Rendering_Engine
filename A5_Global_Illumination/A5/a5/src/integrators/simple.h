/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Simple direct illumination integrator.
 */
struct SimpleIntegrator : Integrator {
    explicit SimpleIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f Li(0.f);
		// TODO: Implement this
		SurfaceInteraction Surfaceinfo;
		SurfaceInteraction secSurface;

		v3f position = scene.getFirstLightPosition();
		v3f intensity = scene.getFirstLightIntensity();

		if (scene.bvh->intersect(ray, Surfaceinfo)) {

			v3f rayDir = glm::normalize(position - Surfaceinfo.p);
			Ray ray = Ray(Surfaceinfo.p, rayDir, 1e-8f, glm::distance(Surfaceinfo.p, position));

			if (!(scene.bvh->intersect(ray, secSurface))) {
				Surfaceinfo.wi = (position - Surfaceinfo.p);

				float R = glm::distance(Surfaceinfo.p, position);
				const BSDF *bsdf = getBSDF(Surfaceinfo);

				Surfaceinfo.wi = glm::normalize(Surfaceinfo.frameNs.toLocal(Surfaceinfo.wi));

				float denom = std::pow(R, 2);

				Li = v3f((intensity / denom) * bsdf->eval(Surfaceinfo));
			}
			else {
				Li = v3f(0.f);
			}
		}
		else {
			Li = v3f(0.f);
		}
		return Li;
    }
};

TR_NAMESPACE_END