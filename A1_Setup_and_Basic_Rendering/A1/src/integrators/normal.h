/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/platform.h>
#include <core/integrator.h>
#include <core/accel.h>

TR_NAMESPACE_BEGIN

/**
 * Surface normal integrator.
 */
struct NormalIntegrator : Integrator {
    explicit NormalIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f color;
		v3f normal;
		SurfaceInteraction info;
		if (scene.bvh->intersect(ray, info)) {
			v3f surfInfo = info.frameNs.n;
			normal = glm::abs(surfInfo);
			color = normal;
		}
		else {
			color = v3f(0.f, 0.f, 0.f);
		}
        return color;
    }
};

TR_NAMESPACE_END