/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/* this determines what type of sampling we're using
0 = sphereical
1 = hemispherical
2 = cosine
*/
int type = 2;

/**
 * Ambient occlusion integrator
 */
struct AOIntegrator : Integrator {
    explicit AOIntegrator(const Scene& scene) : Integrator(scene) { }

	v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f Li(0.f);

		SurfaceInteraction sInfo;
		int numDir = 1;

		bool intersect1 = scene.bvh->intersect(ray, sInfo);
		//check if the ray intersects
		if (intersect1) {

			for (int i = 0; i < numDir; i++) {

				SurfaceInteraction shadow;
				v3f sampledDir;
				p2f p(sampler.next2D());

				//choosing between sampling methods
				if (type == 0) {
					sampledDir = Warp::squareToUniformSphere(p);
				}
				else if (type == 1) {
					sampledDir = Warp::squareToUniformHemisphere(p);
				}
				else if (type == 2) {
					sampledDir = Warp::squareToCosineHemisphere(p);
				}

				Ray ray = Ray(sInfo.p + Epsilon, sInfo.frameNs.toWorld(sampledDir), Epsilon, scene.aabb.getBSphere().radius / 2 - Epsilon);

				bool intersect2 = scene.bvh->intersect(ray, sInfo);

				//checking if the reflected ray intersects
				if (intersect2) {
					Li += v3f(0.f);
				}

				else {

					float albedo = 1;
					sInfo.wi = sampledDir;
					float cosTheta = Frame::cosTheta(sampledDir);

					// avoid cosine clamping
					if (cosTheta > 0) {

						v3f param1 = glm::normalize(sInfo.wi);
						v3f param2 = glm::normalize(sInfo.frameNs.toLocal(sInfo.frameNs.n));
						float x = glm::dot(param1, param2);

						//computing Li for the chosen sampling method
						if (type == 0) {
							Li += albedo * x * INV_PI / (numDir * Warp::squareToUniformSpherePdf());
						}
						else if (type == 1) {
							Li += albedo * x * INV_PI / (numDir * Warp::squareToUniformHemispherePdf(sInfo.wi));
						}
						else if (type == 2) {
							Li += albedo * x * INV_PI / (numDir * Warp::squareToCosineHemispherePdf(sInfo.wi));
						}
					}
				}
			}
		}
		return Li;
	}
};

TR_NAMESPACE_END