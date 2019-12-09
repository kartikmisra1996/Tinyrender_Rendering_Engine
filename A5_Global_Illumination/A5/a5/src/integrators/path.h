/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Path tracer integrator
 */
struct PathTracerIntegrator : Integrator {
    explicit PathTracerIntegrator(const Scene& scene) : Integrator(scene) {
        m_isExplicit = scene.config.integratorSettings.pt.isExplicit;
        m_maxDepth = scene.config.integratorSettings.pt.maxDepth;
        m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
        m_rrProb = scene.config.integratorSettings.pt.rrProb;
    }

	
	v3f renderImplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {		

		v3f Li = v3f(1.f);
		v2f thisSample = sampler.next2D();
		SurfaceInteraction intersect = hit;

		int depth = -1;
		while (depth <= m_maxDepth) {

			depth++;
			v3f emission = getEmission(hit);
			float dotProd = glm::dot(v3f(0, 0, 1), hit.wo);
			bool cond = emission != v3f(0.f) && (dotProd > 0);

			if (cond) {
				Li = Li * emission;
				return Li;
			}

				float pdf;
				float* pdfPointer = &pdf;
				const BSDF* bsdf = getBSDF(hit);

				v3f contVec = bsdf->sample(intersect, thisSample, pdfPointer);

				Li *= contVec; 
				v3f wiW = intersect.frameNs.toWorld(hit.wi);

				wiW = glm::normalize(wiW);
				Ray ray = Ray(intersect.p, wiW); 

				bool noHit = !scene.bvh->intersect(ray, intersect);
				if (noHit) {
					return v3f(0.f);
				}
		}
		return v3f(0.f);
	}

	v3f computDL(SurfaceInteraction& hit, Sampler& sampler) const {

		float emPDF;
		float areaPDF;

		size_t selectedEmitter = selectEmitter(sampler.next(), emPDF);
		const Emitter& emitter = getEmitterByID(selectedEmitter);

		v3f normal;
		v3f position;

		sampleEmitterPosition(sampler, emitter, normal, position, areaPDF); 

		v3f DLightingValue = v3f(0.f);
		v3f wiW = glm::normalize(position - hit.p); 

		hit.wi = hit.frameNs.toLocal(wiW);
		float cosT = glm::dot(-wiW, normal); 

		if (cosT > 0) {
			SurfaceInteraction shadowSI;
			Ray ray = Ray(hit.p, wiW);

			bool intersect = scene.bvh->intersect(ray, shadowSI);
			if (intersect) {

				v3f emission = getEmission(shadowSI);

				bool isEmitting = emission != v3f(0.f);
				if (isEmitting) {
					float distance = glm::distance2(position, hit.p);
					float correction = cosT / distance;

					DLightingValue = emission * getBSDF(hit)->eval(hit) ;
					DLightingValue /= (emPDF*areaPDF)*correction;
				}

			}
			else {
				DLightingValue = v3f(0.f);
			}
		}
		else {
			DLightingValue = v3f(0.f);
		}
		return DLightingValue;
	}

	v3f computeIL(SurfaceInteraction& hit, Sampler& sampler, int depth) const {

		int recD = depth + 1;
		float rusPDF = 1.f;
		if (m_maxDepth == -1) {

			if (recD >= m_rrDepth) {
				if (sampler.next() > m_rrProb) {
					return v3f(0.f);
				}
				else {
					rusPDF = m_rrProb;
				}
			}

		}
		else {

			if (recD >= m_maxDepth) {
				return v3f(0.f);
			}

		}

		v3f bsdfFinal(1.f);

		v3f emission = v3f(1.f);
		SurfaceInteraction ray;

		bool noEm = emission != v3f(0.f);
		while (noEm) {
			float pdf;
			bsdfFinal = getBSDF(hit)->sample(hit, sampler.next2D(), &pdf);
			v3f wiW = hit.frameNs.toWorld(hit.wi);

			wiW = glm::normalize(wiW);
			Ray nextRay = Ray(hit.p, wiW);

			bool intersect = scene.bvh->intersect(nextRay, ray);
			if (intersect) {
				emission = getEmission(ray);
			}
			else {
				return v3f(0.f);
			}
		}

		v3f totalLighting = (computDL(ray, sampler) + computeIL(ray, sampler, recD));
		return 1 / rusPDF * bsdfFinal * totalLighting;

	}

	v3f renderExplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		v3f Li(0.f);


		v3f emission = getEmission(hit);
		if (emission != v3f(0.f)) {
			return emission;
		}

		v3f dir;
		v3f indir;

		if (m_maxDepth == 0) {
			dir = v3f(0.f);

		}
		else {

			dir = computDL(hit, sampler);

		}
		int depth = 0;

		indir = computeIL(hit, sampler, depth);

		Li = dir + indir;
		return Li;
	}

    v3f render(const Ray& ray, Sampler& sampler) const override {
        Ray r = ray;
        SurfaceInteraction hit;

        if (scene.bvh->intersect(r, hit)) {
            if (m_isExplicit)
                return this->renderExplicit(ray, sampler, hit);
            else
                return this->renderImplicit(ray, sampler, hit);
        }
        return v3f(0.0);
    }

    int m_maxDepth;     // Maximum number of bounces
    int m_rrDepth;      // When to start Russian roulette
    float m_rrProb;     // Russian roulette probability
    bool m_isExplicit;  // Implicit or explicit
};

TR_NAMESPACE_END
