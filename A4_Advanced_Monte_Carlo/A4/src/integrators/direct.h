/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
struct DirectIntegrator : Integrator {
    explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
        m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
        m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
        m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
    }

    static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        float f = nf * fPdf, g = ng * gPdf;
        return f / (f + g);
    }

    void sampleSphereByCosineHemisphere(const p2f& sample,
                                        const v3f& n,
                                        const p3f& pShading,
                                        const v3f& emitterCenter,
                                        float emitterRadius,
                                        v3f& wiW,
                                        float& pdf) const {
        // TODO: Implement this
		// Code in renderCosineHemisphere as it was easier to trace when everything was together
    }

    void sampleSphereByArea(const p2f& sample,
                            const p3f& pShading,
                            const v3f& emitterCenter,
                            float emitterRadius,
                            v3f& pos,
                            v3f& ne,
                            v3f& wiW,
                            float& pdf) const {
        // TODO: Implement this
		// Code in renderArea as it was easier to trace when everything was together
    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        // TODO: Implement this
		// Code in renderSolidAngle and renderMIS as it was easier to trace when everything was together
    }

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
		int noOfSamples = scene.config.integratorSettings.di.emitterSamples;
		SurfaceInteraction sInteraction;

		bool isIntersecting = scene.bvh->intersect(ray, sInteraction);
		if (scene.bvh->intersect(ray, sInteraction)) {
			v3f emission = getEmission(sInteraction);

			bool onEmitter = emission != v3f(0.0f);
			if (onEmitter) {
				return emission;
			}
		}

		for (int i = 0; i < noOfSamples; i++) {
			const p2f sample = sampler.next2D();
			const p3f pShading = sInteraction.p;
			float emitterPDF;
			size_t id = selectEmitter(sampler.next(), emitterPDF);
			const Emitter& em = getEmitterByID(id);
			const v3f emitterCenter = scene.getShapeCenter(em.shapeID);
			float emitterRadius = scene.getShapeRadius(em.shapeID);
			v3f pos, ne, wiW;
			float pdf;

			//start of sample sphere by area
			v3f worldCoord = emitterRadius * Warp::squareToUniformSphere(sample) + emitterCenter;
			pos = worldCoord;
			ne = glm::normalize(emitterRadius * Warp::squareToUniformSphere(sample));
			wiW = glm::normalize(worldCoord - pShading); 
			pdf = 1 / (4 * M_PI*emitterRadius*emitterRadius);
			//end of sample sphere by area

			sInteraction.wi = sInteraction.frameNs.toLocal(wiW);
						
			Ray shadow = Ray(sInteraction.p, wiW);
			SurfaceInteraction surfaceShadow;

			if (scene.bvh->intersect(shadow, surfaceShadow)) {

				bool intersectLightSource = getEmission(surfaceShadow) != v3f(0, 0, 0);
				if (intersectLightSource) {
					float alignement = glm::dot(wiW, ne);
					if (glm::dot(emitterCenter - pShading, ne) >= 0) {
						float distance = glm::distance2(pos, pShading);
						pdf = pdf * distance / alignement;
						Lr += getBSDF(sInteraction)->eval(sInteraction) * getEmission(surfaceShadow) / (pdf*emitterPDF);;
					}
				}
			}
		}
		Lr /= noOfSamples;
		return Lr;
    }

    v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
		int noOfSamples = scene.config.integratorSettings.di.emitterSamples;
		SurfaceInteraction sInteraction;

		bool intersection = scene.bvh->intersect(ray, sInteraction);
		if (intersection) {

			bool emIntersect = getEmission(sInteraction) != v3f(0.0f);
			if (emIntersect) {
				return getEmission(sInteraction);
			}
			else {
				for (int i = 0; i < noOfSamples; i++) {

					v3f dir = Warp::squareToCosineHemisphere(sampler.next2D());
					float pdf = Warp::squareToCosineHemispherePdf(dir);

					dir = glm::normalize(sInteraction.frameNs.toWorld(dir));

					Ray shadow = Ray(sInteraction.p, dir);
					SurfaceInteraction shadowSI;

					bool shadowIntersect = scene.bvh->intersect(shadow, shadowSI);
					if (shadowIntersect) {

						bool emIntersect = getEmission(shadowSI) != v3f(0.0f);
						if (emIntersect) {
							v3f tempWI = sInteraction.frameNs.toLocal(dir);
							sInteraction.wi = tempWI;
							v3f BRDF_V = getBSDF(sInteraction)->eval(sInteraction);
							Lr += 1 / pdf * BRDF_V * getEmission(shadowSI);
						}
						else {
							//nothing to do
						}
					}
				}
			}
		}
		Lr = Lr / noOfSamples;
		return Lr;
    }

	v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);

		// TODO: Implement this
		int emitterSamples = scene.config.integratorSettings.di.emitterSamples;
		SurfaceInteraction sInteraction;

		bool isIntersection = scene.bvh->intersect(ray, sInteraction);
		if (isIntersection) {
			
			bool emIntersect = getEmission(sInteraction) != v3f(0.0f);
			if (emIntersect) {
				return getEmission(sInteraction);
			}

			else {
				for (int i = 0; i < emitterSamples; i++) {

					float pdf;
					v3f sample = getBSDF(sInteraction)->sample(sInteraction, sampler.next2D(), &pdf);						

					v3f siwi = sInteraction.wi;
					siwi = glm::normalize(sInteraction.frameNs.toWorld(siwi));

					Ray shadowRay = Ray(sInteraction.p, siwi);
					SurfaceInteraction shadowDir;

					bool isSIntersection = scene.bvh->intersect(shadowRay, shadowDir);
					if (isSIntersection) {

						bool lightSourceIntersection = getEmission(shadowDir) != v3f(0.0f);
						if (lightSourceIntersection){
							v3f shadowEmission = getEmission(shadowDir);
							Lr = Lr + sample * shadowEmission;
						}
						else {
							//nothing to do
						}
					}

				}
			}
		}
		Lr /= emitterSamples;
		return Lr;
	}

    v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
		int noOfSamples = scene.config.integratorSettings.di.emitterSamples;
		SurfaceInteraction sInteraction;

		if (scene.bvh->intersect(ray, sInteraction)) {

			bool intersectEmitter = getEmission(sInteraction) != v3f(0.0f);
			if (getEmission(sInteraction) != v3f(0.0f)) {
				return getEmission(sInteraction);
			}
		}

		for (int i = 0; i < noOfSamples; i++) {

			const p2f s = sampler.next2D();
			const p3f pShading = sInteraction.p;
			float emPdf;
			size_t id = selectEmitter(sampler.next(), emPdf);
			const Emitter& em = getEmitterByID(id);
			const v3f emitterCenter = scene.getShapeCenter(em.shapeID);
			float emitterRadius = scene.getShapeRadius(em.shapeID);
			v3f wiW;
			float pdf;

			//start of sampling sphere by solid angle
			float distance = glm::distance(emitterCenter, pShading);
			float maxCosTheta = distance / sqrt(pow(distance, 2) + pow(emitterRadius, 2));

			float cosT = (1 - s.x) + s.x * maxCosTheta;
			float sinT = glm::sqrt(glm::fmax(0.f, 1 - cosT * cosT));
			float p = s.y * 2 * M_PI;
			float x = sinT * cos(p);
			float y = sinT * sin(p);
			float z = cosT;

			v3f wiL = v3f(x, y, z);

			v3f emitterCenterDirection = glm::normalize(emitterCenter - pShading);
			glm::quat temp = glm::quat(v3f(0, 0, 1), emitterCenterDirection);
			wiW = glm::mat4(temp)*v4f(wiL, 1);
			wiW = glm::normalize(wiW);
			pdf = pdf = 1 / (2 * M_PI*(1 - maxCosTheta));
			//end of sampling phere by solid angle

			sInteraction.wi = sInteraction.frameNs.toLocal(wiW);

			Ray shadow = Ray(sInteraction.p, wiW);
			SurfaceInteraction shadowInteraction;

			//check if intersecting
			bool shIntersect = scene.bvh->intersect(shadow, shadowInteraction);
			if (shIntersect) {

				bool hitLightSource = getEmission(shadowInteraction) != v3f(0.0f);
				if (hitLightSource){
					float temp = pdf * emPdf;
					Lr += getBSDF(sInteraction)->eval(sInteraction) * getEmission(shadowInteraction) / temp;  						
				}
			}
		}
		Lr = Lr / noOfSamples;
		return Lr;
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {
		v3f result = v3f(0.f);
        v3f Lr(0.f);
		int noOfSamplesEm = scene.config.integratorSettings.di.emitterSamples;
		int noOfSamplesBSDF = scene.config.integratorSettings.di.bsdfSamples;

		SurfaceInteraction sInteraction;

		bool intersection = scene.bvh->intersect(ray, sInteraction);
		if (intersection) {
			if (getEmission(sInteraction) != v3f(0.0f)) {
				return  getEmission(sInteraction);
			}
		}

		v3f EmLR = v3f(0.f);
		for (int i = 0; i < noOfSamplesEm; i++) {
			const p2f sample = sampler.next2D();
			const p3f pShading = sInteraction.p;
			float emPdf;
			size_t id = selectEmitter(sampler.next(), emPdf);
			const Emitter& em = getEmitterByID(id);
			const v3f emitterCenter = scene.getShapeCenter(em.shapeID);
			float emitterRadius = scene.getShapeRadius(em.shapeID);
			v3f wiW;
			float pdf;

			//start of sampling sphere by solid angle
			float distance = glm::distance(emitterCenter, pShading);
			float maxCosTheta = distance / sqrt(pow(distance, 2) + pow(emitterRadius, 2));

			float cosT = (1 - sample.x) + sample.x * maxCosTheta;
			float sinT = glm::sqrt(glm::fmax(0.f, 1 - cosT * cosT));
			float p = sample.y * 2 * M_PI;
			float x = sinT * cos(p);
			float y = sinT * sin(p);
			float z = cosT;

			v3f wiL = v3f(x, y, z);

			v3f emitterCenterDirection = glm::normalize(emitterCenter - pShading);
			glm::quat temp = glm::quat(v3f(0, 0, 1), emitterCenterDirection);
			wiW = glm::mat4(temp)*v4f(wiL, 1);
			wiW = glm::normalize(wiW);
			pdf = pdf = 1 / (2 * M_PI*(1 - maxCosTheta));
			//end of sampling sphere by solid angle

			sInteraction.wi = sInteraction.frameNs.toLocal(wiW);					
			Ray shadow = Ray(sInteraction.p, wiW);
			SurfaceInteraction shadowInteraction;

			bool intersection = scene.bvh->intersect(shadow, shadowInteraction);
			if (intersection) {
				v3f emission = getEmission(shadowInteraction);

				bool hitLightSource = emission != v3f(0.0f);
				if (hitLightSource){
					float BSDFFunction = getBSDF(sInteraction)->pdf(sInteraction);
					pdf *= emPdf;
					EmLR += getBSDF(sInteraction)->eval(sInteraction) * emission / pdf * balanceHeuristic(noOfSamplesEm, pdf, noOfSamplesBSDF, BSDFFunction);
				}

				else {
					//nothing to do 
				}
			}
		}

		v3f BSDFLR = v3f(0.f);
		for (int y = 0; y < noOfSamplesBSDF; y++) {

			float pdf;
			v3f sample = getBSDF(sInteraction)->sample(sInteraction, sampler.next2D(), &pdf);
	
			v3f wi = sInteraction.wi;
			wi = glm::normalize(sInteraction.frameNs.toWorld(wi));

			Ray shadow = Ray(sInteraction.p, wi);
			SurfaceInteraction shadowInteraction;

			bool intersect = scene.bvh->intersect(shadow, shadowInteraction);
			if (intersect) {

				bool hitLightSource = getEmission(shadowInteraction) != v3f(0.0f);
				if (hitLightSource){
					const v3f pShading = sInteraction.p;
								
					const Emitter& em = getEmitterByID(getEmitterIDByShapeID(shadowInteraction.shapeID));
					float emitterFunction = getEmitterPdf(em);
					const v3f emitterCenter = scene.getShapeCenter(em.shapeID);
					float emitterRadius = scene.getShapeRadius(em.shapeID);
					float distance = glm::distance(emitterCenter, pShading);
					float maxCosTheta = distance / sqrt(distance*distance + emitterRadius * emitterRadius);
					float solidAngleFunction = pdf = 1 / (2 * M_PI*(1 - maxCosTheta))*emitterFunction;
					float balanceH = balanceHeuristic(noOfSamplesBSDF, pdf, noOfSamplesEm, solidAngleFunction);
					BSDFLR += sample * getEmission(shadowInteraction) * balanceH;
				}

				else {
					//nothing to do
				}
			}
		}

		if (noOfSamplesEm <= 0) {
			result = BSDFLR / noOfSamplesBSDF;
		}
		else if (noOfSamplesBSDF <= 0) {
			result = EmLR / noOfSamplesEm;
		}
		else {
			result = (EmLR / noOfSamplesEm) + (BSDFLR / noOfSamplesBSDF);
		}
		return result;
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        if (m_samplingStrategy == "mis")
            return this->renderMIS(ray, sampler);
        else if (m_samplingStrategy == "area")
            return this->renderArea(ray, sampler);
        else if (m_samplingStrategy == "solidAngle")
            return this->renderSolidAngle(ray, sampler);
        else if (m_samplingStrategy == "cosineHemisphere")
            return this->renderCosineHemisphere(ray, sampler);
        else if (m_samplingStrategy == "bsdf")
            return this->renderBSDF(ray, sampler);
        std::cout << "Error: wrong strategy" << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t m_emitterSamples;     // Number of emitter samples
    size_t m_bsdfSamples;        // Number of BSDF samples
    string m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END