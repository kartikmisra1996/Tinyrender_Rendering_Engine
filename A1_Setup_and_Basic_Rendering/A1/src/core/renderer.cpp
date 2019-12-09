/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/core.h>
#include <core/accel.h>
#include <core/renderer.h>
#include <GL/glew.h>

#ifdef __APPLE__
#include "SDL.h"
#include <OpenGL/gl.h>
#else
#ifdef _WIN32
#include <GL/gl.h>
#include "SDL.h"
#else
#include <GL/gl.h>
#include "SDL2/SDL.h"
#endif
#endif

#include <bsdfs/diffuse.h>

#include <integrators/normal.h>
#include <renderpasses/normal.h>

int const rayPerPixel = 4;

TR_NAMESPACE_BEGIN

Renderer::Renderer(const Config& config) : scene(config) { }

bool Renderer::init(const bool isRealTime, bool nogui) {
    realTime = isRealTime;
    this->nogui = nogui;
    realTimeCameraFree = false;

    if (!scene.load(isRealTime)) return false;

    if (realTime) {
        if (scene.config.renderpass == ENormalRenderPass) {
            renderpass = std::unique_ptr<NormalPass>(new NormalPass(scene));
        } else {
            throw std::runtime_error("Invalid renderpass type");
        }

        bool succ = renderpass.get()->initOpenGL(scene.config.width, scene.config.height);
        if (!succ) return false;

        return renderpass->init(scene.config);
    } else {
        if (scene.config.integrator == ENormalIntegrator) {
            integrator = std::unique_ptr<NormalIntegrator>(new NormalIntegrator(scene));
        } else {
            throw std::runtime_error("Invalid integrator type");
        }

        return integrator->init();
    }
}

v3f computeAvg(int al[rayPerPixel][3]) {
	float x;
	float y;
	float z;
	v3f result;

	for (int i = 0; i < rayPerPixel; ++i) {
		x += (float) al[i][0];
	}
	x = x / rayPerPixel;

	for (int j = 0; j < rayPerPixel; ++j) {
		y += (float)al[j][1];
	}
	y = y / rayPerPixel;

	for (int k = 0; k < rayPerPixel; ++k) {
		z += (float)al[k][2];
	}
	z = z / rayPerPixel;

	result = v3f(x, y, z);
	return result;
}

void Renderer::render() {
    if (realTime) {
        /**
         * 1) Detect and handle the quit event.
         * 2) Call the render function using renderpass->render().
         * 3) Output the rendered image into the GUI window using SDL_GL_SwapWindow(renderpass->window).
         */
		bool run = 1;
		while (run) {
			SDL_Event event;
			SDL_PollEvent(&event);
			while (SDL_PollEvent (&event)) {
				if (event.type == SDL_QUIT) {
					run = 0;
					break;
				}
			}
			renderpass->render();
			SDL_GL_SwapWindow(renderpass->window);
		}

	}
	else {

		// 1) Calculate the camera perspective, the camera-to-world transformation matrix and the aspect ratio.

		// 2) Clear integral RGB buffer.

		// 3) Loop over all pixels on the image plane.

		// 4) Generate a ray through each pixel center.

		// 5) Splat their contribution onto the image plane.

		/* the following renderer has been implemented using spp rays 
		per pixel (bonus questions), in this case, spp = 16 */

		Sampler* sampler = new Sampler(1);

		const v3f eye = scene.config.camera.o;
		const v3f center = scene.config.camera.at;
		const v3f up = scene.config.camera.up;

		glm::mat4 view = glm::lookAt(eye, center, up);

		float w = scene.config.width;
		float h = scene.config.height;

		integrator->rgb->clear();

		int numRays = scene.config.spp;

		for (int i = 0; i < scene.config.width; i++) {
			for (int j = 0; j < scene.config.height; j++) {

				float intermediateX = 0.f;
				float intermediateY = 0.f;
				float intermediateZ = 0.f;

				for (int k = 0; k < numRays; ++k) {

					float pixelX = (2 * (i + sampler->next()) / (scene.config.width) - 1) * (tan(deg2rad*(scene.config.camera.fov) / 2)) * (w/h);
					float pixelY = (1 - 2 * (j + sampler->next()) / (scene.config.height)) * (tan(deg2rad*(scene.config.camera.fov) / 2));

					v3f direction = normalize(v3f(glm::transpose(view) * v4f(pixelX, pixelY, -1, 0)));

					Ray ray = Ray(scene.config.camera.o, direction);

					v3f rgb = integrator->render(ray, *sampler);

					intermediateX += rgb.x;
					intermediateY += rgb.y;
					intermediateZ += rgb.z;
				}

				v3f result = v3f(intermediateX/numRays, intermediateY/numRays, intermediateZ/numRays);

				integrator->rgb->data[j*scene.config.width + i] = result;
			}
		}

		/*Sampler *sampler = new Sampler(260663577);
		float fov = scene.config.camera.fov;
		glm::mat4 inverseView = glm::lookAt(scene.config.camera.o, scene.config.camera.at, scene.config.camera.up);
		integrator->rgb->clear();
		float w = scene.config.width;
		float h = scene.config.height;
		v3f rgb;
		int arrayL[rayPerPixel][3];

		for (int a = 0; a < w; ++a) {
			for (int b = 0; b < h; ++b) {

				//pixel level
				for (int c = 0; c < rayPerPixel; ++c) {
					float f1 = sampler->next();
					float f2 = sampler->next();
					//std::cout << f1 << std::endl;
					float pixelX = (2 * ((a + f1) / scene.config.width - 1) * tan(fov / 2 * M_PI / 180) * (w / h));
					float pixelY = (1 - 2 * ((b + f2) / scene.config.height) * tan(fov / 2 * M_PI / 180));
					v3f directionVector = normalize(v3f(glm::transpose(inverseView)*v4f(pixelX, pixelY, -1, 0)));
					Ray ray = Ray(scene.config.camera.o, directionVector);
					//std::cout << rgb.x << std::endl;
					rgb = integrator->render(ray, *sampler);
					arrayL[c][0] = rgb.x;
					arrayL[c][1] = rgb.y;
					arrayL[c][2] = rgb.z;
				}
				rgb = computeAvg(arrayL);
				integrator->rgb->data[b*scene.config.height + a] = rgb;
			}
		}*/
    }
}

/**
 * Post-rendering step.
 */
void Renderer::cleanUp() {
    if (realTime) {
        renderpass->cleanUp();
    } else {
        integrator->cleanUp();
    }
}

BSDF::BSDF(const WorldData& d, const Config& c, const size_t matID) : worldData(d), config(c) {
    emission = glm::make_vec3(worldData.materials[matID].emission);
}

Scene::Scene(const Config& config) : config(config) { }

bool Scene::load(bool isRealTime) {
    fs::path file(config.objFile);
    bool ret = false;
    std::string err;

    if (!file.is_absolute())
        file = (config.tomlFile.parent_path() / file).make_preferred();

    tinyobj::attrib_t* attrib_ = &worldData.attrib;
    std::vector<tinyobj::shape_t>* shapes_ = &worldData.shapes;
    std::vector<tinyobj::material_t>* materials_ = &worldData.materials;
    std::string* err_ = &err;
    const string filename_ = file.string();
    const string mtl_basedir_ = file.make_preferred().parent_path().string();
    ret = tinyobj::LoadObj(attrib_, shapes_, materials_, err_, filename_.c_str(), mtl_basedir_.c_str(), true);

    if (!err.empty()) { std::cout << "Error: " << err.c_str() << std::endl; }
    if (!ret) {
        std::cout << "Failed to load scene " << config.objFile << " " << std::endl;
        return false;
    }

    // Build list of BSDFs
    bsdfs = std::vector<std::unique_ptr<BSDF>>(worldData.materials.size());
    for (size_t i = 0; i < worldData.materials.size(); i++) {
        if (worldData.materials[i].illum == 7)
            bsdfs[i] = std::unique_ptr<BSDF>(new DiffuseBSDF(worldData, config, i));
    }

    // Build list of emitters (and print what has been loaded)
    std::string nbShapes = worldData.shapes.size() > 1 ? " shapes" : " shape";
    std::cout << "Found " << worldData.shapes.size() << nbShapes << std::endl;
    worldData.shapesCenter.resize(worldData.shapes.size());

    for (size_t i = 0; i < worldData.shapes.size(); i++) {
        const tinyobj::shape_t& shape = worldData.shapes[i];
        const BSDF* bsdf = bsdfs[shape.mesh.material_ids[0]].get();
        std::cout << "Mesh " << i << ": " << shape.name << " ["
                  << shape.mesh.indices.size() / 3 << " primitives | ";

        if (bsdf->isEmissive()) {
            Distribution1D faceAreaDistribution;
            float shapeArea = getShapeArea(i, faceAreaDistribution);
            emitters.emplace_back(Emitter{i, shapeArea, bsdf->emission, faceAreaDistribution});
            std::cout << "Emitter]" << std::endl;
        } else {
            std::cout << bsdf->toString() << "]" << std::endl;
        }

        // Build world AABB and shape centers
        worldData.shapesCenter[i] = v3f(0.0);
        for (auto idx: shape.mesh.indices) {
            v3f p = {worldData.attrib.vertices[3 * idx.vertex_index + 0],
                     worldData.attrib.vertices[3 * idx.vertex_index + 1],
                     worldData.attrib.vertices[3 * idx.vertex_index + 2]};
            worldData.shapesCenter[i] += p;
            aabb.expandBy(p);
        }
        worldData.shapesCenter[i] /= float(shape.mesh.indices.size());
    }

    // Build BVH
    bvh = std::unique_ptr<TinyRender::AcceleratorBVH>(new TinyRender::AcceleratorBVH(this->worldData));

    const clock_t beginBVH = clock();
    bvh->build();
    std::cout << "BVH built in " << float(clock() - beginBVH) / CLOCKS_PER_SEC << "s" << std::endl;

    return true;
}

float Scene::getShapeArea(const size_t shapeID, Distribution1D& faceAreaDistribution) {
    const tinyobj::shape_t& s = worldData.shapes[shapeID];

    for (size_t i = 0; i < s.mesh.indices.size(); i += 3) {
        const int i0 = s.mesh.indices[i + 0].vertex_index;
        const int i1 = s.mesh.indices[i + 1].vertex_index;
        const int i2 = s.mesh.indices[i + 2].vertex_index;
        const v3f v0{worldData.attrib.vertices[3 * i0 + 0], worldData.attrib.vertices[3 * i0 + 1],
                     worldData.attrib.vertices[3 * i0 + 2]};
        const v3f v1{worldData.attrib.vertices[3 * i1 + 0], worldData.attrib.vertices[3 * i1 + 1],
                     worldData.attrib.vertices[3 * i1 + 2]};
        const v3f v2{worldData.attrib.vertices[3 * i2 + 0], worldData.attrib.vertices[3 * i2 + 1],
                     worldData.attrib.vertices[3 * i2 + 2]};

        const v3f e1{v1 - v0};
        const v3f e2{v2 - v0};
        const v3f e3{glm::cross(e1, e2)};
        faceAreaDistribution.add(0.5f * std::sqrt(e3.x * e3.x + e3.y * e3.y + e3.z * e3.z));
    }
    const float area = faceAreaDistribution.cdf.back();
    faceAreaDistribution.normalize();
    return area;
}

TR_NAMESPACE_END
