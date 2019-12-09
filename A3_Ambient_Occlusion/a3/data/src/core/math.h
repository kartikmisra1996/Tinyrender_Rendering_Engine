/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Computes barycentric coordinates.
 */
template<class T>
inline T barycentric(const T& a, const T& b, const T& c, const float u, const float v) {
    return a * (1 - u - v) + b * u + c * v;
}

/**
 * Restricts a value to a given interval.
 */
template<class T>
inline T clamp(T v, T min, T max) {
    return std::min(std::max(v, min), max);
}

/**
 * Checks if vector is zero.
 */
inline bool isZero(const v3f v) {
    return glm::dot(v, v) < Epsilon;
}

/**
 * Generates coordinate system.
 */
inline void coordinateSystem(const v3f& a, v3f& b, v3f& c) {
    if (std::abs(a.x) > std::abs(a.y)) {
        float invLen = 1.f / std::sqrt(a.x * a.x + a.z * a.z);
        c = v3f(a.z * invLen, 0.f, -a.x * invLen);
    } else {
        float invLen = 1.f / std::sqrt(a.y * a.y + a.z * a.z);
        c = v3f(0.f, a.z * invLen, -a.y * invLen);
    }
    b = glm::cross(c, a);
}

/**
 * Converts RGB value to luminance.
 */
inline float getLuminance(const v3f& rgb) {
    return glm::dot(rgb, v3f(0.212671f, 0.715160f, 0.072169f));
}

/**
 * Pseudo-random sampler (Mersenne Twister 19937) structure.
 */
struct Sampler {
    std::mt19937 g;
    std::uniform_real_distribution<float> d;
    explicit Sampler(int seed) {
        g = std::mt19937(seed);
        d = std::uniform_real_distribution<float>(0.f, 1.f);
    }
    float next() { return d(g); }
    p2f next2D() { return {d(g), d(g)}; }
    void setSeed(int seed) {
        g.seed(seed);
        d.reset();
    }
};

/**
 * 1D discrete distribution.
 */
struct Distribution1D {
    std::vector<float> cdf{0};
    bool isNormalized = false;

    inline void add(float pdfVal) {
        cdf.push_back(cdf.back() + pdfVal);
    }

    size_t size() {
        return cdf.size() - 1;
    }

    float normalize() {
        float sum = cdf.back();
        for (float& v : cdf) {
            v /= sum;
        }
        isNormalized = true;
        return sum;
    }

    inline float pdf(size_t i) const {
        assert(isNormalized);
        return cdf[i + 1] - cdf[i];
    }

    int sample(float sample) const {
        assert(isNormalized);
        const auto it = std::upper_bound(cdf.begin(), cdf.end(), sample);
        return clamp(int(distance(cdf.begin(), it)) - 1, 0, int(cdf.size()) - 2);
    }
};


/**
 * Warping functions.
 */
namespace Warp {


	inline v3f squareToUniformSphere(const p2f& sample) {
		v3f result;
		float R = (float)sqrt((1.f - (1.f - 2.f * sample.x) * (1.f - 2.f * sample.x)));
		float p = 2 * M_PI * (sample.y);
		float x = R * cos(p);
		float y = R * sin(p);
		result = v3f(x, y, (1.f - 2.f * sample.x));
		return result;
	}

	inline float squareToUniformSpherePdf() {
		float result = 1.f / (4.f * M_PI);
		return result;
	}

	inline v3f squareToUniformHemisphere(const p2f& sample) {
		v3f result;
		float a = sample.x;
		float b = sample.y;
		float c = a;
		float R = (float)sqrt(std::max(0.f, 1.f - c * c));
		float p = 2 * M_PI * b;
		float x = R * cosf(p);
		float y = R * sinf(p);
		result = v3f(x, y, c);
		return result;
	}

	inline float squareToUniformHemispherePdf(const v3f& v) {
		float result = INV_TWOPI;
		return result;
	}

	inline v2f squareToUniformDiskConcentric(const p2f& sample) {
		v2f result;
		float a = 2.0f*sample.x - 1.0f;
		float b = 2.0f*sample.y - 1.0f;
		float p;
		float R;
		if (a == 0 && b == 0) {
			R = p = 0;
		}
		else if (a*a > b*b) {
			R = a;
			p = (M_PI/4.0f)*(b/a);
		}
		else {
			R = b;
			p = (M_PI/2.0f)-(a/b)*(M_PI/4.0f);
		}
		float x = R*cosf(p);
		float y = R*sinf(p);
		result = v2f(x, y);
		return result;
	}

	inline v3f squareToCosineHemisphere(const p2f& sample) {
		v3f result;
		v2f p = squareToUniformDiskConcentric(sample);
		float z = sqrt(std::max(0.f, 1.0f - p.x*p.x - p.y*p.y));
		result = v3f(p.x, p.y, z);
		return result;
	}

	inline float squareToCosineHemispherePdf(const v3f& v) {
		float result = 0.f;
		result = v.z*INV_PI;
		return result;
	}

	inline v3f squareToPhongLobe(const p2f& sample, int n) {
		v3f result(0.f);
		float theta = acos(pow((1 - sample.x), 1.0f / (n + 2)));
		float p = 2 * M_PI*(sample.y);
		result = v3f(sin(theta) * cos(p), sin(theta) * sin(p), cos(theta));
		return result;
	}

	inline float squareToPhongLobePdf(const v3f& v, int n) {
		float result = 0.f;
		result = (n + 2) * INV_TWOPI * pow(v.z, n);
		return result;
	}

}

TR_NAMESPACE_END