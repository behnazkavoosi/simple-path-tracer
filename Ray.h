#ifndef __RAY_H__
#define __RAY_H__

#include <random>
#include <glm/glm.hpp>
#include <opencv2/opencv.hpp>
#include "iostream"
#include "fstream"
#include "Objects.h"
#include "Light.h"
#include "KdNode.h"

enum RayType { kPrimaryRay, kShadowRay };

static const glm::vec3 kDefaultBackgroundColor = glm::vec3(128, 128, 128);

std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);

struct Options
{
    uint32_t width = 640;
    uint32_t height = 480;
    float fov = 90;
    glm::vec3 backgroundColor = kDefaultBackgroundColor;
    glm::mat4 cameraToWorld;
    float bias = 0.001;
    uint32_t maxDepth = 5;
};

inline
float clamp(const float& lo, const float& hi, const float& v)
{
    return std::max(lo, std::min(hi, v));
}

using namespace glm;
using namespace cv;

void createCoordinateSystem(const glm::vec3& N, glm::vec3& Nt, glm::vec3& Nb)
{
    if (std::fabs(N.x) > std::fabs(N.y))
        Nt = glm::vec3(N.z, 0, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
    else
        Nt = glm::vec3(0, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z);
    Nb = cross(N, Nt);
}

glm::vec3 uniformSampleHemisphere(const float& r1, const float& r2)
{
    float sinTheta = sqrtf(1 - r1 * r1);
    float phi = 2 * M_PI * r2;
    float x = sinTheta * cosf(phi);
    float z = sinTheta * sinf(phi);
    return glm::vec3(x, r1, z);
}

// Compute reflection direction
glm::vec3 reflect(const glm::vec3& I, const glm::vec3& N)
{
    return I - 2 * dot(I, N) * N;
}

// Compute refraction direction
glm::vec3 refract(const glm::vec3& I, const glm::vec3& N, const float& ior)
{
    float cosi = clamp(-1, 1, dot(I, N));
    float etai = 1, etat = ior;
    glm::vec3 n = N;
    if (cosi < 0) { cosi = -cosi; }
    else { std::swap(etai, etat); n = -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? glm::vec3(1, 0, 0) : eta * I + (eta * cosi - sqrtf(k)) * n;
}


// Evaluate Fresnel equation (ration of reflected light for a given incident direction and surface normal)
void fresnel(const glm::vec3& I, const glm::vec3& N, const float& ior, float& kr)
{
    float cosi = clamp(-1, 1, dot(I, N));
    float etai = 1, etat = ior;
    if (cosi > 0) { std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        kr = 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
}

inline float modulo(const float& f)
{
    return f - std::floor(f);
}

inline double Max(glm::vec3 x) {

    if ((x[0] >= x[1]) && (x[0] >= x[2])) {
        return x[0];
    }
    else if ((x[1] >= x[0]) && (x[1] >= x[2])) {
        return x[1];
    }
    else if ((x[2] >= x[0]) && (x[2] >= x[1])) {
        return x[2];
    }
}

// Schlick Fresnel approximation
inline float fresnelConductor(float f0, float vdoth)
{
    // F(v, h) = f₀ + (1 - f₀) * (1 - v·h)⁵

    return f0 + (1 - f0) * pow(1 - vdoth, 5);
}

// Smith GGX geometry factor for either shadowing or masking
inline float ggxNormalDistribution(float alpha, float ndotx)
{
    //
    // G₁(n, x) = (2 n·x) / (n·x + sqrt(α² + (1 - α²)(n·x)²))

    const float aa = alpha * alpha;

    return 2 * ndotx / (ndotx + sqrt(aa + (1 - aa) * ndotx * ndotx));
}

// Smith GGX visibility factor for either shadowing or masking
inline float V1(float alpha, float ndotx)
{
    //
    // V₁(n, x) = 1 / (n·x + sqrt(α² + (1 - α²)(n·x)²)

    const float aa = alpha * alpha;

    return 1 / (ndotx + sqrt(aa + (1 - aa) * ndotx * ndotx));
}

// Smith GGX visibility factor for both shadowing and masking
inline float ggxVisibilityTerm(float alpha, float ndotl, float ndotv)
{
    //
    // V = V₁(n, l) V₁(n, v)

    return V1(alpha, ndotl) * V1(alpha, ndotv);
}


glm::vec3 castRay(
    const glm::vec3& orig, const glm::vec3& dir,
    const KDNode*& Kdtree,
    const std::vector<std::unique_ptr<Light>>& lights,
    const Options& options,
    const uint32_t& depth)
{

    IsectInfo ist;
    if (depth >= options.maxDepth) return options.backgroundColor;
    glm::vec3 hitColor = options.backgroundColor;

    float tNear = kInfinity;
    vec2 uv;
    uint32_t triIndex;
    float t, u, v;

    int count = 0;
    if (!(Kdtree->kIntersect(Kdtree, orig, dir, &ist)))
        return options.backgroundColor;
    else {
        hitColor = glm::vec3(0, 0, 0);

        for (int i = 0; i < Kdtree->meshes->triangles.size(); i++) {

            if (!(rayTriangleIntersect(orig, dir, Kdtree->meshes->triangles[i].vertex[0], Kdtree->meshes->triangles[i].vertex[1], Kdtree->meshes->triangles[i].vertex[2], t, u, v)))
                continue;
            if (t < ist.tNear) {

                ist.tNear = t;
                ist.isEmpty = false;
                ist.position = orig + t * dir;
                ist.normal = normalize(Kdtree->meshes->triangles[i].normal[0] * (1 - u - v) +
                    Kdtree->meshes->triangles[i].normal[1] * u + Kdtree->meshes->triangles[i].normal[2] * v);
                ist.mesh = Kdtree->meshes;
                ist.ior = Kdtree->meshes->triangles[i].ior;
                ist.Kd = Kdtree->meshes->triangles[i].Kd;
                ist.Ks = Kdtree->meshes->triangles[i].Ks;
                ist.Ka = Kdtree->meshes->triangles[i].Ka;
                ist.Ke = Kdtree->meshes->triangles[i].Ke;
                ist.Tr = Kdtree->meshes->triangles[i].Tr;
                ist.n = Kdtree->meshes->triangles[i].n;
                ist.illum = Kdtree->meshes->triangles[i].illum;  
                ist.type = Kdtree->meshes->triangles[i].type;
            }

            /*if ((Max(ist.Ka) > 1) and (Max(ist.Kd) == 0)) {
                return glm::vec3(255,255,255);
            }*/
            if ((Max(ist.Ke) > 1)) {
                return glm::vec3(255, 255, 255);
            }

            if (depth > 5) {
                if (rand() / (double)RAND_MAX < Max(ist.Kd))	
                    ist.Kd = ist.Kd * (1 / Max(ist.Kd));
                else 
                    return glm::vec3(ist.Ka[0], ist.Ka[1], ist.Ka[2]);
            }

            if (depth > 20) {
                return glm::vec3(ist.Ka[0], ist.Ka[1], ist.Ka[2]);
            }

            glm::vec3 normal_real = dot(ist.normal, dir) < 0 ? ist.normal : ist.normal * -1.0;
            switch (ist.type) {

            case Lambertian:
            {
                double r1 = ((float)rand() / RAND_MAX);
                double r2 = ((float)rand() / RAND_MAX);
                float inv_uniform_pdf = 2 * M_PI;
                float inv_cos_pdf = M_PI / r1;
                double r2_sqrt = sqrt(r2);
                glm::vec3 Nt, Nb;
                glm::vec3 w = normal_real;
                glm::vec3 u = (fabs(w.x) > 0.1 ? glm::vec3(0.0, 1.0, 0.0) : normalize(cross(glm::vec3(1.0, 0.0, 0.0), w)));
                glm::vec3 v = cross(w, u);
                glm::vec3 dir1 = normalize((u * cos(2 * M_PI * r1) * r2_sqrt + v * sin(2 * M_PI * r1) * r2_sqrt + w * sqrt(1 - r2)));

                /*if (std::fabs(ist.normal.x) > std::fabs(ist.normal.y))
                    Nt = glm::vec3(ist.normal.z, 0, -ist.normal.x) / sqrtf(ist.normal.x * ist.normal.x + ist.normal.z * ist.normal.z);
                else
                    Nt = glm::vec3(0, -ist.normal.z, ist.normal.y) / sqrtf(ist.normal.y * ist.normal.y + ist.normal.z * ist.normal.z);
                Nb = cross(ist.normal, Nt);

                glm::vec3 sampleWorld(
                    sample.x * Nb.x + sample.y * ist.normal.x + sample.z * Nt.x,
                    sample.x * Nb.y + sample.y * ist.normal.y + sample.z * Nt.y,
                    sample.x * Nb.z + sample.y * ist.normal.z + sample.z * Nt.z);*/

                glm::vec3 xx = castRay(ist.position + dir1 * options.bias, dir1,
                    Kdtree, lights, options, depth + 1);

                hitColor = (glm::vec3(ist.Kd[0] * xx[0] + ist.Ka[0],
                    ist.Kd[1] * xx[1] + ist.Ka[1],
                    ist.Kd[2] * xx[2] + ist.Ka[2]));// *inv_cos_pdf* r1 / M_PI );

                break;

            }
            
            case Mirror:
            {
                glm::vec3 R = normalize(reflect(dir, ist.normal));
                glm::vec3 xx = castRay(ist.position + options.bias * ist.normal, R, Kdtree, lights, options, depth + 1);

                hitColor = glm::vec3(ist.Ka[0] + xx[0], ist.Ka[1] + xx[1], ist.Ka[2] + xx[2]);

                break;
            }

            case Glass:
            {
                glm::vec3 refractionColor = vec3(0), reflectionColor = vec3(0);

                // compute fresnel
                float kr;
                fresnel(dir, ist.normal, ist.ior, kr);
                bool outside = dot(dir, ist.normal) < 0;
                glm::vec3 bias = options.bias * ist.normal;
                // compute refraction if it is not a case of total internal reflection
                if (kr < 1) {
                    glm::vec3 refractionDirection = normalize(refract(dir, ist.normal, ist.ior));
                    glm::vec3 refractionRayOrig = outside ? ist.position - bias : ist.position + bias;
                    refractionColor = castRay(refractionRayOrig, refractionDirection, Kdtree, lights, options, depth + 1);
                }

                glm::vec3 reflectionDirection = normalize(reflect(dir, ist.normal));
                glm::vec3 reflectionRayOrig = outside ? ist.position + bias : ist.position - bias;
                reflectionColor = castRay(reflectionRayOrig, reflectionDirection, Kdtree, lights, options, depth + 1);

                hitColor = reflectionColor * kr + refractionColor * (1 - kr);

                break;   
            }

            case Phong:
            {
                double r1 = ((float)rand() / RAND_MAX);
                double r2 = ((float)rand() / RAND_MAX);
                float inv_uniform_pdf = 2 * M_PI;
                float inv_cos_pdf = M_PI / r1;
                double r2_sqrt = sqrt(r2);
                glm::vec3 w = normal_real;
                glm::vec3 u = (fabs(w.x) > 0.1 ? glm::vec3(0.0, 1.0, 0.0) : normalize(cross(glm::vec3(1.0, 0.0, 0.0), w)));
                glm::vec3 v = cross(w, u);
                glm::vec3 dir1 = normalize((u * cos(2 * M_PI * r1) * r2_sqrt + v * sin(2 * M_PI * r1) * r2_sqrt + w * sqrt(1 - r2)));
                glm::vec3 diffuse;
                glm::vec3 specular;
                IsectInfo isct;

                for (uint32_t i = 0; i < lights.size(); ++i) {

                    glm::vec3 lightDir;
                    glm::vec3 lightIntensity;

                    lights[i]->illuminate(ist.position, lightDir, lightIntensity, isct.tNear);
                    bool visibility = !(Kdtree->kIntersect(Kdtree, ist.position + options.bias * ist.normal, -lightDir, &isct));

                    // compute the diffuse component
                    diffuse += visibility * ist.Ka * lightIntensity * std::max(0.f, dot(ist.normal, -lightDir));

                    // compute the specular component
                    glm::vec3 R = reflect(lightDir, ist.normal);
                    specular += visibility * lightIntensity * std::pow(std::max(0.f, dot(R, -dir)), ist.n);  
                }

                glm::vec3 direct = glm::vec3(diffuse[0] * ist.Kd[0] + specular[0] * ist.Ks[0], diffuse[1] * ist.Kd[1] + specular[1] * ist.Ks[1], diffuse[2] * ist.Kd[2] + specular[2] * ist.Ks[2]);

                glm::vec3 indirect = castRay(ist.position + dir1 * options.bias, dir1, Kdtree, lights, options, depth + 1);

                hitColor = (direct + indirect) * inv_cos_pdf * r1 / M_PI;
                break;
            }

            case Conductor:
            {
                glm::vec3 R = normalize(reflect(dir, ist.normal));
                glm::vec3 xx = castRay(ist.position + options.bias * ist.normal, R, Kdtree, lights, options, depth + 1);

                hitColor = glm::vec3(ist.Kd[0] * xx[0], ist.Kd[1] * xx[1], ist.Kd[2] * xx[2]);

                break;
                
                /*double r1 = ((float)rand() / RAND_MAX);
                double r2 = ((float)rand() / RAND_MAX);
                double r2_sqrt = sqrt(r2);
                glm::vec3 w = normal_real;
                glm::vec3 u = (fabs(w.x) > 0.1 ? glm::vec3(0.0, 1.0, 0.0) : normalize(cross(glm::vec3(1.0, 0.0, 0.0), w)));
                glm::vec3 v = cross(w, u);
                glm::vec3 R = normalize((u * cos(2 * M_PI * r1) * r2_sqrt + v * sin(2 * M_PI * r1) * r2_sqrt + w * sqrt(1 - r2)));
                
                glm::vec3 wI = -dir; // reverse incoming ray, so it points out of surface
                float3 h = normalize(wI + R);
                float NdotO = dot(ist.normal, h);
                float NdotI = dot(ist.normal, wI);
                float HdotO = dot(ist.normal, R);
                float NdotH = dot(ist.normal, h);
                
                    float rs = ist.n * ist.n;
                    float F = fresnelConductor(1, HdotO);
                    float D = ggxNormalDistribution(rs, NdotH);
                    float G = ggxVisibilityTerm(rs, NdotO, NdotI);
                    float bsdf = F * D * G / (4.0 * NdotI);
                    float pdf = D * NdotH / (4.0 * HdotO);
                    float weight = F * ((G * HdotO) / (NdotI * NdotH));
                
                    glm::vec3 xx = castRay(ist.position + options.bias * ist.normal, h, Kdtree, lights, options, depth + 1);
                    //hitColor = glm::vec3(ist.Kd[0] * xx[0], ist.Kd[1] * xx[1], ist.Kd[2] * xx[2]);

                //glm::vec3 indirect = castRay(ist.position + wO * options.bias, wO, Kdtree, lights, options, depth + 1);
                hitColor = xx * bsdf * ist.Kd / pdf;
               // */
                break;
            }

            default:
                break;
            }
        }
    }
            
   
    return hitColor;
}


#endif
