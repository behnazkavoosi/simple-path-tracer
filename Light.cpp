#include "Light.h"

PointLight::PointLight(const glm::vec3& p, const cv::Vec3b& c, const float& i) : Light(pos, c, i)
{
    pos = p;
}

void PointLight::illuminate(const glm::vec3& P, glm::vec3& lightDir, glm::vec3& lightIntensity, float& distance) const
{

    lightDir = normalize(P - pos);
    float r2 = lightDir[0] * lightDir[0] + lightDir[1] * lightDir[1] + lightDir[2] * lightDir[2];
    distance = sqrt(r2);
    lightDir.x /= distance, lightDir.y /= distance, lightDir.z /= distance;
    // avoid division by 0
    lightIntensity = glm::vec3(color[0] * intensity, color[1] * intensity, color[2] * intensity);
}

DistantLight::DistantLight(const glm::vec3& d, const cv::Vec3b& c, const float& i) : Light(dir, c, i)
{
    dir = normalize(d); // in case the matrix scales the light
}

void DistantLight::illuminate(const glm::vec3& P, glm::vec3& lightDir, glm::vec3& lightIntensity, float& distance) const
{
    lightDir = dir;
    lightIntensity = glm::vec3(color[0] * intensity, color[1] * intensity, color[2] * intensity);
    distance = kInfinity;
}