#ifndef __LIGHT_H__
#define __LIGHT_H__

#include <glm/glm.hpp>
#include <opencv2/opencv.hpp>
#include "iostream"
#include "fstream"
#include "Objects.h"


// Light base class
class Light
{
public:
    Light(const glm::vec3& pos, const cv::Vec3b& c = cv::Vec3b(1, 1, 1), const float& i = 1) : position(pos), color(c), intensity(i) {}
    virtual ~Light() {}
    virtual void illuminate(const glm::vec3& P, glm::vec3&, glm::vec3&, float&) const = 0;
    cv::Vec3b color;
    float intensity;
    glm::vec3 position;
};


class PointLight : public Light
{
    
public:
    PointLight(const glm::vec3& p, const cv::Vec3b& c = cv::Vec3b(1, 1, 1), const float& i = 1);
    
    void illuminate(const glm::vec3& P, glm::vec3& lightDir, glm::vec3& lightIntensity, float& distance) const;
    glm::vec3 pos;
};

class DistantLight : public Light
{
    
public:
    DistantLight(const glm::vec3& d, const cv::Vec3b& c = cv::Vec3b(1, 1, 1), const float& i = 1);

    void illuminate(const glm::vec3& P, glm::vec3& lightDir, glm::vec3& lightIntensity, float& distance) const;
    glm::vec3 dir;
};

#endif
