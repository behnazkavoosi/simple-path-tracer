#ifndef __BBOX_H
#define __BBOX_H

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>


static const float kInfinity = std::numeric_limits<float>::max();

class BBox
{
public:
    BBox() = default;

    BBox(glm::vec3 min_, glm::vec3 max_)
    {
        bounds[0] = min_;
        bounds[1] = max_;
    }

    virtual void extendBy(BBox exBox);
    virtual glm::vec3 getMidpoint(const glm::vec3& p);
    virtual int longest_axis() const;
    virtual void computeBBoxFromOriginalPointSet(glm::vec3 vertices[3]);
    virtual void extrameDistanceAlongDir(glm::vec3 dir, glm::vec3 vertices[3], unsigned int* min, unsigned int* max);
    virtual bool bIntersect(const glm::vec3& orig, const glm::vec3& dir, glm::vec3* vcHit) const;

    glm::vec3 min;
    glm::vec3 max;
    glm::vec3 bounds[2] = {glm::vec3(kInfinity, kInfinity, kInfinity), glm::vec3(-kInfinity, -kInfinity, -kInfinity)};
};


#endif

