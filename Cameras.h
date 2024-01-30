#ifndef __CAMERAS_H__
#define __CAMERAS_H__

#include <glm/glm.hpp>
#include <opencv2/opencv.hpp>
#include "Objects.h"

inline
float deg2rad(const float& deg)
{
    return deg * M_PI / 180;
}


class Camera {
public:

    virtual ~Camera() = default;

    virtual glm::vec3 set_camera(size_t i, size_t j) const = 0;
};

class PinholeCamera : public Camera {
public:

    PinholeCamera(const glm::vec3& orig, int width, int height, float fov);
    virtual ~PinholeCamera() = default;
    glm::vec3 set_camera(size_t i, size_t j) const;

    glm::vec3 origin;
    int width;
    int height;
    float fov;
};

#endif