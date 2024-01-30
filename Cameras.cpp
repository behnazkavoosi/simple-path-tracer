#include "Cameras.h"


PinholeCamera::PinholeCamera(const glm::vec3& orig, int width, int height, float fov) : origin(orig), width(width), height(height), fov(fov) {}

glm::vec3 PinholeCamera::set_camera(size_t i, size_t j) const {
    float scale = tan(deg2rad(fov * 0.5));
    float imageAspectRatio = width / (float)height;
    float dir_x = (2 * (i + 0.5) / (float)width - 1) * imageAspectRatio * scale;
    float dir_y = (1 - 2 * (j + 0.5) / (float)height) * scale;
    float dir_z = -height / (2. * tan(fov / 2.));

    return (glm::vec3(dir_x, dir_y, dir_z));
}