#include "Objects.h"

// Compute the roots of a quadratic equation
bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0) return false;
    else if (discr == 0) {
        x0 = x1 = -0.5 * b / a;
    }
    else {
        float q = (b > 0) ?
            -0.5 * (b + sqrt(discr)) :
            -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }

    return true;
}

Sphere::Sphere(const glm::mat4& o2w, const float& r) : Object(), radius(r), radius2(r* r)
{
    center = glm::vec3(o2w * glm::vec4(0, 0, 0, 1));
    std::cout << center[0] << ' ' << center[1] << ' ' << center[2] << std::endl;
}

bool Sphere::intersect(
    const glm::vec3& orig,
    const glm::vec3& dir,
    float& tNear) const    // not used for sphere
{
    float t0, t1; // solutions for t if the ray intersects
    // analytic solution
    glm::vec3 L = orig - center;
    float a = dot(dir, dir);
    float b = 2 * dot(dir, L);
    float c = dot(L, L) - radius2;
    if (!solveQuadratic(a, b, c, t0, t1)) return false;

    if (t0 > t1) std::swap(t0, t1);

    if (t0 < 0) {
        t0 = t1; // if t0 is negative, let's use t1 instead
        if (t0 < 0) return false; // both t0 and t1 are negative
    }

    tNear = t0;

    return true;
}

// Set surface data such as normal and texture coordinates at a given point on the surface
void Sphere::getSurfaceProperties(
    const glm::vec3& hitPoint,
    const glm::vec3& viewDirection,
    glm::vec3& hitNormal) const
{
    hitNormal = hitPoint - center;
    hitNormal = normalize(hitNormal);

}

bool rayTriangleIntersect(
    const glm::vec3& orig, const glm::vec3& dir,
    const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2,
    float& t, float& u, float& v)
{
    auto edge1 = v1 - v0;
    auto edge2 = v2 - v0;
    auto h = cross(dir, edge2);
    auto a = dot(edge1, h);
    if (a > -kEpsilon && a < kEpsilon) {
        return false;  // This ray is parallel to this triangle.
    }
    float f = 1.0f / a;
    auto s = orig - v0;
    u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) {
        return false;
    }
    auto q = cross(s, edge1);
    v = f * dot(dir, q);
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }
    // At this stage we can compute t to find out where the intersection point is on the line.
    t = f * dot(edge2, q);
    if (t > kEpsilon) {
        return true;
    }
    else  // This means that there is a line intersection but not a ray intersection.
        return false;
}


glm::vec3 Triangle::get_midpoint() const {
    return glm::vec3((bbox.max.x + bbox.min.x) / 2, (bbox.max.y + bbox.min.y) / 2, (bbox.max.z + bbox.min.z) / 2);
}

BBox Triangle::get_bounding_box() {
    bbox.computeBBoxFromOriginalPointSet(vertex);
    return bbox;
}

Mesh::Mesh(std::vector<Triangle> const& triangles) : triangles(triangles) {
}

    

   