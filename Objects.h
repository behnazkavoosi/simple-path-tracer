#ifndef __OBJECTS_H__
#define __OBJECTS_H__


#include "iostream"
#include "fstream"
#include "BBox.h"
#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <C:\Users\behna\Documents\Courses\Advanced Image Synthesis and Global Illumination\Topic 5\MCPT-master\MCPT\MCPT\assimp\include\assimp/Importer.hpp>      
#include <C:\Users\behna\Documents\Courses\Advanced Image Synthesis and Global Illumination\Topic 5\MCPT-master\MCPT\MCPT\assimp\include\assimp/scene.h>           
#include <C:\Users\behna\Documents\Courses\Advanced Image Synthesis and Global Illumination\Topic 5\MCPT-master\MCPT\MCPT\assimp\include\assimp/postprocess.h>     
#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/matrix_transform.hpp>

using std::string;
using std::cout;
using std::endl;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const float kEpsilon = 1e-8;


enum MaterialType {
    Lambertian, Glass, Mirror, Phong, Conductor
};

class Object
{
public:
    Object() {

    }
    virtual ~Object() = default;
    bool intersect(const glm::vec3&, const glm::vec3&, float&, uint32_t&, glm::vec2&) const;
    void getSurfaceProperties(const glm::vec3&, const glm::vec3&, const uint32_t&, const glm::vec2&, glm::vec3&, glm::vec2&) const;
    glm::mat4 objectToWorld, worldToObject;
    // BBox<> bbox;
	MaterialType type = Lambertian;
    glm::vec3 albedo = glm::vec3(0, 0, 0);
	glm::vec3 Kd = glm::vec3(0.8); // phong model diffuse weight
	glm::vec3 Ks = glm::vec3(0.2); // phong model specular weight
	glm::vec3 Ka = glm::vec3(0.2); // phong model specular weight
	glm::vec3 Tr = glm::vec3(0.2); // phong model specular weight
    glm::vec3 Ke = glm::vec3(0);
    float n = 1;   // phong specular exponent
    float ior = 1;
	float illum = 1;
};

// Sphere class. A sphere type object
class Sphere : public Object
{
public:
    Sphere(const glm::mat4& o2w, const float& r);

    // Ray-sphere intersection test
    bool intersect(const glm::vec3& orig, const glm::vec3& dir, float& tNear) const;

    // Set surface data such as normal and texture coordinates at a given point on the surface
    void getSurfaceProperties(const glm::vec3& hitPoint, const glm::vec3& viewDirection, glm::vec3& hitNormal) const;

    float radius, radius2;
    glm::vec3 center;
};




class Triangle : public Object {
public:

    glm::vec3 vertex[3];
	glm::vec3 normal[3];
    BBox bbox;
    glm::vec3 get_midpoint() const;
    BBox get_bounding_box();
	int mtl;
};



class Mesh : public Object {
public:
	
    Mesh(std::vector<Triangle> const& triangles);
	Mesh() = default;
    std::vector<Triangle> triangles;

};


bool rayTriangleIntersect(
    const glm::vec3& orig, const glm::vec3& dir,
    const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2,
    float& t, float& u, float& v);

#endif
