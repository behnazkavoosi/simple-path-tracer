#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <chrono>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <opencv2/opencv.hpp>
#include "Objects.h"
#include "Light.h"
#include "Ray.h"
#include "Cameras.h"
#include "BBox.h"
#include "KdNode.h"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include "tiny_obj_loader.h"

std::vector<Triangle> load_mesh(std::string const& inputfile, glm::mat4& o2w) {
    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = "./";  // Path to material files

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(inputfile, reader_config)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning();
    }

    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();
    auto& materials = reader.GetMaterials();

    // Loop over shapes
    int triCounter;

    std::vector<Triangle> tris;
    for (size_t s = 0; s < shapes.size(); s++) {
        if (s > 0) {
            triCounter += shapes[s].mesh.num_face_vertices.size();
        }
        else if (s == 0) {
            triCounter = shapes[s].mesh.num_face_vertices.size();
        }
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {

            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
            if (fv != 3) {
                std::runtime_error("Only triangles allowed!");
            }

            Triangle tri;
            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
                tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
                tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];

                tri.vertex[v] = glm::vec3(o2w * vec4(vx, vy, vz, 1));

                if (idx.normal_index >= 0) {
                    //std::cout<<" + normal"<<std::endl;
                    tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
                    tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
                    tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];

                    tri.normal[v] = glm::vec3(transpose(inverse(o2w)) * vec4(nx, ny, nz, 0));

                }
                else {
                    std::runtime_error("Normals expected!");
                }
            }
            index_offset += fv;
            
            // per-face material
            int materialID = shapes[s].mesh.material_ids[f];

            tri.ior = materials[materialID].ior;
            tri.Kd = glm::vec3(materials[materialID].diffuse[0], materials[materialID].diffuse[1], materials[materialID].diffuse[2]);
            tri.Ks = glm::vec3(materials[materialID].specular[0], materials[materialID].specular[1], materials[materialID].specular[2]);
            tri.Ka = glm::vec3(materials[materialID].ambient[0], materials[materialID].ambient[1], materials[materialID].ambient[2]);
            tri.Tr = glm::vec3(materials[materialID].transmittance[0], materials[materialID].transmittance[1], materials[materialID].transmittance[2]);
            tri.n = materials[materialID].shininess;
            tri.illum = materials[materialID].illum;
            tri.Ke = glm::vec3(materials[materialID].emission[0], materials[materialID].emission[1], materials[materialID].emission[2]);
            /*if (materials[materialID].name == "sphere_mirror")
                tri.type = Mirror;

            else if (materials[materialID].name == "sphere_transparent")
                tri.type = Glass;

            else if (materials[materialID].name == "phong")
                tri.type = Phong;

            else
                tri.type = Lambertian;*/
            if (materials[materialID].shininess > 1) {
                tri.type = Mirror;
                std::cout << "mirror: " << materials[materialID].name << std::endl;
            }
            else if ((materials[materialID].shininess == 1) and (materials[materialID].ior > 1)){
                tri.type = Glass;
                std::cout << "Glass: " << materials[materialID].name << std::endl;
            }
            else if ((materials[materialID].shininess == 1) and (materials[materialID].ior == 1)) {
                tri.type = Lambertian;
                std::cout << "Lambertian: " << materials[materialID].name << std::endl;
            }


            tris.push_back(tri);
        }
    }
    std::cout << "Number of Triangles: " << triCounter << std::endl;
    return tris;
}


void render(
    const Options& options,
    const KDNode*& Kdtree,
    const std::vector<std::unique_ptr<Light>>& lights)
{
    cv::Mat image(options.height, options.width, CV_8UC3);

    glm::vec3 orig;
    orig = glm::vec3(options.cameraToWorld * vec4(0, 0, 0, 1));
    int nsample = 200;
    auto timeStart = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic,1)
    for (uint32_t j = 0; j < options.height; ++j) {
        for (uint32_t i = 0; i < options.width; ++i) {

            glm::vec3 dir = PinholeCamera(glm::vec3(0, 0, 0), options.width, options.height, options.fov).set_camera(i, j);
            glm::vec3 dirT = normalize(glm::vec3(options.cameraToWorld * vec4(dir.x, dir.y, -1, 0)));
            glm::vec3 color = glm::vec3(0, 0, 0);

            #pragma omp parallel for schedule(dynamic,1)
            for (int s = 0; s < nsample; s++) {
                color = color + castRay(orig, dirT, Kdtree, lights, options, 0)/ (double)nsample;
            }

            image.at<cv::Vec3b>(Point(i, j)) = cv::Vec3b(color[2], color[1], color[0]);
        }

        fprintf(stderr, "\r%3d%c", uint32_t(j / (float)options.height * 100), '%');
    }

    auto timeEnd = std::chrono::high_resolution_clock::now();
    auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
    fprintf(stderr, "\rRendering Time: %.2f (sec)\n", passedTime / 1000);

    cv::String windowName = "Rendered Image"; //Name of the window

    cv::namedWindow(windowName); // Create a window

    cv::imshow(windowName, image); // Show our image inside the created window.


    int k = cv::waitKey(0);  // Wait for a keystroke in the window
    if (k == 's') {
        cv::Mat m = cv::Mat::zeros(options.height, options.width, CV_8UC3);
        for (int y = 0; y < options.height; ++y) {
            for (int x = 0; x < options.width; ++x) {
                cv::Vec3f src_pt = image.at<cv::Vec3b>(cv::Point(x, y));
                cv::Vec3b& dst_pt = m.at<cv::Vec3b>(cv::Point(x, y));
                dst_pt[0] = src_pt[0] ;
                dst_pt[1] = src_pt[1] ;
                dst_pt[2] = src_pt[2] ;
            }
        }
        cv::imwrite("cornell_box.png", m);
    }
    cv::destroyWindow(windowName); //destroy the created window
}

int main(int argc, char** argv)
{
    // loading gemetry
    std::vector<std::unique_ptr<KDNode>> kdNodes;
    // lights
    std::vector<std::unique_ptr<Light>> lights;


    Options options;

    options.fov = 52;
    options.width = 700;   //700
    options.height = 600;  //600
    //options.cameraToWorld[3][2] = 15;
    //options.cameraToWorld[3][1] = 5;
    options.cameraToWorld[3][2] = 3;
    options.cameraToWorld[3][1] = 1;
    options.maxDepth = 8;

    std::vector<Mesh*> meshes;
    Mesh* onemesh = new Mesh();

    glm::mat4 cube = glm::mat4(1);
    //cube = glm::translate(cube, glm::vec3(0, -100, 0));
    //cube = glm::rotate(cube, static_cast<float>(M_PI / 4), glm::vec3(0.5, 0.5, 0));
    //cube = glm::scale(cube, glm::vec3(80, 10, 40));
    //KDNode* KDTree = new KDNode();
    Mesh* mesh1 = new Mesh(load_mesh("cornellbox-mixed.obj", cube));

    if (mesh1 != nullptr)
        meshes.push_back(mesh1);

    
     

    Triangle triangle;
    for (int i = 0; i < meshes.size(); i++) {
        for (int j = 0; j < meshes[i]->triangles.size(); j++) {
            

            triangle.vertex[0] = meshes[i]->triangles[j].vertex[0];
            triangle.vertex[1] = meshes[i]->triangles[j].vertex[1];
            triangle.vertex[2] = meshes[i]->triangles[j].vertex[2];

            triangle.normal[0] = meshes[i]->triangles[j].normal[0];
            triangle.normal[1] = meshes[i]->triangles[j].normal[1];
            triangle.normal[2] = meshes[i]->triangles[j].normal[2];

            triangle.albedo = meshes[i]->albedo;
            triangle.ior = meshes[i]->triangles[j].ior;
            triangle.Kd = meshes[i]->triangles[j].Kd;
            triangle.Ks = meshes[i]->triangles[j].Ks;
            triangle.Ka = meshes[i]->triangles[j].Ka * 100;
            triangle.Tr = meshes[i]->triangles[j].Tr;
            triangle.n = meshes[i]->triangles[j].n;
            triangle.illum = meshes[i]->triangles[j].illum;
            triangle.type = meshes[i]->triangles[j].type;
            triangle.Ke = meshes[i]->triangles[j].Ke;

            onemesh->triangles.push_back(triangle);
        }

    }
    auto buildingTimeStart = std::chrono::high_resolution_clock::now();
    const KDNode* KDTree = new KDNode();
    KDTree = KDTree->build(onemesh, 0);

    auto buildingTimeEnd = std::chrono::high_resolution_clock::now();
    auto buildingPassedTime = std::chrono::duration<double, std::milli>(buildingTimeEnd - buildingTimeStart).count();
    fprintf(stderr, "\rKd-Tree Building Time: %.2f (sec)\n", buildingPassedTime / 1000);

    lights.push_back(std::unique_ptr<Light>(new PointLight(glm::vec3(120, 120, 120), cv::Vec3b(255, 255, 255), .5)));
    lights.push_back(std::unique_ptr<Light>(new PointLight(glm::vec3(-120, 120, 120), cv::Vec3b(255, 255, 255), .5)));

    // finally, render
    render(options, KDTree, lights);
    return 0;
}