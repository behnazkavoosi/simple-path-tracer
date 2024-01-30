#include "Objects.h"
#include "BBox.h"


#ifndef KDNODE_H_
#define KDNODE_H_

struct IsectInfo
{
	glm::vec3 position;
	glm::vec3 normal;
	Mesh* mesh;
	float tNear = kInfinity;
	bool isEmpty;
	glm::vec3 albedo;
	glm::vec3 Ka, Kd, Ks, Tr, Ke;
	float n, ior, illum;
	MaterialType type;
	//int type;
};

class KDNode
{
public:
	BBox bbox;
	KDNode* left;
	KDNode* right;
	Mesh* meshes;
	KDNode() = default;
	virtual KDNode* build(Mesh*& meshs, int depth) const;
	virtual bool kIntersect(const KDNode* node, const glm::vec3& orig, const glm::vec3& dir, IsectInfo* isct)const;

};


#endif