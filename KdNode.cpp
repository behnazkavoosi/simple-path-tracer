#include "KdNode.h"

KDNode* KDNode::build(Mesh*& meshs, int depth) const {
	//std::vector<Triangle> triangles1;
	KDNode* node = new KDNode();
	node->meshes = meshs;
	node->left = NULL;
	node->right = NULL;
	node->bbox = BBox();


	if (meshs->triangles.size() == 0) {

		return node;
	}
	if (meshs->triangles.size() == 1) {
		node->bbox = meshs->triangles[0].get_bounding_box();
		node->left = new KDNode();
		node->right = new KDNode();
		node->left->meshes = new Mesh();
		node->right->meshes = new Mesh();
		return node;
	}

	node->bbox = meshs->triangles[0].get_bounding_box();

	for (int i = 1; i < meshs->triangles.size(); i++) {
		//node->bbox.extendBy(meshs->triangles[i].vertex[j]);
		node->bbox.extendBy(meshs->triangles[i].get_bounding_box());
	}

	glm::vec3 midPoint = glm::vec3(0, 0, 0);
	for (int i = 0; i < meshs->triangles.size(); i++) {
		//	for (int j = 0; j < 3; j++) {
		midPoint = midPoint + meshs->triangles[i].get_midpoint();

		//	}
	}

	midPoint /= meshs->triangles.size();
	Mesh* leftMeshes = new Mesh();
	Mesh* rightMeshes = new Mesh();

	int axis = node->bbox.longest_axis();
	for (int i = 0; i < meshs->triangles.size(); i++) {

		switch (axis) {
		case 0:
			midPoint.x >= meshs->triangles[i].get_midpoint().x ? rightMeshes->triangles.push_back(meshs->triangles[i]) : leftMeshes->triangles.push_back(meshs->triangles[i]);
			break;
		case 1:
			midPoint.y >= meshs->triangles[i].get_midpoint().y ? rightMeshes->triangles.push_back(meshs->triangles[i]) : leftMeshes->triangles.push_back(meshs->triangles[i]);
			break;
		case 2:
			midPoint.z >= meshs->triangles[i].get_midpoint().z ? rightMeshes->triangles.push_back(meshs->triangles[i]) : leftMeshes->triangles.push_back(meshs->triangles[i]);
			break;
		}
	}

	if (leftMeshes->triangles.size() == 0 && rightMeshes->triangles.size() > 0) leftMeshes = rightMeshes;
	if (rightMeshes->triangles.size() == 0 && leftMeshes->triangles.size() > 0) rightMeshes = leftMeshes;
	int matches = 0;

	for (int i = 0; i < leftMeshes->triangles.size(); i++) {
		for (int j = 0; j < rightMeshes->triangles.size(); j++) {
			if ((leftMeshes->triangles[i].vertex[0] == rightMeshes->triangles[j].vertex[0]) &&
				(leftMeshes->triangles[i].vertex[1] == rightMeshes->triangles[j].vertex[1]) &&
				(leftMeshes->triangles[i].vertex[2] == rightMeshes->triangles[j].vertex[2]))
			{
				matches++;

			}

		}
	}



	if ((double)matches / leftMeshes->triangles.size() < 0.5 && (double)matches / rightMeshes->triangles.size() < 0.5) {

		node->left = build(leftMeshes, depth + 1);
		node->right = build(rightMeshes, depth + 1);
	}
	else {
		node->left = new KDNode();
		node->right = new KDNode();
		node->left->meshes = new Mesh();
		node->right->meshes = new Mesh();
	}

	return node;
}

bool KDNode::kIntersect(const KDNode* node, const glm::vec3& orig, const glm::vec3& dir, IsectInfo* isct) const{
	bool hit_tri = false;
	float t = kInfinity, u = 0, v = 0;
	uint32_t idx = 0;
	glm::vec3 vchit;
	isct->isEmpty = true;
	IsectInfo isect;

	if (node->bbox.bIntersect(orig, dir, &vchit)) {



		if (node->left->meshes->triangles.size() > 0 || node->right->meshes->triangles.size() > 0) {

			bool hitleft = kIntersect(node->left, orig, dir, isct);
			bool hitright = kIntersect(node->right, orig, dir, isct);
			return hitleft || hitright;
		}

		else {

			for (int i = 0; i < node->meshes->triangles.size(); i++) {
				if (rayTriangleIntersect(orig, dir, node->meshes->triangles[i].vertex[0], node->meshes->triangles[i].vertex[1], node->meshes->triangles[i].vertex[2], t, u, v)) {
					if (t < isect.tNear) {

						return true;
					}
				}
				//idx++;
			}
		}
	}


	return false;
}

