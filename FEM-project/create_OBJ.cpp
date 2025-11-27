#include <iostream>
#include <map>
#include <set>

#include "create_OBJ.h"

void createOBJ(hed::Triangulation& triang, std::ofstream& objfile) {
    const list<hed::Node*>* nodelist = triang.getNodes();
    const list<hed::Edge*>& leading_edges = triang.getLeadingEdges();

    /// Add vertices
    int node_ind = 1;
    std::map<int, int> node_idx;
    std::set<int> idx_set;
    for (hed::Node* node : *nodelist) {
        objfile << "v " << node->x() << " " << node->z() << " " << node->y() << "\n";

        node_idx.insert({node->id(), node_ind});
        node_ind++;

        idx_set.insert(node->id());
    }

    /// Add normals
    for (hed::Edge* edge : leading_edges) {
        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();

        Eigen::Vector3d a(N2->x() - N1->x(), N2->y() - N1->y(), N2->z() - N1->z());
        Eigen::Vector3d b(N3->x() - N1->x(), N3->y() - N1->y(), N3->z() - N1->z());
        Eigen::Vector3d normal = a.cross(b);

        objfile << "vn " << normal.x() << " " << normal.z() << " " << normal.y() << "\n";
    }

    /// Add faces
    int norm_ind = 1; // counts the index of a normal vector
    for (hed::Edge* edge : leading_edges) {
        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();

        objfile << "f "  <<
            node_idx[N1->id()] << "//" << norm_ind << " " <<
            node_idx[N2->id()] << "//" << norm_ind << " " <<
            node_idx[N3->id()] << "//" << norm_ind << "\n";

        norm_ind++;
    }

    // std::cout << "id:\n";
    // for(int id : idx_set){
    //     std::cout << id << ",\n";
    // }

}
