#ifndef CREATE_OBJ_H
#define CREATE_OBJ_H

#include <ttl/halfedge/HeTriang.h>
#include <Eigen/Dense>

#include <fstream>



void createOBJ(hed::Triangulation& triang, std::ofstream& objfile);


#endif // CREATE_OBJ_H
