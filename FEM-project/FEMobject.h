#ifndef FEMOBJECT_H
#define FEMOBJECT_H
#define _USE_MATH_DEFINES

#include <ttl/halfedge/HeTriang.h>
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>
#include <Eigen/Eigen>
//#include <Eigen/Eigenvalues>

#include <cmath>
#include <string>
#include "create_OBJ.h"


class FEMobject {
public:
    FEMobject();

    /// Methods
    void setProblemType(std::string problem);
    void circleMesh(int n, int m, double r);
    void squareMesh(int n, double d, Eigen::Vector2d Op);
    void prepareLaplace(int n = 8, int m = 12);
    void preparePoisson(int n = 8, int m = 12);
    void prepareHelmholtz(int n = 20);
    void prepareEigenvalue(int n = 20);

    void setSolution(Eigen::VectorXd zeta);
    void setSolution(int col=0);
    void solve();

    Eigen::SparseMatrix<double> massMat();      // M
    Eigen::SparseMatrix<double> stiffMat();     // A
    Eigen::VectorXd             loadVect();     // b
    Eigen::SparseMatrix<double> robinMat();     // R
    Eigen::VectorXd             robinVect();    // r

    /// Parameters
    double kappa(double x, double y);   // BC type
    double gN(double x, double y);      // Neumann BC
    double gD(double x, double y);      // Dirichlet BC
    double f(double x, double y);       // function f
    double uExact(double x, double y);  // exact solution

    /// Helper functions
    double triArea(hed::Node* N1, hed::Node* N2, hed::Node* N3);
    Eigen::Vector2<Eigen::Vector3d> gradients(
        Eigen::Vector3d x, Eigen::Vector3d y, double area);
    void localToGlobal();

    /// Getters
    double getError();  //computes the error
    int    getDoFs();   //returns the number of degrees of freedom

    /// Visualization. Create .obj file
    void visualization(std::string filename);
    void visualize_eigen(int i, std::vector<int> columns);

    /// Variables
    std::string problem_type;
    int num_nodes = 0;
    std::map<int, int> global_node_idx; // Node* id -> global id
    hed::Triangulation triang;          // mesh

    Eigen::MatrixXd zeta_mat;
};


#endif // FEMOBJECT_H
