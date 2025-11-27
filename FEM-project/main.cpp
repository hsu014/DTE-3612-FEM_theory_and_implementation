//#define DEBUG_TTL_CONSTR
//#define DEBUG_TTL_CONSTR_PLOT

#include <ttl/halfedge/HeTriang.h>
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>
#include <Eigen/Dense>
#include <iostream>

#include "FEMobject.h"
// #include "test_ttl.h"
// #include "test_eigen.h"

void singleLaplace(FEMobject& fem);
void singlePoisson(FEMobject& fem);
void singleHelmholtz(FEMobject& fem);

// void multipleLaplace(FEMobject& fem, int m);
// void multiplePoisson(FEMobject& fem, int m);
// void multipleHelmholtz(FEMobject& fem);


void singleLaplace(FEMobject &fem) {
    std::cout << "Laplace's equation\n";
    fem.prepareLaplace();
    fem.solve();
    fem.visualization("laplace.obj");
    std::cout << "DoF:   " << fem.getDoFs() << std::endl;
    std::cout << "Error: " << fem.getError() << std::endl;
}

void singlePoisson(FEMobject& fem) {
    std::cout << "Poisson's equation\n";
    fem.preparePoisson(8);
    fem.solve();
    fem.visualization("poisson.obj");
    std::cout << "DoF:   " << fem.getDoFs() << std::endl;
    std::cout << "Error: " << fem.getError() << std::endl;
}

void singleHelmholtz(FEMobject& fem) {
    std::cout << "Helmholtz equation\n";
    fem.prepareHelmholtz();
    fem.solve();
    fem.visualization("helmholtz.obj");
    std::cout << "DoF:   " << fem.getDoFs() << std::endl;
    std::cout << "Error: " << fem.getError() << std::endl;
}



int main() {
    FEMobject fem;
    int m = 10;

    // fem.squareMesh(10, 1.0, {-0.0, -0.5});
    // fem.visualization("test_square_mesh.obj");

    // std::cout << "Laplace's equation\n";
    // std::cout << "DoF\tError" << std::endl;
    // for (int i = 3; i < 20; i+=2) {
    //     fem.prepareLaplace(i, 10);
    //     fem.solve();
    //     std::cout << fem.getDoFs() << "\t" << fem.getError() << std::endl;
    //     fem.visualization(
    //         std::string("laplace_n") + std::to_string(i) +
    //         "_m" + std::to_string(m) + ".obj");
    // }
    // std::cout << "\n";

    // std::cout << "Poisson's equation\n";
    // std::cout << "DoF\tError" << std::endl;
    // for (int i = 3; i < 20; i+=2) {
    //     fem.preparePoisson(i, m);
    //     fem.solve();
    //     std::cout << fem.getDoFs() << "\t" << fem.getError() << std::endl;
    //     fem.visualization(
    //         std::string("poisson_n") + std::to_string(i) +
    //         "_m" + std::to_string(m) + ".obj");
    // }
    // std::cout << "\n";

    // std::cout << "Helmholtz equation\n";
    // std::cout << "DoF\tError" << std::endl;
    // for (int i = 8; i < 48; i+=4) {
    //     fem.prepareHelmholtz(i);
    //     fem.solve();
    //     std::cout << fem.getDoFs() << "\t" << fem.getError() << std::endl;
    //     fem.visualization(
    //         std::string("helmholtz_n") + std::to_string(i) + ".obj");
    // }
    // std::cout << "\n";

    std::vector<int> columns = {
        0, 1, 2, 3, 4
    };

    std::cout << "Eigenvalue problem\n";
    std::cout << "DoF" << std::endl;
    // int i = 7; i <= 13; i+=2
    for (int i = 11; i <= 11; i+=2) {
        fem.prepareEigenvalue(i);
        fem.solve();
        std::cout << fem.getDoFs() << std::endl;
        fem.visualize_eigen(i, columns);
    }
    std::cout << "\n";










    std::cout << "Eigenvalue problem\n";


    return 0;
}
