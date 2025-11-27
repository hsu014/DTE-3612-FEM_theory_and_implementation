#ifndef TEST_EIGEN_H
#define TEST_EIGEN_H

#include <Eigen/Dense>
#include <iostream>

int test_eigen() {
    Eigen::MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    std::cout << m << std::endl;

    return 0;
}

#endif // TEST_EIGEN_H
