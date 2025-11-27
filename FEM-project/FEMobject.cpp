#include "FEMobject.h"
//#include <iostream>

FEMobject::FEMobject() {

}

/// Methods
void FEMobject::setProblemType(std::string problem) {
    this->problem_type = problem;
}

void FEMobject::circleMesh(int n, int m, double r) {
    // n: number of rings
    // m: number of points in first ring
    // r: radius

    std::vector<hed::Node*>* nodes = new std::vector<hed::Node*>;
    hed::Node* p = new hed::Node(0, 0);
    nodes->push_back(p);

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m * (j+1); i++) {
            double dist = (r * (j+1) / n);
            double alpha = (i * 2 * M_PI) / (m * (j+1));

            hed::Node* p = new hed::Node(dist * cos(alpha), dist * sin(alpha));
            nodes->push_back(p);
        }
    }
    this->triang.createDelaunay(nodes->begin(), nodes->end());

    this->num_nodes = getDoFs();
    localToGlobal();
}

void FEMobject::squareMesh(int n, double d, Eigen::Vector2d Op) {
    // n: number of segments per side
    // d: side length

    std::vector<hed::Node*>* nodes = new std::vector<hed::Node*>;

    for (int j = 0; j <= n; j++) {
        for (int i = 0; i <= n; i++) {
            hed::Node* p = new hed::Node(
                Op(0) + (i*d/n), Op(1) + (j*d/n) );

            nodes->push_back(p);
        }
    }
    this->triang.createDelaunay(nodes->begin(), nodes->end());

    this->num_nodes = getDoFs();
    localToGlobal();
}

void FEMobject::prepareLaplace(int n, int m) {
    setProblemType("laplace");
    circleMesh(n, m, 1.0);
}

void FEMobject::preparePoisson(int n, int m) {
    setProblemType("poisson");
    circleMesh(n, m, 1.0);
}

void FEMobject::prepareHelmholtz(int n) {
    setProblemType("helmholtz");
    squareMesh(n, 1.0, {0.0, -0.5});
}

void FEMobject::prepareEigenvalue(int n) {
    setProblemType("eigenvalue");
    squareMesh(n, 1.0, {0.0, -0.5});
}

void FEMobject::setSolution(Eigen::VectorXd zeta){
    const list<hed::Node*>* nodelist = triang.getNodes();

    // loop over nodes
    for (hed::Node* node : *nodelist) {
        int id = this->global_node_idx[node->id()];
        node->init(node->x(), node->y(), zeta[id]);
    }
}

void FEMobject::setSolution(int col) {
    // For eigenvalue problem
    Eigen::VectorXd zeta = this->zeta_mat.col(col);
    const list<hed::Node*>* nodelist = triang.getNodes();

    // loop over nodes
    for (hed::Node* node : *nodelist) {
        int id = this->global_node_idx[node->id()];
        node->init(node->x(), node->y(), zeta[id]);
    }
}

/*
 * Solution
 * ℒ 𝜁 = ℓ         L * zeta = l
 *
 * Eigenvalue solution
 * ℒ 𝜁 = Λ 𝑀 𝜁     L * zeta = lambda * M * zeta
 *
 */

void FEMobject::solve() {
    if (this->problem_type == "laplace") {
        Eigen::SparseMatrix<double> A = stiffMat();
        Eigen::SparseMatrix<double> R = robinMat();

        Eigen::VectorXd r = robinVect();

        Eigen::SparseMatrix<double> L = A + R;
        Eigen::VectorXd l = r;

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        Eigen::VectorXd zeta(this->num_nodes);
        zeta = solver.compute(L).solve(l);

        setSolution(zeta);

        return;
    }
    else if (this->problem_type == "poisson") {
        Eigen::SparseMatrix<double> A = stiffMat();
        Eigen::SparseMatrix<double> R = robinMat();

        Eigen::VectorXd b = loadVect();
        Eigen::VectorXd r = robinVect();

        Eigen::SparseMatrix<double> L = A + R;
        Eigen::VectorXd l = b + r;

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        Eigen::VectorXd zeta(this->num_nodes);
        zeta = solver.compute(L).solve(l);

        setSolution(zeta);

        return;
    }
    else if (this->problem_type == "helmholtz") {
        double lambda = 81;

        Eigen::SparseMatrix<double> A = stiffMat();
        Eigen::SparseMatrix<double> R = robinMat();
        Eigen::SparseMatrix<double> M = massMat();

        Eigen::VectorXd r = robinVect();

        Eigen::SparseMatrix<double> L = A + R - lambda * M;
        Eigen::VectorXd l = r;

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        Eigen::VectorXd zeta(this->num_nodes);
        zeta = solver.compute(L).solve(l);

        setSolution(zeta);

        return;
    }
    else if (this->problem_type == "eigenvalue") {
        Eigen::SparseMatrix<double> M = massMat();
        Eigen::SparseMatrix<double> A = stiffMat();
        Eigen::SparseMatrix<double> R = robinMat();

        Eigen::SparseMatrix<double> L = A + R;

        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver;
        Eigen::MatrixXd Ldense = Eigen::MatrixXd(L);
        Eigen::MatrixXd Mdense = Eigen::MatrixXd(M);
        solver.compute(Mdense, Ldense);
        Eigen::MatrixXcd zetac = solver.eigenvectors(); // here
        Eigen::MatrixXd zeta = zetac.real();

        this->zeta_mat = zeta;

        return;
    }

    std::cout << "Problem type '" << this->problem_type << "' not defined\n";

}

Eigen::SparseMatrix<double> FEMobject::massMat() {     // M
    Eigen::SparseMatrix<double> M(this->num_nodes, this->num_nodes);

    // loop over triangles
    const list<hed::Edge*>& leading_edges = this->triang.getLeadingEdges();
    for (hed::Edge* edge : leading_edges) {
        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();

        // Global indices
        Eigen::Vector3<int> loc2glb = {
            this->global_node_idx[N1->id()],
            this->global_node_idx[N2->id()],
            this->global_node_idx[N3->id()]
        };

        // extract coordinates of the nodes
        Eigen::Vector3d x = {N1->x(), N2->x(), N3->x()};
        Eigen::Vector3d y = {N1->y(), N2->y(), N3->y()};

        double area = triArea(N1, N2, N3);

        // compute MK
        Eigen::Matrix3d mat {
            {2, 1, 1},
            {1, 2, 1},
            {1, 1, 2}
        };
        Eigen::Matrix3d MK = 1.0 / 12.0 * mat * area;

         // insert MK to the global mass matrix M
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                M.coeffRef(loc2glb[i], loc2glb[j]) += MK(i,j);
            }
        }
    }

    return M;
}

Eigen::SparseMatrix<double> FEMobject::stiffMat() {    // A
    Eigen::SparseMatrix<double> A(this->num_nodes, this->num_nodes);

    // loop over triangles
    const list<hed::Edge*>& leading_edges = this->triang.getLeadingEdges();
    for (hed::Edge* edge : leading_edges) {
        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();

        // Global indices
        Eigen::Vector3<int> loc2glb = {
           this->global_node_idx[N1->id()],
           this->global_node_idx[N2->id()],
           this->global_node_idx[N3->id()]
        };

        // extract coordinates of the nodes
        Eigen::Vector3d x = {N1->x(), N2->x(), N3->x()};
        Eigen::Vector3d y = {N1->y(), N2->y(), N3->y()};

        double area = triArea(N1, N2, N3);
        Eigen::Vector2<Eigen::Vector3d> gradphi = gradients(x, y, area);

        // compute AK
        Eigen::Vector3d b = gradphi(0);
        Eigen::Vector3d c = gradphi(1);
        Eigen::Matrix3d AK = (b*b.transpose() + c*c.transpose())*area;

        // insert AK to the global stiffness matrix A
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                A.coeffRef(loc2glb[i], loc2glb[j]) += AK(i,j);
            }
        }
    }

    return A;
}

Eigen::VectorXd FEMobject::loadVect() {    // b
    Eigen::VectorXd b = Eigen::VectorXd::Zero(this->num_nodes);

    // loop over triangles
    const list<hed::Edge*>& leading_edges = this->triang.getLeadingEdges();
    for (hed::Edge* edge : leading_edges) {
        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();

        // Global indices
        Eigen::Vector3<int> loc2glb = {
            this->global_node_idx[N1->id()],
            this->global_node_idx[N2->id()],
            this->global_node_idx[N3->id()]
        };

        // extract coordinates of the nodes
        Eigen::Vector3d x = {N1->x(), N2->x(), N3->x()};
        Eigen::Vector3d y = {N1->y(), N2->y(), N3->y()};
        double area = triArea(N1, N2, N3);

        // compute bK
        double xc = (x(0) + x(1) + x(2)) / 3.0;
        double yc = (y(0) + y(1) + y(2)) / 3.0;

        double F = f(xc, yc);
        Eigen::Vector3d vec = {1, 1, 1};
        Eigen::Vector3d bK = F * area / 3.0 * vec;

        //add bK to the global load vector b
        for (int j = 0; j < 3; j++) {
            b(loc2glb[j]) +=  bK(j);
        }
    }
    return b;
}

Eigen::SparseMatrix<double> FEMobject::robinMat() {    // R
    Eigen::SparseMatrix<double> R(this->num_nodes, this->num_nodes);

    // loop over boundary edges
    hed::Edge* edge = this->triang.getBoundaryEdge();
    hed::Dart b_dart(edge);
    list<hed::Dart> boundary;
    ttl::getBoundary(b_dart, boundary);

    list<hed::Dart>::iterator E;
    for (E = boundary.begin();  E != boundary.end(); E++) {
        hed::Node* N1 = E->getNode();
        hed::Node* N2 = E->getOppositeNode();

        // Global indices
        Eigen::Vector2<int> loc2glb = {
            this->global_node_idx[N1->id()],
            this->global_node_idx[N2->id()]
        };

        // extract coordinates of the nodes
        Eigen::Vector2d x = {N1->x(), N2->x()};
        Eigen::Vector2d y = {N1->y(), N2->y()};

        double len = sqrt((x(0) - x(1)) * (x(0) - x(1)) + (y(0) - y(1)) * (y(0) - y(1)));

        // find edge centroid
        double xc = (x(0) + x(1)) / 2.0;
        double yc = (y(0) + y(1)) / 2.0;
        double k = kappa(xc, yc); // compute the value of κ at centroid

        // compute RE
        Eigen::Matrix2d mat{
            {2, 1},
            {1, 2}
        };

        Eigen::Matrix2d RE = k / 6.0 * mat * len;
        //add RE to the global boundary matrix R
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                R.coeffRef(loc2glb[i], loc2glb[j]) += RE(i,j);
            }
        }
    }

    return R;
}

Eigen::VectorXd FEMobject::robinVect() {   // r
    Eigen::VectorXd r = Eigen::VectorXd::Zero(this->num_nodes);

    // loop over boundary edges
    hed::Edge* edge = this->triang.getBoundaryEdge();
    hed::Dart b_dart(edge);
    list<hed::Dart> boundary;
    ttl::getBoundary(b_dart, boundary);

    list<hed::Dart>::iterator E;
    for (E = boundary.begin();  E != boundary.end(); E++) {
        hed::Node* N1 = E->getNode();
        hed::Node* N2 = E->getOppositeNode();

        // Global indices
        Eigen::Vector2<int> loc2glb = {
            this->global_node_idx[N1->id()],
            this->global_node_idx[N2->id()]
        };

        // extract coordinates of the nodes
        Eigen::Vector2d x = {N1->x(), N2->x()};
        Eigen::Vector2d y = {N1->y(), N2->y()};

        double len = sqrt((x(0) - x(1)) * (x(0) - x(1)) + (y(0) - y(1)) * (y(0) - y(1)));

        // find edge centroid
        double xc = (x(0) + x(1)) / 2.0;
        double yc = (y(0) + y(1)) / 2.0;
        double tmp = kappa(xc, yc) * gD(xc, yc) + gN(xc, yc); // compute the value of the boundary conditions

        Eigen::Vector2d vec = {1, 1};
        Eigen::Vector2d rE = tmp * vec * len / 2.0;
        //add rE to the global boundary vector r
        for (int j = 0; j < 2; j++) {
            r(loc2glb[j]) += rE(j);
        }
    }

    return r;
}

/// Parameters
double FEMobject::kappa(double x, double y) {
    if (this->problem_type == "laplace") {
        return pow(10.0, 6.0);
    }
    else if (this->problem_type == "poisson") {
        return pow(10.0, 6.0);
    }
    else if (this->problem_type == "helmholtz") {
        if (x > 0.0) {
            return 0.0;
        }
        else {
            return pow(10.0, 6.0);
        }
    }
    else if (this->problem_type == "eigenvalue") {
        if (x > 0.0) {
            return 0.0;
        }
        else {
            return pow(10.0, 6.0);
        }
    }
    return 0.0;
}

double FEMobject::gN(double x, double y) {
    if (this->problem_type == "laplace") {
        return 0.0;
    }
    else if (this->problem_type == "poisson") {
        return 0.0;
    }
    else if (this->problem_type == "helmholtz") {
        return 0.0;
    }
    else if (this->problem_type == "eigenvalue") {
        return 0.0;
    }

    return 0.0;
}

double FEMobject::gD(double x, double y) {    
    if (this->problem_type == "laplace") {
        double phi = atan2(y, x);

        return cos(4.0 * phi);
    }
    else if (this->problem_type == "poisson") {
        return pow(y, 2.0) / 2.0;
    }
    else if (this->problem_type == "helmholtz") {
        return 0.25;
    }
    else if (this->problem_type == "eigenvalue") {
        return 0.0;
    }

    return 0.0;
}

double FEMobject::f(double x, double y) {
    if (this->problem_type == "laplace") {
        return 0.0;
    }
    else if (this->problem_type == "poisson") {
        return 1.0;
    }
    else if (this->problem_type == "helmholtz") {
        return 0.0;
    }
    else if (this->problem_type == "eigenvalue") {
        return 0.0;
    }

    return 0.0;
}

double FEMobject::uExact(double x, double y) {
    if (this->problem_type == "laplace") {
        double rho = sqrt(x * x + y * y);
        double phi = atan2(y, x);
        return (pow(rho, 4.0)) * cos(4.0 * phi);
    }
    else if (this->problem_type == "poisson") {
        return ((1.0 - pow(x, 2.0)) / 2.0);
    }
    else if (this->problem_type == "helmholtz") {
        // λ = ?
        double lambda = 81;

        return 0.25 * (cos(sqrt(lambda) * x) + tan(sqrt(lambda)) * sin(sqrt(lambda) * x));
    }
    else if (this->problem_type == "eigenvalue") {
        return 0.0;
    }

    return 0.0;
}

/// Helper functions
double FEMobject::triArea(hed::Node* N1, hed::Node* N2, hed::Node* N3) {
    Eigen::Vector3d a = {N2->x() - N1->x(), N2->y() - N1->y(), 0};
    Eigen::Vector3d b = {N3->x() - N1->x(), N3->y() - N1->y(), 0};
    double area = 0.5 * (b.cross(a)).norm();

    return area;
}

Eigen::Vector2<Eigen::Vector3d> FEMobject::gradients(Eigen::Vector3d x, Eigen::Vector3d y, double area) {
    Eigen::Vector3d b((y(1)-y(2))/(2*area), (y(2)-y(0))/(2*area), (y(0)-y(1))/(2*area));
    Eigen::Vector3d c((x(2)-x(1))/(2*area), (x(0)-x(2))/(2*area), (x(1)-x(0))/(2*area));
    Eigen::Vector2<Eigen::Vector3d> gradphi = {b, c};

    return gradphi;
}

void FEMobject::localToGlobal() {
    // Map Node*->id() to global id
    // std::map<node_id, global_id> global_node_idx
    this->global_node_idx.clear();

    const list<hed::Node*>* nodelist = this->triang.getNodes();
    int cur_id = 0;

    for (hed::Node* node : *nodelist) {
        this->global_node_idx.insert({node->id(), cur_id});
        cur_id++;
    }
}

/// Getters
double FEMobject::getError() {
    double error = 0.0;

    // loop over triangles
    const list<hed::Edge*>& leading_edges = this->triang.getLeadingEdges();
    for (hed::Edge* edge : leading_edges) {
        hed::Node* node_a = edge->getSourceNode();
        hed::Node* node_b = edge->getTargetNode();
        hed::Node* node_c = edge->getNextEdgeInFace()->getTargetNode();

        Eigen::Vector3d x = {node_a->x(), node_b->x(), node_c->x()};
        Eigen::Vector3d y = {node_a->y(), node_b->y(), node_c->y()};

        double area = triArea(node_a, node_b, node_c);

        double xc = (x(0) + x(1) + x(2)) / 3.0;
        double yc = (y(0) + y(1) + y(2)) / 3.0;
        double ubar = uExact(xc,yc);
        double uhbar = (node_a->z() + node_b->z() + node_c->z()) / 3.0;

        double eK = pow((ubar - uhbar), 2.0) *area;
        error += eK;
    }

    return sqrt(error);
}

int FEMobject::getDoFs() {
    return this->triang.getNodes()->size();
}

/// Visualization
void FEMobject::visualization(std::string filename) {
    std::ofstream obj_file(std::string("obj\\") + filename);
    createOBJ(this->triang, obj_file);
}

void FEMobject::visualize_eigen(int i, std::vector<int> columns) {
    if (this->problem_type == "eigenvalue") {
        for (int col : columns) {
            if (col < 0 || col >= (this->num_nodes)) continue;

            localToGlobal();
            setSolution(col);
            visualization(
                std::string("eigenvalue_n") + std::to_string(i) +
                "_c" + std::to_string(col) + ".obj");

        }
    }
}

