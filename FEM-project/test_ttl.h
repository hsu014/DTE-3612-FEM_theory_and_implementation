#ifndef TEST_TTL_H
#define TEST_TTL_H

#include <ttl/halfedge/HeTriang.h>
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>
//using namespace hed; // (to avoid using prefix hed::)

//#include <fstream>
#include <iostream>
#include <algorithm>

#include "create_OBJ.h"



// ------------------------------------------------------------------------------------------------
// Interpret two points as being coincident
inline bool eqPoints(hed::Node*& p1, hed::Node*& p2) {
    double dx = p1->x() - p2->x();
    double dy = p1->y() - p2->y();
    double dist2 = dx*dx + dy*dy;
    const double eps = 1.0e-12;
    if (dist2 < eps)
        return true;

    return false;
}


// ------------------------------------------------------------------------------------------------
// Lexicographically compare two points (2D)
inline bool ltLexPoint(const hed::Node* p1, const hed::Node* p2) {
    return (p1->x() < p2->x()) || (p1->x() == p2->x() && p1->y() < p2->y());
};


// =================================================================================================
// A main program using TTL and the half-edge data structure to create
// a Delaunay triangulation.
// The program also demonstrates how to use other function templates in TTL.
// =================================================================================================

int test_ttl() {

    // ===============================================================
    // CREATE A DELAUNAY TRIANGULATION FROM RANDOM POINTS IN THE PLANE
    // ===============================================================

    // Create random test data.
    int no_of_nodes = 100;
    std::vector<hed::Node*>* nodes = ttl_util::createRandomData<hed::Node>(no_of_nodes);

    // Sort the nodes lexicographically in the plane.
    // This is recommended since the triangulation algorithm will run much faster.
    // (ltLexPoint is defined above)
    std::sort(nodes->begin(), nodes->end(), ltLexPoint);

    // Remove coincident points to avoid degenerate triangles. (eqPoints is defined above)
    std::vector<hed::Node*>::iterator new_end = std::unique(nodes->begin(), nodes->end(), eqPoints);

    // Make the triangulation
    hed::Triangulation triang;
    triang.createDelaunay(nodes->begin(), new_end);

    // <... Print triangulation; see end of file ...>

    // ========================================================
    // SOME EXAMPLES USING TTL (Functions in namespace ttl)
    // ========================================================

    // Insert a new node in the Delaunay triangulation.
    // We need an arbitrary CCW (counterclockwise) dart for TTL.
    // Make the dart from the first edge in the list of leading edges.
    // ( Could also use Triangulation::createDart() )
    const list<hed::Edge*>& l_edges = triang.getLeadingEdges();
    hed::Edge* edge = *l_edges.begin();
    hed::Dart dart(edge);
    hed::Node* point = new hed::Node(0.3, 0.6, 0);
    ttl::insertNode<hed::TTLtraits>(dart, *point);

    // Locate a triangle in the triangulation containing the given point.
    // The given dart will be repositioned to that triangle while maintaining
    // its orientation (CCW or CW).
    // If the given point is outside the triangulation, the dart will be
    // positioned at a boundary edge.
    point->init(0.5, 0.5, 0);
    bool found = ttl::locateTriangle<hed::TTLtraits>(*point, dart);
    if (!found) {
        std::cout << "The given points is outside the triangulation" << std::endl;
        // (and the dart is positioned at a boundary edge)
        exit(-1);
    }

    // The degree (or valency) of a node V in a triangulation, is defined
    // as the number of edges incident with V.
    // Get the degree of the node associated with the dart.
    int degree = ttl::getDegreeOfNode(dart);
    std::cout << "Degree of node = " << degree << std::endl;

    // Check if the edge associated with the dart is at the boundary of the triangulation.
    if (ttl::isBoundaryEdge(dart))
        std::cout << "The edge is at the boundary" << std::endl;

    // Check if the node associated with the dart is at the boundary of the triangulation.
    if (ttl::isBoundaryNode(dart))
        std::cout << "The node is at the boundary" << std::endl;

    // Remove the node associated with the dart used above.
    ttl::removeNode<hed::TTLtraits>(dart);

    // Get the boundary of the triangulation represented as a list of darts.
    // Start with an arbitrary dart at the boundary.
    edge = triang.getBoundaryEdge();
    hed::Dart b_dart(edge);
    list<hed::Dart> boundary;
    ttl::getBoundary(b_dart,boundary);
    std::cout << "No. of edges on boundary = " << boundary.size() << std::endl;

    // Check if the triangulation is Delaunay
    // (This is not a TTL function)
    if (triang.checkDelaunay())
        std::cout << "Triangulation is Delaunay" << std::endl;
    else
        std::cout << "WARNING: Triangulation is not Delaunay" << std::endl;

    // Insert two nodes and then insert a constrained edge between them
    // (Note that this could also be implemented in the code for the data structure.
    //  Here we call ttl directly to demonstrate the generic concept.)
    hed::Dart d1 = triang.createDart(); // an arbitrary CCW dart
    ttl::insertNode<hed::TTLtraits>(d1, *new hed::Node(0.1, 0.25, 0));
    hed::Dart d2 = triang.createDart();
    ttl::insertNode<hed::TTLtraits>(d2, *new hed::Node(0.6, 0.85, 0));
    // (Note that d1 is not necessarily valid after having inserted d2 since insertion
    //  of d2 may affect d1. Here d2 is "far from" d1, so we are relatively safe).
    bool optimizeDelaunay = true; // optimizes to a constrained Delaunay triangulation
    dart = ttl::insertConstraint<hed::TTLtraits>(d1, d2, optimizeDelaunay);

    // Set the edge as constrained (fixed) such that it is not swapped when inserting nodes later
    dart.getEdge()->setConstrained();

    // Insert nodes and a constraint near that above to demonstrate fixed edges
    d1 = triang.createDart();
    ttl::insertNode<hed::TTLtraits>(d1, *new hed::Node(0.35, 0.56, 0));
    d2 = triang.createDart();
    ttl::insertNode<hed::TTLtraits>(d2, *new hed::Node(0.1, 0.9, 0));
    dart = ttl::insertConstraint<hed::TTLtraits>(d1, d2, optimizeDelaunay);
    dart.getEdge()->setConstrained();

    // Check if the boundary is convex (in the plane)
    if (ttl::convexBoundary<hed::TTLtraits>(b_dart))
        std::cout << "\nBoundary is convex:" << std::endl;

    // Print the nodes at the boundary
    list<hed::Dart>::const_iterator it;
    for (it = boundary.begin();  it != boundary.end(); ++it)
        std::cout << it->getNode()->x() << " " << it->getNode()->y() << '\n';

    // Print the triangulation (its edges)  to file
    // (for gnuplot or other line drawing programs)
    std::cout << "\nPrinting edges to file qweEdges.dat..." << std::endl;
    std::cout << "Plot triangulation with: gnuplot qwe.gnu" << std::endl;
    std::cout << "Printin triangulation to file qweEdges.obj" << std::endl;

    std::ofstream ofile("obj\\qweEdges.dat");
    triang.printEdges(ofile);

    std::ofstream obj_file("obj\\qweEdges.obj");
    createOBJ(triang, obj_file);

    return 0;
}


#endif // TEST_TTL_H
