
#include "wd3d/common.hpp"
#include "wd3d/regular_cpu_bw.hpp"
#include <iostream>
using namespace wd3d;

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "need file\n"; return 1; }
    std::vector<Vertex> P;
    if (!load_xyzw(argv[1], P)) { std::cerr << "load failed\n"; return 2; }
    auto rt = regular_triangulation_cpu_bowyer(P, 0);
    if (rt.tets.empty()) rt = regular_triangulation_cpu_bruteforce(P, 0);
    auto edges = edges_from_triangulation(rt);
    if (edges.edges.empty()) { std::cerr << "no edges\n"; return 4; }
    std::cout << "OK " << P.size() << " points, " << rt.tets.size() << " tets, " << edges.edges.size() << " edges\n";
    return 0;
}
