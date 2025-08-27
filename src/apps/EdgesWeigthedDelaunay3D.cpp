
#include "wd3d/common.hpp"
#include "wd3d/regular_cpu_bw.hpp"
#include <iostream>

using namespace wd3d;

static void usage() {
    std::cerr << "EdgesWeigthedDelaunay3D <in.xyzw> <out.edges>\n";
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(); return 1; }
    std::string in = argv[1];
    std::string out = argv[2];
    std::vector<Vertex> P;
    ParseStats st;
    if (!load_xyzw(in, P, &st)) { std::cerr << "Failed to read " << in << "\n"; return 2; }
    RegularTriangulation rt = regular_triangulation_cpu_bowyer(P, 0);
    if (rt.tets.empty()) rt = regular_triangulation_cpu_bruteforce(P, 0);
    EdgeList E = edges_from_triangulation(rt);
    if (!save_edges(out, E)) { std::cerr << "Failed to write " << out << "\n"; return 3; }
    return 0;
}
