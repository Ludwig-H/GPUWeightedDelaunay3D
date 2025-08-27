
#include "wd3d/common.hpp"
#include "wd3d/regular_cpu_bw.hpp"
#include <iostream>

using namespace wd3d;

static void print_usage() {
    std::cout << "WeigthedDelaunay3D\n"
              << "Usage: WeigthedDelaunay3D --input file.xyzw [--out-edges edges.txt] [--out-tets tets.txt] [--cpu] [--verbose]\n";
}

int main(int argc, char** argv) {
    std::string in, out_edges, out_tets;
    bool force_cpu=false;
    int verbose=0;
    for (int i=1;i<argc;i++) {
        std::string a = argv[i];
        if (a=="--input" && i+1<argc) in=argv[++i];
        else if (a=="--out-edges" && i+1<argc) out_edges=argv[++i];
        else if (a=="--out-tets" && i+1<argc) out_tets=argv[++i];
        else if (a=="--cpu") force_cpu=true;
        else if (a=="--verbose") verbose=1;
        else if (a=="-h"||a=="--help") { print_usage(); return 0; }
    }
    if (in.empty()) { print_usage(); return 1; }

    std::vector<Vertex> P;
    ParseStats st;
    if (!load_xyzw(in, P, &st)) {
        std::cerr << "Failed to read " << in << "\n";
        return 2;
    }
    if (verbose) std::cerr << "Loaded " << st.num_points << " points (" << st.num_skipped << " skipped)\n";

    RegularTriangulation rt = regular_triangulation_cpu_bowyer(P, verbose);
    if (rt.tets.empty()) rt = regular_triangulation_cpu_bruteforce(P, verbose);

    if (!out_tets.empty()) {
        if (!save_tets(out_tets, rt)) { std::cerr << "Failed to save " << out_tets << "\n"; return 3; }
    }
    if (!out_edges.empty()) {
        EdgeList E = edges_from_triangulation(rt);
        if (!save_edges(out_edges, E)) { std::cerr << "Failed to save " << out_edges << "\n"; return 4; }
    }
    if (out_tets.empty() && out_edges.empty()) {
         EdgeList E = edges_from_triangulation(rt);
         for (auto e : E.edges) std::cout << e.u << " " << e.v << "\n";
    }
    return 0;
}
