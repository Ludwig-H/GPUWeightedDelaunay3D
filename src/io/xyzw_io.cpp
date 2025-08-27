
#include "wd3d/common.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

namespace wd3d {

bool load_xyzw(const std::string& path, std::vector<Vertex>& out, ParseStats* stats) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    size_t n=0, skipped=0;
    out.clear();
    while (std::getline(in, line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream iss(line);
        double x,y,z,w;
        if (!(iss >> x >> y >> z >> w)) { skipped++; continue; }
        out.push_back(Vertex{x,y,z,w});
        n++;
    }
    if (stats) { stats->num_points = n; stats->num_skipped = skipped; }
    return true;
}

bool save_edges(const std::string& path, const EdgeList& E) {
    std::ofstream out(path);
    if (!out) return false;
    for (auto e : E.edges) out << e.u << " " << e.v << "\n";
    return true;
}

bool save_tets(const std::string& path, const RegularTriangulation& rt) {
    std::ofstream out(path);
    if (!out) return false;
    for (auto t : rt.tets) out << t.a << " " << t.b << " " << t.c << " " << t.d << "\n";
    return true;
}

} // namespace wd3d
