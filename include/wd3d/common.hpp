
#pragma once
#include <vector>
#include <array>
#include <cstdint>
#include <string>
#include <optional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <algorithm>

namespace wd3d {

struct Vertex { double x, y, z, w; inline double lift() const noexcept { return x*x + y*y + z*z - w; } };
struct Tetra { int a, b, c, d; };
struct Edge  { int u, v; };

struct RegularTriangulation { std::vector<Tetra> tets; };
struct EdgeList { std::vector<Edge> edges; };

struct ParseStats { size_t num_points=0, num_skipped=0; };

bool load_xyzw(const std::string& path, std::vector<Vertex>& out, ParseStats* stats=nullptr);
bool save_edges(const std::string& path, const EdgeList& E);
bool save_tets(const std::string& path, const RegularTriangulation& rt);

RegularTriangulation regular_triangulation_cpu_bruteforce(const std::vector<Vertex>& P, int verbose=0);
RegularTriangulation regular_triangulation_cpu_bowyer(const std::vector<Vertex>& P, int verbose=0);
EdgeList edges_from_triangulation(const RegularTriangulation& rt);

inline uint64_t pack_edge(int u, int v) { if (u>v) std::swap(u,v); return (uint64_t(uint32_t(u))<<32) | uint64_t(uint32_t(v)); }

} // namespace wd3d
