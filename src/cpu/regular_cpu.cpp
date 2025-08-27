
#include "wd3d/common.hpp"
#include "wd3d/power_predicates.hpp"
#include <iostream>

namespace wd3d {

RegularTriangulation regular_triangulation_cpu_bruteforce(const std::vector<Vertex>& P, int /*verbose*/) {
    RegularTriangulation rt;
    const int n = (int)P.size();
    if (n < 4) return rt;
    auto eps = 1e-18L;
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) for (int k=j+1;k<n;k++) for (int l=k+1;l<n;l++) {
        const Vertex &a=P[i], &b=P[j], &c=P[k], &d=P[l];
        long double o = orient3d_ld(a,b,c,d);
        if (fabsl(o) < eps) continue;
        bool empty = true;
        for (int m=0;m<n;m++) {
            if (m==i||m==j||m==k||m==l) continue;
            if (in_power_sphere_sign(a,b,c,d,P[m], eps) > 0) { empty=false; break; }
        }
        if (!empty) continue;
        if (o < 0) rt.tets.push_back({i,k,j,l}); else rt.tets.push_back({i,j,k,l});
    }
    return rt;
}

EdgeList edges_from_triangulation(const RegularTriangulation& rt) {
    std::unordered_set<uint64_t> S; S.reserve(rt.tets.size()*6);
    auto add=[&](int u,int v){ S.insert(pack_edge(u,v)); };
    for (auto t : rt.tets) {
        add(t.a,t.b); add(t.a,t.c); add(t.a,t.d);
        add(t.b,t.c); add(t.b,t.d); add(t.c,t.d);
    }
    EdgeList E; E.edges.reserve(S.size());
    for (auto code : S) {
        int u = int(code>>32), v=int(code & 0xffffffffu);
        E.edges.push_back({u,v});
    }
    std::sort(E.edges.begin(), E.edges.end(), [](auto&a, auto&b){ return a.u<b.u || (a.u==b.u && a.v<b.v); });
    return E;
}

} // namespace wd3d
