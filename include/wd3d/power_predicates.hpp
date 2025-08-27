
#pragma once
#include "common.hpp"

namespace wd3d {

inline long double orient3d_ld(const Vertex& a, const Vertex& b, const Vertex& c, const Vertex& d) {
    long double adx = (long double)a.x - d.x;
    long double ady = (long double)a.y - d.y;
    long double adz = (long double)a.z - d.z;
    long double bdx = (long double)b.x - d.x;
    long double bdy = (long double)b.y - d.y;
    long double bdz = (long double)b.z - d.z;
    long double cdx = (long double)c.x - d.x;
    long double cdy = (long double)c.y - d.y;
    long double cdz = (long double)c.z - d.z;
    return adx * (bdy * cdz - bdz * cdy)
         - ady * (bdx * cdz - bdz * cdx)
         + adz * (bdx * cdy - bdy * cdx);
}

inline long double in_power_sphere5_ld(const Vertex& a, const Vertex& b, const Vertex& c, const Vertex& d, const Vertex& e) {
    auto L = [](const Vertex& p)->long double {
        return (long double)p.x*(long double)p.x + (long double)p.y*(long double)p.y + (long double)p.z*(long double)p.z - (long double)p.w;
    };
    long double A[5][5] = {
        {(long double)a.x,(long double)a.y,(long double)a.z,L(a),1.0L},
        {(long double)b.x,(long double)b.y,(long double)b.z,L(b),1.0L},
        {(long double)c.x,(long double)c.y,(long double)c.z,L(c),1.0L},
        {(long double)d.x,(long double)d.y,(long double)d.z,L(d),1.0L},
        {(long double)e.x,(long double)e.y,(long double)e.z,L(e),1.0L}
    };
    int n=5, sign=1;
    long double M[5][5]; for(int i=0;i<n;i++) for(int j=0;j<n;j++) M[i][j]=A[i][j];
    for (int k=0;k<n;k++){
        int piv=k; long double maxa=fabsl(M[k][k]);
        for (int r=k+1;r<n;r++){ long double v=fabsl(M[r][k]); if (v>maxa){ maxa=v; piv=r; } }
        if (maxa==0) return 0;
        if (piv!=k){ for(int c=k;c<n;c++) std::swap(M[k][c],M[piv][c]); sign=-sign; }
        long double pivot=M[k][k];
        for (int r=k+1;r<n;r++){ long double f=M[r][k]/pivot; for(int c=k;c<n;c++) M[r][c]-=f*M[k][c]; }
    }
    long double det=sign; for(int i=0;i<n;i++) det*=M[i][i];
    return det;
}

inline int in_power_sphere_sign(const Vertex& a, const Vertex& b, const Vertex& c, const Vertex& d, const Vertex& e, long double eps=1e-18L) {
    long double o = orient3d_ld(a,b,c,d);
    if (fabsl(o) < eps) return 0;
    long double det5 = in_power_sphere5_ld(a,b,c,d,e);
    long double s = det5 * (o > 0 ? 1.0L : -1.0L);
    if (s > eps) return +1;
    if (s < -eps) return -1;
    return 0;
}

} // namespace wd3d
