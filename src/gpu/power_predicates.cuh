
#pragma once
#ifndef WD3D_NO_CUDA
#include <cuda_runtime.h>
#endif
#include <math.h>
struct float4h { float x,y,z,w; };
#ifndef WD3D_NO_CUDA
__device__ inline float orient3d_f(const float4h& a, const float4h& b, const float4h& c, const float4h& d) {
    float adx = a.x - d.x, ady = a.y - d.y, adz = a.z - d.z;
    float bdx = b.x - d.x, bdy = b.y - d.y, bdz = b.z - d.z;
    float cdx = c.x - d.x, cdy = c.y - d.y, cdz = c.z - d.z;
    return adx * (bdy*cdz - bdz*cdy) - ady*(bdx*cdz - bdz*cdx) + adz*(bdx*cdy - bdy*cdx);
}
__device__ inline float lift_f(const float4h& p) { return p.x*p.x + p.y*p.y + p.z*p.z - p.w; }
__device__ inline float det5(const float A[5][5]) {
    float m[5][5];
    #pragma unroll
    for(int i=0;i<5;i++) for(int j=0;j<5;j++) m[i][j]=A[i][j];
    int sign=1;
    for(int k=0;k<5;k++){
        int piv=k; float maxa=fabsf(m[k][k]);
        for(int r=k+1;r<5;r++){ float v=fabsf(m[r][k]); if(v>maxa){ maxa=v; piv=r; } }
        if(maxa==0) return 0;
        if(piv!=k){ for(int c=k;c<5;c++){ float tmp=m[k][c]; m[k][c]=m[piv][c]; m[piv][c]=tmp; } sign=-sign; }
        float pivot=m[k][k];
        for(int r=k+1;r<5;r++){ float f=m[r][k]/pivot; for(int c=k;c<5;c++) m[r][c]-=f*m[k][c]; }
    }
    float det=(float)sign; for(int i=0;i<5;i++) det*=m[i][i]; return det;
}
__device__ inline int in_power_sphere_sign_gpu(const float4h& a, const float4h& b, const float4h& c, const float4h& d, const float4h& e) {
    float A[5][5] = {
        {a.x,a.y,a.z,lift_f(a),1.0f},
        {b.x,b.y,b.z,lift_f(b),1.0f},
        {c.x,c.y,c.z,lift_f(c),1.0f},
        {d.x,d.y,d.z,lift_f(d),1.0f},
        {e.x,e.y,e.z,lift_f(e),1.0f}
    };
    float o = orient3d_f(a,b,c,d);
    if (o==0) return 0;
    float det = det5(A);
    float s = det * (o>0?1.0f:-1.0f);
    return (s>0) - (s<0);
}
#endif
