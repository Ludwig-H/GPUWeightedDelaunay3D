
#include "wd3d/regular_cpu_bw.hpp"
#include "wd3d/power_predicates.hpp"
#include <queue>
#include <unordered_map>
#include <iostream>

namespace wd3d {

struct Tet { int v[4]; int nbr[4]; Tet(){ v[0]=v[1]=v[2]=v[3]=-1; nbr[0]=nbr[1]=nbr[2]=nbr[3]=-1; } };
static inline void sort3(int& a,int& b,int& c){ if(a>b) std::swap(a,b); if(b>c) std::swap(b,c); if(a>b) std::swap(a,b); }
struct FaceKey { int a,b,c; FaceKey(int i,int j,int k){ a=i;b=j;c=k; sort3(a,b,c);} bool operator==(const FaceKey&o) const noexcept {return a==o.a&&b==o.b&&c==o.c;} };
struct FaceKeyHash { size_t operator()(const FaceKey& f) const noexcept { return (size_t)f.a*73856093u ^ (size_t)f.b*19349663u ^ (size_t)f.c*83492791u; } };
static inline FaceKey face_of(const Tet& t, int fi){ int a=t.v[(fi+1)&3],b=t.v[(fi+2)&3],c=t.v[(fi+3)&3]; return FaceKey(a,b,c); }
static inline int in_power_sphere_sign_idx(const std::vector<Vertex>& P, const Tet& t, int e) {
    return in_power_sphere_sign(P[t.v[0]],P[t.v[1]],P[t.v[2]],P[t.v[3]],P[e]);
}

static void build_super_tetra(const std::vector<Vertex>& P, std::vector<Vertex>& Pext, Tet& super, int idx[4]) {
    double minx=1e300,miny=1e300,minz=1e300,maxx=-1e300,maxy=-1e300,maxz=-1e300, minw=0;
    for (auto& p:P){ minx=std::min(minx,p.x); miny=std::min(miny,p.y); minz=std::min(minz,p.z);
                     maxx=std::max(maxx,p.x); maxy=std::max(maxy,p.y); maxz=std::max(maxz,p.z);
                     minw=std::min(minw,p.w); }
    double span = std::max({maxx-minx, maxy-miny, maxz-minz, 1.0});
    double cx=0.5*(minx+maxx), cy=0.5*(miny+maxy), cz=0.5*(minz+maxz), R=10.0*span;
    Pext=P;
    idx[0]=Pext.size(); Pext.push_back({cx-R,cy,cz,minw});
    idx[1]=Pext.size(); Pext.push_back({cx+R,cy,cz,minw});
    idx[2]=Pext.size(); Pext.push_back({cx,cy-R,cz,minw});
    idx[3]=Pext.size(); Pext.push_back({cx,cy,cz+R,minw});
    super.v[0]=idx[0]; super.v[1]=idx[1]; super.v[2]=idx[2]; super.v[3]=idx[3];
}

RegularTriangulation regular_triangulation_cpu_bowyer(const std::vector<Vertex>& P, int verbose) {
    RegularTriangulation out;
    if (P.size() < 4) return out;

    std::vector<Vertex> Pext; Tet super; int sidx[4]; build_super_tetra(P, Pext, super, sidx);
    const int n = (int)P.size();
    std::vector<Tet> tets; tets.reserve(2*n); tets.push_back(super);
    std::unordered_map<FaceKey,std::pair<int,int>,FaceKeyHash> fmap;

    auto add_face = [&](int tid,int fidx){
        FaceKey fk=face_of(tets[tid],fidx);
        auto it=fmap.find(fk);
        if(it==fmap.end()) fmap.emplace(fk, std::make_pair(tid,fidx));
        else { int tid2=it->second.first, fidx2=it->second.second; tets[tid].nbr[fidx]=tid2; tets[tid2].nbr[fidx2]=tid; fmap.erase(it); }
    };
    auto add_tet = [&](int a,int b,int c,int d)->int{ Tet T; T.v[0]=a;T.v[1]=b;T.v[2]=c;T.v[3]=d; int id=(int)tets.size(); tets.push_back(T); add_face(id,0); add_face(id,1); add_face(id,2); add_face(id,3); return id; };
    auto remove_tet = [&](int tid){ tets[tid].v[0]=-1; };

    std::vector<int> order(n); for(int i=0;i<n;i++) order[i]=i;

    auto point_in_tet = [&](const Vertex& p, const Tet& T)->bool{
        const Vertex &A=Pext[T.v[0]],&B=Pext[T.v[1]],&C=Pext[T.v[2]],&D=Pext[T.v[3]];
        long double o1=orient3d_ld(A,B,C,p), o2=orient3d_ld(A,B,p,D), o3=orient3d_ld(A,p,C,D), o4=orient3d_ld(p,B,C,D);
        bool pos=(o1>0&&o2>0&&o3>0&&o4>0), neg=(o1<0&&o2<0&&o3<0&&o4<0); return pos||neg;
    };

    for (int ii=0; ii<n; ++ii) {
        int pidx=order[ii]; const Vertex& Pcur=Pext[pidx];
        int seed=-1;
        for (int tid=0; tid<(int)tets.size(); ++tid) { if (tets[tid].v[0]<0) continue; if (point_in_tet(Pcur,tets[tid])) { seed=tid; break; } }
        if (seed<0) seed=0;

        std::vector<char> in_cavity(tets.size(),0);
        std::queue<int> q; if (seed>=0){ q.push(seed); in_cavity[seed]=1; }
        while(!q.empty()){
            int tid=q.front(); q.pop();
            Tet& T=tets[tid];
            int sgn=in_power_sphere_sign_idx(Pext,T,pidx);
            if (sgn<=0) { in_cavity[tid]=0; continue; }
            for (int f=0; f<4; ++f){
                int nb=T.nbr[f]; if (nb<0 || nb>=(int)tets.size()) continue; if (tets[nb].v[0]<0) continue;
                if (!in_cavity[nb]){ in_cavity[nb]=1; q.push(nb); }
            }
        }
        bool any=false; for (auto c:in_cavity) if (c){ any=true; break; } if (!any && seed>=0) in_cavity[seed]=1;

        struct BFace{int a,b,c,tet_id,face_idx;}; std::vector<BFace> boundary; boundary.reserve(64);
        for (int tid=0; tid<(int)tets.size(); ++tid) if (in_cavity[tid]){
            Tet& T=tets[tid];
            for (int f=0; f<4; ++f){
                int nb=T.nbr[f]; bool is_boundary = (nb < 0) || (nb >= (int)tets.size()) || !in_cavity[nb];
                if (is_boundary){ int a=T.v[(f+1)&3], b=T.v[(f+2)&3], c=T.v[(f+3)&3]; boundary.push_back({a,b,c,tid,f}); }
            }
        }

        for (int tid=0; tid<(int)tets.size(); ++tid) if (in_cavity[tid]) remove_tet(tid);
        // Purge: invalider tous les voisins pour éviter d’utiliser des pointeurs périmés
        for (int tid=0; tid<(int)tets.size(); ++tid) if (tets[tid].v[0] >= 0) {
            tets[tid].nbr[0] = tets[tid].nbr[1] = tets[tid].nbr[2] = tets[tid].nbr[3] = -1;
        }

        fmap.clear();
        for (int tid=0; tid<(int)tets.size(); ++tid) if (tets[tid].v[0]>=0){ add_face(tid,0); add_face(tid,1); add_face(tid,2); add_face(tid,3); }

        for (auto& bf: boundary){
            Tet T; T.v[0]=bf.a; T.v[1]=bf.b; T.v[2]=bf.c; T.v[3]=pidx;
            long double o=orient3d_ld(Pext[T.v[0]],Pext[T.v[1]],Pext[T.v[2]],Pext[T.v[3]]);
            if (o<0) std::swap(T.v[0],T.v[1]);
            int id=(int)tets.size(); tets.push_back(T);
            add_face(id,0); add_face(id,1); add_face(id,2); add_face(id,3);
        }
    }

    std::vector<Tetra> final_tets; final_tets.reserve(tets.size());
    for (int tid=0; tid<(int)tets.size(); ++tid){
        Tet& T=tets[tid]; if (T.v[0]<0) continue;
        bool touch=false; for (int k=0;k<4;k++){ int v=T.v[k]; if (v==sidx[0]||v==sidx[1]||v==sidx[2]||v==sidx[3]) { touch=true; break; } }
        if (!touch) final_tets.push_back({T.v[0],T.v[1],T.v[2],T.v[3]});
    }
    out.tets=std::move(final_tets);
    return out;
}

} // namespace wd3d
