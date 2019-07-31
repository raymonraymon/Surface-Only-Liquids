//
//  Scenes.cpp
//  MultiTracker
//
//  Created by Fang Da on 15/1/26.
//
//

#include "Scenes.h"
#include <map>
#include "SimOptions.h"
#include "MeshIO.h"
#include "Sim.h"
#include "BPS3D.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Scene-specific initialization
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace
{
    class OrderComp
    {
    public:
        bool operator () (const std::pair<int, double> & o1, const std::pair<int, double> & o2) const { return o1.second < o2.second; }
    };
    
    // subdivide every triangle in the mesh into four triangles, by introducing three new vertices at edge midpoints
    // after subdivision, vs will be expanded (retaining all original vertices), while fs and ls will be replaced (no original faces remain)
    // if r is positive, new vertices will be projected onto the sphere centered at center, with a radius of r (for ICO sphere generation)
    void subdivide(const Vec3d & center, double r, std::vector<Vec3d> & vs, std::vector<Vec3i> & fs, std::vector<Vec2i> & ls)
    {
        std::vector<Vec3i> new_fs;
        std::vector<Vec2i> new_ls;
        std::map<std::pair<int, int>, int> new_vs_map;
        
        // map edge-midpoint to new vertices
        for (size_t j = 0; j < fs.size(); j++)
        {
            Vec3i & f = fs[j];
            int v0 = f[0];
            int v1 = f[1];
            int v2 = f[2];
            
            std::pair<int, int> p01 = (v0 < v1 ? std::pair<int, int>(v0, v1) : std::pair<int, int>(v1, v0));
            std::pair<int, int> p12 = (v1 < v2 ? std::pair<int, int>(v1, v2) : std::pair<int, int>(v2, v1));
            std::pair<int, int> p20 = (v2 < v0 ? std::pair<int, int>(v2, v0) : std::pair<int, int>(v0, v2));
            
            new_vs_map[p01] = 0;
            new_vs_map[p12] = 0;
            new_vs_map[p20] = 0;
        }
        
        // create the new vertices
        for (std::map<std::pair<int, int>, int>::iterator j = new_vs_map.begin(); j != new_vs_map.end(); j++)
        {
            j->second = vs.size();
            if (r > 0)
                vs.push_back(((vs[j->first.first] + vs[j->first.second]) / 2 - center).normalized() * r + center);
            else
                vs.push_back((vs[j->first.first] + vs[j->first.second]) / 2);
        }
        
        // triangulate
        for (size_t j = 0; j < fs.size(); j++)
        {
            Vec3i & f = fs[j];
            int v0 = f[0];
            int v1 = f[1];
            int v2 = f[2];
            
            std::pair<int, int> p01 = (v0 < v1 ? std::pair<int, int>(v0, v1) : std::pair<int, int>(v1, v0));
            std::pair<int, int> p12 = (v1 < v2 ? std::pair<int, int>(v1, v2) : std::pair<int, int>(v2, v1));
            std::pair<int, int> p20 = (v2 < v0 ? std::pair<int, int>(v2, v0) : std::pair<int, int>(v0, v2));
            int nv01 = new_vs_map[p01];
            int nv12 = new_vs_map[p12];
            int nv20 = new_vs_map[p20];
            
            Vec2i & l = ls[j];
            new_fs.push_back(Vec3i(v0, nv01, nv20));   new_ls.push_back(l);
            new_fs.push_back(Vec3i(nv01, v1, nv12));   new_ls.push_back(l);
            new_fs.push_back(Vec3i(nv20, nv12, v2));   new_ls.push_back(l);
            new_fs.push_back(Vec3i(nv12, nv20, nv01)); new_ls.push_back(l);
        }
        
        fs = new_fs;
        ls = new_ls;
    }
    
    void createIcoSphere(const Vec3d & center, double r, int subdivision, std::vector<Vec3d> & vs, std::vector<Vec3i> & fs, std::vector<Vec2i> & ls, const Vec2i & label = Vec2i(1, 0))
    {
        // create the initial icosahedron
        double phi = (1.0 + sqrt(5.0)) / 2.0;
        double len = Vec3d(0, 1, phi).norm();
        
        vs.resize(12);
        vs[ 0] = center + r * Vec3d(0,  1,  phi) / len;
        vs[ 1] = center + r * Vec3d(0, -1,  phi) / len;
        vs[ 2] = center + r * Vec3d(0,  1, -phi) / len;
        vs[ 3] = center + r * Vec3d(0, -1, -phi) / len;
        vs[ 4] = center + r * Vec3d( 1,  phi, 0) / len;
        vs[ 5] = center + r * Vec3d(-1,  phi, 0) / len;
        vs[ 6] = center + r * Vec3d( 1, -phi, 0) / len;
        vs[ 7] = center + r * Vec3d(-1, -phi, 0) / len;
        vs[ 8] = center + r * Vec3d( phi, 0,  1) / len;
        vs[ 9] = center + r * Vec3d( phi, 0, -1) / len;
        vs[10] = center + r * Vec3d(-phi, 0,  1) / len;
        vs[11] = center + r * Vec3d(-phi, 0, -1) / len;
        
        fs.push_back(Vec3i( 0,  1,  8));    ls.push_back(label);
        fs.push_back(Vec3i( 1,  0, 10));    ls.push_back(label);
        fs.push_back(Vec3i( 2,  3, 11));    ls.push_back(label);
        fs.push_back(Vec3i( 3,  2,  9));    ls.push_back(label);
        
        fs.push_back(Vec3i( 4,  5,  0));    ls.push_back(label);
        fs.push_back(Vec3i( 5,  4,  2));    ls.push_back(label);
        fs.push_back(Vec3i( 6,  7,  3));    ls.push_back(label);
        fs.push_back(Vec3i( 7,  6,  1));    ls.push_back(label);
        
        fs.push_back(Vec3i( 8,  9,  4));    ls.push_back(label);
        fs.push_back(Vec3i( 9,  8,  6));    ls.push_back(label);
        fs.push_back(Vec3i(10, 11,  7));    ls.push_back(label);
        fs.push_back(Vec3i(11, 10,  5));    ls.push_back(label);
        
        fs.push_back(Vec3i( 0,  8,  4));    ls.push_back(label);
        fs.push_back(Vec3i( 1,  6,  8));    ls.push_back(label);
        fs.push_back(Vec3i( 0,  5, 10));    ls.push_back(label);
        fs.push_back(Vec3i( 1, 10,  7));    ls.push_back(label);
        
        fs.push_back(Vec3i(11,  3,  7));    ls.push_back(label);
        fs.push_back(Vec3i( 5,  2, 11));    ls.push_back(label);
        fs.push_back(Vec3i( 6,  3,  9));    ls.push_back(label);
        fs.push_back(Vec3i( 9,  2,  4));    ls.push_back(label);
        
        for (int i = 0; i < subdivision; i++)
            subdivide(center, r, vs, fs, ls);
        
    }
    
    void createUVSphere(const Vec3d & center, double r, int N, int M, std::vector<Vec3d> & vs, std::vector<Vec3i> & fs, std::vector<Vec2i> & ls, const Vec2i & label = Vec2i(1, 0))
    {
        vs.push_back(center + r * Vec3d(0, 0, 1));
        for (int i = 1; i < N; i++)
            for (int j = 0; j < M; j++)
                vs.push_back(center + r * Vec3d(sin(i * M_PI / N) * cos(j * 2 * M_PI / M), sin(i * M_PI / N) * sin(j * 2 * M_PI / M), cos(i * M_PI / N)));
        vs.push_back(center + r * Vec3d(0, 0, -1));

        for (int j = 0; j < M; j++)
            fs.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), ls.push_back(label);
        for (int i = 0; i < N - 2; i++)
            for (int j = 0; j < M; j++)
                fs.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     ls.push_back(label),
                fs.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), ls.push_back(label);
        for (int j = 0; j < M; j++)
            fs.push_back(Vec3i((N - 1) * M + 1, (N - 2) * M + (j + 1) % M + 1, (N - 2) * M + j + 1)), ls.push_back(label);
    }
    
}

BPS3D * Scenes::scene(Sim * sim, const std::string & scene, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<Vec3d> & vels, std::vector<char> & solids)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");
    
    if (scene == "sphere")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r = Options::doubleValue("radius");
        createIcoSphere(Vec3d(0, 0, 0), r, N, v, f, l);
//        createUVSphere(Vec3d(0, 0, 0), r, N, M, v, f, l);
        
        for (size_t i = 0; i < v.size(); i++) v[i].z() *= Options::doubleValue("influx");
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(0, 0, 0);
        solids.resize(v.size(), false);
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "donut")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r1 = 1.0;
        double r2 = 0.3;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(Vec3d((r1 + r2 * cos(j * 2 * M_PI / M)) * cos(i * 2 * M_PI / N), (r1 + r2 * cos(j * 2 * M_PI / M)) * sin(i * 2 * M_PI / N), r2 * sin(j * 2 * M_PI / M)));
        
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j, ((i + 1) % N) * M + j, ((i + 1) % N) * M + (j + 1) % M)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i(((i + 1) % N) * M + (j + 1) % M, i * M + (j + 1) % M, i * M + j)),       l.push_back(Vec2i(1, 0));
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
//        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(v[i].y() / pow(v[i].x() * v[i].x() + v[i].y() * v[i].y(), 1.5), -v[i].x() / pow(v[i].x() * v[i].x() + v[i].y() * v[i].y(), 1.5), 0.0);
	for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(v[i].y(), -v[i].x(), 0);
        solids.resize(v.size(), false);
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);

    } else if (scene == "tet")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r = 1.0;
        double len = 1.0;
        v.push_back(Vec3d(-1, -1, -1) * r / len);
        v.push_back(Vec3d(-1,  1,  1) * r / len);
        v.push_back(Vec3d( 1, -1,  1) * r / len);
        v.push_back(Vec3d( 1,  1, -1) * r / len);
        
        f.push_back(Vec3i(0, 2, 1));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(0, 3, 2));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(0, 1, 3));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, 2, 3));    l.push_back(Vec2i(1, 0));
        
        for (int i = 0; i < N; i++)
            subdivide(Vec3d(0, 0, 0), 0, v, f, l);
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(Options::doubleValue("init-vel-x"), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z"));
        solids.resize(v.size(), false);

        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "cube")
    {
        int N = Options::intValue("mesh-size-n");
        
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r = 1.0;
        double len = 1.0;
        v.push_back(Vec3d(-1, -1, -1) * r / len);
        v.push_back(Vec3d( 1, -1, -1) * r / len);
        v.push_back(Vec3d(-1,  1, -1) * r / len);
        v.push_back(Vec3d(-1, -1,  1) * r / len);
        v.push_back(Vec3d(-1,  1,  1) * r / len);
        v.push_back(Vec3d( 1, -1,  1) * r / len);
        v.push_back(Vec3d( 1,  1, -1) * r / len);
        v.push_back(Vec3d( 1,  1,  1) * r / len);
        
        f.push_back(Vec3i(0, 2, 1));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(6, 1, 2));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(0, 3, 2));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(4, 2, 3));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(0, 1, 3));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(5, 3, 1));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, 6, 5));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(5, 6, 7));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(2, 4, 6));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(6, 4, 7));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(3, 5, 4));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(4, 5, 7));    l.push_back(Vec2i(1, 0));
        
        for (int i = 0; i < N; i++)
            subdivide(Vec3d(0, 0, 0), 0, v, f, l);
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(Options::doubleValue("init-vel-x"), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z"));
        solids.resize(v.size(), false);

        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "pilar")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r = 0.2;
        double h = 1.0;
        v.push_back(Vec3d(0, 0, h));
        for (int i = 0; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(Vec3d(r * cos(j * 2 * M_PI / M), r * sin(j * 2 * M_PI / M), h * (1 - i * 2.0 / N)));
        v.push_back(Vec3d(0, 0, -h));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i((N + 1) * M + 1, N * M + (j + 1) % M + 1, N * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(0, 0, 0);
        
        solids.resize(v.size(), false);
        solids[0] = true;
        for (size_t j = 0; j < M; j++)
            solids[1 + j] = true,
            solids[1 + N * M + j] = true;
        solids[(N + 1) * M + 1] = true;
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "dome")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double alpha_end = M_PI / 2;
        double r = Options::doubleValue("radius");
        v.push_back(Vec3d(0, 0, 1));
        for (int i = 1; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(Vec3d(sin(i * alpha_end / N) * cos((j - 0.5 * i) * 2 * M_PI / M), sin(i * alpha_end / N) * sin((j - 0.5 * i) * 2 * M_PI / M), cos(i * alpha_end / N)) * r);
        v.push_back(Vec3d(0, 0, cos(alpha_end) * r));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(N * M + 1, (N - 1) * M + (j + 1) % M + 1, (N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(0, 0, 0);
        
        solids.resize(v.size(), false);
        for (size_t j = 0; j < M; j++)
            solids[(N - 1) * M + j + 1] = true;
        solids[N * M + 1] = true;
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "collision")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r = Options::doubleValue("radius");
        double d = Options::doubleValue("distance");
        double od = Options::doubleValue("offcenter-distance");
        
//        createIcoSphere(Vec3d(-1.1, 0, 0), r, N, v, f, l);
        createUVSphere(Vec3d(-d / 2, 0, -od / 2), r, N, M, v, f, l);
        
        size_t v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
        
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(Vec3d( std::abs(Options::doubleValue("init-vel-x")), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z")));
        
        v.clear();
        f.clear();
        l.clear();
        
//        createIcoSphere(Vec3d( 1.1, 0, 0), r, N, v, f, l);
        createUVSphere(Vec3d(d / 2, 0, od / 2), r, N, M, v, f, l);

        v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
        
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(Vec3d(-std::abs(Options::doubleValue("init-vel-x")), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z")));
        
        solids.resize(vs.size(), false);
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "dripping")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double alpha_end = M_PI / 3;
        double r = Options::doubleValue("radius");
        v.push_back(Vec3d(0, 0, -r));
        for (int i = 1; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(Vec3d(sin(i * alpha_end / N) * cos((j - 0.0 * i) * 2 * M_PI / M), -sin(i * alpha_end / N) * sin((j - 0.0 * i) * 2 * M_PI / M), -cos(i * alpha_end / N)) * r);
        v.push_back(Vec3d(0, 0, -cos(alpha_end) * r));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(N * M + 1, (N - 1) * M + (j + 1) % M + 1, (N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(0, 0, 0);
        
        solids.resize(v.size(), false);
        for (size_t j = 0; j < M; j++)
            solids[(N - 1) * M + j + 1] = true;
        solids[N * M + 1] = true;
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "transport")
    {
        double alpha_end = M_PI / 3;
        double r1 = Options::doubleValue("radius");
        double r2 = Options::doubleValue("radius") / 2;
        int M1 = M;
        int M2 = (int)(M1 * r2 / r1 + 0.5);
        double len = 12.0;
        int Nl = (int)(M1 * len / (2 * M_PI * r1) + 0.5);
        
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        // half dome on one side
        v.push_back(Vec3d(0, len / 2, (1 - cos(alpha_end)) * r1));
        for (int i = 1; i <= N; i++)
            for (int j = 0; j <= M1 / 2; j++)
                v.push_back(Vec3d(sin(i * alpha_end / N) * cos((j - 0.0 * i) * 2 * M_PI / M1), sin(i * alpha_end / N) * sin((j - 0.0 * i) * 2 * M_PI / M1), cos(i * alpha_end / N)) * r1 + Vec3d(0, len / 2, -cos(alpha_end) * r1));
        v.push_back(Vec3d(0, len / 2, 0));

        // half dome on the other side
        v.push_back(Vec3d(0, -len / 2, (1 - cos(alpha_end)) * r2));
        for (int i = 1; i <= N; i++)
            for (int j = M2 / 2; j <= M2; j++)
                v.push_back(Vec3d(sin(i * alpha_end / N) * cos((j - 0.0 * i) * 2 * M_PI / M2), sin(i * alpha_end / N) * sin((j - 0.0 * i) * 2 * M_PI / M2), cos(i * alpha_end / N)) * r2 + Vec3d(0, -len / 2, -cos(alpha_end) * r2));
        v.push_back(Vec3d(0, -len / 2, 0));
        
        // connecting the two domes in the middle
        for (int j = 1; j < Nl; j++)
            v.push_back((1 - (double)j / Nl) * Vec3d(0, len / 2, (1 - cos(alpha_end)) * r1) + ((double)j / Nl) * Vec3d(0, -len / 2, (1 - cos(alpha_end)) * r2));
        for (int i = 1; i <= N; i++)
        {
            for (int j = 1; j < Nl; j++)
                v.push_back((1 - (double)j / Nl) * Vec3d( sin(i * alpha_end / N) * r1, len / 2, (cos(i * alpha_end / N) - cos(alpha_end)) * r1) + ((double)j / Nl) * Vec3d( sin(i * alpha_end / N) * r2, -len / 2, (cos(i * alpha_end / N) - cos(alpha_end)) * r2));
            for (int j = 1; j < Nl; j++)
                v.push_back((1 - (double)j / Nl) * Vec3d(-sin(i * alpha_end / N) * r1, len / 2, (cos(i * alpha_end / N) - cos(alpha_end)) * r1) + ((double)j / Nl) * Vec3d(-sin(i * alpha_end / N) * r2, -len / 2, (cos(i * alpha_end / N) - cos(alpha_end)) * r2));
        }
        for (int j = 1; j < Nl; j++)
            v.push_back((1 - (double)j / Nl) * Vec3d(0, len / 2, 0) + ((double)j / Nl) * Vec3d(0, -len / 2, 0));

        // half dome on one side
        int m1 = M1 / 2 + 1;
        for (int j = 0; j < M1 / 2; j++)
            f.push_back(Vec3i(0, j + 1, j + 2)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M1 / 2; j++)
                f.push_back(Vec3i((i + 0) * m1 + j + 1, (i + 1) * m1 + j + 1, (i + 1) * m1 + j + 2)), l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * m1 + j + 2, (i + 0) * m1 + j + 2, (i + 0) * m1 + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M1 / 2; j++)
            f.push_back(Vec3i(N * m1 + 1, (N - 1) * m1 + j + 2, (N - 1) * m1 + j + 1)), l.push_back(Vec2i(1, 0));
        
        // half dome on the other side
        int m2 = M2 / 2 + 1;
        int n1 = N * m1 + 2;
        for (int j = 0; j < M2 / 2; j++)
            f.push_back(Vec3i(n1, n1 + j + 1, n1 + j + 2)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M2 / 2; j++)
                f.push_back(Vec3i(n1 + (i + 0) * m2 + j + 1, n1 + (i + 1) * m2 + j + 1, n1 + (i + 1) * m2 + j + 2)), l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i(n1 + (i + 1) * m2 + j + 2, n1 + (i + 0) * m2 + j + 2, n1 + (i + 0) * m2 + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M2 / 2; j++)
            f.push_back(Vec3i(n1 + N * m2 + 1, n1 + (N - 1) * m2 + j + 2, n1 + (N - 1) * m2 + j + 1)), l.push_back(Vec2i(1, 0));
        
        // connecting the two domes in the middle
        int n2 = n1 + N * m2 + 2;
        int ml = (Nl - 1) * 2;
        
        f.push_back(Vec3i(0, n2, n2 + (Nl - 1))),               l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(n2 + (Nl - 1), 1, 0)),                l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(0, M1 / 2 + 1, n2 + (Nl - 1) * 2)),   l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(n2 + (Nl - 1) * 2, n2, 0)),           l.push_back(Vec2i(1, 0));
        for (int j = 1; j < Nl - 1; j++)
            f.push_back(Vec3i(n2 + (j - 1), n2 + (j - 0), n2 + (Nl - 1) + (j - 0))),                    l.push_back(Vec2i(1, 0)),
            f.push_back(Vec3i(n2 + (Nl - 1) + (j - 0), n2 + (Nl - 1) + (j - 1), n2 + (j - 1))),         l.push_back(Vec2i(1, 0)),
            f.push_back(Vec3i(n2 + (j - 1), n2 + (Nl - 1) * 2 + (j - 1), n2 + (Nl - 1) * 2 + (j - 0))), l.push_back(Vec2i(1, 0)),
            f.push_back(Vec3i(n2 + (Nl - 1) * 2 + (j - 0), n2 + (j - 0), n2 + (j - 1))),                l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(n2 + Nl - 2, n1, n1 + M2 / 2 + 1)),                     l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(n1 + M2 / 2 + 1, n2 + (Nl - 1) + Nl - 2, n2 + Nl - 2)), l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(n2 + Nl - 2, n2 + (Nl - 1) * 2 + Nl - 2, n1 + 1)),      l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(n1 + 1, n1, n2 + Nl - 2)),                              l.push_back(Vec2i(1, 0));

        for (int i = 0; i < N - 1; i++)
        {
            f.push_back(Vec3i(0, 1, 1) * (n2 + (Nl - 1) + i * ml           ) + Vec3i((i + 0) * m1 + 1, 0, 2 * (Nl - 1) + 0)),                                  l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i(1, 0, 0) * (n2 + (Nl - 1) + i * ml           ) + Vec3i(2 * (Nl - 1) + 0, (i + 1) * m1 + 1, (i + 0) * m1 + 1)),                   l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i(0, 0, 1) * (n2 + (Nl - 1) + i * ml + (Nl - 1)) + Vec3i((i + 0) * m1 + M1 / 2 + 1, (i + 1) * m1 + M1 / 2 + 1, 2 * (Nl - 1) + 0)), l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i(1, 1, 0) * (n2 + (Nl - 1) + i * ml + (Nl - 1)) + Vec3i(2 * (Nl - 1) + 0, 0, (i + 0) * m1 + M1 / 2 + 1)),                         l.push_back(Vec2i(1, 0));
            for (int j = 1; j < Nl - 1; j++)
                f.push_back(Vec3i(1, 1, 1) * (n2 + (Nl - 1) + i * ml           ) + Vec3i(j - 1, j - 0, 2 * (Nl - 1) + j - 0)),                 l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i(1, 1, 1) * (n2 + (Nl - 1) + i * ml           ) + Vec3i(2 * (Nl - 1) + j - 0, 2 * (Nl - 1) + j - 1, j - 1)),  l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i(1, 1, 1) * (n2 + (Nl - 1) + i * ml + (Nl - 1)) + Vec3i(j - 1, 2 * (Nl - 1) + j - 1, 2 * (Nl - 1) + j - 0)),  l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i(1, 1, 1) * (n2 + (Nl - 1) + i * ml + (Nl - 1)) + Vec3i(2 * (Nl - 1) + j - 0, j - 0, j - 1)),                 l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i(1, 0, 0) * (n2 + (Nl - 1) + i * ml           ) + Vec3i(Nl - 2, n1 + (i + 0) * m2 + M2 / 2 + 1, n1 + (i + 1) * m2 + M2 / 2 + 1)), l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i(0, 1, 1) * (n2 + (Nl - 1) + i * ml           ) + Vec3i(n1 + (i + 1) * m2 + M2 / 2 + 1, 2 * (Nl - 1) + Nl - 2, Nl - 2)),          l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i(1, 1, 0) * (n2 + (Nl - 1) + i * ml + (Nl - 1)) + Vec3i(Nl - 2, 2 * (Nl - 1) + Nl - 2, n1 + (i + 1) * m2 + 1)),                   l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i(0, 0, 1) * (n2 + (Nl - 1) + i * ml + (Nl - 1)) + Vec3i(n1 + (i + 1) * m2 + 1, n1 + (i + 0) * m2 + 1, Nl - 2)),                   l.push_back(Vec2i(1, 0));
        }
        
        f.push_back(Vec3i(0, 1, 1) * (n2 + (Nl - 1) + N * ml) + Vec3i(1, 0, 0) * (N * m1 + 1) + Vec3i(-m1 + 0, -ml + 0, 0)),                 l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, 0, 0) * (n2 + (Nl - 1) + N * ml) + Vec3i(0, 1, 1) * (N * m1 + 1) + Vec3i(0, 0, -m1 + 0)),                       l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(0, 0, 1) * (n2 + (Nl - 1) + N * ml) + Vec3i(1, 1, 0) * (N * m1 + 1) + Vec3i(-m1 + M1 / 2, 0, 0)),                  l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, 1, 0) * (n2 + (Nl - 1) + N * ml) + Vec3i(0, 0, 1) * (N * m1 + 1) + Vec3i(0, -ml + (Nl - 1) + 0, -m1 + M1 / 2)), l.push_back(Vec2i(1, 0));
        for (int j = 1; j < Nl - 1; j++)
            f.push_back(Vec3i(1, 1, 1) * (n2 + (Nl - 1) + N * ml) + Vec3i(-ml + (j - 1), -ml + (j - 0), (j - 0))),                       l.push_back(Vec2i(1, 0)),
            f.push_back(Vec3i(1, 1, 1) * (n2 + (Nl - 1) + N * ml) + Vec3i((j - 0), (j - 1), -ml + (j - 1))),                             l.push_back(Vec2i(1, 0)),
            f.push_back(Vec3i(1, 1, 1) * (n2 + (Nl - 1) + N * ml) + Vec3i(-ml + (Nl - 1) + (j - 1), (j - 1), (j - 0))),                  l.push_back(Vec2i(1, 0)),
            f.push_back(Vec3i(1, 1, 1) * (n2 + (Nl - 1) + N * ml) + Vec3i((j - 0), -ml + (Nl - 1) + (j - 0), -ml + (Nl - 1) + (j - 1))), l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, 0, 0) * (n2 + (Nl - 1) + N * ml) + Vec3i(0, 1, 1) * (n1 + N * m2 + 1) + Vec3i(-ml + Nl - 2, -m2 + M2 / 2, 0)),       l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(0, 1, 1) * (n2 + (Nl - 1) + N * ml) + Vec3i(1, 0, 0) * (n1 + N * m2 + 1) + Vec3i(0, Nl - 2, -ml + Nl - 2)),             l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, 1, 0) * (n2 + (Nl - 1) + N * ml) + Vec3i(0, 0, 1) * (n1 + N * m2 + 1) + Vec3i(-ml + (Nl - 1) + Nl - 2, Nl - 2, 0)),  l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(0, 0, 1) * (n2 + (Nl - 1) + N * ml) + Vec3i(1, 1, 0) * (n1 + N * m2 + 1) + Vec3i(0, -m2 + 0, -ml + (Nl - 1) + Nl - 2)), l.push_back(Vec2i(1, 0));
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(0, 0, 0);
        
        solids.resize(v.size(), false);
        for (size_t j = 0; j <= M1 / 2; j++)
            solids[(N - 1) * m1 + j + 1] = true;
        solids[N * m1 + 1] = true;
        for (size_t j = 0; j <= M2 / 2; j++)
            solids[n1 + (N - 1) * m2 + j + 1] = true;
        solids[n1 + N * m2 + 1] = true;
        for (size_t j = 0; j < Nl - 1; j++)
            solids[n2 + (Nl - 1) + (N - 1) * ml + j] = true,
            solids[n2 + (Nl - 1) + (N - 1) * ml + (Nl - 1) + j] = true,
            solids[n2 + (Nl - 1) + N * ml + j] = true;
        
        BPS3D * bps = new BPS3D(sim, vs, fs, ls, vels, solids);
        
        return bps;
        
    } else if (scene == "floorsplash")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r = Options::doubleValue("radius");
        double h = 1.1;
//        createIcoSphere(Vec3d(0, 0, h), r, N, v, f, l);
        createUVSphere(Vec3d(0, 0, h), r, N, M, v, f, l);
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        vels.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) vels[i] = Vec3d(Options::doubleValue("init-vel-x"), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z"));
        solids.resize(v.size(), false);
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "jet")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double nozzle_theta = Options::doubleValue("nozzle-angle") * M_PI / 180;    // direction of the nozzle
        Vec3d nozzle_center(0, 0, 0);
        Mat3d R;
        double c = cos(nozzle_theta), s = sin(nozzle_theta);
        R << c, 0, -s, 0, 1, 0, s, 0, c;
        
        double alpha_end = M_PI / 2;
        double r = Options::doubleValue("radius");
        
        double influx = Options::doubleValue("influx");
        
        v.push_back(nozzle_center + R * Vec3d(0, 0, -r));
        for (int i = 1; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(nozzle_center + R * Vec3d(sin(i * alpha_end / N) * cos((j - 0.0 * i) * 2 * M_PI / M), -sin(i * alpha_end / N) * sin((j - 0.0 * i) * 2 * M_PI / M), -cos(i * alpha_end / N)) * r);
        v.push_back(nozzle_center + R * Vec3d(0, 0, 0));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(N * M + 1, (N - 1) * M + (j + 1) % M + 1, (N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
        
        solids.resize(v.size(), false);
        for (size_t j = 0; j < M; j++)
            solids[(N - 1) * M + j + 1] = true;
        solids[N * M + 1] = true;
        
        Vec3d n = R * Vec3d(0, 0, -1);
        vels.reserve(v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(solids[i] ? Vec3d::Zero() : Vec3d(n * influx));
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);

    } else if (scene == "jetchain")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double nozzle_theta = Options::doubleValue("nozzle-angle") * M_PI / 180;    // direction of the nozzle
        Vec3d nozzle_center(-Options::doubleValue("distance"), 0, 0);
        Mat3d R;
        double c = cos(nozzle_theta), s = sin(nozzle_theta);
        R << c, 0, -s, 0, 1, 0, s, 0, c;
        
        double alpha_end = M_PI / 2;
        double r = Options::doubleValue("radius");
        
        double influx = Options::doubleValue("influx");
        
        int arrangement = Options::intValue("nozzle-arrangement");
        Mat3d arrangementR;
        if (arrangement == 0)       // nozzles in left and right: original geometry
            arrangementR << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        else if (arrangement == 1)  // nozzles in down and up: permute (x,y,z) to (z,x,y)
            arrangementR << 0, 0, 1, 1, 0, 0, 0, 1, 0;
        else if (arrangement == 2)  // nozzles in front and back: permute (x,y,z) to (y,z,x)
            arrangementR << 0, 1, 0, 0, 0, 1, 1, 0, 0;
        
        // left nozzle
        v.push_back(nozzle_center + R * Vec3d(0, 0, -r));
        for (int i = 1; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(nozzle_center + R * Vec3d(sin(i * alpha_end / N) * cos((j - 0.0 * i) * 2 * M_PI / M), -sin(i * alpha_end / N) * sin((j - 0.0 * i) * 2 * M_PI / M), -cos(i * alpha_end / N)) * r);
        v.push_back(nozzle_center + R * Vec3d(0, 0, 0));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(N * M + 1, (N - 1) * M + (j + 1) % M + 1, (N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        for (size_t i = 0; i < v.size(); i++) v[i] = arrangementR * v[i];
        
        size_t v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));

        solids.reserve(solids.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) solids.push_back(false);
        for (size_t j = 0; j < M; j++) solids[v0 + (N - 1) * M + j + 1] = true;
        solids[v0 + N * M + 1] = true;

        Vec3d n = arrangementR * R * Vec3d(0, 0, -1);
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(solids[v0 + i] ? Vec3d::Zero() : Vec3d(n * influx));
        
        // right nozzle
        r *= Options::doubleValue("nozzle-radius-ratio");
        v.clear();
        f.clear();
        l.clear();
        
        nozzle_center = -nozzle_center;
        R << c, 0, s, 0, 1, 0, -s, 0, c;
        v.push_back(nozzle_center + R * Vec3d(0, 0, -r));
        for (int i = 1; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(nozzle_center + R * Vec3d(sin(i * alpha_end / N) * cos((j - 0.0 * i) * 2 * M_PI / M), -sin(i * alpha_end / N) * sin((j - 0.0 * i) * 2 * M_PI / M), -cos(i * alpha_end / N)) * r);
        v.push_back(nozzle_center + R * Vec3d(0, 0, 0));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(N * M + 1, (N - 1) * M + (j + 1) % M + 1, (N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        for (size_t i = 0; i < v.size(); i++) v[i] = arrangementR * v[i];

        v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));

        solids.reserve(solids.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) solids.push_back(false);
        for (size_t j = 0; j < M; j++) solids[v0 + (N - 1) * M + j + 1] = true;
        solids[v0 + N * M + 1] = true;
        
        n = arrangementR * R * Vec3d(0, 0, -1);
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(solids[v0 + i] ? Vec3d::Zero() : Vec3d(n * influx));
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "transportinject")
    {
        // load the equilibrium configuration
        BPS3D * bps = new BPS3D(sim, std::vector<LosTopos::Vec3d>(), std::vector<LosTopos::Vec3st>(), std::vector<LosTopos::Vec2i>(), std::vector<Vec3d>(), std::vector<char>());
        MeshIO::load(*bps, "init/transportinject.rec");
        
        vs = bps->m_st->pm_positions;
        fs = bps->m_st->m_mesh.m_tris;
        ls = bps->m_st->m_mesh.m_triangle_labels;
        vels = bps->vel();
        solids = bps->solid_labels();
        
        for (size_t i = 0; i < vels.size(); i++)    // clear the velocity: the loaded velocity is not desired.
            vels[i].setZero();
        
        delete bps;
        
        bps = new BPS3D(sim, vs, fs, ls, vels, solids);

        return bps;
        
    } else if (scene == "pool")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r_pool = 40.0;
        double d_pool = 40.0;
        
        // pool
        v.push_back(Vec3d(0, 0, 0));
        for (int i = 1; i < N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(Vec3d(r_pool * i / N * cos(j * 2 * M_PI / M), r_pool * i / N * sin(j * 2 * M_PI / M), 0));
        for (int i = 0; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(Vec3d(r_pool * cos(j * 2 * M_PI / M), r_pool * sin(j * 2 * M_PI / M), -d_pool * (double)i / N));
        v.push_back(Vec3d(0, 0, -d_pool));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N * 2 - 1; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i((2 * N) * M + 1, (2 * N - 1) * M + (j + 1) % M + 1, (2 * N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        // introduce a perturbation
        double perturbation = 1.5;
        v.front() += Vec3d(0, 0, 1) * perturbation;
        
        size_t v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
        
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(Vec3d(0, 0, 0));
        
        solids.reserve(solids.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) solids.push_back(true);
        for (int i = 0; i < N - 1; i++) for (int j = 0; j < M; j++) solids[v0 + i * M + j + 1] = false;
        solids[v0] = false;
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);

    } else if (scene == "crownsplash")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
//        double r_pool = 40.0;
//        double d_pool = 40.0;
        double r_droplet = Options::doubleValue("radius");
        double z_droplet = Options::doubleValue("radius") - Options::doubleValue("init-vel-z") * Options::doubleValue("time-step") + 1e-6;
        
        // pool
//        v.push_back(Vec3d(0, 0, 0));
//        for (int i = 1; i < N; i++)
//            for (int j = 0; j < M; j++)
//                v.push_back(Vec3d(r_pool * i / N * cos(j * 2 * M_PI / M), r_pool * i / N * sin(j * 2 * M_PI / M), 0));
//        for (int i = 0; i <= N; i++)
//            for (int j = 0; j < M; j++)
//                v.push_back(Vec3d(r_pool * cos(j * 2 * M_PI / M), r_pool * sin(j * 2 * M_PI / M), -d_pool * (double)i / N));
//        v.push_back(Vec3d(0, 0, -d_pool));
//        
//        for (int j = 0; j < M; j++)
//            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
//        for (int i = 0; i < N * 2 - 1; i++)
//            for (int j = 0; j < M; j++)
//                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
//                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
//        for (int j = 0; j < M; j++)
//            f.push_back(Vec3i((2 * N) * M + 1, (2 * N - 1) * M + (j + 1) % M + 1, (2 * N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
//        
//        size_t v0 = vs.size();
//        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
//        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
//        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
//        
//        vels.reserve(vels.size() + v.size());
//        for (size_t i = 0; i < v.size(); i++) vels.push_back(Vec3d(0, 0, 0));
//        
//        solids.reserve(solids.size() + v.size());
//        for (size_t i = 0; i < v.size(); i++) solids.push_back(true);
//        for (int i = 0; i < N - 1; i++) for (int j = 0; j < M; j++) solids[v0 + i * M + j + 1] = false;
//        solids[v0] = false;
        
        BPS3D * bps = new BPS3D(sim, std::vector<LosTopos::Vec3d>(), std::vector<LosTopos::Vec3st>(), std::vector<LosTopos::Vec2i>(), std::vector<Vec3d>(), std::vector<char>());
        MeshIO::load(*bps, "init/crownsplash.rec");
        
        size_t v0 = vs.size();
        vs.insert(vs.end(), bps->m_st->pm_positions.begin(), bps->m_st->pm_positions.end());
        fs.reserve(fs.size() + bps->m_st->m_mesh.m_tris.size()); for (size_t i = 0; i < bps->m_st->m_mesh.m_tris.size(); i++) fs.push_back(bps->m_st->m_mesh.m_tris[i] + LosTopos::Vec3st(v0, v0, v0));
        ls.insert(ls.end(), bps->m_st->m_mesh.m_triangle_labels.begin(), bps->m_st->m_mesh.m_triangle_labels.end());

        std::vector<Vec3d> bpsvels = bps->vel();       vels.insert(vels.end(), bpsvels.begin(), bpsvels.end());
        std::vector<char> bpssl = bps->solid_labels(); solids.insert(solids.end(), bpssl.begin(), bpssl.end());
        
        for (size_t i = 0; i < vels.size(); i++)    // clear the velocity: the loaded velocity is not desired.
            vels[i].setZero();
        
        delete bps;
        
        // droplet
        v.clear();
        f.clear();
        l.clear();
        
//        createIcoSphere(Vec3d(0, 0, z_droplet), r_droplet, 2, v, f, l);
        createUVSphere(Vec3d(0, 0, z_droplet), r_droplet, 16, 32, v, f, l);
        
        v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
        
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(Vec3d(Options::doubleValue("init-vel-x"), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z")));
        
        solids.reserve(solids.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) solids.push_back(false);

        return new BPS3D(sim, vs, fs, ls, vels, solids);
        
    } else if (scene == "domecollision")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double sigma_sa = Options::doubleValue("sigma-sa");
        double sigma_sl = Options::doubleValue("sigma-sl");
        double contact_angle = acos(std::min(1.0, std::max(-1.0, sigma_sa - sigma_sl)));
        std::cout << "equilibrium contact angle = " << contact_angle * 180 / M_PI << std::endl;
        
        double alpha_end = contact_angle;
        double r = Options::doubleValue("radius");
        double d = Options::doubleValue("distance") / 2;
        
        v.push_back(Vec3d(-d, 0, 0) + Vec3d(0, 0, r - cos(alpha_end) * r));
        for (int i = 1; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(Vec3d(-d, 0, 0) + Vec3d(sin(i * alpha_end / N) * cos((j - 0.5 * i) * 2 * M_PI / M), sin(i * alpha_end / N) * sin((j - 0.5 * i) * 2 * M_PI / M), cos(i * alpha_end / N) - cos(alpha_end)) * r);
        v.push_back(Vec3d(-d, 0, 0) + Vec3d(0, 0, 0));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(N * M + 1, (N - 1) * M + (j + 1) % M + 1, (N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        size_t v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
        
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(Vec3d( std::abs(Options::doubleValue("init-vel-x")), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z")));
        
        solids.reserve(solids.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) solids.push_back(false);
        for (size_t j = 0; j < M; j++)
            solids[v0 + (N - 1) * M + j + 1] = true;
        solids[v0 + N * M + 1] = true;
        
        v.clear();
        f.clear();
        l.clear();

        v.push_back(Vec3d(d, 0, 0) + Vec3d(0, 0, r - cos(alpha_end) * r));
        for (int i = 1; i <= N; i++)
            for (int j = 0; j < M; j++)
                v.push_back(Vec3d(d, 0, 0) + Vec3d(sin(i * alpha_end / N) * cos((j - 0.5 * i) * 2 * M_PI / M), sin(i * alpha_end / N) * sin((j - 0.5 * i) * 2 * M_PI / M), cos(i * alpha_end / N) - cos(alpha_end)) * r);
        v.push_back(Vec3d(d, 0, 0) + Vec3d(0, 0, 0));
        
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(Vec2i(1, 0));
        for (int i = 0; i < N - 1; i++)
            for (int j = 0; j < M; j++)
                f.push_back(Vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)),     l.push_back(Vec2i(1, 0)),
                f.push_back(Vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(Vec2i(1, 0));
        for (int j = 0; j < M; j++)
            f.push_back(Vec3i(N * M + 1, (N - 1) * M + (j + 1) % M + 1, (N - 1) * M + j + 1)), l.push_back(Vec2i(1, 0));
        
        v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
        
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(Vec3d(-std::abs(Options::doubleValue("init-vel-x")), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z")));
        
        solids.reserve(solids.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) solids.push_back(false);
        for (size_t j = 0; j < M; j++)
            solids[v0 + (N - 1) * M + j + 1] = true;
        solids[v0 + N * M + 1] = true;
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
  
    } else if (scene == "puddlesplash")
    {
        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        double r_droplet = Options::doubleValue("radius");
        double z_droplet = Options::doubleValue("radius") + 0.5;
        
        BPS3D * bps = new BPS3D(sim, std::vector<LosTopos::Vec3d>(), std::vector<LosTopos::Vec3st>(), std::vector<LosTopos::Vec2i>(), std::vector<Vec3d>(), std::vector<char>());
        MeshIO::load(*bps, "init/puddle.rec");
        
        size_t v0 = vs.size();
        vs.insert(vs.end(), bps->m_st->pm_positions.begin(), bps->m_st->pm_positions.end());
        fs.reserve(fs.size() + bps->m_st->m_mesh.m_tris.size()); for (size_t i = 0; i < bps->m_st->m_mesh.m_tris.size(); i++) fs.push_back(bps->m_st->m_mesh.m_tris[i] + LosTopos::Vec3st(v0, v0, v0));
        ls.insert(ls.end(), bps->m_st->m_mesh.m_triangle_labels.begin(), bps->m_st->m_mesh.m_triangle_labels.end());
        
        std::vector<Vec3d> bpsvels = bps->vel();       vels.insert(vels.end(), bpsvels.begin(), bpsvels.end());
        std::vector<char> bpssl = bps->solid_labels(); solids.insert(solids.end(), bpssl.begin(), bpssl.end());
        
        for (size_t i = 0; i < vels.size(); i++)    // clear the velocity: the loaded velocity is not desired.
            vels[i].setZero();
        
        delete bps;
        
        // droplet
        v.clear();
        f.clear();
        l.clear();
        
//        createIcoSphere(Vec3d(0, 0, z_droplet), r_droplet, 2, v, f, l);
        createUVSphere(Vec3d(0, 0, z_droplet), r_droplet, 16, 32, v, f, l);
        
        v0 = vs.size();
        vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
        fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
        ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
        
        vels.reserve(vels.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) vels.push_back(Vec3d(Options::doubleValue("init-vel-x"), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z")));
        
        solids.reserve(solids.size() + v.size());
        for (size_t i = 0; i < v.size(); i++) solids.push_back(false);
        
        return new BPS3D(sim, vs, fs, ls, vels, solids);
    
    } else
    {
        std::cout << "Unrecognized scene name: " << scene << std::endl;
    }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Scene-specific time stepping
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Scenes::pre_step(Sim * sim, const std::string & scene, double dt, BPS3D * bps)
{
    if (scene == "transportinject")
    {
        static double s_transportinject_last_inject = -50;
        if (sim->time() - s_transportinject_last_inject > 50)
        {
            s_transportinject_last_inject = sim->time();
            
            // add a ball to be injected
            std::vector<Vec3d> v;
            std::vector<Vec3i> f;
            std::vector<Vec2i> l;
            
            double r = Options::doubleValue("radius");
            double h = 3.0;
            double d = Options::doubleValue("distance");
            
            int N = Options::intValue("mesh-size-n");
            
            createIcoSphere(Vec3d(0, -d, h), r, N, v, f, l);
//            createUVSphere(Vec3d(0, -d, h), r, N, M, v, f, l);
            
            size_t v0 = sim->bps()->nv();
            for (size_t i = 0; i < v.size(); i++)
                sim->bps()->surfTrack()->add_vertex(vc(v[i]), LosTopos::Vec3d(1, 1, 1));
            for (size_t i = 0; i < f.size(); i++)
                sim->bps()->surfTrack()->add_triangle(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]), LosTopos::Vec2i(l[i][0], l[i][1]));
            
            for (size_t i = 0; i < v.size(); i++)
                sim->bps()->vel(v0 + i) = Vec3d(Options::doubleValue("init-vel-x"), Options::doubleValue("init-vel-y"), Options::doubleValue("init-vel-z"));
        }
        
    }
}

void Scenes::post_step(Sim * sim, const std::string & scene, double dt, BPS3D * bps)
{

}


