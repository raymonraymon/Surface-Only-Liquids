//
//  MeshIO.cpp
//
//  Christopher Batty, Fang Da 2014
//
//

#include <fstream>
#include "MeshIO.h"
#include <map>

#define MESHIO_REC_VERSION  (0)

bool MeshIO::save(BPS3D & bps, const std::string & filename, bool binary)
{
    LosTopos::SurfTrack & st = *(bps.surfTrack());
    
    if (binary)
    {
        std::ofstream os(filename.c_str(), std::ios::binary);
        size_t n;
        
        // version number
        n = MESHIO_REC_VERSION;
        os.write((char *)&n , sizeof (size_t));
        
        // write vertices
        n = st.m_mesh.nv();
        os.write((char *)&n, sizeof (size_t));
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x = st.get_position(i);
            os.write((char *)&(x[0]), sizeof (x[0]));
            os.write((char *)&(x[1]), sizeof (x[1]));
            os.write((char *)&(x[2]), sizeof (x[2]));
            
            Vec3d v = bps.vel(i);
            os.write((char *)&(v[0]), sizeof (v[0]));
            os.write((char *)&(v[1]), sizeof (v[1]));
            os.write((char *)&(v[2]), sizeof (v[2]));
        }
        
        // write faces
        n = st.m_mesh.nt();
        os.write((char *)&n, sizeof (size_t));
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t = st.m_mesh.get_triangle(i);
            if (t[0] == t[1])
                continue;
            LosTopos::Vec2i l = st.m_mesh.get_triangle_label(i);
            
            os.write((char *)&(t[0]), sizeof (t[0]));
            os.write((char *)&(t[1]), sizeof (t[1]));
            os.write((char *)&(t[2]), sizeof (t[2]));
            
            os.write((char *)&(l[0]), sizeof (l[0]));
            os.write((char *)&(l[1]), sizeof (l[1]));
        }
        
        // write the solid vertex indices
        n = 0;
        for (size_t i = 0; i < bps.nv(); i++)
            if (bps.vertex_is_solid(i))
                n++;
        os.write((char *)&n, sizeof (size_t));
        
        for (size_t i = 0; i < bps.nv(); i++)
            if (bps.vertex_is_solid(i))
                os.write((char *)&i, sizeof (size_t));
        
        os.close();
        
        return os.good();
    } else
    {
        std::ofstream os(filename.c_str());

        // write the version number
        os << MESHIO_REC_VERSION << std::endl;

        // write the vertices
        os << st.m_mesh.nv() << std::endl;
        for (size_t i = 0; i < st.m_mesh.nv(); i++)
        {
            os << st.get_position(i) << " " << bps.vel(i)[0] << " " << bps.vel(i)[1] << " " << bps.vel(i)[2] << std::endl;
        }
        
        // write the faces
        os << st.m_mesh.nt() << std::endl;
        for (size_t i = 0; i < st.m_mesh.nt(); i++)
        {
            LosTopos::Vec3st t = st.m_mesh.get_triangle(i);
            if (t[0] == t[1])
                continue;
            LosTopos::Vec2i l = st.m_mesh.get_triangle_label(i);
            
            os << t << " " << l << std::endl;
        }
        
        // write the solid vertex indices
        // write the solid vertex indices
        size_t n = 0;
        for (size_t i = 0; i < bps.nv(); i++)
            if (bps.vertex_is_solid(i))
                n++;
        os << n << std::endl;
        
        for (size_t i = 0; i < bps.nv(); i++)
            if (bps.vertex_is_solid(i))
                os << i << std::endl;
        
        os.close();
        
        return os.good();
    }
}

bool MeshIO::load(BPS3D & bps, const std::string & filename, bool binary)
{
    LosTopos::SurfTrack & st = *(bps.surfTrack());

    std::vector<LosTopos::Vec3d> xs;
    std::vector<LosTopos::Vec3d> vels;
    std::vector<LosTopos::Vec3st> fs;
    std::vector<LosTopos::Vec2i> ls;
    std::vector<size_t> solid_vertices;
    
    bool good = loadIntoRaw(xs, fs, ls, vels, solid_vertices, filename, binary);
    if (!good)
        return false;

    // copy the loaded data to the surface tracker
    for (size_t i = 0; i < st.m_mesh.nt(); i++)
    {
        if (st.m_mesh.get_triangle(i)[0] == st.m_mesh.get_triangle(i)[1])
            continue;
        st.remove_triangle(i);
    }
    
    for (size_t i = 0; i < st.m_mesh.nv(); i++)
        st.remove_vertex(i);
    
    size_t nv = xs.size();
    st.m_mesh.set_num_vertices(nv);

    st.m_masses.resize(nv);
    for (size_t i = 0; i < nv; i++)
        st.m_masses[i] = LosTopos::Vec3d(1, 1, 1);
    for (size_t i = 0; i < solid_vertices.size(); i++)
        st.m_masses[solid_vertices[i]] = LosTopos::Vec3d(1, 1, std::numeric_limits<double>::infinity()); //&&&& for now assume that the solid surface is always perpendicular to the z axis. For generality, LosTopos's impulse code will need a major overhaul anyway.
    
    st.pm_positions = xs;
    st.pm_newpositions = xs;
    st.set_all_remesh_velocities(vels);

    st.m_mesh.replace_all_triangles(fs, ls);
    
    assert(nv == st.m_mesh.m_vertex_to_triangle_map.size());
    st.pm_positions.resize(nv);
    st.pm_newpositions.resize(nv);
    st.pm_velocities.resize(nv);
    st.m_velocities.resize(nv);
    
    st.set_all_positions(xs);
    st.set_all_newpositions(xs);
    
    for (size_t i = 0; i < bps.mesh().nv(); i++)
        bps.vel(i) = vc(vels[i]);
    
    return good;
}


bool MeshIO::loadIntoRaw(std::vector<LosTopos::Vec3d> & xs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<LosTopos::Vec3d> & vels, std::vector<size_t> & solid_vertices, const std::string & filename, bool binary)
{
    std::ifstream test(filename.c_str());
    if (!test.is_open())
    {
        std::cout << "[MeshIO::load] Error: file " << filename << " not found." << std::endl;
        return false;
    }
 
    if (binary)
    {
        std::ifstream is(filename.c_str(), std::ios::binary);
        
        // load and confirm the version
        size_t version;
        is.read((char *)&version, sizeof (size_t));
        if (version != MESHIO_REC_VERSION)
        {
            std::cout << "[MeshIO::load] Error: the REC file " << filename << " has version " << version << ", which is different than the implemenation version " << MESHIO_REC_VERSION << std::endl;
            is.close();
            return false;
        }
        
        size_t n;
        
        // load vertices
        is.read((char *)&n, sizeof (size_t));
        xs.resize(n);
        vels.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x;
            is.read((char *)&(x[0]), sizeof (x[0]));
            is.read((char *)&(x[1]), sizeof (x[1]));
            is.read((char *)&(x[2]), sizeof (x[2]));
            xs[i] = x;
            
            LosTopos::Vec3d v;
            is.read((char *)&(v[0]), sizeof (v[0]));
            is.read((char *)&(v[1]), sizeof (v[1]));
            is.read((char *)&(v[2]), sizeof (v[2]));
            vels[i] = v;
        }
        
        // load faces
        is.read((char *)&n, sizeof (size_t));
        fs.resize(n);
        ls.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t;
            is.read((char *)&(t[0]), sizeof (t[0]));
            is.read((char *)&(t[1]), sizeof (t[1]));
            is.read((char *)&(t[2]), sizeof (t[2]));
            fs[i] = t;
            
            LosTopos::Vec2i l;
            is.read((char *)&(l[0]), sizeof (l[0]));
            is.read((char *)&(l[1]), sizeof (l[1]));
            ls[i] = l;
        }
        
        // load the solid vertex indices
        is.read((char *)&n, sizeof (size_t));
        solid_vertices.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            size_t cv;
            is.read((char *)&cv, sizeof (cv));
            solid_vertices[i] = cv;
        }
        
        is.close();
        
        return is.good();
        
    } else
    {
        std::ifstream is(filename.c_str());
        
        // load and confirm the version
        size_t version;
        is >> version;
        if (version != MESHIO_REC_VERSION)
        {
            std::cout << "[MeshIO::load] Error: the REC file " << filename << " has version " << version << ", which is different than the implemenation version " << MESHIO_REC_VERSION << std::endl;
            is.close();
            return false;
        }
        
        size_t n;
        
        // load vertices
        is >> n;
        xs.resize(n);
        vels.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x;
            is >> x[0] >> x[1] >> x[2];
            xs[i] = x;
            
            LosTopos::Vec3d v;
            is >> v[0] >> v[1] >> v[2];
            vels[i] = v;
        }
        
        // load faces
        is >> n;
        fs.resize(n);
        ls.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t;
            is >> t[0] >> t[1] >> t[2];
            fs[i] = t;
            
            LosTopos::Vec2i l;
            is >> l[0] >> l[1];
            ls[i] = l;
        }
        
        // load solid vertex indices
        is >> n;
        solid_vertices.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            size_t cv;
            is >> cv;
            solid_vertices[i] = cv;
        }
        
        is.close();
        
        return is.good();
    }
}

namespace
{
    struct VertexPatchPair
    {
        size_t v;
        Vec2i rp;   // must satisfy rp[0] < rp[1]
        
        bool operator () (const VertexPatchPair & vpp1, const VertexPatchPair & vpp2) const { return (vpp1.v < vpp2.v || (vpp1.v == vpp2.v && (vpp1.rp[0] < vpp2.rp[0] || (vpp1.rp[0] == vpp2.rp[0] && vpp1.rp[1] < vpp2.rp[1])))); }
    };

}

bool MeshIO::saveOBJ(BPS3D & bps, const std::string & filename)
{
    std::cout << "Saving mesh to " << filename << "... ";
    
    // define all the normals (for manifold vertices, average the incident face normals weighted by area; for nonmanifold vertex, one normal instance is created for each manifold patch.)
    std::map<VertexPatchPair, int, VertexPatchPair> vertex_normal_index;
    std::vector<Vec3d> vertex_normals;
    for (size_t i = 0; i < bps.mesh().nv(); i++)
    {
        std::set<Vec2i, Vec2iComp> incident_regions;
        for (size_t j = 0; j < bps.mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            LosTopos::Vec2i l = bps.mesh().get_triangle_label(bps.mesh().m_vertex_to_triangle_map[i][j]);
            Vec2i rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
            incident_regions.insert(rp);
        }
        
        for (std::set<Vec2i, Vec2iComp>::iterator j = incident_regions.begin(); j != incident_regions.end(); j++)
        {
            VertexPatchPair vpp;
            vpp.v = i;
            vpp.rp = *j;
            vertex_normal_index[vpp] = (int)vertex_normals.size();
            Vec3d n(0, 0, 0);
            for (size_t k = 0; k < bps.mesh().m_vertex_to_triangle_map[i].size(); k++)
            {
                LosTopos::Vec2i l = bps.mesh().get_triangle_label(bps.mesh().m_vertex_to_triangle_map[i][k]);
                if ((l[0] == (*j)[0] && l[1] == (*j)[1]) || (l[1] == (*j)[0] && l[0] == (*j)[1]))
                {
                    LosTopos::Vec3st t = bps.mesh().get_triangle(bps.mesh().m_vertex_to_triangle_map[i][k]);
                    Vec3d fn = (bps.pos(t[1]) - bps.pos(t[0])).cross(bps.pos(t[2]) - bps.pos(t[0])) * (l[0] == (*j)[0] ? 1 : -1);
                    n += fn;
                }
            }
            n.normalize();
            vertex_normals.push_back(n);
        }
    }
    
    std::ofstream os(filename.c_str());
    for (size_t i = 0; i < bps.mesh().nv(); i++)
        os << "v " << bps.pos(i).x() << " " << bps.pos(i).y() << " " << bps.pos(i).z() << std::endl;
    
    for (size_t i = 0; i < vertex_normals.size(); i++)
        os << "vn " << vertex_normals[i].x() << " " << vertex_normals[i].y() << " " << vertex_normals[i].z() << std::endl;
    
    for (size_t i = 0; i < bps.mesh().nt(); i++)
    {
        LosTopos::Vec3st t = bps.mesh().get_triangle(i);
        LosTopos::Vec2i l = bps.mesh().get_triangle_label(i);
        VertexPatchPair vpp;
        vpp.rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
        vpp.v = t[0];   int vn0 = vertex_normal_index[vpp];
        vpp.v = t[1];   int vn1 = vertex_normal_index[vpp];
        vpp.v = t[2];   int vn2 = vertex_normal_index[vpp];
        os << "f " << t[0] + 1 << "//" << vn0 + 1 << " " << t[1] + 1 << "//" << vn1 + 1 << " " << t[2] + 1 << "//" << vn2 + 1 << std::endl;
    }
    
    os.close();

    std::cout << "done." << std::endl;
    
    return true;
}

