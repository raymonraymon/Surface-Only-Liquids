//
//  BPS3D.h
//  Droplet3D
//
//  Created by Fang Da on 10/3/15.
//
//

#ifndef __Droplet3D__BPS3D__
#define __Droplet3D__BPS3D__

#include <stdio.h>
#include <iostream>
#include "eigenheaders.h"
#include "surftrack.h"
#include "Timer.h"

class Sim;
class Scenes;

class BPS3D : public LosTopos::SurfTrack::SolidVerticesCallback, public LosTopos::T1Transition::VelocityFieldCallback, public LosTopos::SurfTrack::MeshEventCallback
{
    friend class Sim;
    friend class Scenes;
    
public:
    BPS3D(Sim * sim, const std::vector<LosTopos::Vec3d> & vs, const std::vector<LosTopos::Vec3st> & fs, const std::vector<LosTopos::Vec2i> & ls, const std::vector<Vec3d> & vels, const std::vector<char> & solids);
    ~BPS3D();
    
    class SimOptions
    {
    public:
        bool implicit;
        double smoothing_coef;
        double sigma;
        double sigma_sl;    // solid-liquid interface surface tension coefficient is sigma * sigma_sl
        double sigma_sa;    //    solid-air interface surface tension coefficient is sigma * sigma_sa
        Vec3d gravity;
        double rho;         // liquid density
        
        double mean_edge_length;
        
        SimOptions() : implicit(false), smoothing_coef(0), sigma(1), sigma_sl(1), sigma_sa(1), gravity(Vec3d(0, 0, 0)), rho(1), mean_edge_length(0)
        { }
    };
    
          SimOptions & simOptions()       { return m_sim_options; }
    const SimOptions & simOptions() const { return m_sim_options; }
    
    void setVerbose(bool verbose) { m_verbose = verbose; m_st->m_verbose = verbose; }
    bool verbose() const { return m_verbose; }
    
public:
    const LosTopos::SurfTrack * surfTrack() const { return m_st; }
          LosTopos::SurfTrack * surfTrack()       { return m_st; }
    const LosTopos::NonDestructiveTriMesh & mesh() const { return m_st->m_mesh; }
          LosTopos::NonDestructiveTriMesh & mesh()       { return m_st->m_mesh; }
    const Sim * sim() const { return m_sim; }
          Sim * sim()       { return m_sim; }

    size_t nv() const { return mesh().nv(); }
    size_t nf() const { return mesh().nt(); }
    
    Vec3d pos(size_t vi) const { return vc(m_st->pm_positions[vi]); }
    
    const Vec3d & vel(size_t vi) const { return (*m_v)[vi]; }
          Vec3d & vel(size_t vi)       { return (*m_v)[vi]; }
    std::vector<Vec3d> vel() const { std::vector<Vec3d> vels(nv()); for (size_t i = 0; i < nv(); i++) vels[i] = vel(i); return vels; }
    VecXd             velv() const { VecXd vels = VecXd::Zero(nv() * 3); for (size_t i = 0; i < nv(); i++) vels.segment<3>(i * 3) = vel(i); return vels; }
    
    std::vector<char> solid_labels() const { std::vector<char> solids(nv()); for (size_t i = 0; i < nv(); i++) solids[i] = vertex_is_solid(i); return solids; }
    
    LosTopos::Vec3st getShuffledTriangle(const LosTopos::Vec3st & t, size_t vertex_to_be_front) const;
    Mat3d getVertexPositions(const LosTopos::Vec3st & t) const;
    Mat3d getVertexVelocities(const LosTopos::Vec3st & t) const;
    
    LosTopos::Vec2st getShuffledEdge(const LosTopos::Vec2st & e, size_t vertex_to_be_front) const;
    Mat32d getVertexPositions(const LosTopos::Vec2st & e) const;
    Mat32d getVertexVelocities(const LosTopos::Vec2st & e) const;
    
    int   intermediate_v_selector() const { return m_intermediate_v_selector; }
    int & intermediate_v_selector()       { return m_intermediate_v_selector; }
    void  intermediate_v_select(int inc)  { m_intermediate_v_selector = (m_intermediate_v_selector + inc + m_intermediate_v.size()) % m_intermediate_v.size(); }
    const std::pair<VecXd, std::string> & intermediate_v(int selector) const { return m_intermediate_v[selector]; }
    const std::pair<VecXd, std::string> & intermediate_v()             const { return m_intermediate_v[m_intermediate_v_selector]; }
    
    const std::vector<std::pair<VecXd, std::string> > & intermediate_vs() const { return m_intermediate_v; }
    
    const VecXd & preStepGeometry() const { return m_pre_stepping_geometry; }
    
public:
    double step(double dt, int substep = 0);
    
    double simTime() const;
    
protected:
    struct Partitioning
    {
        std::vector<std::vector<size_t> > p2v;  // a list of vertices for each partition
        std::vector<std::vector<size_t> > p2f;  // a list of faces for each partition
        std::vector<int> v2p;                   // a partition index for each vertex
        std::vector<int> f2p;                   // a partition index for each face
        
        std::vector<size_t> flattened_partition_vertices;   // flatten out all partitions of vertices into one linear array
        std::vector<size_t> indices_in_partitions;          // for each vertex in flattened_partition_vertices, record its index inside its respective partition
    };
    
    void partition_mesh(Partitioning & partitioning);
    
    void step_explicit(double dt, const Partitioning & partitioning);
    void step_implicit(double dt, const Partitioning & partitioning);
    
    void step_HelmholtzDecomposition(double dt, VecXd & v, const Partitioning & partitioning);
    void step_PressureSolve         (double dt, VecXd & v, const Partitioning & partitioning);
    void step_SmoothVelocity        (double dt, VecXd & v, const Partitioning & partitioning, double coef);
    
    void computeCollocCoefficients(std::vector<MatXd> & A_solid, std::vector<MatXd> & A_air, std::vector<MatXd> & B_solid, std::vector<MatXd> & B_air, const Partitioning & partitioning);
    
    void step_PressureSolveIterative(double dt, VecXd & v, const Partitioning & partitioning);
    void step_HD_FMM                (double dt, VecXd & v, const Partitioning & partitioning);
    
public:
    // convenient mesh topology query
    size_t edge_other_vertex(size_t e, size_t v) const { LosTopos::Vec2st edge = mesh().m_edges[e]; return edge[0] == v ? edge[1] : edge[0]; }
    
    // convenient mesh geometry query
    double face_area(size_t f)    const { return m_st->get_triangle_area(f); }
    Vec3d edge_tangent(size_t e)  const { return (pos(mesh().m_edges[e][1]) - pos(mesh().m_edges[e][0])).normalized(); }
    
    // compute edge/vertex neighborhood areas (for manifold regions)
    double edge_area(size_t e)    const { double a = 0; for (size_t i = 0; i < mesh().m_edge_to_triangle_map[e].size(); i++) a += face_area(mesh().m_edge_to_triangle_map[e][i]); return a / 3; }
    double vert_area(size_t v)    const { double a = 0; for (size_t i = 0; i < mesh().m_vertex_to_triangle_map[v].size(); i++) a += face_area(mesh().m_vertex_to_triangle_map[v][i]); return a / 3; }
    
    double edge_length(size_t e)  const { return m_st->get_edge_length(e); }
    
    // compute the vertex angles
    double vert_interior_solid_angle(size_t v) const;
    double vert_interior_solid_angle_tet_decomposition(size_t v) const;
    double vert_interior_solid_angle_dihedral(size_t v) const;
    double face_interior_angle(size_t f, size_t v) const;
    
    // find the outward normal by assuming the liquid interior has a smaller region label than the exterior
    // note that the normal is not necessarily the same as m_st->get_triangle_normal (depending on the region labels)
    Vec3d  face_outward_normal(size_t f)  const { LosTopos::Vec3st t = mesh().m_tris[f]; Vec3d n = (pos(t[1]) - pos(t[0])).cross(pos(t[2]) - pos(t[0])).normalized(); LosTopos::Vec2i l = mesh().get_triangle_label(f); if (l[0] < l[1]) return -n; else return n; }
    Vec3d  vert_outward_normal(size_t v)  const { Vec3d n = Vec3d::Zero(); for (size_t i = 0; i < mesh().m_vertex_to_triangle_map[v].size(); i++) n += face_area(mesh().m_vertex_to_triangle_map[v][i]) * face_outward_normal(mesh().m_vertex_to_triangle_map[v][i]); return n.normalized(); }
    
    // mean curvature (for pressure BC)
    double vert_integral_mean_curvature(size_t v) const;
    
    // special handling for surface tension force on the triple junction: pretend the faces incident on the triple junction has a small (global constant) width away from the triple junction, because the geometry has a sharp crease here
    double triple_junction_virtual_width() const { return simOptions().mean_edge_length * 0.25; }
    
    // solid contact query
    bool vertex_is_solid(size_t v) const { return m_st->vertex_is_any_solid(v); }
    bool face_is_solid(size_t f)   const { LosTopos::Vec3st t = mesh().m_tris[f]; return (vertex_is_solid(t[0]) && vertex_is_solid(t[1]) && vertex_is_solid(t[2])); }
    
    bool vertex_is_on_triple_junction(size_t v) const { if (!vertex_is_solid(v)) return false; for (size_t i = 0; i < mesh().m_vertex_to_triangle_map[v].size(); i++) if (!face_is_solid(mesh().m_vertex_to_triangle_map[v][i])) return true; return false; }
    bool edge_is_on_triple_junction(size_t e) const { bool solid = false; bool air = false; for (size_t i = 0; i < mesh().m_edge_to_triangle_map[e].size(); i++) if (face_is_solid(mesh().m_edge_to_triangle_map[e][i])) solid = true; else air = true; return (solid && air); } // note that just testing if both endpoints of the edge are on the tirple junction won't work.
    
    // volume query
    double liquid_volume() const;
    double liquid_volume_change_rate() const;
    double liquid_volume_change_rate(const VecXd & v) const;
    
    double partition_volume(const std::vector<size_t> & faces) const;
    
protected:
    // SurfTrack::SolidVerticesCallback method
    bool            generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    bool            generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    bool            generate_snapped_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    LosTopos::Vec3c generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    LosTopos::Vec3c generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    LosTopos::Vec3c generate_snapped_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    bool            generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower);
    bool            generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b);
    bool            solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e);
    
    // T1Transition::VelocityFieldCallback methods
    LosTopos::Vec3d sampleVelocity(LosTopos::Vec3d & pos);
    bool sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output);
    
    // SurfTrack::MeshEventCallback
    void pre_collapse(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_collapse(const LosTopos::SurfTrack & st, size_t e, size_t merged_vertex, void * data);
    
    void pre_split(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_split(const LosTopos::SurfTrack & st, size_t e, size_t new_vertex, void * data);
    
    void pre_flip(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_flip(const LosTopos::SurfTrack & st, size_t e, void * data);
    
    void pre_t1(const LosTopos::SurfTrack & st, size_t v, void ** data);
    void post_t1(const LosTopos::SurfTrack & st, size_t v, size_t a, size_t b, void * data);
    
    void pre_facesplit(const LosTopos::SurfTrack & st, size_t f, void ** data);
    void post_facesplit(const LosTopos::SurfTrack & st, size_t f, size_t new_vertex, void * data);
    
    void pre_snap(const LosTopos::SurfTrack & st, size_t v0, size_t v1, void ** data);
    void post_snap(const LosTopos::SurfTrack & st, size_t v_kept, size_t v_deleted, void * data);
    
    void pre_smoothing(const LosTopos::SurfTrack & st, void ** data);
    void post_smoothing(const LosTopos::SurfTrack & st, void * data);
    
    std::ostream & log() { static std::stringstream ss; return (m_verbose ? std::cout : ss); }
    std::ostream & callback_log() { static std::stringstream ss; return (m_verbose ? std::cout : ss); }
    
protected:
    Sim * m_sim;
    
    LosTopos::SurfTrack * m_st;
    
    SimOptions m_sim_options;
    
    // dynamics
    LosTopos::NonDestructiveTriMesh::VertexData<Vec3d> * m_v;

    int m_intermediate_v_selector;
    std::vector<std::pair<VecXd, std::string> > m_intermediate_v;

    VecXd m_pre_stepping_geometry;  // this stores the vertex positions after the remeshing at the beginning of the step, but before the actually time stepping
        
    VecXd m_dbg_v1;
    
    bool m_verbose;
};


#endif /* defined(__Droplet3D__BPS3D__) */
