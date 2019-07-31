//
//  BPS3D.cpp
//  Droplet3D
//
//  Created by Fang Da on 10/3/15.
//
//

#include "BPS3D.h"
#include "SimOptions.h"
#include "LosTopos/LosTopos3D/subdivisionscheme.h"
#include "Sim.h"
#include "FMMTLWrapper.h"

namespace
{
    double angleAroundAxis(const Vec3d & v0, const Vec3d & v1, const Vec3d & a)   // angle from v0 to v1 around axis a
    {
        double asq = a.squaredNorm();
        assert(asq != 0);
        
        Vec3d u = v0 - v0.dot(a) * a / asq;
        assert(u.squaredNorm() != 0);
        u.normalize();
        
        Vec3d v = a.cross(u).normalized();
        
        return atan2(v1.dot(v), v1.dot(u));
    }
    
}

BPS3D::BPS3D(Sim * sim, const std::vector<LosTopos::Vec3d> & vs, const std::vector<LosTopos::Vec3st> & fs, const std::vector<LosTopos::Vec2i> & ls, const std::vector<Vec3d> & vels, const std::vector<char> & solids) :
    m_sim(sim),
    m_verbose(false)
{
    // load sim options
    m_sim_options.implicit = Options::boolValue("implicit-integration");
    m_sim_options.smoothing_coef = Options::doubleValue("smoothing-coef");
    m_sim_options.sigma = Options::doubleValue("sigma");
    m_sim_options.sigma_sl = Options::doubleValue("sigma-sl");
    m_sim_options.sigma_sa = Options::doubleValue("sigma-sa");
    Vec3d gravity_dir = Vec3d(Options::doubleValue("gravity-dir-x"), Options::doubleValue("gravity-dir-y"), Options::doubleValue("gravity-dir-z"));
    m_sim_options.gravity = Options::doubleValue("gravity") * gravity_dir.normalized();
    m_sim_options.rho = Options::doubleValue("rho");
    
    // construct the surface tracker
    double target_edge_len = Options::doubleValue("remeshing-resolution");
    if (target_edge_len == 0)
    {
        double mean_edge_len = 0;
        for (size_t i = 0; i < fs.size(); i++)
        {
            mean_edge_len += mag(vs[fs[i][0]] - vs[fs[i][1]]);
            mean_edge_len += mag(vs[fs[i][1]] - vs[fs[i][2]]);
            mean_edge_len += mag(vs[fs[i][2]] - vs[fs[i][0]]);
        }
        mean_edge_len /= (fs.size() * 3);
        target_edge_len = mean_edge_len;
    }
    m_sim_options.mean_edge_length = target_edge_len;
    
    std::cout << "target edge length = " << target_edge_len << std::endl;
    
    LosTopos::SurfTrackInitializationParameters params;
    params.m_proximity_epsilon = Options::doubleValue("lostopos-collision-epsilon-fraction") * target_edge_len;
    params.m_merge_proximity_epsilon = Options::doubleValue("lostopos-merge-proximity-epsilon-fraction") * target_edge_len;
    params.m_merge_proximity_epsilon_for_liquid_sheet_puncture = params.m_merge_proximity_epsilon * Options::doubleValue("lspmt-ratio");
    params.m_allow_vertex_movement_during_collapse = true;
    params.m_perform_smoothing = Options::boolValue("lostopos-perform-smoothing");
//    params.m_min_to_target_ratio = 0.5;   // default
//    params.m_max_to_target_ratio = 1.5;   // default
    params.m_max_adjacent_target_edge_length_ratio = Options::doubleValue("lostopos-iar-adjacent-ratio");   // default: 1.5
    params.m_target_edge_length_coef_curvature = Options::doubleValue("lostopos-iar-coef-curvature");   // default: 0.05
    params.m_target_edge_length_coef_velocity = Options::doubleValue("lostopos-iar-coef-velocity") / Options::doubleValue("time-step"); // default: 0.01/dt
    params.m_min_target_edge_length = target_edge_len;
    params.m_max_target_edge_length = target_edge_len * Options::doubleValue("remeshing-dr");
    params.m_refine_solid = !Options::boolValue("nrs");
    params.m_refine_triple_junction = !Options::boolValue("nrtp");
    params.m_max_volume_change = Options::doubleValue("lostopos-max-volume-change-fraction") * pow(target_edge_len, 3);
    params.m_min_triangle_angle = Options::doubleValue("lostopos-min-triangle-angle");
    params.m_max_triangle_angle = Options::doubleValue("lostopos-max-triangle-angle");
    params.m_large_triangle_angle_to_split = Options::doubleValue("lostopos-large-triangle-angle-to-split");
    params.m_min_triangle_area = Options::doubleValue("lostopos-min-triangle-area-fraction") * pow(target_edge_len, 2);
    params.m_verbose = false;
    params.m_allow_non_manifold = Options::boolValue("lostopos-allow-non-manifold");
    params.m_allow_topology_changes = Options::boolValue("lostopos-allow-topology-changes");
    params.m_collision_safety = true;
    params.m_remesh_boundaries = true;
    params.m_t1_transition_enabled = Options::boolValue("lostopos-t1-transition-enabled");
    params.m_pull_apart_distance = Options::doubleValue("lostopos-t1-pull-apart-distance-fraction") * target_edge_len;
    
    params.m_velocity_field_callback = NULL;
    
    if (Options::boolValue("lostopos-smooth-subdivision"))
        params.m_subdivision_scheme = new LosTopos::ButterflyScheme();
    else
        params.m_subdivision_scheme = new LosTopos::MidpointScheme();
    
    params.m_use_curvature_when_collapsing = false;
    params.m_use_curvature_when_splitting = false;
    
    std::vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
    for (size_t i = 0; i < vs.size(); i++)  // the solids labels supercede the constrained vertices; the latter has really no effect.
        masses[i] = LosTopos::Vec3d(1, 1, solids[i] ? std::numeric_limits<double>::infinity() : 1); //&&&& for now assume that the solid surface is always perpendicular to the z axis. For generality, LosTopos's impulse code will need a major overhaul anyway.
    m_st = new LosTopos::SurfTrack(vs, fs, ls, masses, params);
    m_st->m_solid_vertices_callback = this;
    m_st->m_mesheventcallback = this;
    
    
    // initialize the dynamics quantities
    m_v = new LosTopos::NonDestructiveTriMesh::VertexData<Vec3d>(&(m_st->m_mesh));
    for (size_t i = 0; i < mesh().nv(); i++)
        (*m_v)[i] = vels[i];

    m_pre_stepping_geometry = VecXd::Zero(nv() * 3);
    for (size_t i = 0; i < nv(); i++)
        m_pre_stepping_geometry.segment<3>(i * 3) = pos(i);
    
    m_intermediate_v_selector = 0;
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(velv(), "Scene initial velocity"));
    
}

BPS3D::~BPS3D()
{
    if (m_st)
        delete m_st;
}

namespace
{
    Mat3d skewSymmetric(const Vec3d & v)
    {
        Mat3d ss = Mat3d::Zero();
        ss(0, 1) = -v(2);
        ss(1, 0) =  v(2);
        ss(0, 2) =  v(1);
        ss(2, 0) = -v(1);
        ss(1, 2) = -v(0);
        ss(2, 1) =  v(0);
        return ss;
    }
}

double BPS3D::step(double dt, int substep)
{
    CSim::TimerMan::timer("Sim.step/BPS.step").start();
    
    m_dbg_v1 = VecXd::Zero(nv());

    double actual_dt = 0;
    if (substep == 0 || substep == 1)
    {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////
    ////    mesh advection
    ////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // prior to integrating, remove the influx velocity from solid vertices, so that their velocities represent the solid velocity, not the liquid velocity
    CSim::TimerMan::timer("Sim.step/BPS.step/solidinflux").start();
    double influx = sim()->currentInflux();
    VecXd nsolid = VecXd::Zero(nv() * 3);
    if (influx != 0)
    {
        for (size_t i = 0; i < nf(); i++)
        {
            if (face_is_solid(i))
            {
                Vec3d n = face_outward_normal(i);
                LosTopos::Vec3st t = mesh().m_tris[i];
                nsolid.segment<3>(t[0] * 3) = Vec3d(nsolid.segment<3>(t[0] * 3) + n).normalized();
                nsolid.segment<3>(t[1] * 3) = Vec3d(nsolid.segment<3>(t[1] * 3) + n).normalized();
                nsolid.segment<3>(t[2] * 3) = Vec3d(nsolid.segment<3>(t[2] * 3) + n).normalized();
            }
        }
        
        for (size_t i = 0; i < nv(); i++)
        {
            if (vertex_is_solid(i))
            {
                Vec3d n = nsolid.segment<3>(i * 3);
                assert(n.squaredNorm() != 0);
                vel(i) = vel(i) - vel(i).dot(n) * n; // now the velocity of solid vertices represents the solid velocity, not the liquid velocity
            }
        }
    }
    CSim::TimerMan::timer("Sim.step/BPS.step/solidinflux").stop();

    CSim::TimerMan::timer("Sim.step/BPS.step/advection").start();
    std::cout << "Advecting mesh..." << std::endl;
    for (size_t i = 0; i < nv(); i++)
        m_st->pm_newpositions[i] = m_st->pm_positions[i] + vc(vel(i)) * dt;

    m_st->rebuild_continuous_broad_phase();
    m_st->integrate(dt, actual_dt); // advection
        if (actual_dt != dt)
            std::cout << "Warning: SurfTrack::integrate() failed to step the full length of the time step!" << std::endl;

    // test for floor contact
    bool solid_integrate_needed = false;
    double floorz = Options::doubleValue("floor-z");
    std::vector<bool> floor_contact(nv(), false);
    for (size_t i = 0; i < nv(); i++)
        if (m_st->pm_newpositions[i][2] <= floorz)
        {
            if (m_st->pm_newpositions[i][2] < floorz)
                solid_integrate_needed = true;  // if a vertex has come strictly below the floor, it'll be moved onto the floor and thus another SurfTrack()::integrate() call will be necessary
            m_st->pm_newpositions[i][2] = floorz;
            m_st->m_masses[i] = LosTopos::Vec3d(1, 1, std::numeric_limits<double>::infinity());
        }
    
    // unlabel any single-vertex-solid-contact as not solid, because the pressure solve assumes any solid contact region will have at least one face.
    for (size_t i = 0; i < nv(); i++)
    {
        bool solid_face = false;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            if (face_is_solid(mesh().m_vertex_to_triangle_map[i][j]))
                solid_face = true;
        if (!solid_face)
            m_st->m_masses[i] = LosTopos::Vec3d(1, 1, 1);   // a solid vertex incident to no solid faces is not a true solid vertex.
    }

    if (solid_integrate_needed)
    {
        m_st->rebuild_continuous_broad_phase();
        m_st->integrate(dt, actual_dt);
    }
    
    for (size_t i = 0; i < nv(); i++)
        m_st->pm_velocities[i] = vc(vel(i));
    
    CSim::TimerMan::timer("Sim.step/BPS.step/advection").stop();

    CSim::TimerMan::timer("Sim.step/BPS.step/solidinflux").start();
    if (influx != 0)
    {
        for (size_t i = 0; i < nv(); i++)
        {
            if (vertex_is_solid(i))
            {
                Vec3d n = nsolid.segment<3>(i * 3);
                assert(n.squaredNorm() != 0);
                vel(i) = vel(i) - vel(i).dot(n) * n - influx * n; // now the velocity of solid vertices represents the liquid velocity, not the solid velocity
            }
        }
    }
    CSim::TimerMan::timer("Sim.step/BPS.step/solidinflux").stop();
    }
    
    
    
    if (substep == 0 || substep == 2)
    {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////
    ////    mesh improvement and defragmentation
    ////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    CSim::TimerMan::timer("Sim.step/BPS.step/remesh").start();
    std::cout << "Remeshing..." << std::endl;
    for (size_t i = 0; i < nv(); i++)
        m_st->pm_velocities[i] = vc(vel(i));
    for(int i = 0; i < Options::intValue("remeshing-iterations"); i++)
    {
        // recompute the target edge lengths (these will be maintained incrementally when performing the remeshing operations, but incremental updates may be too conservative so they are recomputed from scratch here)
        CSim::TimerMan::timer("Sim.step/BPS.step/remesh/target_lengths").start();
        m_st->compute_all_vertex_target_edge_lengths();
        CSim::TimerMan::timer("Sim.step/BPS.step/remesh/target_lengths").stop();
        
        CSim::TimerMan::timer("Sim.step/BPS.step/remesh/topology_changes").start();
        m_st->topology_changes();
        CSim::TimerMan::timer("Sim.step/BPS.step/remesh/topology_changes").stop();
        
        CSim::TimerMan::timer("Sim.step/BPS.step/remesh/target_lengths").start();
        m_st->compute_all_vertex_target_edge_lengths();
        CSim::TimerMan::timer("Sim.step/BPS.step/remesh/target_lengths").stop();

        CSim::TimerMan::timer("Sim.step/BPS.step/remesh/improve_mesh").start();
        m_st->improve_mesh();
        CSim::TimerMan::timer("Sim.step/BPS.step/remesh/improve_mesh").stop();
    }
    
    CSim::TimerMan::timer("Sim.step/BPS.step/remesh/defrag").start();
    std::vector<size_t> dummy;
    m_st->defrag_mesh_from_scratch(dummy);
    CSim::TimerMan::timer("Sim.step/BPS.step/remesh/defrag").start();
    
    // confirm floor contact (there may be cases where a vertex is exactly on the floor but can't be labeled as solid by the floor detection code above integrate() because it'd be a branch from the triple junction curve (i.e. it is incident to no solid face) and the remeshing steps above may have made it okay to label it as solid. if we don't label it as solid in such cases, the pressure solve will have trouble because that the triple junction will have a vertex whose air and solid normals are identical (because the air faces are really already on the floor; they just haven't been labeled as solid yet), which makes the triple junciton tangent NaN and the BEM solve badly conditioned. see the crash in 1449952789)
    double floorz = Options::doubleValue("floor-z");
    for (size_t i = 0; i < nv(); i++)
        if (m_st->pm_positions[i][2] == floorz)
            m_st->m_masses[i] = LosTopos::Vec3d(1, 1, std::numeric_limits<double>::infinity());

    // unlabel any single-vertex-solid-contact as not solid, because the pressure solve assumes any solid contact region will have at least one face.
    for (size_t i = 0; i < nv(); i++)
    {
        bool solid_face = false;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            if (face_is_solid(mesh().m_vertex_to_triangle_map[i][j]))
                solid_face = true;
        if (!solid_face)
            m_st->m_masses[i] = LosTopos::Vec3d(1, 1, 1);   // a solid vertex incident to no solid faces is not a true solid vertex.
    }

    m_pre_stepping_geometry = VecXd::Zero(nv() * 3);
    for (size_t i = 0; i < nv(); i++)
        m_pre_stepping_geometry.segment<3>(i * 3) = pos(i); // store the geometry immediately after the remeshing, for display as comparison to the end of step geometry to debug the remeshing results
    std::cout << "remeshing finished with nv = " << nv() << " nf = " << nf() << std::endl;
    CSim::TimerMan::timer("Sim.step/BPS.step/remesh").stop();
    }
    
    
    
    Partitioning partitioning;
    if (substep == 0 || substep == 3)
    {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////
    ////    partition the mesh into connected components, so that each can be treated separately in HD and BEM (because they have nonlinear complexity)
    ////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////=
    CSim::TimerMan::timer("Sim.step/BPS.step/partition").start();
    partition_mesh(partitioning);
    std::cout << nv() << " vertices, " << nf() << " faces, " << partitioning.p2v.size() << " partitions" << std::endl;
  
    // find the partition with the largest positive volume: it may be used below
    size_t largest_partition = 0;
    double largest_volume = 0;
    for (size_t pi = 0; pi < partitioning.p2v.size(); pi++)
    {
        double volume = partition_volume(partitioning.p2v[pi]);
        if (volume > largest_volume)
            largest_volume = volume,
            largest_partition = pi;
    }
        
    bool deletion_happened = false;
    bool merging_happened = false;
//    Vec3d o(0, 0, 0);
    for (size_t pi = 0; pi < partitioning.p2f.size(); pi++)
    {
        double volume = partition_volume(partitioning.p2f[pi]);
        
        bool outofview = true;
        double viewrange = std::abs(Options::doubleValue("remove-partitions-below-z"));
        for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
            if (pos(partitioning.p2v[pi][i]).norm() < viewrange)
                outofview = false;
        
        if ((Options::boolValue("ibr") && volume < 0) || outofview)
        {
            // IBR
            // remove any closed interior bubbles (isolated closed surfaces with negative liquid volume). They are usually by-products of merging and are created at close to zero volumes, but due to our air boundary condition, our solver doesn't have the ability to handle them in a meaningful way.
            // also remove any droplet that's far away from the center of the scene
            for (size_t i = 0; i < partitioning.p2f[pi].size(); i++)
                surfTrack()->remove_triangle(partitioning.p2f[pi][i]);
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
                surfTrack()->remove_vertex(partitioning.p2v[pi][i]);
            deletion_happened = true;

            // clear this partition
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
                partitioning.v2p[partitioning.p2v[pi][i]] = -1;
            for (size_t i = 0; i < partitioning.p2f[pi].size(); i++)
                partitioning.f2p[partitioning.p2f[pi][i]] = -1;
            partitioning.p2v[pi].clear();
            partitioning.p2f[pi].clear();
            
        } else if (volume < 0)
        {
            // this is the case of an internal bubble, but it hasn't been deleted in the above branch because IBR is turned off.
            // although we don't have a theory for this, we should do it by treating the entire connected liquid body together, i.e. doing HD on the inward facing (internal
            //  bubble) partitions and the outward facing partition together. this is because the outward facing partition's interior velocity field is not zero, and it affects
            //  the internal bubble partitions. at least this will get rid of the sudden change of behavior upon closing of a concavity which forms a bubble.
            // TODO: ideally we should be able to do a point-mesh inside/outside test to figure out which outward facing partition contains this inward facing partition.
            // For now we don't have time to do this. Let's just say the outward facing partition containing this bubble is the one with the largest positive volume.
            merging_happened = true;
            assert(pi != largest_partition);    // a negative volume partition can't be the largest partition.
            
            partitioning.p2v[largest_partition].insert(partitioning.p2v[largest_partition].end(), partitioning.p2v[pi].begin(), partitioning.p2v[pi].end());
            partitioning.p2f[largest_partition].insert(partitioning.p2f[largest_partition].end(), partitioning.p2f[pi].begin(), partitioning.p2f[pi].end());
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
                partitioning.v2p[partitioning.p2v[pi][i]] = largest_partition;
            for (size_t i = 0; i < partitioning.p2f[pi].size(); i++)
                partitioning.f2p[partitioning.p2f[pi][i]] = largest_partition;
            partitioning.p2v[pi].clear();
            partitioning.p2f[pi].clear();
            
        }
    }
    
    if (deletion_happened)
    {
        // if there has been deletion, the mesh needs to be defragged; as a result, the partitioning's indices will be invalidated and need to be mapped/rebuilt
        std::vector<size_t> defrag_vmap;
        defrag_vmap.reserve(nv());
        for (size_t i = 0; i < nv(); i++)
            defrag_vmap.push_back(i);
        m_st->defrag_mesh_from_scratch(defrag_vmap);
        
        // rebuild p2v, p2f, v2p, f2p
        // map the vertex indices
        for (size_t pi = 0; pi < partitioning.p2v.size(); pi++)
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
            {
                partitioning.p2v[pi][i] = defrag_vmap[partitioning.p2v[pi][i]];
                assert(partitioning.p2v[pi][i] >= 0);   // if a vertex is deleted, its partition should have already been emptied, so this shouldn't happen
            }

        // rebuild the inverse maps
        std::vector<int> f2p(nf(), -1);
        std::vector<int> v2p(nf(), -1);
        for (size_t pi = 0; pi < partitioning.p2v.size(); pi++)
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
            {
                v2p[partitioning.p2v[pi][i]] = pi;
                for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[partitioning.p2v[pi][i]].size(); j++)
                    f2p[mesh().m_vertex_to_triangle_map[partitioning.p2v[pi][i]][j]] = pi;
            }
        
        // rebuild the face indices
        for (size_t pi = 0; pi < partitioning.p2f.size(); pi++)
            partitioning.p2f[pi].clear();
        for (size_t i = 0; i < nf(); i++)
        {
            assert(f2p[i] >= 0);    // if a face is deleted, its partition should have alrady been emptied, so this shouldn't happen
            partitioning.p2f[f2p[i]].push_back(i);
        }
        
        partitioning.v2p = v2p;
        partitioning.f2p = f2p;
    }
        
    if (deletion_happened || merging_happened)
    {
        // merging doesn't require rebuilding p2v, p2f, v2p, f2p like deletion does, but both require rebuilding the flattened array (flattened_partition_vertices) and its inverse map (indices_in_partitions)
        partitioning.flattened_partition_vertices.clear();
        partitioning.flattened_partition_vertices.reserve(nv());
        partitioning.indices_in_partitions.clear();
        partitioning.indices_in_partitions.reserve(nv());
        for (size_t pi = 0; pi < partitioning.p2v.size(); pi++)
        {
            partitioning.flattened_partition_vertices.insert(partitioning.flattened_partition_vertices.end(), partitioning.p2v[pi].begin(), partitioning.p2v[pi].end());
            for (size_t pvi = 0; pvi < partitioning.p2v[pi].size(); pvi++)
                partitioning.indices_in_partitions.push_back(pvi);
        }
        assert(partitioning.flattened_partition_vertices.size() == nv());
        assert(partitioning.indices_in_partitions.size() == nv());
    }
    
    std::cout << nv() << " vertices, " << nf() << " faces, " << partitioning.p2v.size() << " partitions" << std::endl;
    CSim::TimerMan::timer("Sim.step/BPS.step/partition").stop();
    }
    
    
    
    m_intermediate_v.clear();
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(velv(), "Beginning of step"));

    if (substep == 0 || substep == 4)
    {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////
    ////    time stepping: from here till the end of the step, the mesh connectivity won't change
    ////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (simOptions().implicit)
        step_implicit(dt, partitioning);
    else
        step_explicit(dt, partitioning);
    }
    
    
    
    if (substep == 0 || substep == 5)
    {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////
    ////    wrap up
    ////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
    
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(velv(), "End of step"));

    CSim::TimerMan::timer("Sim.step/BPS.step").stop();
    return actual_dt;
}

double BPS3D::simTime() const
{
    return sim()->time();
}

LosTopos::Vec3st BPS3D::getShuffledTriangle(const LosTopos::Vec3st & t, size_t vertex_to_be_front) const
{
    assert(t[0] == vertex_to_be_front || t[1] == vertex_to_be_front || t[2] == vertex_to_be_front);
    
    LosTopos::Vec3st nt = t;
    while (nt[0] != vertex_to_be_front)
    {
        size_t tmp = nt[0];
        nt[0] = nt[1];
        nt[1] = nt[2];
        nt[2] = tmp;
    }
    assert(nt[0] == vertex_to_be_front);
    
    return nt;
}

Mat3d BPS3D::getVertexPositions(const LosTopos::Vec3st & t) const
{
    Mat3d x;
    x.col(0) = pos(t[0]);
    x.col(1) = pos(t[1]);
    x.col(2) = pos(t[2]);
    return x;
}

Mat3d BPS3D::getVertexVelocities(const LosTopos::Vec3st & t) const
{
    Mat3d v;
    v.col(0) = vel(t[0]);
    v.col(1) = vel(t[1]);
    v.col(2) = vel(t[2]);
    return v;
}

LosTopos::Vec2st BPS3D::getShuffledEdge(const LosTopos::Vec2st & e, size_t vertex_to_be_front) const
{
    return (e[0] == vertex_to_be_front ? e : LosTopos::Vec2st(e[1], e[0]));
}

Mat32d BPS3D::getVertexPositions(const LosTopos::Vec2st & e) const
{
    Mat32d x;
    x.col(0) = pos(e[0]);
    x.col(1) = pos(e[1]);
    return x;
}

Mat32d BPS3D::getVertexVelocities(const LosTopos::Vec2st & e) const
{
    Mat32d v;
    v.col(0) = vel(e[0]);
    v.col(1) = vel(e[1]);
    return v;
}

double BPS3D::vert_integral_mean_curvature(size_t v) const
{
    double H = 0;
    for (size_t i = 0; i < mesh().m_vertex_to_edge_map[v].size(); i++)
    {
        size_t e = mesh().m_vertex_to_edge_map[v][i];
        assert(mesh().m_edge_to_triangle_map[e].size() == 2);
        size_t f0 = mesh().m_edge_to_triangle_map[e][0];
        size_t f1 = mesh().m_edge_to_triangle_map[e][1];
        bool o0 = (mesh().oriented(mesh().m_edges[e][0], mesh().m_edges[e][1], mesh().m_tris[f0]) == mesh().get_triangle_label(f0)[0] > mesh().get_triangle_label(f0)[1]);
        bool o1 = (mesh().oriented(mesh().m_edges[e][0], mesh().m_edges[e][1], mesh().m_tris[f1]) == mesh().get_triangle_label(f1)[0] > mesh().get_triangle_label(f1)[1]);
        assert(o0 != o1);
        if (o1)
            std::swap(f0, f1),
            std::swap(o0, o1);
        assert(o0);
        assert(!o1);
        Vec3d n0 = face_outward_normal(f0);
        Vec3d n1 = face_outward_normal(f1);
        
        double angle = angleAroundAxis(n0, n1, edge_tangent(e));
        
        H += angle * edge_length(e);
    }
    
    return H / 2;
}

double BPS3D::vert_interior_solid_angle(size_t v) const
{
//    return vert_interior_solid_angle_tet_decomposition(v);
    return vert_interior_solid_angle_dihedral(v);
}

double BPS3D::vert_interior_solid_angle_dihedral(size_t v) const
{
    // compute the solid angle of a valence-n vertex by summing the incident edge dihedral angles and subtracting (n-2) * pi
    // reference: unknown
    double sa = 0;
    for (size_t i = 0; i < mesh().m_vertex_to_edge_map[v].size(); i++)
    {
        size_t ei = mesh().m_vertex_to_edge_map[v][i];
        assert(mesh().m_edge_to_triangle_map[ei].size() == 2);
        LosTopos::Vec2st e = mesh().m_edges[ei];
        Vec3d et = edge_tangent(ei);
        
        size_t f0 = mesh().m_edge_to_triangle_map[ei][0];
        size_t f1 = mesh().m_edge_to_triangle_map[ei][1];
        bool o0 = (mesh().oriented(e[0], e[1], mesh().m_tris[f0]) == (mesh().get_triangle_label(f0)[0] > mesh().get_triangle_label(f0)[1]));
        bool o1 = (mesh().oriented(e[0], e[1], mesh().m_tris[f1]) == (mesh().get_triangle_label(f1)[0] > mesh().get_triangle_label(f1)[1]));
        assert(o0 != o1);
        if (o1)
            std::swap(f0, f1), std::swap(o0, o1);
        assert(o0 && !o1);
        
        Vec3d n0 = face_outward_normal(f0);
        Vec3d n1 = face_outward_normal(f1);
        double dihedral_angle = M_PI + (n0.cross(n1).dot(et) > 0 ? -1 : 1) * acos(std::min(1.0, std::max(-1.0, n0.dot(n1))));
        
        sa += dihedral_angle;
    }
    
    sa -= (mesh().m_vertex_to_edge_map[v].size() - 2) * M_PI;
    
    return sa;
}

double BPS3D::vert_interior_solid_angle_tet_decomposition(size_t v) const
{
    // compute the solid angle at the vertex v by decomposing the 1-ring polygon into triangles (which may subtend a tetrahedron with positive or
    //  negative determinant with repsect to vertex v), and summing the (signed) solid angles subtended by each triangle using the formula by
    //      Oosterom 1983, The Solid Angle of A Plane Triangle, http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=4121581
    //      See also Wikipedia: https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron

    // sort the 1-ring neighbor vertices first (make the vertices consecutive and the order CW when looking at v from outside)
    std::vector<size_t> theonering;

    assert(mesh().m_vertex_to_edge_map[v].size() > 2);
    LosTopos::Vec3st t0 = getShuffledTriangle(mesh().m_tris[mesh().m_vertex_to_triangle_map[v][0]], v);
    if (mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[v][0])[0] > mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[v][0])[1])
        std::swap(t0[1], t0[2]);
    theonering.push_back(t0[1]);
    theonering.push_back(t0[2]);
    while (theonering.size() < mesh().m_vertex_to_edge_map[v].size())
    {
        for (size_t i = 0; i < mesh().m_vertex_to_triangle_map[v].size(); i++)
        {
            LosTopos::Vec3st t = mesh().m_tris[mesh().m_vertex_to_triangle_map[v][i]];
            assert(mesh().triangle_contains_vertex(t, v));
            if (mesh().triangle_contains_vertex(t, theonering.back()) && !mesh().triangle_contains_vertex(t, *(theonering.rbegin() + 1)))
            {
                theonering.push_back(mesh().get_third_vertex(v, theonering.back(), t));
                break;
            }
        }
    }

    for (size_t i = 0; i < theonering.size(); i++)
    {
        assert(mesh().get_triangle_index(v, theonering[i], theonering[(i + 1) % theonering.size()]) < nf());
    }

    // sum the solid angle subtended by the triangle formed theonering[0] and each consecutive pair of vertices in theonering
    double sa = 0;
    Vec3d o = pos(v);
    for (size_t i = 1; i + 1 < theonering.size(); i++)
    {
        Vec3d va = pos(theonering[0])     - o;
        Vec3d vb = pos(theonering[i])     - o;
        Vec3d vc = pos(theonering[i + 1]) - o;
        double a = va.norm();
        double b = vb.norm();
        double c = vc.norm();
        
        // formula from https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
        // since the vertices in 1-ring are ordered CW when looking at v from outside, the triple product a cross b dot c should be positive when the vertex is convex (solid angle < 2 pi)
        // Note that the formula given by Oosterom 1983 only works when the solid angle is smaller than pi because they assume cos(omega/2) >= 0 when deriving Eq. 7 (their application
        //  only requires computing tiny triangles that are usually far away from the origin). The following code tries to account for all the other cases as well.
        double numerator = va.cross(vb).dot(vc);
        double denominator = (a * b * c + va.dot(vb) * c + vb.dot(vc) * a + vc.dot(va) * b);
        double tanhalfomega = numerator / denominator;
        double omega = 0;
        
        // the following handling produces a solid angle between -2 pi and 2 pi, where concave angles are treated as negative (instead of between 2 pi and 4 pi), so that they can be properly summed
        if         (denominator >= 0 && numerator >= 0)  omega =             atan(tanhalfomega) * 2;  // the original case (the triangle is far away): omega between 0 and pi
        else if    (denominator <  0 && numerator >= 0)  omega =  2 * M_PI + atan(tanhalfomega) * 2;  // the tet is still positively oriented, but the triangle is close to the origin, subtending a large solid angle: omega between pi and 2 pi
        else if    (denominator <  0 && numerator <  0)  omega = -2 * M_PI + atan(tanhalfomega) * 2;  // the tet is negatively oriented, and the triangle is close to the origin: omega between -2 pi and -pi
        else assert(denominator >= 0 && numerator <  0), omega =           + atan(tanhalfomega) * 2;  // the tet is negatively oriented, and the triangle is far away: omega between -pi and 0
        if (denominator == 0) omega = (numerator >= 0 ? M_PI : -M_PI);   // the corner cases where cos(omega/2) is zero
        
        sa += omega;
    }
    
    if (sa < 0) // convert negative angles (concave) to the positive counterpart
        sa += 4 * M_PI;
    
    return sa;
}

double BPS3D::face_interior_angle(size_t f, size_t v) const
{
    assert(mesh().triangle_contains_vertex(mesh().m_tris[f], v));
    LosTopos::Vec3st t = getShuffledTriangle(mesh().m_tris[f], v);
    
    double as = (pos(t[1]) - pos(t[0])).squaredNorm();
    double bs = (pos(t[2]) - pos(t[0])).squaredNorm();
    double cs = (pos(t[2]) - pos(t[1])).squaredNorm();
    double a = sqrt(as);
    double b = sqrt(bs);
    double cosa = (as + bs - cs) / (2 * a * b);
    if (cosa < -1) cosa = -1;
    if (cosa > 1)  cosa = 1;
    
    return acos(cosa);
}

double BPS3D::liquid_volume() const
{
    double volume = 0;
    Vec3d xref(0, 0, 0);
    for (size_t i = 0; i < mesh().nt(); i++)
    {
        LosTopos::Vec3st t = mesh().get_triangle(i);
        Vec3d x0 = pos(t[0]);
        Vec3d x1 = pos(t[1]);
        Vec3d x2 = pos(t[2]);
        double v = (x0 - xref).cross(x1 - xref).dot(x2 - xref) / 6;
        
        LosTopos::Vec2i l = mesh().get_triangle_label(i);
        volume += v * (l[0] > l[1] ? 1 : -1);
    }
    
    return volume;
}

double BPS3D::liquid_volume_change_rate() const
{
    return liquid_volume_change_rate(velv());
}

double BPS3D::liquid_volume_change_rate(const VecXd & v) const
{
    double dvdt = 0;
    Vec3d xref(0, 0, 0);
    for (size_t i = 0; i < mesh().nt(); i++)
    {
        LosTopos::Vec3st t = mesh().get_triangle(i);
        double a = face_area(i);
        Vec3d n = face_outward_normal(i);
        
        double dv0dt = v.segment<3>(t[0] * 3).dot(n) * a / 3;
        double dv1dt = v.segment<3>(t[1] * 3).dot(n) * a / 3;
        double dv2dt = v.segment<3>(t[2] * 3).dot(n) * a / 3;
        
        dvdt += dv0dt + dv1dt + dv2dt;
    }
    
    return dvdt;
}

double BPS3D::partition_volume(const std::vector<size_t> & faces) const
{
    double volume = 0;
    Vec3d xref(0, 0, 0);
    for (size_t i = 0; i < faces.size(); i++)
    {
        LosTopos::Vec3st t = mesh().m_tris[faces[i]];
        Vec3d x0 = pos(t[0]);
        Vec3d x1 = pos(t[1]);
        Vec3d x2 = pos(t[2]);
        double v = (x0 - xref).cross(x1 - xref).dot(x2 - xref) / 6;
        
        LosTopos::Vec2i l = mesh().get_triangle_label(faces[i]);
        volume += v * (l[0] > l[1] ? 1 : -1);
    }
    
    return volume;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Callbacks
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool BPS3D::generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos)
{
    if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
    {
        if (vertex_is_on_triple_junction(v0) && vertex_is_on_triple_junction(v1))
        {
            pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        } else if (vertex_is_on_triple_junction(v0))
        {
            pos = st.pm_positions[v0];
        } else if (vertex_is_on_triple_junction(v1))
        {
            pos = st.pm_positions[v1];
        } else
        {
            pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        }
    } else if (st.vertex_is_any_solid(v0))
    {
        pos = st.pm_positions[v0];
    } else if (st.vertex_is_any_solid(v1))
    {
        pos = st.pm_positions[v1];
    } else
    {
        pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    }
    
    return true;
}

bool BPS3D::generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos)
{
    callback_log() << "solid callback: generate split position: " << v0 << " " << v1 << " " << (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1)) << std::endl;
    pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    return true;
}

bool BPS3D::generate_snapped_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos)
{
    bool v0solid = st.vertex_is_any_solid(v0);
    bool v1solid = st.vertex_is_any_solid(v1);

    // scan the neighborhood right now to prune out the case of a singular solid vertex with no incident solid faces. This can happen as T1 pulls apart a solid vertex which is a singular vertex shared by a free surface and a solid surface (e.g. the middle of the thin sheet in floor splash; see 1450021611).
    bool v0_solid_face = false;
    for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[v0].size(); j++)
        if (face_is_solid(mesh().m_vertex_to_triangle_map[v0][j]))
            v0_solid_face = true;
    if (!v0_solid_face)
        v0solid = false;   // a solid vertex incident to no solid faces is not a true solid vertex.
    
    bool v1_solid_face = false;
    for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[v1].size(); j++)
        if (face_is_solid(mesh().m_vertex_to_triangle_map[v1][j]))
            v1_solid_face = true;
    if (!v1_solid_face)
        v1solid = false;   // a solid vertex incident to no solid faces is not a true solid vertex.
    
    if (v0solid && v1solid)
    {
        if (vertex_is_on_triple_junction(v0) && vertex_is_on_triple_junction(v1))
        {
            pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        } else if (vertex_is_on_triple_junction(v0))
        {
            pos = st.pm_positions[v0];
        } else if (vertex_is_on_triple_junction(v1))
        {
            pos = st.pm_positions[v1];
        } else
        {
            pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        }
    } else if (v0solid)
    {
        pos = st.pm_positions[v0];
    } else if (v1solid)
    {
        pos = st.pm_positions[v1];
    } else
    {
        pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    }
    
    return true;
}

LosTopos::Vec3c BPS3D::generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1)
{
    return LosTopos::Vec3c(label0[0] || label1[0], label0[1] || label1[1], label0[2] || label1[2]);
}

LosTopos::Vec3c BPS3D::generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1)
{
    return LosTopos::Vec3c(label0[0] && label1[0], label0[1] && label1[1], label0[2] && label1[2]);
}

LosTopos::Vec3c BPS3D::generate_snapped_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1)
{
    return LosTopos::Vec3c(label0[0] || label1[0], label0[1] || label1[1], label0[2] || label1[2]);
}

bool BPS3D::generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower)
{
    return false;
}

bool BPS3D::generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b)
{
    double floorz = Options::doubleValue("floor-z");
    if (pos_a[2] < floorz) pos_a[2] = floorz;
    if (pos_b[2] < floorz) pos_b[2] = floorz;
    
    return true;
}

bool BPS3D::solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e)
{
    return edge_is_on_triple_junction(e);
}

LosTopos::Vec3d BPS3D::sampleVelocity(LosTopos::Vec3d & pos)
{
    return LosTopos::Vec3d(0, 0, 0);
}

bool BPS3D::sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output)
{
    return false;
}

struct CollapseTempData
{
    size_t v0;
    size_t v1;
    
    Vec3d old_x0;
    Vec3d old_x1;

    Vec3d old_u0;
    Vec3d old_u1;
};

void BPS3D::pre_collapse(const LosTopos::SurfTrack & st, size_t e, void ** data)
{
    CollapseTempData * td = new CollapseTempData;
    td->v0 = st.m_mesh.m_edges[e][0];
    td->v1 = st.m_mesh.m_edges[e][1];
    
    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);
    
    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);
    
    *data = (void *)td;
    callback_log() << "pre collapse: " << e << ": " << td->v0 << " " << td->v1 << std::endl;
}

void BPS3D::post_collapse(const LosTopos::SurfTrack & st, size_t e, size_t merged_vertex, void * data)
{
    CollapseTempData * td = (CollapseTempData *)data;
    callback_log() << "post collapse: " << e << ": " << td->v0 << " " << td->v1 << " => " << merged_vertex << std::endl;
    assert((st.m_mesh.vertex_is_deleted(td->v0) && merged_vertex == td->v1) || (st.m_mesh.vertex_is_deleted(td->v1) && merged_vertex == td->v0));

    Vec3d merged_x = vc(st.pm_positions[merged_vertex]);
    double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();
    if (s > 1) s = 1;
    if (s < 0) s = 0;
    Vec3d new_u = td->old_u0 * (1 - s) + td->old_u1 * s;
    
    vel(merged_vertex) = new_u;
}

struct SplitTempData
{
    size_t v0;
    size_t v1;

    Vec3d old_x0;
    Vec3d old_x1;
    
    Vec3d old_u0;
    Vec3d old_u1;
};

void BPS3D::pre_split(const LosTopos::SurfTrack & st, size_t e, void ** data)
{
    SplitTempData * td = new SplitTempData;
    td->v0 = st.m_mesh.m_edges[e][0];
    td->v1 = st.m_mesh.m_edges[e][1];
    
    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);
    
    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);
    
    *data = (void *)td;
    callback_log() << "pre split: " << e << ": " << td->v0 << " " << td->v1 << std::endl;
}

void BPS3D::post_split(const LosTopos::SurfTrack & st, size_t e, size_t new_vertex, void * data)
{
    SplitTempData * td = (SplitTempData *)data;
    callback_log() << "post split: " << e << ": " << td->v0 << " " << td->v1 << " => " << new_vertex << std::endl;
    
    Vec3d midpoint_x = vc(st.pm_positions[new_vertex]);
    double s = (midpoint_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();
    if (s > 1) s = 1;
    if (s < 0) s = 0;
    Vec3d new_u = td->old_u0 * (1 - s) + td->old_u1 * s;
    
    vel(new_vertex) = new_u;
}

void BPS3D::pre_flip(const LosTopos::SurfTrack & st, size_t e, void ** data)
{
    
}

void BPS3D::post_flip(const LosTopos::SurfTrack & st, size_t e, void * data)
{
    
}

struct T1TempData
{

};

void BPS3D::pre_t1(const LosTopos::SurfTrack & st, size_t v, void ** data)
{
    callback_log() << "pre t1: " << v << std::endl;
}

void BPS3D::post_t1(const LosTopos::SurfTrack & st, size_t v, size_t a, size_t b, void * data)
{
    callback_log() << "post t1: " << v << " => " << a << " " << b << std::endl;
    
    vel(a) = vel(v);
    vel(b) = vel(v);
}

struct FaceSplitTempData
{
    size_t v0;
    size_t v1;
    size_t v2;
    
    Vec3d old_x0;
    Vec3d old_x1;
    Vec3d old_x2;
    
    Vec3d old_u0;
    Vec3d old_u1;
    Vec3d old_u2;
};

void BPS3D::pre_facesplit(const LosTopos::SurfTrack & st, size_t f, void ** data)
{
    FaceSplitTempData * td = new FaceSplitTempData;
    
    td->v0 = st.m_mesh.get_triangle(f)[0];
    td->v1 = st.m_mesh.get_triangle(f)[1];
    td->v2 = st.m_mesh.get_triangle(f)[2];

    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);
    td->old_x2 = vc(st.pm_positions[td->v2]);

    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);
    td->old_u2 = vel(td->v2);
    
    *data = (void *)td;
    callback_log() << "pre facesplit: " << f << ": " << td->v0 << " " << td->v1 << " " << td->v2 << std::endl;
}

void BPS3D::post_facesplit(const LosTopos::SurfTrack & st, size_t f, size_t new_vertex, void * data)
{
    FaceSplitTempData * td = (FaceSplitTempData *)data;
    callback_log() << "post facesplit: " << f << " => " << new_vertex << std::endl;

    Vec3d new_x = vc(st.pm_positions[new_vertex]);
    Vec3d c = Vec3d::Zero();
    Vec3d n = (td->old_x1 - td->old_x0).cross(td->old_x2 - td->old_x0);
    double nsq = n.squaredNorm();
    c[0] = 1 - (new_x - td->old_x0).dot(n.cross(td->old_x1 - td->old_x2)) / nsq;
    c[1] = 1 - (new_x - td->old_x1).dot(n.cross(td->old_x2 - td->old_x0)) / nsq;
    if (c[0] > 1)        c[0] = 1;
    if (c[0] < 0)        c[0] = 0;
    if (c[1] > 1 - c[0]) c[1] = 1 - c[0];
    if (c[1] < 0)        c[1] = 0;
    c[2] = 1 - c[0] - c[1];
    
    vel(new_vertex) = td->old_u0 * c[0] + td->old_u1 * c[1] + td->old_u2 * c[2];
}

struct SnapTempData
{
    size_t v0;
    size_t v1;
    
    Vec3d old_x0;
    Vec3d old_x1;

    Vec3d old_u0;
    Vec3d old_u1;
};

void BPS3D::pre_snap(const LosTopos::SurfTrack & st, size_t v0, size_t v1, void ** data)
{
    SnapTempData * td = new SnapTempData;
    td->v0 = v0;
    td->v1 = v1;
    
    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);
    
    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);
    
    *data = (void *)td;
    callback_log() << "pre snap: " << td->v0 << " " << td->v1 << std::endl;
}

void BPS3D::post_snap(const LosTopos::SurfTrack & st, size_t v_kept, size_t v_deleted, void * data)
{
    SnapTempData * td = (SnapTempData *)data;
    callback_log() << "post snap: " << td->v0 << " " << td->v1 << " => " << v_kept << std::endl;
    assert((td->v0 == v_kept && td->v1 == v_deleted) || (td->v1 == v_kept && td->v0 == v_deleted));
    assert(v_kept != v_deleted);
//    assert(st.m_mesh.vertex_is_deleted(v_deleted));
//    assert(!st.m_mesh.vertex_is_deleted(v_kept));
    
    Vec3d merged_x = vc(st.pm_positions[v_kept]);
    double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();
    if (s > 1) s = 1;
    if (s < 0) s = 0;
    Vec3d new_u = td->old_u0 * (1 - s) + td->old_u1 * s;
    
    vel(v_kept) = new_u;
}

void BPS3D::pre_smoothing(const LosTopos::SurfTrack & st, void ** data)
{
    
}

void BPS3D::post_smoothing(const LosTopos::SurfTrack & st, void * data)
{
    
}
