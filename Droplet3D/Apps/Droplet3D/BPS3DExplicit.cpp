//
//  BPS3DExplicit.cpp
//  Droplet3D
//
//  Created by Fang Da on 10/3/15.
//
//

#include "BPS3D.h"
#include "SimOptions.h"

#ifndef _OPENMP
#define _OPENMP
#endif

#define SUMMATION_TECHNIQUE_DIRECT  0
#define SUMMATION_TECHNIQUE_FMM     1

#define SUMMATION_TECHNIQUE         (SUMMATION_TECHNIQUE_DIRECT)

void BPS3D::step_explicit(double dt, const Partitioning & partitioning)
{
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit").start();
    std::cout << "Explicit time stepping" << std::endl;
    
    
    
    // prepare a convenient data structure, taking advantage of the fact that the mesh does not change within such a step
    VecXd v = VecXd::Zero(nv() * 3);
    for (size_t i = 0; i < nv(); i++)
        v.segment<3>(i * 3) = (*m_v)[i];
    
    
    
    // Helmholtz Decomposition to purify the velocity field
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD").start();
    VecXd oldv = v;
    bool HD = false;
    static int s_counter = 0;
    s_counter++;
    if (Options::intValue("hd-interval") == 0 || simTime() < Options::doubleValue("hd-interval-start") || s_counter == Options::intValue("hd-interval"))
    {
#if (SUMMATION_TECHNIQUE == SUMMATION_TECHNIQUE_DIRECT)
        step_HelmholtzDecomposition(dt, v, partitioning);
#elif (SUMMATION_TECHNIQUE == SUMMATION_TECHNIQUE_FMM)
        step_HD_FMM(dt, v, partitioning);
#endif
        HD = true;
        s_counter = 0;
    }
    assert(v == v); // assert HD doesn't produce NaNs
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD").stop();
    
    
    
    // velocity smoothing (also sharpening)
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/smooth").start();
    double smoothing_coef = simOptions().smoothing_coef * dt;
    if (HD) smoothing_coef -= 0.5;
    step_SmoothVelocity(dt, v, partitioning, smoothing_coef);
    
    // CSOC: clamp the inward normal velocity on concave vertices to the pre-HD level if the inward motion became stronger due to the comebined effect of HD and smoothing.
    if (Options::boolValue("csoc"))
    {
        for (size_t i = 0; i < nv(); i++)
        {
            if (vert_interior_solid_angle(i) > 2.75 * M_PI)
            {
                Vec3d n_i = vert_outward_normal(i);
                Vec3d ovi = oldv.segment<3>(i * 3);
                Vec3d vi = v.segment<3>(i * 3);
                Vec3d dvi = vi - ovi;
                if (vi.dot(n_i) < 0 && dvi.dot(n_i) < 0)    // the velocity is both inward and increasing
                {
//                    double normal_velocity_delta = (std::min(ovi.dot(n_i), 0.0) - vi.dot(n_i));
                    double normal_velocity_delta = (ovi.dot(n_i) - vi.dot(n_i));    // ignore the clamp with zero: if all the neighbors of this concave point is advancing forward, clamping this concave point to have zero normal velocity is not good enough.
                    assert(normal_velocity_delta >= 0);
                    std::cout << "CSOC: vertex " << i << " delta = " << normal_velocity_delta << std::endl;
                    v.segment<3>(i * 3) = vi + normal_velocity_delta * n_i;
                }
            }
        }
    }
    
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v, "smoothing"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v - oldv, "change due to HD + smoothing"));
    assert(v == v); // assert smoothing doesn't produce NaNs
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/smooth").stop();


    
    // Add gravity
    for (size_t i = 0; i < nv(); i++)
        v.segment<3>(i * 3) += simOptions().gravity * dt;
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v, "gravity"));
    
    
    
    // Pressure solve taking into account solid contact and surface tension
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM").start();
#if (SUMMATION_TECHNIQUE == SUMMATION_TECHNIQUE_DIRECT)
    step_PressureSolve(dt, v, partitioning);
#elif (SUMMATION_TECHNIQUE == SUMMATION_TECHNIQUE_FMM)
    step_PressureSolveIterative(dt, v, partitioning);
#endif
    assert(v == v); // assert the pressure solve doesn't produce NaNs
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM").stop();

    
    
    // update the velocity field with the result v
    for (size_t i = 0; i < nv(); i++)
        (*m_v)[i] = v.segment<3>(i * 3);
    
    
    
    std::cout << "Explicit time stepping finished" << std::endl;
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit").stop();
}

void BPS3D::step_SmoothVelocity(double dt, VecXd & v, const Partitioning & partitioning, double coef)
{
    // velocity smoothing (implicit)
    // u_new = u_old + smoothing_coef * dt * Laplacian u_new
    // the equation is constructed as if u_new and u_old are scalar fields, and solved three times, each with a different coordinate component as a different rhs
    // note that the Laplacian discretization is consistent with the smoothing introduced by HD, so the effect should be exactly canceling if smoothing_coef * dt = -0.5 here.
    // note also that although the incoming arguments include a partitioning, the implementation below ignores it and does a global solve including all connected components of the mesh.
    Eigen::SparseMatrix<double> smoothing_A(nv(), nv());
    Eigen::Matrix<double, Eigen::Dynamic, 3> smoothing_rhs(nv(), 3);   // the three columns are for x, y and z components respectively
    smoothing_rhs.setZero();
    assert(smoothing_rhs.rows() == nv());
    assert(smoothing_rhs.cols() == 3);
    
    std::vector<Eigen::Triplet<double> > smoothing_A_triplets;
    for (size_t i = 0; i < nv(); i++)
    {
        double ai = vert_area(i) * 3;
        for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
        {
            size_t e = mesh().m_vertex_to_edge_map[i][j];
            size_t vother = edge_other_vertex(e, i);
            
            assert(mesh().m_edge_to_triangle_map[e].size() == 2);
            double aj = 0;
            for (size_t k = 0; k < mesh().m_edge_to_triangle_map[e].size(); k++)
                aj += face_area(mesh().m_edge_to_triangle_map[e][k]);
            double w = aj / (ai * 2);
            
            smoothing_A_triplets.push_back(Eigen::Triplet<double>(i, vother, -w * coef));
        }
        
        smoothing_A_triplets.push_back(Eigen::Triplet<double>(i, i, 1 + coef));
        
        smoothing_rhs.row(i) = v.segment<3>(i * 3).transpose();
    }
    
    smoothing_A.setFromTriplets(smoothing_A_triplets.begin(), smoothing_A_triplets.end());
    
    VecXd newv = VecXd::Zero(nv() * 3);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > smoothing_solver(smoothing_A);
    VecXd smoothing_solution =  smoothing_solver.solve(smoothing_rhs.col(0));    // the x component of u_new
    for (size_t i = 0; i < nv(); i++) newv[i * 3 + 0] = smoothing_solution[i];
    smoothing_solution =        smoothing_solver.solve(smoothing_rhs.col(1));    // the y component of u_new
    for (size_t i = 0; i < nv(); i++) newv[i * 3 + 1] = smoothing_solution[i];
    smoothing_solution =        smoothing_solver.solve(smoothing_rhs.col(2));    // the y component of u_new
    for (size_t i = 0; i < nv(); i++) newv[i * 3 + 2] = smoothing_solution[i];
    
    v = newv;
}

void BPS3D::partition_mesh(Partitioning & partitioning)
{
    std::vector<std::vector<size_t> > & p2v = partitioning.p2v;
    std::vector<std::vector<size_t> > & p2f = partitioning.p2f;
    std::vector<int>                  & v2p = partitioning.v2p;
    std::vector<int>                  & f2p = partitioning.f2p;
    std::vector<size_t> & flattened_partition_vertices = partitioning.flattened_partition_vertices;
    std::vector<size_t> & indices_in_partitions        = partitioning.indices_in_partitions;

    p2v.clear();
    p2f.clear();
    v2p.assign(nv(), -1);
    f2p.assign(nf(), -1);
    
    // build the partition of vertices first
    while (true)
    {
        size_t seed_vertex = nv();
        for (size_t i = 0; i < nv(); i++)
            if (v2p[i] < 0)
            {
                seed_vertex = i;
                break;
            }
        if (seed_vertex == nv())    // can't find a vertex that doesn't belong to any partition. we're done.
            break;
        
        int new_partition_id = (int)p2v.size();
        p2v.push_back(std::vector<size_t>());
        
        std::vector<size_t> openset;
        openset.push_back(seed_vertex);
        while (openset.size() > 0)    // DPS to traverse this paritition, finding all its vertices
        {
            size_t current = openset.back();
            openset.pop_back();
            
            if (v2p[current] >= 0)
            {
                assert(v2p[current] == new_partition_id);
                continue;
            }
            
            v2p[current] = new_partition_id;
            p2v.back().push_back(current);
            
            for (size_t i = 0; i < mesh().m_vertex_to_edge_map[current].size(); i++)
            {
                LosTopos::Vec2st e = mesh().m_edges[mesh().m_vertex_to_edge_map[current][i]];
                size_t vother = (e[0] == current ? e[1] : e[0]);
                
                openset.push_back(vother);
            }
        }
    }
    
    // sanity check on the vertex counts
    assert(p2v.size() > 0);
    size_t partition_vertex_sum = 0;
    for (size_t i = 0; i < p2v.size(); i++)
        partition_vertex_sum += p2v[i].size();
    assert(partition_vertex_sum == nv());

    // build the partition of faces from the partition of vertices
    p2f.resize(p2v.size());
    for (size_t i = 0; i < nf(); i++)
    {
        LosTopos::Vec3st t = mesh().m_tris[i];
        int p = v2p[t[0]];
        assert(p == v2p[t[1]]);
        assert(p == v2p[t[2]]);
        assert(p >= 0);
        f2p[i] = p;
        p2f[p].push_back(i);
    }
    
    // sanity check on the face counts
    assert(p2f.size() > 0);
    size_t partition_face_sum = 0;
    for (size_t i = 0; i < p2f.size(); i++)
        partition_face_sum += p2f[i].size();
    assert(partition_face_sum == nf());
    
    for (size_t i = 0; i < nv(); i++)
        assert(v2p[i] >= 0);
    for (size_t i = 0; i < nf(); i++)
        assert(f2p[i] >= 0);
    
    // build the flattened array of vertices which helps with parallelization
    flattened_partition_vertices.clear();
    flattened_partition_vertices.reserve(nv());
    indices_in_partitions.clear();
    indices_in_partitions.reserve(nv());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        flattened_partition_vertices.insert(flattened_partition_vertices.end(), p2v[pi].begin(), p2v[pi].end());
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
            indices_in_partitions.push_back(pvi);
    }
    assert(flattened_partition_vertices.size() == nv());
    assert(indices_in_partitions.size() == nv());

}

