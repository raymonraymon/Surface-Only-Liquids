//
//  BPS3DPressureSolve.cpp
//  Droplet3D
//
//  Created by Fang Da on 10/8/15.
//
//
//#include <omp.h>
#include "BPS3D.h"
#include "BoundaryIntegral.h"
#include "SimOptions.h"
#include "Sim.h"

void BPS3D::computeCollocCoefficients(std::vector<MatXd> & A_solid, std::vector<MatXd> & A_air, std::vector<MatXd> & B_solid, std::vector<MatXd> & B_air, const Partitioning & partitioning)
{
    const std::vector<std::vector<size_t> > & p2v = partitioning.p2v;
    const std::vector<std::vector<size_t> > & p2f = partitioning.p2f;
    const std::vector<int>                  & v2p = partitioning.v2p;
    const std::vector<int>                  & f2p = partitioning.f2p;
    const std::vector<size_t> & flattened_partition_vertices = partitioning.flattened_partition_vertices;
    const std::vector<size_t> & indices_in_partitions        = partitioning.indices_in_partitions;
    
    A_solid.resize(p2v.size());
    A_air  .resize(p2v.size());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        A_solid[pi] = MatXd::Zero(p2v[pi].size(), p2v[pi].size());
        A_air  [pi] = MatXd::Zero(p2v[pi].size(), p2v[pi].size());
    }
    
    B_solid.resize(p2v.size());
    B_air  .resize(p2v.size());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        B_solid[pi] = MatXd::Zero(p2v[pi].size(), p2v[pi].size());
        B_air  [pi] = MatXd::Zero(p2v[pi].size(), p2v[pi].size());
    }
    
    const double oneover4pi = 1 / (4 * M_PI);

    const std::vector<Vec3d> quadrature_square =   BoundaryIntegral::quadrature_square();

    std::vector<bool> face_is_solids(nf(), false);
    std::vector<Vec3d> face_outward_normals(nf(), Vec3d::Zero());
    std::vector<double> face_areas(nf(), 0);
    for (size_t i = 0; i < nf(); i++)
    {
        face_is_solids[i] = face_is_solid(i);
        face_outward_normals[i] = face_outward_normal(i);
        face_areas[i] = face_area(i);
    }
    
    std::vector<Vec3d> vert_solid_normals(nv(), Vec3d::Zero());
    std::vector<Vec3d> vert_air_normals(nv(), Vec3d::Zero());
    std::vector<double> vert_solid_areas(nv(), 0);
    std::vector<double> vert_air_areas(nv(), 0);
    for (size_t i = 0; i < nv(); i++)
    {
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][j];
            Vec3d n = face_outward_normals[face_i];
            double a = face_areas[face_i];
            if (face_is_solids[face_i])
                vert_solid_normals[i] += a * n,
                vert_solid_areas[i] += a;   // don't divide by 3 here, because we're looking at the support of the vertex basis function theta_i
            else
                vert_air_normals[i] += a * n,
                vert_air_areas[i] += a;     // don't divide by 3 here, because we're looking at the support of the vertex basis function theta_i
        }
        if (vert_solid_areas[i] > 0)
            vert_solid_normals[i].normalize();
        if (vert_air_areas[i] > 0)
            vert_air_normals[i].normalize();
    }
    
#pragma omp parallel default (none) shared (vert_solid_normals, vert_air_normals, vert_solid_areas, vert_air_areas, face_areas, face_outward_normals, face_is_solids, A_solid, A_air, B_solid, B_air, p2v, p2f, v2p, f2p, flattened_partition_vertices, indices_in_partitions, std::cout)
{
#pragma omp master
    {
//        std::cout << "BEM collocA and collocB loop: spawning " << omp_get_num_threads() << " threads." << std::endl;
    }
        
#pragma omp for schedule (dynamic, 10)
    for (size_t fpi = 0; fpi < nv(); fpi++)   // degree of freedom on vertex i
    {
        size_t i = flattened_partition_vertices[fpi];
        size_t pvi = indices_in_partitions[fpi];
        size_t pi = v2p[i];
        
        double mel = 0;
        for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
            mel += edge_length(mesh().m_vertex_to_edge_map[i][j]);
        mel /= mesh().m_vertex_to_edge_map[i].size();
        
    }
//#pragma omp for

}
//#pragma omp parallel
        
}

void BPS3D::step_PressureSolve(double dt, VecXd & v, const Partitioning & partitioning)
{
    bool TPCF = Options::boolValue("tpcf");
    bool TDMC = Options::boolValue("tdmc");
    
    const std::vector<std::vector<size_t> > & p2v = partitioning.p2v;
    const std::vector<std::vector<size_t> > & p2f = partitioning.p2f;
    const std::vector<int>                  & v2p = partitioning.v2p;
    const std::vector<int>                  & f2p = partitioning.f2p;
    const std::vector<size_t> & flattened_partition_vertices = partitioning.flattened_partition_vertices;
    const std::vector<size_t> & indices_in_partitions        = partitioning.indices_in_partitions;

    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);

    size_t N = nv();
    
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/BC").start();
    // compute the boundary conditions
    // note that both p and dpdn here are time integrated (multiplied by dt)
    VecXd BC_p = VecXd::Zero(N);       // pressure at vertex (set 0 for solid interior vertices, which will not be used anyway)
    VecXd BC_dpdn = VecXd::Zero(N);    // dpdn at vertex (the solid side of the vertex for triple junction vertices, and 0 for air interior vertices which will not be used anyway)
    VecXd BC_dpdna = VecXd::Zero(N);   // dpdn at triple junction vertex on the air side, in the amount just enough to keep the triple junction stationary (dpdna = v_tangent.dot(na) * rho = ((I - ns ns^T) v).dot(na) * rho)
    double influx = sim()->currentInflux();
    for (size_t i = 0; i < nv(); i++)
    {
        Vec3d n_i = vert_outward_normal(i);
        
        if (!vertex_is_solid(i))
        {
            // vertex i is on the interior of the liquid-air interface
            //  BC_p[i] = the pressure jump due to surface tension; BC_dpdn[i] = 0 (unused)
            BC_p[i] = vert_integral_mean_curvature(i) * simOptions().sigma / vert_area(i) * dt;     // multiplied by dt because the pressure is time integrated
            BC_dpdn[i] = 0;
            BC_dpdna[i] = 0;
            assert(BC_p[i] == BC_p[i]);
            
        } else if (vertex_is_solid(i) && !vertex_is_on_triple_junction(i))
        {
            // vertex i is on the interior of the liquid-solid interface
            //  BC_p[i] = 0 (unused); BC_dpdn[i] = the solid boundary condition
            BC_p[i] = 0;
            BC_dpdn[i] = (v.segment<3>(i * 3).dot(n_i) + influx) * simOptions().rho;
            BC_dpdna[i] = 0;
            
        } else
        {
            // vertex i is on triple junction
            //  BC_p[i] = the pressure jump due to the triple junction surface tension; BC_dpdn[i] = the solid side boundary condition
            // the pressure jump comes from two sources: the curvature along the triple junction tangent direction, and the balance of the
            //  three-phase surface tension forces in the plane orthogonal to the triple junction (which can be seen as the curvature along
            //  the direction orthogonal to the triple junction).
            Vec3d n_solid = Vec3d::Zero();
            Vec3d n_air = Vec3d::Zero();
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                size_t fj = mesh().m_vertex_to_triangle_map[i][j];
                if (face_is_solid(fj))  n_solid += face_interior_angle(fj, i) * face_outward_normal(fj);
                else                    n_air   += face_interior_angle(fj, i) * face_outward_normal(fj);
            }
            n_solid.normalize();
            n_air.normalize();

            Vec3d triple_junction_tangent = n_solid.cross(n_air).normalized();
            Vec3d t_solid_outward = triple_junction_tangent.cross(n_solid).normalized();
            Vec3d t_air_outward =  -triple_junction_tangent.cross(n_air).normalized();
            
            // compute the pressure jump in the direction orthogonal to the triple junction (due to combined three-face surface tension forces)
            Vec3d surface_tension_force_sl = -t_solid_outward * simOptions().sigma * simOptions().sigma_sl;
            Vec3d surface_tension_force_sa =  t_solid_outward * simOptions().sigma * simOptions().sigma_sa;
            Vec3d surface_tension_force_la =   -t_air_outward * simOptions().sigma;
            
            Vec3d surface_tension_force_combined = surface_tension_force_sl + surface_tension_force_sa + surface_tension_force_la;
            
            double pressure_jump_normal_to_triple_junction = surface_tension_force_combined.dot(-t_solid_outward) / triple_junction_virtual_width();  // project the combined force in the plane orthogonal to the triple junction in the direction of the solid tangent, because the component normal to the solid surface will be canceled by the solid normal force

            // compute the pressure jump in the direction tangent to the triple junction (due to the curvature of the triple junction curve)
            double triple_junction_solid_angle = 0;
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
                if (face_is_solid(mesh().m_vertex_to_triangle_map[i][j]))
                    triple_junction_solid_angle += face_interior_angle(mesh().m_vertex_to_triangle_map[i][j], i);
            double triple_junction_curve_integral_curvature = M_PI - triple_junction_solid_angle;
            
            std::vector<size_t> triple_junction_edges;
            for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
                if (edge_is_on_triple_junction(mesh().m_vertex_to_edge_map[i][j]))
                    triple_junction_edges.push_back(mesh().m_vertex_to_edge_map[i][j]);
//            assert(triple_junction_edges.size() == 2);
            
            double pressure_jump_tangent_to_triple_junction = triple_junction_curve_integral_curvature * simOptions().sigma / ((edge_length(triple_junction_edges[0]) + edge_length(triple_junction_edges[1])) / 2);
            
            if (triple_junction_edges.size() != 2)  // this can only happen transiently when a droplet is beginning to come into contact with a solid surface but hasn't fully formed a contiguous contact region yet
            {
                pressure_jump_tangent_to_triple_junction = 0;   // don't do anything in this case because the system is lacking DOF and there's no right thing to do. hopefully the inertial motion will resolve this situation very soon.
                pressure_jump_normal_to_triple_junction = 0;
            }
            assert(pressure_jump_normal_to_triple_junction == pressure_jump_normal_to_triple_junction);
            assert(pressure_jump_tangent_to_triple_junction == pressure_jump_tangent_to_triple_junction);
            
            // combined the two pressure jumps to find the pressure Dirichlet BC
            BC_p[i] = (pressure_jump_normal_to_triple_junction + pressure_jump_tangent_to_triple_junction) * dt;
            BC_dpdn[i] = (v.segment<3>(i * 3).dot(n_solid) + influx) * simOptions().rho;
            BC_dpdna[i] = ((Mat3d::Identity() - n_solid * n_solid.transpose()) * v.segment<3>(i * 3)).dot(n_air) * simOptions().rho;
        }

    }
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/BC").stop();
    assert(BC_p == BC_p);
    assert(BC_dpdn == BC_dpdn);
    assert(BC_dpdna == BC_dpdna);
    
    m_dbg_v1 = BC_p;

    // collocation equation:
    //  omega_j p_j = \int p dGdny - \int dpdn G
    //  omega_j p_j = \sum p_i \int theta_i dGdny - \sum dpdn_i \int theta_i G
    //  omega_j p_j = collocA_{ij} p_i - collocB_{ij} dpdn_i
    //  omega_j p_j = colloc_A^T p - collocB^T dpdn
    std::vector<MatXd> A(p2v.size());
    std::vector<VecXd> rhs(p2v.size());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        size_t Npi = p2v[pi].size();
        A[pi] = MatXd::Zero(Npi, Npi);
        rhs[pi] = VecXd::Zero(Npi);
    }
    
    // allocate one VecXd per vertex, so that the OMP loop below doesn't need to access the same variable in different loop iterations which causes a race. these per-vertex rhs are summed into one VecXd (per partition) after the OMP loop.
    std::vector<std::vector<VecXd> > rhs_per_vertex(p2v.size());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        size_t Npi = p2v[pi].size();
        rhs_per_vertex[pi].assign(Npi, VecXd::Zero(Npi));
    }

    // temporary data (code level optimization for the loop below)
    const double oneover4pi = 1 / (4 * M_PI);
    
    const std::vector<Vec3d> quadrature_square =   BoundaryIntegral::quadrature_square();
    
    std::vector<bool> face_is_solids(nf(), false);
    std::vector<Vec3d> face_outward_normals(nf(), Vec3d::Zero());
    std::vector<double> face_areas(nf(), 0);
    for (size_t i = 0; i < nf(); i++)
    {
        face_is_solids[i] = face_is_solid(i);
        face_outward_normals[i] = face_outward_normal(i);
        face_areas[i] = face_area(i);
    }
    
    std::vector<Vec3d> vert_solid_normals(nv(), Vec3d::Zero());
    std::vector<Vec3d> vert_air_normals(nv(), Vec3d::Zero());
    std::vector<double> vert_solid_areas(nv(), 0);
    std::vector<double> vert_air_areas(nv(), 0);
    for (size_t i = 0; i < nv(); i++)
    {
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][j];
            Vec3d n = face_outward_normals[face_i];
            double a = face_areas[face_i];
            if (face_is_solids[face_i])
                vert_solid_normals[i] += a * n,
                vert_solid_areas[i] += a;   // don't divide by 3 here, because we're looking at the support of the vertex basis function theta_i
            else
                vert_air_normals[i] += a * n,
                vert_air_areas[i] += a;     // don't divide by 3 here, because we're looking at the support of the vertex basis function theta_i
        }
        if (vert_solid_areas[i] > 0)
            vert_solid_normals[i].normalize();
        if (vert_air_areas[i] > 0)
            vert_air_normals[i].normalize();
    }
    
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/coefs_and_assembly").start();
#pragma omp parallel default (none) shared (A, rhs, rhs_per_vertex, BC_p, BC_dpdn, BC_dpdna, TPCF, p2v, p2f, v2p, f2p, flattened_partition_vertices, indices_in_partitions, vert_solid_normals, vert_air_normals, vert_solid_areas, vert_air_areas, face_areas, face_outward_normals, face_is_solids, std::cout)
{
#pragma omp master
    {
//        std::cout << "BEM matrix assembly: spawning " << omp_get_num_threads() << " threads." << std::endl;
    }
        
#pragma omp for schedule (dynamic, 10)
    
    for (size_t fpi = 0; fpi < nv(); fpi++)   // dof on vertex i, where the unknown is either p or dpdn depending on where vertex i is sitting
    {
        size_t i = flattened_partition_vertices[fpi];
        size_t pvi = indices_in_partitions[fpi];
        size_t pi = v2p[i];
        size_t Npi = p2v[pi].size();

        double omega = vert_interior_solid_angle(i) / (4 * M_PI);
        
        double mel = 0;
        for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
            mel += edge_length(mesh().m_vertex_to_edge_map[i][j]);
        mel /= mesh().m_vertex_to_edge_map[i].size();

        // one row of the original collocA_solid, collocA_air, collocB_solid and collocB_air matrices, corresponding to DOF i
        // component j in one of these vectors corresponds to collocation point j
        VecXd collocA_solid_rowi = VecXd::Zero(Npi);
        VecXd collocA_air_rowi   = VecXd::Zero(Npi);
        VecXd collocB_solid_rowi = VecXd::Zero(Npi);
        VecXd collocB_air_rowi   = VecXd::Zero(Npi);
        
        for (size_t pvj = 0; pvj < Npi; pvj++)   // the collocation point is x_j
        {
            size_t j = p2v[pi][pvj];
            
            Vec3d x = pos(j);
            
            if ((x - pos(i)).norm() > 5 * mel)
            {
                // for far-field vertices, use one sample point at the vertex to represent the whole integral domain
                Vec3d dx = pos(i) - x;
                double dxn = dx.norm();
                Vec3d dGdy = oneover4pi * dx / (dxn * dxn * dxn);
                double G = -oneover4pi / dxn;
                
                if (vert_solid_areas[i] > 0)
                    collocA_solid_rowi[pvj] = vert_solid_areas[i] * dGdy.dot(vert_solid_normals[i]) / 3,  // the factor 1/3 is the integral of the basis function theta_i over an incident face of unit area
                    collocB_solid_rowi[pvj] = vert_solid_areas[i] * G / 3;
                if (vert_air_areas[i] > 0)
                    collocA_air_rowi[pvj]   = vert_air_areas[i]   * dGdy.dot(vert_air_normals[i]) / 3,    // the factor 1/3 is the integral of the basis function theta_i over an incident face of unit area
                    collocB_air_rowi[pvj]   = vert_air_areas[i]   * G / 3;
            } else
            {
                // for near-field vertices, use Gaussian quadrature (with Duffy transform)
                for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
                {
                    size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];
                    
                    LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
                    Vec3d n_i = face_outward_normals[face_i];
                    double a_i = face_areas[face_i];
                    
                    Vec3d x0 = pos(fi[0]);
                    Vec3d x1 = pos(fi[1]);
                    Vec3d x2 = pos(fi[2]);
                    
                    double Iij = 0;
                    for (size_t qi = 0; qi < quadrature_square.size(); qi++)
                    {
                        Vec2d qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                        double qiw = quadrature_square[qi].z();
                        
                        Vec3d y;
                        double jacobian_i;
                        Vec3d c_i;
//                        BoundaryIntegral::duffyTransform(xis, 0, qik, y, jacobian_i, c_i);
                        
                        Vec2d t((qik.x() + 1) / 2, (qik.x() + 1) / 2 * (qik.y() + 1) / 2);
                        
                        jacobian_i = t.x() * a_i * 0.5;
                        y = x0 * (1 - t.x()) + x1 * (t.x() - t.y()) + x2 * t.y();
                        
                        double theta_i = 1 - t.x();     // the pyramid function for vertex i
                        
                        Vec3d dx = y - x;
                        double dxn = dx.norm();
                        double dGdy = oneover4pi * dx.dot(n_i) / (dxn * dxn * dxn);
                        double G = -oneover4pi / dxn;
                        
                        double Ia = qiw * jacobian_i * theta_i * dGdy;
                        double Ib = qiw * jacobian_i * theta_i * G;
                        
                        if (face_is_solids[face_i])
                            collocA_solid_rowi[pvj] += Ia,
                            collocB_solid_rowi[pvj] += Ib;
                        else
                            collocA_air_rowi[pvj] += Ia,
                            collocB_air_rowi[pvj] += Ib;
                    }
                }
            }
        }
        
            if (!vertex_is_solid(i))
            {
                // vertex i is in the interior of the liquid-air interface
                //  dpdn is the unknown; p is known
                A[pi].col(pvi) = -(collocB_solid_rowi + collocB_air_rowi);
                rhs_per_vertex[pi][pvi] = -(collocA_solid_rowi + collocA_air_rowi) * BC_p[i];
                rhs[pi][pvi] += omega * BC_p[i];
                
            } else if (vertex_is_solid(i) && !vertex_is_on_triple_junction(i))
            {
                // vertex i is in the interior of the liquid-solid interface
                //  p is the unknown; dpdn is known
                A[pi].col(pvi) =  (collocA_solid_rowi + collocA_air_rowi);
                rhs_per_vertex[pi][pvi] = (collocB_solid_rowi + collocB_air_rowi) * BC_dpdn[i];
                A[pi](pvi, pvi) += -omega;
                
            } else
            {
                // vertex i is on triple junction
                if (TPCF)
                {
                    // TPCF: the triple junction has a positional constraint, which makes p unknown but dpdn on the air side is known (whatever's necessary to cancel the triple junction's lateral motion)
                    A[pi].col(pvi) = (collocA_solid_rowi + collocA_air_rowi);
                    rhs_per_vertex[pi][pvi] = collocB_solid_rowi * BC_dpdn[i] + collocB_air_rowi * BC_dpdna[i];
                    A[pi](pvi, pvi) += -omega;
                } else
                {
                    //  dpdn on the air side is the unknown; dpdn on the solid side and p are both known
                    A[pi].col(pvi) = -collocB_air_rowi;
                    rhs_per_vertex[pi][pvi] = -(collocA_solid_rowi + collocA_air_rowi) * BC_p[i] + collocB_solid_rowi * BC_dpdn[i];
                    rhs[pi][pvi] += omega * BC_p[i];
                }
            }
    }
//#pragma omp for
    
}
//#pragma omp parallel
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/coefs_and_assembly").stop();

    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/assembly_rhssum").start();
    for (size_t pi = 0; pi < p2v.size(); pi++)
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
            rhs[pi] += rhs_per_vertex[pi][pvi];
    rhs_per_vertex.clear();
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/assembly_rhssum").stop();
    
    // solve for the unknowns and interpret the result
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/solve").start();
    std::vector<VecXd> solution(p2v.size());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
//        solution[pi] = A[pi].partialPivLu().solve(rhs[pi]);
        Eigen::BiCGSTAB<MatXd> s(A[pi]);
        solution[pi] = s.solve(rhs[pi]);
    }
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/solve").stop();
    assert(solution[0] == solution[0]);
    
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/vel_update").start();
    VecXd p = BC_p;
    VecXd dpdn = BC_dpdn;
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];   // dof on vertex i, where the unknown is either p or dpdn depending on where vertex i is sitting
            
            if (!vertex_is_solid(i))
            {
                // vertex i is in the interior of the liquid-air interface
                //  dpdn is the unknown; p is known
                dpdn[i] = solution[pi][pvi];
                
            } else if (vertex_is_solid(i) && !vertex_is_on_triple_junction(i))
            {
                // vertex i is in the interior of the liquid-solid interface
                //  p is the unknown; dpdn is known
                p[i] = solution[pi][pvi];
                
            } else
            {
                // vertex i is on triple junction
                if (TPCF)
                {
                    // TPCF: for triple junction the dpdn on the air side is in BC_dpdna, not in the solution; the solution for this vertex is going to be a pressure, as it is no longer known
                    dpdn[i] = BC_dpdna[i];
                    p[i] = solution[pi][pvi];
                } else
                {
                    //  dpdn on the air side is the unknown; dpdn on the solid side and p are both known
                    dpdn[i] = solution[pi][pvi];  // this is for the air side; the solid side is stored in BC_dpdn[i]
                }
            }
        }
    }
    assert(p == p);
    assert(dpdn == dpdn);
    
    // compute the velocity update
    VecXd dv = VecXd::Zero(nv() * 3);
    for (size_t i = 0; i < nv(); i++)
    {
        if (!vertex_is_on_triple_junction(i))
        {
            // for smooth regions (away from triple junctions) the averaged normal should work (at least in the limit of mesh refinement)
            Vec3d n_i = vert_outward_normal(i);

            // the normal pressure derivative comes from either the Neumann BC or the solve
            Vec3d dpdn_i = dpdn[i] * n_i;
            
            // the tangential pressure derivatives are computed by finite differencing the p field (averaging the piecewise-constant gradient on each incident face)
            Vec3d dpdt_i = Vec3d::Zero();
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                LosTopos::Vec3st t = mesh().m_tris[mesh().m_vertex_to_triangle_map[i][j]];
                Vec3d x0 = pos(t[0]);
                Vec3d x1 = pos(t[1]);
                Vec3d x2 = pos(t[2]);
                Vec3d n = (x1 - x0).cross(x2 - x0).normalized();
                Vec3d gradp_integral = (p[t[0]] * n.cross(x2 - x1) + p[t[1]] * n.cross(x0 - x2) + p[t[2]] * n.cross(x1 - x0)) / 2;  // grad p times face area
                
                dpdt_i += gradp_integral / 3;   // a third of the face integral goes to vertex i
            }
            dpdt_i /= vert_area(i);
            
            dv.segment<3>(i * 3) = -(dpdn_i + dpdt_i) / simOptions().rho;
            
        } else
        {
            // the triple junction is not smooth
            // use both BC_dpdn[i] (which corresponds to the solid side) and the dpdn value in solution[i] (which corresponds to the air side) to
            //  determine a gradient in the plane normal to the triple junction.
            // note that the tangential derivative of p also contains information about this gradient. not taking into account this information is dangerous
            //  because n_solid and n_air could be near colinear (e.g. in the case of superhydrophobia surfaces) and the solve using those two alone can be
            //  ill conditioned. solving a normal equation that uses both dpdn on both sides as well as dpdt instead is a better approach.
            Vec3d n_solid = Vec3d::Zero();
            Vec3d n_air = Vec3d::Zero();
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                size_t fj = mesh().m_vertex_to_triangle_map[i][j];
                if (face_is_solid(fj))  n_solid += face_interior_angle(fj, i) * face_outward_normal(fj);
                else                    n_air   += face_interior_angle(fj, i) * face_outward_normal(fj);
            }
            n_solid.normalize();
            n_air.normalize();
            
            Vec3d triple_junction_tangent = n_solid.cross(n_air).normalized();
            
            Vec3d t_solid = n_solid.cross(triple_junction_tangent).normalized();
            Vec3d t_air   = n_air  .cross(triple_junction_tangent).normalized();
                        
            // compute the tangential pressure derivatives by finite differencing the p field, for the air and solid sides separately
            Vec3d dpdt_i_solid = Vec3d::Zero();
            Vec3d dpdt_i_air = Vec3d::Zero();
            double vert_area_solid = 0;
            double vert_area_air = 0;
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                size_t fj = mesh().m_vertex_to_triangle_map[i][j];
                LosTopos::Vec3st t = mesh().m_tris[fj];
                Vec3d x0 = pos(t[0]);
                Vec3d x1 = pos(t[1]);
                Vec3d x2 = pos(t[2]);
                Vec3d n = (x1 - x0).cross(x2 - x0).normalized();
                Vec3d gradp_integral = (p[t[0]] * n.cross(x2 - x1) + p[t[1]] * n.cross(x0 - x2) + p[t[2]] * n.cross(x1 - x0)) / 2;  // grad p times face area
                double fa = face_area(fj);
                
                if (face_is_solid(fj))
                {
                    dpdt_i_solid += gradp_integral / 3;   // a third of the face integral goes to vertex i
                    vert_area_solid += fa;
                } else
                {
                    dpdt_i_air   += gradp_integral / 3;   // a third of the face integral goes to vertex i
                    vert_area_air += fa;
                }

            }
            double dpdt_triple_junction = (dpdt_i_solid + dpdt_i_air).dot(triple_junction_tangent) / vert_area(i);
            dpdt_i_solid /= vert_area_solid;
            dpdt_i_air   /= vert_area_air;
            
            
            // solve a normal equation that includes both dpdn and dpdt information
            MatXd gradp_proj = MatXd::Zero(5, 3);
            VecXd gradp_proj_rhs = VecXd::Zero(5);
            gradp_proj.row(0) = n_solid.transpose();        gradp_proj_rhs[0] = BC_dpdn[i];                 // dpdn on the solid side
            gradp_proj.row(1) = n_air.  transpose();        gradp_proj_rhs[1] = dpdn[i];                    // dpdn on the air side
            gradp_proj.row(2) = t_solid.transpose();        gradp_proj_rhs[2] = dpdt_i_solid.dot(t_solid);  // dpdt on the solid side
            gradp_proj.row(3) = t_air.  transpose();        gradp_proj_rhs[3] = dpdt_i_air  .dot(t_air);    // dpdn on the air side
            gradp_proj.row(4) = triple_junction_tangent;    gradp_proj_rhs[4] = dpdt_triple_junction;       // dpdt along the triple junction
            
            Vec3d gradp = (gradp_proj.transpose() * gradp_proj).partialPivLu().solve(gradp_proj.transpose() * gradp_proj_rhs);
            
            // project the pressure gradient to remove the solid normal component, because that component should be the prescribed BC_dpdn[i]
            gradp = (Mat3d::Identity() - n_solid * n_solid.transpose()) * gradp + n_solid * BC_dpdn[i];
            
            dv.segment<3>(i * 3) = -gradp / simOptions().rho;
        }
    }
    assert(dv == dv);
    
    // TPCF: set the tangential velocity of all solid vertices to the negative of current velocity. this is equivalent to a solid tangential force applied explicitly.
    if (TPCF)
    {
        for (size_t i = 0; i < nf(); i++)
        {
            if (face_is_solid(i))
            {
                Vec3d nsolid = face_outward_normal(i);
                LosTopos::Vec3st t = mesh().m_tris[i];
                for (int j = 0; j < 3; j++)
                {
                    Mat3d P = Mat3d::Identity() - nsolid * nsolid.transpose();
                    dv.segment<3>(t[j] * 3) = -P * v.segment<3>(t[j] * 3) + dv.segment<3>(t[j] * 3).dot(nsolid) * nsolid;    // set tangential velocity to cancel the original velocity's tangential component, without changing the normal component
                }
            }
        }
    }

    // TDMCV: manual (pseudo-)momentum conservation for tiny free droplets where the numerics can't hope to be accurate.
    if (TDMC)
    {
        for (size_t pi = 0; pi < p2v.size(); pi++)
        {
            double meanedgelen = 0;
            int pine = 0;
            for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
            {
                for (size_t j = 0; j < mesh().m_vertex_to_edge_map[p2v[pi][pvi]].size(); j++)
                    meanedgelen += edge_length(mesh().m_vertex_to_edge_map[p2v[pi][pvi]][j]);
                pine += mesh().m_vertex_to_edge_map[p2v[pi][pvi]].size();
            }
            meanedgelen /= pine;
            
            double pivolume = partition_volume(p2f[pi]);
            if (pivolume > pow(meanedgelen * 2, 3) * M_PI * 4 / 3)  // this is not a tiny droplet
                continue;
            
            std::cout << "TDMCV: processing tiny droplet (" << p2v[pi].size() << " vertices, " << p2f[pi].size() << " faces)." << std::endl;
            
            bool solid = false;
            for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
                if (vertex_is_solid(p2v[pi][pvi]))
                    solid = true;
            if (solid)      // this is not a free droplet
                continue;
            
            // remove the global rigid translation mode from dv for this droplet
            Vec3d dvmean = Vec3d(0, 0, 0);
            for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
                dvmean += dv.segment<3>(p2v[pi][pvi] * 3);
            dvmean /= p2v[pi].size();
            for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
                dv.segment<3>(p2v[pi][pvi] * 3) -= dvmean;
        }
    }
    assert(dv == dv);
    
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/BEM/vel_update").stop();

    // update the velocity field
    VecXd newv = v + dv;
    assert(v == v);
    assert(newv == newv);

    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v,    "pressure solve: original"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(dv,   "pressure solve: dv"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(newv, "pressure solve: new"));

    v = newv;
    
    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);

}


















