//
//  BPS3DHelmholtzDecomposition.cpp
//  Droplet3D
//
//  Created by Fang Da on 10/8/15.
//
//

//#include <omp.h>
#include "BPS3D.h"
#include "BoundaryIntegral.h"
#include "FMMTLWrapper.h"


void BPS3D::step_HelmholtzDecomposition(double dt, VecXd & v, const Partitioning & partitioning)
{
    const std::vector<std::vector<size_t> > & p2v = partitioning.p2v;
    const std::vector<std::vector<size_t> > & p2f = partitioning.p2f;
    const std::vector<int>                  & v2p = partitioning.v2p;
    const std::vector<int>                  & f2p = partitioning.f2p;
    const std::vector<size_t> & flattened_partition_vertices = partitioning.flattened_partition_vertices;
    const std::vector<size_t> & indices_in_partitions        = partitioning.indices_in_partitions;

    // remove the global (rigid translation) component, improving accuracy by only reconstructing a smaller component
    std::vector<Vec3d> vglobal(p2v.size(), Vec3d::Zero());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        double areasum = 0;
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];
            vglobal[pi] += v.segment<3>(i * 3) * vert_area(i);
            areasum += vert_area(i);
        }
        vglobal[pi] /= areasum;
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];
            v.segment<3>(i * 3) -= vglobal[pi];
        }
    }
    
    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);
    
    

    // cache geometry some data
    std::vector<double> face_areas(nf());
    std::vector<Vec3d> face_outward_normals(nf());
    std::vector<Vec3d> face_centers(nf());
    std::vector<Vec3d> face_center_velocities(nf());
    for (size_t i = 0; i < nf(); i++)
    {
        LosTopos::Vec3st t = mesh().m_tris[i];
        
        face_areas[i] = face_area(i);
        face_outward_normals[i] = face_outward_normal(i);
        face_centers[i] = (pos(t[0]) + pos(t[1]) + pos(t[2])) / 3;
        face_center_velocities[i] = (v.segment<3>(t[0] * 3) + v.segment<3>(t[1] * 3) + v.segment<3>(t[2] * 3)) / 3;
    }
    
    Vec3d y;
    double jacobian_j;
    const double oneover4pi = 1 / (4 * M_PI);
    
    // compute the gradient field component
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/dPhidn").start();
    std::vector<Vec3d> dPhidn(nf(), Vec3d::Zero()); // one integral for each face as a neighbor of each vertex
    
    const std::vector<Vec2d> quadrature_line =     BoundaryIntegral::quadrature_line();
    const std::vector<Vec3d> quadrature_square =   BoundaryIntegral::quadrature_square();
    const std::vector<Vec3d> quadrature_triangle = BoundaryIntegral::quadrature_triangle();
    const Vec3d qt0 = quadrature_triangle[0];
    
#pragma omp parallel default (none) shared (dPhidn, face_outward_normals, face_centers, face_areas, face_center_velocities, p2v, p2f, v2p, f2p, flattened_partition_vertices, indices_in_partitions, std::cout) private (y, jacobian_j)
{
#pragma omp master
    {
//        std::cout << "HD dPhidn loop: spawning " << omp_get_num_threads() << " threads." << std::endl;
    }
        
#pragma omp for schedule (dynamic, 10)
    for (size_t fpi = 0; fpi < nv(); fpi++)   // degree of freedom on vertex i
    {
        size_t i = flattened_partition_vertices[fpi];
        size_t pvi = indices_in_partitions[fpi];
        size_t pi = v2p[i];
        
            for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
            {
                size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];
                
                LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
                Mat3d xis = getVertexPositions(fi);
                Vec3d n_i = face_outward_normals[face_i];
                
                double I = 0;
                for (size_t qi = 0; qi < quadrature_square.size(); qi++)
                {
                    Vec2d qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                    double qiw = quadrature_square[qi].z();
                    
                    Vec3d x;
                    double jacobian_i;
                    Vec3d c_i;
                    BoundaryIntegral::duffyTransform(xis, 0, qik, x, jacobian_i, c_i);
                    
                    double theta_i = c_i[0];    // the pyramid function for vertex i
                    
                    double I_dSLP = 0;
                    for (size_t pfj = 0; pfj < p2f[pi].size(); pfj++)   // inner integral: SLP normal derivative
                    {
                        size_t j = p2f[pi][pfj];
                        
                        if (j == face_i)
                            continue;   // the inner integral is the normal derivative of SLP, which receives no contribution from the face where the evaluation point x is located in

                        Vec3d n_j = face_outward_normals[j];
                        
                        y = face_centers[j];
                        Vec3d dx = x - y;
                        double dxn = dx.norm();
                        double dGdx = dx.dot(n_i) / (dxn * dxn * dxn) * oneover4pi;

                        I_dSLP += face_areas[j] * face_center_velocities[j].dot(n_j) * dGdx;
                    }
                    
                    I += qiw * jacobian_i * theta_i * I_dSLP;
                }
                
                LosTopos::Vec2st dummy;
                dPhidn[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
            }
    }
//#pragma omp for
        
}
//#pragma omp parallel
    
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/dPhidn").stop();

    
    
    // compute the curl field component
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/curl_A").start();
    std::vector<Vec3d> curl_A(nf(), Vec3d::Zero());  // one integral for each face as a neighbor of each vertex
    
#pragma omp parallel default (none) shared (curl_A, face_outward_normals, face_centers, face_areas, face_center_velocities, p2v, p2f, v2p, f2p, flattened_partition_vertices, indices_in_partitions, std::cout) private (y, jacobian_j)
{
#pragma omp master
    {
//        std::cout << "HD curl_A loop: spawning " << omp_get_num_threads() << " threads." << std::endl;
    }
        
#pragma omp for schedule (dynamic, 10)
    // TODO: merge this loop with the one for dPhidn
    for (size_t fpi = 0; fpi < nv(); fpi++)   // degree of freedom on vertex i
    {
        size_t i = flattened_partition_vertices[fpi];
        size_t pvi = indices_in_partitions[fpi];
        size_t pi = v2p[i];
        
            // integrate over the faces
            for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
            {
                size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];
                
                LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
                Mat3d xis = getVertexPositions(fi);
                Mat3d vis = getVertexVelocities(fi);
                Vec3d n_i = face_outward_normal(face_i);
                double area_fi = face_area(face_i);
                
                double I = 0;
                for (size_t qi = 0; qi < quadrature_square.size(); qi++)
                {
                    Vec2d qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                    double qiw = quadrature_square[qi].z();
                    
                    Vec3d x;
                    double jacobian_i;
                    Vec3d c_i;
                    BoundaryIntegral::duffyTransform(xis, 0, qik, x, jacobian_i, c_i);

                    double theta_i = c_i[0];
                    Vec3d grad_theta_i = (xis.col(1) - xis.col(0)).cross(xis.col(2) - xis.col(0)).normalized().cross(xis.col(2) - xis.col(1)) / (area_fi * 2);
                    
                    Vec3d v_x = vis * c_i;  // find the velocity at x by interpolation

                    Vec3d I_SLP = Vec3d::Zero();
                    for (size_t pfj = 0; pfj < p2f[pi].size(); pfj++)   // inner integral: SLP
                    {
                        size_t j = p2f[pi][pfj];
                        
                        if (j == face_i)
                        {
                            // for the coincident case (j == face_i), subdivide the face at position x to create three faces, each integrated by quadrature individually
                            // otherwise, treat the entire face as a whole
                            for (int k = 0; k < 3; k++)
                            {
                                LosTopos::Vec3st fj = mesh().m_tris[j];
                                Mat3d xjs = getVertexPositions(fj);
                                Mat3d vjs = getVertexVelocities(fj);
                                Vec3d n_j = face_outward_normal(j);

                                xjs.col(k) = x;     // override one of the vertices because this triangle is one of the three subdivision triangles of face j
                                vjs.col(k) = v_x;
                                
                                for (size_t qj = 0; qj < quadrature_square.size(); qj++)
                                {
                                    Vec2d qjk = quadrature_square[qj].segment<2>(0);   // coordinates in the square ref domain
                                    double qjw = quadrature_square[qj].z();
                                    
                                    Vec3d y;
                                    double jacobian_j;
                                    Vec3d c_j;
                                    BoundaryIntegral::duffyTransform(xjs, k, qjk, y, jacobian_j, c_j);   // use duffy transform only for the coincident case
                                    
                                    Vec3d v_y = vjs * c_j;  // find the velocity at y by interpolation
                                    
                                    I_SLP += qjw * jacobian_j * n_j.cross(v_y) * BoundaryIntegral::G(x, y);
                                }
                            }
                        } else
                        {
                            Vec3d n_j = face_outward_normals[j];
                            
                            Vec3d y = face_centers[j];
                            Vec3d dx = x - y;
                            double dxn = dx.norm();
                            double G = -oneover4pi / dxn;
                            
                            I_SLP += face_areas[j] * n_j.cross(face_center_velocities[j]) * G;
                        }
                    }
                    
                    I += -qiw * jacobian_i * n_i.dot(grad_theta_i.cross(I_SLP));
                }
                
                LosTopos::Vec2st dummy;
                curl_A[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
            }

    }
//#pragma omp for

}
//#pragma omp parallel
    
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/curl_A").stop();

    
    
    // the jump term (i.e. half the SLP charge in the scalar potential term)
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/charge").start();
    std::vector<Vec3d> charge(nf(), Vec3d::Zero());  // one integral for each face as a neighbor of each vertex
    
    // TODO: merge this loop with the ones for dPhidn and curl_A
    for (size_t i = 0; i < nv(); i++)
    {
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
        {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];
            
            LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
            Mat3d xis = getVertexPositions(fi);
            Mat3d vis = getVertexVelocities(fi);
            Vec3d n_i = face_outward_normal(face_i);
            double area_fi = face_area(face_i);
            
            double I = -(vis.col(0) / 2 + vis.col(1) / 4 + vis.col(2) / 4).dot(n_i) * (area_fi / 3) / 2;   // jump is half the charge
            
            LosTopos::Vec2st dummy;
            charge[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
        }
    }
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/charge").stop();

    
    
    // composition
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/composition").start();
    VecXd v_dPhidn = VecXd::Zero(nv() * 3);
    VecXd v_charge = VecXd::Zero(nv() * 3);
    VecXd v_curl_A = VecXd::Zero(nv() * 3);
    VecXd v_n = VecXd::Zero(nv() * 3);
    VecXd v_t = VecXd::Zero(nv() * 3);
    VecXd newv = v;
    for (size_t i = 0; i < nv(); i++)
    {
        double area_vi = vert_area(i);
        
        Vec3d n_basis = Vec3d::Zero();
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
        {
            size_t fi = mesh().m_vertex_to_triangle_map[i][ii];
            Vec3d fn = face_outward_normal(fi);
            double fa = face_area(fi);
            n_basis += fn * fa / 3; // n_basis is the volume gradient, or equivalently, the area-weighted average of face normals
        }
        // note that n_basis is not normalized here, retaining its magnitude as volume gradient

        double n_component = 0;             // the normal component (scalar) integrated over the vertex neighborhood (weighted by theta_i), interpreted as volume change rate
        Vec3d v_integral = Vec3d::Zero();   // old velocity (full vector) integrated over the vertex neighborhood (weighted by theta_i)
        double nc_dPhidn = 0;
        double nc_charge = 0;
        double nc_curl_A = 0;
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
        {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];

            LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
            Mat3d xis = getVertexPositions(fi);
            Mat3d vis = getVertexVelocities(fi);
            Vec3d n_i = face_outward_normal(face_i);
            double area_fi = face_area(face_i);
            
            LosTopos::Vec2st dummy;
            size_t vi = mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy);
            
//            n_component += (curl_A[face_i][vi] - (dPhidn[face_i][vi] + charge[face_i][vi]));   // the normal component (integral), which is also the volume change rate
            nc_dPhidn += -dPhidn[face_i][vi];
            nc_charge += -charge[face_i][vi];
            nc_curl_A +=  curl_A[face_i][vi];
            
            v_integral += (vis.col(0) / 2 + vis.col(1) / 4 + vis.col(2) / 4) * (area_fi / 3);   // the integral of old velocity times theta_i
        }
        
        n_component = nc_dPhidn + nc_charge + nc_curl_A;
        
        Mat3d P = Mat3d::Identity() - n_basis * n_basis.transpose() / n_basis.squaredNorm();    // tangent projector
        
        Vec3d vt = P * v_integral / area_vi;
        Vec3d vn = n_component * n_basis / n_basis.squaredNorm();
        
        v_dPhidn.segment<3>(i * 3) = nc_dPhidn * n_basis / n_basis.squaredNorm();
        v_charge.segment<3>(i * 3) = nc_charge * n_basis / n_basis.squaredNorm();
        v_curl_A.segment<3>(i * 3) = nc_curl_A * n_basis / n_basis.squaredNorm();
        v_n.segment<3>(i * 3) = vn;
        v_t.segment<3>(i * 3) = vt;

        newv.segment<3>(i * 3) = vt + vn;
    }
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/composition").stop();

    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_dPhidn,  "HD: dPhidn"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_charge,  "HD: charge"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_curl_A,  "HD: curl A"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_n,       "HD: normal"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_t,       "HD: tangential"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v,         "HD: original"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(newv,      "HD: new"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(newv - v,  "HD: change"));
    
    v = newv;
    
    
    
    // added back the previously removed global (rigid translation) component
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];
            v.segment<3>(i * 3) += vglobal[pi];
        }
    }
    
    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);


}

void BPS3D::step_HD_FMM(double dt, VecXd & v, const BPS3D::Partitioning & partitioning)
{
    const std::vector<std::vector<size_t> > & p2v = partitioning.p2v;
    const std::vector<std::vector<size_t> > & p2f = partitioning.p2f;
    const std::vector<int>                  & v2p = partitioning.v2p;
    const std::vector<int>                  & f2p = partitioning.f2p;
    const std::vector<size_t> & flattened_partition_vertices = partitioning.flattened_partition_vertices;
    const std::vector<size_t> & indices_in_partitions        = partitioning.indices_in_partitions;
    
    // remove the global (rigid translation) component, improving accuracy by only reconstructing a smaller component
    std::vector<Vec3d> vglobal(p2v.size(), Vec3d::Zero());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        double areasum = 0;
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];
            vglobal[pi] += v.segment<3>(i * 3) * vert_area(i);
            areasum += vert_area(i);
        }
        vglobal[pi] /= areasum;
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];
            v.segment<3>(i * 3) -= vglobal[pi];
        }
    }
    
    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);
    
    
    
    // cache geometry some data
    std::vector<double> face_areas(nf());
    std::vector<Vec3d> face_outward_normals(nf());
    std::vector<Vec3d> face_centers(nf());
    std::vector<Vec3d> face_center_velocities(nf());
    for (size_t i = 0; i < nf(); i++)
    {
        LosTopos::Vec3st t = mesh().m_tris[i];
        
        face_areas[i] = face_area(i);
        face_outward_normals[i] = face_outward_normal(i);
        face_centers[i] = (pos(t[0]) + pos(t[1]) + pos(t[2])) / 3;
        face_center_velocities[i] = (v.segment<3>(t[0] * 3) + v.segment<3>(t[1] * 3) + v.segment<3>(t[2] * 3)) / 3;
    }
    
    Vec3d y;
    double jacobian_j;
    const double oneover4pi = 1 / (4 * M_PI);
    
    const std::vector<Vec2d> quadrature_line =     BoundaryIntegral::quadrature_line();
    const std::vector<Vec3d> quadrature_square =   BoundaryIntegral::quadrature_square();
    const std::vector<Vec3d> quadrature_triangle = BoundaryIntegral::quadrature_triangle();
    const Vec3d qt0 = quadrature_triangle[0];
    
    std::vector<size_t> v2pi(nv(), 0);  // map from a vertex in mesh to index of this vertex in its partition p2v[pi] (note that this is not indices_in_partition; that is ordered for flattened_partition_vertices, not the ordering of vertices in LosTopos)
    for (size_t pi = 0; pi < p2v.size(); pi++)
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
            v2pi[p2v[pi][pvi]] = pvi;
    std::vector<size_t> f2pi(nf(), 0);  // map from a face in mesh to index of this face in its partition p2f[pi]
    for (size_t pi = 0; pi < p2f.size(); pi++)
        for (size_t pfi = 0; pfi < p2f[pi].size(); pfi++)
            f2pi[p2f[pi][pfi]] = pfi;

    

    std::vector<Vec3d> dPhidn(nf(), Vec3d::Zero()); // one integral for each face as a neighbor of each vertex
    std::vector<Vec3d> curl_A(nf(), Vec3d::Zero()); // one integral for each face as a neighbor of each vertex
    
    // process one partition at a time
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        // FMMTL setup
        CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/fmm_setup").start();
        
        std::vector<Vec3d> sources;
        sources.reserve(p2f[pi].size());    // sources are the face centers
        for (size_t pfi = 0; pfi < p2f[pi].size(); pfi++)
        {
            LosTopos::Vec3st f = mesh().m_tris[p2f[pi][pfi]];
            sources.push_back((pos(f[0]) + pos(f[1]) + pos(f[2])) / 3);
        }
        
        std::vector<Vec3d> targets;
        targets.reserve(p2f[pi].size() * 3 * quadrature_square.size());    // targets are the duffy 2x2 quadrature points for each of the three vertices respectively
        std::vector<double> theta;
        theta.reserve(targets.size());
        std::vector<Vec3d> grad_theta;
        grad_theta.reserve(targets.size());
        std::vector<Vec3d> target_v;
        target_v.reserve(targets.size());
        for (size_t pfi = 0; pfi < p2f[pi].size(); pfi++)
        {
            size_t face_i = p2f[pi][pfi];

            for (int k = 0; k < 3; k++)
            {
                LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], mesh().m_tris[face_i][k]);
                Mat3d xis = getVertexPositions(fi);
                Mat3d vis = getVertexVelocities(fi);
                Vec3d n_i = face_outward_normals[face_i];
                double area_fi = face_areas[face_i];
                
                for (size_t qi = 0; qi < quadrature_square.size(); qi++)
                {
                    Vec2d qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                    double qiw = quadrature_square[qi].z();
                    
                    Vec3d x;
                    double jacobian_i;
                    Vec3d c_i;
                    BoundaryIntegral::duffyTransform(xis, 0, qik, x, jacobian_i, c_i);

                    targets.push_back(x);
                    
                    theta.push_back(c_i[0]);
                    grad_theta.push_back((xis.col(1) - xis.col(0)).cross(xis.col(2) - xis.col(0)).normalized().cross(xis.col(2) - xis.col(1)) / (area_fi * 2));
                    target_v.push_back(vis * c_i);
                }
            }
        }
        
        std::cout << "Creating FMMPlan" << std::endl;
        FMMPlanP2P fmm(sources, targets);
        
        CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/fmm_setup").stop();

    
    
        // compute the dPhidn part
        CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/dPhidn").start();

        // inner integral
        VecXd dphidn_charges = VecXd::Zero(p2f[pi].size());
        for (size_t pfj = 0; pfj < p2f[pi].size(); pfj++)
        {
            size_t j = p2f[pi][pfj];
            Vec3d nj = face_outward_normals[j];
            Vec3d uj = face_center_velocities[j];
            double aj = face_areas[j];
            dphidn_charges[pfj] = nj.dot(uj) * aj;
        }
        
        VecXd dphidn_result = fmm.evaluateGradSLP(dphidn_charges);  // result contains grad Phi_Gamma (3 components consecutively) at each of the 4 quadrature points for each vertex in each face
        
        // outer integral
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];
            
            for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
            {
                size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];
                LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
                Mat3d xis = getVertexPositions(fi);
                Vec3d ni = face_outward_normals[face_i];
                double ai = face_areas[face_i];
                
                double I = 0;
                for (size_t qi = 0; qi < quadrature_square.size(); qi++)
                {
                    Vec2d qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                    double qiw = quadrature_square[qi].z();
                    
                    Vec3d x;
                    double jacobian_i;
                    Vec3d c_i;
                    BoundaryIntegral::duffyTransform(xis, 0, qik, x, jacobian_i, c_i);
                    
                    LosTopos::Vec2st dummy;
                    size_t dphidn_result_idx = f2pi[face_i] * 3 * quadrature_square.size() + mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy) * quadrature_square.size() + qi;
//                    assert(targets[dphidn_result_idx] == x);

                    Vec3d gradPhi = dphidn_result.segment<3>(dphidn_result_idx * 3);
                    double theta_i = theta[dphidn_result_idx];

                    I += qiw * jacobian_i * theta_i * ni.dot(gradPhi);
                }

                LosTopos::Vec2st dummy;
                dPhidn[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
            }
        }
        
        CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/dPhidn").stop();
        
        
        
        // compute the curl_A part
        CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/curl_A").start();
        
        // inner integral
        VecXd curl_A_charges_x = VecXd::Zero(p2f[pi].size());
        VecXd curl_A_charges_y = VecXd::Zero(p2f[pi].size());
        VecXd curl_A_charges_z = VecXd::Zero(p2f[pi].size());
        for (size_t pfj = 0; pfj < p2f[pi].size(); pfj++)
        {
            size_t j = p2f[pi][pfj];
            Vec3d nj = face_outward_normals[j];
            Vec3d uj = face_center_velocities[j];
            double aj = face_areas[j];
            Vec3d ncrossu = nj.cross(uj);
            curl_A_charges_x[pfj] = ncrossu.x() * aj;
            curl_A_charges_y[pfj] = ncrossu.y() * aj;
            curl_A_charges_z[pfj] = ncrossu.z() * aj;
        }

        VecXd curl_A_result_x = fmm.evaluateSLP(curl_A_charges_x);  // result contains A_Gamma at each of the 4 quadrature points for each vertex in each face
        VecXd curl_A_result_y = fmm.evaluateSLP(curl_A_charges_y);
        VecXd curl_A_result_z = fmm.evaluateSLP(curl_A_charges_z);
        
        // manually fix the local neighborhood of every target, using duffy 2x2 quadrature instead
        for (size_t pfi = 0; pfi < p2f[pi].size(); pfi++)
        {
            size_t face_i = p2f[pi][pfi];
            LosTopos::Vec3st fi = mesh().m_tris[face_i];
            Vec3d n_i = face_outward_normals[face_i];
            double a_i = face_areas[face_i];
            
            LosTopos::Vec3st fj = fi;
            Vec3d n_j = n_i;
            
            for (size_t k = 0; k < 3 * quadrature_square.size(); k++)   // iterate through the quadrature points in this face
            {
                Vec3d x = targets[pfi * 3 * quadrature_square.size() + k];
                Vec3d v_x = target_v[pfi * 3 * quadrature_square.size() + k];
                
                // for the coincident case (j == face_i), subdivide the face at position x to create three faces, each integrated by quadrature individually
                // otherwise, treat the entire face as a whole
                Vec3d I_SLP = Vec3d::Zero();
                for (int k = 0; k < 3; k++)
                {
                    Mat3d xjs = getVertexPositions(fj);
                    Mat3d vjs = getVertexVelocities(fj);
                    
                    xjs.col(k) = x;     // override one of the vertices because this triangle is one of the three subdivision triangles of face j
                    vjs.col(k) = v_x;
                    
                    for (size_t qj = 0; qj < quadrature_square.size(); qj++)
                    {
                        Vec2d qjk = quadrature_square[qj].segment<2>(0);   // coordinates in the square ref domain
                        double qjw = quadrature_square[qj].z();
                        
                        Vec3d y;
                        double jacobian_j;
                        Vec3d c_j;
                        BoundaryIntegral::duffyTransform(xjs, k, qjk, y, jacobian_j, c_j);   // use duffy transform only for the coincident case
                        
                        Vec3d v_y = vjs * c_j;  // find the velocity at y by interpolation
                        
                        I_SLP += qjw * jacobian_j * n_j.cross(v_y) * BoundaryIntegral::G(x, y);
                    }
                }

                Vec3d curl_A_result_local_to_add = I_SLP;
                
                Vec3d curl_A_result_local_to_remove = a_i * n_i.cross(face_center_velocities[face_i]) * BoundaryIntegral::G(x, face_centers[face_i]);
                
                Vec3d curl_A_result_local_modification = curl_A_result_local_to_add - curl_A_result_local_to_remove;
                curl_A_result_x[pfi * 3 * quadrature_square.size() + k] += curl_A_result_local_modification.x();
                curl_A_result_y[pfi * 3 * quadrature_square.size() + k] += curl_A_result_local_modification.y();
                curl_A_result_z[pfi * 3 * quadrature_square.size() + k] += curl_A_result_local_modification.z();
            }
        }
        
        // outer integral
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];
            
            for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
            {
                size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];
                LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
                Mat3d xis = getVertexPositions(fi);
                Vec3d ni = face_outward_normals[face_i];
                double ai = face_areas[face_i];
                
                double I = 0;
                for (size_t qi = 0; qi < quadrature_square.size(); qi++)
                {
                    Vec2d qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                    double qiw = quadrature_square[qi].z();
                    
                    Vec3d x;
                    double jacobian_i;
                    Vec3d c_i;
                    BoundaryIntegral::duffyTransform(xis, 0, qik, x, jacobian_i, c_i);
                    
                    LosTopos::Vec2st dummy;
                    size_t curl_A_result_idx = f2pi[face_i] * 3 * quadrature_square.size() + mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy) * quadrature_square.size() + qi;
                    
                    Vec3d A = Vec3d(curl_A_result_x[curl_A_result_idx], curl_A_result_y[curl_A_result_idx], curl_A_result_z[curl_A_result_idx]);
                    double theta_i = theta[curl_A_result_idx];
                    Vec3d grad_theta_i = grad_theta[curl_A_result_idx];
                    
                    I += -qiw * jacobian_i * ni.dot(grad_theta_i.cross(A));
                }
                
                LosTopos::Vec2st dummy;
                curl_A[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
            }
            
        }
        
        CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/curl_A").stop();
    }
    
    
    
    // the jump term (i.e. half the SLP charge in the scalar potential term)
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/charge").start();
    std::vector<Vec3d> charge(nf(), Vec3d::Zero());  // one integral for each face as a neighbor of each vertex
    
    // TODO: merge this loop with the ones for dPhidn and curl_A
    for (size_t i = 0; i < nv(); i++)
    {
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
        {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];
            
            LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
            Mat3d xis = getVertexPositions(fi);
            Mat3d vis = getVertexVelocities(fi);
            Vec3d n_i = face_outward_normal(face_i);
            double area_fi = face_area(face_i);
            
            double I = -(vis.col(0) / 2 + vis.col(1) / 4 + vis.col(2) / 4).dot(n_i) * (area_fi / 3) / 2;   // jump is half the charge
            
            LosTopos::Vec2st dummy;
            charge[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
        }
    }
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/charge").stop();
    
    
    
    // composition
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/composition").start();
    VecXd v_dPhidn = VecXd::Zero(nv() * 3);
    VecXd v_charge = VecXd::Zero(nv() * 3);
    VecXd v_curl_A = VecXd::Zero(nv() * 3);
    VecXd v_n = VecXd::Zero(nv() * 3);
    VecXd v_t = VecXd::Zero(nv() * 3);
    VecXd newv = v;
    for (size_t i = 0; i < nv(); i++)
    {
        double area_vi = vert_area(i);
        
        Vec3d n_basis = Vec3d::Zero();
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
        {
            size_t fi = mesh().m_vertex_to_triangle_map[i][ii];
            Vec3d fn = face_outward_normal(fi);
            double fa = face_area(fi);
            n_basis += fn * fa / 3; // n_basis is the volume gradient, or equivalently, the area-weighted average of face normals
        }
        // note that n_basis is not normalized here, retaining its magnitude as volume gradient
        
        double n_component = 0;             // the normal component (scalar) integrated over the vertex neighborhood (weighted by theta_i), interpreted as volume change rate
        Vec3d v_integral = Vec3d::Zero();   // old velocity (full vector) integrated over the vertex neighborhood (weighted by theta_i)
        double nc_dPhidn = 0;
        double nc_charge = 0;
        double nc_curl_A = 0;
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++)
        {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];
            
            LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
            Mat3d xis = getVertexPositions(fi);
            Mat3d vis = getVertexVelocities(fi);
            Vec3d n_i = face_outward_normal(face_i);
            double area_fi = face_area(face_i);
            
            LosTopos::Vec2st dummy;
            size_t vi = mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy);
            
            nc_dPhidn += -dPhidn[face_i][vi];
            nc_charge += -charge[face_i][vi];
            nc_curl_A +=  curl_A[face_i][vi];
            
            v_integral += (vis.col(0) / 2 + vis.col(1) / 4 + vis.col(2) / 4) * (area_fi / 3);   // the integral of old velocity times theta_i
        }
        
        n_component = nc_dPhidn + nc_charge + nc_curl_A;
        
        Mat3d P = Mat3d::Identity() - n_basis * n_basis.transpose() / n_basis.squaredNorm();    // tangent projector
        
        Vec3d vt = P * v_integral / area_vi;
        Vec3d vn = n_component * n_basis / n_basis.squaredNorm();
        
        v_dPhidn.segment<3>(i * 3) = nc_dPhidn * n_basis / n_basis.squaredNorm();
        v_charge.segment<3>(i * 3) = nc_charge * n_basis / n_basis.squaredNorm();
        v_curl_A.segment<3>(i * 3) = nc_curl_A * n_basis / n_basis.squaredNorm();
        v_n.segment<3>(i * 3) = vn;
        v_t.segment<3>(i * 3) = vt;
        
        newv.segment<3>(i * 3) = vt + vn;
    }
    CSim::TimerMan::timer("Sim.step/BPS.step/explicit/HD/composition").stop();
    
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_dPhidn,  "HD: dPhidn"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_charge,  "HD: charge"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_curl_A,  "HD: curl A"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_n,       "HD: normal"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v_t,       "HD: tangential"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(v,         "HD: original"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(newv,      "HD: new"));
    m_intermediate_v.push_back(std::pair<VecXd, std::string>(newv - v,  "HD: change"));
    
    v = newv;
    
    
    
    // added back the previously removed global (rigid translation) component
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
        {
            size_t i = p2v[pi][pvi];
            v.segment<3>(i * 3) += vglobal[pi];
        }
    }
    
    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);
    

}

