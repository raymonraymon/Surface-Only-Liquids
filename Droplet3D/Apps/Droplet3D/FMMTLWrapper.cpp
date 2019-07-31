//
//  FMMTLWrapper.cpp
//  Droplet3D
//
//  Created by Fang Da on 2/25/16.
//
//

#include "FMMTLWrapper.h"
#include "fmmtl/fmmtl/KernelMatrix.hpp"
//#include "fmmtl/kernel/RMSpherical.hpp"
#include "fmmtl/kernel/BiotSpherical.hpp"
#include "fmmtl/kernel/LaplaceSpherical.hpp"
#include "fmmtl/kernel/DLPSpherical.hpp"
#include "fmmtl/fmmtl/Direct.hpp"
#include "fmmtl/fmmtl/util/Clock.hpp"
#include "BPS3D.h"
#include "BoundaryIntegral.h"
#include "SimOptions.h"
#include "Timer.h"

#define VERBOSE (false)

typedef LaplaceSpherical KernelSLP;
typedef DLPSpherical     KernelDLP;

typedef KernelSLP::point_type  point_type;
typedef KernelSLP::source_type source_type;
typedef KernelSLP::target_type target_type;

typedef KernelSLP::charge_type charge_type_SLP;
typedef KernelSLP::result_type result_type_SLP;
typedef KernelDLP::charge_type charge_type_DLP;
typedef KernelDLP::result_type result_type_DLP;

std::vector<source_type> s_sources;
std::vector<target_type> s_targets;

FMMPlan::FMMPlan(const std::vector<Vec3d> & vertices, const std::vector<Vec3st> & faces, const std::vector<VertexType> & vertex_types, const std::vector<double> & vertex_internal_solid_angles)
{
    // Init the FMM Kernel and options
    FMMOptions opts = get_options(0, NULL);
    
    // Init kernel
    KernelSLP K_SLP;
    KernelDLP K_DLP;
    
    // Copy primary data in mesh
    m_mesh.vertices = vertices;
    m_mesh.faces = faces;
    m_mesh.vertex_types = vertex_types;
    
    // Build derived data in mesh
    m_mesh.face_normals.resize(m_mesh.faces.size(), Vec3d::Zero());
    m_mesh.face_areas.resize(m_mesh.faces.size(), 0);
    m_mesh.v2f.resize(m_mesh.vertices.size(), std::vector<size_t>());
    m_mesh.vertex_internal_solid_angles.resize(m_mesh.vertices.size(), 0);
    
    for (size_t i = 0; i < m_mesh.faces.size(); i++)
    {
        Vec3st t = m_mesh.faces[i];
        
        Vec3d x0 = m_mesh.vertices[t[0]];
        Vec3d x1 = m_mesh.vertices[t[1]];
        Vec3d x2 = m_mesh.vertices[t[2]];
        
        Vec3d n = (x1 - x0).cross(x2 - x0);
        double a = n.norm() / 2;
        n.normalize();
        
        m_mesh.face_normals[i] = n;
        m_mesh.face_areas[i] = a;
        
        m_mesh.v2f[t[0]].push_back(i);
        m_mesh.v2f[t[1]].push_back(i);
        m_mesh.v2f[t[2]].push_back(i);
        
//        std::cout << "face = " << i << " t = " << t.transpose() << " n = " << n.transpose() << " a = " << a << " x0 = " << x0.transpose() << " x1 = " << x1.transpose() << " x2 = " << x2.transpose() << std::endl;
    }
    
    assert(vertex_internal_solid_angles.size() == m_mesh.vertices.size());
    for (size_t i = 0; i < m_mesh.vertices.size(); i++)
        m_mesh.vertex_internal_solid_angles[i] = vertex_internal_solid_angles[i];   // don't do solid angle computation here; just ask the constructor caller to feed it. keep the solid angle computation implementation in BPS3D.
    
    // Init points and charges
    std::vector<source_type> sources;
    std::vector<target_type> targets;
    
    targets.reserve(m_mesh.vertices.size());
    for (size_t i = 0; i < m_mesh.vertices.size(); i++)
    {
        Vec3d x = m_mesh.vertices[i];
        x *= 1.0;
        targets.push_back(Vec<3, double>(x[0], x[1], x[2]));
        if (VERBOSE) std::cout << "target " << i << " = " << x.transpose() << std::endl;
    }
    
    sources.reserve(m_mesh.vertices.size());
    for (size_t j = 0; j < m_mesh.faces.size(); j++)
    {
        Vec3st t = m_mesh.faces[j];
        
        Vec3d x0 = m_mesh.vertices[t[0]];
        Vec3d x1 = m_mesh.vertices[t[1]];
        Vec3d x2 = m_mesh.vertices[t[2]];
        
        Vec3d xp = (x0 + x1 + x2) / 3;
        
        sources.push_back(Vec<3, double>(xp[0], xp[1], xp[2]));
        if (VERBOSE) std::cout << "source " << j << " = " << xp.transpose() << std::endl;
    }
    
    // Build the FMM
    fmmtl::kernel_matrix<KernelSLP> * pA_SLP = new fmmtl::kernel_matrix<KernelSLP>();
    *pA_SLP = K_SLP(targets, sources);
    pA_SLP->set_options(opts);
    m_A_SLP = pA_SLP;
    
    fmmtl::kernel_matrix<KernelDLP> * pA_DLP = new fmmtl::kernel_matrix<KernelDLP>();
    *pA_DLP = K_DLP(targets, sources);
    pA_DLP->set_options(opts);
    m_A_DLP = pA_DLP;

    s_sources = sources;
    s_targets = targets;
    

}

FMMPlan::~FMMPlan()
{
    if (m_A_SLP)
        delete ((fmmtl::kernel_matrix<KernelSLP> *)m_A_SLP);
    if (m_A_DLP)
        delete ((fmmtl::kernel_matrix<KernelDLP> *)m_A_DLP);
}

// result(x) = \int_{\Gamma} charges(y) G(x, y) dA
VecXd FMMPlan::evaluateSLP(const VecXd & face_charges)
{
    assert(face_charges.size() == m_mesh.faces.size());
    assert(m_mesh.faces.size() == m_mesh.face_areas.size());
    assert(m_mesh.faces.size() == m_mesh.face_normals.size());
    
    std::vector<charge_type_SLP> chargesvec(face_charges.size());
    for (size_t i = 0; i < face_charges.size(); i++)
        chargesvec[i] = face_charges[i] * m_mesh.face_areas[i];
    
    fmmtl::kernel_matrix<KernelSLP> * pA_SLP = (fmmtl::kernel_matrix<KernelSLP> *)m_A_SLP;
    std::vector<result_type_SLP> resultsvec = (*pA_SLP) * chargesvec;
    
    VecXd results = VecXd::Zero(resultsvec.size());
    for (size_t i = 0; i < resultsvec.size(); i++)
        results[i] = resultsvec[i][0];
    
    return -results / (4 * M_PI);   // account for the sign and the factor of 4pi: FMMTL defines the Laplacian kernel as G = 1/r, while BEM formulations typically use G = -1/(4 pi r)
}

// result(x) = \int_{\Gamma} charges(y) dG/dny(x, y) dA
// calculation is done as
//       \int_{\Gamma} charges(y) dG/dny(x, y) dA
//    =  \int_{\Gamma} charges(y) ny \cdot \nabla_y G(x, y) dA
//    = -\int_{\Gamma} charges(y) ny \cdot \nabla_x G(x, y) dA
//    = -\int_{\Gamma} charges(y) (ny^0 e^0 + ny^1 e^1 + ny^2 e^2) \cdot \nabla_x G(x, y) dA
//    = -e^0 \cdot \int_{\Gamma} (charges(y) ny^0) \nabla_x G(x, y) dA
//      -e^1 \cdot \int_{\Gamma} (charges(y) ny^1) \nabla_x G(x, y) dA
//      -e^2 \cdot \int_{\Gamma} (charges(y) ny^2) \nabla_x G(x, y) dA
//  where the superscripts 0, 1 and 2 refer to the three cartesian coordinates. after this transformation, the three boundary integrals are each a GradSLP integral.
VecXd FMMPlan::evaluateDLP(const VecXd & face_charges)
{
#if false
    
    assert(face_charges.size() == m_mesh.faces.size());
    assert(m_mesh.faces.size() == m_mesh.face_areas.size());
    assert(m_mesh.faces.size() == m_mesh.face_normals.size());
    
    std::vector<charge_type_SLP> chargesvec0(face_charges.size());
    std::vector<charge_type_SLP> chargesvec1(face_charges.size());
    std::vector<charge_type_SLP> chargesvec2(face_charges.size());
    for (size_t i = 0; i < face_charges.size(); i++)
        chargesvec0[i] = face_charges[i] * m_mesh.face_areas[i] * m_mesh.face_normals[i][0],
        chargesvec1[i] = face_charges[i] * m_mesh.face_areas[i] * m_mesh.face_normals[i][1],
        chargesvec2[i] = face_charges[i] * m_mesh.face_areas[i] * m_mesh.face_normals[i][2];
    
    fmmtl::kernel_matrix<KernelSLP> * pA_SLP = (fmmtl::kernel_matrix<KernelSLP> *)m_A_SLP;
    std::vector<result_type_SLP> resultsvec0 = (*pA_SLP) * chargesvec0;
    std::vector<result_type_SLP> resultsvec1 = (*pA_SLP) * chargesvec1;
    std::vector<result_type_SLP> resultsvec2 = (*pA_SLP) * chargesvec2;
    
    VecXd results = VecXd::Zero(resultsvec0.size());
    for (size_t i = 0; i < resultsvec0.size(); i++)
        results[i] = -(resultsvec0[i][1] + resultsvec1[i][2] + resultsvec2[i][3]);  // there's a minus sign in \int \psi dGdny = -\nabla_x \cdot \int \psi n_y G
    
    return -results / (4 * M_PI);   // account for the sign and the factor of 4pi: FMMTL defines the Laplacian kernel as G = 1/r, while BEM formulations typically use G = -1/(4 pi r)
    
#else
    
    assert(face_charges.size() == m_mesh.faces.size());
    assert(m_mesh.faces.size() == m_mesh.face_areas.size());
    assert(m_mesh.faces.size() == m_mesh.face_normals.size());
    
    std::vector<charge_type_DLP> chargesvec(face_charges.size());
    for (size_t i = 0; i < face_charges.size(); i++)
    {
        Vec3d n = m_mesh.face_normals[i];
        chargesvec[i] = face_charges[i] * m_mesh.face_areas[i] * Vec<3, double>(n[0], n[1], n[2]);
    }
    
    fmmtl::kernel_matrix<KernelDLP> * pA_DLP = (fmmtl::kernel_matrix<KernelDLP> *)m_A_DLP;
    std::vector<result_type_DLP> resultsvec = (*pA_DLP) * chargesvec;
    
    VecXd results = VecXd::Zero(resultsvec.size());
    for (size_t i = 0; i < resultsvec.size(); i++)
        results[i] = resultsvec[i];
    
    return -results / (4 * M_PI);   // account for the sign and the factor of 4pi: FMMTL defines the Laplacian kernel as G = 1/r, while BEM formulations typically use G = -1/(4 pi r)
    
#endif
}

// [result(3x), result(3x+1), result(3x+2)]^T = \nabla \int_{\Gamma} charges(y) G(x, y) dA = \int_{\Gamma} charges(y) \nabla_x G(x, y) dA
VecXd FMMPlan::evaluateGradSLP(const VecXd & face_charges)
{
    assert(face_charges.size() == m_mesh.faces.size());
    assert(m_mesh.faces.size() == m_mesh.face_areas.size());
    assert(m_mesh.faces.size() == m_mesh.face_normals.size());
    
    std::vector<charge_type_SLP> chargesvec(face_charges.size());
    for (size_t i = 0; i < face_charges.size(); i++)
        chargesvec[i] = face_charges[i] * m_mesh.face_areas[i];
    
    fmmtl::kernel_matrix<KernelSLP> * pA_SLP = (fmmtl::kernel_matrix<KernelSLP> *)m_A_SLP;
    std::vector<result_type_SLP> resultsvec = (*pA_SLP) * chargesvec;
    
    VecXd results = VecXd::Zero(resultsvec.size() * 3);
    for (size_t i = 0; i < resultsvec.size(); i++)
        results[3 * i + 0] = resultsvec[i][1],
        results[3 * i + 1] = resultsvec[i][2],
        results[3 * i + 2] = resultsvec[i][3];
    
    return -results / (4 * M_PI);   // account for the sign and the factor of 4pi: FMMTL defines the Laplacian kernel as G = 1/r, while BEM formulations typically use G = -1/(4 pi r)
}

VecXd FMMPlan::evaluateSLPwithVertexCharges(const VecXd & vertex_charges)
{
    if (VERBOSE) std::cout << "slpv: v charges = " << vertex_charges.transpose() << std::endl;
    
    assert(vertex_charges.size() == m_mesh.faces.size() * 3);

    VecXd face_charges = VecXd::Zero(m_mesh.faces.size());
    for (size_t i = 0; i < m_mesh.faces.size(); i++)
    {
        Vec3st t = m_mesh.faces[i];
        face_charges[i] = (vertex_charges[i * 3 + 0] + vertex_charges[i * 3 + 1] + vertex_charges[i * 3 + 2]) / 3;  // single-point quadrature
    }
    
    if (VERBOSE) std::cout << "slpv: f charges = " << face_charges.transpose() << std::endl;

    VecXd result = evaluateSLP(face_charges);
    
    if (VERBOSE) std::cout << "slpv: preliminary result = " << result.transpose() << std::endl;
    
    // remove the S2T contributions from these 1-ring faces that are baked in by FMMTL: they'll be accounted for explicitly below
    for (size_t i = 0; i < m_mesh.vertices.size(); i++)
    {
        Vec3d x = m_mesh.vertices[i];
        for (size_t ii = 0; ii < m_mesh.v2f[i].size(); ii++)
        {
            size_t face_i = m_mesh.v2f[i][ii];
            
            Vec3st t = m_mesh.faces[face_i];
            Vec3d y = (m_mesh.vertices[t[0]] + m_mesh.vertices[t[1]] + m_mesh.vertices[t[2]]) / 3;
            
            double dxn = (y - x).norm();
            double G = -1 / (4 * M_PI * dxn);
            
            result[i] -= G * m_mesh.face_areas[face_i] * face_charges[face_i];
            
        }
    }
    
    if (VERBOSE) std::cout << "slpv: result after removal = " << result.transpose() << std::endl;
    
    // add back the contributions from 1-ring faces, using 4-point duffy quadrature
    const std::vector<Vec3d> quadrature_square = BoundaryIntegral::quadrature_square();
    for (size_t i = 0; i < m_mesh.vertices.size(); i++)
    {
        Vec3d x = m_mesh.vertices[i];
        for (size_t ii = 0; ii < m_mesh.v2f[i].size(); ii++)
        {
            size_t face_i = m_mesh.v2f[i][ii];
            
            Vec3st fi = m_mesh.faces[face_i];
            Vec3st fip = Vec3st(0, 1, 2);   // permutation applied to vertices of fi to achieve the shuffle effect of BPS3D::getShuffledTriangle(fi, i)
            assert(fi[0] == i || fi[1] == i || fi[2] == i);
            while (fi[fip[0]] != i)  // shuffle (BPS3D::getShuffledTriangle(fi, i))
            {
                size_t tmp = fip[0];
                fip[0] = fip[1];
                fip[1] = fip[2];
                fip[2] = tmp;
            }
            assert(fi[fip[0]] == i);
            
            Vec3d n_i = m_mesh.face_normals[face_i];
            double a_i = m_mesh.face_areas[face_i];
            
            Vec3d x0 = m_mesh.vertices[fi[fip[0]]];
            Vec3d x1 = m_mesh.vertices[fi[fip[1]]];
            Vec3d x2 = m_mesh.vertices[fi[fip[2]]];
            
            double Iij = 0;
            for (size_t qi = 0; qi < quadrature_square.size(); qi++)
            {
                Vec2d qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                double qiw = quadrature_square[qi].z();
                
                Vec3d y;
                double jacobian_i;
                Vec3d c_i;
//                BoundaryIntegral::duffyTransform(xis, 0, qik, y, jacobian_i, c_i);
                
                Vec2d t((qik.x() + 1) / 2, (qik.x() + 1) / 2 * (qik.y() + 1) / 2);
                
                jacobian_i = t.x() * a_i * 0.5;
//                y = x0 * (1 - t.x()) + x1 * (t.x() - t.y()) + x2 * t.y();
                
                double theta_i = 1 - t.x();     // the pyramid function for vertex i (i.e. fi[0])
                double theta_1 = t.x() - t.y(); //                      for vertex fi[1]
                double theta_2 = t.y();         //                      for vertex fi[2]
                
                y = x0 * theta_i + x1 * theta_1 + x2 * theta_2;
                double charge = vertex_charges[face_i * 3 + fip[0]] * theta_i + vertex_charges[face_i * 3 + fip[1]] * theta_1 + vertex_charges[face_i * 3 + fip[2]] * theta_2;
                
                double dxn = (y - x).norm();
                double G = -1 / (4 * M_PI * dxn);

                result[i] += G * jacobian_i * qiw * charge;
                
            }
        }
    }
    
    if (VERBOSE) std::cout << "slpv result = " << result.transpose() << std::endl;

    return result;
}

VecXd FMMPlan::evaluateDLPwithVertexCharges(const VecXd & vertex_charges)
{
    if (VERBOSE) std::cout << "dlpv: v charges = " << vertex_charges.transpose() << std::endl;

    assert(vertex_charges.size() == m_mesh.faces.size() * 3);
    
    VecXd face_charges = VecXd::Zero(m_mesh.faces.size());
    for (size_t i = 0; i < m_mesh.faces.size(); i++)
    {
        Vec3st t = m_mesh.faces[i];
        face_charges[i] = (vertex_charges[i * 3 + 0] + vertex_charges[i * 3 + 1] + vertex_charges[i * 3 + 2]) / 3;  // single-point quadrature
    }

    if (VERBOSE) std::cout << "dlpv: f charges = " << face_charges.transpose() << std::endl;

    VecXd result = FMMPlan::evaluateDLP(face_charges);
    
    if (VERBOSE) std::cout << "dlpv result = " << result.transpose() << std::endl;
    
    return result;
}

VecXd FMMPlan::evaluateGreensRepresentation(const VecXd & u, const VecXd & dudn)
{
    CSim::TimerMan::timer("FMM_multiply/GreensRep").start();
    
    assert(u.size() == m_mesh.faces.size() * 3);
    assert(dudn.size() == m_mesh.faces.size() * 3);
    
    VecXd slp = evaluateSLPwithVertexCharges(dudn);
    VecXd dlp = evaluateDLPwithVertexCharges(u);
    
    size_t N = m_mesh.vertices.size();
    VecXd jump = VecXd::Zero(N);
    for (size_t i = 0; i < N; i++)
    {
        double omega = m_mesh.vertex_internal_solid_angles[i];
        double ui = 0;
        for (size_t j = 0; j < m_mesh.v2f[i].size(); j++)
            ui += u[m_mesh.v2f[i][j] * 3 + m_mesh.index_in_face(i, m_mesh.v2f[i][j])];
        ui /= m_mesh.v2f[i].size(); // although we take the average here, in practice the u data is always continuous, so faces incident to a vertex should all have the same value for that vertex.
        jump[i] = ui * omega / (4 * M_PI);
    }
    
    if (VERBOSE) std::cout << "slp = " << slp.transpose() << std::endl << "dlp = " << dlp.transpose() << std::endl << "jump = " << jump.transpose() << std::endl << "result = " << (dlp - slp - jump).transpose() << std::endl;
    
    VecXd result = dlp - slp - jump;
    
    CSim::TimerMan::timer("FMM_multiply/GreensRep").stop();
    
    return result;
}

VecXd FMMPlan::operator * (const VecXd & x)
{
    CSim::TimerMan::timer("FMM_multiply").start();
    
    bool TPCF = Options::boolValue("tpcf");
    if (TPCF)
        assert(!"not implemented");
    
    VecXd p    = VecXd::Zero(m_mesh.faces.size() * 3);   // store the charge for each vertex on each face separately (allowing for discontinuous elements; this helps with handling triple junctions)
    VecXd dpdn = VecXd::Zero(m_mesh.faces.size() * 3);

    // collect data in x (each element may store a p or a dpdn depending on the corresponding element type) into p and dpdn
    // data not carried by x will be left zero. for example, the
    for (size_t i = 0; i < m_mesh.faces.size(); i++)
    {
        Vec3st f = m_mesh.faces[i];
        
        if (m_mesh.face_is_solid(i))
        {
            // for all vertices of a solid face: dpdn is known. p is unknown for non-triple-junction vertices; on the triple junction the air-side dpdn is the unknown.
            for (size_t j = 0; j < 3; j++)
                if (m_mesh.vertex_types[f[j]] != VT_TRIPLE_JUNCTION)
                    p[i * 3 + j] = x[f[j]]; // this element in x stores the p
                else
                    ;                           // this element in x stores the air-side dpdn, which should not be reflected in this solid triangle
        } else
        {
            // for all vertices of an air face: p is known, and dpdn is unknown. on the triple junction the unknown is the air-side dpdn to be precise.
            for (size_t j = 0; j < 3; j++)
                dpdn[i * 3 + j] = x[f[j]];
        }
    }
    
//    if (VERBOSE) std::cout << "x = " << x.transpose() << std::endl << "p = " << p.transpose() << std::endl << "dpdn = " << dpdn.transpose() << std::endl;
    
    VecXd result = evaluateGreensRepresentation(p, dpdn);
    
    CSim::TimerMan::timer("FMM_multiply").stop();
    
    return result;
}




FMMPlanP2P::FMMPlanP2P(const std::vector<Vec3d> & source_points, const std::vector<Vec3d> & target_points)
{
    // Init the FMM Kernel and options
    FMMOptions opts = get_options(0, NULL);
    
    // Init kernel
    KernelSLP K_SLP;
    KernelDLP K_DLP;
    
    // Init points and charges
    std::vector<source_type> sources;
    std::vector<target_type> targets;
    
    targets.reserve(target_points.size());
    for (size_t i = 0; i < target_points.size(); i++)
    {
        Vec3d x = target_points[i];
        x *= 1.0;
        targets.push_back(Vec<3, double>(x[0], x[1], x[2]));
        if (VERBOSE) std::cout << "target " << i << " = " << x.transpose() << std::endl;
    }
    
    sources.reserve(source_points.size());
    for (size_t j = 0; j < source_points.size(); j++)
    {
        Vec3d x = source_points[j];
        
        sources.push_back(Vec<3, double>(x[0], x[1], x[2]));
        if (VERBOSE) std::cout << "source " << j << " = " << x.transpose() << std::endl;
    }
    
    // Build the FMM
    fmmtl::kernel_matrix<KernelSLP> * pA_SLP = new fmmtl::kernel_matrix<KernelSLP>();
    *pA_SLP = K_SLP(targets, sources);
    pA_SLP->set_options(opts);
    m_A_SLP = pA_SLP;
    
    m_sources = source_points;
    m_targets = target_points;
}

FMMPlanP2P::~FMMPlanP2P()
{
    if (m_A_SLP)
        delete ((fmmtl::kernel_matrix<KernelSLP> *)m_A_SLP);
}

// result(x) = \int_{\Gamma} charges(y) G(x, y) dA
VecXd FMMPlanP2P::evaluateSLP(const VecXd & charges)
{
    assert(charges.size() == m_sources.size());
    
    std::vector<charge_type_SLP> chargesvec(charges.size());
    for (size_t i = 0; i < charges.size(); i++)
        chargesvec[i] = charges[i];
    
    fmmtl::kernel_matrix<KernelSLP> * pA_SLP = (fmmtl::kernel_matrix<KernelSLP> *)m_A_SLP;
    std::vector<result_type_SLP> resultsvec = (*pA_SLP) * chargesvec;
    
    VecXd results = VecXd::Zero(resultsvec.size());
    for (size_t i = 0; i < resultsvec.size(); i++)
        results[i] = resultsvec[i][0];
    
    return -results / (4 * M_PI);   // account for the sign and the factor of 4pi: FMMTL defines the Laplacian kernel as G = 1/r, while BEM formulations typically use G = -1/(4 pi r)
}


// [result(3x), result(3x+1), result(3x+2)]^T = \nabla \int_{\Gamma} charges(y) G(x, y) dA = \int_{\Gamma} charges(y) \nabla_x G(x, y) dA
VecXd FMMPlanP2P::evaluateGradSLP(const VecXd & charges)
{
    assert(charges.size() == m_sources.size());
    
    std::vector<charge_type_SLP> chargesvec(charges.size());
    for (size_t i = 0; i < charges.size(); i++)
        chargesvec[i] = charges[i];
    
    fmmtl::kernel_matrix<KernelSLP> * pA_SLP = (fmmtl::kernel_matrix<KernelSLP> *)m_A_SLP;
    std::vector<result_type_SLP> resultsvec = (*pA_SLP) * chargesvec;
    
    VecXd results = VecXd::Zero(resultsvec.size() * 3);
    for (size_t i = 0; i < resultsvec.size(); i++)
        results[3 * i + 0] = resultsvec[i][1],
        results[3 * i + 1] = resultsvec[i][2],
        results[3 * i + 2] = resultsvec[i][3];
    
    return -results / (4 * M_PI);   // account for the sign and the factor of 4pi: FMMTL defines the Laplacian kernel as G = 1/r, while BEM formulations typically use G = -1/(4 pi r)
}
