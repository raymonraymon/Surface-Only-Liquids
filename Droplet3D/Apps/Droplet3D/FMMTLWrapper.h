//
//  FMMTLWrapper.h
//  Droplet3D
//
//  Created by Fang Da on 2/25/16.
//
//

#ifndef __Droplet3D__FMMTLWrapper__
#define __Droplet3D__FMMTLWrapper__

#include <stdio.h>
#include "eigenheaders.h"
#include <vector>

class BPS3D;

// This class encapsulates FMM-accelerated boundary integral evaluation operations. Due to the need to implement quadratures, it is aware of the geometry and mesh connectivity. It does not query LosTopos for this info because of the partitioning. LosTopos stores the entire mesh (all partitions) together in one mesh, while each instance of this class can store a single partition.
class FMMPlan
{
public:
    enum VertexType
    {
        VT_SOLID,
        VT_AIR,
        VT_TRIPLE_JUNCTION,
        
        VT_COUNT
    };
    
    typedef Eigen::Matrix<size_t, 3, 1> Vec3st;
    
    struct Mesh
    {
    public:
        // primary data
        std::vector<Vec3d> vertices;
        std::vector<Vec3st> faces;  // it is assumed that these faces face outwards (i.e. if there were labels, label[1] will be 0 (air))
        std::vector<VertexType> vertex_types;
        
        // derived data (for convenience of access)
        std::vector<std::vector<size_t> > v2f;
        std::vector<Vec3d> face_normals;
        std::vector<double> face_areas;
        std::vector<double> vertex_internal_solid_angles;
        
    public:
        bool vertex_is_solid(size_t v) const { return vertex_types[v] == VT_SOLID || vertex_types[v] == VT_TRIPLE_JUNCTION; }
        bool face_is_solid(size_t f) const { Vec3st t = faces[f]; return (vertex_is_solid(t[0]) && vertex_is_solid(t[1]) && vertex_is_solid(t[2])); }
        int index_in_face(size_t v, size_t f) const { Vec3st t = faces[f]; if (t[0] == v) return 0; else if (t[1] == v) return 1; else if (t[2] == v) return 2; else return -1; }
    };
    
public:
    FMMPlan(const std::vector<Vec3d> & vertices, const std::vector<Vec3st> & faces, const std::vector<VertexType> & vertex_types, const std::vector<double> & vertex_internal_solid_angles);
    virtual ~FMMPlan();
    
          Mesh & mesh()       { return m_mesh; }
    const Mesh & mesh() const { return m_mesh; }
    
public:
    // the following charges should be the DOFs directly, i.e. not having the element areas baked into them.
    virtual VecXd evaluateSLP(const VecXd & face_charges);
    virtual VecXd evaluateDLP(const VecXd & face_charges);
    virtual VecXd evaluateGradSLP(const VecXd & face_charges);
    
    // call FMMPlanF2V::evalauteSLP/DLP() to evaluate SLP/DLP from vertex charges (using a single-point quadrature for all faces other than those in 1-ring, and 4-point duffy for 1-ring faces)
    VecXd evaluateSLPwithVertexCharges(const VecXd & charges);
    VecXd evaluateDLPwithVertexCharges(const VecXd & charges);
    
    // use the SLP/DLP evaluations above to compute Green's representation formula (DLP - SLP - jump)
    // u and dudn are the dirichlet and neumann data on the boundary (the u field is the charge for DLP and jump; the dudn field is the charge for SLP)
    // note that u and dudn are not aligned with vertices. they are 3Nf long (Nf = number of faces) and stores the DOF for each vertex in each face separately, allowing discontinuous elements.
    virtual VecXd evaluateGreensRepresentation(const VecXd & u, const VecXd & dudn);
    
    // this operator is used in the iterative solver to operate on unkonwns. the implementation of this method must be aware of the DOF mapping in the assembly of the linear system.
    VecXd operator * (const VecXd & x);
    
protected:
    Mesh m_mesh;
    
    // the kernel matrix
    void * m_A_SLP; // declared as void * because we don't want to have to include FMMTL headers here, due to the S2T_Compressed duplicate symbol issue. A type cast will be needed whenever it's used.
    void * m_A_DLP;
    
};


// a simpler wrapper, oblivious to the mesh and connectivity, for use with HD
class FMMPlanP2P
{
public:
    FMMPlanP2P(const std::vector<Vec3d> & sources, const std::vector<Vec3d> & targets);
    virtual ~FMMPlanP2P();
    
public:
    // these charges should have element areas baked in, because this wrapper is not aware of the mesh
    virtual VecXd evaluateSLP(const VecXd & charges);
//    virtual VecXd evaluateDLP(const VecXd & charges); // due to lack of mesh info, DLP cannot be implemented; but it's not needed by HD anyway.
    virtual VecXd evaluateGradSLP(const VecXd & charges);
    
protected:
    void * m_A_SLP;
//    void * m_A_DLP;

    std::vector<Vec3d> m_sources;
    std::vector<Vec3d> m_targets;
};


#endif /* defined(__Droplet3D__FMMTLWrapper__) */
