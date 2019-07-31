//
//  Sim.h
//
//  Fang Da 2014
//
//

#ifndef __Sim__
#define __Sim__

#include <iostream>
#include <string>
#include "BPS3D.h"
#include "Scenes.h"

class Sim
{
    friend class Scenes;
public:
    Sim(bool verbose);
    ~Sim();
    
    BPS3D * bps() { return m_bps; }
    
public:
    bool loadOptions(const std::string & option_file, const std::vector<std::pair<std::string, std::string> > & option_overrides);
    bool init();
    
    bool headless() const { return m_headless; }
    
public:
    void step(int substep = 0);
    void stepOutput();
    bool isFinished() const { return m_finished; }
    
    double dt() const { return m_dt; }
    double time() const { return m_time; }
    
    bool incrementTime(int inc = 1, bool loadframe = true);
   
public:
    double currentInflux() const;
    
public:
    enum RenderMode
    {
        RM_TRANSPARENT,
        RM_NONMANIFOLD,
        RM_OPAQUE_FLAT_SHADED,
        RM_OPAQUE_SMOOTH_SHADED,
        
        RM_COUNT
    };
    
    static const int SM_VERTEX = 0x01;
    static const int SM_EDGE   = 0x02;
    static const int SM_FACE   = 0x04;
    void render(RenderMode rm, double field_scale, bool prestepgeometry, const Vec2d & mousepos, int selection_mask = SM_VERTEX | SM_EDGE | SM_FACE);
    
    void showPrimitiveInfo();
    
    void selectPrimitive();
    int nearestVertex() const { return m_nearest_vertex; }
    int nearestEdge()   const { return m_nearest_edge;   }
    int nearestFace()   const { return m_nearest_face;   }
    
    bool isVertexSelected(size_t i) const { return (*m_vertex_selected)[i]; }
    bool isEdgeSelected  (size_t i) const { return (*m_edge_selected)[i]; }
    bool isFaceSelected  (size_t i) const { return (*m_face_selected)[i]; }
    std::vector<size_t> selectedVertices() const { std::vector<size_t> s; for (size_t i = 0; i < m_bps->mesh().nv(); i++) if (isVertexSelected(i)) s.push_back(i); return s; }
    std::vector<size_t> selectedEdges()    const { std::vector<size_t> s; for (size_t i = 0; i < m_bps->mesh().ne(); i++) if (isEdgeSelected(i))   s.push_back(i); return s; }
    std::vector<size_t> selectedFaces()    const { std::vector<size_t> s; for (size_t i = 0; i < m_bps->mesh().nt(); i++) if (isFaceSelected(i))   s.push_back(i); return s; }
    
    void showSelected();
    
public:
    struct Stats
    {
        struct Step
        {
            int nv;
            int ne;
            int nf;
            double total_vol;
            
            double sim_time;
            double real_time;
        };
        
        std::vector<Step> steps;
    };
    
    void addStepToStats(double timecost);
    void statsReport();
    
protected:
    bool m_verbose;
    
    bool m_headless;
    std::string m_scene;
    std::string m_output_directory;
    std::string m_load_directory;
    
    BPS3D * m_bps;

    double m_dt;
    double m_time;
    int m_frameid;
    bool m_finished;

    int m_nearest_vertex;
    int m_nearest_edge;
    int m_nearest_face;
    
    LosTopos::NonDestructiveTriMesh::VertexData<char> * m_vertex_selected;
    LosTopos::NonDestructiveTriMesh::EdgeData<char>   * m_edge_selected;
    LosTopos::NonDestructiveTriMesh::FaceData<char>   * m_face_selected;
  
    Stats m_stats;
};

#endif
