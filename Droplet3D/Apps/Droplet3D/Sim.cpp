//
//  Sim.cpp
//
//  Fang Da 2014
//
//

#include <sstream>
#include <fstream>
#include <iomanip>
#include "Sim.h"
#include "SimOptions.h"
#include "MeshIO.h"
#include <cmath>
#include "eigenheaders.h"
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "YImage.h"
#include <omp.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <GLFW/glfw3.h>

Sim::Sim(bool verbose) :
    m_verbose(verbose),
    m_headless(false),
    m_scene("unspecified"),
    m_output_directory(""),
    m_bps(NULL),
    m_dt(0),
    m_time(0),
    m_frameid(0),
    m_finished(false),
    m_nearest_vertex(-1),
    m_nearest_edge(-1),
    m_nearest_face(-1)
{
    
}

Sim::~Sim()
{

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  General initialization of a simulation
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Sim::loadOptions(const std::string & option_file, const std::vector<std::pair<std::string, std::string> > & option_overrides)
{
    // declare and load the options
    Options::addStringOption ("scene", "sphere");
    Options::addBooleanOption("output", false);
    Options::addBooleanOption("headless", false);
    Options::addStringOption ("load-dir", "");
    Options::addStringOption ("output-dir", "");
    Options::addDoubleOption ("time-step", 0.01);
    Options::addDoubleOption ("simulation-time", 10.0);
    Options::addDoubleOption ("start-time", 0.0);
    Options::addStringOption ("start-rec", "");
    Options::addIntegerOption("num-threads", 0);
    
    Options::addBooleanOption("implicit-integration", false);
    Options::addDoubleOption ("smoothing-coef", 10.0);
    Options::addDoubleOption ("sigma", 1.0);
    Options::addDoubleOption ("sigma-sl", 1.0);
    Options::addDoubleOption ("sigma-sa", 1.0);
    Options::addDoubleOption ("gravity", 0.0);
    Options::addDoubleOption ("rho", 1.0);
    Options::addBooleanOption("fmmtl", false);

    Options::addBooleanOption("tpcf", false);                       // triple junction positional constraint (with post-stabilization that eliminates the tangential velocity on solid vertices)
    Options::addBooleanOption("tdmc", true);                        // remove global translational modes from the pressure solve's velocity update for tiny free droplets (those with a volume smaller than the sphere of a radius twice the average edge length of this droplet, with no solid contact)
    Options::addBooleanOption("nrs", false);                        // not refining solids. set to true to use a large target edge length in solid interior regardless of curvature/velocity (100 times the global max target edge length).
    Options::addBooleanOption("nrtp", false);                       // not refining triple junction. set to true to use the global max target edge length for triple junction vertices curvature/velocity.
    Options::addBooleanOption("csoc", false);                       // concave smoothing overshoot clamp (prevent smoothing from overshooting to create tumors in concave cavities)
    Options::addDoubleOption ("lspmt-ratio", 1.0);                  // ratio between m_merge_proximity_epsilon_for_liquid_sheet_puncture and m_merge_proximity_epsilon in LosTopos
    Options::addIntegerOption("hd-interval", 0) ;                   // how often HD is called (how many frames between two HD calls). 1 (or 0) means calling HD every frame. 
    Options::addDoubleOption ("hd-interval-start", 0);              // the time at which hd-interval starts to take effect. before that time, HD is called every frame (hd-interval = 0)
    Options::addBooleanOption("ibr", true);                         // internal bubble removal. generally it should be turned on since there's no theoretical support for what happens there, but for certain scenes (such as jet bell) it may be desirable to turn it off.

    Options::addDoubleOption ("remove-partitions-below-z", 20);     // remove any closed meshes whose vertices all have z coordinates below this threshold (to help culling far away droplets in scenes like dripping). From 20151215, this is used as a radius to cull droplets far away from the origin, regardless of direction.
    Options::addDoubleOption ("influx", 0.0);                       // influx (normal velocity) through the solid surface
    Options::addDoubleOption ("influx-inc", 0.0);                   // the influx change rate (influx is incremented by this amount times dt after every time step)
    Options::addDoubleOption ("floor-z", -1e+10);                   // floor's z coordinate. vertices lower than this will be considered in contact with the floor.
    Options::addDoubleOption ("init-vel-x", 0.0);                   // initial velocity of droplet for scene tet, cube, floorsplash, transportinject, crownsplash
    Options::addDoubleOption ("init-vel-y", 0.0);
    Options::addDoubleOption ("init-vel-z", -1.0);
    Options::addDoubleOption ("nozzle-angle", 45);                  // nozzle orientation angle (in degrees, not radians; 0 means downward) for scene jet, jetchain
    Options::addDoubleOption ("radius", 1.0);                       // radius of the droplet for scene sphere, collision, dripping, floorsplash, transport (wider end), transportinject, crownsplash, and radius of the nozzle for scene jet and jetchain
    Options::addDoubleOption ("distance", 3.0);                     // distance between nozzle centers in scene jetchain, and position of the injected droplet along the strip in transportinject
    Options::addIntegerOption("nozzle-arrangement", 0);             // for scene jetchain: 0 for left and right (x axis), 1 for front and back (y axis), 2 for down and up (z axis). at nozzle-angle = 0, nozzles will be facing -z, -x, -y direction for nozzle-arrangement = 0, 1, 2 respectively.
    Options::addDoubleOption ("nozzle-radius-ratio", 1.0);          // for scene jetchain: the ratio between the right nozzle and the left nozzle
    Options::addDoubleOption ("offcenter-distance", 0);             // for scene collision: the distance in z direction. using a nonzero offcenter-distance with the default horizontal initial velocity creates an off-center collision.
    Options::addDoubleOption ("gravity-dir-x", 0.0);                // gravity direction (will be normalized before multiplying by "gravity")
    Options::addDoubleOption ("gravity-dir-y", 0.0);
    Options::addDoubleOption ("gravity-dir-z", -1.0);
    
    Options::addBooleanOption("output-png", true);
    Options::addIntegerOption("output-png-every-n-frames", 0);     // 0 means synching with simulation frame rate (equivalent to 1).
    Options::addBooleanOption("output-mesh", true);
    Options::addIntegerOption("output-mesh-every-n-frames", 0);    // 0 means synching with simulation frame rate (equivalent to 1).
    Options::addBooleanOption("output-obj", false);
    Options::addIntegerOption("output-obj-every-n-frames", 0);     // 0 means synching with simulation frame rate (equivalent to 1).
    
    Options::addDoubleOption ("remeshing-resolution", 0);   // 0 means using the initial mean edge length
    Options::addDoubleOption ("remeshing-dr", 100);         // dynamic range of the adaptive remeshing resolution (max target edge length = resolution * dr)
    Options::addIntegerOption("remeshing-iterations", 1);

    Options::addDoubleOption ("lostopos-collision-epsilon-fraction", 1e-4);       // lostopos collision epsilon (fraction of mean edge length)
    Options::addDoubleOption ("lostopos-merge-proximity-epsilon-fraction", 0.1);  // lostopos merge proximity epsilon (fraction of mean edge length)
    Options::addBooleanOption("lostopos-perform-smoothing", false);               // whether or not to perform smoothing
    Options::addDoubleOption ("lostopos-max-volume-change-fraction", 1e-2);       // maximum allowed volume change during a remeshing operation (fraction of mean edge length cubed)
    Options::addDoubleOption ("lostopos-min-triangle-angle", 3.0);                // min triangle angle (in degrees)
    Options::addDoubleOption ("lostopos-max-triangle-angle", 177.0);              // max triangle angle (in degrees)
    Options::addDoubleOption ("lostopos-large-triangle-angle-to-split", 160.0);   // threshold for large angles to be split
    Options::addDoubleOption ("lostopos-min-triangle-area-fraction", 0.02);       // minimum allowed triangle area (fraction of mean edge length squared)
    Options::addBooleanOption("lostopos-t1-transition-enabled", true);            // whether t1 is enabled
    Options::addDoubleOption ("lostopos-t1-pull-apart-distance-fraction", 0.1);   // t1 pull apart distance (fraction of mean edge legnth)
    Options::addBooleanOption("lostopos-smooth-subdivision", false);              // whether to use smooth subdivision during remeshing
    Options::addBooleanOption("lostopos-allow-non-manifold", true);               // whether to allow non-manifold geometry in the mesh
    Options::addBooleanOption("lostopos-allow-topology-changes", true);           // whether to allow topology changes
    Options::addDoubleOption ("lostopos-iar-coef-curvature", 0.05);               // SurfTrack::m_target_edge_length_coef_curvature
    Options::addDoubleOption ("lostopos-iar-coef-velocity", 0.01);                // SurfTrack::m_target_edge_length_coef_velocity
    Options::addDoubleOption ("lostopos-iar-adjacent-ratio", 1.5);                // SurfTrack::m_max_adjacent_target_edge_length_ratio
    
    Options::addIntegerOption("mesh-size-n", 2);
    Options::addIntegerOption("mesh-size-m", 2);
    
    Options::parseOptionFile(option_file, option_overrides, m_verbose);
    
    m_scene = Options::strValue("scene");
    m_headless = Options::boolValue("headless");

    return true;
}

bool Sim::init()
{
    if (Options::boolValue("output"))
    {
        m_output_directory = Options::strValue("output-dir");
        if (m_output_directory == "auto")
        {
            std::stringstream output_dir_ss;
            output_dir_ss << "output_" << ::time(NULL);
            m_output_directory = output_dir_ss.str();
        }
        
        mkdir(m_output_directory.c_str(), 0755);
        std::cout << "Outputing to directory: " << m_output_directory << std::endl;
        
        std::ofstream fo((m_output_directory + "/options.txt").c_str());
        Options::outputOptionValues(fo);
        fo.close();
    }
    
    m_load_directory = Options::strValue("load-dir");
    if (m_scene == "load")  // a replay run: load the first rec file in the load-dir directory
    {
        assert(m_load_directory != "");
        std::cout << "Loading from directory: " << m_load_directory << std::endl;
        
        BPS3D * bps = new BPS3D(this, std::vector<LosTopos::Vec3d>(), std::vector<LosTopos::Vec3st>(), std::vector<LosTopos::Vec2i>(), std::vector<Vec3d>(), std::vector<char>());
        MeshIO::load(*bps, m_load_directory + "/mesh000000.rec");
        
        std::vector<LosTopos::Vec3d> vertices = bps->m_st->pm_positions;
        std::vector<LosTopos::Vec3st> faces = bps->m_st->m_mesh.m_tris;
        std::vector<LosTopos::Vec2i> face_labels = bps->m_st->m_mesh.m_triangle_labels;
        std::vector<Vec3d> vels = bps->vel();
        std::vector<char> solids = bps->solid_labels();
        
        m_bps = new BPS3D(this, vertices, faces, face_labels, vels, solids);
        MeshIO::load(*m_bps, m_load_directory + "/mesh000000.rec");
        
    } else if (Options::strValue("start-rec") != "")    // a continuation run: load the specified rec file as the initial state
    {
        std::string startrec = Options::strValue("start-rec");
        std::cout << "Loading " << startrec << " as the initial state." << std::endl;
        
        BPS3D * bps = new BPS3D(this, std::vector<LosTopos::Vec3d>(), std::vector<LosTopos::Vec3st>(), std::vector<LosTopos::Vec2i>(), std::vector<Vec3d>(), std::vector<char>());
        MeshIO::load(*bps, startrec);
        
        std::vector<LosTopos::Vec3d> vertices = bps->m_st->pm_positions;
        std::vector<LosTopos::Vec3st> faces = bps->m_st->m_mesh.m_tris;
        std::vector<LosTopos::Vec2i> face_labels = bps->m_st->m_mesh.m_triangle_labels;
        std::vector<Vec3d> vels = bps->vel();
        std::vector<char> solids = bps->solid_labels();
        
        m_bps = new BPS3D(this, vertices, faces, face_labels, vels, solids);
        MeshIO::load(*m_bps, startrec);
        
    } else      // a normal run: call scene() to construct the initial state
    {
        std::vector<LosTopos::Vec3d> vertices;
        std::vector<LosTopos::Vec3st> faces;
        std::vector<LosTopos::Vec2i> face_labels;
        std::vector<Vec3d> vels;
        std::vector<char> solids;
        
        m_bps = Scenes::scene(this, m_scene, vertices, faces, face_labels, vels, solids);
        
        std::cout << "nv = " << vertices.size() << " nf = " << faces.size() << std::endl;
    }
    
    m_vertex_selected = new LosTopos::NonDestructiveTriMesh::VertexData<char>(&(m_bps->mesh()));
    m_edge_selected   = new LosTopos::NonDestructiveTriMesh::EdgeData<char>  (&(m_bps->mesh()));
    m_face_selected   = new LosTopos::NonDestructiveTriMesh::FaceData<char>  (&(m_bps->mesh()));
    for (size_t i = 0; i < m_bps->mesh().nv(); i++) (*m_vertex_selected)[i] = 0;
    for (size_t i = 0; i < m_bps->mesh().ne(); i++) (*m_edge_selected)[i] = 0;
    for (size_t i = 0; i < m_bps->mesh().nt(); i++) (*m_face_selected)[i] = 0;
    
    int nt = Options::intValue("num-threads");
//    if (nt > 0) omp_set_num_threads(nt);    // sim option overrides the environment variable OMP_NUM_THREADS

    
    // prepare to start the simulation
    m_time = Options::doubleValue("start-time");
    m_dt = Options::doubleValue("time-step");
    m_finished = false;
    

    // output the initial frame
    if (m_output_directory != "" && Options::boolValue("output-mesh"))
    {
        int frameid = (int)(time() / dt() + 0.5);
        std::stringstream mesh_ss;
        mesh_ss << m_output_directory << "/mesh" << std::setfill('0') << std::setw(6) << frameid << ".rec";
        MeshIO::save(*m_bps, mesh_ss.str());
    }
    
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Time stepping
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Sim::step(int substep)
{
    CSim::TimerMan::timer("Sim.step").start();
    
    assert(m_scene != "unspecified");
    assert(m_bps);
    assert(m_bps->surfTrack());
    std::cout << "Time stepping: t = " << m_time << ", dt = " << m_dt << std::endl;
    
    
    // scene-specific time stepping
    Scenes::pre_step(this, m_scene, m_dt, m_bps);

    
    // general time stepping
    double dt = m_bps->step(m_dt, substep);
    
    
    // compute the volumes
    static double s_init_volume = 0;
    double volume = m_bps->liquid_volume();
    if (s_init_volume == 0) s_init_volume = volume;
    std::cout << "Volume = " << volume << " (" << (volume > s_init_volume ? "+" : "-") << (int)(fabs(volume - s_init_volume) / s_init_volume * 100000) / 1000.0 << "%)" << std::endl;
//    std::cerr << m_time << " " << volume << std::endl;
    
    
    // scene-specific time stepping
    Scenes::post_step(this, m_scene, m_dt, m_bps);

    
    // advance time
    m_frameid++;
    m_time += m_dt;
    if (m_time >= Options::doubleValue("simulation-time"))
        m_finished = true;
    
    
    // check for a termination flag in the file system, since there's no other way to interact with a headless process.
    if (m_headless && m_output_directory != "")
    {
        std::ifstream ftest(std::string(m_output_directory + "/terminate.txt").c_str());
        if (ftest.is_open())
            m_finished = true;
    }
    
    CSim::TimerMan::timer("Sim.step").stop();
    
    addStepToStats(CSim::TimerMan::timer("Sim.step").lasttime());
}

void Sim::stepOutput()
{
    if (m_output_directory != "")
    {
        int frameid = (int)(time() / dt() + 0.5);
        
        int pngfd = Options::intValue("output-png-every-n-frames");
        if ((pngfd == 0 || frameid % pngfd == 0) && Options::boolValue("output-png") && !headless())
        {
            std::stringstream png_ss;
            png_ss << m_output_directory << "/frame" << std::setfill('0') << std::setw(6) << frameid << ".png";

            GLint viewport[4];
            glGetIntegerv(GL_VIEWPORT, viewport);   // the viewport is automatically set by GLFW to the new window client area upon user reshaping.
            int w = viewport[2];
            int h = viewport[3];

            YImage img;
            img.resize(w, h);
            glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char *)(img.data()));
            img.flip();
            img.save(png_ss.str().c_str());
        }
        
        int meshfd = Options::intValue("output-mesh-every-n-frames");
        if ((meshfd == 0 || frameid % meshfd == 0) && Options::boolValue("output-mesh"))
        {
            std::stringstream mesh_ss;
            mesh_ss << m_output_directory << "/mesh" << std::setfill('0') << std::setw(6) << frameid << ".rec";
            MeshIO::save(*m_bps, mesh_ss.str());
        }
        
        int objfd = Options::intValue("output-obj-every-n-frames");
        if ((objfd == 0 || frameid % objfd == 0) && Options::boolValue("output-obj"))
        {
            std::stringstream obj_ss;
            obj_ss << m_output_directory << "/mesh" << std::setfill('0') << std::setw(6) << frameid << ".obj";
            MeshIO::saveOBJ(*m_bps, obj_ss.str());
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Loading saved simulation
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Sim::incrementTime(int inc, bool loadframe)
{
    int current_frame = (int)((m_time + m_dt * 0.5) / m_dt);
    int next_frame = current_frame + inc;
    if (next_frame < 0) next_frame = 0;
    
    if (loadframe)
    {
        std::stringstream ss;
        ss << m_load_directory << "/mesh" << std::setfill('0') << std::setw(6) << next_frame << ".rec";
        if (!MeshIO::load(*m_bps, ss.str()))
        {
            std::cout << "Loading frame " << ss.str() << " unsuccessful." << std::endl;
            return false;
        }
        
        m_bps->m_intermediate_v.clear();
        m_bps->m_intermediate_v.push_back(std::pair<VecXd, std::string>(m_bps->velv(), "Loaded end-of-step velocity"));
        m_bps->m_intermediate_v_selector = 0;
        
        std::cout << "Loaded frame " << ss.str() << ": " << bps()->nv() << " vertices, " << bps()->nf() << " faces, liquid total volume = " << bps()->liquid_volume() << std::endl;
    }
    
    m_time = m_dt * next_frame;
    
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Rendering
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace
{
    bool is_edge_nonmanifold(const LosTopos::SurfTrack & st, size_t e)
    {
        return st.m_mesh.m_edge_to_triangle_map[e].size() != 2;
    }
    
    bool is_vertex_nonmanifold(const LosTopos::SurfTrack & st, size_t v)
    {
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v].size(); i++)
            if (is_edge_nonmanifold(st, st.m_mesh.m_vertex_to_edge_map[v][i]))
                return true;
        return false;
    }
    
    bool is_face_next_to_nonmanifold_vertices(const LosTopos::SurfTrack & st, size_t f)
    {
        return is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][0]) || is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][1]) || is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][2]);
    }
    
    bool is_edge_next_to_nonmanifold_vertices(const LosTopos::SurfTrack & st, size_t e)
    {
        return is_vertex_nonmanifold(st, st.m_mesh.m_edges[e][0]) || is_vertex_nonmanifold(st, st.m_mesh.m_edges[e][1]);
    }
    
    bool is_vertex_next_to_nonmanifold_vertices(const LosTopos::SurfTrack & st, size_t v)
    {
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v].size(); i++)
            if (is_edge_next_to_nonmanifold_vertices(st, st.m_mesh.m_vertex_to_edge_map[v][i]))
                return true;
        return false;
    }
}

void Sim::render(RenderMode rm, double field_scale, bool prestepgeometry, const Vec2d & mousepos, int selection_mask)
{
    std::vector<Vec3d> xs(m_bps->nv());
    if (prestepgeometry)
        for (size_t i = 0; i < m_bps->nv(); i++)
            xs[i] = m_bps->preStepGeometry().segment<3>(i * 3);
    else
        for (size_t i = 0; i < m_bps->nv(); i++)
            xs[i] = m_bps->pos(i);
    
    // find the primitive being picked by cursor
    Mat4d MV;
    Mat4d PJ;
    {
        float mv[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mv);
        float pj[16];
        glGetFloatv(GL_PROJECTION_MATRIX, pj);
        MV << mv[0], mv[4], mv[8], mv[12], mv[1], mv[5], mv[9], mv[13], mv[2], mv[6], mv[10], mv[14], mv[3], mv[7], mv[11], mv[15];
        PJ << pj[0], pj[4], pj[8], pj[12], pj[1], pj[5], pj[9], pj[13], pj[2], pj[6], pj[10], pj[14], pj[3], pj[7], pj[11], pj[15];
    }
    Mat4d MVP = PJ * MV;
    
    double mind = -1;
    m_nearest_vertex = -1;
    m_nearest_edge = -1;
    m_nearest_face = -1;
    if (selection_mask & SM_VERTEX)
    {
        for (size_t i = 0; i < m_bps->mesh().nv(); i++)
        {
            Vec3d pos = xs[i];
            Vec4d scrpos_h = MVP * Vec4d(pos.x(), pos.y(), pos.z(), 1.0);
            Vec2d scrpos = Vec2d(scrpos_h.x(), scrpos_h.y()) / scrpos_h.w();
            
            double distance = (scrpos - mousepos).norm();
            if (distance < mind || mind < 0)
            {
                mind = distance;
                m_nearest_vertex = i;
            }
        }
    }
    
    if (selection_mask & SM_EDGE)
    {
        for (size_t i = 0; i < m_bps->mesh().ne(); i++)
        {
            Vec3d v0 = xs[m_bps->mesh().m_edges[i][0]];
            Vec3d v1 = xs[m_bps->mesh().m_edges[i][1]];
            
            Vec4d scrv0_h = MVP * Vec4d(v0.x(), v0.y(), v0.z(), 1.0);
            Vec2d scrv0 = Vec2d(scrv0_h.x(), scrv0_h.y()) / scrv0_h.w();
            Vec4d scrv1_h = MVP * Vec4d(v1.x(), v1.y(), v1.z(), 1.0);
            Vec2d scrv1 = Vec2d(scrv1_h.x(), scrv1_h.y()) / scrv1_h.w();
            
            double distance = (mousepos - (scrv0 + scrv1) / 2).norm();
//            double distance = (mousepos - scrv0 - (mousepos - scrv0).dot(scrv1 - scrv0) * (scrv1 - scrv0) / (scrv1 - scrv0).squaredNorm()).norm();
            if (distance < mind || mind < 0)
            {
                mind = distance;
                m_nearest_vertex = -1;
                m_nearest_edge = i;
            }
        }
    }
    
    if (selection_mask & SM_FACE)
    {
        for (size_t i = 0; i < m_bps->mesh().nt(); i++)
        {
            const LosTopos::Vec3st & t = m_bps->mesh().get_triangle(i);
            Vec3d v0 = xs[t[0]];
            Vec3d v1 = xs[t[1]];
            Vec3d v2 = xs[t[2]];
            
            Vec4d scrv0_h = MVP * Vec4d(v0.x(), v0.y(), v0.z(), 1.0);
            Vec2d scrv0 = Vec2d(scrv0_h.x(), scrv0_h.y()) / scrv0_h.w();
            Vec4d scrv1_h = MVP * Vec4d(v1.x(), v1.y(), v1.z(), 1.0);
            Vec2d scrv1 = Vec2d(scrv1_h.x(), scrv1_h.y()) / scrv1_h.w();
            Vec4d scrv2_h = MVP * Vec4d(v2.x(), v2.y(), v2.z(), 1.0);
            Vec2d scrv2 = Vec2d(scrv2_h.x(), scrv2_h.y()) / scrv2_h.w();
            
            double distance = (mousepos - (scrv0 + scrv1 + scrv2) / 3).norm();
            if (distance < mind || mind < 0)
            {
                mind = distance;
                m_nearest_vertex = -1;
                m_nearest_edge = -1;
                m_nearest_face = i;
            }
        }
    }
    
    assert(mind >= 0);
    assert(m_nearest_vertex >= 0 || m_nearest_edge >= 0 || m_nearest_face >= 0);
    
    bool truncate = false;
    
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    
    if (rm == RM_TRANSPARENT || rm == RM_NONMANIFOLD)
    {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
    } else
    {
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        
        GLfloat mat_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
        glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
        GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
        GLfloat mat_shininess[] = { 50.0 };
        glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

        GLfloat light_ambient[] = { 0.3, 0.3, 0.3, 1.0 };
        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
        GLfloat light_direction[] = { 1.0, 1.0, 1.0, 0.0 };
        glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light_direction);
    }
    
    // pre-compute vertex normals (area-weighted face normals), for RM_OPAQUE_SMOOTH_SHADED mode
    std::vector<Vec3d> vn(m_bps->mesh().nv(), Vec3d(0, 0, 0));
    for (size_t i = 0; i < m_bps->mesh().nt(); i++)
    {
        LosTopos::Vec3st t = m_bps->surfTrack()->m_mesh.get_triangle(i);
        Vec3d x0 = xs[t[0]];
        Vec3d x1 = xs[t[1]];
        Vec3d x2 = xs[t[2]];
        
        Vec3d nt = (x1 - x0).cross(x2 - x0);
        if (m_bps->surfTrack()->m_mesh.get_triangle_label(i)[0] < m_bps->surfTrack()->m_mesh.get_triangle_label(i)[1])
            nt = -nt;
        
        vn[t[0]] += nt;
        vn[t[1]] += nt;
        vn[t[2]] += nt;
    }
    for (size_t i = 0; i < m_bps->mesh().nv(); i++)
        vn[i].normalize();
    
    if (false)
    {
        glLineWidth(5);
        glBegin(GL_LINES);
        for (size_t i = 0; i < m_bps->mesh().nt(); i++)
        {
            LosTopos::Vec3st t = m_bps->surfTrack()->m_mesh.get_triangle(i);
            Vec3d x0 = xs[t[0]];
            Vec3d x1 = xs[t[1]];
            Vec3d x2 = xs[t[2]];
            Vec3d c = (x0 + x1 + x2) / 3;
            Vec3d n = (x1 - x0).cross(x2 - x0).normalized();
            Vec3d eout = c + n * 0.03;
            Vec3d ein = c - n * 0.03;
            
            LosTopos::Vec2i l = m_bps->mesh().get_triangle_label(i);

            if (l[1] == 0)      glColor3d(1, 0, 0);
            else if (l[1] == 1) glColor3d(0, 1, 0);
            else if (l[1] == 2) glColor3d(0, 0, 1);
            else if (l[1] == 3) glColor3d(0.8, 0.8, 0);
            else if (l[1] == 4) glColor3d(0.8, 0, 0.8);
            else if (l[1] == 5) glColor3d(0, 0.8, 0.8);
            else if (l[1] == 6) glColor3d(0.4, 0.4, 1);
            else if (l[1] == 7) glColor3d(0.4, 1, 0.4);
            else if (l[1] == 8) glColor3d(1, 0.4, 0.4);
            else                glColor3d(0, 0, 0);
            glVertex3d(c[0], c[1], c[2]);
            glVertex3d(eout[0], eout[1], eout[2]);

            if (l[0] == 0)      glColor3d(1, 0, 0);
            else if (l[0] == 1) glColor3d(0, 1, 0);
            else if (l[0] == 2) glColor3d(0, 0, 1);
            else if (l[0] == 3) glColor3d(0.8, 0.8, 0);
            else if (l[0] == 4) glColor3d(0.8, 0, 0.8);
            else if (l[0] == 5) glColor3d(0, 0.8, 0.8);
            else if (l[0] == 6) glColor3d(0.4, 0.4, 1);
            else if (l[0] == 7) glColor3d(0.4, 1, 0.4);
            else if (l[0] == 8) glColor3d(1, 0.4, 0.4);
            else                glColor3d(0, 0, 0);
            glVertex3d(c[0], c[1], c[2]);
            glVertex3d(ein[0], ein[1], ein[2]);
        }
        glEnd();
        glLineWidth(1);

    }
    
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < m_bps->surfTrack()->m_mesh.nt(); i++)
    {
        if (rm == RM_NONMANIFOLD)
        {
            if (!is_face_next_to_nonmanifold_vertices(*m_bps->surfTrack(), i))
                continue;
        }
        
//        double r, g, b;
//        cm.getColorByDensity((gamma(i) - gamma_min) / (gamma_max - gamma_min), r, g, b);
//        glColor3d(r, g, b);

        LosTopos::Vec3st t = m_bps->surfTrack()->m_mesh.get_triangle(i);
        LosTopos::Vec3d x0 = vc(xs[t[0]]);
        LosTopos::Vec3d x1 = vc(xs[t[1]]);
        LosTopos::Vec3d x2 = vc(xs[t[2]]);
        LosTopos::Vec3d c = (x0 + x1 + x2) / 3;
        
        if (truncate)
            if (c[0] > 0.5 || c[0] < -0.5)
                continue;
        
        double shrink = (rm == RM_TRANSPARENT ? 0.05 : 0);
        x0 += (c - x0) * shrink;
        x1 += (c - x1) * shrink;
        x2 += (c - x2) * shrink;
        Vec3d n0, n1, n2;
        
        if (rm == RM_OPAQUE_FLAT_SHADED)
        {
            n0 = vc(cross(x1 - x0, x2 - x0));
            if (m_bps->surfTrack()->m_mesh.get_triangle_label(i)[0] < m_bps->surfTrack()->m_mesh.get_triangle_label(i)[1])
                n0 = -n0;
            n0.normalize();
            n1 = n0;
            n2 = n0;
        } else
        {
            n0 = vn[t[0]];
            n1 = vn[t[1]];
            n2 = vn[t[2]];
        }

        if (m_nearest_face == i)
            glColor4d(0.4, 0.5, 0.6, 0.5);
        else if (isFaceSelected(i))
            glColor4d(0.5, 0.6, 0.7, 0.5);
        else if (m_bps->face_is_solid(i))
            glColor4d(0.5, 0.0, 0.8, 0.2);
        else
            glColor4d(0.7, 0.8, 0.9, 0.2);
        
        glNormal3d(n0[0], n0[1], n0[2]);    glVertex3d(x0[0], x0[1], x0[2]);
        glNormal3d(n1[0], n1[1], n1[2]);    glVertex3d(x1[0], x1[1], x1[2]);
        glNormal3d(n2[0], n2[1], n2[2]);    glVertex3d(x2[0], x2[1], x2[2]);
    }
    glEnd();
    
    if (rm == RM_TRANSPARENT || rm == RM_NONMANIFOLD)
    {
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    } else
    {
        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);
    }
    
    // render edges
    if (rm != RM_OPAQUE_SMOOTH_SHADED)
    {
        glLineWidth(1);
        glColor3d(0, 0, 0);
        glBegin(GL_LINES);
        for (size_t i = 0; i < m_bps->mesh().ne(); i++)
        {
            if (rm == RM_NONMANIFOLD)
            {
//                if (!is_edge_next_to_nonmanifold_vertices(*m_vs->surfTrack(), i))
                if (!is_edge_nonmanifold(*m_bps->surfTrack(), i))
                    continue;
            }
            
            Vec3d x0 = xs[m_bps->mesh().m_edges[i][0]];
            Vec3d x1 = xs[m_bps->mesh().m_edges[i][1]];
            
            if (truncate)
                if (x0.x() + x1.x() > 1 || x0.x() + x1.x() < -1)
                    continue;
            
            glVertex3d(x0[0], x0[1], x0[2]);
            glVertex3d(x1[0], x1[1], x1[2]);
        }
        glEnd();
    }
    
    // render vertices, with mean curvature coloring
    if (rm != RM_OPAQUE_SMOOTH_SHADED)
    {
        glPointSize(2);
        glColor3f(0, 0, 0);
        glBegin(GL_POINTS);
        for (size_t i = 0; i < m_bps->surfTrack()->m_mesh.nv(); i++)
        {
            if (rm == RM_NONMANIFOLD)
            {
//                if (!is_vertex_next_to_nonmanifold_vertices(*m_vs->surfTrack(), i))
                if (!is_vertex_nonmanifold(*m_bps->surfTrack(), i))
                    continue;
            }
            
            Vec3d x = xs[i];
            
            if (truncate)
                if (x[0] > 0.5 || x[0] < -0.5)
                    continue;
            
            glVertex3d(x[0], x[1], x[2]);
        }
        glEnd();
        glPointSize(1);
    }
    
    // render solid vertices
    if (rm != RM_OPAQUE_SMOOTH_SHADED)
    {
        glPointSize(6);
        glColor3f(1, 0, 0.7);
        glBegin(GL_POINTS);
        for (size_t i = 0; i < m_bps->mesh().nv(); i++)
        {
            Vec3d x = xs[i];
            if (m_bps->vertex_is_solid(i))
                glVertex3d(x[0], x[1], x[2]);
        }
        glEnd();
        glPointSize(1);
    }
    
    // render selected vertices
    glPointSize(10);
    glColor3f(0, 0, 0);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < m_bps->mesh().nv(); i++)
    {
        Vec3d x = xs[i];
        if (isVertexSelected(i))
            glVertex3d(x[0], x[1], x[2]);
    }
    glEnd();
    glPointSize(1);
    
    if (m_nearest_vertex >= 0)
    {
        glPointSize(12);
        glColor3f(0, 0, 0);
        glBegin(GL_POINTS);
        Vec3d x = xs[m_nearest_vertex];
        glVertex3d(x[0], x[1], x[2]);
        glEnd();
        glPointSize(1);
    }
    
    // render velocity
    if (field_scale > 0)
    {
        VecXd vs = m_bps->intermediate_v().first;
        std::string info = m_bps->intermediate_v().second;
        
        glColor3d(0, 0.5, 0);
        glLineWidth(3);
        glBegin(GL_LINES);
        for (size_t i = 0; i < m_bps->mesh().nv(); i++)
        {
            Vec3d x = xs[i];
            Vec3d v = vs.segment<3>(i * 3);
            
            if (truncate)
                if (x[0] > 0.5 || x[0] < -0.5)
                    continue;
            
            if (rm == RM_NONMANIFOLD)
            {
                if (!is_vertex_next_to_nonmanifold_vertices(*m_bps->surfTrack(), i))
                    continue;
            }

            double s = 0.1 * field_scale;
            Vec3d e = x + v * s;
            glVertex3d(x[0], x[1], x[2]);
            glVertex3d(e[0], e[1], e[2]);
        }
        glEnd();
        glLineWidth(1);
    }
    
    // render selected edges
    glLineWidth(3);
    glColor3d(0, 0, 0);
    glBegin(GL_LINES);
    for (size_t i = 0; i < m_bps->mesh().ne(); i++)
    {
        Vec3d x0 = xs[m_bps->mesh().m_edges[i][0]];
        Vec3d x1 = xs[m_bps->mesh().m_edges[i][1]];
        if (isEdgeSelected(i))
            glVertex3d(x0[0], x0[1], x0[2]),
            glVertex3d(x1[0], x1[1], x1[2]);
    }
    glEnd();
    glLineWidth(1);
    
    if (m_nearest_edge >= 0)
    {
        glColor3d(0, 0, 0);
        glLineWidth(5);
        glBegin(GL_LINES);
        Vec3d x0 = xs[m_bps->mesh().m_edges[m_nearest_edge][0]];
        Vec3d x1 = xs[m_bps->mesh().m_edges[m_nearest_edge][1]];
        glVertex3d(x0[0], x0[1], x0[2]);
        glVertex3d(x1[0], x1[1], x1[2]);
        glEnd();
        glLineWidth(1);
    }
}

double Sim::currentInflux() const
{
    return Options::doubleValue("influx") + Options::doubleValue("influx-inc") * time();
}

void Sim::showPrimitiveInfo()
{
    if (m_nearest_vertex >= 0)
    {
        std::cout << "Vertex of Interest: " << m_nearest_vertex << " (" << m_bps->pos(m_nearest_vertex).transpose() << ") v = (" << m_bps->vel(m_nearest_vertex).transpose() << ") intermediate_v = (" << m_bps->intermediate_v().first.segment<3>(m_nearest_vertex * 3).transpose() << ") area = " << m_bps->vert_area(m_nearest_vertex) << " solid angle = " << m_bps->vert_interior_solid_angle(m_nearest_vertex) << (m_bps->vertex_is_solid(m_nearest_vertex) ? " solid" : " free") << (m_bps->vertex_is_on_triple_junction(m_nearest_vertex) ? " on triple junction" : "") << std::endl;
        std::cout << "  incident edges:"; for (size_t i = 0; i < m_bps->mesh().m_vertex_to_edge_map    [m_nearest_vertex].size(); i++) std::cout << " " << m_bps->mesh().m_vertex_to_edge_map    [m_nearest_vertex][i]; std::cout << std::endl;
        std::cout << "  incident faces:"; for (size_t i = 0; i < m_bps->mesh().m_vertex_to_triangle_map[m_nearest_vertex].size(); i++) std::cout << " " << m_bps->mesh().m_vertex_to_triangle_map[m_nearest_vertex][i]; std::cout << std::endl;
        if (m_bps->m_dbg_v1.size() == m_bps->nv())
            std::cout << "  m_dbg_v1 = " << m_bps->m_dbg_v1[m_nearest_vertex] << std::endl;
        std::cout << "  remeshing target edge length = " << m_bps->surfTrack()->m_target_edge_lengths[m_nearest_vertex] << std::endl;
    }
    
    if (m_nearest_edge >= 0)
    {
        std::cout << "Edge of Interest: " << m_nearest_edge << ": " << m_bps->mesh().m_edges[m_nearest_edge][0] << " (" << m_bps->pos(m_bps->mesh().m_edges[m_nearest_edge][0]).transpose() << ") - " << m_bps->mesh().m_edges[m_nearest_edge][1] << " (" << m_bps->pos(m_bps->mesh().m_edges[m_nearest_edge][1]).transpose() << ") length = " << m_bps->edge_length(m_nearest_edge) << (m_bps->edge_is_on_triple_junction(m_nearest_edge) ? " on triple junction" : "") << std::endl;
        std::cout << "  incident faces:"; for (size_t i = 0; i < m_bps->mesh().m_edge_to_triangle_map[m_nearest_edge].size(); i++) std::cout << " " << m_bps->mesh().m_edge_to_triangle_map[m_nearest_edge][i]; std::cout << std::endl;
    }
    
    if (m_nearest_face >= 0)
    {
        std::cout << "Face of Interest: " << m_nearest_face << ": " << m_bps->mesh().m_tris[m_nearest_face][0] << " (" << m_bps->pos(m_bps->mesh().m_tris[m_nearest_face][0]).transpose() << "), " << m_bps->mesh().m_tris[m_nearest_face][1] << " (" << m_bps->pos(m_bps->mesh().m_tris[m_nearest_face][1]).transpose() << "), " << m_bps->mesh().m_tris[m_nearest_face][2] << " (" << m_bps->pos(m_bps->mesh().m_tris[m_nearest_face][2]).transpose() << ") area = " << m_bps->face_area(m_nearest_face) << (m_bps->face_is_solid(m_nearest_face) ? " solid" : " free") << std::endl;
        std::cout << "  labels: " << m_bps->mesh().get_triangle_label(m_nearest_face) << std::endl;
        std::cout << "  incident edges:"; for (size_t i = 0; i < 3; i++) std::cout << " " << m_bps->mesh().m_triangle_to_edge_map[m_nearest_face][i]; std::cout << std::endl;
    }
}

void Sim::selectPrimitive()
{
    if (m_nearest_vertex >= 0)
    {
        if (isVertexSelected(m_nearest_vertex))
        {
            (*m_vertex_selected)[m_nearest_vertex] = 0;
            std::cout << "Deselected vertex " << m_nearest_vertex << std::endl;
        } else
        {
            (*m_vertex_selected)[m_nearest_vertex] = 1;
            std::cout << "Selected vertex " << m_nearest_vertex << std::endl;
        }
    }
    if (m_nearest_edge >= 0)
    {
        if (isEdgeSelected(m_nearest_edge))
        {
            (*m_edge_selected)[m_nearest_edge] = 0;
            std::cout << "Deselected edge " << m_nearest_edge << std::endl;
        } else
        {
            (*m_edge_selected)[m_nearest_edge] = 1;
            std::cout << "Selected edge " << m_nearest_edge << std::endl;
        }
    }
    if (m_nearest_face >= 0)
    {
        if (isFaceSelected(m_nearest_face))
        {
            (*m_face_selected)[m_nearest_face] = 0;
            std::cout << "Deselected face " << m_nearest_face << std::endl;
        } else
        {
            (*m_face_selected)[m_nearest_face] = 1;
            std::cout << "Selected face " << m_nearest_face << std::endl;
        }
    }
}

void Sim::showSelected()
{
    std::vector<size_t> s;
    s = selectedVertices();
    std::cout << s.size() << " vertices selected: "; for (size_t i = 0; i < s.size(); i++) std::cout << s[i] << " "; std::cout << std::endl;
    s = selectedEdges();
    std::cout << s.size() << " edges selected:    "; for (size_t i = 0; i < s.size(); i++) std::cout << s[i] << " "; std::cout << std::endl;
    s = selectedFaces();
    std::cout << s.size() << " faces selected:    "; for (size_t i = 0; i < s.size(); i++) std::cout << s[i] << " "; std::cout << std::endl;
}

void Sim::addStepToStats(double timecost)
{
    Stats::Step step;
    step.nv = bps()->nv();
    step.ne = bps()->mesh().ne();
    step.nf = bps()->nf();
    step.total_vol = bps()->liquid_volume();
    step.sim_time = time();
    step.real_time = timecost;
    m_stats.steps.push_back(step);
}

void Sim::statsReport()
{
    double nvmean = 0;
    double nemean = 0;
    double nfmean = 0;
    double volssq = 0;
    double volsum = 0;
    double volmean = 0;
    double timemean = 0;
    for (size_t i = 0; i < m_stats.steps.size(); i++)
    {
        const Stats::Step & s = m_stats.steps[i];
        nvmean += s.nv;
        nemean += s.ne;
        nfmean += s.nf;
        volssq += s.total_vol * s.total_vol;
        volsum += s.total_vol;
        timemean += s.real_time;
    }
    int n = m_stats.steps.size();
    nvmean /= n;
    nemean /= n;
    nfmean /= n;
    double volstd = sqrt((n * volssq - volsum * volsum) / (n * (n - 1)));
    volmean = volsum / n;
    timemean /= n;

    std::cout << "=======================================================" << std::endl;
    std::cout << "                   Sim Stats Report" << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "sim-time       nv      ne      nf        vol    time(s)" << std::endl;
    for (size_t i = 0; i < m_stats.steps.size(); i++)
    {
        const Stats::Step & s = m_stats.steps[i];
        printf("%8.3f %8d%8d%8d %10f %10f\n", s.sim_time, s.nv, s.ne, s.nf, s.total_vol, s.real_time);
    }
    std::cout << "-------------------------------------------------------" << std::endl;
    printf("mean     %8d%8d%8d %10f %10f\n", (int)(nvmean + 0.5), (int)(nemean + 0.5), (int)(nfmean + 0.5), volmean, timemean);
    printf("volume std                        %10f\n", volstd);
    std::cout << "=======================================================" << std::endl;
}

