#include <iostream>
#include <sstream>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <GLFW/glfw3.h>

#include "Sim.h"
#include "MeshIO.h"
#include <iomanip>
#include "YImage.h"
#include "TextRenderer.h"
#include "Timer.h"

Sim g_sim(false);

struct SimControl
{
    bool headless;
    
    GLFWwindow * window;
    int win_w;
    int win_h;
    
    bool step;
    bool run;
    int substep_id;
    bool autoload;
    bool finished;
    
    double view_theta;
    double view_alpha;
    double view_dist;
    Vec3d view_center;
    double view_field_scale;
    
    double mouse_x;
    double mouse_y;
    
    bool ldrag;
    bool sldrag;
    double ldrag_start_x;
    double ldrag_start_y;
    bool rdrag;
    double rdrag_start_x;
    double rdrag_start_y;
    Sim::RenderMode render_mode;
    
    int selection_mode;
    
} g_sc;

void reset_camera()
{
    g_sc.view_theta = 0;
    g_sc.view_alpha = 0;
    g_sc.view_dist = 4;
    g_sc.view_center = Vec3d::Zero();
    g_sc.view_field_scale = 1;
    
    g_sc.mouse_x = 0;
    g_sc.mouse_y = 0;
    
    g_sc.ldrag = false;
    g_sc.sldrag = false;
    g_sc.ldrag_start_x = 0;
    g_sc.ldrag_start_y = 0;
    g_sc.rdrag = false;
    g_sc.rdrag_start_x = 0;
    g_sc.rdrag_start_y = 0;
}

void display()
{
    // object center and zoom
    Vec3d center(0, 0, 0);
    for (size_t i = 0; i < g_sim.bps()->mesh().nv(); i++)
        center += g_sim.bps()->pos(i);
    center /= g_sim.bps()->mesh().nv();
    Vec3d radius(0, 0, 0);
    for (size_t i = 0; i < g_sim.bps()->mesh().nv(); i++)
        for (size_t j = 0; j < 3; j++)
            radius[0] = std::max(radius[0], (g_sim.bps()->pos(i) - center)[0]);
    double min_d = std::max(std::max(radius[0], radius[1]), radius[2]) * 2.2;
//    g_sc.view_dist = std::max(min_d, g_sc.view_dist);
    
    glClearColor(1, 1, 1, 1);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double ar = (double)g_sc.win_w / g_sc.win_h;
    double vfh = 0.01 * 0.4;
    glFrustum(-vfh * ar, vfh * ar, -vfh, vfh, 0.01, 500);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslated(0, 0, -g_sc.view_dist);
    glRotated(-90, 1, 0, 0);
    glRotated(g_sc.view_alpha, 1, 0, 0);
    glRotated(g_sc.view_theta, 0, 0, 1);
    glTranslated(-g_sc.view_center.x(), -g_sc.view_center.y(), -g_sc.view_center.z());
    
    glBegin(GL_LINES);
    glColor3d(1, 0, 0);     glVertex3d(0, 0, 0);    glVertex3d(2, 0, 0);
    glColor3d(0, 1, 0);     glVertex3d(0, 0, 0);    glVertex3d(0, 2, 0);
    glColor3d(0, 0, 1);     glVertex3d(0, 0, 0);    glVertex3d(0, 0, 2);
    glEnd();
    
    g_sim.render(g_sc.render_mode, g_sc.view_field_scale, glfwGetKey(g_sc.window, GLFW_KEY_G) == GLFW_PRESS, Vec2d((double)g_sc.mouse_x / g_sc.win_w * 2 - 1, 1 - (double)g_sc.mouse_y / g_sc.win_h * 2), g_sc.selection_mode);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, g_sc.win_w, 0, g_sc.win_h, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    std::stringstream ss;
    ss << "T = " << g_sim.time();
    std::string s = ss.str();
    glColor4f(0.0, 0.0, 0.0, 1.0);
    TextRenderer::renderString(s, 32, g_sc.win_h - 32, 0, 32);
    
    ss.str("");
    if (g_sim.nearestVertex() >= 0)
        ss << "Picked: Vertex " << g_sim.nearestVertex();
    else if (g_sim.nearestEdge() >= 0)
        ss << "Picked: Edge " << g_sim.nearestEdge();
    else if (g_sim.nearestFace() >= 0)
        ss << "Picked: Face " << g_sim.nearestFace();
    s = ss.str();
    if (s.size() > 0)
        TextRenderer::renderString(s, 32, g_sc.win_h - 32 * 2, 0, 32);
    
    ss.str("");
    double influx = g_sim.currentInflux();
    ss << "Influx: " << influx;
    s = ss.str();
    TextRenderer::renderString(s, 32, g_sc.win_h - 32 * 3, 0, 32);
}

void idle()
{
    if (g_sc.run || g_sc.step)
    {
        g_sim.step(g_sc.substep_id);
        g_sc.step = false;
        g_sc.substep_id = 0;
        std::cout << "Finished step: T = " << g_sim.time() << std::endl;
        g_sim.stepOutput();
        
        if (g_sim.isFinished())
        {
            g_sc.finished = true;
            if (!g_sc.headless)
                glfwSetWindowShouldClose(g_sc.window, 1);
        }
    }
    
    if (g_sc.autoload)
    {
        if (!g_sim.incrementTime(10, true))
            exit(0);
        g_sim.stepOutput();
    }
}

void key(GLFWwindow * window, int k, int sc, int action, int mods)
{
    if (action == GLFW_PRESS)
    {
        if (k == GLFW_KEY_Q || k == GLFW_KEY_ESCAPE)
        {
            glfwSetWindowShouldClose(window, 1);
        } else if (k == GLFW_KEY_SPACE)
        {
            g_sc.run = !g_sc.run;
        } else if (k == GLFW_KEY_S)
        {
            g_sc.step = true;
        } else if (k >= GLFW_KEY_0 && k <= GLFW_KEY_9)
        {
            g_sc.step = true;
            g_sc.substep_id = k - GLFW_KEY_0;
        } else if (k == GLFW_KEY_M)
        {
            g_sc.render_mode = (Sim::RenderMode)(((int)g_sc.render_mode + (mods & GLFW_MOD_SHIFT ? -1 : 1)) % ((int)Sim::RM_COUNT));
            std::cout << "Render mode: " << (int)g_sc.render_mode << std::endl;
        } else if (k == GLFW_KEY_V)
        {
            g_sc.selection_mode = (!(mods & GLFW_MOD_SHIFT) ? (g_sc.selection_mode | Sim::SM_VERTEX) : (g_sc.selection_mode & ~Sim::SM_VERTEX));
            std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        } else if (k == GLFW_KEY_E)
        {
            g_sc.selection_mode = (!(mods & GLFW_MOD_SHIFT) ? (g_sc.selection_mode | Sim::SM_EDGE) : (g_sc.selection_mode & ~Sim::SM_EDGE));
            std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        } else if (k == GLFW_KEY_F)
        {
            g_sc.selection_mode = (!(mods & GLFW_MOD_SHIFT) ? (g_sc.selection_mode | Sim::SM_FACE) : (g_sc.selection_mode & ~Sim::SM_FACE));
            std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        } else if (k == GLFW_KEY_N)
        {
            g_sim.showPrimitiveInfo();
        } else if (k == GLFW_KEY_LEFT_BRACKET)      // [
        {
            int inc = -1;
            if (mods & GLFW_MOD_SHIFT) inc *= 10;
            if (mods & GLFW_MOD_ALT) inc *= 100;
            g_sim.incrementTime(inc, !(mods & GLFW_MOD_CONTROL));   // hold down control to just decrement time, without loading REC files
        } else if (k == GLFW_KEY_RIGHT_BRACKET)     // ]
        {
            int inc = 1;
            if (mods & GLFW_MOD_SHIFT) inc *= 10;
            if (mods & GLFW_MOD_ALT) inc *= 100;
            g_sim.incrementTime(inc, !(mods & GLFW_MOD_CONTROL));   // hold down control to just increment time, without loading REC files
        } else if (k == GLFW_KEY_COMMA)             // ,
        {
            g_sim.bps()->intermediate_v_select(-1);
            std::cout << "Intermediate velocity [" << g_sim.bps()->intermediate_v_selector() << "]: " << g_sim.bps()->intermediate_v().second << std::endl;
        } else if (k == GLFW_KEY_PERIOD)            // .
        {
            g_sim.bps()->intermediate_v_select(1);
            std::cout << "Intermediate velocity [" << g_sim.bps()->intermediate_v_selector() << "]: " << g_sim.bps()->intermediate_v().second << std::endl;
        } else if (k == GLFW_KEY_O)
        {
            std::stringstream obj_ss;
            obj_ss << "mesh_" << ::time(NULL) << ".obj";
            MeshIO::saveOBJ(*(g_sim.bps()), obj_ss.str().c_str());
        } else if (k == GLFW_KEY_P)
        {
            std::stringstream png_ss;
            png_ss << "screenshot_" << ::time(NULL) << ".png";
            
            int w, h;
            glfwGetFramebufferSize(window, &w, &h);
            
            YImage img;
            img.resize(w, h);
            glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char *)(img.data()));
            img.flip();
            img.save(png_ss.str().c_str());
        } else if (k == GLFW_KEY_R)
        {
            std::stringstream rec_ss;
            rec_ss << "mesh_" << ::time(NULL) << ".rec";
            MeshIO::save(*(g_sim.bps()), rec_ss.str().c_str());
        } else if (k == GLFW_KEY_L)
        {
            g_sc.autoload = true;
        } else if (k == GLFW_KEY_T)
        {
            if (mods & GLFW_MOD_SHIFT)
                g_sim.showSelected();
            else
                g_sim.selectPrimitive();
        } else if (k == GLFW_KEY_H)
        {
            reset_camera();
        } else if (k == GLFW_KEY_C)
        {
            if (g_sim.nearestVertex() >= 0)
            {
                g_sc.view_center = g_sim.bps()->pos(g_sim.nearestVertex());
            } else if (g_sim.nearestEdge() >= 0)
            {
                LosTopos::Vec2st e = g_sim.bps()->mesh().m_edges[g_sim.nearestEdge()];
                g_sc.view_center = (g_sim.bps()->pos(e[0]) + g_sim.bps()->pos(e[1])) / 2;
            } else if (g_sim.nearestFace() >= 0)
            {
                LosTopos::Vec3st f = g_sim.bps()->mesh().m_tris[g_sim.nearestFace()];
                g_sc.view_center = (g_sim.bps()->pos(f[0]) + g_sim.bps()->pos(f[1]) + g_sim.bps()->pos(f[2])) / 3;
            }
        } else if (k == GLFW_KEY_B)
        {
            g_sim.bps()->setVerbose(!g_sim.bps()->verbose());
        } else if (k == GLFW_KEY_EQUAL)
        {
            if (mods & GLFW_MOD_SHIFT)
                g_sc.view_field_scale = std::abs(g_sc.view_field_scale);    // turn the field display on (e.g. velocity field)
            else
                g_sc.view_field_scale *= 1.5;
        } else if (k == GLFW_KEY_MINUS)
        {
            if (mods & GLFW_MOD_SHIFT)
                g_sc.view_field_scale = -std::abs(g_sc.view_field_scale);   // turn the field display off (e.g. velocity field)
            else
                g_sc.view_field_scale /= 1.5;
        }
    }
}

void mouse(GLFWwindow * window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        if (mods & GLFW_MOD_SHIFT)
            g_sc.sldrag = true;
        else
            g_sc.ldrag = true;
        glfwGetCursorPos(window, &g_sc.ldrag_start_x, &g_sc.ldrag_start_y);
    } else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        g_sc.ldrag = false;
        g_sc.sldrag = false;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {
        g_sc.rdrag = true;
        glfwGetCursorPos(window, &g_sc.rdrag_start_x, &g_sc.rdrag_start_y);
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE)
    {
        g_sc.rdrag = false;
    }
}

void motion(GLFWwindow * window, double x, double y)
{
    if (g_sc.ldrag || g_sc.sldrag)
    {
        if (g_sc.ldrag)
        {
            g_sc.view_theta += (x - g_sc.ldrag_start_x) * 1.0;
            g_sc.view_alpha += (y - g_sc.ldrag_start_y) * 1.0;
            if (g_sc.view_alpha >  90) g_sc.view_alpha =  90;
            if (g_sc.view_alpha < -90) g_sc.view_alpha = -90;
        } else
        {
            double a = g_sc.view_alpha * M_PI / 180;
            double t = g_sc.view_theta * M_PI / 180;
            Vec3d right = Vec3d(cos(-t), sin(-t), 0);
            Vec3d up = Vec3d(-sin(a) * sin(-t), sin(a) * cos(-t), cos(a));
            g_sc.view_center -= (right * (x - g_sc.ldrag_start_x) - up * (y - g_sc.ldrag_start_y)) * g_sc.view_dist * 0.001;
        }
        
        g_sc.ldrag_start_x = x;
        g_sc.ldrag_start_y = y;
    }
    if (g_sc.rdrag)
    {
        g_sc.view_dist *= pow(2.0, (y - g_sc.rdrag_start_y) * 0.01);
        
        g_sc.rdrag_start_x = x;
        g_sc.rdrag_start_y = y;
    }
    
    g_sc.mouse_x = x;
    g_sc.mouse_y = y;
}

void reshape(GLFWwindow * window, int w, int h)
{
    g_sc.win_w = w;
    g_sc.win_h = h;
    // the viewport is automatically set by GLFW to the new window client area upon user reshaping.
}

void error(int error, const char * description)
{
    std::cerr << "GLFW Error: " << description << std::endl;
}

int main(int argc, char * argv[])
{
    if (argc < 2 || argc % 2 != 0)
    {
        std::cout << "Usage:\n\tDroplet3D option_file [param1 override_value1 [param2 override_value2 [...]]].\n" << std::endl;
        return 0;
    }
    
    // simulation setup
//    CSim::TimerMan::setReport(false);
    
    g_sc.run = false;
    g_sc.step = false;
    g_sc.substep_id = 0;
    g_sc.autoload = false;
    g_sc.finished = false;
    
    g_sc.win_w = 1000;
    g_sc.win_h = 1000;
    
    reset_camera();
    
    g_sc.render_mode = Sim::RM_TRANSPARENT;
    
    g_sc.selection_mode = Sim::SM_VERTEX | Sim::SM_EDGE | Sim::SM_FACE;
    
    std::vector<std::pair<std::string, std::string> > option_overrides;
    for (int i = 2; i < argc; i += 2)
    {
        std::pair<std::string, std::string> o;
        o.first = argv[i];
        assert(i + 1 < argc);
        o.second = argv[i + 1];
        option_overrides.push_back(o);
    }
    
    if (!g_sim.loadOptions(argv[1], option_overrides))
        return 1;
    
    g_sc.headless = g_sim.headless();
    
    if (!g_sc.headless)
    {
        std::cout << "Initializing main window..." << std::endl;
        
        // glut setup
        if (!glfwInit())
            assert(!"GLFW initialization failed.");
        
        glfwSetErrorCallback(error);

        glfwWindowHint(GLFW_VISIBLE, false);
        g_sc.window = glfwCreateWindow(g_sc.win_w, g_sc.win_h, "Droplet3D", NULL, NULL);
        if (!g_sc.window)
        {
            glfwTerminate();
            assert(!"GLFW window creation failed.");
        }
        glfwSetWindowPos(g_sc.window, 0, 0);
        glfwShowWindow(g_sc.window);
        
        glfwSetKeyCallback(g_sc.window, key);
        glfwSetMouseButtonCallback(g_sc.window, mouse);
        glfwSetCursorPosCallback(g_sc.window, motion);
        glfwSetWindowSizeCallback(g_sc.window, reshape);
        
        glfwMakeContextCurrent(g_sc.window);
        
        TextRenderer::initialize();

    }

    if (!g_sim.init())
        return 1;
    
    std::cout << "Initialization complete. Starting the simulation..." << std::endl;

    // main loop
    if (g_sc.headless)
    {
        g_sc.run = true;
        while (!g_sc.finished)
            idle();
    } else
    {
        while (!glfwWindowShouldClose(g_sc.window))
        {
            idle();
            
            display();
            
            glfwSwapBuffers(g_sc.window);
            
            glfwPollEvents();
        }
    }
    
    g_sim.statsReport();
    
    return 0;
}

