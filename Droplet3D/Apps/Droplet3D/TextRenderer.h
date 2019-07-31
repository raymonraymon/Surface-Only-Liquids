//
//  TextRenderer.h
//  Droplet3D
//
//  Created by Fang Da on 10/9/15.
//
//

#ifndef __Droplet3D__TextRenderer__
#define __Droplet3D__TextRenderer__

#include <iostream>
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

class TextRenderer
{
public:
    static bool initialize();
    
public:
    static bool renderString(std::string & s, float x = 0.0f, float y = 0.0f, float z = 0.0f, float size = 1.0f);
    
public:
    struct Type
    {
        enum Font
        {
            FONT1,
            FONT2,
            FONT_COUNT
        };
        
        Font font;
        
        float width;        // relative to height (specified as size)
        float increment;    // relative to height (specified as size)
        bool italic;
        bool bold;
        
        Type() : font(FONT1), width(0.5), increment(0.5), italic(false), bold(false) { }
    };
    
protected:
    static bool loadTexture(GLuint tex, const std::string & filename);
    static void setType(const Type & t);
    
protected:
    static GLuint s_texFont1;
    static GLuint s_texFont2;
    static Type s_type;
    
};

#endif /* defined(__Droplet3D__TextRenderer__) */
