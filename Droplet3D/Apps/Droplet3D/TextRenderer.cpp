//
//  TextRenderer.cpp
//  Droplet3D
//
//  Created by Fang Da on 10/9/15.
//
//

#include "TextRenderer.h"
#include <png.h>
#include <assert.h>

GLuint TextRenderer::s_texFont1;
GLuint TextRenderer::s_texFont2;
TextRenderer::Type TextRenderer::s_type;

bool TextRenderer::initialize()
{
    glGenTextures(1, &s_texFont1);
    glGenTextures(1, &s_texFont2);

    loadTexture(s_texFont1, "assets/fonts/font1.png");
    loadTexture(s_texFont2, "assets/fonts/font2.png");
}

bool TextRenderer::renderString(std::string & s, float x, float y, float z, float size)
{
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, s_texFont2);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);
    glBegin(GL_QUADS);
    for (size_t i = 0; i < s.size(); i++)
    {
        char c = s[i];
        float u0 = (c % 16) / 16.0f;
        float v0 = (15 - c / 16) / 16.0f;
        float u1 = u0 + 1 / 16.0f;
        float v1 = v0 + 1 / 16.0f;
        
        float x0 = x + size * s_type.increment * i;
        float y0 = y;
        float x1 = x0 + size * s_type.width;
        float y1 = y + size;
        
        int rep = (s_type.bold ? 1 : 0);
        float offset = size / 32.0f;
        float topoffset = size * (s_type.italic ? 1 : 0) * 0.5f;
        for (int j = -rep; j <= rep; j++)
        {
            for (int k = -rep; k <= rep; k++)
            {
                glTexCoord2f(u0, v0);   glVertex3f(x0 + offset * j, y0 + offset * k, z);
                glTexCoord2f(u1, v0);   glVertex3f(x1 + offset * j, y0 + offset * k, z);
                glTexCoord2f(u1, v1);   glVertex3f(x1 + offset * j, y1 + offset * k + topoffset, z);
                glTexCoord2f(u0, v1);   glVertex3f(x0 + offset * j, y1 + offset * k + topoffset, z);
            }
        }
    }
    glEnd();
    
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);
}

bool TextRenderer::loadTexture(GLuint tex, const std::string & filename)
{
    //
    // PNG loader code adapted from:
    //  http://zarb.org/~gc/html/libpng.html
    //
    png_byte header[8];    // 8 is the maximum size that can be checked
    
    /* open file and test for it being a png */
    FILE * fp = fopen(filename.c_str(), "rb");
    if (!fp)
    {
        std::cout << "Error: " << "File " << filename << " could not be opened for reading." << std::endl;
        return false;
    }
    
    fread(header, 1, 8, fp);
    if (png_sig_cmp(header, 0, 8))
    {
        std::cout << "Error: " << "File " << filename << " is not recognized as a PNG file." << std::endl;
        return false;
    }
    
    /* initialize stuff */
    png_structp png_ptr;
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!png_ptr)
    {
        std::cout << "Error: " << "png_create_read_struct failed." << std::endl;
        return false;
    }
    
    png_infop info_ptr;
    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        std::cout << "Error: " << "png_create_info_struct failed." << std::endl;
        return false;
    }
    
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        std::cout << "Error: " << "Error during init_io." << std::endl;
        return false;
    }
    
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);
    
    png_read_info(png_ptr, info_ptr);
    
    int w = png_get_image_width(png_ptr, info_ptr);
    int h = png_get_image_height(png_ptr, info_ptr);
    //color_type = png_get_color_type(png_ptr, info_ptr);
    //bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    assert(png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_GRAY);
    
    int rowbytes = png_get_rowbytes(png_ptr, info_ptr);
    
    //number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);
    
    /* read file */
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        std::cout << "Error: " << "Error during read_image." << std::endl;
        return false;
    }
    
    unsigned char * image_data = new unsigned char[w * h];
    
    png_bytep * row_pointers;
    row_pointers = (png_bytep *)malloc(sizeof(png_bytep) * h);
    for (int y = 0; y < h; y++)
        row_pointers[y] = (png_byte *)(image_data + rowbytes * (h - 1 - y));
    
    png_read_image(png_ptr, row_pointers);
    free(row_pointers);
    
    fclose(fp);
    
    // construct the texture data from the image data
    unsigned char * texture_data = new unsigned char[w * h * 4];
    for (int i = 0; i < w * h; i++)
    {
        texture_data[i * 4 + 0] = 0xff;
        texture_data[i * 4 + 1] = 0xff;
        texture_data[i * 4 + 2] = 0xff;
        texture_data[i * 4 + 3] = image_data[i];
    }

    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, texture_data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    
    delete image_data;
    
    return true;
}
