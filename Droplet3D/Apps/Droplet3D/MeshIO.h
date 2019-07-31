//
//  MeshIO.h
//
//  Christopher Batty, Fang Da 2014
//
//

#ifndef __MeshIO__
#define __MeshIO__

#include <iostream>
#include <vector>
#include "surftrack.h"
#include "BPS3D.h"

class MeshIO
{
public:
    static bool save(BPS3D & bps, const std::string & filename, bool binary = true);
    static bool load(BPS3D & bps, const std::string & filename, bool binary = true);
    
    static bool loadIntoRaw(std::vector<LosTopos::Vec3d> & xs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<LosTopos::Vec3d> & vels, std::vector<size_t> & solid_vertices, const std::string & filename, bool binary = true);
    
    static bool saveOBJ(BPS3D & bps, const std::string & filename);
    
};



#endif /* defined(__MeshIO__) */
