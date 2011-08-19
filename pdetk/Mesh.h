#ifndef MESH_H
#define MESH_H

#include "MeshFD.h"

#include "Mapping.h"
namespace Daetk
{

typedef HillSlope1Map MAPPING;
//typedef HillSlope2Map MAPPING;
//typedef YotovMap MAPPING;
//typedef IdentityMap MAPPING;
//typedef RotationMap MAPPING;
typedef MeshFD<MAPPING> MESH;

}//Daetk

#endif
