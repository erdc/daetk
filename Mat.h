#ifndef MAT_H
#define MAT_H

#include "Definitions.h"
#include "mvm.h"

namespace Daetk 
{
#ifndef USE_SINGLE_PRECISION
typedef MV_ColMat Mat;
#else
typedef MV_ColMat Mat;
#endif
}
#endif
