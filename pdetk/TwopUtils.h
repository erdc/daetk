#ifndef TWOP_UTILS_H
#define TWOP_UTILS_H

#include "Definitions.h"
#include "ParameterDatabase.h"
#include "Vec.h"
#include "PetscSecondOrderFd.h"
#include "Psk2p.h"

namespace Daetk 
{
namespace TwoPhaseFlow
{
namespace Utilities
{
//--------------------------------------------------
//global utility functions
//--------------------------------------------------
void useScaleHeterogeneity(Psk2p& psk,
			   ParameterDatabase& pd,
			   Petsc::SecondOrderFd* coarseNode, 
			   const Vec& deltaIn,
			   Petsc::SecondOrderFd& node, 
			   const Vec& vecExample,
			   const real& dx,
			   const real& dy,
			   const real& dz);


inline
int round(double x)
{ int r((int)(x)); if ((x - (real)(r)) >= 0.5) r+=1; return r;} 

}//Utilities
}//TwoPhaseFlow
}//Daetk

#endif
