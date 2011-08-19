#ifndef TEST_POST_PROC_H
#define TEST_POST_PROC_H

#include "DaetkPetscVec.h"
#include "PetscSecondOrderFd.h"
#include "PetscStencilMM.h"

namespace Daetk
{
namespace Petsc
{

//  //try to do this like chris' << routine and ex8 in vec/examples
//  void postProcessFluxesParallel(SecondOrderFd& stencil,
//  			       Vec& Qx,Vec& Qy,Vec& Qz);
//  //try to do this like chris' << routine and ex8 in vec/examples
//  void postProcessFluxesParallel(StencilMM& stencil,
//  			       Vec& Qx,Vec& Qy,Vec& Qz);

//try to just collect L2 error values for fluxes from all of the
//processors?
void printOutFluxError(real L2errQ, real L2exQ); 
void printOutError(real L2errQ, real L2exQ,
		   const char* nameErr, const char* nameEx); 
}//Petsc
}//Daetk
#endif
