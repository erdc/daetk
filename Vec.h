#ifndef VEC_H
#define VEC_H

#include <vector>
#include "Definitions.h"
#include "CMRVec.h"
#include "DaetkPetscVec.h"
#ifndef USE_SINGLE_PRECISION
namespace Daetk 
{
//typedef CMRVec<real> Vec;
  //typedef Petsc::Vec Vec;
  using Petsc::Vec;
  using PetscVecOperators::operator>>;
  using PetscVecOperators::operator<<;
  using PetscVecOperators::nrm2;
  using PetscVecOperators::dot;
  using PetscVecOperators::norm;
  using PetscVecOperators::max;
  using PetscVecOperators::min;
#else
typedef CMRVec<float> Vec;
#endif
//typedef CMRVecIndex VecIndex;
typedef Petsc::VecIndex VecIndex;

 using std::vector; 
typedef vector<Vec> VecVec;
typedef vector<VecVec> VecVecVec;
typedef vector<VecVecVec> VecVecVecVec;
}
#endif
