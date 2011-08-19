#ifndef DAETKPETSCVECBLASD_H
#define DAETKPETSCVECBLASD_H

#include "Definitions.h"
#include "DaetkPetscVec.h"

namespace Daetk 
{
namespace PetscVecBlas 
{
  using Petsc::Vec;
//need to fix the nonlocal stuff

///
void rot(Petsc::Vec& X, Petsc::Vec& Y, real& C, real& S);
/** Apply plane rotation */


///
void swap(Petsc::Vec& X, Petsc::Vec& Y);
/** X <-> Y */

///
void scal(const real& ALPHA, Petsc::Vec& X);
/** X <- alpha * X */

///
void copy(const Petsc::Vec& X, Petsc::Vec& Y);
/** Y <- X */

///
void axpy(const real& ALPHA, const Petsc::Vec& X, Petsc::Vec& Y);
/** Y <- alpha * X + Y*/

///
//real dot(const Petsc::Vec& X,const Petsc::Vec& Y);
/** ddot <- X^T * Y */

///
//real nrm2(const Petsc::Vec& X);
/** dnrm2 <- ||X||_2 */

///
//real asum(const Petsc::Vec& X);
/** dasum <- ||re(X)||_1 + ||im(X)||_1 */

///
//int imax(const Petsc::Vec& X);
/** idamax <- 1st k s.t. X(k) = max { X(i)| X.base() <= i < X.dim() }. */
}//PetscVecBlas
}//Daetk
#endif
