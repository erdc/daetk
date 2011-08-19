#include "DaetkPetscVecBlas.h"
//need to fix the nonlocal stuff
namespace Daetk 
{
namespace PetscVecBlas 
{
  using Petsc::Vec;
extern "C"
{
  void F77NAME(drot)(const int& N, double* X, const int& INCX, double* Y,const int& INCY, double& C, 
            double& S);

  void F77NAME(dswap)(const int& N,  double* X, const int& INCX, double* Y, const int& INCY);

  void F77NAME(dscal)(const int& N,  const double& ALPHA,  double* X,  const int& INCX);

  void F77NAME(dcopy)(const int& N,const double* X, const int& INCX, double* Y, const int& INCY);

  void F77NAME(daxpy)(const int& N, const double& ALPHA, const double* X, const int& INCX, double* Y, 
             const int& INCY);

  double F77NAME(ddot)(const int& N, const double* X, const int& INCX, const double* Y, const int& INCY);

  double F77NAME(dnrm2)(const int& N, const double* X, const int& INCX);

  double F77NAME(dasum)(const int& N, const double* X, const int& INCX);

  int F77NAME(idamax)(const int& N, const double* X, const int& INCX);
}


void rot(Petsc::Vec& X, Petsc::Vec& Y, double& C, double& S)
{
  F77NAME(drot)(X.ldim(),X.castToArray(),X.stride(),Y.castToArray(),Y.stride(),C,S);
}
  
void swap(Petsc::Vec& X, Petsc::Vec& Y)
{
  F77NAME(dswap)(X.ldim(),X.castToArray(), X.stride(), Y.castToArray(), Y.stride());
}

void scal(const double& ALPHA, Petsc::Vec& X)
{
  F77NAME(dscal)(X.ldim(), ALPHA, X.castToArray(), X.stride());
}
  
void copy(const Petsc::Vec& X, Petsc::Vec& Y)
{
  F77NAME(dcopy)(X.ldim(), X.castToConstArray(), X.stride(), Y.castToArray(), 
        Y.stride());
}
 
void axpy(const double& ALPHA, const Petsc::Vec& X, Petsc::Vec& Y)
{
  F77NAME(daxpy)(X.ldim(), ALPHA, X.castToConstArray(), X.stride(), Y.castToArray(), 
        Y.stride());
}

//  double dot(const Petsc::Vec& X,const Petsc::Vec& Y)
//  {
//    return F77NAME(ddot)(X.dim(), X.castToConstArray(), X.stride(), 
//                Y.castToConstArray(), Y.stride());
//  }
 
//  double nrm2(const Petsc::Vec& X)
//  {
//    return F77NAME(dnrm2)(X.dim(), X.castToConstArray(), X.stride());
//  }
 
//  double asum(const Petsc::Vec& X)
//  {
//    return F77NAME(dasum)(X.dim(), X.castToConstArray(), X.stride());
//  }

//  int imax(const Petsc::Vec& X)
//  {
//    return F77NAME(idamax)(X.dim(), X.castToConstArray(), X.stride()) - 1 - X.base();
//  }

}//PetscVecBlas
}//Daetk
