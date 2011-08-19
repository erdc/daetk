#include "CMRVecBlasd.h"
namespace Daetk 
{

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


void rot(CMRVecd& X, CMRVecd& Y, double& C, double& S)
{
  F77NAME(drot)(X.dim(),X.castToArray(),X.stride(),Y.castToArray(),Y.stride(),C,S);
}
  
void swap(CMRVecd& X, CMRVecd& Y)
{
  F77NAME(dswap)(X.dim(),X.castToArray(), X.stride(), Y.castToArray(), Y.stride());
}

void scal(const double& ALPHA, CMRVecd& X)
{
  F77NAME(dscal)(X.dim(), ALPHA, X.castToArray(), X.stride());
}
  
void copy(const CMRVecd& X, CMRVecd& Y)
{
  F77NAME(dcopy)(X.dim(), X.castToConstArray(), X.stride(), Y.castToArray(), 
        Y.stride());
}
 
void axpy(const double& ALPHA, const CMRVecd& X, CMRVecd& Y)
{
  F77NAME(daxpy)(X.dim(), ALPHA, X.castToConstArray(), X.stride(), Y.castToArray(), 
        Y.stride());
}

double dot(const CMRVecd& X,const CMRVecd& Y)
{
  return F77NAME(ddot)(X.dim(), X.castToConstArray(), X.stride(), 
              Y.castToConstArray(), Y.stride());
}
 
double nrm2(const CMRVecd& X)
{
  return F77NAME(dnrm2)(X.dim(), X.castToConstArray(), X.stride());
}
 
double asum(const CMRVecd& X)
{
  return F77NAME(dasum)(X.dim(), X.castToConstArray(), X.stride());
}

int imax(const CMRVecd& X)
{
  return F77NAME(idamax)(X.dim(), X.castToConstArray(), X.stride()) - 1 - X.base();
}

}//Daetk
