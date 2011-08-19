#include "CMRVecBlasd.h"

extern "C"
{
  void drot(const int& N, double* X, const int& INCX, double* Y,const int& INCY, double& C, 
            double& S);

  void dswap(const int& N,  double* X, const int& INCX, double* Y, const int& INCY);

  void dscal(const int& N,  const double& ALPHA,  double* X,  const int& INCX);

  void dcopy(const int& N,const double* X, const int& INCX, double* Y, const int& INCY);

  void daxpy(const int& N, const double& ALPHA, const double* X, const int& INCX, double* Y, 
             const int& INCY);

  double ddot(const int& N, const double* X, const int& INCX, const double* Y, const int& INCY);

  double dnrm2(const int& N, const double* X, const int& INCX);

  double dasum(const int& N, const double* X, const int& INCX);

  int idamax(const int& N, const double* X, const int& INCX);
}


void drot(CMRVecd& X, CMRVecd& Y, double& C, double& S)
{
  drot(X.dim(),X.castToArray(),X.stride(),Y.castToArray(),Y.stride(),C,S);
}
  
void dswap(CMRVecd& X, CMRVecd& Y)
{
  dswap(X.dim(),X.castToArray(), X.stride(), Y.castToArray(), Y.stride());
}

void dscal(const double& ALPHA, CMRVecd& X)
{
  dscal(X.dim(), ALPHA, X.castToArray(), X.stride());
}
  
void dcopy(const CMRVecd& X, CMRVecd& Y)
{
  dcopy(X.dim(), X.castToConstArray(), X.stride(), Y.castToArray(), 
        Y.stride());
}
 
void daxpy(const double& ALPHA, const CMRVecd& X, CMRVecd& Y)
{
  daxpy(X.dim(), ALPHA, X.castToConstArray(), X.stride(), Y.castToArray(), 
        Y.stride());
}

double ddot(const CMRVecd& X,const CMRVecd& Y)
{
  return ddot(X.dim(), X.castToConstArray(), X.stride(), 
              Y.castToConstArray(), Y.stride());
}
 
double dnrm2(const CMRVecd& X)
{
  return dnrm2(X.dim(), X.castToConstArray(), X.stride());
}
 
double dasum(const CMRVecd& X)
{
  return dasum(X.dim(), X.castToConstArray(), X.stride());
}

int idamax(const CMRVecd& X)
{
  return idamax(X.dim(), X.castToConstArray(), X.stride()) - 1 - X.base();
}








