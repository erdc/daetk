#include "CMRVecBlass.h"
namespace Daetk 
{

extern "C"
{
#ifndef CRAYCC
  void F77NAME(srot)(const int * N, float* X, const int* INCX, float* Y,const int* INCY, float* C, 
            float* S);

  void F77NAME(sswap)(const int* N,  float* X, const int* INCX, float* Y, const int* INCY);

  void F77NAME(sscal)(const int* N,  const float* ALPHA,  float* X,  const int* INCX);

  void F77NAME(scopy)(const int* N,const float* X, const int* INCX, float* Y, const int* INCY);

  void F77NAME(saxpy)(const int* N, const float* ALPHA, const float* X, const int* INCX, float* Y, 
             const int* INCY);

  float F77NAME(sdot)(const int* N, const float* X, const int* INCX, const float* Y, const int* INCY);

  float F77NAME(snrm2)(const int* N, const float* X, const int* INCX);

  float F77NAME(sasum)(const int* N, const float* X, const int* INCX);

  int F77NAME(isamax)(const int* N, const float* X, const int* INCX);
#else
  void F77NAME(SROT)(const int * N, float* X, const int* INCX, float* Y,const int* INCY, float* C, 
            float* S);

  void F77NAME(SSWAP)(const int* N,  float* X, const int* INCX, float* Y, const int* INCY);

  void F77NAME(SSCAL)(const int* N,  const float* ALPHA,  float* X,  const int* INCX);

  void F77NAME(SCOPY)(const int* N,const float* X, const int* INCX, float* Y, const int* INCY);

  void F77NAME(SAXPY)(const int* N, const float* ALPHA, const float* X, const int* INCX, float* Y, 
             const int* INCY);

  float F77NAME(SDOT)(const int* N, const float* X, const int* INCX, const float* Y, const int* INCY);

  float F77NAME(SNRM2)(const int* N, const float* X, const int* INCX);

  float F77NAME(SASUM)(const int* N, const float* X, const int* INCX);

  int F77NAME(ISAMAX)(const int* N, const float* X, const int* INCX);
#endif
}


void rot(CMRVecs& X, CMRVecs& Y, float& C, float& S)
{
  int dim(X.dim()),xstride(X.stride()),ystride(Y.stride());
  float *x(X.castToArray()),*y(Y.castToArray());
#ifndef CRAYCC
  F77NAME(srot)(&dim,x,&xstride,y,&ystride,&C,&S);
#else
  F77NAME(SROT)(&dim,x,&xstride,y,&ystride,&C,&S);
#endif
}
  
void swap(CMRVecs& X, CMRVecs& Y)
{
  int dim(X.dim()),xstride(X.stride()),ystride(Y.stride());
  float *x(X.castToArray()),*y(Y.castToArray());
#ifndef CRAYCC
  F77NAME(sswap)(&dim,x, &xstride,y,&ystride);
#else
  F77NAME(SSWAP)(&dim,x, &xstride,y,&ystride);
#endif
}

void scal(const float& ALPHA, CMRVecs& X)
{
  int dim(X.dim()),xstride(X.stride());
  float *x(X.castToArray()),alpha(ALPHA);
#ifndef CRAYCC
  F77NAME(sscal)(&dim, &alpha, x, &xstride);
#else
  F77NAME(SSCAL)(&dim, &alpha, x, &xstride);
#endif
}
  
void copy(const CMRVecs& X, CMRVecs& Y)
{
  int dim(X.dim()),xstride(X.stride()),ystride(Y.stride());
  const float *x(X.castToConstArray());
  float *y(Y.castToArray());
#ifndef CRAYCC
  F77NAME(scopy)(&dim,x, &xstride,y, &ystride);
#else
  F77NAME(SCOPY)(&dim,x, &xstride,y, &ystride);
#endif
}
 
void axpy(const float& ALPHA, const CMRVecs& X, CMRVecs& Y)
{
  int dim(X.dim()),xstride(X.stride()),ystride(Y.stride());
  const float *x(X.castToConstArray()),alpha(ALPHA);
  float *y(Y.castToArray());
#ifndef CRAYCC
  F77NAME(saxpy)(&dim, &alpha, x, &xstride, y, &ystride);
#else
  F77NAME(SAXPY)(&dim, &alpha, x, &xstride, y, &ystride);
#endif
}

float dot(const CMRVecs& X,const CMRVecs& Y)
{
  int dim(X.dim()),xstride(X.stride()),ystride(Y.stride());
  const float *x(X.castToConstArray()),*y(Y.castToConstArray());
#ifndef CRAYCC
  return F77NAME(sdot)(&dim, x, &xstride, y, &ystride);
#else
  return F77NAME(SDOT)(&dim, x, &xstride, y, &ystride);
#endif
}
 
float nrm2(const CMRVecs& X)
{
  int dim(X.dim()),xstride(X.stride());
  const float *x(X.castToConstArray());
#ifndef CRAYCC
  return F77NAME(snrm2)(&dim,x,&xstride);
#else
  return F77NAME(SNRM2)(&dim,x,&xstride);
#endif
}
 
float sum(const CMRVecs& X)
{
  int dim(X.dim()),xstride(X.stride());
  const float *x(X.castToConstArray());
#ifndef CRAYCC
  return F77NAME(sasum)(&dim, x, &xstride);
#else
  return F77NAME(SASUM)(&dim, x, &xstride);
#endif
}

int imax(const CMRVecs& X)
{
  int dim(X.dim()),xstride(X.stride());
  const float *x(X.castToConstArray());
#ifndef CRAYCC
  return F77NAME(isamax)(&dim, x, &xstride) - 1 - X.base();
#else
  return F77NAME(ISAMAX)(&dim, x, &xstride) - 1 - X.base();
#endif
}

}//Daetk
