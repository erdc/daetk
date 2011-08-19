#include "LaFullCholeskySolver.h"

namespace Daetk 
{

using std::cerr;
using std::endl;
  
void LaFullCholeskySolver::storeLower(){uplo='L';}

void LaFullCholeskySolver::storeUpper(){uplo='U';}

LaFullCholeskySolver::LaFullCholeskySolver():
  uplo('L'),
  neq(0),
  LEAD_DIM_STORAGE(0),
  errorFlag(0),
  x(this),
  arrayptr(0),
  M(0)
{}
  
LaFullCholeskySolver::LaFullCholeskySolver(Mat& Min):
  uplo('L'),
  neq(Min.dim(1)),
  LEAD_DIM_STORAGE(Min.dim(0)),
  errorFlag(0),
  x(this,Min.dim(1)),
  arrayptr(Min.castToArray()),
  M(&Min)
{}

LaFullCholeskySolver::~LaFullCholeskySolver(){}


bool LaFullCholeskySolver::prepare()
{
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
  F77NAME(dpotrf)(uplo,neq,arrayptr,LEAD_DIM_STORAGE,errorFlag);
#else
  F77NAME(spotrf)(uplo,neq,arrayptr,LEAD_DIM_STORAGE,errorFlag);
#endif
#else
#ifndef USE_SINGLE_PRECISION
  F77NAME(DPOTRF)(uplo,neq,arrayptr,LEAD_DIM_STORAGE,errorFlag);
#else
  F77NAME(SPOTRF)(uplo,neq,arrayptr,LEAD_DIM_STORAGE,errorFlag);
#endif
#endif
  if (errorFlag!=0)
    {
      cerr<<"error in cholesky factorization, code "<<errorFlag<<endl;
    }
  return errorFlag;
}
  
bool LaFullCholeskySolver::solve(const Vec& bIn,Vec& xIn)
{
  int numberOfRhs(1),ldaRhs(x.v_.dim());

#ifndef USE_BLAS
  xIn=bIn;
#else
  copy(bIn,xIn);
#endif

  x.attachToTarget(xIn);

#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
  F77NAME(dpotrs)(uplo,neq,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,x.v_.castToArray(),ldaRhs,errorFlag);
#else
  F77NAME(spotrs)(uplo,neq,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,x.castToArray(),ldaRhs,errorFlag);
#endif
#else
 #ifndef USE_SINGLE_PRECISION
  F77NAME(DPOTRS)(uplo,neq,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,x.castToArray(),ldaRhs,errorFlag);
#else
  F77NAME(SPOTRS)(uplo,neq,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,x.castToArray(),ldaRhs,errorFlag);
#endif
#endif
  x.restoreToTarget();
  if (errorFlag!=0)
    {
      cerr<<"error in cholesky back substitution, code "<<errorFlag<<endl;
    }
  return errorFlag;
}

}//Daetk
