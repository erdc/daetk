#include "LiBandedDirectSolver.h"

namespace Daetk 
{

LiBandedDirectSolver::LiBandedDirectSolver():
  neq(0),
  kl(0),
  ku(0),
  LEAD_DIM_STORAGE(0),
  errorFlag(true),
  arrayptr(0),
  pivotptr(0),
  x(this),
  BCM(0)
{
  Tracer tr("LiBandedDirectSolver::LiBandedDirectSolver()");
}


LiBandedDirectSolver::LiBandedDirectSolver(BandColMat& BCMin):
  neq(BCMin.getNeq()),
  kl(BCMin.getLowerBandWidth()),
  ku(BCMin.getUpperBandWidth()),
  LEAD_DIM_STORAGE(2*kl+ku+1),
  arrayptr(0),
  x(this,BCMin.getNeq()),
  BCM(&BCMin)
{  
  Tracer tr("LiBandedDirectSolver::LiBandedDirectSolver(BandColMat& BCMin)");
  arrayptr=BCMin.castToArray();
  pivotptr=new int[neq];
}

LiBandedDirectSolver::~LiBandedDirectSolver()
{
  Tracer tr("LiBandedDirectSolver::~LiBandedDirectSolver()");
  delete [] pivotptr;
}

bool LiBandedDirectSolver::prepare()
{ 
  errorFlag = 0;
  /*void dgbfa(double* abd,const int& lda,const int& n,
                      const int& ml,const int& mu,int* ipvt, int& info)*/
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
  F77NAME(dgbfa)(arrayptr,LEAD_DIM_STORAGE,neq,kl,ku,pivotptr,errorFlag);
#else
  F77NAME(sgbfa)(arrayptr,LEAD_DIM_STORAGE,neq,kl,ku,pivotptr,errorFlag);
#endif  
#else
#ifndef USE_SINGLE_PRECISION
  F77NAME(DGBFA)(arrayptr,LEAD_DIM_STORAGE,neq,kl,ku,pivotptr,errorFlag);
#else
  F77NAME(SGBFA)(arrayptr,LEAD_DIM_STORAGE,neq,kl,ku,pivotptr,errorFlag);
#endif
#endif  
  return errorFlag;
}


bool LiBandedDirectSolver::solve(const Vec& bIn,Vec& xIn)
{
  errorFlag = 0;
#ifndef USE_BLAS
  xIn=bIn;
#else
  copy(bIn,xIn);
#endif

  x.attachToTarget(xIn);

  /*void dgbsl(double* abd,const int& lda,const int& n,const int& ml,
                      const int& mu,int* ipvt,double* b,int& job);*/
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
  F77NAME(dgbsl)(arrayptr,LEAD_DIM_STORAGE,neq,kl,ku,pivotptr,x.v_.castToArray(),
        errorFlag);
#else
  F77NAME(sgbsl)(arrayptr,LEAD_DIM_STORAGE,neq,kl,ku,pivotptr,x.v_.castToArray(),
        errorFlag);
#endif
#else
#ifndef USE_SINGLE_PRECISION
  F77NAME(DGBSL)(arrayptr,LEAD_DIM_STORAGE,neq,kl,ku,pivotptr,x.v_.castToArray(),
        errorFlag);
#else
  F77NAME(SGBSL)(arrayptr,LEAD_DIM_STORAGE,neq,kl,ku,pivotptr,x.v_.castToArray(),
        errorFlag);
#endif
#endif
  
  x.restoreToTarget();

  return errorFlag;
}

}//Daetk
