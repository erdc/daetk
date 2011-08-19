#include "LaBandedDirectSolver.h"

namespace Daetk 
{

using std::max;
using std::min;
using std::endl;
using std::cerr;
using std::cout;

void LaBandedDirectSolver::calculateCondition(DataCollector& d)
{
  data = &d;
  CALCULATE_CONDITION=true; 
  rwork = new real[3*neq]; 
  iwork = new int[neq];
}

LaBandedDirectSolver::LaBandedDirectSolver():
  PRINT_MATRICES(0),
  CALCULATE_CONDITION(0),
  x(this,0),
  neq(0),
  kl(0),
  ku(0),
  LEAD_DIM_STORAGE(0),
  errorFlag(true),
  arrayptr(0),
  rwork(0),
  pivotptr(0),
  iwork(0),
  BCM(0)
{
  Tracer tr("LaBandedDirectSolver::LaBandedDirectSolver()");
}

void LaBandedDirectSolver::printMatrices(const char*filename)
{PRINT_MATRICES = true; matOut.open(filename);}

void LaBandedDirectSolver::attachMat(BandColMat& BCMin)
{
  BCM = &BCMin;
  neq = BCM->getNeq();
  x.v_.newsize(neq);
  kl = BCM->getLowerBandWidth();
  ku = BCM->getUpperBandWidth();
  LEAD_DIM_STORAGE= 2*kl+ku+1;
  arrayptr = BCM->castToArray();
  pivotptr = BCM->pivotPtr();
}

LaBandedDirectSolver::LaBandedDirectSolver(BandColMat& BCMin):LinearSolver(),
  PRINT_MATRICES(0),
  CALCULATE_CONDITION(0),
  x(this,BCMin.getNeq()),
  neq(BCMin.getNeq()),
  kl(BCMin.getLowerBandWidth()),
  ku(BCMin.getUpperBandWidth()),
  LEAD_DIM_STORAGE(2*kl+ku+1),
  arrayptr(0),
  rwork(0),
  iwork(0),
  BCM(&BCMin)
{  
  Tracer tr("LaBandedDirectSolver::LaBandedDirectSolver(BandColMat& BCMin)");
  arrayptr=BCMin.castToArray();
  pivotptr=BCMin.pivotPtr();
}

LaBandedDirectSolver::~LaBandedDirectSolver()
{
  Tracer tr("LaBandedDirectSolver::~LaBandedDirectSolver()");
  delete [] rwork; 
  delete [] iwork;
}

bool LaBandedDirectSolver::prepare()
{   
  real anorm=0.0;
  if (CALCULATE_CONDITION)
    {
      real sum=0;
      for (int i=0;i<neq;i++)
        {
          for (int j=max(0,i-kl);j<min(neq-1,i+kl);j++)
            sum+=fabs((*BCM)(i,j));
          if (sum > anorm)
            anorm = sum;
          sum=0.0;
        }
    }
  
  if (PRINT_MATRICES)
    matOut<<*BCM<<endl;
  /*  dgbtrf(const int&,const int&,const int&,const int&,double*,const int&,
         int*,int&);*/
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
  F77NAME(dgbtrf)(neq,neq,kl,ku,arrayptr,LEAD_DIM_STORAGE,pivotptr,errorFlag);
#else
  F77NAME(sgbtrf)(neq,neq,kl,ku,arrayptr,LEAD_DIM_STORAGE,pivotptr,errorFlag);
#endif
#else
#ifndef USE_SINGLE_PRECISION
  F77NAME(DGBTRF)(neq,neq,kl,ku,arrayptr,LEAD_DIM_STORAGE,pivotptr,errorFlag);
#else
  F77NAME(SGBTRF)(neq,neq,kl,ku,arrayptr,LEAD_DIM_STORAGE,pivotptr,errorFlag);
#endif
#endif
  if (errorFlag != 0)
    cerr<<"error in factor "<<errorFlag<<endl;
  if (CALCULATE_CONDITION)
    {
#ifndef CRAYCC
  const char* nrm="1";
#else
  char* one="1";
  _fcd nrm = _cptofcd(one,1);
#endif
      real rcond;
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
      F77NAME(dgbcon)(nrm,neq,kl,ku,arrayptr,LEAD_DIM_STORAGE,pivotptr,&anorm,&rcond,rwork,iwork,errorFlag);
#else
      F77NAME(sgbcon)(nrm,neq,kl,ku,arrayptr,LEAD_DIM_STORAGE,pivotptr,&anorm,&rcond,rwork,iwork,errorFlag);
#endif
#else
#ifndef USE_SINGLE_PRECISION
      F77NAME(DGBCON)(nrm,neq,kl,ku,arrayptr,LEAD_DIM_STORAGE,pivotptr,&anorm,&rcond,rwork,iwork,errorFlag);
#else
      F77NAME(SGBCON)(nrm,neq,kl,ku,arrayptr,LEAD_DIM_STORAGE,pivotptr,&anorm,&rcond,rwork,iwork,errorFlag);
#endif
#endif
      rcond=1./rcond;
      data->conditionNumber(rcond);
    }    

  return errorFlag;
}


bool LaBandedDirectSolver::solve(const Vec& bIn,Vec& xIn)
{
  int numberOfRhs(1);  
#ifndef USE_BLAS
  xIn=bIn;
#else
  copy(bIn,xIn);
#endif
  
  x.attachToTarget(xIn);

#ifndef CRAYCC
  const char* n="N";
#else
  char* N="N";
  _fcd n = _cptofcd(N,1);
#endif
  /* dgbtrs(const char*,const int&,const int&,const int&,const int&,
         double*,const int&,int*,double*,const int&,int&);*/
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
  F77NAME(dgbtrs)(n,neq,kl,ku,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,pivotptr,
         x.v_.castToArray(),neq,errorFlag);
#else
  F77NAME(sgbtrs)(n,neq,kl,ku,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,pivotptr,
         x.v_.castToArray(),neq,errorFlag);
#endif
#else
#ifndef USE_SINGLE_PRECISION
  F77NAME(DGBTRS)(n,neq,kl,ku,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,pivotptr,
         x.v_.castToArray(),neq,errorFlag);
#else
  F77NAME(SGBTRS)(n,neq,kl,ku,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,pivotptr,
         x.v_.castToArray(),neq,errorFlag);
#endif
#endif

  x.restoreToTarget();

  return errorFlag;
}

}//Daetk
