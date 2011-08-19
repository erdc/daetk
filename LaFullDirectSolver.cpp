#include "LaFullDirectSolver.h"

namespace Daetk 
{

void LaFullDirectSolver::calculateCondition(DataCollector& d)
{
  data = &d;
  CALCULATE_CONDITION=true; 
  rwork = new real[3*neq]; 
  iwork = new int[neq];
}

LaFullDirectSolver::LaFullDirectSolver():
  neq(0),
  LEAD_DIM_STORAGE(0),
  errorFlag(true),
  arrayptr(0),
  rwork(0),
  pivotptr(0),
  iwork(0),
  x(this),
  M(0)
{
  Tracer tr("LaFullDirectSolver::LaFullDirectSolver()");
}
 

void LaFullDirectSolver::printMatrices(const char*filename)
{PRINT_MATRICES = true; matOut.open(filename);}


LaFullDirectSolver::LaFullDirectSolver(Mat& Min):
  PRINT_MATRICES(false),
  CALCULATE_CONDITION(false),
  neq(Min.dim(0)),
  LEAD_DIM_STORAGE(Min.dim(0)),
  arrayptr(Min.castToArray()),
  rwork(0),
  iwork(0),
  x(this),
  M(&Min)
{  
  Tracer tr("LaFullDirectSolver::LaFullDirectSolver(Mat& Min)");
  arrayptr=Min.castToArray();
  pivotptr=Min.pivotPtr();
}

void LaFullDirectSolver::attachMat(Mat& Min)
{
  M = &Min;
  arrayptr = Min.castToArray();
  pivotptr = Min.pivotPtr();
  neq = Min.dim(0);
  x.v_.newsize(neq);
  LEAD_DIM_STORAGE = Min.dim(0);
}

LaFullDirectSolver::~LaFullDirectSolver()
{
  Tracer tr("LaFullDirectSolver::~LaFullDirectSolver()");
  delete [] rwork; 
  delete [] iwork;
}

bool LaFullDirectSolver::prepare()
{
  real anorm=0.0;
  if (CALCULATE_CONDITION)
    {
      real sum=0;
      for (int i=0;i<neq;i++)
        {
          for (int j=0;j<neq;j++)
            sum+=fabs((*M)(i,j));
          if (sum > anorm)
            anorm = sum;
          sum=0.0;
        }
    }

  if (PRINT_MATRICES)
    matOut<<*M<<std::endl;

  /*
    void dgetrf(const int& m, const int& n, double* A, const int& lda,
    int* ipiv, int& info);
    */
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
  F77NAME(dgetrf)(neq,neq,arrayptr,LEAD_DIM_STORAGE, pivotptr,errorFlag);
#else
  F77NAME(sgetrf)(neq,neq,arrayptr,LEAD_DIM_STORAGE, pivotptr,errorFlag);
#endif
#else
#ifndef USE_SINGLE_PRECISION
  F77NAME(DGETRF)(&neq,&neq,arrayptr,&LEAD_DIM_STORAGE, pivotptr,&errorFlag);
#else
  F77NAME(SGETRF)(&neq,&neq,arrayptr,&LEAD_DIM_STORAGE, pivotptr,&errorFlag);
#endif
#endif
  if (errorFlag != 0)
    std::cerr<<"error in factor "<<errorFlag<<std::endl;

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
      F77NAME(dgecon)(nrm,neq,arrayptr,LEAD_DIM_STORAGE,&anorm,&rcond,rwork,iwork,errorFlag);
#else
      F77NAME(sgecon)(nrm,neq,arrayptr,LEAD_DIM_STORAGE,&anorm,&rcond,rwork,iwork,errorFlag);
#endif
#else
#ifndef USE_SINGLE_PRECISION
      F77NAME(DGECON)(nrm,neq,arrayptr,LEAD_DIM_STORAGE,&anorm,&rcond,rwork,iwork,errorFlag);
#else
      F77NAME(SGECON)(nrm,neq,arrayptr,LEAD_DIM_STORAGE,&anorm,&rcond,rwork,iwork,errorFlag);
#endif
#endif
      rcond=1./rcond;
      data->conditionNumber(rcond);
    }    
  return errorFlag;
}
 
bool LaFullDirectSolver::solve(const Vec& bIn, Vec& xIn)
{
  int numberOfRhs(1);
#ifndef USE_BLAS
  xIn=bIn;
#else
  copy(bIn,xIn);
#endif

  x.attachToTarget(xIn);
  //opt
  //copy(b,x);
  //
#ifndef CRAYCC
  const char* n="N";
#else
  char* N="N";
  _fcd n = _cptofcd(N,1);
#endif

  /*
    void dgetrs(const char* trans, const int& n, const int& nrhs, 
    double * A, const int& lda, int* ipiv, double* b, 
    const int& ldb, int& info);
    */
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
  F77NAME(dgetrs)(n,neq,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,pivotptr,x.v_.castToArray(),
         neq,errorFlag);
#else
  F77NAME(sgetrs)(n,neq,numberOfRhs,arrayptr,LEAD_DIM_STORAGE,pivotptr,x.v_.castToArray(),
         neq,errorFlag);
#endif
#else
#ifndef USE_SINGLE_PRECISION
  F77NAME(DGETRS)(n,&neq,&numberOfRhs,arrayptr,&LEAD_DIM_STORAGE,pivotptr,x.v_.castToArray(),
         neq,errorFlag);
#else
  F77NAME(SGETRS)(n,&neq,&numberOfRhs,arrayptr,&LEAD_DIM_STORAGE,pivotptr,x.v_.castToArray(),
         &neq,&errorFlag);
#endif
#endif
  x.restoreToTarget();
  return errorFlag;
}

}//Daetk
