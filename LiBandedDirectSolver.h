#ifndef LIBANDEDDIRECTSOLVER_H
#define LIBANDEDDIRECTSOLVER_H

#include "Definitions.h" 
#include "BandColMat.h"
#include "LinearSolver.h"
#include "Vec.h"
#include "VecBlas.h"

namespace Daetk 
{

#ifdef CRAYCC
extern "C" void F77NAME(DGBFA)(double* abd,const int& lda,const int& n,
                      const int& ml,const int& mu,int* ipvt, int& info);

extern "C" void F77NAME(DGBSL)(double* abd,const int& lda,const int& n,const int& ml,
                      const int& mu,int* ipvt,double* b,int& job);
extern "C" void F77NAME(SGBFA)(float* abd,const int& lda,const int& n,
                      const int& ml,const int& mu,int* ipvt, int& info);

extern "C" void F77NAME(SGBSL)(float* abd,const int& lda,const int& n,const int& ml,
                      const int& mu,int* ipvt,float* b,int& job);
#else
extern "C" void F77NAME(dgbfa)(double* abd,const int& lda,const int& n,
                      const int& ml,const int& mu,int* ipvt, int& info);

extern "C" void F77NAME(dgbsl)(double* abd,const int& lda,const int& n,const int& ml,
                      const int& mu,int* ipvt,double* b,int& job);
extern "C" void F77NAME(sgbfa)(float* abd,const int& lda,const int& n,
                      const int& ml,const int& mu,int* ipvt, int& info);

extern "C" void F77NAME(sgbsl)(float* abd,const int& lda,const int& n,const int& ml,
                      const int& mu,int* ipvt,float* b,int& job);
#endif

class LiBandedDirectSolver : public LinearSolver
{ 
public:
  LiBandedDirectSolver();
  LiBandedDirectSolver(BandColMat& BCMin);
  virtual ~LiBandedDirectSolver();
  bool prepare(); 
  bool solve(const Vec& b,Vec& x);
private:
  int neq,kl,ku;
  int LEAD_DIM_STORAGE,errorFlag;
  real* arrayptr;
  int* pivotptr;
  AttacheVec x;
  BandColMat* BCM;
};
}//Daetk  
#endif
