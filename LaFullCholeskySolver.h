#ifndef LAFULLCHOLESKYSOLVER_H
#define LAFULLCHOLESKYSOLVER_H
#include "Definitions.h" 
#include "Mat.h"
#include "LinearSolver.h"
#include "Vec.h"
#include "VecBlas.h"
#include "VectorFunction.h"

namespace Daetk 
{

extern "C" 
{
#ifdef CRAYCC
#include <fortran.h>
void F77NAME(dpotrf)(const char& uplo, const int& n, double* a, const int& lda, int& info);
void F77NAME(dpotrs)(const char& uplo, const int& n, const int& nrhs, double* a, const int& lda, 
	    double* b,const int& ldb, int& info);
void F77NAME(SPOTRF)(const char& uplo, const int& n, float* a, const int& lda, int& info);
void F77NAME(SPOTRS)(const char& uplo, const int& n, const int& nrhs, float* a, const int& lda, 
	    float* b,const int& ldb, int& info);
#else
void F77NAME(dpotrf)(const char& uplo, const int& n, double* a, const int& lda, int& info);
void F77NAME(dpotrs)(const char& uplo, const int& n, const int& nrhs, double* a, const int& lda, 
	    double* b,const int& ldb, int& info);
void F77NAME(spotrf)(const char& uplo, const int& n, float* a, const int& lda, int& info);
void F77NAME(spotrs)(const char& uplo, const int& n, const int& nrhs, float* a, const int& lda, 
	    float* b,const int& ldb, int& info);
#endif
}



class LaFullCholeskySolver : public LinearSolver
{ 
public:
  LaFullCholeskySolver();
  LaFullCholeskySolver(Mat& Min);
  virtual ~LaFullCholeskySolver();
  void storeLower();
  void storeUpper();
  bool prepare(); 
  bool solve(const Vec& b,Vec& x);
private:
  char uplo;
  int neq;
  int LEAD_DIM_STORAGE,errorFlag;
  AttacheVec x;
  real* arrayptr;
  Mat* M;
};
}//Daetk
#endif
