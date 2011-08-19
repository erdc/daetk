#ifndef LAFULLDIRECTSOLVER_H
#define LAFULLDIRECTSOLVER_H

#include <fstream>
#include "Definitions.h"
#include "IntVec.h"
#include "Mat.h"
#include "LinearSolver.h"
#include "DataCollector.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{

#ifdef CRAYCC
#include <fortran.h>
extern "C" void F77NAME(DGETRF)(const int& m, const int* n, double* A, const int* lda,
                       int* ipiv, int* info);

extern "C" void F77NAME(DGETRS)(_fcd trans, const int* n, const int* nrhs, 
                       double * A, const int* lda, int* ipiv, double* b, 
                       const int* ldb, int* info);
extern "C" void F77NAME(DGECON)(_fcd,const int&, 
                                double*, const int&, double*, 
                                double*,double*,int*,const int&);
extern "C" void F77NAME(SGETRF)(const int* m, const int* n, float* A, const int* lda,
                       int* ipiv, int* info);

extern "C" void F77NAME(SGETRS)(_fcd trans, const int* n, const int* nrhs, 
                       float * A, const int* lda, int* ipiv, float* b, 
                       const int* ldb, int* info);
#else
extern "C" void F77NAME(dgetrf)(const int& m, const int& n, double* A, const int& lda,
                       int* ipiv, int& info);

extern "C" void F77NAME(dgetrs)(const char* trans, const int& n, const int& nrhs, 
                       double * A, const int& lda, int* ipiv, double* b, 
                       const int& ldb, int& info);
extern "C" void F77NAME(dgecon)(const char*,const int&, 
                                double*, const int&, double*, 
                                double*,double*,int*,const int&);
extern "C" void F77NAME(sgetrf)(const int& m, const int& n, float* A, const int& lda,
                       int* ipiv, int& info);

extern "C" void F77NAME(sgetrs)(const char* trans, const int& n, const int& nrhs, 
                       float * A, const int& lda, int* ipiv, float* b, 
                       const int& ldb, int& info);
#endif
class LaFullDirectSolver : public LinearSolver
{
public:
  LaFullDirectSolver();
  LaFullDirectSolver(Mat& Min);
  virtual ~LaFullDirectSolver();
  bool prepare();
  bool solve(const Vec& bIn, Vec& xIn);
  void attachMat(Mat& Min);
  void printMatrices(const char* filename);
  void calculateCondition(DataCollector& d);
private:
  bool PRINT_MATRICES,CALCULATE_CONDITION;
  std::ofstream matOut;
  int neq;
  int LEAD_DIM_STORAGE,errorFlag;
  real* arrayptr,*rwork;\
  int* pivotptr,*iwork;
  AttacheVec x;
  Mat* M;
  DataCollector* data;
};
}//Daetk
#endif
