#ifndef LABANDEDDIRECTSOLVER_H
#define LABANDEDDIRECTSOLVER_H

#include <fstream>
#include "Definitions.h"
#include "BandColMat.h"
#include "LinearSolver.h"
#include "DataCollector.h"
#include "Vec.h"
#include "VecBlas.h"
#include "IntVec.h"

namespace Daetk 
{

#ifdef CRAYCC
#include <fortran.h>
extern "C" void F77NAME(DGBTRF)(const int&,const int&,const int&,const int&,
                       double*,const int&,int*,int&);

extern "C" void F77NAME(DGBTRS)(_fcd,const int&,const int&,const int&,
                       const int&,double*,const int&,int*,double*,
                       const int&,int&);

extern "C" void F77NAME(DGBCON)(_fcd,const int&, const int&, const int&, 
                                double*, const int&, int*, double*, 
                                double*,double*,int*,const int&);

extern "C" void F77NAME(SGBTRF)(const int&,const int&,const int&,const int&,
                       float*,const int&,int*,int&);

extern "C" void F77NAME(SGBTRS)(_fcd,const int&,const int&,const int&,
                       const int&,float*,const int&,int*,float*,
                       const int&,int&);
#else
extern "C" void F77NAME(dgbtrf)(const int&,const int&,const int&,const int&,
                       double*,const int&,int*,int&);

extern "C" void F77NAME(dgbtrs)(const char*,const int&,const int&,const int&,
                       const int&,double*,const int&,int*,double*,
                       const int&,int&);

extern "C" void F77NAME(dgbcon)(const char*,const int&, const int&, const int&, 
                                double*, const int&, int*, double*, 
                                double*,double*,int*,const int&);

extern "C" void F77NAME(sgbtrf)(const int&,const int&,const int&,const int&,
                       float*,const int&,int*,int&);

extern "C" void F77NAME(sgbtrs)(const char*,const int&,const int&,const int&,
                       const int&,float*,const int&,int*,float*,
                       const int&,int&);
#endif

class LaBandedDirectSolver : public LinearSolver
{ 
public:
  LaBandedDirectSolver(); 
  LaBandedDirectSolver(BandColMat& BCMin);
  virtual ~LaBandedDirectSolver();
  void attachMat(BandColMat& BCMin);
  bool prepare(); 
  bool solve(const Vec& bIn,Vec& xIn);
  void printMatrices(const char* filename);
  void calculateCondition(DataCollector& d);
private:
  bool PRINT_MATRICES,CALCULATE_CONDITION;
  std::ofstream matOut;
  AttacheVec x;
  int neq,kl,ku;
  int LEAD_DIM_STORAGE,errorFlag;
  real* arrayptr,*rwork;
  int* pivotptr,*iwork;
  BandColMat* BCM;
  DataCollector* data;
};
}//Daetk 
#endif
