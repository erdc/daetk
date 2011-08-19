#ifndef DAFULLDIRECTSOLVER_H
#define DAFULLDIRECTSOLVER_H

#include "Definitions.h"
#include "LinearSolver.h"
#include "Mat.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
/// 
class FullDirectSolver : public LinearSolver 
/** A class for solving the inverse Jacobian vector product arising in 
  the solution of nonlinear equations */ 
{ 
public:
  ///
  FullDirectSolver(Mat& M);
  virtual ~FullDirectSolver(); 
  ///
  bool prepare(); 
  ///
  bool solve(const Vec& bIn,Vec& xIn);
private:
  int neq;
  AttacheVec   x,b;
  Vec scaleFactor;
  int* permutationVector;
  Mat* M;
};
}//Daetk
#endif
