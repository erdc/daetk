#ifndef   LINEARSOLVER_H  
#define   LINEARSOLVER_H

#include "Definitions.h"
#include "Vec.h"
#include "DataCollector.h"

namespace Daetk 
{
class LinearSolver
{
 public:
  LinearSolver();
  virtual ~LinearSolver();
  virtual bool prepare()=0;
  virtual bool solve(const Vec& b,Vec& x)=0;
  virtual void checkTrueResidual(real resTol);
  virtual void printMatrices(const char* filename);
  virtual void calculateCondition(DataCollector& d);
  void solveSubSystem(int start,int end,int stride, int block=1);

  class AttacheVec
  {
  public:
    Vec v_;
    AttacheVec(LinearSolver* lsIn, int dimIn=0);
    void attachToTarget(const Vec& target);
    void restoreToTarget();
  private:
    LinearSolver* this_ls;
    Vec attache;
  };

protected:
  friend class AttacheVec;
  bool SOLVE_SUB;
  int str,blk;
  VecIndex index;
  Vec attache;
};
}//Daetk
#endif
