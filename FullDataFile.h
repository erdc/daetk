#ifndef FULLDATAFILE_H
#define FULLDATAFILE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include "Definitions.h"
#include "Utilities.h"
#include "TimingDataFile.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
class FullDataFile : public TimingDataFile
{
 public:
  FullDataFile(const real& t0=0.0,const char* filename="data.txt");
  virtual ~FullDataFile();  
  void functionEvaluation();
  void jacobianEvaluation();
  void linearSolverIteration();
  void setLinearSolverIterations(int i);
  void nonlinearSolverIteration();
  void lineSearch();
  virtual void stepTaken(int k,real h, real tn, real errorEstimate=0);
  virtual void includeSolution(const real& t,const Vec& y);
  void errorFailure();
  void nonlinearSolverFailure();
  void linearSolverFailure();
  void conditionNumber(real& cond);
  int getLinearSolverIterations();
  int getGlobalLinearSolverIterations();
  int getTotalLinearSolverIterations();
  int getNonlinearSolverIterations();
  int getGlobalNonlinearSolverIterations();
  int getTotalNonlinearSolverIterations();
  int getGlobalFunctionEvaluations ();
  int getTotalFunctionEvaluations ();
  //mwf added for convenience
  int getGlobalJacobianEvaluations();
  int getGlobalStepsTaken();
  int getGlobalStepsFailed();
  int getGlobalLinesearches();
  int getGlobalLinearSolverFailures();
  int getGlobalNonlinearSolverFailures();
  //mwf end additions
  void omitSolution();
  void reset();
protected:
  bool INCLUDE_SOLUTION;
  real tlast, condition, avgCondition;
  int functionEvaluations,
    jacobianEvaluations,
    linearSolverIterations,
    nonlinearSolverIterations,
    lineSearches,
    stepsTaken,
    orderOneStepsTaken,
    orderTwoStepsTaken,
    orderThreeStepsTaken,
    orderFourStepsTaken,
    orderFiveStepsTaken,
    errorFailures,
    nonlinearSolverFailures,
    linearSolverFailures,
    totalJacobianEvaluations,
    totalFunctionEvaluations,
    totalStepsTaken,
    totalStepsFailed,
    totalNonlinearSolverIterations,
    totalNonlinearSolverFailures,
    totalLinearSolverIterations,
    totalLinearSolverFailures,
    totalLinesearches,
    totalErrorFailures,
    globalJacobianEvaluations,
    globalFunctionEvaluations,
    globalStepsTaken,
    globalStepsFailed,
    globalNonlinearSolverIterations,
    globalNonlinearSolverFailures,
    globalLinearSolverIterations,
    globalLinearSolverFailures,
    globalLinesearches,
    globalErrorFailures,
    globalOrderOneStepsTaken,
    globalOrderTwoStepsTaken,
    globalOrderThreeStepsTaken,
    globalOrderFourStepsTaken,
    globalOrderFiveStepsTaken;
};
}//Daetk
#endif
