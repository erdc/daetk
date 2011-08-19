#ifndef DATACOLLECTOR_H
#define DATACOLLECTOR_H

#include "Definitions.h"
#include "Utilities.h"
#include "Vec.h"

namespace Daetk 
{
class DataCollector
{
 public:
  DataCollector();
  virtual ~DataCollector();
  void setAnalyticSolution(bool (*)(const real& t,Vec& solution));
  virtual void startUserStep();
  virtual void endUserStep();
  virtual void functionEvaluation();
  virtual void jacobianEvaluation();
  virtual void linearSolverIteration();
  virtual void setLinearSolverIterations(int i);
  virtual int  getLinearSolverIterations();
  virtual void nonlinearSolverIteration();
  virtual void lineSearch();
  virtual void stepTaken(int,real, real,real errorEstimate=0);
  //mwf for Sarah's richardson extrapolation routines
  virtual void stepTaken(int,real, real,real errorEstimate, real cpu, real errorRichExt, real extraArg=-12345.0);
  virtual void includeSolution(const real&,const Vec&);
  virtual void errorFailure();
  virtual void nonlinearSolverFailure();
  virtual void linearSolverFailure();
  virtual void conditionNumber(real& k);
  virtual void reset();
protected:
  bool (*analyticSolution)(const real& t,Vec& solution);
  bool analyticSolutionAvailable;
};
}//Daetk
#endif
