#include "DataCollector.h"

namespace Daetk 
{

DataCollector::DataCollector():
  analyticSolution(0),
  analyticSolutionAvailable(false)
{
  Tracer("DataCollector::DataCollector()");
}

DataCollector::~DataCollector()
{
  Tracer("DataCollector::~DataCollector()");
}

void DataCollector::setAnalyticSolution(bool (*aS)(const real& t,Vec& solution))
{
  analyticSolutionAvailable = true;
  analyticSolution = aS;
}


void DataCollector::startUserStep(){}
void DataCollector::endUserStep(){}
void DataCollector::functionEvaluation(){}
void DataCollector::jacobianEvaluation(){}
void DataCollector::linearSolverIteration(){}
void DataCollector::setLinearSolverIterations(int i){}
int  DataCollector::getLinearSolverIterations(){return 0;}
void DataCollector::nonlinearSolverIteration(){}
void DataCollector::lineSearch(){}
void DataCollector::stepTaken(int,real, real, real){}
void DataCollector::stepTaken(int,real, real, real, real, real, real){}
void DataCollector::includeSolution(const real&,const Vec&){}
void DataCollector::errorFailure(){}
void DataCollector::nonlinearSolverFailure(){}
void DataCollector::linearSolverFailure(){}
void DataCollector::conditionNumber(real& k){}
void DataCollector::reset(){}

}//Daetk
