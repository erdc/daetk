#include "FullDataFileForRichExtrap.h"

namespace Daetk 
{

using std::setw;
using std::setprecision;
using std::endl;

FullDataFileForRichExtrap::FullDataFileForRichExtrap(const real& t0,const char* filename):
  TimingDataFile(filename),
  INCLUDE_SOLUTION(true),
  tlast(t0),
  condition(0),
  avgCondition(0),
  functionEvaluations(0),
  jacobianEvaluations(0),
  linearSolverIterations(0),
  nonlinearSolverIterations(0),
  lineSearches(0),
  orderOneStepsTaken(0),
  orderTwoStepsTaken(0),
  orderThreeStepsTaken(0),
  orderFourStepsTaken(0),
  orderFiveStepsTaken(0),
  errorFailures(0),
  nonlinearSolverFailures(0),
  linearSolverFailures(0),
  totalJacobianEvaluations(0),
  totalFunctionEvaluations(0),
  totalStepsTaken(0),
  totalStepsFailed(0),
  totalNonlinearSolverIterations(0),
  totalNonlinearSolverFailures(0),
  totalLinearSolverIterations(0),
  totalLinearSolverFailures(0),
  totalLinesearches(0),
  totalErrorFailures(0),
  globalJacobianEvaluations(0),
  globalFunctionEvaluations(0),
  globalStepsTaken(0),
  globalStepsFailed(0),
  globalNonlinearSolverIterations(0),
  globalNonlinearSolverFailures(0),
  globalLinearSolverIterations(0),
  globalLinearSolverFailures(0),
  globalLinesearches(0),
  globalErrorFailures(0),
  globalOrderOneStepsTaken(0),
  globalOrderTwoStepsTaken(0),
  globalOrderThreeStepsTaken(0),
  globalOrderFourStepsTaken(0),
  globalOrderFiveStepsTaken(0)
{
  Tracer tr("FullDataFileForRichExtrap::FullDataFileForRichExtrap(const real& t0,const int Neq)");
  fout<<"======================================================================================================================================"<<endl
      <<" k |  Step Size   |      tn      |errf|nlsf|ilsf|func|jacs|NLit|L_IT|linS|     cond.    |    error     |     cpu      |  errorRichEx |     extra    |"<<endl;
}

FullDataFileForRichExtrap::~FullDataFileForRichExtrap()
{
  Tracer tr("FullDataFileForRichExtrap::~FullDataFileForRichExtrap()");
  fout<<endl
      <<"Data for entire run"<<endl
      <<"------------------------------"<<endl
      <<"Jacobians evaluated:--------"<<globalJacobianEvaluations<<endl
      <<"Function Evaluations:-------"<<globalFunctionEvaluations<<endl
      <<"Steps Taken:----------------"<<globalStepsTaken<<endl
      <<"Steps Failed:---------------"<<globalStepsFailed<<endl
      <<"Error Failures:-------------"<<globalErrorFailures<<endl
      <<"Nonlinear Solver Iterations:"<<globalNonlinearSolverIterations<<endl
      <<"Nonlinear Solver Failures:--"<<globalNonlinearSolverFailures<<endl
      <<"Linear Solver Iterations:---"<<globalLinearSolverIterations<<endl
      <<"Linear Solver Failures:-----"<<globalLinearSolverFailures<<endl
      <<"Line searches:--------------"<<globalLinesearches<<endl
      <<"Order 1 steps:--------------"<<globalOrderOneStepsTaken<<endl
      <<"Order 2 steps:--------------"<<globalOrderTwoStepsTaken<<endl
      <<"Order 3 steps:--------------"<<globalOrderThreeStepsTaken<<endl
      <<"Order 4 steps:--------------"<<globalOrderFourStepsTaken<<endl
      <<"Order 5 steps:--------------"<<globalOrderFiveStepsTaken<<endl
      <<endl;

}

void  FullDataFileForRichExtrap::stepTaken(int k,real h,real tn, real err)
{ stepTaken(k,h,tn,err,-12345.,-12345.); }

void  FullDataFileForRichExtrap::stepTaken(int k,real h,real tn, real err, real cpu, real err2, real extra)
{
  avgCondition += condition;
  totalStepsTaken++;
  totalFunctionEvaluations+=functionEvaluations;
  totalJacobianEvaluations+=jacobianEvaluations;
  totalNonlinearSolverIterations+=nonlinearSolverIterations;
  totalNonlinearSolverFailures+=nonlinearSolverFailures;
  totalLinesearches+=lineSearches;
  totalLinearSolverIterations+=linearSolverIterations;
  totalLinearSolverFailures+=linearSolverFailures;
  totalErrorFailures+=errorFailures;
  totalStepsFailed+=(errorFailures+nonlinearSolverFailures);

  switch(k)
    {
    case 1:
      orderOneStepsTaken++;
      break;
    case 2:
      orderTwoStepsTaken++;
      break;
    case 3:
      orderThreeStepsTaken++;
      break;
    case 4:
      orderFourStepsTaken++;
      break;
    case 5:
      orderFiveStepsTaken++;
      break;
    }

  fout<<setw(3)<<k<<setw(1)<<"|"
      <<setw(14)<<setprecision(9)<<h<<setw(1)<<"|"
      <<setw(14)<<tn<<setw(1)<<"|"
      <<setw(4)<<errorFailures<<setw(1)<<"|"
      <<setw(4)<<nonlinearSolverFailures<<setw(1)<<"|"
      <<setw(4)<<linearSolverFailures<<setw(1)<<"|"
      <<setw(4)<<functionEvaluations<<setw(1)<<"|"
      <<setw(4)<<jacobianEvaluations<<setw(1)<<"|"
      <<setw(4)<<nonlinearSolverIterations<<setw(1)<<"|"
      <<setw(4)<<linearSolverIterations<<setw(1)<<"|"
      <<setw(4)<<lineSearches<<setw(1)<<"|"
      <<setw(14)<<setprecision(9)<<condition<<setw(1)<<"|"
      <<setw(14)<<setprecision(9)<<err<<setw(1)<<"|"
      <<setw(14)<<setprecision(9)<<cpu<<setw(1)<<"|"
      <<setw(14)<<setprecision(9)<<err2<<setw(1)<<"|"
      <<setw(14)<<setprecision(9)<<extra<<setw(1)<<"|"
      <<endl;
  errorFailures=0;
  functionEvaluations=0;  
  jacobianEvaluations=0; 
  nonlinearSolverIterations=0;
  nonlinearSolverFailures=0;
  linearSolverIterations=0; 
  linearSolverFailures=0;
  lineSearches=0;
}

void FullDataFileForRichExtrap::includeSolution(const real& tout,const Vec& solution)
{    
  runTime = clock.elapsed();
  clock.reset();
  if (totalStepsTaken)
    avgCondition/=totalStepsTaken;
  else
    avgCondition=0.0;
  fout<<endl
      <<"Data over the interval ["<<tlast<<","<<tout<<"]"<<endl
      <<"------------------------------"<<endl
      <<"Jacobian Evaluations:-------"<<totalJacobianEvaluations<<endl
      <<"Function Evaluations:-------"<<totalFunctionEvaluations<<endl
      <<"Steps Taken:----------------"<<totalStepsTaken<<endl
      <<"Steps Failed:---------------"<<totalStepsFailed<<endl
      <<"Error Failures:-------------"<<totalErrorFailures<<endl
      <<"Nonlinear Solver Iterations:"<<totalNonlinearSolverIterations<<endl
      <<"Nonlinear Solver Failures:--"<<totalNonlinearSolverFailures<<endl
      <<"Linear Solver Iterations:---"<<totalLinearSolverIterations<<endl
      <<"Linear Solver Failures:-----"<<totalLinearSolverFailures<<endl
      <<"Average Condition Number:---"<<avgCondition<<endl
      <<"Line Searches:--------------"<<totalLinesearches<<endl
      <<"Order 1 steps:--------------"<<orderOneStepsTaken<<endl
      <<"Order 2 steps:--------------"<<orderTwoStepsTaken<<endl
      <<"Order 3 steps:--------------"<<orderThreeStepsTaken<<endl
      <<"Order 4 steps:--------------"<<orderFourStepsTaken<<endl
      <<"Order 5 steps:--------------"<<orderFiveStepsTaken<<endl
      <<"Wall Clock Time:------------"<<runTime<<endl
      <<endl;
  int neq = solution.dim();
  if (analyticSolutionAvailable)
    {
      real rnorm;
      Vec r(solution.dim());
      analyticSolution(tout,r);
      fout<<"Solution at t="<<tout<<"   Analytic Solution"<<endl<<endl;
      for (int i=0;i<neq;i++)
	{
	  fout<<setw(14)<<solution(i)<<setw(3)<<" "<<setw(14)<<r(i)<<endl;
	}
      rnorm = nrm2(r);
      r-=solution;
      fout<<endl<<"Relative error= "<<(norm(r)/rnorm)<<endl<<endl;
    }
  else
    {
      if (INCLUDE_SOLUTION)
        {
          fout<<setw(1)<<"Solution at t="<<tout<<endl
              <<solution<<endl;
        }
    }
  fout<<"======================================================================================================================================"<<endl
      <<" k |  Step Size   |      tn      |errf|nlsf|ilsf|func|jacs|NLit|L_IT|linS|     cond.    |    error     |     cpu      |  errorRichEx |"<<endl;

  globalJacobianEvaluations+=totalJacobianEvaluations;
  globalFunctionEvaluations+=totalFunctionEvaluations;
  globalStepsTaken+=totalStepsTaken;
  globalNonlinearSolverIterations+=totalNonlinearSolverIterations;
  globalNonlinearSolverFailures+=totalNonlinearSolverFailures;
  globalLinearSolverIterations+=totalLinearSolverIterations;
  globalLinearSolverFailures+=totalLinearSolverFailures;
  globalLinesearches+=totalLinesearches;
  globalErrorFailures+=totalErrorFailures;
  globalStepsFailed+=totalStepsFailed;
  globalOrderOneStepsTaken+=orderOneStepsTaken;
  globalOrderTwoStepsTaken+=orderTwoStepsTaken;
  globalOrderThreeStepsTaken+=orderThreeStepsTaken;
  globalOrderFourStepsTaken+=orderFourStepsTaken;
  globalOrderFiveStepsTaken+=orderFiveStepsTaken;

  avgCondition =0;
  totalJacobianEvaluations=0;
  totalFunctionEvaluations=0;
  totalStepsTaken=0;
  totalNonlinearSolverIterations=0;
  totalNonlinearSolverFailures=0;
  totalLinearSolverIterations=0;
  totalLinearSolverFailures=0;
  totalLinesearches=0;
  totalErrorFailures=0;
  totalStepsFailed=0;
  tlast=tout; 
  orderOneStepsTaken=0;
  orderTwoStepsTaken=0;
  orderThreeStepsTaken=0;
  orderFourStepsTaken=0;
  orderFiveStepsTaken=0;
}

  void FullDataFileForRichExtrap::functionEvaluation(){functionEvaluations++;}
  void FullDataFileForRichExtrap::jacobianEvaluation(){jacobianEvaluations++;}
  void FullDataFileForRichExtrap::linearSolverIteration(){linearSolverIterations++;}
  void FullDataFileForRichExtrap::setLinearSolverIterations(int i){linearSolverIterations+=i;}
  void FullDataFileForRichExtrap::nonlinearSolverIteration(){nonlinearSolverIterations++;}
  void FullDataFileForRichExtrap::lineSearch(){lineSearches++;}

  void FullDataFileForRichExtrap::errorFailure(){errorFailures++;}
  void FullDataFileForRichExtrap::nonlinearSolverFailure(){nonlinearSolverFailures++;}
  void FullDataFileForRichExtrap::linearSolverFailure(){linearSolverFailures++;}
  void FullDataFileForRichExtrap::conditionNumber(real& cond){condition = cond;}
  int FullDataFileForRichExtrap::getLinearSolverIterations(){return linearSolverIterations;}
  int FullDataFileForRichExtrap::getGlobalLinearSolverIterations(){return globalLinearSolverIterations;}
  int FullDataFileForRichExtrap::getTotalLinearSolverIterations(){return totalLinearSolverIterations;}
  int FullDataFileForRichExtrap::getNonlinearSolverIterations(){return nonlinearSolverIterations;}
  int FullDataFileForRichExtrap::getGlobalNonlinearSolverIterations(){return globalNonlinearSolverIterations;}
  int FullDataFileForRichExtrap::getTotalNonlinearSolverIterations(){return totalNonlinearSolverIterations;}
  int FullDataFileForRichExtrap::getGlobalFunctionEvaluations (){return globalFunctionEvaluations;}
  int FullDataFileForRichExtrap::getTotalFunctionEvaluations (){return totalFunctionEvaluations;}
  void FullDataFileForRichExtrap::omitSolution(){ INCLUDE_SOLUTION=false;}
  void FullDataFileForRichExtrap::reset()
    {
      globalJacobianEvaluations=0;
    globalFunctionEvaluations=0;
    globalStepsTaken=0;
    globalStepsFailed=0;
    globalNonlinearSolverIterations=0;
    globalNonlinearSolverFailures=0;
    globalLinearSolverIterations=0;
    globalLinearSolverFailures=0;
    globalLinesearches=0;
    globalErrorFailures=0;
    globalOrderOneStepsTaken=0;
    globalOrderTwoStepsTaken=0;
    globalOrderThreeStepsTaken=0;
    globalOrderFourStepsTaken=0;
    globalOrderFiveStepsTaken=0;
    }

//mwf------- added functions for convenience ---------
int FullDataFileForRichExtrap::getGlobalJacobianEvaluations()
{ return globalJacobianEvaluations; }
int FullDataFileForRichExtrap::getGlobalStepsTaken()
{ return globalStepsTaken; }
int FullDataFileForRichExtrap::getGlobalStepsFailed()
{ return globalStepsFailed; }
int FullDataFileForRichExtrap::getGlobalLinesearches()
{ return globalLinesearches; }
int FullDataFileForRichExtrap::getGlobalLinearSolverFailures()
{ return globalLinearSolverFailures; }
int FullDataFileForRichExtrap::getGlobalNonlinearSolverFailures()
{ return globalNonlinearSolverFailures; }


}//Daetk
