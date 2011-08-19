#include "TexDataFile.h"

namespace Daetk 
{

using std::setw;
using std::setprecision;
using std::endl;
  using std::flush;

  TexDataFile::TexDataFile(const real& t0,const char* filename):
    FullDataFile(t0),
    texSummaryOut(filename),
    matlabOut("data.m")
  {
    texSummaryOut<<"\\documentclass[10pt]{article}"<<endl
                 <<"\\begin{document}"<<endl
                 <<"\\begin{tabular}{cccccccccccccccccccc}"<<endl
                 <<"& $t_{m}$  & $t_{p}$ & Jeval  & Feval & Steps & F & EF & NLI & NLF & LI & LF & $\\kappa$ & LS & BDF1 & BDF2 & BDF3 & BDF4 & BDF5 & Clock \\\\"<<endl
                 <<"\\hline"<<flush;
    
    matlabOut<<"r1=["<<flush;
  }
    
TexDataFile::~TexDataFile()
{
  texSummaryOut<<"\\\\"<<endl
               <<"\\end{tabular}"<<endl
               <<"\\end{document}"<<endl;
  matlabOut<<"];"<<endl;
}

  void TexDataFile::stepTaken(int k, real h, real tn, real err)
  {
    matlabOut<<k<<'\t'<<h<<'\t'<<tn<<'\t'<<errorFailures<<'\t'<<nonlinearSolverFailures<<'\t'<<linearSolverFailures<<'\t'<<functionEvaluations<<'\t'<<jacobianEvaluations<<'\t'<<nonlinearSolverIterations<<'\t'<<linearSolverIterations<<'\t'<<lineSearches<<'\t'<<condition<<'\t'<<err<<";"<<endl;
    FullDataFile::stepTaken(k,h,tn,err);
  }

  void TexDataFile::includeSolution(const real& tout,const Vec& solution)
  {
    static int rNumber=1;
    runTime = clock.elapsed();
    clock.reset();
    if (totalStepsTaken)
      avgCondition/=totalStepsTaken;
    else
      avgCondition=0.0;

    texSummaryOut<<"\\\\"<<endl
                 <<"r"<<rNumber
                 <<" & "<<tlast
                 <<" & "<<tout
                 <<" & "<<totalJacobianEvaluations
                 <<" & "<<totalFunctionEvaluations
                 <<" & "<<totalStepsTaken
                 <<" & "<<totalStepsFailed
                 <<" & "<<totalErrorFailures
                 <<" & "<<totalNonlinearSolverIterations
                 <<" & "<<totalNonlinearSolverFailures
                 <<" & "<<totalLinearSolverIterations
                 <<" & "<<totalLinearSolverFailures
                 <<" & "<<avgCondition
                 <<" & "<<totalLinesearches
                 <<" & "<<orderOneStepsTaken
                 <<" & "<<orderTwoStepsTaken
                 <<" & "<<orderThreeStepsTaken
                 <<" & "<<orderFourStepsTaken
                 <<" & "<<orderFiveStepsTaken
                 <<" & "<<runTime<<flush;
    matlabOut<<"];"<<endl
             <<"r"<<rNumber+1<<"=["<<flush;
    FullDataFile::includeSolution(tout,solution);
    rNumber++;
  }


}//Daetk
