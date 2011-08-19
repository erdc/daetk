#include "NewtonOneIt.h"

namespace Daetk 
{

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::max;
using std::min;



NewtonOneIt::NewtonOneIt():
  Newton()
{
  Tracer tr("NewtonOneIt::Newton()");
}

NewtonOneIt::NewtonOneIt(LinearSolver& linearSolverIn, 
			 Jacobian& jacIn,VectorNorm& normIn,
			 DataCollector& dataIn, int neq, real lTol, 
			 real nlTol, int maxit, 
			 LineSearchType lsType):
  Newton(linearSolverIn, jacIn,normIn,
	 dataIn,neq,lTol,nlTol,maxit,lsType)
{
  Tracer tr("NewtonOneIt::NewtonOneIt(LinearSolver& linearSolverIn, Jacobian& jacIn,VectorNorm& W,DataCollector& dataIn, int neq, real lTol,real nlTol, int maxit,LineSearchType lsType)");
}

NewtonOneIt::~NewtonOneIt()
{
  Tracer tr("NewtonOneIt::~NewtonOneIt()");
}


bool NewtonOneIt::solve(Vec& correction,VectorFunction& F)
{
  if (LOG_STEPS)
    data->startUserStep();
  bool evalFailed,linearSolverFailed;
  residual=F.value(evalFailed);
  if (evalFailed)
    {
      cerr<<"Initial guess in Newton is out of range, exiting"<<endl;
      exit(1);
    }

  weightedNorm->setWeight(F.argument());
  normOfInitialGuess =(*weightedNorm)(F.argument());
  
  if (TEST_RESIDUAL)
    {
      r0 = nrm2(residual);
    }

  iterations=0;
  correction=0.0;

  if (LOG_STEPS)
    data->stepTaken(iterations,lin_it,0.0,r0);
  
  roundOffTolerance = 100.0 * normOfInitialGuess * MACHINE_EPSILON;
  //one iteration only
  iterations = 1;
  data->nonlinearSolverIteration();
  lin_it=0;//line searches
  p=0;
  
  tmpArg=F.argument();
  tmpF=residual;
  
  data->jacobianEvaluation();
  F.computeDeltaForJacobian();

  bool jacEvalFailed = jac->evaluate(tmpArg,tmpF);
  if (jacEvalFailed || evalFailed)
    {
      cerr<<"jacobian eval failed in NewtonOneIt"<<endl;
      return jacEvalFailed;
    }

  linearSolverFailed = linearSolver->prepare();
  if (linearSolverFailed)
    {
      if (LOG_STEPS)
	{
	  data->endUserStep();
	  data->linearSolverFailure();
	  FnNew = nrm2(residual); //old residual
	  data->stepTaken(iterations,lin_it,0.0,FnNew);
	}
      else
	data->linearSolverFailure();
      
      cerr<<"linearSolver->prepare failed in Newton"<<endl;
      return linearSolverFailed;
    }
      
  linearSolverFailed=linearSolver->solve(residual,p);
  if (linearSolverFailed)
    {
      if (LOG_STEPS)
	{
	  data->endUserStep();
	  data->linearSolverFailure();
	  FnNew = nrm2(residual); //old residual
	  data->stepTaken(iterations,lin_it,0.0,FnNew);
	}
      else
	data->linearSolverFailure();
      
      cerr<<"linearSolver->solve failed in NewtonOneIt"<<endl;
      return linearSolverFailed;
    }
  p0=p;
  argPrev = F.argument();
  funPrev = residual;
  normOfLastCorrection=normOfCorrection;

  //correctArgument may change pLS so we need to save p in case we do a linesearch
  pLS=p;

  F.correctArgument(pLS);

  normOfCorrection=(*weightedNorm)(pLS);

  //skip this because not needed
  //residual = F.value(evalFailed);

  if (evalFailed)
    {
      cerr <<"evalFailed in NewtonOneIt iteration, quitting"<<endl;
      if (LOG_STEPS)
	{
	  data->endUserStep();
	  data->nonlinearSolverFailure();
	  //don't do this to save an iteration FnNew=nrm2(F.value(evalError));
	  FnNew = nrm2(residual); //old residual
	  data->stepTaken(iterations,lin_it,0.0,FnNew);
	}
      return evalFailed;
    }

  p=pLS;
#ifndef USE_BLAS
  correction-=p;
#else
  axpy(-1.0,p,correction);
#endif

  if (LOG_STEPS)
    {
      //don't do this to save an iteration FnNew=nrm2(F.value(evalError));
      FnNew = nrm2(residual); //old residual
      data->stepTaken(iterations,lin_it,0.0,FnNew);
    }

  return false;
}


}//Daetk
