#include "ThomasGladwellIntegrator.h"

#include <cassert>
#include <iostream>

#include "DaeDefinition.h"
#include "DataCollector.h"

#include "DaetkPetscVecBlas.h"

using Daetk::PetscVecBlas::axpy;
using Daetk::PetscVecBlas::scal;

namespace Daetk
{

//======================================================================
//ThomasGladwellIntegrator
//======================================================================
ThomasGladwellIntegrator::
ThomasGladwellIntegrator(DaeDefinition* daeIn, DataCollector* dataIn,
			 Newton* nlSolverIn, 
			 VectorNorm* errNormIn,
			 TimeIntegrationErrorController* errContrIn):
  ThomasGladwellIntegratorBase<Newton>(daeIn,dataIn,nlSolverIn,
				       errNormIn,errContrIn)
{
  nonlinearSolverStrategy = 0;
}

ThomasGladwellIntegrator::~ThomasGladwellIntegrator()
{

}

//======================================================================
//ThomasGladwellNonIterativeIntegrator
//======================================================================
ThomasGladwellNonIterativeIntegrator::
ThomasGladwellNonIterativeIntegrator(DaeDefinition* daeIn, 
				     DataCollector* dataIn,
				     Newton* nlSolverIn, 
				     VectorNorm* errNormIn,
				     TimeIntegrationErrorController* 
				     errContrIn):
  ThomasGladwellIntegratorBase<Newton>(daeIn,dataIn,nlSolverIn,
				       errNormIn,errContrIn)
{
  assert(dataOk());
  //only allow one nonlinear iteration
  nlSolver->setMaxIterations(1);
  nlSolver->useLineSearch(0);
  nonlinearSolverStrategy = 2;
}

ThomasGladwellNonIterativeIntegrator::~ThomasGladwellNonIterativeIntegrator()
{

}

bool ThomasGladwellNonIterativeIntegrator::ok()
{
  bool isOk = ThomasGladwellIntegratorBase<Newton>::ok();
  
  isOk = isOk && (nlSolver->getMaxIterations() == 1);

  return isOk;
}

//======================================================================
//main routines
//======================================================================
bool ThomasGladwellNonIterativeIntegrator::step(const real& tOut,
						real& tStep,
						Vec& solutionAtTStep, 
						Vec& solutionPrime)
{
  /***********************************************************************
 
    try to take one successful step starting at tn and trying to
    reach tOut. Record the point reached in tStep.

    Assumes that tn has been set correctly and deltaTnp1 has also
    been chosen correctly coming into the call. 
 
    recalculates dt until it can take a successful step

    and selects next dt after a successful step.

    Quit if reach maximum allowed failures of error type or nonlinear
    solver type.
   
  **********************************************************************/
  assert(ok());
  bool stepFailed = false;

  int errorFailures  = 0;
  int stepsFailed    = 0;

  //ynp1 and ynp1Prime hold solution at new time level.
  //the first step is to solve for ynp1 as  the first order
  //backward Euler estimate

  //borrow users storage
  ynp1      = &solutionAtTStep;
  ynp1Prime = &solutionPrime;

  //assumes that yn, tn, ynPrime, deltaTnp1 have all been set
  //by calling routine
  
  while (stepsFailed < maxFailures)
    {
      //predict new solution and corresponding derivative
      stepFailed = calculatePredictor(*ynp1);
      stepFailed = calculateYPrime(*ynp1,deltaTnp1,*ynp1Prime);
#ifdef DEBUG_TG
      std::cout<<"TG NonIt step("<<tOut<<") tn = "<<tn
	       <<" dt = "<<deltaTnp1<<std::endl;
      std::cout<<"predicted value = "<<*ynp1<<std::endl;
#endif

      //pass new values on to dae
      dae->resetFunction();
      dae->yDaeDef     = *ynp1;
      dae->ypDaeDef    = *ynp1Prime;
      dae->alphaDaeDef = 1.0/deltaTnp1;
      dae->betaDaeDef  = yn;
      scal(-1.0/deltaTnp1,dae->betaDaeDef);
      dae->tDaeDef     = tn + deltaTnp1;

      //if want to play games with Jacobian evalutation would have
      //to put that here
      assert(nonlinearSolverStrategy == 2); //NonIterative Newton
      dae->computeOwnJacobianDelta = true;
      //solve nonlinear problem for new time step
      corrLast = 0.0;
      //In non-iterative version we only allow one step
      //of Newton solve 
      //DaeDefinition correction step should take care of
      //ypDaeDef for me
      nlSolver->solve(corrLast,*dae);
  
      //get ynp1, ynp1Prime back from dae
      *ynp1     = dae->yDaeDef;
      *ynp1Prime= dae->ypDaeDef;

      ynp1_1 = *ynp1;   //Backward Euler solution at 
      //assume that ynp1Prime has been set correctly by
      //DaeDefinition correction?
      calculateYPrime(*ynp1,deltaTnp1,*ynp1Prime);
      //get second order TG solution
      calculateThomasGladwellSolution(yn,ynPrime,*ynp1Prime,
				      deltaTnp1,ynp1_2);
      
      //determine if estimated truncation error is sufficiently low
      calculateErrorEstimate(ynp1_1,ynp1_2);
      bool errorOk = false;
      //keep track of last truncation error estimate
      locErrLast = locErr;
      
      real dtNew = controller->estimateStepSize(locErr,deltaTnp1,errorOk);

#ifdef DEBUG_TG
      std::cout<<"TG NonIt step("<<tOut<<") nlsolve to = "<<tn+deltaTnp1
	       <<" succeeded "<<std::endl;
      std::cout<<"locErr = "<<locErr<<" errorOk= "<<errorOk
	       <<" dtNew = "<<dtNew<<std::endl;
#endif

      if (!errorOk) //truncation error estimate too high
	{
	  //have to take another step
	  bool newDtFailed = checkStepSize(tn,tOut,dtNew);
	  if (newDtFailed)
	    {
	      std::cerr<<"TG NonIt step("<<tOut<<") nlsolve to = "
		       <<tn+deltaTnp1
		       <<" dtNew= "<<dtNew<<" failed "<<std::endl;
	      return newDtFailed;
	    }
	  errorFailures++;
	  stepsFailed++;
	  data->errorFailure();
	  deltaTnp1 = dtNew;
	}
      else //successful step
	{
	  bool newDtFailed = checkStepSize(tn+deltaTnp1,tOut,dtNew);
	  if (newDtFailed)
	    {
	      std::cerr<<"TG NonIt step("<<tOut<<") nlsolve to = "
		       <<tn+deltaTnp1
		       <<" dtNew= "<<dtNew<<" failed "<<std::endl;
	      return newDtFailed;
	    }
	  //cycle for next step
	  tn       += deltaTnp1;
	  tStep    = tn;
	  
	  deltaTnm1 = deltaTn;
	  deltaTn   = deltaTnp1;
	  deltaTnp1 = dtNew;
	  
	  ynm1Prime = ynPrime;
	  ynPrime   = *ynp1Prime;
	  
	  //use TG solution as solution for next time step
	  *ynp1     = ynp1_2;
	  yn        = ynp1_2;
	  yn_2      = ynp1_2;
	  yn_1      = ynp1_1;
	  
	  //record step
	  data->stepTaken(2,deltaTn,tn,locErr);
	  dae->stepTaken();
	  
	  firstStep = false;
	  return false;
	}//end successful step

    }//end while number of failures 
  //must have failed if got to this point
  data->stepTaken(2,deltaTnp1,tn);
  return true;

}

//======================================================================
//basic algorithm functions 
//======================================================================
bool ThomasGladwellNonIterativeIntegrator::calculatePredictor(Vec& ypred)
{
  bool  predFailed = false;
  //use different extrapolation formula (13) in Kavetski etal WRR 02
  //y^{n+1,1} = y^{n} + \Delta t^{n+1} \dot{y}^{n}
  
  //y^{n+1,1} = y^{n}
  ypred = yn;
  
  //y^{n+1,1} += \Delta t^{n+1} \dot{y}^{n}
  axpy(deltaTnp1,ynPrime,ypred);

  return predFailed;
}

}//end Daetk
