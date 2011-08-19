#ifndef THOMAS_GLADWELL_INTEGRATOR_H
#define THOMAS_GLADWELL_INTEGRATOR_H

#include "Definitions.h"
#include "Vec.h"

#include "Integrator.h"
#include "DaeDefinition.h"
#include "DataCollector.h"
#include "NonlinearSolver.h"
#include "Newton.h"
#include "TimeIntegrationErrorController.h"
#include "VectorNorm.h"

namespace Daetk
{

template <class NLSOLVER>
class ThomasGladwellIntegratorBase : public Integrator
{
  /***********************************************************************

    The purpose of this class is to implement the core of the
    Thomas-Gladwell time integration scheme presented in
    \cite{Kavetski_Binning_etal_01, Kavetski_Binning_etal_02} for a
    system of implicit ODEs.

    
    \mat{M}(\vec{y})\cdot \od{\vec{y}}{t} + 
        \mat{K}(\vec{y})\cdot \vec{y}    = \vec{F}(\vec{y})          (1)


    The basic approach will be to solve 

       [\mat{M}^{n+1} + \Delta t\mat{K}^{n+1}]\dot{\vec{y}}^{n+1}  
          &=& -\mat{K}^{n+1}\cdot \vec{y}^n + \vec{F}^{n+1}          (2)

       \vec{y}^{n+1}_{(1)} =
                         \vec{y}^{n} + \Delta t  \dot{\vec{y}}^{n+1} (3)
       \vec{y}^{n+1}_{(2)} = \vec{y}^{n} + 
         \frac{\Delta t}{2}(\dot{\vec{y}}^{n} + \dot{\vec{y}}^{n+1}) (4)

    with the coeficients $\mat{X}^{n+1}$ evaluated at $t^{n+1}$ and
    $\vec{y}^{n+1}_{(1)}$ ($\mat{X}=\mat{M},\mat{K},\mat{F}$).

    Here $\vec{y}^{n+1}_{(1)}$ is just the backward Euler solution to
    (1), while $\vec{y}^{n+1}_{(2)}$ is a second order approximation
    from the Thomas-Gladwell formulation.

    I will follow the suggestion in \cite{Kavetski_Binning_etal_02}
    and setup the code to solve for the backward Euler solution
    $\vec{y}^{n+1}_{(1)}$ and then calculate \dot{\vec{y}}^{n+1} and
    \vec{y}^{n+1}_{(2)} from (3) and (4).

    I'll start by trying to reuse the Integrator and DaeDefinition
    interfaces from Chris' daetk library. The Integrator interface is
    pretty mininimal. I also don't think I'll need all of the
    DaeDefinition interface, but one never knows. Hopefully, I'll only
    need to use

   
  **********************************************************************/
public:

  ThomasGladwellIntegratorBase(DaeDefinition* daeIn, 
			       DataCollector* dataIn,
			       NLSOLVER* nlSolverIn,
			       VectorNorm* errNormIn,
			       TimeIntegrationErrorController* errContrIn);

  virtual ~ThomasGladwellIntegratorBase();

  virtual bool ok();

  //---------- Integrator Interface ----------
  //step to tout taking as many substeps as necessary, 
  //load soln and deriv into solution, solutionPrime
  virtual bool calculateSolution(const real& tout,Vec& solution, 
                                 Vec& solutionPrime);
  //take a single successful step to tOut
  //record t value reached in tStep
  virtual bool step(const real& tOut,real& tStep,Vec& solutionAtTStep, 
                    Vec& solutionPrime);
  virtual void reset();

protected:

  //---------- Data Members ----------
  //user-defined problem supposed to be solving
  DaeDefinition* dae;
  //collect statistics about simulation
  DataCollector* data;
  //Nonlinear solver for problem
  NLSOLVER* nlSolver;
  //type of nonlinear solver strategy
  //0 --- Full Newton (solver updates Jacobian each  time step) [default]
  //1 --- modified Newton (1 jacobian update per attempted time step)
  //  --- not implemented yet!!!
  //2 --- Noniterative solver, 1 Newton iteration per step
  int nonlinearSolverStrategy;
  //norm to use in error calculations
  VectorNorm* errNorm;

  //approximate solutions
  Vec* ynp1;       //\vec{y}^{n+1}, storage taken from user
  Vec* ynp1Prime;  //\dot{\vec{y}}^{n+1}, storage taken from user
  Vec ynp1_1;      //\vec{y}^{n+1}_{(1)}, backward Euler solution
  Vec ynp1_2;      //\vec{y}^{n+1}_{(2)}, tg solution

  real deltaTnp1;   //\Delta t^{n+1}_{j+1}

  //history
  Vec yn;           //\vec{y}^{n}
  Vec ynPrime;      //\dot{\vec{y}}^{n}
  Vec ynm1Prime;    //\dot{\vec{y}}^{n-1}
  //not strictly necessary
  Vec yn_1;           //\vec{y}^{n}_{(1)}
  Vec yn_2;           //\vec{y}^{n}_{(2)}
  
  real tn;          //t^{n}
  real deltaTn;     //\Delta t^{n}
  real deltaTnm1;   //\Delta t^{n-1}



  //data necessary for building approximation
  bool firstStep;   //is history valid?

  real timeEps;     //\eps for making sure no divide by zero
  Vec errEst;       //\vec e^{n+1}
  real locErr;      //\|\vec e^{n+1} \| 
  real locErrLast;  //\|\vec e^{n} \| or \|\vec e^{n+1,j} \|

  //put default error controll in a class to play with
  TimeIntegrationErrorController* controller;

  int maxErrorFailures;         //maximum number of trucation error failures
                                //per step
  int maxNLSolverFailures;      //maximum number of nonlinear solver failures
                                //per step
  int maxFailures;              //maximum failures of either type
  real dtRnlFail;               //how much to cut step if there's a nl failure

  //work arrays
  Vec tmpvec,corrLast;

  //---------- Basic Functions For Algorithm ----------
  //set y^{n+1,j+1} = y^{n} + \Delta t^{n+1} \dot{y}^{n} 
  //   + (\Delta t^{n+1})^2(\dot{y}^{n}-\dot{y}^{n-1}})/(2\Delta t^{n})
  virtual bool calculatePredictor(Vec& ypred);

  //set \dot{\vec y}^{n+1} = (\vec y^{n+1}_{(1)} - \vec y^{n})/dt
  virtual bool calculateYPrime(const Vec& y, const real& dt, Vec& yp);

  //calculate \vec{e}^{n+1} = 
  //  \vec{y}^{n+1}_{(1)}-\vec{y}^{n+1}_{(2)}
  //and  \|\vec {e}^{n+1}\| using max weighted error test from
  // Kavetski_Binning_etal_01 (AWR)
  //sets variables 
  // errEst, locErr
  virtual bool calculateErrorEstimate(const Vec& y1, const Vec& y2);

  
  //calculate second order Thomas-Gladwell solution following (16) 
  //from Kavetski_Binning_etal_02
  virtual bool calculateThomasGladwellSolution(const Vec& yOld,  
					       const Vec& ypOld,
					       const Vec& ypNew, 
					       const real& dt,
					       Vec& yNew);

  //look ahead algorithm to make sure \Delta t^{n+1}_{j+1} 
  //doesn't step over t^{out}
  virtual bool checkStepSize(const real& t0, const real& tout, real& dt);
  //pick \Delta t^{0}
  virtual bool calculateInitialDt(const real& t0, const real& tout);

  //---------- utility functions ----------
  //dimension everything
  virtual bool setSizes();
  //check basic data dependencies
  virtual bool dataOk();
  //check that sizes match
  virtual bool sizesOk();
  //check everything including ynp1, ynp1Prime
  virtual bool allDataOk();

};//end ThomasGladwellIntegratorBase

class ThomasGladwellIntegrator: 
public ThomasGladwellIntegratorBase<Newton>
{
  /***********************************************************************
   The default Thomas-Gladwell algorithm that is second order accurate
   and solves nonlinear systems using a full Newton solve. That is, it
   doesn't try to do a noniterative version, or a modified Newton
   solve.

  **********************************************************************/
public: 
  ThomasGladwellIntegrator(DaeDefinition* daeIn, 
			   DataCollector* dataIn,
			   Newton* nlSolverIn,
			   VectorNorm* errNormIn,
			   TimeIntegrationErrorController* errControllerIn);

  virtual ~ThomasGladwellIntegrator();

};

class ThomasGladwellNonIterativeIntegrator: 
public ThomasGladwellIntegratorBase<Newton>
{
  /***********************************************************************
   The default non-iterative Thomas-Gladwell algorithm that is second
   order accurate and solves nonlinear systems using a single Newton
   iteration.

   Should require only changing a few routines: 
      step  :  remove test for nonlinear solver failure
      pred  :  drop (\Delta t^{n+1})^2 term 
               (see eqn 13 in Kavetski etal WRR 02)
      ok    : enforce Newton solve has max its = 1
  **********************************************************************/
public: 
  ThomasGladwellNonIterativeIntegrator(DaeDefinition* daeIn, 
				       DataCollector* dataIn,
				       Newton* nlSolverIn,
				       VectorNorm* errNormIn,
				       TimeIntegrationErrorController* 
				       errControllerIn);

  virtual ~ThomasGladwellNonIterativeIntegrator();

  virtual bool ok();

  //--- main routine ---
  virtual bool step(const real& tOut,real& tStep,Vec& solutionAtTStep, 
                    Vec& solutionPrime);

protected:

  //---------- Basic Functions For Algorithm ----------
  //set y^{n+1,j+1} = y^{n} + \Delta t^{n+1} \dot{y}^{n} 
  // following (13) in Kavetski etal WRR 02
  virtual bool calculatePredictor(Vec& ypred);

};

//======================================================================
//ThomasGladwellIntegratorBase implementation
//======================================================================
template <class NL>
ThomasGladwellIntegratorBase<NL>::
ThomasGladwellIntegratorBase(DaeDefinition* daeIn, DataCollector* dataIn,
			     NL* nlSolverIn, 
			     VectorNorm* errNormIn,
			     TimeIntegrationErrorController* errControllerIn):
  dae(daeIn), data(dataIn), nlSolver(nlSolverIn),
  nonlinearSolverStrategy(0),errNorm(errNormIn),
  ynp1(0),ynp1Prime(0),
  ynp1_1(),ynp1_2(),
  deltaTnp1(-12345.),
  yn(), ynPrime(), ynm1Prime(),
  yn_1(),yn_2(),
  tn(-12345.), deltaTn(-12345.), deltaTnm1(-12345.),
  firstStep(true),
  timeEps(1.0e-10),
  errEst(),locErr(-12345.),locErrLast(-12345.0),
  controller(errControllerIn),
  maxErrorFailures(12),maxNLSolverFailures(12),
  maxFailures(12), dtRnlFail(0.25),
  tmpvec(),corrLast()
{
  assert(dataOk());
  maxErrorFailures = controller->getMaxFailures();
  maxFailures = std::max(maxErrorFailures,maxNLSolverFailures);

  bool resizeFailed = setSizes();
  assert(!resizeFailed);

}

template <class NL>
ThomasGladwellIntegratorBase<NL>::~ThomasGladwellIntegratorBase()
{

}

template <class NL>
bool ThomasGladwellIntegratorBase<NL>::setSizes()
{
  bool resizeFailed = false;
  if (!dataOk())
    {
      resizeFailed = true;
      return resizeFailed;
    }

  int neq = dae->getY0().size();
  ynp1_1.newsize(neq);
  ynp1_2.newsize(neq);
  yn.newsize(neq);
  ynPrime.newsize(neq);
  ynm1Prime.newsize(neq);
  yn_1.newsize(neq);
  yn_2.newsize(neq);
  errEst.newsize(neq);
  tmpvec.newsize(neq);
  corrLast.newsize(neq);

  return false;
}

template <class NL>
void ThomasGladwellIntegratorBase<NL>::reset()
{
  firstStep = true;
  ynp1 = 0;
  ynp1Prime = 0;

  tn = -12345.0; deltaTn = -12345.0; deltaTnm1 = -12345.0;
  deltaTnp1 = -12345.0; locErr = -12345.0; locErrLast = -12345.0;
   
  yn   = -12345.0; ynPrime = -12345.0; ynm1Prime = -12345.0;
  yn_1 = -12345.0; yn_2    = -12345.0;
  errEst = -12345.0;
  tmpvec = 0.0; corrLast = 0.0;

}
//======================================================================
//diagnostic routies and sanity checks
//======================================================================

template <class NL>
bool ThomasGladwellIntegratorBase<NL>::dataOk()
{
  bool isOk = true;

  isOk = isOk && dae;
  isOk = isOk && data;
  isOk = isOk && nlSolver;
  isOk = isOk && errNorm;
  isOk = isOk && controller;

  return isOk;
}

template <class NL>
bool ThomasGladwellIntegratorBase<NL>::sizesOk()
{
  bool isOk = true;
  isOk = isOk && dataOk();
  if (isOk)
    {
      int neq = dae->getY0().size();
      isOk = isOk && ynp1_1.size()    == neq;
      isOk = isOk && ynp1_2.size()    == neq;
      isOk = isOk && yn.size()        == neq;
      isOk = isOk && yn_1.size()      == neq;
      isOk = isOk && yn_2.size()      == neq;
      isOk = isOk && ynPrime.size()   == neq;
      isOk = isOk && ynm1Prime.size() == neq;
      isOk = isOk && errEst.size()    == neq;
      isOk = isOk && tmpvec.size()   == neq;
      isOk = isOk && corrLast.size()  == neq;
    }

  return isOk;
}

template <class NL>
bool ThomasGladwellIntegratorBase<NL>::allDataOk()
{
  bool isOk = dataOk();
 
  isOk = isOk && ynp1;
  isOk = isOk && ynp1Prime;

  return isOk;
}

template <class NL>
bool ThomasGladwellIntegratorBase<NL>::ok()
{
  bool isOk = dataOk();
  
  isOk = isOk && sizesOk();

  return isOk;
}

//======================================================================
//main routines
//======================================================================
template <class NL>
bool ThomasGladwellIntegratorBase<NL>::calculateSolution(const real& tout,
						 Vec& solutionAtTout, 
						 Vec& solutionPrime)
{
  /***********************************************************************
 
    integrate to tOut in as many steps as necessary
    up to predefined number of allowed failures 
   
  **********************************************************************/
  data->startUserStep();

  assert(ok());
  bool stepFailed = false;
  real tCurrent;


  //initialize everything if it's the first step
  if (firstStep)
    {
      yn      = dae->getY0();
      tn      = dae->getT0();
      ynPrime = dae->getY0prime();

      calculateInitialDt(tn,tout);

#ifdef DEBUG_TG
      std::cout<<"TG calcSoln("<<tout<<") dt0 = "<<deltaTnp1<<std::endl;
#endif
    }
  else
    {
      //otherwise assume that all the history is set and
      //that error estimate from previous solution is ok to
      //get step size
      //probably need to make this more robust
      bool lastStepOk = true;
      deltaTnp1 = controller->estimateStepSize(locErrLast,deltaTn,lastStepOk);
#ifdef DEBUG_TG
      std::cout<<"TG calcSoln("<<tout<<") tn = "<<tn<<" dt = "<<deltaTnp1
	       <<" lastStepOk= "<<lastStepOk<<std::endl;
#endif
      assert(lastStepOk);
    }
  tCurrent = tn;

  while (tCurrent < tout)
    {
      stepFailed = step(tout,tCurrent,solutionAtTout,solutionPrime);
      if (stepFailed)
	{
	  std::cerr<<"ThomasGladwell calcSoln ("<<tn<<" --> "<<tCurrent
		   <<") failed "<<std::endl;
          data->includeSolution(-1,solutionAtTout);
	  return stepFailed;
	}
    }

  data->endUserStep();
  data->includeSolution(tout,solutionAtTout);
  
  return false;
}


template <class NL>
bool ThomasGladwellIntegratorBase<NL>::step(const real& tOut,
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
  int nlSolverFailures= 0;
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
      std::cout<<"TG step("<<tOut<<") tn = "<<tn<<" dt = "<<deltaTnp1<<std::endl;
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
      assert(nonlinearSolverStrategy == 0); //Full Newton
      dae->computeOwnJacobianDelta = true;
      //solve nonlinear problem for new time step
      corrLast = 0.0;
      //DaeDefinition correction step should take care of
      //ypDaeDef for me
      bool nlSolveFailed = nlSolver->solve(corrLast,*dae);
  
      if (nlSolveFailed)
	{
	  //pick dt using nonlinear solver failure heuristic
	  std::cerr<<"TG  corrector failed "<<std::endl;
	  nlSolverFailures++;
	  stepsFailed++;

	  deltaTnp1 *= dtRnlFail;
	  //mwf let time step controller know last step history no good
	  controller->invalidateLastEstimate();
	  
	}
      else 
	{
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

	  bool newDtFailed = checkStepSize(tn,tOut,dtNew);


	  if (newDtFailed)
	    return newDtFailed;

#ifdef DEBUG_TG
	  std::cout<<"TG step("<<tOut<<") nlsolve to = "<<tn+deltaTnp1
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

	}//end nonlinear solver failure succeeded
    }//end while number of failures 
  //must have failed if got to this point
  data->stepTaken(2,deltaTnp1,tn);
  return true;

}
//======================================================================
//basic algorithm functions 
//======================================================================
template <class NL>
bool ThomasGladwellIntegratorBase<NL>::calculatePredictor(Vec& ypred)
{
  bool  predFailed = false;
  if (!firstStep)
    {
      //supposed to be second order extrapolation 
      //y^{n+1,1} = y^{n} + \Delta t^{n+1} \dot{y}^{n} + 
      //  (\Delta t^{n+1})^2(\dot{y}^{n}-\dot{y}^{n-1})/(2\Delta t^{n})

      //y^{n+1,1} = y^{n}
      ypred = yn;

      //y^{n+1,1} += \Delta t^{n+1} \dot{y}^{n}
      axpy(deltaTnp1,ynPrime,ypred);

      //mwf debug turn off quadratic predictor
      //(\Delta t^{n+1})^2/(2\Delta t^{n})
      real dtFac = 0.5*deltaTnp1*deltaTnp1/deltaTn;

      //\dot{y}^{n}-\dot{y}^{n-1}
      tmpvec = ynPrime;
      axpy(-1.0,ynm1Prime,tmpvec);

      //y^{n+1,1} += last term
      axpy(dtFac,tmpvec,ypred);

    }
  else
    {
      //use forward euler for predictor
      //y^{1,1} = y^{0} + \Delta t^{1}*\dot{y}^{0}
      ypred = yn;
      //y^{n+1,1} += \Delta t^{n+1} \dot{y}^{n}
      axpy(deltaTnp1,ynPrime,ypred);
    }

  return predFailed;
}

template <class NL>
bool 
ThomasGladwellIntegratorBase<NL>::calculateYPrime(const Vec& y, 
						  const real& dt, 
						  Vec& yp)
{
  /***********************************************************************
     set \dot{\vec y}^{n+1} = \vec y^{n+1}_{(1)} - \vec y^{n}

  **********************************************************************/
  yp = y;
  axpy(-1.0,yn,yp);
  scal(1.0/dt,yp);

  return false;
}

template <class NL>
bool ThomasGladwellIntegratorBase<NL>::calculateErrorEstimate(const Vec& y1,
							      const Vec& y2)
{
  /***************************************************

  calculate \vec{e}^{n+1} = 
    \vec{y}^{n+1}_{(1)}-\vec{y}^{n+1}_{(2)}
  and  \|\vec {e}^{n+1}\| 
  Kavetski_Binning_etal_01 (AWR) use weighted max error test
  Now allow generic error norm by default
  sets variables 
    errEst, locErr

  **************************************************/
  //originally always use use max weighted error test from
  //Kavetski_Binning_etal_01 (AWR)
  //now allow more general error norms

  errEst = y1;
  axpy(-1.0,y2,errEst);

  assert(errNorm);
  tmpvec = y1;
  for (int i=0; i < y1.size(); i++)
    {
      if (std::fabs(yn[i]) > std::fabs(y1[i]))
	tmpvec[i] = yn[i];
    }
  errNorm->setWeight(tmpvec);

  locErr = (*errNorm)(errEst);

  bool errEstFailed = locErr < 0.0;
  return errEstFailed;
}

template <class NL>
bool ThomasGladwellIntegratorBase<NL>::calculateInitialDt(const real& t0,
							  const real& tout)
{

   //follow approach from Shampine_94 pg 379
   //assumes ynPrime is set
   assert(errNorm);
   errNorm->setWeight(yn);

   real normY0Prime = (*errNorm)(ynPrime);
   real denom    = std::max(timeEps,normY0Prime);
   int errOrder  = controller->getOrder();
   real errOrdInv= 1.0/(errOrder+1.0);
   real timeTol = controller->getTolerance();
   real numer = std::pow(timeTol,errOrdInv);
   real dtMin = numer/denom;

  deltaTnp1 = std::min(tout-t0,dtMin);
  bool calcFailed = false;
  calcFailed = t0 + deltaTnp1 > tout;
  
  return calcFailed;
}

template <class NL>
bool 
ThomasGladwellIntegratorBase<NL>::
calculateThomasGladwellSolution(const Vec& yOld,  const Vec& ypOld,
				const Vec& ypNew, const real& dt,
				Vec& yNew)
{
  //calculate second order Thomas-Gladwell solution following (16) 
  //from Kavetski_Binning_etal_02

  tmpvec = ypOld;
  axpy(1.0,ypNew,tmpvec);

  yNew = yOld;
  axpy(0.5*dt,tmpvec,yNew);

  return false;
}

template <class NL>
bool ThomasGladwellIntegratorBase<NL>::checkStepSize(const real& t0, 
						     const real& tout,
						     real& dt)
{
  //assuming have an estimate for \delta t^{n+1} in 
  //dt then check to make sure not going over requested 
  //time level
  
  //follow algorithm 2.7 in Kavetsk_Binning_etal_01 sec 2.7
  if (t0 + dt >= tout)
    {
      dt = tout-t0;
      controller->invalidateLastEstimate();
    }
  else if (t0 + 2.0*dt >= tout)
    {
      dt = 0.5*(tout-t0);
      controller->invalidateLastEstimate();
    }
  bool checkFailed = false;
  checkFailed = t0 + dt > tout || (dt <= 1.0e-14 && t0 + dt < tout);

  return checkFailed;
}


}//Daetk

#endif
