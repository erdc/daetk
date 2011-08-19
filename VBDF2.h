#ifndef VBDF2_INTEGRATOR_H
#define VBDF2_INTEGRATOR_H

#include "Definitions.h"
#include "Vec.h"

#include "Integrator.h"

#include "NonlinearSolver.h"
#include "Newton.h"
#include "TimeIntegrationErrorController.h"
#include "VectorNorm.h"

namespace Daetk
{
//forward declarations
class DaeDefinition;
class DataCollector;


template <class NLSOLVER>
class VBDF2Base : public Integrator
{
  /***********************************************************************

    The purpose of this class is to implement the core of a variable
    coefficient, variable step BDF 1 and BDF 2 integrator for DAE's,
    F(t,\vec y,\vec \dot{y}) = 0. The approximations are built on
    
      \vec \dot{y}^{n+1} \approx 
           \alpha^{n+1}_{k}\vec y^{n+1} + \vec \beta^n_{k}

    and solving G^{n+1}(\vec y) = 0, 
        G^{n+1} = F(t^{n+1},\vec y,\alpha^{n+1}_{k}\vec y + \vec \beta^n_{k}
    
    The integrator uses a predictor corrector formalism with either a linear
    or a quadratic predictor. Initially, I'll use a linear predictor for k=1,2

    The code will use second order unless it's the first step and the
    history is invalid for now. I will look into variable order
    calculations later.

    I'll start by trying to reuse the Integrator and DaeDefinition
    interfaces from Chris' daetk library. The Integrator interface is
    pretty mininimal. I also don't think I'll need all of the
    DaeDefinition interface, but one never knows. Hopefully, I'll only
    need to use

   
  **********************************************************************/
public:

  VBDF2Base(DaeDefinition* daeIn, 
	    DataCollector* dataIn,
	    NLSOLVER* nlSolverIn,
	    VectorNorm* errNormIn,
	    TimeIntegrationErrorController* errContrIn);

  virtual ~VBDF2Base();

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

  real deltaTnp1;   //\Delta t^{n+1}_{j+1}

  //predicted value
  Vec ynp1Pred;    //\vec y^{n+1,p}
  //history
  Vec yn;           //\vec{y}^{n}
  Vec ynPrime;      //\dot{\vec{y}}^{n}
  Vec ynm1;        //\dot{\vec{y}}^{n-1}
  
  real tn;          //t^{n}
  real deltaTn;     //\Delta t^{n}
  real deltaTnm1;   //\Delta t^{n-1}



  //data necessary for building approximation
  real alphaBDF;
  Vec betaBDF;

  bool firstStep;   //is history valid?
  int approxOrder;  //order using for current step 
  int predOrder;    //order to use for predictor

  real timeEps;     //\eps for making sure no divide by zero
  Vec errEst;       //\vec e^{n+1}
  real locErr;      //\|\vec e^{n+1} \| 

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
  //set y^{n+1,j+1} = y^{n} + \Delta t^{n+1} \dot{y}^{n}  for now 
  //
  virtual bool calculatePredictor(Vec& ypred);

  //set \alpha_k and \vec \beta^n 
  //depending on stepsize and order
  virtual bool calculateBDFCoefs();

  //set \dot{\vec y}^{n+1} = \alpha_k\vec y + \vec \beta^n 
  //depending on stepsize and order
  virtual bool calculateYPrime(const Vec& y, const real& dt, Vec& yp);


  //calculate \vec{e}^{n+1} = 
  // errEst, locErr
  virtual bool calculateErrorEstimate(const Vec& yp, const Vec& y);

  virtual bool chooseOrder();

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

};//end VBDF2Base


//======================================================================
//VBDF2IntegratorBase implementation
//======================================================================
template <class NL>
VBDF2Base<NL>::
VBDF2Base(DaeDefinition* daeIn, DataCollector* dataIn,
			     NL* nlSolverIn, 
			     VectorNorm* errNormIn,
			     TimeIntegrationErrorController* errControllerIn):
  dae(daeIn), data(dataIn), nlSolver(nlSolverIn),
  nonlinearSolverStrategy(0),errNorm(errNormIn),
  ynp1(0),ynp1Prime(0),
  deltaTnp1(-12345.),
  ynp1Pred(),
  yn(), ynPrime(), ynm1(),
  tn(-12345.), deltaTn(-12345.), deltaTnm1(-12345.),
  alphaBDF(-12345.0),betaBDF(),
  firstStep(true),
  approxOrder(1),predOrder(1),
  timeEps(1.0e-10),
  errEst(),locErr(-12345.),
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
VBDF2Base<NL>::~VBDF2Base()
{

}

template <class NL>
bool VBDF2Base<NL>::setSizes()
{
  bool resizeFailed = false;
  if (!dataOk())
    {
      resizeFailed = true;
      return resizeFailed;
    }

  int neq = dae->getY0().size();
  ynp1Pred.newsize(neq);
  betaBDF.newsize(neq);
  yn.newsize(neq);
  ynPrime.newsize(neq);
  ynm1.newsize(neq);
  errEst.newsize(neq);
  tmpvec.newsize(neq);
  corrLast.newsize(neq);

  return false;
}

template <class NL>
void VBDF2Base<NL>::reset()
{
  firstStep = true;
  ynp1 = 0;
  ynp1Prime = 0;

  tn = -12345.0; deltaTn = -12345.0; deltaTnm1 = -12345.0;
  alphaBDF  = -12345.0;
  deltaTnp1 = -12345.0; locErr = -12345.0; 
   
  ynp1Pred = -12345.0;
  yn   = -12345.0; ynPrime = -12345.0; ynm1 = -12345.0;
  betaBDF  = -12345.0;
  errEst = -12345.0;
  tmpvec = 0.0; corrLast = 0.0;

}
//======================================================================
//diagnostic routies and sanity checks
//======================================================================

template <class NL>
bool VBDF2Base<NL>::dataOk()
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
bool VBDF2Base<NL>::sizesOk()
{
  bool isOk = true;
  isOk = isOk && dataOk();
  if (isOk)
    {
      int neq = dae->getY0().size();
      isOk = isOk && ynp1Pred.size()  == neq;
      isOk = isOk && betaBDF.size()   == neq;
      isOk = isOk && yn.size()        == neq;
      isOk = isOk && ynPrime.size()   == neq;
      isOk = isOk && ynm1.size()      == neq;
      isOk = isOk && errEst.size()    == neq;
      isOk = isOk && tmpvec.size()    == neq;
      isOk = isOk && corrLast.size()  == neq;
    }

  return isOk;
}

template <class NL>
bool VBDF2Base<NL>::allDataOk()
{
  bool isOk = dataOk();
 
  isOk = isOk && ynp1;
  isOk = isOk && ynp1Prime;

  return isOk;
}

template <class NL>
bool VBDF2Base<NL>::ok()
{
  bool isOk = dataOk();
  
  isOk = isOk && sizesOk();

  return isOk;
}

//======================================================================
//main routines
//======================================================================
template <class NL>
bool VBDF2Base<NL>::calculateSolution(const real& tout,
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
      approxOrder = 1;
      predOrder   = 1;

      yn      = dae->getY0();
      tn      = dae->getT0();
      ynPrime = dae->getY0prime();

      calculateInitialDt(tn,tout);

#ifdef DEBUG_VBDF2
      std::cout<<"VBDF2 calcSoln("<<tout<<") dt0 = "<<deltaTnp1<<std::endl;
#endif
    }
  else
    {
      //otherwise assume that all the history is set and
      //that error estimate from previous solution is ok to
      //get step size
      //probably need to make this more robust
      bool lastStepOk = true;
      controller->setOrder(approxOrder); //\tilde{p} controller uses 
                                         //1/(\tilde{p}+1)
      deltaTnp1 = controller->estimateStepSize(locErr,deltaTn,lastStepOk);
#ifdef DEBUG_VBDF2
      std::cout<<"VBDF2 calcSoln("<<tout<<") tn = "<<tn<<" dt = "<<deltaTnp1
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
	  std::cerr<<"VBDF2 calcSoln ("<<tn<<" --> "<<tCurrent
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
bool VBDF2Base<NL>::step(const real& tOut,
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
  bool stepFailed      = false;

  int errorFailures    = 0;
  int nlSolverFailures = 0;
  int stepsFailed      = 0;

  //ynp1 and ynp1Prime hold solution at new time level.

  //borrow users storage
  ynp1      = &solutionAtTStep;
  ynp1Prime = &solutionPrime;

  //assumes that yn, tn, ynPrime, deltaTnp1 have all been set
  //by calling routine
  
  while (stepsFailed < maxFailures)
    {
      chooseOrder();
      //predict new solution and corresponding derivative
      stepFailed = calculatePredictor(ynp1Pred);
      *ynp1 = ynp1Pred;
      stepFailed = calculateBDFCoefs(); //sets alpha_k and \beta_k
      stepFailed = calculateYPrime(*ynp1,deltaTnp1,*ynp1Prime);
#ifdef DEBUG_VBDF2
      std::cout<<"VBDF2 step("<<tOut<<") tn = "<<tn<<" dt = "
	       <<deltaTnp1<<std::endl;
      std::cout<<"predicted value = "<<*ynp1<<std::endl;
#endif

      //pass new values on to dae
      dae->resetFunction();
      dae->yDaeDef     = *ynp1;
      dae->ypDaeDef    = *ynp1Prime;
      dae->alphaDaeDef = alphaBDF;
      dae->betaDaeDef  = betaBDF;
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
	  std::cerr<<"VBDF2  corrector failed "<<std::endl;
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

	  //make sure this got calculated correctly?
	  calculateYPrime(*ynp1,deltaTnp1,*ynp1Prime);

	  //determine if estimated truncation error is sufficiently low
	  calculateErrorEstimate(ynp1Pred,*ynp1);
	  bool errorOk = false;

	  controller->setOrder(approxOrder); //\tilde{p} controller uses 
                                             //1/(\tilde{p}+1)
	  real dtNew = controller->estimateStepSize(locErr,deltaTnp1,errorOk);

	  bool newDtFailed = checkStepSize(tn,tOut,dtNew);


	  if (newDtFailed)
	    return newDtFailed;

#ifdef DEBUG_VBDF2
	  std::cout<<"VBDF2 step("<<tOut<<") nlsolve to = "<<tn+deltaTnp1
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
		  std::cerr<<"VBDF2 NonIt step("<<tOut<<") nlsolve to = "
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
		  std::cerr<<"VBDF2 NonIt step("<<tOut<<") nlsolve to = "
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

	      ynm1      = yn;
	      ynPrime   = *ynp1Prime;

	      //
	      yn        = *ynp1;
	      
	      //record step
	      data->stepTaken(approxOrder,deltaTn,tn,locErr);
	      dae->stepTaken();

	      firstStep = false;
	      return false;
	    }//end successful step

	}//end nonlinear solver failure succeeded
    }//end while number of failures 
  //must have failed if got to this point
  data->stepTaken(approxOrder,deltaTnp1,tn);
  return true;

}
//======================================================================
//basic algorithm functions 
//======================================================================
template <class NL>
bool VBDF2Base<NL>::calculatePredictor(Vec& ypred)
{
  bool  predFailed = false;
  switch (predOrder)
    {
    case 1:
      ypred = yn;
      //y^{n+1,1} += \Delta t^{n+1} \dot{y}^{n}
      axpy(deltaTnp1,ynPrime,ypred);
      break;
    default:
      std::cerr<<"WARNING predOrder = "<<predOrder<<" not implemented "
	       <<" using 1st order approx "<<std::endl;
      ypred = yn;
      //y^{n+1,1} += \Delta t^{n+1} \dot{y}^{n}
      axpy(deltaTnp1,ynPrime,ypred);
      break;
    }//end switch on order

  return predFailed;
}

template <class NL>
bool 
VBDF2Base<NL>::calculateBDFCoefs()
{
  /***********************************************************************
     set \alpha_k and \vec \beta^{n}
     depending on order and step sizes
  **********************************************************************/
  real dtInv = 1.0/deltaTnp1;
  switch (approxOrder)
    {
    case 2:
      {
	real r = deltaTnp1/deltaTn;
	alphaBDF = (1.0+2.0*r)/(1.0+r)*dtInv;
	real b0  = -(1.0+r)*dtInv;
	real b1  = r*r/(1.0+r)*dtInv;
	betaBDF  = yn;
	scal(b0,betaBDF);
	axpy(b1,ynm1,betaBDF);
	break;
      }
    default: //first order
      {
	alphaBDF = dtInv;
	betaBDF  = yn;
	scal(-dtInv,betaBDF);
	break;
      }
    }//end switch
  return false;
}

template <class NL>
bool 
VBDF2Base<NL>::calculateYPrime(const Vec& y, 
					 const real& dt, 
					 Vec& yp)
{
  /***********************************************************************
     set \dot{\vec y}^{n+1} = \alpha_k\vec y^{n+1}_{(1)} + \vec \beta^{n}
 
     assuming alpha and \beta already determined based on order
  **********************************************************************/
  yp = betaBDF;
  axpy(alphaBDF,y,yp);

  return false;
}

template <class NL>
bool VBDF2Base<NL>::calculateErrorEstimate(const Vec& yp,
						     const Vec& y)
{
  /***************************************************

  calculate \vec{e}^{n+1} 
  depending on order of approximation

  Initially, use (\vec y^{n+1}-\vec y^{n+1,p})/2 for first order
  and 
    r/(1+r)(\vec y^{n+1}-(1+r)\vec y^{n}+r\vec y^{n-1}) for second order
  sets variables 
    errEst, locErr

  **************************************************/
  switch (approxOrder)
    {
    case 2:
      {
	real r = deltaTnp1/deltaTn;
	errEst = y;
	axpy(-(1.0+r),yn,errEst);
	axpy(r,ynm1,errEst);
	scal(r/(1.0+r),errEst);
	break;
      }
    default:
      {
	errEst = y;
	axpy(-1.0,yp,errEst);
	scal(0.5,errEst);
	break;
      }
    }//end switch

  assert(errNorm);
  tmpvec = y;
  for (int i=0; i < y.size(); i++)
    {
      if (std::fabs(yn[i]) > std::fabs(y[i]))
	tmpvec[i] = yn[i];
    }
  errNorm->setWeight(tmpvec);

  locErr = (*errNorm)(errEst);

  bool errEstFailed = locErr < 0.0;
  return errEstFailed;
}
template <class NL>
bool VBDF2Base<NL>::chooseOrder()
{
  if (firstStep)
    {
      approxOrder = 1;
      predOrder = 1;
    }
  else
    {
      approxOrder = 2;
      predOrder = 1; //for now keep linear predictor
    }
  return false;
}

template <class NL>
bool VBDF2Base<NL>::calculateInitialDt(const real& t0,
				       const real& tout)
{

   //follow approach from Shampine_94 pg 379
   //assumes ynPrime is set
   assert(errNorm);
   errNorm->setWeight(yn);

   real normY0Prime = (*errNorm)(ynPrime);
   real denom    = std::max(timeEps,normY0Prime);
   controller->setOrder(approxOrder); //\tilde{p} controller uses 
                                      //1/(\tilde{p}+1)
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
bool VBDF2Base<NL>::checkStepSize(const real& t0, 
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
