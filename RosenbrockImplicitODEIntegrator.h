#ifndef ROSENBROCK_IMPLICIT_ODE_INTEGRATOR_H
#define ROSENBROCK_IMPLICIT_ODE_INTEGRATOR_H

#include <vector>
#include "Definitions.h"
#include "Vec.h"

#include "RosenbrockIntegratorBase.h"

#include "LinearSolver.h"
#include "Jacobian.h"

#include "RosenbrockDaeDefinition.h"
#include "BandColMat.h"
#include "PetscMat.h"
#include "PetscVecBlas.h"


namespace Daetk
{
//forward declarations

template <class MATREP>
class RosenbrockImplicitODEIntegrator : public RosenbrockIntegratorBase
{
  /***********************************************************************

    Implement basic Rosenbrock time integration scheme for implicit
    ODEs as formulated in \cite{Hairer_Wanner_96} and \cite{Lang_01},
    but try to solve in reduced form as suggested in \cite{Lang_01}

    \mat{M(t,\vec y)}\dot{\vec y} &=& \vec {F}(t,\vec y)

    that formally gets transformed into a system 
    \dot{\vec y} = \vec z
               0 = \mat{M}\vec{z} - \vec{F}(t,\vec y)



    The algorithm core is defined in RosenbrockIntegratorBase.

   IMPORTANT: The code assumes that the attached matrix is the same
   one assembled by the user definition class
   
  **********************************************************************/
public:

  RosenbrockImplicitODEIntegrator(RosenbrockDaeDefinition* daeIn, 
				  DataCollector* dataIn,
				  LinearSolver* linSolverIn,
				  Jacobian* jacIn,
				  MATREP* JyMinusMgIn,
				  VectorNorm* errNormIn,
				  VectorNorm* errNormAuxIn,
				  TimeIntegrationErrorController* 
				  errControllerIn,
				  const RosenbrockScheme& schemeIn, 
				  real absTolIn, 
				  real relTolIn,
				  real dt0MaxIn,
				  bool userDefinedMassMatrixIn = false,
				  bool userDefinedDFDtIn = false,
				  int auxErrFlagIn = 0);

  virtual ~RosenbrockImplicitODEIntegrator();

  virtual bool ok();

  virtual void reset();

  //---------- Integrator Interface ----------
  virtual bool calculateSolution(const real& tout,
				 Vec& solutionAtTout, 
				 Vec& solutionPrime);

  //take a single successful step to tOut
  //record t value reached in tStep
  virtual bool step(const real& tOut,real& tStep,Vec& solutionAtTStep, 
                    Vec& solutionPrime);

  //---------- specific to Rosenbrock Integrators ----------
  virtual bool usesAuxiliaryVariable() const
  { return true; }
  virtual bool getAuxiliaryVariable(Vec& aux); 

protected:

  //---------- Data Members ----------
  //linear solver
  LinearSolver* linSolver;
  //Jacobian \vec{F}_{y} for system
  Jacobian* jac;
  //matrix representation for \vec{F}_{y}-\frac{1}{\gamma\Detlta}
  //IMPORTANT: assumes this is the same matrix that is built by
  //Jacobian
  MATREP* JyMinusMg;
  //norm to use in error calculations for auxiliary variable
  VectorNorm* errNormAux;

  //how to incorporate error in auxiliary variable
  //0 --- use solution variable error only 
  //1 --- use auxiliary variable error only 
  //2 --- combine error as L2 or weighted L2 error
  //3 --- take max of error in soln variable and auxiliary variable
  int auxErrControl;

  //auxiliary variable and its error estimate
  Vec zn;      //z^n
  Vec znp1;    //z^{n+1}
  Vec zhatnp1; //\hat{z}^{n+1}
  //stage value for z
  Vec Zi;    //Z_i
  //err estimate for auxiliary variable
  Vec errEstZ;   //\vec e_z = z^{n+1}-\hat{z}^{n+1}
  real locErrCurZ;   //\|e\|
  //stage value for Y
  Vec Yi;    //Y_i


  bool userDefinedMassMatrix;
  bool userDefinedDFDt;

  //go from t^{0} to t^{0} + \Delta t starting with values y^0, yp^0
  virtual
  bool takeOneRosenbrockStep(const real& tIn,
			     const Vec& yIn,
			     const Vec& yInPrime,
			     const real& dt,
			     Vec& yOut,
			     Vec& yErrOut);


  //---------- Basic Functions For Algorithm ----------
  //calculate \vec{e_z}^{n+1} = 
  //  \vec{z}^{n+1}_{(1)}-\hat{\vec{z}}^{n+1}
  //and  \|\vec {e}^{n+1}\| 
  virtual bool calculateAuxiliaryErrorEstimate(const Vec& y1, const Vec& y2,
					       Vec& errZ, real& Dz); 

  //---------- utility functions ----------
  //dimension everything
  virtual bool setSizes();
  //check basic data dependencies
  virtual bool dataOk();
  //check that sizes match
  virtual bool sizesOk();
  //check everything including ynp1, ynp1Prime
  virtual bool allDataOk();


private:

  

};//end RosenbrockImplicitODEIntegrator


//======================================================================
//RosenbrockIntegratorBase implementation
//======================================================================
template <class MAT>
RosenbrockImplicitODEIntegrator<MAT>::
RosenbrockImplicitODEIntegrator(RosenbrockDaeDefinition* daeIn, 
				DataCollector* dataIn,
				LinearSolver* linSolverIn,
				Jacobian* jacIn,
				MAT* JyMinusMgIn,
				VectorNorm* errNormIn,
				VectorNorm* errNormAuxIn,
				TimeIntegrationErrorController* 
				errControllerIn,
				const RosenbrockScheme& schemeIn,
				real absTolIn, 
				real relTolIn,
				real dt0MaxIn,
				bool userDefinedMassMatrixIn,
				bool userDefinedDFDtIn,
				int auxErrFlagIn):
  RosenbrockIntegratorBase(daeIn,dataIn,errNormIn,errControllerIn,
			   schemeIn,absTolIn,relTolIn,dt0MaxIn),
  linSolver(linSolverIn),
  jac(jacIn), JyMinusMg(JyMinusMgIn), 
  errNormAux(errNormAuxIn),
  auxErrControl(auxErrFlagIn),
  zn(),znp1(),zhatnp1(),Zi(),errEstZ(),Yi(),
  userDefinedMassMatrix(userDefinedMassMatrixIn),
  userDefinedDFDt(userDefinedDFDtIn)
{
  assert(dataOk());

  bool resizeFailed = setSizes();
  assert(!resizeFailed);

}

template <class MAT>
RosenbrockImplicitODEIntegrator<MAT>::~RosenbrockImplicitODEIntegrator()
{

}


//======================================================================
//diagnostic routies and sanity checks
//======================================================================

template <class MAT>
bool RosenbrockImplicitODEIntegrator<MAT>::dataOk()
{
  bool isOk = RosenbrockIntegratorBase::dataOk();

  isOk = isOk && linSolver;
  isOk = isOk && jac;
  isOk = isOk && JyMinusMg;

  return isOk;
}
template <class MAT>
bool RosenbrockImplicitODEIntegrator<MAT>::setSizes()
{
  bool resizeFailed = false;
  if (!dataOk())
    {
      resizeFailed = true;
      return resizeFailed;
    }
  resizeFailed = RosenbrockIntegratorBase::setSizes();

  int neq = dae->getY0().size();

  zn.newsize(neq);
  znp1.newsize(neq);
  zhatnp1.newsize(neq);
  Zi.newsize(neq);
  errEstZ.newsize(neq);
  Yi.newsize(neq);

  return false;
}


template <class MAT>
bool RosenbrockImplicitODEIntegrator<MAT>::sizesOk()
{
  bool isOk = true;
  isOk = isOk && dataOk();
  isOk = RosenbrockIntegratorBase::sizesOk();
  if (isOk)
    {
      int neq = dae->getY0().size();
      isOk = isOk && JyMinusMg->dimRange()  == neq;
      isOk = isOk && JyMinusMg->dimDomain() == neq;
      isOk = isOk && zn.size()              == neq;
      isOk = isOk && znp1.size()            == neq;
      isOk = isOk && zhatnp1.size()         == neq;
      isOk = isOk && Zi.size()              == neq;
      isOk = isOk && errEstZ.size()         == neq;
      isOk = isOk && Yi.size()              == neq;
    }

  return isOk;
}

template <class MAT>
bool RosenbrockImplicitODEIntegrator<MAT>::allDataOk()
{
  bool isOk = dataOk();
 
  isOk = isOk && ynp1;
  isOk = isOk && ynp1Prime;

  return isOk;
}

template <class MAT>
bool RosenbrockImplicitODEIntegrator<MAT>::ok()
{
  bool isOk = dataOk();
  
  isOk = isOk && sizesOk();

  return isOk;
}
template <class MAT>
void RosenbrockImplicitODEIntegrator<MAT>::reset()
{

  RosenbrockIntegratorBase::reset();

  locErrCurZ = -12345.0;
  Yi      = -12345.0; 
  Zi      = -12345.0; 
  zn      = -12345.0; 
  znp1    = -12345.0;
  zhatnp1 = -12345.0;
  errEstZ = -12345.0;
}

template <class MAT>
bool RosenbrockImplicitODEIntegrator<MAT>::getAuxiliaryVariable(Vec& aux)
{
  if (aux.size() != znp1.size())
    aux.newsize(znp1.size());
  aux = znp1;
  return false;
}
//======================================================================
//main routines
//======================================================================
template <class MAT>
bool 
RosenbrockImplicitODEIntegrator<MAT>::calculateSolution(const real& tout,
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
      zn      = dae->getY0prime();//should this be Z or F?
      
      dae->rightHandSideValue(tn,yn,ynPrime);
      
      deltaTnp1 = calculateInitialDt(tn,tout);

#ifdef DEBUG_ROS
      std::cout<<"RoseImpOde calcSoln("<<tout<<") dt0 = "<<deltaTnp1<<std::endl;
#endif
    }
  else
    {
      //assumes yn and ynPrime have been set correctly
      //otherwise assume that all the history is set and
      //that error estimate from previous solution is ok to
      //get step size
      //probably need to make this more robust
      bool lastStepOk = true;
      deltaTnp1 = estimateStepSize(locErrCur,deltaTn,lastStepOk);
      assert(lastStepOk);

      bool newDtFailed = checkStepSize(tn,tout,deltaTnp1);
      if (newDtFailed)
	{
	  std::cerr<<"RoseCalcSoln step("<<tn<<") solve to = "
		   <<tout
		   <<" initial deltaTnp1= "<<deltaTnp1<<" failed "<<std::endl;
	  
	  stepFailed = newDtFailed;
	  return stepFailed;
	}

#ifdef DEBUG_ROS
      std::cout<<"In RoseImpOde calcSoln("<<tout<<") tn = "<<tn<<" dt = "
	       <<deltaTnp1<<std::endl;
#endif

    }
  tCurrent = tn;

  while (tCurrent < tout)
    {
      stepFailed = step(tout,tCurrent,solutionAtTout,solutionPrime);
      if (stepFailed)
	{
	  std::cerr<<"RosenbrockBase calcSoln ("<<tn<<" --> "<<tCurrent
		   <<") failed "<<std::endl;
          data->includeSolution(-1,solutionAtTout);
	  return stepFailed;
	}
    }

  data->endUserStep();
  data->includeSolution(tout,solutionAtTout);
  
  return false;
}
template <class MAT>
bool 
RosenbrockImplicitODEIntegrator<MAT>::step(const real& tOut,
					 real& tStep,Vec& solutionAtTStep, 
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

  //borrow users storage
  ynp1      = &solutionAtTStep;
  ynp1Prime = &solutionPrime;

  //assumes that yn, ynPrime, tn, deltaTnp1 have all been set
  //by calling routine
  
  while (stepsFailed < maxFailures)
    {

      //take a step from tn to deltaTnp1 using Rosenbrock algorithm
      //sets value of ynp1 and yerr
      bool roseStepFailed = takeOneRosenbrockStep(tn,yn,ynPrime,deltaTnp1,
						  *ynp1,yerr);

      if (roseStepFailed)
	{
	  std::cerr<<"RoseImpODE step ("<<tn<<" --> "<<tOut
		   <<" takOneRosenbrockStep failed "<<std::endl;
	  stepFailed = roseStepFailed;
	  return stepFailed;
	}

      //determine if estimated truncation error is sufficiently low
      bool estFailed = calculateErrorEstimate(*ynp1,yerr,errEst,locErrCur);
      if (estFailed)
	{
	  std::cerr<<"RoseImp step ("<<tn<<" --> "<<tOut
		   <<" ErrEstp failed "<<std::endl;
	  stepFailed = estFailed;
	  return stepFailed;
	}

      //what about auxiliary variable
      estFailed = calculateAuxiliaryErrorEstimate(znp1,zhatnp1,errEstZ,locErrCurZ);
      if (estFailed)
	{
	  std::cerr<<"RoseImp step ("<<tn<<" --> "<<tOut
		   <<" Aux ErrEstp failed "<<std::endl;
	  stepFailed = estFailed;
	  return stepFailed;
	}

      //how to incorporate error in auxiliary variable
      //0 --- use solution variable error only 
      //1 --- use auxiliary variable error only 
      //2 --- combine error as L2 or weighted L2 error
      //3 --- take max of error in soln variable and auxiliary variable

      switch (auxErrControl)
	{
	case 1:
	  {
	    locErrCur = locErrCurZ;
	    break;
	  }
	case 2:
	  {
	    real locErrComb = 
	      0.5*(locErrCur*locErrCur + locErrCurZ*locErrCurZ);
	    locErrCur = std::sqrt(locErrComb);
	    break;
	  }
	case 3:
	  {
	    real locErrComb = std::max(locErrCur,locErrCurZ);
	    locErrCur = locErrComb;
	    break;
	  }
	default:
	  {
	    //locErrCur already set
	    break;
	  }
	}//end switch
      
      //checks locErrCur
      bool errorOk = false; //mwf 5/12/06 now let error controller decide
      //sets value of locErrLastAcc, or locErrLastRej
      //and lastStepAccepted
      real dtNew       = estimateStepSize(locErrCur,deltaTnp1,errorOk);
#ifdef DEBUG_ROS_REJ
      //mwf test rejection phase
      if (tn > 0.1 && stepsFailed < 1)
	{
	  errorOk = false;
	  std::cerr<<"RoseBase step tn= "<<tn<<" deltaTnp1= "
		   <<deltaTnp1<<" setting errorOk= false= "
		   <<errorOk<<std::endl;
	}
#endif

#ifdef DEBUG_ROS
      std::cout<<"RoseImpODE step("<<tOut<<") solve to = "
	       <<tn+deltaTnp1<<std::endl;
      std::cout<<"locErrCur = "<<locErrCur<<" locErrCurZ= "
	       <<locErrCurZ<<" errorOk= "
	       <<errorOk<<" dtNew = "<<dtNew<<std::endl;
#endif
      
      if (!errorOk) //truncation error estimate too high
	{
	  //have to take another step
	  bool newDtFailed = checkStepSize(tn,tOut,dtNew);
	  errorFailures++;
	  stepsFailed++;
	  data->errorFailure();
	  if (newDtFailed)
	    {
	      std::cerr<<"RoseImpODE step("<<tOut<<") solve to = "
		       <<tn+deltaTnp1
		       <<" dtNew= "<<dtNew<<" failed "<<std::endl;
	      stepFailed = newDtFailed;
	      //record failures
	      data->stepTaken(approxOrder,dtNew,tn);
	      return stepFailed;
	    }
	  deltaTnp1 = dtNew;
	}
      else //successful step
	{
	  bool newDtFailed = checkStepSize(tn+deltaTnp1,tOut,dtNew);
	  if (newDtFailed)
	    {
	      std::cerr<<"RoseImpODE step("<<tOut<<") solve to = "
		       <<tn+deltaTnp1
		       <<" dtNew= "<<dtNew<<" failed "<<std::endl;

	      stepFailed = newDtFailed;
	      return stepFailed;
	    }
	  //cycle for next step
	  tn       += deltaTnp1;
	  //evaluate new yPrime
	  dae->rightHandSideValue(tn,*ynp1,*ynp1Prime);
	  //record jacobian evaluation
	  data->functionEvaluation();

	  tStep    = tn;
	  deltaTn   = deltaTnp1;
	  deltaTnp1 = dtNew;

	  //Y variables	  
	  ynPrime   = *ynp1Prime;
	  yn        = *ynp1;
	  //auxiliary variables
	  zn        = znp1;

	  //record step
	  data->stepTaken(approxOrder,deltaTn,tn,locErrCur);
	  dae->stepTaken();

	  firstStep = false;
	  stepFailed = false;
	  return stepFailed;
	}//end successful step
    }//end while number of failures 
  //must have failed if got to this point
  data->stepTaken(approxOrder,deltaTnp1,tn);
  return true;

}

//======================================================================
//basic algorithm functions 
//======================================================================

template <class MAT>
bool 
RosenbrockImplicitODEIntegrator<MAT>::takeOneRosenbrockStep(const real& tIn,
						     const Vec& yIn,
						     const Vec& yInPrime,
						     const real& dt,
						     Vec& yOut,
						     Vec& yErrOut)
{
  //assumes that tn, yn, ynPrime have been set correctly
  //sets value of ynp1 and yerr
  bool stepFailed = false;

  //reset dae and make sure its holding the values for yn and ynprime
  //that it integrator thinks it should hold
  dae->resetFunction();
  dae->yDaeDef = yIn;
  dae->tDaeDef = tIn;

  //evaluate F_y(t^n,y^n)
  //Important: assumes that JyMinusMg is attached to Jacobian so that 
  //           it gets set here
  //recall this is F_y for implicit ode not whole extended system
#ifdef ROSENBROCK_CALC_JAC_DELTA
  //RosenbrockDaeDefinition has to compute own Jacobian
  computeDeltaForJacobian();
#else
  dae->computeOwnJacobianDelta = true;
  dae->computeDeltaForJacobian();
#endif
  bool jacEvalFailed = jac->evaluate(yIn,yInPrime);
  //record jacobian evaluation
  data->jacobianEvaluation();
  if (jacEvalFailed)
    {
      std::cerr<<"RoseImplicitODE stepOne tn= "<<tn<<" dt= "<<dt
	       <<" JacEval failed "<<std::endl;
      stepFailed = jacEvalFailed;
      return stepFailed;
    }
  //now evaluate derivative of \mat{M(t,\vec y)}\vec z with respect to
  //\vec y and append it to jacobian \mat{F}_{\vec y}
  //\mat{A}^n = \mat{F}_{\vec y} - \pd{}{\vec y}[\mat{M(t,\vec y)}\vec z]
  
  bool dmassFailed = false;
  if (userDefinedMassMatrix)
    {
      dmassFailed = dae->appendMinusDMassDyDotZ(tIn,yIn,zn);
    }
  else
    {
      //assume Mass matrix is constant
    }
  if (dmassFailed)
    {
      std::cerr<<"RoseImplicitODE stepOne tn= "<<tn<<" dt= "<<dt
	       <<" DMassDy failed "<<std::endl;
      stepFailed = dmassFailed;
      return stepFailed;
    }
  
  //now include minus Mass matrix for actual system 
  // -\mat{J} = -\frac{1}{\Delta t \gamma}\mat{M} + \mat{A}^{n}  
  //mwf data collector needs to record this step somehow
  if (userDefinedMassMatrix)
    {
      dae->appendMinusMassMatrix(tIn,yIn,1.0/(gamma[0][0]*dt));
    }
  else
    {
      //assume that Mass Matrix is Identity
      int neq = JyMinusMg->dimRange();
      for (int k = 0; k < neq; k++)
	(*JyMinusMg)(k,k) += -1.0/(gamma[0][0]*dt);
    }

  bool dfdtFailed = false; 
  //evaluate F_t either numerically or analytically
  if (userDefinedDFDt)
    {
      //mwf data collector needs to record this step somehow
      dfdtFailed = dae->evaluateDFDt(tIn,yIn,Jt);
    }
  else
    {
      //numerically for now
      dfdtFailed = calculateDFDt(tIn,yIn,yInPrime,Jt);

    }

  //now have to get \mat{M_t}\dot \vec {z} and append to DFDt
  if (userDefinedMassMatrix)
    {
      dfdtFailed = dae->appendMinusDMassDtDotZ(tIn,yIn,zn,Jt);
    }
  else
    {
      //assume mass matrix is constant with respect to time
    }

  if (dfdtFailed)
    {
      std::cerr<<"RoseImplicitODE stepOne tn= "<<tn<<" dt= "<<dt
	       <<" DFDt failed "<<std::endl;
      stepFailed = dfdtFailed;
      return stepFailed;
    }

  //assumes that linear solver has matrix attached
  bool linearSolverFailed = linSolver->prepare();
      
  if (linearSolverFailed)
    {
      std::cerr<<"RoseImplicitODE oneStep LinearSolver prepare failed :"<<tn
	       <<" --> "<<tn+dt<<std::endl;
      stepFailed = true;
      data->linearSolverFailure();
      return stepFailed;
    }
 

  //loop through stages
  for (int i=0; i < numStages; i++)
    {

      //--- build right hand side of system to be solved
      //    make this negative of what I want see Lang V.5
      //get stage time level and stage Y_i and Z_i
      //t_i = t^n + \alpha_i\Delta t^{n+1} (see Lang pg 49 )
      real ti = tIn + alphaSum[i]*dt;

      //\vec Y_i = \vec y^n + \sum^{i-1}_{j=1}\vec Y_{nj} (see Lang pg 51)
      Yi = yIn;
      for (int j=0; j < i; j++)
	axpy(a[i][j],Yni[j],Yi);

      //\vec Z_i = (1-\sigma_i)\vec z^n + \sum^{i-1}_{j=1}
      //       (s_{ij}/\Delta t \vec Y_{nj})           (see Lang pg 51)
      Zi = zn;
      scal(1.0-sigmaSum[i],Zi);
      for (int j=0; j < i; j++)
	axpy(s[i][j]/dt,Yni[j],Zi);

      //use tmpvec to 
      //accumulate \sum_{j=1}^{i-1}\frac{c_{ij}}{\Delta t}\vec Y_{nj}
      tmpvec = 0.0;
      tmpprod= 0.0;
      for (int j=0; j < i; j++)
	{
	  axpy(c[i][j]/dt,Yni[j],tmpvec);
	}

      //mwf figure out if should count matrix vec multiplies too
      //now calculate 
      //\vec{M}\cdot\sum_{j=1}^{i-1}\frac{c_{ij}}{\Delta t}\vec Y_{nj}
      if (userDefinedMassMatrix)
	{
	  dae->applyMassMatrix(tn,yn,tmpvec,rhs);
	}
      else //assume mass matrix is identity
	rhs = tmpvec;



      //get user to evaluate rhs at stage values
      //F(t_i,\vec Y_i)
      tmpvec = 0.0;
      if (i == 0)
	{
	  tmpvec = yInPrime;
	}
      else
	{
	  dae->rightHandSideValue(ti,Yi,tmpvec);
	  //mwf count function eval
	  data->functionEvaluation();
	}
      //accumulate - F(t_i,Y_i) into rhs
      axpy(-1.0,tmpvec,rhs);
      
      //get negative derivative wrt to time
      //-\gamma_i\Delta t\vec{F}_{t}(t^n,\vec y^n)
      axpy(-dt*gammaSum[i],Jt,rhs);

      //last, have to apply difference in mass matrix at 
      //t^n and t^i. This is zero if first stage
      if (i > 0)
	{
	  tmpprod = 0.0;
	  dae->applyMassMatrix(tn,yn,Zi,tmpprod);
	  axpy(-1.0,tmpprod,rhs);
	  tmpprod = 0.0;
	  dae->applyMassMatrix(ti,Yi,Zi,tmpprod);
	  axpy(1.0,tmpprod,rhs);
	}
      
      //Finally ....
      //solve for stage value \vec Y_{ni} following Lang_01 V.11
      linearSolverFailed = linSolver->solve(rhs,Yni[i]);
      if (linearSolverFailed)
	{
	  std::cerr<<"RosImpODE LinearSolver solve failed stage "<< i<<" : "
		   <<tn<<" --> " <<tn+deltaTnp1<<std::endl;
	  data->linearSolverFailure();
	  return linearSolverFailed;
	}
      
    }//end stage loop

  //now get solution and error estimate at new time level
  yOut    = yIn;
  yErrOut = yIn;
  for (int i=0; i < numStages; i++)
    {
      axpy(m[i],Yni[i],yOut);
      axpy(mhat[i],Yni[i],yErrOut);
    }

  //also get the auxiliary solution at the  new level too  
  znp1   = zn;
  zhatnp1= zn;

  for (int i=0; i < numStages; i++)
    {
      //(sigma_i-1)z^n
      tmpvec = zn;
      scal(sigmaSum[i]-1.0,tmpvec);
      //\sum^{i}_{j=1}(c_{ij}-s_{ij})Y_{nj}
      tmpprod = 0.0;
      for (int j=0; j <= i; j++) //is this j<=i?
	axpy(c[i][j]-s[i][j],Yni[j],tmpprod);
      
      //accumulate sum/dt into tmpvec = (sigma-1)z^n
      axpy(1.0/dt,tmpprod,tmpvec);

      //now accumulate m_i times (sum/dt + (sigma_i z^n)) into 
      //z^{n+1}
      axpy(m[i],tmpvec,znp1);
      //\hat{z}^{n+1}
      axpy(mhat[i],tmpvec,zhatnp1);
    }

  return false;

}

//======================================================================
//basic algorithm functions 
//======================================================================
template <class MAT>
bool 
RosenbrockImplicitODEIntegrator<MAT>::
calculateAuxiliaryErrorEstimate(const Vec& y1,
				const Vec& y2,
				Vec& errZ,
				real& Dz)
{
  /***************************************************

  calculate \vec{e}^{n+1} = 
    \vec{y}^{n+1}_{(1)}-\hat{\vec{y}}^{n+1}
  and  \|\vec {e}^{n+1}\| using WRMS norm or other norm that's passed in

  Now follow RODAS lead and use max(|y_n|,|y^{n+1}|) in weighting
  **************************************************/
  assert(errNormAux);
  tmpvec = y1;
  for (int i=0; i < y1.size(); i++)
    {
      if (std::fabs(zn[i]) > std::fabs(y1[i]))
	tmpvec[i] = zn[i];
    }
  errNormAux->setWeight(tmpvec);


  errZ = y1;
  axpy(-1.0,y2,errZ);

  //follow Lubich and Roche 104 in Lang_01 suggestion
  scal(deltaTnp1,errZ);
  Dz = (*errNormAux)(errZ);

  bool errEstFailed = Dz < 0.0;
  return errEstFailed;
}


}//Daetk

#endif
