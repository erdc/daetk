#include "FLCBDF.h"

namespace Daetk 
{

using std::min;
using std::max;
using std::cerr;
using std::endl;
using std::cout;

  void FLCBDF::useFixedStep(const real step)
  {
    USE_FIXED_STEP=true;
    hFixed=step;
  }
  
  void FLCBDF::useFixedOrder(const int order)
  {
    USE_FIXED_ORDER=true;
    kFixed=order;
  }
  
  void FLCBDF::updateJacobians(){UPDATE_JACOBIANS=true;}
  void FLCBDF::useInterpolant(){STEP_EXACT=false;}
FLCBDF::FLCBDF():
  prepareJacobian(true),                  
  predictorAlreadyUpdated(false),                      
  firstCallToSolver(true),
  jacobianIsOld(true),
  correctorStepFailed(false),
  USE_FIXED_STEP(false),
  USE_FIXED_ORDER(false),
  UPDATE_JACOBIANS(false),
  STEP_EXACT(true),
  neq(0),
  k(1),
  kLast(0),
  constStepsTaken(0),
  stepsFailed(0),
  MAX_FAILURES(11),
  shortCutFactor(1),
  kFixed(1),
  h(0.0),
  hLast(-1.0),
  tn(0.0),
  tnPlusH(0.0),
  errorEstimate(1.0),
  norm_yCminusyP(1.0),
  convergenceFactor(1.0),
  hFixed(0.0),
  alphaS(0.0),
  alphaoNp1(0.0),
  alphaOld(0.0),
  alphaLast(0.0),
  startupType(0),
  stepIncreaseCeilingStartup(2.0),
  stepIncreaseFloorStartup(2.0),
  stepDecreaseCeilingStartup(0.9),
  stepDecreaseFloorStartup(0.5),
  forceRichardsonExtrapolationOrder(2),
  chooseStepSizeFlag(0),
  stepIncreaseCeiling(2.0),
  stepIncreaseFloor(2.0),
  stepDecreaseCeiling(0.9),
  stepDecreaseFloor(0.5),
  dae(0),
  data(0),
  linearSolver(0),
  nonlinearSolver(0),
  jacobian(0),
  weightedNorm(0),
  stepRunTime(0)
{
  //default constructor is called for all vectors
  Tracer tr("FLCBDF::FLCBDF()");
  int i;
  for (i=0;i<7;i++)
    {
      sigmaNp1[i] = 1.0;
      psiNp1[i]=1.0;
    }
  for (i=0;i<6;i++)
    {
      alphaNp1[i]=0.0;
      betaNp1[i]=1.0;
      gammaNp1[i]=0.0;
    }
  alphaNp1[0]=1.0;
}
  
FLCBDF::FLCBDF(DaeDefinition& daeIn,LinearSolver& linearSolverIn,
               NonlinearSolver& nonlinearSolverIn, 
               Jacobian& jacobianIn,VectorNorm& norm,
               DataCollector& dataIn): 
  prepareJacobian(true),                  
  predictorAlreadyUpdated(false),                      
  firstCallToSolver(true),
  jacobianIsOld(true),
  correctorStepFailed(false),
  USE_FIXED_STEP(false),
  USE_FIXED_ORDER(false),
  UPDATE_JACOBIANS(false),
  STEP_EXACT(true),
  neq(daeIn.getY0().dim()),
  k(1),
  kLast(0),
  constStepsTaken(0),
  stepsFailed(0),
  MAX_FAILURES(11),
  shortCutFactor(1),
  kFixed(1),
  h(0.0),
  hLast(-1.0),
  tn(0.0),
  tnPlusH(0.0),
  errorEstimate(1.0),
  norm_yCminusyP(1.0),
  convergenceFactor(1.0),
  hFixed(0.0),
  yn(0),
  ynprime(0),
  yCminusyP(daeIn.getY0().dim(),0.0),
  tempvec(daeIn.getY0().dim(),0.0),
  alphaS(0.0),
  alphaoNp1(0.0),
  alphaOld(0.0),
  alphaLast(0.0),  
  startupType(0),
  stepIncreaseCeilingStartup(2.0),
  stepIncreaseFloorStartup(2.0),
  stepDecreaseCeilingStartup(0.9),
  stepDecreaseFloorStartup(0.5),
  forceRichardsonExtrapolationOrder(2),
  chooseStepSizeFlag(0),
  stepIncreaseCeiling(2.0),
  stepIncreaseFloor(2.0),
  stepDecreaseCeiling(0.9),
  stepDecreaseFloor(0.5),
  norm_yByS(1.0),
  tmp(daeIn.getY0().dim(),0.0),
  yBigStep(daeIn.getY0().dim(),0.0),
  ySmallStep(daeIn.getY0().dim(),0.0),
  yExt(daeIn.getY0().dim(),0.0),
  ySminusyB(daeIn.getY0().dim(),0.0),
  dae(&daeIn),
  data(&dataIn),
  linearSolver(&linearSolverIn),
  nonlinearSolver(&nonlinearSolverIn),
  jacobian(&jacobianIn),
  weightedNorm(&norm),
  stepRunTime(0.0)
{
  Tracer tr("FLCBDF::FLCBDF(DaeDefinition& daeIn,NonlinearSolver& nonlinearSolverIn,BDFJacobian& jacobianIn,VectorNorm& norm,Data& dataIn)");

  //finish allocating vectors and initialize history
  int i;
  for (i=0;i<7;i++)
    {
//        divdifN[i].newsize(daeIn.getY0().dim());
//        divdifN[i] = 0.0;
//        divdifNm1[i].newsize(daeIn.getY0().dim());
//        divdifNm1[i] = 0.0;
      phiN[i].newsize(daeIn.getY0().dim());
      phiN[i] =  0.0;
      psiNp1[i]=1.0;
      sigmaNp1[i]=1.0;
    }
  
  for (i=0;i<6;i++)
    {
      alphaNp1[i]=0.0;
      betaNp1[i]=1.0;
      gammaNp1[i]=0.0;
    }
  alphaNp1[0]=1.0;
}

FLCBDF::~FLCBDF()
{
  Tracer tr("FLCBDF::~FLCBDF()");
}

bool FLCBDF::takeInitialSteps(const real& tout,real& tStep,Vec& 
                              solutionAtTStep, Vec& sp)
{
  //now this is a wrapper routine for calling different alternatives
  switch (startupType)
    {
    case 1:
      return takeInitialStepsWithModifiedCoefs(tout,tStep,solutionAtTStep,sp);
    case 2:
      return takeInitialStepsWithRichExtrap(tout,tStep,solutionAtTStep,sp);
    case 0:
      return takeInitialStepsOriginal(tout,tStep,solutionAtTStep,sp);
    default:
      return takeInitialStepsOriginal(tout,tStep,solutionAtTStep,sp);
    }
}
bool FLCBDF::takeInitialStepsOriginal(const real& tout,real& tStep,Vec& 
				      solutionAtTStep, Vec& sp)
{
  clock.reset();
  clock.start();

  yn = &solutionAtTStep;
  ynprime = &sp;

  tn=dae->getT0();

#ifndef USE_BLAS
  *yn=dae->getY0();
#else
  copy(dae->getY0(),*yn);
#endif

#ifndef USE_BLAS
  *ynprime=dae->getY0prime();
#else
  copy(dae->getY0prime(),*ynprime);
#endif

  weightedNorm->setWeight(*yn);
  //conservative estimate of first step
  if ((*weightedNorm)(*ynprime) > 0)
    h=min((1e-3)*fabs(tout-tn),0.5/(*weightedNorm)(*ynprime));
  else
    h=(1e-3)*fabs(tout-tn);
  h=copysign(h,tout-tn);
  alphaOld = 1e-20; //forces jacobian evaluation on first step
  if (USE_FIXED_STEP && hFixed < h)
    h=hFixed;

  /* Try to run the order up to 5 and double the step size at each
     step. Exit the function and start using normal error control if the
     error estimates look particularly bad.
     */

  int i;
  int maxOrd=5;
  if (USE_FIXED_ORDER)
    maxOrd=kFixed;
  for (i=1;i<maxOrd+1;i++)
    {
      errorFailures=0;
      while ( stepsFailed <= MAX_FAILURES )
        {
          if (i==1)
            {
              constStepsTaken = 0;
              /* Alter the history directly in order to take the first
                 step. Form a history that reflects a step of h*y0prime from a
                 previous value of yn - h*ynprime to get to the current value,
                 yn=y0; In this way the predictor corrector framework is used to take
                 a backward Euler step.
              */
  
              //update scalar coefficients with stepsize h for the "pseudo-step"

              updateCoefficients();

              hLast=0;
              kLast=1;
              k=1;
              

              updateCoefficients();
              //alter the history array directly on the first step.
#ifndef USE_BLAS
              phiN[0] = dae->getY0();
#else
              copy(dae->getY0(),phiN[0]);
#endif
#ifndef USE_BLAS
              phiN[1] = h*(dae->getY0prime());
#else
              copy(dae->getY0prime(),phiN[1]);
              scal(h,phiN[1]);
#endif
              phiN[2] = 0.0;
            }
          else
            updateCoefficients();
          predictor();
          //cout<<"in first call to corrector"<<endl<<flush;
          correctorStepFailed=corrector();
          if (correctorStepFailed) 
            {
              cerr<<"corrector Failed"<<endl;
              stepsFailed+=1;
              data->nonlinearSolverFailure();
              h*=.25;
              if (stepsFailed > 2)
                k=1;
              checkStepSize();
              resetCoefficients();
            }
          else if (errorForStepTooLarge())
            {
              stepsFailed+=1;
              data->errorFailure();
              chooseOrderForStep();
              if (errorFailures==1)
                {
                  h*=max( min(0.9,0.9*errorEstimate), 0.25);
                }
              else if (errorFailures==2)
                {
                  h*=.25;
                }
              else
                {
                  k=1;
                  h*=.25;
                }
              checkStepSize();
              resetCoefficients();
            }
          else 
            {
              tn += h;
              tStep  = tn;
              hLast = h;
              kLast = k;
              weightedNorm->setWeight(*yn);
	      //mwf from seg?
	      stepRunTime=clock.elapsed();
              data->stepTaken(k,h,tn,error,stepRunTime,0.0);
              dae->stepTaken();
              predictorAlreadyUpdated=false;
              chooseOrderForStep();
              if ( k < i )
                {
                  //exit initial phase succesfully but k < maxord
                  chooseStepSizeOriginal(tn,tout);
                  return 0;
                }
	      else if (tStep >= tout)//relies on direction--fix 
                return 0;
              else
                {
                  //take another intial phase step
                  stepsFailed = 0;
                  if (k!=5)
                    {
                      if (USE_FIXED_STEP)
                        chooseStepSizeOriginal(tn,tout);
                      else
                        h*=2.0;
                      checkStepSize();
                      if (USE_FIXED_ORDER)
                        {
                          if (k < kFixed )
                            k++;
                        }
                      else
                        k++;
                    }
                  break;
                }
            }
          if (stepsFailed >=  MAX_FAILURES)
            {
              data->stepTaken(k,h,tn,error);
              return 1;
            }
        }
    }
  /* Exited the loop successfully with k=5. Reset h to it's previous
     value and choose the stepsize according to the error estimates. */
  chooseStepSizeOriginal(tn,tout);
  return 0;
}
bool FLCBDF::takeInitialStepsWithModifiedCoefs(const real& tout,real& tStep,Vec& 
					       solutionAtTStep, Vec& sp)
{
  clock.reset();
  clock.start();

  yn = &solutionAtTStep;
  ynprime = &sp;

  tn=dae->getT0();

#ifndef USE_BLAS
  *yn=dae->getY0();
#else
  copy(dae->getY0(),*yn);
#endif

#ifndef USE_BLAS
  *ynprime=dae->getY0prime();
#else
  copy(dae->getY0prime(),*ynprime);
#endif

  weightedNorm->setWeight(*yn);
  //conservative estimate of first step
  if ((*weightedNorm)(*ynprime) > 0)
    h=min((1e-3)*fabs(tout-tn),0.5/(*weightedNorm)(*ynprime));
  else
    h=(1e-3)*fabs(tout-tn);
  h=copysign(h,tout-tn);
  alphaOld = 1e-20; //forces jacobian evaluation on first step
  if (USE_FIXED_STEP && hFixed < h)
    h=hFixed;

  /* Try to run the order up to 5 and double the step size at each
     step. Exit the function and start using normal error control if the
     error estimates look particularly bad.
     */

  int i;
  int maxOrd=5;
  if (USE_FIXED_ORDER)
    maxOrd=kFixed;
  for (i=1;i<maxOrd+1;i++)
    {
      errorFailures=0;
      while ( stepsFailed <= MAX_FAILURES )
        {
          if (i==1)
            {
              constStepsTaken = 0;
              /* Alter the history directly in order to take the first
                 step. Form a history that reflects a step of h*y0prime from a
                 previous value of yn - h*ynprime to get to the current value,
                 yn=y0; In this way the predictor corrector framework is used to take
                 a backward Euler step.
              */
  
              //update scalar coefficients with stepsize h for the "pseudo-step"

              updateCoefficients();

              hLast=0;
              kLast=1;
              k=1;
              

              updateCoefficients();
              //alter the history array directly on the first step.
#ifndef USE_BLAS
              phiN[0] = dae->getY0();
#else
              copy(dae->getY0(),phiN[0]);
#endif
#ifndef USE_BLAS
              phiN[1] = h*(dae->getY0prime());
#else
              copy(dae->getY0prime(),phiN[1]);
              scal(h,phiN[1]);
#endif
              phiN[2] = 0.0;
            }
          else
            updateCoefficients();
          predictor();
          //cout<<"in first call to corrector"<<endl<<flush;
          correctorStepFailed=corrector();
          if (correctorStepFailed) 
            {
              cerr<<"corrector Failed"<<endl;
              stepsFailed+=1;
              data->nonlinearSolverFailure();
              h*=.25;
              if (stepsFailed > 2)
                k=1;
              checkStepSize();
              resetCoefficients();
            }
          else if (errorForStepTooLarge())
            {
              stepsFailed+=1;
              data->errorFailure();
              chooseOrderForStep();
              if (errorFailures==1)
                {
                  h*=max( min(0.9,0.9*errorEstimate), 0.25);
                }
              else if (errorFailures==2)
                {
                  h*=.25;
                }
              else
                {
                  k=1;
                  h*=.25;
                }
              checkStepSize();
              resetCoefficients();
            }
          else 
            {
              tn += h;
              tStep  = tn;
              hLast = h;
              kLast = k;
              weightedNorm->setWeight(*yn);
	      //mwf from seg?
	      stepRunTime=clock.elapsed();
              data->stepTaken(k,h,tn,error,stepRunTime,0.0);
              //data->stepTaken(k,h,tn,error);
              dae->stepTaken();
              predictorAlreadyUpdated=false;
              chooseOrderForStep();
              if ( k < i )
                {
                  //exit initial phase succesfully but k < maxord
                  chooseStepSizeStartup(tn,tout);
                  return 0;
                }
	      else if (tStep >= tout)//relies on direction--fix 
                return 0;
              else
                {
                  //take another intial phase step
                  stepsFailed = 0;
                  if (k!=5)
                    {
		      chooseStepSizeStartup(tn,tout);
                      checkStepSize();
                      if (USE_FIXED_ORDER)
                        {
                          if (k < kFixed )
                            k++;
                        }
                      else
                        k++;
                    }
                  break;
                }
            }
          if (stepsFailed >=  MAX_FAILURES)
            {
              data->stepTaken(k,h,tn,error);
              return 1;
            }
        }
    }
  /* Exited the loop successfully with k=5. Reset h to it's previous
     value and choose the stepsize according to the error estimates. */
  chooseStepSizeStartup(tn,tout);
  return 0;
}

bool FLCBDF::takeInitialStepsWithRichExtrap(const real& tout,real& tStep,Vec& 
					    solutionAtTStep, Vec& sp)
{
  clock.reset();
  clock.start();

  yn = &solutionAtTStep;
  ynprime = &sp;

  tn=dae->getT0();

#ifndef USE_BLAS
  *yn=dae->getY0();
#else
  copy(dae->getY0(),*yn);
#endif

#ifndef USE_BLAS
  *ynprime=dae->getY0prime();
#else
  copy(dae->getY0prime(),*ynprime);
#endif

  weightedNorm->setWeight(*yn);
  //conservative estimate of first step
  if ((*weightedNorm)(*ynprime) > 0)
    h=min((1e-3)*fabs(tout-tn),0.5/(*weightedNorm)(*ynprime));
  else
    h=(1e-3)*fabs(tout-tn);
  h=copysign(h,tout-tn);
  alphaOld = 1e-20; //forces jacobian evaluation on first step
  if (USE_FIXED_STEP && hFixed < h)
    h=hFixed;

  /* Try to run the order up to 5 and double the step size at each
     step. Exit the function and start using normal error control if the
     error estimates look particularly bad.
     */

  int i;
  int maxOrd=2;//mwf check out this choice, whether or not used unless USE_FIXED_ORDER
  int stepsStabilized = 0;
  bool useExtrapolatedSolutionAsSolution =  false;
  int minStableRichExtrapSteps = 10;
  real richNLfailureFactor = 0.25;
  if (USE_FIXED_ORDER)
    maxOrd=kFixed;

  bool done = false;
  int iTimeStep = 0; //number of global time steps taken so far
  while (!done)
    {
      iTimeStep++;
      errorFailures = 0;
      stepsFailed   = 0; //this means corrector failures for now
      cout<<"RichExtrap starting step number "<<iTimeStep<<" errorFailures = "<<errorFailures<<endl;
      bool richExtrapStepDone = false; 
      
      while (stepsFailed <= MAX_FAILURES && !richExtrapStepDone) //failure if any one of the
					  //Richardson Extrapolation
					  //steps fails
	{
	  //loop through 1 whole step and 2 half steps
	  //   starting with the whole step when ii==1
	  //    continuing for two half steps ii==2, ii==3
	  for (int iiRichExtrapPhase = 1; iiRichExtrapPhase < 4; iiRichExtrapPhase++)
	    {
	      if (iiRichExtrapPhase == 2) h *= 0.5;
	      
	      cout<<"\t  starting phase "<<iiRichExtrapPhase<<" h= "<<h<<endl;
	      
	      if (iTimeStep == 1 && iiRichExtrapPhase != 3)
		{
		  constStepsTaken = 0;
		  /* Alter the history directly in order to take the first
		     step. Form a history that reflects a step of h*y0prime from a
		     previous value of yn - h*ynprime to get to the current value,
		     yn=y0; In this way the predictor corrector framework is used to take
		     a backward Euler step.
		  */
  
		  //update scalar coefficients with stepsize h for the "pseudo-step"
		  
		  updateCoefficients();
		  
		  hLast=0;
		  kLast=1;
		  k=1;
		  
		  updateCoefficients();
		  //alter the history array directly on the first step.
#ifndef USE_BLAS
		  phiN[0] = dae->getY0();
#else
		  copy(dae->getY0(),phiN[0]);
#endif
#ifndef USE_BLAS
		  phiN[1] = h*(dae->getY0prime());
#else
		  copy(dae->getY0prime(),phiN[1]);
		  scal(h,phiN[1]);
#endif
		  phiN[2] = 0.0;
		  
		}//starting new time step from scratch
	      else
		updateCoefficients();
	      
	      predictor();
	      correctorStepFailed=corrector();

	      if (correctorStepFailed)
		{
		  stepsFailed+=1;
		  data->nonlinearSolverFailure();
		  //if first step try to cut back time step and redo
		  if (iiRichExtrapPhase == 1)
		    {
		      h *= richNLfailureFactor;
		    }
		  else //rerun with smaller whole step if failure occurs on a half step
		    {  //this won't go undo first half step though
		      h *= 2.0;
		      h *= richNLfailureFactor;
		      
		    }
		  cerr<<"corrector Failed in Richardson Extrapolation phase h= "<<h<<endl;
		  checkStepSize();
		  resetCoefficients();
		  break;
		}//corrector failed
	      else
		{
		  if (iiRichExtrapPhase == 1)
		    {
		      copy(*yn,yBigStep);
		      resetCoefficients();
		    }
		  else if (iiRichExtrapPhase == 2)
		    {
		      //mwf put in step taken
		      tn += h;
		      tStep = tn;
		      hLast = h;
		      kLast = k;
		      weightedNorm->setWeight(*yn);
                      dae->stepTaken();
                      //seg put in CPU times
		      stepRunTime=clock.elapsed();
                      cout<<"time elapsed = "<<stepRunTime<<endl;
		      //mwf error and errorRichExt from last time step?
	              data->stepTaken(k,h,tn,error,stepRunTime,errorRichExt);
		      predictorAlreadyUpdated=false;
		    }
		  else
		    {
                      //seg put in CPU times
		      stepRunTime=clock.elapsed();
                      cout<<"time elapsed = "<<stepRunTime<<endl;
		      //mwf put in step taken
		      //weightedNorm->setWeight(*yn);
		      copy(*yn,ySmallStep);
		      int ldim=ySmallStep.ldim();
		      for (int jj=0;jj<ldim;jj++)
			{
			  ySminusyB[jj]= fabs(ySmallStep[jj]-yBigStep[jj]);
			}
		      norm_yByS=(*weightedNorm)(ySminusyB);
		      // extrapolate to more accurate solution for yn
		      if (useExtrapolatedSolutionAsSolution)
     			{
		      	  for (int jj=0;jj<ldim;jj++)
			     {
			       yExt[jj] = ySmallStep[jj] + 1.0/(pow(2.0,k)-1.0)*(ySminusyB[jj]);
			     }
			  //could make this different variable
			  norm_yByS *= 1.0/(pow(2.0,k)-1.0);
			  //if going to use extrapolated solution for evolution
			  copy(yExt,*yn);
		        }
		      richExtrapStepDone = true;
		    }//iiRichPhase 
		}//corrector did not fail
	    }//Richardson Extrapolation loop
	}//step failure while
      if (stepsFailed >=  MAX_FAILURES)
	{
	  //data->stepTaken(k,h,tn,error);
	  data->stepTaken(k,h,tn,error,stepRunTime,errorRichExt);
	  return 1;
	}
      //need to add check to see if Richardson step was ok error wise
      //put in steps for FLC eror calculation here to determine its error and errork
      bool flcbdfErrorTooLarge = errorForStepTooLarge();
      //mwf Feb 17 2009 I believe this is error for next time step
      real errork = sigmaNp1[k]*norm_yCminusyP;

      //computed a full time step with Richardson Extrap
      //h here is the small step 
      tn += h;
      tStep = tn;
      hLast = h;
      kLast = k;
      weightedNorm->setWeight(*yn);
 //      data->stepTaken(k,h,tn,error);
      dae->stepTaken();
      predictorAlreadyUpdated=false;
      chooseOrderForStep();
      cout<<"Richardson Extrap error norm = "<<norm_yByS<<", flcbdf errork= "<<errork<<" error= "<<error<<endl;  
      cout<<"RichExtrap step number "<<iTimeStep<<" took step to "<<tn<<" hLast= "<<hLast<<endl;
      if (tStep >= tout - 1.0e-12)//relies on direction--fix 
        {
	  //data->stepTaken(k,h,tn,error);
	  done = true;
	  data->stepTaken(k,h,tn,error,stepRunTime,errorRichExt);
	  //mwf debug
	  std::cout<<"Leaving RichExtrap Init tStep= "<<tStep<<" tout= "<<tout<<std::endl;
	  return 0;
        }
      else //continue on, use classical error control to pick next step size
	{
	  const real safety = 0.975;
	  const real richardsonExtrapCeiling = 5.0;
	  const real richardsonExtrapFloor   = 0.2;
	  errorRichExt = norm_yByS/(pow(2.0,k)-1.0);
	  real fac = pow(1.0/errorRichExt,1.0/(k+1.0));
	  fac = std::min(fac,richardsonExtrapCeiling);
	  fac = std::max(fac,richardsonExtrapFloor);
	  if (fac<2.0 && fac > 0.9)
	    stepsStabilized++;
	  data->stepTaken(k,h,tn,error,stepRunTime,errorRichExt);
	  //put in steps for FLC eror calculation here to determine its error and errork
	  //mwf don't need this again? Feb 17 2009
	  //bool flcbdfErrorTooLarge = errorForStepTooLarge();
	  //mwf or this Feb 17 2009
	  //real errork = sigmaNp1[k]*norm_yCminusyP;
	  // chooseOrderForStep();  seg--don't think it's needed because calculating errork above
          real errorRatio = errorRichExt/errork;
	  cout<<"Richardson Extrap error norm = "<<errorRichExt<<", flcbdf errork= "<<errork<<" error= "<<error<<"errorRatio = "<<errorRatio<<endl;  
	  cout<<" fac stepStabilized= "<<fac<<", "<<stepsStabilized<<endl;
	  //decide if time to quit
          if (stepsStabilized > minStableRichExtrapSteps && errorRatio > safety*2.0/k)
	    {
	      	  done = true;
	    }
	  else //not done, continue for 
	    {
	      h *= 2.0*fac; //recall h is small step size coming in
	      if (forceRichardsonExtrapolationOrder > 0)
		k = forceRichardsonExtrapolationOrder;

	      if (STEP_EXACT)
		{
		  //mwf put this in to see if I can get it to just go to tout                                           
		  if (tn+h > tout)
		    {
		      real htmp = fabs(tout-tn);
		      if (tn+htmp > tn)
			h = htmp;
		      
		    }
		}
 
	      cout<<"Continuing Richardson Extrap step = "<<iTimeStep<<" order= "<<k<<" new h = "<<h<<std::endl; 
	      checkStepSize();
	    }

	}
    }//startup phase end richards extrapolation phases
      
  //seg -- take out chooseStepSize() call, want to send RE step choice to integrator first.  
  // otherwise, this reduces stepsize based on flcbdf error before step() is taken outside of start-up routine.
  //chooseStepSize(tn,tout);  
  
  return 0;
}

bool FLCBDF::step(const real& tout,real& tStep,Vec& solutionAtTStep, Vec& stp)
{
  //clock.start();
  yn = &solutionAtTStep;
  ynprime = &stp;
  errorFailures=0;
  stepsFailed=0;
  while ( stepsFailed <= MAX_FAILURES )
    {
      updateCoefficients();
      predictor();
      correctorStepFailed=corrector();
      if (correctorStepFailed) 
        {
          cerr<<"corrector Failed"<<endl;
          stepsFailed+=1;
          data->nonlinearSolverFailure();
          h*=.25;
          if (stepsFailed>2)
            k=1;
          checkStepSize();
          resetCoefficients();
        }
      else if (errorForStepTooLarge())
        {
          errorFailures++;
          stepsFailed+=1;
          data->errorFailure();
          chooseOrderForStep();
          if (errorFailures==1)
            {
              h*=max( min(0.9,0.9*errorEstimate), 0.25);
            }
          else if (errorFailures==2)
            {
              h*=.25;
            }
          else
            {
              k=1;
              h*=.25;
            }
          checkStepSize();
          resetCoefficients();
        }
      //step was accepted; exit function here 
      else
        {
          tn += h;
          tStep  = tn;
          if (k>kLast)
            kRaisedOnLastStep=true;
          else
            kRaisedOnLastStep=false;
          hLast = h;
          kLast = k;
          //seg put in CPU times
	  stepRunTime=clock.elapsed();
          cout<<"time elapsed = "<<stepRunTime<<endl;
	  data->stepTaken(k,h,tn,error,stepRunTime,0.0);
          // data->stepTaken(k,h,tn,error);
          dae->stepTaken();
          weightedNorm->setWeight(*yn);
          predictorAlreadyUpdated=false;
          chooseOrderForStep();
          chooseStepSize(tn,tout);
          return 0;
        }
    }
  data->stepTaken(k,h,tn);
  return 1;//step failed
}


void FLCBDF::updateCoefficients()
{
  int i;
  real temp1(0),temp2(0);
  //The following coefficients are constants throughout the integration:
  //betaNp1[0]=1.0, sigmaNp1[0]=1.0, gammaNp1[0]=0.0, alphaNp1[0]=1.0;
 
  if (k==kLast && h==hLast)
      constStepsTaken+=1;
  else 
      constStepsTaken = 0;
  
  if (!predictorAlreadyUpdated)
    {
      if(kLast<5)
        {
#ifndef USE_BLAS
          phiN[kLast+1] = yCminusyP;
#else
	  copy(yCminusyP,phiN[kLast+1]);
#endif
        }
#ifndef USE_BLAS
      phiN[kLast] += yCminusyP;
#else
      axpy(1.0,yCminusyP,phiN[kLast]);
#endif
      for (i=kLast-1;i>=0;i--)
        {
#ifndef USE_BLAS
          phiN[i] += phiN[i+1];
#else
	  axpy(1.0,phiN[i+1],phiN[i]);
#endif
        }
    }
  
  if (predictorAlreadyUpdated)
    {
      int i;
      for (i=constStepsTakenLastUpdate+1;i<kLastUpdate+1;i++)  
	{
          temp1 = 1.0/betaNp1[i];
#ifndef USE_BLAS	  
	  phiN[i]*=temp1;
#else
	  scal(temp1,phiN[i]);
#endif
	}
      if (kLastUpdate >= constStepsTakenLastUpdate)
        for (i=1;i<kLastUpdate+1;i++)  
          psiNp1[i-1] = psiNp1[i] - hLastUpdate;
    }

  if (k >= constStepsTaken)
    {
      temp1 = h;
      for (i=1;i<k+1;i++) 
	{
	  temp2 = psiNp1[i-1];
	  psiNp1[i-1] = temp1;
	  betaNp1[i] = betaNp1[i-1]*psiNp1[i-1]/temp2;
	  temp1 = temp2+h;
	  alphaNp1[i] = h/temp1;
	  sigmaNp1[i] = real(i)*sigmaNp1[i-1]*alphaNp1[i];
	  gammaNp1[i] = gammaNp1[i-1]+alphaNp1[i-1]/h;
	}
      psiNp1[k] = temp1;
    }

  hLastUpdate = h;
  kLastUpdate = k;
  constStepsTakenLastUpdate = constStepsTaken;
  predictorAlreadyUpdated=false;
}


void FLCBDF::resetCoefficients()
{
  predictorAlreadyUpdated=true;
}  


void FLCBDF::predictor() 
{
#ifndef USE_BLAS
  *yn=phiN[0];
#else
  copy(phiN[0],*yn);
#endif
  *ynprime=0.0;
  //recall betaNp1[0] = 1.0;
  int i;
  for (i=constStepsTaken+1;i<k+1;i++)
    {
#ifndef USE_BLAS
      phiN[i]*=betaNp1[i];
#else
      scal(betaNp1[i],phiN[i]);
#endif
    }
  for (i=1;i<k+1;i++)  
    {
#ifndef USE_BLAS
      *ynprime+=gammaNp1[i]*phiN[i];
#else
      axpy(gammaNp1[i],phiN[i],*ynprime);
#endif
      
#ifndef USE_BLAS
      *yn+=phiN[i];
#else
      axpy(1.0,phiN[i],*yn);
#endif
    }
}


void FLCBDF::computeDeltaForJacobian()
{
  Vec::UnitStrideIterator deltai,yPprimei,yPi,weighti;
  dae->deltaVF = 0.0;
  deltai = dae->deltaVF.begin();
  yPprimei = ynprime->begin();
  yPi = yn->begin();
  weighti = weightedNorm->getWeightBegin();
  const real* weightEnd =  weightedNorm->getWeightEnd();
  //we only compute delta the length of the weights because
  //we're assuming that everything past the weights is excluded from the jac
  while (weighti < weightEnd)
    {
      *deltai=SQRT_MACHINE_EPSILON*
	max3(fabs(*yPi),fabs(h*(*yPprimei)),1.0/(*weighti));
      *deltai=-copysign(*deltai,h*(*yPprimei));
      //make sure roundoff doesn't change the y
      *deltai=*yPi - (*yPi - *deltai);
      //the minus here is for the correctArgument routine
      ++deltai;++yPprimei;++yPi;++weighti;
   }
}

bool FLCBDF::corrector()
{
  //Determine new leading coefficient

  alphaS = 0.0;
  alphaoNp1 = 0.0;
  for (int i=1;i<=k;i++)
    {
      alphaS-=1.0/(real(i));
      alphaoNp1-=alphaNp1[i-1];
    }

  alphaLast = dae->alphaDaeDef;
  dae->alphaDaeDef=(-alphaS/h);

  if (dae->alphaDaeDef != alphaLast)
    {
      nonlinearSolver->recomputeConvergenceRate();
    }
 
 //determine new nonlinear system of corrector equations

  dae->tDaeDef = tn+h;

#ifndef USE_BLAS
  dae->yDaeDef = *yn;
  dae->ypDaeDef = *ynprime;
  dae->betaDaeDef = *ynprime - dae->alphaDaeDef*(*yn);
#else
  copy(*yn,dae->yDaeDef);
  copy(*ynprime,dae->ypDaeDef);
  copy(*ynprime,dae->betaDaeDef);
  axpy(-dae->alphaDaeDef,*yn,dae->betaDaeDef);
#endif

  dae->resetFunction();

  bool evalError=false;
  //cout<<"in dae->value in corrector"<<endl;
#ifndef USE_BLAS
  tempvec=dae->value(evalError);
#else
  copy(dae->value(evalError),tempvec);
#endif

  if (evalError)
    {
      cerr<<"predictor has S or P out of range"<<endl;
      return evalError;
    }

  /* check to see if the jacobian should be re-evaluated for new equations
     if so reevalute them. */
  
  evalError = analyzeJacobian();
  if (evalError)
    {
      cerr<<"Failed to evaluate Jacobian"<<endl;
      return evalError;
    }
  bool didNotConverge=false;
  didNotConverge=nonlinearSolver->solve(yCminusyP,*dae);
  //cout<<"out of newton solve"<<endl<<flush;
  
  if ( didNotConverge && jacobianIsOld )
    {
#ifndef USE_BLAS
      dae->yDaeDef = *yn;
      dae->ypDaeDef= *ynprime;
      dae->betaDaeDef = *ynprime - dae->alphaDaeDef*(*yn);
#else
      copy(*yn,dae->yDaeDef);
      copy(*ynprime,dae->ypDaeDef);
      copy(*ynprime,dae->betaDaeDef);
      axpy(-dae->alphaDaeDef,*yn,dae->betaDaeDef);
#endif
      dae->initializeFunction(*yn,*ynprime,tempvec);
      evalError = evaluateJacobian();
      if (evalError)
        {
          cerr<<"Failed to evaluate the Jacobian"<<endl;
          return evalError;
        }
      didNotConverge=nonlinearSolver->solve(yCminusyP,*dae);
      //cout<<"out of newton solve"<<endl<<flush;
    }
  
  if (!didNotConverge)
    {
#ifndef USE_BLAS
      *yn = dae->yDaeDef;
      *ynprime = dae->ypDaeDef;
#else
      copy(dae->yDaeDef,*yn);
      copy(dae->ypDaeDef,*ynprime);
#endif
    }
  return didNotConverge;
}

bool FLCBDF::evaluateJacobian()
{
  data->jacobianEvaluation();
  computeDeltaForJacobian();
  bool evalError=false;
  evalError=jacobian->evaluate(*yn,tempvec);
  //cout<<"out of evaluate jacobian"<<endl<<flush;
  if (evalError)
    return evalError;
  convergenceFactor=1.0;
  nonlinearSolver->setConvergenceFactor(convergenceFactor);
  nonlinearSolver->recomputeConvergenceRate();
  alphaOld = dae->alphaDaeDef;// Jac = dF/dy + alpha dF/dy'
  //cout<<"ready for prepare"<<endl<<flush;
  if (linearSolver->prepare())
    cerr<<"error in linearSolver->prepare()"<<endl;
  jacobianIsOld = false;
  dae->initializeFunction(*yn,*ynprime,tempvec);
  return false;
}

bool FLCBDF::analyzeJacobian()
{
  //real lambda,temp1,temp2;
  const real temp1 = (1.0 - .25)/(1.0 + .25),
    temp2=1.0/temp1;

  convergenceFactor=2.0/(1.0 + dae->alphaDaeDef/alphaOld);
  nonlinearSolver->setConvergenceFactor(convergenceFactor);
  jacobianIsOld = true;
  if ( dae->alphaDaeDef/alphaOld < temp1 || dae->alphaDaeDef/alphaOld > temp2 || UPDATE_JACOBIANS)
    {
      //cout<<"evaluating jacobian"<<endl<<flush;
      return evaluateJacobian();
    }
  else
    {
      //cout<<"not evaluating jacobian"<<endl<<flush;
      return false;
    }
}

void FLCBDF::chooseOrderForStep()
{
  /*

    Change this so that you don't multiply and divide by k,k+1 etc 

    */
  //this statement avoids errors in log function

  if (norm_yCminusyP == 0.0 )
    {
      errorEstimate=2.0;//cause step size to be doubled
      return;
    }

#ifndef USE_BLAS
  tempvec=yCminusyP;
#else
  copy(yCminusyP,tempvec);
#endif
  //the following terms will be estimates of the remainder terms of 
  //Taylor series expansions of orders k-2,k-1,k, and k+1

  real termkm2=12345,termkm1=12345,termk=12345,termkp1=12345,errorkm2=12345,errorkm1=12345,errork=12345,errorkp1=12345;

  //compute error estimates
   
  errork = sigmaNp1[k]*norm_yCminusyP;
  termk=(k+1.0)*errork;
  if (USE_FIXED_ORDER && k==kFixed)
    {
       errorEstimate=pow((errork*2.0)+.0001,-1.0/(k+1.0));
       return;
    }

  if (k > 1)
    {
#ifndef USE_BLAS
      tempvec+=phiN[k];                        //tempvec=phiNp1~=tempvec+phiN
#else
      axpy(1.0,phiN[k],tempvec);
#endif
      errorkm1=sigmaNp1[k-1]*(*weightedNorm)(tempvec);
      termkm1=k*errorkm1;

      if ((k <= 2) && (termkm1 <= .5*termk))
        {
          k-=1;
          errorEstimate=pow((errorkm1*2.0)+.0001,-1.0/(k+1.0));
          return;
        }
      if (k > 2)
        {
#ifndef USE_BLAS
          tempvec+=phiN[k-1];
#else
          axpy(1.0,phiN[k-1],tempvec);
#endif
          errorkm2=sigmaNp1[k-2]*(*weightedNorm)(tempvec);
          termkm2=(k-1.)*errorkm2;
          if ((max(termkm1,termkm2)) <= termk)
	    {
	      k-=1;
              errorEstimate=pow((errorkm1*2.0)+.0001,-1.0/(k+1.0));
	      return;
	    }
	}
    }
  if ((k!=5) && (k < constStepsTaken) && (stepsFailed==0) 
      && !kRaisedOnLastStep)
    {
#ifndef USE_BLAS
      tempvec=yCminusyP-phiN[k+1];
#else
      copy(yCminusyP,tempvec);
      axpy(-1.0,phiN[k+1],tempvec);
#endif
      errorkp1=(1.0/(k+2.0))*(*weightedNorm)(tempvec);
      termkp1=errorkp1*(k+2.0);
      if (k>1)
        {   
          if(termkm1 <= min(termk,termkp1))
            {
              k-=1;
              errorEstimate=pow((errorkm1*2.0)+.0001,-1.0/(k+1.0));
              return;
            }
          else if(termkp1 >= termk)
            {
              errorEstimate=pow((errork*2.0)+.0001,-1.0/(k+1.0));
              return;
            }
          else
            {
              k+=1;
              errorEstimate=pow((errorkp1*2.0)+.0001,-1.0/(k+1.0));
	      return;
            }
	}
      //k==1
      else if (termkp1 >= .5*termk)
	{
	  errorEstimate=pow((errork*2.0)+.0001,-1.0/(k+1.0));
	  return;
	} 
      else                          
	{
	  k=2;
	  errorEstimate=pow((errorkp1*2.0)+.0001,-1.0/(k+1.0));
	  return;
	}
    }
  else
    {
      errorEstimate=pow((errork*2.0) + .0001,-1.0/(k+1.0));
    }
}

bool FLCBDF::errorForStepTooLarge()
{
  real M;
  M=max(alphaNp1[k],fabs(alphaNp1[k]+alphaS-alphaoNp1));
  norm_yCminusyP=(*weightedNorm)(yCminusyP);
  error = M*norm_yCminusyP;
  return (error > 1.0);
}

void FLCBDF::chooseStepSize(const real& tin, const real& tout)
{
  if (chooseStepSizeFlag != 1)
    {
      chooseStepSizeOriginal(tin,tout);
      return;
    }
  real r;
  //errorEstimate is really 1/errorEstimate
  r=errorEstimate;
  if (r >= stepIncreaseFloor) //worth increasing
    {
      if (r >= stepIncreaseCeiling)
	{
	  r = stepIncreaseCeiling;
	}
      h *= r;
      if (USE_FIXED_STEP && h > hFixed)
	h = hFixed;
    }
  else if (r < 1)
    {
      if (r >= stepDecreaseCeiling)
	{
	  r = stepDecreaseCeiling;
	}
      else if (r <= stepDecreaseFloor)
	{
	  r = stepDecreaseFloor;
	}
      h *= r;
    }
  else
    {
      r = 1.0;
    }
  if (STEP_EXACT)
    {
      //mwf put this in to see if I can get it to just go to tout
      if (tin+h > tout)
        {
          real htmp = fabs(tout-tin);
          if (tn+htmp > tn)
            h = htmp;
        }
    }
  checkStepSize();
}
void FLCBDF::chooseStepSizeOriginal(const real& tin, const real& tout)
{
  real r;
  //errorEstimate is really 1/errorEstimate
  r=errorEstimate;
  if (r>=2.0 ) 
    {
      r=2.0;
      h*=r;
      if (USE_FIXED_STEP && h > hFixed)
        h=hFixed;
    }
  else if (r>1.0 && r<2.0)
    {
      r=1.0;
      //don't change stepsize
    }
  else if (r<=1.0 && r>.5)
    {
      if(r>.9)
        r=.9;
      else
        r=r;
      h*=r;
    }
  else //r <= .5
    {
      r=.5;
      h*=r;
    }
  if (STEP_EXACT)
    {
      //mwf put this in to see if I can get it to just go to tout
      if (tin+h > tout)
        {
          real htmp = fabs(tout-tin);
          if (tn+htmp > tn)
            h = htmp;
        }
    }
  checkStepSize();
}
void FLCBDF::chooseStepSize()
{
  if (chooseStepSizeFlag != 1)
    {
      chooseStepSizeOriginal();
      return;
    }
  real r;
  //errorEstimate is really 1/errorEstimate
  r=errorEstimate;
  if (r >= stepIncreaseFloor) //worth increasing
    {
      if (r >= stepIncreaseCeiling)
	{
	  r = stepIncreaseCeiling;
	}
      h *= r;
      if (USE_FIXED_STEP && h > hFixed)
	h = hFixed;
    }
  else if (r < 1)
    {
      if (r >= stepDecreaseCeiling)
	{
	  r = stepDecreaseCeiling;
	}
      else if (r <= stepDecreaseFloor)
	{
	  r = stepDecreaseFloor;
	}
      h *= r;
    }
  else
    {
      r = 1.0;
    }
  checkStepSize();
}
void FLCBDF::chooseStepSizeOriginal()
{
  real r;
  //errorEstimate is really 1/errorEstimate
  r=errorEstimate;
  if (r>=2.0 ) 
    {
      r=2.0;
      h*=r;
      if (USE_FIXED_STEP && h > hFixed)
        h=hFixed;
    }
  else if (r>1.0 && r<2.0)
    {
      r=1.0;
      //don't change stepsize
    }
  else if (r<=1.0 && r>.5)
    {
      if(r>.9)
        r=.9;
      else
        r=r;
      h*=r;
    }
  else //r <= .5
    {
      r=.5;
      h*=r;
    }
  checkStepSize();
}
void FLCBDF::chooseStepSizeStartup(const real& tin, const real& tout)
{
  real r;
  //errorEstimate is really 1/errorEstimate
  r=errorEstimate;
  if (r >= stepIncreaseFloorStartup) //worth increasing
    {
      if (r >= stepIncreaseCeilingStartup)
	{
	  r = stepIncreaseCeilingStartup;
	}
      h *= r;
      if (USE_FIXED_STEP && h > hFixed)
	h = hFixed;
    }
  else if (r < 1)
    {
      if (r >= stepDecreaseCeilingStartup)
	{
	  r = stepDecreaseCeilingStartup;
	}
      else if (r <= stepDecreaseFloorStartup)
	{
	  r = stepDecreaseFloorStartup;
	}
      h *= r;
    }
  else
    {
      r = 1.0;
    }
  if (STEP_EXACT)
    {
      //mwf put this in to see if I can get it to just go to tout
      if (tin+h > tout)
        {
          real htmp = fabs(tout-tin);
          if (tn+htmp > tn)
            h = htmp;
        }
    }
  checkStepSize();
}

void FLCBDF::checkStepSize()
{
  if (tn+h == tn)
    {
      cerr<<"FLCBDF:step size has been reduced to machine precision--exiting"<<endl;
      exit(1);
    }
}

void FLCBDF::interpolant(const real& tout,Vec& yAtTout,Vec& yPAtTout)
{
  std::cerr<<"Interpolating *************************************"<<std::endl;
  int i;
  //this is C^0 interpolant see [1,p121]
  //the implementation is a bit tricky because it uses coefficient data
  //to compute t-t_i, i=n,...n-k+2. This eliminates the need to store t_i's
  //updateCoefficients();
  real toutminustn,D,C,gamma,psiN[7];
  for (i=0;i<7;i++)
    psiN[i] = psiNp1[i];
  updateCoefficients();
  toutminustn=tout-tn;
#ifndef USE_BLAS
  yAtTout = phiN[0];
#else
  copy(phiN[0],yAtTout);
#endif

  yPAtTout = 0.0;

  C=1.0;
  D=0.0;
  gamma = toutminustn/psiN[0];
  for (i=1;i<kLast+1;i++)
    {
      D=D*gamma+C/psiN[i-1];
      C=C*gamma;
      gamma=(toutminustn + psiN[i-1])/psiN[i];
#ifndef USE_BLAS
      yAtTout+=C*phiN[i];
#else
      axpy(C,phiN[i],yAtTout);
#endif

#ifndef USE_BLAS
      yPAtTout+=D*phiN[i];
#else
      axpy(D,phiN[i],yPAtTout);
#endif
    }
  resetCoefficients();
  //also need to fix a few things that reset coefficients doesn't know to do
  constStepsTaken-=1;// the interpolant doesn't count as a step so should not modify this
  

}  
bool FLCBDF::calculateSteadyStateSolution(Vec& solutionAtTout, Vec& sp, real sstol)
{
  data->startUserStep();

  bool solverFailed;
  real tCurrent,tout=1.0e6;//got to set tout to something
  real nrm2_sp0 = nrm2(sp);

  if (firstCallToSolver) 
    { 
      if (solverFailed=takeInitialSteps(tout,tCurrent,solutionAtTout,sp))
        {
          return solverFailed;
        }
      firstCallToSolver=false;
    }
  else
    tCurrent=tn;

  while ( nrm2(sp) > nrm2_sp0*sstol + sstol)
    {
      //make sure tout isn't too small...step size will at most be doubled
      if (tCurrent >= 0.5*tout)
        tout = 2.0*tout;
      solverFailed=step(tout,tCurrent,solutionAtTout,sp);
      if (solverFailed)
        {
          data->endUserStep();
          data->includeSolution(-1,solutionAtTout);
          return solverFailed;
        }
    }

  data->endUserStep();
  data->includeSolution(tout,solutionAtTout);
  return solverFailed=false;
}

bool FLCBDF::calculateSolution(const real& tout,Vec& solutionAtTout, Vec& sp)
{
  data->startUserStep();

  bool solverFailed;
  real tCurrent;
  if (firstCallToSolver) 
    { 
      if (solverFailed=takeInitialSteps(tout,tCurrent,solutionAtTout,sp))
        {
          return solverFailed;
        }
      firstCallToSolver=false;
    }
  else
    {
      tCurrent=tn;
      if (STEP_EXACT)
        {
          chooseStepSize(tn,tout);
        }
    }
//   std::cerr<<this<<" calculate solution "<<tCurrent<<'\t'<<tout<<std::endl;
  //mwf not clear that need this if takeInitialSteps goes all the way 
  //assert(tCurrent < tout);
  while ( tCurrent < tout)
    {
//       std::cerr<<this<<" stepping "<<tCurrent<<'\t'<<tout<<std::endl;
      solverFailed=step(tout,tCurrent,solutionAtTout,sp);
      if (solverFailed)
        {
//           std::cerr<<this<<" solver failed"<<std::endl;
          data->includeSolution(-1,solutionAtTout);
          return solverFailed;
        }
//       std::cerr<<this<<" after stepping "<<tCurrent<<'\t'<<tout<<'\t'<<tCurrent + h<<std::endl;
    }
  if (!STEP_EXACT)
    {
//       std::cerr<<this<<" Interpolating "<<std::endl;
      interpolant(tout,solutionAtTout,sp);
    }
  data->endUserStep();
//   std::cerr<<this<<" calling include solution"<<std::endl; 
  data->includeSolution(tout,solutionAtTout);
  return solverFailed=false;
}

bool FLCBDF::stepSpecial(const real& tout,real& tStep, Vec& solutionAtTout, 
                         Vec& sp)
{
  bool solverFailed;
  if (firstCallToSolver) 
    { 
      firstCallToSolver = false;
      yn = &solutionAtTout;
      ynprime = &sp;
      
      tn=dae->getT0();
      
#ifndef USE_BLAS
      *yn=dae->getY0();
#else
      copy(dae->getY0(),*yn);
#endif
      
#ifndef USE_BLAS
      *ynprime=dae->getY0prime();
#else
      copy(dae->getY0prime(),*ynprime);
#endif
      
      weightedNorm->setWeight(*yn);
      
      //conservative estimate of first step
      
      h=min((1e-3)*fabs(tout-tn),0.5/(*weightedNorm)(*ynprime));
      h=copysign(h,tout-tn);
      alphaOld = 1e-20; //forces jacobian evaluation on first step
      
      /* Alter the history directly in order to take the first
         step. Form a history that reflects a step of h*y0prime from a
         previous value of yn - h*ynprime to get to the current value,
         yn=y0; In this way the predictor corrector framework is used to take
         a backward Euler step.
         */
      
      //update scalar coefficients with stepsize h for the "pseudo-step"
      
      updateCoefficients();
      
      /* Now try to run the order up to 5 and double the step size at each
         step. Exit the function and start using normal error control if the
         error estimates look particularly bad.
         */
      
      hLast=0.0;
      kLast=1;
      k=1;
      
      //updateCoefficients for first step
      
      updateCoefficients();
      //alter the history array directly on the first step.
#ifndef USE_BLAS
      phiN[0] = *yn;
#else
      copy(*yn,phiN[0]);
#endif
#ifndef USE_BLAS
      phiN[1] = h*(*ynprime);
#else
      copy(*ynprime,phiN[1]);
      scal(h,phiN[1]);
#endif
      phiN[2] = 0.0;
      
      errorFailures=0;
      while ( stepsFailed <= MAX_FAILURES )
        {
          predictor();
          correctorStepFailed=corrector();
          if (correctorStepFailed) 
            {
              cerr<<"corrector Failed"<<endl;
              stepsFailed+=1;
              data->nonlinearSolverFailure();
              h*=.25;
              checkStepSize();
              resetCoefficients();
              updateCoefficients();
            }
          else if (errorForStepTooLarge())
            {
              stepsFailed+=1;
              data->errorFailure();
              chooseOrderForStep();
              if (errorFailures==1)
                {
                  h*=max( min(0.9,0.9*errorEstimate), 0.25);
                }
              else if (errorFailures==2)
                {
                  h*=.25;
                }
              else
                {
                  k=1;
                  h*=.25;
                }                  
              checkStepSize();
              resetCoefficients();
              updateCoefficients();
            }
          else 
            {
              tn += h;
              tStep  = tn;
              hLast = h;
              kLast = k;
              weightedNorm->setWeight(*yn);
              data->stepTaken(k,h,tn,error);
              dae->stepTaken();
              predictorAlreadyUpdated=false;
              chooseOrderForStep();
              chooseStepSize(tn,tout);
              updateCoefficients();
              break;
            }
        }
    }
  else
    {
      if (solverFailed=step(tout,tStep,solutionAtTout,sp))
        return solverFailed;
    }
  return solverFailed=false;
}


void FLCBDF::reset()
{
  firstCallToSolver=true;
  predictorAlreadyUpdated=false;
  jacobianIsOld=true;
  k=1;
  kLast = 0;
  constStepsTaken=0;
  stepsFailed=0;
  shortCutFactor=1;
  correctorStepFailed=false;
  h=0.0;
  hLast=-1.0;
  tn=0.0;
  tnPlusH=0.0;
  errorEstimate=1.0;
  norm_yCminusyP=1.0;
  convergenceFactor=1.0;  
  yn=0;
  ynprime=0;
  yCminusyP=0.0;
  alphaS=0.0;
  alphaoNp1=0.0;
  alphaOld=0.0;
  int i;

  for (i=0;i<7;i++)
    {
//        divdifN[i] = 0.0;
//        divdifNm1[i] = 0.0;
      phiN[i] =  0.0;
      psiNp1[i]=1.0;
      //psiN[i]=1.0;
      sigmaNp1[i]=1.0;
    }
  
  for (i=0;i<6;i++)
    {
      alphaNp1[i]=0.0;
      betaNp1[i]=1.0;
      gammaNp1[i]=0.0;
    }
  alphaNp1[0]=1.0;  
}
//--- begin utility routines for general stepping procedures ---
void FLCBDF::setStepSizeProcedure(const int& stepFlag)
{
  chooseStepSizeFlag = stepFlag;
}
void FLCBDF::setStepIncreaseFloor(const real& floorIn)
{
  stepIncreaseFloor = floorIn;
}
void FLCBDF::setStepIncreaseCeiling(const real& ceilIn)
{
  stepIncreaseCeiling = ceilIn;
}
void FLCBDF::setStepDecreaseFloor(const real& floorIn)
{
  stepDecreaseFloor = floorIn;
}
void FLCBDF::setStepDecreaseCeiling(const real& ceilIn)
{
  stepDecreaseCeiling = ceilIn;
}
//--- begin utility routines for startup procedures ---
void FLCBDF::setStartupProcedure(const int& startFlag)
{
  startupType = startFlag;
}
void FLCBDF::setStepIncreaseFloorForStartup(const real& floorIn)
{
  stepIncreaseFloorStartup = floorIn;
}
void FLCBDF::setStepIncreaseCeilingForStartup(const real& ceilIn)
{
  stepIncreaseCeilingStartup = ceilIn;
}
void FLCBDF::setStepDecreaseFloorForStartup(const real& floorIn)
{
  stepDecreaseFloorStartup = floorIn;
}
void FLCBDF::setStepDecreaseCeilingForStartup(const real& ceilIn)
{
  stepDecreaseCeilingStartup = ceilIn;
}
void FLCBDF::setRichardsonExtrapolationOrder(const int& richExtrapOrder)
{
  forceRichardsonExtrapolationOrder = richExtrapOrder;
}
}//Daetk
