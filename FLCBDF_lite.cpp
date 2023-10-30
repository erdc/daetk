#include "FLCBDF_lite.h"
#include <iostream>
namespace Daetk 
{

using std::min;
using std::max;
using std::cerr;
using std::endl;
using std::cout;

void FLCBDF_lite::useFixedStep(const real step)
{
  USE_FIXED_STEP=true;
  hFixed=step;
}
  
void FLCBDF_lite::useFixedOrder(const int order)
{
  USE_FIXED_ORDER=true;
  kFixed=order;
}
  
void FLCBDF_lite::updateJacobians()
{
  UPDATE_JACOBIANS=true;
}

void FLCBDF_lite::useInterpolant()
{
  STEP_EXACT=false;
}

  //mwf add for building off diagonal jacobians etc
double FLCBDF_lite::getCurrentAlpha() const
{
  return alpha;
}
FLCBDF_lite::FLCBDF_lite():
  inInitializationPhase(true),
  iStep(1),
  prepareJacobian(true),                  
  predictorAlreadyUpdated(false),                      
  firstCallToSolver(true),
  jacobianIsOld(true),
  correctorStepFailed(false),
  USE_FIXED_STEP(false),
  USE_FIXED_ORDER(false),
  UPDATE_JACOBIANS(false),
  STEP_EXACT(false),
  neq(0),
  k(1),
  kLast(0),
  constStepsTaken(0),
  constStepsTakenSave(0),
  stepsFailed(0),
  MAX_FAILURES(11),
  shortCutFactor(1),
  kFixed(1),
  h(1.0),
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
  data(0),
  weightedNorm(0),
  useOldChooseStepSize(true),
  stepIncreaseCeiling(2.0),
  stepIncreaseFloor(2.0),
  stepDecreaseCeiling(0.9),
  stepDecreaseFloor(0.5) 
{
  Tracer tr("FLCBDF_lite::FLCBDF_lite()");
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
  
FLCBDF_lite::FLCBDF_lite(int dim,
                         VectorNorm& norm,
                         DataCollector& dataIn): 
  inInitializationPhase(true),
  iStep(1),
  firstStep(true),
  prepareJacobian(true),                  
  predictorAlreadyUpdated(false),                      
  firstCallToSolver(true),
  jacobianIsOld(true),
  correctorStepFailed(false),
  kRaisedOnLastStep(true),
  USE_FIXED_STEP(false),
  USE_FIXED_ORDER(false),
  UPDATE_JACOBIANS(false),
  STEP_EXACT(false),
  recomputeNonlinearConvergenceRate(true),
  stepFailed(false),
  neq(dim),
  k(1),
  kLast(0),
  kLastUpdate(0),
  constStepsTaken(0),
  constStepsTakenLastUpdateSave(0),
  stepsFailed(0),
  MAX_FAILURES(11),
  shortCutFactor(1),
  kFixed(1),
  h(1.0),
  hLast(-1.0),
  hLastUpdate(-1.0),
  tn(0.0),
  tnPlusH(0.0),
  errorEstimate(1.0),
  error(1.0),
  alpha(1.0),
  norm_yCminusyP(1.0),
  convergenceFactor(1.0),
  hFixed(0.0),
  yn(dim,0.0),
  ynprime(dim,0.0),
  yCminusyP(dim,0.0),
  tempvec(dim,0.0),
  alphaS(0.0),
  alphaoNp1(0.0),
  alphaOld(0.0),
  alphaLast(0.0),  
  data(&dataIn),
  weightedNorm(&norm),
  useOldChooseStepSize(true),
  stepIncreaseCeiling(10.0),
  stepIncreaseFloor(1.5),
  stepDecreaseCeiling(0.9),
  stepDecreaseFloor(0.5) 
{
  Tracer tr("FLCBDF_lite::FLCBDF_lite(DaeDefinition& daeIn,NonlinearSolver& nonlinearSolverIn,BDFJacobian& jacobianIn,VectorNorm& norm,Data& dataIn)");
  int i;
  for (i=0;i<7;i++)
    {
      phiN[i].newsize(dim);
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

FLCBDF_lite::~FLCBDF_lite()
{
  Tracer tr("FLCBDF_lite::~FLCBDF_lite()");
}

real FLCBDF_lite::chooseInitialStepSize(const real& t, 
                                        const real& tout,
                                        const Vec& y,
                                        const Vec& yPrime)
{
  firstStep=true;
  tn = t;
  weightedNorm->setWeight(y);
  if ((*weightedNorm)(yPrime) > 0)
    h=min((1e-3)*fabs(tout-tn),0.5/(*weightedNorm)(yPrime));
  else
    h=(1e-3)*fabs(tout-tn);
  h=copysign(h,tout-tn);
  alphaOld = 1e-20; //forces jacobian evaluation on first step
  if (USE_FIXED_STEP && hFixed < h)
    h=hFixed;
  return h;
}

real FLCBDF_lite::chooseDT(const real& t, 
                           const real& tout)
{
  if (inInitializationPhase == true)
    {
      chooseStepSize(tn,tout);
      initializationPhaseStep();
    }
  else
    {
      chooseStepSize(tn,tout);
      step();
    }
  return  h;
}

real FLCBDF_lite::setDT(const real& hin)
{
  if (firstStep==true)
    {
      h=hin;
    }
  else if (inInitializationPhase == true)
    {
      h = hin;
      resetCoefficients();
      initializationPhaseStep();
    }
  else
    {
      h = hin;
      resetCoefficients();
      step();
    }
  return  h;
}

real FLCBDF_lite::estimateError(const Vec& y)
{
  firstStep=false;
  if(inInitializationPhase)
    initializationPhaseCheckError(y);
  else
    checkError(y);
  ynprime = y;
  ynprime -= yn;
  ynprime *= alpha;
  ynprime += ynprime;
  yn = y;
}

/* Try to run the order up to 5 and double the step size at each
   step. Exit the function and start using normal error control if the
   error estimates look particularly bad.
*/
bool FLCBDF_lite::initializeTimeHistory(const Vec& yIn, const Vec& yPrimeIn)
{
  firstStep = true;
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
  phiN[0] = yIn;
#else
  copy(yIn,phiN[0]);
#endif
#ifndef USE_BLAS
  phiN[1] = h*(yPrimeIn);
#else
  copy(yPrimeIn,phiN[1]);
  scal(h,phiN[1]);
#endif
  phiN[2] = 0.0;

  predictor();

  //Determine new leading coefficient
  alphaS = 0.0;
  alphaoNp1 = 0.0;
  for (int i=1;i<=k;i++)
    {
      alphaS-=1.0/(real(i));
      alphaoNp1-=alphaNp1[i-1];
    }

  alphaLast = alpha;
  alpha=(-alphaS/h);

  if (alpha != alphaLast)
    recomputeNonlinearConvergenceRate=true;
  else
    recomputeNonlinearConvergenceRate=false;
}

/* Try to run the order up to 5 and double the step size at each
   step. Exit the function and start using normal error control if the
   error estimates look particularly bad.
*/
bool FLCBDF_lite::initializationPhaseStep()
{
  int maxOrd=5;
  if (USE_FIXED_ORDER)
    maxOrd=kFixed;
  updateCoefficients();
  predictor();

  //Determine new leading coefficient
  alphaS = 0.0;
  alphaoNp1 = 0.0;
  for (int i=1;i<=k;i++)
    {
      alphaS-=1.0/(real(i));
      alphaoNp1-=alphaNp1[i-1];
    }

  alphaLast = alpha;
  alpha=(-alphaS/h);

  if (alpha != alphaLast)
    recomputeNonlinearConvergenceRate=true;
  else
    recomputeNonlinearConvergenceRate=false;
}

bool FLCBDF_lite::initializationPhaseSolverFailure(int& stepsFailed)
{
  stepsFailed+=1;
  data->nonlinearSolverFailure();
  h*=.25;
  if (stepsFailed > 2)
    k=1;
  checkStepSize();
  resetCoefficients();
}

bool FLCBDF_lite::initializationPhaseCheckError(const Vec& y)
{
  //not really check error anymore
  int maxOrd=5;
  if (USE_FIXED_ORDER)
    maxOrd=kFixed;
  stepFailed=false;
  tn += h;
  data->stepTaken(k,h,tn);
  hLast = h;
  kLast = k;
  weightedNorm->setWeight(y);
  predictorAlreadyUpdated=false;
  chooseOrderForStep();
  if ( k < iStep )
    {
      inInitializationPhase=false;
      //exit initial phase succesfully but k < maxord
      return 0;
    }
  else
    {
      //take another intial phase step
      stepsFailed = 0;
      if (k!=5)
        {
          //               if (USE_FIXED_STEP)
          //                 chooseStepSize(tn,tout);
          //               else
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
    }
  iStep += 1;
  if (iStep == maxOrd+1)
    {
      //choose step size as usual
      inInitializationPhase=false;
    }
  errorFailures=0;
  stepsFailed=0;
  return 0;
}


bool FLCBDF_lite::step()
{
  updateCoefficients();
  predictor();

  //Determine new leading coefficient
  alphaS = 0.0;
  alphaoNp1 = 0.0;
  for (int i=1;i<=k;i++)
    {
      alphaS-=1.0/(real(i));
      alphaoNp1-=alphaNp1[i-1];
    }

  alphaLast = alpha;
  alpha=(-alphaS/h);

  if (alpha != alphaLast)
    recomputeNonlinearConvergenceRate=true;
  else
    recomputeNonlinearConvergenceRate=false;
}

double FLCBDF_lite::retryStep_solverFailure()
{
  stepsFailed+=1;
  data->nonlinearSolverFailure();
  h*=.25;
  if (stepsFailed>2)
    k=1;
  checkStepSize();
  resetCoefficients();
  if (inInitializationPhase == true)
    {
      initializationPhaseStep();
    }
  else
    {
      step();
    }
  return h;
}

double FLCBDF_lite::retryStep_errorFailure()
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
  if (inInitializationPhase == true)
    {
      initializationPhaseStep();
    }
  else
    {
      step();
    }
  return h;
}

bool FLCBDF_lite::checkError(const Vec& y)
{
  //not really checkError anymore
  tn += h;
  data->stepTaken(k,h,tn,error);
  if (k>kLast)
    kRaisedOnLastStep=true;
  else
    kRaisedOnLastStep=false;
  hLast = h;
  kLast = k;
  weightedNorm->setWeight(y);
  predictorAlreadyUpdated=false;
  chooseOrderForStep();
  errorFailures=0;
  stepsFailed=0;
  return 0;
}


void FLCBDF_lite::updateCoefficients()
{
  int i;
  real temp1(0),temp2(0);
  //The following coefficients are constants throughout the integration:
  //betaNp1[0]=1.0, sigmaNp1[0]=1.0, gammaNp1[0]=0.0, alphaNp1[0]=1.0;
 
  constStepsTakenSave = constStepsTaken;

  if (k==kLast && h==hLast)
      constStepsTaken+=1;
  else 
    {
      constStepsTaken = 0;
    }
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


void FLCBDF_lite::resetCoefficients()
{
  constStepsTaken = constStepsTakenSave;
  //constStepsTakenLastUpdate = constStepsTakenLastUpdateSave;
  predictorAlreadyUpdated=true;
}  


void FLCBDF_lite::predictor() 
{
#ifndef USE_BLAS
  yn=phiN[0];
#else
  copy(phiN[0],yn);
#endif
  ynprime=0.0;
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
      ynprime+=gammaNp1[i]*phiN[i];
#else
      axpy(gammaNp1[i],phiN[i],ynprime);
#endif
      
#ifndef USE_BLAS
      yn+=phiN[i];
#else
      axpy(1.0,phiN[i],yn);
#endif
    }
}


void FLCBDF_lite::computeDeltaForJacobian(Vec& deltaVF)
{
  Vec::UnitStrideIterator deltai,yPprimei,yPi,weighti;
  deltaVF = 0.0;
  deltai = deltaVF.begin();
  yPprimei = ynprime.begin();
  yPi = yn.begin();
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

bool FLCBDF_lite::calculate_yprime(const Vec& y,
                                   const Vec& Dy,
                                   Vec& yprime,
                                   Vec& Dyprime)
{
#ifndef USE_BLAS
  yprime = y;
  yprime -= yn;
  yprime *= alpha;
  yprime += ynprime;
  Dyprime = Dy;
  Dyprime *= alpha;
#else
  copy(y,yprime);
  axpy(-1.0,yn,yprime);
  scal(alpha,yprime);
  axpy(1.0,ynprime,yprime);
  copy(Dy,Dyprime);
  scal(alpha,Dyprime);
#endif
}

bool FLCBDF_lite::analyzeJacobian()
{
  bool updateJacobian=false;
  //real lambda,temp1,temp2;
  const real temp1 = (1.0 - .25)/(1.0 + .25),
    temp2=1.0/temp1;

  convergenceFactor=2.0/(1.0 + alpha/alphaOld);
//   nonlinearSolver->setConvergenceFactor(convergenceFactor);
  jacobianIsOld = true;
  if ( alpha/alphaOld < temp1 || alpha/alphaOld > temp2 || UPDATE_JACOBIANS)
    {
      updateJacobian=true;
    }
  else
    {
      updateJacobian=false;
    }
  return updateJacobian;
}

void FLCBDF_lite::chooseOrderForStep()
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

bool FLCBDF_lite::errorForStepTooLarge(const Vec& y)
{
  yCminusyP = y;
  yCminusyP -= yn;
  real M;
  M=max(alphaNp1[k],fabs(alphaNp1[k]+alphaS-alphaoNp1));
  norm_yCminusyP=(*weightedNorm)(yCminusyP);
  error = M*norm_yCminusyP;
  return (error > 1.0);
}

void FLCBDF_lite::chooseStepSize(const real& tin, const real& tout)
{
  real r;
  //errorEstimate is really 1/errorEstimate
  r=errorEstimate;
  cout<<"error estimate"<<r<<endl;
  if (useOldChooseStepSize)
    {
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
    }
  else
    {
      if (r >= stepIncreaseFloor)
        {
          if (r>=stepIncreaseCeiling ) 
            {
              r=stepIncreaseCeiling;
            }
          h*=r;
          if (USE_FIXED_STEP && h > hFixed)
            h=hFixed;
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
          h*=r;
        }
      else
        {
          r=1.0;
          //don't change stepsize
        }
    }
  if (STEP_EXACT)
    {
      if ( tin + h  > tout - SQRT_MACHINE_EPSILON*tout)
        {
          h = tout - tin;
          cerr<<"FLCBDF setting h to be exact"<<h<<endl;
        }
    }
  cout<<r<<'\t'<<h<<'\t'<<endl;
  checkStepSize();
}

void FLCBDF_lite::checkStepSize()
{
  if (tn+h == tn)
    {
      cerr<<"FLCBDF_lite:step size has been reduced to machine precision--exiting"<<endl;
      exit(1);
    }
}

void FLCBDF_lite::interpolant(const real& tout,Vec& yAtTout,Vec& yPAtTout)
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

void FLCBDF_lite::reset()
{
  inInitializationPhase = true;
  iStep = 1;
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

}//Daetk
