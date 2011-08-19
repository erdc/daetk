#include "RosenbrockIntegratorBase.h"

#include <cassert>
#include <iostream>

namespace Daetk
{

//======================================================================
//RosenbrockIntegratorBase implementation
//======================================================================

RosenbrockIntegratorBase::
RosenbrockIntegratorBase(RosenbrockDaeDefinition* daeIn, 
			 DataCollector* dataIn,
			 VectorNorm* errNormIn,
			 TimeIntegrationErrorController* errControllerIn,
			 const RosenbrockScheme& schemeIn,
			 real absTolIn, 
			 real relTolIn,
			 real dt0MaxIn):
  dae(daeIn), data(dataIn), 
  errNorm(errNormIn),
  scheme(schemeIn),
  numStages(1),
  approxOrder(1),
  absTol(absTolIn),relTol(relTolIn),dt0Max(dt0MaxIn),
  ynp1(0),ynp1Prime(0),
  Yni(),
  deltaTnp1(-12345.),
  Jt(),
  yn(), ynPrime(), yerr(),
  tn(-12345.), deltaTn(-12345.),
  locErrCur(-12345.0),
  errEst(),
  errController(errControllerIn),
  firstStep(true),
  timeTol(1.0),timeEps(1.0e-10),
  maxErrorFailures(12),
  maxFailures(12), 
  tmpvec(),tmpprod(),rhs()
{
  assert(dataOk());

  //make sure timeTol is consistent with error controller
  timeTol = errController->getTolerance();

  bool coefsFailed = setRosenbrockCoefficients(scheme);
  assert(!coefsFailed);
  bool resizeFailed = setSizes();
  assert(!resizeFailed);

}


RosenbrockIntegratorBase::~RosenbrockIntegratorBase()
{

}


bool RosenbrockIntegratorBase::setSizes()
{
  bool resizeFailed = false;
  if (!dataOk())
    {
      resizeFailed = true;
      return resizeFailed;
    }

  Yni.resize(numStages);
  int neq = dae->getY0().size();
  for (int i=0; i < numStages; i++)
    {
      Yni[i].newsize(neq);
    }

  Jt.newsize(neq);
  yn.newsize(neq);
  ynPrime.newsize(neq);
  yerr.newsize(neq);
  errEst.newsize(neq);
  tmpvec.newsize(neq);
  tmpprod.newsize(neq);
  rhs.newsize(neq);

  return false;
}


void RosenbrockIntegratorBase::reset()
{
  firstStep = true;
  ynp1 = 0;
  ynp1Prime = 0;

  deltaTnp1 = -12345.0;
  tn = -12345.0; deltaTn = -12345.0;
  locErrCur = -12345.0;
 
  Jt     = 0.0;
  yn     = -12345.0; 
  ynPrime= -12345.0;
  yerr   = -12345.0;
  errEst = -12345.0;
  tmpvec = 0.0; 
  tmpprod= 0.0; 
  rhs = 0.0;

  errController->reset();

}
//======================================================================
//diagnostic routies and sanity checks
//======================================================================


bool RosenbrockIntegratorBase::dataOk()
{
  bool isOk = true;

  isOk = isOk && dae;
  isOk = isOk && data;
  isOk = isOk && errNorm;
  isOk = isOk && errController;

  return isOk;
}


bool RosenbrockIntegratorBase::sizesOk()
{
  bool isOk = true;
  isOk = isOk && dataOk();
  if (isOk)
    {
      int neq = dae->getY0().size();
      for (int i=0; i < numStages; i++)
	{
	  isOk = isOk && Yni[i].size() == neq;
	}
      isOk = isOk && Jt.size()      == neq;
      isOk = isOk && yn.size()      == neq;
      isOk = isOk && ynPrime.size() == neq;
      isOk = isOk && yerr.size()    == neq;
      isOk = isOk && errEst.size()  == neq;
      isOk = isOk && tmpvec.size()  == neq;
      isOk = isOk && tmpprod.size() == neq;
      isOk = isOk && rhs.size()     == neq;

    }

  return isOk;
}


bool RosenbrockIntegratorBase::allDataOk()
{
  bool isOk = dataOk();
 
  isOk = isOk && ynp1;
  isOk = isOk && ynp1Prime;

  return isOk;
}


bool RosenbrockIntegratorBase::ok()
{
  bool isOk = dataOk();
  
  isOk = isOk && sizesOk();

  return isOk;
}

bool 
RosenbrockIntegratorBase::
setRosenbrockCoefficients(const RosenbrockScheme& schemeType)
{
  //set  numStages, approxOrder, gamma, c, gammaSum, m, mhat 
  //based on scheme trying
  //use tables from Lang_01

  bool setupOk = dataOk();
  if (!setupOk)
    return true;

  //zero everything by default
  for (int i=0; i < MAXSTAGES; i++)
    {
      m[i]        = 0.0;
      mhat[i]     = 0.0;
      gammaSum[i] = 0.0;
      alphaSum[i] = 0.0;
      sigmaSum[i] = 0.0;
      for (int j = 0; j < MAXSTAGES; j++)
	{
	  gamma[i][j] = 0.0;
	  a[i][j]     = 0.0;
	  c[i][j]     = 0.0;
	  s[i][j]     = 0.0;
	}
    }
  if (schemeType == ROS2) 
    {
      //second order 2 stage algorithm, L-stable see also
      //Hundsdorfer_Verwer_03, Formulas in Lang_01 appendix look
      //a little off
      numStages   = 2;
      approxOrder = 2;

      //gamma = 1 + sqrt(2)/2
      //b_1 = 1 - b_2
      //alpha_21 = 1/(2*b_1)
      //gamma_21 = -gamma/b_2

      //here b_i = 1/2
      //\hat{b}_1 = 1, \hat{b}_2 = 0
      //\Gamma^{-1} = 0.58578643762690                  0
      //              1.17157287525381   0.58578643762690
      /*****************************************
        below, I've used the values from Lang_01
      ******************************************/

      //\gamma
      real gconst = 1.707106781186547;
      gamma[0][0] = gconst;
      gamma[1][0] = -gconst/0.5;

      alphaSum[0] = 0.0;
      alphaSum[1] = 1.0;

      //(a_{ij})^{s}_{i,j=1} = (\alpha_{ij})^{s}_{i,j=1}\gvec{\Gamma}^{-1}
      a[0][0]     = 0.0;
      a[1][0]     = 0.5857864376269050;
      a[1][1]     = 0.0;

      //\Gamma^{-1}
      c[0][0]     = 0.5857864376269050;
      c[1][0]     = 1.171572875253810;
      c[1][1]     = 0.5857864376269050;

      //sum^{i}_{j=1} \Gamma_{ij}
      gammaSum[0] = 1.707106781186547;
      gammaSum[1] =-1.707106781186547;

      //\vec b^{T}\Gamma^{-1}
      m[0] = 0.8786796564403575;
      m[1] = 0.2928932188134525;

      //\hat{\vec b}^{T}\Gamma^{-1} 
      mhat[0] = 0.5857864376269050;
      mhat[1] = 0.0;

      //\sigma\Gamma^{-1}
      s[0][0] = 0.000000000000000;
      s[1][0] = 3.431457505076198e-1;
      s[1][1] = 0.000000000000000;

      //\sigma_i
      sigmaSum[0] = 0.000000000000000;
      sigmaSum[1] = 0.5857864376269050;

    }
  else if (schemeType == ROS3P) 
    {
      //third order, three stage algorithm, A-stable R(\infty) \approx 0.73
      //see Lang_Verwer_01
      //
      numStages   = 3;
      approxOrder = 3;

      //gamma = 1/2 + sqrt(3)/6
      //b_1 = 2/3, b_2 = 0, b_3 = 1/3
      //gamma_{21} = -1
      //gamma_{32} = 1/2 - 2gamma, gamma_{31} = -gamma

      //
      //\hat{b}_1 = 1/3, \hat{b}_2 = 1/3, \hat{b}_3 = 1/3


      //\gamma
      real gconst = 7.886751345948129e-1;
      gamma[0][0] = gconst;
      gamma[1][0] = -1.00000000000000e0;
      gamma[2][1] = 0.5 - 2.0*gconst;
      gamma[2][0] = -gconst;

      alphaSum[0] = 0.000000000000000e0;
      alphaSum[1] = 1.000000000000000e0;
      alphaSum[2] = 1.000000000000000e0;

      //(a_{ij})^{s}_{i,j=1} = (\alpha_{ij})^{s}_{i,j=1}\gvec{\Gamma}^{-1}
      a[1][0]     = 1.267949192431123e0;
      a[2][0]     = 1.267949192431123e0;
      a[2][1]     = 0.000000000000000e0;

      //\Gamma^{-1}
      c[0][0]     = 1.267949192431123;
      c[1][0]     = 1.607695154586736e0;
      c[1][1]     = 1.267949192431123;
      c[2][0]     = 3.464101615137755e0;
      c[2][1]     = 1.732050807568877e0; 
      c[2][2]     = 1.267949192431123;

      //sum^{i}_{j=1} \Gamma_{ij}
      gammaSum[0] =  7.886751345948129e-1; 
      gammaSum[1] = -2.113248654051871e-1;
      gammaSum[2] = -1.077350269189626e0;

      //\vec b^{T}\Gamma^{-1}
      m[0] = 2.000000000000000e0;
      m[1] = 5.773502691896258e-1;
      m[2] = 4.226497308103742e-1;

      //\hat{\vec b}^{T}\Gamma^{-1} 
      mhat[0] = 2.113248654051871e0;
      mhat[1] = 1.000000000000000e0;
      mhat[2] = 4.226497308103742e-1;

      //\sigma\Gamma^{-1}
      s[0][0] = 0.000000000000000;
      s[1][0] = 1.607695154586736;
      s[1][1] = 0.000000000000000;
      s[2][0] = 1.607695154586736;
      s[2][1] = 0.000000000000000;
      s[2][2] = 0.000000000000000;

      //\sigma_i
      sigmaSum[0] = 0.000000000000000;
      sigmaSum[1] = 1.267949192431123;
      sigmaSum[2] = 1.267949192431123;

    }
  else if (schemeType == RODAS4) 
    {
      //Rosenbrock solver from Hairer and Wanner. Method is order 4(3)
      //
      numStages   = 6;
      approxOrder = 4;


      //\gamma
      real gconst = 0.2500000000000000;
      gamma[0][0] = gconst;
      //I calculated gamma by backing out from c
      gamma[0][0] =   gconst;
      gamma[1][0] =  -3.543000000000029e-1;
      gamma[1][1] =   gconst;
      gamma[2][0] =  -1.336025052700005e-1;
      gamma[2][1] =  -1.289749473199982e-2;
      gamma[2][2] =   gconst;
      gamma[3][0] =   1.526849173000010e0;
      gamma[3][1] =  -5.336562887500030e-1; 
      gamma[3][2] =  -1.279392884300008e0;
      gamma[3][3] =   gconst;
      gamma[4][0] =   6.981190951800040e0;
      gamma[4][1] =  -2.092930097000013e0;
      gamma[4][2] =  -5.870067663000035e0;
      gamma[4][3] =   7.318068082500054e-1;
      gamma[4][4] =   gconst;
      gamma[5][0] =  -2.080189494200011e0;
      gamma[5][1] =   5.957623556800035e-1; 
      gamma[5][2] =   1.701617798300010e0;
      gamma[5][3] =  -8.851451983600143e-2;
      gamma[5][4] =  -3.786761399300002e-1; 
      gamma[5][5] =  gconst;

      //go ahead and read off \Gamma^{-1} = c directly

      alphaSum[0] = 0.0;
      alphaSum[1] = 0.3860000000000000;
      alphaSum[2] = 0.2100000000000000;
      alphaSum[3] = 0.6300000000000000;
      alphaSum[4] = 1.0000000000000000;
      alphaSum[5] = 1.0000000000000000;

      //(a_{ij})^{s}_{i,j=1} = (\alpha_{ij})^{s}_{i,j=1}\gvec{\Gamma}^{-1}
      a[0][0] =   0.0;
      a[1][0] =   1.5440000000000000;
      a[2][0] =   0.9466785281022800;
      a[2][1] =   0.2557011699000000;
      a[3][0] =   3.3148251870787937;
      a[3][1] =   2.8961240159798773;
      a[3][2] =   0.9986419140000000;
      a[4][0] =   1.2212245090707980;
      a[4][1] =   6.0191344810926299;
      a[4][2] =  12.5370833291149566;
      a[4][3] =  -0.6878860361200000;
      a[5][0] =   1.2212245092986209;
      a[5][1] =   6.0191344813485754;
      a[5][2] =  12.5370833293196604;
      a[5][3] =  -0.6878860360800001;
      a[5][4] =   1.0000000000000000;

      //\Gamma^{-1}
      c[0][0] =   1.0/gconst;
      c[1][0] =   5.6688000000000000;
      c[1][1] =   1.0/gconst;
      c[2][0] =   2.4300933568670464;
      c[2][1] =   0.2063599157120000;
      c[2][2] =   1.0/gconst;
      c[3][0] =   0.1073529065055983;
      c[3][1] =   9.5945622510667228;
      c[3][2] =  20.4702861487999996;
      c[3][3] =   1.0/gconst;
      c[4][0] =  -7.4964433159050206;
      c[4][1] =  10.2468043146053738;
      c[4][2] =  33.9999035259299589;
      c[4][3] = -11.7089089319999999;
      c[4][4] =   1.0/gconst;
      c[5][0] =  -8.0832467990118602;
      c[5][1] =   7.9811329880455499;
      c[5][2] =  31.5215943254324245;
      c[5][3] =  -16.3193054312706352;
      c[5][4] =   6.0588182388799998;
      c[5][5] =   1.0/gconst;

      //sum^{i}_{j=1} \Gamma_{ij}
      gammaSum[0] = 0.2500000000000000;
      gammaSum[1] =-0.1043000000000000;
      gammaSum[2] = 0.1034999999980000;
      gammaSum[3] =-0.0362000000000000;
      gammaSum[4] = 0.0000000000000000;
      gammaSum[5] = 0.0000000000000000;
 
      //\vec b^{T}\Gamma^{-1}
      m[0] =   1.2212245092981915;
      m[1] =   6.0191344813101981;
      m[2] =  12.5370833292377792;
      m[3] =  -0.6878860360960002;
      m[4] =   1.0000000000000000;
      m[5] =   1.0000000000000000;

      //\hat{\vec b}^{T}\Gamma^{-1} 
      mhat[0] =   1.2212245090707956;
      mhat[1] =   6.0191344810926299;
      mhat[2] =  12.5370833291149584;
      mhat[3] =  -0.6878860361200001;
      mhat[4] =   1.0000000000000000;
      mhat[5] =   0.0000000000000000;

      s[0][0] =   0.0;
      s[1][0] =   6.175999999999998e0;
      s[1][1] =   0.0;
      s[2][0] =   3.657022479035840e0;
      s[2][1] =   1.022804679600000e0;
      s[2][2] =   0.0;
      s[3][0] =   1.056512380050520e1;
      s[3][1] =   1.076916010223512e1;
      s[3][2] =   3.994567656000000e0;
      s[3][3] =   0.0;
      s[4][0] =  -6.356365049150499e0;
      s[4][1] =   1.464869134576237e1;
      s[4][2] =   3.881491663021866e1;
      s[4][3] =  -2.751544144480000e0;
      s[4][4] =   0.0;
      s[5][0] =  -7.496443315793021e0;
      s[5][1] =   1.024680431541504e1;
      s[5][2] =   3.399990352740777e1;
      s[5][3] =  -1.170890893184000e1;
      s[5][4] =   4.000000000000000e0;
      s[5][5] =   0.0;

      //\sigma_i
      sigmaSum[0] = 0.0;
      sigmaSum[1] = 1.544000000000000e0;
      sigmaSum[2] = 8.075770916766799e-1;
      sigmaSum[3] = 1.931495303851188e0;
      sigmaSum[4] = 1.000000000000000e0;
      sigmaSum[5] = 1.000000000000000e0;

    }
  else if (schemeType == ROWDAIND2) 
    {
      //Rosenbrock solver from Roche and Lubich and Roche
      //Formulas taken from Lang_01 appendix
      //
      numStages   = 4;
      approxOrder = 3;


      //\gamma
      real gconst = 0.300000000000000;
      gamma[0][0] = gconst;
      //I could calculate gamma by backing out from c
      //go ahead and read off \Gamma^{-1} = c directly

      alphaSum[0] = 0.0;
      alphaSum[1] = 0.5000000000000000;
      alphaSum[2] = 1.0000000000000000;
      alphaSum[3] = 1.0000000000000000;


      //(a_{ij})^{s}_{i,j=1} = (\alpha_{ij})^{s}_{i,j=1}\gvec{\Gamma}^{-1}
      a[0][0] =   0.0;
      a[1][0] =   1.6666666666666667;
      a[1][1] =   0.0;
      a[2][0] =   1.830769230769234;
      a[2][1] =   2.400000000000000;
      a[3][3] =   0.0;
      a[3][0] =   1.830769230769234;
      a[3][1] =   2.400000000000000;
      a[3][2] =   0.0;
      a[3][3] =   0.0;


      //\Gamma^{-1}
      c[0][0] =   1.0/gconst;
      c[1][0] =   1.246438746438751;
      c[1][1] =   1.0/gconst;
      c[2][0] =  -1.226780626780621e1; 
      c[2][1] =   4.266666666666667e1;
      c[2][2] =   1.0/gconst;
      c[3][0] =   5.824628046850726e-2;
      c[3][1] =   3.259259259259259;
      c[3][2] =  -3.703703703703704e-1;
      c[3][3] =   1.0/gconst;

      //sum^{i}_{j=1} \Gamma_{ij}
      gammaSum[0] = 0.3000000000000000;
      gammaSum[1] = 0.1878205128205124;
      gammaSum[2] =-1.0000000000000000;
      gammaSum[3] = 0.0000000000000000;
 
      //\vec b^{T}\Gamma^{-1}
      m[0] =   1.830769230769234;
      m[1] =   2.400000000000000;
      m[2] =   0.000000000000000;
      m[3] =   1.000000000000000;

      //\hat{\vec b}^{T}\Gamma^{-1} 
      mhat[0] =   2.214433650496747;
      mhat[1] =   1.831186394371970;
      mhat[2] =   8.264462809917363e-3;
      mhat[3] =   0.000000000000000;

      s[0][0] =   0.0;
      s[1][0] =   5.555555555555556;
      s[1][1] =   0.0;
      s[2][0] =  -4.239316239316217; 
      s[2][1] =   8.000000000000000;
      s[2][2] =   0.0;
      s[3][0] =  -4.239316239316217;
      s[3][1] =   8.000000000000000;
      s[3][2] =   0.0;
      s[3][3] =   0.0;

      //\sigma_i
      sigmaSum[0] = 0.0;
      sigmaSum[1] = 1.666666666666667;
      sigmaSum[2] = 2.307692307692341e-1;
      sigmaSum[3] = 2.307692307692341e-1;


    }
  else
    {
      std::cerr<<"PROBLEM RoseBase scheme= "<<schemeType<<" not implemented "
	       <<" quitting "<<std::endl;
      assert(0);

    }
  return false;
}

//======================================================================
//main routines
//======================================================================

bool RosenbrockIntegratorBase::calculateSolution(const real& tout,
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

      deltaTnp1 = calculateInitialDt(tn,tout);

#ifdef DEBUG_ROS
      std::cout<<"RoseBase calcSoln("<<tout<<") dt0 = "<<deltaTnp1<<std::endl;
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
      std::cout<<"In RoseBase calcSoln("<<tout<<") tn = "<<tn<<" dt = "
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



bool RosenbrockIntegratorBase::step(const real& tOut,
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
	  std::cerr<<"RoseBase step ("<<tn<<" --> "<<tOut
		   <<" takOneRosenbrockStep failed "<<std::endl;
	  stepFailed = roseStepFailed;
	  return stepFailed;
	}

      //determine if estimated truncation error is sufficiently low
      bool estFailed = calculateErrorEstimate(*ynp1,yerr,errEst,locErrCur);
      if (estFailed)
	{
	  std::cerr<<"RoseBase step ("<<tn<<" --> "<<tOut
		   <<" ErrEstp failed "<<std::endl;
	  stepFailed = estFailed;
	  return stepFailed;
	}
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
      std::cout<<"RoseBase step("<<tOut<<") solve to = "
	       <<tn+deltaTnp1<<std::endl;
      std::cout<<"locErrCur = "<<locErrCur<<" errorOk= "
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
	      std::cerr<<"RoseBase step("<<tOut<<") solve to = "
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
	      std::cerr<<"RoseBase step("<<tOut<<") solve to = "
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

	  
	  ynPrime   = *ynp1Prime;
	  yn        = *ynp1;
	      
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

bool RosenbrockIntegratorBase::calculateErrorEstimate(const Vec& y1,
							   const Vec& y2,
							   Vec& err,
							   real& D)
{
  /***************************************************

  calculate \vec{e}^{n+1} = 
    \vec{y}^{n+1}_{(1)}-\vec{y}^{n+1}_{(2)}
  and  \|\vec {e}^{n+1}\| using WRMS norm or other norm that's passed in

  Now follow RODAS lead and use max(|y_n|,|y^{n+1}|) in weighting
  **************************************************/
  assert(errNorm);
  tmpvec = y1;
  for (int i=0; i < y1.size(); i++)
    {
      if (std::fabs(yn[i]) > std::fabs(y1[i]))
	tmpvec[i] = yn[i];
    }
  errNorm->setWeight(tmpvec);


  err = y1;
  axpy(-1.0,y2,err);

  D = (*errNorm)(err);
#ifdef DEBUG_ERROR_CONTROLLERS
  //std::cout<<"RosBase calcErrEst D= "<<D<<" y1= "<<y1<<" y2= "
  //	   <<y2<<" err = "<<err<<std::endl;
  //   //mwf try a manual calculation ?
  //   real D2 = 0.0;
  //   int neq =  y1.size();
  //   for (int i=0; i < neq; i++)
  //     {
  //       real sk = absTol + 
  //          relTol*std::max(std::fabs(y1[i]),std::fabs(yn[i]));
  //       D2 += err[i]*err[i]/sk/sk;
  //     }
  //   D2 = std::sqrt(D2/neq);
  //   std::cout<<"RosBase rodas wgt D= "<<D2<<" wrms D= "<<D<<std::endl;
#endif

  bool errEstFailed = D < 0.0;
  return errEstFailed;
}




real RosenbrockIntegratorBase::calculateInitialDt(const real& t0,
						       const real& tout)
{
  //follow approach from Shampine_94 pg 379
  //assumes ynPrime is set
  assert(errNorm);
  errNorm->setWeight(yn);

  real normY0Prime = (*errNorm)(ynPrime);
  real denom = std::max(timeEps,normY0Prime);
  int errOrder  = errController->getOrder();
  real errOrdInv= 1.0/(errOrder+1.0);
  real numer = std::pow(timeTol,errOrdInv);
  real dt = std::min(dt0Max,std::min(tout-t0,numer/denom));
#ifdef DEBUG_ERROR_CONTROLLERS
  std::cout<<"RosenbrockBase initial dt = "<<dt
	   <<" normY0Prime= "<<normY0Prime<<std::endl;
#endif
  return dt;
}


bool RosenbrockIntegratorBase::checkStepSize(const real& t0, 
					     const real& tout,
					     real& dt)
{
  //assuming have an estimate for \delta t^{n+1} in 
  //dt then check to make sure not going over requested 
  //time level
  
  //follow algorithm 2.7 in Kavetsk_Binning_etal_01 sec 2.7
  bool stepChanged = false;
  if (t0 + dt >= tout)
    {
      dt = tout-t0;
      stepChanged = true;
    }
  else if (t0 + 2.0*dt >= tout)
    {
      dt = 0.5*(tout-t0);
      stepChanged = true;
    }
  if (stepChanged)
    errController->invalidateLastEstimate();

  bool checkFailed = false;
  checkFailed = t0 + dt > tout || (dt <= 1.0e-14 && t0 + dt < tout);

  return checkFailed;
}


real RosenbrockIntegratorBase::estimateStepSize(const real& err,
						const real& dtLast,
						bool& stepOk)
{
  //use passed in error Controller now
  real dtOut(-1.0);

  dtOut = errController->estimateStepSize(err,dtLast,stepOk);

  return dtOut;
}


bool RosenbrockIntegratorBase::calculateDFDt(const real& tCur,
						  const Vec& yCur,
						  const Vec& Fcur,
						  Vec& DFDt)
{
  //calculate F_t using simple numerical differencing
  //F_t \approx (F(t^n+\delta,\vec y^n)-F(t^n,\vec y))/\delta

  real delta = 1.0e-6; // figure out a good value for this
  bool evalFailed = dae->rightHandSideValue(tCur+delta,yCur,DFDt);

  if (evalFailed)
    return evalFailed;
  

  //F(t+\delta,y)-F(t,y)
  axpy(-1.0,Fcur,DFDt);

  //[F(t+\delta,y)-F(t,y)]/\delta
  scal(1.0/delta,DFDt);

  return false;

}

  //follow approach suggested in Lang_01
void RosenbrockIntegratorBase::computeDeltaForJacobian()
{
  dae->computeOwnJacobianDelta = false;
  dae->deltaVF = 0.0;
  for (int i=0;i<dae->yDaeDef.ldim();i++)
    {
      real ysign = 1.0;
      if (dae->yDaeDef[i] < 0.0)
	ysign = -1.0;
      dae->deltaVF[i] = SQRT_MACHINE_EPSILON*
	ysign*std::max(std::fabs(dae->yDaeDef[i]),absTol);
    }
    
}

}//end Daekt
