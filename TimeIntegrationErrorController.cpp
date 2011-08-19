#include "TimeIntegrationErrorController.h"

#include <cassert>
#include <iostream>
#include <cmath>

namespace Daetk
{

////////////////////////////////////////////////////////////////////////
//DefaultTimeIntegrationErrorController
////////////////////////////////////////////////////////////////////////

DefaultTimeIntegrationErrorController::
DefaultTimeIntegrationErrorController(const real& timeTolIn,
				      const int& orderIn,
				      const real& dtRmaxIn,
				      const real& dtRminIn,
				      const real& safetyIn,
				      const real& timeEpsIn):
  timeTol(timeTolIn),
  order(orderIn),dtRmax(dtRmaxIn),dtRmin(dtRminIn),
  safety(safetyIn),timeEps(timeEpsIn)
{

}

DefaultTimeIntegrationErrorController::
~DefaultTimeIntegrationErrorController()
{}

bool DefaultTimeIntegrationErrorController::reset()
{
  return false;
}

real 
DefaultTimeIntegrationErrorController::estimateStepSize(const real& err,
							const real& dtLast,
							bool& stepOk)
{

  //estimate step size using basic controller from
  //Kavetski_Binning_etal_02 (23--24) essentially the same as
  //Kavetski_Binning_etal_01 (15--16)
  //
  //assumes that err is the local trunction error
  //thas already been calculated for corresponding
  //solution
  //
  //output can either be \Delta t^{n+2}_{1} or \Delta t^{n+1}_{j+1}
  //depending on whether or not step was accepted
  stepOk = err < timeTol;

  real dtOut(-1.0);
  real ordInv = 1.0/(order+1.0);
  real minErr = std::max(err,timeEps);
  real term1 = safety*std::pow(timeTol/minErr,ordInv);

  if (stepOk) //last step was a successful one 
    dtOut = dtLast*std::min(term1,dtRmax);
  else
    dtOut = dtLast*std::max(term1,dtRmin);

  return dtOut;

}//end est
////////////////////////////////////////////////////////////////////////
//RODASerrorController
////////////////////////////////////////////////////////////////////////

RODASerrorController::
RODASerrorController(const bool& usePredControlIn,
		     const real& timeTolIn,
		     const int& orderIn,
		     const real& dtRmaxIn,
		     const real& dtRminIn,
		     const real& safetyIn,
		     const real& timeEpsIn):
  DefaultTimeIntegrationErrorController(timeTolIn,orderIn,dtRmaxIn,
					dtRminIn,safetyIn,timeEpsIn),
  usePred(usePredControlIn),numAccept(0),numReject(0),
  dtLastAcc(-1.0),errLastAcc(-1.0),lastEstimateGood(false)
{

}

RODASerrorController::~RODASerrorController()
{}

bool RODASerrorController::reset()
{
  numAccept  = 0;
  numReject  = 0;
  errLastAcc = -1.0;
  dtLastAcc  = -1.0;
  lastEstimateGood = false;
  return false;
}

real 
RODASerrorController::estimateStepSize(const real& err,
				       const real& dtLast,
				       bool& stepOk)
{
  /***********************************************************************
    Follow approach used in RODAS code. Note this is in terms of
    r= dtOld/dtNew rather than r = dtNew/dtOld

    Get standard estimate with usual floor and ceiling by default.

         r_std = max(1/r_{max},min(1/r_min,D^{1/p+1})/s)

    if (step accepted)
       if (not the first step)

          calculate predictive controller estimate for new step
             r_p = (dt_a/dt)*(D^2/D_a)^{1/p+1}/s
          pick estimate to be min of predictive value and classical value 
             r   = max(r_std,r_p)
             dt_a= dt
             D_a = max(.01,D)
       else
          r = r_std
       endif
    else

       r = r_std

    endif   

    dt^{+} = dt/r
  **********************************************************************/
  stepOk = err < timeTol;
  //1/(p+1)
  real orderInv = 1.0/(order+1.0);
  real dtRminInv = 1.0/dtRmin;
  real dtRmaxInv = 1.0/dtRmax;
  //standard estimate
  real term0 = std::pow(err,orderInv)/safety;
  real term1 = std::min(term0,dtRminInv);
  real rstd  = std::max(dtRmaxInv,term1);
  //by default use standard approx
#ifdef DEBUG_ERROR_CONTROLLERS
  std::cout<<"Rodas contr est err= "<<err
	   <<"\n\t std ratio= "<<rstd<<" dtStd= "<<dtLast/rstd<<std::endl;
#endif
  real rout  = rstd;
  if (stepOk) //accepted
    {
      numAccept++;
      if (usePred)//try Gustafssson
	{
	  if (numAccept > 1 && lastEstimateGood) //not first step
	    {
	      assert(dtLastAcc > 0.0);
	      assert(errLastAcc > 0.0);
	      real term0 = err*err/errLastAcc;
	      real term1 = std::pow(term0,orderInv)/safety;
	      real rpred = dtLastAcc/dtLast*term1;
	      rout       = std::max(rpred,rstd);
#ifdef DEBUG_ERROR_CONTROLLERS
	      std::cout<<"\n\t pred: errRat= "<<term0
		       <<" errPow= "<<term1
		       <<" dtLastAcc= "<<dtLastAcc<<" dtLast= "<<dtLast
		       <<"\n\t dtRat= "<<dtLastAcc/dtLast
		       <<" rpred= "<<rpred<<std::endl;
#endif
	    } //end not first step
	  dtLastAcc = dtLast;
	  errLastAcc= std::max(0.01,err);
	}//end try pred
    }//end accepted
  else
    numReject++;


  //rout should be set correctly now
  real dtOut =  dtLast/rout;

  lastEstimateGood = true;

#ifdef DEBUG_ERROR_CONTROLLERS
  std::cout<<"\n\t chose dt= "<<dtOut<<std::endl;
#endif
  return dtOut;
}//end est

////////////////////////////////////////////////////////////////////////
//SoederlindErrorController
////////////////////////////////////////////////////////////////////////

SoederlindErrorController::
SoederlindErrorController(const real& timeTolIn,
			  const int& orderIn,
			  const real& dtRmaxIn,
			  const real& dtRminIn,
			  const real& safetyIn,
			  const real& timeEpsIn):
  DefaultTimeIntegrationErrorController(timeTolIn,orderIn,dtRmaxIn,
					dtRminIn,safetyIn,timeEpsIn),
  dtPrevCall(-1.0),errPrevCall(-1.0),lastEstimateGood(false)
{

}

SoederlindErrorController::~SoederlindErrorController()
{}

bool SoederlindErrorController::reset()
{
  dtPrevCall = -1.0;
  errPrevCall= -1.0;
  lastEstimateGood = false;
  return false;
}

real 
SoederlindErrorController::estimateStepSize(const real& err,
					    const real& dtLast,
					    bool& stepOk)
{
  /************************************************************************
     See Soederlind_Wang_06
   ************************************************************************/
  //1/(p+1)
  int errOrder = order+1;
  
  real errRatPrev     = timeTol/errPrevCall;
  real errRat         = timeTol/err;
  real stepRatioPrev  = dtLast/dtPrevCall;

  real stepRatTmp = -1.0;
  bool useSecondOrderEst = lastEstimateGood;
#ifdef FORCE_FIRST_ORDER_IN_SOEDERLIND
  useSecondOrderEst = false;
#endif
  //mwf debug useSecondOrderEst = false;
  if (useSecondOrderEst)
    stepRatTmp = filterH211B(errRat,errRatPrev,stepRatioPrev,errOrder);
  else
    stepRatTmp = filterSTD(errRat,errRatPrev,stepRatioPrev,errOrder);
  assert(stepRatTmp > 0.0);

  real stepRatLim = limitRatio(stepRatTmp);
  //should I put safe guards on ratio?
#ifdef USE_SAFEGUARD_IN_SOEDERLIND
  real stepRat = std::max(dtRmin,std::min(stepRatLim,dtRmax));
#else
  real stepRat = stepRatLim;
#endif
  real dtOut  = dtLast*stepRat;
#ifdef DEBUG_ERROR_CONTROLLERS
  std::cout<<"Soeder err= "<<err<<" errPrev= "<<errPrevCall
	   <<"\n\t stepRatPrev= "<<stepRatioPrev
	   <<" new unlim = "<<stepRatTmp<<" new lim= "<<stepRatLim
	   <<"\n\t stepRat= "<<stepRat<<" dtNew = "<<dtOut<<std::endl;
#endif
  //are these updated for failed steps?
  errPrevCall = err;
  dtPrevCall  = dtLast;

  lastEstimateGood = true;

  stepOk = true;
  if (stepRat < 0.9)
    stepOk = false;

  return dtOut;
}//end est
real SoederlindErrorController::limitRatio(const real& ratIn)
{
  real kappa = 1.0;
  real limit = 1.0 + kappa*std::atan((ratIn-1.0)/kappa);
  return limit;
}

real SoederlindErrorController::filterH211B(const real& err, const real& errOld, 
					    const real& ratOld,
					    const int& order)
{
  real b = 4.0;
  assert(b*order > 0.0);
  real bkinv = 1.0/(b*order);
  real binv = 1.0/b;

  real h211b = std::pow(err,bkinv)*std::pow(errOld,bkinv)*std::pow(ratOld,-binv);
  return h211b;
}
real SoederlindErrorController::filterPI42(const real& err, const real& errOld, 
					   const real& ratOld,
					   const int& order)
{
  real k35inv = 3.0/(5.0*order);
  real k15inv = 1.0/(5.0*order);
  real pi42   = std::pow(err,k35inv)*std::pow(errOld,-k15inv);
  return pi42;

}
real SoederlindErrorController::filterSTD(const real& err, const real& errOld, 
					  const real& ratOld,
					  const int& order)
{
  assert(order > 0);
  real ordInv = 1.0/order;

  real stdFilter = std::pow(err,ordInv);
  return stdFilter;
}


}//end Daetk
