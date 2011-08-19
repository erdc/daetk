#ifndef TIME_INTEGRATION_ERROR_CONTROLLER_H
#define TIME_INTEGRATION_ERROR_CONTROLLER_H

#include "Definitions.h"

namespace Daetk
{


  /***********************************************************************

    Try to extract basic interface for some time integration error 
    controllers: 

      "standard" error controller with basic ceiling and floor on
         step ratio

      PC controller as described in Gustafsson_94

    This version will be responsible for determining if a step is
     accepted or not
  **********************************************************************/
class TimeIntegrationErrorController
{
  /***********************************************************************

  **********************************************************************/
public:

  virtual ~TimeIntegrationErrorController() {}

  //if need to keep track of history. E.g., if need to know that
  //current step is the first step
  virtual bool reset() = 0;

  virtual real estimateStepSize(const real& err,const real& dtLast,
				bool& stepOk) = 0;

  virtual int getOrder() const = 0;
  virtual void setOrder(const int& ordIn) = 0;
  virtual real getTolerance() const = 0;

  //tell controller that last estimate was modified
  virtual bool invalidateLastEstimate() = 0;

  virtual int getMaxFailures() const = 0;
};

class DefaultTimeIntegrationErrorController : 
    public TimeIntegrationErrorController
{
  /***********************************************************************

    simple time integration error controller. Basic
    function provided is

    
    real estimateStepSize(const real& err, const real& dtLast, 
	                  const bool& stepOk);


   Does not handle initial step size estimate, if user wants to do 
   something different there.

   Default strategy is one described by Kavetski_Binning_etal_02

  **********************************************************************/
public:
  DefaultTimeIntegrationErrorController(const real& tolIn    = 1.0,
					const int& orderIn   = 1,
					const real& dtRmaxIn = 4.0,
					const real& dtRminIn = 0.1,
					const real& safetyIn = 0.8,
					const real& timeEpsIn= 1.e-10);

  virtual ~DefaultTimeIntegrationErrorController();

  //if need to keep track of history. E.g., if need to know that
  //current step is the first step
  virtual bool reset();

  virtual real estimateStepSize(const real& err,const real& dtLast,
				bool& stepOk);

  virtual int getOrder() const
  { return order; }
  virtual void setOrder(const int& ordIn)
  { order = ordIn; }
  virtual real getTolerance() const
  { return timeTol; }

  //no history in default estimate
  virtual bool invalidateLastEstimate() 
  {return false;}

  virtual int getMaxFailures() const
  { return 12; }

protected:
  /***********************************************************************
    timeTol     --- tolerance for integration 
    order       --- consistency order for estimate
    dtRmax      --- max allowed increase in step size
    dtRmin      --- min allowed decrease in step size
    safety      --- safety factor for truncation error estimate
    timeEps     --- bottom for error estimate to keep from x/0
  **********************************************************************/
  real timeTol;
  int order;
  real dtRmax,dtRmin,safety,timeEps;

};


class RODASerrorController: public DefaultTimeIntegrationErrorController
{
  /***********************************************************************

    Try strategy used in RODAS code see (Hairer_Wanner_96 pg 111) for
    basic strategy. Also mix with predictive controller from
    Gustafsson as is done in actual RODAS code (IWORK(3) = 0 or 1)

    
    
  **********************************************************************/
public:
  RODASerrorController(const bool& usePredControlIn = true,
		       const real& tolIn    = 1.0,
		       const int& orderIn   = 1,
		       const real& dtRmaxIn = 4.0,
		       const real& dtRminIn = 0.1,
		       const real& safetyIn = 0.8,
		       const real& timeEpsIn= 1.e-10);

  virtual ~RODASerrorController();

  //if need to keep track of history. E.g., if need to know that
  //current step is the first step
  virtual bool reset();

  virtual real estimateStepSize(const real& err,const real& dtLast,
				bool& stepOk);
  //no history in default estimate
  virtual bool invalidateLastEstimate() 
  { lastEstimateGood = false; return false;}

protected:
  /***********************************************************************
    usePred    --- try Gustafsson based estimate
    numAccept  --- number of accepted steps
    numReject  --- number of rejected steps
    dtLastAcc  --- \Delta t for last accepted step
    errLastAcc --- truncation error estimate for last step
    lastEstimateGood --- did the integrator use the last estimate 
                         i.e., is the history good
  **********************************************************************/
  bool usePred;
  int numAccept, numReject;
  real dtLastAcc,errLastAcc;
  bool lastEstimateGood;
}; //end RODASerrorController


class SoederlindErrorController: public DefaultTimeIntegrationErrorController
{
  /***********************************************************************

    Try strategy used in Soederlind_Wang_06 
    
    
  **********************************************************************/
public:
  SoederlindErrorController(const real& tolIn    = 1.0,
			    const int& orderIn   = 1,
			    const real& dtRmaxIn = 4.0,
			    const real& dtRminIn = 0.1,
			    const real& safetyIn = 0.8,
			    const real& timeEpsIn= 1.e-10);

  virtual ~SoederlindErrorController();

  //if need to keep track of history. E.g., if need to know that
  //current step is the first step
  virtual bool reset();

  virtual real estimateStepSize(const real& err,const real& dtLast,
				bool& stepOk);
  //no history in default estimate
  virtual bool invalidateLastEstimate() 
  { lastEstimateGood = false; return false;}

  virtual int getMaxFailures() const
  { return 100; }

protected:
  /***********************************************************************
    dtPrevCall  --- \Delta t for previous call
    errPrevCall --- truncation error estimate for last step
    lastEstimateGood --- did the integrator use the last estimate 
                         i.e., is the history good
  **********************************************************************/
  real dtPrevCall,errPrevCall;
  bool lastEstimateGood;
  //mwf added for filter
  real limitRatio(const real& ratIn);
  real filterH211B(const real& err, const real& errOld, const real& ratOld,
		   const int& order);
  real filterPI42(const real& err, const real& errOld, const real& ratOld,
		   const int& order);
  real filterSTD(const real& err, const real& errOld, const real& ratOld,
		 const int& order);
}; //end SoederlindErrorController

}//Daetk

#endif
