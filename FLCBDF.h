#ifndef FLCBDF_H
#define FLCBDF_H

#include <fstream>
#include <iostream>
#include "Definitions.h"
#include "Integrator.h"
#include "Utilities.h"
#include "DaeDefinition.h"
#include "LinearSolver.h"
#include "ModifiedNewton.h"
#include "Jacobian.h"
#include "VectorNorm.h"
#include "DataCollector.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"
#include "Chronograph.h"

namespace Daetk 
{
class FLCBDF : public Integrator
{
public:
  FLCBDF();

  virtual ~FLCBDF();

  FLCBDF(DaeDefinition&,LinearSolver&,NonlinearSolver&,Jacobian&,
         VectorNorm&,DataCollector&);

  bool takeInitialSteps(const real& tout,real& tStep,Vec& solutionAtTStep, 
			Vec& sp);
  bool takeInitialStepsOriginal(const real& tout,real& tStep,Vec& solutionAtTStep, 
				Vec& sp);
  bool takeInitialStepsWithRichExtrap(const real& tout,real& tStep,Vec& solutionAtTStep, 
				      Vec& sp);
  bool takeInitialStepsWithModifiedCoefs(const real& tout,real& tStep,Vec& solutionAtTStep, 
					 Vec& sp);
  
  bool step(const real& tout,real& tStep,Vec& solutionAtTStep, Vec& sp);
  bool stepSpecial(const real& tout,real& tStep,Vec& solutionAtTStep, Vec& sp);

  bool calculateSolution(const real& tout,Vec& solutionAtTout, Vec& sp);
  bool calculateSteadyStateSolution(Vec& solutionAtTout, Vec& sp, real sstol);

  void useFixedStep(const real step);
  void useFixedOrder(const int order);
  void updateJacobians();
  void useInterpolant();
  void reset();
  //for testing startup procedures and step size choice
  void setStepSizeProcedure(const int& stepSizeFlag);
  void setStepIncreaseFloor(const real& floorIn);
  void setStepIncreaseCeiling(const real& ceilIn);
  void setStepDecreaseFloor(const real& floorIn);
  void setStepDecreaseCeiling(const real& ceilIn);

  void setStartupProcedure(const int& startFlag);
  void setStepIncreaseFloorForStartup(const real& floorIn);
  void setStepIncreaseCeilingForStartup(const real& ceilIn);
  void setStepDecreaseFloorForStartup(const real& floorIn);
  void setStepDecreaseCeilingForStartup(const real& ceilIn);
  void setRichardsonExtrapolationOrder(const int& richExtrapOrder);

private:
  
  bool prepareJacobian,                  
    predictorAlreadyUpdated,                      
    firstCallToSolver,
    jacobianIsOld,
    correctorStepFailed,
    kRaisedOnLastStep,
    USE_FIXED_STEP,
    USE_FIXED_ORDER,
    UPDATE_JACOBIANS,
    STEP_EXACT;

  int neq,                               //number of equtions    
    k,                                   //order
    kLast,
    constStepsTaken,
    kLastUpdate,
    constStepsTakenLastUpdate,
    stepsFailed,
    errorFailures,
    MAX_FAILURES,
    shortCutFactor,
    kFixed;
  
  real h,                                 //stepsize
    hLast,
    hLastUpdate,
    tn,                                   //time of previous step
    tnPlusH,
    error,
    errorEstimate,
    errorRichExt,			 //seg for Richardson Extrapolation
    norm_yCminusyP,
    convergenceFactor,                    //for nonlinear solver
    hFixed;
							    
  Vec *yn,*ynprime,                         //solution values at tn
    yCminusyP,tempvec;

  //coefficient data

  real psiNp1[7],
    alphaNp1[6],
    betaNp1[6],
    sigmaNp1[7],
    gammaNp1[6],
    alphaS,
    alphaoNp1,
    alphaOld,
    alphaLast;
  //data members for startup routines
  int startupType; // 0 original
                   // 1 original routine with user dependent
                   //   floor/ceiling on step size growth
                   // 2 Richardson Extrapolation
  //bounds for ratio used in step size increase/decrease for startup
  real stepIncreaseCeilingStartup,
    stepIncreaseFloorStartup,
    stepDecreaseCeilingStartup,
    stepDecreaseFloorStartup;
  //mwf add for Richardson extrapolation configuration
  int forceRichardsonExtrapolationOrder;
  int chooseStepSizeFlag; // 0 original
                          // 1 original routine with user dependent
                          //   floor/ceiling on step size growth
  //bounds for ratio used in step size increase/decrease for general running
  real stepIncreaseCeiling,
    stepIncreaseFloor,
    stepDecreaseCeiling,
    stepDecreaseFloor;

  // seg--timing each step during startup routine 
  Chronograph clock;
  real stepRunTime;
  
  real norm_yByS;			// variables for extrapolation norm
  Vec  tmp,
    yBigStep,ySmallStep,		  // temporary vectors for initialSteps routine -- seg
    yExt,ySminusyB;

  Vec phiN[7];
  
  DaeDefinition* dae;

  DataCollector* data;

  LinearSolver* linearSolver;

  NonlinearSolver* nonlinearSolver;

  Jacobian* jacobian;
 
  VectorNorm* weightedNorm;

  //private member functions 

  bool analyzeJacobian();
  bool evaluateJacobian();
  void computeDeltaForJacobian();
  void resetCoefficients();
  void updateCoefficients();
  void chooseOrderForStep();
  void predictor();
  void  interpolant(const real& T,Vec& yAtT,Vec& yPAtT);
  void chooseStepSize();
  void chooseStepSizeOriginal();
  //mwf added so that I can try and eliminate need to interpolate to
  //mwf time for final soluiton
  void chooseStepSize(const real& tin, const real& tout);
  void chooseStepSizeOriginal(const real& tin, const real& tout);
  void chooseStepSizeStartup(const real& tin, const real& tout);

  bool corrector();
  bool errorForStepTooLarge();
  void checkStepSize();

};
}//Daetk
#endif
