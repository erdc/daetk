#ifndef ROSENBROCK_INTEGRATOR_BASE_H
#define ROSENBROCK_INTEGRATOR_BASE_H

#include <vector>
#include "Definitions.h"
#include "Vec.h"

#include "Integrator.h"
#include "DataCollector.h"
#include "VectorNorm.h"

#include "RosenbrockDaeDefinition.h"
#include "TimeIntegrationErrorController.h"

namespace Daetk
{
//forward declarations

class RosenbrockIntegratorBase : public Integrator
{
  /***********************************************************************

    The purpose of this class is to implement the core of a
    Rosenbrock time integration scheme for ODEs as formulated in
    \cite{Hairer_Wanner_96,Lang_01} 

    \mat{M}\dot{\vec y} &=& \vec {F}(t,\vec y)

    where \mat{M} is constant but potentially singular.

    The basic quantities are (see notes)

    \vec y^{n+1} &=& \vec y^{n} + \Delta t\sum^{s}_{i=1}m_i\vec Y_{ni}

    \hat{\vec y}^{n+1} &=& \vec y^{n} + 
                           \Delta t\sum^{s}_{i=1}\hat{b}_i\vec Y_{ni} 

    with stage values given by

    (\frac{1}{\gamma\Delta t}\mat{M} - \mat{F}_{y}(t^n,\vec y^n))
      \cdot\vec{Y}_{ni} &=& \vec F(t_i,\vec U_i) - 
                 \mat{M}\sum_{j=1}^{i-1}\frac{c_{ij}}{\Delta t}\vec Y_{nj} +
                 \gamma_i\Delta t\vec{F}_{t}(t^n,\vec y^n)
  
    and internal values

    t_i = t^n + \alpha_i\Delta t

    \vec Y_i = \vec y^n + \sum^{i-1}_{j=1}a_{ij}\vec Y_{nj}


    At each time step, the code requires a Jacobian evaluation 
      \vec {F}_{y} = \pd{F}{\vec y}
    and
      \vec {F}_{t} = \pd{F}_{t}
    where the terms are evaluated at t^n and \vec y^n

    Right now, F_{y} can be generated numerically, but the user 
    has to provide F_{t} until I write a numerical difference 
    operator 

    The integrator by default uses a PI controller as described by
    Gustafsson_94 for time step selection
 
    Right now templated on class that is used to represent matrix
    for 1/(\gamma\Delta t)\mat{M} -\vec{F}_{y}.
  
  **********************************************************************/
public:

  enum RosenbrockScheme {ROS2, ROS3P, RODAS3, RODAS4, RODASP, ROWDAIND2};

  RosenbrockIntegratorBase(RosenbrockDaeDefinition* daeIn, 
			   DataCollector* dataIn,
			   VectorNorm* errNormIn,
			   TimeIntegrationErrorController* errControllerIn,
			   const RosenbrockScheme& schemeIn, 
			   real absTolIn, 
			   real relTolIn,
			   real dt0MaxIn);

  virtual ~RosenbrockIntegratorBase();

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

  //---------- specific to Rosenbrock Integrators ----------
  virtual bool usesAuxiliaryVariable() const
  { return false; }
  virtual bool getAuxiliaryVariable(Vec& aux) 
  { return true; }
protected:

  //---------- Data Members ----------
  //user-defined problem supposed to be solving
  RosenbrockDaeDefinition* dae;
  //collect statistics about simulation
  DataCollector* data;
  //norm to use in error calculations
  VectorNorm* errNorm;
  //which Rosenbrock scheme to use
  RosenbrockScheme scheme;
  //number of stages and approximation order of scheme
  int numStages;
  int approxOrder;

  //tolerances in local truncation error estimate
  real absTol; //\eps_{a} 
  real relTol; //\eps_{r}

  real dt0Max; //\Delta t^{0}_{min}

  //approximate solutions
  Vec* ynp1;       //\vec{y}^{n+1}, storage taken from user
  Vec* ynp1Prime;  //\dot{\vec{y}}^{n+1}, storage taken from user

  //--- stage values ---
  //Y_{ni} see Lang_01 V.5
  std::vector<Vec> Yni;

  real deltaTnp1;   //\Delta t^{n+1}_{j+1}

  //F_t at t^n,\vec y^n
  Vec Jt;

  //history
  Vec yn;           //\vec{y}^{n}
  Vec ynPrime;      //\dot{\vec{y}}^{n}, not strictly necessary
  //lower order embedded solution for error estimation
  Vec yerr;         //\hat{\vec y}^{n+1}

  real tn;          //t^{n}
  real deltaTn;     //\Delta t^{n}

  //for step selection 
  real locErrCur;        //D^{c} = \|\vec e^{n+1} \| for current   step
  Vec errEst;            //\vec e^{n+1}

  TimeIntegrationErrorController* errController;

  //data necessary for building approximation
  bool firstStep;         //is history valid?

  real timeTol;     //\tau = 1  for relative norm tests
  real timeEps;     //\eps for making sure no divide by zero

  int maxErrorFailures;         //maximum number of trucation error failures
                                //per step
  int maxFailures;              //maximum failures of any type

  //work arrays
  Vec tmpvec,tmpprod,rhs;

  //---------- coefficients defining actual Rosenbrock scheme ----------
  //           see Lang_01 Appendix C
  static const int MAXSTAGES  = 6;

  //defined in Lang_01 pg 48--50
  real gamma[MAXSTAGES][MAXSTAGES]; //\Gamma
  real a[MAXSTAGES][MAXSTAGES];     //\alpha_{ij} \Gamma^{-1}
  real c[MAXSTAGES][MAXSTAGES];     //\Gamma^{-1}
  real gammaSum[MAXSTAGES];         //\gamma_i
  real alphaSum[MAXSTAGES];         //\alpha_i
  real m[MAXSTAGES];                //b^{T}\Gamma^{-1}
  real mhat[MAXSTAGES];             //\hat{b}^{T}\Gamma^{-1}

  //necessary for implicit ode formulation
  real s[MAXSTAGES][MAXSTAGES];    //\sigma\Gamma^{-1}
  real sigmaSum[MAXSTAGES];        //\sigma_i

  //---------- Basic Functions For Algorithm ----------
  //sets numStages, approxOrder, gamma, c, gammaSum, m, mhat 
  bool setRosenbrockCoefficients(const RosenbrockScheme& schemeType);

  //go from t^{0} to t^{0} + \Delta t starting with values y^0, yp^0
  virtual
  bool takeOneRosenbrockStep(const real& tIn,
			     const Vec& yIn,
			     const Vec& yInPrime,
			     const real& dt,
			     Vec& yOut,
			     Vec& yErrOut) = 0;

  //calculate \vec{e}^{n+1} = 
  //  \vec{y}^{n+1}_{(1)}-\vec{y}^{n+1}_{(2)}
  //and  \|\vec {e}^{n+1}\| 
  virtual bool calculateErrorEstimate(const Vec& y1, const Vec& y2,
				      Vec& err, real& D);

  //estimate step size using PI basic controller from
  //Gustafsson_94 section 5
  //sets values of lastStepAccepted and locErrLastAcc, or locErrLastRej
  //and deltaTlastAcc or deltaTlastRej
  virtual real estimateStepSize(const real& err, const real& dtLast, 
				bool& stepOk);

  //look ahead algorithm to make sure \Delta t^{n+1}_{j+1} 
  //doesn't step over t^{out}
  virtual bool checkStepSize(const real& t0, const real& tout, real& dt);
  //pick \Delta t^{0} following standard approach from Shampine_94
  virtual real calculateInitialDt(const real& t0, const real& tout);

  //calculate F_t using simple numerical differencing
  //F_t \approx (F(t^n+\delta,\vec y^n)-F(t^n,\vec y))/\delta
  virtual bool calculateDFDt(const real& tCur, const Vec& yCur,
			     const Vec& Fcur,
			     Vec& DFDt);

  //compute increment in numerical Jacobian, can also let
  //user code do this
  virtual void computeDeltaForJacobian();

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

  

};//end RosenbrockIntegratorBase



}//Daetk

#endif
