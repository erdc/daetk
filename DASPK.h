#ifndef DASPK_H
#define DASPK_H

#include "Definitions.h"
#include "Integrator.h"
#include "DaeDefinition.h"
#include "DataCollector.h"
#include "Jacobian.h"
#include "Preconditioner.h"
#include "DaspkRes.h"
#include "DaspkJac.h"
#include "DaspkPsol.h"

#include "Definitions.h"
#include "Vec.h"
#include "IntVec.h"



namespace Daetk 
{
class DASPK : public Integrator
{
public:

  DASPK();
  virtual ~DASPK(){}
  DASPK(DaeDefinition&,DataCollector&);
  bool calculateSolution(const real& tout,Vec& solution,Vec& solutionPrime);
  bool step(const real& tOut,real& tStep,Vec& solutionAtTStep, Vec& sP);
  bool computeInitialConditions();

  void setTolerances(real& atol,real& rtol);
  void setTolerances(Vec& atol,Vec& rtol);
  void setTstop(const real& tStop);
  void useAnalyticJacobian(Jacobian& J);
  void usePreconditioner(Preconditioner& P);
  void useFullLinearSolver();
  void useBandedLinearSolver(int kl,int ku);
  void useKrylovSolver();
  void useKrylovSolver(int maxl,int kmp,int nrmax,real& epli);
  void setMaxStepSize(const real& tmax);
  void setInitialStepSize(const real& h0);
  void setMaxOrder(int maxOrder);
  void setInitialConstraints(const IntVec& constraintFlag);
  void enforceNonnegativity();
  void findConsistentY0();
  void findConsistentY0andY0Prime(const IntVec& variableFlag);
  void excludeAlgebraicVariables( const IntVec& variableFlag);
  void setIntialConditionCalcParms(int mxnit,int mxnj,int mxnh,int lsoff, 
                                   const real&  stptol, const real& epinit); 
  void useExtraPrinting(int printLevel);
  void reset();
  void resetTstop();
  void init();
private:

  void collectData();
  int nreOld,
    njeOld,
    netfOld,
    ncfnOld,
    ncflOld,
    nniOld,
    nliOld;
    

  int idid,
    kl,
    ku,
    maxord,
    neq,
    lrw,
    liw,
    maxl,
    kmp,
    nrmax,
    mxnit,
    mxnj,
    mxnh,
    lsoff,
    lenwp,
    leniwp;

  real *atolp,
    *rtolp,
    tstop,
    epli,
    tmax,
    h0,
    stptol,
    epinit,
    t;

  IntVec ipar,
    info,
    iwork,
    cFlag,
    avFlag;
  
 
  Vec *vecAtol,
    *vecRtol,
    rpar,
    rwork,
    y,
    yp;

  DaeDefinition* dae;
  DataCollector* data;
  Jacobian* jacobian;
  Preconditioner* preconditioner;
   
  DaspkResidual DR;
  DaspkJacobian DJ;
  DaspkPsol  DP;
};

typedef void (*GlobDaspkRes)(const real& t,real* y,real* yp,
                         real& cj,real* delta, int* ires, real* rpar, 
                         int *ipar);


typedef void (*GlobDaspkPsol)(const int& neq, const real& t, real* y, 
                          real* yp, real* savr,real* wk,
                          const real& cj,real* wght,real * wp,
                          int* iwp, real* b, const real& eplin, int* ier, 
                          real* rpar,int* ipar);

#ifndef CRAYCC
extern "C" void F77NAME(ddaspk)(GlobDaspkRes res, const int& neq, real& t, 
                  real* y, real* yp, const real& tout, 
                  int* info, real* rtol, real* atol, 
                  int& idid, real* rwork, int& lrw, int* iwork, 
                  int& liw, real* rpar, int* ipar, void* jac, 
                  GlobDaspkPsol psol);
extern "C" void F77NAME(sdaspk)(GlobDaspkRes res, const int& neq, real& t, 
                  real* y, real* yp, const real& tout, 
                  int* info, real* rtol, real* atol, 
                  int& idid, real* rwork, int& lrw, int* iwork, 
                  int& liw, real* rpar, int* ipar, void* jac,
                  GlobDaspkPsol psol);
#else
extern "C" void F77NAME(DDASPK)(GlobDaspkRes res, const int& neq, real& t, 
                  real* y, real* yp, const real& tout, 
                  int* info, real* rtol, real* atol, 
                  int& idid, real* rwork, int& lrw, int* iwork, 
                  int& liw, real* rpar, int* ipar, void* jac,
                  GlobDaspkPsol psol);
extern "C" void F77NAME(SDASPK)(GlobDaspkRes res, const int& neq, real& t, 
                  real* y, real* yp, const real& tout, 
                  int* info, real* rtol, real* atol, 
                  int& idid, real* rwork, int& lrw, int* iwork, 
                  int& liw, real* rpar, int* ipar, void* jac,
                  GlobDaspkPsol psol);
#endif

}//Daetk
#endif








