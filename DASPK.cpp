#include "DASPK.h"

namespace Daetk 
{
using std::cerr;
using std::endl;

DASPK::DASPK():
  maxord(5),
  neq(0),
  nrmax(5),
  lsoff(0),
  lenwp(0),
  leniwp(0),
  atolp(0),
  rtolp(0),
  epli(0.05), 
  stptol(pow(MACHINE_EPSILON,2.0/3.0)),
  epinit(0.01),
  info(20,0,1),
  cFlag(0),
  avFlag(0),
  dae(0),
  data(0),
  jacobian(0),
  preconditioner(0)
{}


DASPK::DASPK(DaeDefinition& daeIn, DataCollector& dataIn):
  maxord(5),
  neq(daeIn.getY0().dim()),
  maxl(std::min(5,int(daeIn.getY0().dim()))),
  kmp(std::min(5,int(daeIn.getY0().dim()))),
  nrmax(5),
  lsoff(0),
  lenwp(0),
  leniwp(0),
  atolp(0),
  rtolp(0),
  epli(0.05), 
  stptol(pow(MACHINE_EPSILON,2.0/3.0)),
  epinit(0.01),
  info(20,0,1),
  cFlag(daeIn.getY0().dim()),
  avFlag(daeIn.getY0().dim()),
  y(daeIn.getY0().dim()),
  yp(daeIn.getY0().dim()),
  dae(&daeIn),
  data(&dataIn),
  jacobian(0),
  preconditioner(0)
{}
  

bool DASPK::calculateSolution(const real& tout,Vec& solution,
                              Vec& solutionPrime)
{
  info(3) = 1;
  if (info(1) == 0) // if first call
    {
      t = dae->getT0();
      y = dae->getY0();
      yp = dae->getY0prime();
    }

  DaspkResidual::theDae = dae;
  DaspkJacobian::theDae = dae;
  DaspkJacobian::thePrec = preconditioner;
  DaspkJacobian::theJac = jacobian;
  DaspkPsol::thePrec  = preconditioner;
  data->startUserStep();
  while (t < tout && idid >= -1)
    {
      if (info(12) == 0)
        {
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
          F77NAME(ddaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#else
          F77NAME(sdaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#endif
#else
#ifndef USE_SINGLE_PRECISION
          F77NAME(DDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#else
          F77NAME(SDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#endif
#endif
        }
      else
        {
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
          F77NAME(ddaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#else
          F77NAME(sdaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#endif
#else
#ifndef USE_SINGLE_PRECISION
          F77NAME(DDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#else
          F77NAME(SDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#endif
#endif
        }
      collectData();
      if (idid >= -1)
        info(1) = 1;
    }
  data->endUserStep();
  data->includeSolution(t,y);
  solution = y;
  solutionPrime = yp;
  if (idid >= -1)
    {
      info(1) = 1;
      return 0;
    }
  else
    return 1;
      
}


bool DASPK::step(const real& tout,real& tStep,Vec& solutionAtTStep, 
                 Vec& solutionPrimeAtTStep)
{
  info(3) = 1;              //this makes ddaspk advance only one step
  
  if (info(1) == 0) // if first call
    {
      t = dae->getT0();
      y = dae->getY0();
      yp = dae->getY0prime();
    }

  DaspkResidual::theDae = dae;
  DaspkJacobian::theJac = jacobian;
  DaspkPsol::thePrec  = preconditioner;
      if (info(12) == 0)
        {
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
          F77NAME(ddaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#else
          F77NAME(sdaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#endif
#else
#ifndef USE_SINGLE_PRECISION
          F77NAME(DDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#else
          F77NAME(SDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#endif
#endif
        }
      else
        {
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
          F77NAME(ddaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#else
          F77NAME(sdaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#endif
#else
#ifndef USE_SINGLE_PRECISION
          F77NAME(DDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#else
          F77NAME(SDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#endif
#endif
        }
  collectData();
  tStep = t;
  solutionAtTStep = y;
  solutionPrimeAtTStep = yp;
  if ( idid >= -1)
    {
      info(1) = 1;
      return 0;
    }
  else
    return 1;
}
 
  
bool DASPK::computeInitialConditions()
{
  info(14) = 1;

  t = dae->getT0();
  y = dae->getY0();
  yp = dae->getY0prime();

  DaspkResidual::theDae = dae;
  DaspkJacobian::theJac = jacobian;
  DaspkPsol::thePrec  = preconditioner;

  real tout;
      if (info(12) == 0)
        {
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
          F77NAME(ddaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#else
          F77NAME(sdaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#endif
#else
#ifndef USE_SINGLE_PRECISION
          F77NAME(DDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#else
          F77NAME(SDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_direct), DaspkPsol::psol);
#endif
#endif
        }
      else
        {
#ifndef CRAYCC
#ifndef USE_SINGLE_PRECISION
          F77NAME(ddaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#else
          F77NAME(sdaspk)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#endif
#else
#ifndef USE_SINGLE_PRECISION
          F77NAME(DDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#else
          F77NAME(SDASPK)(DaspkResidual::res,neq, t ,y.castToArray(), 
                          yp.castToArray(), tout, info.castToArray(), 
                          rtolp, atolp, idid, rwork.castToArray(),lrw, 
                          iwork.castToArray(), liw, rpar.castToArray(), ipar.castToArray(),
                          (void*)(DaspkJacobian::jac_krylov), DaspkPsol::psol);
#endif
#endif
        }
  if (idid == 4)
    {
     info(14) = 0;
     return 0;
    }
  else
    return 1;
}

  
void DASPK::setTolerances(real& atol,real& rtol)
{
  info(2) = 0;
  atolp = &atol;
  rtolp = &rtol;
}

void DASPK::setTolerances(Vec& atol,Vec& rtol)
{
  info(2) = 1;
  atolp = atol.castToArray();
  rtolp = rtol.castToArray();
}

void DASPK::setTstop(const real& tstopIn)
{
  info(4) = 1;
 
  tstop = tstopIn;       //rwork(1)
}

void DASPK::useAnalyticJacobian(Jacobian& J)
{
  info(5) = 1;
  jacobian = &J;
}

void DASPK::usePreconditioner(Preconditioner& P)
{
  info(15) = 1;
  preconditioner = &P;
}
  
void DASPK::useFullLinearSolver()
{
  info(6) = 0;
}

void DASPK::useBandedLinearSolver(int klIn,int kuIn)
{
  info(6) = 1;
  kl = klIn;         //iwork(1)
  ku = kuIn;         //iwork(2)
}

void DASPK::useKrylovSolver()
{
  info(12) = 1;
}
  
void DASPK::useKrylovSolver(int maxlIn,int kmpIn,int nrmaxIn, real& epliIn)
{
  info(12) = 1;
  info(13) = 1;
 
  maxl = maxlIn;               //iwork(24)
  kmp = kmpIn;                 //iwork(25)
  nrmax = nrmaxIn;             //iwork(26)
  epli = epliIn;               //rwork(10)
}
  
void DASPK::setMaxStepSize(const real& tmaxIn)
{
  info(7) = 1;

  tmax = tmaxIn;
}

void DASPK::setInitialStepSize(const real& h0In)
{
  info(8) = 1;
  h0 = h0In;  //rwork(3)
}

void DASPK::setMaxOrder(int maxOrder)
{
  info(9) = 1;
  maxord = maxOrder;
}

void DASPK::setInitialConstraints(const IntVec& constraintFlag)
{
  if (info(10) == 2 || info(10) == 3)
    info(10) = 3;
  else
    info(10) = 1;

  cFlag = constraintFlag;
}

void DASPK::enforceNonnegativity()
{
  if (info(10) == 1 || info(10) == 3)
    info(10) = 3;
  else
    info(10) = 2;
}
  
void DASPK::findConsistentY0()
{
  info(11) = 2;
}
  
void DASPK::findConsistentY0andY0Prime(const IntVec& variableFlag)
{
  info(11) = 1;
  avFlag = variableFlag;
}
  
void DASPK::excludeAlgebraicVariables( const IntVec& variableFlag)
{
  info(16) = 1;

  avFlag = variableFlag;
}
  
void DASPK::setIntialConditionCalcParms(int mxnitIn,int mxnjIn,int mxnhIn,
                                        int lsoffIn, const real&  stptolIn, 
                                        const real& epinitIn)
{
  info(17) = 1;

  mxnit = mxnitIn;        //iwork(32)
  mxnj  = mxnjIn;         //iwork(33)
  mxnh  = mxnhIn;         //iwork(34)
  lsoff = lsoffIn;        //iwork(35)
  stptol = stptolIn;      //rwork(14)
  epinit = epinitIn;      //rwork(15)
} 
  
void DASPK::useExtraPrinting(int printLevel)
{
  info(18) = printLevel;
}


//not finished
void DASPK::reset(){info(1) = 0;}


//not finished
void DASPK::resetTstop(){}

void DASPK::collectData()
{
  int i;
  for (i=0;i< iwork(12) - nreOld;i++)
    data->functionEvaluation();
  nreOld = iwork(12);
  
  for (i=0;i< iwork(13) - njeOld;i++)
    data->jacobianEvaluation();
  njeOld = iwork(13);
  
  for (i=0;i< iwork(14) - netfOld;i++)
    data->errorFailure();
  netfOld = iwork(14);
  
  for (i=0;i< iwork(15) - ncfnOld;i++)
    data->nonlinearSolverFailure();
  ncfnOld = iwork(15);

  for (i=0;i< iwork(16) - ncflOld;i++)
    data->linearSolverFailure();
  ncflOld = iwork(16);
  
  for (i=0;i< iwork(19) - nniOld;i++)
    data->nonlinearSolverIteration();
  nniOld = iwork(19);
  
  for (i=0;i< iwork(20) - nliOld;i++)
    data->linearSolverIteration();
  nliOld = iwork(20);
  
  data->stepTaken(iwork(8),rwork(7),rwork(4));//k,h,tn
}

void DASPK::init()
{
  //compute liw and lrw -- the necessary storage for this configuration
  //remember that vectors are base 1

  if (info(12) == 0) //standard direct method
    {
      //base values
      lrw = 50 + std::max(maxord +4,7)*neq;
      liw = 40 + neq;


      //options

      if ( info(6) == 0)    //dense matrix
        lrw+=neq*neq;
      else                  //banded matrix
        {
          if ( info(5) == 0)       //numerical jac
            lrw+= (2*kl+ku+1)*neq + 2*(neq/(kl+ku+1) + 1);
          else                     //analytical jac
            lrw+= (2*kl+ku+1)*neq;
        }
      if (info(16) == 1)
        lrw+=neq;
    }
  else //krylov method
    {
      lrw = 50 + (maxord +5)*neq + (maxl + 3 + std::max(0,std::min(1,maxl-kmp)))*neq+
        + (maxl+3)*maxl + 1 + lenwp;
      liw = 40 + leniwp;
      if (info(16) == 1)
        lrw+=neq;
    }

  //remainder of cases for liw
  if ( info(10) == 1 || info(10) == 3)
    liw+=neq;
  if ( info(11) == 1)
    liw+=neq;
  

  //allocate storage

  iwork.newsize(liw);
  iwork.setBase(1);
  rwork.newsize(lrw);
  rwork.setBase(1);
  iwork = 0;
  rwork = 0;
  //set user options (some of these will not be used depending on info)

  iwork(1) = kl;
  iwork(2) = ku;
  iwork(24) = maxl;
  iwork(25) = kmp;
  iwork(26) = nrmax;
  iwork(32) = mxnit;
  iwork(33) = mxnj;
  iwork(34) = mxnh;
  iwork(35) = lsoff;

  rwork(1) = tstop;
  rwork(2) = tmax;
  rwork(3) = h0;
  rwork(10) = epli; 
  rwork(14) = stptol;
  rwork(15) = epinit;

  
  // indicate how variables should be constrained if necessary

  if (info(10) == 1 || info(10) == 3)
    {
      IntVec iworkBlock(IntVec::REF,iwork,CMRVecIndex(40+1,40+neq));
      iworkBlock = cFlag;
    }


  // indicate algebraic variables if necessary

  if (info(11) == 1 || info(16) == 1)
    {
      if (info(10) == 0 || info(10) == 2)
        {
          IntVec iworkBlock(IntVec::REF,iwork,CMRVecIndex(40+1,40+neq));
          iworkBlock = avFlag;
        }
      if (info(10) == 1 || info(10) == 3)
        {
          IntVec iworkBlock(IntVec::REF,iwork,CMRVecIndex(40+1,40+2*neq));
          iworkBlock = avFlag;
        }
    }
}

}//Daetk
