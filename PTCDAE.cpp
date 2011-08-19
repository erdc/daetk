#include "PTCDAE.h"

namespace Daetk 
{
using std::min;
using std::max;
using std::cerr;
using std::endl;
using std::cout;
 
  void PTCDAE::linearSolverIsInexact()
{INEXACT_LINEAR_SOLVER=true;}

PTCDAE::PTCDAE():
  STEADYSTATE(12345),
  nIts(0),
  delform(1),
  maxPTCits(1000),
  nCorrectorIts(1),
  err(1),
  evalFailed(1),
  useWRMSstop(0),
  useL2stop(1),
  INEXACT_LINEAR_SOLVER(false),
  dt(12345),
  alpha(12345),
  dyNorm(12345),
  dy2Norm(12345),
  tau(12345),
  ssTol(1.0),
  ss2Tol(1.0e-10),
  dtMax(1000.0),
  dtMin(SQRT_MACHINE_EPSILON),
  t(0.0),
  ypNorm(12345),
  yp2Norm(12345),
  userDt0(0.0),
  ypNormLast(12345),
  yp2NormLast(12345),
  dtLast(12345),
  dtNew(12345),
  linearTolerance(0.001),
  iterout("iterations.grf")
{
  iterout.setf(std::ios::scientific);
  iterout.precision(10);
}
  
PTCDAE::~PTCDAE(){}
   
PTCDAE::PTCDAE(DaeDefinition& daeIn,LinearSolver& lsIn,Jacobian& jIn,
    VectorNorm& nmIn,DataCollector& dataIn, real lTol):
  STEADYSTATE(12345),
  nIts(0),
  delform(1),
  maxPTCits(1000),
  nCorrectorIts(1),
  err(1),
  evalFailed(1),
  useWRMSstop(0),
  useL2stop(1),
  INEXACT_LINEAR_SOLVER(false),
  dt(12345),
  alpha(12345),
  dyNorm(12345),
  dy2Norm(12345),
  tau(12345),
  ssTol(1.0),
  ss2Tol(1.0e-10),
  dtMax(1000.0),
  dtMin(SQRT_MACHINE_EPSILON),
  t(0.0),
  ypNorm(12345),
  yp2Norm(12345),
  userDt0(0.0),
  ypNormLast(12345),
  yp2NormLast(12345),
  dtLast(12345),
  dtNew(12345),
  linearTolerance(lTol),
  residual(daeIn.getY0().dim(),13245),
  residualLast(daeIn.getY0().dim(),13245),
  solutionLast(daeIn.getY0().dim(),13245),
  solutionBeforeLast(daeIn.getY0().dim(),13245),
  solutionPrimeLast(daeIn.getY0().dim(),13245),
  correction(daeIn.getY0().dim(),13245),
  correctionLast(daeIn.getY0().dim(),13245),
  totalCorrection(daeIn.getY0().dim(),13245),
  dyDiff(daeIn.getY0().dim(),13245),
  xtt(daeIn.getY0().dim(),13245),
  xtt_a(daeIn.getY0().dim(),13245),
  iterout("iterations.grf"),
  dae(&daeIn),
  data(&dataIn),
  weightedNorm(&nmIn),
  jacobian(&jIn),
  linearSolver(&lsIn)
{
  std::cerr<<"int here"<<std::endl;
  iterout.setf(std::ios::scientific);
  iterout.precision(10);
}
  
void PTCDAE::readParameters(ParameterDatabase& pd)
{
  delform = pd.i("delform");
  tau = pd.r("tau");
  ssTol=pd.r("ssTol");
  ss2Tol=pd.r("ss2Tol");
  useWRMSstop=pd.b("useWRMSstop");
  useL2stop = !useWRMSstop;
  dt   = pd.r("dtPtcOrig");
  userDt0=dt;
  dtMax= pd.r("dtMax");
  dtMin= pd.r("dtMin");
  maxPTCits = pd.i("maxPTCits");
  nCorrectorIts = pd.i("nCorrectorIts");
}

bool PTCDAE::step(const real& tout,real& tStep,Vec& solution, Vec& solutionPrime)
{

//   cerr<<"step t="<<t<<endl;

//   iterout<<dt<<endl;
  
//   cerr<<"it = "<<nIts<<endl;
  
  
  //for history
  dtLast=dt;
  solutionBeforeLast=solutionLast;
  solutionLast = solution;
  solutionPrimeLast = solutionPrime;
  residualLast = residual;
  correctionLast=correction;
  err=true;

  //try to take a step, watch for jac eval error, linear solver err, and fEval errors
  while (err)
    {  
      correction=0.0;
      totalCorrection=0.0;
      nIts++;
      tStep = t + dt;
      int correctorIts=0;
      err = evaluateJacobian();
      while (correctorIts < nCorrectorIts && !err)
//       while (correctorIts < nCorrectorIts)
        {
          correctorIts++;
          if (correctorIts > 1)
            {
              data->jacobianEvaluation();
              //   computeDeltaForJacobian();
              err=jacobian->evaluate(solution,residual);
              dae->initializeFunction(solution,solutionPrime,residual);
//               err = evaluateJacobian();
            }
          if (err)
            {
              cerr<<"Jacobian evaluation error"<<endl;
            }
          else
            {
              err = linearSolver->prepare();
              if (err)
                {
                  cerr<<"linear solver prepare failure"<<endl;
                }
              else
                {
                  err = linearSolver->solve(residual,correction);
                  if (err)
                    {
                      cerr<<"linear solver failure"<<endl;
                    }
                  else
                    {
                      dae->correctArgument(correction);
                      totalCorrection+=correction;
                      residual = dae->value(err);
                      //residual = Fvalue(err);
                      if (err)
                        {
                          cerr<<"function eval error"<<endl;
                          cerr<<"correctorIts "<<correctorIts<<endl;
                        }
                      else
                        {
                          solution = dae->yDaeDef;
                          solutionPrime = dae->ypDaeDef;
                        }//fEval
                    }//lSolve
                }//lPrepare
            }//jEval
        }//corrector Its
      if (err)
        {
          cerr<<"halving the step size because of function eval erro"<<endl;
          dt*=0.5;
          if (dt < dtMin )
            {
              cout <<" step size dt = "<<dt<<" below minimum allowed = "
                   <<dtMin<<" quitting "<<endl;
              cerr <<" step size dt = "<<dt<<" below minimum allowed = "
                   <<dtMin<<" quitting "<<endl;
              return true;
              //methodFailure = " \\delta \\longrightarrow 0 ";
            }
        }//err
      if (nIts > maxPTCits)
        {
          cerr<<"quitting, exceeded max steps in ptc ="<<nIts<<endl;
          return true;
        }
    }//while
  
  residual = Fvalue(err);

  dy2NormLast = dy2Norm;
  dy2Norm = nrm2(totalCorrection);
// //   cerr <<"dy2Norm = "<<dyNorm<<endl;

// //   dyNormLast = dyNorm;
// //   dyNorm = (*weightedNorm)(correction);
// //   cerr <<"dyNorm = "<<dyNorm<<endl;
  
  yp2NormLast = yp2Norm;
  yp2Norm = nrm2(solutionPrime);
//   cerr <<"yp2Norm = "<<yp2Norm<<endl;

//   ypNormLast = ypNorm;
//   ypNorm = (*weightedNorm)(solutionPrime);
//   cerr <<"ypNorm = "<<ypNorm<<endl;

  res2NormLast = res2Norm;
  res2Norm = nrm2(residual);
//   cerr<<"res2Norm = "<<res2Norm<<endl;

//   weightedNorm->setWeight(solution);
  data->stepTaken(1,dt,t,res2Norm);

  return false;
}

void PTCDAE::chooseDt(const Vec& solution, const Vec& solutionPrime)
{
  //choose new dt
  if (delform==1) // SER
    {
      //use norm of F 
//       dt = dt*yp2NormLast/yp2Norm;
      dt = dt*res2NormLast/res2Norm;
//       cerr <<"SER delta update "<<res2NormLast<<'\t'<<res2Norm<<'\t'<<dt<<endl;
    }
  else if (delform==2) // SSF (Step Size Formula)
    {
      //try norm of step diff
      if (dy2Norm < 1.0e-12)
        dt = dtMax;
      else
        dt = dt/dy2Norm;
//       cerr <<"SSF delta update"<<endl;
    }
  else if (delform==3) // TTE
    {
      //mwf TTE isn't picking up the solution blowing up.
      //mwf it's using dt = 1000
      for (int i=0;i<solution.getLocalHigh();i++)
        {
          xtt(i) = 
            (2.0/(dt+dtLast))*( (solution[i]-solutionLast[i])/dt 
                               -(solutionLast[i]-solutionBeforeLast[i])/dtLast);
        }
      for (int i=0;i<solution.getLocalHigh();i++)
        {
          xtt_a[i] = fabs(xtt[i]) / (2.0*(1.0+fabs(solution[i])));
        }
      
      //mwf Feb 07 try 
//       cerr <<"xtt_a.max() = "<<xtt_a.max()<<endl;
      real l2_infty = max(xtt_a)/xtt_a.dim();
      if ( l2_infty < tau/dtMax )
        { dtNew = dtMax; }
      else
        { dtNew = sqrt(tau/l2_infty); }
      
      dt = dtNew;
      
//       cerr <<"TTE delta update"<<endl;
      
    } //end TTE
  else if (delform==4) //weighted norm SER
    {
      //use norm of F 
      dt = dt*ypNormLast/ypNorm;
      cerr <<"SERW delta update"<<endl;
      //                    dt = dt*(ypNormLast/ypNorm);
      cerr <<"SERW delta update "<<yp2NormLast<<'\t'<<yp2Norm<<'\t'<<dt<<endl;
    }
  else if (delform==5) //weighted norm SSF
    {
      //try norm of step diff
      if (dyNorm < 1)
        dt = dtMax;
      else
        dt = dt/dyNorm;
      cerr <<"SSFW delta update"<<endl;
    }
  else // Unknown
    {
      cerr <<"Unknown delta update formula."<<endl;
    }
  if (dt > dtMax)
    dt = dtMax;
  
//   cerr <<"dt = "<<dt <<endl;
}

bool PTCDAE::calculateSolution(const real& tout,Vec& solution, Vec& solutionPrime)
{
  cerr<<"trying to start PseudoTC"<<endl;
  data->startUserStep();
  correction = 0.0;

  //get initial values and setup DaeDefinition for initial step
  t=dae->getT0();
  solution = dae->getY0();
  solutionPrime = dae->getY0prime();
  dae->yDaeDef=solution;
  residual = Fvalue(evalFailed);
  if (evalFailed)
    {
      cerr <<"evalFailed for initial values, exiting"<<endl;
      exit(1);
    }

  //weightedNorm->setWeight(solution);
  
  
  nIts=0;
  
  
  //norms
  //ypNorm = (*weightedNorm)(solutionPrime);
  yp2Norm = nrm2(solutionPrime);
  res2Norm = nrm2(residual);
  res0 = res2Norm;
  ypNormLast = ypNorm;
  yp2NormLast = yp2Norm;
  res2NormLast = res2Norm;
  dyNorm = 0.0;

  //cerr<<"original yp2Norm  = "<<yp2Norm<<endl;
  //cerr<<"original res2Norm = "<<res2Norm<<endl;

  //pick initial dt
  if (userDt0)
    dt = userDt0;
  else
    {
      if (ypNorm == 0.0)
        {
          cerr<<"Norm of initial yPrime is 0.0, exiting"<<endl;
          exit(1);
        }
      dt=0.5/ypNorm;
    }

  //pretend we did a forward Euler step to get here
  dtLast=dt;
//   solutionLast = solution - dt*solutionPrime;
  solutionLast = solutionPrime;
  scal(-dt,solutionLast);
  axpy(1.0,solution,solutionLast);
//   correction = dt*solutionPrime;
  correction = solutionPrime;
  scal(dt,correction);

  bool stepFailed;
  real tStep;
  while (!converged(solution,solutionPrime))
    {
     stepFailed=step(STEADYSTATE,tStep,solution,solutionPrime);
     if (stepFailed) 
       {
         data->includeSolution(-1,solution);
         return true;
       }
     t=tStep;
     chooseDt(solution,solutionPrime);
    }
  data->endUserStep();
  data->includeSolution(t,solution);
  return false;
}

bool PTCDAE::converged(const Vec& solution, const Vec& solutionPrime)
{
//   if (INEXACT_LINEAR_SOLVER)
//     {
//       if (dyNorm <= linearTolerance)
//         {
          
//           cerr<<"dyNorm is less than the linear iterative solver tolerance, decreasing linearTolerance"<<endl;
//           cerr<<"Old linearTolerance "<<linearTolerance;
//           linearTolerance*=.1;
//           cerr<<" New linearTolerance "<<linearTolerance<<endl;
//         }
//     }
//   else
//     {
//       roundOffTolerance = 100.0 * (*weightedNorm)(solution) * MACHINE_EPSILON;
//       if (dyNorm <= roundOffTolerance)
//         {
//           cerr<<"dyNorm is less than roundOffTolerance, continuing anyway"<<endl;
//         }
//     }
 
//   if (dyNorm <= roundOffTolerance)
//     {
//       cerr<<"dyNorm less than round off or linear solver accuracy, continuing iteration anyway"<<endl;
//       cerr<<"dyNorm = "<<dyNorm<<" round off/linear tolerance = "<<roundOffTolerance<<endl;
//     }
  if (ypNorm <= ssTol && useWRMSstop)
    {
      cerr<<"ypNorm less than tol = "<<ypNorm<<endl;
      cerr<<"while yp2Norm = "<<yp2Norm<<endl;
      return true;
    }
  else if (res2Norm <= res0*ss2Tol+ss2Tol && useL2stop)
    {
      cerr<<"yp2Norm less than tol = "<<yp2Norm<<endl;
      cerr<<"while ypNorm = "<<ypNorm<<endl;
      return true;
    }
  else  
    return false;
}

bool PTCDAE::evaluateJacobian()
{
  //set up DaeDefinition for Backward Euler with initial guess given by the previous solution 
  //i.e. yp=0.0 initially
  dae->alphaDaeDef = 1/dt;
//cek test   dae->alphaDaeDef = -1/dt;
//   dae->betaDaeDef = -dae->alphaDaeDef*solutionLast;
  dae->betaDaeDef = solutionLast;
  scal(-dae->alphaDaeDef,dae->betaDaeDef);
  dae->tDaeDef = t;
  dae->yDaeDef = solutionLast;
  dae->ypDaeDef = 0.0;
  dae->resetFunction();

  data->jacobianEvaluation();
//   computeDeltaForJacobian();
  bool evalError=false;
  evalError=jacobian->evaluate(solutionLast,residualLast);
  dae->ypDaeDef = 0.0;
  dae->initializeFunction(solutionLast,dae->ypDaeDef,residualLast);
  if (evalError)
    return evalError;
  return false;
}

void PTCDAE::computeDeltaForJacobian()
{
  for (int i=0; i<solutionLast.ldim(); i++)
    {
      dae->deltaVF[i]=SQRT_MACHINE_EPSILON*fabs(solutionLast[i])+SQRT_MACHINE_EPSILON;
      dae->deltaVF[i]=copysign(dae->deltaVF[i],correctionLast[i]); //this is just + on the first call
      dae->deltaVF[i]=solutionLast[i] - (solutionLast[i] - dae->deltaVF[i]);
    }
}

const Vec& PTCDAE::Fvalue(bool& errFlag)
{
  dae->alphaDaeDef = 0.0;
  dae->betaDaeDef = 0.0;
  dae->tDaeDef = t;
  dae->ypDaeDef = 0.0; 
  dae->resetFunction();
  return dae->value(errFlag);
}

void PTCDAE::reset()
{
  nIts=0;
}

}//Daetk
