#include "PTC.h"

namespace Daetk 
{
using std::min;
using std::max;
using std::cerr;
using std::endl;
using std::cout;
 
  void PTC::linearSolverIsInexact()
{INEXACT_LINEAR_SOLVER=true;}

PTC::PTC():
  STEADYSTATE(12345),
  nIts(0),
  delform(1),
  maxPTCits(1000),
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
  ypNormOld(12345),
  yp2NormOld(12345),
  dtOld(12345),
  dtNew(12345),  
  linearTolerance(1.0e-4),
  residual(dae->getY0().dim(),13245),
  solutionOld(dae->getY0().dim(),13245),
  solutionPrev(dae->getY0().dim(),13245),
  correction(dae->getY0().dim(),13245),
  solutionPrimeSave(dae->getY0().dim(),13245),
  residualSave(dae->getY0().dim(),13245),
  dyDiff(dae->getY0().dim(),13245),
  xtt(dae->getY0().dim(),13245),
  xtt_a(dae->getY0().dim(),13245),
  iterout("iterations.grf")
{
  iterout.setf(std::ios::scientific);
  iterout.precision(10);
}
  
PTC::~PTC(){}
   
PTC::PTC(DaeDefinition& daeIn,LinearSolver& lsIn,Jacobian& jIn,
    VectorNorm& nmIn,DataCollector& dataIn, real lTol):
  STEADYSTATE(12345),
  nIts(0),
  delform(1),
  maxPTCits(1000),
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
  ypNormOld(12345),
  yp2NormOld(12345),
  dtOld(12345),
  dtNew(12345),
  linearTolerance(lTol),
  residual(daeIn.getY0().dim(),13245),
  solutionOld(daeIn.getY0().dim(),13245),
  solutionPrev(daeIn.getY0().dim(),13245),
  correction(daeIn.getY0().dim(),13245),
  solutionPrimeSave(daeIn.getY0().dim(),13245),
  residualSave(daeIn.getY0().dim(),13245),
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
  iterout.setf(std::ios::scientific);
  iterout.precision(10);
}
  
void PTC::readParameters(ParameterDatabase& pd)
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
}

bool PTC::step(const real& tout,real& tStep,Vec& solution, Vec& solutionPrime)
{

  cerr<<"step t="<<t<<endl;

  iterout<<dt<<endl;
  
  cerr<<"it = "<<nIts<<endl;
  
  
  //for history
  solutionOld = solutionPrev;
  solutionPrev = solution;
  
  residual*=-1.0;

  err=true;

  //try to take a step, watch for jac eval error, linear solver err, and fEval errors
  while (err)
    {  
      nIts++;
      tStep = t + dt;
      alpha = 1.0/dt;
      //set up DaeDefinition
      dae->alphaDaeDef = alpha;
      dae->betaDaeDef = 0.0;
      dae->tDaeDef = t;
      dae->yDaeDef = solution;
      dae->ypDaeDef = solutionPrime;
      //might need to fix numerical jacobian here. The #2 arg is usually res (i.e. the func val)
      err = evaluateJacobian(solution,residual);
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
                  bool err2;
                  //cerr<<"in PTC solution= "<<solution<<endl;
                  //cerr<<"in PTC residual= "<<residual<<endl;
                  //cerr<<"in PTC correction= "<<correction<<endl;
                  solution +=correction;
                  err = dae->yPrimeValue(tStep,solution,solutionPrimeSave);
                  err2 = dae->residual(tStep,solution,solutionPrime,residualSave);
                  err = err || err2;
                  //cerr<<"in PTC solutionPrimeSave= "<<solutionPrimeSave<<endl;
                  if (err)
                    {
                      cerr<<"function eval error"<<endl;
                      solution=solutionPrev;
                    }
                  else
                    {
                      solutionPrime = solutionPrimeSave;
                      residual=residualSave;
                    }//fEval
                }//lSolve
            }//lPrepare
        }//jEval
      if (err)
        {
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
          if (nIts > maxPTCits)
            {
              cerr<<"quitting, exceeded max steps in ptc ="<<nIts<<endl;
              return true;
            }
        }//err                   
    }//while

  data->stepTaken(1,dt,t);

  dy2NormOld = dy2Norm;
  dy2Norm = nrm2(correction);
  cerr <<"dy2Norm = "<<dyNorm<<endl;

  dyNormOld = dyNorm;
  dyNorm = (*weightedNorm)(correction);
  cerr <<"dyNorm = "<<dyNorm<<endl;
  
  yp2NormOld = yp2Norm;
  yp2Norm = nrm2(solutionPrime);
  cerr <<"yp2Norm = "<<yp2Norm<<endl;

  ypNormOld = ypNorm;
  ypNorm = (*weightedNorm)(solutionPrime);
  cerr <<"ypNorm = "<<ypNorm<<endl;

  res2NormOld = res2Norm;
  res2Norm = nrm2(residual);
  cerr<<"res2Norm = "<<res2Norm<<endl;

  weightedNorm->setWeight(solution);

  return false;
}

void PTC::chooseDt(const Vec& solution, const Vec& solutionPrime)
{
  //choose new dt
  if (delform==1) // SER
    {
      //use norm of F 
      dt = dt*yp2NormOld/yp2Norm;
      cerr <<"SER delta update"<<endl;
      //                    dt = dt*(ypNormOld/ypNorm);
      cerr <<"SER delta update "<<yp2NormOld<<'\t'<<yp2Norm<<'\t'<<dt<<endl;
    }
  else if (delform==2) // SSF (Step Size Formula)
    {
      //try norm of step diff
      if (dy2Norm < 1.0e-12)
        dt = dtMax;
      else
        dt = dt/dy2Norm;
      cerr <<"SSF delta update"<<endl;
    }
  else if (delform==3) // TTE
    {
      //mwf TTE isn't picking up the solution blowing up.
      //mwf it's using dt = 1000
      if (nIts == 1) 
        {
          dtOld = dt;
          solutionOld = solution;
        }
      //mwf was solution(i)-solutionOld(i)/dtOld
      for (int i=0;i<solution.getLocalHigh();i++)
        {
          xtt(i) = 
            (2.0/(dt+dtOld))*( (solution[i]-solutionPrev[i])/dt 
                               -(solutionPrev[i]-solutionOld[i])/dtOld);
        }
      for (int i=0;i<solution.getLocalHigh();i++)
        {
          xtt_a[i] = fabs(xtt[i]) / (2.0*(1.0+fabs(solution[i])));
        }
      
      //mwf Feb 07 try 
      cerr <<"xtt_a.max() = "<<xtt_a.max()<<endl;
      if ( max(xtt_a) < tau/dtMax )
        { dtNew = dtMax; }
      else
        { dtNew = sqrt(tau/max(xtt_a)); }
      
      dtOld = dt;
      dt = dtNew;
      
      cerr <<"TTE delta update"<<endl;
      
    } //end TTE
  else if (delform==4) //weighted norm SER
    {
      //use norm of F 
      dt = dt*ypNormOld/ypNorm;
      cerr <<"SERW delta update"<<endl;
      //                    dt = dt*(ypNormOld/ypNorm);
      cerr <<"SERW delta update "<<yp2NormOld<<'\t'<<yp2Norm<<'\t'<<dt<<endl;
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
  
  cerr <<"dt = "<<dt <<endl;
}

bool PTC::calculateSolution(const real& tout,Vec& solution, Vec& solutionPrime)
{
  cerr<<"trying to start PseudoTC"<<endl;
  
  //methodName = " PTC ";
  
  correction = 0.0;
  
  t=dae->getT0();
  solution = dae->getY0();

  cerr <<"dt = "<<dt <<endl;
  
  
  weightedNorm->setWeight(solution);
  
  //for history
  solutionPrev     = solution;
  
  nIts=0;
  
  evalFailed = dae->yPrimeValue(dae->getT0(),solution,solutionPrime);
  bool temp  = dae->residual(t,solution,solutionPrime,residual);

  evalFailed = evalFailed || temp;
  if (evalFailed)
    {
      cerr <<"evalFailed for initial values, exiting"<<endl;
      exit(1);
    }
  
  ypNorm = (*weightedNorm)(solutionPrime);
  yp2Norm = nrm2(solutionPrime);
  res2Norm = nrm2(residual);
  res0 = res2Norm;
  ypNormOld = ypNorm;
  yp2NormOld = yp2Norm;
  res2NormOld = res2Norm;
  dyNorm = 0.0;

  cerr<<"original yp2Norm  = "<<yp2Norm<<endl;
  cerr<<"original res2Norm = "<<res2Norm<<endl;
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
  
  dtOld = 0.0;
  dtNew = 0.0;
  
  residualSave = residual;
  solutionPrimeSave = solutionPrime;

   bool stepFailed;
  real tStep;
  while (!converged(solution,solutionPrime))
    {
     stepFailed=step(STEADYSTATE,tStep,solution,solutionPrime);
      if (stepFailed) return true;
      t=tStep;
      chooseDt(solution,solutionPrime);
    }
  data->includeSolution(t,solution);
  return false;
}

bool PTC::converged(const Vec& solution, const Vec& solutionPrime)
{
  if (INEXACT_LINEAR_SOLVER)
    {
      if (dyNorm <= linearTolerance)
        {
          
          cerr<<"dyNorm is less than the linear iterative solver tolerance, decreasing linearTolerance"<<endl;
          cerr<<"Old linearTolerance "<<linearTolerance;
          linearTolerance*=.1;
          cerr<<" New linearTolerance "<<linearTolerance<<endl;
        }
    }
  else
    {
      roundOffTolerance = 100.0 * (*weightedNorm)(solution) * MACHINE_EPSILON;
      if (dyNorm <= roundOffTolerance)
        {
          cerr<<"dyNorm is less than roundOffTolerance, continuing anyway"<<endl;
        }
    }
 
  if (dyNorm <= roundOffTolerance)
    {
      cerr<<"dyNorm less than round off or linear solver accuracy, continuing iteration anyway"<<endl;
      cerr<<"dyNorm = "<<dyNorm<<" round off/linear tolerance = "<<roundOffTolerance<<endl;
    }
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

bool PTC::evaluateJacobian(const Vec& solution, const Vec& solutionPrime)
{
  data->jacobianEvaluation();
  computeDeltaForJacobian(solution,solutionPrime);
  bool evalError=false;
  evalError=jacobian->evaluate(solution,residual);
  if (evalError)
    return evalError;
  dae->initializeFunction(solution,solutionPrime,residual);
  return false;
}

void PTC::computeDeltaForJacobian(const Vec& solution, const Vec& solutionPrime)
{
  for (int i=0; i<solution.ldim(); i++)
    {
      dae->deltaVF[i]=SQRT_MACHINE_EPSILON*
	max3(fabs(solutionPrime[i]),fabs(dt*(solutionPrime[i])),1.0/(weightedNorm->getWeight()[i]));
      dae->deltaVF[i]=-copysign(dae->deltaVF[i],dt*(solutionPrime[i]));
      //make sure roundoff doesn't change the y
      dae->deltaVF[i]=solutionPrime[i] - (solutionPrime[i] - dae->deltaVF[i]);
      //the minus here is for the correctArgument routine
    }
}

void PTC::reset()
{
  nIts=0;
}

}//Daetk
