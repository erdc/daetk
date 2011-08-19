#include "BiCGstab.h"
namespace Daetk 
{
using std::cerr;
using std::endl;
using std::ostream;
 using std::operator|;
 using std::operator&;
void BiCGstab::useInitialGuess(){USE_GUESS=true;}

void BiCGstab::checkTrueResidual(real resTol)
{CHECK_RESIDUAL = true; trueResTol = resTol;}

BiCGstab::BiCGstab():  
  RESTART(false),
  CHECK_RESIDUAL(false),
  error(true),
  USE_GUESS(false),
  maxIterations(0),
  tol(0.0),
  x(this),
  b(this),
  vectorNorm(0),
  data(0),
  linearOperator(0),
  preconditioner(0)
{
  Tracer tr("BiCGstab::BiCGstab()");
}


BiCGstab::BiCGstab(LinearOperator& linOpIn, 
                   Preconditioner& precIn, 
                   VectorNorm& normIn, 
                   DataCollector& dataIn, 
                   ParameterDatabase& pd):
  RESTART(pd.b("RESTART")),
  CHECK_RESIDUAL(pd.b("CHECK_RESIDUAL")),
  error(pd.b("error")),
  USE_GUESS(false),
  maxIterations(pd.i("maxIterations")),
  tol(pd.r("tol")),
  trueResTol(pd.r("trueResTol")),
  rhoK(pd.r("rhoK")),
  rhoKm1(pd.r("rhoKm1")),
  alpha(pd.r("alpha")),
  omega(pd.r("omega")),
  k(pd.r("k")),
  beta(pd.r("beta")),
  realTemp(pd.r("realTemp")),
  rNorm(pd.r("rNorm")),
  xNorm(pd.r("xNorm")),
  r0(pd.r("r0")),
  r(pd.v("r")),
  rHat(pd.v("rHat")),
  v(pd.v("v")),
  p(pd.v("p")),
  s(pd.v("s")),
  unscaled_p(pd.v("unscaled_p")),
  unscaled_s(pd.v("unscaled_s")),
  t(pd.v("t")),
  vecTemp(pd.v("vecTemp")),
  x(this,vecTemp.dim()),
  b(this,vecTemp.dim()),
  vectorNorm(&normIn),
  data(&dataIn),
  linearOperator(&linOpIn),
  preconditioner(&precIn)
{  
  Tracer tr("BiCGstab::BiCGstab(BandColMat& BCMin)");
}

BiCGstab::BiCGstab(LinearOperator& linOpIn, 
                   Preconditioner& precIn, 
                   VectorNorm& normIn, 
                   DataCollector& dataIn, 
                   real tolIn,
                   int maxIts,int neq):
  RESTART(false),
  CHECK_RESIDUAL(false),
  error(true),
  USE_GUESS(false),
  maxIterations(maxIts),
  tol(tolIn),
  r(neq,12345),
  rHat(neq,12345),
  v(neq,12345),
  p(neq,12345),
  s(neq,12345),
  unscaled_p(neq,12345),
  unscaled_s(neq,12345),
  t(neq,12345),
  vecTemp(neq,12345),
  x(this,neq),
  b(this,neq),
  vectorNorm(&normIn),
  data(&dataIn),
  linearOperator(&linOpIn),
  preconditioner(&precIn)
{  
  Tracer tr("BiCGstab::BiCGstab(BandColMat& BCMin)");
}

BiCGstab::~BiCGstab()
{
  Tracer tr("BiCGstab::~BiCGstab()");
}

bool BiCGstab::prepare()
{
  return preconditioner->prepare();
}

bool BiCGstab::solve(const Vec& bIn, Vec& xIn)
{
  b.attachToTarget(bIn);
  x.attachToTarget(xIn);

  //1. Initialize.  The real system is A_0 x = b_0
  // we solve DMA_0D^-1 z = DMb_0 <=> A_0 x = b_0; A = DMA_0, b = DMb_0, x = D^-1z

  if (!USE_GUESS)
    x.v_=0.0;
  bool evalError=false;
  real tau;
//    VecIndex ind(0,r.dim()-1);
//    bRef.attachToVecMulti(Vec::REF,b,ind);
//    xRef.attachToVecMulti(Vec::REF,x,ind);
//    if (r.dim() < x.dim())
//    {
//      VecIndex indRem(r.dim(),x.dim()-1);
//      bRemRef.attachToVecMulti(Vec::REF,b,indRem);
//      xRemRef.attachToVecMulti(Vec::REF,x,indRem);
    
//  #ifndef USE_BLAS
//      xRemRef = bRemRef;
//  #else
//      copy(bRemRef,xRemRef);
//  #endif
//    }
  
  if (!RESTART)  //first iteration
    {
      if (USE_GUESS)
        {
          evalError = linearOperator->apply(x.v_,vecTemp);               //A_0 x
          while (evalError)
            {
              cerr<<"trying to cut back on x "<<endl;
              x.v_*=0.5;
              evalError = linearOperator->apply(x.v_,vecTemp);               //A_0 x
            }
#ifndef USE_BLAS
          rHat = b.v_ - vecTemp;                             //b - A_0 x  !!rHat is just used as a temp
#else
          copy(b.v_,rHat);
          axpy(-1.0,vecTemp,rHat);
#endif
          if (CHECK_RESIDUAL)
            r0 = nrm2(rHat);
          error = preconditioner->apply(rHat,r);
        }
      else
        {
          if (CHECK_RESIDUAL)
            r0 = nrm2(b.v_);
          error = preconditioner->apply(b.v_,r);
        }
      if (error) return true;
      
  
      vectorNorm->scale(r,rHat);                      //r = DM(b-A_0 x)  = DM(b-A_0D^-1z)
      
#ifndef USE_BLAS
      r = rHat;                                   
#else
      copy(rHat,r);
#endif
//        cout<<"res"<<endl<<r;
//    xNorm = (*vectorNorm)(x);
      //check r 
      rNorm = nrm2(r);
//        if (rNorm <= MACHINE_EPSILON*xNorm)
//          {
//            cerr<<"Early exit for BiCGstab with preconditioned residual = "<<rNorm<<endl;
//            return false;
//          }
      
      rhoK = alpha = omega = 1.0;
      
      p = v = s = t = 0.0;
      
      k = 0;
    }
  
  while((rNorm >= tol && k < maxIterations) || k==0)
    {
      //cout<<rNorm<<endl;
      data->linearSolverIteration();
      ++k;

      rhoKm1 = rhoK;
      rhoK = dot(rHat,r);
      beta = (rhoK/rhoKm1)*(alpha/omega);

#ifndef USE_BLAS
      p = r + beta*(p - omega*v);
#else
      axpy(-omega,v,p);
      scal(beta,p);
      axpy(1.0,r,p);
#endif
      vectorNorm->deScale(p,unscaled_p);                        //D^-1p
      evalError=linearOperator->apply(unscaled_p,v);            //A_0D^-1p
      while(evalError)
        {
          cerr<<"cutting back p in BiCGstab"<<endl;
          p*=0.5;
          vectorNorm->deScale(p,unscaled_p);                    //D^-1p
          evalError=linearOperator->apply(unscaled_p,v);        //A_0D^-1p
        }
      error = preconditioner->apply(v,vecTemp);                 //MA_0D^-1p
      if (error) return true;
      vectorNorm->scale(vecTemp,v);                             //v = DMA_0D^-1p
      //compute alpha
      tau = dot(rHat,v);
//        I commented out the old test an replaced with the simple one--cek
//        tau = dot(rHat,v);
//        if (fabs(tau) < MACHINE_EPSILON*xNorm && 
//  	  fabs(tau) < 100*MACHINE_EPSILON*rhoK ) 
      if (fabs(tau) == 0.0) //this test may need improvement
        {
          //unhappy breakdown
          cerr<<"BiCGstab Failure: v and rHat are orthogonal"<<endl;
//  	  cerr <<"dot(rHat,v)= "<<realTemp<<"  xNorm = "<<xNorm<<endl;
//  	  cerr <<"nrm2(v) = "<<nrm2(v)<<"  rhoK= "<<rhoK<<endl;
	  cerr <<"dot(rHat,v)= "<<tau<<endl;
	  cerr <<"nrm2(v) = "<<nrm2(v)<<"  rhoK= "<<rhoK<<endl;
          return true;
        }
      else
        alpha = rhoK/tau;

#ifndef USE_BLAS
      s = r - alpha*v;
#else
      copy(r,s);
      axpy(-alpha,v,s);
#endif
      vectorNorm->deScale(s,unscaled_s);
      evalError=linearOperator->apply(unscaled_s,t);
      while(evalError)
        {
          cerr<<"cutting back s in BiCGstab"<<endl;
          s*=0.5;
          vectorNorm->deScale(s,unscaled_s);
          evalError=linearOperator->apply(unscaled_s,t);
        }
      error = preconditioner->apply(t,vecTemp);
      if (error) return true;
      vectorNorm->scale(vecTemp,t);
      
      tau = dot(t,t);
      if (fabs(tau) == 0.0) //this test may need improvement
        {
          //happy breakdown
          cerr<<"Early exit from BiCGstab"<<endl;
#ifndef USE_BLAS
          x.v_+=alpha*unscaled_p;
#else
          axpy(alpha,unscaled_p,x.v_);
#endif
          
          if (SOLVE_SUB)
            xIn=bIn;
          x.restoreToTarget();
          return false; 
        }
      omega = dot(t,s)/tau;
      if (fabs(omega) == 0)
        {
          cerr<<"omega==0 in BiCGS: Failure exiting"<<endl;
          return true;
        }
     
#ifndef USE_BLAS
      x.v_+= alpha*unscaled_p + omega*unscaled_s;
      r = s - omega*t;
#else
      axpy(alpha,unscaled_p,x.v_);
      axpy(omega,unscaled_s,x.v_);
      copy(s,r);
      axpy(-omega,t,r);
#endif
//        xNorm = (*vectorNorm)(x);
      rNorm = nrm2(r);
      //cerr<<"nrm2(r)"<<rNorm<<endl;
      if (CHECK_RESIDUAL)
        {
          linearOperator->apply(x.v_,vecTemp);
#ifndef USE_BLAS
          vecTemp-=b.v_;
#else
          axpy(-1.0,b.v_,vecTemp);
#endif
          //we'll assume that the user has incorporated ||b|| into the tolerances
          //cerr<<nrm2(vecTemp)<<endl;
          if (nrm2(vecTemp) <= trueResTol*r0 + trueResTol)
            {
              if (SOLVE_SUB)
                xIn=bIn;
              x.restoreToTarget();
              return false;
            }
          else
            rNorm = tol+1; //force another iteration
        }
    }
  if(k == maxIterations)
    {
      cerr<<"failure in BiCGstab"<<endl;
      return true;
    }
  else
    {
      //cerr<<"rNorm"<<r<<endl;
      if (SOLVE_SUB)
        xIn=bIn;
      x.restoreToTarget();
      return false;
    }
}


void BiCGstab::dump(ostream& os)
{
  os<<"bool"<<'\t'<<"RESTART"<<'\t'<<RESTART<<endl
    <<"bool"<<'\t'<<"CHECK_RESIDUAL"<<'\t'<<CHECK_RESIDUAL<<endl
    <<"bool"<<'\t'<<"error"<<'\t'<<error<<endl
    <<"int"<<'\t'<<"maxIterations"<<'\t'<<maxIterations<<endl
    <<"real"<<'\t'<<"tol"<<'\t'<<tol<<endl
    <<"real"<<'\t'<<"trueResTol"<<'\t'<<trueResTol<<endl
    <<"real"<<'\t'<<"rhoK"<<'\t'<<rhoK<<endl
    <<"real"<<'\t'<<"rhoKm1"<<'\t'<<rhoKm1<<endl
    <<"real"<<'\t'<<"alpha"<<'\t'<<alpha<<endl
    <<"real"<<'\t'<<"omega"<<'\t'<<omega<<endl
    <<"real"<<'\t'<<"k"<<'\t'<<k<<endl
    <<"real"<<'\t'<<"beta"<<'\t'<<beta<<endl
    <<"real"<<'\t'<<"tau"<<'\t'<<realTemp<<endl
    <<"real"<<'\t'<<"rNorm"<<'\t'<<rNorm<<endl
    <<"real"<<'\t'<<"xNorm"<<'\t'<<xNorm<<endl
    <<"real"<<'\t'<<"r0"<<'\t'<<r0<<endl
    <<"Vec"<<'\t'<<"r"<<'\t'<<r.dim()<<endl<<r<<endl
    <<"Vec"<<'\t'<<"rHat"<<'\t'<<rHat.dim()<<endl<<rHat<<endl
    <<"Vec"<<'\t'<<"v"<<'\t'<<v.dim()<<endl<<v<<endl
    <<"Vec"<<'\t'<<"p"<<'\t'<<p.dim()<<endl<<p<<endl
    <<"Vec"<<'\t'<<"s"<<'\t'<<s.dim()<<endl<<s<<endl
    <<"Vec"<<'\t'<<"unscaled_p"<<'\t'<<unscaled_p.dim()<<endl<<unscaled_p<<endl
    <<"Vec"<<'\t'<<"unscaled_s"<<'\t'<<unscaled_s.dim()<<endl<<unscaled_s<<endl
    <<"Vec"<<'\t'<<"t"<<'\t'<<t.dim()<<endl<<t<<endl
    <<"Vec"<<'\t'<<"vecTemp"<<'\t'<<vecTemp.dim()<<endl<<vecTemp<<endl;
}

}//Daetk
