#include "PreCG.h"
namespace Daetk 
{

using std::cerr;
using std::endl;
  
void PreCG::useInitialGuess(){USE_GUESS=true;}

PreCG::PreCG(LinearOperator& linOpIn, 
                   Preconditioner& precIn, 
                   VectorNorm& normIn, 
                   DataCollector& dataIn, 
                   real tolIn,
                   int maxIts):
  error(true),
  USE_GUESS(false),
  maxIterations(maxIts),
  tol(tolIn), 
  rhoK(12345),
  alpha(12345),
  beta(12345),
  tauKm1(12345),
  tauKm2(12345),
  bNorm(12345),
  pDotw(12345),
  zNorm(12345),
  xNorm(12345),
  r(linOpIn.dimDomain(),12345),
  p(linOpIn.dimDomain(),12345),
  Ax(linOpIn.dimDomain(),12345),
  z(linOpIn.dimDomain(),12345),
  w(linOpIn.dimDomain(),12345),
  x(this,linOpIn.dimDomain()),
  b(this,linOpIn.dimDomain()),
  vectorNorm(&normIn),
  data(&dataIn),
  linearOperator(&linOpIn),
  preconditioner(&precIn)
{  
  Tracer tr("PreCG::PreCG(linOp& ...)");
}

PreCG::~PreCG()
{
  Tracer tr("PreCG::~PreCG()");
}

bool PreCG::prepare()
{
  return preconditioner->prepare();
}

//just took simple version from Tim Kelley's book for now
bool PreCG::solve(const Vec& bIn, Vec& xIn)
{ 
  b.attachToTarget(bIn);
  x.attachToTarget(xIn);

  if (!USE_GUESS)
    x.v_=0.0;
  real tolFinal;
  //1. Initialize.  The real system is A_0 x = b_0
  k=0;

  linearOperator->apply(x.v_,Ax);               //A_0 x
#ifndef USE_BLAS
  r = b.v_ - Ax;                                //b - A_0 x  
#else
  copy(b.v_,r);
  axpy(-1.0,Ax,r);
#endif
  rhoK = dot(r,r);                           //||r||_2^2
  tolFinal=sqrt(rhoK)*tol;
  k=0;

  xNorm = (*vectorNorm)(x.v_);

  while(sqrt(rhoK) >= tolFinal && k < maxIterations )
    {
      k++;
      data->linearSolverIteration();

      error = preconditioner->apply(r,z);         //z = Mr
      if (error) return 1;
      
      tauKm2 = tauKm1;
      tauKm1 = dot(z,r);
      if (k==1)
        {
          beta=0.0;
#ifndef USE_BLAS
          p=z;
#else
          copy(z,p);
#endif
        }
      else
        {
          beta = tauKm1/tauKm2;
#ifndef USE_BLAS
          p=z+beta*p;
#else
          scal(beta,p);
          axpy(1.0,z,p);
#endif
        }

      linearOperator->apply(p,w);
      
      alpha = tauKm1/dot(p,w);
      
      if (alpha <= MACHINE_EPSILON*xNorm)
        {
          cerr<<"PreCG: Nonpositive curvature, alpha = "<<alpha<<endl
              <<"PreCG: Linear System may not be positive definite or"<<endl
              <<"PreCG: machine precision has been reached without"<<endl
              <<"PreCG: convergence"<<endl;
          return 1;
        }
#ifndef USE_BLAS
      x.v_  += alpha*p;
      r  -= alpha*w;
#else
      axpy(alpha,p,x.v_);
      axpy(-alpha,w,r);
#endif
      rhoK= dot(r,r);
      //cerr<<sqrt(rhoK)<<endl;
      xNorm=(*vectorNorm)(x.v_);
    }

  if(k == maxIterations)
    {
#ifdef DEBUG
      cerr <<"failed with max PreCG iterations exceeded, used ="
	   <<k<<endl;
#endif
      return 1;
    }
  else
    {
#ifdef DEBUG
      cerr <<"success with  "<<k<<" PreCG iterations"<<endl;
#endif
      if (SOLVE_SUB)
        xIn=bIn;
      x.restoreToTarget();      
      return 0;
    }
}

}//Daetk
