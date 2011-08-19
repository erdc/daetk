#include "DaspkJac.h"

namespace Daetk 
{

Jacobian* DaspkJacobian::theJac = 0;
Preconditioner* DaspkJacobian::thePrec = 0;
DaeDefinition* DaspkJacobian::theDae = 0;

void DaspkJacobian::jac_direct(real& t, real* y, real* yprime, 
                         real* pd, real& cj, real* rpar, int* ipar)
{
  //won't work with difference jacobian until I do something about the pd matrix
  int Neq = theDae->getY0().dim();
  const Vec yVec(Vec::REF,y,Neq,0,1);
  const Vec ypVec(Vec::REF,yprime,Neq,0,1);

  theDae->tDaeDef = t;
  theDae->alphaDaeDef = cj;
#ifndef USE_BLAS
  theDae->yDaeDef = yVec;
  theDae->ypDaeDef = ypVec;
  theDae->betaDaeDef = ypVec - theDae->alphaDaeDef*yVec;
#else
  copy(yVec,theDae->yDaeDef);
  copy(ypVec,theDae->ypDaeDef);
  copy(ypVec,theDae->betaDaeDef);
  axpy(-theDae->alphaDaeDef,yVec,theDae->betaDaeDef);
#endif
  theDae->resetFunction();
  theJac->evaluate(yVec,yVec);  
}

void DaspkJacobian::jac_krylov(GlobDaspkRes res,int* ires,int* neq,  real& t,  real* y,  real* yprime,
                         real* weight,  real* residual,real& work, 
                         real& h,  real& cj, 
                        real* wp, int* iwp,  int& ier, real* rpar,int* ipar)  
{
  int Neq = theDae->getY0().dim();
  const Vec yVec(Vec::REF,y,Neq,0,1);
  const Vec ypVec(Vec::REF,yprime,Neq,0,1);
  const Vec rVec(Vec::REF,residual,Neq,0,1);

  theDae->tDaeDef = t;
  theDae->alphaDaeDef = cj;
#ifndef USE_BLAS
  theDae->yDaeDef = yVec;
  theDae->ypDaeDef = ypVec;
  theDae->betaDaeDef = ypVec - theDae->alphaDaeDef*yVec;
#else
  copy(yVec,theDae->yDaeDef);
  copy(ypVec,theDae->ypDaeDef);
  copy(ypVec,theDae->betaDaeDef);
  axpy(-theDae->alphaDaeDef,yVec,theDae->betaDaeDef);
#endif
  theDae->resetFunction();
  theJac->evaluate(yVec,rVec);  
  thePrec->prepare();
}

}//Daetk
