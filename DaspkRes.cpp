#include "DaspkRes.h"
namespace Daetk 
{

DaeDefinition* DaspkResidual::theDae = 0;

DaspkResidual::DaspkResidual(){}

DaspkResidual::~DaspkResidual(){}

void DaspkResidual::res(const real& t, real* y, 
                                 real* yp, real& cj, real* delta, 
                                 int* ires, real* rpar, int *ipar)
{
  int neq = theDae->getY0().dim();
  const Vec yVec(Vec::REF,y,neq,0,1);
  const Vec ypVec(Vec::REF,yp,neq,0,1);
  Vec resVec(Vec::REF,delta,neq,0,1);
  *ires = theDae->residual(t,yVec,ypVec,resVec);
}

}//Daetk
