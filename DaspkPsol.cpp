#include "DaspkPsol.h"

namespace Daetk 
{
  
DaspkPsol::~DaspkPsol()
{
  bVec.clear();
  b2.clear();
  firstCall=true;
}

Preconditioner* DaspkPsol::thePrec = 0;

Vec DaspkPsol::bVec;
Vec DaspkPsol::b2;

bool DaspkPsol::firstCall(true);

void DaspkPsol::psol(const int& neq, const real& t, real* y, 
                          real* yp, real* savr,real* wk,
                          const real& cj,real* wght,real * wp,
                          int* iwp, real* b, const real& eplin, int* ier, 
                          real* rpar,int* ipar)
{
  if (firstCall)
    {
      b2.newsize(neq);
      firstCall=false;
    }
  bVec.attachToArray(b,neq,0,1);
  b2 = bVec;
  *ier = thePrec->apply(b2,bVec);
}
   
}//Daetk
