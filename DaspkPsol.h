#ifndef DASPKPSOL_H
#define DASPKPSOL_H

#include  "Definitions.h"
#include "Preconditioner.h"

namespace Daetk 
{
class DaspkPsol
{
public:
  DaspkPsol(){}
  virtual ~DaspkPsol();
  static void psol(const int& neq, const real& t, real* y, 
                          real* yp, real* savr,real* wk,
                          const real& cj,real* wght,real * wp,
                          int* iwp, real* b, const real& eplin, int* ier, 
                          real* rpar,int* ipar);

  static Preconditioner* thePrec;
protected:
  static bool firstCall;
  static Vec bVec,b2;
};
}//Daetk
#endif
