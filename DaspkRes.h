#ifndef DASPKRES_H
#define DASPKRES_H

#include "Definitions.h"
#include "DaeDefinition.h"

namespace Daetk 
{
class DaspkResidual
{
public:
  DaspkResidual();
  virtual ~DaspkResidual();
  static void res(const real& t, real* y, real* yp, real& cj,
             real* delta, int* ires, real* rpar, int *ipar);
  static DaeDefinition* theDae;
};
}//Daetk
#endif
