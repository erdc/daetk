#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Definitions.h"
#include "Utilities.h"
#include "Vec.h"

namespace Daetk 
{
class Integrator 
{
public:
  Integrator();

  virtual ~Integrator();

  virtual bool calculateSolution(const real& tout,Vec& solution, 
                                 Vec& solutionPrime)=0;
  virtual bool step(const real& tOut,real& tStep,Vec& solutionAtTStep, 
                    Vec& solutionPrime)=0;
  virtual void reset();
};
}//Daetk
#endif
