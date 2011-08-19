#ifndef DIVERGENCE
#define DIVERGENCE

#include "Definitions.h"
#include "Vec.h"

namespace Daetk
{
class Divergence
{
public:
  virtual const Vec& getFlux_x(int n=0)=0;
  virtual const Vec& getFlux_y(int n=0)=0;
  virtual const Vec& getFlux_z(int n=0)=0;
  virtual Vec& setFlux_x(int n=0)=0;
  virtual Vec& setFlux_y(int n=0)=0;
  virtual Vec& setFlux_z(int n=0)=0;
  virtual void computeInterfaceK(const Vec& K)=0;
};
}//Daetk
#endif
