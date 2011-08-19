#include "GasDensity.h"
namespace Daetk
{
GasDensity::GasDensity(real alpha,real rho, real g, real rhoBase)
{
  if (alpha != 0)
    {
      alpha_=alpha/rhoBase;
      alphaBar_=alpha*g;
      stdPressure_=rho/alpha;
      rho0_=rho/rhoBase;
    }
  else //incompressible
    {
      alpha_=alpha/rhoBase;
      alphaBar_=alpha*g;
      stdPressure_=0.0;
      rho0_=rho/rhoBase;
    }        
}

GasDensity::~GasDensity(){}


  void GasDensity::setHead(const Vec& p, Vec& r, Vec& Dr)
  {
    for (int i=0;i<p.ldim();i++)
      {
        r[i] = rho0_ + alphaBar_*p[i];
        Dr[i] = alphaBar_;
      }
  }
  
  void GasDensity::setHeadD(const Vec& p, Vec& r, Vec& Dr, Vec& DDr)
  {
    for (int i=0;i<p.ldim();i++)
      {
        r[i] = rho0_ + alphaBar_*p[i];
        Dr[i] = alphaBar_;
        DDr[i] = 0.0;
      }
  }

}//Daetk






