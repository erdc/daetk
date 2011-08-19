#ifndef GASDENSITY_H
#define GASDENSITY_H
#include "Definitions.h"
#include "Utilities.h"
#include "ParameterDatabase.h"
#include "Density.h"

namespace Daetk 
{
  class GasDensity : public Density
{
public:
  inline GasDensity(const GasDensity& gd);
  GasDensity(real alpha,real rho, real g, real rhoBase);
 
  virtual ~GasDensity();

  inline void setHead(const real& psi);
  void setHead(const Vec& p, Vec& r, Vec& Dr);
  void setHeadD(const Vec& p, Vec& r, Vec& Dr, Vec& DDr);
  inline void setPressure(const real& p);
  inline const real& getRho();
  inline const real& getDrho();
  inline const real& getDDrho();
  real rho_,Drho_,DDrho_,p_,alpha_,alphaBar_,stdPressure_,rho0_;
};

inline GasDensity::GasDensity(const GasDensity& gd):
  rho_(gd.rho_),          
  Drho_(gd.Drho_),        
  DDrho_(gd.DDrho_),       
  p_(gd.p_),           
  alpha_(gd.alpha_),       
  alphaBar_(gd.alphaBar_),    
  stdPressure_(gd.stdPressure_), 
  rho0_(gd.rho0_)         
  {}
inline void GasDensity::setHead(const real& psi)
{
  //already divided by rhoBase
  rho_ = rho0_ + alphaBar_*psi;
  Drho_ = alphaBar_;
}

inline void GasDensity::setPressure(const real& p)
{
  //already divided by rhoBase
  if (alpha_)
    {
      p_ = p + stdPressure_;
      rho_ = alpha_*p_;
      Drho_ = alpha_;
    }
  else
    {
      rho_ = rho0_;
      Drho_ = 0.0;
    }
}

inline const real& GasDensity::getRho()
{
  return rho_;
}

inline const real& GasDensity::getDrho()
{
  return Drho_;
}

inline const real& GasDensity::getDDrho()
{
  return DDrho_=0.0;
}

}//Daetk
#endif
