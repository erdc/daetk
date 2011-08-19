#ifndef SLCOMPRESS_H
#define SLCOMPRESS_H

#include "Definitions.h"
#include "Utilities.h"
#include "ParameterDatabase.h"
#include "Density.h"

namespace Daetk 
{
  class SlCompress : public Density
  {
  public:
    inline SlCompress(const SlCompress& slc);
    SlCompress(real beta,real rho, real gravity, real rhoBase);
    
    inline void setHead(const real& psi);
    void setHead(const Vec& p, Vec& r, Vec& Dr);
    void setHeadD(const Vec& p, Vec& r, Vec& Dr, Vec& DDr);
    
    inline void setPressure(const real& p);
    
    inline const real& getRho();
    
    inline const real& getDrho();
    
    inline const real& getDDrho();
    
    real p_,rho_,Drho_,DDrho_,beta_,rho0_,g_,brg;
  };
  
  inline SlCompress::SlCompress(const SlCompress& slc):
    p_(slc.p_),
    rho_(slc.rho_),
    Drho_(slc.Drho_),
    DDrho_(slc.DDrho_),
    beta_(slc.beta_),
    rho0_(slc.rho0_),
    g_(slc.g_),
    brg(slc.brg)
  {}
  
  inline void SlCompress::setHead(const real& psi)
  {
    //already dividec by rhoBase
    rho_ = rho0_*exp(brg*psi);
    Drho_ = brg*rho_;
  }
  
  inline void SlCompress::setPressure(const real& p)
  {
    //already divided by rhoBase
    p_ = p;
    rho_ = rho0_*exp(beta_*p_);
    Drho_ = beta_*rho_;
  }
  
  inline const real& SlCompress::getRho()
  {
    return rho_;
  }
  
  inline const real& SlCompress::getDrho()
  {
    return Drho_;
  }
  
  inline const real& SlCompress::getDDrho()
  {
    return DDrho_=brg*Drho_;
  }
  
}//Daetk
#endif

