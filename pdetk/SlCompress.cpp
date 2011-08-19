#include "SlCompress.h"

namespace Daetk
{
  SlCompress::SlCompress(real beta,real rho, real gravity, real rhoBase):
    beta_(beta),
    rho0_(rho/rhoBase),
    g_(gravity),
    brg(beta*rhoBase*gravity)
    {}

 void SlCompress::setHead(const Vec& p, Vec& r, Vec& Dr)
    {
      const int end=p.getLocalHigh();
      for (int i=0;i<end;i++)
	{
	  r[i] = rho0_*exp(brg*p[i]);
	  Dr[i] = brg*r[i];
	}
    }

  void SlCompress::setHeadD(const Vec& p, Vec& r, Vec& Dr, Vec& DDr)
    {
      const int end=p.getLocalHigh();
      for (int i=0;i<end;i++)
	{
	  r[i] = rho0_*exp(brg*p[i]);
	  Dr[i] = brg*r[i];
	  DDr[i]= brg*Dr[i];
	}
    }
}//Daetk


