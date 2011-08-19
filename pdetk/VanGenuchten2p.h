#ifndef VANGENUCHTEN2P_H
#define VANGENUCHTEN2P_H

#include "Definitions.h"
#include "Psk2p.h"
#include "ParameterDatabase.h"
#include "VecBlas.h"
#include "VecOperators.h"
#include <stdlib.h>
#include <math.h>
namespace Daetk 
{
class VanGenuchten2p : public Psk2p
{
public:
  VanGenuchten2p();
  virtual ~VanGenuchten2p();
  void readParameters(ParameterDatabase& pd);
  inline void setCapillaryHead(int node,real psiCIn);
  virtual void setHeads(int node,real psiWIn, real psiNIn=0){setCapillaryHead(node,psiNIn-psiWIn);}
  inline void setVFraction(const real& thetaWIn, int node);
  inline void calculateDerivatives();
  inline void calculateDerivativesVFraction(){calculateDerivatives();}
  void setKWs (Vec &KWsIn);
  void setAlpha(Vec &alphaIn);
  void setN (Vec &nIn);
  void setThetaS (Vec &thetaSIn);
  void setThetaR (Vec &thetaRIn);
  void millerSimilarScaling(Vec& delta);
protected:
  real vBar,uBar,
    alphaPsiC, alphaPsiC_n, alphaPsiC_nM1, alphaPsiC_nM2,
    onePlus_alphaPsiC_n,
    sqrt_sBar, sqrt_1minusSbar,
    sBarByOnePlus_alphaPsiC_n, sBarBy_onePlus_alphaPsiC_n_2;
  Vec global_n,global_m,global_alpha;
  Vec n,m,alpha;
};
  
  inline void VanGenuchten2p::setCapillaryHead(int node,real psiCIn)
{
  psiC = psiCIn;
  i = node; //local node#
  
  alphaPsiC = alpha[i]*psiC;
  alphaPsiC_n = pow(alphaPsiC,n[i]);
  alphaPsiC_nM1 = alphaPsiC_n/alphaPsiC;
  onePlus_alphaPsiC_n = 1.0 + alphaPsiC_n;
  sBar = pow(onePlus_alphaPsiC_n,-m[i]);
  sBarByOnePlus_alphaPsiC_n = sBar/onePlus_alphaPsiC_n;
  DsBar_DpC = -alpha[i]*(n[i]-1.0)*alphaPsiC_nM1 
    *sBarByOnePlus_alphaPsiC_n;
  if(psiC<=0.0) 
    {
      sBar = 1.0;
      DsBar_DpC = 0.0;
    }    
}

inline void VanGenuchten2p::setVFraction(const real& thetaWIn, int node)
{
  i = node;
  sBar = (thetaWIn - thetaR[i])/thetaSR[i];

  sBar= std::min(std::max(sBar,0.0),1.0);
  thetaW = thetaSR[i]*sBar + thetaR[i];

  onePlus_alphaPsiC_n = pow(sBar,1.0/-m[i]);
  alphaPsiC_n = onePlus_alphaPsiC_n - 1.0;
  alphaPsiC = pow(alphaPsiC_n,1.0/n[i]);
  psiC = alphaPsiC/alpha[i];
  
  alphaPsiC_nM1 = alphaPsiC_n/alphaPsiC;
  sBarByOnePlus_alphaPsiC_n = sBar/onePlus_alphaPsiC_n;
  DsBar_DpC = -alpha[i]*(n[i]-1.0)*alphaPsiC_nM1 
    *sBarByOnePlus_alphaPsiC_n;

  if(psiC<=0.0) 
    {
      DsBar_DpC = 0.0;
    }
}


inline void VanGenuchten2p::calculateDerivatives()
{
  alphaPsiC_nM2 =   alphaPsiC_nM1/alphaPsiC;      
  
  sBarBy_onePlus_alphaPsiC_n_2 = sBarByOnePlus_alphaPsiC_n
    /onePlus_alphaPsiC_n;
  DDsBar_DDpC =  alpha[i]*alpha[i]*(n[i]-1)
    *((2*n[i]-1)*alphaPsiC_nM1*alphaPsiC_nM1
      *sBarBy_onePlus_alphaPsiC_n_2
      -
      (n[i]-1)*alphaPsiC_nM2
      *sBarByOnePlus_alphaPsiC_n);

  if (psiC <= 0.0)
    {
      DDsBar_DDpC = 0.0;
    }
}

}//Daetk
#endif





