#ifndef BROOKSCOREY2P_H
#define BROOKSCOREY2P_H

#include "Definitions.h"
#include "Psk2p.h"
#include "Psk2pSpline.h"
#include "ParameterDatabase.h"
#include "PetscSecondOrderFd.h"
#include "VecOperators.h"

namespace Daetk 
{
class BrooksCorey2p : public Psk2p
{
public:
  BrooksCorey2p();
  virtual ~BrooksCorey2p();
  void readParameters(ParameterDatabase& pd);
  void readZones(ParameterDatabase& pd, Petsc::SecondOrderFd& node);
  inline void setHeads(int node,real psiWIn, real psiNIn=0);
  inline void setVFraction(const real& thetaWIn, int node);
  inline void calculateDerivatives();
  inline void calculateDerivativesVFraction();
  virtual void setHeads(const Vec& psiW_vec, const Vec& psiN_vec);
  virtual void setVFraction(const Vec& thetaW_vec);
  virtual void calculateDerivativesHead(const Vec& psiW_vec, const Vec& psiN_vec);
  virtual void calculateDerivativesVFraction(const Vec& thetaW_vec);
  void millerSimilarScaling(Vec& delta);
  friend class Psk2pSpline<BrooksCorey2p>;
private:
  Vec global_psiD,global_lambda;
  Vec psiD,lambda;
};

inline void BrooksCorey2p::setHeads(int node,real psiWIn, real psiNIn)
{
  psiC = psiNIn - psiWIn;
  i = node;
  
  sBar = pow(psiC/psiD[i],-lambda[i]);
  thetaW = thetaSR[i]*sBar + thetaR[i];
  DsBar_DpC = -(lambda[i]/psiD[i])*pow(psiC/psiD[i],-lambda[i]-1);
  DthetaW_DpC = thetaSR[i] * DsBar_DpC; 
  krW = pow(sBar,(2.+3.*lambda[i])/lambda[i]);
  KW = KWs[i]*krW;
  krN = (1.-sBar)*(1-sBar)*(1.-pow(sBar,(2.0 + lambda[i])/lambda[i]));
  KN = KWs[i]*muW_by_muN*krN;
  
  if(psiC < psiD[i]) 
    {
      sBar = 1.0;
      thetaW = thetaS[i];
      DsBar_DpC = 0.0;
      DthetaW_DpC = 0.0;
      krW = 1.0;
      KW = KWs[i];
      krN = 0.0;
      KN = 0.0;
    }    
}

inline void BrooksCorey2p::setVFraction(const real& thetaWIn, int node)
{
  i = node;
  sBar = (thetaWIn - thetaR[i])/thetaSR[i];
  sBar= std::min(std::max(sBar,0.0),1.0);
  thetaW = thetaSR[i]*sBar + thetaR[i];
  psiC = psiD[i]*pow(sBar,-1./lambda[i]);
  //std::cout<<psiC<<'\t'<<psiD[i]<<'\t'<<lambda[i]<<'\t'<<sBar<<'\t'<<i<<std::endl;
  DsBar_DpC = -(lambda[i]/psiD[i])*pow(psiC/psiD[i],-lambda[i]-1.);
  DthetaW_DpC = thetaSR[i] * DsBar_DpC; 
  krW = pow(sBar,(2.+3.*lambda[i])/lambda[i]);
  KW = KWs[i]*krW;
  krN = (1.0-sBar)*(1.0-sBar)*(1.0-pow(sBar,(2.+lambda[i])/lambda[i]));
  KN = KWs[i]*muW_by_muN*krN;
  
  if (psiC < psiD[i])
    {
      DsBar_DpC = 0.0;
      DthetaW_DpC = 0.0;
      krW = 1.0;
      KW = KWs[i];
      krN = 0.0;
      KN = 0.0;
    }
}
  
inline void BrooksCorey2p::calculateDerivatives()
{
  //probably not right
  DDsBar_DDpC = -(lambda[i]/psiD[i])*((-lambda[i]-1.)/psiD[i])*pow(psiC/psiD[i],-lambda[i]-2.);
  DDthetaW_DDpC = thetaSR[i]*DDsBar_DDpC;
  DkrW_DpC = ((2.+3.*lambda[i])/lambda[i])*pow(sBar,(2.+3.*lambda[i])/lambda[i]-1.)*DsBar_DpC;
  DKW_DpC = KWs[i]*DkrW_DpC;
  //mwf should this be different, second term pow should be 2/lambda[i]
  //mwf was
  //DkrN_DpC =-2.0*(1.-sBar)*DsBar_DpC*(1.-pow(sBar,(2.+lambda[i])/lambda[i])) 
  // - (1.-sBar)*(1.-sBar)*pow(sBar,(2.+lambda[i])/lambda[i])*(2.+lambda[i])/lambda[i]*DsBar_DpC;

  //mwf try
  DkrN_DpC =-2.0*(1.-sBar)*DsBar_DpC*(1.-pow(sBar,(2.+lambda[i])/lambda[i])) 
   -(1.-sBar)*(1.-sBar)*pow(sBar,2./lambda[i])*(2.+lambda[i])/lambda[i]
    *DsBar_DpC;

  DKN_DpC = KWs[i]*muW_by_muN*DkrN_DpC;

  if (psiC < psiD[i])
    {
      DDsBar_DDpC = 0.0;
      DDthetaW_DDpC = 0.0;
      DkrW_DpC = 0.0;
      DKW_DpC = 0.0;
      DkrN_DpC = 0.0;
      DKN_DpC = 0.0;
    }
}
inline void BrooksCorey2p::calculateDerivativesVFraction()
{
  calculateDerivatives();
}

}//Daetk
#endif






