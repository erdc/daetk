#ifndef Poisson2p_H
#define Poisson2p_H

#include "Definitions.h"
#include "ParameterDatabase.h"
#include "Psk2p.h"

namespace Daetk 
{
class Poisson2p : public Psk2p
{
public:
  Poisson2p();
  virtual ~Poisson2p();
  void readParameters(ParameterDatabase& pd);
  
  inline void setHeads(int node, const real psiWIn, const real psiNIn=0);
  
  inline void setVFraction(const real& thetaWIn, int node);
  
  inline void calculateDerivatives();

  void millerSimilarScaling(Vec& delta);
};

  inline void Poisson2p::setHeads(int node, const real psiWIn, const real psiNIn)
  {
    psiC = psiNIn - psiWIn;
    sBar = thetaS[node] - 0.5* psiC;
    thetaW = sBar;
    DsBar_DpC = -0.5;
    DthetaW_DpC = -0.5;
    krW = 1.0;
    KW = krW;
    krN = 1.0;
    KN = krN;  
  }
  
inline  void Poisson2p::setVFraction(const real& thetaWIn, int node)
  {
    psiC = 2*(thetaS[node] - thetaWIn);
    sBar = thetaWIn;
    DsBar_DpC = -0.5;
    DthetaW_DpC = -0.5;
    DsBar_DpC = 0.0;
    DthetaW_DpC = 0.0;
    krW = 1.0;
    KW = krW;
    krN = 1.0;
    KN = krN;  
  }  
  
inline void Poisson2p::calculateDerivatives()
  {
      DDsBar_DDpC = 0.0;
      DDthetaW_DDpC = 0.0;
      DkrW_DpC = 0.0;
      DKW_DpC = 0.0;
      DkrN_DpC = 0.0;
      DKN_DpC = 0.0;
    }

}//Daetk
#endif
