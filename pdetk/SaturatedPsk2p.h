#ifndef SATURATEDPSK2P_H
#define SATURATEDPSK2P_H

#include "Definitions.h"
#include "Psk2p.h"
#include "Psk2pSpline.h"
#include "ParameterDatabase.h"
#include "PetscSecondOrderFd.h"
#include "VecOperators.h"

namespace Daetk 
{
class SaturatedPsk2p : public Psk2p
{
public:
  SaturatedPsk2p();
  virtual ~SaturatedPsk2p();
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
};

inline void SaturatedPsk2p::setHeads(int node,real psiWIn, real psiNIn)
{
  psiC = 0.0;
  i = node;
  sBar = 1.0;
  thetaW = thetaS[i];
  DsBar_DpC = 0.0;
  DthetaW_DpC = 0.0;
  krW = 1.0;
  KW = KWs[i];
  krN = 0.0;
  KN = 0.0;
}

inline void SaturatedPsk2p::setVFraction(const real& thetaWIn, int node)
{
  i = node;
  sBar = 1.0;
  psiC = 0.0;
  DsBar_DpC = 0.0;
  DthetaW_DpC = 0.0;
  krW = 1.0;
  KW = KWs[i];
  krN = 0.0;
  KN = 0.0;
}
  
inline void SaturatedPsk2p::calculateDerivatives()
{
  DDsBar_DDpC = 0.0;
  DDthetaW_DDpC = 0.0;
  DkrW_DpC = 0.0;
  DKW_DpC = 0.0;
  DkrN_DpC = 0.0;
  DKN_DpC = 0.0;
}

inline void SaturatedPsk2p::calculateDerivativesVFraction()
{
  calculateDerivatives();
}

}//Daetk
#endif






