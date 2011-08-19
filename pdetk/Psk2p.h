#ifndef PSK2P_H
#define PSK2P_H

#include "Definitions.h"
#include "ParameterDatabase.h"
#include "VecBlas.h"
#include "VecOperators.h"
#include "PetscSecondOrderFd.h"

namespace Daetk 
{
class Psk2p
{
public:
  //mwf added reference to argument
  Psk2p(ParameterDatabase& pd);

  Psk2p();

  virtual ~Psk2p();
  virtual void setInitialConditions(const Vec& pw){}
  virtual void updateHistory(){}
  virtual void readParameters(ParameterDatabase& pd);
  virtual void readZones(ParameterDatabase& pd, Petsc::SecondOrderFd& node){};

  virtual void setHeads(int node,real psiWIn, real psiNIn=0)=0;
  virtual void setVFraction(const real& thetaWIn, int node)=0;
  virtual void calculateDerivatives()=0;
  virtual void calculateDerivativesVFraction()=0;
  
  //*****vector utilities
  
  //use these to tell psk which values to calculate
  void attachPsiC(Vec& psiCIn);

  void attachThetaW(Vec& thetaWIn);
  void attachDthetaW_DpC(Vec& DthetaWIn);
  void attachDDthetaW_DDpC(Vec& DDthetaWIn);

  void attachKrW(Vec& krwIn);
  void attachDKrW_DpC(Vec& DkrwIn);

  void attachKrN(Vec& krnIn);
  void attachDKrN_DpC(Vec& DkrnIn);

  //use these to calculate the attached vectors for the given vector arguments
  //the entire vector is set and all vectors must be the same layout
  virtual void setHeads(const Vec& psiW_vec, const Vec& psiN_vec);
  virtual void setVFraction(const Vec& thetaW_vec);
  virtual void calculateDerivativesHead(const Vec& psiW_vec, const Vec& psiN_vec);
  virtual void calculateDerivativesVFraction(const Vec& thetaW_vec);

  //*****vector utilities

  virtual void setAlpha(Vec &alphaIn){}
  virtual void setN (Vec &nIn){}
  virtual void setThetaS (Vec &thetaSIn);
  virtual void setThetaR (Vec &thetaRIn);
  virtual void setKWs(Vec& KWsIn);
  virtual void millerSimilarScaling(Vec& delta)=0;

  inline const real&   getPsiC(){return psiC;}
  inline const real&   getThetaW(){return thetaW;}
  inline const real&   getDthetaW_DpC(){return DthetaW_DpC;}
  inline const real&   getDDthetaW_DDpC(){return DDthetaW_DDpC;}
  inline const real&   getSbar(){return sBar;}
  inline const real&   getDsBar_DpC(){return DsBar_DpC;}
  inline const real&   getKW(){return KW;}
  inline const real&   getKN(){return KN;}
  inline const real&   getDKW_DpC(){return DKW_DpC;}
  inline const real&   getDKN_DpC(){return DKN_DpC;}
  inline const real&   getKrW(){return krW;}
  inline const real&   getKrN(){return krN;}
  inline const real&   getDkrW_DpC(){return DkrW_DpC;}
  inline const real&   getDkrN_DpC(){return DkrN_DpC;}
  inline const Vec&   getThetaS(){return global_thetaS;}
  inline const Vec&   getThetaR(){return global_thetaR;}
  inline const Vec&   getKWs(){return global_KWs;}
  inline const Vec&   getThetaSR(){return global_thetaSR;}
  inline const Vec&   getLocalThetaS(){return thetaS;}
  inline const Vec&   getLocalThetaR(){return thetaR;}
  inline const Vec&   getLocalKWs(){return KWs;}
  inline const Vec&   getLocalThetaSR(){return thetaSR;}

protected:
  int nNodes,i;
  
  real psiC,
    thetaW,
    DthetaW_DpC,
    DDthetaW_DDpC,
    KW,KN,DKW_DpC,DKN_DpC,
    krW,krN,DkrW_DpC,DkrN_DpC,    
    sBar, DsBar_DpC, DDsBar_DDpC,
    muW_by_muN;

  Vec global_thetaS,global_thetaR,global_KWs,global_thetaSR;
  Vec thetaS,thetaR,KWs,thetaSR;
  Vec *psiC_p, 
    *thetaW_p,*DthetaW_DpC_p,*DDthetaW_DDpC_p,
    *krW_p,*DkrW_DpC_p,
    *krN_p,*DkrN_DpC_p;
};

}//Daetk
#endif

