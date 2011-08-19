#ifndef GERHARD_KUEPER2P_H
#define GERHARD_KUEPER2P_H

#include "Definitions.h"
#include "Psk2p.h"
#include "Psk2pSpline.h"
#include "ParameterDatabase.h"
#include "PetscSecondOrderFd.h"
#include "VecOperators.h"

namespace Daetk 
{
class GerhardKueper2p : public Psk2p
{
  /***********************************************************************
    implement a subset of the Gerhard and Kueper psk relations from
    their WRR 2003 papers. I am interested initially only in getting a
    reasonable estimate of nonwetting phase residual, rather than full
    hysteresis. 


    To incorporate residual formation or "trapping," they use a linear
    relationship between the reversal saturation (from drainage to
    imbibition) and the extinctiction saturation.

      S^X_W = (1-S^{X,\ast}_W)(S^'_S-1) + 1

      S^'_W : current reversal point from drainage to imbibition which
              I will define as the minimum wetting phase saturation
              'seen' by a location up to the current time.

      S^X_W : extinction saturation, wetting phase saturation when
              nonwetting phase becomes immobile.
 
      S^{X,\ast}_W : wetting phase saturation at maximum complete
                     nonwetting phase residual
  

    To implement the residual "trapping," the only extra parameter
    that is needed is S^{X,\ast}_W. I will attempt to lag the
    calculation of S^'_W a time step, so that the constitutive
    relations and derivatives are still essentially
    Brooks-Corey. Specifically, I'll use

      P_C = P_DS_e^{-1/lambda_d}, P_N >= P_D
 
      S_e = (S_W-S_{W,r})/(S^X_W-S_{W,r}) where
         
             S^{X,\ast}_W <= S^X_W <= 1.0 and
             S^X_W = 0 for primary drainage (S^{min}_{W} = 1.0)

      P_D = displacement pressure (entry pressure at S_W=1?)
      \lambda_d = drainage pore-size distribution index


    For secondary imbibition, the basic capillary pressure relationship is

      P_C = P_T(\tilde{S}_e)^{-1/lambda_i}

    for
     
      \tilde{S}_e = (S_W - S^f_W)/(S^X_W-S^f_W)

    where

      P_T  : terminal capillary pressure (at S^X_W). 
       
      lambda_i : imbibition pore-size distribution index

      S^f_W : wetting phase saturation value that the imbition curve
              approaches at infinite capillary pressure. For main
              imbibition, S^f_W = S_r, and S^X_W = S^{X,\ast}_W. In
              general, Gerhard and Kueper define it so that the
              imbibition and drainage curves intersect at a reversal
              point. This requires solving a local nonlinear problem
              according to Gerhard and Kueper, but I think in some
              cases at least you can solve it directly (see ghPSK.py).


    Wetting phase relative permeability is the standard brooks corey
    relationship

      k_{r,W} = S_e^{(2+3\lambda_d)/\lambda_d}

         S_e = (S_W-S^k_r)/(1.0 - S^k_r)

    The nonwetting phase relative permeability is a little more
    complicated. The full Gerhard and Kueper model for drainage is

      k_{r,N} = k^{max}_{r,N}(1-\bar{S}^{\ast,d}_e)^{2\tau_d}
                  (1-\bar{S}^{\ast,d}_e^{(2+\lambda_d)/\lambda_d}

      where
          k^{max}_{r,N} : maximum nonwetting phase permeability

          \tau_d : nonwetting phase relative tortuosity coefficient
          
          \bar{S}^{\ast,d}_e : effective saturation for drainage
                 = (S_W-S^k_r)/[(S^M_W + \Delta S^{\ast,d}_W) - S^k_r]
                   for S^k_r <= S_W <= S^M_W
                 = 1 for S^M_W < S_W <= 1
           
          \S^M_W : emergence saturation, wetting phase saturation at
                   which nonwetting phase flow begins (k_{r,N} >0). 
                   Corresponds to entry pressure
          
          \Delta S^{\ast,d}_W : jump fitting parameter to allow
                   nonwetting phase relative permeability to go from 0
                   to a finite positive value at the emergence
                   saturation, S^{M}_W

     For secondary imbibition, Gerhard and Kueper use

       k_{r,N} = k^'_{r,N}(1-\hat{S}^{\ast,i}_e)^{2\tau_i}
                   (1-\hat{S}^{\ast,i}_e^{(2+\lambda_i)/\lambda_i}

     where
       
        k^'_{r,N} : nonwetting phase relative permeability at reversal 
                    point

        \tau_i : nonwetting phase relative tortuosity coefficient for
             imbibition

        \hat{S}^{\ast,i}_e : effective saturation for imbibition

               = (S_W-S^'_W)/[(S^X_W + \Delta S^{\ast,i}_W) - S^'_W]
                   for S^'_W <= S_W <= S^X_W
               = 1 for S^X_W = S_W

        S^'_W : wetting phase saturation at reversal point from
                drainage to imbibition
    
        \Delta S^{\ast,i}_W : jump fitting parameter to allow
                   nonwetting phase relative permeability to go from 
                   a finite positive value at the extinction
                   saturation, S^{X}_W, to zero    

     Their are a number of new parameters that the model
     introduces. Fortunately, there are defaults that reduce almost
     everything back to Brooks Corey.

       \tau_d = \tau_i = 1.0

       \Delta S^{\ast,d}_W = \Delta S^{\ast,i}_W = 0

       k^{max}_{r,N} = 1.0

       S^M_W = 1.0

     If I ignore hysteresis, I can drop the distinctions between
       \lambda_i --- \lambda_d,  and P_D -- P_T

  ----------------------------------------------------------------------
  To summarize, the initial version of these psk relations will use


     k_{r,W} = S_e^{(2+3\lambda)/\lambda}

         S_e = (S_W-S^k_r)/(1.0 - S^k_r)


     k_{r,N} = k^{max}_{r,N}(1-S^{\ast}_e)^{2}
                (1-S^{\ast}_e^{(2+\lambda)/\lambda}

     S^{\ast}_e = (S_W-S^k_r)/(S^{\ast}_W-S^k_r)

     S^{\ast}_W : wetting phase saturation at which nonwetting phase
                  flow initiates or ceases. In other words, it
                  corresponds to the emergence saturation on drainage
                  and extinction saturation on imbibition

     S^{\ast}_W = (1. - S^{X,\ast}_{W})(S^'_W - 1) + 1

     S^{X,\ast}_W : wetting phase saturation at maximum complete
                    nonwetting phase residual

     S^'_W = min_{a \in [0,t]} S_W

     P_C = P_DS^{\ast}_e^{-1/lambda}, P_N >= P_D

     Note that in general, Gerhard and Kueper include a different
     residual saturation, S_r, for the capillary pressure relationship
 
        S_e = (S_W - S_r)/(S^{\ast}_W - S_r)

     However, I will use S^k_r :=  S_r




  **********************************************************************/
public:
  GerhardKueper2p();
  virtual ~GerhardKueper2p();
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
  friend class Psk2pSpline<GerhardKueper2p>;

  //try to incorporate residual history
  virtual void setThetaS (Vec &thetaSIn);
  //set initial wetting phase volume fraction
  virtual void setInitialConditions(const Vec& local_thw);
  //set historical minimum volume fraction
  virtual void updateHistory();

  //mwf for debugging
  //private:
public:
  Vec global_psiD,global_lambda;

  Vec psiD,lambda;
  //keep track of different effective saturation scaling for wetting phase
  real sBar_krW;
  Vec global_thetaWmaxResidual;
  Vec thetaWmaxResidual;
  //mwf add for residual calculation
  //    leave thetaS as maximum available saturation corresponding to S_W=1
  Vec thetaSeff,thetaSReff;
  //mwf saturation history
  Vec thetaWminSoFar, thetaWsave; //start out at thetaS
};

inline void GerhardKueper2p::setHeads(int node,real psiWIn, real psiNIn)
{
  psiC = psiNIn - psiWIn;
  i = node;
  
  sBar = pow(psiC/psiD[i],-lambda[i]);
  thetaW = thetaSReff[i]*sBar + thetaR[i];
  DsBar_DpC = -(lambda[i]/psiD[i])*pow(psiC/psiD[i],-lambda[i]-1);
  DthetaW_DpC = thetaSReff[i] * DsBar_DpC;
  //scale relative permeability for wetting phase differently 
  sBar_krW = (thetaW-thetaR[i])/thetaSR[i];
  sBar_krW = std::max(0.0,std::min(sBar_krW,1.0));
  krW = pow(sBar_krW,(2.+3.*lambda[i])/lambda[i]);

  KW = KWs[i]*krW;
  krN = (1.-sBar)*(1-sBar)*(1.-pow(sBar,(2.0 + lambda[i])/lambda[i]));
  KN = KWs[i]*muW_by_muN*krN;
  
  if(psiC < psiD[i]) 
    {
      sBar = 1.0;
      sBar_krW = 1.0;
      thetaW = thetaSeff[i];
      DsBar_DpC = 0.0;
      DthetaW_DpC = 0.0;
      krW = 1.0;
      KW = KWs[i];
      krN = 0.0;
      KN = 0.0;
    }   
  //keep most recent thetaW?
  thetaWsave[i] = thetaW;
}

inline void GerhardKueper2p::setVFraction(const real& thetaWIn, int node)
{
  i = node;
  sBar = (thetaWIn - thetaR[i])/thetaSReff[i];
  sBar= std::min(std::max(sBar,0.0),1.0);
  thetaW = thetaSReff[i]*sBar + thetaR[i];
  psiC = psiD[i]*pow(sBar,-1./lambda[i]);
  //std::cout<<psiC<<'\t'<<psiD[i]<<'\t'<<lambda[i]<<'\t'<<sBar<<'\t'<<i<<std::endl;
  DsBar_DpC = -(lambda[i]/psiD[i])*pow(psiC/psiD[i],-lambda[i]-1.);
  DthetaW_DpC = thetaSReff[i] * DsBar_DpC;
  //use different scaling for wetting phase relative perm
  sBar_krW = (thetaWIn - thetaR[i])/thetaSR[i];
  sBar_krW = std::max(0.0,std::min(sBar_krW,1.0));
  //what about DsBar_DpC?
  krW = pow(sBar_krW,(2.+3.*lambda[i])/lambda[i]);
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
  //keep most recent thetaW?
  thetaWsave[i] = thetaW;
}
  
inline void GerhardKueper2p::calculateDerivatives()
{
  //probably not right
  //mwf add some decimals
  DDsBar_DDpC = -(lambda[i]/psiD[i])*((-lambda[i]-1.)/psiD[i])*pow(psiC/psiD[i],-lambda[i]-2.);
  DDthetaW_DDpC = thetaSReff[i]*DDsBar_DDpC;
  //mwf use different scaling for wetting phase relative perm
  real dthetaSREffdthetaSR = thetaSReff[i]/thetaSR[i];
  //real dthetaSREffdthetaSR = 1.0;
  
  DkrW_DpC = ((2.+3.*lambda[i])/lambda[i])*pow(sBar_krW,(2.+3.*lambda[i])/lambda[i]-1.)*
    DsBar_DpC*dthetaSREffdthetaSR;//include scaling to get this to be correct

  DKW_DpC = KWs[i]*DkrW_DpC;

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
  //mwf debug
  //std::cout<<"BC calcDerivs psiC= "<<psiC<<" psiD["<<i<<"]= "<<psiD[i]
  //   <<" DkrN_DpC= "<<DkrN_DpC<<std::endl;
}
inline void GerhardKueper2p::calculateDerivativesVFraction()
{
  calculateDerivatives();
}

}//Daetk
#endif






