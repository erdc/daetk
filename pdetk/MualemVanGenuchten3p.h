#ifndef MUALEMVANGENUCHTEN3P_H
#define MUALEMVANGENUCHTEN3P_H

#include "Definitions.h"
#include "Psk3p.h"
#include "ParameterDatabase.h"
#include "VecOperators.h"
#include <stdlib.h>
#include <math.h>
namespace Daetk 
{
class MualemVanGenuchten3p : public Psk3p
{
public:
  MualemVanGenuchten3p();
  virtual ~MualemVanGenuchten3p();
  void readParameters(ParameterDatabase& pd);
  inline void setHeads(int node,real psiWIn, real psiNIn=0);
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
  
  typedef Psk3pSpline<MualemVanGenuchten3p> MualemVanGenuchten3pSplineTest;

class MualemVanGenuchten3pSpline : public  MualemVanGenuchten3p
{
public:
  virtual ~MualemVanGenuchten3pSpline();
  MualemVanGenuchten3pSpline();
  void readParameters(ParameterDatabase& pd);
  void setHeads(int node,real psiWIn, real psiNIn=0);
  void setVFraction(const real& thetaWIn, int node);
  void calculateDerivatives();
  void calculateDerivativesVFraction();

private:
  int nSplineKnots,index;

  real splineMinPsiC,
    splineMaxPsiC,
    splineMinThetaW,
    splineMaxThetaW,
    splineDx,
    oneOverSplineDx,
    splineDx_t,
    oneOverSplineDx_t,
    xi,xm1,xp1,xp2,xi2m1,x3p1,x3m1,
    N01,N02,N11,N12,
    dN01,dN02,dN11,dN12,
    ddN01,ddN02,ddN11,ddN12,
    ln_of_10,
    xi2,xi_3,xm12,
    oneOverSplineDx_t_2,oneOverSplineDx_t_1p5,splineDx_t_p125;

  real *psiC_spline,
    *thetaW_spline,
    *dthetaW_spline,
    *KW_spline,
    *dKW_spline,
    *KN_spline,
    *dKN_spline;

  real *psiC_spline_t,
    *thetaW_spline_t,
    *dpsiC_spline_t,
    *KW_spline_t,
    *dKW_spline_t,
    *KN_spline_t,
    *dKN_spline_t;

};

inline void MualemVanGenuchten3p::setHeads(int node,real psiWIn, real psiNIn)
{
  psiC = psiNIn - psiWIn;
  i = node; //local node#
  
  alphaPsiC = alpha[i]*psiC;
  alphaPsiC_n = pow(alphaPsiC,n[i]);
  alphaPsiC_nM1 = alphaPsiC_n/alphaPsiC;
  onePlus_alphaPsiC_n = 1.0 + alphaPsiC_n;
  sBar = pow(onePlus_alphaPsiC_n,-m[i]);
  sBarByOnePlus_alphaPsiC_n = sBar/onePlus_alphaPsiC_n;
  sqrt_sBar = sqrt(sBar);
  sqrt_1minusSbar = sqrt(1.0 - sBar);
  thetaW = thetaSR[i]*sBar + thetaR[i];
  DsBar_DpC = -alpha[i]*(n[i]-1.0)*alphaPsiC_nM1 
    *sBarByOnePlus_alphaPsiC_n;
  DthetaW_DpC = thetaSR[i] * DsBar_DpC; 
  vBar = 1.0-alphaPsiC_nM1*sBar;
  uBar = alphaPsiC_nM1*sBar;
  krW = sqrt_sBar*vBar*vBar;
  KW = KWs[i]*krW;
  krN = sqrt_1minusSbar*uBar*uBar;
  KN = KWs[i]*muW_by_muN*krN;

  if(psiC<=0.0) 
    {
      sBar = 1.0;
      thetaW = thetaS[i];
      DsBar_DpC = 0.0;
      DthetaW_DpC = 0.0;
      krW = 1.0;
      KW =  KWs[i];
      krN = 0.0;
      KN = 0.0;
    }    
}

inline void MualemVanGenuchten3p::setVFraction(const real& thetaWIn, int node)
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
  sqrt_sBar = sqrt(sBar);
  sqrt_1minusSbar = sqrt(1.0 - sBar);
  thetaW = thetaSR[i]*sBar + thetaR[i];
  DsBar_DpC = -alpha[i]*(n[i]-1.0)*alphaPsiC_nM1 
    *sBarByOnePlus_alphaPsiC_n;
  DthetaW_DpC = thetaSR[i] * DsBar_DpC; 
  vBar = 1.0-alphaPsiC_nM1*sBar;
  uBar = alphaPsiC_nM1*sBar;
  krW = sqrt_sBar*vBar*vBar;
  KW = KWs[i]*krW;
  krN = sqrt_1minusSbar*uBar*uBar;
  KN = KWs[i]*muW_by_muN*krN;
  if(psiC<=0.0) 
    {
      DsBar_DpC = 0.0;
      DthetaW_DpC = 0.0;
      krW = 1.0;
      KW = KWs[i];
      krN = 0.0;
      KN = 0.0;
    }
}


inline void MualemVanGenuchten3p::calculateDerivatives()
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
  DDthetaW_DDpC = thetaSR[i]*DDsBar_DDpC;

  DkrW_DpC = (0.5/sqrt_sBar)*DsBar_DpC*vBar*vBar
    -
    2.0*sqrt_sBar*vBar*
    (alpha[i]*(n[i]-1.0)*alphaPsiC_nM2*sBar
     + alphaPsiC_nM1 * DsBar_DpC);

  DKW_DpC = KWs[i]*DkrW_DpC;

  DkrN_DpC = -(0.5/sqrt_1minusSbar)*DsBar_DpC*uBar*uBar
    +
    2.0*sqrt_1minusSbar*uBar*
    (alpha[i]*(n[i]-1.0)*alphaPsiC_nM2*sBar
     + alphaPsiC_nM1 * DsBar_DpC);
  
  //recalculate if necessary
  if (sqrt_1minusSbar <= SQRT_MACHINE_EPSILON)
    {
      DkrN_DpC =((1.0 - sBar)/MACHINE_EPSILON)*2.0*SQRT_MACHINE_EPSILON*uBar*
        (alpha[i]*(n[i]-1.0)*alphaPsiC_nM2*sBar
         + alphaPsiC_nM1 * DsBar_DpC)
        - (DsBar_DpC/MACHINE_EPSILON)*SQRT_MACHINE_EPSILON*uBar*uBar;
    }

  DKN_DpC = KWs[i]*muW_by_muN*DkrN_DpC;

  if (psiC <= 0.0)
    {
      DDsBar_DDpC = 0.0;
      DDthetaW_DDpC = 0.0;
      DkrW_DpC = 0.0;
      DKW_DpC = 0.0;
      DkrN_DpC = 0.0;
      DKN_DpC = 0.0;
    }
}

inline void MualemVanGenuchten3pSpline::setHeads(int node,real psiWIn, real psiNIn)
{
  //only for homogeneous right now--we ignore node
  psiC = psiNIn - psiWIn;
  
  index = (int)(oneOverSplineDx*(psiC-splineMinPsiC));
  index = std::max(std::min(index,nSplineKnots-1),0);
  xi = 2.0*oneOverSplineDx*(psiC-psiC_spline[index]) - 1.0;
  xm1 = xi-1.0;
  xp1 = xi+1.0;
  xp2 = xi+2.0;
  xi2m1 = xi*xi-1.0;
  x3p1 = 3.0*xi+1.0;
  x3m1 = 3.0*xi-1.0;
  
  N01 = 0.25*xm1*xm1*xp2;
  N02 = 1.0-N01;
  N11 = 0.125*splineDx*xm1*xm1*xp1;
  N12 = N11*xp1/xm1;
  
  dN01 = 1.5*oneOverSplineDx*xi2m1;
  dN02 = -dN01;
  dN11 = 0.25*x3p1*xm1;
  dN12 = 0.25*x3m1*xp1;
  
  thetaW = N01*thetaW_spline[index] + N02*thetaW_spline[index+1] +
    N11*dthetaW_spline[index] + N12*dthetaW_spline[index+1];
  
  DthetaW_DpC = dN01*thetaW_spline[index] + dN02*thetaW_spline[index+1] +
    dN11*dthetaW_spline[index] + dN12*dthetaW_spline[index+1];
  
  krW = N01*KW_spline[index] + N02*KW_spline[index+1] +
    N11*dKW_spline[index] + N12*dKW_spline[index+1];     
  
  KW = KWs[i]*krW;

  krN = N01*KN_spline[index] + N02*KN_spline[index+1] +
    N11*dKN_spline[index] + N12*dKN_spline[index+1];     

  KN = KWs[i]*muW_by_muN*krN;

  if (psiC < splineMinPsiC || psiC > splineMaxPsiC)
    MualemVanGenuchten3p::setHeads(node,psiWIn,psiNIn);
}

inline void MualemVanGenuchten3pSpline::setVFraction(const real& thetaWIn, int node)
{
  thetaW = thetaWIn;
  
  //only for homogeneous right now--we ignore node
  index = (int)(oneOverSplineDx_t*(thetaW-thetaW_spline_t[0]));
  xi = 2.0*oneOverSplineDx_t*(thetaW-thetaW_spline_t[index]) - 1.0;
  //        std::cout<<index<<'\t'<<xi<<'\t'<<thetaW_spline_t[index]<<'\t'<<'\t'<<thetaW<<'\t'<<thetaW_spline_t[index+1]<<std::endl;
  xm1 = xi-1.0;
  xp1 = xi+1.0;
  xp2 = xi+2.0;
  xi2m1 = xi*xi-1.0;
  x3p1 = 3.0*xi+1.0;
  x3m1 = 3.0*xi-1.0;
  
  N01 = 0.25*xm1*xm1*xp2;
  N02 = 1.0-N01;
  N11 = 0.125*splineDx_t*xm1*xm1*xp1;
  N12 = N11*xp1/xm1;
  
  dN01 = 1.5*oneOverSplineDx_t*xi2m1;
  dN02 = -dN01;
  dN11 = 0.25*x3p1*xm1;
  dN12 = 0.25*x3m1*xp1;
  
  //        // log(psiC) spline
  //        psiC = pow(10.0,N01*psiC_spline_t[index] + N02*psiC_spline_t[index+1] +
  //                   N11*dpsiC_spline_t[index] + N12*dpsiC_spline_t[index+1]);
  //        //actuall DpC_DthetaW
  //        DthetaW_DpC = psiC*ln_of_10*(dN01*psiC_spline_t[index] + dN02*psiC_spline_t[index+1] +
  //                                     dN11*dpsiC_spline_t[index] + dN12*dpsiC_spline_t[index+1]);
  
  // psiC spline
  psiC = N01*psiC_spline_t[index] + N02*psiC_spline_t[index+1] +
    N11*dpsiC_spline_t[index] + N12*dpsiC_spline_t[index+1];
  //actually DpC_DthetaW
  DthetaW_DpC = dN01*psiC_spline_t[index] + dN02*psiC_spline_t[index+1] +
    dN11*dpsiC_spline_t[index] + dN12*dpsiC_spline_t[index+1];
  
  if (fabs(DthetaW_DpC) > 0.0)
    DthetaW_DpC = 1.0/DthetaW_DpC;
  else
    DthetaW_DpC = 0.0;
  
  krW = N01*KW_spline_t[index] + N02*KW_spline_t[index+1] +
    N11*dKW_spline_t[index] + N12*dKW_spline_t[index+1];     
  
  KW = KWs[i]*krW;
  
  krN = N01*KN_spline_t[index] + N02*KN_spline_t[index+1] +
    N11*dKN_spline_t[index] + N12*dKN_spline_t[index+1];     
  
  KN = KWs[i]*muW_by_muN*krN;
  
  if (thetaW <= splineMinThetaW || thetaW >= splineMaxThetaW)
    MualemVanGenuchten3p::setVFraction(thetaW,node);
}

inline void MualemVanGenuchten3pSpline::calculateDerivatives()
{
  ddN01 =  6.0 * oneOverSplineDx*oneOverSplineDx*xi;
  ddN02 = - ddN01;
  ddN11 = oneOverSplineDx*x3m1;
  ddN12 = oneOverSplineDx*x3p1;
  
  DDthetaW_DDpC = ddN01*thetaW_spline[index] + ddN02*thetaW_spline[index+1] +
    ddN11*dthetaW_spline[index] + ddN12*dthetaW_spline[index+1];

  DkrW_DpC = dN01*KW_spline[index] + dN02*KW_spline[index+1] +
    dN11*dKW_spline[index] + dN12*dKW_spline[index+1];

  DKW_DpC = KWs[i]*DkrW_DpC;

  DkrN_DpC = dN01*KN_spline[index] + dN02*KN_spline[index+1] +
    dN11*dKN_spline[index] + dN12*dKN_spline[index+1];

  DKN_DpC = KWs[i]*muW_by_muN*DkrN_DpC;

  if (psiC < splineMinPsiC || psiC > splineMaxPsiC)
    MualemVanGenuchten3p::calculateDerivatives();
}

inline void MualemVanGenuchten3pSpline::calculateDerivativesVFraction()
{
  ddN01 =  6.0 * oneOverSplineDx_t*oneOverSplineDx_t*xi;
  ddN02 = - ddN01;
  ddN11 = oneOverSplineDx_t*x3m1;
  ddN12 = oneOverSplineDx_t*x3p1;

  //really DDpC_DDthetaW
  DDthetaW_DDpC = ddN01*psiC_spline_t[index] + ddN02*psiC_spline_t[index+1] +
    ddN11*dpsiC_spline_t[index] + ddN12*dpsiC_spline_t[index+1];
  //now compute DDthetaW_DDpC
  DDthetaW_DDpC = - DthetaW_DpC*DthetaW_DpC*DDthetaW_DDpC*DthetaW_DpC;
  
  DkrW_DpC = (dN01*KW_spline_t[index] + dN02*KW_spline_t[index+1] +
              dN11*dKW_spline_t[index] + dN12*dKW_spline_t[index+1])*DthetaW_DpC;
  
  DKW_DpC = KWs[i]*DkrW_DpC;
  
  DkrN_DpC = (dN01*KN_spline_t[index] + dN02*KN_spline_t[index+1] +
              dN11*dKN_spline_t[index] + dN12*dKN_spline_t[index+1])*DthetaW_DpC;
  
  DKN_DpC = KWs[i]*muW_by_muN*DkrN_DpC;
  
  if (thetaW <= splineMinThetaW || thetaW >= splineMaxThetaW)
    MualemVanGenuchten3p::calculateDerivativesVFraction();
}

}//Daetk
#endif





