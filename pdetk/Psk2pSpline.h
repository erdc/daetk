#ifndef PSK2PSPLINE_H
#define PSK2PSPLINE_H

#include "Definitions.h"
#include "ParameterDatabase.h"
#include "VecOperators.h"
#include "Psk2p.h"

namespace Daetk 
{
template<class ANALYTIC_PSK> 
class Psk2pSpline : public Psk2p
{
public:
  virtual ~Psk2pSpline();
  Psk2pSpline();
  void readParameters(ParameterDatabase& pd);
  void setHeads(int node,real psiWIn, real psiNIn=0);
  void setVFraction(const real& thetaWIn, int node);
  void calculateDerivatives();
  void calculateDerivativesVFraction();
  
  //vector call sequences
  virtual void setHeads(const Vec& psiW_vec, const Vec& psiN_vec);
  virtual void setVFraction(const Vec& thetaW_vec);
  virtual void calculateDerivativesHead(const Vec& psiW_vec, const Vec& psiN_vec);
  virtual void calculateDerivativesVFraction(const Vec& thetaW_vec);
  
  void millerSimilarScaling(Vec& delta){mypsk.millerSimilarScaling(delta);}
  void setAlpha(Vec &alphaIn){mypsk.setAlpha(alphaIn);}
  void setN (Vec &nIn){mypsk.setN(nIn);}

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

  ANALYTIC_PSK mypsk;
};

template<class ANALYTIC_PSK> 
inline void Psk2pSpline<ANALYTIC_PSK>::setHeads(int node,real psiWIn, real psiNIn)
{
  //only for homogeneous right now--we ignore node
  psiC = psiNIn - psiWIn;
  i=node;

  index = (int)(oneOverSplineDx*(psiC-splineMinPsiC));
  index = std::max(std::min(index,nSplineKnots-2),0);
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
    {
      mypsk.setHeads(node,psiWIn,psiNIn);
      thetaW = mypsk.thetaW;
      DthetaW_DpC = mypsk.DthetaW_DpC;
      krW = mypsk.krW;
      KW = mypsk.KW;
      krN = mypsk.krN;
      KN = mypsk.KN;
    }
}

template<class ANALYTIC_PSK> 
inline void Psk2pSpline<ANALYTIC_PSK>::setVFraction(const real& thetaWIn, int node)
{
  thetaW = thetaWIn;
  i=node;
  //only for homogeneous right now--we ignore node
  index = (int)(oneOverSplineDx_t*(thetaW-thetaW_spline_t[0]));
  index = std::max(std::min(index,nSplineKnots-2),0);
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
    {
      mypsk.setVFraction(thetaW,node);
      psiC = mypsk.psiC;
      DthetaW_DpC = mypsk.DthetaW_DpC;
      krW = mypsk.krW;
      KW = mypsk.KW;
      krN = mypsk.krN;
      KN = mypsk.KN;
    }
}

template<class ANALYTIC_PSK> 
inline void Psk2pSpline<ANALYTIC_PSK>::calculateDerivatives()
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
    {
      mypsk.calculateDerivatives();
      DDthetaW_DDpC = mypsk.DDthetaW_DDpC;
      DkrW_DpC = mypsk.DkrW_DpC;
      DKW_DpC = mypsk.DKW_DpC;
      DkrN_DpC = mypsk.DkrN_DpC;
      DKN_DpC = mypsk.DKN_DpC;
    }
}

template<class ANALYTIC_PSK> 
inline void Psk2pSpline<ANALYTIC_PSK>::calculateDerivativesVFraction()
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
    {
      mypsk.calculateDerivativesVFraction();
      DDthetaW_DDpC = mypsk.DDthetaW_DDpC;
      DkrW_DpC = mypsk.DkrW_DpC;
      DKW_DpC = mypsk.DKW_DpC;
      DkrN_DpC = mypsk.DkrN_DpC;
      DKN_DpC = mypsk.DKN_DpC;
    }
}

template<class ANALYTIC_PSK> 
Psk2pSpline<ANALYTIC_PSK>::~Psk2pSpline()
{
  delete []   psiC_spline;
  delete []   thetaW_spline;
  delete []   dthetaW_spline;
  delete []   KW_spline;
  delete []   dKW_spline;
  delete []   KN_spline;
  delete []   dKN_spline;
  delete []   psiC_spline_t;
  delete []   thetaW_spline_t;
  delete []   dpsiC_spline_t;
  delete []   KW_spline_t;
  delete []   dKW_spline_t;
  delete []   KN_spline_t;
  delete []   dKN_spline_t;
}

template<class ANALYTIC_PSK> 
Psk2pSpline<ANALYTIC_PSK>::Psk2pSpline():
  ln_of_10(log(10.0))
{}

template<class ANALYTIC_PSK> 
void Psk2pSpline<ANALYTIC_PSK>::readParameters(ParameterDatabase& pd)
{
  mypsk.readParameters(pd);
  psiC_spline = new real[pd.i("nSplineKnots")];
  thetaW_spline = new real[pd.i("nSplineKnots")];
  dthetaW_spline = new real[pd.i("nSplineKnots")];
  KW_spline = new real[pd.i("nSplineKnots")];
  dKW_spline = new real[pd.i("nSplineKnots")];
  KN_spline = new real[pd.i("nSplineKnots")];
  dKN_spline = new real[pd.i("nSplineKnots")];
  psiC_spline_t = new real[pd.i("nSplineKnots")];
  thetaW_spline_t = new real[pd.i("nSplineKnots")];
  dpsiC_spline_t = new real[pd.i("nSplineKnots")];
  KW_spline_t = new real[pd.i("nSplineKnots")];
  dKW_spline_t = new real[pd.i("nSplineKnots")];
  KN_spline_t = new real[pd.i("nSplineKnots")];
  dKN_spline_t = new real[pd.i("nSplineKnots")];

  VecIndex all;
  global_thetaS.attachToVec(Vec::REF,mypsk.global_thetaS,all);
  global_thetaR.attachToVec(Vec::REF,mypsk.global_thetaR,all);
  global_KWs.attachToVec(Vec::REF,mypsk.global_KWs,all);
  global_thetaSR.attachToVec(Vec::REF,mypsk.global_thetaSR,all);

  //need to do these as references too, but I don't have that capability in PetscVec yet
  thetaS.newsize(Vec::LOCAL,mypsk.global_thetaS.getDA());
  thetaR.newsize(Vec::LOCAL,mypsk.global_thetaS.getDA());
  KWs.newsize(Vec::LOCAL,mypsk.global_thetaS.getDA());
  thetaSR.newsize(Vec::LOCAL,mypsk.global_thetaS.getDA());
  thetaS = mypsk.thetaS;
  thetaR = mypsk.thetaR;
  KWs = mypsk.KWs;
  thetaSR = mypsk.thetaSR;

  //read in real and integer constants for building table

  splineMinPsiC = pd.r("splineMinPsiC");
  splineMaxPsiC = pd.r("splineMaxPsiC");
  nSplineKnots  = pd.i("nSplineKnots");
  
  splineDx = (splineMaxPsiC - splineMinPsiC) / ( (real) (nSplineKnots-1) );
  oneOverSplineDx = 1.0/splineDx;
  
  psiC_spline[0] = splineMinPsiC;
  for(int i=1; i<nSplineKnots; i++)
    psiC_spline[i] = splineMinPsiC + (real)(i)*splineDx;       
  splineMaxPsiC = psiC_spline[nSplineKnots-1];
  for(int i=0; i<nSplineKnots; i++) 
    {
      //use node 0 to set up the spline table
      //that means the problem MUST be homogeneous

      mypsk.setHeads(0,0.0,psiC_spline[i]);
      mypsk.calculateDerivatives();

      thetaW_spline[i] = mypsk.thetaW;
      dthetaW_spline[i] = mypsk.DthetaW_DpC;

      KW_spline[i] = mypsk.krW;
      dKW_spline[i] = mypsk.DkrW_DpC;

      KN_spline[i] = mypsk.krN;
      dKN_spline[i] = mypsk.DkrN_DpC;
    }
  
  splineMaxThetaW = thetaW_spline[0];
  splineMinThetaW = thetaW_spline[nSplineKnots-1];

  splineDx_t = (splineMaxThetaW - splineMinThetaW) / ( (real) (nSplineKnots-1) );
  oneOverSplineDx_t = 1.0/splineDx_t;
  
  thetaW_spline_t[0] = splineMinThetaW;
  for(int i=1; i<nSplineKnots; i++)
    thetaW_spline_t[i] =  thetaW_spline_t[0] + (real)(i)*splineDx_t;       
  splineMaxThetaW = thetaW_spline_t[nSplineKnots-1];

  for(int i=0; i<nSplineKnots; i++) 
    {
      //use node 0 to set up the spline table
      //that means the problem MUST be homogeneous

      mypsk.setVFraction(thetaW_spline_t[i],0);
      mypsk.calculateDerivativesVFraction();

//        //log(psiC) spline     
//        psiC_spline_t[i] = log10(psiC);
//        if (fabs(ln_of_10*psiC*DthetaW_DpC) > SQRT_MACHINE_EPSILON)
//          dpsiC_spline_t[i] = 1.0/(ln_of_10*psiC*DthetaW_DpC);
//        else
//          dpsiC_spline_t[i] = 0.0;

      real DpC_DthetaW;

      if (fabs(DthetaW_DpC) > 0.0)
        DpC_DthetaW = 1.0/mypsk.DthetaW_DpC;
      else
        DpC_DthetaW = 0.0;
      
      //psiC spline
      psiC_spline_t[i] = mypsk.psiC;
      dpsiC_spline_t[i] = DpC_DthetaW;
          
      KW_spline_t[i] = mypsk.krW;
      dKW_spline_t[i] = mypsk.DkrW_DpC*DpC_DthetaW;
      
      KN_spline_t[i] = mypsk.krN;
      dKN_spline_t[i] = mypsk.DkrN_DpC*DpC_DthetaW;
    }

  //test
//    for (int i=0;i<3*(nSplineKnots-1)+1;i++)
//      {
//        mypsk.setVFraction(thetaW_spline_t[i/3] + i%3*(splineDx_t/3.0),0);
//        mypsk.calculateDerivativesVFraction();
//        std::cout<<psiC<<'\t'<<DthetaW_DpC<<'\t'<<'\t'<<DDthetaW_DDpC<<'\t'<<KW<<'\t'<<DKW_DpC<<'\t'<<KN<<'\t'<<DKN_DpC<<std::endl;
//        setVFraction(thetaW_spline_t[i/3] + i%3*(splineDx_t/3.0),0);
//        calculateDerivativesVFraction();
//        std::cout<<psiC<<'\t'<<DthetaW_DpC<<'\t'<<'\t'<<DDthetaW_DDpC<<'\t'<<KW<<'\t'<<DKW_DpC<<'\t'<<KN<<'\t'<<DKN_DpC<<std::endl;
//      } 

  //std::cout<<splineMaxThetaW<<'\t'<<thetaW_spline_t[nSplineKnots-1]<<std::endl;

  //optimizations
  oneOverSplineDx_t_2 =oneOverSplineDx_t*2.0;
  oneOverSplineDx_t_1p5=oneOverSplineDx_t*1.5;
  splineDx_t_p125=splineDx_t*0.125;
}
  
template<class ANALYTIC_PSK> 
void Psk2pSpline<ANALYTIC_PSK>::setHeads(const Vec& psiW_vec, const Vec& psiN_vec)
{
  for(int i=0;i<psiW_vec.ldim();i++)
    {
      setHeads(i,psiW_vec[i],psiN_vec[i]);
      (*thetaW_p)[i] = thetaW;
      (*krW_p)[i] = krW;
      if (DthetaW_DpC_p)
        (*DthetaW_DpC_p)[i]=DthetaW_DpC;
      if (krN_p)
        (*krN_p)[i]=krN;
    }
}

template<class ANALYTIC_PSK> 
void Psk2pSpline<ANALYTIC_PSK>::setVFraction(const Vec& thetaW_vec)
{
  for(int i=0;i<thetaW_vec.ldim();i++)
    {
      setVFraction(thetaW_vec[i],i);
      (*psiC_p)[i] = psiC;
      (*krW_p)[i] = krW;
      if (DthetaW_DpC_p)
        (*DthetaW_DpC_p)[i]=DthetaW_DpC;
      if (krN_p)
        (*krN_p)[i]=krN;
    }
}

template<class ANALYTIC_PSK> 
void Psk2pSpline<ANALYTIC_PSK>::calculateDerivativesHead(const Vec& psiW_vec, const Vec& psiN_vec)
{
  for(int i=0;i<psiW_vec.ldim();i++)
    {
      setHeads(i,psiW_vec[i],psiN_vec[i]);
      calculateDerivatives();
      (*thetaW_p)[i] = thetaW;
      (*DthetaW_DpC_p)[i] = DthetaW_DpC;
      (*krW_p)[i] = krW;
      (*DkrW_DpC_p)[i] = DkrW_DpC;
      
      if (DDthetaW_DDpC_p)
        (*DDthetaW_DDpC_p)[i]=DDthetaW_DDpC;
      
      if (krN_p)
        {
          (*krN_p)[i]=krN;
          (*DkrN_DpC_p)[i]=DkrN_DpC;
        }
    }
}

template<class ANALYTIC_PSK> 
void Psk2pSpline<ANALYTIC_PSK>::calculateDerivativesVFraction(const Vec& thetaW_vec)
{
  for(int i=0;i<thetaW_vec.ldim();i++)
    {
      setVFraction(thetaW_vec[i],i);
      calculateDerivativesVFraction();
      (*psiC_p)[i] = psiC;
      (*DthetaW_DpC_p)[i] = DthetaW_DpC;
      (*krW_p)[i] = krW;
      (*DkrW_DpC_p)[i] = DkrW_DpC;
      
      if (DDthetaW_DDpC_p)
        (*DDthetaW_DDpC_p)[i]=DDthetaW_DDpC;
      
      if (krN_p)
        {
          (*krN_p)[i]=krN;
          (*DkrN_DpC_p)[i]=DkrN_DpC;
        }
    }
}
}//Daetk
#endif
