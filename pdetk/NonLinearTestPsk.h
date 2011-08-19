#ifndef NONLINEAR_TEST_PSK
#define NONLINEAR_TEST_PSK

#include "Psk2p.h"

namespace Daetk
{
using std::max;
using std::min;

class NonLinearTestPsk : public Psk2p
{
public:
  NonLinearTestPsk(ParameterDatabase& pd):
    Psk2p(pd),nonlin(pd.r("anaNonLin"))
  { 

  }
  
  NonLinearTestPsk():
    Psk2p(),nonlin(1.0)
    {}
  ~NonLinearTestPsk() {}

  virtual void readParameters(ParameterDatabase& pd)
    {
      Psk2p::readParameters(pd);
      nonlin = pd.r("anaNonLin");
    }

  void setAlpha(Vec & tempA) {}
  void setN(Vec & tempN) {}
  void millerSimilarScaling(Vec& delta) {}

  void setHeads(int node, const real& psiWIn, const real& zc=0.0)
  {
    i = node;
    psiC = 0.0;
    sBar = psiWIn;
    thetaW = sBar;
    DsBar_DpC = -1.0;
    DthetaW_DpC = -1.0;
    krW = exp(-nonlin*psiWIn);
    KW = KWs[i]*krW;
    krN = 1.0;
    KN = krN;
  }
  
  void setVFraction(const real& thetaWIn, int node)
  {
//      psiC = 2*(thetaS(node) - thetaWIn);
    i = node;
    psiC = 0.0;
    sBar  = thetaWIn;
    thetaW= sBar;
//      DsBar_DpC = -0.5;
//      DthetaW_DpC = -0.5;
    DsBar_DpC = -1.0;
    DthetaW_DpC = -1.0;
    krW = exp(-nonlin*thetaWIn);
    KW = KWs[i]*krW;
    krN = 1.0;
    KN = krN;  
  }  
  
  void calculateDerivatives()
  {
      DDsBar_DDpC = 0.0;
      DDthetaW_DDpC = 0.0;
      //these get used as negative of value (actual value is -nonlin*krW)
      DkrW_DpC = nonlin*krW;
      DKW_DpC  = nonlin*KW;
      DkrN_DpC = 0.0;
      DKN_DpC = 0.0;
    }

  //
  real nonlin;

};

class NonLinearTestPsk2 : public NonLinearTestPsk
{
public:
  NonLinearTestPsk2(ParameterDatabase pd):
    NonLinearTestPsk(pd),krWoffset(1.0)
  { 
    krWoffset = pd.r("relPermOffset");
  }
  
  NonLinearTestPsk2():
    NonLinearTestPsk(), krWoffset(1.0)
    {}
  ~NonLinearTestPsk2() {}

  virtual void readParameters(ParameterDatabase& pd)
    {
      NonLinearTestPsk::readParameters(pd);
      krWoffset = pd.r("relPermOffset");
    }
  
  void setAlpha(Vec & tempA) {}
  void setN(Vec & tempN) {}
  void millerSimilarScaling(Vec& delta) {}

  void setHeads(int node, const real& psiWIn, const real& zc=0.0)
  {
    i = node;
    //cerr<<"in NLPSK 2 zc = "<<zc<<endl;
    psiC = 0.0;
    sBar = psiWIn;
    thetaW = sBar;
    DsBar_DpC = -1.0;
    DthetaW_DpC = -1.0;
    krW = nonlin*(psiWIn + zc) + krWoffset;
    assert(krW > 0.0);
    KW =  KWs[i]*krW;
    krN = 1.0;
    KN = krN;
  }
  
  void setVFraction(const real& thetaWIn, int node, const real& zc=0.0)
  {
    i = node;
//      psiC = 2*(thetaS(node) - thetaWIn);
    psiC = 0.0;
    sBar = thetaWIn;
//      DsBar_DpC = -0.5;
//      DthetaW_DpC = -0.5;
    thetaW= sBar;
    DsBar_DpC = -1.0;
    DthetaW_DpC = -1.0;
    krW = nonlin*(thetaWIn + zc) + krWoffset;
    assert(krW > 0.0);
    KW = KWs[i]*krW;
    krN = 1.0;
    KN = krN;  
  }  
  
  void calculateDerivatives()
  {
      DDsBar_DDpC = 0.0;
      DDthetaW_DDpC = 0.0;
      //this has to be negative of derivative
      DkrW_DpC = -1.0*nonlin;
      //hopefully i's been set already
      DKW_DpC = DkrW_DpC*KWs[i];
      DkrN_DpC = 0.0;
      DKN_DpC = 0.0;
    }

  //
  real krWoffset;

};

//mwf try and do Putti's second example psk relation
class PuttiPsk2p : public Psk2p
{
public:
  PuttiPsk2p():
    Psk2p(),
    global_n(),
    global_m(),
    global_alpha(),
    global_krPower(),
    n(),
    m(),
    alpha(),
    krPower(),
    psiWS(0.0)
    {}
  //mwf just added to use with same interface as
  //mwf the simple test psk for analytical sols
  //mwf ask Chris if this is alright
  PuttiPsk2p(ParameterDatabase pd):
    Psk2p(pd),
    global_n(),
    global_m(),
    global_alpha(),
    global_krPower(),
    n(),
    m(),
    alpha(),
    krPower(),
    psiWS(0.0)
   {
      readParameters(pd);
    }

  virtual ~PuttiPsk2p(){}

  void readParameters(ParameterDatabase& pd)
  {
    Psk2p::readParameters(pd);

    global_n.newsize(nNodes);
    global_n=pd.r("n");

    global_m.newsize(nNodes);
    global_m = pd.r("mPutti");

    global_alpha.newsize(nNodes);
    global_alpha = pd.r("alpha");

    global_krPower.newsize(nNodes);
    global_krPower = pd.r("krPower");

    n.newsize(Vec::LOCAL,global_n.getDA());
    n=pd.r("n");
    m.newsize(Vec::LOCAL,global_n.getDA());
    m=pd.r("mPutti");
    alpha.newsize(Vec::LOCAL,global_n.getDA());
    alpha = pd.r("alpha");

    krPower.newsize(Vec::LOCAL,global_n.getDA());
    krPower=pd.r("krPower");

    psiWS = pd.r("puttiPsiWS");
//      n.startSetFromGlobal(global_n);
//      n.endSetFromGlobal(global_n);
//      m.startSetFromGlobal(global_m);
//      m.endSetFromGlobal(global_m);
//      alpha.startSetFromGlobal(global_alpha);
//      alpha.endSetFromGlobal(global_alpha);
  }

  void setHeads(int node,real psiWIn, real psiNIn=0)
  {
    //mwf now put psiC = psiNIn - psiWIn;
    if (fabs(psiNIn) > 1.0e-16)
      psiC = psiNIn - psiWIn;
    else
      psiC = psiWS - psiWIn;

    i = node; //local node#

    if(psiC>0.0) 
      {
	//mwf same calculation for theta as mualem van genuchten class
        alphaPsiC = alpha[i]*psiC;
        alphaPsiC_n = pow(alphaPsiC,n[i]);
        alphaPsiC_nM1 = alphaPsiC_n/alphaPsiC;
        onePlus_alphaPsiC_n = 1.0 + alphaPsiC_n;
        sBar = pow(onePlus_alphaPsiC_n,-m[i]);
	sBar = min(sBar,1.0);
        sBarByOnePlus_alphaPsiC_n = sBar/onePlus_alphaPsiC_n;
        sqrt_sBar = sqrt(sBar);
        sqrt_1minusSbar = sqrt(1.0 - sBar);
        thetaW = thetaSR[i]*sBar + thetaR[i];
	//mwf change this since m and n no longer related by
	//mwf van genuchten mualem relation
        DsBar_DpC = -m[i]*alpha[i]*n[i]*alphaPsiC_nM1 
          *sBarByOnePlus_alphaPsiC_n;
        DthetaW_DpC = thetaSR[i] * DsBar_DpC; 
        vBar = 1.0-alphaPsiC_nM1*sBar;
        uBar = alphaPsiC_nM1*sBar;
	//mwf now do different krw
        krW = pow(sBar,krPower[i]);
        KW = KWs[i]*krW;
        krN = sqrt_1minusSbar*uBar*uBar;
        KN = muW_by_muN*KWs[i]*krN;
      }
    else
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
  
  void setVFraction(const real& thetaWIn, int node)
  {
    i = node;
    sBar = (thetaWIn - thetaR[i])/thetaSR[i];
//      std::cout<<"i "<<i<<" thetaWIn "<<thetaWIn<<" thetaSR "<<thetaSR[i]<<" thetaR "<<thetaR[i]<<std::endl<<std::flush;
//      std::cout<<"sBar"<<sBar<<std::endl<<std::flush;
    if (sBar > 1.0)
      {
        sBar = 1.0;
        thetaW = thetaS[i];
      }
    else if (sBar <= 0.0)
      {
        std::cerr<<"thetaW has been reduced below resdiual; exiting"<<std::endl;
        exit(1);
      }
    else
      thetaW = thetaWIn;
    onePlus_alphaPsiC_n = pow(sBar,1.0/-m[i]);
    alphaPsiC_n = onePlus_alphaPsiC_n - 1.0;
    alphaPsiC = pow(alphaPsiC_n,1.0/n[i]);
    psiC = alphaPsiC/alpha[i];

//      std::cout<<"psiC"<<psiC<<std::endl<<std::flush;
    if(psiC>0.0) 
      {
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
        //krW = sqrt_sBar*vBar*vBar;
	//mwf now do different krw
        krW = pow(sBar,krPower[i]);
        KW = KWs[i]*krW;
        krN = sqrt_1minusSbar*uBar*uBar;
        KN = muW_by_muN*KWs[i]*krN;
      }
    else
      {
        DsBar_DpC = 0.0;
        DthetaW_DpC = 0.0;
        krW = 1.0;
        KW = KWs[i];
        krN = 0.0;
        KN = 0.0;
      }
  }
  
 
  void calculateDerivatives()
    {
      if (psiC > 0.0)
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
	  //mwf now do different calc for putti problem
//            DkrW_DpC = (0.5/sqrt_sBar)*DsBar_DpC*vBar*vBar
//              -
//              2.0*sqrt_sBar*vBar*
//              (alpha[i]*(n[i]-1.0)*alphaPsiC_nM2*sBar
//               + alphaPsiC_nM1 * DsBar_DpC);
	  //
	  DkrW_DpC = krPower[i]*pow(sBar,krPower[i]-1)*DsBar_DpC;
          DKW_DpC = KWs[i]*DkrW_DpC;
          if (sqrt_1minusSbar > SQRT_MACHINE_EPSILON)
            {
              DkrN_DpC = -(0.5/sqrt_1minusSbar)*DsBar_DpC*uBar*uBar
                +
                2.0*sqrt_1minusSbar*uBar*
                (alpha[i]*(n[i]-1.0)*alphaPsiC_nM2*sBar
                 + alphaPsiC_nM1 * DsBar_DpC);
            }
          else
            {
              DkrN_DpC =((1.0 - sBar)/MACHINE_EPSILON)*2.0*SQRT_MACHINE_EPSILON*uBar*
                (alpha[i]*(n[i]-1.0)*alphaPsiC_nM2*sBar
                 + alphaPsiC_nM1 * DsBar_DpC)
                - (DsBar_DpC/MACHINE_EPSILON)*SQRT_MACHINE_EPSILON*uBar*uBar;
            }
          DKN_DpC = muW_by_muN*KWs[i]*DkrN_DpC;
        }
      else
        {
          DDsBar_DDpC = 0.0;
          DDthetaW_DDpC = 0.0;
          DkrW_DpC = 0.0;
          DKW_DpC = 0.0;
          DkrN_DpC = 0.0;
          DKN_DpC = 0.0;
        }
    }

  //mwf add to set parameters for heterogeneous test problem
  void setAlpha(const Vec& alphaIn)
    {
      global_alpha = alphaIn;
      alpha.newsize(Vec::LOCAL,global_alpha.getDA());
      alpha.startSetFromGlobal(global_alpha);
      alpha.endSetFromGlobal(global_alpha);
    }

  void setKWs (Vec &KWsIn)
  {
    global_KWs=KWsIn;
    KWs.startSetFromGlobal(global_KWs);
    KWs.endSetFromGlobal(global_KWs);
  }
  void setAlpha (Vec &alphaIn)
  {
    global_alpha=alphaIn;
    alpha.startSetFromGlobal(global_alpha);
    alpha.endSetFromGlobal(global_alpha);
  }
  void setN (Vec &nIn)
  {
    global_n=nIn;
    n.startSetFromGlobal(global_n);
    n.endSetFromGlobal(global_n);
  }
  void setM (Vec &mIn)
  {
    global_m=mIn;
    m.startSetFromGlobal(global_m);
    m.endSetFromGlobal(global_m);
  }
  void setkrPower (Vec &krPowerIn)
  {
    global_krPower=krPowerIn;
    krPower.startSetFromGlobal(global_krPower);
    krPower.endSetFromGlobal(global_krPower);
  }
  void setThetaS (Vec &thetaSIn)
  {
    global_thetaS=thetaSIn; 
    global_thetaSR = global_thetaS - global_thetaR;
    thetaS.startSetFromGlobal(global_thetaS);
    thetaS.endSetFromGlobal(global_thetaS);
    thetaSR.startSetFromGlobal(global_thetaSR);
    thetaSR.endSetFromGlobal(global_thetaSR);
  }
  void setThetaR (Vec &thetaRIn)
  {
    global_thetaR=thetaRIn;
    global_thetaSR = global_thetaS - global_thetaR;
    thetaR.startSetFromGlobal(global_thetaR);
    thetaR.endSetFromGlobal(global_thetaR);
    thetaSR.startSetFromGlobal(global_thetaSR);
    thetaSR.endSetFromGlobal(global_thetaSR);
  }


  void millerSimilarScaling(Vec& delta)
  {
    for (int i=0;i<KWs.getLocalHigh();i++)
      {
        global_KWs[i]*=delta[i]*delta[i];
        global_alpha[i]*=delta[i];
      }
    KWs.startSetFromGlobal(global_KWs);
    KWs.endSetFromGlobal(global_KWs);
    alpha.startSetFromGlobal(global_alpha);
    alpha.endSetFromGlobal(global_alpha);
  }
protected:
  real vBar,uBar,
    alphaPsiC, alphaPsiC_n, alphaPsiC_nM1, alphaPsiC_nM2,
    onePlus_alphaPsiC_n,
    sqrt_sBar, sqrt_1minusSbar,
    sBarByOnePlus_alphaPsiC_n, sBarBy_onePlus_alphaPsiC_n_2;


  Vec global_n,global_m,global_alpha,global_krPower;
  Vec n,m,alpha,krPower;
  //mwf for 3rd example so I don't have to change code
  real psiWS;
};

}//Daetk
#endif
