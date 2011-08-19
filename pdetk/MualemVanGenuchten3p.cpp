#include "MualemVanGenuchten3p.h"

namespace Daetk 
{
  
MualemVanGenuchten3p::MualemVanGenuchten3p():
  Psk3p(),
  global_n(),
  global_m(),
  global_alpha(),
  n(),
  m(),
  alpha()
{}
  
MualemVanGenuchten3p::~MualemVanGenuchten3p(){}
  
void MualemVanGenuchten3p::readParameters(ParameterDatabase& pd)
{
  Psk3p::readParameters(pd);
  
  global_n.newsize(nNodes);
  global_n=pd.r("n");
  
  global_m.newsize(nNodes);
  global_m = 1.0-1.0/pd.r("n");
  
  global_alpha.newsize(nNodes);
  global_alpha = pd.r("alpha");
  
  n.newsize(Vec::LOCAL,global_n.getDA());
  n=pd.r("n");
  m.newsize(Vec::LOCAL,global_n.getDA());
  m=1.0-1.0/pd.r("n");
  alpha.newsize(Vec::LOCAL,global_n.getDA());
  alpha = pd.r("alpha");
  
}

MualemVanGenuchten3pSpline::~MualemVanGenuchten3pSpline()
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
MualemVanGenuchten3pSpline::MualemVanGenuchten3pSpline():
  ln_of_10(log(10.0))
{}

void MualemVanGenuchten3pSpline::readParameters(ParameterDatabase& pd)
{
  MualemVanGenuchten3p::readParameters(pd);

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

      MualemVanGenuchten3p::setHeads(0,0.0,psiC_spline[i]);
      MualemVanGenuchten3p::calculateDerivatives();

      thetaW_spline[i] = thetaW;
      dthetaW_spline[i] = DthetaW_DpC;

      KW_spline[i] = krW;
      dKW_spline[i] = DkrW_DpC;

      KN_spline[i] = krN;
      dKN_spline[i] = DkrN_DpC;
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

      MualemVanGenuchten3p::setVFraction(thetaW_spline_t[i],0);
      MualemVanGenuchten3p::calculateDerivativesVFraction();

//        //log(psiC) spline     
//        psiC_spline_t[i] = log10(psiC);
//        if (fabs(ln_of_10*psiC*DthetaW_DpC) > SQRT_MACHINE_EPSILON)
//          dpsiC_spline_t[i] = 1.0/(ln_of_10*psiC*DthetaW_DpC);
//        else
//          dpsiC_spline_t[i] = 0.0;

      real DpC_DthetaW;

      if (fabs(DthetaW_DpC) > 0.0)
        DpC_DthetaW = 1.0/DthetaW_DpC;
      else
        DpC_DthetaW = 0.0;
      
      //psiC spline
      psiC_spline_t[i] = psiC;
      dpsiC_spline_t[i] = DpC_DthetaW;
          
      KW_spline_t[i] = krW;
      dKW_spline_t[i] = DkrW_DpC*DpC_DthetaW;
      
      KN_spline_t[i] = krN;
      dKN_spline_t[i] = DkrN_DpC*DpC_DthetaW;
    }

  //test
//    for (int i=0;i<3*(nSplineKnots-1)+1;i++)
//      {
//        MualemVanGenuchten3p::setVFraction(thetaW_spline_t[i/3] + i%3*(splineDx_t/3.0),0);
//        MualemVanGenuchten3p::calculateDerivativesVFraction();
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

void MualemVanGenuchten3p::setKWs (Vec &KWsIn)
{
  global_KWs=KWsIn;
  KWs.startSetFromGlobal(global_KWs);
  KWs.endSetFromGlobal(global_KWs);
}

void MualemVanGenuchten3p::setAlpha (Vec &alphaIn)
{
  global_alpha=alphaIn;
  alpha.startSetFromGlobal(global_alpha);
  alpha.endSetFromGlobal(global_alpha);
}
void MualemVanGenuchten3p::setN (Vec &nIn)
{
  global_n=nIn;
  n.startSetFromGlobal(global_n);
  n.endSetFromGlobal(global_n);
  for (int i=0;i<global_n.ldim();i++)
    global_m[i]=1.0 - 1.0/global_n[i];
  m.startSetFromGlobal(global_m);
  m.endSetFromGlobal(global_m);
}
void MualemVanGenuchten3p::setThetaS (Vec &thetaSIn)
{
  global_thetaS=thetaSIn; 
  global_thetaSR = global_thetaS - global_thetaR;
  thetaS.startSetFromGlobal(global_thetaS);
  thetaS.endSetFromGlobal(global_thetaS);
  thetaSR.startSetFromGlobal(global_thetaSR);
  thetaSR.endSetFromGlobal(global_thetaSR);
}
void MualemVanGenuchten3p::setThetaR (Vec &thetaRIn)
{
  global_thetaR=thetaRIn;
  global_thetaSR = global_thetaS - global_thetaR;
  thetaR.startSetFromGlobal(global_thetaR);
  thetaR.endSetFromGlobal(global_thetaR);
  thetaSR.startSetFromGlobal(global_thetaSR);
  thetaSR.endSetFromGlobal(global_thetaSR);
}

void MualemVanGenuchten3p::millerSimilarScaling(Vec& delta)
{
  assert(global_KWs.checkConformance(delta));
  for (int i=0;i<global_KWs.getLocalHigh();i++)
    {
      global_KWs[i]*=(delta[i]*delta[i]);
      global_alpha[i]*=delta[i];
    }
  KWs.startSetFromGlobal(global_KWs);
  KWs.endSetFromGlobal(global_KWs);
  alpha.startSetFromGlobal(global_alpha);
  alpha.endSetFromGlobal(global_alpha);
}

}//Daetk
