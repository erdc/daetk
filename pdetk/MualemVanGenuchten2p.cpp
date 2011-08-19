#include "MualemVanGenuchten2p.h"

namespace Daetk 
{

  void MualemVanGenuchten2p::readZones(ParameterDatabase& pd,Petsc::SecondOrderFd& node)
  {
    int nZones = pd.i("nZones");
    real dx=0.0,dy=0.0,dz=0.0;
    dx=(pd.r("xRight") - pd.r("xLeft")) / (real(pd.i("nxNodes"))-1);
    if (pd.i("nyNodes") > 1)
      dy=(pd.r("yBack") - pd.r("yFront")) / (real(pd.i("nyNodes"))-1);
    if (pd.i("nzNodes") > 1)
      dz=(pd.r("zTop") - pd.r("zBottom")) / (real(pd.i("nzNodes"))-1);
    
    for (int i=node.local_z0;i<node.local_z0+node.local_nzNodes;i++)
      for (int j=node.local_y0;j<node.local_y0+node.local_nyNodes;j++)
        for (int k=node.local_x0;k<node.local_x0+node.local_nxNodes;k++)
          {
            node(i,j,k);
            real x=k*dx,y=j*dy,z=i*dz;
            for (int n=0;n<nZones;n++)
              {
                if (x >= pd.v("zoneLeft")(n) && x < pd.v("zoneRight")(n) &&
                    y >= pd.v("zoneFront")(n) && y < pd.v("zoneBack")(n) &&
                    z >= pd.v("zoneBottom")(n) && z < pd.v("zoneTop")(n) )
                  {
                    global_alpha(node.center) = pd.v("alphaZone")(n);
                    global_n(node.center) = pd.v("nZone")(n);
                    global_KWs(node.center) = pd.v("kZone")(n);
                    global_thetaS(node.center) = pd.v("thetaS_Zone")(n);
                    global_thetaR(node.center) = pd.v("thetaR_Zone")(n);
                  }
              }
          } 
//     global_thetaSR = global_thetaS - global_thetaR;
    global_thetaSR = global_thetaS;
    axpy(-1.0,global_thetaR,global_thetaSR);
//     std::cout<<global_psiD<<std::endl
//              <<global_lambda<<std::endl
//              <<global_thetaS<<std::endl
//              <<global_KWs<<std::endl;

     
    for (int i=0;i<global_n.ldim();i++)
      global_m[i] = 1.0-1.0/global_n[i];
    alpha.startSetFromGlobal(global_alpha);
    alpha.endSetFromGlobal(global_alpha);
    n.startSetFromGlobal(global_n);
    n.endSetFromGlobal(global_n);
    m.startSetFromGlobal(global_m);
    m.endSetFromGlobal(global_m);
    KWs.startSetFromGlobal(global_KWs);
    KWs.endSetFromGlobal(global_KWs);
    thetaS.startSetFromGlobal(global_thetaS);
    thetaS.endSetFromGlobal(global_thetaS);
    thetaR.startSetFromGlobal(global_thetaR);
    thetaR.endSetFromGlobal(global_thetaR);
    thetaSR.startSetFromGlobal(global_thetaSR);
    thetaSR.endSetFromGlobal(global_thetaSR);
  } 
MualemVanGenuchten2p::MualemVanGenuchten2p():
  Psk2p(),
  global_n(),
  global_m(),
  global_alpha(),
  n(),
  m(),
  alpha()
{}
  
MualemVanGenuchten2p::~MualemVanGenuchten2p(){}
  
void MualemVanGenuchten2p::readParameters(ParameterDatabase& pd)
{
  Psk2p::readParameters(pd);
  ns_del = pd.r("vg_ns_del");

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

MualemVanGenuchten2pSplineOld::~MualemVanGenuchten2pSplineOld()
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
MualemVanGenuchten2pSplineOld::MualemVanGenuchten2pSplineOld():
  ln_of_10(log(10.0))
{}

void MualemVanGenuchten2pSplineOld::readParameters(ParameterDatabase& pd)
{
  MualemVanGenuchten2p::readParameters(pd);

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

      MualemVanGenuchten2p::setHeads(0,0.0,psiC_spline[i]);
      MualemVanGenuchten2p::calculateDerivatives();

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

      MualemVanGenuchten2p::setVFraction(thetaW_spline_t[i],0);
      MualemVanGenuchten2p::calculateDerivativesVFraction();

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
//        MualemVanGenuchten2p::setVFraction(thetaW_spline_t[i/3] + i%3*(splineDx_t/3.0),0);
//        MualemVanGenuchten2p::calculateDerivativesVFraction();
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

void MualemVanGenuchten2p::setAlpha (Vec &alphaIn)
{
  global_alpha=alphaIn;
  alpha.startSetFromGlobal(global_alpha);
  alpha.endSetFromGlobal(global_alpha);
}
void MualemVanGenuchten2p::setN (Vec &nIn)
{
  global_n=nIn;
  n.startSetFromGlobal(global_n);
  n.endSetFromGlobal(global_n);
  for (int i=0;i<global_n.ldim();i++)
    global_m[i]=1.0 - 1.0/global_n[i];
  m.startSetFromGlobal(global_m);
  m.endSetFromGlobal(global_m);
}

void MualemVanGenuchten2p::millerSimilarScaling(Vec& delta)
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

void MualemVanGenuchten2p::setHeads(const Vec& psiW_vec, const Vec& psiN_vec)
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

void MualemVanGenuchten2p::setVFraction(const Vec& thetaW_vec)
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

void MualemVanGenuchten2p::calculateDerivativesHead(const Vec& psiW_vec, const Vec& psiN_vec)
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

void MualemVanGenuchten2p::calculateDerivativesVFraction(const Vec& thetaW_vec)
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
