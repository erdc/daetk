#include "BrooksCorey2p.h"

namespace Daetk
{

  BrooksCorey2p::BrooksCorey2p():
    Psk2p(),
    global_psiD(),
    global_lambda(),
    psiD(),
    lambda()
  {
    Tracer tr("BrooksCorey2p()");
  }
  BrooksCorey2p::~BrooksCorey2p(){}
  void BrooksCorey2p::readParameters(ParameterDatabase& pd)
  {
    Psk2p::readParameters(pd);

    global_psiD.newsize(nNodes);
    global_psiD = pd.r("psiD");

    global_lambda.newsize(nNodes);
    global_lambda = pd.r("lambda");

    psiD.newsize(Vec::LOCAL,global_psiD.getDA());
    psiD = pd.r("psiD");
    lambda.newsize(Vec::LOCAL,global_psiD.getDA());
    lambda = pd.r("lambda");

//      psiD.startSetFromGlobal(global_psiD);
//      psiD.endSetFromGlobal(global_psiD);

//      lambda.startSetFromGlobal(global_lambda);
//      lambda.endSetFromGlobal(global_lambda);
  }

  void BrooksCorey2p::readZones(ParameterDatabase& pd,Petsc::SecondOrderFd& node)
  {
    int nZones = pd.i("nZones");
    real dx=0.0,dy=0.0,dz=0.0;
    dx=(pd.r("xRight") - pd.r("xLeft")) / (real(pd.i("nxNodes"))-1);
    if (pd.i("nyNodes") > 1)
      dy=(pd.r("yBack") - pd.r("yFront")) / (real(pd.i("nyNodes"))-1);
    if (pd.i("nzNodes") > 1)
      dz=(pd.r("zTop") - pd.r("zBottom")) / (real(pd.i("nzNodes"))-1);
    real gravity=pd.r("gravity"), // m/d^2
      density=pd.r("rhoW"), // kg / m^3
      viscosity=pd.r("muW"); // kg /m d
    
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
                    global_psiD(node.center) = pd.v("pdZone")(n);
                    global_lambda(node.center) = pd.v("lambdaZone")(n);
                    global_KWs(node.center) = pd.v("permZone")(n)*gravity*density/viscosity;
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

      
    psiD.startSetFromGlobal(global_psiD);
    psiD.endSetFromGlobal(global_psiD);
    lambda.startSetFromGlobal(global_lambda);
    lambda.endSetFromGlobal(global_lambda);
    KWs.startSetFromGlobal(global_KWs);
    KWs.endSetFromGlobal(global_KWs);
    thetaS.startSetFromGlobal(global_thetaS);
    thetaS.endSetFromGlobal(global_thetaS);
    thetaR.startSetFromGlobal(global_thetaR);
    thetaR.endSetFromGlobal(global_thetaR);
    thetaSR.startSetFromGlobal(global_thetaSR);
    thetaSR.endSetFromGlobal(global_thetaSR);
  }
  
  void BrooksCorey2p::millerSimilarScaling(Vec& delta)
  {
    global_KWs.checkConformance(delta);
    for (int i=0;i<global_KWs.getLocalHigh();i++)
      {
        global_KWs[i]*=(delta[i]*delta[i]);
        global_psiD[i]/=delta[i];
      }
    KWs.startSetFromGlobal(global_KWs);
    KWs.endSetFromGlobal(global_KWs);
    psiD.startSetFromGlobal(global_psiD);
    psiD.endSetFromGlobal(global_psiD);
  }

  void BrooksCorey2p::setHeads(const Vec& psiW_vec, const Vec& psiN_vec)
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

  void BrooksCorey2p::setVFraction(const Vec& thetaW_vec)
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
  
  void BrooksCorey2p::calculateDerivativesHead(const Vec& psiW_vec, const Vec& psiN_vec)
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
  
  void BrooksCorey2p::calculateDerivativesVFraction(const Vec& thetaW_vec)
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
