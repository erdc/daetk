#include "VanGenuchten2p.h"

namespace Daetk 
{
  
VanGenuchten2p::VanGenuchten2p():
  Psk2p(),
  global_n(),
  global_m(),
  global_alpha(),
  n(),
  m(),
  alpha()
{}
  
VanGenuchten2p::~VanGenuchten2p(){}
  
void VanGenuchten2p::readParameters(ParameterDatabase& pd)
{
  Psk2p::readParameters(pd);
  
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

void VanGenuchten2p::setKWs (Vec &KWsIn)
{
  global_KWs=KWsIn;
  KWs.startSetFromGlobal(global_KWs);
  KWs.endSetFromGlobal(global_KWs);
}

void VanGenuchten2p::setAlpha (Vec &alphaIn)
{
  global_alpha=alphaIn;
  alpha.startSetFromGlobal(global_alpha);
  alpha.endSetFromGlobal(global_alpha);
}
void VanGenuchten2p::setN (Vec &nIn)
{
  global_n=nIn;
  n.startSetFromGlobal(global_n);
  n.endSetFromGlobal(global_n);
  for (int i=0;i<global_n.ldim();i++)
    global_m[i]=1.0 - 1.0/global_n[i];
  m.startSetFromGlobal(global_m);
  m.endSetFromGlobal(global_m);
}
void VanGenuchten2p::setThetaS (Vec &thetaSIn)
{
  global_thetaS=thetaSIn; 
//   global_thetaSR = global_thetaS - global_thetaR;
  global_thetaSR = global_thetaS;
  axpy(-1.0,global_thetaR,global_thetaSR);
  thetaS.startSetFromGlobal(global_thetaS);
  thetaS.endSetFromGlobal(global_thetaS);
  thetaSR.startSetFromGlobal(global_thetaSR);
  thetaSR.endSetFromGlobal(global_thetaSR);
}
void VanGenuchten2p::setThetaR (Vec &thetaRIn)
{
  global_thetaR=thetaRIn;
//   global_thetaSR = global_thetaS - global_thetaR;
  global_thetaSR = global_thetaS;
  axpy(-1.0, global_thetaR, global_thetaSR);
  thetaR.startSetFromGlobal(global_thetaR);
  thetaR.endSetFromGlobal(global_thetaR);
  thetaSR.startSetFromGlobal(global_thetaSR);
  thetaSR.endSetFromGlobal(global_thetaSR);
}

void VanGenuchten2p::millerSimilarScaling(Vec& delta)
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
