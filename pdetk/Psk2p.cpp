#include "Psk2p.h"

namespace Daetk
{
  //mwf added reference to argument
  Psk2p::Psk2p(ParameterDatabase& pd):
    nNodes(pd.i("nxNodes")*pd.i("nyNodes")*pd.i("nzNodes")),
    muW_by_muN(pd.r("muW_by_muN")),
    thetaS(nNodes,pd.r("thetaS")),
    thetaR(nNodes,pd.r("thetaR")),
    KWs(nNodes,pd.r("Ks")),
    thetaSR(nNodes,pd.r("thetaS") - pd.r("thetaR")),
    thetaW_p(0),
    DthetaW_DpC_p(0),
    DDthetaW_DDpC_p(0),
    krW_p(0),
    DkrW_DpC_p(0),
    krN_p(0),
    DkrN_DpC_p(0)
  {}

  Psk2p::Psk2p(): 
    nNodes(0),
    muW_by_muN(0),
    global_thetaS(),
    global_thetaR(),
    global_KWs(),
    global_thetaSR(),
    thetaS(),
    thetaR(),
    KWs(),
    thetaSR(),
    thetaW_p(0),
    DthetaW_DpC_p(0),
    DDthetaW_DDpC_p(0),
    krW_p(0),
    DkrW_DpC_p(0),
    krN_p(0),
    DkrN_DpC_p(0)
  {
    Tracer tr("Psk()");
  }

Psk2p::~Psk2p(){}

  void Psk2p::readParameters(ParameterDatabase& pd)
  {
    Petsc::Sys psys;
    nNodes= pd.i("nxNodes")*pd.i("nyNodes")*pd.i("nzNodes");
    muW_by_muN = pd.r("muW_by_muN");

    global_thetaS.newsize(nNodes);
    global_thetaS = pd.r("thetaS");

    global_thetaR.newsize(nNodes);
    global_thetaR = pd.r("thetaR");

    global_KWs.newsize(nNodes);
    global_KWs = pd.r("Ks");

    global_thetaSR.newsize(nNodes);
    global_thetaSR = pd.r("thetaS") -  pd.r("thetaR");
//      std::cout<<"nNodes in Prams "<<nNodes<<std::endl<<std::flush;
//      std::cout<<"thetaSR in read Params "<<std::endl<<global_thetaSR.dim()<<std::endl<<std::flush<<std::endl<<"da "<<global_thetaS.getDA()<<std::endl<<std::flush;
//      for (int i=global_thetaSR.getGlobalLow();i<global_thetaSR.getGlobalHigh();i++)
//        std::cout<<" global thetaSR "<<global_thetaSR(i)<<std::endl<<std::flush;
    psys.barrier();
    thetaS.newsize(Vec::LOCAL,global_thetaS.getDA());
    thetaS = pd.r("thetaS");
    thetaR.newsize(Vec::LOCAL,global_thetaS.getDA());
    thetaR = pd.r("thetaR");
    KWs.newsize(Vec::LOCAL,global_thetaS.getDA());
    KWs = pd.r("Ks");
    thetaSR.newsize(Vec::LOCAL,global_thetaS.getDA());
    thetaSR = pd.r("thetaS") -  pd.r("thetaR");

//      thetaS.startSetFromGlobal(global_thetaS);
//      thetaS.endSetFromGlobal(global_thetaS);

//      thetaR.startSetFromGlobal(global_thetaR);
//      thetaR.endSetFromGlobal(global_thetaR);

//      KWs.startSetFromGlobal(global_KWs);
//      KWs.endSetFromGlobal(global_KWs);

//      thetaSR.startSetFromGlobal(global_thetaSR);
//      thetaSR.endSetFromGlobal(global_thetaSR);

//      std::cout<<"local thetaSR dim "<<std::endl<<thetaSR.dim()<<std::endl<<std::flush;
//      int j=global_thetaSR.getGlobalLow();
//      for (int i=0;i<thetaSR.dim();i++)
//        {
//          std::cout<<" local thetaSR "<<thetaSR[i]<<'\t'<<global_thetaSR(j)<<std::endl<<std::flush;
//          ++j;
//        }
  }
  
  void Psk2p::attachPsiC(Vec& psiCIn){psiC_p=&psiCIn;}

  void Psk2p::attachThetaW(Vec& thetaWIn){thetaW_p=&thetaWIn;}
  void Psk2p::attachDthetaW_DpC(Vec& DthetaWIn){DthetaW_DpC_p=&DthetaWIn;}
  void Psk2p::attachDDthetaW_DDpC(Vec& DDthetaWIn){DDthetaW_DDpC_p=&DDthetaWIn;}

  void Psk2p::attachKrW(Vec& krwIn){krW_p=&krwIn;}
  void Psk2p::attachDKrW_DpC(Vec& DkrwIn){DkrW_DpC_p=&DkrwIn;}

  void Psk2p::attachKrN(Vec& krnIn){krN_p=&krnIn;}
  void Psk2p::attachDKrN_DpC(Vec& DkrnIn){DkrN_DpC_p=&DkrnIn;}

  void Psk2p::setHeads(const Vec& psiW_vec, const Vec& psiN_vec)
  {
    for(int i=0;i<psiW_vec.ldim();i++)
      {
        setHeads(i,psiW_vec[i],psiN_vec[i]);
        (*thetaW_p)[i] = getThetaW();
        (*krW_p)[i] = getKrW();
        if (DthetaW_DpC_p)
          (*DthetaW_DpC_p)[i]=getDthetaW_DpC();
        if (krN_p)
          (*krN_p)[i]=getKrN();
      }
  }

  void Psk2p::setVFraction(const Vec& thetaW_vec)
  {
    for(int i=0;i<thetaW_vec.ldim();i++)
      {
        setVFraction(thetaW_vec[i],i);
        (*psiC_p)[i] = getPsiC();
        (*krW_p)[i] = getKrW();
        if (DthetaW_DpC_p)
          (*DthetaW_DpC_p)[i]=getDthetaW_DpC();
        if (krN_p)
          (*krN_p)[i]=getKrN();
      }
  }

  void Psk2p::calculateDerivativesHead(const Vec& psiW_vec, const Vec& psiN_vec)
  {
    for(int i=0;i<psiW_vec.ldim();i++)
      {
        setHeads(i,psiW_vec[i],psiN_vec[i]);
        calculateDerivatives();
        (*thetaW_p)[i] = getThetaW();
        (*DthetaW_DpC_p)[i] = getDthetaW_DpC();
        (*krW_p)[i] = getKrW();
        (*DkrW_DpC_p)[i] = getDkrW_DpC();
        
        if (DDthetaW_DDpC_p)
          (*DDthetaW_DDpC_p)[i]=getDDthetaW_DDpC();
        
        if (krN_p)
          {
            (*krN_p)[i]=getKrN();
            (*DkrN_DpC_p)[i]=getDkrN_DpC();
          }
      }
  }

  void Psk2p::calculateDerivativesVFraction(const Vec& thetaW_vec)
  {
    for(int i=0;i<thetaW_vec.ldim();i++)
      {
        setVFraction(thetaW_vec[i],i);
        calculateDerivativesVFraction();
        (*psiC_p)[i] = getPsiC();
        (*DthetaW_DpC_p)[i] = getDthetaW_DpC();
        (*krW_p)[i] = getKrW();
        (*DkrW_DpC_p)[i] = getDkrW_DpC();
        
        if (DDthetaW_DDpC_p)
          (*DDthetaW_DDpC_p)[i]=getDDthetaW_DDpC();
        
        if (krN_p)
          {
            (*krN_p)[i]=getKrN();
            (*DkrN_DpC_p)[i]=getDkrN_DpC();
          }
      }
  }    

void Psk2p::setThetaS(Vec &thetaSIn)
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
void Psk2p::setThetaR (Vec &thetaRIn)
{
  global_thetaR=thetaRIn;
//   global_thetaSR = global_thetaS - global_thetaR;
  global_thetaSR = global_thetaS;
  axpy(-1.0,global_thetaR,global_thetaSR);
  thetaR.startSetFromGlobal(global_thetaR);
  thetaR.endSetFromGlobal(global_thetaR);
  thetaSR.startSetFromGlobal(global_thetaSR);
  thetaSR.endSetFromGlobal(global_thetaSR);
}

void Psk2p::setKWs(Vec& KWsIn)
{
      global_KWs = KWsIn;
      KWs.newsize(Vec::LOCAL,global_KWs.getDA());
      KWs.startSetFromGlobal(global_KWs);
      KWs.endSetFromGlobal(global_KWs);
    }

}//Daetk



