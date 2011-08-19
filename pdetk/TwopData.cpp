#include "TwopData.h"

namespace Daetk 
{
namespace TwoPhaseFlow
{

TwopData::TwopData(Petsc::SecondOrderFd& s,ParameterDatabase& pd,
		   int dim):
  node(s),
  nNodes(pd.i("nxNodes")*pd.i("nyNodes")*pd.i("nzNodes")),
  nxNodes(pd.i("nxNodes")),
  nyNodes(pd.i("nyNodes")),
  nzNodes(pd.i("nzNodes")),
  dirScale(1.0),
  tForBCreset(-1.0),
  fwIndex(0,dim-1),
  fnIndex(1,dim-1),
  mwIndex(2,dim-1),
  mnIndex(3,dim-1)
{
  Tracer tr("TwopData()");

  thetaS.newsize(Vec::GLOBAL,node.da_);
  thetaS.setExample();
  thetaS = pd.r("thetaS");
  local_thetaS.newsize(Vec::LOCAL,node.da_);
  local_thetaS = pd.r("thetaS");

  local_y.newsize(Vec::LOCAL,node.dadof_);//note, this one has 2 or 4 dof
  local_y=pd.r("thetaS");

  ygdofExample.newsize(Vec::GLOBAL,node.dadof_);
  ygdofExample.setExample();


  D_Div.resize(2);
  local_DK.resize(2);
  for (int alpha=0;alpha<2;alpha++)
    {
      local_DK[alpha].resize(2);
      D_Div[alpha].resize(2);
      for (int beta=0;beta<2;beta++)
        {
          local_DK[alpha][beta].newsize(Vec::LOCAL,node.da_);
          D_Div[alpha][beta].resize(7);
          for (int s=0;s<7;s++)
            {
              D_Div[alpha][beta][s].newsize(nNodes);
              D_Div[alpha][beta][s]=0.0;
            }
        }
      m[alpha].newsize(nNodes);
      Dm[alpha].newsize(nNodes);
      local_mCurrent[alpha].newsize(Vec::LOCAL,node.da_);
      p[alpha].newsize(nNodes);
      local_p[alpha].newsize(Vec::LOCAL,node.da_);
      Dp[alpha].newsize(nNodes);
      theta[alpha].newsize(nNodes);
      Dtheta[alpha].newsize(nNodes);
      div[alpha].newsize(nNodes);
      local_K[alpha].newsize(Vec::LOCAL,node.da_);
      local_theta[alpha].newsize(Vec::LOCAL,node.da_);
      local_rho[alpha].newsize(Vec::LOCAL,node.da_);
      local_Drho[alpha].newsize(Vec::LOCAL,node.da_);

      local_DDrho[alpha].newsize(Vec::LOCAL,node.da_);
      m[alpha]=12345;
      Dm[alpha]=12345;
      local_mCurrent[alpha]=12345;
      p[alpha]=12345;
      local_p[alpha]=12345;
      Dp[alpha] = 12345;
      theta[alpha] = 12345;
      Dtheta[alpha] = 12345;
      theta[alpha] = 12345;
    }

  local_DpC_DthetaW.newsize(Vec::LOCAL,node.da_);
  local_DthetaW_DpC.newsize(Vec::LOCAL,node.da_);
  local_DDthetaW_DDpC.newsize(Vec::LOCAL,node.da_);
  local_DDpC_DDthetaW.newsize(Vec::LOCAL,node.da_);
  if (nyNodes == 1)
    {
      DIM = ONE_D;
      dx = (pd.r("xRight") - pd.r("xLeft")) / (nxNodes-1);
      dz = dy = dx;
      oneOverdz = oneOverdy = oneOverdx = 1.0 / dx;
      if (nzNodes != 1)
        {
          std::cerr<<"Setting nyNodes to 1 in the input file selects a 1D"
		   <<" simulation."<<std::endl
		   <<"Use y as the second dimension if you want to run a 2D"
		   <<" simulation."<<std::endl;
        }
      dirScale = 2*oneOverdx*oneOverdx;
    }
  else if (nzNodes ==1)
    {
      DIM = TWO_D;
      dx = (pd.r("xRight") - pd.r("xLeft")) / (nxNodes-1);
      oneOverdx = 1.0 / dx;
      dy = (pd.r("yBack") - pd.r("yFront")) / (nyNodes-1);
      oneOverdz = oneOverdy = 1.0 / dy;
      dz = dy;
      dirScale = std::max(2*oneOverdx*oneOverdx,2*oneOverdy*oneOverdy);
    }
  else
    {
      dx = (pd.r("xRight") - pd.r("xLeft")) / (nxNodes-1);
      oneOverdx = 1.0 / dx;
      dy = (pd.r("yBack") - pd.r("yFront")) / (nyNodes-1);
      oneOverdy = 1.0 / dy;
      dz = (pd.r("zTop") - pd.r("zBottom")) / (nzNodes-1);
      oneOverdz = 1.0 / dz;
      DIM = THREE_D;
      dirScale = max3(2*oneOverdx*oneOverdx,2*oneOverdy*oneOverdy,
		      2*oneOverdz*oneOverdz);
    }
}

TwopData::~TwopData() 
{
}

}//TwoPhaseFlow
}//Daetk

