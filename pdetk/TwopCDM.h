#ifndef TWOP_CDM_H
#define TWOP_CDM_H

#include "DaeDefinition.h"
#include "Definitions.h"
#include "Utilities.h"
#include "ParameterDatabase.h"
#include "TwopData.h"
#include "DivKgradUG.h"
#include <fstream>
#include "DaetkPetscSys.h"

namespace Daetk 
{
namespace TwoPhaseFlow 
{

class TwopDaeDef : public TwopData, public DaeDefinition
{
public:
  /***********************************************************************
    base class for TwopNaplWater codes that should take care of 
    TwopData and DaeDefinition interfaces
   ***********************************************************************/

  TwopDaeDef(Petsc::SecondOrderFd& s,ParameterDatabase& pd, 
	     DataCollector& dataIn, int dim);
  virtual ~TwopDaeDef();

  void setInitialConditions(const real& t,const Vec& y) {t0=t; y0 = y;}
  //for implementing discontinuous boundary conditions
  virtual void resetBoundaryConditions(const real& t) 
  { tForBCreset = t;}
  //DaeDefinition interface
  const real& getT0()      {return t0;}
  const Vec&  getY0()      {return y0;}
  const Vec&  getY0prime() {return y0prime;}

  virtual const Vec& getFlux_x(int n) =0;
  virtual const Vec& getFlux_y(int n) =0;
  virtual const Vec& getFlux_z(int n) =0;

protected:
  
  real t0;
  Vec y0,y0prime;
  //data members for DaeDefinition
  const Vec* myY;

  std::ofstream twopOut;
  Petsc::Sys psys;
  
  virtual void setICvars()=0;
  virtual bool calculateCoefficients()=0;
  virtual bool calculateJacCoefficients()=0;

};


template <class PROB, class JAC> 
class TwopCDMSP : public TwopDaeDef
{
public:
  //mwf added to switch out DivKgrad types
  //upwind rel perms
  typedef DivKgradUG<typename PROB::BC,2>   DIVKGRAD;
  typedef DivKgradUG1d<typename PROB::BC,2> DIVKGRAD1d;
  typedef DivKgradUG2d<typename PROB::BC,2> DIVKGRAD2d;
  typedef DivKgradUG3d<typename PROB::BC,2> DIVKGRAD3d;
  //arithmetic avg rel perms
  //typedef DivKgrad<BC,2> DIVKGRAD;
  //typedef DivKgrad1d<BC,2> DIVKGRAD1d;
  //typedef DivKgrad2d<BC,2> DIVKGRAD2d;
  //typedef DivKgrad3d<BC,2> DIVKGRAD3d;

  enum variable {S, P};

  TwopCDMSP(Petsc::SecondOrderFd& s,ParameterDatabase& pd, 
	    DataCollector& dataIn);
  virtual ~TwopCDMSP() {delete dkg;}

  void setICvars();
  bool residual(const real& t,const Vec& y,const Vec& yp, Vec& res);
  bool residual(const VecIndex&,const real& t,const Vec& y,
		const Vec& yp, Vec& res){return true;}  
  bool yPrimeValue(const real& t, const Vec& y, Vec& yp);
  bool evaluateDaeJacobian(const real &t,const  Vec& y,const  Vec& yp,
			   const real& alphaBDF);

  bool jacVec(const Vec& x, Vec& Jx){jacMat->apply(x,Jx); return false;}
  void setJac(JAC& jacRep) {jacMat = &jacRep;}
  //try to put in residual calculation at least
  void stepTaken()         {pc.psk.updateHistory();}


  bool calculateCoefficients();
  bool calculateJacCoefficients();

  virtual const Vec& getFlux_x(int n) { return dkg->getFlux_x(n); }
  virtual const Vec& getFlux_y(int n) { return dkg->getFlux_y(n); }
  virtual const Vec& getFlux_z(int n) { return dkg->getFlux_z(n); }

  PROB pc;
  JAC *jacMat;
  
protected:
  DIVKGRAD* dkg;

private:
  Vec pwtmp,twtmp,pntmp;
};

template<class PROB, class JAC>
TwopCDMSP<PROB,JAC>::TwopCDMSP(Petsc::SecondOrderFd& s,
			       ParameterDatabase& pd,
			       DataCollector& dataIn):
  TwopDaeDef(s,pd,dataIn,
	     2*pd.i("nxNodes")*pd.i("nyNodes")*pd.i("nzNodes")),
  pc(pd),
  jacMat(0),
  dkg(0),
  pwtmp(nNodes),
  twtmp(nNodes),
  pntmp(nNodes)
{
  Tracer tr("TwopCDMSP");
  //setup dkg

  if (DIM == ONE_D)
    {
      dkg = new DIVKGRAD1d(pc.bc,node,nxNodes,pd.r("gx"),oneOverdx);
    }
  else if (DIM == TWO_D)
    {
      dkg = new DIVKGRAD2d(pc.bc,node,nxNodes,nyNodes,
			   pd.r("gx"),pd.r("gy"),
			   oneOverdx,oneOverdy);
    }
  else
    {
      dkg = new DIVKGRAD3d(pc.bc,node,nxNodes,nyNodes,nzNodes,
			   pd.r("gx"),pd.r("gy"),pd.r("gz"),
			   oneOverdx,oneOverdy,oneOverdz);
    }
  assert(dkg);
  dkg->computeInterfaceK(pc.psk.getLocalKWs());

  //If I want to switch approach out to pass in a problem
  //then I'll need to do this stuff outside as well
  //
  pc.setBounds(-10.0*pd.r("m_atol"));//make this more general?

  setICvars();
  pc.setupProblem(pd,*this,&alphaDaeDef);
  //also recalculate interface conductivities
  dkg->computeInterfaceK(pc.psk.getLocalKWs());

//    std::cout<<"thetaW"<<theta[W]<<std::endl;
#ifdef DEBUG
  psys.barrier();
  pc.bc[W].print();
  pc.bc[N].print();
  psys.barrier();
#endif
//    std::cout<<"calc coeff"<<std::endl;
  bool evalError=false;
  psys.barrier();
  evalError = psys.catchError(calculateCoefficients());
  psys.barrier();
  if (evalError)
    {
      std::cerr<<"initial values out of range"<<std::endl<<std::flush;
      exit(1);
    }
  Vec tmp(thetaS);
  tmp.setFromLocal(local_mCurrent[W]);  m[W] = tmp;
  tmp.setFromLocal(local_mCurrent[N]);  m[N] = tmp;
  psys.barrier();
  //    std::cout<<"ypval"<<std::endl;
  yPrimeValue(t0,y0,y0prime);
  psys.barrier();
#ifdef DEBUG
  twopOut<<"y0"<<y0<<std::endl<<std::flush<<"y0'"
	 <<y0prime<<std::flush<<std::endl;
#endif
  betaDaeDef = y0prime;
  Vec temp(y0);
  //    std::cout<<"res cal"<<std::endl;
  psys.barrier();
  residual(t0,y0,y0prime,temp);
  psys.barrier();
#ifdef DEBUG  
  twopOut<<"y0"<<y0<<std::endl<<std::flush<<"y0'"
	 <<y0prime<<std::flush<<std::endl<<"res"<<temp<<std::flush;
#endif
  twopOut<<"Relative Residual of Initial Data: "
	 <<nrm2(temp)/nrm2(y0)<<std::endl;

}

template<class PROB, class JAC>
void TwopCDMSP<PROB,JAC>::setICvars()
{
  myY=&y0;
  theta[W].attachToVecMulti(Vec::REF,y0,fwIndex);
  theta[W].setStrideMulti(2);
  p[W].attachToVecMulti(Vec::REF,y0,fnIndex);
  p[W].setStrideMulti(2);
  Dtheta[W].attachToVecMulti(Vec::REF,y0prime,fwIndex);
  Dtheta[W].setStrideMulti(2);
  Dp[W].attachToVecMulti(Vec::REF,y0prime,fnIndex);
  Dp[W].setStrideMulti(2);
}

template <class PROB, class JAC>
bool TwopCDMSP<PROB,JAC>::residual(const real &t,const  Vec& y,const  Vec& yp, 
				Vec& res)
{  
  myY = &y;

  //split vector into NAPL and water pressures
  theta[W].attachToVecMulti(Vec::REF,y,fwIndex);
  Dtheta[W].attachToVecMulti(Vec::REF,yp,fwIndex);

  p[W].attachToVecMulti(Vec::REF,y,fnIndex);
  Dp[W].attachToVecMulti(Vec::REF,yp,fnIndex);
  
  theta[W].setStrideMulti(2);
  Dtheta[W].setStrideMulti(2);

  p[W].setStrideMulti(2);
  Dp[W].setStrideMulti(2);

  //split res int NAPL and water mass balance residuals

  resF[W].attachToVecMulti(Vec::REF,res,fwIndex);
  resF[N].attachToVecMulti(Vec::REF,res,fnIndex);
  
  resF[W].setStrideMulti(2);
  resF[N].setStrideMulti(2);

  pc.adjustBC(t);
  

  bool evalError=false;
  evalError = psys.catchError(calculateCoefficients());
  if (evalError)
    return evalError;

  // compute res on entire domain
  
  for (int i=0;i<node.local_nzNodes;i++) 
    for (int j=0;j<node.local_nyNodes;j++)
      {
	int k=0;
	node.localIndex(i,j,k);
	for (k=0;k<node.local_nxNodes;k++)
	  {
            resF[W][node.center_noGhost] = 
              (local_rho[W][node.center]*Dtheta[W][node.center_noGhost] + 
               local_Drho[W][node.center]*theta[W][node.center_noGhost]*
               Dp[W][node.center_noGhost]) 
              - div[W][node.center_noGhost];
            //the density derivative is tricky. it IS Dp[W] in the
            //first term since Dp[N]/Dp[W] = +1. Dp[N]/DpC = +1 as well.
            resF[N][node.center_noGhost] =(local_Drho[N][node.center]*
                                           local_theta[N][node.center]*
                                           local_DpC_DthetaW[node.center] - 
                                           local_rho[N][node.center])*
              Dtheta[W][node.center_noGhost] + 
              local_Drho[N][node.center]*local_theta[N][node.center]*
              Dp[W][node.center_noGhost]  
              - div[N][node.center_noGhost];

	    //mwf calculate p[N] for boundary conditions now?
	    p[N][node.center_noGhost] = local_p[N][node.center];
            ++node;
          }
      }
  pc.bc[W].applyDirichletConditions(node,resF[W]);
  pc.bc[N].applyDirichletConditions(node,resF[N]);
  return false;
}

template <class PROB, class JAC>
bool TwopCDMSP<PROB,JAC>::yPrimeValue(const real &t,const  Vec& y,Vec& yp)
{ 
  myY=&y;
  
  pc.adjustBC(t);

  real a,b,c,d,det;
         
  theta[W].attachToVecMulti(Vec::REF,y,fwIndex);
  p[W].attachToVecMulti(Vec::REF,y,fnIndex);
  
  theta[W].setStrideMulti(2);
  p[W].setStrideMulti(2);

  Dtheta[W].attachToVecMulti(Vec::REF,yp,fwIndex);
  Dp[W].attachToVecMulti(Vec::REF,yp,fnIndex);
  
  Dtheta[W].setStrideMulti(2);
  Dp[W].setStrideMulti(2);
  

  
  bool evalError = psys.catchError(calculateCoefficients());
  if (evalError)
    return evalError;
  
  //    compute yp on interior and at left and right boundaries

  int n;
  int end=p[W].getLocalHigh();

  for (n=0;n<end;n++)
    {
      node.localNode(n);
      //solve 2x2 system at each node
      a = local_rho[W][node.center];
      b = local_Drho[W][node.center]*theta[W][n];
      c = (local_Drho[N][node.center]*local_theta[N][node.center]*local_DpC_DthetaW[node.center] - local_rho[N][node.center]);
      d = local_Drho[N][node.center]*local_theta[N][node.center];

      det = a*d - b*c;
      
      if (fabs(det) < MACHINE_EPSILON)
        {
          if (local_DthetaW_DpC[node.center] && 
	      local_Drho[W][node.center]==0 && local_Drho[N][node.center] ==0)
            {
              twopOut<<"system is incompressible and locally unsaturated"<<std::endl;
              Dtheta[W][n] = dkg->getDiv(node,W)/a;
              //set the value of Dp just to give the numerical jac information
              Dp[W][n] = dkg->getDiv(node,W)/(- local_rho[W][node.center]*local_DthetaW_DpC[node.center]);
            }
          else if (!local_DthetaW_DpC[node.center] && local_Drho[W][node.center]==0 && local_Drho[N][node.center] ==0)
            {
              twopOut<<"system is locally saturated and incompressible"<<std::endl;
              Dtheta[W][n] = dkg->getDiv(node,W)/a;
              Dp[W][n] = 0.0;
            }
          else if (local_Drho[W][node.center] && local_Drho[N][node.center]==0) // water phase is compressible
            {
              twopOut<<"water phase is compressible; napl phase is incompressible"<<std::endl;
              Dtheta[W][n] = dkg->getDiv(node,N)/c;
              if (b)
                Dp[W][n] = (dkg->getDiv(node,W) - a*Dtheta[W][n])/b;
              else if (local_DthetaW_DpC[node.center])
                Dp[W][n] = dkg->getDiv(node,W)/(- local_rho[W][node.center]*local_DthetaW_DpC[node.center]);
              else
                Dp[W][n] = 0.0;
            }
          else if (local_Drho[N][node.center] && local_Drho[W][node.center]==0) // napl phase is compressible
            {
              twopOut<<"water phase is incompressible; napl phase is compressible"<<std::endl;
              Dtheta[W][n] = dkg->getDiv(node,W)/a;
              if (d)
                Dp[W][n] = (dkg->getDiv(node,N) - c*Dtheta[W][n])/d;
              else if (local_DthetaW_DpC[node.center])
                Dp[W][n] = dkg->getDiv(node,W)/(- local_rho[W][node.center]*local_DthetaW_DpC[node.center]);
              else
                Dp[W][n] = 0.0;
            }
          else
            {
              twopOut<<"confused by singularity of local system"<<std::endl
                  <<"system appears to be compressible and singular"<<std::endl
                  <<"guessing at the initial derivatives"<<std::endl
                  <<a<<'\t'<<b<<'\t'<<c<<'\t'<<d<<std::endl;
              Dtheta[W][n] = dkg->getDiv(node,W)/a;
              Dp[W][n] = dkg->getDiv(node,W)/b;
            }
        }
      else
        {
          Dtheta[W][n] = d*dkg->getDiv(node,W)/det - b*dkg->getDiv(node,N)/det;
          Dp[W][n] = a*dkg->getDiv(node,N)/det - c*dkg->getDiv(node,W)/det;
        }
      //std::cout<<"ypval "<<Dtheta[W][n]<<'\t'<<Dp[W][n]<<'\t'<<theta[W][n]<<'\t'<<p[W][n]<<std::endl<<std::flush;
    }
  //std::cout<<"D before"<<Dp[W]<<std::endl<<std::flush;
  pc.bc[W].applyDirichletYprime(node);
  pc.bc[N].applyDirichletYprime(node);
  //std::cout<<"D after"<<Dp[W]<<std::endl<<std::flush;
  return false;
}


template <class PROB, class JAC>
bool TwopCDMSP<PROB,JAC>::evaluateDaeJacobian(const real &t,
					   const  Vec& y,
					   const  Vec& yp,
					   const real& alphaBDF) 
{  
  jacMat->zeroAll();
  myY = &y;
  const int W=TwopData::W,N=TwopData::N;
  const int S=TwopCDMSP::S,P=TwopCDMSP::P;
  //split vector into NAPL and water pressures
  theta[W].attachToVecMulti(Vec::REF,y,fwIndex);
  Dtheta[W].attachToVecMulti(Vec::REF,yp,fwIndex);

  p[W].attachToVecMulti(Vec::REF,y,fnIndex);
  Dp[W].attachToVecMulti(Vec::REF,yp,fnIndex);

  theta[W].setStrideMulti(2);
  Dtheta[W].setStrideMulti(2);

  p[W].setStrideMulti(2);
  Dp[W].setStrideMulti(2);

  pc.adjustBC(t);

  bool evalError=false;
  evalError = psys.catchError(calculateJacCoefficients());
  if (evalError)
    return evalError;

  real ypJac[2][2];
  real DypJac[2][2];

  for (int i=node.local_z0;i<node.local_nzNodes+node.local_z0;i++)
    for (int j=node.local_y0;j<node.local_nyNodes+node.local_y0;j++)
      for (int k=node.local_x0;k<node.local_nxNodes+node.local_x0;k++)
        {
          node.localIndex(i-node.local_z0,j-node.local_y0,k-node.local_x0);
          int nl=node.center_noGhost,
            nlg=node.center,
            ghosted_left=node.left,
            ghosted_right=node.right,
            ghosted_front=node.front,
            ghosted_back=node.back,
            ghosted_top=node.top,
            ghosted_bottom=node.bottom;
          node.globalIndex(i,j,k);

          ypJac[W][S] = local_rho[W][nlg];
          ypJac[W][P] = local_Drho[W][nlg]*theta[W][nl];
          ypJac[N][S] = local_Drho[N][nlg]*local_theta[N][nlg]*local_DpC_DthetaW[nlg] - local_rho[N][nlg];
          ypJac[N][P] = local_Drho[N][nlg]*local_theta[N][nlg]; 
          
          DypJac[W][S] = local_Drho[W][nlg]*Dp[W][nl];
          DypJac[W][P] = local_Drho[W][nlg]*Dtheta[W][nl] + theta[W][nl]*local_DDrho[N][nlg]*Dp[W][nl];

          DypJac[N][S] = (local_DDrho[N][nlg]*local_DpC_DthetaW[nlg]*local_theta[N][nlg]*local_DpC_DthetaW[nlg] 
                          -local_Drho[N][nlg]*local_DpC_DthetaW[nlg] 
                          + local_Drho[N][nlg]*local_theta[N][nlg]*local_DDpC_DDthetaW[nlg]
                          - local_Drho[N][nlg]*local_DpC_DthetaW[nlg])*Dtheta[W][nl]
            + (-local_Drho[N][nlg] + local_DDrho[N][nlg]*local_DpC_DthetaW[nlg]*local_theta[N][nlg])*Dp[W][nl];

          DypJac[N][P] = (local_DDrho[N][nlg]*local_theta[N][nlg]*local_DpC_DthetaW[nlg] - local_Drho[N][nlg])*Dtheta[W][nl] 
            + local_DDrho[N][nlg]*local_theta[N][nlg]*Dp[W][nl];
          
          (*jacMat)(2*node.center+W,2*node.center+S)  = DypJac[W][S]
            -D_Div[W][N][DIVKGRAD::CENTER][nl]*local_DpC_DthetaW[nlg]
            + alphaBDF*ypJac[W][S];
              
          (*jacMat)(2*node.center+N,2*node.center+S)  = DypJac[N][S]
            -D_Div[N][N][DIVKGRAD::CENTER][nl]*local_DpC_DthetaW[nlg]
            + alphaBDF*ypJac[N][S];
          
          (*jacMat)(2*node.center+W,2*node.center+P)  = DypJac[W][P]
            -  D_Div[W][W][DIVKGRAD::CENTER][nl]
            -  D_Div[W][N][DIVKGRAD::CENTER][nl]
            + alphaBDF*ypJac[W][P];
          
          (*jacMat)(2*node.center+N,2*node.center+P)  = DypJac[N][P]
            - D_Div[N][W][DIVKGRAD::CENTER][nl]
            - D_Div[N][N][DIVKGRAD::CENTER][nl]
            + alphaBDF*ypJac[N][P];
          
        if (node.anchor->k > 0)
          {
            (*jacMat)(2*node.center+W,2*node.left+S)= -D_Div[W][N][DIVKGRAD::LEFT][nl]*local_DpC_DthetaW[ghosted_left];
            (*jacMat)(2*node.center+N,2*node.left+S)= -D_Div[N][N][DIVKGRAD::LEFT][nl]*local_DpC_DthetaW[ghosted_left];
            (*jacMat)(2*node.center+W,2*node.left+P)= -D_Div[W][W][DIVKGRAD::LEFT][nl]-D_Div[W][N][DIVKGRAD::LEFT][nl];
            (*jacMat)(2*node.center+N,2*node.left+P)= -D_Div[N][W][DIVKGRAD::LEFT][nl]-D_Div[N][N][DIVKGRAD::LEFT][nl];
          }
        if (node.anchor->k < nxNodes-1)
          {
            (*jacMat)(2*node.center+W,2*node.right+S)= -D_Div[W][N][DIVKGRAD::RIGHT][nl]*local_DpC_DthetaW[ghosted_right];
            (*jacMat)(2*node.center+N,2*node.right+S)= -D_Div[N][N][DIVKGRAD::RIGHT][nl]*local_DpC_DthetaW[ghosted_right];
            (*jacMat)(2*node.center+W,2*node.right+P)= -D_Div[W][W][DIVKGRAD::RIGHT][nl]-D_Div[W][N][DIVKGRAD::RIGHT][nl];
            (*jacMat)(2*node.center+N,2*node.right+P)= -D_Div[N][W][DIVKGRAD::RIGHT][nl]-D_Div[N][N][DIVKGRAD::RIGHT][nl];
          }
        
        if (DIM > ONE_D)
          {
            if (node.anchor->j > 0)
              {
                (*jacMat)(2*node.center+W,2*node.front+S)= -D_Div[W][N][DIVKGRAD::FRONT][nl]*local_DpC_DthetaW[ghosted_front];
                (*jacMat)(2*node.center+N,2*node.front+S)= -D_Div[N][N][DIVKGRAD::FRONT][nl]*local_DpC_DthetaW[ghosted_front];
                (*jacMat)(2*node.center+W,2*node.front+P)= -D_Div[W][W][DIVKGRAD::FRONT][nl]-D_Div[W][N][DIVKGRAD::FRONT][nl];
                (*jacMat)(2*node.center+N,2*node.front+P)= -D_Div[N][W][DIVKGRAD::FRONT][nl]-D_Div[N][N][DIVKGRAD::FRONT][nl];
              }
            if (node.anchor->j < nyNodes-1)
              {
                (*jacMat)(2*node.center+W,2*node.back+S)= -D_Div[W][N][DIVKGRAD::BACK][nl]*local_DpC_DthetaW[ghosted_back];
                (*jacMat)(2*node.center+N,2*node.back+S)= -D_Div[N][N][DIVKGRAD::BACK][nl]*local_DpC_DthetaW[ghosted_back];
                (*jacMat)(2*node.center+W,2*node.back+P)= -D_Div[W][W][DIVKGRAD::BACK][nl]-D_Div[W][N][DIVKGRAD::BACK][nl];
                (*jacMat)(2*node.center+N,2*node.back+P)= -D_Div[N][W][DIVKGRAD::BACK][nl]-D_Div[N][N][DIVKGRAD::BACK][nl];
              }
          }
        if (DIM == THREE_D)
          {
            if (node.anchor->i > 0)
              {
                (*jacMat)(2*node.center+W,2*node.bottom+S)= -D_Div[W][N][DIVKGRAD::BOTTOM][nl]*local_DpC_DthetaW[ghosted_bottom];
                (*jacMat)(2*node.center+N,2*node.bottom+S)= -D_Div[N][N][DIVKGRAD::BOTTOM][nl]*local_DpC_DthetaW[ghosted_bottom];
                (*jacMat)(2*node.center+W,2*node.bottom+P)= -D_Div[W][W][DIVKGRAD::BOTTOM][nl]-D_Div[W][N][DIVKGRAD::BOTTOM][nl];
                (*jacMat)(2*node.center+N,2*node.bottom+P)= -D_Div[N][W][DIVKGRAD::BOTTOM][nl]-D_Div[N][N][DIVKGRAD::BOTTOM][nl];
              }
            
            if (node.anchor->i < nzNodes-1)
              {
                (*jacMat)(2*node.center+W,2*node.top+S)= -D_Div[W][N][DIVKGRAD::TOP][nl]*local_DpC_DthetaW[ghosted_top];
                (*jacMat)(2*node.center+N,2*node.top+S)= -D_Div[N][N][DIVKGRAD::TOP][nl]*local_DpC_DthetaW[ghosted_top];
                (*jacMat)(2*node.center+W,2*node.top+P)= -D_Div[W][W][DIVKGRAD::TOP][nl]-D_Div[W][N][DIVKGRAD::TOP][nl];
                (*jacMat)(2*node.center+N,2*node.top+P)= -D_Div[N][W][DIVKGRAD::TOP][nl]-D_Div[N][N][DIVKGRAD::TOP][nl];
              } 
          }
        jacMat->finalizeBlockRow(node.center);
        }
  int dof=2;
  //will need to be fixed or take  a different approach when mixing bc variables.
  pc.bc[W].applyDirichletDerivatives(node,*jacMat,dof,S); 
  pc.bc[N].applyDirichletDerivatives(node,*jacMat,dof,P); 
  //mwf try to include dirichlet conditions for nonwetting phase pressure
  pc.bc[N].adjustMixedDirichletDerivatives(node,*jacMat,dof,P,
					   local_DpC_DthetaW,S);
  jacMat->beginAssembly();
  jacMat->endAssembly();

  return false;
}


template <class PROB, class JAC>
bool TwopCDMSP<PROB,JAC>::calculateCoefficients()
{
  //    std::cout<<"myY address "<<myY
  //<<std::endl<<"myY"<<(*myY)<<std::flush<<std::endl;
  local_y.startSetFromGlobal(*myY);
  local_y.endSetFromGlobal(*myY);
  int n;
  bool evalError=false;
  int end=local_p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      local_theta[W][n] = local_y[2*n];
      local_p[W][n]=local_y[2*n+1];
    }
  for (n=0;n<end;n++)
    { 
      local_theta[N][n] = (local_thetaS[n] - local_theta[W][n]);
      
      if (local_theta[N][n] < pc.mNegative)
        {
          std::cerr<<"1theta["<<n<<"]="<<local_theta[N][n]
		   <<" out of range pc.mNegative="<<pc.mNegative<<std::endl;
          evalError = true;
        }
      else if (local_theta[W][n] <= pc.psk.getLocalThetaR()[n])
        {
          std::cerr<<"2theta["<<n<<"]="<<local_theta[W][n]
		   <<" out of range"<<std::endl;
          evalError =  true;
        }
      else
        pc.psk.setVFraction(local_theta[W][n],n);
    
       if (fabs(local_p[W][n]) > pc.pBig)
        {
          std::cerr<<"pW["<<n<<"]="<<local_p[W][n]
		   <<" out of range"<<std::endl;
          evalError = true;
        }
      else
        {
          local_p[N][n] = pc.psk.getPsiC() + local_p[W][n];
          
          pc.densityW.setHead(local_p[W][n]);
          pc.densityN.setHead(local_p[N][n]);
        }

      local_DthetaW_DpC[n] = pc.psk.getDthetaW_DpC();

      if (fabs(local_DthetaW_DpC[n]) > 0.0)
        local_DpC_DthetaW[n] = 1.0/local_DthetaW_DpC[n];
      else
        local_DpC_DthetaW[n] = 0.0;

      local_K[W][n] = pc.psk.getKrW();
      local_K[N][n] = pc.psk.getKrN();

      local_rho[W][n] = pc.densityW.getRho();
      local_rho[N][n] = pc.densityN.getRho();

      local_Drho[W][n] = pc.densityW.getDrho();
      local_Drho[N][n] = pc.densityN.getDrho();
    }      
  
  dkg->computeDivergence(local_K,local_rho,local_p,div);

  return evalError;
}

template <class PROB, class JAC>
bool TwopCDMSP<PROB,JAC>::calculateJacCoefficients()
{
  local_y.startSetFromGlobal(*myY);
  local_y.endSetFromGlobal(*myY);
  bool evalError=false;
  int n;
  int end=local_p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      local_theta[W][n] = local_y[2*n];
      local_p[W][n]=local_y[2*n+1];
    }
  for (n=0;n<end;n++)
    {
      local_theta[N][n] = (local_thetaS[n] - local_theta[W][n]);
      
      if (local_theta[N][n] < pc.mNegative)
        {
          std::cerr<<"1theta("<<n<<")="<<local_theta[N][n]
		   <<" out of range pc.mNegative="<<pc.mNegative<<std::endl;
          evalError = true;
        }
      else if (local_theta[W][n] <= pc.psk.getLocalThetaR()[n])
        {
          std::cerr<<"2theta("<<n<<")="<<local_theta[W][n]
		   <<" out of range"<<std::endl;
          evalError =  true;
        }
      else
        pc.psk.setVFraction(local_theta[W][n],n);
      
      if (fabs(local_p[W][n]) > pc.pBig)
        {
          std::cerr<<"pW("<<n<<")="<<local_p[W][n]<<" out of range"<<std::endl;
          evalError = true;
        }
      else
        {
          local_p[N][n] = pc.psk.getPsiC() + local_p[W][n];
          
          pc.densityW.setHead(local_p[W][n]);
          pc.densityN.setHead(local_p[N][n]);
        }

      local_DthetaW_DpC[n] = pc.psk.getDthetaW_DpC();

      if (fabs(local_DthetaW_DpC[n]) > 0)
        local_DpC_DthetaW[n] = 1.0/local_DthetaW_DpC[n];
      else
        local_DpC_DthetaW[n] = 0.0;

      local_K[W][n] = pc.psk.getKrW();
      local_K[N][n] = pc.psk.getKrN();

      local_rho[W][n] = pc.densityW.getRho();
      local_rho[N][n] = pc.densityN.getRho();

      local_Drho[W][n] = pc.densityW.getDrho();
      local_Drho[N][n] = pc.densityN.getDrho();


      pc.psk.calculateDerivativesVFraction();
      local_DDthetaW_DDpC[n] = pc.psk.getDDthetaW_DDpC();
      local_DDpC_DDthetaW[n] = -local_DpC_DthetaW[n]*local_DpC_DthetaW[n]*local_DDthetaW_DDpC[n]*local_DpC_DthetaW[n];

      //these are derivatives w.r.t. pw,pn
      local_DK[W][W][n] = -pc.psk.getDkrW_DpC();
      local_DK[W][N][n] = pc.psk.getDkrW_DpC();
      local_DK[N][W][n] = -pc.psk.getDkrN_DpC();
      local_DK[N][N][n] = pc.psk.getDkrN_DpC();
      //mwf debug
      //std::cout<<"TwopCDMSP calcJacCoefs n= "<<n<<" psk.getDkrN_DpC()= "
      //       <<psk.getDkrN_DpC()<<std::endl;
      local_DDrho[W][n] = pc.densityW.getDDrho();
      local_DDrho[N][n] = pc.densityN.getDDrho();
    }

  dkg->computeDivergenceJac(local_K,local_rho,local_p,
                            local_DK,local_Drho,D_Div);
  return evalError;
}


template <class PROB, class JAC> 
class TwopCDMMS : public TwopDaeDef
{
public:
  //mwf added to switch out DivKgrad types
  //upwind rel perms
  typedef DivKgradUG<typename PROB::BC,2>   DIVKGRAD;
  typedef DivKgradUG1d<typename PROB::BC,2> DIVKGRAD1d;
  typedef DivKgradUG2d<typename PROB::BC,2> DIVKGRAD2d;
  typedef DivKgradUG3d<typename PROB::BC,2> DIVKGRAD3d;
  //arithmetic avg rel perms
  //typedef DivKgrad<BC,2> DIVKGRAD;
  //typedef DivKgrad1d<BC,2> DIVKGRAD1d;
  //typedef DivKgrad2d<BC,2> DIVKGRAD2d;
  //typedef DivKgrad3d<BC,2> DIVKGRAD3d;

  enum variable {S, P};

  TwopCDMMS(Petsc::SecondOrderFd& s,ParameterDatabase& pd, 
	    DataCollector& dataIn);
  virtual ~TwopCDMMS() {delete dkg;}
  void setICvars();
  bool residual(const real& t,const Vec& y,const Vec& yp, Vec& res);
  bool residual(const VecIndex&,const real& t,const Vec& y,
                       const Vec& yp, Vec& res){return true;}  
  bool yPrimeValue(const real& t, const Vec& y, Vec& yp);
  bool evaluateDaeJacobian(const real &t,const  Vec& y,const  Vec& yp,
                                    const real& alphaBDF);
  bool jacVec(const Vec& x, Vec& Jx){jacMat->apply(x,Jx); return false;}
  void setJac(JAC& jacRep) {jacMat = &jacRep;}
  //try to put in residual calculation at least
  void stepTaken()         {pc.psk.updateHistory();}

  void correctArgument(Vec& correction);
  bool calculateCoefficients();
  bool calculateJacCoefficients();

  virtual const Vec& getFlux_x(int n) { return dkg->getFlux_x(n); }
  virtual const Vec& getFlux_y(int n) { return dkg->getFlux_y(n); }
  virtual const Vec& getFlux_z(int n) { return dkg->getFlux_z(n); }

  PROB pc;
  JAC *jacMat;

protected:
  DIVKGRAD* dkg;

private:
  Vec beta[2],
    thetaWcor,
    pWcor,
    mWcor,
    mNcor;
};

template<class PROB, class JAC>
TwopCDMMS<PROB,JAC>::TwopCDMMS(Petsc::SecondOrderFd& s,
			       ParameterDatabase& pd,
			       DataCollector& dataIn):
  TwopDaeDef(s,pd,dataIn,
	     2*2*pd.i("nxNodes")*pd.i("nyNodes")*pd.i("nzNodes")),
  pc(pd),
  jacMat(0),
  dkg(0),
  thetaWcor(nNodes),
  pWcor(nNodes),
  mWcor(nNodes),
  mNcor(nNodes)
{
  Tracer tr("TwopCDMMS");
  //setup dkg

  if (DIM == ONE_D)
    {
      dkg = new DIVKGRAD1d(pc.bc,node,nxNodes,pd.r("gx"),oneOverdx);
    }
  else if (DIM == TWO_D)
    {
      dkg = new DIVKGRAD2d(pc.bc,node,nxNodes,nyNodes,
			   pd.r("gx"),pd.r("gy"),
			   oneOverdx,oneOverdy);
    }
  else
    {
      dkg = new DIVKGRAD3d(pc.bc,node,nxNodes,nyNodes,nzNodes,
			   pd.r("gx"),pd.r("gy"),pd.r("gz"),
			   oneOverdx,oneOverdy,oneOverdz);
    }
  assert(dkg);
  dkg->computeInterfaceK(pc.psk.getLocalKWs());

  //If I want to switch approach out to pass in a problem
  //then I'll need to do this stuff outside as well
  //
  pc.setBounds(-10.0*pd.r("m_atol"));//make this more general?

  setICvars();
  pc.setupProblem(pd,*this,&alphaDaeDef);
  //also recalculate interface conductivities
  dkg->computeInterfaceK(pc.psk.getLocalKWs());

//    std::cout<<"thetaW"<<theta[W]<<std::endl;
#ifdef DEBUG
  psys.barrier();
  pc.bc[W].print();
  pc.bc[N].print();
  psys.barrier();
#endif
//    std::cout<<"calc coeff"<<std::endl;
  bool evalError=false;
  psys.barrier();
  evalError = psys.catchError(calculateCoefficients());
  psys.barrier();
  if (evalError)
    {
      std::cerr<<"initial values out of range"<<std::endl<<std::flush;
      exit(1);
    }
  Vec tmp(thetaS);
  tmp.setFromLocal(local_mCurrent[W]);  m[W] = tmp;
  tmp.setFromLocal(local_mCurrent[N]);  m[N] = tmp;
  psys.barrier();
//    std::cout<<"ypval"<<std::endl;
  yPrimeValue(t0,y0,y0prime);
  psys.barrier();
#ifdef DEBUG
  twopOut<<"y0"<<y0<<std::endl<<std::flush<<"y0'"
	 <<y0prime<<std::flush<<std::endl;
#endif
  betaDaeDef = y0prime;
  Vec temp(y0);
  //    std::cout<<"res cal"<<std::endl;
  psys.barrier();
  residual(t0,y0,y0prime,temp);
  psys.barrier();
#ifdef DEBUG  
  twopOut<<"y0"<<y0<<std::endl<<std::flush<<"y0'"
	 <<y0prime<<std::flush<<std::endl<<"res"<<temp<<std::flush;
#endif
  twopOut<<"Relative Residual of Initial Data: "
	 <<nrm2(temp)/nrm2(y0)<<std::endl;


}

template<class PROB, class JAC>
void TwopCDMMS<PROB,JAC>::setICvars()
{
  myY=&y0;
  theta[W].attachToVecMulti(Vec::REF,y0,fwIndex);  theta[W].setStrideMulti(4);
  p[W].attachToVecMulti(Vec::REF,y0,fnIndex);      p[W].setStrideMulti(4);
  m[W].attachToVecMulti(Vec::REF,y0,mwIndex);      m[W].setStrideMulti(4);
  m[N].attachToVecMulti(Vec::REF,y0,mnIndex);      m[N].setStrideMulti(4);

  Dtheta[W].attachToVecMulti(Vec::REF,y0prime,fwIndex);  
  Dtheta[W].setStrideMulti(4);
  Dp[W].attachToVecMulti(Vec::REF,y0prime,fnIndex);      
  Dp[W].setStrideMulti(4);
  Dm[W].attachToVecMulti(Vec::REF,y0prime,mwIndex);      
  Dm[W].setStrideMulti(4);
  Dm[N].attachToVecMulti(Vec::REF,y0prime,mnIndex);     
  Dm[N].setStrideMulti(4);
}


template <class PROB, class JAC>
void TwopCDMMS<PROB,JAC>::correctArgument(Vec& correction)
{
  yLast = yDaeDef;
  ypLast = ypDaeDef;
  Flast = Fcurrent;

  updateJac=true;
  updateF=true;

  //split vector into NAPL and water pressures
  theta[W].attachToVecMulti(Vec::REF,yDaeDef,fwIndex);   
  theta[W].setStrideMulti(4);
  Dtheta[W].attachToVecMulti(Vec::REF,ypDaeDef,fwIndex); 
  Dtheta[W].setStrideMulti(4);

  p[W].attachToVecMulti(Vec::REF,yDaeDef,fnIndex);       
  p[W].setStrideMulti(4);
  Dp[W].attachToVecMulti(Vec::REF,ypDaeDef,fnIndex);     
  Dp[W].setStrideMulti(4);

  m[W].attachToVecMulti(Vec::REF,yDaeDef,mwIndex);       
  m[W].setStrideMulti(4);
  Dm[W].attachToVecMulti(Vec::REF,ypDaeDef,mwIndex);     
  Dm[W].setStrideMulti(4);

  m[N].attachToVecMulti(Vec::REF,yDaeDef,mnIndex);       
  m[N].setStrideMulti(4);
  Dm[N].attachToVecMulti(Vec::REF,ypDaeDef,mnIndex);     
  Dm[N].setStrideMulti(4);
  
  thetaWcor.attachToVecMulti(Vec::REF,correction,fwIndex);           
  thetaWcor.setStrideMulti(4);  
  pWcor.attachToVecMulti(Vec::REF,correction,fnIndex);               
  pWcor.setStrideMulti(4);      
  mWcor.attachToVecMulti(Vec::REF,correction,mwIndex);               
  mWcor.setStrideMulti(4);      
  mNcor.attachToVecMulti(Vec::REF,correction,mnIndex);               
  mNcor.setStrideMulti(4);      

#ifndef USE_BLAS
  theta[W]-=thetaWcor;
  p[W]-=pWcor;
#else
  axpy(-1.0,thetaWcor,theta[W]);
  axpy(-1.0,pWcor,p[W]);
#endif
#ifndef USE_BLAS
  Dtheta[W]-=alphaDaeDef*thetaWcor;
  Dp[W]-=alphaDaeDef*pWcor;
#else
  axpy(-alphaDaeDef,thetaWcor,Dtheta[W]);
  axpy(-alphaDaeDef,pWcor,Dp[W]);
#endif
  int n;
  int end=p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      node.localNode(n);
      mWcor[n] += local_rho[W][node.center]*thetaWcor[n] 
        + local_Drho[W][node.center]*theta[W][n]*pWcor[n];      
      m[W][n]  -= mWcor[n];
      Dm[W][n] -= alphaDaeDef*mWcor[n];

      mNcor[n] += (local_Drho[N][node.center]
		   *local_theta[N][node.center]
		   *local_DpC_DthetaW[node.center] 
		   - local_rho[N][node.center])*thetaWcor[n] 
        + local_Drho[N][node.center]*local_theta[N][node.center]*pWcor[n];
      m[N][n]  -= mNcor[n];
      Dm[N][n] -= alphaDaeDef*mNcor[n];
    }
}

template <class PROB, class JAC>
bool TwopCDMMS<PROB,JAC>::residual(const real &t,const  Vec& y,const  Vec& yp, 
				   Vec& res)
{
  myY=&y;

  theta[W].attachToVecMulti(Vec::REF,y,fwIndex);       theta[W].setStrideMulti(4);
  Dtheta[W].attachToVecMulti(Vec::REF,yp,fwIndex);     Dtheta[W].setStrideMulti(4); 

  p[W].attachToVecMulti(Vec::REF,y,fnIndex);           p[W].setStrideMulti(4);        
  Dp[W].attachToVecMulti(Vec::REF,yp,fnIndex);         Dp[W].setStrideMulti(4);      

  m[W].attachToVecMulti(Vec::REF,y,mwIndex);           m[W].setStrideMulti(4);        
  Dm[W].attachToVecMulti(Vec::REF,yp,mwIndex);         Dm[W].setStrideMulti(4);     
  
  m[N].attachToVecMulti(Vec::REF,y,mnIndex);           m[N].setStrideMulti(4);        
  Dm[N].attachToVecMulti(Vec::REF,yp,mnIndex);         Dm[N].setStrideMulti(4);      

  resF[W].attachToVecMulti(Vec::REF,res,fwIndex);        resF[W].setStrideMulti(4);   
  resF[N].attachToVecMulti(Vec::REF,res,fnIndex);        resF[N].setStrideMulti(4);  

  resM[W].attachToVecMulti(Vec::REF,res,mwIndex);        resM[W].setStrideMulti(4);   
  resM[N].attachToVecMulti(Vec::REF,res,mnIndex);        resM[N].setStrideMulti(4);   
  
  beta[W].attachToVecMulti(Vec::REF,betaDaeDef,mwIndex); beta[W].setStrideMulti(4);   
  beta[N].attachToVecMulti(Vec::REF,betaDaeDef,mnIndex); beta[N].setStrideMulti(4);    

  pc.adjustBC(t);

  bool evalError=false;
  evalError = psys.catchError(calculateCoefficients());
  if (evalError)
    return evalError;

  int n;
  int end=p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      node.localNode(n);
      resM[W][n] = m[W][n] - local_mCurrent[W][node.center];
      resM[N][n] = m[N][n] - local_mCurrent[N][node.center];
      resF[W][n] = alphaDaeDef*local_mCurrent[W][node.center] + beta[W][n] 
	- dkg->getDiv(node,W);
      resF[N][n] = alphaDaeDef*local_mCurrent[N][node.center] + beta[N][n] 
	- dkg->getDiv(node,N);

      //mwf calculate p[N] for boundary conditions now?
      p[N][n] = local_p[N][n];

    }

  pc.bc[W].applyDirichletConditions(node,resF[W]);
  pc.bc[N].applyDirichletConditions(node,resF[N]);

  return false;
}

template <class PROB, class JAC>
bool TwopCDMMS<PROB,JAC>::yPrimeValue(const real &,const  Vec& y,Vec& yp)
{  
  myY=&y;
  real a,b,c,d,det;
             
  theta[W].attachToVecMulti(Vec::REF,y,fwIndex);       theta[W].setStrideMulti(4);
  Dtheta[W].attachToVecMulti(Vec::REF,yp,fwIndex);     Dtheta[W].setStrideMulti(4); 

  p[W].attachToVecMulti(Vec::REF,y,fnIndex);           p[W].setStrideMulti(4);        
  Dp[W].attachToVecMulti(Vec::REF,yp,fnIndex);         Dp[W].setStrideMulti(4);     

  m[W].attachToVecMulti(Vec::REF,y,mwIndex);           m[W].setStrideMulti(4);       
  m[N].attachToVecMulti(Vec::REF,y,mnIndex);           m[N].setStrideMulti(4);        

  Dm[W].attachToVecMulti(Vec::REF,yp,mwIndex);         Dm[W].setStrideMulti(4);     
  Dm[N].attachToVecMulti(Vec::REF,yp,mnIndex);         Dm[N].setStrideMulti(4);      

  bool evalError = psys.catchError(calculateCoefficients());
  if (evalError)
    return evalError;

  //    compute yp on interior and at left and right boundaries

  int n;
  int end=p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      node.localNode(n);
      //solve 2x2 system at each node
      a = local_rho[W][node.center];
      b = local_Drho[W][node.center]*theta[W][n];
      c = (local_Drho[N][node.center]*local_theta[N][node.center]*local_DpC_DthetaW[node.center] - local_rho[N][node.center]);
      d = local_Drho[N][node.center]*local_theta[N][node.center];

      Dm[W][n] = dkg->getDiv(node,W);
      Dm[N][n] = dkg->getDiv(node,N);

      det = a*d - b*c;
      
      if (fabs(det) < MACHINE_EPSILON)
        {
          if (local_DthetaW_DpC[node.center] && local_Drho[W][node.center]==0 && 
	      local_Drho[N][node.center] ==0)
            {
              twopOut<<"system is incompressible and locally unsaturated"<<std::endl;
              Dtheta[W][n] = Dm[W][n]/a;
              //set the value of Dp just to give the numerical jac information
              Dp[W][n] = Dm[W][n]/
		(- local_rho[W][node.center]*local_DthetaW_DpC[node.center]);
            }
          else if (!local_DthetaW_DpC[node.center] && 
		   local_Drho[W][node.center]==0 && local_Drho[N][node.center] ==0)
            {
              twopOut<<"system is locally saturated and incompressible"<<std::endl;
              Dtheta[W][n] = Dm[W][n]/a;
              Dp[W][n] = 0.0;
            }
          else if (local_Drho[W][node.center] && 
		   local_Drho[N][node.center]==0) // water phase is compressible
            {
              twopOut<<"water phase is compressible; napl phase is incompressible"<<std::endl;
              Dtheta[W][n] = Dm[N][n]/c;
              if (b)
                Dp[W][node.center] = (Dm[W][n] - a*Dtheta[W][n])/b;
              else if (local_DthetaW_DpC[node.center])
                Dp[W][n] = Dm[W][n]/(- local_rho[W][node.center]*local_DthetaW_DpC[node.center]);
              else
                Dp[W][n] = 0.0;
            }
          else if (local_Drho[N][node.center] && 
		   local_Drho[W][node.center]==0) // napl phase is compressible
            {
              twopOut<<"water phase is incompressible; napl phase is compressible"<<std::endl;
              Dtheta[W][n] = Dm[W][n]/a;
              if (d)
                Dp[W][n] = (Dm[N][n] - c*Dtheta[W][n])/d;
              else if (local_DthetaW_DpC[n])
                Dp[W][n] = Dm[W][n]/(- local_rho[W][node.center]*local_DthetaW_DpC[node.center]);
              else
                Dp[W][n] = 0.0;
            }
          else
            {
              twopOut<<"confused by singularity of local system"<<std::endl
                  <<"system appears to be compressible and singular"<<std::endl
                  <<"guessing at the initial derivatives"<<std::endl
                  <<a<<'\t'<<b<<'\t'<<c<<'\t'<<d<<std::endl;
              Dtheta[W][n] = Dm[W][n]/a;
              Dp[W][n] = Dm[W][n]/b;
            }
        }
      else
        {
          Dtheta[W][n] = d*Dm[W][n]/det - b*Dm[N][n]/det;
          Dp[W][n] = a*Dm[N][n]/det - c*Dm[W][n]/det;
        }
      //std::cout<<"ypval "<<Dtheta[W][n]<<'\t'<<Dp[W][n]<<'\t'
      //<<theta[W][n]<<'\t'<<p[W][n]<<std::endl<<std::flush;
    }
  //std::cout<<"D before"<<Dp[W]<<std::endl<<std::flush;
  pc.bc[W].applyDirichletYprime(node);
  pc.bc[N].applyDirichletYprime(node);
  //std::cout<<"D after"<<Dp[W]<<std::endl<<std::flush;
  return false;
}


template <class PROB, class JAC>
bool TwopCDMMS<PROB,JAC>::evaluateDaeJacobian(const real &t,
					      const  Vec& y,
					      const  Vec& yp,
					      const real& alphaBDF) 
{  
  jacMat->zeroAll();
  myY = &y;
  const int W=TwopData::W,N=TwopData::N;
  const int S=TwopCDMMS::S,P=TwopCDMMS::P;
  //split vector into NAPL and water pressures
  theta[W].attachToVecMulti(Vec::REF,y,fwIndex);       theta[W].setStrideMulti(4);
  Dtheta[W].attachToVecMulti(Vec::REF,yp,fwIndex);     Dtheta[W].setStrideMulti(4); 

  p[W].attachToVecMulti(Vec::REF,y,fnIndex);           p[W].setStrideMulti(4);        
  Dp[W].attachToVecMulti(Vec::REF,yp,fnIndex);         Dp[W].setStrideMulti(4);      

  m[W].attachToVecMulti(Vec::REF,y,mwIndex);           m[W].setStrideMulti(4);       
  Dm[W].attachToVecMulti(Vec::REF,yp,mwIndex);         Dm[W].setStrideMulti(4);     

  m[N].attachToVecMulti(Vec::REF,y,mnIndex);           m[N].setStrideMulti(4);        
  Dm[N].attachToVecMulti(Vec::REF,yp,mnIndex);         Dm[N].setStrideMulti(4);      

  pc.adjustBC(t);

  bool evalError=false;
  evalError = psys.catchError(calculateJacCoefficients());
  if (evalError)
    return evalError;

  real ypJac[2][2];
  for (int i=node.local_z0;i<node.local_nzNodes+node.local_z0;i++)
    for (int j=node.local_y0;j<node.local_nyNodes+node.local_y0;j++)
      for (int k=node.local_x0;k<node.local_nxNodes+node.local_x0;k++)
        {
          node.localIndex(i-node.local_z0,j-node.local_y0,k-node.local_x0);
          int nlg=node.center,
            ghosted_left=node.left,
            ghosted_right=node.right,
            ghosted_front=node.front,
            ghosted_back=node.back,
            ghosted_top=node.top,
            ghosted_bottom=node.bottom;
          node(i,j,k);
          int nl=node.center-node.globalLow;

          ypJac[W][S] = local_rho[W][nlg];
          ypJac[W][P] = local_Drho[W][nlg]*theta[W][nl];
          ypJac[N][S] = local_Drho[N][nlg]*local_theta[N][nlg]*local_DpC_DthetaW[nl] 
	    - local_rho[N][nlg];
          ypJac[N][P] = local_Drho[N][nlg]*local_theta[N][nlg]; 
          

          (*jacMat)(2*node.center+W,2*node.center+S)  = alphaBDF*ypJac[W][S]
            -D_Div[W][N][DIVKGRAD::CENTER][nl]*local_DpC_DthetaW[nlg];
          
          (*jacMat)(2*node.center+N,2*node.center+S)  = alphaBDF*ypJac[N][S]
            -D_Div[N][N][DIVKGRAD::CENTER][nl]*local_DpC_DthetaW[nlg];

          (*jacMat)(2*node.center+W,2*node.center+P)  = alphaBDF*ypJac[W][P]
            - D_Div[W][W][DIVKGRAD::CENTER][nl]
            - D_Div[W][N][DIVKGRAD::CENTER][nl];
          
          (*jacMat)(2*node.center+N,2*node.center+P)  = alphaBDF*ypJac[N][P]
            - D_Div[N][W][DIVKGRAD::CENTER][nl]
            - D_Div[N][N][DIVKGRAD::CENTER][nl];
        
        if (node.anchor->k > 0)
          {
            (*jacMat)(2*node.center+W,2*node.left+S)= 
	      -D_Div[W][N][DIVKGRAD::LEFT][nl]*local_DpC_DthetaW[ghosted_left];
            (*jacMat)(2*node.center+N,2*node.left+S)= 
	      -D_Div[N][N][DIVKGRAD::LEFT][nl]*local_DpC_DthetaW[ghosted_left];
            (*jacMat)(2*node.center+W,2*node.left+P)= 
	      -D_Div[W][W][DIVKGRAD::LEFT][nl]-D_Div[W][N][DIVKGRAD::LEFT][nl];
            (*jacMat)(2*node.center+N,2*node.left+P)= 
	      -D_Div[N][W][DIVKGRAD::LEFT][nl]-D_Div[N][N][DIVKGRAD::LEFT][nl];
          }
        if (node.anchor->k < nxNodes-1)
          {
            (*jacMat)(2*node.center+W,2*node.right+S)= 
	      -D_Div[W][N][DIVKGRAD::RIGHT][nl]*local_DpC_DthetaW[ghosted_right];
            (*jacMat)(2*node.center+N,2*node.right+S)= 
	      -D_Div[N][N][DIVKGRAD::RIGHT][nl]*local_DpC_DthetaW[ghosted_right];
            (*jacMat)(2*node.center+W,2*node.right+P)= 
	      -D_Div[W][W][DIVKGRAD::RIGHT][nl]-D_Div[W][N][DIVKGRAD::RIGHT][nl];
            (*jacMat)(2*node.center+N,2*node.right+P)= 
	      -D_Div[N][W][DIVKGRAD::RIGHT][nl]-D_Div[N][N][DIVKGRAD::RIGHT][nl];
          }
        
        if (DIM > ONE_D)
          {
            if (node.anchor->j > 0)
              {
                (*jacMat)(2*node.center+W,2*node.front+S)= 
		  -D_Div[W][N][DIVKGRAD::FRONT][nl]*local_DpC_DthetaW[ghosted_front];
                (*jacMat)(2*node.center+N,2*node.front+S)= 
		  -D_Div[N][N][DIVKGRAD::FRONT][nl]*local_DpC_DthetaW[ghosted_front];
                (*jacMat)(2*node.center+W,2*node.front+P)= 
		  -D_Div[W][W][DIVKGRAD::FRONT][nl]-D_Div[W][N][DIVKGRAD::FRONT][nl];
                (*jacMat)(2*node.center+N,2*node.front+P)= 
		  -D_Div[N][W][DIVKGRAD::FRONT][nl]-D_Div[N][N][DIVKGRAD::FRONT][nl];
              }
            if (node.anchor->j < nyNodes-1)
              {
                (*jacMat)(2*node.center+W,2*node.back+S)= 
		  -D_Div[W][N][DIVKGRAD::BACK][nl]*local_DpC_DthetaW[ghosted_back];
                (*jacMat)(2*node.center+N,2*node.back+S)= 
		  -D_Div[N][N][DIVKGRAD::BACK][nl]*local_DpC_DthetaW[ghosted_back];
                (*jacMat)(2*node.center+W,2*node.back+P)= 
		  -D_Div[W][W][DIVKGRAD::BACK][nl]-D_Div[W][N][DIVKGRAD::BACK][nl];
                (*jacMat)(2*node.center+N,2*node.back+P)= 
		  -D_Div[N][W][DIVKGRAD::BACK][nl]-D_Div[N][N][DIVKGRAD::BACK][nl];
              }
          }
        if (DIM == THREE_D)
          {
            if (node.anchor->i > 0)
              {
                (*jacMat)(2*node.center+W,2*node.bottom+S)= 
		  -D_Div[W][N][DIVKGRAD::BOTTOM][nl]*local_DpC_DthetaW[ghosted_bottom];
                (*jacMat)(2*node.center+N,2*node.bottom+S)= 
		  -D_Div[N][N][DIVKGRAD::BOTTOM][nl]*local_DpC_DthetaW[ghosted_bottom];
                (*jacMat)(2*node.center+W,2*node.bottom+P)= 
		  -D_Div[W][W][DIVKGRAD::BOTTOM][nl]-D_Div[W][N][DIVKGRAD::BOTTOM][nl];
                (*jacMat)(2*node.center+N,2*node.bottom+P)= 
		  -D_Div[N][W][DIVKGRAD::BOTTOM][nl]-D_Div[N][N][DIVKGRAD::BOTTOM][nl];
              }
            
            if (node.anchor->i < nzNodes-1)
              {
                (*jacMat)(2*node.center+W,2*node.top+S)= 
		  -D_Div[W][N][DIVKGRAD::TOP][nl]*local_DpC_DthetaW[ghosted_top];
                (*jacMat)(2*node.center+N,2*node.top+S)= 
		  -D_Div[N][N][DIVKGRAD::TOP][nl]*local_DpC_DthetaW[ghosted_top];
                (*jacMat)(2*node.center+W,2*node.top+P)= 
		  -D_Div[W][W][DIVKGRAD::TOP][nl]-D_Div[W][N][DIVKGRAD::TOP][nl];
                (*jacMat)(2*node.center+N,2*node.top+P)= 
		  -D_Div[N][W][DIVKGRAD::TOP][nl]-D_Div[N][N][DIVKGRAD::TOP][nl];
              } 
          }
        jacMat->finalizeBlockRow(node.center);
      }
  
  //will need to be fixed or take  a different approach when mixing bc variables.
  pc.bc[W].applyDirichletDerivatives(node,*jacMat,2,S); 
  pc.bc[N].applyDirichletDerivatives(node,*jacMat,2,P); 
  //mwf try to include dirichlet conditions for nonwetting phase pressure
  pc.bc[N].adjustMixedDirichletDerivatives(node,*jacMat,2,P,
					   local_DpC_DthetaW,S);
  jacMat->beginAssembly();
  jacMat->endAssembly();

  return false;
}


template <class PROB, class JAC>
bool TwopCDMMS<PROB,JAC>::calculateCoefficients()
{
  local_y.startSetFromGlobal(*myY);
  local_y.endSetFromGlobal(*myY);

  bool evalError=false;
  int n;
  int end=local_p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      local_theta[W][n] = local_y[4*n];
      local_p[W][n]=local_y[4*n+1];
    }
  for (n=0;n<end;n++)
    { 
      local_theta[N][n] = (local_thetaS[n] - local_theta[W][n]);
      
      if (local_theta[N][n] < pc.mNegative)
        {
          std::cerr<<"1theta["<<n<<"]="<<local_theta[N][n]<<" out of range"<<std::endl;
          evalError = true;
        }
      else if (local_theta[W][n] <= pc.psk.getLocalThetaR()[n])
        {
          std::cerr<<"2theta["<<n<<"]="<<local_theta[W][n]<<" out of range"<<std::endl;
          evalError =  true;
        }
      else
        pc.psk.setVFraction(local_theta[W][n],n);
    
       if (fabs(local_p[W][n]) > pc.pBig)
        {
          std::cerr<<"pW["<<n<<"]="<<local_p[W][n]<<" out of range"<<std::endl;
          evalError = true;
        }
      else
        {
          local_p[N][n] = pc.psk.getPsiC() + local_p[W][n];
          
          pc.densityW.setHead(local_p[W][n]);
          pc.densityN.setHead(local_p[N][n]);
        }

      local_DthetaW_DpC[n] = pc.psk.getDthetaW_DpC();

      if (fabs(local_DthetaW_DpC[n]) > 0)
        local_DpC_DthetaW[n] = 1.0/local_DthetaW_DpC[n];
      else
        local_DpC_DthetaW[n] = 0.0;

      local_K[W][n] = pc.psk.getKrW();
      local_K[N][n] = pc.psk.getKrN();

      local_rho[W][n] = pc.densityW.getRho();
      local_rho[N][n] = pc.densityN.getRho();

      local_Drho[W][n] = pc.densityW.getDrho();
      local_Drho[N][n] = pc.densityN.getDrho();
      local_mCurrent[W][n] = local_theta[W][n] * local_rho[W][n];
      local_mCurrent[N][n] = local_theta[N][n] * local_rho[N][n];
    }      
  
  dkg->computeDivergence(local_K,local_rho,local_p);

  return evalError;
}

template <class PROB, class JAC>
bool TwopCDMMS<PROB,JAC>::calculateJacCoefficients()
{
  local_y.startSetFromGlobal(*myY);
  local_y.endSetFromGlobal(*myY);
  bool evalError=false;
  int n;
  int end=local_p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      local_theta[W][n] = local_y[4*n];
      local_p[W][n]=local_y[4*n+1];
    }
  for (n=0;n<end;n++)
    {
      local_theta[N][n] = (local_thetaS[n] - local_theta[W][n]);
      
      if (local_theta[N][n] < pc.mNegative)
        {
          std::cerr<<"1theta("<<n<<")="<<local_theta[N][n]<<" out of range"<<std::endl;
          evalError = true;
        }
      else if (local_theta[W][n] <= pc.psk.getLocalThetaR()[n])
        {
          std::cerr<<"2theta("<<n<<")="<<local_theta[W][n]<<" out of range"<<std::endl;
          evalError =  true;
        }
      else
        pc.psk.setVFraction(local_theta[W][n],n);
      
      if (fabs(local_p[W][n]) > pc.pBig)
        {
          std::cerr<<"pW("<<n<<")="<<local_p[W][n]<<" out of range"<<std::endl;
          evalError = true;
        }
      else
        {
          local_p[N][n] = pc.psk.getPsiC() + local_p[W][n];
          
          pc.densityW.setHead(local_p[W][n]);
          pc.densityN.setHead(local_p[N][n]);
        }

      local_DthetaW_DpC[n] = pc.psk.getDthetaW_DpC();

      if (fabs(local_DthetaW_DpC[n]) > 0)
        local_DpC_DthetaW[n] = 1.0/local_DthetaW_DpC[n];
      else
        local_DpC_DthetaW[n] = 0.0;

      local_K[W][n] = pc.psk.getKrW();
      local_K[N][n] = pc.psk.getKrN();

      local_rho[W][n] = pc.densityW.getRho();
      local_rho[N][n] = pc.densityN.getRho();

      local_Drho[W][n] = pc.densityW.getDrho();
      local_Drho[N][n] = pc.densityN.getDrho();


      pc.psk.calculateDerivativesVFraction();
      local_DDthetaW_DDpC[n] = pc.psk.getDDthetaW_DDpC();
      local_DDpC_DDthetaW[n] = -local_DpC_DthetaW[n]*local_DpC_DthetaW[n]
	*local_DDthetaW_DDpC[n]*local_DpC_DthetaW[n];

      local_DK[W][W][n] = -pc.psk.getDkrW_DpC();
      local_DK[W][N][n] = pc.psk.getDkrW_DpC();
      local_DK[N][W][n] = -pc.psk.getDkrN_DpC();
      local_DK[N][N][n] = pc.psk.getDkrN_DpC();

      local_DDrho[W][n] = pc.densityW.getDDrho();
      local_DDrho[N][n] = pc.densityN.getDDrho();
    }

  dkg->computeDivergenceJac(local_K,local_rho,local_p,
                            local_DK,local_Drho,D_Div);
  return evalError;
}

template <class PROB, class JAC> 
class TwopCDMPP : public TwopDaeDef
{
public:
  //mwf added to switch out DivKgrad types
  //upwind rel perms
  typedef DivKgradUG<typename PROB::BC,2>   DIVKGRAD;
  typedef DivKgradUG1d<typename PROB::BC,2> DIVKGRAD1d;
  typedef DivKgradUG2d<typename PROB::BC,2> DIVKGRAD2d;
  typedef DivKgradUG3d<typename PROB::BC,2> DIVKGRAD3d;
  //arithmetic avg rel perms
  //typedef DivKgrad<BC,2> DIVKGRAD;
  //typedef DivKgrad1d<BC,2> DIVKGRAD1d;
  //typedef DivKgrad2d<BC,2> DIVKGRAD2d;
  //typedef DivKgrad3d<BC,2> DIVKGRAD3d;

  enum variable {PW, PN};
  TwopCDMPP(Petsc::SecondOrderFd& s,ParameterDatabase& pd, DataCollector& dataIn);
  virtual ~TwopCDMPP() {delete dkg;}
  void setICvars();
  bool residual(const real& t,const Vec& y,const Vec& yp, Vec& res);
  bool residual(const VecIndex&,const real& t,const Vec& y,
                       const Vec& yp, Vec& res){return true;}  
  bool yPrimeValue(const real& t, const Vec& y, Vec& yp);
  bool evaluateDaeJacobian(const real &t,const  Vec& y,const  Vec& yp,
                                    const real& alphaBDF);
  //try to put in residual calculation at least
  void stepTaken()         {pc.psk.updateHistory();}

  bool jacVec(const Vec& x, Vec& Jx){jacMat->apply(x,Jx); return false;}
  void setJac(JAC& jacRep) {jacMat = &jacRep;}
  bool calculateCoefficients();
  bool calculateJacCoefficients();

  virtual const Vec& getFlux_x(int n) { return dkg->getFlux_x(n); }
  virtual const Vec& getFlux_y(int n) { return dkg->getFlux_y(n); }
  virtual const Vec& getFlux_z(int n) { return dkg->getFlux_z(n); }

  PROB pc;
  JAC *jacMat;

protected:
  DIVKGRAD* dkg;

};

template<class PROB, class JAC>
TwopCDMPP<PROB,JAC>::TwopCDMPP(Petsc::SecondOrderFd& s,ParameterDatabase& pd,
			       DataCollector& dataIn):
  TwopDaeDef(s,pd,dataIn,
	     2*pd.i("nxNodes")*pd.i("nyNodes")*pd.i("nzNodes")),
  pc(pd),
  jacMat(0),
  dkg(0)
{
  Tracer tr("TwopCDMPP");
  //setup dkg

  if (DIM == ONE_D)
    {
      dkg = new DIVKGRAD1d(pc.bc,node,nxNodes,pd.r("gx"),oneOverdx);
    }
  else if (DIM == TWO_D)
    {
      dkg = new DIVKGRAD2d(pc.bc,node,nxNodes,nyNodes,
			   pd.r("gx"),pd.r("gy"),
			   oneOverdx,oneOverdy);
    }
  else
    {
      dkg = new DIVKGRAD3d(pc.bc,node,nxNodes,nyNodes,nzNodes,
			   pd.r("gx"),pd.r("gy"),pd.r("gz"),
			   oneOverdx,oneOverdy,oneOverdz);
    }
  assert(dkg);
  dkg->computeInterfaceK(pc.psk.getLocalKWs());

  //If I want to switch approach out to pass in a problem
  //then I'll need to do this stuff outside as well
  //
  pc.setBounds(-10.0*pd.r("m_atol"));//make this more general?

  setICvars();
  pc.setupProblem(pd,*this,&alphaDaeDef);
  //also recalculate interface conductivities
  dkg->computeInterfaceK(pc.psk.getLocalKWs());

//    std::cout<<"thetaW"<<theta[W]<<std::endl;
#ifdef DEBUG
  psys.barrier();
  pc.bc[W].print();
  pc.bc[N].print();
  psys.barrier();
#endif
//    std::cout<<"calc coeff"<<std::endl;
  bool evalError=false;
  psys.barrier();
  evalError = psys.catchError(calculateCoefficients());
  psys.barrier();
  if (evalError)
    {
      std::cerr<<"initial values out of range"<<std::endl<<std::flush;
      exit(1);
    }
  Vec tmp(thetaS);
  tmp.setFromLocal(local_mCurrent[W]);  m[W] = tmp;
  tmp.setFromLocal(local_mCurrent[N]);  m[N] = tmp;
  psys.barrier();
//    std::cout<<"ypval"<<std::endl;
  yPrimeValue(t0,y0,y0prime);
  psys.barrier();
#ifdef DEBUG
  twopOut<<"y0"<<y0<<std::endl<<std::flush<<"y0'"
	 <<y0prime<<std::flush<<std::endl;
#endif
  betaDaeDef = y0prime;
  Vec temp(y0);
  //    std::cout<<"res cal"<<std::endl;
  psys.barrier();
  residual(t0,y0,y0prime,temp);
  psys.barrier();
#ifdef DEBUG  
  twopOut<<"y0"<<y0<<std::endl<<std::flush<<"y0'"
	 <<y0prime<<std::flush<<std::endl<<"res"<<temp<<std::flush;
#endif
  twopOut<<"Relative Residual of Initial Data: "
	 <<nrm2(temp)/nrm2(y0)<<std::endl;


}


template<class PROB, class JAC>
void TwopCDMPP<PROB,JAC>::setICvars()
{
  myY=&y0;
  p[W].attachToVecMulti(Vec::REF,y0,fwIndex);
  p[W].setStrideMulti(2);
  p[N].attachToVecMulti(Vec::REF,y0,fnIndex);
  p[N].setStrideMulti(2);
  Dp[W].attachToVecMulti(Vec::REF,y0prime,fwIndex);
  Dp[W].setStrideMulti(2);
  Dp[N].attachToVecMulti(Vec::REF,y0prime,fnIndex);
  Dp[N].setStrideMulti(2);
}



template <class PROB, class JAC>
bool TwopCDMPP<PROB,JAC>::residual(const real &t,const  Vec& y,const  Vec& yp, 
				Vec& res)
{
  myY=&y;
  //split vector into NAPL and water pressures
             
  p[W].attachToVecMulti(Vec::REF,y,fwIndex);
  Dp[W].attachToVecMulti(Vec::REF,yp,fwIndex);

  p[N].attachToVecMulti(Vec::REF,y,fnIndex);
  Dp[N].attachToVecMulti(Vec::REF,yp,fnIndex);

  p[W].setStrideMulti(2);
  Dp[W].setStrideMulti(2);

  p[N].setStrideMulti(2);
  Dp[N].setStrideMulti(2);

  //split res int NAPL and water mass balance residuals

  resF[W].attachToVecMulti(Vec::REF,res,fwIndex);
  resF[N].attachToVecMulti(Vec::REF,res,fnIndex);
  
  resF[W].setStrideMulti(2);
  resF[N].setStrideMulti(2);
  
  pc.adjustBC(t);

  bool evalError;
  evalError=psys.catchError(calculateCoefficients());
  if (evalError)
    return evalError;
  // compute res on entire domain
  int n;
  int end=p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      node.localNode(n);
      resF[W][n] = (local_Drho[W][node.center]*local_theta[W][node.center] 
		    - local_rho[W][node.center]*local_DthetaW_DpC[node.center])*Dp[W][n]
        + local_rho[W][node.center]*local_DthetaW_DpC[node.center]*Dp[N][n] 
        - dkg->getDiv(node,W);
      
      resF[N][n] = local_rho[N][node.center]*local_DthetaW_DpC[node.center]*Dp[W][n] 
        + (local_Drho[N][node.center]*local_theta[N][node.center] 
	   - local_rho[N][node.center]*local_DthetaW_DpC[node.center])*Dp[N][n] 
        - dkg->getDiv(node,N);
    }

  
  pc.bc[W].applyDirichletConditions(node,resF[W]);
  pc.bc[N].applyDirichletConditions(node,resF[N]);

  return false;
}

template <class PROB, class JAC>
bool TwopCDMPP<PROB,JAC>::yPrimeValue(const real &t,const  Vec& y,Vec& yp)
{  
  myY=&y;
  pc.adjustBC(t);
  real a(12345),b(12345),c(12345),d(12345),det(12345);
             
  p[W].attachToVecMulti(Vec::REF,y,fwIndex);
  p[N].attachToVecMulti(Vec::REF,y,fnIndex);
  
  p[W].setStrideMulti(2);
  p[N].setStrideMulti(2);

  Dp[W].attachToVecMulti(Vec::REF,yp,fwIndex);
  Dp[N].attachToVecMulti(Vec::REF,yp,fnIndex);
  
  Dp[W].setStrideMulti(2);
  Dp[N].setStrideMulti(2);

  bool evalError = psys.catchError(calculateCoefficients());
  if (evalError)
    return evalError;
  
  //    compute yp on interior and at left and right boundaries
  int n;
  int end=p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      node.localNode(n);
      //solve 2x2 system at each node
      a = local_Drho[W][node.center]*local_theta[W][node.center] - local_rho[W][node.center]*local_DthetaW_DpC[node.center];
      b = local_rho[W][node.center]*local_DthetaW_DpC[node.center];
      c = local_rho[N][node.center]*local_DthetaW_DpC[node.center];
      d = local_Drho[N][node.center]*local_theta[N][node.center] - local_rho[N][node.center]*local_DthetaW_DpC[node.center];

      det = a*d - b*c;

      if (fabs(det) < MACHINE_EPSILON)
        {
          twopOut<<"local system singular n="<<n<<" dkg->getDiv(node,W)="
		 <<dkg->getDiv(node,W)<<" dkg->getDiv(node,N)="
		 <<dkg->getDiv(node,N)<<std::endl;
          real Sw= (local_theta[W][node.center]/thetaS[n]),     
            Sn = (1.0 - local_theta[W][node.center]/thetaS[n]);
          if (local_DthetaW_DpC[node.center]!=0 && 
	      (local_Drho[W][node.center]==0 || 
	       local_Drho[N][node.center] ==0))
            {
              twopOut<<"system is incompressible and locally unsaturated"<<std::endl;
              Dp[W][n] = Sw*dkg->getDiv(node,W)/a;
              Dp[N][n] = Sn*dkg->getDiv(node,W)/b;
            }
          else
            {
              twopOut<<"system is locally saturated and incompressible"<<std::endl;
              Dp[W][n] = 0.0;
              Dp[N][n] = 0.0;
            }
        }
      else
        {
          Dp[W][n] = d*dkg->getDiv(node,W)/det - b*dkg->getDiv(node,N)/det;
          Dp[N][n] = a*dkg->getDiv(node,N)/det - c*dkg->getDiv(node,W)/det;
        }
    }

  pc.bc[W].applyDirichletYprime(node);
  pc.bc[N].applyDirichletYprime(node);

  for (n=0;n<end;n++)
    {
      if (Dp[W][n] == 0.0 && Dtheta[W][n] == 0.0)
        {
          Dp[N][n] = 0.0;
        }
      else if (Dp[N][n] == 0.0 && Dp[W][n] !=0.0)
        {
          Dp[W][n] = dkg->getDiv(node,W)/a;
        }
      else if (Dp[W][n] == 0.0 && Dp[N][n] !=0.0)
        {
          Dp[N][n] = dkg->getDiv(node,W)/b;
        }
    }
  return false;
}


template <class PROB, class JAC>
bool TwopCDMPP<PROB,JAC>::evaluateDaeJacobian(const real &t,const  Vec& y,
					   const  Vec& yp,const real& alphaBDF) 
{  
  jacMat->zeroAll();
  myY = &y;
  const int W=TwopData::W,N=TwopData::N;

  //split vector into NAPL and water pressures
  p[W].attachToVecMulti(Vec::REF,y,fwIndex);
  Dp[W].attachToVecMulti(Vec::REF,yp,fwIndex);

  p[N].attachToVecMulti(Vec::REF,y,fnIndex);
  Dp[N].attachToVecMulti(Vec::REF,yp,fnIndex);

  p[W].setStrideMulti(2);
  Dp[W].setStrideMulti(2);

  p[N].setStrideMulti(2);
  Dp[N].setStrideMulti(2);

  pc.adjustBC(t);

  bool evalError=false;
  evalError = psys.catchError(calculateJacCoefficients());
  if (evalError)
    return evalError;

  real ypJac[2][2];
  real DypJac[2][2];

  for (int i=node.local_z0;i<node.local_nzNodes+node.local_z0;i++)
    for (int j=node.local_y0;j<node.local_nyNodes+node.local_y0;j++)
      for (int k=node.local_x0;k<node.local_nxNodes+node.local_x0;k++)
        {
          node.localIndex(i-node.local_z0,j-node.local_y0,k-node.local_x0);
          int nl=node.center_noGhost,
            nlg=node.center;
          node.globalIndex(i,j,k);
          
          ypJac[W][W] = local_Drho[W][nlg]*local_theta[W][nlg] 
	    - local_rho[W][nlg]*local_DthetaW_DpC[nlg];
          ypJac[W][N] = local_rho[W][nlg]*local_DthetaW_DpC[nlg];                         
          ypJac[N][W] = local_rho[N][nlg]*local_DthetaW_DpC[nlg];                         
          ypJac[N][N] = local_Drho[N][nlg]*local_theta[N][nlg] 
	    - local_rho[N][nlg]*local_DthetaW_DpC[nlg];
          
          DypJac[W][W] = (local_DDrho[W][nlg]*local_theta[W][nlg] 
			  - local_Drho[W][nlg]*local_DthetaW_DpC[nlg]  
                          - local_Drho[W][nlg]*local_DthetaW_DpC[nlg] 
			  + local_rho[W][nlg]*local_DDthetaW_DDpC[nlg])*Dp[W][nl] 
            + (local_Drho[W][nlg]*local_DthetaW_DpC[nlg] 
	       - local_rho[W][nlg]*local_DDthetaW_DDpC[nlg])*Dp[N][nl];
          
          DypJac[W][N] = (local_Drho[W][nlg]*local_DthetaW_DpC[nlg] 
			  - local_rho[W][nlg]*local_DDthetaW_DDpC[nlg] )*Dp[W][nl]
            + local_rho[W][nlg]*local_DDthetaW_DDpC[nlg]*Dp[N][nl];
          
          DypJac[N][W] = (-local_rho[N][nlg]*local_DDthetaW_DDpC[nlg])*Dp[W][nl] 
            + (local_Drho[N][nlg]*local_DthetaW_DpC[nlg] 
	       + local_rho[N][nlg]*local_DDthetaW_DDpC[nlg])*Dp[N][nl];
          
          DypJac[N][N] = (local_Drho[N][nlg]*local_DthetaW_DpC[nlg] 
			  + local_rho[N][nlg]*local_DDthetaW_DDpC[nlg])*Dp[W][nl] 
            + (local_DDrho[N][nlg]*local_theta[N][nlg] 
	       - local_Drho[N][nlg]*local_DthetaW_DpC[nlg] 
	       - local_Drho[N][nlg]*local_DthetaW_DpC[nlg] 
	       - local_rho[N][nlg]*local_DDthetaW_DDpC[nlg])*Dp[N][nl];
          

          (*jacMat)(2*node.center+W,2*node.center+W)  = DypJac[W][W]
            -D_Div[W][W][DIVKGRAD::CENTER][nl]
            + alphaBDF*ypJac[W][W];
          
          (*jacMat)(2*node.center+N,2*node.center+W)  = DypJac[N][W]
            -D_Div[N][W][DIVKGRAD::CENTER][nl]
            + alphaBDF*ypJac[N][W];

          (*jacMat)(2*node.center+W,2*node.center+N)  = DypJac[W][N]
            -  D_Div[W][N][DIVKGRAD::CENTER][nl]
            + alphaBDF*ypJac[W][N];
          
          (*jacMat)(2*node.center+N,2*node.center+N)  = DypJac[N][N]
            - D_Div[N][N][DIVKGRAD::CENTER][nl]
            + alphaBDF*ypJac[N][N];
          
          if (node.anchor->k > 0)
	    {
	      (*jacMat)(2*node.center+W,2*node.left+W)= 
		-D_Div[W][W][DIVKGRAD::LEFT][nl];
	      (*jacMat)(2*node.center+N,2*node.left+W)= 
		-D_Div[N][W][DIVKGRAD::LEFT][nl];
	      (*jacMat)(2*node.center+W,2*node.left+N)= 
		-D_Div[W][N][DIVKGRAD::LEFT][nl];
	      (*jacMat)(2*node.center+N,2*node.left+N)= 
		-D_Div[N][N][DIVKGRAD::LEFT][nl];
          }
        if (node.anchor->k < nxNodes-1)
          {
            (*jacMat)(2*node.center+W,2*node.right+W)= 
	      -D_Div[W][W][DIVKGRAD::RIGHT][nl];
            (*jacMat)(2*node.center+N,2*node.right+W)= 
	      -D_Div[N][W][DIVKGRAD::RIGHT][nl];
            (*jacMat)(2*node.center+W,2*node.right+N)= 
	      -D_Div[W][N][DIVKGRAD::RIGHT][nl];
            (*jacMat)(2*node.center+N,2*node.right+N)= 
	      -D_Div[N][N][DIVKGRAD::RIGHT][nl];
          }
        
        if (DIM > ONE_D)
          {
            if (node.anchor->j > 0)
              {
                (*jacMat)(2*node.center+W,2*node.front+W)= 
		  -D_Div[W][W][DIVKGRAD::FRONT][nl];
                (*jacMat)(2*node.center+N,2*node.front+W)= 
		  -D_Div[N][W][DIVKGRAD::FRONT][nl];
                (*jacMat)(2*node.center+W,2*node.front+N)= 
		  -D_Div[W][N][DIVKGRAD::FRONT][nl];
                (*jacMat)(2*node.center+N,2*node.front+N)= 
		  -D_Div[N][N][DIVKGRAD::FRONT][nl];
              }
            if (node.anchor->j < nyNodes-1)
              {
                (*jacMat)(2*node.center+W,2*node.back+W)= 
		  -D_Div[W][W][DIVKGRAD::BACK][nl];
                (*jacMat)(2*node.center+N,2*node.back+W)= 
		  -D_Div[N][W][DIVKGRAD::BACK][nl];
                (*jacMat)(2*node.center+W,2*node.back+N)= 
		  -D_Div[W][N][DIVKGRAD::BACK][nl];
                (*jacMat)(2*node.center+N,2*node.back+N)= 
		  -D_Div[N][N][DIVKGRAD::BACK][nl];
              }
          }
        if (DIM == THREE_D)
          {
            if (node.anchor->i > 0)
              {
                (*jacMat)(2*node.center+W,2*node.bottom+W)= 
		  -D_Div[W][W][DIVKGRAD::BOTTOM][nl];
                (*jacMat)(2*node.center+N,2*node.bottom+W)= 
		  -D_Div[N][W][DIVKGRAD::BOTTOM][nl];
                (*jacMat)(2*node.center+W,2*node.bottom+N)= 
		  -D_Div[W][N][DIVKGRAD::BOTTOM][nl];
                (*jacMat)(2*node.center+N,2*node.bottom+N)= 
		  -D_Div[N][N][DIVKGRAD::BOTTOM][nl];
              }
            
            if (node.anchor->i < nzNodes-1)
              {
                (*jacMat)(2*node.center+W,2*node.top+W)= 
		  -D_Div[W][W][DIVKGRAD::TOP][nl];
                (*jacMat)(2*node.center+N,2*node.top+W)= 
		  -D_Div[N][W][DIVKGRAD::TOP][nl];
                (*jacMat)(2*node.center+W,2*node.top+N)= 
		  -D_Div[W][N][DIVKGRAD::TOP][nl];
                (*jacMat)(2*node.center+N,2*node.top+N)= 
		  -D_Div[N][N][DIVKGRAD::TOP][nl];
              } 
          }
        jacMat->finalizeBlockRow(node.center);
      }
  
  int dof=2;
  //will need to be fixed or take  a different approach when mixing bc variables.
  pc.bc[W].applyDirichletDerivatives(node,*jacMat,dof,W); 
  pc.bc[N].applyDirichletDerivatives(node,*jacMat,dof,N);

  //fix boundary conditions (set for theta_W and P_N instead of the primary variables for this form)
  //now fix theta_W BC
  pc.bc[W].dit = pc.bc[W].dirichlet.begin();
  while(pc.bc[W].dit != pc.bc[W].dirichlet.end()) 
    { 
      node.localNode(pc.bc[W].dit->n - thetaS.getGlobalLow());
      int nlg=node.center;
      node.globalNode(pc.bc[W].dit->n);
      (*jacMat)(2*node.center+W,2*node.center+W) = 
	-(*(pc.bc[W].dit->scale))*local_DthetaW_DpC[nlg];
      (*jacMat)(2*node.center+W,2*node.center+N) = 
	(*(pc.bc[W].dit->scale))*local_DthetaW_DpC[nlg];
      jacMat->finalizeRow(2*node.center+W);
      ++(pc.bc[W].dit);
    }  
  //now fix theta_W BC
  pc.bc[N].dit = pc.bc[N].dirichlet.begin();
  while(pc.bc[N].dit != pc.bc[N].dirichlet.end()) 
    { 
      node.globalNode(pc.bc[N].dit->n);
      (*jacMat)(2*node.center+N,2*node.center+W) = (*(pc.bc[N].dit->scale));
      (*jacMat)(2*node.center+N,2*node.center+N) = 0.0;
      jacMat->finalizeRow(2*node.center+N);
      ++(pc.bc[N].dit);
    }  
  jacMat->beginAssembly();
  jacMat->endAssembly();

  return false;
}


template <class PROB, class JAC>
bool TwopCDMPP<PROB,JAC>::calculateCoefficients()
{
  local_y.startSetFromGlobal(*myY);
  local_y.endSetFromGlobal(*myY);

  bool evalError=false;
  int n;
  int end=local_p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      local_p[W][n] = local_y[2*n];
      local_p[N][n]=local_y[2*n+1];
    }
  for (n=0;n<end;n++)
    { 
      if (fabs(local_p[W][n]) > pc.pBig || fabs(local_p[N][n]) > pc.pBig)
        {
          if (fabs(local_p[W][n]) > pc.pBig)
            std::cerr<<"local_pW("<<n<<")="<<local_p[W][n]<<" out of range"<<std::endl;
          else 
            std::cerr<<"local_pN("<<n<<")="<<local_p[N][n]<<" out of range"<<std::endl;
          evalError=true;
        }
      pc.densityW.setHead(local_p[W][n]);
      pc.densityN.setHead(local_p[N][n]);

      pc.psk.setHeads(n,local_p[W][n],local_p[N][n]);

      local_theta[W][n] = pc.psk.getThetaW();
      local_theta[N][n] = (local_thetaS[n] - local_theta[W][n]);

      local_DthetaW_DpC[n] = pc.psk.getDthetaW_DpC();

      local_K[W][n] = pc.psk.getKrW();
      local_K[N][n] = pc.psk.getKrN();

      local_rho[W][n] = pc.densityW.getRho();
      local_rho[N][n] = pc.densityN.getRho();

      local_Drho[W][n] = pc.densityW.getDrho();
      local_Drho[N][n] = pc.densityN.getDrho();
    }
  dkg->computeDivergence(local_K,local_rho,local_p);
  return evalError;
}

template <class PROB, class JAC>
bool TwopCDMPP<PROB,JAC>::calculateJacCoefficients()
{
  local_y.startSetFromGlobal(*myY);
  local_y.endSetFromGlobal(*myY);

  bool evalError=false;
  int n;
  int end=local_p[W].getLocalHigh();
  for (n=0;n<end;n++)
    {
      local_p[W][n] = local_y[2*n];
      local_p[N][n]=local_y[2*n+1];
    }
  for (n=0;n<end;n++)
    { 
      if (fabs(local_p[W][n]) > pc.pBig || fabs(local_p[N][n]) > pc.pBig)
        {
          if (fabs(local_p[W][n]) > pc.pBig)
            std::cerr<<"local_pW("<<n<<")="<<local_p[W][n]<<" out of range"<<std::endl;
          else 
            std::cerr<<"local_pN("<<n<<")="<<local_p[N][n]<<" out of range"<<std::endl;
          evalError= true;
        }
      pc.densityW.setHead(local_p[W][n]);
      pc.densityN.setHead(local_p[N][n]);

      pc.psk.setHeads(n,local_p[W][n],local_p[N][n]);

      local_theta[W][n] = pc.psk.getThetaW();
      local_theta[N][n] = (local_thetaS[n] - local_theta[W][n]);

      local_DthetaW_DpC[n] = pc.psk.getDthetaW_DpC();

      local_K[W][n] = pc.psk.getKrW();
      local_K[N][n] = pc.psk.getKrN();

      local_rho[W][n] = pc.densityW.getRho();
      local_rho[N][n] = pc.densityN.getRho();

      local_Drho[W][n] = pc.densityW.getDrho();
      local_Drho[N][n] = pc.densityN.getDrho();

      pc.psk.calculateDerivatives();
      local_DDthetaW_DDpC[n] = pc.psk.getDDthetaW_DDpC();

      //these are derivatives w.r.t. pw,pn
      local_DK[W][W][n] = -pc.psk.getDkrW_DpC();
      local_DK[W][N][n] =  pc.psk.getDkrW_DpC();
      local_DK[N][W][n] = -pc.psk.getDkrN_DpC();
      local_DK[N][N][n] =  pc.psk.getDkrN_DpC();

      local_DDrho[W][n] = pc.densityW.getDDrho();
      local_DDrho[N][n] = pc.densityN.getDDrho();
    }

  dkg->computeDivergenceJac(local_K,local_rho,local_p,
                            local_DK,local_Drho,D_Div);
  return evalError;
}


}//TwoPhaseFlow
}//Daetk
#endif
