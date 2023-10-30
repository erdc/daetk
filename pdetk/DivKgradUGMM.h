#ifndef DIVKGRADUGMM_H
#define DIVKGRADUGMM_H
#include "Definitions.h"
#include "PetscSecondOrderFd.h"
#include "Divergence.h"
#include <vector>

//  FIX DIVERGENCE JAC -- MUST APPLY CORRECT DERIVATIVES FOR NEUMAN BC'S

namespace Daetk 
{

using Petsc::SecondOrderFd;
// typedef CMRVec<Vec> VecVec;
// typedef CMRVec<VecVec> VecVecVec;
// typedef CMRVec<VecVecVec> VecVecVecVec;

template<class BC, int nv>
class DivKgradUGMM : public Divergence
{
public:
  enum JacEntry {CENTER, LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP};
  DivKgradUGMM(BC* bcIn, SecondOrderFd& nodeIn, real* vis_ratio=0);
  
  virtual ~DivKgradUGMM();

  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div)=0;

  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)=0;  
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P)=0;

  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho)=0;  
  
  virtual void computeInterfaceK(const Vec& K)=0;

  inline void setLocal_dy_dxi(Vec& local_dy_dxiIn){local_dy_dxi=&local_dy_dxiIn;}
  inline void setKsx(const Vec& K);
  inline void setKsy(const Vec& K);
  inline void setKsz(const Vec& K);

  inline real Kx(const Vec& K, const real& grad);
  inline real DKx_center(const Vec& DK, const real& grad);
  inline real DKx_right(const Vec& DK, const real& grad);

  inline real Ky(const Vec& K, const real& grad);
  inline real DKy_center(const Vec& DK, const real& grad);
  inline real DKy_back(const Vec& DK, const real& grad);

  inline real Kz(const Vec& K, const real& grad);
  inline real DKz_center(const Vec& DK, const real& grad);
  inline real DKz_top(const Vec& DK, const real& grad);

  inline real Rhox(const Vec& Rho);
  inline real DRhox_center(const Vec& DRho);
  inline real DRhox_right(const Vec& DRho);

  inline real Rhoy(const Vec& Rho);
  inline real DRhoy_center(const Vec& DRho);
  inline real DRhoy_back(const Vec& DRho);

  inline real Rhoz(const Vec& Rho);
  inline real DRhoz_center(const Vec& DRho);
  inline real DRhoz_top(const Vec& DRho);

  inline real DpDx(const Vec& P);
  inline real DpDy(const Vec& P);
  inline real DpDz(const Vec& P);

  inline void setFlux_x(int v,const Vec& K, const Vec& Rho, const Vec& P);

  inline void setFlux_x(const Vec* K, const Vec* Rho, const Vec* P);
  inline void setFlux_y(const Vec* K, const Vec* Rho, const Vec* P);
  inline void setFlux_z(const Vec* K, const Vec* Rho, const Vec* P);

  inline void setDFlux_x(const Vec* K, const Vec* Rho, const Vec* P, 
                         const VecVecVec& DK, const Vec* DRho);
  inline void setDFlux_y(const Vec* K, const Vec* Rho, const Vec* P, 
                         const VecVecVec& DK, const Vec* DRho);
  inline void setDFlux_z(const Vec* K, const Vec* Rho, const Vec* P, 
                         const VecVecVec& DK, const Vec* DRho);


  inline real getDiv(SecondOrderFd& n, int div_i=0);

  inline real getDivJacCenter(SecondOrderFd& n, int div_i=0, int u_j=0);
  inline real getDivJacLeft(SecondOrderFd& n, int div_i=0, int u_j=0);
  inline real getDivJacRight(SecondOrderFd& n, int div_i=0, int u_j=0);
  inline real getDivJacFront(SecondOrderFd& n, int div_i=0, int u_j=0);
  inline real getDivJacBack(SecondOrderFd& n, int div_i=0, int u_j=0);
  inline real getDivJacBottom(SecondOrderFd& n, int div_i=0, int u_j=0);
  inline real getDivJacTop(SecondOrderFd& n, int div_i=0, int u_j=0);

  const Vec& getFlux_x(int n=0);
  const Vec& getFlux_y(int n=0);
  const Vec& getFlux_z(int n=0);
  
  Vec& setFlux_x(int n=0);
  Vec& setFlux_y(int n=0);
  Vec& setFlux_z(int n=0);
  
  int i,j,k,vi,vj,nxNodes,nyNodes,nzNodes,local_nxNodes,local_nyNodes,local_nzNodes;
  real oneOverdx,oneOverdy,oneOverdz,gx,gy,gz; 
  real viscosity_ratio[nv];
  Vec* local_dy_dxi;
  Vec flux_x[nv],flux_y[nv],flux_z[nv],Ksx,Ksy,Ksz;
  VecVecVec Dflux_x_center,Dflux_x_right,Dflux_y_center,Dflux_y_back,Dflux_z_center,Dflux_z_top;
  BC* bc;
  SecondOrderFd& node;
};

template<class BC, int nv>
class DivKgrad1dUGMM : public DivKgradUGMM<BC,nv>
{
public:
  DivKgrad1dUGMM(BC* bcIn, 
             SecondOrderFd& nodeIn,
             int nNodes, 
             real g, 
             real oneOverd, 
             real* vis_ratio=0);

  virtual ~DivKgrad1dUGMM();
  
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div);

  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac);  

  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P);
  
  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho);
  
  virtual void computeInterfaceK(const Vec& K);
}; 

template<class BC, int nv>
class DivKgrad2dUGMM : public DivKgradUGMM<BC,nv>
{
public:
  DivKgrad2dUGMM(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn, 
             real gxIn, real gyIn,
             real oneOverdxIn, real oneOverdyIn, real* vis_ratio=0);
  virtual ~DivKgrad2dUGMM();
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div);

  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac);  
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P);
  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho);
   virtual void computeInterfaceK(const Vec& K);
};  

template<class BC, int nv>
class DivKgrad3dUGMM : public DivKgradUGMM<BC,nv>
{
public:
  DivKgrad3dUGMM(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn,int nzNodesIn, 
             real gxIn, real gyIn, real gzIn,
             real oneOverdxIn, real oneOverdyIn, real oneOverdzIn,
             real* vis_ratio=0);
  virtual ~DivKgrad3dUGMM();
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div);

  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac);  
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P);
  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho);
  virtual void computeInterfaceK(const Vec& K);
};  
  
template<class BC, int nv>
DivKgradUGMM<BC,nv>::DivKgradUGMM(BC* bcIn, SecondOrderFd& nodeIn, real* vis_ratio):
  i(0),
  j(0),
  k(0),
  vi(0),
  vj(0),
  nxNodes(1),
  nyNodes(1),
  nzNodes(1),
  local_nxNodes(1),
  local_nyNodes(1),
  local_nzNodes(1),
  bc(bcIn),
  node(nodeIn)
{
  if (vis_ratio)
    {
      for (int i=0;i<nv;i++)
        viscosity_ratio[i] = vis_ratio[i];
    }
  else
    for (int i=0;i<nv;i++)
      viscosity_ratio[i]=1.0;
  
  Dflux_x_center.resize(nv);
  Dflux_x_right.resize(nv);
  Dflux_y_center.resize(nv);
  Dflux_y_back.resize(nv);
  Dflux_z_center.resize(nv);
  Dflux_z_top.resize(nv);
  Ksx.newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
  Ksy.newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
  Ksz.newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
  Ksx = 0.0;
  Ksy = 0.0;
  Ksz = 0.0;
  for (vi=0;vi<nv;vi++)
    {
      flux_x[vi].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
      flux_y[vi].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
      flux_z[vi].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
      flux_x[vi] = 0.0;
      flux_y[vi] = 0.0;
      flux_z[vi] = 0.0;
      Dflux_x_center[vi].resize(nv);
      Dflux_x_right[vi].resize(nv);
      Dflux_y_center[vi].resize(nv);
      Dflux_y_back[vi].resize(nv);
      Dflux_z_center[vi].resize(nv);
      Dflux_z_top[vi].resize(nv);
      for (vj=0;vj<nv;vj++)
	{
	  Dflux_x_center[vi][vj].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
	  Dflux_x_center[vi][vj] = 0.0;
	  Dflux_x_right[vi][vj].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
	  Dflux_x_right[vi][vj] = 0.0;
	  
	  Dflux_y_center[vi][vj].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
	  Dflux_y_center[vi][vj] = 0.0;
	  Dflux_y_back[vi][vj].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
	  Dflux_y_back[vi][vj] = 0.0;
	  
	  Dflux_z_center[vi][vj].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
	  Dflux_z_center[vi][vj] = 0.0;
	  Dflux_z_top[vi][vj].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
	  Dflux_z_top[vi][vj] = 0.0;
	}
    }
}

template<class BC, int nv>
DivKgradUGMM<BC,nv>::~DivKgradUGMM(){}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setKsx(const Vec& K)
{
  Ksx[node.interRight] = 2.0 / ( 1.0/K[node.center] + 1.0/K[node.right]);
}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setKsy(const Vec& K)
{
  Ksy[node.interBack] = 2.0 / ( 1.0/K[node.center] + 1.0/K[node.back]);
}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setKsz(const Vec& K)
{
  Ksz[node.interTop] = 2.0 / ( 1.0/K[node.center] + 1.0/K[node.top]);
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::Kx(const Vec& K, const real& grad )
{ 
  if (grad > 0.0)
    return K[node.right]*Ksx[node.interRight];
  return K[node.center]*Ksx[node.interRight];
}


template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DKx_center(const Vec& DK, const real& grad)
{
  if (grad > 0.0)
    return 0.0;
  return DK[node.center]*Ksx[node.interRight];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DKx_right(const Vec& DK, const real& grad)
{
  if (grad > 0.0)
    return DK[node.right]*Ksx[node.interRight];
  return 0.0;
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::Ky(const Vec& K, const real& grad)
{ 
  if (grad > 0.0)
    return K[node.back]*Ksy[node.interBack];
  return K[node.center]*Ksy[node.interBack];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DKy_center(const Vec& DK,const real& grad)
{ 
  if (grad > 0.0)
    return 0.0;
  return DK[node.center]*Ksy[node.interBack];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DKy_back(const Vec& DK, const real& grad)
{ 
  if (grad > 0.0)
    return DK[node.back]*Ksy[node.interBack];
  return 0.0;
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::Kz(const Vec& K, const real& grad)
{ 
  if (grad > 0.0)
    return K[node.top]*Ksz[node.interTop];
  return K[node.center]*Ksz[node.interTop];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DKz_center(const Vec& DK, const real& grad)
{ 
  if (grad > 0.0)
    return 0.0;
  return DK[node.center]*Ksz[node.interTop];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DKz_top(const Vec& DK, const real& grad)
{ 
  if (grad > 0.0)
    return DK[node.top]*Ksz[node.interTop];
  return 0.0;
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::Rhox(const Vec& Rho)
{ 
  return 0.5*(Rho[node.center]+Rho[node.right]);
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DRhox_center(const Vec& DRho)
{ 
  return 0.5*DRho[node.center];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DRhox_right(const Vec& DRho)
{ 
  return 0.5*DRho[node.right];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::Rhoy(const Vec& Rho)
{ 
  return 0.5*(Rho[node.center]+Rho[node.back]);
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DRhoy_center(const Vec& DRho)
{ 
  return 0.5*DRho[node.center];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DRhoy_back(const Vec& DRho)
{ 
  return 0.5*DRho[node.back];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::Rhoz(const Vec& Rho)
{ 
  return 0.5*(Rho[node.center]+Rho[node.top]);
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DRhoz_center(const Vec& DRho)
{ 
  return 0.5*DRho[node.center];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DRhoz_top(const Vec& DRho)
{ 
  return 0.5*DRho[node.top];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DpDx(const Vec& P)
{ 
  return (oneOverdx/((*local_dy_dxi)[node.interRight]))*(P[node.right] - P[node.center]);
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DpDy(const Vec& P)
{ 
  return oneOverdy*(P[node.back]  - P[node.center]);
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::DpDz(const Vec& P)
{ 
  return oneOverdz*(P[node.top]   - P[node.center]);
}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setFlux_x(int v,const Vec& K, const Vec& Rho, const Vec& P)
{
  real rhox=Rhox(Rho),
    grad=DpDx(P) + rhox*gx;
  flux_x[v][node.interRight] = -rhox*viscosity_ratio[v]*
    Kx(K,grad) * grad;
}


template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setFlux_x(const Vec* K, const Vec* Rho, const Vec* P)
{	
  for (vi=0;vi<nv;vi++)
    {
      real rhox=Rhox(Rho[vi]),
	grad=DpDx(P[vi]) + rhox*gx;
      flux_x[vi][node.interRight] = -rhox*viscosity_ratio[vi]*
        Kx(K[vi],grad) * grad;
    }
}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setFlux_y(const Vec* K, const Vec* Rho, const Vec* P)
{
  for (vi=0;vi<nv;vi++)
    {
      real grad=DpDy(P[vi]) + Rhoy(Rho[vi])*gy;
      flux_y[vi][node.interBack]  = -Rhoy(Rho[vi])*viscosity_ratio[vi]*
        Ky(K[vi],grad) * grad;
    }
}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setFlux_z(const Vec* K, const Vec* Rho, const Vec* P)
{
  for (vi=0;vi<nv;vi++)
    {
      real grad=DpDz(P[vi]) + Rhoz(Rho[vi])*gz;
      flux_z[vi][node.interTop]	  = -Rhoz(Rho[vi])*viscosity_ratio[vi]*
        Kz(K[vi],grad) * grad;
    }
}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setDFlux_x(const Vec* K, const Vec* Rho, const Vec* P, 
                                        const VecVecVec& DK, const Vec* DRho)
{
//   //it's very important that the compiler be able to unroll these
//   for (vi=0;vi<nv;vi++)
//     {
//       for (vj=0;vj<nv;vj++)
//         {
//           //interior
//           real grad=DpDx(P[vi]) + Rhox(Rho[vi])*gx;
//           Dflux_x_center[vi][vj][node.interRight] = -Rhox(Rho[vi])*
//             DKx_center(DK[vi][vj],grad)*grad;
          
//           Dflux_x_right[vi][vj][node.interRight] = -Rhox(Rho[vi])*
//             DKx_right(DK[vi][vj],grad) * grad;
          
//         }
//     }
  
//   for (vi=0;vi<nv;vi++)
//     {
//       real grad=DpDx(P[vi]) + Rhox(Rho[vi])*gx;
//       Dflux_x_center[vi][vi][node.interRight] += 
//         -Rhox(Rho[vi])*Kx(K[vi],grad)*
//         (-oneOverdx + DRhox_center(DRho[vi])*gx) 
//         -DRhox_center(DRho[vi])*Kx(K[vi],grad)*grad;
      
//       Dflux_x_right[vi][vi][node.interRight] +=
//         -Rhox(Rho[vi])*Kx(K[vi],grad)*
//         (oneOverdx + DRhox_right(DRho[vi])*gx) 
//         -DRhox_right(DRho[vi])*Kx(K[vi],grad)*grad;
//     }

  //it's very important that the compiler be able to unroll these
  for (vi=0;vi<nv;vi++)
    {
      real rhox=Rhox(Rho[vi]),
	grad=DpDx(P[vi]) + rhox*gx,
	kx=viscosity_ratio[vi]*Kx(K[vi],grad),
	kxGrad=kx*grad,
	mrhoxKx=-rhox*kx,
	mrhoxGrad=-rhox*grad,
	drhox_center=DRhox_center(DRho[vi]),
	drhox_right=DRhox_right(DRho[vi]);


      
      for (vj=0;vj<nv;vj++)
	{
          //interior
          Dflux_x_center[vi][vj][node.interRight] = mrhoxGrad*
            DKx_center(DK[vi][vj],grad);
          
          Dflux_x_right[vi][vj][node.interRight] = mrhoxGrad*
            DKx_right(DK[vi][vj],grad);
          
	}

      Dflux_x_center[vi][vi][node.interRight] += 
        mrhoxKx*
        (-(oneOverdx/((*local_dy_dxi)[node.interRight])) + drhox_center*gx) 
        -drhox_center*kxGrad;
      
      Dflux_x_right[vi][vi][node.interRight] +=
        mrhoxKx*
        ((oneOverdx/((*local_dy_dxi)[node.interRight])) + drhox_right*gx) 
        -drhox_right*kxGrad;
    }
}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setDFlux_y(const Vec* K, const Vec* Rho, const Vec* P, 
                                        const VecVecVec& DK, const Vec* DRho)
{
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          real grad=DpDy(P[vi]) + Rhoy(Rho[vi])*gy;
          Dflux_y_center[vi][vj][node.interBack] = -Rhoy(Rho[vi])*
            viscosity_ratio[vi]*DKy_center(DK[vi][vj],grad)*grad;
          
          Dflux_y_back[vi][vj][node.interBack] = -Rhoy(Rho[vi])*
            viscosity_ratio[vi]*DKy_back(DK[vi][vj],grad) * grad;
          
        }
    }
  
  for (vi=0;vi<nv;vi++)
    {
      real grad=DpDy(P[vi]) + Rhoy(Rho[vi])*gy;
      Dflux_y_center[vi][vi][node.interBack] += 
        -Rhoy(Rho[vi])*viscosity_ratio[vi]*Ky(K[vi],grad)*
        (-oneOverdy + DRhoy_center(DRho[vi])*gy) 
        -DRhoy_center(DRho[vi])*viscosity_ratio[vi]*Ky(K[vi],grad)*grad;
      
      Dflux_y_back[vi][vi][node.interBack] +=
        -Rhoy(Rho[vi])*viscosity_ratio[vi]*Ky(K[vi],grad)*
        (oneOverdy + DRhoy_back(DRho[vi])*gy) 
        -DRhoy_back(DRho[vi])*viscosity_ratio[vi]*Ky(K[vi],grad)*grad;
    }
}

template<class BC, int nv>
inline void DivKgradUGMM<BC,nv>::setDFlux_z(const Vec* K, const Vec* Rho, const Vec* P, 
                                        const VecVecVec& DK, const Vec* DRho)
{
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          real grad=DpDz(P[vi]) + Rhoz(Rho[vi])*gz;
          Dflux_z_center[vi][vj][node.interTop] = -Rhoz(Rho[vi])*
            viscosity_ratio[vi]*DKz_center(DK[vi][vj],grad)*grad;
          
          Dflux_z_top[vi][vj][node.interTop] = -Rhoz(Rho[vi])*
            viscosity_ratio[vi]*DKz_top(DK[vi][vj],grad) *grad;
          
        }
    }
  
  for (vi=0;vi<nv;vi++)
    { 
      real grad=DpDz(P[vi]) + Rhoz(Rho[vi])*gz;
      Dflux_z_center[vi][vi][node.interTop] += 
        -Rhoz(Rho[vi])*viscosity_ratio[vi]*Kz(K[vi],grad)*
        (-oneOverdz + DRhoz_center(DRho[vi])*gz) 
        -DRhoz_center(DRho[vi])*viscosity_ratio[vi]*Kz(K[vi],grad)*grad;
      
      Dflux_z_top[vi][vi][node.interTop] +=
        -Rhoz(Rho[vi])*viscosity_ratio[vi]*Kz(K[vi],grad)*
	  (oneOverdz + DRhoz_top(DRho[vi])*gz) 
        -DRhoz_top(DRho[vi])*viscosity_ratio[vi]*Kz(K[vi],grad)*grad;
    }
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::getDiv(SecondOrderFd& n, int div_i)
{
  return oneOverdx*(flux_x[div_i][node.interLeft] - 
                    flux_x[div_i][node.interRight])/
            (0.5*( (*local_dy_dxi)[node.interLeft]) + (*local_dy_dxi)[node.interRight])
          + oneOverdy*(flux_y[div_i][node.interFront] - 
                       flux_y[div_i][node.interBack])
          + oneOverdz*(flux_z[div_i][node.interBottom] - 
                       flux_z[div_i][node.interTop]);
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::getDivJacCenter(SecondOrderFd& n, int div_i, int u_j)
{
  return oneOverdx*
    (Dflux_x_right[div_i][u_j][node.interLeft]
     -Dflux_x_center[div_i][u_j][node.interRight])
    + oneOverdy*
    (Dflux_y_back[div_i][u_j][node.interFront]
     -Dflux_y_center[div_i][u_j][node.interBack])
    + oneOverdz*
    (Dflux_z_top[div_i][u_j][node.interBottom]
     -Dflux_z_center[div_i][u_j][node.interTop]); 
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::getDivJacLeft(SecondOrderFd& n, int div_i, int u_j)
{
  return oneOverdx*
    Dflux_x_center[div_i][u_j][node.interLeft];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::getDivJacRight(SecondOrderFd& n, int div_i, int u_j)
{
  return -oneOverdx*
    Dflux_x_right[div_i][u_j][node.interRight];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::getDivJacFront(SecondOrderFd& n, int div_i, int u_j)
{
  return oneOverdy*
    Dflux_y_center[div_i][u_j][node.interFront];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::getDivJacBack(SecondOrderFd& n, int div_i, int u_j)
{
  return -oneOverdy*
    Dflux_y_back[div_i][u_j][node.interBack];
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::getDivJacBottom(SecondOrderFd& n, int div_i, int u_j)
{
  return oneOverdz*
    Dflux_z_center[div_i][u_j][node.interBottom]; 
}

template<class BC, int nv>
inline real DivKgradUGMM<BC,nv>::getDivJacTop(SecondOrderFd& n, int div_i, int u_j)
{
  return -oneOverdz*
    Dflux_z_top[div_i][u_j][node.interTop];
}


template<class BC, int nv>
const Vec& DivKgradUGMM<BC,nv>::getFlux_x(int n){return flux_x[n];}

template<class BC, int nv>
const Vec& DivKgradUGMM<BC,nv>::getFlux_y(int n){return flux_y[n];}

template<class BC, int nv>
const Vec& DivKgradUGMM<BC,nv>::getFlux_z(int n){return flux_z[n];}

template<class BC, int nv>
Vec& DivKgradUGMM<BC,nv>::setFlux_x(int n){return flux_x[n];}

template<class BC, int nv>
Vec& DivKgradUGMM<BC,nv>::setFlux_y(int n){return flux_y[n];}

template<class BC, int nv>
Vec& DivKgradUGMM<BC,nv>::setFlux_z(int n){return flux_z[n];}

template<class BC, int nv>
DivKgrad1dUGMM<BC,nv>::DivKgrad1dUGMM(BC* bcIn, 
             SecondOrderFd& nodeIn,
             int nNodes, 
             real g, 
             real oneOverd, real* vis_ratio):
   DivKgradUGMM<BC,nv>(bcIn,nodeIn,vis_ratio)
  {
    nxNodes=nNodes;
    local_nxNodes = nodeIn.local_nxNodes;

    gx=g;
    oneOverdx=oneOverd;

    Dflux_x_center.resize(nv);
    Dflux_x_right.resize(nv);
    Ksx.newsize(Vec::LOCAL,local_nxNodes+1);
    Ksx = 0.0;
    for (vi=0;vi<nv;vi++)
      {
        flux_x[vi].newsize(Vec::LOCAL,local_nxNodes+1);
        flux_x[vi] = 0.0;
        Dflux_x_center[vi].resize(nv);
        Dflux_x_right[vi].resize(nv);
        for (vj=0;vj<nv;vj++)
          {
            Dflux_x_center[vi][vj].newsize(Vec::LOCAL,local_nxNodes+1);
            Dflux_x_center[vi][vj] = 0.0;
            Dflux_x_right[vi][vj].newsize(Vec::LOCAL,local_nxNodes+1);
            Dflux_x_right[vi][vj] = 0.0;
          }
      }
  }

template<class BC, int nv>
DivKgrad1dUGMM<BC,nv>::~DivKgrad1dUGMM(){}

template<class BC, int nv>
void DivKgrad1dUGMM<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P)
{  
  for (k= -node.ghost_xOffSet; k< (node.ghost_nxNodes-node.ghost_xOffSet-1); k++)
    {
      node.localIndex(k);
      setFlux_x(K,Rho,P);
    }
  
  for (vi=0; vi<nv; vi++)
    bc[vi].applyNeumannConditions(node,flux_x[vi]);
}  

template<class BC, int nv>
void DivKgrad1dUGMM<BC,nv>::computeInterfaceK(const Vec& K)
{
  for (k= -node.ghost_xOffSet; k< (node.ghost_nxNodes-node.ghost_xOffSet-1); k++)
    {
      node.localIndex(k);
      setKsx(K);
    }
}  

template<class BC, int nv>
void DivKgrad1dUGMM<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, Vec* Div)
{
  computeDivergence(K,Rho,P);
  for (k=0; k<local_nxNodes; k++)
    {
      node.localIndex(k);
      for (vi=0;vi<nv;vi++)
        {
          Div[vi][node.center_noGhost] = oneOverdx*
            (flux_x[vi][node.interLeft] - flux_x[vi][node.interRight])/
            (0.5*( (*local_dy_dxi)[node.interLeft]) + (*local_dy_dxi)[node.interRight]);
        }
    }
}  

template<class BC, int nv>
void DivKgrad1dUGMM<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                             const VecVecVec& DK, const Vec* DRho)
{
  for (k=-node.ghost_xOffSet; k< node.ghost_nxNodes-node.ghost_xOffSet-1; k++)
    {
      node.localIndex(k);
      setDFlux_x(K,Rho,P,DK,DRho);
    }
}  

template<class BC, int nv>
void DivKgrad1dUGMM<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                             const VecVecVec& DK, const Vec* DRho, 
                                             VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Rho,P,DK,DRho);
  for (k=0;k<local_nxNodes;k++)
    {
      node.localIndex(k);
      for (vi=0;vi<nv;vi++)
        {
          for (vj=0;vj<nv;vj++)
            {
              DivJac[vi][vj][LEFT][node.center_noGhost] = oneOverdx*
                Dflux_x_center[vi][vj][node.interLeft]/
            (0.5*( (*local_dy_dxi)[node.interLeft]) + (*local_dy_dxi)[node.interRight]);
              
              DivJac[vi][vj][CENTER][node.center_noGhost] = oneOverdx*
                (Dflux_x_right[vi][vj][node.interLeft] 
                 -Dflux_x_center[vi][vj][node.interRight])/
            (0.5*( (*local_dy_dxi)[node.interLeft]) + (*local_dy_dxi)[node.interRight]); 
              
              DivJac[vi][vj][RIGHT][node.center_noGhost]  = -oneOverdx*
                Dflux_x_right[vi][vj][node.interRight]/
            (0.5*( (*local_dy_dxi)[node.interLeft]) + (*local_dy_dxi)[node.interRight]); 
            }
        }
    }
  
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          bc[vi].applyNeumannDerivatives(node,DivJac[vi][vj],
                                         Dflux_x_center[vi][vj],Dflux_x_right[vi][vj],
                                         oneOverdx);
        }
    }
}  

template<class BC, int nv>
DivKgrad2dUGMM<BC,nv>::DivKgrad2dUGMM(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn, 
             real gxIn, real gyIn,
             real oneOverdxIn, real oneOverdyIn, real* vis_ratio):
    DivKgradUGMM<BC,nv>(bcIn,nodeIn,vis_ratio)
  {
    nxNodes=nxNodesIn;
    local_nxNodes = nodeIn.local_nxNodes;
    gx=gxIn;
    oneOverdx = oneOverdxIn;

    nyNodes=nyNodesIn;
    local_nyNodes = nodeIn.local_nyNodes;
    gy=gyIn;
    oneOverdy = oneOverdyIn;

    Dflux_x_center.resize(nv);
    Dflux_x_right.resize(nv);
    Dflux_y_center.resize(nv);
    Dflux_y_back.resize(nv);
    Ksx.newsize(Vec::LOCAL,local_nyNodes*(local_nxNodes+1));
    Ksy.newsize(Vec::LOCAL,(local_nyNodes+1)*local_nxNodes);
    Ksx = 0.0;
    Ksy = 0.0;
    for (vi=0;vi<nv;vi++)
      {
          flux_x[vi].newsize(Vec::LOCAL,local_nyNodes*(local_nxNodes+1));
          flux_y[vi].newsize(Vec::LOCAL,(local_nyNodes+1)*local_nxNodes);
          flux_x[vi] = 0.0;
          flux_y[vi] = 0.0;
          Dflux_x_center[vi].resize(nv);
          Dflux_x_right[vi].resize(nv);
          Dflux_y_center[vi].resize(nv);
          Dflux_y_back[vi].resize(nv);
          for (vj=0;vj<nv;vj++)
            {
              Dflux_x_center[vi][vj].newsize(Vec::LOCAL,local_nyNodes*(local_nxNodes+1));
              Dflux_x_center[vi][vj] = 0.0;
              Dflux_x_right[vi][vj].newsize(Vec::LOCAL,local_nyNodes*(local_nxNodes+1));
              Dflux_x_right[vi][vj] = 0.0;
              Dflux_y_center[vi][vj].newsize(Vec::LOCAL,(local_nyNodes+1)*local_nxNodes);
              Dflux_y_center[vi][vj] = 0.0;
              Dflux_y_back[vi][vj].newsize(Vec::LOCAL,(local_nyNodes+1)*local_nxNodes);
              Dflux_y_back[vi][vj] = 0.0;
            }
        }
    }

template<class BC, int nv>
DivKgrad2dUGMM<BC,nv>::~DivKgrad2dUGMM(){}

template<class BC, int nv>
void DivKgrad2dUGMM<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P)
{
  //bottom ghost line
  if (node.ghost_yOffSet)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(-node.ghost_yOffSet,k);
        setFlux_y(K,Rho,P);
      }
  //top ghost line
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(local_nyNodes-1,k);
        setFlux_y(K,Rho,P);
      }
  
  //left ghost line
  if (node.ghost_xOffSet)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,-node.ghost_xOffSet);
        setFlux_x(K,Rho,P);
      }
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,local_nxNodes-1);
        setFlux_x(K,Rho,P);
      }
  
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(j,k);
          setFlux_x(K,Rho,P);
          setFlux_y(K,Rho,P);
        }
      //last y flux in row
      node.localIndex(j,k);
      setFlux_y(K,Rho,P);
    }
  //last row of x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setFlux_x(K,Rho,P);
    }
  
  for (vi=0;vi<nv;vi++)
    bc[vi].applyNeumannConditions(node,flux_x[vi], flux_y[vi]);
}

template<class BC, int nv>
void DivKgrad2dUGMM<BC,nv>::computeInterfaceK(const Vec& K)
{
  //bottom ghost line
  if (node.ghost_yOffSet)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(-node.ghost_yOffSet,k);
        setKsy(K);
      }
  //top ghost line
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(local_nyNodes-1,k);
        setKsy(K);
      }
  
  //left ghost line
  if (node.ghost_xOffSet)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,-node.ghost_xOffSet);
        setKsx(K);
      }
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,local_nxNodes-1);
        setKsx(K);
      }
  
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(j,k);
          setKsx(K);
          setKsy(K);
        }
      //last y flux in row
      node.localIndex(j,k);
      setKsy(K);
    }
  //last row of x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setKsx(K);
    }
}

template<class BC, int nv>
void DivKgrad2dUGMM<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, Vec* Div)
{
  computeDivergence(K,Rho,P);
  for (j=0;j<local_nyNodes;j++)
    for (k=0;k<local_nxNodes;k++)
      {
        node.localIndex(j,k);
        for (vi=0;vi<nv;vi++)
          {
            Div[vi][node.center_noGhost] = oneOverdx*
              (flux_x[vi][node.interLeft] - flux_x[vi][node.interRight])
              + oneOverdy*
              (flux_y[vi][node.interFront] - flux_y[vi][node.interBack]);
          }
      }
}

template<class BC, int nv>
void DivKgrad2dUGMM<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                             const VecVecVec& DK, const Vec* DRho)
{
  //bottom ghost line
  if (node.ghost_yOffSet)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(-node.ghost_yOffSet,k);
        setDFlux_y(K,Rho,P,DK,DRho);
      }
  //top ghost line
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(local_nyNodes-1,k);
        setDFlux_y(K,Rho,P,DK,DRho);
      }
  
  //left ghost line
  if (node.ghost_xOffSet)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,-node.ghost_xOffSet);
        setDFlux_x(K,Rho,P,DK,DRho);
      }
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,local_nxNodes-1);
        setDFlux_x(K,Rho,P,DK,DRho);
      }
  
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(j,k);
          setDFlux_x(K,Rho,P,DK,DRho);
          setDFlux_y(K,Rho,P,DK,DRho);
        }
      //last y flux in row
      node.localIndex(j,k);
      setDFlux_y(K,Rho,P,DK,DRho);
    }
  //last row of x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setDFlux_x(K,Rho,P,DK,DRho);
    }
}

template<class BC, int nv>
void DivKgrad2dUGMM<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                             const VecVecVec& DK, const Vec* DRho, 
                                             VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Rho,P,DK,DRho);
  for (j=0;j<local_nyNodes;j++)
    for (k=0;k<local_nxNodes;k++)
      {
        node.localIndex(j,k);
        for (vi=0;vi<nv;vi++)
          {
            for (vj=0;vj<nv;vj++)
              {
                
                DivJac[vi][vj][FRONT][node.center_noGhost] = oneOverdy*
                  Dflux_y_center[vi][vj][node.interFront];
                
                DivJac[vi][vj][LEFT][node.center_noGhost] = oneOverdx*
                  Dflux_x_center[vi][vj][node.interLeft];
                
                DivJac[vi][vj][CENTER][node.center_noGhost] = oneOverdx*
                  (Dflux_x_right[vi][vj][node.interLeft] 
                   -Dflux_x_center[vi][vj][node.interRight]) 
                  + oneOverdy*
                  (Dflux_y_back[vi][vj][node.interFront] 
                   -Dflux_y_center[vi][vj][node.interBack]); 
                
                DivJac[vi][vj][RIGHT][node.center_noGhost]  = -oneOverdx*
                  Dflux_x_right[vi][vj][node.interRight]; 
                
                DivJac[vi][vj][BACK][node.center_noGhost]  = -oneOverdy*
                  Dflux_y_back[vi][vj][node.interBack]; 
                
              }
          }
      }
  
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          bc[vi].applyNeumannDerivatives(node,DivJac[vi][vj],
                                         Dflux_x_center[vi][vj],Dflux_x_right[vi][vj],
                                         Dflux_y_center[vi][vj],Dflux_y_back[vi][vj],
                                         oneOverdx,oneOverdy);
          
        }
    }
}


template<class BC, int nv>
DivKgrad3dUGMM<BC,nv>::DivKgrad3dUGMM(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn,int nzNodesIn, 
             real gxIn, real gyIn, real gzIn,
             real oneOverdxIn, real oneOverdyIn, real oneOverdzIn, real* vis_ratio):
    DivKgradUGMM<BC,nv>(bcIn,nodeIn,vis_ratio)
    {
      nxNodes=nxNodesIn;
      local_nxNodes=nodeIn.local_nxNodes;
      gx=gxIn;
      oneOverdx = oneOverdxIn;

      nyNodes=nyNodesIn;
      local_nyNodes=nodeIn.local_nyNodes;
      gy=gyIn;
      oneOverdy = oneOverdyIn;

      nzNodes=nzNodesIn;
      local_nzNodes=nodeIn.local_nzNodes;
      gz=gzIn;
      oneOverdz = oneOverdzIn;

      Dflux_x_center.resize(nv);
      Dflux_x_right.resize(nv);
      Dflux_y_center.resize(nv);
      Dflux_y_back.resize(nv);
      Dflux_z_center.resize(nv);
      Dflux_z_top.resize(nv);
      Ksx.newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
      Ksy.newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
      Ksz.newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
      Ksx = 0.0;
      Ksy = 0.0;
      Ksz = 0.0;
      for (vi=0;vi<nv;vi++)
        {
          flux_x[vi].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
          flux_y[vi].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
          flux_z[vi].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
          flux_x[vi] = 0.0;
          flux_y[vi] = 0.0;
          flux_z[vi] = 0.0;
          Dflux_x_center[vi].resize(nv);
          Dflux_x_right[vi].resize(nv);
          Dflux_y_center[vi].resize(nv);
          Dflux_y_back[vi].resize(nv);
          Dflux_z_center[vi].resize(nv);
          Dflux_z_top[vi].resize(nv);
          for (vj=0;vj<nv;vj++)
            {
              Dflux_x_center[vi][vj].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
              Dflux_x_center[vi][vj] = 0.0;
              Dflux_x_right[vi][vj].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
              Dflux_x_right[vi][vj] = 0.0;

              Dflux_y_center[vi][vj].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
              Dflux_y_center[vi][vj] = 0.0;
              Dflux_y_back[vi][vj].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
              Dflux_y_back[vi][vj] = 0.0;

              Dflux_z_center[vi][vj].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
              Dflux_z_center[vi][vj] = 0.0;
              Dflux_z_top[vi][vj].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
              Dflux_z_top[vi][vj] = 0.0;
            }
        }
    }

template<class BC, int nv>
DivKgrad3dUGMM<BC,nv>::~DivKgrad3dUGMM(){}

template<class BC, int nv>
void DivKgrad3dUGMM<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P)
{
  //ghost faces
  if (node.ghost_zOffSet)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(-node.ghost_zOffSet,j,k);
          setFlux_z(K,Rho,P);
        }
  if (node.ghost_nzNodes -  node.ghost_zOffSet - node.local_nzNodes)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(local_nzNodes-1,j,k);
          setFlux_z(K,Rho,P);
        }
  if (node.ghost_yOffSet)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,-node.ghost_yOffSet,k);
          setFlux_y(K,Rho,P);
        }
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,local_nyNodes-1,k);
          setFlux_y(K,Rho,P);
        }
  if (node.ghost_xOffSet)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,-node.ghost_xOffSet);
          setFlux_x(K,Rho,P);
        }      
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,local_nxNodes-1);
          setFlux_x(K,Rho,P);
        }
  
  //interior
  for (i=0;i<local_nzNodes-1;i++)
    {
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(i,j,k);
              setFlux_x(K,Rho,P);
              setFlux_y(K,Rho,P);
              setFlux_z(K,Rho,P);
            }
          // last y and z fluxes of x row
          node.localIndex(i,j,k);
          setFlux_y(K,Rho,P);
          setFlux_z(K,Rho,P);
        }
      //last row of x and z fluxes in y plane
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setFlux_x(K,Rho,P);
          setFlux_z(K,Rho,P);
        }
      // last z fluz in x row
      node.localIndex(i,j,k);
      setFlux_z(K,Rho,P);
    }
  //top face of x and y fluxes
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setFlux_x(K,Rho,P);
          setFlux_y(K,Rho,P);
        }
      //last y flux in x row
      node.localIndex(i,j,k);
      setFlux_y(K,Rho,P);
    }
  //last row of x fluxes in top face
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(i,j,k);
      setFlux_x(K,Rho,P);
    }
  
  for (vi=0;vi<nv;vi++)
    bc[vi].applyNeumannConditions(node,flux_x[vi], flux_y[vi], flux_z[vi]);
}

template<class BC, int nv>
void DivKgrad3dUGMM<BC,nv>::computeInterfaceK(const Vec& K)
{
  //ghost faces
  if (node.ghost_zOffSet)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(-node.ghost_zOffSet,j,k);
          setKsz(K);
        }
  if (node.ghost_nzNodes -  node.ghost_zOffSet - node.local_nzNodes)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(local_nzNodes-1,j,k);
          setKsz(K);
        }
  if (node.ghost_yOffSet)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,-node.ghost_yOffSet,k);
          setKsy(K);
        }
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,local_nyNodes-1,k);
          setKsy(K);
        }
  if (node.ghost_xOffSet)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,-node.ghost_xOffSet);
          setKsx(K);
        }      
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,local_nxNodes-1);
          setKsx(K);
        }
  
  //interior
  for (i=0;i<local_nzNodes-1;i++)
    {
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(i,j,k);
              setKsx(K);
              setKsy(K);
              setKsz(K);
            }
          // last y and z fluxes of x row
          node.localIndex(i,j,k);
          setKsy(K);
          setKsz(K);
        }
      //last row of x and z fluxes in y plane
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setKsx(K);
          setKsz(K);
        }
      // last z fluz in x row
      node.localIndex(i,j,k);
      setKsz(K);
    }
  //top face of x and y fluxes
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setKsx(K);
          setKsy(K);
        }
      //last y flux in x row
      node.localIndex(i,j,k);
      setKsy(K);
    }
  //last row of x fluxes in top face
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(i,j,k);
      setKsx(K);
    }
}

template<class BC, int nv>
void DivKgrad3dUGMM<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, Vec* Div)
{
  computeDivergence(K,Rho,P);
  for (i=0;i<local_nzNodes;i++)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,j,k);
          for (vi=0;vi<nv;vi++)
            {
              Div[vi][node.center_noGhost] = 
                oneOverdx*(flux_x[vi][node.interLeft] - 
                           flux_x[vi][node.interRight])
                + oneOverdy*(flux_y[vi][node.interFront] - 
                             flux_y[vi][node.interBack])
                + oneOverdz*(flux_z[vi][node.interBottom] - 
                             flux_z[vi][node.interTop]);
            }
        }
}
  
template<class BC, int nv>
void DivKgrad3dUGMM<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho)
{
  //ghost faces
  if (node.ghost_zOffSet)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(-node.ghost_zOffSet,j,k);
          setDFlux_z(K,Rho,P,DK,DRho);
        }
  if (node.ghost_nzNodes -  node.ghost_zOffSet - node.local_nzNodes)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(local_nzNodes-1,j,k);
          setDFlux_z(K,Rho,P,DK,DRho);
        }
  if (node.ghost_yOffSet)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,-node.ghost_yOffSet,k);
          setDFlux_y(K,Rho,P,DK,DRho);
        }
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,local_nyNodes-1,k);
          setDFlux_y(K,Rho,P,DK,DRho);
        }
  if (node.ghost_xOffSet)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,-node.ghost_xOffSet);
          setDFlux_x(K,Rho,P,DK,DRho);
        }      
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,local_nxNodes-1);
          setDFlux_x(K,Rho,P,DK,DRho);
        }
  
  for (i=0;i<local_nzNodes-1;i++)
    {
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              //interior
              node.localIndex(i,j,k);
              setDFlux_x(K,Rho,P,DK,DRho);
              setDFlux_y(K,Rho,P,DK,DRho);
              setDFlux_z(K,Rho,P,DK,DRho);
            }
          // last y and z fluxes of x row
          node.localIndex(i,j,k);
          setDFlux_y(K,Rho,P,DK,DRho);
          setDFlux_z(K,Rho,P,DK,DRho);
        }
      //last row of x and z fluxes in y plane
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setDFlux_x(K,Rho,P,DK,DRho);
          setDFlux_z(K,Rho,P,DK,DRho);
        }
      // last z fluz in x row
      node.localIndex(i,j,k);
      setDFlux_z(K,Rho,P,DK,DRho);
    }
  //top face of x and y fluxes
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setDFlux_x(K,Rho,P,DK,DRho);
          setDFlux_y(K,Rho,P,DK,DRho);
        }
      //last y flux in x row
      node.localIndex(i,j,k);
      setDFlux_y(K,Rho,P,DK,DRho);
    }
  //last row of x fluxes in top face
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(i,j,k);
      setDFlux_x(K,Rho,P,DK,DRho);
    }
}

template<class BC, int nv>
void DivKgrad3dUGMM<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Rho,P,DK,DRho);
  for (i=0;i<local_nzNodes;i++)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,j,k);
          for (vi=0;vi<nv;vi++)
            {
              for (vj=0;vj<nv;vj++)
                {
                  
                  DivJac[vi][vj][BOTTOM][node.center_noGhost] = oneOverdz*
                    Dflux_z_center[vi][vj][node.interBottom];
                  
                  DivJac[vi][vj][FRONT][node.center_noGhost] = oneOverdy*
                    Dflux_y_center[vi][vj][node.interFront];
                  
                  DivJac[vi][vj][LEFT][node.center_noGhost] = oneOverdx*
                    Dflux_x_center[vi][vj][node.interLeft];
                  
                  DivJac[vi][vj][CENTER][node.center_noGhost] = oneOverdx*
                    (Dflux_x_right[vi][vj][node.interLeft] 
                     -Dflux_x_center[vi][vj][node.interRight]) 
                    + oneOverdy*
                    (Dflux_y_back[vi][vj][node.interFront] 
                     -Dflux_y_center[vi][vj][node.interBack]) 
                    + oneOverdz*
                    (Dflux_z_top[vi][vj][node.interBottom] 
                     -Dflux_z_center[vi][vj][node.interTop]); 
                  
                  DivJac[vi][vj][RIGHT][node.center_noGhost]  = -oneOverdx*
                    Dflux_x_right[vi][vj][node.interRight]; 
                  
                  DivJac[vi][vj][BACK][node.center_noGhost]  = -oneOverdy*
                    Dflux_y_back[vi][vj][node.interBack]; 
                  
                  DivJac[vi][vj][TOP][node.center_noGhost]  = -oneOverdz*
                    Dflux_z_top[vi][vj][node.interTop]; 
                }
            }
        }
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          bc[vi].applyNeumannDerivatives(node,DivJac[vi][vj],
                                         Dflux_x_center[vi][vj],Dflux_x_right[vi][vj],
                                         Dflux_y_center[vi][vj],Dflux_y_back[vi][vj],
                                         Dflux_z_center[vi][vj],Dflux_z_top[vi][vj],
                                         oneOverdx,oneOverdy,oneOverdz);
        }
    }  
}

}//Daetk
#endif
