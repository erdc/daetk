#ifndef DIVKGRAD_H
#define DIVKGRAD_H
#include "Definitions.h"
#include "PetscSecondOrderFd.h"
#include "Divergence.h"
#include <vector>

//  FIX DIVERGENCE JAC -- MUST APPLY CORRECT DERIVATIVES FOR NEUMAN BC'S

namespace Daetk 
{

using Petsc::SecondOrderFd;

template<class BC, int nv>
class DivKgrad : public Divergence
{
public:
  enum JacEntry {CENTER, LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP};
  DivKgrad(BC* bcIn, SecondOrderFd& nodeIn, real* vis_ratio=0);
  
  virtual ~DivKgrad();

  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div)=0;

  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)=0;  
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P)=0;

  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho)=0;  
  
  virtual void computeInterfaceK(const Vec& K)=0;

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
  
  int nxNodes,nyNodes,nzNodes,local_nxNodes,local_nyNodes,local_nzNodes;
  real oneOverdx,oneOverdy,oneOverdz,gx,gy,gz; 
  real viscosity_ratio[nv];
  Vec flux_x[nv],flux_y[nv],flux_z[nv],Ksx,Ksy,Ksz;
  VecVecVec Dflux_x_center,Dflux_x_right,Dflux_y_center,Dflux_y_back,Dflux_z_center,Dflux_z_top;
  BC* bc;
  SecondOrderFd& node;
};

template<class BC, int nv>
class DivKgrad1d : public DivKgrad<BC,nv>
{
public:
  DivKgrad1d(BC* bcIn, 
             SecondOrderFd& nodeIn,
             int nNodes, 
             real g, 
             real oneOverd, 
             real* vis_ratio=0);

  virtual ~DivKgrad1d();
  
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
class DivKgrad2d : public DivKgrad<BC,nv>
{
public:
  DivKgrad2d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn, 
             real gxIn, real gyIn,
             real oneOverdxIn, real oneOverdyIn, real* vis_ratio=0);
  virtual ~DivKgrad2d();
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
class DivKgrad3d : public DivKgrad<BC,nv>
{
public:
  DivKgrad3d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn,int nzNodesIn, 
             real gxIn, real gyIn, real gzIn,
             real oneOverdxIn, real oneOverdyIn, real oneOverdzIn,
             real* vis_ratio=0);
  virtual ~DivKgrad3d();
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
DivKgrad<BC,nv>::DivKgrad(BC* bcIn, SecondOrderFd& nodeIn, real* vis_ratio):
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
        this->viscosity_ratio[i] = vis_ratio[i];
    }
  else
    for (int i=0;i<nv;i++)
      this->viscosity_ratio[i]=1.0;
  
  this->Dflux_x_center.resize(nv);
  this->Dflux_x_right.resize(nv);
  this->Dflux_y_center.resize(nv);
  this->Dflux_y_back.resize(nv);
  this->Dflux_z_center.resize(nv);
  this->Dflux_z_top.resize(nv);
  this->Ksx.newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
  this->Ksy.newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
  this->Ksz.newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
  this->Ksx = 0.0;
  this->Ksy = 0.0;
  this->Ksz = 0.0;
  for (int vi=0;vi<nv;vi++)
    {
      this->flux_x[vi].newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
      this->flux_y[vi].newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
      this->flux_z[vi].newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
      this->flux_x[vi] = 0.0;
      this->flux_y[vi] = 0.0;
      this->flux_z[vi] = 0.0;
      this->Dflux_x_center[vi].resize(nv);
      this->Dflux_x_right[vi].resize(nv);
      this->Dflux_y_center[vi].resize(nv);
      this->Dflux_y_back[vi].resize(nv);
      this->Dflux_z_center[vi].resize(nv);
      this->Dflux_z_top[vi].resize(nv);
      for (int vj=0;vj<nv;vj++)
	{
	  this->Dflux_x_center[vi][vj].newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
	  this->Dflux_x_center[vi][vj] = 0.0;
	  this->Dflux_x_right[vi][vj].newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
	  this->Dflux_x_right[vi][vj] = 0.0;
	  
	  this->Dflux_y_center[vi][vj].newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
	  this->Dflux_y_center[vi][vj] = 0.0;
	  this->Dflux_y_back[vi][vj].newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
	  this->Dflux_y_back[vi][vj] = 0.0;
	  
	  this->Dflux_z_center[vi][vj].newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
	  this->Dflux_z_center[vi][vj] = 0.0;
	  this->Dflux_z_top[vi][vj].newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
	  this->Dflux_z_top[vi][vj] = 0.0;
	}
    }
}

template<class BC, int nv>
DivKgrad<BC,nv>::~DivKgrad(){}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setKsx(const Vec& K)
{
  this->Ksx[this->node.interRight] = 2.0 / ( 1.0/K[this->node.center] + 1.0/K[this->node.right]);
}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setKsy(const Vec& K)
{
  this->Ksy[this->node.interBack] = 2.0 / ( 1.0/K[this->node.center] + 1.0/K[this->node.back]);
}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setKsz(const Vec& K)
{
  this->Ksz[this->node.interTop] = 2.0 / ( 1.0/K[this->node.center] + 1.0/K[this->node.top]);
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::Kx(const Vec& K, const real& grad )
{ 
  return 0.5*(K[this->node.center]+K[this->node.right])*this->Ksx[this->node.interRight];
}


template<class BC, int nv>
inline real DivKgrad<BC,nv>::DKx_center(const Vec& DK, const real& grad)
{
  return 0.5*DK[this->node.center]*this->Ksx[this->node.interRight];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DKx_right(const Vec& DK, const real& grad)
{
  return 0.5*DK[this->node.right]*this->Ksx[this->node.interRight];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::Ky(const Vec& K, const real& grad)
{ 
  return 0.5*(K[this->node.center]+K[this->node.back])*this->Ksy[this->node.interBack];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DKy_center(const Vec& DK,const real& grad)
{ 
  return 0.5*DK[this->node.center]*this->Ksy[this->node.interBack];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DKy_back(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[this->node.back]*this->Ksy[this->node.interBack];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::Kz(const Vec& K, const real& grad)
{ 
  return 0.5*(K[this->node.center]+K[this->node.top])*this->Ksz[this->node.interTop];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DKz_center(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[this->node.center]*this->Ksz[this->node.interTop];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DKz_top(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[this->node.top]*this->Ksz[this->node.interTop];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::Rhox(const Vec& Rho)
{ 
  return 0.5*(Rho[this->node.center]+Rho[this->node.right]);
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DRhox_center(const Vec& DRho)
{ 
  return 0.5*DRho[this->node.center];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DRhox_right(const Vec& DRho)
{ 
  return 0.5*DRho[this->node.right];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::Rhoy(const Vec& Rho)
{ 
  return 0.5*(Rho[this->node.center]+Rho[this->node.back]);
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DRhoy_center(const Vec& DRho)
{ 
  return 0.5*DRho[this->node.center];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DRhoy_back(const Vec& DRho)
{ 
  return 0.5*DRho[this->node.back];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::Rhoz(const Vec& Rho)
{ 
  return 0.5*(Rho[this->node.center]+Rho[this->node.top]);
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DRhoz_center(const Vec& DRho)
{ 
  return 0.5*DRho[this->node.center];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DRhoz_top(const Vec& DRho)
{ 
  return 0.5*DRho[this->node.top];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DpDx(const Vec& P)
{ 
  return this->oneOverdx*(P[this->node.right] - P[this->node.center]);
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DpDy(const Vec& P)
{ 
  return this->oneOverdy*(P[this->node.back]  - P[this->node.center]);
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::DpDz(const Vec& P)
{ 
  return this->oneOverdz*(P[this->node.top]   - P[this->node.center]);
}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setFlux_x(int v,const Vec& K, const Vec& Rho, const Vec& P)
{
  real rhox=Rhox(Rho),
    grad=DpDx(P) + rhox*this->gx;
  this->flux_x[v][this->node.interRight] = -rhox*this->viscosity_ratio[v]*
    Kx(K,grad) * grad;
}


template<class BC, int nv>
inline void DivKgrad<BC,nv>::setFlux_x(const Vec* K, const Vec* Rho, const Vec* P)
{	
  for (int vi=0;vi<nv;vi++)
    {
      real rhox=Rhox(Rho[vi]),
	grad=DpDx(P[vi]) + rhox*this->gx;
      this->flux_x[vi][this->node.interRight] = -rhox*this->viscosity_ratio[vi]*
        Kx(K[vi],grad) * grad;
    }
}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setFlux_y(const Vec* K, const Vec* Rho, const Vec* P)
{
  for (int vi=0;vi<nv;vi++)
    {
      real grad=DpDy(P[vi]) + Rhoy(Rho[vi])*this->gy;
      this->flux_y[vi][this->node.interBack]  = -Rhoy(Rho[vi])*this->viscosity_ratio[vi]*
        Ky(K[vi],grad) * grad;
    }
}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setFlux_z(const Vec* K, const Vec* Rho, const Vec* P)
{
  for (int vi=0;vi<nv;vi++)
    {
      real grad=DpDz(P[vi]) + Rhoz(Rho[vi])*this->gz;
      this->flux_z[vi][this->node.interTop]	  = -Rhoz(Rho[vi])*this->viscosity_ratio[vi]*
        Kz(K[vi],grad) * grad;
    }
}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setDFlux_x(const Vec* K, const Vec* Rho, const Vec* P, 
                                        const VecVecVec& DK, const Vec* DRho)
{
//   //it's very important that the compiler be able to unroll these
//   for (int vi=0;vi<nv;vi++)
//     {
//       for (int vj=0;vj<nv;vj++)
//         {
//           //interior
//           real grad=DpDx(P[vi]) + Rhox(Rho[vi])*this->gx;
//           this->Dflux_x_center[vi][vj][this->node.interRight] = -Rhox(Rho[vi])*
//             DKx_center(DK[vi][vj],grad)*grad;
          
//           this->Dflux_x_right[vi][vj][this->node.interRight] = -Rhox(Rho[vi])*
//             DKx_right(DK[vi][vj],grad) * grad;
          
//         }
//     }
  
//   for (int vi=0;vi<nv;vi++)
//     {
//       real grad=DpDx(P[vi]) + Rhox(Rho[vi])*this->gx;
//       this->Dflux_x_center[vi][vi][this->node.interRight] += 
//         -Rhox(Rho[vi])*Kx(K[vi],grad)*
//         (-this->oneOverdx + DRhox_center(DRho[vi])*this->gx) 
//         -DRhox_center(DRho[vi])*Kx(K[vi],grad)*grad;
      
//       this->Dflux_x_right[vi][vi][this->node.interRight] +=
//         -Rhox(Rho[vi])*Kx(K[vi],grad)*
//         (this->oneOverdx + DRhox_right(DRho[vi])*this->gx) 
//         -DRhox_right(DRho[vi])*Kx(K[vi],grad)*grad;
//     }

  //it's very important that the compiler be able to unroll these
  for (int vi=0;vi<nv;vi++)
    {
      real rhox=Rhox(Rho[vi]),
	grad=DpDx(P[vi]) + rhox*this->gx,
	kx=this->viscosity_ratio[vi]*Kx(K[vi],grad),
	kxGrad=kx*grad,
	mrhoxKx=-rhox*kx,
	mrhoxGrad=-rhox*grad,
	drhox_center=DRhox_center(DRho[vi]),
	drhox_right=DRhox_right(DRho[vi]);


      
      for (int vj=0;vj<nv;vj++)
	{
          //interior
          this->Dflux_x_center[vi][vj][this->node.interRight] = mrhoxGrad*
            DKx_center(DK[vi][vj],grad);
          
          this->Dflux_x_right[vi][vj][this->node.interRight] = mrhoxGrad*
            DKx_right(DK[vi][vj],grad);
          
	}

      this->Dflux_x_center[vi][vi][this->node.interRight] += 
        mrhoxKx*
        (-this->oneOverdx + drhox_center*this->gx) 
        -drhox_center*kxGrad;
      
      this->Dflux_x_right[vi][vi][this->node.interRight] +=
        mrhoxKx*
        (this->oneOverdx + drhox_right*this->gx) 
        -drhox_right*kxGrad;
    }
}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setDFlux_y(const Vec* K, const Vec* Rho, const Vec* P, 
                                        const VecVecVec& DK, const Vec* DRho)
{
  for (int vi=0;vi<nv;vi++)
    {
      for (int vj=0;vj<nv;vj++)
        {
          real grad=DpDy(P[vi]) + Rhoy(Rho[vi])*this->gy;
          this->Dflux_y_center[vi][vj][this->node.interBack] = -Rhoy(Rho[vi])*
            this->viscosity_ratio[vi]*DKy_center(DK[vi][vj],grad)*grad;
          
          this->Dflux_y_back[vi][vj][this->node.interBack] = -Rhoy(Rho[vi])*
            this->viscosity_ratio[vi]*DKy_back(DK[vi][vj],grad) * grad;
          
        }
    }
  
  for (int vi=0;vi<nv;vi++)
    {
      real grad=DpDy(P[vi]) + Rhoy(Rho[vi])*this->gy;
      this->Dflux_y_center[vi][vi][this->node.interBack] += 
        -Rhoy(Rho[vi])*this->viscosity_ratio[vi]*Ky(K[vi],grad)*
        (-this->oneOverdy + DRhoy_center(DRho[vi])*this->gy) 
        -DRhoy_center(DRho[vi])*this->viscosity_ratio[vi]*Ky(K[vi],grad)*grad;
      
      this->Dflux_y_back[vi][vi][this->node.interBack] +=
        -Rhoy(Rho[vi])*this->viscosity_ratio[vi]*Ky(K[vi],grad)*
        (this->oneOverdy + DRhoy_back(DRho[vi])*this->gy) 
        -DRhoy_back(DRho[vi])*this->viscosity_ratio[vi]*Ky(K[vi],grad)*grad;
    }
}

template<class BC, int nv>
inline void DivKgrad<BC,nv>::setDFlux_z(const Vec* K, const Vec* Rho, const Vec* P, 
                                        const VecVecVec& DK, const Vec* DRho)
{
  for (int vi=0;vi<nv;vi++)
    {
      for (int vj=0;vj<nv;vj++)
        {
          real grad=DpDz(P[vi]) + Rhoz(Rho[vi])*this->gz;
          this->Dflux_z_center[vi][vj][this->node.interTop] = -Rhoz(Rho[vi])*
            this->viscosity_ratio[vi]*DKz_center(DK[vi][vj],grad)*grad;
          
          this->Dflux_z_top[vi][vj][this->node.interTop] = -Rhoz(Rho[vi])*
            this->viscosity_ratio[vi]*DKz_top(DK[vi][vj],grad) *grad;
          
        }
    }
  
  for (int vi=0;vi<nv;vi++)
    { 
      real grad=DpDz(P[vi]) + Rhoz(Rho[vi])*this->gz;
      this->Dflux_z_center[vi][vi][this->node.interTop] += 
        -Rhoz(Rho[vi])*this->viscosity_ratio[vi]*Kz(K[vi],grad)*
        (-this->oneOverdz + DRhoz_center(DRho[vi])*this->gz) 
        -DRhoz_center(DRho[vi])*this->viscosity_ratio[vi]*Kz(K[vi],grad)*grad;
      
      this->Dflux_z_top[vi][vi][this->node.interTop] +=
        -Rhoz(Rho[vi])*this->viscosity_ratio[vi]*Kz(K[vi],grad)*
	  (this->oneOverdz + DRhoz_top(DRho[vi])*this->gz) 
        -DRhoz_top(DRho[vi])*this->viscosity_ratio[vi]*Kz(K[vi],grad)*grad;
    }
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::getDiv(SecondOrderFd& n, int div_i)
{
  return this->oneOverdx*(this->flux_x[div_i][this->node.interLeft] - 
                    this->flux_x[div_i][this->node.interRight])
    + this->oneOverdy*(this->flux_y[div_i][this->node.interFront] - 
                 this->flux_y[div_i][this->node.interBack])
    + this->oneOverdz*(this->flux_z[div_i][this->node.interBottom] - 
                 this->flux_z[div_i][this->node.interTop]);
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::getDivJacCenter(SecondOrderFd& n, int div_i, int u_j)
{
  return this->oneOverdx*
    (this->Dflux_x_right[div_i][u_j][this->node.interLeft]
     -this->Dflux_x_center[div_i][u_j][this->node.interRight])
    + this->oneOverdy*
    (this->Dflux_y_back[div_i][u_j][this->node.interFront]
     -this->Dflux_y_center[div_i][u_j][this->node.interBack])
    + this->oneOverdz*
    (this->Dflux_z_top[div_i][u_j][this->node.interBottom]
     -this->Dflux_z_center[div_i][u_j][this->node.interTop]); 
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::getDivJacLeft(SecondOrderFd& n, int div_i, int u_j)
{
  return this->oneOverdx*
    this->Dflux_x_center[div_i][u_j][this->node.interLeft];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::getDivJacRight(SecondOrderFd& n, int div_i, int u_j)
{
  return -this->oneOverdx*
    this->Dflux_x_right[div_i][u_j][this->node.interRight];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::getDivJacFront(SecondOrderFd& n, int div_i, int u_j)
{
  return this->oneOverdy*
    this->Dflux_y_center[div_i][u_j][this->node.interFront];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::getDivJacBack(SecondOrderFd& n, int div_i, int u_j)
{
  return -this->oneOverdy*
    this->Dflux_y_back[div_i][u_j][this->node.interBack];
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::getDivJacBottom(SecondOrderFd& n, int div_i, int u_j)
{
  return this->oneOverdz*
    this->Dflux_z_center[div_i][u_j][this->node.interBottom]; 
}

template<class BC, int nv>
inline real DivKgrad<BC,nv>::getDivJacTop(SecondOrderFd& n, int div_i, int u_j)
{
  return -this->oneOverdz*
    this->Dflux_z_top[div_i][u_j][this->node.interTop];
}


template<class BC, int nv>
const Vec& DivKgrad<BC,nv>::getFlux_x(int n){return this->flux_x[n];}

template<class BC, int nv>
const Vec& DivKgrad<BC,nv>::getFlux_y(int n){return this->flux_y[n];}

template<class BC, int nv>
const Vec& DivKgrad<BC,nv>::getFlux_z(int n){return this->flux_z[n];}

template<class BC, int nv>
Vec& DivKgrad<BC,nv>::setFlux_x(int n){return this->flux_x[n];}

template<class BC, int nv>
Vec& DivKgrad<BC,nv>::setFlux_y(int n){return this->flux_y[n];}

template<class BC, int nv>
Vec& DivKgrad<BC,nv>::setFlux_z(int n){return this->flux_z[n];}

template<class BC, int nv>
DivKgrad1d<BC,nv>::DivKgrad1d(BC* bcIn, 
             SecondOrderFd& nodeIn,
             int nNodes, 
             real gxIn, 
             real oneOverdxIn, real* vis_ratio):
   DivKgrad<BC,nv>(bcIn,nodeIn,vis_ratio)
  {
    this->nxNodes=nNodes;
    this->local_nxNodes = nodeIn.local_nxNodes;

    this->gx=gxIn;
    this->oneOverdx=oneOverdxIn;

    this->Dflux_x_center.resize(nv);
    this->Dflux_x_right.resize(nv);
    this->Ksx.newsize(Vec::LOCAL,this->local_nxNodes+1);
    this->Ksx = 0.0;
    for (int vi=0;vi<nv;vi++)
      {
        this->flux_x[vi].newsize(Vec::LOCAL,this->local_nxNodes+1);
        this->flux_x[vi] = 0.0;
        this->Dflux_x_center[vi].resize(nv);
        this->Dflux_x_right[vi].resize(nv);
        for (int vj=0;vj<nv;vj++)
          {
            this->Dflux_x_center[vi][vj].newsize(Vec::LOCAL,this->local_nxNodes+1);
            this->Dflux_x_center[vi][vj] = 0.0;
            this->Dflux_x_right[vi][vj].newsize(Vec::LOCAL,this->local_nxNodes+1);
            this->Dflux_x_right[vi][vj] = 0.0;
          }
      }
  }

template<class BC, int nv>
DivKgrad1d<BC,nv>::~DivKgrad1d(){}

template<class BC, int nv>
void DivKgrad1d<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P)
{  
  for (int k= -this->node.ghost_xOffSet; k< (this->node.ghost_nxNodes - this->node.ghost_xOffSet-1); k++)
    {
      this->node.localIndex(k);
      DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
    }
  
  for (int vi=0; vi<nv; vi++)
    this->bc[vi].applyNeumannConditions(this->node,this->flux_x[vi]);
}  

template<class BC, int nv>
void DivKgrad1d<BC,nv>::computeInterfaceK(const Vec& K)
{
  for (int k= -this->node.ghost_xOffSet; k< (this->node.ghost_nxNodes - this->node.ghost_xOffSet-1); k++)
    {
      this->node.localIndex(k);
      DivKgrad<BC,nv>::setKsx(K);
    }
}  

template<class BC, int nv>
void DivKgrad1d<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, Vec* Div)
{
  computeDivergence(K,Rho,P);
  for (int k=0; k<this->local_nxNodes; k++)
    {
      this->node.localIndex(k);
      for (int vi=0;vi<nv;vi++)
        {
          Div[vi][this->node.center_noGhost] = this->oneOverdx*
            (this->flux_x[vi][this->node.interLeft] - this->flux_x[vi][this->node.interRight]);
        }
    }
}  

template<class BC, int nv>
void DivKgrad1d<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                             const VecVecVec& DK, const Vec* DRho)
{
  for (int k=-this->node.ghost_xOffSet; k< this->node.ghost_nxNodes - this->node.ghost_xOffSet-1; k++)
    {
      this->node.localIndex(k);
      DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
    }
}  

template<class BC, int nv>
void DivKgrad1d<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                             const VecVecVec& DK, const Vec* DRho, 
                                             VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Rho,P,DK,DRho);
  for (int k=0;k<this->local_nxNodes;k++)
    {
      this->node.localIndex(k);
      for (int vi=0;vi<nv;vi++)
        {
          for (int vj=0;vj<nv;vj++)
            {
              DivJac[vi][vj][this->LEFT][this->node.center_noGhost] = this->oneOverdx*
                this->Dflux_x_center[vi][vj][this->node.interLeft];
              
              DivJac[vi][vj][this->CENTER][this->node.center_noGhost] = this->oneOverdx*
                (this->Dflux_x_right[vi][vj][this->node.interLeft] 
                 -this->Dflux_x_center[vi][vj][this->node.interRight]); 
              
              DivJac[vi][vj][this->RIGHT][this->node.center_noGhost]  = -this->oneOverdx*
                this->Dflux_x_right[vi][vj][this->node.interRight]; 
            }
        }
    }
  
  for (int vi=0;vi<nv;vi++)
    {
      for (int vj=0;vj<nv;vj++)
        {
          this->bc[vi].applyNeumannDerivatives(this->node,DivJac[vi][vj],
                                         this->Dflux_x_center[vi][vj],this->Dflux_x_right[vi][vj],
                                         this->oneOverdx);
        }
    }
}  

template<class BC, int nv>
DivKgrad2d<BC,nv>::DivKgrad2d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn, 
             real gxIn, real gyIn,
             real oneOverdxIn, real oneOverdyIn, real* vis_ratio):
    DivKgrad<BC,nv>(bcIn,nodeIn,vis_ratio)
  {
    this->nxNodes=nxNodesIn;
    this->local_nxNodes = nodeIn.local_nxNodes;
    this->gx=gxIn;
    this->oneOverdx = oneOverdxIn;

    this->nyNodes=nyNodesIn;
    this->local_nyNodes = nodeIn.local_nyNodes;
    this->gy=gyIn;
    this->oneOverdy = oneOverdyIn;

    this->Dflux_x_center.resize(nv);
    this->Dflux_x_right.resize(nv);
    this->Dflux_y_center.resize(nv);
    this->Dflux_y_back.resize(nv);
    this->Ksx.newsize(Vec::LOCAL,this->local_nyNodes*(this->local_nxNodes+1));
    this->Ksy.newsize(Vec::LOCAL,(this->local_nyNodes+1)*this->local_nxNodes);
    this->Ksx = 0.0;
    this->Ksy = 0.0;
    for (int vi=0;vi<nv;vi++)
      {
          this->flux_x[vi].newsize(Vec::LOCAL,this->local_nyNodes*(this->local_nxNodes+1));
          this->flux_y[vi].newsize(Vec::LOCAL,(this->local_nyNodes+1)*this->local_nxNodes);
          this->flux_x[vi] = 0.0;
          this->flux_y[vi] = 0.0;
          this->Dflux_x_center[vi].resize(nv);
          this->Dflux_x_right[vi].resize(nv);
          this->Dflux_y_center[vi].resize(nv);
          this->Dflux_y_back[vi].resize(nv);
          for (int vj=0;vj<nv;vj++)
            {
              this->Dflux_x_center[vi][vj].newsize(Vec::LOCAL,this->local_nyNodes*(this->local_nxNodes+1));
              this->Dflux_x_center[vi][vj] = 0.0;
              this->Dflux_x_right[vi][vj].newsize(Vec::LOCAL,this->local_nyNodes*(this->local_nxNodes+1));
              this->Dflux_x_right[vi][vj] = 0.0;
              this->Dflux_y_center[vi][vj].newsize(Vec::LOCAL,(this->local_nyNodes+1)*this->local_nxNodes);
              this->Dflux_y_center[vi][vj] = 0.0;
              this->Dflux_y_back[vi][vj].newsize(Vec::LOCAL,(this->local_nyNodes+1)*this->local_nxNodes);
              this->Dflux_y_back[vi][vj] = 0.0;
            }
        }
    }

template<class BC, int nv>
DivKgrad2d<BC,nv>::~DivKgrad2d(){}

template<class BC, int nv>
void DivKgrad2d<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P)
{
  //bottom ghost line
  if (this->node.ghost_yOffSet)
    for(int k=0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(-this->node.ghost_yOffSet,k);
        DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
      }
  //top ghost line
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for(int k=0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(this->local_nyNodes-1,k);
        DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
      }
  
  //left ghost line
  if (this->node.ghost_xOffSet)
    for (int j=0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,-this->node.ghost_xOffSet);
        DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
      }
  
  //right ghost line
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int j=0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,this->local_nxNodes-1);
        DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
      }
  
  for (int j=0;j<this->local_nyNodes-1;j++)
    {
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(j,k);
          DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
          DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
        }
      //last y flux in row
      this->node.localIndex(j,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
    }
  //last row of x fluxes
  for (int k=0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nyNodes-1,k);
      DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
    }
  
  for (int vi=0;vi<nv;vi++)
    this->bc[vi].applyNeumannConditions(this->node,this->flux_x[vi], this->flux_y[vi]);
}

template<class BC, int nv>
void DivKgrad2d<BC,nv>::computeInterfaceK(const Vec& K)
{
  //bottom ghost line
  if (this->node.ghost_yOffSet)
    for(int k=0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(-this->node.ghost_yOffSet,k);
        DivKgrad<BC,nv>::setKsy(K);
      }
  //top ghost line
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for(int k=0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(this->local_nyNodes-1,k);
        DivKgrad<BC,nv>::setKsy(K);
      }
  
  //left ghost line
  if (this->node.ghost_xOffSet)
    for (int j=0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,-this->node.ghost_xOffSet);
        DivKgrad<BC,nv>::setKsx(K);
      }
  
  //right ghost line
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int j=0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,this->local_nxNodes-1);
        DivKgrad<BC,nv>::setKsx(K);
      }
  
  for (int j=0;j<this->local_nyNodes-1;j++)
    {
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(j,k);
          DivKgrad<BC,nv>::setKsx(K);
          DivKgrad<BC,nv>::setKsy(K);
        }
      //last y flux in row
      this->node.localIndex(j,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setKsy(K);
    }
  //last row of x fluxes
  for (int k=0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nyNodes-1,k);
      DivKgrad<BC,nv>::setKsx(K);
    }
}

template<class BC, int nv>
void DivKgrad2d<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, Vec* Div)
{
  computeDivergence(K,Rho,P);
  for (int j=0;j<this->local_nyNodes;j++)
    for (int k=0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(j,k);
        for (int vi=0;vi<nv;vi++)
          {
            Div[vi][this->node.center_noGhost] = this->oneOverdx*
              (this->flux_x[vi][this->node.interLeft] - this->flux_x[vi][this->node.interRight])
              + this->oneOverdy*
              (this->flux_y[vi][this->node.interFront] - this->flux_y[vi][this->node.interBack]);
          }
      }
}

template<class BC, int nv>
void DivKgrad2d<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                             const VecVecVec& DK, const Vec* DRho)
{
  //bottom ghost line
  if (this->node.ghost_yOffSet)
    for(int k=0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(-this->node.ghost_yOffSet,k);
        DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
      }
  //top ghost line
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for(int k=0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(this->local_nyNodes-1,k);
        DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
      }
  
  //left ghost line
  if (this->node.ghost_xOffSet)
    for (int j=0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,-this->node.ghost_xOffSet);
        DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
      }
  
  //right ghost line
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int j=0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,this->local_nxNodes-1);
        DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
      }
  
  for (int j=0;j<this->local_nyNodes-1;j++)
    {
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(j,k);
          DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
          DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
        }
      //last y flux in row
      this->node.localIndex(j,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
    }
  //last row of x fluxes
  for (int k=0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nyNodes-1,k);
      DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
    }
}

template<class BC, int nv>
void DivKgrad2d<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                             const VecVecVec& DK, const Vec* DRho, 
                                             VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Rho,P,DK,DRho);
  for (int j=0;j<this->local_nyNodes;j++)
    for (int k=0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(j,k);
        for (int vi=0;vi<nv;vi++)
          {
            for (int vj=0;vj<nv;vj++)
              {
                
                DivJac[vi][vj][this->FRONT][this->node.center_noGhost] = this->oneOverdy*
                  this->Dflux_y_center[vi][vj][this->node.interFront];
                
                DivJac[vi][vj][this->LEFT][this->node.center_noGhost] = this->oneOverdx*
                  this->Dflux_x_center[vi][vj][this->node.interLeft];
                
                DivJac[vi][vj][this->CENTER][this->node.center_noGhost] = this->oneOverdx*
                  (this->Dflux_x_right[vi][vj][this->node.interLeft] 
                   -this->Dflux_x_center[vi][vj][this->node.interRight]) 
                  + this->oneOverdy*
                  (this->Dflux_y_back[vi][vj][this->node.interFront] 
                   -this->Dflux_y_center[vi][vj][this->node.interBack]); 
                
                DivJac[vi][vj][this->RIGHT][this->node.center_noGhost]  = -this->oneOverdx*
                  this->Dflux_x_right[vi][vj][this->node.interRight]; 
                
                DivJac[vi][vj][this->BACK][this->node.center_noGhost]  = -this->oneOverdy*
                  this->Dflux_y_back[vi][vj][this->node.interBack]; 
                
              }
          }
      }
  
  for (int vi=0;vi<nv;vi++)
    {
      for (int vj=0;vj<nv;vj++)
        {
          this->bc[vi].applyNeumannDerivatives(this->node,DivJac[vi][vj],
                                         this->Dflux_x_center[vi][vj],this->Dflux_x_right[vi][vj],
                                         this->Dflux_y_center[vi][vj],this->Dflux_y_back[vi][vj],
                                         this->oneOverdx,this->oneOverdy);
          
        }
    }
}


template<class BC, int nv>
DivKgrad3d<BC,nv>::DivKgrad3d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn,int nzNodesIn, 
             real gxIn, real gyIn, real gzIn,
             real oneOverdxIn, real oneOverdyIn, real oneOverdzIn, real* vis_ratio):
    DivKgrad<BC,nv>(bcIn,nodeIn,vis_ratio)
    {
      this->nxNodes=nxNodesIn;
      this->local_nxNodes=nodeIn.local_nxNodes;
      this->gx=gxIn;
      this->oneOverdx = oneOverdxIn;

      this->nyNodes=nyNodesIn;
      this->local_nyNodes=nodeIn.local_nyNodes;
      this->gy=gyIn;
      this->oneOverdy = oneOverdyIn;

      this->nzNodes=nzNodesIn;
      this->local_nzNodes=nodeIn.local_nzNodes;
      this->gz=gzIn;
      this->oneOverdz = oneOverdzIn;

      this->Dflux_x_center.resize(nv);
      this->Dflux_x_right.resize(nv);
      this->Dflux_y_center.resize(nv);
      this->Dflux_y_back.resize(nv);
      this->Dflux_z_center.resize(nv);
      this->Dflux_z_top.resize(nv);
      this->Ksx.newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
      this->Ksy.newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
      this->Ksz.newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
      this->Ksx = 0.0;
      this->Ksy = 0.0;
      this->Ksz = 0.0;
      for (int vi=0;vi<nv;vi++)
        {
          this->flux_x[vi].newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
          this->flux_y[vi].newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
          this->flux_z[vi].newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
          this->flux_x[vi] = 0.0;
          this->flux_y[vi] = 0.0;
          this->flux_z[vi] = 0.0;
          this->Dflux_x_center[vi].resize(nv);
          this->Dflux_x_right[vi].resize(nv);
          this->Dflux_y_center[vi].resize(nv);
          this->Dflux_y_back[vi].resize(nv);
          this->Dflux_z_center[vi].resize(nv);
          this->Dflux_z_top[vi].resize(nv);
          for (int vj=0;vj<nv;vj++)
            {
              this->Dflux_x_center[vi][vj].newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
              this->Dflux_x_center[vi][vj] = 0.0;
              this->Dflux_x_right[vi][vj].newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
              this->Dflux_x_right[vi][vj] = 0.0;

              this->Dflux_y_center[vi][vj].newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
              this->Dflux_y_center[vi][vj] = 0.0;
              this->Dflux_y_back[vi][vj].newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
              this->Dflux_y_back[vi][vj] = 0.0;

              this->Dflux_z_center[vi][vj].newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
              this->Dflux_z_center[vi][vj] = 0.0;
              this->Dflux_z_top[vi][vj].newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
              this->Dflux_z_top[vi][vj] = 0.0;
            }
        }
    }

template<class BC, int nv>
DivKgrad3d<BC,nv>::~DivKgrad3d(){}

template<class BC, int nv>
void DivKgrad3d<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P)
{
  //ghost faces
  if (this->node.ghost_zOffSet)
    for (int j=0;j<this->local_nyNodes;j++)
      for (int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(-this->node.ghost_zOffSet,j,k);
          DivKgrad<BC,nv>::setFlux_z(K,Rho,P);
        }
  if (this->node.ghost_nzNodes -  this->node.ghost_zOffSet - this->node.local_nzNodes)
    for (int j=0;j<this->local_nyNodes;j++)
      for (int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          DivKgrad<BC,nv>::setFlux_z(K,Rho,P);
        }
  if (this->node.ghost_yOffSet)
    for (int i=0;i<this->local_nzNodes;i++)
      for(int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,-this->node.ghost_yOffSet,k);
          DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
        }
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for (int i=0;i<this->local_nzNodes;i++)
      for(int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
        }
  if (this->node.ghost_xOffSet)
    for (int i=0;i<this->local_nzNodes;i++)
      for (int j=0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,-this->node.ghost_xOffSet);
          DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
        }      
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int i=0;i<this->local_nzNodes;i++)
      for (int j=0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,this->local_nxNodes-1);
          DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
        }
  
  //interior
  for (int i=0;i<this->local_nzNodes-1;i++)
    {
      for (int j=0;j<this->local_nyNodes-1;j++)
        {
          for (int k=0;k<this->local_nxNodes-1;k++)
            {
              this->node.localIndex(i,j,k);
              DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
              DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
              DivKgrad<BC,nv>::setFlux_z(K,Rho,P);
            }
          // last y and z fluxes of x row
          this->node.localIndex(i,j,this->local_nxNodes-1);
          DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
          DivKgrad<BC,nv>::setFlux_z(K,Rho,P);
        }
      //last row of x and z fluxes in y plane
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
          DivKgrad<BC,nv>::setFlux_z(K,Rho,P);
        }
      // last z fluz in x row
      this->node.localIndex(i,this->local_nyNodes-1,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setFlux_z(K,Rho,P);
    }
  //top face of x and y fluxes
  for (int j=0;j<this->local_nyNodes-1;j++)
    {
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
          DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
        }
      //last y flux in x row
      this->node.localIndex(this->local_nzNodes-1,j,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setFlux_y(K,Rho,P);
    }
  //last row of x fluxes in top face
  for (int k=0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nzNodes-1,this->local_nyNodes-1,k);
      DivKgrad<BC,nv>::setFlux_x(K,Rho,P);
    }
  
  for (int vi=0;vi<nv;vi++)
    this->bc[vi].applyNeumannConditions(this->node,this->flux_x[vi], this->flux_y[vi], this->flux_z[vi]);
}

template<class BC, int nv>
void DivKgrad3d<BC,nv>::computeInterfaceK(const Vec& K)
{
  //ghost faces
  if (this->node.ghost_zOffSet)
    for (int j=0;j<this->local_nyNodes;j++)
      for (int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(-this->node.ghost_zOffSet,j,k);
          DivKgrad<BC,nv>::setKsz(K);
        }
  if (this->node.ghost_nzNodes -  this->node.ghost_zOffSet - this->node.local_nzNodes)
    for (int j=0;j<this->local_nyNodes;j++)
      for (int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          DivKgrad<BC,nv>::setKsz(K);
        }
  if (this->node.ghost_yOffSet)
    for (int i=0;i<this->local_nzNodes;i++)
      for(int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,-this->node.ghost_yOffSet,k);
          DivKgrad<BC,nv>::setKsy(K);
        }
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for (int i=0;i<this->local_nzNodes;i++)
      for(int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          DivKgrad<BC,nv>::setKsy(K);
        }
  if (this->node.ghost_xOffSet)
    for (int i=0;i<this->local_nzNodes;i++)
      for (int j=0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,-this->node.ghost_xOffSet);
          DivKgrad<BC,nv>::setKsx(K);
        }      
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int i=0;i<this->local_nzNodes;i++)
      for (int j=0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,this->local_nxNodes-1);
          DivKgrad<BC,nv>::setKsx(K);
        }
  
  //interior
  for (int i=0;i<this->local_nzNodes-1;i++)
    {
      for (int j=0;j<this->local_nyNodes-1;j++)
        {
          for (int k=0;k<this->local_nxNodes-1;k++)
            {
              this->node.localIndex(i,j,k);
              DivKgrad<BC,nv>::setKsx(K);
              DivKgrad<BC,nv>::setKsy(K);
              DivKgrad<BC,nv>::setKsz(K);
            }
          // last y and z fluxes of x row
          this->node.localIndex(i,j,this->local_nxNodes-1);
          DivKgrad<BC,nv>::setKsy(K);
          DivKgrad<BC,nv>::setKsz(K);
        }
      //last row of x and z fluxes in y plane
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          DivKgrad<BC,nv>::setKsx(K);
          DivKgrad<BC,nv>::setKsz(K);
        }
      // last z fluz in x row
      this->node.localIndex(i,this->local_nyNodes-1,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setKsz(K);
    }
  //top face of x and y fluxes
  for (int j=0;j<this->local_nyNodes-1;j++)
    {
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          DivKgrad<BC,nv>::setKsx(K);
          DivKgrad<BC,nv>::setKsy(K);
        }
      //last y flux in x row
      this->node.localIndex(this->local_nzNodes-1,j,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setKsy(K);
    }
  //last row of x fluxes in top face
  for (int k=0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nzNodes-1,this->local_nyNodes-1,k);
      DivKgrad<BC,nv>::setKsx(K);
    }
}

template<class BC, int nv>
void DivKgrad3d<BC,nv>::computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, Vec* Div)
{
  computeDivergence(K,Rho,P);
  for (int i=0;i<this->local_nzNodes;i++)
    for (int j=0;j<this->local_nyNodes;j++)
      for (int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,j,k);
          for (int vi=0;vi<nv;vi++)
            {
              Div[vi][this->node.center_noGhost] = 
                this->oneOverdx*(this->flux_x[vi][this->node.interLeft] - 
                           this->flux_x[vi][this->node.interRight])
                + this->oneOverdy*(this->flux_y[vi][this->node.interFront] - 
                             this->flux_y[vi][this->node.interBack])
                + this->oneOverdz*(this->flux_z[vi][this->node.interBottom] - 
                             this->flux_z[vi][this->node.interTop]);
            }
        }
}
  
template<class BC, int nv>
void DivKgrad3d<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho)
{
  //ghost faces
  if (this->node.ghost_zOffSet)
    for (int j=0;j<this->local_nyNodes;j++)
      for (int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(-this->node.ghost_zOffSet,j,k);
          DivKgrad<BC,nv>::setDFlux_z(K,Rho,P,DK,DRho);
        }
  if (this->node.ghost_nzNodes -  this->node.ghost_zOffSet - this->node.local_nzNodes)
    for (int j=0;j<this->local_nyNodes;j++)
      for (int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          DivKgrad<BC,nv>::setDFlux_z(K,Rho,P,DK,DRho);
        }
  if (this->node.ghost_yOffSet)
    for (int i=0;i<this->local_nzNodes;i++)
      for(int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,-this->node.ghost_yOffSet,k);
          DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
        }
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for (int i=0;i<this->local_nzNodes;i++)
      for(int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
        }
  if (this->node.ghost_xOffSet)
    for (int i=0;i<this->local_nzNodes;i++)
      for (int j=0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,-this->node.ghost_xOffSet);
          DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
        }      
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int i=0;i<this->local_nzNodes;i++)
      for (int j=0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,this->local_nxNodes-1);
          DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
        }
  
  for (int i=0;i<this->local_nzNodes-1;i++)
    {
      for (int j=0;j<this->local_nyNodes-1;j++)
        {
          for (int k=0;k<this->local_nxNodes-1;k++)
            {
              //interior
              this->node.localIndex(i,j,k);
              DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
              DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
              DivKgrad<BC,nv>::setDFlux_z(K,Rho,P,DK,DRho);
            }
          // last y and z fluxes of x row
          this->node.localIndex(i,j,this->local_nxNodes-1);
          DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
          DivKgrad<BC,nv>::setDFlux_z(K,Rho,P,DK,DRho);
        }
      //last row of x and z fluxes in y plane
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
          DivKgrad<BC,nv>::setDFlux_z(K,Rho,P,DK,DRho);
        }
      // last z fluz in x row
      this->node.localIndex(i,this->local_nyNodes-1,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setDFlux_z(K,Rho,P,DK,DRho);
    }
  //top face of x and y fluxes
  for (int j=0;j<this->local_nyNodes-1;j++)
    {
      for (int k=0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
          DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
        }
      //last y flux in x row
      this->node.localIndex(this->local_nzNodes-1,j,this->local_nxNodes-1);
      DivKgrad<BC,nv>::setDFlux_y(K,Rho,P,DK,DRho);
    }
  //last row of x fluxes in top face
  for (int k=0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nzNodes-1,this->local_nyNodes-1,k);
      DivKgrad<BC,nv>::setDFlux_x(K,Rho,P,DK,DRho);
    }
}

template<class BC, int nv>
void DivKgrad3d<BC,nv>::computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Rho,P,DK,DRho);
  for (int i=0;i<this->local_nzNodes;i++)
    for (int j=0;j<this->local_nyNodes;j++)
      for (int k=0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,j,k);
          for (int vi=0;vi<nv;vi++)
            {
              for (int vj=0;vj<nv;vj++)
                {
                  
                  DivJac[vi][vj][this->BOTTOM][this->node.center_noGhost] = this->oneOverdz*
                    this->Dflux_z_center[vi][vj][this->node.interBottom];
                  
                  DivJac[vi][vj][this->FRONT][this->node.center_noGhost] = this->oneOverdy*
                    this->Dflux_y_center[vi][vj][this->node.interFront];
                  
                  DivJac[vi][vj][this->LEFT][this->node.center_noGhost] = this->oneOverdx*
                    this->Dflux_x_center[vi][vj][this->node.interLeft];
                  
                  DivJac[vi][vj][this->CENTER][this->node.center_noGhost] = this->oneOverdx*
                    (this->Dflux_x_right[vi][vj][this->node.interLeft] 
                     -this->Dflux_x_center[vi][vj][this->node.interRight]) 
                    + this->oneOverdy*
                    (this->Dflux_y_back[vi][vj][this->node.interFront] 
                     -this->Dflux_y_center[vi][vj][this->node.interBack]) 
                    + this->oneOverdz*
                    (this->Dflux_z_top[vi][vj][this->node.interBottom] 
                     -this->Dflux_z_center[vi][vj][this->node.interTop]); 
                  
                  DivJac[vi][vj][this->RIGHT][this->node.center_noGhost]  = -this->oneOverdx*
                    this->Dflux_x_right[vi][vj][this->node.interRight]; 
                  
                  DivJac[vi][vj][this->BACK][this->node.center_noGhost]  = -this->oneOverdy*
                    this->Dflux_y_back[vi][vj][this->node.interBack]; 
                  
                  DivJac[vi][vj][this->TOP][this->node.center_noGhost]  = -this->oneOverdz*
                    this->Dflux_z_top[vi][vj][this->node.interTop]; 
                }
            }
        }
  for (int vi=0;vi<nv;vi++)
    {
      for (int vj=0;vj<nv;vj++)
        {
          this->bc[vi].applyNeumannDerivatives(this->node,DivJac[vi][vj],
                                         this->Dflux_x_center[vi][vj],this->Dflux_x_right[vi][vj],
                                         this->Dflux_y_center[vi][vj],this->Dflux_y_back[vi][vj],
                                         this->Dflux_z_center[vi][vj],this->Dflux_z_top[vi][vj],
                                         this->oneOverdx,this->oneOverdy,this->oneOverdz);
        }
    }  
}

}//Daetk
#endif
