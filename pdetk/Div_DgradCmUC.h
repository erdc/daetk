#ifndef DIV_DGRADCMUC_H
#define DIV_DGRADCMUC_H
#include "Definitions.h"
#include "PetscSecondOrderFd.h"


//  FIX DIVERGENCE JAC -- MUST APPLY CORRECT DERIVATIVES FOR NEUMAN BC'S

namespace Daetk 
{

using Petsc::SecondOrderFd;
typedef CMRVec<Vec> VecVec;
typedef CMRVec<VecVec> VecVecVec;
typedef CMRVec<VecVecVec> VecVecVecVec;

template<class BC, int nv>
class Div_DgradCmUC
{
public:
  enum JacEntry {CENTER, LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP};

  Div_DgradCmUC(BC* bcIn, SecondOrderFd& nodeIn);
  
  virtual ~Div_DgradCmUC();

  virtual void computeDivergence(const Vec* D, const Vec* U, const Vec* C, 
                                 Vec* Div)=0;

  virtual void computeDivergenceJac(const Vec* D, const Vec* U, const Vec* C,
                                    const VecVecVec& DD, const Vec* DU, 
                                    VecVecVecVec& DivJac)=0;  
  
  inline real Dx(const Vec& D, const Vec& C);
  inline real DDx_center(const Vec& DD, const Vec& C);
  inline real DDx_right(const Vec& DD, const Vec& C);

  inline real Dy(const Vec& D, const Vec& C);
  inline real DDy_center(const Vec& DD, const Vec& C);
  inline real DDy_back(const Vec& DD, const Vec& C);

  inline real Dz(const Vec& D, const Vec& C);
  inline real DDz_center(const Vec& DD, const Vec& C);
  inline real DDz_top(const Vec& DD, const Vec& C);

  inline real Fx(const Vec& U, const Vec& C);
  inline real DFx_center(const Vec& U);
  inline real DFx_right(const Vec& U);

  inline real Fy(const Vec& U, const Vec& C);
  inline real DFy_center(const Vec& U);
  inline real DFy_back(const Vec& U);

  inline real Fz(const Vec& U, const Vec& C);
  inline real DFz_center(const Vec& U);
  inline real DFy_top(const Vec& U);

  inline real DcDx(const Vec& C);
  inline real DcDy(const Vec& C);
  inline real DcDz(const Vec& C);

  inline void setFlux_x(const Vec* D, const Vec* U, const Vec* C);
  inline void setFlux_y(const Vec* D, const Vec* U, const Vec* C);
  inline void setFlux_z(const Vec* D, const Vec* U, const Vec* C);

  inline void setDFlux_x(const Vec* D, const Vec* U, const Vec* C, 
                         const VecVecVec& DD, const Vec* DU);
  inline void setDFlux_y(const Vec* D, const Vec* U, const Vec* C, 
                         const VecVecVec& DD, const Vec* DU);
  inline void setDFlux_z(const Vec* D, const Vec* U, const Vec* C, 
                         const VecVecVec& DD, const Vec* DU);

  const Vec& getFlux_x(int n=0);
  const Vec& getFlux_y(int n=0);
  const Vec& getFlux_z(int n=0);
  
protected:
  int i,j,k,vi,vj,nxNodes,nyNodes,nzNodes,local_nxNodes,local_nyNodes,local_nzNodes;
  real oneOverdx,oneOverdy,oneOverdz,gx,gy,gz;
  Vec flux_x[nv],flux_y[nv],flux_z[nv];
  VecVecVec Dflux_x_center,Dflux_x_right,Dflux_y_center,Dflux_y_back,Dflux_z_center,Dflux_z_top;
  BC* bc;
  SecondOrderFd& node;
};

template<class BC, int nv>
class Div_DgradCmUC1d : public Div_DgradCmUC<BC,nv>
{
public:
  Div_DgradCmUC1d(BC* bcIn, 
                  SecondOrderFd& nodeIn,
                  int nNodes, 
                  real g, 
                  real oneOverd);  

  virtual ~Div_DgradCmUC1d();
  
  virtual void computeDivergence(const Vec* D, const Vec* U, const Vec* C, 
                                 Vec* Div);
  
  virtual void computeDivergenceJac(const Vec* D, const Vec* U, const Vec* C,
                                    const VecVecVec& DD, const Vec* DU, 
                                    VecVecVecVec& DivJac);
}; 

template<class BC, int nv>
class Div_DgradCmUC2d : public Div_DgradCmUC<BC,nv>
{
public:
  Div_DgradCmUC2d(BC* bcIn,
                  SecondOrderFd& nodeIn,
                  int nxNodesIn,int nyNodesIn, 
                  real gxIn, real gyIn,
                  real oneOverdxIn, real oneOverdyIn);  
  
  virtual ~Div_DgradCmUC2d();

  virtual void computeDivergence(const Vec* D, const Vec* U, const Vec* C, 
                                 Vec* Div);
  virtual void computeDivergenceJac(const Vec* D, const Vec* U, const Vec* C, 
                                    const VecVecVec& DD, const Vec* DU, 
                                    VecVecVecVec& DivJac);
};  

template<class BC, int nv>
class Div_DgradCmUC3d : public Div_DgradCmUC<BC,nv>
{
public:
  Div_DgradCmUC3d(BC* bcIn,
                  SecondOrderFd& nodeIn,
                  int nxNodesIn,int nyNodesIn,int nzNodesIn, 
                  real gxIn, real gyIn, real gzIn,
                  real oneOverdxIn, real oneOverdyIn, real oneOverdzIn);
  virtual ~Div_DgradCmUC3d();
  virtual void computeDivergence(const Vec* D, const Vec* U, const Vec* C, 
                                 Vec* Div);
  
  virtual void computeDivergenceJac(const Vec* D, const Vec* U, const Vec* C, 
                                    const VecVecVec& DD, const Vec* DU, 
                                    VecVecVecVec& DivJac);
};  
  
template<class BC, int nv>
Div_DgradCmU<BC,nv>::Div_DgradCmUCC(BC* bcIn, SecondOrderFd& nodeIn):bc(bcIn),node(nodeIn){}
  
template<class BC, int nv>
virtual Div_DgradCmU<BC,nv>::~Div_DgradCmUC(){}

template<class BC, int nv>
inline real  Div_DgradCmU<BC,nv>::Dx(const Vec& D, const Vec& C)
{ 
  return 0.5*(D[node.center]+D[node.right]);
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DDx_center(const Vec& DD, const Vec& C)
{
  return 0.5*DD[node.center];
}

inline real Div_DgradCmU<BC,nv>::DDx_right(const Vec& DD, const Vec& C)
  {
    return 0.5*DD[node.right];
  }
  
template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::Dy(const Vec& D, const Vec& C)
{ 
  return 0.5*(D[node.center]+D[node.back]);
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DDy_center(const Vec& DD, const Vec& C)
{ 
  return 0.5*DD[node.center];
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DDy_back(const Vec& DD, const Vec& C)
{ 
  return 0.5*DD[node.back];
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::Dz(const Vec& D, const Vec& C)
{ 
  return 0.5*(D[node.center]+D[node.top]);
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DDz_center(const Vec& DD, const Vec& C)
{ 
  return 0.5*DD[node.center];
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DDz_top(const Vec& DD, const Vec& C)
{ 
  return 0.5*DD[node.top];
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::Fx(const Vec& U, const Vec& C)
{ 
  if (U[node.interRight] >= 0)
    return U[node.interRight]*C[node.center];
  else
    return U[node.interRight]*C[node.right];
}
  
template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DFx_center(const Vec& U)
{ 
  if (U[node.interRight] >= 0)
    return U[node.interRight];
  else
    return 0.0;
}
  
template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DFx_right(const Vec& U)
{ 
  if (U[node.interRight] >= 0)
    return 0.0;
  else
    return U[node.interRight];
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::Fy(const Vec& U, const Vec& C)
{ 
  if (U[node.interBack] >= 0)
    return U[node.interBack]*C[node.center];
  else
    return U[node.interBack]*C[node.back];
}
  
template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DFy_center(const Vec& U)
{ 
  if (U[node.interBack] >= 0)
    return U[node.interBack];
  else
    return 0.0;
}
  
template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DFy_back(const Vec& U)
{ 
  if (U[node.interBack] >= 0)
    return 0.0;
  else
    return U[node.interBack];
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::Fz(const Vec& U, const Vec& C)
{ 
  if (U[node.interTop] >= 0)
    return U[node.interTop]*C[node.center];
  else
    return U[node.interTop]*C[node.top];
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DFz_center(const Vec& U)
{ 
  if (U[node.interTop] >= 0)
    return U[node.interTop];
  else
    return 0.0;
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DFy_top(const Vec& U)
{ 
  if (U[node.interTop] >= 0)
    return 0.0;
  else
    return U[node.interTop];
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DcDx(const Vec& C)
{ 
  return oneOverdx*(C[node.right] - C[node.center]);
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DcDy(const Vec& C)
{ 
  return oneOverdy*(C[node.back]  - C[node.center]);
}

template<class BC, int nv>
inline real Div_DgradCmU<BC,nv>::DcDz(const Vec& C)
{ 
  return oneOverdz*(C[node.top]   - C[node.center]);
}

template<class BC, int nv>
inline void Div_DgradCmU<BC,nv>::setFlux_x(const Vec* D, const Vec* U, const Vec* C)
{                  
  for (vi=0;vi<nv;vi++)
    flux_x[vi][node.interRight] =-Dx(D[vi],C[vi]) * DcDx(C[vi]) + Fx(U[vi],C[vi]);
}

template<class BC, int nv>
inline void Div_DgradCmU<BC,nv>::setFlux_y(const Vec* D, const Vec* U, const Vec* C)
{
  for (vi=0;vi<nv;vi++)
    flux_y[vi][node.interBack] = -Dy(D[vi],C[vi]) * DcDy(C[vi]) + Fy(U[vi],C[vi]);
}

template<class BC, int nv>
inline void Div_DgradCmU<BC,nv>::setFlux_z(const Vec* D, const Vec* U, const Vec* C)
{
  for (vi=0;vi<nv;vi++)
    flux_z[vi][node.interTop] = -Dz(D[vi],C[vi]) * DcDz(C[vi]) + Fz(U[vi],C[vi]);
}

template<class BC, int nv>
inline void Div_DgradCmU<BC,nv>::setDFlux_x(const Vec* D, const Vec* U, const Vec* C, 
                                            const VecVecVec& DD, const Vec* DU)
{
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          //interior
          Dflux_x_center[vi][vj][node.interRight] = -DDx_center(DD[vi][vj],C[vi])*DcDx(C[vi]);
          
          Dflux_x_right[vi][vj][node.interRight]  = -DDx_right(DD[vi][vj],C[vi]) * DcDx(C[vi]);
          
          if (vi == vj)
            { 
                Dflux_x_center[vi][vj][node.interRight] += 
                  -Dx(D[vi],C[vi])*(-oneOverdx) 
                  +DFx_center(U[vi]);
                
                Dflux_x_right[vi][vj][node.interRight] +=
                  -Dx(D[vi],C[vi])*(oneOverdx)
                  +DFx_right(U[vi]);
            }
        }
    }
}

template<class BC, int nv>
inline void Div_DgradCmU<BC,nv>::setDFlux_y(const Vec* D, const Vec* U, const Vec* C, 
                                            const VecVecVec& DD, const Vec* DU)
{
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          //interior
          Dflux_y_center[vi][vj][node.interBack] = -DDy_center(DD[vi][vj],C[vi])*DcDy(C[vi]);
          
          Dflux_y_back[vi][vj][node.interBack]  = -DDy_back(DD[vi][vj],C[vi]) * DcDy(C[vi]);
          
            if (vi == vj)
              { 
                Dflux_y_center[vi][vj][node.interBack] += 
                  -Dy(D[vi],C[vi])*(-oneOverdy) 
                  +DFy_center(U[vi]);
                
                Dflux_y_back[vi][vj][node.interBack] +=
                  -Dy(D[vi],C[vi])*(oneOverdy)
                  +DFy_back(U[vi]);
              }
        }
    }
}
  
template<class BC, int nv>
inline void Div_DgradCmU<BC,nv>::setDFlux_z(const Vec* D, const Vec* U, const Vec* C, 
                                            const VecVecVec& DD, const Vec* DU)
{
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          //interior
          Dflux_z_center[vi][vj][node.interTop] = -DDz_center(DD[vi][vj],C[vi])*DcDz(C[vi]);
          
          Dflux_z_top[vi][vj][node.interTop]  = -DDz_top(DD[vi][vj],C[vi]) * DcDz(C[vi]);
          
          if (vi == vj)
            { 
              Dflux_z_center[vi][vj][node.interTop] += 
                -Dz(D[vi],C[vi])*(-oneOverdz) 
                +DFz_center(U[vi]);
              
              Dflux_z_top[vi][vj][node.interTop] +=
                -Dz(D[vi],C[vi])*(oneOverdz)
                +DFz_top(U[vi]);
            }
        }
    }
}

template<class BC, int nv>
const Vec& Div_DgradCmU<BC,nv>::getFlux_x(int n=0){return flux_x[n];}

template<class BC, int nv>
const Vec& Div_DgradCmU<BC,nv>::getFlux_y(int n=0){return flux_y[n];}

template<class BC, int nv>
const Vec& Div_DgradCmU<BC,nv>::getFlux_z(int n=0){return flux_z[n];}

template<class BC, int nv>
Div_DgradCmUC1d<BC,nv>::Div_DgradCmUC1d(BC* bcIn, 
                                        SecondOrderFd& nodeIn,
                                        int nNodes, 
                                        real g, 
                                        real oneOverd):
  Div_DgradCmUC<BC,nv>(bcIn,nodeIn)
{
  nxNodes=nNodes;
  local_nxNodes = nodeIn.local_nxNodes;
  
  gx=g;
  oneOverdx=oneOverd;
  
  Dflux_x_center.newsize(nv);
  Dflux_x_right.newsize(nv);
  for (vi=0;vi<nv;vi++)
    {
      flux_x[vi].newsize(Vec::LOCAL,local_nxNodes+1);
        flux_x[vi] = 0.0;
        Dflux_x_center[vi].newsize(nv);
        Dflux_x_right[vi].newsize(nv);
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
Div_DgradCmUC1d<BC,nv>::~Div_DgradCmUC1d(){}

template<class BC, int nv>
void computeDivergence(const Vec* D, const Vec* U, const Vec* C, 
                       Vec* Div)
{
  for (k= -node.ghost_xOffSet; k< (node.ghost_nxNodes-node.ghost_xOffSet-1); k++)
    {
      node.localIndex(k);
      setFlux_x(D,U,C);
    }
  
  for (vi=0; vi<nv; vi++)
    bc[vi].applyNeumannConditions(node,flux_x[vi]);
  
  for (k=0; k<local_nxNodes; k++)
    {
      node.localIndex(k);
      for (vi=0;vi<nv;vi++)
        {
          Div[vi][node.center_noGhost] = oneOverdx*
            (flux_x[vi][node.interLeft] - flux_x[vi][node.interRight]);
        }
    }
}  

template<class BC, int nv>
void Div_DgradCmUC1d<BC,nv>::computeDivergenceJac(const Vec* D, const Vec* U, const Vec* C,
                                                  const VecVecVec& DD, const Vec* DU, 
                                                  VecVecVecVec& DivJac)
{
  for (k=-node.ghost_xOffSet; k< node.ghost_nxNodes-node.ghost_xOffSet-1; k++)
    {
      node.localIndex(k);
      setDFlux_x(D,U,C,DD,DU);
    }
  
  
  for (k=0;k<local_nxNodes;k++)
    {
      node.localIndex(k);
      for (vi=0;vi<nv;vi++)
        {
          for (vj=0;vj<nv;vj++)
            {
              DivJac[vi][vj][LEFT][node.center_noGhost] = oneOverdx*
                Dflux_x_center[vi][vj][node.interLeft];
              
              DivJac[vi][vj][CENTER][node.center_noGhost] = oneOverdx*
                  (Dflux_x_right[vi][vj][node.interLeft] 
                   -Dflux_x_center[vi][vj][node.interRight]); 
              
              DivJac[vi][vj][RIGHT][node.center_noGhost]  = -oneOverdx*
                Dflux_x_right[vi][vj][node.interRight]; 
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
Div_DgradCmUC2d<BC,nv>::Div_DgradCmUC2d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn, 
             real gxIn, real gyIn,
             real oneOverdxIn, real oneOverdyIn):
  Div_DgradCmUC<BC,nv>(bcIn,nodeIn)
{
  nxNodes=nxNodesIn;
  local_nxNodes = nodeIn.local_nxNodes;
  gx=gxIn;
  oneOverdx = oneOverdxIn;
  
  nyNodes=nyNodesIn;
  local_nyNodes = nodeIn.local_nyNodes;
  gy=gyIn;
  oneOverdy = oneOverdyIn;
  
  Dflux_x_center.newsize(nv);
  Dflux_x_right.newsize(nv);
  Dflux_y_center.newsize(nv);
  Dflux_y_back.newsize(nv);
  for (vi=0;vi<nv;vi++)
    {
      flux_x[vi].newsize(Vec::LOCAL,local_nyNodes*(local_nxNodes+1));
      flux_y[vi].newsize(Vec::LOCAL,(local_nyNodes+1)*local_nxNodes);
      flux_x[vi] = 0.0;
      flux_y[vi] = 0.0;
      Dflux_x_center[vi].newsize(nv);
      Dflux_x_right[vi].newsize(nv);
      Dflux_y_center[vi].newsize(nv);
      Dflux_y_back[vi].newsize(nv);
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
Div_DgradCmUC2d<BC,nv>::~Div_DgradCmUC2d(){}

template<class BC, int nv>
void Div_DgradCmUC2d<BC,nv>::computeDivergence(const Vec* D, const Vec* U, const Vec* C, 
                                               Vec* Div)
{
  //bottom ghost line
  if (node.ghost_yOffSet)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(-node.ghost_yOffSet,k);
        setFlux_y(D,U,C);
      }
  //top ghost line
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(local_nyNodes-1,k);
        setFlux_y(D,U,C);
      }
  
  //left ghost line
  if (node.ghost_xOffSet)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,-node.ghost_xOffSet);
        setFlux_x(D,U,C);
      }
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,local_nxNodes-1);
        setFlux_x(D,U,C);
      }

  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(j,k);
          setFlux_x(D,U,C);
          setFlux_y(D,U,C);
            }
      //last y flux in row
      node.localIndex(j,k);
      setFlux_y(D,U,C);
    }
  //last row of x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setFlux_x(D,U,C);
    }
  
  for (vi=0;vi<nv;vi++)
    bc[vi].applyNeumannConditions(node,flux_x[vi], flux_y[vi]);
  
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
void Div_DgradCmUC2d<BC,nv>::computeDivergenceJac(const Vec* D, const Vec* U, const Vec* C, 
                                                  const VecVecVec& DD, const Vec* DU, 
                                                  VecVecVecVec& DivJac)
{
  //bottom ghost line
  if (node.ghost_yOffSet)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(-node.ghost_yOffSet,k);
        setDFlux_y(D,U,C,DD,DU);
      }
  //top ghost line
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(local_nyNodes-1,k);
        setDFlux_y(D,U,C,DD,DU);
      }
  
  //left ghost line
  if (node.ghost_xOffSet)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,-node.ghost_xOffSet);
        setDFlux_x(D,U,C,DD,DU);
      }
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,local_nxNodes-1);
        setDFlux_x(D,U,C,DD,DU);
      }
  
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(j,k);
          setDFlux_x(D,U,C,DD,DU);
          setDFlux_y(D,U,C,DD,DU);
        }
      //last y flux in row
      node.localIndex(j,k);
      setDFlux_y(D,U,C,DD,DU);
    }
  //last row of x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setDFlux_x(D,U,C,DD,DU);
    }
  
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
Div_DgradCmUC3d<BC,nv>::Div_DgradCmUC3d(BC* bcIn,
                                        SecondOrderFd& nodeIn,
                                        int nxNodesIn,int nyNodesIn,int nzNodesIn, 
                                        real gxIn, real gyIn, real gzIn,
                                        real oneOverdxIn, real oneOverdyIn, real oneOverdzIn):
  Div_DgradCmUC<BC,nv>(bcIn,nodeIn)
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
  
  Dflux_x_center.newsize(nv);
  Dflux_x_right.newsize(nv);
  Dflux_y_center.newsize(nv);
  Dflux_y_back.newsize(nv);
  Dflux_z_center.newsize(nv);
  Dflux_z_top.newsize(nv);
  for (vi=0;vi<nv;vi++)
    {
      flux_x[vi].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
      flux_y[vi].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
      flux_z[vi].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
      flux_x[vi] = 0.0;
      flux_y[vi] = 0.0;
      flux_z[vi] = 0.0;
      Dflux_x_center[vi].newsize(nv);
      Dflux_x_right[vi].newsize(nv);
      Dflux_y_center[vi].newsize(nv);
      Dflux_y_back[vi].newsize(nv);
      Dflux_z_center[vi].newsize(nv);
      Dflux_z_top[vi].newsize(nv);
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
Div_DgradCmUC3d<BC,nv>::~Div_DgradCmUC3d(){}

template<class BC, int nv>
void Div_DgradCmUC3d<BC,nv>::computeDivergence(const Vec* D, const Vec* U, const Vec* C, 
                                 Vec* Div)
    {
      //ghost faces
      if (node.ghost_zOffSet)
        for (j=0;j<local_nyNodes;j++)
          for (k=0;k<local_nxNodes;k++)
            {
              node.localIndex(-node.ghost_zOffSet,j,k);
              setFlux_z(D,U,C);
            }
      if (node.ghost_nzNodes -  node.ghost_zOffSet - node.local_nzNodes)
        for (j=0;j<local_nyNodes;j++)
          for (k=0;k<local_nxNodes;k++)
            {
                node.localIndex(local_nzNodes-1,j,k);
                setFlux_z(D,U,C);
            }
      if (node.ghost_yOffSet)
        for (i=0;i<local_nzNodes;i++)
          for(k=0;k<local_nxNodes;k++)
            {
              node.localIndex(i,-node.ghost_yOffSet,k);
              setFlux_y(D,U,C);
            }
      if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
        for (i=0;i<local_nzNodes;i++)
          for(k=0;k<local_nxNodes;k++)
            {
              node.localIndex(i,local_nyNodes-1,k);
              setFlux_y(D,U,C);
            }
      if (node.ghost_xOffSet)
        for (i=0;i<local_nzNodes;i++)
          for (j=0;j<local_nyNodes;j++)
            {
              node.localIndex(i,j,-node.ghost_xOffSet);
              setFlux_x(D,U,C);
            }      
      if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
        for (i=0;i<local_nzNodes;i++)
          for (j=0;j<local_nyNodes;j++)
            {
              node.localIndex(i,j,local_nxNodes-1);
              setFlux_x(D,U,C);
            }

      //interior
      for (i=0;i<local_nzNodes-1;i++)
        {
          for (j=0;j<local_nyNodes-1;j++)
            {
              for (k=0;k<local_nxNodes-1;k++)
                {
                  node.localIndex(i,j,k);
                  setFlux_x(D,U,C);
                  setFlux_y(D,U,C);
                  setFlux_z(D,U,C);
                }
              // last y and z fluxes of x row
              node.localIndex(i,j,k);
              setFlux_y(D,U,C);
              setFlux_z(D,U,C);
            }
          //last row of x and z fluxes in y plane
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(i,j,k);
              setFlux_x(D,U,C);
              setFlux_z(D,U,C);
            }
          // last z fluz in x row
          node.localIndex(i,j,k);
          setFlux_z(D,U,C);
        }
      //top face of x and y fluxes
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(i,j,k);
              setFlux_x(D,U,C);
              setFlux_y(D,U,C);
            }
          //last y flux in x row
          node.localIndex(i,j,k);
          setFlux_y(D,U,C);
        }
      //last row of x fluxes in top face
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setFlux_x(D,U,C);
        }

//        std::cout<<"in bc's"<<std::endl;
      for (vi=0;vi<nv;vi++)
        bc[vi].applyNeumannConditions(node,flux_x[vi], flux_y[vi], flux_z[vi]);
      
//        std::cout<<"in div"<<std::endl;


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
//        std::cout<<"out of dkg"<<std::endl;
    }


template<class BC, int nv>
void Div_DgradCmUC3d<BC,nv>::computeDivergenceJac(const Vec* D, const Vec* U, const Vec* C, 
                                    const VecVecVec& DD, const Vec* DU, 
                                    VecVecVecVec& DivJac)
  {
    //ghost faces
    if (node.ghost_zOffSet)
      for (j=0;j<local_nyNodes;j++)
        for (k=0;k<local_nxNodes;k++)
          {
            node.localIndex(-node.ghost_zOffSet,j,k);
            setDFlux_z(D,U,C,DD,DU);
          }
    if (node.ghost_nzNodes -  node.ghost_zOffSet - node.local_nzNodes)
      for (j=0;j<local_nyNodes;j++)
        for (k=0;k<local_nxNodes;k++)
          {
            node.localIndex(local_nzNodes-1,j,k);
            setDFlux_z(D,U,C,DD,DU);
          }
    if (node.ghost_yOffSet)
      for (i=0;i<local_nzNodes;i++)
        for(k=0;k<local_nxNodes;k++)
          {
            node.localIndex(i,-node.ghost_yOffSet,k);
            setDFlux_y(D,U,C,DD,DU);
          }
    if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
      for (i=0;i<local_nzNodes;i++)
        for(k=0;k<local_nxNodes;k++)
          {
            node.localIndex(i,local_nyNodes-1,k);
            setDFlux_y(D,U,C,DD,DU);
          }
    if (node.ghost_xOffSet)
      for (i=0;i<local_nzNodes;i++)
        for (j=0;j<local_nyNodes;j++)
          {
            node.localIndex(i,j,-node.ghost_xOffSet);
            setDFlux_x(D,U,C,DD,DU);
          }      
    if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
      for (i=0;i<local_nzNodes;i++)
        for (j=0;j<local_nyNodes;j++)
          {
            node.localIndex(i,j,local_nxNodes-1);
              setDFlux_x(D,U,C,DD,DU);
          }
    
      for (i=0;i<local_nzNodes-1;i++)
        {
          for (j=0;j<local_nyNodes-1;j++)
            {
              for (k=0;k<local_nxNodes-1;k++)
                {
                //interior
                  node.localIndex(i,j,k);
                  setDFlux_x(D,U,C,DD,DU);
                  setDFlux_y(D,U,C,DD,DU);
                  setDFlux_z(D,U,C,DD,DU);
                }
              // last y and z fluxes of x row
              node.localIndex(i,j,k);
              setDFlux_y(D,U,C,DD,DU);
              setDFlux_z(D,U,C,DD,DU);
            }
          //last row of x and z fluxes in y plane
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(i,j,k);
              setDFlux_x(D,U,C,DD,DU);
              setDFlux_z(D,U,C,DD,DU);
            }
          // last z fluz in x row
          node.localIndex(i,j,k);
          setDFlux_z(D,U,C,DD,DU);
        }
      //top face of x and y fluxes
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(i,j,k);
              setDFlux_x(D,U,C,DD,DU);
              setDFlux_y(D,U,C,DD,DU);
            }
          //last y flux in x row
          node.localIndex(i,j,k);
          setDFlux_y(D,U,C,DD,DU);
        }
      //last row of x fluxes in top face
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setDFlux_x(D,U,C,DD,DU);
        }
      

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


