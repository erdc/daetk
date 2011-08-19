#ifndef DIVKGRAD_PICARD_H
#define DIVKGRAD_PICARD_H
#include "PetscSecondOrderFd.h"
#include <vector>

//  FIX DIVERGENCE JAC -- MUST APPLY CORRECT DERIVATIVES FOR NEUMAN BC'S

namespace Daetk 
{

using Petsc::SecondOrderFd;


//  FIX DIVERGENCE JAC -- MUST APPLY CORRECT DERIVATIVES FOR NEUMAN BC'S

template<class BC, int nv>
class DivKgradPicard
{
public:
  enum JacEntry {CENTER, LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP};
  DivKgradPicard(BC* bcIn, SecondOrderFd& nodeIn):bc(bcIn),node(nodeIn),
    picardMultiple(1.0){}
  
  virtual ~DivKgradPicard(){}

  void usePicardIteration(bool usePic)
    { 
      if (usePic) 
	picardMultiple=1.0;
      else
	picardMultiple=0.0;
    }
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div)=0;

  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)=0;  
  
  inline real Kx(const Vec& K, const Vec& P)
    { 
      
      return 0.5*(K[node.center]+K[node.right]);

//        return (P(node.center) >= P(node.right)) ? K(node.center) : K(node.right);

//        if ((DpDx(P) + gx) > 0)
//          return K(node.right);
//        else
//          return K(node.center);

//        if (K[node.center]==0.0 || K[node.right]==0.0)
//          return 0.0;
//        else
//          return 2.0/(1.0/K[node.center]  + 1.0/K[node.right]);
      
//        if (P(node.center) == P(node.right))
//          return sqrt(K(node.center));
//        else
//          return fabs(K(node.right)-K(node.center));
    }
  inline real DKx_center(const Vec& DK, const Vec& P)
    {
      return 0.5*DK[node.center]*picardMultiple;
      //return (P[node.center] >= P[node.right]) ? DK[node.center] : 0.0;
    }
  inline real DKx_right(const Vec& DK, const Vec& P)
    {
      return 0.5*DK[node.right]*picardMultiple;
      //return (P[node.center] >= P[node.right]) ? 0.0 : DK[node.right];
    }

  inline real Ky(const Vec& K, const Vec& P)
    { 
      return 0.5*(K[node.center]+K[node.back]);
      //return (P[node.center] >= P[node.back]) ? K[node.center] : K[node.back];
    }
  inline real DKy_center(const Vec& DK, const Vec& P)
    { 
      return 0.5*DK[node.center]*picardMultiple;
      //return (P[node.center] >= P[node.back]) ? DK[node.center] : 0.0;
    }
  inline real DKy_back(const Vec& DK, const Vec& P)
    { 
      return 0.5*DK[node.back]*picardMultiple;
      //return (P[node.center] >= P[node.back]) ? 0.0 : DK[node.back];
    }

  inline real Kz(const Vec& K, const Vec& P)
    { 
      return 0.5*(K[node.center]+K[node.top]);
      //return (P[node.center] >= P[node.top]) ? K[node.center] : K[node.top];
    }
  inline real DKz_center(const Vec& DK, const Vec& P)
    { 
      return 0.5*DK[node.center]*picardMultiple;
      //return (P[node.center] >= P[node.top]) ? DK[node.center] : 0.0;
    }
  inline real DKz_top(const Vec& DK, const Vec& P)
    { 
      return 0.5*DK[node.top]*picardMultiple;
      //return (P[node.center] >= P[node.top]) ? 0.0 : DK[node.top];
    }

  inline real Rhox(const Vec& Rho)
    { return 0.5*(Rho[node.center]+Rho[node.right]);}
  inline real DRhox_center(const Vec& DRho)
    { return 0.5*DRho[node.center]*picardMultiple;}
  inline real DRhox_right(const Vec& DRho)
    { return 0.5*DRho[node.right]*picardMultiple;}

  inline real Rhoy(const Vec& Rho)
    { return 0.5*(Rho[node.center]+Rho[node.back]);}
  inline real DRhoy_center(const Vec& DRho)
    { return 0.5*DRho[node.center]*picardMultiple;}
  inline real DRhoy_back(const Vec& DRho)
    { return 0.5*DRho[node.back]*picardMultiple;}

  inline real Rhoz(const Vec& Rho)
    { return 0.5*(Rho[node.center]+Rho[node.top]);}
  inline real DRhoz_center(const Vec& DRho)
    { return 0.5*DRho[node.center]*picardMultiple;}
  inline real DRhoz_top(const Vec& DRho)
    { return 0.5*DRho[node.top]*picardMultiple;}

  inline real DpDx(const Vec& P)
    { return oneOverdx*(P[node.right] - P[node.center]);}
  inline real DpDy(const Vec& P)
  { 
//      cout<<"DpDy "<<oneOverdy<<'\t'<<P[node.back]<<'\t'<<P[node.center]<<endl;
    return oneOverdy*(P[node.back]  - P[node.center]);
  }
  inline real DpDz(const Vec& P)
    { return oneOverdz*(P[node.top]   - P[node.center]);}


  inline void setFlux_x(const Vec* K, const Vec* Rho, const Vec* P)
  {                  
    //std::cout<<node.center<<'\t'<<node.interRight<<'\t'<<node.right<<std::endl<<std::flush;
    //std::cout<<flux_x[0][node.interRight]<<std::endl<<std::flush;
    //std::cout <<Rhox(Rho[0])<<std::endl<<std::flush;
    for (vi=0;vi<nv;vi++)
      flux_x[vi][node.interRight] = -Rhox(Rho[vi])*
        Kx(K[vi],P[vi]) * (DpDx(P[vi]) + Rhox(Rho[vi])*gx);
  }
  
  inline void setFlux_y(const Vec* K, const Vec* Rho, const Vec* P)
  {
    for (vi=0;vi<nv;vi++)
      flux_y[vi][node.interBack]  = -Rhoy(Rho[vi])*
        Ky(K[vi],P[vi]) * (DpDy(P[vi]) + Rhoy(Rho[vi])*gy);
//      cout<<"y "<<node.center<<'\t'<<node.interBack<<'\t'<<node.back<<endl<<flush;
//      cout<<flux_y[0][node.interBack]<<endl<<flush;
//      cout <<Rhoy(Rho[0])<<'\t'<<DpDy(P[vi])<<endl<<flush;
  }
  
  inline void setFlux_z(const Vec* K, const Vec* Rho, const Vec* P)
  {
//      cout<<"z "<<node.center<<'\t'<<node.interTop<<'\t'<<node.top<<endl<<flush;
//      cout<<flux_z[0][node.interTop]<<endl<<flush;
//      cout <<Rhoz(Rho[0])<<endl<<flush;
    for (vi=0;vi<nv;vi++)
      flux_z[vi][node.interTop]   = -Rhoz(Rho[vi])*
        Kz(K[vi],P[vi]) * (DpDz(P[vi]) + Rhoz(Rho[vi])*gz);
  }

  inline void setDFlux_x(const Vec* K, const Vec* Rho, const Vec* P, 
                         const VecVecVec& DK, const Vec* DRho)
  {
    for (vi=0;vi<nv;vi++)
      {
        for (vj=0;vj<nv;vj++)
          {
            //interior
            Dflux_x_center[vi][vj][node.interRight] = -Rhox(Rho[vi])*
              DKx_center(DK[vi][vj],P[vi])*
              (DpDx(P[vi]) + Rhox(Rho[vi])*gx);
            
            Dflux_x_right[vi][vj][node.interRight] = -Rhox(Rho[vi])*
              DKx_right(DK[vi][vj],P[vi]) * 
              (DpDx(P[vi]) + Rhox(Rho[vi])*gx);
            
            if (vi == vj)
              { 
                Dflux_x_center[vi][vj][node.interRight] += 
                  -Rhox(Rho[vi])*Kx(K[vi],P[vi])*
                  (-oneOverdx + DRhox_center(DRho[vi])*gx) 
                  -DRhox_center(DRho[vi])*Kx(K[vi],P[vi])*
                  (DpDx(P[vi]) + Rhox(Rho[vi])*gx);
                
                Dflux_x_right[vi][vj][node.interRight] +=
                  -Rhox(Rho[vi])*Kx(K[vi],P[vi])*
                  (oneOverdx + DRhox_right(DRho[vi])*gx) 
                  -DRhox_right(DRho[vi])*Kx(K[vi],P[vi])*
                  (DpDx(P[vi]) + Rhox(Rho[vi])*gx);
              }
          }
      }
  }

  inline void setDFlux_y(const Vec* K, const Vec* Rho, const Vec* P, 
                         const VecVecVec& DK, const Vec* DRho)
  {
    for (vi=0;vi<nv;vi++)
      {
        for (vj=0;vj<nv;vj++)
          {
            //interior
            Dflux_y_center[vi][vj][node.interBack] = -Rhoy(Rho[vi])*
              DKy_center(DK[vi][vj],P[vi])*
              (DpDy(P[vi]) + Rhoy(Rho[vi])*gy);
            
            Dflux_y_back[vi][vj][node.interBack] = -Rhoy(Rho[vi])*
              DKy_back(DK[vi][vj],P[vi]) * 
              (DpDy(P[vi]) + Rhoy(Rho[vi])*gy);
            
            if (vi == vj)
              { 
                Dflux_y_center[vi][vj][node.interBack] += 
                  -Rhoy(Rho[vi])*Ky(K[vi],P[vi])*
                  (-oneOverdy + DRhoy_center(DRho[vi])*gy) 
                  -DRhoy_center(DRho[vi])*Ky(K[vi],P[vi])*
                  (DpDy(P[vi]) + Rhoy(Rho[vi])*gy);
                
                Dflux_y_back[vi][vj][node.interBack] +=
                  -Rhoy(Rho[vi])*Ky(K[vi],P[vi])*
                  (oneOverdy + DRhoy_back(DRho[vi])*gy) 
                  -DRhoy_back(DRho[vi])*Ky(K[vi],P[vi])*
                  (DpDy(P[vi]) + Rhoy(Rho[vi])*gy);
              }
          }
      }
  }

  inline void setDFlux_z(const Vec* K, const Vec* Rho, const Vec* P, 
                         const VecVecVec& DK, const Vec* DRho)
  {
    for (vi=0;vi<nv;vi++)
      {
        for (vj=0;vj<nv;vj++)
          {
            //interior
            Dflux_z_center[vi][vj][node.interTop] = -Rhoz(Rho[vi])*
              DKz_center(DK[vi][vj],P[vi])*
              (DpDz(P[vi]) + Rhoz(Rho[vi])*gz);
            
            Dflux_z_top[vi][vj][node.interTop] = -Rhoz(Rho[vi])*
              DKz_top(DK[vi][vj],P[vi]) * 
              (DpDz(P[vi]) + Rhoz(Rho[vi])*gz);
            
            if (vi == vj)
              { 
                Dflux_z_center[vi][vj][node.interTop] += 
                  -Rhoz(Rho[vi])*Kz(K[vi],P[vi])*
                  (-oneOverdz + DRhoz_center(DRho[vi])*gz) 
                  -DRhoz_center(DRho[vi])*Kz(K[vi],P[vi])*
                  (DpDz(P[vi]) + Rhoz(Rho[vi])*gz);
                
                Dflux_z_top[vi][vj][node.interTop] +=
                  -Rhoz(Rho[vi])*Kz(K[vi],P[vi])*
                  (oneOverdz + DRhoz_top(DRho[vi])*gz) 
                  -DRhoz_top(DRho[vi])*Kz(K[vi],P[vi])*
                  (DpDz(P[vi]) + Rhoz(Rho[vi])*gz);
              }
          }
      }
  }
protected:
  int nxNodes,nyNodes,nzNodes,local_nxNodes,local_nyNodes,local_nzNodes;
  real oneOverdx,oneOverdy,oneOverdz,gx,gy,gz;
  Vec flux_x[nv],flux_y[nv],flux_z[nv];
  VecVecVec Dflux_x_center,Dflux_x_right,Dflux_y_center,Dflux_y_back,Dflux_z_center,Dflux_z_top;
  BC* bc;
  SecondOrderFd& node;
  //mwf added for picard iteration
  real picardMultiple;
};

template<class BC, int nv>
class DivKgradPicard1d : public DivKgradPicard<BC,nv>
{
public:
  DivKgradPicard1d(BC* bcIn, 
             SecondOrderFd& nodeIn,
             int nNodes, 
             real g, 
             real oneOverd):
    DivKgradPicard<BC,nv>(bcIn,nodeIn)
  {
    nxNodes=nNodes;
    local_nxNodes = nodeIn.local_nxNodes;

    gx=g;
    oneOverdx=oneOverd;

    Dflux_x_center.resize(nv);
    Dflux_x_right.resize(nv);
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

  virtual ~DivKgradPicard1d(){}
  
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div)
  {
    for (k= -node.ghost_xOffSet; k< (node.ghost_nxNodes-node.ghost_xOffSet-1); k++)
      {
        node.localIndex(k);
        setFlux_x(K,Rho,P);
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
  
  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
  {
    for (k=-node.ghost_xOffSet; k< node.ghost_nxNodes-node.ghost_xOffSet-1; k++)
      {
        node.localIndex(k);
        setDFlux_x(K,Rho,P,DK,DRho);
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
  
}; 

template<class BC, int nv>
class DivKgradPicard2d : public DivKgradPicard<BC,nv>
{
public:
  DivKgradPicard2d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn, 
             real gxIn, real gyIn,
             real oneOverdxIn, real oneOverdyIn):
    DivKgradPicard<BC,nv>(bcIn,nodeIn)
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
  virtual ~DivKgradPicard2d(){}
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div)
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
  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
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
};  

template<class BC, int nv>
class DivKgradPicard3d : public DivKgradPicard<BC,nv>
{
public:
  DivKgradPicard3d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn,int nzNodesIn, 
             real gxIn, real gyIn, real gzIn,
             real oneOverdxIn, real oneOverdyIn, real oneOverdzIn):
    DivKgradPicard<BC,nv>(bcIn,nodeIn)
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
  virtual ~DivKgradPicard3d(){}
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div)
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
  
  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
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
};  

}//Daetk
#endif









