#ifndef DIV_KGRADPLUSFMM_H
#define DIV_KGRADPLUSFMM_H
#include "Definitions.h"
#include "PetscSecondOrderFd.h"
#include "Divergence.h"
#include "Div_KgradPlusFMM.h"

#include <vector>

//  FIX DIVERGENCE JAC -- MUST APPLY CORRECT DERIVATIVES FOR NEUMAN BC'S

namespace Daetk 
{
template<class FLUX_K, class FLUX_F, class BC, int nv>
class Div_KgradPlusFMM : public Divergence
{
public:
  enum JacEntry {CENTER, LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP};
  Div_KgradPlusFMM(BC* bcIn, SecondOrderFd& nodeIn);
  
  virtual ~Div_KgradPlusFMM();

  virtual void computeDivergence(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P,  
                                 Vec* Div)=0;

  virtual void computeDivergenceJac(const Vec* K, 
                                    const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                    const Vec* P, 
                                    const VecVecVec& DK, 
                                    const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz, 
                                    VecVecVecVec& DivJac)=0;  
  virtual void computeDivergence(const Vec* K, 
                                 const Vec* Fx,const Vec* Fy,const Vec* Fz, 
                                 const Vec* P)=0;

  virtual void computeDivergenceJac(const Vec* K, 
                                    const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                    const Vec* P,
                                    const VecVecVec& DK, 
                                    const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz)=0;  
  

  inline void setLocal_dy_dxi(Vec& local_dy_dxiIn){local_dy_dxi=&local_dy_dxiIn;}

  virtual void computeInterfaceK(const Vec& K)=0;
  inline void setKsx(const Vec& K);
  inline void setKsy(const Vec& K);
  inline void setKsz(const Vec& K);

  virtual void computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec& Fz)=0;
  inline void setFsx(const Vec& F);
  inline void setFsy(const Vec& F);
  inline void setFsz(const Vec& F);

  inline real ilKx(const Vec& K, const real& grad);
  inline real ilDKx_center(const Vec& DK, const real& grad);
  inline real ilDKx_right(const Vec& DK, const real& grad);

  inline real ilKy(const Vec& K, const real& grad);
  inline real ilDKy_center(const Vec& DK, const real& grad);
  inline real ilDKy_back(const Vec& DK, const real& grad);

  inline real ilKz(const Vec& K, const real& grad);
  inline real ilDKz_center(const Vec& DK, const real& grad);
  inline real ilDKz_top(const Vec& DK, const real& grad);

  inline real ilFx(const Vec& F, const Vec& P);
  inline real ilDFx_center(const Vec& DF, const Vec& F, const Vec& P);
  inline real ilDFx_right(const Vec& DF, const Vec& F, const Vec& P);

  inline real ilFy(const Vec& F);
  inline real ilDFy_center(const Vec& DF);
  inline real ilDFy_back(const Vec& DF);

  inline real ilFz(const Vec& F);
  inline real ilDFz_center(const Vec& DF);
  inline real ilDFz_top(const Vec& DF);

  inline real DpDx(const Vec& P);
  inline real DpDy(const Vec& P);
  inline real DpDz(const Vec& P);

  inline void setFlux_x(int v,const Vec& K, const Vec& F, const Vec& P);

  inline void setFlux_x(const Vec* K, const Vec* F, const Vec* P);
  inline void setFlux_y(const Vec* K, const Vec* F, const Vec* P);
  inline void setFlux_z(const Vec* K, const Vec* F, const Vec* P);

  inline void setDFlux_x(const Vec* K, const Vec* F, const Vec* P, 
                         const VecVecVec& DK, const VecVecVec& DF);
  inline void setDFlux_y(const Vec* K, const Vec* F, const Vec* P, 
                         const VecVecVec& DK, const VecVecVec& DF);
  inline void setDFlux_z(const Vec* K, const Vec* F, const Vec* P, 
                         const VecVecVec& DK, const VecVecVec& DF);


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
  real oneOverdx,oneOverdy,oneOverdz; 
  Vec flux_x[nv],flux_y[nv],flux_z[nv],Ksx,Ksy,Ksz,Fsx,Fsy,Fsz;
  Vec* local_dy_dxi;
  VecVecVec Dflux_x_center,Dflux_x_right,Dflux_y_center,Dflux_y_back,Dflux_z_center,Dflux_z_top;
  BC* bc;
  SecondOrderFd& node;
  FLUX_K flk;
  FLUX_F flf;
};

template<class FLUX_K, class FLUX_F, class BC, int nv>
class Div_KgradPlusFMM1d : public Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>
{
public:
  Div_KgradPlusFMM1d(BC* bcIn, 
             SecondOrderFd& nodeIn,
             int nNodes, 
             real oneOverd);

  virtual ~Div_KgradPlusFMM1d();
  
  virtual void computeDivergence(const Vec* K, 
                                 const Vec* Fx,const Vec* Fy,const Vec* Fz, 
                                 const Vec* P, 
                                 Vec* Div);

  virtual void computeDivergenceJac(const Vec* K, 
                                    const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                    const Vec* P,
                                    const VecVecVec& DK, 
                                    const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz, 
                                    VecVecVecVec& DivJac);  

  virtual void computeDivergence(const Vec* K, 
                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                 const Vec* P);
  
  virtual void computeDivergenceJac(const Vec* K, 
                                    const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                    const Vec* P,
                                    const VecVecVec& DK, 
                                    const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz);
  
  virtual void computeInterfaceK(const Vec& K);
  virtual void computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec& Fz);
}; 

template<class FLUX_K, class FLUX_F, class BC, int nv>
class Div_KgradPlusFMM2d : public Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>
{
public:
  Div_KgradPlusFMM2d(BC* bcIn,
                   SecondOrderFd& nodeIn,
                   int nxNodesIn,int nyNodesIn,
                   real oneOverdxIn, real oneOverdyIn);
  virtual ~Div_KgradPlusFMM2d();
  virtual void computeDivergence(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P, 
                                 Vec* Div);

  virtual void computeDivergenceJac(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P,
                                    const VecVecVec& DK, const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz, 
                                    VecVecVecVec& DivJac);  
  virtual void computeDivergence(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P);
  virtual void computeDivergenceJac(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P, 
                                    const VecVecVec& DK, const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz);
  virtual void computeInterfaceK(const Vec& K);
  virtual void computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec& Fz);
};  

template<class FLUX_K, class FLUX_F, class BC, int nv>
class Div_KgradPlusFMM3d : public Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>
{
public:
  Div_KgradPlusFMM3d(BC* bcIn,
                   SecondOrderFd& nodeIn,
                   int nxNodesIn,int nyNodesIn,int nzNodesIn,
                   real oneOverdxIn, real oneOverdyIn, real oneOverdzIn);
  virtual ~Div_KgradPlusFMM3d();
  virtual void computeDivergence(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P, 
                                 Vec* Div);

  virtual void computeDivergenceJac(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P,
                                    const VecVecVec& DK, const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz, 
                                    VecVecVecVec& DivJac);  
  virtual void computeDivergence(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P);
  virtual void computeDivergenceJac(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P, 
                                    const VecVecVec& DK, const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz);
  virtual void computeInterfaceK(const Vec& K);
  virtual void computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec& Fz);
};  
  
template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::Div_KgradPlusFMM(BC* bcIn, SecondOrderFd& nodeIn):
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
  Fsx.newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
  Fsy.newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
  Fsz.newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
  Fsx = 0.0;
  Fsy = 0.0;
  Fsz = 0.0;
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::~Div_KgradPlusFMM(){}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setKsx(const Vec& K)
{
  Ksx[node.interRight] = 2.0 / ( 1.0/K[node.center] + 1.0/K[node.right]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setKsy(const Vec& K)
{
  Ksy[node.interBack] = 2.0 / ( 1.0/K[node.center] + 1.0/K[node.back]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setKsz(const Vec& K)
{
  Ksz[node.interTop] = 2.0 / ( 1.0/K[node.center] + 1.0/K[node.top]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFsx(const Vec& F)
{
  Fsx[node.interRight] = 2.0 / ( 1.0/F[node.center] + 1.0/F[node.right]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFsy(const Vec& F)
{
  Fsy[node.interBack] = 2.0 / ( 1.0/F[node.center] + 1.0/F[node.back]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFsz(const Vec& F)
{
  Fsz[node.interTop] = 2.0 / ( 1.0/F[node.center] + 1.0/F[node.top]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilKx(const Vec& K, const real& grad )
{ 
  return 0.5*(K[node.center]+K[node.right])*Ksx[node.interRight];
}


template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDKx_center(const Vec& DK, const real& grad)
{
  return 0.5*DK[node.center]*Ksx[node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDKx_right(const Vec& DK, const real& grad)
{
  return 0.5*DK[node.right]*Ksx[node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilKy(const Vec& K, const real& grad)
{ 
  return 0.5*(K[node.center]+K[node.back])*Ksy[node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDKy_center(const Vec& DK,const real& grad)
{ 
  return 0.5*DK[node.center]*Ksy[node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDKy_back(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[node.back]*Ksy[node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilKz(const Vec& K, const real& grad)
{ 
  return 0.5*(K[node.center]+K[node.top])*Ksz[node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDKz_center(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[node.center]*Ksz[node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDKz_top(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[node.top]*Ksz[node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilFx(const Vec& F, const Vec& P)
{ 
//   std::cout<<"flf.flux "<<flf.flux(F[node.center],F[node.right],
//                                 P[node.center],P[node.right])<<'\t'<<F[node.center]<<std::endl;
  return flf.flux(F[node.center],F[node.right],
                  P[node.center],P[node.right])*Fsx[node.interRight];
  //upwind by hand for v>0
  //return F[node.center]*Fsx[node.interRight];;
  //return 0.5*(F[node.center]+F[node.right])*Fsx[node.interRight];;
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDFx_center(const Vec& DF, const Vec& F, const Vec& P)
{ 
//   std::cout<<"flf.Dflux_l "<<flf.Dflux_l(DF[node.center],DF[node.right],
//                    F[node.center],F[node.right],
//                                  P[node.center],P[node.right])<<'\t'<<DF[node.center]<<std::endl;
  return flf.Dflux_l(DF[node.center],DF[node.right],
                   F[node.center],F[node.right],
                   P[node.center],P[node.right])*Fsx[node.interRight];
  //upwind by hand for v>0
  //return DF[node.center]*Fsx[node.interRight];
  //return 0.5*DF[node.center]*Fsx[node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDFx_right(const Vec& DF, const Vec& F, const Vec& P)
{ 
//   std::cout<<"flf.Dflux_r "<<flf.Dflux_r(DF[node.center],DF[node.right],
//                    F[node.center],F[node.right],
//                                  P[node.center],P[node.right])<<'\t'<<0.0<<std::endl;
  return flf.Dflux_r(DF[node.center],DF[node.right],
                   F[node.center],F[node.right],
                   P[node.center],P[node.right])*Fsx[node.interRight];
  //upwind by hand for v>0
  //return 0.0;
  //return 0.5*DF[node.right]*Fsx[node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilFy(const Vec& F)
{ 
  return 0.5*(F[node.center]+F[node.back])*Fsy[node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDFy_center(const Vec& DF)
{ 
  return 0.5*DF[node.center]*Fsy[node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDFy_back(const Vec& DF)
{ 
  return 0.5*DF[node.back]*Fsy[node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilFz(const Vec& F)
{ 
  return 0.5*(F[node.center]+F[node.top])*Fsz[node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDFz_center(const Vec& DF)
{ 
  return 0.5*DF[node.center]*Fsz[node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::ilDFz_top(const Vec& DF)
{ 
  return 0.5*DF[node.top]*Fsz[node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::DpDx(const Vec& P)
{ 
  return (oneOverdx/((*local_dy_dxi)[node.interRight]))*(P[node.right] - P[node.center]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::DpDy(const Vec& P)
{ 
  return oneOverdy*(P[node.back]  - P[node.center]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::DpDz(const Vec& P)
{ 
  return oneOverdz*(P[node.top]   - P[node.center]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFlux_x(int v,const Vec& K, const Vec& F, const Vec& P)
{
  real fx=ilFx(F,P),
    grad=DpDx(P);
  flux_x[v][node.interRight] = -ilKx(K,grad) * grad + fx;
}


template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFlux_x(const Vec* K, 
                                                           const Vec* F, 
                                                           const Vec* P)
{	
  for (vi=0;vi<nv;vi++)
    {
      real fx=ilFx(F[vi],P[vi]),
        grad=DpDx(P[vi]);
      flux_x[vi][node.interRight] = -ilKx(K[vi],grad) * grad + fx;
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFlux_y(const Vec* K, 
                                                           const Vec* F, 
                                                           const Vec* P)
{
  for (vi=0;vi<nv;vi++)
    {
      real grad=DpDy(P[vi]),
        fy=ilFy(F[vi]);
      flux_y[vi][node.interBack]  = -ilKy(K[vi],grad) * grad + fy;
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFlux_z(const Vec* K, 
                                                           const Vec* F,
                                                           const Vec* P)
{
  for (vi=0;vi<nv;vi++)
    {
      real grad=DpDz(P[vi]),
        fz=ilFz(F[vi]);
      flux_z[vi][node.interTop]	  = -ilKz(K[vi],grad) * grad + fz;
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(const Vec* K, 
                                                            const Vec* F, 
                                                            const Vec* P, 
                                                            const VecVecVec& DK, 
                                                            const VecVecVec& DF)
{
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          real grad=DpDx(P[vi]);
          Dflux_x_center[vi][vj][node.interRight] = -ilDKx_center(DK[vi][vj],grad)*grad + ilDFx_center(DF[vi][vj],F[vi],P[vi]);
          
          Dflux_x_right[vi][vj][node.interRight] = -ilDKx_right(DK[vi][vj],grad)*grad + ilDFx_right(DF[vi][vj],F[vi],P[vi]);
          
        }
    }
  
  for (vi=0;vi<nv;vi++)
    {
      real grad=DpDx(P[vi]);
      Dflux_x_center[vi][vi][node.interRight] += 
        -ilKx(K[vi],grad)*(-oneOverdx/((*local_dy_dxi)[node.interRight]));
      
      Dflux_x_right[vi][vi][node.interRight] +=
        -ilKx(K[vi],grad)*(oneOverdx/((*local_dy_dxi)[node.interRight]));
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(const Vec* K, 
                                                            const Vec* F, 
                                                            const Vec* P, 
                                                            const VecVecVec& DK, 
                                                            const VecVecVec& DF)
{
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          real grad=DpDy(P[vi]);
          Dflux_y_center[vi][vj][node.interBack] = -ilDKy_center(DK[vi][vj],grad)*grad + ilDFy_center(DF[vi][vj]);
          
          Dflux_y_back[vi][vj][node.interBack] = -ilDKy_back(DK[vi][vj],grad) * grad + ilDFy_back(DF[vi][vj]);
          
        }
    }
  
  for (vi=0;vi<nv;vi++)
    {
      real grad=DpDy(P[vi]);
      Dflux_y_center[vi][vi][node.interBack] += 
        -ilKy(K[vi],grad)*(-oneOverdy);
      
      Dflux_y_back[vi][vi][node.interBack] +=
        -ilKy(K[vi],grad)*(oneOverdy);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setDFlux_z(const Vec* K, 
                                                            const Vec* F,
                                                            const Vec* P, 
                                                            const VecVecVec& DK, 
                                                            const VecVecVec& DF)
{
  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          real grad=DpDz(P[vi]);
          Dflux_z_center[vi][vj][node.interTop] = -ilDKz_center(DK[vi][vj],grad)*grad + ilDFz_center(DF[vi][vj]);
          
          Dflux_z_top[vi][vj][node.interTop] = -ilDKz_top(DK[vi][vj],grad) *grad + ilDFz_top(DF[vi][vj]);
          
        }
    }
  
  for (vi=0;vi<nv;vi++)
    { 
      real grad=DpDz(P[vi]);
      Dflux_z_center[vi][vi][node.interTop] += 
        -ilKz(K[vi],grad)*(-oneOverdz);
      
      Dflux_z_top[vi][vi][node.interTop] +=
        -ilKz(K[vi],grad)*(oneOverdz);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getDiv(SecondOrderFd& n, int div_i)
{
  return oneOverdx*(flux_x[div_i][node.interLeft] - 
                    flux_x[div_i][node.interRight])/
    (0.5*((*local_dy_dxi)[node.interLeft] + (*local_dy_dxi)[node.interRight]))
    + oneOverdy*(flux_y[div_i][node.interFront] - 
                 flux_y[div_i][node.interBack])
    + oneOverdz*(flux_z[div_i][node.interBottom] - 
                 flux_z[div_i][node.interTop]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getDivJacCenter(SecondOrderFd& n, int div_i, int u_j)
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getDivJacLeft(SecondOrderFd& n, int div_i, int u_j)
{
  return oneOverdx*
    Dflux_x_center[div_i][u_j][node.interLeft];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getDivJacRight(SecondOrderFd& n, int div_i, int u_j)
{
  return -oneOverdx*
    Dflux_x_right[div_i][u_j][node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getDivJacFront(SecondOrderFd& n, int div_i, int u_j)
{
  return oneOverdy*
    Dflux_y_center[div_i][u_j][node.interFront];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getDivJacBack(SecondOrderFd& n, int div_i, int u_j)
{
  return -oneOverdy*
    Dflux_y_back[div_i][u_j][node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getDivJacBottom(SecondOrderFd& n, int div_i, int u_j)
{
  return oneOverdz*
    Dflux_z_center[div_i][u_j][node.interBottom]; 
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getDivJacTop(SecondOrderFd& n, int div_i, int u_j)
{
  return -oneOverdz*
    Dflux_z_top[div_i][u_j][node.interTop];
}


template<class FLUX_K, class FLUX_F, class BC, int nv>
const Vec& Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getFlux_x(int n){return flux_x[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
const Vec& Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getFlux_y(int n){return flux_y[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
const Vec& Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::getFlux_z(int n){return flux_z[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
Vec& Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFlux_x(int n){return flux_x[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
Vec& Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFlux_y(int n){return flux_y[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
Vec& Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>::setFlux_z(int n){return flux_z[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusFMM1d<FLUX_K,FLUX_F,BC,nv>::Div_KgradPlusFMM1d(BC* bcIn, 
                                          SecondOrderFd& nodeIn,
                                          int nNodes,
                                          real oneOverd):
   Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>(bcIn,nodeIn)
  {
    nxNodes=nNodes;
    local_nxNodes = nodeIn.local_nxNodes;

    oneOverdx=oneOverd;

    Dflux_x_center.resize(nv);
    Dflux_x_right.resize(nv);
    Ksx.newsize(Vec::LOCAL,local_nxNodes+1);
    Ksx = 0.0;
    Fsx.newsize(Vec::LOCAL,local_nxNodes+1);
    Fsx = 0.0;
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusFMM1d<FLUX_K,FLUX_F,BC,nv>::~Div_KgradPlusFMM1d(){}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM1d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec*, const Vec*, 
                                                              const Vec* P)
{  
  for (k= -node.ghost_xOffSet; k< (node.ghost_nxNodes-node.ghost_xOffSet-1); k++)
    {
      node.localIndex(k);
      setFlux_x(K,Fx,P);
    }
  
  for (vi=0; vi<nv; vi++)
    bc[vi].applyNeumannConditions(node,flux_x[vi]);
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM1d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceK(const Vec& K)
{
  for (k= -node.ghost_xOffSet; k< (node.ghost_nxNodes-node.ghost_xOffSet-1); k++)
    {
      node.localIndex(k);
      setKsx(K);
    }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM1d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceF(const Vec& Fx, const Vec&, const Vec&)
{
  for (k= -node.ghost_xOffSet; k< (node.ghost_nxNodes-node.ghost_xOffSet-1); k++)
    {
      node.localIndex(k);
      setFsx(Fx);
    }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM1d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                              const Vec* P, 
                                                              Vec* Div)
{
  computeDivergence(K,Fx,Fy,Fz,P);
  for (vi=0;vi<nv;vi++)
    {
      //      node.localIndex(0);
      for (k=0; k<local_nxNodes; k++)
	{
	  node.localIndex(k);
	  Div[vi][node.center_noGhost] = oneOverdx*
	    (flux_x[vi][node.interLeft] - flux_x[vi][node.interRight])/
            (0.5*( (*local_dy_dxi)[node.interLeft] + (*local_dy_dxi)[node.interRight]));
	}
    }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM1d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec*, const Vec*, 
                                                                 const Vec* P,
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec&, const VecVecVec&)
{
  //   node.localIndex(-node.ghost_xOffSet);
   for (k=-node.ghost_xOffSet; k< node.ghost_nxNodes-node.ghost_xOffSet-1; k++)
     {
       node.localIndex(k);
       setDFlux_x(K,Fx,P,DK,DFx);
     }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM1d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                                 const Vec* P,
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz,
                                                                 VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Fx,Fy,Fz,P,DK,DFx,DFy,DFz);

  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
	{
	  node.localIndex(0);
	  for (k=0;k<local_nxNodes;k++)
	    {
	      DivJac[vi][vj][LEFT][node.center_noGhost] = oneOverdx*
		Dflux_x_center[vi][vj][node.interLeft]/
            (0.5*( (*local_dy_dxi)[node.interLeft] + (*local_dy_dxi)[node.interRight]));
	      
              DivJac[vi][vj][CENTER][node.center_noGhost] = oneOverdx*
                (Dflux_x_right[vi][vj][node.interLeft] 
                 -Dflux_x_center[vi][vj][node.interRight])/
            (0.5*( (*local_dy_dxi)[node.interLeft] + (*local_dy_dxi)[node.interRight])); 
              
              DivJac[vi][vj][RIGHT][node.center_noGhost]  = -oneOverdx*
                Dflux_x_right[vi][vj][node.interRight]/
            (0.5*( (*local_dy_dxi)[node.interLeft] + (*local_dy_dxi)[node.interRight])); 
	      ++node.center_noGhost;++node.interLeft;++node.interRight;
            }

          bc[vi].applyNeumannDerivatives(node,DivJac[vi][vj],
                                         Dflux_x_center[vi][vj],Dflux_x_right[vi][vj],
                                         oneOverdx);
        }
    }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusFMM2d<FLUX_K,FLUX_F,BC,nv>::Div_KgradPlusFMM2d(BC* bcIn,
                                          SecondOrderFd& nodeIn,
                                          int nxNodesIn,int nyNodesIn,
                                          real oneOverdxIn, real oneOverdyIn):
    Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>(bcIn,nodeIn)
  {
    nxNodes=nxNodesIn;
    local_nxNodes = nodeIn.local_nxNodes;
    oneOverdx = oneOverdxIn;

    nyNodes=nyNodesIn;
    local_nyNodes = nodeIn.local_nyNodes;
    oneOverdy = oneOverdyIn;

    Dflux_x_center.resize(nv);
    Dflux_x_right.resize(nv);
    Dflux_y_center.resize(nv);
    Dflux_y_back.resize(nv);
    Ksx.newsize(Vec::LOCAL,local_nyNodes*(local_nxNodes+1));
    Ksy.newsize(Vec::LOCAL,(local_nyNodes+1)*local_nxNodes);
    Ksx = 0.0;
    Ksy = 0.0;
    Fsx.newsize(Vec::LOCAL,local_nyNodes*(local_nxNodes+1));
    Fsy.newsize(Vec::LOCAL,(local_nyNodes+1)*local_nxNodes);
    Fsx = 0.0;
    Fsy = 0.0;
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusFMM2d<FLUX_K,FLUX_F,BC,nv>::~Div_KgradPlusFMM2d(){}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM2d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec*, 
                                                              const Vec* P)
{
  //bottom ghost line
  if (node.ghost_yOffSet)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(-node.ghost_yOffSet,k);
        setFlux_y(K,Fy,P);
      }
  //top ghost line
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(local_nyNodes-1,k);
        setFlux_y(K,Fy,P);
      }
  
  //left ghost line
  if (node.ghost_xOffSet)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,-node.ghost_xOffSet);
        setFlux_x(K,Fx,P);
      }
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,local_nxNodes-1);
        setFlux_x(K,Fx,P);
      }
  
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(j,k);
          setFlux_x(K,Fx,P);
          setFlux_y(K,Fy,P);
        }
      //last y flux in row
      node.localIndex(j,k);
      setFlux_y(K,Fy,P);
    }
  //last row of x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setFlux_x(K,Fx,P);
    }
  
  for (vi=0;vi<nv;vi++)
    bc[vi].applyNeumannConditions(node,flux_x[vi], flux_y[vi]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM2d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceK(const Vec& K)
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM2d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec&)
{
  //bottom ghost line
  if (node.ghost_yOffSet)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(-node.ghost_yOffSet,k);
        setFsy(Fy);
      }
  //top ghost line
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(local_nyNodes-1,k);
        setFsy(Fy);
      }
  
  //left ghost line
  if (node.ghost_xOffSet)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,-node.ghost_xOffSet);
        setFsx(Fx);
      }
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,local_nxNodes-1);
        setFsx(Fx);
      }
  
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(j,k);
          setFsx(Fx);
          setFsy(Fy);
        }
      //last y flux in row
      node.localIndex(j,k);
      setFsy(Fy);
    }
  //last row of x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setFsx(Fx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM2d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                              const Vec* P, 
                                                              Vec* Div)
{
  computeDivergence(K,Fx,Fy,Fz,P);
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM2d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec*, 
                                                                 const Vec* P, 
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec&)
{
  //bottom ghost line
  if (node.ghost_yOffSet)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(-node.ghost_yOffSet,k);
        setDFlux_y(K,Fy,P,DK,DFy);
      }
  //top ghost line
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for(k=0;k<local_nxNodes;k++)
      {
        node.localIndex(local_nyNodes-1,k);
        setDFlux_y(K,Fy,P,DK,DFy);
      }
  
  //left ghost line
  if (node.ghost_xOffSet)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,-node.ghost_xOffSet);
        setDFlux_x(K,Fx,P,DK,DFx);
      }
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (j=0;j<local_nyNodes;j++)
      {
        node.localIndex(j,local_nxNodes-1);
        setDFlux_x(K,Fx,P,DK,DFx);
      }
  
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(j,k);
          setDFlux_x(K,Fx,P,DK,DFx);
          setDFlux_y(K,Fy,P,DK,DFy);
        }
      //last y flux in row
      node.localIndex(j,k);
      setDFlux_y(K,Fy,P,DK,DFy);
    }
  //last row of x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setDFlux_x(K,Fx,P,DK,DFx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM2d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                                 const Vec* P, 
                                                                 const VecVecVec& DK, const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz, 
                                                                 VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Fx,Fy,Fz,P,DK,DFx,DFy,DFz);
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


template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusFMM3d<FLUX_K,FLUX_F,BC,nv>::Div_KgradPlusFMM3d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn,int nzNodesIn, 
             real oneOverdxIn, real oneOverdyIn, real oneOverdzIn):
    Div_KgradPlusFMM<FLUX_K,FLUX_F,BC,nv>(bcIn,nodeIn)
    {
      nxNodes=nxNodesIn;
      local_nxNodes=nodeIn.local_nxNodes;
      oneOverdx = oneOverdxIn;

      nyNodes=nyNodesIn;
      local_nyNodes=nodeIn.local_nyNodes;
      oneOverdy = oneOverdyIn;

      nzNodes=nzNodesIn;
      local_nzNodes=nodeIn.local_nzNodes;
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
      Fsx.newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
      Fsy.newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
      Fsz.newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
      Fsx = 0.0;
      Fsy = 0.0;
      Fsz = 0.0;
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusFMM3d<FLUX_K,FLUX_F,BC,nv>::~Div_KgradPlusFMM3d(){}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM3d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                              const Vec* P)
{
  //ghost faces
  if (node.ghost_zOffSet)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(-node.ghost_zOffSet,j,k);
          setFlux_z(K,Fz,P);
        }
  if (node.ghost_nzNodes -  node.ghost_zOffSet - node.local_nzNodes)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(local_nzNodes-1,j,k);
          setFlux_z(K,Fz,P);
        }
  if (node.ghost_yOffSet)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,-node.ghost_yOffSet,k);
          setFlux_y(K,Fy,P);
        }
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,local_nyNodes-1,k);
          setFlux_y(K,Fy,P);
        }
  if (node.ghost_xOffSet)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,-node.ghost_xOffSet);
          setFlux_x(K,Fx,P);
        }      
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,local_nxNodes-1);
          setFlux_x(K,Fx,P);
        }
  
  //interior
  for (i=0;i<local_nzNodes-1;i++)
    {
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(i,j,k);
              setFlux_x(K,Fx,P);
              setFlux_y(K,Fy,P);
              setFlux_z(K,Fz,P);
            }
          // last y and z fluxes of x row
          node.localIndex(i,j,k);
          setFlux_y(K,Fy,P);
          setFlux_z(K,Fz,P);
        }
      //last row of x and z fluxes in y plane
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setFlux_x(K,Fx,P);
          setFlux_z(K,Fz,P);
        }
      // last z fluz in x row
      node.localIndex(i,j,k);
      setFlux_z(K,Fz,P);
    }
  //top face of x and y fluxes
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setFlux_x(K,Fx,P);
          setFlux_y(K,Fy,P);
        }
      //last y flux in x row
      node.localIndex(i,j,k);
      setFlux_y(K,Fy,P);
    }
  //last row of x fluxes in top face
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(i,j,k);
      setFlux_x(K,Fx,P);
    }
  
  for (vi=0;vi<nv;vi++)
    bc[vi].applyNeumannConditions(node,flux_x[vi], flux_y[vi], flux_z[vi]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM3d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceK(const Vec& K)
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM3d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec& Fz)
{
  //ghost faces
  if (node.ghost_zOffSet)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(-node.ghost_zOffSet,j,k);
          setFsz(Fz);
        }
  if (node.ghost_nzNodes -  node.ghost_zOffSet - node.local_nzNodes)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(local_nzNodes-1,j,k);
          setFsz(Fz);
        }
  if (node.ghost_yOffSet)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,-node.ghost_yOffSet,k);
          setFsy(Fy);
        }
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,local_nyNodes-1,k);
          setFsy(Fy);
        }
  if (node.ghost_xOffSet)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,-node.ghost_xOffSet);
          setFsx(Fx);
        }      
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,local_nxNodes-1);
          setFsx(Fx);
        }
  
  //interior
  for (i=0;i<local_nzNodes-1;i++)
    {
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(i,j,k);
              setFsx(Fx);
              setFsy(Fy);
              setFsz(Fz);
            }
          // last y and z fluxes of x row
          node.localIndex(i,j,k);
          setFsy(Fy);
          setFsz(Fz);
        }
      //last row of x and z fluxes in y plane
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setFsx(Fx);
          setFsz(Fz);
        }
      // last z fluz in x row
      node.localIndex(i,j,k);
      setFsz(Fz);
    }
  //top face of x and y fluxes
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setFsx(Fx);
          setFsy(Fy);
        }
      //last y flux in x row
      node.localIndex(i,j,k);
      setFsy(Fy);
    }
  //last row of x fluxes in top face
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(i,j,k);
      setFsx(Fx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM3d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                              const Vec* P, 
                                                              Vec* Div)
{
  computeDivergence(K,Fx,Fy,Fz,P);
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
  
template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM3d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                                 const Vec* P, 
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz)
{
  //ghost faces
  if (node.ghost_zOffSet)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(-node.ghost_zOffSet,j,k);
          setDFlux_z(K,Fz,P,DK,DFz);
        }
  if (node.ghost_nzNodes -  node.ghost_zOffSet - node.local_nzNodes)
    for (j=0;j<local_nyNodes;j++)
      for (k=0;k<local_nxNodes;k++)
        {
          node.localIndex(local_nzNodes-1,j,k);
          setDFlux_z(K,Fz,P,DK,DFz);
        }
  if (node.ghost_yOffSet)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,-node.ghost_yOffSet,k);
          setDFlux_y(K,Fy,P,DK,DFy);
        }
  if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
    for (i=0;i<local_nzNodes;i++)
      for(k=0;k<local_nxNodes;k++)
        {
          node.localIndex(i,local_nyNodes-1,k);
          setDFlux_y(K,Fy,P,DK,DFy);
        }
  if (node.ghost_xOffSet)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,-node.ghost_xOffSet);
          setDFlux_x(K,Fx,P,DK,DFx);
        }      
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    for (i=0;i<local_nzNodes;i++)
      for (j=0;j<local_nyNodes;j++)
        {
          node.localIndex(i,j,local_nxNodes-1);
          setDFlux_x(K,Fx,P,DK,DFx);
        }
  
  for (i=0;i<local_nzNodes-1;i++)
    {
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              //interior
              node.localIndex(i,j,k);
              setDFlux_x(K,Fx,P,DK,DFx);
              setDFlux_y(K,Fy,P,DK,DFy);
              setDFlux_z(K,Fz,P,DK,DFz);
            }
          // last y and z fluxes of x row
          node.localIndex(i,j,k);
          setDFlux_y(K,Fy,P,DK,DFy);
          setDFlux_z(K,Fz,P,DK,DFz);
        }
      //last row of x and z fluxes in y plane
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setDFlux_x(K,Fx,P,DK,DFx);
          setDFlux_z(K,Fz,P,DK,DFz);
        }
      // last z fluz in x row
      node.localIndex(i,j,k);
      setDFlux_z(K,Fz,P,DK,DFz);
    }
  //top face of x and y fluxes
  for (j=0;j<local_nyNodes-1;j++)
    {
      for (k=0;k<local_nxNodes-1;k++)
        {
          node.localIndex(i,j,k);
          setDFlux_x(K,Fx,P,DK,DFx);
          setDFlux_y(K,Fy,P,DK,DFy);
        }
      //last y flux in x row
      node.localIndex(i,j,k);
      setDFlux_y(K,Fy,P,DK,DFy);
    }
  //last row of x fluxes in top face
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(i,j,k);
      setDFlux_x(K,Fx,P,DK,DFx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusFMM3d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                                 const Vec* P, 
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz,
                                                                 VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Fx,Fy,Fz,P,DK,DFx,DFy,DFz);
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
