#ifndef DIV_KGRADPLUSF_H
#define DIV_KGRADPLUSF_H
#include "Definitions.h"
#include "PetscSecondOrderFd.h"
#include "Divergence.h"

#include <vector>

//  FIX DIVERGENCE JAC -- MUST APPLY CORRECT DERIVATIVES FOR NEUMAN BC'S

namespace Daetk 
{

using Petsc::SecondOrderFd;

  class SimpleUpwindFlux
  {
  public:
    SimpleUpwindFlux(){}
    ~SimpleUpwindFlux(){}
    inline real flux(const real& f_l, const real& f_r,
                              const real& u_l, const real& u_r)
    {
      if ( u_r == u_l)
        return f_l;
      else if ( (f_r - f_l) / (u_r - u_l) > 0 )
        return f_l;
      else
        return f_r;
    }
    
    inline real Dflux_l(const real& df_l, const real& df_r,
                               const real& f_l, const real& f_r,
                               const real& u_l, const real& u_r)
    {
      if ( u_r == u_l)
        return df_l;
      else if ( (f_r - f_l) / (u_r - u_l) > 0 )
        return df_l;
      else
        return 0.0;
    }
    
    inline real Dflux_r(const real& df_l, const real& df_r,
                               const real& f_l, const real& f_r,
                               const real& u_l, const real& u_r)
    {
      if ( u_r == u_l)
        return 0.0;
      else if ( (f_r - f_l) / (u_r - u_l) > 0 )
        return 0.0;
      else
        return df_r;
    }
  };

  class MidpointFlux
  {
  public:
    MidpointFlux(){}
    ~MidpointFlux(){}
    inline real flux(const real& f_l, const real& f_r,
                              const real& u_l, const real& u_r)
    {
      return 0.5*(f_r + f_l);
    }
    inline real Dflux_l(const real& df_l, const real& df_r,
                              const real& f_l, const real& f_r,
                              const real& u_l, const real& u_r)
    {
      return 0.5*df_l;
    }
    inline real Dflux_r(const real& df_l, const real& df_r,
                              const real& f_l, const real& f_r,
                              const real& u_l, const real& u_r)
    {
      return 0.5*df_r;
    }
  };

  
template<class FLUX_K, class FLUX_F, class BC, int nv>
class Div_KgradPlusF : public Divergence
{
public:
  enum JacEntry {CENTER, LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP};
  Div_KgradPlusF(BC* bcIn, SecondOrderFd& nodeIn);
  
  virtual ~Div_KgradPlusF();

  virtual void computeDivergence(const Vec* K, const Vec* Fx, const Vec* Fy, const Vec* Fz, const Vec* P,  
                                 Vec* Div)= 0;

  virtual void computeDivergenceJac(const Vec* K, 
                                    const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                    const Vec* P, 
                                    const VecVecVec& DK, 
                                    const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz, 
                                    VecVecVecVec& DivJac)= 0;  
  virtual void computeDivergence(const Vec* K, 
                                 const Vec* Fx,const Vec* Fy,const Vec* Fz, 
                                 const Vec* P)= 0;

  virtual void computeDivergenceJac(const Vec* K, 
                                    const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                    const Vec* P,
                                    const VecVecVec& DK, 
                                    const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz)= 0;  
  
  virtual void computeInterfaceK(const Vec& K)= 0;
  inline void setKsx(const Vec& K);
  inline void setKsy(const Vec& K);
  inline void setKsz(const Vec& K);

  virtual void computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec& Fz)= 0;
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


  inline real getDiv(SecondOrderFd& n, int div_i= 0);

  inline real getDivJacCenter(SecondOrderFd& n, int div_i= 0, int u_j= 0);
  inline real getDivJacLeft(SecondOrderFd& n, int div_i= 0, int u_j= 0);
  inline real getDivJacRight(SecondOrderFd& n, int div_i= 0, int u_j= 0);
  inline real getDivJacFront(SecondOrderFd& n, int div_i= 0, int u_j= 0);
  inline real getDivJacBack(SecondOrderFd& n, int div_i= 0, int u_j= 0);
  inline real getDivJacBottom(SecondOrderFd& n, int div_i= 0, int u_j= 0);
  inline real getDivJacTop(SecondOrderFd& n, int div_i= 0, int u_j= 0);

  const Vec& getFlux_x(int n= 0);
  const Vec& getFlux_y(int n= 0);
  const Vec& getFlux_z(int n= 0);
  
  Vec& setFlux_x(int n= 0);
  Vec& setFlux_y(int n= 0);
  Vec& setFlux_z(int n= 0);
  
  int nxNodes,nyNodes,nzNodes,local_nxNodes,local_nyNodes,local_nzNodes;
  real oneOverdx,oneOverdy,oneOverdz; 
  Vec flux_x[nv],flux_y[nv],flux_z[nv],Ksx,Ksy,Ksz,Fsx,Fsy,Fsz;
  VecVecVec Dflux_x_center,Dflux_x_right,Dflux_y_center,Dflux_y_back,Dflux_z_center,Dflux_z_top;
  BC* bc;
  SecondOrderFd& node;
  FLUX_K flk;
  FLUX_F flf;
};

template<class FLUX_K, class FLUX_F, class BC, int nv>
class Div_KgradPlusF1d : public Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>
{
public:
  Div_KgradPlusF1d(BC* bcIn, 
             SecondOrderFd& nodeIn,
             int nNodes, 
             real oneOverd);

  virtual ~Div_KgradPlusF1d();
  
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
class Div_KgradPlusF2d : public Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>
{
public:
  Div_KgradPlusF2d(BC* bcIn,
                   SecondOrderFd& nodeIn,
                   int nxNodesIn,int nyNodesIn,
                   real oneOverdxIn, real oneOverdyIn);
  virtual ~Div_KgradPlusF2d();
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
class Div_KgradPlusF3d : public Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>
{
public:
  Div_KgradPlusF3d(BC* bcIn,
                   SecondOrderFd& nodeIn,
                   int nxNodesIn,int nyNodesIn,int nzNodesIn,
                   real oneOverdxIn, real oneOverdyIn, real oneOverdzIn);
  virtual ~Div_KgradPlusF3d();
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
Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::Div_KgradPlusF(BC* bcIn, SecondOrderFd& nodeIn):
  nxNodes(1),
  nyNodes(1),
  nzNodes(1),
  local_nxNodes(1),
  local_nyNodes(1),
  local_nzNodes(1),
  bc(bcIn),
  node(nodeIn)
{
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
  this->Fsx.newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
  this->Fsy.newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
  this->Fsz.newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
  this->Fsx = 0.0;
  this->Fsy = 0.0;
  this->Fsz = 0.0;
  for (int vi= 0;vi<nv;vi++)
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
      for (int vj= 0;vj<nv;vj++)
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::~Div_KgradPlusF(){}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(const Vec& K)
{
  this->Ksx[this->node.interRight] = 2.0 / ( 1.0/K[this->node.center] + 1.0/K[this->node.right]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(const Vec& K)
{
  this->Ksy[this->node.interBack] = 2.0 / ( 1.0/K[this->node.center] + 1.0/K[this->node.back]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsz(const Vec& K)
{
  this->Ksz[this->node.interTop] = 2.0 / ( 1.0/K[this->node.center] + 1.0/K[this->node.top]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(const Vec& F)
{
  this->Fsx[this->node.interRight] = 2.0 / ( 1.0/F[this->node.center] + 1.0/F[this->node.right]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(const Vec& F)
{
  this->Fsy[this->node.interBack] = 2.0 / ( 1.0/F[this->node.center] + 1.0/F[this->node.back]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsz(const Vec& F)
{
  this->Fsz[this->node.interTop] = 2.0 / ( 1.0/F[this->node.center] + 1.0/F[this->node.top]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKx(const Vec& K, const real& grad )
{ 
  return 0.5*(K[this->node.center]+K[this->node.right])*this->Ksx[this->node.interRight];
}


template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKx_center(const Vec& DK, const real& grad)
{
  return 0.5*DK[this->node.center]*this->Ksx[this->node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKx_right(const Vec& DK, const real& grad)
{
  return 0.5*DK[this->node.right]*this->Ksx[this->node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKy(const Vec& K, const real& grad)
{ 
  return 0.5*(K[this->node.center]+K[this->node.back])*this->Ksy[this->node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKy_center(const Vec& DK,const real& grad)
{ 
  return 0.5*DK[this->node.center]*this->Ksy[this->node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKy_back(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[this->node.back]*this->Ksy[this->node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKz(const Vec& K, const real& grad)
{ 
  return 0.5*(K[this->node.center]+K[this->node.top])*this->Ksz[this->node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKz_center(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[this->node.center]*this->Ksz[this->node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKz_top(const Vec& DK, const real& grad)
{ 
  return 0.5*DK[this->node.top]*this->Ksz[this->node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilFx(const Vec& F, const Vec& P)
{ 
//   std::cout<<"this->flf.flux "<<this->flf.flux(F[this->node.center],F[this->node.right],
//                                 P[this->node.center],P[this->node.right])<<'\t'<<F[this->node.center]<<std::endl;
  return this->flf.flux(F[this->node.center],F[this->node.right],
                  P[this->node.center],P[this->node.right])*this->Fsx[this->node.interRight];
  //upwind by hand for v>0
  //return F[this->node.center]*this->Fsx[this->node.interRight];;
  //return 0.5*(F[this->node.center]+F[this->node.right])*this->Fsx[this->node.interRight];;
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFx_center(const Vec& DF, const Vec& F, const Vec& P)
{ 
//   std::cout<<"this->flf.Dflux_l "<<this->flf.Dflux_l(DF[this->node.center],DF[this->node.right],
//                    F[this->node.center],F[this->node.right],
//                                  P[this->node.center],P[this->node.right])<<'\t'<<DF[this->node.center]<<std::endl;
  return this->flf.Dflux_l(DF[this->node.center],DF[this->node.right],
                   F[this->node.center],F[this->node.right],
                   P[this->node.center],P[this->node.right])*this->Fsx[this->node.interRight];
  //upwind by hand for v>0
  //return DF[this->node.center]*this->Fsx[this->node.interRight];
  //return 0.5*DF[this->node.center]*this->Fsx[this->node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFx_right(const Vec& DF, const Vec& F, const Vec& P)
{ 
//   std::cout<<"this->flf.Dflux_r "<<this->flf.Dflux_r(DF[this->node.center],DF[this->node.right],
//                    F[this->node.center],F[this->node.right],
//                                  P[this->node.center],P[this->node.right])<<'\t'<<0.0<<std::endl;
  return this->flf.Dflux_r(DF[this->node.center],DF[this->node.right],
                   F[this->node.center],F[this->node.right],
                   P[this->node.center],P[this->node.right])*this->Fsx[this->node.interRight];
  //upwind by hand for v>0
  //return 0.0;
  //return 0.5*DF[this->node.right]*this->Fsx[this->node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilFy(const Vec& F)
{ 
  return 0.5*(F[this->node.center]+F[this->node.back])*this->Fsy[this->node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFy_center(const Vec& DF)
{ 
  return 0.5*DF[this->node.center]*this->Fsy[this->node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFy_back(const Vec& DF)
{ 
  return 0.5*DF[this->node.back]*this->Fsy[this->node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilFz(const Vec& F)
{ 
  return 0.5*(F[this->node.center]+F[this->node.top])*this->Fsz[this->node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFz_center(const Vec& DF)
{ 
  return 0.5*DF[this->node.center]*this->Fsz[this->node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFz_top(const Vec& DF)
{ 
  return 0.5*DF[this->node.top]*this->Fsz[this->node.interTop];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDx(const Vec& P)
{ 
  return this->oneOverdx*(P[this->node.right] - P[this->node.center]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDy(const Vec& P)
{ 
  return this->oneOverdy*(P[this->node.back]  - P[this->node.center]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDz(const Vec& P)
{ 
  return this->oneOverdz*(P[this->node.top]   - P[this->node.center]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(int v,const Vec& K, const Vec& F, const Vec& P)
{
  real fx=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilFx(F,P),
    grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDx(P);
  this->flux_x[v][this->node.interRight] = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKx(K,grad) * grad + fx;
}


template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(const Vec* K, 
                                                           const Vec* F, 
                                                           const Vec* P)
{	
  for (int vi= 0;vi<nv;vi++)
    {
      real fx=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilFx(F[vi],P[vi]),
        grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDx(P[vi]);
      this->flux_x[vi][this->node.interRight] = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKx(K[vi],grad) * grad + fx;
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(const Vec* K, 
                                                           const Vec* F, 
                                                           const Vec* P)
{
  for (int vi= 0;vi<nv;vi++)
    {
      real grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDy(P[vi]),
        fy=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilFy(F[vi]);
      this->flux_y[vi][this->node.interBack]  = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKy(K[vi],grad) * grad + fy;
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_z(const Vec* K, 
                                                           const Vec* F,
                                                           const Vec* P)
{
  for (int vi= 0;vi<nv;vi++)
    {
      real grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDz(P[vi]),
        fz=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilFz(F[vi]);
      this->flux_z[vi][this->node.interTop]	  = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKz(K[vi],grad) * grad + fz;
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(const Vec* K, 
                                                            const Vec* F, 
                                                            const Vec* P, 
                                                            const VecVecVec& DK, 
                                                            const VecVecVec& DF)
{
  for (int vi= 0;vi<nv;vi++)
    {
      for (int vj= 0;vj<nv;vj++)
        {
          real grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDx(P[vi]);
          this->Dflux_x_center[vi][vj][this->node.interRight] = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKx_center(DK[vi][vj],grad)*grad + Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFx_center(DF[vi][vj],F[vi],P[vi]);
          
          this->Dflux_x_right[vi][vj][this->node.interRight] = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKx_right(DK[vi][vj],grad)*grad + Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFx_right(DF[vi][vj],F[vi],P[vi]);
          
        }
    }
  
  for (int vi= 0;vi<nv;vi++)
    {
      real grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDx(P[vi]);
      this->Dflux_x_center[vi][vi][this->node.interRight] += 
        -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKx(K[vi],grad)*(-this->oneOverdx);
      
      this->Dflux_x_right[vi][vi][this->node.interRight] +=
        -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKx(K[vi],grad)*(this->oneOverdx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(const Vec* K, 
                                                            const Vec* F, 
                                                            const Vec* P, 
                                                            const VecVecVec& DK, 
                                                            const VecVecVec& DF)
{
  for (int vi= 0;vi<nv;vi++)
    {
      for (int vj= 0;vj<nv;vj++)
        {
          real grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDy(P[vi]);
          this->Dflux_y_center[vi][vj][this->node.interBack] = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKy_center(DK[vi][vj],grad)*grad + Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFy_center(DF[vi][vj]);
          
          this->Dflux_y_back[vi][vj][this->node.interBack] = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKy_back(DK[vi][vj],grad) * grad + Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFy_back(DF[vi][vj]);
          
        }
    }
  
  for (int vi= 0;vi<nv;vi++)
    {
      real grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDy(P[vi]);
      this->Dflux_y_center[vi][vi][this->node.interBack] += 
        -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKy(K[vi],grad)*(-this->oneOverdy);
      
      this->Dflux_y_back[vi][vi][this->node.interBack] +=
        -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKy(K[vi],grad)*(this->oneOverdy);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline void Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_z(const Vec* K, 
                                                            const Vec* F,
                                                            const Vec* P, 
                                                            const VecVecVec& DK, 
                                                            const VecVecVec& DF)
{
  for (int vi= 0;vi<nv;vi++)
    {
      for (int vj= 0;vj<nv;vj++)
        {
          real grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDz(P[vi]);
          this->Dflux_z_center[vi][vj][this->node.interTop] = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKz_center(DK[vi][vj],grad)*grad + Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFz_center(DF[vi][vj]);
          
          this->Dflux_z_top[vi][vj][this->node.interTop] = -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDKz_top(DK[vi][vj],grad) *grad + Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilDFz_top(DF[vi][vj]);
          
        }
    }
  
  for (int vi= 0;vi<nv;vi++)
    { 
      real grad=Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::DpDz(P[vi]);
      this->Dflux_z_center[vi][vi][this->node.interTop] += 
        -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKz(K[vi],grad)*(-this->oneOverdz);
      
      this->Dflux_z_top[vi][vi][this->node.interTop] +=
        -Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::ilKz(K[vi],grad)*(this->oneOverdz);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getDiv(SecondOrderFd& n, int div_i)
{
  return this->oneOverdx*(this->flux_x[div_i][this->node.interLeft] - 
                    this->flux_x[div_i][this->node.interRight])
    + this->oneOverdy*(this->flux_y[div_i][this->node.interFront] - 
                 this->flux_y[div_i][this->node.interBack])
    + this->oneOverdz*(this->flux_z[div_i][this->node.interBottom] - 
                 this->flux_z[div_i][this->node.interTop]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getDivJacCenter(SecondOrderFd& n, int div_i, int u_j)
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getDivJacLeft(SecondOrderFd& n, int div_i, int u_j)
{
  return this->oneOverdx*
    this->Dflux_x_center[div_i][u_j][this->node.interLeft];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getDivJacRight(SecondOrderFd& n, int div_i, int u_j)
{
  return -this->oneOverdx*
    this->Dflux_x_right[div_i][u_j][this->node.interRight];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getDivJacFront(SecondOrderFd& n, int div_i, int u_j)
{
  return this->oneOverdy*
    this->Dflux_y_center[div_i][u_j][this->node.interFront];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getDivJacBack(SecondOrderFd& n, int div_i, int u_j)
{
  return -this->oneOverdy*
    this->Dflux_y_back[div_i][u_j][this->node.interBack];
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getDivJacBottom(SecondOrderFd& n, int div_i, int u_j)
{
  return this->oneOverdz*
    this->Dflux_z_center[div_i][u_j][this->node.interBottom]; 
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
inline real Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getDivJacTop(SecondOrderFd& n, int div_i, int u_j)
{
  return -this->oneOverdz*
    this->Dflux_z_top[div_i][u_j][this->node.interTop];
}


template<class FLUX_K, class FLUX_F, class BC, int nv>
const Vec& Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getFlux_x(int n){return this->flux_x[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
const Vec& Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getFlux_y(int n){return this->flux_y[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
const Vec& Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::getFlux_z(int n){return this->flux_z[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
Vec& Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(int n){return this->flux_x[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
Vec& Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(int n){return this->flux_y[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
Vec& Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_z(int n){return this->flux_z[n];}

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusF1d<FLUX_K,FLUX_F,BC,nv>::Div_KgradPlusF1d(BC* bcIn, 
                                          SecondOrderFd& nodeIn,
                                          int nNodes,
                                          real oneOverdxIn):
   Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>(bcIn,nodeIn)
  {
    this->nxNodes=nNodes;
    this->local_nxNodes = nodeIn.local_nxNodes;

    this->oneOverdx=oneOverdxIn;

    this->Dflux_x_center.resize(nv);
    this->Dflux_x_right.resize(nv);
    this->Ksx.newsize(Vec::LOCAL,this->local_nxNodes+1);
    this->Ksx = 0.0;
    this->Fsx.newsize(Vec::LOCAL,this->local_nxNodes+1);
    this->Fsx = 0.0;
    for (int vi= 0;vi<nv;vi++)
      {
        this->flux_x[vi].newsize(Vec::LOCAL,this->local_nxNodes+1);
        this->flux_x[vi] = 0.0;
        this->Dflux_x_center[vi].resize(nv);
        this->Dflux_x_right[vi].resize(nv);
        for (int vj= 0;vj<nv;vj++)
          {
            this->Dflux_x_center[vi][vj].newsize(Vec::LOCAL,this->local_nxNodes+1);
            this->Dflux_x_center[vi][vj] = 0.0;
            this->Dflux_x_right[vi][vj].newsize(Vec::LOCAL,this->local_nxNodes+1);
            this->Dflux_x_right[vi][vj] = 0.0;
          }
      }
  }

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusF1d<FLUX_K,FLUX_F,BC,nv>::~Div_KgradPlusF1d(){}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF1d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec*, const Vec*, 
                                                              const Vec* P)
{  
  for (int k = -this->node.ghost_xOffSet; k< (this->node.ghost_nxNodes-this->node.ghost_xOffSet-1); k++)
    {
      this->node.localIndex(k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
    }
  
  for (int vi= 0; vi<nv; vi++)
    this->bc[vi].applyNeumannConditions(this->node,this->flux_x[vi]);
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF1d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceK(const Vec& K)
{
  for (int k = -this->node.ghost_xOffSet; k< (this->node.ghost_nxNodes-this->node.ghost_xOffSet-1); k++)
    {
      this->node.localIndex(k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
    }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF1d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceF(const Vec& Fx, const Vec&, const Vec&)
{
  for (int k = -this->node.ghost_xOffSet; k< (this->node.ghost_nxNodes-this->node.ghost_xOffSet-1); k++)
    {
      this->node.localIndex(k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
    }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF1d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                              const Vec* P, 
                                                              Vec* Div)
{
  computeDivergence(K,Fx,Fy,Fz,P);
  for (int vi= 0;vi<nv;vi++)
    {
      //      this->node.localIndex(0);
      for (int k = 0; k<this->local_nxNodes; k++)
	{
	  this->node.localIndex(k);
	  Div[vi][this->node.center_noGhost] = this->oneOverdx*
	    (this->flux_x[vi][this->node.interLeft] - this->flux_x[vi][this->node.interRight]);
	}
    }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF1d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec*, const Vec*, 
                                                                 const Vec* P,
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec&, const VecVecVec&)
{
  //   this->node.localIndex(-this->node.ghost_xOffSet);
   for (int k =-this->node.ghost_xOffSet; k< this->node.ghost_nxNodes-this->node.ghost_xOffSet-1; k++)
     {
       this->node.localIndex(k);
       Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
     }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF1d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                                 const Vec* P,
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz,
                                                                 VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Fx,Fy,Fz,P,DK,DFx,DFy,DFz);

  for (int vi= 0;vi<nv;vi++)
    {
      for (int vj= 0;vj<nv;vj++)
	{
	  this->node.localIndex(0);
	  for (int k = 0;k<this->local_nxNodes;k++)
	    {
	      DivJac[vi][vj][this->LEFT][this->node.center_noGhost] = this->oneOverdx*
		this->Dflux_x_center[vi][vj][this->node.interLeft];
	      
              DivJac[vi][vj][this->CENTER][this->node.center_noGhost] = this->oneOverdx*
                (this->Dflux_x_right[vi][vj][this->node.interLeft] 
                 -this->Dflux_x_center[vi][vj][this->node.interRight]); 
              
              DivJac[vi][vj][this->RIGHT][this->node.center_noGhost]  = -this->oneOverdx*
                this->Dflux_x_right[vi][vj][this->node.interRight]; 
	      ++this->node.center_noGhost;++this->node.interLeft;++this->node.interRight;
            }

          this->bc[vi].applyNeumannDerivatives(this->node,DivJac[vi][vj],
                                         this->Dflux_x_center[vi][vj],this->Dflux_x_right[vi][vj],
                                         this->oneOverdx);
        }
    }
}  

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusF2d<FLUX_K,FLUX_F,BC,nv>::Div_KgradPlusF2d(BC* bcIn,
                                          SecondOrderFd& nodeIn,
                                          int nxNodesIn,int nyNodesIn,
                                          real oneOverdxIn, real oneOverdyIn):
    Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>(bcIn,nodeIn)
  {
    this->nxNodes=nxNodesIn;
    this->local_nxNodes = nodeIn.local_nxNodes;
    this->oneOverdx = oneOverdxIn;

    this->nyNodes=nyNodesIn;
    this->local_nyNodes = nodeIn.local_nyNodes;
    this->oneOverdy = oneOverdyIn;

    this->Dflux_x_center.resize(nv);
    this->Dflux_x_right.resize(nv);
    this->Dflux_y_center.resize(nv);
    this->Dflux_y_back.resize(nv);
    this->Ksx.newsize(Vec::LOCAL,this->local_nyNodes*(this->local_nxNodes+1));
    this->Ksy.newsize(Vec::LOCAL,(this->local_nyNodes+1)*this->local_nxNodes);
    this->Ksx = 0.0;
    this->Ksy = 0.0;
    this->Fsx.newsize(Vec::LOCAL,this->local_nyNodes*(this->local_nxNodes+1));
    this->Fsy.newsize(Vec::LOCAL,(this->local_nyNodes+1)*this->local_nxNodes);
    this->Fsx = 0.0;
    this->Fsy = 0.0;
    for (int vi= 0;vi<nv;vi++)
      {
          this->flux_x[vi].newsize(Vec::LOCAL,this->local_nyNodes*(this->local_nxNodes+1));
          this->flux_y[vi].newsize(Vec::LOCAL,(this->local_nyNodes+1)*this->local_nxNodes);
          this->flux_x[vi] = 0.0;
          this->flux_y[vi] = 0.0;
          this->Dflux_x_center[vi].resize(nv);
          this->Dflux_x_right[vi].resize(nv);
          this->Dflux_y_center[vi].resize(nv);
          this->Dflux_y_back[vi].resize(nv);
          for (int vj= 0;vj<nv;vj++)
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusF2d<FLUX_K,FLUX_F,BC,nv>::~Div_KgradPlusF2d(){}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF2d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec*, 
                                                              const Vec* P)
{
  //bottom ghost line
  if (this->node.ghost_yOffSet)
    for(int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(-this->node.ghost_yOffSet,k);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
      }
  //top ghost line
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for(int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(this->local_nyNodes-1,k);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
      }
  
  //left ghost line
  if (this->node.ghost_xOffSet)
    for (int j= 0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,-this->node.ghost_xOffSet);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
      }
  
  //right ghost line
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int j= 0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,this->local_nxNodes-1);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
      }
  
  for (int j= 0;j<this->local_nyNodes-1;j++)
    {
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
        }
      //last y flux in row
      this->node.localIndex(j,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
    }
  //last row of x fluxes
  for (int k = 0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nyNodes-1,k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
    }
  
  for (int vi= 0;vi<nv;vi++)
    this->bc[vi].applyNeumannConditions(this->node,this->flux_x[vi], this->flux_y[vi]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF2d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceK(const Vec& K)
{
  //bottom ghost line
  if (this->node.ghost_yOffSet)
    for(int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(-this->node.ghost_yOffSet,k);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
      }
  //top ghost line
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for(int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(this->local_nyNodes-1,k);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
      }
  
  //left ghost line
  if (this->node.ghost_xOffSet)
    for (int j= 0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,-this->node.ghost_xOffSet);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
      }
  
  //right ghost line
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int j= 0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,this->local_nxNodes-1);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
      }
  
  for (int j= 0;j<this->local_nyNodes-1;j++)
    {
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
        }
      //last y flux in row
      this->node.localIndex(j,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
    }
  //last row of x fluxes
  for (int k = 0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nyNodes-1,k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF2d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec&)
{
  //bottom ghost line
  if (this->node.ghost_yOffSet)
    for(int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(-this->node.ghost_yOffSet,k);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
      }
  //top ghost line
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for(int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(this->local_nyNodes-1,k);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
      }
  
  //left ghost line
  if (this->node.ghost_xOffSet)
    for (int j= 0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,-this->node.ghost_xOffSet);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
      }
  
  //right ghost line
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int j= 0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,this->local_nxNodes-1);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
      }
  
  for (int j= 0;j<this->local_nyNodes-1;j++)
    {
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
        }
      //last y flux in row
      this->node.localIndex(j,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
    }
  //last row of x fluxes
  for (int k = 0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nyNodes-1,k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF2d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                              const Vec* P, 
                                                              Vec* Div)
{
  computeDivergence(K,Fx,Fy,Fz,P);
  for (int j= 0;j<this->local_nyNodes;j++)
    for (int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(j,k);
        for (int vi= 0;vi<nv;vi++)
          {
            Div[vi][this->node.center_noGhost] = this->oneOverdx*
              (this->flux_x[vi][this->node.interLeft] - this->flux_x[vi][this->node.interRight])
              + this->oneOverdy*
              (this->flux_y[vi][this->node.interFront] - this->flux_y[vi][this->node.interBack]);
          }
      }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF2d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec*, 
                                                                 const Vec* P, 
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec&)
{
  //bottom ghost line
  if (this->node.ghost_yOffSet)
    for(int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(-this->node.ghost_yOffSet,k);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
      }
  //top ghost line
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for(int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(this->local_nyNodes-1,k);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
      }
  
  //left ghost line
  if (this->node.ghost_xOffSet)
    for (int j= 0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,-this->node.ghost_xOffSet);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
      }
  
  //right ghost line
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int j= 0;j<this->local_nyNodes;j++)
      {
        this->node.localIndex(j,this->local_nxNodes-1);
        Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
      }
  
  for (int j= 0;j<this->local_nyNodes-1;j++)
    {
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
        }
      //last y flux in row
      this->node.localIndex(j,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
    }
  //last row of x fluxes
  for (int k = 0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nyNodes-1,k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF2d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                                 const Vec* P, 
                                                                 const VecVecVec& DK, const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz, 
                                                                 VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Fx,Fy,Fz,P,DK,DFx,DFy,DFz);
  for (int j= 0;j<this->local_nyNodes;j++)
    for (int k = 0;k<this->local_nxNodes;k++)
      {
        this->node.localIndex(j,k);
        for (int vi= 0;vi<nv;vi++)
          {
            for (int vj= 0;vj<nv;vj++)
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
  
  for (int vi= 0;vi<nv;vi++)
    {
      for (int vj= 0;vj<nv;vj++)
        {
          this->bc[vi].applyNeumannDerivatives(this->node,DivJac[vi][vj],
                                         this->Dflux_x_center[vi][vj],this->Dflux_x_right[vi][vj],
                                         this->Dflux_y_center[vi][vj],this->Dflux_y_back[vi][vj],
                                         this->oneOverdx,this->oneOverdy);
          
        }
    }
}


template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusF3d<FLUX_K,FLUX_F,BC,nv>::Div_KgradPlusF3d(BC* bcIn,
             SecondOrderFd& nodeIn,
             int nxNodesIn,int nyNodesIn,int nzNodesIn, 
             real oneOverdxIn, real oneOverdyIn, real oneOverdzIn):
    Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>(bcIn,nodeIn)
    {
      this->nxNodes=nxNodesIn;
      this->local_nxNodes=nodeIn.local_nxNodes;
      this->oneOverdx = oneOverdxIn;

      this->nyNodes=nyNodesIn;
      this->local_nyNodes=nodeIn.local_nyNodes;
      this->oneOverdy = oneOverdyIn;

      this->nzNodes=nzNodesIn;
      this->local_nzNodes=nodeIn.local_nzNodes;
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
      this->Fsx.newsize(Vec::LOCAL,this->local_nzNodes*this->local_nyNodes*(this->local_nxNodes+1));
      this->Fsy.newsize(Vec::LOCAL,this->local_nzNodes*(this->local_nyNodes+1)*this->local_nxNodes);
      this->Fsz.newsize(Vec::LOCAL,(this->local_nzNodes+1)*this->local_nyNodes*this->local_nxNodes);
      this->Fsx = 0.0;
      this->Fsy = 0.0;
      this->Fsz = 0.0;
      for (int vi= 0;vi<nv;vi++)
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
          for (int vj= 0;vj<nv;vj++)
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

template<class FLUX_K, class FLUX_F, class BC, int nv>
Div_KgradPlusF3d<FLUX_K,FLUX_F,BC,nv>::~Div_KgradPlusF3d(){}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF3d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                              const Vec* P)
{
  //ghost faces
  if (this->node.ghost_zOffSet)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(-this->node.ghost_zOffSet,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_z(K,Fz,P);
        }
  if (this->node.ghost_nzNodes -  this->node.ghost_zOffSet - this->node.local_nzNodes)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_z(K,Fz,P);
        }
  if (this->node.ghost_yOffSet)
    for (int i= 0;i<this->local_nzNodes;i++)
      for(int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,-this->node.ghost_yOffSet,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
        }
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for (int i= 0;i<this->local_nzNodes;i++)
      for(int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
        }
  if (this->node.ghost_xOffSet)
    for (int i= 0;i<this->local_nzNodes;i++)
      for (int j= 0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,-this->node.ghost_xOffSet);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
        }      
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int i= 0;i<this->local_nzNodes;i++)
      for (int j= 0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,this->local_nxNodes-1);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
        }
  
  //interior
  for (int i= 0;i<this->local_nzNodes-1;i++)
    {
      for (int j= 0;j<this->local_nyNodes-1;j++)
        {
          for (int k = 0;k<this->local_nxNodes-1;k++)
            {
              this->node.localIndex(i,j,k);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_z(K,Fz,P);
            }
          // last y and z fluxes of x row
          this->node.localIndex(i,j,this->local_nxNodes-1);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_z(K,Fz,P);
        }
      //last row of x and z fluxes in y plane
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_z(K,Fz,P);
        }
      // last z fluz in x row
      this->node.localIndex(i,this->local_nyNodes-1,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_z(K,Fz,P);
    }
  //top face of x and y fluxes
  for (int j= 0;j<this->local_nyNodes-1;j++)
    {
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
        }
      //last y flux in x row
      this->node.localIndex(this->local_nzNodes-1,j,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_y(K,Fy,P);
    }
  //last row of x fluxes in top face
  for (int k = 0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nzNodes-1,this->local_nyNodes-1,k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFlux_x(K,Fx,P);
    }
  
  for (int vi= 0;vi<nv;vi++)
    this->bc[vi].applyNeumannConditions(this->node,this->flux_x[vi], this->flux_y[vi], this->flux_z[vi]);
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF3d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceK(const Vec& K)
{
  //ghost faces
  if (this->node.ghost_zOffSet)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(-this->node.ghost_zOffSet,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsz(K);
        }
  if (this->node.ghost_nzNodes -  this->node.ghost_zOffSet - this->node.local_nzNodes)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsz(K);
        }
  if (this->node.ghost_yOffSet)
    for (int i= 0;i<this->local_nzNodes;i++)
      for(int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,-this->node.ghost_yOffSet,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
        }
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for (int i= 0;i<this->local_nzNodes;i++)
      for(int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
        }
  if (this->node.ghost_xOffSet)
    for (int i= 0;i<this->local_nzNodes;i++)
      for (int j= 0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,-this->node.ghost_xOffSet);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
        }      
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int i= 0;i<this->local_nzNodes;i++)
      for (int j= 0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,this->local_nxNodes-1);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
        }
  
  //interior
  for (int i= 0;i<this->local_nzNodes-1;i++)
    {
      for (int j= 0;j<this->local_nyNodes-1;j++)
        {
          for (int k = 0;k<this->local_nxNodes-1;k++)
            {
              this->node.localIndex(i,j,k);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsz(K);
            }
          // last y and z fluxes of x row
          this->node.localIndex(i,j,this->local_nxNodes-1);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsz(K);
        }
      //last row of x and z fluxes in y plane
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsz(K);
        }
      // last z fluz in x row
      this->node.localIndex(i,this->local_nyNodes-1,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsz(K);
    }
  //top face of x and y fluxes
  for (int j= 0;j<this->local_nyNodes-1;j++)
    {
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
        }
      //last y flux in x row
      this->node.localIndex(this->local_nzNodes-1,j,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsy(K);
    }
  //last row of x fluxes in top face
  for (int k = 0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nzNodes-1,this->local_nyNodes-1,k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setKsx(K);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF3d<FLUX_K,FLUX_F,BC,nv>::computeInterfaceF(const Vec& Fx, const Vec& Fy, const Vec& Fz)
{
  //ghost faces
  if (this->node.ghost_zOffSet)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(-this->node.ghost_zOffSet,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsz(Fz);
        }
  if (this->node.ghost_nzNodes -  this->node.ghost_zOffSet - this->node.local_nzNodes)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsz(Fz);
        }
  if (this->node.ghost_yOffSet)
    for (int i= 0;i<this->local_nzNodes;i++)
      for(int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,-this->node.ghost_yOffSet,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
        }
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for (int i= 0;i<this->local_nzNodes;i++)
      for(int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
        }
  if (this->node.ghost_xOffSet)
    for (int i= 0;i<this->local_nzNodes;i++)
      for (int j= 0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,-this->node.ghost_xOffSet);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
        }      
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int i= 0;i<this->local_nzNodes;i++)
      for (int j= 0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,this->local_nxNodes-1);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
        }
  
  //interior
  for (int i= 0;i<this->local_nzNodes-1;i++)
    {
      for (int j= 0;j<this->local_nyNodes-1;j++)
        {
          for (int k = 0;k<this->local_nxNodes-1;k++)
            {
              this->node.localIndex(i,j,this->local_nxNodes-1);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsz(Fz);
            }
          // last y and z fluxes of x row
          this->node.localIndex(i,j,this->local_nxNodes-1);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsz(Fz);
        }
      //last row of x and z fluxes in y plane
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsz(Fz);
        }
      // last z fluz in x row
      this->node.localIndex(i,this->local_nyNodes-1,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsz(Fz);
    }
  //top face of x and y fluxes
  for (int j= 0;j<this->local_nyNodes-1;j++)
    {
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
        }
      //last y flux in x row
      this->node.localIndex(this->local_nzNodes-1,j,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsy(Fy);
    }
  //last row of x fluxes in top face
  for (int k = 0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nzNodes-1,this->local_nyNodes-1,k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setFsx(Fx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF3d<FLUX_K,FLUX_F,BC,nv>::computeDivergence(const Vec* K, 
                                                              const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                              const Vec* P, 
                                                              Vec* Div)
{
  computeDivergence(K,Fx,Fy,Fz,P);
  for (int i= 0;i<this->local_nzNodes;i++)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,j,k);
          for (int vi= 0;vi<nv;vi++)
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
  
template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF3d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                                 const Vec* P, 
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz)
{
  //ghost faces
  if (this->node.ghost_zOffSet)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(-this->node.ghost_zOffSet,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_z(K,Fz,P,DK,DFz);
        }
  if (this->node.ghost_nzNodes -  this->node.ghost_zOffSet - this->node.local_nzNodes)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_z(K,Fz,P,DK,DFz);
        }
  if (this->node.ghost_yOffSet)
    for (int i= 0;i<this->local_nzNodes;i++)
      for(int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,-this->node.ghost_yOffSet,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
        }
  if (this->node.ghost_nyNodes - this->node.ghost_yOffSet - this->node.local_nyNodes)
    for (int i= 0;i<this->local_nzNodes;i++)
      for(int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
        }
  if (this->node.ghost_xOffSet)
    for (int i= 0;i<this->local_nzNodes;i++)
      for (int j= 0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,-this->node.ghost_xOffSet);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
        }      
  if (this->node.ghost_nxNodes - this->node.ghost_xOffSet - this->node.local_nxNodes)
    for (int i= 0;i<this->local_nzNodes;i++)
      for (int j= 0;j<this->local_nyNodes;j++)
        {
          this->node.localIndex(i,j,this->local_nxNodes-1);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
        }
  
  for (int i= 0;i<this->local_nzNodes-1;i++)
    {
      for (int j= 0;j<this->local_nyNodes-1;j++)
        {
          for (int k = 0;k<this->local_nxNodes-1;k++)
            {
              //interior
              this->node.localIndex(i,j,k);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
              Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_z(K,Fz,P,DK,DFz);
            }
          // last y and z fluxes of x row
          this->node.localIndex(i,j,this->local_nxNodes-1);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_z(K,Fz,P,DK,DFz);
        }
      //last row of x and z fluxes in y plane
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(i,this->local_nyNodes-1,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_z(K,Fz,P,DK,DFz);
        }
      // last z fluz in x row
      this->node.localIndex(i,this->local_nyNodes-1,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_z(K,Fz,P,DK,DFz);
    }
  //top face of x and y fluxes
  for (int j= 0;j<this->local_nyNodes-1;j++)
    {
      for (int k = 0;k<this->local_nxNodes-1;k++)
        {
          this->node.localIndex(this->local_nzNodes-1,j,k);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
          Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
        }
      //last y flux in x row
      this->node.localIndex(this->local_nzNodes-1,j,this->local_nxNodes-1);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_y(K,Fy,P,DK,DFy);
    }
  //last row of x fluxes in top face
  for (int k = 0;k<this->local_nxNodes-1;k++)
    {
      this->node.localIndex(this->local_nzNodes-1,this->local_nyNodes-1,k);
      Div_KgradPlusF<FLUX_K,FLUX_F,BC,nv>::setDFlux_x(K,Fx,P,DK,DFx);
    }
}

template<class FLUX_K, class FLUX_F, class BC, int nv>
void Div_KgradPlusF3d<FLUX_K,FLUX_F,BC,nv>::computeDivergenceJac(const Vec* K, 
                                                                 const Vec* Fx, const Vec* Fy, const Vec* Fz, 
                                                                 const Vec* P, 
                                                                 const VecVecVec& DK, 
                                                                 const VecVecVec& DFx, const VecVecVec& DFy, const VecVecVec& DFz,
                                                                 VecVecVecVec& DivJac)
{
  computeDivergenceJac(K,Fx,Fy,Fz,P,DK,DFx,DFy,DFz);
  for (int i= 0;i<this->local_nzNodes;i++)
    for (int j= 0;j<this->local_nyNodes;j++)
      for (int k = 0;k<this->local_nxNodes;k++)
        {
          this->node.localIndex(i,j,k);
          for (int vi= 0;vi<nv;vi++)
            {
              for (int vj= 0;vj<nv;vj++)
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
  for (int vi= 0;vi<nv;vi++)
    {
      for (int vj= 0;vj<nv;vj++)
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
