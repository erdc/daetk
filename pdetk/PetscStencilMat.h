#ifndef PETSCSTENCILMAT_H
#define PETSCSTENCILMAT_H

#include "Definitions.h"
#include "LinearOperator.h"
#include "Vec.h"

namespace Daetk 
{
namespace Petsc 
{
typedef CMRVec<Vec> VecVec;

//stuff in Stencil base class
//  struct Location 
//  {
//    int i,j,k,globalNodeNumber,point;
//  }
//  typedef iterator  //iterator for a container of Locations
//  iterator begin(); //return the beginning of Loc. container
//  iterator end();   //return 1 past end of the Loc. container

//  Required STENCIL  Functionality
//  int size() //return the maximum number of points in the stencil
//  int nNodes() //return the number of nodes/unknowns in the domain
//  int globalNode(int i) //send stencil to global node i
//  enum Point {} //an enum of the stencil points

template <class STENCIL, int nv>
class StencilMat;

template <class STENCIL, int nv>
std::ostream&  operator<<(std::ostream& s, const StencilMat<STENCIL,nv>& m);

template <class STENCIL, int nv>
class StencilMat : public Daetk::LinearOperator
{
public:
  StencilMat(STENCIL& s);
  
  StencilMat(int,int,int dim,STENCIL& s);
  
  virtual ~StencilMat();
  
  StencilMat& operator=(const real& scalar);

  real& operator()(int iNode, int jNode, int vi=0,int vj=0);
  
  //this one assumes a global node number and 
  //and a point on stencil(i)
  //does not assume lexicogrphic ordering of global nodes
  inline real& operator()(int iNode, typename STENCIL::Point p, int vj=0,int vk=0);
  inline Vec& operator[](int i);
  inline void zeroRow(int i);
  inline real& operator()(typename STENCIL::Location& si, 
                          typename STENCIL::Location& sj, int vi=0, int vj=0);

  bool apply(const Vec& x, Vec& Ax);
  friend std::ostream&  operator<<  <> (std::ostream& s, const StencilMat& m);

  inline int dim(int i=0);
  inline int stencilSize();
  inline STENCIL& getStencil();
private:
  int nPoints,nNodes;
  STENCIL& stencil;
  VecVec mat;
};

template <class STENCIL, int nv>
inline std::ostream&  operator<<(std::ostream& s, const Petsc::StencilMat<STENCIL,nv>& m)
{
  s.setf(std::ios::scientific);
  s.precision(17);
  for (int i=0;i<m.nNodes;i++)
    {
      m.stencil.globalNode(i);
      typename STENCIL::iterator sit = m.stencil.begin();
      const typename STENCIL::iterator end=m.stencil.end();
      while (sit != end)
        {
          for (int vi=0; vi < nv; vi++)
            for (int vj=0; vj < nv; vj++)
              s<<"("<<i*nv + vi<<","<<sit->second.globalNodeNumber*nv+vj<<")="<<m.mat[i*nv+vi][sit->second.point*nv+vj]<<'\t';
          ++sit;
        }
      s<<std::endl;
    }
  return s;
}

template <class STENCIL, int nv>
StencilMat<STENCIL,nv>::StencilMat(STENCIL& s):
  Daetk::LinearOperator(s.nNodes()*nv,s.nNodes()*nv),
         nPoints(7),
         nNodes(s.nNodes()),
         stencil(s),
         mat(s.nNodes()*nv)
{
  for (int i=0;i<nNodes;i++)
    for (int vi=0; vi < nv ; vi++)
      {
        mat[i*nv + vi].newsize(nPoints*nv);
        mat[i*nv + vi]=0.0;
      }
}

template <class STENCIL, int nv>
StencilMat<STENCIL,nv>::StencilMat(int,int,int dim,STENCIL& s):
  Daetk:: LinearOperator(s.nNodes()*nv,s.nNodes()*nv),
          nPoints(7),
          nNodes(s.nNodes()),
          stencil(s),
          mat(dim*nv)
{
  for (int i=0;i<nNodes;i++)
    for (int vi=0; vi < nv ; vi++)
      {
        mat[i*nv + vi].newsize(nPoints*nv);
        mat[i*nv + vi]=0.0;
      }
}

template <class STENCIL, int nv>
StencilMat<STENCIL,nv>::~StencilMat(){}

template <class STENCIL, int nv>
StencilMat<STENCIL,nv>& StencilMat<STENCIL,nv>::operator=(const real& scalar)
{
  for (int i=0;i<nNodes;i++)
    for (int vi=0; vi < nv ; vi++)
      mat[i*nv + vi]=scalar;
  return *this;
}

template <class STENCIL, int nv>
inline real& StencilMat<STENCIL,nv>::operator()(int iNode, int jNode, int vi,int vj)
{
  //index checking should be put first, if at all
  stencil.globalNode(iNode);
  typename STENCIL::iterator sit = stencil.begin();
  const typename STENCIL::iterator end=stencil.end();
  while (sit != end)
    {
      //            cout<<"here"<<std::endl;
      if (sit->second.globalNodeNumber == jNode)
        return mat[iNode*nv + vi][sit->second.point*nv+vj];
      ++sit;
    }
  //        cout<<"shouldn't be here"<<std::endl;
  static real ZERO=0.0;
  return ZERO=0.0;
}

template <class STENCIL, int nv>
inline real& StencilMat<STENCIL,nv>::operator()(int iNode, typename STENCIL::Point p, int vj,int vk)
{
  return mat[iNode*nv + vj][p*nv+vk];
}

template <class STENCIL, int nv>
inline Vec& StencilMat<STENCIL,nv>::operator[](int i)
{
  return mat[i];
}

template <class STENCIL, int nv>
inline void StencilMat<STENCIL,nv>::zeroRow(int i){mat[i] = 0.0;}

template <class STENCIL, int nv>
inline real& StencilMat<STENCIL,nv>::operator()(typename STENCIL::Location& si, 
                   typename STENCIL::Location& sj, int vi, int vj)
{
  return mat[si.globalNodeNumber*nv + vi][sj.point*nv+vj];
}

template <class STENCIL, int nv>
inline bool StencilMat<STENCIL,nv>::apply(const Vec& x, Vec& Ax)
{
  Ax=0.0;
  //we should swap these loops, but it's clearer this way
  //we could also let the stencil  do the row mult
  for (int i=0;i<nNodes;i++)
    {
      stencil.globalNode(i);
      typename STENCIL::iterator sit = stencil.begin();
      const typename STENCIL::iterator end=stencil.end();
      while (sit != end)
        {
          for (int vi=0;vi<nv;vi++)
            for (int vj=0;vj<nv;vj++)
              Ax(i*nv + vi) += mat[i*nv+vi][sit->second.point*nv+vj]*x(sit->second.globalNodeNumber*nv + vj);
          ++sit;
        }
    }
  return false;
}

template <class STENCIL, int nv>
inline int StencilMat<STENCIL,nv>::dim(int i){ return nNodes*nv;}

template <class STENCIL, int nv>
inline int StencilMat<STENCIL,nv>::stencilSize(){return nPoints;}

template <class STENCIL, int nv>
inline STENCIL& StencilMat<STENCIL,nv>::getStencil(){return stencil;}

}//Petsc
}//Daetk
#endif








