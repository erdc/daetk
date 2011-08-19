#ifndef STENCILNUMERICALJACOBIAN_H
#define STENCILNUMERICALJACOBIAN_H

#include "Definitions.h"
#include "IntVec.h"
#include "Utilities.h"
#include "NumericalJacobian.h"
#include "StencilMat.h"

namespace Daetk 
{
template<class STENCIL, int nv>
class StencilNumericalJacobian : public NumericalJacobian
{
 public:
  virtual ~StencilNumericalJacobian();
  StencilNumericalJacobian(StencilMat<STENCIL,nv>& M,VectorFunction& F);

  bool evaluate(const Vec& x,const Vec& F);
private:
  int nNodes,neq;
  IntVec pFlag,xFlag,xDoneFlag;
  Vec tempDelta;
  STENCIL& stencil;
  StencilMat<STENCIL,nv>& matrix; 
  inline void markPoints(STENCIL& s, IntVec& p);
};

template <class STENCIL,int nv>
inline  void StencilNumericalJacobian<STENCIL,nv>::markPoints(STENCIL& s, IntVec& p)
{        
  typename STENCIL::iterator sit=stencil.begin();
  const typename STENCIL::iterator end=stencil.end();
  while (sit != end)
    {
      for (int vj=0;vj<nv;vj++)
        p(sit->globalNodeNumber*nv + vj) = 1;
      ++sit;
    }
}

template <class STENCIL,int nv>
bool StencilNumericalJacobian<STENCIL,nv>::evaluate(const Vec& x,const Vec& Fatx)
{
  matrix=0.0;
  int checkSum=0;
  int j;
  xDoneFlag=0.0;
  while(checkSum<neq)
    {
      tempDelta=0.0;
      xFlag=0;
      pFlag=0;
      for (j=0;j<neq;j++)
        {
          if (xDoneFlag(j) == 0 && pFlag(j) == 0)
            {
              stencil.globalNode(j/nv);
              checkSum++;
              xFlag(j) = 1;
              xDoneFlag(j) = 1;

              //if x_j affects the equation at another node mark off all the x's that would affect 
              //that equation (including the x's that affect this equation)
              markPoints(stencil,pFlag);
              //I reposition the stencil inside the loop
              //so I have to save my end value (stencil.end() will
              //change inside the loop
              typename STENCIL::iterator sit=stencil.begin();
              const typename STENCIL::iterator end=stencil.end();
              while (sit != end)
                {
                  //cek changed from () to globalNode. I don't think it should be () (that's a 1d index)
                  stencil.globalNode(sit->globalNodeNumber);
                  markPoints(stencil,pFlag);
                  ++sit;
                }
              tempDelta(j) = Fp->deltaVF(j);
            }
        }
      Fp->correctArgument(tempDelta);
      bool evalError=false;
      FatxPdelta = Fp->value(evalError);
      if (evalError)
        {
          Fp->unCorrect();
          cerr<<"Delta in numerical jacobian is causing S or P to be out of range"<<endl;
          return true;
        }
      for (j=0;j<neq;j++)
	{
	  real tempDelInverse;
	  if (xFlag(j) == 1)
	    {
              stencil.globalNode(j/nv);
              int vj=j%nv;
              //minus is because deltap is really in the other direction
              tempDelInverse = -1.0/tempDelta(j);
              typename STENCIL::iterator sit=stencil.begin();
              const typename STENCIL::iterator end=stencil.end();
              while (sit != end)
                {
                  for (int vi=0;vi<nv;vi++) 
                    {
                      matrix(sit->globalNodeNumber,j/nv,vi,vj)=
                        tempDelInverse*(FatxPdelta(sit->globalNodeNumber*nv+vi)-Fatx(sit->globalNodeNumber*nv+vi));
//                        if (j/nv == 0 && j%nv == 1)
//                          {
//                              cout<<(sit->globalNodeNumber*nv+vi)<<'\t'<<vi<<'\t'<<vj<<'\t'<<j/nv<<'\t'<<tempDelInverse*(FatxPdelta(sit->globalNodeNumber*nv+vi)-Fatx(sit->globalNodeNumber*nv+vi))<<'\t'<<matrix(sit->globalNodeNumber,j/nv,vi,vj)<<endl;
//                          }
//                        if (j==2)
//                          {
//                            cout<<sit->globalNodeNumber<<'\t'<<j/nv<<'\t'<<vi<<'\t'<<vj<<endl;
//                            cout<<(sit->globalNodeNumber*nv+vi)<<'\t'<<vi<<'\t'<<vj<<'\t'<<j/nv<<'\t'<<tempDelInverse*(FatxPdelta(sit->globalNodeNumber*nv+vi)-Fatx(sit->globalNodeNumber*nv+vi))<<'\t'<<matrix(sit->globalNodeNumber,j/nv,vi,vj)<<endl;
//                          }
                    }
                  ++sit;
                }
//                cout<<j<<'\t'<<matrix(0,0,1,1)<<endl;
	    }
	}
      Fp->unCorrect();
    }
  return false;
}
template <class STENCIL,int nv>
StencilNumericalJacobian<STENCIL,nv>::~StencilNumericalJacobian(){}

template <class STENCIL,int nv>
StencilNumericalJacobian<STENCIL,nv>::StencilNumericalJacobian(StencilMat<STENCIL,nv>& M,VectorFunction& F):
    NumericalJacobian(F),
    nNodes(M.getStencil().nNodes()),
    neq(M.dim()),
    pFlag(M.dim()),
    xFlag(M.dim()),
    xDoneFlag(M.dim()),
    tempDelta(F.dimDomain()),
    stencil(M.getStencil()),
    matrix(M)
    {};

}//Daetk
#endif
