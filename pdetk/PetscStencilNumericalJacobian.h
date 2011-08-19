#ifndef PETSCSTENCILNUMERICALJACOBIAN_H
#define PETSCSTENCILNUMERICALJACOBIAN_H

#include "Definitions.h"
#include "IntVec.h"
#include "Utilities.h"
#include "NumericalJacobian.h"
#include "PetscStencilMat.h"

namespace Daetk 
{
namespace Petsc
{
template<class STENCIL, int nv>
class StencilNumericalJacobian : public Daetk::NumericalJacobian
{
 public:
  StencilNumericalJacobian(Petsc::StencilMat<STENCIL,nv>& M,VectorFunction& F, STENCIL& stencilIn);

  bool evaluate(const Vec& x,const Vec& F);
private:
  int nNodes,neq;
  IntVec pFlag,xFlag,xDoneFlag;
  Vec tempDelta,localDelta,localF,localFpd;
  STENCIL& stencil;
  Petsc::StencilMat<STENCIL,nv>& matrix; 
  inline void markPoints(STENCIL& s, IntVec& p);
};

template <class STENCIL,int nv>
inline void Petsc::StencilNumericalJacobian<STENCIL,nv>::markPoints(STENCIL& s, IntVec& p)
{        
  typename STENCIL::iterator sit=stencil.begin();
  const typename STENCIL::iterator end=stencil.end();
  while (sit != end)
    {
      for (int vj=0;vj<nv;vj++)
        p(sit->second.globalNodeNumber*nv + vj) = 1;
      ++sit;
    }
}

template <class STENCIL,int nv>
bool Petsc::StencilNumericalJacobian<STENCIL,nv>::evaluate(const Vec& x,const Vec& Fatx)
{
  localF.startSetFromGlobal(Fatx);
  localF.endSetFromGlobal(Fatx);
  //  matrix=0.0;
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
          stencil.globalNode(j/nv);    
          int petscJ=stencil.center*nv + j%nv;              
          if (xDoneFlag(petscJ) == 0 && pFlag(petscJ) == 0)
            {
              checkSum++;
              xFlag(petscJ) = 1;
              xDoneFlag(petscJ) = 1;

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
                  //cek changed to globalNode() from () don't think it should be a 1d index
                  stencil.globalNode(sit->second.globalNodeNumber);
                  markPoints(stencil,pFlag);
                  ++sit;
                }
              if (petscJ >= tempDelta.getGlobalLow() && petscJ < tempDelta.getGlobalHigh())
                tempDelta(petscJ) = Fp->deltaVF(petscJ);
            }
        }
      Fp->correctArgument(tempDelta);
      bool evalError=false;
      FatxPdelta = Fp->value(evalError);
      localFpd.startSetFromGlobal(FatxPdelta);
      localFpd.endSetFromGlobal(FatxPdelta);
      
      if (evalError)
        {
          Fp->unCorrect();
          cerr<<"Delta in numerical jacobian is causing S or P to be out of range"<<endl;
          return true;
        }
      int petscJ;
      for (j=0;j<neq;j++)
	{
          int oldj;
          stencil.globalNode(j/nv);
          oldj = stencil.center;
          petscJ=stencil.center*nv + j%nv;              
          int vj=petscJ%nv;

          //cout<<"in set"<<endl;
          real tempDelInverse;
          if (xFlag(petscJ) == 1 && petscJ >= tempDelta.getGlobalLow() && petscJ < tempDelta.getGlobalHigh())
            {
              tempDelta(petscJ) = Fp->deltaVF(petscJ);
              //minus is because deltap is really in the other direction
              tempDelInverse = -1.0/tempDelta(petscJ);
              typename STENCIL::iterator sit=stencil.begin();
              const typename STENCIL::iterator end=stencil.end();
              while (sit != end)
                {
                  for (int vi=0;vi<nv;vi++) 
                    {
//                        int i=sit->second.globalNodeNumber*nv+vi;
 //                       stencil.localIndex(sit->second.i-stencil.local_z0,
//                                           sit->second.j-stencil.local_y0,
//                                           sit->second.k-stencil.local_x0);
//                        cout<<sit->first<<'\t'<<nv<<'\t'<<i<<'\t'<<petscJ<<endl<<flush;
                      matrix(sit->second.globalNodeNumber,petscJ/nv,vi,vj) =
                        tempDelInverse*(FatxPdelta(sit->second.globalNodeNumber*nv+vi)-Fatx(sit->second.globalNodeNumber*nv+vi));
                    }
                  ++sit;
                }
            }
          //cout<<"out of set"<<endl;
	}
      Fp->unCorrect();
    }
  Petsc::Sys psys;
  //cout<<matrix<<endl;
  return false;
}


template <class STENCIL,int nv>
StencilNumericalJacobian<STENCIL,nv>::StencilNumericalJacobian(Petsc::StencilMat<STENCIL,nv>& M,VectorFunction& F, STENCIL& stencilIn):
  Daetk::NumericalJacobian(F),
         nNodes(stencilIn.nNodes()),
         neq(M.dimDomain()),
         pFlag(M.dimDomain()),
         xFlag(M.dimDomain()),
         xDoneFlag(M.dimDomain()),
         tempDelta(F.dimDomain()),
         localDelta(Vec::LOCAL,stencilIn.dadof_),
         localF(Vec::LOCAL,stencilIn.dadof_),
         localFpd(Vec::LOCAL,stencilIn.dadof_),
         stencil(stencilIn),
         matrix(M)
{};

}//Petsc
}//Daetk
#endif
