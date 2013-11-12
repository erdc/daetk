#ifndef DAETKPETSCNUMERICALJACOBIAN_H
#define DAETKPETSCNUMERICALJACOBIAN_H

#include "Definitions.h"
#include "IntVec.h"
#include "Utilities.h"
#include "NumericalJacobian.h"
#include "DaetkPetscMat.h"
#include "DaetkPetscSys.h"
#include "DaetkPetscAnalyticalJacobian.h"
#include <vector>

namespace Daetk 
{
namespace Petsc
{
  namespace cc
    {
      struct _p_Mat;
      struct _p_Vec;
    }
class JacobianBase : public Daetk::NumericalJacobian
{
public:
  JacobianBase(Petsc::Mat& M, VectorFunction& F);
  virtual ~JacobianBase();
  virtual Petsc::Mat& getMatShell();
  static int petscMatVec(cc::_p_Mat* A, cc::_p_Vec* x, cc::_p_Vec* Ax);
  typedef int (*PetscMatVecType)(cc::_p_Mat* A, cc::_p_Vec* x, cc::_p_Vec* Ax);
  static JacobianBase* theJacVec;
  Petsc::Mat matShell;
  Vec xVec,AxVec;
  Petsc::AnalyticalJacobian pajac;
};

template<class STENCIL, int nv>
class NumericalJacobian : public JacobianBase
{
 public:
  NumericalJacobian(Petsc::Mat& M,VectorFunction& F, STENCIL& stencilIn);
  virtual ~NumericalJacobian();
  bool evaluate(const Vec& x,const Vec& F);
  virtual void attachToSubSystem(VectorFunction& F, const Vec& Fatx);
private:
  int nNodes,neq,nColors,globalLow,globalHigh;

  IntVec pFlag,xFlag,xDoneFlag;
  Vec localDelta,localF,localFpd;
  STENCIL& stencil;
  Err ierr;
  Petsc::Mat& matrix; 
  inline void markPoints(STENCIL& s, IntVec& p);

  bool hasColoring;
  std::vector< std::vector<int> > coloringVector;
  std::vector< std::vector<int> > coloringVectorWithGhost;
  Vec localDeltaAttache;
};
  
template <class STENCIL,int nv>
NumericalJacobian<STENCIL,nv>::NumericalJacobian(Petsc::Mat& M,VectorFunction& F, STENCIL& stencilIn):
    JacobianBase(M,F),
    nNodes(stencilIn.nNodes()),
    neq(M.dimDomain()),
    nColors(0),
    pFlag(M.dimDomain()),
    xFlag(M.dimDomain()),
    xDoneFlag(M.dimDomain()),
    localDelta(Vec::LOCAL,stencilIn.dadof_),
    localF(Vec::LOCAL,stencilIn.dadof_),
    localFpd(Vec::LOCAL,stencilIn.dadof_),
    stencil(stencilIn),
    matrix(M),
    hasColoring(false)
  {}
template <class STENCIL,int nv>
void NumericalJacobian<STENCIL,nv>::attachToSubSystem(VectorFunction& F, const Vec& Fatx)
{
  Daetk::NumericalJacobian::attachToSubSystem(F,Fatx);
  if(!SOLVE_SUB)
    {
      VecIndex indexAll;
      localDeltaAttache.attachToVecMulti(Vec::REF,localDelta,indexAll);
    }
  else
    {
      localDeltaAttache.attachToVecMulti(Vec::REF,localDelta,index);
      localDeltaAttache.setStrideMulti(str);
    }
}

template <class STENCIL,int nv>
NumericalJacobian<STENCIL,nv>::~NumericalJacobian(){}

template <class STENCIL,int nv>
inline void NumericalJacobian<STENCIL,nv>::markPoints(STENCIL& s, IntVec& p)
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
bool NumericalJacobian<STENCIL,nv>::evaluate(const Vec& x,const Vec& Fatx)
{
  if (USE_ANALYTICAL_JACOBIAN)
    return pajac.evaluate(x,Fatx);
  matrix.zeroAll();
  attachToSubSystem(*Fp,Fatx);
  localDelta.startSetFromGlobalMulti(Fp->deltaVF);

  //these are for the vectors, not necessarily the grid (it could have m components)
  globalLow = x.getGlobalLow();
  globalHigh = x.getGlobalHigh();
  
  if (!hasColoring)
    {  
      xDoneFlag=0;
      int checkSum=0;
      while(checkSum<neq)
        {
          nColors++;
          coloringVector.resize(nColors);
          std::vector<int>& color = coloringVector.back();

          coloringVectorWithGhost.resize(nColors);
          std::vector<int>& colorWithGhost = coloringVectorWithGhost.back();

          //          tempDelta=0.0;
          xFlag=0;
          pFlag=0;
          for (int j=0;j<neq;j++)
            {
              stencil.setGlobalPointList(j/nv);    
              int petscJ=stencil.anchor->globalNodeNumber*nv + j%nv;              
              if (xDoneFlag(petscJ) == 0 && pFlag(petscJ) == 0)
                {
                  checkSum++;
                  xFlag(petscJ) = 1;
                  xDoneFlag(petscJ) = 1;
                  
                  //if x_j affects the equation at another node mark off all the x's that would affect 
                  //that equation (including the x's that affect this equation)
                  //markPoints(stencil,pFlag);
                  //I reposition the stencil inside the loop
                  //so I have to save my end value (stencil.end() will
                  //change inside the loop
                  typename STENCIL::iterator sit=stencil.begin();
                  const typename STENCIL::iterator end=stencil.end();
                  while (sit != end)
                    {
                      if ((sit->globalNodeNumber >= stencil.globalLow) && 
                          (sit->globalNodeNumber < stencil.globalHigh))
                        colorWithGhost.push_back(petscJ);
                      
                      stencil.setGlobalPointList(sit->i,sit->j,sit->k);
                      markPoints(stencil,pFlag);
                      ++sit;
                    }
                  if ((petscJ >= globalLow) && (petscJ < globalHigh))
                    //tempDelta(petscJ) = Fp->deltaVF(petscJ);
                    color.push_back( petscJ);
                }
            } 
        }
      hasColoring = true;
    }

  std::vector<int>::iterator index,indexEnd;
  localDelta.endSetFromGlobalMulti(Fp->deltaVF);
  
  int i,j,local_i,local_j;
  
  for (int c=0;c<nColors;c++)
    {
      tempDelta=0.0;

      index=coloringVector[c].begin();
      indexEnd=coloringVector[c].end();

      while (index < indexEnd)
        {
          //use local indexing because it's faster
          tempDeltaAttache[*index - globalLow] = deltaAttache[*index - globalLow];
          ++index;
        }
      
      Fp->correctArgument(tempDelta);
      
      bool evalError=false;
      FatxPdelta = Fp->value(evalError);
      
      int trys=0;
      while (evalError && trys < 10)
        {
          trys++;
          Fp->unCorrect();
          std::cerr<<"Delta in numerical jacobian is causing S or P to be out of range"<<std::endl;
          (localDelta)*=0.1;
          (tempDelta)*=0.1;
          (Fp->deltaVF)*=0.1;
          Fp->correctArgument(tempDelta);
          
         FatxPdelta = Fp->value(evalError);
        }
      if (evalError)
        return evalError;

      int globalNode;
      real tempDelInverse;

      index=coloringVectorWithGhost[c].begin();
      indexEnd = coloringVectorWithGhost[c].end();
      while (index < indexEnd)      
        {
          j = (*index);

          globalNode = stencil.petscToGlobal( j/nv );
          stencil.setGlobalPointList( globalNode );
          //local_j needs to be used on a local vector with ghost regions
          local_j = stencil.globalToLocal(globalNode)*nv + j%nv;

          tempDelInverse = -1.0/localDeltaAttache[local_j];
          
          typename STENCIL::iterator sit=stencil.begin();
          const typename STENCIL::iterator end=stencil.end();
          while (sit != end)
            {
              if ( (sit->globalNodeNumber >= stencil.globalLow) &&
                   (sit->globalNodeNumber < stencil.globalHigh ) )
                {
                  for (int vi=0;vi<nv;vi++) 
                    {
                      i=sit->globalNodeNumber*nv+vi;
                      //local_i is just a local index for a slice out of a global vector
                      local_i = i - globalLow;

                      matrix(i,j) = tempDelInverse*
                        ( FatxPdeltaAttache[local_i] - FatxAttache[local_i] );
                      matrix.finalizeRow(i);//real slow
                    }
                }
              ++sit;
            }
          ++index;
        }	
      Fp->unCorrect();
    }
  
  matrix.beginAssembly();
  matrix.endAssembly();
  return false;
}

}//Petsc
}//Daetk
#endif
