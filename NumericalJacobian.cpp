#include "NumericalJacobian.h"

namespace Daetk
{
  NumericalJacobian::NumericalJacobian(VectorFunction& F):
    SOLVE_SUB(false),
    Jacobian(F),
    USE_ANALYTICAL_JACOBIAN(false),
    tempDelta(F.dimDomain()),
    FatxPdelta(F.dimRange()),
    ajac(F)
    {}

  void NumericalJacobian::setFunction(VectorFunction& F)
    {
      Fp=&F;
      FatxPdelta.newsize(F.dimRange());
    }

  NumericalJacobian::~NumericalJacobian(){}
  
  bool NumericalJacobian::apply(const Vec& x, Vec& Ax)
    {
      if (USE_ANALYTICAL_JACOBIAN)
        return ajac.apply(x,Ax);
      return Fp->numericalJacVec(x,Ax);
    }

  void NumericalJacobian::solveSubSystem(int start,int end,int stride,
                                         int dimLS)
  {SOLVE_SUB=true; VecIndex i(start,end); index=i; str=stride;}
  
  void NumericalJacobian::attachToSubSystem(VectorFunction& F, const Vec& Fatx)
  {
    if(!SOLVE_SUB)
      {
        VecIndex indexAll;
        tempDeltaAttache.attachToVecMulti(Vec::REF,tempDelta,indexAll);
        deltaAttache.attachToVecMulti(Vec::REF,Fp->deltaVF,indexAll);
        FatxPdeltaAttache.attachToVecMulti(Vec::REF,FatxPdelta,indexAll);
        FatxAttache.attachToVecMulti(Vec::REF,Fatx,indexAll);
      }
    else
      {
        tempDeltaAttache.attachToVecMulti(Vec::REF,tempDelta,index);
        tempDeltaAttache.setStrideMulti(str);
        
        deltaAttache.attachToVecMulti(Vec::REF,Fp->deltaVF,index);
        deltaAttache.setStrideMulti(str);
        
        FatxPdeltaAttache.attachToVecMulti(Vec::REF,FatxPdelta,index);
        FatxPdeltaAttache.setStrideMulti(str);
        
        FatxAttache.attachToVecMulti(Vec::REF,Fatx,index);
        FatxAttache.setStrideMulti(str);
      }
  }
  
void NumericalJacobian::useAnalyticalJacobian(bool flag)
{
  USE_ANALYTICAL_JACOBIAN=flag;
}

}//Daetk
