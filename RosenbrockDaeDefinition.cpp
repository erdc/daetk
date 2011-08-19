#include "RosenbrockDaeDefinition.h"
#include <iostream>

#include "DataCollector.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk
{

RosenbrockDaeDefinition::RosenbrockDaeDefinition(const int& neqIn,
						 DataCollector* dataIn):
  VectorFunction(neqIn,neqIn),
  data(dataIn),
  tDaeDef(0.),
  yDaeDef(neqIn,12345.),
  Fcurrent(neqIn,12345.),
  yLast(neqIn,12345.),
  Flast(neqIn,12345.),
  updateF(true),
  updateJac(true),
  updateJacT(true),
  computeOwnJacobianDelta(false),
  del(neqIn,12345.)
{

}

RosenbrockDaeDefinition::~RosenbrockDaeDefinition()
{
}

bool RosenbrockDaeDefinition::ok()
{
  bool isOk = true;
  isOk = isOk && data;

  return isOk;
}
//======================================================================
//Basic methods for use with Jacobians and other solvers
//======================================================================
bool RosenbrockDaeDefinition::resetFunction()
{
  updateF   = true;
  updateJac = true;
  updateJacT= true;

  return false;

}
//======================================================================
//Default methods for use with Jacobians and other solvers
//======================================================================
//calculate \mat{M}
//It is very important that this not overwrite storage
//should set jac --> jac - gamDtInv*M
bool 
RosenbrockDaeDefinition::appendMinusMassMatrix(const real& t, const Vec& y,
					       const real& gamDtInv)
{

  std::cerr<<"You must override appendMinusMassMatrix to use a"
	   <<" non-identity mass matrix"<<std::endl;
  return true;

}

bool 
RosenbrockDaeDefinition::appendMinusDMassDyDotZ(const real& t, const Vec& y,
				       const Vec& z)
{

  std::cerr<<"You must override appendMinusDMassDyDotZ to use a"
	   <<" non-constant mass matrix"<<std::endl;
  return true;

}
bool 
RosenbrockDaeDefinition::appendMinusDMassDtDotZ(const real& t, const Vec& y,
				       const Vec& z, Vec& Jt)
{

  std::cerr<<"You must override appendMinusDMassDtDotZ to use a"
	   <<" non-constant in time mass matrix"<<std::endl;
  return true;

}



bool 
RosenbrockDaeDefinition::applyMassMatrix(const real& t, const Vec& y, 
					 const Vec& x, Vec& Mx)
{
  std::cerr<<"You must override applyMassMatrix to use a"
	   <<" non-identity mass matrix"<<std::endl;
  return true;

}

bool RosenbrockDaeDefinition::evaluateDFDy(const real& t,const Vec& y)
{
  std::cerr<<"You must override evaluateDFDy to use an"
	   <<" analytical jacobian"<<std::endl;
  return true;

}
bool RosenbrockDaeDefinition::evaluateDFDt(const real& t,const Vec& y,
					   Vec& Jt)
{
  std::cerr<<"You must override evaluateDFDt to use an"
	   <<" analytical jacobian with respect to t"<<std::endl;
  return true;

}
bool RosenbrockDaeDefinition::jacVec(const Vec& x, Vec& Jx)
{
  std::cerr<<"You must override jacVec to use an analytical jacobian"
	   <<" with an iterative method"<<std::endl;
  
  return true;
}

bool 
RosenbrockDaeDefinition::constraintValue(const real& t, const Vec& y,
					 const Vec& yaux, Vec& cval)
{

  std::cerr<<"You must override constraintValue to use the simple "
	   <<" constraint formulation"<<std::endl;
  return true;
}

bool 
RosenbrockDaeDefinition::evaluateDConstraintDt(const real& t, const Vec& y,
					 const Vec& yaux, Vec& JCt)
{

  std::cerr<<"You must override evaluateDConstraintDt to use the simple "
	   <<" constraint formulation"<<std::endl;
  return true;
}

bool 
RosenbrockDaeDefinition::applyDConstraintDy(const real& t, const Vec& y,
					    const Vec& x, Vec& Mx)
{

  std::cerr<<"You must override applyDConstraintDy to use the simple "
	   <<" constraint formulation"<<std::endl;

  return true;
}

bool 
RosenbrockDaeDefinition::
appendMinusMassMatrixDConstraintDy(const real& t,
				   const Vec& y,
				   const real& gamDtInv)
{

  std::cerr<<"You must override appendMinusMassMatrixDConstraintDy "
	   <<" to use the simple "
	   <<" constraint formulation"<<std::endl;

  return true;
}
//======================================================================
//VectorFunction Interface
//======================================================================
const Vec& RosenbrockDaeDefinition::argument()
{
  return yDaeDef;
}

const Vec& RosenbrockDaeDefinition::value(bool& evalError)
{
  evalError = false;
  if (updateF)
    {
      data->functionEvaluation();
      evalError = rightHandSideValue(tDaeDef,yDaeDef,Fcurrent);
      if (!evalError)
	updateF=false;
    }
  return Fcurrent;
}

void RosenbrockDaeDefinition::correctArgument(Vec& correction)
{
  yLast = yDaeDef;
  Flast = Fcurrent;

  updateJac=true;
  updateF=true;
#ifndef USE_BLAS
  yDaeDef-=correction;
#else
  axpy(-1.0,correction,yDaeDef);
#endif
}

void RosenbrockDaeDefinition::unCorrect()
{
  updateF = false;
  yDaeDef = yLast;
  Fcurrent = Flast;
}

bool RosenbrockDaeDefinition::evaluateAnalyticalJacobian()
{
  return evaluateDFDy(tDaeDef,yDaeDef);
}

bool RosenbrockDaeDefinition::numericalJacVec(const Vec& v, Vec& Jv)
{
  bool evalError=false;
  real delFac=1.0e-5;
  value(evalError);//sets Fcurrent
  if (evalError)
    return evalError;
  else
    {
//cek       del = -delFac*v;
      del = v;
      scal(-delFac,del);
      //end cek
      correctArgument(del);//sets Flast to Fcurrent -- 
                          //remember correction is -del = delFac*v
      value(evalError);
      while (evalError)
        {
	  std::cerr<<"cutting back on del in numerical jac vec"<<std::endl;
          unCorrect();//sets Flast to Fcurrent
          delFac*=0.1;
//cek           del = -delFac*v;
          del = v;
          scal(-delFac,del);
          //end cek
          correctArgument(del);
          value(evalError);//sets Fcurrent
        }
//cek      Jv = (Fcurrent - Flast)/delFac;
      Jv = Fcurrent;
      axpy(-1.0,Flast,Jv);
      scal(delFac,Jv);
      //end cek
      unCorrect();
    }
  return evalError;
}

bool RosenbrockDaeDefinition::analyticalJacVec(const Vec& v, Vec& Jv)
{
  bool evalError = false;
  if (updateJac)
    {
      evalError = evaluateAnalyticalJacobian();
      if (evalError)
        return evalError;
      else
        updateJac=false;
    }
  jacVec(v,Jv);
  return evalError;
}

//mwf added to allow other integrators that don't compute delta
void RosenbrockDaeDefinition::computeDeltaForJacobian()
{

  if (computeOwnJacobianDelta)
    {
      deltaVF = 0.0;
      for (int i=0;i<yDaeDef.ldim();i++)
	{
	  //mwf this is what it should be 
	  deltaVF[i] = 
	    SQRT_MACHINE_EPSILON*fabs(yDaeDef[i])+SQRT_MACHINE_EPSILON;
	  deltaVF[i]= yDaeDef[i] - (yDaeDef[i] - deltaVF[i]);
	}
    }
}

}//end Daetk
