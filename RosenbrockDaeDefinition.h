#ifndef ROSENBROCK_DAE_DEFINITION_H
#define ROSENBROCK_DAE_DEFINITION_H

#include "Definitions.h"
#include "Vec.h"
#include "VectorFunction.h"

namespace Daetk
{
class DataCollector;

class RosenbrockDaeDefinition: public VectorFunction
{
  /***********************************************************************
  
     Initial effort to provide an interface for solving systems including
     
       \mat{M}\dot{\vec y} = \vec F(t,\vec y)

     Using Rosenbrock time integration methods. There are three basic
     formulations depending on the nature of \mat{M}.

     I.   \mat{M} is constant but maybe singular (linearly implicit ODE)
    
     II.  \mat{M} = \mat{M}(t,\vec y) (implicit ODE)

     III. \mat{M} is constant, singular with special form
          corresponding to our "extended" mass formulations where we
          have a conservation equation supplemented by a simple
          nonlinear constraint

          \pd{z}{t} = g(x,t,u)
                  z = \mu(u)

          I'll try to express the corresponding MOL'd DAE system as
 
          \mat{M}\dot{\vec y} = \vec F(t,\vec y)
     
          \vec y   = [\vec y_1, \vec y_2]^T, 
          \vec y_1 = [...,z_i,...]^T
          \vec y_2 = [...,u_i,...]^T

          \vec F   = [\vec f_1, \vec f_2]^T

          \vec f_1 = [..., \mu(u_i)-z_i, ...]^T
          \vec f_2 = [..., G(t,\vec y_2)_i, ...]^T


          \mat{M} = | 0            0|
                    |\mat{M}_{21}  0|
          
          \mat{M}_{21} constant

          \mat{F}_y = |-\mat{I}  \mat{f}_{1,y_2} |
                      |   0      \mat{f}_{2,y_2} |


     I will try to stay as close
     to the original DaeDefinition  as possible to see if I can get 
     back to a similar interface

     Initially, the main thing is that I want to setup the inherited
     VectorFunction interface so that I can use the jacobian
     interface to calculate \vec F_{\vec y}
 
     I will include an interface for calculating \vec F_{t}

     For formulation III, the routines for evaluating quantities in
     formulations I and II (like evaluateDFDt) apply to \vec f_2 and
     \vec y_2. The additional routines for the constraint and
     auxixiliary variable include 
 
       constraintValue(t,y,yaux,val)
       evaluateDminusConstraintDt(t,y,yaux,JCt)
       appendMinusMassMatrixDConstraintDy(t,y,gammInv,A)
       applyDConstraintDy(t,y,x,Mx)

  **********************************************************************/
public:
  RosenbrockDaeDefinition(const int& neqIn, DataCollector* data);
  virtual ~RosenbrockDaeDefinition();

  //---------- Basic Interface ----------
  virtual bool ok();

  virtual const real& getT0()     =0;     //t^0
  virtual const Vec& getY0()      =0;     //\vec y^0
  virtual const Vec& getY0prime() =0;     //\dot{\vec y}^{0}
  
  //calculate \mat{M}
  //It is very important that this not overwrite storage
  //should set jac --> jac - gamDtInv*M
  virtual bool appendMinusMassMatrix(const real& t, const Vec& y,
				     const real& gamDtInv);
  //perform operation \mat{M}.\vec{x}
  virtual bool applyMassMatrix(const real& t, const Vec& y, 
			       const Vec& x, Vec& Mx);


  //\vec F(t,\vec y)
  virtual bool rightHandSideValue(const real& t, const Vec& y, Vec& rhs) = 0;
  //\vec F_{\vec y}(t,\vec y)
  virtual bool evaluateDFDy(const real& t,const Vec& y);
  //\vec F_{t}(t,\vec y)
  virtual bool evaluateDFDt(const real& t,const Vec& y, Vec& Jt);

  //---------- for implicit ODE formulation ---------- 
  //calculate \pd{}{\vec y}[\mat{M(t,\vec y)}\vec z] for implicit
  //ode formulation and append it to the jacobian matrix
  //It is very important that this not overwrite storage
   virtual bool appendMinusDMassDyDotZ(const real& t, const Vec& y,
				       const Vec& z);
  //calculate \pd{}{t}[\mat{M(t,\vec y)}\vec z] for implicit
  //ode formulation and append it to the F_t vector
  //It is very important that this not overwrite storage
   virtual bool appendMinusDMassDtDotZ(const real& t, const Vec& y,
				       const Vec& z, Vec& Jt);


  //---------- for system with simple constraints ----------
  //\vec f_1 = \mu(y_2)-y_1
  virtual bool constraintValue(const real& t, const Vec& y,
				    const Vec& yaux, Vec& cval);

  virtual bool evaluateDConstraintDt(const real& t, const Vec& y,
				     const Vec& yaux, Vec& JCt);

  //apply \mat{f}_{1,y_2} = \pd{f_1}{\vec y_2} to x
  virtual bool applyDConstraintDy(const real& t, const Vec& y,
				  const Vec& x, Vec& Mx);

  //calculate \mat{M}_{21}.\pd{f_1}{\vec y_2} and append it to the
  //jacobian, i.e. should set jac --> jac - gamDtInv*M*F1y2
  virtual bool appendMinusMassMatrixDConstraintDy(const real& t,
						  const Vec& y,
						  const real& gamDtInv);


  virtual void stepTaken(){}
  virtual bool resetFunction();
  virtual bool jacVec(const Vec& x, Vec& Jx);

  virtual void adjustVectorTimeTolerances(const Vec& atol,
					  const Vec& rtol) {}

  virtual void resetForDiscontinuity(real tIn, const Vec& yIn) {}
  //---------- VectorFunction Interface ----------
  virtual const Vec& argument();                //\vec y^{c}
  virtual const Vec& value(bool& evalError);    //F(t^{c},\vec y^{c})
  virtual void correctArgument(Vec& correction);//y -= correction
  virtual void unCorrect();
  virtual bool evaluateAnalyticalJacobian();
  virtual bool numericalJacVec(const Vec& v, Vec& Jv);
  virtual bool analyticalJacVec(const Vec& v, Vec& Jv);
  //mwf added so that can compute own Jacobian for other integrators
  virtual void computeDeltaForJacobian();


  //---------- data members ----------
  //record function evaluations
  DataCollector* data;

  //current time
  real tDaeDef;
  //current solution
  Vec yDaeDef;
  //current value of right hand side
  Vec Fcurrent;
  //last solution and value of right hand side
  Vec yLast,Flast;

  //is Fcurrent or Jacobian out of date?
  bool updateF,updateJac,updateJacT;

  //allow object to compute its own delta for jacobian 
  //by default if integrator hasn't done it
  bool computeOwnJacobianDelta;
  //for jacobian calculations
  Vec del;
};


}//end Daetk

#endif
