#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "Definitions.h"
#include "DaeDefinition.h"

namespace Daetk
{
template <class STENCIL> 
class BoundaryCondition
{
public:
  enum face_type {CENTER, LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP,
		  INTERIOR_X,INTERIOR_Y,DUMMY_FACE}; 
  //{LEFT,RIGHT,FRONT,BACK,BOTTOM,TOP};
  enum condition_type {DIRICHLET,NEUMANN,ROBBINS,INTERIOR,DUMMY,GRADIENT};

  virtual ~BoundaryCondition(){}

//    virtual void setBoundaryConditions(face_type ft, 
//                               condition_type ct, 
//                               STENCIL& node, 
//                               real bc_value, Vec* var=0, Vec* varprime=0, real* scale=0)=0;

  virtual void setBoundaryConditions(face_type ft, 
				     condition_type ct, 
				     STENCIL& node, 
				     real bc_value, Vec* var=0, 
				     Vec* varprime=0, 
                                     real* value1=0, real* value2=0)=0;
  //mwf add if need to set a boundary condition that has a mixed type
  virtual void setBoundaryConditions(face_type ft, 
				     condition_type ct, 
				     STENCIL& node, 
				     real bc_value, bool mixedDirBC, 
				     Vec* var=0, Vec* varprime=0, 
				     real *value1=0, real *value2=0)
  { setBoundaryConditions(ft,ct,node,bc_value,var,varprime, 
			  value1,value2); }
 
  virtual void applyDirichletConditionsRHS(STENCIL&node, Vec& rhs)=0;
  virtual void print()=0;
  virtual void clear()=0;

  struct Dirichlet
  {
    face_type face;
    int n;
    real value;
    //mwf add for evaluating psk on boundary
    real *scale;
    Vec *var,*varprime;
    //mwf added for doing boundary conditions for non solution variables
    bool mixedBC;
  };
  
  struct Neumann
  {
    face_type face;
    int n;
    real value;
  };
  
  struct InteriorBC
  {
    face_type face;
    int n;
  };

  struct DummyBC
  {
    face_type face;
    int n;
    real value;
    Vec *var,*varprime;
  };

  struct Robbins
  {
    face_type face;
    int n;
    real value;
    //mwf add for evaluating psk on boundary
    //mwf switched back to reals
    real dirMult,neuMult;
    Vec *var,*varprime;
  };

  struct Gradient
  {
    face_type face;
    int n;
    real value;
  };
 
  std::vector<Dirichlet> dirichlet;
  typename std::vector<Dirichlet>::iterator dit;
  
  std::vector<Neumann> neumann;
  typename std::vector<Neumann>::iterator nit;

  std::vector<InteriorBC> interior;
  typename std::vector<InteriorBC>::iterator init;

  std::vector<DummyBC> dummy;
  typename std::vector<DummyBC>::iterator duit;

  std::vector<Robbins> robbins;
  typename std::vector<Robbins>::iterator ruit;

  std::vector<Gradient> noflow;
  typename std::vector<Gradient>::iterator noit;

  void reserveDirichlet(int n){dirichlet.reserve(n);}
  void reserveNeumann(int n){neumann.reserve(n);}
  void reserveRobbins(int n){robbins.reserve(n);}
  void reserveInterior(int n){interior.reserve(n);}
  void reserveDummy(int n){dummy.reserve(n);}
  void reserveGradient(int n){noflow.reserve(n);}
  virtual void deleteDirichlet(typename std::vector<Dirichlet>::iterator& dit_in)=0;
};
}//Daetk
#endif
