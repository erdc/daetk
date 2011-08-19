#ifndef CONSTBC_H
#define CONSTBC_H
#include <vector>
#include <map>
#include <iostream>
#include "Definitions.h"
#include "Vec.h"
#include "BandColMat.h"
#include "PetscSecondOrderFd.h"
#include "BoundaryCondition.h"

//  NEED TO ADD NEUMANN DERIVATIVES--ZERO FLUX DERIVATIVES FOR DIVKGRAD
namespace Daetk 
{

class ConstBC : public BoundaryCondition<Petsc::SecondOrderFd>
{
public:
  ConstBC();

  void applyDirichletConditions(Petsc::SecondOrderFd& node, Vec& res);
  void applyDirichletYprime(Petsc::SecondOrderFd& node);
  void applyDirichletConditionsRHS(Petsc::SecondOrderFd& node, Vec& rhs);

  template <class JAC>
  void applyDirichletDerivatives(Petsc::SecondOrderFd& node, JAC& DresDy,
				 int dof=1,int var=0);
  //mwf added for mixed bc's  
  template <class JAC>
  void adjustMixedDirichletDerivatives(Petsc::SecondOrderFd& node, 
				       JAC& DresDy,int dof,int var,
				       const Vec& dbcVarDsolVar,
				       int solvar);


  void applyNeumannConditions(Petsc::SecondOrderFd& node, Vec& flux_x);
  void applyNeumannConditions(Petsc::SecondOrderFd& node, Vec& flux_x, Vec& flux_y);
  void applyNeumannConditions(Petsc::SecondOrderFd& node, Vec& flux_x, Vec& flux_y, Vec& flux_z);

  typedef std::vector<Vec> VecVec;
  
  void applyNeumannDerivatives(Petsc::SecondOrderFd& node, VecVec& jac, 
                                      Vec& Dflux_x_center, Vec& Dflux_x_right, 
                                      real oneOverdx);

  void applyNeumannDerivatives(Petsc::SecondOrderFd& node, VecVec& jac,
                                      Vec& Dflux_x_center, Vec& Dflux_x_right, 
                                      Vec& Dflux_y_center, Vec& Dflux_y_back,
                                      real oneOverdx, real oneOverdy);

  void applyNeumannDerivatives(Petsc::SecondOrderFd& node, VecVec& jac,
                                      Vec& Dflux_x_center, Vec& Dflux_x_right, 
                                      Vec& Dflux_y_center, Vec& Dflux_y_back, 
                                      Vec& Dflux_z_center, Vec& Dflux_z_top,
                                      real oneOverdx, real oneOverdy, real oneOverdz);

  template <class JAC, class DKG>
  void applyNeumannDer(Petsc::SecondOrderFd& node, JAC& jac, DKG& dkg, const Vec& diag);


  //mwf had tried to put in
//    void setBoundaryConditions(face_type ft, 
//                               condition_type ct, 
//                               Petsc::SecondOrderFd& node, 
//                               real bc_value, Vec* var=0, Vec* varprime=0, 
//  			     real* scale=0);

  void setBoundaryConditions(face_type ft, 
                             condition_type ct, 
                             Petsc::SecondOrderFd& node, 
                             real bc_value, Vec* var=0, Vec* varprime=0, 
			     real *value1=0, real *value2=0);
  void setBoundaryConditions(face_type ft, 
                             condition_type ct, 
                             Petsc::SecondOrderFd& node, 
                             real bc_value, bool mixedDirBC, 
			     Vec* var=0, Vec* varprime=0, 
			     real *value1=0, real *value2=0);
  
  void print();
  void reserveDirichlet(int n){dirichlet.reserve(n);}
  void reserveNeumann(int n){neumann.reserve(n);}
  void clear(){dirichlet.clear(); neumann.clear();}
  void deleteDirichlet(std::vector<Dirichlet>::iterator& dit_in)
  {
    dirichlet.erase(dit_in);
  }

};

template <class JAC>
void ConstBC::applyDirichletDerivatives(Petsc::SecondOrderFd& node, JAC& DresDy,int dof,int var)
{
  //remember this just sets the scaling of the bc variable, not the solution variable
  //the actual jacobian entry should be scaling*Dbc_var/Dsol_var
  //jacobian is set for  neumann
  //DresDy is a row of the matrix with possible non-unit stride
  //Dres_depDy are the other non-unit stride row vectors
  dit = dirichlet.begin();
  while(dit != dirichlet.end()) 
    { 
      node.globalNode(dit->n);
      Petsc::SecondOrderFd::iterator sit = node.begin();
      while (sit != node.end())
        {
          for (int vj=0;vj<dof;vj++)
            DresDy(dof*node.center+var,dof*sit->globalNodeNumber+vj) = 0.0;
          ++sit;
        }
      if (dit->scale)
        DresDy(dof*node.center+var,dof*node.center+var) = (*dit->scale);
      else
        DresDy(dof*node.center+var,dof*node.center+var) = 1.0;
      DresDy.finalizeRow(dof*node.center+var);
       ++dit;
    }
}

template <class JAC>
void ConstBC::adjustMixedDirichletDerivatives(Petsc::SecondOrderFd& node, 
					      JAC& DresDy,int dof,int var,
					      const Vec& dbcVarDsolVar,
					      int solvar)
{
  //mwf now try to adjust Dirichlet boundary conditions to account for
  //the boundary condition variable's dependence on the solution variable

  dit = dirichlet.begin();
  while(dit != dirichlet.end())
    { 
      if (dit->mixedBC)
	{
	  node.globalNode(dit->n);
	  if (node.isLocal)
	    {
	      //tricky part is getting local indeces for derivative
	      //values
	      int globRow = dof*node.center+var;
	      int globCol = dof*node.center+solvar;
	      //if node.isLocal then i >= node.local_z0 etc
	      node.localIndex(node.anchor->i-node.local_z0,
			      node.anchor->j-node.local_y0,
			      node.anchor->k-node.local_x0);
	      int nlg = node.center;
	      if (dit->scale)
		DresDy(globRow,globCol) = (*dit->scale)*dbcVarDsolVar[nlg];
	      else
		DresDy(globRow,globCol) = dbcVarDsolVar[nlg];
	      DresDy.finalizeRow(globRow);
	    }//end local
	}//end mixed bc
      ++dit;
    }
}

template <class JAC, class DKG>
void ConstBC::applyNeumannDer(Petsc::SecondOrderFd& node, JAC& jac, DKG& dkg, const Vec& diag)
{
  int vi=0,vj=0;
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
//        switch(nit->face)
//          {
//          case LEFT:
//            {
//              jac(node.center,node.right) = -dkg.getDivJacRight(node,vi,vj) + dkg.oneOverdx * dkg.Dflux_x_right[vi][vj][node.interRight];
//              jac(node.center,node.center) = diag(node.center) - dkg.getDivJacCenter(node,vi,vj)+ dkg.oneOverdx * dkg.Dflux_x_center[vi][vj][node.interRight];
//              break;
//            }
//          case RIGHT:
//            {
//              jac(node.center,node.left) = -dkg.getDivJacLeft(node,vi,vj) - dkg.oneOverdx*dkg.Dflux_x_center[vi][vj][node.interLeft];
//              jac(node.center,node.center) = diag(node.center) - dkg.getDivJacCenter(node,vi,vj) - dkg.oneOverdx*dkg.Dflux_x_right[vi][vj][node.interLeft];
//              break;
//            }
//          case FRONT:
//            {
//              jac(node.center,node.back)=-dkg.getDivJacBack(node,vi,vj)+dkg.oneOverdy*dkg.Dflux_y_back[vi][vj][node.interBack];
//              jac(node.center,node.center)=diag(node.center) - dkg.getDivJacCenter(node,vi,vj)+dkg.oneOverdy*dkg.Dflux_y_center[vi][vj][node.interBack];
//              break;
//            }
//          case BACK:
//            {
//              jac(node.center,node.front)=-dkg.getDivJacFront(node,vi,vj)-dkg.oneOverdy*dkg.Dflux_y_center[vi][vj][node.interFront];
//              jac(node.center,node.center)=diag(node.center) - dkg.getDivJacCenter(node,vi,vj)-dkg.oneOverdy*dkg.Dflux_y_back[vi][vj][node.interFront];
//              break;
//            }
//          case BOTTOM:
//            {
//              jac(node.center,node.top)=-dkg.getDivJacTop(node,vi,vj)+dkg.oneOverdz*dkg.Dflux_z_top[vi][vj][node.interTop];
//              jac(node.center,node.center)=diag(node.center) - dkg.getDivJacCenter(node,vi,vj)+dkg.oneOverdz*dkg.Dflux_z_center[vi][vj][node.interTop];
//              break;
//            }
//          case TOP:
//            {
//              jac(node.center,node.bottom)=-dkg.getDivJacBottom(node,vi,vj)-dkg.oneOverdz*dkg.Dflux_z_center[vi][vj][node.interBottom];
//              jac(node.center,node.center)=diag(node.center) - dkg.getDivJacCenter(node,vi,vj)-dkg.oneOverdz*dkg.Dflux_z_top[vi][vj][node.interBottom];
//              break;
//            }
//          }
//        jac.finalizeRow(node.center);
      switch(nit->face)
        {
        case LEFT:
          {
            jac(node.center,node.right) +=dkg.oneOverdx * dkg.Dflux_x_right[vi][vj][node.interRight];
            jac(node.center,node.center) +=dkg.oneOverdx * dkg.Dflux_x_center[vi][vj][node.interRight];
            break;
          }
        case RIGHT:
          {
            jac(node.center,node.left) -= dkg.oneOverdx*dkg.Dflux_x_center[vi][vj][node.interLeft];
            jac(node.center,node.center) -= dkg.oneOverdx*dkg.Dflux_x_right[vi][vj][node.interLeft];
            break;
          }
        case FRONT:
          {
            jac(node.center,node.back)+=dkg.oneOverdy*dkg.Dflux_y_back[vi][vj][node.interBack];
            jac(node.center,node.center)+=dkg.oneOverdy*dkg.Dflux_y_center[vi][vj][node.interBack];
            break;
          }
        case BACK:
          {
            jac(node.center,node.front)-=dkg.oneOverdy*dkg.Dflux_y_center[vi][vj][node.interFront];
            jac(node.center,node.center)-=dkg.oneOverdy*dkg.Dflux_y_back[vi][vj][node.interFront];
            break;
          }
        case BOTTOM:
          {
            jac(node.center,node.top)+=dkg.oneOverdz*dkg.Dflux_z_top[vi][vj][node.interTop];
            jac(node.center,node.center)+=dkg.oneOverdz*dkg.Dflux_z_center[vi][vj][node.interTop];
            break;
          }
        case TOP:
          {
            jac(node.center,node.bottom)-=dkg.oneOverdz*dkg.Dflux_z_center[vi][vj][node.interBottom];
            jac(node.center,node.center)-=dkg.oneOverdz*dkg.Dflux_z_top[vi][vj][node.interBottom];
            break;
          }
        default:
          {
            std::cerr<<"ConstBC.h:Face Type "<<nit->face<<" not implemented"<<std::endl;
            exit (1);
          }
        }
      jac.finalizeAddRow(node.center);
      ++nit;
    }
}

}//Daetk
#endif
