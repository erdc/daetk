#include "ConstBC.h"

namespace Daetk 
{

using std::cerr;
using std::cout;
using std::endl;

ConstBC::ConstBC()
{}

void ConstBC::setBoundaryConditions(face_type ft, 
                                    condition_type ct, 
                                    Petsc::SecondOrderFd& node,
                                    real bc_value, 
                                    Vec* var, 
                                    Vec* varprime,
                                    real* value1,
                                    real* value2) 
{
  if (node.isLocal)
    {
      switch (ct)
        {
        case DIRICHLET:
          {
            assert(var != 0);
            assert(varprime != 0);
            Dirichlet d={ft,node.petscToGlobal(node.center),bc_value,value1,var,varprime,false};
            dirichlet.push_back(d);
            break;
          }
        case NEUMANN:
          {
            Neumann n={ft,node.petscToGlobal(node.center),bc_value};
            neumann.push_back(n);
            break;
          }
        default:
          {
            cerr<<"Boundary condition type "<<ft<<"not implemented"<<endl;
            exit(1);
          }
        }
    }
}

void ConstBC::setBoundaryConditions(face_type ft, 
                                    condition_type ct, 
                                    Petsc::SecondOrderFd& node,
                                    real bc_value, 
				    bool mixedDirBC,
                                    Vec* var, 
                                    Vec* varprime,
                                    real* value1,
                                    real* value2) 
{
  if (mixedDirBC && ct == DIRICHLET)
    {
      if (node.isLocal)
	{
	  assert(var != 0);
	  assert(varprime != 0);
	  Dirichlet d={ft,node.petscToGlobal(node.center),
		       bc_value,value1,var,varprime,true};
	  dirichlet.push_back(d);
	}
    }
  else
    {

      setBoundaryConditions(ft,ct,node,bc_value,var,varprime,
			    value1,value2);
    }
}

void ConstBC::print()
{
  dit=dirichlet.begin();
  while(dit!=dirichlet.end())
    {
      cout<<"D"<<'\t'<<dit->face<<'\t'<<dit->n<<'\t'<<dit->value<<endl;
      ++dit;
    }
  
  nit=neumann.begin();
  while(nit!=neumann.end())
    {
      cout<<"N"<<'\t'<<nit->face<<'\t'<<nit->n<<'\t'<<nit->value<<endl;
      ++nit;
    }
}

void ConstBC::applyNeumannConditions(Petsc::SecondOrderFd& node, Vec& flux_x, Vec& flux_y)
{
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      switch(nit->face)
        {
        case LEFT:
          {
            flux_x[node.interLeft] = - 2.0*nit->value - flux_x[node.interRight];
            break;
          }
        case RIGHT:
          {
            flux_x[node.interRight] = 2.0*nit->value - flux_x[node.interLeft];
            break;
          }
        case FRONT:
          {
            flux_y[node.interFront] = - 2.0*nit->value - flux_y[node.interBack];
            break;
          }
        case BACK:
          {
            flux_y[node.interBack] = 2.0*nit->value - flux_y[node.interFront];
            break;
          }
        default:
          {
            std::cerr<<"ConstBC.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
            exit (1);
          }
        }
      ++nit;
    }
}

void ConstBC::applyNeumannConditions(Petsc::SecondOrderFd& node, Vec& flux_x, Vec& flux_y, Vec& flux_z)
{
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      switch(nit->face)
        {
        case LEFT:
          {
            flux_x[node.interLeft] = - 2.0*nit->value  - flux_x[node.interRight];
            break;
          }
        case RIGHT:
          {
            flux_x[node.interRight] = 2.0*nit->value  - flux_x[node.interLeft];
            break;
          }
        case FRONT:
          {
            flux_y[node.interFront] = - 2.0*nit->value  - flux_y[node.interBack];
            break;
          }
        case BACK:
          {
            flux_y[node.interBack] = 2.0*nit->value - flux_y[node.interFront];
            break;
          }
        case BOTTOM:
          {
            flux_z[node.interBottom] = - 2.0*nit->value - flux_z[node.interTop];
            break;
          }
        case TOP:
          {
            flux_z[node.interTop] = 2.0*nit->value - flux_z[node.interBottom];
            break;
          }
//          case LEFT:
//            {
//              flux_x[node.interLeft] = - nit->value;
//              break;
//            }
//          case RIGHT:
//            {
//              flux_x[node.interRight] = nit->value;
//              break;
//            }
//          case FRONT:
//            {
//              flux_y[node.interFront] = - nit->value;
//              break;
//            }
//          case BACK:
//            {
//              flux_y[node.interBack] = nit->value;
//              break;
//            }
//          case BOTTOM:
//            {
//              flux_z[node.interBottom] = - nit->value;
//              break;
//            }
//          case TOP:
//            {
//              flux_z[node.interTop] = nit->value;
//              break;
//            }
        default:
          {
            std::cerr<<"ConstBC.h:Face Type "<<nit->face<<" not implemented"<<std::endl;
            exit (1);
          }
        }
      ++nit;
    }
}

void ConstBC::applyNeumannConditions(Petsc::SecondOrderFd& node, Vec& flux_x)
{
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      switch(nit->face)
        {
        case LEFT:
          {
            flux_x[node.interLeft] = - 2.0*nit->value - flux_x[node.interRight];
            break;
          }
        case RIGHT:
          {
            flux_x[node.interRight] = 2.0*nit->value - flux_x[node.interLeft];
            break;
          }
        default:
          {
            std::cerr<<"ConstBC.h:Face Type "<<nit->face<<" not implemented"<<std::endl;
            exit (1);
          }
        }
      ++nit;
    }
}

void ConstBC::applyNeumannDerivatives(Petsc::SecondOrderFd& node, VecVec& jac,
                                             Vec& Dflux_x_center, Vec& Dflux_x_right, 
                                             Vec& Dflux_y_center, Vec& Dflux_y_back,
                                             real oneOverdx, real oneOverdy)
{
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      switch(nit->face)
        {
        case LEFT:
          {
            jac[LEFT](node.center)=0.0;
            jac[RIGHT](node.center)-=oneOverdx*Dflux_x_right[node.interRight];
            jac[CENTER](node.center)-=oneOverdx*Dflux_x_center[node.interRight];
            break;
          }
        case RIGHT:
          {
            jac[LEFT](node.center)+=oneOverdx*Dflux_x_center[node.interLeft];
            jac[RIGHT](node.center)=0.0;
            jac[CENTER](node.center)+=oneOverdx*Dflux_x_right[node.interLeft];
            break;
          }
        case FRONT:
          {
            jac[FRONT](node.center)=0.0;
            jac[BACK](node.center)-=oneOverdy*Dflux_y_back[node.interBack];
            jac[CENTER](node.center)-=oneOverdy*Dflux_y_center[node.interBack];
            break;
          }
        case BACK:
          {
            jac[FRONT](node.center)+=oneOverdy*Dflux_y_center[node.interFront];
            jac[BACK](node.center)=0.0;
            jac[CENTER](node.center)+=oneOverdy*Dflux_y_back[node.interFront];
            break;
          }
        default:
          {
            std::cerr<<"ConstBC.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
            exit (1);
          }
        }
      ++nit;
    }
}
  

void ConstBC::applyNeumannDerivatives(Petsc::SecondOrderFd& node, VecVec& jac,
                                             Vec& Dflux_x_center, Vec& Dflux_x_right, 
                                             Vec& Dflux_y_center, Vec& Dflux_y_back, 
                                             Vec& Dflux_z_center, Vec& Dflux_z_top,
                                             real oneOverdx, real oneOverdy, real oneOverdz)
{
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      switch(nit->face)
        {
        case LEFT:
          {
            jac[LEFT](node.center)=0.0;
            jac[RIGHT](node.center)-=oneOverdx*Dflux_x_right[node.interRight];
            jac[CENTER](node.center)-=oneOverdx*Dflux_x_center[node.interRight];
            break;
          }
        case RIGHT:
          {
            jac[LEFT](node.center)+=oneOverdx*Dflux_x_center[node.interLeft];
            jac[RIGHT](node.center)=0.0;
            jac[CENTER](node.center)+=oneOverdx*Dflux_x_right[node.interLeft];
            break;
          }
        case FRONT:
          {
            jac[FRONT](node.center)=0.0;
            jac[BACK](node.center)-=oneOverdy*Dflux_y_back[node.interBack];
            jac[CENTER](node.center)-=oneOverdy*Dflux_y_center[node.interBack];
            break;
          }
        case BACK:
          {
            jac[FRONT](node.center)+=oneOverdy*Dflux_y_center[node.interFront];
            jac[BACK](node.center)=0.0;
            jac[CENTER](node.center)+=oneOverdy*Dflux_y_back[node.interFront];
            break;
          }
        case BOTTOM:
          {
            jac[BOTTOM](node.center)=0.0;
            jac[TOP](node.center)-=oneOverdz*Dflux_z_top[node.interTop];
            jac[CENTER](node.center)-=oneOverdz*Dflux_z_center[node.interTop];
            break;
          }
        case TOP:
          {
            jac[BOTTOM](node.center)+=oneOverdz*Dflux_z_center[node.interBottom];
            jac[TOP](node.center)=0.0;
            jac[CENTER](node.center)+=oneOverdz*Dflux_z_top[node.interBottom];
            break;
          }
        default:
          {
            std::cerr<<"ConstBC.h:Face Type "<<nit->face<<" not implemented"<<std::endl;
            exit (1);
          }
        }
      ++nit;
    }
}

void ConstBC::applyNeumannDerivatives(Petsc::SecondOrderFd& node, VecVec& jac, 
                                             Vec& Dflux_x_center, Vec& Dflux_x_right,
                                             real oneOverdx)
{
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      switch(nit->face)
        {
        case LEFT:
          {
            jac[LEFT](node.center)=0.0;
            jac[RIGHT](node.center)-=oneOverdx*Dflux_x_right[node.interRight];
            jac[CENTER](node.center)-=oneOverdx*Dflux_x_center[node.interRight];
            break;
          }
        case RIGHT:
          {
            jac[LEFT](node.center)+=oneOverdx*Dflux_x_center[node.interLeft];
            jac[RIGHT](node.center)=0.0;
            jac[CENTER](node.center)+=oneOverdx*Dflux_x_right[node.interLeft];
            break;
          }
        default:
          {
            std::cerr<<"ConstBC.h:Face Type "<<nit->face<<" not implemented"<<std::endl;
            exit (1);
          }
        }
      ++nit;
    }
}

void ConstBC::applyDirichletConditions(Petsc::SecondOrderFd& node, Vec& res)
{
  dit = dirichlet.begin();
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      if (dit->scale)
        res(node.center) = ((*dit->var)(node.center) - dit->value)*(*dit->scale);
      else
        res(node.center) = ((*dit->var)(node.center) - dit->value);
      ++dit;
    }
}

void ConstBC::applyDirichletConditionsRHS(Petsc::SecondOrderFd& node, Vec& rhs)
{
  dit = dirichlet.begin();
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      if (dit->scale)
        rhs(node.center) = dit->value*(*dit->scale);
      else
        rhs(node.center) = dit->value;
      ++dit;
    }
}

void ConstBC::applyDirichletYprime(Petsc::SecondOrderFd& node)
{
  //yprime is set to the yprime for neumann 
  dit = dirichlet.begin();
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      (*dit->varprime)(node.center) = 0.0;
      ++dit;
    }
}


}//Daetk
