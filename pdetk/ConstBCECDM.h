#ifndef ConstBCECDM_H
#define ConstBCECDM_H
#include <vector>
#include <map>
#include <iostream>
#include "Definitions.h"
#include "Vec.h"
#include "BandColMat.h"
#include "BoundaryCondition.h"
#include "PetscStencilMM.h"

//mwf can I include DivKgradECDM class into this definition
#include "DivKgradECDM.h"
//#define USE_CONT_INTERFACE_TERMS_ECDM

namespace Daetk 
{
using std::cerr;
using std::endl;

//  NEED TO ADD NEUMANN DERIVATIVES--ZERO FLUX DERIVATIVES FOR DIVKGRAD

class ConstBCECDM : public BoundaryCondition<Petsc::StencilMM>
{
public:
  void applyDirichletConditionsRHS(Petsc::StencilMM&, Daetk::Petsc::Vec&)
  {std::cerr<<"not yet implemented for ECDM"<<std::endl;}

  inline void setCurrentDivKgrad(DivKgradECDM<ConstBCECDM,1>*  div)
  {
    assert(div);
    currentDiv = div;
  }

  ConstBCECDM():currentDiv(0)
    {}

  //==================== interface terms====================

  template <int nv>
  inline 
  void applyInteriorConditionToConstRelation(Petsc::StencilMM& 
					     node,
					     DivKgradECDM<ConstBCECDM,nv>*  div,
					     Vec& Kr,Vec& Rho,
					     const Vec& pVal);
   
  void applyInteriorConditionToConstRelationDeriv(Petsc::StencilMM& 
						  node,
						  Vec& DKr,Vec& DRho);

  //=== 1d ===
  void applyInteriorConditionToInterfaceValues(Petsc::StencilMM& 
					       node,
					       const Vec* Kr,
					       const Vec* Rho,
					       Vec& KrX,Vec& RhoX);
  //=== 2d ===
  void applyInteriorConditionToInterfaceValues(Petsc::StencilMM& 
					       node,
					       const Vec* Kr,
					       const Vec* Rho,
					       Vec& KrX,Vec& KrY,
					       Vec& RhoX,Vec& RhoY);

  template <int nv>
  inline 
  void applyInteriorConditionToInterfaceValueDerivs(Petsc::StencilMM& 
						    node,
						    DivKgradECDM<ConstBCECDM,nv>*
						    div,
						    const VecVecVec& DKr,
						    const Vec* DRho);


  //ECDM has explicit dirichlet neumann and dummy conditions 
  //these are for setting the residual and jacobian directly
  //==================== dirichlet ====================
  inline void applyDirichletConditionsToRes(Petsc::StencilMM& node, 
					    Vec& res, real scaling=1.0);
  inline void applyDirichletYprime(Petsc::StencilMM& node);

  template <class JAC>
  inline void applyDirichletDerivatives(Petsc::StencilMM& node, JAC& DresDy,
					real scaling=1.0,int dof=1,int var=0);  

  //==================== robbins ====================

  //make global residual be flux-bc
  template <int nv>
  inline void applyRobbinsConditionsToRes(Petsc::StencilMM& node, Vec& res,
					  DivKgradECDM<ConstBCECDM,nv>* div,
					  real scaling=1.0);

  inline void applyRobbinsConditionsYprime(Petsc::StencilMM& node, 
					   Vec& yPrime);  
  template <class JAC>
  inline void applyRobbinsDerivativesToRes(Petsc::StencilMM& node, JAC& DresDy,
					   VecVec& Div_D, int DIM, 
					   real scaling =1.0);


   //==================== dummy ====================
  inline void applyDummyConditionsToRes(Petsc::StencilMM& node, 
					Vec& res, real scaling=1.0);
  inline void applyDummyYprime(Petsc::StencilMM& node);

  template <class JAC>
  inline void applyDummyDerivatives(Petsc::StencilMM& node, JAC& DresDy,
				    real scaling=1.0,int dof=1,int var=0);  


  //==================== neumann ====================

  //make global residual be flux-bc
  template <int nv>
  inline void applyNeumannConditionsToRes(Petsc::StencilMM& node, Vec& res,
					  DivKgradECDM<ConstBCECDM,nv>* div,
					  real scaling=1.0);

  inline void applyNeumannConditionsYprime(Petsc::StencilMM& node, 
					   Vec& yPrime);  
  template <class JAC>
  inline void applyNeumannDerivativesToRes(Petsc::StencilMM& node, JAC& DresDy,
					   VecVec& Div_D, int DIM, 
					   real scaling =1.0);
  //====== these use mass balance and neumann bc
  template <int nv>
  inline void applyNeumannConditionsYprime(Petsc::StencilMM& node, 
					   DivKgradECDM<ConstBCECDM,nv>* div,
					   Vec& yPrime);  
  //=== pressure head version
  //mwf now try and put in approx. using mass bal
  template <int nv>
  inline void applyNeumannConditionsToRes(Petsc::StencilMM& node, Vec& res,
					  DivKgradECDM<ConstBCECDM,nv>* div,
					  const Vec& rho, const Vec& theta,
					  const Vec& Drho, const Vec& Dtheta,
					  const Vec& detMap, const Vec& Dp,
					  real scaling=1.0);


  //mwf now try and put in approx. using mass bal
  template <class JAC, int nv>
  inline void applyNeumannDerivativesToRes(Petsc::StencilMM& node, JAC& DresDy,
					   VecVec& Div_D, int DIM, 
					   DivKgradECDM<ConstBCECDM,nv>* div,
					   const Vec& rho, const Vec& theta,
					   const Vec& Drho, const Vec& Dtheta,
					   const Vec& detMap, const Vec& Dp, 
					   real alphaBDF, const Vec& DDrho, 
					   const Vec& DDtheta,
					   const Vec& hXcen, const Vec& hYcen,
					   real scaling =1.0);
  //=== mass conservative version
  //mwf now try and put in approx. using mass bal
  template <int nv>
  inline void applyNeumannConditionsToRes(Petsc::StencilMM& node, Vec& res,
					  DivKgradECDM<ConstBCECDM,nv>* div,
					  const Vec& mCurrent, 
					  const Vec& betaDaeDef,
					  const Vec& detMap,
					  real alphaBDF,
					  real scaling=1.0);


  //mwf now try and put in approx. using mass bal
  template <class JAC, int nv>
  inline void applyNeumannDerivativesToRes(Petsc::StencilMM& node, JAC& DresDy,
					   VecVec& Div_D, int DIM, 
					   DivKgradECDM<ConstBCECDM,nv>* div,
					   const Vec& rho, const Vec& theta,
					   const Vec& Drho, const Vec& Dtheta,
					   const Vec& detMap, 
					   real alphaBDF, 
					   real scaling =1.0);

  //==================== interior ====================
  inline void applyInteriorConditionsYprime(Petsc::StencilMM& node, 
					    Vec& yPrime);  
//make global residual be fluxLeft-fluxRight
  template <int nv>
  inline void applyInteriorConditionsToRes(Petsc::StencilMM& node, Vec& res,
					   DivKgradECDM<ConstBCECDM,nv>* div,
					   real scaling=1.0);


  template <class JAC>
  inline void applyInteriorDerivativesToRes(Petsc::StencilMM& node, JAC& DresDy,
					    VecVec& Div_D, int DIM, 
					    real scaling =1.0);

 
 //this one updates derivatives straight into jacobian for left and right terms  
  template <class JAC, int nv>
  inline 
  void updateInteriorDerivativesToRes(Petsc::StencilMM& node, 
				      JAC& DresDy,
				      const VecVec* K, const Vec* Kr,
				      const Vec* Rho, const Vec* P,
				      const VecVecVec& DKr, 
				      const Vec* DRho,
				      VecVec& D_Div, int DIM,
				      DivKgradECDM<ConstBCECDM,nv>* div,
				      real scaling=1.0);
  


  //==================== no flow ====================

  //make global residual be flux-bc
  template <int nv>
  inline void applyNoFlowConditionsToRes(Petsc::StencilMM& node, Vec& res,
					 DivKgradECDM<ConstBCECDM,nv>* div,
					 real scaling=1.0);

  inline void applyNoFlowConditionsYprime(Petsc::StencilMM& node, 
					  Vec& yPrime);  
  template <class JAC, int nv>
  inline void applyNoFlowDerivativesToRes(Petsc::StencilMM& node, JAC& DresDy,
					  DivKgradECDM<ConstBCECDM,nv>* div,
					  int DIM, 
					  real scaling =1.0); 
  

   //==================== all ====================
  //account for 3 exterior boundary conditions in fluxes and
  //pressure differences
  //1d
  void 
  applyBoundaryConditionsToPressureDifferences(Petsc::StencilMM& node,
					       Vec& pDiff_x);

  void 
  applyBoundaryConditionsToPressureDifferences(Petsc::StencilMM& node,
					       Vec& pDiff_x,
					       Vec& pDiff_y);

  //=== 1d ===
  void 
  applyBoundaryConditionsToFluxes(Petsc::StencilMM& node,
				  Vec& flux_x);

  //=== 2d ===
  void 
  applyBoundaryConditionsToFluxes(Petsc::StencilMM& node,
				  Vec& flux_x,
				  Vec& flux_y);

  //=== 1d ====
  //zero everything except interior pressure difference term
  void 
  applyBoundaryConditionsToPressureDerivatives(Petsc::StencilMM& node,
					       Vec& DpDiff_x_center,
					       Vec& DpDiff_x_right);
  //=== 2d ====
  //zero everything except interior pressure difference term
  void 
  applyBoundaryConditionsToPressureDerivatives(Petsc::StencilMM& node,
					       Vec& DpDiff_x_center,
					       Vec& DpDiff_x_right,
					       Vec& DpDiff_y_center,
					       Vec& DpDiff_y_back);

  //=== 1d ====
  //zero everything except derivatives evaluated at interior point
  void applyBoundaryConditionDerivativesECDM(Petsc::StencilMM& node,
					     Vec& Dflux_x_center, 
					     Vec& Dflux_x_right);
  //=== 2d ====
  //zero everything except derivatives evaluated at interior point
  void applyBoundaryConditionDerivativesECDM(Petsc::StencilMM& node,
					     Vec& Dflux_x_center, 
					     Vec& Dflux_x_right,
					     Vec& Dflux_x_back, 
					     Vec& Dflux_x_front,
					     Vec& Dflux_x_rightBack, 
					     Vec& Dflux_x_rightFront,
					     Vec& Dflux_y_center, 
					     Vec& Dflux_y_back, 
					     Vec& Dflux_y_right, 
					     Vec& Dflux_y_left,
					     Vec& Dflux_y_backRight, 
					     Vec& Dflux_y_backLeft);
  
  //////////////////////////

  //this is for time integration
  void adjustTimeIntegrationTolerance(Petsc::StencilMM& node,
				      Vec &  tol, real newTol);

  void setBoundaryConditions(face_type ft, 
                             condition_type ct, 
                             Petsc::StencilMM& node, 
                             real bc_value, Vec* var=0, Vec* varprime=0,
			     real *value1 = 0,
			     real *value2 = 0);

  //mwf try and make sure lists don't share values?
  void checkForDuplicateBCs();

  void print();
  void clear()
  {
    dirichlet.clear(); neumann.clear(); robbins.clear();
    interior.clear(); dummy.clear(); noflow.clear();
  }

  void deleteDirichlet(std::vector<Dirichlet>::iterator& dit_in)
  {
    dirichlet.erase(dit_in);
  }

  //mwf added to simplify interface, may not really need
  DivKgradECDM<ConstBCECDM,1>*  currentDiv;
  

};
/////////////////////////////////////////
//==================== interface conditions for const relations ==========
//try to put in either arithmetic or upwind average
//into constitutive relation value at internal interface
template <int nv>
inline void 
ConstBCECDM::applyInteriorConditionToConstRelation(Petsc::StencilMM& 
						   node,
						   DivKgradECDM<ConstBCECDM,nv>*
						   div,
						   Vec& Kr,Vec& Rho,
						   const  Vec& pVal)
{
  init = interior.begin();
  //std::cout<<init<<'\t'<<interior.end()<<std::endl<<std::flush;
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {

		Rho[node.center] = div->avgDensityAcrossInteriorBoundaryX(Rho,pVal);
		Kr[node.center]  = div->avgRelPermAcrossInteriorBoundaryX(Kr,Rho,pVal);

		break; 
	      }
	    case INTERIOR_Y:
	      {
		//arithmetic average
		Rho[node.center] = div->avgDensityAcrossInteriorBoundaryY(Rho,pVal);
		Kr[node.center]  = div->avgRelPermAcrossInteriorBoundaryY(Kr,Rho,pVal);
		
		break;
	      }
	    default:
	      {
		cerr<<"wrong face = "<<init->face
		    <<" in applyInteriorConditionToConstRelation "<<endl;
		assert(0);
		break;
	      }

	    }//end switch
	}//end if local
      
      ++init;
    }//end while
}//end function


template <int nv>
inline 
void 
ConstBCECDM::applyInteriorConditionToInterfaceValueDerivs(Petsc::StencilMM& 
							  node,
							  DivKgradECDM<ConstBCECDM,nv>*
							  div,
							  const VecVecVec& DKr,
							  const Vec* DRho)
{
  init = interior.begin();
  //std::cout<<init<<'\t'<<interior.end()<<std::endl<<std::flush;
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf only set if node is local to this processor?
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
		div->setDKrXinterface(DKr);
		div->setDRhoXinterface(DRho);
		break; 
	      }
	    case INTERIOR_Y:
	      {
		div->setDKrYinterface(DKr);
		div->setDRhoYinterface(DRho);
		break;
	      }
	    default:
	      {
		cerr<<"wrong face = "<<init->face
		    <<" in applyInteriorConditionToConstRelation "<<endl;
		assert(0);
		break;
	      }

	    }//end switch
	}//end if local
      
      ++init;
    }//end while
}//end function

//////////////////////

//=========================dirichlet====================
inline void ConstBCECDM::applyDirichletConditionsToRes(Petsc::StencilMM& node, 
						       Vec& res, real scaling)
{
  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      if (node.isLocal)
	res(node.center) = ((*dit->var)(node.center) - dit->value)*scaling;
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    }
}

inline void ConstBCECDM::applyDirichletYprime(Petsc::StencilMM& node)
{
  //yprime is set to the yprime for neumann 
  dit = dirichlet.begin();
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      if (node.isLocal)
	(*dit->varprime)(node.center) = 0.0;
      ++dit;
    }
}

template <class JAC>
inline void ConstBCECDM::applyDirichletDerivatives(Petsc::StencilMM& node, 
						   JAC& DresDy,
						   real scaling,int dof,
						   int var)
{
  //remember this just sets the scaling of the bc variable, not the solution variable
  //the actual jacobian entry should be scaling*Dbc_var/Dsol_var
  //jacobian is set for  neumann
  //DresDy is a row of the matrix with possible non-unit stride
  //Dres_depDy are the other non-unit stride row vectors
  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl;
  while(dit != dirichlet.end()) 
    { 
      node.globalNode(dit->n);
      //mwf added isLocal here 

      if (node.isLocal)
	{
	  Petsc::StencilMM::iterator sit = node.begin();
	  while (sit != node.end())
	    {
	      //DresDy(node.center,sit->second.globalNodeNumber) = 0.0;
	      for (int vj=0;vj<dof;vj++)
		DresDy(dof*node.center+var,dof*sit->globalNodeNumber+vj) = 0.0;
	      ++sit;
	    }
	  //DresDy.zeroRow(node.center);
	  DresDy(dof*node.center+var,dof*node.center+var) = scaling;
	  DresDy.finalizeRow(dof*node.center+var);
	}
      ++dit;
    }
}

//====================dummy===========================

inline void ConstBCECDM::applyDummyConditionsToRes(Petsc::StencilMM& node, 
						   Vec& res, real scaling)
{
  duit = dummy.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      if (node.isLocal)
	res(node.center) = ((*duit->var)(node.center) - duit->value)*scaling;
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++duit;
    }
}

inline void ConstBCECDM::applyDummyYprime(Petsc::StencilMM& node)
{

  duit = dummy.begin();
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      if (node.isLocal)
	(*duit->varprime)(node.center) = 0.0;
      ++duit;
    }
}

template <class JAC>
inline void ConstBCECDM::applyDummyDerivatives(Petsc::StencilMM& node, 
					       JAC& DresDy,
					       real scaling,int dof,int var)
{
  //remember this just sets the scaling of the bc variable, not the solution variable
  //the actual jacobian entry should be scaling*Dbc_var/Dsol_var
  //jacobian is set for  neumann
  //DresDy is a row of the matrix with possible non-unit stride
  //Dres_depDy are the other non-unit stride row vectors
  duit = dummy.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl;
  while(duit != dummy.end()) 
    { 
      node.globalNode(duit->n);
      //mwf added isLocal here 
      if (node.isLocal)
	{
	  Petsc::StencilMM::iterator sit = node.begin();
	  while (sit != node.end())
	    {
	      //DresDy(node.center,sit->second.globalNodeNumber) = 0.0;
	      for (int vj=0;vj<dof;vj++)
		DresDy(dof*node.center+var,dof*sit->globalNodeNumber+vj) = 0.0;
	      ++sit;
	    }
	  //DresDy.zeroRow(node.center);
	  DresDy(dof*node.center+var,dof*node.center+var) = scaling;
	  DresDy.finalizeRow(dof*node.center+var);
	}
       ++duit;
    }
}


//======================neumann=========================
//===== these use qval-bc
//mwf I don't know yet what to do about neumann and interior bc's
//as far as time integration goes
inline
void  ConstBCECDM::applyNeumannConditionsYprime(Petsc::StencilMM& node, 
						Vec& yPrime)
{

  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  yPrime(globalCenter) = 0.0;
	}
      ++nit;
    }
}

//make global residual be fluxInterior-qVal
template <int nv>
inline 
void  ConstBCECDM::applyNeumannConditionsToRes(Petsc::StencilMM& node, 
					       Vec& res,
					       DivKgradECDM<ConstBCECDM,nv>* div,
					       real scaling)
{
  assert(nv == 1);
  assert(div);
  const Vec & flux_x = div->getFlux_x(0);
  const Vec & flux_y = div->getFlux_y(0);

  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(nit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		//using interior flux for the Dirichlet face
		//and zeroing the exterior flux
		//remember sign of normal for minus side 
		res(globalCenter) = (flux_x[node.interRight] + nit->value)*scaling;
		if (fabs(res(globalCenter)) > 1.0e4)
		  {
		    cerr<<"in ConstBCECDM::neumannBCsToRes, node.center= "
			<<node.center<<endl;
		    cerr<<"\t flux_x["<<node.interRight<<"]= "
			<<flux_x[node.interRight]<<endl;
		    cerr<<"\t scaling = "<<scaling<<" nit->value= "<<nit->value<<endl;
		  } 
		break;
	      }
	    case RIGHT:
	      {
		res(globalCenter) = (flux_x[node.interLeft] - nit->value)*scaling;
		if (fabs(res(globalCenter)) > 1.0e4)
		  {
		    cerr<<"in ConstBCECDM::neumannBCsToRes, node.center= "
			<<node.center<<endl;
		    cerr<<"\t flux_x["<<node.interLeft<<"]= "
			<<flux_x[node.interLeft]<<endl;
		    cerr<<"\t scaling = "<<scaling<<" nit->value= "<<nit->value<<endl;
		  } 
		break;
	      }
	    case FRONT:
	      {
		//remeber sign for negative boundary
		res(globalCenter) = (flux_y[node.interBack] +nit->value)*scaling;
		//              cout<<"FF"<<'\t'<<flux_y[node.interFront]<<'\t'<<nit->value <<'\t'<< flux_y[node.interBack]<<endl;
		break;
	      }
	    case BACK:
	      {
		res(globalCenter) = (flux_y[node.interFront] -nit->value)*scaling;
		std::cout<<"FB"<<'\t'<<flux_y[node.interFront] <<'\t'<<nit->value<<'\t'<< res(globalCenter)<<endl;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }
}

//make global residual be fluxInterior-qVal
template <class JAC>
inline 
void  ConstBCECDM::applyNeumannDerivativesToRes(Petsc::StencilMM& node, 
						JAC& DresDy,
						VecVec& D_Div, int DIM,
						real scaling)
{

  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      int nl=node.center-node.globalLow;
      real signOfDiv(1.0);
      if (node.isLocal)
	{
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		signOfDiv = -1.0;
		break;
	      }
	    case RIGHT:
	      {
		signOfDiv = 1.0;
		break;
	      }
	    case FRONT:
	      {
		signOfDiv = -1.0;
		break;
	      }
	    case BACK:
	      {
		signOfDiv = 1.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	  signOfDiv *= scaling;
	      
          DresDy(node.center,node.center)  = 
	    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::CENTER][nl];
          
          if (node.anchor->k > 0)
            DresDy(node.center,node.left)  = 
	     signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT][nl];

          if (node.anchor->k < node.nxNodes-1)
            DresDy(node.center,node.right) = 
	      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT][nl];
          
          if (DIM > 0)
            {
              if (node.anchor->j > 0)
		{
		  DresDy(node.center,node.front)= 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::FRONT][nl];
		  //(j-1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(node.center,node.leftFront)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_FRONT][nl];
		  //(j-1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(node.center,node.rightFront) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_FRONT][nl];

		}
              if (node.anchor->j < node.nyNodes-1)
                {
		  DresDy(node.center,node.back) = 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BACK][nl];
		  //(j+1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(node.center,node.leftBack)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_BACK][nl];
		  //(j+1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(node.center,node.rightBack) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_BACK][nl];
		  
		}
	    }//end DIM==TWO_D
          if (DIM > 1) //DIM==THREE_D
            {
              if (node.anchor->i > 0)
                DresDy(node.center,node.bottom)= 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BOTTOM][nl];
              if (node.anchor->i < node.nzNodes-1)
                DresDy(node.center,node.top) = 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::TOP][nl];
            }
	  DresDy.finalizeRow(node.center);
        }//end isLocal
      ++nit;
    }//end while
}//end function

//=============== these use mass balance
template <int nv>
inline
void  ConstBCECDM::applyNeumannConditionsYprime(Petsc::StencilMM& node, 
						DivKgradECDM<ConstBCECDM,nv>* 
						div,
						Vec& yPrime)
{
  assert(nv == 1);
  assert(div);

  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		int globalRight = node.right;

		yPrime(globalCenter) = yPrime(globalRight);

		break;
	      }
	    case RIGHT:
	      {
		int globalLeft = node.left;


		yPrime(globalCenter) = yPrime(globalLeft);

		break;
	      }
	    case FRONT:
	      {

		int globalBack = node.back;
		yPrime(globalCenter) = yPrime(globalBack);

		break;
	      }
	    case BACK:
	      {
		int globalFront = node.front;
		yPrime(globalCenter) = yPrime(globalFront);
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }
}

//========= pressure head form 
//mwf now try and put in approx. using mass bal. only 2d for now
template <int nv>
inline void 
ConstBCECDM::applyNeumannConditionsToRes(Petsc::StencilMM& node, Vec& res,
					 DivKgradECDM<ConstBCECDM,nv>* div,
					 const Vec& rho, const Vec& theta,
					 const Vec& Drho, const Vec& Dtheta,
					 const Vec& detMap, const Vec& Dp,
					 real scaling)
{
  assert(nv == 1);
  assert(div);
  const Vec & flux_x = div->getFlux_x(0);
  const Vec & flux_y = div->getFlux_y(0);
  const Vec & hXcen  = div->getSpatialIncrement(0);
  const Vec & hYcen  = div->getSpatialIncrement(1);

  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      int globalCenter= node.center;
      int globalLeft  = node.left;
      int globalRight = node.right;
      int globalFront = node.front;
      int globalBack  = node.back;

      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(nit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);

	  switch(nit->face)
	    {
	    case LEFT:
	      {
		//remember sign of normal for minus side 
		//res(globalCenter) = (flux_x[node.interRight] + nit->value)*scaling;

		int globalNeighbor=globalRight;
		int localNeighbor = node.right;

		real accumulation=(rho(globalNeighbor)*Dtheta(globalNeighbor) 
				   + Drho(globalNeighbor)*theta(globalNeighbor));


		//use Dp as global variable
		accumulation *= Dp(globalNeighbor);

		real volume = detMap(globalNeighbor)*hYcen[localNeighbor]
		  *hXcen[localNeighbor];

		real fluxCalc = flux_x[node.interLeft];


		//now switch to localNeighbor to get divergence
		node.localNode(localNeighbor);

		//normal res calc uses -div from DivKgradECDM
		res(globalCenter) = volume*accumulation
		  + flux_x[node.interRight] - fluxCalc  
		  + flux_y[node.interBack] - flux_y[node.interFront];


		if (fabs(res(globalCenter)) > 1.0e4)
		  {
		    cerr<<"in ConstBCECDM::neumannBCsToRes, node.center= "
			<<node.center<<endl;
		    cerr<<"\t flux_x["<<node.interRight<<"]= "
			<<flux_x[node.interRight]<<endl;
		    cerr<<"\t scaling = "<<scaling<<" nit->value= "<<nit->value<<endl;
		  } 
		break;
	      }
	    case RIGHT:
	      {
		//  		res(globalCenter) = (flux_x[node.interLeft] - nit->value)*scaling;
		int globalNeighbor=globalLeft;
		int localNeighbor = node.left;

		real accumulation=(rho(globalNeighbor)*Dtheta(globalNeighbor) 
				   + Drho(globalNeighbor)*theta(globalNeighbor));


		//use Dp as global variable
		accumulation *= Dp(globalNeighbor);

		real volume = detMap(globalNeighbor)*hYcen[localNeighbor]
		  *hXcen[localNeighbor];

		real fluxCalc = flux_x[node.interRight];

		//now switch to localNeighbor to get divergence
		node.localNode(localNeighbor);

		//normal res calc uses -div from DivKgradECDM
		res(globalCenter) = volume*accumulation
		  + fluxCalc - flux_x[node.interLeft] 
		  + flux_y[node.interBack] - flux_y[node.interFront];

		if (fabs(res(globalCenter)) > 1.0e4)
		  {
		    cerr<<"in ConstBCECDM::neumannBCsToRes, node.center= "
			<<node.center<<endl;
		    cerr<<"\t flux_x["<<node.interLeft<<"]= "
			<<flux_x[node.interLeft]<<endl;
		    cerr<<"\t scaling = "<<scaling<<" nit->value= "<<nit->value<<endl;
		  } 
		break;
	      }
	    case FRONT:
	      {
		//remeber sign for negative boundary
//  		res(globalCenter) = (flux_y[node.interBack] +nit->value)*scaling;


		int globalNeighbor=globalBack;
		int localNeighbor = node.back;

		real accumulation=(rho(globalNeighbor)*Dtheta(globalNeighbor) 
				   + Drho(globalNeighbor)*theta(globalNeighbor));


		//use Dp as global variable
		accumulation *= Dp(globalNeighbor);

		real volume = detMap(globalNeighbor)*hYcen[localNeighbor]
		  *hXcen[localNeighbor];

		real fluxCalc = flux_y[node.interFront];

		//now switch to localNeighbor to get divergence
		node.localNode(localNeighbor);

		//normal res calc uses -div from DivKgradECDM
		res(globalCenter) = volume*accumulation
		  + flux_x[node.interRight] - flux_x[node.interLeft] 
		  + flux_y[node.interBack] -fluxCalc;



		//              cout<<"FF"<<'\t'<<flux_y[node.interFront]<<'\t'<<nit->value <<'\t'<< flux_y[node.interBack]<<endl;
		break;
	      }
	    case BACK:
	      {
		//  		res(globalCenter) = (flux_y[node.interFront] -nit->value)*scaling;

		int globalNeighbor=globalFront;
		int localNeighbor = node.front;

		real accumulation=(rho(globalNeighbor)*Dtheta(globalNeighbor) 
				   + Drho(globalNeighbor)*theta(globalNeighbor));


		//use Dp as global variable
		accumulation *= Dp(globalNeighbor);

		real volume = detMap(globalNeighbor)*hYcen[localNeighbor]
		  *hXcen[localNeighbor];

		real fluxCalc = flux_y[node.interBack];

		std::cout<<" in appNeuRes glbno = "<<globalCenter<<" glbnei= "
			 <<globalNeighbor<<" vol= "<<volume
			 <<" accu = "<<accumulation
			 <<" fluxCalc= "<<fluxCalc
			 <<" res("<<globalNeighbor<<")= "
			 <<res(globalNeighbor)<<std::endl;

		//now switch to localNeighbor to get divergence
		node.localNode(localNeighbor);

		//normal res calc uses -div from DivKgradECDM
		res(globalCenter) = volume*accumulation
		  + flux_x[node.interRight] - flux_x[node.interLeft] 
		  + fluxCalc - flux_y[node.interFront];

		std::cout<<" res("<<globalCenter<<")= "<<res(globalCenter)
			 <<std::endl;
		break;

	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch

	  res(globalCenter) *= scaling;
	}//end if isLocal
      ++nit;
    }
}

template <class JAC, int nv>
inline void 
ConstBCECDM::applyNeumannDerivativesToRes(Petsc::StencilMM& node, JAC& DresDy,
					  VecVec& D_Div, int DIM,
					  DivKgradECDM<ConstBCECDM,nv>* div,
					  const Vec& rho, const Vec& theta,
					  const Vec& Drho, const Vec& Dtheta,
					  const Vec& detMap, const Vec& Dp, 
					  real alphaBDF, const Vec& DDrho, 
					  const Vec& DDtheta,
					  const Vec& hXcen, const Vec& hYcen,
					  real scaling)

{

  nit = neumann.begin();
  while(nit != neumann.end()) 
    {

      node.globalNode(nit->n);
      int globalCenter = node.center;
      int globalLeft  = node.left;
      int globalRight = node.right;
      int globalFront = node.front;
      int globalBack  = node.back;

      int globalNeighbor(-1),localNeighbor(-1);

      if (node.isLocal)
	{
	  real boundJacTerm(0.0);
	  int  offset = node.center-node.globalLow;

	  int localNodeNumber = node.globalToLocal(nit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);

	  switch(nit->face)
	    {
	    case LEFT:
	      {
		globalNeighbor = globalRight;
		localNeighbor  = node.right;

		break;
	      }
	    case RIGHT:
	      {
		globalNeighbor = globalLeft;
		localNeighbor  = node.left;

		break;
	      }
	    case FRONT:
	      {
		globalNeighbor = globalBack;
		localNeighbor  = node.back;

		break;
	      }
	    case BACK:
	      {
		globalNeighbor = globalFront;
		localNeighbor  = node.front;
		boundJacTerm = 	 
		  D_Div[DivKgradECDM<ConstBCECDM,1>::BACK][offset];

		break;
	      }
	    }//end first switch for local
	
	  //first get jacobian term for actual boundary cell

	  //mwf now go back and center stencil on interior cell
	  node.globalNode(globalNeighbor);
	  int nl=node.center-node.globalLow;
	  //use same sign convention for divergence as in regular res
	  real signOfDiv(-1.0);

	  real volume = hXcen[nl]*hYcen[nl]*fabs(detMap[nl]);

	  //put in scaling here?
	  signOfDiv *= scaling;
	  volume    *= scaling;

	  real ypJac = volume*(rho[nl]*Dtheta[nl] + Drho[nl] * theta[nl]);
	  real Daccum= volume*(rho[nl]*DDtheta[nl] +
			       Drho[nl]*Dtheta[nl] +
			       DDrho[nl]*theta[nl] +
			       Drho[nl]*Dtheta[nl])* Dp[nl];

	  //go through and evaluate jacobian. The stencil is now centered
	  //at interior neighbor

          DresDy(globalCenter,node.center)  = Daccum + alphaBDF*ypJac + 
	    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::CENTER][nl];
          
          if (node.anchor->k > 0)
            DresDy(globalCenter,globalLeft)  = 
	     signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT][nl];

          if (node.anchor->k < node.nxNodes-1)
            DresDy(globalCenter,node.right) = 
	      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT][nl];
          
          if (DIM > 0)
            {
              if (node.anchor->j > 0)
		{
		  DresDy(globalCenter,node.front)= 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::FRONT][nl];
		  std::cerr<<" in appNeu DresDy("<<globalCenter<<","
			   <<node.front<<")= "<<signOfDiv
		    *D_Div[DivKgradECDM<ConstBCECDM,1>::FRONT][nl]<<std::endl;

		  //(j-1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(globalCenter,node.leftFront)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_FRONT][nl];
		  //(j-1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(globalCenter,node.rightFront) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_FRONT][nl];

		}
              if (node.anchor->j < node.nyNodes-1)
                {
		  DresDy(globalCenter,node.back) = 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BACK][nl];

		  std::cerr<<" in appNeu DresDy("<<globalCenter<<","
			   <<node.back<<")= "<<signOfDiv
		    *D_Div[DivKgradECDM<ConstBCECDM,1>::FRONT][nl]<<std::endl;

		  //(j+1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(globalCenter,node.leftBack)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_BACK][nl];
		  //(j+1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(globalCenter,node.rightBack) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_BACK][nl];
		  
		}
	    }//end DIM==TWO_D
          if (DIM > 1) //DIM==THREE_D
            {
              if (node.anchor->i > 0)
                DresDy(globalCenter,node.bottom)= 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BOTTOM][nl];
              if (node.anchor->i < node.nzNodes-1)
                DresDy(globalCenter,node.top) = 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::TOP][nl];
            }


	  DresDy(globalCenter,globalCenter)  = 
	    signOfDiv*boundJacTerm;

	  DresDy.finalizeRow(globalCenter);

	  DresDy(globalNeighbor,globalCenter)  = 0.0;

	  //no dependence on actual boundary node?
	  //how about zeroing neighbors instead?

	  DresDy.finalizeRow(globalNeighbor);

        }//end isLocal
      ++nit;
    }//end while
}//end function

//========= mass conservative form 
//mwf now try and put in approx. using mass bal. only 2d for now
template <int nv>
inline void 
ConstBCECDM::applyNeumannConditionsToRes(Petsc::StencilMM& node, Vec& res,
					 DivKgradECDM<ConstBCECDM,nv>* div,
					 const Vec& mCurrent, 
					 const Vec& betaDaeDef, 
					 const Vec& detMap,
					 real alphaBDF,
					 real scaling)
{
  assert(nv == 1);
  assert(div);
  const Vec & flux_x = div->getFlux_x(0);
  const Vec & flux_y = div->getFlux_y(0);
  const Vec & hXcen  = div->getSpatialIncrement(0);
  const Vec & hYcen  = div->getSpatialIncrement(1);

  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(nit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  //mwf this looks like it might be a problem
	  real accumulation=alphaBDF*mCurrent(globalCenter) + betaDaeDef(2*globalCenter+1);
	    
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		//remember sign of normal for minus side 
		//res(globalCenter) = (flux_x[node.interRight] + nit->value)*scaling;
		//hXcen should be zero I believe
		real volume = detMap(globalCenter)*hYcen[node.center]
		  *0.5*(hXcen[node.center]+hXcen[node.right]);

		res(globalCenter) = volume*accumulation
		  + flux_x[node.interRight] + nit->value;

		//mwf no scaling for now

		if (fabs(res(globalCenter)) > 1.0e4)
		  {
		    cerr<<"in ConstBCECDM::neumannBCsToRes, node.center= "
			<<node.center<<endl;
		    cerr<<"\t flux_x["<<node.interRight<<"]= "
			<<flux_x[node.interRight]<<endl;
		    cerr<<"\t scaling = "<<scaling<<" nit->value= "<<nit->value<<endl;
		  } 
		break;
	      }
	    case RIGHT:
	      {
//  		res(globalCenter) = (flux_x[node.interLeft] - nit->value)*scaling;
		//hXcen should be zero I believe
		real volume = detMap(globalCenter)*hYcen[node.center]
		  *0.5*(hXcen[node.center]+hXcen[node.left]);

		res(globalCenter) = volume*accumulation
		  - flux_x[node.interLeft] + nit->value;


		//mwf no scaling for now
		if (fabs(res(globalCenter)) > 1.0e4)
		  {
		    cerr<<"in ConstBCECDM::neumannBCsToRes, node.center= "
			<<node.center<<endl;
		    cerr<<"\t flux_x["<<node.interLeft<<"]= "
			<<flux_x[node.interLeft]<<endl;
		    cerr<<"\t scaling = "<<scaling<<" nit->value= "<<nit->value<<endl;
		  } 
		break;
	      }
	    case FRONT:
	      {
		//remeber sign for negative boundary
//  		res(globalCenter) = (flux_y[node.interBack] +nit->value)*scaling;

		//hYcen should be zero I believe
		real volume = detMap(globalCenter)*hXcen[node.center]
		  *0.5*(hYcen[node.center]+hYcen[node.back]);

		res(globalCenter) = volume*accumulation
		  + flux_y[node.interBack] + nit->value;


		//mwf no scaling for now

		//              cout<<"FF"<<'\t'<<flux_y[node.interFront]<<'\t'<<nit->value <<'\t'<< flux_y[node.interBack]<<endl;
		break;
	      }
	    case BACK:
	      {
//  		res(globalCenter) = (flux_y[node.interFront] -nit->value)*scaling;
		//hYcen should be zero I believe
		real volume = detMap(globalCenter)*hXcen[node.center]
		  *0.5*(hYcen[node.center]+hYcen[node.front]);

		res(globalCenter) = volume*accumulation  
		  -  flux_y[node.interFront] + nit->value;


		//mwf no scaling for now
		//              cout<<"FB"<<'\t'<<flux_y[node.interBack] <<'\t'<<nit->value<<'\t'<< flux_y[node.interFront]<<endl;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }
}

template <class JAC, int nv>
inline void 
ConstBCECDM::applyNeumannDerivativesToRes(Petsc::StencilMM& node, JAC& DresDy,
					  VecVec& D_Div, int DIM,
					  DivKgradECDM<ConstBCECDM,nv>* div,
					  const Vec& rho, const Vec& theta,
					  const Vec& Drho, const Vec& Dtheta,
					  const Vec& detMap, 
					  real alphaBDF, 
					  real scaling)

{
  const Vec & hXcen  = div->getSpatialIncrement(0);
  const Vec & hYcen  = div->getSpatialIncrement(1);


  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      real volume(1.0);

      node.globalNode(nit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(nit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);

	  switch(nit->face)
	    {
	    case LEFT:
	      {
		volume = detMap(globalCenter)*hYcen[node.center]
		  *0.5*(hXcen[node.center]+hXcen[node.right]);
		break;
	      }
	    case RIGHT:
	      {
		volume = detMap(globalCenter)*hYcen[node.center]
		  *0.5*(hXcen[node.center]+hXcen[node.left]);
		break;
	      }
	    case FRONT:
	      {
		volume = detMap(globalCenter)*hXcen[node.center]
		  *0.5*(hYcen[node.center]+hYcen[node.back]);
		break;
	      }
	    case BACK:
	      {
		volume = detMap(globalCenter)*hXcen[node.center]
		  *0.5*(hYcen[node.center]+hYcen[node.front]);
		break;
	      }
	    }//end first switch for local
	}//end first isLocal

      //mwf now go back through like usual
      node.globalNode(nit->n);
      int nl=node.center-node.globalLow;
      real signOfDiv(1.0);
      if (node.isLocal)
	{
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		signOfDiv = -1.0;
		break;
	      }
	    case RIGHT:
	      {
		signOfDiv = -1.0;
		break;
	      }
	    case FRONT:
	      {
		signOfDiv = -1.0;
		break;
	      }
	    case BACK:
	      {
		signOfDiv = -1.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	  //mwf no scaling for now signOfDiv *= scaling;
	  real ypJac = volume*(rho[nl]*Dtheta[nl] + 
				Drho[nl] * theta[nl]);

          DresDy(node.center,node.center)  = alphaBDF*ypJac + 
	    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::CENTER][nl];
          
          if (node.anchor->k > 0)
            DresDy(node.center,node.left)  = 
	     signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT][nl];

          if (node.anchor->k < node.nxNodes-1)
            DresDy(node.center,node.right) = 
	      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT][nl];
          
          if (DIM > 0)
            {
              if (node.anchor->j > 0)
		{
		  DresDy(node.center,node.front)= 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::FRONT][nl];
		  //(j-1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(node.center,node.leftFront)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_FRONT][nl];
		  //(j-1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(node.center,node.rightFront) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_FRONT][nl];

		}
              if (node.anchor->j < node.nyNodes-1)
                {
		  DresDy(node.center,node.back) = 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BACK][nl];
		  //(j+1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(node.center,node.leftBack)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_BACK][nl];
		  //(j+1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(node.center,node.rightBack) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_BACK][nl];
		  
		}
	    }//end DIM==TWO_D
          if (DIM > 1) //DIM==THREE_D
            {
              if (node.anchor->i > 0)
                DresDy(node.center,node.bottom)= 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BOTTOM][nl];
              if (node.anchor->i < node.nzNodes-1)
                DresDy(node.center,node.top) = 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::TOP][nl];
            }
	  DresDy.finalizeRow(node.center);
        }//end isLocal
      ++nit;
    }//end while
}//end function
//====================robbins===========================

//mwf I don't know yet what to do about robbins and interior bc's
//as far as time integration goes
inline
void  ConstBCECDM::applyRobbinsConditionsYprime(Petsc::StencilMM& node, 
						Vec& yPrime)
{

  ruit = robbins.begin();
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  yPrime(globalCenter) = 0.0;
	}
      ++ruit;
    }
}

//make global residual be dirMult*p + neuMult*fluxInterior-rVal
template <int nv>
inline 
void  ConstBCECDM::applyRobbinsConditionsToRes(Petsc::StencilMM& node, 
					       Vec& res,
					       DivKgradECDM<ConstBCECDM,nv>* 
					       div,
					       real scaling)
{
  assert(nv == 1);
  assert(div);
  const Vec & flux_x = div->getFlux_x(0);
  const Vec & flux_y = div->getFlux_y(0);

  ruit = robbins.begin();
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      int globalCenter = node.center;
      real resTerm(0.0);

      if (node.isLocal)
	{
	  //do dirichlet part first
	  resTerm = (*ruit->var)(globalCenter)*(ruit->dirMult);

	  int localNodeNumber = node.globalToLocal(ruit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  //now go through and get face values for neumann part
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		//using interior flux for the Dirichlet face
		//and zeroing the exterior flux
		//remember sign of normal for minus side 
		resTerm -= (ruit->neuMult)*flux_x[node.interRight];
		break;
	      }
	    case RIGHT:
	      {
		resTerm += (ruit->neuMult)*flux_x[node.interLeft];
		break;
	      }
	    case FRONT:
	      {
		//remeber sign for negative boundary
		resTerm -= (ruit->neuMult)*flux_y[node.interBack];
		break;
	      }
	    case BACK:
	      {
		resTerm += (ruit->neuMult)*flux_y[node.interFront];
		//if (globalCenter >= 1069 && globalCenter <= 1082)
		//  {
		//    std::cout<<"in BCrob value= "<<ruit->value
		//	     <<" flux_y["<<node.interFront<<"]= "
		//	     <<flux_y[node.interFront]<<" dirMult,neuMult= "
		//	     <<ruit->dirMult<<" "<<ruit->neuMult<<std::endl;
		//  }//end if for debugging
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	  //now load in specified value
	  res(globalCenter) = (resTerm - ruit->value)*scaling;

	}//end if isLocal
      ++ruit;
    }
}

//make global residual be dirMult*p + neuMult*fluxInterior-rVal
template <class JAC>
inline 
void  ConstBCECDM::applyRobbinsDerivativesToRes(Petsc::StencilMM& node, 
						JAC& DresDy,
						VecVec& D_Div, int DIM,
						real scaling)
{

  ruit = robbins.begin();
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      int nl=node.center-node.globalLow;
      real signOfDiv(1.0);
      if (node.isLocal)
	{
	  //do dirichlet part first 
	  real dirichletDeriv = (ruit->dirMult)*scaling;
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		signOfDiv = 1.0;
		break;
	      }
	    case RIGHT:
	      {
		signOfDiv = 1.0;
		break;
	      }
	    case FRONT:
	      {
		signOfDiv = 1.0;
		break;
	      }
	    case BACK:
	      {
		signOfDiv = 1.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	  signOfDiv *= (ruit->neuMult)*scaling;
	      
          DresDy(node.center,node.center)  = dirichletDeriv + 
	    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::CENTER][nl];
          
          if (node.anchor->k > 0)
            DresDy(node.center,node.left)  = 
	     signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT][nl];

          if (node.anchor->k < node.nxNodes-1)
            DresDy(node.center,node.right) = 
	      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT][nl];
          
          if (DIM > 0)
            {
              if (node.anchor->j > 0)
		{
		  DresDy(node.center,node.front)= 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::FRONT][nl];
		  //(j-1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(node.center,node.leftFront)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_FRONT][nl];
		  //(j-1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(node.center,node.rightFront) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_FRONT][nl];

		}
              if (node.anchor->j < node.nyNodes-1)
                {
		  DresDy(node.center,node.back) = 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BACK][nl];
		  //(j+1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(node.center,node.leftBack)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_BACK][nl];
		  //(j+1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(node.center,node.rightBack) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_BACK][nl];
		  
		}
	    }//end DIM==TWO_D
          if (DIM > 1) //DIM==THREE_D
            {
              if (node.anchor->i > 0)
                DresDy(node.center,node.bottom)= 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BOTTOM][nl];
              if (node.anchor->i < node.nzNodes-1)
                DresDy(node.center,node.top) = 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::TOP][nl];
            }
	  DresDy.finalizeRow(node.center);
        }//end isLocal
      ++ruit;
    }//end while
}//end function

//====================no flow===========================

//mwf I don't know yet what to do about neumann  bc's
//as far as time integration goes
inline
void  ConstBCECDM::applyNoFlowConditionsYprime(Petsc::StencilMM& node, 
					       Vec& yPrime)
{

  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  yPrime(globalCenter) = 0.0;
	}
      ++noit;
    }
}

//make global residual be dirMult*p + neuMult*fluxInterior-rVal
template <int nv>
inline 
void  ConstBCECDM::applyNoFlowConditionsToRes(Petsc::StencilMM& node, 
					      Vec& res,
					      DivKgradECDM<ConstBCECDM,nv>* 
					      div,
					      real scaling)
{
  assert(nv == 1);
  assert(div);

  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      int globalCenter = node.center;

      if (node.isLocal)
	{

	  real fluxVal(0.0);
	  int localNodeNumber = node.globalToLocal(noit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  //now go through and get face values for neumann part
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		//using single phase flux for gradient
		//remember sign of normal for minus side 
		fluxVal = div->getSinglePhaseFlux_xLeftBnd();
		res(globalCenter)= (fluxVal + noit->value)*scaling;
		break;
	      }
	    case RIGHT:
	      {
		fluxVal = div->getSinglePhaseFlux_xRightBnd();
		res(globalCenter)= (fluxVal - noit->value)*scaling;
		break;
	      }
	    case FRONT:
	      {
		//remeber sign for negative boundary
		fluxVal = div->getSinglePhaseFlux_yFrontBnd();
		res(globalCenter)= (fluxVal + noit->value)*scaling;
		break;
	      }
	    case BACK:
	      {
		fluxVal = div->getSinglePhaseFlux_yBackBnd();
		res(globalCenter)= (fluxVal - noit->value)*scaling;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++noit;
    }//end switch
}//end function

//make global residual be dirMult*p + neuMult*fluxInterior-rVal
template <class JAC, int nv>
inline 
void  ConstBCECDM::applyNoFlowDerivativesToRes(Petsc::StencilMM& node, 
					       JAC& DresDy,
					       DivKgradECDM<ConstBCECDM,nv>* 
					       div, int DIM,
					       real scaling)
{

  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      //int nl=node.center-node.globalLow;

      if (node.isLocal)
	{
	  //save these for use with jacobian
	  int globalCenter    = node.center;
	  int globalLeft      = node.left;
	  int globalRight     = node.right;
	  int globalFront     = node.front;
	  int globalBack      = node.back;
	  int globalLeftBack  = node.leftBack;
	  int globalRightBack = node.rightBack;
	  int globalLeftFront = node.leftFront;
	  int globalRightFront= node.rightFront;
	  int globalBackLeft  = node.leftBack;
	  int globalBackRight = node.rightBack;
	  int globalFrontLeft = node.leftFront;
	  int globalFrontRight= node.rightFront;

	  int localNodeNumber = node.globalToLocal(noit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);

	  switch(noit->face)
	    {
	    case LEFT:
	      {
		real DfluxVal(0.0);
		DfluxVal = div->getDSinglePhaseFlux_xLeftBnd_center();
		DresDy(globalCenter,globalCenter) = scaling*DfluxVal;
		
		DfluxVal = div->getDSinglePhaseFlux_xLeftBnd_right();
		DresDy(globalCenter,globalRight) = scaling*DfluxVal;
		          
		if (DIM > 0)
		  {

		    DfluxVal = div->getDSinglePhaseFlux_xLeftBnd_rightBack();
		    DresDy(globalCenter,globalRightBack) = scaling*DfluxVal;

		    DfluxVal = div->getDSinglePhaseFlux_xLeftBnd_rightFront();
		    DresDy(globalCenter,globalRightFront) = scaling*DfluxVal;
		  }
		DresDy.finalizeRow(globalCenter);
		break;
	      }
	    case RIGHT:
	      {
		real DfluxVal(0.0);
		DfluxVal = div->getDSinglePhaseFlux_xRightBnd_center();
		DresDy(globalCenter,globalCenter) = scaling*DfluxVal;
		
		DfluxVal = div->getDSinglePhaseFlux_xRightBnd_left();
		DresDy(globalCenter,globalLeft) = scaling*DfluxVal;
		
		          
		if (DIM > 0)
		  {
		    
		    DfluxVal = div->getDSinglePhaseFlux_xRightBnd_leftBack();
		    DresDy(globalCenter,globalLeftBack) = scaling*DfluxVal;
		    
		    DfluxVal = div->getDSinglePhaseFlux_xRightBnd_leftFront();
		    DresDy(globalCenter,globalLeftFront) = scaling*DfluxVal;
		  }

		DresDy.finalizeRow(globalCenter);

		break;
	      }
	    case FRONT:
	      {
		real DfluxVal(0.0);
		DfluxVal = div->getDSinglePhaseFlux_yFrontBnd_center();
		DresDy(globalCenter,globalCenter) = scaling*DfluxVal;

		DfluxVal = div->getDSinglePhaseFlux_yFrontBnd_back();
		DresDy(globalCenter,globalBack) = scaling*DfluxVal;

		DfluxVal = div->getDSinglePhaseFlux_yFrontBnd_backRight();
		DresDy(globalCenter,globalBackRight) = scaling*DfluxVal;

		DfluxVal = div->getDSinglePhaseFlux_yFrontBnd_backLeft();
		DresDy(globalCenter,globalBackLeft) = scaling*DfluxVal;

		DresDy.finalizeRow(globalCenter);

		break;
	      }
	    case BACK:
	      {
		real DfluxVal(0.0);
		DfluxVal = div->getDSinglePhaseFlux_yBackBnd_center();
		DresDy(globalCenter,globalCenter) = scaling*DfluxVal;

		DfluxVal = div->getDSinglePhaseFlux_yBackBnd_front();
		DresDy(globalCenter,globalFront) = scaling*DfluxVal;

		DfluxVal = div->getDSinglePhaseFlux_yBackBnd_frontRight();
		DresDy(globalCenter,globalFrontRight) = scaling*DfluxVal;

		DfluxVal = div->getDSinglePhaseFlux_yBackBnd_frontLeft();
		DresDy(globalCenter,globalFrontLeft) = scaling*DfluxVal;

		DresDy.finalizeRow(globalCenter);

		break;

	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++noit;
    }//end switch


}//end function

//====================interior===========================

//make global residual be fluxLeft-fluxRight
template <int nv>
inline 
void  ConstBCECDM::applyInteriorConditionsToRes(Petsc::StencilMM& node, 
						Vec& res,
						DivKgradECDM<ConstBCECDM,nv>* div,
						real scaling)
{
  assert(nv == 1);
  assert(div);

  const   Vec & flux_x = div->getFlux_x(0);
  const   Vec & flux_y = div->getFlux_y(0);
  init = interior.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      int globalCenter = node.center;
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
		//right flux is headed in negative direction
		res(globalCenter) = (flux_x[node.interLeft]
				     -flux_x[node.interRight])*scaling;
		break;
	      }
	    case INTERIOR_Y:
	      {
		//right flux is headed in negative direction
		res(globalCenter) = (flux_y[node.interFront]
				     -flux_y[node.interBack])*scaling;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<init->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++init;
    } //end init loop

}//end function


//make global residual be fluxInterior-qVal
template <class JAC>
inline 
void  ConstBCECDM::applyInteriorDerivativesToRes(Petsc::StencilMM& node, 
						 JAC& DresDy,
						 VecVec& D_Div, int DIM,
						 real scaling)
{

  init = interior.begin();
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      int nl=node.center-node.globalLow;
      real signOfDiv(1.0);
      signOfDiv *= scaling;
      if (node.isLocal)
	{

	  DresDy(node.center,node.center)  = 
	    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::CENTER][nl];
	  
	  if (node.anchor->k > 0)
	    DresDy(node.center,node.left)  = 
	      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT][nl];
	  
	  if (node.anchor->k < node.nxNodes-1)
	    DresDy(node.center,node.right) = 
	  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT][nl];
	  
	  if (DIM > 0)
	    {
	      if (node.anchor->j > 0)
		{
		  DresDy(node.center,node.front)= 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::FRONT][nl];
		  //(j-1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(node.center,node.leftFront)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_FRONT][nl];
		  //(j-1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(node.center,node.rightFront) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_FRONT][nl];
		  
		}
	      if (node.anchor->j < node.nyNodes-1)
		{
		  DresDy(node.center,node.back) = 
		    signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BACK][nl];
		  //(j+1,k-1)
		  if (node.anchor->k > 0)
		    DresDy(node.center,node.leftBack)=
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::LEFT_BACK][nl];
		  //(j+1,k+1)
		  if (node.anchor->k < node.nxNodes-1)
		    DresDy(node.center,node.rightBack) = 
		      signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::RIGHT_BACK][nl];
		  
		}
	    }//end DIM==TWO_D
	  if (DIM > 1) //DIM==THREE_D
	    {
	      if (node.anchor->i > 0)
		DresDy(node.center,node.bottom)= 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::BOTTOM][nl];
	      if (node.anchor->i < node.nzNodes-1)
		DresDy(node.center,node.top) = 
		  signOfDiv*D_Div[DivKgradECDM<ConstBCECDM,1>::TOP][nl];
	    }//end 3d
	  DresDy.finalizeRow(node.center);
	}//end isLocal
      ++init;
   }//end while
}//end function
//this one updates derivatives straight into jacobian for left and right terms  
template <class JAC, int nv>
inline 
void  
ConstBCECDM::updateInteriorDerivativesToRes(Petsc::StencilMM& node, 
					    JAC& DresDy,
					    const VecVec* K, const Vec* Kr,
					    const Vec* Rho, const Vec* P,
					    const VecVecVec& DKr, 
					    const Vec* DRho,
					    VecVec& D_Div, int DIM,
					    DivKgradECDM<ConstBCECDM,nv>* div,
					    real scaling)
{

  init = interior.begin();
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //int nl=node.center-node.globalLow;
      real signOfDiv(1.0);
      signOfDiv *= scaling;

      if (node.isLocal)
	{
	  //save these for use with jacobian
	  //int globalCenter    = node.center;
	  int globalLeft      = node.left;
	  int globalRight     = node.right;
	  int globalFront     = node.front;
	  int globalBack      = node.back;

	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);


	  //set jacobian for neighbors on 
	  //either side of interface. They each need an extra
	  //term to express dependence of rel perm coefficient
	  //on neighbors
	  //this goes into the equations for the left and right neighbors
	  //or the back and front (Y)
	  if (init->face == INTERIOR_X)
	    {	  
	      real fluxDeriv(0.0);
	      int face(0);
	      //interface coefficient is upwind(Kr_{left},Kr_{right})*
	      //1/2(rho_{left}+rho_{right})
	      
	      //from the right
	      //flux gets used as left flux in right mass eqn
	      face = 1;
	      fluxDeriv = div->getInteriorFluxDerivTerm_x(face);

	      //only way get this derivative
	      DresDy(globalRight,globalLeft) = -fluxDeriv;
	      DresDy.finalizeRow(globalRight);

	      //from the left
	      face = -1;
	      //flux gets used as right flux in right mass eqn
	      fluxDeriv = div->getInteriorFluxDerivTerm_x(face);

	      DresDy(globalLeft,globalRight) = fluxDeriv;
	      DresDy.finalizeRow(globalLeft);
	    }//end x face
	  if (init->face == INTERIOR_Y)
	    {	  
	      //set boundary values with div routine, don't use rel perm or density
	      //since they're continuous anyway
	      real fluxDeriv(0.0);
	      int face(0);

	      //interface coefficient is upwind(Kr_{back},Kr_{front})*
	      //1/2(rho_{back}+rho_{front})

	      //from the back
	      //flux gets used as front flux in back mass eqn
	      face = 1;
	      fluxDeriv = div->getInteriorFluxDerivTerm_y(face);

	      DresDy(globalBack,globalFront) = -fluxDeriv;
	      DresDy.finalizeRow(globalBack);

	      //from the front
	      face = -1;
	      fluxDeriv = div->getInteriorFluxDerivTerm_y(face);

	      DresDy(globalFront,globalBack) = fluxDeriv;
	      DresDy.finalizeRow(globalFront);
	    }//end Y face
	      
	}//end isLocal
      ++init;
   }//end while
}//end function

//mwf I don't know yet what to do about neumann and interior bc's
//as far as time integration goes
inline
void  ConstBCECDM::applyInteriorConditionsYprime(Petsc::StencilMM& node, 
						 Vec& yPrime)
{

  init = interior.begin();
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  yPrime(globalCenter) = 0.0;
	}
      ++init;
    }
}

//////////////////////
//mwf now I have to use the stupid STL algorithms I guess
template <class BcType>
class BCnodeIsTheSame
{
public:
  typedef typename std::vector<BcType>::iterator first_argument_type ;
  typedef int second_argument_type ;
  typedef bool result_type ;

  BCnodeIsTheSame(int num = 0):gnum(num) {}

  bool operator()(BcType& iterVal)
    { return iterVal.n == gnum; }

  int gnum;
};


}//Daetk

#endif











