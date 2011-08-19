#ifndef DIVKGRAD_ECDM_H
#define DIVKGRAD_ECDM_H
#include "Divergence.h"
#include "PetscStencilMM.h"

namespace Daetk 
{
using Petsc::StencilMM;

//  typedef CMRVec<Vec> VecVec;
//  typedef CMRVec<VecVec> VecVecVec;
//  typedef CMRVec<VecVecVec> VecVecVecVec;
typedef std::vector<Vec> VecVec;
typedef std::vector<VecVec> VecVecVec;
typedef std::vector<VecVecVec> VecVecVecVec;


//  FIX DIVERGENCE JAC -- MUST APPLY CORRECT DERIVATIVES FOR NEUMAN BC'S
//#define USE_ARITHMETIC_AVG_ECDM
//#define USE_CONT_INTERFACE_TERMS_ECDM

template<class BC, int nv>
class DivKgradECDM : public Divergence
{
public:
  void computeInterfaceK(const Vec&){std::cout<<"DivKgradECDM::computeInterfaceK has not been implemented"<<std::endl;}
 
  enum JacEntry {CENTER, LEFT, RIGHT, FRONT, BACK, 
		 RIGHT_FRONT, RIGHT_BACK, LEFT_FRONT, LEFT_BACK,
		 BOTTOM, TOP};
  //for permeability tensor
  enum KmatEntry {XX,YY,XY,YX,ZZ,XZ,YZ};
  enum KmatCorner{MINUS,PLUS};

  DivKgradECDM(BC* bcIn, StencilMM& nodeIn):bc(bcIn),node(nodeIn) {}
  
  virtual ~DivKgradECDM(){}

  virtual void computeDivergence(const Vec* K, 
				 const Vec* Rho, 
				 const Vec* P, 
				 Vec* Div)=0;

  virtual void computeDivergence(const VecVec* K, const Vec* Kr,
				 const Vec* Rho, 
				 const Vec* P, 
				 Vec* Div)
    {
      Vec Kscalar[nv];
      for (vi=0; vi < nv; vi++);
      {
	Kscalar[vi] = K[vi][XX];
      }
      computeDivergence(Kscalar,Rho,P,Div); 
    }

  virtual void computeDivergenceJac(const Vec* K,
				    const Vec* Rho, 
				    const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)=0;  
  
  virtual void computeDivergenceJac(const VecVec* K, const Vec* Kr,
				    const Vec* Rho, 
				    const Vec* P,
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
    {
      Vec Kscalar[nv];
      for (vi=0; vi < nv; vi++);
      {
	Kscalar[vi] = K[vi][XX];
      }
      computeDivergenceJac(Kscalar,Rho,P,DK,DRho,DivJac) ; 
    }
  //arithmetic average
  inline real avgA(const real& a, const real& b)
    { return 0.5*(a+b); }
  //harmonic average
  inline real avgH(const real& a, const real& b)
    { 
      if (a == 0.0 && b == 0.0)
	return 0.0;
      else
	return 2.0*a*b/(a+b); 
    }
  //upwinding (a if uVar >= 0) -- trick: what's uVar? 
  inline real avgU(const real& a, const real& b,const real& uVar)
    { 
      real eps = 1.0e-32;
      real uVarAbs = fabs(uVar);
      if (uVar >= eps*uVarAbs) 
	return a;
      else
	return b;
    }
  //weighted average (using cell increments)
  //note wb comes before wa
  //this one gives you value on this side
  //of interior boundary
  inline real avgWA2(const real& a, const real& b,
		    const real& wb, const real& wa)
    { 
      //mwf just test this to see 
      real eps = 1.0e-32;
      if (fabs(wa) <= eps)
	return a;
      if (fabs(wb) <= eps)
	return b;
      else
	return (a*wa + b*wb)/(wb + wa); 
    }

  //put these in for intrinsic perm
  inline real KIx(const Vec& K)
    { 
//        real harmAvg = avgWH(K[node.center],K[node.right],
//  			   hXcen[node.center],hXcen[node.right]);
      real aritAvg = avgWA2(K[node.center],K[node.right],
			    hXcen[node.center],hXcen[node.right]);
      return aritAvg;
      //return harmAvg;
      //return avgH(K[node.center],K[node.right]);
    }
  //put these in for intrinsic perm
  inline real KIy(const Vec& K)
    { 
//        real harmAvg = avgWH(K[node.center],K[node.back],
//  			   hYcen[node.center],hYcen[node.back]);
      real aritAvg = avgWA2(K[node.center],K[node.back],
			   hYcen[node.center],hYcen[node.back]);
      return aritAvg;
//        return harmAvg;
      //return avgH(K[node.center],K[node.back]);
    }

  //relative permeability for interior boundary
  //can't use pDiff here because not calculated with latest P yet
  inline real avgRelPermAcrossInteriorBoundaryX(const Vec& Kr,const Vec& Rho,
						const Vec& P) 
    {  
#ifdef USE_ARITHMETIC_AVG_ECDM
      return  avgA(Kr[node.left],Kr[node.right]);
#else
      //mwf have to use whole driving force to do upwinding
      //mwf use these for now even though they're out of date?
      real driveForce = getSinglePhaseFlux_x();
     //this worked for the most part
//        real driveForce = hYcen[node.center]*(-DpDx(P))
//  	+ hYcen[node.center]*Rho[node.center]*gravX[node.interRightGhost];
      return  avgU(Kr[node.left],Kr[node.right],driveForce);

#endif
    }
  inline real avgRelPermAcrossInteriorBoundaryY(const Vec& Kr,const Vec& Rho,
						const Vec& P) 
    {  

#ifdef USE_ARITHMETIC_AVG_ECDM
      return  avgA(Kr[node.front],Kr[node.back]);
#else
      //mwf use these for now even though they're out of date?
      real driveForce = getSinglePhaseFlux_y();
//        real driveForce =  hXcen[node.center]*(-DpDy(P))
//  	+ hXcen[node.center]*Rho[node.center]*gravY[node.interBackGhost];

      return avgU(Kr[node.front],Kr[node.back],driveForce);
#endif
    }

  inline real avgDensityAcrossInteriorBoundaryX(const Vec& Rho,const Vec& P) 
    {  
      return  avgA(Rho[node.left],Rho[node.right]);
    }
  inline real avgDensityAcrossInteriorBoundaryY(const Vec& Rho,const Vec& P) 
    {  
      return  avgA(Rho[node.front],Rho[node.back]);
    }

  inline real DpDx(const Vec& P)
    { 
      return 2.0*(P[node.right] - P[node.center])
	/(hXcen[node.center]+hXcen[node.right]);
    }
  inline real DpDy(const Vec& P)
  { 
    //cout<<"DpDy "<<hYint[node.interBack]<<'\t'<<P[node.back]
    //	<<'\t'<<P[node.center]<<endl;
    return 2.0*(P[node.back] - P[node.center])
      /(hYcen[node.center] + hYcen[node.back]);
  }
  inline real DpDz(const Vec& P)
    { 
      return 2.0*(P[node.top] - P[node.center])
	/(hZcen[node.center]+hZcen[node.front]);
    }

  //assumes that RhoX has been set
  inline virtual void setPdiff_x(const Vec* Rho, const Vec* P)
    {                  
      pDiff_x[0][node.interRightGhost] = 
	hYcen[node.center]*(-DpDx(P[0]))
	+ hYcen[node.center]
	*RhoX[node.interRightGhost]*gravX[node.interRightGhost];
    }

  //assumes that RhoY has been set
  inline virtual void setPdiff_y(const Vec* Rho, const Vec* P)
    {
      pDiff_y[0][node.interBackGhost] = 
	  hXcen[node.center]*(-DpDy(P[0]))
	+ hXcen[node.center]
	*RhoY[node.interBackGhost]*gravY[node.interBackGhost];

//      cout<<"y "<<node.center<<'\t'<<node.interBack<<'\t'
//  	<<node.back<<endl<<flush;
//      cout<<pDiff_y[0][node.interBackGhost]<<endl<<flush;
//      cout <<Rhoy(Rho[0])<<'\t'<<DpDy(P[0])<<endl<<flush;
//      cout <<"p "<<P[0][node.center]<<'\t'<<P[0][node.back]<<endl<<flush;
//      cout<<"dy "<<hYcen[node.center]<<'\t'<<hYcen[node.back]<<'\t'
//  	<<endl<<flush;
    }
  //assumes that DRhoX_center and DRhoX_right have been set
  //mwf now do derivatives of pressure differences
  inline virtual void setDPdiff_x(const Vec* Rho, const Vec* DRho, const Vec* P)
    {                  

      DpDiff_x_center[0][node.interRightGhost] = 
	2.0*hYcen[node.center]/(hXcen[node.center]+hXcen[node.right])
	+ hYcen[node.center]*DRhoX_center[node.interRightGhost]
	*gravX[node.interRightGhost];

      DpDiff_x_right[0][node.interRightGhost] = 
	-2.0*hYcen[node.center]/(hXcen[node.center]+hXcen[node.right])
	+ hYcen[node.center]*DRhoX_right[node.interRightGhost]
	*gravX[node.interRightGhost];
    }

  //assumes that DRhoY_center and DRhoY_back have been set
  inline virtual void setDPdiff_y(const Vec* Rho, const Vec* DRho,  const Vec* P)
    {
      
      DpDiff_y_center[0][node.interBackGhost] = 
	2.0*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back])
	+ hXcen[node.center]*DRhoY_center[node.interBackGhost]
	*gravY[node.interBackGhost];

      DpDiff_y_back[0][node.interBackGhost] = 
	-2.0*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back])
	+ hXcen[node.center]*DRhoY_back[node.interBackGhost]
	*gravY[node.interBackGhost];
      
    }


  //now do interface terms
  //assumes that pDiff_x has been set for this interface
  inline void setKrX(const Vec* Kr, const Vec* P)
    { 
#ifdef USE_ARITHMETIC_AVG_ECDM
      KrX[node.interRightGhost] = 
	avgA(Kr[0][node.center],Kr[0][node.right]);
#else
      //mwf have to use whole driving force to do upwinding
      real driveForce = getSinglePhaseFlux_x();
      KrX[node.interRightGhost] = 
	avgU(Kr[0][node.center],Kr[0][node.right],driveForce);
#endif
    }
  //assumes that pDiff_x has been set for this interface
  inline void setKrY(const Vec* Kr, const Vec* P)
    {                  
#ifdef USE_ARITHMETIC_AVG_ECDM
      KrY[node.interBackGhost] = 
	avgA(Kr[0][node.center],Kr[0][node.back]);
#else
      //mwf have to use whole driving force to do upwinding
      real driveForce = getSinglePhaseFlux_y();
      KrY[node.interBackGhost] = 
	avgU(Kr[0][node.center],Kr[0][node.back],
	     driveForce);
#endif
    }

  inline void setDKrX(const VecVecVec& DKr, const Vec* P)
    {                  
#ifdef USE_ARITHMETIC_AVG_ECDM
      DKrX_center[node.interRightGhost] = 
	avgA(DKr[0][0][node.center],0.0);
      
      DKrX_right[node.interRightGhost] = 
	avgA(0.0,DKr[0][0][node.right]);
#else
      real driveForce = getSinglePhaseFlux_x();
      DKrX_center[node.interRightGhost] = 
	avgU(DKr[0][0][node.center],0.0,driveForce);
      
      DKrX_right[node.interRightGhost] = 
	avgU(0.0,DKr[0][0][node.right],driveForce);
#endif
      
    }
  
  inline void setDKrY(const VecVecVec& DKr, const Vec* P)
    {                  
#ifdef USE_ARITHMETIC_AVG_ECDM
      DKrY_center[node.interBackGhost] = 
	avgA(DKr[0][0][node.center],0.0);
      
      DKrY_back[node.interBackGhost] = 
	avgA(0.0,DKr[0][0][node.back]);
#else
      real driveForce = getSinglePhaseFlux_y();
      DKrY_center[node.interBackGhost] = 
	avgU(DKr[0][0][node.center],0.0,driveForce);
      
      DKrY_back[node.interBackGhost] = 
	avgU(0.0,DKr[0][0][node.back],driveForce);
#endif
      
    }
  //try and set up stupid interface terms 
  inline void setDKrXinterface(const VecVecVec& DKr)
    {                  
      DKrX_center[node.interRightGhost] = 0.0;
      DKrX_right[node.interLeftGhost]   = 0.0;

#ifdef USE_ARITHMETIC_AVG_ECDM
      DKrX_right[node.interRightGhost] = 
	avgA(0.0,DKr[0][0][node.right]);

      DKrX_center[node.interLeftGhost] = 
	avgA(DKr[0][0][node.left],0.0);
#else
      //mwf have to use whole driving force to do upwinding
      real driveForce = getSinglePhaseFlux_x();
      DKrX_right[node.interRightGhost] = 
	avgU(0.0,DKr[0][0][node.right],driveForce);

      DKrX_center[node.interLeftGhost] = 
	avgU(DKr[0][0][node.left],0.0,driveForce);
#endif
    }

  //try and set up stupid interface terms 
  inline void setDKrYinterface(const VecVecVec& DKr)
    {                  
      DKrY_center[node.interBackGhost] = 0.0;
      DKrY_back[node.interFrontGhost]   = 0.0;

#ifdef USE_ARITHMETIC_AVG_ECDM
      DKrY_back[node.interBackGhost] = 
	avgA(0.0,DKr[0][0][node.back]);

      DKrY_center[node.interFrontGhost] = 
	avgA(DKr[0][0][node.front],0.0);
#else
      real driveForce = getSinglePhaseFlux_y();
      DKrY_back[node.interBackGhost] = 
	avgU(0.0,DKr[0][0][node.back],driveForce);

      DKrY_center[node.interFrontGhost] = 
	avgU(DKr[0][0][node.front],0.0,driveForce);
#endif
    }

  //assumes that pDiff_x has been set for this interface
  inline void setRhoX(const Vec* Rho, const Vec* P)
    {                  
      RhoX[node.interRightGhost] = 
	avgA(Rho[0][node.center],Rho[0][node.right]);
      
    }
  inline void setRhoY(const Vec* Rho, const Vec* P)
    {                  
      RhoY[node.interBackGhost] = 
	avgA(Rho[0][node.center],Rho[0][node.back]);

    }


  inline void setDRhoX(const Vec* DRho, const Vec* P)
    {                  
      DRhoX_center[node.interRightGhost] = 
	avgA(DRho[0][node.center],0.0);

      DRhoX_right[node.interRightGhost] = 
	avgA(0.0,DRho[0][node.right]);

    }
  
  inline void setDRhoY(const Vec* DRho, const Vec* P)
    {                  
      DRhoY_center[node.interBackGhost] = 
	avgA(DRho[0][node.center],0.0);

      DRhoY_back[node.interBackGhost] = 
	avgA(0.0,DRho[0][node.back]);

    }

  inline void setDRhoXinterface(const Vec* DRho)
    {                  
      //always zero lagrange face term
      DRhoX_center[node.interRightGhost] = 0.0;
      DRhoX_right[node.interLeftGhost]   = 0.0;

      DRhoX_right[node.interRightGhost] = 
	avgA(0.0,DRho[0][node.right]);
      DRhoX_center[node.interLeftGhost] = 
	avgA(DRho[0][node.left],0.0);
    }

  inline void setDRhoYinterface(const Vec* DRho)
    {                  
      //always zero lagrange face term
      DRhoY_center[node.interBackGhost] = 0.0;
      DRhoY_back[node.interFrontGhost]  = 0.0;

      DRhoY_back[node.interBackGhost] = 
	avgA(0.0,DRho[0][node.back]);
      DRhoY_center[node.interFrontGhost] = 
	avgA(DRho[0][node.front],0.0);
    }
  
  //////////////////////////////////////////
  //get flux deriv for extra conductivity term
  inline virtual real getInteriorFluxDerivTerm_x(int face)
  {   
    real fluxDeriv(0.0),gradVal(0.0),coefDeriv(0.0);

    //this is like normal flux calc on right side
    //except derivative term is wrt left variable
    if (face == 1)
      {
	//now get corner terms at x+1/2
	real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
				(*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
    
	real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
	real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];
    
	real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
	
	gradVal = Ki11Avg*pDiff_x[0][node.interRightGhost];

	gradVal += hX12inv
	  *(Ki12front*(pDiff_y[0][node.interFrontGhost]
		       +pDiff_y[0][node.interRightFrontGhost])
	    +Ki12back*(pDiff_y[0][node.interBackGhost]
		       +pDiff_y[0][node.interRightBackGhost]));
	
	//at interface KrX[interRightGhost]=KrX[interLeftGhost]
	coefDeriv = RhoX[node.interLeftGhost]*DKrX_center[node.interLeftGhost]
	  +DRhoX_center[node.interLeftGhost]*KrX[node.interLeftGhost];

	fluxDeriv = gradVal*coefDeriv;

      }
    else
      {
	//this is a gradient on the left side of the interface
	//except derivative term is wrt right variable
	assert(face == -1);
	//now get corner terms at x+1/2
	real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interLeftGhost],
				(*KWsMatCorner)[XX][PLUS][node.interLeftGhost]);
    
	real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interLeftGhost];
	real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interLeftGhost];
    
	real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.left] + hXcen[node.center]);
	
	gradVal = Ki11Avg*pDiff_x[0][node.interLeftGhost];

	gradVal += hX12inv
	  *(Ki12front*(pDiff_y[0][node.interLeftFrontGhost]
		       +pDiff_y[0][node.interFrontGhost])
	    +Ki12back*(pDiff_y[0][node.interLeftBackGhost]
		       +pDiff_y[0][node.interBackGhost]));
	
	//at interface KrX[interRightGhost]=KrX[interLeftGhost]
	coefDeriv = RhoX[node.interLeftGhost]*DKrX_right[node.interLeftGhost]
	  +DRhoX_right[node.interLeftGhost]*KrX[node.interLeftGhost];

	fluxDeriv = gradVal*coefDeriv;
      }//end else
    return fluxDeriv;
  }

  inline virtual real getInteriorFluxDerivTerm_y(int face)
  {   
    real fluxDeriv(0.0),gradVal(0.0),coefDeriv(0.0);

    //this is like normal flux calc on back side
    //except derivative term is wrt left variable
    if (face == 1)
      {
	//now get corner terms at y+/12
	real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
			     (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
	real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
	real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
	
	real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);
 
	gradVal = Ki22Avg*pDiff_y[0][node.interBackGhost];

	gradVal += hY12inv
	  *(Ki12left*(pDiff_x[0][node.interLeftGhost]
		      +pDiff_x[0][node.interBackLeftGhost])
	    +Ki12right*(pDiff_x[0][node.interRightGhost]
			+pDiff_x[0][node.interBackRightGhost]));


	//at interface KrY[interBackGhost]=KrY[interFrontGhost]
	coefDeriv = RhoY[node.interFrontGhost]*DKrY_center[node.interFrontGhost]
	  +DRhoY_center[node.interFrontGhost]*KrY[node.interFrontGhost];

	fluxDeriv = gradVal*coefDeriv;
	
      }
    else
      {
	//this is a gradient on the front side of the interface
	assert(face == -1);
	//now get corner terms at y+/12
	real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interFrontGhost],
			     (*KWsMatCorner)[YY][PLUS][node.interFrontGhost]);
	real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interFrontGhost];
	real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interFrontGhost];
	
	real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.front] + hYcen[node.center]);
 
	gradVal = Ki22Avg*pDiff_y[0][node.interFrontGhost];

	gradVal += hY12inv
	  *(Ki12left*(pDiff_x[0][node.interFrontLeftGhost]
		      +pDiff_x[0][node.interLeftGhost])
	    +Ki12right*(pDiff_x[0][node.interFrontRightGhost]
			+pDiff_x[0][node.interRightGhost]));


	//at interface KrY[interBackGhost]=KrY[interFrontGhost]
	coefDeriv = RhoY[node.interBackGhost]*DKrY_back[node.interBackGhost]
	  +DRhoY_back[node.interBackGhost]*KrY[node.interBackGhost];

	fluxDeriv = gradVal*coefDeriv;
		
      }//end else
    return fluxDeriv;
  }

  //////////////////////////////////////////
  //use this for upwinding
  inline virtual real getSinglePhaseFlux_x()
    {
      real fluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
      
      fluxVal = Ki11Avg*pDiff_x[0][node.interRightGhost];

      fluxVal += hX12inv
	*(Ki12front*(pDiff_y[0][node.interFrontGhost]
		     +pDiff_y[0][node.interRightFrontGhost])
	  +Ki12back*(pDiff_y[0][node.interBackGhost]
		     +pDiff_y[0][node.interRightBackGhost]));
      
      return fluxVal;
    }
  inline virtual real getSinglePhaseFlux_y()
    {
      real fluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);
      
      fluxVal  =Ki22Avg*pDiff_y[0][node.interBackGhost];

      fluxVal += hY12inv
	*(Ki12left*(pDiff_x[0][node.interLeftGhost]
		      +pDiff_x[0][node.interBackLeftGhost])
	  +Ki12right*(pDiff_x[0][node.interRightGhost]
		      +pDiff_x[0][node.interBackRightGhost]));

      return fluxVal;
    }
  inline virtual void setFluxNew_x(const VecVec* K, const Vec* Kr,
				   const Vec* Rho, const Vec* P)
  {                  
    for (vi=0;vi<nv;vi++)
      {
	real rhoAvgKrAvg =  RhoX[node.interRightGhost]*
	  KrX[node.interRightGhost];
	//now get corner terms at x+1/2
	real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
				(*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
	
	real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
	real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];

	//mwf put in correct off diagonal scaling terms
	real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
 
	flux_x[vi][node.interRight] =
	  rhoAvgKrAvg*Ki11Avg*pDiff_x[vi][node.interRightGhost];

	flux_x[vi][node.interRight] += hX12inv*rhoAvgKrAvg
	  *(Ki12front*(pDiff_y[vi][node.interFrontGhost]
		       +pDiff_y[vi][node.interRightFrontGhost])
	    +Ki12back*(pDiff_y[vi][node.interBackGhost]
		       +pDiff_y[vi][node.interRightBackGhost]));

//  	//mwf add for debugging
//  	if (fabs(flux_x[vi][node.interRight]) > 1.0e4)
//  	  {
//  	    cerr<<"in setFluxNew_x, flux_x[vi]["<<node.interRight
//  		<<"]= "<<flux_x[vi][node.interRight]<<endl;
//  	    cerr<<"\t rhoAvgKrAvg= "<<rhoAvgKrAvg<<endl;
//  	    cerr<<"\t P["<< node.center<<"]= "<<P[vi][node.center]<<endl;
//  	    cerr<<"\t pDiff_y[node.interBackGhost]= "
//  		<<pDiff_y[vi][node.interBackGhost]<<endl;
//  	    cerr<<"\t pDiff_y[node.interFrontGhost]= "
//  		<<pDiff_y[vi][node.interFrontGhost]<<endl;
//  	    cerr<<"\t pDiff_x[node.interRightGhost]= "
//  		<<pDiff_x[vi][node.interRightGhost]<<endl;
//  	  }
      }//end vi loop
  }


  inline virtual void setFluxNew_y(const VecVec* K, const Vec* Kr,
				   const Vec* Rho, const Vec* P)
  {                  
    for (vi=0;vi<nv;vi++)
      {
	real rhoAvgKrAvg =  RhoY[node.interBackGhost]*
	  KrY[node.interBackGhost];
	//now get corner terms at y+/12
	real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
			     (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
	real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
	real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
	
	//mwf put in correct off diagonal scaling term
	real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);
 
	flux_y[vi][node.interBack] =
	  rhoAvgKrAvg*Ki22Avg*pDiff_y[vi][node.interBackGhost];

//  	std::cerr<<"in DifKgradECDM:: flux_y["<<node.interBack<<"]= "
//  		 <<flux_y[vi][node.interBack]
//  		 <<" pDiff_y[vi]["<<node.interBackGhost<<"]= "
//  		 <<pDiff_y[vi][node.interBackGhost]<<endl;
//  	std::cerr<<" rhoAvgKrAvg= "<<rhoAvgKrAvg<<" Ki22Avg = "
//  		 <<Ki22Avg<<endl;
	  
	flux_y[vi][node.interBack] += hY12inv*rhoAvgKrAvg
	  *(Ki12left*(pDiff_x[vi][node.interLeftGhost]
		      +pDiff_x[vi][node.interBackLeftGhost])
	    +Ki12right*(pDiff_x[vi][node.interRightGhost]
			+pDiff_x[vi][node.interBackRightGhost]));

      }
  }

  //mwf add new functions for doing analytical jacobians
  inline virtual void setDFluxNew_x(const VecVec* K, const Vec* Kr,
				    const Vec* Rho, const Vec* P, 
				    const VecVecVec& DKr, const Vec* DRho)
    {

      for (vi=0;vi<nv;vi++)
	{
	  for (vj=0;vj<nv;vj++)
	    {
	      //interior
	      //mwf from flux calculation
	      real rhoAvgKrAvg=  RhoX[node.interRightGhost]
		*KrX[node.interRightGhost];
	      real Ki11Avg  = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
				   (*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
	      
	      real Ki12front= (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
	      real Ki12back = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];


		//mwf put in correct off diagonal scaling terms
	      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
	      //mwf derivative terms
	      real DRhoAvgKrAvg_center= DRhoX_center[node.interRightGhost]
		*KrX[node.interRightGhost];
	      real DRhoAvgKrAvg_right = DRhoX_right[node.interRightGhost]
		*KrX[node.interRightGhost];
	      
	      real rhoAvgDKrAvg_center= RhoX[node.interRightGhost] 
		*DKrX_center[node.interRightGhost];

	      real rhoAvgDKrAvg_right =  RhoX[node.interRightGhost] 
		*DKrX_right[node.interRightGhost];
	      
	      real transTerm = hX12inv
		*(Ki12front*(pDiff_y[vi][node.interFrontGhost]
			     +pDiff_y[vi][node.interRightFrontGhost])
		  +Ki12back*(pDiff_y[vi][node.interBackGhost]
			     +pDiff_y[vi][node.interRightBackGhost]));

	      //do derivatives of coefficients first 
	      //mwf coefficients only a function of center and right (for now)
	      //mwf relative perm can be a function of the other phase
	      Dflux_x_center[vi][vj][node.interRight] = rhoAvgDKrAvg_center*
		Ki11Avg*pDiff_x[vi][node.interRightGhost];
	      
	      Dflux_x_center[vi][vj][node.interRight] += rhoAvgDKrAvg_center*transTerm;
	      
	      Dflux_x_right[vi][vj][node.interRight] = rhoAvgDKrAvg_right*
		Ki11Avg*pDiff_x[vi][node.interRightGhost];
	      
	      Dflux_x_right[vi][vj][node.interRight] += rhoAvgDKrAvg_right*transTerm;

	      if (vi == vj)
		{ 
		  //take care of DRho terms first
		  Dflux_x_center[vi][vj][node.interRight] += DRhoAvgKrAvg_center*
		    Ki11Avg*pDiff_x[vi][node.interRightGhost];
		  
		  Dflux_x_center[vi][vj][node.interRight] += 
		    DRhoAvgKrAvg_center*transTerm;

		  Dflux_x_right[vi][vj][node.interRight] += DRhoAvgKrAvg_right*
		    Ki11Avg*pDiff_x[vi][node.interRightGhost];
		  
		  Dflux_x_right[vi][vj][node.interRight] += 
		    DRhoAvgKrAvg_right*transTerm;
		  
		  //now start getting derivative terms from appearance of 
		  //actual pressure
		  //finish out center and right
		  Dflux_x_center[vi][vj][node.interRight] += rhoAvgKrAvg*
		    Ki11Avg*DpDiff_x_center[vi][node.interRightGhost];
	
		  Dflux_x_center[vi][vj][node.interRight] += 
		    hX12inv*rhoAvgKrAvg
		    *(Ki12back*DpDiff_y_center[vi][node.interBackGhost]
		      +Ki12front*DpDiff_y_back[vi][node.interFrontGhost]);
		  
		  Dflux_x_right[vi][vj][node.interRight] += rhoAvgKrAvg*
		    Ki11Avg*DpDiff_x_right[vi][node.interRightGhost];
		  
		  Dflux_x_right[vi][vj][node.interRight] += 
		    hX12inv*rhoAvgKrAvg
		    *(Ki12back*DpDiff_y_center[vi][node.interRightBackGhost]
		      +Ki12front*DpDiff_y_back[vi][node.interRightFrontGhost]);
		  
		  //now do up and down terms (i,j+1)
		  Dflux_x_back[vi][vj][node.interRight] = 
		    hX12inv*rhoAvgKrAvg
		    *Ki12back*DpDiff_y_back[vi][node.interBackGhost];
		  //now do up and down terms (i,j-1)
		  Dflux_x_front[vi][vj][node.interRight] = 
		    hX12inv*rhoAvgKrAvg
		    *Ki12front*DpDiff_y_center[vi][node.interFrontGhost];
		  
		  //now do corner terms (i+1,j+1)
		  Dflux_x_rightBack[vi][vj][node.interRight] = 
		    hX12inv*rhoAvgKrAvg
		    *Ki12back*DpDiff_y_back[vi][node.interRightBackGhost];
		  //now do corner terms (i+1,j-1)
		  Dflux_x_rightFront[vi][vj][node.interRight] = 
		    hX12inv*rhoAvgKrAvg
		    *Ki12front*DpDiff_y_center[vi][node.interRightFrontGhost];
		  
		}//end if (vi==vj)
//  	      //mwf add for debugging
//  	      if (abs(Dflux_x_center[vi][vi][node.interRight]) < 1.0e-18)
//  		{
//  		  cerr<<"in setDFluxNew_x, Dflux_x_center[vi][vi]["<<node.interRight
//  		      <<"]= "<<Dflux_x_center[vi][vi][node.interRight]<<endl;
//  		  cerr<<"\t rhoAvgKrAvg= "<<rhoAvgKrAvg<<endl;
//  		  cerr<<"\t P["<< node.center<<"]= "<<P[vi][node.center]<<endl;
//  		  cerr<<"\t DpDiff_y_back[node.interRightBackGhost]= "
//  		      <<DpDiff_y_back[vi][node.interRightBackGhost]<<endl;
//  		  cerr<<"\t DpDiff_y_center[node.interRightFrontGhost]= "
//  		      <<DpDiff_y_center[vi][node.interRightFrontGhost]<<endl;
//  		  cerr<<"\t Dflux_x_rightBack[vi][vj][node.interRight]  = "
//  		      <<Dflux_x_rightBack[vi][vj][node.interRight] <<endl;
//  		}

	    }//end vj				 
	}//end vi
    }//end function

  //mwf add new functions for doing analytical jacobians
  inline virtual void setDFluxNew_y(const VecVec* K, const Vec* Kr,
				    const Vec* Rho, const Vec* P, 
				    const VecVecVec& DKr, const Vec* DRho)
    {

      for (vi=0;vi<nv;vi++)
	{
	  for (vj=0;vj<nv;vj++)
	    {
	      //these terms are from flux calculation
	      real rhoAvgKrAvg =  RhoY[node.interBackGhost]
		*KrY[node.interBackGhost];

	      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
				   (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
	      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
	      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
	
	      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);
	      
	      //now get terms for derivatve of coefficients
	      real DRhoAvgKrAvg_center= DRhoY_center[node.interBackGhost]
		*KrY[node.interBackGhost];
	      real DRhoAvgKrAvg_back  = DRhoY_back[node.interBackGhost]
		*KrY[node.interBackGhost];
	      
	      real rhoAvgDKrAvg_center=  
		RhoY[node.interBackGhost]*DKrY_center[node.interBackGhost];
	      real rhoAvgDKrAvg_back =  
		RhoY[node.interBackGhost]*DKrY_back[node.interBackGhost];
	      
	      real transTerm = hY12inv*
		(Ki12left*(pDiff_x[vi][node.interLeftGhost]
			   +pDiff_x[vi][node.interBackLeftGhost])
		 +Ki12right*(pDiff_x[vi][node.interRightGhost]
			     +pDiff_x[vi][node.interBackRightGhost]));

	      //do derivatives of coefficients first 
	      //mwf coefficients only a function of center and back (for now)
	      //mwf relative perm can be a function of the other phase
	      
	      Dflux_y_center[vi][vj][node.interBack] =
		rhoAvgDKrAvg_center*Ki22Avg*pDiff_y[vi][node.interBackGhost];
	      
	      Dflux_y_center[vi][vj][node.interBack] += rhoAvgDKrAvg_center*transTerm;
	      
	      Dflux_y_back[vi][vj][node.interBack] =
		rhoAvgDKrAvg_back*Ki22Avg*pDiff_y[vi][node.interBackGhost];
	      
	      Dflux_y_back[vi][vj][node.interBack] += rhoAvgDKrAvg_back*transTerm;
//  	      if (node.center==26)
//  		{
//  		  cerr<<"in DivKgradECDM setDFlux_y["<<node.interBack
//  		      <<"] rhoAvgDKrAvg_center= "<<rhoAvgDKrAvg_center<<endl;
//  		  cerr<<"   DRhoAvgKrAvg_center= "<<DRhoAvgKrAvg_center<<endl;
//  		  cerr<<"   Kr["<<node.center<<"]= "<<Kr[vi][node.center]<<endl;
//  		  cerr<<"   Kr["<<node.back<<"]= "<<Kr[vi][node.back]<<endl;
//  		  cerr<<"   Ki22Avg= "<<Ki22Avg<<" pDiff_y["<<node.interBackGhost
//  		      <<"]= "<<pDiff_y[vi][node.interBackGhost]<<endl;
//  		  cerr<<"  Ki12 = "<<Ki12<<" Ki12back= "<<Ki12back<<endl;
//  		}
	  //now go on to same phase for rest of terms
	      if (vi == vj)
		{
		  //do density terms
		  Dflux_y_center[vi][vj][node.interBack] +=
		    DRhoAvgKrAvg_center*Ki22Avg*pDiff_y[vi][node.interBackGhost];
		  
		  Dflux_y_center[vi][vj][node.interBack] += 
		    DRhoAvgKrAvg_center*transTerm;
		  
		  Dflux_y_right[vi][vj][node.interBack] +=
		    DRhoAvgKrAvg_back*Ki22Avg*pDiff_y[vi][node.interBackGhost];
		  
		  Dflux_y_right[vi][vj][node.interBack] += 
		    DRhoAvgKrAvg_back*transTerm;
		  
		  //now go on to derivatives for terms with pressure explicit
		  //finish off center and back first
		  Dflux_y_center[vi][vj][node.interBack] += rhoAvgKrAvg*
		    Ki22Avg*DpDiff_y_center[vi][node.interBackGhost];
		  
		  Dflux_y_center[vi][vj][node.interBack] += 
		    hY12inv*rhoAvgKrAvg*
		    (Ki12right*DpDiff_x_center[vi][node.interRightGhost]
		     +Ki12left*DpDiff_x_right[vi][node.interLeftGhost]);
		  
		  Dflux_y_back[vi][vj][node.interBack] += rhoAvgKrAvg
		    *Ki22Avg*DpDiff_y_back[vi][node.interBackGhost];
		  
		  Dflux_y_back[vi][vj][node.interBack] += 
		    hY12inv*rhoAvgKrAvg*
		    (Ki12right*DpDiff_x_center[vi][node.interBackRightGhost]
		     +Ki12left*DpDiff_x_right[vi][node.interBackLeftGhost]);
		  
		  //now do left and right terms
		  //(i+1,j)
		  Dflux_y_right[vi][vj][node.interBack]= 
		    hY12inv*rhoAvgKrAvg
		    *Ki12right*DpDiff_x_right[vi][node.interRightGhost];
		  
		  //(i-1,j)
		  Dflux_y_left[vi][vj][node.interBack]= 
		    hY12inv*rhoAvgKrAvg
		    *Ki12left*DpDiff_x_center[vi][node.interLeftGhost];

		  //now do corner terms
		  //(i+1,j+1)
		  Dflux_y_backRight[vi][vj][node.interBack]= 
		    hY12inv*rhoAvgKrAvg
		    *Ki12right*DpDiff_x_right[vi][node.interBackRightGhost];
		    
		  //i-1,j+1
		  Dflux_y_backLeft[vi][vj][node.interBack]= 
		    hY12inv*rhoAvgKrAvg
		    *Ki12left*DpDiff_x_center[vi][node.interBackLeftGhost];
		  

		}//end if vi==vj
//  	      //mwf add for debugging
//  	      if (abs(Dflux_y_center[vi][vi][node.interBack]) < 1.0e-18)
//  		{
//  		  cerr<<"in setDFluxNew_y, Dflux_y_center[vi][vi]["<<node.interBack
//  		      <<"]= "<<Dflux_y_center[vi][vi][node.interBack]<<endl;
//  		  cerr<<"\t rhoAvgKrAvg= "<<rhoAvgKrAvg<<endl;
//  		  cerr<<"\t P["<< node.center<<"]= "<<P[vi][node.center]<<endl;
//  		  cerr<<"\t DpDiff_x_right[node.interBackRightGhost]= "
//  		      <<DpDiff_x_right[vi][node.interBackRightGhost]<<endl;
//  		  cerr<<"\t DpDiff_x_center[node.interBackLeftGhost]= "
//  		      <<DpDiff_x_center[vi][node.interBackLeftGhost]<<endl;
//  		  cerr<<"\t Dflux_y_backRight[vi][vj][node.interBack]  = "
//  		      <<Dflux_y_backRight[vi][vj][node.interBack] <<endl;
//  		}
	    }//end loop for vj
	}//end loop for vi
    }//end function
  ////////
  //====== gradient derivative terms
  //try and use these for no flow or gradient conditions 
  inline virtual real getSinglePhaseFlux_xRightBnd()
    {
      //have to use values at one cell in?
      real fluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interLeftGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interLeftGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interLeftGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interLeftGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.left]);
      
      fluxVal = Ki11Avg*pDiff_x[0][node.interLeftGhost];

      fluxVal += hX12inv
	*(Ki12front*(pDiff_y[0][node.interLeftFrontGhost]
		     +pDiff_y[0][node.interFrontGhost])
	  +Ki12back*(pDiff_y[0][node.interLeftBackGhost]
		     +pDiff_y[0][node.interBackGhost]));
      
      return fluxVal;
    }
  //same as right flux for left boundary?
  inline virtual real getSinglePhaseFlux_xLeftBnd()
    {
      real fluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
      
      fluxVal = Ki11Avg*pDiff_x[0][node.interRightGhost];

      fluxVal += hX12inv
	*(Ki12front*(pDiff_y[0][node.interFrontGhost]
		     +pDiff_y[0][node.interRightFrontGhost])
	  +Ki12back*(pDiff_y[0][node.interBackGhost]
		     +pDiff_y[0][node.interRightBackGhost]));
      
      return fluxVal;
    }
  inline virtual real getSinglePhaseFlux_yBackBnd()
    {
      real fluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interFrontGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interFrontGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interFrontGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interFrontGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.front]);
      
      fluxVal  =Ki22Avg*pDiff_y[0][node.interFrontGhost];

      fluxVal += hY12inv
	*(Ki12left*(pDiff_x[0][node.interFrontLeftGhost]
		      +pDiff_x[0][node.interLeftGhost])
	  +Ki12right*(pDiff_x[0][node.interFrontRightGhost]
		      +pDiff_x[0][node.interRightGhost]));

      return fluxVal;
    }

  inline virtual real getSinglePhaseFlux_yFrontBnd()
    {
      real fluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);
      
      fluxVal  =Ki22Avg*pDiff_y[0][node.interBackGhost];

      fluxVal += hY12inv
	*(Ki12left*(pDiff_x[0][node.interLeftGhost]
		      +pDiff_x[0][node.interBackLeftGhost])
	  +Ki12right*(pDiff_x[0][node.interRightGhost]
		      +pDiff_x[0][node.interBackRightGhost]));

      return fluxVal;
    }

  //mwf now use this for derivatives of no flow condition?
  //======xRightBnd derivatives=====
  //try and use these for no flow or gradient conditions 
  inline virtual real getDSinglePhaseFlux_xRightBnd_center()
    {
      //have to use values at one cell in?
      real DfluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interLeftGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interLeftGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interLeftGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interLeftGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.left]);
      
      DfluxVal = Ki11Avg*DpDiff_x_right[0][node.interLeftGhost];

      //these should be zero
      DfluxVal += hX12inv
	*(Ki12front*DpDiff_y_back[0][node.interFrontGhost]
	  +Ki12back*DpDiff_y_center[0][node.interBackGhost]);
      
      return DfluxVal;
    }
  //try and use these for no flow or gradient conditions 
  inline virtual real getDSinglePhaseFlux_xRightBnd_left()
    {
      //have to use values at one cell in?
      real DfluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interLeftGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interLeftGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interLeftGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interLeftGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.left]);
      
      DfluxVal = Ki11Avg*DpDiff_x_center[0][node.interLeftGhost];

      DfluxVal += hX12inv
	*(Ki12front*DpDiff_y_back[0][node.interLeftFrontGhost]
	  +Ki12back*DpDiff_y_center[0][node.interLeftBackGhost]);
      
      return DfluxVal;
    }
  //try and use these for no flow or gradient conditions 
  inline virtual real getDSinglePhaseFlux_xRightBnd_leftBack()
    {
      //have to use values at one cell in?
      real DfluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interLeftGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interLeftGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interLeftGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interLeftGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.left]);
      
      DfluxVal = hX12inv*Ki12back
	*DpDiff_y_back[0][node.interLeftBackGhost];

      
      return DfluxVal;
    }
   //try and use these for no flow or gradient conditions 
  inline virtual real getDSinglePhaseFlux_xRightBnd_leftFront()
    {
      //have to use values at one cell in?
      real DfluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interLeftGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interLeftGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interLeftGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interLeftGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.left]);
      
      DfluxVal = hX12inv*Ki12front
	*DpDiff_y_center[0][node.interLeftFrontGhost];

      return DfluxVal;
    }
  //====== xLeftBnd derivatives =====
  //same as right flux for left boundary?
  inline virtual real getDSinglePhaseFlux_xLeftBnd_center()
    {
      real DfluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
      
      DfluxVal = Ki11Avg*DpDiff_x_center[0][node.interRightGhost];

      //this should be zero
      DfluxVal += hX12inv
	*(Ki12front*DpDiff_y_back[0][node.interFrontGhost]
	  +Ki12back*DpDiff_y_center[0][node.interBackGhost]);
      
      return DfluxVal;
    }
  inline virtual real getDSinglePhaseFlux_xLeftBnd_right()
    {
      real DfluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
      
      DfluxVal = Ki11Avg*DpDiff_x_right[0][node.interRightGhost];

      DfluxVal += hX12inv
	*(Ki12front*DpDiff_y_back[0][node.interRightFrontGhost]
	  +Ki12back*DpDiff_y_center[0][node.interRightBackGhost]);
      
      return DfluxVal;
    }
  inline virtual real getDSinglePhaseFlux_xLeftBnd_rightBack()
    {
      real DfluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
      
      DfluxVal = hX12inv*Ki12back
	*DpDiff_y_back[0][node.interRightBackGhost];
      
      return DfluxVal;
    }
   //same as right flux for left boundary?
  inline virtual real getDSinglePhaseFlux_xLeftBnd_rightFront()
    {  
      real DfluxVal(0.0);
      //now get corner terms at x+1/2
      real Ki11Avg     = avgA((*KWsMatCorner)[XX][MINUS][node.interRightGhost],
			      (*KWsMatCorner)[XX][PLUS][node.interRightGhost]);
      
      real Ki12front = (*KWsMatCorner)[XY][MINUS][node.interRightGhost];
      real Ki12back  = (*KWsMatCorner)[XY][PLUS][node.interRightGhost];
      
      real hX12inv     =  0.5*hYcen[node.center]/(hXcen[node.center] + hXcen[node.right]);
      

      DfluxVal = hX12inv*Ki12front
	*DpDiff_y_center[0][node.interRightFrontGhost];
      
      return DfluxVal;
    }
  //============ yBackBnd derivatives
  inline virtual real getDSinglePhaseFlux_yBackBnd_center()
    {
      real DfluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interFrontGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interFrontGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interFrontGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interFrontGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.front]);
      
      DfluxVal  =Ki22Avg*DpDiff_y_back[0][node.interFrontGhost];

      //mwf these should be zero
      DfluxVal += hY12inv
	*(Ki12left*DpDiff_x_right[0][node.interLeftGhost]
	  +Ki12right*DpDiff_x_center[0][node.interRightGhost]);

      return DfluxVal;
    }
  inline virtual real getDSinglePhaseFlux_yBackBnd_front()
    {
      real DfluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interFrontGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interFrontGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interFrontGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interFrontGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.front]);
      
      DfluxVal  =Ki22Avg*DpDiff_y_center[0][node.interFrontGhost];
      
      DfluxVal += hY12inv
	*(Ki12left*DpDiff_x_right[0][node.interFrontLeftGhost]
	  +Ki12right*DpDiff_x_center[0][node.interFrontRightGhost]);
		

      return DfluxVal;
    }

  inline virtual real getDSinglePhaseFlux_yBackBnd_frontRight()
    {
      real DfluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interFrontGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interFrontGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interFrontGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interFrontGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.front]);

      DfluxVal = hY12inv*Ki12right
	*DpDiff_x_right[0][node.interFrontRightGhost];

      return DfluxVal;
    }

  inline virtual real getDSinglePhaseFlux_yBackBnd_frontLeft()
    {
      real DfluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interFrontGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interFrontGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interFrontGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interFrontGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.front]);

      DfluxVal = hY12inv*Ki12left
	*DpDiff_x_center[0][node.interFrontLeftGhost];

      return DfluxVal;
    }
  //============ yFrontBnd derivatives
  inline virtual real getDSinglePhaseFlux_yFrontBnd_center()
    {
      real DfluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);
      
      DfluxVal  =Ki22Avg*DpDiff_y_center[0][node.interBackGhost];

      //these should be zero
      DfluxVal += hY12inv
	*(Ki12left*DpDiff_x_right[0][node.interLeftGhost]
	  +Ki12right*DpDiff_x_center[0][node.interRightGhost]);

      return DfluxVal;
    }
  inline virtual real getDSinglePhaseFlux_yFrontBnd_back()
    {
      real DfluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);

      DfluxVal  =Ki22Avg*DpDiff_y_back[0][node.interBackGhost];

      DfluxVal += hY12inv
	*(Ki12left*DpDiff_x_right[0][node.interBackLeftGhost]
	  +Ki12right*DpDiff_x_center[0][node.interBackRightGhost]);

      return DfluxVal;
    }

  inline virtual real getDSinglePhaseFlux_yFrontBnd_backRight()
    {
      real DfluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);

      DfluxVal = hY12inv*Ki12right
	*DpDiff_x_right[0][node.interBackRightGhost];

      return DfluxVal;
    }

  inline virtual real getDSinglePhaseFlux_yFrontBnd_backLeft()
    {
      real DfluxVal(0.0);

      //now get corner terms at y+/12
      real Ki22Avg  = avgA((*KWsMatCorner)[YY][MINUS][node.interBackGhost],
			   (*KWsMatCorner)[YY][PLUS][node.interBackGhost]);
      real Ki12left =  (*KWsMatCorner)[YX][MINUS][node.interBackGhost];
      real Ki12right=  (*KWsMatCorner)[YX][PLUS][node.interBackGhost];
      
      real hY12inv     =  0.5*hXcen[node.center]/(hYcen[node.center] + hYcen[node.back]);

      DfluxVal = hY12inv*Ki12left
	*DpDiff_x_center[0][node.interBackLeftGhost];

      return DfluxVal;
    }
  //mwf add for getting fluxes out for postProcessing??
  const Vec& getFlux_x(int vi) { return flux_x[vi]; }
  const Vec& getFlux_y(int vi) { return flux_y[vi]; }
  const Vec& getFlux_z(int vi) { return flux_z[vi]; }
  Vec& setFlux_x(int n=0){ return flux_x[vi]; };
  Vec& setFlux_y(int n=0){ return flux_y[vi]; };
  Vec& setFlux_z(int n=0){ return flux_z[vi]; };

  //have Rich class pass these in
  void setSpatialIncrement(int dim,const Vec& hVec);
  const Vec& getSpatialIncrement(int dim);

  //setup terms for irregular grid?
  void setReferenceGravity(int dim, const Vec& refG);
  const Vec& getReferenceGravity(int dim);

  void setIntrinsicPermeability(VecVecVec* kwsInP)
    { assert(kwsInP); KWsMatCorner = kwsInP; }
  VecVecVec* getReferenceGravity()
    { return KWsMatCorner; }

protected:
  int i,j,k,vi,vj,nxNodes,nyNodes,nzNodes,local_nxNodes,local_nyNodes,local_nzNodes;
  real oneOverdx,oneOverdy,oneOverdz,gx,gy,gz;
  Vec flux_x[nv],flux_y[nv],flux_z[nv];
  VecVecVec Dflux_x_center,Dflux_x_right,Dflux_y_center,Dflux_y_back,Dflux_z_center,Dflux_z_top;
  BC* bc;
  StencilMM& node;
  //mwf add for MM
  //spatial increments
  Vec hXcen,hYcen,hZcen;

  //differences in pressure (\Qtilde^{x}_{i+1/2} and \Qtilde^{y}_{ij+1/2})
  Vec pDiff_x[nv],pDiff_y[nv],pDiff_z[nv];
  //gravity terms
  Vec gravX,gravY,gravZ;

  //derivatives of differences in pressure
  Vec DpDiff_x_center[nv],DpDiff_x_right[nv],
    DpDiff_y_center[nv],DpDiff_y_back[nv];
  //have to add extra arrays for flux derivatives wrt corner points
  VecVecVec Dflux_x_back,Dflux_x_front,Dflux_x_rightBack,Dflux_x_rightFront;
  VecVecVec Dflux_y_right,Dflux_y_left,Dflux_y_backRight,Dflux_y_backLeft;

  //try and lump together intrinsic permeabilities here?
  VecVecVec* KWsMatCorner;
  //hold interface rel perms and densities here now?
  //dimensioned like 
  Vec KrX,KrY,RhoX,RhoY;
  //hold the derivatives of these quantities
  Vec DKrX_center,DKrX_right,DKrY_center,DKrY_back;
  Vec DRhoX_center,DRhoX_right,DRhoY_center,DRhoY_back;

};

//have Rich class pass these in
template<class BC, int nv>
void DivKgradECDM<BC,nv>::setSpatialIncrement(int dim,const Vec& hVec)
{ 
  switch (dim)
    {
    case 0:
      {
	hXcen = hVec;
	break;
      }
    case 1:
      {
	hYcen = hVec;
	break;
      }
    case 2:
      {
	hZcen = hVec;
	break;
      }
    default:
      {
	std::cerr<<"incorrect dimension in setSpatial increment"<<std::endl;
	assert(0);
	break;
      }
    }
}

//have Rich class pass these in
template<class BC, int nv>
const Vec& DivKgradECDM<BC,nv>::getSpatialIncrement(int dim)
{ 
  switch (dim)
    {
    case 0:
      {
	return hXcen;
	break;
      }
    case 1:
      {
	return hYcen;
	break;
      }
    case 2:
      {
	return hZcen;
	break;
      }
    default:
      {
	std::cerr<<"incorrect dimension in setSpatial increment"<<std::endl;
	assert(0);
	return hXcen;
	break;
      }
    }
}
//have Rich class pass these in
template<class BC, int nv>
void DivKgradECDM<BC,nv>::setReferenceGravity(int dim,const Vec& refG)
{ 
  switch (dim)
    {
    case 0:
      {
	gravX = refG;
	break;
      }
    case 1:
      {
	gravY = refG;
	break;
      }
    case 2:
      {
	gravZ = refG;
	break;
      }
    default:
      {
	std::cerr<<"incorrect dimension in setSpatial increment"<<std::endl;
	assert(0);
      }
    }
}

//have Rich class pass these in
template<class BC, int nv>
const Vec& DivKgradECDM<BC,nv>::getReferenceGravity(int dim)
{ 
  switch (dim)
    {
    case 0:
      {
	return gravX;
	break;
      }
    case 1:
      {
	return gravY;
	break;
      }
    case 2:
      {
	return gravZ;
	break;
      }
    default:
      {
	std::cerr<<"incorrect dimension in setSpatial increment"<<std::endl;
	assert(0);
      }
    }
}

template<class BC, int nv>
class DivKgradECDM1d : public DivKgradECDM<BC,nv>
{
public:
  DivKgradECDM1d(BC* bcIn, 
                 StencilMM& nodeIn,
                 int nNodes, 
                 real g, 
                 real oneOverd);
  virtual ~DivKgradECDM1d(){}
  inline void setPdiff_x(const Vec* Rho, const Vec* P);
  inline void setDPdiff_x(const Vec* Rho, const Vec* DRho, const Vec* P);
  inline real getInteriorFluxDerivTerm_x(int face);
  inline real getSinglePhaseFlux_x();
  inline void setFluxNew_x(const VecVec* K, const Vec* Kr,
                           const Vec* Rho, const Vec* P);

  inline void setDFluxNew_x(const VecVec* K, const Vec* Kr,
                                  const Vec* Rho, const Vec* P, 
                            const VecVecVec& DKr, const Vec* DRho);
  inline real getSinglePhaseFlux_xRightBnd();
  inline real getSinglePhaseFlux_xLeftBnd();
  inline real getDSinglePhaseFlux_xRightBnd_center();
  inline real getDSinglePhaseFlux_xRightBnd_left();
  inline real getDSinglePhaseFlux_xRightBnd_leftBack();
  inline real getDSinglePhaseFlux_xRightBnd_leftFront();
  inline real getDSinglePhaseFlux_xRightBnd_leftFront();
  inline real getDSinglePhaseFlux_xLeftBnd_center();
  inline real getDSinglePhaseFlux_xLeftBnd_right();
  inline real getDSinglePhaseFlux_xLeftBnd_rightBack();
  inline real getDSinglePhaseFlux_xLeftBnd_rightFront();
  virtual void computeDivergence(const VecVec* K, const Vec* Kr,
                                 const Vec* Rho, 
                                 const Vec* P, 
                                 Vec* Div);
virtual void computeDivergenceJac(const VecVec* K, const Vec* Kr,
                                  const Vec* Rho, 
                                  const Vec* P, 
                                  const VecVecVec& DKr, const Vec* DRho, 
                                  VecVecVecVec& DivJac);
virtual void computeDivergence(const Vec* K, 
                               const Vec* Rho, 
                               const Vec* P, 
                               Vec* Div);
virtual void computeDivergenceJac(const Vec* K,
                                  const Vec* Rho, 
                                  const Vec* P,
                                  const VecVecVec& DK, const Vec* DRho, 
                                  VecVecVecVec& DivJac);
}; 

template<class BC, int nv>
DivKgradECDM1d<BC,nv>::DivKgradECDM1d(BC* bcIn, 
                                      StencilMM& nodeIn,
                                      int nNodes, 
                                      real g, 
                                      real oneOverd):
  DivKgradECDM<BC,nv>(bcIn,nodeIn)
{
  this->nxNodes=nNodes;
  local_nxNodes = nodeIn.local_nxNodes;
  
  gx=g;
  oneOverdx=oneOverd;
  
  //setup gravity term for ghost interfaces
  gravX.newsize(Vec::LOCAL,node.ghost_nxNodes+1);
  gravX = gx;
  //what do I need to do here?
  hXcen.newsize(Vec::LOCAL,node.da_);
  hYcen.newsize(Vec::LOCAL,node.da_);
  hZcen.newsize(Vec::LOCAL,node.da_);

  hXcen = 1.0/oneOverdx;
  hYcen = 1.0;
  hZcen = 1.0;

  Dflux_x_center.resize(nv);
  Dflux_x_right.resize(nv);
  for (vi=0;vi<nv;vi++)
    {
      flux_x[vi].newsize(Vec::LOCAL,local_nxNodes+1);
      flux_x[vi] = 0.0;
      Dflux_x_center[vi].resize(nv);
      Dflux_x_right[vi].resize(nv);
      for (vj=0;vj<nv;vj++)
        {
          Dflux_x_center[vi][vj].newsize(Vec::LOCAL,local_nxNodes+1);
          Dflux_x_center[vi][vj] = 0.0;
          Dflux_x_right[vi][vj].newsize(Vec::LOCAL,local_nxNodes+1);
          Dflux_x_right[vi][vj] = 0.0;
        }
      //mwf initialize difference increments
      pDiff_x[vi].newsize(Vec::LOCAL,node.ghost_nxNodes+1);
      pDiff_x[vi] = 0.0;
      DpDiff_x_center[vi].newsize(Vec::LOCAL,node.ghost_nxNodes+1);
      DpDiff_x_right[vi].newsize(Vec::LOCAL,node.ghost_nxNodes+1);
      DpDiff_x_right[vi] = 0.0;

    }

  KrX.newsize(Vec::LOCAL,node.ghost_nxNodes+1);
  RhoX.newsize(Vec::LOCAL,node.ghost_nxNodes+1);
  DKrX_center.newsize(Vec::LOCAL,node.ghost_nxNodes+1);
  DKrX_right.newsize(Vec::LOCAL,node.ghost_nxNodes+1);
  DRhoX_center.newsize(Vec::LOCAL,node.ghost_nxNodes+1);
  DRhoX_right.newsize(Vec::LOCAL,node.ghost_nxNodes+1);
    
  DKrX_center = 0.0;
  DKrX_right  = 0.0;
  DRhoX_center= 0.0;
  DRhoX_right = 0.0;

}

template<class BC, int nv>
inline void DivKgradECDM1d<BC,nv>::setPdiff_x(const Vec* Rho, const Vec* P)
{                  
  pDiff_x[0][node.interRightGhost] = 
    (-DpDx(P[0]))+
    RhoX[node.interRightGhost]*gravX[node.interRightGhost];
}

template<class BC, int nv>
inline void DivKgradECDM1d<BC,nv>::setDPdiff_x(const Vec* Rho, const Vec* DRho, const Vec* P)
{                  

  DpDiff_x_center[0][node.interRightGhost] = 
    2.0/(hXcen[node.center]+hXcen[node.right])
    + DRhoX_center[node.interRightGhost]*gravX[node.interRightGhost];

  DpDiff_x_right[0][node.interRightGhost] = 
    -2.0/(hXcen[node.center]+hXcen[node.right])
    + DRhoX_right[node.interRightGhost]*gravX[node.interRightGhost];
}

//get flux deriv for extra conductivity term
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getInteriorFluxDerivTerm_x(int face)
{   
  real fluxDeriv(0.0),gradVal(0.0),coefDeriv(0.0);

  //this is like normal flux calc on right side
  //except derivative term is wrt left variable
  if (face == 1)
    {
      //now get corner terms at x+1/2
      real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interRightGhost];
    
	
      gradVal = Ki11Avg*pDiff_x[0][node.interRightGhost];

      //at interface KrX[interRightGhost]=KrX[interLeftGhost]
      coefDeriv = RhoX[node.interLeftGhost]*DKrX_center[node.interLeftGhost]
        +DRhoX_center[node.interLeftGhost]*KrX[node.interLeftGhost];

      fluxDeriv = gradVal*coefDeriv;

    }
  else
    {
      //this is a gradient on the left side of the interface
      //except derivative term is wrt right variable
      assert(face == -1);
      //now get corner terms at x+1/2
      real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interLeftGhost];

      gradVal = Ki11Avg*pDiff_x[0][node.interLeftGhost];

      //at interface KrX[interRightGhost]=KrX[interLeftGhost]
      coefDeriv = RhoX[node.interLeftGhost]*DKrX_right[node.interLeftGhost]
        +DRhoX_right[node.interLeftGhost]*KrX[node.interLeftGhost];

      fluxDeriv = gradVal*coefDeriv;
    }//end else
  return fluxDeriv;
}

//use this for upwinding
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getSinglePhaseFlux_x()
{
  real fluxVal(0.0);
  //now get corner terms at x+1/2
  real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interRightGhost];

      
  fluxVal = Ki11Avg*pDiff_x[0][node.interRightGhost];

  return fluxVal;
}

//mwf do a separate one for 1d for now
template<class BC, int nv>
inline void DivKgradECDM1d<BC,nv>::setFluxNew_x(const VecVec* K, const Vec* Kr,
                                 const Vec* Rho, const Vec* P)
{                  
  for (vi=0;vi<nv;vi++)
    {
      real rhoAvgKrAvg =  RhoX[node.interRightGhost]*
        KrX[node.interRightGhost];
      //now get corner terms at x+1/2
      real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interRightGhost];

      flux_x[vi][node.interRight] =
        rhoAvgKrAvg*Ki11Avg*pDiff_x[vi][node.interRightGhost];

    }//end vi loop
}

//mwf put this one in for 1d for now
template<class BC, int nv>
inline void DivKgradECDM1d<BC,nv>::setDFluxNew_x(const VecVec* K, const Vec* Kr,
                                  const Vec* Rho, const Vec* P, 
                                  const VecVecVec& DKr, const Vec* DRho)
{

  for (vi=0;vi<nv;vi++)
    {
      for (vj=0;vj<nv;vj++)
        {
          //interior
          //mwf from flux calculation
          real rhoAvgKrAvg=  RhoX[node.interRightGhost]
            *KrX[node.interRightGhost];
          real Ki11Avg  = (*KWsMatCorner)[XX][MINUS][node.interRightGhost];

          //mwf derivative terms
          real DRhoAvgKrAvg_center= DRhoX_center[node.interRightGhost]
            *KrX[node.interRightGhost];
          real DRhoAvgKrAvg_right = DRhoX_right[node.interRightGhost]
            *KrX[node.interRightGhost];
	      
          real rhoAvgDKrAvg_center= RhoX[node.interRightGhost] 
            *DKrX_center[node.interRightGhost];

          real rhoAvgDKrAvg_right =  RhoX[node.interRightGhost] 
            *DKrX_right[node.interRightGhost];
	      
          //do derivatives of coefficients first 
          //mwf coefficients only a function of center and right (for now)
          //mwf relative perm can be a function of the other phase
          Dflux_x_center[vi][vj][node.interRight] = rhoAvgDKrAvg_center*
            Ki11Avg*pDiff_x[vi][node.interRightGhost];
	      
	      
          Dflux_x_right[vi][vj][node.interRight] = rhoAvgDKrAvg_right*
            Ki11Avg*pDiff_x[vi][node.interRightGhost];
	      

          if (vi == vj)
            { 
              //take care of DRho terms first
              Dflux_x_center[vi][vj][node.interRight] += DRhoAvgKrAvg_center*
                Ki11Avg*pDiff_x[vi][node.interRightGhost];
		  

              Dflux_x_right[vi][vj][node.interRight] += DRhoAvgKrAvg_right*
                Ki11Avg*pDiff_x[vi][node.interRightGhost];
		  
		  
              //now start getting derivative terms from appearance of 
              //actual pressure
              //finish out center and right
              Dflux_x_center[vi][vj][node.interRight] += rhoAvgKrAvg*
                Ki11Avg*DpDiff_x_center[vi][node.interRightGhost];
	
		  
              Dflux_x_right[vi][vj][node.interRight] += rhoAvgKrAvg*
                Ki11Avg*DpDiff_x_right[vi][node.interRightGhost];
		  
            }//end if (vi==vj)
        }//end vj				 
    }//end vi
}//end function

//====== gradient derivative terms
//try and use these for no flow or gradient conditions 
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getSinglePhaseFlux_xRightBnd()
{
  //have to use values at one cell in?
  real fluxVal(0.0);
  //now get corner terms at x+1/2
  real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interLeftGhost];

  fluxVal = Ki11Avg*pDiff_x[0][node.interLeftGhost];

  return fluxVal;
}
  
//same as right flux for left boundary?
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getSinglePhaseFlux_xLeftBnd()
{
  real fluxVal(0.0);
  //now get corner terms at x+1/2
  real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interRightGhost];

  fluxVal = Ki11Avg*pDiff_x[0][node.interRightGhost];

  return fluxVal;
}

//======xRightBnd derivatives=====
//try and use these for no flow or gradient conditions 
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getDSinglePhaseFlux_xRightBnd_center()
{
  //have to use values at one cell in?
  real DfluxVal(0.0);
  //now get corner terms at x+1/2
  real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interLeftGhost];

      
  DfluxVal = Ki11Avg*DpDiff_x_right[0][node.interLeftGhost];

  return DfluxVal;
}
//try and use these for no flow or gradient conditions 
template<class BC, int nv>
inline virtual real DivKgradECDM1d<BC,nv>::getDSinglePhaseFlux_xRightBnd_left()
{
  //have to use values at one cell in?
  real DfluxVal(0.0);
  //now get corner terms at x+1/2
  real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interLeftGhost];
      
  DfluxVal = Ki11Avg*DpDiff_x_center[0][node.interLeftGhost];

  return DfluxVal;
}
//try and use these for no flow or gradient conditions 
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getDSinglePhaseFlux_xRightBnd_leftBack()
{
  return 0.0;
}
//try and use these for no flow or gradient conditions 
template<class BC, int nv>
inline real  DivKgradECDM1d<BC,nv>::getDSinglePhaseFlux_xRightBnd_leftFront()
{
  return 0.0;
}

//====== xLeftBnd derivatives =====
//same as right flux for left boundary?
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getDSinglePhaseFlux_xLeftBnd_center()
{
  real DfluxVal(0.0);
  //now get corner terms at x+1/2
  real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interRightGhost];
      
  DfluxVal = Ki11Avg*DpDiff_x_center[0][node.interRightGhost];

  return DfluxVal;
}
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getDSinglePhaseFlux_xLeftBnd_right()
{
  real DfluxVal(0.0);
  //now get corner terms at x+1/2
  real Ki11Avg     = (*KWsMatCorner)[XX][MINUS][node.interRightGhost];
      
  DfluxVal = Ki11Avg*DpDiff_x_right[0][node.interRightGhost];

  return DfluxVal;
}
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getDSinglePhaseFlux_xLeftBnd_rightBack()
{
  return 0.0;
}
//same as right flux for left boundary?
template<class BC, int nv>
inline real DivKgradECDM1d<BC,nv>::getDSinglePhaseFlux_xLeftBnd_rightFront()
{  
  return 0.0;
}

//try and do ECDM method here
template<class BC, int nv>
virtual void DivKgradECDM1d<BC,nv>::computeDivergence(const VecVec* K, const Vec* Kr,
                               const Vec* Rho, 
                               const Vec* P, 
                               Vec* Div)
{
  //first have to set pDiff_x now
  //also RhoX, and KrX
  //order RhoX, pDiff_X, KrX
    
  //left ghost line
  if (node.ghost_xOffSet)
    {
      int kStart = max(-1,-node.ghost_xOffSet);
	
      node.localIndex(kStart);
      setRhoX(Rho,P);
      setPdiff_x(Rho,P);
	
      //mwf only get KrX at interfaces I'll need (not x direction
      //mwf along y ghost row)  because have to use perm arrays for
      //mwf driveForce
      //if (j >= 0 && j < local_nyNodes)
      setKrX(Kr,P);
    }
	  
  
  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    {
      int kEnd   =  local_nxNodes-1;
	
      node.localIndex(kEnd);
      setRhoX(Rho,P);
      setPdiff_x(Rho,P);
	
      //mwf only get KrX at interfaces I'll need (not x direction
      //mwf along y ghost row)  because have to use perm arrays for
      //mwf driveForce
      //if (j >= 0 && j < local_nyNodes)
      setKrX(Kr,P);
	
    }
  //now do interior for X
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(k);
      setRhoX(Rho,P);
      setPdiff_x(Rho,P);
	
      //mwf only get KrX at interfaces I'll need (not x direction
      //mwf along y ghost row)  because have to use perm arrays for
      //mwf driveForce
      //if (j >= 0 && j < local_nyNodes)
      setKrX(Kr,P);
	      
    }

  //k=0(interLeft) or k=local_nxNodes-1 (interRight) 
  //are done either by ghost boundary or BC class

  // mwf now account for Dirichlet BC's in pressure differences
      
  for (vi=0;vi<nv;vi++)
    {
      bc[vi].applyBoundaryConditionsToPressureDifferences(node,pDiff_x[vi]);
    }//end loop for applying BC's to pressure differences

  //now account for interior conditions for const. rels
#ifdef USE_CONT_INTERFACE_TERMS_ECDM
  for (vi=0;vi<nv;vi++)
    {
      bc[vi].setCurrentDivKgrad(this);
      bc[vi].applyInteriorConditionToInterfaceValues(node,Kr,Rho,
                                                     KrX,RhoX);
    }//end loop for applying BC's to const. rels
#endif
  //try and look at rel perms now
  //        ofstream krxout("KrX.grf");
  //        krxout <<KrX<<endl;

  //mwf now maybe I can set the fluxes
  //////////
  //left ghost line
  if (node.ghost_xOffSet)
    {
      int kStart = max(-1,-node.ghost_xOffSet);
      //mwf changed this to be -1 only 
      node.localIndex(kStart);
      setFluxNew_x(K,Kr,Rho,P);
    }

  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    {
      node.localIndex(local_nxNodes-1);
      setFluxNew_x(K,Kr,Rho,P);
    }

  //interior for x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(k);
      setFluxNew_x(K,Kr,Rho,P);
    }
	
  //mwf now account for lagrange multipliers in boundary conditions
      
  for (vi=0;vi<nv;vi++)
    {
      bc[vi].applyBoundaryConditionsToFluxes(node,flux_x[vi]);

    }//end loop for applying bc's to fluxes

  for (k=0;k<local_nxNodes;k++)
    {
      node.localIndex(k);
      //mwf this is minus the divergence!
      for (vi=0;vi<nv;vi++)
        {
          Div[vi][node.center_noGhost] = flux_x[vi][node.interLeft] 
            - flux_x[vi][node.interRight];
        }
    }
    
}  
  
////////////////////////////////////////
template<class BC, int nv>
virtual void DivKgradECDM1d<BC,nv>::DcomputeDivergenceJac(const VecVec* K, const Vec* Kr,
                                  const Vec* Rho, 
                                  const Vec* P, 
                                  const VecVecVec& DKr, const Vec* DRho, 
                                  VecVecVecVec& DivJac)
{
  //first have to set DpDiff_x, now
  //also DRhoX, and DKrX
  //order DRhoX, DpDiff_X, DKrX
  //left ghost line
  if (node.ghost_xOffSet)
    {
      int kStart = max(-1,-node.ghost_xOffSet);

      node.localIndex(kStart);
      setDRhoX(DRho,P);
      setDPdiff_x(Rho,DRho,P);
      //mwf only get KrX at interfaces I'll need (not x direction
      //mwf along y ghost row)  because have to use perm arrays for
      //mwf driveForce
      //if (j >= 0 && j < local_nyNodes)
      setDKrX(DKr,P);
    }

  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    {
      int kEnd   =  local_nxNodes-1;

      node.localIndex(kEnd);
      setDRhoX(DRho,P);
      setDPdiff_x(Rho,DRho,P);
      //mwf only get KrX at interfaces I'll need (not x direction
      //mwf along y ghost row)  because have to use perm arrays for
      //mwf driveForce
      //if (j >= 0 && j < local_nyNodes)
      setDKrX(DKr,P);

    }

  //now do interior for X
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(k);
      setDRhoX(DRho,P);
      setDPdiff_x(Rho,DRho,P);

      //mwf only get KrX at interfaces I'll need (not x direction
      //mwf along y ghost row)  because have to use perm arrays for
      //mwf driveForce
      //if (j >= 0 && j < local_nyNodes)
      setDKrX(DKr,P);

    }
        
  //k=0(interLeft) or k=local_nxNodes-1 (interRight) 
  //are done either by ghost boundary or BC class

  for (vi=0;vi<nv;vi++)
    {
      bc[vi].applyBoundaryConditionsToPressureDerivatives(node,
                                                          DpDiff_x_center[vi],
                                                          DpDiff_x_right[vi]);
    }


#ifdef USE_CONT_INTERFACE_TERMS_ECDM
  //now account for interior conditions for const. rels
  for (vi=0;vi<nv;vi++)
    {
      bc[vi].setCurrentDivKgrad(this);
      bc[vi].applyInteriorConditionToInterfaceValueDerivs(node,
                                                          this,
                                                          DKr,DRho);
    }//end loop for applying BC's to const. rels
#endif

  //mwf now maybe I can set the fluxes
  //////////
  //left ghost line
  if (node.ghost_xOffSet)
    {
      int kStart = max(-1,-node.ghost_xOffSet);
      node.localIndex(kStart);
      setDFluxNew_x(K,Kr,Rho,P,DKr,DRho);
    }

  //right ghost line
  if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
    {
      node.localIndex(local_nxNodes-1);
      setDFluxNew_x(K,Kr,Rho,P,DKr,DRho);
    }

  //interior for x fluxes
  for (k=0;k<local_nxNodes-1;k++)
    {
      node.localIndex(j,k);
      setDFluxNew_x(K,Kr,Rho,P,DKr,DRho);
    }
	

  //mwf now account for Dirichlet BC's in fluxes by maybe upwinding
  //mwf using dirichlet values
      
  //here since should be all zero?
  for (vi=0;vi<nv;vi++)
    {
      bc[vi].applyBoundaryConditionDerivativesECDM(node,
                                                   Dflux_x_center[vi][vi], 
                                                   Dflux_x_right[vi][vi]);
    }
      
  for (k=0;k<local_nxNodes;k++)
    {
      node.localIndex(k);
      for (vi=0;vi<nv;vi++)
        {
          for (vj=0;vj<nv;vj++)
            {
              DivJac[vi][vj][LEFT][node.center_noGhost] = 
                Dflux_x_center[vi][vj][node.interLeft];
		    

              DivJac[vi][vj][CENTER][node.center_noGhost] = 
                Dflux_x_right[vi][vj][node.interLeft] 
                -Dflux_x_center[vi][vj][node.interRight]; 
		     
                    
              DivJac[vi][vj][RIGHT][node.center_noGhost]  = 
                -Dflux_x_right[vi][vj][node.interRight];
            }//vj
        }//vi
    }//end loop through for jac

}//end function
  
template<class BC, int nv>
virtual void DivKgradECDM1d<BC,nv>::computeDivergence(const Vec* K, 
                               const Vec* Rho, 
                               const Vec* P, 
                               Vec* Div)
{ assert(0); }


template<class BC, int nv>
virtual void DivKgradECDM1d<BC,nv>::computeDivergenceJac(const Vec* K,
                                  const Vec* Rho, 
                                  const Vec* P,
                                  const VecVecVec& DK, const Vec* DRho, 
                                  VecVecVecVec& DivJac)
{ assert(0); }   

template<class BC, int nv>
class DivKgradECDM2d : public DivKgradECDM<BC,nv>
{
public:
  DivKgradECDM2d(BC* bcIn,
             StencilMM& nodeIn,
             int nxNodesIn,int nyNodesIn, 
             real gxIn, real gyIn,
             real oneOverdxIn, real oneOverdyIn):
    DivKgradECDM<BC,nv>(bcIn,nodeIn)
  {
    nxNodes=nxNodesIn;
    local_nxNodes = nodeIn.local_nxNodes;
    gx=gxIn;
    oneOverdx = oneOverdxIn;

    nyNodes=nyNodesIn;
    local_nyNodes = nodeIn.local_nyNodes;
    gy=gyIn;
    oneOverdy = oneOverdyIn;
   
    //what do I need to do here?
    hXcen.newsize(Vec::LOCAL,node.da_);
    hYcen.newsize(Vec::LOCAL,node.da_);
    hZcen.newsize(Vec::LOCAL,node.da_);

    hXcen = 1.0/oneOverdx;
    hYcen = 1.0/oneOverdy;
    hZcen = 1.0;

    //setup gravity term for ghost interfaces
    gravX.newsize(Vec::LOCAL,
		  node.ghost_nyNodes*(node.ghost_nxNodes+1));
    gravY.newsize(Vec::LOCAL,
		  (node.ghost_nyNodes+1)*node.ghost_nxNodes);
    gravX = gx;
    gravY = gy;

    Dflux_x_center.resize(nv);
    Dflux_x_right.resize(nv);
    Dflux_y_center.resize(nv);
    Dflux_y_back.resize(nv);
    //mwf now add more terms for cross derivatives
    Dflux_x_front.resize(nv);
    Dflux_x_back.resize(nv);
    Dflux_x_rightFront.resize(nv);
    Dflux_x_rightBack.resize(nv);
    //
    Dflux_y_right.resize(nv);
    Dflux_y_left.resize(nv);
    Dflux_y_backRight.resize(nv);
    Dflux_y_backLeft.resize(nv);

    for (vi=0;vi<nv;vi++)
      {
	flux_x[vi].newsize(Vec::LOCAL,local_nyNodes*(local_nxNodes+1));
	flux_y[vi].newsize(Vec::LOCAL,(local_nyNodes+1)*local_nxNodes);
	flux_x[vi] = 0.0;
	flux_y[vi] = 0.0;
	Dflux_x_center[vi].resize(nv);
	Dflux_x_right[vi].resize(nv);
	Dflux_y_center[vi].resize(nv);
	Dflux_y_back[vi].resize(nv);

	//mwf now add more terms for cross derivatives
	Dflux_x_front[vi].resize(nv);
	Dflux_x_back[vi].resize(nv);
	Dflux_x_rightFront[vi].resize(nv);
	Dflux_x_rightBack[vi].resize(nv);
	//
	Dflux_y_right[vi].resize(nv);
	Dflux_y_left[vi].resize(nv);
	Dflux_y_backRight[vi].resize(nv);
	Dflux_y_backLeft[vi].resize(nv);


	for (vj=0;vj<nv;vj++)
	  {
	    Dflux_x_center[vi][vj].newsize(Vec::LOCAL,
					   local_nyNodes*(local_nxNodes+1));
	    Dflux_x_center[vi][vj] = 0.0;
	    Dflux_x_right[vi][vj].newsize(Vec::LOCAL,
					  local_nyNodes*(local_nxNodes+1));
	    Dflux_x_right[vi][vj] = 0.0;
	    Dflux_y_center[vi][vj].newsize(Vec::LOCAL,
					   (local_nyNodes+1)*local_nxNodes);
	    Dflux_y_center[vi][vj] = 0.0;
	    Dflux_y_back[vi][vj].newsize(Vec::LOCAL,
					 (local_nyNodes+1)*local_nxNodes);
	    Dflux_y_back[vi][vj] = 0.0;

	    //mwf now do cross derivative terms
	    Dflux_x_front[vi][vj].newsize(Vec::LOCAL,
					   local_nyNodes*(local_nxNodes+1));
	    Dflux_x_front[vi][vj] = 0.0;

	    Dflux_x_back[vi][vj].newsize(Vec::LOCAL,
					   local_nyNodes*(local_nxNodes+1));
	    Dflux_x_back[vi][vj] = 0.0;

	    Dflux_x_rightFront[vi][vj].newsize(Vec::LOCAL,
					   local_nyNodes*(local_nxNodes+1));
	    Dflux_x_rightFront[vi][vj] = 0.0;

	    Dflux_x_rightBack[vi][vj].newsize(Vec::LOCAL,
					   local_nyNodes*(local_nxNodes+1));
	    Dflux_x_rightBack[vi][vj] = 0.0;
	    
	    //
	    Dflux_y_right[vi][vj].newsize(Vec::LOCAL,
					   (local_nyNodes+1)*local_nxNodes);
	    Dflux_y_right[vi][vj] = 0.0;

	    Dflux_y_left[vi][vj].newsize(Vec::LOCAL,
					 (local_nyNodes+1)*local_nxNodes);
	    Dflux_y_left[vi][vj] = 0.0;

	    Dflux_y_backRight[vi][vj].newsize(Vec::LOCAL,
					      (local_nyNodes+1)*local_nxNodes);
	    Dflux_y_backRight[vi][vj] = 0.0;

	    Dflux_y_backLeft[vi][vj].newsize(Vec::LOCAL,
					     (local_nyNodes+1)*local_nxNodes);
	    Dflux_y_backLeft[vi][vj] = 0.0;

	    
	  }
	//mwf initialize difference increments
	pDiff_x[vi].newsize(Vec::LOCAL,
			    node.ghost_nyNodes*(node.ghost_nxNodes+1));
	pDiff_y[vi].newsize(Vec::LOCAL,
			    (node.ghost_nyNodes+1)*node.ghost_nxNodes);
	pDiff_x[vi] = 0.0;
	pDiff_y[vi] = 0.0;
	
	//mwf now do derivative terms
	DpDiff_x_center[vi].newsize(Vec::LOCAL,
				    node.ghost_nyNodes*(node.ghost_nxNodes+1));
	DpDiff_x_right[vi].newsize(Vec::LOCAL,
				   node.ghost_nyNodes*(node.ghost_nxNodes+1));
	DpDiff_y_center[vi].newsize(Vec::LOCAL,
				    (node.ghost_nyNodes+1)*node.ghost_nxNodes);
	DpDiff_y_back[vi].newsize(Vec::LOCAL,
				  (node.ghost_nyNodes+1)*node.ghost_nxNodes);
	DpDiff_x_center[vi] = 0.0;
	DpDiff_x_right[vi]  = 0.0;
	DpDiff_y_center[vi] = 0.0;
	DpDiff_y_back[vi]   = 0.0;
      }
	
    KrX.newsize(Vec::LOCAL,
		node.ghost_nyNodes*(node.ghost_nxNodes+1));
    RhoX.newsize(Vec::LOCAL,
			node.ghost_nyNodes*(node.ghost_nxNodes+1));
    KrY.newsize(Vec::LOCAL,
		(node.ghost_nyNodes+1)*node.ghost_nxNodes);
    RhoY.newsize(Vec::LOCAL,
		 (node.ghost_nyNodes+1)*node.ghost_nxNodes);
    KrX = 0.0;
    RhoX = 0.0;
    KrY = 0.0;
    RhoY = 0.0;

    DKrX_center.newsize(Vec::LOCAL,
			node.ghost_nyNodes*(node.ghost_nxNodes+1));
    DKrX_right.newsize(Vec::LOCAL,
		       node.ghost_nyNodes*(node.ghost_nxNodes+1));
    DRhoX_center.newsize(Vec::LOCAL,
			 node.ghost_nyNodes*(node.ghost_nxNodes+1));
    DRhoX_right.newsize(Vec::LOCAL,
			node.ghost_nyNodes*(node.ghost_nxNodes+1));
    DKrY_center.newsize(Vec::LOCAL,
			(node.ghost_nyNodes+1)*node.ghost_nxNodes);
    DKrY_back.newsize(Vec::LOCAL,
		      (node.ghost_nyNodes+1)*node.ghost_nxNodes);
    DRhoY_center.newsize(Vec::LOCAL,
			(node.ghost_nyNodes+1)*node.ghost_nxNodes);
    DRhoY_back.newsize(Vec::LOCAL,
		       (node.ghost_nyNodes+1)*node.ghost_nxNodes);

    DKrX_center = 0.0;
    DKrX_right  = 0.0;
    DRhoX_center= 0.0;
    DRhoX_right = 0.0;
    DKrY_center = 0.0;
    DKrY_back   = 0.0;
    DRhoY_center= 0.0;
    DRhoY_back  = 0.0;

  
  }
  virtual ~DivKgradECDM2d(){}

  //try and do ECDM method here
  virtual void computeDivergence(const VecVec* K, const Vec* Kr,
				 const Vec* Rho, 
				 const Vec* P, 
                                 Vec* Div)
    {


      //first have to set pDiff_x,pDiff_y now
      //also RhoX, and KrX
      //order RhoX, pDiff_X, KrX

      //bottom ghost line
      if (node.ghost_yOffSet)
	{
	  int kStart = max(-1,-node.ghost_xOffSet);
	  int kEnd   = max(local_nxNodes,
			   node.ghost_nxNodes-node.ghost_xOffSet);
	  int jStart = max(-1,-node.ghost_yOffSet);
	  for(k = kStart; k< kEnd ;k++)
	    {
	      node.localIndex(jStart,k);
	      setRhoY(Rho,P);
	      setPdiff_y(Rho,P);

	      //mwf only get KrY at interfaces I'll need (not y direction
	      //mwf along x ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (k >= 0 && k < local_nxNodes)
	      setKrY(Kr,P);

	    }
	}
      //top ghost line
      if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
	{
	  int kStart = max(-1,-node.ghost_xOffSet);
	  int kEnd   = max(local_nxNodes,
			   node.ghost_nxNodes-node.ghost_xOffSet);
	  int jEnd   = local_nyNodes-1;

	  for(k = kStart;
	      k< kEnd ;k++)
	    {
	      node.localIndex(jEnd,k);
	      setRhoY(Rho,P);
	      setPdiff_y(Rho,P);
	      
	      //mwf only get KrY at interfaces I'll need (not y direction
	      //mwf along x ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (k >= 0 && k < local_nxNodes)
	      setKrY(Kr,P);
	    }
	}
      //left ghost line
      if (node.ghost_xOffSet)
	{
	  int jStart = max(-1,-node.ghost_yOffSet);
	  int jEnd   = max(local_nyNodes,
			   node.ghost_nyNodes-node.ghost_yOffSet);
	  int kStart = max(-1,-node.ghost_xOffSet);

	  for (j = jStart;
	       j< jEnd ;j++)
	    {
	      node.localIndex(j,kStart);
	      setRhoX(Rho,P);
	      setPdiff_x(Rho,P);

	      //mwf only get KrX at interfaces I'll need (not x direction
	      //mwf along y ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (j >= 0 && j < local_nyNodes)
	      setKrX(Kr,P);
	    }
	  
	}
      //right ghost line
      if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
	{
	  int jStart = max(-1,-node.ghost_yOffSet);
	  int jEnd   = max(local_nyNodes,
			   node.ghost_nyNodes-node.ghost_yOffSet);
	  int kEnd   =  local_nxNodes-1;

	  for (j = jStart; j< jEnd ;j++)
	    {
	      node.localIndex(j,kEnd);
	      setRhoX(Rho,P);
	      setPdiff_x(Rho,P);

	      //mwf only get KrX at interfaces I'll need (not x direction
	      //mwf along y ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (j >= 0 && j < local_nyNodes)
	      setKrX(Kr,P);
	    }
	}
      //now do interior for X
      //mwf June 24, should this be node.ghost_nyNodes-ghost_yOffset?
      //mwf June 24 was node.ghost_nyNodes-1
      int jLoopMax = max(local_nyNodes,
			 node.ghost_nyNodes-node.ghost_yOffSet);
      int jStart   = max(-1,-node.ghost_yOffSet);
      for (j = jStart; j< jLoopMax  ;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(j,k);
	      setRhoX(Rho,P);
              setPdiff_x(Rho,P);

	      //mwf only get KrX at interfaces I'll need (not x direction
	      //mwf along y ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (j >= 0 && j < local_nyNodes)
	      setKrX(Kr,P);
	      
            }
        }
      //k=0(interLeft) or k=local_nxNodes-1 (interRight) 
      //are done either by ghost boundary or BC class

      //now do interior for Y
      //mwf June 24, should this be node.ghost_nxNodes-ghost_xOffset?
      //mwf June 24 was node.ghost_nxNodes-1
      int kLoopMax = max(local_nxNodes,
			 node.ghost_nxNodes-node.ghost_xOffSet);
      int kStart   = max(-1,-node.ghost_xOffSet);
      for (j=0; j < local_nyNodes-1; j++)
	{
	  for(k = kStart;
	      k< kLoopMax ;k++)
	    {
	      node.localIndex(j,k);
	      setRhoY(Rho,P);
	      setPdiff_y(Rho,P);
	      //mwf only get KrY at interfaces I'll need (not y direction
	      //mwf along x ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (k >= 0 && k < local_nxNodes)
	      setKrY(Kr,P);

	    }
	}
      //j=0(interFront) or j=local_nyNodes-1 (interBack) 
      //are done either by ghost boundary or BC class

      // mwf now account for Dirichlet BC's in pressure differences
      
      for (vi=0;vi<nv;vi++)
        {
	  bc[vi].applyBoundaryConditionsToPressureDifferences(node,pDiff_x[vi],
							      pDiff_y[vi]);
	}//end loop for applying BC's to pressure differences

      //now account for interior conditions for const. rels
#ifdef USE_CONT_INTERFACE_TERMS_ECDM
      for (vi=0;vi<nv;vi++)
        {
	  bc[vi].setCurrentDivKgrad(this);
	  bc[vi].applyInteriorConditionToInterfaceValues(node,Kr,Rho,
							 KrX,KrY,
							 RhoX,RhoY);
	}//end loop for applying BC's to const. rels
#endif
      //try and look at rel perms now
//        ofstream krxout("KrX.grf");
//        krxout <<KrX<<endl;
//        ofstream kryout("KrY.grf");
//        kryout <<KrY<<endl;

      //mwf now maybe I can set the fluxes
      //////////
      //bottom ghost line
      if (node.ghost_yOffSet)
	{
	  int jStart   = max(-1,-node.ghost_yOffSet);
	  for(k=0;k<local_nxNodes;k++)
	    {
	      //mwf changed this to be -1 only
	      node.localIndex(jStart,k);
	      setFluxNew_y(K,Kr,Rho,P);
	    }
	}
      //top ghost line
      if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
        for(k=0;k<local_nxNodes;k++)
          {
            node.localIndex(local_nyNodes-1,k);
            setFluxNew_y(K,Kr,Rho,P);
          }
      
      //left ghost line
      if (node.ghost_xOffSet)
        for (j=0;j<local_nyNodes;j++)
          {
	    int kStart = max(-1,-node.ghost_xOffSet);
	    //mwf changed this to be -1 only 
            node.localIndex(j,kStart);
            setFluxNew_x(K,Kr,Rho,P);
          }

      //right ghost line
      if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
        for (j=0;j<local_nyNodes;j++)
          {
            node.localIndex(j,local_nxNodes-1);
            setFluxNew_x(K,Kr,Rho,P);
          }

      //interior for x fluxes
      for (j=0;j<local_nyNodes;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(j,k);
              setFluxNew_x(K,Kr,Rho,P);
	    }
	}

      //interior for y fluxes
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes;k++)
            {
              node.localIndex(j,k);
              setFluxNew_y(K,Kr,Rho,P);
	    }
	}

      //mwf now account for lagrange multipliers in boundary conditions
      
      for (vi=0;vi<nv;vi++)
        {
	  bc[vi].applyBoundaryConditionsToFluxes(node,flux_x[vi],
						 flux_y[vi]);

	}//end loop for applying bc's to fluxes

      for (j=0;j<local_nyNodes;j++)
        for (k=0;k<local_nxNodes;k++)
          {
            node.localIndex(j,k);
	    //mwf this is minus the divergence!
            for (vi=0;vi<nv;vi++)
              {
                Div[vi][node.center_noGhost] = flux_x[vi][node.interLeft] 
		  - flux_x[vi][node.interRight]
                  +flux_y[vi][node.interFront] 
		  - flux_y[vi][node.interBack];
              }
          }
    }
  ////////////////////////////////////////
  virtual void computeDivergenceJac(const VecVec* K, const Vec* Kr,
				    const Vec* Rho, 
				    const Vec* P, 
                                    const VecVecVec& DKr, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
    {


      //first have to set DpDiff_x,DpDiff_y now
      //also DRhoX, and DKrX
      //order DRhoX, DpDiff_X, DKrX
      //bottom ghost line
      if (node.ghost_yOffSet)
	{
	  int kStart = max(-1,-node.ghost_xOffSet);
	  int kEnd   = max(local_nxNodes,
			   node.ghost_nxNodes-node.ghost_xOffSet);
	  int jStart = max(-1,-node.ghost_yOffSet);
	  for(k = kStart;
	      k< kEnd ;k++)
	    {
	      node.localIndex(jStart,k);
	      setDRhoY(DRho,P);
	      setDPdiff_y(Rho,DRho,P);
	      //mwf only get KrY at interfaces I'll need (not y direction
	      //mwf along x ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (k >= 0 && k < local_nxNodes)
	      setDKrY(DKr,P);
	    }
	}
      //top ghost line
      if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
	{
	  int kStart = max(-1,-node.ghost_xOffSet);
	  int kEnd   = max(local_nxNodes,
			   node.ghost_nxNodes-node.ghost_xOffSet);
	  int jEnd   = local_nyNodes-1;

	  for(k = kStart;
	      k< kEnd ;k++)
	    {
	      node.localIndex(jEnd,k);
	      setDRhoY(DRho,P);
	      setDPdiff_y(Rho,DRho,P);
	      //mwf only get KrY at interfaces I'll need (not y direction
	      //mwf along x ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (k >= 0 && k < local_nxNodes)
	      setDKrY(DKr,P);

	    }
	}
      //left ghost line
      if (node.ghost_xOffSet)
	{
	  int jStart = max(-1,-node.ghost_yOffSet);
	  int jEnd   = max(local_nyNodes,
			   node.ghost_nyNodes-node.ghost_yOffSet);
	  int kStart = max(-1,-node.ghost_xOffSet);

	  for (j = jStart;
	       j< jEnd ;j++)
	    {
	      node.localIndex(j,kStart);
	      setDRhoX(DRho,P);
	      setDPdiff_x(Rho,DRho,P);
	      //mwf only get KrX at interfaces I'll need (not x direction
	      //mwf along y ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (j >= 0 && j < local_nyNodes)
	      setDKrX(DKr,P);
	    }
	}
      //right ghost line
      if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
	{
	  int jStart = max(-1,-node.ghost_yOffSet);
	  int jEnd   = max(local_nyNodes,
			   node.ghost_nyNodes-node.ghost_yOffSet);
	  int kEnd   =  local_nxNodes-1;

	  for (j = jStart;
	       j< jEnd ;j++)
	    {
	      node.localIndex(j,kEnd);
	      setDRhoX(DRho,P);
	      setDPdiff_x(Rho,DRho,P);
	      //mwf only get KrX at interfaces I'll need (not x direction
	      //mwf along y ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (j >= 0 && j < local_nyNodes)
	      setDKrX(DKr,P);

	    }
	}
      //now do interior for X
      //mwf June 24, should this be node.ghost_nyNodes-ghost_yOffset?
      //mwf June 24 was node.ghost_nyNodes-1
      int jLoopMax = max(local_nyNodes,
			 node.ghost_nyNodes-node.ghost_yOffSet);
      int jStart   = max(-1,-node.ghost_yOffSet);
      for (j = jStart;
	   j< jLoopMax  ;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(j,k);
	      setDRhoX(DRho,P);
              setDPdiff_x(Rho,DRho,P);

	      //mwf only get KrX at interfaces I'll need (not x direction
	      //mwf along y ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (j >= 0 && j < local_nyNodes)
	      setDKrX(DKr,P);

            }
        }
      //k=0(interLeft) or k=local_nxNodes-1 (interRight) 
      //are done either by ghost boundary or BC class

      //now do interior for Y
      //mwf June 24, should this be node.ghost_nxNodes-ghost_xOffset?
      //mwf June 24 was node.ghost_nxNodes-1
      int kLoopMax = max(local_nxNodes,
			 node.ghost_nxNodes-node.ghost_xOffSet);
      int kStart   = max(-1,-node.ghost_xOffSet);
      for (j=0; j < local_nyNodes-1; j++)
	{
	  for(k = kStart;
	      k< kLoopMax ;k++)
	    {
	      node.localIndex(j,k);
	      setDRhoY(DRho,P);
	      setDPdiff_y(Rho,DRho,P);
	      
	      //mwf only get KrY at interfaces I'll need (not y direction
	      //mwf along x ghost row)  because have to use perm arrays for
	      //mwf driveForce
	      //if (k >= 0 && k < local_nxNodes)
	      setDKrY(DKr,P);

	    }
	}
      //j=0(interFront) or j=local_nyNodes-1 (interBack) 
      //are done either by ghost boundary or BC class
      // mwf now account for Dirichlet BC's in pressure differences
      
      for (vi=0;vi<nv;vi++)
	{
	  bc[vi].applyBoundaryConditionsToPressureDerivatives(node,
							      DpDiff_x_center[vi],
							      DpDiff_x_right[vi],
							      DpDiff_y_center[vi],
							      DpDiff_y_back[vi]);
	}


#ifdef USE_CONT_INTERFACE_TERMS_ECDM
      //now account for interior conditions for const. rels
      for (vi=0;vi<nv;vi++)
        {
	  bc[vi].setCurrentDivKgrad(this);
	  bc[vi].applyInteriorConditionToInterfaceValueDerivs(node,
							      this,
							      DKr,DRho);
	}//end loop for applying BC's to const. rels
#endif

      //mwf now maybe I can set the fluxes
      //////////
      //bottom ghost line
      if (node.ghost_yOffSet)
	{
	  int jStart   = max(-1,-node.ghost_yOffSet);
	  for(k=0;k<local_nxNodes;k++)
	    {
	      //mwf changed this to -1
	      node.localIndex(jStart,k);
	      setDFluxNew_y(K,Kr,Rho,P,DKr,DRho);
	    }
	}
      //top ghost line
      if (node.ghost_nyNodes - node.ghost_yOffSet - node.local_nyNodes)
        for(k=0;k<local_nxNodes;k++)
          {
            node.localIndex(local_nyNodes-1,k);
            setDFluxNew_y(K,Kr,Rho,P,DKr,DRho);
          }
      
      //left ghost line
      if (node.ghost_xOffSet)
        for (j=0;j<local_nyNodes;j++)
          {
	    int kStart = max(-1,-node.ghost_xOffSet);
            node.localIndex(j,kStart);
            setDFluxNew_x(K,Kr,Rho,P,DKr,DRho);
          }

      //right ghost line
      if (node.ghost_nxNodes - node.ghost_xOffSet - node.local_nxNodes)
        for (j=0;j<local_nyNodes;j++)
          {
            node.localIndex(j,local_nxNodes-1);
            setDFluxNew_x(K,Kr,Rho,P,DKr,DRho);
          }

      //interior for x fluxes
      for (j=0;j<local_nyNodes;j++)
        {
          for (k=0;k<local_nxNodes-1;k++)
            {
              node.localIndex(j,k);
              setDFluxNew_x(K,Kr,Rho,P,DKr,DRho);
	    }
	}

      //interior for y fluxes
      for (j=0;j<local_nyNodes-1;j++)
        {
          for (k=0;k<local_nxNodes;k++)
            {
              node.localIndex(j,k);
              setDFluxNew_y(K,Kr,Rho,P,DKr,DRho);
	    }
	}

      //mwf now account for Dirichlet BC's in fluxes by maybe upwinding
      //mwf using dirichlet values
      
      //here since should be all zero?
      for (vi=0;vi<nv;vi++)
	{
	  bc[vi].applyBoundaryConditionDerivativesECDM(node,
						       Dflux_x_center[vi][vi], 
						       Dflux_x_right[vi][vi],
						       Dflux_x_back[vi][vi], 
						       Dflux_x_front[vi][vi],
						       Dflux_x_rightBack[vi][vi], 
						       Dflux_x_rightFront[vi][vi],
						       Dflux_y_center[vi][vi], 
						       Dflux_y_back[vi][vi], 
						       Dflux_y_right[vi][vi], 
						       Dflux_y_left[vi][vi],
						       Dflux_y_backRight[vi][vi], 
						       Dflux_y_backLeft[vi][vi]);

	}
      
      for (j=0;j<local_nyNodes;j++)
        for (k=0;k<local_nxNodes;k++)
          {
            node.localIndex(j,k);
            for (vi=0;vi<nv;vi++)
              {
                for (vj=0;vj<nv;vj++)
                  {
                    
                    DivJac[vi][vj][FRONT][node.center_noGhost] = 
                      Dflux_y_center[vi][vj][node.interFront]
		      -Dflux_x_front[vi][vj][node.interRight]
		      +Dflux_x_rightFront[vi][vj][node.interLeft];
                    
                    DivJac[vi][vj][LEFT][node.center_noGhost] = 
                      Dflux_x_center[vi][vj][node.interLeft]
		      -Dflux_y_left[vi][vj][node.interBack]
		      +Dflux_y_backLeft[vi][vj][node.interFront];

                    DivJac[vi][vj][CENTER][node.center_noGhost] = 
                      Dflux_x_right[vi][vj][node.interLeft] 
		      -Dflux_x_center[vi][vj][node.interRight] 
		      +Dflux_y_back[vi][vj][node.interFront] 
		      -Dflux_y_center[vi][vj][node.interBack]; 
                    
                    DivJac[vi][vj][RIGHT][node.center_noGhost]  = 
                      -Dflux_x_right[vi][vj][node.interRight]
		      -Dflux_y_right[vi][vj][node.interBack]
		      +Dflux_y_backRight[vi][vj][node.interFront];

                    DivJac[vi][vj][BACK][node.center_noGhost]  = 
                      -Dflux_y_back[vi][vj][node.interBack]
		      -Dflux_x_back[vi][vj][node.interRight]
		      +Dflux_x_rightBack[vi][vj][node.interLeft];
		    
		    //now do corner terms
		    //(i+1,j+1)
		    DivJac[vi][vj][RIGHT_BACK][node.center_noGhost]  =
		      -Dflux_x_rightBack[vi][vj][node.interRight]
		      -Dflux_y_backRight[vi][vj][node.interBack];

		    //(i+1,j-1)
		    DivJac[vi][vj][RIGHT_FRONT][node.center_noGhost]  =
		      -Dflux_x_rightFront[vi][vj][node.interRight]
		      +Dflux_y_right[vi][vj][node.interFront];

		    //(i-1,j+1)
		    DivJac[vi][vj][LEFT_BACK][node.center_noGhost]  =
		      Dflux_x_back[vi][vj][node.interLeft]
		      -Dflux_y_backLeft[vi][vj][node.interBack];

		    //(i-1,j-1)
		    DivJac[vi][vj][LEFT_FRONT][node.center_noGhost]  =
		      Dflux_x_front[vi][vj][node.interLeft]
		      +Dflux_y_left[vi][vj][node.interFront];

                  }//vj
              }//vi
          }//end loop through for jac

    }//end function

  ///////////////////////////////////////////
  //mwf not correct, maybe works like standard finite difference still
  virtual void computeDivergence(const Vec* K, const Vec* Rho, const Vec* P, 
                                 Vec* Div)
    {
      assert(0);
    }
  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, 
				    const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
    {
      assert(0);
    }
};  

template<class BC, int nv>
class DivKgradECDM3d : public DivKgradECDM<BC,nv>
{
public:
  DivKgradECDM3d(BC* bcIn,
             StencilMM& nodeIn,
             int nxNodesIn,int nyNodesIn,int nzNodesIn, 
             real gxIn, real gyIn, real gzIn,
             real oneOverdxIn, real oneOverdyIn, real oneOverdzIn):
    DivKgradECDM<BC,nv>(bcIn,nodeIn)
    {
      nxNodes=nxNodesIn;
      local_nxNodes=nodeIn.local_nxNodes;
      gx=gxIn;
      oneOverdx = oneOverdxIn;

      nyNodes=nyNodesIn;
      local_nyNodes=nodeIn.local_nyNodes;
      gy=gyIn;
      oneOverdy = oneOverdyIn;

      nzNodes=nzNodesIn;
      local_nzNodes=nodeIn.local_nzNodes;
      gz=gzIn;
      oneOverdz = oneOverdzIn;

      //what do I need to do here?
      hXcen.newsize(Vec::LOCAL,node.da_);
      hYcen.newsize(Vec::LOCAL,node.da_);
      hZcen.newsize(Vec::LOCAL,node.da_);
      
      hXcen = 1.0/oneOverdx;
      hYcen = 1.0/oneOverdy;
      hZcen = 1.0/oneOverdz;

      Dflux_x_center.resize(nv);
      Dflux_x_right.resize(nv);
      Dflux_y_center.resize(nv);
      Dflux_y_back.resize(nv);
      Dflux_z_center.resize(nv);
      Dflux_z_top.resize(nv);
      for (vi=0;vi<nv;vi++)
        {
          flux_x[vi].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes
			     *(local_nxNodes+1));
          flux_y[vi].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)
			     *local_nxNodes);
          flux_z[vi].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes
			     *local_nxNodes);
          flux_x[vi] = 0.0;
          flux_y[vi] = 0.0;
          flux_z[vi] = 0.0;
          Dflux_x_center[vi].resize(nv);
          Dflux_x_right[vi].resize(nv);
          Dflux_y_center[vi].resize(nv);
          Dflux_y_back[vi].resize(nv);
          Dflux_z_center[vi].resize(nv);
          Dflux_z_top[vi].resize(nv);
          for (vj=0;vj<nv;vj++)
            {
              Dflux_x_center[vi][vj].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
              Dflux_x_center[vi][vj] = 0.0;
              Dflux_x_right[vi][vj].newsize(Vec::LOCAL,local_nzNodes*local_nyNodes*(local_nxNodes+1));
              Dflux_x_right[vi][vj] = 0.0;

              Dflux_y_center[vi][vj].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
              Dflux_y_center[vi][vj] = 0.0;
              Dflux_y_back[vi][vj].newsize(Vec::LOCAL,local_nzNodes*(local_nyNodes+1)*local_nxNodes);
              Dflux_y_back[vi][vj] = 0.0;

              Dflux_z_center[vi][vj].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
              Dflux_z_center[vi][vj] = 0.0;
              Dflux_z_top[vi][vj].newsize(Vec::LOCAL,(local_nzNodes+1)*local_nyNodes*local_nxNodes);
              Dflux_z_top[vi][vj] = 0.0;
            }
	  pDiff_x[vi].newsize(Vec::LOCAL,node.ghost_nzNodes*node.ghost_nyNodes
			      *(node.ghost_nxNodes+1));
          pDiff_y[vi].newsize(Vec::LOCAL,
			      node.ghost_nzNodes*(node.ghost_nyNodes+1)
			      *node.ghost_nxNodes);
          pDiff_z[vi].newsize(Vec::LOCAL,
			      (node.ghost_nzNodes+1)*node.ghost_nyNodes
			      *node.ghost_nxNodes);
          pDiff_x[vi] = 0.0;
          pDiff_y[vi] = 0.0;
          pDiff_z[vi] = 0.0;

        }
    }
  virtual ~DivKgradECDM3d(){}
  //not correct, maybe works like standard fd still  
  virtual void computeDivergence(const Vec* K, const Vec* Rho, 
				 const Vec* P, 
                                 Vec* Div)
    {
      assert(0);
    }
  //not correct, maybe works like standard fd still  
  virtual void computeDivergenceJac(const Vec* K, const Vec* Rho, 
				    const Vec* P, 
                                    const VecVecVec& DK, const Vec* DRho, 
                                    VecVecVecVec& DivJac)
    {
      assert(0);
    }
};  


}//Daetk
#endif









