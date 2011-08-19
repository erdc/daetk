#ifndef MESH_FD_H
#define MESH_FD_H

#include "ParameterDatabase.h"
#include "PetscSecondOrderFd.h"
#include "PetscStencilMM.h"
#include "IntVec.h"
//#include "Mapping.h"

//just a simple class for keeping track of mesh coordinates
//for mhfem test problems, requires a Stencil though

//templatized on stencil but going to do 3d together with
//2d
//#define IDENTITY_MAP
//#define ROTATION_MAP
//#define YOTOV_MAP

namespace Daetk
{

template <class MAPPING>
class MeshFD
{
public:
  MeshFD(ParameterDatabase& pd);

  ~MeshFD() {}

  //set current Node to NodeNum
  template <class STENCIL>
  void globalNode(const STENCIL& stencil);
  //
  template <class STENCIL>
  void localNode(const  STENCIL& stencil);

  template <class STENCIL>
  void localNode(const STENCIL& stencil,
		 int k);

  template <class STENCIL>
  void localNode(const STENCIL& stencil,
		 int j,int k);

  template <class STENCIL>
  void localNode(const STENCIL& stencil,
		 int i,int j,int k);

  void globalNode(int i,int j,int k);


  void setCurrentFace(int im);

  void setCurrentCorner(int im);

  //how about trying to print coordinates for now
  template <class STENCIL>
  void printGrid(STENCIL& stencil, bool printRef=true);

  //given global nodes, tells you if point
  //is dummy or not?
  inline bool isDummyNode(int j, int k);
  inline bool isDummyNode(int i,int j, int k);
  //data members
  inline bool isCellCenter(int i,int j, int k);  
  //
  bool usePointDistributedGrid;
  //
  bool useEnhancedCellCenteredGrid,useMM;
  //just basic mesh information
  int nxNodes,nyNodes,nzNodes;
  //starting points of mesh
  real xLeft,yFront,zBottom;
  //ending points of mesh
  real xRight,yBack,zTop;
  //
  real dX,dY,dZ;
  //coordinates of cell center
  real Xc,Yc,Zc;
  //coordinates of current Face
  real Xf,Yf,Zf;
  //stenicl anchor's current i,j,k
  int ic,jc,kc;
  //I don't know how this will work in 
  //parallel
  //is the current point on a face?
  IntVec isXface,isYface,isZface;
  //location of current point
  Vec Xvec,Yvec,Zvec;
  //this is the number of subdomains in each direction
  int NX,NY,NZ;
  //these are the number of cells 
  //per subdomain in each direction
  IntVec nxSub,nySub,nzSub;
  //number of cells in each direction
  int nxCells,nyCells,nzCells;
  //mwf these are for mapping reference grid
  //to physical domain. Call Map F
  // All other variables
  //are reference domain
  //mwf Two D only!
  //Fx(Xc,Yc),Fy(Xc,Yc)
  real XcP,YcP,ZcP;
  //Fx(Xf,Yf),Fy(Xf,Yf)
  real XfP,YfP,ZfP;
  //dF_x/dX^{hat},dF_x/dY^{hat} eval at (Xc,Yc)
  real dXdXcHat,dXdYcHat;
  //dF_y/dX^{hat},dF_y/dY^{hat} eval at (Xc,Yc)
  real dYdXcHat,dYdYcHat;
  //dF_x/dX^{hat},dF_x/dY^{hat} eval at (Xf,Yf)
  real dXdXfHat,dXdYfHat;
  //dF_y/dX^{hat},dF_y/dY^{hat} eval at (Xf,Yf)
  real dYdXfHat,dYdYfHat;
  //determinant of DF
  real jacC;
  //matrix entries for G(Xc,Yc)=1/jacC*DF
  //write DF as function of ref variables Xc,Yc
  real Gxx,Gxy,Gyx,Gyy;
  //Boundary integral terms for Neumann BC's
  //evaluated at current Xf,Yf
  //X same for left and right
  //Y same for front and back
  real sideJacXf,sideJacYf;
  //just for stupid irregular mesh definition
  real slopeConst,slopeDenom;
  //set up points on physical grid corresponding
  //to current Xc,Yc, and Xf,Yf
  inline void mapCoordinates();
  //set up points and derivatives on physical grid corresponding
  //to current Xc,Yc, and Xf,Yf
  inline void mapDerivatives();
  //map vector from reference to physical domain
  inline real mapVectorX(real Vx,real Vy);
  inline real mapVectorY(real Vx,real Vy);
  //map point from physical domain to referece domain
  inline real mapPointToRefX(real x, real y);
  inline real mapPointToRefY(real x, real y);

  MAPPING mapping;
  //these are the actual functions used 
  //for given application
  inline real Fx(real x,real y);
  inline real Fy(real x,real y);
  inline real dFxdX(real x, real y);
  inline real dFxdY(real x, real y);
  inline real dFydX(real x, real y);
  inline real dFydY(real x, real y);
  //inverse of F
  inline real FxInv(real x,real y);
  inline real FyInv(real x,real y);
protected:
  //I should try and find a good place to keep this
  enum Point { pXminus = 0, 
	       pXplus  = 1,
	       pYminus = 2,
	       pYplus  = 3,
	       pZminus = 4,
	       pZplus  = 5};

};

template <class MAPPING>
inline
void MeshFD<MAPPING>::globalNode(int i,int j,int k)
{
  kc = k;
  jc = j;
  ic = i;

  Xc = Xvec[k];
  Yc = Yvec[j];
  Zc = Zvec[i];

  Xf = Xc;  Yf = Yc;  Zf = Zc;

}

template <class MAPPING>
template <class STENCIL>
inline
void MeshFD<MAPPING>::globalNode(const STENCIL& stencil)
{
  //assumes that  stencil.globalNode(n); or (i,j,k)  
  //has been called
  
  ic = stencil.anchor->i;
  jc = stencil.anchor->j;
  kc = stencil.anchor->k;


  Xc = Xvec[kc];
  Yc = Yvec[jc];
  Zc = Zvec[ic];


  Xf = Xc;  Yf = Yc;  Zf = Zc;

}

template <class MAPPING>
template <class STENCIL>
inline
void MeshFD<MAPPING>::localNode(const STENCIL& stencil)
{
  //assumes that  stencil.localNode(n);
  //has been called
  
  ic = stencil.local_z0+stencil.anchor->i;
  jc = stencil.local_y0+stencil.anchor->j;
  kc = stencil.local_x0+stencil.anchor->k;


  Xc = Xvec[kc];
  Yc = Yvec[jc];
  Zc = Zvec[ic];

  Xf = Xc;  Yf = Yc;  Zf = Zc;

}

template <class MAPPING>
template <class STENCIL>
inline
void MeshFD<MAPPING>::localNode(const STENCIL& stencil,
		       int i,int j,int k)
{
  //assumes that  stencil.localNode(n);
  //has been called
  kc = stencil.local_x0+k;
  jc = stencil.local_y0+j;
  ic = stencil.local_z0+i;

  Xc = Xvec[kc];
  Yc = Yvec[jc];
  Zc = Zvec[ic];
  
  Xf = Xc;  Yf = Yc;  Zf = Zc;

}

template <class MAPPING>
template <class STENCIL>
inline
void MeshFD<MAPPING>::localNode(const STENCIL& stencil,
		       int j,int k)
{
  //assumes that  stencil.localNode(n);
  //has been called
  kc = stencil.local_x0+k;
  jc = stencil.local_y0+j;
  ic =0;

  Xc = Xvec[kc];
  Yc = Yvec[jc];
  Zc = Zvec[ic];

  Xf = Xc;  Yf = Yc;  Zf = Zc;

}

template <class MAPPING>
template <class STENCIL>
inline
void MeshFD<MAPPING>::localNode(const STENCIL& stencil,
		       int k)
{
  //assumes that  stencil.localNode(n);
  //has been called
  kc = stencil.local_x0+k;
  jc = 0;
  ic = 0;

  Xc = Xvec[kc];
  Yc = Yvec[jc];
  Zc = Zvec[ic];

  Xf = Xc;  Yf = Yc;  Zf = Zc;

}


//assumes that globalCell has been called
//mwf now see if FD classes need boundar fluxes
//mwf to be xLeft-0.5*dX actually
template <class MAPPING>
inline
void MeshFD<MAPPING>::setCurrentFace(int p)
{
  assert(0<= p <= 5);
  switch (p)
    {
      case (pXminus):
	{
	  Xf = Xc - 0.5*dX;  Yf = Yc; Zf = Zc;
	  if (isXface(kc)) 
	    {
	      //actually an x face
	      Xf = Xc;
	    }
	  
	  break;
	}
    case (pXplus):
	{
	  Xf = Xc + 0.5*dX; Yf = Yc; Zf = Zc;
	  if (isXface(kc)) 
	    {
	      //actually a face
	      Xf = Xc;
	    }

	  break;
	}
    case (pYminus):
	{
	  Yf = Yc - 0.5*dY; Xf = Xc; Zf = Zc;
	  if (isYface(jc)) 
	    {
	      //actually a face
	      Yf = Yc;
	    }
	  break;
	}
    case (pYplus):
	{
	  Yf = Yc + 0.5*dY; Xf = Xc; Zf = Zc;
	  if (isYface(jc)) 
	    {
	      //actually a face
	      Yf = Yc;
	    }
	  break;
	}
    case (pZminus):
	{
	  Zf = Zc - 0.5*dZ; Xf = Xc; Yf = Yc;
	  if (isZface(ic)) 
	    {
	      //actually a face
	      Zf = Zc;
	    }
	  break;
	}
    case (pZplus):
	{
	  Zf = Zc + 0.5*dZ; Xf = Xc; Yf = Yc;
	  if (isZface(ic)) 
	    {
	      //actually a face
	      Zf = Zc;
	    }
	  break;
	}
    default:
      {
	assert(0);
	break;
      }
    }//end switch?
}

//only do 2d for now
template <class MAPPING>
inline
void MeshFD<MAPPING>::setCurrentCorner(int p)
{
  assert(0<= p <= 4);
  switch (p)
    {
      case (0):
	{
	  //SW corner
	  Xf = Xc - 0.5*dX;  Yf = Yc-0.5*dY; Zf = Zc;
	  if (isXface(kc)) 
	    {
	      //actually an x face
	      Xf = Xc;
	    }
	  if (isYface(jc)) 
	    {
	      //actually an y face
	      Yf = Yc;
	    }
	  
	  break;
	}
    case (1):
	{
	  //SE corner
	  Xf = Xc + 0.5*dX; Yf = Yc-0.5*dY; Zf = Zc;
	  if (isXface(kc)) 
	    {
	      //actually a face
	      Xf = Xc;
	    }
	  if (isYface(jc)) 
	    {
	      //actually an y face
	      Yf = Yc;
	    }

	  break;
	}
    case (2):
	{
	  //NW corner
	  Yf = Yc + 0.5*dY; Xf = Xc-0.5*dX; Zf = Zc;
	  if (isYface(jc)) 
	    {
	      //actually a face
	      Yf = Yc;
	    }
	  if (isXface(kc)) 
	    {
	      //actually a face
	      Xf = Xc;
	    }
	  break;
	}
    case (3):
	{
	  //NE corner
	  Yf = Yc + 0.5*dY; Xf = Xc + 0.5*dX; Zf = Zc;
	  if (isYface(jc)) 
	    {
	      //actually a face
	      Yf = Yc;
	    }
	  if (isXface(kc)) 
	    {
	      //actually a face
	      Xf = Xc;
	    }
	  break;
	}
    default:
      {
	assert(0);
	break;
      }
    }//end switch?
}

//assumes that j,k are global
template <class MAPPING>
inline 
bool MeshFD<MAPPING>::isDummyNode(int i, int j, int k)
{
  bool dummy=false;
  
  if (useEnhancedCellCenteredGrid ) // || useMM)
    {
      if (i==0 || i==nzNodes-1)
	{
	  if (j==0 || j==nyNodes-1)
	    {
	      if (k==0 || k==nxNodes-1)
		dummy =  true;
	    }//j
	}//i
    }//end if on useEnhanced
  return dummy;
}

//assumes that j,k are global
template <class MAPPING>
inline 
bool MeshFD<MAPPING>::isDummyNode(int j, int k)
{
  bool dummy=false;
  
  if (useEnhancedCellCenteredGrid) // || useMM)
    {
      if (j==0 || j==nyNodes-1)
	{
	  if (k==0 || k==nxNodes-1)
	    dummy =  true;
	}
    }//end if on useEnhanced
  return dummy;
}

//assumes i,j,k global indeces
template <class MAPPING>
inline 
bool MeshFD<MAPPING>::isCellCenter(int i,int j, int k)
{
  //ignore i for now
  bool cellCenter=true;
  
  if (useEnhancedCellCenteredGrid) //  || useMM)
    {
      if (isXface(k) || isYface(j))
	cellCenter= false;

    }//end if on useEnhanced
  return cellCenter;
}


template <class MAPPING>
template <class STENCIL>
void MeshFD<MAPPING>::printGrid(STENCIL& stencil,bool printRefCoordinates)
{
  //make global vectors for x,y,z coordinates 
  Vec Xvec(nxNodes*nyNodes*nzNodes),Yvec(nxNodes*nyNodes*nzNodes),
    Zvec(nxNodes*nyNodes*nzNodes);
  //make global vectors for global I,J,K indeces
  
  Vec Kvec(nxNodes*nyNodes*nzNodes),Jvec(nxNodes*nyNodes*nzNodes),
    Ivec(nxNodes*nyNodes*nzNodes);

  //record which unknowns are actually faces
  Vec KfaceVec(nxNodes*nyNodes*nzNodes),JfaceVec(nxNodes*nyNodes*nzNodes),
    IfaceVec(nxNodes*nyNodes*nzNodes);

  //now go through nodes on this processor and load into vectors
  for (int i=stencil.local_z0;
       i<stencil.local_z0+stencil.local_nzNodes;i++)
    {
      for (int j=stencil.local_y0;
	   j<stencil.local_y0+stencil.local_nyNodes;j++)
	{
	  for (int k=stencil.local_x0;
	       k<stencil.local_x0+stencil.local_nxNodes;k++)
            {
	      stencil(i,j,k);
              globalNode(i,j,k);
	      if (printRefCoordinates)
		{
		  Zvec(stencil.center) = Zc;
		  Yvec(stencil.center) = Yc;
		  Xvec(stencil.center) = Xc;
		}
	      else
		{
		  mapCoordinates();
		  Zvec(stencil.center) = ZcP;
		  Yvec(stencil.center) = YcP;
		  Xvec(stencil.center) = XcP;
		}
	      Ivec(stencil.center) = i;
	      Jvec(stencil.center) = j;
	      Kvec(stencil.center) = k;

	      if (isXface(k))
		KfaceVec(stencil.center) = -1;
	      else
		KfaceVec(stencil.center) = k;
	      if (isYface(j))
		JfaceVec(stencil.center) = -1;
	      else
		JfaceVec(stencil.center) = j;
	      if (isZface(i))
		IfaceVec(stencil.center) = -1;
	      else
		IfaceVec(stencil.center) = i;

	    }//end loop on k
	}//end loop for j
    }//end loop for i
  using namespace std;
  ofstream xout("X.grf");
  xout.setf(ios::scientific);
  xout.precision(10);
  ofstream yout("Y.grf");
  yout.setf(ios::scientific);
  yout.precision(10);
  ofstream zout("Z.grf");
  zout.setf(ios::scientific);
  zout.precision(10);
  
  ofstream iout("I.grf");
  iout.setf(ios::scientific);
  iout.precision(10);
  ofstream jout("J.grf");
  jout.setf(ios::scientific);
  jout.precision(10);
  ofstream kout("K.grf");
  kout.setf(ios::scientific);
  kout.precision(10);
  
  ofstream ifout("Iface.grf");
  ifout.setf(ios::scientific);
  ifout.precision(10);
  ofstream jfout("Jface.grf");
  jfout.setf(ios::scientific);
  jfout.precision(10);
  ofstream kfout("Kface.grf");
  kfout.setf(ios::scientific);
  kfout.precision(10);
  
  xout <<Xvec;
  yout <<Yvec;
  zout <<Zvec;

  iout <<Ivec;
  jout <<Jvec;
  kout <<Kvec;

  ifout <<IfaceVec;
  jfout <<JfaceVec;
  kfout <<KfaceVec;

  //now print out subdomain information
  ofstream nzout("nzSub.grf");
  ofstream nyout("nySub.grf");
  ofstream nxout("nxSub.grf");
  
  nzout<<nzSub<<endl;
  nyout<<nySub<<endl;
  nxout<<nxSub<<endl;
 
  //what about trying to create similar vectors for flux variables as well?
  Vec XeVec,YeVec,ZeVec;
  XeVec.newsize(Vec::LOCAL,stencil.local_nyNodes*(stencil.local_nxNodes+1));
  YeVec.newsize(Vec::LOCAL,stencil.local_nxNodes*(stencil.local_nyNodes+1));
  //y coordinate at x faces
  Vec  XceVec,YceVec,ZceVec;
  YceVec.newsize(Vec::LOCAL,stencil.local_nyNodes*(stencil.local_nxNodes+1));
  //x coordinate at y faces
  XceVec.newsize(Vec::LOCAL,stencil.local_nxNodes*(stencil.local_nyNodes+1));

  Vec KeVec,JeVec,IeVec;
  KeVec.newsize(Vec::LOCAL,stencil.local_nyNodes*(stencil.local_nxNodes+1));
  JeVec.newsize(Vec::LOCAL,stencil.local_nxNodes*(stencil.local_nyNodes+1));

  //loop through including left/front/bottom?
  int iglb(0),jglb(0),kglb(0);
  int i(0),j(0),k(0);

  //x flux
  for (i=0; i < stencil.local_nzNodes; i++)
    {
      for (j=0; j < stencil.local_nyNodes; j++)
	{
	  //do left face too
	  k = 0;
	  stencil.localIndex(i,j,k);
	  iglb = i + stencil.local_z0;
	  jglb = j + stencil.local_y0;
	  kglb = k + stencil.local_x0;

	  KeVec[stencil.interLeft] = iglb*(stencil.nxyNodes+stencil.nyNodes)
	    + jglb*(stencil.nxNodes+1)+kglb;

	  globalNode(iglb,jglb,kglb);
	  //left face
	  setCurrentFace(0);
	  if (printRefCoordinates)
	    {
	      XeVec[stencil.interLeft] = Xf;
	      YceVec[stencil.interLeft]= Yf;
	    }
	  else
	    {
	      mapCoordinates();
	      XeVec[stencil.interLeft] = XfP;
	      YceVec[stencil.interLeft]= YfP;
	    }

	  for (k=0; k < stencil.local_nxNodes; k++)
	    {

	      stencil.localIndex(i,j,k);
	      iglb = i + stencil.local_z0;
	      jglb = j + stencil.local_y0;
	      kglb = k + stencil.local_x0;

	      KeVec[stencil.interRight] = iglb*(stencil.nxyNodes+stencil.nyNodes)
		+ jglb*(stencil.nxNodes+1)+kglb+1;

	      globalNode(iglb,jglb,kglb);
	      //right face
	      setCurrentFace(1);
	      if (printRefCoordinates)
		{
		  XeVec[stencil.interRight] = Xf;
		  YceVec[stencil.interRight]= Yf;
		}
	      else
		{
		  mapCoordinates();
		  XeVec[stencil.interRight] = XfP;
		  YceVec[stencil.interRight]= YfP;
		}


	    }
	}
    }//end loop for Qx
  //y flux
  for (i=0; i < stencil.local_nzNodes; i++)
    {
      for (k=0; k < stencil.local_nxNodes; k++)
	{
	  //do front face too
	  j = 0;
	  stencil.localIndex(i,j,k);
	  iglb = i + stencil.local_z0;
	  jglb = j + stencil.local_y0;
	  kglb = k + stencil.local_x0;

	  JeVec[stencil.interFront] = iglb*(stencil.nxyNodes+stencil.nxNodes)
	    + jglb*(stencil.nxNodes)+kglb;


	  globalNode(iglb,jglb,kglb);
	  //front face
	  setCurrentFace(2);
	  if (printRefCoordinates)
	    {
	      YeVec[stencil.interFront] = Yf;
	      XceVec[stencil.interFront] = Xf;
	    }
	  else
	    {
	      mapCoordinates();
	      YeVec[stencil.interFront] = YfP;
	      XceVec[stencil.interFront] = XfP;
	    }

	  for (j=0; j < stencil.local_nyNodes; j++)
	    {

	      stencil.localIndex(i,j,k);
	      iglb = i + stencil.local_z0;
	      jglb = j + stencil.local_y0;
	      kglb = k + stencil.local_x0;

	      JeVec[stencil.interBack] = iglb*(stencil.nxyNodes+stencil.nxNodes)
		+ (jglb+1)*(stencil.nxNodes)+kglb;
  

	      globalNode(iglb,jglb,kglb);
	      //back face
	      setCurrentFace(3);
	      if (printRefCoordinates)
		{
		  YeVec[stencil.interBack] = Yf;
		  XceVec[stencil.interBack]= Xf;
		}
	      else
		{
		  mapCoordinates();
		  YeVec[stencil.interBack] = YfP;
		  XceVec[stencil.interBack]= XfP;
		}

	    }
	}
    }//end loop for Qy

  //print out x indeces
  ofstream keout("Ke.grf");
  keout.setf(ios::scientific);
  keout.precision(10);

  keout <<KeVec;

  //print out y indeces
  ofstream jeout("Je.grf");
  jeout.setf(ios::scientific);
  jeout.precision(10);

  jeout <<JeVec;

  //print out x face coords
  ofstream xeout("Xe.grf");
  xeout.setf(ios::scientific);
  xeout.precision(10);

  xeout <<XeVec;

  //print out y coords at x face
  ofstream yceeout("Yce.grf");
  yceeout.setf(ios::scientific);
  yceeout.precision(10);

  yceeout <<YceVec;

  //print out yf coords
  ofstream yeout("Ye.grf");
  yeout.setf(ios::scientific);
  yeout.precision(10);

  yeout <<YeVec;

  //print out x coords at y face
  ofstream xceeout("Xce.grf");
  xceeout.setf(ios::scientific);
  xceeout.precision(10);

  xceeout <<XceVec;

}//end routine

///////////////////////// for mapping to physical domain
//assumes globalNode or something similar has been called
template <class MAPPING>
inline 
void MeshFD<MAPPING>::mapCoordinates()
{
  XcP = mapping.Fx(Xc,Yc);  YcP = mapping.Fy(Xc,Yc); 
  XfP = mapping.Fx(Xf,Yf);  YfP = mapping.Fy(Xf,Yf);
  //don't do 3d for now
  ZcP = Zc;   ZfP = Zf;
}

template <class MAPPING>
inline 
void MeshFD<MAPPING>::mapDerivatives()
{
  mapCoordinates();
  //calculate first order partial derivs of F
  //at cell center
  dXdXcHat = mapping.dFxdX(Xc,Yc);  dXdYcHat = mapping.dFxdY(Xc,Yc); 
  dYdXcHat = mapping.dFydX(Xc,Yc);  dYdYcHat = mapping.dFydY(Xc,Yc); 
  //at current cell face
  dXdXfHat = mapping.dFxdX(Xf,Yf);  dXdYfHat = mapping.dFxdY(Xf,Yf); 
  dYdXfHat = mapping.dFydX(Xf,Yf);  dYdYfHat = mapping.dFydY(Xf,Yf); 
  
  //determinant of DF
  jacC = dXdXcHat*dYdYcHat - dYdXcHat*dXdYcHat;
  assert(fabs(jacC) > 0.0);

  //piola transform for mapping vector valued functions
  //evaluated at center
  Gxx = dXdXcHat/jacC;  Gxy = dXdYcHat/jacC;
  Gyx = dYdXcHat/jacC;  Gyy = dYdYcHat/jacC;

  //norms of "normal part" of DF for doing boundary integral
  //calculations. this is translation of arclength parameter 
  //to reference coordinates
  sideJacXf = sqrt(dXdYfHat*dXdYfHat + dYdYfHat*dYdYfHat);
  sideJacYf = sqrt(dYdXfHat*dYdXfHat + dXdXfHat*dXdXfHat);

  //mwf debug
  //sideJacXf = 1.0;
  //sideJacYf = 1.0;

}
//calculate x coordinate of G*V
template <class MAPPING>
inline 
real MeshFD<MAPPING>::mapVectorX(real Vx, real Vy)
{
  return Gxx*Vx + Gxy*Vy;
}

//calculate y coordinate of G*V
template <class MAPPING>
inline 
real MeshFD<MAPPING>::mapVectorY(real Vx, real Vy)
{
  return Gyx*Vx + Gyy*Vy;
}
template <class MAPPING>
inline 
real MeshFD<MAPPING>::mapPointToRefX(real x, real y)
{
  return mapping.FxInv(x,y);
}

//calculate y coordinate of G*V
template <class MAPPING>
inline 
real MeshFD<MAPPING>::mapPointToRefY(real x, real y)
{
  return mapping.FyInv(x,y);
}



using std::max;
using std::min;

template <class MAPPING>
MeshFD<MAPPING>::MeshFD(ParameterDatabase& pd):
  usePointDistributedGrid(pd.b("usePointDistributedGrid")),
  useEnhancedCellCenteredGrid(pd.b("useECDM")),
  useMM(pd.b("useMM")),
  nxNodes(pd.i("nxNodes")),
  nyNodes(pd.i("nyNodes")),
  nzNodes(pd.i("nzNodes")),
  xLeft(pd.r("xLeft")),
  yFront(pd.r("yFront")),
  zBottom(pd.r("zBottom")),
  xRight(pd.r("xRight")),
  yBack(pd.r("yBack")),
  zTop(pd.r("zTop")),
  dX(1.0),dY(1.0),dZ(1.0),
  Xc(0.0),Yc(0.0),Zc(0.0),
  Xf(0.0),Yf(0.0),Zf(0.0),
  ic(0),jc(0),kc(0),
  isXface(pd.i("nxNodes"),12345),
  isYface(pd.i("nyNodes"),12345),
  isZface(pd.i("nzNodes"),12345),
  Xvec(),
  Yvec(),
  Zvec(),
  NX(pd.i("NX")),
  NY(pd.i("NY")),
  NZ(pd.i("NZ")),
  nxSub(pd.i("NX"),1),
  nySub(pd.i("NY"),1),
  nzSub(pd.i("NZ"),1),
  nxCells(pd.i("nxNodes")),
  nyCells(pd.i("nyNodes")),
  nzCells(pd.i("nzNodes")),
  mapping(pd)
{
  Petsc::Sys psys;

  assert(nxNodes >= 1);
  assert(nyNodes >= 1);
  assert(nzNodes >= 1);

  dX = dY = dZ = 1.0;
  //mwf go ahead and set base spatial increments for
  //enhanced cell centered differences
  if (nxNodes > 1)
    {
      if (usePointDistributedGrid)
	dX = (xRight-xLeft)/(nxNodes-1);
      else
	dX = (xRight-xLeft)/(nxNodes);
    }
  if (nyNodes > 1)
    {
      if (usePointDistributedGrid)
	dY = (yBack-yFront)/(nyNodes-1);
      else
	dY = (yBack-yFront)/(nyNodes);
    }
  if (nzNodes > 1)
    {
      if (usePointDistributedGrid)
	dZ = (zTop-zBottom)/(nzNodes-1);
      else
	dZ = (zTop-zBottom)/(nzNodes);
    }
  if (!usePointDistributedGrid 
      && !useEnhancedCellCenteredGrid) // && !useMM)
    {
      xLeft  += 0.5*dX;
      yFront += 0.5*dY;
      zBottom+= 0.5*dZ;
    }



  //adjust number of cells for useEnhancedCellCenteredGrid
  //or useMM
  if (useEnhancedCellCenteredGrid) // || useMM)
    {
      nxCells = nxNodes-2 - (NX-1);
      if (nyNodes > 2)
	nyCells = nyNodes-2 - (NY-1);
      else
	nyCells = 1;
      if (nzNodes > 2)
	nzCells = nzNodes-2 - (NZ-1);
      else
	nzCells = 1;
      //account for maybe only in 2d 
      nxCells = max(1,nxCells);
      nyCells = max(1,nyCells);
      nzCells = max(1,nzCells);

      dX = (xRight-xLeft)/(nxCells);
      if (nyCells > 1)
	dY = (yBack-yFront)/(nyCells);
      else
	dY = 1.0;
      if (nzCells > 1)
	dZ = (zTop-zBottom)/(nzCells);
      else
	dZ = 1.0;
    }
  //set up subdomain information?
  if (pd.b("setSubDomainSizeDirectly"))
    {
       ParameterDatabase pds("subDomains.txt");

      IntVec nxSubTmp(NX);
      IntVec nySubTmp(NY);

	   
      nxSubTmp = pds.iv("nxSubDomain");
      nySubTmp = pds.iv("nySubDomain");
	
      assert(nxSubTmp.size()==NX);
      assert(nySubTmp.size()==NY);

      for (int m=0; m < NX ; m++)
	nxSub(m) = nxSubTmp(m);
      for (int m=0; m < NY ; m++)
	nySub(m) = nySubTmp(m);
    
      int nSum(0);
      for (int m=0; m < NX ; m++)
	nSum += nxSub(m);

      assert(nSum == nxCells);
      
      nSum = 0;

      for (int m=0; m < NY ; m++)
	nSum += nySub(m);

      assert(nSum == nyCells);

    }
  else
    {
      for (int ns=0; ns < NX; ns++)
	nxSub(ns) = nxCells/NX;
      for (int ns=0; ns < NY; ns++)
	nySub(ns) = nyCells/NY;
      for (int ns=0; ns < NZ; ns++)
	nzSub(ns) = nzCells/NZ;
      
      //last subdomain in each direction might get extra cells
      nxSub(NX-1) += nxCells % NX;
      nySub(NY-1) += nyCells % NY;
      nzSub(NZ-1) += nzCells % NZ;
    }
  //these should be local vectors of global length
  Xvec.newsize(Vec::LOCAL,nxNodes);
  Yvec.newsize(Vec::LOCAL,nyNodes);
  Zvec.newsize(Vec::LOCAL,nzNodes);

  //I'll have to set the face information from the calling class
  //set up a default configuration
  if (!useEnhancedCellCenteredGrid) // && !useMM)
    {
      for (int k=0; k < nxNodes; k++)
	{
	  isXface(k) = 0;
	  Xvec[k] = xLeft + k*dX;
	}
      for (int j=0; j < nyNodes; j++)
	{
	  isYface(j) = 0;
	  Yvec[j] = yFront + j*dY;
	}
      for (int i=0; i < nzNodes; i++)
	{
	  isZface(i) = 0;
	  Zvec[i]  = zBottom + i*dZ;
	}
    }//end !useEnhanced && !useMM
  else
    {
      //global index 
      int ig(-1);
      //global location
      real xLoc(xLeft);
      //go through sub domains in x direction
      for (int ns=0; ns < NX; ns++)
	{
	  int nc = nxSub(ns);
	  //should be first unknown in subdomain
	  ig++;
	  Xvec[ig] = xLoc;
	  isXface(ig)=1;
	  //try and count through the cells in the
	  //current direction
	  for (int ic = 0; ic < nc; ic++)
	    {
	      ig++;
	      xLoc += dX;
	      isXface(ig)=0;
	      Xvec[ig] = xLoc-0.5*dX;
	    }
	}//end loop through subdomains
      //now check to see that ig is nxNodes-2
      assert(ig==nxNodes-2);
      ig++;
      Xvec[ig] = xRight;
      isXface(ig)=1;


      //now do y direction
      if (nyNodes < 3)
	{
	  Yvec[0] = yFront;
	  isYface(0)=1;
	  Yvec[nyNodes-1]=yBack;
	  isYface(nyNodes-1)=1;
	}
      else
	{
	  ig = -1;
	  real yLoc(yFront);

	  //go through sub domains in y direction
	  for (int ns=0; ns < NY; ns++)
	    {
	      int nc = nySub(ns);
	      //should be first unknown in subdomain
	      ig++;
	      Yvec[ig] = yLoc;
	      isYface(ig)=1;
	      //try and count through the cells in the
	      //current direction
	      for (int ic = 0; ic < nc; ic++)
		{
		  ig++;
		  yLoc += dY;
		  isYface(ig)=0;
		  Yvec[ig] = yLoc-0.5*dY;
		}
	    }//end loop through subdomains
	  //now check to see that ig is nxNodes-2
	  assert(ig==nyNodes-2);
	  ig++;
	  Yvec[ig] = yBack;
	  isYface(ig)=1;
	  
	}//end else on nyNodes >= 3

      //now do y direction
      if (nzNodes < 3)
	{
	  Zvec[0] = zBottom;
	  isZface(0)=1;
	  Zvec[nzNodes-1]=zTop;
	  isZface(nzNodes-1)=1;
	}
      else
	{
	  ig = -1;
	  real zLoc(zBottom);

	  //go through sub domains in y direction
	  for (int ns=0; ns < NZ; ns++)
	    {
	      int nc = nzSub(ns);
	      //should be first unknown in subdomain
	      ig++;
	      Zvec[ig] = zLoc;
	      isZface(ig)=1;
	      //try and count through the cells in the
	      //current direction
	      for (int ic = 0; ic < nc; ic++)
		{
		  ig++;
		  zLoc += dZ;
		  isZface(ig)=0;
		  Zvec[ig] = zLoc-0.5*dZ;
		}
	    }//end loop through subdomains
	  //now check to see that ig is nxNodes-2
	  assert(ig==nzNodes-2);
	  ig++;
	  Zvec[ig] = zTop;
	  isZface(ig)=1;
	  
	}//end else on nyNodes > 3

    }//end else on useEnhanced or useMM

}


}//Daetk
#endif

