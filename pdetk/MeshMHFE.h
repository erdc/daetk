#ifndef MESH_MHFE_H
#define MESH_MHFE_H

#include "ParameterDatabase.h"
#include "StencilMHFEre.h"
//#include "StencilMHFEre.h"
//just a simple class for keeping track of mesh coordinates
//for mhfe test problems, requires a Stencil though

//templatized on stencil but going to do 3d together with
//2d

class MeshMHFE
{
public:
  MeshMHFE(ParameterDatabase& pd);

  ~MeshMHFE() {}

  //set current cell to cellNum
  void globalCell(const Petsc::StencilMHFEreBase& stencil);
  //same as globalCell for now
  void localCell(const Petsc::StencilMHFEreBase& stencil)
    { globalCell(stencil); }
  //set face of current cell to im
  void setCurrentFace(int im);

  //data members
  //just basic mesh information
  int nCellsX,nCellsY,nCellsZ;
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

  
protected:
  //I should try and find a good place to keep this
  enum Point { pXminus = 0, 
	       pXplus  = 1,
	       pYminus = 2,
	       pYplus  = 3,
	       pZminus = 4,
	       pZplus  = 5};

};

inline
void MeshMHFE::globalCell(const Petsc::StencilMHFEreBase& stencil)
{
  //assumes that  stencil.globalCell(cellNumber); 
  //has been called

  Xc = xLeft  + 0.5*dX + dX*stencil.getCurrentCell_K();
  Yc = yFront + 0.5*dY + dY*stencil.getCurrentCell_J();
  Zc = zBottom+ 0.5*dZ + dZ*stencil.getCurrentCell_I();

  Xf = Xc;  Yf = Yc;  Zf = Zc;

}


//assumes that globalCell has been called
inline
void MeshMHFE::setCurrentFace(int p)
{
  assert(0<= p <= 5);
  switch (p)
    {
      case (pXminus):
	{
	  Xf = Xc - 0.5*dX; Yf = Yc; Zf = Zc;
	  break;
	}
    case (pXplus):
	{
	  Xf = Xc + 0.5*dX; Yf = Yc; Zf = Zc;
	  break;
	}
    case (pYminus):
	{
	  Yf = Yc - 0.5*dY; Xf = Xc; Zf = Zc;
	  break;
	}
    case (pYplus):
	{
	  Yf = Yc + 0.5*dY; Xf = Xc; Zf = Zc;
	  break;
	}
    case (pZminus):
	{
	  Zf = Zc - 0.5*dZ; Xf = Xc; Yf = Yc;
	  break;
	}
    case (pZplus):
	{
	  Zf = Zc + 0.5*dZ; Xf = Xc; Yf = Yc;
	  break;
	}
    default:
      {
	assert(0);
	break;
      }
    }//end switch?
}



#endif
