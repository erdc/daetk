#ifndef TWOP_DATA_H
#define TWOP_DATA_H

#include "Definitions.h"
#include "Utilities.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"
#include "ParameterDatabase.h"
#include "PetscSecondOrderFd.h"

namespace Daetk 
{
namespace TwoPhaseFlow
{

class TwopData
{
  /***********************************************************************
   put as much data as possible here, 
  **********************************************************************/
 public:

  enum phase {W,N};
  enum dimension {ONE_D,TWO_D,THREE_D};

  TwopData(Petsc::SecondOrderFd& s,
	   ParameterDatabase& pd, int dim);

  virtual ~TwopData();


  //--------------------------------------------------
  //data members
  //--------------------------------------------------
  Petsc::SecondOrderFd& node;

  dimension DIM;
  int nNodes, nxNodes, nyNodes, nzNodes;

  real dx, dy, dz,
    oneOverdx, oneOverdy, oneOverdz;

  real dirScale;
  //for keeping track of boundary condition changes
  real tForBCreset;

  //example vector for degrees of freedom holding global unknowns
  Vec ygdofExample;

  Vec m[2],Dm[2],
    theta[2],
    Dtheta[2],
    thetaS,
    p[2],
    Dp[2],
    resF[2],
    resM[2],
    div[2];

  Vec local_mCurrent[2],local_p[2],local_K[2],local_rho[2], 
    local_Drho[2], local_theta[2],
    local_DthetaW_DpC,local_DpC_DthetaW,
    local_thetaS,local_y,local_DDthetaW_DDpC,
    local_DDpC_DDthetaW,local_DDrho[2];

  VecIndex fwIndex,fnIndex,mwIndex,mnIndex;

  VecVecVec DK;
  VecVecVec local_DK;
  VecVecVecVec D_Div;

};//end TwopData

}//TwoPhaseFlow
}//Daetk

#endif
