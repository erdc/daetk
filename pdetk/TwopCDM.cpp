#include "TwopCDM.h"

namespace Daetk 
{
namespace TwoPhaseFlow 
{

TwopDaeDef::TwopDaeDef(Petsc::SecondOrderFd& s,
		       ParameterDatabase& pd, 
		       DataCollector& dataIn, 
		       int dim):
  TwopData(s,pd,dim),DaeDefinition(dim,&dataIn),t0(0.0),myY(0),
  twopOut("twop.out"),psys()
{
  //examples should be set in TwopData
  y0.newsize(dim);
  y0prime.newsize(dim);

}

TwopDaeDef::~TwopDaeDef()
{}

}//TwoPhaseFlow
}//Daetk
