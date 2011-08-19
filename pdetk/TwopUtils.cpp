#include "TwopUtils.h"

namespace Daetk 
{
namespace TwoPhaseFlow
{
namespace Utilities
{

//a utilitity so that I don't have to reimplement
//this function all the time? 
void useScaleHeterogeneityUtil(Psk2p& psk,
			       ParameterDatabase& pd,
			       Petsc::SecondOrderFd* coarseNode, 
			       const Vec& deltaIn,
			       Petsc::SecondOrderFd& node, 
			       const Vec& vecExample,
			       const real& dx,
			       const real& dy,
			       const real& dz)
{
  int i,j,k,ci,cj,ck;
  real lx,ly,lz,DX=pd.r("scalDX"),DY=pd.r("scalDY"),DZ=pd.r("scalDZ");

  assert(coarseNode); 
  Vec delta(vecExample);

  for (i=node.local_z0;i<node.local_z0+node.local_nzNodes;i++)
    for (j=node.local_y0;j<node.local_y0+node.local_nyNodes;j++)
      for (k=node.local_x0;k<node.local_x0+node.local_nxNodes;k++)
        {
          lz = i*dz;
          ly = j*dy;
          lx = k*dx;
          ci = std::max(round(lz/DZ),coarseNode->local_z0);
          cj = std::max(round(ly/DY),coarseNode->local_y0);
          ck = std::max(round(lx/DX),coarseNode->local_x0);
          if (ci < coarseNode->local_z0 || 
	      ci >= coarseNode->local_z0+coarseNode->local_nzNodes)
            std::cerr<<"z dir"<<ci<<'\t'<< coarseNode->local_z0<<'\t'
		     <<coarseNode->local_z0+coarseNode->local_nzNodes
		     <<std::endl<<std::flush;
          if (cj < coarseNode->local_y0 || 
	      cj >= coarseNode->local_y0+coarseNode->local_nyNodes)
            std::cerr<<"y dir"<<cj<<'\t'<< coarseNode->local_y0<<'\t'
		     <<coarseNode->local_y0+coarseNode->local_nyNodes
		     <<std::endl<<std::flush;
          if (ck < coarseNode->local_x0 || 
	      ck >= coarseNode->local_x0+coarseNode->local_nxNodes)
            std::cerr<<"x dir"<<ck<<'\t'<< coarseNode->local_x0<<'\t'
		     <<coarseNode->local_x0+coarseNode->local_nxNodes
		     <<std::endl<<std::flush;
          //should probably do something about the case lz%DZ = DZ/2.0;
          (*coarseNode)(ci,cj,ck);
          node(i,j,k);
          delta(node.center) = deltaIn(coarseNode->center);
        }
  psk.millerSimilarScaling(delta);
  
}

}//Utilities
}//TwoPhaseFlow
}//Daetk
