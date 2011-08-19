#include "SecondOrderFd.h"

namespace Daetk
{
  SecondOrderFd::SecondOrderFd(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1):
    nxNodes(nxNodesIn),
    nyNodes(nyNodesIn),
    nzNodes(nzNodesIn),
    nxyNodes(nxNodesIn*nyNodesIn),
    nxyzNodes(nxyNodes*nzNodes),
//      pointMapField(nxyzNodes)
    pointListField(nxyzNodes)
  {}
  
  SecondOrderFd::~SecondOrderFd(){}

  void SecondOrderFd::setDomainSize(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1)
  {
    nxNodes = nxNodesIn;
    nyNodes = nyNodesIn;
    nzNodes = nzNodesIn;
    nxyNodes =nxNodesIn*nyNodesIn;
    nxyzNodes = nxyNodes*nzNodes;
  }
 
  SecondOrderFd1d::SecondOrderFd1d(int nxNodesIn, int nyNodes=1, int nzNodes =1):
    SecondOrderFd(nxNodesIn,1,1)
    {
      assert(nyNodes ==1);
      assert(nzNodes ==1);
      int i=0,j=0;
      for (int k=0;k<nxNodes;k++)
        {
          (*this)(k);
          pointList->push_back(Location(i,j,k,center,CENTER));
          if (k-1 >=0)
            pointList->push_back(Location(i,j,k-1,left,LEFT));
          if (k+1 < nxNodes)
            pointList->push_back(Location(i,j,k+1,right,RIGHT));
        }          
    }
  SecondOrderFd2d::SecondOrderFd2d(int nxNodesIn, int nyNodesIn, int nzNodes=1):
    SecondOrderFd(nxNodesIn,nyNodesIn,1)
    {
      assert(nzNodes==1);
      int i=0;
      for (int j=0;j<nyNodes;j++)
        for (int k=0;k<nxNodes;k++)
          {
            (*this)(j,k);
            pointList->push_back(Location(i,j,k,center,CENTER));
            if (k-1 >=0)
              pointList->push_back(Location(i,j,k-1,left,LEFT));
            if (k+1 < nxNodes)
              pointList->push_back(Location(i,j,k+1,right,RIGHT));
            if (j-1 >=0)
              pointList->push_back(Location(i,j-1,k,front,FRONT));
            if (j+1 < nyNodes)
              pointList->push_back(Location(i,j+1,k,back,BACK));
          }
    }
  SecondOrderFd3d::SecondOrderFd3d(int nxNodesIn, int nyNodesIn, int nzNodesIn):
    SecondOrderFd(nxNodesIn,nyNodesIn,nzNodesIn)
    {      
      for (int i=0;i<nzNodes;i++)
        for (int j=0;j<nyNodes;j++)
          for (int k=0;k<nxNodes;k++)
            {
              (*this)(i,j,k);
              pointList->push_back(Location(i,j,k,center,CENTER));
              if (k-1 >=0)
                pointList->push_back(Location(i,j,k-1,left,LEFT));
              if (k+1 < nxNodes)
                pointList->push_back(Location(i,j,k+1,right,RIGHT));
              if (j-1 >=0)
                pointList->push_back(Location(i,j-1,k,front,FRONT));
              if (j+1 < nyNodes)
                pointList->push_back(Location(i,j+1,k,back,BACK));
              if (i-1 >=0)
                pointList->push_back(Location(i-1,j,k,bottom,BOTTOM));
              if (i+1 < nzNodes)
                pointList->push_back(Location(i+1,j,k,top,TOP));
            }
    }

}//Daetk
