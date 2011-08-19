#include "DtHistoryDataFile.h"

namespace Daetk
{

DtHistoryDataFile::DtHistoryDataFile(const real& t0, const char* filename):
  TexDataFile(t0,filename),
  dtOut("dtSim.grf"),orderOut("orderSim.grf")
{


}

DtHistoryDataFile::~DtHistoryDataFile()
{



}

void DtHistoryDataFile::stepTaken(int k,real h, real tn, 
				  real errorEstimate)
{
  TexDataFile::stepTaken(k,h,tn,errorEstimate);

  dtOut<<tn<<"\t"<<h<<std::endl;
  orderOut<<tn<<"\t"<<k<<std::endl;
}

}//end Daetk
