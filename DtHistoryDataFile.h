#ifndef DT_HISTORY_DATAFILE_H
#define DT_HISTORY_DATAFILE_H

#include "TexDataFile.h"

namespace Daetk
{

class DtHistoryDataFile: public TexDataFile
{
  /***********************************************************************
    record all of the things in TexDataFile (and hence FullDataFile)
    but also record evolution of dt and order over simulation.

    writes file dtSim.grf that records simulation's time steps
       t^1    \Delta t^1
       ...
       t^n    \Delta t^n

    and ordSim.grf that records simulation's temporal approx. order
       t^1      k^1
       ...
       t^n      k^n
   
  **********************************************************************/
public:

  DtHistoryDataFile(const real& t0, const char* filename="data.tex");
  virtual ~DtHistoryDataFile();

  virtual void stepTaken(int k,real h, real tn, real errorEstimate=0);
protected:
  std::ofstream dtOut,orderOut;

};//DtHistoryDataFile

}//end Daetk


#endif
