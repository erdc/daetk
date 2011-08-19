#ifndef TIMINGDATAFILE_H
#define TIMINGDATAFILE_H

#include "Definitions.h"
#include "DataFile.h"
#include "Chronograph.h"

namespace Daetk 
{
class TimingDataFile : public DataFile
{
public:
  TimingDataFile(const char* filename="data.txt");
  virtual ~TimingDataFile();
  virtual void startUserStep();
  virtual void endUserStep();
  virtual void includeSolution(const real& t,const Vec& s);
  virtual double& getRunTime();
protected:
  Chronograph clock;
  double runTime;
};
}//Daetk
#endif
