#include "TimingDataFile.h"

namespace Daetk
{
  TimingDataFile::TimingDataFile(const char* filename):DataFile(filename){}
   
  TimingDataFile::~TimingDataFile(){}
  
  void TimingDataFile::startUserStep(){ clock.start();}
   
  void TimingDataFile::endUserStep(){ clock.stop();}

  void TimingDataFile::includeSolution(const real& t,const Vec& s)
  {
    runTime = clock.elapsed();
    clock.reset();
    fout<<"Time for step = "<<runTime<<std::endl
        <<"Solution at t = "<<t<<std::endl
        <<s<<std::endl;
  }
   
  double& TimingDataFile::getRunTime(){ return runTime;}

}//Daetk
