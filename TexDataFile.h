#include "FullDataFile.h"

namespace Daetk
{
  class TexDataFile : public FullDataFile
  {
  public:
    TexDataFile(const real& t0=0.0,const char* filename="data.tex");
    virtual ~TexDataFile();
    virtual void stepTaken(int k,real h, real tn, real errorEstimate=0);
    virtual void includeSolution(const real& t,const Vec& y);
  protected:
    std::ofstream texSummaryOut,matlabOut;
  };

}//namespace Daetk
