#ifndef DATAFILE_H
#define DATAFILE_H

#include <fstream>
#include "Definitions.h"
#include "Utilities.h"
#include "Vec.h"
#include "DataCollector.h"


namespace Daetk 
{
class DataFile : public DataCollector
{
 public:
  DataFile(const char* filename="data.txt");
  virtual ~DataFile();
  void setOutputFile(const char* filename);
protected:
  std::ofstream fout;
};
}//Daetk
#endif
