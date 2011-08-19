#include "DataFile.h"

namespace Daetk 
{

DataFile::DataFile(const char* filename):
  DataCollector(),
  fout(filename)
{
  Tracer tr("DataFile::DataFile()");
  //mwf debug
  //std::cout<<"In DataFile ctor filename= "<<filename<<std::endl;
  //fout<<"FOOOOOOO"<<std::endl;
  //fout.flush();
}

DataFile::~DataFile()
{
  Tracer tr("DataFile::~DataFile()");
  //mwf debug
  //std::cout<<"In DataFile dtor "<<std::endl;
}

void DataFile::setOutputFile(const char* filename)
{
  fout.close();
  fout.open(filename); 
}
 
}//Daetk
