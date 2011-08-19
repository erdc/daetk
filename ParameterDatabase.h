#ifndef PARAMETERDATABASE_H
#define PARAMETERDATABASE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>

#ifndef HAVE_OLD_STD_LIB
#include <sstream>
#else
#include <strstream>
#endif

#include <vector>
#include <map>

#include "Definitions.h"
#include "Utilities.h"
#include "IntVec.h"
#include "Vec.h"
#include "DaetkPetscSys.h"

namespace Daetk 
{
//mwf
class MenuStatus
{
public:
  MenuStatus() : start(false),quit(false),exitLevel(false) {}
  ~MenuStatus() {}
  //
  bool start,quit,exitLevel;
};



class ParameterDatabase
{
public:
  
  enum Type {DEFAULT, BOOL, INT, REAL, INTVEC, VEC, PARAMSTRING};
  struct Parameter 
  {
    Type t;
    union 
    {
      bool b;
      int  i;
      real r;
    };
    IntVec iv;
    CMRVec<real> v;
    std::string str;
  };

  ParameterDatabase();
  ParameterDatabase(std::istream& s);
  //mwf see if I can make this const
  ParameterDatabase(const  char* filename);
  virtual ~ParameterDatabase();
  ParameterDatabase& clear();
  unsigned int readAll(std::istream& s);
  std::ostream& writeAll(std::ostream& s);
  ParameterDatabase& insert(std::string str, Parameter& v);
  Parameter& operator()(std::string strKey, Type t=DEFAULT);
  Parameter& operator[](std::string strKey);
  bool& b(std::string strKey);
  int& i(std::string strKey);
  real& r(std::string strKey);
  IntVec& iv(std::string strKey);
  CMRVec<real>& v(std::string strKey);
  std::string& s(std::string strKey);
  MenuStatus menu(std::istream& fin);
  bool quit();
  //mwf try to insert sub database?
  ParameterDatabase& insertParameterDatabase(std::string str,
					     ParameterDatabase* sdb);
  //mwf try to remove sub database?
  int removeParameterDatabase(std::string str);
  void broadcast();
protected:
  bool QUIT;
  std::ofstream missingOut;
  typedef std::map<std::string,Parameter> ParameterMap ;
  ParameterMap pmap;
  bool  findInSubDatabase(std::string strKey,ParameterMap::iterator& ip);

  //mwf
  //try to allow adding "sub databases to parameter database"
  typedef std::map<std::string,ParameterDatabase*> SubDatabaseMap;
  //does not do any memory allocation for now
  SubDatabaseMap subParameterDatabase;
  Petsc::Sys pSys;
};
}//Daetk
#endif
