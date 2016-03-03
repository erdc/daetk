#include "ParameterDatabase.h"
#include "mpi.h"

#include <cassert>

namespace Daetk 
{
namespace Petsc
{
  namespace cc
  {
#include "petsc.h"
  }
}
//mwf put in because my gcc doesn't seem to have this ios_base
#ifndef HAVE_OLD_STD_LIB
using std::ios_base;
#endif
using std::ios;
using std::istream;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::setw;
using std::flush;
using std::setiosflags;
//using std::
  
  using namespace Daetk::Petsc::cc;
bool& ParameterDatabase::b(std::string strKey) 
{ return (*this)(strKey,BOOL).b;}

int& ParameterDatabase::i(std::string strKey) 
{ return (*this)(strKey,INT).i;}

real& ParameterDatabase::r(std::string strKey)
{ return (*this)(strKey,REAL).r;}

IntVec& ParameterDatabase::iv(std::string strKey) 
{ return (*this)(strKey,INTVEC).iv;}

CMRVec<real>& ParameterDatabase::v(std::string strKey) 
{ return (*this)(strKey,VEC).v;}

std::string& ParameterDatabase::s(std::string strKey) 
{ return (*this)(strKey,PARAMSTRING).str;}

bool ParameterDatabase::quit(){return QUIT;}

MenuStatus  ParameterDatabase::menu(istream& fin)
{
  MenuStatus status;
  status.start     = false;
  status.quit      = false;
  status.exitLevel = false;
  if (pSys.master())
    {
      fin.setf(ios::skipws);
      while (!status.start && !status.quit)
	{
	  string name;
	  cout<<"Enter a parameter or (m)enu/(e)xit menu/(s)tart/(q)uit: "<<flush;
	  fin>>name;
	  if (!isalpha(name[0]))
	    fin.ignore(80,'\n');
	  else if (name == "start" || name == "s")
	    status.start=true;
	  else if (name == "quit" || name == "q")
	    {
	      status.quit = true;
	      QUIT = true;
	    }
	  else if (name == "menu" || name == "m")
	    {
	      cout<<endl<<endl;
	      writeAll(cout);
	      cout<<endl<<endl;
	    }
	  else if (name == "exit" || name == "e")
	    {
	      cout <<"exiting sub menu"<<endl;
	      status.quit = true;
	      QUIT        = true;
	      status.exitLevel = true;
	    }
	  else
	    {  
	      bool foundInLocalDB = true;
	      ParameterMap::iterator ip = pmap.find(name);
	      bool foundInSubDB   = false;
	      SubDatabaseMap::iterator 
		ispd = subParameterDatabase.find(name);
	      if (ip == pmap.end()) 
		{
		  foundInLocalDB = false;
		}
	      if (!foundInLocalDB)
		{
		  if (ispd != subParameterDatabase.end()) 
		    foundInSubDB = true;
		}
	      if (!foundInLocalDB && !foundInSubDB)
		{
		  cout<<name<<" not found"<<endl;
		}
	      if (foundInLocalDB)
		{
		  cout<<name<<" = "<<flush;
		  switch (ip->second.t)
		    {
		    case DEFAULT:
		      {
			fin>>ip->second.r;
			break;
		      }
		    case BOOL:
		      {
			fin>>ip->second.b;
			break;
		      }
		    case INT:
		      {
			fin>>ip->second.i;
			break;
		      }
		    case REAL:
		      {
			fin>>ip->second.r;
			break;
		      }
		    case INTVEC:
		      {
			fin>>ip->second.iv;
			break;
		      }
		    case VEC:
		      {
			fin>>ip->second.iv;
			break;
		      }
		    case PARAMSTRING:
		      {
			fin>>ip->second.iv;
			break;
		      }
		    }
		}//end localDB 
	      //don't search sub data base if found in local db 
	      else if(foundInSubDB)
		{
		  cout <<"calling menu for "<<name<<endl;
		  assert(ispd->second);
		  MenuStatus subStatus = ispd->second->menu(fin);
		  status = subStatus;
		  QUIT = status.quit;
		  if (subStatus.exitLevel == true)
		    {
		      status.quit     = false;
		      QUIT            = false;
		      status.exitLevel= false;
		    }
		}//end subDB
	    }//end else 
	}//end not quit
      QUIT = status.quit;
    }//end if master

  broadcast();

  Petsc::Err ierr;

  int  st=status.start,qu=status.quit,el=status.exitLevel;

  ierr = MPI_Bcast(&st,1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
  ierr = MPI_Bcast(&qu,1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
  ierr = MPI_Bcast(&el,1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);

  QUIT = qu;
  status.start=st;
  status.quit = qu;
  status.exitLevel = el;
  return status;
}

ostream&  operator<<(ostream& s,  ParameterDatabase& D)
{
  D.writeAll(s);
  return s;
}


ParameterDatabase::ParameterDatabase():
  QUIT(false),
#ifdef HAVE_OLD_STD_LIB
  missingOut("missingParameters.txt",ios::app)
#else
  missingOut("missingParameters.txt",ios_base::app)
#endif
{
  Tracer tr("ParameterDatabase::ParameterDatabase()");
}

ParameterDatabase::ParameterDatabase(istream& s):
  QUIT(false)
{
  Tracer tr("ParameterDatabase::ParameterDatabase()");

  if (pSys.master())
    readAll(s);
  broadcast();
}

ParameterDatabase::ParameterDatabase(const char* filename):
  QUIT(false),
#ifdef HAVE_OLD_STD_LIB
  missingOut("missingParameters.txt",ios::app)
#else
  missingOut("missingParameters.txt",ios_base::app)
#endif
{
  Tracer tr("ParameterDatabase::ParameterDatabase()");
  if (pSys.master())
    {
      ifstream temp(filename);
      readAll(temp);
    }   
  broadcast();
}

ParameterDatabase::~ParameterDatabase()
{ Tracer tr("ParameterDatabase::~ParameterDatabase()"); }

unsigned int ParameterDatabase::readAll(istream& s)
{
  string line,type,name;
  
  while (s) 
    {
      line = "";
      name = "";
      type = "";
      s>>type;
      Parameter tmp;
      //mwf try to not read comments?
      bool insertValue(true);
      //        if (istr && (type.size() > 0) && isalpha(type[0])) 
      if (s && type[0] == '/')
	{
	  s.ignore(256,'\n');
	}

      if (s && (type.size() > 0) && isalpha(type[0])) 
	{
          s>>name;
          if (type == "bool")
            {
              s>>tmp.b;
              tmp.t = BOOL;
	      insertValue = true;
            }
          else if (type == "int")
            {
              s>>tmp.i;
              tmp.t = INT;
	      insertValue = true;
            }
          else if (type == "real")
            {
              s>>tmp.r;
              tmp.t = REAL;
	      insertValue = true;
            }
          else if (type == "IntVec")
            {
              s>>tmp.iv;
              tmp.t = INTVEC;
 	      insertValue = true;
           }
          else if (type == "Vec")
            {
              s>>tmp.v;
              tmp.t = VEC;
	      insertValue = true;
            }
	  else if (type == "string")
            {
              s>>tmp.str;
              tmp.t = PARAMSTRING;
	      insertValue = true;
            }
          else 
            {
              cerr<<"Type "<<type<<" not recognized by ParamterMap"<<endl;
              cerr<<"Storing value as type real"<<endl;
              s>>tmp.r;
              tmp.t = REAL;
	      insertValue = false;
            }
	  if (insertValue)
	    pmap.insert(ParameterMap::value_type(name,tmp));
        }
    }
  return pmap.size();
}

ParameterDatabase&
ParameterDatabase::clear()
{
  pmap.clear();
  return *this;
}

ParameterDatabase&
ParameterDatabase::insert(string str,  Parameter& val )
{
  pmap.insert(ParameterMap::value_type(str,val));
  return *this;
}
//mwf
ParameterDatabase&
ParameterDatabase::insertParameterDatabase(string str,  
					   ParameterDatabase* sdb )
{
  subParameterDatabase.insert(SubDatabaseMap::value_type(str,sdb));
  return *this;
}
//mwf?
int 
ParameterDatabase::removeParameterDatabase(string str)
{
  int count = subParameterDatabase.erase(str);
  return count;
}

string typeString(ParameterDatabase::Type t)
{
  switch (t)
    {
    case ParameterDatabase::DEFAULT:
      return "default";
    case ParameterDatabase::BOOL:
      return "bool";
    case ParameterDatabase::INT:
      return "int";
    case ParameterDatabase::REAL:
      return "real";
    case ParameterDatabase::INTVEC:
      return "IntVec";
    case ParameterDatabase::VEC:
      return "Vec";
    case ParameterDatabase::PARAMSTRING:
      return "string";
    default:
      return "Missing";
    }
}
ParameterDatabase::Parameter&
ParameterDatabase::operator()(string strKey, Type t) 
{ 
   Parameter& p((*this)[strKey]);
  if (p.t != t && t!=DEFAULT)
    {
      cerr<<"ParameterDatabase: Requested parameter "<<strKey<<" of type "<<typeString(t)
          <<" is of type "<<typeString(p.t)<<endl
          <<"ParameterDatabase: Changing the type."<<endl
          <<"ParameterDatabase: The following appended to missingParameters.txt:"<<endl
          <<typeString(t)<<'\t'<<strKey<<'\t'<<0<<endl;
      missingOut<<typeString(t)<<'\t'<<strKey<<'\t'<<0<<endl;
      p.t = t;
    }
  return p;
}

ParameterDatabase::Parameter&
ParameterDatabase::operator[](string strKey) 
{ 
  //first search local database
  bool foundInLocalDB = true;
  ParameterMap::iterator ip = pmap.find(strKey);
  bool foundInSubDB   = false;
  if (ip == pmap.end()) 
    {
      foundInLocalDB = false;
    }
  if (!foundInLocalDB)
    {
      SubDatabaseMap::iterator 
	ispd = subParameterDatabase.begin();
      while (ispd != subParameterDatabase.end() && !foundInSubDB)
	{
	  assert(ispd->second);
	  bool subStatus = ispd->second->findInSubDatabase(strKey,ip);
	  if (subStatus == true)
	    {
	      foundInSubDB = true;
	    }
	  ispd++;
	}
    }
  bool success = foundInSubDB || foundInLocalDB;
  if (!success)
    {
      cerr<<"ParameterDatabase: "<<strKey<<" not found."<<endl
          <<"ParameterDatabase: Inserting "<<strKey
	  <<" into database (uninitialized)."<<endl
          <<"ParameterDatabase: You should fix the database"<<endl;
      Parameter tmp;
      pmap.insert(ParameterMap::value_type(strKey,tmp));
      ip = pmap.find(strKey);
      return ip->second;
      // exit(1);
    }
  else
    {
      return ip->second;
    }
}
bool ParameterDatabase::findInSubDatabase(string strKey,
ParameterMap::iterator& ip) 
{ 
  //first search local database
  bool foundInLocalDB = true;
  //ParameterMap::iterator ip = pmap.find(strKey);
  ip = pmap.find(strKey);
  bool foundInSubDB   = false;
  if (ip == pmap.end()) 
    {
      foundInLocalDB = false;
    }
  if (!foundInLocalDB)
    {
      SubDatabaseMap::iterator 
	ispd = subParameterDatabase.begin();
      while (ispd != subParameterDatabase.end() && !foundInSubDB)
	{
	  assert(ispd->second);
	  bool subStatus = ispd->second->findInSubDatabase(strKey,ip);
	  if (subStatus == true)
	    {
	      foundInSubDB = true;
	    }
	  ispd++;
	}
    }
  bool success = foundInSubDB || foundInLocalDB;
//    if (!success)
//      {
//        cerr<<"ParameterDatabase: "<<strKey<<" not found"<<endl;
//        exit(1);
//      }
//    else
//      {
//        //p = ip->second;
//      }
  return success;
}

ostream&
ParameterDatabase::writeAll(ostream& s) 
{
  //  if (pSys.master())
    {
  s<<setiosflags(ios::left);
  int col = 1;
  ParameterMap::iterator it = pmap.begin();
  while (it != pmap.end()) 
    {
      col = (++col)%2;
      s<<setw(15)<<it->first.c_str();
      switch (it->second.t)
        {
          case DEFAULT:
            {
              s<<setw(15)<<it->second.r;
              break;
            }
          case BOOL:
            {
              s<<setw(15)<<it->second.b;
              break;
              }
          case INT:
            {
              s<<setw(15)<<it->second.i;
              break;
            }
          case REAL:
            {
              s<<setw(15)<<it->second.r;
              break;
              }
          case INTVEC:
            {
              s<<setw(15)<<it->second.iv;
              break;
            }
          case VEC:
            {
              s<<setw(15)<<it->second.v;
              break;
            }
          case PARAMSTRING:
            {
              s<<it->second.str;
              break;
            }
        }
      if (col == 1)
        s<<endl;
      ++it;
    }
  s<<endl;
  //now write sub databases?
  col = 1;
  SubDatabaseMap::iterator 
    ispd = subParameterDatabase.begin();
  while (ispd != subParameterDatabase.end())
    {
      col = (++col)%2;
      s<<setw(15)<<"subDatabase  "
       <<setw(15)<<ispd->first.c_str() <<endl;
      ispd++;
      if (col == 1)
        s<<endl;
    }
  s<<endl;
    }
  return s;
}

void
ParameterDatabase::broadcast() 
{
  Petsc::Err ierr;
  bool master=pSys.master();
  int size;
  if (master)
    size=pmap.size();
  else
    {
      pmap.clear();
      size=0;
    }
  ierr = MPI_Bcast(&size,1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
  ParameterMap::iterator it = pmap.begin();
  for(int i=0;i<size;i++)
    {      
      Parameter tmp;

      int length=0;
      if (master)
        length=it->first.length();
      ierr = MPI_Bcast(&length,1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
      char* cstr = new char[length];
      if (master)
        it->first.copy(cstr,length);
      ierr = MPI_Bcast(cstr,length,MPI_CHAR,0,Petsc::cc::PETSC_COMM_WORLD);
      std::string name(cstr,length);
      delete [] cstr;
      if (master)
        assert(name == it->first);
      if (master)
        tmp.t=it->second.t;
      ierr = MPI_Bcast(&(tmp.t),1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
      
      switch (tmp.t)
        {
        case DEFAULT:
          {
            if (master)
              tmp.r=it->second.r;
            ierr = MPI_Bcast(&(tmp.r),1,MPI_DOUBLE,0,Petsc::cc::PETSC_COMM_WORLD);
            break;
          }
        case BOOL:
          {
            if (master)
              tmp.b=it->second.b;
            int ibool=tmp.b;
            ierr = MPI_Bcast(&ibool,1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
            tmp.b=ibool;
            break;
          }
        case INT:
          {
            if (master)
              tmp.i=it->second.i;
             ierr = MPI_Bcast(&(tmp.i),1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
            break;
          }
        case REAL:
          {
            if (master)
              tmp.r=it->second.r;
            ierr = MPI_Bcast(&(tmp.r),1,MPI_DOUBLE,0,Petsc::cc::PETSC_COMM_WORLD);
            break;
          }
        case INTVEC: //I don't do the right thing here, but can be fixed
          {
            int size;
            if (master)
              {
                tmp.iv=it->second.iv;
                size=tmp.iv.storageLength();
              }
            ierr = MPI_Bcast(&(size),1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
        
            if (!master)
              tmp.iv.newsize(size);
            
            ierr = MPI_Bcast(tmp.iv.castToArray(),size,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
  
            break;
          }
        case VEC:
          {
            int size;
            if (master)
              {
                tmp.v=it->second.v;
                size=tmp.v.storageLength();
              }
            ierr = MPI_Bcast(&(size),1,MPI_INT,0,Petsc::cc::PETSC_COMM_WORLD);
            
            if (!master)
              tmp.v.newsize(size);
            
            ierr = MPI_Bcast(tmp.v.castToArray(),size,MPI_DOUBLE,0,Petsc::cc::PETSC_COMM_WORLD);
  
            break;
          }
	case PARAMSTRING:
	  {
	    int valuelen = 0;
	    if (master)
	      valuelen = it->second.str.length();
	    //broadcast size
	    ierr = MPI_Bcast(&valuelen,1,MPI_INT,0,
			     Petsc::cc::PETSC_COMM_WORLD);
	    char * cvalue = new char[valuelen];
	    if (master)
	      it->second.str.copy(cvalue,valuelen);
	    
	    ierr =  MPI_Bcast(cvalue,valuelen,MPI_CHAR,0,
			      Petsc::cc::PETSC_COMM_WORLD);
	    std::string vname(cvalue,valuelen);
	    delete [] cvalue;
	    if (master)
	      assert(vname == it->second.str);

	    tmp.str = vname;

	    break;
	  }//end string
        }//end switch
      if (!master)
          pmap.insert(ParameterMap::value_type(name,tmp));
      else
        ++it;
    }

  SubDatabaseMap::iterator 
    ispd = subParameterDatabase.begin();
  while (ispd != subParameterDatabase.end())
    {
      ispd->second->broadcast();
      ++ispd;
    }
}

}//Daetk
