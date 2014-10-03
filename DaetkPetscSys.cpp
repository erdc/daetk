#include "DaetkPetscSys.h"
#include <cstdio>
namespace Daetk 
{
namespace Petsc
{
  namespace cc
  {
    extern "C"
    {
#ifndef DAETK_DEF_CPLUSPLUS_FOR_PETSC_H
#undef __cplusplus
#endif
#include "petsc.h"
#ifndef DAETK_DEF_CPLUSPLUS_FOR_PETSC_H
#define __cplusplus
#endif
    }
  }
  
bool Sys::master(){return !rank;}

bool Sys::initialized=false;

int Sys::rank=0;

int Sys::size=0;
  
Err::Err(){}

Err::Err(int i)
{
  //  CHKERRQ(i);
}

int Err::operator=(int i)
{
  using namespace cc;
  int line=0;
  if (i)
    PetscError(PETSC_COMM_SELF,line,"unknown Daetk function","unknown file","daetk source",i,PETSC_ERROR_INITIAL,"a funtion in the petsc library threw an error which the enclosing daetk function can't handle");
  return i;
}

Sys::Sys():
  commCreator(false)
{}

Sys::Sys(int& argc, char **argv,char* help, char* file)
{
  using namespace cc;
  using std::cerr;
  using std::endl;
  using std::flush;
  using std::cout;
  if (initialized)
    {
      PetscBool isInitialized;
      Err ierr=PetscInitialized(&isInitialized);
      assert(isInitialized);
      commCreator = false;
      ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
      ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    }
  else
    {
      if (!help)
	  help = (char*)(PETSC_NULL);
      if (!file)
	  file = (char*)(PETSC_NULL);
      commCreator=true;
      Err ierr;
      initialized=true;
      ierr = PetscInitialize(&argc,&argv,file,help);
      ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
      ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    }
}

Sys::~Sys()
{
  using namespace cc;
  if (commCreator)
    {
      Err ierr;
      PetscBool finalizedAlready(PETSC_FALSE);
      ierr = cc::PetscFinalized(&finalizedAlready);
      if (!finalizedAlready)
	ierr = cc::PetscFinalize();
      initialized=false;
    }
}

void Sys::barrier()
{
  using namespace cc;
  Err ierr;
  ierr = MPI_Barrier(PETSC_COMM_WORLD);
}


bool Sys::isInitialized()
{
  return initialized;
}

void Sys::beginSequential(int ng)
{
  using namespace cc;
  ierr =  PetscSequentialPhaseBegin(PETSC_COMM_WORLD,ng);
}

void Sys::endSequential(int ng)
{
  using namespace cc;
  ierr =  PetscSequentialPhaseEnd(PETSC_COMM_WORLD,ng);
}

bool Sys::catchError(bool error)
{
  using namespace cc;
  int thisVal=error;
  int result=error;
  MPI_Allreduce(&thisVal,&result,1,MPI_INT,MPI_LOR,PETSC_COMM_WORLD);
  return result;
}

int Sys::getRank(){return rank;}
int Sys::getSize(){return size;}

}//Petsc
}//Daetk
