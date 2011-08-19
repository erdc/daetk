#include "TestPostProc.h"

namespace Daetk
{
namespace Petsc
{
  namespace cc
  {
    extern "C"
    {
#include "petsc.h"
#include "petscvec.h"
#include "petscda.h"
    }
  }
  //}

  //using namespace Petsc;
using namespace Petsc::cc;
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::istream;

//  #define TRY_NEW_GLOBAL_VEC 
//  //try to do this like chris' << routine and ex8 in vec/examples
//  void Petsc::postProcessFluxesParallel(Petsc::SecondOrderFd& stencil,
//  			      Vec& Qx,Vec& Qy,Vec& Qz)
//  {

//    Err ierr;
//    Petsc::Sys pSys;
 
//    int rank(0),size(0);
//    MPI_Status        status;
//    //get number of processors and local processor number
//    rank = pSys.getRank();
//    size = pSys.getSize();

//    //do Qx first
//    int local_nxNodes = stencil.local_nxNodes;
//    int local_nyNodes = stencil.local_nyNodes;
//    int local_nzNodes = stencil.local_nzNodes;

//  #ifndef TRY_NEW_GLOBAL_VEC  
//    PetscViewer viewer;

//    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"qxTest.grf",&viewer);
//    ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX);
//    ierr = VecView(Qx.castToPetsc(),viewer);
//    ierr = PetscViewerDestroy(viewer);
//  #else  
//    Petsc::cc::_p_Vec*  pQx_;
//    int local_nQx = (local_nxNodes+1)*(local_nyNodes);
//    ierr = VecCreate(PETSC_COMM_WORLD,local_nQx,PETSC_DETERMINE,&pQx_);
//    ierr = VecSetFromOptions(pQx_);

 
//    //create set of global indices corresponding to this processors local
//    //nodes

//    int * globalIndicesForLocalQx = new int[local_nQx];
//    int * globalIndicesCheck= new int[local_nQx];
//    int * indicesForLocalQx = new int[local_nQx];

//    //also get array holding local values
//    real * localQxValues = new real[local_nQx];

//    int viewRank = 1;
//    if (rank == viewRank)
//      cout <<"in postProcessFluxesParallel on proc "<<rank<<endl;

//    //now loop through local indeces and set local flux vector using Qx
//    for (int j=0;j<local_nyNodes;j++)
//      for (int k=0; k<local_nxNodes;k++)
//        {
//  	stencil.localIndex(j,k);
//  	int interRight = stencil.interRight;
//  	localQxValues[interRight] = Qx[interRight];
//  	indicesForLocalQx[interRight] = interRight;
//  	//need to put in global flux index into stencil
//  	int jglb = j + stencil.local_y0;
//  	int kglb = k + stencil.local_x0;
//  	stencil(jglb,kglb);
//  	globalIndicesForLocalQx[interRight] = jglb*(stencil.nxNodes+1) + kglb + 1;

//  	if (rank == viewRank)
//  	  {
//  	    //cout<<"(j,k)= "<<j<<","<<k<<" (jglb,kglb)= "<<jglb<<","<<kglb;
//  	    cout<<endl<<" indicesForLocalQx["<<interRight<<"]= "
//  		<<indicesForLocalQx[interRight];
//  	    cout<<" localQxValues["<<interRight<<"]= "
//  		<<localQxValues[interRight];
//  	    cout<<" globalIndicesForLocalQx["<<interRight<<"]= "
//  		<<globalIndicesForLocalQx[interRight];
//  	    cout<<" petscToGlobal("<<stencil.center<<")= "
//  		<<stencil.petscToGlobal(stencil.center);
//  	  }
	    
//        }
  
//    int k = 0;
//    for (int j=0;j<local_nyNodes;j++)
//      {
//        stencil.localIndex(j,k);
//        int interLeft = stencil.interLeft;
//        if (stencil.local_x0 > 0)
//  	{
//  	  localQxValues[interLeft] = -123456.7;
//  	  indicesForLocalQx[interLeft] = -1;
//  	  globalIndicesForLocalQx[interLeft] = -1;
//  	}
//        else
//  	{
//  	  localQxValues[interLeft] = Qx[interLeft];
//  	  indicesForLocalQx[interLeft] = interLeft;
//  	  //need to put in global flux index into stencil
//  	  int jglb = j + stencil.local_y0;
//  	  int kglb = k + stencil.local_x0;
//  	  stencil(jglb,kglb);
//  	  globalIndicesForLocalQx[interLeft] = jglb*(stencil.nxNodes+1) + kglb;

//  	  if (rank == viewRank)
//  	    {
//  	      //cout<<"(j,k)= "<<j<<","<<k<<" (jglb,kglb)= "<<jglb<<","<<kglb;
//  	      cout<<endl<<" indicesForLocalQx["<<interLeft<<"]= "
//  		  <<indicesForLocalQx[interLeft];
//  	      cout<<" localQxValues["<<interLeft<<"]= "
//  		<<localQxValues[interLeft];
//  	      cout<<" globalIndicesForLocalQx["<<interLeft<<"]= "
//  		<<globalIndicesForLocalQx[interLeft];
//  	      cout<<" petscToGlobal("<<stencil.center<<")= "
//  		  <<stencil.petscToGlobal(stencil.center);
//  	    }

//  	}      
//      }
//    //create mapping from these local indeces to the global indeces

//    {//local scope
//      ISLocalToGlobalMapping ltog;
//      ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_SELF,local_nQx,
//  					globalIndicesForLocalQx,
//  					&ltog);
//      ierr = VecSetLocalToGlobalMapping(pQx_,ltog);
//      //mwf check some more localToGlobal mapping stuff
//      ierr = ISLocalToGlobalMappingApply(ltog,local_nQx,indicesForLocalQx,globalIndicesCheck);
//      //cout<<" now printing out ltog"<<endl;
//      //ierr = ISLocalToGlobalMappingView(ltog,PETSC_VIEWER_STDOUT_SELF);
//      ierr = ISLocalToGlobalMappingDestroy(ltog);

//    }//end local scope
  
//    ierr = PetscIntView(local_nQx,globalIndicesCheck,PETSC_VIEWER_STDOUT_SELF);

//    //insert values from local_Qx
//    ierr = VecSetValuesLocal(pQx_,local_nQx,indicesForLocalQx,
//  			   localQxValues,INSERT_VALUES);


//    //clean up
//    delete [] globalIndicesForLocalQx;
//    delete [] localQxValues;
//    delete [] indicesForLocalQx;
//    delete [] globalIndicesCheck;

//    ierr = VecAssemblyBegin(pQx_);
//    ierr = VecAssemblyEnd(pQx_);

//    /*
//        View the vector; then destroy it.
//    */
//    ierr = VecView(pQx_,PETSC_VIEWER_STDOUT_WORLD);
//    ierr = VecDestroy(pQx_);
//  #endif  
//  }

//  //try to do this like chris' << routine and ex8 in vec/examples
//  void Petsc::postProcessFluxesParallel(Petsc::StencilMM& stencil,
//  				      Vec& Qx,Vec& Qy,Vec& Qz)
//  {
//    Err ierr;
//    Petsc::Sys pSys;
 
//    int rank(0),size(0);
//    MPI_Status        status;
//    //get number of processors and local processor number
//    rank = pSys.getRank();
//    size = pSys.getSize();

//    //do Qx first
//    int local_nxNodes = stencil.local_nxNodes;
//    int local_nyNodes = stencil.local_nyNodes;
//    int local_nzNodes = stencil.local_nzNodes;


//  #ifndef TRY_NEW_GLOBAL_VEC  
//    PetscViewer viewer;

//    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"qxTest.grf",&viewer);
//    ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX);
//    ierr = VecView(Qx.castToPetsc(),viewer);
//    ierr = PetscViewerDestroy(viewer);
//  #else  
//    Petsc::cc::_p_Vec*  pQx_;
//    int local_nQx = (local_nxNodes+1)*(local_nyNodes);
//    ierr = VecCreate(PETSC_COMM_WORLD,local_nQx,PETSC_DETERMINE,&pQx_);
//    ierr = VecSetFromOptions(pQx_);

 
//    //create set of global indices corresponding to this processors local
//    //nodes

//    int * globalIndicesForLocalQx = new int[local_nQx];
//    int * globalIndicesCheck= new int[local_nQx];
//    int * indicesForLocalQx = new int[local_nQx];

//    //also get array holding local values
//    real * localQxValues = new real[local_nQx];

//    int viewRank = 1;
//    if (rank == viewRank)
//      cout <<"in postProcessFluxesParallel on proc "<<rank<<endl;

//    //now loop through local indeces and set local flux vector using Qx
//    for (int j=0;j<local_nyNodes;j++)
//      for (int k=0; k<local_nxNodes;k++)
//        {
//  	stencil.localIndex(j,k);
//  	int interRight = stencil.interRight;
//  	localQxValues[interRight] = Qx[interRight];
//  	indicesForLocalQx[interRight] = interRight;
//  	//need to put in global flux index into stencil
//  	int jglb = j + stencil.local_y0;
//  	int kglb = k + stencil.local_x0;
//  	stencil(jglb,kglb);
//  	globalIndicesForLocalQx[interRight] = jglb*(stencil.nxNodes+1) + kglb + 1;

//  	if (rank == viewRank)
//  	  {
//  	    //cout<<"(j,k)= "<<j<<","<<k<<" (jglb,kglb)= "<<jglb<<","<<kglb;
//  	    cout<<endl<<" indicesForLocalQx["<<interRight<<"]= "
//  		<<indicesForLocalQx[interRight];
//  	    cout<<" localQxValues["<<interRight<<"]= "
//  		<<localQxValues[interRight];
//  	    cout<<" globalIndicesForLocalQx["<<interRight<<"]= "
//  		<<globalIndicesForLocalQx[interRight];
//  	    cout<<" petscToGlobal("<<stencil.center<<")= "
//  		<<stencil.petscToGlobal(stencil.center);
//  	  }
	    
//        }
  
//    int k = 0;
//    for (int j=0;j<local_nyNodes;j++)
//      {
//        stencil.localIndex(j,k);
//        int interLeft = stencil.interLeft;
//        if (stencil.local_x0 > 0)
//  	{
//  	  localQxValues[interLeft] = -123456.7;
//  	  indicesForLocalQx[interLeft] = -1;
//  	  globalIndicesForLocalQx[interLeft] = -1;
//  	}
//        else
//  	{
//  	  localQxValues[interLeft] = Qx[interLeft];
//  	  indicesForLocalQx[interLeft] = interLeft;
//  	  //need to put in global flux index into stencil
//  	  int jglb = j + stencil.local_y0;
//  	  int kglb = k + stencil.local_x0;
//  	  stencil(jglb,kglb);
//  	  globalIndicesForLocalQx[interLeft] = jglb*(stencil.nxNodes+1) + kglb;

//  	  if (rank == viewRank)
//  	    {
//  	      //cout<<"(j,k)= "<<j<<","<<k<<" (jglb,kglb)= "<<jglb<<","<<kglb;
//  	      cout<<endl<<" indicesForLocalQx["<<interLeft<<"]= "
//  		  <<indicesForLocalQx[interLeft];
//  	      cout<<" localQxValues["<<interLeft<<"]= "
//  		<<localQxValues[interLeft];
//  	      cout<<" globalIndicesForLocalQx["<<interLeft<<"]= "
//  		<<globalIndicesForLocalQx[interLeft];
//  	      cout<<" petscToGlobal("<<stencil.center<<")= "
//  		  <<stencil.petscToGlobal(stencil.center);
//  	    }

//  	}      
//      }
//    //create mapping from these local indeces to the global indeces

//    {//local scope
//      ISLocalToGlobalMapping ltog;
//      ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_SELF,local_nQx,
//  					globalIndicesForLocalQx,
//  					&ltog);
//      ierr = VecSetLocalToGlobalMapping(pQx_,ltog);
//      //mwf check some more localToGlobal mapping stuff
//      ierr = ISLocalToGlobalMappingApply(ltog,local_nQx,indicesForLocalQx,globalIndicesCheck);
//      //cout<<" now printing out ltog"<<endl;
//      //ierr = ISLocalToGlobalMappingView(ltog,PETSC_VIEWER_STDOUT_SELF);
//      ierr = ISLocalToGlobalMappingDestroy(ltog);

//    }//end local scope
  
//    ierr = PetscIntView(local_nQx,globalIndicesCheck,PETSC_VIEWER_STDOUT_SELF);

//    //insert values from local_Qx
//    ierr = VecSetValuesLocal(pQx_,local_nQx,indicesForLocalQx,
//  			   localQxValues,INSERT_VALUES);


//    //clean up
//    delete [] globalIndicesForLocalQx;
//    delete [] localQxValues;
//    delete [] indicesForLocalQx;
//    delete [] globalIndicesCheck;

//    ierr = VecAssemblyBegin(pQx_);
//    ierr = VecAssemblyEnd(pQx_);

//    /*
//        View the vector; then destroy it.
//    */
//    ierr = VecView(pQx_,PETSC_VIEWER_STDOUT_WORLD);
//    ierr = VecDestroy(pQx_);
//  #endif  
  
//  }

//try to do this like chris' << routine 
void printOutFluxError(real L2errQ, real L2exQ)
{
  Err ierr;
  Petsc::Sys pSys;
  int tag=10;
  int               rank,len,n,size;//work,
  MPI_Status        status;
  real  values[2];
  real  array[2];
  real finalValues[2];

  /* determine maximum message to arrive */
  rank = pSys.getRank();
  size = pSys.getSize();
  //work = 2;
  len  = 2;

  //load in error values into array
  array[0] = L2errQ;
  array[1] = L2exQ;

  if (!rank) 
    {
      finalValues[0] = array[0];
      finalValues[1] = array[1];

      /* receive and print messages */

      for (int j=1; j<size; j++) 
        {
          ierr = MPI_Recv(values,len,MPIU_SCALAR,j,tag,
			  PETSC_COMM_WORLD,&status);
          ierr = MPI_Get_count(&status,MPIU_SCALAR,&n);         
	  assert(n == 2);
          for (int i=0; i<n; i++) 
            finalValues[i] += values[i];
        } 

      //mwf now print out to 
      cout<<" on processor 0, L2errQ= "<<sqrt(finalValues[0])
	  <<" L2exQ  = "<<sqrt(finalValues[1])<<endl;
    } 
  else 
    {
    /* send values */
      ierr = MPI_Send(array,2,MPIU_SCALAR,0,tag,PETSC_COMM_WORLD);
    }
  
}

//try to do this like chris' << routine 
void printOutError(real L2errQ, real L2exQ,
			  const char * nameErr,
			  const char * nameEx)
{
  Err ierr;
  Petsc::Sys pSys;
  int tag=10;
  int               rank,len,n,size;//work,
  MPI_Status        status;
  real  values[2];
  real  array[2];
  real finalValues[2];

  /* determine maximum message to arrive */
  rank = pSys.getRank();
  size = pSys.getSize();
  //work = 2;
  len  = 2;

  //load in error values into array
  array[0] = L2errQ;
  array[1] = L2exQ;

  if (!rank) 
    {
      finalValues[0] = array[0];
      finalValues[1] = array[1];

      /* receive and print messages */

      for (int j=1; j<size; j++) 
        {
          ierr = MPI_Recv(values,len,MPIU_SCALAR,j,tag,
			  PETSC_COMM_WORLD,&status);
          ierr = MPI_Get_count(&status,MPIU_SCALAR,&n);         
	  assert(n == 2);
          for (int i=0; i<n; i++) 
            finalValues[i] += values[i];
        } 

      //mwf now print out to 
      cout<<" on processor 0, "<<nameErr<<" = "<<sqrt(finalValues[0])
	  <<" "<<nameEx<<"  = "<<sqrt(finalValues[1])<<endl;
    } 
  else 
    {
    /* send values */
      ierr = MPI_Send(array,2,MPIU_SCALAR,0,tag,PETSC_COMM_WORLD);
    }
  
}

}//Petsc
}//Daetk
