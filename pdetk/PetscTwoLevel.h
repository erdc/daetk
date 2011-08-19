#ifndef PETSCTWOLEVEL_H
#define PETSCTWOLEVEL_H

#include <fstream>
#include "Definitions.h"
#include "IntVec.h"
#include "BandColMat.h"
#include "LinearSolver.h"
#include "LaFullDirectSolver.h"
#include "LaBandedDirectSolver.h"
#include "DataCollector.h"
#include "Definitions.h" 
#include "WeightedRMSNorm.h"
#include "VectorFunction.h"
#include "Preconditioner.h"
#include "LinearOperator.h"
#include "PetscStencilMat.h"
#include "DaetkPetscVec.h"
#include "DaetkPetscMat.h"
#include "DaetkPetscLinearSolver.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
namespace Petsc
{

template <class STENCIL, int nv>
class TwoLevel : public Daetk::Preconditioner
{ 
public:
  TwoLevel(Petsc::Mat& mIn, 
           LinearOperator& linOpIn,
           int NXIn,
           int NYIn,
           int NZIn, 
           STENCIL& stencilIn, DataCollector& dataIn, bool myASM=false);
  virtual ~TwoLevel();

  bool prepare(); 

  void coarseToFine(const Vec& cVec, Vec& fVec);
  void fineToCoarse(const Vec& fVec, Vec& cVec);
  void fineToSubdomain(const Vec& fVec, Vec* subDVec);
  void subdomainToFine(const Vec* subDVec, Vec& fVec);

  bool apply(const Vec& p, Vec& Mp);

  void setPrec(bool use_coarse_flag,
               bool do_fine_last_flag = false,
               bool use_multiplicative_flag = false);


  void setCoarseGridTolerances(real rtol, real atol, real dtol, int maxit);
  void printMatrices(const char* filename);

  void doNothing();
  void useMyAdditiveSchwarz();
private:
  bool USE_COARSE, DO_FINE_LAST, MULTIPLICATIVE,
    PRINT_MATRICES,DO_NOTHING,USE_MY_ADDITIVE_SCHWARZ,error;

  int nNodes,nx,NX,NX_Local,ny,NY,NY_Local,nz,NZ,NZ_Local,NXY,NXY_Local,nxExtra,nyExtra,nzExtra,nsd,nsdLocal, coarseGlobalLow;
  real wIJ,wIJP1,wIJM1,wIP1J,wIM1J;
  Vec xCoarse,bCoarse,coarseCorrection,fineCorrection;
  Vec  *xSubD,*bSubD;
  Petsc::Mat& fineMat;
  Petsc::Mat *subDomainMat,coarseMat;
  IntVec nxSubD,nySubD,nzSubD,nNodesSubD;
  std::vector< std::vector<int> > subdomainIndeces;
  Petsc::LinearSolver coarseSolver,*subDomainSolver;
  Petsc::LinearSolver subdomainSolver2;
  std::ofstream matOut;

  bool solveCoarsePrec(const Vec& v, Vec& pv);
  bool solveFinePrec(const Vec& v, Vec& pv);

  Norm2 vnorm;
  LinearOperator* linearOperator;
  DataCollector emptyData;
  STENCIL& stencil,coarseStencil;
  struct Location { int i,j,k,n;};
  CMRVec<Location> coarseLocation,subdomainLocation;
  void calculateWeights(Location& subD,Location& coarse);

  inline int computeOffPCoarseNode(int i,int j,int k);

  inline void computeNodeIndeces(int i,int j,int k, Location& coarse, Location& subdomain);
  
  template <int nvf>
  friend void finalize(int i,Petsc::Mat& m);
  template <int nvf>
  friend void finalizeAdd(int i,Petsc::Mat& m);
};


template <int nvf>
void finalize(int i,Petsc::Mat& m)
{
  m.finalizeBlockRow(i);
}

template <int nvf>
void finalizeAdd(int i,Petsc::Mat& m)
{
  m.finalizeAddBlockRow(i); 
}

template <>
void finalize<1>(int i,Petsc::Mat& m)
{
  m.finalizeRow(i);
}

template <>
void finalizeAdd<1>(int i,Petsc::Mat& m)
{
  m.finalizeAddRow(i); 
}

template <class STENCIL, int nv>
int TwoLevel<STENCIL,nv>::computeOffPCoarseNode(int i,int j,int k)
  {
    int I,J,K;
    i-=stencil.local_z0;
    j-=stencil.local_y0;
    k-=stencil.local_x0;
    
    if (i < 0)
      I = coarseStencil.local_z0-1;
    else if (i >= stencil.local_nzNodes)
      I = coarseStencil.local_z0 + coarseStencil.local_nzNodes;
    else
      I =  std::min(i/nzSubD[0],NZ_Local-1) + coarseStencil.local_z0;

    if (j < 0)
      J = coarseStencil.local_y0-1;
    else if (j >= stencil.local_nyNodes)
      J = coarseStencil.local_y0 + coarseStencil.local_nyNodes;
    else
      J =  std::min(j/nySubD[0],NY_Local-1)+ coarseStencil.local_y0;

    if (k < 0)
      K = coarseStencil.local_x0-1;
    else if  (k >= stencil.local_nxNodes)
      K = coarseStencil.local_x0 + coarseStencil.local_nxNodes;
    else
      K =  std::min(k/nxSubD[0],NX_Local-1)+ coarseStencil.local_x0;

    coarseStencil(I,J,K);
    return coarseStencil.center;
  }

template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::computeNodeIndeces(int i,int j,int k, Location& coarse, Location& subdomain)
    {
      i-=stencil.local_z0;
      j-=stencil.local_y0;
      k-=stencil.local_x0;
      coarse.i = std::min(i/nzSubD[0],NZ_Local-1);
      coarse.j = std::min(j/nySubD[0],NY_Local-1); //catches edge subdomains with slop
      coarse.k = std::min(k/nxSubD[0],NX_Local-1); // ditto
      coarse.n = coarse.i*NXY_Local + coarse.j*NX_Local + coarse.k;
      subdomain.i = i%nzSubD[0] + (std::max(i/nzSubD[0],NZ_Local-1) - NZ_Local+1)*nzSubD[0];//catches edge's again
      subdomain.j = j%nySubD[0] + (std::max(j/nySubD[0],NY_Local-1) - NY_Local+1)*nySubD[0];
      subdomain.k = k%nxSubD[0] + (std::max(k/nxSubD[0],NX_Local-1) - NX_Local+1)*nxSubD[0];
      subdomain.n = subdomain.i*nySubD[coarse.j]*nxSubD[coarse.k]+subdomain.j*nxSubD[coarse.k]+subdomain.k;
    }


template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::coarseToFine(const Vec& cVec, Vec& fVec)
{
  const int end=stencil.local_nxyzNodes;
  for (int n=0;n<end;n++)
    {
      const Location &nodeCoarse(coarseLocation[n]);
      for (int vj=0;vj<nv;vj++)
        fVec[nv*n+vj] = cVec[nv*nodeCoarse.n+vj];
    }
}

template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::fineToCoarse(const Vec& fVec, Vec& cVec)
{
  const int end=stencil.local_nxyzNodes;
  for (int n=0;n<end;n++)
    {
      const Location &nodeCoarse(coarseLocation[n]);
      for (int vj=0;vj<nv;vj++)
        cVec[nv*nodeCoarse.n+vj] += fVec[nv*n+vj];
    }
}

template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::fineToSubdomain(const Vec& fVec, Vec* subDVec)
{
  const int end=stencil.local_nxyzNodes;
  for (int n=0;n<end;n++)
    {
      const Location &nodeCoarse(coarseLocation[n]), &nodeSubdomain(subdomainLocation[n]);
      for (int vj=0;vj<nv;vj++)
        subDVec[nodeCoarse.n][nv*nodeSubdomain.n+vj] = fVec[nv*n+vj];
    }
}

template  <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::subdomainToFine(const Vec* subDVec, Vec& fVec)
{
  const int end=stencil.local_nxyzNodes;
  for (int n=0;n<end;n++)
    {
      const Location &nodeCoarse(coarseLocation[n]), &nodeSubdomain(subdomainLocation[n]);
      for (int vj=0;vj<nv;vj++)
        fVec[nv*n+vj] = subDVec[nodeCoarse.n][nv*nodeSubdomain.n+vj];
    }  
}

template<class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::setCoarseGridTolerances(real rtol, real atol, real dtol, int maxits)
{
  coarseSolver.setTolerances(rtol,atol,dtol,maxits);
}
  
template<class STENCIL, int nv>
TwoLevel<STENCIL,nv>::TwoLevel(Petsc::Mat& mIn, 
                                      LinearOperator& linOpIn,
                                      int NXIn,
                                      int NYIn,
                                      int NZIn, 
                                      STENCIL& stencilIn, 
                                      DataCollector& dataIn,
                                      bool myASM):
  USE_COARSE(true),
  DO_FINE_LAST(false),
  MULTIPLICATIVE(true),
  PRINT_MATRICES(false),
  DO_NOTHING(false),
  USE_MY_ADDITIVE_SCHWARZ(myASM),
  error(true),
  nNodes((stencilIn.nxNodes)*(stencilIn.nyNodes)*(stencilIn.nzNodes)),
  nx(stencilIn.nxNodes),
  NX(NXIn),
  ny(stencilIn.nyNodes),
  NY(NYIn),
  nz(stencilIn.nzNodes),
  NZ(NZIn),
  NXY(NXIn*NYIn),
  nsd(NXIn*NYIn*NZIn),
  coarseCorrection(nv*nNodes,0.0),
  fineCorrection(nv*nNodes,0.0),
  xSubD(0),
  bSubD(0),
  fineMat(mIn),
  subDomainMat(0),
  coarseSolver(emptyData),
  matOut(),
  linearOperator(&linOpIn),
  stencil(stencilIn),
  coarseStencil(stencilIn,NXIn,NYIn,NZIn,nv),
  coarseLocation(stencilIn.local_nxyzNodes),
  subdomainLocation(stencilIn.local_nxyzNodes)
{
  //make sure user hasn't specified an odd decomposition
  if( (NX % stencil.npx) != 0)
    {
      std::cerr<<"Subdomains must be evenly divided among processors"<<std::endl
          <<"Reseting input value NX ="<<NX<<" to ";
      NX = (NX/stencil.npx) * stencil.npx;
      std::cerr<<NX<<" for npx = "<<stencil.npx<<" processors"<<std::endl;
      exit(1);
    }

  if( (NY % stencil.npy) != 0)
    {
      std::cerr<<"Subdomains must be evenly divided among processors"<<std::endl
          <<"Reseting input value NY ="<<NY<<" to ";
      NY = (NY/stencil.npy) * stencil.npy;
      std::cerr<<NY<<" for npy = "<<stencil.npy<<" processors"<<std::endl;
      exit(1);
    }

  if( (NZ % stencil.npz) != 0)
    {
      std::cerr<<"Subdomains must be evenly divided among processors"<<std::endl
          <<"Reseting input value NZ ="<<NZ<<" to ";
      NZ = (NZ/stencil.npz) * stencil.npz;
      std::cerr<<NZ<<" for npz = "<<stencil.npz<<" processors"<<std::endl;
      exit(1);
    }

  NX_Local = NX/stencil.npx;
  NY_Local = NY/stencil.npy;
  NZ_Local = NZ/stencil.npz;
  NXY_Local = NX_Local*NY_Local;
  
  nxSubD.newsize(NX_Local); nxSubD = stencil.local_nxNodes/NX_Local;
  nySubD.newsize(NY_Local); nySubD = stencil.local_nyNodes/NY_Local;
  nzSubD.newsize(NZ_Local); nzSubD = stencil.local_nzNodes/NZ_Local;

  nsdLocal = NX_Local*NY_Local*NZ_Local;
  subdomainIndeces.resize(nsdLocal);
  nNodesSubD.newsize(nsdLocal);

  xCoarse.newsize(Vec::GLOBAL,coarseStencil.dadof_);
  xCoarse.setExample();
  assert(xCoarse.dim() == nsd*nv);
  bCoarse.newsize(nsd*nv);
  coarseGlobalLow = coarseStencil.globalLow;

  if (USE_MY_ADDITIVE_SCHWARZ)
    {
      subDomainMat = new Petsc::Mat[nsdLocal];
      subDomainSolver = new Petsc::LinearSolver[nsdLocal];

      xSubD = new Vec[nsdLocal];
      bSubD = new Vec[nsdLocal];
    }
  //the last subdomain in each direction will have the extra nodes

  nxExtra = stencil.local_nxNodes%NX_Local;
  nyExtra = stencil.local_nyNodes%NY_Local;
  nzExtra = stencil.local_nzNodes%NZ_Local;

  nxSubD[NX_Local - 1] += nxExtra;
  nySubD[NY_Local - 1] += nyExtra;
  nzSubD[NZ_Local - 1] += nzExtra;

  bool oneD=false,twoD=false;
  if (ny == 1)
    oneD=true;
  else if (nz == 1)
    twoD=true;

  int offDiagonals;
  if (oneD)
    offDiagonals=0;
  else if (twoD)
    offDiagonals = 2;
  else
    offDiagonals = 4;
  
  //need to create with local coarse dim & possibly also another with local storage + coarse solvers
  bool isSymmetric=false; //will work for symmetric too.

  coarseMat.newsize(nv*nsd,3,offDiagonals,isSymmetric,nv);
  coarseSolver.attachMat(coarseMat);
  coarseSolver.attachPreconditioner(coarseMat);
  coarseSolver.attachData(emptyData);
  coarseSolver.setMethod(Petsc::LinearSolver::PREONLY);
  if (stencil.npx*stencil.npy*stencil.npz == 1)
    {
      coarseSolver.setPreconditioner(Petsc::LinearSolver::LU);
    }
  else
    {
      coarseSolver.setPreconditioner(Petsc::LinearSolver::REDUNDANT);
      coarseSolver.setRedundantPreconditioner(Petsc::LinearSolver::LU);
    }
  for (int N=0;N<nsdLocal;N++)
    {
      int I=N/NXY_Local,
        J=(N/NX_Local)%NY_Local,
        K=N%NX_Local;
      
      nNodesSubD[N] = nxSubD[K] * nySubD[J]* nzSubD[I];

      if (USE_MY_ADDITIVE_SCHWARZ)
        {
          subDomainMat[N].newsizeSequential(nv*nNodesSubD[N],3,offDiagonals,isSymmetric,nv);
          subDomainSolver[N].attachSerialMat(subDomainMat[N]);
          subDomainSolver[N].attachPreconditioner(subDomainMat[N]);
          subDomainSolver[N].attachData(emptyData);
          subDomainSolver[N].setMethod(Petsc::LinearSolver::PREONLY);
          subDomainSolver[N].setPreconditioner(Petsc::LinearSolver::LU);
          
//            subDomainSolver[N].setMethod(Petsc::LinearSolver::GMRES);
//            subDomainSolver[N].setPreconditioner(Petsc::LinearSolver::NONE);
          xSubD[N].newsizeSerial(nv*nNodesSubD[N]);
          bSubD[N].newsizeSerial(nv*nNodesSubD[N]);
        }
    }

  for (int n=stencil.globalLow;n<stencil.globalHigh;n++)
    {
      stencil.globalNode(stencil.petscToGlobal(n));
      assert(n==stencil.center);
      computeNodeIndeces(stencil.anchor->i,
                         stencil.anchor->j,
                         stencil.anchor->k,
                         coarseLocation[n - stencil.globalLow],
                         subdomainLocation[n - stencil.globalLow]);
      subdomainIndeces[ coarseLocation[n-stencil.globalLow].n ].push_back(n);
    }

  if (!USE_MY_ADDITIVE_SCHWARZ)
    {
      bool clear=true;
      subdomainSolver2.attachMat(fineMat,clear);
      subdomainSolver2.attachPreconditioner(fineMat);
      subdomainSolver2.attachData(emptyData);
      subdomainSolver2.setMethod(Petsc::LinearSolver::PREONLY);
      subdomainSolver2.setPreconditioner(Petsc::LinearSolver::ASM);
      subdomainSolver2.additiveSchwarz.setType(Petsc::LinearSolver::AdditiveSchwarz::BASIC);
      subdomainSolver2.additiveSchwarz.setSubdomainIndexSets(subdomainIndeces,nv);
      subdomainSolver2.additiveSchwarz.setOverlap(0);
    }
  Tracer tr("Petsc::TwoLevel<STENCIL>::Petsc::TwoLevel<STENCIL>(...)");
}

template <class STENCIL, int nv>
TwoLevel<STENCIL,nv>::~TwoLevel()
{
  Tracer tr("TwoLevel<STENCIL>::~TwoLevel<STENCIL>()");
  if (USE_MY_ADDITIVE_SCHWARZ)
    {
      delete [] subDomainSolver;
      delete [] subDomainMat;
      delete [] xSubD;
      delete [] bSubD;
    }
}

template <class STENCIL, int nv>
bool TwoLevel<STENCIL,nv>::prepare()
{ 
  Petsc::Sys psys;
  coarseMat.zeroAll();

  if (USE_MY_ADDITIVE_SCHWARZ)
    {
      for (int sd=0;sd<nsdLocal;sd++)
        {
          subDomainMat[sd].zeroAll();
        }
    }

  if(PRINT_MATRICES)
    fineMat.print();
  
  for (int n=stencil.globalLow;n<stencil.globalHigh;n++)
    { 
      int CNbar;
      stencil.setGlobalPointList(stencil.petscToGlobal(n));
      int nLocal = n-stencil.globalLow;

      const Location &nodeCoarse(coarseLocation[nLocal]),&nodeSubdomain(subdomainLocation[nLocal]);
      CNbar = nodeCoarse.n + coarseGlobalLow;
      typename STENCIL::iterator sit = stencil.begin();
      const typename STENCIL::iterator end=stencil.end();

      while (sit != end)
        {
          int CN,sn;

          //column entry is local
          if ((sit->globalNodeNumber >= stencil.globalLow) && (sit->globalNodeNumber < stencil.globalHigh))
            {
              int nLocalIt = sit->globalNodeNumber - stencil.globalLow;
              const Location &pointCoarse(coarseLocation[nLocalIt]),&pointSubdomain(subdomainLocation[nLocalIt]);
              CN = pointCoarse.n+coarseGlobalLow;
              sn = pointSubdomain.n;
            }
          else //column entry is off proc
            {
              CN = computeOffPCoarseNode(sit->i,sit->j,sit->k);
              if (CN==CNbar)
		{
		  std::cout<<"node dependency "<<sit->globalNodeNumber<<std::endl
		      <<"Is being located in coarse cell"<<CN<<std::endl
		      <<"which is in this nodes coasrse cell"<<CNbar<<std::endl
		      <<"This node is "<<n<<std::endl;
		  exit(1);
		}
              sn=0;
            }

          if (USE_MY_ADDITIVE_SCHWARZ)
            {
              if (CN == CNbar)
                {
                  for (int vk=0;vk<nv;vk++)
                    for (int vj=0;vj<nv;vj++)
                      {
                        subDomainMat[nodeCoarse.n](nv*nodeSubdomain.n+vj,nv*sn+vk) 
                          = fineMat.get(n*nv+vj,sit->globalNodeNumber*nv+vk);
                      }
                  finalize<nv>(nodeSubdomain.n,subDomainMat[nodeCoarse.n]);
//                    subDomainMat[nodeCoarse.n].finalizeBlockRow(nodeSubdomain.n);
                }
            }

          if (USE_COARSE )   
            {
              for (int vk=0;vk<nv;vk++)
                for (int vj=0;vj<nv;vj++)
                  coarseMat(nv*CNbar + vj,nv*CN+vk) 
                    += fineMat.get(n*nv+vj,sit->globalNodeNumber*nv+vk);
              finalizeAdd<nv>(CNbar,coarseMat);
//                coarseMat.finalizeAddBlockRow(CNbar);
            }
          ++sit;
        }
    }
  if (USE_COARSE)
    coarseMat.beginAssembly();
  
  if (USE_MY_ADDITIVE_SCHWARZ)
    {
      for (int sd=0;sd<nsdLocal;sd++)
        {
          subDomainMat[sd].beginAssembly();  
          subDomainMat[sd].endAssembly();
        }
    

      for (int sd=0;sd<nsdLocal;sd++)
        {
          if(PRINT_MATRICES)
            subDomainMat[sd].print();
          error = subDomainSolver[sd].prepare();
          if (error) 
            {
              std::cerr<<"Subdomain Matrix is singular"<<std::endl;
              return(1);
            }
        }
    }
  else
    {
      subdomainSolver2.prepare();
      subdomainSolver2.additiveSchwarz.setSubdomainSolvers();
      //sets to PREONLY and LU prec.
    }
  if (USE_COARSE)
    {
      coarseMat.endAssembly();
      if(PRINT_MATRICES)
        coarseMat.print();      
      error=coarseSolver.prepare();
      if (error) 
        {
          std::cerr<<"Coarse Matrix is singular"<<std::endl;
          return(1);
        }
    }
  return 0;
}



template <class STENCIL, int nv>
bool TwoLevel<STENCIL,nv>::solveCoarsePrec(const Vec& v, Vec& pv)
{  
  pv=0.0;
  bCoarse=0.0;

  fineToCoarse(v,bCoarse);

  error=coarseSolver.apply(bCoarse,xCoarse);
  if (error)
    {
      std::cout<<"coarse solver failure"<<std::endl;
      pv=0.0;
    } 
  else
    {
      coarseToFine(xCoarse,pv);
    }
  
  return 0;
}



template <class STENCIL, int nv>
bool TwoLevel<STENCIL,nv>::solveFinePrec(const Vec& v, Vec& pv)
{
  fineToSubdomain(v,bSubD);

  for (int sd=0;sd<nsdLocal;sd++)
    {
      error=subDomainSolver[sd].apply(bSubD[sd],xSubD[sd]);
      if (error)
        {
          std::cerr<<"Fine grid solve failed"<<std::endl;
        }
    }
 
  subdomainToFine(xSubD,pv);

  return 0;
}

template <class STENCIL, int nv>
bool TwoLevel<STENCIL,nv>::apply(const Vec& v, Vec& pv)
{ 
  //  Vec fineCorrection2(fineCorrection),pv2(pv);

  if (DO_NOTHING)
    {
#ifndef USE_BLAS
    pv=v;
#else
    copy(v,pv);
#endif
    }
  else
    {
#ifndef USE_BLAS
      coarseCorrection = v;
      fineCorrection = v;
      pv = v;
#else
      copy(v,coarseCorrection);
      copy(v,fineCorrection);
      copy(v,pv);
#endif
      if ((!DO_FINE_LAST))
        {
          if (USE_MY_ADDITIVE_SCHWARZ)
            error=solveFinePrec(v,fineCorrection);
          else
            error = subdomainSolver2.apply(v,fineCorrection);
          if (error)
            {
              std::cerr<<"Error solving fine prec"<<std::endl;
              return 1;
            }
          if (USE_COARSE && MULTIPLICATIVE)
            {
              bool evalError=false;
              evalError=linearOperator->apply(fineCorrection,coarseCorrection);
              while(evalError)
                {
                  std::cerr<<"cutting back on fineCorrection in TwoLevel"<<std::endl;
                  fineCorrection*=0.5;
                  evalError=linearOperator->apply(fineCorrection,coarseCorrection);
                }
#ifndef USE_BLAS
              pv -= coarseCorrection;
#else
              axpy(-1.0,coarseCorrection,pv);
#endif
            }
        }
      
      if (USE_COARSE)
        {
          error=solveCoarsePrec(pv,coarseCorrection);
          if (error)
            {
              std::cerr<<"Error solving coarse matrix"<<std::endl;
              return 1;
            }
        }
      
      if (DO_FINE_LAST)
        {
          if (USE_COARSE && MULTIPLICATIVE)
            {
              bool evalError=false;
              evalError=linearOperator->apply(coarseCorrection,fineCorrection);
              while (evalError)
                {
                  std::cerr<<"cutting back coarse Correction in TwoLevel"<<std::endl;
                  coarseCorrection*=0.5;
                  evalError=linearOperator->apply(coarseCorrection,fineCorrection);                  
                }
#ifndef USE_BLAS
              pv = v - fineCorrection;
#else
              copy(v,pv);
              axpy(-1.0,fineCorrection,pv);
#endif
            }
          if (USE_MY_ADDITIVE_SCHWARZ)
            error=solveFinePrec(pv,fineCorrection);
          else
            error = subdomainSolver2.apply(pv,fineCorrection);

          if (error)
            {
              std::cerr<<"Error solving fine prec"<<std::endl;
              return 1;
            }
        }
      if (USE_COARSE)
        {
#ifndef USE_BLAS
          pv = fineCorrection+coarseCorrection;
#else
          copy(fineCorrection,pv);
          axpy(1.0,coarseCorrection,pv);
#endif
        }
      else
        {
#ifndef USE_BLAS
        pv = fineCorrection;
#else
        copy(fineCorrection,pv);
#endif
        }
    }
  return 0;
}

template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::setPrec(bool use_coarse_flag,
                                   bool do_fine_last_flag,
                                   bool use_multiplicative_flag)
{
  USE_COARSE = use_coarse_flag;
  DO_FINE_LAST = do_fine_last_flag;
  MULTIPLICATIVE = use_multiplicative_flag;
} 

template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::printMatrices(const char*filename)
{PRINT_MATRICES = true; matOut.open(filename);}


template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::calculateWeights(Location& subD,Location& coarse)
{
  if (subD.k < nxSubD(0)/2)
    {
      if (coarse.k > 0)
        {
          wIJ = (real(subD.k) + 0.5)/real(nxSubD(0)) + 0.5;
          wIJM1 = 1.0 - wIJ;
          wIJP1 = 0.0;
        }
      else
        {
          wIJ = real(subD.k + 0.5)/real(nxSubD(0));
          wIJM1 = 0.0;
          wIJP1 = 0.0;
        }
    }
  else
    {
      if (coarse.k < NX-1)
        {
          wIJ = 1.5 - (real(subD.k) + 0.5)/real(nxSubD(0));
          wIJM1 = 1.0 - wIJ;
          wIJP1 = 0.0;
        }
      else
        {
          wIJ = 1.0 - (real(subD.k) + 0.5)/real(nxSubD(0));
          wIJM1 = 0.0;
          wIJP1 = 0.0;
        }
    }
  if (subD.j < nySubD(0)/2)
    {
      if (coarse.j > 0)
        {
          wIJ *= (real(subD.j) + 0.5)/real(nySubD(0)) + 0.5;
          wIM1J = 0.5 - (real(subD.j) + 0.5)/real(nySubD(0));
          wIP1J = 0.0;
        }
      else
        {
          wIJ *= real(subD.j + 0.5)/real(nySubD(0));
          wIM1J = 0.0;
          wIP1J = 0.0;
        }
    }
  else
    {
      if (coarse.j < NY-1)
        {
          wIJ = 1.5 - (real(subD.j) + 0.5)/real(nySubD(0));
          wIM1J = (real(subD.j) + 0.5)/real(nySubD(0)) - 0.5;
          wIP1J = 0.0;
        }
      else
        {
          wIJ = 1.0 - (real(subD.j) + 0.5)/real(nxSubD(0));
          wIM1J = 0.0;
          wIP1J = 0.0;
        }
    }
}
 
template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::doNothing(){DO_NOTHING=true;}
  
template <class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::useMyAdditiveSchwarz(){USE_MY_ADDITIVE_SCHWARZ=true;}

}//Petsc
}//Daetk
#endif
