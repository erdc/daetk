#ifndef TWOLEVEL_H
#define TWOLEVEL_H

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
#include "StencilMat.h"
#include "Vec.h"
#include "VecOperators.h"
#include "VecBlas.h"

namespace Daetk 
{
template <class STENCIL, int nv>
class TwoLevel : public Preconditioner
{ 
public:
  TwoLevel(StencilMat<STENCIL,nv>& Min, 
           LinearOperator& linOpIn,
           int nxIn,int NXIn,
           int nyIn=1,int NYIn=1, 
           int nzIn=1,int NZIn=1);
  virtual ~TwoLevel();
  bool prepare(); 
  bool apply(const Vec& p, Vec& Mp);
  void setPrec(bool use_coarse_flag,
               bool do_fine_last_flag = false,
               bool use_multiplicative_flag = false);
  void printMatrices(const char* filename);
  void doNothing();
private:
  bool USE_COARSE, DO_FINE_LAST, MULTIPLICATIVE,
    PRINT_MATRICES,DO_NOTHING,error;
  int nNodes,nx,NX,ny,NY,nz,NZ,NXY,nxExtra,nyExtra,nzExtra,nsd;
  real wIJ,wIJP1,wIJM1,wIP1J,wIM1J;
  Vec xCoarse,bCoarse,coarseCorrection,fineCorrection;
  Vec *xSubD,*bSubD;
  StencilMat<STENCIL,nv>& stencilMatrix;
  BandColMat *subDomain,coarseMat;
  IntVec *subDomainPivots,coarsePivots,nxSubD,nySubD,nzSubD,nNodesSubD;
  LaFullDirectSolver fullSolver;
  LaBandedDirectSolver bandSolver;
  std::ofstream matOut;
  bool solveCoarsePrec(const Vec& v, Vec& pv);
  bool solveFinePrec(const Vec& v, Vec& pv);
  LinearOperator* linearOperator;
  STENCIL& stencil;
  struct Location { int i,j,k,n;};
  CMRVec<Location> coarseLocation,subdomainLocation;
  void calculateWeights(Location& subD,Location& coarse);
  inline void computeNodeIndeces(int i,int j,int k, Location& coarse, Location& subdomain)
    {
      coarse.i = min(i/nzSubD(0),NZ-1);
      coarse.j = min(j/nySubD(0),NY-1); //catches edge subdomains with slop
      coarse.k = min(k/nxSubD(0),NX-1); // ditto
      coarse.n = coarse.i*NXY + coarse.j*NX + coarse.k;
      subdomain.i = i%nzSubD(0) + (max(i/nzSubD(0),NZ-1) - NZ+1)*nzSubD(0);
      subdomain.j = j%nySubD(0) + (max(j/nySubD(0),NY-1) - NY+1)*nySubD(0);
      subdomain.k = k%nxSubD(0) + (max(k/nxSubD(0),NX-1) - NX+1)*nxSubD(0);
      subdomain.n = subdomain.i*nySubD(coarse.j)*nxSubD(coarse.k)+subdomain.j*nxSubD(coarse.k)+subdomain.k;
    }
};
  
template<class STENCIL, int nv>
void TwoLevel<STENCIL,nv>::doNothing(){DO_NOTHING=true;}

template<class STENCIL, int nv>
TwoLevel<STENCIL,nv>::TwoLevel(StencilMat<STENCIL,nv>& mIn, 
                               LinearOperator& linOpIn,
                               int nxIn,int NXIn,
                               int nyIn,int NYIn,
                               int nzIn,int NZIn):
  USE_COARSE(true),
  DO_FINE_LAST(false),
  MULTIPLICATIVE(true),
  PRINT_MATRICES(false),
  DO_NOTHING(false),
  error(true),
  nNodes(nxIn*nyIn*nzIn),
  nx(nxIn),
  NX(NXIn),
  ny(nyIn),
  NY(NYIn),
  nz(nzIn),
  NZ(NZIn),
  NXY(NXIn*NYIn),
  nxExtra(nxIn%NXIn),
  nyExtra(nyIn%NYIn),
  nzExtra(nzIn%NZIn),
  nsd(NXIn*NYIn*NZIn),
  xCoarse(nv*nsd,0.0),
  bCoarse(nv*nsd,0.0),
  coarseCorrection(nv*nNodes,0.0),
  fineCorrection(nv*nNodes,0.0),
  xSubD(0),
  bSubD(0),
  stencilMatrix(mIn),
  subDomain(0),
  subDomainPivots(0),
  coarsePivots(nv*nsd,0.0),
  nxSubD(NXIn,nxIn/NXIn),
  nySubD(NYIn,nyIn/NYIn),
  nzSubD(NZIn,nzIn/NZIn),
  nNodesSubD(NXIn*NYIn*NZIn,(nxIn/NXIn)*(nyIn/NYIn)*(nzIn/NZIn)),
  fullSolver(),
  bandSolver(),
  matOut(),
  linearOperator(&linOpIn),
  stencil(stencilMatrix.getStencil()),
  coarseLocation(nNodes),
  subdomainLocation(nNodes)
{  
  subDomain = new BandColMat[nsd];
  subDomainPivots = new IntVec[nsd];
  xSubD = new Vec[nsd];
  bSubD = new Vec[nsd];
  //the last subdomain in each direction will have the extra nodes
  nxSubD(NX - 1) += nxExtra;
  nySubD(NY - 1) += nyExtra;
  nzSubD(NZ - 1) += nzExtra;

  int coarseUpperBandWidth, subDomainUpperBandWidth; //= lowerBandWidth
  bool oneD=false,twoD=false;
  if (ny == 1)
    oneD=true;
  else if (nz == 1)
    twoD=true;

  if (oneD)
    coarseUpperBandWidth = nv*1 + nv-1;
  else if (twoD)
    coarseUpperBandWidth = nv*NX + nv - 1;
  else
    coarseUpperBandWidth = nv*NXY + nv - 1;
  
  coarseMat.newsize(coarseUpperBandWidth,coarseUpperBandWidth,nv*nsd);
  
  for (int N=0;N<nsd;N++)
    {
      int I=N/NXY,
        J=(N/NX)%NY,
        K=N%NX;
      
      nNodesSubD(N) = nxSubD(K) * nySubD(J)* nzSubD(I);
      if (oneD)
        subDomainUpperBandWidth = nv*1 + nv-1;
      else if (twoD)
        subDomainUpperBandWidth = nv*nxSubD(K) + nv - 1;
      else
        subDomainUpperBandWidth = nv*nxSubD(K)*nySubD(J) + nv - 1;
      subDomain[N].newsize(subDomainUpperBandWidth,subDomainUpperBandWidth,nv*nNodesSubD(N));
      subDomain[N] = 0.0;
      subDomainPivots[N].newsize(nv*nNodesSubD(N));
      subDomainPivots[N] = 0;
      xSubD[N].newsize(nv*nNodesSubD(N));
      bSubD[N].newsize(nv*nNodesSubD(N));
    }
  for (int n=0;n<nNodes;n++)
    {
      stencil.globalNode(n);
      computeNodeIndeces(stencil.anchor->i,
                         stencil.anchor->j,
                         stencil.anchor->k,
                         coarseLocation(n),
                         subdomainLocation(n));
    }
  Tracer tr("TwoLevel<STENCIL>::TwoLevel<STENCIL>(...)");
}

template <class STENCIL, int nv>
TwoLevel<STENCIL,nv>::~TwoLevel()
{
  Tracer tr("TwoLevel<STENCIL>::~TwoLevel<STENCIL>()");
  delete [] subDomain;
  delete [] subDomainPivots;
  delete [] xSubD;
  delete [] bSubD;
}

template <class STENCIL, int nv>
bool TwoLevel<STENCIL,nv>::prepare()
{ 
  coarseMat = 0.0;
  for (int sd=0;sd<nsd;sd++)
    {
      subDomain[sd] = 0.0;
    }

  if(PRINT_MATRICES)
    matOut<<"Jacobian"<<endl<<stencilMatrix<<endl;

  
  //opt Location nodeCoarse,nodeSubdomain,pointCoarse,pointSubdomain;
  //real wRowIJ,wRowIJP1,wRowIJM1,wRowIP1J,wRowIM1J;
  for (int n=0;n<nNodes;n++)
    {          
      stencil.globalNode(n);
      /* opt computeNodeIndeces(stencil.anchor->i,
                         stencil.anchor->j,
                         stencil.anchor->k,
                         nodeCoarse,nodeSubdomain);*/
      const Location &nodeCoarse(coarseLocation(n)),&nodeSubdomain(subdomainLocation(n));
      //calculateWeights(nodeSubdomain,nodeCoarse);
      //wRowIJ=wIJ;wRowIJP1=wIJP1;wRowIJM1=wIJM1;wRowIP1J=wIP1J;wRowIM1J=wIM1J;
      typename STENCIL::iterator sit = stencil.begin();
      const typename STENCIL::iterator end=stencil.end();
      while (sit != end)
        {
          //opt computeNodeIndeces(sit->i,sit->j,sit->k,pointCoarse,pointSubdomain);
          const Location &pointCoarse(coarseLocation(sit.second->globalNodeNumber)),&pointSubdomain(subdomainLocation(sit.second->globalNodeNumber));
          //          calculateWeights(pointSubdomain,pointCoarse);
          //if the nonzero element is in the same subdomain insert it into subdomain matrix
          if (pointCoarse.n == nodeCoarse.n)
            for (int vj=0;vj<nv;vj++)
              for (int vk=0;vk<nv;vk++)
               {
                 subDomain[nodeCoarse.n](nv*nodeSubdomain.n+vj,nv*pointSubdomain.n+vk) 
                   = stencilMatrix(n,(typename STENCIL::Point)(sit.second->point),vj,vk);
               }
          //add the nonyero element to the coarse mat
          if (USE_COARSE )   
            for (int vj=0;vj<nv;vj++)
              for (int vk=0;vk<nv;vk++)
                coarseMat(nv*nodeCoarse.n + vj,nv*pointCoarse.n+vk) 
                  += stencilMatrix(n,(typename STENCIL::Point)(sit.second->point),vj,vk);
//              for (int vj=0;vj<nv;vj++)
//                for (int vk=0;vk<nv;vk++)
//                  {
//                    coarseMat(nv*nodeCoarse.n + vj,nv*pointCoarse.n+vk) 
//                      += (wIJ/wRowIJ)*stencilMatrix(n,(typename STENCIL::Point)(sit->point),vj,vk);
//                    if (pointCoarse.j < NX-1)
//                      coarseMat(nv*nodeCoarse.n + vj,nv*pointCoarse.n+1+vk) 
//                        += (wIJP1/wRowIJ)*stencilMatrix(n,(typename STENCIL::Point)(sit->point),vj,vk);
//                    if (pointCoarse.j > 0)
//                      coarseMat(nv*nodeCoarse.n + vj,nv*pointCoarse.n-1+vk) 
//                        += (wIJM1/wRowIJ)*stencilMatrix(n,(typename STENCIL::Point)(sit->point),vj,vk);
//                    if (pointCoarse.i < NY-1)
//                      coarseMat(nv*nodeCoarse.n + vj,nv*pointCoarse.n+NX+vk) 
//                        += (wIP1J/wRowIJ)*stencilMatrix(n,(typename STENCIL::Point)(sit->point),vj,vk);
//                    if (pointCoarse.i > 0)
//                      coarseMat(nv*nodeCoarse.n + vj,nv*pointCoarse.n-NX+vk) 
//                        += (wIM1J/wRowIJ)*stencilMatrix(n,(typename STENCIL::Point)(sit->point),vj,vk);
//                  }
          ++sit;
        }
    }
  
  for (int sd=0;sd<nsd;sd++)
    {
      if(PRINT_MATRICES)
        matOut<<"subDomain["<<sd<<"]"<<endl<<subDomain[sd];
      bandSolver.attachMat(subDomain[sd]);
      error = bandSolver.prepare();
      if (error) 
        {
          //cout<<"sten"<<stencilMatrix<<endl;
          //cout<<"subD"<<subDomain[sd]<<endl;
          cerr<<"Subdomain Matrix is singular"<<endl;
          return(1);
        }
    }

  if(PRINT_MATRICES)
    matOut<<"coarseMat"<<endl<<coarseMat<<endl;

  if (USE_COARSE)
    {
      coarsePivots = 0;
      bandSolver.attachMat(coarseMat);
      error=bandSolver.prepare();
      if (error) 
        {
          cerr<<"Coarse Matrix is singular"<<endl;
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
  //opt Location nodeCoarse,nodeSubdomain;
  for (int n=0;n<nNodes;n++)
    {
      stencil.globalNode(n);
      /*computeNodeIndeces(stencil.anchor->i,
                         stencil.anchor->j,
                         stencil.anchor->k,
                         nodeCoarse,nodeSubdomain);*/
      const Location &nodeCoarse(coarseLocation(n));
      for (int vj=0;vj<nv;vj++)
        bCoarse(nv*nodeCoarse.n+vj) += v(nv*n+vj);
//        calculateWeights(nodeSubdomain,nodeCoarse);
//        for (int vj=0;vj<nv;vj++)
//          {
//            bCoarse(nv*nodeCoarse.n+vj)      += v(nv*n+vj)/wIJ;
//            if (nodeCoarse.k < NX-1)
//              bCoarse(nv*nodeCoarse.n+1+vj)  += v(nv*n+vj)/wIJP1;
//            if (nodeCoarse.k > 0)
//              bCoarse(nv*nodeCoarse.n-1+vj)  += v(nv*n+vj)/wIJM1;
//            if (nodeCoarse.j < NY-1)
//              bCoarse(nv*nodeCoarse.n+NX+vj) += v(nv*n+vj)/wIP1J;
//            if (nodeCoarse.j > 0)
//              bCoarse(nv*nodeCoarse.n-NX+vj) += v(nv*n+vj)/wIM1J;
//          }
    }
  

  bandSolver.attachMat(coarseMat);
  error=bandSolver.solve(bCoarse,xCoarse);
  if (error)
    {
      cerr<<"Coarse Grid solve failed"<<endl;
      return(1);
    } 

  for (int n=0;n<nNodes;n++)
    {
      stencil.globalNode(n);
      /* opt computeNodeIndeces(stencil.anchor->i,
                         stencil.anchor->j,
                         stencil.anchor->k,
                         nodeCoarse,nodeSubdomain);*/
      const Location &nodeCoarse(coarseLocation(n));
      
      for (int vj=0;vj<nv;vj++)
        pv(nv*n+vj) = xCoarse(nv*nodeCoarse.n+vj);
//        calculateWeights(nodeSubdomain,nodeCoarse);
//        for (int vj=0;vj<nv;vj++)
//          {
//            pv(nv*n+vj) = xCoarse(nv*nodeCoarse.n+vj)*wIJ;
//            if (nodeCoarse.k < NX-1)
//              pv(nv*n+vj) += xCoarse(nv*nodeCoarse.n+1+vj)*wIJP1;
//            if (nodeCoarse.k > 0)
//              pv(nv*n+vj) += xCoarse(nv*nodeCoarse.n-1+vj)*wIJM1;
//            if (nodeCoarse.j < NY-1)
//              pv(nv*n+vj) +=xCoarse(nv*nodeCoarse.n+NX+vj)*wIP1J;
//            if (nodeCoarse.j > 0)
//              pv(nv*n+vj) += xCoarse(nv*nodeCoarse.n-NX+vj)*wIM1J;
//          }  
    }
  return 0;
}

template <class STENCIL, int nv>
bool TwoLevel<STENCIL,nv>::solveFinePrec(const Vec& v, Vec& pv)
{
  //opt  Location nodeCoarse,nodeSubdomain;
  for (int n=0;n<nNodes;n++)
    {
      stencil.globalNode(n);
      /* opt computeNodeIndeces(stencil.anchor->i,
                         stencil.anchor->j,
                         stencil.anchor->k,
                         nodeCoarse,nodeSubdomain);*/
      const Location &nodeCoarse(coarseLocation(n)), &nodeSubdomain(subdomainLocation(n));
       for (int vj=0;vj<nv;vj++)
        bSubD[nodeCoarse.n](nv*nodeSubdomain.n+vj) = v(nv*n+vj);
    }

  for (int sd=0;sd<nsd;sd++)
    {
      bandSolver.attachMat(subDomain[sd]);
      error=bandSolver.solve(bSubD[sd],xSubD[sd]);
      if (error)
        {
          cerr<<"Fine grid solve failed"<<endl;
        }
    }
  

  for (int n=0;n<nNodes;n++)
    {
      stencil.globalNode(n);
      /*opt      computeNodeIndeces(stencil.anchor->i,
                         stencil.anchor->j,
                         stencil.anchor->k,
                         nodeCoarse,nodeSubdomain);*/
      const Location &nodeCoarse(coarseLocation(n)), &nodeSubdomain(subdomainLocation(n));
      for (int vj=0;vj<nv;vj++)
        pv(nv*n+vj) = xSubD[nodeCoarse.n](nv*nodeSubdomain.n+vj);
    }
  return 0;
}

template <class STENCIL, int nv>
bool TwoLevel<STENCIL,nv>::apply(const Vec& v, Vec& pv)
{ 
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
          error=solveFinePrec(v,fineCorrection);
          if (error)
            {
              cerr<<"Error solving fine prec"<<endl;
              return 1;
            }
          if (USE_COARSE && MULTIPLICATIVE)
            {
              bool evalError=false;
              evalError=linearOperator->apply(fineCorrection,coarseCorrection);
              while(evalError)
                {
                  cerr<<"cutting back on fineCorrection in TwoLevel"<<endl;
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
              cerr<<"Error solving coarse matrix"<<endl;
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
                  cerr<<"cutting back coarse Correction in TwoLevel"<<endl;
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
          error=solveFinePrec(pv,fineCorrection);
          if (error)
            {
              cerr<<"Error solving fine prec"<<endl;
              return 1;
            }
        }
      //cout<<fineCorrection<<coarseCorrection;
      if (USE_COARSE)
        {
          //cout<<fineCorrection<<coarseCorrection;
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
  //cout<<wIM1J<<'\t'<<wIJM1<<'\t'<<wIJ<<'\t'<<wIJP1<<'\t'<<wIP1J<<endl;
}
}//Daetk
#endif
