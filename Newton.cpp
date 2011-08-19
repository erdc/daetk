#include "Newton.h"

namespace Daetk 
{
// cek stuff I cut out of the public archive's version of ModifiedNewton
// ///////////////////////////// start line search stuff
// //mwf now put these in base class?
// bool ModifiedNewton::lineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
// 				  bool& evalFailed,
// 				  Vec& corr,
// 				  VectorFunction& F)
// {
//   bool LSFailed(false);
//   switch (lineSearchMethod)
//     {
//     case cubicLS:
//       {
// 	LSFailed = cubicLineSearch(yp,ppLS,Fp,Fnew,evalFailed,corr,F);
// 	break;
//       }
//     case armijoLS:
//       {
// 	LSFailed = armijoLineSearch(yp,ppLS,Fp,Fnew,evalFailed,corr,F);
// 	break;
//       }
//     case simpleLS:
//       {
// 	LSFailed = simpleLineSearch(yp,ppLS,Fp,Fnew,evalFailed,corr,F);
// 	break;
//       }
//     default:
//       {
// 	LSFailed = poorMansLineSearch(yp,ppLS,Fp,Fnew,evalFailed,corr,F);
// 	break;
//       }
//     }
//   return LSFailed;
// }
// //mwf try to use num rec. line search
// bool ModifiedNewton::cubicLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
// 				     bool& evalFailed,
// 				     Vec& corr,
// 				     VectorFunction& F)
// {

// #ifndef DEBUG
// #define DEBUG
// #define DEBUG_LOCAL
// #endif
//   /////line search stuff 
//   //already comes in with full Newton step attempted
// #ifdef DEBUG
//   cout<<"entering cubic line search, evalFailed = "<<evalFailed<<endl;
//   cerr<<"entering cubic line search, evalFailed = "<<evalFailed<<endl;
// #endif
//   if (USE_LINE_SEARCH)
//     { 
//       //mwf try and force first step to take line search for now?
//       //Fp is F(yp)
//       //go ahead and try to scale step?
// //        real stepMax = 1.0;

//       real normCorrection = nrm2(ppLS);
//       ppLS *=1.0/normCorrection;

//       real yTol=1.0e-7;
//       real alpha=1.0e-4;
      
//       //slope is grad(F*F/2)*p but should be -Fp*Fp
      
//       real FnOld = 0.5*dot(Fp,Fp);
//       real slope = -2.0*FnOld;

//       //compute minimum step factor
//       real test(0.0),temp(0.0);
//       int nn = ppLS.size();
//       for (int i=0; i < nn; i++)
// 	{
// 	  temp = fabs(ppLS(i))/std::max(fabs(yp(i)),1.0);
// 	  if (temp > test)
// 	    test = temp;
// 	}
//       real lamin = yTol/test;
//       real lam   = 1.0;
// #ifdef DEBUG
//       //  cout <<"in line search, lamin="<<lamin<<endl;
//       //  cerr <<"in line search, lamin="<<lamin<<endl;
// #endif

//       //temporaries for line search      
//       int lin_it=0;
//       real tmplam(0.0),a(0.0),b(0.0),rhs1(0.0),rhs2(0.0);
//       real FnOld2(0.0),FnNew2(0.0),lam2(0.0),disc(0.0);

//       //mwf wastes an evaluation but lets me know if 
//       real FnNew= 0.5*dot(Fnew,Fnew);

// #ifdef DEBUG
//       cout<<"entering cubic line search, FnOld="<<FnOld<<endl;
//       cerr<<"entering cubic line search, FnOld="<<FnOld<<endl;
//       cout<<"entering cubic line search, FnNew="<<FnNew<<endl;
//       cerr<<"entering cubic line search, FnNew="<<FnNew<<endl;
// #endif
//       while (lin_it < nLineSearches &&  lam >= lamin && 
// 	     (FnNew > FnOld + alpha*lam*slope || evalFailed ))
// 	{
// 	  //should I uncorrect argument before I start
// 	  //new line search?
// 	  F.unCorrect();

// 	  if (lam >= 1.0)
// 	    {
// 	      //first step
// 	      tmplam= -slope/(2.0*(FnNew-FnOld-slope));
// 	    }
// 	  else
// 	    {
// 	      rhs1 = FnNew-FnOld-lam*slope;
// 	      rhs2 = FnNew2-FnOld2-lam2*slope;
// 	      a=(rhs1/(lam*lam)-rhs2/(lam2*lam2))/(lam-lam2);
// 	      b=(-lam2*rhs1/(lam*lam)+lam*rhs2/(lam2*lam2))/(lam-lam2);
// 	      if (a==0.0)
// 		{
// 		  tmplam=-slope/(2.0*b);
// 		}
// 	      else
// 		{
// 		  disc=b*b-3.0*a*slope;
// 		  if (disc < 0.0)
// 		    {
// 		      cerr<<"roundoff problem in num rec. line search"<<endl;
// 		      return 1;
// 		    }
// 		  else
// 		    {
// 		      tmplam=(-b+sqrt(disc))/(3.0*a);
// 		    }
// 		}
// 	      if (tmplam > 0.5*lam)
// 		tmplam = 0.5*lam;
// 	    }//end portion for steps after first

// 	  lam2  =lam;
// 	  FnNew2=FnNew;
// 	  //should I recalculate FnOld?
// 	  //is FnOld2 to be set to FnNew or something?
// 	  FnOld2=FnOld; 
// #ifdef DEBUG
// 	  cout <<"in cubic line search, tmplam="<<tmplam<<endl;
// 	  cerr <<"in cubic line search, tmplam="<<tmplam<<endl;
// #endif
// 	  //mwf try something more aggressive?
// 	  //mwf changed to 1.0e-3 from 1.0e-1
// 	  lam=max(tmplam,1.0e-1*lam);
	  

// #ifdef DEBUG
// 	  cout <<"in cubic line search, lam="<<lam<<endl;
// 	  cerr <<"in cubic line search, lam="<<lam<<endl;
// #endif
// 	  ppLS*=lam;
// 	  F.correctArgument(ppLS);
// #ifndef USE_BLAS
// 	  corr-=ppLS;
// #else
// 	  axpy(-1.0,ppLS,corr);
// #endif

// 	  Fnew = F.value(evalFailed);
// 	  if (!evalFailed)
// 	    {	   
// 	      FnNew = 0.5*dot(Fnew,Fnew);
// #ifdef DEBUG
// 	      cout<<"in cubic line search, FnNew="<<FnNew<<endl;
// 	      cerr<<"in cubic line search, FnNew="<<FnNew<<endl;
// #endif
// 	      cerr<<"+"<<endl;
// 	    }
// 	  else
// 	    {
// 	      //mwf what to do if evaluation fails? 
// 	      //mwf want to reject new iterate
// 	      cerr<<"evalFailed in line search"<<endl;
// 	      cerr<<"-"<<endl;
// 	    }

	  
	  
// 	  data->lineSearch();
// 	  lin_it++;

// 	} // end while 
//       //Fp is residual in solve 
//       Fp = Fnew;
//       if (lin_it == nLineSearches)
// 	{
// 	  cerr <<"line search reached max"<<endl;
// 	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
// 	  //still try another newton iteration?
// 	  return true;
	  
// 	}
//       if (lam <= lamin)
// 	{
// 	  cerr <<"line search lambda = "<<lam<<" below min= "
// 	       <<lamin<<endl;
// 	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
// 	  return true;
// 	}
      
//     }//end if for line search
//   ////// end line search stuff
// #ifdef DEBUG
//       cout<<"successfully completed cubic line search?"<<endl;
//       cerr<<"successfully completed cubic line search?"<<endl;
// #endif
// #ifdef DEBUG_LOCAL
// #undef DEBUG
// #endif
//   return false;
// }

// //mwf try to use poor man's  line search again
// bool ModifiedNewton::simpleLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
// 				      bool& evalFailed,
// 				      Vec& corr,
// 				      VectorFunction& F)
// {
//   if (USE_LINE_SEARCH)
//     {
//       real FnOld = nrm2(Fp);
//       real FnNew = FnOld;

//       if (evalFailed)
// 	cerr<<"p or S out of range in Newton Iteration simpLS"<<endl;
//       else
// 	FnNew = nrm2(Fnew);


//       std::cout<<" in simpLS "
// 	       <<" FnNew = "<<FnNew<<" FnOld= "<<FnOld<<std::endl;

//       int lin_it=0;

//       if ((FnNew > FnOld) || evalFailed)
// 	{

// #ifndef USE_BLAS
// 	  corr+=ppLS;
// #else
// 	  axpy(1.0,p,corr);
// #endif         

// 	  while (((FnNew > FnOld) || evalFailed) && lin_it < nLineSearches)
// 	    {
// 	      F.unCorrect();
// 	      data->lineSearch();

// 	      //only appLSly linesearch to subsystem correction
// 	      if(!SOLVE_SUB)
// 		ppLS*=lsRedFact;
// 	      else
// 		{
// 		  attache.attachToVecMulti(Vec::REF,ppLS,index);attache.setStrideMulti(str);
// 		  attache*=lsRedFact;
// 		}
                  
// 	      lin_it++;
// 	      F.correctArgument(ppLS);
// 	      residual = F.value(evalFailed);
// 	      if (evalFailed)
// 		cerr<<"p or S out of range in line search"<<endl;
// 	      else
// 		FnNew = nrm2(residual);
// 	      //cout<<FnNew<<endl;
// 	    }//end while (FnNew > FnOld)
// 	  FnOld = FnNew;
// #ifndef USE_BLAS
// 	  corr-=ppLS;
// #else
// 	  axpy(-1.0,ppLS,corr);
              
// #endif         
// 	  std::cout<<lin_it<<" line searches"<<std::endl;
// 	  if (lin_it == nLineSearches)
// 	    {
// 	      F.unCorrect();
// 	      //mwf reset boundary conditions based on FnOld
// 	      residual = F.value(evalFailed);
// 	      cerr<<"Max # of line searches exceeded"
// 		  <<" last FnNew=  "<<FnNew<<endl;
// 	      return true;
// 	    }
// 	} //end original if (FnNew > FnOld)
//     }// end USE_LINE_SEARCH


//   return false;

// }
// //mwf try to use poor man's  line search again
// bool ModifiedNewton::poorMansLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
// 					  bool& evalFailed,
// 					  Vec& corr,
// 					  VectorFunction& F)
// {
//   bool localEvalError(false);
//   /////line search stuff 
//   //mwf turn on debugging?
// #ifndef DEBUG
// #define DEBUG
// #define DEBUG_LOCAL
// #endif
//   //already comes in with full Newton step attempted
// #ifdef DEBUG
//     cout<<"entering poor man's line search, evalFailed = "<<evalFailed<<endl;
//     cerr<<"entering poor man's line search, evalFailed = "<<evalFailed<<endl;
// #endif
//   if (USE_LINE_SEARCH)
//     { 

//       real alpha=1.0e-4;
//       real lamin=1.0e-12;
//       //real lamax=0.5;

//       real FnOld = norm(Fp);
     

//       real lam   = 1.0;

//       //temporaries for line search      
//       int lin_it=0;

//       //mwf wastes an evaluation but lets me know
//       Fnew = F.value(localEvalError);

//       real FnNew(100.0*FnOld);
//       if (!evalFailed)
// 	FnNew= norm(Fnew);

// #ifdef DEBUG
//       cout<<"entering po line search, localEvalError = "<<localEvalError<<endl;
//       cout<<"entering po line search, evalFailed = "<<evalFailed<<endl;
//       cerr<<"entering po line search, localEvalError = "<<localEvalError<<endl;
//       cerr<<"entering po line search, evalFailed = "<<evalFailed<<endl;
//       cout<<"entering po line search, FnOld="<<FnOld<<endl;
//       cerr<<"entering po line search, FnOld="<<FnOld<<endl;
//       cout<<"entering po line search, FnNew="<<FnNew<<endl;
//       cerr<<"entering po line search, FnNew="<<FnNew<<endl;
// #endif


//       while (lin_it < nLineSearches &&  lam >= lamin && 
// 	     (FnNew >= (1.0-alpha*lam)*FnOld   || evalFailed) )
// 	{
// 	  //uncorrect argument before I start
// 	  //new line search
// 	  F.unCorrect();

// 	  //mwf should I check that FnOld is same as value now
// #ifdef DEBUG
// 	  //cout<<"after uncorrecting norm(F.value())= "
// 	  //  <<norm(F.value(localEvalError))<<endl;
// #endif

// 	  lam = lsRedFact*lam;

// #ifdef DEBUG
// 	  //cout <<"in line search, lam="<<lam<<endl;
// 	  //cerr <<"in line search, lam="<<lam<<endl;
// #endif
// 	  ppLS*=lam;
// 	  F.correctArgument(ppLS);
// 	  //mwf do I still need this?
// #ifndef USE_BLAS
// 	  corr-=ppLS;
// #else
// 	  axpy(-1.0,ppLS,corr);
// #endif

// 	  Fnew = F.value(evalFailed);
// 	  if (!evalFailed)
// 	    {	   
// 	      FnNew = norm(Fnew);
// #ifdef DEBUG
// 	      //cout<<"in line search, FnNew="<<FnNew<<endl;
// 	      //cerr<<"in line search, FnNew="<<FnNew<<endl;
// #endif
// 	      cerr<<"+"<<endl;
// 	    }
// 	  else
// 	    {
// 	      //mwf what to do if evaluation fails? 
// 	      //mwf want to reject new iterate
// 	      //FnNew = 100.0*FnOld;
// 	      cerr<<"evalFailed in line search"<<endl;
// 	      cerr<<"-"<<endl;
// 	    }

// #ifdef DEBUG
// 	  //cerr <<"suff decrease = "<<FnNew - (1.0-alpha*lam)*FnOld<<endl;
// 	  //cout <<"suff decrease = "<<FnNew - (1.0-alpha*lam)*FnOld<<endl;
// #endif
  
	  
// 	  data->lineSearch();
// 	  lin_it++;
// 	} // end while 
//       //Fp is residual in solve 
//       Fp = Fnew;
//       if (lin_it == nLineSearches)
// 	{
// 	  cerr <<"line search reached max"<<endl;
// 	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
// 	  return true;
	  
// 	}
//       if (lam <= lamin)
// 	{
// 	  cerr <<"line search lambda = "<<lam<<" below min= "
// 	       <<lamin<<endl;
// 	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
// 	  return true;
// 	}
      
//     }//end if for line search
//   ////// end line search stuff
// #ifdef DEBUG
//   //cerr <<"leaving poor man's line search "<<endl;
//   //cout <<"leaving poor man's line search "<<endl;
// #endif

// #ifdef DEBUG_LOCAL
// #undef DEBUG
// #endif
//   return false;
// }

// //mwf try to armijo  line search 
// bool ModifiedNewton::armijoLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
// 					bool& evalFailed,
// 					Vec& corr,
// 					VectorFunction& F)
// {
// //    bool localEvalError(false);
//   /////line search stuff 
//   //mwf turn on debugging?
// #ifndef DEBUG
// #define DEBUG
// #define DEBUG_LOCAL
// #endif
//   //already comes in with full Newton step attempted
// #ifdef DEBUG
//     cout<<"entering armijo routine, evalFailed = "<<evalFailed<<endl;
//     cerr<<"entering armijo routine, evalFailed = "<<evalFailed<<endl;
// #endif
//   if (USE_LINE_SEARCH)
//     { 
//       //
//       real suffDecrParam=1.0e-4;
//       //
//       real initialStepSize=1.0;
//       //
//       real scaleFactor=lsRedFact*initialStepSize;

//       real scaleMin=1.0e-12;
//       int lin_it=0;
//       //original Newton Step?
//       const Vec origStep(ppLS);

//       //now minimize 0.5*F^{T}*F
//       real FnOld = 0.5*dot(Fp,Fp);
     

//       //mwf wastes an evaluation but lets me know if 
//       real FnNew= 0.5*dot(Fnew,Fnew);

// #ifdef DEBUG
//       //cout<<"entering line search, FnOld="<<FnOld<<endl;
//       cerr<<"entering line search, FnOld="<<FnOld<<endl;
//       //cout<<"entering line search, FnNew="<<FnNew<<endl;
//       cerr<<"entering line search, FnNew="<<FnNew<<endl;
// #endif

//       while (lin_it < nLineSearches && scaleFactor >= scaleMin && 
// 	     (FnNew >= (1.0-2.0*suffDecrParam*scaleFactor)*FnOld  
// 	      || evalFailed ))
// 	{
// 	  //uncorrect argument before I start
// 	  //new line search
// 	  F.unCorrect();

// 	  //mwf should I check that FnOld is same as value now
// #ifdef DEBUG
// //  	  cout<<"after uncorrecting norm(F.value())= "
// //  	      <<norm(F.value(localEvalError))<<endl;
// #endif

// 	  scaleFactor *= lsRedFact;
// #ifdef DEBUG
// //  	  cout <<"in line search, scaleFactor="<<scaleFactor<<endl;
// //  	  cerr <<"in line search, scaleFactor="<<scaleFactor<<endl;
// #endif


// #ifdef DEBUG
// //  	  cout <<"in line search, scaleFactor="<<scaleFactor<<endl;
// //  	  cerr <<"in line search, scaleFactor="<<scaleFactor<<endl;
// #endif
// 	  ppLS = origStep;
// 	  ppLS*= scaleFactor;
// 	  F.correctArgument(ppLS);
// #ifndef USE_BLAS
// 	  corr-=ppLS;
// #else
// 	  axpy(-1.0,ppLS,corr);
// #endif

// 	  Fnew = F.value(evalFailed);
// 	  if (!evalFailed)
// 	    {	   
// 	      FnNew = 0.5*dot(Fnew,Fnew);
// #ifdef DEBUG
//   	      //cout<<"in line search, FnNew="<<FnNew<<endl;
//   	      cerr<<"in line search, FnNew="<<FnNew<<endl;
// #endif
// 	      cerr<<"+"<<endl;
// 	    }
// 	  else
// 	    {
// 	      //mwf what to do if evaluation fails? 
// 	      //mwf want to reject new iterate
// 	      FnNew = 100.0*FnOld;
// 	      cerr<<"evalFailed in armijo line search"<<endl;
// 	      cerr<<"-"<<endl;
// 	    }

	  
	  
// 	  data->lineSearch();
// 	  lin_it++;
// 	} // end while 
//       //Fp is residual in solve 
//       Fp = Fnew;
//       if (lin_it == nLineSearches)
// 	{
// 	  cerr <<"line search reached max"<<endl;
// 	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
// 	  return true;
	  
// 	}
//       if (scaleFactor <= scaleMin)
// 	{
// 	  cerr <<"line search scaleFactor = "<<scaleFactor<<" below min= "
// 	       <<scaleMin<<endl;
// 	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
// 	  return true;
// 	}
      
//     }//end if for line search
//   ////// end line search stuff
// #ifdef DEBUG
//   //cerr <<"leaving armijo line search "<<endl;
//   //cout <<"leaving armijo line search "<<endl;
// #endif

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::max;
using std::min;

void Newton::testResidual(real tol)
{TEST_RESIDUAL = true;resTol=tol;}
  
void Newton::testConvergenceRate(bool flag)
{TEST_RATE = flag;}
  
void Newton::linearSolverIsInexact()
{INEXACT_LINEAR_SOLVER=true;}
  
void Newton::useLineSearch(int nls) 
{ nLineSearches=nls; USE_LINE_SEARCH = true; }
  
void Newton::setLineSearchFact(real lsRedIn)
{ lsRedFact= lsRedIn; }
  
void Newton::setLineSearchMethod(LineSearchType lsType)
{
  lineSearchMethod = lsType;
}


Newton::Newton():
  USE_LINE_SEARCH(false),
  LOG_STEPS(true),
  USE_PICARD_HYBRID(false),
  lineSearchMethod(poorMansLS),
  nLineSearches(0),
  lsRedFact(0.5),
  INEXACT_LINEAR_SOLVER(false),
  TEST_RESIDUAL(false),
  TEST_RATE(true),
  maxIterations(20),
  s(100.0),
  nonlinearTolerance(.33),
  p(),
  p0(),
  tmpArg(),
  tmpF(),
  argPrev(),
  funPrev(),
  pLS(),
  residual(),
  weightedNorm(0),
  data(0),
  linearSolver(0),
  mnout("modNewt.txt")
{
  Tracer tr("Newton::Newton()");
}

Newton::Newton(LinearSolver& linearSolverIn, Jacobian& jacIn,VectorNorm& normIn,
                   DataCollector& dataIn, int neq, real lTol, 
 		   real nlTol, int maxit, 
		   LineSearchType lsType):
  USE_LINE_SEARCH(false),
  LOG_STEPS(true),
  USE_PICARD_HYBRID(false),
  lineSearchMethod(poorMansLS),
  nLineSearches(0),
  lsRedFact(0.5),
  INEXACT_LINEAR_SOLVER(false),
  TEST_RESIDUAL(false),
  TEST_RATE(true),
  maxIterations(maxit),
  s(100.0), 
  nonlinearTolerance(nlTol),
  linearTolerance(lTol),
  p(neq),
  pLS(neq),
  p0(neq),
  tmpArg(neq),
  tmpF(neq),
  argPrev(neq),
  funPrev(neq),
  residual(neq),
  weightedNorm(&normIn),
  data(&dataIn),
  linearSolver(&linearSolverIn),
  mnout("modNewt.txt"),
  jac(&jacIn)
{
  Tracer tr("Newton::Newton(LinearSolver& linearSolverIn, Jacobian& jacIn,VectorNorm& W,DataCollector& dataIn, int neq, real lTol,real nlTol, int maxit,LineSearchType lsType)");
  lineSearchMethod = lsType; 
}

Newton::~Newton()
{
  Tracer tr("Newton::~Newton()");
}


void Newton::computeRate()
{
  
  if (iterations > 1)
    {
      rate = normOfCorrection/normOfLastCorrection;
      if (rate < 1.0)
        s=rate/(1.0-rate);
      else
        {
          rate = 1.0;
          s=100;
        }
    }
  else
    {
      rate=1.0;
      s=100;
    }
}

bool Newton::lineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
				  bool& evalFailed,
				  VectorFunction& F)
{
  bool LSFailed(false);
  switch (lineSearchMethod)
    {
    case cubicLS:
      {
	LSFailed = cubicLineSearch(yp,ppLS,pp,Fp,Fnew,evalFailed,F);
	break;
      }
    case armijoLS:
      {
	LSFailed = armijoLineSearch(yp,ppLS,pp,Fp,Fnew,evalFailed,F);
	break;
      }
    default:
      {
	LSFailed = poorMansLineSearch(yp,ppLS,pp,Fp,Fnew,evalFailed,F);
	break;
      }
    }
  return LSFailed;
}

//Line search from numerical recipes
bool Newton::cubicLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
                             bool& evalFailed,
                             VectorFunction& F)
{
  real normCorrection = nrm2(pp);
//   ppLS =1.0/normCorrection*pp;
  ppLS =pp;
  scal(1.0/normCorrection,ppLS);
  
  real yTol=1.0e-7;
  real alpha=1.0e-4;
  
  //slope is grad(F*F/2)*p but should be -Fp*Fp
  
  FnOld = 0.5*dot(Fp,Fp);
  real slope = -2.0*FnOld;
  
  //compute minimum step factor
  real test(0.0),temp(0.0);
  int nn = ppLS.size();
  for (int i=0; i < nn; i++)
    {
      temp = fabs(ppLS(i))/std::max(fabs(yp(i)),1.0);
      if (temp > test)
        test = temp;
    }
  real lamin = yTol/test;
  real lam   = 1.0;
#ifdef DEBUG
  //  cout <<"in line search, lamin="<<lamin<<endl;
  //  cerr <<"in line search, lamin="<<lamin<<endl;
#endif
  
  //temporaries for line search      
  lin_it=0;
  real tmplam(0.0),a(0.0),b(0.0),rhs1(0.0),rhs2(0.0);
  real FnOld2(0.0),FnNew2(0.0),lam2(0.0),disc(0.0);
  
  //mwf wastes an evaluation but lets me know if 
  FnNew= 0.5*dot(Fnew,Fnew);
  
#ifdef DEBUG
  cout<<"entering cubic line search, FnOld="<<FnOld<<endl;
  cerr<<"entering cubic line search, FnOld="<<FnOld<<endl;
  cout<<"entering cubic line search, FnNew="<<FnNew<<endl;
  cerr<<"entering cubic line search, FnNew="<<FnNew<<endl;
#endif
  while (lin_it < nLineSearches &&  
         lam >= lamin && 
         (FnNew > (FnOld + alpha*lam*slope) 
          || evalFailed ))
    {
      F.unCorrect();
      
      if (lam >= 1.0)
        {
          //first step
          tmplam= -slope/(2.0*(FnNew-FnOld-slope));
        }
      else
        {
          rhs1 = FnNew-FnOld-lam*slope;
          rhs2 = FnNew2-FnOld2-lam2*slope;
          a=(rhs1/(lam*lam)-rhs2/(lam2*lam2))/(lam-lam2);
          b=(-lam2*rhs1/(lam*lam)+lam*rhs2/(lam2*lam2))/(lam-lam2);
          if (a==0.0)
            {
              tmplam=-slope/(2.0*b);
            }
          else
            {
              disc=b*b-3.0*a*slope;
              if (disc < 0.0)
                {
                  cerr<<"roundoff problem in num rec. line search"<<endl;
                  return 1;
                }
              else
                {
                  tmplam=(-b+sqrt(disc))/(3.0*a);
                }
            }
          if (tmplam > 0.5*lam)
            tmplam = 0.5*lam;
        }//end portion for steps after first
      
      lam2  =lam;
      FnNew2=FnNew;
      FnOld2=FnOld; 

      lam=max(tmplam,1.0e-1*lam);
//       ppLS=lam*pp;
      ppLS=pp;
      scal(lam,ppLS);
      F.correctArgument(ppLS);
      Fnew = F.value(evalFailed);
      if (!evalFailed)
        {	   
          FnNew = 0.5*dot(Fnew,Fnew);
        }
      else
        {
//           cerr<<"evalFailed in line search"<<endl;
//           cerr<<"-"<<endl;
        }
      
      data->lineSearch();
      lin_it++;
      
    }

  Fp = Fnew;
  if (lin_it == nLineSearches)
    {
//       cerr <<"line search reached max"<<endl;
      return true;
      
    }
  if (lam <= lamin)
    {
//       cerr <<"line search lambda = "<<lam<<" below min= "
//            <<lamin<<endl;
      return true;
    }
  return false;
}
//mwf try to use poor man's  line search again
bool Newton::poorMansLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
					  bool& evalFailed,
					  VectorFunction& F)
{
  lin_it=0;
  real lam=1.0;
  FnOld = nrm2(Fp);
//   cout<<"FnOld***************** "<<FnOld<<endl;
  
  Fnew = F.value(evalError);
  FnNew=100.0*FnOld;
  if (!evalFailed)
    FnNew= nrm2(Fnew);
//   cout<<"FnNew***************** "<<FnNew<<endl;
  
  while (lin_it < nLineSearches &&  
         (FnNew >= FnOld   || 
          evalFailed) )
    {
      F.unCorrect();
      lam*=lsRedFact;
//       ppLS=lam*pp;
      ppLS=pp;
      scal(lam,ppLS);
      F.correctArgument(ppLS);
      Fnew = F.value(evalFailed);
      if (!evalFailed)
        {	   
          FnNew = nrm2(Fnew);
        }
      else
        {
//           cerr<<"evalFailed in line search"<<endl;
//           cerr<<"-"<<endl;
        }
      
//       cout<<"FnNew***************** "<<FnNew<<'\t'<<" FnOld*********** "<<FnOld<<endl;

      data->lineSearch();
      lin_it++;
    } 
  
  Fp = Fnew;
  if (lin_it == nLineSearches)
    {
//       cerr <<"line search reached max"<<endl;
      return true;
      
    }
  return false;
}

//mwf try to armijo  line search 
bool Newton::armijoLineSearch(Vec& yp,Vec& ppLS, const Vec& pp,Vec& Fp,Vec& Fnew,
					bool& evalFailed,
					VectorFunction& F)
{
  real alpha_ls=1.0e-4,
    lambda_c=1.0,
    lambda_m=lambda_c,
    sigma_c=lsRedFact,
    sigma_0=0.1,
    sigma_1=lsRedFact,
    f0,fc,fm;

  lin_it=0;

  FnOld = nrm2(Fp);
  f0=FnOld*FnOld;
//   cerr<<"FnOld***************** "<<FnOld<<endl;
  
  Fnew = F.value(evalFailed);
  FnNew=(100.0*FnOld);
  if (!evalFailed)
    FnNew= nrm2(Fnew);
  fc=FnNew*FnNew;
//   cerr<<"FnNew***************** "<<FnNew<<endl;
  
  real FnLast=FnNew;

  while (lin_it < nLineSearches && 
         (FnNew >= (1.0-alpha_ls*lambda_c)*FnOld  || 
          evalFailed ))
    {
//       cerr<<"FnNew********** "<<FnNew<<" FnOld **************** "<<FnOld<<endl;
      F.unCorrect();
      
      //use three-point parabolic model to compute the new sigma (from CTK's book)
      if (FnLast == FnNew || (lambda_c-lambda_m) == 0.0 )
        sigma_c=sigma_1;
      else if (evalFailed)
        sigma_c=sigma_1;
      else
        {
          if (lambda_c == 0.0 || lambda_m == 0.0)
            {
              std::cerr<<"Newton, problem in armijo line search, lambda=0, exit(1)"<<std::endl;
              exit(1);
            }
          real dp0=(1.0/(lambda_c - lambda_m))*(lambda_c*(fm-f0)/lambda_m - lambda_m*(fc-f0)/lambda_c),
            ddp0= (2.0/(lambda_c*lambda_m*(lambda_c-lambda_m)))*
            (lambda_m*(fc - f0) - lambda_c*(fm - f0));
          if (ddp0>0)
            {
              //lambda_c_new = -dp0/ddp0 = lambda_c_old *sigma_c therefore
              sigma_c = -dp0/(ddp0*lambda_c);
              if (sigma_c < sigma_0)
                sigma_c = sigma_0;
              else if (sigma_c > sigma_1)
                sigma_c = sigma_1;
            }
          else
            sigma_c=sigma_1;
        }
      lambda_m = lambda_c;
      lambda_c*=sigma_c;
      
//       std::cerr<<"sigma_c "<<sigma_c<<" lambda_c "<<lambda_c<<std::endl;

      //ppLS*=sigma_c;
//       ppLS= lambda_c*pp; //=lambda_c*pp_orig or sigma_c*ppLS
      ppLS= pp;
      scal(lambda_c,ppLS);

      F.correctArgument(ppLS);
      FnLast=FnNew;
      fm=fc;
      Fnew = F.value(evalFailed);
      if (!evalFailed)
        {	   
          FnNew = nrm2(Fnew);
        }
      else
        {
          FnNew = 100.0*FnOld;
//           cerr<<"evalFailed in armijo line search"<<endl;
//           cerr<<"-"<<endl;
        }
      fc=FnNew*FnNew;

//       cout<<"FnNew***************** "<<FnNew<<'\t'<<" FnOld*********** "<<FnOld<<endl;

      data->lineSearch();
      lin_it++;
    }
  
  Fp = Fnew;
  if (lin_it == nLineSearches)
    {
      cerr <<"line search reached max"<<endl;
      return true;
      
    }
  return false;
}

bool Newton::solve(Vec& correction,VectorFunction& F)
{
  if (LOG_STEPS)
    data->startUserStep();
  bool evalFailed,lineSearchFailed,linearSolverFailed;
  residual=F.value(evalFailed);
  if (evalFailed)
    {
      cerr<<"Initial guess in Newton is out of range, exiting"<<endl;
      exit(1);
    }

  weightedNorm->setWeight(F.argument());
  normOfInitialGuess =(*weightedNorm)(F.argument());
  
  if (TEST_RESIDUAL)
    {
      r0 = nrm2(residual);
    }

  iterations=0;
  correction=0.0;

  if (LOG_STEPS)
    data->stepTaken(iterations,lin_it,0.0,r0);
  
//   if (INEXACT_LINEAR_SOLVER)
//     roundOffTolerance = linearTolerance;
//   else
    roundOffTolerance = 100.0 * normOfInitialGuess * MACHINE_EPSILON;
  
  while ( iterations <  maxIterations )
    {      
      iterations++;
      data->nonlinearSolverIteration();
      lin_it=0;//line searches
      p=0;

      tmpArg=F.argument();
      tmpF=residual;

      data->jacobianEvaluation();
      F.computeDeltaForJacobian();

      bool jacEvalFailed = jac->evaluate(tmpArg,tmpF);
      if (jacEvalFailed || evalFailed)
        {
          cerr<<"jacobian eval failed in Newton"<<endl;
          return jacEvalFailed;
        }

      linearSolverFailed = linearSolver->prepare();
      if (linearSolverFailed)
        {
          if (LOG_STEPS)
            {
              data->endUserStep();
              data->linearSolverFailure();
              data->stepTaken(iterations,lin_it,0.0,FnNew);
            }
          else
            data->linearSolverFailure();

          cerr<<"linearSolver->prepare failed in Newton"<<endl;
          return linearSolverFailed;
        }
      
      linearSolverFailed=linearSolver->solve(residual,p);
      if (linearSolverFailed)
	{
          if (LOG_STEPS)
            {
              data->endUserStep();
              data->linearSolverFailure();
              data->stepTaken(iterations,lin_it,0.0,FnNew);
            }
          else
            data->linearSolverFailure();
            
	  cerr<<"linearSolver->solve failed in Newton"<<endl;
	  return linearSolverFailed;
	}
      if (iterations ==1) 
        p0=p;
      argPrev = F.argument();
      funPrev = residual;
      normOfLastCorrection=normOfCorrection;

      //correctArgument may change pLS so we need to save p in case we do a linesearch
      pLS=p;

      F.correctArgument(pLS);

      normOfCorrection=(*weightedNorm)(pLS);
      residual = F.value(evalFailed);

      if (USE_PICARD_HYBRID && iterations > 1 && nrm2(residual) < 1.0e-2*r0)
        {
          std::cerr<<"switching to Newton at iteration "<<iterations<<std::endl;
          USE_PICARD_HYBRID=false;
          F.usePicardApproximation(false);
          useLineSearch(1000);
          setLineSearchFact(0.5);
          setLineSearchMethod(armijoLS);
       }

      if (!evalFailed)
        {
          if (normOfCorrection <= roundOffTolerance)
            {
#ifndef USE_BLAS
              correction-=pLS;
#else
              axpy(-1.0,pLS,correction);
#endif
              cerr <<"Newton: normOfCorrection <= roundOffTolerance, exiting"<<endl
                   <<"Newton: normOfCorrection = "<<normOfCorrection<<" roundOffTolerance = "<<roundOffTolerance<<endl;
              FnNew=nrm2(F.value(evalError));
              if (LOG_STEPS)
                {
                  data->stepTaken(iterations,lin_it,0.0,FnNew);
                  data->endUserStep();
                }              
              return false;
            }
          computeRate();
          if ( TEST_RATE && rate > 0.9 )
            {  
              cerr<<"convergence is too slow in Newton, exiting with rate = "<<rate<<endl;
              FnNew=nrm2(F.value(evalError));
              if (LOG_STEPS)
                {
                  data->endUserStep();
                  data->nonlinearSolverFailure();
                  data->stepTaken(iterations,lin_it,0.0,FnNew);
                }
              else
                data->nonlinearSolverFailure();
              
              return 1;
            }
          else if ( converged(F) )
            {
#ifndef USE_BLAS
              correction-=pLS;
#else
              axpy(-1.0,pLS,correction);
#endif
              FnNew=nrm2(F.value(evalError));
              if (LOG_STEPS)
                {
                  data->stepTaken(iterations,lin_it,0.0,FnNew);
                  data->endUserStep();
                }
              return false;
            }
        }
      else
	{
//  	  cerr <<"evalFailed in Newton iteration"<<endl;
	}
      
      if(USE_LINE_SEARCH)
        {
          lineSearchFailed = lineSearch(argPrev,pLS,p,funPrev,residual,evalFailed,F);
          if (lineSearchFailed)
            {
              cerr<<"lineSearch Failed"<<endl;
              if (LOG_STEPS)
                {
                  data->endUserStep();
                  data->nonlinearSolverFailure();
                  data->stepTaken(iterations,lin_it,0.0,FnNew);
                }
              else
                data->nonlinearSolverFailure();
              return 1;
            }
        }
      else
        {
          if (evalFailed)
            {
              real lam=1.0;
              int lastChanceIts=0;
              while(evalFailed && lastChanceIts < 100)
                {
                  lastChanceIts++;
                  F.unCorrect();
                  lam*=0.5;
//                   pLS=lam*p;
                  pLS=p;
                  scal(lam,pLS);
                  F.correctArgument(pLS);
                  F.value(evalFailed);
                }
              lin_it=lastChanceIts;
              if(!evalFailed)
                {
                  residual=F.value(evalFailed);
                  FnNew=nrm2(F.value(evalError));
                }
              else
                {
                  cerr<<"eval failed in newton iteration and no line search"<<endl;
                  if (LOG_STEPS)
                    {
                      data->endUserStep();
                      data->nonlinearSolverFailure();
                      data->stepTaken(iterations,lin_it,0.0,FnNew);
                    }
                  else
                    data->nonlinearSolverFailure();
                  return 1;
                }
            }  
        }
      
      p=pLS;
#ifndef USE_BLAS
      correction-=p;
#else
      axpy(-1.0,p,correction);
#endif

      if (LOG_STEPS)
        {
          FnNew=nrm2(F.value(evalError));
          data->stepTaken(iterations,lin_it,0.0,FnNew);
        }
    }
  data->nonlinearSolverFailure();
  if (LOG_STEPS)
    {
      FnNew=nrm2(F.value(evalError));
      data->stepTaken(iterations,lin_it,0.0,FnNew);
    }
  cerr<<"nonlinear solver exceeded max iterations "<<iterations<<endl;
  return true;
}

bool Newton::converged(VectorFunction& F)
{
   if (TEST_RESIDUAL)
    {
      real r=nrm2(F.value(evalError));
      if (evalError)
        {
//           cerr<<"Newton: eval error in ::converged, shouldn't be here"<<endl;
          return  false;
        }
//       cerr<<"Newton: norm2(residual) after Newton correction = "<<r<<endl;
      return r <= (r0*resTol + resTol);
    }
  else if (iterations > 1)
    {
      return s*normOfCorrection <= nonlinearTolerance;
    }
  else
    return false;
}

}//Daetk
