#include "ModifiedNewton.h"

namespace Daetk 
{

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::max;
using std::min;

void ModifiedNewton::doChordIteration()
{CHORD_ITERATION = true;}

void ModifiedNewton::doFullNewton()
{CHORD_ITERATION = false;}
  
void ModifiedNewton::testResidual(real tol)
{TEST_RESIDUAL = true;resTol=tol;}
  
void ModifiedNewton::testConvergenceRate(bool flag)
{TEST_RATE = flag;}
  
void ModifiedNewton::linearSolverIsInexact()
{INEXACT_LINEAR_SOLVER=true;}
  
void ModifiedNewton::useLineSearch(int nls) 
{ nLineSearches=nls; USE_LINE_SEARCH = true; }
  
void ModifiedNewton::setLineSearchFact(real lsRedIn)
{ lsRedFact= lsRedIn; }
  
void ModifiedNewton::setLineSearchMethod(LineSearchType lsType)
{
  lineSearchMethod = lsType;
}

void ModifiedNewton::solveSubSystem(int start,int end,int stride,
				    int dimLS)
{SOLVE_SUB=true; VecIndex i(start,end); index=i; str=stride;}


ModifiedNewton::ModifiedNewton():
  SOLVE_SUB(false),
  USE_LINE_SEARCH(false),
  lineSearchMethod(poorMansLS),
  nLineSearches(0),
  lsRedFact(0.5),
  CHORD_ITERATION(true),
  INEXACT_LINEAR_SOLVER(false),
  TEST_RESIDUAL(false),
  TEST_RATE(true),
  RECOMPUTE_RATE(true),
  maxIterations(4),
  convergenceFactorp(0),
  s(100.0),
  nonlinearTolerance(.33),
  p(),
  residual(),
  weightedNorm(0),
  data(0),
  linearSolver(0),
  mnout("modNewt.txt")
{
  Tracer tr("ModifiedNewton::ModifiedNewton()");
}

ModifiedNewton::ModifiedNewton(LinearSolver& linearSolverIn,
                               VectorNorm& normIn, 
                               DataCollector& dataIn, int neq, 
                               real linearTol,
                               real nonlinearTol, int maxit):
  SOLVE_SUB(false),
  USE_LINE_SEARCH(false),
  lineSearchMethod(poorMansLS),
  nLineSearches(0),
  lsRedFact(0.5),
  CHORD_ITERATION(true),
  INEXACT_LINEAR_SOLVER(false),
  TEST_RESIDUAL(false),
  TEST_RATE(true),
  RECOMPUTE_RATE(true),
  maxIterations(maxit),
  convergenceFactorp(0),
  s(100.0), 
  nonlinearTolerance(nonlinearTol),
  roundOffTolerance(linearTol),
  p(neq),
  residual(neq),
  weightedNorm(&normIn),
  data(&dataIn),
  linearSolver(&linearSolverIn),
  mnout("modNewt.txt")
{
  Tracer tr("ModifiedNewton::ModifiedNewton(LinearSolver& linearSolverIn, WeightedVectorNorm& normIn, Data& dataIn, int neq)");
}

ModifiedNewton::~ModifiedNewton()
{
  Tracer tr("ModifiedNewton::~ModifiedNewton()");
}

void ModifiedNewton::setConvergenceFactor(const real& cf)
{
  convergenceFactorp = &cf;
}

void ModifiedNewton::computeRate()
{
  rate = pow(normOfCorrection/normOfFirstCorrection,1.0/real(iterations-1));
  if (rate < 1.0)
    s=rate/(1.0-rate);
  else
    {
      rate = 1.0;
      s=100;
    }
  RECOMPUTE_RATE = false;
}

void ModifiedNewton::recomputeConvergenceRate()
{
  RECOMPUTE_RATE=true;
}

bool ModifiedNewton::converged(VectorFunction& F)
{
  bool evalError=false, conv=false;
  if (TEST_RESIDUAL)
    {
      //      cout<<nrm2(F.value(evalError))<<'\t'<<r0<<endl;
      conv = (nrm2(F.value(evalError)) <= r0*resTol + resTol);
      if (evalError) //should do a line search or something here
        return false;
      else 
        return conv;
    }
  else if (RECOMPUTE_RATE)
    return false;
  else
    return s*normOfCorrection <= nonlinearTolerance;
}

bool ModifiedNewton::solve(Vec& correction,VectorFunction& F)
{

  //cout<<"in newton solve"<<endl<<flush;
  //This function will not exit in one iteration if RECOMPUTE_RATE is true, unlike DASPK's nonlinear
  //solver which forces the rate to be recomputed by setting it (the inverse rate) to a large value (100).
  //It is sometimes the case that even with s=100 the first correction is small enough to give convergence
  //while still not being below roundoff. This causes a slight difference in the way the two codes run.
 
  bool lsFailure,evalError=false;
  residual=F.value(evalError);

  if (evalError)
    {
      cerr<<"Predictor caused evaluation error; nonlinear solver returning failure"<<endl
          <<"The calling routine should have caught this"<<endl;
      return evalError;
    }
  if (TEST_RESIDUAL)
    r0 = nrm2(residual);

  iterations=0;
  correction=0.0;

  normOfPredictor=(*weightedNorm)(F.argument());

  if (!INEXACT_LINEAR_SOLVER)
    roundOffTolerance = 100.0 * normOfPredictor * MACHINE_EPSILON;
  //roundOffTolerance = linearTolerance otherwise

  //take one Newton Step and check to see if correction is extremely small

  data->nonlinearSolverIteration();
  ++iterations;
  if (*convergenceFactorp !=1.0 && CHORD_ITERATION)
    {
      //don't apply convergence factor to the entire vector
      //if we're solving a subsystem of the equations
      if(!SOLVE_SUB)
        residual*=(*convergenceFactorp);
      else
        {
          attache.attachToVecMulti(Vec::REF,residual,index);attache.setStrideMulti(str);
          attache*=(*convergenceFactorp);
        }
    }
  lsFailure=linearSolver->solve(residual,p);
  if (lsFailure)
    {
//       cerr<<"Linear Solver solve failure in Newton Iteration"<<endl;
      data->linearSolverFailure();
      return lsFailure;
    }

  F.correctArgument(p);


  normOfCorrection=(*weightedNorm)(p);

//uncomment for an alternavie subsystem solver...only check the norm of the subsystem correction  
//    if(!SOLVE_SUB)
//      normOfCorrection=(*weightedNorm)(p);
//    else
//      {
//        real sum=0.0;
//        const Vec& scaling(weightedNorm->getScaling());
//        for (int i=0;i<scaling.ldim()/2;i++)
//          sum+=scaling[str*i]*p[str*i]*scaling[str*i]*p[str*i];
//        normOfCorrection=sqrt(sum);
//      }
  
 
#ifndef USE_BLAS
  correction-=p;
#else
  axpy(-1.0,p,correction);
#endif
  if (normOfCorrection <= roundOffTolerance )
    {
      return false;
    }
  else
    {
      normOfFirstCorrection = normOfCorrection;
      if (!CHORD_ITERATION) // if this is full or inexact newton we can't use a rate from the past
        recomputeConvergenceRate();
      if ( converged(F) )
        return false;
    }
  real normOfOldCorrection=0.0, FnOld=r0; //force line search check on first iteration
  while ( iterations <  maxIterations )
    {      
      residual = F.value(evalError);
      //mnout<<residual<<F.argument();
      if (evalError)
        {
          cerr<<"feval in ModifiedNewton Iteration, exiting with failure"<<endl;
          return true;
        }
//       if ((normOfCorrection >= normOfOldCorrection && USE_LINE_SEARCH) || evalError)
//         {
//           real FnNew=FnOld;
//           if (evalError)
//               cerr<<"p or S out of range in Newton Iteration"<<endl;
//           else
//             FnNew = nrm2(F.value(evalError));
          
//           if (((FnNew > FnOld) && USE_LINE_SEARCH) || evalError)
//             {
//               int lin_it=0;
// #ifndef USE_BLAS
//               correction+=p;
// #else
//               axpy(1.0,p,correction);
// #endif         
              
//               while (( ((FnNew > FnOld) && USE_LINE_SEARCH) || evalError) && lin_it < 20)
//                 {
//                   std::cerr<<"In Line Search--not safe for Mass formulations yet"<<std::endl;
//                   F.unCorrect();
//                   data->lineSearch();

//                   //only apply linesearch to subsystem correction
//                   if(!SOLVE_SUB)
//                     p*=0.01;
//                   else
//                     {
//                       attache.attachToVecMulti(Vec::REF,p,index);attache.setStrideMulti(str);
//                       attache*=0.01;
//                     }
                  
//                   lin_it++;
//                   F.correctArgument(p);
//                   residual = F.value(evalError);
//                   if (evalError)
//                     cerr<<"p or S out of range in line search"<<endl;
//                   else
//                     FnNew = nrm2(residual);
//                   //cout<<FnNew<<endl;
//                 }
//               FnOld = FnNew;
// #ifndef USE_BLAS
//               correction-=p;
// #else
//               axpy(-1.0,p,correction);
// #endif         
//               //cerr<<lin_it<<" line searches"<<endl;
//               if (lin_it == 20)
//                 {
//                   F.unCorrect();
//                   cerr<<"Max # of line searches exceeded"<<endl;
//                   return true;
//                 }
//             }
//        }
      data->nonlinearSolverIteration();
      ++iterations;
      if (*convergenceFactorp !=1.0 && CHORD_ITERATION)
	{
          if(!SOLVE_SUB)
            residual*=(*convergenceFactorp);
          else
            {
              attache.attachToVecMulti(Vec::REF,residual,index);attache.setStrideMulti(str);
              attache*=(*convergenceFactorp);
            }
        }
      lsFailure=linearSolver->solve(residual,p);

      if (lsFailure)
	{
//           cerr<<"Linear Solver solve failure in Newton iteration"<<endl;
	  data->linearSolverFailure();
	  return lsFailure;
	}
      F.correctArgument(p);
      normOfOldCorrection=normOfCorrection;
      normOfCorrection=(*weightedNorm)(p);

//    //testing different method ---- REPLACE with normOfCorrection=(*weightedNorm)(p);
//    //hack
//    if(!SOLVE_SUB)
//      normOfCorrection=(*weightedNorm)(p);
//    else
//      {
//        real sum=0.0;
//        const Vec& scaling(weightedNorm->getScaling());
//        for (int i=0;i<scaling.ldim()/2;i++)
//          sum+=scaling[str*i]*p[str*i]*scaling[str*i]*p[str*i];
//        normOfCorrection=sqrt(sum);
//      }


#ifndef USE_BLAS
      correction-=p;
#else
      axpy(-1.0,p,correction);
#endif

      //decide whether to exit solver
      if (normOfCorrection <= roundOffTolerance )
        {
          return false;
        }
      computeRate();
      if ( TEST_RATE && rate > 0.9 )
        {  
//           cerr<<"nonlinear solver failed due to poor convergence rate "<<rate<<endl;
          return true;//convergence is too slow return failure
        }
      else if ( converged(F) )
        {
          return false;
        }
    }
//   cerr<<"nonlinear solver exceeded max iterations "<<iterations<<endl;
  return true;
}

//mwf===============================================
//probably don't really need this, but had some changes for ECDM
//==================================================
// cek stuff I cut out of the public archive's version of ModifiedNewton
///////////////////////////// start line search stuff
//mwf now put these in base class?
bool ModifiedNewtonMM::lineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
				  bool& evalFailed,
				  Vec& corr,
				  VectorFunction& F)
{
  bool LSFailed(false);
  switch (lineSearchMethod)
    {
    case cubicLS:
      {
	LSFailed = cubicLineSearch(yp,pp,Fp,Fnew,evalFailed,corr,F);
	break;
      }
    case armijoLS:
      {
	LSFailed = armijoLineSearch(yp,pp,Fp,Fnew,evalFailed,corr,F);
	break;
      }
    case simpleLS:
      {
	LSFailed = simpleLineSearch(yp,pp,Fp,Fnew,evalFailed,corr,F);
	break;
      }
    default:
      {
	LSFailed = poorMansLineSearch(yp,pp,Fp,Fnew,evalFailed,corr,F);
	break;
      }
    }
  return LSFailed;
}
//mwf try to use num rec. line search
bool ModifiedNewtonMM::cubicLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
				     bool& evalFailed,
				     Vec& corr,
				     VectorFunction& F)
{

#ifndef DEBUG
#define DEBUG
#define DEBUG_LOCAL
#endif
  /////line search stuff 
  //already comes in with full Newton step attempted
#ifdef DEBUG
  cout<<"entering cubic line search, evalFailed = "<<evalFailed<<endl;
  cerr<<"entering cubic line search, evalFailed = "<<evalFailed<<endl;
#endif
  if (USE_LINE_SEARCH)
    { 
      //mwf try and force first step to take line search for now?
      //Fp is F(yp)
      //go ahead and try to scale step?
//        real stepMax = 1.0;

      real normCorrection = nrm2(pp);
      pp *=1.0/normCorrection;

      real yTol=1.0e-7;
      real alpha=1.0e-4;
      
      //slope is grad(F*F/2)*p but should be -Fp*Fp
      
      real FnOld = 0.5*dot(Fp,Fp);
      real slope = -2.0*FnOld;

      //compute minimum step factor
      real test(0.0),temp(0.0);
      int nn = pp.size();
      for (int i=0; i < nn; i++)
	{
	  temp = fabs(pp(i))/std::max(fabs(yp(i)),1.0);
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
      int lin_it=0;
      real tmplam(0.0),a(0.0),b(0.0),rhs1(0.0),rhs2(0.0);
      real FnOld2(0.0),FnNew2(0.0),lam2(0.0),disc(0.0);

      //mwf wastes an evaluation but lets me know if 
      real FnNew= 0.5*dot(Fnew,Fnew);

#ifdef DEBUG
      cout<<"entering cubic line search, FnOld="<<FnOld<<endl;
      cerr<<"entering cubic line search, FnOld="<<FnOld<<endl;
      cout<<"entering cubic line search, FnNew="<<FnNew<<endl;
      cerr<<"entering cubic line search, FnNew="<<FnNew<<endl;
#endif
      while (lin_it < nLineSearches &&  lam >= lamin && 
	     (FnNew > FnOld + alpha*lam*slope || evalFailed ))
	{
	  //should I uncorrect argument before I start
	  //new line search?
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
	  //should I recalculate FnOld?
	  //is FnOld2 to be set to FnNew or something?
	  FnOld2=FnOld; 
#ifdef DEBUG
	  cout <<"in cubic line search, tmplam="<<tmplam<<endl;
	  cerr <<"in cubic line search, tmplam="<<tmplam<<endl;
#endif
	  //mwf try something more aggressive?
	  //mwf changed to 1.0e-3 from 1.0e-1
	  lam=max(tmplam,1.0e-1*lam);
	  

#ifdef DEBUG
	  cout <<"in cubic line search, lam="<<lam<<endl;
	  cerr <<"in cubic line search, lam="<<lam<<endl;
#endif
	  pp*=lam;
	  F.correctArgument(pp);
#ifndef USE_BLAS
	  corr-=pp;
#else
	  axpy(-1.0,pp,corr);
#endif

	  Fnew = F.value(evalFailed);
	  if (!evalFailed)
	    {	   
	      FnNew = 0.5*dot(Fnew,Fnew);
#ifdef DEBUG
	      cout<<"in cubic line search, FnNew="<<FnNew<<endl;
	      cerr<<"in cubic line search, FnNew="<<FnNew<<endl;
#endif
	      cerr<<"+"<<endl;
	    }
	  else
	    {
	      //mwf what to do if evaluation fails? 
	      //mwf want to reject new iterate
	      cerr<<"evalFailed in line search"<<endl;
	      cerr<<"-"<<endl;
	    }

	  
	  
	  data->lineSearch();
	  lin_it++;

	} // end while 
      //Fp is residual in solve 
      Fp = Fnew;
      if (lin_it == nLineSearches)
	{
	  cerr <<"line search reached max"<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  //still try another newton iteration?
	  return true;
	  
	}
      if (lam <= lamin)
	{
	  cerr <<"line search lambda = "<<lam<<" below min= "
	       <<lamin<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	}
      
    }//end if for line search
  ////// end line search stuff
#ifdef DEBUG
      cout<<"successfully completed cubic line search?"<<endl;
      cerr<<"successfully completed cubic line search?"<<endl;
#endif
#ifdef DEBUG_LOCAL
#undef DEBUG
#endif
  return false;
}

//mwf try to use poor man's  line search again
bool ModifiedNewtonMM::simpleLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
				      bool& evalFailed,
				      Vec& corr,
				      VectorFunction& F)
{
  if (USE_LINE_SEARCH)
    {
      real FnOld = nrm2(Fp);
      real FnNew = FnOld;

      if (evalFailed)
	cerr<<"p or S out of range in Newton Iteration simpLS"<<endl;
      else
	FnNew = nrm2(Fnew);


      std::cout<<" in simpLS "
	       <<" FnNew = "<<FnNew<<" FnOld= "<<FnOld<<std::endl;

      int lin_it=0;

      if ((FnNew > FnOld) || evalFailed)
	{

#ifndef USE_BLAS
	  corr+=pp;
#else
	  axpy(1.0,p,corr);
#endif         

	  while (((FnNew > FnOld) || evalFailed) && lin_it < nLineSearches)
	    {
	      F.unCorrect();
	      data->lineSearch();

	      //only apply linesearch to subsystem correction
	      if(!SOLVE_SUB)
		pp*=lsRedFact;
	      else
		{
		  attache.attachToVecMulti(Vec::REF,pp,index);attache.setStrideMulti(str);
		  attache*=lsRedFact;
		}
                  
	      lin_it++;
	      F.correctArgument(pp);
	      residual = F.value(evalFailed);
	      if (evalFailed)
		cerr<<"p or S out of range in line search"<<endl;
	      else
		FnNew = nrm2(residual);
	      //cout<<FnNew<<endl;
	    }//end while (FnNew > FnOld)
	  FnOld = FnNew;
#ifndef USE_BLAS
	  corr-=pp;
#else
	  axpy(-1.0,pp,corr);
              
#endif         
	  std::cout<<lin_it<<" line searches"<<std::endl;
	  if (lin_it == nLineSearches)
	    {
	      F.unCorrect();
	      //mwf reset boundary conditions based on FnOld
	      residual = F.value(evalFailed);
	      cerr<<"Max # of line searches exceeded"
		  <<" last FnNew=  "<<FnNew<<endl;
	      return true;
	    }
	} //end original if (FnNew > FnOld)
    }// end USE_LINE_SEARCH


  return false;

}
//mwf try to use poor man's  line search again
bool ModifiedNewtonMM::poorMansLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
					  bool& evalFailed,
					  Vec& corr,
					  VectorFunction& F)
{
  bool localEvalError(false);
  /////line search stuff 
  //mwf turn on debugging?
#ifndef DEBUG
#define DEBUG
#define DEBUG_LOCAL
#endif
  //already comes in with full Newton step attempted
#ifdef DEBUG
    cout<<"entering poor man's line search, evalFailed = "<<evalFailed<<endl;
    cerr<<"entering poor man's line search, evalFailed = "<<evalFailed<<endl;
#endif
  if (USE_LINE_SEARCH)
    { 

      real alpha=1.0e-4;
      real lamin=1.0e-12;
      //real lamax=0.5;

      real FnOld = norm(Fp);
     

      real lam   = 1.0;

      //temporaries for line search      
      int lin_it=0;

      //mwf wastes an evaluation but lets me know
      Fnew = F.value(localEvalError);

      real FnNew(100.0*FnOld);
      if (!evalFailed)
	FnNew= norm(Fnew);

#ifdef DEBUG
      cout<<"entering po line search, localEvalError = "<<localEvalError<<endl;
      cout<<"entering po line search, evalFailed = "<<evalFailed<<endl;
      cerr<<"entering po line search, localEvalError = "<<localEvalError<<endl;
      cerr<<"entering po line search, evalFailed = "<<evalFailed<<endl;
      cout<<"entering po line search, FnOld="<<FnOld<<endl;
      cerr<<"entering po line search, FnOld="<<FnOld<<endl;
      cout<<"entering po line search, FnNew="<<FnNew<<endl;
      cerr<<"entering po line search, FnNew="<<FnNew<<endl;
#endif


      while (lin_it < nLineSearches &&  lam >= lamin && 
	     (FnNew >= (1.0-alpha*lam)*FnOld   || evalFailed) )
	{
	  //uncorrect argument before I start
	  //new line search
	  F.unCorrect();

	  //mwf should I check that FnOld is same as value now
#ifdef DEBUG
	  //cout<<"after uncorrecting norm(F.value())= "
	  //  <<norm(F.value(localEvalError))<<endl;
#endif

	  lam = lsRedFact*lam;

#ifdef DEBUG
	  //cout <<"in line search, lam="<<lam<<endl;
	  //cerr <<"in line search, lam="<<lam<<endl;
#endif
	  pp*=lam;
	  F.correctArgument(pp);
	  //mwf do I still need this?
#ifndef USE_BLAS
	  corr-=pp;
#else
	  axpy(-1.0,pp,corr);
#endif

	  Fnew = F.value(evalFailed);
	  if (!evalFailed)
	    {	   
	      FnNew = norm(Fnew);
#ifdef DEBUG
	      //cout<<"in line search, FnNew="<<FnNew<<endl;
	      //cerr<<"in line search, FnNew="<<FnNew<<endl;
#endif
	      cerr<<"+"<<endl;
	    }
	  else
	    {
	      //mwf what to do if evaluation fails? 
	      //mwf want to reject new iterate
	      //FnNew = 100.0*FnOld;
	      cerr<<"evalFailed in line search"<<endl;
	      cerr<<"-"<<endl;
	    }

#ifdef DEBUG
	  //cerr <<"suff decrease = "<<FnNew - (1.0-alpha*lam)*FnOld<<endl;
	  //cout <<"suff decrease = "<<FnNew - (1.0-alpha*lam)*FnOld<<endl;
#endif
  
	  
	  data->lineSearch();
	  lin_it++;
	} // end while 
      //Fp is residual in solve 
      Fp = Fnew;
      if (lin_it == nLineSearches)
	{
	  cerr <<"line search reached max"<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	  
	}
      if (lam <= lamin)
	{
	  cerr <<"line search lambda = "<<lam<<" below min= "
	       <<lamin<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	}
      
    }//end if for line search
  ////// end line search stuff
#ifdef DEBUG
  //cerr <<"leaving poor man's line search "<<endl;
  //cout <<"leaving poor man's line search "<<endl;
#endif

#ifdef DEBUG_LOCAL
#undef DEBUG
#endif
  return false;
}

//mwf try to armijo  line search 
bool ModifiedNewtonMM::armijoLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
					bool& evalFailed,
					Vec& corr,
					VectorFunction& F)
{
//    bool localEvalError(false);
  /////line search stuff 
  //mwf turn on debugging?
#ifndef DEBUG
#define DEBUG
#define DEBUG_LOCAL
#endif
  //already comes in with full Newton step attempted
#ifdef DEBUG
    cout<<"entering armijo routine, evalFailed = "<<evalFailed<<endl;
    cerr<<"entering armijo routine, evalFailed = "<<evalFailed<<endl;
#endif
  if (USE_LINE_SEARCH)
    { 
      //
      real suffDecrParam=1.0e-4;
      //
      real initialStepSize=1.0;
      //
      real scaleFactor=lsRedFact*initialStepSize;

      real scaleMin=1.0e-12;
      int lin_it=0;
      //original Newton Step?
      const Vec origStep(pp);

      //now minimize 0.5*F^{T}*F
      real FnOld = 0.5*dot(Fp,Fp);
     

      //mwf wastes an evaluation but lets me know if 
      real FnNew= 0.5*dot(Fnew,Fnew);

#ifdef DEBUG
      //cout<<"entering line search, FnOld="<<FnOld<<endl;
      cerr<<"entering line search, FnOld="<<FnOld<<endl;
      //cout<<"entering line search, FnNew="<<FnNew<<endl;
      cerr<<"entering line search, FnNew="<<FnNew<<endl;
#endif

      while (lin_it < nLineSearches && scaleFactor >= scaleMin && 
	     (FnNew >= (1.0-2.0*suffDecrParam*scaleFactor)*FnOld  
	      || evalFailed ))
	{
	  //uncorrect argument before I start
	  //new line search
	  F.unCorrect();

	  //mwf should I check that FnOld is same as value now
#ifdef DEBUG
//  	  cout<<"after uncorrecting norm(F.value())= "
//  	      <<norm(F.value(localEvalError))<<endl;
#endif

	  scaleFactor *= lsRedFact;
#ifdef DEBUG
//  	  cout <<"in line search, scaleFactor="<<scaleFactor<<endl;
//  	  cerr <<"in line search, scaleFactor="<<scaleFactor<<endl;
#endif


#ifdef DEBUG
//  	  cout <<"in line search, scaleFactor="<<scaleFactor<<endl;
//  	  cerr <<"in line search, scaleFactor="<<scaleFactor<<endl;
#endif
	  pp = origStep;
	  pp*= scaleFactor;
	  F.correctArgument(pp);
#ifndef USE_BLAS
	  corr-=pp;
#else
	  axpy(-1.0,pp,corr);
#endif

	  Fnew = F.value(evalFailed);
	  if (!evalFailed)
	    {	   
	      FnNew = 0.5*dot(Fnew,Fnew);
#ifdef DEBUG
  	      //cout<<"in line search, FnNew="<<FnNew<<endl;
  	      cerr<<"in line search, FnNew="<<FnNew<<endl;
#endif
	      cerr<<"+"<<endl;
	    }
	  else
	    {
	      //mwf what to do if evaluation fails? 
	      //mwf want to reject new iterate
	      FnNew = 100.0*FnOld;
	      cerr<<"evalFailed in armijo line search"<<endl;
	      cerr<<"-"<<endl;
	    }

	  
	  
	  data->lineSearch();
	  lin_it++;
	} // end while 
      //Fp is residual in solve 
      Fp = Fnew;
      if (lin_it == nLineSearches)
	{
	  cerr <<"line search reached max"<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	  
	}
      if (scaleFactor <= scaleMin)
	{
	  cerr <<"line search scaleFactor = "<<scaleFactor<<" below min= "
	       <<scaleMin<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	}
      
    }//end if for line search
  ////// end line search stuff
#ifdef DEBUG
  //cerr <<"leaving armijo line search "<<endl;
  //cout <<"leaving armijo line search "<<endl;
#endif
  return false;
}

ModifiedNewtonMM::ModifiedNewtonMM() : ModifiedNewton() 
{}

ModifiedNewtonMM::ModifiedNewtonMM(LinearSolver& linearSolverIn, VectorNorm& W,
				   DataCollector& dataIn, int neq, 
				   real lTol,real nlTol, int maxit):
  ModifiedNewton(linearSolverIn,W,
		 dataIn,neq,lTol,nlTol,maxit),
  attache2(),attache3(),resForLSolve(),pForLSolve()
{}

ModifiedNewtonMM::ModifiedNewtonMM(LinearSolver& linearSolverIn, VectorNorm& W,
				   DataCollector& dataIn, int neq,
				   VectorNorm& Wlin,
				   real lTol,real nlTol, int maxit):
  ModifiedNewton(linearSolverIn,W,
		   dataIn,neq,lTol,nlTol,maxit),
    attache2(),attache3(),resForLSolve(),pForLSolve(),
    weightedNormLinSys(&Wlin)
    {}

ModifiedNewtonMM::~ModifiedNewtonMM() 
{}

void ModifiedNewtonMM::solveSubSystem(int start,int end,int stride,
				      int dimLS)
{
  SOLVE_SUB=true; VecIndex i(start,end); index=i; str=stride;
  if (dimLS > 0)
    {
      resForLSolve.newsize(dimLS);
      pForLSolve.newsize(dimLS);
    }
}


//mwf try a different way to switch out right hand side for subsystems
bool ModifiedNewtonMM::solve(Vec& correction,VectorFunction& F)
{

  //cout<<"in newton solve"<<endl<<flush;
  //This function will not exit in one iteration if RECOMPUTE_RATE is true, unlike DASPK's nonlinear
  //solver which forces the rate to be recomputed by setting it (the inverse rate) to a large value (100).
  //It is sometimes the case that even with s=100 the first correction is small enough to give convergence
  //while still not being below roundoff. This causes a slight difference in the way the two codes run.

  //mwf try and get vectors for linear solves
  //mwf index may not be set for normal system?

  bool lsFailure,evalError=false;
  residual=F.value(evalError);

  if (evalError)
    {
      cerr<<"Predictor caused evaluation error; nonlinear solver returning failure"<<endl
          <<"The calling routine should have caught this"<<endl;
      return evalError;
    }
  if (TEST_RESIDUAL)
    r0 = nrm2(residual);

  iterations=0;
  correction=0.0;

  //mwf now see if I need to set weighted norm?
  Vec dummy = F.argument();
  weightedNorm->setWeight(dummy);
  normOfPredictor=(*weightedNorm)(F.argument());

  if (!INEXACT_LINEAR_SOLVER)
    roundOffTolerance = 100.0 * normOfPredictor * MACHINE_EPSILON;
  //roundOffTolerance = linearTolerance otherwise

  //take one Newton Step and check to see if correction is extremely small

  data->nonlinearSolverIteration();
  ++iterations;
  if (*convergenceFactorp !=1.0 && CHORD_ITERATION)
    {
      //don't apply convergence factor to the entire vector
      //if we're solving a subsystem of the equations
      if(!SOLVE_SUB)
        residual*=(*convergenceFactorp);
      else
        {
          attache.attachToVecMulti(Vec::REF,residual,index);attache.setStrideMulti(str);
          attache*=(*convergenceFactorp);
        }
    }
  //mwf see if Chris' new LinearSolver class fixes this  
  //mwf seems to work ok
    lsFailure=linearSolver->solve(residual,p);
//   if(!SOLVE_SUB)
//     lsFailure=linearSolver->solve(residual,p);
//   else
//     {
//       //mwf now try and set weight for linear solver too
//       Vec tempArg = F.argument();
//       attache3.attachToVecMulti(Vec::REF,tempArg,index);
//       attache3.setStrideMulti(str);
//       assert(weightedNormLinSys);
//       weightedNormLinSys->setWeight(attache3);

//       //mwf copy over residual to correction first
//       p = residual;
//       attache.attachToVecMulti(Vec::REF,residual,index);
//       attache.setStrideMulti(str);

//       for (int i=0;i<resForLSolve.getLocalHigh();i++)
//         //resForLSolve[i]=attache[str*i];
// 	resForLSolve[i]=attache[i];


//       attache2.attachToVecMulti(Vec::REF,p,index);
//       attache2.setStrideMulti(str);

//       for (int i=0;i<pForLSolve.getLocalHigh();i++)
//         //pForLSolve[i]=attache2[str*i];
// 	pForLSolve[i]=attache2[i];

//       lsFailure=linearSolver->solve(resForLSolve,pForLSolve);

//       //now load back?
//       for (int i=0;i<resForLSolve.getLocalHigh();i++)
//         //attache[str*i]=resForLSolve[i];
//         attache[i]=resForLSolve[i];

//       for (int i=0;i<pForLSolve.getLocalHigh();i++)
//         //attache2[str*i]=pForLSolve[i];
//         attache2[i]=pForLSolve[i];
//     }
       
  if (lsFailure)
    {
      cerr<<"Linear Solver solve failure in Newton Iteration"<<endl;
      data->linearSolverFailure();
      return lsFailure;
    }

  F.correctArgument(p);


  normOfCorrection=(*weightedNorm)(p);

//uncomment for an alternavie subsystem solver...only check the norm of the subsystem correction  
//    if(!SOLVE_SUB)
//      normOfCorrection=(*weightedNorm)(p);
//    else
//      {
//        real sum=0.0;
//        const Vec& scaling(weightedNorm->getScaling());
//        for (int i=0;i<scaling.ldim()/2;i++)
//          sum+=scaling[str*i]*p[str*i]*scaling[str*i]*p[str*i];
//        normOfCorrection=sqrt(sum);
//      }
  
 
#ifndef USE_BLAS
  correction-=p;
#else
  axpy(-1.0,p,correction);
#endif
  if (normOfCorrection <= roundOffTolerance )
    {
      return false;
    }
  else
    {
      normOfFirstCorrection = normOfCorrection;
      if (!CHORD_ITERATION) // if this is full or inexact newton we can't use a rate from the past
        recomputeConvergenceRate();
      if ( converged(F) )
        return false;
    }
  real normOfOldCorrection=0.0, FnOld=r0; //force line search check on first iteration
  while ( iterations <  maxIterations )
    {      
      residual = F.value(evalError);
      //mnout<<residual<<F.argument();
      if ((normOfCorrection >= normOfOldCorrection && USE_LINE_SEARCH) || evalError)
        {
          real FnNew=FnOld;
          if (evalError)
            cerr<<"p or S out of range in Newton Iteration"<<endl;
          else
            FnNew = nrm2(F.value(evalError));
          
          if ((FnNew > FnOld) || evalError)
            {
              int lin_it=0;
#ifndef USE_BLAS
              correction+=p;
#else
              axpy(1.0,p,correction);
#endif         
              
              while (((FnNew > FnOld) || evalError) && lin_it < 5)
                {
                  F.unCorrect();
                  data->lineSearch();

                  //only apply linesearch to subsystem correction
                  if(!SOLVE_SUB)
                    p*=0.01;
                  else
                    {
                      attache.attachToVecMulti(Vec::REF,p,index);attache.setStrideMulti(str);
                      attache*=0.01;
                    }
                  
                  lin_it++;
                  F.correctArgument(p);
                  residual = F.value(evalError);
                  if (evalError)
                    cerr<<"p or S out of range in line search"<<endl;
                  else
                    FnNew = nrm2(residual);
                  //cout<<FnNew<<endl;
                }
              FnOld = FnNew;
#ifndef USE_BLAS
	      correction-=p;
#else
	      axpy(-1.0,p,correction);
              
#endif         
              //cerr<<lin_it<<" line searches"<<endl;
              if (lin_it == 5)
                {
                  F.unCorrect();
                  cerr<<"Max # of line searches exceeded"<<endl;
                  return true;
                }
            }
        }
      data->nonlinearSolverIteration();
      ++iterations;
      if (*convergenceFactorp !=1.0 && CHORD_ITERATION)
	{
          if(!SOLVE_SUB)
            residual*=(*convergenceFactorp);
          else
            {
              attache.attachToVecMulti(Vec::REF,residual,index);attache.setStrideMulti(str);
              attache*=(*convergenceFactorp);
            }
        }

      //mwf check to see if chris new linear solver fixes this
      //mwf seems to work ok
      lsFailure=linearSolver->solve(residual,p);
//       if(!SOLVE_SUB)
// 	lsFailure=linearSolver->solve(residual,p);
//       else
// 	{
// 	  //mwf copy over residual to correction first
// 	  p = residual;

// 	  attache.attachToVecMulti(Vec::REF,residual,index);
// 	  attache.setStrideMulti(str);
// 	  for (int i=0;i<resForLSolve.getLocalHigh();i++)
// 	    //resForLSolve[i]=attache[str*i];
// 	    resForLSolve[i]=attache[i];
	  
	  
// 	  attache2.attachToVecMulti(Vec::REF,p,index);
// 	  attache2.setStrideMulti(str);
	  
// 	  for (int i=0;i<pForLSolve.getLocalHigh();i++)
//   	    //pForLSolve[i]=attache2[str*i];
// 	    pForLSolve[i]=attache2[i];
	  
// 	  lsFailure=linearSolver->solve(resForLSolve,pForLSolve);
	  
// 	  //now load back?
// 	  for (int i=0;i<resForLSolve.getLocalHigh();i++)
//   	    //attache[str*i]=resForLSolve[i];
// 	    attache[i]=resForLSolve[i];
	  
// 	  for (int i=0;i<pForLSolve.getLocalHigh();i++)
//   	    //attache2[str*i]=pForLSolve[i];
// 	    attache2[i]=pForLSolve[i];
	  
// 	}
      
      

      if (lsFailure)
	{
          cerr<<"Linear Solver solve failure in Newton iteration"<<endl;
	  data->linearSolverFailure();
	  return lsFailure;
	}
      F.correctArgument(p);
      normOfOldCorrection=normOfCorrection;
      normOfCorrection=(*weightedNorm)(p);

//    //testing different method ---- REPLACE with normOfCorrection=(*weightedNorm)(p);
//    //hack
//    if(!SOLVE_SUB)
//      normOfCorrection=(*weightedNorm)(p);
//    else
//      {
//        real sum=0.0;
//        const Vec& scaling(weightedNorm->getScaling());
//        for (int i=0;i<scaling.ldim()/2;i++)
//          sum+=scaling[str*i]*p[str*i]*scaling[str*i]*p[str*i];
//        normOfCorrection=sqrt(sum);
//      }


#ifndef USE_BLAS
      correction-=p;
#else
      axpy(-1.0,p,correction);
#endif

      //decide whether to exit solver
      if (normOfCorrection <= roundOffTolerance )
        {
          return false;
        }
      computeRate();
      if ( TEST_RATE && rate > 0.9 )
        {  
          //cerr<<"nonlinear solver failed due to poor convergence rate "<<rate<<endl;
          return true;//convergence is too slow return failure
        }
      else if ( converged(F) )
        {
          return false;
        }
    }
  //  cerr<<"nonlinear solver exceeded max iterations "<<iterations<<endl;
  return true;
}

//mwf==================================================
//mwf try to switch to see what types of line search to use
//cek moved some stuff from the declaration in order to prevent implicit inlining
 
ModifiedNewtonSS::ModifiedNewtonSS() : ModifiedNewton() {}

ModifiedNewtonSS::ModifiedNewtonSS(LinearSolver& linearSolverIn, Jacobian& jacIn,VectorNorm& W,
                   DataCollector& dataIn, int neq, real lTol, 
 		   real nlTol, int maxit, 
		   LineSearchType lsType):
  ModifiedNewton(linearSolverIn,W,dataIn,neq,lTol,nlTol,maxit),
  jac(&jacIn)
{ lineSearchMethod = lsType; }
  
ModifiedNewtonSS::~ModifiedNewtonSS() {}

bool ModifiedNewtonSS::lineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
				  bool& evalFailed,
				  Vec& corr,
				  VectorFunction& F)
{
  bool LSFailed(false);
  switch (lineSearchMethod)
    {
    case cubicLS:
      {
	LSFailed = cubicLineSearch(yp,pp,Fp,Fnew,evalFailed,corr,F);
	break;
      }
    case armijoLS:
      {
	LSFailed = armijoLineSearch(yp,pp,Fp,Fnew,evalFailed,corr,F);
	break;
      }
    default:
      {
	LSFailed = poorMansLineSearch(yp,pp,Fp,Fnew,evalFailed,corr,F);
	break;
      }
    }
  return LSFailed;
}
//mwf try to use num rec. line search
bool ModifiedNewtonSS::cubicLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
				       bool& evalFailed,
				       Vec& corr,
				       VectorFunction& F)
{
#ifndef DEBUG
#define DEBUG
#define DEBUG_LOCAL
#endif
  /////line search stuff 
  //already comes in with full Newton step attempted
#ifdef DEBUG
  cout<<"entering cubic line search, evalFailed = "<<evalFailed<<endl;
  cerr<<"entering cubic line search, evalFailed = "<<evalFailed<<endl;
#endif
  if (USE_LINE_SEARCH)
    { 
      //mwf try and force first step to take line search for now?
      //Fp is F(yp)
      //go ahead and try to scale step?
//        real stepMax = 1.0;

      real normCorrection = nrm2(pp);
      pp *=1.0/normCorrection;

      real yTol=1.0e-7;
      real alpha=1.0e-4;
      
      //slope is grad(F*F/2)*p but should be -Fp*Fp
      
      real FnOld = 0.5*dot(Fp,Fp);
      real slope = -2.0*FnOld;

      //compute minimum step factor
      real test(0.0),temp(0.0);
      int nn = pp.size();
      for (int i=0; i < nn; i++)
	{
	  temp = fabs(pp(i))/std::max(fabs(yp(i)),1.0);
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
      int lin_it=0;
      real tmplam(0.0),a(0.0),b(0.0),rhs1(0.0),rhs2(0.0);
      real FnOld2(0.0),FnNew2(0.0),lam2(0.0),disc(0.0);

      //mwf wastes an evaluation but lets me know if 
      real FnNew= 0.5*dot(Fnew,Fnew);

#ifdef DEBUG
      cout<<"entering cubic line search, FnOld="<<FnOld<<endl;
      cerr<<"entering cubic line search, FnOld="<<FnOld<<endl;
      cout<<"entering cubic line search, FnNew="<<FnNew<<endl;
      cerr<<"entering cubic line search, FnNew="<<FnNew<<endl;
#endif
      while (lin_it < nLineSearches &&  lam >= lamin && 
	     (FnNew > FnOld + alpha*lam*slope || evalFailed ))
	{
	  //should I uncorrect argument before I start
	  //new line search?
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
	  //should I recalculate FnOld?
	  //is FnOld2 to be set to FnNew or something?
	  FnOld2=FnOld; 
#ifdef DEBUG
	  cout <<"in cubic line search, tmplam="<<tmplam<<endl;
	  cerr <<"in cubic line search, tmplam="<<tmplam<<endl;
#endif
	  //mwf try something more aggressive?
	  //mwf changed to 1.0e-3 from 1.0e-1
	  lam=max(tmplam,1.0e-1*lam);
	  

#ifdef DEBUG
	  cout <<"in cubic line search, lam="<<lam<<endl;
	  cerr <<"in cubic line search, lam="<<lam<<endl;
#endif
	  pp*=lam;
	  F.correctArgument(pp);
#ifndef USE_BLAS
	  corr-=pp;
#else
	  axpy(-1.0,pp,corr);
#endif

	  Fnew = F.value(evalFailed);
	  if (!evalFailed)
	    {	   
	      FnNew = 0.5*dot(Fnew,Fnew);
#ifdef DEBUG
	      cout<<"in cubic line search, FnNew="<<FnNew<<endl;
	      cerr<<"in cubic line search, FnNew="<<FnNew<<endl;
#endif
	      cerr<<"+"<<endl;
	    }
	  else
	    {
	      //mwf what to do if evaluation fails? 
	      //mwf want to reject new iterate
	      cerr<<"evalFailed in line search"<<endl;
	      cerr<<"-"<<endl;
	    }

	  
	  
	  data->lineSearch();
	  lin_it++;

	} // end while 
      //Fp is residual in solve 
      Fp = Fnew;
      if (lin_it == nLineSearches)
	{
	  cerr <<"line search reached max"<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  //still try another newton iteration?
	  return true;
	  
	}
      if (lam <= lamin)
	{
	  cerr <<"line search lambda = "<<lam<<" below min= "
	       <<lamin<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	}
      
    }//end if for line search
  ////// end line search stuff
#ifdef DEBUG
      cout<<"successfully completed cubic line search?"<<endl;
      cerr<<"successfully completed cubic line search?"<<endl;
#endif
#ifdef DEBUG_LOCAL
#undef DEBUG
#endif
  return false;
}
//mwf try to use poor man's  line search again
bool ModifiedNewtonSS::poorMansLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
					  bool& evalFailed,
					  Vec& corr,
					  VectorFunction& F)
{
  bool localEvalError(false);
  /////line search stuff 
  //mwf turn on debugging?
#ifndef DEBUG
#define DEBUG
#define DEBUG_LOCAL
#endif
  //already comes in with full Newton step attempted
#ifdef DEBUG
    cout<<"entering poor man's line search, evalFailed = "<<evalFailed<<endl;
    cerr<<"entering poor man's line search, evalFailed = "<<evalFailed<<endl;
#endif
  if (USE_LINE_SEARCH)
    { 

      real alpha=1.0e-4;
      real lamin=1.0e-12;
      //real lamax=0.5;

      real FnOld = norm(Fp);
     

      real lam   = 1.0;

      //temporaries for line search      
      int lin_it=0;

      //mwf wastes an evaluation but lets me know
      Fnew = F.value(localEvalError);

      real FnNew(100.0*FnOld);
      if (!evalFailed)
	FnNew= norm(Fnew);

#ifdef DEBUG
      cout<<"entering po line search, localEvalError = "<<localEvalError<<endl;
      cout<<"entering po line search, evalFailed = "<<evalFailed<<endl;
      cerr<<"entering po line search, localEvalError = "<<localEvalError<<endl;
      cerr<<"entering po line search, evalFailed = "<<evalFailed<<endl;
      cout<<"entering po line search, FnOld="<<FnOld<<endl;
      cerr<<"entering po line search, FnOld="<<FnOld<<endl;
      cout<<"entering po line search, FnNew="<<FnNew<<endl;
      cerr<<"entering po line search, FnNew="<<FnNew<<endl;
#endif


      while (lin_it < nLineSearches &&  lam >= lamin && 
	     (FnNew >= (1.0-alpha*lam)*FnOld   || evalFailed) )
	{
	  //uncorrect argument before I start
	  //new line search
	  F.unCorrect();

	  //mwf should I check that FnOld is same as value now
#ifdef DEBUG
	  //cout<<"after uncorrecting norm(F.value())= "
	  //  <<norm(F.value(localEvalError))<<endl;
#endif

	  lam = lsRedFact*lam;

#ifdef DEBUG
	  //cout <<"in line search, lam="<<lam<<endl;
	  //cerr <<"in line search, lam="<<lam<<endl;
#endif
	  pp*=lam;
	  F.correctArgument(pp);
	  //mwf do I still need this?
#ifndef USE_BLAS
	  corr-=pp;
#else
	  axpy(-1.0,pp,corr);
#endif

	  Fnew = F.value(evalFailed);
	  if (!evalFailed)
	    {	   
	      FnNew = norm(Fnew);
#ifdef DEBUG
	      //cout<<"in line search, FnNew="<<FnNew<<endl;
	      //cerr<<"in line search, FnNew="<<FnNew<<endl;
#endif
	      cerr<<"+"<<endl;
	    }
	  else
	    {
	      //mwf what to do if evaluation fails? 
	      //mwf want to reject new iterate
	      //FnNew = 100.0*FnOld;
	      cerr<<"evalFailed in line search"<<endl;
	      cerr<<"-"<<endl;
	    }

#ifdef DEBUG
	  //cerr <<"suff decrease = "<<FnNew - (1.0-alpha*lam)*FnOld<<endl;
	  //cout <<"suff decrease = "<<FnNew - (1.0-alpha*lam)*FnOld<<endl;
#endif
  
	  
	  data->lineSearch();
	  lin_it++;
	} // end while 
      //Fp is residual in solve 
      Fp = Fnew;
      if (lin_it == nLineSearches)
	{
	  cerr <<"line search reached max"<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	  
	}
      if (lam <= lamin)
	{
	  cerr <<"line search lambda = "<<lam<<" below min= "
	       <<lamin<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	}
      
    }//end if for line search
  ////// end line search stuff
#ifdef DEBUG
  //cerr <<"leaving poor man's line search "<<endl;
  //cout <<"leaving poor man's line search "<<endl;
#endif

#ifdef DEBUG_LOCAL
#undef DEBUG
#endif
  return false;
}

//mwf try to armijo  line search 
bool ModifiedNewtonSS::armijoLineSearch(Vec& yp,Vec& pp,Vec& Fp,Vec& Fnew,
					bool& evalFailed,
					Vec& corr,
					VectorFunction& F)
{
//    bool localEvalError(false);
  /////line search stuff 
  //mwf turn on debugging?
#ifndef DEBUG
#define DEBUG
#define DEBUG_LOCAL
#endif
  //already comes in with full Newton step attempted
#ifdef DEBUG
    cout<<"entering armijo routine, evalFailed = "<<evalFailed<<endl;
    cerr<<"entering armijo routine, evalFailed = "<<evalFailed<<endl;
#endif
  if (USE_LINE_SEARCH)
    { 
      //
      real suffDecrParam=1.0e-4;
      //
      real initialStepSize=1.0;
      //
      real scaleFactor=lsRedFact*initialStepSize;

      real scaleMin=1.0e-12;
      int lin_it=0;
      //original Newton Step?
      const Vec origStep(pp);

      //now minimize 0.5*F^{T}*F
      real FnOld = 0.5*dot(Fp,Fp);
     

      //mwf wastes an evaluation but lets me know if 
      real FnNew= 0.5*dot(Fnew,Fnew);

#ifdef DEBUG
      //cout<<"entering line search, FnOld="<<FnOld<<endl;
      cerr<<"entering line search, FnOld="<<FnOld<<endl;
      //cout<<"entering line search, FnNew="<<FnNew<<endl;
      cerr<<"entering line search, FnNew="<<FnNew<<endl;
#endif

      while (lin_it < nLineSearches && scaleFactor >= scaleMin && 
	     (FnNew >= (1.0-2.0*suffDecrParam*scaleFactor)*FnOld  
	      || evalFailed ))
	{
	  //uncorrect argument before I start
	  //new line search
	  F.unCorrect();

	  //mwf should I check that FnOld is same as value now
#ifdef DEBUG
//  	  cout<<"after uncorrecting norm(F.value())= "
//  	      <<norm(F.value(localEvalError))<<endl;
#endif

	  scaleFactor *= lsRedFact;
#ifdef DEBUG
//  	  cout <<"in line search, scaleFactor="<<scaleFactor<<endl;
//  	  cerr <<"in line search, scaleFactor="<<scaleFactor<<endl;
#endif


#ifdef DEBUG
//  	  cout <<"in line search, scaleFactor="<<scaleFactor<<endl;
//  	  cerr <<"in line search, scaleFactor="<<scaleFactor<<endl;
#endif
	  pp = origStep;
	  pp*= scaleFactor;
	  F.correctArgument(pp);
#ifndef USE_BLAS
	  corr-=pp;
#else
	  axpy(-1.0,pp,corr);
#endif

	  Fnew = F.value(evalFailed);
	  if (!evalFailed)
	    {	   
	      FnNew = 0.5*dot(Fnew,Fnew);
#ifdef DEBUG
  	      //cout<<"in line search, FnNew="<<FnNew<<endl;
  	      cerr<<"in line search, FnNew="<<FnNew<<endl;
#endif
	      cerr<<"+"<<endl;
	    }
	  else
	    {
	      //mwf what to do if evaluation fails? 
	      //mwf want to reject new iterate
	      FnNew = 100.0*FnOld;
	      cerr<<"evalFailed in armijo line search"<<endl;
	      cerr<<"-"<<endl;
	    }

	  
	  
	  data->lineSearch();
	  lin_it++;
	} // end while 
      //Fp is residual in solve 
      Fp = Fnew;
      if (lin_it == nLineSearches)
	{
	  cerr <<"line search reached max"<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	  
	}
      if (scaleFactor <= scaleMin)
	{
	  cerr <<"line search scaleFactor = "<<scaleFactor<<" below min= "
	       <<scaleMin<<endl;
	  cerr <<"norm of correction= "<<nrm2(corr)<<endl;
	  return true;
	}
      
    }//end if for line search
  ////// end line search stuff
#ifdef DEBUG
  //cerr <<"leaving armijo line search "<<endl;
  //cout <<"leaving armijo line search "<<endl;
#endif

#ifdef DEBUG_LOCAL
#undef DEBUG
#endif
  return false;
}

//mwf just made a few changes here to use with vector function in elliptic
//mwf problem
bool ModifiedNewtonSS::solve(Vec& correction,VectorFunction& F)
{
  //This function will not exit in one iteration if RECOMPUTE_RATE is true, unlike DASPK's nonlinear
  //solver which forces the rate to be recomputed by setting it (the inverse rate) to a large value (100).
  //It is sometimes the case that even with s=100 the first correction is small enough to give convergence
  //while still not being below roundoff. This causes a slight difference in the way the two codes run.
 
  bool lsFailure;

  bool evalFailed(false);

  //assume that original value is ok?
  residual=F.value(evalFailed);
  //#ifdef DEBUG
  cerr<<"initial residual in ModifiedNewtonSS:nrm2(residual)= "
      <<nrm2(residual)<<endl;
  //#endif

  //for holding solution value for numerical jacobian
  Vec tmpVec(residual);
  Vec tmpArg= F.argument();

  if (evalFailed)
    {
      cerr <<"original evaluation of VectorFunction failed in solve"<<endl;
      return 1;
    }
  if (TEST_RESIDUAL)
    r0 = nrm2(residual);
  iterations=0;
  correction=0.0;

  normOfPredictor=(*weightedNorm)(F.argument());

  if (!INEXACT_LINEAR_SOLVER)
    roundOffTolerance = 100.0 * normOfPredictor * MACHINE_EPSILON;
  //roundOffTolerance = linearTolerance otherwise

  //take one Newton Step and check to see if correction is extremely small

  data->nonlinearSolverIteration();
  ++iterations;
//    if (*convergenceFactorp != 1.0 && CHORD_ITERATION)
//      {
//  #ifndef USE_BLAS
//        residual*=(*convergenceFactorp);
//  #else
//        scal(*convergenceFactorp,residual);
//  #endif
//      }
  p=0;

  //mwf try to see what's going on here
  //mwf need to be able to turn this on and off?
  if (!CHORD_ITERATION)
    {
      //mwf try to do numerical jacobian sometime
      //mwf this should be like FLCBDF::evaluatJacobian
      data->jacobianEvaluation();
      //mwf need to find a way to do this only if actually
      //mwf using numerical jacobian
      F.computeDeltaForJacobian();
      tmpVec = F.value(evalFailed);
      tmpArg = F.argument();
      bool jacEvalFailed = 
	jac->evaluate(tmpArg,tmpVec);
      if (jacEvalFailed || evalFailed)
	{
	  cerr<<"jacobian eval failed in ModifiedNewtonSS"<<endl;
	  return jacEvalFailed;
	}
      bool linSolveFailed = linearSolver->prepare();
      if (linSolveFailed)
	{
	  cerr<<"linearSolver->prepare failed in ModifiedNewtonSS"<<endl;
	  return linSolveFailed;
	}
      
    }
//    if (!CHORD_ITERATION)
//       {//mwf try to do numerical jacobian sometime
//        jac->evaluate(F.argument(),F.value(evalFailed));
//        linearSolver->prepare();
//      }

  lsFailure=linearSolver->solve(residual,p);
  if (lsFailure)
    {
      data->linearSolverFailure();
      cout <<"linear solver failure"<<endl;
      cerr <<"linear solver failure"<<endl;
      return lsFailure;
    }
  //mwf need this to try and do num rec. line search
  Vec argPrev = F.argument();
  Vec funPrev = residual;

#ifdef DEBUG
//    cout <<"newton step, norm initial residual="<<norm(residual)<<endl;
//    cerr <<"newton step, norm initial residual="<<norm(residual)<<endl;
#endif
  F.correctArgument(p);

  normOfCorrection=(*weightedNorm)(p);
  
 
#ifndef USE_BLAS
  correction-=p;
#else
  axpy(-1.0,p,correction);
#endif

  
  if (normOfCorrection <= roundOffTolerance )
    {
      //mwf 
#ifdef DEBUG
      cout <<"in ModifiedNewtonMM::solve, \n normOfCorrection= "
	   <<normOfCorrection<<endl;
      cout <<"in ModifiedNewtonMM::solve, \n l2 norm of Correction= "
	   <<norm(p)<<endl;
      cout <<"normOfCorrection <= roundOffTolerance"<<endl;
      cerr <<"normOfCorrection <= roundOffTolerance"<<endl;
#endif
      return 0;
    }
  else
    {
      normOfFirstCorrection = normOfCorrection;
      if (!CHORD_ITERATION) // if this is full or inexact newton we can't use a rate from the past
        recomputeConvergenceRate();
      if ( converged(F) )
	return 0;
      
    }
  //  real normOfOldCorrection=normOfCorrection;

  while ( iterations <  maxIterations )
    {      
      //this should be the value of F at the previous iterate
      funPrev = residual;
      //this should be value at next newton iterate
      //mwf if evalFails, take care of this in line search?
      residual = F.value(evalFailed);

      if (evalFailed)
	{
#ifdef DEBUG
	  cout <<"evalFailed in Newton iteration"<<endl;
	  cerr <<"evalFailed in Newton iteration"<<endl;
#endif
	  //should I uncorrect argument here?
	  //line search alg does this for me
	}
      bool lineSearchFailed = lineSearch(argPrev,p,funPrev,
					 residual,evalFailed,
					 correction,F);

      if (lineSearchFailed)
	{
	  cout<<"lineSearch Failed"<<endl;
	  cerr<<"lineSearch Failed"<<endl;
	  return 1;
	}

      if (evalFailed && !USE_LINE_SEARCH)
	{
	  cout<<"eval failed in newton iteration and no line search"<<endl;
	  cerr<<"eval failed in newton iteration and no line search"<<endl;
	  return 1;
	}  
      data->nonlinearSolverIteration();
      ++iterations;
//        if (*convergenceFactorp !=1.0 && CHORD_ITERATION)
//  	{
//            //tricky for mass method
//            Vec::UnitStrideIterator weighti=weightedNorm->getWeightBegin();
//            const real* weightEnd =  weightedNorm->getWeightEnd();
//            int i=0;
//            while (weighti < weightEnd)
//              {
//                residual(i)*=(*convergenceFactorp);
//                ++i; ++weighti;
//              }
//          }
      p=0;
      //mwf try to see what's going on here
      if (!CHORD_ITERATION)
	{
	  //mwf try to do numerical jacobian sometime
	  //mwf this should be like FLCBDF::evaluatJacobian
	  data->jacobianEvaluation();
	  //need to find a way to do this only if doing 
	  //numerical calculation
	  F.computeDeltaForJacobian();
	  tmpVec = F.value(evalFailed);
	  tmpArg = F.argument();

	  bool jacEvalFailed = 
	    jac->evaluate(tmpArg,tmpVec);
	  if (jacEvalFailed || evalFailed)
	    {
	      cerr<<"jacobian eval failed in ModifiedNewtonSS"<<endl;
	      return jacEvalFailed;
	    }
	  bool linSolveFailed = linearSolver->prepare();
	  if (linSolveFailed)
	    {
	      cerr<<"linearSolver->prepare failed in ModifiedNewtonSS"<<endl;
	      return linSolveFailed;
	    }
	}
//        //mwf try to see what's going on here
//        if (!CHORD_ITERATION)
//  	{
//            jac->evaluate(F.argument(),F.value(evalFailed));
//            linearSolver->prepare();
//  	}
      lsFailure=linearSolver->solve(residual,p);

      if (lsFailure)
	{
	  data->linearSolverFailure();
	  cout <<"linear solver failure"<<endl;
	  cerr <<"linear solver failure"<<endl;
	  return lsFailure;
	}
      //mwf try and hold onto previous argument for line search
      argPrev = F.argument();

      F.correctArgument(p);
      //normOfOldCorrection=normOfCorrection;
      normOfCorrection=(*weightedNorm)(p);


#ifndef USE_BLAS
      correction-=p;
#else
      axpy(-1.0,p,correction);
#endif

      //decide whether to exit solver
      if (normOfCorrection <= roundOffTolerance )
        {
#ifdef DEBUG
	  cout <<"in ModifiedNewtonMM::solve, \n normOfCorrection= "
	       <<normOfCorrection<<endl;
	  cout <<"in ModifiedNewtonMM::solve, \n l2 norm of Correction= "
	       <<norm(p)<<endl;
	  cout <<"normOfCorrection <= roundOffTolerance"<<endl;
	  cerr <<"normOfCorrection <= roundOffTolerance"<<endl;
#endif
          return false;
        }
      computeRate();
      if ( TEST_RATE && rate > 0.9 )
        {  
	  //mwf
	  cerr<<"convergence test failed"<<endl;
          return 1;//convergence is too slow return failure
        }
      else if ( converged(F) )
        {
          return false;
        }
    }
  cerr<<"nonlinear solver exceeded max iterations "<<iterations<<endl;
  return true;
}

bool ModifiedNewtonSS::converged(VectorFunction& F)
{
  bool evalError=false, conv=false;
  if (TEST_RESIDUAL)
    {
      
      Vec tmp = F.value(evalError);
      //mwf was this
//        conv = (nrm2(tmp) <= r0*resTol);
//        //#ifdef DEBUG
//        cerr <<"in ModifiedNewtonSS::converged"<<endl;
//        cerr <<"nrm2(tmp)= "<<nrm2(tmp)<<"  "
      conv = (nrm2(tmp) <= r0*resTol + resTol);
      //#ifdef DEBUG
      cerr <<"in ModifiedNewtonSS::converged"<<endl;
      cerr <<"nrm2(tmp)= "<<nrm2(tmp)<<"  "
	   <<"resTol= "<<resTol<<"  "
	   <<"r0= "<<r0<<endl;
      //#endif
      if (!evalError)
	return  conv;
      else
	{
	  cout <<"evalError in converged = "<<evalError<<endl;
	  cerr <<"evalError in converged = "<<evalError<<endl;
	  return false;
	}
    }
  else if (RECOMPUTE_RATE)
    return false;
  else
    {
      //#ifdef DEBUG
      cerr <<"in ModifiedNewtonSS::converged"<<endl;
      cerr <<"s= "<<s<<"  "
	   <<"normOfCorrection= "<<normOfCorrection<<"  "
	   <<"nonlinearTolerance= "<<nonlinearTolerance<<endl;
      //#endif
      return s*normOfCorrection <= nonlinearTolerance;
    }
}

}//Daetk
