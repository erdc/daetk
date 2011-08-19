#include "LinearSolver.h"

namespace Daetk
{

LinearSolver::LinearSolver():SOLVE_SUB(false){Tracer tr("LinearSolver::LinearSolver()");}

LinearSolver::~LinearSolver(){Tracer tr("LinearSolver::~LinearSolver()");}

void LinearSolver::checkTrueResidual(real resTol){}

void LinearSolver::printMatrices(const char* filename){}

  void LinearSolver::calculateCondition(DataCollector& d){std::cout<<"calCond not implemented"<<std::endl;}
  
void LinearSolver::solveSubSystem(int start,int end,int stride, int block){SOLVE_SUB=true; VecIndex i(start,end); index=i; str=stride;blk=block;}


LinearSolver::AttacheVec::AttacheVec(LinearSolver* lsIn, int dimIn):
  v_(dimIn),
  this_ls(lsIn)
{}

void LinearSolver::AttacheVec::attachToTarget(const Vec& target)
{
  if (this_ls->SOLVE_SUB)
    {
      attache.attachToVecMulti(Vec::REF,target,this_ls->index);
      attache.setBlockStrideMulti(this_ls->str,this_ls->blk);
      v_ = attache;
    }
  else
    {
      VecIndex all;
      v_.attachToVecMulti(Vec::REF,target,all);
    }
}

void LinearSolver::AttacheVec::restoreToTarget()
{
  if (this_ls->SOLVE_SUB)
    attache=v_;
}

}//Daetk

