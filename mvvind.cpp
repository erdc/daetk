#include "mvvind.h"

namespace Daetk
{
        MV_VecIndex::MV_VecIndex() : start_(0), end_(0), all_(1) {}

        MV_VecIndex::MV_VecIndex(unsigned int i1) :start_(i1), end_(i1), all_(0) {}

        MV_VecIndex::MV_VecIndex(unsigned int i1, unsigned int i2): start_(i1), end_(i2),
            all_(0)
        {
            assert(i1 <= i2);
        }

        MV_VecIndex::MV_VecIndex(const MV_VecIndex &s) : start_(s.start_), end_(s.end_), 
            all_(s.all_){}


        int MV_VecIndex::start() const { return (all_==1) ? 0 : start_;}

        int MV_VecIndex::end() const { return (all_ ==1) ? 0 : end_;}

        int MV_VecIndex::length() const { 
            return (all_==1) ? 0 : (end_-start_+1);}

        int MV_VecIndex::all() const { return all_; }

        MV_VecIndex& MV_VecIndex::operator=(const MV_VecIndex& I)
            { start_=I.start_; end_ = I.end_; return *this;}

        MV_VecIndex MV_VecIndex::operator+(int i)
            { return MV_VecIndex(start_ +i, end_ +i); }

        MV_VecIndex& MV_VecIndex::operator+=(int i)
            { start_ += i; end_ += i; return *this; }

        MV_VecIndex MV_VecIndex::operator-(int i)
            { return MV_VecIndex(start_ -i, end_ -i); }

        MV_VecIndex& MV_VecIndex::operator-=(int i)
            { start_ -= i; end_ -= i; return *this; }


}//Daetk

