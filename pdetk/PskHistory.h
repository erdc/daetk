#ifndef PSKHISTORY_H
#define PSKHISTORY_H

#include <vector>
#include "Definitions.h"

namespace Daetk
{
class PskHistory
{
public:
  PskHistory(){}
  ~PskHistory(){}

  //Saturations
  
  inline real& lastWaterSat();
  inline real& lastNaplSat();
  inline real& lastAirSat();
  inline real& lastTotalSat();
  
  //Water Equivalent Capillary Heads
  
  inline real& lastAirWaterHead();
  inline real& lastAirNaplHead();
  inline real& lastNaplWaterHead();
  inline real& airWaterHead();
  inline real& airNaplHead();
  inline real& naplWaterHead();
  
  //Apparent Saturations
  
  inline real& waterAppSat();
  inline real& totalAppSat();
  
  inline real& lastWaterAppSat();
  inline real& lastTotalAppSat();
  
  //Important Historical Values
  
  inline real& histMin2pWaterAppSat(); 
  inline real& histMin3pWaterAppSat();
  inline real& histMinTotalAppSat();
  
  //Maximum Residual Effective Saturations of Nonwetting Phases
  
  inline real& maxResNaplEffSat();
  inline real& maxResAirTrappedByWaterEffSat();
  inline real& maxResAirTrappedByNaplEffSat();

  inline real& maxResAirOnMainNaplImbibition();
  inline real& maxResAirOnMainWaterImbibition();
  inline real& maxResNaplOnMainWaterImbibition();

  //irreducible water saturation
  inline real& irreducibleWaterSat(){ return *_irreducibleWaterSat;}
  //Reversal Points
  
  inline void pushDtoIReversalPoint_total(const real& S,const real& H,const real& S_d,const real& S_i);
  inline void pushDtoIReversalPoint_water(const real& S,const real& H,const real& S_d,const real& S_i);
  inline void pushItoDReversalPoint_total(const real& S,const real& H,const real& S_d,const real& S_i);
  inline void pushItoDReversalPoint_water(const real& S,const real& H,const real& S_d,const real& S_i);
  
  inline void popDtoIReversalPoint_total();
  inline void popDtoIReversalPoint_water();
  inline void popItoDReversalPoint_total();
  inline void popItoDReversalPoint_water();
  
  inline void lastDtoIReversalPoint_total(real& S, real& H, real& S_d, real& S_i);
  inline void lastDtoIReversalPoint_water(real& S, real& H, real& S_d, real& S_i);
  inline void lastItoDReversalPoint_total(real& S, real& H, real& S_d, real& S_i);
  inline void lastItoDReversalPoint_water(real& S, real& H, real& S_d, real& S_i);
  
  inline void nextToLastDtoIReversalPoint_total(real& S, real& H, real& S_d, real& S_i);
  inline void nextToLastDtoIReversalPoint_water(real& S, real& H, real& S_d, real& S_i);
  inline void nextToLastItoDReversalPoint_total(real& S, real& H, real& S_d, real& S_i);
  inline void nextToLastItoDReversalPoint_water(real& S, real& H, real& S_d, real& S_i);

  inline void getPreviousDtoIReversalPoint_total(int level,real& S, real& H, real& S_d, real& S_i);
  inline void getPreviousDtoIReversalPoint_water(int level,real& S, real& H, real& S_d, real& S_i);
  inline void getPreviousItoDReversalPoint_total(int level,real& S, real& H, real& S_d, real& S_i);
  inline void getPreviousItoDReversalPoint_water(int level,real& S, real& H, real& S_d, real& S_i);
  
  
  //Current Direction of Interface Movement (imbibing or draining)
  
  enum SaturationPath {DRAINAGE,IMBIBITION};
  
  SaturationPath& waterSaturationPath();
  SaturationPath& totalSaturationPath();
  
  //Flag for two or three phase system
  
  enum SystemState {TWOPHASE,THREEPHASE};
  
  SystemState& state(){return *_state;}
  SystemState& lastState(){return *_lastState;}
protected:
  real *_lastWaterSat,
    *_lastNaplSat,
    *_lastAirSat,
    *_lastTotalSat,
    *_lastAirWaterHead,
    *_lastAirNaplHead,
    *_lastNaplWaterHead,
    *_waterAppSat,
    *_totalAppSat,
    *_airWaterHead,
    *_airNaplHead,
    *_naplWaterHead,
    *_lastWaterAppSat,
    *_lastTotalAppSat,
    *_histMin2pWaterAppSat,
    *_histMin3pWaterAppSat,
    *_histMinTotalAppSat,
    *_maxResNaplEffSat,
    *_maxResAirByWaterEffSat,
    *_maxResAirByNaplEffSat,
    *_maxResAirOnMainNaplImbibition,
    *_maxResAirOnMainWaterImbibition,
    *_maxResNaplOnMainWaterImbibition,
    *_irreducibleWaterSat;

  struct ReversalPoint
  {
    real S,
      H,
      S_d,
      S_i;
  };
      
  std::vector<ReversalPoint> *_ItoDStack_water,
    *_DtoIStack_water,
    *_ItoDStack_total,
    *_DtoIStack_total;
  
  SaturationPath *_waterPath,*_totalPath;
  SystemState *_state,*_lastState;  
};



inline real& PskHistory::lastWaterSat()
{
  return *_lastWaterSat;
}
 
inline real& PskHistory::lastNaplSat()
{
  return *_lastNaplSat;
}

inline real& PskHistory::lastAirSat()
{
  return *_lastAirSat;
}

inline real& PskHistory::lastTotalSat()
{
  return *_lastTotalSat;
}

inline real& PskHistory::lastAirWaterHead()
{
  return *_lastAirWaterHead;
}

inline real& PskHistory::lastAirNaplHead()
{
  return *_lastAirNaplHead;
}

inline real& PskHistory::lastNaplWaterHead()
{
  return *_lastNaplWaterHead;
}

inline real& PskHistory::airWaterHead()
{
  return *_airWaterHead;
}

inline real& PskHistory::airNaplHead()
{
  return *_airNaplHead;
}

inline real& PskHistory::naplWaterHead()
{
  return *_naplWaterHead;
}

inline real& PskHistory::waterAppSat()
{
  return *_waterAppSat;
}

inline real& PskHistory::totalAppSat()
{
  return *_totalAppSat;
}

inline real& PskHistory::lastWaterAppSat()
{
  return *_lastWaterAppSat;
}

inline real& PskHistory::lastTotalAppSat()
{
  return *_lastTotalAppSat;
}

inline real& PskHistory::histMin2pWaterAppSat()
{
  return *_histMin2pWaterAppSat;
} 

inline real& PskHistory::histMin3pWaterAppSat()
{
  return *_histMin3pWaterAppSat;
}

inline real& PskHistory::histMinTotalAppSat()
{
  return *_histMinTotalAppSat;
}

inline real& PskHistory::maxResNaplEffSat()
{
  return *_maxResNaplEffSat;
}

inline real& PskHistory::maxResAirTrappedByWaterEffSat()
{
  return *_maxResAirByWaterEffSat;
}

inline real& PskHistory::maxResAirTrappedByNaplEffSat()
{
return *_maxResAirByNaplEffSat;
}
  
inline real& PskHistory::maxResAirOnMainNaplImbibition()
{
  return *_maxResAirOnMainNaplImbibition;
}
inline real& PskHistory::maxResAirOnMainWaterImbibition()
{
  return *_maxResAirOnMainWaterImbibition;
}

inline real& PskHistory::maxResNaplOnMainWaterImbibition()
{
  return *_maxResNaplOnMainWaterImbibition;
}

inline void PskHistory::pushDtoIReversalPoint_total(const real& S,const real& H, const real& S_d, const real& S_i)
{
  ReversalPoint p={S,H,S_d,S_i};
  _DtoIStack_total->push_back(p);
}

inline void PskHistory::pushDtoIReversalPoint_water(const real& S,const real& H, const real& S_d, const real& S_i)
{
  ReversalPoint p={S,H,S_d,S_i};
  _DtoIStack_water->push_back(p);
}

inline void PskHistory::pushItoDReversalPoint_total(const real& S,const real& H, const real& S_d, const real& S_i)
{
  ReversalPoint p={S,H,S_d,S_i};
  _ItoDStack_total->push_back(p);
}

inline void PskHistory::pushItoDReversalPoint_water(const real& S,const real& H, const real& S_d, const real& S_i)
{
  ReversalPoint p={S,H,S_d,S_i};
  _ItoDStack_water->push_back(p);
}

inline void PskHistory::popDtoIReversalPoint_total()
{
  if (_DtoIStack_total->size() > 1)
    _DtoIStack_total->pop_back();
}

inline void PskHistory::popDtoIReversalPoint_water()
{
  if (_DtoIStack_water->size() > 1)
    _DtoIStack_water->pop_back();
}

inline void PskHistory::popItoDReversalPoint_total()
{
  if (_ItoDStack_total->size() > 1)
    _ItoDStack_total->pop_back();
}

inline void PskHistory::popItoDReversalPoint_water()
{
  if (_ItoDStack_water->size() > 1)
    _ItoDStack_water->pop_back();
}

inline void PskHistory::lastDtoIReversalPoint_total(real& S, real& H, real& S_d, real& S_i)
{
  ReversalPoint p=_DtoIStack_total->back();
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::lastDtoIReversalPoint_water(real& S, real& H, real& S_d, real& S_i)
{
  ReversalPoint p = _DtoIStack_water->back();
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::lastItoDReversalPoint_total(real& S, real& H, real& S_d, real& S_i)
{  
  ReversalPoint p = _ItoDStack_total->back();
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::lastItoDReversalPoint_water(real& S, real& H, real& S_d, real& S_i)
{  
  ReversalPoint p = _ItoDStack_water->back();
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::nextToLastDtoIReversalPoint_total(real& S, real& H, real& S_d, real& S_i)
{
  ReversalPoint p = _DtoIStack_total->at(std::max(0,int(_DtoIStack_total->size())-2));
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::nextToLastDtoIReversalPoint_water(real& S, real& H, real& S_d, real& S_i)
{
  ReversalPoint p = _DtoIStack_water->at(std::max(0,int(_DtoIStack_water->size())-2));
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::nextToLastItoDReversalPoint_total(real& S, real& H, real& S_d, real& S_i)
{  
  ReversalPoint p = _ItoDStack_total->at(std::max(0,int(_ItoDStack_total->size())-2));
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::nextToLastItoDReversalPoint_water(real& S, real& H, real& S_d, real& S_i)
{  
  ReversalPoint p = _ItoDStack_water->at(std::max(0,int(_ItoDStack_water->size())-2));
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}
  
inline void PskHistory::getPreviousDtoIReversalPoint_total(int level,real& S, real& H, real& S_d, real& S_i)
{
  ReversalPoint p = _DtoIStack_total->at(std::max(0,int(_DtoIStack_total->size())-1-level));
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::getPreviousDtoIReversalPoint_water(int level,real& S, real& H, real& S_d, real& S_i)
{
  ReversalPoint p = _DtoIStack_water->at(std::max(0,int(_DtoIStack_water->size())-1-level));
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::getPreviousItoDReversalPoint_total(int level,real& S, real& H, real& S_d, real& S_i)
{  
  ReversalPoint p = _ItoDStack_total->at(std::max(0,int(_ItoDStack_total->size())-1-level));
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline void PskHistory::getPreviousItoDReversalPoint_water(int level,real& S, real& H, real& S_d, real& S_i)
{  
  ReversalPoint p = _ItoDStack_water->at(std::max(0,int(_ItoDStack_water->size())-1-level));
  S = p.S;
  H = p.H;
  S_d = p.S_d;
  S_i = p.S_i;
}

inline PskHistory::SaturationPath& PskHistory::waterSaturationPath()
{
  return *_waterPath;
}
inline PskHistory::SaturationPath& PskHistory::totalSaturationPath()
{
  return *_totalPath;
}
}//Daetk
#endif
