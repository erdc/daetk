#ifndef PSHYSTERESIS_H
#define PSHYSTERESIS_H

#include "Psk2p.h"
#include "MualemVanGenuchten2p.h"
#include "VanGenuchten2p.h"
#include "PskHistory_Vec.h"

namespace Daetk 
{
class ParkerHysteresis : public MualemVanGenuchten2p
{
public:
  ParkerHysteresis():_naplTolerance(1e-4){};
  void setInitialConditions(const Vec& p_w)
  {
    for(i=0;i<m.dim();i++)
      {
        h.node(i);
        h.maxResAirOnMainWaterImbibition() = .20;
        h.irreducibleWaterSat() = thetaR[i]/thetaS[i];
        
        h.lastAirWaterHead() = -p_w[i];
        
        mainDrainage_water2p.setCapillaryHead(i,h.lastAirWaterHead());
        h.lastWaterSat() = mainDrainage_water2p.getSbar();
        h.lastWaterAppSat() = mainDrainage_water2p.getSbar();
        h.lastAirSat() = 1.0 - h.lastWaterSat();
        h.state() = PskHistory::TWOPHASE;
        h.lastState() = PskHistory::TWOPHASE;

        h.waterSaturationPath() = PskHistory::DRAINAGE;

        h.histMin2pWaterAppSat() = h.lastWaterSat();

        h.pushItoDReversalPoint_water(1.0,0.0,1.0,1.0);
        h.pushDtoIReversalPoint_water(0.0,1000,0.0,0.0);
              
        h.maxResAirTrappedByWaterEffSat() = (1.0 - h.histMin2pWaterAppSat()) / 
          ( 1.0 + (1.0/h.maxResAirOnMainWaterImbibition() - 1.0) * 
            (1.0 - h.histMin2pWaterAppSat()));
        
      }
    
  }
  void readParameters(ParameterDatabase& pd)
  {
    MualemVanGenuchten2p::readParameters(pd);
    h.newsize(m.dim());
    
    ParameterDatabase mit("imbTot.txt"),
      mdt("drainTot.txt"),
      miw2p("imbWater2p.txt"),
      mdw2p("drainWater2p.txt"),
      miw3p("imbWater3p.txt"),
      mdw3p("drainWater3p.txt");
    
    mainDrainage_total.readParameters(mdt);
    mainImbibition_total.readParameters(mit);
    mainDrainage_water2p.readParameters(mdw2p);
    mainImbibition_water2p.readParameters(miw2p);
    mainDrainage_water3p.readParameters(mdw3p);
    mainImbibition_water3p.readParameters(miw3p);
  }

  ~ParkerHysteresis(){};

  void setHeads(int node,real psiWIn, real psiNIn=0);
  void computeSaturationsFromHeads(const real& waterH,const real& naplH,
                                   const real& airH, real& waterS,real& naplS, 
                                   real& airS);

  inline const real& dWaterSat_dWaterHead();
  inline const real& dWaterSat_dNaplHead();
  inline const real& dWaterSat_dAirHead();
  inline const real& dNaplSat_dWaterHead();
  inline const real& dNaplSat_dNaplHead();
  inline const real& dNaplSat_dAirHead();
  inline const real& dAirSat_dWaterHead();
  inline const real& dAirSat_dNaplHead();
  inline const real& dAirSat_dAirHead();

  void updateHistory();

  VanGenuchten2p mainDrainage_total,
    mainImbibition_total,
    mainDrainage_water2p,
    mainImbibition_water2p,
    mainDrainage_water3p,
    mainImbibition_water3p;

private:	
  void computeTotalAppSat();
  void compute3pWaterAppSat();
  void compute2pWaterAppSat();
  void compute3pTrappedAirSat();
  void compute2pTrappedAirSat();
  void computeTrappedNaplSat();
  inline real& saturationFromDrainageScanningCurve();
  inline real& saturationFromImbibitionScanningCurve();
  inline real& dAppSat_dS_d();
  inline real& dAppSat_dS_i();
  real _S_i_of_h,
    _S_d_of_h,
    _S_i_of_h_id,
    _S_i_of_h_di,
    _S_d_of_h_di,
    _S_d_of_h_id,
    _totalAppSat,
    _waterAppSat,
    _naplTolerance,
    _airTrappedByWaterEffSat,
    _airTrappedByNaplEffSat,
    _trappedAirEffSat,
    _trappedNaplEffSat,
    _airTrappedInNaplByNaplEffSat,
    _airTrappedInNaplByWaterEffSat,
    _waterEffSat,
    _naplEffSat,
    _airEffSat,
    _totalTrappedNaplEffSat,
    _airWaterHead,
    _airNaplHead,
    _naplWaterHead,
    _h_id,
    _h_di,
    _S_id,
    _S_di,          
    _dTotalAppSat_dSatCurve,
    _dWaterAppSat_dSatCurve,
    _total_dSatCurve_dHead,
    _water_dSatCurve_dHead,
    _dAirTrappedByWaterEffSat_dWaterAppSat,
    _dTrappedAirEffSat_dTotalAppSat,
    _dAirTrappedByNaplEffSat_dWaterAppSat,
    _dAirTrappedByNaplEffSat_dTotalAppSat,
    _dTotalTrappedNaplEffSat_dWaterAppSat,
    _dTrappedNaplEffSat_dWaterAppSat,
    _dAirTrappedInNaplByWaterEffSat_dWaterAppSat,
    _dAirTrappedInNaplByNaplEffSat_dWaterAppSat,
    _dWaterSat_dAirHead,
    _dWaterSat_dWaterHead,
    _dWaterSat_dNaplHead,
    _dNaplSat_dAirHead,
    _dNaplSat_dWaterHead,
    _dNaplSat_dNaplHead,
    _dAirSat_dAirHead,
    _dAirSat_dWaterHead,
    _dAirSat_dNaplHead;

  PskHistory_Vec h;
};

inline real& ParkerHysteresis::saturationFromDrainageScanningCurve()
{
  static real saturation;
  if (_h_di == 1000 || _S_d_of_h_di == 1.0)//on main drainage or still saturated
    saturation = _S_d_of_h;
  else if (std::fabs(_S_d_of_h_id - _S_d_of_h_di) < SQRT_MACHINE_EPSILON)
    {
      saturation = _S_id;
      std::cout<<"problem in sat from drainage scan "<<_S_d_of_h_id<<'\t'<<_S_d_of_h_di<<std::endl;
    }
  else
    {
      saturation = (( _S_d_of_h - _S_d_of_h_di ) * (_S_id - _S_di)) /
        (_S_d_of_h_id - _S_d_of_h_di) + _S_di;
    }
  return saturation;
}

inline real& ParkerHysteresis::saturationFromImbibitionScanningCurve()
{
//    std::cout<<"scan "<<_S_i_of_h<<'\t'<<_S_i_of_h_id<<'\t'<<_S_di<<'\t'<<_S_id<<'\t'<<_S_i_of_h_di<<std::endl;
  static real saturation;
  if (_h_di == 1000 || _S_i_of_h_di == 1.0) //on main imbibition or already saturated on imbibition
    saturation = _S_i_of_h;
  else if (std::fabs(_S_i_of_h_di - _S_i_of_h_id) < SQRT_MACHINE_EPSILON)
    {
      saturation = _S_di;
      std::cout<<"problem in sat from imb scan "<<_S_i_of_h_di<<'\t'<<_S_i_of_h_id<<std::endl;
    }
  else
    {
      saturation = (( _S_i_of_h - _S_i_of_h_id ) * (_S_di - _S_id)) /
        (_S_i_of_h_di - _S_i_of_h_id) + _S_id;
    }
  return saturation;
}

inline real& ParkerHysteresis::dAppSat_dS_d()
{
  static real dSdS_d;
  if (_h_di == 1000 || _S_d_of_h_di == 1.0)//on main drainage or still saturated
    dSdS_d =  1.0;
  else if (std::fabs(_S_d_of_h_id - _S_d_of_h_di) < SQRT_MACHINE_EPSILON)
    {
      dSdS_d =  1.0;
      std::cout<<"problem in dSdS_d "<<_S_d_of_h_id<<'\t'<<_S_d_of_h_di<<std::endl;
    }
  else
    {
      dSdS_d = (_S_id - _S_di) / (_S_d_of_h_id - _S_d_of_h_di);
    }
  return dSdS_d;
}

inline real& ParkerHysteresis::dAppSat_dS_i()
{
  static real dSdS_i;
  if (_h_di == 1000 || _S_i_of_h_di == 1.0) //on main imbibition or already saturated on imbibition
    dSdS_i = 1.0;
  else if (std::fabs(_S_i_of_h_di - _S_i_of_h_id) < SQRT_MACHINE_EPSILON)
    {
      dSdS_i = 1.0;
      std::cout<<"problem in dSdS_i"<<_S_i_of_h_id<<'\t'<<_S_i_of_h_di<<std::endl;
    }
  else
    {
      dSdS_i = (_S_di - _S_id) / (_S_i_of_h_di - _S_i_of_h_id);
    }
  return dSdS_i;
}

inline const real& ParkerHysteresis::dWaterSat_dWaterHead()
{
  return _dWaterSat_dWaterHead;
}
 
inline const real& ParkerHysteresis::dWaterSat_dNaplHead()
{
  return _dWaterSat_dNaplHead;
}
  
inline const real& ParkerHysteresis::dWaterSat_dAirHead()
{
  return _dWaterSat_dAirHead;
}
  
inline const real& ParkerHysteresis::dNaplSat_dWaterHead()
{
  return _dNaplSat_dWaterHead;
}
  
inline const real& ParkerHysteresis::dNaplSat_dNaplHead()
{
  return _dNaplSat_dNaplHead;
}
  
inline const real& ParkerHysteresis::dNaplSat_dAirHead()
{
  return _dNaplSat_dAirHead;
}
  
inline const real& ParkerHysteresis::dAirSat_dWaterHead()
{
  return _dAirSat_dWaterHead;
}
  
inline const real& ParkerHysteresis::dAirSat_dNaplHead()
{
  return _dAirSat_dNaplHead;
}
  
inline const real& ParkerHysteresis::dAirSat_dAirHead()
{
  return _dAirSat_dAirHead;
}
}//Daetk  
#endif






