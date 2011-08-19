#include "ParkerHysteresis.h"
namespace Daetk
{
  using std::cout;
  using std::endl;

void  ParkerHysteresis::computeTotalAppSat()
{
//    if (h.totalSaturationPath() == PskHistory::IMBIBITION)
//      {
//        if (_airNaplHead > h.lastAirNaplHead()) //reversal: i to d
//  	{
//  	  _h_id = h.lastAirNaplHead();
//  	  _S_id = h.lastTotalAppSat();
//  	  h.lastDtoIReversalPoint_total(_S_di,_h_di);
	  
//  	  mainDrainage_total.setCapillaryHead(i,_h_id);
//            _S_d_of_h_id=mainDrainage_total.getSbar();
//  	  mainDrainage_total.setCapillaryHead(i,_h_di);
//            _S_d_of_h_di=mainDrainage_total.getSbar();
//  	  mainDrainage_total.setCapillaryHead(i,_airNaplHead);
//            _S_d_of_h=mainDrainage_total.getSbar();
//            _total_dSatCurve_dHead=mainDrainage_total.getDsBar_DpC();
	  
//  	  _totalAppSat = saturationFromDrainageScanningCurve();
//            _dTotalAppSat_dSatCurve = dAppSat_dS_d(); 
//  	}
//        else if (_airNaplHead < h.lastAirNaplHead()) //still imbibing
//  	{
//  	  h.lastItoDReversalPoint_total(_S_id,_h_id);
//            h.lastDtoIReversalPoint_total(_S_di,_h_di);
//  	  if ( _airNaplHead < _h_id ) //get on previous i scanning curve
//              {
//                h.nextToLastItoDReversalPoint_total(_S_id,_h_id);
//                h.nextToLastDtoIReversalPoint_total(_S_di,_h_di);
//              }
	  
//  	  mainImbibition_total.setCapillaryHead(i,_h_di);
//            _S_i_of_h_di=mainImbibition_total.getSbar();
//  	  mainImbibition_total.setCapillaryHead(i,_h_id);
//            _S_i_of_h_id=mainImbibition_total.getSbar();
//  	  mainImbibition_total.setCapillaryHead(i,_airNaplHead);
//            _S_i_of_h=mainImbibition_total.getSbar();
//            _total_dSatCurve_dHead=mainImbibition_total.getDsBar_DpC();
	  
//  	  _totalAppSat = saturationFromImbibitionScanningCurve();
//            _dTotalAppSat_dSatCurve = dAppSat_dS_i(); 
//  	}
//        else //airNaplHead == _lastAirNaplHead
//          {
//  	  mainImbibition_total.setCapillaryHead(i,_airNaplHead);
//            _S_i_of_h=mainImbibition_total.getSbar();
//            _total_dSatCurve_dHead=mainImbibition_total.getDsBar_DpC();
//            _totalAppSat = h.lastTotalAppSat();
//            _dTotalAppSat_dSatCurve = dAppSat_dS_i(); 
//          }
//      }
//    else	//drainage
//      {
//        if (_airNaplHead < h.lastAirNaplHead()) //reversal: d to i
//  	{
//  	  _h_di = h.lastAirNaplHead();
//  	  _S_di = h.lastTotalAppSat();
//  	  h.lastItoDReversalPoint_total(_S_id,_h_id);
	  
//  	  mainImbibition_total.setCapillaryHead(i,_h_di);
//            _S_i_of_h_di=mainImbibition_total.getSbar();
//  	  mainImbibition_total.setCapillaryHead(i,_h_id);
//            _S_i_of_h_id=mainImbibition_total.getSbar();
//  	  mainImbibition_total.setCapillaryHead(i,_airNaplHead);
//            _S_i_of_h=mainImbibition_total.getSbar();
//            _total_dSatCurve_dHead=mainImbibition_total.getDsBar_DpC();
	  
//  	  _totalAppSat = saturationFromImbibitionScanningCurve();
//            _dTotalAppSat_dSatCurve = dAppSat_dS_i(); 
//  	}
//        else if (_airNaplHead > h.lastAirNaplHead()) //still draining
//  	{
//  	  h.lastDtoIReversalPoint_total(_S_di,_h_di);
//            h.lastItoDReversalPoint_total(_S_id,_h_id);
//  	  if ( _airNaplHead > _h_di ) //get on previous d scanning curve
//              {
//                h.nextToLastDtoIReversalPoint_total(_S_di,_h_di);
//                h.nextToLastItoDReversalPoint_total(_S_id,_h_id);
//              }
	  
//  	  mainDrainage_total.setCapillaryHead(i,_h_di);
//            _S_d_of_h_di=mainDrainage_total.getSbar();
//  	  mainDrainage_total.setCapillaryHead(i,_h_id);
//            _S_d_of_h_id=mainDrainage_total.getSbar();
//  	  mainDrainage_total.setCapillaryHead(i,_airNaplHead);
//            _S_d_of_h=mainDrainage_total.getSbar();
//            _total_dSatCurve_dHead=mainDrainage_total.getDsBar_DpC();
	  
//  	  _totalAppSat = saturationFromDrainageScanningCurve();
//            _dTotalAppSat_dSatCurve = dAppSat_dS_d(); 
//  	}
//        else //airNaplHead == _lastAirNaplHead
//  	{
//  	  mainDrainage_total.setCapillaryHead(i,_airNaplHead);
//            _S_d_of_h=mainDrainage_total.getSbar();
//            _total_dSatCurve_dHead=mainDrainage_total.getDsBar_DpC();
//            _totalAppSat = h.lastTotalAppSat();
//            _dTotalAppSat_dSatCurve = dAppSat_dS_d(); 
//          }
//      }
}

void ParkerHysteresis::compute3pWaterAppSat()
{		
//    if (h.waterSaturationPath() == PskHistory::IMBIBITION)
//      {
//        if (_naplWaterHead  > h.lastNaplWaterHead() )	//reversal: i to d
//  	{
//  	  _h_id = h.lastNaplWaterHead();
//  	  _S_id = h.lastWaterAppSat();
//  	  h.lastDtoIReversalPoint_water(_S_di,_h_di);
	  
//  	  mainDrainage_water3p.setCapillaryHead(i,_h_id);
//            _S_d_of_h_id =mainDrainage_water3p.getSbar();
//  	  mainDrainage_water3p.setCapillaryHead(i,_h_di);
//            _S_d_of_h_di=mainDrainage_water3p.getSbar();
//  	  mainDrainage_water3p.setCapillaryHead(i,_naplWaterHead);
//            _S_d_of_h=mainDrainage_water3p.getSbar();
//            _water_dSatCurve_dHead = mainDrainage_water3p.getDsBar_DpC();
	  
//  	  _waterAppSat = saturationFromDrainageScanningCurve();
//            _dWaterAppSat_dSatCurve = dAppSat_dS_d(); 
//          }
//        else if (_naplWaterHead < h.lastNaplWaterHead()) //still imbibiing
//  	{
//  	  h.lastItoDReversalPoint_water(_S_id,_h_id);
//            h.lastDtoIReversalPoint_water(_S_di,_h_di);
//  	  if( _naplWaterHead < _h_id) //get on previous i scanning curve
//  	     {
//                 h.nextToLastItoDReversalPoint_water(_S_id,_h_id);
//                 h.nextToLastDtoIReversalPoint_water(_S_di,_h_di);
//               }

//  	  mainImbibition_water3p.setCapillaryHead(i,_h_di);
//            _S_i_of_h_di=mainImbibition_water3p.getSbar();
//  	  mainImbibition_water3p.setCapillaryHead(i,_h_id);
//            _S_i_of_h_id=mainImbibition_water3p.getSbar();
//  	  mainImbibition_water3p.setCapillaryHead(i,_naplWaterHead);
//            _S_i_of_h=mainImbibition_water3p.getSbar();
//            _water_dSatCurve_dHead = mainImbibition_water3p.getDsBar_DpC();

//  	  _waterAppSat = saturationFromImbibitionScanningCurve();
//            _dWaterAppSat_dSatCurve = dAppSat_dS_i(); 
//  	}
//        else //naplWaterHead == _lastNaplWaterHead
//  	{
//  	  mainImbibition_water3p.setCapillaryHead(i,_naplWaterHead);
//            _S_i_of_h=mainImbibition_water3p.getSbar();
//            _water_dSatCurve_dHead = mainImbibition_water3p.getDsBar_DpC();
//            _waterAppSat = h.lastWaterAppSat();
//            _dWaterAppSat_dSatCurve = dAppSat_dS_i(); 
//          }

//      }
//    else //draining
//      {		
//        if (_naplWaterHead < h.lastNaplWaterHead()) //reversal: d to i
//  	{
//  	  _h_di = h.lastNaplWaterHead();
//  	  _S_di = h.lastWaterAppSat();
//  	  h.lastItoDReversalPoint_water(_S_id,_h_id);
	  
//  	  mainImbibition_water3p.setCapillaryHead(i,_h_id);
//            _S_i_of_h_id=mainImbibition_water3p.getSbar();
//  	  mainImbibition_water3p.setCapillaryHead(i,_h_di);
//            _S_i_of_h_di=mainImbibition_water3p.getSbar();
//  	  mainImbibition_water3p.setCapillaryHead(i,_naplWaterHead);
//            _S_i_of_h=mainImbibition_water3p.getSbar();
//            _water_dSatCurve_dHead=mainImbibition_water3p.getDsBar_DpC();
	  
//  	  _waterAppSat = saturationFromImbibitionScanningCurve();
//            _dWaterAppSat_dSatCurve = dAppSat_dS_i(); 
//  	}
//        else if (_naplWaterHead > h.lastNaplWaterHead()) //still draining
//  	{
//  	  h.lastDtoIReversalPoint_water(_S_di,_h_di);
//            h.lastItoDReversalPoint_water(_S_id,_h_id);
//  	  if ( _naplWaterHead > _h_di ) //get on previous d scanning curve
//              {
//                h.nextToLastDtoIReversalPoint_water(_S_di,_h_di);
//                h.nextToLastItoDReversalPoint_water(_S_id,_h_id);
//              }

	  
//  	  mainDrainage_water3p.setCapillaryHead(i,_h_di);
//            _S_d_of_h_di=mainDrainage_water3p.getSbar();
//  	  mainDrainage_water3p.setCapillaryHead(i,_h_id);
//            _S_d_of_h_id=mainDrainage_water3p.getSbar();
//  	  mainDrainage_water3p.setCapillaryHead(i,_naplWaterHead);
//            _S_d_of_h=mainDrainage_water3p.getSbar();
//  	  _water_dSatCurve_dHead=mainDrainage_water3p.getDsBar_DpC();

//  	  _waterAppSat = saturationFromDrainageScanningCurve();
//            _dWaterAppSat_dSatCurve = dAppSat_dS_d(); 
//  	}
//        else //naplWaterHead == _lastNaplWaterHead
//  	{
//  	  mainDrainage_water3p.setCapillaryHead(i,_naplWaterHead);
//            _S_d_of_h=mainDrainage_water3p.getSbar();
//  	  _water_dSatCurve_dHead=mainDrainage_water3p.getDsBar_DpC();
//            _waterAppSat = h.lastWaterAppSat();          
//            _dWaterAppSat_dSatCurve = dAppSat_dS_d(); 
//          }
//      }
}

void ParkerHysteresis::compute2pWaterAppSat()
{		
  if (_airWaterHead <= 0.0 )
    {
      _water_dSatCurve_dHead = 0.0;
      _waterAppSat = 1.0;
      _dWaterAppSat_dSatCurve = 1.0;
    }
  else 
    {
      h.lastItoDReversalPoint_water(_S_id,_h_id,_S_d_of_h_id,_S_i_of_h_id);
      h.lastDtoIReversalPoint_water(_S_di,_h_di,_S_d_of_h_di,_S_i_of_h_di); 

      if (h.waterSaturationPath() == PskHistory::IMBIBITION) 
        { 
          if (_airWaterHead > h.lastAirWaterHead() 
              && (h.lastWaterAppSat() - _S_di) > SQRT_MACHINE_EPSILON) //reversal: i to d 
            { 
//                  std::cout<<"reversal i to d "<<i<<std::endl;
              
              _h_id = h.lastAirWaterHead(); 
              _S_id = h.lastWaterAppSat(); 
  
              
              mainDrainage_water2p.setCapillaryHead(i,_h_id);
              _S_d_of_h_id=mainDrainage_water2p.getSbar();

              mainDrainage_water2p.setCapillaryHead(i,_airWaterHead);
              _S_d_of_h=mainDrainage_water2p.getSbar(); 
              _water_dSatCurve_dHead=mainDrainage_water2p.getDsBar_DpC();
              
              _waterAppSat = saturationFromDrainageScanningCurve(); 
              _dWaterAppSat_dSatCurve = dAppSat_dS_d(); 
            } 
          else
            { 
//                std::cout<<"remaining on i "<<i<<std::endl;
              if(_airWaterHead < _h_id) //get on previous i scanning curve 
                {
//                    std::cout<<"getting on previous scanning curve "<<i<<std::endl;
                  int level=1;
                  h.getPreviousItoDReversalPoint_water(level,_S_id,_h_id,_S_d_of_h_id,_S_i_of_h_id); 
                  h.getPreviousDtoIReversalPoint_water(level,_S_di,_h_di,_S_d_of_h_di,_S_i_of_h_di);
                  while (_airWaterHead < _h_id ) 
                    {
                      level++;
                      h.getPreviousItoDReversalPoint_water(level,_S_id,_h_id,_S_d_of_h_id,_S_i_of_h_id); 
                      h.getPreviousDtoIReversalPoint_water(level,_S_di,_h_di,_S_d_of_h_di,_S_i_of_h_di);
                    }
                }
              
              mainImbibition_water2p.setCapillaryHead(i,_airWaterHead);
              _S_i_of_h=mainImbibition_water2p.getSbar();
              _water_dSatCurve_dHead=mainImbibition_water2p.getDsBar_DpC();
              
              _waterAppSat = saturationFromImbibitionScanningCurve();
              _dWaterAppSat_dSatCurve = dAppSat_dS_i(); 
            } 
        } 
      else //draining 
        { 
          if (_airWaterHead < h.lastAirWaterHead() 
              && (_S_id - h.lastWaterAppSat()) > SQRT_MACHINE_EPSILON ) //reversal: d to i 
            { 
//                std::cout<<"reversal d to i "<<i<<std::endl;
              _h_di = h.lastAirWaterHead(); 
              _S_di = h.lastWaterAppSat(); 
              
              mainImbibition_water2p.setCapillaryHead(i,_h_di);
              _S_i_of_h_di=mainImbibition_water2p.getSbar();

              mainImbibition_water2p.setCapillaryHead(i,_airWaterHead);
              _S_i_of_h=mainImbibition_water2p.getSbar();
              _water_dSatCurve_dHead=mainImbibition_water2p.getDsBar_DpC();
              
              _waterAppSat = saturationFromImbibitionScanningCurve();
              _dWaterAppSat_dSatCurve = dAppSat_dS_i(); 
            }
          else
            { 
//                std::cout<<"remaining on d "<<i<<std::endl;
              if (_airWaterHead > _h_di ) //get on previous d scanning curve 
                {
//                    std::cout<<"getting on previous scanning curve "<<i<<std::endl;
                 int level=1;
                  h.getPreviousItoDReversalPoint_water(level,_S_id,_h_id,_S_d_of_h_id,_S_i_of_h_id); 
                  h.getPreviousDtoIReversalPoint_water(level,_S_di,_h_di,_S_d_of_h_di,_S_i_of_h_di);
                  while (_airWaterHead > _h_di ) 
                    {
                      level++;
                      h.getPreviousItoDReversalPoint_water(level,_S_id,_h_id,_S_d_of_h_id,_S_i_of_h_id); 
                      h.getPreviousDtoIReversalPoint_water(level,_S_di,_h_di,_S_d_of_h_di,_S_i_of_h_di);
                    }
                }
              

              mainDrainage_water2p.setCapillaryHead(i,_airWaterHead);
              _S_d_of_h=mainDrainage_water2p.getSbar();
              _water_dSatCurve_dHead=mainDrainage_water2p.getDsBar_DpC();
              
              _waterAppSat = saturationFromDrainageScanningCurve();
              _dWaterAppSat_dSatCurve = dAppSat_dS_d(); 
            }
        }
    }
}

void ParkerHysteresis::compute2pTrappedAirSat()
{
  if (_waterAppSat <= h.histMin2pWaterAppSat())
    {
      _airTrappedByWaterEffSat = 0.0;
      _dAirTrappedByWaterEffSat_dWaterAppSat = 0.0;
    }
  else
    {
      _airTrappedByWaterEffSat = h.maxResAirTrappedByWaterEffSat() * 
        ( _waterAppSat - h.histMin2pWaterAppSat()) /
        (1.0 - h.histMin2pWaterAppSat() );
      
      _dAirTrappedByWaterEffSat_dWaterAppSat = h.maxResAirTrappedByWaterEffSat()/
        (1.0 - h.histMin2pWaterAppSat());
    }
}

void ParkerHysteresis::compute3pTrappedAirSat()
{
  if ( h.histMinTotalAppSat() > h.histMin2pWaterAppSat())
    {
      _trappedAirEffSat = h.maxResAirTrappedByWaterEffSat() * 
        ( ( h.histMinTotalAppSat() - h.histMin2pWaterAppSat() ) /
          (1.0 - h.histMin2pWaterAppSat() ) ) + 
        h.maxResAirTrappedByNaplEffSat() * 
        ((_totalAppSat - h.histMinTotalAppSat() ) /
         (1.0 - h.histMinTotalAppSat() ) ); 
      _dTrappedAirEffSat_dTotalAppSat = h.maxResAirTrappedByNaplEffSat()/
        (1.0 - h.histMinTotalAppSat() );
      if (_waterAppSat >= h.histMinTotalAppSat() )
        {
          _airTrappedByWaterEffSat = h.maxResAirTrappedByWaterEffSat() 
            * ( ( h.histMinTotalAppSat() - h.histMin2pWaterAppSat() ) /
                (1.0 - h.histMin2pWaterAppSat() ) ) 
            + h.maxResAirTrappedByNaplEffSat() 
            * ((_waterAppSat - h.histMinTotalAppSat() ) / (1.0 - h.histMinTotalAppSat() ) );          
          _dAirTrappedByWaterEffSat_dWaterAppSat = h.maxResAirTrappedByNaplEffSat()/ (1.0 - h.histMinTotalAppSat()); 
        }
      else if (_waterAppSat > h.histMin2pWaterAppSat()) //_histMinTotalAppSat > _waterAppSat > _histMin2pWaterAppSat 
        {
          _airTrappedByWaterEffSat = h.maxResAirTrappedByWaterEffSat() * 
            ( ( _waterAppSat - h.histMin2pWaterAppSat() ) /
              (1.0 - h.histMin2pWaterAppSat() ) );
          _dAirTrappedByWaterEffSat_dWaterAppSat = h.maxResAirTrappedByWaterEffSat()/
            (1.0 - h.histMin2pWaterAppSat() );
        }
      else //_waterAppSat <= _histMin2pWaterAppSat
        {
          _airTrappedByWaterEffSat = 0.0;
          _dAirTrappedByWaterEffSat_dWaterAppSat = 0.0;
        }
    }
  else //_histMinTotalAppSat <= _histMin2pWaterAppSat
    {
      _trappedAirEffSat = h.maxResAirTrappedByNaplEffSat() * 
        ((_totalAppSat - h.histMinTotalAppSat() ) /
         (1.0 - h.histMinTotalAppSat() ) );
      _dTrappedAirEffSat_dTotalAppSat = h.maxResAirTrappedByNaplEffSat()/
        (1.0 - h.histMinTotalAppSat() );
      if (_waterAppSat > h.histMinTotalAppSat() )
        {
          _airTrappedByWaterEffSat = h.maxResAirTrappedByNaplEffSat() 
            * ((_waterAppSat - h.histMinTotalAppSat() ) /
               (1.0 - h.histMinTotalAppSat() ) );
          _dAirTrappedByWaterEffSat_dWaterAppSat = h.maxResAirTrappedByNaplEffSat()/
            (1.0 - h.histMinTotalAppSat() ) ;
        }
      else // _waterAppSat <= h.histMinTotalAppSat
        {
          _airTrappedByWaterEffSat = 0.0;
          _dAirTrappedByWaterEffSat_dWaterAppSat = 0.0;
        }
    }
  _airTrappedByNaplEffSat = _trappedAirEffSat - _airTrappedByWaterEffSat;
  _dAirTrappedByNaplEffSat_dTotalAppSat = _dTrappedAirEffSat_dTotalAppSat;
  _dAirTrappedByNaplEffSat_dWaterAppSat = -_dAirTrappedByWaterEffSat_dWaterAppSat;
}

void ParkerHysteresis::computeTrappedNaplSat()
{
  //first compute the total trapped napl

  if (_waterAppSat > h.histMin3pWaterAppSat())
    {
      _totalTrappedNaplEffSat = h.maxResNaplEffSat() * 
        ((_waterAppSat - h.histMin3pWaterAppSat()) / 
         ( 1.0 - h.histMin3pWaterAppSat()));
      _dTotalTrappedNaplEffSat_dWaterAppSat = h.maxResNaplEffSat()/ 
        ( 1.0 - h.histMin3pWaterAppSat());
    }
  else //_waterAppSat <= h.histMin3pWaterAppSat
    {
      _totalTrappedNaplEffSat = 0.0;
      _dTotalTrappedNaplEffSat_dWaterAppSat = 0.0;
    }
  //breakdown of total trapped napl
  
  if ( _waterAppSat > h.histMinTotalAppSat() && h.histMinTotalAppSat() > h.histMin2pWaterAppSat())
    {
      _airTrappedInNaplByNaplEffSat = (h.maxResAirTrappedByNaplEffSat() * 
                                       h.maxResNaplEffSat() * (_waterAppSat - h.histMinTotalAppSat()))/
        ((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMinTotalAppSat()));
      _dAirTrappedInNaplByNaplEffSat_dWaterAppSat =  (h.maxResAirTrappedByNaplEffSat() * h.maxResNaplEffSat())
        /((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMinTotalAppSat()));
      if (h.histMin3pWaterAppSat() > h.histMin2pWaterAppSat())
        {
          _airTrappedInNaplByWaterEffSat = (h.maxResAirTrappedByWaterEffSat() * h.maxResNaplEffSat() * 
                                            (h.histMinTotalAppSat() - h.histMin3pWaterAppSat()))/
            ((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMin2pWaterAppSat()));
          _dAirTrappedInNaplByWaterEffSat_dWaterAppSat = 0.0; 
        }
      else //h.histMin3pWaterAppSat <= h.histMin2pWaterAppSat
        {
          _airTrappedInNaplByWaterEffSat = (h.maxResAirTrappedByWaterEffSat() * 
                                            h.maxResNaplEffSat() * 
                                            (h.histMinTotalAppSat() - h.histMin2pWaterAppSat()))/
            ((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMin2pWaterAppSat()));
          _dAirTrappedInNaplByWaterEffSat_dWaterAppSat = 0.0; 
        }
    }
  else if (_waterAppSat <= h.histMinTotalAppSat() &&  h.histMinTotalAppSat() > h.histMin2pWaterAppSat())
    {
      _airTrappedInNaplByNaplEffSat = 0.0;
      _dAirTrappedInNaplByNaplEffSat_dWaterAppSat = 0.0;
      if ( h.histMin3pWaterAppSat() > h.histMin2pWaterAppSat())
        {
          _airTrappedInNaplByWaterEffSat = (h.maxResAirTrappedByWaterEffSat() * 
                                            h.maxResNaplEffSat() * (_waterAppSat - h.histMin3pWaterAppSat()))/
            ((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMin2pWaterAppSat()));
          _dAirTrappedInNaplByWaterEffSat_dWaterAppSat = h.maxResAirTrappedByWaterEffSat() * 
            h.maxResNaplEffSat()/((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMin2pWaterAppSat()));
        }
      else //h.histMin3pWaterAppSat <= h.histMin2pWaterAppSat
        {
          _airTrappedInNaplByWaterEffSat = (h.maxResAirTrappedByWaterEffSat() * 
                                            h.maxResNaplEffSat() * (_waterAppSat - h.histMin2pWaterAppSat()))/
            ((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMin2pWaterAppSat()));
          _dAirTrappedInNaplByWaterEffSat_dWaterAppSat = h.maxResAirTrappedByWaterEffSat() * 
            h.maxResNaplEffSat()/((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMin2pWaterAppSat()));
        }
    }
  else if ( _waterAppSat > h.histMinTotalAppSat() &&  h.histMinTotalAppSat() <= h.histMin2pWaterAppSat())
    {
      _airTrappedInNaplByNaplEffSat = (h.maxResAirTrappedByNaplEffSat() * 
                                       h.maxResNaplEffSat() * (_waterAppSat - h.histMinTotalAppSat()))/
        ((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMinTotalAppSat()));
      _dAirTrappedInNaplByNaplEffSat_dWaterAppSat = h.maxResAirTrappedByNaplEffSat() * 
        h.maxResNaplEffSat() / ((1.0 - h.histMin3pWaterAppSat())*(1.0 - h.histMinTotalAppSat()));
      _airTrappedInNaplByWaterEffSat = 0.0;
      _dAirTrappedInNaplByWaterEffSat_dWaterAppSat = 0.0;
    }
  else // _waterAppSat <= h.histMinTotalAppSat AND h.histMinTotalAppSat <=  h.histMin2pWaterAppSat
    {
      _airTrappedInNaplByNaplEffSat = 0.0;
      _dAirTrappedInNaplByNaplEffSat_dWaterAppSat = 0.0;
      _airTrappedInNaplByWaterEffSat = 0.0;
      _dAirTrappedInNaplByWaterEffSat_dWaterAppSat = 0.0;
    }
  _trappedNaplEffSat = _totalTrappedNaplEffSat - _airTrappedInNaplByWaterEffSat - _airTrappedInNaplByNaplEffSat;
  _dTrappedNaplEffSat_dWaterAppSat =  _dTotalTrappedNaplEffSat_dWaterAppSat - _dAirTrappedInNaplByWaterEffSat_dWaterAppSat 
    - _dAirTrappedInNaplByNaplEffSat_dWaterAppSat;
}

void ParkerHysteresis::setHeads(int node,real psiWIn, real psiNIn)
{
  //this is for twophase flow so psiN is considered the air head and we set naplH=psiW so that total sat == water sat
  i=node;
//    MualemVanGenuchten2p::setHeads(i,psiWIn);
//    std::cout<<'\t'<<thetaW<<'\t'<<DthetaW_DpC<<'\t'<<krW<<'\t'<<krN<<std::endl;
//    std::cout<<'\t'<<sBar<<'\t'<<DsBar_DpC<<std::endl;
//    mainImbibition_water2p.setCapillaryHead(i,-psiWIn);
//    std::cout<<'\t'<<mainImbibition_water2p.getSbar()<<'\t'<<mainImbibition_water2p.getDsBar_DpC()<<std::endl;
//    mainDrainage_water2p.setCapillaryHead(i,-psiWIn);
//    std::cout<<'\t'<<mainDrainage_water2p.getSbar()<<'\t'<<mainDrainage_water2p.getDsBar_DpC()<<std::endl;
  h.node(i);
  real waterSat;//,airSat;
  _airWaterHead = psiNIn - psiWIn;
  h.airWaterHead() = psiNIn - psiWIn;

  compute2pWaterAppSat();
//    mainDrainage_water2p.setCapillaryHead(i,_airWaterHead);
//    std::cout<<_waterAppSat<<'\t'<<mainDrainage_water2p.getSbar()<<'\t'<<_water_dSatCurve_dHead<<'\t'<<mainDrainage_water2p.getDsBar_DpC()<<std::endl;
//    _waterAppSat=mainDrainage_water2p.getSbar(); 
//    _dWaterAppSat_dSatCurve=1.0;
//    _water_dSatCurve_dHead=mainDrainage_water2p.getDsBar_DpC();
  
  //now we have apparent saturations and we know whether it's a two phase or three phase case
  
  compute2pTrappedAirSat();

  _waterEffSat = _waterAppSat - _airTrappedByWaterEffSat;
  _airEffSat = 1.0 - _waterEffSat;
//    _dAirTrappedByWaterEffSat_dWaterAppSat=0.0;
//    _waterEffSat = _waterAppSat;

  waterSat = _waterEffSat * (1.0 - h.irreducibleWaterSat()) + h.irreducibleWaterSat();
  //airSat = 1.0 - waterSat;
  _dWaterSat_dAirHead   = (1.0 - h.irreducibleWaterSat())* (1.0 - _dAirTrappedByWaterEffSat_dWaterAppSat)
    *_dWaterAppSat_dSatCurve*_water_dSatCurve_dHead;
  _dWaterSat_dWaterHead = -_dWaterSat_dAirHead;
  _dWaterSat_dNaplHead  = 0.0;
  _dAirSat_dAirHead    = - _dWaterSat_dAirHead;
  _dAirSat_dWaterHead  = _dWaterSat_dAirHead;
  _dAirSat_dNaplHead   =0.0;
  
  //set all the values that are necessary into the history
  h.waterAppSat() = _waterAppSat;

  //get in terms of theta
  thetaW = waterSat * thetaS[i];
  DthetaW_DpC = _dWaterSat_dAirHead*thetaS[i];
  
  real corrFac;
  if (h.maxResAirTrappedByWaterEffSat() == 0.0 || h.histMin2pWaterAppSat() == 0.0)
    corrFac = 0.0;
  else
    corrFac = h.maxResAirTrappedByWaterEffSat() / ( 1.0 - h.histMin2pWaterAppSat());
  
  krW = sqrt(_waterEffSat)*
    pow(1.0 
        -
        (1.0 - corrFac)*
        pow(1.0 - pow(_waterAppSat,1.0/m[i]),m[i]) 
        -
        corrFac*
        pow(1.0 - pow(h.histMin2pWaterAppSat(),1.0/m[i]),m[i]),2);
  
  KW = KWs[i]*krW;
  
  krN = sqrt(_airEffSat - _airTrappedByWaterEffSat)*
    pow(1.0 - pow(_waterAppSat,1.0/m[i]),2*m[i]);
  
  KN = KWs[i]*muW_by_muN*krN;
  //    std::cout<<i<<'\t'<<thetaW<<'\t'<<DthetaW_DpC<<'\t'<<krW<<'\t'<<krN<<std::endl;
}

void ParkerHysteresis::computeSaturationsFromHeads(const real& waterH,const real& naplH,
					       const real& airH, real& waterSat,real& naplSat, 
					       real& airSat)
{		
  _airWaterHead = airH - waterH;
  _airNaplHead = airH - naplH;
  _naplWaterHead = naplH - waterH;
  h.airWaterHead() = airH - waterH;
  h.airNaplHead() = airH - naplH;
  h.naplWaterHead() = naplH - waterH;
  
  //compute apparent saturations for the three phase case
  
  computeTotalAppSat();
  compute3pWaterAppSat();
  
  
  // check to see if the system is still air water (two phase) and act accordingly 
  if ( h.lastState() == PskHistory::TWOPHASE && ((_totalAppSat - _waterAppSat) < _naplTolerance) )
    {
      compute2pWaterAppSat(); //still twophase; get water saturation from two phase relation
    }
  else
    {
      cout<<_waterAppSat<<'\t'<<_naplWaterHead<<endl;
      cout<<"3p"<<_totalAppSat<<'\t'<<_airNaplHead<<endl;
      h.state() = PskHistory::THREEPHASE;
    }
  //now we have apparent saturations and we know whether it's a two phase or three phase case
  
  if ( h.state() == PskHistory::TWOPHASE)
    {
      compute2pTrappedAirSat();
      _waterEffSat = _waterAppSat - _airTrappedByWaterEffSat;
      waterSat = _waterEffSat * (1.0 - h.irreducibleWaterSat()) + h.irreducibleWaterSat();
      naplSat = 0.0;
      airSat = 1.0 - waterSat;
      _dWaterSat_dAirHead   = (1.0 - h.irreducibleWaterSat())* (1.0 - _dAirTrappedByWaterEffSat_dWaterAppSat)
        *_dWaterAppSat_dSatCurve*_water_dSatCurve_dHead;
      _dWaterSat_dWaterHead = -_dWaterSat_dAirHead;
      _dWaterSat_dNaplHead  = 0.0;
      _dNaplSat_dAirHead    =0.0;
      _dNaplSat_dWaterHead  =0.0;
      _dNaplSat_dNaplHead   =0.0;
      _dAirSat_dAirHead    = - _dWaterSat_dAirHead;
      _dAirSat_dWaterHead  = _dWaterSat_dAirHead;
      _dAirSat_dNaplHead   =0.0;
    }
  else //three phase system
    {
      compute3pTrappedAirSat();
      computeTrappedNaplSat();

      //compute effective saturations
      _waterEffSat = _waterAppSat - _trappedNaplEffSat - _airTrappedByWaterEffSat;
      _naplEffSat = _totalAppSat - _waterEffSat - _trappedAirEffSat;
      
      //compute actual saturations
      waterSat = _waterEffSat * (1.0 - h.irreducibleWaterSat()) + h.irreducibleWaterSat();
      naplSat = _naplEffSat * (1.0 - h.irreducibleWaterSat());
      airSat = 1.0 - waterSat - naplSat;

      _dWaterSat_dAirHead   = 0.0;

      _dWaterSat_dNaplHead  = (1.0 - _dTrappedNaplEffSat_dWaterAppSat - _dAirTrappedByWaterEffSat_dWaterAppSat)*
        _dWaterAppSat_dSatCurve*_water_dSatCurve_dHead* (1.0 - h.irreducibleWaterSat());
      
      _dWaterSat_dWaterHead = -_dWaterSat_dNaplHead;
      
      _dNaplSat_dAirHead    =((1.0 - _dTrappedAirEffSat_dTotalAppSat) *_dTotalAppSat_dSatCurve*_total_dSatCurve_dHead 
                             - (1.0 - _dTrappedNaplEffSat_dWaterAppSat - _dAirTrappedByWaterEffSat_dWaterAppSat)*
                             _dWaterAppSat_dSatCurve*_water_dSatCurve_dHead)
        *(1.0 - h.irreducibleWaterSat());
      
      _dNaplSat_dWaterHead  = ((1.0 - _dTrappedNaplEffSat_dWaterAppSat - _dAirTrappedByWaterEffSat_dWaterAppSat)*
                              _dWaterAppSat_dSatCurve*_water_dSatCurve_dHead)
        *(1.0 - h.irreducibleWaterSat());
      
      _dNaplSat_dNaplHead   = -(1.0 - _dTrappedAirEffSat_dTotalAppSat) *_dTotalAppSat_dSatCurve*_total_dSatCurve_dHead
        *(1.0 - h.irreducibleWaterSat());

      _dAirSat_dAirHead    = -_dWaterSat_dAirHead - _dNaplSat_dAirHead;
      _dAirSat_dWaterHead  = -_dWaterSat_dWaterHead - _dNaplSat_dWaterHead;
      _dAirSat_dNaplHead   = -_dWaterSat_dNaplHead - _dNaplSat_dNaplHead;
    }
  
  //set all the values that are necessary into the history
  h.waterAppSat() = _waterAppSat;
  h.totalAppSat() = _totalAppSat;
}

void ParkerHysteresis::updateHistory()
{
  std::cout<<"******************update***************"<<std::endl;
  for(int i=0;i<m.dim();i++)
    {
      h.node(i);
      //update reversal points
  
      h.lastState() = h.state();
      real S,H;
      
      if (h.state() == PskHistory::TWOPHASE)
        {
          if (h.airWaterHead() > 0)
            {
              h.lastDtoIReversalPoint_water(_S_di,_h_di,_S_d_of_h_id,_S_i_of_h_id);
              h.lastItoDReversalPoint_water(_S_id,_h_id,_S_d_of_h_di,_S_i_of_h_di);
              if ( h.waterSaturationPath() == PskHistory::IMBIBITION)
                {
                  if ( h.airWaterHead() > h.lastAirWaterHead() 
                       && (h.lastWaterAppSat() - _S_di) > SQRT_MACHINE_EPSILON)
                    {
                      std::cout<<"checking in reversal point i to d "<<i<<std::endl;
                      h.waterSaturationPath() = PskHistory::DRAINAGE;

                      mainDrainage_water2p.setCapillaryHead(i,h.lastAirWaterHead());
                      _S_d_of_h_id=mainDrainage_water2p.getSbar();

                      mainImbibition_water2p.setCapillaryHead(i,h.lastAirWaterHead());
                      _S_i_of_h_id=mainImbibition_water2p.getSbar();
                      
                      h.pushItoDReversalPoint_water(h.lastWaterAppSat(),
                                                    h.lastAirWaterHead(),
                                                    _S_d_of_h_id,
                                                    _S_i_of_h_id);

                      h.lastDtoIReversalPoint_water(S,H,_S_d_of_h_di,_S_i_of_h_di);
                      while (h.waterAppSat() < S)
                        {
                          std::cout<<"popping reversal points on drain scan after reversal"<<i<<std::endl;
                          h.popDtoIReversalPoint_water();
                          h.popItoDReversalPoint_water();
                          h.lastDtoIReversalPoint_water(S,H,_S_d_of_h_di,_S_i_of_h_di);
                        }          
                    }
                  else
                    {
                      h.lastItoDReversalPoint_water(S,H,_S_d_of_h_id,_S_i_of_h_id);
                      while (h.waterAppSat() > S)
                        {
                          std::cout<<"popping reversal points on imb scan "<<i<<std::endl;
                          h.popItoDReversalPoint_water();
                          h.popDtoIReversalPoint_water();
                          h.lastItoDReversalPoint_water(S,H,_S_d_of_h_id,_S_i_of_h_id);
                        }          
                    }
                }
              else
                {
                  if ( h.airWaterHead() < h.lastAirWaterHead() 
                       && (_S_id - h.lastWaterAppSat()) > SQRT_MACHINE_EPSILON)
                    {
                      std::cout<<"checking in reversal point d to i "<<i<<std::endl;
                      h.waterSaturationPath() = PskHistory::IMBIBITION;

                      mainDrainage_water2p.setCapillaryHead(i,h.lastAirWaterHead());
                      _S_d_of_h_di=mainDrainage_water2p.getSbar();

                      mainImbibition_water2p.setCapillaryHead(i,h.lastAirWaterHead());
                      _S_i_of_h_di=mainImbibition_water2p.getSbar();

                      h.pushDtoIReversalPoint_water(h.lastWaterAppSat(),
                                                    h.lastAirWaterHead(),
                                                    _S_d_of_h_di,
                                                    _S_i_of_h_di);

                      h.lastItoDReversalPoint_water(S,H,_S_d_of_h_id,_S_i_of_h_id);
                      while (h.waterAppSat() > S)
                        {
                          std::cout<<"popping reversal points on imb scan after reversal"<<i<<std::endl;
                          h.popItoDReversalPoint_water();
                          h.popDtoIReversalPoint_water();
                          h.lastItoDReversalPoint_water(S,H,_S_d_of_h_id,_S_i_of_h_id);
                        }          
                    }
                  else 
                    {
                      h.lastDtoIReversalPoint_water(S,H,_S_d_of_h_id,_S_i_of_h_id);
                      while (h.waterAppSat() < S)
                        {
                          std::cout<<"popping reversal points on drain scan "<<i<<std::endl;
                          h.popDtoIReversalPoint_water();
                          h.popItoDReversalPoint_water();
                          h.lastDtoIReversalPoint_water(S,H,_S_d_of_h_id,_S_i_of_h_id);
                        }          
                    }
                }
            }
        }
      else
        {
//            if ( h.waterSaturationPath() == PskHistory::IMBIBITION && 
//                 (h.waterAppSat() <= h.lastWaterAppSat() && h.airWaterHead() > h.lastAirWaterHead()) )
//              {
//                h.waterSaturationPath() = PskHistory::DRAINAGE;
//                h.pushItoDReversalPoint_water(h.lastWaterAppSat(),
//                                              h.lastNaplWaterHead());
//              }
//            else if ( h.waterSaturationPath() == PskHistory::IMBIBITION && 
//                      (h.waterAppSat() >= h.lastWaterAppSat() && h.airWaterHead() < h.lastAirWaterHead()) )
//              {
//                h.lastItoDReversalPoint_water(S,H);
//                if (h.waterAppSat() > S)
//                  {
//                    h.popDtoIReversalPoint_water();
//                    h.popItoDReversalPoint_water();
//                  }          
//              }
//            else if (h.waterSaturationPath() == PskHistory::DRAINAGE && 
//                     (h.waterAppSat() >= h.lastWaterAppSat() && h.airWaterHead() < h.lastAirWaterHead()))
//              {
//                h.waterSaturationPath() = PskHistory::IMBIBITION;
//                h.pushDtoIReversalPoint_water(h.lastWaterAppSat(),
//                                              h.lastNaplWaterHead());
//              }
//            else if (h.waterSaturationPath() == PskHistory::DRAINAGE && 
//                     (h.waterAppSat() <= h.lastWaterAppSat() && h.airWaterHead() > h.lastAirWaterHead()))
//              {
//                h.lastDtoIReversalPoint_water(S,H);
//                if (h.waterAppSat() < S)
//                  {
//                    h.popDtoIReversalPoint_water();
//                    h.popItoDReversalPoint_water();
//                  }          
//              }
//            if ( h.totalSaturationPath() == PskHistory::IMBIBITION && 
//                 h.totalAppSat() < h.lastTotalAppSat())
//              {
//                h.totalSaturationPath() = PskHistory::DRAINAGE;
//                h.pushItoDReversalPoint_total(h.lastTotalAppSat(),
//                                              h.lastAirNaplHead());
//              }
//            else if ( h.totalSaturationPath() == PskHistory::IMBIBITION && 
//                      h.totalAppSat() > h.lastTotalAppSat())
//              {
//                h.lastItoDReversalPoint_total(S,H);
//                if (h.totalAppSat() > S)
//                  {
//                    h.popDtoIReversalPoint_total();
//                    h.popItoDReversalPoint_total();
//                  }          
//              }
//            else if (h.totalSaturationPath() == PskHistory::DRAINAGE && 
//                     h.totalAppSat() > h.lastTotalAppSat())
//              {
//                h.totalSaturationPath() = PskHistory::IMBIBITION;
//                h.pushDtoIReversalPoint_total(h.lastTotalAppSat(),
//                                              h.lastAirNaplHead());
//              }
//            else if (h.totalSaturationPath() == PskHistory::DRAINAGE && 
//                     h.totalAppSat() < h.lastTotalAppSat())
//              {
//                h.lastDtoIReversalPoint_total(S,H);
//                if (h.totalAppSat() < S)
//                  {
//                    h.popDtoIReversalPoint_total();
//                    h.popItoDReversalPoint_total();
//                  }          
//              }
        }
      
      //update last saturations and heads
      h.lastWaterAppSat() = h.waterAppSat();
      h.lastTotalAppSat() = h.totalAppSat();
      h.lastAirWaterHead() = h.airWaterHead();
      h.lastAirNaplHead() = h.airNaplHead();
      h.lastNaplWaterHead() = h.naplWaterHead();
      
      //update historic minimums and maximum residuals if necessary
      
      if (h.state() == PskHistory::TWOPHASE)
        {
          if (h.waterAppSat() < h.histMin2pWaterAppSat())
            {
              h.histMin2pWaterAppSat() = h.waterAppSat();
              h.maxResAirTrappedByWaterEffSat() = (1.0 - h.histMin2pWaterAppSat()) / 
                ( 1.0 + (1.0/h.maxResAirOnMainWaterImbibition() - 1.0) * 
                  (1.0 - h.histMin2pWaterAppSat()));
              h.histMin3pWaterAppSat() = h.histMin2pWaterAppSat();
              h.histMinTotalAppSat() = h.histMin2pWaterAppSat();
              h.maxResAirTrappedByNaplEffSat() = ((1.0 - h.histMinTotalAppSat()) / 
                                                  (1.0 + 
                                                   (1.0/h.maxResAirOnMainNaplImbibition()  
                                                    - 1.0) *
                                                   (1.0 - h.histMinTotalAppSat())));
            }
        }
      else
        {
          if (h.waterAppSat() < h.histMin3pWaterAppSat())
            {
              h.histMin3pWaterAppSat() = h.waterAppSat();
              h.maxResAirTrappedByWaterEffSat() = (1.0 - h.histMin3pWaterAppSat()) / 
                ( 1.0 + (1.0/h.maxResAirOnMainWaterImbibition() - 1.0) * 
                  (1.0 - h.histMin3pWaterAppSat()));
              h.maxResNaplEffSat() = (1.0 - h.histMin3pWaterAppSat()) / 
                ( 1.0 + (1.0/h.maxResNaplOnMainWaterImbibition() - 1.0) * 
                  (1.0 - h.histMin3pWaterAppSat()));
            }		
          if (h.totalAppSat() < h.histMinTotalAppSat())
            {
              h.histMinTotalAppSat() = h.totalAppSat();
              h.maxResAirTrappedByNaplEffSat() = ((1.0 - h.histMinTotalAppSat()) / 
                                                  (1.0 + 
                                                   (1.0/h.maxResAirOnMainNaplImbibition()  
                                                    - 1.0) *
                                                   (1.0 - h.histMinTotalAppSat())));
            }
        }
    }
}
}//Daetk









