/*******************************************************************************
 * Ech2o, a spatially-distributed, ecohydrologic simulator
 * Copyright (c) 2016 Marco Maneta <marco.maneta@umontana.edu>
 *
 *     This file is part of ech2o, a hydrologic model developed at the 
 *     University of Montana.
 *
 *     Ech2o is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     Ech2o is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Ech2o.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *    Marco Maneta, Sylvain Kuppel
 *******************************************************************************/
/*
 * MixingV_plusExtraGW.cpp
 *
 *  Created on: May 14, 2020
 *      Author: Xiaoqiang Yang 
 */

#include "Basin.h"

void Tracking::MixingV_plusExtraGW(Basin &bsn, Control &ctrl, 
			     double &d1, double &d2, double &d3, double &fc,
			     double &Qk1, double &dtdx, double &dx, double &dt, int r, int c)
{
  int mixmod = ctrl.toggle_mix;

  // Soil state before routing
  double pond_old = bsn.getPondingWater()->matrix[r][c];//without including exfiltration, so it's "old"
  double theta1_old = bsn.getSoilMoist1()->matrix[r][c];
  double theta2_old = bsn.getSoilMoist2()->matrix[r][c];
  double theta3_old = bsn.getSoilMoist3()->matrix[r][c];

  // Vertical fluxes
  double L1toSrf = bsn.getFluxExfilt()->matrix[r][c];
  double L2toL1 = bsn.getFluxL2toL1()->matrix[r][c];
  double L3toL2 = bsn.getFluxL3toL2()->matrix[r][c];
  double GWtoChn = bsn.getFluxGWtoChn()->matrix[r][c];
  double SnowtoSrf = bsn.getFluxSnowtoSrf()->matrix[r][c];
  double SrftoChn = bsn.getFluxSrftoChn()->matrix[r][c];														
  // Lateral out
  double GWtoLat = bsn.getFluxGWtoLat()->matrix[r][c];
  double ChntoLat = Qk1*dtdx/dx;
  double SrftoLat = bsn.getFluxSrftoLat()->matrix[r][c];
  // Lateral in
  double LattoGW = bsn.getFluxLattoGW()->matrix[r][c]; 
  double LattoChn = bsn.getFluxLattoChn()->matrix[r][c]; 
  double LattoSrf = bsn.getFluxLattoSrf()->matrix[r][c]; 

  // Extra GW yangx 2020-05
  double ExtraGWtoChn = bsn.getFluxExtraGWtoChn()->matrix[r][c];
  double ExtraGWtoLat = bsn.getFluxExtraGWtoLat()->matrix[r][c];
  double LattoExtraGW = bsn.getFluxLattoExtraGW()->matrix[r][c];
  double ExtraGWtoLat_old = bsn.getFluxExtraGWtoLat_old()->matrix[r][c];
  double ExtraGWtoChn_old = bsn.getFluxExtraGWtoChn_old()->matrix[r][c];
  double LattoExtraGW_old = bsn.getFluxLattoExtraGW_old()->matrix[r][c];
  double ExtraGW = bsn.getExtraGW()->matrix[r][c];
  double Leakage = bsn.getBedrockLeakage()->matrix[r][c] * dt; 
  double Leakage_old = bsn.getBedrockLeakage_old()->matrix[r][c] * dt; 
  //leak is stored as m/s in ln 155 SoilWaterRedistribution.cpp

  // For GW and surface (pond+channel), equivalent lateral inputs values
  double FinSrf = SrftoChn + GWtoChn + LattoChn + ExtraGWtoChn; //Extra GW --yangx
  double FinSrf2 = L1toSrf ; //sum of ponding inputs that needed to be mixed, the rest has already handled in mixingV_down.cpp
  double channel_old = bsn.getChannel_old()->matrix[r][c];//getChannelWater()																		   
  double d2Hin = 0;
  double d18Oin = 0;
  double Agein= 0;
  double d2HinChn = 0;
  double d18OinChn = 0;
  double AgeinChn= 0; 
  double ex_d2Hin = 0;
  double ex_d2Hin_old = 0;
  double ex_d18Oin = 0;
  double ex_d18Oin_old = 0;
  double ex_Agein = 0;
  double ex_Agein_old = 0;  
  
  double exq_out = 0;
  double exq_out_old = 0;

  // Two-pore stuff
  double theta_MW1 = 0;
  double theta_MW2 = 0;
  double theta_r = 0;
  //double porosity = 0;
  double L3toTB2 = 0;
  double L3toMW2 = 0;
  double MW2toTB1 = 0;
  double MW2toMW1 = 0;
  double MW1toSrf = 0;

  double d_old = 0;

  //initialize _old variables if "current_ts_count = 1"
  if(ctrl.current_ts_count == 1){
	LattoExtraGW_old = LattoExtraGW;
	ExtraGWtoLat_old = ExtraGWtoLat;
	ExtraGWtoChn_old = ExtraGWtoChn;
	//_Fd2HLattoExtraGW_old->matrix[r][c] = _Fd2HLattoExtraGW->matrix[r][c];
  }


  if(ctrl.sw_TPD){
    theta_r = bsn.getSoilMoistR()->matrix[r][c];
    theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
    theta_MW2 = bsn.getMoistureMW2()->matrix[r][c];
    // Return flow to L2: weighted between TB2 (if there's deficit there) and MW2
    L3toTB2 = std::min<double>(L3toL2,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d2*(theta_MW2-theta2_old)));
    L3toMW2 = std::max<double>(0, L3toL2 - L3toTB2);
    // Return flow to L1: : weighted between TB1 (if there's deficit there) and MW1
    MW2toTB1 = std::min<double>(L2toL1,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d1*(theta_MW1-theta1_old)));
    MW2toMW1 = std::max<double>(0,L2toL1 - MW2toTB1);
    // Return flow to surface: : only from mobile water in L1
    MW1toSrf = L1toSrf;
  }

  // Extra GW --yangx 2020-05
  // update tracer signatures in deepGW flux 
  if(LattoExtraGW > RNDOFFERR){

	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){
		exq_out = ExtraGWtoLat + ExtraGWtoChn;
        exq_out_old = ExtraGWtoLat_old + ExtraGWtoChn_old;
    } else {
		exq_out = ExtraGWtoLat;
        exq_out_old = ExtraGWtoLat_old;
    }	  
    if(ctrl.sw_2H){
	//up-slope grid input, might be receving more than one grid cell
	  ex_d2Hin = _Fd2HLattoExtraGW->matrix[r][c] / LattoExtraGW;
	  ex_d2Hin_old = _Fd2HLattoExtraGW_old->matrix[r][c] / LattoExtraGW_old;
	  _d2HExtraGWtoLat->matrix[r][c]  = DeepGWMixing(Leakage, _d2Hleakage->matrix[r][c],
	                     LattoExtraGW, ex_d2Hin, exq_out,
						 Leakage_old, _d2Hleakage_old->matrix[r][c],
	                     LattoExtraGW_old, ex_d2Hin_old, exq_out_old,
						 _d2HExtraGWtoLat->matrix[r][c], ExtraGW);
    }
    if(ctrl.sw_18O){
	  //up-slope grid input, might be receving more than one grid cell
	  ex_d18Oin = _Fd18OLattoExtraGW->matrix[r][c] / LattoExtraGW;
	  ex_d18Oin_old = _Fd18OLattoExtraGW_old->matrix[r][c] / LattoExtraGW_old;
	  _d18OExtraGWtoLat->matrix[r][c]  = DeepGWMixing(Leakage, _d18Oleakage->matrix[r][c],
	                     LattoExtraGW, ex_d18Oin, exq_out,
						 Leakage_old, _d18Oleakage_old->matrix[r][c],
	                     LattoExtraGW_old, ex_d18Oin_old, exq_out_old,
						 _d18OExtraGWtoLat->matrix[r][c], ExtraGW);
    }  
    if(ctrl.sw_Age){
	//up-slope grid input, might be receving more than one grid cell
	  ex_Agein = _FAgeLattoExtraGW->matrix[r][c] / LattoExtraGW;
	  ex_Agein_old = _FAgeLattoExtraGW_old->matrix[r][c] / LattoExtraGW_old;
	  _AgeExtraGWtoLat->matrix[r][c]  = DeepGWMixing(Leakage, _Ageleakage->matrix[r][c],
	                     LattoExtraGW, ex_Agein, exq_out,
						 Leakage_old, _Ageleakage_old->matrix[r][c],
	                     LattoExtraGW_old, ex_Agein_old, exq_out_old,
						 _AgeExtraGWtoLat->matrix[r][c], ExtraGW);
	//
	  if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	    _AgeExtraGWtoChn->matrix[r][c] = ExtraGWtoChn > RNDOFFERR ? _AgeExtraGWtoLat->matrix[r][c] : 0.0;
	  }
    }
  }
  // Layer 3 (GW included) --------------------------------------------------------------------

  if(LattoGW > RNDOFFERR){  
    if(ctrl.sw_2H){
      // Equivalent input signature for GW
      d2Hin = _Fd2HLattoGW->matrix[r][c] / LattoGW ;
      _d2Hsoil3->matrix[r][c] = InOutMix(theta3_old*d3, _d2Hsoil3->matrix[r][c],
					 LattoGW, d2Hin, L3toL2+GWtoLat, mixmod);
      _d2Hgroundwater->matrix[r][c] = _d2Hsoil3->matrix[r][c];
    }
    if(ctrl.sw_18O){
      d18Oin = _Fd18OLattoGW->matrix[r][c] / LattoGW ;
      _d18Osoil3->matrix[r][c] = InOutMix(theta3_old*d3, _d18Osoil3->matrix[r][c],
					  LattoGW, d18Oin, L3toL2+GWtoLat, mixmod);
      _d18Ogroundwater->matrix[r][c] = _d18Osoil3->matrix[r][c];
    }
    if(ctrl.sw_Age){
      Agein = _FAgeLattoGW->matrix[r][c] / LattoGW ;
      _Agesoil3->matrix[r][c] = InOutMix(theta3_old*d3, _Agesoil3->matrix[r][c],
					 LattoGW, Agein, L3toL2+GWtoLat, mixmod);
      _Agegroundwater->matrix[r][c] = _Agesoil3->matrix[r][c];
    }
  }

  // Layer 2 ------------------------------------------------------------------------

  // If two-pore domain activated
  if(ctrl.sw_TPD and L3toL2 > RNDOFFERR){
    // Tightly-bound
    if(ctrl.sw_2H)
      _d2H_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2, 
					_d2H_TB2->matrix[r][c],
					L3toTB2, _d2Hsoil3->matrix[r][c]);
    if(ctrl.sw_18O)
      _d18O_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2, 
					 _d18O_TB2->matrix[r][c],
					 L3toTB2, _d18Osoil3->matrix[r][c]);
    if(ctrl.sw_Age)
      _Age_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2, 
					_Age_TB2->matrix[r][c],
					L3toTB2, _Agesoil3->matrix[r][c]);
    
    // Mobile water
    if(ctrl.sw_2H)
      _d2H_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _d2H_MW2->matrix[r][c], 
		 L3toMW2, _d2Hsoil3->matrix[r][c], L2toL1, mixmod) : _d2Hsoil3->matrix[r][c] ;

    if(ctrl.sw_18O)
      _d18O_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _d18O_MW2->matrix[r][c], 
		 L3toMW2, _d18Osoil3->matrix[r][c], L2toL1, mixmod) : _d18Osoil3->matrix[r][c] ;
    
    if(ctrl.sw_Age)
      _Age_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _Age_MW2->matrix[r][c], 
		 L3toMW2, _Agesoil3->matrix[r][c], L2toL1, mixmod) : _Agesoil3->matrix[r][c] ;

  } else if (L3toL2 > RNDOFFERR) { // Soil-averaged values
    if(ctrl.sw_2H)
      _d2Hsoil2->matrix[r][c] = InOutMix(theta2_old*d2, _d2Hsoil2->matrix[r][c],
					 L3toL2, _d2Hsoil3->matrix[r][c], L2toL1, mixmod);
    if(ctrl.sw_18O)
      _d18Osoil2->matrix[r][c] = InOutMix(theta2_old*d2, _d18Osoil2->matrix[r][c],
					  L3toL2, _d18Osoil3->matrix[r][c], L2toL1, mixmod);
    
    if(ctrl.sw_Age)
      _Agesoil2->matrix[r][c] = InOutMix(theta2_old*d2, _Agesoil2->matrix[r][c],
					 L3toL2, _Agesoil3->matrix[r][c], L2toL1, mixmod);
  }

  // Layer 1 ------------------------------------------------------------------------

  // If two-pore domain activated: return flow only from MW2
  if(ctrl.sw_TPD and L2toL1 > RNDOFFERR){
    // Tightly-bound
    if(ctrl.sw_2H){
      d_old = _d2H_TB1->matrix[r][c];      
      _d2H_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, 
					_d2H_TB1->matrix[r][c],	
					MW2toTB1, _d2H_MW2->matrix[r][c]);

      if(abs(_d2H_TB1->matrix[r][c])>100)
	cout << r << " " << c << "| d2H " << 
	  "| dTB1_new:" << _d2H_TB1->matrix[r][c] << "| dTB1_old:" << d_old << 
	  "| dMW2:" << _d2H_MW2->matrix[r][c] << "| MW2toTB1:" << MW2toTB1 << endl;
    }
    if(ctrl.sw_18O){
      d_old = _d18O_TB1->matrix[r][c];      
      _d18O_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, 
					 _d18O_TB1->matrix[r][c],
					 MW2toTB1, _d18O_MW2->matrix[r][c]);
      
      if(abs(_d18O_TB1->matrix[r][c])>100)
	cout << r << " " << c << "| d18O " << 
	  "| dTB1_new:" << _d18O_TB1->matrix[r][c] << "| dTB1_old:" << d_old << 
	  "| dMW2:" << _d18O_MW2->matrix[r][c] << "| MW2toTB1:" << MW2toTB1 << endl;
    }
    if(ctrl.sw_Age)
      _Age_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, 
					_Age_TB1->matrix[r][c],
					MW2toTB1, _Age_MW2->matrix[r][c]);

    // Mobile water
    if(ctrl.sw_2H){
      d_old = _d2H_MW1->matrix[r][c];
      _d2H_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _d2H_MW1->matrix[r][c], 
		 MW2toMW1, _d2H_MW2->matrix[r][c], L1toSrf, mixmod) : _d2H_MW2->matrix[r][c] ;

      if(abs(_d2H_MW1->matrix[r][c])>100 ) //or (r==80 and c==119))
	cout << r << " " << c << "| d2H " << 
	  "| dMW1_new:" << _d2H_MW1->matrix[r][c] << "| dMW1_old:" << d_old << 
	  "| dMW2:" << _d2H_MW2->matrix[r][c] <<
	  "| MW2toMW1:" << MW2toMW1 << "| MW1toSrf:" << L1toSrf << endl;
    }
    if(ctrl.sw_18O){
      d_old = _d18O_MW1->matrix[r][c];
      _d18O_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _d18O_MW1->matrix[r][c], 
		 MW2toMW1, _d18O_MW2->matrix[r][c], L1toSrf, mixmod) : _d18O_MW2->matrix[r][c] ;

      if(abs(_d18O_MW1->matrix[r][c])>100)// or (r==80 and c==119))
	cout << r << " " << c << "| d2H " << 
	  "| dMW1_new:" << _d18O_MW1->matrix[r][c] << "| dMW1_old:" << d_old << 
	  "| dMW2:" << _d18O_MW2->matrix[r][c] <<
	  "| MW2toMW1:" << MW2toMW1 << "| MW1toSrf:" << L1toSrf << endl;
    }

    if(ctrl.sw_Age)
      _Age_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _Age_MW1->matrix[r][c], 
		 MW2toMW1, _Age_MW2->matrix[r][c], L1toSrf, mixmod) : _Age_MW2->matrix[r][c] ;

  } else if (L2toL1 > RNDOFFERR) { // Soil-averaged values    
    if(ctrl.sw_2H)
      _d2Hsoil1->matrix[r][c] = InOutMix(theta1_old*d1, _d2Hsoil1->matrix[r][c],
					 L2toL1, _d2Hsoil2->matrix[r][c], L1toSrf, mixmod);
    if(ctrl.sw_18O)
      _d18Osoil1->matrix[r][c] = InOutMix(theta1_old*d1, _d18Osoil1->matrix[r][c],
					  L2toL1, _d18Osoil2->matrix[r][c], L1toSrf, mixmod);
    if(ctrl.sw_Age)
      _Agesoil1->matrix[r][c] = InOutMix(theta1_old*d1, _Agesoil1->matrix[r][c],
					 L2toL1, _Agesoil2->matrix[r][c], L1toSrf, mixmod);
  }

  // Surface --------------------------------------------------------------------------------
  
  // If non-channel, return flow + runon + runoff.
  // If channel cell add input discharge, seepage, river outflow.
  // If surface storage initially null, propagated signature is that of inputs
  // Extra GW by yangx 2020-05
  // yangx 2020-11  separate channel and ponding water storage
  if(FinSrf > RNDOFFERR) {
    // If two-pore domain activated: return flow only from MW1
    if(ctrl.sw_TPD){
      if(ctrl.sw_2H){
/* 	d_old = _d2Hsurface->matrix[r][c];
	// Equivalent input signature for surface inputs
	d2Hin = (L1toSrf*_d2H_MW1->matrix[r][c] + GWtoChn *_d2Hgroundwater->matrix[r][c] + 
		 ExtraGWtoChn *_d2HExtraGWtoLat->matrix[r][c] + _Fd2HLattoChn->matrix[r][c]) / FinSrf ;
	_d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  InOutMix(pond_old, _d2Hsurface->matrix[r][c], FinSrf, d2Hin, 
		   ChntoLat+SrftoLat, mixmod): d2Hin; */
	    if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){

		//for channel water, to update inputs: GWtoChn, SrftoChn, and LattoChn
	      d2HinChn = (GWtoChn *_d2Hgroundwater->matrix[r][c] + SrftoChn*_d2Hsurface->matrix[r][c] +
		    ExtraGWtoChn *_d2HExtraGWtoLat->matrix[r][c] + _Fd2HLattoChn->matrix[r][c]) / FinSrf ;
	      _d2Hchan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _d2Hchan->matrix[r][c], FinSrf, d2HinChn, 
		     ChntoLat, mixmod): d2HinChn;          
		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat + SrftoChn  
          d2Hin = _d2H_MW1->matrix[r][c];//d2h for L1toSrf
	      _d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d2Hsurface->matrix[r][c], L1toSrf, d2Hin, 
		     SrftoChn+SrftoLat, mixmod): d2Hin;
			 
        } else {  //or without channel
        //if it's not a channel cell, only have to update 
		//input (exfiltration L1toSrf) d2Hin = _d2Hsoil1->matrix[r][c]
		//output(SrftoLat):
          d2Hin = _d2H_MW1->matrix[r][c];
	      _d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d2Hsurface->matrix[r][c], L1toSrf, d2Hin, 
		     SrftoLat, mixmod): d2Hin;
			 
          _d2Hchan->matrix[r][c] = 0.0;			 
        }
		

	/* if(abs(_d2Hsurface->matrix[r][c])>100 or abs(d_old)>100 )
	//if(bsn.getChannelWidth()->matrix[r][c] > 0)
	cout << r << " " << c << "| d2H-MixingRouting2 " << "| pond_old:" << pond_old << 
	"| L1toSrf:" << L1toSrf << "| GWtoChn:" << GWtoChn << "| LattoChn:" << LattoChn << 
	"| ChntoLat:" << ChntoLat << "| SrftoLat:" << SrftoLat << "| pond_new:" << pond_new <<
	"| dSrf_new:" << _d2Hsurface->matrix[r][c] << "|dSrf_old:" << d_old << 
	"| d2Hchn_lat:" << _FCd2HLattoChn->matrix[r][c]/LattoChn  << 
	"|d2Hgw_tmp:" << d2Hgw_tmp << 
	"|d2Hin:" << d2Hin << "| ischannel(>0):" << bsn.getChannelWidth()->matrix[r][c] << 
	"| mixmod: " << mixmod << endl; */
      }

      if(ctrl.sw_18O){
/* 	d_old = _d18Osurface->matrix[r][c];
	// Equivalent input signature for surface inputs
	d18Oin = (L1toSrf*_d18O_MW1->matrix[r][c] + GWtoChn *_d18Ogroundwater->matrix[r][c] +
		  ExtraGWtoChn *_d18OExtraGWtoLat->matrix[r][c] + _Fd18OLattoChn->matrix[r][c]) / FinSrf ;
	_d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  InOutMix(pond_old, _d18Osurface->matrix[r][c], FinSrf, d18Oin, 
		   ChntoLat+SrftoLat, mixmod): d18Oin; */
	    if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){

		//for channel water, to update inputs: GWtoChn, SrftoChn, and LattoChn
	      d18OinChn = (GWtoChn *_d18Ogroundwater->matrix[r][c] + SrftoChn*_d18Osurface->matrix[r][c] +
		    ExtraGWtoChn *_d18OExtraGWtoLat->matrix[r][c] + _Fd18OLattoChn->matrix[r][c]) / FinSrf ;
	      _d18Ochan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _d18Ochan->matrix[r][c], FinSrf, d18OinChn, 
		     ChntoLat, mixmod): d18OinChn;          

		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat + SrftoChn  
          d18Oin = _d18O_MW1->matrix[r][c];//d2h for L1toSrf
	      _d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d18Osurface->matrix[r][c], L1toSrf, d18Oin, 
		     SrftoChn+SrftoLat, mixmod): d18Oin;
			 
        } else {  //or without channel
        //if it's not a channel cell, only have to update 
		//input (exfiltration L1toSrf)
		//output(SrftoLat):
          d18Oin = _d18O_MW1->matrix[r][c];
	      _d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d18Osurface->matrix[r][c], L1toSrf, d18Oin, 
		     SrftoLat, mixmod): d18Oin;
			 
          _d18Ochan->matrix[r][c] = 0.0;			 
        }

	/* if(abs(_d18Osurface->matrix[r][c])>100 or abs(d_old)>100 )
	//if(bsn.getChannelWidth()->matrix[r][c] > 0)
	cout << r << " " << c << "| d18O-MixingRouting2 " << "| pond_old:" << pond_old << 
	"| L1toSrf:" << L1toSrf << "| GWtoChn:" << GWtoChn << "| LattoChn:" << LattoChn << 
	"| ChntoLat:" << ChntoLat << "| SrftoLat:" << SrftoLat << "| pond_new:" << pond_new <<
	"| dSrf_new:" << _d18Osurface->matrix[r][c] << "|dSrf_old:" << d_old << 
	"| d18Ochn_lat:" << _FCd18OLattoChn->matrix[r][c]/LattoChn  << 
	"|d18Ogw_tmp:" << d18Ogw_tmp << 
	"|d18Oin:" << d18Oin << "| ischannel(>0):" << bsn.getChannelWidth()->matrix[r][c] << 
	"| mixmod: " << mixmod << endl; */
      }

      if(ctrl.sw_Age){
/* 	// Equivalent input signature for surface inputs
	Agein = (L1toSrf*_Age_MW1->matrix[r][c] + GWtoChn *_Agegroundwater->matrix[r][c] +
		 ExtraGWtoChn *_AgeExtraGWtoLat->matrix[r][c] + _FAgeLattoChn->matrix[r][c]) / FinSrf ;
	_Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  InOutMix(pond_old, _Agesurface->matrix[r][c], FinSrf, Agein, 
		   ChntoLat+SrftoLat, mixmod): Agein; */
	    if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){

		//for channel water, to update inputs: GWtoChn, SrftoChn, and LattoChn
	      AgeinChn = (GWtoChn *_Agegroundwater->matrix[r][c] + SrftoChn*_Agesurface->matrix[r][c] +
		    ExtraGWtoChn *_AgeExtraGWtoLat->matrix[r][c] + _FAgeLattoChn->matrix[r][c]) / FinSrf ;
	      _Agechan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _Agechan->matrix[r][c], FinSrf, AgeinChn, 
		     ChntoLat, mixmod): AgeinChn;          
          _AgeGWtoChn->matrix[r][c] = _Agegroundwater->matrix[r][c];
          _AgeSrftoChn->matrix[r][c] = pond_old > RNDOFFERR ? _Agesurface->matrix[r][c] : 0.0; 
		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat + SrftoChn  
          Agein = _Age_MW1->matrix[r][c];//d2h for L1toSrf
	      _Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _Agesurface->matrix[r][c], L1toSrf, Agein, 
		     SrftoChn+SrftoLat, mixmod): Agein;

        } else {  //or without channel
        //if it's not a channel cell, only have to update 
		//input (exfiltration L1toSrf)
		//output(SrftoLat):
          Agein = _Age_MW1->matrix[r][c];
	      _Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _Agesurface->matrix[r][c], L1toSrf, Agein, 
		     SrftoLat, mixmod): Agein;
			 
          _Agechan->matrix[r][c] = 0.0;			 
        }

	// If it's a channel cell, there's no reinfiltration, so that all lateral+snow+exfilt
	// inputs ultimately feed the channel
/* 	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	  _AgeGWtoChn->matrix[r][c] = GWtoChn > RNDOFFERR ? _Agegroundwater->matrix[r][c] : 0.0;
	  _AgeSrftoChn->matrix[r][c] = FinSrf2 > RNDOFFERR ? 
	    (L1toSrf*_Age_MW1->matrix[r][c] + _FAgeLattoSrf->matrix[r][c] + 
	     SnowtoSrf*_Agesnowmelt->matrix[r][c]) / FinSrf2 : 0.0;
	} */
      }
      
    } else {      
      if(ctrl.sw_2H){
	// d2Hin = (L1toSrf*_d2Hsoil1->matrix[r][c] + GWtoChn *_d2Hgroundwater->matrix[r][c] +
		 // ExtraGWtoChn *_d2HExtraGWtoLat->matrix[r][c] + _Fd2HLattoChn->matrix[r][c]) / FinSrf ;
	// _d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  // InOutMix(pond_old, _d2Hsurface->matrix[r][c], FinSrf, d2Hin, 
		   // ChntoLat+SrftoLat, mixmod): d2Hin;
      // }
	    if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){
		//for channel water, to update inputs: GWtoChn, SrftoChn, and LattoChn
	      d2HinChn = (GWtoChn *_d2Hgroundwater->matrix[r][c] + SrftoChn*_d2Hsurface->matrix[r][c] +
		    ExtraGWtoChn *_d2HExtraGWtoLat->matrix[r][c] + _Fd2HLattoChn->matrix[r][c]) / FinSrf ;
	      _d2Hchan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _d2Hchan->matrix[r][c], FinSrf, d2HinChn, 
		     ChntoLat, mixmod): d2HinChn;          
		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat + SrftoChn  
          d2Hin = _d2Hsoil1->matrix[r][c];//d2h for L1toSrf
	      _d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d2Hsurface->matrix[r][c], L1toSrf, d2Hin, 
		     SrftoChn+SrftoLat, mixmod): d2Hin;
        } else {  //or without channel
        //if it's not a channel cell, only have to update 
		//input (exfiltration L1toSrf) d2Hin = _d2Hsoil1->matrix[r][c]
		//output(SrftoLat):
          d2Hin = _d2Hsoil1->matrix[r][c];
	      _d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d2Hsurface->matrix[r][c], L1toSrf, d2Hin, 
		     SrftoLat, mixmod): d2Hin;
			 
          _d2Hchan->matrix[r][c] = 0.0;			 
        }
	  }		
      if(ctrl.sw_18O){
	    if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){
		//for channel water, to update inputs: GWtoChn, SrftoChn, and LattoChn
	      d18OinChn = (GWtoChn *_d18Ogroundwater->matrix[r][c] + SrftoChn*_d18Osurface->matrix[r][c] +
		    ExtraGWtoChn *_d18OExtraGWtoLat->matrix[r][c] + _Fd18OLattoChn->matrix[r][c]) / FinSrf ;
	      _d18Ochan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _d18Ochan->matrix[r][c], FinSrf, d18OinChn, 
		     ChntoLat, mixmod): d18OinChn;          
		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat + SrftoChn  
          d18Oin = _d18Osoil1->matrix[r][c];//d2h for L1toSrf
	      _d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d18Osurface->matrix[r][c], L1toSrf, d18Oin, 
		     SrftoChn+SrftoLat, mixmod): d18Oin;
        } else {  //or without channel
        //if it's not a channel cell, only have to update 
		//input (exfiltration L1toSrf)
		//output(SrftoLat):
          d18Oin = _d18Osoil1->matrix[r][c];
	      _d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d18Osurface->matrix[r][c], L1toSrf, d18Oin, 
		     SrftoLat, mixmod): d18Oin;
			 
          _d18Ochan->matrix[r][c] = 0.0;			 
        }	
	// d18Oin = (L1toSrf*_d18Osoil1->matrix[r][c] + GWtoChn *_d18Ogroundwater->matrix[r][c] +
		  // ExtraGWtoChn *_d18OExtraGWtoLat->matrix[r][c] + _Fd18OLattoChn->matrix[r][c]) / FinSrf ;
	// _d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  // InOutMix(pond_old, _d18Osurface->matrix[r][c], FinSrf, d18Oin, 
		   // ChntoLat+SrftoLat, mixmod): d18Oin;
      // }
      }
      if(ctrl.sw_Age) {
	    if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){

		//for channel water, to update inputs: GWtoChn, SrftoChn, and LattoChn
	      AgeinChn = (GWtoChn *_Agegroundwater->matrix[r][c] + SrftoChn*_Agesurface->matrix[r][c] +
		    ExtraGWtoChn *_AgeExtraGWtoLat->matrix[r][c] + _FAgeLattoChn->matrix[r][c]) / FinSrf ;
	      _Agechan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _Agechan->matrix[r][c], FinSrf, AgeinChn, 
		     ChntoLat, mixmod): AgeinChn;          
          _AgeGWtoChn->matrix[r][c] = _Agegroundwater->matrix[r][c];
          _AgeSrftoChn->matrix[r][c] = pond_old > RNDOFFERR ? _Agesurface->matrix[r][c] : 0.0; 
		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat + SrftoChn  
          Agein = _Agesoil1->matrix[r][c];//d2h for L1toSrf
	      _Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _Agesurface->matrix[r][c], L1toSrf, Agein, 
		     SrftoChn+SrftoLat, mixmod): Agein;
        } else {  //or without channel
        //if it's not a channel cell, only have to update 
		//input (exfiltration L1toSrf)
		//output(SrftoLat):
          Agein = _Agesoil1->matrix[r][c];
	      _Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _Agesurface->matrix[r][c], L1toSrf, Agein, 
		     SrftoLat, mixmod): Agein;
			 
          _Agechan->matrix[r][c] = 0.0;			 
        }
/* 	Agein = (L1toSrf*_Agesoil1->matrix[r][c] + GWtoChn *_Agegroundwater->matrix[r][c] +
		 ExtraGWtoChn *_AgeExtraGWtoLat->matrix[r][c] + _FAgeLattoChn->matrix[r][c]) / FinSrf ;
	_Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  InOutMix(pond_old, _Agesurface->matrix[r][c], FinSrf, Agein, 
		   ChntoLat+SrftoLat, mixmod): Agein;

	// If it's a channel cell, there's no reinfiltration, so that all lateral+snow+exfilt
	// inputs ultimately feed the channel
	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	  _AgeGWtoChn->matrix[r][c] = GWtoChn > RNDOFFERR ? _Agegroundwater->matrix[r][c] : 0.0;
	  _AgeSrftoChn->matrix[r][c] = FinSrf2 > RNDOFFERR ? 
	    (L1toSrf*_Agesoil1->matrix[r][c] + _FAgeLattoSrf->matrix[r][c] + 
	     SnowtoSrf*_Agesnowmelt->matrix[r][c]) / FinSrf2 : 0.0;
	} */
      }
    }
  }
  // -------------------------------------------------------------------------------
}


