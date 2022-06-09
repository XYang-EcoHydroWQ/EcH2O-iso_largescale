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
 * GWrouting.cpp
 *
 *  Created on: Dec 2, 2010
 *      Author: Marco.Maneta
 */

#include"Basin.h"

int Basin::DailyGWRouting(Atmosphere &atm, Control &ctrl, Tracking &trck) {

	int r, c, d;
	int rr, cc;
	bool lat_ok;
	REAL8 dtdx;
	REAL8 alpha;
	REAL8 qj1i;
	REAL8 hji1;
	REAL8 hj1i1;
	REAL8 R;
	REAL8 dt = ctrl.dt;
	REAL8 poros1, poros2, poros3; //porosity
	REAL8 soildepth, d1, d2, d3; //guess
	REAL8 fc; //field capacity
	REAL8 deficit; //soil water deficit to reach field capacity in m
    REAL8 actArea = 0; //actual area [m2] for each cell yangx 2020-11
    REAL8 Outletpond = 0; //[m3/s] the remaining ponding water at the outlet
	//surface routing parameters
    //REAL8 srfrnoff = 0; //the direct overland surface flow to channel
	REAL8 ponding = 0;
    REAL8 chan_store = 0; //[m3]
    REAL8 benthic_area = 0; // [m2]
    //REAL8 paved_ponding = 0; //surface runoff
    REAL8 pro_srfrnoff = 0; 
	REAL8 theta1 = 0;
	REAL8 theta2 = 0;
	REAL8 theta3 = 0;
	REAL8 f = 0;
	REAL8 F = 0;
	//REAL8 ca = 0; //catchment area
	REAL8 gw = 0; //gravitational water
	REAL8 returnflow = 0; //flow from gw in excess of the available soil storage
	//REAL8 maxR = 0; //maximum gravitational water possible
	REAL8 qc = 0; // water transfered from the subsurface system to the channel
	REAL8 qall = 0; //lateral inflows to channel
	REAL8 Qij1 = 0; //new discharge from the upstream boundary
	REAL8 Qk1 = 0; //new discharge out of the cell
	REAL8 Si1j1 = 0; //storage in channel at the end of time step

	REAL8 leak = 0;
	
	REAL8 extragw_all = 0; //storage in extra deep GW--yangx 2020-05
	REAL8 ex_Fhydro = 0; 
	REAL8 extragw = 0;  //hydrologically active part of deep extra GW storage
	REAL8 ex_qc = 0;
	REAL8 ex_hj1i1 = 0;
	REAL8 ex_alpha, ex_qj1i, ex_R;

	dtdx = dt / _dx;

	// Reinitialize to zero the fluxes modified earlier / in the previous time step
	_FluxExfilt->reset();
	_FluxL2toL1->reset();
	_FluxL3toL2->reset();
	if(ctrl.sw_trck){
	  _FluxSrftoL1->reset();
	  _FluxL1toL2->reset();
	  _FluxL2toL3->reset();
	  if(ctrl.sw_2H)
	    trck.resetFd2HLat();
	  if(ctrl.sw_18O)
	    trck.resetFd18OLat();
	  if(ctrl.sw_Age)
	    trck.resetFAgeLat();
	}
	// Initialize others
	_FluxLattoSrf->reset();
	_FluxLattoChn->reset();
	_FluxLattoGW->reset();
    // Extra GW
	if(ctrl.sw_extraGW){
	  _FluxExtraGWtoLat->reset();
	  _FluxLattoExtraGW->reset();
	  _FluxExtraGWtoChn->reset();
	  _dailyExtraGwtrOutput.cells.clear();
	  _ExtraGWupstreamBC->reset();
	  if(ctrl.sw_2H)
	  //_Fd2HLattoExtraGW->reset();
	    trck.resetFd2HExtra();
      if(ctrl.sw_18O)
	  //_Fd18OLattoExtraGW->reset();
	    trck.resetFd18OExtra();
      if(ctrl.sw_Age)
	  //_FAgeLattoExtraGW->reset();
	    trck.resetFAgeExtra();
    }
  if(ctrl.sw_chan_evap)
    _FTemp_w->reset();
	// --------------------------------------------------------------------------------------
	for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
	  r = _vSortedGrid.cells[j].row;
	  c = _vSortedGrid.cells[j].col;
	  d = _vSortedGrid.cells[j].dir;

	  //surface routing stuff
	  //srfrnoff = 0;
	  returnflow = 0;
	  Qij1 = _Disch_upstreamBC->matrix[r][c];
	  qall = 0;
      Qk1 = 0;
      Outletpond = 0;
	  ponding = _ponding->matrix[r][c];
	  theta1 = _soilmoist1->matrix[r][c];
	  theta2 = _soilmoist2->matrix[r][c];
	  theta3 = _soilmoist3->matrix[r][c];
	  //ca = _catcharea->matrix[r][c];
	  gw = _GravityWater->matrix[r][c];
	  
	  fc = _fieldcapL3->matrix[r][c];
	  soildepth = _soildepth->matrix[r][c];
	  d1 = _depth_layer1->matrix[r][c];
	  d2 = _depth_layer2->matrix[r][c];
	  d3 = soildepth - d1 - d2;
	  //yangx 2020-11
      //get the actual area of specific grid cell
      actArea = _dx * _dx * _ttarea->matrix[r][c]; //[m2]
      pro_srfrnoff = 0; 
      //paved_ponding = 0;
      chan_store = 0;
      benthic_area = 0;
      //a proportion of surface flow goes to river, estimated by the relative lengths of channel and grid size
      if (ctrl.sw_channel && _channellength->matrix[r][c] > 0) {
        pro_srfrnoff = _channellength->matrix[r][c] / (_dx + _channellength->matrix[r][c]);
        //srfrnoff = ponding * pro_srfrnoff; //get the direct overland flows
        benthic_area = _channellength->matrix[r][c] * _channelwidth->matrix[r][c];
      }
	  //if reinfiltration switch is on is not a channel cell or the channel switch is off
	  //if (ctrl.sw_reinfilt && !(ctrl.sw_channel && _channelwidth->matrix[r][c] > 0)){
      if (ctrl.sw_reinfilt){
        //if its the channel cells, the part of ponding storage that goes to river dosen't reinfiltrate again
        //if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0){
		//  ponding = ponding*(1-pro_srfrnoff);		
	    Infilt_GreenAmpt(ctrl, f, F, theta1, theta2, theta3, ponding, gw, dt, r, c);
	    // Store time-step- and cumulated- infitlration flux (only if there is reinfiltration)
	    _FluxInfilt->matrix[r][c] += _FluxSrftoL1->matrix[r][c];
	    _AccInfilt->matrix[r][c] += _FluxSrftoL1->matrix[r][c];

     	// Tracking
	    if(ctrl.sw_trck){
	    // Mixing across the profile, accounting for snowmelt and lateral input
	      trck.MixingV_down(*this, ctrl, d1, d2, d3, fc, leak, r, c, 1);

	    // Back up soil moisture before vertical redistrib
	      _soilmoist1->matrix[r][c] = theta1; //soil moisture at t=t+1
	      _soilmoist2->matrix[r][c] = theta2;
	      _soilmoist3->matrix[r][c] = theta3;
	    // Back up ponding before GW seepage (for rivers)
	      _ponding->matrix[r][c] = ponding;
	    }
	  }	  
	  // For the rest of the routine, theta3 is only the content of the non-saturated
	  // part of L3
	  if (theta3 > fc) {
	    gw = (theta3 - fc) * d3;
	    theta3 = fc;
	  } else
	    gw = 0;
	  

	  if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) { 
	    //if this is a channel cell and channels are activated
	    //maxR = ( _porosity->matrix[r][c] - fc ) * soildepth; 
	    //calculates the maximum gravitational water that can go

	    qc = _KsatL3->matrix[r][c] * gw
	      * (1 - expl(-_chGWparam->matrix[r][c] * gw));
        // to update gw
	    gw -= qc * dtdx * _channellength->matrix[r][c]/_dx; 
		
		
        // Darcy's flux Q = qc*_chanlength
        // to update gw = gw - Q/(grid_area: _dx*_dx) * dt
	  }

	  _GravityWater->matrix[r][c] = gw;
	  
	  //enter groundwater
	  poros1 = _porosityL1->matrix[r][c];
	  poros2 = _porosityL2->matrix[r][c];
	  poros3 = _porosityL3->matrix[r][c];
	  alpha = _KsatL3->matrix[r][c] * sin(atan(_slope->matrix[r][c]));
	  
	  deficit = 0;
	  if (fabs(fc - theta3) > RNDOFFERR) {
	    deficit = (fc - theta3) * d3;
	  }
	  
	  // discharge (j is timestep) so j1i is total upstream discharge per unit width at t+1
	  qj1i = _GWupstreamBC->matrix[r][c];	
	  //Not used since local GW head is embedded in the updated theta3 portion, becoming R	 
	  hji1 = 0;
	  //recharge to the groundwater system at the end of the time step in meters
	  R = _GravityWater->matrix[r][c]; 
	  //gravity water becomes groundwater
	  _GravityWater->matrix[r][c] = 0; 
	  
	  // Solution of the kinematic wave (hj1i1 = head "ready to go" downstream)
	  hj1i1 = (dtdx * qj1i + hji1 + R - returnflow - deficit)
	    / (1 + alpha * dtdx); //R is in meters so no need to multiply times dt here
	  
	  
	  // If there's deficit and negative head -> capillary flow to L3, no GW outflow
	  if (deficit > 0 && hj1i1 < 0) {
	    theta3 += (dtdx * qj1i + hji1 + R - returnflow) / d3;
	    hj1i1 = 0;
	    // If there's deficit and positive head -> capillary flow to L3
	  } else if (deficit > 0 && hj1i1 >= 0){
	    theta3 += (dtdx * qj1i + hji1 + R - returnflow
		       - hj1i1 * (1 + alpha * dtdx)) / d3;
	  }
	  
	  // If the new amount of water in the cell is larger than the soil storage:
	  // --> Solve for return flow
	  if (((poros3 - theta3) * d3) < hj1i1) {
	    returnflow = -(poros3 - theta3) * d3 * (1 + alpha * dtdx)	\
	      + (dtdx * qj1i) + hji1 + R - deficit;
	    theta2 += returnflow / d2;
	    
	    // Tracking
	    //if(ctrl.sw_trck)
	    _FluxL3toL2->matrix[r][c] = returnflow ;
	    
	    if (theta2 > poros2) {
	      theta1 += (theta2 - poros2) * d2 / d1;
	      // Tracking
	      //if(ctrl.sw_trck)
		_FluxL2toL1->matrix[r][c] = (theta2 - poros2) * d2 ;
	      theta2 = poros2;
	    }
	    if (theta1 > poros1) {
	      ponding += (theta1 - poros1) * d1;
	      _FluxExfilt->matrix[r][c] = (theta1 - poros1) * d1 ;
	      theta1 = poros1;
	    }
	    
	    //ponding += returnflow;
	    // Update head
	    hj1i1 = (poros3 - theta3) * d3;
	  }
	  
	  //Extra GW terrestrial dynamics yangx 2020-05
	  if (ctrl.sw_extraGW){
		extragw_all = _ExtraGW->matrix[r][c];
		ex_Fhydro = _Hydrofrac_ExtraGW->matrix[r][c];
		//hydrologically active part only
		extragw = extragw_all * ex_Fhydro;
	    //ExtraGWDynamics(ctrl, extragw, ex_qc, ex_hj1i1, dt, r, c);
		//ex_qc = 0;   //[m2/s]
		//ex_hj1i1 = 0;


		//get the global variable  
		leak = _BedrockLeakageFlux->matrix[r][c]; //stored as m/s in in ln 155 SoilWaterRedistribution.cpp
		leak = leak * dt; //[m]
		
		//update storage
		extragw += leak;
		//baseflow generation if it is channel cell
		if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) { 
		  //if this is a channel cell and channels are activated
		  ex_qc = _KsatL3->matrix[r][c] * extragw * (1 - expl(-_chExGWparam->matrix[r][c] * extragw));//_chExGWparam->matrix[r][c];
		  //similarily as qc yangx 2020-11
		  //extragw -= ex_qc * dtdx;
          extragw -= ex_qc * dtdx * _channellength->matrix[r][c]/_dx;
		}
		
		//terrestrial routing factor, the same as original GW 
		ex_alpha = _KsatL3->matrix[r][c] * sin(atan(_slope->matrix[r][c]));	
		
		//input from up-slope grid cell
		ex_qj1i = _ExtraGWupstreamBC->matrix[r][c]; //[m2/s]
		// DeepGW head
		ex_R = extragw;
		//
		ex_hj1i1 = (dtdx*ex_qj1i + ex_R)/(1+ex_alpha*dtdx); //[m]
			
	  }
	  //if there is paved area presents, precip in this part
      // otherwise no surface lateral exchange nor surface runoff generation in channel cells
      //yangx 2020-11
      //if (_sealedarea->matrix[r][c] > 0.000001){
      //  paved_ponding = _sealedarea->matrix[r][c] * atm.getPrecipitation()->matrix[r][c];
      //  ponding = ponding * (1- _sealedarea->matrix[r][c] ) + paved_ponding; //update ponding storage
      //}//else{
      //  srfflw = 0.0;
      //}
	  //channel routing
	  if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) {
	    //ponding += qc * dtdx;
        //convert global var. [m] to [m3] to ensuring mass balance yangx 2020-11
		//link to line 489, here dx*dx is just for unit transform
        chan_store = _chan_store->matrix[r][c] * (_dx * _dx);
        //the overland flow part of the ponding water will go to channel directly
        //add surface runoff
        chan_store +=  ponding * pro_srfrnoff * actArea;
        //add subsurface runoff (gw)		
        chan_store += qc * dt * _ttarea->matrix[r][c] * _channellength->matrix[r][c];  //[m3]
		//also add baseflow --Extra GW yangx
		if (ctrl.sw_extraGW)
		  //ponding += ex_qc * dtdx;
          chan_store += ex_qc * dt * _ttarea->matrix[r][c] * _channellength->matrix[r][c];
        //Changed by yangx 2020-11
        //ALSO should be adjusted in KenematicWave.cpp	    
	    //qall = ponding * _dx / dt;
	    qall = chan_store / dt; //[m3/s]
		
	    KinematicWave(Qk1, Si1j1, Qij1, qall, dt, r, c);
	    //Qk1 = ponding * _dx*_dx/dt  + Qij1 ; // oooold (before kinematic wave)

	    // Not all ponding water get routed once in the channel..
	    // what is the actual amount of run-off that contributes to streamflow then?
	    // For lack of a better solution, it is the amount of ponding that is effectively routed
	    // AFTER all groundwater and upstream streamflow has been used
	    // (since groundwater effectively enters the stream)
	    _FluxGWtoChn->matrix[r][c] = qc*dtdx * _channellength->matrix[r][c]/_dx; //accounting for river length	    
	    //yangx 2020-11
		//_FluxSrftoChn->matrix[r][c] = std::max<double>(0.0,(Qk1 - (Qij1 + qc*_dx))*dtdx/_dx);
	    _FluxSrftoChn->matrix[r][c] = ponding * pro_srfrnoff; //[m]

        //also baseflow --Extra GW yangx
		if (ctrl.sw_extraGW){
		  _FluxExtraGWtoChn->matrix[r][c] = ex_qc*dtdx*_channellength->matrix[r][c]/_dx;
          _AccExtraGWtoChn->matrix[r][c] += _FluxExtraGWtoChn->matrix[r][c];
		  //update surface to channel flow
		  //_FluxSrftoChn->matrix[r][c] = std::max<double>(0.0,(Qk1 - (Qij1 + qc*_dx + ex_qc*_dx))*dtdx/_dx);
	    }
		
	    // Accumulated fluxes
	    _AccGWtoChn->matrix[r][c] += _FluxGWtoChn->matrix[r][c];
	    _AccSrftoChn->matrix[r][c] += _FluxSrftoChn->matrix[r][c];	
		
		
        chan_store = 0;
	    
	  }	  

	  // Locate downstream cell (if it exists)
	  lat_ok = 0;
	  switch (d) 
	    {
	    case 1:
	      rr = r+1;
	      cc = c-1;
	      lat_ok = 1;
	      break;
	    case 2:
	      rr = r+1;
	      cc = c;
	      lat_ok = 1;
	      break;
	    case 3:
	      rr = r+1;
	      cc = c+1;
	      lat_ok = 1;
	      break;
	    case 4:
	      rr = r;
	      cc = c-1;
	      lat_ok = 1;
	      break;
	    case 5: //if it is an outlet store the outflow m3s-1
	      _dailyGwtrOutput.cells.push_back(cell(r, c, (alpha * hj1i1 * _dx)));
	      //yangx 2020-11
		  //_dailyOvlndOutput.cells.push_back(cell(r, c, Qk1+ponding * _dx *_dx / dt));
          Outletpond = ponding*(1-pro_srfrnoff)*actArea/ dt; //[m3/s]directly put into final outputs
	      _dailyOvlndOutput.cells.push_back(cell(r, c, Qk1+Outletpond));
	      //second term needed to account for outer at outlets with no channel
          if (ctrl.sw_extraGW)
            _dailyExtraGwtrOutput.cells.push_back(cell(r, c, (alpha * ex_hj1i1 * _dx))); //[m3/s]			  
	      
	      break;
	    case 6:
	      rr = r;
	      cc = c+1;
	      lat_ok = 1;
	      break;
	    case 7:
	      rr = r-1;
	      cc = c-1;
	      lat_ok = 1;
	      break;
	    case 8:
	      rr = r-1;
	      cc = c;
	      lat_ok = 1;
	      break;
	    case 9:
	      rr = r-1;
	      cc = c+1;
	      lat_ok = 1;
	      break;
	    default:
	      return -1;
	    }
	  
	  // Check there is downstream cell
	  if(lat_ok){

	    // Input water for downstream cells (additive)
	    _FluxLattoSrf->matrix[rr][cc] += ponding * (1 - pro_srfrnoff) ;
	    _FluxLattoChn->matrix[rr][cc] += Qk1*dtdx/_dx;
	    _FluxLattoGW->matrix[rr][cc] += hj1i1 * alpha * dtdx;
	    // Accumulated fluxes
	    _AccLattoSrf->matrix[rr][cc] += ponding * (1 - pro_srfrnoff) ; //ponding
	    _AccLattoChn->matrix[rr][cc] += Qk1*dtdx/_dx ;
	    _AccLattoGW->matrix[rr][cc] += hj1i1 * alpha * dtdx;
	    
		// Extra GW --yangx
		// we use the same equation for terrestrial routing, so alpha = ex_alpha
		if (ctrl.sw_extraGW){
		  _FluxLattoExtraGW->matrix[rr][cc] += ex_hj1i1 * alpha * dtdx; //[m]
	      _AccLattoExtraGW->matrix[rr][cc] += ex_hj1i1 * alpha * dtdx;
		  _ExtraGWupstreamBC->matrix[rr][cc] += ex_hj1i1 * alpha; //[m2/s]
	    }

	    // Add the previously calculated *discharge* (not elevation) to the downstream cell
	    _GWupstreamBC->matrix[rr][cc] += hj1i1 * alpha;
	    _Disch_upstreamBC->matrix[rr][cc] += Qk1;  //includes newly added baseflow
	    _ponding->matrix[rr][cc] += ponding * (1 - pro_srfrnoff) ;     	

	  }

	  // Outgoing water (outside of lat_ok because can be 0)
	  //if(ctrl.sw_trck){
	  _FluxSrftoLat->matrix[r][c] = lat_ok != 0 ? ponding * (1 - pro_srfrnoff) : 0.0;
	  _FluxGWtoLat->matrix[r][c] = lat_ok != 0 ? hj1i1 * alpha * dtdx : 0.0;
      //Aaron Smith separated channel store
      _FluxChntoLat->matrix[r][c] = lat_ok != 0 ? Qk1*dtdx / _dx : 0.0;	  
	  // Accumulated fluxes
	  _AccSrftoLat->matrix[r][c] += _FluxSrftoLat->matrix[r][c];
	  _AccGWtoLat->matrix[r][c] += _FluxGWtoLat->matrix[r][c];
      //Aaron Smith separated channel store
      _AccChntoLat->matrix[r][c] += _FluxChntoLat->matrix[r][c];
      //outgoing water from Extra GW --yangx
	  if (ctrl.sw_extraGW){
		_FluxExtraGWtoLat->matrix[r][c] = lat_ok != 0 ? ex_hj1i1 * alpha * dtdx : 0.0; //[m]
	    _AccLattoExtraGW->matrix[r][c] += _FluxExtraGWtoLat->matrix[r][c];
	  }      	  

	  // Tracking of lateral in/out + return + seepage
	  // modified by yangx for extra GW
	  if(ctrl.sw_trck){
	    if (!ctrl.sw_extraGW){
		  trck.MixingV_latup(*this, ctrl, d1, d2, d3, fc, 
			         Qk1, dtdx, _dx, r, c);
	      // Tracking lateral inputs to the downstream cell
	      // Summed tracking contribution downstream cells (for mixing)
	      if(lat_ok == 1)
	        trck.FCdownstream(*this, ctrl, Qk1, dtdx, _dx, r, c, rr, cc);
	      else
	        // Catchment outlets' values
	        trck.OutletVals(*this, ctrl, 1,  Qk1, Outletpond, r, c);
		}else{  
		  //extra GW --yangx 2020-05
	      trck.MixingV_plusExtraGW(*this, ctrl, d1, d2, d3, fc,Qk1, 
		                          dtdx, _dx, dt, r, c);
		  if(lat_ok == 1){
			trck.FCdownstream_plusExtraGW(*this, ctrl, Qk1, dtdx, _dx, r, c, rr, cc);
		  }else{
			trck.OutletVals_plusExtraGW(*this, ctrl, 1, Qk1, Outletpond, r, c);
		  }
		}
	  }

	  // Update ponding and water contents
	  if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0){
        _ponding->matrix[r][c] = 0.0;
	    _chan_store->matrix[r][c] = Si1j1/(_dx * _dx); // link to ln314, just unit transform and store as [m]
        if (ctrl.sw_chan_evap){
	      SolveChannelEnergyBalance(atm,ctrl,trck,0,Qk1,r,c);
	      if (lat_ok){
	        _FTemp_w->matrix[rr][cc] += Qk1*dtdx/_dx * _Temp_w->matrix[r][c];
	      }
        } else {
	      _chan_evap->matrix[r][c] = 0.0; // no channel evaporation
	      if (ctrl.sw_trck)
	      trck.MixingV_evapW(atm,*this,ctrl,atm.getTemperature()->matrix[r][c],
			     _chan_evap->matrix[r][c],_chan_store->matrix[r][c],r,c);//to get a value for evap isotopes (no frac)
        }
	  }else{
	    _ponding->matrix[r][c] = 0.0;
        _chan_store->matrix[r][c] = 0.0;
        _chan_evap->matrix[r][c] = 0.0;
        if (ctrl.sw_chan_evap){
	       _Temp_w->matrix[r][c] = 0.0;
        }
      }
	  _soilmoist1->matrix[r][c] = theta1;
	  _soilmoist2->matrix[r][c] = theta2;
	  _soilmoist3->matrix[r][c] = theta3 + hj1i1 / d3;
	  _GrndWater->matrix[r][c] = hj1i1;
	  if (ctrl.sw_extraGW)
        _ExtraGW->matrix[r][c] = ex_hj1i1 + extragw_all*(1-ex_Fhydro);
	  // Save river discharge
	  _Disch_old->matrix[r][c] = Qk1;
	  Qk1 = 0;

	  // Accumulated fluxes
	  _AccExfilt->matrix[r][c] += _FluxExfilt->matrix[r][c];
	  
/* 	  // Deep groundwater subroutines for considering baseflow
	  if (ctrl.sw_extraGW)
	    ExtraGWDynamics(ctrl, trck, lat_ok, dt, r, c)
 */
	}
	
	// Save previous GW and surface state
	*_GrndWater_old = *_GrndWater;
	*_ponding_old = *_ponding;
	*_chan_store_old = *_chan_store;
	// Save previous Extra GW states
	if (ctrl.sw_extraGW){
	  *_FluxLattoExtraGW_old = *_FluxLattoExtraGW;
	  *_FluxExtraGWtoLat_old = *_FluxExtraGWtoLat;
	  *_FluxExtraGWtoChn_old = *_FluxExtraGWtoChn;
	  *_BedrockLeakageFlux_old = *_BedrockLeakageFlux;
      //to backup tracer fluxes: puls "_old"
	  if(ctrl.sw_trck){
	    trck.PrefluxTrck_plusExtraGW(ctrl);
	  }
	}	

	return EXIT_SUCCESS;
}
