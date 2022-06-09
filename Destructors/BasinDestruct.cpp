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
 * BasinDestruct.cpp
 *
 *  Created on: Oct 13, 2009
 *      Author: Marco Maneta
 */

#include "Basin.h"

Basin::~Basin(){

	if(_DEM)
		delete _DEM;
	if(_ldd)
		delete _ldd;
	if(_snow)
		delete _snow;
	if(_Rn)
		delete _Rn;
	if(_Rn_sum)
		delete _Rn_sum;
	if(_latheat)
		delete _latheat;
	if(_sensheat)
		delete _sensheat;
	if(_grndheat)
		delete _grndheat;
	if(_snwheat)
		delete _snwheat;
	if(_Temp_s)
		delete _Temp_s;
	if(_Temp_s_old)
		delete _Temp_s_old;
	if(_albedo)
		delete _albedo;
	if(_emiss_surf)
		delete _emiss_surf;
	if(_soil_dry_heatcap)
		delete _soil_dry_heatcap;
	if(_soil_dry_thermcond)
		delete _soil_dry_thermcond;
	if(_dampdepth)
		delete _dampdepth;
	if(_Temp_d)
		delete _Temp_d;
	if(_ponding)
		delete _ponding;
	if(_Ksat0)
		delete _Ksat0;
	if(_kKsat)
		delete _kKsat;
	if(_KsatL1)
		delete _KsatL1;
	if(_KsatL2)
		delete _KsatL2;
	if(_KsatL3)
		delete _KsatL3;
	if(_KvKs)
		delete _KvKs;
	if(_random_roughness)
		delete _random_roughness;
	if(_slope)
		delete _slope;
	if(_porosity0)
		delete _porosity0;
	if(_kporos)
		delete _kporos;
	if(_porosityL1)
		delete _porosityL1;
	if(_porosityL2)
		delete _porosityL2;
	if(_porosityL3)
		delete _porosityL3;
	if(_psi_ae)
		delete _psi_ae;
	if(_BClambda)
		delete _BClambda;
	if(_theta_r)
		delete _theta_r;
	if(_infilt_cap)
		delete _infilt_cap;
	if(_soilmoist1)
		delete _soilmoist1;
	if(_AccumInfilt)
		delete _AccumInfilt;
	if(_soildepth)
		delete _soildepth;
	if(_depth_layer1)
		delete _depth_layer1;
	if(_depth_layer2)
		delete _depth_layer2;
	//if(_rootfrac1)
	//	delete _rootfrac1;
	//if(_rootfrac2)
	//      delete _rootfrac2;
	//if(_Kroot)
	//	delete _Kroot;
	if(_fieldcapL1)
		delete _fieldcapL1;
	if(_fieldcapL2)
		delete _fieldcapL2;
	if(_fieldcapL3)
		delete _fieldcapL3;
	if(_paramWc)
		delete _paramWc;
	if(_paramWp)
		delete _paramWp;
	if(_meltCoeff)
		delete _meltCoeff;
	if(_Evaporation)
		delete _Evaporation;
	if(_BedrockLeakageFlux)
		delete _BedrockLeakageFlux;
	if(_SoilWaterDepth)
		delete _SoilWaterDepth;
	if(_SoilWaterDepthL1)
		delete _SoilWaterDepthL1;
	if(_SoilWaterDepthL2)
		delete _SoilWaterDepthL2;
	if(_SoilWaterDepthL3)
		delete _SoilWaterDepthL3;
	if(_WaterTableDepth)
		delete _WaterTableDepth;
	if(_SoilSatDeficit)
		delete _SoilSatDeficit;
	if(_CanopyStorage)
		delete _CanopyStorage;
	if(_Disch_old)
		delete _Disch_old;
	if(_Disch_upstreamBC)
		delete _Disch_upstreamBC;
	if(_catcharea)
		delete _catcharea;
	if(_GravityWater)
		delete _GravityWater;
	if(_GrndWater_old)
		delete _GrndWater_old;
	if(_GrndWater)
		delete _GrndWater;
	if(_GWupstreamBC)
		delete _GWupstreamBC;
	if(_channelwidth)
		delete _channelwidth;
	if(_chGWparam)
		delete _chGWparam;
	if(_Manningn)
		delete _Manningn;
	if(_soilmoist_av)
		delete _soilmoist_av;
	if(_soilmoist_12)
		delete _soilmoist_12;
	if(_soilmoist2)
		delete _soilmoist2;
	if(_soilmoist3)
		delete _soilmoist3;
	if(_bedrock_leak)
		delete _bedrock_leak;
	if(_IsSaturated)
		delete _IsSaturated;
	if(_EvaporationS_all)
	  delete _EvaporationS_all;
	if(_EvaporationI_all)
	  delete _EvaporationI_all;
	if(_Transpiration_all)
	  delete _Transpiration_all;

    if(_Temp_w)
      delete _Temp_w;
    if(_FTemp_w)
      delete _FTemp_w;
    if(_chan_evap)
      delete _chan_evap;
    if(_chan_roughness);
      delete _chan_roughness;

	if(fForest)
		delete fForest;

	if(_ponding_old)
	  delete _ponding_old;
	if(_FluxCnptoSrf)
		delete _FluxCnptoSrf;
	if(_FluxCnptoSnow)
		delete _FluxCnptoSnow;
	if(_FluxSnowtoSrf)
		delete _FluxSnowtoSrf;
	if(_FluxSrftoL1)
		delete _FluxSrftoL1;
	if(_FluxInfilt)
		delete _FluxInfilt;
	if(_FluxExfilt)
		delete _FluxExfilt;
	if(_FluxRecharge)
		delete _FluxRecharge;
	if(_FluxLattoSrf)
		delete _FluxLattoSrf;
	if(_FluxLattoGW)
		delete _FluxLattoGW;
	if(_FluxLattoChn)
		delete _FluxLattoChn;
	if(_FluxSrftoLat)
		delete _FluxSrftoLat;
	if(_FluxGWtoLat)
		delete _FluxGWtoLat;
	if(_FluxGWtoChn)
		delete _FluxGWtoChn;
	if(_FluxSrftoChn)
		delete _FluxSrftoChn;
	if(_AccInfilt)
		delete _AccInfilt;
	if(_AccExfilt)
		delete _AccExfilt;
	if(_AccLattoGW)
		delete _AccLattoGW;
	if(_AccLattoSrf)
		delete _AccLattoSrf;
	if(_AccLattoChn)
		delete _AccLattoChn;
	if(_AccSrftoLat)
		delete _AccSrftoLat;
	if(_AccGWtoLat)
		delete _AccGWtoLat;
	if(_AccGWtoChn)
		delete _AccGWtoChn;
	if(_AccSrftoChn)
		delete _AccSrftoChn;

	// Tracking
	if(_psi_MW)
	  delete _psi_MW;
	if(_moist_MW1)
	  delete _moist_MW1;
	if(_moist_MW2)
	  delete _moist_MW2;
	if(_fracMW1)
	  delete _fracMW1;
	if(_fracMW2)
	  delete _fracMW2;
	if(_fracMW12)
	  delete _fracMW12;
	if(_FluxL1toL2)
		delete _FluxL1toL2;
	if(_FluxL2toL3)
		delete _FluxL2toL3;
	//if(_FluxL2toGW)
	//		delete _FluxL2toGW;
	//if(_FluxL3toGW)
	//		delete _FluxL3toGW;
	if(_FluxL2toL1)
		delete _FluxL2toL1;
	if(_FluxL3toL2)
		delete _FluxL3toL2;
	//if(_FluxGWtoL2)
	//	delete _FluxGWtoL2;
	//	if(_FluxGWtoL3)
	//	delete _FluxGWtoL3;

}
