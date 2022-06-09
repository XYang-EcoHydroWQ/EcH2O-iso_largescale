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
 * Basin.h
 *
 *  Created on: Oct 9, 2009
 *      Author: Marco Maneta
 *
 */

#ifndef BASIN_H_
#define BASIN_H_

#include "Atmosphere.h"
#include "ConstAndFuncs.h"
#include "Forest.h"
#include "Grid.h"
#include "InitConf.h"
#include "SortGrid.h"
#include "Tracking.h"

#include <omp.h>
#include <errno.h>

using namespace std;
class Forest;
class Tracking;
class Basin {

  /*Properties of the simulation*/

  UINT4 _NRows;
  UINT4 _NCols;
  REAL8 _dx;
  REAL8 _north,_south,_west,_east;
  REAL8 _nodata;

  //DEM map, which acts as the base map setting the geometry of the basin
  grid *_DEM;

  vectCells _vSortedGrid;

  Forest *fForest;

  /*Spatial properties of the basin*/
  grid *_porosity0; // Reference maximum volumetric water content (m3.m-3)
  // Layer-integrated maximum volumetric water content in L1, L2, and L3 (m3.m-3)
  grid *_porosityL1, *_porosityL2, *_porosityL3 ; // poros = poros0*exp(-z/kporos)
  grid *_kporos; // coefficient for exponential profile of porosity, in m 
  grid *_Ksat0; // Reference saturated hydraulic conductivity in m s-1
  // Layer-integrated saturated hydraulic conductivity in L1, L2, and L3 (m.s-1)
  grid *_KsatL1, *_KsatL2, *_KsatL3 ; // Ksat = Ksat0*exp(-z/kKsat)
  grid *_kKsat; // coefficient for exponential profile of Ksat, in m
  grid *_KvKs; //vertical to horiontal Ks anisotropy ratio
  grid *_random_roughness; //terrain roughness to calculate aerodynamic resistance (m)
  grid *_theta_r; //residual soil moisture volumetric
  grid *_psi_ae; // soil air entry pressure in m
  grid *_BClambda; //Brooks and Corey lambda parameter
  grid *_soildepth; //soil depth m
  grid *_depth_layer1; //depth of layer 1. 0.1 m by default
  grid *_depth_layer2; //depth of layer 2. Depth of layer 3 is calculated form depth
  //grid *_Kroot; // exponential root profile shape (m-1)
  //grid *_rootfrac1; //fraction of roots in soil layer 1
  //grid *_rootfrac2; //fraction of roots in soil layer 2. For layer three it is calculated from layer 1 and 2
  grid *_fieldcapL1, *_fieldcapL2, *_fieldcapL3; //field capacity layer 1 (volumetric)
  grid *_paramWc; //empirical parameter in water efficiency function for GPP calculation (see Landsberg and Waring, 1997 or TRIPLEX paper)
  grid *_paramWp; ////empirical parameter in water efficiency function for GPP calculation (see Landsberg and Waring, 1997 or TRIPLEX paper)
  grid *_meltCoeff; //snowmelt coefficient m s-1 C-1
  grid *_channelwidth; //width of channel in m. 0 if no channel
  grid *_chGWparam; //subsurface to channel water transfer parameter [dimensionless]

  grid *_Manningn; //manning's n

  grid *_ldd; //local drain direction (steepest 8 neighbor algorithm)
  grid *_catcharea; //catchment area (m2)
  grid *_slope; //slope in mm-1
  //yangx 2020-11
  grid *_ttarea; // areal proportion for each grid cell 
  grid *_channellength; //channel length in m. 0 if no channel
  grid *_chan_store; //channel storage (m3)
  grid *_chan_store_old; //channel storage for previous step (m3)
  grid *_sealedarea; //threshold of ponding water for surface flow generation (m)
  grid *_dampdepth; // soil depth at which there is no diurnal temperature variation
  grid *_Temp_d; //temperature at damping depth

  grid *_Temp_w;
  grid *_FTemp_w;
  grid *_chan_evap;
  grid *_chan_roughness; //channel roughness to calculate aerodynamic resistance (m)  
  /*State variables*/

  grid *_albedo; //surface albedo, no units
  grid *_emiss_surf; //emissivity of surface, no units
  grid *_soil_dry_heatcap; //dry soil heat capacity Jm-3C-1
  grid *_soil_dry_thermcond; //dry soil thermal conductivity Wm-1C-1
  grid *_snow; //snow water equivalent in m
  grid *_ponding; // water ponding on the soil surface in m
  grid *_infilt_cap; //infilt capacity m s-1
  grid *_IsSaturated; // Saturated locations
  grid *_soilmoist1; // average volumetric soil moisture over layer 1 or over entire soil profile
  grid *_soilmoist2; //average volumetric soil moisture of the second soil layer
  grid *_soilmoist3; //average volumetric soil moisture of the bottom soil layer
  grid *_soilmoist_12; //average volumetric soil moisture of the two upper layers
  grid *_soilmoist_av; //average volumetric soil moisture of the entire soil profile
  grid *_SoilWaterDepth; //soil moisture depth (m) for entire soil profile
  grid *_SoilWaterDepthL1; //soil moisture depth (m) for the first layer
  grid *_SoilWaterDepthL2; //soil moisture depth (m) for the second layer
  grid *_SoilWaterDepthL3; //soil moisture depth (m) non-gravitational in the third layer
  grid *_WaterTableDepth; //reconstructed WTD (ignoring perched aquifers)
  grid *_SoilSatDeficit; //soil saturation deficit (1 full deficit - 0 saturation)
  grid *_AccumInfilt; //Accumulated infiltration m
  grid *_Evaporation; //actual evaporation and transpiration in m s-1
  grid *_EvaporationS_all; //actual soil evaporation in m s-1
  grid *_EvaporationI_all; //actual evaporation from summed vegetation in m s-1
  grid *_Transpiration_all; //transpiration from canopy in m s-1

  grid *_BedrockLeakageFlux; //water flux down the bottom of the soil in m s-1
  grid *_CanopyStorage; //current water stored in the canopy (m)
  grid *_GravityWater; //current water stored in the soil beyond field capacity (m) (percolation or water traveling in the vadose zone)
  grid *_ponding_old; //water stored at the surface at the beginning of the time step (0 if not channel, m)
  grid *_GrndWater_old; //water stored in the gw system at the beginning of the time step (m)
  grid *_GrndWater; //water stored in the gw system at the end of the time step (m)
  grid *_GWupstreamBC; //gw flux upstream boundary condition (m2.s-1)

  grid *_Rn; //Net radiation for the soil surface Wm-2
  grid *_Rn_sum; //Net radiation summed at the top of the canopy Wm-2
  grid *_latheat; //latent heat flux into the atmosphere Wm-2
  grid *_sensheat; //sensible heat into the atmosphere Wm-2
  grid *_grndheat; //ground heat flux Wm-2
  grid *_snwheat; //snow heat flux Wm-2
  grid *_Temp_s; //temperature of the surface in dg C
  grid *_Temp_s_old; //temperature of the surface in dg C in the previous time step

  grid *_Disch_old;// streamflow out of each channel cell at the beginning of hte time step (m3 s-1)
  grid *_Disch_upstreamBC; //upstream boundary condition (m3s-1)
  vectCells _dailyOvlndOutput; //vector containing water output for each cell with no drainage (ldd value of 5). The vectCell structure contains the row and col
  //of the cell with no output and the area draining to that cell m3s-1
  vectCells _dailyGwtrOutput; //vector containing water output for each cell with no drainage (ldd value of 5). The vectCell structure contains the row and col
  //of the cell with no output and the area draining to that cell m3s-1 ???

  grid *_bedrock_leak;

  // ---------------------------------------------------------------------------------------
  // Head and moisture defining tightly-bound -> mobile transition
  grid *_psi_MW ; // pressure in meters of head (1 m = 9804.139432 Pa)
  grid *_moist_MW1, *_moist_MW2 ;
  // Mobile water water fractions
  grid *_fracMW1; // of _soilmoist1
  grid *_fracMW2; // of _soilmoist2
  grid *_fracMW12; // average over 2 upper layers

  // Vertical (intra-cell) fluxes downwards
  grid *_FluxCnptoSrf; // canopy/sky to surface
  grid *_FluxCnptoSnow; // canopy/sky to snowpack
  grid *_FluxSnowtoSrf; // snowmelt to surface
  grid *_FluxSrftoL1; // infiltration from first layer (changes twice during time step)
  grid *_FluxInfilt; // timestep-summed infiltration from first layer
  grid *_FluxL1toL2; // percolation L1 to L2
  grid *_FluxL2toL3; // percolation L2 to L3
  //grid *_FluxL2toGW; // recharge L2 to groundwater
  //grid *_FluxL3toGW; // recharge L3 to groundwater
  grid *_FluxRecharge; // recharge to L3
  // Vertical (intra-cell) fluxes upwards
  grid *_FluxExfilt; // first layer to surface (return flow)
  grid *_FluxL2toL1; // capillary + return flow, L2 to L1
  grid *_FluxL3toL2; // capillary + return flow, L2 to L3
  grid *_FluxGWtoL2; // return flow, groundwater to L2
  grid *_FluxGWtoL3; // discharge, groundwater to L2 unsaturated
  // Other intra-cell fluxes
  grid *_FluxGWtoChn; // discharge, groundwater to channel
  grid *_FluxSrftoChn; // overland flow to channel
  // Lateral (inter-cell), for report only
  grid *_FluxLattoChn; //Channel input
  grid *_FluxLattoSrf; //Lateral overland input
  grid *_FluxLattoGW; //Lateral GW input 
  grid *_FluxSrftoLat; //Lateral overland output
  grid *_FluxGWtoLat; //Lateral overland output 
  grid *_FluxChntoLat; //channel transport
  // Cumulated fluxes, for report only
  grid *_AccInfilt; // Infiltration from surface to first layer
  grid *_AccExfilt; // Exfiltration from first layer to surface
  grid *_AccLattoChn; //Channel input
  grid *_AccLattoSrf; //Lateral overland input
  grid *_AccLattoGW; //Lateral GW input 
  grid *_AccSrftoLat; //Lateral overland output
  grid *_AccGWtoLat; //Lateral groundwater output
  grid *_AccGWtoChn; //Groundwater seepage to channel
  grid *_AccSrftoChn; //Overland to channel
  grid *_AccChntoLat; //channel transport
  //for extra GW yangx 2020-05
  grid *_chExGWparam; //baseflow to channel water transfer parameter [dimensionless]
  grid *_ExtraGW; //Extra GW storage for each grid cell
  grid *_Hydrofrac_ExtraGW; //hydrologically active fraction of the Extra GW

  grid *_FluxExtraGWtoChn; // Deep baseflow from the extra GW storge
  grid *_FluxLattoExtraGW; // Lateral deep GW input
  grid *_FluxExtraGWtoLat; // Lateral deep GW output
  grid *_AccExtraGWtoChn; // Accumulated baseflow export to channel
  grid *_AccLattoExtraGW; // Accumulated baseflow input
  grid *_AccExtraGWtoLat; // Accumulated baseflow output
  grid *_ExtraGWupstreamBC; //ExtraGW flux upstream boundary condition (m2.s-1) 

  grid *_FluxLattoExtraGW_old; // store value of previous step
  grid *_FluxExtraGWtoLat_old; // store value of previous step
  grid *_FluxExtraGWtoChn_old; // store value of previous step
  grid *_BedrockLeakageFlux_old; // 

  vectCells _dailyExtraGwtrOutput; //vector containing baseflow output for each cell with no drainage (ldd value of 5). The vectCell structure contains the row and col
  //of the cell with no output 

  // --------------------------------------------------------------------------------------

  vectCells SortGridLDD();

  void CheckMaps(Control &ctrl); 
  // Check maps mainly to make sure no nodata values are in the domain. 
  // Also sets slope to MIN_SLOPE for pixels with 0 slope.

  int CalcCatchArea();
  int CalcKsatLayers(Control &ctrl);
  int CalcPorosLayers(Control &ctrl);
  int CalcFieldCapacity();
  int CalcInitialStreamStorage();
  int CalcRootDistrib();
  int CalcTPDMoisture(Control &ctrl);

  double NetRad(Atmosphere &atm, const double &Ts, REAL8 Kbeers, REAL8 lai,
		REAL8 ec, REAL8 Tc, int row, int col);
		
  double NetRad_water(Atmosphere &atm, const double &Tw, int row, int col);

  double LatHeat(Atmosphere &atm, double beta, double ra, double rs,
		 double rc, const double &Ts, int row, int col);
  double SensHeat(Atmosphere &atm, double ra, const double &Ts, int row,
		  int col);
  double GrndHeat(Atmosphere &atm, Control &ctrl, const double &theta,
		  const double &Ts, const double &Td, int row, int col);
  double SnowHeat(Atmosphere &atm, Control &ctrl, const double &Ts, int row,
		  int col);
  double MeltHeat(Atmosphere &atm, Control &ctrl, const double &Ts,
		  const double &swe, const double &M, int row, int col);
  double RainHeat(Atmosphere &atm, double R, int row, int col);
  double SnowOutput(Atmosphere &atm, Control &ctrl, Tracking &trck,
		    const double &meltheat, int row, int col);

  REAL8 CalcAerodynResist(REAL8 u_za, REAL8 z_a, REAL8 z_0u, REAL8 z_du,
			  REAL8 z_0o, REAL8 z_do, REAL8 Ht, REAL8 LAI, REAL8 Ts, REAL8 Ta,
			  INT4 option, bool surface);
  REAL8 CalcSoilResist(double &theta, int row, int col, UINT4 option);

  double InOutMix_Temperature(double hold, double Told, double qin, double Tin, double qout, int mode);

  //Hydrologic processes

  void Infilt_GreenAmpt(Control &ctrl, double &f, double &F, double &theta, double &theta2,
			double &theta3, double &pond, double &gw, double dt, int r, int c);
  void SoilWaterRedistribution(Control &ctrl, const double &F, double &theta1,
			       double &theta2, double &theta3, double &pond, double &gw, double &leak, double dt,
			       int r, int c);
  void Infilt_Richards(Control &ctrl, double &f, double &F, double &theta1,
		       double &theta2, double &theta3, double &leak, double &pond, double &percolat,
		       double dt, int r, int c, int flowdir);
  int SolveSurfaceEnergyBalance(Atmosphere &atm, Control &ctrl, Tracking &trck,
				REAL8 ra,	REAL8 rs, REAL8 rc, REAL8 Kbeers, REAL8 lai, REAL8 emis_can,
				REAL8 Temp_can, REAL8 &nrad, REAL8 &latheat, REAL8 &sensheat,
				REAL8 &grndheat, REAL8 &snowheat, REAL8 &meltheat, REAL8 &Tsold,
				REAL8 &etp, REAL8 &pond, REAL8 &theta, REAL8 &Ts1, REAL8 &Tdold,
				REAL8 p, UINT4 r, UINT4 c, UINT4 s);

  int SolveChannelEnergyBalance(Atmosphere &atm, Control &ctrl, Tracking &trck,
				REAL8 chan_prec, REAL8 qout, UINT4 r, UINT4 c);

  void ChannelEvaporation(REAL8 LE, REAL8 lambda, REAL8 rs, REAL8 &evap, REAL8 &chan_store, REAL8 dt, UINT4 r, UINT4 c);


  //This functions updates soil moisture by solving the local soil water balance
  void SoilEvapotranspiration(REAL8 LE, //input latent heat
			      REAL8 Ts, //input surface temperature
			      REAL8 lambda, //input the latent heat (either latent heat of vaporization or of sublimation)
			      REAL8 rs, // input the potential exfiltration capacity
			      REAL8 &etp, //output updated evapotranspiration
			      REAL8 &theta, //output updated soil moisture
			      REAL8 dt, //time step
			      UINT4 r, UINT4 c);

  REAL8 ExfiltrationCapacity(REAL8 theta, //soil moisture
			     REAL8 dt, UINT4 r, UINT4 c);

  void CalcSoilMoistureProfile(Atmosphere &atm, Control &ctrl, REAL8 theta,
			       UINT4 row, UINT4 col);

  void KinematicWave(REAL8 &Qk1,  REAL8 &S,   REAL8 &Qij1,  REAL8 &qall,  REAL8 dt, UINT4 r, UINT4 c);

  //extra GW yangx 2020-05
  void ExtraGWDynamics(Control &ctrl, double extragw,
			  double ex_qc, double ex_hj1il,
		      double dt,int r, int c);	 

  /*	//This function updates _psi with the soil tension corresponding to the current soil moisture status
	void UpdateSoilWaterPotential() {
	int r, c;
	double S;
        #pragma omp parallel for\
	default(none) private(r,c,S)
	for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
	r = _vSortedGrid.cells[j].row;
	c = _vSortedGrid.cells[j].col;
	S = (_soilmoist1->matrix[r][c] - _theta_r->matrix[r][c])
	/ (_porosity->matrix[r][c] - _theta_r->matrix[r][c]);

	_psi->matrix[r][c] = fabs(_psi_ae->matrix[r][c])
	/ pow(S, _BClambda->matrix[r][c]);
	}
	}*/

 public:

  //Constructors
  Basin();
  Basin(Control &ctrl);

  //dtor
  ~Basin();

  int UpdateSnowPack(Atmosphere &atm, Control &ctrl);

  int SolveCanopyFluxes(Atmosphere &atm, Control &ctrl, Tracking &trck);
  int SolveSurfaceFluxes(Atmosphere &atm, Control &ctrl, Tracking &trck);
  int CalculateGrowForest(const Atmosphere &atm, const Control &ctrl);

  //int DailySurfaceRouting(Atmosphere &atm, Control &ctrl);
  int DailyGWRouting(Atmosphere &atm, Control &ctrl, Tracking &trck);
  int CalculateSatArea(Atmosphere &atm, Control &ctrl);

  int CalcFracMobileWater();

  //Getters

  REAL8 getCellSize() const {
    return _dx;
  }

  UINT4 getNumRows() const {
    return _NRows;
  }

  UINT4 getNumCols() const {
    return _NCols;
  }
  
  REAL8 getSouthCoord() const {
    return _south;
  }

  REAL8 getWestCoord() const {
    return _west;
  }

  REAL8 getModelnodata() const {
    return _nodata;
  }
  
  const vectCells &getSortedGrid() const {
    return _vSortedGrid;
  }

  const vectCells *getDailyOvlndOutput() const {
    return &_dailyOvlndOutput;
  }

  const vectCells *getDailyGwtrOutput() const {
    return &_dailyGwtrOutput;
  }

  grid *getDEM() const {
    return _DEM;
  }
  grid *getTTarea() const {
    return _ttarea;
  }  
  grid *getSoilDepth() const {
    return _soildepth;
  }
  grid *getSoilDepth1() const {
    return _depth_layer1;
  }
  grid *getSoilDepth2() const {
    return _depth_layer2;
  }

  grid *getNetRad() const {
    return _Rn;
  }

  grid *getNetRad_sum() const {
    return _Rn_sum;
  }

  grid *getLatheat() const {
    return _latheat;
  }

  grid *getSensHeat() const {
    return _sensheat;
  }

  grid *getGrndHeat() const {
    return _grndheat;
  }

  grid *getSnwHeat() const {
    return _snwheat;
  }

  grid *getSoilTemp() const {
    //return _Temp_s;
    return _Temp_d;
  }

  grid *getSkinTemp() const {
    return _Temp_s;
  }

  grid *getCanopyStorage() const {
    return _CanopyStorage;
  }
  grid *getSnowWaterEquiv() const {
    return _snow;
  }

  grid *getPondingWater() const {
    return _ponding;
  }

  grid *getChannelWater() const {
    return _chan_store;
  }
  
  grid *getChannel_old() const {
    return _chan_store_old;
  }
  grid *getTemp_w() const {
    return _Temp_w;
  }
  grid *getChanEvap() const {
    return _chan_evap;
  }

  grid *getPonding_old() const {
    return _ponding_old;
  }
  grid *getStreamflow() const {
    return _Disch_old;
  }
  grid *getSatArea() const {
    return _IsSaturated;
  }
  grid *getSoilMoist1() const {
    return _soilmoist1;
  }
  grid *getSoilMoist2() const {
    return _soilmoist2;
  }
  grid *getSoilMoist3() const {
    return _soilmoist3;
  }

  grid *getSoilMoist_av() const {

    double depth;
    double d1, d2, d3;
    int r, c;
#pragma omp parallel for			\
  default(none) private(r,c,depth, d1, d2, d3)
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      depth = _soildepth->matrix[r][c];
      d1 = _depth_layer1->matrix[r][c];
      d2 = _depth_layer2->matrix[r][c];
      d3 = depth - d1 - d2;
      _soilmoist_av->matrix[r][c] = (_soilmoist1->matrix[r][c] * d1
				     + _soilmoist2->matrix[r][c] * d2
				     + _soilmoist3->matrix[r][c] * d3) / depth;
    }

    return _soilmoist_av;
  }

  grid *getSoilMoist_12() const {

    double d1, d2;
    int r, c;
#pragma omp parallel for			\
  default(none) private(r,c, d1, d2)
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      d1 = _depth_layer1->matrix[r][c];
      d2 = _depth_layer2->matrix[r][c];
      _soilmoist_12->matrix[r][c] = (_soilmoist1->matrix[r][c] * d1
				     + _soilmoist2->matrix[r][c] * d2) / (d1+d2);
    }

    return _soilmoist_12;
  }

  // --- Two-domain stuff ------------------------------------
  grid *getMoistureMW1() const {
    return _moist_MW1;
  }
  grid *getMoistureMW2() const {
    return _moist_MW2;
  }
  grid *getFracMW1() const {
    return _fracMW1;
  }
  grid *getFracMW2() const {
    return _fracMW2;
  }
  grid *getFracMW12() const {

    double d1, d2;
    int r, c;
#pragma omp parallel for			\
  default(none) private(r,c, d1, d2)
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      d1 = _depth_layer1->matrix[r][c];
      d2 = _depth_layer2->matrix[r][c];
      _fracMW12->matrix[r][c] = (_fracMW1->matrix[r][c]*d1 + _fracMW2->matrix[r][c]*d2) / (d1+d2);
    }

    return _fracMW12;
  }
  // ----------------------------------------------------------


  grid *getSoilWaterDepth() const {
    int r, c;
    double depth;
    double d1, d2, d3;
#pragma omp parallel for			\
  default(none) private(r,c,depth, d1, d2, d3)
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      depth = _soildepth->matrix[r][c];
      d1 = _depth_layer1->matrix[r][c];
      d2 = _depth_layer2->matrix[r][c];
      d3 = depth - d1 - d2;
      _SoilWaterDepth->matrix[r][c] = (_soilmoist1->matrix[r][c] * d1
				       + _soilmoist2->matrix[r][c] * d2
				       + _soilmoist3->matrix[r][c] * d3);
    }

    return _SoilWaterDepth;
  }

  grid *getSoilWaterDepthL1() const {
    int r, c;
#pragma omp parallel for			\
  default(none) private(r,c)
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      _SoilWaterDepthL1->matrix[r][c] = _soilmoist1->matrix[r][c] * _depth_layer1->matrix[r][c];
    }	  
    return _SoilWaterDepthL1;
  }

  grid *getSoilWaterDepthL2() const {
    int r, c;
#pragma omp parallel for			\
  default(none) private(r,c)
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      _SoilWaterDepthL2->matrix[r][c] = _soilmoist2->matrix[r][c] * _depth_layer2->matrix[r][c];
    }	  
    return _SoilWaterDepthL2;
  }

  // Non-gravitational water depth in third layer
  grid *getSoilWaterDepthL3() const {
    int r, c;
    double fc3, d3;
#pragma omp parallel for			\
  default(none) private(r,c, fc3, d3)
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      fc3 = _fieldcapL3->matrix[r][c];
      d3 = _soildepth->matrix[r][c]-_depth_layer1->matrix[r][c]-_depth_layer2->matrix[r][c];
      _SoilWaterDepthL3->matrix[r][c] = min<double>(_soilmoist3->matrix[r][c],fc3) * d3;
    }	  
    return _SoilWaterDepthL3;
  }

  grid *getInitGroundwater() const {
    int r, c;
    double fc3, d3;
#pragma omp parallel for			\
  default(none) private(r,c, fc3, d3)
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      fc3 = _fieldcapL3->matrix[r][c];
      d3 = _soildepth->matrix[r][c]-_depth_layer1->matrix[r][c]-_depth_layer2->matrix[r][c];
      _GrndWater->matrix[r][c] = max<double>(0,_soilmoist3->matrix[r][c] - fc3) * d3;
    }	  
    return _GrndWater;
  }

  grid *getWaterTableDepth() const {
    /*int r, c;
      double depth, fc, eta;
      double d1, d2, d3;
      #pragma omp parallel for						\
      default(none) private(r, c, fc, eta, depth, d1, d2, d3)
      for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      fc = _fieldcap->matrix[r][c];
      eta = _porosity->matrix[r][c];
      depth = _soildepth->matrix[r][c];
      // If the theta3 is not above field cap, then no water table within the profile
      if(fc - _soilmoist3->matrix[r][c] > 0)
      _WaterTableDepth->matrix[r][c] = depth;
      else{
      d1 = _depth_layer1->matrix[r][c];
      d2 = _depth_layer2->matrix[r][c];
      d3 = depth - d1 - d2;
      // If theta3 below porosity, water table within third layer
      if(fabs(eta - _soilmoist3->matrix[r][c]) > RNDOFFERR)
      _WaterTableDepth->matrix[r][c] =
      d1 + d2 + d3*(eta - _soilmoist3->matrix[r][c])/(eta - fc);
      // If layer 3 saturated...
      else{
      // If theta2 below porosity, water table within third layer
      if(fabs(eta - _soilmoist2->matrix[r][c]) > RNDOFFERR)
      _WaterTableDepth->matrix[r][c] =
      d1 + d2*(eta - _soilmoist2->matrix[r][c])/(eta - fc);
      // If layer 2 saturated, water table within first layer
      else
      _WaterTableDepth->matrix[r][c] =
      d1*(eta - _soilmoist1->matrix[r][c])/(eta - fc);
      }
      }
      }*/
	  
    return _WaterTableDepth;
  }

  grid * getSaturationDeficit() const {

    int r, c;
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      _SoilSatDeficit->matrix[r][c] = 1 - ((_soilmoist1->matrix[r][c] - _theta_r->matrix[r][c])
	   / (_porosityL1->matrix[r][c] - _theta_r->matrix[r][c]));
    }

    return _SoilSatDeficit;

  }

  /*	grid *getSoilWaterPotential() {
	if (!_depth_layer1) //if we are using the lumped soil hydrology, update the soil water potential
	UpdateSoilWaterPotential();

	return _psi;
	}*/
  grid *getInfiltCap() const {
    return _infilt_cap;
  }
	
  //grid *getRootFrac1() const {
  //	return _rootfrac1;
  //}
  //grid *getRootFrac2() const {
  //	return _rootfrac2;
  //}
	
  grid *getEvaporation() const {
    return _Evaporation;
  }
  grid *getEvaporationS_all() const {
    return _EvaporationS_all;
  }
  grid *getEvaporationI_all() const {
    return _EvaporationI_all;
  }
  grid *getTranspiration_all() const {
    return _Transpiration_all;
  }
  grid *getBedrockLeakage() const{
    return _BedrockLeakageFlux;
  }

  grid *getParamWc() const {
    return _paramWc;
  }

  grid *getParamWp() const {
    return _paramWp;
  }

  grid *getPorosityL1() const {
    return _porosityL1;
  }
  grid *getPorosityL2() const {
    return _porosityL2;
  }
  grid *getPorosityL3() const {
    return _porosityL3;
  }
  grid *getSoilMoistR() const {
    return _theta_r;
  }
  grid *getPsiAE() const {
    return _psi_ae;
  }
  grid *getBClambda() const {
    return _BClambda;
  }

  grid *getFieldCapacityL1() const {
    return _fieldcapL1;
  }
  grid *getFieldCapacityL2() const {
    return _fieldcapL2;
  }
  grid *getFieldCapacityL3() const {
    return _fieldcapL3;
  }

  grid *getCatchmentArea() const {
    return _catcharea;
  }

  grid *getGravityWater() const {
    return _GravityWater;
  }

  grid *getGrndWater() const {
    return _GrndWater;
  }
  grid *getGrndWater_old() const {
    return _GrndWater_old;
  }

  // Internal fluxes
  grid *getFluxSrftoL1() const {
    return _FluxSrftoL1;
  }
  grid *getFluxInfilt() const {
    return _FluxInfilt;
  }
  grid *getFluxExfilt() const {
    return _FluxExfilt;
  }
  grid *getFluxRecharge() const {
    return _FluxRecharge;
  }
  grid *getFluxLattoGW() const {
    return _FluxLattoGW;
  }
  grid *getFluxLattoSrf() const {
    return _FluxLattoSrf;
  }
  grid *getFluxLattoChn() const {
    return _FluxLattoChn;
  }
  grid *getFluxGWtoLat() const {
    return _FluxGWtoLat;
  }
  grid *getFluxSrftoLat() const {
    return _FluxSrftoLat;
  }
  grid *getFluxGWtoChn() const {
    return _FluxGWtoChn;
  }
  grid *getFluxSrftoChn() const {
    return _FluxSrftoChn;
  }
  grid *getFluxChntoLat() const {
    return _FluxChntoLat;
  }
  grid *getAccChntoLat() const {
    return _AccChntoLat;
  }

  grid *getAccInfilt() const {
    return _AccInfilt;
  }
  grid *getAccExfilt() const {
    return _AccExfilt;
  }
  grid *getAccLattoGW() const {
    return _AccLattoGW;
  }
  grid *getAccLattoSrf() const {
    return _AccLattoSrf;
  }
  grid *getAccLattoChn() const {
    return _AccLattoChn;
  }
  grid *getAccSrftoLat() const {
    return _AccSrftoLat;
  }
  grid *getAccGWtoLat() const {
    return _AccGWtoLat;
  }
  grid *getAccGWtoChn() const {
    return _AccGWtoChn;
  }
  grid *getAccSrftoChn() const {
    return _AccSrftoChn;
  }

  // Addition tracking
  // ---------------------------------
  grid *getChannelWidth() const {
    return _channelwidth;
  }
  grid *getFluxCnptoSrf() const {
    return _FluxCnptoSrf;
  }
  grid *getFluxCnptoSnow() const {
    return _FluxCnptoSnow;
  }
  grid *getFluxSnowtoSrf() const {
    return _FluxSnowtoSrf;
  }
  grid *getFluxL1toL2() const {
    return _FluxL1toL2;
  }
  grid *getFluxL2toL3() const {
    return _FluxL2toL3;
  }
  //grid *getFluxL2toGW() const {
  //	return _FluxL2toGW;
  //}
  //grid *getFluxL3toGW() const {
  //  return _FluxL3toGW;
  //}
  grid *getFluxL2toL1() const {
    return _FluxL2toL1;
  }
  grid *getFluxL3toL2() const {
    return _FluxL3toL2;
  }
  grid *getFluxGWtoL2() const {
    return _FluxGWtoL2;
  }
  grid *getFluxGWtoL3() const {
    return _FluxGWtoL3;
  }
  
  // Extra GW related yangx 2020-05
  grid *getExtraGW() const {
    return _ExtraGW;
  }
  grid *getFracHydroExtraGW() const {
    return _Hydrofrac_ExtraGW;
  }
  grid *getFluxExtraGWtoChn() const {
    return _FluxExtraGWtoChn;
  }
  grid *getAccExtraGWtoChn() const {
    return _AccExtraGWtoChn;
  }
  grid *getFluxExtraGWtoLat() const {
    return _FluxExtraGWtoLat;
  }
  grid *getAccExtraGWtoLat() const {
    return _AccExtraGWtoLat;
  }
  grid *getFluxLattoExtraGW() const {
    return _FluxLattoExtraGW;
  }
  grid *getAccLattoExtraGW() const {
    return _AccLattoExtraGW;
  }  
  grid *getFluxExtraGWtoLat_old() const {
    return _FluxExtraGWtoLat_old;
  }
  grid *getFluxExtraGWtoChn_old() const {
    return _FluxExtraGWtoChn_old;
  }
  grid *getFluxLattoExtraGW_old() const {
    return _FluxLattoExtraGW_old;
  }
  grid *getBedrockLeakage_old() const{
    return _BedrockLeakageFlux_old;
  }
  //budget
  const vectCells *getDailyExtraGwtrOutput() const {
    return &_dailyExtraGwtrOutput;
  }
  
  // --------------------------------------------------------------------------------------
  // -- Getters of fForest getters

  UINT4 getNumSpecies() const ;

  grid *getVegetFrac(UINT4 n) const;

  grid *getLAI(UINT4 n) const;

  grid *getStemDensity(UINT4 n) const;

  grid *getStandAge(UINT4 n) const;

  grid *getCanopyCond(UINT4 n) const;

  grid *getGPP(UINT4 n) const;

  grid *getNPP(UINT4 n) const;

  grid *getBasalArea(UINT4 n) const;

  grid *getTreeHeight(UINT4 n) const;

  grid *getRootMass(UINT4 n) const;

  grid *getCanopyTemp(UINT4 n) const;

  grid *getCanopyNetRad(UINT4 n) const;

  grid *getCanopyLatHeatE(UINT4 n) const;

  grid *getCanopyLatHeatT(UINT4 n) const;

  grid *getCanopySensHeat(UINT4 n) const;

  grid *getCanopyWaterStor(UINT4 n) const;

  grid *getETspecies(UINT4 n) const;

  grid *getTranspiration(UINT4 n) const;

  grid *getEinterception(UINT4 n) const;

  grid *getEsoil(UINT4 n) const;

  grid *getLeafWaterPotential(UINT4 n) const;

  grid *getRootFrac1(UINT4 n) const;
  grid *getRootFrac2(UINT4 n) const;

  // Addition tracking
  // ---------------------------------
  // -- Getters of fTracking getters
  // 
  grid *getd2Hcanopy(UINT4 n) const ;
  grid *getd2HevapI(UINT4 n) const ;
  grid *getd2HevapT(UINT4 n) const ;
  grid *getd2HevapS(UINT4 n) const ;

  // 18O
  grid *getd18Ocanopy(UINT4 n) const ;
  grid *getd18OevapI(UINT4 n) const ;
  grid *getd18OevapT(UINT4 n) const ;
  grid *getd18OevapS(UINT4 n) const ;

  // Age
  grid *getAgecanopy(UINT4 n) const ;
  grid *getAgeevapI(UINT4 n) const ;
  grid *getAgeevapT(UINT4 n) const ;
  grid *getAgeevapS(UINT4 n) const ;
  //setters

  void setAgecanopy(UINT4 n, UINT4 r, UINT4 c, REAL8 value) const;
  void setAgeevapS(UINT4 n, UINT4 r, UINT4 c, REAL8 value) const;
  void setAgeevapI(UINT4 n, UINT4 r, UINT4 c, REAL8 value) const;
  void setAgeevapT(UINT4 n, UINT4 r, UINT4 c, REAL8 value) const;
};

#endif /* BASIN_H_ */
