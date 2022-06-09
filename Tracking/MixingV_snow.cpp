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
 * MixingV_snow.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingV_snow(Atmosphere &atm, Basin &bsn, Control &ctrl,
			    double &h, double &dh, int r, int c) //time step
{
  
  double h_eff; // Effective SWE used for mixing
  double snow_in = bsn.getFluxCnptoSnow()->matrix[r][c];

  // - in snowpack (snowfall in + snowmelt out), 
  // considering that snowmelt "flushes" the most recent snowfall first, without mixing
  if(abs(h) < RNDOFFERR){
    if(ctrl.sw_2H)
      _d2Hsnowpack->matrix[r][c] = -1000;
    if(ctrl.sw_18O)
      _d18Osnowpack->matrix[r][c] = -1000;
    if(ctrl.sw_Age)
      _Agesnowpack->matrix[r][c] = 0.0;
  }

  // Case where there is more snowfall than snowmelt: 
  // snowpack mixed, snowmelt has snowfall signature
  if(h > RNDOFFERR and snow_in > dh){

    h_eff = snow_in - dh;

    if(ctrl.sw_2H){
      // Snowpack: last (same timestep) in, first melt
      _d2Hsnowpack->matrix[r][c] = InputMix(h, _d2Hsnowpack->matrix[r][c],
					    h_eff, atm.getd2Hprecip()->matrix[r][c]);
      // Snowmelt: snowfall (=rain) signature
      _d2Hsnowmelt->matrix[r][c] = atm.getd2Hprecip()->matrix[r][c];
    }

    if(ctrl.sw_18O){
      // Snowpack: last (same timestep) in, first melt
      _d18Osnowpack->matrix[r][c] = InputMix(h, _d18Osnowpack->matrix[r][c],
					     h_eff, atm.getd18Oprecip()->matrix[r][c]);
      // Snowmelt: snowfall (=rain) signature
      _d18Osnowmelt->matrix[r][c] = atm.getd18Oprecip()->matrix[r][c];
    }

    if(ctrl.sw_Age){
      // Snowpack: last (same timestep) in, first melt
      _Agesnowpack->matrix[r][c] = InputMix(h, _Agesnowpack->matrix[r][c], h_eff, 0.0);
      // Snowmelt: age 0
      _Agesnowmelt->matrix[r][c] = 0.0;
    }

  } else {

    // Case where there is more snowmelt than snowfall: 
    // no mixing in snowpack, snowmelt has mixed signature
    h_eff = dh - snow_in;

    if(ctrl.sw_2H){
      // Snowpack: no change (all recent snow has melted)
      // Snowmelt: mix of snowpack and throughfall
      _d2Hsnowmelt->matrix[r][c] = InputMix(h_eff, _d2Hsnowpack->matrix[r][c],
					    snow_in, atm.getd2Hprecip()->matrix[r][c]);
    }    
    if(ctrl.sw_18O){
      // Snowpack: no change (all recent snow has melted)
      // Snowmelt: mix of snowpack and throughfall
      _d18Osnowmelt->matrix[r][c] = InputMix(h_eff, _d18Osnowpack->matrix[r][c],
					     snow_in, atm.getd18Oprecip()->matrix[r][c]);
    }
    
    if(ctrl.sw_Age){
      // Snowpack: no change (all snow in melted)
      // Snowmelt: mix of snowpack and throughfall
      _Agesnowmelt->matrix[r][c] = InputMix(h_eff, _Agesnowpack->matrix[r][c],
					    snow_in, 0.0);
      
    }
  }
  
}


