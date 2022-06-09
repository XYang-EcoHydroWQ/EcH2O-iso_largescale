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
 * PrefluxTrck_plusExtraGW.cpp
 *
 *  Created on: May 14, 2020
 *      Author: Xiaoqiang Yang
 *  This subroutine is created for backing up tracer states/fluxes 
 *  for previous time step (used for updating Extra GW outgoing signatures)
 */

#include "Basin.h"

void Tracking::PrefluxTrck_plusExtraGW(Control &ctrl)
{
  // Deuterium
  if(ctrl.sw_2H){
    *_Fd2HLattoExtraGW_old = *_Fd2HLattoExtraGW;
	*_d2Hleakage_old = *_d2Hleakage;
  }
  // Oxygen 18
  if(ctrl.sw_18O){
    *_Fd18OLattoExtraGW_old = *_Fd18OLattoExtraGW;	
	*_d18Oleakage_old = *_d18Oleakage;
  }
  
  // Water age
  if(ctrl.sw_Age){
    *_FAgeLattoExtraGW_old = *_FAgeLattoExtraGW;
	*_Ageleakage_old = *_Ageleakage;
  }
  
}


