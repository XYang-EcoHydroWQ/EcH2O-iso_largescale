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
 * OutletVals.cpp
 *
 *  Created on: Mar 1, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::OutletVals(Basin &bsn, Control &ctrl, int mode, 
                 double Qk1, double Outletpond, int r, int c)
{

  double trnew = 0;
    
  if (mode == 0){ // Initialization
    if(ctrl.sw_2H){
      _d2HOvlndOutput.cells.clear();
      _d2HGwtrOutput.cells.clear();
    }
    if(ctrl.sw_18O){
      _d18OOvlndOutput.cells.clear();
      _d18OGwtrOutput.cells.clear();
    }
    if(ctrl.sw_Age){
      _AgeOvlndOutput.cells.clear();
      _AgeGwtrOutput.cells.clear();
    }
  }

  else if (mode == 1){ // Active mode

    // Deuterium
    if(ctrl.sw_2H){
      _d2HGwtrOutput.cells.push_back(cell(r, c, _d2Hgroundwater->matrix[r][c]));
	  //yangx 2020-11 mixing of surface ponding and channel storage
      //_d2HOvlndOutput.cells.push_back(cell(r, c, _d2Hsurface->matrix[r][c]));
      //here just need the relative weights of channel and ponding waters, no need for unit transformation	  
	  trnew = (Qk1+Outletpond) > RNDOFFERR ? (Qk1*_d2Hchan->matrix[r][c] +
           	   Outletpond*_d2Hsurface->matrix[r][c])/(Qk1+Outletpond) : _d2Hsurface->matrix[r][c];
      _d2HOvlndOutput.cells.push_back(cell(r, c, trnew));
    }
    // Oxygen 18
    if(ctrl.sw_18O){
      _d18OGwtrOutput.cells.push_back(cell(r, c, _d18Ogroundwater->matrix[r][c]));
      //_d18OOvlndOutput.cells.push_back(cell(r, c, _d18Osurface->matrix[r][c])); 
	  trnew = (Qk1+Outletpond) > RNDOFFERR ? (Qk1*_d18Ochan->matrix[r][c] +
           	   Outletpond*_d18Osurface->matrix[r][c])/(Qk1+Outletpond) : _d18Osurface->matrix[r][c];
      _d18OOvlndOutput.cells.push_back(cell(r, c, trnew));
    }
    // Age
    if(ctrl.sw_Age){
      _AgeGwtrOutput.cells.push_back(cell(r, c, _Agegroundwater->matrix[r][c]));
      //_AgeOvlndOutput.cells.push_back(cell(r, c, _Agesurface->matrix[r][c])); 
	  trnew = (Qk1+Outletpond) > RNDOFFERR ? (Qk1*_Agechan->matrix[r][c] +
           	   Outletpond*_Agesurface->matrix[r][c])/(Qk1+Outletpond) : _Agesurface->matrix[r][c];
      _AgeOvlndOutput.cells.push_back(cell(r, c, trnew));
    }
  }
}


