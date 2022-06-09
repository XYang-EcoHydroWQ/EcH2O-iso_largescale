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
 * ExtraGWDynamics.cpp
 *
 *  Created on: May 11, 2020
 *      Author: Xiaoqiang Yang
 */

#define ARMA_NO_DEBUG //disables armadillo bound checks for speed optimization
//#include <armadillo>
#include"Basin.h"

//using namespace arma;

void Basin::ExtraGWDynamics(Control &ctrl, double extragw,
			  double ex_qc, double ex_hj1i1,
		      double dt,int r, int c) {

	REAL8 leak = 0; //lekage

	REAL8 ex_alpha = 0;
    REAL8 ex_qj1i = 0;
    REAL8 ex_R = 0;
	REAL8 dtdx = 0;

	
	dtdx = dt/_dx;
	ex_qc = 0;   //[m2/s]
	ex_hj1i1 = 0;


    //get the global variable  
	leak = _BedrockLeakageFlux->matrix[r][c]; //stored as m/s in in ln 155 SoilWaterRedistribution.cpp
    leak = leak * dt;
	
	//update storage
	extragw += leak;
	//baseflow generation if it is channel cell
	if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) { 
	//if this is a channel cell and channels are activated
	//maxR = ( _porosity->matrix[r][c] - fc ) * soildepth; 
	//calculates the maximum gravitational water that can go
	  ex_qc = _KsatL3->matrix[r][c] * extragw * _chExGWparam->matrix[r][c];
	  extragw -= ex_qc * dtdx;
    }
	
	//terrestrial routing factor, the same as original subsurface flow "alpha"
	ex_alpha = _KsatL3->matrix[r][c] * sin(atan(_slope->matrix[r][c]));	
	
	//input from up-slope grid cell
	ex_qj1i = _ExtraGWupstreamBC->matrix[r][c]; //[m2/s]
	// DeepGW head
	ex_R = extragw;
	//
	ex_hj1i1 = (dtdx*ex_qj1i + ex_R)/(1+ex_alpha*dtdx); //[m]
	
	
}

