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
 * CalcPropLayers.cpp
 *
 *  Created on: Jun 22nd, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"


// Layer-integrated saturated hydraulic conductivities based on exponentially-decreasing profile
int Basin::CalcKsatLayers(Control &ctrl){

  UINT4 r, c;
  REAL8 k, K0, d1, d2, d;
  
#pragma omp parallel default(none)		\
  private(r,c,k,K0,d,d1,d2) shared(ctrl)			
  { 
#pragma omp for nowait

    for (UINT4 j = 0; j < _vSortedGrid.cells.size(); j++)
      {
	r = _vSortedGrid.cells[j].row;
	c = _vSortedGrid.cells[j].col;

	K0 = _Ksat0->matrix[r][c];

	// Check if profile option is activated: exponential profile
	if(ctrl.sw_expKsat){

	  d = _soildepth->matrix[r][c];
	  d1 = _depth_layer1->matrix[r][c];
	  d2 = _depth_layer2->matrix[r][c];
	  k = _kKsat->matrix[r][c];
	  
	  _KsatL1->matrix[r][c] = k*K0 * (1 - expl(-d1/k))/ d1;
	  _KsatL2->matrix[r][c] = k*K0 * (expl(-d1/k) - expl(-(d1+d2)/k)) / d2;
	  _KsatL3->matrix[r][c] = k*K0 * (expl(-(d1+d2)/k) - expl(-d/k)) / (d-d1-d2);

	} else {
	  _KsatL1->matrix[r][c] = K0 ;
	  _KsatL2->matrix[r][c] = K0 ;
	  _KsatL3->matrix[r][c] = K0 ;
	}
	
      } // for
  } //end omp parallel block
  return EXIT_SUCCESS;

}

// Layer-integrated porosities based on exponentially-decreasing profile
int Basin::CalcPorosLayers(Control &ctrl){

  UINT4 r, c;
  REAL8 k, phi0, d1, d2, d;
  
#pragma omp parallel default(none)		\
  private(r,c,k,phi0,d,d1,d2) shared(ctrl)			
  { 
#pragma omp for nowait

    for (UINT4 j = 0; j < _vSortedGrid.cells.size(); j++)
      {
	r = _vSortedGrid.cells[j].row;
	c = _vSortedGrid.cells[j].col;

	phi0 = _porosity0->matrix[r][c];

	// Check if profile option is activated: exponential profile
	if(ctrl.sw_expPoros){

	  d = _soildepth->matrix[r][c];
	  d1 = _depth_layer1->matrix[r][c];
	  d2 = _depth_layer2->matrix[r][c];
	  k = _kporos->matrix[r][c];
	  
	  _porosityL1->matrix[r][c] = k*phi0 * (1 - expl(-d1/k)) / d1;
	  _porosityL2->matrix[r][c] = k*phi0 * (expl(-d1/k) - expl(-(d1+d2)/k)) / d2;
	  _porosityL3->matrix[r][c] = k*phi0 * (expl(-(d1+d2)/k) - expl(-d/k)) / (d-d1-d2);

	} else {
	  _porosityL1->matrix[r][c] = phi0 ;
	  _porosityL2->matrix[r][c] = phi0 ;
	  _porosityL3->matrix[r][c] = phi0 ;
	}
	
      } // for
  } //end omp parallel block
  return EXIT_SUCCESS;

}

// Layer-dependent field capacity
int Basin::CalcFieldCapacity(){

	UINT4 r, c;
	UINT4 length = _vSortedGrid.cells.size();

#pragma omp parallel default(none)\
  private(  r,c), shared(length)
	for (UINT4 j = 0; j < length ; j++)
	  {
	    r = _vSortedGrid.cells[j].row;
	    c = _vSortedGrid.cells[j].col;
	    
	    _fieldcapL1->matrix[r][c] =
	      powl(_psi_ae->matrix[r][c] / 3.36 ,1/_BClambda->matrix[r][c])
	      * (_porosityL1->matrix[r][c] - _theta_r->matrix[r][c])
	      + _theta_r->matrix[r][c];

	    _fieldcapL2->matrix[r][c] =
	      powl(_psi_ae->matrix[r][c] / 3.36 ,1/_BClambda->matrix[r][c])
	      * (_porosityL2->matrix[r][c] - _theta_r->matrix[r][c])
	      + _theta_r->matrix[r][c];

	    _fieldcapL3->matrix[r][c] =
	      powl(_psi_ae->matrix[r][c] / 3.36 ,1/_BClambda->matrix[r][c])
	      * (_porosityL3->matrix[r][c] - _theta_r->matrix[r][c])
	      + _theta_r->matrix[r][c];

	  }
		
	return EXIT_SUCCESS;
}
