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
 * CalcRootDistrib.cpp
 *
 *  Created on: Feb 6th, 2017
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"


int Basin::CalcRootDistrib(){

  UINT4 r, c;
  UINT4 s, nsp;
  REAL8 frac1, frac2;
  REAL8 k, d1, d2, d;
  
  nsp = fForest->getNumSpecies();

#pragma omp parallel default(none)		\
  private(s,r,c,k,d,d1,d2,frac1,frac2) \
  shared(nsp)
  { 
#pragma omp for nowait

    for (UINT4 j = 0; j < _vSortedGrid.cells.size(); j++)
      {
	r = _vSortedGrid.cells[j].row;
	c = _vSortedGrid.cells[j].col;

	d = _soildepth->matrix[r][c];
	d1 = _depth_layer1->matrix[r][c];
	d2 = _depth_layer2->matrix[r][c];
	
	for (s = 0; s < nsp; s++) {

	  if (s == nsp - 1) { //if this is bare ground set fracs to 0
	    frac1 = 0;
	    frac2 = 0;
	  } 
	  else {
	    // use exponential profile
	    k = fForest->getKRoot(s);
	    frac1 = (1 - expl(-k*d1))/(1-expl(-k*d));
	    frac2 = (expl(-k*d1) - expl(-k*(d1+d2)))/(1-expl(-k*d));
	  }

	  fForest->setRootFrac1Species(s, r, c, frac1);
	  fForest->setRootFrac2Species(s, r, c, frac2);
	} 
	/*
	k = _Kroot->matrix[r][c];
	_rootfrac1->matrix[r][c] = (1 - expl(-k*d1))/(1-expl(-k*d));
	_rootfrac2->matrix[r][c] = (expl(-k*d1) - expl(-k*(d1+d2)))/(1-expl(-k*d));
	*/
      } // for
  } //end omp parallel block
  return EXIT_SUCCESS;

}
