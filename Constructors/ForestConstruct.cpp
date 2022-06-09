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
 * ForestConstruct.cpp
 *
 *  Created on: May 5, 2010
 *      Author: Marco.Maneta
 */

#include "Forest.h"

Forest::Forest(Control &ctrl)
{
	try{

	//Read the base map and writes the dimensions of the grid
		_patches = new grid(ctrl.path_BasinFolder + ctrl.fn_patches, ctrl.MapType);
		_NRows = _patches->r;
		_NCols = _patches->c;
		_dx = _patches->dx;
		_Nsp = ctrl.NumSpecs + 1; //num of species plus bare soil

		/*sorts the basin with data cells according
		 * to the ldd after _DEM and _ldd have been created*/
		_vSortedGrid = SortGrid();

		_species = new Grove[_Nsp]; //creates Grove array with default constructor

		for (UINT4 i = 0; i < _Nsp; i++){ //initializes the grids in each Grove object
			_species[i].CreateGrids(_patches);
			// Tracking: vegetation-dependent maps
			if(ctrl.sw_trck && ctrl.sw_2H)
				_species[i].CreateGridsd2H(_patches);
			if(ctrl.sw_trck && ctrl.sw_18O)
				_species[i].CreateGridsd18O(_patches);
			if(ctrl.sw_trck && ctrl.sw_Age)
				_species[i].CreateGridsAge(_patches);
		}

		SetSpeciesParameters(ctrl);

		if(!ctrl.ForestStateVarsInputType.compare("tables"))
			SetStateVarsTabs(ctrl);
		else if(!ctrl.ForestStateVarsInputType.compare("maps"))
			SetStateVarsMaps(ctrl);
		else{
			cerr << "Illegal type " << ctrl.ForestStateVarsInputType << endl;
			throw;
		}


		//SetSpeciesParameters(ctrl);
        checkForestDatabase(); //check the sanity of the database

	}catch (std::bad_alloc &)
	  { cerr << "Cleaning up the forest..." << "\n";
		if(_patches)
			delete _patches;
		if(_species)
			delete[] _species;

	  throw;
	  }

}
