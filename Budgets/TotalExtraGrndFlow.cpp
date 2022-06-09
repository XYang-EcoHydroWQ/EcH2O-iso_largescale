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
 *    Marco Maneta
 *******************************************************************************/
/*
 * TotalExtraGrndFlow.cpp
 *
 *  Created on: May 14, 2020
 *      Author: Xiaoqiang Yang
 *  Following TotalGrndFlow.cpp
 */

#include "Budget.h"

void Budget::TotalExtraGrndFlow(const vectCells *timeseries, const Basin *b)
{
	ex_gwtrflow += AccountFluxes(timeseries, b); //ex_gwtrflow should also be defined in Budget.h!!
}

void Budget::TotalExtraGrndFlow_d2H(const vectCells* timeseries1, const vectCells* timeseries2)
{
  ex_gwtrflow_d2H += AccountTrckFluxes(timeseries1, timeseries2);
  //gwtrflow_d2H = AccountTrckFluxes(timeseries1, timeseries2);
}

void Budget::TotalExtraGrndFlow_d18O(const vectCells* timeseries1, const vectCells* timeseries2)
{
  ex_gwtrflow_d18O += AccountTrckFluxes(timeseries1, timeseries2);
  //gwtrflow_d18O = AccountTrckFluxes(timeseries1, timeseries2);
}

// the water that already left is kept in the balance and "aging" as well
void Budget::TotalExtraGrndFlow_Age(const vectCells* timeseries1, const vectCells* timeseries2)
{
  ex_gwtrflow_Age += ex_gwtrflow * dt / 86400 + AccountTrckFluxes(timeseries1, timeseries2);  //?don't understand...
  //gwtrflow_Age = AccountTrckFluxes(timeseries1, timeseries2);
}

// Instantaneous age reporting
void Budget::InstExtraGrndFlow_Age(const vectCells* timeseries1, const vectCells* timeseries2)
{
  ex_AgeGWOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
