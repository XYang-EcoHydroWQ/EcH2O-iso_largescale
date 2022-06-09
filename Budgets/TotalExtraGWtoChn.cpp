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
 *    Marco Maneta, Sylvain Kuppl
 *******************************************************************************/
/*
 * TotalExtraGWtoChn.cpp
 *
 *  Created on: May 14, 2020
 *      Author: Xiaoqiang Yang
 *  Following TotalGWtoChn.cpp
 */

#include "Budget.h"

void Budget::TotalExtraGWtoChn(const grid* map1, const grid* map2, const Basin *b)
{
  ex_gwtochn += AccountStorages(map1, map2, b);
  // AccountStorages is used because FluxextraGWtoChn is already in m/tstep
  // (no need to multiply by dt)
}
void Budget::TotalExtraGWtoChn_d2H(const grid* map1, const grid* map2, const grid* map3, const Basin *b)
{
  ex_gwtochn_d2H += AccountTrckFluxes(map1, map2, map3, b);

}

void Budget::TotalExtraGWtoChn_d18O(const grid* map1, const grid* map2, const grid* map3, const Basin *b)
{
  ex_gwtochn_d18O += AccountTrckFluxes(map1, map2, map3, b);

}

// the water that already left is kept in the balance and "aging" as well
void Budget::TotalExtraGWtoChn_Age(const grid* map1, const grid* map2, const grid* map3, const Basin *b)
{
  ex_gwtochn_Age += ex_gwtochn * dt / 86400 + AccountTrckFluxes(map1, map2, map3, b);

}


// Instantaneous age reporting
void Budget::InstExtraGWtoChn_Age(const grid *map1, const grid *map2, const Basin *b)
{
  ex_Agegwtochn = AccountTrckFluxes2(map1, map2, b);
}
