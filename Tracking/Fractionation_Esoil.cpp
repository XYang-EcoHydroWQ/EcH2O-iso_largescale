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
 * Fractionation_Esoil.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: Sylvain Kuppel
 */

#include "Tracking.h"

int Tracking::Frac_Esoil(Atmosphere &atm, Basin &bsn, Control &ctrl,
			 REAL8 V_old, REAL8 V_new, REAL8 &beta,
			 REAL8 &di_old, REAL8 &di_new, REAL8 &di_evap,
			 REAL8 &Ts, int r, int c, int iso){
  
  REAL8 Ta = atm.getTemperature()->matrix[r][c] + 273.15 ; // Air temperature (K)
  REAL8 ha = atm.getRelativeHumidty()->matrix[r][c]; // Atmospheric relative humidity (fraction)
  REAL8 es_s = SatVaporPressure(Ts-273.15) ; //saturated vapor pressure in the atmosphere (Pa)
  REAL8 ea_s = SatVaporPressure(Ta-273.15) ; // saturated vapor pressure at the surface (Pa)
  
  REAL8 ha_p; // Corrected relatvie air humidity above the surface (fraction)
  REAL8 hs; // Soil vapor saturation at the surface (fraction)
  REAL8 f; // Water loss fraction after evaporation (fraction)
  REAL8 alpha_p = 0; // equilibrium isotope fractionation factor (fraction)
  REAL8 eps_p; // equilibrium isotope fractionation factor (per mil)
  REAL8 eps_k = 0; // kinetic isotope fractionation factor (per mil)
  REAL8 eps; // total isotope fractionation factor (per mil)
  REAL8 di_atm = 0; // Isotopic signatures (permil)
  REAL8 di_s; // Limiting isotopic composition (per mil)
  REAL8 m; // Calculation factor (-)
  REAL8 n; // Parameter translating dominant water transport mode (-)
  REAL8 th_r, th_s, psiae, bclambda;
  REAL8 d1 = bsn.getSoilDepth1()->matrix[r][c];
  
  // Corrected relative humidity at the surface
  ha_p = ha*ea_s/es_s;
  
  // Relative humidity in the soil
  if(ctrl.toggle_hs == 0) {
    // Just 1
    hs = 1;
  } else if (ctrl.toggle_hs == 1) {
    // Corrected h for surface temperature and beta
    hs = beta + (1-beta)*ha_p;
  } else if (ctrl.toggle_hs == 2) {
    // Sorderberg et al. (2012), orginially for dry soils 
    // (departs from EcH2O's evaporation conceptualization of evap!
    th_r = bsn.getSoilMoistR()->matrix[r][c];
    th_s = bsn.getPorosityL1()->matrix[r][c];
    psiae = bsn.getPsiAE()->matrix[r][c];
    bclambda = bsn.getBClambda()->matrix[r][c];
    hs = expl(psiae*powl((V_old/d1-th_r)/(th_s-th_r),-bclambda)*18.0145/(8.3145*Ts*1000)); // waterpot*molarmass/(R*Ts*waterdensity)
  } else {
    std::cout << "Wrong option in the surface humidity toggle switch for fractionation. Please verify the configuration file" << std::endl;
    exit(EXIT_FAILURE);  
  }
  
  // Horita and Wesolowski (1994)
  if(iso == 0) // deuterium
    alpha_p = expl((1158.8*powl(Ta,3)*1e-9 - 1620.1*powl(Ta,2)*1e-6 + 794.84*Ta*1e-3 - \
		    161.04 + 2.9992*1e9/powl(Ta,3))/1000);
  else if (iso == 1) // oxygen 18
    alpha_p = expl((-7.685 + 6.7123*1000/Ta - 1.6664*1e6/powl(Ta,2) +
		    0.35041*1e9/powl(Ta,3))/1000);
  
  // Skrzypek et al. (2015)
  eps_p = (1 - 1/alpha_p)*1000;
  
  // (Gat, 1995) + (Gibson and Reid, 2014)
  if(iso == 0) // deuterium
    di_atm = (atm.getd2Hprecip()->matrix[r][c] - eps_p)/ alpha_p;
  else if(iso ==1) // oxygen 18
    di_atm = (atm.getd18Oprecip()->matrix[r][c] - eps_p)/ alpha_p;

  // Water transport mode: from diffusive (=1, dry soil) to turbulent (=0.5, water body)
  if(ctrl.toggle_n == 0) {
    n = 1;
  } else if (ctrl.toggle_n == 1) {
    // Mathieu and Bariac (1996) + Braud et al. (2005)
    n = 1 - 0.5*(V_old/d1-bsn.getSoilMoistR()->matrix[r][c])/
      (bsn.getPorosityL1()->matrix[r][c]-bsn.getSoilMoistR()->matrix[r][c]);
  } else {
    std::cout << "Wrong option in the soil fractionation n toggle switch. Please verify the configuration file" << std::endl;
    exit(EXIT_FAILURE);  
  }
  
  // Kinetic fractionation factor epsilon_k
  if(ctrl.toggle_ek == 0) {
    // Value of Di/D from Merlivat (1978)
    if(iso==0) 
      eps_k = (hs-ha_p)*(1-0.9757)*1000*n;
    else if(iso==1) 
      eps_k = (hs-ha_p)*(1-0.9727)*1000*n;

  } else if (ctrl.toggle_ek == 1) {
    // Value of Di/D from Vogt (1976)
    if(iso==0) // deuterium
      eps_k = (hs-ha_p)*(1-0.9877)*1000*n; 
    else if(iso==1) // oxygen 18
      eps_k = (hs-ha_p)*(1-0.9859)*1000*n; 

  } else if (ctrl.toggle_ek == 2) {
    REAL8 u = atm.getWindSpeed()->matrix[r][c]; // wind speed (m.s-1)
    // From Merlivat and Jouzel (1979) adapted by Haese et al. (2013)
    if(iso==0) // deuterium
      eps_k = u > 7.0 ? (hs-ha_p)*0.88*(0.285*u+0.82) : (hs-ha_p)*0.88*6;
    else if (iso==1) // oxygen 18
      eps_k = u > 7.0 ? (hs-ha_p)*(0.285*u+0.82) : (hs-ha_p)*6;
  }

  // --- Generalized following Good et al. (2014) -------------
  // Gibson and Reid (2010)
  eps = hs*eps_p + eps_k;
  // (Gat and Levy, 1978) + (Gat, 1981)
  di_s = (ha_p*di_atm + eps) / (ha_p - eps/1000);
  // (Welhan and Fritz, 1977) + (Allison and Leaney, 1982)
  m = (ha_p - eps/1000) / (hs - ha_p + eps_k/1000);
  // ----------------------------------------------------------
  
  // Evaporative loss fraction
  f = V_new/V_old;
  
  // (Hamilton et al., 2005)
  // New isotopic signature in topsoil
  di_new = di_s - (di_s - di_old) * powl(f,m);
  
  // Isotopic signature of evaporated water
  di_evap = std::max<double>(-1000,(hs*alpha_p*di_new - ha_p*di_atm - eps)/ (hs - ha_p + eps_k/1000));

  // if(abs(di_new)>100 or abs(di_evap)>1e3)
  //    cout << r << " " << c << "| iso:" << iso << "| evapS:" << V_old-V_new << "| disoil_new:" << di_new << "|di_old:" << di_old <<
  //      "| di_star:" << di_s << "| di_evap:" << di_evap << 
  //      "| f:" << f <<"| m:" << m <<"| ha_p:" << ha_p <<"| hs:" << hs << endl;
  
  return EXIT_SUCCESS;
  
}
