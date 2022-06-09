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
 * Report2screen.cpp
 *
 *  Created on: Aug 2, 2010
 *      Author: Marco.Maneta
 */

#include "Sativa.h"

int Report2Screen(){

  /*
  // A few tracking reports here, so that 
  // the age mass balance check uses beg-of-time-step values, more simple!
  if(oControl->sw_trck){
    if(oControl->sw_2H){
      if(oControl->Rep_d2HsoilUp || oControl->RepTs_d2HsoilUp)
	oTracking->Calcd2Hsoil_12(*oBasin);
      if(oControl->Rep_d2HsoilAv || oControl->RepTs_d2HsoilAv)
	oTracking->Calcd2Hsoil_Av(*oBasin);
    }
    if(oControl->sw_18O){
      if(oControl->Rep_d18OsoilUp || oControl->RepTs_d18OsoilUp)
	oTracking->Calcd18Osoil_12(*oBasin);
      if(oControl->Rep_d18OsoilAv || oControl->RepTs_d18OsoilAv)
	oTracking->Calcd18Osoil_Av(*oBasin);
    }
    
    if(oControl->sw_Age){
      // Increment age by one time step duration
      oTracking->IncrementAge(*oBasin, *oControl);
      // Reported quantities
      if(oControl->Rep_AgesoilUp || oControl->RepTs_AgesoilUp)
	oTracking->CalcAgesoil_12(*oBasin);
      if(oControl->Rep_AgesoilAv || oControl->RepTs_AgesoilAv)
	oTracking->CalcAgesoil_Av(*oBasin);
    }
    }*/

  // ==== BasinSummary.txt --------------------------------------------------------
  // -----------------------------------------------------
  printf("\nTotal Precipitation (m3):  %.2f \t", oBudget->precipitation);
  ofSummary << oBudget->precipitation << "\t";

  printf("SWE (m3): %.2f \n", oBudget->snowpack);
  ofSummary << oBudget->snowpack << "\t";

  printf("Canopy Storage (m3): %.2f \t", oBudget->canopy);
  ofSummary << oBudget->canopy << "\t";

  printf("Ponding (m3): %.2f \n", oBudget->ponding);
  ofSummary << oBudget->ponding << "\t";

  printf("Soil water (m3): %.2f \t", oBudget->vadose);
  ofSummary << oBudget->vadose << "\t";

  printf("of which Groundwater (m3): %.2f \n", oBudget->grndwater);
  ofSummary << oBudget->grndwater << "\t";

  printf("Total Evapotranspiration (m3): %.2f \t", oBudget->evaporation);
  ofSummary << oBudget->evaporation << "\t";

  printf("Total Soil Evaporation (m3): %.2f \n", oBudget->evaporationS);
  ofSummary << oBudget->evaporationS << "\t";

  printf("Total Canopy Evaporation (m3): %.2f \t", oBudget->evaporationI);
  ofSummary << oBudget->evaporationI << "\t";

  printf("Total Transpiration (m3): %.2f \n", oBudget->transpiration);
  ofSummary << oBudget->transpiration << "\t";

  printf("Bedrock Leak (m3): %.2f \n", oBudget->leakage);
  ofSummary << oBudget->leakage << "\t";

  printf("Total OvlndFlow output (m3): %.2f \t", oBudget->ovlndflow);
  ofSummary << oBudget->ovlndflow << "\t";

  printf("Total GWFlow output (m3): %.2f \n", oBudget->gwtrflow);
  ofSummary << oBudget->gwtrflow << "\t";

  printf("Run-off to channel (m3): %.2f \t", oBudget->srftochn);
  ofSummary << oBudget->srftochn << "\t";

  printf("GW to channel (m3): %.2f \n", oBudget->gwtochn);
  ofSummary << oBudget->gwtochn << "\t";

  printf("Recharge (to layer 3, m3): %.2f \t", oBudget->recharge);
  ofSummary << oBudget->recharge << "\t";

  // Saturated area (% of the catchment)
  printf("Saturated area fraction: %.2f \n", oBudget->satarea);
  ofSummary << oBudget->satarea << "\t";

  printf("Mass Balance Error (%): %e \n", oBudget->MBErr);
  ofSummary << oBudget->MBErr ;//<< "\t";

  if(oControl->sw_trck and oControl->sw_2H){
    //printf("Deuterium Mass Balance Error (%): %e \n", oBudget->MBErr_d2H);
    ofSummary << "\t" << oBudget->MBErr_d2H ;
  }

  if(oControl->sw_trck and oControl->sw_18O){
    //printf("Oxygen 18 Mass Balance Error (%): %e \n", oBudget->MBErr_d18O);
    ofSummary << "\t" << oBudget->MBErr_d18O ;
  }

  if(oControl->sw_trck and oControl->sw_Age){
    //printf("Age Mass Balance Error (%): %e \n", oBudget->MBErr_Age);
    ofSummary << "\t" << oBudget->MBErr_Age ;
  }

  ofSummary << "\n";

  // ==== BasinAgeSummary.txt --------------------------------------------------------
  // -----------------------------------------------------
  if(oControl->sw_trck and oControl->sw_Age){
    ofAgeSummary << oBudget->AgeTot << "\t";
    ofAgeSummary << oBudget->Agesnowpack << "\t";
    ofAgeSummary << oBudget->Agecanopy << "\t";
    ofAgeSummary << oBudget->Ageponding << "\t";
    ofAgeSummary << oBudget->Agevadose << "\t";
    ofAgeSummary << oBudget->AgesoilL1 << "\t";
    ofAgeSummary << oBudget->AgesoilL2 << "\t";
    ofAgeSummary << oBudget->AgesoilL3 << "\t";
    ofAgeSummary << oBudget->Agegrndwater << "\t";
    ofAgeSummary << oBudget->AgeET << "\t";
    ofAgeSummary << oBudget->AgeevapS << "\t";
    ofAgeSummary << oBudget->AgeevapI << "\t";
    ofAgeSummary << oBudget->AgeevapT << "\t";
    ofAgeSummary << oBudget->Ageleakage << "\t";
    ofAgeSummary << oBudget->AgeOvlndOut << "\t";
    ofAgeSummary << oBudget->AgeGWOut << "\t";
    ofAgeSummary << oBudget->AgeOut << "\t";
    ofAgeSummary << oBudget->Agesrftochn << "\t";
    ofAgeSummary << oBudget->Agegwtochn << "\t";   
    ofAgeSummary << oBudget->Agerecharge << "\n";

  }

  return EXIT_SUCCESS;
}
