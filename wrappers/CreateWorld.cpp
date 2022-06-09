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
 * CreateWorld.cpp
 *
 *  Created on: Jul 30, 2010
 *      Author: Marco.Maneta
 */

#include "Sativa.h"

int CreateWorld(char* argv[]){

  oControl = new Control;
  cout << "Control created ok... " << "\n";
  oControl->ReadConfigFile(argv[1]);
  cout << "Config.ini read ok... " << "\n";

  oBasin = new Basin(*oControl);
  cout << "Basin created ok... " << "\n";

  oAtmosphere = new Atmosphere(*oControl);
  cout << "Atmosphere created ok... " << "\n";

  oReport = new Report(*oControl);
  cout << "Report created ok... " << "\n";

  oTracking = new Tracking(*oControl, *oBasin);
  if(oControl->sw_trck)
    cout << "Isotope module created ok... " << "\n";

  oBudget = new Budget(oBasin, oControl, oTracking);
  cout << "Budget created ok... " << "\n";

  // == Basin Summary ==========================================================
  // ---
  try{
    ofSummary.open((oControl->path_ResultsFolder + "BasinSummary.txt").c_str());
    if(!ofSummary)
      throw std::ios::failure("Error opening BasinSummary.txt buffer\n");

  }catch(const std::exception &e){
    cout << e.what() << endl;
    throw;
  }
  // Headers for BasinSummary
  ofSummary << "Precip\t";
  ofSummary << "SWE\t";
  ofSummary << "Intrcp\t";
  ofSummary << "Surface\t";
  ofSummary << "SoilW\t";
  ofSummary << "GW\t";
  ofSummary << "ET\t";
  ofSummary << "EvapS\t";
  ofSummary << "EvapI\t";
  ofSummary << "EvapT\t";
  ofSummary << "Leakage\t";
  ofSummary << "SrfOut\t";
  ofSummary << "GWOut\t";
  ofSummary << "SrftoCh\t";
  ofSummary << "GWtoCh\t";
  ofSummary << "Rchrge\t";
  ofSummary << "SatExt\t";
  ofSummary << "MBErr";
  if(oControl->sw_trck and oControl->sw_2H)
    ofSummary << "\t2H_MBE";
  if(oControl->sw_trck and oControl->sw_18O)
    ofSummary << "\t18O_MBE";
  if(oControl->sw_trck and oControl->sw_Age)
    ofSummary << "\tAge_MBE";
  ofSummary << "\n";

  // == Age Summary ==========================================================
  // ---  
  if(oControl->sw_trck and oControl->sw_Age){
    try{
      ofAgeSummary.open((oControl->path_ResultsFolder + "BasinAgeSummary.txt").c_str());
      if(!ofAgeSummary)
	throw std::ios::failure("Error opening BasinAgeSummary.txt buffer\n");
      
    }catch(const std::exception &e){
      cout << e.what() << endl;
      throw;
    }
    // Headers for BasinSummary
    ofAgeSummary << "StorAll\t";
    ofAgeSummary << "Snow\t";
    ofAgeSummary << "Intercp\t";
    ofAgeSummary << "Surface\t";
    ofAgeSummary << "Soil\t";
    ofAgeSummary << "SoilL1\t";
    ofAgeSummary << "SoilL2\t";
    ofAgeSummary << "SoilL3\t";
    ofAgeSummary << "GW\t";
    ofAgeSummary << "ET\t";
    ofAgeSummary << "EvapS\t";
    ofAgeSummary << "EvapI\t";
    ofAgeSummary << "EvapT\t";
    ofAgeSummary << "Leakage\t";
    ofAgeSummary << "SrfOut\t";
    ofAgeSummary << "GWOut\t";
    ofAgeSummary << "AllOut\t";
    ofAgeSummary << "SrftoCh\t";
    ofAgeSummary << "GWtoCh\t";
    ofAgeSummary << "Rchrge\n";

  }

  
  return EXIT_SUCCESS;
}
