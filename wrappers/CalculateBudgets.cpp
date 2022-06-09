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
 * CalculateBudgets.cpp
 *
 *  Created on: Aug 2, 2010
 *      Author: Marco.Maneta
 */

#include "Sativa.h"

int CalculateBudgets(){

  //REAL8 output;

  oBudget->TotalPrecipitation(oAtmosphere->getPrecipitation(), oBasin->getTTarea(), oAtmosphere);
  oBudget->TotalEvaporation(oBasin->getEvaporation(), oBasin->getTTarea(), oBasin);
  oBudget->TotalEvaporationS(oBasin->getEvaporationS_all(), oBasin->getTTarea(), oBasin);

  oBudget->TotalEvaporationI(oBasin->getEvaporationI_all(), oBasin->getTTarea(), oBasin);
  oBudget->TotalTranspiration(oBasin->getTranspiration_all(), oBasin->getTTarea(), oBasin);
  oBudget->TotalBedrockLeakage(oBasin->getBedrockLeakage(), oBasin->getTTarea(), oBasin);
  oBudget->TotalOvlndFlow(oBasin->getDailyOvlndOutput(), oBasin);
  oBudget->TotalGrndFlow(oBasin->getDailyGwtrOutput(), oBasin);
  oBudget->TotalStorage(oBasin->getCanopyStorage(),
			oBasin->getSnowWaterEquiv(),
            oBasin->getPondingWater(),
			oBasin->getChannelWater(),
			//oBasin->getSoilWaterDepth(),
			oBasin->getSoilWaterDepthL1(),
			oBasin->getSoilWaterDepthL2(),
			oBasin->getSoilWaterDepthL3(),
			//oBasin->getGravityWater(),
			oBasin->getGrndWater(),
            oBasin->getTTarea(),
			oBasin);
  oBudget->TotalSrftoChn(oBasin->getFluxSrftoChn(), oBasin->getTTarea(), oBasin);
  oBudget->TotalGWtoChn(oBasin->getFluxGWtoChn(), oBasin->getTTarea(), oBasin);
  oBudget->TotalRecharge(oBasin->getFluxRecharge(), oBasin->getTTarea(), oBasin);
  oBudget->TotalSaturationArea(oBasin->getSatArea(), oBasin->getTTarea(), oBasin);
  
  //budget for Extra GW yangx 2020-05
  if(oControl->sw_extraGW){
    oBudget->TotalExtraGrndFlow(oBasin->getDailyExtraGwtrOutput(), oBasin);
    oBudget->TotalExtraGWtoChn(oBasin->getFluxExtraGWtoChn(), oBasin->getTTarea(), oBasin);
  }
  // ---------------------------------------------------------------------------------------------

  // Tracking -------------------------------------------------------------------------------------
  if(oControl->sw_trck){
    // Deuterium
    if (oControl->sw_2H){
      oBudget->TotalPrecipitation_d2H(oAtmosphere->getPrecipitation(), oAtmosphere->getd2Hprecip(), 
				     oBasin->getTTarea(), oAtmosphere);
      oBudget->TotalEvaporationS_d2H(oBasin->getEvaporationS_all(), oTracking->getd2HevapS_sum(), 
				    oBasin->getTTarea(), oBasin);
      oBudget->TotalEvaporationC_d2H(oBasin->getChanEvap(), oTracking->getd2HevapC_sum(), 
                    oBasin->getTTarea(), oBasin);
      oBudget->TotalEvaporationI_d2H(oBasin->getEvaporationI_all(), oTracking->getd2HevapI_sum(), 
				    oBasin->getTTarea(), oBasin);
      oBudget->TotalTranspiration_d2H(oBasin->getTranspiration_all(), oTracking->getd2HevapT_sum(), 
				    oBasin->getTTarea(), oBasin);
      oBudget->TotalBedrockLeakage_d2H(oBasin->getBedrockLeakage(), oTracking->getd2Hleakage(),
                    oBasin->getTTarea(), oBasin);
      oBudget->TotalOvlndFlow_d2H(oBasin->getDailyOvlndOutput(), oTracking->getd2HOvlndOutput());
      oBudget->TotalGrndFlow_d2H(oBasin->getDailyGwtrOutput(), oTracking->getd2HGwtrOutput());
      oBudget->TotalStorage_d2H(oBasin->getCanopyStorage(), oTracking->getd2Hcanopy_sum(),
			       oBasin->getSnowWaterEquiv(),  oTracking->getd2Hsnowpack(),
			       oBasin->getPondingWater(), oTracking->getd2Hsurface(),
			       oBasin->getChannelWater(),oTracking->getd2Hchan(),
			       oBasin->getSoilWaterDepthL1(), oTracking->getd2Hsoil1(),
			       oBasin->getSoilWaterDepthL2(), oTracking->getd2Hsoil2(),
			       oBasin->getSoilWaterDepthL3(), oTracking->getd2Hsoil3(),
			       oBasin->getGrndWater(), oTracking->getd2Hgroundwater(),
                   oBasin->getTTarea(), 
			       oBasin);//, oControl);
      //budget for Extra GW yangx 2020-05
      if(oControl->sw_extraGW){
        oBudget->TotalExtraGrndFlow_d2H(oBasin->getDailyExtraGwtrOutput(), 
                                        oTracking->getd2HExtraGWtrOutput());
        //d2HExtraGWtoChn = d2HExtraGWtoLat 
        oBudget->TotalExtraGWtoChn_d2H(oBasin->getFluxExtraGWtoChn(), oTracking->getd2HExtraGWtoLat(),
                                       oBasin->getTTarea(), oBasin);
      }	
    }

    // Oxygen 18
    if (oControl->sw_18O){
      oBudget->TotalPrecipitation_d18O(oAtmosphere->getPrecipitation(), oAtmosphere->getd18Oprecip(), 
				       oBasin->getTTarea(), oAtmosphere);
      oBudget->TotalEvaporationS_d18O(oBasin->getEvaporationS_all(), oTracking->getd18OevapS_sum(), 
					oBasin->getTTarea(), oBasin);
      oBudget->TotalEvaporationC_d18O(oBasin->getChanEvap(), oTracking->getd18OevapC_sum(), 
					oBasin->getTTarea(), oBasin);
      oBudget->TotalEvaporationI_d18O(oBasin->getEvaporationI_all(), oTracking->getd18OevapI_sum(), 
				      oBasin->getTTarea(), oBasin);
      oBudget->TotalTranspiration_d18O(oBasin->getTranspiration_all(), oTracking->getd18OevapT_sum(), 
				       oBasin->getTTarea(), oBasin);
      oBudget->TotalBedrockLeakage_d18O(oBasin->getBedrockLeakage(), oTracking->getd18Oleakage(), 
					oBasin->getTTarea(), oBasin);
      oBudget->TotalOvlndFlow_d18O(oBasin->getDailyOvlndOutput(), oTracking->getd18OOvlndOutput());
      oBudget->TotalGrndFlow_d18O(oBasin->getDailyGwtrOutput(), oTracking->getd18OGwtrOutput());
      oBudget->TotalStorage_d18O(oBasin->getCanopyStorage(), oTracking->getd18Ocanopy_sum(),
				 oBasin->getSnowWaterEquiv(),  oTracking->getd18Osnowpack(),
				 oBasin->getPondingWater(), oTracking->getd18Osurface(),
			     oBasin->getChannelWater(),oTracking->getd18Ochan(),
				 oBasin->getSoilWaterDepthL1(), oTracking->getd18Osoil1(),
				 oBasin->getSoilWaterDepthL2(), oTracking->getd18Osoil2(),
				 oBasin->getSoilWaterDepthL3(), oTracking->getd18Osoil3(),
				 oBasin->getGrndWater(), oTracking->getd18Ogroundwater(),
                 oBasin->getTTarea(), 
				 oBasin);//, oControl);
	      //budget for Extra GW yangx 2020-05
      if(oControl->sw_extraGW){
        oBudget->TotalExtraGrndFlow_d18O(oBasin->getDailyExtraGwtrOutput(), 
                                         oTracking->getd18OExtraGWtrOutput());
        oBudget->TotalExtraGWtoChn_d2H(oBasin->getFluxExtraGWtoChn(), oTracking->getd18OExtraGWtoLat(),
                                       oBasin->getTTarea(), oBasin);
      }
    }
    // Age
    if (oControl->sw_Age){
     // Age for mass balance calculation
      oBudget->TotalPrecipitation_Age();
      //oAtmosphere->getPrecipitation(), oAtmosphere);//, oControl);
      //oBudget->precipitation_Age += oBudget
      oBudget->TotalEvaporationS_Age(oBasin->getEvaporationS_all(), oTracking->getAgeevapS_sum(), 
				     oBasin->getTTarea(), oBasin);
      oBudget->TotalEvaporationC_Age(oBasin->getChanEvap(), oTracking->getAgeevapC_sum(), 
				     oBasin->getTTarea(),oBasin);
      oBudget->TotalEvaporationI_Age(oBasin->getEvaporationI_all(), oTracking->getAgeevapI_sum(), 
				     oBasin->getTTarea(), oBasin);
      oBudget->TotalTranspiration_Age(oBasin->getTranspiration_all(), oTracking->getAgeevapT_sum(), 
				      oBasin->getTTarea(), oBasin);
      oBudget->TotalBedrockLeakage_Age(oBasin->getBedrockLeakage(), oTracking->getAgeleakage(), 
				       oBasin->getTTarea(), oBasin);
      oBudget->TotalOvlndFlow_Age(oBasin->getDailyOvlndOutput(), oTracking->getAgeOvlndOutput());
      oBudget->TotalGrndFlow_Age(oBasin->getDailyGwtrOutput(), oTracking->getAgeGwtrOutput());

      //output = oTracking->getAgesnowpack()->maxi();
      //std::cout << output << endl;
      //output = oTracking->getAgesnowpack()->mini();
      //std::cout << output << endl;
      oBudget->TotalStorage_Age(oBasin->getCanopyStorage(), oTracking->getAgecanopy_sum(),
				oBasin->getSnowWaterEquiv(),  oTracking->getAgesnowpack(),
				oBasin->getPondingWater(), oTracking->getAgesurface(),
                oBasin->getChannelWater(),oTracking->getAgechan(),
				oBasin->getSoilWaterDepthL1(), oTracking->getAgesoil1(),
				oBasin->getSoilWaterDepthL2(), oTracking->getAgesoil2(),
				oBasin->getSoilWaterDepthL3(), oTracking->getAgesoil3(),
				oBasin->getGrndWater(), oTracking->getAgegroundwater(),
                oBasin->getTTarea(), 
				oBasin);//, oControl);
				
	      //budget for Extra GW yangx 2020-05
      if(oControl->sw_extraGW){
        oBudget->TotalExtraGrndFlow_Age(oBasin->getDailyExtraGwtrOutput(),
                                        oTracking->getAgeExtraGWtrOutput());
        oBudget->TotalExtraGWtoChn_d2H(oBasin->getFluxExtraGWtoChn(), oTracking->getAgeExtraGWtoChn(),
                                       oBasin->getTTarea(), oBasin);
      }
      // Age for BasinAgeSummary.txt
      oBudget->InstEvaporation_Age(oBasin->getEvaporationS_all(), oTracking->getAgeevapS_sum(), 
				   oBasin->getEvaporationI_all(), oTracking->getAgeevapI_sum(), 
				   oBasin->getTranspiration_all(), oTracking->getAgeevapT_sum(), 
				   oBasin);
      oBudget->InstEvaporationS_Age(oBasin->getEvaporationS_all(), oTracking->getAgeevapS_sum(), 
				     oBasin);
      oBudget->InstEvaporationI_Age(oBasin->getEvaporationI_all(), oTracking->getAgeevapI_sum(), 
				     oBasin);
      oBudget->InstTranspiration_Age(oBasin->getTranspiration_all(), oTracking->getAgeevapT_sum(), 
				      oBasin);
      oBudget->InstBedrockLeakage_Age(oBasin->getBedrockLeakage(), oTracking->getAgeleakage(), 
				       oBasin);

      oBudget->InstOvlndFlow_Age(oBasin->getDailyOvlndOutput(), oTracking->getAgeOvlndOutput());
      oBudget->InstGrndFlow_Age(oBasin->getDailyGwtrOutput(), oTracking->getAgeGwtrOutput());

      oBudget->InstOut_Age(oBasin->getEvaporationS_all(), oTracking->getAgeevapS_sum(), 
			   oBasin->getEvaporationI_all(), oTracking->getAgeevapI_sum(), 
			   oBasin->getTranspiration_all(), oTracking->getAgeevapT_sum(), 
			   oBasin->getBedrockLeakage(), oTracking->getAgeleakage(), 
			   oBasin->getDailyOvlndOutput(), oTracking->getAgeOvlndOutput(),
			   oBasin->getDailyGwtrOutput(), oTracking->getAgeGwtrOutput(),
               oBasin->getTTarea(),
			   oBasin);

      oBudget->InstSrftoChn_Age(oBasin->getFluxSrftoChn(), oTracking->getAgeSrftoChn(), oBasin);
      oBudget->InstGWtoChn_Age(oBasin->getFluxGWtoChn(), oTracking->getAgeGWtoChn(), oBasin);
      oBudget->InstRecharge_Age(oBasin->getFluxRecharge(), oTracking->getAgeRecharge(), oBasin);     
      //Extra GW
      if(oControl->sw_extraGW){
		oBudget->InstExtraGrndFlow_Age(oBasin->getDailyExtraGwtrOutput(), oTracking->getAgeExtraGWtrOutput());
        oBudget->InstExtraGWtoChn_Age(oBasin->getFluxExtraGWtoChn(), oTracking->getAgeExtraGWtoChn(), oBasin);
      }
    }

  } // ---------------------------------------------------------------------------------------------
  // Mass balance check
  oBudget->TrckBalanceError(oControl);
  //oBudget->MassBalanceError();

  return EXIT_SUCCESS;
}
