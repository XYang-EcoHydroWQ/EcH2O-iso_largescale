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
 * SolveCanopyFluxes.cpp
 *
 *  Created on: Jul 9, 2010
 *      Author: Marco.Maneta
 */

#include"Basin.h"

int Basin::SolveCanopyFluxes(Atmosphere &atm, Control &ctrl, Tracking &trck) {
  
  UINT4 r, c;
  REAL8 dt = ctrl.dt;
  
  REAL8 Tp = 0;
  REAL8 maxTp = 0;
  REAL8 minTp = 0;
  REAL8 snow = 0; //amount of snow reaching the ground ms-1
  REAL8 rain = 0;//amount of rain reaching the ground ms-1
  REAL8 sno_rain_thres = 0; //temperature threshold for snow rain transition, degC

  REAL8 ra; //soil aerodynamic resistance


  REAL8 evap = 0; //evaporation for the tree groves
  REAL8 transp = 0; //transpiratin for the tree groves
  REAL8 netR = 0; //net radiation for the tree groves
  REAL8 evap_f = 0; //total evaporation for the entire cell
  REAL8 transp_f = 0; //total transpiration for the entire cell
  //REAL8 ETP;

  //canopy storage parameters
  REAL8 D = 0; //canopy trascolation (amount of water that actually reach the ground)
  REAL8 DelCanStor = 0; //Canopy Storage

  //soil parameters
  REAL8 rootdepth;
  REAL8 thetar;
  REAL8 fc3;
  REAL8 poros1, poros2, poros3;
  REAL8 psi_ae;
  REAL8 bclambda;

  REAL8 froot1;
  REAL8 froot2;
  REAL8 froot3;
  //soil layer depths
  REAL8 d1;
  REAL8 d2;
  REAL8 d3;

  //aerodynamic resistance parameters
  REAL8 za; //height of wind speed measurements
  REAL8 z0o; // roughness length
  REAL8 zdo; //zero plane displacement
  REAL8 wind; //wind speed
  REAL8 treeheight;

  REAL8 theta;// = 0;
  REAL8 theta2;// = 0;
  REAL8 theta3;// = 0;
  REAL8 theta_available=0; //water available to roots

  UINT4 nsp;
  REAL8 p; //fraction of species s
  REAL8 veg_p; //summed fraction of vegetation

  //unsigned int j;
  UINT4 s;
  int  thre=0;

  // Tracking
  REAL8 pTrp1, pTrp2, pTrp3;
  REAL8 d2HevapT_f, d18OevapT_f, AgeevapT_f;
  REAL8 d2HevapI_f, d18OevapI_f, AgeevapI_f;
  REAL8 d2Hcanopy_f, d18Ocanopy_f, Agecanopy_f;
  REAL8 kTB_L1, kTB_L2;
  REAL8 kMW_L1, kMW_L2;
  REAL8 dth1, dth2, dth3;
  REAL8 theta_MW1, theta_MW2;
  // Initialize to zero
  _Rn_sum->reset();
  if(ctrl.sw_trck){
    //For tracking
    _FluxCnptoSrf->reset(); // canopy/sky to surface
    _FluxCnptoSnow->reset(); // canopy/sky to snowpack
  }

#pragma omp parallel default(none)					\
  private( s, r,c, p,  treeheight, wind, za, z0o, zdo,			\
	   Tp, maxTp, minTp, snow, rain, sno_rain_thres, evap,		\
	   transp, netR, evap_f, transp_f, D, DelCanStor, theta, theta2, theta3, theta_available, ra, \
	   psi_ae, bclambda, rootdepth, froot1, froot2, froot3, d1, d2, d3, thetar, fc3, \
	   pTrp1, pTrp2, pTrp3, veg_p, d2HevapT_f, d18OevapT_f, AgeevapT_f, \
	   d2HevapI_f, d18OevapI_f, AgeevapI_f, d2Hcanopy_f, d18Ocanopy_f, Agecanopy_f, \
	   poros1, poros2, poros3,					\
	   kTB_L1, kMW_L1, kTB_L2, kMW_L2, \
	   theta_MW1, theta_MW2, dth1, dth2, dth3)		\
  shared(nsp, atm, ctrl, trck, dt, thre, std::cout)

  {
    thre = omp_get_num_threads();
#pragma omp single
    printf("\nnum threads %d: \n", thre);
#pragma omp for nowait
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {


      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;

      /*--------*/
      nsp = fForest->getNumSpecies();

      veg_p = 0;
      treeheight = 0;
      evap_f = 0;
      transp_f = 0;

      theta = _soilmoist1->matrix[r][c]; //soil moisture at time t
      theta2 = _soilmoist2->matrix[r][c];
      theta3 = _soilmoist3->matrix[r][c];
      thetar = _theta_r->matrix[r][c];
      fc3 = _fieldcapL3->matrix[r][c];


      psi_ae = _psi_ae->matrix[r][c];
      bclambda = _BClambda->matrix[r][c];

      d1 = _depth_layer1->matrix[r][c];
      d2 = _depth_layer2->matrix[r][c];
      d3 = _soildepth->matrix[r][c] - d1 - d2;

      poros1 = _porosityL1->matrix[r][c];
      poros2 = _porosityL2->matrix[r][c];
      poros3 = _porosityL3->matrix[r][c];
      
      // Tracking: initialize summed values
      d2HevapT_f = 0;
      d18OevapT_f = 0;
      AgeevapT_f = 0;
      d2HevapI_f = 0;
      d18OevapI_f = 0;
      AgeevapI_f = 0;
      d2Hcanopy_f = 0;
      d18Ocanopy_f = 0;
      Agecanopy_f = 0;
      kTB_L1 = 0;
      kTB_L2 = 0;
      kMW_L1 = 0;
      kMW_L2 = 0;
      theta_MW1 = 0;
      theta_MW2 = 0;

      if(ctrl.sw_trck and ctrl.sw_TPD){
	theta_MW1 = _moist_MW1->matrix[r][c];
	theta_MW2 = _moist_MW2->matrix[r][c];
      }    
      // --- LOOP ON SPECIES --------------------------------------
      for (s = 0; s < nsp; s++) {
     
	p = fForest->getPropSpecies(s, r, c);
	if (p == 0)
	  continue; //if no species j present, continue
     
	DelCanStor = 0;
	D = 0;
	evap = 0;
	transp = 0;
	netR = 0;
     
	if (s == nsp - 1) { //if this is bare ground set D to precip and skip the tree stuff
       
	  D = atm.getPrecipitation()->matrix[r][c];

	} else {
       
	  wind = atm.getWindSpeed()->matrix[r][c];
       
	  veg_p += p;
	  treeheight = max<REAL8>(0.01, fForest->getTreeHeight(s, r, c));
       
	  /*TODO: Tentative relationship between forest height and wind velocity profile parameters*/
	  za = treeheight + 2;
	  z0o = powl(treeheight, 1.19) * 0.057544; //powl( 10, -1.24+1.19*log10l(treeheight) );     //treeheight > 1 ? 0.1 : treeheight * 0.1;
	  zdo = powl(treeheight, 0.98) * 0.707946; //powl( 10, 0.98*log10l(treeheight)-0.15); //treeheight > 1 ? 0.1 : treeheight * 0.7;
       
	  rootdepth = _soildepth->matrix[r][c];
	  theta = _soilmoist1->matrix[r][c]; //soil moisture at time t
	  theta2 = _soilmoist2->matrix[r][c];
	  theta3 = _soilmoist3->matrix[r][c];
	  froot1 = fForest->getRootFrac1(s)->matrix[r][c];
	  froot2 = fForest->getRootFrac2(s)->matrix[r][c];
	  //froot1 = _rootfrac1->matrix[r][c];
	  //froot2 = _rootfrac2->matrix[r][c];
	  froot3 = 1 - froot1 - froot2;

	  theta_available = (theta-thetar) * froot1 + 
	    (theta2-thetar) * froot2 + (theta3-thetar) * froot3;

	  //cout << "theta_a : " << theta_available << ", theta1 : " << theta <<
	  //  ", theta2 : " << theta2 << ", theta3 : " << theta3 << endl;

	  //root depth is the depth of layers that contain 95% of roots
	  if (froot1 > 0.95)
	    rootdepth = froot1 > d1;
	  else if ((froot1 + froot2) > 0.95)
	    rootdepth = d1 + d2;
       
	  ra = CalcAerodynResist(wind, za, 0, 0, z0o, zdo, treeheight,
				 fForest->getLAISpecies(s, r, c),
				 getCanopyTemp(s)->matrix[r][c],
				 atm.getTemperature()->matrix[r][c], ctrl.toggle_ra,
				 false);
       
	  fForest->CanopyInterception(atm, ctrl, DelCanStor, D, s, r, c); //calculates canopy interception and trascolation

	  // Tracking: canopy storage
	  if (ctrl.sw_trck and DelCanStor > RNDOFFERR){
	    if(ctrl.sw_2H)
	      fForest->setd2Hcanopy(s, r, c, 
				    fForest->getIntercWater(s, r, c) <= RNDOFFERR ? -1000 :
				    trck.InputMix(fForest->getIntercWater(s, r, c) - DelCanStor,
						  fForest->getd2Hcanopy(s)->matrix[r][c],
						  DelCanStor, 
						  atm.getd2Hprecip()->matrix[r][c]));
	 
	    if(ctrl.sw_18O)
	      fForest->setd18Ocanopy(s, r, c, 
				    fForest->getIntercWater(s, r, c) <= RNDOFFERR ? -1000 :
				     trck.InputMix(fForest->getIntercWater(s, r, c) - DelCanStor,
						   fForest->getd18Ocanopy(s)->matrix[r][c],
						   DelCanStor, 
						   atm.getd18Oprecip()->matrix[r][c]));
	 
	    if(ctrl.sw_Age)
	      fForest->setAgecanopy(s, r, c, 
				    fForest->getIntercWater(s, r, c) <= RNDOFFERR ? 0 :
				    trck.InputMix(fForest->getIntercWater(s, r, c) - DelCanStor,
						  fForest->getAgecanopy(s)->matrix[r][c],
						  DelCanStor,0));
	 
	  }
	  // -----------------------------------
 
	  // TODO : Implement fractionation in canopy interception in SolveCanopyEnergyBalance
	  fForest->SolveCanopyEnergyBalance(*this, atm, ctrl, 
					    thetar, rootdepth, psi_ae, bclambda, 
					    ra, DelCanStor, evap, transp, netR, 
					    s, r, c);
       
	  // Canopy evap-related update etc.
	  _CanopyStorage->matrix[r][c] += DelCanStor * p;
       
	  _Rn_sum->matrix[r][c] += netR * p ;
       
	  if (_CanopyStorage->matrix[r][c] < RNDOFFERR)
	    _CanopyStorage->matrix[r][c] = 0.0;
       
	  evap_f += evap * p; //evaporation at t=t+1
	  // ---------------------------
	  if(ctrl.sw_trck){
	    if(ctrl.sw_2H)
	      d2Hcanopy_f += fForest->getd2Hcanopy(s)->matrix[r][c]*p*fForest->getIntercWater(s, r, c);
	    if(ctrl.sw_18O)
	      d18Ocanopy_f += fForest->getd18Ocanopy(s)->matrix[r][c]*p*fForest->getIntercWater(s, r, c);
	    if(ctrl.sw_Age)
	      Agecanopy_f += fForest->getAgecanopy(s)->matrix[r][c]*p*fForest->getIntercWater(s, r, c);
	  }

      
	  // Tracking (evapI) TODO : fractionation in canopy (does not affect soil)
	  if(ctrl.sw_trck){
	    if(ctrl.sw_2H){
	      fForest->setd2HevapI(s, r, c, 
				   evap <= RNDOFFERR ? -1000 : fForest->getd2Hcanopy(s)->matrix[r][c]);
	      d2HevapI_f += fForest->getd2HevapI(s)->matrix[r][c] * p * evap ;
	    }
	    if(ctrl.sw_18O){
	      fForest->setd18OevapI(s, r, c, 
				    evap <= RNDOFFERR ? -1000 : fForest->getd18Ocanopy(s)->matrix[r][c]);
	      d18OevapI_f += fForest->getd18OevapI(s)->matrix[r][c] * p * evap ;
	    }
	    if(ctrl.sw_Age){
	      fForest->setAgeevapI(s, r, c, 
				   evap <= RNDOFFERR ? 0 : fForest->getAgecanopy(s)->matrix[r][c]);
	      AgeevapI_f += fForest->getAgeevapI(s)->matrix[r][c] * p * evap ;
	    }
	  }
       
	  // Transpiration-related update etc.
	  transp_f += transp * p;
       
	  pTrp1 = ((theta-thetar)*froot1) / theta_available;
	  pTrp2 = ((theta2-thetar)*froot2) / theta_available;
	  pTrp3 = ((theta3-thetar)*froot3) / theta_available;
       
	  dth1 = transp * p * dt * pTrp1 /d1;
	  dth2 = transp * p * dt * pTrp2 /d2;
	  dth3 = transp * p * dt * pTrp3 /d3;

	  theta  -= dth1 ; //soil moisture at t=t+1
	  theta2 -= dth2 ; //soil moisture at t=t+1
	  theta3 -= dth3 ; //soil moisture at t=t+1

	  // Tracking (evapT):
       
	  if(ctrl.sw_trck){
	    // If two-pore domain activated: weighted contribution using
	    // relative pore volume + transfer-mixing to tightly-bound domain if needed		
	    if(ctrl.sw_TPD){
	      //k_MW = (poros - theta_MW) / (poros - thetar);
	      //k_MW3 = (fc - theta_MW) / (fc - thetar);
	      // Check if the water content in mobile water is limiting,
	      // i.e. as compared to the potential root uptakes: k_MW*(theta_old-theta_new)
	      // Note that here the tightly-bound domain cannot be limiting, since tightly-bound has to be
	      // full to have fast (so at a minimum transpi takes fully from tightly-bound)
	      if(dth1 > RNDOFFERR){
		kMW_L1 = max<double>(0,(_soilmoist1->matrix[r][c] - theta_MW1) /
				     (_soilmoist1->matrix[r][c] - thetar));
		kTB_L1 = max<double>(0,1 - kMW_L1);
	      }
	      if(dth2 > RNDOFFERR) {
		kMW_L2 = min<double>(0, (_soilmoist2->matrix[r][c] - theta_MW2) /
				     (_soilmoist2->matrix[r][c] - thetar));
		kTB_L2 = max<double>(0,1 - kMW_L2);
	      }
	      if(ctrl.sw_2H){
		fForest->setd2HevapT(s, r, c,
				     transp <= RNDOFFERR ? -1000 :
				     pTrp1*(kTB_L1*trck.getd2H_TB1()->matrix[r][c]+
					    kMW_L1*trck.getd2H_MW1()->matrix[r][c])+
				     pTrp2*(kTB_L2*trck.getd2H_TB2()->matrix[r][c]+
					    kMW_L2*trck.getd2H_MW2()->matrix[r][c])+
				     pTrp3*trck.getd2Hsoil3()->matrix[r][c]);
		d2HevapT_f += fForest->getd2HevapT(s)->matrix[r][c] * p * transp ;		       
	      }
	      if(ctrl.sw_18O){
		fForest->setd18OevapT(s, r, c,
				      transp <= RNDOFFERR ? -1000 :
				      pTrp1*(kTB_L1*trck.getd18O_TB1()->matrix[r][c]+
					     kMW_L1*trck.getd18O_MW1()->matrix[r][c])+
				      pTrp2*(kTB_L2*trck.getd18O_TB2()->matrix[r][c]+
					     kMW_L2*trck.getd18O_MW2()->matrix[r][c])+
				      pTrp3*trck.getd18Osoil3()->matrix[r][c]);
		d18OevapT_f += fForest->getd18OevapT(s)->matrix[r][c] * p * transp ;
	      }
	      if(ctrl.sw_Age){
		fForest->setAgeevapT(s, r, c,
				     transp <= RNDOFFERR ? 0 :
				     pTrp1*(kTB_L1*trck.getAge_TB1()->matrix[r][c]+
					    kMW_L1*trck.getAge_MW1()->matrix[r][c])+
				     pTrp2*(kTB_L2*trck.getAge_TB2()->matrix[r][c]+
					    kMW_L2*trck.getAge_MW2()->matrix[r][c])+
				     pTrp3*trck.getAgesoil3()->matrix[r][c]);
		AgeevapT_f += fForest->getAgeevapT(s)->matrix[r][c] * p * transp ;
	      }
	      // If uptake from tightly-bound domain, refill from mobile water (potentially emptying
	      // the latter), no change for layer hydro but mixing!
	      trck.MixingTPD_postET(*this, ctrl, dth1, dth2,
				    kTB_L1, kTB_L2, kMW_L1, kMW_L2, r, c);
	   
	   
	    } else {
	      // Else, it assumes total mixing of the water pulled from each soil layer	   
	      if(ctrl.sw_2H){
		fForest->setd2HevapT(s, r, c,
				     transp <= RNDOFFERR ? -1000 :
				     pTrp1*trck.getd2Hsoil1()->matrix[r][c]+
				     pTrp2*trck.getd2Hsoil2()->matrix[r][c]+
				     pTrp3*trck.getd2Hsoil3()->matrix[r][c]);
		d2HevapT_f += fForest->getd2HevapT(s)->matrix[r][c] * p * transp ;
	      }
	      if(ctrl.sw_18O){
		fForest->setd18OevapT(s, r, c,
				      transp <= RNDOFFERR ? -1000 :
				      pTrp1*trck.getd18Osoil1()->matrix[r][c]+
				      pTrp2*trck.getd18Osoil2()->matrix[r][c]+
				      pTrp3*trck.getd18Osoil3()->matrix[r][c]);
		d18OevapT_f += fForest->getd18OevapT(s)->matrix[r][c] * p * transp ;
	      }
	      if(ctrl.sw_Age){
		fForest->setAgeevapT(s, r, c,
				     transp <= RNDOFFERR ? 0 :
				     pTrp1*trck.getAgesoil1()->matrix[r][c]+
				     pTrp2*trck.getAgesoil2()->matrix[r][c]+
				     pTrp3*trck.getAgesoil3()->matrix[r][c]);
		AgeevapT_f += fForest->getAgeevapT(s)->matrix[r][c] * p * transp ;
	      }
	    }
	  }
	  // Update soil moisture objects
	  _soilmoist1->matrix[r][c] = theta;
	  _soilmoist2->matrix[r][c] = theta2;
	  _soilmoist3->matrix[r][c] = theta3;
	      
	} // end if bare / veg
     
	Tp = atm.getTemperature()->matrix[r][c];
	maxTp = atm.getMaxTemperature()->matrix[r][c];
	minTp = atm.getMinTemperature()->matrix[r][c];
	sno_rain_thres = atm.getRainSnowTempThreshold();
     
	if(maxTp <= sno_rain_thres){
	  snow = D;
	  rain = 0;
	}
	else if(minTp > sno_rain_thres){
	  rain = D;
	  snow = 0;
	}
	else{
	  snow = D * max<REAL8>(0.0, (sno_rain_thres - minTp) /(maxTp - minTp));
	  rain = D - snow;
	}
     
	// Water tracking
	if(ctrl.sw_trck){
	  // Used later in SolveSurfaceFluxes
	  _FluxCnptoSrf->matrix[r][c] += rain * p * dt;
	  _FluxCnptoSnow->matrix[r][c] += snow * p * dt;
       
	  // Mix in surface pool from rain throughfall
	  trck.MixingV_through(atm, *this, ctrl, rain, p, r, c);
	}
     
	_snow->matrix[r][c] +=  snow * dt * p;
	_ponding->matrix[r][c] += rain * dt * p;
     
      } //end for over species
   
      _Evaporation->matrix[r][c] = evap_f + transp_f; //total evaporation for the entire cell
      // Vegetation-summed values
      _Transpiration_all->matrix[r][c]  = transp_f ;
      _EvaporationI_all->matrix[r][c] = evap_f ;
      // applies for tracking as well
      // summed evapT 2H, 18O and Age weighted by flux magnitude and cover ONLY over vegetated fraction
      // (otherwise results makes no sense whre bare frac is significant)
      if(ctrl.sw_trck && ctrl.sw_2H){
	trck.setd2HevapT_sum(r, c, transp_f > RNDOFFERR ? d2HevapT_f / transp_f : -1000); 
	trck.setd2HevapI_sum(r, c, evap_f > RNDOFFERR ? d2HevapI_f / evap_f : -1000); 
	trck.setd2Hcanopy_sum(r, c,_CanopyStorage->matrix[r][c] > RNDOFFERR ? 
			      d2Hcanopy_f / _CanopyStorage->matrix[r][c] : -1000); 
      }
      if(ctrl.sw_trck && ctrl.sw_18O){
	trck.setd18OevapT_sum(r, c, transp_f > RNDOFFERR ? d18OevapT_f / transp_f : -1000);
	trck.setd18OevapI_sum(r, c, evap_f > RNDOFFERR ? d18OevapI_f / evap_f : -1000); 
	trck.setd18Ocanopy_sum(r, c, _CanopyStorage->matrix[r][c] > RNDOFFERR ? 
			       d18Ocanopy_f / _CanopyStorage->matrix[r][c] : -1000); 
      }
      if(ctrl.sw_trck && ctrl.sw_Age){
	trck.setAgeevapT_sum(r, c, transp_f > RNDOFFERR ? AgeevapT_f / transp_f : 0);
	trck.setAgeevapI_sum(r, c, evap_f > RNDOFFERR ? AgeevapI_f / evap_f : 0); 
	trck.setAgecanopy_sum(r, c, _CanopyStorage->matrix[r][c] > RNDOFFERR ? 
			      Agecanopy_f / _CanopyStorage->matrix[r][c] : 0); 
      }
   
    }//end for
  }//end omp parallel
   
   
  return EXIT_SUCCESS;
}
