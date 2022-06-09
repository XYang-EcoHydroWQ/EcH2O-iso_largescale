#!/bin/bash

COUNT=1
NITER=11

while [ $COUNT -lt $NITER ]
do
	echo "Running iteration ${COUNT}"
	xterm -e ech2o config.ini
	
	echo "finished run, copying files"
	cp -f ./Results/root0_00.365 ./Spatial/root_0.map
	cp -f ./Results/p0_00000.365 ./Spatial/p_0.map
	cp -f ./Results/ntr0_000.365 ./Spatial/ntr_0.map
	cp -f ./Results/lai0_000.365 ./Spatial/lai_0.map
	cp -f ./Results/hgt0_000.365 ./Spatial/hgt_0.map
	cp -f ./Results/bas0_000.365 ./Spatial/bas_0.map
	cp -f ./Results/age0_000.365 ./Spatial/age_0.map

	cp -f ./Results/SWE00000.365 ./Spatial/SWE.map
	cp -f ./Results/SWC1_000.365 ./Spatial/Soil_moisture_1.map
	cp -f ./Results/SWC2_000.365 ./Spatial/Soil_moisture_2.map
	cp -f ./Results/SWC3_000.365 ./Spatial/Soil_moisture_3.map
	cp -f ./Results/Ts000000.365 ./Spatial/soiltemp.map
	cp -f ./Results/Q0000000.365 ./Spatial/streamflow.map

	cat ./Results/lai_0.tab >> ./Results/laiaccum.txt
	cat ./Results/NPP_0.tab >> ./Results/NPPaccum.txt
	cat ./Results/SoilMoistureAv.tab >> ./Results/SWCaccum.txt


	let "COUNT++"
done
	echo "Finished simulation\n"


