################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Hydro/CalcCatchArea.cpp \
../Hydro/CalcFracMobileWater.cpp \
../Hydro/CalcInitialStreamStorage.cpp \
../Hydro/CalcSoilResist.cpp \
../Hydro/CalcTPDMoisture.cpp \
../Hydro/CalcRootDistrib.cpp \
../Hydro/CalcSoilMoistureProfile.cpp \
../Hydro/CalcPropLayers.cpp \
../Hydro/CalculateForestGrowth.cpp \
../Hydro/CalculateSatArea.cpp \
../Hydro/CanopyInterception.cpp \
../Hydro/GWrouting.cpp \
../Hydro/Infilt_GreenAmpt.cpp \
../Hydro/ExtraGWDynamics.cpp \
../Hydro/KinematicWave.cpp \
../Hydro/SnowOutputPhase.cpp \
../Hydro/SoilEvapotranspiration.cpp \
../Hydro/SoilWaterRedistribution.cpp \
../Hydro/SolveCanopyFluxes.cpp \
../Hydro/SolveSurfaceEnergyBalance.cpp \
../Hydro/UpdateSnowPack.cpp \
../Hydro/SolveChannelEnergyBalance.cpp \
../Hydro/InOutMix_Temperature.cpp \
../Hydro/ChannelEvaporation.cpp

OBJS += \
./Hydro/CalcCatchArea.o \
./Hydro/CalcFracMobileWater.o \
./Hydro/CalcInitialStreamStorage.o \
./Hydro/CalcTPDMoisture.o \
./Hydro/CalcSoilResist.o \
./Hydro/CalcRootDistrib.o \
./Hydro/CalculateForestGrowth.o \
./Hydro/CalcSoilMoistureProfile.o \
./Hydro/CalcPropLayers.o \
./Hydro/CalculateSatArea.o \
./Hydro/CanopyInterception.o \
./Hydro/GWrouting.o \
./Hydro/Infilt_GreenAmpt.o \
./Hydro/ExtraGWDynamics.o \
./Hydro/KinematicWave.o \
./Hydro/SnowOutputPhase.o \
./Hydro/SoilEvapotranspiration.o \
./Hydro/SoilWaterRedistribution.o \
./Hydro/SolveCanopyFluxes.o \
./Hydro/SolveSurfaceEnergyBalance.o \
./Hydro/UpdateSnowPack.o \
./Hydro/SolveChannelEnergyBalance.o \
./Hydro/InOutMix_Temperature.o \
./Hydro/ChannelEvaporation.o

CPP_DEPS += \
./Hydro/CalcCatchArea.d \
./Hydro/CalcFracMobileWater.d \
./Hydro/CalcInitialStreamStorage.d \
./Hydro/CalcTPDMoisture.d \
./Hydro/CalcSoilResist.d \
./Hydro/CalcRootDistrib.d \
./Hydro/CalcSoilMoistureProfile.d \
./Hydro/CalcPropLayers.d \
./Hydro/CalculateForestGrowth.d \
./Hydro/CalculateSatArea.d \
./Hydro/CanopyInterception.d \
./Hydro/GWrouting.d \
./Hydro/Infilt_GreenAmpt.d \
./Hydro/ExtraGWDynamics.d \
./Hydro/KinematicWave.d \
./Hydro/SnowOutputPhase.d \
./Hydro/SoilEvapotranspiration.d \
./Hydro/SoilWaterRedistribution.d \
./Hydro/SolveCanopyFluxes.d \
./Hydro/SolveSurfaceEnergyBalance.d \
./Hydro/UpdateSnowPack.d \
./Hydro/SolveChannelEnergyBalance.d \
./Hydro/InOutMix_Temperature.d \
./Hydro/ChannelEvaporation.d


# Each subdirectory must supply rules for building sources it contributes
Hydro/%.o: ../Hydro/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


