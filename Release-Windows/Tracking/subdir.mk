################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tracking/CalcTrck_L1L2.cpp \
../Tracking/CalcTrck_SoilAv.cpp \
../Tracking/CheckMapsTrck.cpp \
../Tracking/Fractionation_Esoil.cpp \
../Tracking/IncrementAge.cpp \
../Tracking/MixingTPD_postET.cpp \
../Tracking/MixingH.cpp \
../Tracking/MixingV_down.cpp \
../Tracking/MixingV_evapS.cpp \
../Tracking/MixingV_seep.cpp \
../Tracking/MixingV_snow.cpp \
../Tracking/MixingV_through.cpp \
../Tracking/MixingV_up.cpp \
../Tracking/OutletVals.cpp 

OBJS += \
./Tracking/CalcTrck_L1L2.o \
./Tracking/CalcTrck_SoilAv.o \
./Tracking/CheckMapsTrck.o\
./Tracking/Fractionation_Esoil.o \
./Tracking/IncrementAge.o \
./Tracking/MixingTPD_postET.o \
./Tracking/MixingH.o \
./Tracking/MixingV_down.o \
./Tracking/MixingV_evapS.o \
./Tracking/MixingV_seep.o \
./Tracking/MixingV_snow.o \
./Tracking/MixingV_through.o \
./Tracking/MixingV_up.o \
./Tracking/OutletVals.o

CPP_DEPS += \
./Tracking/CalcTrck_L1L2.d \
./Tracking/CalcTrck_SoilAv.d \
./Tracking/CheckMapsTrck.d \
./Tracking/Fractionation_Esoil.d \
./Tracking/IncrementAge.d \
./Tracking/MixingTPD_postET.d \
./Tracking/MixingH.d \
./Tracking/MixingV_down.d \
./Tracking/MixingV_evapS.d \
./Tracking/MixingV_seep.d \
./Tracking/MixingV_snow.d \
./Tracking/MixingV_through.d \
./Tracking/MixingV_up.d  \
./Tracking/OutletVals.d

# Each subdirectory must supply rules for building sources it contributes
Tracking/%.o: ../Tracking/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	x86_64-w64-mingw32-gcc -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


