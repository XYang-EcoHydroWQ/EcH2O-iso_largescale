################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Budgets/AccountFluxes.cpp \
../Budgets/AccountRelArea.cpp \
../Budgets/AccountStorages.cpp \
../Budgets/MassBalanceError.cpp \
../Budgets/TotalEvaporation.cpp \
../Budgets/TotalEvaporationS.cpp \
../Budgets/TotalEvaporationI.cpp \
../Budgets/TotalGWtoChn.cpp \
../Budgets/TotalTranspiration.cpp \
../Budgets/TotalLeakage.cpp \
../Budgets/TotalOvlndFlow.cpp \
../Budgets/TotalPrecipitation.cpp \
../Budgets/TotalSaturationArea.cpp \
../Budgets/TotalSrftoChn.cpp \
../Budgets/TotalStorage.cpp \
../Budgets/totalGrndFlow.cpp \
../Budgets/TrckBalanceError.cpp

OBJS += \
./Budgets/AccountFluxes.o \
./Budgets/AccountRelArea.o \
./Budgets/AccountStorages.o \
./Budgets/MassBalanceError.o \
./Budgets/TotalEvaporation.o \
./Budgets/TotalEvaporationS.o \
./Budgets/TotalEvaporationI.o \
./Budgets/TotalGWtoChn.o \
./Budgets/TotalTranspiration.o \
./Budgets/TotalLeakage.o \
./Budgets/TotalOvlndFlow.o \
./Budgets/TotalPrecipitation.o \
./Budgets/TotalSaturationArea.o \
./Budgets/TotalSrftoChn.o \
./Budgets/TotalStorage.o \
./Budgets/totalGrndFlow.o \
./Budgets/TrckBalanceError.o

CPP_DEPS += \
./Budgets/AccountFluxes.d \
./Budgets/AccountRelArea.d \
./Budgets/AccountStorages.d \
./Budgets/MassBalanceError.d \
./Budgets/TotalEvaporation.d \
./Budgets/TotalEvaporationS.d \
./Budgets/TotalEvaporationI.d \
./Budgets/TotalGWtoChn.d \
./Budgets/TotalTranspiration.d \
./Budgets/TotalLeakage.d \
./Budgets/TotalOvlndFlow.d \
./Budgets/TotalPrecipitation.d \
./Budgets/TotalSaturationArea.d \
./Budgets/TotalSrftoChn.d \
./Budgets/TotalStorage.d \
./Budgets/totalGrndFlow.d \
./Budgets/TrckBalanceError.d


# Each subdirectory must supply rules for building sources it contributes
Budgets/%.o: ../Budgets/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	x86_64-w64-mingw32-g++ -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


