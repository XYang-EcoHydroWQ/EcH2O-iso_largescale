################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Constructors/AtmosphConstruct.cpp \
../Constructors/BasinConstruct.cpp \
../Constructors/CheckMaps.cpp \
../Constructors/ForestConstruct.cpp \
../Constructors/GroveConstruct.cpp \
../Constructors/ReportConstruct.cpp \
../Constructors/TrackingConstruct.cpp \
../Constructors/checkForestDatabase.cpp 

OBJS += \
./Constructors/AtmosphConstruct.o \
./Constructors/BasinConstruct.o \
./Constructors/CheckMaps.o \
./Constructors/ForestConstruct.o \
./Constructors/GroveConstruct.o \
./Constructors/ReportConstruct.o \
./Constructors/TrackingConstruct.o \
./Constructors/checkForestDatabase.o 

CPP_DEPS += \
./Constructors/AtmosphConstruct.d \
./Constructors/BasinConstruct.d \
./Constructors/CheckMaps.d \
./Constructors/ForestConstruct.d \
./Constructors/GroveConstruct.d \
./Constructors/ReportConstruct.d \
./Constructors/TrackingConstruct.d \
./Constructors/checkForestDatabase.d 


# Each subdirectory must supply rules for building sources it contributes
Constructors/%.o: ../Constructors/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	x86_64-w64-mingw32-g++ -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


