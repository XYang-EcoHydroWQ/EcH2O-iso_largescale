################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Grid/SortGrid.cpp \
../Grid/SortGridForest.cpp \
../Grid/SortGridLDD.cpp \
../Grid/grid.cpp 

OBJS += \
./Grid/SortGrid.o \
./Grid/SortGridForest.o \
./Grid/SortGridLDD.o \
./Grid/grid.o 

CPP_DEPS += \
./Grid/SortGrid.d \
./Grid/SortGridForest.d \
./Grid/SortGridLDD.d \
./Grid/grid.d 


# Each subdirectory must supply rules for building sources it contributes
Grid/%.o: ../Grid/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	x86_64-w64-mingw32-g++ -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


