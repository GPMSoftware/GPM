################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BulkEditorDialog.cpp \
../src/CpuGpuPanel.cpp \
../src/DatabasePanel.cpp \
../src/InputPanel.cpp \
../src/MatchPanel.cpp \
../src/MegaChoiceDialog.cpp \
../src/MetadataDialog.cpp \
../src/MyFrame.cpp \
../src/MyGrid.cpp \
../src/ResultsPanel.cpp \
../src/SelectProfilesPanel.cpp \
../src/SpmPanel.cpp \
../src/gmatch.cpp 

OBJS += \
./src/BulkEditorDialog.o \
./src/CpuGpuPanel.o \
./src/DatabasePanel.o \
./src/InputPanel.o \
./src/MatchPanel.o \
./src/MegaChoiceDialog.o \
./src/MetadataDialog.o \
./src/MyFrame.o \
./src/MyGrid.o \
./src/ResultsPanel.o \
./src/SelectProfilesPanel.o \
./src/SpmPanel.o \
./src/gmatch.o 

CPP_DEPS += \
./src/BulkEditorDialog.d \
./src/CpuGpuPanel.d \
./src/DatabasePanel.d \
./src/InputPanel.d \
./src/MatchPanel.d \
./src/MegaChoiceDialog.d \
./src/MetadataDialog.d \
./src/MyFrame.d \
./src/MyGrid.d \
./src/ResultsPanel.d \
./src/SelectProfilesPanel.d \
./src/SpmPanel.d \
./src/gmatch.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/cuda/include -I"/home/gareth/GPMSoftware/eclipse_workspace" -I"/usr/local/lib/wx/include/gtk2-unicode-2.9" -I/usr/local/include/wx-2.9 -I/usr/include/mysql -I/usr/local/geotrans/include -O0 -g3 -Wall -c -fmessage-length=0 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -D__WXGTK__ -pthread -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


