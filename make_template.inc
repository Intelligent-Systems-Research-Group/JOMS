CXX = /clang/install/bin/clang++
LINK =/clang/install bin/clang++
NVCC = /usr/local/cuda/bin/nvcc
SRC ?= src
#CUDAHOME += /usr/local/cuda
#CUDAHOME += -L/usr/local/cuda-9.2/targets/x86_64-linux/lib/ -lcudart
#/release
TERRAHOME ?= /terra/install
TERRA ?= /terra/install/bin/terra
#TERRAHOME ?= /home/dl/projects/Optlang/terra_working
#TERRA ?= /home/dl/projects/Optlang/terra_working

ifeq ($(UNAME), Linux)
  S4 = 64
endif
# -g -O3 -DNDEBUG -ltcmalloc -march=native
FLAGS += -O3 -march=native -std=c++11 -I$(SRC) -I$(SRC)/cutil/inc -I../Optlang/API/release/include -I$(TERRAHOME)/include -I/usr/local/cuda-9.2/targets/x86_64-linux/include
FLAGS += -I/fbxsdk/include -I/opensubdiv/install/include 
FLAGS += -I/nanoflann/install/include -I/embree/install/include
FLAGS += -I/cxxopts/install/include -I/json/install/include  
NVCCFLAGS += -arch sm_61 $(FLAGS)
LFLAGS += -std=c++11 -L/usr/local/cuda-9.2/targets/x86_64-linux/lib/ -L/fbxsdk/lib/gcc4/x64/release -lfbxsdk -lcudart 
LFLAGS += -L/terra/install/lib -L/opensubdiv/install/lib 
LFLAGS += -L/embree/install/lib 
LFLAGS += -lterra -ldl -ltinfo -lz -lembree3 
LFLAGS += -L/clang/install/lib -lLLVMLTO -lLLVMPasses -lLLVMObjCARCOpts -lLLVMSymbolize -lLLVMDebugInfoPDB -lLLVMDebugInfoDWARF -lLLVMMIRParser -lLLVMFuzzMutate -lLLVMCoverage -lLLVMTableGen -lLLVMDlltoolDriver -lLLVMOrcJIT -lLLVMXCoreDisassembler -lLLVMXCoreCodeGen -lLLVMXCoreDesc -lLLVMXCoreInfo -lLLVMXCoreAsmPrinter -lLLVMSystemZDisassembler -lLLVMSystemZCodeGen -lLLVMSystemZAsmParser -lLLVMSystemZDesc -lLLVMSystemZInfo -lLLVMSystemZAsmPrinter -lLLVMSparcDisassembler -lLLVMSparcCodeGen -lLLVMSparcAsmParser -lLLVMSparcDesc -lLLVMSparcInfo -lLLVMSparcAsmPrinter -lLLVMPowerPCDisassembler -lLLVMPowerPCCodeGen -lLLVMPowerPCAsmParser -lLLVMPowerPCDesc -lLLVMPowerPCInfo -lLLVMPowerPCAsmPrinter -lLLVMNVPTXCodeGen -lLLVMNVPTXDesc -lLLVMNVPTXInfo -lLLVMNVPTXAsmPrinter -lLLVMMSP430CodeGen -lLLVMMSP430Desc -lLLVMMSP430Info -lLLVMMSP430AsmPrinter -lLLVMMipsDisassembler -lLLVMMipsCodeGen -lLLVMMipsAsmParser -lLLVMMipsDesc
LFLAGS += -lLLVMMipsInfo -lLLVMMipsAsmPrinter -lLLVMLanaiDisassembler -lLLVMLanaiCodeGen -lLLVMLanaiAsmParser -lLLVMLanaiDesc -lLLVMLanaiAsmPrinter -lLLVMLanaiInfo -lLLVMHexagonDisassembler -lLLVMHexagonCodeGen -lLLVMHexagonAsmParser -lLLVMHexagonDesc -lLLVMHexagonInfo -lLLVMBPFDisassembler -lLLVMBPFCodeGen -lLLVMBPFAsmParser -lLLVMBPFDesc -lLLVMBPFInfo -lLLVMBPFAsmPrinter -lLLVMARMDisassembler -lLLVMARMCodeGen -lLLVMARMAsmParser -lLLVMARMDesc -lLLVMARMInfo -lLLVMARMAsmPrinter -lLLVMARMUtils -lLLVMAMDGPUDisassembler -lLLVMAMDGPUCodeGen -lLLVMAMDGPUAsmParser -lLLVMAMDGPUDesc -lLLVMAMDGPUInfo -lLLVMAMDGPUAsmPrinter -lLLVMAMDGPUUtils -lLLVMAArch64Disassembler -lLLVMAArch64CodeGen -lLLVMAArch64AsmParser -lLLVMAArch64Desc -lLLVMAArch64Info -lLLVMAArch64AsmPrinter -lLLVMAArch64Utils -lLLVMObjectYAML -lLLVMLibDriver -lLLVMOption -lLLVMWindowsManifest -lLLVMX86Disassembler -lLLVMX86AsmParser -lLLVMX86CodeGen -lLLVMGlobalISel -lLLVMSelectionDAG -lLLVMAsmPrinter -lLLVMX86Desc -lLLVMMCDisassembler -lLLVMX86Info -lLLVMX86AsmPrinter -lLLVMX86Utils -lLLVMMCJIT -lLLVMLineEditor -lLLVMInterpreter -lLLVMExecutionEngine -lLLVMRuntimeDyld -lLLVMCodeGen -lLLVMTarget -lLLVMCoroutines -lLLVMipo -lLLVMInstrumentation -lLLVMVectorize -lLLVMScalarOpts -lLLVMLinker -lLLVMIRReader -lLLVMAsmParser -lLLVMInstCombine -lLLVMBitWriter -lLLVMAggressiveInstCombine -lLLVMTransformUtils -lLLVMAnalysis -lLLVMProfileData -lLLVMObject -lLLVMMCParser -lLLVMMC -lLLVMDebugInfoCodeView -lLLVMDebugInfoMSF -lLLVMBitReader -lLLVMCore -lLLVMBinaryFormat -lLLVMSupport -lLLVMDemangle
LFLAGS += -lrt -ldl -lpthread -lm 
UNAME := $(shell uname)

#-L/home/caduser/samuel/fbxsdk/lib/gcc4/x64/release
#-I../fbxsdk/include
ifeq ($(UNAME), Darwin)
  LFLAGS += -pagezero_size 10000 -image_base 100000000
endif

ifeq ($(UNAME), Linux)
#  LFLAGS += -Wl,-rpath,$(CUDAHOME)/lib$(S4) -pthread
  LFLAGS += -Wl,-rpath,/usr/local/cuda-9.2/targets/x86_64-linux/lib/ -pthread
  S4 = 64
endif

.SUFFIXES:
.PHONY:	all purge clean

all: $(EXECUTABLE)

OPTLIB = /Optlang/API/release/lib/libOptDev.a

build/%.o:	$(SRC)/%.cpp
	@mkdir -p build
	$(CXX) $(FLAGS) -MP -MD $< -c -o $@
	#$(CXX) --cuda-gpu-arch=sm_61 $(FLAGS) -MP -MD $< -c -o $@
build/%.o:	$(SRC)/%.cu
	@mkdir -p build
	$(NVCC) $(NVCCFLAGS) $< -c -o $@

$(EXECUTABLE):	$(OBJS) $(OPTLIB)
	$(CXX) $^ -o $@ $(LFLAGS)

clean:
	rm -rf build/*.o build/*.d
	rm -rf $(EXECUTABLE)

purge: clean
	rm -rf build/*

# dependency rules
include $(wildcard build/*.d)
