EXECUTABLE = trainer
OBJS = build/main.o build/corpora.o build/MeshTopology.o build/rotations.o build/skeleton.o build/VertexGroupLoader.o  build/SubdivEvaluator.o build/particles.o build/correspondences.o build/calculus.o build/model.o
OBJS += build/raytracer.o  build/Camera.o
OBJS += build/INIReader.o build/Pose2d.o build/graph_def.o
OBJS += build/ini.o build/Icp.o build/InputData.o build/TermBuilder.o
OBJS += build/exporter_special.o build/cost.o
OBJS += build/GeometryUtility.o build/AnimationUtility.o build/Common.o
#build/InputData.o
SRC = src/

#CXXFLAGS += -ggp
#CCFLAGS += -ggp

FLAGS += -I/usr/include/eigen3 -I/usr/include/fbxsdk -I/home/dl/projects/nanoflann/include
FLAGS += -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include

LFLAGS += -L/usr/lib/eigen3 -L/usr/lib/gcc4/x64/release -losdCPU 
LFLAGS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_cpp -lhdf5
#-lfbxsdk
include make_template.inc
