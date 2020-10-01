FROM nvidia/cuda:9.2-devel
SHELL ["/bin/bash", "-c"]
ENV BASH_ENV "/root/.bashrc"
RUN apt-get update && apt-get install -y cmake g++ libeigen3-dev \
 libhdf5-dev libtinfo-dev  \
 libncurses5-dev libedit-dev zlib1g-dev \
 libtbb-dev nano gdb \
&& echo "export CC=/clang/install/bin/clang" >> /bashrc \
&& echo "export CXX=/clang/install/bin/clang++" >> /bashrc \
&& echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/clang/install/lib/:/terra/install/lib/" \
>> /bashrc
COPY --from=joms/clang:latest /clang/install /clang/install
COPY --from=joms/terra:latest /terra /terra
COPY --from=joms/optlang:latest /Optlang /Optlang
COPY --from=joms/json:latest /json/install /json/install
COPY --from=joms/nanoflann:latest /nanoflann/install /nanoflann/install
COPY --from=joms/cxxopts:latest /cxxopts/install /cxxopts/install
COPY --from=joms/fbxsdk:latest /fbxsdk /fbxsdk
COPY --from=joms/embree:latest /embree/install /embree/install
COPY --from=joms/opensubdiv:latest /opensubdiv/install /opensubdiv/install
ADD . /human-model-trainer
WORKDIR /human-model-trainer
RUN echo "alias llvm-config='/clang/install/bin/llvm-config'" >> /root/bashrc && \
echo "alias clang='/clang/install/bin/clang'" >> /root/bashrc && \
echo "alias clang++='/clang/install/bin/clang++'" >> /root/bashrc
RUN echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/opensubdiv/install/lib" >> bashrc && \
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/clang/install/lib" >> bashrc && \
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/terra/install/lib" >> bashrc && \
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/fbxsdk/lib/gcc4/x64/release" >> bashrc && \
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/embree/install/lib" >> bashrc && \
cat /bashrc >> /human-model-trainer/bashrc && \
cat /human-model-trainer/bashrc >> ~/.bashrc
RUN bash -c "source /human-model-trainer/bashrc && make -j16"
