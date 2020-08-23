FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
	build-essential \
	cmake	\
	git \
	wget \
	bash \
	ncurses-bin \
	ncurses-base \
	gfortran \
	python3.8 \
	zlib1g-dev \
	libncurses5-dev \
	python3-pip

RUN pip3 install --upgrade pip

#ENV TERM=xterm-256color
ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++
ENV FC=/usr/bin/gfortran
RUN ${CC} --version

ENV TUT_WORKDIR=/wlsTutorial
ENV TUTTPLS=${TUT_WORKDIR}/tpls
RUN mkdir -p ${TUTTPLS}

RUN mkdir -p ${TUTTPLS}/openblas \
	&& cd ${TUTTPLS}/openblas \
	&& wget https://github.com/xianyi/OpenBLAS/archive/v0.3.10.tar.gz \
	&& tar zxf v0.3.10.tar.gz \
	&& cd OpenBLAS-0.3.10 \
	&& make BINARY=64 HOSTCC=${CC} \
	&& make PREFIX=${PWD}/../install install

# EIGEN
ENV EIGENVERSION=3.3.7
ENV EIGENURL=https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
ENV EIGENUNPACKEDDIRNAME=eigen-3.3.7
RUN mkdir -p ${TUTTPLS}/eigen \
	&& cd ${TUTTPLS}/eigen \
	&& wget ${EIGENURL} \
	&& tar xf eigen-${EIGENVERSION}.tar.gz \
	&& mv ${EIGENUNPACKEDDIRNAME} eigen

# PRESSIO
RUN mkdir -p ${TUTTPLS}/pressio \
	&& cd ${TUTTPLS}/pressio \
	&& git clone https://github.com/Pressio/pressio.git \
	&& cd pressio \
	&& git checkout rom/wls/hyperreduction \
	&& cd .. \
	&& mkdir build && cd build \
	&& cmake -D CMAKE_INSTALL_PREFIX=../install \
	-D PRESSIO_ENABLE_TPL_MPI=Off \
	-D PRESSIO_ENABLE_TPL_TRILINOS=On \
	-D PRESSIO_ENABLE_TPL_KOKKOS=On \
	-D PRESSIO_ENABLE_TPL_BLAS=On \
	-D PRESSIO_ENABLE_TPL_LAPACK=On \
	../pressio \
	&& make install

# TRILINOS
ENV TRILPFX=${TUTTPLS}/trilinos/install
ENV BLASPFX=${TUTTPLS}/openblas/install
ENV LAPACKPFX=${TUTTPLS}/openblas/install

RUN mkdir -p ${TUTTPLS}/trilinos \
	&& cd ${TUTTPLS}/trilinos \
	&& git clone https://github.com/trilinos/Trilinos.git \
	&& cd Trilinos \
	&& git checkout d5abc894f1052682c4933b9b6951d904e74aa7fe \
	&& cd .. \
	&& mkdir build && cd build \
	&& cmake \
	-D CMAKE_VERBOSE_MAKEFILE=TRUE \
	-D Trilinos_VERBOSE_CONFIGURE=FALSE \
	-D CMAKE_INSTALL_PREFIX:PATH=${TRILPFX} \
	-D CMAKE_BUILD_TYPE:STRING=Release \
	-D BUILD_SHARED_LIBS=On \
	-D Trilinos_CXX11_FLAGS=-std=c++11 \
	-D TPL_ENABLE_MPI:BOOL=OFF \
	-D CMAKE_C_COMPILER=$CC \
	-D CMAKE_CXX_COMPILER=$CXX \
	-D CMAKE_Fortran_COMPILER=$FC \
	-D Trilinos_ENABLE_ALL_PACKAGES=OFF \
	-D Trilinos_ENABLE_TESTS=OFF \
	-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=OFF \
	-D Trilinos_ENABLE_Teuchos:BOOL=ON \
	-D Trilinos_ENABLE_Tpetra:BOOL=ON \
	-D Tpetra_ENABLE_DEPRECATED_CODE:BOOL=OFF \
	-D Kokkos_ENABLE_SERIAL:BOOL=ON \
	-D Kokkos_ENABLE_OPENMP:BOOL=OFF \
	-D Kokkos_ENABLE_DEPRECATED_CODE=OFF \
	-D TPL_ENABLE_BLAS=ON \
	-D BLAS_LIBRARY_DIRS:PATH=${BLASPFX}/lib \
	-D BLAS_LIBRARY_NAMES=openblas \
	-D TPL_ENABLE_LAPACK=ON \
	-D LAPACK_LIBRARY_DIRS:PATH=${LAPACKPFX}/lib \
	-D LAPACK_LIBRARY_NAMES=openblas \
	-D Trilinos_EXTRA_LINK_FLAGS:STRING=-ldl \
	../Trilinos \
	&& make -j 4 \
	&& make install
