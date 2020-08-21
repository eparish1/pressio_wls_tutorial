#!/usr/bin/env bash

function build_pressio(){
    local PWD=`pwd`
    local PARENTDIR=$PWD

    # get pressio
    if [ ! -d pressio ]; then
	mkdir -p pressio
	cd pressio
	cp ${TOPDIR}/tpls/pressio.tar.gz .
	tar zxf pressio.tar.gz
	cd ..
	#cd pressio
	#git clone https://github.com/Pressio/pressio.git
	#git clone git@github.com:Pressio/pressio.git
	#cd ..
    fi
    cd pressio

    #cd pressio
    #git checkout rom/wls/hyperreduction
    #cd ..

    [[ ! -d build ]] && mkdir build

    cd build
    cmake -D CMAKE_INSTALL_PREFIX=../install \
	  -D PRESSIO_ENABLE_TPL_MPI=Off \
	  -D PRESSIO_ENABLE_TPL_TRILINOS=On \
	  -D PRESSIO_ENABLE_TPL_KOKKOS=On \
	  -D PRESSIO_ENABLE_TPL_BLAS=On \
	  -D PRESSIO_ENABLE_TPL_LAPACK=On \
	  ../pressio > ../config.txt

    make install > ../install.txt
    echo "${fggreen}Pressio installed successfully${fgrst}"

    cd ${PARENTDIR}
}
