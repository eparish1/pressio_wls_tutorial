#!/usr/bin/env bash

function build_openblas(){
    local PWD=`pwd`
    local PARENTDIR=$PWD

    #-----------------------------------
    # create dir
    [[ ! -d openblas ]] && mkdir openblas
    cd openblas

    if [ ! -f v0.3.10.tar.gz ]; then
	#wget https://github.com/xianyi/OpenBLAS/archive/v0.3.10.tar.gz
	cp ${TOPDIR}/tpls/v0.3.10.tar.gz .
    fi

    if [ ! -d OpenBLAS-0.3.10 ]; then
	tar zxf v0.3.10.tar.gz
    fi

    echo ${PWD}
    if [ ! -d install ]; then
	cd OpenBLAS-0.3.10

	BFName="${PWD}/../build.txt"
	IFName="${PWD}/../install.txt"

	echo "${fgyellow}Building OpenBLAS ${fgrst}"
	make BINARY=64 HOSTCC=$CC > ${BFName} 2>&1
	echo "Build output written to ${BFName}"
	if grep -q "OpenBLAS build complete" "${BFName}"; then
	    echo "${fggreen}OpenBLAS built successfully${fgrst}"
	else
	    echo "${fgred}OpenBLAS built unsuccessfully${fgrst}"
	    exit 44
	fi

	echo "${fgyellow}Installing OpenBLAS ${fgrst}"
	make PREFIX=${PWD}/../install install > ${IFName} 2>&1
	echo "Install output written to ${IFName}"
    fi

    cd ${PARENTDIR}
}
