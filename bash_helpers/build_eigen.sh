#!/usr/bin/env bash

function build_eigen(){
    local PWD=`pwd`
    local PARENTDIR=$PWD

    # get eigen
    if [ ! -d eigen ]; then
	EIGENVERSION=3.3.7
	#EIGENTARURL=https://gitlab.com/libeigen/eigen/-/archive/${EIGENVERSION}/eigen-${EIGENVERSION}.tar.gz
	EIGENUNPACKEDDIRNAME=eigen-${EIGENVERSION}

	mkdir -p ${WORKINGDIR}/tpls/eigen
	cd ${WORKINGDIR}/tpls/eigen
	#wget ${EIGENTARURL}
	cp ${TOPDIR}/tpls/eigen-${EIGENVERSION}.tar.gz .
	tar zxf eigen-${EIGENVERSION}.tar.gz
	mv ${EIGENUNPACKEDDIRNAME} eigen
	cd ${WORKINGDIR}
    fi

    echo "${fggreen}Eigen installed successfully${fgrst}"

    cd ${PARENTDIR}
}
