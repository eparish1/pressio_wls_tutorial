#!/usr/bin/env bash

function build_eigen(){
    local PWD=`pwd`
    local PARENTDIR=$PWD

    # get eigen
    if [ ! -d eigen ]; then
  EIGENVERSION=3.3.7
  EIGENURL=https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
  EIGENUNPACKEDDIRNAME=eigen-3.3.7
  mkdir -p ${WORKINGDIR}/tpls/eigen
	cd ${WORKINGDIR}/tpls/eigen
	wget ${EIGENURL}
	tar xf eigen-${EIGENVERSION}.tar.gz
	mv ${EIGENUNPACKEDDIRNAME} eigen

	cd ${WORKINGDIR}
    fi

    echo "${fggreen}Eigen installed successfully${fgrst}"

    cd ${PARENTDIR}
}
