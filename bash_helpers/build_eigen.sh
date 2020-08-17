#!/usr/bin/env bash

function build_eigen(){
    local PWD=`pwd`
    local PARENTDIR=$PWD

    # get eigen
    if [ ! -d eigen ]; then
	EIGENVERSION=3.3.7
	EIGENURL=http://bitbucket.org/eigen/eigen/get/${EIGENVERSION}
	EIGENUNPACKEDDIRNAME=eigen-eigen-323c052e1731

	mkdir -p ${WORKINGDIR}/tpls/eigen
	cd ${WORKINGDIR}/tpls/eigen
	wget ${EIGENURL}.tar.bz2
	tar xf ${EIGENVERSION}.tar.bz2
	mv ${EIGENUNPACKEDDIRNAME} eigen
	cd ${WORKINGDIR}
    fi

    echo "${fggreen}Eigen installed successfully${fgrst}"

    cd ${PARENTDIR}
}
