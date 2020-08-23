#!/usr/bin/env bash

function version {
    echo "$@" | awk -F. '{ printf("%d%02d%02d\n", $1,$2,$3); }';
}

function detect_cmake(){
    # check if cmake is found in PATH
    if [ -x "$(command -v cmake)" ]; then
	CMAKEVERSIONDETECTED="$(cmake --version | perl -pe '($_)=/([0-9]+([.][0-9]+)+)/')"
    fi
}

function have_admissible_cmake(){
    # detect cmake
    detect_cmake
    echo ""
    echo "${fgyellow}+++ Detecting which cmake you have +++${fgrst}"
    echo "CMAKEVERSIONDETECTED=$CMAKEVERSIONDETECTED"

    if [ -z $CMAKEVERSIONDETECTED ]; then
	# if cmakeversion is empty
	echo "${fgred}No Cmake found. Terminating. ${fgrst}"
	exit 15
    else
	# if cmakeversion is NOT empty but wrong version
    	if [ $(version $CMAKEVERSIONDETECTED) -lt $(version "$CMAKEVERSIONMIN") ]; then
	    echo "${fgred}You have cmake ${CMAKEVERSIONDETECTED} ${fgrst}"
	    echo "${fgred}while I need >=${CMAKEVERSIONMIN}. Terminate.${fgrst}"
    	    exit 15
    	else
	    return 0
    	fi
    fi
}
