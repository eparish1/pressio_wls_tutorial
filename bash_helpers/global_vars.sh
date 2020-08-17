#!/bin/bash

# top dir where this lives or is sourced
TOPDIR=${PWD}

# source
CPPSRC=${TOPDIR}

# the working dir
WORKINGDIR=

# yes/no wipe existing data content of target directory
WIPEEXISTING=no

# var to detect the os type [linux or mac]
ARCH=
if [[ $OSTYPE == *"darwin"* ]]; then
    ARCH=mac
else
    ARCH=linux
fi

function print_global_vars(){
    echo "CPPSRC         = $CPPSRC"
    echo "WORKINGDIR     = $WORKINGDIR"
    echo "WIPEEXISTING   = ${WIPEEXISTING}"
    echo "ARCH           = $ARCH"
}

function check_minimum_vars_set(){
    if [[ -z $WORKINGDIR ]]; then
	echo "--working-dir is empty, must be set: exiting"
	exit 1
    fi
}
