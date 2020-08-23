#!/bin/bash

set -e

source ./bash_helpers/colors.sh
source ./bash_helpers/global_vars.sh
source ./bash_helpers/cmd_line_options.sh

# check that all basic variables are set, otherwise leave
check_minimum_vars_set

echo ""
echo "--------------------------------------------"
echo " current setting is: "
echo ""
print_global_vars
echo ""
echo "--------------------------------------------"

# detect valid cmake
source ./bash_helpers/detect_cmake.sh
have_admissible_cmake
echo "${fggreen}Valid cmake found: ok! ${fgrst}"
echo ""

# create working dir if not existing
[[ ! -d ${WORKINGDIR} ]] && mkdir -p ${WORKINGDIR}

# create tpls dir if not existing
[[ ! -d ${WORKINGDIR}/tpls ]] && mkdir -p ${WORKINGDIR}/tpls

# # wipe everything if set to 1
# [[ $WIPEEXISTING == yes ]] && wipe_existing_data_in_target_dir

#---------------------------
# go to working dir
cd ${WORKINGDIR}

# do tpls
cd tpls
source ${TOPDIR}/bash_helpers/build_openblas.sh && build_openblas
source ${TOPDIR}/bash_helpers/build_trilinos.sh && build_trilinos
source ${TOPDIR}/bash_helpers/build_eigen.sh && build_eigen
source ${TOPDIR}/bash_helpers/build_pressio.sh && build_pressio
cd ..

# do SWE
bdirname=build
[[ ! -d ${bdirname} ]] && mkdir ${bdirname}

# enter
rm -rf ${bdirname}/*
cd ${bdirname}
cmake -DCMAKE_CXX_COMPILER=${CXX} \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
      -DCMAKE_BUILD_TYPE=Debug \
      \
      -DEIGEN_INCLUDE_DIR=${WORKINGDIR}/tpls/eigen/eigen \
      -DPRESSIO_INCLUDE_DIR=${WORKINGDIR}/tpls/pressio/install/include \
      \
      -DOpenBLAS_DIR=${WORKINGDIR}/tpls/openblas/install/lib/cmake/openblas \
      -DTRILINOS_INC_DIR=${WORKINGDIR}/tpls/trilinos/install/include \
      -DTRILINOS_LIB_DIR=${WORKINGDIR}/tpls/trilinos/install/lib \
      \
      ${CPPSRC}
make
cd ..

# go back where we started
cd ${TOPDIR}
