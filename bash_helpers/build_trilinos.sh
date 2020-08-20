#!/usr/bin/env bash

function build_trilinos(){
    local PWD=`pwd`
    local PARENTDIR=$PWD
    local TPLname=trilinos
    local nCoreMake=6
    #-----------------------------------

    # create dir
    [[ ! -d trilinos ]] && mkdir trilinos
    cd trilinos

    # clone repo
    if [ ! -d Trilinos ]; then
	git clone https://github.com/trilinos/Trilinos.git
	#git clone git@github.com:trilinos/Trilinos.git
	cd Trilinos
	git checkout master
	git checkout d5abc894f1052682c4933b9b6951d904e74aa7fe
	cd ..
    fi

    local PFX=${PWD}/install
    local BLASPFX=${PWD}/../openblas/install
    local LAPACKPFX=${PWD}/../openblas/install

    # create build
    if [ ! -d build ]; then
	mkdir build
    fi
    cd build

    CMAKELINE=""
    CMAKELINE+="-D CMAKE_VERBOSE_MAKEFILE=TRUE "
    CMAKELINE+="-D Trilinos_VERBOSE_CONFIGURE=FALSE "
    CMAKELINE+="-D CMAKE_INSTALL_PREFIX:PATH=${PFX} "
    CMAKELINE+="-D CMAKE_BUILD_TYPE:STRING=Release "
    CMAKELINE+="-D BUILD_SHARED_LIBS=On "
    CMAKELINE+="-D Trilinos_CXX11_FLAGS=-std=c++11 "

    CMAKELINE+="-D TPL_ENABLE_MPI:BOOL=OFF "
    CMAKELINE+="-D CMAKE_C_COMPILER=$CC "
    CMAKELINE+="-D CMAKE_CXX_COMPILER=$CXX "
    CMAKELINE+="-D CMAKE_Fortran_COMPILER=$FC "
    # CMAKELINE+="-D TPL_ENABLE_MPI:BOOL=OFF "
    # CMAKELINE+="-D MPI_C_COMPILER:FILEPATH=${CC} "
    # CMAKELINE+="-D MPI_CXX_COMPILER:FILEPATH=${CXX} "
    # CMAKELINE+="-D Trilinos_ENABLE_Fortran:BOOL=ON "
    # CMAKELINE+="-D MPI_Fortran_COMPILER:FILEPATH=${FC} "
    # CMAKELINE+="-D MPI_EXEC:FILEPATH=mpirun "
    # CMAKELINE+="-D MPI_USE_COMPILER_WRAPPERS:BOOL=ON "

    CMAKELINE+="-D Trilinos_ENABLE_ALL_PACKAGES=OFF "
    CMAKELINE+="-D Trilinos_ENABLE_TESTS=OFF "
    CMAKELINE+="-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=OFF "
    CMAKELINE+="-D Trilinos_ENABLE_Teuchos:BOOL=ON "
    CMAKELINE+="-D Trilinos_ENABLE_Tpetra:BOOL=ON "
    CMAKELINE+="-D Tpetra_ENABLE_DEPRECATED_CODE:BOOL=OFF "
    #CMAKELINE+="-D Tpetra_ENABLE_TSQR:BOOL=ON "
    #CMAKELINE+="-D Trilinos_ENABLE_AztecOO:BOOL=ON "
    #CMAKELINE+="-D Trilinos_ENABLE_Ifpack:BOOL=ON "
    #CMAKELINE+="-D Trilinos_ENABLE_Ifpack2:BOOL=ON "
    #CMAKELINE+="-D Trilinos_ENABLE_ROL:BOOL=ON "

    CMAKELINE+="-D Kokkos_ENABLE_SERIAL:BOOL=ON "
    CMAKELINE+="-D Kokkos_ENABLE_OPENMP:BOOL=OFF "
    CMAKELINE+="-D Kokkos_ENABLE_DEPRECATED_CODE=OFF "
    CMAKELINE+="-D TPL_ENABLE_BLAS=ON "
    CMAKELINE+="-D BLAS_LIBRARY_DIRS:PATH=${BLASPFX}/lib "
    CMAKELINE+="-D BLAS_LIBRARY_NAMES=openblas "
    CMAKELINE+="-D TPL_ENABLE_LAPACK=ON "
    CMAKELINE+="-D LAPACK_LIBRARY_DIRS:PATH=${LAPACKPFX}/lib "
    CMAKELINE+="-D LAPACK_LIBRARY_NAMES=openblas "
    CMAKELINE+="-D Trilinos_EXTRA_LINK_FLAGS:STRING=-ldl "
    CMAKELINE+="../Trilinos"

    echo "${fgyellow}Configuring ${TPLname} ${fgrst}"
    CFName="${PWD}/../config.txt"
    cmake eval ${CMAKELINE} > ${CFName} 2>&1
    echo "Config output written to ${CFName}"

    echo "${fgyellow}Building ${TPLname} ${fgrst}"
    BFName="${PWD}/../build.txt"
    make -j ${nCoreMake} > ${BFName} 2>&1
    echo "Build output written to ${BFName}"

    echo "${fgyellow}Installing ${TPLname} ${fgrst}"
    IFName="${PWD}/../install.txt"
    make install > ${IFName} 2>&1
    echo "Install output written to ${IFName}"

    cd ${PARENTDIR}
}
