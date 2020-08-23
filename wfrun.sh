#!/bin/bash

function version {
    echo "$@" | awk -F. '{ printf("%d%02d%02d\n", $1,$2,$3); }';
}

function have_required_python_package(){
    echo ""
    echo "${fgyellow}--- Detecting if Python supports required package $1 ---${fgrst}"
    if python -c "import $1" &> /dev/null; then
	echo "${fggreen}Found $1 installed. ${fgrst}"
    else
	echo "${fgred} You don't have $1 installed. Please install it. Terminating. ${fgrst}"
	exit 11
    fi
}

if [ -x "$(command -v python3)" ]; then
    PYVERS="$(python3 --version | perl -pe '($_)=/([0-9]+([.][0-9]+)+)/')"
    echo "Found Python with version ${PYVERS}"
    if [ $(version $PYVERS) -lt $(version "3.6.10") ]; then
    	echo "${fgred}You have Python ${PYVERS} ${fgrst}"
    	echo "${fgred}while I need >=3.6.10. Terminate.${fgrst}"
    	exit 15
    fi

    have_required_python_package matplotlib
    have_required_python_package numpy
    have_required_python_package scipy
else
    echo "Python not found, terminating"
    exit 16
fi

# run everything
OMP_PLACES=cores; OMP_PROC_BIND=true; ./sweFom
python3 make_hyper_basis.py
./sweRom
python3 postProcess.py
