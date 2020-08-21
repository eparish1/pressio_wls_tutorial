#!/bin/bash

OMP_PLACES=cores; OMP_PROC_BIND=true; ./sweFom
python make_hyper_basis.py
./sweRom
python postProcess.py
