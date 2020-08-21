#!/bin/bash

./sweFom
python3 make_hyper_basis.py
./sweRom
python3 postProcess.py
