# pressio_wls_tutorial

## Building

In order to build the tutorial, you need:

* Environment variables `CC`, `CXX` and `FC` to point to your serial compilers.
This has been tested wiht GCC-8.4.0

* CMake >= 3.11.0;

* Bash >= 3.2.57.

(a) To make things easier and cleaner, create the follwoing environment variable:
```bash
export TUT_WORKDIR=<somewhere-where-you-want-to-build-tutorial>
```

From the top-level directory, execute:
```bash
./do_build.sh -working-dir=${TUT_WORKDIR}
```

If the process succeeds, you should have the a build subdirectory inside the ${TUT_WORKDIR} which contains the tutorial executables.

## Running

Proceed as follows:
```bash
cd ${TUT_WORKDIR}/build
./sweFom
```
this takes some time to run, and should create the following files:
```bash
solution0.bin
solution1.bin
...
solution9.bin
```
then generate the POD modes by running
```bash
python make_basis_hyper.py
```
and then run the WLS ROM with
```bash
./sweRom
```
and view the results by running
```bash
postProcess.py
```
