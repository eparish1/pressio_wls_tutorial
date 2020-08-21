# pressio_wls_tutorial

## Prerequisites

In order to successfully complete the tutorial, you need:

* Environment variables `CC`, `CXX` and `FC` to point to your serial compilers.
This has been tested wiht GCC-8.4.0

* CMake >= 3.11.0;

* Bash >= 3.2.57.

* Python>3.6 with numpy, scipy and matplotlib


## Building

(a) To make things easier/cleaner, create the follwoing environment variable:
```bash
export TUT_WORKDIR=<somewhere-where-you-want-to-build-tutorial>
```

From the top-level directory of this repo, do:
```bash
./do_build.sh -working-dir=${TUT_WORKDIR}
```

If the process succeeds, you should have the a build subdirectory inside
the ${TUT_WORKDIR} which contains the tutorial executables.


## Running
Proceed as follows:
```bash
cd ${TUT_WORKDIR}/build
./wfrun.sh
```
this takes a couple minutes and should generate a plot: `result.png'.
