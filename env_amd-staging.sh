#!/bin/bash -l
module load cmake/3.27.4 curl/7.61.1
module load amd-staging
module list

export CC=$(which amdclang)
export CXX=$(which amdclang++)
export FC=$(which amdflang)
