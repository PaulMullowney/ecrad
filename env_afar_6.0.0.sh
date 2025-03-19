#!/bin/bash -l
module load cmake/3.27.4 curl/7.61.1
module load rocm-afar-drop/5.3.0
module list

CC=$(which amdclang)
CXX=$(which amdclang++)
FC=$(which amdflang-new)
