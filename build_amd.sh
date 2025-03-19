#!/bin/bash -l
#sh build.sh -r=8 -rocm=5.4.3 -cuda=12.5 -offload=omp
staging=
cmake=
arch=gfx90a
MPI=OFF
for i in "$@"; do
    case "$1" in
	-r=*|--ranks=*)
	    ranks="${i#*=}"
	    shift # past argument=value
	    ;;
	-arch=*|--arch=*)
	    arch="${i#*=}"
	    shift # past argument=value
	    ;;	
	-staging|--staging)
   	    staging=YES
	    shift # past argument=value
	    ;;
	-cmake|--cmake)
            cmake=YES
	    shift # past argument=value
	    ;;
	-branch=*|--branch=*)
	    branch="${i#*=}"
	    shift # past argument=value
	    ;;
	-mpi|--mpi)
   	    MPI=ON
	    shift # past argument=value
	    ;;
	--)
	    shift
	    break
	    ;;
    esac
done
OMP=ON
ACC=OFF

if [[ $clear_cmakecache ]];
then
    if test -f ./CMakeCache.txt; then
	rm CMakeCache.txt
    fi
fi

module load cmake/3.27.4 curl/7.61.1
if [[ $staging ]]; then
    module load amd-staging
else
    module load rocm-afar-drop/6.1.0
fi
module list

if [[ $staging ]]; then
    BUILD_DIR=build.amd-staging.Release.${arch}
else
    BUILD_DIR=build.rocm-afar-6.1.0.Release.${arch}
fi

if [[ $MPI == "ON" ]]; then
    mpi_ext=""
else
    mpi_ext="-serial"
fi

if [[ $cmake ]]; then
    rm ${BUILD_DIR}/CMakeCache.txt
    FIAT_INSTALL_PATH=${HOME}/ecmwf/install/fiat/develop/
    if [[ $staging ]]; then
	#  -DDEBUG_CORRECTNESS_RADIATION
	cmake -B ${BUILD_DIR} -S . -DENABLE_BITIDENTITY_TESTING=ON -DCMAKE_BUILD_TYPE=RELEASE -DENABLE_ACC=${ACC} -DENABLE_OMP=${OMP} -DENABLE_GPU=ON -DENABLE_ROCTX=OFF -DECBUILD_Fortran_FLAGS="-O3 --offload-arch=${arch} -fopenmp -fPIC -DDEBUG_WARNING -DWORKAROUND_NESTED_PARALLEL_MCICA_OMP_LW -DWORKAROUND_NESTED_PARALLEL_MCICA_OMP_SW -DWORKAROUND_NESTED_PARALLEL_LW_DERIVATIVES -DWORKAROUND_NESTED_PARALLEL_FLUX" -DENABLE_SINGLE_PRECISION=OFF -DENABLE_DOUBLE_PRECISION=ON -DNETCDF_ROOT=$HOME/ecmwf/install/netcdf/amd-staging -Decbuild_ROOT=${ECBUILD_PATH} -Dfiat_ROOT=${FIAT_INSTALL_PATH}/amd-staging${mpi_ext}/ -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_EXE_LINKER_FLAGS_INIT="-lflang_rt.hostdevice --offload-arch=${arch} -fopenmp -O3" -DBUILD_SHARED_LIBS=OFF
    else
	#  -DDEBUG_CORRECTNESS_RADIATION
	cmake -B ${BUILD_DIR} -S . -DENABLE_BITIDENTITY_TESTING=ON -DCMAKE_BUILD_TYPE=RELEASE -DENABLE_ACC=${ACC} -DENABLE_OMP=${OMP} -DENABLE_GPU=ON -DENABLE_ROCTX=ON -DECBUILD_Fortran_FLAGS="-O3 --offload-arch=${arch} -fopenmp -DDEBUG_WARNING -DWORKAROUND_NESTED_PARALLEL_MCICA_OMP_LW -DWORKAROUND_NESTED_PARALLEL_MCICA_OMP_SW -DWORKAROUND_NESTED_PARALLEL_LW_DERIVATIVES -DWORKAROUND_NESTED_PARALLEL_FLUX -UDEBUG_CORRECTNESS_RADIATION" -DENABLE_SINGLE_PRECISION=ON -DENABLE_DOUBLE_PRECISION=ON -DNETCDF_ROOT=$HOME/ecmwf/install/netcdf/rocm-afar-6.1.0 -Decbuild_ROOT=${ECBUILD_PATH} -Dfiat_ROOT=${FIAT_INSTALL_PATH}/rocm-afar-6.1.0${mpi_ext}/ -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_EXE_LINKER_FLAGS_INIT="-lflang_rt.hostdevice --offload-arch=${arch} -fopenmp -O3" -DBUILD_SHARED_LIBS=OFF
    fi
fi
cmake --build ${BUILD_DIR} -- -j ${ranks}
