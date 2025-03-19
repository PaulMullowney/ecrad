#!/bin/bash -l
#sh build.sh -r=8 -rocm=5.4.3 -cuda=12.5 -offload=omp
device="A100"
MPI=OFF
cmake=
for i in "$@"; do
    case "$1" in
	-r=*|--ranks=*)
	    ranks="${i#*=}"
	    shift # past argument=value
	    ;;
	-device=*|--device=*)
	    device="${i#*=}"
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
	-cmake|--cmake)
            cmake=YES
	    shift # past argument=value
	    ;;
	--)
	    shift
	    break
	    ;;
    esac
done

#module load openmpi/4.1.5-ucx1.13.1-cuda12.0 CUDA/${cuda}
#module list
module load cmake/3.27.4 curl/7.61.1
module list

if [[ $device == "A100" ]]; then
    arch="80"
    hpctk_ver=24.9
    hpcdr_ver=12.6
fi
if [[ $device == "H100" ]]; then
    arch="90"
    hpctk_ver=24.5
    hpcdr_ver=12.4
fi

NVHPC_CUDA_HOME=/opt/nvidia/hpc_sdk
TARGET=Linux_x86_64
VERSION=${hpctk_ver}
DRIVER_VER=${hpcdr_ver}
OMP_VER=openmpi4

CUDA_DIR=${NVHPC_CUDA_HOME}/${TARGET}/${VERSION}/cuda
nvcompdir=${NVHPC_CUDA_HOME}/${TARGET}/${VERSION}/compilers
nvmathdir=${NVHPC_CUDA_HOME}/${TARGET}/${VERSION}/math_libs
nvcommdir=${NVHPC_CUDA_HOME}/${TARGET}/${VERSION}/comm_libs/${DRIVER_VER}
nvompidir=${nvcommdir}/${OMP_VER}
# nsys_profilers_dir=${NVHPC_CUDA_HOME}/${TARGET}/${VERSION}/profilers/Nsight_Systems
# ncu_profilers_dir=${NVHPC_CUDA_HOME}/${TARGET}/${VERSION}/profilers/Nsight_Compute

export NVHPC=${NVHPC_CUDA_HOME}
export CC=${nvcompdir}/bin/nvc
export CXX=${nvcompdir}/bin/nvc++
export FC=${nvcompdir}/bin/nvfortran
export F90=${nvcompdir}/bin/nvfortran
export F77=${nvcompdir}/bin/nvfortran
export NVHPC_CMAKE=${NVHPC_CUDA_HOME}/${TARGET}/${VERSION}/cmake

export OPAL_PREFIX=${nvcommdir}/hpcx/latest/ompi
# export OPAL_PREFIX=$nvcommdir/openmpi4/openmpi-4.0.5
export PATH=${nvcompdir}/bin:${PATH}
export PATH=${CUDA_DIR}/bin:${PATH}
export PATH=${OPAL_PREFIX}/bin:${PATH}
# export PATH=${nvompidir}/bin:${PATH}
#export PATH=${nsys_profilers_dir}/bin:${PATH}
#export PATH=${ncu_profilers_dir}:${PATH}

export LD_LIBRARY_PATH=${nvcompdir}/lib/:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${CUDA_DIR}/lib64/:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${nvmathdir}/lib64/:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${OPAL_PREFIX}/lib/:${LD_LIBRARY_PATH}
# export LD_LIBRARY_PATH=${nvompidir}/lib/:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${nvcommdir}/nccl/lib/:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${nvcommdir}/nvshmem/lib:${LD_LIBRARY_PATH}

export LIBRARY_PATH=${nvcompdir}/lib/:${LIBRARY_PATH}
export LIBRARY_PATH=${CUDA_DIR}/lib64/:${LIBRARY_PATH}
export LIBRARY_PATH=${nvmathdir}/lib64/:${LIBRARY_PATH}
export LIBRARY_PATH=${OPAL_PREFIX}/lib/:${LIBRARY_PATH}
# export LIBRARY_PATH=${nvompidir}/lib/:${LIBRARY_PATH}
export LIBRARY_PATH=${nvcommdir}/nccl/lib/:${LIBRARY_PATH}
export LIBRARY_PATH=${nvcommdir}/nvshmem/lib:${LIBRARY_PATH}

# export NVHPC_CUDA_HOME=${CUDA_DIR}
export FFTW3_ROOT=${HOME}/workspace/FFT/FFTW3_INS
export C_INCLUDE_PATH=${CUDA_DIR}/include:${C_INCLUDE_PATH}
export C_INCLUDE_PATH=${nvmathdir}/include:${C_INCLUDE_PATH}

export CPLUS_INCLUDE_PATH=${CUDA_DIR}/include:${CPLUS_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${nvmathdir}/include:${CPLUS_INCLUDE_PATH}

export CPATH=${CUDA_DIR}/include:${CPATH}
export CPATH=${nvmathdir}/include:${CPATH}

export INCLUDE=${CUDA_DIR}/include:${INCLUDE}
export INCLUDE=${nvmathdir}/include:${INCLUDE}

#export CUBLAS_INCLUDE_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/math_libs/12.4/targets/x86_64-linux/
#export CUBLAS_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/math_libs/12.4/targets/x86_64-linux/

#export CUDA_cublas_ROOT=$CUDA_DIR/math_libs/12.4/targets/x86_64-linux/
#export CUDA_cublas_LIBRARY=${CUDA_cublas_ROOT}/lib/
#export CUDA_cufft_ROOT=$CUDA_DIR/math_libs/12.4/targets/x86_64-linux/
#export CUDA_cufft_LIBRARY=${CUDA_cufft_ROOT}/lib/

FIAT_INSTALL_PATH=${HOME}/ecmwf/install/fiat/develop/
ECBUILD_PATH=${HOME}/ecmwf/ecbuild

if [[ $MPI == "ON" ]]; then
    mpi_ext=""
else
    mpi_ext="-serial"
fi

if [[ $cmake == "YES" ]]; then
    if [[ $device == "A100" ]]; then
	cmake -B build.nv.${hpctk_ver}.Release -S . -DENABLE_BITIDENTITY_TESTING=ON -DCMAKE_BUILD_TYPE=Release -DENABLE_ACC=ON -DENABLE_OMP=OFF -DENABLE_GPU=ON -DENABLE_NVTX=ON -DENABLE_SINGLE_PRECISION=OFF -DENABLE_DOUBLE_PRECISION=ON -DECBUILD_Fortran_FLAGS="-UDEBUG_CORRECTNESS_RADIATION" -DOpenACC_Fortran_FLAGS="-acc=gpu -gpu=cc80,fastmath" -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_Fortran_COMPILER=nvfortran -DNETCDF_ROOT=$HOME/ecmwf/install/netcdf/${hpctk_ver} -Decbuild_ROOT=${ECBUILD_PATH} -Dfiat_ROOT=${FIAT_INSTALL_PATH}/nvhpc-${hpctk_ver}${mpi_ext}/
    fi

    if [[ $device == "H100" ]]; then
	cmake -B build.nv.${hpctk_ver}.Release -S . -DENABLE_BITIDENTITY_TESTING=ON -DCMAKE_BUILD_TYPE=Release -DENABLE_ACC=ON -DENABLE_OMP=OFF -DENABLE_GPU=ON -DENABLE_NVTX=ON -DENABLE_SINGLE_PRECISION=OFF -DENABLE_DOUBLE_PRECISION=ON -DECBUILD_Fortran_FLAGS="-UDEBUG_CORRECTNESS_RADIATION" -DOpenACC_Fortran_FLAGS="-acc=gpu -gpu=cc90,fastmath" -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_Fortran_COMPILER=nvfortran -DNETCDF_ROOT=$HOME/ecmwf/install/netcdf/${hpctk_ver} -Decbuild_ROOT=${ECBUILD_PATH} -Dfiat_ROOT=${FIAT_INSTALL_PATH}/nvhpc-${hpctk_ver}${mpi_ext}/
    fi
fi
cmake --build build.nv.${hpctk_ver}.Release -- -j ${ranks}
