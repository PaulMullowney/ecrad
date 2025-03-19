#!/usr/bin/env bash

version=4.9.2

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export NETCDF_INSTALL_DIR=$(pwd)/netcdf-install
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export NETCDF_INSTALL_DIR="$2"; shift
        ;;
    "--tmpdir")
        TEMPORARY_FILES="$2"; shift
        ;;
    "--version")
        version="$2"; shift
        ;;
    "--hdf5-root")
        export HDF5_ROOT="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

NETCDF_MIRROR=https://downloads.unidata.ucar.edu/netcdf-c/
NETCDF_VERSION=${version}

URL=${NETCDF_MIRROR}/${NETCDF_VERSION}/netcdf-c-${NETCDF_VERSION}.tar.gz
FOLDER=netcdf-c-${NETCDF_VERSION}

if [ ! -d "${TEMPORARY_FILES}/${FOLDER}" ]; then
  echo "Downloading ${TEMPORARY_FILES}/${FOLDER} from URL [${URL}]"
  mkdir -p ${TEMPORARY_FILES}
  curl --location \
       "${URL}" | tar zx -C "${TEMPORARY_FILES}"
else
   echo "Download already present in ${TEMPORARY_FILES}/${FOLDER}"
fi

mkdir -p ${TEMPORARY_FILES}/build-${FOLDER} && cd ${TEMPORARY_FILES}/build-${FOLDER}
rm -rf ./*
cmake ${TEMPORARY_FILES}/${FOLDER} \
    -DHDF5_DIR=${HDF5_ROOT}/cmake -DCMAKE_INSTALL_PREFIX="${NETCDF_INSTALL_DIR}" \
    -DENABLE_TESTS=OFF  -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC -DCURL_LIBRARY="-lcurl" -DCMAKE_BUILD_TYPE=RELEASE
cmake --build . --config Release
cmake --install .
