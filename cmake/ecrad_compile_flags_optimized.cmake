# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# OPTIMIZED COMPILER FLAGS FOR ECRAD PERFORMANCE
# This file contains enhanced compiler flags for better performance

# Capture ecbuild defaults and/or flags set by a toolchain
set( ${PNAME}_Fortran_FLAGS "${${PNAME}_Fortran_FLAGS} ${ECBUILD_Fortran_FLAGS}" )
set( ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} ${ECBUILD_Fortran_FLAGS_BIT}" )
set( ${PNAME}_Fortran_FLAGS_DEBUG "${${PNAME}_Fortran_FLAGS_DEBUG} ${ECBUILD_Fortran_FLAGS_DEBUG}" )

# Enhanced optimization flags for different compilers
if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set(checkbounds_flags   "-Rb")
  set(fpe_flags           "-Ktrap=fp")
  set(initsnan_flags      "-ei")
  set(enhanced_flags      "-O3 -hflex_mp=conservative -haggress -hnopattern")
  set(vectorization_flags "-hvector3 -hfp3")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(checkbounds_flags   "-fcheck=bounds")
  set(fpe_flags           "-ffpe-trap=invalid,zero,overflow")
  set(initsnan_flags      "-finit-real=snan")
  # Enhanced GCC optimization flags
  set(enhanced_flags      "-O3 -march=native -mtune=native -funroll-loops")
  set(vectorization_flags "-ftree-vectorize -fvect-cost-model=unlimited -ftree-loop-vectorize")
  set(fast_math_flags     "-ffast-math -fno-math-errno -ffinite-math-only")
  set(lto_flags           "-flto=auto -ffat-lto-objects")
  set(cache_flags         "-floop-nest-optimize -ftree-loop-distribution")
  set(simd_flags          "-fopenmp-simd")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(align_flags         "-align array64byte")
  set(alloc_flags         "-heap-arrays 32")
  set(checkbounds_flags   "-check bounds")
  set(initsnan_flags      "-init=snan")
  set(inline_flags        "-finline-functions -finline-limit=1500 -Winline")
  set(vectorization_flags "-assume byterecl,realloc_lhs")
  set(fpmodel_flags       "-fpe0 -fp-model precise -fp-speculation=safe -ftz")
  set(transcendentals_flags "-fast-transcendentals")
  # Enhanced Intel optimization flags
  set(enhanced_flags      "-O3 -xHost -ipo -no-prec-div")
  set(fast_math_flags     "-fp-model fast=2")
  set(parallel_flags      "-parallel -par-report3")
  set(vec_report_flags    "-vec-report3 -opt-report3")
  set(lto_flags           "-ipo")
  set(cache_flags         "-opt-prefetch -opt-mem-layout-trans")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  set(fpe_flags           "-Ktrap=fp")
  set(vectorization_flags "-O3 -fast")
  string(REPLACE "-O2" "" ${PNAME}_Fortran_FLAGS_BIT ${${PNAME}_Fortran_FLAGS_BIT})
  set(checkbounds_flags   "-Mbounds")
  # Enhanced PGI/NVHPC optimization flags
  set(enhanced_flags      "-O3 -fast -Munroll=c:4 -Mvect=simd")
  set(gpu_flags           "-acc -gpu=managed,lineinfo")
  set(parallel_flags      "-mp -Minfo=mp")
  set(cache_flags         "-Mcache_align -Mprefetch")

endif()

# Performance-oriented flag combinations
ecbuild_add_fortran_flags( "-g -O0"   NAME base_debug BUILD DEBUG )

# Apply enhanced flags in release builds
if(CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo")
  if( DEFINED enhanced_flags )
    ecbuild_add_fortran_flags( "${enhanced_flags}" NAME enhanced_opt BUILD BIT )
  endif()
  if( DEFINED fast_math_flags )
    ecbuild_add_fortran_flags( "${fast_math_flags}" NAME fast_math BUILD BIT )
  endif()
  if( DEFINED lto_flags )
    ecbuild_add_fortran_flags( "${lto_flags}" NAME lto BUILD BIT )
  endif()
  if( DEFINED cache_flags )
    ecbuild_add_fortran_flags( "${cache_flags}" NAME cache_opt BUILD BIT )
  endif()
  if( DEFINED parallel_flags )
    ecbuild_add_fortran_flags( "${parallel_flags}" NAME parallel_opt BUILD BIT )
  endif()
  if( DEFINED simd_flags )
    ecbuild_add_fortran_flags( "${simd_flags}" NAME simd_opt BUILD BIT )
  endif()
endif()

# Standard flags (unchanged from original)
if( DEFINED align_flags )
  ecbuild_add_fortran_flags( "${align_flags}"   NAME align )
endif()
if( DEFINED alloc_flags )
  ecbuild_add_fortran_flags( "${alloc_flags}"   NAME alloc )
endif()
if( DEFINED convert_flags )
  ecbuild_add_fortran_flags( "${convert_flags}"   NAME convert )
endif()
if( DEFINED vectorization_flags )
  # vectorization flags must be per-sourcefile overrideable, so are set via ${PNAME}_Fortran_FLAGS
  set( ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} ${vectorization_flags}" )
endif()
if( DEFINED fpmodel_flags )
  ecbuild_add_fortran_flags( "${fpmodel_flags}"   NAME fpmodel BUILD BIT )
endif()
if( DEFINED transcendentals_flags )
  ecbuild_add_fortran_flags( "${transcendentals_flags}"   NAME transcendentals BUILD BIT )
endif()
if( DEFINED inline_flags )
  ecbuild_add_fortran_flags( "${inline_flags}"   NAME inline BUILD BIT )
endif()

# GPU-specific optimizations
if( DEFINED gpu_flags AND HAVE_ACC )
  ecbuild_add_fortran_flags( "${gpu_flags}" NAME gpu_opt BUILD BIT )
endif()

# Compiler-specific reporting flags for performance analysis
if( DEFINED vec_report_flags )
  ecbuild_add_fortran_flags( "${vec_report_flags}" NAME vec_report BUILD BIT )
endif()

# Debug flags (unchanged from original)
if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  foreach( debug_flag    fpe initsnan checkbounds )
    if( ${debug_flag}_flags )
      set( ${PNAME}_Fortran_FLAGS_DEBUG "${${PNAME}_Fortran_FLAGS_DEBUG} ${${debug_flag}_flags}" )
    endif()
  endforeach()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # In case '-check all' has been added, we need to remove the '-check arg_temp_created' warnings
    set( ${PNAME}_Fortran_FLAGS_DEBUG "${${PNAME}_Fortran_FLAGS_DEBUG} -check noarg_temp_created" )
  endif()
endif()

# Performance validation flags
if( ENABLE_PERFORMANCE_PROFILING )
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set( ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} -prof-gen=srcpos" )
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set( ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} -pg -fprofile-arcs" )
  endif()
endif()

# Memory optimization flags
if( ENABLE_MEMORY_OPTIMIZATION )
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set( ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} -opt-mem-layout-trans=4" )
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set( ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} -floop-nest-optimize" )
  endif()
endif()

# Print optimization summary
message(STATUS "ECRAD Optimization Summary:")
message(STATUS "  Compiler: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo")
  message(STATUS "  Enhanced Optimizations: ENABLED")
  if( DEFINED enhanced_flags )
    message(STATUS "    Enhanced flags: ${enhanced_flags}")
  endif()
  if( DEFINED lto_flags )
    message(STATUS "    LTO flags: ${lto_flags}")
  endif()
  if( DEFINED vectorization_flags )
    message(STATUS "    Vectorization flags: ${vectorization_flags}")
  endif()
else()
  message(STATUS "  Enhanced Optimizations: DISABLED (Debug build)")
endif()