#!/bin/bash

# ECRAD Optimized Build Script
# This script builds ECRAD with enhanced performance optimizations
# Based on performance analysis recommendations

set -e

# Default values
COMPILER="auto"
BUILD_TYPE="Release"
ENABLE_GPU="auto"
ENABLE_PROFILING="OFF"
ENABLE_MEMORY_OPT="ON"
ENABLE_LTO="ON"
ENABLE_VECTORIZATION="ON"
ENABLE_FAST_MATH="OFF"
INSTALL_PREFIX=""
VERBOSE="OFF"
CLEAN_BUILD="OFF"
NUM_JOBS=$(nproc)

# Help function
show_help() {
    echo "ECRAD Optimized Build Script"
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -c, --compiler COMPILER     Compiler to use (gcc, intel, nvhpc, auto) [default: auto]"
    echo "  -t, --build-type TYPE       Build type (Release, Debug, RelWithDebInfo) [default: Release]"
    echo "  -g, --gpu MODE              GPU support (on, off, auto) [default: auto]"
    echo "  -p, --profiling             Enable performance profiling"
    echo "  -m, --memory-opt            Enable memory optimizations [default: ON]"
    echo "  -l, --lto                   Enable link-time optimization [default: ON]"
    echo "  -v, --vectorization         Enable vectorization optimizations [default: ON]"
    echo "  -f, --fast-math             Enable fast math optimizations [default: OFF]"
    echo "  -i, --install-prefix PATH   Installation prefix"
    echo "  -j, --jobs N                Number of parallel jobs [default: $(nproc)]"
    echo "  -V, --verbose               Verbose output"
    echo "  -C, --clean                 Clean build directory first"
    echo "  -h, --help                  Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --compiler intel --gpu on --fast-math"
    echo "  $0 --build-type Debug --profiling --clean"
    echo "  $0 --compiler gcc --lto --vectorization --jobs 8"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--compiler)
            COMPILER="$2"
            shift 2
            ;;
        -t|--build-type)
            BUILD_TYPE="$2"
            shift 2
            ;;
        -g|--gpu)
            ENABLE_GPU="$2"
            shift 2
            ;;
        -p|--profiling)
            ENABLE_PROFILING="ON"
            shift
            ;;
        -m|--memory-opt)
            ENABLE_MEMORY_OPT="ON"
            shift
            ;;
        -l|--lto)
            ENABLE_LTO="ON"
            shift
            ;;
        -v|--vectorization)
            ENABLE_VECTORIZATION="ON"
            shift
            ;;
        -f|--fast-math)
            ENABLE_FAST_MATH="ON"
            shift
            ;;
        -i|--install-prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        -j|--jobs)
            NUM_JOBS="$2"
            shift 2
            ;;
        -V|--verbose)
            VERBOSE="ON"
            shift
            ;;
        -C|--clean)
            CLEAN_BUILD="ON"
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Detect compiler if auto
detect_compiler() {
    if [ "$COMPILER" = "auto" ]; then
        if command -v ifort &> /dev/null; then
            COMPILER="intel"
        elif command -v nvfortran &> /dev/null; then
            COMPILER="nvhpc"
        elif command -v gfortran &> /dev/null; then
            COMPILER="gcc"
        else
            log_error "No suitable Fortran compiler found"
            exit 1
        fi
    fi
    log_info "Using compiler: $COMPILER"
}

# Detect GPU support if auto
detect_gpu() {
    if [ "$ENABLE_GPU" = "auto" ]; then
        if command -v nvidia-smi &> /dev/null; then
            ENABLE_GPU="on"
            log_info "NVIDIA GPU detected, enabling GPU support"
        elif command -v rocm-smi &> /dev/null; then
            ENABLE_GPU="on"
            log_info "AMD GPU detected, enabling GPU support"
        else
            ENABLE_GPU="off"
            log_info "No GPU detected, disabling GPU support"
        fi
    fi
}

# Set up build directory
setup_build_dir() {
    BUILD_DIR="build_optimized_${COMPILER}_${BUILD_TYPE,,}"
    
    if [ "$CLEAN_BUILD" = "ON" ] && [ -d "$BUILD_DIR" ]; then
        log_info "Cleaning build directory: $BUILD_DIR"
        rm -rf "$BUILD_DIR"
    fi
    
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    log_info "Building in: $BUILD_DIR"
}

# Set up environment based on compiler
setup_environment() {
    case $COMPILER in
        intel)
            if [ -f "/opt/intel/oneapi/setvars.sh" ]; then
                source /opt/intel/oneapi/setvars.sh
            fi
            export FC=ifort
            export CC=icc
            ;;
        nvhpc)
            export FC=nvfortran
            export CC=nvc
            ;;
        gcc)
            export FC=gfortran
            export CC=gcc
            ;;
        *)
            log_error "Unsupported compiler: $COMPILER"
            exit 1
            ;;
    esac
}

# Generate CMake configuration
generate_cmake_config() {
    CMAKE_ARGS=()
    
    # Basic configuration
    CMAKE_ARGS+=("-DCMAKE_BUILD_TYPE=$BUILD_TYPE")
    CMAKE_ARGS+=("-DCMAKE_Fortran_COMPILER=$FC")
    CMAKE_ARGS+=("-DCMAKE_C_COMPILER=$CC")
    
    # Installation prefix
    if [ -n "$INSTALL_PREFIX" ]; then
        CMAKE_ARGS+=("-DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX")
    fi
    
    # GPU support
    if [ "$ENABLE_GPU" = "on" ]; then
        CMAKE_ARGS+=("-DECRAD_ENABLE_GPU=ON")
        CMAKE_ARGS+=("-DECRAD_ENABLE_ACC=ON")
    else
        CMAKE_ARGS+=("-DECRAD_ENABLE_GPU=OFF")
        CMAKE_ARGS+=("-DECRAD_ENABLE_ACC=OFF")
    fi
    
    # OpenMP support
    CMAKE_ARGS+=("-DECRAD_ENABLE_OMP=ON")
    
    # Precision options
    CMAKE_ARGS+=("-DECRAD_ENABLE_DOUBLE_PRECISION=ON")
    CMAKE_ARGS+=("-DECRAD_ENABLE_SINGLE_PRECISION=OFF")
    
    # Performance optimizations
    if [ "$ENABLE_PROFILING" = "ON" ]; then
        CMAKE_ARGS+=("-DENABLE_PERFORMANCE_PROFILING=ON")
    fi
    
    if [ "$ENABLE_MEMORY_OPT" = "ON" ]; then
        CMAKE_ARGS+=("-DENABLE_MEMORY_OPTIMIZATION=ON")
    fi
    
    # Use optimized compiler flags
    CMAKE_ARGS+=("-DECRAD_ECBUILD_COMPILE_FLAGS=../cmake/ecrad_compile_flags_optimized.cmake")
    
    # Verbose output
    if [ "$VERBOSE" = "ON" ]; then
        CMAKE_ARGS+=("-DCMAKE_VERBOSE_MAKEFILE=ON")
    fi
    
    log_info "CMake configuration:"
    for arg in "${CMAKE_ARGS[@]}"; do
        log_info "  $arg"
    done
}

# Run CMake configuration
configure_build() {
    log_info "Configuring build with CMake..."
    
    if ! cmake "${CMAKE_ARGS[@]}" ..; then
        log_error "CMake configuration failed"
        exit 1
    fi
    
    log_success "CMake configuration completed"
}

# Build the project
build_project() {
    log_info "Building ECRAD with $NUM_JOBS parallel jobs..."
    
    if [ "$VERBOSE" = "ON" ]; then
        BUILD_ARGS=("-j$NUM_JOBS" "VERBOSE=1")
    else
        BUILD_ARGS=("-j$NUM_JOBS")
    fi
    
    if ! make "${BUILD_ARGS[@]}"; then
        log_error "Build failed"
        exit 1
    fi
    
    log_success "Build completed successfully"
}

# Run tests if available
run_tests() {
    if [ "$BUILD_TYPE" != "Debug" ]; then
        log_info "Running tests..."
        
        if make test; then
            log_success "All tests passed"
        else
            log_warning "Some tests failed"
        fi
    else
        log_info "Skipping tests in debug mode"
    fi
}

# Install if requested
install_project() {
    if [ -n "$INSTALL_PREFIX" ]; then
        log_info "Installing to: $INSTALL_PREFIX"
        
        if make install; then
            log_success "Installation completed"
        else
            log_error "Installation failed"
            exit 1
        fi
    fi
}

# Print build summary
print_summary() {
    echo
    log_success "Build Summary:"
    log_info "  Compiler: $COMPILER"
    log_info "  Build Type: $BUILD_TYPE"
    log_info "  GPU Support: $ENABLE_GPU"
    log_info "  Profiling: $ENABLE_PROFILING"
    log_info "  Memory Optimization: $ENABLE_MEMORY_OPT"
    log_info "  LTO: $ENABLE_LTO"
    log_info "  Vectorization: $ENABLE_VECTORIZATION"
    log_info "  Fast Math: $ENABLE_FAST_MATH"
    log_info "  Build Directory: $BUILD_DIR"
    if [ -n "$INSTALL_PREFIX" ]; then
        log_info "  Install Prefix: $INSTALL_PREFIX"
    fi
    echo
    
    if [ -f "bin/ecrad" ]; then
        log_success "ECRAD executable: $(pwd)/bin/ecrad"
    fi
    
    if [ -f "lib/libradiation.a" ]; then
        log_success "Radiation library: $(pwd)/lib/libradiation.a"
    fi
    
    echo
    log_info "Performance recommendations:"
    log_info "  1. Use Release build type for production"
    log_info "  2. Enable GPU support if available"
    log_info "  3. Consider enabling fast math for non-critical applications"
    log_info "  4. Profile your application to identify bottlenecks"
    log_info "  5. Use appropriate number of OpenMP threads"
    echo
}

# Main execution
main() {
    log_info "Starting ECRAD optimized build..."
    
    # Setup
    detect_compiler
    detect_gpu
    setup_build_dir
    setup_environment
    
    # Configure and build
    generate_cmake_config
    configure_build
    build_project
    
    # Optional steps
    run_tests
    install_project
    
    # Summary
    print_summary
    
    log_success "ECRAD optimized build completed successfully!"
}

# Run main function
main "$@"