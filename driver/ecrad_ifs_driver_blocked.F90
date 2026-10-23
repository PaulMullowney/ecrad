! ecrad_ifs_driver_blocked.F90 - Driver for offline ECRAD radiation scheme
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! ECRAD is the radiation scheme used in the ECMWF Integrated
! Forecasting System in cycle 43R3 and later. Several solvers are
! available, including McICA, Tripleclouds and SPARTACUS (the Speedy
! Algorithm for Radiative Transfer through Cloud Sides, a modification
! of the two-stream formulation of shortwave and longwave radiative
! transfer to account for 3D radiative effects). Gas optical
! properties are provided by the RRTM-G gas optics scheme.

! This program takes three arguments:
! 1) Namelist file to configure the radiation calculation, but note
!    that only the radiation_config group is read
! 2) Name of a NetCDF file containing one or more atmospheric profiles
! 3) Name of output NetCDF file
!
! This version uses the infrastructure of the IFS, such as computing
! effective radius and cloud overlap from latitude and other
! variables. To configure ecRad in this version you need to edit
! ifs/yoerad.F90 in the ecRad package, but these options can be
! overridden with the "radiation" namelist. This file requires the
! input data to have compatible settings, e.g. the right number of
! aerosol variables, and surface albedo/emissivity bands; a test file
! satisfying this requirement is test/ifs/ecrad_meridian.nc in the
! ecRad package.
!
! Note that the purpose of this file is simply to demonstrate the use
! of the setup_radiation_scheme and radiation_scheme routines as well
! as the use of a blocked memory layout to improve cache efficiency;
! all the rest is using the offline ecRad driver containers to read
! a NetCDF file to memory and pass it into these routines.

program ecrad_ifs_driver

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------
  use parkind1,                 only : jprb, jprd ! Working/double precision
  use yomhook,                  only : dr_hook_init
#ifdef HAVE_FIAT
  use mpl_module,               only : mpl_init, mpl_end
#endif

  use radiation_io,             only : nulout
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, IMassMixingRatio, &
       &   IH2O, ICO2, IO3, IN2O, INO2, ICO, ICH4, IO2, ICFC11, ICFC12, &
       &   IHCFC22, ICCl4
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use radiation_flux,           only : flux_type
  use radiation_save,           only : save_net_fluxes
  use radiation_setup,          only : tradiation, setup_radiation_scheme
  use radiation_constants,      only : Pi
  use ecrad_driver_config,      only : driver_config_type
  use ecrad_driver_read_input,  only : read_input
  use easy_netcdf
  use ifs_blocking
  use iso_fortran_env,          only : int64
#ifdef HAVE_NVTX
  use nvtx
#endif
#ifdef HAVE_ROCTX
  use roctx_profiling, only: roctxstartrange, roctxendrange
  use iso_c_binding, only: c_null_char
#endif

  implicit none

#include "radiation_scheme.intfb.h"

  ! The NetCDF file containing the input profiles
  type(netcdf_file)         :: file

  ! Configuration for the radiation scheme, IFS style
  type(tradiation)          :: yradiation

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol

  ! Configuration specific to this driver
  type(driver_config_type)  :: driver_config

  ! Derived type containing outputs from the radiation scheme
  type(flux_type)           :: flux

  ! Additional arrays passed to radiation_scheme
  real(jprb), allocatable, dimension(:) :: ccn_land, ccn_sea, sin_latitude, longitude_rad, land_frac
  real(jprb), allocatable, dimension(:,:) :: pressure_fl, temperature_fl
  real(jprb), allocatable, dimension(:) :: flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
       &  emissivity_out
  real(jprb), allocatable, dimension(:,:) :: flux_diffuse_band, flux_direct_band

  ! Bespoke data types to set-up the blocked memory layout
  type(ifs_config_type)        :: ifs_config
  real(kind=jprb), allocatable :: zrgp(:,:,:) ! monolithic IFS data structure
#ifdef BITIDENTITY_TESTING
  integer, allocatable         :: iseed(:,:) ! Seed for random number generator
#endif

  integer :: ncol, nlev         ! Number of columns and levels
  integer :: nproma             ! block size

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

#ifndef NO_OPENMP
  ! OpenMP functions
  integer, external :: omp_get_thread_num
  real(kind=jprd), external :: omp_get_wtime
  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop
#endif

  ! For demonstration of get_sw_weights later on
  ! Ultraviolet weightings
  !integer    :: nweight_uv
  !integer    :: iband_uv(100)
  !real(jprb) :: weight_uv(100)
  ! Photosynthetically active radiation weightings
  !integer    :: nweight_par
  !integer    :: iband_par(100)
  !real(jprb) :: weight_par(100)

  ! Loop index for repeats (for benchmarking)
  integer :: jrepeat

  integer(kind=int64) :: count_rate, t(4)
  real(kind=jprd) :: total_dt, kernel_dt

  ! Loop index
  integer :: jrl, ibeg, iend, il, ib

  ! If enabled, COPY_ASYNC stages copy in, compute and copy out of ZRGP
  ! in different streams, overlapping the copy-in of the next block with the
  ! compute of the current block, and the copy-out of the processed block
  ! with the compute of the next block. However, this does not seem to
  ! improve overall performance and it is therefore disabled by default.
! #define COPY_ASYNC
#ifdef COPY_ASYNC
  integer :: next_il
#endif

  ! Execution tile (OMPGPU only, enabled by ECRAD_FULL_BATCH_SIZE=K > 0):
  ! zrgp stays resident on the device with its (nproma,field,block) layout,
  ! and K consecutive blocks are gathered on the device into one contiguous
  ! (K*nproma,field) array so that each radiation_scheme call processes
  ! K*nproma columns. Each column keeps the seed of its original block.
  integer :: full_batch_size, full_width, full_active
  integer :: jblock, jcolumn, jfield, icol, iblk
  integer :: iinend, ioutbeg, ioutend, ifldstot, iaprs
  real(kind=jprb), allocatable :: full_fields(:,:)
  integer, allocatable :: full_seed(:)

  ! Multi-process partitioning (ECRAD_NPARTS=N, ECRAD_PART=r with 0 <= r < N):
  ! blocks are dealt round-robin in chunks of ECRAD_PART_CHUNK consecutive
  ! blocks, chunk c going to part mod(c,N), and this process computes only its
  ! own blocks. Run one process per GPU; driver/merge_parts.py must use the same
  ! N and chunk to recombine the outputs.
  integer :: nparts, ipart, part_chunk, nblocks, nown, il_last
  integer :: jtile, ntile_blocks
  integer, allocatable :: own_blocks(:)
  logical :: tile_overlap
  logical, allocatable :: is_own_block(:)

  ! Are any variables out of bounds?
  logical :: is_out_of_bounds

!  integer    :: iband(20), nweights
!  real(jprb) :: weight(20)

  ! Initialise MPI if not done yet
#ifdef HAVE_FIAT
  call mpl_init
#endif

  call dr_hook_init()
  call system_clock(count_rate=count_rate)

  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: ecrad config.nam input_file.nc output_file.nc'
  end if

  ! Use namelist to configure the radiation calculation
  call get_command_argument(1, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of namelist file as string of length < 512'
  end if

  ! Read "radiation_driver" namelist into radiation driver config type
  call driver_config%read(file_name)
  nproma = driver_config%nblocksize

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '-------------------------- OFFLINE ECRAD RADIATION SCHEME --------------------------'
    write(nulout,'(a)') 'Copyright (C) 2014- ECMWF'
    write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
#ifdef PARKIND1_SINGLE
    write(nulout,'(a)') 'Floating-point precision: single'
#else
    write(nulout,'(a)') 'Floating-point precision: double'
#endif
  end if

  ! Albedo/emissivity intervals may be specified like this
  !call config%define_sw_albedo_intervals(6, &
  !     &  [0.25e-6_jprb, 0.44e-6_jprb, 0.69e-6_jprb, &
  !     &     1.19_jprb, 2.38e-6_jprb], [1,2,3,4,5,6], &
  !     &   do_nearest=.false.)
  !call config%define_lw_emiss_intervals(3, &
  !     &  [8.0e-6_jprb, 13.0e-6_jprb], [1,2,1], &
  !     &   do_nearest=.false.)

  ! If monochromatic aerosol properties are required, then the
  ! wavelengths can be specified (in metres) as follows - these can be
  ! whatever you like for the general aerosol optics, but must match
  ! the monochromatic values in the aerosol input file for the older
  ! aerosol optics
  !call config%set_aerosol_wavelength_mono( &
  !     &  [3.4e-07_jprb, 3.55e-07_jprb, 3.8e-07_jprb, 4.0e-07_jprb, 4.4e-07_jprb, &
  !     &   4.69e-07_jprb, 5.0e-07_jprb, 5.32e-07_jprb, 5.5e-07_jprb, 6.45e-07_jprb, &
  !     &   6.7e-07_jprb, 8.0e-07_jprb, 8.58e-07_jprb, 8.65e-07_jprb, 1.02e-06_jprb, &
  !     &   1.064e-06_jprb, 1.24e-06_jprb, 1.64e-06_jprb, 2.13e-06_jprb, 1.0e-05_jprb])

  call yradiation%rad_config%read(file_name=file_name)

  ! Setup aerosols
  if (yradiation%rad_config%use_aerosols) then
    yradiation%yrerad%naermacc = 1 ! MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
  else
    yradiation%yrerad%naermacc = 0
  endif

  ! Setup the radiation scheme: load the coefficients for gas and
  ! cloud optics, currently from RRTMG
  call setup_radiation_scheme(yradiation, .true., file_name=file_name)
  ! Or call without specifying the namelist filename, in which case
  ! the default settings are from yoerad.F90
  !call setup_radiation_scheme(yradiation, .true.)

  ! Demonstration of how to get weights for UV and PAR fluxes
  !if (config%do_sw) then
  !  call config%get_sw_weights(0.2e-6_jprb, 0.4415e-6_jprb,&
  !       &  nweight_uv, iband_uv, weight_uv,&
  !       &  'ultraviolet')
  !  call config%get_sw_weights(0.4e-6_jprb, 0.7e-6_jprb,&
  !       &  nweight_par, iband_par, weight_par,&
  !       &  'photosynthetically active radiation, PAR')
  !end if

  ! --------------------------------------------------------
  ! Section 3: Read input data file
  ! --------------------------------------------------------

  ! Get NetCDF input file name
  call get_command_argument(2, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of input NetCDF file as string of length < 512'
  end if

  ! Open the file and configure the way it is read
  call file%open(trim(file_name), iverbose=driver_config%iverbose)

  ! Get NetCDF output file name
  call get_command_argument(3, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of output NetCDF file as string of length < 512'
  end if

  ! 2D arrays are assumed to be stored in the file with height varying
  ! more rapidly than column index. Specifying "true" here transposes
  ! all 2D arrays so that the column index varies fastest within the
  ! program.
  call file%transpose_matrices(.true.)

  ! Read input variables from NetCDF file, noting that cloud overlap
  ! and effective radius are ignored
  call read_input(file, yradiation%rad_config, driver_config, ncol, nlev, &
       &          single_level, thermodynamics, &
       &          gas, cloud, aerosol)

  ! Latitude is used for cloud overlap and ice effective radius
  if (file%exists('lat')) then
    call file%get('lat', sin_latitude)
    sin_latitude = sin(sin_latitude * Pi/180.0_jprb)
  else
    allocate(sin_latitude(ncol))
    sin_latitude = 0.0_jprb
  end if

  if (file%exists('lon')) then
    call file%get('lon', longitude_rad)
    longitude_rad = longitude_rad * Pi/180.0_jprb
  else
    allocate(longitude_rad(ncol))
    longitude_rad = 0.0_jprb
  end if

  ! Close input file
  call file%close()

  ! Convert gas units to mass-mixing ratio
  call gas%set_units(IMassMixingRatio, lacc=.false.)

  ! Compute seed from skin temperature residual
  !  single_level%iseed = int(1.0e9*(single_level%skin_temperature &
  !       &                            -int(single_level%skin_temperature)))

  ! Set first and last columns to process
  if (driver_config%iendcol < 1 .or. driver_config%iendcol > ncol) then
    driver_config%iendcol = ncol
  end if

  if (driver_config%istartcol > driver_config%iendcol) then
    write(nulout,'(a,i0,a,i0,a,i0,a)') '*** Error: requested column range (', &
         &  driver_config%istartcol, &
         &  ' to ', driver_config%iendcol, ') is out of the range in the data (1 to ', &
         &  ncol, ')'
    stop 1
  end if

  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  ! Compute saturation with respect to liquid (needed for aerosol
  ! hydration) call
  call thermodynamics%calc_saturation_wrt_liquid(driver_config%istartcol,driver_config%iendcol)

  ! Check inputs are within physical bounds, printing message if not
  is_out_of_bounds =     gas%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.   single_level%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or. thermodynamics%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.          cloud%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.        aerosol%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs)

  ! Allocate memory for the flux profiles, which may include arrays
  ! of dimension n_bands_sw/n_bands_lw, so must be called after
  ! setup_radiation
  call flux%allocate(yradiation%rad_config, 1, ncol, nlev)

  ! set relevant fluxes to zero
  flux%lw_up(:,:) = 0._jprb
  flux%lw_dn(:,:) = 0._jprb
  flux%sw_up(:,:) = 0._jprb
  flux%sw_dn(:,:) = 0._jprb
  flux%sw_dn_direct(:,:) = 0._jprb
  flux%lw_up_clear(:,:) = 0._jprb
  flux%lw_dn_clear(:,:) = 0._jprb
  flux%sw_up_clear(:,:) = 0._jprb
  flux%sw_dn_clear(:,:) = 0._jprb
  flux%sw_dn_direct_clear(:,:) = 0._jprb

  flux%lw_dn_surf_canopy(:,:) = 0._jprb
  flux%sw_dn_diffuse_surf_canopy(:,:) = 0._jprb
  flux%sw_dn_direct_surf_canopy(:,:) = 0._jprb
  flux%lw_derivatives(:,:) = 0._jprb

  ! Allocate memory for additional arrays
  allocate(ccn_land(ncol))
  allocate(ccn_sea(ncol))
  allocate(land_frac(ncol))
  allocate(pressure_fl(ncol,nlev))
  allocate(temperature_fl(ncol,nlev))
  allocate(flux_sw_direct_normal(ncol))
  allocate(flux_uv(ncol))
  allocate(flux_par(ncol))
  allocate(flux_par_clear(ncol))
  allocate(emissivity_out(ncol))
  allocate(flux_diffuse_band(ncol,yradiation%yrerad%nsw))
  allocate(flux_direct_band(ncol,yradiation%yrerad%nsw))

  ccn_land = yradiation%yrerad%rccnlnd
  ccn_sea = yradiation%yrerad%rccnsea
  pressure_fl = 0.5_jprb * (thermodynamics%pressure_hl(:,1:nlev)+thermodynamics%pressure_hl(:,2:nlev+1))
  temperature_fl = 0.5_jprb * (thermodynamics%temperature_hl(:,1:nlev)+thermodynamics%temperature_hl(:,2:nlev+1))

  ! --------------------------------------------------------
  ! Section 4a: Reshuffle into blocked memory layout
  ! --------------------------------------------------------

  call ifs_setup_indices(driver_config, ifs_config, yradiation, nlev)
  call ifs_copy_inputs_to_blocked(driver_config, ifs_config, yradiation,&
        & ncol, nlev, single_level, thermodynamics, gas, cloud, aerosol,&
        & sin_latitude, longitude_rad, land_frac, pressure_fl, temperature_fl,&
        & zrgp &
#ifdef BITIDENTITY_TESTING
        &, iseed=iseed &
#endif
        & )

  ! --------------------------------------------------------
  ! Section 4b: Call radiation_scheme with blocked memory data
  ! --------------------------------------------------------

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)')  'Performing radiative transfer calculations'
  end if

  full_batch_size = env_int('ECRAD_FULL_BATCH_SIZE', 0)
  if (full_batch_size < 0) stop 'ECRAD_FULL_BATCH_SIZE must be non-negative'
  tile_overlap = env_int('ECRAD_TILE_OVERLAP', 0) /= 0
#if !defined(OMPGPU) || defined(COPY_ASYNC)
  if (full_batch_size > 0) stop 'ECRAD_FULL_BATCH_SIZE requires an OMPGPU build without COPY_ASYNC'
#endif

  nparts     = env_int('ECRAD_NPARTS', 1)
  ipart      = env_int('ECRAD_PART', 0)
  part_chunk = env_int('ECRAD_PART_CHUNK', 8)
  if (nparts < 1 .or. ipart < 0 .or. ipart >= nparts .or. part_chunk < 1) then
    stop 'Need ECRAD_NPARTS >= 1, 0 <= ECRAD_PART < ECRAD_NPARTS and ECRAD_PART_CHUNK >= 1'
  end if
  nblocks = (ncol+nproma-1)/nproma
  il_last = ncol - (nblocks-1)*nproma
  allocate(is_own_block(nblocks), own_blocks(nblocks))
  nown = 0
  do ib = 1, nblocks
    is_own_block(ib) = mod((ib-1)/part_chunk, nparts) == ipart
    if (is_own_block(ib)) then
      nown = nown + 1
      own_blocks(nown) = ib
    end if
  end do
  if (nparts > 1) then
    write(nulout,'(a,i0,a,i0,a,i0,a,i0,a)') 'Part ', ipart, ' of ', nparts, ': ', &
         &  nown, ' of ', nblocks, ' blocks'
  end if

  if (nown == 0) full_batch_size = 0

  if (full_batch_size > 0) then
    full_width = min(full_batch_size, nown)*nproma
    iinend   = ifs_config%iinend
    ioutbeg  = ifs_config%ioutbeg
    ioutend  = ifs_config%ioutend
    ifldstot = ifs_config%ifldstot
    iaprs    = ifs_config%iaprs
    if (driver_config%iverbose >= 2) then
      write(nulout,'(a,i0,a,i0,a)') 'Execution tile: ', full_width/nproma, &
           &  ' blocks of ', nproma, ' columns per radiation_scheme call'
    end if
    allocate(full_fields(full_width,ifldstot))
    allocate(full_seed(full_width))
    full_fields = 0._jprb
    ! Partial last block: keep the padding rows defined, as the per-block path does
    if (il_last < nproma) then
      zrgp(il_last+1:nproma, :, nblocks) = 0._jprb
    end if
    !$OMP TARGET ENTER DATA MAP(ALLOC: zrgp, full_fields, full_seed)
    !$OMP TARGET ENTER DATA MAP(TO: own_blocks)
#ifdef BITIDENTITY_TESTING
    !$OMP TARGET ENTER DATA MAP(TO: iseed)
#endif
  end if

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  do jrepeat = 1,driver_config%nrepeat
     total_dt = 0.0_jprd
     kernel_dt = 0.0_jprd

#ifdef HAVE_NVTX
     call nvtxStartRange("ecrad_it")
#endif
#ifdef HAVE_ROCTX
    call roctxStartRange("ecrad_it"//c_null_char)
#endif

#ifdef COPY_ASYNC
#if defined(OMPGPU)
#error "NOT IMPLEMENTED!"
#endif
#if defined(_OPENACC)
    !$acc enter data create(zrgp(:,:,1)) &
#ifdef BITIDENTITY_TESTING
    !$acc&     copyin(iseed(:,1)) &
#endif
    !$acc&     async(2)

    next_il = min(nproma,ncol)
    !$acc update device(zrgp(:,ifs_config%iinbeg:ifs_config%iinend,1), &
    !$acc&              zrgp(:,ifs_config%ioutend+1:ifs_config%ifldstot,1)) async(2)
#endif
#endif !! COPY_ASYNC

#ifndef NO_OPENMP
    if (jrepeat == driver_config%nwarmup + 1) then
      tstart = omp_get_wtime()
    end if
#endif

!    if (driver_config%do_parallel) then
      ! Run radiation scheme over blocks of columns in parallel

#if defined(OMPGPU) && !defined(COPY_ASYNC)
    if (full_batch_size > 0) then
      call system_clock(count=t(1))
      ! Each tile moves only its own blocks in and its own outputs out. With
      ! overlap, tile k+1's blocks are copied in while tile k runs.
      if (tile_overlap) call update_tile(1, min(full_batch_size, nown), .false.)

      do jtile = 1, nown, full_batch_size
        ntile_blocks = min(full_batch_size, nown-jtile+1)
        if (.not. tile_overlap) then
          call update_tile(jtile, ntile_blocks, .false.)
        else if (jtile+full_batch_size <= nown) then
          call update_tile(jtile+full_batch_size, &
               &           min(full_batch_size, nown-jtile-full_batch_size+1), .true.)
        end if
        full_active = (ntile_blocks-1)*nproma
        if (own_blocks(jtile+ntile_blocks-1) == nblocks) then
          full_active = full_active + il_last
        else
          full_active = full_active + nproma
        end if

        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(jblock, jcolumn, iblk)
        do jfield = 1, ifldstot
          do icol = 1, full_width
            if (icol <= full_active .and. (jfield <= iinend .or. jfield > ioutend)) then
              jblock  = (icol-1)/nproma
              jcolumn = icol - jblock*nproma
              iblk    = own_blocks(jtile+jblock)
              full_fields(icol,jfield) = zrgp(jcolumn,jfield,iblk)
            else
              full_fields(icol,jfield) = 0._jprb
            end if
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(jblock, jcolumn, iblk)
        do icol = 1, full_width
          jblock  = (icol-1)/nproma
          jcolumn = icol - jblock*nproma
#ifdef BITIDENTITY_TESTING
          if (icol <= full_active) then
            iblk = own_blocks(jtile+jblock)
            full_seed(icol) = iseed(jcolumn,iblk)
          else
            full_seed(icol) = 0
          end if
#else
          ! Matches init_seed_simple applied to each original block
          full_seed(icol) = jcolumn
#endif
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

        ! radiation_scheme reads the top two half-level pressures on the host
        do jfield = iaprs, iaprs+1
          do icol = 1, full_active
            jblock  = (icol-1)/nproma
            jcolumn = icol - jblock*nproma
            full_fields(icol,jfield) = zrgp(jcolumn,jfield,own_blocks(jtile+jblock))
          end do
        end do

        call system_clock(count=t(3))
        call run_radiation_tile(full_width, full_active, full_fields, full_seed)
        call system_clock(count=t(4))
        kernel_dt = kernel_dt + (t(4)-t(3))/dble(count_rate)

        ! The tile's outputs are contiguous on the device: copy them back in one
        ! transfer and scatter them into the blocked host array
        !$OMP TARGET UPDATE FROM(full_fields(:,ioutbeg:ioutend))
        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(jblock, iblk, il)
        do jfield = ioutbeg, ioutend
          do jblock = 0, ntile_blocks-1
            iblk = own_blocks(jtile+jblock)
            il   = min(nproma, full_active - jblock*nproma)
            zrgp(1:il,jfield,iblk) = full_fields(jblock*nproma+1:jblock*nproma+il,jfield)
          end do
        end do
        !$OMP END PARALLEL DO

        if (tile_overlap) then
          ! Next tile's inputs
          !$OMP TASKWAIT
        end if
      end do
      call system_clock(count=t(2))
      total_dt = total_dt + (t(2)-t(1))/dble(count_rate)
    else
#endif

#ifndef OMPGPU
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) REDUCTION(+:total_dt,kernel_dt)&
      !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB,t)
#endif
      do jrl=1,ncol,nproma
        if (.not. is_own_block((jrl-1)/nproma+1)) cycle
        call system_clock(count=t(1))
        ibeg=jrl
        iend=min(ibeg+nproma-1,ncol)
        il=iend-ibeg+1
        ib=(jrl-1)/nproma+1

        if (driver_config%iverbose >= 3) then
#ifndef NO_OPENMP
          write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
               &  ' processing columns ', ibeg, '-', iend
#else
          write(nulout,'(a,i0,a,i0)')  'Processing columns ', ibeg, '-', iend
#endif
        end if

#ifdef COPY_ASYNC
        if (iend < ncol) then
#if defined(OMPGPU)
#endif
#if defined(_OPENACC)
          !$acc enter data create(zrgp(:,:,ib+1)) &
#ifdef BITIDENTITY_TESTING
          !$acc&     copyin(iseed(:,ib+1)) &
#endif
          !$acc&     async(2) wait(2)
#endif
          next_il = min(iend+nproma,ncol) - (ibeg+nproma) + 1
#if defined(OMPGPU)
#endif
#if defined(_OPENACC)
          !$acc update device(zrgp(:,ifs_config%iinbeg:ifs_config%iinend,ib+1), &
          !$acc&              zrgp(:,ifs_config%ioutend+1:ifs_config%ifldstot,ib+1)) async(2)
#endif
        endif

#else  /* COPY_ASYNC */

        ! Pad unused column rows with zero on host so full nproma-sized
        ! OMPGPU transfers are safe when il < nproma (partial last block).
#if defined(OMPGPU)
        if (il < nproma) then
          zrgp(il+1:nproma, ifs_config%iinbeg:ifs_config%iinend, ib) = 0._jprb
          zrgp(il+1:nproma, ifs_config%ioutend+1:ifs_config%ifldstot, ib) = 0._jprb
        endif
        !$OMP TARGET ENTER DATA MAP(ALLOC:zrgp(:,:,ib))
#ifdef BITIDENTITY_TESTING
        !$OMP TARGET UPDATE TO(iseed(:,ib))
#endif
        !$OMP TARGET UPDATE TO(zrgp(:,ifs_config%iinbeg:ifs_config%iinend,ib), &
        !$OMP&                 zrgp(:,ifs_config%ioutend+1:ifs_config%ifldstot,ib))
#endif
#if defined(_OPENACC)
        !$acc data create(zrgp(:,:,ib)) &
#ifdef BITIDENTITY_TESTING
        !$acc&     copyin(iseed(:,ib)) &
#endif
        !$acc&

        !$acc update device(zrgp(1:il,ifs_config%iinbeg:ifs_config%iinend,ib), &
        !$acc&            zrgp(1:il,ifs_config%ioutend+1:ifs_config%ifldstot,ib)) &
        !$acc& async(1)
#endif
#endif /* COPY_ASYNC */

        call system_clock(count=t(3))

        ! Call the ECRAD radiation scheme
        call radiation_scheme &
             & (yradiation, &
             &  1, il, nproma, &                       ! startcol, endcol, ncol
             &  nlev, size(aerosol%mixing_ratio,3), &    ! nlev, naerosols
             &  single_level%solar_irradiance, &                               ! solar_irrad
             ! array inputs
             &  zrgp(:,ifs_config%iamu0,ib), zrgp(:,ifs_config%its,ib), &    ! mu0, skintemp
             &  zrgp(:,ifs_config%iald,ib) , zrgp(:,ifs_config%ialp,ib), &    ! albedo_dif, albedo_dir
             &  zrgp(:,ifs_config%iemiss,ib), &                   ! spectral emissivity
             &  zrgp(:,ifs_config%iccnl,ib), zrgp(:,ifs_config%iccno,ib) ,&  ! CCN concentration, land and sea
             &  zrgp(:,ifs_config%igelam,ib),zrgp(:,ifs_config%igemu,ib), &  ! longitude, sine of latitude
             &  zrgp(:,ifs_config%islm,ib), &                     ! land sea mask
             &  zrgp(:,ifs_config%ipr,ib),   zrgp(:,ifs_config%iti,ib),  &   ! full level pressure and temperature
             &  zrgp(:,ifs_config%iaprs,ib), zrgp(:,ifs_config%ihti,ib), &   ! half-level pressure and temperature
             &  zrgp(:,ifs_config%iwv,ib),   zrgp(:,ifs_config%iico2,ib), &
             &  zrgp(:,ifs_config%iich4,ib), zrgp(:,ifs_config%iin2o,ib), &
             &  zrgp(:,ifs_config%ino2,ib),  zrgp(:,ifs_config%ic11,ib), &
             &  zrgp(:,ifs_config%ic12,ib),  zrgp(:,ifs_config%ic22,ib), &
             &  zrgp(:,ifs_config%icl4,ib),  zrgp(:,ifs_config%ioz,ib), &
             &  zrgp(:,ifs_config%iclc,ib),  zrgp(:,ifs_config%ilwa,ib), &
             &  zrgp(:,ifs_config%iiwa,ib),  zrgp(:,ifs_config%irwa,ib), &
             &  zrgp(:,ifs_config%iswa,ib), &
             &  zrgp(:,ifs_config%iaer,ib),  zrgp(:,ifs_config%iaero,ib), &
             ! flux outputs
             &  zrgp(:,ifs_config%ifrso,ib), zrgp(:,ifs_config%ifrth,ib), &
             &  zrgp(:,ifs_config%iswfc,ib), zrgp(:,ifs_config%ilwfc,ib),&
             &  zrgp(:,ifs_config%ifrsod,ib),zrgp(:,ifs_config%ifrted,ib), &
             &  zrgp(:,ifs_config%ifrsodc,ib),zrgp(:,ifs_config%ifrtedc,ib),&
             &  zrgp(:,ifs_config%ifdir,ib), zrgp(:,ifs_config%icdir,ib), &
             &  zrgp(:,ifs_config%isudu,ib), &
             &  zrgp(:,ifs_config%iuvdf,ib), zrgp(:,ifs_config%iparf,ib), &
             &  zrgp(:,ifs_config%iparcf,ib),zrgp(:,ifs_config%itincf,ib), &
             &  zrgp(:,ifs_config%iemit,ib) ,zrgp(:,ifs_config%ilwderivative,ib), &
             &  zrgp(:,ifs_config%iswdiffuseband,ib), zrgp(:,ifs_config%iswdirectband,ib)&
#ifdef BITIDENTITY_TESTING
            ! To validate results against standalone ecrad, we overwrite effective
            ! radii, cloud overlap and seed with input values
             &  ,pre_liq=zrgp(:,ifs_config%ire_liq,ib), &
             &  pre_ice=zrgp(:,ifs_config%ire_ice,ib), &
             &  pcloud_overlap=zrgp(:,ifs_config%ioverlap,ib), &
             &  iseed=iseed(:,ib) &
#endif
             & )

        call system_clock(count=t(4))

#if defined(OMPGPU)
#ifdef COPY_ASYNC
#else
        ! Flang omptarget mishandles partial row slices (1:il,:); use the
        ! same full nproma output-field transfer as full blocks. Only
        ! rows 1:il are read when scattering fluxes back to global arrays.
        !$OMP TARGET UPDATE FROM(zrgp(:,ifs_config%ioutbeg:ifs_config%ioutend,ib))
        !$OMP TARGET EXIT DATA MAP(DELETE:zrgp(:,:,ib))
#endif
#endif
#if defined(_OPENACC)
#ifdef COPY_ASYNC
          !$acc update host(zrgp(:,ifs_config%ioutbeg:ifs_config%ioutend,ib)) async(3) wait(1)
          !$acc exit data delete(zrgp(:,:,ib)) async(3)
#else
          !$acc update host(zrgp(1:il,ifs_config%ioutbeg:ifs_config%ioutend,ib)) async(1)
          !$acc end data
#endif
#endif
        call system_clock(count=t(2))
        total_dt = total_dt + (t(2)-t(1))/dble(count_rate)
        kernel_dt = kernel_dt + (t(4)-t(3))/dble(count_rate)
      end do
#ifndef OMPGPU
      !$OMP END PARALLEL DO
#endif
#if defined(OMPGPU) && !defined(COPY_ASYNC)
    end if
#endif
!    else
      ! Run radiation scheme serially
!      if (driver_config%iverbose >= 3) then
!        write(nulout,'(a,i0,a)')  'Processing ', ncol, ' columns'
!      end if

      ! Call the ECRAD radiation scheme
!      call radiation_scheme(ncol, nlev, driver_config%istartcol, driver_config%iendcol, &
!           &  config, single_level, thermodynamics, gas, cloud, aerosol, flux)

!    end if
#ifdef HAVE_NVTX
     call nvtxEndRange
#endif
#ifdef HAVE_ROCTX
     call roctxEndRange
#endif

     write(nulout, '(a,g12.5,a)') 'time elapsed in radiative transfer: ', &
          &                         total_dt, ' seconds'
     write(nulout, '(a,g12.5,a)') 'time elapsed in radiative transfer kernel: ', &
          &                         kernel_dt, ' seconds'
     write(nulout, '(a,i0)') 'Columns/s : ', int((ncol)/total_dt)
     write(nulout, '(a,i0)') 'Columns/s : ', int((ncol)/kernel_dt)

  end do

#if defined(OMPGPU)
#endif
#if defined(_OPENACC)
!$acc wait
#endif

  if (full_batch_size > 0) then
    !$OMP TARGET EXIT DATA MAP(DELETE: zrgp, full_fields, full_seed, own_blocks)
#ifdef BITIDENTITY_TESTING
    !$OMP TARGET EXIT DATA MAP(DELETE: iseed)
#endif
    deallocate(full_fields, full_seed)
  end if

  ! --------------------------------------------------------
  ! Section 4c: Copy fluxes from blocked memory data
  ! --------------------------------------------------------

  call ifs_copy_fluxes_from_blocked(driver_config, ifs_config, yradiation, ncol, nlev,&
          & zrgp, flux, flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
          & emissivity_out, flux_diffuse_band, flux_direct_band)

  ! "up" fluxes are actually net fluxes at this point - we modify the
  ! upwelling flux so that net=dn-up, while the TOA and surface
  ! downwelling fluxes are correct.
  flux%sw_up = -flux%sw_up
  flux%sw_up(:,1) = flux%sw_up(:,1)+flux%sw_dn(:,1)
  flux%sw_up(:,nlev+1) = flux%sw_up(:,nlev+1)+flux%sw_dn(:,nlev+1)

  flux%lw_up = -flux%lw_up
  flux%lw_up(:,1) = flux%lw_up(:,1)+flux%lw_dn(:,1)
  flux%lw_up(:,nlev+1) = flux%lw_up(:,nlev+1)+flux%lw_dn(:,nlev+1)

  flux%sw_up_clear = -flux%sw_up_clear
  flux%sw_up_clear(:,1) = flux%sw_up_clear(:,1)+flux%sw_dn_clear(:,1)
  flux%sw_up_clear(:,nlev+1) = flux%sw_up_clear(:,nlev+1)+flux%sw_dn_clear(:,nlev+1)

  flux%lw_up_clear = -flux%lw_up_clear
  flux%lw_up_clear(:,1) = flux%lw_up_clear(:,1)+flux%lw_dn_clear(:,1)
  flux%lw_up_clear(:,nlev+1) = flux%lw_up_clear(:,nlev+1)+flux%lw_dn_clear(:,nlev+1)

  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------

  ! This is unreliable because only the net fluxes are valid:
  !is_out_of_bounds = flux%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol)

  ! Store the fluxes in the output file
  yradiation%rad_config%do_surface_sw_spectral_flux = .false.
  yradiation%rad_config%do_canopy_fluxes_sw = .false.
  yradiation%rad_config%do_canopy_fluxes_lw = .false.

  call save_net_fluxes(file_name, yradiation%rad_config, thermodynamics, flux, &
       &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
       &   experiment_name=driver_config%experiment_name, &
       &   is_double_precision=driver_config%do_write_double_precision)

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '------------------------------------------------------------------------------------'
  end if

  ! Finalise MPI if not done yet
#ifdef HAVE_FIAT
  call mpl_end(ldmeminfo=.false.)
#endif

contains

  integer function env_int(name, default)
    character(len=*), intent(in) :: name
    integer, intent(in)          :: default
    character(len=32) :: value
    integer :: status

    env_int = default
    call get_environment_variable(name, value, status=status)
    if (status == 0 .and. len_trim(value) > 0) read(value,*) env_int
  end function env_int

  ! Copy the blocks own_blocks(jfirst:jfirst+nblk-1) to the device, one
  ! contiguous run of blocks at a time. Sections strided across blocks are
  ! mishandled by this runtime, so each run is copied whole.
  subroutine update_tile(jfirst, nblk, async)
    integer, intent(in) :: jfirst, nblk
    logical, intent(in) :: async
    integer :: jpos, jbeg, ib1, ib2

    jbeg = jfirst
    do jpos = jfirst, jfirst+nblk-1
      if (jpos == jfirst+nblk-1 .or. own_blocks(min(jpos+1,nown)) /= own_blocks(jpos)+1) then
        ib1 = own_blocks(jbeg)
        ib2 = own_blocks(jpos)
        if (async) then
          !$OMP TARGET UPDATE TO(zrgp(:,:,ib1:ib2)) NOWAIT
        else
          !$OMP TARGET UPDATE TO(zrgp(:,:,ib1:ib2))
        end if
        jbeg = jpos + 1
      end if
    end do
  end subroutine update_tile

  subroutine run_radiation_tile(kcol, kend, pfields, kseed)
    integer, intent(in)            :: kcol, kend
    real(kind=jprb), intent(inout) :: pfields(kcol, ifs_config%ifldstot)
    integer, intent(in)            :: kseed(kcol)

    call radiation_scheme &
         & (yradiation, 1, kend, kcol, nlev, size(aerosol%mixing_ratio,3), &
         &  single_level%solar_irradiance, &
         &  pfields(:,ifs_config%iamu0), pfields(:,ifs_config%its), &
         &  pfields(:,ifs_config%iald),  pfields(:,ifs_config%ialp), &
         &  pfields(:,ifs_config%iemiss), &
         &  pfields(:,ifs_config%iccnl), pfields(:,ifs_config%iccno), &
         &  pfields(:,ifs_config%igelam),pfields(:,ifs_config%igemu), &
         &  pfields(:,ifs_config%islm), &
         &  pfields(:,ifs_config%ipr),   pfields(:,ifs_config%iti), &
         &  pfields(:,ifs_config%iaprs), pfields(:,ifs_config%ihti), &
         &  pfields(:,ifs_config%iwv),   pfields(:,ifs_config%iico2), &
         &  pfields(:,ifs_config%iich4), pfields(:,ifs_config%iin2o), &
         &  pfields(:,ifs_config%ino2),  pfields(:,ifs_config%ic11), &
         &  pfields(:,ifs_config%ic12),  pfields(:,ifs_config%ic22), &
         &  pfields(:,ifs_config%icl4),  pfields(:,ifs_config%ioz), &
         &  pfields(:,ifs_config%iclc),  pfields(:,ifs_config%ilwa), &
         &  pfields(:,ifs_config%iiwa),  pfields(:,ifs_config%irwa), &
         &  pfields(:,ifs_config%iswa), &
         &  pfields(:,ifs_config%iaer),  pfields(:,ifs_config%iaero), &
         &  pfields(:,ifs_config%ifrso), pfields(:,ifs_config%ifrth), &
         &  pfields(:,ifs_config%iswfc), pfields(:,ifs_config%ilwfc), &
         &  pfields(:,ifs_config%ifrsod),pfields(:,ifs_config%ifrted), &
         &  pfields(:,ifs_config%ifrsodc),pfields(:,ifs_config%ifrtedc), &
         &  pfields(:,ifs_config%ifdir), pfields(:,ifs_config%icdir), &
         &  pfields(:,ifs_config%isudu), &
         &  pfields(:,ifs_config%iuvdf), pfields(:,ifs_config%iparf), &
         &  pfields(:,ifs_config%iparcf),pfields(:,ifs_config%itincf), &
         &  pfields(:,ifs_config%iemit), pfields(:,ifs_config%ilwderivative), &
         &  pfields(:,ifs_config%iswdiffuseband), pfields(:,ifs_config%iswdirectband) &
#ifdef BITIDENTITY_TESTING
         &  ,pre_liq=pfields(:,ifs_config%ire_liq), &
         &  pre_ice=pfields(:,ifs_config%ire_ice), &
         &  pcloud_overlap=pfields(:,ifs_config%ioverlap) &
#endif
         &  ,iseed=kseed)
  end subroutine run_radiation_tile

end program ecrad_ifs_driver
