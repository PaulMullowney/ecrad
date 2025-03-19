! radiation_flux.F90 - Derived type to store the output fluxes
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
! Modifications
!   2017-09-08  R. Hogan  Store g-point fluxes
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2019-01-08  R. Hogan  Added "indexed_sum_profile"
!   2019-01-14  R. Hogan  out_of_physical_bounds calls routine in radiation_config
!   2021-01-20  R. Hogan  Added heating_rate_out_of_physical_bounds function
!   2022-12-07  R. Hogan  Added top-of-atmosphere spectral output

#include "ecrad_config.h"

module radiation_flux

  use parkind1,     only : jprb
  use radiation_io, only : nulout, nulerr

  implicit none
  public

  !---------------------------------------------------------------------
  ! This derived type contains the output from the radiation
  ! calculation.  Currently this is solely flux profiles, but in
  ! future surface fluxes in each band may be stored in order that the
  ! calling program can compute surface-radiation such as
  ! photosynthetically active radiation and UV index.
  type flux_type
     ! All the following are broad-band fluxes in W m-2 with
     ! dimensions (ncol,nlev+1).  Note that only those fluxes that are
     ! requested will be used, so clear-sky and direct-beam arrays may
     ! not be allocated
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up, lw_dn, &   ! Upwelling and downwelling longwave
          &  sw_up, sw_dn, &   ! Upwelling and downwelling shortwave
          &  sw_dn_direct, &   ! Direct-beam shortwave into a horizontal plane
          &  lw_up_clear, lw_dn_clear, & ! Clear-sky quantities...
          &  sw_up_clear, sw_dn_clear, &
          &  sw_dn_direct_clear
     ! As above but fluxes in each spectral band in W m-2 with
     ! dimensions (nband,ncol,nlev+1).  These are only allocated if
     ! config%do_save_spectral_flux==.true.
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  lw_up_band, lw_dn_band, &   ! Upwelling and downwelling longwave
          &  sw_up_band, sw_dn_band, &   ! Upwelling and downwelling shortwave
          &  sw_dn_direct_band, &        ! Direct-beam shortwave
          &  lw_up_clear_band, lw_dn_clear_band, & ! Clear-sky quantities...
          &  sw_up_clear_band, sw_dn_clear_band, &
          &  sw_dn_direct_clear_band
     ! Surface downwelling quantities at each g point, dimensioned
     ! (ng,ncol), that are always saved by the solver, except for the
     ! clear-sky ones that are only produced if
     ! config%do_clear==.true.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_surf_g, lw_dn_surf_clear_g, &
          &  sw_dn_diffuse_surf_g, sw_dn_direct_surf_g, &
          &  sw_dn_diffuse_surf_clear_g, sw_dn_direct_surf_clear_g
     ! Top-of-atmosphere quantities at each g point, dimensioned
     ! (ng,ncol), that are always saved by the solver, except for the
     ! clear-sky ones that are only produced if
     ! config%do_clear==.true.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up_toa_g, lw_up_toa_clear_g, &
          &  sw_dn_toa_g, sw_up_toa_g, sw_up_toa_clear_g
     ! Shortwave downwelling spectral fluxes in W m-2 at the surface,
     ! from which quantities such as photosynthetically active and UV
     ! radiation can be computed. Only allocated if
     ! config%do_surface_sw_spectral_flux==.true.  Note that the
     ! clear-sky quantities are only computed if
     ! config%do_clear==.true., but direct fluxes are computed whether
     ! or not do_direct==.true.. The dimensions are (nband,ncol).
     real(jprb), allocatable, dimension(:,:) :: &
          &  sw_dn_surf_band, sw_dn_direct_surf_band, &
          &  sw_dn_surf_clear_band, sw_dn_direct_surf_clear_band
     ! Top-of-atmosphere spectral fluxes in W m-2. Only allocated if
     ! config%do_toa_spectral_flux=.true.. Note that the clear-sky
     ! quantities are only computed if config%do_clear==.true.. The
     ! dimensions are (nband,ncol).
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up_toa_band, lw_up_toa_clear_band, &
          &  sw_dn_toa_band, sw_up_toa_band, sw_up_toa_clear_band
     ! Surface downwelling fluxes in W m-2 at the spectral resolution
     ! needed by any subsequent canopy radiative transfer.  If
     ! config%use_canopy_full_spectrum_[sw|lw] then these will be at
     ! g-point resolution; otherwise they will be at
     ! config%n_albedo_bands and config%n_emiss_bands resolution.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_surf_canopy, &
          &  sw_dn_diffuse_surf_canopy, sw_dn_direct_surf_canopy

     ! Diagnosed cloud cover from the short- and long-wave solvers
     real(jprb), allocatable, dimension(:) :: &
          &  cloud_cover_lw, cloud_cover_sw
     ! Longwave derivatives needed by Hogan and Bozzo (2015) method
     ! for approximate longwave updates in between the full radiation
     ! calls: rate of change of upwelling broad-band flux with respect
     ! to surface value, dimensioned (ncol,nlev+1)
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_derivatives

   contains
    procedure :: allocate   => allocate_flux_type
    procedure :: deallocate => deallocate_flux_type
    procedure, nopass :: calc_surface_spectral
    procedure :: calc_toa_spectral
    procedure :: out_of_physical_bounds
    procedure :: heating_rate_out_of_physical_bounds
#if defined(_OPENACC)  || defined(OMPGPU)
    procedure, nopass :: create_device => create_device_flux
    procedure, nopass :: update_host   => update_host_flux
    procedure, nopass :: update_device => update_device_flux
    procedure, nopass :: delete_device => delete_device_flux
#endif
  end type flux_type

! Added for DWD (2020)
#ifdef DWD_VECTOR_OPTIMIZATIONS
  logical, parameter :: use_indexed_sum_vec = .true.
#else
  logical, parameter :: use_indexed_sum_vec = .false.
#endif

  !$omp declare target(indexed_sum)

contains

  !---------------------------------------------------------------------
  ! Allocate arrays for flux profiles, using config to define which
  ! fluxes are needed.  The arrays are dimensioned for columns between
  ! istartcol, iendcol and levels from 1 to nlev+1
  subroutine allocate_flux_type(this, config, istartcol, iendcol, nlev)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_io,     only : nulerr, radiation_abort
    use radiation_config, only : config_type

    integer, intent(in)             :: istartcol, iendcol, nlev
    class(flux_type), intent(inout) :: this
    type(config_type), intent(in)   :: config

    integer                         :: jcol
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_flux:allocate',0,hook_handle)

    ! Allocate longwave arrays
    if (config%do_lw) then
      allocate(this%lw_up(istartcol:iendcol,nlev+1))
      allocate(this%lw_dn(istartcol:iendcol,nlev+1))
      if (config%do_clear) then
        allocate(this%lw_up_clear(istartcol:iendcol,nlev+1))
        allocate(this%lw_dn_clear(istartcol:iendcol,nlev+1))
      end if

      if (config%do_save_spectral_flux) then
        if (config%n_spec_lw == 0) then
          write(nulerr,'(a)') '*** Error: number of LW spectral points to save not yet defined ' &
               & // 'so cannot allocate spectral flux arrays'
          call radiation_abort()
        end if

        allocate(this%lw_up_band(config%n_spec_lw,istartcol:iendcol,nlev+1))
        allocate(this%lw_dn_band(config%n_spec_lw,istartcol:iendcol,nlev+1))
        if (config%do_clear) then
          allocate(this%lw_up_clear_band(config%n_spec_lw, &
               &                         istartcol:iendcol,nlev+1))
          allocate(this%lw_dn_clear_band(config%n_spec_lw, &
               &                         istartcol:iendcol,nlev+1))
        end if
      end if

      if (config%do_lw_derivatives) then
        allocate(this%lw_derivatives(istartcol:iendcol,nlev+1))
      end if

      if (config%do_toa_spectral_flux) then
        if (config%n_bands_lw == 0) then
          write(nulerr,'(a)') '*** Error: number of LW bands not yet defined ' &
               & // 'so cannot allocate TOA spectral flux arrays'
          call radiation_abort()
        end if
        allocate(this%lw_up_toa_band(config%n_bands_lw, istartcol:iendcol))
        if (config%do_clear) then
          allocate(this%lw_up_toa_clear_band(config%n_bands_lw, istartcol:iendcol))
        end if
      end if

      ! Allocate g-point downwelling fluxes at surface, and TOA fluxes
      allocate(this%lw_dn_surf_g(config%n_g_lw,istartcol:iendcol))
      allocate(this%lw_up_toa_g (config%n_g_lw,istartcol:iendcol))
      if (config%do_clear) then
        allocate(this%lw_dn_surf_clear_g(config%n_g_lw,istartcol:iendcol))
        allocate(this%lw_up_toa_clear_g (config%n_g_lw,istartcol:iendcol))
      end if

      if (config%do_canopy_fluxes_lw) then
        ! Downward fluxes at top of canopy at the spectral resolution
        ! used in the canopy radiative transfer scheme
        allocate(this%lw_dn_surf_canopy(config%n_canopy_bands_lw,istartcol:iendcol))
      end if
    end if

    ! Allocate shortwave arrays
    if (config%do_sw) then
      allocate(this%sw_up(istartcol:iendcol,nlev+1))
      allocate(this%sw_dn(istartcol:iendcol,nlev+1))
      if (config%do_sw_direct) then
        allocate(this%sw_dn_direct(istartcol:iendcol,nlev+1))
      end if
      if (config%do_clear) then
        allocate(this%sw_up_clear(istartcol:iendcol,nlev+1))
        allocate(this%sw_dn_clear(istartcol:iendcol,nlev+1))
        if (config%do_sw_direct) then
          allocate(this%sw_dn_direct_clear(istartcol:iendcol,nlev+1))
        end if
      end if

      if (config%do_save_spectral_flux) then
        if (config%n_spec_sw == 0) then
          write(nulerr,'(a)') '*** Error: number of SW spectral points to save not yet defined ' &
               & // 'so cannot allocate spectral flux arrays'
          call radiation_abort()
        end if

        allocate(this%sw_up_band(config%n_spec_sw,istartcol:iendcol,nlev+1))
        allocate(this%sw_dn_band(config%n_spec_sw,istartcol:iendcol,nlev+1))

        if (config%do_sw_direct) then
          allocate(this%sw_dn_direct_band(config%n_spec_sw, &
               &                          istartcol:iendcol,nlev+1))
        end if
        if (config%do_clear) then
          allocate(this%sw_up_clear_band(config%n_spec_sw, &
               &                         istartcol:iendcol,nlev+1))
          allocate(this%sw_dn_clear_band(config%n_spec_sw, &
               &                         istartcol:iendcol,nlev+1))
          if (config%do_sw_direct) then
            allocate(this%sw_dn_direct_clear_band(config%n_spec_sw, &
                 &                                istartcol:iendcol, nlev+1))
          end if
        end if
      end if

      if (config%do_surface_sw_spectral_flux) then
        if (config%n_bands_sw == 0) then
          write(nulerr,'(a)') '*** Error: number of SW bands not yet defined ' &
               & // 'so cannot allocate TOA spectral flux arrays'
          call radiation_abort()
        end if
        allocate(this%sw_dn_surf_band(config%n_bands_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_band(config%n_bands_sw,istartcol:iendcol))
        if (config%do_clear) then
          allocate(this%sw_dn_surf_clear_band(config%n_bands_sw, &
               &                              istartcol:iendcol))
          allocate(this%sw_dn_direct_surf_clear_band(config%n_bands_sw, &
               &                                     istartcol:iendcol))
        end if
      end if

      if (config%do_toa_spectral_flux) then
        if (config%n_bands_sw == 0) then
          write(nulerr,'(a)') '*** Error: number of SW bands not yet defined ' &
               & // 'so cannot allocate surface spectral flux arrays'
          call radiation_abort()
        end if
        allocate(this%sw_dn_toa_band(config%n_bands_sw, istartcol:iendcol))
        allocate(this%sw_up_toa_band(config%n_bands_sw, istartcol:iendcol))
        if (config%do_clear) then
          allocate(this%sw_up_toa_clear_band(config%n_bands_sw, istartcol:iendcol))
        end if
      end if

      ! Allocate g-point downwelling fluxes at surface, and TOA fluxes
      allocate(this%sw_dn_diffuse_surf_g(config%n_g_sw,istartcol:iendcol))
      allocate(this%sw_dn_direct_surf_g (config%n_g_sw,istartcol:iendcol))
      allocate(this%sw_dn_toa_g         (config%n_g_sw,istartcol:iendcol))
      allocate(this%sw_up_toa_g         (config%n_g_sw,istartcol:iendcol))
      if (config%do_clear) then
        allocate(this%sw_dn_diffuse_surf_clear_g(config%n_g_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_clear_g (config%n_g_sw,istartcol:iendcol))
        allocate(this%sw_up_toa_clear_g         (config%n_g_sw,istartcol:iendcol))
      end if

      if (config%do_canopy_fluxes_sw) then
        ! Downward fluxes at top of canopy at the spectral resolution
        ! used in the canopy radiative transfer scheme
        allocate(this%sw_dn_diffuse_surf_canopy(config%n_canopy_bands_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_canopy (config%n_canopy_bands_sw,istartcol:iendcol))
      end if
    end if

    ! Allocate cloud cover arrays
    allocate(this%cloud_cover_lw(istartcol:iendcol))
    allocate(this%cloud_cover_sw(istartcol:iendcol))

    ! Some solvers may not write to cloud cover, so we initialize to
    ! an unphysical value
    do jcol = istartcol,iendcol
      this%cloud_cover_lw(jcol) = -1.0_jprb
      this%cloud_cover_sw(jcol) = -1.0_jprb
    end do

    if (lhook) call dr_hook('radiation_flux:allocate',1,hook_handle)

  end subroutine allocate_flux_type


  !---------------------------------------------------------------------
  ! Deallocate flux arrays
  subroutine deallocate_flux_type(this)

    use yomhook,          only : lhook, dr_hook, jphook

    class(flux_type), intent(inout) :: this
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_flux:deallocate',0,hook_handle)

    if (allocated(this%lw_up)) then
      deallocate(this%lw_up)
      if (allocated(this%lw_dn))       deallocate(this%lw_dn)
      if (allocated(this%lw_up_clear)) deallocate(this%lw_up_clear)
      if (allocated(this%lw_dn_clear)) deallocate(this%lw_dn_clear)
    end if

    if (allocated(this%sw_up)) then
      deallocate(this%sw_up)
      if (allocated(this%sw_dn))        deallocate(this%sw_dn)
      if (allocated(this%sw_up_clear))  deallocate(this%sw_up_clear)
      if (allocated(this%sw_dn_clear))  deallocate(this%sw_dn_clear)
      if (allocated(this%sw_dn_direct)) deallocate(this%sw_dn_direct)
      if (allocated(this%sw_dn_direct_clear)) &
           &   deallocate(this%sw_dn_direct_clear)
    end if

    if (allocated(this%lw_up_band)) then
      deallocate(this%lw_up_band)
      if (allocated(this%lw_dn_band))       deallocate(this%lw_dn_band)
      if (allocated(this%lw_up_clear_band)) deallocate(this%lw_up_clear_band)
      if (allocated(this%lw_dn_clear_band)) deallocate(this%lw_dn_clear_band)
    end if

    if (allocated(this%sw_up_band)) then
      deallocate(this%sw_up_band)
      if (allocated(this%sw_dn_band))        deallocate(this%sw_dn_band)
      if (allocated(this%sw_up_clear_band))  deallocate(this%sw_up_clear_band)
      if (allocated(this%sw_dn_clear_band))  deallocate(this%sw_dn_clear_band)
      if (allocated(this%sw_dn_direct_band)) deallocate(this%sw_dn_direct_band)
      if (allocated(this%sw_dn_direct_clear_band)) &
           &   deallocate(this%sw_dn_direct_clear_band)
    end if

    if (allocated(this%sw_dn_surf_band)) then
      deallocate(this%sw_dn_surf_band)
      deallocate(this%sw_dn_direct_surf_band)
    end if
    if (allocated(this%sw_dn_surf_clear_band)) then
      deallocate(this%sw_dn_surf_clear_band)
      deallocate(this%sw_dn_direct_surf_clear_band)
    end if

    if (allocated(this%lw_dn_surf_canopy)) deallocate(this%lw_dn_surf_canopy)
    if (allocated(this%sw_dn_diffuse_surf_canopy)) deallocate(this%sw_dn_diffuse_surf_canopy)
    if (allocated(this%sw_dn_direct_surf_canopy)) deallocate(this%sw_dn_direct_surf_canopy)

    if (allocated(this%cloud_cover_sw)) then
      deallocate(this%cloud_cover_sw)
    end if
    if (allocated(this%cloud_cover_lw)) then
      deallocate(this%cloud_cover_lw)
    end if

    if (allocated(this%lw_derivatives)) then
      deallocate(this%lw_derivatives)
    end if

    if (allocated(this%lw_dn_surf_g))               deallocate(this%lw_dn_surf_g)
    if (allocated(this%lw_dn_surf_clear_g))         deallocate(this%lw_dn_surf_clear_g)
    if (allocated(this%sw_dn_diffuse_surf_g))       deallocate(this%sw_dn_diffuse_surf_g)
    if (allocated(this%sw_dn_direct_surf_g))        deallocate(this%sw_dn_direct_surf_g)
    if (allocated(this%sw_dn_diffuse_surf_clear_g)) deallocate(this%sw_dn_diffuse_surf_clear_g)
    if (allocated(this%sw_dn_direct_surf_clear_g))  deallocate(this%sw_dn_direct_surf_clear_g)

    ! !$ACC EXIT DATA DELETE(this%lw_up_toa_g) ASYNC(1) IF(allocated(this%lw_up_toa_g))
    ! !$ACC EXIT DATA DELETE(this%sw_up_toa_g) ASYNC(1) IF(allocated(this%sw_up_toa_g))
    ! !$ACC EXIT DATA DELETE(this%sw_dn_toa_g) ASYNC(1) IF(allocated(this%sw_dn_toa_g))
    ! !$ACC EXIT DATA DELETE(this%lw_up_toa_clear_g) ASYNC(1) IF(allocated(this%lw_up_toa_clear_g))
    ! !$ACC EXIT DATA DELETE(this%sw_up_toa_clear_g) ASYNC(1) IF(allocated(this%sw_up_toa_clear_g))
    ! !$ACC WAIT
    if (allocated(this%lw_up_toa_g))                deallocate(this%lw_up_toa_g)
    if (allocated(this%sw_up_toa_g))                deallocate(this%sw_up_toa_g)
    if (allocated(this%sw_dn_toa_g))                deallocate(this%sw_dn_toa_g)
    if (allocated(this%lw_up_toa_clear_g))          deallocate(this%lw_up_toa_clear_g)
    if (allocated(this%sw_up_toa_clear_g))          deallocate(this%sw_up_toa_clear_g)

    if (lhook) call dr_hook('radiation_flux:deallocate',1,hook_handle)

  end subroutine deallocate_flux_type


  !---------------------------------------------------------------------
  ! Calculate surface downwelling fluxes in each band using the
  ! downwelling surface fluxes at each g point
  subroutine calc_surface_spectral(this, config, istartcol, iendcol)

    use yomhook,          only : lhook, dr_hook, jphook
#if defined(_OPENACC) || defined(OMPGPU)
    use radiation_io,     only : nulerr, radiation_abort
#endif
    use radiation_config, only : config_type

#ifdef HAVE_NVTX
    use nvtx
#endif
#ifdef HAVE_ROCTX
    use roctx_profiling, only: roctxmarka, roctxrangepusha, roctxrangepop
    use iso_c_binding, only: c_null_char
    integer(4) :: roctx_ret
#endif

    type(flux_type),  intent(inout) :: this
    type(config_type), intent(in)    :: config
    integer,           intent(in)    :: istartcol, iendcol

    integer :: jcol, jband, jalbedoband, nalbedoband, jg

    ! Longwave surface downwelling in each band needed to compute
    ! canopy fluxes
    real(jprb) :: lw_dn_surf_band(config%n_bands_lw,istartcol:iendcol)

    real(jphook) :: hook_handle
#if defined(OMPGPU)
    integer :: istart, iend, ig
    real(jprb) :: s1, s2
#endif
    
#ifdef HAVE_NVTX
      call nvtxStartRange("radiation_flux::calc_surface_spectral")
#endif
#ifdef HAVE_ROCTX
      roctx_ret = roctxRangePushA("radiation_flux::calc_surface_spectral"//c_null_char)
#endif

    if (lhook) call dr_hook('radiation_flux:calc_surface_spectral',0,hook_handle)

#if defined(_OPENACC) || defined(OMPGPU)
    if (use_indexed_sum_vec) then
      write(nulerr,'(a)') '*** Error: radiation_flux:calc_surface_spectral use_indexed_sum_vec==.true. not ported to GPU'
      call radiation_abort()
    endif
#endif

    if (config%do_sw .and. config%do_surface_sw_spectral_flux) then

      if (use_indexed_sum_vec) then
        call indexed_sum_vec(this%sw_dn_direct_surf_g, &
             &               config%i_band_from_reordered_g_sw, &
             &               this%sw_dn_direct_surf_band, istartcol, iendcol)
        call indexed_sum_vec(this%sw_dn_diffuse_surf_g, &
             &               config%i_band_from_reordered_g_sw, &
             &               this%sw_dn_surf_band, istartcol, iendcol)
        do jcol = istartcol,iendcol
          this%sw_dn_surf_band(:,jcol) &
               &  = this%sw_dn_surf_band(:,jcol) &
               &  + this%sw_dn_direct_surf_band(:,jcol)
        end do
     else

#if defined(OMPGPU)
        istart = lbound(this%sw_dn_surf_band,1)
        iend = ubound(this%sw_dn_surf_band,1)
#endif
        !!$OMP TARGET DATA MAP(PRESENT, ALLOC: this%sw_dn_direct_surf_g,  config%i_band_from_reordered_g_sw, &
        !!$OMP   this%sw_dn_direct_surf_band, this%sw_dn_diffuse_surf_g, this%sw_dn_surf_band)

#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
        !$OMP TARGET TEAMS DISTRIBUTE
#endif
        !$ACC PARALLEL DEFAULT(PRESENT) NUM_GANGS(iendcol-istartcol+1) NUM_WORKERS(1) &
        !$ACC   VECTOR_LENGTH(32*((config%n_g_sw-1)/32+1)) ASYNC(1)
        !$ACC LOOP GANG
        do jcol = istartcol,iendcol
          call indexed_sum(this%sw_dn_direct_surf_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_dn_direct_surf_band(:,jcol))
          call indexed_sum(this%sw_dn_diffuse_surf_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_dn_surf_band(:,jcol))
#if defined(OMPGPU)
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP PARALLEL DO
#endif
          do ig = istart, iend
             this%sw_dn_surf_band(ig,jcol) &
                  &  = this%sw_dn_surf_band(ig,jcol) &
                  &  + this%sw_dn_direct_surf_band(ig,jcol)
          end do
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP END PARALLEL DO
#endif
#else
          this%sw_dn_surf_band(:,jcol) &
               &  = this%sw_dn_surf_band(:,jcol) &
               &  + this%sw_dn_direct_surf_band(:,jcol)
#endif
        end do
        !$ACC END PARALLEL
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
        !$OMP END TARGET TEAMS DISTRIBUTE
#endif
        !!$OMP END TARGET DATA
      end if

      if (config%do_clear) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%sw_dn_direct_surf_clear_g, &
               &               config%i_band_from_reordered_g_sw, &
               &               this%sw_dn_direct_surf_clear_band, istartcol, iendcol)
          call indexed_sum_vec(this%sw_dn_diffuse_surf_clear_g, &
               &               config%i_band_from_reordered_g_sw, &
               &               this%sw_dn_surf_clear_band, istartcol, iendcol)
          do jcol = istartcol,iendcol
            this%sw_dn_surf_clear_band(:,jcol) &
                 &  = this%sw_dn_surf_clear_band(:,jcol) &
                 &  + this%sw_dn_direct_surf_clear_band(:,jcol)
          end do
       else
#if defined(OMPGPU)
          istart = lbound(this%sw_dn_surf_clear_band,1)
          iend = ubound(this%sw_dn_surf_clear_band,1)
#endif
          !!$OMP TARGET DATA MAP(PRESENT, ALLOC: &
          !!$OMP   this%sw_dn_direct_surf_clear_g, config%i_band_from_reordered_g_sw, this%sw_dn_direct_surf_clear_band,  &
          !!$OMP   this%sw_dn_diffuse_surf_clear_g, this%sw_dn_surf_clear_band)
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
          !$OMP TARGET TEAMS DISTRIBUTE          
#endif
          !$ACC PARALLEL DEFAULT(PRESENT) NUM_GANGS(iendcol-istartcol+1) NUM_WORKERS(1) &
          !$ACC   VECTOR_LENGTH(32*(config%n_g_sw-1)/32+1) ASYNC(1)
          !$ACC LOOP GANG
          do jcol = istartcol,iendcol
            call indexed_sum(this%sw_dn_direct_surf_clear_g(:,jcol), &
                 &           config%i_band_from_reordered_g_sw, &
                 &           this%sw_dn_direct_surf_clear_band(:,jcol))
            call indexed_sum(this%sw_dn_diffuse_surf_clear_g(:,jcol), &
                 &           config%i_band_from_reordered_g_sw, &
                 &           this%sw_dn_surf_clear_band(:,jcol))
#if defined(OMPGPU)
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
            !$OMP PARALLEL DO
#endif
            do ig = istart, iend
               this%sw_dn_surf_clear_band(ig,jcol) &
                    &  = this%sw_dn_surf_clear_band(ig,jcol) &
                    &  + this%sw_dn_direct_surf_clear_band(ig,jcol)
            end do
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
            !$OMP END PARALLEL DO
#endif
#else
            this%sw_dn_surf_clear_band(:,jcol) &
                 &  = this%sw_dn_surf_clear_band(:,jcol) &
                 &  + this%sw_dn_direct_surf_clear_band(:,jcol)
#endif
          end do
          !$ACC END PARALLEL
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
          !$OMP END TARGET TEAMS DISTRIBUTE
#endif
          !!$OMP END TARGET DATA
       end if
      end if

    end if ! do_surface_sw_spectral_flux

    ! Fluxes in bands required for canopy radiative transfer
    if (config%do_sw .and. config%do_canopy_fluxes_sw) then
      if (config%use_canopy_full_spectrum_sw) then
        !$OMP TARGET DATA MAP(PRESENT, ALLOC: this%sw_dn_diffuse_surf_canopy, this%sw_dn_diffuse_surf_g, &
        !$OMP   this%sw_dn_direct_surf_canopy, this%sw_dn_direct_surf_g)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        do jcol = istartcol,iendcol
          do jg = 1,config%n_g_sw
            this%sw_dn_diffuse_surf_canopy(jg,jcol) = this%sw_dn_diffuse_surf_g(jg,jcol)
            this%sw_dn_direct_surf_canopy (jg,jcol) = this%sw_dn_direct_surf_g (jg,jcol)
          end do
        end do
        !$ACC END PARALLEL
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        !$OMP END TARGET DATA
      else if (config%do_nearest_spectral_sw_albedo) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%sw_dn_direct_surf_g, &
               &               config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
               &               this%sw_dn_direct_surf_canopy, istartcol, iendcol)
          call indexed_sum_vec(this%sw_dn_diffuse_surf_g, &
               &               config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
               &               this%sw_dn_diffuse_surf_canopy, istartcol, iendcol)
        else
          !$OMP TARGET DATA MAP(PRESENT, ALLOC: &
          !$OMP   this%sw_dn_direct_surf_g, config%i_albedo_from_band_sw, config%i_band_from_reordered_g_sw, &
          !$OMP   this%sw_dn_direct_surf_canopy, this%sw_dn_diffuse_surf_g, this%sw_dn_diffuse_surf_canopy)
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
          !$OMP TARGET TEAMS DISTRIBUTE
#endif
          !$ACC PARALLEL DEFAULT(PRESENT) NUM_GANGS(iendcol-istartcol+1) NUM_WORKERS(1) &
          !$ACC   VECTOR_LENGTH(32*(config%n_g_sw-1)/32+1) ASYNC(1)
          !$ACC LOOP GANG
          do jcol = istartcol,iendcol
            call indexed_sum(this%sw_dn_direct_surf_g(:,jcol), &
                 &           config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
                 &           this%sw_dn_direct_surf_canopy(:,jcol))
            call indexed_sum(this%sw_dn_diffuse_surf_g(:,jcol), &
                 &           config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
                 &           this%sw_dn_diffuse_surf_canopy(:,jcol))
          end do
          !$ACC END PARALLEL
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
          !$OMP END TARGET TEAMS DISTRIBUTE
#endif
          !$OMP END TARGET DATA
        end if
      else
        ! More accurate calculations using weights, but requires
        ! this%sw_dn_[direct_]surf_band to be defined, i.e.
        ! config%do_surface_sw_spectral_flux == .true.
        nalbedoband = size(config%sw_albedo_weights,1)
        !!$OMP TARGET DATA MAP(PRESENT, ALLOC: this%sw_dn_diffuse_surf_canopy,  this%sw_dn_direct_surf_canopy)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) &
        !$ACC     PRESENT(this%sw_dn_diffuse_surf_canopy, this%sw_dn_direct_surf_canopy, &
        !$ACC             config%sw_albedo_weights)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        do jcol = istartcol,iendcol
          do jalbedoband = 1,nalbedoband
            this%sw_dn_diffuse_surf_canopy(jalbedoband,jcol) = 0.0_jprb
            this%sw_dn_direct_surf_canopy (jalbedoband,jcol) = 0.0_jprb
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        !!$OMP END TARGET DATA

        !!$OMP TARGET DATA MAP(PRESENT, ALLOC: &
        !!$OMP   config%sw_albedo_weights, this%sw_dn_diffuse_surf_canopy, this%sw_dn_diffuse_surf_canopy, this%sw_dn_surf_band, &
        !!$OMP   this%sw_dn_direct_surf_canopy, this%sw_dn_direct_surf_canopy, this%sw_dn_direct_surf_band)
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
#else
        !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2)
#endif
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        do jcol = istartcol, iendcol
          do jalbedoband = 1,nalbedoband
#if defined(OMPGPU)
            s1 = 0 
            s2 = 0
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
            !$OMP PARALLEL DO REDUCTION(+ : s1, s2)
#endif
            do jband = 1,config%n_bands_sw
              if (config%sw_albedo_weights(jalbedoband,jband) /= 0.0_jprb) then
                ! Initially, "diffuse" is actually "total"
                s1 = s1 + config%sw_albedo_weights(jalbedoband,jband) &
                    &    * this%sw_dn_surf_band(jband,jcol)
                s2 = s2 + config%sw_albedo_weights(jalbedoband,jband) &
                   &    * this%sw_dn_direct_surf_band(jband,jcol)
              end if
            end do
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
            !$OMP END PARALLEL DO
#endif
            this%sw_dn_diffuse_surf_canopy(jalbedoband,jcol) = this%sw_dn_diffuse_surf_canopy(jalbedoband,jcol) + s1
            this%sw_dn_direct_surf_canopy(jalbedoband,jcol) = this%sw_dn_direct_surf_canopy(jalbedoband,jcol) + s2
#else
            !$ACC LOOP SEQ
            do jband = 1,config%n_bands_sw
              if (config%sw_albedo_weights(jalbedoband,jband) /= 0.0_jprb) then
                ! Initially, "diffuse" is actually "total"
                this%sw_dn_diffuse_surf_canopy(jalbedoband,jcol) &
                    &  = this%sw_dn_diffuse_surf_canopy(jalbedoband,jcol) &
                    &  + config%sw_albedo_weights(jalbedoband,jband) &
                    &    * this%sw_dn_surf_band(jband,jcol)
                this%sw_dn_direct_surf_canopy(jalbedoband,jcol) &
                    &  = this%sw_dn_direct_surf_canopy(jalbedoband,jcol) &
                    &  + config%sw_albedo_weights(jalbedoband,jband) &
                    &    * this%sw_dn_direct_surf_band(jband,jcol)
              end if
            end do
#endif 
          end do
        end do
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
        !$OMP END TARGET TEAMS DISTRIBUTE
#endif
        !!$OMP END TARGET DATA

        !!$OMP TARGET DATA MAP(PRESENT, ALLOC: this%sw_dn_diffuse_surf_canopy, this%sw_dn_direct_surf_canopy)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        do jcol = istartcol,iendcol
          do jalbedoband = 1,nalbedoband
            ! Subtract the direct from total to get diffuse
            this%sw_dn_diffuse_surf_canopy(jalbedoband,jcol) &
                &  = this%sw_dn_diffuse_surf_canopy(jalbedoband,jcol) &
                &  - this%sw_dn_direct_surf_canopy(jalbedoband,jcol)
          end do
        end do
        !$ACC END PARALLEL
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        !!$OMP END TARGET DATA
      end if

    end if ! do_canopy_fluxes_sw

    if (config%do_lw .and. config%do_canopy_fluxes_lw) then
      if (config%use_canopy_full_spectrum_lw) then
        !!$OMP TARGET DATA MAP(PRESENT, ALLOC: this%lw_dn_surf_canopy, this%lw_dn_surf_g)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        do jcol = istartcol,iendcol
          do jg = 1,config%n_g_lw
            this%lw_dn_surf_canopy(jg,jcol) = this%lw_dn_surf_g(jg,jcol)
          end do
        end do
        !$ACC END PARALLEL
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        !!$OMP END TARGET DATA
      else if (config%do_nearest_spectral_lw_emiss) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%lw_dn_surf_g, &
               &               config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw), &
               &               this%lw_dn_surf_canopy, istartcol, iendcol)
        else
#ifdef _OPENACC
          !!$OMP TARGET DATA MAP(PRESENT, ALLOC: this%lw_dn_surf_canopy, config%i_emiss_from_band_lw, config%i_band_from_reordered_g_lw, &
          !!$OMP   this%lw_dn_surf_g)
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
          !$OMP TARGET TEAMS DISTRIBUTE
#endif
          !$ACC PARALLEL DEFAULT(PRESENT) &
          !$ACC   NUM_GANGS(iendcol-istartcol+1) NUM_WORKERS(1) &
          !$ACC   VECTOR_LENGTH(32*(config%n_g_lw-1)/32+1) ASYNC(1)
          !$ACC LOOP GANG
          do jcol = istartcol,iendcol
            ! Inlined indexed_sum because of an on-device segfault due to the nested
            ! array index-subscript passed as `ind` to indexed_sum
            this%lw_dn_surf_canopy(:,jcol) = 0.0
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
            !$OMP PARALLEL DO PRIVATE(jband)
#endif
            !$ACC LOOP VECTOR
            do jg = 1,config%n_g_lw
              jband = config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw(jg))
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
              !$OMP ATOMIC UPDATE
#endif
              !$ACC ATOMIC UPDATE
              this%lw_dn_surf_canopy(jband,jcol) = this%lw_dn_surf_canopy(jband,jcol) + this%lw_dn_surf_g(jg,jcol)
              !$ACC END ATOMIC
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
              !$OMP END ATOMIC
#endif
            end do
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
            !$OMP END PARALLEL DO
#endif
          end do
          !$ACC END PARALLEL
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
          !$OMP END TARGET TEAMS DISTRIBUTE
#endif
          !!$OMP END TARGET DATA
#else
          do jcol = istartcol,iendcol
            call indexed_sum(this%lw_dn_surf_g(:,jcol), &
                 &           config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw), &
                 &           this%lw_dn_surf_canopy(:,jcol))
          end do
#endif
        end if
      else
        !$OMP TARGET ENTER DATA MAP(ALLOC: lw_dn_surf_band)
        !$ACC DATA CREATE(lw_dn_surf_band) ASYNC(1)
        ! Compute fluxes in each longwave emissivity interval using
        ! weights; first sum over g points to get the values in bands
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%lw_dn_surf_g, &
               &               config%i_band_from_reordered_g_lw, &
               &               lw_dn_surf_band, istartcol, iendcol)
        else
          !!$OMP TARGET DATA MAP(PRESENT, ALLOC: this%lw_dn_surf_g, config%i_band_from_reordered_g_lw)
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
          !$OMP TARGET TEAMS DISTRIBUTE
#endif
          !$ACC PARALLEL DEFAULT(PRESENT) NUM_GANGS(iendcol-istartcol+1) NUM_WORKERS(1) &
          !$ACC   VECTOR_LENGTH(32*(config%n_g_lw-1)/32+1) ASYNC(1)
          !$ACC LOOP GANG
          do jcol = istartcol,iendcol
            call indexed_sum(this%lw_dn_surf_g(:,jcol), &
                 &           config%i_band_from_reordered_g_lw, &
                 &           lw_dn_surf_band(:,jcol))
          end do
          !$ACC END PARALLEL
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
          !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
          !$OMP END TARGET TEAMS DISTRIBUTE
#endif
          !!$OMP END TARGET DATA
        end if
        nalbedoband = size(config%lw_emiss_weights,1)
        !!$OMP TARGET DATA MAP(PRESENT, ALLOC: this%lw_dn_surf_canopy)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
        !$ACC PARALLEL DEFAULT(PRESENT) &
        !$ACC   PRESENT(this%lw_dn_surf_canopy, config%lw_emiss_weights) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        do jcol = istartcol,iendcol
          do jalbedoband = 1,nalbedoband
            this%lw_dn_surf_canopy(jalbedoband,jcol) = 0.0_jprb
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        !!$OMP END TARGET DATA

        !!$OMP TARGET DATA MAP(PRESENT, ALLOC: config%lw_emiss_weights, this%lw_dn_surf_canopy, this%lw_dn_surf_canopy, config%lw_emiss_weights)
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
#else
        !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2)
#endif
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        do jcol = istartcol,iendcol
          do jalbedoband = 1,nalbedoband
#if defined(OMPGPU)
            s1 = 0
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
            !$OMP PARALLEL DO REDUCTION(+:s1)
#endif
            do jband = 1,config%n_bands_lw
               if (config%lw_emiss_weights(jalbedoband,jband) /= 0.0_jprb) then
                  s1 = s1 + config%lw_emiss_weights(jalbedoband,jband) &
                   &    * lw_dn_surf_band(jband,jcol)
               end if
            end do
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
            !$OMP END PARALLEL DO
#endif
            this%lw_dn_surf_canopy(jalbedoband,jcol) &
                 &  = this%lw_dn_surf_canopy(jalbedoband,jcol) + s1
#else
            !$ACC LOOP SEQ
            do jband = 1,config%n_bands_lw
               if (config%lw_emiss_weights(jalbedoband,jband) /= 0.0_jprb) then
                  this%lw_dn_surf_canopy(jalbedoband,jcol) &
                       &  = this%lw_dn_surf_canopy(jalbedoband,jcol) &
                       &  + config%lw_emiss_weights(jalbedoband,jband) &
                       &    * lw_dn_surf_band(jband,jcol)
               end if
            end do
#endif
        end do
        end do
        !$ACC END PARALLEL
#ifdef WORKAROUND_NESTED_PARALLEL_FLUX
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
        !$OMP END TARGET TEAMS DISTRIBUTE
#endif
        !!$OMP END TARGET DATA

        !$ACC END DATA
        !$OMP TARGET EXIT DATA MAP(DELETE: lw_dn_surf_band)
      end if
    end if

    if (lhook) call dr_hook('radiation_flux:calc_surface_spectral',1,hook_handle)

#ifdef HAVE_NVTX
    !$ACC WAIT(1)
    call nvtxEndRange
#endif
#ifdef HAVE_ROCTX
    call roctxRangePop()
    call roctxMarkA("radiation_flux:calc_surface_spectra"//c_null_char)
#endif
  end subroutine calc_surface_spectral


  !---------------------------------------------------------------------
  ! Calculate top-of-atmosphere fluxes in each band using the fluxes
  ! at each g point
  subroutine calc_toa_spectral(this, config, istartcol, iendcol)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_config, only : config_type
#if defined(_OPENACC) || defined(OMPGPU)
    use radiation_io,     only : nulerr, radiation_abort
#endif


    class(flux_type),  intent(inout) :: this
    type(config_type), intent(in)    :: config
    integer,           intent(in)    :: istartcol, iendcol

    integer :: jcol, jband

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_flux:calc_toa_spectral',0,hook_handle)

#if defined(_OPENACC) || defined(OMPGPU)
    if (config%do_toa_spectral_flux) then
      write(nulerr,'(a)') '*** Error: radiation_flux:calc_toa_spectral not ported to GPU.'
      call radiation_abort()
    end if
#endif

    if (config%do_sw .and. config%do_toa_spectral_flux) then
      if (use_indexed_sum_vec) then
        call indexed_sum_vec(this%sw_dn_toa_g, &
             &               config%i_band_from_reordered_g_sw, &
             &               this%sw_dn_toa_band, istartcol, iendcol)
        call indexed_sum_vec(this%sw_up_toa_g, &
             &               config%i_band_from_reordered_g_sw, &
             &               this%sw_up_toa_band, istartcol, iendcol)
      else
        do jcol = istartcol,iendcol
          call indexed_sum(this%sw_dn_toa_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_dn_toa_band(:,jcol))
          call indexed_sum(this%sw_up_toa_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_up_toa_band(:,jcol))
        end do
      end if

      if (config%do_clear) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%sw_up_toa_clear_g, &
               &               config%i_band_from_reordered_g_sw, &
               &               this%sw_up_toa_clear_band, istartcol, iendcol)
        else
          do jcol = istartcol,iendcol
            call indexed_sum(this%sw_up_toa_clear_g(:,jcol), &
                 &               config%i_band_from_reordered_g_sw, &
                 &               this%sw_up_toa_clear_band(:,jcol))
          end do
        end if
      end if
    end if

    if (config%do_lw .and. config%do_toa_spectral_flux) then
      if (use_indexed_sum_vec) then
        call indexed_sum_vec(this%lw_up_toa_g, &
             &               config%i_band_from_reordered_g_lw, &
             &               this%lw_up_toa_band, istartcol, iendcol)
      else
        do jcol = istartcol,iendcol
          call indexed_sum(this%lw_up_toa_g(:,jcol), &
               &           config%i_band_from_reordered_g_lw, &
               &           this%lw_up_toa_band(:,jcol))
        end do
      end if

      if (config%do_clear) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%lw_up_toa_clear_g, &
               &               config%i_band_from_reordered_g_lw, &
               &               this%lw_up_toa_clear_band, istartcol, iendcol)
        else
          do jcol = istartcol,iendcol
            call indexed_sum(this%lw_up_toa_clear_g(:,jcol), &
                 &               config%i_band_from_reordered_g_lw, &
                 &               this%lw_up_toa_clear_band(:,jcol))
          end do
        end if
      end if
    end if

    if (lhook) call dr_hook('radiation_flux:calc_toa_spectral',1,hook_handle)

  end subroutine calc_toa_spectral


  !---------------------------------------------------------------------
  ! Return .true. if the most important flux variables are out of a
  ! physically sensible range, optionally only considering columns
  ! between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol) result(is_bad)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_2d

    class(flux_type), intent(inout) :: this
    integer, optional,intent(in) :: istartcol, iendcol
    logical                      :: is_bad

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_flux:out_of_physical_bounds',0,hook_handle)

    is_bad =    out_of_bounds_2d(this%lw_up, 'lw_up', 10.0_jprb, 900.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%lw_dn, 'lw_dn', 0.0_jprb,  800.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_up, 'sw_up', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn, 'sw_dn', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_direct, 'sw_dn_direct', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%lw_derivatives, 'lw_derivatives', 0.0_jprb, 1.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_surf_band, 'sw_dn_surf_band', 0.0_jprb, 1500.0_jprb, &
         &                       .false., j1=istartcol, j2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_surf_clear_band, 'sw_dn_surf_clear_band', 0.0_jprb, 1500.0_jprb, &
         &                       .false., j1=istartcol, j2=iendcol)

    if (lhook) call dr_hook('radiation_flux:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds

  !---------------------------------------------------------------------
  ! Return .true. if the heating rates are out of a physically
  ! sensible range, optionally only considering columns between
  ! istartcol and iendcol. This function allocates and deallocates
  ! memory due to the requirements for inputs of out_of_bounds_2d.
  function heating_rate_out_of_physical_bounds(this, nlev, istartcol, iendcol, pressure_hl) result(is_bad)

    use radiation_check, only : out_of_bounds_2d
    use radiation_constants, only : AccelDueToGravity

    ! "Cp" (J kg-1 K-1)
    real(jprb), parameter :: SpecificHeatDryAir = 1004.0

    class(flux_type), intent(inout) :: this
    integer, intent(in) :: istartcol, iendcol, nlev
    logical                      :: is_bad

    real(jprb), intent(in) :: pressure_hl(:,:)

    real(jprb), allocatable :: hr_K_day(:,:)

    real(jprb) :: scaling(istartcol:iendcol,nlev)

    allocate(hr_K_day(istartcol:iendcol,nlev))

    scaling = -(24.0_jprb * 3600.0_jprb * AccelDueToGravity / SpecificHeatDryAir) &
         &  / (pressure_hl(istartcol:iendcol,2:nlev+1) - pressure_hl(istartcol:iendcol,1:nlev))
    ! Shortwave
    hr_K_day = scaling * (this%sw_dn(istartcol:iendcol,2:nlev+1) - this%sw_up(istartcol:iendcol,2:nlev+1) &
         &               -this%sw_dn(istartcol:iendcol,1:nlev)   + this%sw_up(istartcol:iendcol,1:nlev))
    is_bad = out_of_bounds_2d(hr_K_day, 'sw_heating_rate_K_day', 0.0_jprb, 200.0_jprb, &
         &                    .false., i1=istartcol, i2=iendcol)

    ! Longwave
    hr_K_day = scaling * (this%lw_dn(istartcol:iendcol,2:nlev+1) - this%lw_up(istartcol:iendcol,2:nlev+1) &
         &               -this%lw_dn(istartcol:iendcol,1:nlev)   + this%lw_up(istartcol:iendcol,1:nlev))
    is_bad = is_bad .or. out_of_bounds_2d(hr_K_day, 'lw_heating_rate_K_day', -250.0_jprb, 150.0_jprb, &
         &                                .false., i1=istartcol, i2=iendcol)

    deallocate(hr_K_day)

  end function heating_rate_out_of_physical_bounds


  !---------------------------------------------------------------------
  ! Sum elements of "source" into "dest" according to index "ind".
  ! "source" and "ind" should have the same size and bounds, and no
  ! element of "ind" should refer outside the bounds of "dest".  This
  ! version increments existing contents of "dest".
  pure subroutine add_indexed_sum(source, ind, dest)

    real(jprb), intent(in)    :: source(:)
    integer,    intent(in)    :: ind(:)
    real(jprb), intent(inout) :: dest(:)

    integer :: ig, jg, istart, iend

    istart = lbound(source,1)
    iend   = ubound(source,1)

    do jg = istart, iend
      ig = ind(jg)
      dest(ig) = dest(ig) + source(jg)
    end do

  end subroutine add_indexed_sum


  !---------------------------------------------------------------------
  ! As "add_indexed_sum" but this version overwrites existing contents
  ! of "dest"
  pure subroutine indexed_sum(source, ind, dest)

    real(jprb), intent(in)  :: source(:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:)

    integer :: ig, jg, istart, iend

    !$ACC ROUTINE VECTOR

    istart = lbound(dest,1)
    iend   = ubound(dest,1)

#if defined(OMPGPU)
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
    !$OMP PARALLEL DO
#endif
    do jg = istart, iend
       dest(jg) = 0.0
    end do
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
    !$OMP END PARALLEL DO
#endif
#else
    dest = 0.0
#endif

    istart = lbound(source,1)
    iend   = ubound(source,1)
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
    !$OMP PARALLEL DO PRIVATE(ig)
#endif
    !$ACC LOOP VECTOR
    do jg = istart, iend
      ig = ind(jg)
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
      !$OMP ATOMIC UPDATE 
#endif
      !$ACC ATOMIC UPDATE
      dest(ig) = dest(ig) + source(jg)
      !$ACC END ATOMIC
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
      !$END OMP ATOMIC
#endif
   end do
#ifndef WORKAROUND_NESTED_PARALLEL_FLUX
    !$OMP END PARALLEL DO
#endif
  end subroutine indexed_sum

  !---------------------------------------------------------------------
  ! Vectorized version of "add_indexed_sum"
  subroutine indexed_sum_vec(source, ind, dest, ist, iend)

    real(jprb), intent(in)  :: source(:,:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:,:)
    integer,    intent(in)  :: ist, iend

    integer :: ig, jg, jc

    dest = 0.0

    do jg = lbound(source,1), ubound(source,1)
      ig = ind(jg)
      do jc = ist, iend
        dest(ig,jc) = dest(ig,jc) + source(jg,jc)
      end do
    end do

  end subroutine indexed_sum_vec

  !---------------------------------------------------------------------
  ! As "add_indexed_sum" but a whole vertical profiles
  pure subroutine add_indexed_sum_profile(source, ind, dest)

    real(jprb), intent(in)  :: source(:,:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:,:)

    integer :: ig, jg, istart, iend, jlev, nlev

    istart = lbound(source,1)
    iend   = ubound(source,1)
    nlev   = size(source,2)

    do jlev = 1,nlev
      do jg = istart, iend
        ig = ind(jg)
        dest(ig,jlev) = dest(ig,jlev) + source(jg,jlev)
      end do
    end do

  end subroutine add_indexed_sum_profile


  !---------------------------------------------------------------------
  ! As "indexed_sum" but a whole vertical profiles
  pure subroutine indexed_sum_profile(source, ind, dest)

    real(jprb), intent(in)  :: source(:,:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:,:)

    integer :: ig, jg, istart, iend, jlev, nlev

    dest = 0.0

    istart = lbound(source,1)
    iend   = ubound(source,1)
    nlev   = size(source,2)

    do jlev = 1,nlev
      do jg = istart, iend
        ig = ind(jg)
        dest(ig,jlev) = dest(ig,jlev) + source(jg,jlev)
      end do
    end do

  end subroutine indexed_sum_profile

#if defined(_OPENACC)  || defined(OMPGPU)
  !---------------------------------------------------------------------
  ! debug dump
  subroutine debug_dump(this,istartcol,iendcol,k)
    use radiation_io,     only : nulout
    type(flux_type), intent(in) :: this
    integer, intent(in) :: istartcol, iendcol,k
    integer :: file_idx, fidx1, fidx2,shift
    character*100 :: FNAME
    !$OMP TARGET UPDATE FROM(this%sw_up) IF(allocated(this%sw_up))
    !$OMP TARGET UPDATE FROM(this%sw_dn) IF(allocated(this%sw_dn))
    !$OMP TARGET UPDATE FROM(this%sw_up_clear) IF(allocated(this%sw_up_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_clear) IF(allocated(this%sw_dn_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct) IF(allocated(this%sw_dn_direct))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy))
    !$OMP TARGET UPDATE FROM(this%lw_up) IF(allocated(this%lw_up))
    !$OMP TARGET UPDATE FROM(this%lw_dn) IF(allocated(this%lw_dn))
    !$OMP TARGET UPDATE FROM(this%lw_up_clear) IF(allocated(this%lw_up_clear))
    !$OMP TARGET UPDATE FROM(this%lw_dn_clear) IF(allocated(this%lw_dn_clear))
    !$OMP TARGET UPDATE FROM(this%lw_derivatives) IF(allocated(this%lw_derivatives))
    !$OMP TARGET UPDATE FROM(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band))

    !$ACC WAIT(1)
    !$ACC UPDATE HOST(this%lw_up) IF(allocated(this%lw_up)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn) IF(allocated(this%lw_dn)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_clear) IF(allocated(this%lw_up_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_clear) IF(allocated(this%lw_dn_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_derivatives) IF(allocated(this%lw_derivatives)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up) IF(allocated(this%sw_up)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn) IF(allocated(this%sw_dn)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_clear) IF(allocated(this%sw_up_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_clear) IF(allocated(this%sw_dn_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct) IF(allocated(this%sw_dn_direct)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band)) ASYNC(1)
    !$ACC WAIT(1)
    shift = k-1
    write(*,*) "shift=",shift

    if (allocated(this%sw_dn_surf_band)) then
       FNAME='sw_dn_surf_band.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%sw_dn_surf_band, dim=1), ubound(this%sw_dn_surf_band, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2+shift,",",this%sw_dn_surf_band(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_dn_surf_clear_band)) then
       FNAME='sw_dn_surf_clear_band.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%sw_dn_surf_clear_band, dim=1), ubound(this%sw_dn_surf_clear_band, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2+shift,",",this%sw_dn_surf_clear_band(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_dn_diffuse_surf_canopy)) then
       FNAME='sw_dn_diffuse_surf_canopy.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%sw_dn_diffuse_surf_canopy, dim=1), ubound(this%sw_dn_diffuse_surf_canopy, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2+shift,",",this%sw_dn_diffuse_surf_canopy(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_dn_direct_surf_canopy)) then
       FNAME='sw_dn_direct_surf_canopy.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%sw_dn_direct_surf_canopy, dim=1), ubound(this%sw_dn_direct_surf_canopy, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2+shift,",",this%sw_dn_direct_surf_canopy(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_up)) then
       FNAME='sw_up.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%sw_up, dim=2), ubound(this%sw_up, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%sw_up(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_dn)) then
       FNAME='sw_dn.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%sw_dn, dim=2), ubound(this%sw_dn, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%sw_dn(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_up_clear)) then
       FNAME='sw_up_clear.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%sw_up_clear, dim=2), ubound(this%sw_up_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%sw_up_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%sw_dn_clear)) then
       FNAME='sw_dn_clear.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%sw_dn_clear, dim=2), ubound(this%sw_dn_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%sw_dn_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%sw_dn_direct_clear)) then
       FNAME='sw_dn_direct_clear.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%sw_dn_direct_clear, dim=2), ubound(this%sw_dn_direct_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%sw_dn_direct_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_dn_direct)) then
       FNAME='sw_dn_direct.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%sw_dn_direct, dim=2), ubound(this%sw_dn_direct, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%sw_dn_direct(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%lw_up)) then
       FNAME='lw_up.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%lw_up, dim=2), ubound(this%lw_up, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%lw_up(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%lw_dn)) then
       FNAME='lw_dn.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%lw_dn, dim=2), ubound(this%lw_dn, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%lw_dn(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%lw_up_clear)) then
       FNAME='lw_up_clear.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%lw_up_clear, dim=2), ubound(this%lw_up_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%lw_up_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%lw_dn_clear)) then
       FNAME='lw_dn_clear.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%lw_dn_clear, dim=2), ubound(this%lw_dn_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%lw_dn_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%lw_derivatives)) then
       FNAME='lw_derivatives.txt'
       if (k==1) then
          open(newunit=file_idx, file=FNAME, action='write',status='replace')
       else
          open(newunit=file_idx, file=FNAME, action='write',position='append')
       endif

       do fidx2 = lbound(this%lw_derivatives, dim=2), ubound(this%lw_derivatives, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1+shift,",",fidx2,",",this%lw_derivatives(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

  end subroutine debug_dump
  
  !---------------------------------------------------------------------
  ! debug dump
  subroutine debug_dump_sw(this,FILE,LINE,istartcol,iendcol)
    use radiation_io,     only : nulout
    type(flux_type), intent(in) :: this
    character, intent(in) :: FILE
    integer, intent(in) :: LINE
    integer, intent(in) :: istartcol, iendcol
    integer :: file_idx, fidx1, fidx2, fidx3
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g))
    !$OMP TARGET UPDATE FROM(this%sw_up_clear) IF(allocated(this%sw_up_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_clear) IF(allocated(this%sw_dn_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct) IF(allocated(this%sw_dn_direct))
    !$OMP TARGET UPDATE FROM(this%sw_up) IF(allocated(this%sw_up))
    !$OMP TARGET UPDATE FROM(this%sw_dn) IF(allocated(this%sw_dn))

    !$ACC UPDATE HOST(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g))  ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_clear) IF(allocated(this%sw_up_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_clear) IF(allocated(this%sw_dn_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct) IF(allocated(this%sw_dn_direct)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up) IF(allocated(this%sw_up)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn) IF(allocated(this%sw_dn)) ASYNC(1)

    !$ACC WAIT(1)
    write(nulout,'(a,a,a,i0)')         "    DUMPING DEVICE ARRAYS FROM ", FILE, " : LINE = ", LINE
    write(nulout,'(a,a,a,i0)')         "                            IN ", __FILE__, " : LINE = ", __LINE__

    if (allocated(this%cloud_cover_sw)) then
       open(newunit=file_idx, file='cloud_cover_sw.txt', status='replace')
       do fidx1 = istartcol, iendcol
          write(file_idx, '(i0,a,es13.6)') fidx1,",",this%cloud_cover_sw(fidx1)
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_dn_diffuse_surf_clear_g)) then
       open(newunit=file_idx, file='sw_dn_diffuse_surf_clear_g.txt', status='replace')
       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%sw_dn_diffuse_surf_clear_g, dim=1), ubound(this%sw_dn_diffuse_surf_clear_g, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_dn_diffuse_surf_clear_g(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%sw_dn_direct_surf_clear_g)) then
       open(newunit=file_idx, file='sw_dn_direct_surf_clear_g.txt', status='replace')
       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%sw_dn_direct_surf_clear_g, dim=1), ubound(this%sw_dn_direct_surf_clear_g, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_dn_direct_surf_clear_g(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_dn_diffuse_surf_g)) then
       open(newunit=file_idx, file='sw_dn_diffuse_surf_g.txt', status='replace')
       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%sw_dn_diffuse_surf_g, dim=1), ubound(this%sw_dn_diffuse_surf_g, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_dn_diffuse_surf_g(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%sw_dn_direct_surf_g)) then
       open(newunit=file_idx, file='sw_dn_direct_surf_g.txt', status='replace')
       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%sw_dn_direct_surf_g, dim=1), ubound(this%sw_dn_direct_surf_g, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_dn_direct_surf_g(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_up_clear)) then
       open(newunit=file_idx, file='sw_up_clear.txt', status='replace')
       do fidx2 = lbound(this%sw_up_clear, dim=2), ubound(this%sw_up_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_up_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%sw_dn_clear)) then
       open(newunit=file_idx, file='sw_dn_clear.txt', status='replace')
       do fidx2 = lbound(this%sw_dn_clear, dim=2), ubound(this%sw_dn_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_dn_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%sw_up)) then
       open(newunit=file_idx, file='sw_up.txt', status='replace')
       do fidx2 = lbound(this%sw_up, dim=2), ubound(this%sw_up, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_up(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%sw_dn)) then
       open(newunit=file_idx, file='sw_dn.txt', status='replace')
       do fidx2 = lbound(this%sw_dn, dim=2), ubound(this%sw_dn, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_dn(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%sw_dn_direct_clear)) then
       open(newunit=file_idx, file='sw_dn_direct_clear.txt', status='replace')
       do fidx2 = lbound(this%sw_dn_direct_clear, dim=2), ubound(this%sw_dn_direct_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_dn_direct_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%sw_dn_direct)) then
       open(newunit=file_idx, file='sw_dn_direct.txt', status='replace')
       do fidx2 = lbound(this%sw_dn_direct, dim=2), ubound(this%sw_dn_direct, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%sw_dn_direct(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

  end subroutine debug_dump_sw
      
  !---------------------------------------------------------------------
  ! debug dump
  subroutine debug_dump_lw(this,FILE,LINE,istartcol,iendcol)
    use radiation_io,     only : nulout
    type(flux_type), intent(in) :: this
    character, intent(in) :: FILE
    integer, intent(in) :: LINE
    integer, intent(in) :: istartcol, iendcol
    integer :: file_idx, fidx1, fidx2, fidx3
    !$OMP TARGET UPDATE FROM(this%lw_up) IF(allocated(this%lw_up))
    !$OMP TARGET UPDATE FROM(this%lw_dn) IF(allocated(this%lw_dn))
    !$OMP TARGET UPDATE FROM(this%lw_up_clear) IF(allocated(this%lw_up_clear))
    !$OMP TARGET UPDATE FROM(this%lw_dn_clear) IF(allocated(this%lw_dn_clear))
    !$OMP TARGET UPDATE FROM(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw))
    !$OMP TARGET UPDATE FROM(this%lw_derivatives) IF(allocated(this%lw_derivatives))
    !$OMP TARGET UPDATE FROM(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g))
    !$OMP TARGET UPDATE FROM(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g))

    !$ACC UPDATE HOST(this%lw_up) IF(allocated(this%lw_up)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn) IF(allocated(this%lw_dn)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_clear) IF(allocated(this%lw_up_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_clear) IF(allocated(this%lw_dn_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_derivatives) IF(allocated(this%lw_derivatives)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g)) ASYNC(1)

    !$ACC WAIT(1)
    write(nulout,'(a,a,a,i0)')         "    DUMPING DEVICE ARRAYS FROM ", FILE, " : LINE = ", LINE
    write(nulout,'(a,a,a,i0)')         "                            IN ", __FILE__, " : LINE = ", __LINE__

    if (allocated(this%lw_up)) then
       open(newunit=file_idx, file='lw_up.txt', status='replace')
       do fidx2 = lbound(this%lw_up, dim=2), ubound(this%lw_up, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%lw_up(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%lw_dn)) then
       open(newunit=file_idx, file='lw_dn.txt', status='replace')
       do fidx2 = lbound(this%lw_dn, dim=2), ubound(this%lw_dn, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%lw_dn(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%lw_up_clear)) then
       open(newunit=file_idx, file='lw_up_clear.txt', status='replace')
       do fidx2 = lbound(this%lw_up_clear, dim=2), ubound(this%lw_up_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%lw_up_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    
    if (allocated(this%lw_dn_clear)) then
       open(newunit=file_idx, file='lw_dn_clear.txt', status='replace')
       do fidx2 = lbound(this%lw_dn_clear, dim=2), ubound(this%lw_dn_clear, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%lw_dn_clear(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%lw_dn_surf_g)) then
       open(newunit=file_idx, file='lw_dn_surf_g.txt', status='replace')
       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%lw_dn_surf_g, dim=1), ubound(this%lw_dn_surf_g, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%lw_dn_surf_g(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if
    if (allocated(this%lw_dn_surf_clear_g)) then
       open(newunit=file_idx, file='lw_dn_surf_clear_g.txt', status='replace')
       do fidx2 = istartcol, iendcol
          do fidx1 = lbound(this%lw_dn_surf_clear_g, dim=1), ubound(this%lw_dn_surf_clear_g, dim=1)
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%lw_dn_surf_clear_g(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

    if (allocated(this%cloud_cover_lw)) then
       open(newunit=file_idx, file='cloud_cover_lw.txt', status='replace')
       do fidx1 = istartcol, iendcol
          write(file_idx, '(i0,a,es13.6)') fidx1,",",this%cloud_cover_lw(fidx1)
       enddo
       close(file_idx)
    end if
    if (allocated(this%lw_derivatives)) then
       open(newunit=file_idx, file='lw_derivatives.txt', status='replace')
       do fidx2 = lbound(this%lw_derivatives, dim=2), ubound(this%lw_derivatives, dim=2)
          do fidx1 = istartcol, iendcol
             write(file_idx, '(i0,a,i0,a,es13.6)') fidx1,",",fidx2,",",this%lw_derivatives(fidx1,fidx2)
          enddo
       enddo
       close(file_idx)
    end if

  end subroutine debug_dump_lw

  !---------------------------------------------------------------------
  ! debug print
  subroutine debug_print_flux(this,FILE,LINE,istartcol,iendcol)
    use radiation_io,     only : nulout
    type(flux_type), intent(in) :: this
    character, intent(in) :: FILE
    integer, intent(in) :: LINE
    integer, intent(in) :: istartcol, iendcol
    integer :: s1, s2, s3

    !$OMP TARGET UPDATE FROM(this%lw_up) IF(allocated(this%lw_up))
    !$OMP TARGET UPDATE FROM(this%lw_dn) IF(allocated(this%lw_dn))
    !$OMP TARGET UPDATE FROM(this%lw_up_clear) IF(allocated(this%lw_up_clear))
    !$OMP TARGET UPDATE FROM(this%lw_dn_clear) IF(allocated(this%lw_dn_clear))
    !$OMP TARGET UPDATE FROM(this%sw_up) IF(allocated(this%sw_up))
    !$OMP TARGET UPDATE FROM(this%sw_dn) IF(allocated(this%sw_dn))
    !$OMP TARGET UPDATE FROM(this%sw_up_clear) IF(allocated(this%sw_up_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_clear) IF(allocated(this%sw_dn_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct) IF(allocated(this%sw_dn_direct))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear))
    !$OMP TARGET UPDATE FROM(this%lw_up_band) IF(allocated(this%lw_up_band))
    !$OMP TARGET UPDATE FROM(this%lw_dn_band) IF(allocated(this%lw_dn_band))
    !$OMP TARGET UPDATE FROM(this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band))
    !$OMP TARGET UPDATE FROM(this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_up_band) IF(allocated(this%sw_up_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_band) IF(allocated(this%sw_dn_band))
    !$OMP TARGET UPDATE FROM(this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band))
    !$OMP TARGET UPDATE FROM(this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy))
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy))
    !$OMP TARGET UPDATE FROM(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw))
    !$OMP TARGET UPDATE FROM(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw))
    !$OMP TARGET UPDATE FROM(this%lw_derivatives) IF(allocated(this%lw_derivatives))
    !$OMP TARGET UPDATE FROM(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g))
    !$OMP TARGET UPDATE FROM(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g))
    !$OMP TARGET UPDATE FROM(this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g))
    !$OMP TARGET UPDATE FROM(this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g))
    !$OMP TARGET UPDATE FROM(this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g))
    !$OMP TARGET UPDATE FROM(this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g))
    !$OMP TARGET UPDATE FROM(this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band))
    !$OMP TARGET UPDATE FROM(this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band))
    !$OMP TARGET UPDATE FROM(this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band))
    !$OMP TARGET UPDATE FROM(this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band))

    !$ACC UPDATE HOST(this%lw_up) IF(allocated(this%lw_up)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn) IF(allocated(this%lw_dn)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_clear) IF(allocated(this%lw_up_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_clear) IF(allocated(this%lw_dn_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up) IF(allocated(this%sw_up)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn) IF(allocated(this%sw_dn)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_clear) IF(allocated(this%sw_up_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_clear) IF(allocated(this%sw_dn_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct) IF(allocated(this%sw_dn_direct)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_band) IF(allocated(this%lw_up_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_band) IF(allocated(this%lw_dn_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_band) IF(allocated(this%sw_up_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_band) IF(allocated(this%sw_dn_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy)) ASYNC(1)
    !$ACC UPDATE HOST(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw)) ASYNC(1)
    !$ACC UPDATE HOST(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_derivatives) IF(allocated(this%lw_derivatives)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band)) ASYNC(1)

    !$ACC WAIT(1)

    write(nulout,'(a,a,a,i0)')         "    DEVICE ARRAY DEBUG FROM ", FILE, " : LINE = ", LINE
    write(nulout,'(a,a,a,i0)')         "                         IN ", __FILE__, " : LINE = ", __LINE__

    if (allocated(this%lw_up)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%lw_up,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_up=", this%lw_up(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_dn)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%lw_dn,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_dn=", this%lw_dn(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_up_clear)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%lw_up_clear,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_up_clear=", this%lw_up_clear(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_dn_clear)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%lw_dn_clear,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_dn_clear=", this%lw_dn_clear(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_up)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%sw_up,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_up=", this%sw_up(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%sw_dn,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn=", this%sw_dn(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_up_clear)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%sw_up_clear,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_up_clear=", this%sw_up_clear(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_clear)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%sw_dn_clear,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_clear=", this%sw_dn_clear(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_direct)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%sw_dn_direct,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_direct=", this%sw_dn_direct(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_direct_clear)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%sw_dn_direct_clear,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_direct_clear=", this%sw_dn_direct_clear(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_up_band)) then
       s1 = size(this%lw_up_band,1)/2
       s2 = size(this%lw_up_band,2)/2
       s3 = size(this%lw_up_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%lw_up_band=", this%lw_up_band(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%lw_dn_band)) then
       s1 = size(this%lw_dn_band,1)/2
       s2 = size(this%lw_dn_band,2)/2
       s3 = size(this%lw_dn_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%lw_dn_band=", this%lw_dn_band(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%sw_up_band)) then
       s1 = size(this%sw_up_band,1)/2
       s2 = size(this%sw_up_band,2)/2
       s3 = size(this%sw_up_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%sw_up_band=", this%sw_up_band(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%sw_dn_band)) then
       s1 = size(this%sw_dn_band,1)/2
       s2 = size(this%sw_dn_band,2)/2
       s3 = size(this%sw_dn_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%sw_dn_band=", this%sw_dn_band(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%sw_dn_direct_band)) then
       s1 = size(this%sw_dn_direct_band,1)/2
       s2 = size(this%sw_dn_direct_band,2)/2
       s3 = size(this%sw_dn_direct_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%sw_dn_direct_band=", this%sw_dn_direct_band(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%lw_up_clear_band)) then
       s1 = size(this%lw_up_clear_band,1)/2
       s2 = size(this%lw_up_clear_band,2)/2
       s3 = size(this%lw_up_clear_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%lw_up_clear_band=", this%lw_up_clear_band(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%lw_dn_clear_band)) then
       s1 = size(this%lw_dn_clear_band,1)/2
       s2 = size(this%lw_dn_clear_band,2)/2
       s3 = size(this%lw_dn_clear_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%lw_dn_clear_band=", this%lw_dn_clear_band(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%sw_dn_direct_clear_band)) then
       s1 = size(this%sw_dn_direct_clear_band,1)/2
       s2 = size(this%sw_dn_direct_clear_band,2)/2
       s3 = size(this%sw_dn_direct_clear_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%sw_dn_direct_clear_band=", this%sw_dn_direct_clear_band(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%lw_dn_surf_g)) then
       s1 = size(this%lw_dn_surf_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_dn_surf_g=", this%lw_dn_surf_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_dn_surf_clear_g)) then
       s1 = size(this%lw_dn_surf_clear_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_dn_surf_clear_g=", this%lw_dn_surf_clear_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_diffuse_surf_g)) then
       s1 = size(this%sw_dn_diffuse_surf_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_diffuse_surf_g=", this%sw_dn_diffuse_surf_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_direct_surf_g)) then
       s1 = size(this%sw_dn_direct_surf_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_direct_surf_g=", this%sw_dn_direct_surf_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_diffuse_surf_clear_g)) then
       s1 = size(this%sw_dn_diffuse_surf_clear_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_diffuse_surf_clear_g=", this%sw_dn_diffuse_surf_clear_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_direct_surf_clear_g)) then
       s1 = size(this%sw_dn_direct_surf_clear_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_direct_surf_clear_g=", this%sw_dn_direct_surf_clear_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_toa_g)) then
       s1 = size(this%sw_dn_toa_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_toa_g=", this%sw_dn_toa_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_up_toa_g)) then
       s1 = size(this%lw_up_toa_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_up_toa_g=", this%lw_up_toa_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_up_toa_clear_g)) then
       s1 = size(this%lw_up_toa_clear_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_up_toa_clear_g=", this%lw_up_toa_clear_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_up_toa_g)) then
       s1 = size(this%sw_up_toa_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_up_toa_g=", this%sw_up_toa_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_up_toa_clear_g)) then
       s1 = size(this%sw_up_toa_clear_g,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_up_toa_clear_g=", this%sw_up_toa_clear_g(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_surf_band)) then
       s1 = size(this%sw_dn_surf_band,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_surf_band=", this%sw_dn_surf_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_direct_surf_band)) then
       s1 = size(this%sw_dn_direct_surf_band,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_direct_surf_band=", this%sw_dn_direct_surf_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_surf_clear_band)) then
       s1 = size(this%sw_dn_surf_clear_band,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_surf_clear_band=", this%sw_dn_surf_clear_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_direct_surf_clear_band)) then
       s1 = size(this%sw_dn_direct_surf_clear_band,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_direct_surf_clear_band=", this%sw_dn_direct_surf_clear_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_toa_band)) then
       s1 = size(this%sw_dn_toa_band,1)/2
       s2 = size(this%sw_dn_toa_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_toa_band=", this%sw_dn_toa_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_up_toa_band)) then
       s1 = size(this%lw_up_toa_band,1)/2
       s2 = size(this%lw_up_toa_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_up_toa_band=", this%lw_up_toa_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_up_toa_clear_band)) then
       s1 = size(this%lw_up_toa_clear_band,1)/2
       s2 = size(this%lw_up_toa_clear_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_up_toa_clear_band=", this%lw_up_toa_clear_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_up_toa_band)) then
       s1 = size(this%sw_up_toa_band,1)/2
       s2 = size(this%sw_up_toa_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_up_toa_band=", this%sw_up_toa_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_up_toa_clear_band)) then
       s1 = size(this%lw_up_toa_clear_band,1)/2
       s2 = size(this%lw_up_toa_clear_band,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_up_toa_clear_band=", this%lw_up_toa_clear_band(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%cloud_cover_lw)) then
       s1 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0)') "this%cloud_cover_lw=", this%cloud_cover_lw(s1), " at indices ", s1
    end if
    if (allocated(this%cloud_cover_sw)) then
       s1 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0)') "this%cloud_cover_sw=", this%cloud_cover_sw(s1), " at indices ", s1
    end if
    if (allocated(this%lw_dn_surf_canopy)) then
       s1 = size(this%lw_dn_surf_canopy,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_dn_surf_canopy=", this%lw_dn_surf_canopy(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_diffuse_surf_canopy)) then
       s1 = size(this%sw_dn_diffuse_surf_canopy,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_diffuse_surf_canopy=", this%sw_dn_diffuse_surf_canopy(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%sw_dn_direct_surf_canopy)) then
       s1 = size(this%sw_dn_direct_surf_canopy,1)/2
       s2 = (iendcol-istartcol)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%sw_dn_direct_surf_canopy=", this%sw_dn_direct_surf_canopy(s1,s2), " at indices ", s1, " ", s2
    end if
    if (allocated(this%lw_derivatives)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%lw_derivatives,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0)') "this%lw_derivatives=", this%lw_derivatives(s1,s2), " at indices ", s1, " ", s2
    end if
  end subroutine debug_print_flux

  !---------------------------------------------------------------------
  ! Creates fields on device
  subroutine create_device_flux(this)

    type(flux_type), intent(inout) :: this
#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_up) IF(allocated(this%lw_up))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_dn) IF(allocated(this%lw_dn))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_up_clear) IF(allocated(this%lw_up_clear))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_dn_clear) IF(allocated(this%lw_dn_clear))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_up) IF(allocated(this%sw_up))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn) IF(allocated(this%sw_dn))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_up_clear) IF(allocated(this%sw_up_clear))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_clear) IF(allocated(this%sw_dn_clear))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct) IF(allocated(this%sw_dn_direct))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_up_band) IF(allocated(this%lw_up_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_dn_band) IF(allocated(this%lw_dn_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_up_band) IF(allocated(this%sw_up_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_band) IF(allocated(this%sw_dn_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_derivatives) IF(allocated(this%lw_derivatives))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band))

    !$ACC ENTER DATA CREATE(this%lw_up) IF(allocated(this%lw_up)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_dn) IF(allocated(this%lw_dn)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_up_clear) IF(allocated(this%lw_up_clear)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_dn_clear) IF(allocated(this%lw_dn_clear)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_up) IF(allocated(this%sw_up)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn) IF(allocated(this%sw_dn)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_up_clear) IF(allocated(this%sw_up_clear)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_clear) IF(allocated(this%sw_dn_clear)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct) IF(allocated(this%sw_dn_direct)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_up_band) IF(allocated(this%lw_up_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_dn_band) IF(allocated(this%lw_dn_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_up_band) IF(allocated(this%sw_up_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_band) IF(allocated(this%sw_dn_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_derivatives) IF(allocated(this%lw_derivatives)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band)) ASYNC(1)

  end subroutine create_device_flux

  !---------------------------------------------------------------------
  ! updates fields on host
  subroutine update_host_flux(this)

    type(flux_type), intent(inout) :: this
#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET UPDATE FROM(this%lw_up) IF(allocated(this%lw_up))
    !$OMP TARGET UPDATE FROM(this%lw_dn) IF(allocated(this%lw_dn))
    !$OMP TARGET UPDATE FROM(this%lw_up_clear) IF(allocated(this%lw_up_clear))
    !$OMP TARGET UPDATE FROM(this%lw_dn_clear) IF(allocated(this%lw_dn_clear))
    !$OMP TARGET UPDATE FROM(this%sw_up) IF(allocated(this%sw_up))
    !$OMP TARGET UPDATE FROM(this%sw_dn) IF(allocated(this%sw_dn))
    !$OMP TARGET UPDATE FROM(this%sw_up_clear) IF(allocated(this%sw_up_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_clear) IF(allocated(this%sw_dn_clear))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct) IF(allocated(this%sw_dn_direct))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear))
    !$OMP TARGET UPDATE FROM(this%lw_up_band) IF(allocated(this%lw_up_band))
    !$OMP TARGET UPDATE FROM(this%lw_dn_band) IF(allocated(this%lw_dn_band))
    !$OMP TARGET UPDATE FROM(this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band))
    !$OMP TARGET UPDATE FROM(this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_up_band) IF(allocated(this%sw_up_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_band) IF(allocated(this%sw_dn_band))
    !$OMP TARGET UPDATE FROM(this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band))
    !$OMP TARGET UPDATE FROM(this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy))
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy))
    !$OMP TARGET UPDATE FROM(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw))
    !$OMP TARGET UPDATE FROM(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw))
    !$OMP TARGET UPDATE FROM(this%lw_derivatives) IF(allocated(this%lw_derivatives))
    !$OMP TARGET UPDATE FROM(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g))
    !$OMP TARGET UPDATE FROM(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g))
    !$OMP TARGET UPDATE FROM(this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g))
    !$OMP TARGET UPDATE FROM(this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g))
    !$OMP TARGET UPDATE FROM(this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g))
    !$OMP TARGET UPDATE FROM(this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g))
    !$OMP TARGET UPDATE FROM(this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g))
    !$OMP TARGET UPDATE FROM(this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band))
    !$OMP TARGET UPDATE FROM(this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band))
    !$OMP TARGET UPDATE FROM(this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band))
    !$OMP TARGET UPDATE FROM(this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band))
    !$OMP TARGET UPDATE FROM(this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band))

    !$ACC UPDATE HOST(this%lw_up) IF(allocated(this%lw_up)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn) IF(allocated(this%lw_dn)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_clear) IF(allocated(this%lw_up_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_clear) IF(allocated(this%lw_dn_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up) IF(allocated(this%sw_up)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn) IF(allocated(this%sw_dn)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_clear) IF(allocated(this%sw_up_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_clear) IF(allocated(this%sw_dn_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct) IF(allocated(this%sw_dn_direct)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_band) IF(allocated(this%lw_up_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_band) IF(allocated(this%lw_dn_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_band) IF(allocated(this%sw_up_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_band) IF(allocated(this%sw_dn_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy)) ASYNC(1)
    !$ACC UPDATE HOST(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw)) ASYNC(1)
    !$ACC UPDATE HOST(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_derivatives) IF(allocated(this%lw_derivatives)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band)) ASYNC(1)
    !$ACC UPDATE HOST(this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band)) ASYNC(1)

  end subroutine update_host_flux

  !---------------------------------------------------------------------
  ! updates fields on device
  subroutine update_device_flux(this)

    type(flux_type), intent(inout) :: this
#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET UPDATE TO(this%lw_up) IF(allocated(this%lw_up))
    !$OMP TARGET UPDATE TO(this%lw_dn) IF(allocated(this%lw_dn))
    !$OMP TARGET UPDATE TO(this%lw_up_clear) IF(allocated(this%lw_up_clear))
    !$OMP TARGET UPDATE TO(this%lw_dn_clear) IF(allocated(this%lw_dn_clear))
    !$OMP TARGET UPDATE TO(this%sw_up) IF(allocated(this%sw_up))
    !$OMP TARGET UPDATE TO(this%sw_dn) IF(allocated(this%sw_dn))
    !$OMP TARGET UPDATE TO(this%sw_up_clear) IF(allocated(this%sw_up_clear))
    !$OMP TARGET UPDATE TO(this%sw_dn_clear) IF(allocated(this%sw_dn_clear))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct) IF(allocated(this%sw_dn_direct))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear))
    !$OMP TARGET UPDATE TO(this%lw_up_band) IF(allocated(this%lw_up_band))
    !$OMP TARGET UPDATE TO(this%lw_dn_band) IF(allocated(this%lw_dn_band))
    !$OMP TARGET UPDATE TO(this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band))
    !$OMP TARGET UPDATE TO(this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band))
    !$OMP TARGET UPDATE TO(this%sw_up_band) IF(allocated(this%sw_up_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_band) IF(allocated(this%sw_dn_band))
    !$OMP TARGET UPDATE TO(this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band))
    !$OMP TARGET UPDATE TO(this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy))
    !$OMP TARGET UPDATE TO(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy))
    !$OMP TARGET UPDATE TO(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw))
    !$OMP TARGET UPDATE TO(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw))
    !$OMP TARGET UPDATE TO(this%lw_derivatives) IF(allocated(this%lw_derivatives))
    !$OMP TARGET UPDATE TO(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g))
    !$OMP TARGET UPDATE TO(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g))
    !$OMP TARGET UPDATE TO(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g))
    !$OMP TARGET UPDATE TO(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g))
    !$OMP TARGET UPDATE TO(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g))
    !$OMP TARGET UPDATE TO(this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g))
    !$OMP TARGET UPDATE TO(this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g))
    !$OMP TARGET UPDATE TO(this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g))
    !$OMP TARGET UPDATE TO(this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g))
    !$OMP TARGET UPDATE TO(this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g))
    !$OMP TARGET UPDATE TO(this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band))
    !$OMP TARGET UPDATE TO(this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band))
    !$OMP TARGET UPDATE TO(this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band))
    !$OMP TARGET UPDATE TO(this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band))
    !$OMP TARGET UPDATE TO(this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band))

    !$ACC UPDATE DEVICE(this%lw_up) IF(allocated(this%lw_up)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_dn) IF(allocated(this%lw_dn)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_up_clear) IF(allocated(this%lw_up_clear)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_dn_clear) IF(allocated(this%lw_dn_clear)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_up) IF(allocated(this%sw_up)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn) IF(allocated(this%sw_dn)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_up_clear) IF(allocated(this%sw_up_clear)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_clear) IF(allocated(this%sw_dn_clear)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct) IF(allocated(this%sw_dn_direct)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_up_band) IF(allocated(this%lw_up_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_dn_band) IF(allocated(this%lw_dn_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_up_band) IF(allocated(this%sw_up_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_band) IF(allocated(this%sw_dn_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_derivatives) IF(allocated(this%lw_derivatives)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band)) ASYNC(1)

  end subroutine update_device_flux

  !---------------------------------------------------------------------
  ! Deletes fields on device
  subroutine delete_device_flux(this)

    type(flux_type), intent(inout) :: this
#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_up) IF(allocated(this%lw_up))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_dn) IF(allocated(this%lw_dn))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_up_clear) IF(allocated(this%lw_up_clear))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_dn_clear) IF(allocated(this%lw_dn_clear))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_up) IF(allocated(this%sw_up))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn) IF(allocated(this%sw_dn))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_up_clear) IF(allocated(this%sw_up_clear))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_clear) IF(allocated(this%sw_dn_clear))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct) IF(allocated(this%sw_dn_direct))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_up_band) IF(allocated(this%lw_up_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_dn_band) IF(allocated(this%lw_dn_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_up_band) IF(allocated(this%sw_up_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_band) IF(allocated(this%sw_dn_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_derivatives) IF(allocated(this%lw_derivatives))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band))

    !$ACC EXIT DATA DELETE(this%lw_up) IF(allocated(this%lw_up)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_dn) IF(allocated(this%lw_dn)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_up_clear) IF(allocated(this%lw_up_clear)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_dn_clear) IF(allocated(this%lw_dn_clear)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_up) IF(allocated(this%sw_up)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn) IF(allocated(this%sw_dn)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_up_clear) IF(allocated(this%sw_up_clear)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_clear) IF(allocated(this%sw_dn_clear)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct) IF(allocated(this%sw_dn_direct)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct_clear) IF(allocated(this%sw_dn_direct_clear)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_up_band) IF(allocated(this%lw_up_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_dn_band) IF(allocated(this%lw_dn_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_up_clear_band) IF(allocated(this%lw_up_clear_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_dn_clear_band) IF(allocated(this%lw_dn_clear_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_up_band) IF(allocated(this%sw_up_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_band) IF(allocated(this%sw_dn_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_up_clear_band) IF(allocated(this%sw_up_clear_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_clear_band) IF(allocated(this%sw_dn_clear_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct_band) IF(allocated(this%sw_dn_direct_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct_clear_band) IF(allocated(this%sw_dn_direct_clear_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_surf_band) IF(allocated(this%sw_dn_surf_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct_surf_band) IF(allocated(this%sw_dn_direct_surf_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_surf_clear_band) IF(allocated(this%sw_dn_surf_clear_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct_surf_clear_band) IF(allocated(this%sw_dn_direct_surf_clear_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_dn_surf_canopy) IF(allocated(this%lw_dn_surf_canopy)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_diffuse_surf_canopy) IF(allocated(this%sw_dn_diffuse_surf_canopy)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct_surf_canopy) IF(allocated(this%sw_dn_direct_surf_canopy)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%cloud_cover_sw) IF(allocated(this%cloud_cover_sw)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%cloud_cover_lw) IF(allocated(this%cloud_cover_lw)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_derivatives) IF(allocated(this%lw_derivatives)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_dn_surf_g) IF(allocated(this%lw_dn_surf_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_dn_surf_clear_g) IF(allocated(this%lw_dn_surf_clear_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_diffuse_surf_g) IF(allocated(this%sw_dn_diffuse_surf_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct_surf_g) IF(allocated(this%sw_dn_direct_surf_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_diffuse_surf_clear_g) IF(allocated(this%sw_dn_diffuse_surf_clear_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_direct_surf_clear_g) IF(allocated(this%sw_dn_direct_surf_clear_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_up_toa_g) IF(allocated(this%lw_up_toa_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_up_toa_clear_g) IF(allocated(this%lw_up_toa_clear_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_toa_g) IF(allocated(this%sw_dn_toa_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_up_toa_g) IF(allocated(this%sw_up_toa_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_up_toa_clear_g) IF(allocated(this%sw_up_toa_clear_g)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_up_toa_band) IF(allocated(this%lw_up_toa_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%lw_up_toa_clear_band) IF(allocated(this%lw_up_toa_clear_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_dn_toa_band) IF(allocated(this%sw_dn_toa_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_up_toa_band) IF(allocated(this%sw_up_toa_band)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%sw_up_toa_clear_band) IF(allocated(this%sw_up_toa_clear_band)) ASYNC(1)

  end subroutine delete_device_flux
#endif

end module radiation_flux
