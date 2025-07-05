! radiation_memory_pool.F90 - Memory pool for radiation calculations
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
! Author:  Performance Optimization Team
! Email:   performance@ecmwf.int
!
! This module provides memory pooling for radiation calculations to reduce
! memory allocation overhead and improve cache locality.

module radiation_memory_pool

  use parkind1, only : jprb, jpim
  use yomhook,  only : lhook, dr_hook, jphook

  implicit none

  private
  public :: memory_pool_type, initialize_memory_pool, finalize_memory_pool

  ! Memory pool for radiation optical properties
  type :: memory_pool_type
    ! Pool dimensions
    integer :: n_g_lw, n_g_sw, n_bands_lw, n_bands_sw
    integer :: nlev, ncol_max
    logical :: initialized = .false.
    
    ! Longwave optical properties pool
    real(jprb), allocatable :: od_lw_pool(:,:,:)
    real(jprb), allocatable :: ssa_lw_pool(:,:,:)
    real(jprb), allocatable :: g_lw_pool(:,:,:)
    real(jprb), allocatable :: od_lw_cloud_pool(:,:,:)
    real(jprb), allocatable :: ssa_lw_cloud_pool(:,:,:)
    real(jprb), allocatable :: g_lw_cloud_pool(:,:,:)
    
    ! Shortwave optical properties pool
    real(jprb), allocatable :: od_sw_pool(:,:,:)
    real(jprb), allocatable :: ssa_sw_pool(:,:,:)
    real(jprb), allocatable :: g_sw_pool(:,:,:)
    real(jprb), allocatable :: od_sw_cloud_pool(:,:,:)
    real(jprb), allocatable :: ssa_sw_cloud_pool(:,:,:)
    real(jprb), allocatable :: g_sw_cloud_pool(:,:,:)
    
    ! Planck function and flux arrays
    real(jprb), allocatable :: planck_hl_pool(:,:,:)
    real(jprb), allocatable :: flux_up_pool(:,:,:)
    real(jprb), allocatable :: flux_dn_pool(:,:,:)
    
    ! Surface properties
    real(jprb), allocatable :: lw_emission_pool(:,:)
    real(jprb), allocatable :: lw_albedo_pool(:,:)
    real(jprb), allocatable :: sw_albedo_direct_pool(:,:)
    real(jprb), allocatable :: sw_albedo_diffuse_pool(:,:)
    real(jprb), allocatable :: incoming_sw_pool(:,:)
    
    ! GPU device pointers for OpenACC
    logical :: on_device = .false.
    
  contains
    procedure :: allocate_pool
    procedure :: deallocate_pool
    procedure :: resize_if_needed
    procedure :: get_lw_arrays
    procedure :: get_sw_arrays
    procedure :: get_flux_arrays
    procedure :: get_surface_arrays
    procedure :: create_device_arrays
    procedure :: delete_device_arrays
    procedure :: print_memory_usage
  end type memory_pool_type

  ! Global memory pool instance
  type(memory_pool_type), save :: global_pool

contains

  !---------------------------------------------------------------------
  ! Initialize the global memory pool
  subroutine initialize_memory_pool(n_g_lw, n_g_sw, n_bands_lw, n_bands_sw, &
       &                           nlev, ncol_max)
    
    integer, intent(in) :: n_g_lw, n_g_sw, n_bands_lw, n_bands_sw
    integer, intent(in) :: nlev, ncol_max
    
    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_memory_pool:initialize_memory_pool',0,hook_handle)
    
    call global_pool%allocate_pool(n_g_lw, n_g_sw, n_bands_lw, n_bands_sw, &
         &                        nlev, ncol_max)
    
    if (lhook) call dr_hook('radiation_memory_pool:initialize_memory_pool',1,hook_handle)
    
  end subroutine initialize_memory_pool

  !---------------------------------------------------------------------
  ! Finalize the global memory pool
  subroutine finalize_memory_pool()
    
    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_memory_pool:finalize_memory_pool',0,hook_handle)
    
    call global_pool%deallocate_pool()
    
    if (lhook) call dr_hook('radiation_memory_pool:finalize_memory_pool',1,hook_handle)
    
  end subroutine finalize_memory_pool

  !---------------------------------------------------------------------
  ! Allocate memory pool arrays
  subroutine allocate_pool(this, n_g_lw, n_g_sw, n_bands_lw, n_bands_sw, &
       &                  nlev, ncol_max)
    
    class(memory_pool_type), intent(inout) :: this
    integer, intent(in) :: n_g_lw, n_g_sw, n_bands_lw, n_bands_sw
    integer, intent(in) :: nlev, ncol_max
    
    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_memory_pool:allocate_pool',0,hook_handle)
    
    ! Store dimensions
    this%n_g_lw = n_g_lw
    this%n_g_sw = n_g_sw
    this%n_bands_lw = n_bands_lw
    this%n_bands_sw = n_bands_sw
    this%nlev = nlev
    this%ncol_max = ncol_max
    
    ! Allocate longwave arrays
    if (n_g_lw > 0) then
      allocate(this%od_lw_pool(n_g_lw, nlev, ncol_max))
      allocate(this%ssa_lw_pool(n_g_lw, nlev, ncol_max))
      allocate(this%g_lw_pool(n_g_lw, nlev, ncol_max))
      allocate(this%planck_hl_pool(n_g_lw, nlev+1, ncol_max))
      allocate(this%lw_emission_pool(n_g_lw, ncol_max))
      allocate(this%lw_albedo_pool(n_g_lw, ncol_max))
    end if
    
    if (n_bands_lw > 0) then
      allocate(this%od_lw_cloud_pool(n_bands_lw, nlev, ncol_max))
      allocate(this%ssa_lw_cloud_pool(n_bands_lw, nlev, ncol_max))
      allocate(this%g_lw_cloud_pool(n_bands_lw, nlev, ncol_max))
    end if
    
    ! Allocate shortwave arrays
    if (n_g_sw > 0) then
      allocate(this%od_sw_pool(n_g_sw, nlev, ncol_max))
      allocate(this%ssa_sw_pool(n_g_sw, nlev, ncol_max))
      allocate(this%g_sw_pool(n_g_sw, nlev, ncol_max))
      allocate(this%flux_up_pool(n_g_sw, nlev+1, ncol_max))
      allocate(this%flux_dn_pool(n_g_sw, nlev+1, ncol_max))
      allocate(this%sw_albedo_direct_pool(n_g_sw, ncol_max))
      allocate(this%sw_albedo_diffuse_pool(n_g_sw, ncol_max))
      allocate(this%incoming_sw_pool(n_g_sw, ncol_max))
    end if
    
    if (n_bands_sw > 0) then
      allocate(this%od_sw_cloud_pool(n_bands_sw, nlev, ncol_max))
      allocate(this%ssa_sw_cloud_pool(n_bands_sw, nlev, ncol_max))
      allocate(this%g_sw_cloud_pool(n_bands_sw, nlev, ncol_max))
    end if
    
    this%initialized = .true.
    
    ! Print memory usage information
    call this%print_memory_usage()
    
    if (lhook) call dr_hook('radiation_memory_pool:allocate_pool',1,hook_handle)
    
  end subroutine allocate_pool

  !---------------------------------------------------------------------
  ! Deallocate memory pool arrays
  subroutine deallocate_pool(this)
    
    class(memory_pool_type), intent(inout) :: this
    
    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_memory_pool:deallocate_pool',0,hook_handle)
    
    if (.not. this%initialized) return
    
    ! Delete device arrays if they exist
    if (this%on_device) then
      call this%delete_device_arrays()
    end if
    
    ! Deallocate longwave arrays
    if (allocated(this%od_lw_pool)) deallocate(this%od_lw_pool)
    if (allocated(this%ssa_lw_pool)) deallocate(this%ssa_lw_pool)
    if (allocated(this%g_lw_pool)) deallocate(this%g_lw_pool)
    if (allocated(this%od_lw_cloud_pool)) deallocate(this%od_lw_cloud_pool)
    if (allocated(this%ssa_lw_cloud_pool)) deallocate(this%ssa_lw_cloud_pool)
    if (allocated(this%g_lw_cloud_pool)) deallocate(this%g_lw_cloud_pool)
    if (allocated(this%planck_hl_pool)) deallocate(this%planck_hl_pool)
    if (allocated(this%lw_emission_pool)) deallocate(this%lw_emission_pool)
    if (allocated(this%lw_albedo_pool)) deallocate(this%lw_albedo_pool)
    
    ! Deallocate shortwave arrays
    if (allocated(this%od_sw_pool)) deallocate(this%od_sw_pool)
    if (allocated(this%ssa_sw_pool)) deallocate(this%ssa_sw_pool)
    if (allocated(this%g_sw_pool)) deallocate(this%g_sw_pool)
    if (allocated(this%od_sw_cloud_pool)) deallocate(this%od_sw_cloud_pool)
    if (allocated(this%ssa_sw_cloud_pool)) deallocate(this%ssa_sw_cloud_pool)
    if (allocated(this%g_sw_cloud_pool)) deallocate(this%g_sw_cloud_pool)
    if (allocated(this%flux_up_pool)) deallocate(this%flux_up_pool)
    if (allocated(this%flux_dn_pool)) deallocate(this%flux_dn_pool)
    if (allocated(this%sw_albedo_direct_pool)) deallocate(this%sw_albedo_direct_pool)
    if (allocated(this%sw_albedo_diffuse_pool)) deallocate(this%sw_albedo_diffuse_pool)
    if (allocated(this%incoming_sw_pool)) deallocate(this%incoming_sw_pool)
    
    this%initialized = .false.
    
    if (lhook) call dr_hook('radiation_memory_pool:deallocate_pool',1,hook_handle)
    
  end subroutine deallocate_pool

  !---------------------------------------------------------------------
  ! Resize pool if needed (grow only for performance)
  subroutine resize_if_needed(this, nlev, ncol_max)
    
    class(memory_pool_type), intent(inout) :: this
    integer, intent(in) :: nlev, ncol_max
    
    if (.not. this%initialized) return
    
    ! Only resize if we need more space
    if (nlev > this%nlev .or. ncol_max > this%ncol_max) then
      ! Deallocate and reallocate with new dimensions
      call this%deallocate_pool()
      call this%allocate_pool(this%n_g_lw, this%n_g_sw, this%n_bands_lw, &
           &                  this%n_bands_sw, max(nlev, this%nlev), &
           &                  max(ncol_max, this%ncol_max))
    end if
    
  end subroutine resize_if_needed

  !---------------------------------------------------------------------
  ! Get longwave arrays from pool
  subroutine get_lw_arrays(this, istartcol, iendcol, nlev, &
       &                   od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &                   planck_hl, lw_emission, lw_albedo)
    
    class(memory_pool_type), intent(inout) :: this
    integer, intent(in) :: istartcol, iendcol, nlev
    real(jprb), pointer, intent(out) :: od_lw(:,:,:)
    real(jprb), pointer, intent(out) :: ssa_lw(:,:,:)
    real(jprb), pointer, intent(out) :: g_lw(:,:,:)
    real(jprb), pointer, intent(out) :: od_lw_cloud(:,:,:)
    real(jprb), pointer, intent(out) :: ssa_lw_cloud(:,:,:)
    real(jprb), pointer, intent(out) :: g_lw_cloud(:,:,:)
    real(jprb), pointer, intent(out) :: planck_hl(:,:,:)
    real(jprb), pointer, intent(out) :: lw_emission(:,:)
    real(jprb), pointer, intent(out) :: lw_albedo(:,:)
    
    ! Ensure pool is large enough
    call this%resize_if_needed(nlev, iendcol)
    
    ! Return pointers to pool arrays
    od_lw => this%od_lw_pool(:, 1:nlev, istartcol:iendcol)
    ssa_lw => this%ssa_lw_pool(:, 1:nlev, istartcol:iendcol)
    g_lw => this%g_lw_pool(:, 1:nlev, istartcol:iendcol)
    od_lw_cloud => this%od_lw_cloud_pool(:, 1:nlev, istartcol:iendcol)
    ssa_lw_cloud => this%ssa_lw_cloud_pool(:, 1:nlev, istartcol:iendcol)
    g_lw_cloud => this%g_lw_cloud_pool(:, 1:nlev, istartcol:iendcol)
    planck_hl => this%planck_hl_pool(:, 1:nlev+1, istartcol:iendcol)
    lw_emission => this%lw_emission_pool(:, istartcol:iendcol)
    lw_albedo => this%lw_albedo_pool(:, istartcol:iendcol)
    
  end subroutine get_lw_arrays

  !---------------------------------------------------------------------
  ! Get shortwave arrays from pool
  subroutine get_sw_arrays(this, istartcol, iendcol, nlev, &
       &                   od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, g_sw_cloud, &
       &                   sw_albedo_direct, sw_albedo_diffuse, incoming_sw)
    
    class(memory_pool_type), intent(inout) :: this
    integer, intent(in) :: istartcol, iendcol, nlev
    real(jprb), pointer, intent(out) :: od_sw(:,:,:)
    real(jprb), pointer, intent(out) :: ssa_sw(:,:,:)
    real(jprb), pointer, intent(out) :: g_sw(:,:,:)
    real(jprb), pointer, intent(out) :: od_sw_cloud(:,:,:)
    real(jprb), pointer, intent(out) :: ssa_sw_cloud(:,:,:)
    real(jprb), pointer, intent(out) :: g_sw_cloud(:,:,:)
    real(jprb), pointer, intent(out) :: sw_albedo_direct(:,:)
    real(jprb), pointer, intent(out) :: sw_albedo_diffuse(:,:)
    real(jprb), pointer, intent(out) :: incoming_sw(:,:)
    
    ! Ensure pool is large enough
    call this%resize_if_needed(nlev, iendcol)
    
    ! Return pointers to pool arrays
    od_sw => this%od_sw_pool(:, 1:nlev, istartcol:iendcol)
    ssa_sw => this%ssa_sw_pool(:, 1:nlev, istartcol:iendcol)
    g_sw => this%g_sw_pool(:, 1:nlev, istartcol:iendcol)
    od_sw_cloud => this%od_sw_cloud_pool(:, 1:nlev, istartcol:iendcol)
    ssa_sw_cloud => this%ssa_sw_cloud_pool(:, 1:nlev, istartcol:iendcol)
    g_sw_cloud => this%g_sw_cloud_pool(:, 1:nlev, istartcol:iendcol)
    sw_albedo_direct => this%sw_albedo_direct_pool(:, istartcol:iendcol)
    sw_albedo_diffuse => this%sw_albedo_diffuse_pool(:, istartcol:iendcol)
    incoming_sw => this%incoming_sw_pool(:, istartcol:iendcol)
    
  end subroutine get_sw_arrays

  !---------------------------------------------------------------------
  ! Get flux arrays from pool
  subroutine get_flux_arrays(this, istartcol, iendcol, nlev, &
       &                     flux_up, flux_dn)
    
    class(memory_pool_type), intent(inout) :: this
    integer, intent(in) :: istartcol, iendcol, nlev
    real(jprb), pointer, intent(out) :: flux_up(:,:,:)
    real(jprb), pointer, intent(out) :: flux_dn(:,:,:)
    
    ! Ensure pool is large enough
    call this%resize_if_needed(nlev, iendcol)
    
    ! Return pointers to pool arrays
    flux_up => this%flux_up_pool(:, 1:nlev+1, istartcol:iendcol)
    flux_dn => this%flux_dn_pool(:, 1:nlev+1, istartcol:iendcol)
    
  end subroutine get_flux_arrays

  !---------------------------------------------------------------------
  ! Get surface arrays from pool
  subroutine get_surface_arrays(this, istartcol, iendcol, &
       &                        lw_emission, lw_albedo, sw_albedo_direct, &
       &                        sw_albedo_diffuse, incoming_sw)
    
    class(memory_pool_type), intent(inout) :: this
    integer, intent(in) :: istartcol, iendcol
    real(jprb), pointer, intent(out) :: lw_emission(:,:)
    real(jprb), pointer, intent(out) :: lw_albedo(:,:)
    real(jprb), pointer, intent(out) :: sw_albedo_direct(:,:)
    real(jprb), pointer, intent(out) :: sw_albedo_diffuse(:,:)
    real(jprb), pointer, intent(out) :: incoming_sw(:,:)
    
    ! Ensure pool is large enough
    call this%resize_if_needed(1, iendcol)
    
    ! Return pointers to pool arrays
    lw_emission => this%lw_emission_pool(:, istartcol:iendcol)
    lw_albedo => this%lw_albedo_pool(:, istartcol:iendcol)
    sw_albedo_direct => this%sw_albedo_direct_pool(:, istartcol:iendcol)
    sw_albedo_diffuse => this%sw_albedo_diffuse_pool(:, istartcol:iendcol)
    incoming_sw => this%incoming_sw_pool(:, istartcol:iendcol)
    
  end subroutine get_surface_arrays

  !---------------------------------------------------------------------
  ! Create device arrays for GPU computation
  subroutine create_device_arrays(this)
    
    class(memory_pool_type), intent(inout) :: this
    
    if (.not. this%initialized) return
    if (this%on_device) return
    
#ifdef _OPENACC
    ! Create device arrays for OpenACC
    !$ACC ENTER DATA CREATE(this%od_lw_pool, this%ssa_lw_pool, this%g_lw_pool)
    !$ACC ENTER DATA CREATE(this%od_lw_cloud_pool, this%ssa_lw_cloud_pool, this%g_lw_cloud_pool)
    !$ACC ENTER DATA CREATE(this%od_sw_pool, this%ssa_sw_pool, this%g_sw_pool)
    !$ACC ENTER DATA CREATE(this%od_sw_cloud_pool, this%ssa_sw_cloud_pool, this%g_sw_cloud_pool)
    !$ACC ENTER DATA CREATE(this%planck_hl_pool, this%flux_up_pool, this%flux_dn_pool)
    !$ACC ENTER DATA CREATE(this%lw_emission_pool, this%lw_albedo_pool)
    !$ACC ENTER DATA CREATE(this%sw_albedo_direct_pool, this%sw_albedo_diffuse_pool)
    !$ACC ENTER DATA CREATE(this%incoming_sw_pool)
#endif
    
    this%on_device = .true.
    
  end subroutine create_device_arrays

  !---------------------------------------------------------------------
  ! Delete device arrays
  subroutine delete_device_arrays(this)
    
    class(memory_pool_type), intent(inout) :: this
    
    if (.not. this%on_device) return
    
#ifdef _OPENACC
    ! Delete device arrays for OpenACC
    !$ACC EXIT DATA DELETE(this%od_lw_pool, this%ssa_lw_pool, this%g_lw_pool)
    !$ACC EXIT DATA DELETE(this%od_lw_cloud_pool, this%ssa_lw_cloud_pool, this%g_lw_cloud_pool)
    !$ACC EXIT DATA DELETE(this%od_sw_pool, this%ssa_sw_pool, this%g_sw_pool)
    !$ACC EXIT DATA DELETE(this%od_sw_cloud_pool, this%ssa_sw_cloud_pool, this%g_sw_cloud_pool)
    !$ACC EXIT DATA DELETE(this%planck_hl_pool, this%flux_up_pool, this%flux_dn_pool)
    !$ACC EXIT DATA DELETE(this%lw_emission_pool, this%lw_albedo_pool)
    !$ACC EXIT DATA DELETE(this%sw_albedo_direct_pool, this%sw_albedo_diffuse_pool)
    !$ACC EXIT DATA DELETE(this%incoming_sw_pool)
#endif
    
    this%on_device = .false.
    
  end subroutine delete_device_arrays

  !---------------------------------------------------------------------
  ! Print memory usage information
  subroutine print_memory_usage(this)
    
    use radiation_io, only : nulout
    
    class(memory_pool_type), intent(in) :: this
    
    real(jprb) :: total_memory_mb
    integer :: total_elements
    
    if (.not. this%initialized) return
    
    ! Calculate total memory usage
    total_elements = 0
    
    ! Longwave arrays
    if (this%n_g_lw > 0) then
      total_elements = total_elements + this%n_g_lw * this%nlev * this%ncol_max * 3  ! od, ssa, g
      total_elements = total_elements + this%n_g_lw * (this%nlev + 1) * this%ncol_max  ! planck_hl
      total_elements = total_elements + this%n_g_lw * this%ncol_max * 2  ! emission, albedo
    end if
    
    if (this%n_bands_lw > 0) then
      total_elements = total_elements + this%n_bands_lw * this%nlev * this%ncol_max * 3  ! cloud arrays
    end if
    
    ! Shortwave arrays
    if (this%n_g_sw > 0) then
      total_elements = total_elements + this%n_g_sw * this%nlev * this%ncol_max * 3  ! od, ssa, g
      total_elements = total_elements + this%n_g_sw * (this%nlev + 1) * this%ncol_max * 2  ! flux arrays
      total_elements = total_elements + this%n_g_sw * this%ncol_max * 3  ! albedo, incoming
    end if
    
    if (this%n_bands_sw > 0) then
      total_elements = total_elements + this%n_bands_sw * this%nlev * this%ncol_max * 3  ! cloud arrays
    end if
    
    ! Convert to MB (assuming 8 bytes per real(jprb))
    total_memory_mb = real(total_elements * 8) / (1024.0_jprb * 1024.0_jprb)
    
    write(nulout, '(a)') 'ECRAD Memory Pool Summary:'
    write(nulout, '(a,i0,a,i0,a,i0,a,i0)') '  Dimensions: n_g_lw=', this%n_g_lw, &
         &  ', n_g_sw=', this%n_g_sw, ', nlev=', this%nlev, ', ncol_max=', this%ncol_max
    write(nulout, '(a,i0)') '  Total elements: ', total_elements
    write(nulout, '(a,f0.2,a)') '  Total memory: ', total_memory_mb, ' MB'
    write(nulout, '(a,l1)') '  On device: ', this%on_device
    
  end subroutine print_memory_usage

end module radiation_memory_pool