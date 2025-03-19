! radiation_aerosol.F90 - Derived type describing aerosol
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
!   2018-04-15  R. Hogan  Add "direct" option
!   2019-01-14  R. Hogan  Added out_of_physical_bounds routine

module radiation_aerosol

  use parkind1, only : jprb
  use radiation_io, only : nulerr, radiation_abort, nulout

  implicit none
  public

  !---------------------------------------------------------------------
  ! Type describing the aerosol content in the atmosphere
  type aerosol_type
     ! The mass mixing ratio of config%n_aerosol_types different
     ! aerosol types dimensioned
     ! (ncol,istartlev:iendlev,config%n_aerosol_types), where ncol is
     ! the number of columns, istartlev:iendlev is the range of model
     ! levels where aerosols are present
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mixing_ratio  ! mass mixing ratio (kg/kg)

     ! Alternatively, if is_direct=true, the optical properties are
     ! provided directly and are dimensioned
     ! (nband,istartlev:iendlev,ncol)
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  od_sw, ssa_sw, g_sw, & ! Shortwave optical properties
          &  od_lw, ssa_lw, g_lw    ! Longwave optical properties

     ! Range of levels in which the aerosol properties are provided
     integer :: istartlev, iendlev

     ! Are the optical properties going to be provided directly by the
     ! user?
     logical :: is_direct = .false.

   contains
      procedure :: allocate        => allocate_aerosol_arrays
      procedure :: allocate_direct => allocate_aerosol_arrays_direct
      procedure :: deallocate      => deallocate_aerosol_arrays
      procedure :: out_of_physical_bounds
#if defined(_OPENACC)  || defined(OMPGPU)
      procedure, nopass :: create_device => create_device_aerosol
      procedure, nopass :: update_host   => update_host_aerosol
      procedure, nopass :: update_device => update_device_aerosol
      procedure, nopass :: delete_device => delete_device_aerosol
#endif
  end type aerosol_type

contains

  !---------------------------------------------------------------------
  ! Allocate array for describing aerosols, although in the offline
  ! code these are allocated when they are read from the NetCDF file
  subroutine allocate_aerosol_arrays(this, ncol, istartlev, iendlev, ntype)

    use yomhook,     only : lhook, dr_hook, jphook

    class(aerosol_type), intent(inout) :: this
    integer, intent(in)                :: ncol  ! Number of columns
    integer, intent(in)                :: istartlev, iendlev ! Level range
    integer, intent(in)                :: ntype ! Number of aerosol types
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:allocate',0,hook_handle)

    allocate(this%mixing_ratio(ncol,istartlev:iendlev,ntype))
    this%is_direct = .false.
    this%istartlev = istartlev
    this%iendlev   = iendlev


    if (lhook) call dr_hook('radiation_aerosol:allocate',1,hook_handle)

  end subroutine allocate_aerosol_arrays


  !---------------------------------------------------------------------
  ! Allocate arrays for describing aerosol optical properties
  subroutine allocate_aerosol_arrays_direct(this, config, &
       &                                    ncol, istartlev, iendlev)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_config, only : config_type

    class(aerosol_type), intent(inout) :: this
    type(config_type),   intent(in)    :: config
    integer, intent(in)                :: ncol  ! Number of columns
    integer, intent(in)                :: istartlev, iendlev ! Level range
    integer                            :: jband, jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:allocate_direct',0,hook_handle)

    this%is_direct = .true.
    this%istartlev = istartlev
    this%iendlev   = iendlev

    if (config%do_sw) then
      allocate(this%od_sw (config%n_bands_sw,istartlev:iendlev,ncol))
      allocate(this%ssa_sw(config%n_bands_sw,istartlev:iendlev,ncol))
      allocate(this%g_sw  (config%n_bands_sw,istartlev:iendlev,ncol))
    end if

    if (config%do_lw) then
      allocate(this%od_lw (config%n_bands_lw,istartlev:iendlev,ncol))
      allocate(this%ssa_lw(config%n_bands_lw,istartlev:iendlev,ncol))
      allocate(this%g_lw  (config%n_bands_lw,istartlev:iendlev,ncol))

      ! for openacc, this is done during create_device
#if defined(_OPENACC) || defined(OMPGPU)
#else
      ! If longwave scattering by aerosol is not to be represented,
      ! then the user may wish to just provide absorption optical
      ! depth in od_lw, in which case we must set the following two
      ! variables to zero

      ! !$ACC WAIT ! ACCWA (nvhpc 22.7) crashes otherwise

      ! !$ACC PARALLEL DEFAULT(NONE) PRESENT(this, config) ASYNC(1)
      ! !$ACC LOOP GANG VECTOR COLLAPSE(3)
      do jcol = 1,ncol
        do jlev = istartlev,iendlev
          do jband = 1,config%n_bands_lw
            this%ssa_lw(jband,jlev,jcol) = 0.0_jprb
            this%g_lw(jband,jlev,jcol) = 0.0_jprb
          end do
        end do
      end do
      ! !$ACC END PARALLEL
#endif
    end if

    if (lhook) call dr_hook('radiation_aerosol:allocate_direct',1,hook_handle)

  end subroutine allocate_aerosol_arrays_direct


  !---------------------------------------------------------------------
  ! Deallocate array
  subroutine deallocate_aerosol_arrays(this)

    use yomhook,     only : lhook, dr_hook, jphook

    class(aerosol_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:deallocate',0,hook_handle)

    if (allocated(this%mixing_ratio)) deallocate(this%mixing_ratio)
    if (allocated(this%od_sw))        deallocate(this%od_sw)
    if (allocated(this%ssa_sw))       deallocate(this%ssa_sw)
    if (allocated(this%g_sw))         deallocate(this%g_sw)
    if (allocated(this%od_lw))        deallocate(this%od_lw)
    if (allocated(this%ssa_lw))       deallocate(this%ssa_lw)
    if (allocated(this%g_lw))         deallocate(this%g_lw)

    if (lhook) call dr_hook('radiation_aerosol:deallocate',1,hook_handle)

  end subroutine deallocate_aerosol_arrays


  !---------------------------------------------------------------------
  ! Return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol, do_fix) result(is_bad)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_3d

    class(aerosol_type),   intent(inout) :: this
    integer,      optional,intent(in) :: istartcol, iendcol
    logical,      optional,intent(in) :: do_fix
    logical                           :: is_bad

    logical    :: do_fix_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:out_of_physical_bounds',0,hook_handle)

    if (present(do_fix)) then
      do_fix_local = do_fix
    else
      do_fix_local = .false.
    end if

    is_bad =    out_of_bounds_3d(this%mixing_ratio, 'aerosol%mixing_ratio', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_3d(this%od_sw, 'aerosol%od_sw', &
         &                       0.0_jprb, 100.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%od_lw, 'aerosol%od_lw', &
         &                       0.0_jprb, 100.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%ssa_sw, 'aerosol%ssa_sw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%ssa_lw, 'aerosol%ssa_lw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%g_sw, 'aerosol%g_sw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%g_lw, 'aerosol%g_lw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol)

    if (lhook) call dr_hook('radiation_aerosol:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds

#if defined(_OPENACC) || defined(OMPGPU)
  !---------------------------------------------------------------------
  ! debug print
  subroutine debug_print_aerosol(this,FILE,LINE,istartcol,iendcol)
    use radiation_io,     only : nulout
    type(aerosol_type), intent(in) :: this
    character, intent(in) :: FILE
    integer, intent(in) :: LINE
    integer, intent(in) :: istartcol, iendcol
    integer :: s1, s2, s3
    
    !$OMP TARGET UPDATE FROM(this%mixing_ratio, this%od_sw, this%ssa_sw, this%g_sw, this%od_lw, this%ssa_lw, this%g_lw)

    !$ACC UPDATE HOST(this%mixing_ratio) ASYNC(1) IF(allocated(this%mixing_ratio))
    !$ACC UPDATE HOST(this%od_sw) ASYNC(1) IF(allocated(this%od_sw))
    !$ACC UPDATE HOST(this%ssa_sw) ASYNC(1) IF(allocated(this%ssa_sw))
    !$ACC UPDATE HOST(this%g_sw) ASYNC(1) IF(allocated(this%g_sw))
    !$ACC UPDATE HOST(this%od_lw) ASYNC(1) IF(allocated(this%od_lw))
    !$ACC UPDATE HOST(this%ssa_lw) ASYNC(1) IF(allocated(this%ssa_lw))
    !$ACC UPDATE HOST(this%g_lw) ASYNC(1) IF(allocated(this%g_lw))

    !$ACC WAIT(1)

    write(nulout,'(a,a,a,i0)')         "    DEVICE ARRAY DEBUG FROM ", FILE, " : LINE = ", LINE
    write(nulout,'(a,a,a,i0)')         "                         IN ", __FILE__, " : LINE = ", __LINE__
    
    if (allocated(this%mixing_ratio)) then
       s1 = (iendcol-istartcol)/2
       s2 = size(this%mixing_ratio,2)/2
       s3 = size(this%mixing_ratio,3)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%mixing_ratio=", this%mixing_ratio(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    s3 = iendcol
    if (allocated(this%od_sw)) then
       s1 = size(this%od_sw,1)/2
       s2 = size(this%od_sw,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%od_sw=", this%od_sw(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%ssa_sw)) then
       s1 = size(this%ssa_sw,1)/2
       s2 = size(this%ssa_sw,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%ssa_sw=", this%ssa_sw(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%g_sw)) then
       s1 = size(this%g_sw,1)/2
       s2 = size(this%g_sw,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%g_sw=", this%g_sw(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%od_lw)) then
       s1 = size(this%od_lw,1)/2
       s2 = size(this%od_lw,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%od_lw=", this%od_lw(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%ssa_lw)) then
       s1 = size(this%ssa_lw,1)/2
       s2 = size(this%ssa_lw,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%ssa_lw=", this%ssa_lw(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
    if (allocated(this%g_lw)) then
       s1 = size(this%g_lw,1)/2
       s2 = size(this%g_lw,2)/2
       write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "this%g_lw=", this%g_lw(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    end if
  end subroutine debug_print_aerosol

  !---------------------------------------------------------------------
  ! Creates fields on device
  subroutine create_device_aerosol(this)

    type(aerosol_type), intent(inout) :: this
#if defined(OMPGPU)
    integer :: i,j,k
#endif
#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET ENTER DATA MAP(ALLOC:this%mixing_ratio) IF(allocated(this%mixing_ratio))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%od_sw) IF(allocated(this%od_sw))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%ssa_sw) IF(allocated(this%ssa_sw))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%g_sw) IF(allocated(this%g_sw))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%od_lw) IF(allocated(this%od_lw))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%ssa_lw) IF(allocated(this%ssa_lw))
    !$OMP TARGET ENTER DATA MAP(ALLOC:this%g_lw) IF(allocated(this%g_lw))

    !$ACC ENTER DATA CREATE(this%mixing_ratio) IF(allocated(this%mixing_ratio)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%od_sw) IF(allocated(this%od_sw)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%ssa_sw) IF(allocated(this%ssa_sw)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%g_sw) IF(allocated(this%g_sw)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%od_lw) IF(allocated(this%od_lw)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%ssa_lw) IF(allocated(this%ssa_lw)) ASYNC(1)
    !$ACC ENTER DATA CREATE(this%g_lw) IF(allocated(this%g_lw)) ASYNC(1)

    ! (note that both should always be allocated)
    if (allocated(this%ssa_lw) .or. allocated(this%g_lw)) then
      ! If longwave scattering by aerosol is not to be represented,
      ! then the user may wish to just provide absorption optical
      ! depth in od_lw, in which case we must set the following two
      ! variables to zero
#if defined(_OPENACC)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
      this%ssa_lw(:,:,:) = 0.0_jprb
      this%g_lw(:,:,:) = 0.0_jprb
      !$ACC END KERNELS
#endif
#if defined(OMPGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3)
      do k = 1,SIZE(this%ssa_lw, 3)
        do j = 1,SIZE(this%ssa_lw, 2)
          do i = 1,SIZE(this%ssa_lw, 1)
            this%ssa_lw(i,j,k) = 0.0_jprb
            this%g_lw(i,j,k) = 0.0_jprb
          end do
        end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#endif
    endif

  end subroutine create_device_aerosol

  !---------------------------------------------------------------------
  ! updates fields on host
  subroutine update_host_aerosol(this)

    type(aerosol_type), intent(inout) :: this
#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET UPDATE FROM(this%mixing_ratio) IF(allocated(this%mixing_ratio))
    !$OMP TARGET UPDATE FROM(this%od_sw) IF(allocated(this%od_sw))
    !$OMP TARGET UPDATE FROM(this%ssa_sw) IF(allocated(this%ssa_sw))
    !$OMP TARGET UPDATE FROM(this%g_sw) IF(allocated(this%g_sw))
    !$OMP TARGET UPDATE FROM(this%od_lw) IF(allocated(this%od_lw))
    !$OMP TARGET UPDATE FROM(this%ssa_lw) IF(allocated(this%ssa_lw))
    !$OMP TARGET UPDATE FROM(this%g_lw) IF(allocated(this%g_lw))

    !$ACC UPDATE HOST(this%mixing_ratio) IF(allocated(this%mixing_ratio)) ASYNC(1)
    !$ACC UPDATE HOST(this%od_sw) IF(allocated(this%od_sw)) ASYNC(1)
    !$ACC UPDATE HOST(this%ssa_sw) IF(allocated(this%ssa_sw)) ASYNC(1)
    !$ACC UPDATE HOST(this%g_sw) IF(allocated(this%g_sw)) ASYNC(1)
    !$ACC UPDATE HOST(this%od_lw) IF(allocated(this%od_lw)) ASYNC(1)
    !$ACC UPDATE HOST(this%ssa_lw) IF(allocated(this%ssa_lw)) ASYNC(1)
    !$ACC UPDATE HOST(this%g_lw) IF(allocated(this%g_lw)) ASYNC(1)

  end subroutine update_host_aerosol

  !---------------------------------------------------------------------
  ! updates fields on device
  subroutine update_device_aerosol(this)

    type(aerosol_type), intent(inout) :: this
#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET UPDATE TO(this%mixing_ratio) IF(allocated(this%mixing_ratio))
    !$OMP TARGET UPDATE TO(this%od_sw) IF(allocated(this%od_sw))
    !$OMP TARGET UPDATE TO(this%ssa_sw) IF(allocated(this%ssa_sw))
    !$OMP TARGET UPDATE TO(this%g_sw) IF(allocated(this%g_sw))
    !$OMP TARGET UPDATE TO(this%od_lw) IF(allocated(this%od_lw))
    !$OMP TARGET UPDATE TO(this%ssa_lw) IF(allocated(this%ssa_lw))
    !$OMP TARGET UPDATE TO(this%g_lw) IF(allocated(this%g_lw))

    !$ACC UPDATE DEVICE(this%mixing_ratio) IF(allocated(this%mixing_ratio)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%od_sw) IF(allocated(this%od_sw)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%ssa_sw) IF(allocated(this%ssa_sw)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%g_sw) IF(allocated(this%g_sw)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%od_lw) IF(allocated(this%od_lw)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%ssa_lw) IF(allocated(this%ssa_lw)) ASYNC(1)
    !$ACC UPDATE DEVICE(this%g_lw) IF(allocated(this%g_lw)) ASYNC(1)

  end subroutine update_device_aerosol

  !---------------------------------------------------------------------
  ! Deletes fields on device
  subroutine delete_device_aerosol(this)

    type(aerosol_type), intent(inout) :: this
#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET EXIT DATA MAP(DELETE:this%mixing_ratio) IF(allocated(this%mixing_ratio))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%od_sw) IF(allocated(this%od_sw))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%ssa_sw) IF(allocated(this%ssa_sw))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%g_sw) IF(allocated(this%g_sw))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%od_lw) IF(allocated(this%od_lw))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%ssa_lw) IF(allocated(this%ssa_lw))
    !$OMP TARGET EXIT DATA MAP(DELETE:this%g_lw) IF(allocated(this%g_lw))

    !$ACC EXIT DATA DELETE(this%mixing_ratio) IF(allocated(this%mixing_ratio)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%od_sw) IF(allocated(this%od_sw)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%ssa_sw) IF(allocated(this%ssa_sw)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%g_sw) IF(allocated(this%g_sw)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%od_lw) IF(allocated(this%od_lw)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%ssa_lw) IF(allocated(this%ssa_lw)) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%g_lw) IF(allocated(this%g_lw)) ASYNC(1)

  end subroutine delete_device_aerosol
#endif

end module radiation_aerosol
