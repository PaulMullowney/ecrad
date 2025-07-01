! radiation_mcica_omp_lw.F90 - Monte-Carlo Independent Column Approximation longtwave solver
!
! (C) Copyright 2015- ECMWF.
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
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-07-12  R. Hogan  Call fast adding method if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

#include "ecrad_config.h"

module radiation_mcica_omp_lw

  use radiation_io, only : nulout
  public

contains

  !---------------------------------------------------------------------
  ! Longwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_omp_lw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, &
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, debug_dump_lw
    use radiation_two_stream, only     : calc_no_scattering_transmittance_lw, calc_ref_trans_acc_lw, calc_ref_trans_lw_OMP_single_level, &
         &                               calc_no_scattering_transmittance_lw_OMP, calc_no_scattering_transmittance_lw_OMP_single_cell
    use radiation_adding_ica_lw, only  : adding_ica_lw, fast_adding_ica_lw, fast_adding_ica_lw_OMP, &
         &                               calc_fluxes_no_scattering_lw, calc_fluxes_no_scattering_lw_OMP
    
#if defined(OMPGPU)
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica, &
              & calc_lw_derivatives_ica_omp, modify_lw_derivatives_ica_omp, &
              & calc_lw_derivatives_ica_omp_gm, modify_lw_derivatives_ica_omp_gm
#else
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
#endif
    use radiation_cloud_generator_acc, only: cloud_generator_acc, cloud_generator_OMP
    use radiation_cloud_cover, only    : beta2alpha, MaxCloudFrac

#ifdef HAVE_ROCTX
    use roctx_profiling, only: roctxmarka, roctxrangepusha, roctxrangepop
    use iso_c_binding, only: c_null_char
#endif

    implicit none

#ifdef HAVE_ROCTX
    integer(4) :: roctx_ret
#endif

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables : Mapped into Global Memory
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up, flux_dn
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up_clear, flux_dn_clear

    ! Identify clear-sky layers
    logical :: is_clear_sky_layer(nlev, istartcol:iendcol)

    ! workaround that allows inling of cloud generator
    real(jprb), dimension(config%pdf_sampler%ncdf, config%pdf_sampler%nfsd)  :: sample_val

    ! temporary arrays to increase performance
    real(jprb), dimension(nlev, istartcol:iendcol) :: frac, frac_std
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: overlap_param

    ! temporary arrays
    real(jprb), dimension(nlev, istartcol:iendcol) :: cum_cloud_cover
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: pair_cloud_cover

    ! Cumulative product needed in computation of total_cloud_cover
    real(jprb) :: cum_product(istartcol:iendcol)

    ! First and last cloudy layers
    integer :: ibegin(istartcol:iendcol), iend(istartcol:iendcol)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables : Was stack and is now in Global Memory
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: trans_clear, reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: source_up, source_dn

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_lw,nlev, istartcol:iendcol) :: od_scaling
    
        ! Temporary working array
    real(jprb), dimension(config%n_g_lw,nlev+1, istartcol:iendcol) :: tmp_work_source
    real(jprb), dimension(config%n_g_lw, istartcol:iendcol) :: tmp_derivatives

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables : Stack
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb) :: od_cloud_new
    
    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb) :: od_total, ssa_total, g_total

    ! Combined scattering optical depth
    real(jprb) :: scat_od, scat_od_total(config%n_g_lw)

    ! Total cloud cover output from the cloud generator
    real(jprb) :: total_cloud_cover

    ! "Alpha" overlap parameter
    real(jprb) :: overlap_alpha

    ! Temporary storage for more efficient summation
    real(jprb) :: sum_up, sum_dn, sum_up_clr, sum_dn_clr

    ! Index of the highest cloudy layer
    integer :: i_cloud_top

    ! Number of g points
    integer :: ng

    ! Loop indices for level, column and g point
    integer :: jlev, jcol, jg

    real(jphook) :: hook_handle
    real(jprb)  totalMem
    integer :: file_idx, fidx1, fidx2, fidx3
#ifdef DEBUG_CORRECTNESS_RADIATION1
    integer :: s1, s2, s3
#endif

#ifdef HAVE_ROCTX
    roctx_ret = roctxRangePushA("radiation::mcica_omp_lw"//c_null_char)
#endif

    if (lhook) call dr_hook('radiation_mcica_omp_lw:solver_mcica_omp_lw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: longwave McICA OMP requires clear-sky calculation to be performed'
      call radiation_abort()
    end if

    ng = config%n_g_lw

    totalMem = 6*ng * nlev
    totalMem = totalMem+5*(ng)*(nlev+1) !flux_up,dn,up_clear,dn_clear,source
    totalmem = totalMem*(iendcol-istartcol)*SIZEOF((real(jprb)))/1.e9
    write(nulout,'(a,a,i0,a,g0.5)') __FILE__, " : LINE = ", __LINE__, " total_memory=",totalMem

    !$OMP TARGET ENTER DATA MAP(ALLOC: flux_up, flux_dn, flux_up_clear, flux_dn_clear, &
    !$OMP             is_clear_sky_layer, &
    !$OMP             sample_val, frac, frac_std, overlap_param, cum_cloud_cover, &
    !$OMP             pair_cloud_cover, cum_product, ibegin, iend)
    !$OMP TARGET DATA MAP(PRESENT, ALLOC: config, single_level, cloud, od, ssa, g, od_cloud, ssa_cloud, &
    !$OMP             g_cloud, planck_hl, emission, albedo, flux)

    !
    ! These allocations need to moved into a TEAMS PRIVATE clause, I think.
    !

    !$OMP TARGET ENTER DATA MAP(ALLOC: trans_clear, od_scaling, &
    !$OMP   reflectance, transmittance, source_up, source_dn, tmp_work_source, &
    !$OMP   tmp_derivatives)


#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
    do jlev = 1,config%pdf_sampler%nfsd
      do jcol = 1,config%pdf_sampler%ncdf
        sample_val(jcol,jlev) = config%pdf_sampler%val(jcol,jlev)
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev
        frac(jlev, jcol) = cloud%fraction(jcol,jlev)
        frac_std(jlev, jcol) = cloud%fractional_std(jcol,jlev)
        is_clear_sky_layer(jlev,jcol) = .true.
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev-1
        overlap_param(jlev, jcol) = cloud%overlap_param(jcol,jlev)
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(overlap_alpha)
    do jcol = istartcol,iendcol
      !---------------------------------------------------------------------
      ! manual inline from cum_cloud_cover_exp_ran >>>>>>>>>>>>>>>>>>>>>>>>
      ! Loop to compute total cloud cover and the cumulative cloud cover
      ! down to the base of each layer
      do jlev = 1,nlev-1
        ! Convert to "alpha" overlap parameter if necessary
        if (config%use_beta_overlap) then
          overlap_alpha = beta2alpha(overlap_param(jlev,jcol), &
                &                     frac(jlev,jcol), frac(jlev+1,jcol))
        else
          overlap_alpha = overlap_param(jlev,jcol)
        end if

        ! Compute the combined cloud cover of layers jlev and jlev+1
        pair_cloud_cover(jlev, jcol) = overlap_alpha*max(frac(jlev,jcol),frac(jlev+1,jcol)) &
              &  + (1.0_jprb - overlap_alpha) &
              &  * (frac(jlev,jcol)+frac(jlev+1,jcol)-frac(jlev,jcol)*frac(jlev+1,jcol))
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

#ifdef DEBUG
    write(nulout,'(a,a,a,i0)') "    ", __FILE__, " : LINE = ", __LINE__
#endif

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    do jcol = istartcol,iendcol
      cum_cloud_cover(1, jcol) = frac(1,jcol)
      cum_product(jcol) = 1.0_jprb - frac(1,jcol)
      do jlev = 1,nlev-1
        if (frac(jlev,jcol) >= MaxCloudFrac) then
          ! Cloud cover has reached one
          cum_product(jcol) = 0.0_jprb
        else
          cum_product(jcol) = cum_product(jcol) * (1.0_jprb - pair_cloud_cover(jlev, jcol)) &
                &  / (1.0_jprb - frac(jlev,jcol))
        end if
        cum_cloud_cover(jlev+1, jcol) = 1.0_jprb - cum_product(jcol)
      end do
      flux%cloud_cover_lw(jcol) = cum_cloud_cover(nlev,jcol);
      if (flux%cloud_cover_lw(jcol) < config%cloud_fraction_threshold) then
        ! Treat column as clear sky: calling function therefore will not
        ! use od_scaling so we don't need to calculate it
        flux%cloud_cover_lw(jcol) = 0.0_jprb
      end if
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    do jcol = istartcol,iendcol
      if (flux%cloud_cover_lw(jcol) >= config%cloud_fraction_threshold) then
        ! Cloud is present: need to calculate od_scaling
        ! Find range of cloudy layers
        ibegin(jcol) = nlev
        do jlev = 1, nlev
          if( frac(jlev,jcol) > 0.0_jprb ) then
            ibegin(jcol) = min(jlev, ibegin(jcol))
          end if
        end do

        iend(jcol) = ibegin(jcol)
        do jlev = ibegin(jcol)+1,nlev
          if (frac(jlev,jcol) > 0.0_jprb) then
            iend(jcol) = max(jlev, iend(jcol))
          end if
        end do
      end if
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev+1
         do jg = 1, ng
            flux_up_clear(jg,jlev,jcol) = 0.0_jprb
            flux_dn_clear(jg,jlev,jcol) = 0.0_jprb
            flux_up(jg,jlev,jcol) = 0.0_jprb
            flux_dn(jg,jlev,jcol) = 0.0_jprb
         end do
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !
    ! This kernel does band independent computations. Some computation is done across veritical levels
    !
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) THREAD_LIMIT(1024)
    do jcol = istartcol,iendcol
       do jg = 1, ng
          call cloud_generator_OMP(jg, ng, nlev, &
               &  single_level%iseed(jcol) + 997, &
               &  config%cloud_fraction_threshold, &
               &  frac(:,jcol), overlap_param(:,jcol), &
               &  config%cloud_inhom_decorr_scaling, frac_std(:,jcol), &
               &  config%pdf_sampler%ncdf, config%pdf_sampler%nfsd, &
               &  config%pdf_sampler%fsd1, config%pdf_sampler%inv_fsd_interval, &
               &  sample_val, &
               &  od_scaling(:,:,jcol), flux%cloud_cover_lw(jcol)+0.0_jprb, & ! Workaround for nvhpc-24.1
               &  ibegin(jcol), iend(jcol), &
               &  cum_cloud_cover=cum_cloud_cover(:,jcol), &
               &  pair_cloud_cover=pair_cloud_cover(:,jcol))
       enddo
    enddo

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(od_cloud_new, od_total, ssa_total, g_total) THREAD_LIMIT(1024)
    do jcol = istartcol,iendcol
       do jg = 1, ng
          ! Clear-sky calculation
          ! Non-scattering case: use simpler functions for
          ! transmission and emission
          call calc_no_scattering_transmittance_lw_OMP(jg, ng, nlev, od(:,:,jcol), &
               &  planck_hl(:,:,jcol), trans_clear(:,:,jcol), source_up(:,:,jcol), &
               &  source_dn(:,:,jcol))

          ! Simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw_OMP(jg, ng, nlev, &
               &  trans_clear(:,:,jcol), source_up(:,:,jcol), source_dn(:,:,jcol), &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up_clear(:,:,jcol), flux_dn_clear(:,:,jcol))

          ! Store surface spectral downwelling fluxes
          flux%lw_dn_surf_clear_g(jg,jcol) = flux_dn_clear(jg,nlev+1,jcol)
          ! Do cloudy-sky calculation; add a prime number to the seed in
          ! the longwave

          if (flux%cloud_cover_lw(jcol) >= config%cloud_fraction_threshold) then
             ! Total-sky calculation
             i_cloud_top = nlev+1

             do jlev = 1,nlev
                ! Compute combined gas+aerosol+cloud optical properties
                if (frac(jlev,jcol) >= config%cloud_fraction_threshold) then
                   is_clear_sky_layer(jlev,jcol) = .false.
                   ! Get index to the first cloudy layer from the top
                   if (i_cloud_top > jlev) then
                      i_cloud_top = jlev
                   end if
                   
                   od_cloud_new = od_scaling(jg,jlev,jcol) &
                        &  * od_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol)
                   od_total  = od(jg,jlev,jcol) + od_cloud_new
                   ssa_total = 0.0_jprb
                   g_total   = 0.0_jprb

                   if (config%do_lw_cloud_scattering) then
                      ! Scattering case: calculate reflectance and
                      ! transmittance at each model level
                      if (od_total > 0.0_jprb) then
                         scat_od = ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                              &     * od_cloud_new
                         ssa_total = scat_od / od_total
                         if (scat_od > 0.0_jprb) then
                            g_total = g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                                 &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                                 &     *  od_cloud_new / scat_od
                         end if
                      end if

                      ! Compute cloudy-sky reflectance, transmittance etc at
                      ! each model level
                      call calc_ref_trans_lw_OMP_single_level(jg, ng, &
                           &  od_total, ssa_total, g_total, &
                           &  planck_hl(:,jlev:jlev+1,jcol), &
                           &  reflectance(:,jlev,jcol), transmittance(:,jlev,jcol), &
                           &  source_up(:,jlev,jcol), source_dn(:,jlev,jcol))
                   else
                      ! No-scattering case: use simpler functions for
                      ! transmission and emission
                      call calc_no_scattering_transmittance_lw_OMP_single_cell(od_total, &
                           &  planck_hl(jg,jlev,jcol), planck_hl(jg,jlev+1,jcol), transmittance(jg,jlev,jcol), &
                           &  source_up(jg,jlev,jcol), source_dn(jg,jlev,jcol))
                   end if

                else
                   ! Clear-sky layer: copy over clear-sky values
                   reflectance(jg,jlev,jcol) = 0.0
                   transmittance(jg,jlev,jcol) = trans_clear(jg,jlev,jcol)
                end if
             end do

             if(config%do_lw_cloud_scattering) then
                ! Use adding method to compute fluxes but optimize for the
                ! presence of clear-sky layers
                call fast_adding_ica_lw_OMP(jg, ng, nlev, reflectance(:,:,jcol), transmittance(:,:,jcol), &
                     &  source_up(:,:,jcol), source_dn(:,:,jcol), &
                     &  emission(:,jcol), albedo(:,jcol), &
                     &  is_clear_sky_layer(:,jcol), i_cloud_top, flux_dn_clear(:,:,jcol), &
                     &  flux_up(:,:,jcol), flux_dn(:,:,jcol), &
                     &  source=tmp_work_source(:,:,jcol))
             else
                ! Simpler down-then-up method to compute fluxes
                call calc_fluxes_no_scattering_lw_OMP(jg, ng, nlev, &
                     &  transmittance(:,:,jcol), source_up(:,:,jcol), source_dn(:,:,jcol), &
                     &  emission(:,jcol), albedo(:,jcol), &
                     &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
             end if

             ! Cloudy flux profiles currently assume completely overcast
             ! skies; perform weighted average with clear-sky profile
             ! Store surface spectral downwelling fluxes
             flux%lw_dn_surf_g(jg,jcol) = flux%cloud_cover_lw(jcol)*flux_dn(jg,nlev+1,jcol) &
                  &  + (1.0_jprb - flux%cloud_cover_lw(jcol))*flux%lw_dn_surf_clear_g(jg,jcol)
          else
             ! No cloud in profile and clear-sky fluxes already
             ! calculated: copy them over
             flux%lw_dn_surf_g(jg,jcol) = flux%lw_dn_surf_clear_g(jg,jcol)
          end if
       end do
       !!$OMP END PARALLEL DO
    end do
    !!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !
    ! Separate out the derivative computation, which has reductions across spectral bands
    !
    if (config%do_lw_derivatives) then
#ifdef WORKAROUND_NESTED_PARALLEL_LW_DERIVATIVES
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
       !$OMP TARGET TEAMS DISTRIBUTE
#endif
       do jcol = istartcol,iendcol
          if (flux%cloud_cover_lw(jcol) >= config%cloud_fraction_threshold) then
             call calc_lw_derivatives_ica_omp_gm(ng, nlev, jcol, transmittance(:,:,jcol), flux_up(:,nlev+1,jcol), &
                  &                       flux%lw_derivatives, tmp_derivatives(:,jcol))
             if (flux%cloud_cover_lw(jcol) < 1.0_jprb - config%cloud_fraction_threshold) then
                ! Modify the existing derivative with the contribution from the clear sky
                call modify_lw_derivatives_ica_omp_gm(ng, nlev, jcol, trans_clear(:,:,jcol), flux_up_clear(:,nlev+1,jcol), &
                     &                             1.0_jprb-flux%cloud_cover_lw(jcol), flux%lw_derivatives, tmp_derivatives(:,jcol))
             end if
          else
             call calc_lw_derivatives_ica_omp_gm(ng, nlev, jcol, trans_clear(:,:,jcol), flux_up_clear(:,nlev+1,jcol), &
                  &                       flux%lw_derivatives, tmp_derivatives(:,jcol))
          endif
       enddo
#ifdef WORKAROUND_NESTED_PARALLEL_LW_DERIVATIVES
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
       !$OMP END TARGET TEAMS DISTRIBUTE
#endif
    endif

#ifdef DEBUG_WARNING
    write(nulout,'(a,a,a,i0,a)') "    ", __FILE__, " : LINE = ", __LINE__, " OMP Parallel reductions aren't working right. This needs more attention."
#endif

    ! Loop through columns
#ifdef WORKAROUND_NESTED_PARALLEL_MCICA_OMP_LW
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(total_cloud_cover, sum_up, sum_dn, sum_up_clr, sum_dn_clr)
#else
    !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) PRIVATE(total_cloud_cover, sum_up, sum_dn, sum_up_clr, sum_dn_clr) THREAD_LIMIT(128)
#endif
    do jcol = istartcol,iendcol
       do jlev = 1,nlev+1

        sum_up_clr = 0._jprb
        sum_dn_clr = 0._jprb
#ifndef WORKAROUND_NESTED_PARALLEL_MCICA_OMP_LW
        !$OMP PARALLEL DO REDUCTION(+:sum_up_clr, sum_dn_clr)
#endif
        do jg = 1,ng
          sum_up_clr = sum_up_clr + flux_up_clear(jg,jlev,jcol)
          sum_dn_clr = sum_dn_clr + flux_dn_clear(jg,jlev,jcol)
        end do
#ifndef WORKAROUND_NESTED_PARALLEL_MCICA_OMP_LW
        !$OMP END PARALLEL DO
#endif
        flux%lw_up_clear(jcol,jlev) = sum_up_clr
        flux%lw_dn_clear(jcol,jlev) = sum_dn_clr

        total_cloud_cover = flux%cloud_cover_lw(jcol)
        if (total_cloud_cover >= config%cloud_fraction_threshold) then

          ! Store overcast broadband fluxes
          sum_up = 0._jprb
          sum_dn = 0._jprb
#ifndef WORKAROUND_NESTED_PARALLEL_MCICA_OMP_LW
          !$OMP PARALLEL DO REDUCTION(+:sum_up_clr, sum_dn_clr)
#endif
          do jg = 1,ng
            sum_up = sum_up + flux_up(jg,jlev,jcol)
            sum_dn = sum_dn + flux_dn(jg,jlev,jcol)
          end do
#ifndef WORKAROUND_NESTED_PARALLEL_MCICA_OMP_LW
          !$OMP END PARALLEL DO
#endif
          flux%lw_up(jcol,jlev) = total_cloud_cover*sum_up + (1.0_jprb - total_cloud_cover)*sum_up_clr
          flux%lw_dn(jcol,jlev) = total_cloud_cover*sum_dn + (1.0_jprb - total_cloud_cover)*sum_dn_clr

        else

          flux%lw_up(jcol,jlev) = sum_up_clr
          flux%lw_dn(jcol,jlev) = sum_dn_clr

        end if
      end do
    end do
#ifdef WORKAROUND_NESTED_PARALLEL_MCICA_OMP_LW
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
    !$OMP END TARGET TEAMS DISTRIBUTE
#endif

#ifdef DEBUG_CORRECTNESS_RADIATION1
    write(nulout,'(a)') "*******************************************************************"
    write(nulout,'(a,a,a,i0,a,g0.5)') "Correctness Check : ", __FILE__, " : LINE = ", __LINE__, " config%cloud_fraction_threshold=", config%cloud_fraction_threshold
    !$OMP TARGET UPDATE FROM(flux_up_clear,flux_dn_clear,flux_up,flux_dn,emission,albedo,planck_hl,od,flux%lw_up,flux%lw_dn,flux%cloud_cover_lw,flux%lw_up_clear,flux%lw_dn_clear,flux%lw_up,flux%lw_dn,source_up_clear,source_dn_clear,trans_clear,od_scaling,pair_cloud_cover)
    s1 = lbound(flux%cloud_cover_lw,1)
    write(nulout,'(a,g0.5,a,i0)') "flux%cloud_cover_lw=", flux%cloud_cover_lw(s1), " at indices ", s1
    s1 = iendcol/2
    write(nulout,'(a,g0.5,a,i0)') "flux%cloud_cover_lw=", flux%cloud_cover_lw(s1), " at indices ", s1
    s1 = iendcol
    write(nulout,'(a,g0.5,a,i0)') "flux%cloud_cover_lw=", flux%cloud_cover_lw(s1), " at indices ", s1
    
    s1 = lbound(flux%lw_up,1)
    s2 = lbound(flux%lw_up,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up=", flux%lw_up(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol
    s2 = lbound(flux%lw_up,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up=", flux%lw_up(s1,s2), " at indices ", s1, " ", s2
    s1 = lbound(flux%lw_up,1)
    s2 = ubound(flux%lw_up,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up=", flux%lw_up(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol/2
    s2 = size(flux%lw_up,2)/2
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up=", flux%lw_up(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol
    s2 = ubound(flux%lw_up,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up=", flux%lw_up(s1,s2), " at indices ", s1, " ", s2

    s1 = lbound(flux%lw_dn,1)
    s2 = lbound(flux%lw_dn,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn=", flux%lw_dn(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol
    s2 = lbound(flux%lw_dn,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn=", flux%lw_dn(s1,s2), " at indices ", s1, " ", s2
    s1 = lbound(flux%lw_dn,1)
    s2 = ubound(flux%lw_dn,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn=", flux%lw_dn(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol/2
    s2 = size(flux%lw_dn,2)/2
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn=", flux%lw_dn(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol
    s2 = ubound(flux%lw_dn,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn=", flux%lw_dn(s1,s2), " at indices ", s1, " ", s2
    
    s1 = lbound(flux%lw_up_clear,1)
    s2 = lbound(flux%lw_up_clear,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up_clear=", flux%lw_up_clear(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol
    s2 = lbound(flux%lw_up_clear,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up_clear=", flux%lw_up_clear(s1,s2), " at indices ", s1, " ", s2
    s1 = lbound(flux%lw_up_clear,1)
    s2 = ubound(flux%lw_up_clear,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up_clear=", flux%lw_up_clear(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol/2
    s2 = size(flux%lw_up_clear,2)/2
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up_clear=", flux%lw_up_clear(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol
    s2 = ubound(flux%lw_up_clear,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_up_clear=", flux%lw_up_clear(s1,s2), " at indices ", s1, " ", s2
    
    s1 = lbound(flux%lw_dn_clear,1)
    s2 = lbound(flux%lw_dn_clear,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn_clear=", flux%lw_dn_clear(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol
    s2 = lbound(flux%lw_dn_clear,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn_clear=", flux%lw_dn_clear(s1,s2), " at indices ", s1, " ", s2
    s1 = lbound(flux%lw_dn_clear,1)
    s2 = ubound(flux%lw_dn_clear,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn_clear=", flux%lw_dn_clear(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol/2
    s2 = size(flux%lw_dn_clear,2)/2
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn_clear=", flux%lw_dn_clear(s1,s2), " at indices ", s1, " ", s2
    s1 = iendcol
    s2 = ubound(flux%lw_dn_clear,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "flux%lw_dn_clear=", flux%lw_dn_clear(s1,s2), " at indices ", s1, " ", s2

    s1 = lbound(flux_up,1)
    s2 = lbound(flux_up,2)
    s3 = lbound(flux_up,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up=", flux_up(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(flux_up,1)
    s2 = ubound(flux_up,2)
    s3 = lbound(flux_up,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up=", flux_up(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(flux_up,1)
    s2 = lbound(flux_up,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up=", flux_up(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(flux_up,1)/2
    s2 = size(flux_up,2)/2
    s3 = iendcol/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up=", flux_up(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(flux_up,1)
    s2 = ubound(flux_up,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up=", flux_up(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3

    s1 = lbound(flux_dn,1)
    s2 = lbound(flux_dn,2)
    s3 = lbound(flux_dn,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn=", flux_dn(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(flux_dn,1)
    s2 = ubound(flux_dn,2)
    s3 = lbound(flux_dn,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn=", flux_dn(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(flux_dn,1)
    s2 = lbound(flux_dn,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn=", flux_dn(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(flux_dn,1)/2
    s2 = size(flux_dn,2)/2
    s3 = iendcol/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn=", flux_dn(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(flux_dn,1)
    s2 = ubound(flux_dn,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn=", flux_dn(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3

    s1 = lbound(flux_up_clear,1)
    s2 = lbound(flux_up_clear,2)
    s3 = lbound(flux_up_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up_clear=", flux_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(flux_up_clear,1)
    s2 = ubound(flux_up_clear,2)
    s3 = lbound(flux_up_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up_clear=", flux_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(flux_up_clear,1)
    s2 = lbound(flux_up_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up_clear=", flux_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(flux_up_clear,1)/2
    s2 = size(flux_up_clear,2)/2
    s3 = iendcol/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up_clear=", flux_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(flux_up_clear,1)
    s2 = ubound(flux_up_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_up_clear=", flux_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3

    s1 = lbound(flux_dn_clear,1)
    s2 = lbound(flux_dn_clear,2)
    s3 = lbound(flux_dn_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn_clear=", flux_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(flux_dn_clear,1)
    s2 = ubound(flux_dn_clear,2)
    s3 = lbound(flux_dn_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn_clear=", flux_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(flux_dn_clear,1)
    s2 = lbound(flux_dn_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn_clear=", flux_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(flux_dn_clear,1)/2
    s2 = size(flux_dn_clear,2)/2
    s3 = iendcol/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn_clear=", flux_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(flux_dn_clear,1)
    s2 = ubound(flux_dn_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "flux_dn_clear=", flux_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3

    s1 = lbound(trans_clear,1)
    s2 = lbound(trans_clear,2)
    s3 = lbound(trans_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "trans_clear=", trans_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(trans_clear,1)
    s2 = ubound(trans_clear,2)
    s3 = lbound(trans_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "trans_clear=", trans_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(trans_clear,1)
    s2 = lbound(trans_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "trans_clear=", trans_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(trans_clear,1)/2
    s2 = size(trans_clear,2)/2
    s3 = iendcol/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "trans_clear=", trans_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(trans_clear,1)
    s2 = ubound(trans_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "trans_clear=", trans_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3

    s1 = lbound(source_up_clear,1)
    s2 = lbound(source_up_clear,2)
    s3 = lbound(source_up_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_up_clear=", source_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(source_up_clear,1)
    s2 = ubound(source_up_clear,2)
    s3 = lbound(source_up_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_up_clear=", source_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(source_up_clear,1)
    s2 = lbound(source_up_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_up_clear=", source_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(source_up_clear,1)/2
    s2 = size(source_up_clear,2)/2
    s3 = iendcol/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_up_clear=", source_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(source_up_clear,1)
    s2 = ubound(source_up_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_up_clear=", source_up_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3

    s1 = lbound(source_dn_clear,1)
    s2 = lbound(source_dn_clear,2)
    s3 = lbound(source_dn_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_dn_clear=", source_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(source_dn_clear,1)
    s2 = ubound(source_dn_clear,2)
    s3 = lbound(source_dn_clear,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_dn_clear=", source_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(source_dn_clear,1)
    s2 = lbound(source_dn_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_dn_clear=", source_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(source_dn_clear,1)/2
    s2 = size(source_dn_clear,2)/2
    s3 = iendcol/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_dn_clear=", source_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(source_dn_clear,1)
    s2 = ubound(source_dn_clear,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "source_dn_clear=", source_dn_clear(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3

    s1 = lbound(od_scaling,1)
    s2 = lbound(od_scaling,2)
    s3 = lbound(od_scaling,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od_scaling=", od_scaling(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(od_scaling,1)
    s2 = ubound(od_scaling,2)
    s3 = lbound(od_scaling,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od_scaling=", od_scaling(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(od_scaling,1)
    s2 = lbound(od_scaling,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od_scaling=", od_scaling(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(od_scaling,1)/2
    s2 = size(od_scaling,2)/2
    s3 = iendcol/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od_scaling=", od_scaling(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(od_scaling,1)
    s2 = ubound(od_scaling,2)
    s3 = iendcol
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od_scaling=", od_scaling(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3

    s1 = lbound(pair_cloud_cover,1)
    s2 = lbound(pair_cloud_cover,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "pair_cloud_cover=", pair_cloud_cover(s1,s2), " at indices ", s1, " ", s2
    s1 = ubound(pair_cloud_cover,1)
    s2 = lbound(pair_cloud_cover,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "pair_cloud_cover=", pair_cloud_cover(s1,s2), " at indices ", s1, " ", s2
    s1 = lbound(pair_cloud_cover,1)
    s2 = ubound(pair_cloud_cover,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "pair_cloud_cover=", pair_cloud_cover(s1,s2), " at indices ", s1, " ", s2
    s1 = size(pair_cloud_cover,1)/2
    s2 = size(pair_cloud_cover,2)/2
    write(nulout,'(a,g0.5,a,i0,a,i0)') "pair_cloud_cover=", pair_cloud_cover(s1,s2), " at indices ", s1, " ", s2
    s1 = ubound(pair_cloud_cover,1)
    s2 = ubound(pair_cloud_cover,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "pair_cloud_cover=", pair_cloud_cover(s1,s2), " at indices ", s1, " ", s2

    s1 = lbound(emission,1)
    s2 = lbound(emission,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "emission=", emission(s1,s2), " at indices ", s1, " ", s2
    s1 = ubound(emission,1)
    s2 = lbound(emission,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "emission=", emission(s1,s2), " at indices ", s1, " ", s2
    s1 = lbound(emission,1)
    s2 = ubound(emission,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "emission=", emission(s1,s2), " at indices ", s1, " ", s2
    s1 = size(emission,1)/2
    s2 = size(emission,2)/2
    write(nulout,'(a,g0.5,a,i0,a,i0)') "emission=", emission(s1,s2), " at indices ", s1, " ", s2
    s1 = ubound(emission,1)
    s2 = ubound(emission,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "emission=", emission(s1,s2), " at indices ", s1, " ", s2
    
    s1 = lbound(albedo,1)
    s2 = lbound(albedo,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "albedo=", albedo(s1,s2), " at indices ", s1, " ", s2
    s1 = ubound(albedo,1)
    s2 = lbound(albedo,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "albedo=", albedo(s1,s2), " at indices ", s1, " ", s2
    s1 = lbound(albedo,1)
    s2 = ubound(albedo,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "albedo=", albedo(s1,s2), " at indices ", s1, " ", s2
    s1 = size(albedo,1)/2
    s2 = size(albedo,2)/2
    write(nulout,'(a,g0.5,a,i0,a,i0)') "albedo=", albedo(s1,s2), " at indices ", s1, " ", s2
    s1 = ubound(albedo,1)
    s2 = ubound(albedo,2)
    write(nulout,'(a,g0.5,a,i0,a,i0)') "albedo=", albedo(s1,s2), " at indices ", s1, " ", s2
    

    s1 = lbound(planck_hl,1)
    s2 = lbound(planck_hl,2)
    s3 = lbound(planck_hl,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "planck_hl=", planck_hl(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(planck_hl,1)
    s2 = lbound(planck_hl,2)+1
    s3 = lbound(planck_hl,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "planck_hl=", planck_hl(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(planck_hl,1)
    s2 = ubound(planck_hl,2)-1
    s3 = lbound(planck_hl,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "planck_hl=", planck_hl(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(planck_hl,1)
    s2 = ubound(planck_hl,2)
    s3 = lbound(planck_hl,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "planck_hl=", planck_hl(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(planck_hl,1)/2
    s2 = size(planck_hl,2)/2
    s3 = size(planck_hl,3)/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "planck_hl=", planck_hl(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(planck_hl,1)
    s2 = ubound(planck_hl,2)
    s3 = ubound(planck_hl,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "planck_hl=", planck_hl(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(od,1)
    s2 = lbound(od,2)
    s3 = lbound(od,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od=", od(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = lbound(od,1)
    s2 = lbound(od,2)+1
    s3 = lbound(od,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od=", od(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(od,1)/2
    s2 = size(od,2)/2
    s3 = size(od,3)/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od=", od(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = ubound(od,1)
    s2 = ubound(od,2)
    s3 = ubound(od,3)
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "od=", od(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    write(nulout,'(a)') "*******************************************************************"
#endif
    !$OMP TARGET EXIT DATA MAP(DELETE: trans_clear, od_scaling, &
    !$OMP   reflectance, transmittance, source_up, source_dn, tmp_work_source, &
    !$OMP   tmp_derivatives)

    !$OMP TARGET EXIT DATA MAP(DELETE: flux_up, flux_dn, flux_up_clear, flux_dn_clear, &
    !$OMP             is_clear_sky_layer, &
    !$OMP             sample_val, frac, frac_std, overlap_param, cum_cloud_cover, &
    !$OMP             pair_cloud_cover, cum_product, ibegin, iend)
    !$OMP END TARGET DATA

#ifdef HAVE_ROCTX
    call roctxRangePop()
    call roctxMarkA("radiation::radiation::mcica_omp_lw"//c_null_char)
#endif
    if (lhook) call dr_hook('radiation_mcica_omp_lw:solver_mcica_omp_lw',1,hook_handle)

  end subroutine solver_mcica_omp_lw

end module radiation_mcica_omp_lw
