SUBROUTINE SRTM_TAUMOL25 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1,&
 & P_COLH2O  , P_COLMOL , P_COLO3,&
 & K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0   &
 & , laytrop_min, laytrop_max)  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 25:  16000-22650 cm-1 (low - H2O; high - nothing)

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG25
USE YOESRTA25, ONLY : ABSA, SFLUXREFC, ABSO3AC, ABSO3BC, RAYLC, LAYREFFR  
USE YOESRTWN , ONLY : NSPA
use radiation_io, only : nulout

IMPLICIT NONE

!-- Output
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLH2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLO3(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 

REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_SFLUXZEN(KIDIA:KFDIA,JPG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUG(KIDIA:KFDIA,KLEV,JPG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUR(KIDIA:KFDIA,KLEV,JPG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, IND0, IND1, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON
INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max
INTEGER(KIND=JPIM) :: I_LAY_NEXT

REAL(KIND=JPRB) ::  &
 & Z_TAURAY  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#ifdef DEBUG_CORRECTNESS_IFSRRTM
    integer :: s1, s2, s3
#endif

#ifdef DEBUG
    write(nulout,'(a,a,a,i0,a)') "    ", __FILE__, " : LINE = ", __LINE__, " Begin"
#endif
    !$ACC DATA CREATE(I_LAYSOLFR) &
    !$ACC     PRESENT(P_FAC00, P_FAC01, P_FAC10, P_FAC11, K_JP, K_JT, K_JT1, &
    !$ACC             P_COLH2O, P_COLMOL, P_COLO3, K_LAYTROP, P_SFLUXZEN, &
    !$ACC             P_TAUG, P_TAUR, PRMU0)
    !$OMP TARGET ENTER DATA MAP(ALLOC: I_LAYSOLFR)
    !$OMP TARGET DATA MAP(PRESENT, ALLOC: P_FAC00, P_FAC01, P_FAC10, P_FAC11, K_JP, K_JT, K_JT1, &
    !$OMP             P_COLH2O, P_COLMOL, P_COLO3, K_LAYTROP, P_SFLUXZEN, &
    !$OMP             P_TAUG, P_TAUR, PRMU0)
    
    !$OMP TARGET DATA MAP(PRESENT, ALLOC: laytrop_min, laytrop_max)
    
    if (.not. present(laytrop_min) .and. .not. present(laytrop_max)) then
#if defined(_OPENACC) || defined(OMPGPU)
    laytrop_min = HUGE(laytrop_min) 
    laytrop_max = -HUGE(laytrop_max)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO REDUCTION(min:laytrop_min) REDUCTION(max:laytrop_max)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR REDUCTION(min:laytrop_min) REDUCTION(max:laytrop_max)
    do iplon = KIDIA,KFDIA
      laytrop_min = MIN(laytrop_min, k_laytrop(iplon))
      laytrop_max = MAX(laytrop_max, k_laytrop(iplon))
    end do
    !$ACC END PARALLEL
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
    laytrop_min = MINVAL(k_laytrop(KIDIA:KFDIA))
    laytrop_max = MAXVAL(k_laytrop(KIDIA:KFDIA))
#endif
    endif

    i_nlayers = klev
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    !$ACC WAIT
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG(STATIC:1) VECTOR
    DO iplon = KIDIA, KFDIA
      i_laysolfr(iplon) = k_laytrop(iplon)
    ENDDO
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

#if defined(OMPGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    DO iplon = KIDIA, KFDIA
       DO i_lay = 1, laytrop_min
         IF (k_jp(iplon,i_lay) < layreffr .AND.   &
              &    k_jp(iplon,i_lay+1) >= layreffr) &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
       ENDDO
    ENDDO
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ind0, ind1, z_tauray)
    DO i_lay = 1, laytrop_min
      DO iplon = KIDIA, KFDIA
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(25) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(25) + 1
!$NEC unroll(NG25)
         DO ig = 1 , ng25
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absa(ind0,ig)   + &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig)  + &
                & p_fac01(iplon,i_lay) * absa(ind1,ig)    + &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) + &
                & p_colo3(iplon,i_lay) * abso3ac(ig)
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
      ENDDO
    ENDDO
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
    !$ACC LOOP SEQ
    DO i_lay = 1, laytrop_min
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ind0, ind1, z_tauray)
      !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1)
      DO iplon = KIDIA, KFDIA
         IF (k_jp(iplon,i_lay) < layreffr .AND.   &
              &    k_jp(iplon,i_lay+1) >= layreffr) &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(25) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(25) + 1
         !$ACC LOOP SEQ PRIVATE(z_tauray)
!$NEC unroll(NG25)
         DO ig = 1 , ng25
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absa(ind0,ig)   + &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig)  + &
                & p_fac01(iplon,i_lay) * absa(ind1,ig)    + &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) + &
                & p_colo3(iplon,i_lay) * abso3ac(ig)
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
      ENDDO
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    ENDDO
#endif

    !$ACC LOOP SEQ
    DO i_lay = laytrop_min+1, laytrop_max
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ind0, ind1, z_tauray)
      !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1)
      DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr .AND.   &
                 &    k_jp(iplon,i_lay+1) >= layreffr) &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(25) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(25) + 1
!$NEC unroll(NG25)
            !$ACC LOOP SEQ PRIVATE(z_tauray)
            DO ig = 1 , ng25
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absa(ind0,ig)   + &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig)  + &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig)    + &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) + &
                   & p_colo3(iplon,i_lay) * abso3ac(ig)
              IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
!$NEC unroll(NG25)
            !$ACC LOOP SEQ PRIVATE(z_tauray)
            DO ig = 1 , ng25
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * abso3bc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
      ENDDO
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    ENDDO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(z_tauray)
    !$ACC LOOP SEQ
    DO ig = 1 , ng25
      !$ACC LOOP SEQ
      DO i_lay = laytrop_max+1, i_nlayers
        !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(z_tauray)
        DO iplon = KIDIA, KFDIA
          z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
          p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * abso3bc(ig)
          p_taur(iplon,i_lay,ig) = z_tauray
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$ACC WAIT
    !$ACC END DATA

    !$OMP TARGET EXIT DATA MAP(DELETE: I_LAYSOLFR)
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA

#ifdef DEBUG_CORRECTNESS_IFSRRTM
    write(nulout,'(a)') "*******************************************************************"
    write(nulout,'(a,a,a,i0)') "Correctness Check : ", __FILE__, " : LINE = ", __LINE__
    !$OMP TARGET UPDATE FROM(P_SFLUXZEN, P_TAUR, P_TAUG)
    !$ACC UPDATE HOST(P_SFLUXZEN) ASYNC(1)
    !$ACC UPDATE HOST(P_TAUR) ASYNC(1)
    !$ACC UPDATE HOST(P_TAUG) ASYNC(1)
    !$ACC WAIT(1)
    s1 = size(P_TAUG,1)/2
    s2 = size(P_TAUG,2)/2
    s3 = size(P_TAUG,3)/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "P_TAUG=", P_TAUG(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(P_TAUR,1)/2
    s2 = size(P_TAUR,2)/2
    s3 = size(P_TAUR,3)/2
    write(nulout,'(a,g0.5,a,i0,a,i0,a,i0)') "P_TAUR=", P_TAUR(s1,s2,s3), " at indices ", s1, " ", s2, " ", s3
    s1 = size(P_SFLUXZEN,1)/2
    s2 = size(P_SFLUXZEN,2)/2
    write(nulout,'(a,g0.5,a,i0,a,i0)') "P_SFLUXZEN=", P_SFLUXZEN(s1,s2), " at indices ", s1, " ", s2
    write(nulout,'(a)') "*******************************************************************"
#endif

#ifdef DEBUG
    write(nulout,'(a,a,a,i0,a)') "    ", __FILE__, " : LINE = ", __LINE__, " Done"
#endif

END SUBROUTINE SRTM_TAUMOL25
