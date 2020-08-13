MODULE mod_util_pdaf
!$AGRIF_DO_NOT_TREAT

CONTAINS

  !$Id: init_info_pdaf.F90 1589 2015-06-12 11:57:58Z lnerger $
  !BOP
  !
  ! !ROUTINE: init_info_pdaf - Screen output on assimilation configuration
  !
  ! !INTERFACE:
  SUBROUTINE init_info_pdaf()

    ! !DESCRIPTION:
    ! This routine performs a model-sided screen output about
    ! the coniguration of the data assimilation system.
    ! Using this output is optional. Most of the information
    ! is also displayed by PDAF itself when it is initialized
    ! in PDAF_init. Not displayed by PDAF is the assimilation
    ! interval (delt_obs), which is unknown to PDAF.
    !
    ! !REVISION HISTORY:
    ! 2011-05 - Lars Nerger - Initial code extracted from init_pdaf
    ! Later revisions - see svn log
    !
    ! !USES:
    USE mod_assimilation_pdaf, & ! Variables for assimilation
         ONLY: filtertype, subtype, dim_ens, delt_obs, model_error, &
         model_err_amp, forget, rank_analysis_enkf, int_rediag

    IMPLICIT NONE

    ! !CALLING SEQUENCE:
    ! Called by: init_pdaf
    !EOP


    ! *****************************
    ! *** Initial Screen output ***
    ! *****************************

    IF (filtertype == 0) THEN
       WRITE (*, '(/21x, a)') 'Filter: SEEK'
       IF (subtype == 2) THEN
          WRITE (*, '(6x, a)') '-- fixed basis filter with update of matrix U'
          WRITE (*, '(6x, a)') '-- no re-diagonalization of VUV^T'
       ELSE IF (subtype == 3) THEN
          WRITE (*, '(6x, a)') '-- fixed basis filter & no update of matrix U'
          WRITE (*, '(6x, a)') '-- no re-diagonalization of VUV^T'
       ELSE IF (subtype == 5) THEN
          WRITE (*, '(6x, a)') '-- Offline mode'
       END IF
       WRITE (*, '(13x, a, i5)') 'number of EOFs:', dim_ens
       IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
       WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
       IF (subtype /= 5) THEN
          IF ((int_rediag > 0) .AND. ((subtype /= 2) .OR. (subtype /= 3))) &
               WRITE (*, '(10x, a, i4, a)') &
               'Re-diag each ', int_rediag, '-th analysis step'
       ELSE
          IF (int_rediag == 1) THEN
             WRITE (*, '(10x, a)') 'Perform re-diagonalization'
          ELSE
             WRITE (*, '(10x, a)') 'No re-diagonalization'
          END IF
       END IF
    ELSE IF (filtertype == 1) THEN
       WRITE (*, '(21x, a)') 'Filter: SEIK'
       IF (subtype == 2) THEN
          WRITE (*, '(6x, a)') '-- fixed error-space basis'
       ELSE IF (subtype == 3) THEN
          WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
       ELSE IF (subtype == 4) THEN
          WRITE (*, '(6x, a)') '-- use ensemble transformation'
       ELSE IF (subtype == 5) THEN
          WRITE (*, '(6x, a)') '-- Offline mode'
       END IF
       WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
       IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
       WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
       IF (model_error) THEN
          WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
       END IF
    ELSE IF (filtertype == 2) THEN
       WRITE (*, '(21x, a)') 'Filter: EnKF'
       IF (subtype == 5) THEN
          WRITE (*, '(6x, a)') '-- Offline mode'
       END IF
       WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
       IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
       WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
       IF (model_error) THEN
          WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
       END IF
       IF (rank_analysis_enkf > 0) THEN
          WRITE (*, '(6x, a, i5)') &
               'analysis with pseudo-inverse of HPH, rank:', rank_analysis_enkf
       END IF
    ELSE IF (filtertype == 3) THEN
       WRITE (*, '(21x, a)') 'Filter: LSEIK'
       IF (subtype == 2) THEN
          WRITE (*, '(6x, a)') '-- fixed error-space basis'
       ELSE IF (subtype == 3) THEN
          WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
       ELSE IF (subtype == 4) THEN
          WRITE (*, '(6x, a)') '-- use ensemble transformation'
       ELSE IF (subtype == 5) THEN
          WRITE (*, '(6x, a)') '-- Offline mode'
       END IF
       WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
       IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
       WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
       IF (model_error) THEN
          WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
       END IF
    ELSE IF (filtertype == 4) THEN
       WRITE (*, '(21x, a)') 'Filter: ETKF'
       IF (subtype == 0) THEN
          WRITE (*, '(6x, a)') '-- Variant using T-matrix'
       ELSE IF (subtype == 1) THEN
          WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
       ELSE IF (subtype == 5) THEN
          WRITE (*, '(6x, a)') '-- Offline mode'
       END IF
       WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
       IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
       WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
       IF (model_error) THEN
          WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
       END IF
    ELSE IF (filtertype == 5) THEN
       WRITE (*, '(21x, a)') 'Filter: LETKF'
       IF (subtype == 0) THEN
          WRITE (*, '(6x, a)') '-- Variant using T-matrix'
       ELSE IF (subtype == 1) THEN
          WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
       ELSE IF (subtype == 5) THEN
          WRITE (*, '(6x, a)') '-- Offline mode'
       END IF
       WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
       IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
       WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
       IF (model_error) THEN
          WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
       END IF
    ELSE IF (filtertype == 6) THEN
       WRITE (*, '(21x, a)') 'Filter: ESTKF'
       IF (subtype == 0) THEN
          WRITE (*, '(6x, a)') '-- Standard mode'
       ELSE IF (subtype == 5) THEN
          WRITE (*, '(6x, a)') '-- Offline mode'
       END IF
       WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
       IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
       WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
       IF (model_error) THEN
          WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
       END IF
    ELSE IF (filtertype == 7) THEN
       WRITE (*, '(21x, a)') 'Filter: LESTKF'
       IF (subtype == 0) THEN
          WRITE (*, '(6x, a)') '-- Standard mode'
       ELSE IF (subtype == 5) THEN
          WRITE (*, '(6x, a)') '-- Offline mode'
       END IF
       WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
       IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
       WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
       IF (model_error) THEN
          WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
       END IF
    END IF

  END SUBROUTINE init_info_pdaf

  !$Id: read_config_pdaf.F90 3 2013-09-05 10:28:51Z lnerger $
  !BOP
  !
  ! !ROUTINE: read_config_pdaf - Read configuration for PDAF
  !
  ! !INTERFACE:
  SUBROUTINE read_config_pdaf()

    ! !DESCRIPTION:
    ! This routine read the namelist file with
    ! parameters controlling data assimilation with
    ! PDAF.

    ! !USES:
    USE mod_parallel_pdaf, &
         ONLY: mype_ens
    USE mod_assimilation_pdaf, &
         ONLY: filtertype, subtype, dim_ens, delt_obs, child_dt_fac, &
         screen, forget, local_range_par, local_range_child,&
         locweight, srange_par, srange_child, istate_t_par,&
         istate_u_par, istate_v_par, istate_t_child,&
         istate_u_child, istate_v_child
    USE mod_agrif_pdaf, &
         ONLY: lowlim_ssh_par, upplim_ssh_par, lowlim_sal_par, &
         upplim_sal_par, lowlim_temp_par, upplim_temp_par, &
         lowlim_uvel_par, upplim_uvel_par, lowlim_vvel_par, &
         upplim_vvel_par, lowlim_ssh_child, upplim_ssh_child, &
         lowlim_sal_child, upplim_sal_child, lowlim_temp_child, &
         upplim_temp_child, lowlim_uvel_child, upplim_uvel_child, &
         lowlim_vvel_child, upplim_vvel_child
    USE mod_obs_ssh_par_pdafomi, &
         ONLY: assim_ssh_par, rms_ssh_par, file_ssh_par, &
         twin_exp_ssh_par, noise_amp_ssh_par
    USE mod_obs_ssh_child_pdafomi, &
         ONLY: assim_ssh_child, rms_ssh_child, file_ssh_child, &
         twin_exp_ssh_child, noise_amp_ssh_child
    USE mod_obs_fake_ssh_par_pdafomi, &
         ONLY: assim_fake_ssh_par, rms_fake_ssh_par, file_fake_ssh_par, &
         twin_exp_fake_ssh_par, noise_amp_fake_ssh_par
    USE mod_obs_fake_ssh_child_pdafomi, &
         ONLY: assim_fake_ssh_child, rms_fake_ssh_child, file_fake_ssh_child, &
         twin_exp_fake_ssh_child, noise_amp_fake_ssh_child

    IMPLICIT NONE
    !EOP

    ! Local variables
    CHARACTER(len=100) :: nmlfile   ! name of namelist file

    NAMELIST /pdaf_nml/ filtertype, subtype, dim_ens,&
         delt_obs, child_dt_fac, screen, forget,&
         local_range_par, local_range_child, locweight, &
         srange_par, srange_child, istate_t_par,&
         istate_u_par, istate_v_par, istate_t_child,&
         istate_u_child, istate_v_child, assim_ssh_par, &
         rms_ssh_par, file_ssh_par, twin_exp_ssh_par, &
         noise_amp_ssh_par, lowlim_ssh_par, upplim_ssh_par, &
         lowlim_sal_par, upplim_sal_par, lowlim_temp_par, &
         upplim_temp_par, lowlim_uvel_par, upplim_uvel_par, &
         lowlim_vvel_par, upplim_vvel_par, assim_ssh_child, &
         rms_ssh_child, file_ssh_child, twin_exp_ssh_child, &
         noise_amp_ssh_child, lowlim_ssh_child, upplim_ssh_child, &
         lowlim_sal_child, upplim_sal_child, lowlim_temp_child, &
         upplim_temp_child, lowlim_uvel_child, upplim_uvel_child, &
         lowlim_vvel_child, upplim_vvel_child, assim_fake_ssh_par, &
         rms_fake_ssh_par, file_fake_ssh_par, twin_exp_fake_ssh_par, &
         noise_amp_fake_ssh_par, assim_fake_ssh_child, &
         rms_fake_ssh_child, file_fake_ssh_child, twin_exp_fake_ssh_child, &
         noise_amp_fake_ssh_child

    ! ****************************************************
    ! ***   Initialize PDAF parameters from namelist   ***
    ! ****************************************************

    ! *** Read namelist file ***
    nmlfile ='namelist.pdaf'

    OPEN (20,file=nmlfile)
    READ (20,NML=pdaf_nml)
    CLOSE(20)

    ! *** Print configuration variables ***
    showconf: IF (mype_ens .EQ. 0) THEN

       WRITE (*,'(/1x,a)') '-- Overview of PDAF configuration --'
       WRITE (*,'(3x,a)') 'PDAF [pdaf_nml]:'
       WRITE (*,'(5x,a,i10)')    'filtertype   ', filtertype
       WRITE (*,'(5x,a,i10)')    'subtype      ', subtype
       WRITE (*,'(5x,a,i10)')    'dim_ens      ', dim_ens
       WRITE (*,'(5x,a,i10)')    'delt_obs     ', delt_obs
       WRITE (*,'(5x,a,i10)')    'child_dt_fac ', child_dt_fac
       WRITE (*,'(5x,a,i10)')    'screen       ', screen
       WRITE (*,'(5x,a,f10.2)')  'forget       ', forget
       WRITE (*,'(5x,a,es10.2)') 'local_range_par ', local_range_par
       WRITE (*,'(5x,a,es10.2)') 'local_range_child ', local_range_child
       WRITE (*,'(5x,a,i10)')    'locweight    ', locweight
       WRITE (*,'(5x,a,es10.2)') 'srange_par  ', srange_par
       WRITE (*,'(5x,a,es10.2)') 'srange_child ', srange_child
       WRITE (*,'(5x,a,a)')  'istate_t_par   ', istate_t_par
       WRITE (*,'(5x,a,a)')  'istate_u_par   ', istate_u_par
       WRITE (*,'(5x,a,a)')  'istate_v_par   ', istate_v_par
       WRITE (*,'(5x,a,a)')  'istate_t_child  ', istate_t_child
       WRITE (*,'(5x,a,a)')  'istate_u_child  ', istate_u_child
       WRITE (*,'(5x,a,a)')  'istate_v_child  ', istate_v_child
       WRITE (*,'(5x,a,l1)') 'assim_ssh_par   ', assim_ssh_par
       WRITE (*,'(5x,a,es10.2)') 'rms_ssh_par ', rms_ssh_par
       WRITE (*,'(5x,a,a)')  'file_ssh_par    ', file_ssh_par
       WRITE (*,'(5x,a,l1)')     'twin_exp_ssh_par  ', twin_exp_ssh_par
       WRITE (*,'(5x,a,es10.2)') 'noise_amp_ssh_par ', noise_amp_ssh_par
       WRITE (*,'(5x,a,es10.2)') 'lowlim_ssh_par ', lowlim_ssh_par
       WRITE (*,'(5x,a,es10.2)') 'upplim_ssh_par ', upplim_ssh_par
       WRITE (*,'(5x,a,es10.2)') 'lowlim_sal_par ', lowlim_sal_par
       WRITE (*,'(5x,a,es10.2)') 'upplim_sal_par ', upplim_sal_par
       WRITE (*,'(5x,a,es10.2)') 'lowlim_temp_par ', lowlim_temp_par
       WRITE (*,'(5x,a,es10.2)') 'upplim_temp_par ', upplim_temp_par
       WRITE (*,'(5x,a,es10.2)') 'lowlim_uvel_par ', lowlim_uvel_par
       WRITE (*,'(5x,a,es10.2)') 'upplim_uvel_par ', upplim_uvel_par
       WRITE (*,'(5x,a,es10.2)') 'lowlim_vvel_par ', lowlim_vvel_par
       WRITE (*,'(5x,a,es10.2)') 'upplim_vvel_par ', upplim_vvel_par
       WRITE (*,'(5x,a,l1)') 'assim_ssh_child   ', assim_ssh_child
       WRITE (*,'(5x,a,es10.2)') 'rms_ssh_child ', rms_ssh_child
       WRITE (*,'(5x,a,a)')  'file_ssh_child    ', file_ssh_child
       WRITE (*,'(5x,a,l1)')     'twin_exp_ssh_child  ', twin_exp_ssh_child
       WRITE (*,'(5x,a,es10.2)') 'noise_amp_ssh_child ', noise_amp_ssh_child
       WRITE (*,'(5x,a,es10.2)') 'lowlim_ssh_child ', lowlim_ssh_child
       WRITE (*,'(5x,a,es10.2)') 'upplim_ssh_child ', upplim_ssh_child
       WRITE (*,'(5x,a,es10.2)') 'lowlim_sal_child ', lowlim_sal_child
       WRITE (*,'(5x,a,es10.2)') 'upplim_sal_child ', upplim_sal_child
       WRITE (*,'(5x,a,es10.2)') 'lowlim_temp_child ', lowlim_temp_child
       WRITE (*,'(5x,a,es10.2)') 'upplim_temp_child ', upplim_temp_child
       WRITE (*,'(5x,a,es10.2)') 'lowlim_uvel_child ', lowlim_uvel_child
       WRITE (*,'(5x,a,es10.2)') 'upplim_uvel_child ', upplim_uvel_child
       WRITE (*,'(5x,a,es10.2)') 'lowlim_vvel_child ', lowlim_vvel_child
       WRITE (*,'(5x,a,es10.2)') 'upplim_vvel_child ', upplim_vvel_child
       WRITE (*,'(5x,a,l1)') 'assim_fake_ssh_par   ', assim_fake_ssh_par
       WRITE (*,'(5x,a,es10.2)') 'rms_fake_ssh_par ', rms_fake_ssh_par
       WRITE (*,'(5x,a,a)')  'file_fake_ssh_par    ', file_fake_ssh_par
       WRITE (*,'(5x,a,l1)')     'twin_exp_fake_ssh_par  ', twin_exp_fake_ssh_par
       WRITE (*,'(5x,a,es10.2)') 'noise_amp_fake_ssh_par ', noise_amp_fake_ssh_par
       WRITE (*,'(5x,a,l1)') 'assim_fake_ssh_child   ', assim_fake_ssh_child
       WRITE (*,'(5x,a,es10.2)') 'rms_fake_ssh_child ', rms_fake_ssh_child
       WRITE (*,'(5x,a,a)')  'file_fake_ssh_child    ', file_fake_ssh_child
       WRITE (*,'(5x,a,l1)')     'twin_exp_fake_ssh_child  ', twin_exp_fake_ssh_child
       WRITE (*,'(5x,a,es10.2)') 'noise_amp_fake_ssh_child ', noise_amp_fake_ssh_child
       WRITE (*,'(1x,a)') '-- End of PDAF configuration overview --'

    END IF showconf

  END SUBROUTINE read_config_pdaf

  SUBROUTINE cleanup_pdaf()

    ! !DESCRIPTION:
    ! Clean-up at end of PDAF.

#if defined key_agrif
    USE mod_agrif_pdaf, &
         ONLY: nimppt_child, nimppt_par, njmppt_child, njmppt_par, &
         glamt_par, glamt_child, gphit_par, gphit_child, tmask_par, &
         tmask_child
    USE mod_statevector_pdaf, &
         ONLY: halo_2d_par, halo_2d_child
#else
    USE mod_agrif_pdaf, &
         ONLY: nimppt_par, njmppt_par, glamt_par, gphit_par, tmask_par
    USE mod_statevector_pdaf, &
         ONLY: halo_2d_par
#endif

    IMPLICIT NONE
    !EOP


    ! Clean-up from mod_agrif_pdaf
#if defined key_agrif
    IF (ALLOCATED(nimppt_par)) DEALLOCATE(nimppt_par)
    IF (ALLOCATED(njmppt_par)) DEALLOCATE(njmppt_par)
    IF (ALLOCATED(gphit_par)) DEALLOCATE(gphit_par)
    IF (ALLOCATED(glamt_par)) DEALLOCATE(glamt_par)
    IF (ALLOCATED(tmask_par)) DEALLOCATE(tmask_par)
    IF (ALLOCATED(halo_2d_par)) DEALLOCATE(halo_2d_par)
    IF (ALLOCATED(nimppt_child)) DEALLOCATE(nimppt_child)
    IF (ALLOCATED(njmppt_child)) DEALLOCATE(njmppt_child)
    IF (ALLOCATED(gphit_child)) DEALLOCATE(gphit_child)
    IF (ALLOCATED(glamt_child)) DEALLOCATE(glamt_child)
    IF (ALLOCATED(tmask_child)) DEALLOCATE(tmask_child)
    IF (ALLOCATED(halo_2d_child)) DEALLOCATE(halo_2d_child)
#else
    IF (ALLOCATED(nimppt_par)) DEALLOCATE(nimppt_par)
    IF (ALLOCATED(njmppt_par)) DEALLOCATE(njmppt_par)
    IF (ALLOCATED(gphit_par)) DEALLOCATE(gphit_par)
    IF (ALLOCATED(glamt_par)) DEALLOCATE(glamt_par)
    IF (ALLOCATED(tmask_par)) DEALLOCATE(tmask_par)
    IF (ALLOCATED(halo_2d_par)) DEALLOCATE(halo_2d_par)
#endif

  END SUBROUTINE cleanup_pdaf

  !$Id: finalize_pdaf.F90 1857 2017-12-14 18:26:01Z lnerger $
  !BOP
  !
  ! !ROUTINE: finalize_pdaf --- Finalize PDAF
  !
  ! !INTERFACE:
  SUBROUTINE finalize_pdaf()

    ! !DESCRIPTION:
    ! Timing and clean-up of PDAF
    !
    ! !REVISION HISTORY:
    ! 2004-11 - Lars Nerger - Initial code
    ! Later revisions - see svn log
    !
    ! !USES:
    USE mod_parallel_pdaf, &
         ONLY: mype_ens

    IMPLICIT NONE

    ! !CALLING SEQUENCE:
    ! Called by: main program
    !EOP

    ! *** Show allocated memory for PDAF ***
    IF (mype_ens==0) CALL PDAF_print_info(2)

    ! *** Print PDAF timings onto screen ***
    IF (mype_ens==0) CALL PDAF_print_info(1)

    ! *** Deallocate PDAF arrays
    CALL PDAF_deallocate()

  END SUBROUTINE finalize_pdaf

!$AGRIF_END_DO_NOT_TREAT
END MODULE mod_util_pdaf
