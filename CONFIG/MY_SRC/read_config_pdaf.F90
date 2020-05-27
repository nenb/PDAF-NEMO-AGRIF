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
  USE mod_assimilation, &
       ONLY: filtertype, subtype, dim_ens, delt_obs, screen,&
       incremental, type_forget, forget, local_range,&
       locweight, srange, rms_obs,type_trans,&
       type_sqrt, covartype, rank_analysis_enkf, istate_fname_t,&
       istate_fname_u, istate_fname_v, istate_fname_w

  IMPLICIT NONE
!EOP

! Local variables
  CHARACTER(len=100) :: nmlfile   ! name of namelist file

  NAMELIST /pdaf_nml/ filtertype, subtype, dim_ens,& 
       delt_obs, screen,incremental, type_forget, forget,&
       local_range, locweight, srange, rms_obs,type_trans,&
       type_sqrt, covartype, rank_analysis_enkf, istate_fname_t,&
       istate_fname_u, istate_fname_v, istate_fname_w

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
    WRITE (*,'(5x,a,i10)')    'screen       ', screen
    WRITE (*,'(5x,a,i10)')    'incremental  ', incremental
    WRITE (*,'(5x,a,i10)')    'type_forget  ', type_forget
    WRITE (*,'(5x,a,f10.2)')  'forget       ', forget
    WRITE (*,'(5x,a,i10)')    'type_trans   ', type_trans
    WRITE (*,'(5x,a,es10.2)') 'local_range  ', local_range
    WRITE (*,'(5x,a,i10)')    'locweight    ', locweight
    WRITE (*,'(5x,a,es10.2)') 'srange       ', srange
    WRITE (*,'(5x,a,es10.2)') 'rms_obs_ssh  ', rms_obs
    WRITE (*,'(5x,a,a)')  'istate_fname_t   ', istate_fname_t
    WRITE (*,'(5x,a,a)')  'istate_fname_u   ', istate_fname_u
    WRITE (*,'(5x,a,a)')  'istate_fname_v   ', istate_fname_v
    WRITE (*,'(5x,a,a)')  'istate_fname_w   ', istate_fname_w
    WRITE (*,'(1x,a)') '-- End of PDAF configuration overview --'

  END IF showconf

END SUBROUTINE read_config_pdaf
