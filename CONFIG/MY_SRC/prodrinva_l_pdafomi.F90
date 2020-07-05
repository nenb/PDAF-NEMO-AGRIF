!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA_l
!!
!! This routine calls the routine PDAFomi_prodRinvA_l
!! for each observation type
!!
SUBROUTINE prodrinva_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
!$AGRIF_DO_NOT_TREAT
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_prodRinvA_l
  ! Include observation types (rename generic name)
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

  ! Include variables for localization
  USE mod_assimilation_pdaf, ONLY: local_range, locweight, srange
  ! Include filter process rank
  USE mod_parallel_pdaf, ONLY: mype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p          !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l         !< Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              !< Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_l(dim_obs_l)  !< Local vector of observations
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix

! *** local variables ***
  INTEGER :: verbose                 ! Verbosity flag
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index


! **********************
! *** INITIALIZATION ***
! **********************

  ! Set verbosity flag (Screen output for first analysis domain)
  IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
     verbose = 1
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p


! *************************************
! *** Compute                       ***
! ***                  -1           ***
! ***           C = W R   A         ***
! ***                               ***
! *** where W are the localization  ***
! *** weights.                      ***
! *************************************

  CALL PDAFomi_prodRinvA_l(obs_A_l, obs_A, dim_obs_l, rank, &
       locweight, local_range, srange, A_l, C_l, verbose)
  CALL PDAFomi_prodRinvA_l(obs_B_l, obs_B, dim_obs_l, rank, &
       locweight, local_range, srange, A_l, C_l, verbose)
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE prodrinva_l_pdafomi
