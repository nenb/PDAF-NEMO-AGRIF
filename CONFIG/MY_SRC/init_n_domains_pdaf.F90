!$Id: init_n_domains_pdaf.F90 343 2020-01-21 14:21:42Z lnerger $
!> Set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF
!!
!! The routine is called in PDAF_X_update 
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains. 
!! It has to set the number of local analysis 
!! domains for the PE-local domain.
!!
!! Implementation for the 2D online example
!! without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)
!$AGRIF_DO_NOT_TREAT
  USE mod_assimilation_pdaf, &   ! Assimilation variables
       ONLY: dim_state_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************
  
  ! Here simply the state dimension
  n_domains_p = dim_state_p
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_n_domains_pdaf
