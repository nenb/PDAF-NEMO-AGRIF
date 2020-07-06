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
  USE mod_statevector_pdaf, &   ! Assimilation variables
       ONLY: mpi_subd_lon_child, mpi_subd_lon_par, mpi_subd_lat_child, &
       mpi_subd_lat_par

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************

  ! Use horizontal localization, hence number of domains corresponds to
  ! number of (x,y) points on parent grid plus child grid (if AGRIF used).
  n_domains_p = (mpi_subd_lon_par*mpi_subd_lat_par)
#if defined key_agrif
  n_domains_p = n_domains_p +  (mpi_subd_lon_child*mpi_subd_lat_child)
#endif

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_n_domains_pdaf
