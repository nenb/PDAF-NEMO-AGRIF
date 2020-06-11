!$Id: collect_state_pdaf.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X or PDAF\_assimilate\_X
! after the propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code based on online implementation
! Later revisions - see svn log
!
! !USES:
  USE kind_pdaf
  USE mod_statevector, &
       ONLY: fill2d_statevector, fill3d_statevector
  
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL(pwp), INTENT(inout) :: state_p(dim_p)  ! local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X    (as U_coll_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
! Call Xfill_statevector
!EOP

! *** local variables ***


! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  CALL fill2d_statevector(dim_p, state_p)
  CALL fill3d_statevector(dim_p, state_p)

END SUBROUTINE collect_state_pdaf
