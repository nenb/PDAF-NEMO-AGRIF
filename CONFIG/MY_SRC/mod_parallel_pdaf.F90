!$Id: mod_parallel_pdaf.F90 1415 2013-09-25 14:33:26Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_parallel_pdaf

! !DESCRIPTION:
! This modules provides variables for the MPI parallelization
! to be shared between model-related routines. The are variables
! that are used in the model, even without PDAF and additional
! variables that are only used, if data assimialtion with PDAF
! is performed.
! In addition methods to initialize and finalize MPI are provided.
! The initialization routine is only for the model itself, the 
! more complex initialization of communicators for xecution with
! PDAF is peformed in init\_parallel\_pdaf.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_kind_pdaf

  IMPLICIT NONE
  SAVE 

  INCLUDE 'mpif.h'

! !PUBLIC DATA MEMBERS:
  ! Basic variables for model state integrations
  INTEGER :: COMM_model  ! MPI communicator for model tasks
  INTEGER :: mype_model  ! Number of PEs in COMM_model
  INTEGER :: npes_model  ! PE rank in COMM_model
  INTEGER :: COMM_ensemble      ! Communicator of all PEs doing model tasks
  INTEGER :: mype_ens, npes_ens ! rank and size in COMM_ensemble

  ! Additional variables for use with PDAF
  INTEGER :: n_modeltasks = 1         ! Number of parallel model tasks
  INTEGER :: n_filterpes  = 1         ! Number of PEs for filter analysis
  INTEGER :: npes_world, mype_world   ! # PEs and PE rank in MPI_COMM_WORLD
  INTEGER :: COMM_filter ! MPI communicator for filter PEs 
  INTEGER :: mype_filter, npes_filter ! # PEs and PE rank in COMM_filter
  INTEGER :: COMM_couple ! MPI communicator for coupling filter and model
  LOGICAL :: modelpe     ! Whether we are on a PE in a COMM_model
  LOGICAL :: filterpe    ! Whether we are on a PE in a COMM_filter
  INTEGER :: task_id     ! Index of my model task (1,...,n_modeltasks)
  INTEGER :: MPIerr      ! Error flag for MPI
  INTEGER :: MPIstatus(MPI_STATUS_SIZE)       ! Status array for MPI
  INTEGER :: screen_parallel = 1       ! Option for screen output
  INTEGER :: ens_test_parallel = 0     ! Ensemble size test - NOT USED
  INTEGER, ALLOCATABLE :: local_npes_model(:) ! # PEs per ensemble
!EOP
  
CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_parallel - Initialize MPI
!
! !INTERFACE:
  SUBROUTINE init_parallel()

! !DESCRIPTION:
! Routine to initialize MPI, the number of PEs
! (npes\_world) and the rank of a PE (mype\_world).
! The model is executed within the scope of the
! communicator Comm_model. It is also initialized
! here together with its size (npes\_model) and 
! the rank of a PE (mype\_model) within Comm_model.
!EOP

    IMPLICIT NONE

    INTEGER :: i
  
    CALL MPI_Init(i);
    CALL MPI_Comm_Size(MPI_COMM_WORLD,npes_world,i)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,i)
   
  END SUBROUTINE init_parallel
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: finalize_parallel - Finalize MPI
!
! !INTERFACE:
  SUBROUTINE finalize_parallel()

! !DESCRIPTION:
! Routine to finalize MPI
!EOP

    IMPLICIT NONE
    
    CALL  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
    CALL  MPI_Finalize(MPIerr)

  END SUBROUTINE finalize_parallel
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: abort_parallel - Abort MPI
!
! !INTERFACE:
  SUBROUTINE abort_parallel()

! !DESCRIPTION:
! Routine to abort MPI program
!EOP

    IMPLICIT NONE
    
    CALL  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  END SUBROUTINE abort_parallel
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_parallel_pdaf --- Initialize communicators for PDAF
!
! !INTERFACE:
  SUBROUTINE init_parallel_pdaf(dim_ens, screen, mpi_comm)

! !DESCRIPTION:
! Parallelization routine for a model with 
! attached PDAF. The subroutine is called in 
! the main program subsequently to the 
! initialization of MPI. It initializes
! MPI communicators for the model tasks, filter 
! tasks and the coupling between model and
! filter tasks. In addition some other variables 
! for the parallelization are initialized.
! The communicators and variables are handed
! over to PDAF in the call to 
! PDAF\_filter\_init.
!
! 3 Communicators are generated:\\
! - COMM\_filter: Communicator in which the
!   filter itself operates\\
! - COMM\_model: Communicators for parallel
!   model forecasts\\
! - COMM\_couple: Communicator for coupling
!   between models and filter\\
! Other variables that have to be initialized are:\\
! - filterpe - Logical: Does the PE execute the 
! filter?\\
! - my\_ensemble - Integer: The index of the PE's 
! model task\\
! - local\_npes\_model - Integer array holding 
! numbers of PEs per model task
!
! For COMM\_filter and COMM\_model also
! the size of the communicators (npes\_filter and 
! npes\_model) and the rank of each PE 
! (mype\_filter, mype\_model) are initialized. 
! These variables can be used in the model part 
! of the program, but are not handed over to PDAF.
!
! This variant is for a domain decomposed 
! model.
!
! This is a template that is expected to work 
! with many domain-decomposed models. However, 
! it might be necessary to adapt the routine 
! for a particular model. Inportant is that the
! communicator COMM_model equals the communicator 
! used in the model. If one plans to run the
! ensemble forecast in parallel COMM_model cannot 
! be MPI_COMM_WORLD! Thus, if the model uses 
! MPI_COMM_WORLD it has to be replaced by an 
! alternative communicator named, e.g., COMM_model.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    USE par_oce, ONLY: jpnij ! Number of model MPI processes

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(inout) :: dim_ens ! Ensemble size or number of EOFs (only SEEK)
    ! Often dim_ens=0 when calling this routine, because the real ensemble size
    ! is initialized later in the program. For dim_ens=0 no consistency check
    ! for ensemble size with number of model tasks is performed.
    INTEGER, INTENT(in)    :: screen ! Whether screen information is shown
    INTEGER, INTENT(inout) :: mpi_comm ! MPI communicator after XIOS splitting

! !CALLING SEQUENCE:
! Called by: main program
! Calls: MPI_Comm_size
! Calls: MPI_Comm_rank
! Calls: MPI_Comm_split
! Calls: MPI_Barrier
!EOP

! ! local variables
    INTEGER :: i, j               ! Counters
    INTEGER :: mype_couple, npes_couple ! Rank and size in COMM_couple
    INTEGER :: pe_index           ! Index of PE
    INTEGER :: my_color, color_couple ! Variables for communicator-splitting 
    CHARACTER(len=100) :: nmlfile   ! name of namelist file
    INTEGER :: tasks  ! Number of model tasks

    NAMELIST /tasks_nml/ tasks


    ! *** Read namelist for number of model tasks ***
    nmlfile ='namelist.pdaf'

    OPEN (20,file=nmlfile)
    READ (20,NML=tasks_nml)
    CLOSE (20)

    n_modeltasks=tasks

    ! ***              COMM_ENSEMBLE                ***
    ! *** Generate communicator for ensemble runs   ***
    ! *** only used to generate model communicators ***
    COMM_ensemble = mpi_comm
    CALL MPI_Comm_Size(COMM_ensemble, npes_ens, MPIerr)
    CALL MPI_Comm_Rank(COMM_ensemble, mype_ens, MPIerr)

    ! *** Initialize communicators for ensemble evaluations ***
    IF (mype_ens == 0) &
         WRITE (*, '(/1x, a)') 'Initialize communicators for assimilation with PDAF'

    ! *** Check consistency of number of parallel ensemble tasks ***
    consist1: IF (n_modeltasks*jpnij /= npes_ens) THEN
       WRITE (*,'(5x, a, i, i, i)') 'n_modeltasks, jpnij, npes_ens =',&
            n_modeltasks,jpnij,npes_ens
       WRITE (*,'(5x, a)') 'ERROR: Total number of processes is not consistent.'
       CALL abort_parallel()
    END IF consist1

    ! *** Store # PEs per ensemble                 ***
    ! *** used for info on PE 0 and for generation ***
    ! *** of model communicators on other Pes      ***
    ALLOCATE(local_npes_model(n_modeltasks))

    local_npes_model = FLOOR(REAL(npes_ens) / REAL(n_modeltasks))

    DO i = 1, (npes_world - n_modeltasks * local_npes_model(1))
       local_npes_model(i) = local_npes_model(i) + 1
    END DO


    ! ***              COMM_MODEL               ***
    ! *** Generate communicators for model runs ***
    ! *** (Split COMM_ENSEMBLE)                 ***
    pe_index = 0
    doens1: DO i = 1, n_modeltasks
       DO j = 1, local_npes_model(i)
          IF (mype_ens == pe_index) THEN
             task_id = i
             EXIT doens1
          END IF
          pe_index = pe_index + 1
       END DO
    END DO doens1

    CALL MPI_Comm_split(COMM_ensemble, task_id, mype_ens, &
         COMM_model, MPIerr)

    ! *** Re-initialize PE informations   ***
    ! *** according to model communicator ***
    CALL MPI_Comm_Size(COMM_model, npes_model, MPIerr)
    CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

    IF (screen > 1) then
       WRITE (*,*) 'MODEL: mype(w)= ', mype_ens, '; model task: ', task_id, &
            '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
    END IF

    ! Init flag FILTERPE (all PEs of model task 1)
    IF (task_id == 1) THEN
       filterpe = .TRUE.
    ELSE
       filterpe = .FALSE.
    END IF

    ! ***         COMM_FILTER                 ***
    ! *** Generate communicator for filter    ***
    IF (filterpe) THEN
       my_color = task_id
    ELSE
       my_color = MPI_UNDEFINED
    ENDIF

    CALL MPI_Comm_split(COMM_ensemble, my_color, mype_ens, &
         COMM_filter, MPIerr)

    ! *** Initialize PE informations         ***
    ! *** according to coupling communicator ***
    IF (filterpe) THEN
       CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
       CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)
    ENDIF


    ! ***              COMM_COUPLE                 ***
    ! *** Generate communicators for communication ***
    ! *** between model and filter PEs             ***
    ! *** (Split COMM_ENSEMBLE)                    ***

    color_couple = mype_model + 1

    CALL MPI_Comm_split(COMM_ensemble, color_couple, mype_ens, &
         COMM_couple, MPIerr)

    ! *** Initialize PE informations         ***
    ! *** according to coupling communicator ***
    CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
    CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

    IF (screen > 0) THEN
       IF (mype_ens == 0) THEN
          WRITE (*, '(/18x, a)') 'PE configuration:'
          WRITE (*, '(2x, a6, a9, a10, a14, a13, /2x, a5, a9, a7, a7, a7, a7, a7, /2x, a)') &
               'world', 'filter', 'model', 'couple', 'filterPE', &
               'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
               '----------------------------------------------------------'
       END IF
       CALL MPI_Barrier(COMM_ensemble, MPIerr)
       IF (task_id == 1) THEN
          WRITE (*, '(2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
               mype_ens, mype_filter, task_id, mype_model, color_couple, &
               mype_couple, filterpe
       ENDIF
       IF (task_id > 1) THEN
          WRITE (*,'(2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
               mype_ens, task_id, mype_model, color_couple, mype_couple, filterpe
       END IF
       CALL MPI_Barrier(COMM_ensemble, MPIerr)

       IF (mype_ens == 0) WRITE (*, '(/a)') ''

    END IF

    ! ******************************************************************************
    ! *** Initialize model equivalents to COMM_model, npes_model, and mype_model ***
    ! ******************************************************************************

    ! If the names of the variables for COMM_model, npes_model, and 
    ! mype_model are different in the numerical model, the 
    ! model-internal variables should be initialized at this point.

    mpi_comm = COMM_model

  END SUBROUTINE init_parallel_pdaf

  SUBROUTINE calc_global_dim(dim_local, dim_global)

    ! !DESCRIPTION:
    !  Computes the dimension of a state variable field for the global state vector.
    
    ! !USES:
    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in)    :: dim_local  ! Local dimension of state variable
    INTEGER, INTENT(inout) :: dim_global ! Global dimension of state variable

    ! ! local variables
    INTEGER :: domain                       ! Counters
    INTEGER :: MPIerr                       ! MPI error status
    INTEGER :: MPIstatus(MPI_STATUS_SIZE)   ! MPI status array
    INTEGER, ALLOCATABLE :: dim_p_array(:)  ! Temporary array for local dimensions

     IF (mype_filter /= 0) THEN
        ! *** Send dim_p on filter-PEs with rank > 0 ***
        CALL MPI_Send(dim_local, 1, MPI_INT, 0, 1, COMM_filter, MPIerr)
     ELSE
        ! Initialize array of dim_p values
        ALLOCATE(dim_p_array(npes_filter))
        dim_p_array(1) = dim_local
        DO domain = 2, npes_filter
           ! Receive dim_p on filter-PE with rank = 0
           CALL MPI_Recv(dim_p_array(domain), 1, MPI_INT, MPI_ANY_SOURCE, 1, &
                COMM_filter, MPIstatus, MPIerr)
        END DO
        dim_global = SUM(dim_p_array)
        DEALLOCATE(dim_p_array)
     END IF

     ! Send dim_state from filter-PE with rank = 0 to all ranks.
     CALL MPI_Bcast(dim_global, 1, MPI_INT, 0, COMM_filter, MPIerr)

   END SUBROUTINE calc_global_dim

   SUBROUTINE gather_ens(rank, dimx_local, dimy_local, dimz_local, dimx_global, &
        dimy_global, dimz_global, dim_local, offset, dim_ens, ens_local, ens)

!  !DESCRIPTION:
!  Gathers local ensemble arrays on PE = rank, then assembles global
!  ensemble array on PE = rank. This global ensemble is reshaped i.e.
!  it is no longer in 1D statevector format, but rather in its output
!  format of (x,y,member) or (x,y,z,member).

     ! !USES:
     IMPLICIT NONE

     ! !ARGUMENTS:
     INTEGER, INTENT(in) :: rank        ! Rank of PE to gather ens on
     INTEGER, INTENT(in) :: dimx_local  ! Dimension of mpi domain: x axis
     INTEGER, INTENT(in) :: dimy_local  ! Dimension of mpi domain: y axis
     INTEGER, INTENT(in) :: dimz_local  ! Dimension of mpi domain: z axis
     INTEGER, INTENT(in) :: dimx_global ! Dimension of global domain: x axis
     INTEGER, INTENT(in) :: dimy_global ! Dimension of global domain: y axis
     INTEGER, INTENT(in) :: dimz_global ! Dimension of global domain: z axis
     INTEGER, INTENT(in) :: dim_local   ! Dimension of local statevector
     INTEGER, INTENT(in) :: offset      ! State variable offset in local statevector
     INTEGER, INTENT(in) :: dim_ens     ! Dimension of ensemble
     REAL(pwp), INTENT(in)    :: ens_local(:,:) ! Local ensemble
     REAL(pwp), INTENT(inout) :: ens(:,:,:,:)   ! Global ensemble (reshaped)

     ! local variables
     INTEGER :: domain, member, i, j, k      ! Counters
     INTEGER :: MPIerr                       ! MPI error status
     INTEGER :: MPIstatus(MPI_STATUS_SIZE)   ! MPI status array
     INTEGER :: dimx_tmp                     ! Dimension of mpi domain: x axis
     INTEGER :: dimy_tmp                     ! Dimension of mpi domain: x axis
     INTEGER :: dimz_tmp                     ! Dimension of mpi domain: x axis
     INTEGER :: offx                         ! Offset of local in global: x axis
     INTEGER :: offy                         ! Offset of local in global: y axis
     INTEGER :: dim_msg                      ! Size of MPI msg
     REAL(pwp), ALLOCATABLE :: flat_ens_local(:)  ! Flattened ens_local array


     mype: IF (mype_filter /= rank) THEN
        ! Send arrays on filter-PEs /= rank to PE = rank.
        ! Flatten arrays and also 'attach' local x/y/z dimensions as values
        ! at the end of the flattened arrays (an ugly hack). These 'attached'
        ! values are needed to reassemble the array on PE = rank.
        ALLOCATE(flat_ens_local(dim_ens*dim_local + 3))
        DO member = 1, dim_ens
           DO i = 1, dim_local
              flat_ens_local(i+dim_local*(member-1) ) = &
                   ens_local(i+offset, member)
           END DO
        END DO
        flat_ens_local(dim_ens*dim_local + 1) = REAL(dimx_local)
        flat_ens_local(dim_ens*dim_local + 2) = REAL(dimy_local)
        flat_ens_local(dim_ens*dim_local + 3) = REAL(dimz_local)
        CALL MPI_Send(flat_ens_local(1), dim_ens*dim_local + 3, &
             MPI_DOUBLE_PRECISION, rank, 1, COMM_filter, MPIerr)
        DEALLOCATE(flat_ens_local)
     ELSE mype ! Assemble each local ensemble in global ensemble array.
        ! Offset for local ensemble in global ensemble.
        offx=0
        offy=0
        DO domain = 1, npes_filter
           ! No receive necessary on PE = rank (domain = rank+1)
           IF (domain == rank+1) THEN
              DO member = 1, dim_ens
                 DO k = 1, dimz_local
                    DO j = 1, dimy_local
                       DO i = 1, dimx_local
                          ens(i+offx,j+offy,k,member) = ens_local(offset+i+&
                               (j-1)*dimx_local+(k-1)*dimx_local*dimy_local,member)
                       END DO
                    END DO
                 END DO
              END DO
              ! Increment offset for local in global.
              IF (offx+dimx_local >= dimx_global) offy=offy+dimy_local
              offx=MOD(offx+dimx_local,dimx_global)
           ELSE ! Receive all messages on PE = rank
              CALL MPI_Probe(domain - 1, 1, COMM_filter, MPIstatus, MPIerr)
              CALL MPI_Get_count(MPIstatus, MPI_DOUBLE_PRECISION, dim_msg, MPIerr)
              ALLOCATE(flat_ens_local(dim_msg))
              CALL MPI_Recv(flat_ens_local(1), dim_msg, MPI_DOUBLE_PRECISION, &
                   domain - 1, 1, COMM_filter, MPIstatus, MPIerr)
              dimx_tmp = INT(flat_ens_local(dim_msg-2))
              dimy_tmp = INT(flat_ens_local(dim_msg-1))
              dimz_tmp = INT(flat_ens_local(dim_msg))
              DO member = 1, dim_ens
                 DO k = 1, dimz_tmp
                    DO j = 1, dimy_tmp
                       DO i = 1, dimx_tmp
                          ens(i+offx,j+offy,k,member) = flat_ens_local(i+ &
                               (j-1)*dimx_tmp+(k-1)*dimx_tmp*dimy_tmp+ &
                               (member-1)*dimx_tmp*dimy_tmp*dimz_tmp)
                       END DO
                    END DO
                 END DO
              END DO
              DEALLOCATE(flat_ens_local)
              ! Increment offset for local in global
              IF (offx+dimx_tmp >= dimx_global) offy=offy+dimy_tmp
              offx=MOD(offx+dimx_tmp,dimx_global)
           END IF
        END DO
     END IF mype

   END SUBROUTINE gather_ens

END MODULE mod_parallel_pdaf
