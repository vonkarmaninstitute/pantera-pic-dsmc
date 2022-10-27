PROGRAM PANTERA 

   USE particle
   USE mpi_common
   USE global
   USE screen
   USE tools
   USE initialization
   USE timecycle
   USE grid_and_partition
   USE mUMFPACK
   USE fields
   USE fully_implicit

   IMPLICIT NONE
   
   ! ========= Init MPI environment ========================
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_THREADS, ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_ID, ierr)
   CALL NEWTYPE ! Define new "particle" type for mpi


   CALL ONLYMASTERPRINT1(PROC_ID, '> TESTING PETSC INIT...')
   !CALL PETSC_INITIAL_TEST
   CALL PETSC_INIT
   CALL ONLYMASTERPRINT1(PROC_ID, '> PETSC INIT DONE!')
   ! ========= Print header (important things first) =======
   CALL PRINTTITLE(PROC_ID)

   ! ========= Read input file and init variables ==========
   CALL ONLYMASTERPRINT1(PROC_ID, '> READING INPUT DATA...')
   CALL READINPUT          ! Read input file
   CALL INITVARIOUS        ! Initialize some additional variables
   CALL INITINJECTION      ! Initialize variables for injection
   CALL INITCOLLISIONS     ! Initialize variables for collisions
   CALL INITREACTIONS      ! Initialize variables for reactions
   CALL INITFIELDS         ! Initialize electromagnetic fields
   CALL COMPUTE_B_FIELD_FROM_SOLENOIDS

   IF (PIC_TYPE .NE. NONE) CALL ASSEMBLE_POISSON
   ! ========= Initial particles seed ======================
   IF (RESTART_TIMESTEP > 0) THEN
      CALL READ_PARTICLES_FILE(RESTART_TIMESTEP)
      tID = RESTART_TIMESTEP
   ELSE
      tID = 0
      RESTART_TIMESTEP = 0
      CALL INITIAL_SEED
   END IF
   ! CALL DUMP_PARTICLES_SCREEN

   ! ========= Time loop ===================================
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
   CALL ONLYMASTERPRINT1(PROC_ID, '> STARTING TIME LOOP...')
   CALL TIME_LOOP

   ! CALL DUMP_PARTICLES_SCREEN

   ! ========== Close MPI environment ======================
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

   CALL ONLYMASTERPRINT1(PROC_ID, '> R O A R !')
   CALL MPI_FINALIZE(ierr)

END PROGRAM PANTERA
