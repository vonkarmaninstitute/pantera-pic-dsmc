PROGRAM PANTERA 

   USE particle
   USE mpi_common
   USE global
   USE screen
   USE tools
   USE initialization
   USE timecycle
   USE grid_and_partition

   IMPLICIT NONE
   integer :: ip
   ! ========= Init MPI environment ========================
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_THREADS, ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_ID, ierr)
   CALL NEWTYPE ! Define new "particle" type for mpi

   ! ========= Print header (important things first) =======
   CALL PRINTTITLE(PROC_ID)

   ! ========= Read input file and init variables ==========
   CALL ONLYMASTERPRINT1(PROC_ID, '> READING INPUT DATA...')
   CALL READINPUT          ! Read input file
   CALL INITVARIOUS        ! Initialize some additional variables
   CALL INITINJECTION      ! Initialize variables for injection
   CALL INITCOLLISIONS     ! Initialize variables for collisions
   CALL INITREACTIONS      ! Initialize variables for reactions

   ! ========= Initial particles seed ======================
   IF (BOOL_INITIAL_SEED)   CALL INITIAL_SEED
   ! CALL DUMP_PARTICLES_SCREEN
   do ip = 1,np_proc
   if (particles(ip)%x .ge. walls(1)%x1 .and. particles(ip)%x .le. walls(1)%x2 .and. particles(ip)%y .ge. walls(4)%y2 .and. &
   particles(ip)%y .le. walls(4)%y1) then
   call remove_particle_array(ip, particles, np_proc)
   end if
   end do
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
