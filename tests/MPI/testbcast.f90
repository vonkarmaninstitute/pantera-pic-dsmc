PROGRAM testalltoall

  USE mpi

  IMPLICIT NONE

  INTEGER :: PROC_ID, N_MPI_THREADS, STATUS(MPI_STATUS_SIZE), ierr

  REAL(KIND=8) :: numberSND, numberRCV

  WRITE(*,*) 'Ciao'


  ! ========= Init MPI environment ========================
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_THREADS, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_ID, ierr)

  if (PROC_ID == 0) numberSND = 1936.27d0

  WRITE(*,*) 'My proc ID is: ', PROC_ID, ' and I say: ', numberSND

  call MPI_BCAST (numberSND, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  WRITE(*,*) 'My proc ID is: ', PROC_ID, ' and I say: ', numberSND

  ! Finalize MPI environment
  CALL MPI_FINALIZE(ierr)

END PROGRAM testalltoall
