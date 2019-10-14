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

  WRITE(*,*) 'My proc ID is: ', PROC_ID

  ! ======== SEND ========
  IF (PROC_ID == 0) THEN ! Master create a message

    numberSND = 12345.44

    ! CALL MPI_SEND (data_to_send, send_count, send_type, destination_ID, tag, comm, ierr) 
    CALL MPI_SEND(numberSND, 1, MPI_DOUBLE_PRECISION, 1, 88, MPI_COMM_WORLD, ierr) 
    
    WRITE(*,*) 'Processor 0 just sent a message: ', numberSND

  ELSE IF (PROC_ID == 1) THEN ! A particular slave receives message
    CALL MPI_RECV (numberRCV, 1, MPI_DOUBLE_PRECISION, 0, 88, MPI_COMM_WORLD, STATUS, ierr) 
    WRITE(*,*) 'Processor 1 just received a message: ', numberRCV

  END IF

  ! Finalize MPI environment
  CALL MPI_FINALIZE(ierr)

END PROGRAM testalltoall
