PROGRAM testalltoall

  USE mpi

  IMPLICIT NONE

  INTEGER :: PROC_ID, N_MPI_THREADS, STATUS(MPI_STATUS_SIZE), ierr

  REAL(KIND=8) :: numberSND, numberRCV

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: snd_vect, rcv_vect
  INTEGER :: num_of_numbers, ii

  WRITE(*,*) 'Ciao'


  ! ========= Init MPI environment ========================
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_THREADS, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_ID, ierr)

  
  ! Master 
  if (PROC_ID == 0) then

    num_of_numbers = 5 ! I will send 5 numbers to everyone

    allocate(snd_vect(N_MPI_THREADS*num_of_numbers))

    do ii = 1,num_of_numbers*N_MPI_THREADS
      snd_vect(ii) = ii
    end do
  end if

  ! Master broadcasts the number of elements that he will send
  call MPI_BCAST(num_of_numbers, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)

  WRITE(*,*) 'Hey, I am process ', PROC_ID, ' and I expect ', num_of_numbers, ' elements' 
  allocate(rcv_vect(num_of_numbers))

  ! ------------
  call MPI_SCATTER(snd_vect, num_of_numbers, MPI_DOUBLE_PRECISION,  &
                   rcv_vect, num_of_numbers, MPI_DOUBLE_PRECISION,  &
                   0, MPI_COMM_WORLD, ierr) 

  do ii = 1, num_of_numbers
    write(*,*) 'Proc ', PROC_ID, ' says: ', ii, rcv_vect(ii)
  end do

  ! Finalize MPI environment
  CALL MPI_FINALIZE(ierr)

END PROGRAM testalltoall
