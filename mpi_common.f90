! This module contains some MPI variables. Not sure it's needed...

MODULE mpi_common

   USE mpi
   
   IMPLICIT NONE
   
   INTEGER :: PROC_ID, N_MPI_THREADS, STATUS(MPI_STATUS_SIZE), ierr 
 
END MODULE mpi_common
