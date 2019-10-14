MODULE grid_and_partition 

   USE global
   USE mpi_common

   IMPLICIT NONE
 
   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CELL_FROM_POSITION -> finds the ID of a grid cell from particle position !!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE CELL_FROM_POSITION(XP,YP,  IDCELL)

      ! Note: first cell is denoted as zero

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: XP, YP ! Location of particle
      INTEGER, INTENT(OUT)     :: IDCELL ! ID of cell to which the particle belongs

      REAL(KIND=8) :: DX, DY

      ! Cartesian grid with equally spaced cells
      DX = (XMAX - XMIN)/NX
      DY = (YMAX - YMIN)/NY

      IDCELL = INT((XP-XMIN)/DX) + NX*INT((YP-YMIN)/DY)
  
   END SUBROUTINE CELL_FROM_POSITION

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE PROC_FROM_POSITION -> finds the process ID from particle position !!!!!!!!!!!!!!
   ! Note this depends on the parallelization strategy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE PROC_FROM_POSITION(XP,YP,  IDPROC)

      ! Note: processes go from 0 (usually termed the Master) to MPI_N_THREADS - 1
      ! If the number of MPI processes is only 1, then IDPROC is 0, the one and only process.
      ! Otherwise, check the partition style (variable "DOMPART_TYPE")

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: XP, YP
      INTEGER, INTENT(OUT)     :: IDPROC

      INTEGER :: NCELLSPP
      INTEGER :: IDCELL

      IF (N_MPI_THREADS == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Serial operation

        IDPROC = 0      

      ELSE IF (DOMPART_TYPE == 0) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ S T R I P S
      !
      ! Domain partition type: strips
      !
      ! This domain partition assigns the same number of grid cells to every process.
      ! The domain partition is thus in strips along the "x" axis, as the cells are assumed 
      ! to be numbered along x.
 
         ! 1) Find ID of cell where the particle is
         CALL CELL_FROM_POSITION(XP, YP, IDCELL)

         ! 2) Find number of cells for each process (NX*NY*NZ = number of cells)
         !    Exceed a bit, so the last processor has slightly less cells, if number
         !    of cells is not divisible by the MPI_threads
         NCELLSPP = CEILING(REAL(NX*NY)/REAL(N_MPI_THREADS))

         ! 3) Here is the process ID 
         IDPROC   = INT(IDCELL/NCELLSPP)

      ELSE IF (DOMPART_TYPE == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ B L O C K S
      !
      ! Domain partition type: blocks
      !
      ! This domain partition style divides the domain into blocks
      ! Note that blocks in this definition are independent from cells! IT IS POSSIBLE TO
      ! GENERALIZE IT (and we could and should.. but for PIC I don't really care).

         IDPROC = INT((XP-XMIN)/DX_BLOCKS) + N_BLOCKS_X*INT((YP-YMIN)/DY_BLOCKS)

      ELSE ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ E R R O R

         WRITE(*,*) 'ATTENTION! Value for variable DOMPART_TYPE = ', DOMPART_TYPE
         WRITE(*,*) ' not recognized! Check input file!! In: PROC_FROM_POSITION() ABORTING!'
         STOP

      END IF

   END SUBROUTINE PROC_FROM_POSITION


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE COUNTING_SORT -> Sorts particles using the "counting sort"  !!
   ! sorting algorithm. Requires "PROC_FROM_POSITION" subroutine.           !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE COUNTING_SORT(PARTICLES_ARRAY, PARTICLES_COUNT)
 
   TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), INTENT(INOUT) :: PARTICLES_ARRAY ! array
   INTEGER                      , DIMENSION(:), INTENT(IN)    :: PARTICLES_COUNT ! freq
   INTEGER                      , DIMENSION(:), ALLOCATABLE   :: disp1
   TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE   :: sorted
   INTEGER                                                    :: i1
   INTEGER                                                    :: IPROC
 
   ALLOCATE(sorted(SIZE(PARTICLES_ARRAY)))
   ALLOCATE(disp1(0:SIZE(PARTICLES_COUNT)-1))
 
   disp1 = PARTICLES_COUNT
 
   DO i1 = 1, N_MPI_THREADS - 1
      disp1(i1) = disp1(i1) + disp1(i1-1)
   END DO
 
   DO i1 = SIZE(PARTICLES_ARRAY), 1, -1
 
      ! Find position of particle
      CALL PROC_FROM_POSITION(PARTICLES_ARRAY(i1)%X, PARTICLES_ARRAY(i1)%Y,  IPROC)
 
      sorted(disp1(IPROC)) = PARTICLES_ARRAY(i1)
      disp1(IPROC)         = disp1(IPROC) - 1
 
   END DO
 
   PARTICLES_ARRAY = sorted
 
   DEALLOCATE(sorted)
   DEALLOCATE(disp1)
 
   RETURN
 
   END SUBROUTINE COUNTING_SORT


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE EXCHANGE -> Exchanges particles between processes !!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE EXCHANGE

   ! This subroutine exchanges particles among processes.
   ! 1) It counts how many particles, in the vector of particles, are to be sent
   ! 2) Copies them in the sendbuf vector and removes them from the "particles" array
   ! 3) Sends around the information on how many they are
   ! 4) Sends them with MPI_ALLTOALLV
   ! 5) Copies them at the bottom of the "particles" array

      IMPLICIT NONE

      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: sendbuf, recvbuf

      INTEGER, DIMENSION(:), ALLOCATABLE :: sendcount, recvcount
      INTEGER, DIMENSION(:), ALLOCATABLE :: senddispl, recvdispl

      INTEGER :: i, IP, JP, JPROC, NP_RECV

      ALLOCATE(sendbuf(NP_PROC)) ! Allocate some space for the send buffer (max NP_PROC particles will be exchanged)

      ALLOCATE(sendcount(N_MPI_THREADS)) ! Allocate other vectors
      ALLOCATE(recvcount(N_MPI_THREADS))
      ALLOCATE(senddispl(N_MPI_THREADS))
      ALLOCATE(recvdispl(N_MPI_THREADS))

      DO i = 1,N_MPI_THREADS
         sendcount(i) = 0
         recvcount(i) = 0
         senddispl(i) = 0
         recvdispl(i) = 0
      END DO

      ! Loop on particles, from the last one to the first one and check if they belong to current process
      
      IP = NP_PROC ! Init 
      JP = 0       ! Init

      DO WHILE( IP .GE. 1 )

         CALL PROC_FROM_POSITION(particles(IP)%X, particles(IP)%Y, JPROC) ! Find which processor the particle belongs to

         IF (JPROC .NE. PROC_ID) THEN ! I shall send it to processor JPROC

            ! Increment the number of particles that shall be sent to processor JPROC 
            sendcount(JPROC+1) = sendcount(JPROC+1) + 1 ! Note: processors ID start from 0. Index ID from 1

            ! Copy particle in the send buffer
            JP = JP + 1
            sendbuf(JP) = particles(IP) 

            ! Remove particle from current array
            CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)

         END IF

         IP = IP - 1

      END DO

      CALL COUNTING_SORT(sendbuf(1:JP), sendcount) ! Reorder send buffer by process ID

      ! ~~~~~~ At this point, exchange particles among processes ~~~~~~

      ! Say how many particles are to be received
      CALL MPI_ALLTOALL(sendcount, 1, MPI_INTEGER, recvcount, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      NP_RECV = SUM(recvcount)
      ALLOCATE(recvbuf(NP_RECV))

      ! Compute position of particle chunks to be sent & received (incremental)
      DO i = 1, N_MPI_THREADS-1
         senddispl(i+1) = senddispl(i) + sendcount(i)
         recvdispl(i+1) = recvdispl(i) + recvcount(i)
      END DO

      ! Compute position of particle chunks to be sent & received (incremental), and send them 
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      CALL MPI_ALLTOALLV(sendbuf, sendcount, senddispl, MPI_PARTICLE_DATA_STRUCTURE, &
                         recvbuf, recvcount, recvdispl, MPI_PARTICLE_DATA_STRUCTURE, &
                         MPI_COMM_WORLD, ierr)

      ! At this point, write the particles received in the particles vector. Note that among these received,
      ! the processor also has its own particles.

      DO IP = 1, NP_RECV
         CALL ADD_PARTICLE_ARRAY(recvbuf(IP), NP_PROC, particles)
      END DO

      ! ~~~~~~~~ Done. Now deallocate stuff that would waste memory ~~~~~~~
      DEALLOCATE(sendbuf) 
      DEALLOCATE(recvbuf) 
      DEALLOCATE(sendcount)
      DEALLOCATE(recvcount)
      DEALLOCATE(senddispl)
      DEALLOCATE(recvdispl)
      
   END SUBROUTINE EXCHANGE





 
!!! ??????   !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ??????   !      ! SUBROUTINE PARTITION_FROM_CELL_ID -> finds the ID of a grid cell from particle position !!!
!!! ??????   !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ??????   !    
!!! ??????   !      SUBROUTINE PARTITION_FROM_CELL_ID(IDCELL, IDPART)
!!! ??????   !    
!!! ??????   !         IMPLICIT NONE
!!! ??????   !   
!!! ??????   !         INTEGER, INTENT(IN)  :: IDCELL ! ID of cell
!!! ??????   !         INTEGER, INTENT(OUT) :: IDPART ! ID of processor partition to which the cell belong
!!! ??????   !   
!!! ??????   !         INTEGER :: NCELLS, NCELL_PER_PART
!!! ??????   !   
!!! ??????   !         ! Check (but remove this in the future)
!!! ??????   !         IF(NCELLS > N_MPI_THREADS) THEN
!!! ??????   !            WRITE(*,*) 'ATTENTION! MORE MPI THREADS THAN PARTICLES. THIS IS NOT ALLOWED!'
!!! ??????   !            EXIT
!!! ??????   !         END IF
!!! ??????   !   
!!! ??????   !         NCELLS         = NX*NY ! Total number of cells
!!! ??????   !         NCELL_PER_PART = INT(NCELLS/N_MPI_THREADS)
!!! ??????   !   
!!! ??????   !         IDPART = INT(IDCELL/NCELLS) + 1
!!! ??????   !   
!!! ??????   !   !       INTEGER :: XID, YID
!!! ??????   !   ! 
!!! ??????   !   !       ! Partition by rows
!!! ??????   !   ! 
!!! ??????   !   !       ! Partition by columns
!!! ??????   !    
!!! ??????   !      END SUBROUTINE CELL_ID_FROM_POSITION


END MODULE grid_and_partition
