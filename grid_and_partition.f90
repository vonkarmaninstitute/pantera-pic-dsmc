MODULE grid_and_partition 

   USE global
   USE mpi_common
   USE screen
   USE tools

   IMPLICIT NONE
 
   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CELL_FROM_POSITION -> finds the ID of a grid cell from particle position      !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE CELL_FROM_POSITION(XP,YP,  IDCELL)

      ! Note: first cell index is 1.

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: XP, YP ! Location of particle
      INTEGER, INTENT(OUT)     :: IDCELL ! ID of cell to which the particle belongs

      INTEGER      :: XCELL, YCELL
      REAL(KIND=8) :: DX, DY

      IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
         ! Cartesian grid with equally spaced cells
         
         DX = (XMAX - XMIN)/NX
         DY = (YMAX - YMIN)/NY

         !WRITE(*,*) 'XP = ', XP, 'XMIN = ', XMIN, 'DX = ', DX
         XCELL = INT((XP-XMIN)/DX)
         YCELL = INT((YP-YMIN)/DY)

         ! IF (XCELL .GT. (NX-1)) THEN 
         !    XCELL = NX-1
         !    WRITE(*,*) 'Particle out of bound xhi!'
         ! ELSE IF (XCELL .LT. 0) THEN
         !    XCELL = 0
         !    WRITE(*,*) 'Particle out of bound xlo!'
         ! END IF

         ! IF (YCELL .GT. (NY-1)) THEN 
         !    YCELL = NY-1
         !    WRITE(*,*) 'Particle out of bound yhi!'
         ! ELSE IF (YCELL .LT. 0) THEN
         !    YCELL = 0
         !    WRITE(*,*) 'Particle out of bound ylo!'
         ! END IF

         IDCELL = XCELL + NX*YCELL + 1
      ELSE IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) THEN
         XCELL = BINARY_SEARCH(XP, XCOORD)
         YCELL = BINARY_SEARCH(YP, YCOORD)

         ! IF (XCELL .GT. (NX)) THEN 
         !    XCELL = NX
         !    WRITE(*,*) 'Particle out of bound xhi!', XCELL, NX
         ! ELSE IF (XCELL .LT. 1) THEN
         !    XCELL = 1
         !    WRITE(*,*) 'Particle out of bound xlo!'
         ! END IF

         ! IF (YCELL .GT. (NY)) THEN 
         !    YCELL = NY
         !    WRITE(*,*) 'Particle out of bound yhi!'
         ! ELSE IF (YCELL .LT. 1) THEN
         !    YCELL = 1
         !    WRITE(*,*) 'Particle out of bound ylo!'
         ! END IF

         IDCELL = XCELL + NX*(YCELL-1)
      END IF

   END SUBROUTINE CELL_FROM_POSITION

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE BINARY_SEARCH -> Finds the INDEX such that                           !
   ! ARRAY(INDEX) < VALUE < ARRAY(INDEX+1), where ARRAY is monotonically increasing. !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION BINARY_SEARCH(VALUE, ARRAY) RESULT(INDEX)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: VALUE
      INTEGER :: L, R
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ARRAY
      INTEGER :: INDEX

      L = LBOUND(ARRAY, DIM=1)
      R = UBOUND(ARRAY, DIM=1)

      INDEX = -1
      IF (VALUE .LT. ARRAY(L) .OR. VALUE .GT. ARRAY(R)) THEN
         WRITE(*,*) 'Particle out of bounds!'
         RETURN
      ELSE IF (R == L+1) THEN
         INDEX = L
         RETURN
      ELSE
         DO
            INDEX = (L+R)/2
            IF (ARRAY(INDEX) .LE. VALUE) THEN
               IF (ARRAY(INDEX+1) .GT. VALUE) RETURN
               L = INDEX
            ELSE
               IF (ARRAY(INDEX-1) .LE. VALUE) THEN
                  INDEX = INDEX-1
                  RETURN
               END IF
               R = INDEX
            END IF
         END DO
         RETURN
      END IF

   END FUNCTION BINARY_SEARCH



   SUBROUTINE PROC_FROM_CELL(IDCELL, IDPROC)

      ! Note: processes go from 0 (usually termed the Master) to MPI_N_THREADS - 1
      ! If the number of MPI processes is only 1, then IDPROC is 0, the one and only process.
      ! Otherwise, check the partition style (variable "DOMPART_TYPE")

      IMPLICIT NONE

      INTEGER, INTENT(IN)     :: IDCELL
      INTEGER, INTENT(OUT)     :: IDPROC
      
      IF (N_MPI_THREADS == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Serial operation

         IDPROC = 0
   
      ELSE

         IDPROC = CELL_PROCS(IDCELL)

      END IF

   END SUBROUTINE PROC_FROM_CELL


   SUBROUTINE ASSIGN_CELLS_TO_PROCS


      INTEGER :: IDCELL
      INTEGER :: NCELLSPP
      INTEGER :: I, J, IPROC, IP
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CENTROID, ARRAY_TEST
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: WEIGHT
      INTEGER, DIMENSION(:), ALLOCATABLE :: NP_CELLS
      REAL(KIND=8) :: COORDMAX, COORDMIN, WEIGHT_PER_PROC, CUMULATIVE_WEIGHT
      INTEGER, DIMENSION(:), ALLOCATABLE :: ORDER, ORDER_TEST

      ! Weight is a generic measure that is ideally directly proportional to the computational cost for each cell.

      IF (.NOT. ALLOCATED(CELL_PROCS)) ALLOCATE(CELL_PROCS(NCELLS))

      IF (N_MPI_THREADS == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Serial operation

         CELL_PROCS = 0
      
      ELSE IF (GRID_TYPE == UNSTRUCTURED) THEN
         IF (.FALSE. .OR. DIMS .NE. 2) THEN

            ALLOCATE(CENTROID(NCELLS))

            IF (DIMS == 3) THEN

               DO I = 1, NCELLS
                  CENTROID(I) = (U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,1), 2) &
                              +  U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,2), 2) &
                              +  U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,3), 2) &
                              +  U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,4), 2)) / 4.
               END DO
         
               COORDMAX = MAXVAL(U3D_GRID%NODE_COORDS(:,2))
               COORDMIN = MINVAL(U3D_GRID%NODE_COORDS(:,2))
               
            ELSE IF (DIMS == 2) THEN

               DO I = 1, NCELLS
                  CENTROID(I) = (U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), 2) &
                              +  U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,2), 2) &
                              +  U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,3), 2)) / 3.
               END DO
         
               COORDMAX = MAXVAL(U2D_GRID%NODE_COORDS(:,2))
               COORDMIN = MINVAL(U2D_GRID%NODE_COORDS(:,2))
         
            END IF   

            DO I = 1, NCELLS
               CELL_PROCS(I) = INT((CENTROID(I)-COORDMIN)/(COORDMAX-COORDMIN)*REAL(N_MPI_THREADS))
            END DO

            DEALLOCATE(CENTROID)

         ELSE

            ALLOCATE(CENTROID(NCELLS))
            ALLOCATE(WEIGHT(NCELLS))
            ALLOCATE(ORDER(NCELLS))
            ALLOCATE(NP_CELLS(NCELLS))
            ! Assign weight
            NP_CELLS = 1
            DO IP = 1, NP_PROC
               NP_CELLS(particles(IP)%IC) = NP_CELLS(particles(IP)%IC) + 1
            END DO

            IF (PROC_ID .EQ. 0) THEN
               CALL MPI_REDUCE(MPI_IN_PLACE, NP_CELLS, NCELLS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            ELSE
               CALL MPI_REDUCE(NP_CELLS,     NP_CELLS, NCELLS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            END IF
      
            CALL MPI_BCAST(NP_CELLS, NCELLS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            WEIGHT = DBLE(NP_CELLS) + SUM(NP_CELLS)/10.
            
            WEIGHT_PER_PROC = SUM(WEIGHT) / DBLE(N_MPI_THREADS)

            IF (DIMS == 2) THEN

               DO I = 1, NCELLS
                  CENTROID(I) = (U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), 2) &
                              +  U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,2), 2) &
                              +  U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,3), 2)) / 3.
                  ORDER(I) = I
               END DO
         
            END IF

            WRITE(*,*) 'Started sorting'
            IF (PROC_ID == 0) THEN
               ALLOCATE(ARRAY_TEST(10))
               ALLOCATE(ORDER_TEST(10))
               ARRAY_TEST = [2.3, 4.4, 0.2, 1.0, 1.12, 7.9, 101.1, 2.1, 67.5, 11.1]
               WRITE(*,*) ARRAY_TEST
               ORDER_TEST = [1,2,3,4,5,6,7,8,9,10]
               CALL QUICKSORT(ARRAY_TEST, ORDER_TEST, 10)
               WRITE(*,*) ARRAY_TEST
               WRITE(*,*) ORDER_TEST
            END IF

            !CALL QUICKSORT(CENTROID, ORDER, NCELLS)
            WRITE(*,*) 'Ended sorting'

            IPROC = 0
            DO I = 1, NCELLS
               CUMULATIVE_WEIGHT = CUMULATIVE_WEIGHT + WEIGHT(I)
               CELL_PROCS(I) = IPROC

               IF (CUMULATIVE_WEIGHT > WEIGHT_PER_PROC) THEN
                  IPROC = IPROC + 1
                  CUMULATIVE_WEIGHT = 0
               END IF

            END DO

            DEALLOCATE(CENTROID)
            DEALLOCATE(WEIGHT)
            DEALLOCATE(ORDER)
            DEALLOCATE(NP_CELLS)

         END IF

      ELSE IF (DOMPART_TYPE == 0) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ S T R I P S
      !
      ! Domain partition type: strips
      !
      ! This domain partition assigns the same number of grid cells to every process.
      ! The domain partition is thus in strips along the "x" axis, as the cells are assumed 
      ! to be numbered along x.

         NCELLSPP = CEILING(REAL(NCELLS)/REAL(N_MPI_THREADS))

         DO IDCELL = 1, NCELLS
            ! 3) Here is the process ID 
            CELL_PROCS(IDCELL) = INT((IDCELL-1)/NCELLSPP) ! Before was giving wrong result for peculiar combinations of NCELLS and NCELLSPP
         END DO

      ELSE IF (DOMPART_TYPE == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ B L O C K S
      !
      ! Domain partition type: blocks
      !
      ! This domain partition style divides the domain into blocks
      ! Note that blocks in this definition are independent from cells! IT IS POSSIBLE TO
      ! GENERALIZE IT.

         DO IDCELL = 1, NCELLS
            I = MOD(IDCELL - 1, NX)
            J = (IDCELL - 1)/NX
            CELL_PROCS(IDCELL) = I*N_BLOCKS_X/NX + N_BLOCKS_X*(J*N_BLOCKS_Y/NY)
         END DO

      ELSE ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ E R R O R

         WRITE(*,*) 'ATTENTION! Value for variable DOMPART_TYPE = ', DOMPART_TYPE
         WRITE(*,*) ' not recognized! Check input file!! In: PROC_FROM_POSITION() ABORTING!'
         STOP
      END IF

   END SUBROUTINE ASSIGN_CELLS_TO_PROCS

   ! SUBROUTINE ASSIGN_CELLS_TO_PROCS

   !    INTEGER, DIMENSION(:), ALLOCATABLE :: NPC, NPCMOD, NPP
   !    INTEGER :: JP, IC, NPPDESIRED, IPROC
      
   !    ALLOCATE(NPC(NCELLS))
   !    NPC = 0

   !    DO JP = 1, NP_PROC
   !       IC = particles(JP)%IC
   !       NPC(IC) = NPC(IC) + 1
   !    END DO

   !    IF (PROC_ID .EQ. 0) THEN
   !       CALL MPI_REDUCE(MPI_IN_PLACE, NPC, NCELLS, MPI_INTEGER,  MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   !    ELSE
   !       CALL MPI_REDUCE(NPC,          NPC, NCELLS, MPI_INTEGER,  MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   !    END IF

   !    CALL MPI_BCAST(NPC, NCELLS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


   !    NPCMOD = NPC + SUM(NPC)/NCELLS/10 + 1

   !    NPPDESIRED = SUM(NPCMOD) / N_MPI_THREADS

   !    ALLOCATE(NPP(N_MPI_THREADS))
   !    NPP = 0

   !    IF (.NOT. ALLOCATED(CELL_PROCS)) ALLOCATE(CELL_PROCS(NCELLS))
   !    CELL_PROCS = 0
   !    IPROC = 0
   !    DO IC = 1, NCELLS
   !       IF (NPP(IPROC) < NPPDESIRED .OR. IPROC == N_MPI_THREADS - 1) THEN
   !          NPP(IPROC) = NPP(IPROC) + NPCMOD(IC)
   !          CELL_PROCS(IC) = IPROC
   !       ELSE
   !          IPROC = IPROC + 1
   !          NPP(IPROC) = NPP(IPROC) + NPCMOD(IC)
   !          CELL_PROCS(IC) = IPROC
   !       END IF
   !    END DO

   !    DEALLOCATE(NPC)
   !    DEALLOCATE(NPP)


   ! END SUBROUTINE ASSIGN_CELLS_TO_PROCS



   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ! SUBROUTINE PROC_FROM_POSITION -> finds the process ID from particle position              !
   ! ! Note this depends on the parallelization strategy                                         !
   ! ! This should never be called when an unstructured grid is used.                            !
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! SUBROUTINE PROC_FROM_POSITION(XP,YP,  IDPROC)

   !    ! Note: processes go from 0 (usually termed the Master) to MPI_N_THREADS - 1
   !    ! If the number of MPI processes is only 1, then IDPROC is 0, the one and only process.
   !    ! Otherwise, check the partition style (variable "DOMPART_TYPE")

   !    IMPLICIT NONE

   !    REAL(KIND=8), INTENT(IN) :: XP, YP
   !    INTEGER, INTENT(OUT)     :: IDPROC

   !    INTEGER :: NCELLSPP
   !    INTEGER :: IDCELL

   !    IF (N_MPI_THREADS == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Serial operation

   !      IDPROC = 0      

   !    ELSE IF (DOMPART_TYPE == 0) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ S T R I P S
   !    !
   !    ! Domain partition type: strips
   !    !
   !    ! This domain partition assigns the same number of grid cells to every process.
   !    ! The domain partition is thus in strips along the "x" axis, as the cells are assumed 
   !    ! to be numbered along x.
 
   !       ! 1) Find ID of cell where the particle is
   !       CALL CELL_FROM_POSITION(XP, YP, IDCELL)
   !       IF (IDCELL .GT. NX*NY .OR. IDCELL .LT. 1) THEN
   !          WRITE(*,*) 'Error! CELL_FROM_POSITION returned cell:', IDCELL, 'Particle position: ', XP, ', ', YP
   !       END IF

   !       ! 2) Find number of cells for each process (NX*NY*NZ = number of cells)
   !       !    Exceed a bit, so the last processor has slightly less cells, if number
   !       !    of cells is not divisible by the MPI_threads
   !       NCELLSPP = CEILING(REAL(NX*NY)/REAL(N_MPI_THREADS))

   !       ! 3) Here is the process ID 
   !       IDPROC   = INT((IDCELL-1)/NCELLSPP) ! Before was giving wrong result for peculiar combinations of NCELLS and NCELLSPP
   !       IF (IDPROC .GT. N_MPI_THREADS-1 .OR. IDPROC .LT. 0) THEN
   !          WRITE(*,*) 'Error! PROC_FROM_POSITION returned proc:', IDPROC
   !       END IF


   !    ELSE IF (DOMPART_TYPE == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ B L O C K S
   !    !
   !    ! Domain partition type: blocks
   !    !
   !    ! This domain partition style divides the domain into blocks
   !    ! Note that blocks in this definition are independent from cells! IT IS POSSIBLE TO
   !    ! GENERALIZE IT (and we could and should.. but for PIC I don't really care).

   !       IDPROC = INT((XP-XMIN)/DX_BLOCKS) + N_BLOCKS_X*INT((YP-YMIN)/DY_BLOCKS)

   !    ELSE ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ E R R O R

   !       WRITE(*,*) 'ATTENTION! Value for variable DOMPART_TYPE = ', DOMPART_TYPE
   !       WRITE(*,*) ' not recognized! Check input file!! In: PROC_FROM_POSITION() ABORTING!'
   !       STOP

   !    END IF

   ! END SUBROUTINE PROC_FROM_POSITION


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
         !CALL PROC_FROM_POSITION(PARTICLES_ARRAY(i1)%X, PARTICLES_ARRAY(i1)%Y,  IPROC)
         CALL PROC_FROM_CELL(particles_array(i1)%IC, IPROC)
   
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

         !CALL PROC_FROM_POSITION(particles(IP)%X, particles(IP)%Y, JPROC) ! Find which processor the particle belongs to
         CALL PROC_FROM_CELL(particles(IP)%IC, JPROC)

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


   SUBROUTINE READ_UNSTRUCTURED_GRID(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE
      
      INTEGER, PARAMETER :: in5 = 2384
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER            :: NUM, I, J, FOUND, V1, V2, INDEX, ELEM_TYPE, NUM_GROUPS
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_LINE_NODES
      INTEGER, DIMENSION(:), ALLOCATABLE        :: TEMP_POINT_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS
      INTEGER :: NUM_CELLS, NUM_LINES, NUM_POINTS

      INTEGER, DIMENSION(2) :: WHICH1, WHICH2

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: STRARRAY

      INTEGER :: NCELLSPP

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. Gmsh v2 file format (*.msh)
      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Reading grid file.'
      WRITE(*,*) '==========================================='
      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE
         WRITE(*,*) 'Read line:', LINE
         IF (ReasonEOF < 0) EXIT 
         IF (LINE == '$Nodes') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM
            ALLOCATE(U2D_GRID%NODE_COORDS(NUM,3))
            ALLOCATE(U2D_GRID%CELL_NODES(NUM,3))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) INDEX, XYZ
               XYZ(3) = 0.d0 ! Stay in the x-y plane
               U2D_GRID%NODE_COORDS(INDEX,:) = XYZ
            END DO
            U2D_GRID%NUM_NODES = NUM
         ELSE IF (LINE == '$Elements') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM
            ALLOCATE(TEMP_POINT_NODES(NUM))
            ALLOCATE(TEMP_LINE_NODES(NUM,2))
            ALLOCATE(TEMP_CELL_NODES(NUM,3))
            NUM_CELLS = 0
            NUM_LINES = 0
            NUM_POINTS = 0
            DO I = 1, NUM
               READ(in5,'(A)', IOSTAT=ReasonEOF) LINE
               CALL SPLIT_STR(LINE, ' ', STRARRAY, N_STR)

               READ(STRARRAY(1),*) INDEX
               READ(STRARRAY(2),*) ELEM_TYPE
               READ(STRARRAY(3),*) NUM_GROUPS

               IF (ELEM_TYPE == 15) THEN
                  NUM_POINTS = NUM_POINTS + 1
                  READ(STRARRAY(4+NUM_GROUPS),*) TEMP_POINT_NODES(NUM_POINTS)
               ELSE IF (ELEM_TYPE == 1) THEN
                  NUM_LINES = NUM_LINES + 1
                  READ(STRARRAY(4+NUM_GROUPS),*) TEMP_LINE_NODES(NUM_LINES,1)
                  READ(STRARRAY(5+NUM_GROUPS),*) TEMP_LINE_NODES(NUM_LINES,2)
               ELSE IF (ELEM_TYPE == 2) THEN
                  NUM_CELLS = NUM_CELLS + 1
                  READ(STRARRAY(4+NUM_GROUPS),*) TEMP_CELL_NODES(NUM_CELLS,1)
                  READ(STRARRAY(5+NUM_GROUPS),*) TEMP_CELL_NODES(NUM_CELLS,2)
                  READ(STRARRAY(6+NUM_GROUPS),*) TEMP_CELL_NODES(NUM_CELLS,3)
               END IF
               DEALLOCATE(STRARRAY)
            END DO
            U2D_GRID%NUM_POINTS = NUM_POINTS
            U2D_GRID%NUM_LINES = NUM_LINES
            U2D_GRID%NUM_CELLS = NUM_CELLS

            U2D_GRID%POINT_NODES = TEMP_POINT_NODES(1:NUM_POINTS)
            U2D_GRID%LINE_NODES = TEMP_LINE_NODES(1:NUM_LINES,:)
            U2D_GRID%CELL_NODES = TEMP_CELL_NODES(1:NUM_CELLS,:)
         END IF
      END DO

      ! Done reading
      CLOSE(in5)



      WRITE(*,*) 'Read grid file. It contains ', U2D_GRID%NUM_POINTS, &
                 'points, ', U2D_GRID%NUM_LINES, 'lines, and ', &
                 U2D_GRID%NUM_CELLS, 'cells.'

      ! Process the mesh: generate connectivity, normals and such...
      !XMIN, XMAX,...

      DO I = 1, NUM_CELLS
         WRITE(*,*) U2D_GRID%CELL_NODES(I,:)
      END DO

      ! Compute cell volumes
      ALLOCATE(CELL_VOLUMES(NUM_CELLS))
      DO I = 1, NUM_CELLS
         A = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), :)
         B = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,2), :)
         C = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,3), :)

         IF (DIMS == 2 .AND. .NOT. AXI) THEN
            CELL_VOLUMES(I) = 0.5*ABS(A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))
            WRITE(*,*) CELL_VOLUMES(I)
         END IF
      END DO

      ! Find cell connectivity
      ALLOCATE(TEMP_CELL_NEIGHBORS(NUM_CELLS, 3))
      TEMP_CELL_NEIGHBORS = -1
      DO I = 1, NUM_CELLS
         DO J = I, NUM_CELLS
            IF (I == J) CYCLE
            FOUND = 0
            DO V1 = 1, 3
               DO V2 = 1, 3
                  IF (U2D_GRID%CELL_NODES(I,V1) == U2D_GRID%CELL_NODES(J,V2)) THEN
                     FOUND = FOUND + 1
                     IF (FOUND .GT. 2) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                     WHICH1(FOUND) = V1
                     WHICH2(FOUND) = V2
                  END IF
               END DO
            END DO
            IF (FOUND == 2) THEN
               IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) TEMP_CELL_NEIGHBORS(I, 1) = J
               IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) TEMP_CELL_NEIGHBORS(I, 2) = J
               IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) TEMP_CELL_NEIGHBORS(I, 3) = J
               IF (ANY(WHICH2 == 1) .AND. ANY(WHICH2 == 2)) TEMP_CELL_NEIGHBORS(J, 1) = I
               IF (ANY(WHICH2 == 2) .AND. ANY(WHICH2 == 3)) TEMP_CELL_NEIGHBORS(J, 2) = I
               IF (ANY(WHICH2 == 3) .AND. ANY(WHICH2 == 1)) TEMP_CELL_NEIGHBORS(J, 3) = I
            END IF
         END DO
      END DO

      U2D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS


      WRITE(*,*) 'Generated grid connectivity. '
      DO I = 1, NUM_CELLS
         WRITE(*,*) 'Cell ', I, ' neighbors cells ', TEMP_CELL_NEIGHBORS(I, :)
      END DO

      ! Distributing cells over MPI processes
      NCELLSPP = CEILING(REAL(NUM_CELLS)/REAL(N_MPI_THREADS))
      ALLOCATE(CELL_PROCS(NUM_CELLS))
      DO I = 1, NUM_CELLS
         CELL_PROCS(I) = INT((I-1)/NCELLSPP)
      END DO

      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Done reading grid file.'
      WRITE(*,*) '==========================================='

   END SUBROUTINE READ_UNSTRUCTURED_GRID





   SUBROUTINE READ_2D_UNSTRUCTURED_GRID_SU2(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER            :: NUM, I, J, FOUND, V1, V2, V3, ELEM_TYPE, NUMELEMS
      INTEGER, DIMENSION(3,2) :: IND
      REAL(KIND=8)       :: X1, X2, X3, Y1, Y2, Y3, LEN, RAD
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(2) :: VLIST, WHICH1, WHICH2

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. SU2 file format (*.su2)
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Reading grid file in SU2 format.'
         WRITE(*,*) '==========================================='
      END IF
      
      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            ALLOCATE(U2D_GRID%NODE_COORDS(NUM,3))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) XYZ(1:2)
               XYZ(3) = 0.d0 ! Stay in the x-y plane.
               U2D_GRID%NODE_COORDS(I,:) = XYZ
            END DO
            U2D_GRID%NUM_NODES = NUM
         ELSE IF (LINE == 'NELEM=') THEN
            ALLOCATE(U2D_GRID%CELL_NODES(NUM,3))

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, U2D_GRID%CELL_NODES(I,:)
               !WRITE(*,*) 'I read element ', I, ' has nodes ', U2D_GRID%CELL_NODES(I,:)
               IF (ELEM_TYPE .NE. 5) WRITE(*,*) 'Element type was not 5!'
            END DO
            U2D_GRID%CELL_NODES = U2D_GRID%CELL_NODES + 1 ! Start indexing from 1.

            U2D_GRID%NUM_CELLS = NUM

            ALLOCATE(U2D_GRID%CELL_EDGES_PG(U2D_GRID%NUM_CELLS, 3))
            U2D_GRID%CELL_EDGES_PG = -1
      
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cell edges.
            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
               IF (LINE .NE. 'MARKER_TAG=') THEN
                  WRITE(*,*) 'Error! did not find marker name.'
               ELSE
                  WRITE(*,*) 'Found marker tag, with groupname: ', GROUPNAME
               END IF
         
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS
               IF (LINE .NE. 'MARKER_ELEMS=') THEN
                  WRITE(*,*) 'Error! did not find marker elements.'
               ELSE
                  WRITE(*,*) 'Found marker elements, with number of elements: ', NUMELEMS
               END IF

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
               END DO
            END DO

         END IF
      END DO

      REWIND(in5)

      ALLOCATE(N_CELLS_WITH_NODE(U2D_GRID%NUM_NODES))
      ALLOCATE(IOF(U2D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U2D_GRID%NUM_CELLS
         DO V1 = 1, 3
            JN = U2D_GRID%CELL_NODES(I,V1)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U2D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U2D_GRID%NUM_CELLS
         DO V1 = 1, 3
            JN = U2D_GRID%CELL_NODES(I,V1)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NELEM=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cell edges.

            ALLOCATE(GRID_BC(NUM)) ! Append the physical group to the list
            N_GRID_BC = NUM

            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
      
               GRID_BC(I)%PHYSICAL_GROUP_NAME = GROUPNAME
         
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, VLIST
                  IF (ELEM_TYPE .NE. 3) WRITE(*,*) 'Error! element type was not line.'
                  VLIST = VLIST + 1

                  JN = VLIST(1)
                  IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                     DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                        JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                        FOUND = 0
                        DO V1 = 1, 3
                           IF (ANY(VLIST == U2D_GRID%CELL_NODES(JC1,V1))) THEN
                              FOUND = FOUND + 1
                              WHICH1(FOUND) = V1
                           END IF
                        END DO
         
                        IF (FOUND == 2) THEN
                           IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) THEN
                              U2D_GRID%CELL_EDGES_PG(JC1, 1) = I
                           ELSE IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) THEN
                              U2D_GRID%CELL_EDGES_PG(JC1, 2) = I
                           ELSE IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) THEN
                              U2D_GRID%CELL_EDGES_PG(JC1, 3) = I
                           END IF
                        END IF
                     END DO
                  END IF

               END DO
            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)

      !WRITE(*,*) 'Read grid file. It contains ', U2D_GRID%NUM_NODES, &
      !           'points, and ', U2D_GRID%NUM_CELLS, 'cells.'

      ! Process the mesh: generate connectivity, normals and such...
      !XMIN, XMAX,...

      !DO I = 1, U2D_GRID%NUM_CELLS
      !   WRITE(*,*) U2D_GRID%CELL_NODES(I,:)
      !END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing cell volumes.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell volumes
      ALLOCATE(CELL_AREAS(U2D_GRID%NUM_CELLS))
      ALLOCATE(CELL_VOLUMES(U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         A = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), :)
         B = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,2), :)
         C = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,3), :)

         CELL_AREAS(I) = 0.5*ABS(A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))
         IF (DIMS == 2 .AND. .NOT. AXI) THEN
            CELL_VOLUMES(I) = CELL_AREAS(I) * (ZMAX-ZMIN)
            !WRITE(*,*) CELL_VOLUMES(I)
         END IF
         IF (DIMS == 2 .AND. AXI) THEN
            RAD = (A(2)+B(2)+C(2))/3.
            CELL_VOLUMES(I) = CELL_AREAS(I) * (ZMAX-ZMIN)*RAD
         END IF
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing grid connectivity.'
         WRITE(*,*) '==========================================='
      END IF

      ! Find cell connectivity
      ALLOCATE(TEMP_CELL_NEIGHBORS(U2D_GRID%NUM_CELLS, 3))
      TEMP_CELL_NEIGHBORS = -1



      DO JN = 1, U2D_GRID%NUM_NODES
         !IF (PROC_ID == 0) WRITE(*,*) 'Checking node ', JN, ' of ',  U2D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)


                  FOUND = 0
                  DO V1 = 1, 3
                     DO V2 = 1, 3
                        IF (U2D_GRID%CELL_NODES(JC1,V1) == U2D_GRID%CELL_NODES(JC2,V2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 2) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           WHICH1(FOUND) = V1
                           WHICH2(FOUND) = V2
                        END IF
                     END DO
                  END DO

                  IF (FOUND == 2) THEN
      
                     IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) THEN
                        TEMP_CELL_NEIGHBORS(JC1, 1) = JC2
                     ELSE IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) THEN
                        TEMP_CELL_NEIGHBORS(JC1, 2) = JC2
                     ELSE IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) THEN
                        TEMP_CELL_NEIGHBORS(JC1, 3) = JC2
                     END IF

                     IF (ANY(WHICH2 == 1) .AND. ANY(WHICH2 == 2)) THEN
                        TEMP_CELL_NEIGHBORS(JC2, 1) = JC1
                     ELSE IF (ANY(WHICH2 == 2) .AND. ANY(WHICH2 == 3)) THEN
                        TEMP_CELL_NEIGHBORS(JC2, 2) = JC1
                     ELSE IF (ANY(WHICH2 == 3) .AND. ANY(WHICH2 == 1)) THEN
                        TEMP_CELL_NEIGHBORS(JC2, 3) = JC1
                     END IF
      
                  END IF


               END DO
            END DO
         END IF
      END DO

      U2D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS



      !WRITE(*,*) 'Generated grid connectivity. '
      !DO I = 1, U2D_GRID%NUM_CELLS
      !   WRITE(*,*) 'Cell ', I, ' neighbors cells ', TEMP_CELL_NEIGHBORS(I, :)
      !END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing face normals.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell edge normals
      IND(1,:) = [1,2]
      IND(2,:) = [2,3]
      IND(3,:) = [3,1]
      ALLOCATE(U2D_GRID%EDGE_NORMAL(U2D_GRID%NUM_CELLS, 3, 3))
      ALLOCATE(U2D_GRID%CELL_EDGES_LEN(U2D_GRID%NUM_CELLS, 3))
      DO I = 1, U2D_GRID%NUM_CELLS
         DO J = 1, 3
            X1 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,IND(J,1)), 1)
            X2 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,IND(J,2)), 1)
            Y1 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,IND(J,1)), 2)
            Y2 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,IND(J,2)), 2)
            LEN = SQRT((Y2-Y1)*(Y2-Y1) + (X2-X1)*(X2-X1))
            U2D_GRID%CELL_EDGES_LEN(I,J) = LEN
            U2D_GRID%EDGE_NORMAL(I,J,1) = (Y2-Y1)/LEN
            U2D_GRID%EDGE_NORMAL(I,J,2) = (X1-X2)/LEN
            U2D_GRID%EDGE_NORMAL(I,J,3) = 0.d0
         END DO
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Checking ordering.'
         WRITE(*,*) '==========================================='
      END IF

      DO I = 1, U2D_GRID%NUM_CELLS
         X1 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,2), 1) &
            - U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), 1)
         X2 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,3), 1) &
            - U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), 1)
         Y1 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,2), 2) &
            - U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), 2)
         Y2 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,3), 2) &
            - U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), 2)

         IF (X1*Y2-X2*Y1 < 0) CALL ERROR_ABORT('2D mesh triangles have negative z-normal.')

      END DO

      NCELLS = U2D_GRID%NUM_CELLS
      NNODES = U2D_GRID%NUM_NODES



      ALLOCATE(U2D_GRID%BASIS_COEFFS(NCELLS,3,3))

      DO I = 1, NCELLS
         V1 = U2D_GRID%CELL_NODES(I,1)
         V2 = U2D_GRID%CELL_NODES(I,2)
         V3 = U2D_GRID%CELL_NODES(I,3)

         X1 = U2D_GRID%NODE_COORDS(V1, 1)
         X2 = U2D_GRID%NODE_COORDS(V2, 1)
         X3 = U2D_GRID%NODE_COORDS(V3, 1)
         Y1 = U2D_GRID%NODE_COORDS(V1, 2)
         Y2 = U2D_GRID%NODE_COORDS(V2, 2)
         Y3 = U2D_GRID%NODE_COORDS(V3, 2)

         ! These are such that PSI_i = SUM_j [ x_j * BASIS_COEFFS(IC,i,j) ] + BASIS_COEFFS(IC,i,3)

         U2D_GRID%BASIS_COEFFS(I,1,1) =  Y2-Y3
         U2D_GRID%BASIS_COEFFS(I,1,2) = -(X2-X3)
         U2D_GRID%BASIS_COEFFS(I,1,3) =  X2*Y3 - X3*Y2

         U2D_GRID%BASIS_COEFFS(I,2,1) = -(Y1-Y3)
         U2D_GRID%BASIS_COEFFS(I,2,2) =  X1-X3
         U2D_GRID%BASIS_COEFFS(I,2,3) =  X3*Y1 - X1*Y3

         U2D_GRID%BASIS_COEFFS(I,3,1) = -(Y2-Y1)
         U2D_GRID%BASIS_COEFFS(I,3,2) =  X2-X1
         U2D_GRID%BASIS_COEFFS(I,3,3) =  X1*Y2 - X2*Y1

         U2D_GRID%BASIS_COEFFS(I,:,:) = 0.5*U2D_GRID%BASIS_COEFFS(I,:,:)/CELL_AREAS(I)

      END DO


      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Done reading grid file.'
         WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
         WRITE(*,*) '==========================================='
      END IF

   END SUBROUTINE READ_2D_UNSTRUCTURED_GRID_SU2








   SUBROUTINE READ_3D_UNSTRUCTURED_GRID_SU2(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER            :: NUM, I, J, FOUND, V1, V2, V3, V4, ELEM_TYPE, NUMELEMS
      INTEGER, DIMENSION(4,3) :: IND
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C, CROSSP

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(3) :: VLIST, WHICH1, WHICH2

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2

      REAL(KIND=8) :: VOLUME, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. SU2 file format (*.su2)
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Reading grid file in SU2 format.'
         WRITE(*,*) '==========================================='
      END IF
      
      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            ALLOCATE(U3D_GRID%NODE_COORDS(NUM,3))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) XYZ
               U3D_GRID%NODE_COORDS(I,:) = XYZ
            END DO
            U3D_GRID%NUM_NODES = NUM
         ELSE IF (LINE == 'NELEM=') THEN
            ALLOCATE(U3D_GRID%CELL_NODES(NUM,4))

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, U3D_GRID%CELL_NODES(I,:)
               !WRITE(*,*) 'I read element ', I, ' has nodes ', U3D_GRID%CELL_NODES(I,:)
               IF (ELEM_TYPE .NE. 10) CALL ERROR_ABORT('Reading 3D grid found element type was not tetrahedron (type 10).')
            END DO
            U3D_GRID%CELL_NODES = U3D_GRID%CELL_NODES + 1 ! Start indexing from 1.

            U3D_GRID%NUM_CELLS = NUM

            ALLOCATE(U3D_GRID%CELL_FACES_PG(U3D_GRID%NUM_CELLS, 4))
            U3D_GRID%CELL_FACES_PG = -1
         ELSE IF (LINE == 'NMARK=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
               IF (LINE .NE. 'MARKER_TAG=') THEN
                  CALL ERROR_ABORT('Error! did not find marker name.')
               ELSE
                  IF (PROC_ID == 0) WRITE(*,*) 'Found marker tag, with groupname: ', GROUPNAME
               END IF

               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS
               IF (LINE .NE. 'MARKER_ELEMS=') THEN
                  CALL ERROR_ABORT('Error! did not find marker elements.')
               ELSE
                  IF (PROC_ID == 0) WRITE(*,*) 'Found marker elements, with number of elements: ', NUMELEMS
               END IF

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
               END DO
            END DO
         END IF
      END DO

      REWIND(in5)



      ALLOCATE(N_CELLS_WITH_NODE(U3D_GRID%NUM_NODES))
      ALLOCATE(IOF(U3D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U3D_GRID%NUM_CELLS
         DO V1 = 1, 4
            JN = U3D_GRID%CELL_NODES(I,V1)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U3D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U3D_GRID%NUM_CELLS
         DO V1 = 1, 4
            JN = U3D_GRID%CELL_NODES(I,V1)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NELEM=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cell edges.

            ALLOCATE(GRID_BC(NUM)) ! Append the physical group to the list
            N_GRID_BC = NUM

            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
      
               GRID_BC(I)%PHYSICAL_GROUP_NAME = GROUPNAME
         
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, VLIST
                  IF (ELEM_TYPE .NE. 5) CALL ERROR_ABORT('Error! element type was not triangle.')
                  VLIST = VLIST + 1

                  JN = VLIST(1)
                  IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                     DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                        JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                        FOUND = 0
                        DO V1 = 1, 4
                           IF (ANY(VLIST == U3D_GRID%CELL_NODES(JC1,V1))) THEN
                              FOUND = FOUND + 1
                              WHICH1(FOUND) = V1
                           END IF
                        END DO
         
                        IF (FOUND == 3) THEN
            
                           IF (ANY(WHICH1 == 1)) THEN
                              IF (ANY(WHICH1 == 2)) THEN
                                 IF (ANY(WHICH1 == 3))  THEN
                                    U3D_GRID%CELL_FACES_PG(JC1, 1) = I
                                 ELSE IF (ANY(WHICH1 == 4)) THEN
                                    U3D_GRID%CELL_FACES_PG(JC1, 2) = I
                                 END IF
                              ELSE IF (ANY(WHICH1 == 3)) THEN
                                 IF (ANY(WHICH1 == 4)) U3D_GRID%CELL_FACES_PG(JC1, 4) = I
                              END IF
                           ELSE IF (ANY(WHICH1 == 2)) THEN
                              IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 4)) U3D_GRID%CELL_FACES_PG(JC1, 3) = I
                           END IF

                        END IF
                     END DO
                  END IF

               END DO
            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)

      !WRITE(*,*) 'Read grid file. It contains ', U3D_GRID%NUM_NODES, &
      !           'points, and ', U3D_GRID%NUM_CELLS, 'cells.'

      ! Process the mesh: generate connectivity, normals and such...
      !XMIN, XMAX,...

      !DO I = 1, U3D_GRID%NUM_CELLS
      !   WRITE(*,*) U3D_GRID%CELL_NODES(I,:)
      !END DO
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing cell volumes.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell volumes
      ALLOCATE(CELL_VOLUMES(U3D_GRID%NUM_CELLS))
      DO I = 1, U3D_GRID%NUM_CELLS
         A = U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,2), :) - U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,1), :)
         B = U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,3), :) - U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,1), :)
         C = U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,4), :) - U3D_GRID%NODE_COORDS(U3D_GRID%CELL_NODES(I,1), :)

         CELL_VOLUMES(I) = ABS(C(1)*(A(2)*B(3)-A(3)*B(2)) + C(2)*(A(3)*B(1)-A(1)*B(3)) + C(3)*(A(1)*B(2)-A(2)*B(1))) / 6.
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing grid connectivity.'
         WRITE(*,*) '==========================================='
      END IF


      ALLOCATE(TEMP_CELL_NEIGHBORS(U3D_GRID%NUM_CELLS, 4))
      TEMP_CELL_NEIGHBORS = -1

      DO JN = 1, U3D_GRID%NUM_NODES
         !IF (PROC_ID == 0) WRITE(*,*) 'Checking node ', JN, ' of ',  U3D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)


                  FOUND = 0
                  DO V1 = 1, 4
                     DO V2 = 1, 4
                        IF (U3D_GRID%CELL_NODES(JC1,V1) == U3D_GRID%CELL_NODES(JC2,V2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 3) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           WHICH1(FOUND) = V1
                           WHICH2(FOUND) = V2
                        END IF
                     END DO
                  END DO

                  IF (FOUND == 3) THEN
      
                     IF (ANY(WHICH1 == 1)) THEN
                        IF (ANY(WHICH1 == 2)) THEN
                           IF (ANY(WHICH1 == 3))  THEN
                              TEMP_CELL_NEIGHBORS(JC1, 1) = JC2
                           ELSE IF (ANY(WHICH1 == 4)) THEN
                              TEMP_CELL_NEIGHBORS(JC1, 2) = JC2
                           END IF
                        ELSE IF (ANY(WHICH1 == 3)) THEN
                           IF (ANY(WHICH1 == 4)) TEMP_CELL_NEIGHBORS(JC1, 4) = JC2
                        END IF
                     ELSE IF (ANY(WHICH1 == 2)) THEN
                        IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 4)) TEMP_CELL_NEIGHBORS(JC1, 3) = JC2
                     END IF
      
      
                     IF (ANY(WHICH2 == 1)) THEN
                        IF (ANY(WHICH2 == 2)) THEN
                           IF (ANY(WHICH2 == 3))  THEN
                              TEMP_CELL_NEIGHBORS(JC2, 1) = JC1
                           ELSE IF (ANY(WHICH2 == 4)) THEN
                              TEMP_CELL_NEIGHBORS(JC2, 2) = JC1
                           END IF
                        ELSE IF (ANY(WHICH2 == 3)) THEN
                           IF (ANY(WHICH2 == 4)) TEMP_CELL_NEIGHBORS(JC2, 4) = JC1
                        END IF
                     ELSE IF (ANY(WHICH2 == 2)) THEN
                        IF (ANY(WHICH2 == 3) .AND. ANY(WHICH2 == 4)) TEMP_CELL_NEIGHBORS(JC2, 3) = JC1
                     END IF
      
                  END IF


               END DO
            END DO
         END IF
      END DO

      U3D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing face normals.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell edge normals
      IND(1,:) = [1,3,2]
      IND(2,:) = [1,2,4]
      IND(3,:) = [2,3,4]
      IND(4,:) = [1,4,3]
      ALLOCATE(U3D_GRID%FACE_NORMAL(U3D_GRID%NUM_CELLS, 4, 3))
      ALLOCATE(U3D_GRID%FACE_TANG1(U3D_GRID%NUM_CELLS, 4, 3))
      ALLOCATE(U3D_GRID%FACE_TANG2(U3D_GRID%NUM_CELLS, 4, 3))
      ALLOCATE(U3D_GRID%FACE_NODES(U3D_GRID%NUM_CELLS, 4, 3))
      ALLOCATE(U3D_GRID%CELL_FACES_COEFFS(U3D_GRID%NUM_CELLS, 4, 4))
      ALLOCATE(U3D_GRID%FACE_AREA(U3D_GRID%NUM_CELLS, 4))

      DO I = 1, U3D_GRID%NUM_CELLS
         DO J = 1, 4
            V1 = U3D_GRID%CELL_NODES(I,IND(J,1))
            V2 = U3D_GRID%CELL_NODES(I,IND(J,2))
            V3 = U3D_GRID%CELL_NODES(I,IND(J,3))

            U3D_GRID%FACE_NODES(I,J,1) = V1
            U3D_GRID%FACE_NODES(I,J,2) = V2
            U3D_GRID%FACE_NODES(I,J,3) = V3

            A = U3D_GRID%NODE_COORDS(V1,:)
            B = U3D_GRID%NODE_COORDS(V2,:)
            C = U3D_GRID%NODE_COORDS(V3,:)
            
            CROSSP = CROSS(B-A,C-A)
            U3D_GRID%FACE_NORMAL(I,J,:) = CROSSP/NORM2(CROSSP)
            U3D_GRID%FACE_AREA(I,J) = 0.5*NORM2(CROSSP)
            
            U3D_GRID%FACE_TANG1(I,J,:) = (B-A)/NORM2(B-A)
            U3D_GRID%FACE_TANG2(I,J,:) = CROSS(U3D_GRID%FACE_NORMAL(I,J,:), U3D_GRID%FACE_TANG1(I,J,:))


            ! The coefficients (a,b,c,d) of a*x + b*y + c*z + d = 0
            U3D_GRID%CELL_FACES_COEFFS(I,J,1) =  A(2)*B(3)-B(2)*A(3) &
                                                +B(2)*C(3)-C(2)*B(3) &
                                                +C(2)*A(3)-A(2)*C(3)
            U3D_GRID%CELL_FACES_COEFFS(I,J,2) = -A(1)*B(3)+B(1)*A(3) &
                                                -B(1)*C(3)+C(1)*B(3) &
                                                -C(1)*A(3)+A(1)*C(3)
            U3D_GRID%CELL_FACES_COEFFS(I,J,3) =  A(1)*B(2)-B(1)*A(2) &
                                                +B(1)*C(2)-C(1)*B(2) &
                                                +C(1)*A(2)-A(1)*C(2)
            U3D_GRID%CELL_FACES_COEFFS(I,J,4) = -A(1)*B(2)*C(3) &
                                                +A(1)*C(2)*B(3) &
                                                +B(1)*A(2)*C(3) &
                                                -B(1)*C(2)*A(3) &
                                                -C(1)*A(2)*B(3) &
                                                +C(1)*B(2)*A(3)

         END DO
      END DO



      NCELLS = U3D_GRID%NUM_CELLS
      NNODES = U3D_GRID%NUM_NODES



      ALLOCATE(U3D_GRID%BASIS_COEFFS(NCELLS,4,4))

      DO I = 1, NCELLS
         VOLUME = CELL_VOLUMES(I)
         V1 = U3D_GRID%CELL_NODES(I,1)
         V2 = U3D_GRID%CELL_NODES(I,2)
         V3 = U3D_GRID%CELL_NODES(I,3)
         V4 = U3D_GRID%CELL_NODES(I,4)

         X1 = U3D_GRID%NODE_COORDS(V1, 1)
         X2 = U3D_GRID%NODE_COORDS(V2, 1)
         X3 = U3D_GRID%NODE_COORDS(V3, 1)
         X4 = U3D_GRID%NODE_COORDS(V4, 1)
         Y1 = U3D_GRID%NODE_COORDS(V1, 2)
         Y2 = U3D_GRID%NODE_COORDS(V2, 2)
         Y3 = U3D_GRID%NODE_COORDS(V3, 2)
         Y4 = U3D_GRID%NODE_COORDS(V4, 2)
         Z1 = U3D_GRID%NODE_COORDS(V1, 3)
         Z2 = U3D_GRID%NODE_COORDS(V2, 3)
         Z3 = U3D_GRID%NODE_COORDS(V3, 3)
         Z4 = U3D_GRID%NODE_COORDS(V4, 3)


         ! These are such that PSI_i = SUM_j [ x_j * BASIS_COEFFS(IC,i,j) ] + BASIS_COEFFS(IC,i,4)

         U3D_GRID%BASIS_COEFFS(I,1,1) =  Y2*Z3-Y3*Z2 -Y2*Z4+Y4*Z2 +Y3*Z4-Y4*Z3
         U3D_GRID%BASIS_COEFFS(I,1,2) = -X2*Z3+X3*Z2 +X2*Z4-X4*Z2 -X3*Z4+X4*Z3
         U3D_GRID%BASIS_COEFFS(I,1,3) =  X2*Y3-X3*Y2 -X2*Y4+X4*Y2 +X3*Y4-X4*Y3
         U3D_GRID%BASIS_COEFFS(I,1,4) = -X2*Y3*Z4 +X3*Y2*Z4 +X2*Y4*Z3 -X4*Y2*Z3 -X3*Y4*Z2 +X4*Y3*Z2

         U3D_GRID%BASIS_COEFFS(I,2,1) = -Y1*Z3+Y3*Z1 +Y1*Z4-Y4*Z1 -Y3*Z4+Y4*Z3
         U3D_GRID%BASIS_COEFFS(I,2,2) =  X1*Z3-X3*Z1 -X1*Z4+X4*Z1 +X3*Z4-X4*Z3
         U3D_GRID%BASIS_COEFFS(I,2,3) = -X1*Y3+X3*Y1 +X1*Y4-X4*Y1 -X3*Y4+X4*Y3
         U3D_GRID%BASIS_COEFFS(I,2,4) =  X1*Y3*Z4 -X3*Y1*Z4 -X1*Y4*Z3 +X4*Y1*Z3 +X3*Y4*Z1 -X4*Y3*Z1

         U3D_GRID%BASIS_COEFFS(I,3,1) =  Y1*Z2-Y2*Z1 -Y1*Z4+Y4*Z1 +Y2*Z4-Y4*Z2
         U3D_GRID%BASIS_COEFFS(I,3,2) = -X1*Z2+X2*Z1 +X1*Z4-X4*Z1 -X2*Z4+X4*Z2
         U3D_GRID%BASIS_COEFFS(I,3,3) =  X1*Y2-X2*Y1 -X1*Y4+X4*Y1 +X2*Y4-X4*Y2
         U3D_GRID%BASIS_COEFFS(I,3,4) = -X1*Y2*Z4 +X2*Y1*Z4 +X1*Y4*Z2 -X4*Y1*Z2 -X2*Y4*Z1 +X4*Y2*Z1

         U3D_GRID%BASIS_COEFFS(I,4,1) = -Y1*Z2+Y2*Z1 +Y1*Z3-Y3*Z1 -Y2*Z3+Y3*Z2
         U3D_GRID%BASIS_COEFFS(I,4,2) =  X1*Z2-X2*Z1 -X1*Z3+X3*Z1 +X2*Z3-X3*Z2
         U3D_GRID%BASIS_COEFFS(I,4,3) = -X1*Y2+X2*Y1 +X1*Y3-X3*Y1 -X2*Y3+X3*Y2
         U3D_GRID%BASIS_COEFFS(I,4,4) =  X1*Y2*Z3 -X2*Y1*Z3 -X1*Y3*Z2 +X3*Y1*Z2 +X2*Y3*Z1 -X3*Y2*Z1

         U3D_GRID%BASIS_COEFFS(I,:,:) = -U3D_GRID%BASIS_COEFFS(I,:,:)/6./VOLUME

      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Done reading grid file.'
         WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
         WRITE(*,*) '==========================================='
      END IF

   END SUBROUTINE READ_3D_UNSTRUCTURED_GRID_SU2


END MODULE grid_and_partition
