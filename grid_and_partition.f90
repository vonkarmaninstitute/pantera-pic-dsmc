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

         IF (XCELL .GT. (NX-1)) THEN 
            XCELL = NX-1
            WRITE(*,*) 'Particle out of bound xhi!'
         ELSE IF (XCELL .LT. 0) THEN
            XCELL = 0
            WRITE(*,*) 'Particle out of bound xlo!'
         END IF

         IF (YCELL .GT. (NY-1)) THEN 
            YCELL = NY-1
            WRITE(*,*) 'Particle out of bound yhi!'
         ELSE IF (YCELL .LT. 0) THEN
            YCELL = 0
            WRITE(*,*) 'Particle out of bound ylo!'
         END IF

         IDCELL = XCELL + NX*YCELL + 1
      ELSE IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) THEN
         XCELL = BINARY_SEARCH(XP, XCOORD)
         YCELL = BINARY_SEARCH(YP, YCOORD)

         IF (XCELL .GT. (NX)) THEN 
            XCELL = NX
            WRITE(*,*) 'Particle out of bound xhi!', XCELL, NX
         ELSE IF (XCELL .LT. 1) THEN
            XCELL = 1
            WRITE(*,*) 'Particle out of bound xlo!'
         END IF

         IF (YCELL .GT. (NY)) THEN 
            YCELL = NY
            WRITE(*,*) 'Particle out of bound yhi!'
         ELSE IF (YCELL .LT. 1) THEN
            YCELL = 1
            WRITE(*,*) 'Particle out of bound ylo!'
         END IF

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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE PROC_FROM_CELL -> Establishes to which processor IDPROC the cell     !
   ! IDCELL belongs to. IDCELL starts from 1, IDPROC starts from 0                   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE PROC_FROM_CELL(IDCELL, IDPROC)

      ! Note: processes go from 0 (usually termed the Master) to MPI_N_THREADS - 1
      ! If the number of MPI processes is only 1, then IDPROC is 0, the one and only process.
      ! Otherwise, check the partition style (variable "DOMPART_TYPE")

      IMPLICIT NONE

      INTEGER, INTENT(IN)     :: IDCELL
      INTEGER, INTENT(OUT)     :: IDPROC

      INTEGER :: NCELLS, NCELLSPP
      
      IF (N_MPI_THREADS == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Serial operation

         IDPROC = 0      

      ELSE IF (DOMPART_TYPE == 0) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ S T R I P S
      !
      ! Domain partition type: strips
      !
      ! This domain partition assigns the same number of grid cells to every process.
      ! The domain partition is thus in strips along the "x" axis, as the cells are assumed 
      ! to be numbered along x.
         IF (GRID_TYPE == UNSTRUCTURED) THEN
            NCELLS = U2D_GRID%NUM_CELLS
         ELSE
            NCELLS = NX*NY
         END IF
         NCELLSPP = CEILING(REAL(NCELLS)/REAL(N_MPI_THREADS))

         ! 3) Here is the process ID 
         IDPROC   = INT((IDCELL-1)/NCELLSPP) ! Before was giving wrong result for peculiar combinations of NCELLS and NCELLSPP
         IF (IDPROC .GT. N_MPI_THREADS-1 .OR. IDPROC .LT. 0) THEN
            WRITE(*,*) 'Error! PROC_FROM_CELL returned proc:', IDPROC
         END IF


      ELSE IF (DOMPART_TYPE == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ B L O C K S
      !
      ! Domain partition type: blocks
      !
      ! This domain partition style divides the domain into blocks
      ! Note that blocks in this definition are independent from cells! IT IS POSSIBLE TO
      ! GENERALIZE IT (and we could and should.. but for PIC I don't really care).

         IDPROC = -1 ! Fix this!

      ELSE ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ E R R O R

         WRITE(*,*) 'ATTENTION! Value for variable DOMPART_TYPE = ', DOMPART_TYPE
         WRITE(*,*) ' not recognized! Check input file!! In: PROC_FROM_POSITION() ABORTING!'
         STOP

      END IF

   END SUBROUTINE PROC_FROM_CELL


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE PROC_FROM_POSITION -> finds the process ID from particle position              !
   ! Note this depends on the parallelization strategy                                         !
   ! This should never be called when an unstructured grid is used.                            !
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
         IF (IDCELL .GT. NX*NY .OR. IDCELL .LT. 1) THEN
            WRITE(*,*) 'Error! CELL_FROM_POSITION returned cell:', IDCELL, 'Particle position: ', XP, ', ', YP
         END IF

         ! 2) Find number of cells for each process (NX*NY*NZ = number of cells)
         !    Exceed a bit, so the last processor has slightly less cells, if number
         !    of cells is not divisible by the MPI_threads
         NCELLSPP = CEILING(REAL(NX*NY)/REAL(N_MPI_THREADS))

         ! 3) Here is the process ID 
         IDPROC   = INT((IDCELL-1)/NCELLSPP) ! Before was giving wrong result for peculiar combinations of NCELLS and NCELLSPP
         IF (IDPROC .GT. N_MPI_THREADS-1 .OR. IDPROC .LT. 0) THEN
            WRITE(*,*) 'Error! PROC_FROM_POSITION returned proc:', IDPROC
         END IF


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
      ALLOCATE(U2D_GRID%CELL_PROC(NUM_CELLS))
      DO I = 1, NUM_CELLS
         U2D_GRID%CELL_PROC(I) = INT((I-1)/NCELLSPP)
      END DO

      GRID_TYPE = UNSTRUCTURED

      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Done reading grid file.'
      WRITE(*,*) '==========================================='

   END SUBROUTINE READ_UNSTRUCTURED_GRID





   SUBROUTINE READ_UNSTRUCTURED_GRID_SU2(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER            :: NUM, I, J, FOUND, V1, V2, ELEM_TYPE, NUMELEMS, IC
      INTEGER, DIMENSION(3,2) :: IND
      REAL(KIND=8)       :: X1, X2, Y1, Y2, LEN
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(2) :: WHICH1, WHICH2

      INTEGER :: NCELLSPP

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. SU2 file format (*.su2)
      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Reading grid file in SU2 format.'
      WRITE(*,*) '==========================================='
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
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cell edges.
            ALLOCATE(U2D_GRID%CELL_EDGES_PG(U2D_GRID%NUM_CELLS, 3))
            U2D_GRID%CELL_EDGES_PG = -1

            ALLOCATE(GRID_BC(NUM)) ! Append the physical group to the list
            N_GRID_BC = NUM

            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
               IF (LINE .NE. 'MARKER_TAG=') THEN
                  WRITE(*,*) 'Error! did not find marker name.'
               ELSE
                  WRITE(*,*) 'Found marker tag, with groupname: ', GROUPNAME
               END IF
      
               GRID_BC(I)%PHYSICAL_GROUP_NAME = GROUPNAME
         
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS
               IF (LINE .NE. 'MARKER_ELEMS=') THEN
                  WRITE(*,*) 'Error! did not find marker elements.'
               ELSE
                  WRITE(*,*) 'Found marker elements, with number of elements: ', NUMELEMS
               END IF

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, V1, V2
                  IF (ELEM_TYPE .NE. 3) WRITE(*,*) 'Error! element type was not line.'

                  V1 = V1 + 1
                  V2 = V2 + 1
                  DO IC = 1, U2D_GRID%NUM_CELLS
                     IF ((U2D_GRID%CELL_NODES(IC, 1) == V1 .AND. U2D_GRID%CELL_NODES(IC, 2) == V2) .OR. &
                         (U2D_GRID%CELL_NODES(IC, 1) == V2 .AND. U2D_GRID%CELL_NODES(IC, 2) == V1)) THEN
                        U2D_GRID%CELL_EDGES_PG(IC, 1) = I
                     ELSE IF ((U2D_GRID%CELL_NODES(IC, 2) == V1 .AND. U2D_GRID%CELL_NODES(IC, 3) == V2) .OR. &
                              (U2D_GRID%CELL_NODES(IC, 2) == V2 .AND. U2D_GRID%CELL_NODES(IC, 3) == V1)) THEN
                        U2D_GRID%CELL_EDGES_PG(IC, 2) = I
                     ELSE IF ((U2D_GRID%CELL_NODES(IC, 3) == V1 .AND. U2D_GRID%CELL_NODES(IC, 1) == V2) .OR. &
                              (U2D_GRID%CELL_NODES(IC, 3) == V2 .AND. U2D_GRID%CELL_NODES(IC, 1) == V1)) THEN
                        U2D_GRID%CELL_EDGES_PG(IC, 3) = I
                     END IF
                  END DO

               END DO
            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)

      WRITE(*,*) 'Read grid file. It contains ', U2D_GRID%NUM_NODES, &
                 'points, and ', U2D_GRID%NUM_CELLS, 'cells.'

      ! Process the mesh: generate connectivity, normals and such...
      !XMIN, XMAX,...

      !DO I = 1, U2D_GRID%NUM_CELLS
      !   WRITE(*,*) U2D_GRID%CELL_NODES(I,:)
      !END DO

      ! Compute cell volumes
      ALLOCATE(CELL_VOLUMES(U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         A = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,1), :)
         B = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,2), :)
         C = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,3), :)

         IF (DIMS == 2 .AND. .NOT. AXI) THEN
            CELL_VOLUMES(I) = 0.5*ABS(A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))
            !WRITE(*,*) CELL_VOLUMES(I)
         END IF
      END DO

      ! Find cell connectivity
      ALLOCATE(TEMP_CELL_NEIGHBORS(U2D_GRID%NUM_CELLS, 3))
      TEMP_CELL_NEIGHBORS = -1
      DO I = 1, U2D_GRID%NUM_CELLS
         DO J = I, U2D_GRID%NUM_CELLS
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


      !WRITE(*,*) 'Generated grid connectivity. '
      !DO I = 1, U2D_GRID%NUM_CELLS
      !   WRITE(*,*) 'Cell ', I, ' neighbors cells ', TEMP_CELL_NEIGHBORS(I, :)
      !END DO

      ! Distributing cells over MPI processes
      NCELLSPP = CEILING(REAL(U2D_GRID%NUM_CELLS)/REAL(N_MPI_THREADS))
      ALLOCATE(U2D_GRID%CELL_PROC(U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         U2D_GRID%CELL_PROC(I) = INT((I-1)/NCELLSPP)
      END DO

      ! Compute cell edge normals
      IND(1,:) = [1,2]
      IND(2,:) = [2,3]
      IND(3,:) = [3,1]
      ALLOCATE(U2D_GRID%EDGE_NORMAL(U2D_GRID%NUM_CELLS, 3, 2))
      DO I = 1, U2D_GRID%NUM_CELLS
         DO J = 1, 3
            X1 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,IND(J,1)), 1)
            X2 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,IND(J,2)), 1)
            Y1 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,IND(J,1)), 2)
            Y2 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(I,IND(J,2)), 2)
            LEN = SQRT((Y2-Y1)*(Y2-Y1) + (X2-X1)*(X2-X1))
            U2D_GRID%EDGE_NORMAL(I,J,1) = (Y2-Y1)/LEN
            U2D_GRID%EDGE_NORMAL(I,J,2) = (X1-X2)/LEN
         END DO
      END DO

      GRID_TYPE = UNSTRUCTURED

      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Done reading grid file.'
      WRITE(*,*) '==========================================='

      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Checking ordering.'
      WRITE(*,*) '==========================================='

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



   END SUBROUTINE READ_UNSTRUCTURED_GRID_SU2



END MODULE grid_and_partition
