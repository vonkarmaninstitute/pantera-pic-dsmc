! Contains Electromagnetic fields related code

MODULE fields

   USE mUMFPACK
   USE global
   USE screen

   IMPLICIT NONE

   CONTAINS

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ASSEMBLE_POISSON -> Prepares the linear system for the solution  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ASSEMBLE_POISSON

      INTEGER :: STATUS
      INTEGER :: I, J
      INTEGER :: IC, IS, IW, IN, IE
      INTEGER :: MAXNNZ, SIZE
      REAL(KIND=8) :: HX, HY
      REAL(KIND=8) :: AX, AY, BX, BY, CX, CY, H1X, H2X, H1Y, H2Y, R
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, K11, K22, K33, K12, K23, K13, AREA, EDGELENGTH
      INTEGER :: V1, V2, V3
      INTEGER :: EDGE_PG
      LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_UNUSED

      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZE = U2D_GRID%NUM_NODES
      ELSE
         SIZE = NPX*NPY
      END IF

      IF (.NOT. ALLOCATED(RHS)) ALLOCATE(RHS(0:SIZE-1))
      RHS = 0.d0

      IF (.NOT. ALLOCATED(PHI_FIELD)) ALLOCATE(PHI_FIELD(0:SIZE-1))
      PHI_FIELD = 0.d0

      IF (.NOT. ALLOCATED(DIRICHLET)) ALLOCATE(DIRICHLET(0:SIZE-1))
      IF (.NOT. ALLOCATED(IS_DIRICHLET)) ALLOCATE(IS_DIRICHLET(0:SIZE-1))
      IS_DIRICHLET = .FALSE.

      IF (.NOT. ALLOCATED(NEUMANN)) ALLOCATE(NEUMANN(0:SIZE-1))
      NEUMANN = 0.d0
      IF (.NOT. ALLOCATED(IS_NEUMANN)) ALLOCATE(IS_NEUMANN(0:SIZE-1))
      IS_NEUMANN = .FALSE.


      ! Create the matrix in Sparse Triplet format.
      ! Entries can be input in random order, but not duplicated.
      ! We have to allocate the ST matrix, but in this case
      ! it is not a problem, since each point takes at most 5 entries,
      ! we can use this as an upper limit.
      ! (We use a 5 point stencil for interior points)
      MAXNNZ = 10*SIZE
      CALL ST_MATRIX_ALLOCATE(A_ST, MAXNNZ)

      ! At this point, populate the matrix
      IF (GRID_TYPE == UNSTRUCTURED) THEN
         ALLOCATE(IS_UNUSED(0:SIZE-1))
         IS_UNUSED = .TRUE.

         DO I = 1, U2D_GRID%NUM_CELLS
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)
            IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.; IS_UNUSED(V3-1) = .FALSE.
            DO J = 1, 3
               EDGE_PG = U2D_GRID%CELL_EDGES_PG(I, J)
               IF (EDGE_PG .NE. -1) THEN
                  IF (GRID_BC(EDGE_PG)%FIELD_BC == DIRICHLET_BC) THEN
                     IF (J==1) THEN
                        DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        IS_DIRICHLET(V1-1) = .TRUE.; IS_DIRICHLET(V2-1) = .TRUE.
                     ELSE IF (J==2) THEN
                        DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        IS_DIRICHLET(V2-1) = .TRUE.; IS_DIRICHLET(V3-1) = .TRUE.
                     ELSE
                        DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        IS_DIRICHLET(V3-1) = .TRUE.; IS_DIRICHLET(V1-1) = .TRUE.
                     END IF
                  END IF

                  IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN
                     IF (J==1) THEN
                        IS_NEUMANN(V1-1) = .TRUE.; IS_NEUMANN(V2-1) = .TRUE.
                     ELSE IF (J==2) THEN
                        IS_NEUMANN(V2-1) = .TRUE.; IS_NEUMANN(V3-1) = .TRUE.
                     ELSE
                        IS_NEUMANN(V3-1) = .TRUE.; IS_NEUMANN(V1-1) = .TRUE.
                     END IF
                  END IF
               END IF
            END DO
         END DO

         !IS_DIRICHLET(0) = .TRUE.    ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDDBDBDBDBDBDBDB
         !DIRICHLET(0) = 0.d0

         DO I = 1, U2D_GRID%NUM_CELLS
            AREA = CELL_AREAS(I)
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)            
            X1 = U2D_GRID%NODE_COORDS(V1, 1)
            X2 = U2D_GRID%NODE_COORDS(V2, 1)
            X3 = U2D_GRID%NODE_COORDS(V3, 1)
            Y1 = U2D_GRID%NODE_COORDS(V1, 2)
            Y2 = U2D_GRID%NODE_COORDS(V2, 2)
            Y3 = U2D_GRID%NODE_COORDS(V3, 2)
            K11 = 0.25*((Y2-Y3)**2 + (X2-X3)**2)/AREA
            K22 = 0.25*((Y1-Y3)**2 + (X1-X3)**2)/AREA
            K33 = 0.25*((Y2-Y1)**2 + (X2-X1)**2)/AREA
            K12 =-0.25*((Y2-Y3)*(Y1-Y3) + (X2-X3)*(X1-X3))/AREA
            K23 = 0.25*((Y1-Y3)*(Y2-Y1) + (X1-X3)*(X2-X1))/AREA
            K13 =-0.25*((Y2-Y3)*(Y2-Y1) + (X2-X3)*(X2-X1))/AREA

            ! We need to ADD to a sparse matrix entry.
            IF (IS_DIRICHLET(V1-1)) THEN
               CALL ST_MATRIX_SET(A_ST, V1-1, V1-1, 1.d0)
            ELSE !IF (.NOT. IS_NEUMANN(V1-1)) THEN
               CALL ST_MATRIX_ADD(A_ST, V1-1, V1-1, K11)
               CALL ST_MATRIX_ADD(A_ST, V1-1, V2-1, K12)
               CALL ST_MATRIX_ADD(A_ST, V1-1, V3-1, K13)
            END IF
            IF (IS_DIRICHLET(V2-1)) THEN
               CALL ST_MATRIX_SET(A_ST, V2-1, V2-1, 1.d0)
            ELSE !IF (.NOT. IS_NEUMANN(V2-1)) THEN
               CALL ST_MATRIX_ADD(A_ST, V2-1, V1-1, K12)
               CALL ST_MATRIX_ADD(A_ST, V2-1, V3-1, K23)
               CALL ST_MATRIX_ADD(A_ST, V2-1, V2-1, K22)
            END IF
            IF (IS_DIRICHLET(V3-1)) THEN
               CALL ST_MATRIX_SET(A_ST, V3-1, V3-1, 1.d0)
            ELSE !IF (.NOT. IS_NEUMANN(V3-1)) THEN
               CALL ST_MATRIX_ADD(A_ST, V3-1, V1-1, K13)
               CALL ST_MATRIX_ADD(A_ST, V3-1, V2-1, K23)
               CALL ST_MATRIX_ADD(A_ST, V3-1, V3-1, K33)
            END IF

            DO J = 1, 3
               EDGE_PG = U2D_GRID%CELL_EDGES_PG(I, J)
               IF (EDGE_PG == -1) CYCLE
               IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN

                  IF (J==1) THEN
                     EDGELENGTH = SQRT((X2-X1)**2 + (Y2-Y1)**2)
                     IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                        NEUMANN(V1-1) = NEUMANN(V1-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
                     END IF
                     IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                        NEUMANN(V2-1) = NEUMANN(V2-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
                     END IF
                  ELSE IF (J==2) THEN
                     EDGELENGTH = SQRT((X3-X2)**2 + (Y3-Y2)**2)
                     IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                        NEUMANN(V2-1) = NEUMANN(V2-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
                     END IF
                     IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                        NEUMANN(V3-1) = NEUMANN(V3-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
                     END IF
                  ELSE
                     EDGELENGTH = SQRT((X1-X3)**2 + (Y1-Y3)**2)
                     IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                        NEUMANN(V1-1) = NEUMANN(V1-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
                     END IF
                     IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                        NEUMANN(V3-1) = NEUMANN(V3-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
                     END IF
                  END IF

               END IF
            END DO

         END DO

         DO I = 1, U2D_GRID%NUM_NODES
            IF (IS_UNUSED(I-1)) THEN
               CALL ST_MATRIX_SET(A_ST, I-1, I-1, 1.d0)
               IS_DIRICHLET(I-1) = .TRUE.
               DIRICHLET(I-1) = 0.d0
            END IF
         END DO
         DEALLOCATE(IS_UNUSED)

      ELSE
         DO I = 0, NPX-1
            DO J = 0, NPY-1


               IC = I+(NPX)*J
               IE = I+(NPX)*J+1
               IW = I+(NPX)*J-1
               IN = I+(NPX)*(J+1)
               IS = I+(NPX)*(J-1)

               HX = (XMAX-XMIN)/DBLE(NX)
               AX = 1./(HX*HX)


               ! Various configurations go here.
               
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Plume 2d cartesian, full domain !
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! IF ((I == 0) .OR. (I == NPX-1) .OR. (J == 0) .OR. (J == NPY-1) .OR. &
               ! ((I .GE. 10) .AND. (I .LE. 50) .AND. (J .GE. 50) .AND. (J .LE. 150)) ) THEN
               !    ! Boundary point
               !    CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)

               !    IF ((I == 50) .AND. (J .GE. 84) .AND. (J .LE. 116)) THEN
               !       ! On the inlet surface.
               !       DIRICHLET(IC) = 2.35d0
               !    ELSE
               !       ! On the boundary or on the rest of the PFG.
               !       DIRICHLET(IC) = 0.d0
               !    END IF
               !    IS_DIRICHLET(IC) = .TRUE.
               ! ELSE



               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Plume 2d axisymmetric (half domain) !
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               IF (J == 0) THEN
                  ! Axis: zero radial gradient.
                  CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)
                  DIRICHLET(IC) = 0.d0
                  CALL ST_MATRIX_SET(A_ST, IC, IN, -1.d0)
                  IS_DIRICHLET(IC) = .TRUE.
               ELSE IF ((I == 0) .OR. (I == NPX-1) .OR. (J == NPY-1) .OR. &
                  ((I .GE. 10) .AND. (I .LE. 30) .AND. (J .LE. 100)) ) THEN
                  ! Boundary point.
                  CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)
                  IF ((I == 30) .AND. (J .LE. 50)) THEN
                  ! On the inlet surface.
                     DIRICHLET(IC) = 2.35d0
                  ELSE
                     ! On the boundary or on the rest of the PFG.
                     DIRICHLET(IC) = 0.d0
                  END IF
                  IS_DIRICHLET(IC) = .TRUE.






               ELSE
                  ! Interior point.

                  IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
                     HX = (XMAX-XMIN)/DBLE(NX)
                     HY = (YMAX-YMIN)/DBLE(NY)
                     AX = 1./(HX*HX)
                     CX = AX
                     BX = -AX-CX
                     AY = 1./(HY*HY)
                     CY = AY
                     BY = -AY-CY
                  ELSE
                     H1X = XSIZE(I)
                     H2X = XSIZE(I+1)
                     H1Y = YSIZE(J)
                     H2Y = YSIZE(J+1)

                     AX = 2./(H1X*(H1X+H2X))
                     CX = 2./(H2X*(H1X+H2X))
                     BX = -AX-CX

                     AY = 2./(H1Y*(H1Y+H2Y))
                     CY = 2./(H2Y*(H1Y+H2Y))
                     BY = -AY-CY
                  END IF
      
                  IF (AXI) THEN
                     IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
                        R = YMIN + J*(YMAX-YMIN)/DBLE(NY)
                     ELSE
                        R = YCOORD(J+1)
                     END IF
                     AY = AY - 1./(R*(H1Y+H2Y))
                     CY = CY + 1./(R*(H1Y+H2Y))
                  END IF
                  
                  IF (DIMS == 2) THEN
                     CALL ST_MATRIX_SET(A_ST, IC, IC, BX+BY)
                     CALL ST_MATRIX_SET(A_ST, IC, IN, CY)
                     CALL ST_MATRIX_SET(A_ST, IC, IS, AY)
                     CALL ST_MATRIX_SET(A_ST, IC, IE, CX)
                     CALL ST_MATRIX_SET(A_ST, IC, IW, AX)
                  ELSE
                     CALL ST_MATRIX_SET(A_ST, IC, IC, BX)
                     CALL ST_MATRIX_SET(A_ST, IC, IE, CX)
                     CALL ST_MATRIX_SET(A_ST, IC, IW, AX)
                  END IF
               END IF


            END DO
         END DO
      END IF

      ! Factorize the matrix for later solution
      !IF (PROC_ID .EQ. 0) CALL ST_MATRIX_PRINT(A_ST, SIZE)

      ! IF (PROC_ID == 0) THEN
      !    DO J = 0, A_ST%NNZ-1
      !       WRITE(*,*) 'R ',  A_ST%RIDX(J), ' C ', A_ST%CIDX(J), ' V ', A_ST%VALUE(J)
      !    END DO
      ! ELSE
      !    CALL SLEEP(5)
      ! END IF
      IF (PROC_ID .EQ. 0) THEN
         CALL ST_MATRIX_TO_CC(A_ST, SIZE, A_CC)
         CALL S_UMFPACK_SYMBOLIC(SIZE, SIZE, A_CC%AP, A_CC%AI, A_CC%AX, STATUS = STATUS)
         WRITE(*,*) 'Status is = ', STATUS
         CALL S_UMFPACK_NUMERIC(A_CC%AP, A_CC%AI, A_CC%AX, STATUS = STATUS)
         WRITE(*,*) 'Status 1 is = ', STATUS
         CALL S_UMFPACK_FREE_SYMBOLIC
      END IF

   END SUBROUTINE ASSEMBLE_POISSON



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE SOLVE_POISSON -> Solves the Poisson equation with the RHS RHS !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE SOLVE_POISSON

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X
      REAL(KIND=8) :: HX, HY
      INTEGER :: I, J
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, AREA
      INTEGER :: V1, V2, V3, SIZE
      INTEGER :: IC, IN, IS, IE, IW

      HX = (XMAX-XMIN)/DBLE(NX)
      HY = (YMAX-YMIN)/DBLE(NY)


      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZE = U2D_GRID%NUM_NODES
      ELSE
         SIZE = NPX*NPY
      END IF

      ALLOCATE(X(0:SIZE-1))

      ! Solve the linear system
      IF (PROC_ID .EQ. 0) THEN
         !WRITE(*,*) 'Solving poisson'
         !WRITE(*,*) RHS
         CALL S_UMFPACK_SOLVE(UMFPACK_A, A_CC%AP, A_CC%AI, A_CC%AX, X, RHS)
         !WRITE(*,*) 'Solution:'
         !WRITE(*,*) X
         !X = 0.d0
      END IF
 
      CALL MPI_BCAST(X, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      PHI_FIELD = X
      PHIBAR_FIELD = X
      DEALLOCATE(X)

      ! Reshape the linear array into a 2D array
      IF (GRID_TYPE == UNSTRUCTURED) THEN

         ! Compute the electric field at grid points
         DO I = 1, U2D_GRID%NUM_CELLS
            AREA = CELL_AREAS(I)
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)            
            X1 = U2D_GRID%NODE_COORDS(V1, 1)
            X2 = U2D_GRID%NODE_COORDS(V2, 1)
            X3 = U2D_GRID%NODE_COORDS(V3, 1)
            Y1 = U2D_GRID%NODE_COORDS(V1, 2)
            Y2 = U2D_GRID%NODE_COORDS(V2, 2)
            Y3 = U2D_GRID%NODE_COORDS(V3, 2)

            E_FIELD(I,1,1) = -0.5/AREA*(  PHI_FIELD(V1-1)*(Y2-Y3) &
                                        - PHI_FIELD(V2-1)*(Y1-Y3) &
                                        - PHI_FIELD(V3-1)*(Y2-Y1))
            E_FIELD(I,1,2) = -0.5/AREA*(- PHI_FIELD(V1-1)*(X2-X3) &
                                        + PHI_FIELD(V2-1)*(X1-X3) &
                                        + PHI_FIELD(V3-1)*(X2-X1))
            E_FIELD(I,1,3) = 0.d0
         END DO

      ELSE
         
         ! Compute the electric field at grid points
         DO I = 0, NPX-1
            DO J = 0, NPY-1

               IC = I+1 + NPX*(J+1)
               IE = I+2 + NPX*(J+1)
               IW = I   + NPX*(J+1)
               IN = I+1 + NPX*(J+2)
               IS = I+1 + NPX*J

               IF (I == 0) THEN ! Left boundary
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = XSIZE(1)
                  E_FIELD(I,J,1) = (PHI_FIELD(IC)-PHI_FIELD(IE))/HX
               ELSE IF (I == 30 .AND. (J .LE. 50)) THEN ! Right side of PFG
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = XSIZE(I+1)
                  E_FIELD(I,J,1) = (PHI_FIELD(IC)-PHI_FIELD(IE))/HX
               ELSE IF (I == NPX-1) THEN ! Right boundary
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = XSIZE(NPX-1)
                  E_FIELD(I,J,1) = (PHI_FIELD(IW)-PHI_FIELD(IC))/HX
               ELSE IF (I == 10 .AND. (J .LE. 50)) THEN ! Left side of PFG
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = XSIZE(I)
                  E_FIELD(I,J,1) = (PHI_FIELD(IW)-PHI_FIELD(IC))/HX
               ELSE ! Interior point
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = 0.5*(XSIZE(I)+XSIZE(I+1))
                  E_FIELD(I,J,1) = 0.5*(PHI_FIELD(IW)-PHI_FIELD(IE))/HX
               END IF
               IF (DIMS == 2) THEN
                  IF (J == 0) THEN ! Bottom boundary
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = YSIZE(1)
                     E_FIELD(I,J,2) = (PHI_FIELD(IC)-PHI_FIELD(IN))/HY
                  ELSE IF (J == 50 .AND. ((I.GE. 10) .AND. (I .LE. 30))) THEN ! Top side of PFG
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = YSIZE(J+1)
                     E_FIELD(I,J,2) = (PHI_FIELD(IC)-PHI_FIELD(IN))/HY
                  ELSE IF (J == NPY-1) THEN ! Top boundary
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = YSIZE(NPY-1)
                     E_FIELD(I,J,2) = (PHI_FIELD(IS)-PHI_FIELD(IC))/HY
                  ! ELSE IF (J == 50 .AND. ((I.GE. 10) .AND. (I .LE. 50))) THEN ! Bottom side of PFG
                  !    IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = YSIZE(J)
                  !    E_FIELD(I,J,2) = (PHI_FIELD(IS)-PHI_FIELD(IC))/HY
                  ELSE ! Interior point
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = 0.5*(YSIZE(J)+YSIZE(J+1))
                     E_FIELD(I,J,2) = 0.5*(PHI_FIELD(IS)-PHI_FIELD(IN))/HY
                  END IF
               ELSE
                  E_FIELD(I,J,2) = 0.d0
               END IF
               E_FIELD(I,J,3) = 0.d0
            END DO
         END DO
         
      END IF

   END SUBROUTINE SOLVE_POISSON


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DEPOSIT_CHARGE -> Deposits the charge of particles on grid points !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DEPOSIT_CHARGE
      
      IMPLICIT NONE

      INTEGER :: JP, I, IC

      REAL(KIND=8) :: K, RHO_Q, CHARGE
      REAL(KIND=8) :: VOL, CFNUM
      REAL(KIND=8), DIMENSION(4) :: WEIGHTS
      INTEGER, DIMENSION(4) :: INDICES, INDI, INDJ
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, XP, YP
      INTEGER :: V1, V2, V3, SIZE
      REAL(KIND=8) :: PSI1, PSI2, PSI3

      K = QE/(EPS0*EPS_SCALING**2) ! [V m] Elementary charge / Dielectric constant of vacuum

      RHS = 0.d0

      DO JP = 1, NP_PROC
         CHARGE = SPECIES(particles(JP)%S_ID)%CHARGE
         IF (ABS(CHARGE) .LT. 1.d-6) CYCLE

         IF (GRID_TYPE == UNSTRUCTURED) THEN 
            IC = particles(JP)%IC
            VOL = CELL_VOLUMES(IC)
            V1 = U2D_GRID%CELL_NODES(IC,1)
            V2 = U2D_GRID%CELL_NODES(IC,2)
            V3 = U2D_GRID%CELL_NODES(IC,3)            
            X1 = U2D_GRID%NODE_COORDS(V1, 1)
            X2 = U2D_GRID%NODE_COORDS(V2, 1)
            X3 = U2D_GRID%NODE_COORDS(V3, 1)
            Y1 = U2D_GRID%NODE_COORDS(V1, 2)
            Y2 = U2D_GRID%NODE_COORDS(V2, 2)
            Y3 = U2D_GRID%NODE_COORDS(V3, 2)
            XP = particles(JP)%X
            YP = particles(JP)%Y
            PSI1 = 0.5*( (Y2-Y3)*(XP-X3) - (X2-X3)*(YP-Y3))/VOL
            PSI2 = 0.5*(-(Y1-Y3)*(XP-X3) + (X1-X3)*(YP-Y3))/VOL
            PSI3 = 0.5*(-(Y2-Y1)*(XP-X1) + (X2-X1)*(YP-Y1))/VOL
            
            RHO_Q = K*CHARGE*FNUM
            RHS(V1-1) = RHS(V1-1) + RHO_Q*PSI1
            RHS(V2-1) = RHS(V2-1) + RHO_Q*PSI2
            RHS(V3-1) = RHS(V3-1) + RHO_Q*PSI3

         ELSE

            CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)

            IF (GRID_TYPE == RECTILINEAR_UNIFORM .AND. .NOT. AXI) THEN
               VOL = CELL_VOL
            ELSE
               VOL = CELL_VOLUMES(particles(JP)%IC)
            END IF

            CFNUM = FNUM
            IF (BOOL_RADIAL_WEIGHTING) CFNUM = CELL_FNUM(particles(JP)%IC)         

            RHO_Q = -K*CHARGE*CFNUM/VOL

            IF (DIMS == 2) THEN
               RHS(INDICES(1)) = RHS(INDICES(1)) + RHO_Q * WEIGHTS(1)
               RHS(INDICES(2)) = RHS(INDICES(2)) + RHO_Q * WEIGHTS(2)
               RHS(INDICES(3)) = RHS(INDICES(3)) + RHO_Q * WEIGHTS(3)
               RHS(INDICES(4)) = RHS(INDICES(4)) + RHO_Q * WEIGHTS(4)            
            ELSE
               RHS(INDICES(1)) = RHS(INDICES(1)) + RHO_Q * WEIGHTS(1)
               RHS(INDICES(2)) = RHS(INDICES(2)) + RHO_Q * WEIGHTS(2)
            END IF

         END IF
      END DO


      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZE = U2D_GRID%NUM_NODES
      ELSE
         SIZE = NPX*NPY
      END IF

      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, RHS, SIZE, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ELSE
         CALL MPI_REDUCE(RHS,          RHS, SIZE, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF


      DO I = 0, SIZE-1
         IF (IS_DIRICHLET(I)) THEN
            RHS(I)=DIRICHLET(I)
         ELSE IF (IS_NEUMANN(I)) THEN
            RHS(I)=RHS(I) + NEUMANN(I)
         END IF
      END DO


      ! RHS = 0
      ! ALLOCATE(Q_TEMP(NPX,NPY))
      ! Q_TEMP = 0
      ! DO I = NPX/4, 3*NPX/4
      !    DO J = NPY/4, 3*NPY/4
      !       Q_TEMP(I,J) = 1.d0
      !    END DO
      ! END DO
      ! RHS = PACK(Q_TEMP, .TRUE.)

   END SUBROUTINE DEPOSIT_CHARGE



   SUBROUTINE SOLVE_POISSON_NONLINEAR
      
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AX, YK, FK, BK, DBDPHI, PHIK
      REAL(KIND=8) :: RESIDUAL
      TYPE(ST_MATRIX) :: JAC_ST
      TYPE(CC_MATRIX) :: JAC_CC
      INTEGER :: J, IC, IN, SIZE, ITERATION, STATUS
      INTEGER :: I
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, AREA
      INTEGER :: V1, V2, V3

      IF (PROC_ID .EQ. 0) THEN

         SIZE = U2D_GRID%NUM_NODES

         ALLOCATE(PHIK, SOURCE = PHI_FIELD)
         ALLOCATE(BK(0:SIZE-1))
         ALLOCATE(DBDPHI(0:SIZE-1))
         ALLOCATE(YK(0:SIZE-1))
         IF (.NOT. ALLOCATED(BOLTZ_NRHOE)) ALLOCATE(BOLTZ_NRHOE(SIZE))

         DO ITERATION = 1, 10
            WRITE(*,*) '---> Starting iteration ', ITERATION
            
            WRITE(*,*) PHIK
            DO IN = 1, SIZE
               BOLTZ_NRHOE(IN) = BOLTZ_N0*EXP( QE*(PHIK(IN-1) - BOLTZ_PHI0) / (KB*BOLTZ_TE) )
            END DO

            BK = 0.d0
            DBDPHI = 0.d0
            DO IC = 1, U2D_GRID%NUM_CELLS
               DO J = 1, 3
                  IN = U2D_GRID%CELL_NODES(IC,J)
                  IF (.NOT. IS_DIRICHLET(IN-1)) THEN
                     BK(IN-1) = BK(IN-1) + QE/EPS0*CELL_AREAS(IC)/3*BOLTZ_NRHOE(IN)
                     DBDPHI(IN-1) = DBDPHI(IN-1) + QE*QE/(EPS0*KB*BOLTZ_TE)*CELL_AREAS(IC)/3*BOLTZ_NRHOE(IN)
                  END IF
               END DO
            END DO

            !CALL ST_MATRIX_MULT(A_ST, PHIK, AX)
            ALLOCATE( AX(0:SIZE-1) )
            AX = 0.d0
      
            DO J = 0, A_ST%NNZ-1
               AX(A_ST%RIDX(J)) = AX(A_ST%RIDX(J)) + A_ST%VALUE(J)*PHIK(A_ST%CIDX(J))
            END DO

            ALLOCATE(FK(0:SIZE-1))
            FK = AX - RHS - BK
            DEALLOCATE(AX)

            JAC_ST%NNZ = A_ST%NNZ
            ALLOCATE(JAC_ST%VALUE, SOURCE = A_ST%VALUE)
            ALLOCATE(JAC_ST%RIDX, SOURCE = A_ST%RIDX)
            ALLOCATE(JAC_ST%CIDX, SOURCE = A_ST%CIDX)

            DO IN = 1, SIZE
               CALL ST_MATRIX_ADD(JAC_ST, IN-1, IN-1, -DBDPHI(IN-1))
            END DO

            CALL ST_MATRIX_TO_CC(JAC_ST, SIZE, JAC_CC)
            CALL S_UMFPACK_SYMBOLIC(SIZE, SIZE, JAC_CC%AP, JAC_CC%AI, JAC_CC%AX, STATUS = STATUS)
            CALL S_UMFPACK_NUMERIC(JAC_CC%AP, JAC_CC%AI, JAC_CC%AX, STATUS = STATUS)
            CALL S_UMFPACK_FREE_SYMBOLIC
            CALL S_UMFPACK_SOLVE(UMFPACK_A, JAC_CC%AP, JAC_CC%AI, JAC_CC%AX, YK, FK)
            DEALLOCATE(FK)
            PHIK = PHIK - YK
            CALL S_UMFPACK_FREE_NUMERIC

            CALL ST_MATRIX_DEALLOCATE(JAC_ST)
            CALL CC_MATRIX_DEALLOCATE(JAC_CC)

            RESIDUAL = NORM2(YK)
            WRITE(*,*) 'Residual of this iteration is: ', RESIDUAL
            IF (RESIDUAL < 1.d-6) EXIT
         END DO
         PHI_FIELD = PHIK

         DEALLOCATE(BK)
         DEALLOCATE(DBDPHI)
         DEALLOCATE(PHIK)
         DEALLOCATE(YK)

         PHIBAR_FIELD = PHI_FIELD
         ! Compute the electric field at grid points
         DO I = 1, U2D_GRID%NUM_CELLS
            AREA = CELL_AREAS(I)
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)            
            X1 = U2D_GRID%NODE_COORDS(V1, 1)
            X2 = U2D_GRID%NODE_COORDS(V2, 1)
            X3 = U2D_GRID%NODE_COORDS(V3, 1)
            Y1 = U2D_GRID%NODE_COORDS(V1, 2)
            Y2 = U2D_GRID%NODE_COORDS(V2, 2)
            Y3 = U2D_GRID%NODE_COORDS(V3, 2)

            E_FIELD(I,1,1) = -0.5/AREA*(  PHI_FIELD(V1-1)*(Y2-Y3) &
                                        - PHI_FIELD(V2-1)*(Y1-Y3) &
                                        - PHI_FIELD(V3-1)*(Y2-Y1))
            E_FIELD(I,1,2) = -0.5/AREA*(- PHI_FIELD(V1-1)*(X2-X3) &
                                        + PHI_FIELD(V2-1)*(X1-X3) &
                                        + PHI_FIELD(V3-1)*(X2-X1))
            E_FIELD(I,1,3) = 0.d0
         END DO

      END IF

      CALL MPI_BCAST(PHI_FIELD, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   END SUBROUTINE SOLVE_POISSON_NONLINEAR


   SUBROUTINE ASSEMBLE_AMPERE

      INTEGER :: STATUS
      INTEGER :: I, J
      INTEGER :: MAXNNZ, SIZE
      TYPE(ST_MATRIX) :: A_ST
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, K11, K22, K33, K12, K23, K13, AREA
      INTEGER :: V1, V2, V3
      INTEGER :: EDGE_PG
      LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_UNUSED
      REAL(KIND=8) :: MM

      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZE = U2D_GRID%NUM_NODES
      ELSE
         SIZE = NPX*NPY
      END IF

      IF (.NOT. ALLOCATED(RHS)) ALLOCATE(RHS(0:SIZE-1))
      RHS = 0.d0

      IF (.NOT. ALLOCATED(DIRICHLET)) ALLOCATE(DIRICHLET(0:SIZE-1))
      IF (.NOT. ALLOCATED(IS_DIRICHLET)) ALLOCATE(IS_DIRICHLET(0:SIZE-1))
      IS_DIRICHLET = .FALSE.

      IF (.NOT. ALLOCATED(NEUMANN)) ALLOCATE(NEUMANN(0:SIZE-1))
      NEUMANN = 0.d0
      IF (.NOT. ALLOCATED(IS_NEUMANN)) ALLOCATE(IS_NEUMANN(0:SIZE-1))
      IS_NEUMANN = .FALSE.


      ! Create the matrix in Sparse Triplet format.
  
      CALL CC_MATRIX_DEALLOCATE(A_CC)
      MAXNNZ = 10*SIZE
      CALL ST_MATRIX_ALLOCATE(A_ST, MAXNNZ)

      ! At this point, populate the matrix
      IF (GRID_TYPE == UNSTRUCTURED) THEN
         ALLOCATE(IS_UNUSED(0:SIZE-1))
         IS_UNUSED = .TRUE.

         DO I = 1, U2D_GRID%NUM_CELLS
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)
            IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.; IS_UNUSED(V3-1) = .FALSE.
            DO J = 1, 3
               EDGE_PG = U2D_GRID%CELL_EDGES_PG(I, J)
               IF (EDGE_PG .NE. -1) THEN
                  IF (GRID_BC(EDGE_PG)%FIELD_BC == DIRICHLET_BC) THEN
                     IF (J==1) THEN
                        DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        IS_DIRICHLET(V1-1) = .TRUE.; IS_DIRICHLET(V2-1) = .TRUE.
                     ELSE IF (J==2) THEN
                        DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        IS_DIRICHLET(V2-1) = .TRUE.; IS_DIRICHLET(V3-1) = .TRUE.
                     ELSE
                        DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                        IS_DIRICHLET(V3-1) = .TRUE.; IS_DIRICHLET(V1-1) = .TRUE.
                     END IF
                  END IF

                  IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN
                     IF (J==1) THEN
                        IS_NEUMANN(V1-1) = .TRUE.; IS_NEUMANN(V2-1) = .TRUE.
                     ELSE IF (J==2) THEN
                        IS_NEUMANN(V2-1) = .TRUE.; IS_NEUMANN(V3-1) = .TRUE.
                     ELSE
                        IS_NEUMANN(V3-1) = .TRUE.; IS_NEUMANN(V1-1) = .TRUE.
                     END IF
                  END IF
               END IF
            END DO
         END DO

         !IS_DIRICHLET(0) = .TRUE. ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDDBDBDBDBDBDBDB
         !DIRICHLET(0) = 0.d0

         DO I = 1, U2D_GRID%NUM_CELLS
            AREA = CELL_AREAS(I)
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)            
            X1 = U2D_GRID%NODE_COORDS(V1, 1)
            X2 = U2D_GRID%NODE_COORDS(V2, 1)
            X3 = U2D_GRID%NODE_COORDS(V3, 1)
            Y1 = U2D_GRID%NODE_COORDS(V1, 2)
            Y2 = U2D_GRID%NODE_COORDS(V2, 2)
            Y3 = U2D_GRID%NODE_COORDS(V3, 2)
            K11 = 0.25*((Y2-Y3)**2 + (X2-X3)**2)/AREA
            K22 = 0.25*((Y1-Y3)**2 + (X1-X3)**2)/AREA
            K33 = 0.25*((Y2-Y1)**2 + (X2-X1)**2)/AREA
            K12 =-0.25*((Y2-Y3)*(Y1-Y3) + (X2-X3)*(X1-X3))/AREA
            K23 = 0.25*((Y1-Y3)*(Y2-Y1) + (X1-X3)*(X2-X1))/AREA
            K13 =-0.25*((Y2-Y3)*(Y2-Y1) + (X2-X3)*(X2-X1))/AREA
            MM = MASS_MATRIX(I)+1.

            ! We need to ADD to a sparse matrix entry.
            IF (IS_DIRICHLET(V1-1)) THEN
               CALL ST_MATRIX_SET(A_ST, V1-1, V1-1, 1.d0)
            ELSE !IF (.NOT. IS_NEUMANN(V1-1)) THEN
               CALL ST_MATRIX_ADD(A_ST, V1-1, V1-1, MM*K11)
               CALL ST_MATRIX_ADD(A_ST, V1-1, V2-1, MM*K12)
               CALL ST_MATRIX_ADD(A_ST, V1-1, V3-1, MM*K13)
               RHS(V1-1) = RHS(V1-1) + PHI_FIELD(V1-1)*K11+PHI_FIELD(V2-1)*K12+PHI_FIELD(V3-1)*K13
            END IF
            IF (IS_DIRICHLET(V2-1)) THEN
               CALL ST_MATRIX_SET(A_ST, V2-1, V2-1, 1.d0)
            ELSE !IF (.NOT. IS_NEUMANN(V2-1)) THEN
               CALL ST_MATRIX_ADD(A_ST, V2-1, V1-1, MM*K12)
               CALL ST_MATRIX_ADD(A_ST, V2-1, V3-1, MM*K23)
               CALL ST_MATRIX_ADD(A_ST, V2-1, V2-1, MM*K22)
               RHS(V2-1) = RHS(V2-1) + PHI_FIELD(V1-1)*K12+PHI_FIELD(V2-1)*K22+PHI_FIELD(V3-1)*K23
            END IF
            IF (IS_DIRICHLET(V3-1)) THEN
               CALL ST_MATRIX_SET(A_ST, V3-1, V3-1, 1.d0)
            ELSE !IF (.NOT. IS_NEUMANN(V3-1)) THEN
               CALL ST_MATRIX_ADD(A_ST, V3-1, V1-1, MM*K13)
               CALL ST_MATRIX_ADD(A_ST, V3-1, V2-1, MM*K23)
               CALL ST_MATRIX_ADD(A_ST, V3-1, V3-1, MM*K33)
               RHS(V3-1) = RHS(V3-1) + PHI_FIELD(V1-1)*K13+PHI_FIELD(V2-1)*K23+PHI_FIELD(V3-1)*K33
            END IF

         ! Neumann part has to be included only if the derivative changes in time.
         !    DO J = 1, 3
         !       EDGE_PG = U2D_GRID%CELL_EDGES_PG(I, J)
         !       IF (EDGE_PG == -1) CYCLE
         !       IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN

         !          IF (J==1) THEN
         !             EDGELENGTH = SQRT((X2-X1)**2 + (Y2-Y1)**2)
         !             IF (.NOT. IS_DIRICHLET(V1-1)) THEN
         !                NEUMANN(V1-1) = NEUMANN(V1-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
         !             END IF
         !             IF (.NOT. IS_DIRICHLET(V2-1)) THEN
         !                NEUMANN(V2-1) = NEUMANN(V2-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
         !             END IF
         !          ELSE IF (J==2) THEN
         !             EDGELENGTH = SQRT((X3-X2)**2 + (Y3-Y2)**2)
         !             IF (.NOT. IS_DIRICHLET(V2-1)) THEN
         !                NEUMANN(V2-1) = NEUMANN(V2-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
         !             END IF
         !             IF (.NOT. IS_DIRICHLET(V3-1)) THEN
         !                NEUMANN(V3-1) = NEUMANN(V3-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
         !             END IF
         !          ELSE
         !             EDGELENGTH = SQRT((X1-X3)**2 + (Y1-Y3)**2)
         !             IF (.NOT. IS_DIRICHLET(V1-1)) THEN
         !                NEUMANN(V1-1) = NEUMANN(V1-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
         !             END IF
         !             IF (.NOT. IS_DIRICHLET(V3-1)) THEN
         !                NEUMANN(V3-1) = NEUMANN(V3-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * 0.5*EDGELENGTH
         !             END IF
         !          END IF

         !       END IF
         !    END DO

         END DO

         DO I = 0, U2D_GRID%NUM_NODES-1
            IF (IS_UNUSED(I)) THEN
               CALL ST_MATRIX_SET(A_ST, I, I, 1.d0)
               IS_DIRICHLET(I) = .TRUE.
               DIRICHLET(I) = 0.d0
            ELSE IF (.NOT. IS_DIRICHLET(I) ) THEN
               RHS(I) = RHS(I) + 0.5*DT/EPS0*J_FIELD(I)
            END IF
         END DO
         DEALLOCATE(IS_UNUSED)

      ELSE
         CALL ERROR_ABORT('Implicit with cartesian grid not implemented!')
      END IF

      ! Factorize the matrix for later solution
      !IF (PROC_ID .EQ. 0) CALL ST_MATRIX_PRINT(A_ST, SIZE)

      ! IF (PROC_ID == 0) THEN
      !    DO J = 0, A_ST%NNZ-1
      !       WRITE(*,*) 'R ',  A_ST%RIDX(J), ' C ', A_ST%CIDX(J), ' V ', A_ST%VALUE(J)
      !    END DO
      ! ELSE
      !    CALL SLEEP(5)
      ! END IF

      DO I = 0, SIZE-1
         IF (IS_DIRICHLET(I)) THEN
            RHS(I)=DIRICHLET(I)
         ELSE IF (IS_NEUMANN(I)) THEN
            RHS(I)=RHS(I) + NEUMANN(I)
         END IF
      END DO

      IF (PROC_ID .EQ. 0) THEN
         CALL ST_MATRIX_TO_CC(A_ST, SIZE, A_CC)
         CALL S_UMFPACK_SYMBOLIC(SIZE, SIZE, A_CC%AP, A_CC%AI, A_CC%AX, STATUS = STATUS)
         WRITE(*,*) 'Status is = ', STATUS
         CALL S_UMFPACK_NUMERIC(A_CC%AP, A_CC%AI, A_CC%AX, STATUS = STATUS)
         WRITE(*,*) 'Status 1 is = ', STATUS
         CALL S_UMFPACK_FREE_SYMBOLIC
      END IF

   END SUBROUTINE ASSEMBLE_AMPERE


   SUBROUTINE DEPOSIT_CURRENT
      
      IMPLICIT NONE

      INTEGER :: JP, IC

      REAL(KIND=8) :: CHARGE
      REAL(KIND=8) :: VOL
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3
      INTEGER :: V1, V2, V3, SIZE, SIZEC
      REAL(KIND=8) :: DPSI1DX, DPSI2DX, DPSI3DX, DPSI1DY, DPSI2DY, DPSI3DY


      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZE = U2D_GRID%NUM_NODES
         SIZEC = U2D_GRID%NUM_CELLS
      ELSE
         CALL ERROR_ABORT('Not implemented.')
      END IF

      IF (.NOT. ALLOCATED(J_FIELD)) ALLOCATE(J_FIELD(0:SIZE-1))
      J_FIELD = 0.d0
      IF (.NOT. ALLOCATED(MASS_MATRIX)) ALLOCATE(MASS_MATRIX(SIZEC))
      MASS_MATRIX = 0.d0

      DO JP = 1, NP_PROC
         CHARGE = SPECIES(particles(JP)%S_ID)%CHARGE
         IF (ABS(CHARGE) .LT. 1.d-6) CYCLE

         IF (GRID_TYPE == UNSTRUCTURED) THEN 
            IC = particles(JP)%IC
            VOL = CELL_VOLUMES(IC)
            V1 = U2D_GRID%CELL_NODES(IC,1)
            V2 = U2D_GRID%CELL_NODES(IC,2)
            V3 = U2D_GRID%CELL_NODES(IC,3)            
            X1 = U2D_GRID%NODE_COORDS(V1, 1)
            X2 = U2D_GRID%NODE_COORDS(V2, 1)
            X3 = U2D_GRID%NODE_COORDS(V3, 1)
            Y1 = U2D_GRID%NODE_COORDS(V1, 2)
            Y2 = U2D_GRID%NODE_COORDS(V2, 2)
            Y3 = U2D_GRID%NODE_COORDS(V3, 2)
            DPSI1DX =  0.5*(Y2-Y3)/VOL
            DPSI2DX = -0.5*(Y1-Y3)/VOL
            DPSI3DX = -0.5*(Y2-Y1)/VOL
            DPSI1DY = -0.5*(X2-X3)/VOL
            DPSI2DY =  0.5*(X1-X3)/VOL
            DPSI3DY =  0.5*(X2-X1)/VOL
            
            J_FIELD(V1-1) = J_FIELD(V1-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI1DX + particles(JP)%VY*DPSI1DY)
            J_FIELD(V2-1) = J_FIELD(V2-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI2DX + particles(JP)%VY*DPSI2DY)
            J_FIELD(V3-1) = J_FIELD(V3-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI3DX + particles(JP)%VY*DPSI3DY)

            MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*particles(JP)%DTRIM/EPS0/VOL*FNUM &
                              * (QE*CHARGE)**2/SPECIES(particles(JP)%S_ID)%MOLECULAR_MASS
         ELSE

            CALL ERROR_ABORT('Not implemented.')

         END IF
      END DO


      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, J_FIELD, SIZE, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ELSE
         CALL MPI_REDUCE(J_FIELD,      J_FIELD, SIZE, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF


      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, MASS_MATRIX, SIZEC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ELSE
         CALL MPI_REDUCE(MASS_MATRIX,  MASS_MATRIX, SIZEC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF


   END SUBROUTINE DEPOSIT_CURRENT


   SUBROUTINE SOLVE_AMPERE

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X
      REAL(KIND=8) :: HX, HY
      INTEGER :: I
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, AREA
      INTEGER :: V1, V2, V3, SIZE

      HX = (XMAX-XMIN)/DBLE(NX)
      HY = (YMAX-YMIN)/DBLE(NY)


      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZE = U2D_GRID%NUM_NODES
      ELSE
         SIZE = NPX*NPY
      END IF

      ALLOCATE(X(0:SIZE-1))

      ! Solve the linear system
      IF (PROC_ID .EQ. 0) THEN
         !WRITE(*,*) 'Solving poisson'
         !WRITE(*,*) RHS
         CALL S_UMFPACK_SOLVE(UMFPACK_A, A_CC%AP, A_CC%AI, A_CC%AX, X, RHS)
         CALL S_UMFPACK_FREE_NUMERIC
         !WRITE(*,*) 'Solution:'
         !WRITE(*,*) X
         !X = 0.d0
      END IF
 
      CALL MPI_BCAST(X, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      PHIBAR_FIELD = X
      DEALLOCATE(X)

      PHI_FIELD = 2*PHIBAR_FIELD-PHI_FIELD

      ! Reshape the linear array into a 2D array
      IF (GRID_TYPE == UNSTRUCTURED) THEN

         ! Compute the electric field at grid points
         DO I = 1, U2D_GRID%NUM_CELLS
            AREA = CELL_AREAS(I)
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)            
            X1 = U2D_GRID%NODE_COORDS(V1, 1)
            X2 = U2D_GRID%NODE_COORDS(V2, 1)
            X3 = U2D_GRID%NODE_COORDS(V3, 1)
            Y1 = U2D_GRID%NODE_COORDS(V1, 2)
            Y2 = U2D_GRID%NODE_COORDS(V2, 2)
            Y3 = U2D_GRID%NODE_COORDS(V3, 2)

            E_FIELD(I,1,1) = -0.5/AREA*(  PHI_FIELD(V1-1)*(Y2-Y3) &
                                        - PHI_FIELD(V2-1)*(Y1-Y3) &
                                        - PHI_FIELD(V3-1)*(Y2-Y1))
            E_FIELD(I,1,2) = -0.5/AREA*(- PHI_FIELD(V1-1)*(X2-X3) &
                                        + PHI_FIELD(V2-1)*(X1-X3) &
                                        + PHI_FIELD(V3-1)*(X2-X1))
            E_FIELD(I,1,3) = 0.d0

            EBAR_FIELD(I,1,1) = -0.5/AREA*(  PHIBAR_FIELD(V1-1)*(Y2-Y3) &
                                           - PHIBAR_FIELD(V2-1)*(Y1-Y3) &
                                           - PHIBAR_FIELD(V3-1)*(Y2-Y1))
            EBAR_FIELD(I,1,2) = -0.5/AREA*(- PHIBAR_FIELD(V1-1)*(X2-X3) &
                                           + PHIBAR_FIELD(V2-1)*(X1-X3) &
                                           + PHIBAR_FIELD(V3-1)*(X2-X1))
            EBAR_FIELD(I,1,3) = 0.d0
         END DO

      ELSE
         
         CALL ERROR_ABORT('Not implemented.')
         
      END IF

   END SUBROUTINE SOLVE_AMPERE



   SUBROUTINE COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: JP
      REAL(KIND=8), DIMENSION(4), INTENT(OUT) :: WEIGHTS
      INTEGER, DIMENSION(4), INTENT(OUT) :: INDICES, INDI, INDJ
      REAL(KIND=8) :: XP, YP, CELL_XMIN, CELL_YMIN, FX, FY, HX, HY
      INTEGER :: IC, I, J

      XP = particles(JP)%X
      YP = particles(JP)%Y
      ! This is where we need the cell properties based on particle position
      ! With rectilinear non-uniform grid this has to be changed.
      IC = particles(JP)%IC
      IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
         HX = (XMAX-XMIN)/DBLE(NX)
         HY = (YMAX-YMIN)/DBLE(NY)
         I = MOD(IC-1, NX)
         J = (IC-1)/NX
         CELL_XMIN = I * HX + XMIN
         CELL_YMIN = J * HY + YMIN
      ELSE
         I = MOD(IC-1, NX)
         J = (IC-1)/NX
         HX = XSIZE(I+1)
         HY = YSIZE(J+1)
         CELL_XMIN = XCOORD(I+1)
         CELL_YMIN = YCOORD(J+1)
      END IF
      FX = (particles(JP)%X - CELL_XMIN)/HX
      FY = (particles(JP)%Y - CELL_YMIN)/HY

      !WRITE(*,*) 'FX=', FX, ', FY=', FY ! DBDBDBDBDBDBDBDDBDBDBDBD
      !IF ((FX .LT. 0.) .OR. (FX .GT. 1.) .OR. (FY .LT. 0.) .OR. (FY .GT. 1.)) WRITE(*,*) 'Whoops!'

      ! Nodes numbering
      !     _________
      !  3 |         | 4   ^ Y
      !    |   CELL  |     |
      !    |    IC   |     |
      !  1 |_________| 2   |------> X
      !

      IF (DIMS == 2) THEN
         ! Row and column indices, starting from 0
         INDI(1) = I
         INDI(2) = I + 1
         INDI(3) = I
         INDI(4) = I + 1

         INDJ(1) = J
         INDJ(2) = J
         INDJ(3) = J + 1
         INDJ(4) = J + 1
         
         INDICES(1) = IC+J-1
         INDICES(2) = INDICES(1) + 1
         INDICES(3) = INDICES(1) + NX + 1
         INDICES(4) = INDICES(3) + 1

         WEIGHTS(1) = (1.-FX)*(1.-FY)
         WEIGHTS(2) = FX*(1.-FY)
         WEIGHTS(3) = (1.-FX)*FY
         WEIGHTS(4) = FX*FY

      ELSE
         ! Row and column indices, starting from 0
         INDI(1) = I
         INDI(2) = I + 1
         INDI(3) = -1
         INDI(4) = -1

         INDJ(1) = 0
         INDJ(2) = 0
         INDJ(3) = -1
         INDJ(4) = -1
      
         ! Left, right
         INDICES(1) = IC
         INDICES(2) = IC + 1
         INDICES(4) = -1
         INDICES(3) = -1

         ! Left, right
         WEIGHTS(1) = (1.-FX)
         WEIGHTS(2) = FX
         WEIGHTS(3) = -1
         WEIGHTS(4) = -1
      END IF
   END SUBROUTINE COMPUTE_WEIGHTS


   SUBROUTINE APPLY_E_FIELD(JP, E)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: E
      INTEGER, INTENT(IN) :: JP
      REAL(KIND=8), DIMENSION(4) :: WEIGHTS
      INTEGER, DIMENSION(4) :: INDICES, INDI, INDJ

      IF (GRID_TYPE == UNSTRUCTURED) THEN
         IF (BOOL_PIC_IMPLICIT) THEN
            E = EBAR_FIELD(particles(JP)%IC, 1, :)
         ELSE
            E = E_FIELD(particles(JP)%IC, 1, :)
         END IF
      ELSE

         CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
         IF (DIMS == 2) THEN
            E = WEIGHTS(1)*E_FIELD(INDI(1), INDJ(1), :) + &
            WEIGHTS(2)*E_FIELD(INDI(2), INDJ(2), :) + &
            WEIGHTS(3)*E_FIELD(INDI(3), INDJ(3), :) + &
            WEIGHTS(4)*E_FIELD(INDI(4), INDJ(4), :)
         ELSE
            E = WEIGHTS(1)*E_FIELD(INDI(1), INDJ(1), :) + &
            WEIGHTS(2)*E_FIELD(INDI(2), INDJ(2), :)
         END IF

      END IF
   END SUBROUTINE APPLY_E_FIELD


   SUBROUTINE APPLY_POTENTIAL(JP, PHI)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(OUT) :: PHI
      INTEGER, INTENT(IN) :: JP
      REAL(KIND=8), DIMENSION(4) :: WEIGHTS
      INTEGER, DIMENSION(4) :: INDICES, INDI, INDJ

      CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
      IF (DIMS == 2) THEN
         PHI = WEIGHTS(1)*PHI_FIELD(INDI(1)+1+NPX*(INDJ(1)+1)) + &
         WEIGHTS(2)*PHI_FIELD(INDI(2)+1 +NPX*(INDJ(2)+1)) + &
         WEIGHTS(3)*PHI_FIELD(INDI(3)+1 +NPX*(INDJ(3)+1)) + &
         WEIGHTS(4)*PHI_FIELD(INDI(4)+1 +NPX*(INDJ(4)+1))
      ELSE
         PHI = WEIGHTS(1)*PHI_FIELD(INDI(1)+1 +NPX*(INDJ(1)+1)) + &
         WEIGHTS(2)*PHI_FIELD(INDI(2)+1 +NPX*(INDJ(2)+1))
      END IF
      
   END SUBROUTINE APPLY_POTENTIAL


   SUBROUTINE ST_MATRIX_PRINT(THIS, SIZE)

      TYPE(ST_MATRIX) :: THIS
      INTEGER, INTENT(IN) :: SIZE
      INTEGER :: I, J, IDX
      LOGICAL :: FOUND

      WRITE(*,*) 'Printing matrix with ', THIS%NNZ, ' non-zero entries.'
      DO I = 0, (SIZE-1)
         DO J = 0, (SIZE-1)
            FOUND = .FALSE.
            DO IDX = 0, THIS%NNZ-1
               IF (THIS%RIDX(IDX) == I .AND. THIS%CIDX(IDX) == J) THEN
                  WRITE(*,'(ES10.3)', ADVANCE = 'no') THIS%VALUE(IDX)
                  FOUND = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT. FOUND)  WRITE(*,'(ES10.3)', ADVANCE = 'no') 0.d0
         END DO
         WRITE(*,*)
      END DO

   END SUBROUTINE ST_MATRIX_PRINT

END MODULE fields
