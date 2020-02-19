! Contains Electromagnetic fields related code

MODULE fields

   USE mUMFPACK
   USE global

   IMPLICIT NONE

   CONTAINS

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ASSEMBLE_POISSON -> Prepares the linear system for the solution !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ASSEMBLE_POISSON

      INTEGER :: STATUS
      INTEGER :: I, J
      INTEGER :: IC, IS, IW, IN, IE
      INTEGER :: MAXNNZ, SIZE
      TYPE(ST_MATRIX) :: A_ST
      REAL(KIND=8) :: HX, HY, DX, DY
      REAL(KIND=8) :: PI

      PI   = 3.141593

      HX = (XMAX-XMIN)/DBLE(NX)
      HY = (YMAX-YMIN)/DBLE(NY)
      DX = 1./(HX*HX)
      DY = 1./(HY*HY)

      SIZE = NPX*NPY

      ALLOCATE(Q_FIELD(0:SIZE-1))
      ALLOCATE(Q_BC(0:SIZE-1))
      Q_BC = .FALSE.

      ! Create the matrix in Sparse Triplet format.
      ! Entries can be input in random order, but not duplicated.
      ! We have to allocate the ST matrix, but in this case
      ! it is not a problem, since each point takes at most 5 entries,
      ! we can use this as an upper limit.
      ! (We use a 5 point stencil for interior points)
      MAXNNZ = 5*SIZE
      CALL ST_MATRIX_ALLOCATE(A_ST, MAXNNZ)


      ! At this point, populate the matrix
      DO I = 0, NPX-1
         DO J = 0, NPY-1
            IC = I+(NPX)*J
            IE = I+(NPX)*J+1
            IW = I+(NPX)*J-1
            IN = I+(NPX)*(J+1)
            IS = I+(NPX)*(J-1)

            ! Simple Dirichlet on the boundary.
            IF (I == 0 .OR. I == NPX-1 .OR. J == 0 .OR. J == NPY-1) THEN
               ! Boundary point
               CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)
               Q_FIELD(IC) = 0.d0
               Q_BC(IC) = .TRUE.
            ELSE
               ! Interior point
               CALL ST_MATRIX_SET(A_ST, IC, IC, -2.*DX-2.*DY)
               CALL ST_MATRIX_SET(A_ST, IC, IN, DY)
               CALL ST_MATRIX_SET(A_ST, IC, IS, DY)
               CALL ST_MATRIX_SET(A_ST, IC, IE, DX)
               CALL ST_MATRIX_SET(A_ST, IC, IW, DX)
            END IF

            ! Example with different BCs.
            ! IF (I == 0) THEN
            !    ! IF (J == NPY-1) THEN
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IC, -1.d0)
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IE, -0.d5)
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IS, -0.d5)
            !    !    B(IC) = -HX*10.
            !    ! ELSE IF (J == 0) THEN
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IE, -0.d5)
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IN, -0.d5)
            !    !    B(IC) = -HX*10.
            !    ! ELSE
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IC, 2.d0)
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IE, -1.d0)
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IN, -0.d5)
            !    !    CALL ST_MATRIX_SET(A_ST, IC, IS, -0.d5)
            !    !    B(IC) = -HX*10. ! This is O(h**2)
            !    ! END IF
            
            !    CALL ST_MATRIX_SET(A_ST, IC, IC, -1.d0)
            !    CALL ST_MATRIX_SET(A_ST, IC, IE, 1.d0)
            !    B(IC) = HX*10.

            ! ELSE IF (I == NPX-1) THEN
            !    ! Example of Dirichlet on boundary  
            !    CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)
            !    B(IC) = 2*SIN(J/(NPY-1.)*2.*PI)
            ! ELSE IF (J == 0) THEN
            !    ! Example of Dirichlet on boundary
            !    CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)
            !    B(IC) = SIN(I/(NPX-1.)*PI)
            ! ELSE IF (J == NPY-1) THEN
            !    ! Example of Dirichlet on boundary
            !    CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)
            !    B(IC) = SIN(I/(NPX-1.)*5.*PI)
            ! ELSE IF ((I-NPX/2.)**2.+(J-3.*NPY/4.)**2 .LT. NPX/4.) THEN
            !    ! Example of Dirichlet in interior
            !    CALL ST_MATRIX_SET(A_ST, IC, IC, 1.d0)
            !    B(IC) = 3.
            ! ELSE
            !    ! Interior point
            !    CALL ST_MATRIX_SET(A_ST, IC, IC, 4.d0)
               
            !    IF (J < NPY) CALL ST_MATRIX_SET(A_ST, IC, IN, -1.d0)
            !    IF (I < NPX) CALL ST_MATRIX_SET(A_ST, IC, IE, -1.d0)
            !    IF (J > 0)   CALL ST_MATRIX_SET(A_ST, IC, IS, -1.d0)
            !    IF (I > 0)   CALL ST_MATRIX_SET(A_ST, IC, IW, -1.d0)
   
            ! END IF
            ! ! Source (RHS) in interior
            ! IF ((I-NPX/2.)**2+(J-NPY/4.)**2 .LT. NPX/4.) B(IC) = -0.1
            ! IF ((I-NPX/2.)**2+(J-NPY/4.)**2 .LT. NPX/4.) B(IC) = -0.1
         END DO
      END DO


      ! Factorize the matrix for later solution
      !IF (PROC_ID .EQ. 0) CALL ST_MATRIX_PRINT(A_ST, NPX, NPY)
      CALL ST_MATRIX_TO_CC(A_ST, SIZE, A_CC)
      CALL S_UMFPACK_SYMBOLIC(SIZE, SIZE, A_CC%AP, A_CC%AI, A_CC%AX, STATUS = STATUS)
      CALL S_UMFPACK_NUMERIC(A_CC%AP, A_CC%AI, A_CC%AX)
      CALL S_UMFPACK_FREE_SYMBOLIC

   END SUBROUTINE ASSEMBLE_POISSON



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE SOLVE_POISSON -> Solves the Poisson equation with the Q_FIELD RHS !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE SOLVE_POISSON

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X
      REAL(KIND=8) :: HX, HY
      INTEGER :: I, J

      HX = (XMAX-XMIN)/DBLE(NX)
      HY = (YMAX-YMIN)/DBLE(NY)

      ALLOCATE(X(0:NPX*NPY-1))

      ! Solve the linear system
      IF (PROC_ID .EQ. 0) THEN
         WRITE(*,*) 'Solving poisson'
         CALL S_UMFPACK_SOLVE(UMFPACK_A, A_CC%AP, A_CC%AI, A_CC%AX, X, Q_FIELD)
      END IF

      CALL MPI_BCAST(X, NPX*NPY, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      
      ! Reshape the linear array into a 2D array
      PHI_FIELD = RESHAPE(X, [NPX, NPY])
      DEALLOCATE(X)


      DO I = 0, NPX-1
         DO J = 0, NPY-1
            IF (I == 0) THEN
               E_FIELD(I,J,1) = (PHI_FIELD(I+1,J+1)-PHI_FIELD(I+2,J+1))/HX
            ELSE IF (I == NPX-1) THEN
               E_FIELD(I,J,1) = (PHI_FIELD(I,J+1)-PHI_FIELD(I+1,J+1))/HX
            ELSE
               E_FIELD(I,J,1) = 0.5*(PHI_FIELD(I,J+1)-PHI_FIELD(I+2,J+1))/HX
            END IF
            IF (J == 0) THEN
               E_FIELD(I,J,2) = (PHI_FIELD(I+1,J+1)-PHI_FIELD(I+1,J+2))/HY
            ELSE IF (J == NPY-1) THEN
               E_FIELD(I,J,2) = (PHI_FIELD(I+1,J)-PHI_FIELD(I+1,J+1))/HY
            ELSE
               E_FIELD(I,J,2) = 0.5*(PHI_FIELD(I+1,J)-PHI_FIELD(I+1,J+2))/HY
            END IF
            E_FIELD(I,J,3) = 0.d0
         END DO
      END DO

   END SUBROUTINE SOLVE_POISSON


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DEPOSIT_CHARGE -> Deposits the charge of particles on grid points !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DEPOSIT_CHARGE
      
      IMPLICIT NONE

      INTEGER :: JP

      REAL(KIND=8) :: K, RHO_Q, CHARGE
      REAL(KIND=8) :: HX, HY
      REAL(KIND=8), DIMENSION(4) :: WEIGHTS
      INTEGER, DIMENSION(4) :: INDICES, INDI, INDJ

      K = 1.8095128E-8 ! [V m] Elementary charge / Dielectric constant of vacuum

      HX = (XMAX-XMIN)/DBLE(NX)
      HY = (YMAX-YMIN)/DBLE(NY)

      Q_FIELD = 0

      DO JP = 1, NP_PROC
         CHARGE = SPECIES(particles(JP)%S_ID)%CHARGE
         IF (CHARGE .EQ. 0.d0) CYCLE

         CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
         RHO_Q = -K*CHARGE*FNUM/(HX*HY)

         IF (Q_BC(INDICES(3)) .EQV. .FALSE.) THEN
            Q_FIELD(INDICES(3)) = Q_FIELD(INDICES(3)) + RHO_Q * WEIGHTS(3)
         END IF
         IF (Q_BC(INDICES(4)) .EQV. .FALSE.) THEN
            Q_FIELD(INDICES(4)) = Q_FIELD(INDICES(4)) + RHO_Q * WEIGHTS(4)
         END IF
         IF (Q_BC(INDICES(2)) .EQV. .FALSE.) THEN
            Q_FIELD(INDICES(2)) = Q_FIELD(INDICES(2)) + RHO_Q * WEIGHTS(2)
         END IF
         IF (Q_BC(INDICES(1)) .EQV. .FALSE.) THEN
            Q_FIELD(INDICES(1)) = Q_FIELD(INDICES(1)) + RHO_Q * WEIGHTS(1)
         END IF

      END DO


      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, Q_FIELD, NPX*NPY, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ELSE
         CALL MPI_REDUCE(Q_FIELD,      Q_FIELD, NPX*NPY, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF

      ! Q_FIELD = 0
      ! ALLOCATE(Q_TEMP(NPX,NPY))
      ! Q_TEMP = 0
      ! DO I = NPX/4, 3*NPX/4
      !    DO J = NPY/4, 3*NPY/4
      !       Q_TEMP(I,J) = 1.d0
      !    END DO
      ! END DO
      ! Q_FIELD = PACK(Q_TEMP, .TRUE.)

   END SUBROUTINE DEPOSIT_CHARGE


   SUBROUTINE COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: JP
      REAL(KIND=8), DIMENSION(4), INTENT(OUT) :: WEIGHTS
      INTEGER, DIMENSION(4), INTENT(OUT) :: INDICES, INDI, INDJ
      REAL(KIND=8) :: XP, YP, CELL_XMIN, CELL_YMIN, FX, FY, HX, HY
      INTEGER :: IC, JCOL, JROW

      XP = particles(JP)%X
      YP = particles(JP)%Y
      ! This is where we need the cell properties based on particle position
      ! With rectilinear non-uniform grid this has to be changed.
      IC = particles(JP)%IC
      HX = (XMAX-XMIN)/DBLE(NX)
      HY = (YMAX-YMIN)/DBLE(NY)
      JROW = MOD(IC-1, NX)
      JCOL = (IC-1)/NX
      CELL_XMIN = JROW * HX + XMIN
      CELL_YMIN = JCOL * HY + YMIN
      FX = (particles(JP)%X - CELL_XMIN)/HX
      FY = (particles(JP)%Y - CELL_YMIN)/HY


      ! Row and column indices, starting from 0
      INDI(1) = JROW
      INDI(2) = JROW
      INDI(3) = JROW + 1
      INDI(4) = JROW + 1

      INDJ(1) = JCOL
      INDJ(2) = JCOL
      INDJ(3) = JCOL + 1
      INDJ(4) = JCOL + 1

      ! CCW starting from SW
      INDICES(1) = IC+JCOL-1
      INDICES(2) = INDICES(1) + 1
      INDICES(4) = INDICES(1) + NX + 1
      INDICES(3) = INDICES(4) + 1

      ! CCW starting from SW
      WEIGHTS(1) = (1.-FX)*(1.-FY)
      WEIGHTS(2) = FX*(1.-FY)
      WEIGHTS(3) = FX*FY
      WEIGHTS(4) = (1.-FX)*FY

   END SUBROUTINE COMPUTE_WEIGHTS


   SUBROUTINE APPLY_E_FIELD(JP, E)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: E
      INTEGER, INTENT(IN) :: JP
      REAL(KIND=8), DIMENSION(4) :: WEIGHTS
      INTEGER, DIMENSION(4) :: INDICES, INDI, INDJ

      CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
      E = WEIGHTS(1)*E_FIELD(INDI(1), INDJ(1), :) + &
      WEIGHTS(2)*E_FIELD(INDI(2), INDJ(2), :) + &
      WEIGHTS(3)*E_FIELD(INDI(3), INDJ(3), :) + &
      WEIGHTS(4)*E_FIELD(INDI(4), INDJ(4), :)

      
   END SUBROUTINE APPLY_E_FIELD

   SUBROUTINE ST_MATRIX_PRINT(THIS, NROWS, NCOLS)

      TYPE(ST_MATRIX) :: THIS
      INTEGER, INTENT(IN) :: NROWS, NCOLS
      INTEGER :: I, J, IDX
      LOGICAL :: FOUND

      WRITE(*,*) 'Printing matrix with ', THIS%NNZ, ' non-zero entries.'
      DO I = 0, (NROWS*NCOLS-1)
         DO J = 0, (NROWS*NCOLS-1)
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