! Contains Electromagnetic fields related code

MODULE fields

#include <petsc/finclude/petscksp.h>

   USE mUMFPACK
   USE global
   USE screen
   USE tools


   
   USE petscksp

   IMPLICIT NONE

   Mat Amat
   Vec bvec, xvec, X_SEQ
   VecScatter ctx
   KSP ksp
   PetscInt one,f9,f6,f30
   PetscInt Istart, Iend
   PetscReal val
   PetscScalar, POINTER :: PHI_FIELD(:)
   PetscScalar, POINTER :: PHIBAR_FIELD(:)
   KSPConvergedReason reason
   PC pc

   CONTAINS
   


   SUBROUTINE PETSC_INITIAL_TEST

      IMPLICIT NONE
!
!  This example demonstrates basic use of the PETSc Fortran interface
!  to vectors.
!
      PetscInt  n
      PetscErrorCode ierr
      PetscBool  flg
      PetscScalar      one,two,three,dot
      PetscReal        norm,rdot
      Vec              x,y,w
      PetscOptions     options

      n     = 20
      one   = 1.0
      two   = 2.0
      three = 3.0

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'Unable to initialize PETSc'
        stop
      endif
      call PetscOptionsCreate(options,ierr)
      call PetscOptionsGetInt(options,PETSC_NULL_CHARACTER,                  &
    &                        '-n',n,flg,ierr)
      call PetscOptionsDestroy(options,ierr)

! Create a vector, then duplicate it
      call VecCreate(PETSC_COMM_WORLD,x,ierr)
      call VecSetSizes(x,PETSC_DECIDE,n,ierr)
      call VecSetFromOptions(x,ierr)
      call VecDuplicate(x,y,ierr)
      call VecDuplicate(x,w,ierr)

      call VecSet(x,one,ierr)
      call VecSet(y,two,ierr)

      call VecDot(x,y,dot,ierr)
      rdot = PetscRealPart(dot)
      write(6,100) rdot
 100  format('Result of inner product ',f10.4)

      call VecScale(x,two,ierr)
      call VecNorm(x,NORM_2,norm,ierr)
      write(6,110) norm
 110  format('Result of scaling ',f10.4)

      call VecCopy(x,w,ierr)
      call VecNorm(w,NORM_2,norm,ierr)
      write(6,120) norm
 120  format('Result of copy ',f10.4)

      call VecAXPY(y,three,x,ierr)
      call VecNorm(y,NORM_2,norm,ierr)
      write(6,130) norm
 130  format('Result of axpy ',f10.4)

      call VecDestroy(x,ierr)
      call VecDestroy(y,ierr)
      call VecDestroy(w,ierr)
      call PetscFinalize(ierr)



   END SUBROUTINE PETSC_INITIAL_TEST


   SUBROUTINE PETSC_INIT

      CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      !PetscCallMPIA(MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)) ---> N_MPI_THREADS
      !PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)) ---> PROC_ID

   END SUBROUTINE PETSC_INIT


   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ASSEMBLE_POISSON -> Prepares the linear system for the solution  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ASSEMBLE_POISSON

      IMPLICIT NONE

      INTEGER :: STATUS
      INTEGER :: I, J
      INTEGER :: ICENTER, INORTH, IEAST, ISOUTH, IWEST
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

      !IF (.NOT. ALLOCATED(PHI_FIELD)) ALLOCATE(PHI_FIELD(0:SIZE-1))
      !PHI_FIELD = 0.d0

      IF (.NOT. ALLOCATED(DIRICHLET)) ALLOCATE(DIRICHLET(0:SIZE-1))
      IF (.NOT. ALLOCATED(IS_DIRICHLET)) ALLOCATE(IS_DIRICHLET(0:SIZE-1))
      IS_DIRICHLET = .FALSE.

      IF (.NOT. ALLOCATED(NEUMANN)) ALLOCATE(NEUMANN(0:SIZE-1))
      NEUMANN = 0.d0
      IF (.NOT. ALLOCATED(IS_NEUMANN)) ALLOCATE(IS_NEUMANN(0:SIZE-1))
      IS_NEUMANN = .FALSE.

      f6 = 6
      f9 = 9   
      one = 1
      f30 = 30

      CALL MatCreate(PETSC_COMM_WORLD,Amat,ierr)
      CALL MatSetSizes( Amat,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
      CALL MatSetType( Amat, MATAIJ, ierr)
      CALL MatSetOption(Amat,MAT_SPD,PETSC_TRUE,ierr)
      IF (N_MPI_THREADS == 1) THEN
         CALL MatSetType( Amat, MATAIJ, ierr)
      ELSE
         CALL MatSetType( Amat, MATMPIAIJ, ierr)
      END IF
      CALL MatMPIAIJSetPreallocation(Amat,f30,PETSC_NULL_INTEGER,f30,PETSC_NULL_INTEGER, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
      CALL MatSetFromOptions( Amat, ierr)
      CALL MatSetUp( Amat, ierr)
      CALL MatGetOwnershipRange( Amat, Istart, Iend, ierr)

      CALL MatCreateVecs( Amat, PETSC_NULL_VEC, xvec, ierr)
      CALL VecSetFromOptions( xvec, ierr)
      CALL VecDuplicate( xvec, bvec, ierr)

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
            IF (AXI) THEN
               K11 = K11*(Y1+Y2+Y3)/3.
               K22 = K22*(Y1+Y2+Y3)/3.
               K33 = K33*(Y1+Y2+Y3)/3.
               K12 = K12*(Y1+Y2+Y3)/3.
               K23 = K23*(Y1+Y2+Y3)/3.
               K13 = K13*(Y1+Y2+Y3)/3.
            END IF

            ! We need to ADD to a sparse matrix entry.
            IF (V1-1 >= Istart .AND. V1-1 < Iend) THEN
               IF (IS_DIRICHLET(V1-1)) THEN
                  CALL MatSetValues(Amat,one,V1-1,one,V1-1,1.d0,ADD_VALUES,ierr)
                  !CALL ST_MATRIX_SET(A_ST, V1-1, V1-1, 1.d0)
               ELSE !IF (.NOT. IS_NEUMANN(V1-1)) THEN
                  CALL MatSetValues(Amat,one,V1-1,one,V1-1,K11,ADD_VALUES,ierr)
                  CALL MatSetValues(Amat,one,V1-1,one,V2-1,K12,ADD_VALUES,ierr)
                  CALL MatSetValues(Amat,one,V1-1,one,V3-1,K13,ADD_VALUES,ierr)

                  !CALL ST_MATRIX_ADD(A_ST, V1-1, V1-1, K11)
                  !CALL ST_MATRIX_ADD(A_ST, V1-1, V2-1, K12)
                  !CALL ST_MATRIX_ADD(A_ST, V1-1, V3-1, K13)
               END IF
            END IF
            IF (V2-1 >= Istart .AND. V2-1 < Iend) THEN
               IF (IS_DIRICHLET(V2-1)) THEN
                  CALL MatSetValues(Amat,one,V2-1,one,V2-1,1.d0,ADD_VALUES,ierr)
                  !CALL ST_MATRIX_SET(A_ST, V2-1, V2-1, 1.d0)
               ELSE !IF (.NOT. IS_NEUMANN(V2-1)) THEN
                  CALL MatSetValues(Amat,one,V2-1,one,V1-1,K12,ADD_VALUES,ierr)
                  CALL MatSetValues(Amat,one,V2-1,one,V3-1,K23,ADD_VALUES,ierr)
                  CALL MatSetValues(Amat,one,V2-1,one,V2-1,K22,ADD_VALUES,ierr)
                  
                  !CALL ST_MATRIX_ADD(A_ST, V2-1, V1-1, K12)
                  !CALL ST_MATRIX_ADD(A_ST, V2-1, V3-1, K23)
                  !CALL ST_MATRIX_ADD(A_ST, V2-1, V2-1, K22)
               END IF
            END IF
            IF (V3-1 >= Istart .AND. V3-1 < Iend) THEN
               IF (IS_DIRICHLET(V3-1)) THEN
                  CALL MatSetValues(Amat,one,V3-1,one,V3-1,1.d0,ADD_VALUES,ierr)
                  !CALL ST_MATRIX_SET(A_ST, V3-1, V3-1, 1.d0)
               ELSE !IF (.NOT. IS_NEUMANN(V3-1)) THEN
                  CALL MatSetValues(Amat,one,V3-1,one,V1-1,K13,ADD_VALUES,ierr)
                  CALL MatSetValues(Amat,one,V3-1,one,V2-1,K23,ADD_VALUES,ierr)
                  CALL MatSetValues(Amat,one,V3-1,one,V3-1,K33,ADD_VALUES,ierr)
                  
                  !CALL ST_MATRIX_ADD(A_ST, V3-1, V1-1, K13)
                  !CALL ST_MATRIX_ADD(A_ST, V3-1, V2-1, K23)
                  !CALL ST_MATRIX_ADD(A_ST, V3-1, V3-1, K33)
               END IF
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
               CALL MatSetValues(Amat,one,I-1,one,I-1,1.d0,ADD_VALUES,ierr)
               IS_DIRICHLET(I-1) = .TRUE.
               DIRICHLET(I-1) = 0.d0
            END IF
         END DO
         DEALLOCATE(IS_UNUSED)

      ELSE
         DO I = 0, NPX-1
            DO J = 0, NPY-1


               ICENTER = I+(NPX)*J
               IF (ICENTER < Istart .OR. ICENTER >= Iend) CYCLE ! Each proc populates part of the matrix.
               IEAST = I+(NPX)*J+1
               IWEST = I+(NPX)*J-1
               INORTH = I+(NPX)*(J+1)
               ISOUTH = I+(NPX)*(J-1)

               HX = (XMAX-XMIN)/DBLE(NX)
               AX = 1./(HX*HX)


               ! Various configurations go here.
               
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Plume 2d cartesian, full domain !
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               IF ((I == 0) .OR. (I == NPX-1) .OR. (J == 0) .OR. (J == NPY-1)) THEN
                  ! Boundary point
                  CALL MatSetValues(Amat,one,ICENTER,one,ICENTER,1.d0,INSERT_VALUES,ierr)
                  ! On the boundary or on the rest of the PFG.
                  DIRICHLET(ICENTER) = 0.d0
                  IS_DIRICHLET(ICENTER) = .TRUE.


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
                     CALL MatSetValues(Amat,one,ICENTER,one,ICENTER,BX+BY,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,INORTH,CY,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,ISOUTH,AY,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,IEAST,CX,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,IWEST,AX,INSERT_VALUES,ierr)
                  ELSE
                     CALL MatSetValues(Amat,one,ICENTER,one,ICENTER,BX,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,IEAST,CX,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,IWEST,AX,INSERT_VALUES,ierr)
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

      CALL MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)

      ! IF (PROC_ID .EQ. 0) THEN
      !    CALL ST_MATRIX_TO_CC(A_ST, SIZE, A_CC)
      !    CALL S_UMFPACK_SYMBOLIC(SIZE, SIZE, A_CC%AP, A_CC%AI, A_CC%AX, STATUS = STATUS)
      !    WRITE(*,*) 'Status is = ', STATUS
      !    CALL S_UMFPACK_NUMERIC(A_CC%AP, A_CC%AI, A_CC%AX, STATUS = STATUS)
      !    WRITE(*,*) 'Status 1 is = ', STATUS
      !    CALL S_UMFPACK_FREE_SYMBOLIC
      ! END IF

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
      INTEGER :: ICENTER, INORTH, ISOUTH, IEAST, IWEST

      HX = (XMAX-XMIN)/DBLE(NX)
      HY = (YMAX-YMIN)/DBLE(NY)


      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZE = U2D_GRID%NUM_NODES
      ELSE
         SIZE = NPX*NPY
      END IF

      ! Solve the linear system

      CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
      CALL KSPSetOperators(ksp,Amat,Amat,ierr)

      ! Set solver type to direct (LU). Could be MUMPS but must install with --download-mumps
      !CALL KSPGetPC(ksp,pc,ierr)
      !CALL PCSetType(pc,PCLU,ierr)
      !CALL KSPSetFromOptions(ksp,ierr)

      CALL KSPSolve(ksp,bvec,xvec,ierr)

      CALL KSPGetConvergedReason(ksp,reason,ierr)
      IF (PROC_ID == 0) WRITE(*,*) 'KSPConvergedReason = ', reason

      CALL VecScatterCreateToAll(xvec,ctx,X_SEQ,ierr)
      CALL VecScatterBegin(ctx,xvec,X_SEQ,INSERT_VALUES,SCATTER_FORWARD)
      CALL VecScatterEnd(ctx,xvec,X_SEQ,INSERT_VALUES,SCATTER_FORWARD)

      CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD,ierr)

      PHIBAR_FIELD => PHI_FIELD

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

            E_FIELD(I,1,1) = -0.5/AREA*(  PHI_FIELD(V1)*(Y2-Y3) &
                                        - PHI_FIELD(V2)*(Y1-Y3) &
                                        - PHI_FIELD(V3)*(Y2-Y1))
            E_FIELD(I,1,2) = -0.5/AREA*(- PHI_FIELD(V1)*(X2-X3) &
                                        + PHI_FIELD(V2)*(X1-X3) &
                                        + PHI_FIELD(V3)*(X2-X1))
            E_FIELD(I,1,3) = 0.d0
         END DO

      ELSE
         ! Compute the electric field at grid points
         DO I = 0, NPX-1
            DO J = 0, NPY-1

               ICENTER = I+1 + NPX*J
               IEAST = I+2 + NPX*J
               IWEST = I   + NPX*J
               INORTH = I+1 + NPX*(J+1)
               ISOUTH = I+1 + NPX*(J-1)

               IF (I == 0) THEN ! Left boundary
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = XSIZE(1)
                  E_FIELD(I,J,1) = (PHI_FIELD(ICENTER)-PHI_FIELD(IEAST))/HX
               ELSE IF (I == NPX-1) THEN ! Right boundary
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = XSIZE(NPX-1)
                  E_FIELD(I,J,1) = (PHI_FIELD(IWEST)-PHI_FIELD(ICENTER))/HX
               ELSE ! Interior point
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = 0.5*(XSIZE(I)+XSIZE(I+1))
                  E_FIELD(I,J,1) = 0.5*(PHI_FIELD(IWEST)-PHI_FIELD(IEAST))/HX
               END IF
               IF (DIMS == 2) THEN
                  IF (J == 0) THEN ! Bottom boundary
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = YSIZE(1)
                     E_FIELD(I,J,2) = (PHI_FIELD(ICENTER)-PHI_FIELD(INORTH))/HY
                  ELSE IF (J == NPY-1) THEN ! Top boundary
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = YSIZE(NPY-1)
                     E_FIELD(I,J,2) = (PHI_FIELD(ISOUTH)-PHI_FIELD(ICENTER))/HY
                  ELSE ! Interior point
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = 0.5*(YSIZE(J)+YSIZE(J+1))
                     E_FIELD(I,J,2) = 0.5*(PHI_FIELD(ISOUTH)-PHI_FIELD(INORTH))/HY
                  END IF
               ELSE
                  E_FIELD(I,J,2) = 0.d0
               END IF
               E_FIELD(I,J,3) = 0.d0
            END DO
         END DO
         
      END IF

      !CALL VecRestoreArrayReadF90(X_SEQ,PHI_FIELD,ierr)

   END SUBROUTINE SOLVE_POISSON


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DEPOSIT_CHARGE -> Deposits the charge of particles on grid points !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DEPOSIT_CHARGE
      
      IMPLICIT NONE

      INTEGER :: JP, I, IC

      REAL(KIND=8) :: K, RHO_Q, CHARGE
      REAL(KIND=8) :: VOL, CFNUM, AREA
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
            AREA = CELL_AREAS(IC)
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
            PSI1 = 0.5*( (Y2-Y3)*(XP-X3) - (X2-X3)*(YP-Y3))/AREA/(ZMAX-ZMIN)
            PSI2 = 0.5*(-(Y1-Y3)*(XP-X3) + (X1-X3)*(YP-Y3))/AREA/(ZMAX-ZMIN)
            PSI3 = 0.5*(-(Y2-Y1)*(XP-X1) + (X2-X1)*(YP-Y1))/AREA/(ZMAX-ZMIN)
            
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

      CALL MPI_BCAST(RHS, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)


      DO I = Istart, Iend-1
         IF (IS_DIRICHLET(I)) THEN
            val = DIRICHLET(I)
         ELSE IF (IS_NEUMANN(I)) THEN
            val = RHS(I) + NEUMANN(I)
         ELSE
            val = RHS(I)
         END IF

         CALL VecSetValues(bvec,one,I,val,INSERT_VALUES,ierr)
      END DO


      CALL VecAssemblyBegin(bvec,ierr)
      CALL VecAssemblyEnd(bvec,ierr)


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

      SIZE = U2D_GRID%NUM_NODES

      IF (PROC_ID .EQ. 0) THEN

         ALLOCATE(PHIK, SOURCE = PHI_FIELD)
         ALLOCATE(BK(0:SIZE-1))
         ALLOCATE(DBDPHI(0:SIZE-1))
         ALLOCATE(YK(0:SIZE-1))
         IF (.NOT. ALLOCATED(BOLTZ_NRHOE)) ALLOCATE(BOLTZ_NRHOE(SIZE))

         DO ITERATION = 1, 10
            WRITE(*,*) '---> Starting iteration ', ITERATION
            
            !WRITE(*,*) PHIK
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
                     !DBDPHI(IN-1) = DBDPHI(IN-1) + QE*QE/(EPS0*KB*BOLTZ_TE)*CELL_AREAS(IC)/3*BOLTZ_N0
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

      END IF

      CALL MPI_BCAST(PHI_FIELD, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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


   END SUBROUTINE SOLVE_POISSON_NONLINEAR


   SUBROUTINE ASSEMBLE_AMPERE

      INTEGER :: STATUS
      INTEGER :: I, J
      INTEGER :: MAXNNZ, SIZE
      TYPE(ST_MATRIX) :: A_ST
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, K11, K22, K33, K12, K23, K13, AREA
      REAL(KIND=8) :: K11TILDE, K22TILDE, K33TILDE, K12TILDE, K23TILDE, K13TILDE
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

      IF (.NOT. ALLOCATED(DIRICHLET)) ALLOCATE(DIRICHLET(0:SIZE-1))
      IF (.NOT. ALLOCATED(IS_DIRICHLET)) ALLOCATE(IS_DIRICHLET(0:SIZE-1))
      IS_DIRICHLET = .FALSE.

      IF (.NOT. ALLOCATED(NEUMANN)) ALLOCATE(NEUMANN(0:SIZE-1))
      NEUMANN = 0.d0
      IF (.NOT. ALLOCATED(IS_NEUMANN)) ALLOCATE(IS_NEUMANN(0:SIZE-1))
      IS_NEUMANN = .FALSE.


      ! Create the matrix in Sparse Triplet format.

      CALL MatDestroy(Amat,ierr)
      CALL MatCreate(PETSC_COMM_WORLD,Amat,ierr)
      CALL MatSetSizes( Amat,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
      CALL MatSetType( Amat, MATAIJ, ierr)
      CALL MatSetOption(Amat,MAT_SPD,PETSC_TRUE,ierr)
      IF (N_MPI_THREADS == 1) THEN
         CALL MatSetType( Amat, MATAIJ, ierr)
      ELSE
         CALL MatSetType( Amat, MATMPIAIJ, ierr)
      END IF
      CALL MatMPIAIJSetPreallocation(Amat,f9,PETSC_NULL_INTEGER,f6,PETSC_NULL_INTEGER, ierr)
      CALL MatSetFromOptions( Amat, ierr)
      CALL MatSetUp( Amat, ierr)
      CALL MatGetOwnershipRange( Amat, Istart, Iend, ierr)

      CALL MatCreateVecs( Amat, PETSC_NULL_VEC, xvec, ierr)
      CALL VecSetFromOptions( xvec, ierr)
      CALL VecDuplicate( xvec, bvec, ierr)


      ! At this point, populate the matrix
      IF (GRID_TYPE == UNSTRUCTURED) THEN
         ALLOCATE(IS_UNUSED(0:SIZE-1))
         IS_UNUSED = .TRUE.

         DO I = 1, U2D_GRID%NUM_CELLS
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)

            IF ((V1-1 < Istart .OR. V1-1 > Iend) .AND. &
                (V2-1 < Istart .OR. V2-1 > Iend) .AND. &
                (V3-1 < Istart .OR. V3-1 > Iend)) CYCLE

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
            V1 = U2D_GRID%CELL_NODES(I,1)
            V2 = U2D_GRID%CELL_NODES(I,2)
            V3 = U2D_GRID%CELL_NODES(I,3)
            
            IF ((V1-1 < Istart .OR. V1-1 > Iend) .AND. &
                (V2-1 < Istart .OR. V2-1 > Iend) .AND. &
                (V3-1 < Istart .OR. V3-1 > Iend)) CYCLE

            AREA = CELL_AREAS(I)
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
            K11TILDE = K11*(Y1+Y2+Y3)/3.
            K22TILDE = K22*(Y1+Y2+Y3)/3.
            K33TILDE = K33*(Y1+Y2+Y3)/3.
            K12TILDE = K12*(Y1+Y2+Y3)/3.
            K23TILDE = K23*(Y1+Y2+Y3)/3.
            K13TILDE = K13*(Y1+Y2+Y3)/3.

            ! We need to ADD to a sparse matrix entry.
            IF (V1-1 >= Istart .AND. V1-1 <= Iend) THEN
               IF (IS_DIRICHLET(V1-1)) THEN
                  CALL MatSetValues(Amat,one,V1-1,one,V1-1,1.d0,INSERT_VALUES,ierr)
               ELSE !IF (.NOT. IS_NEUMANN(V1-1)) THEN
                  IF (AXI) THEN
                     CALL MatSetValues(Amat,one,V1-1,one,V1-1,K11TILDE + MASS_MATRIX(I)*K11,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V1-1,one,V2-1,K12TILDE + MASS_MATRIX(I)*K12,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V1-1,one,V3-1,K13TILDE + MASS_MATRIX(I)*K13,ADD_VALUES,ierr)
                     val = PHI_FIELD(V1-1)*K11TILDE+PHI_FIELD(V2-1)*K12TILDE+PHI_FIELD(V3-1)*K13TILDE
                     CALL VecSetValues(bvec,one,V1-1,val,ADD_VALUES,ierr)
                  ELSE
                     CALL MatSetValues(Amat,one,V1-1,one,V1-1,(MASS_MATRIX(I)+1.)*K11,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V1-1,one,V2-1,(MASS_MATRIX(I)+1.)*K12,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V1-1,one,V3-1,(MASS_MATRIX(I)+1.)*K13,ADD_VALUES,ierr)
                     val = PHI_FIELD(V1-1)*K11+PHI_FIELD(V2-1)*K12+PHI_FIELD(V3-1)*K13
                     CALL VecSetValues(bvec,one,V1-1,val,ADD_VALUES,ierr)
                  END IF
               END IF
            END IF
            IF (V2-1 >= Istart .AND. V2-1 <= Iend) THEN
               IF (IS_DIRICHLET(V2-1)) THEN
                  CALL MatSetValues(Amat,one,V2-1,one,V2-1,1.d0,INSERT_VALUES,ierr)
               ELSE !IF (.NOT. IS_NEUMANN(V2-1)) THEN
                  IF (AXI) THEN
                     CALL MatSetValues(Amat,one,V2-1,one,V1-1,K12TILDE + MASS_MATRIX(I)*K12,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V2-1,one,V3-1,K23TILDE + MASS_MATRIX(I)*K23,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V2-1,one,V2-1,K22TILDE + MASS_MATRIX(I)*K22,ADD_VALUES,ierr)
                     val = PHI_FIELD(V1-1)*K12TILDE+PHI_FIELD(V2-1)*K22TILDE+PHI_FIELD(V3-1)*K23TILDE
                     CALL VecSetValues(bvec,one,V2-1,val,ADD_VALUES,ierr)
                  ELSE
                     CALL MatSetValues(Amat,one,V2-1,one,V1-1,(MASS_MATRIX(I)+1.)*K12,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V2-1,one,V3-1,(MASS_MATRIX(I)+1.)*K23,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V2-1,one,V2-1,(MASS_MATRIX(I)+1.)*K22,ADD_VALUES,ierr)
                     val = PHI_FIELD(V1-1)*K12+PHI_FIELD(V2-1)*K22+PHI_FIELD(V3-1)*K23
                     CALL VecSetValues(bvec,one,V2-1,val,ADD_VALUES,ierr)
                  END IF
               END IF
            END IF
            IF (V3-1 >= Istart .AND. V3-1 <= Iend) THEN
               IF (IS_DIRICHLET(V3-1)) THEN
                  CALL MatSetValues(Amat,one,V3-1,one,V3-1,1.d0,INSERT_VALUES,ierr)
               ELSE !IF (.NOT. IS_NEUMANN(V3-1)) THEN
                  IF (AXI) THEN
                     CALL MatSetValues(Amat,one,V3-1,one,V1-1,K13TILDE + MASS_MATRIX(I)*K13,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V3-1,one,V2-1,K23TILDE + MASS_MATRIX(I)*K23,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V3-1,one,V3-1,K33TILDE + MASS_MATRIX(I)*K33,ADD_VALUES,ierr)
                     val = PHI_FIELD(V1-1)*K13TILDE+PHI_FIELD(V2-1)*K23TILDE+PHI_FIELD(V3-1)*K33TILDE
                     CALL VecSetValues(bvec,one,V3-1,val,ADD_VALUES,ierr)
                  ELSE

                     CALL MatSetValues(Amat,one,V3-1,one,V1-1,(MASS_MATRIX(I)+1.)*K13,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V3-1,one,V2-1,(MASS_MATRIX(I)+1.)*K23,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V3-1,one,V3-1,(MASS_MATRIX(I)+1.)*K33,ADD_VALUES,ierr)
                     val = PHI_FIELD(V1-1)*K13+PHI_FIELD(V2-1)*K23+PHI_FIELD(V3-1)*K33
                     CALL VecSetValues(bvec,one,V3-1,val,ADD_VALUES,ierr)
                  END IF
               END IF
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

         DO I = Istart, Iend-1
            IF (IS_UNUSED(I)) THEN
               CALL MatSetValues(Amat,one,I,one,I,1.d0,INSERT_VALUES,ierr)
               IS_DIRICHLET(I) = .TRUE.
               DIRICHLET(I) = 0.d0
            ELSE IF (.NOT. IS_DIRICHLET(I) ) THEN
               val = 0.5*DT/EPS0*J_FIELD(I)
               CALL VecSetValues(bvec,one,I,val,ADD_VALUES,ierr)
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

      DO I = Istart, Iend-1
         IF (IS_DIRICHLET(I)) THEN
            val = DIRICHLET(I)
            CALL VecSetValues(bvec,one,I,val,INSERT_VALUES,ierr)
         ELSE IF (IS_NEUMANN(I)) THEN
            val = NEUMANN(I)
            CALL VecSetValues(bvec,one,I,val,ADD_VALUES,ierr)
         END IF
      END DO


      CALL MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
      CALL VecAssemblyBegin(bvec,ierr)
      CALL VecAssemblyEnd(bvec,ierr)

   END SUBROUTINE ASSEMBLE_AMPERE


   SUBROUTINE DEPOSIT_CURRENT
      
      IMPLICIT NONE

      INTEGER :: JP, IC

      REAL(KIND=8) :: CHARGE
      REAL(KIND=8) :: AREA
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
            AREA = CELL_AREAS(IC)
            V1 = U2D_GRID%CELL_NODES(IC,1)
            V2 = U2D_GRID%CELL_NODES(IC,2)
            V3 = U2D_GRID%CELL_NODES(IC,3)            
            X1 = U2D_GRID%NODE_COORDS(V1, 1)
            X2 = U2D_GRID%NODE_COORDS(V2, 1)
            X3 = U2D_GRID%NODE_COORDS(V3, 1)
            Y1 = U2D_GRID%NODE_COORDS(V1, 2)
            Y2 = U2D_GRID%NODE_COORDS(V2, 2)
            Y3 = U2D_GRID%NODE_COORDS(V3, 2)
            DPSI1DX =  0.5*(Y2-Y3)/AREA
            DPSI2DX = -0.5*(Y1-Y3)/AREA
            DPSI3DX = -0.5*(Y2-Y1)/AREA
            DPSI1DY = -0.5*(X2-X3)/AREA
            DPSI2DY =  0.5*(X1-X3)/AREA
            DPSI3DY =  0.5*(X2-X1)/AREA
            
            IF (AXI) THEN
               J_FIELD(V1-1) = J_FIELD(V1-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI1DX + particles(JP)%VY*DPSI1DY)/(ZMAX-ZMIN)
               J_FIELD(V2-1) = J_FIELD(V2-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI2DX + particles(JP)%VY*DPSI2DY)/(ZMAX-ZMIN)
               J_FIELD(V3-1) = J_FIELD(V3-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI3DX + particles(JP)%VY*DPSI3DY)/(ZMAX-ZMIN)

               MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*particles(JP)%DTRIM/EPS0/AREA/(ZMAX-ZMIN)*FNUM &
                                 * (QE*CHARGE)**2/SPECIES(particles(JP)%S_ID)%MOLECULAR_MASS
            ELSE
               J_FIELD(V1-1) = J_FIELD(V1-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI1DX + particles(JP)%VY*DPSI1DY)/(ZMAX-ZMIN)
               J_FIELD(V2-1) = J_FIELD(V2-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI2DX + particles(JP)%VY*DPSI2DY)/(ZMAX-ZMIN)
               J_FIELD(V3-1) = J_FIELD(V3-1) + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI3DX + particles(JP)%VY*DPSI3DY)/(ZMAX-ZMIN)

               MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*particles(JP)%DTRIM/EPS0/AREA/(ZMAX-ZMIN)*FNUM &
                                 * (QE*CHARGE)**2/SPECIES(particles(JP)%S_ID)%MOLECULAR_MASS
            END IF
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


      CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
      CALL KSPSetOperators(ksp,Amat,Amat,ierr)

      CALL KSPSolve(ksp,bvec,xvec,ierr)

      CALL VecScatterCreateToAll(xvec,ctx,X_SEQ,ierr)
      CALL VecScatterBegin(ctx,xvec,X_SEQ,INSERT_VALUES,SCATTER_FORWARD)
      CALL VecScatterEnd(ctx,xvec,X_SEQ,INSERT_VALUES,SCATTER_FORWARD)

      CALL VecGetArrayF90(X_SEQ,PHIBAR_FIELD,ierr)

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



   SUBROUTINE COMPUTE_B_FIELD_FROM_SOLENOIDS()

      IMPLICIT NONE

      REAL(KIND=8) :: DTHETA, THETA, WIRER, WIREZ
      REAL(KIND=8), DIMENSION(3) :: POINT, DL, RPRIME
      INTEGER :: ICOIL, IN, IX, IY, ITHETA, NTHETA

      NTHETA = 100

      DTHETA = 2*PI/REAL(NTHETA)

      IF (GRID_TYPE == UNSTRUCTURED) THEN

         DO ICOIL = 1, N_SOLENOIDS
            DO IN = 1, U2D_GRID%NUM_NODES
               POINT = U2D_GRID%NODE_COORDS(IN, :)
               DO IX = 1, SOLENOIDS(ICOIL)%N_WIRES_X
                  
                  IF (SOLENOIDS(ICOIL)%N_WIRES_X == 1) THEN
                     WIREZ = SOLENOIDS(ICOIL)%X1
                  ELSE
                     WIREZ = SOLENOIDS(ICOIL)%X1 + (SOLENOIDS(ICOIL)%X2-SOLENOIDS(ICOIL)%X1)/ &
                     REAL(SOLENOIDS(ICOIL)%N_WIRES_X-1)*REAL(IX-1)
                  END IF

                  DO IY = 1, SOLENOIDS(ICOIL)%N_WIRES_Y
                     
                     IF (SOLENOIDS(ICOIL)%N_WIRES_Y == 1) THEN
                        WIRER = SOLENOIDS(ICOIL)%Y1
                     ELSE
                        WIRER = SOLENOIDS(ICOIL)%Y1 + (SOLENOIDS(ICOIL)%Y2-SOLENOIDS(ICOIL)%Y1)/ &
                        REAL(SOLENOIDS(ICOIL)%N_WIRES_Y-1)*REAL(IY-1)
                     END IF

                     DO ITHETA = 1, NTHETA
                        
                        THETA = 2*PI/REAL(NTHETA)*REAL(ITHETA-1)
                        
                        DL(1) =  0.d0
                        DL(2) = -WIRER*SIN(THETA)*DTHETA
                        DL(3) =  WIRER*COS(THETA)*DTHETA
                        
                        RPRIME(1) = POINT(1) - WIREZ
                        RPRIME(2) = POINT(2) - WIRER*COS(THETA)
                        RPRIME(3) = POINT(3) - WIRER*SIN(THETA)
                        
                        B_FIELD(IN, 1, :) = B_FIELD(IN, 1, :) + MU0/(4*PI)*CROSS(DL, RPRIME) / &
                        (RPRIME(1)*RPRIME(1) + RPRIME(2)*RPRIME(2) + RPRIME(3)*RPRIME(3))**1.5
                     
                     END DO
                  END DO
               END DO
            END DO
         END DO

      ELSE IF (N_SOLENOIDS > 0) THEN
         CALL ERROR_ABORT('Error! B field from solenoids is only implemented for unstructured grids!')
      END IF

   END SUBROUTINE COMPUTE_B_FIELD_FROM_SOLENOIDS



   SUBROUTINE APPLY_RF_E_FIELD(JP, E)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(3), INTENT(INOUT) :: E
      INTEGER, INTENT(IN) :: JP
      REAL(KIND=8) :: RF_XMIN, RF_XMAX, RF_YMAX
      REAL(KIND=8) :: RF_FREQ, NOVERL, COIL_CURRENT

      RF_XMIN = -0.08d0
      RF_XMAX = -0.055d0
      RF_YMAX = 0.015d0

      RF_FREQ = 13.56d6
      NOVERL = 250.d0
      COIL_CURRENT = 1.d0

      IF (particles(JP)%X > RF_XMIN .AND. particles(JP)%X < RF_XMAX .AND. particles(JP)%Y < RF_YMAX) THEN
         E(3) = E(3) - MU0*PI*RF_FREQ*NOVERL*COIL_CURRENT * particles(JP)%Y * SIN(2*PI*RF_FREQ*tID*DT)
      END IF

   END SUBROUTINE APPLY_RF_E_FIELD

END MODULE fields
