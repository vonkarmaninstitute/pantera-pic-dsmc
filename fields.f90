! Contains Electromagnetic fields related code

MODULE fields

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>

   USE global
   USE screen
   USE tools
   USE fully_implicit


   USE petsc
   USE petscksp

   IMPLICIT NONE

   CONTAINS



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

      INTEGER :: I, J
      INTEGER :: ICENTER, INORTH, IEAST, ISOUTH, IWEST
      INTEGER :: SIZE
      REAL(KIND=8) :: HX, HY
      REAL(KIND=8) :: AX, AY, BX, BY, CX, CY, H1X, H2X, H1Y, H2Y, R
      REAL(KIND=8) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4
      REAL(KIND=8) :: K11, K22, K33, K12, K23, K13, AREA, EDGELENGTH, FACEAREA
      INTEGER :: V1, V2, V3, V4
      INTEGER :: P, Q, VP, VQ
      REAL(KIND=8) :: KIJ, VOLUME
      INTEGER :: EDGE_PG
      LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_UNUSED

      SIZE = NNODES

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

      one = 1

      CALL MatCreate(PETSC_COMM_WORLD,Amat,ierr)
      CALL MatSetSizes( Amat,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
      CALL MatSetType( Amat, MATMPIAIJ, ierr)
      !CALL MatSetOption(Amat,MAT_SPD,PETSC_TRUE,ierr)
      CALL MatMPIAIJSetPreallocation(Amat,100,PETSC_NULL_INTEGER,100,PETSC_NULL_INTEGER, ierr)
      CALL MatSetFromOptions( Amat, ierr)
      CALL MatSetUp( Amat, ierr)
      CALL MatGetOwnershipRange( Amat, Istart, Iend, ierr)

      CALL MatCreateVecs( Amat, PETSC_NULL_VEC, xvec, ierr)
      CALL VecSetFromOptions( xvec, ierr)
      CALL VecDuplicate( xvec, bvec, ierr)


      ! At this point, populate the matrix
      IF (GRID_TYPE == UNSTRUCTURED) THEN

         IF (DIMS == 2) THEN

            ALLOCATE(IS_UNUSED(0:SIZE-1))
            IS_UNUSED = .TRUE.

            DO I = 1, NCELLS
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

            DO I = 1, NCELLS
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
                  IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                     CALL MatSetValues(Amat,one,V1-1,one,V1-1,K11,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V1-1,one,V2-1,K12,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V1-1,one,V3-1,K13,ADD_VALUES,ierr)
                  END IF
               END IF
               IF (V2-1 >= Istart .AND. V2-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                     CALL MatSetValues(Amat,one,V2-1,one,V1-1,K12,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V2-1,one,V3-1,K23,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V2-1,one,V2-1,K22,ADD_VALUES,ierr)
                  END IF
               END IF
               IF (V3-1 >= Istart .AND. V3-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                     CALL MatSetValues(Amat,one,V3-1,one,V1-1,K13,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V3-1,one,V2-1,K23,ADD_VALUES,ierr)
                     CALL MatSetValues(Amat,one,V3-1,one,V3-1,K33,ADD_VALUES,ierr)
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

         ELSE IF (DIMS == 3) THEN

            ALLOCATE(IS_UNUSED(0:SIZE-1))
            IS_UNUSED = .TRUE.

            DO I = 1, NCELLS
               V1 = U3D_GRID%CELL_NODES(I,1)
               V2 = U3D_GRID%CELL_NODES(I,2)
               V3 = U3D_GRID%CELL_NODES(I,3)
               V4 = U3D_GRID%CELL_NODES(I,4)
               IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.
               IS_UNUSED(V3-1) = .FALSE.; IS_UNUSED(V4-1) = .FALSE.

               DO J = 1, 4
                  EDGE_PG = U3D_GRID%CELL_FACES_PG(I, J)
                  IF (EDGE_PG .NE. -1) THEN
                     IF (GRID_BC(EDGE_PG)%FIELD_BC == DIRICHLET_BC) THEN
                        IF (J==1) THEN
                           DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V1-1) = .TRUE.
                           IS_DIRICHLET(V3-1) = .TRUE.
                           IS_DIRICHLET(V2-1) = .TRUE.
                        ELSE IF (J==2) THEN
                           DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V4-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V1-1) = .TRUE.
                           IS_DIRICHLET(V2-1) = .TRUE.
                           IS_DIRICHLET(V4-1) = .TRUE.
                        ELSE IF (J == 3) THEN
                           DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V4-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V2-1) = .TRUE.
                           IS_DIRICHLET(V3-1) = .TRUE.
                           IS_DIRICHLET(V4-1) = .TRUE.
                        ELSE
                           DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V4-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V1-1) = .TRUE.
                           IS_DIRICHLET(V4-1) = .TRUE.
                           IS_DIRICHLET(V3-1) = .TRUE.
                        END IF
                     END IF

                     IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN
                        IF (J==1) THEN
                           IS_NEUMANN(V1-1) = .TRUE.
                           IS_NEUMANN(V3-1) = .TRUE.
                           IS_NEUMANN(V2-1) = .TRUE.
                        ELSE IF (J==2) THEN
                           IS_NEUMANN(V1-1) = .TRUE.
                           IS_NEUMANN(V2-1) = .TRUE.
                           IS_NEUMANN(V4-1) = .TRUE.
                        ELSE IF (J==3) THEN
                           IS_NEUMANN(V2-1) = .TRUE.
                           IS_NEUMANN(V3-1) = .TRUE.
                           IS_NEUMANN(V4-1) = .TRUE.
                        ELSE
                           IS_NEUMANN(V1-1) = .TRUE.
                           IS_NEUMANN(V4-1) = .TRUE.
                           IS_NEUMANN(V3-1) = .TRUE.
                        END IF
                     END IF
                  END IF
               END DO
            END DO

            !IS_DIRICHLET(0) = .TRUE.    ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDDBDBDBDBDBDBDB
            !DIRICHLET(0) = 0.d0

            DO I = 1, NCELLS

               VOLUME = CELL_VOLUMES(I)

               DO P = 1, 4
                  VP = U3D_GRID%CELL_NODES(I,P) - 1
                  IF (VP >= Istart .AND. VP < Iend) THEN
                     IF (.NOT. IS_DIRICHLET(VP)) THEN
                        DO Q = 1, 4
                           VQ = U3D_GRID%CELL_NODES(I,Q) - 1
                           KIJ = VOLUME*(U3D_GRID%BASIS_COEFFS(I,P,1)*U3D_GRID%BASIS_COEFFS(I,Q,1) &
                                       + U3D_GRID%BASIS_COEFFS(I,P,2)*U3D_GRID%BASIS_COEFFS(I,Q,2) &
                                       + U3D_GRID%BASIS_COEFFS(I,P,3)*U3D_GRID%BASIS_COEFFS(I,Q,3))
                           CALL MatSetValues(Amat,one,VP,one,VQ,KIJ,ADD_VALUES,ierr)
                        END DO
                     END IF
                  END IF
               END DO

               DO J = 1, 4
                  EDGE_PG = U3D_GRID%CELL_FACES_PG(I, J)
                  IF (EDGE_PG == -1) CYCLE
                  IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN

                     FACEAREA = U3D_GRID%FACE_AREA(I, J)
                     IF (J==1) THEN
                        IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                           NEUMANN(V1-1) = NEUMANN(V1-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                        IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                           NEUMANN(V3-1) = NEUMANN(V3-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                        IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                           NEUMANN(V2-1) = NEUMANN(V2-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                     ELSE IF (J==2) THEN
                        IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                           NEUMANN(V1-1) = NEUMANN(V1-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                        IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                           NEUMANN(V2-1) = NEUMANN(V2-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                        IF (.NOT. IS_DIRICHLET(V4-1)) THEN
                           NEUMANN(V4-1) = NEUMANN(V4-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                     ELSE IF (J==3) THEN
                        IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                           NEUMANN(V2-1) = NEUMANN(V2-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                        IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                           NEUMANN(V3-1) = NEUMANN(V3-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                        IF (.NOT. IS_DIRICHLET(V4-1)) THEN
                           NEUMANN(V4-1) = NEUMANN(V4-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                     ELSE
                        IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                           NEUMANN(V1-1) = NEUMANN(V1-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                        IF (.NOT. IS_DIRICHLET(V4-1)) THEN
                           NEUMANN(V4-1) = NEUMANN(V4-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                        IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                           NEUMANN(V3-1) = NEUMANN(V3-1) + GRID_BC(EDGE_PG)%WALL_EFIELD * FACEAREA/3.
                        END IF
                     END IF

                  END IF
               END DO

            END DO

         END IF

         DO I = 1, NNODES
            IF (IS_UNUSED(I-1)) THEN
               !CALL MatSetValues(Amat,one,I-1,one,I-1,1.d0,ADD_VALUES,ierr)
               IS_DIRICHLET(I-1) = .TRUE.
               DIRICHLET(I-1) = 0.d0
            END IF
         END DO
         DEALLOCATE(IS_UNUSED)

      ELSE
         IF (DIMS == 1) THEN
            DO I = 0, NPX-1

               IF (I < Istart .OR. I >= Iend) CYCLE ! Each proc populates part of the matrix.
               ICENTER = I
               IEAST = I+1
               IWEST = I-1

               IF ((I == 0) .OR. (I == NPX-1)) THEN
                  ! Boundary point
                  DIRICHLET(ICENTER) = 0.d0
                  IS_DIRICHLET(ICENTER) = .TRUE.

               ELSE
                  ! Interior point.

                  IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
                     HX = (XMAX-XMIN)/DBLE(NX)
                     AX = 1./(HX*HX)
                     CX = AX
                     BX = -AX-CX
                  ELSE
                     H1X = XSIZE(I)
                     H2X = XSIZE(I+1)

                     AX = 2./(H1X*(H1X+H2X))
                     CX = 2./(H2X*(H1X+H2X))
                     BX = -AX-CX
                  END IF
                        
                  CALL MatSetValues(Amat,one,ICENTER,one,ICENTER,BX,INSERT_VALUES,ierr)
                  CALL MatSetValues(Amat,one,ICENTER,one,IEAST,CX,INSERT_VALUES,ierr)
                  CALL MatSetValues(Amat,one,ICENTER,one,IWEST,AX,INSERT_VALUES,ierr)

               END IF

            END DO
         ELSE IF (DIMS == 2) THEN
            DO I = 0, NPX-1
               DO J = 0, NPY-1


                  ICENTER = I+(NPX)*J
                  IF (ICENTER < Istart .OR. ICENTER >= Iend) CYCLE ! Each proc populates part of the matrix.
                  IEAST = I+(NPX)*J+1
                  IWEST = I+(NPX)*J-1
                  INORTH = I+(NPX)*(J+1)
                  ISOUTH = I+(NPX)*(J-1)


                  ! Various configurations go here.
                  
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! Plume 2d cartesian, full domain !
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  IF ((I == 0) .OR. (I == NPX-1) .OR. (J == 0) .OR. (J == NPY-1)) THEN
                     ! Boundary point
                     ! CALL MatSetValues(Amat,one,ICENTER,one,ICENTER,1.d0,INSERT_VALUES,ierr)
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
                     
                     CALL MatSetValues(Amat,one,ICENTER,one,ICENTER,BX+BY,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,INORTH,CY,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,ISOUTH,AY,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,IEAST,CX,INSERT_VALUES,ierr)
                     CALL MatSetValues(Amat,one,ICENTER,one,IWEST,AX,INSERT_VALUES,ierr)

                  END IF


               END DO
            END DO
         END IF
      END IF


      CALL MatAssemblyBegin(Amat,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Amat,MAT_FLUSH_ASSEMBLY,ierr)

      DO I = Istart, Iend-1
         IF (IS_DIRICHLET(I)) CALL MatSetValues(Amat,one,I,one,I,1.d0,INSERT_VALUES,ierr)
      END DO

      CALL MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)

      !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !CALL MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr)


   END SUBROUTINE ASSEMBLE_POISSON




   SUBROUTINE ASSEMBLE_AMPERE

      INTEGER :: I, J
      INTEGER :: MAXNNZ, SIZE
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, K11, K22, K33, K12, K23, K13, AREA, VOLUME
      REAL(KIND=8) :: K11TILDE, K22TILDE, K33TILDE, K12TILDE, K23TILDE, K13TILDE, KIJ
      INTEGER :: V1, V2, V3, V4, P, Q, VP, VQ
      INTEGER :: EDGE_PG
      LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_UNUSED

      SIZE = NNODES

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
      CALL MatSetType( Amat, MATMPIAIJ, ierr)
      !CALL MatSetOption(Amat,MAT_SPD,PETSC_TRUE,ierr)
      CALL MatMPIAIJSetPreallocation(Amat,30,PETSC_NULL_INTEGER,30,PETSC_NULL_INTEGER, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
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

         IF (DIMS == 2) THEN
            DO I = 1, NCELLS
               V1 = U2D_GRID%CELL_NODES(I,1)
               V2 = U2D_GRID%CELL_NODES(I,2)
               V3 = U2D_GRID%CELL_NODES(I,3)

               !IF ((V1-1 < Istart .OR. V1-1 >= Iend) .AND. &
               !    (V2-1 < Istart .OR. V2-1 >= Iend) .AND. &
               !    (V3-1 < Istart .OR. V3-1 >= Iend)) CYCLE

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

            DO I = 1, NCELLS
               V1 = U2D_GRID%CELL_NODES(I,1)
               V2 = U2D_GRID%CELL_NODES(I,2)
               V3 = U2D_GRID%CELL_NODES(I,3)
               
               !IF ((V1-1 < Istart .OR. V1-1 >= Iend) .AND. &
               !    (V2-1 < Istart .OR. V2-1 >= Iend) .AND. &
               !    (V3-1 < Istart .OR. V3-1 >= Iend)) CYCLE

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
               IF (V1-1 >= Istart .AND. V1-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V1-1)) THEN
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
                        val = PHI_FIELD(V1)*K11+PHI_FIELD(V2)*K12+PHI_FIELD(V3)*K13
                        CALL VecSetValues(bvec,one,V1-1,val,ADD_VALUES,ierr)
                     END IF
                  END IF
               END IF
               IF (V2-1 >= Istart .AND. V2-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                     IF (AXI) THEN
                        CALL MatSetValues(Amat,one,V2-1,one,V1-1,K12TILDE + MASS_MATRIX(I)*K12,ADD_VALUES,ierr)
                        CALL MatSetValues(Amat,one,V2-1,one,V3-1,K23TILDE + MASS_MATRIX(I)*K23,ADD_VALUES,ierr)
                        CALL MatSetValues(Amat,one,V2-1,one,V2-1,K22TILDE + MASS_MATRIX(I)*K22,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K12TILDE+PHI_FIELD(V2)*K22TILDE+PHI_FIELD(V3)*K23TILDE
                        CALL VecSetValues(bvec,one,V2-1,val,ADD_VALUES,ierr)
                     ELSE
                        CALL MatSetValues(Amat,one,V2-1,one,V1-1,(MASS_MATRIX(I)+1.)*K12,ADD_VALUES,ierr)
                        CALL MatSetValues(Amat,one,V2-1,one,V3-1,(MASS_MATRIX(I)+1.)*K23,ADD_VALUES,ierr)
                        CALL MatSetValues(Amat,one,V2-1,one,V2-1,(MASS_MATRIX(I)+1.)*K22,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K12+PHI_FIELD(V2)*K22+PHI_FIELD(V3)*K23
                        CALL VecSetValues(bvec,one,V2-1,val,ADD_VALUES,ierr)
                     END IF
                  END IF
               END IF
               IF (V3-1 >= Istart .AND. V3-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                     IF (AXI) THEN
                        CALL MatSetValues(Amat,one,V3-1,one,V1-1,K13TILDE + MASS_MATRIX(I)*K13,ADD_VALUES,ierr)
                        CALL MatSetValues(Amat,one,V3-1,one,V2-1,K23TILDE + MASS_MATRIX(I)*K23,ADD_VALUES,ierr)
                        CALL MatSetValues(Amat,one,V3-1,one,V3-1,K33TILDE + MASS_MATRIX(I)*K33,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K13TILDE+PHI_FIELD(V2)*K23TILDE+PHI_FIELD(V3)*K33TILDE
                        CALL VecSetValues(bvec,one,V3-1,val,ADD_VALUES,ierr)
                     ELSE

                        CALL MatSetValues(Amat,one,V3-1,one,V1-1,(MASS_MATRIX(I)+1.)*K13,ADD_VALUES,ierr)
                        CALL MatSetValues(Amat,one,V3-1,one,V2-1,(MASS_MATRIX(I)+1.)*K23,ADD_VALUES,ierr)
                        CALL MatSetValues(Amat,one,V3-1,one,V3-1,(MASS_MATRIX(I)+1.)*K33,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K13+PHI_FIELD(V2)*K23+PHI_FIELD(V3)*K33
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
         ELSE IF (DIMS == 3) THEN

            DO I = 1, NCELLS
               V1 = U3D_GRID%CELL_NODES(I,1)
               V2 = U3D_GRID%CELL_NODES(I,2)
               V3 = U3D_GRID%CELL_NODES(I,3)
               V4 = U3D_GRID%CELL_NODES(I,4)
               IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.
               IS_UNUSED(V3-1) = .FALSE.; IS_UNUSED(V4-1) = .FALSE.

               DO J = 1, 4
                  EDGE_PG = U3D_GRID%CELL_FACES_PG(I, J)
                  IF (EDGE_PG .NE. -1) THEN
                     IF (GRID_BC(EDGE_PG)%FIELD_BC == DIRICHLET_BC) THEN
                        IF (J==1) THEN
                           DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V1-1) = .TRUE.
                           IS_DIRICHLET(V3-1) = .TRUE.
                           IS_DIRICHLET(V2-1) = .TRUE.
                        ELSE IF (J==2) THEN
                           DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V4-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V1-1) = .TRUE.
                           IS_DIRICHLET(V2-1) = .TRUE.
                           IS_DIRICHLET(V4-1) = .TRUE.
                        ELSE IF (J == 3) THEN
                           DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V4-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V2-1) = .TRUE.
                           IS_DIRICHLET(V3-1) = .TRUE.
                           IS_DIRICHLET(V4-1) = .TRUE.
                        ELSE
                           DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V4-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           DIRICHLET(V3-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V1-1) = .TRUE.
                           IS_DIRICHLET(V4-1) = .TRUE.
                           IS_DIRICHLET(V3-1) = .TRUE.
                        END IF
                     END IF

                     IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN
                        IF (J==1) THEN
                           IS_NEUMANN(V1-1) = .TRUE.
                           IS_NEUMANN(V3-1) = .TRUE.
                           IS_NEUMANN(V2-1) = .TRUE.
                        ELSE IF (J==2) THEN
                           IS_NEUMANN(V1-1) = .TRUE.
                           IS_NEUMANN(V2-1) = .TRUE.
                           IS_NEUMANN(V4-1) = .TRUE.
                        ELSE IF (J==3) THEN
                           IS_NEUMANN(V2-1) = .TRUE.
                           IS_NEUMANN(V3-1) = .TRUE.
                           IS_NEUMANN(V4-1) = .TRUE.
                        ELSE
                           IS_NEUMANN(V1-1) = .TRUE.
                           IS_NEUMANN(V4-1) = .TRUE.
                           IS_NEUMANN(V3-1) = .TRUE.
                        END IF
                     END IF
                  END IF
               END DO
            END DO

            !IS_DIRICHLET(0) = .TRUE.    ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDDBDBDBDBDBDBDB
            !DIRICHLET(0) = 0.d0
            DO I = 1, NCELLS

               VOLUME = CELL_VOLUMES(I)

               DO P = 1, 4
                  VP = U3D_GRID%CELL_NODES(I,P) - 1
                  IF (VP >= Istart .AND. VP < Iend) THEN
                     IF (.NOT. IS_DIRICHLET(VP)) THEN
                        DO Q = 1, 4
                           VQ = U3D_GRID%CELL_NODES(I,Q) - 1
                           KIJ = VOLUME*(U3D_GRID%BASIS_COEFFS(I,P,1)*U3D_GRID%BASIS_COEFFS(I,Q,1) &
                                       + U3D_GRID%BASIS_COEFFS(I,P,2)*U3D_GRID%BASIS_COEFFS(I,Q,2) &
                                       + U3D_GRID%BASIS_COEFFS(I,P,3)*U3D_GRID%BASIS_COEFFS(I,Q,3))
                           CALL MatSetValues(Amat,one,VP,one,VQ,(MASS_MATRIX(I)+1.)*KIJ,ADD_VALUES,ierr)
                           val = PHI_FIELD(VQ+1)*KIJ
                           CALL VecSetValues(bvec,one,VP,val,ADD_VALUES,ierr)
                        END DO
                     END IF
                  END IF
               END DO
            
            END DO

         END IF

         DO I = Istart, Iend-1
            IF (IS_UNUSED(I)) THEN
               !CALL MatSetValues(Amat,one,I,one,I,1.d0,ADD_VALUES,ierr)
               IS_DIRICHLET(I) = .TRUE.
               DIRICHLET(I) = 0.d0
            ELSE IF (.NOT. IS_DIRICHLET(I) ) THEN
               val = 0.5/EPS0*J_FIELD(I)
               CALL VecSetValues(bvec,one,I,val,ADD_VALUES,ierr)
            END IF
         END DO
         DEALLOCATE(IS_UNUSED)

      ELSE
         CALL ERROR_ABORT('Semi-implicit with cartesian grid not implemented!')
      END IF


      CALL MatAssemblyBegin(Amat,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Amat,MAT_FLUSH_ASSEMBLY,ierr)

      DO I = Istart, Iend-1
         IF (IS_DIRICHLET(I)) THEN
            CALL MatSetValues(Amat,one,I,one,I,1.d0,INSERT_VALUES,ierr)

            val = DIRICHLET(I)
            CALL VecSetValues(bvec,one,I,val,ADD_VALUES,ierr)
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
      INTEGER :: V1, V2, V3, V4, SIZE, SIZEC, P, VP
      REAL(KIND=8) :: DPSI1DX, DPSI2DX, DPSI3DX, DPSI1DY, DPSI2DY, DPSI3DY


      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZE = NNODES
         SIZEC = NCELLS
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

            IF (DIMS == 2) THEN
               
               AREA = CELL_AREAS(IC)

               V1 = U2D_GRID%CELL_NODES(IC,1)
               V2 = U2D_GRID%CELL_NODES(IC,2)
               V3 = U2D_GRID%CELL_NODES(IC,3)            

               DPSI1DX = U2D_GRID%BASIS_COEFFS(IC,1,1)
               DPSI2DX = U2D_GRID%BASIS_COEFFS(IC,2,1)
               DPSI3DX = U2D_GRID%BASIS_COEFFS(IC,3,1)
               DPSI1DY = U2D_GRID%BASIS_COEFFS(IC,1,2)
               DPSI2DY = U2D_GRID%BASIS_COEFFS(IC,2,2)
               DPSI3DY = U2D_GRID%BASIS_COEFFS(IC,3,2)
               
               IF (AXI) THEN
                  J_FIELD(V1-1) = J_FIELD(V1-1) &
                  + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI1DX + particles(JP)%VY*DPSI1DY)/(ZMAX-ZMIN)*particles(JP)%DTRIM
                  J_FIELD(V2-1) = J_FIELD(V2-1) &
                  + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI2DX + particles(JP)%VY*DPSI2DY)/(ZMAX-ZMIN)*particles(JP)%DTRIM
                  J_FIELD(V3-1) = J_FIELD(V3-1) &
                  + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI3DX + particles(JP)%VY*DPSI3DY)/(ZMAX-ZMIN)*particles(JP)%DTRIM

                  MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*particles(JP)%DTRIM/EPS0/AREA/(ZMAX-ZMIN)*FNUM &
                                    * (QE*CHARGE)**2/SPECIES(particles(JP)%S_ID)%MOLECULAR_MASS
               ELSE
                  J_FIELD(V1-1) = J_FIELD(V1-1) &
                  + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI1DX + particles(JP)%VY*DPSI1DY)/(ZMAX-ZMIN)*particles(JP)%DTRIM
                  J_FIELD(V2-1) = J_FIELD(V2-1) &
                  + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI2DX + particles(JP)%VY*DPSI2DY)/(ZMAX-ZMIN)*particles(JP)%DTRIM
                  J_FIELD(V3-1) = J_FIELD(V3-1) &
                  + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI3DX + particles(JP)%VY*DPSI3DY)/(ZMAX-ZMIN)*particles(JP)%DTRIM

                  MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*particles(JP)%DTRIM/EPS0/AREA/(ZMAX-ZMIN)*FNUM &
                                    * (QE*CHARGE)**2/SPECIES(particles(JP)%S_ID)%MOLECULAR_MASS
               END IF

            ELSE IF (DIMS == 3) THEN
               
               DO P = 1, 4
                  VP = U3D_GRID%CELL_NODES(IC,P) - 1
                  J_FIELD(VP) = J_FIELD(VP) + FNUM*QE*CHARGE*(particles(JP)%VX*U3D_GRID%BASIS_COEFFS(IC,P,1) &
                                                            + particles(JP)%VY*U3D_GRID%BASIS_COEFFS(IC,P,2) &
                                                            + particles(JP)%VZ*U3D_GRID%BASIS_COEFFS(IC,P,3))*particles(JP)%DTRIM
               END DO
               MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*particles(JP)%DTRIM/EPS0/CELL_VOLUMES(IC)*FNUM &
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

      CALL MPI_BCAST(J_FIELD, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, MASS_MATRIX, SIZEC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ELSE
         CALL MPI_REDUCE(MASS_MATRIX,  MASS_MATRIX, SIZEC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF

      CALL MPI_BCAST(MASS_MATRIX, SIZEC, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   END SUBROUTINE DEPOSIT_CURRENT


   SUBROUTINE SOLVE_AMPERE

      IMPLICIT NONE

      CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
      CALL KSPSetOperators(ksp,Amat,Amat,ierr)

      !CALL KSPSetTolerances(ksp,1.d-11,1.d-11,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)

      CALL KSPSolve(ksp,bvec,xvec,ierr)
      
      CALL KSPGetConvergedReason(ksp,reason,ierr)
      IF (PROC_ID == 0) WRITE(*,*) 'KSPConvergedReason = ', reason
      
      CALL VecScatterCreateToAll(xvec,ctx,X_SEQ,ierr)
      CALL VecScatterBegin(ctx,xvec,X_SEQ,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterEnd(ctx,xvec,X_SEQ,INSERT_VALUES,SCATTER_FORWARD,ierr)

      CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD_TEMP,ierr)
      DEALLOCATE(PHIBAR_FIELD)
      ALLOCATE(PHIBAR_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(X_SEQ,PHI_FIELD_TEMP,ierr)

      PHI_FIELD = 2*PHIBAR_FIELD-PHI_FIELD

      CALL KSPDestroy(ksp,ierr)

   END SUBROUTINE SOLVE_AMPERE





   SUBROUTINE APPLY_POTENTIAL(JP, PHI)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(OUT) :: PHI
      INTEGER, INTENT(IN) :: JP
      REAL(KIND=8), DIMENSION(4) :: WEIGHTS
      INTEGER, DIMENSION(4) :: INDICES, INDI, INDJ
      INTEGER :: IC, V1, V2, V3
      REAL(KIND=8) :: AREA, X1, X2, X3, Y1, Y2, Y3, XP, YP, PSI1, PSI2, PSI3

      IF (DIMS == 2) THEN
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

            PHI = PSI1*PHI_FIELD(V1) + PSI2*PHI_FIELD(V2) + PSI3*PHI_FIELD(V3)
         ELSE
            CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
            PHI = WEIGHTS(1)*PHI_FIELD(INDICES(1)+1) + &
                  WEIGHTS(2)*PHI_FIELD(INDICES(2)+1) + &
                  WEIGHTS(3)*PHI_FIELD(INDICES(3)+1) + &
                  WEIGHTS(4)*PHI_FIELD(INDICES(4)+1)
         END IF
      ELSE
         CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
         PHI = WEIGHTS(1)*PHI_FIELD(INDICES(1)+1) + &
               WEIGHTS(2)*PHI_FIELD(INDICES(2)+1)
      END IF
      
   END SUBROUTINE APPLY_POTENTIAL




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
                        
                        B_FIELD(IN, 1, :) = B_FIELD(IN, 1, :) + SOLENOIDS(ICOIL)%WIRE_CURRENT*MU0/(4*PI)*CROSS(DL, RPRIME) / &
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
