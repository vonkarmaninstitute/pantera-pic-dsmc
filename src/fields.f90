! Copyright (C) 2024 von Karman Institute for Fluid Dynamics (VKI)
!
! This file is part of PANTERA PIC-DSMC, a software for the simulation
! of rarefied gases and plasmas using particles.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.PANTERA PIC-DSMC

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

      CALL PetscPopSignalHandler(ierr)

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
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3
      REAL(KIND=8) :: K11, K22, K33, K12, K23, K13, AREA, EDGELENGTH, FACEAREA
      INTEGER :: V1, V2, V3, V4
      INTEGER :: P, Q, VP, VQ
      REAL(KIND=8) :: KIJ, VOLUME, LENGTH
      INTEGER :: EDGE_PG
      LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_UNUSED
      REAL(KIND=8) :: EPS_REL

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

      CALL MatDestroy(Amat,ierr)
      CALL MatCreate(PETSC_COMM_WORLD,Amat,ierr)
      CALL MatSetSizes( Amat,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
      CALL MatSetType( Amat, MATMPIAIJ, ierr)
      !CALL MatSetOption(Amat,MAT_SPD,PETSC_TRUE,ierr)
      CALL MatMPIAIJSetPreallocation(Amat,100,PETSC_NULL_INTEGER_ARRAY,100,PETSC_NULL_INTEGER_ARRAY, ierr)
      CALL MatSetFromOptions( Amat, ierr)
      CALL MatSetUp( Amat, ierr)
      CALL MatGetOwnershipRange( Amat, Istart, Iend, ierr)

      CALL VecDestroy(xvec,ierr)
      CALL MatCreateVecs( Amat, PETSC_NULL_VEC, xvec, ierr)
      CALL VecSetFromOptions( xvec, ierr)
      CALL VecDestroy(bvec,ierr)
      CALL VecDuplicate( xvec, bvec, ierr)


      ! At this point, populate the matrix
      IF (GRID_TYPE == UNSTRUCTURED) THEN
         IF (DIMS == 1) THEN

            ALLOCATE(IS_UNUSED(0:SIZE-1))
            IS_UNUSED = .TRUE.

            DO I = 1, NCELLS
               V1 = U1D_GRID%CELL_NODES(1,I)
               V2 = U1D_GRID%CELL_NODES(2,I)
               IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.
               DO J = 1, 2
                  EDGE_PG = U1D_GRID%CELL_EDGES_PG(J, I)
                  IF (EDGE_PG .NE. -1) THEN
                     IF (GRID_BC(EDGE_PG)%FIELD_BC == DIRICHLET_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == RF_VOLTAGE_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == DECOUPLED_RF_VOLTAGE_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == CONDUCTIVE_BC) THEN
                        IF (J==1) THEN
                           DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V1-1) = .TRUE.
                        ELSE IF (J==2) THEN
                           DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V2-1) = .TRUE.
                        END IF
                     END IF

                     IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN
                        IF (J==1) THEN
                           IS_NEUMANN(V1-1) = .TRUE.
                        ELSE IF (J==2) THEN
                           IS_NEUMANN(V2-1) = .TRUE.
                        END IF
                     END IF
                  END IF
               END DO
            END DO


            IF (PIC_TYPE == EXPLICITLIMITED) THEN
               CALL COMPUTE_DENSITY_TEMPERATURE(particles)
               IF (.NOT. ALLOCATED(DXLDRATIO)) ALLOCATE(DXLDRATIO(NCELLS))
               DO I = 1, NCELLS
                  LENGTH = U1D_GRID%SEGMENT_LENGTHS(I)
                  IF (CELL_TE(I) > 0) THEN
                     DXLDRATIO(I) = LENGTH**2*(CELL_NE(I)*QE**2/(EPS0*KB*CELL_TE(I)))
                     !WRITE(*,*) 'Cell ', I, ' T= ', CELL_TE(I), ' n= ', CELL_NE(I), ' dx^2/lD^2= ', RATIO
                  ELSE
                     DXLDRATIO(I) = 0
                  END IF
               END DO
            END IF

            !IS_DIRICHLET(0) = .TRUE.    ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDDBDBDBDBDBDBDB
            !DIRICHLET(0) = 0.d0

            DO I = 1, NCELLS
               IF (U1D_GRID%CELL_PG(I) == -1) THEN
                  EPS_REL = 1.d0
               ELSE
                  EPS_REL = GRID_BC(U1D_GRID%CELL_PG(I))%EPS_REL
               END IF
               LENGTH = U1D_GRID%SEGMENT_LENGTHS(I)
               V1 = U1D_GRID%CELL_NODES(1,I)
               V2 = U1D_GRID%CELL_NODES(2,I)


               K11 = 1.0/LENGTH*EPS_REL
               K22 = 1.0/LENGTH*EPS_REL
               K12 =-1.0/LENGTH*EPS_REL


               IF (PIC_TYPE == EXPLICITLIMITED) THEN
                  K11 = K11 * (1. + DXLDRATIO(I))
                  K22 = K22 * (1. + DXLDRATIO(I))
                  K12 = K12 * (1. + DXLDRATIO(I))
               END IF

               ! We need to ADD to a sparse matrix entry.
               IF (V1-1 >= Istart .AND. V1-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                     CALL MatSetValue(Amat,V1-1,V1-1,K11,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V1-1,V2-1,K12,ADD_VALUES,ierr)
                  END IF
               END IF
               IF (V2-1 >= Istart .AND. V2-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                     CALL MatSetValue(Amat,V2-1,V1-1,K12,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V2-1,V2-1,K22,ADD_VALUES,ierr)
                  END IF
               END IF

               DO J = 1, 2
                  EDGE_PG = U1D_GRID%CELL_EDGES_PG(J, I)
                  IF (EDGE_PG == -1) CYCLE
                  IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN

                     IF (J==1) THEN
                        IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                           NEUMANN(V1-1) = NEUMANN(V1-1) + GRID_BC(EDGE_PG)%WALL_EFIELD
                        END IF
                     ELSE IF (J==2) THEN
                        IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                           NEUMANN(V2-1) = NEUMANN(V2-1) + GRID_BC(EDGE_PG)%WALL_EFIELD
                        END IF
                     END IF

                  END IF
               END DO

            END DO

         ELSE IF (DIMS == 2) THEN

            ALLOCATE(IS_UNUSED(0:SIZE-1))
            IS_UNUSED = .TRUE.

            DO I = 1, NCELLS
               V1 = U2D_GRID%CELL_NODES(1,I)
               V2 = U2D_GRID%CELL_NODES(2,I)
               V3 = U2D_GRID%CELL_NODES(3,I)
               IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.; IS_UNUSED(V3-1) = .FALSE.
               DO J = 1, 3
                  EDGE_PG = U2D_GRID%CELL_EDGES_PG(J, I)
                  IF (EDGE_PG .NE. -1) THEN
                     IF (GRID_BC(EDGE_PG)%FIELD_BC == DIRICHLET_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == RF_VOLTAGE_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == DECOUPLED_RF_VOLTAGE_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == CONDUCTIVE_BC) THEN
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


            IF (PIC_TYPE == EXPLICITLIMITED) THEN
               CALL COMPUTE_DENSITY_TEMPERATURE(particles)
               IF (.NOT. ALLOCATED(DXLDRATIO)) ALLOCATE(DXLDRATIO(NCELLS))
               DO I = 1, NCELLS
                  LENGTH = U1D_GRID%SEGMENT_LENGTHS(I)
                  IF (CELL_TE(I) > 0) THEN
                     DXLDRATIO(I) = LENGTH**2*(CELL_NE(I)*QE**2/(EPS0*KB*CELL_TE(I)))
                     !WRITE(*,*) 'Cell ', I, ' T= ', CELL_TE(I), ' n= ', CELL_NE(I), ' dx^2/lD^2= ', RATIO
                  ELSE
                     DXLDRATIO(I) = 0
                  END IF
               END DO
            END IF

            !IS_DIRICHLET(0) = .TRUE.    ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDDBDBDBDBDBDBDB
            !DIRICHLET(0) = 0.d0

            DO I = 1, NCELLS
               IF (U2D_GRID%CELL_PG(I) == -1) THEN
                  EPS_REL = 1.d0
               ELSE
                  EPS_REL = GRID_BC(U2D_GRID%CELL_PG(I))%EPS_REL
               END IF
               AREA = U2D_GRID%CELL_AREAS(I)
               V1 = U2D_GRID%CELL_NODES(1,I)
               V2 = U2D_GRID%CELL_NODES(2,I)
               V3 = U2D_GRID%CELL_NODES(3,I)            
               X1 = U2D_GRID%NODE_COORDS(1, V1)
               X2 = U2D_GRID%NODE_COORDS(1, V2)
               X3 = U2D_GRID%NODE_COORDS(1, V3)
               Y1 = U2D_GRID%NODE_COORDS(2, V1)
               Y2 = U2D_GRID%NODE_COORDS(2, V2)
               Y3 = U2D_GRID%NODE_COORDS(2, V3)
               K11 = 0.25*((Y2-Y3)**2 + (X2-X3)**2)/AREA*EPS_REL
               K22 = 0.25*((Y1-Y3)**2 + (X1-X3)**2)/AREA*EPS_REL
               K33 = 0.25*((Y2-Y1)**2 + (X2-X1)**2)/AREA*EPS_REL
               K12 =-0.25*((Y2-Y3)*(Y1-Y3) + (X2-X3)*(X1-X3))/AREA*EPS_REL
               K23 = 0.25*((Y1-Y3)*(Y2-Y1) + (X1-X3)*(X2-X1))/AREA*EPS_REL
               K13 =-0.25*((Y2-Y3)*(Y2-Y1) + (X2-X3)*(X2-X1))/AREA*EPS_REL
               IF (AXI) THEN
                  K11 = K11*(Y1+Y2+Y3)/3.
                  K22 = K22*(Y1+Y2+Y3)/3.
                  K33 = K33*(Y1+Y2+Y3)/3.
                  K12 = K12*(Y1+Y2+Y3)/3.
                  K23 = K23*(Y1+Y2+Y3)/3.
                  K13 = K13*(Y1+Y2+Y3)/3.
               END IF

               IF (PIC_TYPE == EXPLICITLIMITED) THEN
                  K11 = K11 * (1. + DXLDRATIO(I))
                  K22 = K22 * (1. + DXLDRATIO(I))
                  K33 = K33 * (1. + DXLDRATIO(I))
                  K12 = K12 * (1. + DXLDRATIO(I))
                  K23 = K23 * (1. + DXLDRATIO(I))
                  K13 = K13 * (1. + DXLDRATIO(I))
               END IF

               ! We need to ADD to a sparse matrix entry.
               IF (V1-1 >= Istart .AND. V1-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                     CALL MatSetValue(Amat,V1-1,V1-1,K11,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V1-1,V2-1,K12,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V1-1,V3-1,K13,ADD_VALUES,ierr)
                  END IF
               END IF
               IF (V2-1 >= Istart .AND. V2-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                     CALL MatSetValue(Amat,V2-1,V1-1,K12,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V2-1,V3-1,K23,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V2-1,V2-1,K22,ADD_VALUES,ierr)
                  END IF
               END IF
               IF (V3-1 >= Istart .AND. V3-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                     CALL MatSetValue(Amat,V3-1,V1-1,K13,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V3-1,V2-1,K23,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V3-1,V3-1,K33,ADD_VALUES,ierr)
                  END IF
               END IF

               DO J = 1, 3
                  EDGE_PG = U2D_GRID%CELL_EDGES_PG(J, I)
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
               V1 = U3D_GRID%CELL_NODES(1,I)
               V2 = U3D_GRID%CELL_NODES(2,I)
               V3 = U3D_GRID%CELL_NODES(3,I)
               V4 = U3D_GRID%CELL_NODES(4,I)
               IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.
               IS_UNUSED(V3-1) = .FALSE.; IS_UNUSED(V4-1) = .FALSE.

               DO J = 1, 4
                  EDGE_PG = U3D_GRID%CELL_FACES_PG(J, I)
                  IF (EDGE_PG .NE. -1) THEN
                     IF (GRID_BC(EDGE_PG)%FIELD_BC == DIRICHLET_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == RF_VOLTAGE_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == DECOUPLED_RF_VOLTAGE_BC &
                         .OR. GRID_BC(EDGE_PG)%FIELD_BC == CONDUCTIVE_BC) THEN
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
            IF (PIC_TYPE == EXPLICITLIMITED) THEN
               CALL COMPUTE_DENSITY_TEMPERATURE(particles)
               IF (.NOT. ALLOCATED(DXLDRATIO)) ALLOCATE(DXLDRATIO(NCELLS))
               DO I = 1, NCELLS
                  IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
                     VOLUME = U1D_GRID%CELL_VOLUMES(I)
                  ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
                     VOLUME = U2D_GRID%CELL_VOLUMES(I)
                  ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
                     VOLUME = U3D_GRID%CELL_VOLUMES(I)
                  END IF

                  IF (CELL_TE(I) > 0) THEN
                     DXLDRATIO(I) = VOLUME**(2./3.)*(CELL_NE(I)*QE**2/(EPS0*KB*CELL_TE(I)))
                     !WRITE(*,*) 'Cell ', I, ' T= ', CELL_TE(I), ' n= ', CELL_NE(I), ' dx^2/lD^2= ', RATIO
                  ELSE
                     DXLDRATIO(I) = 0
                  END IF
               END DO
            END IF

            !WRITE(*,*) 'PROC= ', PROC_ID, 'Max DXLDRATIO= ', MAXVAL(DXLDRATIO), 'Min DXLDRATIO= ', MINVAL(DXLDRATIO)

            DO I = 1, NCELLS
               IF (U3D_GRID%CELL_PG(I) == -1) THEN
                  EPS_REL = 1.d0
               ELSE
                  EPS_REL = GRID_BC(U3D_GRID%CELL_PG(I))%EPS_REL
               END IF

               IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
                  VOLUME = U1D_GRID%CELL_VOLUMES(I)
               ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
                  VOLUME = U2D_GRID%CELL_VOLUMES(I)
               ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
                  VOLUME = U3D_GRID%CELL_VOLUMES(I)
               END IF

               DO P = 1, 4

                  VP = U3D_GRID%CELL_NODES(P,I) - 1
                  IF (VP >= Istart .AND. VP < Iend) THEN
                     IF (.NOT. IS_DIRICHLET(VP)) THEN
                        DO Q = 1, 4
                           VQ = U3D_GRID%CELL_NODES(Q,I) - 1
                           KIJ = VOLUME*(U3D_GRID%BASIS_COEFFS(1,P,I)*U3D_GRID%BASIS_COEFFS(1,Q,I) &
                                       + U3D_GRID%BASIS_COEFFS(2,P,I)*U3D_GRID%BASIS_COEFFS(2,Q,I) &
                                       + U3D_GRID%BASIS_COEFFS(3,P,I)*U3D_GRID%BASIS_COEFFS(3,Q,I)) * EPS_REL
                           IF (PIC_TYPE == EXPLICITLIMITED) KIJ = KIJ * (1. + DXLDRATIO(I))
                           CALL MatSetValue(Amat,VP,VQ,KIJ,ADD_VALUES,ierr)
                        END DO
                     END IF
                  END IF
               END DO

               DO J = 1, 4
                  EDGE_PG = U3D_GRID%CELL_FACES_PG(J, I)
                  IF (EDGE_PG == -1) CYCLE
                  IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN

                     FACEAREA = U3D_GRID%FACE_AREA(J, I)
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
               !CALL MatSetValue(Amat,I-1,I-1,1.d0,ADD_VALUES,ierr)
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
                        
                  CALL MatSetValue(Amat,ICENTER,ICENTER,BX,INSERT_VALUES,ierr)
                  CALL MatSetValue(Amat,ICENTER,IEAST,CX,INSERT_VALUES,ierr)
                  CALL MatSetValue(Amat,ICENTER,IWEST,AX,INSERT_VALUES,ierr)

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
                     ! CALL MatSetValue(Amat,ICENTER,ICENTER,1.d0,INSERT_VALUES,ierr)
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
                     
                     CALL MatSetValue(Amat,ICENTER,ICENTER,BX+BY,INSERT_VALUES,ierr)
                     CALL MatSetValue(Amat,ICENTER,INORTH,CY,INSERT_VALUES,ierr)
                     CALL MatSetValue(Amat,ICENTER,ISOUTH,AY,INSERT_VALUES,ierr)
                     CALL MatSetValue(Amat,ICENTER,IEAST,CX,INSERT_VALUES,ierr)
                     CALL MatSetValue(Amat,ICENTER,IWEST,AX,INSERT_VALUES,ierr)

                  END IF


               END DO
            END DO
         END IF
      END IF


      CALL MatAssemblyBegin(Amat,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Amat,MAT_FLUSH_ASSEMBLY,ierr)

      DO I = Istart, Iend-1
         IF (IS_DIRICHLET(I)) CALL MatSetValue(Amat,I,I,1.d0,INSERT_VALUES,ierr)
      END DO

      CALL MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)

      !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !CALL MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr)


   END SUBROUTINE ASSEMBLE_POISSON




   SUBROUTINE ASSEMBLE_AMPERE

      INTEGER :: I, J
      INTEGER :: SIZE
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, K11, K22, K33, K12, K23, K13, AREA, VOLUME, LENGTH
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
      CALL MatMPIAIJSetPreallocation(Amat,30,PETSC_NULL_INTEGER_ARRAY,30,PETSC_NULL_INTEGER_ARRAY, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
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
         IF (DIMS == 1) THEN
            DO I = 1, NCELLS
               V1 = U2D_GRID%CELL_NODES(1,I)
               V2 = U2D_GRID%CELL_NODES(2,I)

               IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.
               DO J = 1, 2
                  EDGE_PG = U1D_GRID%CELL_EDGES_PG(J, I)
                  IF (EDGE_PG .NE. -1) THEN
                     IF (GRID_BC(EDGE_PG)%FIELD_BC == DIRICHLET_BC) THEN
                        IF (J==1) THEN
                           DIRICHLET(V1-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V1-1) = .TRUE.
                        ELSE IF (J==2) THEN
                           DIRICHLET(V2-1) = GRID_BC(EDGE_PG)%WALL_POTENTIAL
                           IS_DIRICHLET(V2-1) = .TRUE.
                        END IF
                     END IF

                     IF (GRID_BC(EDGE_PG)%FIELD_BC == NEUMANN_BC) THEN
                        IF (J==1) THEN
                           IS_NEUMANN(V1-1) = .TRUE.
                        ELSE IF (J==2) THEN
                           IS_NEUMANN(V2-1) = .TRUE.
                        END IF
                     END IF
                  END IF
               END DO
            END DO

            !IS_DIRICHLET(0) = .TRUE. ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDDBDBDBDBDBDBDB
            !DIRICHLET(0) = 0.d0

            DO I = 1, NCELLS
               LENGTH = U1D_GRID%SEGMENT_LENGTHS(I)
               V1 = U1D_GRID%CELL_NODES(1,I)
               V2 = U1D_GRID%CELL_NODES(2,I)


               K11 = 1.0/LENGTH
               K22 = 1.0/LENGTH
               K12 =-1.0/LENGTH


               ! We need to ADD to a sparse matrix entry.
               IF (V1-1 >= Istart .AND. V1-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V1-1)) THEN
                     CALL MatSetValue(Amat,V1-1,V1-1,(MASS_MATRIX(I)+1.)*K11,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V1-1,V2-1,(MASS_MATRIX(I)+1.)*K12,ADD_VALUES,ierr)
                     val = PHI_FIELD(V1)*K11+PHI_FIELD(V2)*K12
                     CALL VecSetValue(bvec,V1-1,val,ADD_VALUES,ierr)
                  END IF
               END IF
               IF (V2-1 >= Istart .AND. V2-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                     CALL MatSetValue(Amat,V2-1,V1-1,(MASS_MATRIX(I)+1.)*K12,ADD_VALUES,ierr)
                     CALL MatSetValue(Amat,V2-1,V2-1,(MASS_MATRIX(I)+1.)*K22,ADD_VALUES,ierr)
                     val = PHI_FIELD(V1)*K12+PHI_FIELD(V2)*K22
                     CALL VecSetValue(bvec,V2-1,val,ADD_VALUES,ierr)
                  END IF
               END IF

            END DO
         ELSE IF (DIMS == 2) THEN
            DO I = 1, NCELLS
               V1 = U2D_GRID%CELL_NODES(1,I)
               V2 = U2D_GRID%CELL_NODES(2,I)
               V3 = U2D_GRID%CELL_NODES(3,I)

               !IF ((V1-1 < Istart .OR. V1-1 >= Iend) .AND. &
               !    (V2-1 < Istart .OR. V2-1 >= Iend) .AND. &
               !    (V3-1 < Istart .OR. V3-1 >= Iend)) CYCLE

               IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.; IS_UNUSED(V3-1) = .FALSE.
               DO J = 1, 3
                  EDGE_PG = U2D_GRID%CELL_EDGES_PG(J, I)
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
               V1 = U2D_GRID%CELL_NODES(1,I)
               V2 = U2D_GRID%CELL_NODES(2,I)
               V3 = U2D_GRID%CELL_NODES(3,I)
               
               !IF ((V1-1 < Istart .OR. V1-1 >= Iend) .AND. &
               !    (V2-1 < Istart .OR. V2-1 >= Iend) .AND. &
               !    (V3-1 < Istart .OR. V3-1 >= Iend)) CYCLE

               AREA = U2D_GRID%CELL_AREAS(I)
               X1 = U2D_GRID%NODE_COORDS(1, V1)
               X2 = U2D_GRID%NODE_COORDS(1, V2)
               X3 = U2D_GRID%NODE_COORDS(1, V3)
               Y1 = U2D_GRID%NODE_COORDS(2, V1)
               Y2 = U2D_GRID%NODE_COORDS(2, V2)
               Y3 = U2D_GRID%NODE_COORDS(2, V3)
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
                        CALL MatSetValue(Amat,V1-1,V1-1,K11TILDE + MASS_MATRIX(I)*K11,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V1-1,V2-1,K12TILDE + MASS_MATRIX(I)*K12,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V1-1,V3-1,K13TILDE + MASS_MATRIX(I)*K13,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1-1)*K11TILDE+PHI_FIELD(V2-1)*K12TILDE+PHI_FIELD(V3-1)*K13TILDE
                        CALL VecSetValue(bvec,V1-1,val,ADD_VALUES,ierr)
                     ELSE
                        CALL MatSetValue(Amat,V1-1,V1-1,(MASS_MATRIX(I)+1.)*K11,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V1-1,V2-1,(MASS_MATRIX(I)+1.)*K12,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V1-1,V3-1,(MASS_MATRIX(I)+1.)*K13,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K11+PHI_FIELD(V2)*K12+PHI_FIELD(V3)*K13
                        CALL VecSetValue(bvec,V1-1,val,ADD_VALUES,ierr)
                     END IF
                  END IF
               END IF
               IF (V2-1 >= Istart .AND. V2-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V2-1)) THEN
                     IF (AXI) THEN
                        CALL MatSetValue(Amat,V2-1,V1-1,K12TILDE + MASS_MATRIX(I)*K12,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V2-1,V3-1,K23TILDE + MASS_MATRIX(I)*K23,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V2-1,V2-1,K22TILDE + MASS_MATRIX(I)*K22,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K12TILDE+PHI_FIELD(V2)*K22TILDE+PHI_FIELD(V3)*K23TILDE
                        CALL VecSetValue(bvec,V2-1,val,ADD_VALUES,ierr)
                     ELSE
                        CALL MatSetValue(Amat,V2-1,V1-1,(MASS_MATRIX(I)+1.)*K12,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V2-1,V3-1,(MASS_MATRIX(I)+1.)*K23,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V2-1,V2-1,(MASS_MATRIX(I)+1.)*K22,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K12+PHI_FIELD(V2)*K22+PHI_FIELD(V3)*K23
                        CALL VecSetValue(bvec,V2-1,val,ADD_VALUES,ierr)
                     END IF
                  END IF
               END IF
               IF (V3-1 >= Istart .AND. V3-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(V3-1)) THEN
                     IF (AXI) THEN
                        CALL MatSetValue(Amat,V3-1,V1-1,K13TILDE + MASS_MATRIX(I)*K13,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V3-1,V2-1,K23TILDE + MASS_MATRIX(I)*K23,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V3-1,V3-1,K33TILDE + MASS_MATRIX(I)*K33,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K13TILDE+PHI_FIELD(V2)*K23TILDE+PHI_FIELD(V3)*K33TILDE
                        CALL VecSetValue(bvec,V3-1,val,ADD_VALUES,ierr)
                     ELSE

                        CALL MatSetValue(Amat,V3-1,V1-1,(MASS_MATRIX(I)+1.)*K13,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V3-1,V2-1,(MASS_MATRIX(I)+1.)*K23,ADD_VALUES,ierr)
                        CALL MatSetValue(Amat,V3-1,V3-1,(MASS_MATRIX(I)+1.)*K33,ADD_VALUES,ierr)
                        val = PHI_FIELD(V1)*K13+PHI_FIELD(V2)*K23+PHI_FIELD(V3)*K33
                        CALL VecSetValue(bvec,V3-1,val,ADD_VALUES,ierr)
                     END IF
                  END IF
               END IF
            ! Neumann part has to be included only if the derivative changes in time.
            !    DO J = 1, 3
            !       EDGE_PG = U2D_GRID%CELL_EDGES_PG(J, I)
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
               V1 = U3D_GRID%CELL_NODES(1,I)
               V2 = U3D_GRID%CELL_NODES(2,I)
               V3 = U3D_GRID%CELL_NODES(3,I)
               V4 = U3D_GRID%CELL_NODES(4,I)
               IS_UNUSED(V1-1) = .FALSE.; IS_UNUSED(V2-1) = .FALSE.
               IS_UNUSED(V3-1) = .FALSE.; IS_UNUSED(V4-1) = .FALSE.

               DO J = 1, 4
                  EDGE_PG = U3D_GRID%CELL_FACES_PG(J, I)
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

               IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
                  VOLUME = U2D_GRID%CELL_VOLUMES(I)
               ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
                  VOLUME = U3D_GRID%CELL_VOLUMES(I)
               END IF

               DO P = 1, 4
                  VP = U3D_GRID%CELL_NODES(P,I) - 1
                  IF (VP >= Istart .AND. VP < Iend) THEN
                     IF (.NOT. IS_DIRICHLET(VP)) THEN
                        DO Q = 1, 4
                           VQ = U3D_GRID%CELL_NODES(Q,I) - 1
                           KIJ = VOLUME*(U3D_GRID%BASIS_COEFFS(1,P,I)*U3D_GRID%BASIS_COEFFS(1,Q,I) &
                                       + U3D_GRID%BASIS_COEFFS(2,P,I)*U3D_GRID%BASIS_COEFFS(2,Q,I) &
                                       + U3D_GRID%BASIS_COEFFS(3,P,I)*U3D_GRID%BASIS_COEFFS(3,Q,I))
                           CALL MatSetValue(Amat,VP,VQ,(MASS_MATRIX(I)+1.)*KIJ,ADD_VALUES,ierr)
                           val = PHI_FIELD(VQ+1)*KIJ
                           CALL VecSetValue(bvec,VP,val,ADD_VALUES,ierr)
                        END DO
                     END IF
                  END IF
               END DO
            
            END DO

         END IF

         DO I = Istart, Iend-1
            IF (IS_UNUSED(I)) THEN
               !CALL MatSetValue(Amat,I,I,1.d0,ADD_VALUES,ierr)
               IS_DIRICHLET(I) = .TRUE.
               DIRICHLET(I) = 0.d0
            ELSE IF (.NOT. IS_DIRICHLET(I) ) THEN
               val = 0.5/EPS0*J_FIELD(I)
               CALL VecSetValue(bvec,I,val,ADD_VALUES,ierr)
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
            CALL MatSetValue(Amat,I,I,1.d0,INSERT_VALUES,ierr)

            val = DIRICHLET(I)
            CALL VecSetValue(bvec,I,val,ADD_VALUES,ierr)
         ELSE IF (IS_NEUMANN(I)) THEN
            val = NEUMANN(I)
            CALL VecSetValue(bvec,I,val,ADD_VALUES,ierr)
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
      INTEGER :: V1, V2, V3, SIZE, SIZEC, P, VP
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

            IF (DIMS == 1) THEN

               V1 = U1D_GRID%CELL_NODES(1,IC)
               V2 = U1D_GRID%CELL_NODES(2,IC)

               DPSI1DX = U1D_GRID%BASIS_COEFFS(1,1,IC)
               DPSI2DX = U1D_GRID%BASIS_COEFFS(1,2,IC)
               
               J_FIELD(V1-1) = J_FIELD(V1-1) &
               + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI1DX)/(YMAX-YMIN)/(ZMAX-ZMIN)*particles(JP)%DTRIM
               J_FIELD(V2-1) = J_FIELD(V2-1) &
               + FNUM*QE*CHARGE*(particles(JP)%VX*DPSI2DX)/(YMAX-YMIN)/(ZMAX-ZMIN)*particles(JP)%DTRIM
              
               MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*particles(JP)%DTRIM/EPS0/U1D_GRID%CELL_VOLUMES(IC)*FNUM &
                                 * (QE*CHARGE)**2/SPECIES(particles(JP)%S_ID)%MOLECULAR_MASS

            ELSE IF (DIMS == 2) THEN
               
               AREA = U2D_GRID%CELL_AREAS(IC)

               V1 = U2D_GRID%CELL_NODES(1,IC)
               V2 = U2D_GRID%CELL_NODES(2,IC)
               V3 = U2D_GRID%CELL_NODES(3,IC)            

               DPSI1DX = U2D_GRID%BASIS_COEFFS(1,1,IC)
               DPSI2DX = U2D_GRID%BASIS_COEFFS(1,2,IC)
               DPSI3DX = U2D_GRID%BASIS_COEFFS(1,3,IC)
               DPSI1DY = U2D_GRID%BASIS_COEFFS(2,1,IC)
               DPSI2DY = U2D_GRID%BASIS_COEFFS(2,2,IC)
               DPSI3DY = U2D_GRID%BASIS_COEFFS(2,3,IC)
               
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
                  VP = U3D_GRID%CELL_NODES(P,IC) - 1
                  J_FIELD(VP) = J_FIELD(VP) + FNUM*QE*CHARGE*(particles(JP)%VX*U3D_GRID%BASIS_COEFFS(1,P,IC) &
                                                            + particles(JP)%VY*U3D_GRID%BASIS_COEFFS(2,P,IC) &
                                                            + particles(JP)%VZ*U3D_GRID%BASIS_COEFFS(3,P,IC))*particles(JP)%DTRIM
               END DO
               MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*particles(JP)%DTRIM/EPS0/U3D_GRID%CELL_VOLUMES(IC)*FNUM &
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
      INTEGER :: IC, VI, I
      REAL(KIND=8) :: XP, YP, ZP, PSII

      
      IF (GRID_TYPE == UNSTRUCTURED) THEN
         IF (DIMS == 1) THEN

            IC = particles(JP)%IC

            XP = particles(JP)%X

            PHI = 0.d0
            DO I = 1, 2
               VI = U1D_GRID%CELL_NODES(I,IC)
               PSII = XP*U1D_GRID%BASIS_COEFFS(1,I,IC) + U1D_GRID%BASIS_COEFFS(2,I,IC)
               PHI = PHI + PSII*PHI_FIELD(VI)
            END DO

         ELSE IF (DIMS == 2) THEN

            IC = particles(JP)%IC

            XP = particles(JP)%X
            YP = particles(JP)%Y

            PHI = 0.d0
            DO I = 1, 3
               VI = U2D_GRID%CELL_NODES(I,IC)
               PSII = XP*U2D_GRID%BASIS_COEFFS(1,I,IC) + YP*U2D_GRID%BASIS_COEFFS(2,I,IC) + U2D_GRID%BASIS_COEFFS(3,I,IC)
               PHI = PHI + PSII*PHI_FIELD(VI)
            END DO

         ELSE IF (DIMS == 3) THEN

            IC = particles(JP)%IC

            XP = particles(JP)%X
            YP = particles(JP)%Y
            ZP = particles(JP)%Z

            PHI = 0.d0
            DO I = 1, 4
               VI = U3D_GRID%CELL_NODES(I,IC)
               PSII = XP*U3D_GRID%BASIS_COEFFS(1,I,IC) &
                    + YP*U3D_GRID%BASIS_COEFFS(2,I,IC) &
                    + ZP*U3D_GRID%BASIS_COEFFS(3,I,IC) &
                       + U3D_GRID%BASIS_COEFFS(4,I,IC)
               PHI = PHI + PSII*PHI_FIELD(VI)
            END DO

         END IF
      ELSE
         IF (DIMS == 2) THEN
            CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
            PHI = WEIGHTS(1)*PHI_FIELD(INDICES(1)+1) + &
                  WEIGHTS(2)*PHI_FIELD(INDICES(2)+1) + &
                  WEIGHTS(3)*PHI_FIELD(INDICES(3)+1) + &
                  WEIGHTS(4)*PHI_FIELD(INDICES(4)+1)
         ELSE
            CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
            PHI = WEIGHTS(1)*PHI_FIELD(INDICES(1)+1) + &
                  WEIGHTS(2)*PHI_FIELD(INDICES(2)+1)
         END IF
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
               POINT = U2D_GRID%NODE_COORDS(:, IN)
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
                        
                        B_FIELD(:, 1, IN) = B_FIELD(:, 1, IN) + SOLENOIDS(ICOIL)%WIRE_CURRENT*MU0/(4*PI)*CROSS(DL, RPRIME) / &
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




END MODULE fields
