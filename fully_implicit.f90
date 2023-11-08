
MODULE fully_implicit

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
   
   USE global
   USE screen
   USE tools
   USE grid_and_partition

   
   USE petsc
   USE petscksp

   IMPLICIT NONE

   Mat Amat, Jmat, Jmfmat, Qmat, Rmat, Qmatfact, Pmat
   Vec bvec, xvec, x_seq, solvec_seq, xvec_seq, rvec, solvec
   VecScatter ctx
   KSP ksp, kspnk
   PetscInt one, maxit, maxf
   PetscInt Istart, Iend
   PetscReal val, norm, f0, stol, rtol, abstol
   PetscScalar, POINTER :: PHI_FIELD_TEMP(:)
   !PetscScalar, POINTER :: PHIBAR_FIELD(:)
   !PetscScalar, POINTER :: PHI_FIELD_OLD(:)
   KSPConvergedReason reason
   SNESConvergedReason snesreason
   PC pc, pcnk
   SNES snes
   SNESLineSearch linesearch
   PetscViewer viewer
   IS rowperm, colperm
   MatFactorInfo  info(MAT_FACTORINFO_SIZE)

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_NE, CELL_TE

   CONTAINS
   

   ! Reminder of JACOBIAN_TYPE:
   ! 0 --> Jacobian containing only the FEM-related derivatives (from Poisson's equation and boundary conditions.), no particle term.
   ! 1 --> Exact Jacobian from following all particles and taking all partial derivatives
   ! 2 --> Approximate Jacobian neglecting d(Dt) terms
   ! 3 --> ECSIM-like mass matrices. Accumulation is done with particles before advection
   ! 4 --> ECSIM-like mass matrices. Accumulation is done with particles after advection
   ! 5 --> The Chacon suggested one. Particle moments taken before advection. Requires a matrix inversion that at the moment does not work.
   ! 6 --> Boltzmann's relation (with fixed density and temperature), currently not working.


   SUBROUTINE SOLVE_POISSON_FULLY_IMPLICIT

      LOGICAL :: SET_SOLVEC_TO_ZERO = .FALSE.
      PetscScalar, POINTER :: solvec_l(:)
      PetscBool :: flg

      CALL SNESCreate(PETSC_COMM_WORLD,snes,ierr)
      !CALL SNESSetType(snes, SNESNEWTONLS ,ierr) ! defalut
      !CALL SNESGetLineSearch(snes,linesearch,ierr)
      !CALL SNESLineSearchSetType(linesearch,SNESLINESEARCHBT,ierr)
      !CALL SNESLineSearchSetOrder(linesearch,SNES_LINESEARCH_ORDER_QUADRATIC,ierr)

      CALL VecCreate(PETSC_COMM_WORLD,rvec,ierr)
      CALL VecSetSizes(rvec,PETSC_DECIDE,NNODES,ierr)
      CALL VecSetFromOptions(rvec, ierr)
      CALL VecDuplicate(rvec, solvec, ierr)


      IF (RESIDUAL_AND_JACOBIAN_COMBINED) THEN
         CALL SNESSetFunction(snes,rvec,FormFunctionAndJacobian,0,ierr) ! Function and Jacobian computed together
      ELSE
         CALL SNESSetFunction(snes,rvec,FormFunction,0,ierr) ! Function and Jacobian computed separately
      END IF



      !CALL MatSetOption(Jmat,MAT_SPD,PETSC_TRUE,ierr)
      !CALL MatMPIAIJSetPreallocation(Jmat,30,PETSC_NULL_INTEGER,30,PETSC_NULL_INTEGER,ierr) ! DBDBDBDBDBDB Large preallocation!
      !CALL MatSetFromOptions(Jmat,ierr)
      !CALL MatSetUp(Jmat,ierr)

      !CALL MatCreate(PETSC_COMM_WORLD,Jmfmat,ierr)
      !CALL MatSetSizes(Jmfmat,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
      !CALL MatSetType(Jmfmat,MATSHELL,ierr)
      !CALL MatSetUp(Jmfmat,ierr)


      CALL MatCreate(PETSC_COMM_WORLD,Jmat,ierr)
      CALL MatSetSizes(Jmat,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
      CALL MatSetType(Jmat, MATMPIAIJ, ierr)

      CALL MatCreate(PETSC_COMM_WORLD,Pmat,ierr)
      CALL MatSetSizes(Pmat,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
      CALL MatSetType(Pmat, MATMPIAIJ, ierr)

      IF (RESIDUAL_AND_JACOBIAN_COMBINED) THEN
         CALL SNESSetJacobian(snes,Jmat,Jmat,PETSC_NULL_FUNCTION,0,ierr)   ! Use this for the combined computation of residual and Jacobian.
      ELSE
         flg  = PETSC_FALSE
         CALL PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-snes_mf_operator",flg,ierr)
         IF (flg) THEN  ! We want only the preconditioner to be filled. The Jacobian is computed from finite differencing.
            CALL SNESSetJacobian(snes,Jmat,Pmat,FormJacobian,0,ierr) ! The expensive but safe one. Jacobian computed independently.
         ELSE
            CALL SNESSetJacobian(snes,Jmat,Jmat,FormJacobian,0,ierr) ! The expensive but safe one. Jacobian computed independently.
         END IF
      END IF


      !abstol = 1.d-5
      !rtol = 1.d-5
      !stol = 1.d0
      !maxit = 10000 !10
      !maxf = 10000 !30
      !CALL SNESSetTolerances(snes, abstol, rtol, stol, maxit, maxf, ierr)
      !CALL SNESSetTolerances(snes, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, maxit, maxf, ierr)
      !CALL SNESSetTolerances(snes, PETSC_DEFAULT_REAL, SNES_RTOL, PETSC_DEFAULT_REAL, maxit, maxf, ierr)

      !CALL SNESGetKSP(snes,kspnk,ierr)
      !CALL KSPGetPC(kspnk,pcnk,ierr)
      !CALL PCSetType(pcnk,PCLU,ierr)

      !CALL KSPSetTolerances(kspnk,PETSC_DEFAULT_REAL,1.d-3,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
      !CALL KSPSetTolerances(ksp,1.e-4,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,20,ierr)

      CALL SNESSetFromOptions(snes,ierr)


      !  Note: The user should initialize the vector, x, with the initial guess
      !  for the nonlinear solver prior to calling SNESSolve().  In particular,
      !  to employ an initial guess of zero, the user should explicitly set
      !  this vector to zero by calling VecSet().

      IF (SET_SOLVEC_TO_ZERO) THEN
         CALL VecSet(solvec,0.d0,ierr)
      ELSE
         CALL VecGetOwnershipRange(solvec,Istart,Iend,ierr)
         CALL VecGetArrayF90(solvec,solvec_l,ierr)
         solvec_l = PHI_FIELD(Istart+1:Iend)
         CALL VecRestoreArrayF90(solvec,solvec_l,ierr)
      END IF
      
      ! Test the jacobian
      !CALL MatCreate(PETSC_COMM_WORLD,testJ,ierr)
      !CALL MatSetSizes(testJ,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
      !CALL MatSetFromOptions(testJ,ierr)
      !CALL MatSetUp(testJ,ierr)
      !CALL SNESComputeJacobianDefault(snes,solvec,testJ,testJ,ierr)
      !CALL MatView(testJ,PETSC_VIEWER_STDOUT_WORLD,ierr)

      IF (JACOBIAN_TYPE == 3) CALL COMPUTE_MASS_MATRICES(particles)
      IF (JACOBIAN_TYPE == 5) CALL COMPUTE_DENSITY_TEMPERATURE(particles)

      CALL SNESSolve(snes,PETSC_NULL_VEC,solvec,ierr)
      CALL SNESGetConvergedReason(snes,snesreason,ierr)
      IF (PROC_ID == 0) WRITE(*,*) 'SNESConvergedReason = ', snesreason
      !CALL VecView(solvec,PETSC_VIEWER_STDOUT_WORLD,ierr)
      !IF (PROC_ID == 0) WRITE(*,*) 'PHI_FIELD was: ', PHI_FIELD

      CALL VecScatterCreateToAll(solvec,ctx,solvec_seq,ierr)
      CALL VecScatterBegin(ctx,solvec,solvec_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterEnd(ctx,solvec,solvec_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterDestroy(ctx, ierr)

      CALL VecGetArrayReadF90(solvec_seq,PHI_FIELD_TEMP,ierr)
      IF (ALLOCATED(PHIBAR_FIELD)) DEALLOCATE(PHIBAR_FIELD)
      ALLOCATE(PHIBAR_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(solvec_seq,PHI_FIELD_TEMP,ierr)

      PHI_FIELD = 2.*PHIBAR_FIELD - PHI_FIELD


      ! Cleanup.
      CALL VecDestroy(solvec, ierr)
      CALL VecDestroy(rvec, ierr)
      CALL VecDestroy(solvec_seq,ierr)
      CALL SNESDestroy(snes, ierr)
      CALL MatDestroy(Jmat,ierr)
      CALL MatDestroy(Pmat,ierr)

   END SUBROUTINE SOLVE_POISSON_FULLY_IMPLICIT


   SUBROUTINE FormJacobianVectorProduct(jac, x, y)

      Mat jac
      Vec x, y

      IF (PROC_ID == 0) THEN
         WRITE(*,*) 'FormJacobianVectorProduct Called'
      END IF



   END SUBROUTINE FormJacobianVectorProduct


   SUBROUTINE FormJacobianVectorProduct_old(jac, x, y)


      Mat jac
      Vec x, y
      INTEGER :: I

      INTEGER :: V1, V2, V3
      INTEGER :: P, Q, VP, VQ
      REAL(KIND=8) :: KPQ, VOLUME, AREA

      IF (PROC_ID == 0) THEN
         WRITE(*,*) 'FormJacobianVectorProduct Called'
      END IF

      CALL MatCreate(PETSC_COMM_WORLD,jac,ierr)
      CALL MatSetSizes(jac,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
      CALL MatSetType(jac, MATMPIAIJ, ierr)

      CALL COMPUTE_MASS_MATRICES(particles)

      CALL MatGetOwnershipRange( jac, Istart, Iend, ierr)

      ! Accumulate Jacobian. All this is in principle not needed since we already have Amat.

      IF (DIMS == 2) THEN
         DO I = 1, NCELLS
            AREA = CELL_AREAS(I)
            
            ! We need to ADD to a sparse matrix entry.
            DO P = 1, 3
               VP = U2D_GRID%CELL_NODES(P,I)
               IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(VP-1)) THEN
                     DO Q = 1, 3
                        VQ = U2D_GRID%CELL_NODES(Q,I)
                        KPQ = AREA*(U2D_GRID%BASIS_COEFFS(1,P,I)*U2D_GRID%BASIS_COEFFS(1,Q,I) &
                                  + U2D_GRID%BASIS_COEFFS(2,P,I)*U2D_GRID%BASIS_COEFFS(2,Q,I))
                        
                        IF (AXI) THEN
                           V1 = U2D_GRID%CELL_NODES(1,I)
                           V2 = U2D_GRID%CELL_NODES(2,I)
                           V3 = U2D_GRID%CELL_NODES(3,I)  
                           KPQ = KPQ*(U2D_GRID%NODE_COORDS(2, V1) &
                                    + U2D_GRID%NODE_COORDS(2, V2) &
                                    + U2D_GRID%NODE_COORDS(2, V3))/3.
                        END IF

                        IF (JACOBIAN_TYPE == 3 .OR. JACOBIAN_TYPE == 4) THEN
                           KPQ = KPQ * (MASS_MATRIX(I)+1)
                        END IF

                        CALL MatSetValues(jac,one,VP-1,one,VQ-1,-2.*KPQ,ADD_VALUES,ierr)
                     END DO
                  END IF
               END IF
            END DO

         END DO
      ELSE IF (DIMS == 3) THEN
         DO I = 1, NCELLS
            VOLUME = CELL_VOLUMES(I)
            
            ! We need to ADD to a sparse matrix entry.
            DO P = 1, 4
               VP = U3D_GRID%CELL_NODES(P,I)
               IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(VP-1)) THEN
                     DO Q = 1, 4
                        VQ = U3D_GRID%CELL_NODES(Q,I)
                        KPQ = VOLUME*(U3D_GRID%BASIS_COEFFS(1,P,I)*U3D_GRID%BASIS_COEFFS(1,Q,I) &
                                    + U3D_GRID%BASIS_COEFFS(2,P,I)*U3D_GRID%BASIS_COEFFS(2,Q,I) &
                                    + U3D_GRID%BASIS_COEFFS(3,P,I)*U3D_GRID%BASIS_COEFFS(3,Q,I))

                        IF (JACOBIAN_TYPE == 3 .OR. JACOBIAN_TYPE == 4) THEN
                           KPQ = KPQ * (MASS_MATRIX(I)+1)
                        END IF

                        CALL MatSetValues(jac,one,VP-1,one,VQ-1,-2.*KPQ,ADD_VALUES,ierr)
                     END DO
                  END IF
               END IF
            END DO

         END DO
      END IF

      CALL MatAssemblyBegin(jac,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(jac,MAT_FLUSH_ASSEMBLY,ierr)

      DO I = Istart, Iend-1
         IF (IS_DIRICHLET(I)) CALL MatSetValues(jac,one,I,one,I,-2.d0,INSERT_VALUES,ierr)
      END DO

      CALL MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)

      CALL MatMult(jac, x, y, ierr)

   END SUBROUTINE FormJacobianVectorProduct_old


   SUBROUTINE FormFunction(snes,x,f,dummy,ierr_l)

      IMPLICIT NONE

      SNES snes
      Vec x,f
      PetscErrorCode ierr_l
      INTEGER dummy(*)
      PetscScalar, POINTER :: RESIDUAL(:)
      CHARACTER(LEN=512)  :: filename

      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, K11, K22, K33, K12, K23, K13, AREA
      INTEGER :: V1, V2, V3, I
      INTEGER :: P, Q, VP, VQ
      REAL(KIND=8) :: KPQ, VOLUME

      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: part_adv

      !CALL VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)

      IF (PROC_ID == 0) THEN
         WRITE(*,*) 'FormFunction Called'
      END IF

      CALL VecScatterCreateToAll(x,ctx,x_seq,ierr)
      CALL VecScatterBegin(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterEnd(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterDestroy(ctx, ierr)

      ! CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD,ierr)
      CALL VecGetArrayReadF90(x_seq,PHI_FIELD_TEMP,ierr)
      IF (ALLOCATED(PHIBAR_FIELD)) DEALLOCATE(PHIBAR_FIELD)
      ALLOCATE(PHIBAR_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(x_seq,PHI_FIELD_TEMP,ierr)
      CALL VecDestroy(x_seq,ierr)

      ! Compute the RHS corresponding to PHIBAR_FIELD
      PHI_FIELD_NEW = 2.*PHIBAR_FIELD - PHI_FIELD

      ALLOCATE(RHS_NEW, SOURCE = RHS)
      RHS_NEW = 0.d0

      IF (DIMS == 2) THEN
         DO I = 1, NCELLS
            AREA = CELL_AREAS(I)
            V1 = U2D_GRID%CELL_NODES(1,I)
            V2 = U2D_GRID%CELL_NODES(2,I)
            V3 = U2D_GRID%CELL_NODES(3,I)            
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
            IF (AXI) THEN
               K11 = K11*(Y1+Y2+Y3)/3.
               K22 = K22*(Y1+Y2+Y3)/3.
               K33 = K33*(Y1+Y2+Y3)/3.
               K12 = K12*(Y1+Y2+Y3)/3.
               K23 = K23*(Y1+Y2+Y3)/3.
               K13 = K13*(Y1+Y2+Y3)/3.
            END IF

            RHS_NEW(V1-1) = RHS_NEW(V1-1) + K11*PHI_FIELD_NEW(V1) + K12*PHI_FIELD_NEW(V2) + K13*PHI_FIELD_NEW(V3)
            RHS_NEW(V2-1) = RHS_NEW(V2-1) + K12*PHI_FIELD_NEW(V1) + K22*PHI_FIELD_NEW(V2) + K23*PHI_FIELD_NEW(V3)
            RHS_NEW(V3-1) = RHS_NEW(V3-1) + K13*PHI_FIELD_NEW(V1) + K23*PHI_FIELD_NEW(V2) + K33*PHI_FIELD_NEW(V3)
         END DO
      ELSE IF (DIMS == 3) THEN
         DO I = 1, NCELLS
            VOLUME = CELL_VOLUMES(I)
            DO P = 1, 4
               VP = U3D_GRID%CELL_NODES(P,I)
               DO Q = 1, 4
                  VQ = U3D_GRID%CELL_NODES(Q,I)
                  KPQ = VOLUME*(U3D_GRID%BASIS_COEFFS(1,P,I)*U3D_GRID%BASIS_COEFFS(1,Q,I) &
                              + U3D_GRID%BASIS_COEFFS(2,P,I)*U3D_GRID%BASIS_COEFFS(2,Q,I) &
                              + U3D_GRID%BASIS_COEFFS(3,P,I)*U3D_GRID%BASIS_COEFFS(3,Q,I))
                  RHS_NEW(VQ-1) = RHS_NEW(VQ-1) + KPQ*PHI_FIELD_NEW(VP)
               END DO
            END DO
         END DO 
      END IF

      DO I = 0, NNODES-1
         IF (IS_DIRICHLET(I)) THEN
            !RHS_NEW(I) = DIRICHLET(I)
            RHS_NEW(I) = PHI_FIELD_NEW(I+1)
         ELSE IF (IS_NEUMANN(I)) THEN
            RHS_NEW(I) = RHS_NEW(I) + NEUMANN(I)
         END IF
      END DO



      CALL COMPUTE_E_FIELD

      !  Advect the particles using the guessed new potential
      ALLOCATE(part_adv, SOURCE = particles)
      CALL TIMER_START(3)
      CALL ADVECT_CN(part_adv, .FALSE., .FALSE., Jmat)
      CALL TIMER_STOP(3)
      CALL DEPOSIT_CHARGE(part_adv)
      DEALLOCATE(part_adv)

      ! Compute the new potential from the charge distribution
      !CALL SOLVE_POISSON

      ! Compute the residual
      CALL VecGetOwnershipRange(f,Istart,Iend,ierr)

      CALL VecGetArrayF90(f,RESIDUAL,ierr_l)
      RESIDUAL = RHS(Istart:Iend-1) - RHS_NEW(Istart:Iend-1)
      !RESIDUAL = PHI_FIELD(Istart+1:Iend) + PHI_FIELD_OLD(Istart+1:Iend) - 2.*PHIBAR_FIELD(Istart+1:Iend)
      CALL VecRestoreArrayF90(f,RESIDUAL,ierr_l)


      CALL VecNorm(f,NORM_2,norm,ierr)
      IF (PROC_ID == 0) THEN
         !WRITE(*,*) ' PHI_FIELD = ', PHI_FIELD
         WRITE(*,*) '||RESIDUAL|| = ', norm !, ' with potential ', PHIBAR_FIELD

         WRITE(filename, "(A,A)") TRIM(ADJUSTL(RESIDUAL_SAVE_PATH)), "residuals" ! Compose filename   
         OPEN(66331, FILE=filename, POSITION='append', STATUS='unknown', ACTION='write')
         WRITE(66331,*) tID, norm
         CLOSE(66331)
      END IF

      DEALLOCATE(RHS_NEW)

   END SUBROUTINE FormFunction



   SUBROUTINE FormJacobian(snes,x,jactofill,precond,dummy,ierr_l)

      IMPLICIT NONE

      SNES snes
      Vec x
      Mat  jac, prec, jactofill, precond
      PetscErrorCode ierr_l
      PetscBool      flg
      INTEGER dummy(*)
      INTEGER I

      REAL(KIND=8) :: AREA
      INTEGER :: V1, V2, V3
      INTEGER :: P, Q, VP, VQ
      REAL(KIND=8) :: KPQ, VOLUME, VALUETOADD
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: part_adv
      REAL(KIND=8) :: MPQ, FACTOR1, FACTOR2, ME
      INTEGER :: ELECTRON_S_ID

      IF (PROC_ID == 0) THEN
         WRITE(*,*) 'FormJacobian Called'
      END IF


      flg  = PETSC_FALSE
      CALL PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-snes_mf_operator",flg,ierr)
      IF (flg) THEN  ! We want only the preconditioner to be filled. The Jacobian is computed from finite differencing.
         jac  = precond
         prec = jactofill
      ELSE ! We want the Jacobian to be filled. The preconditioner should be identically filled.
         jac  = jactofill
         prec = precond
      END IF

      CALL MatMPIAIJSetPreallocation(jac,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER,ierr) ! DBDBDBDBDBDB Large preallocation!
      CALL MatSetFromOptions(jac,ierr)
      CALL MatSetUp(jac,ierr)


      !CALL MatConvert(Amat,MATSAME,MAT_INITIAL_MATRIX,jac,ierr)
      !mult = -2.d0
      !CALL MatScale(jac,mult,ierr)
      CALL MatZeroEntries(jac,ierr)
      !CALL MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      !CALL MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
      !WRITE(*,*) 'We got here 1 on proc ', PROC_ID
      !CALL MatAXPY(jac,-2.d0,Amat,UNKNOWN_NONZERO_PATTERN,ierr)
      !WRITE(*,*) 'We got here 2 on proc ', PROC_ID
      !CALL MatDuplicate(Amat,MAT_COPY_VALUES,jac,ierr)
      !CALL MatScale(jac,mult,ierr)
      CALL MatAssemblyBegin(jac,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(jac,MAT_FLUSH_ASSEMBLY,ierr)

      CALL VecScatterCreateToAll(x,ctx,x_seq,ierr)
      CALL VecScatterBegin(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterEnd(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterDestroy(ctx, ierr)

      ! CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD,ierr)
      CALL VecGetArrayReadF90(x_seq,PHI_FIELD_TEMP,ierr)
      IF (ALLOCATED(PHIBAR_FIELD)) DEALLOCATE(PHIBAR_FIELD)
      ALLOCATE(PHIBAR_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(x_seq,PHI_FIELD_TEMP,ierr)
      CALL VecDestroy(x_seq,ierr)

      CALL COMPUTE_E_FIELD

      !CALL MatView(jac,PETSC_VIEWER_STDOUT_WORLD,ierr)
      !  Advect the particles using the guessed new potential
      ALLOCATE(part_adv, SOURCE = particles)
      CALL TIMER_START(3)
      CALL ADVECT_CN(part_adv, .FALSE., .TRUE., jac)
      CALL TIMER_STOP(3)
      IF (JACOBIAN_TYPE == 4) CALL COMPUTE_MASS_MATRICES(part_adv)
      IF (JACOBIAN_TYPE == 6) CALL COMPUTE_DENSITY_TEMPERATURE(part_adv)
      IF (JACOBIAN_TYPE == 7) CALL COMPUTE_DENSITY_TEMPERATURE(part_adv)
      DEALLOCATE(part_adv)


      CALL MatGetOwnershipRange( jac, Istart, Iend, ierr)
      IF (JACOBIAN_TYPE == 6) THEN
         DO I = 1, NCELLS
            AREA = CELL_AREAS(I)
            DO P = 1, 3
               VP = U2D_GRID%CELL_NODES(P,I)
               IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(VP-1)) THEN
                     VALUETOADD = -QE*QE*5.d11/(EPS0*KB*11600.d0)*AREA/3.
                     CALL MatSetValues(jac,one,VP-1,one,VP-1,VALUETOADD,ADD_VALUES,ierr)
                  END IF
               END IF
            END DO
         END DO
      END IF

      IF (JACOBIAN_TYPE == 7) THEN
         DO I = 1, NCELLS
            AREA = CELL_AREAS(I)
            DO P = 1, 3
               VP = U2D_GRID%CELL_NODES(P,I)
               IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(VP-1)) THEN
                     VALUETOADD = -QE*QE*CELL_NE(I)/(EPS0*KB*CELL_TE(I))*AREA/3.
                     CALL MatSetValues(jac,one,VP-1,one,VP-1,VALUETOADD,ADD_VALUES,ierr)
                  END IF
               END IF
            END DO
         END DO
      END IF
      !CALL MatAssemblyBegin(jac,MAT_FLUSH_ASSEMBLY,ierr)
      !CALL MatAssemblyEnd(jac,MAT_FLUSH_ASSEMBLY,ierr)

      !DO I = 0, U2D_GRID%NUM_NODES - 1
      !   IF (IS_DIRICHLET(I)) THEN
      !      CALL MatSetValue(jac,I,I,-2.d0,INSERT_VALUES,ierr)
      !   END IF
      !END DO

      
      CALL MatAssemblyBegin(jac,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(jac,MAT_FLUSH_ASSEMBLY,ierr)

      !CALL MatSetValues(jac,1,0,1,9,1.234d0,ADD_VALUES,ierr)

      !CALL MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      !CALL MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)









      CALL MatGetOwnershipRange( jac, Istart, Iend, ierr)




      IF (JACOBIAN_TYPE == 5) THEN

         CALL MatCreate(PETSC_COMM_WORLD,Qmat,ierr)
         CALL MatSetSizes(Qmat,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
         CALL MatSetType(Qmat, MATMPIAIJ, ierr)
         CALL MatMPIAIJSetPreallocation(Qmat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER,ierr) ! DBDBDBDBDBDB Large preallocation!
         CALL MatSetFromOptions(Qmat,ierr)
         CALL MatSetUp(Qmat,ierr)

         CALL MatCreate(PETSC_COMM_WORLD,Rmat,ierr)
         CALL MatSetSizes(Rmat,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
         CALL MatSetType(Rmat, MATMPIAIJ, ierr)
         CALL MatMPIAIJSetPreallocation(Rmat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER,ierr) ! DBDBDBDBDBDB Large preallocation!
         CALL MatSetFromOptions(Rmat,ierr)
         CALL MatSetUp(Rmat,ierr)

         ELECTRON_S_ID = SPECIES_NAME_TO_ID('e')
         ME = SPECIES(ELECTRON_S_ID)%MOLECULAR_MASS
         IF (DIMS == 2) THEN
            DO I = 1, NCELLS
               AREA = CELL_AREAS(I)
               FACTOR1 = 0.25*DT*DT*KB*CELL_TE(I)/ME
               FACTOR2 = 0.25*DT*DT*QE*QE*CELL_NE(I)/EPS0/ME
               ! We need to ADD to a sparse matrix entry.
               DO P = 1, 3
                  VP = U2D_GRID%CELL_NODES(P,I)
                  IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                     DO Q = 1, 3
                        VQ = U2D_GRID%CELL_NODES(Q,I)
                        IF (Q==P) THEN
                           MPQ = AREA/6.
                        ELSE
                           MPQ = AREA/12.
                        END IF
                        KPQ = AREA*(U2D_GRID%BASIS_COEFFS(1,P,I)*U2D_GRID%BASIS_COEFFS(1,Q,I) &
                                  + U2D_GRID%BASIS_COEFFS(2,P,I)*U2D_GRID%BASIS_COEFFS(2,Q,I))
                        
                        CALL MatSetValues(Qmat,one,VP-1,one,VQ-1,MPQ,ADD_VALUES,ierr)
                        CALL MatSetValues(Qmat,one,VP-1,one,VQ-1,KPQ*FACTOR1,ADD_VALUES,ierr)
                        CALL MatSetValues(Rmat,one,VP-1,one,VQ-1,KPQ*FACTOR2,ADD_VALUES,ierr)
                     END DO
                  END IF
               END DO
   
            END DO
         ELSE
            CALL ERROR_ABORT('Not implemented!')
         END IF

         CALL MatAssemblyBegin(Qmat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(Qmat,MAT_FINAL_ASSEMBLY,ierr)

         CALL MatAssemblyBegin(Rmat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(Rmat,MAT_FINAL_ASSEMBLY,ierr)

         ! CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'QMatrix', FILE_MODE_WRITE, viewer, ierr)
         ! CALL MatView(Qmat,viewer,ierr)
         ! CALL PetscViewerDestroy(viewer,ierr)
         ! WRITE(*,*) 'Q written to file.'

         ! CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'RMatrix', FILE_MODE_WRITE, viewer, ierr)
         ! CALL MatView(Rmat,viewer,ierr)
         ! CALL PetscViewerDestroy(viewer,ierr)
         ! WRITE(*,*) 'R written to file.'

         ! CALL SLEEP(10)

      END IF


      ! Accumulate Jacobian. All this is in principle not needed since we already have Amat.

      IF (DIMS == 2) THEN
         DO I = 1, NCELLS
            AREA = CELL_AREAS(I)
            
            ! We need to ADD to a sparse matrix entry.
            DO P = 1, 3
               VP = U2D_GRID%CELL_NODES(P,I)
               IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(VP-1)) THEN
                     DO Q = 1, 3
                        VQ = U2D_GRID%CELL_NODES(Q,I)
                        KPQ = AREA*(U2D_GRID%BASIS_COEFFS(1,P,I)*U2D_GRID%BASIS_COEFFS(1,Q,I) &
                                  + U2D_GRID%BASIS_COEFFS(2,P,I)*U2D_GRID%BASIS_COEFFS(2,Q,I))
                        
                        IF (AXI) THEN
                           V1 = U2D_GRID%CELL_NODES(1,I)
                           V2 = U2D_GRID%CELL_NODES(2,I)
                           V3 = U2D_GRID%CELL_NODES(3,I)  
                           KPQ = KPQ*(U2D_GRID%NODE_COORDS(2, V1) &
                                    + U2D_GRID%NODE_COORDS(2, V2) &
                                    + U2D_GRID%NODE_COORDS(2, V3))/3.
                        END IF

                        IF (JACOBIAN_TYPE == 3 .OR. JACOBIAN_TYPE == 4) THEN
                           KPQ = KPQ * (MASS_MATRIX(I)+1)
                        END IF

                        CALL MatSetValues(jac,one,VP-1,one,VQ-1,-2.*KPQ,ADD_VALUES,ierr)
                     END DO
                  END IF
               END IF
            END DO

         END DO
      ELSE IF (DIMS == 3) THEN
         DO I = 1, NCELLS
            VOLUME = CELL_VOLUMES(I)
            
            ! We need to ADD to a sparse matrix entry.
            DO P = 1, 4
               VP = U3D_GRID%CELL_NODES(P,I)
               IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(VP-1)) THEN
                     DO Q = 1, 4
                        VQ = U3D_GRID%CELL_NODES(Q,I)
                        KPQ = VOLUME*(U3D_GRID%BASIS_COEFFS(1,P,I)*U3D_GRID%BASIS_COEFFS(1,Q,I) &
                                    + U3D_GRID%BASIS_COEFFS(2,P,I)*U3D_GRID%BASIS_COEFFS(2,Q,I) &
                                    + U3D_GRID%BASIS_COEFFS(3,P,I)*U3D_GRID%BASIS_COEFFS(3,Q,I))

                        IF (JACOBIAN_TYPE == 3 .OR. JACOBIAN_TYPE == 4) THEN
                           KPQ = KPQ * (MASS_MATRIX(I)+1)
                        END IF

                        CALL MatSetValues(jac,one,VP-1,one,VQ-1,-2.*KPQ,ADD_VALUES,ierr)
                     END DO
                  END IF
               END IF
            END DO

         END DO
      END IF

      CALL MatAssemblyBegin(jac,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(jac,MAT_FLUSH_ASSEMBLY,ierr)

      DO I = Istart, Iend-1
         IF (IS_DIRICHLET(I)) CALL MatSetValues(jac,one,I,one,I,-2.d0,INSERT_VALUES,ierr)
      END DO





      !CALL MatAXPY(jac,mult,Amat,DIFFERENT_NONZERO_PATTERN,ierr)
      !IF (PROC_ID == 0) WRITE(*,*) 'jac matrix'
      !CALL MatView(jac,PETSC_VIEWER_STDOUT_WORLD,ierr)
      !IF (PROC_ID == 0) WRITE(*,*) 'Amat matrix'
      !CALL MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr)
      CALL MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)

      flg  = PETSC_FALSE
      CALL PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-snes_mf_operator",flg,ierr)
      IF (flg) THEN  ! We want only the preconditioner to be filled. The Jacobian is computed from finite differencing.
         CALL MatSetFromOptions(prec,ierr)
         CALL MatSetUp(prec,ierr)
         CALL MatAssemblyBegin(prec,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(prec,MAT_FINAL_ASSEMBLY,ierr)
      END IF



      !jac = prec
 
      ! CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'JacobianMatrix', FILE_MODE_WRITE, viewer, ierr)
      ! CALL PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER, &
      !                          PETSC_NULL_CHARACTER,0,0,684,684,viewer,ierr)
      ! CALL MatView(prec,viewer,ierr)
      ! CALL PetscSleep(5.d0,ierr)
      ! CALL PetscViewerDestroy(viewer,ierr)

      ! CALL PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER, &
      ! PETSC_NULL_CHARACTER,0,0,684,684,viewer,ierr)
      ! CALL MatView(jac,viewer,ierr)
      ! CALL PetscSleep(5.d0,ierr)
      ! CALL PetscViewerDestroy(viewer,ierr)
      ! WRITE(*,*) 'Jacobian written to file.'
      ! CALL SLEEP(10)


   END SUBROUTINE FormJacobian




   SUBROUTINE FormFunctionAndJacobian(snes,x,f,dummy,ierr_l)

      ! SNES variable x is the mid-step potential \phi^{n+1/2}.

      IMPLICIT NONE

      SNES snes
      Vec x,f


      PetscErrorCode ierr_l
      INTEGER dummy(*)
      PetscScalar, POINTER :: RESIDUAL(:)
      CHARACTER(LEN=512)  :: filename

      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, K11, K22, K33, K12, K23, K13, AREA
      INTEGER :: V1, V2, V3, I
      INTEGER :: P, Q, VP, VQ
      REAL(KIND=8) :: KPQ, VOLUME
      REAL(KIND=8) :: MPQ, FACTOR1, FACTOR2, ME
      INTEGER :: ELECTRON_S_ID

      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: part_adv

      !CALL VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)

      IF (PROC_ID == 0) THEN
         WRITE(*,*) 'FormFunctionAndJacobian Called'
      END IF



      CALL MatMPIAIJSetPreallocation(Jmat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER,ierr) ! DBDBDBDBDBDB Large preallocation!
      CALL MatSetFromOptions(Jmat,ierr)
      CALL MatSetUp(Jmat,ierr)
      CALL MatZeroEntries(Jmat,ierr)
      CALL MatAssemblyBegin(Jmat,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Jmat,MAT_FLUSH_ASSEMBLY,ierr)



      CALL VecScatterCreateToAll(x,ctx,x_seq,ierr)
      CALL VecScatterBegin(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterEnd(ctx,x,x_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterDestroy(ctx, ierr)

      ! CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD,ierr)
      CALL VecGetArrayReadF90(x_seq,PHI_FIELD_TEMP,ierr)
      IF (ALLOCATED(PHIBAR_FIELD)) DEALLOCATE(PHIBAR_FIELD)
      ALLOCATE(PHIBAR_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(x_seq,PHI_FIELD_TEMP,ierr)
      CALL VecDestroy(x_seq,ierr)

      ! Compute the RHS corresponding to PHIBAR_FIELD
      PHI_FIELD_NEW = 2.*PHIBAR_FIELD - PHI_FIELD

      ALLOCATE(RHS_NEW, SOURCE = RHS)
      RHS_NEW = 0.d0

      IF (DIMS == 2) THEN
         DO I = 1, NCELLS
            AREA = CELL_AREAS(I)
            V1 = U2D_GRID%CELL_NODES(1,I)
            V2 = U2D_GRID%CELL_NODES(2,I)
            V3 = U2D_GRID%CELL_NODES(3,I)            
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
            IF (AXI) THEN
               K11 = K11*(Y1+Y2+Y3)/3.
               K22 = K22*(Y1+Y2+Y3)/3.
               K33 = K33*(Y1+Y2+Y3)/3.
               K12 = K12*(Y1+Y2+Y3)/3.
               K23 = K23*(Y1+Y2+Y3)/3.
               K13 = K13*(Y1+Y2+Y3)/3.
            END IF

            RHS_NEW(V1-1) = RHS_NEW(V1-1) + K11*PHI_FIELD_NEW(V1) + K12*PHI_FIELD_NEW(V2) + K13*PHI_FIELD_NEW(V3)
            RHS_NEW(V2-1) = RHS_NEW(V2-1) + K12*PHI_FIELD_NEW(V1) + K22*PHI_FIELD_NEW(V2) + K23*PHI_FIELD_NEW(V3)
            RHS_NEW(V3-1) = RHS_NEW(V3-1) + K13*PHI_FIELD_NEW(V1) + K23*PHI_FIELD_NEW(V2) + K33*PHI_FIELD_NEW(V3)
         END DO
      ELSE IF (DIMS == 3) THEN
         DO I = 1, NCELLS
            VOLUME = CELL_VOLUMES(I)
            DO P = 1, 4
               VP = U3D_GRID%CELL_NODES(P,I)
               DO Q = 1, 4
                  VQ = U3D_GRID%CELL_NODES(Q,I)
                  KPQ = VOLUME*(U3D_GRID%BASIS_COEFFS(1,P,I)*U3D_GRID%BASIS_COEFFS(1,Q,I) &
                              + U3D_GRID%BASIS_COEFFS(2,P,I)*U3D_GRID%BASIS_COEFFS(2,Q,I) &
                              + U3D_GRID%BASIS_COEFFS(3,P,I)*U3D_GRID%BASIS_COEFFS(3,Q,I))
                  RHS_NEW(VQ-1) = RHS_NEW(VQ-1) + KPQ*PHI_FIELD_NEW(VP)
               END DO
            END DO
         END DO 
      END IF

      DO I = 0, NNODES-1
         IF (IS_DIRICHLET(I)) THEN
            !RHS_NEW(I) = DIRICHLET(I)
            RHS_NEW(I) = PHI_FIELD_NEW(I+1)
         ELSE IF (IS_NEUMANN(I)) THEN
            RHS_NEW(I) = RHS_NEW(I) + NEUMANN(I)
         END IF
      END DO



      CALL COMPUTE_E_FIELD

      !  Advect the particles using the guessed new potential
      ALLOCATE(part_adv, SOURCE = particles)
      CALL TIMER_START(3)
      CALL ADVECT_CN(part_adv, .FALSE., .TRUE., Jmat)
      CALL TIMER_STOP(3)
      CALL DEPOSIT_CHARGE(part_adv)
      IF (JACOBIAN_TYPE == 4) CALL COMPUTE_MASS_MATRICES(part_adv)
      DEALLOCATE(part_adv)

      ! Compute the new potential from the charge distribution
      !CALL SOLVE_POISSON

      ! Compute the residual
      CALL VecGetOwnershipRange(f,Istart,Iend,ierr)

      CALL VecGetArrayF90(f,RESIDUAL,ierr_l)
      RESIDUAL = RHS(Istart:Iend-1) - RHS_NEW(Istart:Iend-1)
      !RESIDUAL = PHI_FIELD(Istart+1:Iend) + PHI_FIELD_OLD(Istart+1:Iend) - 2.*PHIBAR_FIELD(Istart+1:Iend)
      CALL VecRestoreArrayF90(f,RESIDUAL,ierr_l)


      CALL VecNorm(f,NORM_2,norm,ierr)
      IF (PROC_ID == 0) THEN
         !WRITE(*,*) ' PHI_FIELD = ', PHI_FIELD
         WRITE(*,*) '||RESIDUAL|| = ', norm !, ' with potential ', PHIBAR_FIELD
   
         WRITE(filename, "(A,A)") TRIM(ADJUSTL(RESIDUAL_SAVE_PATH)), "residuals" ! Compose filename   
         OPEN(66331, FILE=filename, POSITION='append', STATUS='unknown', ACTION='write')
         WRITE(66331,*) tID, norm
         CLOSE(66331)
      END IF

      DEALLOCATE(RHS_NEW)






      
      CALL MatAssemblyBegin(Jmat,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Jmat,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatGetOwnershipRange( Jmat, Istart, Iend, ierr)



      IF (JACOBIAN_TYPE == 5) THEN

         CALL MatCreate(PETSC_COMM_WORLD,Qmat,ierr)
         CALL MatSetSizes(Qmat,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
         CALL MatSetType(Qmat, MATMPIAIJ, ierr)
         CALL MatMPIAIJSetPreallocation(Qmat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER,ierr) ! DBDBDBDBDBDB Large preallocation!
         CALL MatSetFromOptions(Qmat,ierr)
         CALL MatSetUp(Qmat,ierr)

         CALL MatCreate(PETSC_COMM_WORLD,Rmat,ierr)
         CALL MatSetSizes(Rmat,PETSC_DECIDE,PETSC_DECIDE,NNODES,NNODES,ierr)
         CALL MatSetType(Rmat, MATMPIAIJ, ierr)
         CALL MatMPIAIJSetPreallocation(Rmat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER,ierr) ! DBDBDBDBDBDB Large preallocation!
         CALL MatSetFromOptions(Rmat,ierr)
         CALL MatSetUp(Rmat,ierr)

         ELECTRON_S_ID = SPECIES_NAME_TO_ID('e')
         ME = SPECIES(ELECTRON_S_ID)%MOLECULAR_MASS

         IF (DIMS == 2) THEN
            DO I = 1, NCELLS
               AREA = CELL_AREAS(I)
               FACTOR1 = 0.25*DT*DT*KB*CELL_TE(I)/ME
               FACTOR2 = -0.25*DT*DT*QE*QE*CELL_NE(I)/EPS0/ME
               ! We need to ADD to a sparse matrix entry.
               DO P = 1, 3
                  VP = U2D_GRID%CELL_NODES(P,I)
                  IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                     DO Q = 1, 3
                        VQ = U2D_GRID%CELL_NODES(Q,I)
                        IF (Q==P) THEN
                           MPQ = AREA/6.
                        ELSE
                           MPQ = AREA/12.
                        END IF
                        KPQ = AREA*(U2D_GRID%BASIS_COEFFS(1,P,I)*U2D_GRID%BASIS_COEFFS(1,Q,I) &
                                  + U2D_GRID%BASIS_COEFFS(2,P,I)*U2D_GRID%BASIS_COEFFS(2,Q,I))
                        
                        CALL MatSetValues(Qmat,one,VP-1,one,VQ-1,MPQ,ADD_VALUES,ierr)
                        CALL MatSetValues(Qmat,one,VP-1,one,VQ-1,KPQ*FACTOR1,ADD_VALUES,ierr)
                        CALL MatSetValues(Rmat,one,VP-1,one,VQ-1,KPQ*FACTOR2,ADD_VALUES,ierr)
                     END DO
                  END IF
               END DO
   
            END DO
         ELSE
            CALL ERROR_ABORT('Not implemented!')
         END IF

         CALL MatAssemblyBegin(Qmat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(Qmat,MAT_FINAL_ASSEMBLY,ierr)

         CALL MatAssemblyBegin(Rmat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(Rmat,MAT_FINAL_ASSEMBLY,ierr)

         CALL MatGetOrdering(Qmat, MATORDERINGNATURAL, rowperm, colperm, ierr)

         WRITE(*,*) 'Solved MatGetOrdering'

         CALL MatGetFactor(Qmat, MATSOLVERMUMPS, MAT_FACTOR_LU, Qmatfact,ierr)
         WRITE(*,*) 'Solved MatGetFactor'

         CALL MatLUFactorSymbolic(Qmatfact, Qmat, rowperm, colperm, info, ierr)
         WRITE(*,*) 'Solved MatLUFactorSymbolic'
         CALL MatLUFactorNumeric(Qmatfact, Qmat, info,ierr)
         WRITE(*,*) 'Solved MatLUFactorNumeric'
         CALL MatMatSolve(Qmatfact, Rmat, Jmat, ierr)

      END IF


      ! Accumulate Jacobian. All this is in principle not needed since we already have Amat.

      IF (DIMS == 2) THEN
         DO I = 1, NCELLS
            AREA = CELL_AREAS(I)
            
            ! We need to ADD to a sparse matrix entry.
            DO P = 1, 3
               VP = U2D_GRID%CELL_NODES(P,I)
               IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(VP-1)) THEN
                     DO Q = 1, 3
                        VQ = U2D_GRID%CELL_NODES(Q,I)
                        KPQ = AREA*(U2D_GRID%BASIS_COEFFS(1,P,I)*U2D_GRID%BASIS_COEFFS(1,Q,I) &
                                  + U2D_GRID%BASIS_COEFFS(2,P,I)*U2D_GRID%BASIS_COEFFS(2,Q,I))
                        
                        IF (AXI) THEN
                           V1 = U2D_GRID%CELL_NODES(1,I)
                           V2 = U2D_GRID%CELL_NODES(2,I)
                           V3 = U2D_GRID%CELL_NODES(3,I)  
                           KPQ = KPQ*(U2D_GRID%NODE_COORDS(2, V1) &
                                    + U2D_GRID%NODE_COORDS(2, V2) &
                                    + U2D_GRID%NODE_COORDS(2, V3))/3.
                        END IF

                        IF (JACOBIAN_TYPE == 3 .OR. JACOBIAN_TYPE == 4) THEN
                           KPQ = KPQ * (MASS_MATRIX(I)+1)
                        END IF

                        CALL MatSetValues(Jmat,one,VP-1,one,VQ-1,-2.*KPQ,ADD_VALUES,ierr)
                     END DO
                  END IF
               END IF
            END DO

         END DO
      ELSE IF (DIMS == 3) THEN
         DO I = 1, NCELLS
            VOLUME = CELL_VOLUMES(I)
            
            ! We need to ADD to a sparse matrix entry.
            DO P = 1, 4
               VP = U3D_GRID%CELL_NODES(P,I)
               IF (VP-1 >= Istart .AND. VP-1 < Iend) THEN
                  IF (.NOT. IS_DIRICHLET(VP-1)) THEN
                     DO Q = 1, 4
                        VQ = U3D_GRID%CELL_NODES(Q,I)
                        KPQ = VOLUME*(U3D_GRID%BASIS_COEFFS(1,P,I)*U3D_GRID%BASIS_COEFFS(1,Q,I) &
                                    + U3D_GRID%BASIS_COEFFS(2,P,I)*U3D_GRID%BASIS_COEFFS(2,Q,I) &
                                    + U3D_GRID%BASIS_COEFFS(3,P,I)*U3D_GRID%BASIS_COEFFS(3,Q,I))
                        
                        IF (JACOBIAN_TYPE == 3 .OR. JACOBIAN_TYPE == 4) THEN
                           KPQ = KPQ * (MASS_MATRIX(I)+1)
                        END IF

                        CALL MatSetValues(Jmat,one,VP-1,one,VQ-1,-2.*KPQ,ADD_VALUES,ierr)
                     END DO
                  END IF
               END IF
            END DO

         END DO
      END IF

      CALL MatAssemblyBegin(Jmat,MAT_FLUSH_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Jmat,MAT_FLUSH_ASSEMBLY,ierr)

      DO I = Istart, Iend-1
         IF (IS_DIRICHLET(I)) CALL MatSetValues(Jmat,one,I,one,I,-2.d0,INSERT_VALUES,ierr)
      END DO

      CALL MatAssemblyBegin(Jmat,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Jmat,MAT_FINAL_ASSEMBLY,ierr)

      ! CALL MatView(Jmat,PETSC_VIEWER_DRAW_WORLD,ierr)

   END SUBROUTINE FormFunctionAndJacobian




   SUBROUTINE COMPUTE_MASS_MATRICES(part_to_deposit)

      IMPLICIT NONE

      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: part_to_deposit
      INTEGER :: JP, IC

      REAL(KIND=8) :: CHARGE
      INTEGER :: SIZEC



      IF (GRID_TYPE == UNSTRUCTURED) THEN
         SIZEC = NCELLS
      ELSE
         CALL ERROR_ABORT('Not implemented.')
      END IF

      IF (.NOT. ALLOCATED(MASS_MATRIX)) ALLOCATE(MASS_MATRIX(SIZEC))
      MASS_MATRIX = 0.d0

      DO JP = 1, NP_PROC
         CHARGE = SPECIES(part_to_deposit(JP)%S_ID)%CHARGE
         IF (ABS(CHARGE) .LT. 1.d-6) CYCLE

         IF (GRID_TYPE == UNSTRUCTURED) THEN
            IC = part_to_deposit(JP)%IC

            IF (DIMS == 2) THEN
               
               MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*part_to_deposit(JP)%DTRIM/EPS0/CELL_AREAS(IC)/(ZMAX-ZMIN)*FNUM &
                                 * (QE*CHARGE)**2/SPECIES(part_to_deposit(JP)%S_ID)%MOLECULAR_MASS

            ELSE IF (DIMS == 3) THEN

               MASS_MATRIX(IC) = MASS_MATRIX(IC) + 0.25*DT*part_to_deposit(JP)%DTRIM/EPS0/CELL_VOLUMES(IC)*FNUM &
                                 * (QE*CHARGE)**2/SPECIES(part_to_deposit(JP)%S_ID)%MOLECULAR_MASS

            END IF
         ELSE

            CALL ERROR_ABORT('Not implemented.')

         END IF
      END DO

      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, MASS_MATRIX, SIZEC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ELSE
         CALL MPI_REDUCE(MASS_MATRIX,  MASS_MATRIX, SIZEC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF

      CALL MPI_BCAST(MASS_MATRIX, SIZEC, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   END SUBROUTINE COMPUTE_MASS_MATRICES




   SUBROUTINE COMPUTE_DENSITY_TEMPERATURE(part_to_deposit)

      IMPLICIT NONE

      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: part_to_deposit

      INTEGER, DIMENSION(:), ALLOCATABLE      :: TIMESTEP_NP
   
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VY
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VZ

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VX2
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VY2
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VZ2

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TTRX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TTRY
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TTRZ
      
      INTEGER :: JP, IC, ELECTRON_S_ID
      REAL(KIND=8) :: ME, NUMPART

      ELECTRON_S_ID = SPECIES_NAME_TO_ID('e')
      ME = SPECIES(ELECTRON_S_ID)%MOLECULAR_MASS

      IF (.NOT. GRID_TYPE == UNSTRUCTURED) CALL ERROR_ABORT('Not implemented.')


      IF (.NOT. ALLOCATED(CELL_NE)) ALLOCATE(CELL_NE(NCELLS))
      CELL_NE = 0
      IF (.NOT. ALLOCATED(CELL_TE)) ALLOCATE(CELL_TE(NCELLS))
      CELL_TE = 0

      ALLOCATE(TIMESTEP_NP(NCELLS))

      ALLOCATE(TIMESTEP_VX(NCELLS))
      ALLOCATE(TIMESTEP_VY(NCELLS))
      ALLOCATE(TIMESTEP_VZ(NCELLS))

      ALLOCATE(TIMESTEP_VX2(NCELLS))
      ALLOCATE(TIMESTEP_VY2(NCELLS))
      ALLOCATE(TIMESTEP_VZ2(NCELLS))

      ALLOCATE(TIMESTEP_TTRX(NCELLS))
      ALLOCATE(TIMESTEP_TTRY(NCELLS))
      ALLOCATE(TIMESTEP_TTRZ(NCELLS))

      TIMESTEP_NP = 0

      TIMESTEP_VX = 0
      TIMESTEP_VY = 0
      TIMESTEP_VZ = 0

      TIMESTEP_VX2 = 0
      TIMESTEP_VY2 = 0
      TIMESTEP_VZ2 = 0

      TIMESTEP_TTRX = 0
      TIMESTEP_TTRY = 0
      TIMESTEP_TTRZ = 0

      DO JP = 1, NP_PROC
         IF (part_to_deposit(JP)%S_ID .NE. ELECTRON_S_ID) CYCLE

         IC = part_to_deposit(JP)%IC

         TIMESTEP_NP(IC) = TIMESTEP_NP(IC) + 1
         TIMESTEP_VX(IC) = TIMESTEP_VX(IC) + part_to_deposit(JP)%VX
         TIMESTEP_VY(IC) = TIMESTEP_VY(IC) + part_to_deposit(JP)%VY
         TIMESTEP_VZ(IC) = TIMESTEP_VZ(IC) + part_to_deposit(JP)%VZ
         TIMESTEP_VX2(IC) = TIMESTEP_VX2(IC) + part_to_deposit(JP)%VX*part_to_deposit(JP)%VX
         TIMESTEP_VY2(IC) = TIMESTEP_VY2(IC) + part_to_deposit(JP)%VY*part_to_deposit(JP)%VY
         TIMESTEP_VZ2(IC) = TIMESTEP_VZ2(IC) + part_to_deposit(JP)%VZ*part_to_deposit(JP)%VZ

      END DO

      ! Collect data from all the processes
      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_NP,  NCELLS, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VX,  NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VY,  NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VZ,  NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VX2, NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VY2, NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VZ2, NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

         DO IC = 1, NCELLS
            IF (TIMESTEP_NP(IC) .GT. 1) THEN
               NUMPART = DBLE(TIMESTEP_NP(IC))
               
               TIMESTEP_VX(IC) = TIMESTEP_VX(IC) / NUMPART
               TIMESTEP_VY(IC) = TIMESTEP_VY(IC) / NUMPART
               TIMESTEP_VZ(IC) = TIMESTEP_VZ(IC) / NUMPART
               TIMESTEP_TTRX(IC) = ME / KB &
                  * (TIMESTEP_VX2(IC)/NUMPART - TIMESTEP_VX(IC)*TIMESTEP_VX(IC))
               TIMESTEP_TTRY(IC) = ME / KB &
                  * (TIMESTEP_VY2(IC)/NUMPART - TIMESTEP_VY(IC)*TIMESTEP_VY(IC))
               TIMESTEP_TTRZ(IC) = ME / KB &
                  * (TIMESTEP_VZ2(IC)/NUMPART - TIMESTEP_VZ(IC)*TIMESTEP_VZ(IC))
            END IF
         END DO
   
         CELL_NE = TIMESTEP_NP*FNUM / CELL_VOLUMES
         CELL_TE = (TIMESTEP_TTRX + TIMESTEP_TTRY + TIMESTEP_TTRZ) / 3.

      ELSE
         CALL MPI_REDUCE(TIMESTEP_NP,   TIMESTEP_NP,  NCELLS, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VX,   TIMESTEP_VX,  NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VY,   TIMESTEP_VY,  NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VZ,   TIMESTEP_VZ,  NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VX2,  TIMESTEP_VX2, NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VY2,  TIMESTEP_VY2, NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VZ2,  TIMESTEP_VZ2, NCELLS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF

      CALL MPI_BCAST(CELL_NE, NCELLS, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(CELL_TE, NCELLS, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      DEALLOCATE(TIMESTEP_NP)

      DEALLOCATE(TIMESTEP_VX)
      DEALLOCATE(TIMESTEP_VY)
      DEALLOCATE(TIMESTEP_VZ)

      DEALLOCATE(TIMESTEP_VX2)
      DEALLOCATE(TIMESTEP_VY2)
      DEALLOCATE(TIMESTEP_VZ2)

      DEALLOCATE(TIMESTEP_TTRX)
      DEALLOCATE(TIMESTEP_TTRY)
      DEALLOCATE(TIMESTEP_TTRZ)

   END SUBROUTINE COMPUTE_DENSITY_TEMPERATURE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ADVECT_CN -> Advects particles in the domain using a Crank-Nicholson scheme for fully implicit !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ADVECT_CN(part_adv, FINAL, COMPUTE_JACOBIAN, jac)

      ! This subroutine advects the particles belonging to the processor.
      ! Particles are advected for a timestep DT. A temporary variable DTRIM is
      ! used: if a collision with a wall happens, the wall_collision subroutine 
      ! is called, the free flight time before collision is computed, and then 
      ! advection is performed for the remaining time. ! TO BE IMPLEMENTED!
      !
      ! After the advection, periodic BCs are applied, so that the periodic 
      ! coordinates stay into the domain boundaries.
      !
      ! Finally, if particles are out of the domain (and it happens only if 
      ! the periodicity wasn't active for some direction), remove them. 

      IMPLICIT NONE

      Mat jac
      Mat dxde, dxdexmat, dxdeymat, dydexmat, dydeymat
      !PetscInt row
      PetscInt ncols
      PetscInt cols(2000)
      PetscScalar dxdexvals(2000), dxdeyvals(2000), dydexvals(2000), dydeyvals(2000), vals(2000)
      PetscInt first_row, last_row

      !PetscViewer  viewer

      LOGICAL, INTENT(IN) :: FINAL, COMPUTE_JACOBIAN
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: part_adv

      INTEGER      :: IP, IC, I, J, OLD_IC, JJ, SIZE
      INTEGER      :: BOUNDCOLL, FACE_PG, NEIGHBOR
      REAL(KIND=8) :: DTCOLL, TOTDTCOLL, rfp, QOM
      REAL(KIND=8) :: COEFA, COEFB, COEFC, DELTA, ALPHA, BETA, GAMMA, DTADV
      !REAL(KIND=8), DIMENSION(2) :: TEST
      REAL(KIND=8), DIMENSION(4) :: NORMX, NORMY, XW, YW, ZW, BOUNDPOS
      ! REAL(KIND=8) :: XCOLL, YCOLL, ZCOLL
      LOGICAL, DIMENSION(:), ALLOCATABLE :: REMOVE_PART
      REAL(KIND=8), DIMENSION(3) :: E, B
      REAL(KIND=8) ::V_PERP, VDUMMY, EROT, EVIB, VDOTN, WALL_TEMP
      INTEGER :: S_ID
      LOGICAL :: HASCOLLIDED
      REAL(KIND=8) :: EDGE_X1, EDGE_Y1, EDGE_X2, EDGE_Y2
      INTEGER, DIMENSION(:), ALLOCATABLE :: LOCAL_BOUNDARY_COLL_COUNT, LOCAL_WALL_COLL_COUNT
      REAL(KIND=8) :: WEIGHT_RATIO
      TYPE(PARTICLE_DATA_STRUCTURE) :: NEWparticle
      INTEGER, DIMENSION(3) :: NEXTVERT
      INTEGER, DIMENSION(8) :: EDGEINDEX
      REAL(KIND=8), DIMENSION(8) :: COLLTIMES
      INTEGER, DIMENSION(2000) :: EBARIDX
      REAL(KIND=8), DIMENSION(2000) :: DVXDEX, DVYDEX, DVXDEY, DVYDEY, DXDEX, DYDEX, DXDEY, DYDEY, DTDEX, DTDEY, &
      DTDEXPREC, DTDEYPREC, DALPHADEX, DALPHADEY, DBETADEX, DBETADEY, DGAMMADEX, DGAMMADEY
      REAL(KIND=8) :: NEWVX, NEWVY
      INTEGER :: LOC, NUMSTEPS, MAXNUMSTEPS
      REAL(KIND=8) :: VALXX, VALXY, VALYX, VALYY
      REAL(KIND=8) :: DPSI1DX, DPSI1DY, DPSI2DX, DPSI2DY, DPSI3DX, DPSI3DY, DPSJ1DX, DPSJ1DY, DPSJ2DX, DPSJ2DY, DPSJ3DX, DPSJ3DY
      INTEGER :: V1I, V2I, V3I, V1J, V2J, V3J
      REAL(KIND=8) :: X_TEMP, Y_TEMP, Z_TEMP, VX_TEMP, VY_TEMP, VZ_TEMP
      REAL(KIND=8), DIMENSION(3) :: FACE_NORMAL, FACE_TANG1, FACE_TANG2
      REAL(KIND=8) :: V_TANG1, V_TANG2
      REAL(KIND=8), DIMENSION(3) :: TANG1, TANG2
      REAL(KIND=8) :: VDOTTANG1, VRM, RN, R1, R2, THETA1, THETA2, DOT_NORM, VTANGENT
      INTEGER :: NI, NJ, VNI, VNJ

      INTEGER :: NCROSSINGS

      INTEGER :: NCOLSXX, NCOLSXY, NCOLSYX, NCOLSYY
      INTEGER, DIMENSION(:), ALLOCATABLE :: COLSXX, COLSXY, COLSYX, COLSYY
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: DXDEXV, DXDEYV, DYDEXV, DYDEYV

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: DXDEBAR, DVDEBAR

      LOGICAL :: FLUIDBOUNDARY
      INTEGER :: NEIGHBORPG
      REAL(KIND=8) :: CHARGE, K, PSIP, RHO_Q
      INTEGER :: VP

      !REAL(KIND=8) :: CHECKVALUE

      ! IF (tID == 4 .AND. FINAL) THEN
      !    OPEN(66332, FILE='crossings', POSITION='append', STATUS='unknown', ACTION='write')
      ! END IF

      IF (JACOBIAN_TYPE .LT. 0 .OR. JACOBIAN_TYPE .GT. 7) CALL ERROR_ABORT('Jacobian type not supported.')

      SIZE = NCELLS

      MAXNUMSTEPS = 0
      IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 1) THEN
         CALL MatCreate(PETSC_COMM_WORLD,dxde,ierr)
         CALL MatSetSizes(dxde,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
         CALL MatSetType(dxde, MATMPIAIJ, ierr)
         !CALL MatSetOption(dxde,MAT_SPD,PETSC_TRUE,ierr)
         CALL MatMPIAIJSetPreallocation(dxde,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
         CALL MatSetFromOptions(dxde, ierr)
         CALL MatSetUp(dxde,ierr)

         CALL MatCreate(PETSC_COMM_WORLD,dxdexmat,ierr)
         CALL MatSetSizes(dxdexmat,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
         CALL MatSetType(dxdexmat, MATMPIAIJ, ierr)
         !CALL MatSetOption(dxdexmat,MAT_SPD,PETSC_TRUE,ierr)
         CALL MatMPIAIJSetPreallocation(dxdexmat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
         CALL MatSetFromOptions(dxdexmat, ierr)
         CALL MatSetUp(dxdexmat,ierr)

         CALL MatCreate(PETSC_COMM_WORLD,dxdeymat,ierr)
         CALL MatSetSizes(dxdeymat,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
         CALL MatSetType(dxdeymat, MATMPIAIJ, ierr)
         !CALL MatSetOption(dxdeymat,MAT_SPD,PETSC_TRUE,ierr)
         CALL MatMPIAIJSetPreallocation(dxdeymat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
         CALL MatSetFromOptions(dxdeymat, ierr)
         CALL MatSetUp(dxdeymat,ierr)

         CALL MatCreate(PETSC_COMM_WORLD,dydexmat,ierr)
         CALL MatSetSizes(dydexmat,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
         CALL MatSetType(dydexmat, MATMPIAIJ, ierr)
         !CALL MatSetOption(dydexmat,MAT_SPD,PETSC_TRUE,ierr)
         CALL MatMPIAIJSetPreallocation(dydexmat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
         CALL MatSetFromOptions(dydexmat, ierr)
         CALL MatSetUp(dydexmat,ierr)

         CALL MatCreate(PETSC_COMM_WORLD,dydeymat,ierr)
         CALL MatSetSizes(dydeymat,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
         CALL MatSetType(dydeymat, MATMPIAIJ, ierr)
         !CALL MatSetOption(dydeymat,MAT_SPD,PETSC_TRUE,ierr)
         CALL MatMPIAIJSetPreallocation(dydeymat,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
         CALL MatSetFromOptions(dydeymat, ierr)
         CALL MatSetUp(dydeymat,ierr)

      ELSE IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 2) THEN
         
         CALL MatCreate(PETSC_COMM_WORLD,dxde,ierr)
         CALL MatSetSizes(dxde,PETSC_DECIDE, PETSC_DECIDE, SIZE, SIZE, ierr)
         CALL MatSetType(dxde, MATMPIAIJ, ierr)
         !CALL MatSetOption(dxde,MAT_SPD,PETSC_TRUE,ierr)
         CALL MatMPIAIJSetPreallocation(dxde,2000,PETSC_NULL_INTEGER,2000,PETSC_NULL_INTEGER, ierr) !! DBDBDBDBDBDBDBDBDDBDB Large preallocation!
         CALL MatSetFromOptions(dxde, ierr)
         CALL MatSetUp(dxde,ierr)

         IF (ALLOCATED(DVDEBAR)) DEALLOCATE(DVDEBAR)
         ALLOCATE(DVDEBAR(SIZE))
         
         IF (ALLOCATED(DXDEBAR)) DEALLOCATE(DXDEBAR)
         ALLOCATE(DXDEBAR(SIZE))
      
      END IF



      !E = [0.d0, 0.d0, 0.d0]
      B = [0.d0, 0.d0, 0.d0]

      NEXTVERT = [2,3,1]
      EDGEINDEX = [1,1,2,2,3,3,4,4]

      NORMX = [ 1., -1., 0., 0. ]
      NORMY = [ 0., 0., 1., -1. ]

      XW = [ XMIN, XMAX, XMIN, XMAX ]
      YW = [ YMAX, YMIN, YMIN, YMAX ]
      ZW = [ 0., 0., 0., 0. ]

      BOUNDPOS = [ XMIN, XMAX, YMIN, YMAX ]

      ALLOCATE(REMOVE_PART(NP_PROC))
      
      ALLOCATE(LOCAL_BOUNDARY_COLL_COUNT(4*N_SPECIES))
      LOCAL_BOUNDARY_COLL_COUNT = 0
      ALLOCATE(LOCAL_WALL_COLL_COUNT(N_WALLS*N_SPECIES))
      LOCAL_WALL_COLL_COUNT = 0

      DO IP = 1, NP_PROC
         NCROSSINGS = 0
         REMOVE_PART(IP) = .FALSE.
         IC = part_adv(IP)%IC

         IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 1) THEN
            DVXDEX = 0.d0
            DVYDEX = 0.d0
            DVXDEY = 0.d0
            DVYDEY = 0.d0

            DXDEX = 0.d0
            DYDEX = 0.d0
            DXDEY = 0.d0
            DYDEY = 0.d0

            DTDEX = 0.d0
            DTDEY = 0.d0

            DTDEXPREC = 0.d0
            DTDEYPREC = 0.d0

            EBARIDX = -1
            NUMSTEPS = 0
         ELSE IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 2) THEN
            DXDEBAR = 0.d0
            DVDEBAR = 0.d0
            EBARIDX = -1
            NUMSTEPS = 0
         END IF


         ! ! Update velocity

         ! V_NEW(1) = part_adv(IP)%VX
         ! V_NEW(2) = part_adv(IP)%VY
         ! V_NEW(3) = part_adv(IP)%VZ

         ! IF (PIC_TYPE .NE. NONE) THEN
         !    !PHIBAR_FIELD = 0.d0
         !    !EBAR_FIELD = 0.d0
         !    CALL APPLY_E_FIELD(IP, E)
         !    !CALL APPLY_RF_E_FIELD(IP, E)
            


         !    V_OLD(1) = part_adv(IP)%VX
         !    V_OLD(2) = part_adv(IP)%VY
         !    V_OLD(3) = part_adv(IP)%VZ

         !    CALL UPDATE_VELOCITY_BORIS(part_adv(IP)%DTRIM, V_OLD, V_NEW, &
         !    SPECIES(part_adv(IP)%S_ID)%CHARGE, SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS, &
         !    E, B)

         !    ! Assign v^n to the particle, for simplicity
         !    !part_adv(IP)%VX = 0.5*(V_OLD(1) + V_NEW(1))
         !    !part_adv(IP)%VY = 0.5*(V_OLD(2) + V_NEW(2))
         !    !part_adv(IP)%VZ = 0.5*(V_OLD(3) + V_NEW(3))
         !    part_adv(IP)%VX = V_NEW(1)
         !    part_adv(IP)%VY = V_NEW(2)
         !    part_adv(IP)%VZ = V_NEW(3)
         ! END IF



         ! ! Forced electric field
         ! E = [3000.d0, 0.d0, 0.d0]
            


         ! V_OLD(1) = part_adv(IP)%VX
         ! V_OLD(2) = part_adv(IP)%VY
         ! V_OLD(3) = part_adv(IP)%VZ

         ! CALL UPDATE_VELOCITY_BORIS(DT, V_OLD, V_NEW, &
         ! SPECIES(part_adv(IP)%S_ID)%CHARGE, SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS, &
         ! E, B)

         ! ! Assign v^n to the particle, for simplicity
         ! !part_adv(IP)%VX = 0.5*(V_OLD(1) + V_NEW(1))
         ! !part_adv(IP)%VY = 0.5*(V_OLD(2) + V_NEW(2))
         ! !part_adv(IP)%VZ = 0.5*(V_OLD(3) + V_NEW(3))
         ! part_adv(IP)%VX = V_NEW(1)
         ! part_adv(IP)%VY = V_NEW(2)
         ! part_adv(IP)%VZ = V_NEW(3)

         ! ! End forced electric field.



         HASCOLLIDED = .FALSE.
         TOTDTCOLL = 0.
         DO WHILE (part_adv(IP)%DTRIM .GT. 0.) ! Repeat the procedure until step is done

            DTCOLL = part_adv(IP)%DTRIM ! Looking for collisions within the remaining time
            ! ______ ADVECTION ______
            IF (PIC_TYPE .NE. NONE) THEN
               E = EBAR_FIELD(:, 1, IC) ! CALL APPLY_E_FIELD(IP, E)
            ELSE
               E = 0.d0
            END IF
            !WRITE(*,*) 'Field on particle ', IP, ' on proc ', PROC_ID, ' in cell ', IC, ' is ', E

            IF (GRID_TYPE == UNSTRUCTURED) THEN
               !WRITE(*,*) 'Moving particle ', IP, ' for ', DTCOLL, ' s. Position: ', part_adv(IP)%X, part_adv(IP)%Y,&
               !' velocity: ', part_adv(IP)%VX, part_adv(IP)%VY
               ! For unstructured, we only need to check the boundaries of the cell.

               ! The new C-N procedure with uniform E field.
               IF (DIMS == 2) THEN
                  COLLTIMES = -1.d0
                  DO I = 1, 3
                     J = NEXTVERT(I)
                     EDGE_X1 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(I,IC))
                     EDGE_Y1 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(I,IC))
                     EDGE_X2 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(J,IC))
                     EDGE_Y2 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(J,IC))
                     !IF (EDGE_X2 == EDGE_X1) THEN
                     !   COEFB = 0; COEFA = 1; COEFC = - EDGE_X1
                     !ELSE IF (EDGE_Y2 == EDGE_Y1) THEN
                     !   COEFA = 0; COEFB = 1; COEFC = - EDGE_Y1
                     !ELSE
                     !   COEFA = (EDGE_Y2-EDGE_Y1)
                     !   COEFB = (EDGE_X1-EDGE_X2)
                     !   COEFC = (EDGE_X2 - EDGE_X1) * EDGE_Y1 - (EDGE_Y2 - EDGE_Y1) * EDGE_X1
                     !   !COEFA = 1
                     !   !COEFB = (EDGE_X1-EDGE_X2)/(EDGE_Y2-EDGE_Y1)
                     !   !COEFC = - COEFB * EDGE_Y1 - EDGE_X1
                     !END IF

                     COEFA = (EDGE_Y2-EDGE_Y1)
                     COEFB = (EDGE_X1-EDGE_X2)
                     !COEFC = (EDGE_X2 - EDGE_X1) * EDGE_Y1 - (EDGE_Y2 - EDGE_Y1) * EDGE_X1
                     COEFC = EDGE_Y1*EDGE_X2 - EDGE_X1*EDGE_Y2

                     IF (PIC_TYPE .NE. NONE) THEN
                        ALPHA = 0.5*SPECIES(part_adv(IP)%S_ID)%CHARGE*QE/SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS* &
                              (COEFA*E(1) + COEFB*E(2))
                     ELSE
                        ALPHA = 0.d0
                     END IF

                     BETA = COEFA*part_adv(IP)%VX + COEFB*part_adv(IP)%VY
                     GAMMA = COEFA*part_adv(IP)%X + COEFB*part_adv(IP)%Y + COEFC

                     IF (ALPHA == 0.d0) THEN
                        COLLTIMES(2*(I-1) + 1) = - GAMMA/BETA
                     ELSE
                        DELTA = BETA*BETA - 4.0*ALPHA*GAMMA
                        IF (DELTA >= 0.d0) THEN
                           COLLTIMES(2*(I-1) + 1) = 0.5*(-BETA - SQRT(DELTA))/ALPHA
                           COLLTIMES(2*I) =         0.5*(-BETA + SQRT(DELTA))/ALPHA
                        END IF
                     END IF
                  END DO

                  ! Find which collision happens first.
                  DTCOLL = part_adv(IP)%DTRIM
                  BOUNDCOLL = -1
                  DO I = 1, 6
                     IF (COLLTIMES(I) > 0 .AND. COLLTIMES(I) < DTCOLL) THEN


                        X_TEMP = part_adv(IP)%X
                        Y_TEMP = part_adv(IP)%Y
                        Z_TEMP = part_adv(IP)%Z

                        VX_TEMP = part_adv(IP)%VX
                        VY_TEMP = part_adv(IP)%VY
                        VZ_TEMP = part_adv(IP)%VZ

                        CALL MOVE_PARTICLE_CN(part_adv, IP, E, COLLTIMES(I))
                        J = EDGEINDEX(I)


                        ! ! A small check that we actually found an intersection.

                        ! JJ = NEXTVERT(J)
                        ! EDGE_X1 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(J,IC))
                        ! EDGE_Y1 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(J,IC))
                        ! EDGE_X2 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(JJ,IC))
                        ! EDGE_Y2 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(JJ,IC))

                        ! COEFA = (EDGE_Y2-EDGE_Y1)
                        ! COEFB = (EDGE_X1-EDGE_X2)
                        ! COEFC = (EDGE_X2 - EDGE_X1) * EDGE_Y1 - (EDGE_Y2 - EDGE_Y1) * EDGE_X1
                        ! WRITE(*,*) (COEFA*part_adv(IP)%X + COEFB*part_adv(IP)%Y + COEFC)/(ABS(COEFA)+ABS(COEFB)+ABS(COEFC))
                        ! ! End of the small check.
                        
                        IF ((part_adv(IP)%VX*U2D_GRID%EDGE_NORMAL(1,J,IC) + part_adv(IP)%VY*U2D_GRID%EDGE_NORMAL(2,J,IC)) > 0) THEN
                           ! Velocity pre-rotation points outwards.
                           DTCOLL = COLLTIMES(I)
                           BOUNDCOLL = J
                           IF (AXI) THEN
                              CALL AXI_ROTATE_VELOCITY(part_adv, IP)
                              IF ((part_adv(IP)%VX*U2D_GRID%EDGE_NORMAL(1,J,IC) &
                                 + part_adv(IP)%VY*U2D_GRID%EDGE_NORMAL(2,J,IC)) < 0) THEN
                                 ! After rotation the velocity points back into the cell.
                                 ! We should move without applying the boundary conditions, and start again from there.
                                 BOUNDCOLL = 0
                                 !WRITE(*,*) 'This happened to particle with id ', part_adv(IP)%ID, ' and V dot N was ', &
                                 !(part_adv(IP)%VX*U2D_GRID%EDGE_NORMAL(1,J,IC) + part_adv(IP)%VY*U2D_GRID%EDGE_NORMAL(2,J,IC))
                                 ! We should move without applying the boundary conditions, and start again from there.
                              END IF
                           END IF
                        ELSE
                           IF (AXI) THEN
                              CALL AXI_ROTATE_VELOCITY(part_adv, IP)
                              IF ((part_adv(IP)%VX*U2D_GRID%EDGE_NORMAL(1,J,IC) &
                                 + part_adv(IP)%VY*U2D_GRID%EDGE_NORMAL(2,J,IC)) > 0) THEN
                                 ! After rotation the velocity points out of the cell.
                                 DTCOLL = COLLTIMES(I)
                                 BOUNDCOLL = J
                                 !WRITE(*,*) 'The opposite happened to particle with id ', part_adv(IP)%ID, ' and V dot N was ', &
                                 !(part_adv(IP)%VX*U2D_GRID%EDGE_NORMAL(1,J,IC) + part_adv(IP)%VY*U2D_GRID%EDGE_NORMAL(2,J,IC))
                                 ! We should move without applying the boundary conditions, and start again from there.
                              END IF
                           END IF
                        END IF

                        part_adv(IP)%X = X_TEMP
                        part_adv(IP)%Y = Y_TEMP
                        part_adv(IP)%Z = Z_TEMP

                        part_adv(IP)%VX = VX_TEMP
                        part_adv(IP)%VY = VY_TEMP
                        part_adv(IP)%VZ = VZ_TEMP

                     END IF
                  END DO
               ELSE IF (DIMS == 3) THEN
                  COLLTIMES = -1.d0
                  DO I = 1, 4
                     IF (PIC_TYPE .NE. NONE) THEN
                        ALPHA = 0.5*SPECIES(part_adv(IP)%S_ID)%CHARGE*QE/SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS* &
                              (U3D_GRID%CELL_FACES_COEFFS(1,I,IC)*E(1) &
                             + U3D_GRID%CELL_FACES_COEFFS(2,I,IC)*E(2) &
                             + U3D_GRID%CELL_FACES_COEFFS(3,I,IC)*E(3))
                     ELSE
                        ALPHA = 0.d0
                     END IF

                     BETA = U3D_GRID%CELL_FACES_COEFFS(1,I,IC)*part_adv(IP)%VX &
                          + U3D_GRID%CELL_FACES_COEFFS(2,I,IC)*part_adv(IP)%VY &
                          + U3D_GRID%CELL_FACES_COEFFS(3,I,IC)*part_adv(IP)%VZ

                     GAMMA = U3D_GRID%CELL_FACES_COEFFS(1,I,IC)*part_adv(IP)%X &
                           + U3D_GRID%CELL_FACES_COEFFS(2,I,IC)*part_adv(IP)%Y &
                           + U3D_GRID%CELL_FACES_COEFFS(3,I,IC)*part_adv(IP)%Z &
                           + U3D_GRID%CELL_FACES_COEFFS(4,I,IC)

                     IF (ALPHA == 0.d0) THEN
                        COLLTIMES(2*(I-1) + 1) = - GAMMA/BETA
                     ELSE
                        DELTA = BETA*BETA - 4.0*ALPHA*GAMMA
                        IF (DELTA >= 0.d0) THEN
                           COLLTIMES(2*(I-1) + 1) = 0.5*(-BETA - SQRT(DELTA))/ALPHA
                           COLLTIMES(2*I) =         0.5*(-BETA + SQRT(DELTA))/ALPHA
                        END IF
                     END IF
                  END DO

                  ! Find which collision happens first.
                  DTCOLL = part_adv(IP)%DTRIM
                  BOUNDCOLL = -1
                  DO I = 1, 8
                     IF (COLLTIMES(I) > 0 .AND. COLLTIMES(I) <= DTCOLL) THEN


                        X_TEMP = part_adv(IP)%X
                        Y_TEMP = part_adv(IP)%Y
                        Z_TEMP = part_adv(IP)%Z

                        VX_TEMP = part_adv(IP)%VX
                        VY_TEMP = part_adv(IP)%VY
                        VZ_TEMP = part_adv(IP)%VZ

                        CALL MOVE_PARTICLE_CN(part_adv, IP, E, COLLTIMES(I))
                        J = EDGEINDEX(I)

                        ! ! A small check that we actually found an intersection.
                        ! CHECKVALUE = U3D_GRID%CELL_FACES_COEFFS(1,J,IC)*part_adv(IP)%X &
                        ! + U3D_GRID%CELL_FACES_COEFFS(2,J,IC)*part_adv(IP)%Y &
                        ! + U3D_GRID%CELL_FACES_COEFFS(3,J,IC)*part_adv(IP)%Z &
                        ! + U3D_GRID%CELL_FACES_COEFFS(4,J,IC)
                        ! IF (ABS(CHECKVALUE) > 1.0d-12) WRITE(*,*) 'Checkvalue too large! = ', CHECKVALUE
                        ! ! End of the small check.
                        
                        IF ((part_adv(IP)%VX*U3D_GRID%FACE_NORMAL(1,J,IC) &
                           + part_adv(IP)%VY*U3D_GRID%FACE_NORMAL(2,J,IC) &
                           + part_adv(IP)%VZ*U3D_GRID%FACE_NORMAL(3,J,IC)) > 0) THEN
                           DTCOLL = COLLTIMES(I)
                           BOUNDCOLL = J
                        END IF

                        part_adv(IP)%X = X_TEMP
                        part_adv(IP)%Y = Y_TEMP
                        part_adv(IP)%Z = Z_TEMP

                        part_adv(IP)%VX = VX_TEMP
                        part_adv(IP)%VY = VY_TEMP
                        part_adv(IP)%VZ = VZ_TEMP

                     END IF
                  END DO
               END IF


               IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 1) THEN
                  ! Accumulate exact Jacobian of the motion
                  ! The particle hit the boundary of a cell.
                  ! It may cross to the neighboring cell, or if this is a domain boundary be reflected or absorbed.
                  ! The total time of the motion is less than dt (or the particle's initial dtrim) and depends on the path taken.

                  LOC = NUMSTEPS + 1
                  DO I = 1, NUMSTEPS
                     IF (EBARIDX(I) == IC-1) THEN
                        LOC = I
                        EXIT
                     END IF
                  END DO
                  IF (LOC == NUMSTEPS + 1) THEN
                     NUMSTEPS = LOC
                     EBARIDX(LOC) = IC-1
                  END IF

                  IF (NUMSTEPS > MAXNUMSTEPS) MAXNUMSTEPS = NUMSTEPS

                  QOM = SPECIES(part_adv(IP)%S_ID)%CHARGE*QE/SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS

                  IF (BOUNDCOLL > 0) THEN
                     DTADV = DTCOLL

                     I = BOUNDCOLL
                     J = NEXTVERT(I)
                     EDGE_X1 = U2D_GRID%NODE_COORDS(1,U2D_GRID%CELL_NODES(I,IC))
                     EDGE_Y1 = U2D_GRID%NODE_COORDS(2,U2D_GRID%CELL_NODES(I,IC))
                     EDGE_X2 = U2D_GRID%NODE_COORDS(1,U2D_GRID%CELL_NODES(J,IC))
                     EDGE_Y2 = U2D_GRID%NODE_COORDS(2,U2D_GRID%CELL_NODES(J,IC))

                     COEFA = (EDGE_Y2-EDGE_Y1)
                     COEFB = (EDGE_X1-EDGE_X2)
                     COEFC = EDGE_Y1*EDGE_X2 - EDGE_X1*EDGE_Y2

                     ALPHA = 0.5*QOM*(COEFA*E(1) + COEFB*E(2))
                     BETA = COEFA*part_adv(IP)%VX + COEFB*part_adv(IP)%VY
                     GAMMA = COEFA*part_adv(IP)%X + COEFB*part_adv(IP)%Y + COEFC
                     DALPHADEX = 0.d0
                     DALPHADEY = 0.d0
                     DALPHADEX(LOC) = 0.5*QOM*COEFA
                     DALPHADEY(LOC) = 0.5*QOM*COEFB
                     DBETADEX = COEFA*DVXDEX+COEFB*DVYDEX
                     DBETADEY = COEFA*DVXDEY+COEFB*DVYDEY
                     DGAMMADEX = COEFA*DXDEX+COEFB*DYDEX
                     DGAMMADEY = COEFA*DXDEY+COEFB*DYDEY
                     DTDEX = -(DALPHADEX*DTADV*DTADV + DBETADEX*DTADV + DGAMMADEX)/(2.*ALPHA*DTADV + BETA)
                     DTDEY = -(DALPHADEY*DTADV*DTADV + DBETADEY*DTADV + DGAMMADEY)/(2.*ALPHA*DTADV + BETA)
                  ELSE
                     DTADV = part_adv(IP)%DTRIM
                  END IF

                  DXDEX = DXDEX + DVXDEX*DTADV
                  DXDEY = DXDEY + DVXDEY*DTADV
                  DYDEX = DYDEX + DVYDEX*DTADV
                  DYDEY = DYDEY + DVYDEY*DTADV

                  DXDEX(LOC) = DXDEX(LOC) + 0.5*QOM*DTADV*DTADV
                  DYDEY(LOC) = DYDEY(LOC) + 0.5*QOM*DTADV*DTADV

                  NEWVX = part_adv(IP)%VX + QOM*E(1)*DTADV
                  NEWVY = part_adv(IP)%VY + QOM*E(2)*DTADV
                  IF (BOUNDCOLL > 0) THEN
                     DXDEX = DXDEX + NEWVX*DTDEX
                     DXDEY = DXDEY + NEWVX*DTDEY
                     DYDEX = DYDEX + NEWVY*DTDEX
                     DYDEY = DYDEY + NEWVY*DTDEY
                  ELSE
                     DXDEX = DXDEX - NEWVX*DTDEXPREC
                     DXDEY = DXDEY - NEWVX*DTDEYPREC
                     DYDEX = DYDEX - NEWVY*DTDEXPREC
                     DYDEY = DYDEY - NEWVY*DTDEYPREC
                  END IF

                  DTDEXPREC = DTDEXPREC + DTDEX
                  DTDEYPREC = DTDEYPREC + DTDEY

                  DVXDEX = DVXDEX + QOM*E(1)*DTDEX
                  DVYDEX = DVYDEX + QOM*E(2)*DTDEX
                  DVXDEY = DVXDEY + QOM*E(1)*DTDEY
                  DVYDEY = DVYDEY + QOM*E(2)*DTDEY

                  DVXDEX(LOC) = DVXDEX(LOC) + QOM*DTADV
                  DVYDEY(LOC) = DVYDEY(LOC) + QOM*DTADV

               ELSE IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 2) THEN
                  ! Accumulate approximate Jacobian of the motion
                  LOC = NUMSTEPS + 1
                  DO I = 1, NUMSTEPS
                     IF (EBARIDX(I) == IC-1) THEN
                        LOC = I
                        EXIT
                     END IF
                  END DO
                  IF (LOC == NUMSTEPS + 1) THEN
                     NUMSTEPS = LOC
                     EBARIDX(LOC) = IC-1
                  END IF
                  
                  IF (NUMSTEPS > MAXNUMSTEPS) MAXNUMSTEPS = NUMSTEPS

                  IF (BOUNDCOLL > 0) THEN
                     DTADV = DTCOLL
                  ELSE
                     DTADV = part_adv(IP)%DTRIM
                  END IF

                  DXDEBAR = DXDEBAR + DVDEBAR*DTCOLL
                  DXDEBAR(LOC) = DXDEBAR(LOC) &
                  + 0.5*SPECIES(part_adv(IP)%S_ID)%CHARGE*QE/SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS*DTADV*DTADV
                  DVDEBAR(LOC) = DVDEBAR(LOC) &
                  + SPECIES(part_adv(IP)%S_ID)%CHARGE*QE/SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS*DTADV
               END IF


               ! Do the advection
               IF (BOUNDCOLL > 0) THEN
                  CALL MOVE_PARTICLE_CN(part_adv, IP, E, DTCOLL)
                  IF (AXI .AND. DIMS == 2) CALL AXI_ROTATE_VELOCITY(part_adv, IP)
                  part_adv(IP)%DTRIM = part_adv(IP)%DTRIM - DTCOLL

                  FLUIDBOUNDARY = .FALSE.
                  IF (DIMS == 2) THEN
                     NEIGHBOR = U2D_GRID%CELL_NEIGHBORS(BOUNDCOLL, IC)
                     IF (NEIGHBOR == -1) THEN
                        FLUIDBOUNDARY = .TRUE.
                     ELSE
                        NEIGHBORPG = U2D_GRID%CELL_PG(NEIGHBOR)
                        IF (NEIGHBORPG .NE. -1) THEN
                           IF (GRID_BC(NEIGHBORPG)%VOLUME_BC == SOLID) FLUIDBOUNDARY = .TRUE.
                        END IF
                     END IF

                     FACE_PG = U2D_GRID%CELL_EDGES_PG(BOUNDCOLL, IC)
                     FACE_NORMAL = U2D_GRID%EDGE_NORMAL(:,BOUNDCOLL,IC)
                     FACE_TANG1 = [-FACE_NORMAL(2), FACE_NORMAL(1), 0.d0]
                     FACE_TANG2 = [0.d0, 0.d0, 1.d0]
                  ELSE
                     NEIGHBOR = U3D_GRID%CELL_NEIGHBORS(BOUNDCOLL, IC)
                     IF (NEIGHBOR == -1) THEN
                        FLUIDBOUNDARY = .TRUE.
                     ELSE
                        NEIGHBORPG = U3D_GRID%CELL_PG(NEIGHBOR)
                        IF (NEIGHBORPG .NE. -1) THEN
                           IF (GRID_BC(NEIGHBORPG)%VOLUME_BC == SOLID) FLUIDBOUNDARY = .TRUE.
                        END IF
                     END IF
                     
                     FACE_PG = U3D_GRID%CELL_FACES_PG(BOUNDCOLL, IC)
                     FACE_NORMAL = U3D_GRID%FACE_NORMAL(:,BOUNDCOLL,IC)
                     FACE_TANG1 = U3D_GRID%FACE_TANG1(:,BOUNDCOLL,IC)
                     FACE_TANG2 = U3D_GRID%FACE_TANG2(:,BOUNDCOLL,IC)
                  END IF

                  ! Move to new cell
                  IF (FLUIDBOUNDARY) THEN
                     ! The particle is at the boundary of the domain
                     IF (FACE_PG .NE. -1) THEN

                        CHARGE = SPECIES(part_adv(IP)%S_ID)%CHARGE
                        IF (GRID_BC(FACE_PG)%FIELD_BC == DIELECTRIC_BC .AND. ABS(CHARGE) .GE. 1.d-6 .AND. FINAL) THEN
                           K = QE/(EPS0*EPS_SCALING**2)
                           IF (DIMS == 2) THEN
                              RHO_Q = K*CHARGE*FNUM/(ZMAX-ZMIN)
                              DO I = 1, 3
                                 VP = U2D_GRID%CELL_NODES(I,IC)
                                 PSIP = U2D_GRID%BASIS_COEFFS(1,I,IC)*part_adv(IP)%X &
                                      + U2D_GRID%BASIS_COEFFS(2,I,IC)*part_adv(IP)%Y &
                                      + U2D_GRID%BASIS_COEFFS(3,I,IC)
                                 SURFACE_CHARGE(VP) = SURFACE_CHARGE(VP) + RHO_Q*PSIP
                              END DO
                           ELSE IF (DIMS == 3) THEN
                              RHO_Q = K*CHARGE*FNUM
                              DO I = 1, 4
                                 VP = U3D_GRID%CELL_NODES(I,IC)
                                 PSIP = U3D_GRID%BASIS_COEFFS(1,I,IC)*part_adv(IP)%X &
                                      + U3D_GRID%BASIS_COEFFS(2,I,IC)*part_adv(IP)%Y &
                                      + U3D_GRID%BASIS_COEFFS(3,I,IC)*part_adv(IP)%Z &
                                      + U3D_GRID%BASIS_COEFFS(4,I,IC)
                                 SURFACE_CHARGE(VP) = SURFACE_CHARGE(VP) + RHO_Q*PSIP
                              END DO
                           END IF
                        END IF

                        ! Apply particle boundary condition
                        IF (GRID_BC(FACE_PG)%PARTICLE_BC == SPECULAR) THEN
                           IF (GRID_BC(FACE_PG)%REACT) THEN
                              CALL WALL_REACT(part_adv, IP, REMOVE_PART(IP))
                           END IF
                           
                           VDOTN = part_adv(IP)%VX*FACE_NORMAL(1) &
                                 + part_adv(IP)%VY*FACE_NORMAL(2) &
                                 + part_adv(IP)%VZ*FACE_NORMAL(3)
                           part_adv(IP)%VX = part_adv(IP)%VX - 2.*VDOTN*FACE_NORMAL(1)
                           part_adv(IP)%VY = part_adv(IP)%VY - 2.*VDOTN*FACE_NORMAL(2)
                           part_adv(IP)%VZ = part_adv(IP)%VZ - 2.*VDOTN*FACE_NORMAL(3)
                        ELSE IF (GRID_BC(FACE_PG)%PARTICLE_BC == DIFFUSE) THEN
                           IF (GRID_BC(FACE_PG)%REACT) THEN
                              CALL WALL_REACT(part_adv, IP, REMOVE_PART(IP))
                           END IF

                           S_ID = part_adv(IP)%S_ID
                           WALL_TEMP = GRID_BC(FACE_PG)%WALL_TEMP
                           CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                           WALL_TEMP, WALL_TEMP, WALL_TEMP, &
                           VDUMMY, V_TANG1, V_TANG2, SPECIES(S_ID)%MOLECULAR_MASS)

                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, WALL_TEMP, EROT)
                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, WALL_TEMP, EVIB)
                                          
                           V_PERP = FLX(0.d0, WALL_TEMP, SPECIES(S_ID)%MOLECULAR_MASS)


                           part_adv(IP)%VX = - V_PERP*FACE_NORMAL(1) &
                                             - V_TANG1*FACE_TANG1(1) &
                                             - V_TANG2*FACE_TANG2(1)
                           part_adv(IP)%VY = - V_PERP*FACE_NORMAL(2) &
                                             - V_TANG1*FACE_TANG1(2) &
                                             - V_TANG2*FACE_TANG2(2)
                           part_adv(IP)%VZ = - V_PERP*FACE_NORMAL(3) &
                                             - V_TANG1*FACE_TANG1(3) &
                                             - V_TANG2*FACE_TANG2(3)
                           
                           part_adv(IP)%EROT = EROT
                           part_adv(IP)%EVIB = EVIB

                        ELSE IF (GRID_BC(FACE_PG)%PARTICLE_BC == CLL) THEN
                           IF (GRID_BC(FACE_PG)%REACT) THEN
                              CALL WALL_REACT(part_adv, IP, REMOVE_PART(IP))
                           END IF

                           VDOTN = part_adv(IP)%VX*FACE_NORMAL(1) &
                                 + part_adv(IP)%VY*FACE_NORMAL(2) &
                                 + part_adv(IP)%VZ*FACE_NORMAL(3)

                           TANG1(1) = part_adv(IP)%VX - VDOTN*FACE_NORMAL(1)
                           TANG1(2) = part_adv(IP)%VY - VDOTN*FACE_NORMAL(2)
                           TANG1(3) = part_adv(IP)%VZ - VDOTN*FACE_NORMAL(3)

                           TANG1 = TANG1/NORM2(TANG1)

                           TANG2 = CROSS(FACE_NORMAL, TANG1)

                           VDOTTANG1 = part_adv(IP)%VX*TANG1(1) &
                                     + part_adv(IP)%VY*TANG1(2) &
                                     + part_adv(IP)%VZ*TANG1(3)

                           S_ID = part_adv(IP)%S_ID
                           WALL_TEMP = GRID_BC(FACE_PG)%WALL_TEMP

                           VRM = SQRT(2.*KB*WALL_TEMP/SPECIES(S_ID)%MOLECULAR_MASS)

                           ! Normal velocity for the CLL model
                           RN = rf()
                           DO WHILE (RN < 1.0D-13)
                              RN = rf()
                           END DO
                           R1 = SQRT(-GRID_BC(FACE_PG)%ACC_N*LOG(RN))
                           THETA1 = 2.*PI*rf()
                           DOT_NORM = VDOTN/VRM * SQRT(1 - GRID_BC(FACE_PG)%ACC_N)
                           V_PERP = VRM * SQRT(R1*R1 + DOT_NORM*DOT_NORM + 2.*R1*DOT_NORM*COS(THETA1))

                           ! Tangential velocity for the CLL model
                           RN = rf()
                           DO WHILE (RN < 1.0D-13)
                              RN = rf()
                           END DO         
                           R2 = SQRT(-GRID_BC(FACE_PG)%ACC_T*LOG(RN))
                           THETA2 = 2.*PI*rf()
                           VTANGENT = VDOTTANG1/VRM * SQRT(1 - GRID_BC(FACE_PG)%ACC_T)
                           V_TANG1 = VRM * (VTANGENT + R2 * COS(THETA2))
                           V_TANG2 = VRM * R2 * SIN(THETA2)

                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, WALL_TEMP, EROT)
                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, WALL_TEMP, EVIB)

                           part_adv(IP)%VX = V_PERP*FACE_NORMAL(1) &
                                           + V_TANG1*TANG1(1) &
                                           + V_TANG2*TANG2(1)
                           part_adv(IP)%VY = V_PERP*FACE_NORMAL(2) &
                                           + V_TANG1*TANG1(2) &
                                           + V_TANG2*TANG2(2)
                           part_adv(IP)%VZ = V_PERP*FACE_NORMAL(3) &
                                           + V_TANG1*TANG1(3) &
                                           + V_TANG2*TANG2(3)
                           
                           part_adv(IP)%EROT = EROT
                           part_adv(IP)%EVIB = EVIB

                        ELSE
                           REMOVE_PART(IP) = .TRUE.
                           part_adv(IP)%DTRIM = 0.d0
                        END IF
                     ELSE
                        REMOVE_PART(IP) = .TRUE.
                        part_adv(IP)%DTRIM = 0.d0
                     END IF
                  ELSE
                     ! The particle is crossing to another cell
                     part_adv(IP)%IC = NEIGHBOR
                     IC = NEIGHBOR
                     NCROSSINGS = NCROSSINGS + 1
                     !WRITE(*,*) 'moved particle to cell: ', IC
                  END IF
               ELSE IF (BOUNDCOLL == 0) THEN
                  CALL MOVE_PARTICLE_CN(part_adv, IP, E, DTCOLL)
                  IF (AXI .AND. DIMS == 2) CALL AXI_ROTATE_VELOCITY(part_adv, IP)
                  part_adv(IP)%DTRIM = part_adv(IP)%DTRIM - DTCOLL
               ELSE
                  ! The particle stays within the current cell. End of the motion.
                  CALL MOVE_PARTICLE_CN(part_adv, IP, E, part_adv(IP)%DTRIM)
                  IF (AXI .AND. DIMS == 2) CALL AXI_ROTATE_VELOCITY(part_adv, IP)
                  part_adv(IP)%DTRIM = 0.d0
               END IF
            
            ELSE
               CALL ERROR_ABORT('Error! C-N mover only implemented for unstructured grids.')
            END IF
            
            ! IF (tID == 4 .AND. FINAL) THEN
            !    WRITE(66332,*) NCROSSINGS, NUMSTEPS
            ! END IF

         END DO ! WHILE (DTRIM .GT. 0.)

            ! _______ APPLY Z PERIODIC BCs ________

            ! Z we should be able to do this a posteriori (remember it could be more than one depth length out!)
         IF (BOOL_Z_PERIODIC) THEN
         
            DO WHILE (part_adv(IP)%Z .GT. ZMAX)
               part_adv(IP)%Z = ZMIN + (part_adv(IP)%Z - ZMAX)
            END DO
            DO WHILE (part_adv(IP)%Z .LT. ZMIN) 
               part_adv(IP)%Z = ZMAX + (part_adv(IP)%Z - ZMIN)
            END DO
      
         END IF 

         part_adv(IP)%DTRIM = DT ! For the next timestep.


         ! IF (COMPUTE_JACOBIAN) THEN

         !    DXDEBAR = DXDEBAR * SPECIES(part_adv(IP)%S_ID)%CHARGE*QE
         !    COUNTER = 0
         !    DO I = 1, U2D_GRID%NUM_CELLS
         !       IF (DXDEBAR(I) .NE. 0.d0) COUNTER = COUNTER + 1
         !    END DO
         !    IF (ALLOCATED(COLINDICES)) DEALLOCATE(COLINDICES)
         !    IF (ALLOCATED(DXDEBAR_NZ)) DEALLOCATE(DXDEBAR_NZ)
         !    ALLOCATE(COLINDICES(COUNTER))
         !    ALLOCATE(DXDEBAR_NZ(COUNTER))
         !    J = 1
         !    DO I = 1, U2D_GRID%NUM_CELLS
         !       IF (DXDEBAR(I) .NE. 0.d0) THEN
         !          COLINDICES(J) = I - 1
         !          DXDEBAR_NZ(J) = DXDEBAR(I)
         !          J = J + 1
         !       END IF
         !    END DO
         !    IF (COUNTER > 0) THEN
         !       CALL MatSetValues(dxde,1,IC-1,COUNTER,COLINDICES,DXDEBAR_NZ,ADD_VALUES,ierr)
         !    END IF
         !    ! Tutto qui. Il resto fuori dal loop particelle.
         ! END IF


         IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 1) THEN

            DXDEX = DXDEX * SPECIES(part_adv(IP)%S_ID)%CHARGE*QE
            DXDEY = DXDEY * SPECIES(part_adv(IP)%S_ID)%CHARGE*QE
            DYDEX = DYDEX * SPECIES(part_adv(IP)%S_ID)%CHARGE*QE
            DYDEY = DYDEY * SPECIES(part_adv(IP)%S_ID)%CHARGE*QE
            !CALL MatSetValues(jac,1,IC-1,NUMSTEPS,EBARIDX,DXDEBAR,ADD_VALUES,ierr)
            DO I = 1, NUMSTEPS
               CALL MatSetValue(dxdexmat,IC-1,EBARIDX(I),DXDEX(I),ADD_VALUES,ierr)
               CALL MatSetValue(dxdeymat,IC-1,EBARIDX(I),DXDEY(I),ADD_VALUES,ierr)
               CALL MatSetValue(dydexmat,IC-1,EBARIDX(I),DYDEX(I),ADD_VALUES,ierr)
               CALL MatSetValue(dydeymat,IC-1,EBARIDX(I),DYDEY(I),ADD_VALUES,ierr)
            END DO
            ! Tutto qui. Il resto fuori dal loop particelle.
         ELSE IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 2) THEN
            DXDEBAR = DXDEBAR * SPECIES(part_adv(IP)%S_ID)%CHARGE*QE
            !CALL MatSetValues(jac,1,IC-1,NUMSTEPS,EBARIDX,DXDEBAR,ADD_VALUES,ierr)
            DO I = 1, NUMSTEPS
               CALL MatSetValue(dxde,IC-1,EBARIDX(I),DXDEBAR(I),ADD_VALUES,ierr)
            END DO

         END IF

      END DO ! End loop: DO IP = 1,NP_PROC

      !WRITE(*,*) 'Max number of steps on proc ', PROC_ID, ' = ', MAXNUMSTEPS


      IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 1) THEN

         CALL MatAssemblyBegin(dxdexmat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyBegin(dxdeymat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyBegin(dydexmat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyBegin(dydeymat,MAT_FINAL_ASSEMBLY,ierr)

         CALL MatAssemblyEnd(dxdexmat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(dxdeymat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(dydexmat,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(dydeymat,MAT_FINAL_ASSEMBLY,ierr)

         !CALL PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER, &
         !                         PETSC_NULL_CHARACTER,0,0,684,684,viewer,ierr)
         !CALL MatView(dxde,viewer,ierr)
         !CALL PetscSleep(5.d0,ierr)

         !CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,'mat.output',viewer,ierr)
         !CALL MatView(dxde,viewer,ierr)
         !CALL PetscSleep(5.d0,ierr)
         !CALL PetscViewerDestroy(viewer,ierr)
   
         CALL MatGetOwnershipRange(dxdexmat,first_row,last_row,ierr)

         IF (DIMS == 2) THEN
            DO I = first_row, last_row-1
               CALL MatGetRow(dxdexmat,I,ncols,cols,dxdexvals,ierr)
               NCOLSXX = ncols
               ALLOCATE(COLSXX, SOURCE=cols)
               ALLOCATE(DXDEXV, SOURCE=dxdexvals)
               CALL MatRestoreRow(dxdexmat,I,ncols,cols,dxdexvals,ierr)
               
               CALL MatGetRow(dxdeymat,I,ncols,cols,dxdeyvals,ierr)
               NCOLSXY = ncols
               ALLOCATE(COLSXY, SOURCE=cols)
               ALLOCATE(DXDEYV, SOURCE=dxdeyvals)
               CALL MatRestoreRow(dxdeymat,I,ncols,cols,dxdeyvals,ierr)

               CALL MatGetRow(dydexmat,I,ncols,cols,dydexvals,ierr)
               NCOLSYX = ncols
               ALLOCATE(COLSYX, SOURCE=cols)
               ALLOCATE(DYDEXV, SOURCE=dydexvals)
               CALL MatRestoreRow(dydexmat,I,ncols,cols,dydexvals,ierr)

               CALL MatGetRow(dydeymat,I,ncols,cols,dydeyvals,ierr)
               NCOLSYY = ncols
               ALLOCATE(COLSYY, SOURCE=cols)
               ALLOCATE(DYDEYV, SOURCE=dydeyvals)
               CALL MatRestoreRow(dydeymat,I,ncols,cols,dydeyvals,ierr)

               
               !WRITE(*,*) 'Row ', I, ' ncols ', ncols, ' cols ', cols, ' vals ', vals, ' ierr ', ierr
               ! MatGetRow(Mat A,PetscInt row, PetscInt *ncols,const PetscInt (*cols)[],const PetscScalar (*vals)[]);

               V1I = U2D_GRID%CELL_NODES(1,I+1)
               V2I = U2D_GRID%CELL_NODES(2,I+1)
               V3I = U2D_GRID%CELL_NODES(3,I+1)

               DPSI1DX = U2D_GRID%BASIS_COEFFS(1,1,I+1)
               DPSI2DX = U2D_GRID%BASIS_COEFFS(1,2,I+1)
               DPSI3DX = U2D_GRID%BASIS_COEFFS(1,3,I+1)
               DPSI1DY = U2D_GRID%BASIS_COEFFS(2,1,I+1)
               DPSI2DY = U2D_GRID%BASIS_COEFFS(2,2,I+1)
               DPSI3DY = U2D_GRID%BASIS_COEFFS(2,3,I+1)

               DO JJ = 1, NCOLSXX
                  J = COLSXX(JJ)

                  V1J = U2D_GRID%CELL_NODES(1,J+1)
                  V2J = U2D_GRID%CELL_NODES(2,J+1)
                  V3J = U2D_GRID%CELL_NODES(3,J+1)

                  DPSJ1DX = U2D_GRID%BASIS_COEFFS(1,1,J+1)
                  DPSJ2DX = U2D_GRID%BASIS_COEFFS(1,2,J+1)
                  DPSJ3DX = U2D_GRID%BASIS_COEFFS(1,3,J+1)
                  DPSJ1DY = U2D_GRID%BASIS_COEFFS(2,1,J+1)
                  DPSJ2DY = U2D_GRID%BASIS_COEFFS(2,2,J+1)
                  DPSJ3DY = U2D_GRID%BASIS_COEFFS(2,3,J+1)
   
                  VALXX = - DXDEXV(JJ)/EPS0/(ZMAX-ZMIN)*FNUM
                  VALXY = - DXDEYV(JJ)/EPS0/(ZMAX-ZMIN)*FNUM
                  VALYX = - DYDEXV(JJ)/EPS0/(ZMAX-ZMIN)*FNUM
                  VALYY = - DYDEYV(JJ)/EPS0/(ZMAX-ZMIN)*FNUM

                  IF (.NOT. IS_DIRICHLET(V1I-1)) THEN
                     CALL MatSetValues(jac,1,V1I-1,1,V1J-1,VALXX*DPSI1DX*DPSJ1DX + VALXY*DPSI1DX*DPSJ1DY + &
                                                           VALYX*DPSI1DY*DPSJ1DX + VALYY*DPSI1DY*DPSJ1DY,ADD_VALUES,ierr)
                     CALL MatSetValues(jac,1,V1I-1,1,V2J-1,VALXX*DPSI1DX*DPSJ2DX + VALXY*DPSI1DX*DPSJ2DY + &
                                                           VALYX*DPSI1DY*DPSJ2DX + VALYY*DPSI1DY*DPSJ2DY,ADD_VALUES,ierr)
                     CALL MatSetValues(jac,1,V1I-1,1,V3J-1,VALXX*DPSI1DX*DPSJ3DX + VALXY*DPSI1DX*DPSJ3DY + &
                                                           VALYX*DPSI1DY*DPSJ3DX + VALYY*DPSI1DY*DPSJ3DY,ADD_VALUES,ierr)
                  END IF
                  IF (.NOT. IS_DIRICHLET(V2I-1)) THEN
                     CALL MatSetValues(jac,1,V2I-1,1,V1J-1,VALXX*DPSI2DX*DPSJ1DX + VALXY*DPSI2DX*DPSJ1DY + &
                                                           VALYX*DPSI2DY*DPSJ1DX + VALYY*DPSI2DY*DPSJ1DY,ADD_VALUES,ierr)
                     CALL MatSetValues(jac,1,V2I-1,1,V2J-1,VALXX*DPSI2DX*DPSJ2DX + VALXY*DPSI2DX*DPSJ2DY + &
                                                           VALYX*DPSI2DY*DPSJ2DX + VALYY*DPSI2DY*DPSJ2DY,ADD_VALUES,ierr)
                     CALL MatSetValues(jac,1,V2I-1,1,V3J-1,VALXX*DPSI2DX*DPSJ3DX + VALXY*DPSI2DX*DPSJ3DY + &
                                                           VALYX*DPSI2DY*DPSJ3DX + VALYY*DPSI2DY*DPSJ3DY,ADD_VALUES,ierr)
                  END IF
                  IF (.NOT. IS_DIRICHLET(V3I-1)) THEN
                     CALL MatSetValues(jac,1,V3I-1,1,V1J-1,VALXX*DPSI3DX*DPSJ1DX + VALXY*DPSI3DX*DPSJ1DY + &
                                                           VALYX*DPSI3DY*DPSJ1DX + VALYY*DPSI3DY*DPSJ1DY,ADD_VALUES,ierr)
                     CALL MatSetValues(jac,1,V3I-1,1,V2J-1,VALXX*DPSI3DX*DPSJ2DX + VALXY*DPSI3DX*DPSJ2DY + &
                                                           VALYX*DPSI3DY*DPSJ2DX + VALYY*DPSI3DY*DPSJ2DY,ADD_VALUES,ierr)
                     CALL MatSetValues(jac,1,V3I-1,1,V3J-1,VALXX*DPSI3DX*DPSJ3DX + VALXY*DPSI3DX*DPSJ3DY + &
                                                           VALYX*DPSI3DY*DPSJ3DX + VALYY*DPSI3DY*DPSJ3DY,ADD_VALUES,ierr)
                  END IF

               END DO

               DEALLOCATE(COLSXX)
               DEALLOCATE(COLSXY)
               DEALLOCATE(COLSYX)
               DEALLOCATE(COLSYY)
               DEALLOCATE(DXDEXV)
               DEALLOCATE(DXDEYV)
               DEALLOCATE(DYDEXV)
               DEALLOCATE(DYDEYV)

            END DO
         ELSE IF (DIMS == 3) THEN
            CALL ERROR_ABORT('Not supported!')   
         END IF

         CALL MatDestroy(dxdexmat,ierr)
         CALL MatDestroy(dxdeymat,ierr)
         CALL MatDestroy(dydexmat,ierr)
         CALL MatDestroy(dydeymat,ierr)

      ELSE IF (COMPUTE_JACOBIAN .AND. JACOBIAN_TYPE == 2) THEN


         CALL MatAssemblyBegin(dxde,MAT_FINAL_ASSEMBLY,ierr)
         CALL MatAssemblyEnd(dxde,MAT_FINAL_ASSEMBLY,ierr)

         !CALL PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER, &
         !                         PETSC_NULL_CHARACTER,0,0,684,684,viewer,ierr)
         !CALL MatView(dxde,viewer,ierr)
         !CALL PetscSleep(5.d0,ierr)

         !CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,'mat.output',viewer,ierr)
         !CALL MatView(dxde,viewer,ierr)
         !CALL PetscSleep(5.d0,ierr)
         !CALL PetscViewerDestroy(viewer,ierr)
   
         CALL MatGetOwnershipRange(dxde,first_row,last_row,ierr)

         IF (DIMS == 2) THEN
            DO I = first_row, last_row-1
               CALL MatGetRow(dxde,I,ncols,cols,vals,ierr)

               DO JJ = 1, ncols
                  J = cols(JJ)

                  VAL = - vals(JJ)/EPS0/(ZMAX-ZMIN)*FNUM
                  DO NI = 1, 3
                     VNI = U2D_GRID%CELL_NODES(NI,I+1)
                     IF (.NOT. IS_DIRICHLET(VNI - 1)) THEN
                        DO NJ = 1, 3
                           VNJ = U2D_GRID%CELL_NODES(NJ,J+1)
                           CALL MatSetValues(jac,1,VNI-1,1,VNJ-1,VAL* &
                            (U2D_GRID%BASIS_COEFFS(1,NI,I+1)*U2D_GRID%BASIS_COEFFS(1,NJ,J+1) &
                           + U2D_GRID%BASIS_COEFFS(2,NI,I+1)*U2D_GRID%BASIS_COEFFS(2,NJ,J+1)),ADD_VALUES,ierr)
                        END DO
                     END IF
                  END DO

               END DO

               CALL MatRestoreRow(dxde,I,ncols,cols,vals,ierr)

            END DO
         ELSE IF (DIMS == 3) THEN
            DO I = first_row, last_row-1
               CALL MatGetRow(dxde,I,ncols,cols,vals,ierr)
   
               DO JJ = 1, ncols
                  J = cols(JJ)

                  VAL = - vals(JJ)/EPS0*FNUM
                  DO NI = 1, 4
                     VNI = U3D_GRID%CELL_NODES(NI,I+1)
                     IF (.NOT. IS_DIRICHLET(VNI - 1)) THEN
                        DO NJ = 1, 4
                           VNJ = U3D_GRID%CELL_NODES(NJ,J+1)
                           CALL MatSetValues(jac,1,VNI-1,1,VNJ-1,VAL* &
                            (U3D_GRID%BASIS_COEFFS(1,NI,I+1)*U3D_GRID%BASIS_COEFFS(1,NJ,J+1) &
                           + U3D_GRID%BASIS_COEFFS(2,NI,I+1)*U3D_GRID%BASIS_COEFFS(2,NJ,J+1) &
                           + U3D_GRID%BASIS_COEFFS(3,NI,I+1)*U3D_GRID%BASIS_COEFFS(3,NJ,J+1)),ADD_VALUES,ierr)
                        END DO
                     END IF
                  END DO
   
               END DO
   
               CALL MatRestoreRow(dxde,I,ncols,cols,vals,ierr)
   
            END DO
   
         END IF

         ! CALL MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
         ! CALL MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)

         !CALL MatMatMult(drhodxx,dxde,MAT_INITIAL_MATRIX,PETSC_DEFAULT,drhodex,ierr)
         !CALL MatMatMult(drhodex,dexdphi,MAT_INITIAL_MATRIX,PETSC_DEFAULT,drhodphi,ierr)

         CALL MatDestroy(dxde,ierr)

      END IF

      ! ==============
      ! ========= Now remove part_adv that are out of the domain and reorder the array. 
      ! ========= Do this in the end of the advection step, since reordering the array 
      ! ========= mixes the part_adv.
      ! ==============
      IF (FINAL) THEN
         IP = NP_PROC
         DO WHILE (IP .GE. 1)

            ! Is particle IP out of the domain? Then remove it!
            IF (REMOVE_PART(IP)) THEN
               CALL REMOVE_PARTICLE_ARRAY(IP, part_adv, NP_PROC)
            ELSE IF (GRID_TYPE .NE. UNSTRUCTURED) THEN
               CALL CELL_FROM_POSITION(part_adv(IP)%X, part_adv(IP)%Y, IC)
               OLD_IC = part_adv(IP)%IC
               part_adv(IP)%IC = IC
               ! If cell-based particle weight has changed,
               IF (BOOL_RADIAL_WEIGHTING .AND. (IC .NE. OLD_IC)) THEN
                  WEIGHT_RATIO = CELL_FNUM(OLD_IC)/CELL_FNUM(IC)
                  ! If the weight has increased, the particle may be removed
                  IF (WEIGHT_RATIO .LT. 1.d0) THEN
                     rfp = rf()
                     IF (WEIGHT_RATIO .LT. rfp) CALL REMOVE_PARTICLE_ARRAY(IP, part_adv, NP_PROC)
                  ELSE IF (WEIGHT_RATIO .GT. 1.d0) THEN
                  ! If the weight has decreased, more part_adv may be added
                     rfp = rf()
                     DO WHILE ((WEIGHT_RATIO - 1.) .GT. rfp)
                        CALL INIT_PARTICLE(part_adv(IP)%X, part_adv(IP)%Y, part_adv(IP)%Z, &
                        part_adv(IP)%VX, part_adv(IP)%VY,part_adv(IP)%VZ, &
                        part_adv(IP)%EROT,part_adv(IP)%EVIB,part_adv(IP)%S_ID,IC,DT,NEWparticle)
                        CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, part_adv)
                        WEIGHT_RATIO = WEIGHT_RATIO - 1.
                     END DO
                  END IF
               END IF
            END IF

            IP = IP - 1

         END DO
      END IF
      DEALLOCATE(REMOVE_PART)


      IF (FINAL) THEN
         ! Tally surface collisions from all processes
         CALL MPI_REDUCE(LOCAL_BOUNDARY_COLL_COUNT, BOUNDARY_COLL_COUNT, 4*N_SPECIES, MPI_INTEGER, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)

         CALL MPI_REDUCE(LOCAL_WALL_COLL_COUNT, WALL_COLL_COUNT, N_WALLS*N_SPECIES, MPI_INTEGER, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
      END IF

      DEALLOCATE(LOCAL_BOUNDARY_COLL_COUNT)
      DEALLOCATE(LOCAL_WALL_COLL_COUNT)

      ! IF (tID == 4 .AND. FINAL) THEN
      !    CLOSE(66332)
      ! END IF
      
      !CALL VecDestroy(dxdexvals,ierr)
      !CALL VecDestroy(dxdeyvals,ierr)
      !CALL VecDestroy(dydexvals,ierr)
      !CALL VecDestroy(dydeyvals,ierr)
      !CALL VecDestroy(vals,ierr)
      !CALL VecDestroy(cols,ierr)


   END SUBROUTINE ADVECT_CN


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MOVE_PARTICLE_CN -> Move the particles !!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE MOVE_PARTICLE_CN(part_adv, IP, E, TIME)

      ! Moves particle with index IP for time TIME.
      ! For now simply rectilinear movement, will have to include Lorentz's force

      IMPLICIT NONE

      INTEGER, INTENT(IN)      :: IP
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: E
      REAL(KIND=8), INTENT(IN) :: TIME
      REAL(KIND=8), DIMENSION(3) :: ACC
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: part_adv


      ! IF (part_adv(IP)%Z .NE. 0.d0) WRITE(*,*) 'Particle Z coordinate was not zero but ', part_adv(IP)%Z

      ACC = E*SPECIES(part_adv(IP)%S_ID)%CHARGE*QE/SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS

      part_adv(IP)%X = part_adv(IP)%X + part_adv(IP)%VX*TIME + 0.5*ACC(1)*TIME*TIME
      part_adv(IP)%Y = part_adv(IP)%Y + part_adv(IP)%VY*TIME + 0.5*ACC(2)*TIME*TIME
      part_adv(IP)%Z = part_adv(IP)%Z + part_adv(IP)%VZ*TIME + 0.5*ACC(3)*TIME*TIME

      part_adv(IP)%VX = part_adv(IP)%VX + ACC(1) * TIME
      part_adv(IP)%VY = part_adv(IP)%VY + ACC(2) * TIME
      part_adv(IP)%VZ = part_adv(IP)%VZ + ACC(3) * TIME

   END SUBROUTINE MOVE_PARTICLE_CN


   SUBROUTINE AXI_ROTATE_VELOCITY(part_adv, IP)

      IMPLICIT NONE

      INTEGER, INTENT(IN)      :: IP
      REAL(KIND=8) :: R, COSTHETA, SINTHETA, VY, VZ
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: part_adv
               
      R = SQRT((part_adv(IP)%Y)**2 + (part_adv(IP)%Z)**2)
      ! Rotate velocity vector back to x-y plane.
      SINTHETA = part_adv(IP)%Z / R
      COSTHETA = SIGN(SQRT(1.-SINTHETA*SINTHETA), part_adv(IP)%Y)

      part_adv(IP)%Z = 0.d0

      VZ = part_adv(IP)%VZ
      VY = part_adv(IP)%VY
      part_adv(IP)%VZ = COSTHETA*VZ - SINTHETA*VY
      part_adv(IP)%VY = SINTHETA*VZ + COSTHETA*VY

   END SUBROUTINE AXI_ROTATE_VELOCITY
 

   SUBROUTINE APPLY_E_FIELD(JP, E)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: E
      INTEGER, INTENT(IN) :: JP
      REAL(KIND=8), DIMENSION(4) :: WEIGHTS
      INTEGER, DIMENSION(4) :: INDICES, INDI, INDJ

      IF (GRID_TYPE == UNSTRUCTURED) THEN
         IF (PIC_TYPE == SEMIIMPLICIT .OR. PIC_TYPE == FULLYIMPLICIT) THEN
            E = EBAR_FIELD(:, 1, particles(JP)%IC)
         ELSE
            E = E_FIELD(:, 1, particles(JP)%IC)
         END IF
      ELSE

         CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)
         IF (DIMS == 2) THEN
            E = WEIGHTS(1)*E_FIELD(:, INDJ(1), INDI(1)) + &
                WEIGHTS(2)*E_FIELD(:, INDJ(2), INDI(2)) + &
                WEIGHTS(3)*E_FIELD(:, INDJ(3), INDI(3)) + &
                WEIGHTS(4)*E_FIELD(:, INDJ(4), INDI(4))
         ELSE
            E = WEIGHTS(1)*E_FIELD(:, INDJ(1), INDI(1)) + &
                WEIGHTS(2)*E_FIELD(:, INDJ(2), INDI(2))
         END IF

      END IF
   END SUBROUTINE APPLY_E_FIELD



   SUBROUTINE APPLY_B_FIELD(JP, B)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: B
      INTEGER, INTENT(IN) :: JP
      INTEGER :: J, VJ, IC
      REAL(KIND=8) :: PSIJ

      IC = particles(JP)%IC
      B = 0.d0
      IF (GRID_TYPE == UNSTRUCTURED) THEN
         DO J = 1, 3
            VJ = U2D_GRID%CELL_NODES(J, IC)
            PSIJ = particles(JP)%X*U2D_GRID%BASIS_COEFFS(1,J,IC) &
                 + particles(JP)%Y*U2D_GRID%BASIS_COEFFS(2,J,IC) &
                 + U2D_GRID%BASIS_COEFFS(3,J,IC)
            B = B + B_FIELD(:, 1, VJ)*PSIJ
         END DO
      END IF


   END SUBROUTINE APPLY_B_FIELD


   SUBROUTINE WALL_REACT(part_adv, IP, REMOVE)
      
      IMPLICIT NONE

      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: part_adv
      INTEGER, INTENT(IN) :: IP
      LOGICAL, INTENT(OUT) :: REMOVE
      INTEGER :: JS, JR, JP
      REAL(KIND=8) :: PROB_SCALE, VEL_SCALE

      JS = part_adv(IP)%S_ID
      PROB_SCALE = 1.
      REMOVE = .FALSE.
      DO JR = 1, N_WALL_REACTIONS
         IF (WALL_REACTIONS(JR)%R_SP_ID == JS) THEN
            IF ( rf() .LE. WALL_REACTIONS(JR)%PROB/PROB_SCALE ) THEN
               
               IF (WALL_REACTIONS(JR)%N_PROD == 0) THEN
                  REMOVE = .TRUE.
                  part_adv(IP)%DTRIM = 0.
               ELSE IF (WALL_REACTIONS(JR)%N_PROD == 1) THEN
                  JP = WALL_REACTIONS(JR)%P1_SP_ID
                  part_adv(IP)%S_ID = JP
                  VEL_SCALE = SPECIES(JS)%MOLECULAR_MASS/SPECIES(JP)%MOLECULAR_MASS
                  part_adv(IP)%VX = part_adv(IP)%VX*VEL_SCALE
                  part_adv(IP)%VY = part_adv(IP)%VY*VEL_SCALE
                  part_adv(IP)%VZ = part_adv(IP)%VZ*VEL_SCALE
               ELSE
                  CALL ERROR_ABORT('Number of products in wall reaction not supported.')
               END IF

            ELSE
               PROB_SCALE = PROB_SCALE - WALL_REACTIONS(JR)%PROB
            END IF

         END IF
      END DO


   END SUBROUTINE WALL_REACT


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE SETUP_POISSON -> Solves the Poisson equation with the RHS RHS !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE SETUP_POISSON

      IMPLICIT NONE

      CALL KSPDestroy(ksp,ierr)
      CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
      CALL KSPSetOperators(ksp,Amat,Amat,ierr)

      ! Set solver type to direct (LU). Could be MUMPS but must install with --download-mumps
      !CALL KSPGetPC(ksp,pc,ierr)
      !CALL PCSetType(pc,PCLU,ierr)
      CALL KSPSetFromOptions(ksp,ierr)

      !abstol = 1.d-15
      !rtol = 1.d-15
      !maxit = 50
      !CALL KSPSetTolerances(ksp,rtol,abstol,PETSC_DEFAULT_REAL,maxit,ierr)

   END SUBROUTINE SETUP_POISSON


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE SOLVE_POISSON -> Solves the Poisson equation with the RHS RHS !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE SOLVE_POISSON

      IMPLICIT NONE

      ! Solve the linear system  (Amat*PHI_FIELD = bvec)
      CALL KSPSolve(ksp,bvec,xvec,ierr)

      CALL KSPGetConvergedReason(ksp,reason,ierr)
      IF (PROC_ID == 0) WRITE(*,*) 'KSPConvergedReason = ', reason

      CALL VecScatterCreateToAll(xvec,ctx,xvec_seq,ierr)
      CALL VecScatterBegin(ctx,xvec,xvec_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterEnd(ctx,xvec,xvec_seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterDestroy(ctx, ierr)

      CALL VecGetArrayReadF90(xvec_seq,PHI_FIELD_TEMP,ierr)
      IF (ALLOCATED(PHI_FIELD)) DEALLOCATE(PHI_FIELD)
      ALLOCATE(PHI_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(xvec_seq,PHI_FIELD_TEMP,ierr)

      CALL VecScatterDestroy(ctx,ierr)
      CALL VecDestroy(xvec_seq,ierr)

      !!CALL KSPDestroy(ksp,ierr)

      !CALL VecRestoreArrayReadF90(X_SEQ,PHI_FIELD,ierr)

      !IF (PROC_ID == 0) WRITE(*,*) 'Max PHI_FIELD= ', MAXVAL(PHI_FIELD), 'Min PHI_FIELD= ', MINVAL(PHI_FIELD)

   END SUBROUTINE SOLVE_POISSON


   SUBROUTINE COMPUTE_E_FIELD
      
      IMPLICIT NONE

      REAL(KIND=8) :: HX, HY
      INTEGER :: I, J, P, VP, ID
      INTEGER :: SIZE
      INTEGER :: ICENTER, INORTH, ISOUTH, IEAST, IWEST

      HX = (XMAX-XMIN)/DBLE(NX)
      HY = (YMAX-YMIN)/DBLE(NY)

      SIZE = NNODES

      ! Reshape the linear array into a 2D array
      IF (GRID_TYPE == UNSTRUCTURED) THEN

         ! Compute the electric field at grid points
         IF (DIMS == 2) THEN

            E_FIELD = 0.d0
            EBAR_FIELD = 0.d0
            DO I = 1, NCELLS
               DO ID = 1, 2
                  DO P = 1, 3
                     VP = U2D_GRID%CELL_NODES(P,I)
                     E_FIELD(ID,1,I) =    E_FIELD(ID,1,I)    - PHI_FIELD(VP)   *U2D_GRID%BASIS_COEFFS(ID,P,I)
                     EBAR_FIELD(ID,1,I) = EBAR_FIELD(ID,1,I) - PHIBAR_FIELD(VP)*U2D_GRID%BASIS_COEFFS(ID,P,I)
                  END DO
               END DO
            END DO

         ELSE IF (DIMS == 3) THEN

            E_FIELD = 0.d0
            EBAR_FIELD = 0.d0
            DO I = 1, NCELLS
               DO ID = 1, 3
                  DO P = 1, 4
                     VP = U3D_GRID%CELL_NODES(P,I)
                     E_FIELD(ID,1,I) =    E_FIELD(ID,1,I)    - PHI_FIELD(VP)   *U3D_GRID%BASIS_COEFFS(ID,P,I)
                     EBAR_FIELD(ID,1,I) = EBAR_FIELD(ID,1,I) - PHIBAR_FIELD(VP)*U3D_GRID%BASIS_COEFFS(ID,P,I)
                  END DO
               END DO
            END DO

         END IF

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
                  E_FIELD(1,J,I) = (PHI_FIELD(ICENTER)-PHI_FIELD(IEAST))/HX
               ELSE IF (I == NPX-1) THEN ! Right boundary
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = XSIZE(NPX-1)
                  E_FIELD(1,J,I) = (PHI_FIELD(IWEST)-PHI_FIELD(ICENTER))/HX
               ELSE ! Interior point
                  IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HX = 0.5*(XSIZE(I)+XSIZE(I+1))
                  E_FIELD(1,J,I) = 0.5*(PHI_FIELD(IWEST)-PHI_FIELD(IEAST))/HX
               END IF
               IF (DIMS == 2) THEN
                  IF (J == 0) THEN ! Bottom boundary
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = YSIZE(1)
                     E_FIELD(2,J,I) = (PHI_FIELD(ICENTER)-PHI_FIELD(INORTH))/HY
                  ELSE IF (J == NPY-1) THEN ! Top boundary
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = YSIZE(NPY-1)
                     E_FIELD(2,J,I) = (PHI_FIELD(ISOUTH)-PHI_FIELD(ICENTER))/HY
                  ELSE ! Interior point
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) HY = 0.5*(YSIZE(J)+YSIZE(J+1))
                     E_FIELD(2,J,I) = 0.5*(PHI_FIELD(ISOUTH)-PHI_FIELD(INORTH))/HY
                  END IF
               ELSE
                  E_FIELD(2,J,I) = 0.d0
               END IF
               E_FIELD(3,J,I) = 0.d0
            END DO
         END DO
         
      END IF

   END SUBROUTINE COMPUTE_E_FIELD


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
         INDICES(1) = IC-1
         INDICES(2) = IC
         INDICES(4) = -1
         INDICES(3) = -1

         ! Left, right
         WEIGHTS(1) = (1.-FX)
         WEIGHTS(2) = FX
         WEIGHTS(3) = -1
         WEIGHTS(4) = -1
      END IF
   END SUBROUTINE COMPUTE_WEIGHTS




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DEPOSIT_CHARGE -> Deposits the charge of particles on grid points !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DEPOSIT_CHARGE(part_adv)
      
      IMPLICIT NONE
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), INTENT(IN) :: part_adv
      INTEGER :: JP, I, IC

      REAL(KIND=8) :: K, RHO_Q, CHARGE
      REAL(KIND=8) :: VOL, CFNUM
      REAL(KIND=8), DIMENSION(4) :: WEIGHTS
      INTEGER, DIMENSION(4) :: INDICES, INDI, INDJ
      INTEGER :: SIZE, P, VP
      REAL(KIND=8) :: PSIP

      K = QE/(EPS0*EPS_SCALING**2) ! [V m] Elementary charge / Dielectric constant of vacuum

      RHS = 0.d0

      DO JP = 1, NP_PROC
         CHARGE = SPECIES(part_adv(JP)%S_ID)%CHARGE
         IF (ABS(CHARGE) .LT. 1.d-6) CYCLE

         IF (GRID_TYPE == UNSTRUCTURED) THEN 
            IC = part_adv(JP)%IC
            IF (DIMS == 2) THEN
               RHO_Q = K*CHARGE*FNUM/(ZMAX-ZMIN)
               DO P = 1, 3
                  VP = U2D_GRID%CELL_NODES(P,IC) - 1
                  PSIP = U2D_GRID%BASIS_COEFFS(1,P,IC)*part_adv(JP)%X &
                       + U2D_GRID%BASIS_COEFFS(2,P,IC)*part_adv(JP)%Y &
                       + U2D_GRID%BASIS_COEFFS(3,P,IC)
                  RHS(VP) = RHS(VP) + RHO_Q*PSIP
               END DO
            ELSE IF (DIMS == 3) THEN
               RHO_Q = K*CHARGE*FNUM
               DO P = 1, 4
                  VP = U3D_GRID%CELL_NODES(P,IC) - 1
                  PSIP = U3D_GRID%BASIS_COEFFS(1,P,IC)*part_adv(JP)%X &
                       + U3D_GRID%BASIS_COEFFS(2,P,IC)*part_adv(JP)%Y &
                       + U3D_GRID%BASIS_COEFFS(3,P,IC)*part_adv(JP)%Z &
                       + U3D_GRID%BASIS_COEFFS(4,P,IC)
                  RHS(VP) = RHS(VP) + RHO_Q*PSIP
               END DO
            END IF

         ELSE

            CALL COMPUTE_WEIGHTS(JP, WEIGHTS, INDICES, INDI, INDJ)

            IF (GRID_TYPE == RECTILINEAR_UNIFORM .AND. .NOT. AXI) THEN
               VOL = CELL_VOL
            ELSE
               VOL = CELL_VOLUMES(part_adv(JP)%IC)
            END IF

            CFNUM = FNUM
            IF (BOOL_RADIAL_WEIGHTING) CFNUM = CELL_FNUM(part_adv(JP)%IC)         

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

      RHS = RHS + SURFACE_CHARGE

      SIZE = NNODES

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


      DO I = 0, NNODES-1
         IF (IS_DIRICHLET(I)) THEN
            RHS(I) = DIRICHLET(I)
         ELSE IF (IS_NEUMANN(I)) THEN
            RHS(I) = RHS(I) + NEUMANN(I)
         END IF
      END DO

   END SUBROUTINE DEPOSIT_CHARGE


END MODULE fully_implicit
