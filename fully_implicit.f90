
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

   Mat Amat
   Vec bvec, xvec, X_SEQ, rvec
   VecScatter ctx
   KSP ksp, kspnk
   PetscInt one,f9,f6,f30
   PetscInt Istart, Iend
   PetscReal val, norm
   PetscScalar, POINTER :: PHI_FIELD_TEMP(:)
   !PetscScalar, POINTER :: PHIBAR_FIELD(:)
   !PetscScalar, POINTER :: PHI_FIELD_OLD(:)
   KSPConvergedReason reason
   PC pc, pcnk
   SNES snes


   CONTAINS
   

   SUBROUTINE DEPOSIT_CHARGE_CN()
      
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
         CHARGE = SPECIES(part_adv(JP)%S_ID)%CHARGE
         IF (ABS(CHARGE) .LT. 1.d-6) CYCLE

         IF (GRID_TYPE == UNSTRUCTURED) THEN 
            IC = part_adv(JP)%IC
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
            XP = part_adv(JP)%X
            YP = part_adv(JP)%Y
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

   END SUBROUTINE DEPOSIT_CHARGE_CN



   SUBROUTINE SOLVE_POISSON_FULLY_IMPLICIT

      PetscScalar, POINTER :: xvec_l(:)

      CALL SNESCreate(PETSC_COMM_WORLD,snes,ierr)
      CALL VecCreate(PETSC_COMM_WORLD,rvec,ierr)
      CALL VecSetSizes(rvec,PETSC_DECIDE,U2D_GRID%NUM_NODES,ierr)
      CALL VecSetFromOptions(rvec, ierr)
      CALL VecDuplicate(rvec, xvec, ierr)

      CALL SNESSetFunction(snes,rvec,FormFunction,0,ierr)


      CALL SNESGetKSP(snes,kspnk,ierr)
      CALL KSPGetPC(kspnk,pcnk,ierr)
      CALL PCSetType(pcnk,PCNONE,ierr)
      ! tol = 1.e-4
      ! CALL KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,i20,ierr)
      !CALL SNESSetFromOptions(snes,ierr)

      !  Note: The user should initialize the vector, x, with the initial guess
      !  for the nonlinear solver prior to calling SNESSolve().  In particular,
      !  to employ an initial guess of zero, the user should explicitly set
      !  this vector to zero by calling VecSet().
      CALL DEPOSIT_CHARGE
      CALL SOLVE_POISSON
      CALL VecGetOwnershipRange(xvec,Istart,Iend,ierr)
      CALL VecGetArrayF90(xvec,xvec_l,ierr)
      WRITE(*,*) 'ON PROC ', PROC_ID, ' Istart ', Istart, ' Iend= ', Iend
      xvec_l = PHI_FIELD(Istart+1:Iend)
      CALL VecRestoreArrayF90(xvec,xvec_l,ierr)

      CALL SNESSolve(snes,PETSC_NULL_VEC,xvec,ierr)

      CALL VecScatterCreateToAll(xvec,ctx,X_SEQ,ierr)
      CALL VecScatterBegin(ctx,xvec,X_SEQ,INSERT_VALUES,SCATTER_FORWARD)
      CALL VecScatterEnd(ctx,xvec,X_SEQ,INSERT_VALUES,SCATTER_FORWARD)

      CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD_TEMP,ierr)
      IF (ALLOCATED(PHI_FIELD)) DEALLOCATE(PHI_FIELD)
      ALLOCATE(PHI_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(X_SEQ,PHI_FIELD_TEMP,ierr)

      !CALL VecGetArrayReadF90(xvec,xvec_l,ierr)
      !PHI_FIELD = xvec_l
      !CALL VecRestoreArrayReadF90(xvec,xvec_l,ierr)

      CALL VecDestroy(xvec, ierr)
      CALL VecDestroy(rvec, ierr)
      CALL SNESDestroy(snes, ierr)

      WRITE(*,*) 'We get here on PROC ', PROC_ID
   END SUBROUTINE SOLVE_POISSON_FULLY_IMPLICIT




   SUBROUTINE FormFunction(snes,x,f,dummy,ierr_l)

      IMPLICIT NONE

      SNES snes
      Vec x,f
      PetscErrorCode ierr_l
      INTEGER dummy(*)
      PetscScalar, POINTER :: RESIDUAL(:)

      CALL VecScatterCreateToAll(x,ctx,X_SEQ,ierr)
      CALL VecScatterBegin(ctx,x,X_SEQ,INSERT_VALUES,SCATTER_FORWARD)
      CALL VecScatterEnd(ctx,x,X_SEQ,INSERT_VALUES,SCATTER_FORWARD)

      ! CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD,ierr)
      CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD_TEMP,ierr)
      IF (ALLOCATED(PHI_FIELD)) DEALLOCATE(PHI_FIELD)
      ALLOCATE(PHI_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(X_SEQ,PHI_FIELD_TEMP,ierr)

      !  Advect the particles using the guessed new potential
      ALLOCATE(part_adv, SOURCE = particles)
      CALL ADVECT_CN(.FALSE.)
      CALL DEPOSIT_CHARGE_CN
      DEALLOCATE(part_adv)

      ALLOCATE(PHI_FIELD_OLD, SOURCE = PHI_FIELD)
      ! Compute the new potential from the charge distribution
      CALL SOLVE_POISSON

      ! Compute the residual
      CALL VecGetOwnershipRange(f,Istart,Iend,ierr)

      CALL VecGetArrayF90(f,RESIDUAL,ierr_l)
      RESIDUAL = PHI_FIELD_OLD(Istart+1:Iend) - PHI_FIELD(Istart+1:Iend)
      CALL VecRestoreArrayF90(f,RESIDUAL,ierr_l)

      WRITE(*,*) 'On PROC ', PROC_ID, ' PHI_FIELD = ', PHI_FIELD

      CALL VecNorm(f,NORM_2,norm,ierr)
      WRITE(*,*) 'FormFunction Called, ||RESIDUAL|| = ', norm

      DEALLOCATE(PHI_FIELD_OLD)

   END SUBROUTINE FormFunction







   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ADVECT_CN -> Advects particles in the domain using a Crank-Nicholson scheme for fully implicit !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ADVECT_CN(FINAL)

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

      LOGICAL, INTENT(IN) :: FINAL

      INTEGER      :: IP, IC, I, J, SOL, OLD_IC
      INTEGER      :: BOUNDCOLL, WALLCOLL, GOODSOL, EDGE_PG
      REAL(KIND=8) :: DTCOLL, TOTDTCOLL, CANDIDATE_DTCOLL, rfp
      REAL(KIND=8) :: COEFA, COEFB, COEFC, DELTA, SOL1, SOL2, ALPHA, BETA, GAMMA
      REAL(KIND=8), DIMENSION(2) :: TEST
      REAL(KIND=8), DIMENSION(4) :: NORMX, NORMY, XW, YW, ZW, BOUNDPOS
      ! REAL(KIND=8) :: XCOLL, YCOLL, ZCOLL
      REAL(KIND=8) :: VN, DX
      LOGICAL, DIMENSION(:), ALLOCATABLE :: REMOVE_PART
      REAL(KIND=8), DIMENSION(3) :: V_OLD, V_NEW
      REAL(KIND=8), DIMENSION(3) :: E, B
      REAL(KIND=8) :: V_NORM, V_PARA, V_PERP, VZ, VDUMMY, EROT, EVIB, VDOTN, WALL_TEMP
      INTEGER :: S_ID
      LOGICAL :: HASCOLLIDED
      REAL(KIND=8) :: XCOLL, YCOLL, COLLDIST, EDGE_X1, EDGE_Y1, EDGE_X2, EDGE_Y2
      INTEGER, DIMENSION(:), ALLOCATABLE :: LOCAL_BOUNDARY_COLL_COUNT, LOCAL_WALL_COLL_COUNT
      REAL(KIND=8) :: WEIGHT_RATIO
      TYPE(PARTICLE_DATA_STRUCTURE) :: NEWparticle
      INTEGER, DIMENSION(3) :: NEXTVERT
      INTEGER, DIMENSION(6) :: EDGEINDEX
      REAL(KIND=8), DIMENSION(6) :: COLLTIMES

      REAL(KIND=8) :: TOL

      TOL = 1e-15

      !E = [0.d0, 0.d0, 0.d0]
      B = [0.d0, 0.d0, 0.d0]

      NEXTVERT = [2,3,1]
      EDGEINDEX = [1,1,2,2,3,3]

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
         REMOVE_PART(IP) = .FALSE.
         IC = part_adv(IP)%IC

         ! ! Update velocity

         ! V_NEW(1) = part_adv(IP)%VX
         ! V_NEW(2) = part_adv(IP)%VY
         ! V_NEW(3) = part_adv(IP)%VZ

         ! IF (BOOL_PIC) THEN
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
            IF (GRID_TYPE == UNSTRUCTURED) THEN
               !WRITE(*,*) 'Moving particle ', IP, ' for ', DTCOLL, ' s. Position: ', part_adv(IP)%X, part_adv(IP)%Y,&
               !' velocity: ', part_adv(IP)%VX, part_adv(IP)%VY
               ! For unstructured, we only need to check the boundaries of the cell.

               ! The new C-N procedure with uniform E field.
               DO I = 1, 3
                  J = NEXTVERT(I)
                  EDGE_X1 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(IC,I),1)
                  EDGE_Y1 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(IC,I),2)
                  EDGE_X2 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(IC,J),1)
                  EDGE_Y2 = U2D_GRID%NODE_COORDS(U2D_GRID%CELL_NODES(IC,J),2)
                  IF (EDGE_X2 == EDGE_X1) THEN
                     COEFB = 0; COEFA = 1; COEFC = - EDGE_X1
                  ELSE IF (EDGE_Y2 == EDGE_Y1) THEN
                     COEFA = 0; COEFB = 1; COEFC = - EDGE_Y1
                  ELSE
                     COEFA = 1
                     COEFB = (EDGE_X1-EDGE_X2)/(EDGE_Y2-EDGE_Y1)
                     COEFC = - COEFB * EDGE_Y1 - EDGE_X1
                  END IF

                  IF (BOOL_PIC) THEN
                     CALL APPLY_E_FIELD(IP, E)
                     ALPHA = 0.5*SPECIES(part_adv(IP)%S_ID)%CHARGE*QE/SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS* &
                             (COEFA*E(1) + COEFB*E(2))
                  ELSE
                     ALPHA = 0.d0
                  END IF

                  BETA = COEFA*part_adv(IP)%VX + COEFB*part_adv(IP)%VY
                  GAMMA = COEFA*part_adv(IP)%X + COEFB*part_adv(IP)%Y + COEFC

                  IF (ALPHA == 0.d0) THEN
                     COLLTIMES(2*(I-1) + 1) = GAMMA/BETA
                     COLLTIMES(2*I) = - 1.d0
                  ELSE
                     DELTA = BETA*BETA - 4.0*ALPHA*GAMMA
                     IF (DELTA >= 0.d0) THEN
                        COLLTIMES(2*(I-1) + 1) = 0.5*(-BETA - SQRT(DELTA))/ALPHA
                        COLLTIMES(2*I) =         0.5*(-BETA + SQRT(DELTA))/ALPHA
                     ELSE
                        COLLTIMES(2*(I-1) + 1) = -1.d0
                        COLLTIMES(2*I) = - 1.d0
                     END IF
                  END IF
               END DO

               ! Find which collision happens first.
               DTCOLL = part_adv(IP)%DTRIM
               BOUNDCOLL = -1
               DO I = 1, 6
                  IF (COLLTIMES(I) > 0 .AND. COLLTIMES(I) < DTCOLL) THEN
                     CALL MOVE_PARTICLE_CN(IP, E, COLLTIMES(I))
                     J = EDGEINDEX(I)
                     IF ((part_adv(IP)%VX*U2D_GRID%EDGE_NORMAL(IC,J,1) + part_adv(IP)%VY*U2D_GRID%EDGE_NORMAL(IC,J,2)) > 0) THEN
                        DTCOLL = COLLTIMES(I)
                        BOUNDCOLL = J
                     END IF
                     CALL MOVE_PARTICLE_CN(IP, E, -COLLTIMES(I))
                  END IF
               END DO

               ! Do the advection
               IF (BOUNDCOLL .NE. -1) THEN
                  CALL MOVE_PARTICLE_CN(IP, E, DTCOLL)
                  part_adv(IP)%DTRIM = part_adv(IP)%DTRIM - DTCOLL

                  ! Move to new cell
                  IF (U2D_GRID%CELL_NEIGHBORS(IC, BOUNDCOLL) == -1) THEN
                     EDGE_PG = U2D_GRID%CELL_EDGES_PG(IC, BOUNDCOLL)
                     IF (EDGE_PG .NE. -1) THEN
                        ! Apply particle boundary condition
                        IF (GRID_BC(EDGE_PG)%PARTICLE_BC == SPECULAR) THEN
                           IF (GRID_BC(EDGE_PG)%REACT) THEN
                              CALL WALL_REACT(IP, REMOVE_PART(IP))
                           END IF
                           
                           VDOTN = part_adv(IP)%VX*U2D_GRID%EDGE_NORMAL(IC,BOUNDCOLL,1) &
                                 + part_adv(IP)%VY*U2D_GRID%EDGE_NORMAL(IC,BOUNDCOLL,2)
                           part_adv(IP)%VX = part_adv(IP)%VX - 2.*VDOTN*U2D_GRID%EDGE_NORMAL(IC,BOUNDCOLL,1)
                           part_adv(IP)%VY = part_adv(IP)%VY - 2.*VDOTN*U2D_GRID%EDGE_NORMAL(IC,BOUNDCOLL,2)
                        ELSE IF (GRID_BC(EDGE_PG)%PARTICLE_BC == DIFFUSE) THEN
                           IF (GRID_BC(EDGE_PG)%REACT) THEN
                              CALL WALL_REACT(IP, REMOVE_PART(IP))
                           END IF

                           S_ID = part_adv(IP)%S_ID
                           WALL_TEMP = GRID_BC(EDGE_PG)%WALL_TEMP
                           CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                           WALL_TEMP, WALL_TEMP, WALL_TEMP, &
                           VDUMMY, V_PARA, VZ, SPECIES(S_ID)%MOLECULAR_MASS)

                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, WALL_TEMP, EROT)
                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, WALL_TEMP, EVIB)
                                          
                           V_PERP = FLX(0.d0, WALL_TEMP, SPECIES(S_ID)%MOLECULAR_MASS)

                           part_adv(IP)%VX = - V_PERP*U2D_GRID%EDGE_NORMAL(IC,BOUNDCOLL,1) &
                                            - V_PARA*U2D_GRID%EDGE_NORMAL(IC,BOUNDCOLL,2)
                           part_adv(IP)%VY = - V_PERP*U2D_GRID%EDGE_NORMAL(IC,BOUNDCOLL,2) &
                                            + V_PARA*U2D_GRID%EDGE_NORMAL(IC,BOUNDCOLL,1)
                           part_adv(IP)%VZ = VZ
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
                     part_adv(IP)%IC = U2D_GRID%CELL_NEIGHBORS(IC, BOUNDCOLL)
                     IC = part_adv(IP)%IC
                     !WRITE(*,*) 'moved particle to cell: ', IC
                  END IF
               ELSE
                  CALL MOVE_PARTICLE_CN(IP, E, part_adv(IP)%DTRIM)
                  part_adv(IP)%DTRIM = 0.d0
               END IF
            
            ELSE
               CALL ERROR_ABORT('Error! C-N mover only implemented for unstructured grids.')
            END IF

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

      END DO ! End loop: DO IP = 1,NP_PROC


      ! Before removing the particles we record the charge displacement.
      CALL DEPOSIT_CHARGE_CN


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

   END SUBROUTINE ADVECT_CN


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MOVE_PARTICLE_CN -> Move the particles !!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE MOVE_PARTICLE_CN(IP, E, TIME)

      ! Moves particle with index IP for time TIME.
      ! For now simply rectilinear movement, will have to include Lorentz's force

      IMPLICIT NONE

      INTEGER, INTENT(IN)      :: IP
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: E
      REAL(KIND=8), INTENT(IN) :: TIME
      REAL(KIND=8), DIMENSION(3) :: ACC

      ACC = E*SPECIES(part_adv(IP)%S_ID)%CHARGE*QE/SPECIES(part_adv(IP)%S_ID)%MOLECULAR_MASS

      part_adv(IP)%VX = part_adv(IP)%VX + ACC(1) * TIME
      part_adv(IP)%VY = part_adv(IP)%VY + ACC(2) * TIME
      part_adv(IP)%VZ = part_adv(IP)%VZ + ACC(3) * TIME

      part_adv(IP)%X = part_adv(IP)%X + part_adv(IP)%VX*TIME + 0.5*ACC(1)*TIME*TIME
      part_adv(IP)%Y = part_adv(IP)%Y + part_adv(IP)%VY*TIME + 0.5*ACC(2)*TIME*TIME
      part_adv(IP)%Z = part_adv(IP)%Z + part_adv(IP)%VZ*TIME + 0.5*ACC(3)*TIME*TIME

   END SUBROUTINE MOVE_PARTICLE_CN

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


   SUBROUTINE WALL_REACT(IP, REMOVE)
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IP
      LOGICAL, INTENT(OUT) :: REMOVE
      INTEGER :: JS, JR, JP
      REAL(KIND=8) :: PROB_SCALE, VEL_SCALE

      JS = particles(IP)%S_ID
      PROB_SCALE = 1.
      REMOVE = .FALSE.
      DO JR = 1, N_WALL_REACTIONS
         IF (WALL_REACTIONS(JR)%R_SP_ID == JS) THEN
            IF ( rf() .LE. WALL_REACTIONS(JR)%PROB/PROB_SCALE ) THEN
               
               IF (WALL_REACTIONS(JR)%N_PROD == 0) THEN
                  REMOVE = .TRUE.
                  particles(IP)%DTRIM = 0.
               ELSE IF (WALL_REACTIONS(JR)%N_PROD == 1) THEN
                  JP = WALL_REACTIONS(JR)%P1_SP_ID
                  particles(IP)%S_ID = JP
                  VEL_SCALE = SPECIES(JS)%MOLECULAR_MASS/SPECIES(JP)%MOLECULAR_MASS
                  particles(IP)%VX = particles(IP)%VX*VEL_SCALE
                  particles(IP)%VY = particles(IP)%VY*VEL_SCALE
                  particles(IP)%VZ = particles(IP)%VZ*VEL_SCALE
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

      CALL VecGetArrayReadF90(X_SEQ,PHI_FIELD_TEMP,ierr)
      IF (ALLOCATED(PHI_FIELD)) DEALLOCATE(PHI_FIELD)
      ALLOCATE(PHI_FIELD, SOURCE = PHI_FIELD_TEMP)
      CALL VecRestoreArrayReadF90(X_SEQ,PHI_FIELD_TEMP,ierr)

      PHIBAR_FIELD = PHI_FIELD

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

   END SUBROUTINE DEPOSIT_CHARGE

END MODULE fully_implicit
