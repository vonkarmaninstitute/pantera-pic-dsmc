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

MODULE timecycle

   USE global
   USE mpi_common
   USE screen
   USE tools
   USE collisions
   USE postprocess
   USE fields
   USE fully_implicit
   USE washboard

   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE TIME_LOOP                                     !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE TIME_LOOP

      IMPLICIT NONE
   
      INTEGER :: NP_TOT, NCOLL_TOT, NREAC_TOT
      REAL(KIND=8) :: FIELD_POWER_TOT
      REAL(KIND=8) :: CURRENT_TIME, CURRENT_CPU_TIME, EST_TIME
      INTEGER :: EST_TIME_H, EST_TIME_M

      CHARACTER(len=512) :: stringTMP

      ! Init variables
      NP_TOT = 0
      CALL INIT_POSTPROCESS
      CALL INIT_BOUNDARY_POSTPROCESS
      ! CALL BOUNDARY_RESET

      ALLOCATE(FIELD_POWER_AVG(FIELD_POWER_NUMAVG))
      FIELD_POWER_AVG = FIELD_POWER_TARGET

      ! ########### Compute poisson ##########################################

      IF (PIC_TYPE .NE. NONE) THEN

         CALL SET_WALL_POTENTIAL

         CALL DEPOSIT_CHARGE(particles)
         IF (PIC_TYPE == HYBRID) THEN           
            CALL SOLVE_BOLTZMANN
         ELSE
            CALL SETUP_POISSON
            CALL SOLVE_POISSON
         END IF
         PHIBAR_FIELD = PHI_FIELD
         CALL COMPUTE_E_FIELD

         !IF (PIC_TYPE == SEMIIMPLICIT) THEN
         !   CALL DEPOSIT_CURRENT
         !   CALL ASSEMBLE_AMPERE
         !   CALL SOLVE_AMPERE
         !ELSE IF (PIC_TYPE == FULLYIMPLICIT) THEN
         !   !CALL SOLVE_POISSON_FULLY_IMPLICIT
         !   ! We get here no problem.
         !END IF
      END IF
      
      ! ########### Dump particles and flowfield after the initial seeding ##########################################
      IF ((tID .GE. DUMP_PART_START) .AND. (tID .NE. RESTART_TIMESTEP)) THEN
         IF (MOD(tID-DUMP_PART_START, DUMP_PART_EVERY) .EQ. 0) CALL DUMP_PARTICLES_FILE(tID)
      END IF

      IF ((tID .GE. DUMP_GRID_START) .AND. (tID .NE. RESTART_TIMESTEP)) THEN
         ! If we are in the grid save timestep, average, then dump the cumulated averages
         IF (MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY*DUMP_GRID_N_AVG) .EQ. 0) THEN
            CALL GRID_AVG
            CALL GRID_SAVE
            CALL GRID_RESET
         END IF
      END IF

      IF ((tID .GE. DUMP_BOUND_START) .AND. (tID .NE. RESTART_TIMESTEP)) THEN
         ! If we are in the grid save timestep, average, then dump the cumulated averages
         IF (MOD(tID-DUMP_BOUND_START, DUMP_BOUND_AVG_EVERY*DUMP_BOUND_N_AVG) .EQ. 0) THEN
            CALL BOUNDARY_GATHER
            CALL BOUNDARY_SAVE
            CALL BOUNDARY_RESET
         END IF
      END IF

      ! ########### Perform the conservation checks ###################################
      IF (PERFORM_CHECKS .AND. MOD(tID, CHECKS_EVERY) .EQ. 0) CALL CHECKS


      ! ########### Start the time loop #################################
      tID = tID + 1
      CALL CPU_TIME(START_CPU_TIME)
      DO WHILE (tID .LE. NT)


         CALL MPI_REDUCE(FIELD_POWER, FIELD_POWER_TOT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(FIELD_POWER_TOT, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         !FIELD_POWER_AVG(MOD(tID, FIELD_POWER_NUMAVG) + 1) = FIELD_POWER_TOT
         !IF (MOD(tID, FIELD_POWER_NUMAVG) .EQ. 0) THEN
         !   WRITE(*,*) 'Adjusting coil current.'
         !   COIL_CURRENT = COIL_CURRENT*SQRT(FIELD_POWER_TARGET/(SUM(FIELD_POWER_AVG)/DBLE(FIELD_POWER_NUMAVG)))
         !END IF



         ! IF (tID == 40001) THEN
         !    DO IP = 1, NP_PROC
         !       IF (rf() < 1.d-4) THEN
         !          particles(IP)%DUMP_TRAJ = .TRUE.
         !       END IF
         !    END DO
         ! END IF

         ! ########### Print simulation info #######################################

         CURRENT_TIME = tID*DT
         CALL CPU_TIME(CURRENT_CPU_TIME)
         CURRENT_CPU_TIME = CURRENT_CPU_TIME - START_CPU_TIME 
         EST_TIME = DBLE(NT-tID)/DBLE(tID-RESTART_TIMESTEP)*CURRENT_CPU_TIME
         EST_TIME_H = INT(EST_TIME/3600.d0)
         EST_TIME_M = INT((EST_TIME - EST_TIME_H*3600.d0)/60.d0)

         IF (MOD(tID, STATS_EVERY) .EQ. 0) THEN
            CALL MPI_REDUCE(NP_PROC, NP_TOT, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(TIMESTEP_COLL, NCOLL_TOT, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(TIMESTEP_REAC, NREAC_TOT, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(FIELD_POWER, FIELD_POWER_TOT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            

            WRITE(stringTMP, '(A13,I8,A4,I8,A9,ES14.3,A17,F10.1,A27,I5,A5,I2,A4,A24,I10,A25,I10,A24,I10)') &
                           '   Timestep: ', tID, ' of ', NT, &
                           ' - time: ', CURRENT_TIME, ' [s] - CPU time: ', CURRENT_CPU_TIME, &
                           ' [s] - est. time required: ', EST_TIME_H, ' [h] ' , EST_TIME_M, ' [m]',&
                           ' - number of particles: ', NP_TOT, &
                           ' - number of collisions: ', NCOLL_TOT, &
                           ' - number of reactions: ', NREAC_TOT

            ! Use this tho have the coil current and field power output to console.
            ! WRITE(stringTMP, '(A13,I8,A4,I8,A9,ES14.3,A17,F10.1,A27,I5,A5,I2,A4,A24,I10, &
            ! A25,I10,A24,I10,A23,ES14.3,A4,A17,ES14.3,A4)') &
            !                '   Timestep: ', tID, ' of ', NT, &
            !                ' - time: ', CURRENT_TIME, ' [s] - CPU time: ', CURRENT_CPU_TIME, &
            !                ' [s] - est. time required: ', EST_TIME_H, ' [h] ' , EST_TIME_M, ' [m]',&
            !                ' - number of particles: ', NP_TOT, &
            !                ' - number of collisions: ', NCOLL_TOT, &
            !                ' - number of reactions: ', NREAC_TOT, &
            !                ' - rf power deposited: ', FIELD_POWER_TOT, ' [W]', &
            !                ' - coil current: ', COIL_CURRENT, ' [A]'

            CALL ONLYMASTERPRINT1(PROC_ID, TRIM(stringTMP))

         END IF


         IF (tID .GT. 0 .AND. MOD(tID, TIMING_STATS_EVERY) .EQ. 0) THEN
            CALL TIMER_SUMMARY
         END IF
         
         !FLUSH(101)

         ! ########### Inject particles from boundaries/injection sources ##########

         CALL BOUNDARIES_INJECT
         CALL LINE_SOURCE_INJECT
         CALL BOUNDARIES_EMIT
         CALL VOLUME_INJECT

         !CALL FIXED_IONIZATION

         ! ########### Perform load balancing ###################################
         IF (LOAD_BALANCE) THEN
            IF (tID .GT. 0 .AND. MOD(tID, LOAD_BALANCE_EVERY) .EQ. 0) THEN
               CALL ASSIGN_CELLS_TO_PROCS
            END IF
         END IF

         ! ########### Exchange particles among processes ##########################

         CALL TIMER_START(5)
         CALL EXCHANGE
         CALL TIMER_STOP(5)

         ! ########### Advect particles and update field ############################################

         IF (PIC_TYPE == EXPLICIT) THEN
            CALL SET_WALL_POTENTIAL
            CALL DEPOSIT_CHARGE(particles)

            CALL TIMER_START(2)
            CALL SOLVE_POISSON
            CALL COMPUTE_E_FIELD
            CALL TIMER_STOP(2)

            CALL TIMER_START(3)
            CALL ADVECT
            CALL TIMER_STOP(3)

         ELSE IF (PIC_TYPE == EXPLICITLIMITED) THEN   
           
            CALL ASSEMBLE_POISSON
            CALL SET_WALL_POTENTIAL
            CALL DEPOSIT_CHARGE(particles)

            CALL TIMER_START(2)
            CALL SETUP_POISSON
            CALL SOLVE_POISSON
            CALL COMPUTE_E_FIELD
            CALL TIMER_STOP(2)

            CALL TIMER_START(3)
            CALL ADVECT
            CALL TIMER_STOP(3)
         
         ELSE IF (PIC_TYPE == SEMIIMPLICIT) THEN
         
            CALL DEPOSIT_CURRENT

            CALL TIMER_START(2)
            CALL ASSEMBLE_AMPERE
            CALL SOLVE_AMPERE
            CALL COMPUTE_E_FIELD
            CALL TIMER_STOP(2)

            CALL TIMER_START(3)
            CALL ADVECT
            CALL TIMER_STOP(3)

         ELSE IF (PIC_TYPE == FULLYIMPLICIT) THEN

            CALL TIMER_START(2)
            CALL SET_WALL_POTENTIAL
            CALL SOLVE_POISSON_FULLY_IMPLICIT
            CALL COMPUTE_E_FIELD
            CALL TIMER_STOP(2)
            
            CALL TIMER_START(3)
            CALL ADVECT_CN_B(particles, .TRUE., .FALSE., Jmat)
            CALL TIMER_STOP(3)

         ELSE IF (PIC_TYPE == HYBRID) THEN

            CALL SET_WALL_POTENTIAL
            CALL DEPOSIT_CHARGE(particles)

            CALL TIMER_START(2)
            CALL SOLVE_BOLTZMANN
            CALL COMPUTE_E_FIELD
            CALL TIMER_STOP(2)

            CALL TIMER_START(3)
            CALL ADVECT
            CALL FLUID_ELECTRONS_SURFACE_CHARGING
            CALL TIMER_STOP(3)
            
         ELSE

            CALL TIMER_START(3)
            CALL ADVECT
            CALL TIMER_STOP(3)
            
         END IF

         ! ########### Exchange particles among processes ##########################

         CALL TIMER_START(5)
         CALL EXCHANGE
         CALL TIMER_STOP(5)
         

         ! ########### Perform collisions ##########################################

         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         CALL TIMER_START(6)
         IF (COLLISION_TYPE == MCC)  CALL MCC_COLLISIONS
         IF (COLLISION_TYPE == MCC_VAHEDI)  CALL MCC_COLLISIONS_VAHEDI


         IF (COLLISION_TYPE == DSMC .OR. COLLISION_TYPE == DSMC_VAHEDI) THEN
            CALL DSMC_COLLISIONS
         END IF

         IF (COLLISION_TYPE == BGK) CALL BGK_COLLISIONS

         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         CALL TIMER_STOP(6)

         IF (BOOL_THERMAL_BATH) CALL THERMAL_BATH




         CALL TIMER_START(4)
         ! ########### Dump particles ##############################################
         IF (tID .GE. DUMP_PART_START) THEN
            IF (MOD(tID-DUMP_PART_START, DUMP_PART_EVERY) .EQ. 0) CALL DUMP_PARTICLES_FILE(tID)
         END IF

         ! ########### Dump flowfield ##############################################
         IF (tID .GT. DUMP_GRID_START) THEN
            ! If we are in the grid save timestep, average, then dump the cumulated averages
            IF (MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY*DUMP_GRID_N_AVG) .EQ. 0) THEN
               CALL GRID_AVG
               CALL GRID_SAVE
               CALL GRID_RESET
            ! If we are just in a grid average timestep, compute the grid average
            ELSE IF (MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY) .EQ. 0) THEN
               CALL GRID_AVG
            END IF
         END IF

         ! ########### Dump boundary ##############################################
         IF (tID .GT. DUMP_BOUND_START) THEN
            ! If we are in the grid save timestep, average, then dump the cumulated averages
            IF (MOD(tID-DUMP_BOUND_START, DUMP_BOUND_AVG_EVERY*DUMP_BOUND_N_AVG) .EQ. 0) THEN
               CALL BOUNDARY_GATHER
               CALL BOUNDARY_SAVE
               CALL BOUNDARY_RESET
            ! If we are just in a grid average timestep, compute the grid average
            ELSE IF (MOD(tID-DUMP_BOUND_START, DUMP_BOUND_AVG_EVERY) .EQ. 0) THEN
               CALL BOUNDARY_GATHER
            END IF
         END IF

         ! ########### Dump particle fluxes #######################################
         IF (BOOL_DUMP_FLUXES) CALL DUMP_FLUXES_FILE(tID)

         ! ########### Dump individual particle ###################################
         CALL DUMP_TRAJECTORY_FILE(tID)
         CALL TIMER_STOP(4)


         ! ########### Perform the conservation checks ###################################
         IF (PERFORM_CHECKS .AND. MOD(tID, CHECKS_EVERY) .EQ. 0) CALL CHECKS

         IF (REMOVE_MIX .NE. -1) CALL REMOVE_PARTICLES_IN_MIXTURE(REMOVE_MIX)


         !IF (PERFORM_CHECKS .AND. MOD(tID, CHECKS_EVERY) .EQ. 0) THEN
         !   CALL ONLYMASTERPRINT1(PROC_ID, '---> Checking if particles are in the correct cells.')
         !   CALL REASSIGN_PARTICLES_TO_CELLS_2D
         !END IF

         ! ~~~~~ Hmm that's it! ~~~~~

         
         tID = tID + 1

      END DO

   END SUBROUTINE TIME_LOOP
 


   SUBROUTINE BOUNDARIES_EMIT
  
      IMPLICIT NONE
   
      INTEGER      :: IP, IC, IS, NFS, ITASK
      REAL(KIND=8) :: DTFRAC, Vdummy, V_NORM, V_TANG1, V_TANG2, BETA, BETA_E, X1, X2, Y1, Y2, R, P, Q, S, T
      REAL(KIND=8) :: X, Y, Z, VX, VY, VZ, EROT, EVIB 
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
      REAL(KIND=8), DIMENSION(3) :: FACE_NORMAL, FACE_TANG1, FACE_TANG2, V1, V2, V3
   
      INTEGER :: S_ID, ELECTRON_S_ID
      REAL(KIND=8) :: M

      TYPE(EMIT_TASK_DATA_STRUCTURE) :: EMIT_TASK

      IF (COLOCATED_ELECTRONS) ELECTRON_S_ID = SPECIES_NAME_TO_ID('e')

      DO ITASK = 1, N_EMIT_TASKS
         EMIT_TASK = EMIT_TASKS(ITASK)

         IC = EMIT_TASK%IC

         IF (DIMS == 1) THEN
            X1 = U1D_GRID%NODE_COORDS(1, EMIT_TASK%IV1)
            X2 = U1D_GRID%NODE_COORDS(1, EMIT_TASK%IV2)

            FACE_NORMAL = U1D_GRID%EDGE_NORMAL(:,EMIT_TASK%IFACE,IC)
            FACE_TANG1 = [0.d0, FACE_NORMAL(1), 0.d0]
            FACE_TANG2 = [0.d0, 0.d0, 1.d0]
         ELSE IF (DIMS == 2) THEN
            X1 = U2D_GRID%NODE_COORDS(1, EMIT_TASK%IV1)
            Y1 = U2D_GRID%NODE_COORDS(2, EMIT_TASK%IV1)
            X2 = U2D_GRID%NODE_COORDS(1, EMIT_TASK%IV2)
            Y2 = U2D_GRID%NODE_COORDS(2, EMIT_TASK%IV2)

            FACE_NORMAL = U2D_GRID%EDGE_NORMAL(:,EMIT_TASK%IFACE,IC)
            FACE_TANG1 = [-FACE_NORMAL(2), FACE_NORMAL(1), 0.d0]
            FACE_TANG2 = [0.d0, 0.d0, 1.d0]
         ELSE IF (DIMS == 3) THEN
            FACE_NORMAL = U3D_GRID%FACE_NORMAL(:,EMIT_TASK%IFACE,IC)
            FACE_TANG1 = U3D_GRID%FACE_TANG1(:,EMIT_TASK%IFACE,IC)
            FACE_TANG2 = U3D_GRID%FACE_TANG2(:,EMIT_TASK%IFACE,IC)

            V1 = U3D_GRID%NODE_COORDS(:,U3D_GRID%FACE_NODES(1,EMIT_TASK%IFACE,IC))
            V2 = U3D_GRID%NODE_COORDS(:,U3D_GRID%FACE_NODES(2,EMIT_TASK%IFACE,IC))
            V3 = U3D_GRID%NODE_COORDS(:,U3D_GRID%FACE_NODES(3,EMIT_TASK%IFACE,IC))
         END IF

         DO IS = 1, MIXTURES(EMIT_TASK%MIX_ID)%N_COMPONENTS ! Loop on mixture components

            S_ID = MIXTURES(EMIT_TASK%MIX_ID)%COMPONENTS(IS)%ID

            M = SPECIES(S_ID)%MOLECULAR_MASS
            NFS = FLOOR(EMIT_TASK%NFS(IS))
            IF (EMIT_TASK%NFS(IS)-REAL(NFS, KIND=8) .GE. rf()) THEN ! Same as SPARTA's perspeciess
               NFS = NFS + 1
            END IF

            IF (EMIT_TASK%TTRA == 0) THEN
               BETA = 0
            ELSE
               BETA = 1./SQRT(2.*KB/M*EMIT_TASK%TTRA)
               ! IF (BOOL_KAPPA_DISTRIBUTION)  BETA = 1./SQRT(2.*KB/M*EMIT_TASK%TTRA*(KAPPA_C-3./2.))
            END IF

            IF (COLOCATED_ELECTRONS) BETA_E = 1./SQRT(2.*KB/SPECIES(ELECTRON_S_ID)%MOLECULAR_MASS*COLOCATED_ELECTRONS_TTRA)

            DO IP = 1, NFS ! Loop on particles to be injected

               CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                           EMIT_TASK%TTRA, EMIT_TASK%TTRA, EMIT_TASK%TTRA, &
                           Vdummy, V_TANG1, V_TANG2, M)

               CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, EMIT_TASK%TROT, EROT)
               CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, EMIT_TASK%TVIB, EVIB)

               IF (EMIT_TASK%TTRA == 0) THEN
                  V_NORM = 0
               ELSE
                  V_NORM = FLX(EMIT_TASK%U_NORM*BETA, EMIT_TASK%TTRA, M)
               END IF

               VX = EMIT_TASK%UX &
                  - V_NORM*FACE_NORMAL(1) &
                  - V_TANG1*FACE_TANG1(1) &
                  - V_TANG2*FACE_TANG2(1)
               VY = EMIT_TASK%UY &
                  - V_NORM*FACE_NORMAL(2) &
                  - V_TANG1*FACE_TANG1(2) &
                  - V_TANG2*FACE_TANG2(2)
               VZ = EMIT_TASK%UZ &
                  - V_NORM*FACE_NORMAL(3) &
                  - V_TANG1*FACE_TANG1(3) &
                  - V_TANG2*FACE_TANG2(3)

               DTFRAC = rf()*DT


               IF (DIMS == 1) THEN
                  R = rf()

                  X = X1 + R*(X2-X1)
                  Y = YMIN + (YMAX-YMIN)*rf()
                  Z = ZMIN + (ZMAX-ZMIN)*rf()
               ELSE IF (DIMS == 2) THEN
                  R = rf()
                  IF (AXI) THEN
                     IF (Y1 == Y2) THEN
                        Y = Y1
                        X = X1 + R*(X2-X1)
                        Z = 0
                     ELSE IF (Y1 < Y2) THEN
                        Y = SQRT(Y1*Y1 + R*(Y2*Y2 - Y1*Y1))
                        X = X1 + (Y-Y1)/(Y2-Y1)*(X2-X1)
                        Z = 0
                     ELSE
                        Y = SQRT(Y2*Y2 + R*(Y1*Y1 - Y2*Y2))
                        X = X2 + (Y-Y2)/(Y1-Y2)*(X1-X2)
                        Z = 0
                     END IF
                  ELSE
                     X = X1 + R*(X2-X1)
                     Y = Y1 + R*(Y2-Y1)
                     Z = ZMIN + (ZMAX-ZMIN)*rf()
                  END IF
               ELSE IF (DIMS == 3) THEN
                  S = rf()
                  T = rf()
                  
                  IF (S > 1-T) THEN
                     Q = 1-T
                     P = 1-S
                  ELSE
                     Q = T
                     P = S
                  END IF

                  X = V1(1) + (V2(1)-V1(1))*P + (V3(1)-V1(1))*Q
                  Y = V1(2) + (V2(2)-V1(2))*P + (V3(2)-V1(2))*Q
                  Z = V1(3) + (V2(3)-V1(3))*P + (V3(3)-V1(3))*Q
               END IF

               ! Init a particle object and assign it to the local vector of particles
               CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTFRAC,  particleNOW)
               CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles)

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Modification to inject electrons colocated with ions.
               IF (COLOCATED_ELECTRONS .AND. SPECIES(S_ID)%CHARGE == 1) THEN
                        
                  CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                  COLOCATED_ELECTRONS_TTRA, COLOCATED_ELECTRONS_TTRA, COLOCATED_ELECTRONS_TTRA, &
                  Vdummy, V_TANG1, V_TANG2, SPECIES(ELECTRON_S_ID)%MOLECULAR_MASS)
   
                  V_NORM = FLX(EMIT_TASK%U_NORM*BETA_E, COLOCATED_ELECTRONS_TTRA, SPECIES(ELECTRON_S_ID)%MOLECULAR_MASS)
   
                  VX = EMIT_TASK%UX &
                     - V_NORM*FACE_NORMAL(1) &
                     - V_TANG1*FACE_TANG1(1) &
                     - V_TANG2*FACE_TANG2(1)
                  VY = EMIT_TASK%UY &
                     - V_NORM*FACE_NORMAL(2) &
                     - V_TANG1*FACE_TANG1(2) &
                     - V_TANG2*FACE_TANG2(2)
                  VZ = EMIT_TASK%UZ &
                     - V_NORM*FACE_NORMAL(3) &
                     - V_TANG1*FACE_TANG1(3) &
                     - V_TANG2*FACE_TANG2(3)

                  CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,0.d0,0.d0,ELECTRON_S_ID,IC,DTFRAC, particleNOW) ! Save in particle
                  CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
               END IF
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            END DO

         END DO
      END DO

   END SUBROUTINE BOUNDARIES_EMIT





   SUBROUTINE VOLUME_INJECT
      
      IMPLICIT NONE

      INTEGER            :: IP, NP_INIT, IC
      REAL(KIND=8)       :: Xp, Yp, Zp, VXp, VYp, VZp, EROT, EVIB, DOMAIN_VOLUME, VOL
      INTEGER            :: CID

      REAL(KIND=8), DIMENSION(3) :: V1, V2, V3, V4
      REAL(KIND=8)       :: P, Q, R, S, T, U

      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
      REAL(KIND=8)  :: M
      INTEGER       :: S_ID, I, J, N_INJECT, ITASK
     
      ! Compute number of particles to be seeded
      !NP_INIT = NINT(NRHO_INIT/FNUM*(XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN))

      !CALL ONLYMASTERPRINT2(PROC_ID, '   Particles to be seeded', REAL(NP_INIT, KIND=8))
 
      ! Compute number of particles to be seeded by every process
      !NPPP_INIT = INT(NP_INIT/N_MPI_THREADS) 

      DO ITASK = 1, N_VOLUME_INJECT_TASKS
         ! Create particles in the domain
         DO I = 1, MIXTURES(VOLUME_INJECT_TASKS(ITASK)%MIX_ID)%N_COMPONENTS
            
            S_ID = MIXTURES(VOLUME_INJECT_TASKS(ITASK)%MIX_ID)%COMPONENTS(i)%ID
            IF (GRID_TYPE == UNSTRUCTURED) THEN

               DO IC = 1, NCELLS
                  ! Compute number of particles of this species per process to be created in this cell.
                  IF (DIMS == 1) THEN
                     VOL = U1D_GRID%CELL_VOLUMES(IC)   
                  ELSE IF (DIMS == 2) THEN
                     VOL = U2D_GRID%CELL_VOLUMES(IC)
                  ELSE IF (DIMS == 3) THEN
                     VOL = U3D_GRID%CELL_VOLUMES(IC)
                  END IF
                  NP_INIT = RANDINT(DT*VOLUME_INJECT_TASKS(ITASK)%NRHODOT/(FNUM*SPECIES(S_ID)%SPWT)*VOL* &
                              MIXTURES(VOLUME_INJECT_TASKS(ITASK)%MIX_ID)%COMPONENTS(i)%MOLFRAC/N_MPI_THREADS)
                  IF (NP_INIT == 0) CYCLE

                  IF (DIMS == 1) THEN
                     V1 = U1D_GRID%NODE_COORDS(:,U2D_GRID%CELL_NODES(1,IC))
                     V2 = U1D_GRID%NODE_COORDS(:,U2D_GRID%CELL_NODES(2,IC))
                     V3 = 0
                     V4 = 0
                  ELSE IF (DIMS == 2) THEN
                     V1 = U2D_GRID%NODE_COORDS(:,U2D_GRID%CELL_NODES(1,IC))
                     V2 = U2D_GRID%NODE_COORDS(:,U2D_GRID%CELL_NODES(2,IC))
                     V3 = U2D_GRID%NODE_COORDS(:,U2D_GRID%CELL_NODES(3,IC))
                     V4 = 0
                  ELSE IF (DIMS == 3) THEN
                     V1 = U3D_GRID%NODE_COORDS(:,U3D_GRID%CELL_NODES(1,IC))
                     V2 = U3D_GRID%NODE_COORDS(:,U3D_GRID%CELL_NODES(2,IC))
                     V3 = U3D_GRID%NODE_COORDS(:,U3D_GRID%CELL_NODES(3,IC))
                     V4 = U3D_GRID%NODE_COORDS(:,U3D_GRID%CELL_NODES(4,IC))
                  END IF

                  DO IP = 1, NP_INIT

                     ! Create particle position randomly in the cell
                     S = rf()
                     T = rf()
                     U = rf()
                     
                     IF (S+T > 1) THEN
                        S = 1-S; T = 1-T
                     END IF

                     IF (DIMS == 1) THEN
                        XP = V1(1) + (V2(1)-V1(1))*S
                        YP = YMIN + (YMAX-YMIN)*T
                        ZP = ZMIN + (ZMAX-ZMIN)*U
                     ELSE IF (DIMS == 2) THEN
                        XP = V1(1) + (V2(1)-V1(1))*S + (V3(1)-V1(1))*T
                        YP = V1(2) + (V2(2)-V1(2))*S + (V3(2)-V1(2))*T
                        IF (AXI) THEN
                           ZP = 0.d0
                        ELSE
                           ZP = ZMIN + (ZMAX-ZMIN)*U
                        END IF
                     ELSE IF (DIMS == 3) THEN
                        ! http://vcg.isti.cnr.it/publications/papers/rndtetra_a.pdf

                        IF (S+T+U <= 1) THEN
                           P = S; Q = T; R = U
                        ELSE
                           IF (T+U > 1) THEN
                              P = S; Q = 1-U; R = 1-S-T
                           ELSE
                              P = 1-T-U; Q = T; R = S+T+U-1
                           END IF
                        END IF

                        XP = V1(1) + (V2(1)-V1(1))*P + (V3(1)-V1(1))*Q + (V4(1)-V1(1))*R
                        YP = V1(2) + (V2(2)-V1(2))*P + (V3(2)-V1(2))*Q + (V4(2)-V1(2))*R
                        ZP = V1(3) + (V2(3)-V1(3))*P + (V3(3)-V1(3))*Q + (V4(3)-V1(3))*R
                     END IF

                     ! Assign velocity and energy following a Boltzmann distribution
                     M = SPECIES(S_ID)%MOLECULAR_MASS
                     CALL MAXWELL(VOLUME_INJECT_TASKS(ITASK)%UX, &
                                 VOLUME_INJECT_TASKS(ITASK)%UY, &
                                 VOLUME_INJECT_TASKS(ITASK)%UZ, &
                                 VOLUME_INJECT_TASKS(ITASK)%TTRAX, &
                                 VOLUME_INJECT_TASKS(ITASK)%TTRAY, &
                                 VOLUME_INJECT_TASKS(ITASK)%TTRAZ, &
                                 VXP, VYP, VZP, M)

                     CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, VOLUME_INJECT_TASKS(ITASK)%TROT, EROT)
                     CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, VOLUME_INJECT_TASKS(ITASK)%TVIB, EVIB)

                     CALL INIT_PARTICLE(XP,YP,ZP,VXP,VYP,VZP,EROT,EVIB,S_ID,IC,DT, particleNOW) ! Save in particle
                     CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
                  END DO
               END DO
               
            ELSE ! Structured grid
               ! Compute number of particles of this species per process to be created.
               IF (AXI) THEN
                  DOMAIN_VOLUME = 0.5*(XMAX-XMIN)*(YMAX**2-YMIN**2)*(ZMAX-ZMIN)
               ELSE
                  DOMAIN_VOLUME = (XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN)
               END IF

               NP_INIT = RANDINT(DT*VOLUME_INJECT_TASKS(ITASK)%NRHODOT/(FNUM*SPECIES(S_ID)%SPWT)*DOMAIN_VOLUME* &
                           MIXTURES(VOLUME_INJECT_TASKS(ITASK)%MIX_ID)%COMPONENTS(i)%MOLFRAC/N_MPI_THREADS)
               IF (NP_INIT == 0) CYCLE
               DO IP = 1, NP_INIT

                  ! Create particle position randomly in the domain
                  XP = XMIN + (XMAX-XMIN)*rf()

                  IF (AXI .AND. (.NOT. BOOL_RADIAL_WEIGHTING)) THEN
                     YP = SQRT(YMIN*YMIN + rf()*(YMAX*YMAX - YMIN*YMIN))
                     ZP = 0
                  ELSE
                     YP = YMIN + (YMAX-YMIN)*rf()
                     ZP = ZMIN + (ZMAX-ZMIN)*rf()
                  END IF


               
                  ! Assign velocity and energy following a Boltzmann distribution
                  M = SPECIES(S_ID)%MOLECULAR_MASS
                  CALL MAXWELL(VOLUME_INJECT_TASKS(ITASK)%UX, &
                               VOLUME_INJECT_TASKS(ITASK)%UY, &
                               VOLUME_INJECT_TASKS(ITASK)%UZ, &
                               VOLUME_INJECT_TASKS(ITASK)%TTRAX, &
                               VOLUME_INJECT_TASKS(ITASK)%TTRAY, &
                               VOLUME_INJECT_TASKS(ITASK)%TTRAZ, &
                               VXP, VYP, VZP, M)

                  CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, VOLUME_INJECT_TASKS(ITASK)%TROT, EROT)
                  CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, VOLUME_INJECT_TASKS(ITASK)%TVIB, EVIB)

                  CALL CELL_FROM_POSITION(XP,YP,  CID) ! Find cell containing particle

                  CALL INIT_PARTICLE(XP,YP,ZP,VXP,VYP,VZP,EROT,EVIB,S_ID,CID,DT, particleNOW) ! Save in particle
                  CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
               END DO
            END IF
         END DO
      END DO





      IF (BOOL_INJECT_FROM_FILE) THEN
         DO I = 1, NP_INJECT_PROC
            N_INJECT = INT(INJECT_PROBABILITY + rf())
            particleNOW = part_inject(I)
            DO J = 1, N_INJECT
               particleNOW%DTRIM = rf()*DT
               CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
            END DO
         END DO
      END IF

      ! ~~~~~~ At this point, exchange particles among processes ~~~~~~

   END SUBROUTINE VOLUME_INJECT




   ! Should not be used in unstructured.
   SUBROUTINE LINE_SOURCE_INJECT
  
      IMPLICIT NONE
   
      INTEGER      :: IP, IC, IS, NFS, ILINE, DUMP_COUNTER
      REAL(KIND=8) :: DTFRAC, Vdummy, V_NORM, V_PERP, X1, X2, Y1, Y2, R
      REAL(KIND=8) :: X, Y, Z, VX, VY, VZ, EROT, EVIB 
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
   
      INTEGER :: S_ID
      REAL(KIND=8) :: M

      INTEGER, DIMENSION(:), ALLOCATABLE :: LOCAL_LINE_EMIT_COUNT

      ALLOCATE(LOCAL_LINE_EMIT_COUNT(N_LINESOURCES*N_SPECIES))
      LOCAL_LINE_EMIT_COUNT = 0

   
      ! Linesource
      DO ILINE = 1, N_LINESOURCES
         
         DO IS = 1, MIXTURES(LINESOURCES(ILINE)%MIX_ID)%N_COMPONENTS ! Loop on mixture components

            S_ID = MIXTURES(LINESOURCES(ILINE)%MIX_ID)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            NFS = FLOOR(LINESOURCES(ILINE)%nfs(IS))
            IF (LINESOURCES(ILINE)%nfs(IS)-REAL(NFS, KIND=8) .GE. rf()) THEN ! Same as SPARTA's perspeciess
               NFS = NFS + 1
            END IF
            DUMP_COUNTER = 0
            DO IP = 1, NFS ! Loop on particles to be injected
               
               CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                            LINESOURCES(ILINE)%TTRA, LINESOURCES(ILINE)%TTRA, LINESOURCES(ILINE)%TTRA, &
                            Vdummy, V_PERP, VZ, M)

               CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, LINESOURCES(ILINE)%TROT, EROT)
               CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, LINESOURCES(ILINE)%TVIB, EVIB)
                            
               V_NORM = FLX(LINESOURCES(ILINE)%S_NORM, LINESOURCES(ILINE)%TTRA, M) !!AAAAAAAAAAA

               VX = V_NORM*LINESOURCES(ILINE)%NORMX - V_PERP*LINESOURCES(ILINE)%NORMY + LINESOURCES(ILINE)%UX
               VY = V_PERP*LINESOURCES(ILINE)%NORMX + V_NORM*LINESOURCES(ILINE)%NORMY + LINESOURCES(ILINE)%UY
               VZ = VZ + LINESOURCES(ILINE)%UZ

               DTFRAC = rf()*DT

               R = rf()
               IF (DIMS == 2 .AND. AXI .AND. (.NOT. BOOL_RADIAL_WEIGHTING)) THEN
                  ! Correct but not really robust when Y1 \approx Y2
                  IF (LINESOURCES(ILINE)%Y1 == LINESOURCES(ILINE)%Y2) THEN
                     Y = LINESOURCES(ILINE)%Y1
                     X = LINESOURCES(ILINE)%X1 + R*(LINESOURCES(ILINE)%X2-LINESOURCES(ILINE)%X1)
                     Z = 0
                  ELSE
                     IF (LINESOURCES(ILINE)%Y1 < LINESOURCES(ILINE)%Y2) THEN
                        Y1 = LINESOURCES(ILINE)%Y1
                        Y2 = LINESOURCES(ILINE)%Y2
                        X1 = LINESOURCES(ILINE)%X1
                        X2 = LINESOURCES(ILINE)%X2
                     ELSE
                        Y1 = LINESOURCES(ILINE)%Y2
                        Y2 = LINESOURCES(ILINE)%Y1
                        X1 = LINESOURCES(ILINE)%X2
                        X2 = LINESOURCES(ILINE)%X1
                     END IF
                     Y = SQRT(Y1*Y1 + R*(Y2*Y2 - Y1*Y1))
                     X = X1 + (Y-Y1)/(Y2-Y1)*(X2-X1)
                     Z = 0
                  END IF
               ELSE
                  X = LINESOURCES(ILINE)%X1 + R*(LINESOURCES(ILINE)%X2-LINESOURCES(ILINE)%X1)
                  Y = LINESOURCES(ILINE)%Y1 + R*(LINESOURCES(ILINE)%Y2-LINESOURCES(ILINE)%Y1)
                  Z = ZMIN + (ZMAX-ZMIN)*rf()
               END IF
   
               CALL CELL_FROM_POSITION(X,Y,  IC)
   
               ! Init a particle object and assign it to the local vector of particles
               CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTFRAC,  particleNOW)
               IF ((tID == DUMP_TRAJECTORY_START) .AND. (DUMP_COUNTER < DUMP_TRAJECTORY_NUMBER)) THEN
                  particleNOW%DUMP_TRAJ = .TRUE.
                  DUMP_COUNTER = DUMP_COUNTER + 1
               END IF
               CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles)
   
            END DO

            LOCAL_LINE_EMIT_COUNT(ILINE+N_LINESOURCES*(S_ID-1)) = LOCAL_LINE_EMIT_COUNT(ILINE+N_LINESOURCES*(S_ID-1)) + NFS

         END DO
   
      END DO


      ! Compute total of particles injected
      CALL MPI_REDUCE(LOCAL_LINE_EMIT_COUNT, LINE_EMIT_COUNT, N_LINESOURCES*N_SPECIES, MPI_INTEGER, MPI_SUM, &
      0, MPI_COMM_WORLD, ierr)

      DEALLOCATE(LOCAL_LINE_EMIT_COUNT)
         
   END SUBROUTINE LINE_SOURCE_INJECT


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE BOUNDARIES_INJECT -> Injects particles from domain borders !!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Deprecated.
   SUBROUTINE BOUNDARIES_INJECT
   
      IMPLICIT NONE

      INTEGER      :: IP, IC, IS, NFS
      REAL(KIND=8) :: DTFRAC, Vdummy 
      REAL(KIND=8) :: X, Y, Z, VX, VY, VZ, EROT, EVIB 
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW

      INTEGER :: S_ID
      REAL(KIND=8) :: M


      ! Low X boundary (left)
      IF (BOOL_INJ_XMIN) THEN

         DO IS = 1, MIXTURES(MIX_BOUNDINJECT)%N_COMPONENTS ! Loop on mixture components
            S_ID = MIXTURES(MIX_BOUNDINJECT)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            NFS = FLOOR(nfs_XMIN(IS))
            IF (nfs_XMIN(IS)-REAL(NFS, KIND=8) .GE. rf()) THEN ! Same as SPARTA's perspeciess
               NFS = NFS + 1
            END IF

            DO IP = 1, NFS ! Loop on particles to be injected

               CALL MAXWELL(UX_BOUND,UY_BOUND,UZ_BOUND, &
                           TTRA_BOUND,TTRA_BOUND,TTRA_BOUND, &
                           Vdummy,VY,VZ,M)
               
               CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, TROT_BOUND, EROT)
               CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, TVIB_BOUND, EVIB)

               VX = UX_BOUND + FLX(S_NORM_XMIN,TTRA_BOUND,M) !!AAAAAAAAAAA

               DTFRAC = rf()*DT
               X = XMIN
               Y = YMIN + (YMAX-YMIN)*rf()
               Z = ZMIN + (ZMAX-ZMIN)*rf()

               CALL CELL_FROM_POSITION(X,Y,  IC)

               ! Init a particle object and assign it to the local vector of particles
               CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTFRAC,  particleNOW)
               CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles)

            END DO
         END DO

      END IF
         
      ! High X boundary (right)
      IF (BOOL_INJ_XMAX) THEN


         DO IS = 1, MIXTURES(MIX_BOUNDINJECT)%N_COMPONENTS ! Loop on mixture components
            S_ID = MIXTURES(MIX_BOUNDINJECT)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            NFS = FLOOR(nfs_XMAX(IS))
            IF (nfs_XMAX(IS)-REAL(NFS, KIND=8) .GE. rf()) THEN
               NFS = NFS + 1
            END IF

            DO IP = 1, NFS ! Loop on particles to be injected

               CALL MAXWELL(UX_BOUND,UY_BOUND,UZ_BOUND, &
                           TTRA_BOUND,TTRA_BOUND,TTRA_BOUND, &
                           Vdummy,VY,VZ,M)

               CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, TROT_BOUND, EROT)
               CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, TVIB_BOUND, EVIB)

               VX = UX_BOUND - FLX(S_NORM_XMAX,TTRA_BOUND,M) !! I think something was wrong here with the sign

               DTFRAC = rf()*DT
               X = XMAX
               Y = YMIN + (YMAX-YMIN)*rf()
               Z = ZMIN + (ZMAX-ZMIN)*rf()
               
               CALL CELL_FROM_POSITION(X,Y,  IC)

               ! Init a particle object and assign it to vector of particles
               CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTFRAC,  particleNOW)
               CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles)
               
            END DO
         END DO

      END IF

      ! Low Y boundary (bottom)
      IF (BOOL_INJ_YMIN) THEN

         DO IS = 1, MIXTURES(MIX_BOUNDINJECT)%N_COMPONENTS ! Loop on mixture components
            S_ID = MIXTURES(MIX_BOUNDINJECT)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            NFS = FLOOR(nfs_YMIN(IS))
            IF (nfs_YMIN(IS)-REAL(NFS, KIND=8) .GE. rf()) THEN
               NFS = NFS + 1
            END IF

            DO IP = 1, NFS ! Loop on particles to be injected

               CALL MAXWELL(UX_BOUND,UY_BOUND,UZ_BOUND, &
                           TTRA_BOUND,TTRA_BOUND,TTRA_BOUND, &
                           VX,Vdummy,VZ,M)

               CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, TROT_BOUND, EROT)
               CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, TVIB_BOUND, EVIB)

               VY = UY_BOUND + FLX(S_NORM_YMIN,TTRA_BOUND,M) !!AAAAAAAAAAA

               DTFRAC = rf()*DT
               X = XMIN + (XMAX-XMIN)*rf()
               Y = YMIN
               Z = ZMIN + (ZMAX-ZMIN)*rf()

               CALL CELL_FROM_POSITION(X,Y,  IC)

               ! Init a particle object and assign it to the local vector of particles
               CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTFRAC,  particleNOW)
               CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles)
            END DO
         END DO

      END IF
         
      ! High Y boundary (top)
      IF (BOOL_INJ_YMAX) THEN

         DO IS = 1, MIXTURES(MIX_BOUNDINJECT)%N_COMPONENTS ! Loop on mixture components
            S_ID = MIXTURES(MIX_BOUNDINJECT)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            NFS = FLOOR(nfs_YMAX(IS))
            IF (nfs_YMAX(IS)-REAL(NFS, KIND=8) .GE. rf()) THEN
               NFS = NFS + 1
            END IF

            DO IP = 1, NFS ! Loop on particles to be injected

               CALL MAXWELL(UX_BOUND,UY_BOUND,UZ_BOUND, &
                           TTRA_BOUND,TTRA_BOUND,TTRA_BOUND, &
                           VX,Vdummy,VZ,M)

               CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, TROT_BOUND, EROT)
               CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, TVIB_BOUND, EVIB)

               VY = UY_BOUND - FLX(S_NORM_YMAX,TTRA_BOUND,M) !!AAAAAAAAAAA

               DTFRAC = rf()*DT
               X = XMIN + (XMAX-XMIN)*rf() ! There was a bug here!
               Y = YMAX
               Z = ZMIN + (ZMAX-ZMIN)*rf()

               CALL CELL_FROM_POSITION(X,Y,  IC)

               ! Init a particle object and assign it to vector of particles
               CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTFRAC,  particleNOW)
               CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles)
            END DO
         END DO

      END IF

   END SUBROUTINE BOUNDARIES_INJECT
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ADVECT -> Advects particles in the domain     !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ADVECT

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

      INTEGER      :: IP, IC, i, SOL, OLD_IC
      INTEGER      :: BOUNDCOLL, WALLCOLL, GOODSOL, FACE_PG, NEIGHBOR
      REAL(KIND=8) :: DTCOLL, TOTDTCOLL, CANDIDATE_DTCOLL, rfp
      REAL(KIND=8) :: COEFA, COEFB, COEFC, DELTA, SOL1, SOL2, ALPHA, BETA
      REAL(KIND=8), DIMENSION(2) :: TEST
      REAL(KIND=8), DIMENSION(4) :: NORMX, NORMY, XW, YW, ZW, BOUNDPOS
      ! REAL(KIND=8) :: XCOLL, YCOLL, ZCOLL
      REAL(KIND=8) :: VN, DX
      LOGICAL, DIMENSION(:), ALLOCATABLE :: REMOVE_PART
      REAL(KIND=8), DIMENSION(3) :: V_OLD, V_NEW
      REAL(KIND=8), DIMENSION(3) :: E, B
      REAL(KIND=8), DIMENSION(3) :: FACE_NORMAL, FACE_TANG1, FACE_TANG2
      REAL(KIND=8), DIMENSION(3) :: TANG1, TANG2
      REAL(KIND=8) :: VDOTTANG1, VRM, RN, R1, R2, THETA1, THETA2, DOT_NORM, VTANGENT
      REAL(KIND=8) :: V_NORM, V_TANG1, V_TANG2, V_PERP, VZ, VDUMMY, EROT, EVIB, VDOTN, WALL_TEMP
      INTEGER :: S_ID
      LOGICAL :: HASCOLLIDED
      REAL(KIND=8) :: XCOLL, YCOLL, COLLDIST, EDGE_X1, EDGE_Y1
      INTEGER, DIMENSION(:), ALLOCATABLE :: LOCAL_BOUNDARY_COLL_COUNT, LOCAL_WALL_COLL_COUNT
      REAL(KIND=8) :: WEIGHT_RATIO
      TYPE(PARTICLE_DATA_STRUCTURE) :: NEWparticle, particleNOW
      LOGICAL :: FLUIDBOUNDARY
      INTEGER :: NEIGHBORPG
      REAL(KIND=8) :: CHARGE, K, PSIP, RHO_Q
      INTEGER :: VP

      REAL(KIND=8) :: VXPRE, VYPRE, VZPRE

      
      REAL(KIND=8) :: TOL = 1.0d-15

      E = [0.d0, 0.d0, 0.d0]
      B = [0.d0, 0.d0, 0.d0]


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

      FIELD_POWER = 0

      ! OPEN(66341, FILE='washboarddump', POSITION='append', STATUS='unknown', ACTION='write')


      DO IP = 1, NP_PROC
         REMOVE_PART(IP) = .FALSE.

         ! Update velocity
         IC = particles(IP)%IC

         V_NEW(1) = particles(IP)%VX
         V_NEW(2) = particles(IP)%VY
         V_NEW(3) = particles(IP)%VZ

         IF (PIC_TYPE .NE. NONE) THEN
            !PHIBAR_FIELD = 0.d0
            !EBAR_FIELD = 0.d0
            CALL APPLY_E_FIELD(IP, E)
            B = 0
            IF (N_SOLENOIDS > 0) CALL APPLY_B_FIELD(IP, B)
            IF (BOOL_MAGNETIC_DIPOLE) CALL APPLY_B_DIPOLE_FIELD(IP, B)
            B = B + EXTERNAL_B_FIELD
            ! CALL APPLY_RF_EB_FIELD(particles, IP, E, B)
            


            V_OLD(1) = particles(IP)%VX
            V_OLD(2) = particles(IP)%VY
            V_OLD(3) = particles(IP)%VZ

            CALL UPDATE_VELOCITY_BORIS(particles(IP)%DTRIM, V_OLD, V_NEW, &
            SPECIES(particles(IP)%S_ID)%CHARGE, SPECIES(particles(IP)%S_ID)%MOLECULAR_MASS, &
            E, B)

            ! Assign v^n to the particle, for simplicity
            !particles(IP)%VX = 0.5*(V_OLD(1) + V_NEW(1))
            !particles(IP)%VY = 0.5*(V_OLD(2) + V_NEW(2))
            !particles(IP)%VZ = 0.5*(V_OLD(3) + V_NEW(3))
            particles(IP)%VX = V_NEW(1)
            particles(IP)%VY = V_NEW(2)
            particles(IP)%VZ = V_NEW(3)
         END IF



         ! ! Forced electric field
         ! E = [3000.d0, 0.d0, 0.d0]
            


         ! V_OLD(1) = particles(IP)%VX
         ! V_OLD(2) = particles(IP)%VY
         ! V_OLD(3) = particles(IP)%VZ

         ! CALL UPDATE_VELOCITY_BORIS(DT, V_OLD, V_NEW, &
         ! SPECIES(particles(IP)%S_ID)%CHARGE, SPECIES(particles(IP)%S_ID)%MOLECULAR_MASS, &
         ! E, B)

         ! ! Assign v^n to the particle, for simplicity
         ! !particles(IP)%VX = 0.5*(V_OLD(1) + V_NEW(1))
         ! !particles(IP)%VY = 0.5*(V_OLD(2) + V_NEW(2))
         ! !particles(IP)%VZ = 0.5*(V_OLD(3) + V_NEW(3))
         ! particles(IP)%VX = V_NEW(1)
         ! particles(IP)%VY = V_NEW(2)
         ! particles(IP)%VZ = V_NEW(3)

         ! ! End forced electric field.



         HASCOLLIDED = .FALSE.
         TOTDTCOLL = 0.
         DO WHILE (particles(IP)%DTRIM .GT. 0.) ! Repeat the procedure until step is done

            DTCOLL = particles(IP)%DTRIM ! Looking for collisions within the remaining time
            ! ______ ADVECTION ______
            IF (GRID_TYPE == UNSTRUCTURED) THEN
               !WRITE(*,*) 'Moving particle ', IP, ' for ', DTCOLL, ' s. Position: ', particles(IP)%X, particles(IP)%Y,&
               !' velocity: ', particles(IP)%VX, particles(IP)%VY
               ! For unstructured, we only need to check the boundaries of the cell.
               BOUNDCOLL = -1
               
                  
               IF (AXI) THEN
                  DO I = 1, 3
                     !!!! DBDBDBDBDDBBDBDBDBDBDBDBDBDBDBDDB work on this. To be adapted to the three sides of the cell.
                     ! We are axisymmetric
                     CANDIDATE_DTCOLL = DTCOLL
                     ! Compute auxiliary parameters
                     EDGE_X1 = U2D_GRID%NODE_COORDS(1,U2D_GRID%CELL_NODES(I,IC))
                     EDGE_Y1 = U2D_GRID%NODE_COORDS(2,U2D_GRID%CELL_NODES(I,IC))

                     IF (ABS(U2D_GRID%EDGE_NORMAL(2,I,IC)) < 1.d-10 ) THEN
                        ! Vertical wall
                        SOL1 = (EDGE_X1 - particles(IP)%X)/particles(IP)%VX
                        IF (SOL1 >= 0 .AND. SOL1 < DTCOLL) THEN
                           GOODSOL = 1
                           TEST(1) = SOL1
                        ELSE
                           GOODSOL = 0
                        END IF
                     ELSE
                        ! Non-vertical wall
                        ALPHA = -U2D_GRID%EDGE_NORMAL(1,I,IC)/U2D_GRID%EDGE_NORMAL(2,I,IC)
                        BETA  = (particles(IP)%X - EDGE_X1)*ALPHA

                        COEFA = particles(IP)%VY**2 + particles(IP)%VZ**2 - ALPHA*ALPHA*particles(IP)%VX**2
                        COEFB = particles(IP)%Y * particles(IP)%VY - (BETA + EDGE_Y1)*ALPHA*particles(IP)%VX
                        COEFC = particles(IP)%Y**2 - (BETA + EDGE_Y1)**2
                        DELTA = COEFB*COEFB-COEFA*COEFC
                        IF (DELTA .GE. 0) THEN
                           ! Compute the solutions, check if they are any good and, in case, order them to be checked further.
                           SOL1 = (-COEFB - SQRT(DELTA))/COEFA
                           SOL2 = (-COEFB + SQRT(DELTA))/COEFA
                           IF (SOL1 >= -TOL .AND. SOL1 < DTCOLL .AND. BETA+ALPHA*particles(IP)%VX*SOL1+EDGE_Y1 >= 0) THEN
                              IF (SOL2 >= -TOL .AND. SOL2 < DTCOLL .AND. BETA+ALPHA*particles(IP)%VX*SOL2+EDGE_Y1 >= 0) THEN
                                 GOODSOL = 2
                                 IF (SOL1 <= SOL2) THEN
                                    TEST(1) = SOL1
                                    TEST(2) = SOL2
                                 ELSE
                                    TEST(1) = SOL2
                                    TEST(2) = SOL1
                                 END IF
                              ELSE
                                 GOODSOL = 1
                                 TEST(1) = SOL1
                              ENDIF
                           ELSE
                              IF (SOL2 >= -TOL .AND. SOL2 < DTCOLL .AND. BETA+ALPHA*particles(IP)%VX*SOL2+EDGE_Y1 >= 0) THEN
                                 GOODSOL = 1
                                 TEST(1) = SOL2
                              ELSE
                                 GOODSOL = 0
                              END IF 
                           END IF
                        ELSE
                           GOODSOL = 0
                        END IF
                     END IF

                     ! Now further checks for each of the good solutions:
                     ! - if the velocity at collision time is actually towards the surface
                     ! - if the collision point is within the extremes of the segment
                     ! Note: this requires moving the particle to the collision point and then moving back.
                     IF (GOODSOL /= 0) THEN
                        DO SOL = 1, GOODSOL
                           CALL MOVE_PARTICLE(IP, TEST(SOL))
                           IF ((particles(IP)%VX*U2D_GRID%EDGE_NORMAL(1,I,IC) &
                              + particles(IP)%VY*U2D_GRID%EDGE_NORMAL(2,I,IC)) > 0) THEN

                              COLLDIST = -U2D_GRID%EDGE_NORMAL(2,I,IC)*(particles(IP)%X-EDGE_X1)/U2D_GRID%CELL_EDGES_LEN(I,IC) &
                                         +U2D_GRID%EDGE_NORMAL(1,I,IC)*(particles(IP)%Y-EDGE_Y1)/U2D_GRID%CELL_EDGES_LEN(I,IC)
                              IF ((COLLDIST .GE. 0) .AND. (COLLDIST .LE. 1)) THEN
                                 ! Collision happens!
                                 DTCOLL = TEST(SOL)
                                 BOUNDCOLL = I  
                                 TOTDTCOLL = TOTDTCOLL + DTCOLL  
                                 HASCOLLIDED = .TRUE.
                              END IF
                           END IF
                           CALL MOVE_PARTICLE(IP, -TEST(SOL))
                        END DO
                     END IF
                  END DO
               ELSE
                  ! We are not axisymmetric (valid also for simplified axisymmetric procedure)
               
                  IF (DIMS == 1) THEN
                     DO I = 1, 2
                        VN = particles(IP)%VX * U1D_GRID%EDGE_NORMAL(1,I,IC)
                        ! Compute the distance from the boundary
                        DX = (U1D_GRID%NODE_COORDS(1,U1D_GRID%CELL_NODES(I,IC)) - particles(IP)%X) * U1D_GRID%EDGE_NORMAL(1,I,IC)

                        ! Check if a collision happens (sooner than previously calculated)
                        IF (VN .GE. 0. .AND. VN * DTCOLL .GE. DX) THEN
                           DTCOLL = DX/VN
                           BOUNDCOLL = I
                           TOTDTCOLL = TOTDTCOLL + DTCOLL  
                           HASCOLLIDED = .TRUE.     
                        END IF
                     END DO
                  ELSE IF (DIMS == 2) THEN
                     DO I = 1, 3
                        VN = particles(IP)%VX * U2D_GRID%EDGE_NORMAL(1,I,IC) &
                           + particles(IP)%VY * U2D_GRID%EDGE_NORMAL(2,I,IC)
                        ! Compute the distance from the boundary
                        DX = (U2D_GRID%NODE_COORDS(1,U2D_GRID%CELL_NODES(I,IC)) - particles(IP)%X) * U2D_GRID%EDGE_NORMAL(1,I,IC) &
                           + (U2D_GRID%NODE_COORDS(2,U2D_GRID%CELL_NODES(I,IC)) - particles(IP)%Y) * U2D_GRID%EDGE_NORMAL(2,I,IC)

                        ! Check if a collision happens (sooner than previously calculated)
                        IF (VN .GE. 0. .AND. VN * DTCOLL .GE. DX) THEN
                           DTCOLL = DX/VN
                           BOUNDCOLL = I
                           TOTDTCOLL = TOTDTCOLL + DTCOLL  
                           HASCOLLIDED = .TRUE.     
                        END IF
                     END DO
                  ELSE IF (DIMS == 3) THEN
                     DO I = 1, 4
                        VN = particles(IP)%VX * U3D_GRID%FACE_NORMAL(1,I,IC) &
                           + particles(IP)%VY * U3D_GRID%FACE_NORMAL(2,I,IC) &
                           + particles(IP)%VZ * U3D_GRID%FACE_NORMAL(3,I,IC)
                        ! Compute the distance from the boundary
                        DX = (U3D_GRID%NODE_COORDS(1,U3D_GRID%CELL_NODES(I,IC)) - particles(IP)%X) * U3D_GRID%FACE_NORMAL(1,I,IC) &
                           + (U3D_GRID%NODE_COORDS(2,U3D_GRID%CELL_NODES(I,IC)) - particles(IP)%Y) * U3D_GRID%FACE_NORMAL(2,I,IC) &
                           + (U3D_GRID%NODE_COORDS(3,U3D_GRID%CELL_NODES(I,IC)) - particles(IP)%Z) * U3D_GRID%FACE_NORMAL(3,I,IC)

                        ! Check if a collision happens (sooner than previously calculated)
                        IF (VN .GE. 0. .AND. VN * DTCOLL .GE. DX) THEN
                           DTCOLL = DX/VN
                           BOUNDCOLL = I
                           TOTDTCOLL = TOTDTCOLL + DTCOLL  
                           HASCOLLIDED = .TRUE.     
                        END IF
                     END DO
                  END IF

               END IF

               IF (BOUNDCOLL .NE. -1) THEN
                  IF (particles(IP)%Y == 0.d0 .AND. DTCOLL == 0.d0) WRITE(*,*) 'Here 1.'
                  CALL MOVE_PARTICLE(IP, DTCOLL)
                  particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

                  FLUIDBOUNDARY = .FALSE.
                  IF (DIMS == 1) THEN
                     NEIGHBOR = U1D_GRID%CELL_NEIGHBORS(BOUNDCOLL, IC)
                     IF (NEIGHBOR == -1) THEN
                        FLUIDBOUNDARY = .TRUE.
                     ELSE
                        NEIGHBORPG = U1D_GRID%CELL_PG(NEIGHBOR)
                        IF (NEIGHBORPG .NE. -1) THEN
                           IF (GRID_BC(NEIGHBORPG)%VOLUME_BC == SOLID) FLUIDBOUNDARY = .TRUE.
                        END IF
                     END IF

                     FACE_PG = U1D_GRID%CELL_EDGES_PG(BOUNDCOLL, IC)
                     FACE_NORMAL = -U1D_GRID%EDGE_NORMAL(:,BOUNDCOLL,IC)
                     FACE_TANG1 = [0.d0, FACE_NORMAL(1), 0.d0]
                     FACE_TANG2 = [0.d0, 0.d0, 1.d0]
                  ELSE IF (DIMS == 2) THEN
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
                     FACE_NORMAL = -U2D_GRID%EDGE_NORMAL(:,BOUNDCOLL,IC)
                     FACE_TANG1 = [-FACE_NORMAL(2), FACE_NORMAL(1), 0.d0]
                     FACE_TANG2 = [0.d0, 0.d0, 1.d0]
                  ELSE IF (DIMS == 3) THEN
                     NEIGHBOR = U3D_GRID%CELL_NEIGHBORS(BOUNDCOLL, IC)
                     IF (NEIGHBOR == -1) THEN
                        FLUIDBOUNDARY = .TRUE.
                     ELSE
                        NEIGHBORPG = U3D_GRID%CELL_PG(NEIGHBOR)
                        IF (NEIGHBORPG .NE. -1) THEN
                           IF (GRID_BC(NEIGHBORPG)%VOLUME_BC == SOLID) FLUIDBOUNDARY = .TRUE.
                        END IF
                     END IF

                     FACE_PG = U3D_GRID%CELL_FACES_PG(BOUNDCOLL,IC)
                     FACE_NORMAL = -U3D_GRID%FACE_NORMAL(:,BOUNDCOLL,IC)
                     FACE_TANG1 = U3D_GRID%FACE_TANG1(:,BOUNDCOLL,IC)
                     FACE_TANG2 = U3D_GRID%FACE_TANG2(:,BOUNDCOLL,IC)
                  END IF


                  ! Apply boundary conditions
                  IF (FLUIDBOUNDARY) THEN

                     IF (FACE_PG .NE. -1) THEN

                        IF (GRID_BC(FACE_PG)%DUMP_FLUXES .AND. (tID .GE. DUMP_PART_BOUND_START)) THEN
                           particleNOW = particles(IP)
                           CALL ADD_PARTICLE_ARRAY(particleNOW, NP_DUMP_PROC, part_dump)
                        END IF


                        ! Tally incident particle fluxes to boundary
                        IF ((tID .GE. DUMP_GRID_START) .AND. (tID .NE. RESTART_TIMESTEP)) THEN
                           IF (MOD(tID-DUMP_BOUND_START, DUMP_BOUND_AVG_EVERY) .EQ. 0) THEN
                              CALL TALLY_PARTICLE_TO_BOUNDARY(.FALSE., particles(IP), IC, BOUNDCOLL)
                           END IF
                        END IF

                        
                        CHARGE = SPECIES(particles(IP)%S_ID)%CHARGE
                        IF ( (GRID_BC(FACE_PG)%FIELD_BC == DIELECTRIC_BC &
                           .OR.GRID_BC(FACE_PG)%FIELD_BC == CONDUCTIVE_BC) .AND. ABS(CHARGE) .GE. 1.d-6) THEN
                           K = QE/(EPS0*EPS_SCALING**2)
                           IF (DIMS == 1) THEN
                              RHO_Q = K*CHARGE*FNUM/(YMAX-YMIN)/(ZMAX-ZMIN)
                              DO I = 1, 2
                                 VP = U1D_GRID%CELL_NODES(I,IC)
                                 PSIP = U1D_GRID%BASIS_COEFFS(1,I,IC)*particles(IP)%X &
                                      + U1D_GRID%BASIS_COEFFS(3,I,IC)
                                 SURFACE_CHARGE(VP) = SURFACE_CHARGE(VP) + RHO_Q*PSIP
                              END DO
                           ELSE IF (DIMS == 2) THEN
                              RHO_Q = K*CHARGE*FNUM/(ZMAX-ZMIN)
                              DO I = 1, 3
                                 VP = U2D_GRID%CELL_NODES(I,IC)
                                 PSIP = U2D_GRID%BASIS_COEFFS(1,I,IC)*particles(IP)%X &
                                      + U2D_GRID%BASIS_COEFFS(2,I,IC)*particles(IP)%Y &
                                      + U2D_GRID%BASIS_COEFFS(3,I,IC)
                                 SURFACE_CHARGE(VP) = SURFACE_CHARGE(VP) + RHO_Q*PSIP
                              END DO
                           ELSE IF (DIMS == 3) THEN
                              RHO_Q = K*CHARGE*FNUM
                              DO I = 1, 4
                                 VP = U3D_GRID%CELL_NODES(I,IC)
                                 PSIP = U3D_GRID%BASIS_COEFFS(1,I,IC)*particles(IP)%X &
                                      + U3D_GRID%BASIS_COEFFS(2,I,IC)*particles(IP)%Y &
                                      + U3D_GRID%BASIS_COEFFS(3,I,IC)*particles(IP)%Z &
                                      + U3D_GRID%BASIS_COEFFS(4,I,IC)
                                 SURFACE_CHARGE(VP) = SURFACE_CHARGE(VP) + RHO_Q*PSIP
                              END DO
                           END IF
                        END IF

                        ! Apply particle boundary condition
                        IF (GRID_BC(FACE_PG)%PARTICLE_BC == SPECULAR) THEN
                           IF (GRID_BC(FACE_PG)%REACT) THEN
                              CALL WALL_REACT(particles, IP, REMOVE_PART(IP))
                           END IF
                           
                           VDOTN = particles(IP)%VX*FACE_NORMAL(1) &
                                 + particles(IP)%VY*FACE_NORMAL(2) &
                                 + particles(IP)%VZ*FACE_NORMAL(3)
                           particles(IP)%VX = particles(IP)%VX - 2.*VDOTN*FACE_NORMAL(1)
                           particles(IP)%VY = particles(IP)%VY - 2.*VDOTN*FACE_NORMAL(2)
                           particles(IP)%VZ = particles(IP)%VZ - 2.*VDOTN*FACE_NORMAL(3)

                        ELSE IF (GRID_BC(FACE_PG)%PARTICLE_BC == DIFFUSE) THEN
                           IF (GRID_BC(FACE_PG)%REACT) THEN
                              CALL WALL_REACT(particles, IP, REMOVE_PART(IP))
                           END IF
                           
                           VXPRE = particles(IP)%VX
                           VYPRE = particles(IP)%VY
                           VZPRE = particles(IP)%VZ

                           S_ID = particles(IP)%S_ID
                           WALL_TEMP = GRID_BC(FACE_PG)%WALL_TEMP
                           CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                           WALL_TEMP, WALL_TEMP, WALL_TEMP, &
                           VDUMMY, V_TANG1, V_TANG2, SPECIES(S_ID)%MOLECULAR_MASS)

                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, WALL_TEMP, EROT)
                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, WALL_TEMP, EVIB)
                                          
                           V_PERP = FLX(0.d0, WALL_TEMP, SPECIES(S_ID)%MOLECULAR_MASS)


                           particles(IP)%VX = V_PERP*FACE_NORMAL(1) &
                                            + V_TANG1*FACE_TANG1(1) &
                                            + V_TANG2*FACE_TANG2(1)
                           particles(IP)%VY = V_PERP*FACE_NORMAL(2) &
                                            + V_TANG1*FACE_TANG1(2) &
                                            + V_TANG2*FACE_TANG2(2)
                           particles(IP)%VZ = V_PERP*FACE_NORMAL(3) &
                                            + V_TANG1*FACE_TANG1(3) &
                                            + V_TANG2*FACE_TANG2(3)
                           
                           particles(IP)%EROT = EROT
                           particles(IP)%EVIB = EVIB

                           ! IF (rf() < 0.01) THEN
                           !    WRITE(66341,*) VXPRE, ',', VYPRE, ',', VZPRE, ',', &
                           !    particles(IP)%VX, ',', particles(IP)%VY, ',', particles(IP)%VZ, ',',&
                           !    FACE_NORMAL(1), ',', FACE_NORMAL(2), ',', FACE_NORMAL(3)
                           ! END IF

                        ELSE IF (GRID_BC(FACE_PG)%PARTICLE_BC == CLL) THEN
                           IF (GRID_BC(FACE_PG)%REACT) THEN
                              CALL WALL_REACT(particles, IP, REMOVE_PART(IP))
                           END IF

                           !VXPRE = particles(IP)%VX
                           !VYPRE = particles(IP)%VY
                           !VZPRE = particles(IP)%VZ

                           VDOTN = particles(IP)%VX*FACE_NORMAL(1) &
                                 + particles(IP)%VY*FACE_NORMAL(2) &
                                 + particles(IP)%VZ*FACE_NORMAL(3)

                           TANG1(1) = particles(IP)%VX - VDOTN*FACE_NORMAL(1)
                           TANG1(2) = particles(IP)%VY - VDOTN*FACE_NORMAL(2)
                           TANG1(3) = particles(IP)%VZ - VDOTN*FACE_NORMAL(3)

                           TANG1 = TANG1/NORM2(TANG1)

                           TANG2 = CROSS(FACE_NORMAL, TANG1)

                           VDOTTANG1 = particles(IP)%VX*TANG1(1) &
                                     + particles(IP)%VY*TANG1(2) &
                                     + particles(IP)%VZ*TANG1(3)

                           S_ID = particles(IP)%S_ID
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

                           particles(IP)%VX = V_PERP*FACE_NORMAL(1) &
                                            + V_TANG1*TANG1(1) &
                                            + V_TANG2*TANG2(1)
                           particles(IP)%VY = V_PERP*FACE_NORMAL(2) &
                                            + V_TANG1*TANG1(2) &
                                            + V_TANG2*TANG2(2)
                           particles(IP)%VZ = V_PERP*FACE_NORMAL(3) &
                                            + V_TANG1*TANG1(3) &
                                            + V_TANG2*TANG2(3)
                           
                           particles(IP)%EROT = EROT
                           particles(IP)%EVIB = EVIB

                           !OPEN(66341, FILE='clldump', POSITION='append', STATUS='unknown', ACTION='write')
                           !WRITE(66341,*) VXPRE, ', ', VYPRE, ', ', VZPRE, ', ', &
                           !particles(IP)%VX, ', ', particles(IP)%VY, ', ', particles(IP)%VZ
                           !CLOSE(66341)
                        ELSE IF (GRID_BC(FACE_PG)%PARTICLE_BC == WB_BC) THEN
                           IF (GRID_BC(FACE_PG)%REACT) THEN
                              CALL WALL_REACT(particles, IP, REMOVE_PART(IP))
                           END IF

                           VXPRE = particles(IP)%VX
                           VYPRE = particles(IP)%VY
                           VZPRE = particles(IP)%VZ

                           VDOTN = particles(IP)%VX*FACE_NORMAL(1) &
                                 + particles(IP)%VY*FACE_NORMAL(2) &
                                 + particles(IP)%VZ*FACE_NORMAL(3)

                           TANG1(1) = particles(IP)%VX - VDOTN*FACE_NORMAL(1)
                           TANG1(2) = particles(IP)%VY - VDOTN*FACE_NORMAL(2)
                           TANG1(3) = particles(IP)%VZ - VDOTN*FACE_NORMAL(3)

                           TANG1 = TANG1/NORM2(TANG1)

                           TANG2 = CROSS(FACE_NORMAL, TANG1)

                           VDOTTANG1 = particles(IP)%VX*TANG1(1) &
                                     + particles(IP)%VY*TANG1(2) &
                                     + particles(IP)%VZ*TANG1(3)


                           V_TANG1 = VDOTTANG1
                           V_PERP = VDOTN
                           CALL WB_SCATTER(FACE_PG, SPECIES(particles(IP)%S_ID)%MOLECULAR_MASS, &
                           V_TANG1, V_TANG2, V_PERP)
                           S_ID = particles(IP)%S_ID
                           WALL_TEMP = GRID_BC(FACE_PG)%WALL_TEMP

                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, WALL_TEMP, EROT)
                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, WALL_TEMP, EVIB)

                           particles(IP)%VX = V_PERP*FACE_NORMAL(1) &
                                            + V_TANG1*TANG1(1) &
                                            + V_TANG2*TANG2(1)
                           particles(IP)%VY = V_PERP*FACE_NORMAL(2) &
                                            + V_TANG1*TANG1(2) &
                                            + V_TANG2*TANG2(2)
                           particles(IP)%VZ = V_PERP*FACE_NORMAL(3) &
                                            + V_TANG1*TANG1(3) &
                                            + V_TANG2*TANG2(3)
                           
                           particles(IP)%EROT = EROT
                           particles(IP)%EVIB = EVIB

                           ! IF (rf() < 0.01) THEN
                           !    WRITE(66341,*) VXPRE, ',', VYPRE, ',', VZPRE, ',', &
                           !    particles(IP)%VX, ',', particles(IP)%VY, ',', particles(IP)%VZ, ',',&
                           !    FACE_NORMAL(1), ',', FACE_NORMAL(2), ',', FACE_NORMAL(3)
                           ! END IF
                           
                        ELSE
                           REMOVE_PART(IP) = .TRUE.
                           particles(IP)%DTRIM = 0.d0
                        END IF


                        ! Tally reflected particle fluxes to boundary
                        IF (.NOT. REMOVE_PART(IP)) THEN
                           IF ((tID .GE. DUMP_GRID_START) .AND. (tID .NE. RESTART_TIMESTEP)) THEN
                              IF (MOD(tID-DUMP_BOUND_START, DUMP_BOUND_AVG_EVERY) .EQ. 0) THEN
                                 CALL TALLY_PARTICLE_TO_BOUNDARY(.TRUE., particles(IP), IC, BOUNDCOLL)
                              END IF
                           END IF
                        END IF
                     ELSE
                        REMOVE_PART(IP) = .TRUE.
                        particles(IP)%DTRIM = 0.d0
                     END IF
                  ELSE
                     particles(IP)%IC = NEIGHBOR
                     IC = NEIGHBOR
                     !WRITE(*,*) 'moved particle to cell: ', IC
                  END IF
               ELSE
                  IF (particles(IP)%Y == 0.d0 .AND. DTCOLL == 0.d0) WRITE(*,*) 'Here 2.'
                  CALL MOVE_PARTICLE(IP, particles(IP)%DTRIM)
                  particles(IP)%DTRIM = 0.d0
               END IF

            ELSE ! Grid type is rectilinear
               BOUNDCOLL = -1
               DO I = 1, 2*DIMS ! Check collisions with boundaries (xmin, xmax, ymin, ymax)
                  IF (DIMS == 2 .AND. AXI) THEN
                     ! We are axisymmetric. Check is more complicated.
                     CANDIDATE_DTCOLL = DTCOLL
                     
                     IF (i==1 .OR. i==2) THEN
                        ! Vertical wall
                        SOL1 = (BOUNDPOS(i)-particles(IP)%X)/particles(IP)%VX
                        IF (SOL1 >= 0 .AND. SOL1 < DTCOLL) THEN
                           GOODSOL = 1
                           TEST(1) = SOL1
                        ELSE
                           GOODSOL = 0
                        END IF
                     ELSE
                        ! Horizontal boundary
                        COEFA = particles(IP)%VY**2 + particles(IP)%VZ**2
                        COEFB = particles(IP)%Y * particles(IP)%VY
                        COEFC = particles(IP)%Y**2 - BOUNDPOS(i)**2
                        DELTA = COEFB*COEFB-COEFA*COEFC
                        IF (DELTA .GE. 0) THEN
                        ! Compute the solutions, check if they are any good and, in case, order them to be checked further.
                           SOL1 = (-COEFB - SQRT(DELTA))/COEFA
                           SOL2 = (-COEFB + SQRT(DELTA))/COEFA
                           IF (SOL1 >= 0 .AND. SOL1 < DTCOLL) THEN
                              IF (SOL2 >= 0 .AND. SOL2 < DTCOLL) THEN
                                 GOODSOL = 2
                                 IF (SOL1 <= SOL2) THEN
                                    TEST(1) = SOL1
                                    TEST(2) = SOL2
                                 ELSE
                                    TEST(1) = SOL2
                                    TEST(2) = SOL1
                                 END IF
                              ELSE
                                 GOODSOL = 1
                                 TEST(1) = SOL1
                              ENDIF
                           ELSE
                              IF (SOL2 >= 0 .AND. SOL2 < DTCOLL) THEN
                                 GOODSOL = 1
                                 TEST(1) = SOL2
                              ELSE
                                 GOODSOL = 0
                              END IF 
                           END IF
                        ELSE
                           GOODSOL = 0
                        END IF
                     END IF

                     ! Now further checks for each of the good solutions:
                     ! - if the velocity at collision time is actually towards the surface
                     ! - if the collision point is within the extremes of the segment
                     ! Note: this requires moving the particle to the collision point and then moving back.
                     IF (GOODSOL /= 0) THEN
                        DO SOL = 1, GOODSOL
                           CALL MOVE_PARTICLE(IP, TEST(SOL))
                           IF ((particles(IP)%VX*NORMX(i) + particles(IP)%VY*NORMY(i)) < 0) THEN
                              ! Collision happens!
                              DTCOLL = TEST(SOL)
                              BOUNDCOLL = i    
                              TOTDTCOLL = TOTDTCOLL + DTCOLL  
                              HASCOLLIDED = .TRUE.
                              IF (BOUNDCOLL == 3) WRITE(*,*) "Collision with axis!"
                           END IF
                           CALL MOVE_PARTICLE(IP, -TEST(SOL))
                        END DO
                     END IF

                     

                  ELSE
                     ! We are not 2d axisymmetric
                     ! Compute the velocity normal to the boundary
                     VN = -(particles(IP)%VX * NORMX(i) + particles(IP)%VY * NORMY(i))
                     ! Compute the distance from the boundary
                     DX = (particles(IP)%X - XW(i)) * NORMX(i) + (particles(IP)%Y - YW(i)) * NORMY(i)
                     ! Check if a collision happens (sooner than previously calculated)
                     IF (VN .GE. 0. .AND. VN * DTCOLL .GT. DX) THEN
                        
                        ! Find candidate collision point (don't need this for boundaries)
                        ! XCOLL = particles(IP)%X + particles(IP)%VX * DTCOLL 
                        ! YCOLL = particles(IP)%Y + particles(IP)%VY * DTCOLL
                        ! ZCOLL = particles(IP)%Z + particles(IP)%VZ * DTCOLL

                        DTCOLL = DX/VN
                        BOUNDCOLL = i    
                        TOTDTCOLL = TOTDTCOLL + DTCOLL  
                        HASCOLLIDED = .TRUE.               

                     END IF
                  END IF
               END DO

               

               ! Check collisions with surfaces
               ! If earlier than dtcoll remember to set BOUNDCOLL to -1 and the new dtcoll
               WALLCOLL = -1
               DO i = 1, N_WALLS
                  IF (AXI) THEN
                     ! We are axisymmetric
                     CANDIDATE_DTCOLL = DTCOLL
                     ! Compute auxiliary parameters
                     IF (WALLS(i)%X1 == WALLS(i)%X2) THEN
                        ! Vertical wall
                        SOL1 = (WALLS(i)%X1-particles(IP)%X)/particles(IP)%VX
                        IF (SOL1 >= 0 .AND. SOL1 < DTCOLL) THEN
                           GOODSOL = 1
                           TEST(1) = SOL1
                        ELSE
                           GOODSOL = 0
                        END IF
                     ELSE
                        ! Non-vertical wall
                        ALPHA = (WALLS(i)%Y2-WALLS(i)%Y1)/(WALLS(i)%X2-WALLS(i)%X1)
                        BETA  = (particles(IP)%X - WALLS(i)%X1)*ALPHA

                        COEFA = particles(IP)%VY**2 + particles(IP)%VZ**2 - ALPHA*ALPHA*particles(IP)%VX**2
                        COEFB = particles(IP)%Y * particles(IP)%VY - (BETA + WALLS(i)%Y1)*ALPHA*particles(IP)%VX
                        COEFC = particles(IP)%Y**2 - (BETA + WALLS(i)%Y1)**2
                        DELTA = COEFB*COEFB-COEFA*COEFC
                        IF (DELTA .GE. 0) THEN
                           ! Compute the solutions, check if they are any good and, in case, order them to be checked further.
                           SOL1 = (-COEFB - SQRT(DELTA))/COEFA
                           SOL2 = (-COEFB + SQRT(DELTA))/COEFA
                           IF (SOL1 >= 0 .AND. SOL1 < DTCOLL .AND. BETA+ALPHA*particles(IP)%VX*SOL1+WALLS(i)%Y1 >= 0) THEN
                              IF (SOL2 >= 0 .AND. SOL2 < DTCOLL .AND. BETA+ALPHA*particles(IP)%VX*SOL2+WALLS(i)%Y1 >= 0) THEN
                                 GOODSOL = 2
                                 IF (SOL1 <= SOL2) THEN
                                    TEST(1) = SOL1
                                    TEST(2) = SOL2
                                 ELSE
                                    TEST(1) = SOL2
                                    TEST(2) = SOL1
                                 END IF
                              ELSE
                                 GOODSOL = 1
                                 TEST(1) = SOL1
                              ENDIF
                           ELSE
                              IF (SOL2 >= 0 .AND. SOL2 < DTCOLL .AND. BETA+ALPHA*particles(IP)%VX*SOL2+WALLS(i)%Y1 >= 0) THEN
                                 GOODSOL = 1
                                 TEST(1) = SOL2
                              ELSE
                                 GOODSOL = 0
                              END IF 
                           END IF
                        ELSE
                           GOODSOL = 0
                        END IF
                     END IF

                     ! Now further checks for each of the good solutions:
                     ! - if the velocity at collision time is actually towards the surface
                     ! - if the collision point is within the extremes of the segment
                     ! Note: this requires moving the particle to the collision point and then moving back.
                     IF (GOODSOL /= 0) THEN
                        DO SOL = 1, GOODSOL
                           CALL MOVE_PARTICLE(IP, TEST(SOL))
                           IF ((particles(IP)%VX*WALLS(i)%NORMX + particles(IP)%VY*WALLS(i)%NORMY) < 0) THEN

                              COLLDIST = ((particles(IP)%X-WALLS(i)%X1)*(WALLS(i)%X2-WALLS(i)%X1) &
                                       + (particles(IP)%Y-WALLS(i)%Y1)*(WALLS(i)%Y2-WALLS(i)%Y1)) &
                                       / ((WALLS(i)%X2-WALLS(i)%X1)**2 + (WALLS(i)%Y2-WALLS(i)%Y1)**2)
                              IF ((COLLDIST .GE. 0) .AND. (COLLDIST .LE. 1)) THEN
                                 ! Collision happens!
                                 DTCOLL = TEST(SOL)
                                 BOUNDCOLL = -1
                                 WALLCOLL = i   
                                 TOTDTCOLL = TOTDTCOLL + DTCOLL  
                                 HASCOLLIDED = .TRUE.
                              END IF
                           END IF
                           CALL MOVE_PARTICLE(IP, -TEST(SOL))
                        END DO
                     END IF


                  ELSE
                     ! We are not axisymmetric
                     ! Compute the velocity normal to the surface
                     VN = -(particles(IP)%VX * WALLS(i)%NORMX + particles(IP)%VY * WALLS(i)%NORMY)
                     ! Compute the distance from the boundary
                     DX = (particles(IP)%X - WALLS(i)%X1) * WALLS(i)%NORMX + (particles(IP)%Y - WALLS(i)%Y1) * WALLS(i)%NORMY
                     ! Check if a collision happens (sooner than previously calculated)
                     IF (DX .GE. 0. .AND. VN * DTCOLL .GE. DX) THEN
                        
                        CANDIDATE_DTCOLL = DX/VN
                        ! Find candidate collision point
                        XCOLL = particles(IP)%X + particles(IP)%VX * CANDIDATE_DTCOLL 
                        YCOLL = particles(IP)%Y + particles(IP)%VY * CANDIDATE_DTCOLL
                        COLLDIST = ((XCOLL-WALLS(i)%X1)*(WALLS(i)%X2-WALLS(i)%X1)+(YCOLL-WALLS(i)%Y1)*(WALLS(i)%Y2-WALLS(i)%Y1))/ &
                        ((WALLS(i)%X2-WALLS(i)%X1)**2 + (WALLS(i)%Y2-WALLS(i)%Y1)**2)
                        IF ((COLLDIST .GE. 0) .AND. (COLLDIST .LE. 1)) THEN
                           BOUNDCOLL = -1
                           WALLCOLL = i
                           DTCOLL = CANDIDATE_DTCOLL
                           TOTDTCOLL = TOTDTCOLL + DTCOLL
                           HASCOLLIDED = .TRUE.
                        END IF
                     END IF
                  END IF
               END DO


               IF (BOUNDCOLL .NE. -1) THEN
                  ! Collision with a boundary
                  ! Tally the collision
                  LOCAL_BOUNDARY_COLL_COUNT(BOUNDCOLL+4*(particles(IP)%S_ID-1)) = &
                  LOCAL_BOUNDARY_COLL_COUNT(BOUNDCOLL+4*(particles(IP)%S_ID-1)) + 1

                  IF (BOOL_PERIODIC(BOUNDCOLL)) THEN

                     CALL MOVE_PARTICLE(IP, DTCOLL)
                     SELECT CASE (BOUNDCOLL)
                        CASE (1)
                           particles(IP)%X = particles(IP)%X + XMAX - XMIN
                        CASE (2)
                           particles(IP)%X = particles(IP)%X - XMAX + XMIN
                        CASE (3)
                           particles(IP)%Y = particles(IP)%Y + YMAX - YMIN
                        CASE (4)
                           particles(IP)%Y = particles(IP)%Y - YMAX + YMIN
                     END SELECT
                     particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

                  ELSE IF (BOOL_SPECULAR(BOUNDCOLL)) THEN
                  
                     CALL MOVE_PARTICLE(IP, DTCOLL)
                     SELECT CASE (BOUNDCOLL)
                        CASE (1)
                           particles(IP)%VX = - particles(IP)%VX
                        CASE (2)
                           particles(IP)%VX = - particles(IP)%VX
                        CASE (3)
                           particles(IP)%VY = - particles(IP)%VY
                        CASE (4)
                           particles(IP)%VY = - particles(IP)%VY
                     END SELECT
                     particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

                     IF (BOOL_REACT(BOUNDCOLL)) THEN
                        CALL WALL_REACT(particles, IP, REMOVE_PART(IP))
                     END IF

                  ELSE IF (BOOL_DIFFUSE(BOUNDCOLL)) THEN

                     CALL MOVE_PARTICLE(IP, DTCOLL)
                     particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
                     IF (BOOL_REACT(BOUNDCOLL)) THEN
                        CALL WALL_REACT(particles, IP, REMOVE_PART(IP))
                     END IF
                     IF (.NOT. REMOVE_PART(IP)) THEN
                        S_ID = particles(IP)%S_ID
                        CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                        BOUNDTEMP, BOUNDTEMP, BOUNDTEMP, &
                        VDUMMY, V_PERP, VZ, SPECIES(S_ID)%MOLECULAR_MASS)

                        CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, BOUNDTEMP, EROT)
                        CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, BOUNDTEMP, EVIB)
                                       
                        V_NORM = FLX(0.d0, BOUNDTEMP, SPECIES(S_ID)%MOLECULAR_MASS)

                        particles(IP)%VX = V_NORM*NORMX(BOUNDCOLL) - V_PERP*NORMY(BOUNDCOLL)
                        particles(IP)%VY = V_PERP*NORMX(BOUNDCOLL) + V_NORM*NORMY(BOUNDCOLL)
                        particles(IP)%VZ = VZ
                        particles(IP)%EROT = EROT
                        particles(IP)%EVIB = EVIB

                     END IF

                  ELSE

                     REMOVE_PART(IP) = .TRUE.
                     particles(IP)%DTRIM = 0.

                  END IF
               ELSE IF (WALLCOLL .NE. -1) THEN
                  ! Collision with a wall
                  ! Tally the collision
                  LOCAL_WALL_COLL_COUNT(WALLCOLL+N_WALLS*(particles(IP)%S_ID-1)) = &
                  LOCAL_WALL_COLL_COUNT(WALLCOLL+N_WALLS*(particles(IP)%S_ID-1)) + 1

                  IF (WALLS(WALLCOLL)%SPECULAR) THEN
                     rfp = rf()
                     IF (WALLS(WALLCOLL)%POROUS .AND. rfp .LE. WALLS(WALLCOLL)%TRANSMISSIVITY) THEN

                        REMOVE_PART(IP) = .TRUE.
                        particles(IP)%DTRIM = 0.

                     ELSE

                        CALL MOVE_PARTICLE(IP, DTCOLL)
                        !WRITE(*,*) "Specular reflection! Velocity before: ",  particles(IP)%VX, ", ", particles(IP)%VY
                        VDOTN = particles(IP)%VX*WALLS(WALLCOLL)%NORMX+particles(IP)%VY*WALLS(WALLCOLL)%NORMY
                        particles(IP)%VX = particles(IP)%VX - 2.*VDOTN*WALLS(WALLCOLL)%NORMX
                        particles(IP)%VY = particles(IP)%VY - 2.*VDOTN*WALLS(WALLCOLL)%NORMY
                        !WRITE(*,*) "Normal:  ",  WALLS(WALLCOLL)%NORMX, ", ", WALLS(WALLCOLL)%NORMY
                        !WRITE(*,*) "Specular reflection! Velocity after:  ",  particles(IP)%VX, ", ", particles(IP)%VY
                        particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

                        IF (WALLS(WALLCOLL)%REACT) THEN
                           CALL WALL_REACT(particles, IP, REMOVE_PART(IP))
                        END IF

                     END IF
                  ELSE IF (WALLS(WALLCOLL)%DIFFUSE) THEN
                     rfp = rf()
                     IF (WALLS(WALLCOLL)%POROUS .AND. rfp .LE. WALLS(WALLCOLL)%TRANSMISSIVITY) THEN

                        REMOVE_PART(IP) = .TRUE.
                        particles(IP)%DTRIM = 0.
                        
                     ELSE
                           
                        CALL MOVE_PARTICLE(IP, DTCOLL)
                        particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
                        IF (WALLS(WALLCOLL)%REACT) THEN
                           CALL WALL_REACT(particles, IP, REMOVE_PART(IP))
                        END IF
                        IF (.NOT. REMOVE_PART(IP)) THEN
                           S_ID = particles(IP)%S_ID
                           CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                           WALLS(WALLCOLL)%TEMP, WALLS(WALLCOLL)%TEMP, WALLS(WALLCOLL)%TEMP, &
                           VDUMMY, V_PERP, VZ, SPECIES(S_ID)%MOLECULAR_MASS)

                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, WALLS(WALLCOLL)%TEMP, EROT)
                           CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, WALLS(WALLCOLL)%TEMP, EVIB)
                                          
                           V_NORM = FLX(0.d0, WALLS(WALLCOLL)%TEMP, SPECIES(S_ID)%MOLECULAR_MASS)

                           particles(IP)%VX = V_NORM*WALLS(WALLCOLL)%NORMX - V_PERP*WALLS(WALLCOLL)%NORMY
                           particles(IP)%VY = V_PERP*WALLS(WALLCOLL)%NORMX + V_NORM*WALLS(WALLCOLL)%NORMY
                           particles(IP)%VZ = VZ
                           particles(IP)%EROT = EROT
                           particles(IP)%EVIB = EVIB

                        END IF
                     END IF

                  ELSE

                     REMOVE_PART(IP) = .TRUE.
                     particles(IP)%DTRIM = 0.

                  END IF

               ELSE
                  ! No collision in this timestep
                  CALL MOVE_PARTICLE(IP, particles(IP)%DTRIM)
                  particles(IP)%DTRIM = 0.

               END IF


               ! Axisymmetric velocity rotation already performed in MOVE_PARTICLE()

         
               ! _______ CHECK WHERE PARTICLE ENDED UP _______

               ! ++++++++ Check if particle is still in the domain ++++++++++++
               !IF (particles(IP)%X .GE. XMIN .AND. particles(IP)%X .LE. XMAX .AND. & 
               !    particles(IP)%Y .GE. YMIN .AND. particles(IP)%Y .LE. YMAX .AND. &
               !    particles(IP)%Z .GE. ZMIN .AND. particles(IP)%Z .LE. ZMAX) THEN   ! Check that the particle is still in the domain
         
                  ! Compute the index of the cell in which the particle ended up
                  !CALL CELL_FROM_POSITION(particles(IP)%X, particles(IP)%Y, IC)
                  !particles(IP)%IC = IC
         
                  !DTRIM = 0.E0 ! Timestep is over.
      
               ! +++++++++ Particle crossed domain boundaries ++++++++
               !ELSE

                  ! CALL IMPACT_BOUNDARY(IP, DTRIM) 
               !   DTRIM = 0.E0 ! TMP TMP TMP TMP TMP TMP TMP TMP TMP
               !END IF

            END IF

         END DO ! WHILE (DTRIM .GT. 0.)

            ! _______ APPLY Z PERIODIC BCs ________

            ! Z we should be able to do this a posteriori (remember it could be more than one depth length out!)
         IF (BOOL_Z_PERIODIC) THEN
         
            DO WHILE (particles(IP)%Z .GT. ZMAX)
               particles(IP)%Z = ZMIN + (particles(IP)%Z - ZMAX)
            END DO
            DO WHILE (particles(IP)%Z .LT. ZMIN) 
               particles(IP)%Z = ZMAX + (particles(IP)%Z - ZMIN)
            END DO
      
         END IF

         IF (DIMS == 1 .AND. BOOL_Y_PERIODIC) THEN
         
            DO WHILE (particles(IP)%Y .GT. YMAX)
               particles(IP)%Y = YMIN + (particles(IP)%Y - YMAX)
            END DO
            DO WHILE (particles(IP)%Y .LT. YMIN) 
               particles(IP)%Y = YMAX + (particles(IP)%Y - YMIN)
            END DO
      
         END IF 

         particles(IP)%DTRIM = DT ! For the next timestep.

         !IF (HASCOLLIDED) THEN
         !   V_OLD(1) = particles(IP)%VX
         !   V_OLD(2) = particles(IP)%VY
         !   V_OLD(3) = particles(IP)%VZ
            !Propagate v_coll to v^{n+1/2}
            !Consistent: does not change for neutral particles.
            !Also does not change is collision happens at DT/2
            !CALL UPDATE_VELOCITY_BORIS(DT/2., V_OLD, V_NEW, &
            !SPECIES(particles(IP)%S_ID)%CHARGE, SPECIES(particles(IP)%S_ID)%MOLECULAR_MASS, &
            !E, B)
         !END IF
         !particles(IP)%VX = V_NEW(1)
         !particles(IP)%VY = V_NEW(2)
         !particles(IP)%VZ = V_NEW(3)


      END DO ! End loop: DO IP = 1,NP_PROC


      

      ! ==============
      ! ========= Now remove particles that are out of the domain and reorder the array. 
      ! ========= Do this in the end of the advection step, since reordering the array 
      ! ========= mixes the particles.
      ! ==============

      IP = NP_PROC
      DO WHILE (IP .GE. 1)

         ! Is particle IP out of the domain? Then remove it!
         IF (REMOVE_PART(IP)) THEN
            CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)
         ELSE IF (GRID_TYPE .NE. UNSTRUCTURED) THEN
            CALL CELL_FROM_POSITION(particles(IP)%X, particles(IP)%Y, IC)
            OLD_IC = particles(IP)%IC
            particles(IP)%IC = IC
            ! If cell-based particle weight has changed,
            IF (BOOL_RADIAL_WEIGHTING .AND. (IC .NE. OLD_IC)) THEN
               WEIGHT_RATIO = CELL_FNUM(OLD_IC)/CELL_FNUM(IC)
               ! If the weight has increased, the particle may be removed
               IF (WEIGHT_RATIO .LT. 1.d0) THEN
                  rfp = rf()
                  IF (WEIGHT_RATIO .LT. rfp) CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)
               ELSE IF (WEIGHT_RATIO .GT. 1.d0) THEN
               ! If the weight has decreased, more particles may be added
                  rfp = rf()
                  DO WHILE ((WEIGHT_RATIO - 1.) .GT. rfp)
                     CALL INIT_PARTICLE(particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
                     particles(IP)%VX, particles(IP)%VY,particles(IP)%VZ, &
                     particles(IP)%EROT,particles(IP)%EVIB,particles(IP)%S_ID,IC,DT,NEWparticle)
                     CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)
                     WEIGHT_RATIO = WEIGHT_RATIO - 1.
                  END DO
               END IF
            END IF
         END IF

         IP = IP - 1

      END DO

      DEALLOCATE(REMOVE_PART)


      ! Tally surface collisions from all processes
      CALL MPI_REDUCE(LOCAL_BOUNDARY_COLL_COUNT, BOUNDARY_COLL_COUNT, 4*N_SPECIES, MPI_INTEGER, MPI_SUM, 0, &
      MPI_COMM_WORLD, ierr)

      DEALLOCATE(LOCAL_BOUNDARY_COLL_COUNT)
      
      CALL MPI_REDUCE(LOCAL_WALL_COLL_COUNT, WALL_COLL_COUNT, N_WALLS*N_SPECIES, MPI_INTEGER, MPI_SUM, 0, &
      MPI_COMM_WORLD, ierr)

      DEALLOCATE(LOCAL_WALL_COLL_COUNT)


      IF ((tID .GT. DUMP_PART_BOUND_START) .AND. (tID .NE. RESTART_TIMESTEP)) THEN
         IF (MOD(tID-DUMP_PART_BOUND_START, DUMP_PART_BOUND_EVERY) .EQ. 0) THEN
            CALL DUMP_BOUNDARY_PARTICLES_FILE(tID)
            IF (ALLOCATED(part_dump)) DEALLOCATE(part_dump)
            NP_DUMP_PROC = 0
         END IF
      END IF


      !CLOSE(66341)

   END SUBROUTINE ADVECT


   SUBROUTINE UPDATE_VELOCITY_BORIS(DTIME, V_OLD, V_NEW, CHARGE, MASS, E, B)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: V_NEW
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: V_OLD, E, B
      REAL(KIND=8), DIMENSION(3) :: V_MINUS, V_PLUS, V_PRIME, T, S
      REAL(KIND=8), INTENT(IN) :: DTIME, CHARGE, MASS
      REAL(KIND=8) :: COULOMBCHARGE

      COULOMBCHARGE = CHARGE * QE
      V_MINUS = V_OLD + 0.5*COULOMBCHARGE*E/MASS*DTIME

      T = 0.5*COULOMBCHARGE*B/MASS*DTIME
      V_PRIME = V_MINUS + CROSS(V_MINUS, T)
      S = 2.*T/(1.+( T(1)*T(1) + T(2)*T(2) + T(3)*T(3) ))
      V_PLUS = V_MINUS + CROSS(V_PRIME, S)

      V_NEW = V_PLUS + 0.5*COULOMBCHARGE*E/MASS*DTIME

   END SUBROUTINE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MOVE_PARTICLE -> Move the particles !!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE MOVE_PARTICLE(IP, TIME)

      ! Moves particle with index IP for time TIME.
      ! For now simply rectilinear movement, will have to include Lorentz's force

      IMPLICIT NONE

      INTEGER, INTENT(IN)      :: IP
      REAL(KIND=8), INTENT(IN) :: TIME
      REAL(KIND=8)             :: SINTHETA, COSTHETA, VZ, VY, R

      ! WRITE(*,*) "Moving particle ", IP, " for time interval ", TIME

      IF (AXI) THEN
         IF (particles(IP)%VZ == 0.d0) WRITE(*,*) 'Particle has exactly zero vz velocity!'
         particles(IP)%X = particles(IP)%X + particles(IP)%VX * TIME
         particles(IP)%Z = 0.d0
         R = SQRT((particles(IP)%Y + particles(IP)%VY*TIME)**2 + (particles(IP)%VZ*TIME)**2)
         ! Rotate velocity vector back to x-y plane.
         SINTHETA = particles(IP)%VZ*TIME / R
         COSTHETA = SIGN(SQRT(1.-SINTHETA*SINTHETA), particles(IP)%Y + particles(IP)%VY*TIME)

         particles(IP)%Y = R

         VZ = particles(IP)%VZ
         VY = particles(IP)%VY
         particles(IP)%VZ = COSTHETA*VZ - SINTHETA*VY
         particles(IP)%VY = SINTHETA*VZ + COSTHETA*VY
         
      ELSE
         particles(IP)%X = particles(IP)%X + particles(IP)%VX * TIME
         particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * TIME
         particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * TIME
      END IF

   END SUBROUTINE MOVE_PARTICLE







   SUBROUTINE REMOVE_PARTICLES_OF_SPECIES(SP_ID)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: SP_ID
      INTEGER :: IP

      IP = NP_PROC
      DO WHILE (IP .GE. 1)
         IF (particles(IP)%S_ID == SP_ID) CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)
         IP = IP - 1
      END DO

   END SUBROUTINE REMOVE_PARTICLES_OF_SPECIES


   SUBROUTINE REMOVE_PARTICLES_IN_MIXTURE(MIX_ID)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MIX_ID
      INTEGER :: I

      DO I = 1, MIXTURES(MIX_ID)%N_COMPONENTS
         CALL REMOVE_PARTICLES_OF_SPECIES(MIXTURES(MIX_ID)%COMPONENTS(I)%ID)
      END DO

   END SUBROUTINE REMOVE_PARTICLES_IN_MIXTURE


END MODULE timecycle
