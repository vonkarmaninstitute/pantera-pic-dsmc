MODULE timecycle

   USE global
   USE mpi_common
   USE screen
   USE tools
   USE collisions
   USE postprocess
   USE fields

   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE TIME_LOOP                                     !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE TIME_LOOP

   IMPLICIT NONE
 
   INTEGER :: NP_TOT, NCOLL_TOT
   REAL(KIND=8) :: CURRENT_TIME, CURRENT_CPU_TIME

   CHARACTER(len=512) :: stringTMP

   ! Init variables
   NP_TOT = 0
   CALL INIT_POSTPROCESS

   ! ########### Compute poisson ##########################################

   IF (BOOL_PIC) THEN
      CALL DEPOSIT_CHARGE
      CALL SOLVE_POISSON
   END IF
   
   ! Dump particles and flowfield before the first time step, but after the initial seeding
   IF (DUMP_START .EQ. 0) THEN
      CALL DUMP_PARTICLES_FILE(0)
   END IF

   IF (DUMP_GRID_START .EQ. 0 .AND. DUMP_GRID_N_AVG .EQ. 1) THEN
      CALL GRID_AVG
      CALL GRID_SAVE
      CALL GRID_RESET
   END IF


   tID = 1
   CALL CPU_TIME(START_CPU_TIME)
   DO WHILE (tID .LE. NT)

      ! ########### Print simulation info #######################################

      CURRENT_TIME = tID*DT
      CALL CPU_TIME(CURRENT_CPU_TIME)
      CURRENT_CPU_TIME = CURRENT_CPU_TIME - START_CPU_TIME 

      IF (MOD(tID, STATS_EVERY) .EQ. 0) THEN
         CALL MPI_REDUCE(NP_PROC, NP_TOT, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_COLL, NCOLL_TOT, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

         WRITE(stringTMP, '(A13,I8,A4,I8,A9,ES14.3,A17,F10.1,A28,I10,A25,I10)') '   Timestep: ', tID, ' of ', NT, &
                          ' - time: ', CURRENT_TIME, ' [s] - CPU time: ', CURRENT_CPU_TIME, &
                          ' [s] - number of particles: ', NP_TOT, &
                          ' - number of collisions: ', NCOLL_TOT

         CALL ONLYMASTERPRINT1(PROC_ID, TRIM(stringTMP))
      END IF
      
      ! ########### Inject particles from boundaries/injection sources ##########

      CALL BOUNDARIES_INJECT
      CALL LINE_SOURCE_INJECT

      ! ########### Advect particles ############################################

      CALL ADVECT 

      ! ########### Exchange particles among processes ##########################

      CALL EXCHANGE

      ! ########### Perform collisions ##########################################

      IF (BOOL_MCC)  CALL MCC_COLLISIONS

      IF (BOOL_DSMC) CALL DSMC_COLLISIONS


      IF (BOOL_THERMAL_BATH) CALL THERMAL_BATH

      ! ########### Compute poisson ##########################################

      IF (BOOL_PIC) THEN
         CALL DEPOSIT_CHARGE
         CALL SOLVE_POISSON
      END IF

      ! ########### Dump particles ##############################################
      IF (tID .GT. DUMP_START) THEN
         IF (MOD(tID-DUMP_START, DUMP_EVERY) .EQ. 0) CALL DUMP_PARTICLES_FILE(tID)
      END IF
      ! ########### Dump flowfield ##############################################

      IF (tID .GT. DUMP_GRID_START) THEN
         ! If we are in the grid save timestep, average, then dump the cumulated averages
         IF (tID .NE. DUMP_GRID_START .AND. MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY*DUMP_GRID_N_AVG) .EQ. 0) THEN
            CALL GRID_AVG
            CALL GRID_SAVE
            CALL GRID_RESET
         ! If we are just in a grid average timestep, compute the grid average
         ELSE IF (MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY) .EQ. 0) THEN
            CALL GRID_AVG
         END IF
      END IF
      ! ~~~~~ Hmm that's it! ~~~~~

      ! Perform the conservation checks
      IF (PERFORM_CHECKS .AND. MOD(tID, CHECKS_EVERY) .EQ. 0) CALL CHECKS

      tID = tID + 1
   END DO

   END SUBROUTINE TIME_LOOP
 


   SUBROUTINE LINE_SOURCE_INJECT
  
      IMPLICIT NONE
   
      INTEGER      :: IP, IC, IS, NFS, ILINE
      REAL(KIND=8) :: DTFRAC, Vdummy, V_NORM, V_PERP, POS
      REAL(KIND=8) :: X, Y, Z, VX, VY, VZ, EROT, EVIB 
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
   
      INTEGER :: S_ID
      REAL(KIND=8) :: M
   
   
      ! Linesource
      DO ILINE = 1, N_LINESOURCES
   
         DO IS = 1, MIXTURES(LINESOURCES(ILINE)%MIX_ID)%N_COMPONENTS ! Loop on mixture components

            S_ID = MIXTURES(LINESOURCES(ILINE)%MIX_ID)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            NFS = FLOOR(LINESOURCES(ILINE)%nfs(IS))
            IF (LINESOURCES(ILINE)%nfs(IS)-REAL(NFS, KIND=8) .GE. rf()) THEN ! Same as SPARTA's perspeciess
               NFS = NFS + 1
            END IF
   
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
               POS = rf()
               X = LINESOURCES(ILINE)%CX + (0.5-POS)*LINESOURCES(ILINE)%DX
               Y = LINESOURCES(ILINE)%CY + (0.5-POS)*LINESOURCES(ILINE)%DY
               Z = ZMIN + (ZMAX-ZMIN)*rf()
   
               CALL CELL_FROM_POSITION(X,Y,  IC)
   
               ! Init a particle object and assign it to the local vector of particles
               CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTFRAC,  particleNOW)
               CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles)
   
            END DO
         END DO
   
      END DO
         
   END SUBROUTINE LINE_SOURCE_INJECT


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE BOUNDARIES_INJECT -> Injects particles from domain borders !!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
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

   INTEGER      :: IP, IC, i
   INTEGER      :: BOUNDCOLL, WALLCOLL
   REAL(KIND=8) :: DTCOLL, TOTDTCOLL, CANDIDATE_DTCOLL
   REAL(KIND=8), DIMENSION(4) :: NORMX, NORMY, XW, YW, ZW
   ! REAL(KIND=8) :: XCOLL, YCOLL, ZCOLL
   REAL(KIND=8) :: VN, DX
   LOGICAL, DIMENSION(:), ALLOCATABLE :: REMOVE_PART
   REAL(KIND=8), DIMENSION(3) :: V_OLD, V_NEW
   REAL(KIND=8), DIMENSION(3) :: E, B
   REAL(KIND=8) :: V_NORM, V_PERP, VZ, VDUMMY, EROT, EVIB, VDOTN
   INTEGER :: S_ID
   LOGICAL :: HASCOLLIDED
   REAL(KIND=8) :: XCOLL, YCOLL, COLLDIST

   !E = [0.d0, 0.d0, 0.d0]
   B = [0.d0, 0.d0, 0.d0]


   NORMX = [ 1., -1., 0., 0. ]
   NORMY = [ 0., 0., 1., -1. ]

   XW = [ XMIN, XMAX, XMIN, XMAX ]
   YW = [ YMAX, YMIN, YMIN, YMAX ]
   ZW = [ 0., 0., 0., 0. ]

   ALLOCATE(REMOVE_PART(NP_PROC))

   DO IP = 1, NP_PROC
      REMOVE_PART(IP) = .FALSE.


      ! Update velocity
      IC = particles(IP)%IC

      IF (BOOL_PIC) THEN
         CALL APPLY_E_FIELD(IP, E)
         


         V_OLD(1) = particles(IP)%VX
         V_OLD(2) = particles(IP)%VY
         V_OLD(3) = particles(IP)%VZ

         CALL UPDATE_VELOCITY_BORIS(DT, V_OLD, V_NEW, &
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


      HASCOLLIDED = .FALSE.
      TOTDTCOLL = 0.
      DO WHILE (particles(IP)%DTRIM .GT. 0.) ! Repeat the procedure until step is done

         DTCOLL = particles(IP)%DTRIM ! Looking for collisions within the remaining time
         ! ______ ADVECTION ______
         BOUNDCOLL = -1
         DO i = 1, 4 ! Check collisions with boundaries (xmin, xmax, ymin, ymax)
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
         END DO


         ! Check collisions with surfaces
         ! If earlier than dtcoll remember to set BOUNDCOLL to -1 and the new dtcoll
         WALLCOLL = -1
         DO i = 1, N_WALLS
            ! Compute the velocity normal to the surface
            VN = -(particles(IP)%VX * WALLS(i)%NORMX + particles(IP)%VY * WALLS(i)%NORMY)
            ! Compute the distance from the boundary
            DX = (particles(IP)%X - WALLS(i)%CX) * WALLS(i)%NORMX + (particles(IP)%Y - WALLS(i)%CY) * WALLS(i)%NORMY
            ! Check if a collision happens (sooner than previously calculated)
            IF (DX .GE. 0. .AND. VN * DTCOLL .GE. DX) THEN
               
               CANDIDATE_DTCOLL = DX/VN
               ! Find candidate collision point
               XCOLL = particles(IP)%X + particles(IP)%VX * CANDIDATE_DTCOLL 
               YCOLL = particles(IP)%Y + particles(IP)%VY * CANDIDATE_DTCOLL
               COLLDIST = (XCOLL-WALLS(i)%CX)**2+(YCOLL-WALLS(i)%CY)**2
               IF (COLLDIST .LE. (WALLS(i)%DX**2+WALLS(i)%DY**2)/4.) THEN
                  BOUNDCOLL = -1
                  WALLCOLL = i
                  DTCOLL = CANDIDATE_DTCOLL
                  TOTDTCOLL = TOTDTCOLL + DTCOLL
                  HASCOLLIDED = .TRUE.
               END IF
            END IF
         END DO


         IF (BOUNDCOLL .NE. -1) THEN
            ! Collision with a boundary
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
                  CALL WALL_REACT(IP, REMOVE_PART(IP))
               END IF

            ELSE IF (BOOL_DIFFUSE(BOUNDCOLL)) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
               IF (BOOL_REACT(BOUNDCOLL)) THEN
                  CALL WALL_REACT(IP, REMOVE_PART(IP))
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

            IF (WALLS(WALLCOLL)%SPECULAR) THEN

               IF (WALLS(WALLCOLL)%POROUS .AND. rf() .LE. WALLS(WALLCOLL)%TRANSMISSIVITY) THEN

                  REMOVE_PART(IP) = .TRUE.
                  particles(IP)%DTRIM = 0.

               ELSE

                  CALL MOVE_PARTICLE(IP, DTCOLL)
                  VDOTN = particles(IP)%VX*WALLS(WALLCOLL)%NORMX+particles(IP)%VY*WALLS(WALLCOLL)%NORMY
                  particles(IP)%VX = particles(IP)%VX + 2.*VDOTN*WALLS(WALLCOLL)%NORMX
                  particles(IP)%VY = particles(IP)%VY + 2.*VDOTN*WALLS(WALLCOLL)%NORMY
                  
                  particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

                  IF (WALLS(WALLCOLL)%REACT) THEN
                     CALL WALL_REACT(IP, REMOVE_PART(IP))
                  END IF

               END IF
            ELSE IF (WALLS(WALLCOLL)%DIFFUSE) THEN

               IF (WALLS(WALLCOLL)%POROUS .AND. rf() .LE. WALLS(WALLCOLL)%TRANSMISSIVITY) THEN

                  REMOVE_PART(IP) = .TRUE.
                  particles(IP)%DTRIM = 0.
                  
               ELSE
                     
                  CALL MOVE_PARTICLE(IP, DTCOLL)
                  particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
                  IF (WALLS(WALLCOLL)%REACT) THEN
                     CALL WALL_REACT(IP, REMOVE_PART(IP))
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


      

         ! Axisymmetric domain
         IF (BOOL_AXI .eqv. .TRUE.) THEN 
   
            ! Transform position, velocities etc
   
         END IF

   
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
      ELSE
         CALL CELL_FROM_POSITION(particles(IP)%X, particles(IP)%Y, IC)
         particles(IP)%IC = IC
      END IF

      IP = IP - 1

   END DO

   DEALLOCATE(REMOVE_PART)

   END SUBROUTINE ADVECT


   SUBROUTINE UPDATE_VELOCITY_BORIS(DT, V_OLD, V_NEW, CHARGE, MASS, E, B)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: V_NEW
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: V_OLD, E, B
      REAL(KIND=8), DIMENSION(3) :: V_MINUS, V_PLUS, V_PRIME, T, S
      REAL(KIND=8), INTENT(IN) :: DT, CHARGE, MASS
      REAL(KIND=8) :: COULOMBCHARGE

      COULOMBCHARGE = CHARGE * 1.602176634e-19
      V_MINUS = V_OLD + 0.5*COULOMBCHARGE*E/MASS*DT

      T = 0.5*COULOMBCHARGE*B/MASS*DT
      V_PRIME = V_MINUS + CROSS(V_MINUS, T)
      S = 2.*T/(1.+( T(1)*T(1) + T(2)*T(2) + T(3)*T(3) ))
      V_PLUS = V_MINUS + CROSS(V_PRIME, S)

      V_NEW = V_PLUS + 0.5*COULOMBCHARGE*E/MASS*DT

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

      particles(IP)%X = particles(IP)%X + particles(IP)%VX * TIME
      particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * TIME
      particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * TIME

   END SUBROUTINE MOVE_PARTICLE


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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MCC_COLLISIONS -> perform MCC collisions !!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE MCC_COLLISIONS

   ! Performs MCC collisions: for each particle in the current process ID, computes the probability of 
   ! colliding with a background species, and (in case) performs the collision.
   ! Collisions are isotropic in the center of mass frame, and assume that the current particle 
   ! collides with a much heavier one (electron-neutral collisions). 
   ! 
   ! There is no need to consider cells for these collisions. The neutrals density is for the
   ! moment uniform in space. In the future, we may associate it to the cells...
   
   IMPLICIT NONE

   INTEGER      :: IP
   REAL(KIND=8) :: V_NOW, P_COLL, CHI, COS_TH, SIN_TH
   REAL(KIND=8) :: PI,PI2

   PI  = 3.141593
   PI2 = 2*PI

   DO IP = 1,NP_PROC

      ! particle velocity
      V_NOW = SQRT(particles(IP)%VX**2 + particles(IP)%VY**2+particles(IP)%VZ**2)

      ! Compute collision probability
      P_COLL = 1 - exp(-DT*MCC_BG_DENS*MCC_SIGMA*V_NOW)

      ! Try the collision
      IF (rf() < P_COLL) THEN ! Collision happens

         ! Sample angles randomly (uniform distribution on a sphere)
         CHI    = PI2*rf()
         COS_TH = 2*rf() - 1
         SIN_TH = SQRT(1 - COS_TH**2)

         ! And compute new velocity (this should be in the center-of-mass frame, but m_electron << m_neutral)
         particles(IP)%VX = V_NOW*SIN_TH*COS(CHI)
         particles(IP)%VY = V_NOW*SIN_TH*SIN(CHI)
         particles(IP)%VZ = V_NOW*COS_TH

      END IF

   END DO

   END SUBROUTINE MCC_COLLISIONS

END MODULE timecycle
