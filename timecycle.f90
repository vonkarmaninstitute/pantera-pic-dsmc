MODULE timecycle

   USE global
   USE mpi_common
   USE screen
   USE tools
   USE collisions
   USE postprocess

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
   
   ! Dump particles before the first time step, but after the initial seeding
   CALL DUMP_PARTICLES_FILE(0)
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

      ! ########### Dump particles ##############################################

      IF (MOD(tID, DUMP_EVERY) .EQ. 0) CALL DUMP_PARTICLES_FILE(tID)

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
      IF (PERFORM_CHECKS .AND. MOD(tID, 100) .EQ. 0) CALL CHECKS

      tID = tID + 1
   END DO

   END SUBROUTINE TIME_LOOP
 


   SUBROUTINE LINE_SOURCE_INJECT
  
      IMPLICIT NONE
   
      INTEGER      :: IP, IC, IS, NFS, ILINE
      REAL(KIND=8) :: DTFRAC, Vdummy, V_NORM, V_PERP, POS
      REAL(KIND=8) :: X, Y, Z, VX, VY, VZ, EI 
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
   
      INTEGER :: S_ID
      REAL(KIND=8) :: M
      REAL(KIND=8) :: ZERO = 0.0
   
   
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
   
               CALL MAXWELL(ZERO, ZERO, ZERO, &
                            LINESOURCES(ILINE)%TTRA, LINESOURCES(ILINE)%TTRA, LINESOURCES(ILINE)%TTRA, &
                            LINESOURCES(ILINE)%TROT, &
                            Vdummy, V_PERP, VZ, EI, M)
                            
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
               CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EI,S_ID,IC,DTFRAC,  particleNOW)
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
   REAL(KIND=8) :: X, Y, Z, VX, VY, VZ, EI 
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
                        TTRA_BOUND,TTRA_BOUND,TTRA_BOUND,TROT_BOUND, &
                        Vdummy,VY,VZ,EI,M)
            
            VX = UX_BOUND + FLX(S_NORM_XMIN,TTRA_BOUND,M) !!AAAAAAAAAAA

            DTFRAC = rf()*DT
            X = XMIN
            Y = YMIN + (YMAX-YMIN)*rf()
            Z = ZMIN + (ZMAX-ZMIN)*rf()

            CALL CELL_FROM_POSITION(X,Y,  IC)

            ! Init a particle object and assign it to the local vector of particles
            CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EI,S_ID,IC,DTFRAC,  particleNOW)
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
                        TTRA_BOUND,TTRA_BOUND,TTRA_BOUND,TROT_BOUND, &
                        Vdummy,VY,VZ,EI,M)

            VX = UX_BOUND - FLX(S_NORM_XMAX,TTRA_BOUND,M) !! I think something was wrong here with the sign

            DTFRAC = rf()*DT
            X = XMAX
            Y = YMIN + (YMAX-YMIN)*rf()
            Z = ZMIN + (ZMAX-ZMIN)*rf()
            
            CALL CELL_FROM_POSITION(X,Y,  IC)

            ! Init a particle object and assign it to vector of particles
            CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EI,S_ID,IC,DTFRAC,  particleNOW)
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
                        TTRA_BOUND,TTRA_BOUND,TTRA_BOUND,TROT_BOUND, &
                        VX,Vdummy,VZ,EI,M)

            VY = UY_BOUND + FLX(S_NORM_YMIN,TTRA_BOUND,M) !!AAAAAAAAAAA

            DTFRAC = rf()*DT
            X = XMIN + (XMAX-XMIN)*rf()
            Y = YMIN
            Z = ZMIN + (ZMAX-ZMIN)*rf()

            CALL CELL_FROM_POSITION(X,Y,  IC)

            ! Init a particle object and assign it to the local vector of particles
            CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EI,S_ID,IC,DTFRAC,  particleNOW)
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
                        TTRA_BOUND,TTRA_BOUND,TTRA_BOUND,TROT_BOUND, &
                        VX,Vdummy,VZ,EI,M)

            VY = UY_BOUND - FLX(S_NORM_YMAX,TTRA_BOUND,M) !!AAAAAAAAAAA

            DTFRAC = rf()*DT
            X = XMIN + (XMAX-XMIN)*rf() ! There was a bug here!
            Y = YMAX
            Z = ZMIN + (ZMAX-ZMIN)*rf()

            CALL CELL_FROM_POSITION(X,Y,  IC)

            ! Init a particle object and assign it to vector of particles
            CALL INIT_PARTICLE(X,Y,Z,VX,VY,VZ,EI,S_ID,IC,DTFRAC,  particleNOW)
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
   INTEGER      :: BOUNDCOLL
   REAL(KIND=8) :: DTCOLL
   REAL(KIND=8), DIMENSION(4) :: NX, NY, NZ, XW, YW, ZW
   ! REAL(KIND=8) :: XCOLL, YCOLL, ZCOLL
   REAL(KIND=8) :: VN, DX
   LOGICAL, DIMENSION(:), ALLOCATABLE :: REMOVE_PART

   NX = (/ 1., -1., 0., 0. /)
   NY = (/ 0., 0., 1., -1. /)
   NZ = (/ 0., 0., 0., 0. /)

   XW = (/ XMIN, XMAX, XMIN, XMAX /)
   YW = (/ YMAX, YMIN, YMIN, YMAX /)
   ZW = (/ 0., 0., 0., 0. /)

   ALLOCATE(REMOVE_PART(NP_PROC))

   DO IP = 1, NP_PROC
      REMOVE_PART(IP) = .FALSE.

      
      DO WHILE (particles(IP)%DTRIM .GT. 0.) ! Repeat the procedure until step is done

         DTCOLL = particles(IP)%DTRIM ! Looking for collisions within the remaining time
         ! ______ ADVECTION ______

         BOUNDCOLL = -1
         DO i = 1, 4 ! Check collisions with boundaries (xmin, xmax, ymin, ymax)
            ! Compute he velocity normal to the boundary
            VN = particles(IP)%VX * NX(i) + particles(IP)%VY * NY(i) + particles(IP)%VZ * NZ(i)
            ! Compute the distance from the boundary
            DX = (XW(i) - particles(IP)%X) * NX(i) + (YW(i) - particles(IP)%Y) * NY(i) + (ZW(i) - particles(IP)%Z) * NZ(i)
            ! Check if a collision happens (sooner than previously calculated)
            IF (VN * DTCOLL .LE. DX) THEN
               
               ! Find candidate collision point (don't need this for boundaries)
               ! XCOLL = particles(IP)%X + particles(IP)%VX * DTCOLL 
               ! YCOLL = particles(IP)%Y + particles(IP)%VY * DTCOLL
               ! ZCOLL = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               DTCOLL = DX/VN
               BOUNDCOLL = i                                 

            END IF
         END DO


         ! Check collisions with surfaces
         ! If earlier than dtcoll remember to set BOUNDCOLL to -1 and the new dtcoll


         IF (BOUNDCOLL == 1) THEN ! Collision with XMIN
            ! Tally boundary properties
            IF (BOOL_X_PERIODIC) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%X = particles(IP)%X + XMAX - XMIN
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

            ELSE IF (BOOL_XMIN_SPECULAR) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%VX = - particles(IP)%VX
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

            ELSE

               REMOVE_PART(IP) = .TRUE.
               particles(IP)%DTRIM = 0.

            END IF
         ELSE IF (BOUNDCOLL == 2) THEN ! Collision with XMAX
            ! Tally boundary properties
            IF (BOOL_X_PERIODIC) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%X = particles(IP)%X - XMAX + XMIN
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

            ELSE IF (BOOL_XMAX_SPECULAR) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%VX = - particles(IP)%VX
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

            ELSE

               REMOVE_PART(IP) = .TRUE.
               particles(IP)%DTRIM = 0.

            END IF
         ELSE IF (BOUNDCOLL == 3) THEN ! Collision with YMIN
            ! Tally boundary properties
            IF (BOOL_Y_PERIODIC) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%Y = particles(IP)%Y + YMAX - YMIN
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

            ELSE IF (BOOL_YMIN_SPECULAR) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%VY = - particles(IP)%VY
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

            ELSE

               REMOVE_PART(IP) = .TRUE.
               particles(IP)%DTRIM = 0.

            END IF
         ELSE IF (BOUNDCOLL == 4) THEN ! Collision with YMAX
            ! Tally boundary properties
            IF (BOOL_Y_PERIODIC) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%Y = particles(IP)%Y - YMAX + YMIN
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

            ELSE IF (BOOL_YMAX_SPECULAR) THEN

               CALL MOVE_PARTICLE(IP, DTCOLL)
               particles(IP)%VY = - particles(IP)%VY
               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL

            ELSE

               REMOVE_PART(IP) = .TRUE.
               particles(IP)%DTRIM = 0.

            END IF
         ELSE

            CALL MOVE_PARTICLE(IP, particles(IP)%DTRIM)
            particles(IP)%DTRIM = 0.

         END IF
         
            
         ! Advect particle in velocity - leapfrog method?
         ! particles(IP)%VX = ...
         ! particles(IP)%VY = ...
         ! particles(IP)%VZ = ...

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

   END SUBROUTINE ADVECT

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
