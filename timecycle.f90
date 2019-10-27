MODULE timecycle

   USE global
   USE mpi_common
   USE screen
   USE tools

   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE TIME_LOOP                                     !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE TIME_LOOP

   IMPLICIT NONE

   INTEGER :: tID 
   INTEGER :: NP_TOT
   REAL(KIND=8) :: CURRENT_TIME

   CHARACTER(len=512) :: stringTMP

   ! Init variables
   NP_TOT = 0

   DO tID = 1, NT

      ! ########### Print simulation info #######################################

      CURRENT_TIME = tID*DT

      CALL MPI_REDUCE(NP_PROC, NP_TOT, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      WRITE(stringTMP, '(A13,I5,A4,I8,A9,ES14.3,A28,I10)') '   Timestep: ', tID, ' of ', NT, &
                       ' - time: ', CURRENT_TIME, ' [s] - number of particles: ', NP_TOT

      CALL ONLYMASTERPRINT1(PROC_ID, TRIM(stringTMP))

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

      CALL DUMP_PARTICLES_FILE(tID)

      ! ~~~~~ Hmm that's it! ~~~~~

   END DO

   END SUBROUTINE TIME_LOOP


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE LINE_SOURCE_INJECT -> Injects particles from line source !!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE LINE_SOURCE_INJECT

   IMPLICIT NONE

   INTEGER      :: IP, IC
   REAL(KIND=8) :: DTFRAC
   REAL(KIND=8) :: X, Y, Z, VX, VY, VZ, EI, Vdummy
   TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW

   ! DDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDB
   REAL(KIND=8) :: RGAS = 200.0d0 ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBD
   INTEGER :: S_ID_DUMMY = -1 ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDDBDBDB
   ! DDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDBDBDDBDBDBDB

   IF (BOOL_LINESOURCE) THEN

      DO IP = 1,nfs_LINESOURCE 

         ! Create particle velocity, internal energy, position
         CALL MAXWELL(UX_LINESOURCE, UY_LINESOURCE, UZ_LINESOURCE,                        &
                      TTRA_LINESOURCE, TTRA_LINESOURCE, TTRA_LINESOURCE, TROT_LINESOURCE, &
                      Vdummy,VY,VZ,EI,RGAS)

         VX = UX_LINESOURCE + FLX(S_NORM_LINESOURCE,TTRA_LINESOURCE,RGAS)
 
         DTFRAC = rf()*DT
   
         X = X_LINESOURCE                     - VX*DTFRAC
         Y = Y_LINESOURCE + L_LINESOURCE*rf() - VY*DTFRAC
         Z = ZMIN         + (ZMAX-ZMIN)*rf()  - VZ*DTFRAC
         
         CALL CELL_FROM_POSITION(X,Y,  IC)

         ! Create particle
         CALL INIT_PARTICLE(X, Y, Z, VX, VY, VZ, EI, S_ID_DUMMY, IC, DTFRAC, particleNOW)

         ! Assign particle to particles array
         CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles)

      END DO

   END IF

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
         M = SPECIES(S_ID)%MOLMASS
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
         M = SPECIES(S_ID)%MOLMASS
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
         M = SPECIES(S_ID)%MOLMASS
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
         M = SPECIES(S_ID)%MOLMASS
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

   INTEGER      :: IP, IC
   REAL(KIND=8) :: DTRIM

   DO IP = 1, NP_PROC

      DTRIM = particles(IP)%DTRIM

      DO WHILE (DTRIM .GT. 0.) ! Repeat the procedure until step is done

         ! ______ ADVECTION ______

         ! Advect particle in space (forward Euler integration)
         particles(IP)%X = particles(IP)%X + DTRIM*particles(IP)%VX
         particles(IP)%Y = particles(IP)%Y + DTRIM*particles(IP)%VY
         particles(IP)%Z = particles(IP)%Z + DTRIM*particles(IP)%VZ
   
         ! Advect particle in velocity
         ! particles(IP)%VX = ...
         ! particles(IP)%VY = ...
         ! particles(IP)%VZ = ...

         ! Axisymmetric domain
         IF (BOOL_AXI .eqv. .TRUE.) THEN 
   
            ! Transform position, velocities etc
   
         END IF

         ! _______ APPLY PERIODIC BCs IF REQUESTED ________

         ! If the periodic BC is active, check if the particle exits one border or the other,
         ! and if so, bring it back in the domain.
       
         ! -------
         IF (BOOL_X_PERIODIC) THEN
      
            IF (particles(IP)%X .GT. XMAX) particles(IP)%X = XMIN + (particles(IP)%X - XMAX)
            IF (particles(IP)%X .LT. XMIN) particles(IP)%X = XMAX + (particles(IP)%X - XMIN)
      
         END IF
      
         ! -------
         IF (BOOL_Y_PERIODIC) THEN
      
            IF (particles(IP)%Y .GT. YMAX) particles(IP)%Y = YMIN + (particles(IP)%Y - YMAX)
            IF (particles(IP)%Y .LT. YMIN) particles(IP)%Y = YMAX + (particles(IP)%Y - YMIN)
      
         END IF 
      
         ! -------
         IF (BOOL_Z_PERIODIC) THEN
      
            IF (particles(IP)%Z .GT. ZMAX) particles(IP)%Z = ZMIN + (particles(IP)%Z - ZMAX)
            IF (particles(IP)%Z .LT. ZMIN) particles(IP)%Z = ZMAX + (particles(IP)%Z - ZMIN)
      
         END IF 
   
         ! _______ CHECK WHERE PARTICLE ENDED UP _______

         ! ++++++++ Check if particle is still in the domain ++++++++++++
         IF (particles(IP)%X .GE. XMIN .AND. particles(IP)%X .LE. XMAX .AND. & 
             particles(IP)%Y .GE. YMIN .AND. particles(IP)%Y .LE. YMAX .AND. &
             particles(IP)%Z .GE. ZMIN .AND. particles(IP)%Z .LE. ZMAX) THEN   ! Check that the particle is still in the domain
   
            ! Compute the index of the cell in which the particle ended up
            CALL CELL_FROM_POSITION(particles(IP)%X, particles(IP)%Y, IC)
            particles(IP)%IC = IC
   
            DTRIM = 0.E0 ! Timestep is over.
 
         ! +++++++++ Particle crossed domain boundaries ++++++++
         ELSE

            ! CALL IMPACT_BOUNDARY(particles(IP)%X, DTRIM) 
            DTRIM = 0.E0 ! TMP TMP TMP TMP TMP TMP TMP TMP TMP
         END IF

      END DO

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
      IF (particles(IP)%X .LT. XMIN .OR. particles(IP)%X .GT. XMAX .OR. & 
          particles(IP)%Y .LT. YMIN .OR. particles(IP)%Y .GT. YMAX .OR. &
          particles(IP)%Z .LT. ZMIN .OR. particles(IP)%Z .GT. ZMAX) THEN
   
         CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)

      END IF

      IP = IP - 1

   END DO

   END SUBROUTINE ADVECT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DSMC_COLLISIONS -> perform DSMC collisions !!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DSMC_COLLISIONS

   IMPLICIT NONE

   CALL ONLYMASTERPRINT1(PROC_ID, "DSMC COLLISIONS NOT IMPLEMENTED YET! IGNORING THIS!")

   END SUBROUTINE DSMC_COLLISIONS

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
