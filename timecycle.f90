MODULE timecycle

   USE global
   USE mpi_common
   USE screen
   USE tools
   USE collisions
   USE postprocess
   USE EM_fields

   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE TIME_LOOP                                     !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE TIME_LOOP

   IMPLICIT NONE
 
   INTEGER :: NP_TOT, NCOLL_TOT
   REAL(KIND=8) :: CURRENT_TIME

   CHARACTER(len=512) :: stringTMP

   ! Init variables
   NP_TOT = 0
   CALL INIT_POSTPROCESS
   
   ! Dump particles before the first time step, but after the initial seeding
   ! CALL DUMP_PARTICLES_FILE(0)

   CALL DUMP_GLOBAL_MOMENTS_FILE(0) ! Compute and dump moments at initial timestep

   tID = 1
   DO WHILE (tID .LE. NT)

      ! ########### Print simulation info #######################################

      CURRENT_TIME = tID*DT ! Current time is a variable of the "global.f90" module

      CALL MPI_REDUCE(NP_PROC, NP_TOT, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TIMESTEP_COLL, NCOLL_TOT, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      WRITE(stringTMP, '(A13,I9,A4,I8,A9,ES14.3,A28,I10,A25,I10)') '   Timestep: ', tID, ' of ', NT, &
                       ' - time: ', CURRENT_TIME, ' [s] - number of particles: ', NP_TOT, &
                       ' - number of collisions: ', NCOLL_TOT

      CALL ONLYMASTERPRINT1(PROC_ID, TRIM(stringTMP))

      ! ########### Inject particles from boundaries/injection sources ##########

      CALL BOUNDARIES_INJECT
      CALL LINE_SOURCE_INJECT

      ! ########### Advect particles ############################################

      ! CALL ADVECT 

      ! PRINT*, "           > ATTENTION! I'm using a weird advection! Check timecycle.f90 ~~~~~~~~~~~~~~~"
      ! CALL ADVECT_0D_ExB_ANALYTICAL(dt) ! Advects only velocity

      PRINT*, "           > ATTENTION! BCs NOT IMPLEMENTED! Check timecycle.f90"
      CALL ADVECT_BORIS

      ! ########### Exchange particles among processes ##########################

      CALL EXCHANGE

      ! ########### Perform collisions ##########################################

      IF (BOOL_MCC)  CALL MCC_COLLISIONS

      IF (BOOL_DSMC) CALL DSMC_COLLISIONS

      IF (BOOL_BGK)  CALL BGK_COLLISIONS

      IF (BOOL_CUSTOM_COLL) CALL CUSTOM_COLLISIONS

      ! ########### Dump quantities ##############################################

      ! ++++ Dump particles

      IF (DUMP_PART_EVERY .NE. -1) THEN ! Ok, dump particles
         IF (MOD(tID, DUMP_PART_EVERY) .EQ. 0) CALL DUMP_PARTICLES_FILE(tID)
      END IF

      ! ++++ Dump flowfield to VTK 

       IF (DUMP_GRID_START .NE. -1) THEN ! Ok, you can dump stuff
          IF (tID .GE. DUMP_GRID_START) THEN
             IF (MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY) .EQ. 0) CALL GRID_AVG_MOMENTS
             IF (MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY*DUMP_GRID_N_AVG) .EQ. 0) THEN
                CALL GRID_AVG_MOMENTS
                CALL GRID_SAVE_MOMENTS
                CALL GRID_RESET_MOMENTS
             END IF
          END IF
       END IF


!       IF (DUMP_GRID_START .NE. -1) THEN ! Ok, you can dump stuff
!          IF (tID .GE. DUMP_GRID_START) THEN
!             IF (MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY) .EQ. 0) CALL GRID_AVG
!             IF (MOD(tID-DUMP_GRID_START, DUMP_GRID_AVG_EVERY*DUMP_GRID_N_AVG) .EQ. 0) THEN
!                CALL GRID_AVG
!                CALL GRID_SAVE
!                CALL GRID_RESET
!             END IF
!          END IF
!       END IF

      ! ++++ Dump moments to file

      IF (DUMP_GLOB_MOM_EVERY .NE. -1) THEN ! Ok, you can dump moments
         IF (MOD(tID, DUMP_GLOB_MOM_EVERY) .EQ. 0) CALL DUMP_GLOBAL_MOMENTS_FILE(tID)
      END IF

      ! ~~~~~ Hmm that's it! ~~~~~

      tID = tID + 1 ! Don't forget this, at the end of the loop!
   END DO

   END SUBROUTINE TIME_LOOP
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE LINE_SOURCE_INJECT -> Injects particles from lines in the domain !!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE LINE_SOURCE_INJECT
  
      ! This subroutine injects particles from all defined linesources.
      ! The number of particles to be injected during a timestep DT was previously computed
      ! and stored in LINESOURCES(i)%nfs.
      ! Note that this is a floating point number! The actual number of particles to be 
      ! injected will be incremented stochastically, if a random number is lower than the
      ! decimal part of LINESOURCES(i)%nfs.

      IMPLICIT NONE
   
      INTEGER      :: IP, IC, IS, NFS_INT, ILINE
      REAL(KIND=8) :: DTFRAC, POS
      REAL(KIND=8) :: UX_LSRC, UY_LSRC, UZ_LSRC, TTRA_LSRC, TROT_LSRC ! Working variables
      REAL(KIND=8) :: X, Y, Z, VX, VY, VZ, EI 
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
   
      INTEGER :: S_ID
      REAL(KIND=8) :: M
   
      ! Linesource
      DO ILINE = 1, N_LINESOURCES
        
         UX_LSRC   = LINESOURCES(ILINE)%UX
         UY_LSRC   = LINESOURCES(ILINE)%UY
         UZ_LSRC   = LINESOURCES(ILINE)%UZ
         TTRA_LSRC = LINESOURCES(ILINE)%TTRA
         TROT_LSRC = LINESOURCES(ILINE)%TROT
         S_ID      = LINESOURCES(ILINE)%S_ID

         M = SPECIES(S_ID)%MOLMASS

         ! Number of particles to be injected is the lower integer of nfs, and we add or not one more particle
         ! taking the remainder.
         ! Example: if nfs = 44.7, we inject 44 particles, plus another one with probability 0.7.

         NFS_INT = FLOOR(LINESOURCES(ILINE)%nfs)
         IF (LINESOURCES(ILINE)%nfs-REAL(NFS_INT, KIND=8) .GE. rf()) THEN ! Same as SPARTA's perspeciess
            NFS_INT = NFS_INT + 1
         END IF
   
         ! Loop on particles and inject them
         DO IP = 1, NFS_INT 
   
            CALL MAXWELL(UX_LSRC, UY_LSRC, UZ_LSRC, &
                         TTRA_LSRC, TTRA_LSRC, TTRA_LSRC, TROT_LSRC, &
                         VX, VY, VZ, EI, M)

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

   INTEGER      :: IP, IC, i
   INTEGER      :: BOUNDCOLL
   REAL(KIND=8) :: DTCOLL
   REAL(KIND=8), DIMENSION(4) :: NX, NY, NZ, XW, YW, ZW
   REAL(KIND=8) :: XCOLL, YCOLL, ZCOLL
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

         DTCOLL = particles(IP)%DTRIM ! Start by looking for collisions within the remaining time

         ! ______ ADVECTION ______

         BOUNDCOLL = -1
         DO i = 1, 4 ! Check collisions with boundaries
            VN = particles(IP)%VX * NX(i) + particles(IP)%VY * NY(i) + particles(IP)%VZ * NZ(i) ! The velocity normal to the boundary
            DX = (XW(i) - particles(IP)%X) * NX(i) + (YW(i) - particles(IP)%Y) * NY(i) + (ZW(i) - particles(IP)%Z) * NZ(i) ! The distance from the boundary
            IF (VN * DTCOLL .LE. DX) THEN
               
               XCOLL = particles(IP)%X + particles(IP)%VX * DTCOLL ! Find candidate collision point
               YCOLL = particles(IP)%Y + particles(IP)%VY * DTCOLL
               ZCOLL = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               DTCOLL = DX/VN
               BOUNDCOLL = i                                 

            END IF
         END DO

         ! Check collisions with surfaces
         ! If earlier than dtcoll remember to set BOUNDCOLL to -1 and the new dtcoll


         IF (BOUNDCOLL == 1) THEN ! Collision with XMIN
            ! Tally boundary properties
            IF (BOOL_X_PERIODIC) THEN
               particles(IP)%X = particles(IP)%X + particles(IP)%VX * DTCOLL + XMAX - XMIN
               particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * DTCOLL
               particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
            ELSE IF (BOOL_XMIN_SPECULAR) THEN
               particles(IP)%X = particles(IP)%X + particles(IP)%VX * DTCOLL
               particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * DTCOLL
               particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               particles(IP)%VX = - particles(IP)%VX

               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
            ELSE
               REMOVE_PART(IP) = .TRUE.

               particles(IP)%DTRIM = 0.
            END IF
         ELSE IF (BOUNDCOLL == 2) THEN ! Collision with XMAX
            ! Tally boundary properties
            IF (BOOL_X_PERIODIC) THEN
               particles(IP)%X = particles(IP)%X + particles(IP)%VX * DTCOLL - XMAX + XMIN
               particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * DTCOLL
               particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
            ELSE IF (BOOL_XMAX_SPECULAR) THEN
               particles(IP)%X = particles(IP)%X + particles(IP)%VX * DTCOLL
               particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * DTCOLL
               particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               particles(IP)%VX = - particles(IP)%VX

               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
            ELSE
               REMOVE_PART(IP) = .TRUE.

               particles(IP)%DTRIM = 0.
            END IF
         ELSE IF (BOUNDCOLL == 3) THEN ! Collision with YMIN
            ! Tally boundary properties
            IF (BOOL_Y_PERIODIC) THEN
               particles(IP)%X = particles(IP)%X + particles(IP)%VX * DTCOLL
               particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * DTCOLL + YMAX - YMIN
               particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
            ELSE IF (BOOL_YMIN_SPECULAR) THEN
               particles(IP)%X = particles(IP)%X + particles(IP)%VX * DTCOLL
               particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * DTCOLL
               particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               particles(IP)%VY = - particles(IP)%VY

               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
            ELSE
               REMOVE_PART(IP) = .TRUE.

               particles(IP)%DTRIM = 0.
            END IF
         ELSE IF (BOUNDCOLL == 4) THEN ! Collision with YMIN
            ! Tally boundary properties
            IF (BOOL_Y_PERIODIC) THEN
               particles(IP)%X = particles(IP)%X + particles(IP)%VX * DTCOLL
               particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * DTCOLL - YMAX + YMIN
               particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
            ELSE IF (BOOL_YMAX_SPECULAR) THEN
               particles(IP)%X = particles(IP)%X + particles(IP)%VX * DTCOLL
               particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * DTCOLL
               particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * DTCOLL

               particles(IP)%VY = - particles(IP)%VY

               particles(IP)%DTRIM = particles(IP)%DTRIM - DTCOLL
            ELSE
               REMOVE_PART(IP) = .TRUE.

               particles(IP)%DTRIM = 0.
            END IF
         ELSE
            particles(IP)%X = particles(IP)%X + particles(IP)%VX * particles(IP)%DTRIM
            particles(IP)%Y = particles(IP)%Y + particles(IP)%VY * particles(IP)%DTRIM
            particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ * particles(IP)%DTRIM
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

   TIMESTEP_COLL = 0 ! Initialize number of collisions done during the timestep.

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

         ! Update number of collisions performed
         TIMESTEP_COLL = TIMESTEP_COLL + 1
      END IF

   END DO

   END SUBROUTINE MCC_COLLISIONS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ADVECT_0D_ExB_ANALYTICAL -> Advects particles in the domain     !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ADVECT_0D_ExB_ANALYTICAL(dt)

      ! This subroutine only changes the velocity (0D simulation), and uses an analytical law for
      ! the motion in the ExB field.
      ! E is assumed to be along the z direction and B is assumed to be along x.

      INTEGER :: IP
      REAL(KIND=8), INTENT(IN) :: dt
      REAL(KIND=8) :: E, B !Electric and magnetic fields
      REAL(KIND=8) :: VY0, VZ0, M, Q, omega_cycl

      E = 20000
      B = 0.01

      DO IP = 1, NP_PROC ! Loop on particles

         VY0 = particles(IP)%VY
         VZ0 = particles(IP)%VZ

         M = SPECIES(particles(IP)%S_ID)%MOLMASS
         Q = SPECIES(particles(IP)%S_ID)%CHARGE

         omega_cycl = Q*B/M

         ! This assumes electric field is along z and magnetic field along x
         particles(IP)%VY =  (VY0 - E/B)*cos(omega_cycl*dt) + VZ0*sin(omega_cycl*dt) + E/B;
         particles(IP)%VZ = -(VY0 - E/B)*sin(omega_cycl*dt) + VZ0*cos(omega_cycl*dt);

      END DO

   END SUBROUTINE ADVECT_0D_ExB_ANALYTICAL

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ADVECT_BORIS -> Advects particles in the domain using Boris scheme !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ADVECT_BORIS

      ! This subroutine advects particles using the Boris pusher.
      ! E is assumed to be along the x direction and B is assumed to be along z.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS FUNCTION IS TO BE VERIFIED!!!!!!!!!!!!!!!!!!!!!!!! 
! THIS FUNCTION IS TO BE VERIFIED!!!!!!!!!!!!!!!!!!!!!!!!
! THIS FUNCTION IS TO BE VERIFIED!!!!!!!!!!!!!!!!!!!!!!!!
! THIS FUNCTION IS TO BE VERIFIED!!!!!!!!!!!!!!!!!!!!!!!!
! THIS FUNCTION IS TO BE VERIFIED!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROPER TREATMENT FOR BOUNDARIES TO BE FORMULATED!!!!!!!
! PROPER TREATMENT FOR BOUNDARIES TO BE FORMULATED!!!!!!!
! PROPER TREATMENT FOR BOUNDARIES TO BE FORMULATED!!!!!!!
! PROPER TREATMENT FOR BOUNDARIES TO BE FORMULATED!!!!!!!
! PROPER TREATMENT FOR BOUNDARIES TO BE FORMULATED!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER :: IP, IC
      REAL(KIND=8) :: E, B !Electric and magnetic fields
      REAL(KIND=8) :: M, Q
      REAL(KIND=8), DIMENSION(3) :: E_vec, B_vec, V_now, V_minus, V_plus, V_prime, t_vec, s_vec
      REAL(KIND=8) :: DTRIM, Xold, Yold, Zold

      ! Loop on particles from the LAST ONE to the FIRST ONE. 
      ! In this way, if a particle exits the domain, you can just remove it, and this doesn't
      ! mix the position of the previous particles.
      DO IP = NP_PROC, 1, -1

         DTRIM = particles(IP)%DTRIM ! Remaining time

         ! Compute fields at the particle position
         CALL GET_E_B_POSITION_1D(particles(IP)%X, E, B)
         E_vec = (/ E,     0.0d0, 0.0d0 /)
         B_vec = (/ 0.0d0, 0.0d0, B /)

         ! Advect particle in space (using the current velocity, which BTW should be staggered)
         Xold = particles(IP)%X
         Yold = particles(IP)%Y
         Zold = particles(IP)%Z

         particles(IP)%X = particles(IP)%X + particles(IP)%VX*DTRIM 
         particles(IP)%Y = particles(IP)%Y + particles(IP)%VY*DTRIM
         particles(IP)%Z = particles(IP)%Z + particles(IP)%VZ*DTRIM

         ! Find new particle velocity
         M = SPECIES(particles(IP)%S_ID)%MOLMASS
         Q = SPECIES(particles(IP)%S_ID)%CHARGE

         V_now   = (/particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ/)
         V_minus = V_now + Q/M*E_vec*DTRIM/2
         t_vec   = Q/M*B_vec*DTRIM/2
         V_prime = V_minus + CROSS(V_minus, t_vec)
         s_vec   = 2*t_vec/(1 + DOT(t_vec, t_vec))
         V_plus  = V_minus + CROSS(V_prime, s_vec)

         particles(IP)%VX = V_plus(1) + Q/M*DTRIM/2*E_vec(1)
         particles(IP)%VY = V_plus(2) + Q/M*DTRIM/2*E_vec(2)
         particles(IP)%VZ = V_plus(3) + Q/M*DTRIM/2*E_vec(3)

         ! Apply boundary condition
          CALL IMPOSE_BOUNDARIES_TEST(IP, Xold, Yold, Zold)

         ! Convection is over, put DTRIM = DT for the next timestep
         particles(IP)%DTRIM = DT 

      END DO

      ! Finally, after particles have even collided with boundaries, find their new cell ID.
      ! Note that NP_PROC is automatically updated if you remove particles via "REMOVE_PARTICLE_..."

      DO IP = 1, NP_PROC
         CALL CELL_FROM_POSITION(particles(IP)%X, particles(IP)%Y, IC)
         particles(IP)%IC = IC
      END DO

   END SUBROUTINE ADVECT_BORIS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE IMPOSE_BOUNDARIES_TEST -> Imposes the boundaries !!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE IMPOSE_BOUNDARIES_TEST(IP, Xold, Yold, Zold)

      ! Checks the position of particle.
      ! IP: particle identifier in the "particles" vector
      ! Xold, Yold, Zold: position before advection took place (we need it for the BCs)
      ! 
      ! ATTENTION!!! The MOD() fortran operation does not work well with the left boundary!
      ! Therefore, I implemented it by hand. 
 
      INTEGER, INTENT(IN) :: IP
      REAL(KIND=8), INTENT(IN) :: Xold, Yold, Zold

      ! =========== Are periodic BCs requested? =========
      ! Then if the particle is outside the domain, put it back in.
      ! Do it recursively, since a very fast particle could be "very very" outside.
      ! Another way (probably more efficient): define how many domains it's "out", and 
      ! put it back by this amount.
      ! 
      ! BEWARE THE MOD() FUNCTION IN F90!! DOESN'T WORK WELL FOR THE LEFT INTERFACE!!

      IF (BOOL_X_PERIODIC) THEN ! Periodicity along X
 
         DO WHILE (particles(IP)%X < XMIN)
            particles(IP)%X = particles(IP)%X + (XMAX - XMIN)
         END DO

         DO WHILE (particles(IP)%X > XMAX)
            particles(IP)%X = particles(IP)%X - (XMAX - XMIN)
         END DO

      END IF
 
      IF (BOOL_Y_PERIODIC) THEN ! Periodicity along y

         DO WHILE (particles(IP)%Y < YMIN)
            particles(IP)%Y = particles(IP)%Y + (YMAX - YMIN)
         END DO

         DO WHILE (particles(IP)%Y > YMAX)
            particles(IP)%Y = particles(IP)%Y - (YMAX - YMIN)
         END DO

      END IF

      IF (BOOL_Z_PERIODIC) THEN ! Periodicity along z

         DO WHILE (particles(IP)%Z < ZMIN)
            particles(IP)%Z = particles(IP)%Z + (ZMAX - ZMIN)
         END DO

         DO WHILE (particles(IP)%Z > ZMAX)
            particles(IP)%Z = particles(IP)%Z - (ZMAX - ZMIN)
         END DO

      END IF

      ! ========= Remove particles outside of the domain =========
      ! At this point, if a particle is still out of the domain 
      ! it means that no condition was imposed -> remove it.

      IF ( (particles(IP)%X .LT. XMIN) .OR. (particles(IP)%X .GT. XMAX) .OR. &
           (particles(IP)%Y .LT. YMIN) .OR. (particles(IP)%Y .GT. YMAX) .OR. &
           (particles(IP)%Z .LT. ZMIN) .OR. (particles(IP)%Z .GT. ZMAX)) THEN

         CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)

      END IF

   END SUBROUTINE IMPOSE_BOUNDARIES_TEST
 

END MODULE timecycle
