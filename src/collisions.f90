! Contains subroutines related to particle collisions

MODULE collisions

   USE global
   USE mpi_common
   USE screen
   USE tools
   
   IMPLICIT NONE
   
   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DSMC_COLLISIONS -> Computes DSMC collisions on the grid  !
   ! Called by TIME_LOOP if DSMC collisions are active                   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DSMC_COLLISIONS 

      ! First, reorder particles
      ! IND => index of particle in particle array, indexed by particles I own and sorted by cells
      ! NPC => number of particles in cell, indexed by cell index (a block length for IND) 
      ! could this be indexed by all cells? for those that I don't own, NPC will turn out to be 0.
      ! + would simplify significantly,
      ! - could be up to about 50 times bigger
      ! + it's just a list of integers
      ! --->  I am gonna try and in case revert back.
      ! IOF => index (for IND) of first particle in a cell (an offset in IND)
      
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NPC, IOF
      INTEGER, DIMENSION(:), ALLOCATABLE :: NPCALL, IOFALL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IND, INDALL
      INTEGER                            :: JP, JS, JC, IDX, IDXALL
      INTEGER                            :: NCOLLREAL

      ALLOCATE(NPC(N_SPECIES,NCELLS))
      ALLOCATE(IOF(N_SPECIES,NCELLS))

      ALLOCATE(NPCALL(NCELLS))
      ALLOCATE(IOFALL(NCELLS))
      ! Count the number of particles in each cell of each species, to allocate arrays later
      NPC = 0
      NPCALL = 0
      DO JP = 1, NP_PROC
         JC = particles(JP)%IC
         JS = particles(JP)%S_ID
         NPC(JS,JC) = NPC(JS,JC) + 1
         NPCALL(JC) = NPCALL(JC) + 1
      END DO

      ! Fill the array of offsets (IOF). IOF(IS,IC) is the the index in IND of the first
      ! particle in cell IC of species IS
      IOF = -1
      IOFALL = -1
      IDX = 1
      IDXALL = 1
      ! DO JC = 1, NCELLS
      !    IF (NPCALL(JC) .NE. 0) THEN
      !       IOFALL(JC) = IDXALL
      !       IDXALL = IDXALL + NPCALL(JC)
      !    END IF
      !    DO JS = 1, N_SPECIES
      !       IF (NPC(JS,JC) .NE. 0) THEN
      !          IOF(JS,JC) = IDX
      !          IDX = IDX + NPC(JS,JC)
      !       END IF
      !    END DO
      ! END DO
      DO JC = 1, NCELLS
         IOFALL(JC) = IDXALL
         IDXALL = IDXALL + NPCALL(JC)
         DO JS = 1, N_SPECIES
            IOF(JS,JC) = IDX
            IDX = IDX + NPC(JS,JC)
         END DO
      END DO

      ! Allocate and fill the array of indices (IND). IND(IP) is the particle index (for the "particles" array)
      ! but ordered by cells, with offsets given by IND(IC) and stride length NPC(IC)
      ALLOCATE(IND(NP_PROC))
      ALLOCATE(INDALL(NP_PROC))
   
      NPC = 0
      NPCALL = 0
      DO JP = 1, NP_PROC
         JC = particles(JP)%IC
         JS = particles(JP)%S_ID
         IND(IOF(JS,JC) + NPC(JS,JC)) = JP
         NPC(JS,JC) = NPC(JS,JC) + 1

         INDALL(IOFALL(JC) + NPCALL(JC)) = JP
         NPCALL(JC) = NPCALL(JC) + 1
      END DO

      ! Verify: ! DBDBDBDBDBDDBDBDBDBDBDBDDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBD
      DO JC = 1, NCELLS
         DO JS = 1, N_SPECIES
            DO JP = IOF(JS,JC), IOF(JS,JC)+NPC(JS,JC)-1
               IF (particles(IND(JP))%S_ID .NE. JS) WRITE(*,*) 'Error in sorting!'
            END DO
         END DO
      END DO
   
      ! Compute collisions between particles
      TIMESTEP_COLL = 0
      TIMESTEP_REAC = 0
      DO JC = 1, NCELLS
         IF (SUM(NPC(:,JC)) .GT. 1) THEN
            ! For cells where there is at least two particles, call the collision procedure.
            IF (COLLISION_TYPE == DSMC) THEN
               ! DSMC temporarily broken because now arrays are per-species.
               CALL VSS_COLLIS(JC, NPCALL, IOFALL, INDALL, NCOLLREAL)
            ELSE IF (COLLISION_TYPE == DSMC_VAHEDI) THEN
               CALL VAHEDI_COLLIS(JC, NPC, IOF, IND, NCOLLREAL)
            END IF
            ! Add to the total number of collisions for this process
            TIMESTEP_COLL = TIMESTEP_COLL + NCOLLREAL
         END IF
      END DO
   
      !WRITE(*,*) 'Number of real collisions: ', TIMESTEP_COLL
      DEALLOCATE(NPC)
      DEALLOCATE(IOF)
      DEALLOCATE(IND)

      DEALLOCATE(NPCALL)
      DEALLOCATE(IOFALL)
      DEALLOCATE(INDALL)
   
   END SUBROUTINE DSMC_COLLISIONS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE VSS_COLLIS -> Compute collisions with VSS model !!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE VSS_COLLIS(JC,NPC,IOF,IND, NCOLLREAL)

      ! Computes the collisions using the VSS (or VHS, HS, depending on parameters)
      ! in cell JC. Needs the particles to be sorted by cell. This is done by the calling
      ! function DSMC_COLLISIONS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: JC
      INTEGER, DIMENSION(:), INTENT(IN) :: NPC, IOF
      INTEGER, DIMENSION(:), INTENT(IN) :: IND
      INTEGER, INTENT(OUT) :: NCOLLREAL


      INTEGER      :: IOFJ,IOLJ,NPCJ,JP1,JP2,JCOL, JP, INDJ, I, JR
      INTEGER      :: SP_ID1, SP_ID2, P1_SP_ID, P2_SP_ID, P3_SP_ID
      INTEGER      :: NCOLL,NCOLLMAX_INT
      REAL(KIND=8) :: SIGMA, COLLSIGMA, ALPHA, NCOLLMAX,FCORR,VR,VR2,rfp,SIGMA_R
      REAL(KIND=8) :: OMEGA, CREF, MRED, COLLPROB, PTCE
      REAL(KIND=8) :: TRDOF, PROT1, PROT2, PVIB1, PVIB2, EI, ETR, ECOLL, TOTDOF, EA, EROT, EVIB
      !REAL(KIND=8) :: B,C,EINT,ETOT,ETR,PHI,SITETA,VRX,VRY,VRZ
      REAL(KIND=8) :: VXMAX,VXMIN,VYMAX,VYMIN,VZMAX,VZMIN,VRMAX
      REAL(KIND=8) :: PI2
      REAL(KIND=8), DIMENSION(3) :: C1, C2, GREL, W
      REAL(KIND=8) :: GX, GY, GZ, G
      REAL(KIND=8) :: COSCHI, SINCHI, THETA, COSTHETA, SINTHETA
      REAL(KIND=8) :: M1, M2, COSA, SINA, BB
      REAL(KIND=8) :: CFNUM, VOL
      LOGICAL      :: SKIP

      TYPE(PARTICLE_DATA_STRUCTURE) :: NEWparticle

      PI2  = 2.*PI

      IOFJ  = IOF(JC)               ! First particle in the cell
      IOLJ  = IOF(JC) + NPC(JC) - 1 ! Last particle in the cell
      NPCJ = NPC(JC)                ! Number of particles in the cell

      IF (BOOL_RADIAL_WEIGHTING) THEN
         CFNUM = CELL_FNUM(JC)
      ELSE
         CFNUM = FNUM
      END IF
      ! Step 1. Compute the number of pairs to test for collision
      ! This in not very efficient for multispecies, could be improved by grouping species.
      ! For now consider all as one group (just one NCOLLMAX)

      VXMAX = -1.D+38
      VXMIN =  1.D+38
      VYMAX = -1.D+38
      VYMIN =  1.D+38
      VZMAX = -1.D+38
      VZMIN =  1.D+38 

      DO JP = IOFJ,IOLJ ! Find velocity envelope
         INDJ = IND(JP)
         VXMIN = DMIN1(VXMIN,particles(INDJ)%VX)
         VXMAX = DMAX1(VXMAX,particles(INDJ)%VX)
         VYMIN = DMIN1(VYMIN,particles(INDJ)%VY)
         VYMAX = DMAX1(VYMAX,particles(INDJ)%VY)
         VZMIN = DMIN1(VZMIN,particles(INDJ)%VZ)
         VZMAX = DMAX1(VZMAX,particles(INDJ)%VZ)
      END DO

      ! The maximum cross section of all the species involved in the collisions is set when the collision file is read.

      ! Compute the "worst case scenario" relative velocity
      VRMAX = SQRT((VXMAX-VXMIN)**2 + (VYMAX-VYMIN)**2 + (VZMAX-VZMIN)**2)
      ! Compute the maximum expected number of collisions
      IF (GRID_TYPE == RECTILINEAR_UNIFORM .AND. DIMS == 2) THEN
         VOL = CELL_VOL
      ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
         VOL = U2D_GRID%CELL_VOLUMES(JC)
      ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
         VOL = U3D_GRID%CELL_VOLUMES(JC)
      END IF

      NCOLLMAX = 0.5*NPC(JC)*(NPC(JC)-1)*SIGMAMAX*VRMAX*CFNUM*DT/VOL

      NCOLLMAX_INT = FLOOR(NCOLLMAX+0.5)

      
      ! Step 2. Compute the number of pairs (real+artificial) => add constraints (min/max) 

      IF (NCOLLMAX_INT .LT. 1) THEN

         NCOLL = 1

      ELSE IF (NCOLLMAX_INT .GT. FLOOR(0.5*NPC(JC))) THEN

         NCOLL = FLOOR(0.5*NPC(JC))

      ELSE 

         NCOLL = NCOLLMAX_INT

      END IF

      FCORR = NCOLLMAX/NCOLL
      !WRITE(*,*) 'Ncollmax_int', NCOLLMAX_INT, 'ncoll:', NCOLL, 'fcorr:', FCORR
      NCOLLREAL = 0

      ! Step 3. Perform the collision => actual probability correct via FCORR

      DO JCOL = 1, NCOLL
         SKIP = .FALSE.
         ! Select first collision partner randomly

         JP1 = IOFJ + INT(NPCJ*rf())
         JP1 = IND(JP1)

         ! Select second collision partner randomly (shouldn't be JP1)
         JP2 = JP1
         DO WHILE (JP2 == JP1)
            JP2 = IOFJ + INT(NPCJ*rf())
            JP2 = IND(JP2)
         END DO

         ! Compute the relative velocity
         VR2 = (particles(JP2)%VX - particles(JP1)%VX)**2 + &
               (particles(JP2)%VY - particles(JP1)%VY)**2 + &
               (particles(JP2)%VZ - particles(JP1)%VZ)**2 

         VR = SQRT(VR2)

         ! Get actual collision parameters for seleced pair
         SP_ID1 = particles(JP1)%S_ID
         SP_ID2 = particles(JP2)%S_ID

         !SIGMA = PI * (0.5 * (SPECIES(SP_ID1)%DIAM + SPECIES(SP_ID2)%DIAM))**2
         !OMEGA    = 0.5 * (SPECIES(SP_ID1)%OMEGA + SPECIES(SP_ID2)%OMEGA)
         !TREF = 0.5 * (SPECIES(SP_ID1)%TREF + SPECIES(SP_ID2)%TREF)
         CREF  = VSS_GREFS(SP_ID1, SP_ID2)
         SIGMA = VSS_SIGMAS(SP_ID1, SP_ID2)
         OMEGA = VSS_OMEGAS(SP_ID1, SP_ID2)
         ALPHA = VSS_ALPHAS(SP_ID1, SP_ID2)

         M1    = SPECIES(SP_ID1)%MOLECULAR_MASS
         M2    = SPECIES(SP_ID2)%MOLECULAR_MASS
         MRED  = M1*M2/(M1+M2)

         ETR = 0.5*MRED*VR2
         
         C1(1) = particles(JP1)%VX
         C1(2) = particles(JP1)%VY
         C1(3) = particles(JP1)%VZ

         C2(1) = particles(JP2)%VX
         C2(2) = particles(JP2)%VY
         C2(3) = particles(JP2)%VZ

         ! Check if collision happens
         rfp = rf()

         ! Compute (simulated) collision probability
         IF (ABS(OMEGA-0.5) .LT. 1.d-6) THEN
            COLLSIGMA = SIGMA
         ELSE
            IF (SIGMAMAX*VRMAX .LT. 1e-20) CALL ERROR_ABORT('The product is zero!')
            IF (VR .LT. 1e-20) CYCLE !CALL ERROR_ABORT('VR is zero!')
            COLLSIGMA = SIGMA*(VR/CREF)**(1.-2.*OMEGA)
         END IF
         COLLPROB = FCORR/(SIGMAMAX*VRMAX)*VR*COLLSIGMA


         IF ( rfp .LT. COLLPROB ) THEN

            ! Rimuovere commento per avere avviso
            IF (COLLPROB .GT. 1.) THEN
               WRITE(*,*) 'Attention => this was a bad DSMC collision!'
            END IF

            NCOLLREAL = NCOLLREAL + 1


            ! Test for chemical reaction with TCE model
            ECOLL = ETR + particles(JP1)%EROT + particles(JP2)%EROT + &
            particles(JP1)%EVIB + particles(JP2)%EVIB
            DO JR = 1, N_REACTIONS
               ! Check if species in collision belong to this reaction.
               ! If they don't, check the next reaction.
               ! Not very efficient with many reactions
               IF (REACTIONS(JR)%R1_SP_ID .NE. SP_ID1) THEN
                  IF (REACTIONS(JR)%R1_SP_ID .NE. SP_ID2) THEN
                     CYCLE
                  ELSE
                     IF (REACTIONS(JR)%R2_SP_ID .NE. SP_ID1) CYCLE
                  END IF
               ELSE
                  IF (REACTIONS(JR)%R2_SP_ID .NE. SP_ID2) CYCLE
               END IF

               EA = REACTIONS(JR)%EA

               IF (REACTIONS(JR)%TYPE == TCE) THEN
                  IF (ECOLL .LE. EA) CYCLE
                  PTCE = REACTIONS(JR)%C1 * (ECOLL-EA)**REACTIONS(JR)%C2 * (1.-EA/ECOLL)**REACTIONS(JR)%C3
               ELSE IF (REACTIONS(JR)%TYPE == LXCAT) THEN
                  SIGMA_R = INTERP_CS(ETR, REACTIONS(JR)%TABLE_ENERGY, REACTIONS(JR)%TABLE_CS)
                  PTCE = SIGMA_R / (SIGMA*(VR/CREF)**(1.-2.*OMEGA))
               ELSE
                  PTCE = 0
               END IF
               ! Here we suppose that the probability is so low that we can test sequentially with acceptable error
               rfp = rf()
               IF (rfp .LT. PTCE) THEN
                  SKIP = .TRUE.
                  TIMESTEP_REAC = TIMESTEP_REAC + 1
                  ! React
                  ECOLL = ECOLL - EA

                  ! Use R1 as P1 and assign internal energy to it.
                  ! Use R2 as P2 (or P2+P3) (Set R2 as the molecule that dissociates in P2+P3)
                  P1_SP_ID = REACTIONS(JR)%P1_SP_ID
                  P2_SP_ID = REACTIONS(JR)%P2_SP_ID

                  particles(JP1)%S_ID = REACTIONS(JR)%P1_SP_ID
                  particles(JP2)%S_ID = REACTIONS(JR)%P2_SP_ID

                  TOTDOF = 3. + SPECIES(P2_SP_ID)%ROTDOF + SPECIES(P2_SP_ID)%VIBDOF + &
                        SPECIES(P1_SP_ID)%ROTDOF
                  IF (REACTIONS(JR)%N_PROD == 3) THEN
                     P3_SP_ID = REACTIONS(JR)%P3_SP_ID
                     TOTDOF = TOTDOF + 3. + SPECIES(P3_SP_ID)%ROTDOF + SPECIES(P3_SP_ID)%VIBDOF
                  END IF

                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P1_SP_ID)%VIBDOF)
                  particles(JP1)%EVIB = EI
                  ECOLL = ECOLL - EI

                  TOTDOF = TOTDOF - SPECIES(P1_SP_ID)%ROTDOF
                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P1_SP_ID)%ROTDOF)
                  particles(JP1)%EROT = EI
                  ECOLL = ECOLL - EI

                  TOTDOF = TOTDOF - SPECIES(P2_SP_ID)%VIBDOF
                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P2_SP_ID)%VIBDOF)
                  particles(JP2)%EVIB = EI
                  ECOLL = ECOLL - EI

                  TOTDOF = TOTDOF - SPECIES(P2_SP_ID)%ROTDOF
                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P2_SP_ID)%ROTDOF)
                  particles(JP2)%EROT = EI
                  ECOLL = ECOLL - EI

                  M1 = SPECIES(P1_SP_ID)%MOLECULAR_MASS
                  IF (REACTIONS(JR)%N_PROD == 2) THEN
                     M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                  ELSE IF (REACTIONS(JR)%N_PROD == 3) THEN
                     M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS + SPECIES(P3_SP_ID)%MOLECULAR_MASS
                  END IF
                  TOTDOF = TOTDOF - 3.
                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, 3)
                  CALL HS_SCATTER(EI, M1, M2, C1, C2)
                  ECOLL = ECOLL - EI

                  
                  particles(JP1)%VX = C1(1)
                  particles(JP1)%VY = C1(2)
                  particles(JP1)%VZ = C1(3)

                  IF (REACTIONS(JR)%N_PROD == 2) THEN
                     particles(JP2)%VX = C2(1)
                     particles(JP2)%VY = C2(2)
                     particles(JP2)%VZ = C2(3)
                  ELSE IF (REACTIONS(JR)%N_PROD == 3) THEN
                     TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%VIBDOF
                     EVIB = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%VIBDOF)
                     ECOLL = ECOLL - EVIB
         
                     TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%ROTDOF
                     EROT = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%ROTDOF)
                     ECOLL = ECOLL - EROT

                     M1 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                     M2 = SPECIES(P3_SP_ID)%MOLECULAR_MASS
                     C1 = C2
                     
                     CALL HS_SCATTER(ECOLL, M1, M2, C1, C2)

                     particles(JP2)%VX = C1(1)
                     particles(JP2)%VY = C1(2)
                     particles(JP2)%VZ = C1(3)
                     

                     CALL INIT_PARTICLE(particles(JP2)%X,particles(JP2)%Y,particles(JP2)%Z, &
                     C2(1),C2(2),C2(3),EROT,EVIB,P3_SP_ID,JC,DT, NEWparticle)
                     !WRITE(*,*) 'Should be adding particle!'
                     CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)
                     
                  END IF

               EXIT
               END IF
            END DO

            IF (SKIP) CYCLE
            ! Perform elastic/inelastic collision
            ! Test for inelastic collision and TR/TV exchange
            TRDOF = 3.

            !TRDOF = 5. -2.*OMEGA why not?

            PROT1 = (TRDOF + SPECIES(SP_ID1)%ROTDOF)/TRDOF * SPECIES(SP_ID1)%ROTREL
            PROT2 = (TRDOF + SPECIES(SP_ID2)%ROTDOF)/TRDOF * SPECIES(SP_ID2)%ROTREL
            PVIB1 = (TRDOF + SPECIES(SP_ID1)%VIBDOF)/TRDOF * SPECIES(SP_ID1)%VIBREL
            PVIB2 = (TRDOF + SPECIES(SP_ID2)%VIBDOF)/TRDOF * SPECIES(SP_ID2)%VIBREL

            TRDOF = 5. -2.*OMEGA

            rfp = rf()
            IF (rfp .LT. PROT1) THEN
               ! Particle 1 selected for rotational energy exchange
               ECOLL = ETR + particles(JP1)%EROT
               particles(JP1)%EROT = COLL_INTERNAL_ENERGY(ECOLL, TRDOF, SPECIES(SP_ID1)%ROTDOF)
               ETR = ECOLL - particles(JP1)%EROT

            ELSE IF (rfp .LT. PROT1 + PROT2) THEN
               ! Particle 2 selected for rotational energy exchange
               ECOLL = ETR + particles(JP2)%EROT
               particles(JP2)%EROT = COLL_INTERNAL_ENERGY(ECOLL, TRDOF, SPECIES(SP_ID2)%ROTDOF)
               ETR = ECOLL - particles(JP2)%EROT

            ELSE IF (rfp .LT. PROT1 + PROT2 + PVIB1) THEN
               ! Particle 1 selected for vibrational energy exchange
               ECOLL = ETR + particles(JP1)%EVIB
               particles(JP1)%EVIB = COLL_INTERNAL_ENERGY(ECOLL, TRDOF, SPECIES(SP_ID1)%VIBDOF)
               ETR = ECOLL - particles(JP1)%EVIB

            ELSE IF (rfp .LT. PROT1 + PROT2 + PVIB1 + PVIB2) THEN
               ! Particle 2 selected for vibrational energy exchange
               ECOLL = ETR + particles(JP2)%EVIB
               particles(JP2)%EVIB = COLL_INTERNAL_ENERGY(ECOLL, TRDOF, SPECIES(SP_ID2)%VIBDOF)
               ETR = ECOLL - particles(JP2)%EVIB

            END IF

            ! Compute post-collision velocities
            ! VSS Collisions

            COSCHI = 2.*rf()**(1./ALPHA) - 1.
            SINCHI = SQRT(1.-COSCHI*COSCHI)
            THETA = PI2*rf()
            COSTHETA = COS(THETA)
            SINTHETA = SIN(THETA)


            IF (rfp .LT. PROT1 + PROT2 + PVIB1 + PVIB2) THEN
               ! Inelastic collision occurred => Randomize relative velocity vector
               G = SQRT(2.*ETR/MRED)

               COSA = 2.*rf()-1.0
               SINA = SQRT(1-COSA*COSA)
               BB = PI2 * rf()

               GX = G*SINA*COS(BB)
               GY = G*SINA*SIN(BB)
               GZ = G*COSA
            
            ELSE
               ! No inelastic collision occured => Take actual velocity vector
               GX = C1(1) - C2(1)
               GY = C1(2) - C2(2)
               GZ = C1(3) - C2(3)
               G = SQRT(GX*GX + GY*GY + GZ*GZ)

            END IF

            ! Compute post-collision velocities

            GREL(1) = GX*COSCHI + SQRT(GY*GY+GZ*GZ)*SINTHETA*SINCHI
            GREL(2) = GY*COSCHI + (G*GZ*COSTHETA - GX*GY*SINTHETA)/SQRT(GY*GY+GZ*GZ)*SINCHI
            GREL(3) = GZ*COSCHI - (G*GY*COSTHETA + GX*GZ*SINTHETA)/SQRT(GY*GY+GZ*GZ)*SINCHI

            ! Compute center of mass velocity vector
            DO I = 1, 3
               W(I) = M1/(M1+M2)*C1(I) + M2/(M1+M2)*C2(I)
            END DO

            ! Compute absolute velocity vector of the two particles
            DO I = 1, 3
               C1(I) = W(I) + M2/(M1+M2)*GREL(I)
               C2(I) = W(I) - M1/(M1+M2)*GREL(I)
            END DO

            particles(JP1)%VX = C1(1)
            particles(JP1)%VY = C1(2)
            particles(JP1)%VZ = C1(3)
               
            particles(JP2)%VX = C2(1)
            particles(JP2)%VY = C2(2)
            particles(JP2)%VZ = C2(3)
            

         END IF

      END DO  

      !WRITE(*,*) 'Actually performed:', NCOLLREAL
      !WRITE(*,*) NCOLL/(DT*NPC(JC))/MCRVHS, NCOLLREAL/(DT*NPC(JC))/MCRVHS 

      RETURN
         
   END SUBROUTINE VSS_COLLIS


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE VAHEDI_COLLIS -> Compute particle-particle collisions !!!
   ! with Vahedi's algorithm and tabluated cross-sections             !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE VAHEDI_COLLIS(JC,NPC,IOF,IND, NCOLLREAL)

      ! Computes the collisions using the VSS (or VHS, HS, depending on parameters)
      ! in cell JC. Needs the particles to be sorted by cell. This is done by the calling
      ! function DSMC_COLLISIONS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: JC
      INTEGER, DIMENSION(:,:), INTENT(IN) :: NPC, IOF
      INTEGER, DIMENSION(:), INTENT(IN) :: IND
      INTEGER, INTENT(OUT) :: NCOLLREAL


      INTEGER      :: JP1,JP2,JCOL, JP, INDJ, JR,IND1,IND2
      INTEGER      :: SP_ID1, SP_ID2, P1_SP_ID, P2_SP_ID, P3_SP_ID
      INTEGER      :: NCOLL,NCOLLMAX_INT
      REAL(KIND=8) :: NCOLLMAX,FCORR,VR,VR2,SIGMA_R,MAX_SIGMA
      REAL(KIND=8) :: MRED
      REAL(KIND=8) :: EI, ETR, ECOLL, TOTDOF, EA, EROT, EVIB
      !REAL(KIND=8) :: B,C,EINT,ETOT,ETR,PHI,SITETA,VRX,VRY,VRZ
      REAL(KIND=8) :: VXMAX,VXMIN,VYMAX,VYMIN,VZMAX,VZMIN,VRMAX
      REAL(KIND=8) :: PI2
      REAL(KIND=8), DIMENSION(3) :: C1, C2
      REAL(KIND=8) :: M1, M2
      REAL(KIND=8) :: CFNUM, VOL
      REAL(KIND=8) :: P_REACT, FACTOR
      LOGICAL, DIMENSION(:), ALLOCATABLE :: HAS_REACTED

      TYPE(PARTICLE_DATA_STRUCTURE) :: NEWparticle

      PI2  = 2.*PI

      IF (BOOL_RADIAL_WEIGHTING) THEN
         CFNUM = CELL_FNUM(JC)
      ELSE
         CFNUM = FNUM
      END IF

      ALLOCATE(HAS_REACTED(SUM(NPC(:,JC))))
      HAS_REACTED = .FALSE.

      DO JR = 1, N_REACTIONS

         REACTIONS(JR)%COUNTS = 0

         SP_ID1 = REACTIONS(JR)%R1_SP_ID
         SP_ID2 = REACTIONS(JR)%R2_SP_ID

         ! Test if there are enough particles of reactants
         IF ( (NPC(SP_ID1,JC) .LT. 1) .OR. (NPC(SP_ID2,JC) .LT. 1) ) EXIT

         ! Step 1. Compute the number of pairs to test for collision

         VXMAX = -1.D+38
         VXMIN =  1.D+38
         VYMAX = -1.D+38
         VYMIN =  1.D+38
         VZMAX = -1.D+38
         VZMIN =  1.D+38 

         DO JP = IOF(SP_ID1,JC), IOF(SP_ID1,JC)+NPC(SP_ID1,JC)-1 ! Find velocity envelope
            INDJ = IND(JP)
            VXMIN = DMIN1(VXMIN,particles(INDJ)%VX)
            VXMAX = DMAX1(VXMAX,particles(INDJ)%VX)
            VYMIN = DMIN1(VYMIN,particles(INDJ)%VY)
            VYMAX = DMAX1(VYMAX,particles(INDJ)%VY)
            VZMIN = DMIN1(VZMIN,particles(INDJ)%VZ)
            VZMAX = DMAX1(VZMAX,particles(INDJ)%VZ)
         END DO

         DO JP = IOF(SP_ID2,JC), IOF(SP_ID2,JC)+NPC(SP_ID2,JC)-1 ! Find velocity envelope
            INDJ = IND(JP)
            VXMIN = DMIN1(VXMIN,particles(INDJ)%VX)
            VXMAX = DMAX1(VXMAX,particles(INDJ)%VX)
            VYMIN = DMIN1(VYMIN,particles(INDJ)%VY)
            VYMAX = DMAX1(VYMAX,particles(INDJ)%VY)
            VZMIN = DMIN1(VZMIN,particles(INDJ)%VZ)
            VZMAX = DMAX1(VZMAX,particles(INDJ)%VZ)
         END DO


         ! Find the maximum cross section of all the species involved in the collisions
         !NSP = MIXTURES(DSMC_COLL_MIX)%N_COMPONENTS
         !SIGMAMAX = 0
         !DO JSP = 1, NSP
         !  SIGMA = SPECIES( MIXTURES(DSMC_COLL_MIX)%COMPONENTS(JSP)%ID )%SIGMA
         !  IF (SIGMA.GT.SIGMAMAX) SIGMAMAX = SIGMA
         !END DO

         ! Compute the "worst case scenario" relative velocity
         VRMAX = SQRT((VXMAX-VXMIN)**2 + (VYMAX-VYMIN)**2 + (VZMAX-VZMIN)**2)

         MAX_SIGMA = REACTIONS(JR)%MAX_SIGMA

         ! Compute the maximum expected number of collisions
         FACTOR = 1
         IF (SP_ID1 == SP_ID2) FACTOR = 0.5

         IF (GRID_TYPE == RECTILINEAR_UNIFORM .AND. DIMS == 2) THEN
            VOL = CELL_VOL
         ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
            VOL = U2D_GRID%CELL_VOLUMES(JC)
         ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
            VOL = U3D_GRID%CELL_VOLUMES(JC)
         END IF

         NCOLLMAX = FACTOR*NPC(SP_ID1,JC)*NPC(SP_ID2,JC)*MAX_SIGMA*VRMAX*CFNUM*DT/VOL

         NCOLLMAX_INT = FLOOR(NCOLLMAX+0.5)


         ! Step 2. Compute the number of pairs (real+artificial) => add constraints (min/max) 

         IF (NCOLLMAX_INT .LT. 1) THEN
            NCOLL = 1
         ELSE IF (NCOLLMAX_INT .GT. FLOOR(FACTOR*NPC(SP_ID1,JC))) THEN
            NCOLL = FLOOR(0.5*NPC(SP_ID1,JC))
         ELSE IF (NCOLLMAX_INT .GT. FLOOR(FACTOR*NPC(SP_ID2,JC))) THEN
            NCOLL = FLOOR(0.5*NPC(SP_ID2,JC))
         ELSE 
            NCOLL = NCOLLMAX_INT
         END IF

         FCORR = NCOLLMAX/NCOLL
         !WRITE(*,*) 'Ncollmax_int', NCOLLMAX_INT, 'ncoll:', NCOLL, 'fcorr:', FCORR
         NCOLLREAL = 0

         ! Step 3. Perform the collision => actual probability correct via FCORR
         WRITE(*,*) 'Testing ', NCOLL, ' collision pairs for reaction ', JR, ' with correction ', FCORR

         DO JCOL = 1, NCOLL

            ! Select a particle pair randomly.

            ! DO
            !    IND1 = IOF(SP_ID1,JC) + INT(NPC(SP_ID1,JC)*rf())
            !    JP1 = IND(IND1)
            !    IND1 = IND1 - IOF(1,JC) + 1
            ! END DO

            !!!!!!!! TODO: check if the IF condition for CYCLE Is correct
            IND1 = IOF(SP_ID1,JC) - IOF(1,JC) + INT(NPC(SP_ID1,JC)*rf()) + 1
            JP1 = IND(IND1 + IOF(1,JC) - 1)
            IF (HAS_REACTED(IND1)) CYCLE

            ! Select second collision partner randomly (shouldn't be JP1)
            ! DO
            !    IND2 = IOF(SP_ID2,JC) + INT(NPC(SP_ID2,JC)*rf())
            !    JP2 = IND(IND2)
            !    IND2 = IND2 - IOF(SP_ID2,JC) + 1
            !    IF ((.NOT. HAS_REACTED(IND2)) .AND. (JP2 .NE. JP1)) EXIT
            ! END DO
            IND2 = IOF(SP_ID2,JC) - IOF(1,JC) + INT(NPC(SP_ID2,JC)*rf()) + 1
            JP2 = IND(IND2 + IOF(1,JC) - 1)
            IF (HAS_REACTED(IND2)) CYCLE
            !!!!!!!!!!!!!!!!!!!!

            ! Compute the relative velocity
            VR2 = (particles(JP2)%VX - particles(JP1)%VX)**2 + &
                  (particles(JP2)%VY - particles(JP1)%VY)**2 + &
                  (particles(JP2)%VZ - particles(JP1)%VZ)**2 

            VR = SQRT(VR2)

            ! Get actual collision parameters for selected pair
            IF ((SP_ID1 .NE. particles(JP1)%S_ID) .OR. (SP_ID2 .NE. particles(JP2)%S_ID)) WRITE(*,*) 'Error in particle selection.'

            M1    = SPECIES(SP_ID1)%MOLECULAR_MASS
            M2    = SPECIES(SP_ID2)%MOLECULAR_MASS
            MRED  = M1*M2/(M1+M2)

            ETR = 0.5*MRED*VR2
            
            C1(1) = particles(JP1)%VX
            C1(2) = particles(JP1)%VY
            C1(3) = particles(JP1)%VZ

            C2(1) = particles(JP2)%VX
            C2(2) = particles(JP2)%VY
            C2(3) = particles(JP2)%VZ

            EA = REACTIONS(JR)%EA
            IF (ETR .LE. EA) CYCLE

            IF (REACTIONS(JR)%TYPE == LXCAT) THEN
               SIGMA_R = INTERP_CS(ETR, REACTIONS(JR)%TABLE_ENERGY, REACTIONS(JR)%TABLE_CS)
               P_REACT = FCORR/(MAX_SIGMA*VRMAX)*VR*SIGMA_R
            ELSE
               CYCLE
            END IF

            ! Try the reaction
            IF (rf() < P_REACT) THEN ! Collision happens

               TIMESTEP_COLL = TIMESTEP_COLL + 1
               REACTIONS(JR)%COUNTS = REACTIONS(JR)%COUNTS + 1
               HAS_REACTED(IND1) = .TRUE.
               HAS_REACTED(IND2) = .TRUE.

               ! Rimuovere commento per avere avviso
               IF (P_REACT .GT. 1.) THEN
                  !WRITE(*,*) 'Attention => this was a bad MCC collision! P_REACT = ', P_REACT
                  TIMESTEP_REAC = TIMESTEP_REAC + 1
               END IF


               !WRITE(*,*) 'Reacting!'
               ! React
               ECOLL = ETR - EA

               ! Use R1 as P1 and assign internal energy to it.
               ! Use R2 as P2 (or P2+P3) (Set R2 as the molecule that dissociates in P2+P3)
               P1_SP_ID = REACTIONS(JR)%P1_SP_ID
               P2_SP_ID = REACTIONS(JR)%P2_SP_ID

               !WRITE(*,*) 'SP_ID changed from ', particles(JP1)%S_ID, particles(JP2)%S_ID, &
               !' to ', REACTIONS(JR)%P1_SP_ID, REACTIONS(JR)%P2_SP_ID

               particles(JP1)%S_ID = REACTIONS(JR)%P1_SP_ID
               particles(JP2)%S_ID = REACTIONS(JR)%P2_SP_ID


               IF (.NOT. REACTIONS(JR)%IS_CEX) THEN
                  
                  TOTDOF = 3. + SPECIES(P2_SP_ID)%ROTDOF + SPECIES(P2_SP_ID)%VIBDOF + &
                         SPECIES(P1_SP_ID)%ROTDOF
                  IF (REACTIONS(JR)%N_PROD == 3) THEN
                     P3_SP_ID = REACTIONS(JR)%P3_SP_ID
                     TOTDOF = TOTDOF + 3. + SPECIES(P3_SP_ID)%ROTDOF + SPECIES(P3_SP_ID)%VIBDOF
                  END IF

                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P1_SP_ID)%VIBDOF)
                  particles(JP1)%EVIB = EI
                  ECOLL = ECOLL - EI

                  TOTDOF = TOTDOF - SPECIES(P1_SP_ID)%ROTDOF
                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P1_SP_ID)%ROTDOF)
                  particles(JP1)%EROT = EI
                  ECOLL = ECOLL - EI

                  TOTDOF = TOTDOF - SPECIES(P2_SP_ID)%VIBDOF
                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P2_SP_ID)%VIBDOF)
                  particles(JP2)%EVIB = EI
                  ECOLL = ECOLL - EI

                  TOTDOF = TOTDOF - SPECIES(P2_SP_ID)%ROTDOF
                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P2_SP_ID)%ROTDOF)
                  particles(JP2)%EROT = EI
                  ECOLL = ECOLL - EI

                  M1 = SPECIES(P1_SP_ID)%MOLECULAR_MASS
                  IF (REACTIONS(JR)%N_PROD == 2) THEN
                     M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                  ELSE IF (REACTIONS(JR)%N_PROD == 3) THEN
                     M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS + SPECIES(P3_SP_ID)%MOLECULAR_MASS
                  END IF
                  TOTDOF = TOTDOF - 3.
                  EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, 3)
                  !EI = ECOLL
                  !IF (TIMESTEP_COLL < 10) WRITE(*,*) 'Colliding particles ', JP1, JP2, ' with mass ', M1, ' and ', M2
                  !IF (TIMESTEP_COLL < 10) WRITE(*,*) 'Pre collision velocities  ', C1, ' and ', C2

                  CALL HS_SCATTER(EI, M1, M2, C1, C2)
                  !IF (TIMESTEP_COLL < 10) WRITE(*,*) 'Post collision velocities ', C1, ' and ', C2
                  ECOLL = ECOLL - EI

                  particles(JP1)%VX = C1(1)
                  particles(JP1)%VY = C1(2)
                  particles(JP1)%VZ = C1(3)


                  IF (REACTIONS(JR)%N_PROD == 2) THEN
                     particles(JP2)%VX = C2(1)
                     particles(JP2)%VY = C2(2)
                     particles(JP2)%VZ = C2(3)
                  ELSE IF (REACTIONS(JR)%N_PROD == 3) THEN
                     TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%VIBDOF
                     EVIB = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%VIBDOF)
                     ECOLL = ECOLL - EVIB
         
                     TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%ROTDOF
                     EROT = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%ROTDOF)
                     ECOLL = ECOLL - EROT
      
                     M1 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                     M2 = SPECIES(P3_SP_ID)%MOLECULAR_MASS
                     C1 = C2
                     
                     CALL HS_SCATTER(ECOLL, M1, M2, C1, C2)
                           
                     particles(JP2)%VX = C1(1)
                     particles(JP2)%VY = C1(2)
                     particles(JP2)%VZ = C1(3)


                     CALL INIT_PARTICLE(particles(JP2)%X,particles(JP2)%Y,particles(JP2)%Z, &
                     C2(1),C2(2),C2(3),EROT,EVIB,P3_SP_ID,particles(JP2)%IC,DT, NEWparticle)
                     !WRITE(*,*) 'Should be adding particle!'
                     CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)


                  END IF
               END IF

            END IF

         END DO

      END DO

      !WRITE(*,*) 'Actually performed:', NCOLLREAL
      !WRITE(*,*) NCOLL/(DT*NPC(JC))/MCRVHS, NCOLLREAL/(DT*NPC(JC))/MCRVHS 

      DEALLOCATE(HAS_REACTED)
         
   END SUBROUTINE VAHEDI_COLLIS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION FUNCTION_I -> Function used for energy redistribution     !
   ! between translational and internal degrees of freedom.             !
   ! Eq. (6.108) of Boyd and Schwartzentruber                           !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   REAL(KIND=8) FUNCTION FUNCTION_I(ECOLL, EI, TRDOF, IDOF)
      
      IMPLICIT NONE

      INTEGER, INTENT(IN)       :: IDOF
      REAL(KIND=8), INTENT(IN)  :: EI, ECOLL, TRDOF
      REAL(KIND=8)              :: RESULT

      IF (IDOF .EQ. 2) THEN
         RESULT = (1.-EI/ECOLL)**(TRDOF/2.-1.)
      ELSE
         RESULT = ((TRDOF+IDOF-4.)/(IDOF-2.))**(IDOF/2.-1.) * ((TRDOF+IDOF-4.)/(TRDOF-2.))**(TRDOF/2.-1.) * &
                     (EI/ECOLL)**(IDOF/2.-1.) * (1.-EI/ECOLL)**(TRDOF/2.-1.)
      END IF
      FUNCTION_I = RESULT

   END FUNCTION FUNCTION_I


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION COLL_INTERNAL_ENERGY -> Borgnakke-Larsen energy sampling  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   REAL(KIND=8) FUNCTION COLL_INTERNAL_ENERGY(ECOLL, TRDOF, INTDOF)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: INTDOF
      REAL(KIND=8), INTENT(IN)  :: ECOLL, TRDOF
      REAL(KIND=8) :: EI, R1, R2

      IF (TRDOF .LT. 1.d-6) THEN
         EI = ECOLL
      ELSE IF (INTDOF == 0) THEN
         EI = 0.d0
      ELSE
         DO
            R1 = rf()
            EI = R1*ECOLL
            R2 = rf()
            IF (R2 .LE. FUNCTION_I(ECOLL, EI, TRDOF, INTDOF)) EXIT
         END DO
      END IF
      COLL_INTERNAL_ENERGY = EI

   END FUNCTION COLL_INTERNAL_ENERGY


   REAL(KIND=8) FUNCTION COLL_INTERNAL_ENERGY_EQUAL(ECOLL, TRDOF, INTDOF)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: INTDOF
      REAL(KIND=8), INTENT(IN)  :: ECOLL, TRDOF

      COLL_INTERNAL_ENERGY_EQUAL = 0.5*ECOLL

   END FUNCTION COLL_INTERNAL_ENERGY_EQUAL

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE HS_SCATTER -> Computes velocities after a hard-sphere   !
   ! i.e. isotropic collision                                           !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE HS_SCATTER(ECOLL, M1, M2, C1, C2)
      
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: ECOLL, M1, M2
      REAL(KIND=8), INTENT(INOUT), DIMENSION(3) :: C1, C2
      REAL(KIND=8), DIMENSION(3) :: W, GREL
      REAL(KIND=8) :: PI2
      REAL(KIND=8) :: MRED, G, COSPHI, SINPHI, THETA
      INTEGER :: I

      PI2  = 2.*PI

      MRED  = M1*M2/(M1+M2)
      ! Relative velocity from given translational energy
      G = SQRT(2.*ECOLL/MRED)

      COSPHI = 2.*rf()-1.
      SINPHI = SQRT(1-COSPHI*COSPHI)
      THETA = PI2 * rf()

      GREL(1) = G*SINPHI*COS(THETA)
      GREL(2) = G*SINPHI*SIN(THETA)
      GREL(3) = G*COSPHI

      ! Compute center of mass velocity vector
      DO I = 1, 3
         W(I) = M1/(M1+M2)*C1(I) + M2/(M1+M2)*C2(I)
      END DO

      ! Compute absolute velocity vector of the two particles
      DO I = 1, 3
         C1(I) = W(I) + M2/(M1+M2)*GREL(I)
         C2(I) = W(I) - M1/(M1+M2)*GREL(I)
      END DO

   END SUBROUTINE HS_SCATTER

     
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! SUBROUTINE BGK_COLLISIONS -> Standalone procedure for       !
   ! simple BGK collisions                                       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   SUBROUTINE BGK_COLLISIONS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Implements BGK model for the collision operator.           !!
   ! Allows two types of relaxation, following the internal     !!
   ! variable "BGK_MODEL_TYPE_INT" (see INITCOLLISIONS          !!
   ! subroutine):                                               !!
   !                                                            !!
   ! BGK_MODEL_TYPE == 0                                        !!
   !   Classical BGK, performs relaxation towards a local       !!
   !   Maxwellian.                                              !!
   !                                                            !!
   ! BGK_MODEL_TYPE == 1                                        !!
   !   Light particles (such as electrons) colliding with       !!
   !   background medium.                                       !!
   !                                                            !!
   ! ATTENTION! WORKS ONLY FOR SINGLE SPECIES GAS MIXTURE!!!    !!
   ! The temperature of the local Maxwellian is found from the  !!
   ! internal energy (translational PLUS rotational, saved as   !!
   ! "EI" in the particle object).                              !!
   ! The local Maxwellian has a rotational temperature equal to !!
   ! the translational one.                                     !!
   !                                                            !!
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   ! First, create vectors holding the position of particles in the "particles" array,
   ! so that it's practical to extract them.
   ! Particles in the "particles" array are not actually moved aroud, in order to save computational
   ! time. Instead, the particle array is scanned so that:
   ! 1) the number of particles in each cell is computed
   ! 2) the ID of each particle is copied in the array IND, ordered by cell: first the particles
   !    belonging to the first cell, then these belonging to the second one etc.
   ! The IDs in the IND array allow to pick the particle quickly from the "particles" array.
   !
   !
   ! NPC => number of particles in each cell. Array of size NX*NY (total cells of the simulation)
   !        In the cells belonging to the current processor, there will be some particle, in other
   !        cells it will be zero.
   !
   ! IND => array of length NP_PROC, containing the ID of particles, ordered per cell.
   ! 
   ! IOF => index of first particle in a cell (it indicates the offset in IND). 
   !        Conceptually: IOF = cumsum(NPC).

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPC, IOF
      INTEGER, DIMENSION(:), ALLOCATABLE :: IND
      INTEGER                            :: JP, JC, IDX
   
      ! For actual collisions
      INTEGER :: IDp, N_DOF_TOT, JP_START, JP_END
      REAL(KIND=8) :: VX_CELL, VY_CELL, VZ_CELL, P_CELL, T_CELL, Ttr_CELL, Ttot_CELL, n_CELL, MASS, VOL
      REAL(KIND=8) :: nu, V_TH
      REAL(KIND=8) :: VX_NOW, VY_NOW, VZ_NOW, V2_NOW, EI_NOW
      REAL(KIND=8) :: Etot_cell, Etot_tra_cell, Ekin_cell, E_SUM, Etr_SUM
      REAL(KIND=8) :: Q, CHI, COS_TH, SIN_TH, COS_CHI, SIN_CHI

      ! Check that the gas is composed by only one species
      IF (SIZE(SPECIES) .NE. 1) THEN ! Print an error 
         PRINT*
         PRINT*, "  ATTENTION! BGK collisions work for single species gas only. Sorry about that."
         PRINT*, "  You may have one only species in your actual particle vector, but I see the SPECIES "
         PRINT*, "  structure has more entries, so I will stop here just to be sure. "
         PRINT*, "  Use only one species, or generalize this BGK function."
         PRINT*
         PRINT*, "  ABORTING!"
         PRINT*
         STOP
      END IF
   
      ! =========== Here, create vectors of particle indices ===========


      ALLOCATE(NPC(NCELLS))
      ALLOCATE(IOF(NCELLS))



      NPC = 0
      DO JP = 1, NP_PROC
         JC = particles(JP)%IC
         NPC(JC) = NPC(JC) + 1
      END DO
   

   
      IOF = -1
      IDX = 1
      DO JC = 1, NCELLS
         IF (NPC(JC) .NE. 0) THEN
            IOF(JC) = IDX
            IDX = IDX + NPC(JC)
         END IF
      END DO
   
      ALLOCATE(IND(NP_PROC))
      
      NPC = 0
      DO JP = 1, NP_PROC
         JC = particles(JP)%IC
         IND(IOF(JC) + NPC(JC)) = JP
         NPC(JC) = NPC(JC) + 1
      END DO
      
      ! =========== LOOP ON CELLS AND COMPUTE BGK COLLISIONS =========
   
      TIMESTEP_COLL = 0 ! Init number of collisions that happened
   
      DO JC = 1, NCELLS
   
         ! Note that in the current formulation, some cells may be not owned by me. 
         ! These will have zero particles inside, so I will skip them.
         IF (NPC(JC) == 0) CYCLE

         JP_START = IOF(JC)                ! First particle in the cell 
         JP_END   = IOF(JC) + NPC(JC) - 1  ! Last particle in the cell
   
         ! ++++++++ Compute average quantities in cell +++++++++++
         ! First, put internal quantities to zero
         VX_CELL = 0
         VY_CELL = 0
         VZ_CELL = 0
   
         E_SUM   = 0
         Etr_SUM = 0
   
         V_TH   = 0
         P_CELL = 0
   
         MASS = SPECIES(1)%MOLECULAR_MASS  ! ONLY 1 SPECIES FOR BGK!
   
         ! Compute average velocities and total energy
         DO JP = JP_START, JP_END ! Loop on particles in cell
            
            IDp = IND(JP) ! ID of current particle
   
            ! Extract velocities and internal energy
            VX_NOW   = particles(IDp)%VX
            VY_NOW   = particles(IDp)%VY
            VZ_NOW   = particles(IDp)%VZ
            EI_NOW   = particles(IDp)%EROT + particles(IDp)%EVIB
   
            ! Average velocity in the cell
            VX_CELL = VX_CELL + VX_NOW/NPC(JC)
            VY_CELL = VY_CELL + VY_NOW/NPC(JC)
            VZ_CELL = VZ_CELL + VZ_NOW/NPC(JC)
   
            ! Total energy, sum over the particles
            E_SUM   = E_SUM   + MASS*(VX_NOW**2 + VY_NOW**2 + VZ_NOW**2)/2 + MASS*EI_NOW ! [J] Total energy
            Etr_SUM = Etr_SUM + MASS*(VX_NOW**2 + VY_NOW**2 + VZ_NOW**2)/2               ! [J] Translational energy
   
         END DO
   
         IF (GRID_TYPE == RECTILINEAR_UNIFORM .AND. DIMS == 2) THEN
            VOL = CELL_VOL
         ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
            VOL = U2D_GRID%CELL_VOLUMES(JC)
         ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
            VOL = U3D_GRID%CELL_VOLUMES(JC)
         END IF


         n_CELL   = FNUM*NPC(JC)/VOL ! Number density in the cell
   
         Etot_cell     = n_CELL*E_SUM/NPC(JC)                                 ! [J/m3] total energy
         Etot_tra_cell = n_CELL*Etr_SUM/NPC(JC)                               ! [J/m3] total translational energy
         Ekin_cell     = MASS*n_CELL*(VX_CELL**2 + VY_CELL**2 + VZ_CELL**2)/2 ! [J/m3] (ordered) kinetic energy in cell
      
         N_DOF_TOT = 3 + SPECIES(1)%ROTDOF + SPECIES(1)%VIBDOF       ! Total number of DOFs. BGK works only for 1 species!!!
         T_CELL    = 2/(N_DOF_TOT*kB*n_CELL)*(Etot_cell - Ekin_cell) ! [K] Compute equilibr temp from internal energy
         Ttr_CELL  = 2/(3*kB*n_CELL)*(Etot_tra_cell - Ekin_cell)     ! [K] Compute translational temperature

         V_TH = sqrt(8*kB*MAX(Ttr_CELL, 0.0d0)/(MASS*pi)) ! [m/s] Thermal speed (from translational energy only!)
   
         Ttot_CELL = Ttr_CELL+MASS*(VX_CELL**2+VY_CELL**2+VZ_CELL**2)/(3.0d0*kB) ! [K] Temperature from total transl ener.

         ! +++++++ Compute collision frequency according to the model ++++++++
         nu = 0
         IF (BGK_MODEL_TYPE_INT == 0)  nu = SQRT(2.d0)*n_CELL*BGK_SIGMA*V_TH    ! [1/s]  ! Classical BGK
         IF (BGK_MODEL_TYPE_INT == 1)  nu = BGK_BG_DENS*BGK_SIGMA*V_TH          ! [1/s]  ! Collisions with background
         IF (BGK_MODEL_TYPE_INT == 2)  nu = 1.0d5                               ! [1/s]  ! HARD-CODED

         ! ++++++++++ Test particles for collisions +++++++++++
         DO JP = JP_START, JP_END

            IDp = IND(JP) ! ID of current particle

            ! Try collision - note that in this BGK model I use a collision frequency function of T:
            !                 collisions happen only depending on \nu and on a random number, independently
            !                 from the actual particle velocity.

            IF (rf() .LT. (1 - EXP(-nu*particles(IDp)%DTRIM))) THEN ! PROBABILITY OF COLLISION = 1 - exp(-nu*dt)
   
               IF (BGK_MODEL_TYPE_INT == 0) THEN ! Classical collisions, just pick from a Maxwellian at the  cell equil temp

                  ! Sample from Maxwellian (isotropic) at local velocity and temperature 
                  CALL MAXWELL(VX_CELL, VY_CELL, VZ_CELL, T_CELL, T_CELL, T_CELL, &
                              VX_NOW, VY_NOW, VZ_NOW, MASS)

               ELSE IF (BGK_MODEL_TYPE_INT == 1) THEN ! Collisions with background (heavy) particle

                  ! New particle velocity, obtained reducing the energy by the average factor 2 M/M_BG
                  V2_NOW = ( particles(IDp)%VX**2 + particles(IDp)%VY**2 + particles(IDp)%VZ**2 ) ! Initial V2
                  MASS   = SPECIES(particles(IDp)%S_ID)%MOLECULAR_MASS
                  V2_NOW = V2_NOW*(1 - 2*MASS/BGK_BG_MASS) ! Rescale it
   
                  ! Now sample new velocity on a sphere (see Bird for example)
                  Q      = 2*rf() - 1
                  COS_TH = Q
                  SIN_TH = SQRT(1 - Q**2)

                  CHI     = 2*pi*rf()
                  COS_CHI = COS(CHI)
                  SIN_CHI = SIN(CHI)

                  VX_NOW = SIN_TH*COS_CHI*SQRT(V2_NOW)
                  VY_NOW = SIN_TH*SIN_CHI*SQRT(V2_NOW)
                  VZ_NOW = COS_TH*SQRT(V2_NOW)

               ELSE IF (BGK_MODEL_TYPE_INT == 2) THEN ! Collisions with background (heavy) particles, relaxation to  
                                                      ! a Maxwellian with T_tot = T + m*u^2/(3*kB) 
                                                      ! (see Boccelli et al, Phys Plasmas, 2020 - submitted)

                  ! Sample from Maxwellian (isotropic) at zero and "total" temperature
                  CALL MAXWELL(0.0d0, 0.0d0, 0.0d0, Ttot_CELL, Ttot_CELL, Ttot_CELL, &
                              VX_NOW, VY_NOW, VZ_NOW, MASS)

               END IF

               ! Assign new velocities and internal energy to particle
               particles(IDp)%VX = VX_NOW
               particles(IDp)%VY = VY_NOW
               particles(IDp)%VZ = VZ_NOW
   
               CALL INTERNAL_ENERGY(SPECIES(1)%ROTDOF, Ttot_CELL, EI_NOW)
               particles(IDp)%EROT = EI_NOW

               CALL INTERNAL_ENERGY(SPECIES(1)%VIBDOF, Ttot_CELL, EI_NOW)
               particles(IDp)%EVIB = EI_NOW
   
               ! Update number of collisions happened
               TIMESTEP_COLL = TIMESTEP_COLL + 1
      
            END IF

         END DO

      
      END DO
   
      DEALLOCATE(IND)
      DEALLOCATE(NPC)
      DEALLOCATE(IOF)
   
   END SUBROUTINE BGK_COLLISIONS



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MCC_COLLISIONS -> perform MCC collisions !!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE MCC_COLLISIONS

      ! Performs MCC collisions: for each particle in the current process ID, computes the probability of 
      ! colliding with a background species, and (in case) performs the collision. 
      ! 
      ! There is no need to consider cells for these collisions. The neutrals density is for the
      ! moment uniform in space. In the future, we may associate it to the cells...
      ! NOTE: Eventually this procedure could be merged with similar procedures in collisions.f90

      IMPLICIT NONE

      INTEGER      :: JP1, JP2, JR, I, J, SP_ID1, SP_ID2, P1_SP_ID, P2_SP_ID, P3_SP_ID, NP_PROC_INITIAL
      REAL(KIND=8) :: P_COLL, PTCE, rfp, BG_NRHO
      REAL(KIND=8) :: SIGMA, OMEGA, CREF, ALPHA, FRAC, SIGMA_R
      REAL(KIND=8) :: PI2
      REAL(KIND=8), DIMENSION(3) :: C1, C2, GREL, W
      REAL(KIND=8) :: GX, GY, GZ, G
      REAL(KIND=8) :: COSCHI, SINCHI, THETA, COSTHETA, SINTHETA
      REAL(KIND=8) :: VR2, VR, MRED, M1, M2, COSA, SINA, BB
      REAL(KIND=8) :: TRDOF, PROT1, PROT2, PVIB1, PVIB2, EI, ETR, ECOLL, TOTDOF, EA, EROT, EVIB
      LOGICAL :: SKIP
      TYPE(PARTICLE_DATA_STRUCTURE) :: NEWparticle

      PI2 = 2*PI

      TIMESTEP_COLL = 0
      TIMESTEP_REAC = 0

      NP_PROC_INITIAL = NP_PROC
      DO JP1 = 1,NP_PROC_INITIAL
         SP_ID1 = particles(JP1)%S_ID

         DO J = 1, MIXTURES(MCC_BG_MIX)%N_COMPONENTS
            SP_ID2 = MIXTURES(MCC_BG_MIX)%COMPONENTS(J)%ID

            IF (BOOL_BG_DENSITY_FILE) THEN
               BG_NRHO = MCC_BG_CELL_NRHO(SP_ID2, particles(JP1)%IC)
            ELSE
               FRAC = MIXTURES(MCC_BG_MIX)%COMPONENTS(J)%MOLFRAC
               BG_NRHO = FRAC*MCC_BG_DENS
            END IF
            IF (BG_NRHO == 0) CYCLE
            
            !SIGMA = PI * (0.5 * (SPECIES(SP_ID1)%DIAM + SPECIES(SP_ID2)%DIAM))**2
            !OMEGA = 0.5 * (SPECIES(SP_ID1)%OMEGA + SPECIES(SP_ID2)%OMEGA)
            CREF = VSS_GREFS(SP_ID1, SP_ID2)
            !ALPHA = 0.5 * (SPECIES(SP_ID1)%ALPHA + SPECIES(SP_ID2)%ALPHA)
            SIGMA = VSS_SIGMAS(SP_ID1, SP_ID2)
            OMEGA = VSS_OMEGAS(SP_ID1, SP_ID2)
            ALPHA = VSS_ALPHAS(SP_ID1, SP_ID2)
            
            ! Sample the velocity of second collision partner, that would be in the MCC backgorund
            ! This may bias collisions towards higher temperatures of particle 2!
            CALL MAXWELL(0.d0, 0.d0, 0.d0, &
            MCC_BG_TTRA, MCC_BG_TTRA, MCC_BG_TTRA, &
            C2(1), C2(2), C2(3), SPECIES(SP_ID2)%MOLECULAR_MASS)

            ! Compute the real relative velocity
            VR2 = (C2(1) - particles(JP1)%VX)**2 + &
            (C2(2) - particles(JP1)%VY)**2 + &
            (C2(3) - particles(JP1)%VZ)**2 

            VR = SQRT(VR2)
            

            ! Compute collision probability
            P_COLL = 1 - exp(-DT*BG_NRHO*VR*SIGMA*(VR/CREF)**(1.-2.*OMEGA))

            ! Try the collision
            IF (rf() < P_COLL) THEN ! Collision happens

               TIMESTEP_COLL = TIMESTEP_COLL + 1

               ! Rimuovere commento per avere avviso
               IF (P_COLL .GT. 1.) THEN
                  WRITE(*,*) 'Attention => this was a bad MCC collision!!!'
               END IF

               C1(1) = particles(JP1)%VX
               C1(2) = particles(JP1)%VY
               C1(3) = particles(JP1)%VZ

               ! Actually create the second collision partner
               CALL INTERNAL_ENERGY(SPECIES(SP_ID2)%ROTDOF, MCC_BG_TTRA, EROT)
               CALL INTERNAL_ENERGY(SPECIES(SP_ID2)%VIBDOF, MCC_BG_TTRA, EVIB)
               CALL INIT_PARTICLE(particles(JP1)%X,particles(JP1)%Y,particles(JP1)%Z, &
               C2(1),C2(2),C2(3),EROT,EVIB,SP_ID2,particles(JP1)%IC,DT, NEWparticle)
               !WRITE(*,*) 'Should be adding particle!'
               CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)
               
               JP2 = NP_PROC

               M1    = SPECIES(SP_ID1)%MOLECULAR_MASS
               M2    = SPECIES(SP_ID2)%MOLECULAR_MASS
               MRED  = M1*M2/(M1+M2)
      
               ETR = 0.5*MRED*VR2

               ! Test for chemical reaction with TCE model
               ECOLL = ETR + particles(JP1)%EROT + particles(JP2)%EROT + &
               particles(JP1)%EVIB + particles(JP2)%EVIB
               SKIP = .FALSE.
               DO JR = 1, N_REACTIONS
                  !WRITE(*,*) 'Checking reaction with R1 = ', REACTIONS(JR)%R1_SP_ID, 'R2 = ', REACTIONS(JR)%R2_SP_ID, &
                  !' against particles with ids ',  SP_ID1, SP_ID2
                  ! Check if species in collision belong to this reaction.
                  ! If they don't, check the next reaction.
                  ! Not very efficient with many reactions
                  IF (REACTIONS(JR)%R1_SP_ID .NE. SP_ID1) THEN
                     IF (REACTIONS(JR)%R1_SP_ID .NE. SP_ID2) THEN
                        CYCLE
                     ELSE
                        IF (REACTIONS(JR)%R2_SP_ID .NE. SP_ID1) CYCLE
                     END IF
                  ELSE
                     IF (REACTIONS(JR)%R2_SP_ID .NE. SP_ID2) CYCLE
                  END IF
      
                  EA = REACTIONS(JR)%EA
                  IF (ECOLL .LE. EA) CYCLE

                  IF (REACTIONS(JR)%TYPE == TCE) THEN
                     IF (ECOLL .LE. EA) CYCLE
                     PTCE = REACTIONS(JR)%C1 * (ECOLL-EA)**REACTIONS(JR)%C2 * (1.-EA/ECOLL)**REACTIONS(JR)%C3
                  ELSE IF (REACTIONS(JR)%TYPE == LXCAT) THEN
                     SIGMA_R = INTERP_CS(ETR, REACTIONS(JR)%TABLE_ENERGY, REACTIONS(JR)%TABLE_CS)
                     PTCE = SIGMA_R / (SIGMA*(VR/CREF)**(1.-2.*OMEGA))
                  ELSE
                     PTCE = 0
                  END IF
                  rfp = rf()
                  ! Here we suppose that the probability is so low that we can test sequentially with acceptable error
                  IF (rfp .LT. PTCE) THEN
                     SKIP = .TRUE.
                     TIMESTEP_REAC = TIMESTEP_REAC + 1
                     !WRITE(*,*) 'Reacting!'
                     ! React
                     ECOLL = ECOLL - EA

                     ! Use R1 as P1 and assign internal energy to it.
                     ! Use R2 as P2 (or P2+P3) (Set R2 as the molecule that dissociates in P2+P3)
                     P1_SP_ID = REACTIONS(JR)%P1_SP_ID
                     P2_SP_ID = REACTIONS(JR)%P2_SP_ID

                     particles(JP1)%S_ID = REACTIONS(JR)%P1_SP_ID
                     particles(JP2)%S_ID = REACTIONS(JR)%P2_SP_ID

                     IF (.NOT. REACTIONS(JR)%IS_CEX) THEN
         
                        TOTDOF = 3. + SPECIES(P2_SP_ID)%ROTDOF + SPECIES(P2_SP_ID)%VIBDOF + &
                              SPECIES(P1_SP_ID)%ROTDOF
                        IF (REACTIONS(JR)%N_PROD == 3) THEN
                           P3_SP_ID = REACTIONS(JR)%P3_SP_ID
                           TOTDOF = TOTDOF + 3. + SPECIES(P3_SP_ID)%ROTDOF + SPECIES(P3_SP_ID)%VIBDOF
                        END IF
         
                        EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P1_SP_ID)%VIBDOF)
                        particles(JP1)%EVIB = EI
                        ECOLL = ECOLL - EI
         
                        TOTDOF = TOTDOF - SPECIES(P1_SP_ID)%ROTDOF
                        EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P1_SP_ID)%ROTDOF)
                        particles(JP1)%EROT = EI
                        ECOLL = ECOLL - EI
         
                        TOTDOF = TOTDOF - SPECIES(P2_SP_ID)%VIBDOF
                        EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P2_SP_ID)%VIBDOF)
                        particles(JP2)%EVIB = EI
                        ECOLL = ECOLL - EI
         
                        TOTDOF = TOTDOF - SPECIES(P2_SP_ID)%ROTDOF
                        EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P2_SP_ID)%ROTDOF)
                        particles(JP2)%EROT = EI
                        ECOLL = ECOLL - EI
         
                        M1 = SPECIES(P1_SP_ID)%MOLECULAR_MASS
                        IF (REACTIONS(JR)%N_PROD == 2) THEN
                           M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                        ELSE IF (REACTIONS(JR)%N_PROD == 3) THEN
                           M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS + SPECIES(P3_SP_ID)%MOLECULAR_MASS
                        END IF
                        TOTDOF = TOTDOF - 3.
                        EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, 3)
                        CALL HS_SCATTER(EI, M1, M2, C1, C2)
                        ECOLL = ECOLL - EI
         
                        particles(JP1)%VX = C1(1)
                        particles(JP1)%VY = C1(2)
                        particles(JP1)%VZ = C1(3)

         
                        IF (REACTIONS(JR)%N_PROD == 2) THEN
                           particles(JP2)%VX = C2(1)
                           particles(JP2)%VY = C2(2)
                           particles(JP2)%VZ = C2(3)
                        ELSE IF (REACTIONS(JR)%N_PROD == 3) THEN
                           TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%VIBDOF
                           EVIB = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%VIBDOF)
                           ECOLL = ECOLL - EVIB
               
                           TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%ROTDOF
                           EROT = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%ROTDOF)
                           ECOLL = ECOLL - EROT
            
                           M1 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                           M2 = SPECIES(P3_SP_ID)%MOLECULAR_MASS
                           C1 = C2
                           
                           CALL HS_SCATTER(ECOLL, M1, M2, C1, C2)
                                 
                           particles(JP2)%VX = C1(1)
                           particles(JP2)%VY = C1(2)
                           particles(JP2)%VZ = C1(3)
                                 
                           
                           CALL INIT_PARTICLE(particles(JP2)%X,particles(JP2)%Y,particles(JP2)%Z, &
                           C2(1),C2(2),C2(3),EROT,EVIB,P3_SP_ID,particles(JP2)%IC,DT, NEWparticle)
                           !WRITE(*,*) 'Should be adding particle!'
                           CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)
                           
                        
                        END IF
                     END IF

                     ! If this pair had a chemical reaction, exit and don't test any other reaction.
                     EXIT

                  END IF
               END DO
      
               IF (SKIP) CYCLE
               ! Perform elastic/inelastic collision
               ! Test for inelastic collision and TR/TV exchange
               TRDOF = 3.
      
               !TRDOF = 5. -2.*OMEGA why not?
      
               PROT1 = (TRDOF + SPECIES(SP_ID1)%ROTDOF)/TRDOF * SPECIES(SP_ID1)%ROTREL
               PROT2 = (TRDOF + SPECIES(SP_ID2)%ROTDOF)/TRDOF * SPECIES(SP_ID2)%ROTREL
               PVIB1 = (TRDOF + SPECIES(SP_ID1)%VIBDOF)/TRDOF * SPECIES(SP_ID1)%VIBREL
               PVIB2 = (TRDOF + SPECIES(SP_ID2)%VIBDOF)/TRDOF * SPECIES(SP_ID2)%VIBREL
      
               TRDOF = 5. -2.*OMEGA
      
               rfp = rf()
               IF (rfp .LT. PROT1) THEN
                  ! Particle 1 selected for rotational energy exchange
                  ECOLL = ETR + particles(JP1)%EROT
                  particles(JP1)%EROT = COLL_INTERNAL_ENERGY(ECOLL, TRDOF, SPECIES(SP_ID1)%ROTDOF)
                  ETR = ECOLL - particles(JP1)%EROT
      
               ELSE IF (rfp .LT. PROT1 + PROT2) THEN
                  ! Particle 2 selected for rotational energy exchange
                  ECOLL = ETR + particles(JP2)%EROT
                  particles(JP2)%EROT = COLL_INTERNAL_ENERGY(ECOLL, TRDOF, SPECIES(SP_ID2)%ROTDOF)
                  ETR = ECOLL - particles(JP2)%EROT
      
               ELSE IF (rfp .LT. PROT1 + PROT2 + PVIB1) THEN
                  ! Particle 1 selected for vibrational energy exchange
                  ECOLL = ETR + particles(JP1)%EVIB
                  particles(JP1)%EVIB = COLL_INTERNAL_ENERGY(ECOLL, TRDOF, SPECIES(SP_ID1)%VIBDOF)
                  ETR = ECOLL - particles(JP1)%EVIB
      
               ELSE IF (rfp .LT. PROT1 + PROT2 + PVIB1 + PVIB2) THEN
                  ! Particle 2 selected for vibrational energy exchange
                  ECOLL = ETR + particles(JP2)%EVIB
                  particles(JP2)%EVIB = COLL_INTERNAL_ENERGY(ECOLL, TRDOF, SPECIES(SP_ID2)%VIBDOF)
                  ETR = ECOLL - particles(JP2)%EVIB
      
               END IF
      
               ! Compute post-collision velocities
               ! VSS Collisions
      
               COSCHI = 2.*rf()**(1./ALPHA) - 1.
               SINCHI = SQRT(1.-COSCHI*COSCHI)
               THETA = PI2*rf()
               COSTHETA = COS(THETA)
               SINTHETA = SIN(THETA)
      
      
               IF (rfp .LT. PROT1 + PROT2 + PVIB1 + PVIB2) THEN
      
                  ! Inelastic collision occurred => Randomize relative velocity vector
                  G = SQRT(2.*ETR/MRED)
                  
                  COSA = 2.*rf()-1.0
                  SINA = SQRT(1-COSA*COSA)
                  BB = PI2 * rf()
      
                  GX = G*SINA*COS(BB)
                  GY = G*SINA*SIN(BB)
                  GZ = G*COSA
               
               ELSE
      
                  ! No inelastic collision occured => Take actual velocity vector
                  GX = C1(1) - C2(1)
                  GY = C1(2) - C2(2)
                  GZ = C1(3) - C2(3)
                  G = SQRT(GX*GX + GY*GY + GZ*GZ)
      
               END IF
      
               ! Compute post-collision velocities
      
               GREL(1) = GX*COSCHI + SQRT(GY*GY+GZ*GZ)*SINTHETA*SINCHI
               GREL(2) = GY*COSCHI + (G*GZ*COSTHETA - GX*GY*SINTHETA)/SQRT(GY*GY+GZ*GZ)*SINCHI
               GREL(3) = GZ*COSCHI - (G*GY*COSTHETA + GX*GZ*SINTHETA)/SQRT(GY*GY+GZ*GZ)*SINCHI
      
               ! Compute center of mass velocity vector
               DO I = 1, 3
                  W(I) = M1/(M1+M2)*C1(I) + M2/(M1+M2)*C2(I)
               END DO
      
               ! Compute absolute velocity vector of the two particles
               DO I = 1, 3
                  C1(I) = W(I) + M2/(M1+M2)*GREL(I)
                  C2(I) = W(I) - M1/(M1+M2)*GREL(I)
               END DO
      
               
               particles(JP1)%VX = C1(1)
               particles(JP1)%VY = C1(2)
               particles(JP1)%VZ = C1(3)
                           
               particles(JP2)%VX = C2(1)
               particles(JP2)%VY = C2(2)
               particles(JP2)%VZ = C2(3)
               
            END IF
         END DO
      END DO

   END SUBROUTINE MCC_COLLISIONS





   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MCC_COLLISIONS_VAHEDI -> perform MCC collisions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Implemented in the style of Vahedi and Surendra, Computer Physics Communications 87 (1995) 179-198 !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE MCC_COLLISIONS_VAHEDI

      ! Performs MCC collisions: for each particle in the current process ID, computes the probability of 
      ! colliding with a background species, and (in case) performs the collision. 
      ! 
      ! There is no need to consider cells for these collisions. The neutrals density is for the
      ! moment uniform in space. In the future, we may associate it to the cells...
      ! NOTE: Eventually this procedure could be merged with similar procedures in collisions.f90

      IMPLICIT NONE

      INTEGER      :: JP1, JP2, JP3, JR, J, SP_ID1, SP_ID2, P1_SP_ID, P2_SP_ID, P3_SP_ID, P4_SP_ID, NP_PROC_INITIAL
      REAL(KIND=8) :: BG_NRHO
      REAL(KIND=8) :: FRAC, SIGMA_R
      REAL(KIND=8) :: PI2
      REAL(KIND=8), DIMENSION(3) :: C1, C2
      REAL(KIND=8) :: VR2, VR, MRED, M1, M2
      REAL(KIND=8) :: EI, ETR, ECOLL, TOTDOF, EA, EROT, EVIB
      TYPE(PARTICLE_DATA_STRUCTURE) :: NEWparticle
      REAL(KIND=8) :: NULL_COLL_FREQ, P_NULL, R_SELECT, P_CUMULATED
      LOGICAL :: HAS_REACTED

      PI2 = 2*PI

      TIMESTEP_COLL = 0
      TIMESTEP_REAC = 0

      DO JR = 1, N_REACTIONS
         REACTIONS(JR)%COUNTS = 0
      END DO

      NULL_COLL_FREQ = MCC_BG_DENS*1e-12
      P_NULL = 1 - EXP(-DT*NULL_COLL_FREQ)
      !P_NULL = DT*NULL_COLL_FREQ
      IF (P_NULL > 1) THEN
         WRITE(*,*) 'Caution! P_NULL is greater than 1! P_NULL = ', P_NULL, 'Setting it to 1 and adjusting the collision frequency.'
         P_NULL = 1.
         NULL_COLL_FREQ = 1./DT
      END IF

      NP_PROC_INITIAL = NP_PROC
      DO JP1 = 1, NP_PROC_INITIAL
         HAS_REACTED = .FALSE.
         IF (rf() > P_NULL) CYCLE
         TIMESTEP_COLL = TIMESTEP_COLL + 1

         SP_ID1 = particles(JP1)%S_ID

         DO J = 1, MIXTURES(MCC_BG_MIX)%N_COMPONENTS
            SP_ID2 = MIXTURES(MCC_BG_MIX)%COMPONENTS(J)%ID

            FRAC = MIXTURES(MCC_BG_MIX)%COMPONENTS(J)%MOLFRAC
            BG_NRHO = FRAC*MCC_BG_DENS

            IF (BG_NRHO == 0) CYCLE
            
            ! Sample the velocity of second collision partner, that would be in the MCC backgorund
            ! This may bias collisions towards higher temperatures of particle 2!
            CALL MAXWELL(0.d0, 0.d0, 0.d0, &
            MCC_BG_TTRA, MCC_BG_TTRA, MCC_BG_TTRA, &
            C2(1), C2(2), C2(3), SPECIES(SP_ID2)%MOLECULAR_MASS)

            ! Compute the real relative velocity
            VR2 = (C2(1) - particles(JP1)%VX)**2 + &
                  (C2(2) - particles(JP1)%VY)**2 + &
                  (C2(3) - particles(JP1)%VZ)**2

            VR = SQRT(VR2)
            

            ! Compute collision probability
            R_SELECT = rf()

            C1(1) = particles(JP1)%VX
            C1(2) = particles(JP1)%VY
            C1(3) = particles(JP1)%VZ

            M1    = SPECIES(SP_ID1)%MOLECULAR_MASS
            M2    = SPECIES(SP_ID2)%MOLECULAR_MASS
            MRED  = M1*M2/(M1+M2)
   
            ! Compute the kinetic energy in the center-of-mass frame (aka collision energy)
            ETR = 0.5*MRED*VR2

            P_CUMULATED = 0

            DO JR = 1, N_REACTIONS
               !WRITE(*,*) 'Checking reaction with R1 = ', REACTIONS(JR)%R1_SP_ID, 'R2 = ', REACTIONS(JR)%R2_SP_ID, &
               !' against particles with ids ',  SP_ID1, SP_ID2
               ! Check if species in collision belong to this reaction.
               ! If they don't, check the next reaction.
               ! Not very efficient with many reactions
               IF ((REACTIONS(JR)%R1_SP_ID .NE. SP_ID1) .OR. (REACTIONS(JR)%R2_SP_ID .NE. SP_ID2)) CYCLE
   
               EA = REACTIONS(JR)%EA
               IF (ETR .LE. EA) CYCLE

               IF (REACTIONS(JR)%TYPE == LXCAT) THEN
                  SIGMA_R = INTERP_CS(ETR, REACTIONS(JR)%TABLE_ENERGY, REACTIONS(JR)%TABLE_CS)
                  !P_CUMULATED = P_CUMULATED + BG_NRHO*SIGMA_R*VR / NULL_COLL_FREQ
                  P_CUMULATED = P_CUMULATED + (1. -EXP(-BG_NRHO*SIGMA_R*VR*DT)) / P_NULL
               ELSE
                  CYCLE
               END IF

               ! Try the reaction
               IF (R_SELECT < P_CUMULATED) THEN ! Collision happens
                  REACTIONS(JR)%COUNTS = REACTIONS(JR)%COUNTS + 1

                  ! Actually create the second collision partner
                  CALL INTERNAL_ENERGY(SPECIES(SP_ID2)%ROTDOF, MCC_BG_TTRA, EROT)
                  CALL INTERNAL_ENERGY(SPECIES(SP_ID2)%VIBDOF, MCC_BG_TTRA, EVIB)
                  CALL INIT_PARTICLE(particles(JP1)%X,particles(JP1)%Y,particles(JP1)%Z, &
                  C2(1),C2(2),C2(3),EROT,EVIB,SP_ID2,particles(JP1)%IC,DT, NEWparticle)

                  ! Add particle to array!
                  CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)
                  JP2 = NP_PROC


                  TIMESTEP_REAC = TIMESTEP_REAC + 1

                  ! Rimuovere commento per avere avviso
                  IF (P_CUMULATED .GT. 1.) THEN
                     WRITE(*,*) 'Attention => this was a bad MCC collision!!!'
                  END IF

                  !WRITE(*,*) 'Reacting!'
                  ! React
                  ECOLL = ETR - EA

                  ! Use R1 as P1 and assign internal energy to it.
                  ! Use R2 as P2 (or P2+P3) (Set R2 as the molecule that dissociates in P2+P3)
                  P1_SP_ID = REACTIONS(JR)%P1_SP_ID
                  P2_SP_ID = REACTIONS(JR)%P2_SP_ID

                  !WRITE(*,*) 'SP_ID changed from ', particles(JP1)%S_ID, particles(JP2)%S_ID, &
                  !' to ', REACTIONS(JR)%P1_SP_ID, REACTIONS(JR)%P2_SP_ID

                  particles(JP1)%S_ID = REACTIONS(JR)%P1_SP_ID
                  particles(JP2)%S_ID = REACTIONS(JR)%P2_SP_ID


                  IF (.NOT. REACTIONS(JR)%IS_CEX) THEN
      
                     TOTDOF = 3. + SPECIES(P2_SP_ID)%ROTDOF + SPECIES(P2_SP_ID)%VIBDOF + &
                           SPECIES(P1_SP_ID)%ROTDOF
                     IF (REACTIONS(JR)%N_PROD .GE. 3) THEN
                        P3_SP_ID = REACTIONS(JR)%P3_SP_ID
                        TOTDOF = TOTDOF + 3. + SPECIES(P3_SP_ID)%ROTDOF + SPECIES(P3_SP_ID)%VIBDOF
                     END IF
                     IF (REACTIONS(JR)%N_PROD .GE. 4) THEN
                        P4_SP_ID = REACTIONS(JR)%P4_SP_ID
                        TOTDOF = TOTDOF + 3. + SPECIES(P4_SP_ID)%ROTDOF + SPECIES(P4_SP_ID)%VIBDOF
                     END IF
      
                     EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P1_SP_ID)%VIBDOF)
                     particles(JP1)%EVIB = EI
                     ECOLL = ECOLL - EI
      
                     TOTDOF = TOTDOF - SPECIES(P1_SP_ID)%ROTDOF
                     EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P1_SP_ID)%ROTDOF)
                     particles(JP1)%EROT = EI
                     ECOLL = ECOLL - EI
      
                     TOTDOF = TOTDOF - SPECIES(P2_SP_ID)%VIBDOF
                     EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P2_SP_ID)%VIBDOF)
                     particles(JP2)%EVIB = EI
                     ECOLL = ECOLL - EI
      
                     TOTDOF = TOTDOF - SPECIES(P2_SP_ID)%ROTDOF
                     EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P2_SP_ID)%ROTDOF)
                     particles(JP2)%EROT = EI
                     ECOLL = ECOLL - EI
      
                     M1 = SPECIES(P1_SP_ID)%MOLECULAR_MASS
                     IF (REACTIONS(JR)%N_PROD == 2) THEN
                        M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                     ELSE IF (REACTIONS(JR)%N_PROD == 3) THEN
                        M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS + SPECIES(P3_SP_ID)%MOLECULAR_MASS
                     ELSE IF (REACTIONS(JR)%N_PROD == 4) THEN
                        M2 = SPECIES(P2_SP_ID)%MOLECULAR_MASS + SPECIES(P3_SP_ID)%MOLECULAR_MASS &
                           + SPECIES(P4_SP_ID)%MOLECULAR_MASS
                     END IF
                     TOTDOF = TOTDOF - 3.
                     EI = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, 3)
                     CALL HS_SCATTER(EI, M1, M2, C1, C2)
                     ECOLL = ECOLL - EI
      
                     particles(JP1)%VX = C1(1)
                     particles(JP1)%VY = C1(2)
                     particles(JP1)%VZ = C1(3)

      
                     IF (REACTIONS(JR)%N_PROD == 2) THEN
                        particles(JP2)%VX = C2(1)
                        particles(JP2)%VY = C2(2)
                        particles(JP2)%VZ = C2(3)
                     ELSE IF (REACTIONS(JR)%N_PROD == 3) THEN
                        TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%VIBDOF
                        EVIB = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%VIBDOF)
                        ECOLL = ECOLL - EVIB
            
                        TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%ROTDOF
                        EROT = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%ROTDOF)
                        ECOLL = ECOLL - EROT
         
                        M1 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                        M2 = SPECIES(P3_SP_ID)%MOLECULAR_MASS
                        C1 = C2
                        
                        CALL HS_SCATTER(ECOLL, M1, M2, C1, C2)
                              
                        particles(JP2)%VX = C1(1)
                        particles(JP2)%VY = C1(2)
                        particles(JP2)%VZ = C1(3)
                              
                        
                        CALL INIT_PARTICLE(particles(JP2)%X,particles(JP2)%Y,particles(JP2)%Z, &
                        C2(1),C2(2),C2(3),EROT,EVIB,P3_SP_ID,particles(JP2)%IC,DT, NEWparticle)
                        !WRITE(*,*) 'Should be adding particle!'
                        CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)
                        
                     ELSE IF (REACTIONS(JR)%N_PROD == 4) THEN
                        TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%VIBDOF
                        EVIB = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%VIBDOF)
                        ECOLL = ECOLL - EVIB
            
                        TOTDOF = TOTDOF - SPECIES(P3_SP_ID)%ROTDOF
                        EROT = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P3_SP_ID)%ROTDOF)
                        ECOLL = ECOLL - EROT
         
                        M1 = SPECIES(P2_SP_ID)%MOLECULAR_MASS
                        M2 = SPECIES(P3_SP_ID)%MOLECULAR_MASS + SPECIES(P4_SP_ID)%MOLECULAR_MASS
                        C1 = C2
                        
                        CALL HS_SCATTER(ECOLL, M1, M2, C1, C2)
                              
                        particles(JP2)%VX = C1(1)
                        particles(JP2)%VY = C1(2)
                        particles(JP2)%VZ = C1(3)
                              
                        
                        CALL INIT_PARTICLE(particles(JP2)%X,particles(JP2)%Y,particles(JP2)%Z, &
                        C2(1),C2(2),C2(3),EROT,EVIB,P3_SP_ID,particles(JP2)%IC,DT, NEWparticle)
                        !WRITE(*,*) 'Should be adding particle!'
                        CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)
                        JP3 = NP_PROC

                        TOTDOF = TOTDOF - SPECIES(P4_SP_ID)%VIBDOF
                        EVIB = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P4_SP_ID)%VIBDOF)
                        ECOLL = ECOLL - EVIB
            
                        TOTDOF = TOTDOF - SPECIES(P4_SP_ID)%ROTDOF
                        EROT = COLL_INTERNAL_ENERGY(ECOLL, TOTDOF, SPECIES(P4_SP_ID)%ROTDOF)
                        ECOLL = ECOLL - EROT
         
                        M1 = SPECIES(P3_SP_ID)%MOLECULAR_MASS
                        M2 = SPECIES(P4_SP_ID)%MOLECULAR_MASS
                        C1 = C2
                        
                        CALL HS_SCATTER(ECOLL, M1, M2, C1, C2)
                              
                        particles(JP3)%VX = C1(1)
                        particles(JP3)%VY = C1(2)
                        particles(JP3)%VZ = C1(3)
                              
                        
                        CALL INIT_PARTICLE(particles(JP3)%X,particles(JP3)%Y,particles(JP3)%Z, &
                        C2(1),C2(2),C2(3),EROT,EVIB,P3_SP_ID,particles(JP3)%IC,DT, NEWparticle)
                        !WRITE(*,*) 'Should be adding particle!'
                        CALL ADD_PARTICLE_ARRAY(NEWparticle, NP_PROC, particles)
                        
                     END IF
                  END IF

                  ! If this pair had a chemical reaction, exit and don't test any other reaction.
                  HAS_REACTED = .TRUE.
                  EXIT
               END IF
            END DO
            IF (HAS_REACTED) EXIT
         END DO

      END DO

   END SUBROUTINE MCC_COLLISIONS_VAHEDI


END MODULE collisions
