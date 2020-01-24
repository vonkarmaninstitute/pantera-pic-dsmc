  ! Contains subroutines related to particle collisions

  MODULE collisions

  USE global
  USE mpi_common
  USE screen
  USE tools
  
  IMPLICIT NONE
  
  CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE DSMC_COLLISIONS -> Calcola collisioni !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSMC_COLLISIONS 

   ! Reorder particles
   ! IND => index of particle in particle array, indexed by particles I own and sorted by cells
   ! NPC => number of particles in cell, indexed by cell index (a block length for IND) 
   ! could this be indexed by all cells? for those that I don't own, NPC will turn out to be 0.
   ! + would simplify significantly,
   ! - could be up to about 50 times bigger
   ! + it's just a list of integers
   ! --->  I am gonna try and in case revert back.
   ! IOF => index (for IND) of first particle in a cell (an offset in IND)
   
   INTEGER, DIMENSION(NX*NY)          :: NPC, IOF
   INTEGER, DIMENSION(:), ALLOCATABLE :: IND
   INTEGER                            :: JP, JC, NCELLS, IDX
   INTEGER                            :: NCOLLREAL

   ! Count the number of particles in each cell, to allocate arrays later
   NPC = 0
   DO JP = 1, NP_PROC
      JC = particles(JP)%IC
      NPC(JC) = NPC(JC) + 1
   END DO

   NCELLS = NX*NY

   ! Fill the array of offsets (IOF). IOF(IC) is the the index of the first
   ! particle in cell IP
   IOF = -1
   IDX = 1
   DO JC = 1, NCELLS
      IF (NPC(JC) .NE. 0) THEN
         IOF(JC) = IDX
         IDX = IDX + NPC(JC)
      END IF
   END DO

   ! Allocate and fill the array of indices (IND). IND(IP) is the particle index (for the "particles" array)
   ! but ordered by cells, with offsets given by IND(IC) and stride length NPC(IC)
   ALLOCATE(IND(NP_PROC))
 
   NPC = 0
   DO JP = 1, NP_PROC
      JC = particles(JP)%IC
      IND(IOF(JC) + NPC(JC)) = JP
      NPC(JC) = NPC(JC) + 1
   END DO
 
   ! Compute collisions between particles
   TIMESTEP_COLL = 0
   DO JC = 1, NCELLS
      IF (NPC(JC) .GT. 1) THEN
         CALL VSS_COLLIS(JC, NPC, IOF, IND, NCOLLREAL)
         ! Add to the total number of collisions for this process
         TIMESTEP_COLL = TIMESTEP_COLL + NCOLLREAL
      END IF
   END DO
 
   !WRITE(*,*) 'Number of real collisions: ', TIMESTEP_COLL

   DEALLOCATE(IND)
 
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


  INTEGER      :: IOFJ,IOLJ,NPCJ,JP1,JP2,JCOL, JP, INDJ, NSP, JSP, I
  INTEGER      :: SP_ID1, SP_ID2
  INTEGER      :: NCOLL,NCOLLMAX_INT
  REAL(KIND=8) :: SIGMA, ALPHA, SIGMAMAX, NCOLLMAX,FCORR,VR,VR2,rfp
  REAL(KIND=8) :: OMEGA, CREF, TREF, ZETA, MRED, COLLPROB
  REAL(KIND=8) :: TRDOF, PROT1, PROT2, PVIB1, PVIB2, EI, ETR, ECOLL, R1, R2
  !REAL(KIND=8) :: B,C,EINT,ETOT,ETR,PHI,SITETA,VRX,VRY,VRZ
  REAL(KIND=8) :: VXMAX,VXMIN,VYMAX,VYMIN,VZMAX,VZMIN,VRMAX
  REAL(KIND=8) :: PI, PI2, KB
  REAL(KIND=8), DIMENSION(3) :: C1, C2, GREL, W
  REAL(KIND=8) :: GX, GY, GZ, G
  REAL(KIND=8) :: COSCHI, SINCHI, THETA, COSTHETA, SINTHETA
  REAL(KIND=8) :: M1, M2, COSA, SINA, BB
  
  PI   = 3.141593
  PI2  = 2.*PI
  
  KB = 1.38064852E-23

  IOFJ  = IOF(JC)               ! First particle in the cell
  IOLJ  = IOF(JC) + NPC(JC) - 1 ! Last particle in the cell
  NPCJ = NPC(JC)                ! Number of particles in the cell

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

  ! Find the maximum cross section of all the species involved in the collisions
  NSP = MIXTURES(DSMC_COLL_MIX)%N_COMPONENTS
  SIGMAMAX = 0
  DO JSP = 1, NSP
    SIGMA = SPECIES( MIXTURES(DSMC_COLL_MIX)%COMPONENTS(JSP)%ID )%SIGMA
    IF (SIGMA.GT.SIGMAMAX) SIGMAMAX = SIGMA
  END DO

  ! Compute the "worst case scenario" relative velocity
  VRMAX    = SQRT((VXMAX-VXMIN)**2 + (VYMAX-VYMIN)**2 + (VZMAX-VZMIN)**2)
  ! Compute the maximum expected number of collisions
  NCOLLMAX = 0.5*NPC(JC)*(NPC(JC)-1)*SIGMAMAX*VRMAX*FNUM*DT/CELL_VOL
  
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

     SIGMA = PI * (0.5 * (SPECIES(SP_ID1)%DIAM + SPECIES(SP_ID2)%DIAM))**2
     OMEGA    = 0.5 * (SPECIES(SP_ID1)%OMEGA + SPECIES(SP_ID2)%OMEGA)
     TREF = 0.5 * (SPECIES(SP_ID1)%TREF + SPECIES(SP_ID2)%TREF)
     CREF = GREFS(SP_ID1, SP_ID2)
     
     ALPHA = 0.5 * (SPECIES(SP_ID1)%ALPHA + SPECIES(SP_ID2)%ALPHA)
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
     ZETA = 0.

     ! Compute (simulated) collision probability
     IF (OMEGA .EQ. 0.5) THEN
       COLLPROB = FCORR/(SIGMAMAX*VRMAX)*VR*SIGMA
     ELSE
       COLLPROB = FCORR/(SIGMAMAX*VRMAX)*VR*SIGMA*(VR/CREF)**(1.-2.*OMEGA)
     END IF


     IF ( rfp .LT. COLLPROB ) THEN

        ! Rimuovere commento per avere avviso
        IF (COLLPROB .GT. 1.) THEN
           WRITE(*,*) 'Attention => this was a bad collision!!!'
        END IF

        NCOLLREAL = NCOLLREAL + 1



        ! Select elastic/inelastic collision 
        TRDOF = 3.

        !TRDOF = 5. -2.*OMEGA

        PROT1 = (TRDOF + SPECIES(SP_ID1)%ROTDOF)/TRDOF * SPECIES(SP_ID1)%ROTREL
        PROT2 = (TRDOF + SPECIES(SP_ID2)%ROTDOF)/TRDOF * SPECIES(SP_ID2)%ROTREL
        PVIB1 = (TRDOF + SPECIES(SP_ID1)%VIBDOF)/TRDOF * SPECIES(SP_ID1)%VIBREL
        PVIB2 = (TRDOF + SPECIES(SP_ID2)%VIBDOF)/TRDOF * SPECIES(SP_ID2)%VIBREL

        TRDOF = 5. -2.*OMEGA
        rfp = rf()

        IF (rfp .LT. PROT1) THEN
           ! Particle 1 selected for rotational energy exchange
           ECOLL = ETR + particles(JP1)%EROT
           DO
              R1 = rf()
              EI = R1*ECOLL
              R2 = rf()
              IF (R2 .LE. FUNCTION_I(ECOLL, EI, TRDOF, SPECIES(SP_ID1)%ROTDOF)) EXIT
           END DO
           particles(JP1)%EROT = EI
           ETR = ECOLL - EI

        ELSE IF (rfp .LT. PROT1 + PROT2) THEN
         ! Particle 2 selected for rotational energy exchange
         ECOLL = ETR + particles(JP2)%EROT
           DO
              R1 = rf()
              EI = R1*ECOLL
              R2 = rf()
              IF (R2 .LE. FUNCTION_I(ECOLL, EI, TRDOF, SPECIES(SP_ID2)%ROTDOF)) EXIT
           END DO
           particles(JP2)%EROT = EI
           ETR = ECOLL - EI

        ELSE IF (rfp .LT. PROT1 + PROT2 + PVIB1) THEN
         ! Particle 1 selected for vibrational energy exchange
         ECOLL = ETR + particles(JP1)%EVIB
         DO
            R1 = rf()
            EI = R1*ECOLL
            R2 = rf()
            IF (R2 .LE. FUNCTION_I(ECOLL, EI, TRDOF, SPECIES(SP_ID1)%VIBDOF)) EXIT
         END DO
         particles(JP1)%EVIB = EI
         ETR = ECOLL - EI

        ELSE IF (rfp .LT. PROT1 + PROT2 + PVIB1 + PVIB2) THEN
         ! Particle 2 selected for vibrational energy exchange
         ECOLL = ETR + particles(JP2)%EVIB
         DO
            R1 = rf()
            EI = R1*ECOLL
            R2 = rf()
            IF (R2 .LE. FUNCTION_I(ECOLL, EI, TRDOF, SPECIES(SP_ID2)%VIBDOF)) EXIT
         END DO
         particles(JP2)%EVIB = EI
         ETR = ECOLL - EI

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

  
  END MODULE collisions

