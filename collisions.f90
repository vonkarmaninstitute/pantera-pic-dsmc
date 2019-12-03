  ! Contains subroutines related to particle collisions

  MODULE collisions

  USE global
  USE mpi_common
  USE screen
  USE tools
  
  IMPLICIT NONE
  
  CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE DSMC_COLLISIONS -> Calcola collisioni !!!!!!!!!!!!!!!!!!!
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

   NPC = 0
   DO JP = 1, NP_PROC
      JC = particles(JP)%IC
      NPC(JC) = NPC(JC) + 1
   END DO

   NCELLS = NX*NY

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
 
   ! Calcola collisioni fra particelle
 
   TIMESTEP_COLL = 0
   DO JC = 1, NCELLS
      IF (NPC(JC) .GT. 1) THEN
         CALL COLLIS(JC, NPC, IOF, IND, NCOLLREAL)
         TIMESTEP_COLL = TIMESTEP_COLL + NCOLLREAL
      END IF
   END DO
 
   !WRITE(*,*) 'Number of real collisions: ', TIMESTEP_COLL

   DEALLOCATE(IND)
 
   END SUBROUTINE DSMC_COLLISIONS
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE COLLIS -> Calcola le collisioni tra le particelle !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COLLIS(JC,NPC,IOF,IND, NCOLLREAL)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: JC
  INTEGER, DIMENSION(:), INTENT(IN) :: NPC, IOF
  INTEGER, DIMENSION(:), INTENT(IN) :: IND
  INTEGER, INTENT(OUT) :: NCOLLREAL


  INTEGER      :: IOFJ,IOLJ,NPCJ,JP1,JP2,JCOL, JP, INDJ, NSP, JSP, I
  INTEGER      :: SP_ID1, SP_ID2
  INTEGER      :: NCOLL,NCOLLMAX_INT
  REAL(KIND=8) :: SIGMA, ALPHA, SIGMAMAX, NCOLLMAX,FCORR,VR,VR2,rfp
  REAL(KIND=8) :: OMEGA, CREF, ZETA, MASS
  REAL(KIND=8) :: PPMAX
  REAL(KIND=8) :: B,C,EINT,ETOT,ETR,PHI,SITETA,VRX,VRY,VRZ
  REAL(KIND=8) :: VXMAX,VXMIN,VYMAX,VYMIN,VZMAX,VZMIN,VRMAX
  REAL(KIND=8) :: PI,PI2
  REAL(KIND=8), DIMENSION(3) :: C1, C2, GREL, W
  REAL(KIND=8) :: GX, GY, GZ, G
  REAL(KIND=8) :: COSCHI, SINCHI, THETA, COSTHETA, SINTHETA
  REAL(KIND=8) :: M1, M2, MTOT, MR, COSA, SINA, BB
  
  PI   = 3.141593
  PI2  = 2.*PI

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

  NSP = MIXTURES(DSMC_COLL_MIX)%N_COMPONENTS
  SIGMAMAX = 0
  DO JSP = 1, NSP
    SIGMA = SPECIES( MIXTURES(DSMC_COLL_MIX)%COMPONENTS(JSP)%ID )%SIGMA
    IF (SIGMA.GT.SIGMAMAX) SIGMAMAX = SIGMA
  END DO

  VRMAX    = SQRT((VXMAX-VXMIN)**2 + (VYMAX-VYMIN)**2 + (VZMAX-VZMIN)**2)

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

  DO JCOL = 1,NCOLL
  
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
     CREF  = 0.5 * (SPECIES(SP_ID1)%CREF + SPECIES(SP_ID2)%CREF) ! Not correct DBDBDBDBDBDBDDBDBDB
     
     ALPHA = 0.5 * (SPECIES(SP_ID1)%ALPHA + SPECIES(SP_ID2)%ALPHA)
     M1    = SPECIES(SP_ID1)%MOLMASS
     M2    = SPECIES(SP_ID2)%MOLMASS
     MTOT  = M1 + M2
     MASS  = 0.5 * MTOT
     
     C1(1) = particles(JP1)%VX
     C1(2) = particles(JP1)%VY
     C1(3) = particles(JP1)%VZ

     C2(1) = particles(JP2)%VX
     C2(2) = particles(JP2)%VY
     C2(3) = particles(JP2)%VZ

     ! Check if collision happens
     rfp = rf()
     ZETA = 0.


     !WRITE(*,*) SIGMA*(VR/CREF)**(1.-2.*OMEGA)/5.3e-19
     !WRITE(*,*) 'Collision probability:', (FCORR/(SIGMAMAX*VRMAX)*SIGMA*(VR/CREF)**(1.-2.*NU))
     IF ( rfp .LT. (FCORR/(SIGMAMAX*VRMAX)*VR*SIGMA*(VR/CREF)**(1.-2.*OMEGA)) ) THEN

        ! Rimuovere commento per avere avviso
        IF ((FCORR/(SIGMAMAX*VRMAX)*VR*SIGMA*(VR/CREF)**(1.-2.*OMEGA)) .GT. 1.) THEN
           WRITE(*,*) 'Attention => this was a bad collision!!!'
        END IF

        NCOLLREAL = NCOLLREAL + 1

        ! Select elastic/inelastic collision 

        rfp = rf()

        IF ( rfp .LT. ZETA ) THEN
               
           ! Inelastic collision - Larsen/Borgnakke model

           ETOT = MASS*VR2/4. + particles(JP1)%EI + particles(JP2)%EI

           ! Step 1. Compute post-collision translational energy

10         rfp = rf()

           PPMAX = ((2.5-OMEGA)/(1.5-OMEGA)*rfp)**(1.5-OMEGA) * (2.5-OMEGA)*(1.-rfp)

           IF ( PPMAX .GE. rf() ) THEN
              ETR = ETOT * rfp
           ELSE 
              GO TO 10                 
           END IF

           ! Step 2. Compute post-collision internal energy

           EINT = ETOT - ETR
           
           particles(JP1)%EI = EINT*rf()
           particles(JP2)%EI = EINT - particles(JP1)%EI

           ! Step 3. Compute post-collision velocities

           VR = 2.*SQRT(ETR/MASS)

           B = 1. - 2.*rf()
           SITETA = SQRT(1.-B*B)

           PHI = PI2*rf()

           VRZ = VR*B
           VRX = VR*SITETA*COS(PHI)
           VRY = VR*SITETA*SIN(PHI)

           C = particles(JP1)%VX + particles(JP2)%VX
           particles(JP1)%VX = (C-VRX)*0.5
           particles(JP2)%VX = (C+VRX)*0.5

           C = particles(JP1)%VY + particles(JP2)%VY
           particles(JP1)%VY = (C-VRY)*0.5
           particles(JP2)%VY = (C+VRY)*0.5

           C = particles(JP1)%VZ + particles(JP2)%VZ
           particles(JP1)%VZ = (C-VRZ)*0.5
           particles(JP2)%VZ = (C+VRZ)*0.5

        ELSE

         ! Elastic collision => Compute post-collision velocities
         ! HS Collisions
         !   B = 1. - 2.*rf() ! COSTHETA
         !   SITETA = SQRT(1.-B*B) ! SINTHETA

         !   PHI = PI2*rf()

         !   VRZ = VR*B
         !   VRX = VR*SITETA*COS(PHI)
         !   VRY = VR*SITETA*SIN(PHI)

         !   C = particles(JP1)%VX + particles(JP2)%VX
         !   particles(JP1)%VX = (C-VRX)*0.5
         !   particles(JP2)%VX = (C+VRX)*0.5

         !   C = particles(JP1)%VY + particles(JP2)%VY
         !   particles(JP1)%VY = (C-VRY)*0.5
         !   particles(JP2)%VY = (C+VRY)*0.5

         !   C = particles(JP1)%VZ + particles(JP2)%VZ
         !   particles(JP1)%VZ = (C-VRZ)*0.5
         !   particles(JP2)%VZ = (C+VRZ)*0.5
         ! End of HS collisions

         ! VSS Collisions
         !ALPHA = 1.0 ! DBDBDBDBDBDBDBDBDBDBBDBDBDBBDBDBDBBDBDBBDBDBDBBDDBBDBDBDBDBBDBDBBDBD
         COSCHI = 2.*rf()**(1./ALPHA) - 1.
         SINCHI = SQRT(1.-COSCHI*COSCHI)
         THETA = PI2*rf()
         COSTHETA = COS(THETA)
         SINTHETA = SIN(THETA)

         MR = M1*M2/(M1+M2)         
         ETR = 0.5*MR*VR**2 ! Will be different for inelastic collisions.
         G = SQRT(2.*ETR/MR)

         ! Randomize relative velocity vector
         COSA = 2.*rf()-1.0
         SINA = SQRT(1-COSA*COSA)
         BB = PI2 * rf()


         GX = G*SINA*COS(BB)
         GY = G*SINA*SIN(BB)
         GZ = G*COSA


         !GX = G*SINCHI*COSTHETA
         !GY = G*SINCHI*SINTHETA
         !GZ = G*COSCHI
         
         !GX = particles(JP1)%VX - particles(JP2)%VX
         !GY = particles(JP1)%VY - particles(JP2)%VY
         !GZ = particles(JP1)%VZ - particles(JP2)%VZ

         GREL(1) = GX*COSCHI + SQRT(GY*GY+GZ*GZ)*SINTHETA*SINCHI
         GREL(2) = GY*COSCHI + (G*GZ*COSTHETA - GX*GY*SINTHETA)/SQRT(GY*GY+GZ*GZ)*SINCHI
         GREL(3) = GZ*COSCHI - (G*GY*COSTHETA + GX*GZ*SINTHETA)/SQRT(GY*GY+GZ*GZ)*SINCHI

         !GREL(1) = GX
         !GREL(2) = GY
         !GREL(3) = GZ

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

     END IF

  END DO  

  !WRITE(*,*) 'Actually performed:', NCOLLREAL
  !WRITE(*,*) NCOLL/(DT*NPC(JC))/MCRVHS, NCOLLREAL/(DT*NPC(JC))/MCRVHS 

  RETURN
      
  END SUBROUTINE COLLIS
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! SUBROUTINE BGK_COLLISIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Implements BGK model for the collision operator.           !!
  ! Performs relaxation towards a local Maxwellian.            !!
  !                                                            !!
  ! ATTENTION! WORKS ONLY FOR SINGLE SPECIES GAS MIXTURE!!!    !!
  ! The temperature of the local Maxwellian is found from the  !!
  ! internal energy (translational PLUS rotational). The local !!
  ! Maxwellian has a rotational temperature equal to the       !!
  ! translational one.                                         !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  SUBROUTINE BGK_COLLISIONS

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

!!!!!!! 
!!!!!!! ! index of the particle in the "particle" array, indexed by particles I own and sorted by cells
!!!!!!!   ! NPC => number of particles in each cell, indexed by cell index (a block length for IND) 
!!!!!!! 
!!!!!!!    ! Reorder particles
!!!!!!!    ! IND => index of particle in particle array, indexed by particles I own and sorted by cells
!!!!!!!    ! NPC => number of particles in cell, indexed by cell index (a block length for IND) 
!!!!!!!    ! could this be indexed by all cells? for those that I don't own, NPC will turn out to be 0.
!!!!!!!    ! + would simplify significantly,
!!!!!!!    ! - could be up to about 50 times bigger
!!!!!!!    ! + it's just a list of integers
!!!!!!!    ! --->  I am gonna try and in case revert back.
!!!!!!!    ! IOF => index (for IND) of first particle in a cell (an offset in IND)
   
   INTEGER, DIMENSION(NX*NY)          :: NPC, IOF
   INTEGER, DIMENSION(:), ALLOCATABLE :: IND
   INTEGER                            :: JP, JC, NCELLS, IDX
   INTEGER                            :: NCOLLREAL

   ! For actual collisions
   INTEGER :: IDp 
   REAL(KIND=8) :: VX_CELL, VY_CELL, VZ_CELL, P_CELL, T_CELL, n_CELL, CX, CY, CZ, C2, MASS
   REAL(KIND=8) :: nu, V_TH
   REAL(KIND=8) :: VX_NEW, VY_NEW, VZ_NEW, EI_NEW


   INTEGER :: JP_START, JP_END
   REAL(KIND=8) :: kB = 1.38064852E-23
   REAL(KIND=8) :: pi = 2.0*ASIN(1.0d0)

   ! Check that the gas is composed by only one species
   IF (SIZE(SPECIES) .NE. 1) THEN
      PRINT*
      PRINT*, "  ATTENTION! BGK collisions work for single species gas only. Sorry about that."
      PRINT*, "  You may have one only species in your actual particle vector, but I see the SPECIES "
      PRINT*, "  structure has more entries, so I will stop here just to be sure. Use one only species!"
      PRINT*
      PRINT*, "  ABORTING!"
      PRINT*
      STOP
   END IF

   PRINT*, "DBDBDB ", SIZE(SPECIES)

   NPC = 0
   DO JP = 1, NP_PROC
      JC = particles(JP)%IC
      NPC(JC) = NPC(JC) + 1
   END DO

   NCELLS = NX*NY

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
 
   ! Calcola collisioni fra particelle
 
   TIMESTEP_COLL = 0 ! Init number of collisions that happened

   DO JC = 1, NCELLS

      JP_START = IOF(JC)                ! First particle in the cell 
      JP_END   = IOF(JC) + NPC(JC) - 1  ! Last particle in the cell

      ! ++++++++ Compute average quantities in cell +++++++++++
      ! First, put internal quantities to zero
      VX_CELL = 0
      VY_CELL = 0
      VZ_CELL = 0

      V_TH   = 0
      P_CELL = 0

      ! Compute average velocities
      DO JP = JP_START, JP_END ! Loop on particles in cell
        
         IDp = IND(JP) ! ID of current particle
         VX_CELL = VX_CELL + particles(IDp)%VX/NPC(JC)
         VY_CELL = VY_CELL + particles(IDp)%VY/NPC(JC)
         VZ_CELL = VZ_CELL + particles(IDp)%VZ/NPC(JC)

      END DO

      ! Pressure, temperature (and internal temperature?)
      DO JP = JP_START, JP_END ! Loop on particles in cell

         IDp = IND(JP) ! ID of current particle
         CX  = particles(IDp)%VX - VX_CELL
         CY  = particles(IDp)%VY - VY_CELL
         CZ  = particles(IDp)%VZ - VZ_CELL

         C2  = CX**2 + CY**2 + CZ**2

         V_TH = V_TH + SQRT(C2)/NPC(JC) ! Compute thermal velocity (based on translational energy)

         MASS = SPECIES(particles(IDp)%S_ID)%MOLMASS 

         P_CELL = P_CELL + FNUM*MASS*(C2)/3/CELL_VOL ! Pressure in the cell

         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"
         PRINT*, "THIS IS A WORK IN PROGRESS!!!!!! INTRODUCE ROTATIONAL ENERGY!!!!!!"

         STOP

      END DO

      n_CELL = FNUM*NPC(JC)/CELL_VOL
      T_CELL = P_CELL/(n_CELL*kB)

      nu = n_CELL*BGK_SIGMA*V_TH ! Compute collision frequency (CONSTANT SIGMA FOR NOW)

      ! ++++++++++ Test particles for collisions +++++++++++
      DO JP = JP_START, JP_END

         IDp = IND(JP) ! ID of current particle

         ! Try collision - note that in this BGK model I use a collision frequency function of T:
         !                 collisions may or may not happen according to \nu, independently to 
         !                 the actual particle velocity.

         IF (rf() .LT. (1 - EXP(-nu*DT))) THEN ! Pcoll = 1 - exp(-nu*dt)

            ! Sample from Maxwellian (isotropic) at local velocity and temperature 
            CALL MAXWELL(VX_CELL, VY_CELL, VZ_CELL, T_CELL, T_CELL, T_CELL, T_CELL, &
                         VX_NEW, VY_NEW, VZ_NEW, EI_NEW, MASS)

            ! Assign new velocities and internal energy to particle
            particles(IDp)%VX = VX_NEW
            particles(IDp)%VY = VY_NEW
            particles(IDp)%VZ = VZ_NEW

            particles(IDp)%EI = EI_NEW

            ! Update number of collisions happened
            TIMESTEP_COLL = TIMESTEP_COLL + 1

         END IF
      END DO

   END DO

   DEALLOCATE(IND)

  END SUBROUTINE BGK_COLLISIONS

  
 
  END MODULE collisions

