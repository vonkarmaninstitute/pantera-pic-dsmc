! This module has the routines for initializing the program

MODULE initialization

   USE global
   USE screen
   USE mpi_common
   USE tools
   USE grid_and_partition

   IMPLICIT NONE

   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBORUTINE READINPUT -> reads input file and initializes variables !!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE READINPUT

      IMPLICIT NONE

      INTEGER, PARAMETER :: in1 = 444
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      ! Open input file for reading
      OPEN(UNIT=in1,FILE='input', STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         PRINT*
         WRITE(*,*)'  Attention, "input" file not found! ABORTING.'
         PRINT*
         STOP
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in1,'(A)', IOSTAT=ReasonEOF) line ! Read line
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         IF (line=='Axial_symmetry_bool:')     READ(in1,'(L6)') BOOL_AXI
         IF (line=='Domain_limits:')           READ(in1,*) XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
         IF (line=='Domain_periodicity:')      READ(in1,*) BOOL_X_PERIODIC, BOOL_Y_PERIODIC, BOOL_Z_PERIODIC 
         IF (line=='Number_of_cells:')         READ(in1,*) NX, NY, NZ

         ! ~~~~~~~~~~~~~  Numerical settings  ~~~~~~~~~~~~~~~~~
         IF (line=='Fnum:')                    READ(in1,*) FNUM
         IF (line=='Timestep:')                READ(in1,*) DT
         IF (line=='Number_of_timesteps:')     READ(in1,*) NT
         IF (line=='RNG_seed:')                READ(in1,*) RNG_SEED_GLOBAL

         ! ~~~~~~~~~~~~~  Collision type  ~~~~~~~~~~~~~~~~~
         IF (line=='Collision_type:')          READ(in1,*) COLLISION_TYPE
         IF (line=='MCC_background_dens:')     READ(in1,*) MCC_BG_DENS
         IF (line=='MCC_cross_section:')       READ(in1,*) MCC_SIGMA

         ! ~~~~~~~~~~~~~  Initial particles seeding  ~~~~~~~~~~~~~~~~~
         IF (line=='Initial_particles_bool:')  READ(in1,*) BOOL_INITIAL_SEED
         IF (line=='Initial_particles_dens:')  READ(in1,*) NRHO_INIT
         IF (line=='Initial_particles_vel:')   READ(in1,*) UX_INIT, UY_INIT, UZ_INIT
         IF (line=='Initial_particles_Ttra:')  READ(in1,*) TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT
         IF (line=='Initial_particles_Trot:')  READ(in1,*) TROT_INIT

         ! ~~~~~~~~~~~~~  Particle injection at boundaries ~~~~~~~~~~~
         IF (line=='Boundaries_inject_bool:')      READ(in1,*)BOOL_BOUNDINJECT
         IF (line=='Boundaries_inject_which_bool:')READ(in1,*)BOOL_INJ_XMIN,BOOL_INJ_XMAX,BOOL_INJ_YMIN, BOOL_INJ_YMAX
         IF (line=='Boundaries_inject_dens:')      READ(in1,*)NRHO_BOUNDINJECT
         IF (line=='Boundaries_inject_vel:')       READ(in1,*)UX_BOUND, UY_BOUND, UZ_BOUND
         IF (line=='Boundaries_inject_Ttra:')      READ(in1,*)TTRA_BOUND
         IF (line=='Boundaries_inject_Trot:')      READ(in1,*)TROT_BOUND

         ! ~~~~~~~~~~~~~  Particle injection at line source ~~~~~~~~~~~
         IF (line=='Linesource_bool:')      READ(in1,*)BOOL_LINESOURCE
         IF (line=='Linesource_center:')    READ(in1,*)X_LINESOURCE, Y_LINESOURCE
         IF (line=='Linesource_length_y:')  READ(in1,*)L_LINESOURCE
         IF (line=='Linesource_dens:')      READ(in1,*)NRHO_LINESOURCE
         IF (line=='Linesource_vel:')       READ(in1,*)UX_LINESOURCE, UY_LINESOURCE, UZ_LINESOURCE
         IF (line=='Linesource_Ttra:')      READ(in1,*)TTRA_LINESOURCE
         IF (line=='Linesource_Trot:')      READ(in1,*)TROT_LINESOURCE 

         ! ~~~~~~~~~~~~~  MPI parallelization settings ~~~~~~~~~~~~~~~
         IF (line=='Partition_style:')         READ(in1,*) DOMPART_TYPE
         IF (line=='Partition_num_blocks:')    READ(in1,*) N_BLOCKS_X, N_BLOCKS_Y

      END DO ! Loop for reading input file

      CLOSE(in1) ! Close input file

      CALL INPUT_DATA_SANITY_CHECK ! Check the values that were read
      CALL PRINTINPUT              ! Print values

   END SUBROUTINE READINPUT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE PRINTINPUT -> Asks master to print the data foud in the input file !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE PRINTINPUT

      IMPLICIT NONE
   
      CHARACTER(LEN=512) string

      IF (PROC_ID == 0) THEN ! Only master prints

         ! ~~~~ Geometry and domain ~~~~   
         WRITE(*,*) '  =========== Geometry and domain =================================='
         string = 'Axisymmetric geometry bool [T/F]: '
         WRITE(*,'(A5, A50,L)') '     ', string, BOOL_AXI
 
         string = 'Domain limits:'
         WRITE(*,'(A5,A50)') '    ', string
         WRITE(*,*)          '      X [m]:', XMIN, XMAX
         WRITE(*,*)          '      Y [m]:', YMIN, YMAX
         WRITE(*,*)          '      Z [m]:', ZMIN, ZMAX

         string = 'Domain periodicity along X, Y, Z [T/F]: '
         WRITE(*,'(A5,A50,L4,L4,L4)')'    ',string,BOOL_X_PERIODIC,BOOL_Y_PERIODIC,BOOL_Z_PERIODIC

         string = 'Number of cells:'
         WRITE(*,'(A5,A50, I8, I8, I8)') '    ', string, NX, NY, NZ

  
         ! ~~~~ Numerical settings ~~~~
         WRITE(*,*) '  =========== Numerical settings ==================================='
         string = 'Fnum:'
         WRITE(*,'(A5,A50,ES14.3)') '     ', string, FNUM

         string = 'Timestep [s]:'
         WRITE(*,'(A5,A50,ES14.3)') '     ', string, DT
 
         string = 'Number of timesteps:'
         WRITE(*,'(A5,A50,I9)') '     ',     string, NT

         string = 'RNG seed (global):'
         WRITE(*,'(A5,A50,I9)') '     ', string, RNG_SEED_GLOBAL

         ! ~~~~ Collisions ~~~~
         WRITE(*,*) '  =========== Collisions ============================'
         string = 'Collision type:'
         WRITE(*,'(A5,A50,A64)') '     ', string, COLLISION_TYPE 

         IF (COLLISION_TYPE == 'MCC') THEN  ! Only print this if MCC collisions are ON

            string = 'Background number density [1/m^3]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, MCC_BG_DENS

            string = 'Collision cross-section [m^2]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, MCC_SIGMA 

         END IF

         ! ~~~~ Initial particles seeding ~~~~
         WRITE(*,*) '  =========== Initial particles seeding ============================'
         string = 'Initial particle seeding bool [T/F]:'
         WRITE(*,'(A5,A50,L)') '     ', string, BOOL_INITIAL_SEED

         IF (BOOL_INITIAL_SEED) THEN ! Only print if this is ON

            string = 'Initial particles number density [1/m^3]:'
            WRITE(*,'(A5,A50,ES14.3)') '    ', string, NRHO_INIT
   
            string = 'Initial particles average velocities [m/s]:'
            WRITE(*,'(A5,A50,ES14.3,ES14.3,ES14.3)') '    ', string, UX_INIT, UY_INIT, UZ_INIT
   
            string = 'Initial particles transl temperatures [K]:'
            WRITE(*,'(A5,A50,ES14.3,ES14.3,ES14.3)') '    ', string, TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT
   
            string = 'Initial particles rotational temp [K]:'
            WRITE(*,'(A5,A50,ES14.3)') '    ', string, TROT_INIT

         END IF

         ! ~~~~ Injection at boundaries ~~~~
         WRITE(*,*) '  =========== Particle injection at boundaries ====================='

         string = 'Particle injection at boundaries bool [T/F]:'
         WRITE(*,'(A5,A50,L)') '     ', string, BOOL_BOUNDINJECT

         IF (BOOL_BOUNDINJECT) THEN ! Only print if this is ON 

            string = 'Injection at XMIN XMAX YMIN YMAX [T/F]:'
            WRITE(*,'(A5,A50,L,L,L,L)') '     ',string,BOOL_INJ_XMIN,BOOL_INJ_XMAX,BOOL_INJ_YMIN,BOOL_INJ_YMAX
   
            string = 'Boundary injection number density [1/m^3]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, NRHO_BOUNDINJECT
   
            string = 'Boundary injection average velocities [m/s]:'
            WRITE(*,'(A5,A50,ES14.3,ES14.3,ES14.3)') '     ', string, UX_BOUND, UY_BOUND, UZ_BOUND
   
            string = 'Boundary injection transl and rot Temps [K]:'
            WRITE(*,'(A5,A50,ES14.3,ES14.3)') '     ', string, TTRA_BOUND, TROT_BOUND
   
         END IF

         ! ~~~~ Injection from line source ~~~~
         WRITE(*,*) '  =========== Particle injection from line source ====================='

         string = 'Particle injection from line source bool [T/F]:'
         WRITE(*,'(A5,A50,L)') '     ', string, BOOL_LINESOURCE

         IF (BOOL_LINESOURCE) THEN ! Only print if this is ON
 
            string = 'Line source center position [m]:'
            WRITE(*,'(A5,A50,ES14.3,ES14.3)') '     ', string, X_LINESOURCE, Y_LINESOURCE
   
            string = 'Line source y-length [m]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, L_LINESOURCE
   
            string = 'Line source injection number density [1/m^3]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, NRHO_LINESOURCE
   
            string = 'Line source injection average velocities [m/s]:'
            WRITE(*,'(A5,A50,ES14.3,ES14.3,ES14.3)')'     ',string,UX_LINESOURCE,UY_LINESOURCE,UZ_LINESOURCE
   
            string = 'Line source transl and rot Temps [K]:'
            WRITE(*,'(A5,A50,ES14.3,ES14.3)') '     ', string, TTRA_LINESOURCE, TROT_LINESOURCE

         END IF

         ! ~~~~ MPI parallelization ~~~~
         WRITE(*,*) '  =========== MPI parallelization settings ========================='

         string = 'Number of MPI processes:'
         WRITE(*,'(A5,A50,I9)') '    ', string, N_MPI_THREADS     

         IF (DOMPART_TYPE .EQ. 0) THEN ! Strips
             string = 'Partition style: strips'
             WRITE(*,'(A5,A50)') '    ', string     

         ELSE IF (DOMPART_TYPE .EQ. 1) THEN ! Blocks
             string = 'Partition style: blocks'
             WRITE(*,'(A5,A50)') '    ', string

             string = 'Number of partitions along X and Y:'
             WRITE(*,'(A5,A50,I9,I9)') '    ', string, N_BLOCKS_X, N_BLOCKS_Y
         END IF



      END IF

   END SUBROUTINE PRINTINPUT


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INITIAL_SEED -> seeds particles at the beginning of !!
   ! the simulation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE INITIAL_SEED
      
      ! This subroutines seeds particles at the beginning of the simulation.
      ! 
      ! Every processor creates some particles in ALL the domain. These particles are then
      ! sent to the respective processors by a call to the "PARTICLE_EXCHANGE" subroutine.
      ! (Alternatively we could ask to one only process to make the particles and send them to the
      ! proper processes, OR ask to the processes to make particles in their own boundaries only,
      ! but you know, this doesn't matter that much. And, this is similar to what is done in the 
      ! time loop, where we use the MPI_ALLTOALL subroutine).

      IMPLICIT NONE

      INTEGER            :: IP, NP_INIT, NPPP_INIT
      REAL(KIND=8)       :: Xp, Yp, Zp, VXp, VYp, VZp, EIp
      INTEGER            :: CID

      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
      REAL(KIND=8)  :: RGAS ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDB
      INTEGER       :: idummy

      ! Print message 
      CALL ONLYMASTERPRINT1(PROC_ID, '> SEEDING INITIAL PARTICLES IN THE DOMAIN...')
     
      ! Compute number of particles to be seeded
      NP_INIT = NINT(NRHO_INIT/FNUM*(XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN))

      CALL ONLYMASTERPRINT2(PROC_ID, '   Particles to be seeded', REAL(NP_INIT, KIND=8))
 
      ! Compute number of particles to be seeded by every process
      NPPP_INIT = INT(NP_INIT/N_MPI_THREADS) 

      ! Create particles in the domain
      DO IP = 1, NPPP_INIT

         ! Create particle position randomly in the domain
         XP = XMIN + (XMAX-XMIN)*rand()
         YP = YMIN + (YMAX-YMIN)*rand()
         ZP = ZMIN + (ZMAX-ZMIN)*rand()

         ! Assign velocity and energy following a Boltzmann distribution
         RGAS = 200. ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDB
         CALL MAXWELL(UX_INIT, UY_INIT, UZ_INIT, &
                      TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT, TROT_INIT, &
                      VXP, VYP, VZP, EIP, RGAS)

         CALL CELL_FROM_POSITION(XP,YP,  CID) ! Find cell containing particle

         idummy = 0
         CALL INIT_PARTICLE(XP,YP,ZP,VXP,VYP,VZP,EIP,idummy,CID,  particleNOW) ! Save in particle
         CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
      END DO

      ! ~~~~~~ At this point, exchange particles among processes ~~~~~~
      CALL EXCHANGE

   END SUBROUTINE INITIAL_SEED

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INPUT_DATA_SANITY_CHECK -> checks validity of some input data !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE INPUT_DATA_SANITY_CHECK

      IMPLICIT NONE
 
      ! ------------ Check values for number of cells ------------
      IF ( (NX < 1) .OR. (NY < 1) ) THEN
         WRITE(*,*)
         WRITE(*,*) 'ERROR! Number of cells along X or Y smaller than one. ABORTING!'
         WRITE(*,*)
         STOP
      END IF 

      IF (NZ /= 1) THEN
         WRITE(*,*)
         WRITE(*,*) 'ERROR! Number of cells along Z different than 1 is not supported. ABORTING!'
         WRITE(*,*)
         STOP
      END IF 
     
      ! ------------ Check geometrical symmetries and flags ------
      ! For axisymmetric simulations, Y cannot be periodic (YMIN is the symmetry axis)
      IF (BOOL_AXI) THEN 
         IF (BOOL_Y_PERIODIC) THEN
            WRITE(*,*)
            WRITE(*,*) 'ERROR! For axisymmetric simulations, Y cannot be periodic:'
            WRITE(*,*) 'YMIN is the symmetry axis! Check the input file. ABORTING!'
            WRITE(*,*)
            STOP
         END IF
      END IF

      ! ------------ MPI settings --------------------------------
      IF (DOMPART_TYPE == -1) THEN ! It's the initialization value. It was not recognized.

         IF (N_MPI_THREADS == 1) THEN ! Serial running, no problem.
            ! Then chill, it's fine, don't do anything.

         ELSE 
            WRITE(*,*)
            WRITE(*,*) '****************************************************************************'
            WRITE(*,*) '* ERROR! No partitioning strategy specified in the input file.              '
            WRITE(*,*) '* Use keyword "Partition_style:" ABORTING!'
            WRITE(*,*) '****************************************************************************'
            WRITE(*,*)
            STOP

         END IF

      END IF

      IF (DOMPART_TYPE == 1) THEN ! "blocks" domain partition

         IF (N_MPI_THREADS /= N_BLOCKS_X*N_BLOCKS_Y) THEN
            WRITE(*,*)
            WRITE(*,*) '****************************************************************************'
            WRITE(*,*) '* ERROR! Number of MPI threads is different than number of blocks'
            WRITE(*,*) '* specified in the input file! ABORTING!'
            WRITE(*,*) '****************************************************************************'
            WRITE(*,*)
            STOP
         END IF

      END IF


   END SUBROUTINE INPUT_DATA_SANITY_CHECK


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INITVARIOUS -> initializes various other variables !!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE INITVARIOUS 

      IMPLICIT NONE

      INTEGER         :: I
      REAL(KIND=8)    :: dummy
 
      ! ~~~~~~~~ Initial allocation for vector of particles in each process ~~~~~~~~
      NP_PROC = 0
      ALLOCATE(particles(0))

      ! ~~~~~~~~ Initialize RNG (random number generator) ~~~~~~~~~~

      ! Create local seed for RNG (different for each processor)
      RNG_SEED_LOCAL = RNG_SEED_GLOBAL*(PROC_ID + 1) 

      CALL SRAND(RNG_SEED_LOCAL) ! Seed RNG with local quantity
      ! And call it a number of times, so that initial correlation between the similar seeds 
      ! is lost (apparently, some nonnegligible correlation is retained!!!!!)
      DO I = 1,100
        dummy = RAND()
      END DO

      ! ~~~~~~~~ Additional variables for MPI domain decomposition ~~~~~~
      IF (DOMPART_TYPE == 1) THEN ! "blocks" domain partition
         ! Compute DX_BLOCKS, DY_BLOCKS
         DX_BLOCKS = (XMAX - XMIN)/N_BLOCKS_X
         DY_BLOCKS = (YMAX - YMIN)/N_BLOCKS_Y
      END IF

   END SUBROUTINE INITVARIOUS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INITINJECTION -> initializes variables for particles injection !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE INITINJECTION

      ! Initializes the variables needed for particle injection, either from 
      ! domain boundaries or from a line source.
      ! 
      ! ~~~~ Boundary injection: ~~~~
      !  
      ! 1) Compute molecular speed ratio "S" (kind of Mach number, it's the standard 
      !    deviation of the Maxwellian) with velocity perpendicular to considered surf;
      !
      ! 2) Compute flux entering through the wall from a Maxwellian (formulas can be found
      !    in Bird's "Molecular Gas Dynamics and Direct Simul..." (1994), pag. 82, Eq. 4.22,
      !    section "Fluxal properties");
      ! 
      ! 2.1) If the flux is exiting the surface and is hypersonic (U > 3 S means roughly Mach 6)
      !      then print a warning, as almost all particles emitted will be out of the domain
      !      and computing them will be worthless
      ! 
      ! 3) Compute how many particles "nfs_..." shall be injected by every processor
      ! 
      ! 4) Compute parameter "KAPPA", later used by the FLX function
      ! 
      ! ~~~~ Line source injection: ~~~~
      !
      ! Exactly like for the boundaries.

      IMPLICIT NONE

      REAL(KIND=8) :: BETA, FLUXBOUND, FLUXLINESOURCE, NtotINJECT, ACCA, Snow
      REAL(KIND=8) :: PI, PI2  
 
      REAL(KIND=8)  :: RGAS = 200. ! DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDB
    
      PI   = 3.141593
      PI2  = 2.*PI

      ! Injection from boundaries
      IF (BOOL_BOUNDINJECT) THEN 

         BETA = 1./SQRT(2.*RGAS*TTRA_BOUND) ! sqrt(M/(2*kB*T)), it's the Maxwellian std dev

         ! +++++++++++ Lower x ++++++++++++
         IF (BOOL_INJ_XMIN) THEN

            S_NORM_XMIN = UX_BOUND*BETA ! Molecular speed ratio normal to boundary
            Snow        = S_NORM_XMIN   ! temp variable

            FLUXBOUND = NRHO_BOUNDINJECT/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) & 
                                      + SQRT(PI)*Snow*(1.+ERF1(Snow)))    ! Tot number flux emitted [1/(s m^2)]
            NtotINJECT = FLUXBOUND*(YMAX-YMIN)*(ZMAX-ZMIN)*DT/FNUM        ! Tot num of particles to be injected
            nfs_XMIN = FLOOR(NtotINJECT/REAL(N_MPI_THREADS,KIND=8)+rf())  ! Particles injected by each proc
            ACCA = SQRT(Snow**2+2.)                                       ! Tmp variable
            KAPPA_XMIN = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA)) ! global variable

            ! Check: if we are hypersonic and stuff exits domain print a warning
            IF (UX_BOUND+3./BETA .LE. 0.) THEN 
               CALL ONLYMASTERPRINT1(PROC_ID, '$$$ Warning! Hypersonic boundary! Almost no &
                                                   &particles will be emitted at XMIN.')
            END IF

         END IF 

         ! +++++++++++ Higher x +++++++++++
         IF (BOOL_INJ_XMAX) THEN

            S_NORM_XMAX = - UX_BOUND*BETA ! Molecular speed ratio normal to boundary (inside)
            Snow        = S_NORM_XMAX     ! temp variable
            FLUXBOUND   = NRHO_BOUNDINJECT/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) &
                                      + SQRT(PI)*Snow*(1.+ERF1(Snow)))      ! Tot number flux emitted 
            NtotINJECT  = FLUXBOUND*(YMAX-YMIN)*(ZMAX-ZMIN)*DT/FNUM         ! Tot num of particles to be injected
            nfs_XMAX    = FLOOR(NtotINJECT/REAL(N_MPI_THREADS,KIND=8)+rf()) ! Particles injected by each proc
            ACCA        = SQRT(Snow**2+2.)                                  ! Tmp variable
            KAPPA_XMAX  = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA))  ! global variable

            ! Check: if we are hypersonic and stuff exits domain print a warning
            IF (UX_BOUND-3./BETA .GE. 0.) THEN 
               CALL ONLYMASTERPRINT1(PROC_ID, '$$$ Warning! Hypersonic boundary! Almost no &
                                                   &particles will be emitted at XMAX.')
            END IF

         END IF 

         ! +++++++++++ Lower y ++++++++++++
         IF (BOOL_INJ_YMIN) THEN

            S_NORM_YMIN = UY_BOUND*BETA ! Molecular speed ratio normal to boundary
            Snow        = S_NORM_YMIN   ! temp variable
            FLUXBOUND   = NRHO_BOUNDINJECT/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) &
                                      + SQRT(PI)*Snow*(1.+ERF1(Snow)))      ! Tot number flux emitted
            NtotINJECT  = FLUXBOUND*(XMAX-XMIN)*(ZMAX-ZMIN)*DT/FNUM         ! Tot num of particles to be injected
            nfs_YMIN    = FLOOR(NtotINJECT/REAL(N_MPI_THREADS,KIND=8)+rf()) ! Particles injected by each proc
            ACCA        = SQRT(Snow**2+2.)                                  ! Tmp variable
            KAPPA_YMIN  = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA))  ! global variable

            ! Check: if we are hypersonic and stuff exits domain print a warning
            IF (UY_BOUND+3./BETA .LE. 0.) THEN 
               CALL ONLYMASTERPRINT1(PROC_ID, '$$$ Warning! Hypersonic boundary! Almost no &
                                                   &particles will be emitted at YMIN.')
            END IF

         END IF 

         ! +++++++++++ Higher y +++++++++++
         IF (BOOL_INJ_YMAX) THEN

            S_NORM_YMAX = - UY_BOUND*BETA ! Molecular speed ratio normal to boundary (inside)
            Snow        = S_NORM_YMAX     ! temp variable
            FLUXBOUND   = NRHO_BOUNDINJECT/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) &
                                      + SQRT(PI)*Snow*(1.+ERF1(Snow)))      ! Tot number flux emitted 
            NtotINJECT  = FLUXBOUND*(XMAX-XMIN)*(ZMAX-ZMIN)*DT/FNUM         ! Tot num of particles to be injected
            nfs_YMAX    = FLOOR(NtotINJECT/REAL(N_MPI_THREADS,KIND=8)+rf()) ! Particles injected by each proc
            ACCA        = SQRT(Snow**2+2.)                                  ! Tmp variable
            KAPPA_YMAX  = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA))  ! global variable

            ! Check: if we are hypersonic and stuff exits domain print a warning
            IF (UY_BOUND-3./BETA .GE. 0.) THEN 
               CALL ONLYMASTERPRINT1(PROC_ID, '$$$ Warning! Hypersonic boundary! Almost no &
                                                   &particles will be emitted at YMAX.')
            END IF

         END IF

      END IF

      ! =====================================================
      ! Injection from line source
      IF (BOOL_LINESOURCE) THEN 

         BETA = 1./SQRT(2.*RGAS*TTRA_LINESOURCE) ! sqrt(M/(2*kB*T)), it's the Maxwellian std dev

         S_NORM_LINESOURCE = UX_LINESOURCE*BETA ! Molecular speed ratio normal to line source (which is y-aligned)
         Snow              = S_NORM_LINESOURCE  ! temp variable

         FLUXLINESOURCE = NRHO_LINESOURCE/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) & 
                                   + SQRT(PI)*Snow*(1.+ERF1(Snow)))          ! Tot number flux emitted [1/(s m^2)]
         NtotINJECT = FLUXLINESOURCE*L_LINESOURCE*(ZMAX-ZMIN)*DT/FNUM        ! Tot num of particles to be injected
         nfs_LINESOURCE = FLOOR(NtotINJECT/REAL(N_MPI_THREADS,KIND=8)+rf())  ! Particles injected by each proc
         ACCA = SQRT(Snow**2+2.)                                             ! Tmp variable
         KAPPA_LINESOURCE = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA)) ! global variable

      END IF 

   END SUBROUTINE INITINJECTION

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INITCOLLISIONS -> Initializes variables for the collisions !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE INITCOLLISIONS

   IMPLICIT NONE

   ! Set logicals
   IF (COLLISION_TYPE == "NONE") THEN
      BOOL_DSMC = .FALSE.
      BOOL_MCC  = .FALSE.

   ELSE IF (COLLISION_TYPE == "DSMC") THEN
      BOOL_DSMC = .TRUE.
      BOOL_MCC  = .FALSE.

   ELSE IF (COLLISION_TYPE == "MCC") THEN
      BOOL_DSMC = .FALSE.
      BOOL_MCC  = .TRUE.

   ELSE ! ERROR!
      CALL ONLYMASTERPRINT1(PROC_ID, '$$$ ATTENTION! Collision type in input file not recognized! ABORTING!')
      STOP

   END IF 

   END SUBROUTINE INITCOLLISIONS



END MODULE initialization
