! This module has the routines for initializing the program

MODULE initialization

   USE global
   USE screen
   USE mpi_common
   USE tools
   USE grid_and_partition
   USE mt19937_64

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

      CHARACTER*512      :: MIXTURE_DEFINITION, VSS_PARAMS_FILENAME, LINESOURCE_DEFINITION, WALL_DEFINITION
      CHARACTER*64       :: MIX_INIT_NAME, MIX_BOUNDINJECT_NAME, DSMC_COLL_MIX_NAME, MCC_BG_MIX_NAME

      ! Open input file for reading
      OPEN(UNIT=in1,FILE='input', STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, "input" file not found! ABORTING.')
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
         IF (line=='Dimensions:')              READ(in1,*) DIMS
         IF (line=='Domain_periodicity:') THEN
            READ(in1,*) BOOL_X_PERIODIC, BOOL_Y_PERIODIC, BOOL_Z_PERIODIC
            BOOL_PERIODIC(1) = BOOL_X_PERIODIC
            BOOL_PERIODIC(2) = BOOL_X_PERIODIC
            BOOL_PERIODIC(3) = BOOL_Y_PERIODIC
            BOOL_PERIODIC(4) = BOOL_Y_PERIODIC
         END IF
         IF (line=='Domain_specular:')         READ(in1,*) BOOL_SPECULAR
         IF (line=='Domain_diffuse:')          READ(in1,*) BOOL_DIFFUSE
         IF (line=='Boundary_temperature:')    READ(in1,*) BOUNDTEMP
         IF (line=='Domain_react:')            READ(in1,*) BOOL_REACT
         IF (line=='Wall_reactions_file:') THEN
            READ(in1,*) WALL_REACTIONS_FILENAME
            CALL READ_WALL_REACTIONS(WALL_REACTIONS_FILENAME)
         END IF
         IF (line=='Number_of_cells:')         READ(in1,*) NX, NY, NZ

         IF (line=='Quadtree_root_cells:') THEN
            READ(in1,*) NX, NY, NZ
            GRID_TYPE = QUADTREE
            CALL INIT_QUADTREE
         END IF

         IF (line=='Grid_file:') THEN
            READ(in1,*) GRID_FILENAME
            CALL READ_GRID_FILE(GRID_FILENAME)
            GRID_TYPE = RECTILINEAR_NONUNIFORM
         END IF
         ! ~~~~~~~~~~~~~  Numerical settings  ~~~~~~~~~~~~~~~~~
         IF (line=='Fnum:')                    READ(in1,*) FNUM
         IF (line=='Timestep:')                READ(in1,*) DT
         IF (line=='Number_of_timesteps:')     READ(in1,*) NT
         IF (line=='RNG_seed:')                READ(in1,*) RNG_SEED_GLOBAL
         IF (line=='Dump_part_every:')         READ(in1,*) DUMP_EVERY
         IF (line=='Dump_part_start:')         READ(in1,*) DUMP_START
         IF (line=='Dump_grid_avgevery:')      READ(in1,*) DUMP_GRID_AVG_EVERY
         IF (line=='Dump_grid_start:')         READ(in1,*) DUMP_GRID_START
         IF (line=='Dump_grid_numavgs:')       READ(in1,*) DUMP_GRID_N_AVG
         IF (line=='Perform_checks:')          READ(in1,*) PERFORM_CHECKS
         IF (line=='Checks_every:')            READ(in1,*) CHECKS_EVERY
         IF (line=='Stats_every:')             READ(in1,*) STATS_EVERY

         IF (line=='Bool_PIC:')                READ(in1,*) BOOL_PIC


         ! ~~~~~~~~~~~~~  Multispecies ~~~~~~~~~~~~~~~
         IF (line=='Species_file:') THEN
            READ(in1,*) SPECIES_FILENAME
            CALL READ_SPECIES
         END IF

         IF (line=='Def_mixture:') THEN
            READ(in1,'(A)') MIXTURE_DEFINITION
            CALL DEF_MIXTURE(MIXTURE_DEFINITION)
         END IF

         ! ~~~~~~~~~~~~~  Collision type  ~~~~~~~~~~~~~~~~~
         IF (line=='Collision_type:')          READ(in1,*) COLLISION_TYPE
         IF (line=='MCC_background_dens:')     READ(in1,*) MCC_BG_DENS
         IF (line=='MCC_background_Ttra:')     READ(in1,*) MCC_BG_TTRA
         IF (line=='MCC_background_mixture:') THEN
            READ(in1,*) MCC_BG_MIX_NAME
            MCC_BG_MIX = MIXTURE_NAME_TO_ID(MCC_BG_MIX_NAME)
         END IF
         IF (line=='VSS_parameters_file:')     THEN
            READ(in1,*) VSS_PARAMS_FILENAME
            CALL READ_VSS(VSS_PARAMS_FILENAME)
         END IF
         IF (line=='DSMC_collisions_mixture:') THEN
            READ(in1,*) DSMC_COLL_MIX_NAME
            DSMC_COLL_MIX = MIXTURE_NAME_TO_ID(DSMC_COLL_MIX_NAME)
         END IF
         
        ! ~~~~~~~~~~~~~  Reactions ~~~~~~~~~~~~~~~
         IF (line=='Reactions_file:') THEN
            READ(in1,*) REACTIONS_FILENAME
            CALL READ_REACTIONS(REACTIONS_FILENAME)
         END IF

         ! ~~~~~~~~~~~~~  Initial particles seeding  ~~~~~~~~~~~~~~~~~
         IF (line=='Initial_particles_bool:')  READ(in1,*) BOOL_INITIAL_SEED
         IF (line=='Initial_particles_dens:')  READ(in1,*) NRHO_INIT
         IF (line=='Initial_particles_vel:')   READ(in1,*) UX_INIT, UY_INIT, UZ_INIT
         IF (line=='Initial_particles_Ttra:')  READ(in1,*) TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT
         IF (line=='Initial_particles_Trot:')  READ(in1,*) TROT_INIT
         IF (line=='Initial_particles_Tvib:')  READ(in1,*) TVIB_INIT
         IF (line=='Initial_particles_mixture:') THEN
            READ(in1,*) MIX_INIT_NAME
            MIX_INIT = MIXTURE_NAME_TO_ID(MIX_INIT_NAME)
         END IF


         ! ~~~~~~~~~~~~~  Thermal bath  ~~~~~~~~~~~~~~~~~
         IF (line=='Thermal_bath_bool:')  READ(in1,*) BOOL_THERMAL_BATH
         IF (line=='Thermal_bath_Ttr:')  READ(in1,*) TBATH
         

         ! ~~~~~~~~~~~~~  Particle injection at boundaries ~~~~~~~~~~~
         IF (line=='Boundaries_inject_bool:')       READ(in1,*) BOOL_BOUNDINJECT
         IF (line=='Boundaries_inject_which_bool:') READ(in1,*) BOOL_INJ_XMIN,BOOL_INJ_XMAX,BOOL_INJ_YMIN, BOOL_INJ_YMAX
         IF (line=='Boundaries_inject_dens:')       READ(in1,*) NRHO_BOUNDINJECT
         IF (line=='Boundaries_inject_vel:')        READ(in1,*) UX_BOUND, UY_BOUND, UZ_BOUND
         IF (line=='Boundaries_inject_Ttra:')       READ(in1,*) TTRA_BOUND
         IF (line=='Boundaries_inject_Trot:')       READ(in1,*) TROT_BOUND
         IF (line=='Boundaries_inject_Tvib:')       READ(in1,*) TVIB_BOUND
         IF (line=='Boundaries_inject_mixture:') THEN
            READ(in1,*) MIX_BOUNDINJECT_NAME
            MIX_BOUNDINJECT = MIXTURE_NAME_TO_ID(MIX_BOUNDINJECT_NAME)
         END IF

         ! ~~~~~~~~~~~~~  Particle injection at line source ~~~~~~~~~~~
         IF (line=='Linesource_bool:')      READ(in1,*)BOOL_LINESOURCE
         IF (line=='Linesource_center:')    READ(in1,*)X_LINESOURCE, Y_LINESOURCE
         IF (line=='Linesource_length_y:')  READ(in1,*)L_LINESOURCE
         IF (line=='Linesource_dens:')      READ(in1,*)NRHO_LINESOURCE
         IF (line=='Linesource_vel:')       READ(in1,*)UX_LINESOURCE, UY_LINESOURCE, UZ_LINESOURCE
         IF (line=='Linesource_Ttra:')      READ(in1,*)TTRA_LINESOURCE
         IF (line=='Linesource_Trot:')      READ(in1,*)TROT_LINESOURCE 

         IF (line=='Linesource:') THEN
            READ(in1,'(A)') LINESOURCE_DEFINITION
            CALL DEF_LINESOURCE(LINESOURCE_DEFINITION)
         END IF


         ! ~~~~~~~~~~~~~  Walls ~~~~~~~~~~~~~~~

         IF (line=='Wall:') THEN
            READ(in1,'(A)') WALL_DEFINITION
            CALL DEF_WALL(WALL_DEFINITION)
         END IF

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
      INTEGER :: j, k
      CHARACTER(LEN=512) string

      IF (PROC_ID == 0) THEN ! Only master prints

         ! ~~~~ Geometry and domain ~~~~   
         WRITE(*,*) '  =========== Geometry and domain =================================='
         string = 'Dimensions of the simulation: '
         WRITE(*,'(A5, A50,I8)') '     ', string, DIMS
         
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

            string = 'Background translational temperature [K]::'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, MCC_BG_TTRA

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

            string = 'Initial particles vibrational temp [K]:'
            WRITE(*,'(A5,A50,ES14.3)') '    ', string, TVIB_INIT

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

         ! ~~~~ Multispecies ~~~~
         WRITE(*,*) '  =========== Species ========================='
         WRITE(*,*) '    ','File read'
         WRITE(*,'(A5,A30,I9)') '    ','Number of species present: ', N_SPECIES


         WRITE(*,*) '  =========== Mixtures ========================='
         DO j = 1, N_MIXTURES
            WRITE(*,*)'    ','Contents of mixture ', MIXTURES(j)%NAME, ' with ', MIXTURES(j)%N_COMPONENTS, ' components:'
            DO k = 1, MIXTURES(j)%N_COMPONENTS
               WRITE(*,*) '    ','Mixture component ', MIXTURES(j)%COMPONENTS(k)%NAME, &
                        ' with fraction ', MIXTURES(j)%COMPONENTS(k)%MOLFRAC, &
                        ' and ID ', MIXTURES(j)%COMPONENTS(k)%ID
            END DO
         END DO

         
      END IF

   END SUBROUTINE PRINTINPUT


   SUBROUTINE READ_SPECIES
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: in2 = 445
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER*10 :: NAME
      REAL(KIND=8) :: MOLWT
      REAL(KIND=8) :: MOLECULAR_MASS
      INTEGER      :: ROTDOF
      REAL(KIND=8) :: ROTREL
      INTEGER      :: VIBDOF
      REAL(KIND=8) :: VIBREL
      REAL(KIND=8) :: VIBTEMP
      REAL(KIND=8) :: SPWT
      REAL(KIND=8) :: CHARGE

   
      ! Open input file for reading
      OPEN(UNIT=in2,FILE=SPECIES_FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, species definition file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      N_SPECIES = 0
      DO

         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         
         !READ(line,'(A2, ES14.3, ES14.3, I1, ES14.3, I1, ES14.3, ES14.3, ES14.3, ES14.3)') &
         READ(line,*) &
         NAME, MOLWT, MOLECULAR_MASS, ROTDOF, ROTREL, VIBDOF, VIBREL, VIBTEMP, SPWT, CHARGE
      
         ALLOCATE(TEMP_SPECIES(N_SPECIES+1)) ! Append the species to the list
         TEMP_SPECIES(1:N_SPECIES) = SPECIES
         CALL MOVE_ALLOC(TEMP_SPECIES, SPECIES)
         N_SPECIES = N_SPECIES + 1

         SPECIES(N_SPECIES)%NAME = NAME
         SPECIES(N_SPECIES)%MOLWT = MOLWT
         SPECIES(N_SPECIES)%MOLECULAR_MASS = MOLECULAR_MASS
         SPECIES(N_SPECIES)%ROTDOF = ROTDOF
         SPECIES(N_SPECIES)%ROTREL = ROTREL
         SPECIES(N_SPECIES)%VIBDOF = VIBDOF
         SPECIES(N_SPECIES)%VIBREL = VIBREL
         SPECIES(N_SPECIES)%VIBTEMP = VIBTEMP
         SPECIES(N_SPECIES)%SPWT = SPWT
         SPECIES(N_SPECIES)%INVSPWT = 1./SPWT
         SPECIES(N_SPECIES)%CHARGE = CHARGE

         IF (SPWT .LT. 1.d0) THEN
            CALL ERROR_ABORT('Attention, Species specific particle weight must be >= 1! ABORTING.')
         ENDIF

      END DO

      CLOSE(in2) ! Close input file

   END SUBROUTINE READ_SPECIES

   SUBROUTINE READ_VSS(FILENAME)

      IMPLICIT NONE

      CHARACTER*64, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in2 = 456
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER*64 :: SP_NAME
      INTEGER      :: SP_ID
      REAL(KIND=8) :: DIAM
      REAL(KIND=8) :: OMEGA
      REAL(KIND=8) :: TREF
      REAL(KIND=8) :: ALPHA
      REAL(KIND=8) :: PI, KB

      INTEGER      :: IS, JS
      REAL(KIND=8) :: M1, M2, MRED

      PI   = 3.141593
      KB = 1.38064852E-23

      ! Open input file for reading
      OPEN(UNIT=in2,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, VSS parameters definition file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         
         !READ(line,'(A2, ES14.3, ES14.3, I1, ES14.3, I1, ES14.3, ES14.3, ES14.3, ES14.3)') &
         READ(line,*) SP_NAME, DIAM, OMEGA, TREF, ALPHA
      
         SP_ID = SPECIES_NAME_TO_ID(SP_NAME)

         SPECIES(SP_ID)%DIAM  = DIAM
         SPECIES(SP_ID)%OMEGA = OMEGA
         SPECIES(SP_ID)%TREF  = TREF
         SPECIES(SP_ID)%ALPHA = ALPHA

         SPECIES(SP_ID)%SIGMA = PI*DIAM**2
         SPECIES(SP_ID)%NU    = OMEGA - 0.5


         IF (ABS(OMEGA-0.5) .LT. 1.d-6) THEN
            SPECIES(SP_ID)%CREF  = (4.*KB*TREF/SPECIES(SP_ID)%MOLECULAR_MASS)**0.5 * 1.2354
         ELSE
            SPECIES(SP_ID)%CREF  = (4.*KB*TREF/SPECIES(SP_ID)%MOLECULAR_MASS)**0.5 * (GAMMA(2.5-OMEGA))**(-0.5/(OMEGA-0.5))
         END IF
         
         
         
         
      END DO

      ALLOCATE(GREFS(N_SPECIES, N_SPECIES))
      DO IS = 1, N_SPECIES
         DO JS = 1, N_SPECIES
            M1    = SPECIES(IS)%MOLECULAR_MASS
            M2    = SPECIES(JS)%MOLECULAR_MASS
            MRED  = M1*M2/(M1+M2)
            OMEGA    = 0.5 * (SPECIES(IS)%OMEGA + SPECIES(JS)%OMEGA)
            TREF = 0.5 * (SPECIES(IS)%TREF + SPECIES(JS)%TREF)
            
            IF (ABS(OMEGA-0.5) .LT. 1.d-6) THEN
               GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * 1.2354
            ELSE
               GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * (GAMMA(2.5-OMEGA))**(-0.5/(OMEGA-0.5))
            END IF
         END DO
      END DO

      CLOSE(in2) ! Close input file


   END SUBROUTINE READ_VSS


   SUBROUTINE DEF_MIXTURE(DEFINITION)
      
      IMPLICIT NONE

      CHARACTER*10 :: MIX_NAME
      CHARACTER*10 :: COMP_NAME
      REAL(KIND=8) :: MOLFRAC

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION
      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      INTEGER :: i, N_COMP
      !TYPE(MIXTURE_COMPONENT) NEW_MIXTURE
      TYPE(MIXTURE), DIMENSION(:), ALLOCATABLE :: TEMP_MIXTURES
      TYPE(MIXTURE_COMPONENT), DIMENSION(:), ALLOCATABLE :: TEMP_COMPONENTS


      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)
      
      READ(STRARRAY(1),'(A10)') MIX_NAME
      N_COMP = (N_STR-1)/2
      ALLOCATE(TEMP_COMPONENTS(N_COMP))
      DO i = 1, N_COMP
         READ(STRARRAY(2*i),'(A10)') COMP_NAME
         READ(STRARRAY(2*i+1),'(ES14.3)') MOLFRAC

         TEMP_COMPONENTS(i)%NAME = COMP_NAME
         TEMP_COMPONENTS(i)%MOLFRAC = MOLFRAC
         TEMP_COMPONENTS(i)%ID = SPECIES_NAME_TO_ID(COMP_NAME)
      END DO


      ALLOCATE(TEMP_MIXTURES(N_MIXTURES+1)) ! Append the mixture to the list
      TEMP_MIXTURES(1:N_MIXTURES) = MIXTURES
      CALL MOVE_ALLOC(TEMP_MIXTURES, MIXTURES)
      N_MIXTURES = N_MIXTURES + 1

      MIXTURES(N_MIXTURES)%NAME = MIX_NAME
      MIXTURES(N_MIXTURES)%N_COMPONENTS = N_COMP

      CALL MOVE_ALLOC(TEMP_COMPONENTS, MIXTURES(N_MIXTURES)%COMPONENTS)


   END SUBROUTINE DEF_MIXTURE


   
   SUBROUTINE DEF_LINESOURCE(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      CHARACTER*64 :: MIX_NAME
      INTEGER      :: MIX_ID
      TYPE(LINESOURCE), DIMENSION(:), ALLOCATABLE :: TEMP_LINESOURCES

      
      ALLOCATE(TEMP_LINESOURCES(N_LINESOURCES+1)) ! Append the mixture to the list
      TEMP_LINESOURCES(1:N_LINESOURCES) = LINESOURCES
      CALL MOVE_ALLOC(TEMP_LINESOURCES, LINESOURCES)
      N_LINESOURCES = N_LINESOURCES + 1



      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      READ(STRARRAY(1), '(ES14.3)') LINESOURCES(N_LINESOURCES)%CX
      READ(STRARRAY(2), '(ES14.3)') LINESOURCES(N_LINESOURCES)%CY
      READ(STRARRAY(3), '(ES14.3)') LINESOURCES(N_LINESOURCES)%DX
      READ(STRARRAY(4), '(ES14.3)') LINESOURCES(N_LINESOURCES)%DY
      READ(STRARRAY(5), '(ES14.3)') LINESOURCES(N_LINESOURCES)%NRHO
      READ(STRARRAY(6), '(ES14.3)') LINESOURCES(N_LINESOURCES)%UX
      READ(STRARRAY(7), '(ES14.3)') LINESOURCES(N_LINESOURCES)%UY
      READ(STRARRAY(8), '(ES14.3)') LINESOURCES(N_LINESOURCES)%UZ
      READ(STRARRAY(9), '(ES14.3)') LINESOURCES(N_LINESOURCES)%TTRA
      READ(STRARRAY(10),'(ES14.3)') LINESOURCES(N_LINESOURCES)%TROT
      READ(STRARRAY(10),'(ES14.3)') LINESOURCES(N_LINESOURCES)%TVIB
      READ(STRARRAY(11),'(A10)') MIX_NAME

      MIX_ID = MIXTURE_NAME_TO_ID(MIX_NAME)
      LINESOURCES(N_LINESOURCES)%MIX_ID = MIX_ID

   END SUBROUTINE DEF_LINESOURCE



   SUBROUTINE DEF_WALL(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      REAL(KIND=8) :: LENGTH

      TYPE(WALL), DIMENSION(:), ALLOCATABLE :: TEMP_WALLS

      
      ALLOCATE(TEMP_WALLS(N_WALLS+1))
      TEMP_WALLS(1:N_WALLS) = WALLS
      CALL MOVE_ALLOC(TEMP_WALLS, WALLS)
      N_WALLS = N_WALLS + 1



      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      READ(STRARRAY(1), '(ES14.3)') WALLS(N_WALLS)%CX
      READ(STRARRAY(2), '(ES14.3)') WALLS(N_WALLS)%CY
      READ(STRARRAY(3), '(ES14.3)') WALLS(N_WALLS)%DX
      READ(STRARRAY(4), '(ES14.3)') WALLS(N_WALLS)%DY
      READ(STRARRAY(5), '(ES14.3)') WALLS(N_WALLS)%TEMP
      READ(STRARRAY(6), *) WALLS(N_WALLS)%SPECULAR
      READ(STRARRAY(7), *) WALLS(N_WALLS)%DIFFUSE
      READ(STRARRAY(8), *) WALLS(N_WALLS)%POROUS
      READ(STRARRAY(9), '(ES14.3)') WALLS(N_WALLS)%TRANSMISSIVITY
      READ(STRARRAY(10), *) WALLS(N_WALLS)%REACT

      LENGTH = SQRT(WALLS(N_WALLS)%DX**2 + WALLS(N_WALLS)%DY**2)

      WALLS(N_WALLS)%NORMX = WALLS(N_WALLS)%DY/LENGTH
      WALLS(N_WALLS)%NORMY = - WALLS(N_WALLS)%DX/LENGTH

   END SUBROUTINE DEF_WALL


   SUBROUTINE READ_REACTIONS(FILENAME)

      IMPLICIT NONE

      CHARACTER*64, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in3 = 457
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER(LEN=80)   :: DEFINITION
      INTEGER :: N_STR, index
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      TYPE(REACTIONS_DATA_STRUCTURE) :: NEW_REACTION
      
      !CHARACTER*64 :: SP_NAME
      !INTEGER      :: SP_ID
      !REAL(KIND=8) :: DIAM
      !REAL(KIND=8) :: OMEGA
      !REAL(KIND=8) :: TREF
      !REAL(KIND=8) :: ALPHA
      !REAL(KIND=8) :: PI, KB

      !INTEGER      :: IS, JS
      !REAL(KIND=8) :: M1, M2, MRED

      !PI   = 3.141593
      !KB = 1.38064852E-23

      ! Open input file for reading
      OPEN(UNIT=in3,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, reactions definition file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in3,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction components line         
         CALL STRIP_COMMENTS(DEFINITION, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

         IF (N_STR == 0) CYCLE ! This is an empty line
         IF (STRARRAY(2) .NE. '+' .OR. STRARRAY(4) .NE. '-->' .OR. STRARRAY(6) .NE. '+') THEN
            CALL ERROR_ABORT('Attention, format is not respected in reactions file.')
         END IF

         IF (N_STR .GE. 7) THEN
            NEW_REACTION%N_PROD = 2
            NEW_REACTION%R1_SP_ID = REACTION_SPECIES_NAME_TO_ID(STRARRAY(1))
            NEW_REACTION%R2_SP_ID = REACTION_SPECIES_NAME_TO_ID(STRARRAY(3))
            NEW_REACTION%P1_SP_ID = REACTION_SPECIES_NAME_TO_ID(STRARRAY(5))
            NEW_REACTION%P2_SP_ID = REACTION_SPECIES_NAME_TO_ID(STRARRAY(7))
         END IF

         IF (N_STR .GE. 9) THEN
            IF (STRARRAY(8) .NE. '+') THEN
               CALL ERROR_ABORT('Attention, format is not respected in reactions file.')
            END IF
            NEW_REACTION%N_PROD = 3
            NEW_REACTION%P3_SP_ID = REACTION_SPECIES_NAME_TO_ID(STRARRAY(9))
         END IF
      
         READ(in3,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction parameters line
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)
         READ(STRARRAY(1), '(ES14.3)') NEW_REACTION%A
         READ(STRARRAY(2), '(ES14.3)') NEW_REACTION%N
         READ(STRARRAY(3), '(ES14.3)') NEW_REACTION%EA

         IF (ReasonEOF < 0) EXIT ! End of file reached
         
         ALLOCATE(TEMP_REACTIONS(N_REACTIONS+1)) ! Append the mixture to the list
         TEMP_REACTIONS(1:N_REACTIONS) = REACTIONS
         CALL MOVE_ALLOC(TEMP_REACTIONS, REACTIONS)
         N_REACTIONS = N_REACTIONS + 1
         REACTIONS(N_REACTIONS) = NEW_REACTION

         DEALLOCATE(STRARRAY)

      END DO
      
      CLOSE(in3) ! Close input file

      IF (PROC_ID == 0) THEN
         DO index = 1, N_REACTIONS
            IF (REACTIONS(index)%N_PROD == 2) THEN
               WRITE(*,*) 'Reaction ', index, 'Has 2 reactants with ids:', REACTIONS(index)%R1_SP_ID, ' and ', & 
               REACTIONS(index)%R2_SP_ID
               WRITE(*,*) REACTIONS(index)%N_PROD, 'products with ids:',  REACTIONS(index)%P1_SP_ID, ' and ', &
               REACTIONS(index)%P2_SP_ID
               WRITE(*,*) 'Parameters:', REACTIONS(index)%A, REACTIONS(index)%N, REACTIONS(index)%EA
            ELSE IF (REACTIONS(index)%N_PROD == 3) THEN
               WRITE(*,*) 'Reaction ', index, 'Has 2 reactants with ids:', REACTIONS(index)%R1_SP_ID, ' and ', &
               REACTIONS(index)%R2_SP_ID
               WRITE(*,*) REACTIONS(index)%N_PROD, 'products with ids:',  REACTIONS(index)%P1_SP_ID, ', ', &
               REACTIONS(index)%P2_SP_ID, ' and ', REACTIONS(index)%P3_SP_ID
               WRITE(*,*) 'Parameters:', REACTIONS(index)%A, REACTIONS(index)%N, REACTIONS(index)%EA
            END IF
         END DO
      END IF

   END SUBROUTINE READ_REACTIONS



   SUBROUTINE READ_WALL_REACTIONS(FILENAME)

      IMPLICIT NONE

      CHARACTER*64, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in4 = 5657
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER(LEN=80)   :: DEFINITION
      INTEGER :: N_STR, I
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      TYPE(WALL_REACTIONS_DATA_STRUCTURE) :: NEW_REACTION
      
      ! Open input file for reading
      OPEN(UNIT=in4,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, wall reactions definition file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in4,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction components line         
         CALL STRIP_COMMENTS(DEFINITION, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

         IF (N_STR == 0) CYCLE ! This is an empty line
         IF (STRARRAY(2) .NE. '-->' ) THEN
            CALL ERROR_ABORT('Attention, format is not respected in wall reactions file.')
         END IF

         NEW_REACTION%R_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1))
         IF (STRARRAY(3) .EQ. 'none' ) THEN
            NEW_REACTION%N_PROD = 0
            IF (N_STR .NE. 3) THEN
               CALL ERROR_ABORT('Attention, format is not respected in wall reactions file.')
            END IF
         ELSE
            NEW_REACTION%N_PROD = (N_STR-1)/2
            DO I = 1, NEW_REACTION%N_PROD
               IF (I == 1) THEN
                  NEW_REACTION%P1_SP_ID = SPECIES_NAME_TO_ID( STRARRAY(3) )
               ELSE IF (I == 2) THEN
                  NEW_REACTION%P2_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(5))
                  IF (STRARRAY(4) .NE. '+' ) THEN
                     CALL ERROR_ABORT('Attention, format is not respected in wall reactions file.')
                  END IF
               ELSE
                  ! Only acceps up to two products.
                  CALL ERROR_ABORT('Attention, format is not respected in wall reactions file.')
               END IF
            END DO
         END IF
      
         READ(in4,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction parameters line
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)
         READ(STRARRAY(1), '(ES14.3)') NEW_REACTION%PROB

         IF (ReasonEOF < 0) EXIT ! End of file reached
         
         ALLOCATE(TEMP_WALL_REACTIONS(N_WALL_REACTIONS+1)) ! Append the mixture to the list
         TEMP_WALL_REACTIONS(1:N_WALL_REACTIONS) = WALL_REACTIONS
         CALL MOVE_ALLOC(TEMP_WALL_REACTIONS, WALL_REACTIONS)
         N_WALL_REACTIONS = N_WALL_REACTIONS + 1
         WALL_REACTIONS(N_WALL_REACTIONS) = NEW_REACTION

         DEALLOCATE(STRARRAY)

      END DO
      
      CLOSE(in4) ! Close input file

      ! DO I = 1, N_WALL_REACTIONS
      !    WRITE(*,*) 'Wall reaction ', I, 'Has 1 reactant with id:', WALL_REACTIONS(I)%R_SP_ID
      !    WRITE(*,*) WALL_REACTIONS(I)%N_PROD, 'products with ids:',  WALL_REACTIONS(I)%P1_SP_ID, ' and ', WALL_REACTIONS(I)%P2_SP_ID
      !    WRITE(*,*) 'Parameters:', WALL_REACTIONS(I)%PROB
      ! END DO

   END SUBROUTINE READ_WALL_REACTIONS


   SUBROUTINE READ_GRID_FILE(FILENAME)

      IMPLICIT NONE

      CHARACTER*64, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in5 = 5557
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER :: I, J
      
      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, grid definition file not found! ABORTING.')
      ENDIF


      READ(in5,*, IOSTAT=ReasonEOF) NX ! Read number of x cells
      IF (ReasonEOF < 0) CALL ERROR_ABORT('Attention, format error in grid definition file! ABORTING.')
      ALLOCATE(XCOORD(NX+1))
      READ(in5,*, IOSTAT=ReasonEOF) XCOORD ! Read coords of x cells boundaries
      XMIN = XCOORD(1)
      XMAX = XCOORD(NX+1)
      ALLOCATE(XSIZE(NX))
      DO I = 1, NX
         XSIZE(I) = XCOORD(I+1)-XCOORD(I)
      END DO
      IF (ReasonEOF < 0) CALL ERROR_ABORT('Attention, format error in grid definition file! ABORTING.')

      READ(in5,*, IOSTAT=ReasonEOF) NY ! Read number of y cells
      IF (ReasonEOF < 0) CALL ERROR_ABORT('Attention, format error in grid definition file! ABORTING.')
      ALLOCATE(YCOORD(NY+1))
      READ(in5,*, IOSTAT=ReasonEOF) YCOORD ! Read coords of y cells boundaries
      YMIN = YCOORD(1)
      YMAX = YCOORD(NY+1)
      ALLOCATE(YSIZE(NY))
      DO I = 1, NY
         YSIZE(I) = YCOORD(I+1)-YCOORD(I)
      END DO
      IF (ReasonEOF < 0) CALL ERROR_ABORT('Attention, format error in grid definition file! ABORTING.')

      READ(in5,*, IOSTAT=ReasonEOF) NZ ! Read number of z cells
      IF (ReasonEOF < 0) CALL ERROR_ABORT('Attention, format error in grid definition file! ABORTING.')
      READ(in5,*, IOSTAT=ReasonEOF) ZMIN, ZMAX ! Read coords of z cells boundaries
      IF (ReasonEOF < 0) CALL ERROR_ABORT('Attention, format error in grid definition file! ABORTING.')
   
      ALLOCATE(CELL_VOLUMES(NX*NY))
      DO I = 1, NX
         DO J = 1, NY
            CELL_VOLUMES(I + NX*(J-1)) = XSIZE(I)*YSIZE(J)*(ZMAX-ZMIN)
         END DO
      END DO
      ! Perform all the checks on these arrays!
      CLOSE(in5)

   END SUBROUTINE READ_GRID_FILE


   SUBROUTINE INIT_QUADTREE
      
      IMPLICIT NONE

      INTEGER :: JC

      IF (NZ .NE. 1) CALL ERROR_ABORT('Error! Number of cells along Z must be 1.')
      ALLOCATE(QUADTREE_ROOT(NX*NY))
      DO JC = 1, NX*NY
         QUADTREE_ROOT(JC)%LEVEL = 1
         QUADTREE_ROOT(JC)%ISLEAF = .TRUE.
         !QUADTREE_ROOT(JC)%PARENT => QUADTREE_ROOT
      END DO

   END SUBROUTINE INIT_QUADTREE

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

      INTEGER            :: IP, NP_INIT
      REAL(KIND=8)       :: Xp, Yp, Zp, VXp, VYp, VZp, EROT, EVIB
      INTEGER            :: CID

      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
      REAL(KIND=8)  :: M
      INTEGER       :: S_ID, i
            
      ! Print message 
      CALL ONLYMASTERPRINT1(PROC_ID, '> SEEDING INITIAL PARTICLES IN THE DOMAIN...')
     
      ! Compute number of particles to be seeded
      !NP_INIT = NINT(NRHO_INIT/FNUM*(XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN))

      !CALL ONLYMASTERPRINT2(PROC_ID, '   Particles to be seeded', REAL(NP_INIT, KIND=8))
 
      ! Compute number of particles to be seeded by every process
      !NPPP_INIT = INT(NP_INIT/N_MPI_THREADS) 

      ! Create particles in the domain
      DO i = 1, MIXTURES(MIX_INIT)%N_COMPONENTS
         
         S_ID = MIXTURES(MIX_INIT)%COMPONENTS(i)%ID
         ! Compute number of particles of this species per process to be created.
         NP_INIT = NINT(NRHO_INIT/(FNUM*SPECIES(S_ID)%SPWT)*(XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN)* &
                       MIXTURES(MIX_INIT)%COMPONENTS(i)%MOLFRAC/N_MPI_THREADS)
         IF (NP_INIT == 0) CYCLE
         DO IP = 1, NP_INIT

            ! Create particle position randomly in the domain
            XP = XMIN + (XMAX-XMIN)*rf()
            !XP = XMIN + 0.5*(XMAX-XMIN)*rf()
            !IF (XP .GT. 0.) CYCLE
            YP = YMIN + (YMAX-YMIN)*rf()
            ZP = ZMIN + (ZMAX-ZMIN)*rf()

            ! Chose particle species based on mixture specifications
            ! RANVAR = rf()
            ! DO i = 1, MIXTURES(MIX_INIT)%N_COMPONENTS
            !    RANVAR = RANVAR - MIXTURES(MIX_INIT)%COMPONENTS(i)%MOLFRAC
            !    IF (RANVAR < 0.) THEN
            !       S_ID = MIXTURES(MIX_INIT)%COMPONENTS(i)%ID
            !       EXIT
            !    END IF
            ! END DO
         
            ! Assign velocity and energy following a Boltzmann distribution
            M = SPECIES(S_ID)%MOLECULAR_MASS
            IF (S_ID == 15) THEN
               CALL MAXWELL(UX_INIT, UY_INIT, UZ_INIT, &
                        3.d-2, 3.d-2, 3.d-2, &
                        VXP, VYP, VZP, M)
            ELSE
               CALL MAXWELL(UX_INIT, UY_INIT, UZ_INIT, &
                           TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT, &
                           VXP, VYP, VZP, M)
            END IF

            CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, TROT_INIT, EROT)
            CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, TVIB_INIT, EVIB)

            CALL CELL_FROM_POSITION(XP,YP,  CID) ! Find cell containing particle

            CALL INIT_PARTICLE(XP,YP,ZP,VXP,VYP,VZP,EROT,EVIB,S_ID,CID,DT, particleNOW) ! Save in particle
            CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
         END DO
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
         CALL ERROR_ABORT('ERROR! Number of cells along X or Y smaller than one. ABORTING!')
      END IF 

      IF (NZ /= 1) THEN
         CALL ERROR_ABORT('ERROR! Number of cells along Z different than 1 is not supported. ABORTING!')
      END IF 
     
      IF (DIMS == 1 .AND. NY /= 1) THEN
         CALL ERROR_ABORT('ERROR! Number of cells along Y different than 1 for a 2d simulation. ABORTING!')
      END IF
      ! ------------ Check geometrical symmetries and flags ------
      ! For axisymmetric simulations, Y cannot be periodic (YMIN is the symmetry axis)
      IF (BOOL_AXI) THEN 
         IF (BOOL_Y_PERIODIC) THEN
            CALL ERROR_ABORT('ERROR! For axisymmetric simulations, Y cannot be periodic.')
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
 
      ! ~~~~~~~~ Initial allocation for vector of particles in each process ~~~~~~~~
      NP_PROC = 0
      ALLOCATE(particles(0))

      ! ~~~~~~~~ Initialize RNG (random number generator) ~~~~~~~~~~

      ! Create local seed for RNG (different for each processor)
      RNG_SEED_LOCAL = RNG_SEED_GLOBAL*(PROC_ID + 1) 

      CALL init_genrand64(RNG_SEED_LOCAL) ! Seed Mersenne Twister PRNG with local quantity

      ! ~~~~~~~~ Additional variables for MPI domain decomposition ~~~~~~
      IF (DOMPART_TYPE == 1) THEN ! "blocks" domain partition
         ! Compute DX_BLOCKS, DY_BLOCKS
         DX_BLOCKS = (XMAX - XMIN)/N_BLOCKS_X
         DY_BLOCKS = (YMAX - YMIN)/N_BLOCKS_Y
      END IF

      CELL_VOL = (XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN)/(NX*NY)

   END SUBROUTINE INITVARIOUS

   
   SUBROUTINE INITFIELDS

      IMPLICIT NONE

      NPX = NX + 1
      IF (DIMS == 2) THEN
         NPY = NY + 1
      ELSE
         NPY = 1
      END IF
      ALLOCATE(E_FIELD(0:NPX-1, 0:NPY-1, 3))

   END SUBROUTINE INITFIELDS

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

      REAL(KIND=8) :: BETA, FLUXBOUND, NtotINJECT, Snow
      REAL(KIND=8) :: U_NORM, S_NORM, FLUXLINESOURCE, LINELENGTH, NORMX, NORMY
      REAL(KIND=8) :: PI, PI2  

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: nfs_LINE
 
      REAL(KIND=8)  :: M, FRAC, KB
      INTEGER :: N_COMP, IS, S_ID, ILINE
    
      PI   = 3.141593
      PI2  = 2.*PI

      KB = 1.38064852E-23

      ! Injection from boundaries
      IF (BOOL_BOUNDINJECT) THEN 
         N_COMP = MIXTURES(MIX_BOUNDINJECT)%N_COMPONENTS

         ALLOCATE(nfs_XMIN(N_COMP))
         ALLOCATE(nfs_XMAX(N_COMP))
         ALLOCATE(nfs_YMIN(N_COMP))
         ALLOCATE(nfs_YMAX(N_COMP))

         DO IS = 1, N_COMP ! Loop on mixture components
            S_ID = MIXTURES(MIX_BOUNDINJECT)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            FRAC = MIXTURES(MIX_BOUNDINJECT)%COMPONENTS(IS)%MOLFRAC
            BETA = 1./SQRT(2.*KB/M*TTRA_BOUND) ! sqrt(M/(2*kB*T)), it's the Maxwellian std dev

            ! +++++++++++ Lower x ++++++++++++
            IF (BOOL_INJ_XMIN) THEN

               S_NORM_XMIN = UX_BOUND*BETA ! Molecular speed ratio normal to boundary
               Snow        = S_NORM_XMIN   ! temp variable

               FLUXBOUND = NRHO_BOUNDINJECT*FRAC/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) & 
                                       + SQRT(PI)*Snow*(1.+ERF1(Snow)))    ! Tot number flux emitted [1/(s m^2)]
               NtotINJECT = FLUXBOUND*(YMAX-YMIN)*(ZMAX-ZMIN)*DT/(FNUM*SPECIES(S_ID)%SPWT)        ! Tot num of particles to be injected
               nfs_XMIN(IS) = NtotINJECT/REAL(N_MPI_THREADS,KIND=8)  ! Particles injected by each proc
               WRITE(*,*)'ON XLO ', nfs_XMIN(1)
               !ACCA = SQRT(Snow**2+2.)                                       ! Tmp variable
               !KAPPA_XMIN = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA)) ! global variable

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

               FLUXBOUND   = NRHO_BOUNDINJECT*FRAC/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) &
                                       + SQRT(PI)*Snow*(1.+ERF1(Snow)))      ! Tot number flux emitted 
               NtotINJECT  = FLUXBOUND*(YMAX-YMIN)*(ZMAX-ZMIN)*DT/(FNUM*SPECIES(S_ID)%SPWT)         ! Tot num of particles to be injected
               nfs_XMAX(IS)    = NtotINJECT/REAL(N_MPI_THREADS,KIND=8) ! Particles injected by each proc
               WRITE(*,*)'ON XHI ', nfs_XMAX(1)
               !ACCA        = SQRT(Snow**2+2.)                                  ! Tmp variable
               !KAPPA_XMAX  = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA))  ! global variable

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

               FLUXBOUND   = NRHO_BOUNDINJECT*FRAC/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) &
                                       + SQRT(PI)*Snow*(1.+ERF1(Snow)))      ! Tot number flux emitted
               NtotINJECT  = FLUXBOUND*(XMAX-XMIN)*(ZMAX-ZMIN)*DT/(FNUM*SPECIES(S_ID)%SPWT)         ! Tot num of particles to be injected
               nfs_YMIN(IS)    = NtotINJECT/REAL(N_MPI_THREADS,KIND=8) ! Particles injected by each proc
               WRITE(*,*)'ON YLO ', nfs_YMIN(1)
               !ACCA        = SQRT(Snow**2+2.)                                  ! Tmp variable
               !KAPPA_YMIN  = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA))  ! global variable

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

               FLUXBOUND   = NRHO_BOUNDINJECT*FRAC/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) &
                                       + SQRT(PI)*Snow*(1.+ERF1(Snow)))      ! Tot number flux emitted 
               NtotINJECT  = FLUXBOUND*(XMAX-XMIN)*(ZMAX-ZMIN)*DT/(FNUM*SPECIES(S_ID)%SPWT)         ! Tot num of particles to be injected
               nfs_YMAX(IS)    = NtotINJECT/REAL(N_MPI_THREADS,KIND=8) ! Particles injected by each proc
               WRITE(*,*)'ON YHI ', nfs_YMAX(1)
               !ACCA        = SQRT(Snow**2+2.)                                  ! Tmp variable
               !KAPPA_YMAX  = 2./(Snow+ACCA) * EXP(0.5 + 0.5*Snow*(Snow-ACCA))  ! global variable

               ! Check: if we are hypersonic and stuff exits domain print a warning
               IF (UY_BOUND-3./BETA .GE. 0.) THEN 
                  CALL ONLYMASTERPRINT1(PROC_ID, '$$$ Warning! Hypersonic boundary! Almost no &
                                                      &particles will be emitted at YMAX.')
               END IF

            END IF
         END DO
      END IF

      ! =====================================================
      ! Injection from line sources
      ! Can have multiple ones and can be used as boundary injection

      DO ILINE = 1, N_LINESOURCES ! Loop on line sources
         
         N_COMP = MIXTURES(LINESOURCES(ILINE)%MIX_ID)%N_COMPONENTS
         
         ALLOCATE(nfs_LINE(N_COMP))

         DO IS = 1, N_COMP ! Loop on mixture components
            ! The species ID of the component
            S_ID = MIXTURES(LINESOURCES(ILINE)%MIX_ID)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            FRAC = MIXTURES(LINESOURCES(ILINE)%MIX_ID)%COMPONENTS(IS)%MOLFRAC
            BETA = 1./SQRT(2.*KB/M*LINESOURCES(ILINE)%TTRA) ! sqrt(M/(2*kB*T)), it's the Maxwellian std dev
         
            LINELENGTH = SQRT(LINESOURCES(ILINE)%DX**2 + LINESOURCES(ILINE)%DY**2)

            NORMX = LINESOURCES(ILINE)%DY/LINELENGTH
            LINESOURCES(ILINE)%NORMX = NORMX
            NORMY = - LINESOURCES(ILINE)%DX/LINELENGTH
            LINESOURCES(ILINE)%NORMY = NORMY

            U_NORM = LINESOURCES(ILINE)%UX*NORMX + LINESOURCES(ILINE)%UY*NORMY ! Molecular speed ratio normal to boundary
            S_NORM = U_NORM*BETA
            LINESOURCES(ILINE)%S_NORM = S_NORM
            Snow   = S_NORM     ! temp variable

            FLUXLINESOURCE   = LINESOURCES(ILINE)%NRHO*FRAC/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) &
                                    + SQRT(PI)*Snow*(1.+ERF1(Snow)))      ! Tot number flux emitted 
            NtotINJECT  = FLUXLINESOURCE*LINELENGTH*(ZMAX-ZMIN)*DT/(FNUM*SPECIES(S_ID)%SPWT)         ! Tot num of particles to be injected

            nfs_LINE(IS)    = NtotINJECT/REAL(N_MPI_THREADS,KIND=8) ! Particles injected by each proc
            
            ! Check: if we are hypersonic and stuff exits domain print a warning
            IF (S_NORM .LE. -3.) THEN 
               CALL ONLYMASTERPRINT1(PROC_ID, '$$$ Warning! Hypersonic boundary! Almost no &
                                                   &particles will be emitted at line source.')
            END IF
         END DO

         CALL MOVE_ALLOC(nfs_LINE, LINESOURCES(ILINE)%nfs)
      
      END DO


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


   SUBROUTINE INITREACTIONS
   ! Here we precompute the TCE coefficients, with the input collision parameters

      IMPLICIT NONE

      REAL(KIND=8) :: OMEGA, TREF, SIGMA, ZETA, ETA, LAMBDA, EPS, MR, EA, M1, M2
      REAL(KIND=8) :: C1, C2, C3
      REAL(KIND=8) :: PI, KB
      INTEGER      :: JR, R1_SP_ID, R2_SP_ID

      PI   = 3.141593
      KB = 1.38064852E-23

      DO JR = 1, N_REACTIONS
         R1_SP_ID = REACTIONS(JR)%R1_SP_ID
         R2_SP_ID = REACTIONS(JR)%R2_SP_ID

         ZETA  = 0.5*(3. + SPECIES(R1_SP_ID)%VIBDOF + 3. + SPECIES(R2_SP_ID)%VIBDOF)
         OMEGA = 0.5*(SPECIES(R1_SP_ID)%OMEGA + SPECIES(R2_SP_ID)%OMEGA)
         SIGMA = 0.25*PI*(SPECIES(R1_SP_ID)%DIAM + SPECIES(R2_SP_ID)%DIAM)**2
         TREF  = 0.5*(SPECIES(R1_SP_ID)%TREF + SPECIES(R2_SP_ID)%TREF)
         M1 = SPECIES(R1_SP_ID)%MOLECULAR_MASS
         M2 = SPECIES(R2_SP_ID)%MOLECULAR_MASS
         MR    = M1*M2/(M1+M2)
         LAMBDA = REACTIONS(JR)%A
         ETA = REACTIONS(JR)%N
         EA = REACTIONS(JR)%EA

         IF (R1_SP_ID == R1_SP_ID) THEN
            EPS = 2.
         ELSE
            EPS = 1.
         END IF

         C1 = SQRT(PI)*EPS*LAMBDA/(2.*SIGMA) * GAMMA(ZETA+2.5-OMEGA)/GAMMA(ZETA+ETA+1.5) &
              * SQRT(MR/(2.*KB*TREF)) * TREF**(1.-OMEGA) / KB**(ETA-1.+OMEGA)
         C2 = ETA - 1. + OMEGA
         C3 = ZETA + 1.5 - OMEGA

         REACTIONS(JR)%C1 = C1
         REACTIONS(JR)%C2 = C2
         REACTIONS(JR)%C3 = C3
      END DO

   END SUBROUTINE INITREACTIONS


   INTEGER FUNCTION SPECIES_NAME_TO_ID(NAME)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER                       :: INDEX, MATCH
      MATCH = -1
      DO INDEX = 1, N_SPECIES
         IF (SPECIES(INDEX)%NAME == NAME) MATCH = INDEX
      END DO


      SPECIES_NAME_TO_ID = MATCH

   END FUNCTION SPECIES_NAME_TO_ID



   INTEGER FUNCTION REACTION_SPECIES_NAME_TO_ID(NAME)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER                       :: INDEX, MATCH
      MATCH = -1
      DO INDEX = 1, N_SPECIES
         IF (SPECIES(INDEX)%NAME == NAME) MATCH = INDEX
      END DO
      IF (NAME == 'M') THEN
         MATCH = 0
      END IF


      REACTION_SPECIES_NAME_TO_ID = MATCH

   END FUNCTION REACTION_SPECIES_NAME_TO_ID



   INTEGER FUNCTION MIXTURE_NAME_TO_ID(NAME)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER                       :: INDEX, MATCH
      MATCH = -1
      DO INDEX = 1, N_MIXTURES
         IF (MIXTURES(INDEX)%NAME == NAME) MATCH = INDEX
      END DO

      IF (MATCH .EQ. -1) THEN
         WRITE(*,*) 'Error! Mixture name not found.'
      END IF

      MIXTURE_NAME_TO_ID = MATCH

   END FUNCTION MIXTURE_NAME_TO_ID


END MODULE initialization
