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

      CHARACTER*512      :: bufferchar ! working variable
      CHARACTER*512      :: MIXTURE_DEFINITION, VSS_PARAMS_FILENAME, LINESOURCE_DEFINITION
      CHARACTER*64       :: MIX_INIT_NAME, MIX_BOUNDINJECT_NAME, DSMC_COLL_MIX_NAME

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
         IF (line=='Domain_specular:')         READ(in1,*) BOOL_XMIN_SPECULAR, BOOL_XMAX_SPECULAR, &
                                                           BOOL_YMIN_SPECULAR, BOOL_YMAX_SPECULAR
         IF (line=='Number_of_cells:')         READ(in1,*) NX, NY, NZ

         ! ~~~~~~~~~~~~~  Numerical settings  ~~~~~~~~~~~~~~~~~
         IF (line=='Fnum:')                    READ(in1,*) FNUM
         IF (line=='Timestep:')                READ(in1,*) DT
         IF (line=='Number_of_timesteps:')     READ(in1,*) NT
         IF (line=='RNG_seed:')                READ(in1,*) RNG_SEED_GLOBAL

         ! ~~~~~~~~~~~~~  Dump quantities  ~~~~~~~~~~~~~~~~
         ! Dump grid-averaged quantities
         IF (line=='Dump_grid_avgevery:')      READ(in1,*) DUMP_GRID_AVG_EVERY
         IF (line=='Dump_grid_start:')         READ(in1,*) DUMP_GRID_START
         IF (line=='Dump_grid_numavgs:')       READ(in1,*) DUMP_GRID_N_AVG

         ! Dump global moments
         IF (line=='Dump_glob_mom_every:')     READ(in1,*) DUMP_GLOB_MOM_EVERY

         ! Dump particles
         IF (line=='Dump_part_every:')         READ(in1,*) DUMP_PART_EVERY

         IF (line=='Dump_part_directory:')     THEN
            bufferchar = ''
            READ(in1,'(A)') bufferchar
            DUMP_PART_PATH = TRIM(ADJUSTL(bufferchar))
         END IF

         ! ~~~~~~~~~~~~~  Electric and magnetic fields ~~~~~~~~~~~~~

         IF (line=='B_field_type:')  READ(in1,*) B_FIELD_TYPE
         IF (line=='E_field_type:')  READ(in1,*) E_FIELD_TYPE
         IF (line=='B_field_params_uniform:')  READ(in1,*) B_UNIFORM_VAL 
         IF (line=='E_field_params_uniform:')  READ(in1,*) E_UNIFORM_VAL
         IF (line=='B_field_params_sin:')  READ(in1,*) B_SIN_AMPL, B_SIN_k, B_SIN_omega, B_SIN_PHASE_rad 
         IF (line=='E_field_params_sin:')  READ(in1,*) E_SIN_AMPL, E_SIN_k, E_SIN_omega, E_SIN_PHASE_rad 

         ! ~~~~~~~~~~~~~  Multispecies ~~~~~~~~~~~~~~~
         IF (line=='Species_file:') THEN ! Attention, must deal with paths, "./" and remove leading spaces
            bufferchar = '' ! reset
            READ(in1,'(A)') bufferchar
            SPECIES_FILENAME = TRIM(ADJUSTL(bufferchar))
            CALL READ_SPECIES
         END IF

         IF (line=='Def_mixture:') THEN
            READ(in1,'(A)') MIXTURE_DEFINITION
            CALL DEF_MIXTURE(MIXTURE_DEFINITION)
         END IF

         ! ~~~~~~~~~~~~~  Collision type  ~~~~~~~~~~~~~~~~~
         IF (line=='Collision_type:')          READ(in1,*) COLLISION_TYPE

         ! MCC collisions
         IF (line=='MCC_background_dens:')     READ(in1,*) MCC_BG_DENS
         IF (line=='MCC_cross_section:')       READ(in1,*) MCC_SIGMA

         ! BGK collisions
         IF (line=='BGK_cross_section:')       READ(in1,*) BGK_SIGMA
         IF (line=='BGK_background_mass:')     READ(in1,*) BGK_BG_MASS
         IF (line=='BGK_background_dens:')     READ(in1,*) BGK_BG_DENS
          
         IF (line=='BGK_model_type:')          THEN
            bufferchar = '' ! reset
            READ(in1,'(A)') bufferchar 
            BGK_MODEL_TYPE = TRIM(ADJUSTL(bufferchar))
         END IF            

         ! DSMC collisions 
         IF (line=='VSS_parameters_file:')     THEN
            bufferchar = '' ! reset
            READ(in1,'(A)') bufferchar
            VSS_PARAMS_FILENAME = TRIM(ADJUSTL(bufferchar))
            CALL READ_VSS(VSS_PARAMS_FILENAME)
         END IF

         IF (line=='DSMC_collisions_mixture:') THEN
            bufferchar = '' ! reset
            READ(in1,'(A)') bufferchar
            DSMC_COLL_MIX_NAME = TRIM(ADJUSTL(bufferchar))
            DSMC_COLL_MIX = MIXTURE_NAME_TO_ID(DSMC_COLL_MIX_NAME)
         END IF

         ! CUSTOM collisions
         IF (line=='CUSTOM_background_mass:')        READ(in1,*) CUSTOM_BG_MASS
         IF (line=='CUSTOM_background_dens:')        READ(in1,*) CUSTOM_BG_DENS
         IF (line=='CUSTOM_excitation_energy_eV:')   READ(in1,*) DeltaE_exc_eV
         IF (line=='CUSTOM_ionization_energy_eV:')   READ(in1,*) DeltaE_iz_eV

         ! ~~~~~~~~~~~~~  Initial particles seeding  ~~~~~~~~~~~~~~~~~
         IF (line=='Initial_particles_bool:')  READ(in1,*) BOOL_INITIAL_SEED
         IF (line=='Initial_particles_dens:')  READ(in1,*) NRHO_INIT
         IF (line=='Initial_particles_vel:')   READ(in1,*) UX_INIT, UY_INIT, UZ_INIT
         IF (line=='Initial_particles_Ttra:')  READ(in1,*) TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT
         IF (line=='Initial_particles_Trot:')  READ(in1,*) TROT_INIT
         IF (line=='Initial_particles_mixture:') THEN
            READ(in1,*) MIX_INIT_NAME
            MIX_INIT = MIXTURE_NAME_TO_ID(MIX_INIT_NAME)
         END IF

         ! ~~~~~~~~~~~~~  Particle injection at boundaries ~~~~~~~~~~~
         IF (line=='Boundaries_inject_bool:')       READ(in1,*)BOOL_BOUNDINJECT
         IF (line=='Boundaries_inject_which_bool:') READ(in1,*)BOOL_INJ_XMIN,BOOL_INJ_XMAX,BOOL_INJ_YMIN, BOOL_INJ_YMAX
         IF (line=='Boundaries_inject_dens:')       READ(in1,*)NRHO_BOUNDINJECT
         IF (line=='Boundaries_inject_vel:')        READ(in1,*)UX_BOUND, UY_BOUND, UZ_BOUND
         IF (line=='Boundaries_inject_Ttra:')       READ(in1,*)TTRA_BOUND
         IF (line=='Boundaries_inject_Trot:')       READ(in1,*)TROT_BOUND
         IF (line=='Boundaries_inject_mixture:') THEN
            READ(in1,*) MIX_BOUNDINJECT_NAME
            MIX_BOUNDINJECT = MIXTURE_NAME_TO_ID(MIX_BOUNDINJECT_NAME)
         END IF

         ! ~~~~~~~~~~~~~  Particle injection at line source ~~~~~~~~~~~
         IF (line=='Linesource:') THEN
            READ(in1,'(A)') LINESOURCE_DEFINITION
            CALL DEF_LINESOURCE(LINESOURCE_DEFINITION)
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

         ! ~~~~ Electric and magnetic fields ~~~~
         WRITE(*,*) '  =========== Electric and magnetic fields =========='
         string = 'Electric field type:'
         WRITE(*,'(A5,A50,A64)') '     ', string, E_FIELD_TYPE
     
         ! If sine, plot parameters 
         IF (E_FIELD_TYPE .EQ. "UNIFORM") THEN

            WRITE(*,'(A5, ES14.7)') '     ', E_UNIFORM_VAL

         ELSE IF (E_FIELD_TYPE .EQ. "SINE") THEN

            string = 'Parameters:'
            WRITE(*,'(A5,ES14.7,A5,ES14.7, A5, ES14.7, A4, ES14.7, A2)') '     ', & 
            E_SIN_AMPL, ' SIN(', E_SIN_k, ' x - ', E_SIN_omega, ' t + ', E_SIN_PHASE_rad , ' )'

         END IF

         string = 'Magnetic field type:'
         WRITE(*,'(A5,A50,A64)') '     ', string, B_FIELD_TYPE

         ! If sine, plot parameters 
         IF (B_FIELD_TYPE .EQ. "UNIFORM") THEN

            WRITE(*,'(A5, ES14.7, A4)') '     ', B_UNIFORM_VAL

         ELSE IF (B_FIELD_TYPE .EQ. "SINE") THEN

            string = 'Parameters:'
            WRITE(*,'(A5,ES14.7,A5,ES14.7, A5, ES14.7, A4, ES14.7, A2)') '     ', & 
            B_SIN_AMPL, ' SIN(', B_SIN_k, ' x - ', B_SIN_omega, ' t + ', B_SIN_PHASE_rad , ' )'

         END IF

         ! ~~~~ Collisions ~~~~
         WRITE(*,*) '  =========== Collisions ============================'
         string = 'Collision type:'
         WRITE(*,'(A5,A50,A64)') '     ', string, COLLISION_TYPE 

         IF (COLLISION_TYPE == 'MCC') THEN  ! Only print this if MCC collisions are ON

            string = 'Background number density [1/m^3]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, MCC_BG_DENS

            string = 'Collision cross-section [m^2]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, MCC_SIGMA 

         ELSE IF (COLLISION_TYPE == 'BGK') THEN ! Only print this if BGK collisions are ON

            string = 'BGK collision model type:'
            WRITE(*,'(A5,A50,A20)') '     ', string, BGK_MODEL_TYPE

            string = 'Collision cross-section [m^2]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, BGK_SIGMA 

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

         ! ~~~~ Multispecies ~~~~
         WRITE(*,*) '  =========== Species ========================='
         WRITE(*,*) 'File read'
         WRITE(*,'(A30,I9)') 'Number of species present: ', N_SPECIES


         WRITE(*,*) '  =========== Mixtures ========================='
         DO j = 1, N_MIXTURES
            WRITE(*,*)'Contents of mixture ', MIXTURES(j)%NAME, 'with ', MIXTURES(j)%N_COMPONENTS, ' components'
            DO k = 1, MIXTURES(j)%N_COMPONENTS
               WRITE(*,*) 'Mix component ', MIXTURES(j)%COMPONENTS(k)%NAME, &
                        ' with frac ', MIXTURES(j)%COMPONENTS(k)%MOLFRAC, &
                        ' and ID ', MIXTURES(j)%COMPONENTS(k)%ID
            END DO
         END DO
         
      END IF

   END SUBROUTINE PRINTINPUT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE READ_SPECIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE READ_SPECIES
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: in2 = 445
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER*10 :: NAME
      REAL(KIND=8) :: MOLWT
      REAL(KIND=8) :: MOLMASS
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
         PRINT*
         WRITE(*,*)'  Attention, species definition file "', SPECIES_FILENAME, '" not found! ABORTING.'
         PRINT*
         STOP
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      N_SPECIES = 0
      DO

         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line until EoF
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         IF (LEN_TRIM(line) .NE. 0) THEN ! Skip empty lines (commented ones for example, that became empty)

            !READ(line,'(A2, ES14.3, ES14.3, I1, ES14.3, I1, ES14.3, ES14.3, ES14.3, ES14.3)') &
            READ(line,*) &
            NAME, MOLWT, MOLMASS, ROTDOF, ROTREL, VIBDOF, VIBREL, VIBTEMP, SPWT, CHARGE
         
            ALLOCATE(TEMP_SPECIES(N_SPECIES+1)) ! Append the species to the list
            TEMP_SPECIES(1:N_SPECIES) = SPECIES
            CALL MOVE_ALLOC(TEMP_SPECIES, SPECIES)
            N_SPECIES = N_SPECIES + 1
   
            SPECIES(N_SPECIES)%NAME = NAME
            SPECIES(N_SPECIES)%MOLWT = MOLWT
            SPECIES(N_SPECIES)%MOLMASS = MOLMASS
            SPECIES(N_SPECIES)%ROTDOF = ROTDOF
            SPECIES(N_SPECIES)%ROTREL = ROTREL
            SPECIES(N_SPECIES)%VIBDOF = VIBDOF
            SPECIES(N_SPECIES)%VIBREL = VIBREL
            SPECIES(N_SPECIES)%VIBTEMP = VIBTEMP
            SPECIES(N_SPECIES)%SPWT = SPWT
            SPECIES(N_SPECIES)%CHARGE = CHARGE
         
         END IF

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

      PI   = 3.141593
      KB = 1.38064852E-23

      ! Open input file for reading
      OPEN(UNIT=in2,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         PRINT*
         WRITE(*,*)'  Attention, VSS parameters definition file not found! ABORTING.'
         PRINT*
         STOP
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         !READ(line,'(A2, ES14.3, ES14.3, I1, ES14.3, I1, ES14.3, ES14.3, ES14.3, ES14.3)') &
         READ(line,*) SP_NAME, DIAM, OMEGA, TREF, ALPHA
      
         SP_ID = SPECIES_NAME_TO_ID(SP_NAME)

         SPECIES(SP_ID)%DIAM  = DIAM
         SPECIES(SP_ID)%OMEGA = OMEGA
         SPECIES(SP_ID)%TREF  = TREF
         SPECIES(SP_ID)%ALPHA = ALPHA

         SPECIES(SP_ID)%SIGMA = PI*DIAM**2
         SPECIES(SP_ID)%NU    = OMEGA - 0.5
         SPECIES(SP_ID)%CREF  = SQRT(3.*KB*TREF / SPECIES(SP_ID)%MOLMASS)
         
      END DO

      CLOSE(in2) ! Close input file


   END SUBROUTINE READ_VSS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DEF_MIXTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DEF_MIXTURE(DEFINITION)
      
      IMPLICIT NONE

      CHARACTER*10 :: MIX_NAME
      CHARACTER*10 :: COMP_NAME
      REAL(KIND=8) :: MOLFRAC

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION
      INTEGER :: N_STR
      CHARACTER(LEN=80), allocatable :: STRARRAY(:)
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


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DEF_LINESOURCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE DEF_LINESOURCE(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      CHARACTER*64 :: S_NAME
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
      READ(STRARRAY(11),'(A200)') S_NAME ! Species name

      LINESOURCES(N_LINESOURCES)%S_ID = SPECIES_NAME_TO_ID(S_NAME)

   END SUBROUTINE

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
      REAL(KIND=8)  :: M
      INTEGER       :: S_ID, i
      REAL(KIND=8) :: RANVAR
            
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

         ! Chose particle species based on mixture specifications
         RANVAR = rand()
         DO i = 1, MIXTURES(MIX_INIT)%N_COMPONENTS
            RANVAR = RANVAR - MIXTURES(MIX_INIT)%COMPONENTS(i)%MOLFRAC
            IF (RANVAR < 0.) THEN
               S_ID = MIXTURES(MIX_INIT)%COMPONENTS(i)%ID
               EXIT
            END IF
         END DO
      
         ! Assign velocity and energy following a Boltzmann distribution
         M = SPECIES(S_ID)%MOLMASS
         CALL MAXWELL(UX_INIT, UY_INIT, UZ_INIT, &
                      TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT, TROT_INIT, &
                      VXP, VYP, VZP, EIP, M)

!          PRINT*, "DBDBDB PART: ", VXP, VYP, VZP

         CALL CELL_FROM_POSITION(XP,YP,  CID) ! Find cell containing particle

         CALL INIT_PARTICLE(XP,YP,ZP,VXP,VYP,VZP,EIP,S_ID,CID,DT, particleNOW) ! Save in particle
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

      ! Check that periodic faces were not also set as specular
      IF (BOOL_X_PERIODIC) THEN
        IF (BOOL_XMIN_SPECULAR .OR. BOOL_XMAX_SPECULAR) THEN
            WRITE(*,*)
            WRITE(*,*) 'ERROR! X faces were specified as periodic, but "specular reflection" was also imposed.' 
            WRITE(*,*) 'Pick only one of the two! Check the input file. ABORTING!'
            WRITE(*,*)
            STOP
        END IF
      END IF

      ! Same thing for Y
      IF (BOOL_Y_PERIODIC) THEN
        IF (BOOL_YMIN_SPECULAR .OR. BOOL_YMAX_SPECULAR) THEN
            WRITE(*,*)
            WRITE(*,*) 'ERROR! Y faces were specified as periodic, but "specular reflection" was also imposed.' 
            WRITE(*,*) 'Pick only one of the two! Check the input file. ABORTING!'
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

      ! --------------- BGK collisions -------------------
      IF (COLLISION_TYPE == "BGK") THEN

         ! For BGK you must specify what kind of BGK (classical / background)
         ! and the cross section (constant).
         ! For "background" also need to specify the background density

         ! Cross section is initialized to negative value to generate this check. 
         IF (BGK_SIGMA .LT. 0.0) THEN
            WRITE(*,*) '*************************************************************************'
            WRITE(*,*) '* ERROR! You did not specify the cross section in the input file, or'
            WRITE(*,*) '* you used a negative value! Use input file keyword "BGK_cross_section:"'
            WRITE(*,*) '* Aborting!'
            WRITE(*,*) '*************************************************************************'
            STOP
         END IF

         ! BGK_MODEL_TYPE (initialized as "NONE" to generate this check)
         IF (BGK_MODEL_TYPE == "NONE") THEN
            WRITE(*,*) '*************************************************************************'
            WRITE(*,*) '* ERROR! Must specify a BGK type! Use "BGK_model_type:" in input file.'
            WRITE(*,*) '* Aborting!'
            WRITE(*,*) '*************************************************************************'
            STOP
         END IF

         IF (BGK_MODEL_TYPE == "background") THEN

            IF (BGK_BG_MASS .LT. 0.0) THEN 
               WRITE(*,*) '************************************************************'
               WRITE(*,*) '* ERROR! BGK background mass not specified or negative! '
               WRITE(*,*) '* Check input file! Keyword: "BGK_background_mass:"'
               WRITE(*,*) '************************************************************'
               STOP
            END IF

            IF (BGK_BG_DENS .LT. 0.0) THEN 
               WRITE(*,*) '************************************************************'
               WRITE(*,*) '* ERROR! BGK background number density not specified or negative! '
               WRITE(*,*) '* Check input file! Keyword: "BGK_background_dens:"'
               WRITE(*,*) '************************************************************'
               STOP
            END IF


         END IF
      END IF

      ! CUSTOM collisions
      IF (COLLISION_TYPE == "CUSTOM_LOAD_RATES") THEN

            IF (CUSTOM_BG_MASS .LT. 0.0) THEN 
               WRITE(*,*) '************************************************************'
               WRITE(*,*) '* ERROR! CUSTOM collision background mass not specified or negative! '
               WRITE(*,*) '* Check input file! Keyword: "CUSTOM_background_mass:"'
               WRITE(*,*) '************************************************************'
               STOP
            END IF

            IF (CUSTOM_BG_DENS .LT. 0.0) THEN 
               WRITE(*,*) '************************************************************'
               WRITE(*,*) '* ERROR! CUSTOM collision background number density not specified or negative! '
               WRITE(*,*) '* Check input file! Keyword: "CUSTOM_background_dens:"'
               WRITE(*,*) '************************************************************'
               STOP
            END IF

            IF (DeltaE_exc_eV .LT. 0.0) THEN 
               WRITE(*,*) '************************************************************'
               WRITE(*,*) '* ERROR! CUSTOM collision excitation energy loss not specified or negative! '
               WRITE(*,*) '* Check input file! Keyword: "CUSTOM_excitation_energy_eV:"'
               WRITE(*,*) '************************************************************'
               STOP
            END IF

            IF (DeltaE_iz_eV .LT. 0.0) THEN 
               WRITE(*,*) '************************************************************'
               WRITE(*,*) '* ERROR! CUSTOM collision excitation energy loss not specified or negative! '
               WRITE(*,*) '* Check input file! Keyword: "CUSTOM_ionization_energy_eV:"'
               WRITE(*,*) '************************************************************'
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

      CELL_VOL = (XMAX-XMIN)*(YMAX-YMIN)/(NX*NY)

      ! ~~~~~~~~ ELECTRIC AND MAGNETIC FIELDS ~~~~~~~~~~
      IF (B_FIELD_TYPE .NE. "NONE") B_FIELD_BOOL = .TRUE.
      IF (E_FIELD_TYPE .NE. "NONE") E_FIELD_BOOL = .TRUE.

      ! ~~~~~~~~ DUMPS ~~~~~~~

      ! Master should make sure that there are the appropriate folders for dumping, and initialize file
      IF (PROC_ID .EQ. 0) THEN 

         CALL system("mkdir -p dumps") ! Make sure the dumps directory exists (only works in unix)
  
         ! ++++ Inits file to dump moments, if required ++++
         IF (DUMP_GLOB_MOM_EVERY .NE. -1) THEN ! "-1" is default, standing for "don't do this."

            OPEN(12, FILE=DUMP_GLOB_MOM_FILENAME) ! Create new file, or empty it if it exists
            WRITE(12, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            WRITE(12, *) "% Every line represents the global moments of the system, considering all particles"
            WRITE(12, *) "% in all the domain, at a certain time. In the order (SI units):"
            WRITE(12, *) "%"
            WRITE(12, *) "% Time [s]"
            WRITE(12, *) "% rho" 
            WRITE(12, *) "% ux  uy  uz" 
            WRITE(12, *) "% Pxx  Pxy  Pxz  Pyy  Pyz  Pzz" 
            WRITE(12, *) "% qx  qy  qz" 
            WRITE(12, *) "% Qxxx  Qxxy  Qxyy  Qyyy  Qyyz  Qyzz  Qzzz  Qxxz  Qxzz  Qxyz" 
            WRITE(12, *) "% Riijj" 
            WRITE(12, *) "% Rxxjj  Rxyjj  Rxzjj  Ryyjj  Ryzjj  Rzzjj" 
            WRITE(12, *) "% Sxiijj  Syiijj  Sziijj"
            WRITE(12, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            CLOSE(12)

         END IF

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

      REAL(KIND=8) :: BETA, FLUXBOUND, NtotINJECT, Snow
      ! REAL(KIND=8) :: U_NORM, S_NORM, FLUXLINESOURCE, NORMX, NORMY
      REAL(KIND=8) :: PI, PI2  

      REAL(KIND=8) :: Vth, Ainject, Linelength, FlxTot


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
            M = SPECIES(S_ID)%MOLMASS
            FRAC = MIXTURES(MIX_BOUNDINJECT)%COMPONENTS(IS)%MOLFRAC
            BETA = 1./SQRT(2.*KB/M*TTRA_BOUND) ! sqrt(M/(2*kB*T)), it's the Maxwellian std dev

            ! +++++++++++ Lower x ++++++++++++
            IF (BOOL_INJ_XMIN) THEN

               S_NORM_XMIN = UX_BOUND*BETA ! Molecular speed ratio normal to boundary
               Snow        = S_NORM_XMIN   ! temp variable

               FLUXBOUND = NRHO_BOUNDINJECT*FRAC/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) & 
                                       + SQRT(PI)*Snow*(1.+ERF1(Snow)))    ! Tot number flux emitted [1/(s m^2)]
               NtotINJECT = FLUXBOUND*(YMAX-YMIN)*(ZMAX-ZMIN)*DT/FNUM        ! Tot num of particles to be injected
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
               NtotINJECT  = FLUXBOUND*(YMAX-YMIN)*(ZMAX-ZMIN)*DT/FNUM         ! Tot num of particles to be injected
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
               NtotINJECT  = FLUXBOUND*(XMAX-XMIN)*(ZMAX-ZMIN)*DT/FNUM         ! Tot num of particles to be injected
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
               NtotINJECT  = FLUXBOUND*(XMAX-XMIN)*(ZMAX-ZMIN)*DT/FNUM         ! Tot num of particles to be injected
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
      ! ========= Injection from line sources ===============
      !
      ! We define a line source as a line in the X-Y plane (a plane in the X-Y-Z volume)
      ! AT which particles are injected following a distribution.
      ! Every MPI process injects particles, which will be exchanged afterwards.
      ! Particles are injected at both sides of the linesource, so the computation of 
      ! the number of particles to be injected is simple:
      !
      ! Gamma = n v_th / 2 = flux of particles through a UNIT surface, both going right and left 
      ! 
      ! Ndot_real = Gamma*Asurf = flux of particles trough surface of area Asurf 
      !                           In our case, Asurf = LINELENGTH * Zdomain
      ! 
      ! Ndot_simulated = Ndot_real/Fnum = particles to be emitted per unit time
      ! 
      ! NtotINJECT = Ndot_simluated*dt
      !
      ! NprocINJECT = NtotINJECT/N_processors 
      
      DO ILINE = 1, N_LINESOURCES ! Loop on line sources

         M                      = SPECIES(  LINESOURCES(ILINE)%S_ID  )%MOLMASS
         Vth                    = SQRT(8*KB*LINESOURCES(ILINE)%TTRA/(PI*M)) ! Thermal velocity
         Linelength             = SQRT(LINESOURCES(ILINE)%DX**2 + LINESOURCES(ILINE)%DY**2)
         Ainject                = Linelength*(ZMAX-ZMIN)
         FlxTot                 = LINESOURCES(ILINE)%NRHO*Vth/2
         NtotINJECT             = Ainject*FlxTot*DT/Fnum

         LINESOURCES(ILINE)%nfs  = NtotINJECT/REAL(N_MPI_THREADS,KIND=8) ! Particles injected by each proc (variable type: double!)

      END DO


   END SUBROUTINE INITINJECTION

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INITCOLLISIONS -> Initializes variables for the collisions !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE INITCOLLISIONS

      IMPLICIT NONE

      ! Set logicals
      IF (COLLISION_TYPE == "NONE") THEN
         BOOL_DSMC         = .FALSE.
         BOOL_MCC          = .FALSE.
         BOOL_BGK          = .FALSE.
         BOOL_CUSTOM_COLL  = .FALSE.

      ELSE IF (COLLISION_TYPE == "DSMC") THEN
         BOOL_DSMC         = .TRUE.
         BOOL_MCC          = .FALSE.
         BOOL_BGK          = .FALSE.
         BOOL_CUSTOM_COLL  = .FALSE.

      ELSE IF (COLLISION_TYPE == "MCC") THEN
         BOOL_DSMC         = .FALSE.
         BOOL_MCC          = .TRUE.
         BOOL_BGK          = .FALSE.
         BOOL_CUSTOM_COLL  = .FALSE.

      ELSE IF (COLLISION_TYPE == "BGK") THEN
         BOOL_DSMC         = .FALSE.
         BOOL_MCC          = .FALSE.
         BOOL_BGK          = .TRUE.
         BOOL_CUSTOM_COLL  = .FALSE.

      ELSE IF (COLLISION_TYPE == "CUSTOM_LOAD_RATES") THEN
         BOOL_DSMC         = .FALSE.
         BOOL_MCC          = .FALSE.
         BOOL_BGK          = .FALSE.
         BOOL_CUSTOM_COLL  = .TRUE.

      ELSE ! ERROR!
         CALL ONLYMASTERPRINT1(PROC_ID, '$$$ ATTENTION! Collision type in input file not recognized! ABORTING!')
         STOP

      END IF

      ! Initialize internal values for BGK collisions 
      IF (BOOL_BGK) THEN
         IF (BGK_MODEL_TYPE .EQ. "classical") THEN
            BGK_MODEL_TYPE_INT = 0
         ELSE IF (BGK_MODEL_TYPE .EQ. "background") THEN
            BGK_MODEL_TYPE_INT = 1
         END IF
      END IF

      ! Read rate files for CUSTOM_COLLISION type
      IF (BOOL_CUSTOM_COLL) THEN
         ! Initialize internal variables
         DeltaE_exc_J = 1.60217662E-19*DeltaE_exc_eV
         DeltaE_iz_J  = 1.60217662E-19*DeltaE_iz_eV

         ! Read rates
         CALL INIT_CUSTOM_COLL_RATES
      END IF

   END SUBROUTINE INITCOLLISIONS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION SPECIES_NAME_TO_ID -> Maps a species name (string) to its ID !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MIXTURE_NAME_TO_ID -> Maps a mixture name (string) to its ID !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INIT_CUSTOM_COLL_RATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE INIT_CUSTOM_COLL_RATES
 
   ! This subroutine reads files and initializes the vectors containing collision rates
 
   INTEGER            :: N_lines, io, ii
   REAL(KIND=8)       :: T_NOW, k_el_NOW, k_exc_NOW, k_iz_NOW
   CHARACTER(LEN=512) :: dummy

   WRITE(*,*) 'CUSTOM collision model: reading rates file...'
 
   ! First, count lines in the file
 
   N_lines = 0 ! Init
 
   OPEN(UNIT=99, FILE='custom_rates.dat', STATUS='old', ACTION='read')
 
   io = 0 ! init
   DO WHILE (io .EQ. 0)
       READ(99,*,IOSTAT=io)  dummy 
       IF (io > 0) THEN
          WRITE(*,*) 'Check reaction rates file.  Something went wrong. ABORTING!'
          STOP
       ELSE IF (io < 0) THEN
          WRITE(*,*)  ' Collision rates contain ', N_lines, ' lines'
       ELSE
          N_lines = N_lines + 1
       END IF
 
   END DO

   CLOSE(99)

   ! Init arrays with this length
   ALLOCATE(T_rates(N_lines))
   ALLOCATE(k_el(N_lines))
   ALLOCATE(k_exc(N_lines))
   ALLOCATE(k_iz(N_lines))
 
   ! Now again, read file and save its values
   OPEN(UNIT=99, FILE='custom_rates.dat', STATUS='old', ACTION='read')
 
   ! Read until EOF
   ii = 0
   DO
      READ(99,*,IOSTAT=io)  T_NOW,  k_el_NOW,  k_exc_NOW,  k_iz_NOW
  
      IF (io > 0) THEN

         WRITE(*,*) 'Check input.  Something was wrong'
         EXIT

      ELSE IF (io < 0) THEN ! EOF reached, quit.
    
         CLOSE(99)
         EXIT

      ELSE ! Proper values were read

         ii = ii + 1

         T_rates(ii) = T_NOW
         k_el(ii)    = k_el_NOW
         k_exc(ii)   = k_exc_NOW
         k_iz(ii)    = k_iz_NOW

      END IF
   END DO
 
  END SUBROUTINE INIT_CUSTOM_COLL_RATES

END MODULE initialization
