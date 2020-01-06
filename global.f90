! This module holds global variables

MODULE global

   USE particle
   USE mpi_common

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Particle variables and arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER :: NP_PROC ! Number of particles in local MPI process
   INTEGER :: MPI_PARTICLE_DATA_STRUCTURE ! We need this for MPI
   TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: particles

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Geometry, domain and grid !!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL      :: BOOL_AXI = .FALSE. ! Assign default value!
   INTEGER      :: NX, NY, NZ
   REAL(KIND=8) :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX ! Dimensions of domain
   REAL(KIND=8) :: CELL_VOL
   LOGICAL      :: BOOL_X_PERIODIC = .FALSE. ! Assign default value!
   LOGICAL      :: BOOL_Y_PERIODIC = .FALSE.
   LOGICAL      :: BOOL_Z_PERIODIC = .FALSE.
   LOGICAL      :: BOOL_XMIN_SPECULAR = .FALSE.
   LOGICAL      :: BOOL_XMAX_SPECULAR = .FALSE.
   LOGICAL      :: BOOL_YMIN_SPECULAR = .FALSE.
   LOGICAL      :: BOOL_YMAX_SPECULAR = .FALSE.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Numerical settings !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER         :: NT, tID
   REAL(KIND=8)    :: FNUM, DT, CURRENT_TIME
   INTEGER(KIND=8) :: RNG_SEED_GLOBAL, RNG_SEED_LOCAL

   ! ++++++ Dump quantities. Init to negative 1: if user does not define them, 
   ! ++++++ dump is not performed.

   ! Dump particles
   INTEGER        :: DUMP_PART_EVERY = -1
   CHARACTER(512) :: DUMP_PART_PATH  = "./dumps/" ! Default position

   ! Dump grid-averaged quantities
   INTEGER :: DUMP_GRID_AVG_EVERY  = -1
   INTEGER :: DUMP_GRID_START      = -1
   INTEGER :: DUMP_GRID_N_AVG      = -1

   ! Dump global moments
   INTEGER :: DUMP_GLOB_MOM_EVERY  = -1
   CHARACTER(26) :: DUMP_GLOB_MOM_FILENAME = "./dumps/global_moments.dat" ! Hard-coded

 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Electrical and magnetic fields !!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL :: B_FIELD_BOOL = .FALSE. ! Init fields boolean value
   LOGICAL :: E_FIELD_BOOL = .FALSE. ! Init fields boolean value

   CHARACTER(LEN=64) :: B_FIELD_TYPE = "NONE" ! Init field to none
   CHARACTER(LEN=64) :: E_FIELD_TYPE = "NONE" ! Init field to none

   ! Parameters for E and B fields
   REAL(KIND=8) :: B_UNIFORM_VAL = 0
   REAL(KIND=8) :: E_UNIFORM_VAL = 0

   REAL(KIND=8) :: B_SIN_AMPL      = 0
   REAL(KIND=8) :: B_SIN_k         = 0
   REAL(KIND=8) :: B_SIN_omega     = 0
   REAL(KIND=8) :: B_SIN_PHASE_rad = 0

   REAL(KIND=8) :: E_SIN_AMPL      = 0
   REAL(KIND=8) :: E_SIN_k         = 0
   REAL(KIND=8) :: E_SIN_omega     = 0
   REAL(KIND=8) :: E_SIN_PHASE_rad = 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Collisions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CHARACTER(LEN=64) :: COLLISION_TYPE = "NONE" ! Init collisions to none

   ! Variables for DSMC
   LOGICAL           :: BOOL_DSMC = .FALSE.
   INTEGER           :: DSMC_COLL_MIX
   INTEGER           :: TIMESTEP_COLL

   ! Variables for BGK
   LOGICAL           :: BOOL_BGK  = .FALSE.
   REAL(KIND=8)      :: BGK_SIGMA      = -1     ! Init value
   REAL(KIND=8)      :: BGK_BG_MASS    = -1     ! Init value
   REAL(KIND=8)      :: BGK_BG_DENS    = -1     ! Init value
   CHARACTER(LEN=64) :: BGK_MODEL_TYPE = "NONE" ! Init 
   INTEGER           :: BGK_MODEL_TYPE_INT      ! INT = "integer" or "internal variable"

   ! Variables for MCC
   LOGICAL           :: BOOL_MCC  = .FALSE.
   REAL(KIND=8)      :: MCC_BG_DENS, MCC_SIGMA

   ! Variables for custom collision model
   LOGICAL           :: BOOL_CUSTOM_COLL = .FALSE.
   REAL(KIND=8)      :: CUSTOM_BG_MASS = -1
   REAL(KIND=8)      :: CUSTOM_BG_DENS = -1
   REAL(KIND=8)      :: DeltaE_exc_J  = -1
   REAL(KIND=8)      :: DeltaE_exc_eV = -1
   REAL(KIND=8)      :: DeltaE_iz_J   = -1
   REAL(KIND=8)      :: DeltaE_iz_eV  = -1
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: T_rates, k_el, k_exc, k_iz

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Initial particles seed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL      :: BOOL_INITIAL_SEED = .FALSE. ! Assign default value!
   REAL(KIND=8) :: NRHO_INIT
   REAL(KIND=8) :: UX_INIT, UY_INIT, UZ_INIT
   REAL(KIND=8) :: TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT, TROT_INIT
   INTEGER      :: MIX_INIT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Particles injection from boundaries !!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL      :: BOOL_BOUNDINJECT = .FALSE. ! Assign default value!
   LOGICAL      :: BOOL_INJ_XMIN=.FALSE.,BOOL_INJ_XMAX=.FALSE. ! Assign default value!
   LOGICAL      :: BOOL_INJ_YMIN=.FALSE.,BOOL_INJ_YMAX=.FALSE. ! Assign default value!
   REAL(KIND=8) :: NRHO_BOUNDINJECT
   REAL(KIND=8) :: UX_BOUND, UY_BOUND, UZ_BOUND
   REAL(KIND=8) :: TTRA_BOUND, TROT_BOUND ! No TTRAX_ TTRAY_ TTRAZ_ for now! IMPLEMENT IT!
   INTEGER      :: MIX_BOUNDINJECT

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: nfs_XMIN, nfs_XMAX, nfs_YMIN, nfs_YMAX
   !INTEGER      :: nfs_XMIN, nfs_XMAX, nfs_YMIN, nfs_YMAX
   REAL(KIND=8) :: KAPPA_XMIN, KAPPA_XMAX, KAPPA_YMIN, KAPPA_YMAX
   REAL(KIND=8) :: S_NORM_XMIN, S_NORM_XMAX, S_NORM_YMIN, S_NORM_YMAX

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Particles injection from line source !!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Not used anymore
   LOGICAL      :: BOOL_LINESOURCE = .FALSE. ! Assign default value!
  
   REAL(KIND=8) :: X_LINESOURCE, Y_LINESOURCE, L_LINESOURCE
   REAL(KIND=8) :: NRHO_LINESOURCE
   REAL(KIND=8) :: UX_LINESOURCE, UY_LINESOURCE, UZ_LINESOURCE
   REAL(KIND=8) :: TTRA_LINESOURCE, TROT_LINESOURCE ! No TTRAX_ TTRAY_ TTRAZ_ for now! IMPLEMENT IT!

   INTEGER      :: nfs_LINESOURCE
   REAL(KIND=8) :: KAPPA_LINESOURCE
   REAL(KIND=8) :: S_NORM_LINESOURCE

   ! This is used instead.
   INTEGER         :: N_LINESOURCES = 0

   TYPE LINESOURCE
      REAL(KIND=8) :: CX, CY, DX, DY
      REAL(KIND=8) :: NRHO
      REAL(KIND=8) :: UX, UY, UZ
      REAL(KIND=8) :: TTRA, TROT ! No TTRAX_ TTRAY_ TTRAZ_ for now! IMPLEMENT IT!
      INTEGER      :: MIX_ID

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: nfs
      REAL(KIND=8) :: KAPPA
      REAL(KIND=8) :: U_NORM
      REAL(KIND=8) :: NORMX, NORMY
   END TYPE LINESOURCE

   TYPE(LINESOURCE), DIMENSION(:), ALLOCATABLE :: LINESOURCES

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! MPI parallelization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER      :: DOMPART_TYPE = -1  ! Type of domain partition (init to -1 to recognize)
   INTEGER      :: N_BLOCKS_X, N_BLOCKS_Y ! Used for "block" partitioning
   REAL(KIND=8) :: DX_BLOCKS, DY_BLOCKS   ! Used for "block" partitioning


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Multispecies !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CHARACTER(LEN=256) :: SPECIES_FILENAME
   INTEGER            :: N_SPECIES = 0
   
   TYPE SPECIES_DATA_STRUCTURE
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
      REAL(KIND=8) :: DIAM
      REAL(KIND=8) :: OMEGA
      REAL(KIND=8) :: TREF
      REAL(KIND=8) :: ALPHA
      REAL(KIND=8) :: SIGMA
      REAL(KIND=8) :: NU
      REAL(KIND=8) :: CREF
      
   END TYPE SPECIES_DATA_STRUCTURE

   TYPE(SPECIES_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: SPECIES, TEMP_SPECIES


   INTEGER            :: N_MIXTURES = 0

   TYPE MIXTURE_COMPONENT
      INTEGER :: ID
      CHARACTER*64 :: NAME
      REAL(KIND=8) :: MOLFRAC
   END TYPE MIXTURE_COMPONENT

   TYPE MIXTURE
      CHARACTER*64 :: NAME
      INTEGER      :: N_COMPONENTS
      TYPE(MIXTURE_COMPONENT), DIMENSION(:), ALLOCATABLE :: COMPONENTS
   END TYPE MIXTURE

   TYPE(MIXTURE), DIMENSION(:), ALLOCATABLE :: MIXTURES
 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Average flowfield !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER, DIMENSION(:), ALLOCATABLE      :: AVG_NP

   INTEGER, DIMENSION(:), ALLOCATABLE      :: AVG_N

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VX
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VY
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VZ

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRX
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRY
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRZ
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTR

   INTEGER                                 :: AVG_CUMULATED


   CONTAINS  ! @@@@@@@@@@@@@@@@@@@@@ SUBROUTINES @@@@@@@@@@@@@@@@@@@@@@@@

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE NEWTYPE -> defines a new type needed by MPI !!!!!!!!!!!!!
   ! to package messages in the "particle" format !!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE NEWTYPE
   
      INTEGER :: ii, extent_dpr, extent_int
      INTEGER, DIMENSION(10) :: blocklengths, oldtypes, offsets
     
      CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent_dpr, ierr)  
      CALL MPI_TYPE_EXTENT(MPI_INTEGER,          extent_int, ierr)
           
      blocklengths = 1
     
      oldtypes(1:8) = MPI_DOUBLE_PRECISION  
      oldtypes(9:10) = MPI_INTEGER
          
      offsets(1) = 0  
      DO ii = 2, 9
         offsets(ii) = offsets(ii - 1) + extent_dpr * blocklengths(ii - 1)
      END DO
      offsets(10) = offsets(9) + extent_int * blocklengths(9)
      
      CALL MPI_TYPE_STRUCT(10, blocklengths, offsets, oldtypes, MPI_PARTICLE_DATA_STRUCTURE, ierr)  
      CALL MPI_TYPE_COMMIT(MPI_PARTICLE_DATA_STRUCTURE, ierr)   
   
   END SUBROUTINE NEWTYPE


END MODULE global
