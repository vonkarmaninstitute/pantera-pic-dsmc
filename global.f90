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
   INTEGER      :: DIMS = 2
   INTEGER      :: NX, NY, NZ
   REAL(KIND=8) :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
   REAL(KIND=8) :: CELL_VOL
   LOGICAL      :: BOOL_X_PERIODIC = .FALSE. ! Assign default value!
   LOGICAL      :: BOOL_Y_PERIODIC = .FALSE.
   LOGICAL      :: BOOL_Z_PERIODIC = .FALSE.
   LOGICAL      :: BOOL_XMIN_SPECULAR = .FALSE.
   LOGICAL      :: BOOL_XMAX_SPECULAR = .FALSE.
   LOGICAL      :: BOOL_YMIN_SPECULAR = .FALSE.
   LOGICAL      :: BOOL_YMAX_SPECULAR = .FALSE.
   LOGICAL      :: BOOL_XMIN_REACT = .FALSE.
   LOGICAL      :: BOOL_XMAX_REACT = .FALSE.
   LOGICAL      :: BOOL_YMIN_REACT = .FALSE.
   LOGICAL      :: BOOL_YMAX_REACT = .FALSE.
   LOGICAL      :: BOOL_XMIN_DIFFUSE = .FALSE.
   LOGICAL      :: BOOL_XMAX_DIFFUSE = .FALSE.
   LOGICAL      :: BOOL_YMIN_DIFFUSE = .FALSE.
   LOGICAL      :: BOOL_YMAX_DIFFUSE = .FALSE.
   LOGICAL, DIMENSION(4) :: BOOL_PERIODIC = .FALSE.
   LOGICAL, DIMENSION(4) :: BOOL_SPECULAR = .FALSE.
   LOGICAL, DIMENSION(4) :: BOOL_REACT    = .FALSE.
   LOGICAL, DIMENSION(4) :: BOOL_DIFFUSE  = .FALSE.
   REAL(KIND=8) :: BOUNDTEMP
   CHARACTER(LEN=256) :: GRID_FILENAME
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: XCOORD, YCOORD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: XSIZE, YSIZE
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_VOLUMES
   
   ENUM, BIND(C)
      ENUMERATOR RECTILINEAR_UNIFORM, RECTILINEAR_NONUNIFORM, QUADTREE
   END ENUM
   INTEGER(KIND(RECTILINEAR_UNIFORM)) :: GRID_TYPE = RECTILINEAR_UNIFORM


   TYPE QUADTREE_CELL
      INTEGER :: LEVEL
      LOGICAL :: ISLEAF
      TYPE(QUADTREE_CELL), DIMENSION(:), ALLOCATABLE :: CHILDREN
      !TYPE(QUADTREE_CELL), POINTER :: PARENT
      REAL(KIND=8) :: XSIZE, YSIZE
   END TYPE QUADTREE_CELL

   TYPE(QUADTREE_CELL), DIMENSION(:), ALLOCATABLE :: QUADTREE_ROOT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Electromagnetic fields !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL :: BOOL_PIC = .FALSE.

   INTEGER :: NPX, NPY
   REAL(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: E_FIELD
   REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: PHI_FIELD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Q_FIELD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: DIRICHLET
   LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_DIRICHLET
   

   TYPE ST_MATRIX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: VALUE
      INTEGER, DIMENSION(:), ALLOCATABLE      :: RIDX, CIDX
      INTEGER :: NNZ
   END TYPE ST_MATRIX

   TYPE CC_MATRIX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AX
      INTEGER, DIMENSION(:), ALLOCATABLE      :: AP, AI
      INTEGER :: NNZ
   END TYPE CC_MATRIX

   TYPE(CC_MATRIX) :: A_CC

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Numerical settings !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER      :: NT, tID
   REAL(KIND=8) :: FNUM, DT, START_CPU_TIME
   INTEGER(KIND=8) :: RNG_SEED_GLOBAL, RNG_SEED_LOCAL
   INTEGER      :: DUMP_EVERY = 1
   INTEGER      :: DUMP_START = 0
   INTEGER      :: DUMP_GRID_AVG_EVERY = 1
   INTEGER      :: DUMP_GRID_START = 0
   INTEGER      :: DUMP_GRID_N_AVG = 1
   LOGICAL      :: PERFORM_CHECKS = .FALSE.
   INTEGER      :: CHECKS_EVERY = 1
   INTEGER      :: STATS_EVERY = 1

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Collisions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CHARACTER(LEN=64) :: COLLISION_TYPE
   LOGICAL           :: BOOL_MCC = .FALSE., BOOL_DSMC = .FALSE.
   REAL(KIND=8)      :: MCC_BG_DENS, MCC_BG_TTRA
   INTEGER           :: MCC_BG_MIX
   INTEGER           :: DSMC_COLL_MIX
   INTEGER           :: TIMESTEP_COLL
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_GREFS ! Matrix of reference relative velocities for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_SIGMAS ! Matrix of reference cross sections for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_ALPHAS ! Matrix of reference scattering coeff. for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_OMEGAS ! Matrix of reference temperature exponent for VSS
   REAL(KIND=8) :: SIGMAMAX

   LOGICAL           :: BOOL_THERMAL_BATH = .FALSE.
   REAL(KIND=8)      :: TBATH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Initial particles seed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL      :: BOOL_INITIAL_SEED = .FALSE. ! Assign default value!
   REAL(KIND=8) :: NRHO_INIT
   REAL(KIND=8) :: UX_INIT, UY_INIT, UZ_INIT
   REAL(KIND=8) :: TTRAX_INIT, TTRAY_INIT, TTRAZ_INIT, TROT_INIT, TVIB_INIT
   INTEGER      :: MIX_INIT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Particles injection from boundaries !!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL      :: BOOL_BOUNDINJECT = .FALSE. ! Assign default value!
   LOGICAL      :: BOOL_INJ_XMIN=.FALSE.,BOOL_INJ_XMAX=.FALSE. ! Assign default value!
   LOGICAL      :: BOOL_INJ_YMIN=.FALSE.,BOOL_INJ_YMAX=.FALSE. ! Assign default value!
   REAL(KIND=8) :: NRHO_BOUNDINJECT
   REAL(KIND=8) :: UX_BOUND, UY_BOUND, UZ_BOUND
   REAL(KIND=8) :: TTRA_BOUND, TROT_BOUND, TVIB_BOUND ! No TTRAX_ TTRAY_ TTRAZ_ for now! IMPLEMENT IT!
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
   REAL(KIND=8) :: TTRA_LINESOURCE, TROT_LINESOURCE, TVIB_LINESOURCE ! No TTRAX_ TTRAY_ TTRAZ_ for now! IMPLEMENT IT!

   INTEGER      :: nfs_LINESOURCE
   REAL(KIND=8) :: KAPPA_LINESOURCE
   REAL(KIND=8) :: S_NORM_LINESOURCE

   ! This is used instead.
   INTEGER         :: N_LINESOURCES = 0

   TYPE LINESOURCE
      REAL(KIND=8) :: CX, CY, DX, DY
      REAL(KIND=8) :: NRHO
      REAL(KIND=8) :: UX, UY, UZ
      REAL(KIND=8) :: TTRA, TROT, TVIB ! No TTRAX_ TTRAY_ TTRAZ_ for now! IMPLEMENT IT!
      INTEGER      :: MIX_ID

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: nfs
      REAL(KIND=8) :: KAPPA
      REAL(KIND=8) :: S_NORM
      REAL(KIND=8) :: NORMX, NORMY
   END TYPE LINESOURCE

   TYPE(LINESOURCE), DIMENSION(:), ALLOCATABLE :: LINESOURCES

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Walls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER         :: N_WALLS = 0

   TYPE WALL
      REAL(KIND=8) :: CX, CY, DX, DY
      REAL(KIND=8) :: TEMP
      LOGICAL      :: SPECULAR, DIFFUSE, POROUS, REACT
      REAL(KIND=8) :: TRANSMISSIVITY
      REAL(KIND=8) :: NORMX, NORMY
   END TYPE WALL

   TYPE(WALL), DIMENSION(:), ALLOCATABLE :: WALLS

   
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
      REAL(KIND=8) :: MOLECULAR_MASS
      INTEGER      :: ROTDOF
      REAL(KIND=8) :: ROTREL
      INTEGER      :: VIBDOF
      REAL(KIND=8) :: VIBREL
      REAL(KIND=8) :: VIBTEMP
      REAL(KIND=8) :: SPWT
      REAL(KIND=8) :: INVSPWT
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
   !!!!!!!!! Chemical reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER(LEN=256) :: REACTIONS_FILENAME
   INTEGER :: N_REACTIONS = 0
   
   TYPE REACTIONS_DATA_STRUCTURE
      INTEGER :: R1_SP_ID
      INTEGER :: R2_SP_ID
      INTEGER :: P1_SP_ID
      INTEGER :: P2_SP_ID
      INTEGER :: P3_SP_ID
      REAL(KIND=8) :: A, N, EA
      REAL(KIND=8) :: C1, C2, C3
      INTEGER :: N_PROD
      LOGICAL :: IS_CEX
   END TYPE REACTIONS_DATA_STRUCTURE

   TYPE(REACTIONS_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: REACTIONS, TEMP_REACTIONS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Wall reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER(LEN=256) :: WALL_REACTIONS_FILENAME
   INTEGER :: N_WALL_REACTIONS = 0
   
   TYPE WALL_REACTIONS_DATA_STRUCTURE
      INTEGER :: R_SP_ID
      INTEGER :: P1_SP_ID
      INTEGER :: P2_SP_ID
      REAL(KIND=8) :: PROB
      INTEGER :: N_PROD
   END TYPE WALL_REACTIONS_DATA_STRUCTURE

   TYPE(WALL_REACTIONS_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: WALL_REACTIONS, TEMP_WALL_REACTIONS


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Average flowfield !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_N

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_NP

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VX
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VY
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VZ

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRX
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRY
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRZ
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTR

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TROT
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TVIB

   INTEGER                                 :: AVG_CUMULATED
   INTEGER, DIMENSION(:), ALLOCATABLE      :: AVG_CUMULATED_INTENSIVE_ONE
   INTEGER, DIMENSION(:), ALLOCATABLE      :: AVG_CUMULATED_INTENSIVE_TWO


   CONTAINS  ! @@@@@@@@@@@@@@@@@@@@@ SUBROUTINES @@@@@@@@@@@@@@@@@@@@@@@@

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE NEWTYPE -> defines a new type needed by MPI !!!!!!!!!!!!!
   ! to package messages in the "particle" format !!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE NEWTYPE
   
      INTEGER :: ii, extent_dpr, extent_int
      INTEGER, DIMENSION(11) :: blocklengths, oldtypes, offsets
     
      CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent_dpr, ierr)  
      CALL MPI_TYPE_EXTENT(MPI_INTEGER,          extent_int, ierr)
           
      blocklengths = 1
     
      oldtypes(1:9) = MPI_DOUBLE_PRECISION  
      oldtypes(10:11) = MPI_INTEGER
          
      offsets(1) = 0  
      DO ii = 2, 10
         offsets(ii) = offsets(ii - 1) + extent_dpr * blocklengths(ii - 1)
      END DO
      offsets(11) = offsets(10) + extent_int * blocklengths(10)
      
      CALL MPI_TYPE_STRUCT(11, blocklengths, offsets, oldtypes, MPI_PARTICLE_DATA_STRUCTURE, ierr)  
      CALL MPI_TYPE_COMMIT(MPI_PARTICLE_DATA_STRUCTURE, ierr)   
   
   END SUBROUTINE NEWTYPE


   SUBROUTINE ST_MATRIX_ALLOCATE(THIS, MAXDIM)

      TYPE(ST_MATRIX), INTENT(INOUT) :: THIS
      INTEGER, INTENT(IN) :: MAXDIM

      ALLOCATE(THIS%RIDX(0:MAXDIM-1))
      ALLOCATE(THIS%CIDX(0:MAXDIM-1))
      ALLOCATE(THIS%VALUE(0:MAXDIM-1))
      THIS%NNZ = 0

   END SUBROUTINE ST_MATRIX_ALLOCATE


   SUBROUTINE ST_MATRIX_SET(THIS, ROWIDX, COLIDX, VAL)

      TYPE(ST_MATRIX), INTENT(INOUT) :: THIS
      INTEGER, INTENT(IN) :: ROWIDX, COLIDX
      REAL(KIND=8), INTENT(IN) :: VAL

      THIS%RIDX(THIS%NNZ) = ROWIDX
      THIS%CIDX(THIS%NNZ) = COLIDX
      THIS%VALUE(THIS%NNZ) = VAL

      THIS%NNZ = THIS%NNZ + 1

   END SUBROUTINE ST_MATRIX_SET


   SUBROUTINE ST_MATRIX_TO_CC(THIS, NCOLS, CC)

      IMPLICIT NONE

      TYPE(ST_MATRIX), INTENT(IN) :: THIS
      TYPE(CC_MATRIX), INTENT(OUT) :: CC
      INTEGER, INTENT(IN) :: NCOLS
      INTEGER :: NNZ, I, J, K, L
      INTEGER :: AI_TEMP
      REAL(KIND=8) :: AX_TEMP
      INTEGER, DIMENSION(:), ALLOCATABLE :: NPC

      NNZ = THIS%NNZ
      CC%NNZ = NNZ
      ALLOCATE(CC%AP(0:NCOLS))
      ALLOCATE(NPC(0:NCOLS))
      ALLOCATE(CC%AI(0:NNZ-1))
      ALLOCATE(CC%AX(0:NNZ-1))

      ! Compute the offsets to the first element of a column
      ! AP(I) is the offset of column I in rows AI with values AX
      ! i.e. the matrix entry is a(AI(J), I) = AX(J), with J = AP(I)+K, J<AP(I+1)
      CC%AP = 0
      DO I = 0, NNZ-1
         CC%AP( THIS%CIDX(I)+1 ) = CC%AP( THIS%CIDX(I)+1 ) + 1
      END DO
      DO I = 1, NCOLS-1
         CC%AP(I+1) =  CC%AP(I+1) + CC%AP(I)
      END DO
      IF (CC%AP(NCOLS) .NE. NNZ) WRITE(*,*) 'I believe we ve had a problem here.'

      ! Now populate AI and AX at the right offsets (but not sorted!)
      NPC = 0
      DO I = 0, NNZ-1
         J = THIS%CIDX(I)
         CC%AI(CC%AP(J) + NPC(J)) = THIS%RIDX(I)
         CC%AX(CC%AP(J) + NPC(J)) = THIS%VALUE(I) 
         NPC(J) = NPC(J) + 1
      END DO

      ! Now sort AI and AX by increasing AI in each column
      DO J = 0, NCOLS-1
         ! A classic sort routine, within the correct bounds
         DO K = CC%AP(J), CC%AP(J+1)-2
            DO L = K+1, CC%AP(J+1)-1
               IF (CC%AI(K) .GT. CC%AI(L)) THEN
                  AI_TEMP = CC%AI(K)
                  AX_TEMP = CC%AX(K)
                  
                  CC%AI(K) = CC%AI(L)
                  CC%AX(K) = CC%AX(L)

                  CC%AI(L) = AI_TEMP
                  CC%AX(L) = AX_TEMP
               END IF
            END DO
         END DO

      END DO


   END SUBROUTINE ST_MATRIX_TO_CC

END MODULE global
