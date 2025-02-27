! Copyright (C) 2024 von Karman Institute for Fluid Dynamics (VKI)
!
! This file is part of PANTERA PIC-DSMC, a software for the simulation
! of rarefied gases and plasmas using particles.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.PANTERA PIC-DSMC

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

      CHARACTER*512      :: MIXTURE_DEFINITION, VSS_PARAMS_FILENAME, LINESOURCE_DEFINITION, WALL_DEFINITION, MCC_BG_FILENAME
      CHARACTER*512      :: BC_DEFINITION, SOLENOID_DEFINITION
      CHARACTER*64       :: MIX_BOUNDINJECT_NAME, DSMC_COLL_MIX_NAME, MCC_BG_MIX_NAME, PIC_TYPE_STRING, PARTITION_STYLE_STRING, &
      COLLISION_TYPE_STRING, REMOVE_MIX_NAME

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
         IF (line=='Restart:')                 READ(in1,*) RESTART_TIMESTEP

         IF (line=='Axisymmetric:')            READ(in1,*) AXI
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

         IF (line=='Grid_file:') THEN
            READ(in1,*) GRID_FILENAME
            CALL READ_GRID_FILE(GRID_FILENAME)
            GRID_TYPE = RECTILINEAR_NONUNIFORM
         END IF

         IF (line=='Mesh_file_SU2:') THEN
            READ(in1,*) MESH_FILENAME
            IF (DIMS == 1) THEN
               CALL READ_1D_UNSTRUCTURED_GRID_SU2(MESH_FILENAME)
            ELSE IF (DIMS == 2) THEN
               CALL READ_2D_UNSTRUCTURED_GRID_SU2(MESH_FILENAME)
            ELSE IF (DIMS == 3) THEN
               CALL READ_3D_UNSTRUCTURED_GRID_SU2(MESH_FILENAME)
            ELSE
               CALL ERROR_ABORT('Asked to read SU2 grid file but dimensions are not set to 2 or 3.')
            END IF
            GRID_TYPE = UNSTRUCTURED
         END IF

         IF (line=='Boundary_condition:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_BOUNDARY_CONDITION(BC_DEFINITION)
         END IF

         IF (line=='Domain_type:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_DOMAIN_TYPE(BC_DEFINITION)
         END IF

         IF (line=='Boundary_emit:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_BOUNDARY_EMIT(BC_DEFINITION)
         END IF

         IF (line=='Solenoid_field:') THEN
            READ(in1,'(A)') SOLENOID_DEFINITION
            CALL DEF_SOLENOID(SOLENOID_DEFINITION)
         END IF

         IF (line=='External_B_field:') READ(in1,*) EXTERNAL_B_FIELD(1), EXTERNAL_B_FIELD(2), EXTERNAL_B_FIELD(3)

         IF (line=='Magnetic_dipole:') THEN
            BOOL_MAGNETIC_DIPOLE = .TRUE.
            READ(in1,*) DIPOLE_POSITION(:), DIPOLE_ORIENTATION(:), MAGNETIC_MOMENT
         END IF

         IF (line=='Boundary_dump_fluxes:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_BOUNDARY_DUMP_FLUXES(BC_DEFINITION)
         END IF


         ! ~~~~~~~~~~~~~  Numerical settings  ~~~~~~~~~~~~~~~~~
         IF (line=='Fnum:')                    READ(in1,*) FNUM
         IF (line=='Bool_radial_weighting:')   READ(in1,*) BOOL_RADIAL_WEIGHTING
         IF (line=='Timestep:')                READ(in1,*) DT
         IF (line=='Number_of_timesteps:')     READ(in1,*) NT
         IF (line=='RNG_seed:')                READ(in1,*) RNG_SEED_GLOBAL
         IF (line=='Perform_checks:')          READ(in1,*) PERFORM_CHECKS
         IF (line=='Checks_every:')            READ(in1,*) CHECKS_EVERY
         IF (line=='Stats_every:')             READ(in1,*) STATS_EVERY
         IF (line=='Epsilon_scaling:')         READ(in1,*) EPS_SCALING
         IF (line=='PIC_type:') THEN
            READ(in1,*) PIC_TYPE_STRING
            IF (PIC_TYPE_STRING == "none") THEN
               PIC_TYPE = NONE
            ELSE IF (PIC_TYPE_STRING == "explicit") THEN
               PIC_TYPE = EXPLICIT
            ELSE IF (PIC_TYPE_STRING == "fullyimplicit") THEN
               PIC_TYPE = FULLYIMPLICIT
            ELSE IF (PIC_TYPE_STRING == "semiimplicit") THEN
               PIC_TYPE = SEMIIMPLICIT
            ELSE IF (PIC_TYPE_STRING == "explicitlimited") THEN
               PIC_TYPE = EXPLICITLIMITED
            ELSE IF (PIC_TYPE_STRING == "hybrid") THEN
               PIC_TYPE = HYBRID
            ELSE
               CALL ERROR_ABORT('PIC type is not specified correctly in input script.')
            END IF
         END IF
         
         ! Fluid electrons setting
         IF (line=='Fluid_electrons_n0:')      READ(in1,*) BOLTZ_N0
         IF (line=='Fluid_electrons_phi0:')    READ(in1,*) BOLTZ_PHI0
         IF (line=='Fluid_electrons_Te:')      READ(in1,*) BOLTZ_TE

         IF (line=='Bool_kappa_fluid:')        READ(in1,*) BOOL_KAPPA_FLUID
         IF (line=='Kappa_constant_fluid:')    READ(in1,*) KAPPA_FLUID_C

         IF (line=='Jacobian_type:')           READ(in1,*) JACOBIAN_TYPE
         IF (line=='Residual_and_jacobian_combined:') READ(in1,*) RESIDUAL_AND_JACOBIAN_COMBINED
         IF (line=='Colocated_electrons:') THEN
            READ(in1,*) COLOCATED_ELECTRONS_TTRA
            COLOCATED_ELECTRONS = .TRUE.
         END IF
         IF (line=='SNES_rtol:')               READ(in1,*) SNES_RTOL
         
         ! ~~~~~~~~~~~~~  File output ~~~~~~~~~~~~~~~

         IF (line=='Flowfield_output:')        READ(in1,*) FLOWFIELD_SAVE_PATH
         IF (line=='Boundary_output:')         READ(in1,*) BOUNDARY_SAVE_PATH
         IF (line=='Particle_dump_output:')    READ(in1,*) PARTDUMP_SAVE_PATH
         IF (line=='Trajectory_dump_output:')  READ(in1,*) TRAJDUMP_SAVE_PATH
         IF (line=='Fluxes_dump_output:')      READ(in1,*) FLUXDUMP_SAVE_PATH
         IF (line=='Checks_output:')           READ(in1,*) CHECKS_SAVE_PATH
         IF (line=='All_output_path:') THEN
            READ(in1,*) FLOWFIELD_SAVE_PATH
            CHECKS_SAVE_PATH = FLOWFIELD_SAVE_PATH
            FLUXDUMP_SAVE_PATH = FLOWFIELD_SAVE_PATH
            TRAJDUMP_SAVE_PATH = FLOWFIELD_SAVE_PATH
            PARTDUMP_SAVE_PATH = FLOWFIELD_SAVE_PATH
            RESIDUAL_SAVE_PATH = FLOWFIELD_SAVE_PATH
            BOUNDARY_SAVE_PATH = FLOWFIELD_SAVE_PATH
         END IF
         IF (line=='Binary_output:')           READ(in1,*) BOOL_BINARY_OUTPUT
         IF (line=='Dump_part_every:')         READ(in1,*) DUMP_PART_EVERY
         IF (line=='Dump_part_start:')         READ(in1,*) DUMP_PART_START
         IF (line=='Dump_part_fracsample:')    READ(in1,*) PARTDUMP_FRACSAMPLE
         IF (line=='Dump_part_bound_every:')   READ(in1,*) DUMP_PART_BOUND_EVERY
         IF (line=='Dump_part_bound_start:')   READ(in1,*) DUMP_PART_BOUND_START
         IF (line=='Dump_bound_avgevery:')     READ(in1,*) DUMP_BOUND_AVG_EVERY
         IF (line=='Dump_bound_start:')        READ(in1,*) DUMP_BOUND_START
         IF (line=='Dump_bound_numavgs:')      READ(in1,*) DUMP_BOUND_N_AVG
         IF (line=='Dump_grid_avgevery:')      READ(in1,*) DUMP_GRID_AVG_EVERY
         IF (line=='Dump_grid_start:')         READ(in1,*) DUMP_GRID_START
         IF (line=='Dump_grid_numavgs:')       READ(in1,*) DUMP_GRID_N_AVG
         IF (line=='Bool_dump_moments:')       READ(in1,*) BOOL_DUMP_MOMENTS
         IF (line=='Bool_dump_fluxes:')        READ(in1,*) BOOL_DUMP_FLUXES
         IF (line=='Dump_traj_start:')         READ(in1,*) DUMP_TRAJECTORY_START
         IF (line=='Dump_traj_number:')        READ(in1,*) DUMP_TRAJECTORY_NUMBER


         IF (line=='Inject_from_file:') THEN
            BOOL_INJECT_FROM_FILE = .TRUE.
            READ(in1,*) INJECT_FILENAME
            CALL READ_INJECT_FILE
         END IF

         IF (line=='Inject_probability:')      READ(in1,*) INJECT_PROBABILITY

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
         IF (line=='Collision_type:') THEN
            READ(in1,*) COLLISION_TYPE_STRING
            IF (COLLISION_TYPE_STRING == "NONE") THEN
               COLLISION_TYPE = NO_COLL
            ELSE IF (COLLISION_TYPE_STRING == "DSMC") THEN
               COLLISION_TYPE = DSMC
            ELSE IF (COLLISION_TYPE_STRING == "MCC") THEN
               COLLISION_TYPE = MCC
            ELSE IF (COLLISION_TYPE_STRING == "VAHEDI_MCC") THEN
               COLLISION_TYPE = MCC_VAHEDI
            ELSE IF (COLLISION_TYPE_STRING == "VAHEDI_DSMC") THEN
               COLLISION_TYPE = DSMC_VAHEDI
            ELSE IF (COLLISION_TYPE_STRING == "BGK") THEN
               COLLISION_TYPE = BGK
            ELSE
               CALL ERROR_ABORT('Collision type is not specified correctly in input script.')
            END IF
         END IF
         IF (line=='MCC_background_dens:')     READ(in1,*) MCC_BG_DENS
         IF (line=='MCC_background_Ttra:')     READ(in1,*) MCC_BG_TTRA
         IF (line=='MCC_background_mixture:') THEN
            READ(in1,*) MCC_BG_MIX_NAME
            MCC_BG_MIX = MIXTURE_NAME_TO_ID(MCC_BG_MIX_NAME)
         END IF
         IF (line=='MCC_background_file:') THEN
            READ(in1,*) MCC_BG_FILENAME
            CALL READ_MCC_BACKGROUND_FILE(MCC_BG_FILENAME)
         END IF
         IF (line=='VSS_parameters_file:')     THEN
            READ(in1,*) VSS_PARAMS_FILENAME
            CALL READ_VSS(VSS_PARAMS_FILENAME)
         END IF
         IF (line=='VSS_parameters_binary_file:')     THEN
            READ(in1,*) VSS_PARAMS_FILENAME
            CALL READ_VSS_BINARY(VSS_PARAMS_FILENAME)
         END IF
         IF (line=='DSMC_collisions_mixture:') THEN
            READ(in1,*) DSMC_COLL_MIX_NAME
            DSMC_COLL_MIX = MIXTURE_NAME_TO_ID(DSMC_COLL_MIX_NAME)
         END IF
         IF (line=='BGK_sigma:')     READ(in1,*) BGK_SIGMA
         
        ! ~~~~~~~~~~~~~  Reactions ~~~~~~~~~~~~~~~
         IF (line=='Reactions_file:') THEN
            READ(in1,*) REACTIONS_FILENAME
            CALL READ_REACTIONS(REACTIONS_FILENAME)
         END IF

         ! ~~~~~~~~~~~~~  Initial particles seeding  ~~~~~~~~~~~~~~~~~
         IF (line=='Initial_particles:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_INITIAL_PARTICLES(BC_DEFINITION)
         END IF

         ! ~~~~~~~~~~~~~  Continuous injection in the volume  ~~~~~~~~~~~~~~~~~
         IF (line=='Volume_inject:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_VOLUME_INJECT(BC_DEFINITION)
         END IF

         ! ~~~~~~~~~~~~~  Particle deletion  ~~~~~~~~~~~~~~~~~
         IF (line=='Remove_particles_in_mixture:') THEN
            READ(in1,*) REMOVE_MIX_NAME
            REMOVE_MIX = MIXTURE_NAME_TO_ID(REMOVE_MIX_NAME)
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
         IF (line=='Partition_style:') THEN
            READ(in1,*) PARTITION_STYLE_STRING
            IF (PARTITION_STYLE_STRING == "stripsx") THEN
               PARTITION_STYLE = STRIPSX
            ELSE IF (PARTITION_STYLE_STRING == "stripsy") THEN
               PARTITION_STYLE = STRIPSY
            ELSE IF (PARTITION_STYLE_STRING == "stripsz") THEN
               PARTITION_STYLE = STRIPSZ
            ELSE IF (PARTITION_STYLE_STRING == "slicesx") THEN
               PARTITION_STYLE = SLICESX
            ELSE IF (PARTITION_STYLE_STRING == "slicesz") THEN
               PARTITION_STYLE = SLICESZ
            ELSE
               CALL ERROR_ABORT('Partition style is not specified correctly in input script.')
            END IF
         END IF

         IF (line=='Load_balance:')    READ(in1,*) LOAD_BALANCE
         IF (line=='Load_balance_every:')    READ(in1,*) LOAD_BALANCE_EVERY
         
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
         WRITE(*,'(A5, A50,L)') '     ', string, AXI
 
         string = 'Domain limits:'
         WRITE(*,'(A5,A50)') '    ', string
         WRITE(*,'(A13,ES14.3,ES14.3)')          '      X [m]: ', XMIN, XMAX
         WRITE(*,'(A13,ES14.3,ES14.3)')          '      Y [m]: ', YMIN, YMAX
         WRITE(*,'(A13,ES14.3,ES14.3)')          '      Z [m]: ', ZMIN, ZMAX

         IF (GRID_TYPE == UNSTRUCTURED) THEN
            string = 'Unstructured grid with number of cells, nodes: '
            WRITE(*,'(A5,A50,I8,I8)')'    ',string,NCELLS,NNODES
         ELSE
            string = 'Structured grid with number of cells: '
            WRITE(*,'(A5,A50, I8, I8, I8)') '    ', string, NX, NY, NZ

            string = 'Domain periodicity along X, Y, Z [T/F]: '
            WRITE(*,'(A5,A50,L4,L4,L4)')'    ',string,BOOL_X_PERIODIC,BOOL_Y_PERIODIC,BOOL_Z_PERIODIC
         END IF
  
         ! ~~~~ Numerical settings ~~~~
         WRITE(*,*) '  =========== Numerical settings ==================================='
         string = 'Fnum:'
         WRITE(*,'(A5,A50,ES14.3)') '     ', string, FNUM

         string = 'Timestep [s]:'
         WRITE(*,'(A5,A50,ES14.3)') '     ', string, DT
 
         string = 'Initial timestep:'
         WRITE(*,'(A5,A50,I9)') '     ',     string, RESTART_TIMESTEP

         string = 'Final timestep:'
         WRITE(*,'(A5,A50,I9)') '     ',     string, NT

         string = 'RNG seed (global):'
         WRITE(*,'(A5,A50,I9)') '     ', string, RNG_SEED_GLOBAL

         ! ~~~~ Collisions ~~~~
         WRITE(*,*) '  =========== Collisions ============================'
         string = 'Collision type:'
         WRITE(*,'(A5,A50,I8)') '     ', string, COLLISION_TYPE

         IF (COLLISION_TYPE == MCC) THEN  ! Only print this if MCC collisions are ON

            string = 'Background number density [1/m^3]:'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, MCC_BG_DENS

            string = 'Background translational temperature [K]::'
            WRITE(*,'(A5,A50,ES14.3)') '     ', string, MCC_BG_TTRA

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

         ! ~~~~ Multispecies ~~~~
         WRITE(*,*) '  =========== Species ========================='
         WRITE(*,'(A5,A30,I9)') '    ','Number of species present: ', N_SPECIES


         WRITE(*,*) '  =========== Mixtures ========================='
         string = 'Number of mixtures:'
         WRITE(*,'(A5,A50,I9)') '    ', string, N_MIXTURES
         DO j = 1, N_MIXTURES

            WRITE(*,*)'    ','Mixture ', j, ' named ', TRIM(MIXTURES(j)%NAME), ' has ', MIXTURES(j)%N_COMPONENTS, ' components:'
            DO k = 1, MIXTURES(j)%N_COMPONENTS
               WRITE(*,*) '    ','Mixture component ', k, ' is ', TRIM(MIXTURES(j)%COMPONENTS(k)%NAME), &
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
      
         IF (ALLOCATED(SPECIES)) THEN
            ALLOCATE(TEMP_SPECIES(N_SPECIES+1)) ! Append the species to the list
            TEMP_SPECIES(1:N_SPECIES) = SPECIES(1:N_SPECIES)
            CALL MOVE_ALLOC(TEMP_SPECIES, SPECIES)
         ELSE
            ALLOCATE(SPECIES(1))
         END IF
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

      INTEGER      :: IS, JS
      REAL(KIND=8) :: M1, M2, MRED

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


         IF (ABS(OMEGA-0.5) .LT. 1.d-6) THEN
            SPECIES(SP_ID)%CREF  = (4.*KB*TREF/SPECIES(SP_ID)%MOLECULAR_MASS)**0.5 * 1.2354
         ELSE
            SPECIES(SP_ID)%CREF  = (4.*KB*TREF/SPECIES(SP_ID)%MOLECULAR_MASS)**0.5 * (GAMMA(2.5-OMEGA))**(-0.5/(OMEGA-0.5))
         END IF
         
         
         
         
      END DO

      ALLOCATE(VSS_GREFS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_SIGMAS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_ALPHAS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_OMEGAS(N_SPECIES, N_SPECIES))


      SIGMAMAX = 0
      DO IS = 1, N_SPECIES
         DO JS = 1, N_SPECIES
            M1    = SPECIES(IS)%MOLECULAR_MASS
            M2    = SPECIES(JS)%MOLECULAR_MASS
            MRED  = M1*M2/(M1+M2)
            OMEGA    = 0.5 * (SPECIES(IS)%OMEGA + SPECIES(JS)%OMEGA)
            TREF = 0.5 * (SPECIES(IS)%TREF + SPECIES(JS)%TREF)
            
            IF (ABS(OMEGA-0.5) .LT. 1.d-6) THEN
               VSS_GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * 1.2354
            ELSE
               VSS_GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * (GAMMA(2.5-OMEGA))**(-0.5/(OMEGA-0.5))
            END IF

            VSS_SIGMAS(IS, JS) = PI * (0.5 * (SPECIES(IS)%DIAM + SPECIES(JS)%DIAM))**2
            VSS_ALPHAS(IS, JS) = 0.5 * (SPECIES(IS)%ALPHA + SPECIES(JS)%ALPHA)
            VSS_OMEGAS(IS, JS) = 0.5 * (SPECIES(IS)%OMEGA + SPECIES(JS)%OMEGA)

            IF (VSS_SIGMAS(IS, JS) .GT. SIGMAMAX) SIGMAMAX = VSS_SIGMAS(IS, JS)
         END DO
      END DO

      CLOSE(in2) ! Close input file


   END SUBROUTINE READ_VSS




   SUBROUTINE READ_VSS_BINARY(FILENAME)

      IMPLICIT NONE

      CHARACTER*64, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in2 = 457
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      INTEGER, ALLOCATABLE :: SP_IDS(:)
      REAL(KIND=8) :: OMEGA
      REAL(KIND=8) :: TREF
      REAL(KIND=8) :: READ_VALUE

      INTEGER      :: IS, JS
      REAL(KIND=8) :: M1, M2, MRED

      ALLOCATE(VSS_GREFS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_SIGMAS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_ALPHAS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_OMEGAS(N_SPECIES, N_SPECIES))

      ! Open input file for reading
      OPEN(UNIT=in2,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, VSS parameters definition file not found! ABORTING.')
      ENDIF

      ! Read Tref
      line = '' ! Init empty
      DO WHILE (line == '')
         
         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) THEN
            CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
         END IF
         
      END DO

      READ(line,*) TREF
      
      ! Read species ordering
      N_STR = 0
      DO WHILE (N_STR == 0)
         line = '' ! Init empty
         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) THEN
            CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
         END IF

         CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
      END DO

      IF (N_STR .NE. N_SPECIES) CALL ERROR_ABORT('Attention, incorrect number of species in VSS definition! ABORTING.')

      ALLOCATE(SP_IDS(N_SPECIES))
      DO IS = 1, N_SPECIES
         SP_IDS(IS) = SPECIES_NAME_TO_ID(STRARRAY(IS))
      END DO


      ! Read molecular diameters
      DO JS = 1, N_SPECIES
         N_STR = 0
         DO WHILE (N_STR == 0)
            line = '' ! Init empty
            READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
            CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

            IF (ReasonEOF < 0) THEN
               CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
            END IF

            CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
         END DO
         
         IF (N_STR .NE. N_SPECIES+1-JS) CALL ERROR_ABORT('Attention, VSS params matrix format wrong! ABORTING.')
         

         DO IS = JS, N_SPECIES
            READ(STRARRAY(1+IS-JS),*) READ_VALUE
            VSS_SIGMAS(SP_IDS(IS), SP_IDS(JS)) = PI*READ_VALUE**2
            VSS_SIGMAS(SP_IDS(JS), SP_IDS(IS)) = VSS_SIGMAS(SP_IDS(IS), SP_IDS(JS))
         END DO
      END DO
         

      ! Read omegas
      DO JS = 1, N_SPECIES
         N_STR = 0
         DO WHILE (N_STR == 0)
            line = '' ! Init empty
            READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
            CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

            IF (ReasonEOF < 0) THEN
               CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
            END IF

            CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
         END DO

         DO IS = JS, N_SPECIES
            READ(STRARRAY(1+IS-JS),*) READ_VALUE
            VSS_OMEGAS(SP_IDS(IS), SP_IDS(JS)) = READ_VALUE
            VSS_OMEGAS(SP_IDS(JS), SP_IDS(IS)) = VSS_OMEGAS(SP_IDS(IS), SP_IDS(JS))
         END DO
      END DO


      ! Read alphas
      DO JS = 1, N_SPECIES
         N_STR = 0
         DO WHILE (N_STR == 0)
            line = '' ! Init empty
            READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
            CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

            IF (ReasonEOF < 0) THEN
               CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
            END IF

            CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
         END DO

         DO IS = JS, N_SPECIES
            READ(STRARRAY(1+IS-JS),*) READ_VALUE
            VSS_ALPHAS(SP_IDS(IS), SP_IDS(JS)) = READ_VALUE
            VSS_ALPHAS(SP_IDS(JS), SP_IDS(IS)) = VSS_ALPHAS(SP_IDS(IS), SP_IDS(JS))
         END DO
      END DO
 

      SIGMAMAX = 0
      DO IS = 1, N_SPECIES
         DO JS = 1, N_SPECIES
            M1    = SPECIES(IS)%MOLECULAR_MASS
            M2    = SPECIES(JS)%MOLECULAR_MASS
            MRED  = M1*M2/(M1+M2)
            OMEGA = VSS_OMEGAS(IS, JS)
                        
            IF (ABS(OMEGA-0.5) .LT. 1.d-6) THEN
               VSS_GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * 1.2354
            ELSE
               VSS_GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * (GAMMA(2.5-OMEGA))**(-0.5/(OMEGA-0.5))
            END IF

            IF (VSS_SIGMAS(IS, JS) .GT. SIGMAMAX) SIGMAMAX = VSS_SIGMAS(IS, JS)
         END DO
      END DO

      CLOSE(in2) ! Close input file

   

      ! WRITE(*,*) 'Tref = ', TREF
      ! WRITE(*,*) 'sigma = '
      ! DO IS = 1, N_SPECIES
      !    WRITE(*,*) VSS_SIGMAS(IS,:)
      ! END DO
      ! WRITE(*,*) 'omega = '
      ! DO IS = 1, N_SPECIES
      !    WRITE(*,*) VSS_OMEGAS(IS,:)
      ! END DO
      ! WRITE(*,*) 'alpha = '
      ! DO IS = 1, N_SPECIES
      !    WRITE(*,*) VSS_ALPHAS(IS,:)
      ! END DO
      ! WRITE(*,*) 'gref = '
      ! DO IS = 1, N_SPECIES
      !    WRITE(*,*) VSS_GREFS(IS,:)
      ! END DO

   END SUBROUTINE READ_VSS_BINARY


   SUBROUTINE DEF_MIXTURE(DEFINITION)
      
      IMPLICIT NONE

      CHARACTER*64 :: MIX_NAME
      CHARACTER*64 :: COMP_NAME
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


      IF (ALLOCATED(MIXTURES)) THEN
         ALLOCATE(TEMP_MIXTURES(N_MIXTURES+1)) ! Append the mixture to the list
         TEMP_MIXTURES(1:N_MIXTURES) = MIXTURES(1:N_MIXTURES)
         CALL MOVE_ALLOC(TEMP_MIXTURES, MIXTURES)
      ELSE
         ALLOCATE(MIXTURES(1))
      END IF
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

      IF (ALLOCATED(LINESOURCES)) THEN
         ALLOCATE(TEMP_LINESOURCES(N_LINESOURCES+1)) ! Append the mixture to the list
         TEMP_LINESOURCES(1:N_LINESOURCES) = LINESOURCES(1:N_LINESOURCES)
         CALL MOVE_ALLOC(TEMP_LINESOURCES, LINESOURCES)
      ELSE
         ALLOCATE(LINESOURCES(1))
      END IF
      N_LINESOURCES = N_LINESOURCES + 1



      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      READ(STRARRAY(1), '(ES14.0)') LINESOURCES(N_LINESOURCES)%X1
      READ(STRARRAY(2), '(ES14.0)') LINESOURCES(N_LINESOURCES)%Y1
      READ(STRARRAY(3), '(ES14.0)') LINESOURCES(N_LINESOURCES)%X2
      READ(STRARRAY(4), '(ES14.0)') LINESOURCES(N_LINESOURCES)%Y2
      READ(STRARRAY(5), '(ES14.0)') LINESOURCES(N_LINESOURCES)%NRHO
      READ(STRARRAY(6), '(ES14.0)') LINESOURCES(N_LINESOURCES)%UX
      READ(STRARRAY(7), '(ES14.0)') LINESOURCES(N_LINESOURCES)%UY
      READ(STRARRAY(8), '(ES14.0)') LINESOURCES(N_LINESOURCES)%UZ
      READ(STRARRAY(9), '(ES14.0)') LINESOURCES(N_LINESOURCES)%TTRA
      READ(STRARRAY(10),'(ES14.0)') LINESOURCES(N_LINESOURCES)%TROT
      READ(STRARRAY(10),'(ES14.0)') LINESOURCES(N_LINESOURCES)%TVIB
      READ(STRARRAY(11),'(A10)') MIX_NAME
      ! WRITE(*,*) LINESOURCES(N_LINESOURCES)%UX
      ! WRITE(*,*) LINESOURCES(N_LINESOURCES)%NRHO
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

      IF (ALLOCATED(WALLS)) THEN
         ALLOCATE(TEMP_WALLS(N_WALLS+1))
         TEMP_WALLS(1:N_WALLS) = WALLS(1:N_WALLS)
         CALL MOVE_ALLOC(TEMP_WALLS, WALLS)
      ELSE
         ALLOCATE(WALLS(1))
      END IF
      N_WALLS = N_WALLS + 1



      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      READ(STRARRAY(1), '(ES14.0)') WALLS(N_WALLS)%X1
      READ(STRARRAY(2), '(ES14.0)') WALLS(N_WALLS)%Y1
      READ(STRARRAY(3), '(ES14.0)') WALLS(N_WALLS)%X2
      READ(STRARRAY(4), '(ES14.0)') WALLS(N_WALLS)%Y2
      READ(STRARRAY(5), '(ES14.0)') WALLS(N_WALLS)%TEMP
      READ(STRARRAY(6), *) WALLS(N_WALLS)%SPECULAR
      READ(STRARRAY(7), *) WALLS(N_WALLS)%DIFFUSE
      READ(STRARRAY(8), *) WALLS(N_WALLS)%POROUS
      READ(STRARRAY(9), '(ES14.0)') WALLS(N_WALLS)%TRANSMISSIVITY
      READ(STRARRAY(10), *) WALLS(N_WALLS)%REACT

      LENGTH = SQRT((WALLS(N_WALLS)%X2 - WALLS(N_WALLS)%X1)**2 + (WALLS(N_WALLS)%Y2 - WALLS(N_WALLS)%Y1)**2)

      WALLS(N_WALLS)%NORMX = (WALLS(N_WALLS)%Y2 - WALLS(N_WALLS)%Y1)/LENGTH
      WALLS(N_WALLS)%NORMY = - (WALLS(N_WALLS)%X2 - WALLS(N_WALLS)%X1)/LENGTH

   END SUBROUTINE DEF_WALL


   SUBROUTINE DEF_BOUNDARY_CONDITION(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR, I, IPG
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)


      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      ! phys_group type parameters
      IPG = -1
      DO I = 1, N_GRID_BC
         IF (GRID_BC(I)%PHYSICAL_GROUP_NAME == STRARRAY(1)) IPG = I
      END DO
      IF (IPG == -1) THEN
         WRITE(*,*) 'Group ', STRARRAY(1), ' not found.'
         CALL ERROR_ABORT('Error in boundary condition definition. Group name not found.')
      END IF
      
      IF (STRARRAY(2) == 'vacuum') THEN
         GRID_BC(IPG)%PARTICLE_BC = VACUUM
      ELSE IF (STRARRAY(2) == 'specular') THEN
         GRID_BC(IPG)%PARTICLE_BC = SPECULAR
      ELSE IF (STRARRAY(2) == 'diffuse') THEN
         GRID_BC(IPG)%PARTICLE_BC = DIFFUSE
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%WALL_TEMP
      ELSE IF (STRARRAY(2) == 'cll') THEN
         GRID_BC(IPG)%PARTICLE_BC = CLL
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%WALL_TEMP
         READ(STRARRAY(4), '(ES14.0)') GRID_BC(IPG)%ACC_N
         READ(STRARRAY(5), '(ES14.0)') GRID_BC(IPG)%ACC_T
      ELSE IF (STRARRAY(2) == 'react') THEN
         GRID_BC(IPG)%REACT = .TRUE.
      ELSE IF (STRARRAY(2) == 'washboard') THEN
         GRID_BC(IPG)%PARTICLE_BC = WB_BC
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%A
         READ(STRARRAY(4), '(ES14.0)') GRID_BC(IPG)%B
         READ(STRARRAY(5), '(ES14.0)') GRID_BC(IPG)%W
         READ(STRARRAY(6), '(ES14.0)') GRID_BC(IPG)%WALL_TEMP
         READ(STRARRAY(7), '(ES14.0)') GRID_BC(IPG)%ACC_N
         GRID_BC(IPG)%ACC_T = 0
         CALL READ_WASHBOARD_TABLES(IPG)
      ELSE IF (STRARRAY(2) == 'axis') THEN
         GRID_BC(IPG)%PARTICLE_BC = AXIS
         ! Needs more info
      ELSE IF (STRARRAY(2) == 'emit') THEN
         GRID_BC(IPG)%PARTICLE_BC = EMIT
         ! Needs more info
      !!! Field BCs
      ELSE IF (STRARRAY(2) == 'dirichlet') THEN
         GRID_BC(IPG)%FIELD_BC = DIRICHLET_BC
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%WALL_POTENTIAL
      ELSE IF (STRARRAY(2) == 'neumann') THEN
         GRID_BC(IPG)%FIELD_BC = NEUMANN_BC
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%WALL_EFIELD
      ELSE IF (STRARRAY(2) == 'dielectric') THEN
         GRID_BC(IPG)%FIELD_BC = DIELECTRIC_BC
      ELSE IF (STRARRAY(2) == 'conductive') THEN
         GRID_BC(IPG)%FIELD_BC = CONDUCTIVE_BC
         BOOL_CONDUCTIVE_BC = .TRUE.
      !!! BCs for both particles and field
      ELSE IF (STRARRAY(2) == 'periodic_master') THEN
         GRID_BC(IPG)%PARTICLE_BC = PERIODIC_MASTER
         GRID_BC(IPG)%FIELD_BC = PERIODIC_SLAVE_BC
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%TRANSLATEVEC(1)
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%TRANSLATEVEC(2)

         ! DO J = 1, U2D_GRID%NUM_CELLS
         !    DO K = 1, 3
         !       IF (U2D_GRID%CELL_EDGES_PG(K, J) == IPG) THEN
         !          ! This is the periodic edge. Find its corresponding neighbor.
         !       END IF
         !    END DO
         ! END DO
      ELSE IF (STRARRAY(2) == 'rf_voltage') THEN
         GRID_BC(IPG)%FIELD_BC = RF_VOLTAGE_BC
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%WALL_RF_POTENTIAL
         READ(STRARRAY(4), '(ES14.0)') GRID_BC(IPG)%WALL_POTENTIAL
         READ(STRARRAY(5), '(ES14.0)') GRID_BC(IPG)%RF_FREQUENCY
      ELSE IF (STRARRAY(2) == 'decoupled_rf_voltage') THEN
         GRID_BC(IPG)%FIELD_BC = DECOUPLED_RF_VOLTAGE_BC
         READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%WALL_RF_POTENTIAL
         READ(STRARRAY(4), '(ES14.0)') GRID_BC(IPG)%WALL_POTENTIAL
         READ(STRARRAY(5), '(ES14.0)') GRID_BC(IPG)%RF_FREQUENCY
         READ(STRARRAY(6), '(ES14.0)') GRID_BC(IPG)%CAPACITANCE
      ELSE
         CALL ERROR_ABORT('Error in boundary condition definition.')
      END IF


   END SUBROUTINE DEF_BOUNDARY_CONDITION


   SUBROUTINE READ_WASHBOARD_TABLES(IPG)

      IMPLICIT NONE
      
      INTEGER :: IPG ! Physical group index
      INTEGER, PARAMETER :: file_unit = 938
      INTEGER :: I, J, K, IDX
      INTEGER :: STATUS

      real, dimension(45, 1600) :: matrix, matrix_neg, matrix_pup

      ALLOCATE(GRID_BC(IPG)%MAX_P_DN(45,40,40))
      ALLOCATE(GRID_BC(IPG)%MAX_P_UP(45,40,40))
      ALLOCATE(GRID_BC(IPG)%P_COLL_UP(45,40,40))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!   Max_P   !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      ! Open the file
      open(unit=file_unit, file='Max_P.csv', status='old', action='read', &
           form='formatted', iostat=status)
      if (status /= 0) then
         WRITE(*,*) 'Error opening file Max_P.csv'
         stop
      endif
      ! Read data from the file into the array
      do i = 1, 45
         read(file_unit, *) (matrix(i, j), j = 1, 1600)
      end do
      ! Close the file
      close(unit=file_unit)
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!   Max_P_neg  !!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Open the file
      open(unit=file_unit, file='Max_neg_P.csv', status='old', action='read', &
           form='formatted', iostat=status)
      if (status /= 0) then
         WRITE(*,*) 'Error opening file neg Max_neg_P.csv'
         stop
      endif
      ! Read data from the file into the array
      do i = 1, 45
          read(file_unit, *) (matrix_neg(i, j), j = 1, 1600)
      end do
      close(unit=file_unit)
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!   P_up_col   !!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Open the file
      open(unit=file_unit, file='Pup_write.csv', status='old', action='read', &
           form='formatted', iostat=status)
      if (status /= 0) then
         WRITE(*,*) 'Error opening file neg Pup_write.csv'
         stop
      endif
      ! Read data from the file into the array
      do i = 1, 45
         read(file_unit, *) (matrix_pup(i, j), j = 1, 1600)
      end do
      close(unit=file_unit)
  
      ! Reshape the matrix
      do i = 1, 45
         idx = 1
         do j = 1, 40
            do k = 1, 40
               GRID_BC(IPG)%MAX_P_DN(i, j, k) = matrix(i, idx)
               GRID_BC(IPG)%MAX_P_UP(i, j, k) = matrix_neg(i, idx)
               GRID_BC(IPG)%P_COLL_UP(i, j, k) = matrix_pup(i, idx)
               idx = idx + 1
            end do
         end do
      end do

   END SUBROUTINE READ_WASHBOARD_TABLES


   SUBROUTINE DEF_DOMAIN_TYPE(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR, I, IPG
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)


      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      ! phys_group type parameters
      IPG = -1
      DO I = 1, N_GRID_BC
         IF (GRID_BC(I)%PHYSICAL_GROUP_NAME == STRARRAY(1)) IPG = I
      END DO
      IF (IPG == -1) THEN
         WRITE(*,*) 'Group ', STRARRAY(1), ' not found.'
         CALL ERROR_ABORT('Error in domain type definition. Group name not found.')
      END IF
      
      IF (STRARRAY(2) == 'fluid') THEN
         GRID_BC(IPG)%VOLUME_BC = FLUID
      ELSE IF (STRARRAY(2) == 'solid') THEN
         GRID_BC(IPG)%VOLUME_BC = SOLID
      ELSE
         CALL ERROR_ABORT('Error in boundary condition definition.')
      END IF

      READ(STRARRAY(3), '(ES14.0)') GRID_BC(IPG)%EPS_REL


   END SUBROUTINE DEF_DOMAIN_TYPE



   SUBROUTINE DEF_BOUNDARY_DUMP_FLUXES(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR, I, IPG
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)


      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      IF (N_STR .NE. 1) CALL ERROR_ABORT('Error in dump fluxes definition.')

      ! phys_group type parameters
      IPG = -1
      DO I = 1, N_GRID_BC
         IF (GRID_BC(I)%PHYSICAL_GROUP_NAME == STRARRAY(1)) IPG = I
      END DO
      IF (IPG == -1) CALL ERROR_ABORT('Error in boundary condition definition. Group name not found.')

      GRID_BC(IPG)%DUMP_FLUXES = .TRUE.


   END SUBROUTINE DEF_BOUNDARY_DUMP_FLUXES



   SUBROUTINE DEF_BOUNDARY_EMIT(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      CHARACTER*64 :: MIX_NAME
      INTEGER      :: MIX_ID, I, IC, IPG
      TYPE(EMIT_TASK_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: TEMP_EMIT_TASKS

      REAL(KIND=8) :: NRHO, UX, UY, UZ, TTRA, TROT, TVIB

      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      IPG = -1
      DO I = 1, N_GRID_BC
         IF (GRID_BC(I)%PHYSICAL_GROUP_NAME == STRARRAY(1)) IPG = I
      END DO
      IF (IPG == -1) CALL ERROR_ABORT('Error in boundary emit definition. Group name not found.')


      READ(STRARRAY(2),'(A10)') MIX_NAME
      READ(STRARRAY(3), '(ES14.0)') NRHO
      READ(STRARRAY(4), '(ES14.0)') UX
      READ(STRARRAY(5), '(ES14.0)') UY
      READ(STRARRAY(6), '(ES14.0)') UZ
      READ(STRARRAY(7), '(ES14.0)') TTRA
      READ(STRARRAY(8),'(ES14.0)') TROT
      READ(STRARRAY(9),'(ES14.0)') TVIB


      MIX_ID = MIXTURE_NAME_TO_ID(MIX_NAME)

      ! WRITE(*,*) 'Read boundary emit definition. Parameters: ', IPG, ', ', MIX_NAME, ', ', MIX_ID, ', ', NRHO, ', ',&
      !  UX, ', ', UY, ', ', UZ, ', ', TTRA, ', ', TROT, ', ', TVIB
      IF (DIMS == 1) THEN 
         DO IC = 1, NCELLS
            DO I = 1, 2
               IF (U1D_GRID%CELL_EDGES_PG(I,IC) == IPG) THEN
                  
                  IF (ALLOCATED(EMIT_TASKS)) THEN
                     ALLOCATE(TEMP_EMIT_TASKS(N_EMIT_TASKS+1)) ! Append the mixture to the list
                     TEMP_EMIT_TASKS(1:N_EMIT_TASKS) = EMIT_TASKS(1:N_EMIT_TASKS)
                     CALL MOVE_ALLOC(TEMP_EMIT_TASKS, EMIT_TASKS)
                  ELSE
                     ALLOCATE(EMIT_TASKS(1))
                  END IF
                  N_EMIT_TASKS = N_EMIT_TASKS + 1
   
                  EMIT_TASKS(N_EMIT_TASKS)%NRHO = NRHO
                  EMIT_TASKS(N_EMIT_TASKS)%UX = UX
                  EMIT_TASKS(N_EMIT_TASKS)%UY = UY
                  EMIT_TASKS(N_EMIT_TASKS)%UZ = UZ
                  EMIT_TASKS(N_EMIT_TASKS)%TTRA = TTRA
                  EMIT_TASKS(N_EMIT_TASKS)%TROT = TROT
                  EMIT_TASKS(N_EMIT_TASKS)%TVIB = TVIB
                  EMIT_TASKS(N_EMIT_TASKS)%MIX_ID = MIX_ID
                  EMIT_TASKS(N_EMIT_TASKS)%IC = IC
                  EMIT_TASKS(N_EMIT_TASKS)%IFACE = I
   
                  EMIT_TASKS(N_EMIT_TASKS)%IV1 = U1D_GRID%CELL_NODES(I,IC)
   
                  ! NFS WILL BE INITIALIZED LATER.
               END IF
            END DO
         END DO
      ELSE IF (DIMS == 2) THEN
         DO IC = 1, NCELLS
            DO I = 1, 3
               IF (U2D_GRID%CELL_EDGES_PG(I,IC) == IPG) THEN
                  
                  IF (ALLOCATED(EMIT_TASKS)) THEN
                     ALLOCATE(TEMP_EMIT_TASKS(N_EMIT_TASKS+1)) ! Append the mixture to the list
                     TEMP_EMIT_TASKS(1:N_EMIT_TASKS) = EMIT_TASKS(1:N_EMIT_TASKS)
                     CALL MOVE_ALLOC(TEMP_EMIT_TASKS, EMIT_TASKS)
                  ELSE
                     ALLOCATE(EMIT_TASKS(1))
                  END IF
                  N_EMIT_TASKS = N_EMIT_TASKS + 1
   
                  EMIT_TASKS(N_EMIT_TASKS)%NRHO = NRHO
                  EMIT_TASKS(N_EMIT_TASKS)%UX = UX
                  EMIT_TASKS(N_EMIT_TASKS)%UY = UY
                  EMIT_TASKS(N_EMIT_TASKS)%UZ = UZ
                  EMIT_TASKS(N_EMIT_TASKS)%TTRA = TTRA
                  EMIT_TASKS(N_EMIT_TASKS)%TROT = TROT
                  EMIT_TASKS(N_EMIT_TASKS)%TVIB = TVIB
                  EMIT_TASKS(N_EMIT_TASKS)%MIX_ID = MIX_ID
                  EMIT_TASKS(N_EMIT_TASKS)%IC = IC
                  EMIT_TASKS(N_EMIT_TASKS)%IFACE = I
   
                  
                  IF (I == 1) THEN
                     EMIT_TASKS(N_EMIT_TASKS)%IV1 = U2D_GRID%CELL_NODES(1,IC)
                     EMIT_TASKS(N_EMIT_TASKS)%IV2 = U2D_GRID%CELL_NODES(2,IC)
                  ELSE IF (I == 2) THEN
                     EMIT_TASKS(N_EMIT_TASKS)%IV1 = U2D_GRID%CELL_NODES(2,IC)
                     EMIT_TASKS(N_EMIT_TASKS)%IV2 = U2D_GRID%CELL_NODES(3,IC)
                  ELSE
                     EMIT_TASKS(N_EMIT_TASKS)%IV1 = U2D_GRID%CELL_NODES(3,IC)
                     EMIT_TASKS(N_EMIT_TASKS)%IV2 = U2D_GRID%CELL_NODES(1,IC)
                  END IF
   
                  ! NFS WILL BE INITIALIZED LATER.
               END IF
            END DO
         END DO
      ELSE IF (DIMS == 3) THEN
         DO IC = 1, NCELLS
            DO I = 1, 4
               IF (U3D_GRID%CELL_FACES_PG(I,IC) == IPG) THEN
                  IF (ALLOCATED(EMIT_TASKS)) THEN
                     ALLOCATE(TEMP_EMIT_TASKS(N_EMIT_TASKS+1)) ! Append the mixture to the list
                     TEMP_EMIT_TASKS(1:N_EMIT_TASKS) = EMIT_TASKS(1:N_EMIT_TASKS)
                     CALL MOVE_ALLOC(TEMP_EMIT_TASKS, EMIT_TASKS)
                  ELSE
                     ALLOCATE(EMIT_TASKS(1))
                  END IF
                  N_EMIT_TASKS = N_EMIT_TASKS + 1
   
                  EMIT_TASKS(N_EMIT_TASKS)%NRHO = NRHO
                  EMIT_TASKS(N_EMIT_TASKS)%UX = UX
                  EMIT_TASKS(N_EMIT_TASKS)%UY = UY
                  EMIT_TASKS(N_EMIT_TASKS)%UZ = UZ
                  EMIT_TASKS(N_EMIT_TASKS)%TTRA = TTRA
                  EMIT_TASKS(N_EMIT_TASKS)%TROT = TROT
                  EMIT_TASKS(N_EMIT_TASKS)%TVIB = TVIB
                  EMIT_TASKS(N_EMIT_TASKS)%MIX_ID = MIX_ID
                  EMIT_TASKS(N_EMIT_TASKS)%IC = IC
                  EMIT_TASKS(N_EMIT_TASKS)%IFACE = I
   
                  ! NFS WILL BE INITIALIZED LATER.
               END IF
            END DO
         END DO
      END IF

      


   END SUBROUTINE DEF_BOUNDARY_EMIT





   SUBROUTINE DEF_INITIAL_PARTICLES(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      CHARACTER*64 :: MIX_NAME
      INTEGER      :: MIX_ID
      TYPE(INITIAL_PARTICLES_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: TEMP_INITIAL_PARTICLES_TASK

      REAL(KIND=8) :: NRHO, UX, UY, UZ, TTRAX, TTRAY, TTRAZ, TROT, TVIB

      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)


      READ(STRARRAY(1),'(A10)') MIX_NAME
      READ(STRARRAY(2), '(ES14.0)') NRHO
      READ(STRARRAY(3), '(ES14.0)') UX
      READ(STRARRAY(4), '(ES14.0)') UY
      READ(STRARRAY(5), '(ES14.0)') UZ
      READ(STRARRAY(6), '(ES14.0)') TTRAX
      READ(STRARRAY(7), '(ES14.0)') TTRAY
      READ(STRARRAY(8), '(ES14.0)') TTRAZ
      READ(STRARRAY(9), '(ES14.0)') TROT
      READ(STRARRAY(10),'(ES14.0)') TVIB

      MIX_ID = MIXTURE_NAME_TO_ID(MIX_NAME)

      ! WRITE(*,*) 'Read initial particle seed definition. Parameters: ', MIX_NAME, ', ', MIX_ID, ', ', NRHO, ', ',&
      !  UX, ', ', UY, ', ', UZ, ', ', TTRAX, ', ', TTRAY, ', ', TTRAZ, ', ', TROT, ', ', TVIB

      IF (ALLOCATED(INITIAL_PARTICLES_TASKS)) THEN
         ALLOCATE(TEMP_INITIAL_PARTICLES_TASK(N_INITIAL_PARTICLES_TASKS+1)) ! Append the mixture to the list
         TEMP_INITIAL_PARTICLES_TASK(1:N_INITIAL_PARTICLES_TASKS) = INITIAL_PARTICLES_TASKS(1:N_INITIAL_PARTICLES_TASKS)
         CALL MOVE_ALLOC(TEMP_INITIAL_PARTICLES_TASK, INITIAL_PARTICLES_TASKS)
      ELSE
         ALLOCATE(INITIAL_PARTICLES_TASKS(1))
      END IF
      N_INITIAL_PARTICLES_TASKS = N_INITIAL_PARTICLES_TASKS + 1

      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%NRHO = NRHO
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%UX = UX
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%UY = UY
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%UZ = UZ
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%TTRAX = TTRAX
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%TTRAY = TTRAY
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%TTRAZ = TTRAZ
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%TROT = TROT
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%TVIB = TVIB
      INITIAL_PARTICLES_TASKS(N_INITIAL_PARTICLES_TASKS)%MIX_ID = MIX_ID

   END SUBROUTINE DEF_INITIAL_PARTICLES


   SUBROUTINE DEF_VOLUME_INJECT(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      CHARACTER*64 :: MIX_NAME
      INTEGER      :: MIX_ID
      TYPE(VOLUME_INJECT_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: TEMP_VOLUME_INJECT_TASK

      REAL(KIND=8) :: NRHODOT, UX, UY, UZ, TTRAX, TTRAY, TTRAZ, TROT, TVIB

      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)


      READ(STRARRAY(1),'(A10)') MIX_NAME
      READ(STRARRAY(2), '(ES14.0)') NRHODOT
      READ(STRARRAY(3), '(ES14.0)') UX
      READ(STRARRAY(4), '(ES14.0)') UY
      READ(STRARRAY(5), '(ES14.0)') UZ
      READ(STRARRAY(6), '(ES14.0)') TTRAX
      READ(STRARRAY(7), '(ES14.0)') TTRAY
      READ(STRARRAY(8), '(ES14.0)') TTRAZ
      READ(STRARRAY(9), '(ES14.0)') TROT
      READ(STRARRAY(10),'(ES14.0)') TVIB

      MIX_ID = MIXTURE_NAME_TO_ID(MIX_NAME)

      ! WRITE(*,*) 'Read volume inject definition. Parameters: ', MIX_NAME, ', ', MIX_ID, ', ', NRHODOT, ', ',&
      !  UX, ', ', UY, ', ', UZ, ', ', TTRAX, ', ', TTRAY, ', ', TTRAZ, ', ', TROT, ', ', TVIB

      IF (ALLOCATED(VOLUME_INJECT_TASKS)) THEN
         ALLOCATE(TEMP_VOLUME_INJECT_TASK(N_VOLUME_INJECT_TASKS+1)) ! Append the mixture to the list
         TEMP_VOLUME_INJECT_TASK(1:N_VOLUME_INJECT_TASKS) = VOLUME_INJECT_TASKS(1:N_VOLUME_INJECT_TASKS)
         CALL MOVE_ALLOC(TEMP_VOLUME_INJECT_TASK, VOLUME_INJECT_TASKS)
      ELSE
         ALLOCATE(VOLUME_INJECT_TASKS(1))
      END IF
      N_VOLUME_INJECT_TASKS = N_VOLUME_INJECT_TASKS + 1

      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%NRHODOT = NRHODOT
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%UX = UX
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%UY = UY
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%UZ = UZ
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%TTRAX = TTRAX
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%TTRAY = TTRAY
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%TTRAZ = TTRAZ
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%TROT = TROT
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%TVIB = TVIB
      VOLUME_INJECT_TASKS(N_VOLUME_INJECT_TASKS)%MIX_ID = MIX_ID

   END SUBROUTINE DEF_VOLUME_INJECT


   SUBROUTINE DEF_SOLENOID(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      TYPE(SOLENOID), DIMENSION(:), ALLOCATABLE :: TEMP_SOLENOIDS
      
      IF (ALLOCATED(SOLENOIDS)) THEN
         ALLOCATE(TEMP_SOLENOIDS(N_SOLENOIDS+1))
         TEMP_SOLENOIDS(1:N_SOLENOIDS) = SOLENOIDS(1:N_SOLENOIDS)
         CALL MOVE_ALLOC(TEMP_SOLENOIDS, SOLENOIDS)
      ELSE
         ALLOCATE(SOLENOIDS(1))
      END IF
      N_SOLENOIDS = N_SOLENOIDS + 1

      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      READ(STRARRAY(1), '(ES14.0)') SOLENOIDS(N_SOLENOIDS)%X1
      READ(STRARRAY(2), '(ES14.0)') SOLENOIDS(N_SOLENOIDS)%Y1
      READ(STRARRAY(3), '(ES14.0)') SOLENOIDS(N_SOLENOIDS)%X2
      READ(STRARRAY(4), '(ES14.0)') SOLENOIDS(N_SOLENOIDS)%Y2
      READ(STRARRAY(5), '(I3)') SOLENOIDS(N_SOLENOIDS)%N_WIRES_X
      READ(STRARRAY(6), '(I3)') SOLENOIDS(N_SOLENOIDS)%N_WIRES_Y
      READ(STRARRAY(7), '(ES14.0)') SOLENOIDS(N_SOLENOIDS)%WIRE_CURRENT
      
      

   END SUBROUTINE DEF_SOLENOID



   SUBROUTINE READ_REACTIONS(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in3 = 457
      INTEGER, PARAMETER :: in4 = 754
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF, ReasonEOFCS

      CHARACTER(LEN=80)   :: DEFINITION
      INTEGER :: N_STR, index
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      TYPE(REACTIONS_DATA_STRUCTURE) :: NEW_REACTION
      INTEGER :: I, NROWS
      CHARACTER*512      :: LINECS
      CHARACTER*256      :: REACTION_FILENAME
      
      !CHARACTER*64 :: SP_NAME
      !INTEGER      :: SP_ID
      !REAL(KIND=8) :: DIAM
      !REAL(KIND=8) :: OMEGA
      !REAL(KIND=8) :: TREF
      !REAL(KIND=8) :: ALPHA

      !INTEGER      :: IS, JS
      !REAL(KIND=8) :: M1, M2, MRED


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
         IF (STRARRAY(2) .NE. '+' .OR. (STRARRAY(4) .NE. '-->' .AND. STRARRAY(4) .NE. '-CEX->') .OR. STRARRAY(6) .NE. '+') THEN
            WRITE(*,*) 'Line in reactions file:', DEFINITION
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
               WRITE(*,*) 'Line in reactions file:', DEFINITION
               CALL ERROR_ABORT('Attention, format is not respected in reactions file.')
            END IF
            NEW_REACTION%N_PROD = 3
            NEW_REACTION%P3_SP_ID = REACTION_SPECIES_NAME_TO_ID(STRARRAY(9))
         END IF

         IF (N_STR .GE. 11) THEN
            IF (STRARRAY(10) .NE. '+') THEN
               WRITE(*,*) 'Line in reactions file:', DEFINITION
               CALL ERROR_ABORT('Attention, format is not respected in reactions file.')
            END IF
            NEW_REACTION%N_PROD = 4
            NEW_REACTION%P4_SP_ID = REACTION_SPECIES_NAME_TO_ID(STRARRAY(11))
         END IF
      
         IF (STRARRAY(4) == '-CEX->') THEN
            NEW_REACTION%IS_CEX = .TRUE.
         ELSE
            NEW_REACTION%IS_CEX = .FALSE.
         END IF

         READ(in3,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction parameters line
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)
         IF (STRARRAY(1) == 'tce') THEN
            NEW_REACTION%TYPE = TCE
            IF (NEW_REACTION%IS_CEX) THEN
               NEW_REACTION%EA = 0.d0
               READ(STRARRAY(2), '(ES14.0)') NEW_REACTION%C1
               READ(STRARRAY(3), '(ES14.0)') NEW_REACTION%C2
            ELSE
               READ(STRARRAY(2), '(ES14.0)') NEW_REACTION%A
               READ(STRARRAY(3), '(ES14.0)') NEW_REACTION%N
               READ(STRARRAY(4), '(ES14.0)') NEW_REACTION%EA
            END IF
         ELSE IF (STRARRAY(1) == 'lxcat') THEN
            NEW_REACTION%TYPE = LXCAT
            READ(STRARRAY(2), *) NEW_REACTION%EA
            READ(STRARRAY(3), *) REACTION_FILENAME


            OPEN(UNIT=in4,FILE=REACTION_FILENAME, STATUS='old',IOSTAT=ios)

            IF (ios .NE. 0) THEN
               CALL ERROR_ABORT('Attention, reactions cross section file not found! ABORTING.')
            ENDIF
      
            LINECS = '' ! Init empty

            DO
               READ(in4,'(A)', IOSTAT=ReasonEOFCS) LINECS
               IF (ReasonEOFCS < 0) CALL ERROR_ABORT('Attention, reactions cross section file format error! ABORTING.')
                  
               IF (LINECS(1:5) == '-----') EXIT
            END DO
            NROWS = 0
            DO
               READ(in4,'(A)', IOSTAT=ReasonEOFCS) LINECS  
               IF (ReasonEOFCS < 0) CALL ERROR_ABORT('Attention, reactions cross section file format error! ABORTING.')
                  
               IF (LINECS(1:5) == '-----') EXIT
               NROWS = NROWS + 1
            END DO
            REWIND(in4)
            IF (ALLOCATED(NEW_REACTION%TABLE_ENERGY)) DEALLOCATE(NEW_REACTION%TABLE_ENERGY)
            IF (ALLOCATED(NEW_REACTION%TABLE_CS)) DEALLOCATE(NEW_REACTION%TABLE_CS)
            ALLOCATE(NEW_REACTION%TABLE_ENERGY(NROWS))
            ALLOCATE(NEW_REACTION%TABLE_CS(NROWS))
            DO
               READ(in4,'(A)', IOSTAT=ReasonEOFCS) LINECS
               IF (ReasonEOFCS < 0) CALL ERROR_ABORT('Attention, reactions cross section file format error! ABORTING.')
                  
               IF (LINECS(1:5) == '-----') EXIT
            END DO
            NEW_REACTION%MAX_SIGMA = 0
            DO I = 1, NROWS
               READ(in4,'(A)', IOSTAT=ReasonEOFCS) LINECS
               IF (ReasonEOFCS < 0) CALL ERROR_ABORT('Attention, reactions cross section file format error! ABORTING.')
               READ(LINECS, *) NEW_REACTION%TABLE_ENERGY(I), NEW_REACTION%TABLE_CS(I)
               IF (NEW_REACTION%TABLE_CS(I) > SIGMAMAX) SIGMAMAX = NEW_REACTION%TABLE_CS(I) ! Global maximum gross section
               IF (NEW_REACTION%TABLE_CS(I) > NEW_REACTION%MAX_SIGMA) NEW_REACTION%MAX_SIGMA = NEW_REACTION%TABLE_CS(I) ! Reaction maximum cross section
            END DO

            CLOSE(in4)

            NEW_REACTION%TABLE_ENERGY = NEW_REACTION%TABLE_ENERGY * QE ! Table is in [eV] but we need [J].

         END IF

         NEW_REACTION%EA = NEW_REACTION%EA * QE
         NEW_REACTION%COUNTS = 0

         IF (ReasonEOF < 0) EXIT ! End of file reached
         
         IF (ALLOCATED(REACTIONS)) THEN
            ALLOCATE(TEMP_REACTIONS(N_REACTIONS+1)) ! Append the reaction to the list
            TEMP_REACTIONS(1:N_REACTIONS) = REACTIONS(1:N_REACTIONS)
            CALL MOVE_ALLOC(TEMP_REACTIONS, REACTIONS)
         ELSE
            ALLOCATE(REACTIONS(1))
         END IF
         N_REACTIONS = N_REACTIONS + 1
         REACTIONS(N_REACTIONS) = NEW_REACTION

         DEALLOCATE(STRARRAY)

      END DO
      
      CLOSE(in3) ! Close input file

      ! IF (PROC_ID == 0) THEN
      !    DO index = 1, N_REACTIONS
      !       WRITE(*,*) 'Reaction ', index, 'Has 2 reactants with ids:', REACTIONS(index)%R1_SP_ID, ' and ', &
      !       REACTIONS(index)%R2_SP_ID
      !       IF (REACTIONS(index)%N_PROD == 2) THEN
      !          WRITE(*,*) REACTIONS(index)%N_PROD, 'products with ids:',  REACTIONS(index)%P1_SP_ID, ' and ', &
      !          REACTIONS(index)%P2_SP_ID
      !          IF (REACTIONS(index)%IS_CEX) WRITE(*,*) 'This is a CEX reaction'
      !       ELSE IF (REACTIONS(index)%N_PROD == 3) THEN
      !          WRITE(*,*) REACTIONS(index)%N_PROD, 'products with ids:',  REACTIONS(index)%P1_SP_ID, ', ', &
      !          REACTIONS(index)%P2_SP_ID, ' and ', REACTIONS(index)%P3_SP_ID
      !       END IF
      !       IF (REACTIONS(index)%TYPE == TCE) THEN
      !          WRITE(*,*) 'Parameters:', REACTIONS(index)%A, REACTIONS(index)%N, REACTIONS(index)%EA
      !       ELSE IF (REACTIONS(index)%TYPE == LXCAT) THEN
      !          WRITE(*,*) 'Reaction ', index, ' is from tabulated data in LxCat format.'
      !          WRITE(*,*) 'Activation energy: ', REACTIONS(index)%EA
      !          WRITE(*,*) 'Here are the energies (eV): ', REACTIONS(index)%TABLE_ENERGY
      !          WRITE(*,*) 'And here are the cross sections (m^2): ', REACTIONS(index)%TABLE_CS
      !       ELSE
      !          WRITE(*,*) 'Reaction ', index, ' is not defined!'
      !       END IF
      !    END DO
      ! END IF

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
         READ(STRARRAY(1), '(ES14.0)') NEW_REACTION%PROB

         IF (ReasonEOF < 0) EXIT ! End of file reached
         
         IF (ALLOCATED(WALL_REACTIONS)) THEN
            ALLOCATE(TEMP_WALL_REACTIONS(N_WALL_REACTIONS+1)) ! Append the mixture to the list
            TEMP_WALL_REACTIONS(1:N_WALL_REACTIONS) = WALL_REACTIONS(1:N_WALL_REACTIONS)
            CALL MOVE_ALLOC(TEMP_WALL_REACTIONS, WALL_REACTIONS)
         ELSE
            ALLOCATE(WALL_REACTIONS(1))
         END IF
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

      IF (DIMS == 1) THEN
         DO I = 1, NX
            CELL_VOLUMES(I) = XSIZE(I)*(YMAX-YMIN)*(ZMAX-ZMIN)
         END DO
      ELSE IF (DIMS == 2) THEN
         DO I = 1, NX
            DO J = 1, NY
               IF (.NOT. AXI) THEN
                  CELL_VOLUMES(I + NX*(J-1)) = XSIZE(I)*YSIZE(J)*(ZMAX-ZMIN)
               ELSE
                  ! 2D Axisymmetric geometry: volume is a cylindrical annulus. Coordinate Z is in radians.
                  CELL_VOLUMES(I + NX*(J-1)) = 0.5*XSIZE(I)*(YCOORD(J+1)**2 - YCOORD(J)**2)*(ZMAX-ZMIN)
               END IF
            END DO
         END DO
      END IF
      ! Perform all the checks on these arrays!
      CLOSE(in5)

   END SUBROUTINE READ_GRID_FILE


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

      INTEGER            :: IP, NP_INIT, IC
      REAL(KIND=8)       :: Xp, Yp, Zp, VXp, VYp, VZp, EROT, EVIB, DOMAIN_VOLUME, VOL
      INTEGER            :: CID

      REAL(KIND=8), DIMENSION(3) :: V1, V2, V3, V4
      REAL(KIND=8)       :: P, Q, R, S, T, U

      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
      REAL(KIND=8)  :: M
      INTEGER       :: S_ID, ELECTRON_S_ID, i, ITASK, CELL_PG

      ! Print message 
      CALL ONLYMASTERPRINT1(PROC_ID, '> SEEDING INITIAL PARTICLES IN THE DOMAIN...')
     
      IF (COLOCATED_ELECTRONS) ELECTRON_S_ID = SPECIES_NAME_TO_ID('e')

      ! Compute number of particles to be seeded
      !NP_INIT = NINT(NRHO_INIT/FNUM*(XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN))

      !CALL ONLYMASTERPRINT2(PROC_ID, '   Particles to be seeded', REAL(NP_INIT, KIND=8))
 
      ! Compute number of particles to be seeded by every process
      !NPPP_INIT = INT(NP_INIT/N_MPI_THREADS) 

      DO ITASK = 1, N_INITIAL_PARTICLES_TASKS
         ! Create particles in the domain
         DO i = 1, MIXTURES(INITIAL_PARTICLES_TASKS(ITASK)%MIX_ID)%N_COMPONENTS
            
            S_ID = MIXTURES(INITIAL_PARTICLES_TASKS(ITASK)%MIX_ID)%COMPONENTS(i)%ID

            !IF (COLOCATED_ELECTRONS .AND. S_ID == ELECTRON_S_ID) CYCLE

            IF (GRID_TYPE == UNSTRUCTURED) THEN
               DO IC = 1, NCELLS
                  IF (DIMS == 1) THEN
                     IF (CELL_PROCS(IC) .NE. PROC_ID) CYCLE
                     CELL_PG = U1D_GRID%CELL_PG(IC)
                  ELSE IF (DIMS == 2) THEN
                     IF (CELL_PROCS(IC) .NE. PROC_ID) CYCLE
                     CELL_PG = U2D_GRID%CELL_PG(IC)
                  ELSE IF (DIMS == 3) THEN
                     IF (CELL_PROCS(IC) .NE. PROC_ID) CYCLE
                     CELL_PG = U3D_GRID%CELL_PG(IC)
                  END IF
                  IF (CELL_PG .NE. -1) THEN
                     IF (GRID_BC(CELL_PG)%VOLUME_BC == SOLID) CYCLE
                  END IF


                  ! Compute number of particles of this species to be created in this cell.
                  IF (DIMS == 1) THEN
                     VOL = U1D_GRID%CELL_VOLUMES(IC)
                  ELSE IF (DIMS == 2) THEN
                     VOL = U2D_GRID%CELL_VOLUMES(IC)
                  ELSE IF (DIMS == 3) THEN
                     VOL = U3D_GRID%CELL_VOLUMES(IC)
                  END IF
                  NP_INIT = RANDINT(INITIAL_PARTICLES_TASKS(ITASK)%NRHO/(FNUM*SPECIES(S_ID)%SPWT)*VOL* &
                              MIXTURES(INITIAL_PARTICLES_TASKS(ITASK)%MIX_ID)%COMPONENTS(i)%MOLFRAC)
                  IF (NP_INIT == 0) CYCLE


                  IF (DIMS == 1) THEN
                     V1 = U1D_GRID%NODE_COORDS(:,U1D_GRID%CELL_NODES(1,IC))
                     V2 = U1D_GRID%NODE_COORDS(:,U1D_GRID%CELL_NODES(2,IC))
                     V3 = 0
                     V4 = 0
                  ELSE IF (DIMS == 2) THEN
                     V1 = U2D_GRID%NODE_COORDS(:,U2D_GRID%CELL_NODES(1,IC))
                     V2 = U2D_GRID%NODE_COORDS(:,U2D_GRID%CELL_NODES(2,IC))
                     V3 = U2D_GRID%NODE_COORDS(:,U2D_GRID%CELL_NODES(3,IC))
                     V4 = 0
                  ELSE IF (DIMS == 3) THEN
                     V1 = U3D_GRID%NODE_COORDS(:,U3D_GRID%CELL_NODES(1,IC))
                     V2 = U3D_GRID%NODE_COORDS(:,U3D_GRID%CELL_NODES(2,IC))
                     V3 = U3D_GRID%NODE_COORDS(:,U3D_GRID%CELL_NODES(3,IC))
                     V4 = U3D_GRID%NODE_COORDS(:,U3D_GRID%CELL_NODES(4,IC))
                  END IF

                  DO IP = 1, NP_INIT

                     ! Create particle position randomly in the cell
                     S = rf()
                     T = rf()
                     U = rf()
                     
                     IF (S+T > 1) THEN
                        S = 1-S; T = 1-T
                     END IF

                     IF (DIMS == 1) THEN
                        XP = V1(1) + (V2(1)-V1(1))*S
                        YP = YMIN + (YMAX-YMIN)*T
                        ZP = ZMIN + (ZMAX-ZMIN)*U
                     ELSE IF (DIMS == 2) THEN
                        XP = V1(1) + (V2(1)-V1(1))*S + (V3(1)-V1(1))*T
                        YP = V1(2) + (V2(2)-V1(2))*S + (V3(2)-V1(2))*T
                        IF (AXI) THEN
                           ZP = 0.d0
                        ELSE
                           ZP = ZMIN + (ZMAX-ZMIN)*U
                        END IF
                     ELSE IF (DIMS == 3) THEN
                        ! http://vcg.isti.cnr.it/publications/papers/rndtetra_a.pdf

                        IF (S+T+U <= 1) THEN
                           P = S; Q = T; R = U
                        ELSE
                           IF (T+U > 1) THEN
                              P = S; Q = 1-U; R = 1-S-T
                           ELSE
                              P = 1-T-U; Q = T; R = S+T+U-1
                           END IF
                        END IF

                        XP = V1(1) + (V2(1)-V1(1))*P + (V3(1)-V1(1))*Q + (V4(1)-V1(1))*R
                        YP = V1(2) + (V2(2)-V1(2))*P + (V3(2)-V1(2))*Q + (V4(2)-V1(2))*R
                        ZP = V1(3) + (V2(3)-V1(3))*P + (V3(3)-V1(3))*Q + (V4(3)-V1(3))*R
                     END IF

                     ! Assign velocity and energy following a Boltzmann distribution
                     M = SPECIES(S_ID)%MOLECULAR_MASS
                     CALL MAXWELL(INITIAL_PARTICLES_TASKS(ITASK)%UX, &
                                  INITIAL_PARTICLES_TASKS(ITASK)%UY, &
                                  INITIAL_PARTICLES_TASKS(ITASK)%UZ, &
                                  INITIAL_PARTICLES_TASKS(ITASK)%TTRAX, &
                                  INITIAL_PARTICLES_TASKS(ITASK)%TTRAY, &
                                  INITIAL_PARTICLES_TASKS(ITASK)%TTRAZ, &
                                  VXP, VYP, VZP, M)

                     CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, INITIAL_PARTICLES_TASKS(ITASK)%TROT, EROT)
                     CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, INITIAL_PARTICLES_TASKS(ITASK)%TVIB, EVIB)

                     CALL INIT_PARTICLE(XP,YP,ZP,VXP,VYP,VZP,EROT,EVIB,S_ID,IC,DT, particleNOW) ! Save in particle
                     CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array



                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     ! Modification to inject electrons colocated with ions.
                     IF (COLOCATED_ELECTRONS .AND. SPECIES(S_ID)%CHARGE == 1) THEN
                        
                        CALL MAXWELL(INITIAL_PARTICLES_TASKS(ITASK)%UX, &
                                    INITIAL_PARTICLES_TASKS(ITASK)%UY, &
                                    INITIAL_PARTICLES_TASKS(ITASK)%UZ, &
                                    COLOCATED_ELECTRONS_TTRA, &
                                    COLOCATED_ELECTRONS_TTRA, &
                                    COLOCATED_ELECTRONS_TTRA, &
                                    VXP, VYP, VZP, SPECIES(ELECTRON_S_ID)%MOLECULAR_MASS)

                        CALL INIT_PARTICLE(XP,YP,ZP,VXP,VYP,VZP,0.d0,0.d0,ELECTRON_S_ID,IC,DT, particleNOW) ! Save in particle
                        CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
                     END IF
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  END DO
               END DO
            ELSE ! Structured grid
               ! Compute number of particles of this species per process to be created.
               IF (AXI) THEN
                  DOMAIN_VOLUME = 0.5*(XMAX-XMIN)*(YMAX**2-YMIN**2)*(ZMAX-ZMIN)
               ELSE
                  DOMAIN_VOLUME = (XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN)
               END IF
               NP_INIT = RANDINT(INITIAL_PARTICLES_TASKS(ITASK)%NRHO/(FNUM*SPECIES(S_ID)%SPWT)*DOMAIN_VOLUME* &
                           MIXTURES(INITIAL_PARTICLES_TASKS(ITASK)%MIX_ID)%COMPONENTS(i)%MOLFRAC/N_MPI_THREADS)
               IF (NP_INIT == 0) CYCLE
               DO IP = 1, NP_INIT

                  ! Create particle position randomly in the domain
                  XP = XMIN + (XMAX-XMIN)*rf()

                  IF (AXI .AND. (.NOT. BOOL_RADIAL_WEIGHTING)) THEN
                     YP = SQRT(YMIN*YMIN + rf()*(YMAX*YMAX - YMIN*YMIN))
                     ZP = 0
                  ELSE
                     YP = YMIN + (YMAX-YMIN)*rf()
                     ZP = ZMIN + (ZMAX-ZMIN)*rf()
                  END IF
               
                  ! Assign velocity and energy following a Boltzmann distribution
                  M = SPECIES(S_ID)%MOLECULAR_MASS
                  CALL MAXWELL(INITIAL_PARTICLES_TASKS(ITASK)%UX, &
                               INITIAL_PARTICLES_TASKS(ITASK)%UY, &
                               INITIAL_PARTICLES_TASKS(ITASK)%UZ, &
                               INITIAL_PARTICLES_TASKS(ITASK)%TTRAX, &
                               INITIAL_PARTICLES_TASKS(ITASK)%TTRAY, &
                               INITIAL_PARTICLES_TASKS(ITASK)%TTRAZ, &
                               VXP, VYP, VZP, M)

                  CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, INITIAL_PARTICLES_TASKS(ITASK)%TROT, EROT)
                  CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, INITIAL_PARTICLES_TASKS(ITASK)%TVIB, EVIB)

                  CALL CELL_FROM_POSITION(XP,YP,  CID) ! Find cell containing particle

                  CALL INIT_PARTICLE(XP,YP,ZP,VXP,VYP,VZP,EROT,EVIB,S_ID,CID,DT, particleNOW) ! Save in particle
                  CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
               END DO
            END IF
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
      IF ( GRID_TYPE .NE. UNSTRUCTURED .AND. ((NX < 1) .OR. (NY < 1)) ) THEN
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
      IF (AXI) THEN 
         IF (BOOL_PERIODIC(3)) THEN
            CALL ERROR_ABORT('ERROR! For axisymmetric simulations, Y cannot be periodic.')
         END IF
         IF (YMIN /= 0) THEN
            CALL ERROR_ABORT('ERROR! For axisymmetric simulations, YMIN must be 0.')
         END IF
         IF (.NOT. BOOL_SPECULAR(3) .OR. BOOL_DIFFUSE(3)) THEN
            CALL ERROR_ABORT('ERROR! For axisymmetric simulations, YMIN must have specular B.C.')
         END IF
         IF (BOOL_REACT(3)) THEN
            CALL ERROR_ABORT('ERROR! For axisymmetric simulations, YMIN cannot be reacting.')
         END IF
      END IF
      ! ---------- Check if fluid electron parameters were set for HYBRID PIC_TYPE ------
      IF (PIC_TYPE == HYBRID .AND. (BOLTZ_N0 == 0. .OR. BOLTZ_TE == 0.)) THEN
         WRITE(*,*) BOLTZ_N0, BOLTZ_TE
         CALL ERROR_ABORT('Hybrid PIC Type was chosen but proper parameter setting for fluid electrons is missing.')
      END IF

   END SUBROUTINE INPUT_DATA_SANITY_CHECK


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INITVARIOUS -> initializes various other variables !!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE INITVARIOUS 

      IMPLICIT NONE
 
      INTEGER :: I, J
      REAL(KIND=8) :: CELL_YMIN, CELL_YMAX

      ! ~~~~~~~~ Initialize RNG (random number generator) ~~~~~~~~~~

      ! Create local seed for RNG (different for each processor)
      RNG_SEED_LOCAL = RNG_SEED_GLOBAL*(PROC_ID + 1) 

      CALL init_genrand64(RNG_SEED_LOCAL) ! Seed Mersenne Twister PRNG with local quantity



      CELL_VOL = (XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN)/(NX*NY)

      IF (GRID_TYPE == RECTILINEAR_UNIFORM .AND. AXI) THEN
         ALLOCATE(CELL_VOLUMES(NX*NY))
         IF (BOOL_RADIAL_WEIGHTING) ALLOCATE(CELL_FNUM(NX*NY))
         DO I = 1, NX
            DO J = 1, NY
               CELL_YMIN = YMIN + (J-1)*(YMAX-YMIN)/NY
               CELL_YMAX = YMIN + J*(YMAX-YMIN)/NY
               CELL_VOLUMES(I + NX*(J-1)) = 0.5*(XMAX-XMIN)/NX*(CELL_YMAX**2 - CELL_YMIN**2)*(ZMAX-ZMIN)
               IF (BOOL_RADIAL_WEIGHTING) CELL_FNUM(I + NX*(J-1)) = FNUM*NY*(CELL_YMAX**2 - CELL_YMIN**2)/(YMAX**2-YMIN**2)
            END DO
         END DO
      END IF


      ! Allocate counters for wall collisions, boundary collisions and particle injection.
      ALLOCATE(WALL_COLL_COUNT(N_WALLS*N_SPECIES))
      WALL_COLL_COUNT = 0
      ALLOCATE(BOUNDARY_COLL_COUNT(4*N_SPECIES))
      BOUNDARY_COLL_COUNT = 0
      ALLOCATE(LINE_EMIT_COUNT(N_LINESOURCES*N_SPECIES))
      LINE_EMIT_COUNT = 0

      ! ~~~~~~~~ Initialize grid number of cells and nodes ~~~~~~~~~~

      IF (.NOT. GRID_TYPE == UNSTRUCTURED) THEN
         NCELLS = NX*NY
         IF (DIMS == 1) THEN
            NNODES = NX+1
         ELSE IF (DIMS == 2) THEN
            NNODES = (NX+1)*(NY+1)
         END IF
      END IF

   END SUBROUTINE INITVARIOUS

   
   SUBROUTINE INITFIELDS

      IMPLICIT NONE

      IF (GRID_TYPE == UNSTRUCTURED) THEN
         ALLOCATE(E_FIELD(3, 1, NCELLS))
         E_FIELD = 0.d0
         ALLOCATE(B_FIELD(3, 1, NNODES))
         B_FIELD = 0.d0
         ALLOCATE(EBAR_FIELD(3, 1, NCELLS))
         EBAR_FIELD = 0.d0
         ALLOCATE(SURFACE_CHARGE(NNODES))
         SURFACE_CHARGE = 0.d0
      ELSE
         NPX = NX + 1
         IF (DIMS == 2) THEN
            NPY = NY + 1
         ELSE
            NPY = 1
         END IF
         ALLOCATE(E_FIELD(3, 0:NPY-1, 0:NPX-1))
         E_FIELD = 0.d0
      END IF

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
      REAL(KIND=8) :: U_NORM, S_NORM, FLUXLINESOURCE, LINELENGTH, NORMX, NORMY, NORMZ, AREA
      REAL(KIND=8) :: PI2  
      REAL(KIND=8) :: INTEG1, INTEG2

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: nfs_LINE
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TASK_NFS
      !INTEGER, DIMENSION(:), ALLOCATABLE :: EMIT_COUNT
 
      REAL(KIND=8)  :: M, FRAC
      INTEGER :: N_COMP, IS, S_ID, ILINE, ITASK, I

      REAL(KIND=8) :: X1, X2, Y1, Y2
    
      PI2  = 2.*PI

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
      ! More compact than one for each boundary, try to use this.

      DO ILINE = 1, N_LINESOURCES ! Loop on line sources
         
         N_COMP = MIXTURES(LINESOURCES(ILINE)%MIX_ID)%N_COMPONENTS
         
         ALLOCATE(nfs_LINE(N_COMP))

         DO IS = 1, N_COMP ! Loop on mixture components
            ! The species ID of the component
            S_ID = MIXTURES(LINESOURCES(ILINE)%MIX_ID)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            FRAC = MIXTURES(LINESOURCES(ILINE)%MIX_ID)%COMPONENTS(IS)%MOLFRAC
            BETA = 1./SQRT(2.*KB/M*LINESOURCES(ILINE)%TTRA) ! sqrt(M/(2*kB*T)), it's the Maxwellian std dev
         
            LINELENGTH = SQRT((LINESOURCES(ILINE)%X2-LINESOURCES(ILINE)%X1)**2 + (LINESOURCES(ILINE)%Y2-LINESOURCES(ILINE)%Y1)**2)

            NORMX = (LINESOURCES(ILINE)%Y2 - LINESOURCES(ILINE)%Y1)/LINELENGTH
            LINESOURCES(ILINE)%NORMX = NORMX
            NORMY = - (LINESOURCES(ILINE)%X2 - LINESOURCES(ILINE)%X1)/LINELENGTH
            LINESOURCES(ILINE)%NORMY = NORMY

            U_NORM = LINESOURCES(ILINE)%UX*NORMX + LINESOURCES(ILINE)%UY*NORMY ! Molecular speed ratio normal to boundary
            S_NORM = U_NORM*BETA
            LINESOURCES(ILINE)%S_NORM = S_NORM
            Snow   = S_NORM     ! temp variable

            FLUXLINESOURCE   = LINESOURCES(ILINE)%NRHO*FRAC/(BETA*2.*SQRT(PI)) * (EXP(-Snow**2) &
                                    + SQRT(PI)*Snow*(1.+ERF1(Snow)))      ! Tot number flux emitted

            IF (AXI) THEN
               AREA = 0.5*(ZMAX-ZMIN)*LINELENGTH*(LINESOURCES(ILINE)%Y1+LINESOURCES(ILINE)%Y2)
            ELSE
               AREA  = LINELENGTH*(ZMAX-ZMIN)
            END IF

            NtotINJECT  = FLUXLINESOURCE*AREA*DT/(FNUM*SPECIES(S_ID)%SPWT)         ! Tot num of particles to be injected

            nfs_LINE(IS)    = NtotINJECT/REAL(N_MPI_THREADS,KIND=8) ! Particles injected by each proc
            
            ! Check: if we are hypersonic and stuff exits domain print a warning
            IF (S_NORM .LE. -3.) THEN 
               CALL ONLYMASTERPRINT1(PROC_ID, '$$$ Warning! Hypersonic boundary! Almost no &
                                                   &particles will be emitted at line source.')
            END IF
         END DO

         CALL MOVE_ALLOC(nfs_LINE, LINESOURCES(ILINE)%nfs)
     
      END DO






      ! =====================================================
      ! Injection from unstructured grid boundaries
      ! WRITE(*,*) 'N emit tasks ', N_EMIT_TASKS
      DO ITASK = 1, N_EMIT_TASKS ! Loop on line sources

         N_COMP = MIXTURES(EMIT_TASKS(ITASK)%MIX_ID)%N_COMPONENTS
         
         ALLOCATE(TASK_NFS(N_COMP))

         IF (DIMS == 1) THEN
            NORMX = -U1D_GRID%EDGE_NORMAL(1, EMIT_TASKS(ITASK)%IFACE, EMIT_TASKS(ITASK)%IC)
            NORMY = 0
            NORMZ = 0

            AREA = (YMAX-YMIN)*(ZMAX-ZMIN)
         ELSE IF (DIMS == 2) THEN

            X1 = U2D_GRID%NODE_COORDS(1, EMIT_TASKS(ITASK)%IV1)
            Y1 = U2D_GRID%NODE_COORDS(2, EMIT_TASKS(ITASK)%IV1)
            X2 = U2D_GRID%NODE_COORDS(1, EMIT_TASKS(ITASK)%IV2)
            Y2 = U2D_GRID%NODE_COORDS(2, EMIT_TASKS(ITASK)%IV2)

            LINELENGTH = SQRT((X2-X1)**2 + (Y2-Y1)**2)

            NORMX = -U2D_GRID%EDGE_NORMAL(1, EMIT_TASKS(ITASK)%IFACE, EMIT_TASKS(ITASK)%IC)
            NORMY = -U2D_GRID%EDGE_NORMAL(2, EMIT_TASKS(ITASK)%IFACE, EMIT_TASKS(ITASK)%IC)
            NORMZ = 0

            IF (AXI) THEN
               AREA = (ZMAX-ZMIN)*0.5*(Y1+Y2)*LINELENGTH
            ELSE
               AREA = LINELENGTH*(ZMAX-ZMIN)
            END IF

         ELSE IF (DIMS == 3) THEN
            NORMX = -U3D_GRID%FACE_NORMAL(1, EMIT_TASKS(ITASK)%IFACE, EMIT_TASKS(ITASK)%IC)
            NORMY = -U3D_GRID%FACE_NORMAL(2, EMIT_TASKS(ITASK)%IFACE, EMIT_TASKS(ITASK)%IC)
            NORMZ = -U3D_GRID%FACE_NORMAL(3, EMIT_TASKS(ITASK)%IFACE, EMIT_TASKS(ITASK)%IC)
            
            AREA = U3D_GRID%FACE_AREA(EMIT_TASKS(ITASK)%IFACE, EMIT_TASKS(ITASK)%IC)
         END IF


         DO IS = 1, N_COMP ! Loop on mixture components
            ! The species ID of the component
            S_ID = MIXTURES(EMIT_TASKS(ITASK)%MIX_ID)%COMPONENTS(IS)%ID
            M = SPECIES(S_ID)%MOLECULAR_MASS
            FRAC = MIXTURES(EMIT_TASKS(ITASK)%MIX_ID)%COMPONENTS(IS)%MOLFRAC
            IF (EMIT_TASKS(ITASK)%TTRA == 0) THEN
               BETA = 0
            ELSE
               BETA = 1./SQRT(2.*KB/M*EMIT_TASKS(ITASK)%TTRA) ! sqrt(M/(2*kB*T)), it's the Maxwellian std dev
            END IF

            U_NORM = EMIT_TASKS(ITASK)%UX*NORMX + EMIT_TASKS(ITASK)%UY*NORMY + EMIT_TASKS(ITASK)%UZ*NORMZ ! Molecular speed ratio normal to boundary
            EMIT_TASKS(ITASK)%U_NORM = U_NORM
            S_NORM = U_NORM*BETA

            IF (EMIT_TASKS(ITASK)%TTRA == 0) THEN
               FLUXLINESOURCE = EMIT_TASKS(ITASK)%NRHO*FRAC*U_NORM      ! Tot number flux emitted
            ELSE
               FLUXLINESOURCE = EMIT_TASKS(ITASK)%NRHO*FRAC/(BETA*2.*SQRT(PI)) * (EXP(-S_NORM**2) &
                              + SQRT(PI)*S_NORM*(1.+ERF1(S_NORM)))      ! Tot number flux emitted
            END IF

            !!!!! KAPPA DISTRIBUTION !!!!! TODO: Obtain from Kappa distr func.
            ! IF (BOOL_KAPPA_DISTRIBUTION) THEN
            !    BETA = 1./SQRT(2.*KB/M*EMIT_TASKS(ITASK)%TTRA*(KAPPA_C-3./2.))
            !    U_NORM = EMIT_TASKS(ITASK)%UX*NORMX + EMIT_TASKS(ITASK)%UY*NORMY + EMIT_TASKS(ITASK)%UZ*NORMZ ! Molecular speed ratio normal to boundary
            !    EMIT_TASKS(ITASK)%U_NORM = U_NORM
            !    S_NORM = U_NORM*BETA

            !    INTEG1 = PI/2.
            !    INTEG2 = ATAN(S_NORM)
            !    DO I = 2, INT(KAPPA_C)
            !       INTEG1 = (2.*I-3.)/(2.*(I-1.))*INTEG1
            !       INTEG2 = S_NORM/(2.*(I-1.)*(1.+(S_NORM)**2)**(I-1.)) + (2.*I-3.)/(2.*(I-1.))*INTEG2
            !    END DO
               
            !    FLUXLINESOURCE = GAMMA(KAPPA_C)/GAMMA(KAPPA_C-1./2.)*EMIT_TASKS(ITASK)%NRHO*FRAC/(BETA*SQRT(PI)) * &
            !    ( 1./(KAPPA_C-1.)/2.*(1+S_NORM**2)**(-KAPPA_C+1.) + S_NORM*(INTEG1 +INTEG2))      ! Tot number flux emitted

            ! END IF


            NtotINJECT = FLUXLINESOURCE*AREA*DT/FNUM         ! Tot num of particles to be injected

            TASK_NFS(IS) = NtotINJECT/REAL(N_MPI_THREADS,KIND=8) ! Particles injected by each proc
            
            ! Check: if we are hypersonic and stuff exits domain print a warning
            IF (S_NORM .LE. -3.) THEN 
               CALL ONLYMASTERPRINT1(PROC_ID, '$$$ Warning! Hypersonic boundary! Almost no &
                                                   &particles will be emitted at boundary.')
            END IF
         END DO
         !WRITE(*,*) 'Task NFS', TASK_NFS
         CALL MOVE_ALLOC(TASK_NFS, EMIT_TASKS(ITASK)%NFS)

      END DO


   END SUBROUTINE INITINJECTION


   SUBROUTINE INITREACTIONS
   ! Here we precompute the TCE coefficients, with the input collision parameters

      IMPLICIT NONE

      REAL(KIND=8) :: OMEGA, TREF, SIGMA, ZETA, ETA, LAMBDA, EPS, MR, EA, M1, M2
      REAL(KIND=8) :: C1, C2, C3
      INTEGER      :: JR, R1_SP_ID, R2_SP_ID


      DO JR = 1, N_REACTIONS
         IF (REACTIONS(JR)%TYPE == TCE .AND. .NOT. REACTIONS(JR)%IS_CEX) THEN
            R1_SP_ID = REACTIONS(JR)%R1_SP_ID
            R2_SP_ID = REACTIONS(JR)%R2_SP_ID

            ZETA  = 0.5*(SPECIES(R1_SP_ID)%ROTDOF + SPECIES(R1_SP_ID)%VIBDOF + SPECIES(R2_SP_ID)%ROTDOF + SPECIES(R2_SP_ID)%VIBDOF)
            OMEGA = 0.5*(SPECIES(R1_SP_ID)%OMEGA + SPECIES(R2_SP_ID)%OMEGA)
            SIGMA = 0.25*PI*(SPECIES(R1_SP_ID)%DIAM + SPECIES(R2_SP_ID)%DIAM)**2
            TREF  = 0.5*(SPECIES(R1_SP_ID)%TREF + SPECIES(R2_SP_ID)%TREF)
            M1 = SPECIES(R1_SP_ID)%MOLECULAR_MASS
            M2 = SPECIES(R2_SP_ID)%MOLECULAR_MASS
            MR    = M1*M2/(M1+M2)
            LAMBDA = REACTIONS(JR)%A
            ETA = REACTIONS(JR)%N
            EA = REACTIONS(JR)%EA

            IF (R1_SP_ID == R2_SP_ID) THEN
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
         END IF
      END DO

   END SUBROUTINE INITREACTIONS


   SUBROUTINE READ_MCC_BACKGROUND_FILE(FILENAME)

      IMPLICIT NONE

      CHARACTER*512, INTENT(IN) :: FILENAME
      CHARACTER*80      :: line
      INTEGER, PARAMETER :: in5 = 5557
      INTEGER            :: ios
      INTEGER            :: ReasonEOF
      INTEGER            :: N_STR, SP_ID, STAT
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      ! Open input file for reading
      !OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(UNIT=in5, FILE=FILENAME, STATUS='OLD', IOSTAT=ios, ACCESS='STREAM', CONVERT='BIG_ENDIAN')
      ELSE
         OPEN(UNIT=in5, FILE=FILENAME, STATUS='OLD', IOSTAT=ios, FORM='FORMATTED')
      END IF

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, vtk file not found! ABORTING.')
      END IF

      ALLOCATE(MCC_BG_CELL_NRHO(N_SPECIES, NCELLS))

      ! ++++++++ Read until the end of file ++++++++
      IF (BOOL_BINARY_OUTPUT) THEN
         DO
            CALL SKIP_TO(in5, 'nrho_mean_', STAT)

            READ(in5, IOSTAT=ReasonEOF) line ! Read line
            IF (ReasonEOF < 0) EXIT ! End of file reached

            SP_ID = SPECIES_NAME_TO_ID(line)
            IF (SP_ID == -1) CALL ERROR_ABORT('Error! Species in MCC background vtk file not found.')

            WRITE(*,*) 'Found data for species ', SP_ID
            CALL SKIP_TO(in5, ACHAR(10), STAT)

            READ(in5, IOSTAT=ReasonEOF) MCC_BG_CELL_NRHO(SP_ID, :)
            IF (ReasonEOF < 0) EXIT ! End of file reached
         END DO ! Loop for reading input file
      ELSE
         DO
            READ(in5, IOSTAT=ReasonEOF) line ! Read line
            CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
            IF (STRARRAY(1)(1:10) == 'nrho_mean_') THEN
               SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1)(11:))
               IF (SP_ID == -1) CALL ERROR_ABORT('Error! Species in MCC background vtk file not found.')

               WRITE(*,*) 'Found data for species ', SP_ID
               READ(in5, IOSTAT=ReasonEOF) MCC_BG_CELL_NRHO(SP_ID, :)
               IF (ReasonEOF < 0) EXIT ! End of file reached
            END IF
         END DO
      END IF

      CLOSE(in5) ! Close input file

      BOOL_BG_DENSITY_FILE = .TRUE.

      WRITE(*,*) 'Done reading.'

   END SUBROUTINE READ_MCC_BACKGROUND_FILE


END MODULE initialization
