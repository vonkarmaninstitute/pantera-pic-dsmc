! PANTERA PIC-DSMC - A software for the simulation of rarefied gases
! and plasmas using particles
! Copyright (C) 2025 von Karman Institute for Fluid Dynamics (VKI)
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

PROGRAM PANTERA 

   USE particle
   USE mpi_common
   USE global
   USE screen
   USE tools
   USE initialization
   USE timecycle
   USE grid_and_partition
   USE fields
   USE washboard

   IMPLICIT NONE

   ! ========= Init MPI environment ========================
   CALL ONLYMASTERPRINT1(PROC_ID, '> INITIALIZING MPI...')
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_THREADS, ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_ID, ierr)
   CALL NEWTYPE ! Define new "particle" type for mpi
   CALL ONLYMASTERPRINT1(PROC_ID, '> MPI INIT DONE!')

   CALL ONLYMASTERPRINT1(PROC_ID, '> INITIALIZING PETSc...')
   CALL PETSC_INIT
   CALL ONLYMASTERPRINT1(PROC_ID, '> PETSc INIT DONE!')
   ! ========= Print header (important things first) =======
   CALL PRINTTITLE(PROC_ID)

   ! ========= Read input file and init variables ==========
   CALL ONLYMASTERPRINT1(PROC_ID, '> READING INPUT DATA...')
   CALL READINPUT          ! Read input file
   
   CALL INITVARIOUS        ! Initialize some additional variables
   CALL INITINJECTION      ! Initialize variables for injection
   CALL INITREACTIONS      ! Initialize variables for reactions
   CALL INITFIELDS         ! Initialize electromagnetic fields
   CALL COMPUTE_B_FIELD_FROM_SOLENOIDS

   ! ========= Initial particles seed ======================
   IF (RESTART_TIMESTEP > 0) THEN
      CALL READ_PARTICLES_FILE(RESTART_TIMESTEP)

      IF (DIMS == 1) THEN 
         CALL REASSIGN_PARTICLES_TO_CELLS_1D
      ELSE IF (DIMS == 2) THEN 
         CALL REASSIGN_PARTICLES_TO_CELLS_2D
      ELSE IF (DIMS == 3) THEN
         CALL REASSIGN_PARTICLES_TO_CELLS_3D
      END IF
      
      CALL ASSIGN_CELLS_TO_PROCS
      tID = RESTART_TIMESTEP
   ELSE
      tID = 0
      RESTART_TIMESTEP = 0
      CALL ASSIGN_CELLS_TO_PROCS
      CALL INITIAL_SEED
   END IF
   ! CALL DUMP_PARTICLES_SCREEN
   
   CALL EXCHANGE

   IF (PIC_TYPE .NE. NONE) CALL ASSEMBLE_POISSON

   ! ========= Time loop ===================================
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
   CALL ONLYMASTERPRINT1(PROC_ID, '> STARTING TIME LOOP...')
   !FLUSH(101)
   CALL TIME_LOOP

   ! CALL DUMP_PARTICLES_SCREEN

   ! ========== Close MPI environment ======================
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

   CALL TIMER_SUMMARY

   CALL ONLYMASTERPRINT1(PROC_ID, '> R O A R !')
   CALL PetscFinalize(ierr)
   CALL MPI_FINALIZE(ierr)

END PROGRAM PANTERA
