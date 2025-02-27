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

MODULE postprocess

   USE global
   USE screen
   USE grid_and_partition
   USE fields

   IMPLICIT NONE

   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GRID_AVG -> Adds timestep to cumulative average !!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GRID_AVG

      IMPLICIT NONE

      INTEGER                            :: JP, JC, JS, INDEX

      REAL(KIND=8) :: DBLE_AVG_CUMULATED, NUMPART, SAMPLEDOF
      REAL(KIND=8) :: CX, CY, CZ, C2, MASS, VOL, CFNUM

   
      INTEGER, DIMENSION(:), ALLOCATABLE      :: TIMESTEP_NP
   
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VY
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VZ

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VX2
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VY2
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_VZ2
   
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TTRX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TTRY
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TTRZ
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TTR

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_EROT
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_EVIB
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TROT
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_TVIB

      LOGICAL, DIMENSION(:), ALLOCATABLE      :: INTENSIVE_AVERAGE_ONE
      LOGICAL, DIMENSION(:), ALLOCATABLE      :: INTENSIVE_AVERAGE_TWO

      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: TIMESTEP_MOMENTS

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PHI

      INTEGER :: LENGTH

      LENGTH = NCELLS * N_SPECIES

      ALLOCATE(TIMESTEP_NP(LENGTH))

      ALLOCATE(TIMESTEP_VX(LENGTH))
      ALLOCATE(TIMESTEP_VY(LENGTH))
      ALLOCATE(TIMESTEP_VZ(LENGTH))

      ALLOCATE(TIMESTEP_VX2(LENGTH))
      ALLOCATE(TIMESTEP_VY2(LENGTH))
      ALLOCATE(TIMESTEP_VZ2(LENGTH))
      
      ALLOCATE(TIMESTEP_TTRX(LENGTH))
      ALLOCATE(TIMESTEP_TTRY(LENGTH))
      ALLOCATE(TIMESTEP_TTRZ(LENGTH))
      ALLOCATE(TIMESTEP_TTR(LENGTH))
      
      ALLOCATE(TIMESTEP_EROT(LENGTH))
      ALLOCATE(TIMESTEP_EVIB(LENGTH))

      ALLOCATE(TIMESTEP_TROT(LENGTH))
      ALLOCATE(TIMESTEP_TVIB(LENGTH))

      ALLOCATE(INTENSIVE_AVERAGE_ONE(LENGTH))
      ALLOCATE(INTENSIVE_AVERAGE_TWO(LENGTH))

      TIMESTEP_NP = 0

      TIMESTEP_VX = 0
      TIMESTEP_VY = 0
      TIMESTEP_VZ = 0

      TIMESTEP_VX2 = 0
      TIMESTEP_VY2 = 0
      TIMESTEP_VZ2 = 0

      TIMESTEP_TTRX = 0
      TIMESTEP_TTRY = 0
      TIMESTEP_TTRZ = 0
      TIMESTEP_TTR = 0

      TIMESTEP_EROT = 0
      TIMESTEP_EVIB = 0

      TIMESTEP_TROT = 0
      TIMESTEP_TVIB = 0

      INTENSIVE_AVERAGE_ONE = .FALSE.
      INTENSIVE_AVERAGE_TWO = .FALSE.

      IF (BOOL_DUMP_MOMENTS) THEN
         ALLOCATE(TIMESTEP_MOMENTS(LENGTH,33))
         TIMESTEP_MOMENTS = 0
      END IF


      ! Compute average values for this timestep on this process

      ! Number of particles
      DO JP = 1, NP_PROC
         JC = particles(JP)%IC
         JS = particles(JP)%S_ID
         INDEX = JC+NCELLS*(JS-1)
         TIMESTEP_NP(INDEX) = TIMESTEP_NP(INDEX) + 1
      END DO

      ! Velocity
      DO JP = 1, NP_PROC
         JC = particles(JP)%IC
         JS = particles(JP)%S_ID
         INDEX = JC+NCELLS*(JS-1)
         TIMESTEP_VX(INDEX) = TIMESTEP_VX(INDEX) + particles(JP)%VX
         TIMESTEP_VY(INDEX) = TIMESTEP_VY(INDEX) + particles(JP)%VY
         TIMESTEP_VZ(INDEX) = TIMESTEP_VZ(INDEX) + particles(JP)%VZ
      END DO


      DO INDEX = 1, NCELLS*N_SPECIES
         IF (TIMESTEP_NP(INDEX) .GT. 0) THEN
            INTENSIVE_AVERAGE_ONE(INDEX) = .TRUE.
            TIMESTEP_VX(INDEX) = TIMESTEP_VX(INDEX) / DBLE(TIMESTEP_NP(INDEX))
            TIMESTEP_VY(INDEX) = TIMESTEP_VY(INDEX) / DBLE(TIMESTEP_NP(INDEX))
            TIMESTEP_VZ(INDEX) = TIMESTEP_VZ(INDEX) / DBLE(TIMESTEP_NP(INDEX))
         END IF
      END DO

      ! Translational Temperature and internal energy
      DO JP = 1, NP_PROC
         JC = particles(JP)%IC
         JS = particles(JP)%S_ID
         INDEX = JC+NCELLS*(JS-1)
         TIMESTEP_VX2(INDEX) = TIMESTEP_VX2(INDEX) + particles(JP)%VX*particles(JP)%VX
         TIMESTEP_VY2(INDEX) = TIMESTEP_VY2(INDEX) + particles(JP)%VY*particles(JP)%VY
         TIMESTEP_VZ2(INDEX) = TIMESTEP_VZ2(INDEX) + particles(JP)%VZ*particles(JP)%VZ

         TIMESTEP_EROT(INDEX) = TIMESTEP_EROT(INDEX) + particles(JP)%EROT
         TIMESTEP_EVIB(INDEX) = TIMESTEP_EVIB(INDEX) + particles(JP)%EVIB
      END DO


      IF (BOOL_DUMP_MOMENTS) THEN
         ! Higher order moments
         DO JP = 1, NP_PROC
            JC = particles(JP)%IC
            JS = particles(JP)%S_ID
            INDEX = JC+NCELLS*(JS-1)

            CX = particles(JP)%VX - TIMESTEP_VX(INDEX)
            CY = particles(JP)%VY - TIMESTEP_VY(INDEX)
            CZ = particles(JP)%VZ - TIMESTEP_VZ(INDEX)

            C2 = CX*CX + CY*CY + CZ*CZ
            TIMESTEP_MOMENTS(INDEX,5)  = TIMESTEP_MOMENTS(INDEX,5)  + CX*CX ! Pxx
            TIMESTEP_MOMENTS(INDEX,6)  = TIMESTEP_MOMENTS(INDEX,6)  + CX*CY ! Pxy
            TIMESTEP_MOMENTS(INDEX,7)  = TIMESTEP_MOMENTS(INDEX,7)  + CX*CZ ! Pxz
            TIMESTEP_MOMENTS(INDEX,8)  = TIMESTEP_MOMENTS(INDEX,8)  + CY*CY ! Pyy
            TIMESTEP_MOMENTS(INDEX,9)  = TIMESTEP_MOMENTS(INDEX,9)  + CY*CZ ! Pyz
            TIMESTEP_MOMENTS(INDEX,10) = TIMESTEP_MOMENTS(INDEX,10) + CZ*CZ ! Pzz

            TIMESTEP_MOMENTS(INDEX,11) = TIMESTEP_MOMENTS(INDEX,11) + CX*C2 ! qx
            TIMESTEP_MOMENTS(INDEX,12) = TIMESTEP_MOMENTS(INDEX,12) + CY*C2 ! qy
            TIMESTEP_MOMENTS(INDEX,13) = TIMESTEP_MOMENTS(INDEX,13) + CZ*C2 ! qz

            TIMESTEP_MOMENTS(INDEX,14) = TIMESTEP_MOMENTS(INDEX,14) + CX*CX*CX ! Qxxx
            TIMESTEP_MOMENTS(INDEX,15) = TIMESTEP_MOMENTS(INDEX,15) + CX*CX*CY ! Qxxy
            TIMESTEP_MOMENTS(INDEX,16) = TIMESTEP_MOMENTS(INDEX,16) + CX*CY*CY ! Qxyy
            TIMESTEP_MOMENTS(INDEX,17) = TIMESTEP_MOMENTS(INDEX,17) + CY*CY*CY ! Qyyy
            TIMESTEP_MOMENTS(INDEX,18) = TIMESTEP_MOMENTS(INDEX,18) + CY*CY*CZ ! Qyyz
            TIMESTEP_MOMENTS(INDEX,19) = TIMESTEP_MOMENTS(INDEX,19) + CY*CZ*CZ ! Qyzz
            TIMESTEP_MOMENTS(INDEX,20) = TIMESTEP_MOMENTS(INDEX,20) + CZ*CZ*CZ ! Qzzz
            TIMESTEP_MOMENTS(INDEX,21) = TIMESTEP_MOMENTS(INDEX,21) + CX*CX*CZ ! Qxxz
            TIMESTEP_MOMENTS(INDEX,22) = TIMESTEP_MOMENTS(INDEX,22) + CX*CZ*CZ ! Qxzz
            TIMESTEP_MOMENTS(INDEX,23) = TIMESTEP_MOMENTS(INDEX,23) + CX*CY*CZ ! Qxyz

            TIMESTEP_MOMENTS(INDEX,24) = TIMESTEP_MOMENTS(INDEX,24) + C2*C2 ! Riijj

            TIMESTEP_MOMENTS(INDEX,25) = TIMESTEP_MOMENTS(INDEX,25) + CX*CX*C2 ! Rxxjj
            TIMESTEP_MOMENTS(INDEX,26) = TIMESTEP_MOMENTS(INDEX,26) + CX*CY*C2 ! Rxyjj
            TIMESTEP_MOMENTS(INDEX,27) = TIMESTEP_MOMENTS(INDEX,27) + CX*CZ*C2 ! Rxzjj
            TIMESTEP_MOMENTS(INDEX,28) = TIMESTEP_MOMENTS(INDEX,28) + CY*CY*C2 ! Ryyjj
            TIMESTEP_MOMENTS(INDEX,29) = TIMESTEP_MOMENTS(INDEX,29) + CY*CZ*C2 ! Ryzjj
            TIMESTEP_MOMENTS(INDEX,30) = TIMESTEP_MOMENTS(INDEX,30) + CZ*CZ*C2 ! Rzzjj

            TIMESTEP_MOMENTS(INDEX,31) = TIMESTEP_MOMENTS(INDEX,31) + CX*C2*C2 ! Sxiijj
            TIMESTEP_MOMENTS(INDEX,32) = TIMESTEP_MOMENTS(INDEX,32) + CY*C2*C2 ! Syiijj
            TIMESTEP_MOMENTS(INDEX,33) = TIMESTEP_MOMENTS(INDEX,33) + CZ*C2*C2 ! Sziijj

         END DO

         DO JC = 1, NCELLS

            IF (GRID_TYPE == RECTILINEAR_UNIFORM .AND. .NOT. AXI) THEN
               VOL = CELL_VOL
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
               VOL = U1D_GRID%CELL_VOLUMES(JC)
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
               VOL = U2D_GRID%CELL_VOLUMES(JC)
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
               VOL = U3D_GRID%CELL_VOLUMES(JC)
            END IF

            IF (BOOL_RADIAL_WEIGHTING) THEN
               CFNUM = CELL_FNUM(particles(JP)%IC)
            ELSE
               CFNUM = FNUM
            END IF

            DO JS = 1, N_SPECIES
               INDEX = JC+NCELLS*(JS-1)
               MASS = SPECIES(JS)%MOLECULAR_MASS
               

               ! rho
               TIMESTEP_MOMENTS(INDEX,1) = MASS*CFNUM/VOL*TIMESTEP_NP(INDEX)
               ! Ux, Uy, Uz
               TIMESTEP_MOMENTS(INDEX,2) = TIMESTEP_VX(INDEX)
               TIMESTEP_MOMENTS(INDEX,3) = TIMESTEP_VY(INDEX)
               TIMESTEP_MOMENTS(INDEX,4) = TIMESTEP_VZ(INDEX)
               ! Higher order moments
               TIMESTEP_MOMENTS(INDEX,5:33) = TIMESTEP_MOMENTS(INDEX,5:33)*MASS*CFNUM/VOL
            END DO
         END DO

      END IF


      DO JC = 1, NCELLS
         DO JS = 1, N_SPECIES
            INDEX = JC+NCELLS*(JS-1)
            IF (TIMESTEP_NP(INDEX) .GT. 0) THEN


               IF (SPECIES(JS)%ROTDOF .EQ. 0) THEN
                  TIMESTEP_TROT(INDEX) = 0.0
               ELSE
                  TIMESTEP_TROT(INDEX) = 2. / SPECIES(JS)%ROTDOF / KB * (TIMESTEP_EROT(INDEX) / DBLE(TIMESTEP_NP(INDEX)) )
               END IF

               IF (SPECIES(JS)%VIBDOF .EQ. 0) THEN
                  TIMESTEP_TVIB(INDEX) = 0.0
               ELSE
                  TIMESTEP_TVIB(INDEX) = 2. / SPECIES(JS)%VIBDOF / KB * (TIMESTEP_EVIB(INDEX) / DBLE(TIMESTEP_NP(INDEX)) )
               END IF




               IF (TIMESTEP_NP(INDEX) .GT. 1) THEN
                  INTENSIVE_AVERAGE_TWO(INDEX) = .TRUE.
                  NUMPART = DBLE(TIMESTEP_NP(INDEX))
                  ! Old method: SAMPLEDOF = NUMPART. Biased when sample is small.
                  ! This implementation (commented) uses a less biased estimation of the standard deviation.
                  ! SAMPLEDOF = NUMPART - 1.5 + 1./(8*(NUMPART-1))
                  SAMPLEDOF = NUMPART
                  TIMESTEP_TTRX(INDEX) = SPECIES(JS)%MOLECULAR_MASS / KB &
                     * (TIMESTEP_VX2(INDEX) - NUMPART*TIMESTEP_VX(INDEX)*TIMESTEP_VX(INDEX)) / SAMPLEDOF
                  TIMESTEP_TTRY(INDEX) = SPECIES(JS)%MOLECULAR_MASS / KB &
                     * (TIMESTEP_VY2(INDEX) - NUMPART*TIMESTEP_VY(INDEX)*TIMESTEP_VY(INDEX)) / SAMPLEDOF
                  TIMESTEP_TTRZ(INDEX) = SPECIES(JS)%MOLECULAR_MASS / KB &
                     * (TIMESTEP_VZ2(INDEX) - NUMPART*TIMESTEP_VZ(INDEX)*TIMESTEP_VZ(INDEX)) / SAMPLEDOF
               END IF


            END IF
         END DO
      END DO



      ! Collect data from all the processes
      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_NP, NCELLS*N_SPECIES, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VX,  NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VY,  NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VZ,  NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_TTRX,  NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_TTRY,  NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_TTRZ,  NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_TROT,  NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_TVIB,  NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  INTENSIVE_AVERAGE_ONE,  NCELLS*N_SPECIES, MPI_LOGICAL, MPI_LOR, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  INTENSIVE_AVERAGE_TWO,  NCELLS*N_SPECIES, MPI_LOGICAL, MPI_LOR, 0, MPI_COMM_WORLD, ierr)
         IF (BOOL_DUMP_MOMENTS) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,TIMESTEP_MOMENTS,NCELLS*N_SPECIES*33,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         END IF
         !TALLY = 0
         !DO INDEX = 1, NCELLS*N_SPECIES
         !  TALLY = TALLY + TIMESTEP_NP(INDEX)
         !END DO
         !WRITE(*,*) 'Total in postprocess', TALLY
      ELSE
         CALL MPI_REDUCE(TIMESTEP_NP,  TIMESTEP_NP,  NCELLS*N_SPECIES, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VX,   TIMESTEP_VX,   NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VY,   TIMESTEP_VY,   NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_VZ,   TIMESTEP_VZ,   NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_TTRX,   TIMESTEP_TTRX,   NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_TTRY,   TIMESTEP_TTRY,   NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_TTRZ,   TIMESTEP_TTRZ,   NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_TROT,   TIMESTEP_TROT,   NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TIMESTEP_TVIB,   TIMESTEP_TVIB,   NCELLS*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(INTENSIVE_AVERAGE_ONE, INTENSIVE_AVERAGE_ONE, NCELLS*N_SPECIES, MPI_LOGICAL, MPI_LOR,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(INTENSIVE_AVERAGE_TWO, INTENSIVE_AVERAGE_TWO, NCELLS*N_SPECIES, MPI_LOGICAL, MPI_LOR,0,MPI_COMM_WORLD,ierr)
         IF (BOOL_DUMP_MOMENTS) THEN
         CALL MPI_REDUCE(TIMESTEP_MOMENTS,TIMESTEP_MOMENTS,NCELLS*N_SPECIES*33,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         END IF
      END IF

      IF (PIC_TYPE .NE. NONE) TIMESTEP_PHI = PHI_FIELD

      ! Add to cumulated average
      DBLE_AVG_CUMULATED = DBLE(AVG_CUMULATED)
      AVG_NP =   (AVG_NP*DBLE_AVG_CUMULATED + DBLE(TIMESTEP_NP))/(AVG_CUMULATED + 1.)
      
      IF (PIC_TYPE .NE. NONE) AVG_PHI = (AVG_PHI*DBLE_AVG_CUMULATED + TIMESTEP_PHI)/(AVG_CUMULATED + 1.)

      AVG_CUMULATED = AVG_CUMULATED + 1

      DO INDEX = 1, NCELLS*N_SPECIES
         IF (INTENSIVE_AVERAGE_ONE(INDEX)) THEN
         AVG_VX(INDEX) = (AVG_VX(INDEX)*AVG_CUMULATED_INTENSIVE_ONE(INDEX) + TIMESTEP_VX(INDEX)) & 
         /(AVG_CUMULATED_INTENSIVE_ONE(INDEX) + 1.)
         AVG_VY(INDEX) = (AVG_VY(INDEX)*AVG_CUMULATED_INTENSIVE_ONE(INDEX) + TIMESTEP_VY(INDEX)) & 
         /(AVG_CUMULATED_INTENSIVE_ONE(INDEX) + 1.)
         AVG_VZ(INDEX) = (AVG_VZ(INDEX)*AVG_CUMULATED_INTENSIVE_ONE(INDEX) + TIMESTEP_VZ(INDEX)) & 
         /(AVG_CUMULATED_INTENSIVE_ONE(INDEX) + 1.)

         AVG_TROT(INDEX) = (AVG_TROT(INDEX)*AVG_CUMULATED_INTENSIVE_ONE(INDEX) + TIMESTEP_TROT(INDEX)) & 
         /(AVG_CUMULATED_INTENSIVE_ONE(INDEX) + 1.)
         AVG_TVIB(INDEX) = (AVG_TVIB(INDEX)*AVG_CUMULATED_INTENSIVE_ONE(INDEX) + TIMESTEP_TVIB(INDEX)) & 
         /(AVG_CUMULATED_INTENSIVE_ONE(INDEX) + 1.)

         AVG_CUMULATED_INTENSIVE_ONE(INDEX) = AVG_CUMULATED_INTENSIVE_ONE(INDEX) + 1
         END IF
         IF (INTENSIVE_AVERAGE_TWO(INDEX)) THEN
         AVG_TTRX(INDEX) = (AVG_TTRX(INDEX)*AVG_CUMULATED_INTENSIVE_TWO(INDEX) + TIMESTEP_TTRX(INDEX)) &
         /(AVG_CUMULATED_INTENSIVE_TWO(INDEX) + 1.)
         AVG_TTRY(INDEX) = (AVG_TTRY(INDEX)*AVG_CUMULATED_INTENSIVE_TWO(INDEX) + TIMESTEP_TTRY(INDEX)) &
         /(AVG_CUMULATED_INTENSIVE_TWO(INDEX) + 1.)
         AVG_TTRZ(INDEX) = (AVG_TTRZ(INDEX)*AVG_CUMULATED_INTENSIVE_TWO(INDEX) + TIMESTEP_TTRZ(INDEX)) &
         /(AVG_CUMULATED_INTENSIVE_TWO(INDEX) + 1.)

         IF (BOOL_DUMP_MOMENTS) THEN
            AVG_MOMENTS(INDEX,:) = (AVG_MOMENTS(INDEX,:)*AVG_CUMULATED_INTENSIVE_TWO(INDEX) + TIMESTEP_MOMENTS(INDEX,:)) &
            /(AVG_CUMULATED_INTENSIVE_TWO(INDEX) + 1.)
         END IF

         AVG_CUMULATED_INTENSIVE_TWO(INDEX) = AVG_CUMULATED_INTENSIVE_TWO(INDEX) + 1
         END IF
      END DO

      AVG_TTR = (AVG_TTRX + AVG_TTRY + AVG_TTRZ) / 3.0



         
      DEALLOCATE(TIMESTEP_NP)

      DEALLOCATE(TIMESTEP_VX)
      DEALLOCATE(TIMESTEP_VY)
      DEALLOCATE(TIMESTEP_VZ)

      DEALLOCATE(TIMESTEP_VX2)
      DEALLOCATE(TIMESTEP_VY2)
      DEALLOCATE(TIMESTEP_VZ2)
      
      DEALLOCATE(TIMESTEP_TTRX)
      DEALLOCATE(TIMESTEP_TTRY)
      DEALLOCATE(TIMESTEP_TTRZ)
      DEALLOCATE(TIMESTEP_TTR)
      
      DEALLOCATE(TIMESTEP_EROT)
      DEALLOCATE(TIMESTEP_EVIB)

      DEALLOCATE(TIMESTEP_TROT)
      DEALLOCATE(TIMESTEP_TVIB)

      DEALLOCATE(INTENSIVE_AVERAGE_ONE)
      DEALLOCATE(INTENSIVE_AVERAGE_TWO)

      IF (BOOL_DUMP_MOMENTS) DEALLOCATE(TIMESTEP_MOMENTS)

      IF (PIC_TYPE .NE. NONE) DEALLOCATE(TIMESTEP_PHI)
      
   END SUBROUTINE GRID_AVG



   FUNCTION ITOA(I) RESULT(RES)
      CHARACTER(:), ALLOCATABLE :: RES
      INTEGER, INTENT(IN) :: I
      CHARACTER(RANGE(I)+2) :: TMP
      WRITE(TMP,'(I0)') I
      RES = TRIM(TMP)
   END FUNCTION



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GRID_SAVE -> Saves cumulated average !!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE GRID_SAVE
         
      IMPLICIT NONE

      CHARACTER*256                      :: string, file_name
      CHARACTER*7, DIMENSION(33)         :: MOMENT_STRING 

      REAL(KIND=8), DIMENSION(NX+1)      :: XNODES
      REAL(KIND=8), DIMENSION(NY+1)      :: YNODES

      INTEGER                            :: I, JS, FIRST, LAST, JPROC, MOM

      INTEGER, DIMENSION(:), ALLOCATABLE :: CELL_PROC_ID


      MOMENT_STRING = ['rho_   ', &
                        'Ux_    ', 'Uy_    ', 'Uz_    ', &
                        'Pxx_   ', 'Pxy_   ', 'Pxz_   ', 'Pyy_   ', 'Pyz_   ', 'Pzz_   ', &
                        'qx_    ', 'qy_    ', 'qz_    ', &
                        'Qxxx_  ', 'Qxxy_  ', 'Qxyy_  ', 'Qyyy_  ', 'Qyyz_  ', &
                        'Qyzz_  ', 'Qzzz_  ', 'Qxxz_  ', 'Qxzz_  ', 'Qxyz_  ', &
                        'Riijj_ ', 'Rxxjj_ ', 'Rxyjj_ ', 'Rxzjj_ ', 'Ryyjj_ ', 'Ryzjj_ ', 'Rzzjj_ ', &
                        'Sxiijj_', 'Syiijj_', 'Sziijj_']

      IF (PROC_ID .EQ. 0) THEN
         IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
            XNODES = 0
            DO I = 0, NX
               XNODES(i+1) = XMIN + (XMAX-XMIN)/NX*i 
            END DO

            YNODES = 0
            DO I = 0, NY
               YNODES(i+1) = YMIN + (YMAX-YMIN)/NY*i
            END DO
         ELSE IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) THEN
            XNODES = XCOORD
            YNODES = YCOORD
         END IF

         ALLOCATE(CELL_PROC_ID(NCELLS))
         DO I = 1, NCELLS
            CALL PROC_FROM_CELL(i, JPROC)
            CELL_PROC_ID(i) = JPROC
         END DO

         ! Write to file.

         !WRITE(file_name,'(A, I0, A)') '/media/pietro/Storage/panteradumps/dsmc_flowfield_', tID, '.vtk'
         WRITE(file_name,'(A, A, I0, A)') TRIM(ADJUSTL(FLOWFIELD_SAVE_PATH)), 'dsmc_flowfield_', tID, '.vtk'

         IF (BOOL_BINARY_OUTPUT) THEN
            OPEN(54321, FILE=file_name, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='NEW', CONVERT='BIG_ENDIAN')

            WRITE(54321) '# vtk DataFile Version 3.0'//ACHAR(10)
            WRITE(54321) 'vtk output'//ACHAR(10)
            WRITE(54321) 'BINARY'//ACHAR(10)


            IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
               WRITE(54321) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

               WRITE(54321) 'POINTS '//ITOA(U1D_GRID%NUM_NODES)//' double'//ACHAR(10)
               DO I = 1, U1D_GRID%NUM_NODES
                  WRITE(54321) U1D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321) 'CELLS '//ITOA(U1D_GRID%NUM_CELLS)//' '//ITOA(3*U1D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U1D_GRID%NUM_CELLS
                  WRITE(54321) 2, (U1D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321) 'CELL_TYPES '//ITOA(U1D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U1D_GRID%NUM_CELLS
                  WRITE(54321) 3
               END DO
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
               WRITE(54321) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

               WRITE(54321) 'POINTS '//ITOA(U2D_GRID%NUM_NODES)//' double'//ACHAR(10)
               DO I = 1, U2D_GRID%NUM_NODES
                  WRITE(54321) U2D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321) 'CELLS '//ITOA(U2D_GRID%NUM_CELLS)//' '//ITOA(4*U2D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U2D_GRID%NUM_CELLS
                  WRITE(54321) 3, (U2D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321) 'CELL_TYPES '//ITOA(U2D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U2D_GRID%NUM_CELLS
                  WRITE(54321) 5
               END DO
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
               WRITE(54321) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

               WRITE(54321) 'POINTS '//ITOA(U3D_GRID%NUM_NODES)//' double'//ACHAR(10)
               DO I = 1, U3D_GRID%NUM_NODES
                  WRITE(54321) U3D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321) 'CELLS '//ITOA(U3D_GRID%NUM_CELLS)//' '//ITOA(5*U3D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U3D_GRID%NUM_CELLS
                  WRITE(54321) 4, (U3D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321) 'CELL_TYPES '//ITOA(U3D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U3D_GRID%NUM_CELLS
                  WRITE(54321) 10
               END DO
            ELSE IF (DIMS == 1) THEN
               WRITE(54321) 'DATASET RECTILINEAR_GRID'//ACHAR(10)
               
               WRITE(54321) 'DIMENSIONS '//ITOA(NX+1)//' '//ITOA(1)//' '//ITOA(1)//ACHAR(10)

               WRITE(54321) 'X_COORDINATES '//ITOA(NX+1)//' double'//ACHAR(10)
               WRITE(54321) XNODES, ACHAR(10)

               WRITE(54321) 'Y_COORDINATES '//ITOA(1)//' double'//ACHAR(10)
               WRITE(54321) 0.d0, ACHAR(10)

               WRITE(54321) 'Z_COORDINATES '//ITOA(1)//' double'//ACHAR(10)
               WRITE(54321) 0.d0, ACHAR(10)
            ELSE IF (DIMS == 2) THEN
               WRITE(54321) 'DATASET RECTILINEAR_GRID'//ACHAR(10)
               
               WRITE(54321) 'DIMENSIONS '//ITOA(NX+1)//' '//ITOA(NY+1)//' '//ITOA(1)//ACHAR(10)

               WRITE(54321) 'X_COORDINATES '//ITOA(NX+1)//' double'//ACHAR(10)
               WRITE(54321) XNODES, ACHAR(10)

               WRITE(54321) 'Y_COORDINATES '//ITOA(NY+1)//' double'//ACHAR(10)
               WRITE(54321) YNODES, ACHAR(10)

               WRITE(54321) 'Z_COORDINATES '//ITOA(1)//' double'//ACHAR(10)
               WRITE(54321) 0.d0, ACHAR(10)
            END IF

            WRITE(54321) 'CELL_DATA '//ITOA(NCELLS)//ACHAR(10)
            IF (BOOL_DUMP_MOMENTS) THEN
               WRITE(54321) 'FIELD FieldData '//ITOA( (13+33)*N_SPECIES+1 )//ACHAR(10)
            ELSE
               WRITE(54321) 'FIELD FieldData '//ITOA( 13*N_SPECIES+1 )//ACHAR(10)
            END IF


            ! Write per-cell value
            WRITE(54321) 'PROC_ID '//ITOA(1)//' '//ITOA(NCELLS)//' integer'//ACHAR(10)
            WRITE(54321) CELL_PROC_ID, ACHAR(10)

            ! Write per-cell, per-species values
            DO JS = 1, N_SPECIES
               FIRST = 1 + (JS-1)*NCELLS
               LAST  = JS*NCELLS
            
               WRITE(string, *) 'number_particles_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_NP(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'nrho_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               IF (GRID_TYPE == RECTILINEAR_NONUNIFORM .OR. GRID_TYPE == UNSTRUCTURED .OR. AXI) THEN
                  IF (BOOL_RADIAL_WEIGHTING) THEN
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) THEN
                        WRITE(54321) CELL_FNUM*AVG_NP(FIRST:LAST)/CELL_VOLUMES, ACHAR(10)
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
                        WRITE(54321) CELL_FNUM*AVG_NP(FIRST:LAST)/U1D_GRID%CELL_VOLUMES, ACHAR(10)
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
                        WRITE(54321) CELL_FNUM*AVG_NP(FIRST:LAST)/U2D_GRID%CELL_VOLUMES, ACHAR(10)
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
                        WRITE(54321) CELL_FNUM*AVG_NP(FIRST:LAST)/U3D_GRID%CELL_VOLUMES, ACHAR(10)
                     END IF
                  ELSE
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) THEN
                        WRITE(54321) FNUM*AVG_NP(FIRST:LAST)/CELL_VOLUMES, ACHAR(10)
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
                        WRITE(54321) FNUM*AVG_NP(FIRST:LAST)/U1D_GRID%CELL_VOLUMES, ACHAR(10)
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
                        WRITE(54321) FNUM*AVG_NP(FIRST:LAST)/U2D_GRID%CELL_VOLUMES, ACHAR(10)
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
                        WRITE(54321) FNUM*AVG_NP(FIRST:LAST)/U3D_GRID%CELL_VOLUMES, ACHAR(10)
                     END IF
                  END IF
               ELSE IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
                  WRITE(54321) FNUM*AVG_NP(FIRST:LAST)/CELL_VOL, ACHAR(10)
               END IF

               WRITE(string, *) 'vx_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_VX(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'vy_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_VY(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'vz_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_VZ(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'Ttrx_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_TTRX(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'Ttry_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_TTRY(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'Ttrz_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_TTRZ(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'Ttr_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_TTR(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'Trot_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_TROT(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'Tvib_mean_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_TVIB(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'num_cumulated_one_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' integer'//ACHAR(10)
               WRITE(54321) AVG_CUMULATED_INTENSIVE_ONE(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'num_cumulated_two_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NCELLS)//' integer'//ACHAR(10)
               WRITE(54321) AVG_CUMULATED_INTENSIVE_TWO(FIRST:LAST), ACHAR(10)

            
               IF (BOOL_DUMP_MOMENTS) THEN
                  ! Now dump the 33 moments
                  DO MOM = 1, 33
                     WRITE(54321) TRIM(MOMENT_STRING(MOM))//SPECIES(JS)%NAME//' '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
                     WRITE(54321) AVG_MOMENTS(FIRST:LAST, MOM), ACHAR(10)
                  END DO
               END IF
         
         
            END DO

            IF (PIC_TYPE .NE. NONE) THEN
               IF (GRID_TYPE == UNSTRUCTURED) THEN
                  IF (PIC_TYPE == EXPLICITLIMITED) THEN
                     WRITE(54321) 'FIELD FieldData '//ITOA(4)//ACHAR(10)
                  ELSE
                     WRITE(54321) 'FIELD FieldData '//ITOA(3)//ACHAR(10)
                  END IF
                  WRITE(54321) 'E_X '//ITOA(1)//' '//ITOA( NCELLS )//' double'//ACHAR(10)
                  WRITE(54321) E_FIELD(1,:,:), ACHAR(10)

                  WRITE(54321) 'E_Y '//ITOA(1)//' '//ITOA( NCELLS )//' double'//ACHAR(10)
                  WRITE(54321) E_FIELD(2,:,:), ACHAR(10)

                  WRITE(54321) 'E_Z '//ITOA(1)//' '//ITOA( NCELLS )//' double'//ACHAR(10)
                  WRITE(54321) E_FIELD(3,:,:), ACHAR(10)

                  IF (PIC_TYPE == EXPLICITLIMITED) THEN
                     WRITE(54321) 'DXLDRATIO '//ITOA(1)//' '//ITOA( NCELLS )//' double'//ACHAR(10)
                     WRITE(54321) DXLDRATIO, ACHAR(10)
                  END IF


                  WRITE(54321) 'POINT_DATA '//ITOA( NNODES )//ACHAR(10)
                  IF (PIC_TYPE == HYBRID) THEN
                     WRITE(54321) 'FIELD FieldData '//ITOA(7)//ACHAR(10)
                  ELSE
                     WRITE(54321) 'FIELD FieldData '//ITOA(6)//ACHAR(10)
                  END IF
               
                  WRITE(54321) 'QRHO '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) RHS, ACHAR(10)

                  WRITE(54321) 'PHI '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) AVG_PHI, ACHAR(10)

                  WRITE(54321) 'PHIBAR '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) PHIBAR_FIELD, ACHAR(10)

                  WRITE(54321) 'B_X '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) B_FIELD(1,:,:), ACHAR(10)

                  WRITE(54321) 'B_Y '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) B_FIELD(2,:,:), ACHAR(10)

                  WRITE(54321) 'B_Z '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) B_FIELD(3,:,:), ACHAR(10)

                  IF (PIC_TYPE == HYBRID) THEN
                     WRITE(54321) 'nrho_e_FLUID '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                     WRITE(54321) BOLTZ_NRHOE, ACHAR(10)
                  END IF
               ELSE
                  WRITE(54321) 'POINT_DATA '//ITOA( NNODES )//ACHAR(10)
                  WRITE(54321) 'FIELD FieldData '//ITOA(5)//ACHAR(10)
               
                  WRITE(54321) 'QRHO '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) RHS, ACHAR(10)

                  WRITE(54321) 'PHI '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) AVG_PHI, ACHAR(10)

                  WRITE(54321) 'E_X '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) E_FIELD(1,:,:), ACHAR(10)

                  WRITE(54321) 'E_Y '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) E_FIELD(2,:,:), ACHAR(10)

                  WRITE(54321) 'E_Z '//ITOA(1)//' '//ITOA( NNODES )//' double'//ACHAR(10)
                  WRITE(54321) E_FIELD(3,:,:), ACHAR(10)
               END IF
            END IF

            CLOSE(54321)
            
         ELSE  ! Write ASCII output.
            OPEN(54321, FILE=file_name, ACCESS='SEQUENTIAL', FORM='FORMATTED', STATUS='NEW')

            WRITE(54321,'(A)') '# vtk DataFile Version 3.0'
            WRITE(54321,'(A)') 'vtk output'
            WRITE(54321,'(A)') 'ASCII'

            IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
               WRITE(54321,'(A)') 'DATASET UNSTRUCTURED_GRID'
               
               WRITE(54321,'(A,I10,A7)') 'POINTS', U1D_GRID%NUM_NODES, 'double'
               DO I = 1, U1D_GRID%NUM_NODES
                  WRITE(54321,*) U1D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321,'(A,I10,I10)') 'CELLS', U1D_GRID%NUM_CELLS, 3*U1D_GRID%NUM_CELLS 
               DO I = 1, U1D_GRID%NUM_CELLS
                  WRITE(54321,*) 2, (U1D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321,'(A,I10)') 'CELL_TYPES', U1D_GRID%NUM_CELLS
               DO I = 1, U1D_GRID%NUM_CELLS
                  WRITE(54321,*) 3
               END DO
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
               WRITE(54321,'(A)') 'DATASET UNSTRUCTURED_GRID'
               
               WRITE(54321,'(A,I10,A7)') 'POINTS', U2D_GRID%NUM_NODES, 'double'
               DO I = 1, U2D_GRID%NUM_NODES
                  WRITE(54321,*) U2D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321,'(A,I10,I10)') 'CELLS', U2D_GRID%NUM_CELLS, 4*U2D_GRID%NUM_CELLS 
               DO I = 1, U2D_GRID%NUM_CELLS
                  WRITE(54321,*) 3, (U2D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321,'(A,I10)') 'CELL_TYPES', U2D_GRID%NUM_CELLS
               DO I = 1, U2D_GRID%NUM_CELLS
                  WRITE(54321,*) 5
               END DO
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
               WRITE(54321,'(A)') 'DATASET UNSTRUCTURED_GRID'
               
               WRITE(54321,'(A,I10,A7)') 'POINTS', U3D_GRID%NUM_NODES, 'double'
               DO I = 1, U3D_GRID%NUM_NODES
                  WRITE(54321,*) U3D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321,'(A,I10,I10)') 'CELLS', U3D_GRID%NUM_CELLS, 5*U3D_GRID%NUM_CELLS 
               DO I = 1, U3D_GRID%NUM_CELLS
                  WRITE(54321,*) 4, (U3D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321,'(A,I10)') 'CELL_TYPES', U3D_GRID%NUM_CELLS
               DO I = 1, U3D_GRID%NUM_CELLS
                  WRITE(54321,*) 10
               END DO
            ELSE IF (DIMS == 1) THEN
               WRITE(54321,'(A)') 'DATASET RECTILINEAR_GRID'
               
               WRITE(54321,'(A,I10,I10,I10)') 'DIMENSIONS', NX+1, NY+1, 1 

               WRITE(54321,'(A,I10,A7)') 'X_COORDINATES', NX+1, 'double'
               WRITE(54321,*) XNODES

               WRITE(54321,'(A,I10,A7)') 'Y_COORDINATES', 1, 'double'
               WRITE(54321,*) 0.

               WRITE(54321,'(A,I10,A7)') 'Z_COORDINATES', 1, 'double'
               WRITE(54321,*) 0.
            ELSE IF (DIMS == 2) THEN
               WRITE(54321,'(A)') 'DATASET RECTILINEAR_GRID'
               
               WRITE(54321,'(A,I10,I10,I10)') 'DIMENSIONS', NX+1, NY+1, 1 

               WRITE(54321,'(A,I10,A7)') 'X_COORDINATES', NX+1, 'double'
               WRITE(54321,*) XNODES

               WRITE(54321,'(A,I10,A7)') 'Y_COORDINATES', NY+1, 'double'
               WRITE(54321,*) YNODES

               WRITE(54321,'(A,I10,A7)') 'Z_COORDINATES', 1, 'double'
               WRITE(54321,*) 0.
            END IF
            
            WRITE(54321,'(A,I10)') 'CELL_DATA', NCELLS
            IF (BOOL_DUMP_MOMENTS) THEN
               WRITE(54321,'(A,I10)') 'FIELD FieldData', (13+33)*N_SPECIES+1
            ELSE
               WRITE(54321,'(A,I10)') 'FIELD FieldData', 13*N_SPECIES+1
            END IF


            ! Write per-cell value
            WRITE(54321,'(A,I10,I10,A8)') 'PROC_ID', 1, NCELLS, 'integer'
            WRITE(54321,*) CELL_PROC_ID

            ! Write per-cell, per-species values
            DO JS = 1, N_SPECIES
               FIRST = 1 + (JS-1)*NCELLS
               LAST  = JS*NCELLS
            
               WRITE(string, *) 'number_particles_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_NP(FIRST:LAST)

               WRITE(string, *) 'nrho_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               IF (GRID_TYPE == RECTILINEAR_NONUNIFORM .OR. GRID_TYPE == UNSTRUCTURED .OR. AXI) THEN
                  IF (BOOL_RADIAL_WEIGHTING) THEN
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) THEN
                        WRITE(54321,*) CELL_FNUM*AVG_NP(FIRST:LAST)/CELL_VOLUMES
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
                        WRITE(54321,*) CELL_FNUM*AVG_NP(FIRST:LAST)/U1D_GRID%CELL_VOLUMES
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
                        WRITE(54321,*) CELL_FNUM*AVG_NP(FIRST:LAST)/U2D_GRID%CELL_VOLUMES
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
                        WRITE(54321,*) CELL_FNUM*AVG_NP(FIRST:LAST)/U3D_GRID%CELL_VOLUMES
                     END IF
                  ELSE
                     IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) THEN
                        WRITE(54321,*) FNUM*AVG_NP(FIRST:LAST)/CELL_VOLUMES
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
                        WRITE(54321,*) FNUM*AVG_NP(FIRST:LAST)/U1D_GRID%CELL_VOLUMES
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
                        WRITE(54321,*) FNUM*AVG_NP(FIRST:LAST)/U2D_GRID%CELL_VOLUMES
                     ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
                        WRITE(54321,*) FNUM*AVG_NP(FIRST:LAST)/U3D_GRID%CELL_VOLUMES
                     END IF
                  END IF
               ELSE IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
                  WRITE(54321,*) FNUM*AVG_NP(FIRST:LAST)/CELL_VOL
               END IF

               WRITE(string, *) 'vx_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_VX(FIRST:LAST)

               WRITE(string, *) 'vy_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_VY(FIRST:LAST)

               WRITE(string, *) 'vz_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_VZ(FIRST:LAST)

               WRITE(string, *) 'Ttrx_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_TTRX(FIRST:LAST)

               WRITE(string, *) 'Ttry_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_TTRY(FIRST:LAST)

               WRITE(string, *) 'Ttrz_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_TTRZ(FIRST:LAST)

               WRITE(string, *) 'Ttr_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_TTR(FIRST:LAST)

               WRITE(string, *) 'Trot_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_TROT(FIRST:LAST)

               WRITE(string, *) 'Tvib_mean_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
               WRITE(54321,*) AVG_TVIB(FIRST:LAST)

               WRITE(string, *) 'num_cumulated_one_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A8)') string, 1, NCELLS, 'integer'
               WRITE(54321,*) AVG_CUMULATED_INTENSIVE_ONE(FIRST:LAST)

               WRITE(string, *) 'num_cumulated_two_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A8)') string, 1, NCELLS, 'integer'
               WRITE(54321,*) AVG_CUMULATED_INTENSIVE_TWO(FIRST:LAST)

            
               IF (BOOL_DUMP_MOMENTS) THEN
                  ! Now dump the 33 moments
                  DO MOM = 1, 33
                     WRITE(54321,'(A,A,I10,I10,A7)') TRIM(MOMENT_STRING(MOM)), SPECIES(JS)%NAME, 1, NCELLS, 'double'
                     WRITE(54321,*) AVG_MOMENTS(FIRST:LAST, MOM)
                  END DO
               END IF
         
         
            END DO

            IF (PIC_TYPE .NE. NONE) THEN
               IF (GRID_TYPE == UNSTRUCTURED) THEN
                  WRITE(54321,'(A,I10)') 'FIELD FieldData', 3

                  WRITE(54321,'(A,I10,I10,A7)') 'E_X', 1, NCELLS, 'double'
                  WRITE(54321,*) E_FIELD(1,:,:)

                  WRITE(54321,'(A,I10,I10,A7)') 'E_Y', 1, NCELLS, 'double'
                  WRITE(54321,*) E_FIELD(2,:,:)

                  WRITE(54321,'(A,I10,I10,A7)') 'E_Z', 1, NCELLS, 'double'
                  WRITE(54321,*) E_FIELD(3,:,:)

                  WRITE(54321,'(A,I10)') 'POINT_DATA', NNODES
                  WRITE(54321,'(A,I10)') 'FIELD FieldData', 6
               
                  WRITE(54321,'(A,I10,I10,A7)') 'QRHO', 1, NNODES, 'double'
                  WRITE(54321,*) RHS

                  WRITE(54321,'(A,I10,I10,A7)') 'PHI', 1, NNODES, 'double'
                  WRITE(54321,*) AVG_PHI

                  WRITE(54321,'(A,I10,I10,A7)') 'PHIBAR', 1, NNODES, 'double'
                  WRITE(54321,*) PHIBAR_FIELD

                  WRITE(54321,'(A,I10,I10,A7)') 'B_X', 1, NNODES, 'double'
                  WRITE(54321,*) B_FIELD(1,:,:)

                  WRITE(54321,'(A,I10,I10,A7)') 'B_Y', 1, NNODES, 'double'
                  WRITE(54321,*) B_FIELD(2,:,:)

                  WRITE(54321,'(A,I10,I10,A7)') 'B_Z', 1, NNODES, 'double'
                  WRITE(54321,*) B_FIELD(3,:,:)


               ELSE
                  WRITE(54321,'(A,I10)') 'POINT_DATA', NNODES
                  WRITE(54321,'(A,I10)') 'FIELD FieldData', 5
               
                  WRITE(54321,'(A,I10,I10,A7)') 'QRHO', 1, NNODES, 'double'
                  WRITE(54321,*) RHS

                  WRITE(54321,'(A,I10,I10,A7)') 'PHI', 1, NNODES, 'double'
                  WRITE(54321,*) AVG_PHI

                  WRITE(54321,'(A,I10,I10,A7)') 'E_X', 1, NNODES, 'double'
                  WRITE(54321,*) E_FIELD(1,:,:)

                  WRITE(54321,'(A,I10,I10,A7)') 'E_Y', 1, NNODES, 'double'
                  WRITE(54321,*) E_FIELD(2,:,:)

                  WRITE(54321,'(A,I10,I10,A7)') 'E_Z', 1, NNODES, 'double'
                  WRITE(54321,*) E_FIELD(3,:,:)
               END IF
            END IF

            CLOSE(54321)
         
         END IF

      END IF 

   END SUBROUTINE GRID_SAVE


   SUBROUTINE INIT_POSTPROCESS

      IMPLICIT NONE

      INTEGER :: LENGTH

      LENGTH = NCELLS * N_SPECIES
      
      ALLOCATE(AVG_NP(LENGTH))

      ALLOCATE(AVG_VX(LENGTH))
      ALLOCATE(AVG_VY(LENGTH))
      ALLOCATE(AVG_VZ(LENGTH))
      
      ALLOCATE(AVG_TTRX(LENGTH))
      ALLOCATE(AVG_TTRY(LENGTH))
      ALLOCATE(AVG_TTRZ(LENGTH))
      ALLOCATE(AVG_TTR(LENGTH))

      ALLOCATE(AVG_TROT(LENGTH))
      ALLOCATE(AVG_TVIB(LENGTH))
      
      ALLOCATE(AVG_CUMULATED_INTENSIVE_ONE(LENGTH))
      ALLOCATE(AVG_CUMULATED_INTENSIVE_TWO(LENGTH))

      AVG_NP = 0

      AVG_VX = 0
      AVG_VY = 0
      AVG_VZ = 0

      AVG_TTRX = 0
      AVG_TTRY = 0
      AVG_TTRZ = 0
      AVG_TTR = 0

      AVG_TROT = 0
      AVG_TVIB = 0


      AVG_CUMULATED = 0
      AVG_CUMULATED_INTENSIVE_ONE = 0
      AVG_CUMULATED_INTENSIVE_TWO = 0

      IF (BOOL_DUMP_MOMENTS) THEN
         ALLOCATE(AVG_MOMENTS(LENGTH,33))
         AVG_MOMENTS = 0
      END IF

      IF (PIC_TYPE .NE. NONE) THEN
         ALLOCATE(AVG_PHI(NNODES))
         AVG_PHI = 0
      END IF

      

   END SUBROUTINE INIT_POSTPROCESS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GRID_RESET -> Resets cumulated average !!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GRID_RESET

      IMPLICIT NONE

      AVG_NP = 0

      AVG_VX = 0
      AVG_VY = 0
      AVG_VZ = 0

      AVG_TTRX = 0
      AVG_TTRY = 0
      AVG_TTRZ = 0
      AVG_TTR = 0

      AVG_TROT = 0
      AVG_TVIB = 0

      AVG_CUMULATED = 0
      AVG_CUMULATED_INTENSIVE_ONE = 0
      AVG_CUMULATED_INTENSIVE_TWO = 0

      IF (BOOL_DUMP_MOMENTS) AVG_MOMENTS = 0

   END SUBROUTINE GRID_RESET



   SUBROUTINE TALLY_PARTICLE_TO_BOUNDARY(REFLECTED, PART, IC, IFACE)

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: REFLECTED
      TYPE(PARTICLE_DATA_STRUCTURE), INTENT(IN) :: PART
      INTEGER, INTENT(IN) :: IC, IFACE
      REAL(KIND=8) :: MOLMASS, K, AREA
      INTEGER :: INDEX

      MOLMASS = SPECIES(PART%S_ID)%MOLECULAR_MASS

      IF (DIMS == 1) THEN
         INDEX = U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX(IFACE, IC) + (PART%S_ID-1)*NBOUNDCELLS
         AREA = U0D_GRID%VERTEX_AREAS(U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX(IFACE, IC))
      ELSE IF (DIMS == 2) THEN
         INDEX = U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(IFACE, IC) + (PART%S_ID-1)*NBOUNDCELLS
         AREA = U1D_GRID%SEGMENT_AREAS(U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(IFACE, IC))
      ELSE IF (DIMS == 3) THEN
         INDEX = U3D_GRID%CELL_FACES_BOUNDARY_INDEX(IFACE, IC) + (PART%S_ID-1)*NBOUNDCELLS
         AREA = U2D_GRID%CELL_AREAS(U3D_GRID%CELL_FACES_BOUNDARY_INDEX(IFACE, IC))
      END IF

      K = FNUM*MOLMASS/AREA/DT

      IF (REFLECTED) THEN
         TIMESTEP_NOUT(INDEX) = TIMESTEP_NOUT(INDEX) - FNUM/AREA/DT
         
         TIMESTEP_PXOUT(INDEX) = TIMESTEP_PXOUT(INDEX) - K * PART%VX
         TIMESTEP_PYOUT(INDEX) = TIMESTEP_PYOUT(INDEX) - K * PART%VY
         TIMESTEP_PZOUT(INDEX) = TIMESTEP_PZOUT(INDEX) - K * PART%VZ
                  
         TIMESTEP_EOUT(INDEX) = TIMESTEP_EOUT(INDEX) - K * &
                                (PART%VX * PART%VX + &
                                 PART%VY * PART%VY + &
                                 PART%VZ * PART%VZ )
      ELSE
         TIMESTEP_NIN(INDEX) = TIMESTEP_NIN(INDEX) + FNUM/AREA/DT
         
         TIMESTEP_PXIN(INDEX) = TIMESTEP_PXIN(INDEX) + K * PART%VX
         TIMESTEP_PYIN(INDEX) = TIMESTEP_PYIN(INDEX) + K * PART%VY
         TIMESTEP_PZIN(INDEX) = TIMESTEP_PZIN(INDEX) + K * PART%VZ
                  
         TIMESTEP_EIN(INDEX) = TIMESTEP_EIN(INDEX) + K * &
                               (PART%VX * PART%VX + &
                                PART%VY * PART%VY + &
                                PART%VZ * PART%VZ )
      END IF

   END SUBROUTINE TALLY_PARTICLE_TO_BOUNDARY


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE BOUNDARY_GATHER -> Adds timestep to cumulative average of boundary !!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE BOUNDARY_GATHER

      IMPLICIT NONE

      INTEGER      :: I, IP, IC, INDEX

      REAL(KIND=8) :: DBLE_BOUNDARY_AVG_CUMULATED, E_MAG2


      IF (PIC_TYPE .NE. NONE) THEN
         DO I = 1, NNODES
            IF (DIMS == 1) THEN
               IF (U1D_GRID%NODES_BOUNDARY_INDEX(I) .NE. -1) THEN
                  TIMESTEP_PHI_BOUND(U1D_GRID%NODES_BOUNDARY_INDEX(I)) = PHI_FIELD(I)
                  TIMESTEP_QRHO_BOUND(U1D_GRID%NODES_BOUNDARY_INDEX(I)) = SURFACE_CHARGE(I)
               END IF
            ELSE IF (DIMS == 2) THEN
               IF (U2D_GRID%NODES_BOUNDARY_INDEX(I) .NE. -1) THEN
                  TIMESTEP_PHI_BOUND(U2D_GRID%NODES_BOUNDARY_INDEX(I)) = PHI_FIELD(I)
                  TIMESTEP_QRHO_BOUND(U2D_GRID%NODES_BOUNDARY_INDEX(I)) = SURFACE_CHARGE(I)
               END IF
            ELSE IF (DIMS == 3) THEN
               IF (U3D_GRID%NODES_BOUNDARY_INDEX(I) .NE. -1) THEN
                  TIMESTEP_PHI_BOUND(U3D_GRID%NODES_BOUNDARY_INDEX(I)) = PHI_FIELD(I)
                  TIMESTEP_QRHO_BOUND(U3D_GRID%NODES_BOUNDARY_INDEX(I)) = SURFACE_CHARGE(I)
               END IF
            END IF
         END DO

         DO IC = 1, NCELLS
            E_MAG2 = MAG(E_FIELD(:,1,IC))**2
            IF (DIMS == 1) THEN 
               DO IP = 1, 2
                  IF (U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX(IP, IC) .NE. -1) THEN
                     INDEX = U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX(IP, IC)

                     TIMESTEP_PXEM(INDEX) = TIMESTEP_PXEM(INDEX) - EPS0*&
                     (E_FIELD(1,1,IC)*E_FIELD(1,1,IC) - 0.5*E_MAG2)*U2D_GRID%EDGE_NORMAL(1,IP,IC)
                  END IF
               END DO
            ELSE IF (DIMS == 2) THEN
               DO IP = 1, 3
                  IF (U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(IP, IC) .NE. -1) THEN
                     INDEX = U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(IP, IC)

                     TIMESTEP_PXEM(INDEX) = TIMESTEP_PXEM(INDEX) - EPS0*&
                     ( (E_FIELD(1,1,IC)*E_FIELD(1,1,IC) - 0.5*E_MAG2)*U2D_GRID%EDGE_NORMAL(1,IP,IC) &
                     +  E_FIELD(1,1,IC)*E_FIELD(2,1,IC)*U2D_GRID%EDGE_NORMAL(2,IP,IC))
                     TIMESTEP_PYEM(INDEX) = TIMESTEP_PYEM(INDEX) - EPS0*&
                     (  E_FIELD(2,1,IC)*E_FIELD(1,1,IC)*U2D_GRID%EDGE_NORMAL(1,IP,IC) &
                     + (E_FIELD(2,1,IC)*E_FIELD(2,1,IC) - 0.5*E_MAG2)*U2D_GRID%EDGE_NORMAL(2,IP,IC))
                  END IF
               END DO
            ELSE IF (DIMS == 3) THEN
               DO IP = 1, 4
                  IF (U3D_GRID%CELL_FACES_BOUNDARY_INDEX(IP, IC) .NE. -1) THEN
                     INDEX = U3D_GRID%CELL_FACES_BOUNDARY_INDEX(IP, IC)

                     TIMESTEP_PXEM(INDEX) = TIMESTEP_PXEM(INDEX) - EPS0*&
                     ( (E_FIELD(1,1,IC)*E_FIELD(1,1,IC) - 0.5*E_MAG2)*U3D_GRID%FACE_NORMAL(1,IP,IC) &
                     +  E_FIELD(1,1,IC)*E_FIELD(2,1,IC)*U3D_GRID%FACE_NORMAL(2,IP,IC) &
                     +  E_FIELD(1,1,IC)*E_FIELD(3,1,IC)*U3D_GRID%FACE_NORMAL(3,IP,IC))
                     TIMESTEP_PYEM(INDEX) = TIMESTEP_PYEM(INDEX) - EPS0*&
                     (  E_FIELD(2,1,IC)*E_FIELD(1,1,IC)*U3D_GRID%FACE_NORMAL(1,IP,IC) &
                     + (E_FIELD(2,1,IC)*E_FIELD(2,1,IC) - 0.5*E_MAG2)*U3D_GRID%FACE_NORMAL(2,IP,IC) &
                     +  E_FIELD(2,1,IC)*E_FIELD(3,1,IC)*U3D_GRID%FACE_NORMAL(3,IP,IC))
                     TIMESTEP_PZEM(INDEX) = TIMESTEP_PZEM(INDEX) - EPS0*&
                     (  E_FIELD(3,1,IC)*E_FIELD(1,1,IC)*U3D_GRID%FACE_NORMAL(1,IP,IC) &
                     +  E_FIELD(3,1,IC)*E_FIELD(2,1,IC)*U3D_GRID%FACE_NORMAL(2,IP,IC) &
                     + (E_FIELD(3,1,IC)*E_FIELD(3,1,IC) - 0.5*E_MAG2)*U3D_GRID%FACE_NORMAL(3,IP,IC))
                  END IF
               END DO
            END IF
         END DO
      END IF

      ! Collect data from all the processes (PHI and PEM does not need reduction)
      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_NIN,   NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_NOUT,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_PXIN,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_PYIN,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_PZIN,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_PXOUT, NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_PYOUT, NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_PZOUT, NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_EIN,   NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_EOUT,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_QRHO_BOUND, NBOUNDNODES,     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      ELSE
         CALL MPI_REDUCE(TIMESTEP_NIN,  TIMESTEP_NIN,   NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_NOUT, TIMESTEP_NOUT,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_PXIN, TIMESTEP_PXIN,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_PYIN, TIMESTEP_PYIN,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_PZIN, TIMESTEP_PZIN,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_PXOUT,TIMESTEP_PXOUT, NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_PYOUT,TIMESTEP_PYOUT, NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_PZOUT,TIMESTEP_PZOUT, NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_EIN,  TIMESTEP_EIN,   NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_EOUT, TIMESTEP_EOUT,  NBOUNDCELLS*N_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(TIMESTEP_QRHO_BOUND,TIMESTEP_QRHO_BOUND, NBOUNDNODES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      END IF

      ! Add to cumulated average
      DBLE_BOUNDARY_AVG_CUMULATED = DBLE(BOUNDARY_AVG_CUMULATED)

      AVG_NIN =   (AVG_NIN*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_NIN)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_NOUT =  (AVG_NOUT*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_NOUT)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_PXIN =  (AVG_PXIN*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PXIN)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_PYIN =  (AVG_PYIN*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PYIN)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_PZIN =  (AVG_PZIN*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PZIN)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_PXOUT = (AVG_PXOUT*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PXOUT)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_PYOUT = (AVG_PYOUT*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PYOUT)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_PZOUT = (AVG_PZOUT*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PZOUT)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_EIN =   (AVG_EIN*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_EIN)/(BOUNDARY_AVG_CUMULATED + 1.)
      AVG_EOUT =  (AVG_EOUT*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_EOUT)/(BOUNDARY_AVG_CUMULATED + 1.)


      IF (PIC_TYPE .NE. NONE) THEN
         AVG_PHI_BOUND = (AVG_PHI_BOUND*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PHI_BOUND)/(BOUNDARY_AVG_CUMULATED + 1.)
         AVG_QRHO_BOUND = (AVG_QRHO_BOUND*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_QRHO_BOUND)/(BOUNDARY_AVG_CUMULATED + 1.)
         AVG_PXEM =  (AVG_PXEM*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PXEM)/(BOUNDARY_AVG_CUMULATED + 1.)
         AVG_PYEM =  (AVG_PYEM*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PYEM)/(BOUNDARY_AVG_CUMULATED + 1.)
         AVG_PZEM =  (AVG_PZEM*DBLE_BOUNDARY_AVG_CUMULATED + TIMESTEP_PZEM)/(BOUNDARY_AVG_CUMULATED + 1.)
      END IF
      
      BOUNDARY_AVG_CUMULATED = BOUNDARY_AVG_CUMULATED + 1


      TIMESTEP_NIN = 0
      TIMESTEP_NOUT = 0

      TIMESTEP_PXIN = 0
      TIMESTEP_PYIN = 0
      TIMESTEP_PZIN = 0
      TIMESTEP_PXOUT = 0
      TIMESTEP_PYOUT = 0
      TIMESTEP_PZOUT = 0

      TIMESTEP_PXEM = 0
      TIMESTEP_PYEM = 0
      TIMESTEP_PZEM = 0

      TIMESTEP_EIN = 0
      TIMESTEP_EOUT = 0

      TIMESTEP_PHI_BOUND = 0
      TIMESTEP_QRHO_BOUND = 0
      
   END SUBROUTINE BOUNDARY_GATHER



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GRID_SAVE -> Saves cumulated average !!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE BOUNDARY_SAVE
         
      IMPLICIT NONE

      CHARACTER*256                      :: string, file_name

      INTEGER                            :: I, JS, FIRST, LAST


      IF (PROC_ID .EQ. 0) THEN
         
         ! Write to file.

         !WRITE(file_name,'(A, I0, A)') '/media/pietro/Storage/panteradumps/dsmc_flowfield_', tID, '.vtk'
         WRITE(file_name,'(A, A, I0, A)') TRIM(ADJUSTL(BOUNDARY_SAVE_PATH)), 'dsmc_boundary_', tID, '.vtk'

         IF (BOOL_BINARY_OUTPUT) THEN
            OPEN(54321, FILE=file_name, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='NEW', CONVERT='BIG_ENDIAN')

            WRITE(54321) '# vtk DataFile Version 3.0'//ACHAR(10)
            WRITE(54321) 'vtk output'//ACHAR(10)
            WRITE(54321) 'BINARY'//ACHAR(10)

            IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
               WRITE(54321) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

               WRITE(54321) 'POINTS '//ITOA(U0D_GRID%NUM_NODES)//' double'//ACHAR(10)
               DO I = 1, U0D_GRID%NUM_NODES
                  WRITE(54321) U0D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321) 'CELLS '//ITOA(U0D_GRID%NUM_POINTS)//' '//ITOA(2*U0D_GRID%NUM_POINTS)//ACHAR(10)
               DO I = 1, U0D_GRID%NUM_POINTS
                  WRITE(54321) 1, (U0D_GRID%POINT_NODES(I) - 1)
               END DO

               WRITE(54321) 'CELL_TYPES '//ITOA(U0D_GRID%NUM_POINTS)//ACHAR(10)
               DO I = 1, U0D_GRID%NUM_POINTS
                  WRITE(54321) 1
               END DO
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
               WRITE(54321) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

               WRITE(54321) 'POINTS '//ITOA(U1D_GRID%NUM_NODES)//' double'//ACHAR(10)
               DO I = 1, U1D_GRID%NUM_NODES
                  WRITE(54321) U1D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321) 'CELLS '//ITOA(U1D_GRID%NUM_CELLS)//' '//ITOA(3*U1D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U1D_GRID%NUM_CELLS
                  WRITE(54321) 2, (U1D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321) 'CELL_TYPES '//ITOA(U1D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U1D_GRID%NUM_CELLS
                  WRITE(54321) 3
               END DO
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
               WRITE(54321) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

               WRITE(54321) 'POINTS '//ITOA(U2D_GRID%NUM_NODES)//' double'//ACHAR(10)
               DO I = 1, U2D_GRID%NUM_NODES
                  WRITE(54321) U2D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321) 'CELLS '//ITOA(U2D_GRID%NUM_CELLS)//' '//ITOA(4*U2D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U2D_GRID%NUM_CELLS
                  WRITE(54321) 3, (U2D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321) 'CELL_TYPES '//ITOA(U2D_GRID%NUM_CELLS)//ACHAR(10)
               DO I = 1, U2D_GRID%NUM_CELLS
                  WRITE(54321) 5
               END DO
            ELSE
               CALL ERROR_ABORT('Boundary dump not implemented on structured grids!')
            END IF

            WRITE(54321) 'CELL_DATA '//ITOA(NBOUNDCELLS)//ACHAR(10)

            WRITE(54321) 'FIELD FieldData '//ITOA( 10*N_SPECIES + 1 )//ACHAR(10)

            ! Write per-cell value
            WRITE(54321) 'PHYSICAL_GROUP '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' integer'//ACHAR(10)
            IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
               WRITE(54321) U0D_GRID%POINT_PG, ACHAR(10)
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
               WRITE(54321) U1D_GRID%CELL_PG, ACHAR(10)
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
               WRITE(54321) U2D_GRID%CELL_PG, ACHAR(10)
            END IF

            ! Write per-cell, per-species values
            DO JS = 1, N_SPECIES
               FIRST = 1 + (JS-1)*NBOUNDCELLS
               LAST  = JS*NBOUNDCELLS
            
               WRITE(string, *) 'flux_in_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_NIN(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'flux_out_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_NOUT(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'mom_x_in_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PXIN(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'mom_y_in_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PYIN(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'mom_z_in_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PZIN(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'mom_x_out_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PXOUT(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'mom_y_out_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PYOUT(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'mom_z_out_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PZOUT(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'ene_in_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_EIN(FIRST:LAST), ACHAR(10)

               WRITE(string, *) 'ene_out_', SPECIES(JS)%NAME
               WRITE(54321) string//' '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_EOUT(FIRST:LAST), ACHAR(10)

            END DO

            IF (PIC_TYPE .NE. NONE) THEN

               WRITE(54321) 'FIELD FieldData '//ITOA(3)//ACHAR(10)

               WRITE(54321) 'pem_x '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PXEM, ACHAR(10)

               WRITE(54321) 'pem_y '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PYEM, ACHAR(10)

               WRITE(54321) 'pem_z '//ITOA(1)//' '//ITOA(NBOUNDCELLS)//' double'//ACHAR(10)
               WRITE(54321) AVG_PZEM, ACHAR(10)


               WRITE(54321) 'POINT_DATA '//ITOA(NBOUNDNODES)//ACHAR(10)

               WRITE(54321) 'FIELD FieldData '//ITOA(2)//ACHAR(10)
            
               WRITE(54321) 'QRHO '//ITOA(1)//' '//ITOA(NBOUNDNODES)//' double'//ACHAR(10)
               WRITE(54321) AVG_QRHO_BOUND, ACHAR(10)

               WRITE(54321) 'PHI '//ITOA(1)//' '//ITOA(NBOUNDNODES)//' double'//ACHAR(10)
               WRITE(54321) AVG_PHI_BOUND, ACHAR(10)

            END IF

            CLOSE(54321)
            
         ELSE  ! Write ASCII output.
            OPEN(54321, FILE=file_name, ACCESS='SEQUENTIAL', FORM='FORMATTED', STATUS='NEW')

            WRITE(54321,'(A)') '# vtk DataFile Version 3.0'
            WRITE(54321,'(A)') 'vtk output'
            WRITE(54321,'(A)') 'ASCII'

            IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
               WRITE(54321,'(A)') 'DATASET UNSTRUCTURED_GRID'
               
               WRITE(54321,'(A,I10,A7)') 'POINTS', U0D_GRID%NUM_NODES, 'double'
               DO I = 1, U0D_GRID%NUM_NODES
                  WRITE(54321,*) U0D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321,'(A,I10,I10)') 'CELLS', U0D_GRID%NUM_POINTS, 2*U0D_GRID%NUM_POINTS 
               DO I = 1, U0D_GRID%NUM_POINTS
                  WRITE(54321,*) 1, (U0D_GRID%POINT_NODES(I) - 1)
               END DO

               WRITE(54321,'(A,I10)') 'CELL_TYPES', U0D_GRID%NUM_POINTS
               DO I = 1, U0D_GRID%NUM_POINTS
                  WRITE(54321,*) 1
               END DO
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
               WRITE(54321,'(A)') 'DATASET UNSTRUCTURED_GRID'
               
               WRITE(54321,'(A,I10,A7)') 'POINTS', U1D_GRID%NUM_NODES, 'double'
               DO I = 1, U1D_GRID%NUM_NODES
                  WRITE(54321,*) U1D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321,'(A,I10,I10)') 'CELLS', U1D_GRID%NUM_CELLS, 3*U1D_GRID%NUM_CELLS 
               DO I = 1, U1D_GRID%NUM_CELLS
                  WRITE(54321,*) 2, (U1D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321,'(A,I10)') 'CELL_TYPES', U1D_GRID%NUM_CELLS
               DO I = 1, U1D_GRID%NUM_CELLS
                  WRITE(54321,*) 3
               END DO
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
               WRITE(54321,'(A)') 'DATASET UNSTRUCTURED_GRID'
               
               WRITE(54321,'(A,I10,A7)') 'POINTS', U2D_GRID%NUM_NODES, 'double'
               DO I = 1, U2D_GRID%NUM_NODES
                  WRITE(54321,*) U2D_GRID%NODE_COORDS(:,I)
               END DO

               WRITE(54321,'(A,I10,I10)') 'CELLS', U2D_GRID%NUM_CELLS, 4*U2D_GRID%NUM_CELLS 
               DO I = 1, U2D_GRID%NUM_CELLS
                  WRITE(54321,*) 3, (U2D_GRID%CELL_NODES(:,I) - 1)
               END DO

               WRITE(54321,'(A,I10)') 'CELL_TYPES', U2D_GRID%NUM_CELLS
               DO I = 1, U2D_GRID%NUM_CELLS
                  WRITE(54321,*) 5
               END DO
            ELSE
               CALL ERROR_ABORT('Boundary dump not implemented on structured grids!')
            END IF
            
            WRITE(54321,'(A,I10)') 'CELL_DATA', NBOUNDCELLS

            WRITE(54321,'(A,I10)') 'FIELD FieldData', 10*N_SPECIES + 1

            ! Write per-cell value
            WRITE(54321,'(A,I10,I10,A8)') 'PHYSICAL_GROUP', 1, NBOUNDCELLS, 'integer'
            IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 1) THEN
               WRITE(54321,*) U0D_GRID%POINT_PG
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN
               WRITE(54321,*) U1D_GRID%CELL_PG
            ELSE IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN
               WRITE(54321,*) U2D_GRID%CELL_PG
            END IF

            ! Write per-cell, per-species values
            DO JS = 1, N_SPECIES
               FIRST = 1 + (JS-1)*NBOUNDCELLS
               LAST  = JS*NBOUNDCELLS
            
               WRITE(string, *) 'flux_in_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_NIN(FIRST:LAST)

               WRITE(string, *) 'flux_out_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_NOUT(FIRST:LAST)

               WRITE(string, *) 'mom_x_in_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_PXIN(FIRST:LAST)

               WRITE(string, *) 'mom_y_in_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_PYIN(FIRST:LAST)

               WRITE(string, *) 'mom_z_in_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_PZIN(FIRST:LAST)

               WRITE(string, *) 'mom_x_out_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_PXOUT(FIRST:LAST)

               WRITE(string, *) 'mom_y_out_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_PYOUT(FIRST:LAST)

               WRITE(string, *) 'mom_z_out_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_PZOUT(FIRST:LAST)

               WRITE(string, *) 'ene_in_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_EIN(FIRST:LAST)

               WRITE(string, *) 'ene_out_', SPECIES(JS)%NAME
               WRITE(54321,'(A,I10,I10,A7)') string, 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) AVG_EOUT(FIRST:LAST)

            END DO

            IF (PIC_TYPE .NE. NONE) THEN

               WRITE(54321,'(A,I10)') 'FIELD FieldData', 3

               WRITE(54321,'(A,I10,I10,A7)') 'pem_x', 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) E_FIELD(1,:,:)

               WRITE(54321,'(A,I10,I10,A7)') 'pem_y', 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) E_FIELD(2,:,:)

               WRITE(54321,'(A,I10,I10,A7)') 'pem_z', 1, NBOUNDCELLS, 'double'
               WRITE(54321,*) E_FIELD(3,:,:)

               WRITE(54321,'(A,I10)') 'POINT_DATA', NBOUNDNODES
               WRITE(54321,'(A,I10)') 'FIELD FieldData', 2
            
               WRITE(54321,'(A,I10,I10,A7)') 'QRHO', 1, NBOUNDNODES, 'double'
               WRITE(54321,*) AVG_QRHO_BOUND

               WRITE(54321,'(A,I10,I10,A7)') 'PHI', 1, NBOUNDNODES, 'double'
               WRITE(54321,*) AVG_PHI_BOUND

            END IF

            CLOSE(54321)
         
         END IF

      END IF 

   END SUBROUTINE BOUNDARY_SAVE


   SUBROUTINE INIT_BOUNDARY_POSTPROCESS

      IMPLICIT NONE

      INTEGER :: LENGTH

      LENGTH = NBOUNDCELLS * N_SPECIES

      ALLOCATE(AVG_NIN(LENGTH))
      ALLOCATE(AVG_NOUT(LENGTH))

      ALLOCATE(AVG_PXIN(LENGTH))
      ALLOCATE(AVG_PYIN(LENGTH))
      ALLOCATE(AVG_PZIN(LENGTH))
      ALLOCATE(AVG_PXOUT(LENGTH))
      ALLOCATE(AVG_PYOUT(LENGTH))
      ALLOCATE(AVG_PZOUT(LENGTH))

      ALLOCATE(AVG_PXEM(NBOUNDCELLS))
      ALLOCATE(AVG_PYEM(NBOUNDCELLS))
      ALLOCATE(AVG_PZEM(NBOUNDCELLS))

      ALLOCATE(AVG_EIN(LENGTH))
      ALLOCATE(AVG_EOUT(LENGTH))

      ALLOCATE(AVG_PHI_BOUND(NBOUNDNODES))
      ALLOCATE(AVG_QRHO_BOUND(NBOUNDNODES))

      AVG_NIN = 0
      AVG_NOUT = 0

      AVG_PXIN = 0
      AVG_PYIN = 0
      AVG_PZIN = 0
      AVG_PXOUT = 0
      AVG_PYOUT = 0
      AVG_PZOUT = 0

      AVG_PXEM = 0
      AVG_PYEM = 0
      AVG_PZEM = 0

      AVG_EIN = 0
      AVG_EOUT = 0

      AVG_PHI_BOUND = 0
      AVG_QRHO_BOUND = 0

      ALLOCATE(TIMESTEP_NIN(LENGTH))
      ALLOCATE(TIMESTEP_NOUT(LENGTH))

      ALLOCATE(TIMESTEP_PXIN(LENGTH))
      ALLOCATE(TIMESTEP_PYIN(LENGTH))
      ALLOCATE(TIMESTEP_PZIN(LENGTH))
      ALLOCATE(TIMESTEP_PXOUT(LENGTH))
      ALLOCATE(TIMESTEP_PYOUT(LENGTH))
      ALLOCATE(TIMESTEP_PZOUT(LENGTH))

      ALLOCATE(TIMESTEP_PXEM(NBOUNDCELLS))
      ALLOCATE(TIMESTEP_PYEM(NBOUNDCELLS))
      ALLOCATE(TIMESTEP_PZEM(NBOUNDCELLS))

      ALLOCATE(TIMESTEP_EIN(LENGTH))
      ALLOCATE(TIMESTEP_EOUT(LENGTH))

      ALLOCATE(TIMESTEP_PHI_BOUND(NBOUNDNODES))
      ALLOCATE(TIMESTEP_QRHO_BOUND(NBOUNDNODES))

      TIMESTEP_NIN = 0
      TIMESTEP_NOUT = 0

      TIMESTEP_PXIN = 0
      TIMESTEP_PYIN = 0
      TIMESTEP_PZIN = 0
      TIMESTEP_PXOUT = 0
      TIMESTEP_PYOUT = 0
      TIMESTEP_PZOUT = 0

      TIMESTEP_PXEM = 0
      TIMESTEP_PYEM = 0
      TIMESTEP_PZEM = 0

      TIMESTEP_EIN = 0
      TIMESTEP_EOUT = 0

      TIMESTEP_PHI_BOUND = 0
      TIMESTEP_QRHO_BOUND = 0

      BOUNDARY_AVG_CUMULATED = 0

   END SUBROUTINE INIT_BOUNDARY_POSTPROCESS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GRID_RESET -> Resets cumulated average !!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE BOUNDARY_RESET

      IMPLICIT NONE

      AVG_NIN = 0
      AVG_NOUT = 0

      AVG_PXIN = 0
      AVG_PYIN = 0
      AVG_PZIN = 0
      AVG_PXOUT = 0
      AVG_PYOUT = 0
      AVG_PZOUT = 0

      AVG_PXEM = 0
      AVG_PYEM = 0
      AVG_PZEM = 0

      AVG_EIN = 0
      AVG_EOUT = 0

      AVG_PHI_BOUND = 0
      AVG_QRHO_BOUND = 0

      BOUNDARY_AVG_CUMULATED = 0

   END SUBROUTINE BOUNDARY_RESET


   SUBROUTINE CHECKS

      IMPLICIT NONE

      INTEGER                            :: JP, JS, JC, JR
   
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TOT_NUM, TOT_REACT_COUNTS
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: TOT_EE_PART, TOT_KE_PART

      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: TOT_MOMENTUM
      REAL(KIND=8)                       :: TOT_KE, TOT_IE, TOT_FE, TOT_EE_FIELD, PHI, CURRENT_TIME, FIELD_POWER_TOT
      REAL(KIND=8)                       :: CFNUM, VOL

      CHARACTER*256                      :: file_name
      CHARACTER*2048                     :: HEADER_STRING
      LOGICAL                            :: FILE_EXISTS

      ALLOCATE(TOT_NUM(N_SPECIES))
      ALLOCATE(TOT_EE_PART(N_SPECIES))
      ALLOCATE(TOT_KE_PART(N_SPECIES))
      ALLOCATE(TOT_MOMENTUM(3,N_SPECIES))
      ALLOCATE(TOT_REACT_COUNTS(N_REACTIONS))

      TOT_NUM = 0
      TOT_MOMENTUM = 0
      TOT_KE = 0
      TOT_IE = 0
      TOT_FE = 0
      TOT_EE_FIELD = 0
      TOT_EE_PART = 0
      TOT_KE_PART = 0
      TOT_REACT_COUNTS = 0

      CURRENT_TIME = tID*DT

      ! Number of particles
      DO JP = 1, NP_PROC

         JS = particles(JP)%S_ID

         TOT_NUM(JS) = TOT_NUM(JS) + 1

         IF (BOOL_RADIAL_WEIGHTING) THEN
            CFNUM = CELL_FNUM(particles(JP)%IC)
         ELSE
            CFNUM = FNUM
         END IF

         ! Momentum
         TOT_MOMENTUM(1,JS) = TOT_MOMENTUM(1,JS) + SPECIES(JS)%MOLECULAR_MASS*particles(JP)%VX * CFNUM
         TOT_MOMENTUM(2,JS) = TOT_MOMENTUM(2,JS) + SPECIES(JS)%MOLECULAR_MASS*particles(JP)%VY * CFNUM
         TOT_MOMENTUM(3,JS) = TOT_MOMENTUM(3,JS) + SPECIES(JS)%MOLECULAR_MASS*particles(JP)%VZ * CFNUM
         

         ! Kinietic energy
         TOT_KE_PART(JS) = TOT_KE_PART(JS) &
         + 0.5*SPECIES(JS)%MOLECULAR_MASS*(particles(JP)%VX**2+particles(JP)%VY**2+particles(JP)%VZ**2) * CFNUM
         TOT_IE = TOT_IE + (particles(JP)%EROT + particles(JP)%EVIB) * CFNUM

         !IF (JS == 4) THEN
         !  TOT_FE = TOT_FE + 15.63e-19/2.
         !END IF
         IF (PIC_TYPE .NE. NONE) THEN
            CALL APPLY_POTENTIAL(JP, PHI)
            TOT_EE_PART(JS)  = TOT_EE_PART(JS) + 0.5*PHI*QE*SPECIES(JS)%CHARGE * CFNUM
         END IF

      END DO

      DO JR = 1, N_REACTIONS
         TOT_REACT_COUNTS(JR) = TOT_REACT_COUNTS(JR) + REACTIONS(JR)%COUNTS
      END DO

      ! Collect data from all the processes and print it
      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_NUM,      N_SPECIES, MPI_INTEGER,  MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_MOMENTUM, 3*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_KE_PART,  N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_IE,       1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_EE_PART,  N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(FIELD_POWER, FIELD_POWER_TOT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE, TOT_REACT_COUNTS, N_REACTIONS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

         IF (GRID_TYPE == UNSTRUCTURED) THEN
            TOT_EE_FIELD = 0.d0
            DO JC = 1, NCELLS
               IF (DIMS == 1) THEN
                  VOL = U1D_GRID%CELL_VOLUMES(JC)
               ELSE IF (DIMS == 2) THEN
                  VOL = U2D_GRID%CELL_VOLUMES(JC)
               ELSE IF (DIMS == 3) THEN
                  VOL = U3D_GRID%CELL_VOLUMES(JC)
               END IF
               TOT_EE_FIELD = TOT_EE_FIELD + (E_FIELD(1, 1, JC)**2 + E_FIELD(2, 1, JC)**2 + E_FIELD(3, 1, JC)**2)*VOL
            END DO
            TOT_EE_FIELD = TOT_EE_FIELD *0.5*EPS0
         END IF
         !CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_FE,       1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

         !HX = (XMAX-XMIN)/DBLE(NX)
         !HY = (YMAX-YMIN)/DBLE(NY)
         !TOT_EE = -HX*HY*8.8541878128E-12*SUM( RHS*PACK(PHI_FIELD, .TRUE.) )/FNUM

         ! WRITE(*,*) ' '
         ! WRITE(*,*) 'Conservation checks:'
         ! WRITE(*,*) 'Number of particles:  ', TOT_NUM
         ! WRITE(*,*) 'Total momentum:       [ ', TOT_MOMENTUM, ' ] [kg m/s]'
         ! WRITE(*,*) 'Total kinetic energy:   ', TOT_KE, ' [J]'
         ! WRITE(*,*) 'Total internal energy:  ', TOT_IE, ' [J]'
         ! WRITE(*,*) 'Total formation energy: ', TOT_FE, ' [J]'
         ! WRITE(*,*) 'Total electrostatic energy: ', TOT_EE, ' [J]'
         ! WRITE(*,*) 'Total energy:          ', TOT_KE+TOT_IE+TOT_FE+TOT_EE, ' [J]'
         ! WRITE(*,*) ' '

         WRITE(file_name,'(A, A)') TRIM(ADJUSTL(CHECKS_SAVE_PATH)), 'conservation_checks'

         INQUIRE(FILE=file_name, EXIST=FILE_EXISTS)

         OPEN(54331, FILE=file_name, POSITION='append', STATUS='unknown', ACTION='write')
         IF (.NOT. FILE_EXISTS) THEN
            HEADER_STRING = ''
            HEADER_STRING = TRIM(HEADER_STRING) // 'time'
            DO JS = 1, N_SPECIES
               HEADER_STRING = TRIM(HEADER_STRING) // ' npart_' // TRIM(SPECIES(JS)%NAME)
            END DO
            DO JS = 1, N_SPECIES
               HEADER_STRING = TRIM(HEADER_STRING) // ' xmom_' // TRIM(SPECIES(JS)%NAME) // ' ' &
                                             // 'ymom_' // TRIM(SPECIES(JS)%NAME) // ' ' &
                                             // 'zmom_' // TRIM(SPECIES(JS)%NAME)
            END DO
            HEADER_STRING = TRIM(HEADER_STRING) // ' totxmom totymom totzmom'
            DO JS = 1, N_SPECIES
               HEADER_STRING = TRIM(HEADER_STRING) // ' ek_' // TRIM(SPECIES(JS)%NAME)
            END DO
            HEADER_STRING = TRIM(HEADER_STRING) // ' totei'
            DO JS = 1, N_SPECIES
               HEADER_STRING = TRIM(HEADER_STRING) // ' eepart_' // TRIM(SPECIES(JS)%NAME)
            END DO
            HEADER_STRING = TRIM(HEADER_STRING) // ' toteefield tote fieldpower coilcurrent'
            DO JS = 1, N_REACTIONS
               HEADER_STRING = TRIM(HEADER_STRING) // ' nreact_' // ITOA(JS)
            END DO

            WRITE(54331,*) TRIM(HEADER_STRING)
         END IF

         WRITE(54331,*) CURRENT_TIME, TOT_NUM, TOT_MOMENTUM, SUM(TOT_MOMENTUM, DIM=2), &
         TOT_KE_PART, TOT_IE, TOT_EE_PART,TOT_EE_FIELD, &
         SUM(TOT_KE_PART) + TOT_IE + TOT_EE_FIELD, &
         FIELD_POWER_TOT, COIL_CURRENT, TOT_REACT_COUNTS !TOT_FE, TOT_EE
         CLOSE(54331)

      ELSE
         CALL MPI_REDUCE(TOT_NUM,       TOT_NUM,      N_SPECIES, MPI_INTEGER,  MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TOT_MOMENTUM,  TOT_MOMENTUM, 3*N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TOT_KE_PART,   TOT_KE_PART,  N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TOT_IE,        TOT_IE,       1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
         CALL MPI_REDUCE(TOT_EE_PART,   TOT_EE_PART,  N_SPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(FIELD_POWER, FIELD_POWER_TOT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(TOT_REACT_COUNTS, TOT_REACT_COUNTS, N_REACTIONS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

         !CALL MPI_REDUCE(TOT_FE,        TOT_FE,       1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF

   END SUBROUTINE CHECKS


   SUBROUTINE DUMP_PARTICLES_VTK(TIMESTEP)

      USE global
      USE mpi_common

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename
      INTEGER :: IP

      WRITE(filename, "(A,A,I0.5,A6,I0.8,A4)") TRIM(ADJUSTL(PARTDUMP_SAVE_PATH)), "proc_", PROC_ID, "_time_", TIMESTEP, ".vtk" ! Compose filename

      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(20, FILE=filename, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='NEW', CONVERT='BIG_ENDIAN')

         WRITE(20) '# vtk DataFile Version 2.0'//ACHAR(10)
         WRITE(20) 'vtk output'//ACHAR(10)
         WRITE(20) 'BINARY'//ACHAR(10)

         WRITE(20) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

         WRITE(20) 'POINTS '//ITOA(NP_PROC)//' double'//ACHAR(10)
         DO IP = 1, NP_PROC
            WRITE(20) particles(IP)%X, particles(IP)%Y, particles(IP)%Z
         END DO

         WRITE(20) 'CELL_TYPES '//ITOA(NP_PROC)//ACHAR(10)
         DO IP = 1, NP_PROC
            WRITE(20) 1
         END DO

         WRITE(20) 'POINT_DATA '//ITOA(NP_PROC)//ACHAR(10)
         WRITE(20) 'VECTORS Velocity double'//ACHAR(10)

         DO IP = 1, NP_PROC
            WRITE(20) particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ
         END DO

      ELSE  ! Write ASCII output.
         OPEN(20, FILE=filename, ACCESS='SEQUENTIAL', FORM='FORMATTED', STATUS='NEW')

         WRITE(20,'(A)') '# vtk DataFile Version 2.0'
         WRITE(20,'(A)') 'vtk output'
         WRITE(20,'(A)') 'ASCII'

         WRITE(20,'(A)') 'DATASET UNSTRUCTURED_GRID'
         
         WRITE(20,'(A,I10,A7)') 'POINTS', NP_PROC, ' double'
         DO IP = 1, NP_PROC
            WRITE(20,*) particles(IP)%X, particles(IP)%Y, particles(IP)%Z
         END DO

         WRITE(20,'(A,I10)') 'CELL_TYPES', NP_PROC
         DO IP = 1, NP_PROC
            WRITE(20,*) 1
         END DO

         WRITE(20,'(A,I10)') 'POINT_DATA', NP_PROC
         WRITE(20,'(A)') 'VECTORS Velocity double'
         DO IP = 1, NP_PROC
            WRITE(20,*) particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ
         END DO

         WRITE(20,'(A,I10)') 'FIELD FieldData', 1

         WRITE(20,'(A,I10,I10,A8)') 'IC', 1, NP_PROC,  ' integer'
         DO IP = 1, NP_PROC
            WRITE(20,*) particles(IP)%IC
         END DO

      END IF

   END SUBROUTINE DUMP_PARTICLES_VTK




END MODULE postprocess
