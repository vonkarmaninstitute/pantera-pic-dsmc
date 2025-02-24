! Copyright (C) 2025 von Karman Institute for Fluid Dynamics (VKI)
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

! This module contains the structure of particles and useful functions
! for handling them

MODULE particle

USE mpi_common

   ! Define particle type and arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   TYPE PARTICLE_DATA_STRUCTURE
      REAL(KIND=8)    :: X, Y, Z        ! position
      REAL(KIND=8)    :: VX, VY, VZ     ! velocities and internal energy
      REAL(KIND=8)    :: EROT, EVIB     ! internal energies
      REAL(KIND=8)    :: DTRIM          ! Remaining time for advection
      INTEGER         :: IC             ! Cell index 
      INTEGER         :: S_ID           ! Species ID
      INTEGER(KIND=8) :: ID             ! Particle identifier
      LOGICAL         :: DUMP_TRAJ      ! The trajectory of this particle should be dumped
   END TYPE PARTICLE_DATA_STRUCTURE
 
   INTEGER(8) :: PARTICLE_ID_COUNTER = 0

   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INIT_PARTICLE -> initializes a particle, given required fields !!!
   ! By using this subroutine, the one does not have to modify everything if a !!!
   ! field is added/removed to/from the particle type.                         !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE INIT_PARTICLE(X, Y, Z, VX, VY, VZ, EROT, EVIB, S_ID, IC, DTRIM, particlept)

      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(IN) :: X, Y, Z, VX, VY, VZ, EROT, EVIB, DTRIM
      INTEGER, INTENT(IN)      :: S_ID, IC

      TYPE(PARTICLE_DATA_STRUCTURE), INTENT(INOUT) :: particlept
     
      particlept%X    = X
      particlept%Y    = Y
      particlept%Z    = Z
 
      particlept%VX   = VX
      particlept%VY   = VY
      particlept%VZ   = VZ

      particlept%EROT   = EROT
      particlept%EVIB   = EVIB
      
 
      particlept%S_ID = S_ID 
      particlept%IC   = IC

      particlept%DTRIM = DTRIM

      particlept%ID = PROC_ID + ISHFT(PARTICLE_ID_COUNTER, 8)
      PARTICLE_ID_COUNTER = PARTICLE_ID_COUNTER + 1

      particlept%DUMP_TRAJ = .FALSE.
      
   END SUBROUTINE INIT_PARTICLE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ADD_PARTICLE_ARRAY -> adds a particle to "particles" array     !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ADD_PARTICLE_ARRAY(particleNOW, NP, particlesARRAY)

      ! This subroutine adds a particle to the particlesARRAY, which contains 
      ! NP particles.
      ! "particlesARRAY" can have a size larger than NP, so that it is not necessary
      ! to change its size every time one more particle is added.
      ! 
      ! For copying a particle, we first check that there is space. If the space is 
      ! over, the vector is enlarged 5 times, and the particle is copied>
      ! 
      ! 1) Increase NP (global variable)
      ! 2) Check if there is enough space in the vector
      !   YES? -> Add particle at new position NP
      !   NO?  -> Make a backup copy of it, make it 5 times larger (reallocate) and then do it.

      IMPLICIT NONE

      TYPE(PARTICLE_DATA_STRUCTURE), INTENT(IN)                                :: particleNOW
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: particlesARRAY
      INTEGER, INTENT(INOUT) :: NP ! Number of particles in particlesARRAY (which is /= of its size!)

      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: COPYARRAY

      IF (NP > 0) THEN
         NP = NP + 1 ! Increase number of particles in array

         IF ( NP > SIZE(particlesARRAY) ) THEN  ! Array too small

         ! Make backup copy of particles array
         ALLOCATE(COPYARRAY(NP-1))
         COPYARRAY(1:NP-1) = particlesARRAY(1:NP-1)

         ! Reallocate particles array to 5 times its size
         DEALLOCATE(particlesARRAY)
         ALLOCATE(particlesARRAY(5*NP))

         ! Re-copy the elements 
         particlesARRAY(1:NP-1) = COPYARRAY(1:NP-1)

         DEALLOCATE(COPYARRAY)
         END IF
         
         ! Add particle to array
         particlesARRAY(NP) = particleNOW
      ELSE
         NP = NP + 1 ! Increase number of particles in array
         IF (.NOT. ALLOCATED(particlesARRAY)) ALLOCATE(particlesARRAY(1))
         particlesARRAY(1) = particleNOW
      END IF

   END SUBROUTINE ADD_PARTICLE_ARRAY

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTOINE REMOVE_PARTICLE_ARRAY -> Removes a particle from the array !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE REMOVE_PARTICLE_ARRAY(ID_REMOVE, particlesARRAY, NP_ARRAY)

      ! This subroutine removes a particle from the particlesARRAY, whose length is NP_ARRAY.
      ! The last particle in the array is copied in place of the particle to be removed "ID_REMOVE".
      ! Then, NP_ARRAY is decremented by 1.

      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: ID_REMOVE
      INTEGER, INTENT(INOUT) :: NP_ARRAY
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), INTENT(INOUT)  :: particlesARRAY
 
      ! First, copy the last particle in the place of the particle to be removed
      !CALL INIT_PARTICLE(particlesARRAY(NP_ARRAY)%X,  particlesARRAY(NP_ARRAY)%Y,    particlesARRAY(NP_ARRAY)%Z,  &
      !                   particlesARRAY(NP_ARRAY)%VX, particlesARRAY(NP_ARRAY)%VY,   particlesARRAY(NP_ARRAY)%VZ, &
      !                   particlesARRAY(NP_ARRAY)%EROT, particlesARRAY(NP_ARRAY)%EVIB, &
      !                   particlesARRAY(NP_ARRAY)%S_ID, particlesARRAY(NP_ARRAY)%IC, &
      !                   particlesARRAY(NP_ARRAY)%DTRIM, particlesARRAY(ID_REMOVE))
      !particlesARRAY(ID_REMOVE)%ID = particlesARRAY(NP_ARRAY)%ID
      ! Just do
      particlesARRAY(ID_REMOVE) = particlesARRAY(NP_ARRAY)
      ! Then put the last place to zeros and decrement counter
      !CALL INIT_PARTICLE(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-1,-1, 0.d0, particlesARRAY(NP_ARRAY)) ! Zero
      ! Just use
      particlesARRAY(NP_ARRAY)%X  = 0.d0
      particlesARRAY(NP_ARRAY)%Y  = 0.d0
      particlesARRAY(NP_ARRAY)%Z  = 0.d0
      particlesARRAY(NP_ARRAY)%VX = 0.d0
      particlesARRAY(NP_ARRAY)%VY = 0.d0
      particlesARRAY(NP_ARRAY)%VZ = 0.d0
      particlesARRAY(NP_ARRAY)%S_ID = -1
      particlesARRAY(NP_ARRAY)%ID = -1
      particlesARRAY(NP_ARRAY)%DUMP_TRAJ = .FALSE.
      

      NP_ARRAY = NP_ARRAY - 1

   END SUBROUTINE REMOVE_PARTICLE_ARRAY

END MODULE particle
