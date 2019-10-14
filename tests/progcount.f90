PROGRAM countprog

  IMPLICIT NONE

  TYPE PARTICLE_DATA_STRUCTURE
     REAL(KIND=8) :: X, Y, Z        ! position
     REAL(KIND=8) :: VX, VY, VZ, EI ! velocities and internal energy
     INTEGER      :: IC             ! Cell index 
     INTEGER      :: S_ID           ! Species ID
  END TYPE PARTICLE_DATA_STRUCTURE

  TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: PARTICLES_ARRAY ! array
  INTEGER                      , DIMENSION(:), ALLOCATABLE :: PARTICLES_COUNT ! freq
  INTEGER                      , DIMENSION(:), ALLOCATABLE   :: disp1
  TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE   :: sorted
  INTEGER                                                    :: i1
  INTEGER                                                    :: IPROC

  ! ===========================================================
  ! Initialize particles array
  ALLOCATE(PARTICLES_ARRAY(12))
  DO i1 = 1, 12
    PARTICLES_ARRAY(i1)%X = 0.1
    PARTICLES_ARRAY(i1)%Y = 0.1
    PARTICLES_ARRAY(i1)%Z = 0.1
    PARTICLES_ARRAY(i1)%VX = 0.1
    PARTICLES_ARRAY(i1)%VY = 0.1
    PARTICLES_ARRAY(i1)%VZ = 0.1
  END DO

  PARTICLES_ARRAY(1)%IC = 1
  PARTICLES_ARRAY(2)%IC = 0
  PARTICLES_ARRAY(3)%IC = 0
  PARTICLES_ARRAY(4)%IC = 1
  PARTICLES_ARRAY(5)%IC = 2
  PARTICLES_ARRAY(6)%IC = 1
  PARTICLES_ARRAY(7)%IC = 2
  PARTICLES_ARRAY(8)%IC = 1
  PARTICLES_ARRAY(9)%IC = 0
  PARTICLES_ARRAY(10)%IC = 0
  PARTICLES_ARRAY(11)%IC = 1
  PARTICLES_ARRAY(12)%IC = 2

  ALLOCATE(PARTICLES_COUNT(12))
  PARTICLES_COUNT(1) = 4
  PARTICLES_COUNT(2) = 5
  PARTICLES_COUNT(3) = 3

  PRINT *
  PRINT *
  DO i1 = 1, 12
    WRITE(*,*) PARTICLES_ARRAY(i1)%IC
  END DO


  ! ===========================================================
  ! SORT ROUTINE

  ALLOCATE(sorted(SIZE(PARTICLES_ARRAY)))
  ALLOCATE(disp1(0:SIZE(PARTICLES_COUNT)-1))

  disp1 = PARTICLES_COUNT

  DO i1 = 1, 3 - 1
     disp1(i1) = disp1(i1) + disp1(i1-1)
  END DO

  DO i1 = SIZE(PARTICLES_ARRAY), 1, -1

     ! Find position of particle
     IPROC = PARTICLES_ARRAY(i1)%IC

     sorted(disp1(IPROC)) = PARTICLES_ARRAY(i1)
     disp1(IPROC)         = disp1(IPROC) - 1

  END DO

  PARTICLES_ARRAY = sorted

  DEALLOCATE(sorted)
  DEALLOCATE(disp1)


  ! PRINT
  PRINT*
  PRINT*
  DO i1 = 1, 12
    WRITE(*,*) PARTICLES_ARRAY(i1)%IC
  END DO

END PROGRAM
