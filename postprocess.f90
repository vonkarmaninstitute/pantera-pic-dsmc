  MODULE postprocess

  USE global

  USE screen

  IMPLICIT NONE

  CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE GRID_AVG -> Adds timestep to cumulative average !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GRID_AVG

    IMPLICIT NONE

    INTEGER                            :: NCELLS, NSPECIES
    INTEGER                            :: JP, JC, JS, INDEX, TALLY

    REAL(KIND=8)                       :: KB

    REAL(KIND=8) :: DBLE_AVG_CUMULATED

 
    INTEGER, DIMENSION(:), ALLOCATABLE      :: TIMESTEP_NP

    INTEGER, DIMENSION(:), ALLOCATABLE      :: TIMESTEP_N
 
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

    INTEGER :: LENGTH

    KB = 1.38064852E-23

    LENGTH = NX*NY * N_SPECIES

    ALLOCATE(TIMESTEP_NP(LENGTH))
    
    ALLOCATE(TIMESTEP_N(LENGTH))

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
    
    TIMESTEP_NP = 0

    TIMESTEP_N = 0

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

    NCELLS = NX*NY
    NSPECIES = N_SPECIES

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


    DO INDEX = 1, NCELLS*NSPECIES
      IF (TIMESTEP_NP(INDEX) .GT. 0) THEN
        TIMESTEP_VX(INDEX) = TIMESTEP_VX(INDEX) / DBLE(TIMESTEP_NP(INDEX))
        TIMESTEP_VY(INDEX) = TIMESTEP_VY(INDEX) / DBLE(TIMESTEP_NP(INDEX))
        TIMESTEP_VZ(INDEX) = TIMESTEP_VZ(INDEX) / DBLE(TIMESTEP_NP(INDEX))
      END IF
    END DO

    ! Translational Temperature
    DO JP = 1, NP_PROC
      JC = particles(JP)%IC
      JS = particles(JP)%S_ID
      INDEX = JC+NCELLS*(JS-1)
      TIMESTEP_VX2(INDEX) = TIMESTEP_VX2(INDEX) + particles(JP)%VX*particles(JP)%VX
      TIMESTEP_VY2(INDEX) = TIMESTEP_VY2(INDEX) + particles(JP)%VY*particles(JP)%VY
      TIMESTEP_VZ2(INDEX) = TIMESTEP_VZ2(INDEX) + particles(JP)%VZ*particles(JP)%VZ
    END DO


    DO JC = 1, NCELLS
      DO JS = 1, NSPECIES
        INDEX = JC+NCELLS*(JS-1)
        IF (TIMESTEP_NP(INDEX) .GT. 0) THEN
          TIMESTEP_TTRX(INDEX) = SPECIES(JS)%MOLECULAR_MASS / KB * (TIMESTEP_VX2(INDEX) / DBLE(TIMESTEP_NP(INDEX)) &
                                - TIMESTEP_VX(INDEX)*TIMESTEP_VX(INDEX))
          TIMESTEP_TTRY(INDEX) = SPECIES(JS)%MOLECULAR_MASS / KB * (TIMESTEP_VY2(INDEX) / DBLE(TIMESTEP_NP(INDEX)) &
                                - TIMESTEP_VY(INDEX)*TIMESTEP_VY(INDEX))
          TIMESTEP_TTRZ(INDEX) = SPECIES(JS)%MOLECULAR_MASS / KB * (TIMESTEP_VZ2(INDEX) / DBLE(TIMESTEP_NP(INDEX)) &
                                - TIMESTEP_VZ(INDEX)*TIMESTEP_VZ(INDEX))
        END IF
      END DO
    END DO



    ! Collect data from all the processes
    IF (PROC_ID .EQ. 0) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_NP, NCELLS*NSPECIES, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VX,  NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VY,  NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_VZ,  NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_TTRX,  NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_TTRY,  NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(MPI_IN_PLACE,  TIMESTEP_TTRZ,  NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      TALLY = 0
      DO INDEX = 1, NCELLS*NSPECIES
        TALLY = TALLY + TIMESTEP_NP(INDEX)
      END DO
      WRITE(*,*) 'Total in postprocess', TALLY
    ELSE
      CALL MPI_REDUCE(TIMESTEP_NP,  TIMESTEP_NP,  NCELLS*NSPECIES, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TIMESTEP_VX,   TIMESTEP_VX,   NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TIMESTEP_VY,   TIMESTEP_VY,   NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TIMESTEP_VZ,   TIMESTEP_VZ,   NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TIMESTEP_TTRX,   TIMESTEP_TTRX,   NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TIMESTEP_TTRY,   TIMESTEP_TTRY,   NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TIMESTEP_TTRZ,   TIMESTEP_TTRZ,   NCELLS*NSPECIES, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    END IF

    ! Add to cumulated average
    DBLE_AVG_CUMULATED = DBLE(AVG_CUMULATED)
    AVG_NP =   (AVG_NP*DBLE_AVG_CUMULATED + DBLE(TIMESTEP_NP))/(AVG_CUMULATED + 1.)
    AVG_VX =   (AVG_VX*DBLE_AVG_CUMULATED + TIMESTEP_VX)/(AVG_CUMULATED + 1.)
    AVG_VY =   (AVG_VY*DBLE_AVG_CUMULATED + TIMESTEP_VY)/(AVG_CUMULATED + 1.)
    AVG_VZ =   (AVG_VZ*DBLE_AVG_CUMULATED + TIMESTEP_VZ)/(AVG_CUMULATED + 1.)
    AVG_TTRX = (AVG_TTRX*DBLE_AVG_CUMULATED + TIMESTEP_TTRX)/(AVG_CUMULATED + 1.)
    AVG_TTRY = (AVG_TTRY*DBLE_AVG_CUMULATED + TIMESTEP_TTRY)/(AVG_CUMULATED + 1.)
    AVG_TTRZ = (AVG_TTRZ*DBLE_AVG_CUMULATED + TIMESTEP_TTRZ)/(AVG_CUMULATED + 1.)
    AVG_TTR = (AVG_TTRX + AVG_TTRY + AVG_TTRZ) / 3.0

    
    AVG_CUMULATED = AVG_CUMULATED + 1

    
  END SUBROUTINE GRID_AVG



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE GRID_SAVE -> Saves cumulated average !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE GRID_SAVE
    
    IMPLICIT NONE

    CHARACTER*30                       :: string, file_name

    REAL(KIND=8), DIMENSION(NX+1)        :: XNODES
    REAL(KIND=8), DIMENSION(NY+1)        :: YNODES


    INTEGER                            :: NCELLS, NSPECIES
    INTEGER                            :: i, JS, FIRST, LAST



    NCELLS = NX*NY
    NSPECIES = N_SPECIES

    IF (PROC_ID .EQ. 0) THEN
      XNODES = 0
      DO i = 0, NX
         XNODES(i+1) = XMIN + (XMAX-XMIN)/NX*i 
      END DO

      YNODES = 0
      DO i = 0, NY
         YNODES(i+1) = YMIN + (YMAX-YMIN)/NY*i
      END DO

      ! DSMC flowfield file

      WRITE(file_name,'(A, I0, A)') 'dsmc_flowfield_', tID, '.vtk'

      OPEN(54321, FILE=file_name, ACCESS='SEQUENTIAL', FORM='formatted', STATUS='new')

      WRITE(54321,'(A)') '# vtk DataFile Version 3.0'
      WRITE(54321,'(A)') 'vtk output'
      WRITE(54321,'(A)') 'ASCII'
      WRITE(54321,'(A)') 'DATASET RECTILINEAR_GRID'
      WRITE(54321,'(A,I10,I10,I10)') 'DIMENSIONS', NX+1, NY+1, 1 

      WRITE(54321,'(A,I10,A7)') 'X_COORDINATES', NX+1, 'double'
      WRITE(54321,*) XNODES

      WRITE(54321,'(A,I10,A7)') 'Y_COORDINATES', NY+1, 'double'
      WRITE(54321,*) YNODES

      WRITE(54321,'(A,I10,A7)') 'Z_COORDINATES', 1, 'double'
      WRITE(54321,*) 0.

      WRITE(54321,'(A,I10)') 'CELL_DATA', NCELLS
      WRITE(54321,'(A,I10)') 'FIELD FieldData', 8*NSPECIES



      DO JS = 1, NSPECIES
        FIRST = 1 + (JS-1)*NCELLS
        LAST  = JS*NCELLS
      
        WRITE(string, *) 'number_particles_', SPECIES(JS)%NAME
        WRITE(54321,'(A,I10,I10,A7)') string, 1, NCELLS, 'double'
        WRITE(54321,*) AVG_NP(FIRST:LAST)

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

      END DO

      !WRITE(54321,'(A,I10)') 'POINT_DATA', (NX+1)*(NY+1)*1

      CLOSE(54321)

   END IF 

  END SUBROUTINE GRID_SAVE


  SUBROUTINE INIT_POSTPROCESS

    IMPLICIT NONE

    INTEGER :: LENGTH

    LENGTH = NX*NY * N_SPECIES
    
    ALLOCATE(AVG_NP(LENGTH))
    
    ALLOCATE(AVG_N(LENGTH))

    ALLOCATE(AVG_VX(LENGTH))
    ALLOCATE(AVG_VY(LENGTH))
    ALLOCATE(AVG_VZ(LENGTH))
    
    ALLOCATE(AVG_TTRX(LENGTH))
    ALLOCATE(AVG_TTRY(LENGTH))
    ALLOCATE(AVG_TTRZ(LENGTH))
    ALLOCATE(AVG_TTR(LENGTH))
    
    AVG_NP = 0

    AVG_N = 0

    AVG_VX = 0
    AVG_VY = 0
    AVG_VZ = 0

    AVG_TTRX = 0
    AVG_TTRY = 0
    AVG_TTRZ = 0
    AVG_TTR = 0

    AVG_CUMULATED = 0

  END SUBROUTINE INIT_POSTPROCESS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE GRID_RESET -> Resets cumulated average !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GRID_RESET

    IMPLICIT NONE

    AVG_NP = 0

    AVG_N = 0

    AVG_VX = 0
    AVG_VY = 0
    AVG_VZ = 0

    AVG_TTRX = 0
    AVG_TTRY = 0
    AVG_TTRZ = 0
    AVG_TTR = 0

    AVG_CUMULATED = 0

  END SUBROUTINE GRID_RESET


  SUBROUTINE CHECKS

    IMPLICIT NONE

    INTEGER                            :: JP, JS
 
    INTEGER                            :: TOT_NUM

    REAL(KIND=8), DIMENSION(3)         :: TOT_MOMENTUM
    REAL(KIND=8)                       :: TOT_KE


    TOT_NUM = 0
    TOT_MOMENTUM = 0
    TOT_KE = 0

    ! Number of particles
    DO JP = 1, NP_PROC
      TOT_NUM = TOT_NUM + 1

      JS = particles(JP)%S_ID

      ! Momentum
      TOT_MOMENTUM(1) = TOT_MOMENTUM(1) + SPECIES(JS)%MOLECULAR_MASS*particles(JP)%VX
      TOT_MOMENTUM(2) = TOT_MOMENTUM(2) + SPECIES(JS)%MOLECULAR_MASS*particles(JP)%VY
      TOT_MOMENTUM(3) = TOT_MOMENTUM(3) + SPECIES(JS)%MOLECULAR_MASS*particles(JP)%VZ

      ! Kinietic energy
      TOT_KE = TOT_KE + 0.5*SPECIES(JS)%MOLECULAR_MASS*(particles(JP)%VX**2+particles(JP)%VY**2+particles(JP)%VZ**2)

    END DO



    ! Collect data from all the processes and print it
    IF (PROC_ID .EQ. 0) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_NUM,      1, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_MOMENTUM, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(MPI_IN_PLACE,  TOT_KE,       1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      WRITE(*,*) ' '
      WRITE(*,*) 'Conservation checks:'
      WRITE(*,*) 'Number of particles:  ', TOT_NUM
      WRITE(*,*) 'Total momentum:       [ ', TOT_MOMENTUM, ' ] [kg m/s]'
      WRITE(*,*) 'Total kinetic energy: ', TOT_KE, ' [J]'
      WRITE(*,*) ' '

    ELSE
      CALL MPI_REDUCE(TOT_NUM,       TOT_NUM,      1, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TOT_MOMENTUM,  TOT_MOMENTUM, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(TOT_KE,        TOT_KE,       1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    END IF

  END SUBROUTINE CHECKS

END MODULE postprocess