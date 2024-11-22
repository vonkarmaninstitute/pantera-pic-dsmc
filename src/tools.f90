MODULE tools
   
USE mpi_common
USE mt19937_64
USE global
USE screen

CONTAINS 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION rf -> Pseudo-random number generator (RNG) !!!!!!!!!!!!!!!!
   ! shorthand for fortran builtin RANDOM_NUMBER function !!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   FUNCTION rf() RESULT(out)

      IMPLICIT NONE

      REAL(KIND=8) :: out

      out = genrand64_real1() ! From the Mersenne Twister module
      
      RETURN

   END FUNCTION rf

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ! FUNCTION rf -> Pseudo-random number generator (RNG) !!!!!!!!!!!!!!!!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   FUNCTION rf(idum) RESULT(out)
!       
!   IMPLICIT NONE
!       
!   INTEGER       :: idum
!   REAL(KIND=8)  :: rand, out  
! 
!   INTEGER,      PARAMETER :: MBIG=1000000000,MSEED=161803398,MZ=0 
!   REAL(KIND=8), PARAMETER :: FAC=1.d0/MBIG
! 
!   INTEGER                      :: i,ii,k
!   INTEGER, SAVE                :: iff,inext,inextp
!   INTEGER                      :: mj,mk
!   INTEGER, DIMENSION(55), SAVE :: ma 
!   
! 321  IF (idum < 0 .OR. iff == 0) THEN
! 
!      iff = 1
!      mj = MSEED - iabs(idum)
!      mj = MOD(mj,MBIG)
!      ma(55) = mj
!      mk = 1
!      
!      DO i = 1,54
! 
!        ii = MOD(21*i,55)
!        ma(ii) = mk
!        mk = mj - mk
! 
!        IF (mk < MZ) mk = mk + MBIG
!          
!        mj = ma(ii)
!         
!      END DO
!      
!      DO  k = 1,4
!         DO  i = 1,55
!            ma(i) = ma(i) - ma(1+MOD(i+30,55))
!            IF (ma(i) < MZ) ma(i) = ma(i) + MBIG
!         END DO
!      END DO
! 
!      inext = 0
!      inextp = 31
!      idum = 1
!   END IF
!   
!   inext = inext + 1
!   
!   IF(inext == 56) inext = 1      
!      
!   inextp = inextp + 1
!    
!   IF(inextp == 56) inextp = 1      
!      
!   mj =ma(inext) - ma(inextp)
!      
!   IF(mj < MZ) mj =mj + MBIG
!   
!   ma(inext) = mj
!   rand = mj * FAC
!   out = rand
! 
!   if (rand.le.1e-16) go to 321
!       
!   RETURN
!       
!   END FUNCTION rf  


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! SUBROUTINE MAXWELL -> Samples velocity and internal energy of      !
   ! particle following Maxwell-Boltzmann distribution at average       ! 
   ! velocity (UX, UY, UZ), translational temperatures (TX, TY, TZ)     !
   ! and rotational temperature TR.                                     ! 
   ! Takes RGAS = kB/Mmolecule[kg] as a parameter.                      !
   !                                                                    !
   ! Rotational energy is continuous and for diatomic species.          !
   ! It works also for monatomic, if you give Trot = 0.                 !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE MAXWELL(UX, UY, UZ, TX, TY, TZ, VX, VY, VZ, M)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN)    :: UX, UY, UZ, TX, TY, TZ
      REAL(KIND=8), INTENT(IN)    :: M ! Molecular mass
      REAL(KIND=8), INTENT(INOUT) :: VX, VY, VZ

      INTEGER                     :: I
      REAL(KIND=8)                :: PI2
      REAL(KIND=8)                :: R, R1, RO, TETA, BETA
      REAL(KIND=8), DIMENSION(3)  :: VEL, TT

      PI2  = 2.*PI

      TT(1) = TX
      TT(2) = TY
      TT(3) = TZ



      DO I = 1,3
         IF (TT(I) == 0.d0) THEN
            VEL(I) = 0
         ELSE

            ! Step 1.
            R = rf()
               
            TETA = PI2*R

            ! Step 2.
            
            BETA = 1./SQRT(2.*KB/M*TT(I))

            ! R goes from 0 to 1 included. Remove the extremes
            ! or the log() will explode
            R1 = rf()
            DO WHILE (R1 < 1.0D-13)
               R1 = rf()
            END DO

            RO = SQRT(-LOG(R1))/BETA ! The random number here shouldn't be correlated to the one for teta!!

            VEL(I) = RO*SIN(TETA)

         END IF
      END DO

      ! Step 3.

      VX = UX + VEL(1)
      VY = UY + VEL(2)
      VZ = UZ + VEL(3)

      RETURN

   END SUBROUTINE MAXWELL


   !!!!! TODO: Implement Kappa distribution for particles
   ! SUBROUTINE KAPPA(UX, UY, UZ, TX, TY, TZ, VX, VY, VZ, M)

   !    IMPLICIT NONE

   !    REAL(KIND=8), INTENT(IN)    :: UX, UY, UZ, TX, TY, TZ
   !    REAL(KIND=8), INTENT(IN)    :: M ! Molecular mass
   !    REAL(KIND=8), INTENT(INOUT) :: VX, VY, VZ

   !    INTEGER                     :: I
   !    REAL(KIND=8)                :: PI2
   !    REAL(KIND=8)                :: R, R1, RO, TETA, BETA
   !    REAL(KIND=8), DIMENSION(3)  :: VEL, TT
   !    REAL(KIND=8) :: DELTA

   !    PI2  = 2.*PI

   !    TT(1) = TX
   !    TT(2) = TY
   !    TT(3) = TZ


   !    !!!!! KAPPA DISTRIBUTION !!!!!

   !    DO I = 1,3
   !       ! Step 1.
   !       R = rf()
             
   !       TETA = PI2*R

   !       ! Step 2.
   !       ! R goes from 0 to 1 included. Remove the extremes
   !       R1 = rf()
   !          DO WHILE (R1 < 1.0D-13)
   !             R1 = rf()
   !          END DO

   !       BETA = 1./SQRT(2.*KB*TT(I)/M*(KAPPA_C-3./2.))
   !       DELTA = (1-GAMMA(KAPPA_C-1./2.)/GAMMA(KAPPA_C+1./2.)*(KAPPA_C-1./2.)*(1-R1))**(-1./(KAPPA_C-1./2.))

   !       RO = SQRT(DELTA-1)/BETA
   !       VEL(I) = RO*SIN(TETA)
   !    END DO

   !    ! Step 3.

   !    VX = UX + VEL(1)
   !    VY = UY + VEL(2)
   !    VZ = UZ + VEL(3)

   !    RETURN

   ! END SUBROUTINE KAPPA


   SUBROUTINE INTERNAL_ENERGY(DOF, TEMP, EI)

      IMPLICIT NONE

      INTEGER, INTENT(IN)       :: DOF
      REAL(KIND=8), INTENT(IN)  :: TEMP
      REAL(KIND=8), INTENT(OUT) :: EI
      REAL(KIND=8)              :: R, X, PARAM

      IF (DOF .EQ. 0) THEN

         EI = 0.
      
      ELSE IF (DOF .EQ. 2) THEN

         R = rf()
         DO WHILE (R < 1.0D-13)
            R = rf()
         END DO
         EI = -LOG(R)*KB*TEMP

      ELSE
         
         IF (DOF .EQ. 1) THEN
            PARAM = 5.0
         ELSE
            PARAM = 0.01
         END IF
         R = rf()
         X = -LOG(R)
         R = rf()
         DO WHILE ( R .GT. (X/PARAM)**(DOF/2.-1.) )
            R = rf()
            X = -LOG(R)
            R = rf()
         END DO
         EI = X*KB*TEMP

      END IF

   END SUBROUTINE INTERNAL_ENERGY



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DUMP_PARTICLES_SCREEN -> dumps particle properties on the screen !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DUMP_PARTICLES_SCREEN(TIMESTEP)

      ! This subroutine dumps particles at the screen. DUMPS MAY OVERWRITE, SINCE STDOUT IS BUFFERED!!!!
      ! Doesn't matter if I loop on processes!!!!
      ! 
      ! OLD OLD OL A loop is performed on the processes, making sure that
      ! OLD OLD OL processes don't do it simultaneously, messing up the dump.

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      INTEGER :: IP, IPROC

      ! Loop on MPI threads
      DO IPROC = 0,N_MPI_THREADS-1

         IF (PROC_ID == IPROC) THEN 
            PRINT*, "I am proc: ", PROC_ID

            DO IP = 1, NP_PROC
               WRITE(*,*) 'DUMP PARTICLES: ', TIMESTEP, particles(IP)%X, particles(IP)%Y, particles(IP)%Z, PROC_ID
            END DO

         END IF

         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) ! Sync
      END DO
   
   END SUBROUTINE DUMP_PARTICLES_SCREEN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DUMP_PARTICLES_FILE -> dumps particle properties to file !!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DUMP_PARTICLES_FILE(TIMESTEP)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename
      INTEGER :: IP

      WRITE(filename, "(A,A,I0.5,A6,I0.8)") TRIM(ADJUSTL(PARTDUMP_SAVE_PATH)), "proc_", PROC_ID, "_time_", TIMESTEP ! Compose filename

      ! Open file for writing
      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(10, FILE=filename, ACCESS='SEQUENTIAL', FORM='UNFORMATTED', STATUS='NEW', CONVERT='BIG_ENDIAN', RECL=80)
         DO IP = 1, NP_PROC
            IF (PARTDUMP_FRACSAMPLE < 1) THEN
               IF (rf() > PARTDUMP_FRACSAMPLE) CYCLE
            END IF
            WRITE(10) particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
            particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ, particles(IP)%EROT, particles(IP)%EVIB, &
            particles(IP)%S_ID, particles(IP)%IC, particles(IP)%DTRIM
         END DO
         CLOSE(10)
      ELSE
         OPEN(10, FILE=filename )
         !WRITE(10,*) '% X | Y | Z | VX | VY | VZ | EROT | EVIB | S_ID | IC | DTRIM'
         DO IP = 1, NP_PROC
            IF (PARTDUMP_FRACSAMPLE < 1) THEN
               IF (rf() > PARTDUMP_FRACSAMPLE) CYCLE
            END IF
            WRITE(10,*) particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
            particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ, particles(IP)%EROT, particles(IP)%EVIB, &
            particles(IP)%S_ID, particles(IP)%IC, particles(IP)%DTRIM
         END DO
         CLOSE(10)
      END IF



   END SUBROUTINE DUMP_PARTICLES_FILE


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DUMP_PARTICLES_FILE -> dumps particle properties to file !!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DUMP_BOUNDARY_PARTICLES_FILE(TIMESTEP)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename
      INTEGER :: IP


      ! Dump particles that hit a boundary to file
      WRITE(filename, "(A,A,I0.5)") TRIM(ADJUSTL(PARTDUMP_SAVE_PATH)), "bound_proc_", PROC_ID ! Compose filename

      ! Open file for writing
      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(28479, FILE=filename, ACCESS='SEQUENTIAL', POSITION='APPEND', FORM='UNFORMATTED', &
         STATUS='UNKNOWN', CONVERT='BIG_ENDIAN', RECL=84)
         DO IP = 1, NP_DUMP_PROC
            WRITE(28479) TIMESTEP, part_dump(IP)%X, part_dump(IP)%Y, part_dump(IP)%Z, &
            part_dump(IP)%VX, part_dump(IP)%VY, part_dump(IP)%VZ, part_dump(IP)%EROT, part_dump(IP)%EVIB, &
            part_dump(IP)%S_ID, part_dump(IP)%IC, part_dump(IP)%DTRIM
         END DO
         CLOSE(28479)
      ELSE
         OPEN(28479, FILE=filename )
         !WRITE(10,*) '% TIMESTEP | X | Y | Z | VX | VY | VZ | EROT | EVIB | S_ID | IPG | DTRIM'
         DO IP = 1, NP_DUMP_PROC
            WRITE(28479,*) TIMESTEP, part_dump(IP)%X, part_dump(IP)%Y, part_dump(IP)%Z, &
            part_dump(IP)%VX, part_dump(IP)%VY, part_dump(IP)%VZ, part_dump(IP)%EROT, part_dump(IP)%EVIB, &
            part_dump(IP)%S_ID, part_dump(IP)%IC, part_dump(IP)%DTRIM
         END DO
         CLOSE(28479)
      END IF

   END SUBROUTINE DUMP_BOUNDARY_PARTICLES_FILE

   SUBROUTINE DUMP_FORCE_FILE(TIMESTEP)

      IMPLICIT NONE

      REAL(KIND=8) :: CURRENT_TIME
      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename

      REAL(KIND=8), DIMENSION(3) :: DUMP_FORCE_DIRECT, DUMP_FORCE_INDIRECT

      DUMP_FORCE_DIRECT = FORCE_DIRECT
      DUMP_FORCE_INDIRECT = FORCE_INDIRECT

      CURRENT_TIME = TIMESTEP*DT


      IF (PROC_ID .EQ. 0) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,  DUMP_FORCE_DIRECT, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(MPI_IN_PLACE,  DUMP_FORCE_INDIRECT, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         
         WRITE(filename, "(A,A,I0.5)") TRIM(ADJUSTL(FLOWFIELD_SAVE_PATH)), "dump_force" ! Compose filename
         OPEN(54331, FILE=filename, POSITION='append', STATUS='unknown', ACTION='write')
         WRITE(54331,*) CURRENT_TIME, DUMP_FORCE_DIRECT, DUMP_FORCE_INDIRECT
         CLOSE(54331)

      ELSE
         CALL MPI_REDUCE(DUMP_FORCE_DIRECT,  DUMP_FORCE_DIRECT, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(DUMP_FORCE_INDIRECT,  DUMP_FORCE_INDIRECT, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF

      FORCE_DIRECT = 0.d0
      FORCE_INDIRECT = 0.d0
   END SUBROUTINE DUMP_FORCE_FILE


   SUBROUTINE READ_PARTICLES_FILE(TIMESTEP)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename
      INTEGER :: ios

      REAL(KIND=8) :: XP, YP, ZP, VX, VY, VZ, EROT, EVIB, DTRIM
      INTEGER      :: S_ID, IC
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
      character(len=100) :: iomsg
      WRITE(filename, "(A,A,I0.5,A6,I0.8)") TRIM(ADJUSTL(PARTDUMP_SAVE_PATH)), "proc_", PROC_ID, "_time_", TIMESTEP ! Compose filename

      ! Open file for reading
      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(1010, FILE=filename, ACCESS='SEQUENTIAL', FORM='UNFORMATTED', STATUS='OLD', &
         CONVERT='BIG_ENDIAN', RECL=80, IOSTAT=ios, IOMSG=iomsg)

         IF (ios .NE. 0) THEN
            WRITE(*,*) 'iomsg was ', iomsg
            CALL ERROR_ABORT('Attention, particle restart file not found! ABORTING.')
         ENDIF

         DO
            READ(1010, IOSTAT=ios) XP, YP, ZP, VX, VY, VZ, EROT, EVIB, S_ID, IC, DTRIM

            IF (ios < 0) EXIT
            CALL INIT_PARTICLE(XP,YP,ZP,VX,VY,VZ,EROT,EVIB,S_ID,IC,DT, particleNOW) ! Save in particle
            CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
         END DO

         CLOSE(1010)
      ELSE
         OPEN(1010, FILE=filename, STATUS='OLD', IOSTAT=ios)
         
         IF (ios .NE. 0) THEN
            CALL ERROR_ABORT('Attention, particle restart file not found! ABORTING.')
         ENDIF

         DO
            READ(1010,*,IOSTAT=ios) XP, YP, ZP, VX, VY, VZ, EROT, EVIB, S_ID, IC, DTRIM
            IF (ios < 0) EXIT
            CALL INIT_PARTICLE(XP,YP,ZP,VX,VY,VZ,EROT,EVIB,S_ID,IC,DT, particleNOW) ! Save in particle
            CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
         END DO

         CLOSE(1010)
      END IF



   END SUBROUTINE READ_PARTICLES_FILE




   SUBROUTINE READ_INJECT_FILE

      IMPLICIT NONE

      CHARACTER(LEN=512)  :: filename
      INTEGER :: ios

      REAL(KIND=8) :: XP, YP, ZP, VX, VY, VZ, EROT, EVIB, DTRIM
      INTEGER      :: S_ID, IC, TIMESTEP
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW
      character(len=100) :: iomsg
      WRITE(filename, "(A,A,A,I0.5)") TRIM(ADJUSTL(PARTDUMP_SAVE_PATH)), TRIM(ADJUSTL(INJECT_FILENAME)), "proc_", PROC_ID ! Compose filename

      ! Open file for reading
      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(1030, FILE=filename, ACCESS='SEQUENTIAL', FORM='UNFORMATTED', STATUS='OLD', &
         CONVERT='BIG_ENDIAN', RECL=84, IOSTAT=ios, IOMSG=iomsg)

         IF (ios .NE. 0) THEN
            WRITE(*,*) 'iomsg was ', iomsg
            CALL ERROR_ABORT('Attention, particle inject file not found! ABORTING.')
         ENDIF

         DO
            READ(1030, IOSTAT=ios) TIMESTEP, XP, YP, ZP, VX, VY, VZ, EROT, EVIB, S_ID, IC, DTRIM

            IF (ios < 0) EXIT
            CALL INIT_PARTICLE(XP,YP,ZP,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTRIM, particleNOW) ! Save in particle
            CALL ADD_PARTICLE_ARRAY(particleNOW, NP_INJECT_PROC, part_inject) ! Add particle to local array
         END DO

         CLOSE(1030)
      ELSE
         OPEN(1030, FILE=filename, STATUS='OLD', IOSTAT=ios)
         
         IF (ios .NE. 0) THEN
            CALL ERROR_ABORT('Attention, particle inject file not found! ABORTING.')
         ENDIF

         DO
            READ(1030,*,IOSTAT=ios) TIMESTEP, XP, YP, ZP, VX, VY, VZ, EROT, EVIB, S_ID, IC, DTRIM
            IF (ios < 0) EXIT
            CALL INIT_PARTICLE(XP,YP,ZP,VX,VY,VZ,EROT,EVIB,S_ID,IC,DTRIM, particleNOW) ! Save in particle
            CALL ADD_PARTICLE_ARRAY(particleNOW, NP_INJECT_PROC, part_inject) ! Add particle to local array
         END DO

         CLOSE(1030)
      END IF



   END SUBROUTINE READ_INJECT_FILE



   SUBROUTINE REASSIGN_PARTICLES_TO_CELLS_3D

      IMPLICIT NONE

      REAL(KIND=8) :: DX, DY, DZ, CELLXMIN, CELLXMAX, CELLYMIN, CELLYMAX, CELLZMIN, CELLZMAX
      REAL(KIND=8) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, XP, YP, ZP, PSIP
      INTEGER :: I, J, K, IC, IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, V1, V2, V3, V4
      INTEGER :: N_PARTITIONS_X, N_PARTITIONS_Y, N_PARTITIONS_Z, IP, VP, CELL, IDX
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NINGRIDCELL, IOFTETRA
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDTETRA
      LOGICAL :: INSIDE, INSIDE_DOMAIN
      !REAL(KIND=8) :: MINABSPSI
      LOGICAL, DIMENSION(:), ALLOCATABLE :: REMOVE_PART

      N_PARTITIONS_X = 100
      N_PARTITIONS_Y = 100
      N_PARTITIONS_Z = 100

      ALLOCATE(REMOVE_PART(NP_PROC))
      REMOVE_PART = .FALSE.

      ! First we assign the unstructured cells to the cells of a cartesian grid, which is easily indexable.
      IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 3) THEN

         DX = (XMAX-XMIN)/N_PARTITIONS_X
         DY = (YMAX-YMIN)/N_PARTITIONS_Y
         DZ = (ZMAX-ZMIN)/N_PARTITIONS_Z
         
         ALLOCATE(NINGRIDCELL(N_PARTITIONS_Z, N_PARTITIONS_Y, N_PARTITIONS_X))
         NINGRIDCELL = 0

         DO IC = 1, NCELLS
            V1 = U3D_GRID%CELL_NODES(1,IC)
            V2 = U3D_GRID%CELL_NODES(2,IC)
            V3 = U3D_GRID%CELL_NODES(3,IC)
            V4 = U3D_GRID%CELL_NODES(4,IC)
   
            X1 = U3D_GRID%NODE_COORDS(1, V1)
            X2 = U3D_GRID%NODE_COORDS(1, V2)
            X3 = U3D_GRID%NODE_COORDS(1, V3)
            X4 = U3D_GRID%NODE_COORDS(1, V4)
            Y1 = U3D_GRID%NODE_COORDS(2, V1)
            Y2 = U3D_GRID%NODE_COORDS(2, V2)
            Y3 = U3D_GRID%NODE_COORDS(2, V3)
            Y4 = U3D_GRID%NODE_COORDS(2, V4)
            Z1 = U3D_GRID%NODE_COORDS(3, V1)
            Z2 = U3D_GRID%NODE_COORDS(3, V2)
            Z3 = U3D_GRID%NODE_COORDS(3, V3)
            Z4 = U3D_GRID%NODE_COORDS(3, V4)

            CELLXMIN = MIN(X1, X2, X3, X4)
            CELLXMAX = MAX(X1, X2, X3, X4)
            CELLYMIN = MIN(Y1, Y2, Y3, Y4)
            CELLYMAX = MAX(Y1, Y2, Y3, Y4)
            CELLZMIN = MIN(Z1, Z2, Z3, Z4)
            CELLZMAX = MAX(Z1, Z2, Z3, Z4)

            IMIN = INT((CELLXMIN - XMIN)/DX) + 1
            IMAX = INT((CELLXMAX - XMIN)/DX) + 1
            JMIN = INT((CELLYMIN - YMIN)/DY) + 1
            JMAX = INT((CELLYMAX - YMIN)/DY) + 1
            KMIN = INT((CELLZMIN - ZMIN)/DZ) + 1
            KMAX = INT((CELLZMAX - ZMIN)/DZ) + 1

            DO I = IMIN, IMAX
               DO J = JMIN, JMAX
                  DO K = KMIN, KMAX
                     NINGRIDCELL(K,J,I) = NINGRIDCELL(K,J,I) + 1
                  END DO
               END DO
            END DO

         END DO

         ALLOCATE(IOFTETRA(N_PARTITIONS_Z, N_PARTITIONS_Y, N_PARTITIONS_X))
         IOFTETRA = -1

         IDX = 1
         DO I = 1, N_PARTITIONS_X
            DO J = 1, N_PARTITIONS_Y
               DO K = 1, N_PARTITIONS_Z
                  IF (NINGRIDCELL(K,J,I) .NE. 0) THEN
                     IOFTETRA(K,J,I) = IDX
                     IDX = IDX + NINGRIDCELL(K,J,I)
                  END IF
               END DO
            END DO
         END DO

         ALLOCATE(INDTETRA(SUM(NINGRIDCELL)))
         NINGRIDCELL = 0

         DO IC = 1, NCELLS
            V1 = U3D_GRID%CELL_NODES(1,IC)
            V2 = U3D_GRID%CELL_NODES(2,IC)
            V3 = U3D_GRID%CELL_NODES(3,IC)
            V4 = U3D_GRID%CELL_NODES(4,IC)
   
            X1 = U3D_GRID%NODE_COORDS(1, V1)
            X2 = U3D_GRID%NODE_COORDS(1, V2)
            X3 = U3D_GRID%NODE_COORDS(1, V3)
            X4 = U3D_GRID%NODE_COORDS(1, V4)
            Y1 = U3D_GRID%NODE_COORDS(2, V1)
            Y2 = U3D_GRID%NODE_COORDS(2, V2)
            Y3 = U3D_GRID%NODE_COORDS(2, V3)
            Y4 = U3D_GRID%NODE_COORDS(2, V4)
            Z1 = U3D_GRID%NODE_COORDS(3, V1)
            Z2 = U3D_GRID%NODE_COORDS(3, V2)
            Z3 = U3D_GRID%NODE_COORDS(3, V3)
            Z4 = U3D_GRID%NODE_COORDS(3, V4)

            CELLXMIN = MIN(X1, X2, X3, X4)
            CELLXMAX = MAX(X1, X2, X3, X4)
            CELLYMIN = MIN(Y1, Y2, Y3, Y4)
            CELLYMAX = MAX(Y1, Y2, Y3, Y4)
            CELLZMIN = MIN(Z1, Z2, Z3, Z4)
            CELLZMAX = MAX(Z1, Z2, Z3, Z4)

            IMIN = INT((CELLXMIN - XMIN)/DX) + 1
            IMAX = INT((CELLXMAX - XMIN)/DX) + 1
            JMIN = INT((CELLYMIN - YMIN)/DY) + 1
            JMAX = INT((CELLYMAX - YMIN)/DY) + 1
            KMIN = INT((CELLZMIN - ZMIN)/DZ) + 1
            KMAX = INT((CELLZMAX - ZMIN)/DZ) + 1


            DO I = IMIN, IMAX
               DO J = JMIN, JMAX
                  DO K = KMIN, KMAX
                     INDTETRA(IOFTETRA(K,J,I) + NINGRIDCELL(K,J,I)) = IC
                     NINGRIDCELL(K,J,I) = NINGRIDCELL(K,J,I) + 1
                  END DO
               END DO
            END DO

         END DO

         DO IP = 1, NP_PROC
            XP = particles(IP)%X
            YP = particles(IP)%Y
            ZP = particles(IP)%Z
            
            IF (XP < XMIN .OR. XP > XMAX .OR. YP < YMIN .OR. YP > YMAX .OR.ZP < ZMIN .OR. ZP > ZMAX) THEN
               WRITE(*,*) 'Error during restart! Particle found ouside the bounds.'
               WRITE(*,*) 'Particle position [x, y, z]=[ ', XP, ', ', YP, ', ', ZP, ']'
               WRITE(*,*) 'Particle species: ', SPECIES(particles(IP)%S_ID)%NAME
               REMOVE_PART(IP) = .TRUE.
               CYCLE
            END IF

            I = INT((XP-XMIN)/DX) + 1
            J = INT((YP-YMIN)/DY) + 1
            K = INT((ZP-ZMIN)/DZ) + 1

            INSIDE_DOMAIN = .FALSE.
            DO CELL = IOFTETRA(K,J,I), IOFTETRA(K,J,I)+NINGRIDCELL(K,J,I)-1
               !MINABSPSI = 1.d100
               IC = INDTETRA(CELL)
               ! Check if particle IP (XP, YP) is in unstructured cell IC.
               INSIDE = .TRUE.
               DO VP = 1, 4
                  PSIP = XP*U3D_GRID%BASIS_COEFFS(1,VP,IC) + YP*U3D_GRID%BASIS_COEFFS(2,VP,IC) &
                       + ZP*U3D_GRID%BASIS_COEFFS(3,VP,IC) + U3D_GRID%BASIS_COEFFS(4,VP,IC)
                  IF (PSIP < 0) INSIDE = .FALSE.
                  !IF (ABS(PSIP) < MINABSPSI) MINABSPSI = ABS(PSIP)
               END DO
               IF (INSIDE) THEN
                  !IF (particles(IP)%IC .NE. IC) WRITE(*,*) 'Particle has been found in a different cell! IC=', IC, &
                  !' instead of ', particles(IP)%IC, '. Minimum ABS(PSI) was ', MINABSPSI
                  INSIDE_DOMAIN = .TRUE.
                  particles(IP)%IC = IC
                  EXIT
               END IF
            END DO
            IF (.NOT. INSIDE_DOMAIN) THEN
               WRITE(*,*) 'Error during restart! Particle found ouside the domain.'
               WRITE(*,*) 'Particle position [x, y, z]=[ ', XP, ', ', YP, ', ', ZP, ']'
               WRITE(*,*) 'Particle species: ', SPECIES(particles(IP)%S_ID)%NAME
               REMOVE_PART(IP) = .TRUE.
            END IF
         END DO

      ELSE
         CALL ERROR_ABORT('Not implemented!')
      END IF




      IP = NP_PROC
      DO WHILE (IP .GE. 1)
         ! Is particle IP out of the domain? Then remove it!
         IF (REMOVE_PART(IP)) THEN
            CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)
         END IF
         IP = IP - 1
      END DO

      DEALLOCATE(REMOVE_PART)


   END SUBROUTINE REASSIGN_PARTICLES_TO_CELLS_3D



   SUBROUTINE REASSIGN_PARTICLES_TO_CELLS_2D

      IMPLICIT NONE

      REAL(KIND=8) :: DX, DY, CELLXMIN, CELLXMAX, CELLYMIN, CELLYMAX
      REAL(KIND=8) :: X1, X2, X3, Y1, Y2, Y3, XP, YP, PSIP
      INTEGER :: I, J, K, IC, IMIN, IMAX, JMIN, JMAX, V1, V2, V3, N_PARTITIONS_X, N_PARTITIONS_Y, IP, VP
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NINGRIDCELL
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: TRISINGRID
      LOGICAL :: INSIDE, CELLFOUND
      REAL(KIND=8) :: MINABSPSI
      LOGICAL, DIMENSION(:), ALLOCATABLE :: REMOVE_PART

      N_PARTITIONS_X = 100
      N_PARTITIONS_Y = 100

      ALLOCATE(REMOVE_PART(NP_PROC))
      REMOVE_PART = .FALSE.

      ! First we assign the unstructured cells to the cells of a cartesian grid, which is easily indexable.
      IF (GRID_TYPE == UNSTRUCTURED .AND. DIMS == 2) THEN

         DX = (XMAX-XMIN)/N_PARTITIONS_X
         DY = (YMAX-YMIN)/N_PARTITIONS_Y
         
         ALLOCATE(NINGRIDCELL(N_PARTITIONS_Y, N_PARTITIONS_X))
         NINGRIDCELL = 0

         DO IC = 1, NCELLS
            V1 = U2D_GRID%CELL_NODES(1,IC)
            V2 = U2D_GRID%CELL_NODES(2,IC)
            V3 = U2D_GRID%CELL_NODES(3,IC)
   
            X1 = U2D_GRID%NODE_COORDS(1, V1)
            X2 = U2D_GRID%NODE_COORDS(1, V2)
            X3 = U2D_GRID%NODE_COORDS(1, V3)
            Y1 = U2D_GRID%NODE_COORDS(2, V1)
            Y2 = U2D_GRID%NODE_COORDS(2, V2)
            Y3 = U2D_GRID%NODE_COORDS(2, V3)

            CELLXMIN = MIN(X1, X2, X3)
            CELLXMAX = MAX(X1, X2, X3)
            CELLYMIN = MIN(Y1, Y2, Y3)
            CELLYMAX = MAX(Y1, Y2, Y3)

            IMIN = INT((CELLXMIN - XMIN)/DX) + 1
            IMAX = INT((CELLXMAX - XMIN)/DX) + 1
            JMIN = INT((CELLYMIN - YMIN)/DY) + 1
            JMAX = INT((CELLYMAX - YMIN)/DY) + 1

            DO I = IMIN, IMAX
               DO J = JMIN, JMAX
                  NINGRIDCELL(J,I) = NINGRIDCELL(J,I) + 1
               END DO
            END DO

         END DO

         ALLOCATE(TRISINGRID(MAXVAL(NINGRIDCELL), N_PARTITIONS_Y, N_PARTITIONS_X))
         NINGRIDCELL = 0

         DO IC = 1, NCELLS
            V1 = U2D_GRID%CELL_NODES(1,IC)
            V2 = U2D_GRID%CELL_NODES(2,IC)
            V3 = U2D_GRID%CELL_NODES(3,IC)
   
            X1 = U2D_GRID%NODE_COORDS(1, V1)
            X2 = U2D_GRID%NODE_COORDS(1, V2)
            X3 = U2D_GRID%NODE_COORDS(1, V3)
            Y1 = U2D_GRID%NODE_COORDS(2, V1)
            Y2 = U2D_GRID%NODE_COORDS(2, V2)
            Y3 = U2D_GRID%NODE_COORDS(2, V3)

            CELLXMIN = MIN(X1, X2, X3)
            CELLXMAX = MAX(X1, X2, X3)
            CELLYMIN = MIN(Y1, Y2, Y3)
            CELLYMAX = MAX(Y1, Y2, Y3)

            IMIN = INT((CELLXMIN - XMIN)/DX) + 1
            IMAX = INT((CELLXMAX - XMIN)/DX) + 1
            JMIN = INT((CELLYMIN - YMIN)/DY) + 1
            JMAX = INT((CELLYMAX - YMIN)/DY) + 1

            DO I = IMIN, IMAX
               DO J = JMIN, JMAX
                  NINGRIDCELL(J,I) = NINGRIDCELL(J,I) + 1
                  TRISINGRID(NINGRIDCELL(J,I), J, I) = IC
               END DO
            END DO

         END DO

         DO IP = 1, NP_PROC
            XP = particles(IP)%X
            YP = particles(IP)%Y
            IF (XP < XMIN .OR. XP > XMAX .OR. YP < YMIN .OR. YP > YMAX) THEN
               WRITE(*,*) 'Particle with ID ', particles(IP)%ID, 'is out of the domain. Deleting it.'
               REMOVE_PART(IP) = .TRUE.
               CYCLE
            END IF
            I = INT((XP-XMIN)/DX) + 1
            J = INT((YP-YMIN)/DY) + 1
            CELLFOUND = .FALSE.
            DO K = 1, NINGRIDCELL(J,I)
               MINABSPSI = 1.d100
               IC = TRISINGRID(K, J, I)
               ! Check if particle IP (XP, YP) is in unstructured cell IC.
               INSIDE = .TRUE.
               DO VP = 1, 3
                  PSIP = XP*U2D_GRID%BASIS_COEFFS(1,VP,IC) + YP*U2D_GRID%BASIS_COEFFS(2,VP,IC) + U2D_GRID%BASIS_COEFFS(3,VP,IC)
                  IF (PSIP < 0) INSIDE = .FALSE.
                  IF (ABS(PSIP) < MINABSPSI) MINABSPSI = ABS(PSIP)
               END DO
               IF (INSIDE) THEN
                  IF (particles(IP)%IC .NE. IC) WRITE(*,*) 'Particle with ID ', particles(IP)%ID, &
                  ' has been found in a different cell! IC=', IC, &
                  ' instead of ', particles(IP)%IC, '. Minimum ABS(PSI) was ', MINABSPSI
                  particles(IP)%IC = IC
                  CELLFOUND = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT. CELLFOUND) THEN
               WRITE(*,*) 'Particle with ID ', particles(IP)%ID, 'is out of the domain. Deleting it.'
               REMOVE_PART(IP) = .TRUE.
            END IF
         END DO

         IP = NP_PROC
         DO WHILE (IP .GE. 1)
            IF (REMOVE_PART(IP)) CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)
            IP = IP - 1
         END DO

      ELSE
         CALL ERROR_ABORT('Not implemented!')
      END IF

      DEALLOCATE(REMOVE_PART)

   END SUBROUTINE REASSIGN_PARTICLES_TO_CELLS_2D



   SUBROUTINE DUPLICATE_PARTICLES

      IMPLICIT NONE

      INTEGER :: IP, NP_PROC_INITIAL
      TYPE(PARTICLE_DATA_STRUCTURE) :: particleNOW

      NP_PROC_INITIAL = NP_PROC
      DO IP = 1, NP_PROC
         CALL INIT_PARTICLE(particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
         particles(IP)%VX + 10.d0*rf()-5.d0, &
         particles(IP)%VY + 10.d0*rf()-5.d0, &
         particles(IP)%VZ + 10.d0*rf()-5.d0, &
         particles(IP)%EROT, particles(IP)%EVIB, particles(IP)%S_ID, &
         particles(IP)%IC, particles(IP)%DTRIM, particleNOW) ! Save in particle
         CALL ADD_PARTICLE_ARRAY(particleNOW, NP_PROC, particles) ! Add particle to local array
      END DO

   END SUBROUTINE DUPLICATE_PARTICLES

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DUMP_TRAJECTORY_FILE -> dumps particle trajectory to file !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DUMP_TRAJECTORY_FILE(TIMESTEP)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename
      INTEGER :: IP

      DO IP = 1, NP_PROC
         IF (particles(IP)%DUMP_TRAJ) THEN
            !WRITE(*,*) 'Writing trajectory file for particle with ID ', particles(IP)%ID
            WRITE(filename, "(A,A,I0.15)") TRIM(ADJUSTL(TRAJDUMP_SAVE_PATH)), "trajectory_", particles(IP)%ID ! Compose filename
            ! Open file for writing
            IF (BOOL_BINARY_OUTPUT) THEN
               OPEN(1610, FILE=filename, ACCESS='SEQUENTIAL', POSITION='APPEND', FORM='UNFORMATTED', &
               STATUS='UNKNOWN', CONVERT='BIG_ENDIAN', RECL=56)
               WRITE(1610) TIMESTEP, particles(IP)%S_ID, particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
               particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ
               CLOSE(1610)
            ELSE
               OPEN(1610, FILE=filename, POSITION='APPEND')
               WRITE(1610,*) TIMESTEP, particles(IP)%S_ID, particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
               particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ
               CLOSE(1610)
            END IF   
         END IF
      END DO

   END SUBROUTINE DUMP_TRAJECTORY_FILE

   SUBROUTINE DUMP_FLUXES_FILE(TIMESTEP)

      ! This subroutine dumps particles at the screen. DUMPS MAY OVERWRITE, SINCE STDOUT IS BUFFERED!!!!
      ! Doesn't matter if I loop on processes!!!!
      ! 
      ! OLD OLD OL A loop is performed on the processes, making sure that
      ! OLD OLD OL processes don't do it simultaneously, messing up the dump.

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename

      IF (PROC_ID .EQ. 0) THEN
         WRITE(filename, "(A,A)") TRIM(ADJUSTL(FLUXDUMP_SAVE_PATH)), "fluxes" ! Compose filename

         ! Open file for writing
         OPEN(15, FILE=filename)

         WRITE(15,*) TIMESTEP, BOUNDARY_COLL_COUNT, WALL_COLL_COUNT, LINE_EMIT_COUNT
      END IF

   END SUBROUTINE DUMP_FLUXES_FILE


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! FUNCTION ERF1 -> Error function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
   FUNCTION ERF1(X) RESULT(out)

   IMPLICIT NONE

   REAL(KIND=8)             :: out, T
   REAL(KIND=8), INTENT(IN) :: X

   REAL(KIND=8), PARAMETER :: P = 0.3275911 ,     A1 = 0.254829592 ,    &
                              A2 = -0.284496736 , A3 = 1.421413741,     &
                              A4 = -1.453152027 , A5 = 1.061405429

   T = 1./(1.+P*ABS(X))
   out = 1 - (A1*T + A2*T**2 + A3*T**3 + A4*T**4 + A5*T**5) * EXP(-X**2)
   out = SIGN(out,X)

   RETURN

   END FUNCTION ERF1

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! FUNCTION FLX -> Function used to sample the velocity orthogonal !!!!
   ! to a surface !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION FLX(SN, TINF, M) RESULT(OUT)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: SN,TINF,M
      REAL(KIND=8) :: OUT
      REAL(KIND=8) :: R1,R2
      REAL(KIND=8) :: y,fM,BETA, KAPPA, ACCA

      BETA = 1./SQRT(2.*KB/M*TINF)
      ! IF (BOOL_KAPPA_DISTRIBUTION) BETA = 1./SQRT(2.*KB/M*TINF*(KAPPA_C-3./2.))

      ACCA = SQRT(SN**2+2.)                              ! Tmp variable
      KAPPA = 2./(SN+ACCA) * EXP(0.5 + 0.5*SN*(SN-ACCA)) ! variable

      ! Step 1.
      DO
         R1 = rf()
         y  = -3.+6.*R1

         ! Step 2.

         R2 = rf()
         fM = KAPPA*(y+sn)*EXP(-y**2)
         ! IF (BOOL_KAPPA_DISTRIBUTION) fM = KAPPA*(y+sn)*GAMMA(KAPPA_C)/GAMMA(KAPPA_C-1./2.)/(1+y**2)**KAPPA_C

         ! Step 3. 

         IF (R2 .LE. fM) THEN
            OUT = y/BETA
            EXIT
         END IF
      END DO

      RETURN

   END FUNCTION FLX


      SUBROUTINE THERMAL_BATH

      IMPLICIT NONE

      REAL(KIND=8)             :: VXP, VYP, VZP, EROT, EVIB, M
      INTEGER                  :: S_ID, JP

      DO JP = 1, NP_PROC
         S_ID = particles(JP)%S_ID
         M = SPECIES(S_ID)%MOLECULAR_MASS
         CALL MAXWELL(0.d0, 0.d0, 0.d0, &
                     TBATH, TBATH, TBATH, &
                     VXP, VYP, VZP, M)
         
         CALL INTERNAL_ENERGY(SPECIES(S_ID)%ROTDOF, TBATH, EROT)
         CALL INTERNAL_ENERGY(SPECIES(S_ID)%VIBDOF, TBATH, EVIB)
         particles(JP)%VX = VXP
         particles(JP)%VY = VYP
         particles(JP)%VZ = VZP
         particles(JP)%EROT = EROT
         particles(JP)%EVIB = EVIB
      END DO

   END SUBROUTINE THERMAL_BATH


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE STRIP_COMMENTS -> strips a part of string containing a comment, !!
   ! and removes white spaces as well.                                          !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE STRIP_COMMENTS(str,c)

      ! This subroutine checks if there is any character of type 'c' in the string str,
      ! and if so, it removes all the remaining of the string, tagged as a comment.
      ! Also, white spaces are removed by a call to the trim() function.

      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(INOUT) :: str
      CHARACTER(LEN=1),INTENT(IN)    :: c !comment character

      CHARACTER(LEN=LEN(str)) :: str_tmp
      INTEGER :: i
      
      ! Check if there is any comment to trim
      i = INDEX(str,c)
      IF (i .GT. 0) THEN
         str_tmp = str(1:i-1)
      ELSE
         str_tmp = str
      END IF
      
      ! Assign str, removing trailing blank spaces if any
      str = TRIM(str_tmp)

   END SUBROUTINE STRIP_COMMENTS


   SUBROUTINE SPLIT_STR(STRING, DELIMITER, STRARRAY, N_STR)
      ! splitstring splits a string to an array of
      ! substrings based on a selected delimiter
      ! note any facing space/blank in substrings will be removed

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE, INTENT(OUT) :: STRARRAY(:)
      CHARACTER(LEN=*), INTENT(IN) :: STRING
      INTEGER   :: n, i, j, idx
      CHARACTER(len=80) :: STRTMP = ''
      CHARACTER, INTENT(IN) :: DELIMITER

      n=LEN(STRING)
      ALLOCATE(STRARRAY(n))
      j = 1
      idx = 0
      DO i=1, n
         
         IF (STRING(i:i) /= DELIMITER) THEN
            STRTMP(j:j) = STRING(i:i)
            j = j + 1
            IF (i==n .OR. STRING(i+1:i+1) == DELIMITER) THEN
               j = 1
               idx = idx + 1
               STRARRAY(idx) = STRTMP
               STRTMP = ''
            END IF
         END IF
      END DO

      N_STR = idx

   END SUBROUTINE SPLIT_STR


   FUNCTION CROSS(A, B)
      REAL(KIND=8), DIMENSION(3) :: CROSS
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: A, B
   
      CROSS(1) = A(2) * B(3) - A(3) * B(2)
      CROSS(2) = A(3) * B(1) - A(1) * B(3)
      CROSS(3) = A(1) * B(2) - A(2) * B(1)

   END FUNCTION CROSS

   FUNCTION DOT(A, B)
      REAL(KIND=8) :: DOT
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: A, B
   
      DOT = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      
   END FUNCTION DOT



   FUNCTION RANDINT(X)
      REAL(KIND=8), INTENT(IN) :: X
      INTEGER :: RANDINT

      RANDINT = INT(X)
      IF (X-RANDINT > rf()) RANDINT = RANDINT + 1

   END FUNCTION RANDINT


   RECURSIVE SUBROUTINE QUICKSORT(ARRAY, ORDER, LO, HI)

      INTEGER, INTENT(IN) :: LO, HI
      REAL(KIND=8), INTENT(INOUT) :: ARRAY(:)
      INTEGER, INTENT(INOUT) :: ORDER(:)
      INTEGER :: I, J, P, TEMPORDER
      REAL(KIND=8) :: PIVOT, TEMPA

      IF (LO >= 1 .AND. HI >= 1 .AND. LO < HI) THEN
         PIVOT = ARRAY(HI)
         I = LO - 1
         DO J = LO, HI-1
            IF (ARRAY(J) <= PIVOT) THEN
               I = I + 1
               TEMPA = ARRAY(I); ARRAY(I) = ARRAY(J); ARRAY(J) = TEMPA
               TEMPORDER = ORDER(I); ORDER(I) = ORDER(J); ORDER(J) = TEMPORDER
            END IF
         END DO
         TEMPA = ARRAY(I+1); ARRAY(I+1) = ARRAY(HI); ARRAY(HI) = TEMPA
         TEMPORDER = ORDER(I+1); ORDER(I+1) = ORDER(HI); ORDER(HI) = TEMPORDER
         P = I + 1

         CALL QUICKSORT(ARRAY, ORDER, LO, P-1)
         CALL QUICKSORT(ARRAY, ORDER, P+1, HI)
      END IF
   END SUBROUTINE QUICKSORT


   SUBROUTINE TIMER_START(SECTION_ID)
      
      INTEGER, INTENT(IN) :: SECTION_ID

      CALL CPU_TIME(TIMERS_START_TIME(SECTION_ID))

   END SUBROUTINE


   SUBROUTINE TIMER_STOP(SECTION_ID)

      INTEGER, INTENT(IN) :: SECTION_ID
      REAL(KIND=8) :: TIMENOW
      
      CALL CPU_TIME(TIMENOW)
      TIMERS_ELAPSED(SECTION_ID) = TIMERS_ELAPSED(SECTION_ID) + (TIMENOW - TIMERS_START_TIME(SECTION_ID))

   END SUBROUTINE


   SUBROUTINE TIMER_SUMMARY

      REAL(KIND=8) :: TOTAL_ELAPSED

      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: ALL_TIMERS_ELAPSED
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: GLOBAL_TIMERS_ELAPSED
      REAL(KIND=8), DIMENSION(6) :: TIMERS_MAX_AMONG_PROC
      REAL(KIND=8), DIMENSION(6) :: TIMERS_MIN_AMONG_PROC
      REAL(KIND=8), DIMENSION(6) :: TIMERS_AVG_AMONG_PROC


      ALLOCATE(ALL_TIMERS_ELAPSED(6*N_MPI_THREADS))
      CALL MPI_GATHER(TIMERS_ELAPSED, 6, MPI_DOUBLE_PRECISION, & 
      ALL_TIMERS_ELAPSED, 6, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      GLOBAL_TIMERS_ELAPSED = RESHAPE(ALL_TIMERS_ELAPSED, [6, N_MPI_THREADS])

      IF (PROC_ID == 0) THEN
         TIMERS_MAX_AMONG_PROC = MAXVAL(GLOBAL_TIMERS_ELAPSED, DIM = 2)
         TIMERS_MIN_AMONG_PROC = MINVAL(GLOBAL_TIMERS_ELAPSED, DIM = 2)
         TIMERS_AVG_AMONG_PROC = SUM(GLOBAL_TIMERS_ELAPSED, DIM = 2)/DBLE(N_MPI_THREADS)
         TOTAL_ELAPSED = SUM(TIMERS_MAX_AMONG_PROC)

         WRITE(*,*) '==================================================================='
         WRITE(*,*) '=======================     TIMING INFO     ======================='
         WRITE(*,*) '==================================================================='
         WRITE(*,*) 'Section             |    MAX    |    MIN    |    AVG    |  MAX/TOT'

         WRITE(*,'(A21,F9.2,A3,F9.2,A3,F9.2,A3,F9.2,A2)') ' Misc:               ',  &
         TIMERS_MAX_AMONG_PROC(1),   ' s ', &
         TIMERS_MIN_AMONG_PROC(1),   ' s ', &
         TIMERS_AVG_AMONG_PROC(1),   ' s ', &
         100*TIMERS_MAX_AMONG_PROC(1)/TOTAL_ELAPSED, '%.'

         WRITE(*,'(A21,F9.2,A3,F9.2,A3,F9.2,A3,F9.2,A2)') ' Field solution:     ',  &
         TIMERS_MAX_AMONG_PROC(2),   ' s ', &
         TIMERS_MIN_AMONG_PROC(2),   ' s ', &
         TIMERS_AVG_AMONG_PROC(2),   ' s ', &
         100*TIMERS_MAX_AMONG_PROC(2)/TOTAL_ELAPSED, '%.'

         WRITE(*,'(A21,F9.2,A3,F9.2,A3,F9.2,A3,F9.2,A2)') ' Particle movement:  ',  &
         TIMERS_MAX_AMONG_PROC(3),   ' s ', &
         TIMERS_MIN_AMONG_PROC(3),   ' s ', &
         TIMERS_AVG_AMONG_PROC(3),   ' s ', &
         100*TIMERS_MAX_AMONG_PROC(3)/TOTAL_ELAPSED, '%.'

         WRITE(*,'(A21,F9.2,A3,F9.2,A3,F9.2,A3,F9.2,A2)') ' File output:        ',  &
         TIMERS_MAX_AMONG_PROC(4),   ' s ', &
         TIMERS_MIN_AMONG_PROC(4),   ' s ', &
         TIMERS_AVG_AMONG_PROC(4),   ' s ', &
         100*TIMERS_MAX_AMONG_PROC(4)/TOTAL_ELAPSED, '%.'

         WRITE(*,'(A21,F9.2,A3,F9.2,A3,F9.2,A3,F9.2,A2)') ' MPI particle comm.: ',  &
         TIMERS_MAX_AMONG_PROC(5),   ' s ', &
         TIMERS_MIN_AMONG_PROC(5),   ' s ', &
         TIMERS_AVG_AMONG_PROC(5),   ' s ', &
         100*TIMERS_MAX_AMONG_PROC(5)/TOTAL_ELAPSED, '%.'

         WRITE(*,'(A21,F9.2,A3,F9.2,A3,F9.2,A3,F9.2,A2)') ' Collisions:         ',  &
         TIMERS_MAX_AMONG_PROC(6),   ' s ', &
         TIMERS_MIN_AMONG_PROC(6),   ' s ', &
         TIMERS_AVG_AMONG_PROC(6),   ' s ', &
         100*TIMERS_MAX_AMONG_PROC(6)/TOTAL_ELAPSED, '%.'

         WRITE(*,*) '===================================================================='
      END IF

      DEALLOCATE(ALL_TIMERS_ELAPSED)

   END SUBROUTINE TIMER_SUMMARY



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



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE BINARY_SEARCH -> Finds the INDEX such that                           !
   ! ARRAY(INDEX) < VALUE < ARRAY(INDEX+1), where ARRAY is monotonically increasing. !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION BINARY_SEARCH(VALUE, ARRAY) RESULT(INDEX)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: VALUE
      INTEGER :: L, R
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ARRAY
      INTEGER :: INDEX

      L = LBOUND(ARRAY, DIM=1)
      R = UBOUND(ARRAY, DIM=1)

      INDEX = -1
      IF (VALUE .LT. ARRAY(L) .OR. VALUE .GT. ARRAY(R)) THEN
         WRITE(*,*) 'Particle out of bounds!'
         RETURN
      ELSE IF (R == L+1) THEN
         INDEX = L
         RETURN
      ELSE
         DO
            INDEX = (L+R)/2
            IF (ARRAY(INDEX) .LE. VALUE) THEN
               IF (ARRAY(INDEX+1) .GT. VALUE) RETURN
               L = INDEX
            ELSE
               IF (ARRAY(INDEX-1) .LE. VALUE) THEN
                  INDEX = INDEX-1
                  RETURN
               END IF
               R = INDEX
            END IF
         END DO
         RETURN
      END IF

   END FUNCTION BINARY_SEARCH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INTERP_CS -> Interpolates a tabulated cross-section given an energy  !
   ! There are a few variations. Here:                                               !
   ! * We interpolate linearly                                                       !
   ! * For energy values out of the table we set the cross section to zero           !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION INTERP_CS(VALUE_EN, TABLE_EN, TABLE_CS) RESULT(VALUE_CS)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: VALUE_EN
      INTEGER :: L, R
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: TABLE_EN, TABLE_CS
      INTEGER :: INDEX
      REAL(KIND=8) :: VALUE_CS

      L = LBOUND(TABLE_EN, DIM=1)
      R = UBOUND(TABLE_EN, DIM=1)

      INDEX = -1
      IF (VALUE_EN .LT. TABLE_EN(L)) THEN
         ! Lower than lower energy value
         !VALUE_CS = TABLE_CS(L)
         VALUE_CS = 0
         RETURN
      ELSE IF (VALUE_EN .GT. TABLE_EN(R)) THEN
         ! Higher than highest energy value
         !VALUE_CS = TABLE_CS(R)
         VALUE_CS = 0
         RETURN
      ELSE IF (R == L+1) THEN
         ! Only two values in the table
         VALUE_CS = TABLE_CS(L) + (TABLE_CS(R)-TABLE_CS(L))*(VALUE_EN-TABLE_EN(L))/(TABLE_EN(R)-TABLE_EN(L))
         RETURN
      ELSE
         DO
            INDEX = (L+R)/2
            IF (TABLE_EN(INDEX) .LE. VALUE_EN) THEN
               IF (TABLE_EN(INDEX+1) .GT. VALUE_EN) EXIT
               L = INDEX
            ELSE
               IF (TABLE_EN(INDEX-1) .LE. VALUE_EN) THEN
                  INDEX = INDEX-1
                  EXIT
               END IF
               R = INDEX
            END IF
         END DO
         ! The value we are looking for is between INDEX and INDEX+1.
         L = INDEX
         R = INDEX+1
         VALUE_CS = TABLE_CS(L) + (TABLE_CS(R)-TABLE_CS(L))*(VALUE_EN-TABLE_EN(L))/(TABLE_EN(R)-TABLE_EN(L))
         RETURN
      END IF

   END FUNCTION INTERP_CS


   SUBROUTINE SKIP_TO(UNIT, STR, STAT)

      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: STR
      INTEGER, INTENT(IN) :: UNIT
      INTEGER, INTENT(OUT) :: STAT
      CHARACTER :: CH
      INTEGER :: IO
    
      DO
         READ(UNIT, IOSTAT=IO) CH

         IF (IO/=0) THEN
            STAT = 1
            RETURN
         END IF
    
         IF (CH==STR(1:1)) THEN
            IF (LEN(STR) == 1) THEN
               STAT = 0
               RETURN
            END IF
            CALL CHECK(UNIT, STR(2:), STAT)
            IF (STAT == 0) RETURN
         END IF
    
      END DO
   END SUBROUTINE

    
   SUBROUTINE CHECK(UNIT, STR, STAT)
      CHARACTER(*), INTENT(IN) :: STR
      INTEGER, INTENT(IN) :: UNIT
      INTEGER, INTENT(OUT) :: STAT
      CHARACTER :: CH
      INTEGER :: I, IO

      STAT = 1
      I = 0

      DO
         I = I + 1

         READ(UNIT, IOSTAT=IO) CH

         IF (IO/=0) RETURN

         IF (CH/=STR(I:I)) RETURN

         IF (I==LEN(STR)) THEN
            STAT = 0
            RETURN
         END IF
      END DO
   END SUBROUTINE CHECK


   FUNCTION SOLVE_QUADRATIC(A, B, C)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: A, B, C
      REAL(KIND=8), DIMENSION(2) :: SOLVE_QUADRATIC
      REAL(KIND=8) :: DELTA, R

      DELTA = B*B-4.0*A*C

      IF (DELTA < 0) THEN
         SOLVE_QUADRATIC(1) = -1
         SOLVE_QUADRATIC(2) = -1
      ELSE IF (DELTA == 0) THEN
         SOLVE_QUADRATIC(1) = B/A
         SOLVE_QUADRATIC(2) = SOLVE_QUADRATIC(1)
      ELSE IF (DELTA > 0) THEN
         IF (A == 0) THEN
            SOLVE_QUADRATIC(1) = 0.5*C/B
            SOLVE_QUADRATIC(2) = SOLVE_QUADRATIC(1)
         ELSE
            R = -B - SIGN(SQRT(DELTA), B)
            SOLVE_QUADRATIC(1) = 2.0*C/R
            SOLVE_QUADRATIC(2) = 0.5*R/A
         END IF
      END IF
   END FUNCTION SOLVE_QUADRATIC



   FUNCTION DQDCRT(A, B, C)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: A, B, C     !! COEFFICIENTS
      REAL(KIND=8), DIMENSION(2) :: DQDCRT   !! REAL COMPONENTS OF ROOTS

      REAL(KIND=8) :: D, R, W

      IF (A == 0) THEN
         !IT IS REALLY A LINEAR EQUATION:
         IF (B == 0) THEN  !DEGENERATE CASE, JUST RETURN ZEROS
            DQDCRT = 0
         ELSE
         !THERE IS ONLY ONE ROOT (REAL), SO JUST DUPLICATE IT:
            DQDCRT = -C/B
         END IF
      ELSE
         IF (C == 0) THEN
            DQDCRT(1) = 0
            DQDCRT(2) = -B/A
         ELSE
            D = B*B - 4.0*A*C
               IF (D < 0) THEN
               !COMPLEX ROOTS
                  DQDCRT = -1
               ELSE
               !DISTINCT REAL ROOTS
               R = SQRT(D)
               IF (B /= 0) THEN
                  W = -(B + SIGN(R, B))
                  DQDCRT(1) = 2.0*C/W
                  DQDCRT(2) = 0.5*W/A
               ELSE
                  DQDCRT(1) = ABS(0.5*R/A)
                  DQDCRT(2) = -DQDCRT(1)
               END IF
            END IF
         END IF
      END IF

   END FUNCTION DQDCRT


END MODULE tools
