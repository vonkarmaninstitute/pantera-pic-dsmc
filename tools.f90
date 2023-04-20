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

      END DO

      ! Step 3.

      VX = UX + VEL(1)
      VY = UY + VEL(2)
      VZ = UZ + VEL(3)

      RETURN

   END SUBROUTINE MAXWELL



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

      USE global
      USE mpi_common

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

      USE global
      USE mpi_common

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename
      INTEGER :: IP

      WRITE(filename, "(A,A,I0.5,A6,I0.8)") TRIM(ADJUSTL(PARTDUMP_SAVE_PATH)), "proc_", PROC_ID, "_time_", TIMESTEP ! Compose filename

      ! Open file for writing
      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(10, FILE=filename, ACCESS='SEQUENTIAL', FORM='UNFORMATTED', STATUS='NEW', CONVERT='BIG_ENDIAN', RECL=80)
         DO IP = 1, NP_PROC
            WRITE(10) particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
            particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ, particles(IP)%EROT, particles(IP)%EVIB, &
            particles(IP)%S_ID, particles(IP)%IC, particles(IP)%DTRIM
         END DO
         CLOSE(10)
      ELSE
         OPEN(10, FILE=filename )
         !WRITE(10,*) '% X | Y | Z | VX | VY | VZ | EROT | EVIB | S_ID | IC | DTRIM'
         DO IP = 1, NP_PROC
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

      USE global
      USE mpi_common

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename
      INTEGER :: IP


      ! Dump particles that hit a boundary to file
      WRITE(filename, "(A,A,I0.5,A6,I0.8)") TRIM(ADJUSTL(PARTDUMP_SAVE_PATH)), "bound_proc_", PROC_ID ! Compose filename

      ! Open file for writing
      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(12, FILE=filename, ACCESS='SEQUENTIAL', FORM='UNFORMATTED', STATUS='OLD', CONVERT='BIG_ENDIAN', RECL=84)
         DO IP = 1, NP_DUMP_PROC
            WRITE(12) TIMESTEP, part_dump(IP)%X, part_dump(IP)%Y, part_dump(IP)%Z, &
            part_dump(IP)%VX, part_dump(IP)%VY, part_dump(IP)%VZ, part_dump(IP)%EROT, part_dump(IP)%EVIB, &
            part_dump(IP)%S_ID, part_dump(IP)%IC, part_dump(IP)%DTRIM
         END DO
         CLOSE(12)
      ELSE
         OPEN(12, FILE=filename )
         !WRITE(10,*) '% X | Y | Z | VX | VY | VZ | EROT | EVIB | S_ID | IPG | DTRIM'
         DO IP = 1, NP_DUMP_PROC
            WRITE(12,*) TIMESTEP, part_dump(IP)%X, part_dump(IP)%Y, part_dump(IP)%Z, &
            part_dump(IP)%VX, part_dump(IP)%VY, part_dump(IP)%VZ, part_dump(IP)%EROT, part_dump(IP)%EVIB, &
            part_dump(IP)%S_ID, part_dump(IP)%IC, part_dump(IP)%DTRIM
         END DO
         CLOSE(12)
      END IF

   END SUBROUTINE DUMP_BOUNDARY_PARTICLES_FILE


   SUBROUTINE READ_PARTICLES_FILE(TIMESTEP)

      USE global
      USE mpi_common

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

      USE global
      USE mpi_common

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TIMESTEP
      CHARACTER(LEN=512)  :: filename
      INTEGER :: IP

      DO IP = 1, NP_PROC
         IF (particles(IP)%DUMP_TRAJ) THEN
            !WRITE(*,*) 'Writing trajectory file for particle with ID ', particles(IP)%ID
            WRITE(filename, "(A,A,I0.15)") TRIM(ADJUSTL(TRAJDUMP_SAVE_PATH)), "trajectory_", particles(IP)%ID ! Compose filename
            ! Open file for writing
            OPEN(1610, FILE=filename, POSITION='APPEND')
            WRITE(1610,*) TIMESTEP, particles(IP)%S_ID, particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
            particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ
            CLOSE(1610)
         END IF
      END DO

   END SUBROUTINE DUMP_TRAJECTORY_FILE

   SUBROUTINE DUMP_FLUXES_FILE(TIMESTEP)

      ! This subroutine dumps particles at the screen. DUMPS MAY OVERWRITE, SINCE STDOUT IS BUFFERED!!!!
      ! Doesn't matter if I loop on processes!!!!
      ! 
      ! OLD OLD OL A loop is performed on the processes, making sure that
      ! OLD OLD OL processes don't do it simultaneously, messing up the dump.

      USE global
      USE mpi_common

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

      ACCA = SQRT(SN**2+2.)                              ! Tmp variable
      KAPPA = 2./(SN+ACCA) * EXP(0.5 + 0.5*SN*(SN-ACCA)) ! variable

      ! Step 1.
      DO
         R1 = rf()
         y  = -3.+6.*R1

         ! Step 2.

         R2 = rf()
         fM = KAPPA*(y+sn)*EXP(-y**2)

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


   FUNCTION RANDINT(X)
      REAL(KIND=8), INTENT(IN) :: X
      INTEGER :: RANDINT

      RANDINT = INT(X)
      IF (X-RANDINT > rf()) RANDINT = RANDINT + 1

   END FUNCTION RANDINT


   RECURSIVE SUBROUTINE QUICKSORT(A,ORDER,NA)

      ! DUMMY ARGUMENTS
      INTEGER, INTENT(IN) :: NA
      INTEGER, DIMENSION(NA), INTENT(IN OUT) :: ORDER
      REAL(KIND=8), DIMENSION(NA), INTENT(IN OUT) :: A

      ! LOCAL VARIABLES
      INTEGER :: LEFT, RIGHT
      !REAL(KIND=8) :: RANDOM
      REAL(KIND=8) :: PIVOT
      INTEGER :: TEMPORDER
      REAL(KIND=8) :: TEMPA
      INTEGER :: MARKER

      IF (NA > 1) THEN

         !RANDOM = rf()
         !PIVOT = A(INT(RANDOM*REAL(NA-1))+1)   ! Choice a random pivot (not best performance, but avoids worst-case)
         PIVOT = A(INT(DBLE(NA-1)/2.)+1)   ! Choice a random pivot (not best performance, but avoids worst-case)
         LEFT = 1
         RIGHT = NA
         ! Partition loop
         DO
            IF (LEFT >= RIGHT) EXIT
            DO
               IF (A(RIGHT) <= PIVOT) EXIT
               RIGHT = RIGHT - 1
            END DO
            DO
               IF (A(LEFT) >= PIVOT) EXIT
               LEFT = LEFT + 1
            END DO
            IF (LEFT < RIGHT) THEN
               TEMPA = A(LEFT); A(LEFT) = A(RIGHT); A(RIGHT) = TEMPA
               TEMPORDER = ORDER(LEFT); ORDER(LEFT) = ORDER(RIGHT); ORDER(RIGHT) = TEMPORDER
            END IF
         END DO

         IF (LEFT == RIGHT) THEN
            MARKER = LEFT + 1
         ELSE
            MARKER = LEFT
         END IF

         CALL QUICKSORT(A(:MARKER-1),ORDER(:MARKER-1),MARKER-1)
         CALL QUICKSORT(A(MARKER:),ORDER(:MARKER-1),NA-MARKER+1)

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

      
      IF (PROC_ID == 0) THEN
         TOTAL_ELAPSED = SUM(TIMERS_ELAPSED)

         WRITE(*,*) '=========================================='
         WRITE(*,*) '==========     TIMING INFO     ==========='
         WRITE(*,*) '=========================================='

         WRITE(*,'(A20,F6.1,A6,F4.1,A2)') 'Initialization:     ',  TIMERS_ELAPSED(1),   ' s  = ', &
                                          100*TIMERS_ELAPSED(1)/TOTAL_ELAPSED, '%.'
         WRITE(*,'(A20,F6.1,A6,F4.1,A2)') 'Field solution:     ',  TIMERS_ELAPSED(2),   ' s  = ', &
                                          100*TIMERS_ELAPSED(2)/TOTAL_ELAPSED, '%.'
         WRITE(*,'(A20,F6.1,A6,F4.1,A2)') 'Particle movement:  ',  TIMERS_ELAPSED(3), ' s  = ', &
                                          100*TIMERS_ELAPSED(3)/TOTAL_ELAPSED, '%.'
         WRITE(*,'(A20,F6.1,A6,F4.1,A2)') 'File output:        ',  TIMERS_ELAPSED(4),   ' s  = ', &
                                          100*TIMERS_ELAPSED(4)/TOTAL_ELAPSED, '%.'
         WRITE(*,'(A20,F6.1,A6,F4.1,A2)') 'MPI particle comm.: ',  TIMERS_ELAPSED(5),   ' s  = ', &
                                          100*TIMERS_ELAPSED(5)/TOTAL_ELAPSED, '%.'

         WRITE(*,*) '=========================================='
      END IF

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





END MODULE tools