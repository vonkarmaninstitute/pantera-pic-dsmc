MODULE tools

USE grid_and_partition
USE mpi_common
USE global

CONTAINS 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION rf -> Pseudo-random number generator (RNG) !!!!!!!!!!!!!!!!
   ! shorthand for fortran builtin RANDOM_NUMBER function !!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   FUNCTION rf() RESULT(out)

      IMPLICIT NONE

      REAL(KIND=8) :: out

      out = RAND()

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

  SUBROUTINE MAXWELL(UX, UY, UZ, TX, TY, TZ, TR, VX, VY, VZ, EI, M)

  IMPLICIT NONE

  REAL(KIND=8), INTENT(IN)    :: UX, UY, UZ, TX, TY, TZ, TR
  !REAL(KIND=8), INTENT(IN)    :: RGAS !RGAS = KB / M
  REAL(KIND=8), INTENT(IN)    :: M ! Molecular mass
  REAL(KIND=8), INTENT(INOUT) :: VX, VY, VZ, EI

  INTEGER                     :: I
  REAL(KIND=8)                :: PI, PI2, KB
  REAL(KIND=8)                :: R, R1, RO, TETA, BETA
  REAL(KIND=8), DIMENSION(3)  :: VEL, TT

  PI   = 3.141593
  PI2  = 2.*PI

  KB = 1.38064852E-23

  TT(1) = TX
  TT(2) = TY
  TT(3) = TZ

  DO I = 1,3

     ! Step 1.
     R = RAND()
      
     TETA = PI2*R

     ! Step 2.
     
     BETA = 1./SQRT(2.*KB/M*TT(I))

     ! R goes from 0 to 1 included. Remove the extremes in a very rough way
     ! or the log() will explode
     !R = rf()*0.99999999999 + 1.0d-13 ! Please fix this!!!!! MAKE THIS ELEGANT!!!
     ! Fixed at the root.
     R1 = RAND()
      DO WHILE (R1 < 1.0D-13)
         R1 = RAND()
      END DO

     RO = SQRT(-LOG(R1))/BETA ! The random number here shouldn't be correlated to the one for teta!!

     VEL(I) = RO*SIN(TETA)

  END DO

  ! Step 3.

  VX = UX + VEL(1)
  VY = UY + VEL(2)
  VZ = UZ + VEL(3)

  ! Step 4. 

  R = RAND()
  DO WHILE (R < 1.0D-13)
    R = RAND()
  END DO

  EI = -LOG(R)*KB*TR

  RETURN

  END SUBROUTINE MAXWELL

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
     INTEGER :: IP

     WRITE(filename, "(A13,I0.5,A6,I0.8)") "./dumps/proc_", PROC_ID, "_time_", TIMESTEP ! Compose filename

     ! Open file for writing
     OPEN(10, FILE=filename )

     WRITE(10,*) '% Timestep | X | Y | Z | VX | VY | VZ | S_ID | proc_ID'
     DO IP = 1, NP_PROC
        WRITE(10,*) TIMESTEP, particles(IP)%X, particles(IP)%Y, particles(IP)%Z, &
        particles(IP)%VX, particles(IP)%VY, particles(IP)%VZ, particles(IP)%S_ID, PROC_ID
     END DO

  END SUBROUTINE DUMP_PARTICLES_FILE


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
  REAL(KIND=8) :: y,fM,BETA,KB, KAPPA, ACCA

  KB = 1.38064852E-23
  BETA = 1./SQRT(2.*KB/M*TINF)

  ACCA = SQRT(SN**2+2.)                              ! Tmp variable
  KAPPA = 2./(SN+ACCA) * EXP(0.5 + 0.5*SN*(SN-ACCA)) ! variable

   ! Step 1.

1  R1 = rf()
   y  = -3.+6.*R1

  ! Step 2.

   R2 = rf()
   fM = KAPPA*(y+sn)*EXP(-y**2)

  ! Step 3. 

   IF (R2 .LE. fM) THEN
      out = y/BETA
   ELSE
      GO TO 1
   END IF

  RETURN

  END FUNCTION FLX

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
  else
     str_tmp = str
  end if
 
  ! Assign str, removing trailing blank spaces if any
  str = TRIM(str_tmp)

  END SUBROUTINE STRIP_COMMENTS


  SUBROUTINE SPLIT_STR(STRING, DELIMITER, STRARRAY, N_STR)
   ! splitstring splits a string to an array of
   ! substrings based on a selected delimiter
   ! note any facing space/blank in substrings will be removed

   IMPLICIT NONE

   INTEGER :: N_STR
   CHARACTER(LEN=80), allocatable :: STRARRAY(:)
   CHARACTER(LEN=*) :: STRING
   INTEGER   :: n, i, j, idx
   CHARACTER(len=80) :: STRTMP = ''
   CHARACTER :: DELIMITER

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

END MODULE tools