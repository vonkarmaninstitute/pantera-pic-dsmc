MODULE tools

USE grid_and_partition
USE mpi_common
USE mt19937_64
USE global

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

  PI   = 2.0D0*ASIN(1.0D0)
  PI2  = 2.0D0*PI

  KB = 1.38064852E-23

  TT(1) = TX
  TT(2) = TY
  TT(3) = TZ

  DO I = 1,3

     ! Step 1.
     R = rf()
      
     TETA = PI2*R

     ! Step 2.
     
     BETA = 1./SQRT(2.*KB/M*TT(I))

     ! R goes from 0 to 1 included. Remove the extremes in a very rough way
     ! or the log() will explode
     !R = rf()*0.99999999999 + 1.0d-13 ! Please fix this!!!!! MAKE THIS ELEGANT!!!
     ! Fixed at the root.
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

  ! Step 4. 

  R = rf()
  DO WHILE (R < 1.0D-13)
    R = rf()
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

     !! old !! WRITE(filename, "(A13,I0.5,A6,I0.8)") "./dumps/proc_", PROC_ID, "_time_", TIMESTEP ! Hard-coded filename
     WRITE(filename, "(A,A6,I0.5,A6,I0.8)") ADJUSTL(TRIM(DUMP_PART_PATH)), &
                                      "/proc_", PROC_ID, "_time_", TIMESTEP ! Compose filename

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE DUMP_GLOBAL_MOMENTS_FILE                                 !!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DUMP_GLOBAL_MOMENTS_FILE(TIMESTEP)

     ! This subroutine computes and dumps global moments (computed from all particles in the domain) 
     ! to a file, appending a line for every timestep.
     ! The moments to be computed are 33, corresponding to the 14 moments equations and their closing
     ! moments. This should be enough for now. For this reason, LOCAL_MOMENTS and MOMENTS arrays have
     ! length 33.
     ! 
     ! TIMESTEP: current timestep of the simulation

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: TIMESTEP

     REAL(KIND=8), DIMENSION(33) :: LOCAL_MOMENTS, MOMENTS ! See description of the function
     INTEGER :: i

     MOMENTS = 0
     LOCAL_MOMENTS = 0

     ! Compute moments at current time. Every process does it.
     ! Moments are all extensive: rho, rho u_i, Pij, etc, so they can be then
     ! added by a reduce command, to get the whole moments.
     CALL COMPUTE_GLOBAL_MOMENTS(LOCAL_MOMENTS)

     ! Now, processors should communicate with each others. All moments are computed as 
     ! extensive quantities, so I can add them. Then I will extract the velocities (intensive).
     !
     ! Note that I used the total volume for the computation, that's why I have to add 
     ! the moments. It's like considering all particles.
     CALL MPI_REDUCE(LOCAL_MOMENTS,  MOMENTS, 33, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

     ! Master writes moments to file
     IF (PROC_ID .EQ. 0) THEN 

        ! First, redefine moments 2, 3 and 4 as velocities, by dividing by the density
        MOMENTS(2) = MOMENTS(2)/(MOMENTS(1) + 1.0e-35)
        MOMENTS(3) = MOMENTS(3)/(MOMENTS(1) + 1.0e-35)
        MOMENTS(4) = MOMENTS(4)/(MOMENTS(1) + 1.0e-35)

        ! Then, write
        OPEN(12, FILE=DUMP_GLOB_MOM_FILENAME, STATUS="old", POSITION="append", ACTION="write") 
        WRITE(12, *) TIMESTEP*DT, (MOMENTS(i), i = 1, SIZE(MOMENTS))
        CLOSE(12)
     END IF

  END SUBROUTINE DUMP_GLOBAL_MOMENTS_FILE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE COMPUTE_GLOBAL_MOMENTS => Computes moments, from all particles in the domain !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_GLOBAL_MOMENTS(MOM_VECT)

     REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: MOM_VECT
     INTEGER :: JP
     REAL(KIND=8) :: RHO, Ux, Uy, Uz, Pxx, Pxy, Pxz, Pyy, Pyz, Pzz
     REAL(KIND=8) :: qx, qy, qz, Qxxx, Qxxy, Qxyy, Qyyy, Qyyz, Qyzz, Qzzz, Qxxz, Qxzz, Qxyz
     REAL(KIND=8) :: Riijj, Rxxjj, Rxyjj, Rxzjj, Ryyjj, Ryzjj, Rzzjj, Sxiijj, Syiijj, Sziijj
     REAL(KIND=8) :: MASS, CX_NOW, CY_NOW, CZ_NOW, C2_NOW, VOL

     ! Check size of MOM_VECT
     IF (SIZE(MOM_VECT) .NE. 33) THEN
        WRITE(*,*) "ERROR! Only 33 moments can be computed for now:"
        WRITE(*,*) "Rho, rho*ui, Pij, qi, Qijk, Riijj, Rijkk, Siijjx (for i,j,k = x,y,z)."
        WRITE(*,*) "Check how this function is called!"
        WRITE(*,*) "In: COMPUTE_GLOBAL_MOMENTS."
        STOP
     END IF

     ! Every processor should perform the operation on its own particles

     VOL = (XMAX - XMIN)*(YMAX - YMIN)*(ZMAX - ZMIN) ! Domain volume

     ! Init moments to zero (I'm sure there is a smarter way.....)
     RHO  = 0

     Ux  = 0; Uy  = 0; Uz  = 0
     Pxx = 0; Pxy = 0; Pxz = 0; Pyy = 0; Pyz = 0; Pzz = 0
     qx  = 0; qy  = 0; qz  = 0

     Qxxx = 0; Qxxy = 0; Qxyy = 0; Qyyy = 0; Qyyz = 0; 
     Qyzz = 0; Qzzz = 0; Qxxz = 0; Qxzz = 0; Qxyz = 0; 

     Riijj = 0
     Rxxjj = 0; Rxyjj = 0; Rxzjj = 0; 
     Ryyjj = 0; Ryzjj = 0; Rzzjj = 0;

     Sxiijj = 0; Syiijj = 0; Sziijj = 0

     ! Compute average velocities and density
     DO JP = 1, NP_PROC ! Loop on the particles that the processor owns

        MASS = SPECIES(particles(JP)%S_ID)%MOLMASS ! [kg]

        RHO = RHO + MASS/VOL*Fnum

        Ux = Ux + particles(JP)%VX/NP_PROC
        Uy = Uy + particles(JP)%VY/NP_PROC
        Uz = Uz + particles(JP)%VZ/NP_PROC

     END DO

     ! Now proceed with central moments
     DO JP = 1, NP_PROC

        MASS = SPECIES(particles(JP)%S_ID)%MOLMASS ! [kg]

        CX_NOW = particles(JP)%VX - Ux
        CY_NOW = particles(JP)%VY - Uy
        CZ_NOW = particles(JP)%VZ - Uz

        C2_NOW = CX_NOW**2 + CY_NOW**2 + CZ_NOW**2

        ! Compute pressure terms
        Pxx = Pxx + MASS*CX_NOW*CX_NOW*Fnum/VOL
        Pxy = Pxy + MASS*CX_NOW*CY_NOW*Fnum/VOL
        Pxz = Pxz + MASS*CX_NOW*CZ_NOW*Fnum/VOL
        Pyy = Pyy + MASS*CY_NOW*CY_NOW*Fnum/VOL
        Pyz = Pyz + MASS*CY_NOW*CZ_NOW*Fnum/VOL
        Pzz = Pzz + MASS*CZ_NOW*CZ_NOW*Fnum/VOL

        ! Heat flux vector components
        qx = qx + MASS*CX_NOW*C2_NOW*Fnum/VOL
        qy = qy + MASS*CY_NOW*C2_NOW*Fnum/VOL
        qz = qz + MASS*CZ_NOW*C2_NOW*Fnum/VOL
 
        ! Heat flux tensor components
        Qxxx = Qxxx + MASS*CX_NOW*CX_NOW*CX_NOW*Fnum/VOL
        Qxxy = Qxxy + MASS*CX_NOW*CX_NOW*CY_NOW*Fnum/VOL
        Qxyy = Qxyy + MASS*CX_NOW*CY_NOW*CY_NOW*Fnum/VOL
        Qyyy = Qyyy + MASS*CY_NOW*CY_NOW*CY_NOW*Fnum/VOL
        Qyyz = Qyyz + MASS*CY_NOW*CY_NOW*CZ_NOW*Fnum/VOL
        Qyzz = Qyzz + MASS*CY_NOW*CZ_NOW*CZ_NOW*Fnum/VOL
        Qzzz = Qzzz + MASS*CZ_NOW*CZ_NOW*CZ_NOW*Fnum/VOL
        Qxxz = Qxxz + MASS*CX_NOW*CX_NOW*CZ_NOW*Fnum/VOL
        Qxzz = Qxzz + MASS*CX_NOW*CZ_NOW*CZ_NOW*Fnum/VOL
        Qxyz = Qxyz + MASS*CX_NOW*CY_NOW*CZ_NOW*Fnum/VOL

        ! Riijj
        Riijj = Riijj + MASS*C2_NOW*C2_NOW*Fnum/VOL

        ! Rijkk
        Rxxjj = Rxxjj + MASS*CX_NOW*CX_NOW*C2_NOW*Fnum/VOL
        Rxyjj = Rxyjj + MASS*CX_NOW*CY_NOW*C2_NOW*Fnum/VOL
        Rxzjj = Rxzjj + MASS*CX_NOW*CZ_NOW*C2_NOW*Fnum/VOL
        Ryyjj = Ryyjj + MASS*CY_NOW*CY_NOW*C2_NOW*Fnum/VOL
        Ryzjj = Ryzjj + MASS*CY_NOW*CZ_NOW*C2_NOW*Fnum/VOL
        Rzzjj = Rzzjj + MASS*CZ_NOW*CZ_NOW*C2_NOW*Fnum/VOL

        ! Sxiijj
        Sxiijj = Sxiijj + MASS*CX_NOW*C2_NOW*C2_NOW*Fnum/VOL
        Syiijj = Syiijj + MASS*CY_NOW*C2_NOW*C2_NOW*Fnum/VOL
        Sziijj = Sziijj + MASS*CZ_NOW*C2_NOW*C2_NOW*Fnum/VOL

     END DO

     ! Assemble moments vector
     MOM_VECT(1)  = RHO
     MOM_VECT(2)  = RHO*Ux ! NOTE HERE, rho Ux, not Ux
     MOM_VECT(3)  = RHO*Uy ! NOTE HERE, rho Uy, not Uy
     MOM_VECT(4)  = RHO*Uz ! NOTE HERE, rho Uz, not Uz
     MOM_VECT(5)  = Pxx
     MOM_VECT(6)  = Pxy
     MOM_VECT(7)  = Pxz
     MOM_VECT(8)  = Pyy
     MOM_VECT(9)  = Pyz
     MOM_VECT(10) = Pzz
     MOM_VECT(11) = qx
     MOM_VECT(12) = qy
     MOM_VECT(13) = qz
     MOM_VECT(13) = qz
     MOM_VECT(14) = Qxxx
     MOM_VECT(15) = Qxxy
     MOM_VECT(16) = Qxyy
     MOM_VECT(17) = Qyyy
     MOM_VECT(18) = Qyyz
     MOM_VECT(19) = Qyzz
     MOM_VECT(20) = Qzzz
     MOM_VECT(21) = Qxxz
     MOM_VECT(22) = Qxzz
     MOM_VECT(23) = Qxyz
     MOM_VECT(24) = Riijj
     MOM_VECT(25) = Rxxjj
     MOM_VECT(26) = Rxyjj
     MOM_VECT(27) = Rxzjj
     MOM_VECT(28) = Ryyjj
     MOM_VECT(29) = Ryzjj
     MOM_VECT(30) = Rzzjj
     MOM_VECT(31) = Sxiijj
     MOM_VECT(32) = Syiijj
     MOM_VECT(33) = Sziijj

  END SUBROUTINE COMPUTE_GLOBAL_MOMENTS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FUNCTION INTERP_VECTOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  FUNCTION INTERP_VECTOR(x_vect, y_vect, x_target)

  ! Very basic interpolation routine. Very inefficient. Use bisection at least!

     REAL(KIND=8) :: INTERP_VECTOR, x_target
     REAL(KIND=8), DIMENSION(:) :: x_vect, y_vect

     INTEGER :: ii
     REAL(KIND=8) :: X1, X2, Y1, Y2

     ! Check that requested element is in same position
     IF (x_target .LE. x_vect(1)) THEN

        WRITE(*,*) 'Attention! In function INTERP_VECTOR(): out of vector range! Using first element.'
        INTERP_VECTOR = y_vect(1)
        RETURN

     ELSE IF (x_target .GE. x_vect(size(x_vect))) THEN

        WRITE(*,*) 'Attention! In function INTERP_VECTOR(): out vector range! Using last element.'
        WRITE(*,*) 'Asked value at: ', x_target, ' Max tabulated value: ', x_vect(size(x_vect))
        INTERP_VECTOR = y_vect(SIZE(y_vect))
        RETURN
       
     END IF

     ! Check that x_vect and y_vect have same size
     IF (SIZE(x_vect) .NE. SIZE(y_vect)) THEN
        WRITE(*,*) 'ERROR! In function INTERP_VECTOR(): number of elements of X and Y does not match. ABORTING!'
        WRITE(*,*) 'Asked value at: ', x_target, ' Min tabulated value: ', x_vect(1)
        STOP
     END IF

     ! Classical interpolation
     DO ii = 1, SIZE(x_vect)
        IF (x_vect(ii) .GT. x_target) THEN ! We just paseed it!

           X1 = x_vect(ii - 1)
           X2 = x_vect(ii)
           Y1 = y_vect(ii - 1)
           Y2 = y_vect(ii)
           
           INTERP_VECTOR = Y1 + (Y2 - Y1)/(X2 - X1)*(x_target - X1)
           RETURN

        END IF
     END DO

     INTERP_VECTOR = -1.0d0 ! Give it a wrong value, if I reach this point, and quit.
     WRITE(*,*) 'ERROR! Apparently the interpolation in function INTERP_VECTOR() did not go through! ABORTING!' 
     STOP

     RETURN
  END FUNCTION

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FUNCTION CROSS(a,b) -> computes cross product !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION CROSS(a,b)

     REAL(KIND=8), DIMENSION(3) :: CROSS
     REAL(KIND=8), DIMENSION(3), INTENT(IN) :: a, b

     CROSS(1) = a(2)*b(3) - a(3)*b(2)
     CROSS(2) = a(3)*b(1) - a(1)*b(3)
     CROSS(3) = a(1)*b(2) - a(2)*b(1)

  END FUNCTION CROSS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FUNCTION DOT(a,b) -> computes dot product !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION DOT(a,b)

     REAL(KIND=8) :: DOT
     REAL(KIND=8), DIMENSION(3), INTENT(IN) :: a, b

     DOT = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

  END FUNCTION DOT


END MODULE tools
