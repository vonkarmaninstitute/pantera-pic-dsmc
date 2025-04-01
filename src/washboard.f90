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

MODULE washboard

   USE global
   USE tools

   CONTAINS
 
   SUBROUTINE WB_SCATTER(IPG, MASS, VX_I, VY_I, VZ_I)

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: IPG
      REAL(KIND=8), INTENT(INOUT) :: VX_I, VZ_I
      REAL(KIND=8), INTENT(IN) :: MASS

      REAL(KIND=8), INTENT(OUT) :: VY_I
      REAL(KIND=8) :: A, G, ALPHA_N, ALPHA_T, T, CI, LIMIT_TRAP, MAXX, MAXX_NEG, PHI_A, R, THETA_A, THETA_LIM, W
      INTEGER :: col_to_trap = 20
      REAL(KIND=8), DIMENSION(3) :: VELOCITIES
      REAL(KIND=8), DIMENSION(2) :: ANGLE_VALUES
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: theta_i_v, theta_r_v, A_v, g_v1, g_v2
      REAL(KIND=8), DIMENSION(40) :: g_v
      INTEGER :: COL, DOWN, I, J, TRAPPING, NDOWN, NUP, TD

      REAL(KIND=8) :: ITPL
      

      TRAPPING = 1
      VY_I = 0.d0
      NDOWN = 0
      NUP = 0
      TD = 0

      T = GRID_BC(IPG)%WALL_TEMP
      A = GRID_BC(IPG)%A
      G = GRID_BC(IPG)%B
      W = GRID_BC(IPG)%W
      ALPHA_N = GRID_BC(IPG)%ACC_N
      ALPHA_T = GRID_BC(IPG)%ACC_T

      !WRITE(*,*) 'A = ', A, ' B = ', G, ' W = ', W, ' T = ', T, ' ALPHA_N = ', ALPHA_N, ' ALPHA_T = ', ALPHA_T
      
      
      CI = SQRT(2*KB*T/MASS)
      LIMIT_TRAP = SQRT(2*W/MASS)

      !Minimum angle for which we check for a collision going up. If the reflected angle is lower than this 
      ! we assume the particle does not collide.
      THETA_LIM = 35./180.*PI


      !Get the grid where the pre-computed values were taken
      THETA_I_V = LINSPACE(PI/2, PI, 45)    ! DBDBDBBDBDBDBDDBBDDB include in the structure
      THETA_R_V = LINSPACE(0.0D0, PI/2, 45)
      !Careful!!! 
      !This is the range of validity of the precomputed values for 
      !the maximums of the probability distribution of the local surface normals.
      ! If you take values for A and g out of this range the sampling procedure is not guaranteed to be correct!
      A_V = LINSPACE(0.5d0, 1.5d0, 40)
      G_V1 = LINSPACE(0.2d0, 1.0d0, 20)       
      G_V2 = LINSPACE(1.05d0, 5.0d0, 20)    
  
      ! Concatenate g_v1 and g_v2 into g_v
      J = 1
      DO I = 1, SIZE(G_V1)
         G_V(J) = G_V1(I)
         J = J + 1
      END DO
      DO I = 1, SIZE(G_V2)
         G_V(J) = G_V2(I)
         J = J + 1
      END DO


      VZ_I = -SQRT(VZ_I*VZ_I + 2*W/MASS)

      THETA_A = ACOS(VZ_I/SQRT(VZ_I*VZ_I + VX_I*VX_I))
      PHI_A = 0



      DO WHILE (VZ_I < 0)
         DOWN = 1
         NDOWN = NDOWN + 1
         !print*, down, vx_i,vy_i,vz_i
         ! Sample local normal
         THETA_A = ACOS(VZ_I/SQRT(VX_I*VX_I + VY_I*VY_I + VZ_I*VZ_I))
         PHI_A = ATAN2(VY_I, VX_I)
         
         MAXX = INTERPOLATE(THETA_I_V, A_V, G_V, GRID_BC(IPG)%MAX_P_DN, THETA_A, A, G)/COS(THETA_A)
         ANGLE_VALUES = ACCEPT_REJECT_SURF(MAXX, THETA_A, PHI_A, A, G,DOWN)
         
         VELOCITIES(:) = CL_KERNEL_LOCAL(VX_I, VY_I, VZ_I, ALPHA_N, ALPHA_T, CI,ANGLE_VALUES(1),ANGLE_VALUES(2))
         VX_I =  VELOCITIES(1) 
         VY_I =  VELOCITIES(2)
         VZ_I =  VELOCITIES(3)

         THETA_A = ACOS(VZ_I/SQRT(VX_I*VX_I + VY_I*VY_I + VZ_I*VZ_I))  
         PHI_A = ATAN2(VY_I, VX_I)
         IF (VZ_I > 0) THEN
            R = rf()
            ITPL = INTERPOLATE(THETA_R_V, A_V, G_V, GRID_BC(IPG)%P_COLL_UP, THETA_A, A, G)
            !WRITE(*,*) 'ITPL = ', ITPL
            IF (R <= ITPL  ) THEN 
               COL = 1
               !WRITE(*,*) 'Upwards collision 1!'
            ELSE 
               COL = 0
               IF (VZ_I < LIMIT_TRAP .AND. TRAPPING == 1) THEN 
                  VZ_I = -VZ_I
               END IF 
            END IF 
            ! col = 0 -> no collision | col = 1 -> collision
            ! Upwards collision loop
            DO WHILE (COL == 1 .AND. THETA_A > THETA_LIM)
               ! If col == 1, sample the next collision site.
               ! I also use a threshold for where the next collision can't occur.
               ! This is when the probability is very low and we would waste significant time finding local normal coordinates.
               ! Exactly the same as before but functions are specific for upwards collisions.
               DOWN = -1
               NUP = NUP + 1
               !print*, down, vx_i,vy_i,vz_i
               MAXX_NEG = INTERPOLATE(THETA_R_V, A_V, G_V, GRID_BC(IPG)%MAX_P_UP, THETA_A, A, G)/COS(THETA_A)    
               ANGLE_VALUES = ACCEPT_REJECT_SURF(MAXX_NEG, THETA_A, PHI_A, A, G,DOWN)
               VELOCITIES(:) = CL_KERNEL_LOCAL(VX_I, VY_I, VZ_I, ALPHA_N, ALPHA_T, CI,ANGLE_VALUES(1),ANGLE_VALUES(2))
               VX_I =  VELOCITIES(1) 
               VY_I =  VELOCITIES(2)
               VZ_I =  VELOCITIES(3)
               THETA_A = ACOS(VZ_I/SQRT(VX_I*VX_I + VY_I*VY_I + VZ_I*VZ_I))    
               PHI_A = ATAN2(VY_I, VX_I)
               ! If after upwards collision the velocity is positive, we check for collision
               IF (VZ_I > 0) THEN
                  R = RF()
                  IF (R <= INTERPOLATE(THETA_R_V, A_V, G_V, GRID_BC(IPG)%P_COLL_UP, THETA_A, A, G)  ) THEN 
                        COL = 1
                        !WRITE(*,*) 'Upwards collision 2!'
                  ELSE 
                        COL = 0
                        IF (VZ_I < LIMIT_TRAP .AND. TRAPPING == 1) THEN 
                           VZ_I = -VZ_I
                        END IF 
                  END IF 
               ELSE
                  COL = 1
               END IF
            END DO
         END IF 
         ! If we reach enough collisions we just assume full accommodation.
         IF (NDOWN + NUP > COL_TO_TRAP .AND. TRAPPING == 1) THEN
               VELOCITIES(:) = CL_KERNEL_LOCAL(VX_I, VY_I, VZ_I, 1.0D0, 1.0D0, CI, 0.0D0, 0.0D0)
               !TD = TD + 1
               EXIT
         END IF 

         IF (VZ_I > 0) THEN
               IF (VZ_I > LIMIT_TRAP) THEN 
                  VELOCITIES(3) = SQRT(VZ_I**2-2*W/MASS)
               ELSE 
                  IF (TRAPPING == 1 .AND. THETA_A < THETA_LIM ) THEN
                     VZ_I = - VZ_I
                  ELSE IF (TRAPPING == 0) THEN
                     VELOCITIES(:) = CL_KERNEL_LOCAL(VX_I, VY_I, VZ_I, 1.0D0, 1.0D0, CI, 0.0D0, 0.0D0)
                     TD = TD + 1
                     EXIT 
                  END IF 
               END IF 
         END IF 
      END DO

      VX_I =  VELOCITIES(1) 
      VY_I =  VELOCITIES(2)
      VZ_I =  VELOCITIES(3)
      !IF (NUP > 0) WRITE(*,*) 'NUP = ', nup

   END SUBROUTINE WB_SCATTER


   !Function analogous to numpy.linspace.
   FUNCTION LINSPACE(START, END, NUM_POINTS) RESULT(RESULT)
      REAL(KIND=8), INTENT(IN) :: START, END
      INTEGER, INTENT(IN) :: NUM_POINTS
      REAL(KIND=8), DIMENSION(NUM_POINTS) :: RESULT
      
      REAL(KIND=8) :: STEP
      INTEGER :: I
      
      STEP = (END - START) / REAL(NUM_POINTS - 1)
      
      DO I = 1, NUM_POINTS
          RESULT(I) = START + REAL(I - 1) * STEP
      END DO
   END FUNCTION LINSPACE



   FUNCTION SURFACE_PROF(THETA_I, PHI, ALPHA, BETA, A, G,DOWN)
      REAL(KIND=8) :: SURFACE_PROF, THETA_I, PHI, ALPHA, BETA, A, G
      REAL(KIND=8) :: SEC, ARG, EE2, DEN, TERM
      INTEGER :: DOWN

      SEC = 1.0 / COS(ALPHA)
      IF (G>=1) THEN
         ARG = 1.0 - 1.0 / G**2
         ! Declare return type of m_ellipE as real when calling it
         EE2 = (M_ELLIPE(ARG))**2

         DEN = PI**2 * A**2 * COS(ALPHA)**3
         TERM = G**2 * EE2 / (A**2 * PI)
         
         SURFACE_PROF = G * EE2 * SIN(ALPHA) / DEN * EXP(-TERM * (SEC**2 - 1.0) * (COS(BETA)**2 + 1.0 / G**2 * SIN(BETA)**2)) 
      ELSE IF (G<1) THEN
         ARG = 1.0 -G**2
         EE2 = (M_ELLIPE(ARG))**2
   
         DEN = G*PI**2*A**2*COS(ALPHA)**3
         TERM = EE2/(A**2*PI)
         
   
         SURFACE_PROF = EE2 * SIN(ALPHA) / DEN * EXP(- TERM * (SEC**2 - 1) * (COS(BETA)**2 + 1.0 / G**2 * SIN(BETA)**2))
      END IF

      IF (DOWN == 1) THEN
         SURFACE_PROF = SURFACE_PROF * MAX(0.0, (TAN(THETA_I) * TAN(ALPHA) * COS(PHI - BETA) + 1.0))
      ELSE IF (DOWN == -1) THEN
         SURFACE_PROF = SURFACE_PROF * MAX(0.0, -(TAN(THETA_I) * TAN(ALPHA) * COS(PHI - BETA) + 1.0))
      END IF
   
   END FUNCTION SURFACE_PROF

   !Approximation of the elliptical integral in the interval 0-1. 
   FUNCTION M_ELLIPE(X)
      REAL(KIND=8) :: M_ELLIPE, X
      M_ELLIPE = (1.56969 - 2.24458*X + 0.728559*X**2)/(1.0 - 1.18949*X + 0.243123*X**2)
   END FUNCTION M_ELLIPE

   !Acceptance rejection algorithm 
   FUNCTION ACCEPT_REJECT_SURF(MAXL, THETA_I, PHI, A, G,DOWN) RESULT(RANDOM_VALUES)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: MAXL, THETA_I, PHI, A, G
      REAL(KIND=8) :: ALPHA_TEST, BETA_TEST, Y
      REAL(KIND=8), DIMENSION(2) :: RANDOM_VALUES
      INTEGER :: DOWN

      !call random_seed()    
      DO
         ! Initialize random number generator

         ! Generate random alpha_test
         ALPHA_TEST = rf()
         ALPHA_TEST = ALPHA_TEST * (PI / 2.0)

         ! Generate random beta_teste
         BETA_TEST = rf()
         BETA_TEST = (2.0 * PI) * (BETA_TEST - 0.5)

         ! Generate random y
         Y = rf()
         Y = Y * MAXL

         IF (Y < SURFACE_PROF(THETA_I, PHI, ALPHA_TEST, BETA_TEST, A, G, DOWN)) THEN 
            EXIT 
         END IF
      END DO 
      RANDOM_VALUES = [ALPHA_TEST, BETA_TEST]
   END FUNCTION ACCEPT_REJECT_SURF

   ! Obtain velocity coordinates in a reference frame whose normal is given by the
   ! polar angle, alpha, and azimuthal angle beta
   ! Align the velocity so that there is only one tangential velocity component.
   FUNCTION LAB_TO_LOCAL_MINE_PHI(VX, VY, VZ, ALPHA, BETA)
      REAL(KIND=8), INTENT(IN) :: VX, VY, VZ, ALPHA, BETA
      REAL(KIND=8) :: COS_ALPHA, SIN_ALPHA, COS_BETA, SIN_BETA
      REAL(KIND=8) :: UX, UY, UZ, PSI
      REAL(KIND=8), DIMENSION(4) :: LAB_TO_LOCAL_MINE_PHI

      UX = SQRT(((VX*SIN(BETA) - VY*COS(BETA))**2 + &
         (VX*COS(ALPHA)*COS(BETA) + VY*SIN(BETA)*COS(ALPHA) - VZ*SIN(ALPHA))**2) / &
         (VX*COS(ALPHA)*COS(BETA) + VY*SIN(BETA)*COS(ALPHA) - VZ*SIN(ALPHA))**2) * &
         (VX*COS(ALPHA)*COS(BETA) + VY*SIN(BETA)*COS(ALPHA) - VZ*SIN(ALPHA))

      UY = 0

      UZ = SIN(ALPHA)*COS(BETA)*VX + SIN(ALPHA)*SIN(BETA)*VY + COS(ALPHA)*VZ

      PSI = ATAN((-VX*SIN(BETA) + VY*COS(BETA)) / &
      (VX*COS(ALPHA)*COS(BETA) + VY*SIN(BETA)*COS(ALPHA) - VZ*SIN(ALPHA)))

      LAB_TO_LOCAL_MINE_PHI = [UX, UY, UZ, PSI]
   END FUNCTION LAB_TO_LOCAL_MINE_PHI

   ! Inverse of the previous expresion.
   FUNCTION LOCAL_TO_LAB_MINE_PHI(UX, UY, UZ, ALPHA, BETA, PSI)
      REAL(KIND=8), INTENT(IN) :: UX, UY, UZ, ALPHA, BETA, PSI
      REAL(KIND=8) :: VX_R, VY_R, VZ_R
      REAL(KIND=8), DIMENSION(3) :: LOCAL_TO_LAB_MINE_PHI


      VX_R = UX*(-SIN(BETA)*SIN(PSI) + COS(ALPHA)*COS(BETA)*COS(PSI)) + &
            UY*(-SIN(BETA)*COS(PSI) - SIN(PSI)*COS(ALPHA)*COS(BETA)) + &
            UZ*SIN(ALPHA)*COS(BETA)

      VY_R = UX*(SIN(BETA)*COS(ALPHA)*COS(PSI) + SIN(PSI)*COS(BETA)) + &
            UY*(-SIN(BETA)*SIN(PSI)*COS(ALPHA) + COS(BETA)*COS(PSI)) + &
            UZ*SIN(ALPHA)*SIN(BETA)
      VZ_R = -UX*SIN(ALPHA)*COS(PSI) + UY*SIN(ALPHA)*SIN(PSI) + UZ*COS(ALPHA)

      LOCAL_TO_LAB_MINE_PHI = [VX_R, VY_R, VZ_R]
   END FUNCTION LOCAL_TO_LAB_MINE_PHI


   !Basic CL Kernel. 
   FUNCTION CL_KERNEL(UN, UT, ALPHA_N, ALPHA_T, CI) !RESULT(KERNEL_VALUES)
      REAL(KIND=8), INTENT(IN) :: UN, UT, ALPHA_N, ALPHA_T, CI
      REAL(KIND=8) :: AL, AM, AN
      REAL(KIND=8) :: RAND1, RAND2, RAND3, RAND4, RAND5, RAND6
      REAL(KIND=8) :: R1, R3, R5, PHI2, PHI4, PHI6, VNM, VTM
      !REAL(KIND=8), dimension(3) :: kernel_values ! Declare the function result type
      REAL(KIND=8), DIMENSION(3) :: CL_KERNEL

      ! Generate random numbers
      RAND1 = rf()
      RAND2 = rf()
      RAND3 = rf()
      RAND4 = rf()
      RAND5 = rf()
      RAND6 = rf()

      ! Calculate r1, r3, r5
      R1 = SQRT(-ALPHA_N * LOG(RAND1))
      R3 = SQRT(-ALPHA_T * LOG(RAND2))
      R5 = SQRT(-ALPHA_T * LOG(RAND3))

      ! Calculate phi2, phi4, phi6
      PHI2 = 2 * PI * RAND4 ! Use the declared constant PI
      PHI4 = 2 * PI * RAND5
      PHI6 = 2 * PI * RAND6

      ! Calculate vnm, vtm
      VNM = UN / CI * SQRT(1 - ALPHA_N)
      VTM = UT / CI * SQRT(1 - ALPHA_T)

      ! Calculate AM, AL, AN
      AM = CI * SQRT(R1**2 + VNM**2 + 2*R1*VNM*COS(PHI2))
      AL = CI * (VTM + R3*COS(PHI4))
      AN = CI * (R5*COS(PHI6))

      CL_KERNEL = [AL, AN, AM]
   END FUNCTION CL_KERNEL

   ! CL Kernel with change to local coordinate frame. A
   ! Alpha and Beta are the polar and azimuthal angles of the local normal
   ! Velocities are changed back to the macroscopic normal frame and returned
   FUNCTION CL_KERNEL_LOCAL(VX, VY, VZ, ALPHA_N, ALPHA_T, CI, ALPHA, BETA) !RESULT(KERNEL_VALUES)
      REAL(KIND=8) :: VX, VY, VZ, ALPHA_N, ALPHA_T, CI, ALPHA, BETA
      REAL(KIND=8) :: UI_V(4), UR_V(3), VR_V(3)
      REAL(KIND=8) :: CL_KERNEL_LOCAL(3)

      UI_V = LAB_TO_LOCAL_MINE_PHI(VX, VY, VZ, ALPHA, BETA)
      UR_V = CL_KERNEL(UI_V(3), UI_V(1), ALPHA_N, ALPHA_T, CI)
      CL_KERNEL_LOCAL = LOCAL_TO_LAB_MINE_PHI(UR_V(1), UR_V(2), UR_V(3), ALPHA, BETA, UI_V(4))
   END FUNCTION CL_KERNEL_LOCAL

    ! Interpolates data that has a 3D input.
    ! x,y,z are the grid coordinates where the data was computed 
    ! data is a 3D matrix 
    ! (x,y,z)_interp are the coordinates where we want to interpolate the existing data.
    ! NOTE: This is not the most efficient, please use an interpolation library 
    ! whose indices and coefficients are computed once before and then just call 
    ! an evaluator when needed. This is not efficient because the interpolation 
    ! coefficients are computed on every call to the function. 
   FUNCTION INTERPOLATE(X, Y, Z, DATA, X_INTERP, Y_INTERP, Z_INTERP) RESULT(RESULT)
      REAL(KIND=8), INTENT(IN) :: X(:), Y(:), Z(:), DATA(:,:,:), X_INTERP, Y_INTERP, Z_INTERP
      REAL(KIND=8) :: RESULT
      REAL(KIND=8) :: F1, F2, F3, F4, F5, F6, F7, F8
      INTEGER :: I, J, K, I1, J1, K1, I2, J2, K2
      REAL(KIND=8) :: X_FRAC, Y_FRAC, Z_FRAC
      
      ! Find the indices surrounding the interpolation point
      CALL FIND_INDICES(X, X_INTERP, I1, I2, X_FRAC)
      CALL FIND_INDICES(Y, Y_INTERP, J1, J2, Y_FRAC)
      CALL FIND_INDICES(Z, Z_INTERP, K1, K2, Z_FRAC)
      ! i1 = 1
      ! i2 = 2
      ! j1 = 1
      ! j2 = 2
      ! k1 = 1
      ! k2 = 2
      ! x_frac = 1
      ! y_frac = 1
      ! z_frac = 1

      ! Perform trilinear interpolation
      F1 = DATA(I1,J1,K1) * (1.0 - X_FRAC) + DATA(I2,J1,K1) * X_FRAC
      F2 = DATA(I1,J2,K1) * (1.0 - X_FRAC) + DATA(I2,J2,K1) * X_FRAC
      F3 = DATA(I1,J1,K2) * (1.0 - X_FRAC) + DATA(I2,J1,K2) * X_FRAC
      F4 = DATA(I1,J2,K2) * (1.0 - X_FRAC) + DATA(I2,J2,K2) * X_FRAC
      
      F5 = F1 * (1.0 - Y_FRAC) + F2 * Y_FRAC
      F6 = F3 * (1.0 - Y_FRAC) + F4 * Y_FRAC
      
      F7 = F5 * (1.0 - Z_FRAC) + F6 * Z_FRAC
      
      RESULT = F7
      
   END FUNCTION INTERPOLATE

   ! ! Find indices and fractional part for linear interpolation
   ! subroutine find_indices(arr, val, idx1, idx2, frac)
   !    REAL(KIND=8), intent(in) :: arr(:), val
   !    integer, intent(out) :: idx1, idx2
   !    REAL(KIND=8), intent(out) :: frac
   !    integer :: i



      
   !    if (val <= arr(1)) then
   !        idx1 = 1
   !        idx2 = 2
   !        frac = 0.0
   !    else if (val >= arr(size(arr))) then
   !        idx1 = size(arr) - 1
   !        idx2 = size(arr)
   !        frac = 1.0
   !    else
   !       !ARRAY(INDEX) < VALUE < ARRAY(INDEX+1)
   !       idx1 = BINARY_SEARCH(val, arr)
   !       idx2 = idx1+1
   !       frac = (val - arr(i)) / (arr(i+1) - arr(i))
   !    end if
   ! end subroutine find_indices

   ! FIND INDICES AND FRACTIONAL PART FOR LINEAR INTERPOLATION
   SUBROUTINE FIND_INDICES(ARR, VAL, IDX1, IDX2, FRAC)
      REAL(KIND=8), INTENT(IN) :: ARR(:), VAL
      INTEGER, INTENT(OUT) :: IDX1, IDX2
      REAL(KIND=8), INTENT(OUT) :: FRAC
      INTEGER :: I
      
      IF (VAL <= ARR(1)) THEN
         IDX1 = 1
         IDX2 = 2
         FRAC = 0.0
         WRITE(*,*) 'Value outside the interpolation range!'
      ELSE IF (VAL >= ARR(SIZE(ARR))) THEN
         IDX1 = SIZE(ARR) - 1
         IDX2 = SIZE(ARR)
         FRAC = 1.0
         WRITE(*,*) 'Value outside the interpolation range!'
      ELSE
         DO I = 1, SIZE(ARR) - 1
            IF (VAL >= ARR(I) .AND. VAL <= ARR(I+1)) THEN
               IDX1 = I
               IDX2 = I + 1
               FRAC = (VAL - ARR(I)) / (ARR(I+1) - ARR(I))
               EXIT
            END IF
         END DO
      END IF
   END SUBROUTINE FIND_INDICES


END MODULE washboard
