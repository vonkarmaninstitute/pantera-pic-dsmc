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
      theta_lim = 35./180.*pi


      !Get the grid where the pre-computed values were taken
      theta_i_v = linspace(pi/2, pi, 45)    ! DBDBDBBDBDBDBDDBBDDB include in the structure
      theta_r_v = linspace(0.0d0, pi/2, 45)
      !Careful!!! 
      !This is the range of validity of the precomputed values for 
      !the maximums of the probability distribution of the local surface normals.
      ! If you take values for A and g out of this range the sampling procedure is not guaranteed to be correct!
      A_v = linspace(0.5d0, 1.5d0, 40)
      g_v1 = linspace(0.2d0, 1.0d0, 20)       
      g_v2 = linspace(1.05d0, 5.0d0, 20)    
  
      ! Concatenate g_v1 and g_v2 into g_v
      j = 1
      do i = 1, size(g_v1)
         g_v(j) = g_v1(i)
         j = j + 1
      end do
      do i = 1, size(g_v2)
         g_v(j) = g_v2(i)
         j = j + 1
      end do


      VZ_I = -SQRT(VZ_I*VZ_I + 2*W/MASS)

      theta_a = ACOS(VZ_I/SQRT(VZ_I*VZ_I + VX_I*VX_I))
      phi_a = 0



      do while (VZ_I < 0)
         down = 1
         ndown = ndown + 1
         !print*, down, vx_i,vy_i,vz_i
         ! Sample local normal
         theta_a = ACOS(VZ_I/SQRT(VX_I*VX_I + VY_I*VY_I + VZ_I*VZ_I))
         phi_a = atan2(vy_i, vx_i)
         
         maxx = interpolate(theta_i_v, A_v, g_v, GRID_BC(IPG)%MAX_P_DN, theta_a, A, g)/cos(theta_a)
         angle_values = accept_reject_surf(maxx, theta_a, phi_a, A, g,down)
         
         velocities(:) = CL_Kernel_local(vx_i, vy_i, vz_i, alpha_n, alpha_t, ci,angle_values(1),angle_values(2))
         vx_i =  velocities(1) 
         vy_i =  velocities(2)
         vz_i =  velocities(3)

         theta_a = ACOS(VZ_I/SQRT(VX_I*VX_I + VY_I*VY_I + VZ_I*VZ_I))  
         phi_a = atan2(vy_i, vx_i)
         if (vz_i > 0) then
               R = rf()
               itpl = interpolate(theta_r_v, A_v, g_v, GRID_BC(IPG)%P_COLL_UP, theta_a, A, g)
               !WRITE(*,*) 'ITPL = ', ITPL
               if (R <= itpl  ) then 
                  col = 1
                  !WRITE(*,*) 'Upwards collision 1!'
               else 
                  col = 0
                  if (vz_i < limit_trap .And. trapping == 1) then 
                     vz_i = -vz_i
                  end if 
               end if 
               ! col = 0 -> no collision | col = 1 -> collision
               ! Upwards collision loop
               do while (col == 1 .AND. theta_a > theta_lim)
                  ! If col == 1, sample the next collision site.
                  ! I also use a threshold for where the next collision can't occur.
                  ! This is when the probability is very low and we would waste significant time finding local normal coordinates.
                  ! Exactly the same as before but functions are specific for upwards collisions.
                  down = -1
                  nup = nup + 1
                  !print*, down, vx_i,vy_i,vz_i
                  maxx_neg = interpolate(theta_r_v, A_v, g_v, GRID_BC(IPG)%MAX_P_UP, theta_a, A, g)/cos(theta_a)    
                  angle_values = accept_reject_surf(maxx_neg, theta_a, phi_a, A, g,down)
                  velocities(:) = CL_Kernel_local(vx_i, vy_i, vz_i, alpha_n, alpha_t, ci,angle_values(1),angle_values(2))
                  vx_i =  velocities(1) 
                  vy_i =  velocities(2)
                  vz_i =  velocities(3)
                  theta_a = ACOS(VZ_I/SQRT(VX_I*VX_I + VY_I*VY_I + VZ_I*VZ_I))    
                  phi_a = atan2(vy_i, vx_i)
                  ! If after upwards collision the velocity is positive, we check for collision
                  if (vz_i > 0) then
                     R = rf()
                     if (R <= interpolate(theta_r_v, A_v, g_v, GRID_BC(IPG)%P_COLL_UP, theta_a, A, g)  ) then 
                           col = 1
                           !WRITE(*,*) 'Upwards collision 2!'
                     else 
                           col = 0
                           if (vz_i < limit_trap .And. trapping == 1) then 
                              vz_i = -vz_i
                           end if 
                     end if 
                  else
                     col = 1
                  end if
               end do
         end if 
         ! If we reach enough collisions we just assume full accommodation.
         if (ndown + nup > col_to_trap .AND. trapping == 1) then
               velocities(:) = CL_Kernel_local(vx_i, vy_i, vz_i, 1.0d0, 1.0d0, ci, 0.0d0, 0.0d0)
               !TD = TD + 1
               exit
         end if 

         if (vz_i > 0) then
               if (vz_i > limit_trap) then 
                  velocities(3) = sqrt(vz_i**2-2*W/MASS)
               else 
                  if (trapping == 1 .AND. theta_a < theta_lim ) then
                     vz_i = - vz_i
                  else if (trapping == 0) then
                     velocities(:) = CL_Kernel_local(vx_i, vy_i, vz_i, 1.0d0, 1.0d0, ci, 0.0d0, 0.0d0)
                     TD = TD + 1
                     exit 
                  end if 
               end if 
         end if 
      end do

      vx_i =  velocities(1) 
      vy_i =  velocities(2)
      vz_i =  velocities(3)
      !IF (NUP > 0) WRITE(*,*) 'NUP = ', nup

   END SUBROUTINE WB_SCATTER


   !Function analogous to numpy.linspace.
   function linspace(start, end, num_points) result(result)
      REAL(KIND=8), intent(in) :: start, end
      integer, intent(in) :: num_points
      REAL(KIND=8), dimension(num_points) :: result
      
      REAL(KIND=8) :: step
      integer :: i
      
      step = (end - start) / real(num_points - 1)
      
      do i = 1, num_points
          result(i) = start + real(i - 1) * step
      end do
   end function linspace



   function surface_prof(theta_i, phi, alpha, beta, A, g,down)
      REAL(KIND=8) :: surface_prof, theta_i, phi, alpha, beta, A, g
      REAL(KIND=8) :: sec, arg, Ee2, den, term
      integer :: down

      sec = 1.0 / cos(alpha)
      if (g>=1) then
         arg = 1.0 - 1.0 / g**2
         ! Declare return type of m_ellipE as real when calling it
         Ee2 = (m_ellipE(arg))**2

         den = pi**2 * A**2 * cos(alpha)**3
         term = g**2 * Ee2 / (A**2 * pi)
         
         surface_prof = g * Ee2 * sin(alpha) / den * exp(-term * (sec**2 - 1.0) * (cos(beta)**2 + 1.0 / g**2 * sin(beta)**2)) 
      else if (g<1) then
         arg = 1.0 -g**2
         Ee2 = (m_ellipE(arg))**2
   
         den = g*pi**2*A**2*cos(alpha)**3
         term = Ee2/(A**2*pi)
         
   
         surface_prof = Ee2 * sin(alpha) / den * exp(- term * (sec**2 - 1) * (cos(beta)**2 + 1.0 / g**2 * sin(beta)**2))
      end if

      if (down == 1) then
         surface_prof = surface_prof * max(0.0, (tan(theta_i) * tan(alpha) * cos(phi - beta) + 1.0))
      else if (down == -1) then
         surface_prof = surface_prof * max(0.0, -(tan(theta_i) * tan(alpha) * cos(phi - beta) + 1.0))
      end if
   
   end function surface_prof

  !Approximation of the elliptical integral in the interval 0-1. 
  function m_ellipE(x)
      REAL(KIND=8) :: m_ellipE, x
      m_ellipE = (1.56969 - 2.24458*x + 0.728559*x**2)/(1.0 - 1.18949*x + 0.243123*x**2)
  end function m_ellipE

  !Acceptance rejection algorithm 
  function accept_reject_surf(maxl, theta_i, phi, A, g,down) result(random_values)
      implicit none
      REAL(KIND=8), intent(in) :: maxl, theta_i, phi, A, g
      REAL(KIND=8) :: alpha_test, beta_test, y
      REAL(KIND=8), dimension(2) :: random_values
      integer :: down

      !call random_seed()    
      do 
         ! Initialize random number generator

         ! Generate random alpha_test
         alpha_test = rf()
         alpha_test = alpha_test * (pi / 2.0)

         ! Generate random beta_teste
         beta_test = rf()
         beta_test = (2.0 * pi) * (beta_test - 0.5)

         ! Generate random y
         y = rf()
         y = y * maxl

         if (y < surface_prof(theta_i, phi, alpha_test, beta_test, A, g,down)) then 
            exit 
         end if
      end do 
      random_values = [alpha_test, beta_test]
  end function accept_reject_surf

      ! Obtain velocity coordinates in a reference frame whose normal is given by the
    ! polar angle, alpha, and azimuthal angle beta
    ! Align the velocity so that there is only one tangential velocity component.
  function Lab_to_local_mine_phi(vx, vy, vz, alpha, beta)
   REAL(KIND=8), intent(in) :: vx, vy, vz, alpha, beta
   REAL(KIND=8) :: cos_alpha, sin_alpha, cos_beta, sin_beta
   REAL(KIND=8) :: ux, uy, uz, psi
   REAL(KIND=8), dimension(4) :: Lab_to_local_mine_phi

   ux = SQRT(((vx*SIN(beta) - vy*COS(beta))**2 + &
      (vx*COS(alpha)*COS(beta) + vy*SIN(beta)*COS(alpha) - vz*SIN(alpha))**2) / &
      (vx*COS(alpha)*COS(beta) + vy*SIN(beta)*COS(alpha) - vz*SIN(alpha))**2) * &
      (vx*COS(alpha)*COS(beta) + vy*SIN(beta)*COS(alpha) - vz*SIN(alpha))

   uy = 0

   uz = SIN(alpha)*COS(beta)*vx + SIN(alpha)*SIN(beta)*vy + COS(alpha)*vz

   psi = ATAN((-vx*SIN(beta) + vy*COS(beta)) / &
(vx*COS(alpha)*COS(beta) + vy*SIN(beta)*COS(alpha) - vz*SIN(alpha)))

   Lab_to_local_mine_phi = [ux, uy, uz, psi]
end function Lab_to_local_mine_phi

! Inverse of the previous expresion.
function Local_to_lab_mine_phi(ux, uy, uz, alpha, beta, psi)
   REAL(KIND=8), intent(in) :: ux, uy, uz, alpha, beta, psi
   REAL(KIND=8) :: vx_r, vy_r, vz_r
   REAL(KIND=8), dimension(3) :: Local_to_lab_mine_phi


   vx_r = ux*(-SIN(beta)*SIN(psi) + COS(alpha)*COS(beta)*COS(psi)) + &
          uy*(-SIN(beta)*COS(psi) - SIN(psi)*COS(alpha)*COS(beta)) + &
          uz*SIN(alpha)*COS(beta)

   vy_r = ux*(SIN(beta)*COS(alpha)*COS(psi) + SIN(psi)*COS(beta)) + &
          uy*(-SIN(beta)*SIN(psi)*COS(alpha) + COS(beta)*COS(psi)) + &
          uz*SIN(alpha)*SIN(beta)
   vz_r = -ux*SIN(alpha)*COS(psi) + uy*SIN(alpha)*SIN(psi) + uz*COS(alpha)

   Local_to_lab_mine_phi = [vx_r, vy_r, vz_r]
end function Local_to_lab_mine_phi


!Basic CL Kernel. 
function CL_Kernel(un, ut, alpha_n, alpha_t, ci) !result(kernel_values)
   REAL(KIND=8), intent(in) :: un, ut, alpha_n, alpha_t, ci
   REAL(KIND=8) :: AL, AM, AN
   REAL(KIND=8) :: rand1, rand2, rand3, rand4, rand5, rand6
   REAL(KIND=8) :: r1, r3, r5, phi2, phi4, phi6, vnm, vtm
   !REAL(KIND=8), dimension(3) :: kernel_values ! Declare the function result type
   REAL(KIND=8), dimension(3) :: CL_kernel

   ! Generate random numbers
   rand1 = rf()
   rand2 = rf()
   rand3 = rf()
   rand4 = rf()
   rand5 = rf()
   rand6 = rf()

   ! Calculate r1, r3, r5
   r1 = sqrt(-alpha_n * log(rand1))
   r3 = sqrt(-alpha_t * log(rand2))
   r5 = sqrt(-alpha_t * log(rand3))

   ! Calculate phi2, phi4, phi6
   phi2 = 2 * PI * rand4 ! Use the declared constant PI
   phi4 = 2 * PI * rand5
   phi6 = 2 * PI * rand6

   ! Calculate vnm, vtm
   vnm = un / ci * sqrt(1 - alpha_n)
   vtm = ut / ci * sqrt(1 - alpha_t)

   ! Calculate AM, AL, AN
   AM = ci * sqrt(r1**2 + vnm**2 + 2*r1*vnm*cos(phi2))
   AL = ci * (vtm + r3*cos(phi4))
   AN = ci * (r5*cos(phi6))

   CL_Kernel = [ AL,  AN, AM ]
end function CL_Kernel

! CL Kernel with change to local coordinate frame. A
! Alpha and Beta are the polar and azimuthal angles of the local normal
! Velocities are changed back to the macroscopic normal frame and returned
function CL_Kernel_local(vx, vy, vz, alpha_n, alpha_t, ci, alpha, beta) !result(kernel_values)
   REAL(KIND=8) :: vx, vy, vz, alpha_n, alpha_t, ci, alpha, beta
   REAL(KIND=8) :: ui_v(4), ur_v(3), vr_v(3)
   REAL(KIND=8) :: CL_kernel_local(3)

   ui_v = Lab_to_local_mine_phi(vx, vy, vz, alpha, beta)
   ur_v = CL_Kernel(ui_v(3), ui_v(1), alpha_n, alpha_t, ci)
   CL_Kernel_local = Local_to_lab_mine_phi(ur_v(1), ur_v(2), ur_v(3), alpha, beta, ui_v(4))
end function CL_Kernel_local

    ! Interpolates data that has a 3D input.
    ! x,y,z are the grid coordinates where the data was computed 
    ! data is a 3D matrix 
    ! (x,y,z)_interp are the coordinates where we want to interpolate the existing data.
    ! NOTE: This is not the most efficient, please use an interpolation library 
    ! whose indices and coefficients are computed once before and then just call 
    ! an evaluator when needed. This is not efficient because the interpolation 
    ! coefficients are computed on every call to the function. 
function interpolate(x, y, z, data, x_interp, y_interp, z_interp) result(result)
   REAL(KIND=8), intent(in) :: x(:), y(:), z(:), data(:,:,:), x_interp, y_interp, z_interp
   REAL(KIND=8) :: result
   REAL(KIND=8) :: f1, f2, f3, f4, f5, f6, f7, f8
   INTEGER :: i, j, k, i1, j1, k1, i2, j2, k2
   REAL(KIND=8) :: x_frac, y_frac, z_frac
   
   ! Find the indices surrounding the interpolation point
   call find_indices(x, x_interp, i1, i2, x_frac)
   call find_indices(y, y_interp, j1, j2, y_frac) ! these are slow af
   call find_indices(z, z_interp, k1, k2, z_frac)
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
   f1 = data(i1,j1,k1) * (1.0 - x_frac) + data(i2,j1,k1) * x_frac
   f2 = data(i1,j2,k1) * (1.0 - x_frac) + data(i2,j2,k1) * x_frac
   f3 = data(i1,j1,k2) * (1.0 - x_frac) + data(i2,j1,k2) * x_frac
   f4 = data(i1,j2,k2) * (1.0 - x_frac) + data(i2,j2,k2) * x_frac
   
   f5 = f1 * (1.0 - y_frac) + f2 * y_frac
   f6 = f3 * (1.0 - y_frac) + f4 * y_frac
   
   f7 = f5 * (1.0 - z_frac) + f6 * z_frac
   
   result = f7
   
end function interpolate

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

    ! Find indices and fractional part for linear interpolation
subroutine find_indices(arr, val, idx1, idx2, frac)
   real(KIND=8), intent(in) :: arr(:), val
   integer, intent(out) :: idx1, idx2
   real(KIND=8), intent(out) :: frac
   integer :: i
   
   if (val <= arr(1)) then
      idx1 = 1
      idx2 = 2
      frac = 0.0
      WRITE(*,*) 'Value outside the interpolation range!'
   else if (val >= arr(size(arr))) then
      idx1 = size(arr) - 1
      idx2 = size(arr)
      frac = 1.0
      WRITE(*,*) 'Value outside the interpolation range!'
   else
      do i = 1, size(arr) - 1
         if (val >= arr(i) .and. val <= arr(i+1)) then
            idx1 = i
            idx2 = i + 1
            frac = (val - arr(i)) / (arr(i+1) - arr(i))
            exit
         end if
      end do
   end if
end subroutine find_indices


END MODULE washboard