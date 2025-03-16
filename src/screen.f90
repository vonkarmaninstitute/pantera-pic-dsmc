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

! This module contains functions to conveniently print stuff on the screen

MODULE screen

   IMPLICIT NONE 

   CONTAINS 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE PRINTTITLE(P_ID) !!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE PRINTTITLE(P_ID)

      INTEGER, INTENT(IN) :: P_ID ! ID of current process ID

      IF (P_ID .EQ. 0) THEN ! Only master prints
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)  'PANTERA PIC-DSMC  Copyright (C) 2025  von Karman Institute for Fluid Dynamics (VKI)'
         WRITE(*,*)  'This program comes with ABSOLUTELY NO WARRANTY.                                    '
         WRITE(*,*)  'This is free software, and you are welcome to redistribute it                      '
         WRITE(*,*)  'under certain conditions. See https://www.gnu.org/licenses/ for details.           '
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)  '                                                            ___.,.   .~-.,    '
         WRITE(*,*)  '                                              .___    .,//.~      `\/     \   '
         WRITE(*,*)  '                                             / `\ \.,/             /       |  '
         WRITE(*,*)  '                                             |   )               ./  (     )  '
         WRITE(*,*)  '                                              \ /                     \,` /   '
         WRITE(*,*)  '                                               )                          )   '
         WRITE(*,*)  '                                               (                              '
         WRITE(*,*)  '                                               /     ~.     .`                '
         WRITE(*,*)  '                                              (  __.  `\   /  ,__,            '
         WRITE(*,*)  '                                               )(o_\.   `,`  ./o_)            '
         WRITE(*,*)  '     _________                                (     )  / /  //                '
         WRITE(*,*)  '     \   ___  \                                \   /    `                     '
         WRITE(*,*)  '     /  /___\ /       __                        \_/                           '
         WRITE(*,*)  '    /  ______/_____  / /____  _________ _        (   /`/ ,``, )/______        '
         WRITE(*,*)  '   /  / / __``/ __ \/ __/ _ \/ ___/ __ `/    ____( .____ `   /\______         '
         WRITE(*,*)  '  /  / / /_/ / / / / /_/  __/ /  / /_/ /    ______\(,  ,)   /V \              '
         WRITE(*,*)  ' /___| \__,_/_/ /_/\__/\___/_/   \__,_/        ____\\__/__~`V   ) ./          '
         WRITE(*,*)  '                                  V  V              V_M._V.M/\/` /            '
         WRITE(*,*)  ' PArticle Numerical Tool for non-Equilibrium         \          /             '
         WRITE(*,*)  '                            Reacting Aerodynamics     (__.,-_.,/              '
         WRITE(*,*)
         WRITE(*,*)  '      ____     ____   ______           ____    _____    __  ___   ______      '
         WRITE(*,*)  '     / __ \   /  _/  / ____/          / __ \  / ___/   /  |/  /  / ____/      '
         WRITE(*,*)  '    / /_/ /   / /   / /      ______  / / / /  \__ \   / /|_/ /  / /           '
         WRITE(*,*)  '   / ____/  _/ /   / /___   /_____/ / /_/ /  ___/ /  / /  / /  / /___         '
         WRITE(*,*)  '  /_/      /___/   \____/          /_____/  /____/  /_/  /_/   \____/         '
         WRITE(*,*)
         WRITE(*,*)
      END IF

   END SUBROUTINE PRINTTITLE


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ONLYMASTERPRINT1 -> Prints a string (master only) !!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE ONLYMASTERPRINT1(P_ID, string)
 
      IMPLICIT NONE
    
      INTEGER, INTENT(IN)           :: P_ID
      CHARACTER(LEN=*), INTENT(IN)  :: string
    
      IF (P_ID .EQ. 0) THEN
         WRITE(*,*) string
      END IF
 
   END SUBROUTINE ONLYMASTERPRINT1

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ONLYMASTERPRINT2 -> Prints string and a number (master only) !!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE ONLYMASTERPRINT2(P_ID, string, num)
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN)           :: P_ID
      CHARACTER(LEN=*), INTENT(IN)  :: string
      REAL(KIND=8), INTENT(IN)      :: num
 
      IF (P_ID .EQ. 0) THEN
        WRITE(*,'(A60,F15.3)') string, num
      END IF
 
   END SUBROUTINE ONLYMASTERPRINT2


   SUBROUTINE ERROR_ABORT(string)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: string

      PRINT*
      WRITE(*,*) '   ', string
      PRINT*
      STOP
      
   END SUBROUTINE ERROR_ABORT


END MODULE screen
