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

MODULE velocity_distribution

    USE mpi_common
    USE particle
    USE global
    USE tools

    IMPLICIT NONE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!    MAXWELL-BOLTZMANN VDF   !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define Maxwell-Boltzmann VDF type

    TYPE, EXTENDS(VELOCITY_DISTRIBUTION_STRUCTURE) :: MAXWELL_VDF
        CONTAINS
        PROCEDURE :: SAMPLE_VELOCITY => SAMPLE_VELOCITY_MAXWELL
        PROCEDURE :: FLX => FLX_MAXWELL
        PROCEDURE :: FLUXSOURCE => FLUXSOURCE_MAXWELL
        PROCEDURE :: BETA => BETA_MAXWELL
    END TYPE MAXWELL_VDF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!          KAPPA VDF         !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define Kappa VDF type
    TYPE, EXTENDS(VELOCITY_DISTRIBUTION_STRUCTURE) :: KAPPA_VDF
        REAL(KIND=8) :: KAPPA = 3.d0
        CONTAINS
        PROCEDURE :: SAMPLE_VELOCITY => SAMPLE_VELOCITY_KAPPA
        PROCEDURE :: FLX => FLX_KAPPA
        PROCEDURE :: FLUXSOURCE => FLUXSOURCE_KAPPA
        PROCEDURE :: BETA => BETA_KAPPA
    END TYPE KAPPA_VDF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!         MANUAL VDF         !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! =====  TO BE FINISHED ========
    ! This VDF should be manually loaded from a text file

    CONTAINS

    SUBROUTINE ASSIGN_VDF(VDF, VDF_NAME)
        CLASS(VELOCITY_DISTRIBUTION_STRUCTURE), ALLOCATABLE, INTENT(INOUT) :: VDF
        CHARACTER(*), INTENT(IN) :: VDF_NAME

        SELECT CASE (VDF_NAME)
        CASE('Maxwell')
          ALLOCATE(MAXWELL_VDF :: VDF)

        CASE('Kappa')
          ALLOCATE(KAPPA_VDF :: VDF)
          
        CASE DEFAULT
          ALLOCATE(MAXWELL_VDF :: VDF)
          CALL ONLYMASTERPRINT1(PROC_ID, 'Velocity distribution not specified for this task! Assuming Maxwell-Boltzmann VDF...')
      END SELECT
    END SUBROUTINE ASSIGN_VDF


    SUBROUTINE SAMPLE_VELOCITY_MAXWELL(THIS,UX, UY, UZ, TX, TY, TZ, VX, VY, VZ, M)

        IMPLICIT NONE

        CLASS(MAXWELL_VDF), INTENT(IN) :: THIS
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

    END SUBROUTINE SAMPLE_VELOCITY_MAXWELL

    FUNCTION FLX_MAXWELL(THIS, SN, TINF, M) RESULT(OUT)

        IMPLICIT NONE

        CLASS(MAXWELL_VDF), INTENT(IN) :: THIS
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

    END FUNCTION FLX_MAXWELL

    FUNCTION FLUXSOURCE_MAXWELL(THIS,S_NORM, U_NORM, TTRA, M) RESULT(OUT)
        CLASS(MAXWELL_VDF), INTENT(IN) :: THIS
        REAL(KIND=8), INTENT(INOUT) :: S_NORM
        REAL(KIND=8), INTENT(IN) :: U_NORM,TTRA,M
        REAL(KIND=8) :: OUT

        REAL(KIND=8) :: BETA

        BETA = 1./SQRT(2.*KB/M*TTRA)
        S_NORM = U_NORM*BETA
        OUT = 1.d0/(BETA*2.*SQRT(PI)) * (EXP(-S_NORM**2) &
                                + SQRT(PI)*S_NORM*(1.+ERF1(S_NORM)))  
    END FUNCTION FLUXSOURCE_MAXWELL

    FUNCTION BETA_MAXWELL(THIS, TTRA, M) RESULT(OUT)
        CLASS(MAXWELL_VDF), INTENT(IN) :: THIS

        REAL(KIND=8), INTENT(IN) :: TTRA, m
        REAL(KIND=8) :: OUT

        OUT = 1./SQRT(2.*KB*TTRA/M)

     END FUNCTION BETA_MAXWELL




    SUBROUTINE SAMPLE_VELOCITY_KAPPA(THIS,UX, UY, UZ, TX, TY, TZ, VX, VY, VZ, M)
  
        IMPLICIT NONE
     
        CLASS(KAPPA_VDF), INTENT(IN) :: THIS
        REAL(KIND=8), INTENT(IN)    :: UX, UY, UZ, TX, TY, TZ
        REAL(KIND=8), INTENT(IN)    :: M ! Molecular mass
        REAL(KIND=8), INTENT(INOUT) :: VX, VY, VZ
     
        INTEGER                     :: I
        REAL(KIND=8)                :: PI2
        REAL(KIND=8)                :: R, R1, RO, TETA, BETA
        REAL(KIND=8), DIMENSION(3)  :: VEL, TT
        REAL(KIND=8) :: DELTA, KAPPA_C
     
        KAPPA_C = THIS%KAPPA
     
        PI2  = 2.*PI
     
        TT(1) = TX
        TT(2) = TY
        TT(3) = TZ
     
     
        !!!!! KAPPA DISTRIBUTION !!!!!
     
        DO I = 1,3
           ! Step 1.
           R = rf()
                 
           TETA = PI2*R
     
           ! Step 2.
           ! R goes from 0 to 1 included. Remove the extremes
           R1 = rf()
              DO WHILE (R1 < 1.0D-13)
                 R1 = rf()
              END DO
     
           BETA = 1./SQRT(2.*KB*TT(I)/M*(KAPPA_C-3./2.))
           DELTA = (1-GAMMA(KAPPA_C-1./2.)/GAMMA(KAPPA_C+1./2.)*(KAPPA_C-1./2.)*(1-R1))**(-1./(KAPPA_C-1./2.))
     
           RO = SQRT(DELTA-1)/BETA
           VEL(I) = RO*SIN(TETA)
        END DO
     
        ! Step 3.
     
        VX = UX + VEL(1)
        VY = UY + VEL(2)
        VZ = UZ + VEL(3)
     
        RETURN
     
        END SUBROUTINE SAMPLE_VELOCITY_KAPPA


        FUNCTION FLX_KAPPA(THIS, SN, TINF, M) RESULT(OUT)

            IMPLICIT NONE
    
            CLASS(KAPPA_VDF), INTENT(IN) :: THIS
            REAL(KIND=8), INTENT(IN) :: SN,TINF,M
            REAL(KIND=8) :: OUT
            REAL(KIND=8) :: R1,R2
            REAL(KIND=8) :: y,fM,BETA, KAPPA, ACCA, KAPPA_C

            KAPPA_C = THIS%KAPPA
    
            BETA = 1./SQRT(2.*KB/M*TINF*(KAPPA_C-3./2.))
    
            ACCA = SQRT(SN**2+2.)                              ! Tmp variable
            KAPPA = 2./(SN+ACCA) * EXP(0.5 + 0.5*SN*(SN-ACCA)) ! variable
            !!! TODO: THIS SHOULD BE CHANGED TO POWER LAW AS WELL PROBABLY
            
    
            ! Step 1.
            DO
            R1 = rf()
            y  = -3.+6.*R1
    
            ! Step 2.
    
            R2 = rf()
            fM = KAPPA*(y+sn)*GAMMA(KAPPA_C)/GAMMA(KAPPA_C-1./2.)/(1+y**2)**KAPPA_C
    
            ! Step 3. 
    
            IF (R2 .LE. fM) THEN
                OUT = y/BETA
                EXIT
            END IF
            END DO
    
            RETURN
    
        END FUNCTION FLX_KAPPA
    
        FUNCTION FLUXSOURCE_KAPPA(THIS,S_NORM, U_NORM, TTRA, M) RESULT(OUT)
            CLASS(KAPPA_VDF), INTENT(IN) :: THIS
            REAL(KIND=8), INTENT(INOUT) :: S_NORM
            REAL(KIND=8), INTENT(IN) :: U_NORM,TTRA,M
            REAL(KIND=8) :: OUT
    
            REAL(KIND=8) :: BETA, KAPPA_C
            REAL(KIND=8) :: INTEG1, INTEG2
            INTEGER :: I

            KAPPA_C = THIS%KAPPA
    
            BETA = 1./SQRT(2.*KB/M*TTRA*(KAPPA_C-3./2.))
            S_NORM = U_NORM*BETA

            INTEG1 = PI/2.
            INTEG2 = ATAN(S_NORM)
            DO I = 2, INT(KAPPA_C)
                INTEG1 = (2.*I-3.)/(2.*(I-1.))*INTEG1
                INTEG2 = S_NORM/(2.*(I-1.)*(1.+(S_NORM)**2)**(I-1.)) + (2.*I-3.)/(2.*(I-1.))*INTEG2
            END DO
            OUT = GAMMA(KAPPA_C)/GAMMA(KAPPA_C-1./2.)/(BETA*SQRT(PI)) * &
                ( 1./(KAPPA_C-1.)/2.*(1+S_NORM**2)**(-KAPPA_C+1.) + S_NORM*(INTEG1 +INTEG2))

        END FUNCTION FLUXSOURCE_KAPPA

        FUNCTION BETA_KAPPA(THIS, TTRA, M) RESULT(OUT)
            CLASS(KAPPA_VDF), INTENT(IN) :: THIS
    
            REAL(KIND=8), INTENT(IN) :: TTRA, m
            REAL(KIND=8) :: OUT

            REAL(KIND=8) :: KAPPA_C
            KAPPA_C = THIS%KAPPA
    
            OUT = 1./SQRT(2.*KB/M*TTRA*(KAPPA_C-3./2.))
    
         END FUNCTION BETA_KAPPA

END MODULE velocity_distribution