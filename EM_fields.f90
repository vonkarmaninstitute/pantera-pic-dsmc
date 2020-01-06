! This module contains routines for elecromagnetic-fields computation, sampling etc.

MODULE EM_fields

   USE global
   USE mpi_common

   IMPLICIT NONE

   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBORUTINE GET_E_B_POSITION_1D -> returns the electric and magnetic fields at a position !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GET_E_B_POSITION_1D(x,  E, B)

      REAL(KIND=8), INTENT(IN)  ::  x
      REAL(KIND=8), INTENT(OUT) ::  E, B

      ! +++++ Magnetic field +++++
 
      IF (B_FIELD_BOOL .EQV. .FALSE.) THEN ! No B field for the simulation: put it to zero
         B = 0
      ELSEIF (B_FIELD_TYPE .EQ. "UNIFORM") THEN ! Uniform field
         B = B_UNIFORM_VAL
      ELSEIF (B_FIELD_TYPE .EQ. "SINE") THEN ! Compute using sin() function -> FIND A BETTER WAY TO DO THIS! NO CHAR COMPARISON
         B = B_SIN_AMPL * SIN( B_SIN_k*x - B_SIN_omega*CURRENT_TIME + B_SIN_PHASE_rad )
      ENDIF
 
      ! +++++ Electric field +++++
 
      IF (E_FIELD_BOOL .EQV. .FALSE.) THEN ! No E field for the simulation: put it to zero
         E = 0
      ELSEIF (E_FIELD_TYPE .EQ. "UNIFORM") THEN ! Uniform field
         E = E_UNIFORM_VAL
      ELSEIF (E_FIELD_TYPE .EQ. "SINE") THEN ! Compute using sin() function -> FIND A BETTER WAY TO DO THIS! NO CHAR COMPARISON
         E = E_SIN_AMPL * SIN( E_SIN_k*x - E_SIN_omega*CURRENT_TIME + E_SIN_PHASE_rad )
      ENDIF     

   END SUBROUTINE GET_E_B_POSITION_1D

END MODULE EM_fields

