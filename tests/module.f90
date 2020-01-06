PROGRAM MOD90

  IMPLICIT NONE

  REAL(KIND=8) :: xp, yp, Lx, Ly, xMin, yMin, xMax, yMax

  xMin = -0.9221;
  xMax = 0.823;

  yMin = 0.1;
  yMax = 2.2;

  Lx = xMax - xMin;
  Ly = yMax - yMin;

  PRINT*, xMin + MOD(-1.031 - xMin, xMax - xMin)

END PROGRAM
