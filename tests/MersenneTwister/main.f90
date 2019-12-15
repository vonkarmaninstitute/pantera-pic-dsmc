PROGRAM test

  USE mt19937_64

  INTEGER(KIND=8) :: seed
  REAL(KIND=8)    :: randnum

  ! Initialize generator with seed
  seed = 1234
  CALL init_genrand64(seed)

  ! Find some numbers
  PRINT*, genrand64_real1()
  PRINT*, genrand64_real1()
  PRINT*, genrand64_real1()
  PRINT*, genrand64_real1()



END PROGRAM test
