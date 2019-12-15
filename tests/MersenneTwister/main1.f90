PROGRAM test

  USE mt19937_64

  IMPLICIT NONE

  INTEGER :: II
  INTEGER(KIND=8) :: seed_MT
  INTEGER(KIND=4) :: seed_rand, seed_BASE, seed_RANDOM
  REAL(KIND=8)    :: rr

  ! ============ Select global seed ========
  seed_BASE = 1234

  ! ============ First, use MT ===========
  seed_MT = seed_BASE
  CALL init_genrand64(seed_MT)
  DO II = 1, 500
    PRINT*, "DBDB MT rand1: ",  genrand64_real1()
  END DO

  seed_MT = 2*seed_BASE
  CALL init_genrand64(seed_MT)
  DO II = 1, 500
    PRINT*, "DBDB MT rand2: ",  genrand64_real1()
  END DO

  seed_MT = 3*seed_BASE
  CALL init_genrand64(seed_MT)
  DO II = 1, 500
    PRINT*, "DBDB MT rand3: ",  genrand64_real1()
  END DO

  ! ============ Now a call to rand() ===========
  seed_rand = seed_BASE
PRINT*, seed_rand
CALL SLEEP(1)
  CALL srand(seed_rand)
PRINT*, seed_rand
CALL SLEEP(1)
  DO II = 1, 500
    PRINT*, "DBDB other1: ",  rand()
  END DO

  seed_rand = 2*seed_BASE
PRINT*, seed_rand
CALL SLEEP(1)
  CALL srand(seed_rand)
PRINT*, seed_rand
CALL SLEEP(1)
  DO II = 1, 500
    PRINT*, "DBDB other2: ",  rand()
  END DO

  seed_rand = 3*seed_BASE
  CALL srand(seed_rand)
  DO II = 1, 500
    PRINT*, "DBDB other3: ",  rand()
  END DO

  ! =========== Now try RANDOM_NUMBER() ==============
  seed_RANDOM = seed_BASE
  CALL RANDOM_SEED(seed_RANDOM)
  DO II = 1, 500
    CALL RANDOM_NUMBER(rr)
    PRINT*, "DBDB RANDOM_NUMBER1: ",  rr
  END DO

  seed_RANDOM = 2*seed_BASE
  CALL RANDOM_SEED(seed_RANDOM)
  DO II = 1, 500
    CALL RANDOM_NUMBER(rr)
    PRINT*, "DBDB RANDOM_NUMBER2: ",  rr
  END DO

  seed_RANDOM = 3*seed_BASE
  CALL RANDOM_SEED(seed_RANDOM)
  DO II = 1, 500
    CALL RANDOM_NUMBER(rr)
    PRINT*, "DBDB RANDOM_NUMBER3: ",  rr
  END DO

END PROGRAM test
