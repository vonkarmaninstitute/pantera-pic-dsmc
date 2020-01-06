
PROGRAM PANTERA 

   USE particle
   USE mpi_common
   USE global
   USE screen
   USE tools
   USE initialization
   USE EM_fields
   USE timecycle
   USE grid_and_partition

   IMPLICIT NONE

   ! ========= Init MPI environment ========================
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_THREADS, ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROC_ID, ierr)
   CALL NEWTYPE ! Define new "particle" type for mpi

   ! ========= Print header (important things first) =======
   CALL PRINTTITLE(PROC_ID)

   ! ========= Read input file and init variables ==========
   CALL ONLYMASTERPRINT1(PROC_ID, '> READING INPUT DATA...')
   CALL READINPUT          ! Read input file
   CALL INITVARIOUS        ! Initialize some additional variables
   CALL INITINJECTION      ! Initialize variables for injection
   CALL INITCOLLISIONS     ! Initialize variables for collisions

   ! ========= Initial particles seed ======================
   IF (BOOL_INITIAL_SEED)   CALL INITIAL_SEED
   ! CALL DUMP_PARTICLES_SCREEN

   ! ========= Time loop ===================================
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
   CALL ONLYMASTERPRINT1(PROC_ID, '> STARTING TIME LOOP...')
   CALL TIME_LOOP

   ! CALL DUMP_PARTICLES_SCREEN

   ! ========== Close MPI environment ======================
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

   CALL ONLYMASTERPRINT1(PROC_ID, '> R O A R !')
   CALL MPI_FINALIZE(ierr)


!   ! Inizializza simulazione e dominio
! 
!   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!   stringa = '0) READ INPUT DATA...'
!   CALL ONLYMASTERPRINT1(PE_ID, vuoto)
!   CALL ONLYMASTERPRINT1(PE_ID, stringa)
!   CALL LETTURA
! 
!   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!   stringa = '2) DSMC INITIALIZATION...'
!   CALL ONLYMASTERPRINT1(PE_ID, vuoto)
!   CALL ONLYMASTERPRINT1(PE_ID, stringa)
!   CALL INIZIO
! 
!   ! Ciclo temporale
! 
!   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!   stringa = '3) TIME CYCLE...'
!   CALL ONLYMASTERPRINT1(PE_ID, vuoto)
!   CALL ONLYMASTERPRINT1(PE_ID, stringa) 
!   CALL ONLYMASTERPRINT1(PE_ID, vuoto)
! 
!   CALL CPU_TIME(tin)
!   CALL CICLO
!   CALL CPU_TIME(tfin)
! 
!   ! Calcola e stampa informazioni
! 
!   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!   CALL ONLYMASTERPRINT1(PE_ID, vuoto)
!   stringa = '4) STATISTICS OF THE RUN...'
!   CALL ONLYMASTERPRINT1(PE_ID, stringa)
!   CALL STATISTICHE
! 
!   ! Libera la memoria
! 
!   CALL DEALLOCA2
!   CLOSE(PE_ID)
! 
!   ! Finalizzazione ambiente parallelo
! 
!   stringa = ' Total time of the simulation [s]:'
!   CALL ONLYMASTERPRINT1(PE_ID, vuoto)
!   CALL ONLYMASTERPRINT2(PE_ID, stringa, tfin-tin)  
!   CALL ONLYMASTERPRINT1(PE_ID, vuoto) 
!   stringa = '   L E S   J E U X   S O N T   F A I T E S !!!'
!   CALL ONLYMASTERPRINT1(PE_ID, stringa) 
!   CALL ONLYMASTERPRINT1(PE_ID, vuoto) 
!   CALL MPI_FINALIZE(ierr)
! 
!   CONTAINS
! 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ! SUBROUTINE NEWTYPE -> Definisce un nuovo tipo necessario ad MPI !!!!
!   ! per impacchettare i messaggi nel formato "particella" !!!!!!!!!!!!!!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   SUBROUTINE NEWTYPE
! 
!   INTEGER ii
! 
!   CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent_dpr, ierr) 
!   CALL MPI_TYPE_EXTENT(MPI_INTEGER,          extent_int, ierr)
!       
!   blocklengths = 1
! 
!   oldtypes(1:7) = MPI_DOUBLE_PRECISION 
!   oldtypes(8:9) = MPI_INTEGER
!      
!   offsets(1) = 0 
!   DO ii = 2, 8
!      offsets(ii) = offsets(ii - 1) + extent_dpr * blocklengths(ii - 1)
!   END DO
!   offsets(9) = offsets(8) + extent_int * blocklengths(8)
!  
!   CALL MPI_TYPE_STRUCT(9, blocklengths, offsets, oldtypes, MPI_PARTICLE_DATA_STRUCTURE, ierr) 
!   CALL MPI_TYPE_COMMIT(MPI_PARTICLE_DATA_STRUCTURE, ierr)  
! 
!   END SUBROUTINE NEWTYPE
! 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ! SUBROUTINE DEALLOCA2 -> Dealloca memoria per simulazione !!!!!!!!!!!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   SUBROUTINE DEALLOCA2 
! 
!   DEALLOCATE(particle)
! 
!   DEALLOCATE(NPC) 
!   DEALLOCATE(IOF) 
! 
!   DEALLOCATE(DMED)   
!   DEALLOCATE(VXMED) 
!   DEALLOCATE(VYMED) 
!   DEALLOCATE(VZMED)  
!   DEALLOCATE(TMED)  
!   DEALLOCATE(TRMED) 
!   DEALLOCATE(TTMED)  
!   DEALLOCATE(AVCORA)
! 
!   DEALLOCATE(DMED1)  
!   DEALLOCATE(NMED1)   
!   DEALLOCATE(VXMED1)  
!   DEALLOCATE(VYMED1) 
!   DEALLOCATE(VZMED1)  
!   DEALLOCATE(TMED1)   
!   DEALLOCATE(TRMED1) 
!   DEALLOCATE(TTMED1) 
!   DEALLOCATE(TXMED1)
!   DEALLOCATE(TYMED1)  
!   DEALLOCATE(TZMED1) 
! 
!   DEALLOCATE(DMEAN)
!   DEALLOCATE(NMEAN)
!   DEALLOCATE(TMEAN)
!   DEALLOCATE(TTMEAN)
!   DEALLOCATE(TRMEAN)
!   DEALLOCATE(VXMEAN)
!   DEALLOCATE(VYMEAN)
!   DEALLOCATE(VZMEAN)
! 
!   DEALLOCATE(CELL_PARTITION) 
! 
!   DEALLOCATE(sendcount)
!   DEALLOCATE(recvcount)
!   DEALLOCATE(senddispl)
!   DEALLOCATE(recvdispl)
! 
!   DEALLOCATE(XNODES)
!   DEALLOCATE(YNODES)
! 
!   DEALLOCATE(DM_LIM1)
!   DEALLOCATE(DM_LIM2)
! 
!   END SUBROUTINE DEALLOCA2

  END PROGRAM PANTERA
