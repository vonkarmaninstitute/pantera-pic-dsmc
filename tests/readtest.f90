program read

      INTEGER, PARAMETER :: in1 = 444
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER*512      :: MIXTURE_DEFINITION, VSS_PARAMS_FILENAME, LINESOURCE_DEFINITION
      CHARACTER*64       :: MIX_INIT_NAME, MIX_BOUNDINJECT_NAME, DSMC_COLL_MIX_NAME

      CHARACTER*512 :: filename

      ! Open input file for reading
      OPEN(UNIT=in1,FILE='input', STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         PRINT*
         WRITE(*,*)'  Attention, "input" file not found! ABORTING.'
         PRINT*
         STOP
      ENDIF


      line = '' ! Init empty
      filename = ''

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in1,'(A)', IOSTAT=ReasonEOF) line ! Read line

         PRINT*, "Current line: ", line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         IF (line=='READ_HERE:')     READ(in1,'(A)') FILENAME
         ! IF (line=='Number_of_cells:')         READ(in1,*) NX, NY, NZ

       END DO

       PRINT*, filename


end program read

