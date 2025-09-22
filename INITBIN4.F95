!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE INITBIN4
!
!**********************************************************************C
!
! M.R. MORTON    02 JUN 1999
!
! **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
!
!----------------------------------------------------------------------C
!
!
!**********************************************************************C
!
! INITIALIZES BINARY FILE FOR EFDC OUTPUT.  PLACES CONTROL
! PARAMETERS FOR POST-PROCESSOR IN HEADER SECTION OF BINARY
! FILE WQSDTS.BIN FOR BENTHIC FLUX RATES.
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'

      PARAMETER(MXPARM=30)
      REAL TEND
      DIMENSION XLON(LCM), YLAT(LCM)
      INTEGER NPARM, NCELLS
      LOGICAL FEXIST, IS1OPEN, IS2OPEN
      CHARACTER*20 WQNAME(MXPARM)
      CHARACTER*10 WQUNITS(MXPARM)
      CHARACTER*3  WQCODE(MXPARM)
!
!**********************************************************************C
!
! THE FOLLOWING PARAMETERS ARE SPECIFIED IN EFDC.INP AND WQ3DSD.INP:
! KCSD     = NUMBER OF VERTICAL LAYERS (FORCED TO 1 FOR BENTHIC FLUX FILE)
! ISMTSDT  = NUMBER OF TIME STEPS PER DATA DUMP OF BENTHIC FLUXES
! DT       = TIME STEP OF EFDC MODEL IN SECONDS
! LA       = NUMBER OF ACTIVE CELLS + 1 IN MODEL
! TBEGAN   = BEGINNING TIME OF RUN IN DAYS
!
! THE PARAMETER NPARM MUST BE CHANGED IF THE OUTPUT DATA
! IS CHANGED IN SUBROUTINE WSMTSBIN:
! NPARM   = NUMBER OF PARAMETERS WRITTEN TO BINARY FILE
!
! NREC4   = NUMBER OF RECORDS WRITTEN TO BINARY FILE (ONE RECORD
!           IS A COMPLETE DATA DUMP FOR TIME INTERVAL IWQDIUDT)
!
      NPARM = 8
      NCELLS = LA-1
      NREC4 = 0
      TEND = TBEGIN
      KCSD = 1
      MAXRECL4 = 32
      IF(NPARM .GE. 8)THEN
        MAXRECL4 = NPARM*4
      ENDIF
!
! THE FOLLOWING WATER QUALITY NAMES, UNITS, AND 3-CHARACTER CODES
! SHOULD BE MODIFIED TO MATCH THE PARAMETERS WRITTEN TO THE BINARY
! FILE IN SUBROUTINE WSMTSBIN.  THE CHARACTER STRINGS MUST BE
! EXACTLY THE LENGTH SPECIFIED BELOW IN ORDER FOR THE POST-PROCESSOR
! TO WORK CORRECTLY.
!
! BE SURE WQNAME STRINGS ARE EXACTLY 20-CHARACTERS LONG:
!------------------'         1         2'
!------------------'12345678901234567890'
      WQNAME( 1) = 'SOD_BENTHIC_FLUX    '
      WQNAME( 2) = 'NH4_BENTHIC_FLUX    '
      WQNAME( 3) = 'NO3_BENTHIC_FLUX    '
      WQNAME( 4) = 'PO4D_BENTHIC_FLUX   '
      WQNAME( 5) = 'SAD_BENTHIC_FLUX    '
      WQNAME( 6) = 'COD_BENTHIC_FLUX    '
      WQNAME( 7) = 'SEDIMENT_TEMPERATURE'
      WQNAME( 8) = 'BENTHIC_STRESS      '
!
! BE SURE WQUNITS STRINGS ARE EXACTLY 10-CHARACTERS LONG:
!-------------------'         1'
!-------------------'1234567890'
      WQUNITS( 1) = 'G/M2/DAY  '
      WQUNITS( 2) = 'G/M2/DAY  '
      WQUNITS( 3) = 'G/M2/DAY  '
      WQUNITS( 4) = 'G/M2/DAY  '
      WQUNITS( 5) = 'G/M2/DAY  '
      WQUNITS( 6) = 'G/M2/DAY  '
      WQUNITS( 7) = 'DEGC      '
      WQUNITS( 8) = 'DAYS      '
!
! BE SURE WQCODE STRINGS ARE EXACTLY 3-CHARACTERS LONG:
!
!------------------'123'
      WQCODE( 1) = 'SOD'
      WQCODE( 2) = 'FNH'
      WQCODE( 3) = 'FNO'
      WQCODE( 4) = 'FP4'
      WQCODE( 5) = 'FSA'
      WQCODE( 6) = 'FCO'
      WQCODE( 7) = 'SMT'
      WQCODE( 8) = 'BST'
!
!---------------------------------------------------------
!
! IF WQSDTS.BIN ALREADY EXISTS, OPEN FOR APPENDING HERE.
!
      IF(ISSDBIN .EQ. 2)THEN
        IO = 1
5       IO = IO+1
        IF(IO .GT. 99)THEN
          WRITE(0,*) ' NO AVAILABLE IO UNITS ... IO > 99'
          STOP ' EFDC HALTED IN SUBROUTINE INITBIN4'
        ENDIF
        INQUIRE(UNIT=IO, OPENED=IS2OPEN)
        IF(IS2OPEN) GOTO 5
        INQUIRE(FILE='WQSDTS.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=IO, FILE='WQSDTS.BIN', ACCESS='DIRECT',FORM='UNFORMATTED', STATUS='UNKNOWN', &
               RECL=MAXRECL4)
          WRITE(0,*) 'OLD FILE WQSDTS.BIN FOUND...OPENING FOR APPEND'
          READ(IO, REC=1) NREC4, TBEGAN, TEND, DT, ISMTSDT, NPARM,NCELLS, KCSD
          NR6 = 1 + NPARM*3 + NCELLS*4 + (NCELLS*KCSD+1)*NREC4 + 1
          CLOSE(IO)
        ELSE
          ISSDBIN=1
        ENDIF
      ENDIF
!
!-------------------------------------------------------------------
!
! IF WQSDTS.BIN ALREADY EXISTS, DELETE IT HERE.
!
      IF(ISSDBIN .EQ. 1)THEN
        TBEGAN = TBEGIN
        IO = 1
10      IO = IO+1
        IF(IO .GT. 99)THEN
          WRITE(0,*) ' NO AVAILABLE IO UNITS ... IO > 99'
          STOP ' EFDC HALTED IN SUBROUTINE INITBIN4'
        ENDIF
        INQUIRE(UNIT=IO, OPENED=IS2OPEN)
        IF(IS2OPEN) GOTO 10
        INQUIRE(FILE='WQSDTS.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=IO, FILE='WQSDTS.BIN')
          CLOSE(UNIT=IO, STATUS='DELETE')
          WRITE(0,*) 'OLD FILE WQSDTS.BIN DELETED...'
        ENDIF

        OPEN(UNIT=IO, FILE='WQSDTS.BIN', ACCESS='DIRECT',FORM='UNFORMATTED', STATUS='UNKNOWN', &
             RECL=MAXRECL4)
!--------------------------------------------------------------------
! WRITE CONTROL PARAMETERS FOR POST-PROCESSOR TO HEADER
! SECTION OF THE WQSDTS.BIN BINARY FILE:
!
        WRITE(IO) NREC4, TBEGAN, TEND, DT, ISMTSDT, NPARM, NCELLS, KCSD
        DO I=1,NPARM
          WRITE(IO) WQNAME(I)
        ENDDO
        DO I=1,NPARM
          WRITE(IO) WQUNITS(I)
        ENDDO
        DO I=1,NPARM
          WRITE(IO) WQCODE(I)
        ENDDO
!
! WRITE CELL I,J MAPPING REFERENCE TO HEADER SECTION OF BINARY FILE:
!
        DO L=2,LA
          WRITE(IO) IL(L)
        ENDDO
        DO L=2,LA
          WRITE(IO) JL(L)
        ENDDO
!
! **  READ IN XLON AND YLAT OR UTME AND UTMN OF CELL CENTERS OF
! **  CURVILINEAR PORTION OF THE GRID FROM FILE LXLY.INP:
!
        IO1 = 0
20      IO1 = IO1+1
        IF(IO1 .GT. 99)THEN
          WRITE(0,*) ' NO AVAILABLE IO UNITS ... IO1 > 99'
          STOP ' EFDC HALTED IN SUBROUTINE INITBIN4'
        ENDIF
        INQUIRE(UNIT=IO1, OPENED=IS1OPEN)
        IF(IS1OPEN) GOTO 20
        OPEN(IO1,FILE='LXLY.INP',STATUS='UNKNOWN')
!
        DO NS=1,4
          READ(IO1,1111)
        ENDDO
 1111   FORMAT(80X)
!
        DO LL=1,LVC
          READ(IO1,*) I,J,XUTME,YUTMN
          L=LIJ(I,J)
          XLON(L)=XUTME
          YLAT(L)=YUTMN
        ENDDO
        CLOSE(IO1)
!
! WRITE XLON AND YLAT OF CELL CENTERS TO HEADER SECTION OF
! BINARY OUTPUT FILE:
!
        DO L=2,LA
          WRITE(IO) XLON(L)
        ENDDO
        DO L=2,LA
          WRITE(IO) YLAT(L)
        ENDDO

        INQUIRE(UNIT=IO, NEXTREC=NR6)
        CLOSE(IO)
      ENDIF

      RETURN
      END