!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE INITBIN3
!
!**********************************************************************C
!
! M.R. MORTON    29 APR 1999
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
! FILE WQDOCOMP.BIN FOR D.O. COMPONENT ANALYSIS.
!
!**********************************************************************C
!
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'

      PARAMETER(MXPARM=33)
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
! THE FOLLOWING PARAMETERS ARE SPECIFIED IN EFDC.INP AND WQ3DWC.INP:
! KC       = NUMBER OF VERTICAL LAYERS
! IWQTSDT  = NUMBER OF TIME STEPS PER DATA DUMP
! DT       = TIME STEP OF EFDC MODEL IN SECONDS
! LA       = NUMBER OF ACTIVE CELLS + 1 IN MODEL
! TBEGAN   = BEGINNING TIME OF RUN IN DAYS
!
! THE PARAMETER NPARM MUST BE CHANGED IF THE OUTPUT DATA
! IS CHANGED IN SUBROUTINE WWQTSBIN:
! NPARM   = NUMBER OF PARAMETERS WRITTEN TO BINARY FILE
!
! NREC3   = NUMBER OF RECORDS WRITTEN TO BINARY FILE (ONE RECORD
!           IS A COMPLETE DATA DUMP FOR TIME INTERVAL IWQDIUDT)
!
      NPARM = 33
      NCELLS = LA-1
      NREC3 = 0
      TEND = TBEGIN
      MAXRECL3 = 32
      IF(NPARM .GE. 8)THEN
        MAXRECL3 = NPARM*4
      ENDIF
!
! THE FOLLOWING WATER QUALITY NAMES, UNITS, AND 3-CHARACTER CODES
! SHOULD BE MODIFIED TO MATCH THE PARAMETERS WRITTEN TO THE BINARY
! FILE IN SUBROUTINE WWQTSBIN.  THE CHARACTER STRINGS MUST BE
! EXACTLY THE LENGTH SPECIFIED BELOW IN ORDER FOR THE POST-PROCESSOR
! TO WORK CORRECTLY.
!
! BE SURE WQNAME STRINGS ARE EXACTLY 20-CHARACTERS LONG:
!------------------'         1         2'
!------------------'12345678901234567890'
      WQNAME( 1) = 'NITROGEN_LIMIT_CYA  '
      WQNAME( 2) = 'NITROGEN_LIMIT_DIA  '
      WQNAME( 3) = 'NITROGEN_LIMIT_GRN  '
      WQNAME( 4) = 'NITROGEN_LIMIT_MAC  '
      WQNAME( 5) = 'PHOSPHORUS_LIMIT_CYA'
      WQNAME( 6) = 'PHOSPHORUS_LIMIT_DIA'
      WQNAME( 7) = 'PHOSPHORUS_LIMIT_GRN'
      WQNAME( 8) = 'PHOSPHORUS_LIMIT_MAC'
      WQNAME( 9) = 'LIGHT_LIMIT_CYA     '
      WQNAME(10) = 'LIGHT_LIMIT_DIA     '
      WQNAME(11) = 'LIGHT_LIMIT_GRN     '
      WQNAME(12) = 'LIGHT_LIMIT_MAC     '
      WQNAME(13) = 'TEMP_LIMIT_CYA      '
      WQNAME(14) = 'TEMP_LIMIT_DIA      '
      WQNAME(15) = 'TEMP_LIMIT_GRN      '
      WQNAME(16) = 'TEMP_LIMIT_MAC      '
      WQNAME(17) = 'VELOCITY_LIMIT_MAC  '
      WQNAME(18) = 'DENSITY_LIMIT_MAC   '
      WQNAME(19) = 'DO_SATURATION       '
      WQNAME(20) = 'DO_POINT_SOURCES    '
      WQNAME(21) = 'DO_SED_OXYGEN_DEMAND'
      WQNAME(22) = 'DO_REAERATION       '
      WQNAME(23) = 'DO_DOC_DECAY        '
      WQNAME(24) = 'DO_NH4_NITRIFICATION'
      WQNAME(25) = 'DO_COD_OXIDATION    '
      WQNAME(26) = 'DO_PHOTOSYNTH_CHL   '
      WQNAME(27) = 'DO_RESPIRATION_CHL  '
      WQNAME(28) = 'DO_PHOTOSYNTH_MAC   '
      WQNAME(29) = 'DO_RESPIRATION_MAC  '
      WQNAME(30) = 'DO_DEFICIT          '
      WQNAME(31) = 'DO_TRANSPORT        '
      WQNAME(32) = 'DO_ALL_COMPONENTS   '
      WQNAME(33) = 'LAYER_THICKNESS     '
!
! BE SURE WQUNITS STRINGS ARE EXACTLY 10-CHARACTERS LONG:
!-------------------'         1'
!-------------------'1234567890'
      WQUNITS( 1) = 'UNITLESS  '
      WQUNITS( 2) = 'UNITLESS  '
      WQUNITS( 3) = 'UNITLESS  '
      WQUNITS( 4) = 'UNITLESS  '
      WQUNITS( 5) = 'UNITLESS  '
      WQUNITS( 6) = 'UNITLESS  '
      WQUNITS( 7) = 'UNITLESS  '
      WQUNITS( 8) = 'UNITLESS  '
      WQUNITS( 9) = 'UNITLESS  '
      WQUNITS(10) = 'UNITLESS  '
      WQUNITS(11) = 'UNITLESS  '
      WQUNITS(12) = 'UNITLESS  '
      WQUNITS(13) = 'UNITLESS  '
      WQUNITS(14) = 'UNITLESS  '
      WQUNITS(15) = 'UNITLESS  '
      WQUNITS(16) = 'UNITLESS  '
      WQUNITS(17) = 'UNITLESS  '
      WQUNITS(18) = 'UNITLESS  '
      WQUNITS(19) = 'MG/L/DAY  '
      WQUNITS(20) = 'MG/L/DAY  '
      WQUNITS(21) = 'MG/L/DAY  '
      WQUNITS(22) = 'MG/L/DAY  '
      WQUNITS(23) = 'MG/L/DAY  '
      WQUNITS(24) = 'MG/L/DAY  '
      WQUNITS(25) = 'MG/L/DAY  '
      WQUNITS(26) = 'MG/L/DAY  '
      WQUNITS(27) = 'MG/L/DAY  '
      WQUNITS(28) = 'MG/L/DAY  '
      WQUNITS(29) = 'MG/L/DAY  '
      WQUNITS(30) = 'MG/L/DAY  '
      WQUNITS(31) = 'MG/L/DAY  '
      WQUNITS(32) = 'MG/L/DAY  '
      WQUNITS(33) = 'METERS    '
!
! BE SURE WQCODE STRINGS ARE EXACTLY 3-CHARACTERS LONG:
!
!------------------'123'
      WQCODE( 1) = 'NLC'
      WQCODE( 2) = 'NLD'
      WQCODE( 3) = 'NLG'
      WQCODE( 4) = 'NLM'
      WQCODE( 5) = 'PLC'
      WQCODE( 6) = 'PLD'
      WQCODE( 7) = 'PLG'
      WQCODE( 8) = 'PLM'
      WQCODE( 9) = 'LLC'
      WQCODE(10) = 'LLD'
      WQCODE(11) = 'LLG'
      WQCODE(12) = 'LLM'
      WQCODE(13) = 'TLC'
      WQCODE(14) = 'TLD'
      WQCODE(15) = 'TLG'
      WQCODE(16) = 'TLM'
      WQCODE(17) = 'VLM'
      WQCODE(18) = 'DLM'
      WQCODE(19) = 'DCS'
      WQCODE(20) = 'DPS'
      WQCODE(21) = 'DSO'
      WQCODE(22) = 'DKA'
      WQCODE(23) = 'DCA'
      WQCODE(24) = 'DNH'
      WQCODE(25) = 'DCO'
      WQCODE(26) = 'DPC'
      WQCODE(27) = 'DRC'
      WQCODE(28) = 'DPM'
      WQCODE(29) = 'DRM'
      WQCODE(30) = 'DEF'
      WQCODE(31) = 'DTR'
      WQCODE(32) = 'DAL'
      WQCODE(33) = 'DZZ'
!
!---------------------------------------------------------
!
! IF WQDOCOMP.BIN ALREADY EXISTS, OPEN FOR APPENDING HERE.
!
      IF(ISCOMP .EQ. 2)THEN
        IO = 1
5       IO = IO+1
        IF(IO .GT. 99)THEN
          WRITE(0,*) ' NO AVAILABLE IO UNITS ... IO > 99'
          STOP ' EFDC HALTED IN SUBROUTINE INITBIN3'
        ENDIF
        INQUIRE(UNIT=IO, OPENED=IS2OPEN)
        IF(IS2OPEN) GOTO 5
        INQUIRE(FILE='WQDOCOMP.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=IO, FILE='WQDOCOMP.BIN', ACCESS='DIRECT',FORM='UNFORMATTED', STATUS='UNKNOWN',&
               RECL=MAXRECL3)
          WRITE(0,*) 'OLD FILE WQDOCOMP.BIN FOUND...OPENING FOR APPEND'
          READ(IO, REC=1) NREC3, TBEGAN, TEND, DT, IWQTSDT, NPARM,NCELLS, KC
          NR5 = 1 + NPARM*3 + NCELLS*4 + (NCELLS*KC+1)*NREC3 + 1
          CLOSE(IO)
        ELSE
          ISCOMP=1
        ENDIF
      ENDIF
!
!-------------------------------------------------------------------
!
! IF WQDOCOMP.BIN ALREADY EXISTS, DELETE IT HERE.
!
      IF(ISCOMP .EQ. 1)THEN
        TBEGAN = TBEGIN
        IO = 1
10      IO = IO+1
        IF(IO .GT. 99)THEN
          WRITE(0,*) ' NO AVAILABLE IO UNITS ... IO > 99'
          STOP ' EFDC HALTED IN SUBROUTINE INITBIN3'
        ENDIF
        INQUIRE(UNIT=IO, OPENED=IS2OPEN)
        IF(IS2OPEN) GOTO 10
        INQUIRE(FILE='WQDOCOMP.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=IO, FILE='WQDOCOMP.BIN')
          CLOSE(UNIT=IO, STATUS='DELETE')
          WRITE(0,*) 'OLD FILE WQDOCOMP.BIN DELETED...'
        ENDIF

        OPEN(UNIT=IO, FILE='WQDOCOMP.BIN', ACCESS='DIRECT',FORM='UNFORMATTED', STATUS='UNKNOWN',  &
             RECL=MAXRECL3)
!--------------------------------------------------------------------
! WRITE CONTROL PARAMETERS FOR POST-PROCESSOR TO HEADER
! SECTION OF THE WQDOCOMP.BIN BINARY FILE:
!
        WRITE(IO) NREC3, TBEGAN, TEND, DT, IWQTSDT, NPARM, NCELLS, KC
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
          STOP ' EFDC HALTED IN SUBROUTINE INITBIN3'
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

        INQUIRE(UNIT=IO, NEXTREC=NR5)
        CLOSE(IO)
      ENDIF

      RETURN
      END