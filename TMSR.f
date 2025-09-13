!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE TMSR
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 02/19/2002        john hamrick       02/19/2002       john hamrick
!  changed tox bed output
!
! 4/10/06           Hugo N Rodriguez   increase output resolution
! 12/12/2020        DRBC LZ - ADD OUTPUT OF TURBULENCE DISSIPATION RATE
!
!----------------------------------------------------------------------C
!
! **  SUBROUTINE TMSR WRITES TIME SERIES FILES FOR SURFACE ELEVATON,
! **  VELOCITY, CONCENTRATION, AND VOLUME SOURCES AT SPECIFIED
! **  (I,J) POINTS
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
!     DIMENSION ATMP(KCM),BTMP(KCM)
!
      DIMENSION SEDSND(KCM),TXWF(KCM),TXWC(KCM),TXWP(KCM),TXBT(KBM),SEDSNDB(KBM),TXBF(KBM),TXBC(KBM),TXBP(KBM),PORH(KBM)
	DIMENSION QCHANUIJ(LCM),QCHANVIJ(LCM)
!
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7,            !hnr
     &         TITLE11,TITLE12,TITLE13,TITLE14,TITLE15,TITLE16,TITLE17,
     &         TITLE18,TITLE19,TFNTXWT,TFNTXWF,TFNTXWC,TFNTXWP,
     &         TFNTXBT,TFNTXBF,TFNTXBC,TFNTXBP,TITLE20                          
      CHARACTER*10 CTUNIT                                                           !hnr
      CHARACTER*1  CZTT1(0:9)                                                       !hnr
      CHARACTER*2  CZTT2(0:99)                                                      !hnr
      CHARACTER*1  CCHTMF,CCHTML                                                    !hnr
      CHARACTER*2  CCHTMS                                                           !hnr
      CHARACTER*2  CNTOX(25),CNSND(6),CNSBL(6)                                      !hnr
!
!**********************************************************************C
!
      IF(JSTMSR.NE.1) GOTO 300
!
!----------------------------------------------------------------------C
!
      IF(MLTMSR.GT.MLTMSRM)THEN
        WRITE (6,600)
        STOP
      ENDIF
      IF(MLTMSR.GT.999)THEN
        WRITE (6,601)
        STOP
      ENDIF
!
  600 FORMAT(' NUMBER OF TIME SER LOC, MLTMSR, EXCEEDS DIM, MLTMSRM')
  601 FORMAT(' NUMBER OF TIME SERIES LOCATIONS EXCEED 999')
!
            TAUW1=0.
            TAUW2=0.
!
      CZTT1(0)='0'
      CZTT1(1)='1'
      CZTT1(2)='2'
      CZTT1(3)='3'
      CZTT1(4)='4'
      CZTT1(5)='5'
      CZTT1(6)='6'
      CZTT1(7)='7'
      CZTT1(8)='8'
      CZTT1(9)='9'
!
      CZTT2( 0)='00'
      CZTT2( 1)='01'
      CZTT2( 2)='02'
      CZTT2( 3)='03'
      CZTT2( 4)='04'
      CZTT2( 5)='05'
      CZTT2( 6)='06'
      CZTT2( 7)='07'
      CZTT2( 8)='08'
      CZTT2( 9)='09'
      CZTT2(10)='10'
      CZTT2(11)='11'
      CZTT2(12)='12'
      CZTT2(13)='13'
      CZTT2(14)='14'
      CZTT2(15)='15'
      CZTT2(16)='16'
      CZTT2(17)='17'
      CZTT2(18)='18'
      CZTT2(19)='19'
      CZTT2(20)='20'
      CZTT2(21)='21'
      CZTT2(22)='22'
      CZTT2(23)='23'
      CZTT2(24)='24'
      CZTT2(25)='25'
      CZTT2(26)='26'
      CZTT2(27)='27'
      CZTT2(28)='28'
      CZTT2(29)='29'
      CZTT2(30)='30'
      CZTT2(31)='31'
      CZTT2(32)='32'
      CZTT2(33)='33'
      CZTT2(34)='34'
      CZTT2(35)='35'
      CZTT2(36)='36'
      CZTT2(37)='37'
      CZTT2(38)='38'
      CZTT2(39)='39'
      CZTT2(40)='40'
      CZTT2(41)='41'
      CZTT2(42)='42'
      CZTT2(43)='43'
      CZTT2(44)='44'
      CZTT2(45)='45'
      CZTT2(46)='46'
      CZTT2(47)='47'
      CZTT2(48)='48'
      CZTT2(49)='49'
      CZTT2(50)='50'
      CZTT2(51)='51'
      CZTT2(52)='52'
      CZTT2(53)='53'
      CZTT2(54)='54'
      CZTT2(55)='55'
      CZTT2(56)='56'
      CZTT2(57)='57'
      CZTT2(58)='58'
      CZTT2(59)='59'
      CZTT2(60)='60'
      CZTT2(61)='61'
      CZTT2(62)='62'
      CZTT2(63)='63'
      CZTT2(64)='64'
      CZTT2(65)='65'
      CZTT2(66)='66'
      CZTT2(67)='67'
      CZTT2(68)='68'
      CZTT2(69)='69'
      CZTT2(70)='70'
      CZTT2(71)='71'
      CZTT2(72)='72'
      CZTT2(73)='73'
      CZTT2(74)='74'
      CZTT2(75)='75'
      CZTT2(76)='76'
      CZTT2(77)='77'
      CZTT2(78)='78'
      CZTT2(79)='79'
      CZTT2(80)='80'
      CZTT2(81)='81'
      CZTT2(82)='82'
      CZTT2(83)='83'
      CZTT2(84)='84'
      CZTT2(85)='85'
      CZTT2(86)='86'
      CZTT2(87)='87'
      CZTT2(88)='88'
      CZTT2(89)='89'
      CZTT2(90)='90'
      CZTT2(91)='91'
      CZTT2(92)='92'
      CZTT2(93)='93'
      CZTT2(94)='94'
      CZTT2(95)='95'
      CZTT2(96)='96'
      CZTT2(97)='97'
      CZTT2(98)='98'
      CZTT2(99)='99'
!
!99TS      DO MLTM=1,MLTMSR                             MLTM=5   MLTM=98   MLTM=101   MLTM=549
!99TS      MSDIG=MOD(MLTM,10)                             5        8            1
!99TS      MTMP=MLTM-MSDIG                                0        90          100
!99TS      MFDIG=MTMP/10                                  0        9           10
!99TS      CCHTMF=CZTT(MFDIG)
!99TS      CCHTMS=CZTT(MSDIG)
!99TS      CNTMSR(MLTM)= CCHTMF // CCHTMS // CCHTML
!99TS      ENDDO
!
      DO MLTM=1,MLTMSR                     !      1      11     99    101    111     199    210   211
	MSDIG=MOD(MLTM,100)                  !      1      11     99      1     11      99     10    11
	MTMP=MLTM-MSDIG                      !      0       0      0    100    100     100    200   200
	MFDIG=MTMP/100                       !      0       0      0      1      1       1      2     2
      CCHTMF=CZTT1(MFDIG)                  !     '0'     '0'    '0'    '1'    '1'     '1'    '2'   '2'
      CCHTMS=CZTT2(MSDIG)                  !    '01'    '11'   '99'   '01'   '11'    '99'   '10'  '11'
      CNTMSR(MLTM)= CCHTMF // CCHTMS       !   '001'   '011'  '099'  '101'  '111'   '199'  '210' '211'
!	WRITE(8,888)MLTM,CNTMSR(MLTM)
      ENDDO
!
  888 FORMAT(I5,2X,A3)
!
      IF(TCTMSR.EQ.1.) CTUNIT='SECONDS'
      IF(TCTMSR.EQ.60.) CTUNIT='MINUTES'
      IF(TCTMSR.EQ.3600.) CTUNIT='HOURS'
      IF(TCTMSR.EQ.86400.) CTUNIT='DAYS'
!
       CNTOX( 1)= '01'
       CNTOX( 2)= '02'
       CNTOX( 3)= '03'
       CNTOX( 4)= '04'
       CNTOX( 5)= '05'
       CNTOX( 6)= '06'
       CNTOX( 7)= '07'
       CNTOX( 8)= '08'
       CNTOX( 9)= '09'
       CNTOX(10)= '10'
       CNTOX(11)= '11'
       CNTOX(12)= '12'
       CNTOX(13)= '13'
       CNTOX(14)= '14'
       CNTOX(15)= '15'
       CNTOX(16)= '16'
       CNTOX(17)= '17'
       CNTOX(18)= '18'
       CNTOX(19)= '19'
       CNTOX(20)= '20'
       CNTOX(21)= '21'
       CNTOX(22)= '22'
       CNTOX(23)= '23'
       CNTOX(24)= '24'
       CNTOX(25)= '25'
!
       CNSND( 1)= '01'
       CNSND( 2)= '02'
       CNSND( 3)= '03'
       CNSND( 4)= '04'
       CNSND( 5)= '05'
       CNSND( 6)= '06'
!
       CNSBL( 1)= '01'
       CNSBL( 2)= '02'
       CNSBL( 3)= '03'
       CNSBL( 4)= '04'
       CNSBL( 5)= '05'
       CNSBL( 6)= '06'
!

! **  WRITE HEADINGS
!
      TITLE1=' SALINITY (PSU) TIME SERIES, K=1,KC'
      TITLE2=' TEMPERATURE (DEG C) TIME SERIES, K=1,KC'
      TITLE3=' DYE CONC (KG/M**3) TIME SERIES, K=1,KC'
      TITLE4=' SED CONC (MG/LITER) TIME SERIES, K=1,KC'
      TITLE5=' TOXIC CONC (M/TOT VOL - UG/LITER) 1-4 BED,5-8 WC'
      TITLE6=' VISCOSITY (CM**2/S) TIME SERIES, K=1,KS'
      TITLE7=' DIFFUSIVITY (CM**2/S) TIME SERIES, K=1,KS'  
      TITLE11=' SURFACE ELEVATION & DEPTH (METERS) TIME SERIES'
      TITLE12=' EXT MODE E,N VEL (CM/S) TBX TBY TB TBCG TBNG (CM/S)**2'
      TITLE13=' EXT MODE U,V TRANSPORT (M**3/S) TIME SERIES'
      TITLE14=' INT MODE EAST VEL (CM/S) TIME SERIES, K=1,KC'
      TITLE15=' INT MODE NORTH VEL (CM/S) TIME SERIES, K=1,KC'
      TITLE16=' EXT MODE VOLUME S/S (M**3/S) TIME SERIES'
      TITLE17=' INT MODE VOL S/S (M**3/S) TIME SERIES, K=1,KC'
      TITLE18=' SED BED LOAD QSX QSY (GM/S) CQSX CQSY (MG/L) '
      TITLE19=' BED TOP  KBT HBED(KBT) HBED(KBT-1) VOIDR FRAC SED/SND'
      TITLE20=' DISSIPATION (W/KG) TIME SERIES, K=1,KC'          
      TFNTXWT=' TOTAL TOXIC CONC WATER COL (M/TOT VOL), UG/LITER'
      TFNTXWF=' FREE DIS TOXIC CONC WATER COL (M/TOT VOL), UG/LITER'
      TFNTXWC=' DOC COMP TOXIC CONC WATER COL (M/TOT VOL), UG/LITER'
      TFNTXWP=' TOT PART TOXIC CONC WATER COL (M/M), UG/GM'
      TFNTXBT=' TOTAL TOXIC CONC SED BED (M/TOT VOL), UG/LITER'
      TFNTXBF=' FREE DIS TOXIC CONC SED BED (M/PORE VOL), UG/LITER'
      TFNTXBC=' DOC COMP TOXIC CONC SED BED (M/PORE VOL), UG/LITER'
      TFNTXBP=' TOT PART TOXIC CONC SED BED (M/M), UG/GM'
!
      IF(ISTMSR.EQ.2)THEN
      DO MLTM=1,MLTMSR
        IF(MTMSRC(MLTM).GE.1)THEN
          IF(ISTRAN(1).GE.1)THEN
            FNSAL(MLTM)='SALTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(2).GE.1)THEN
            FNTEM(MLTM)='TEMTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(3).GE.1)THEN
            FNDYE(MLTM)='DYETS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(4).GE.1)THEN
            FNSFL(MLTM)='SFLTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(6).GE.1)THEN
            FNSED(MLTM)='SEDTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
!          IF(ISTRAN(7).GE.1)THEN
!            FNSND(MLTM)='SNDTS' // CNTMSR(MLTM) // '.OUT'
!          ENDIF
          IF(ISTRAN(7).GE.1)THEN
            DO NX=1,NSND
            FNSND(MLTM,NX)='SND' // CNSND(NX) // 'TS' // CNTMSR(MLTM) // '.OUT'
            FNSBL(MLTM,NX)='SBL' // CNSBL(NX) // 'TS' // CNTMSR(MLTM) // '.OUT'
            ENDDO
          ENDIF
          IF(ISTRAN(8).GE.1)THEN
            FNDOX(MLTM)='DOXTS' // CNTMSR(MLTM) // '.OUT'
            FNTOC(MLTM)='TOCTS' // CNTMSR(MLTM) // '.OUT'
            FNNHX(MLTM)='NHXTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(5).GE.1)THEN
            DO NT=1,NTOX
            FNTXWT(MLTM,NT)='TXWT' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
            FNTXBT(MLTM,NT)='TXBT' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
            IF(MTMSRC(MLTM).EQ.2)THEN
              FNTXWF(MLTM,NT)='TXWF' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
              FNTXWC(MLTM,NT)='TXWC' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
              FNTXWP(MLTM,NT)='TXWP' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
              FNTXBF(MLTM,NT)='TXBF' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
              FNTXBC(MLTM,NT)='TXBC' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
              FNTXBP(MLTM,NT)='TXBP' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
            ENDIF
            ENDDO
          ENDIF
        ENDIF
        IF(MTMSRA(MLTM).EQ.1)THEN
          FNAVV(MLTM)='AVVTS' // CNTMSR(MLTM) // '.OUT'
          FNAVB(MLTM)='AVBTS' // CNTMSR(MLTM) // '.OUT'
          FNDIS(MLTM)='DISPR' // CNTMSR(MLTM) // '.OUT'   !DRBC, LZ, 12/14/2020, ADD TKE DISSIPATION RATE OUTPUT       
        ENDIF
        IF(MTMSRP(MLTM).EQ.1)THEN
          FNSEL(MLTM)='SELTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRUE(MLTM).EQ.1)THEN
          FNUVE(MLTM)='UVETS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRUT(MLTM).EQ.1)THEN
          FNUVT(MLTM)='UVTTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRU(MLTM).GE.1)THEN
          FNU3D(MLTM)='U3DTS' // CNTMSR(MLTM) // '.OUT'
          FNV3D(MLTM)='V3DTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRQE(MLTM).EQ.1)THEN
          FNQQE(MLTM)='QQETS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRQ(MLTM).EQ.1)THEN
          FNQ3D(MLTM)='Q3DTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
      ENDDO
      JSTMSR=0
      ENDIF
!
      IF(JSTMSR.EQ.0) GOTO 300
!
      DO MLTM=1,MLTMSR
        IF(MTMSRC(MLTM).GE.1)THEN
          IF(ISTRAN(1).GE.1)THEN
            FNSAL(MLTM)='SALTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(11,FILE=FNSAL(MLTM),STATUS='UNKNOWN')
            CLOSE(11,STATUS='DELETE')
            OPEN(11,FILE=FNSAL(MLTM),STATUS='UNKNOWN')
            WRITE (11,100) TITLE1
            WRITE (11,101) CLTMSR(MLTM)
            WRITE (11,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (11,102) CTUNIT
            CLOSE(11)
          ENDIF
          IF(ISTRAN(2).GE.1)THEN
            FNTEM(MLTM)='TEMTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(21,FILE=FNTEM(MLTM),STATUS='UNKNOWN')
            CLOSE(21,STATUS='DELETE')
            OPEN(21,FILE=FNTEM(MLTM),STATUS='UNKNOWN')
            WRITE (21,100) TITLE2
            WRITE (21,101) CLTMSR(MLTM)
            WRITE (21,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (21,102) CTUNIT
            CLOSE(21)
          ENDIF
          IF(ISTRAN(3).GE.1)THEN
            FNDYE(MLTM)='DYETS' // CNTMSR(MLTM) // '.OUT'
            OPEN(31,FILE=FNDYE(MLTM),STATUS='UNKNOWN')
            CLOSE(31,STATUS='DELETE')
            OPEN(31,FILE=FNDYE(MLTM),STATUS='UNKNOWN')
            WRITE (31,100) TITLE3
            WRITE (31,101) CLTMSR(MLTM)
            WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (31,102) CTUNIT
            CLOSE(31)
          ENDIF
          IF(ISTRAN(4).GE.1)THEN
            FNDYE(MLTM)='SFLTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(31,FILE=FNSFL(MLTM),STATUS='UNKNOWN')
            CLOSE(31,STATUS='DELETE')
            OPEN(31,FILE=FNSFL(MLTM),STATUS='UNKNOWN')
            WRITE (31,100) TITLE3
            WRITE (31,101) CLTMSR(MLTM)
            WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (31,102) CTUNIT
            CLOSE(31)
          ENDIF
          IF(ISTRAN(6).GE.1)THEN
            FNSED(MLTM)='SEDTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(41,FILE=FNSED(MLTM),STATUS='UNKNOWN')
            CLOSE(41,STATUS='DELETE')
            OPEN(41,FILE=FNSED(MLTM),STATUS='UNKNOWN')
            WRITE (41,100) TITLE4
            WRITE (41,101) CLTMSR(MLTM)
            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (41,102) CTUNIT
            CLOSE(41)
          ENDIF
          IF(ISTRAN(7).GE.1)THEN
            DO NX=1,NSND
            FNSND(MLTM,NX)='SND'// CNSND(NX) // 'TS' // CNTMSR(MLTM) // '.OUT'
            OPEN(41,FILE=FNSND(MLTM,NX),STATUS='UNKNOWN')
            CLOSE(41,STATUS='DELETE')
            OPEN(41,FILE=FNSND(MLTM,NX),STATUS='UNKNOWN')
            WRITE (41,100) TITLE4
            WRITE (41,101) CLTMSR(MLTM)
            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (41,102) CTUNIT
            CLOSE(41)
            ENDDO
            DO NX=1,NSND
            FNSBL(MLTM,NX)='SBL'// CNSBL(NX) // 'TS' // CNTMSR(MLTM) // '.OUT'
            OPEN(41,FILE=FNSBL(MLTM,NX),STATUS='UNKNOWN')
            CLOSE(41,STATUS='DELETE')
            OPEN(41,FILE=FNSBL(MLTM,NX),STATUS='UNKNOWN')
            WRITE (41,100) TITLE18
            WRITE (41,101) CLTMSR(MLTM)
            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (41,102) CTUNIT
            CLOSE(41)
            ENDDO
!            FNSND(MLTM)='SNDTS' // CNTMSR(MLTM) // '.OUT'
!            OPEN(41,FILE=FNSND(MLTM),STATUS='UNKNOWN')
!            CLOSE(41,STATUS='DELETE')
!            OPEN(41,FILE=FNSND(MLTM),STATUS='UNKNOWN')
!            WRITE (41,100) TITLE4
!            WRITE (41,101) CLTMSR(MLTM)
!            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
!            WRITE (41,102) CTUNIT
!            CLOSE(41)
          ENDIF
          IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
            FNBED(MLTM)='BEDTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(41,FILE=FNBED(MLTM),STATUS='UNKNOWN')
            CLOSE(41,STATUS='DELETE')
            OPEN(41,FILE=FNBED(MLTM),STATUS='UNKNOWN')
            WRITE (41,100) TITLE19
            WRITE (41,101) CLTMSR(MLTM)
            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (41,102) CTUNIT
            CLOSE(41)
          ENDIF
          IF(ISTRAN(8).GE.1)THEN
            FNDOX(MLTM)='DOXTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(41,FILE=FNDOX(MLTM),STATUS='UNKNOWN')
            CLOSE(41,STATUS='DELETE')
            OPEN(41,FILE=FNDOX(MLTM),STATUS='UNKNOWN')
            WRITE (41,100) TITLE4
            WRITE (41,101) CLTMSR(MLTM)
            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (41,102) CTUNIT
            CLOSE(41)
            FNTOC(MLTM)='TOCTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(42,FILE=FNTOC(MLTM),STATUS='UNKNOWN')
            CLOSE(42,STATUS='DELETE')
            OPEN(42,FILE=FNTOC(MLTM),STATUS='UNKNOWN')
            WRITE (42,100) TITLE4
            WRITE (42,101) CLTMSR(MLTM)
            WRITE (42,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (42,102) CTUNIT
            CLOSE(42)
            FNNHX(MLTM)='NHXTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(43,FILE=FNNHX(MLTM),STATUS='UNKNOWN')
            CLOSE(43,STATUS='DELETE')
            OPEN(43,FILE=FNNHX(MLTM),STATUS='UNKNOWN')
            WRITE (43,100) TITLE4
            WRITE (43,101) CLTMSR(MLTM)
            WRITE (43,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (43,102) CTUNIT
            CLOSE(43)
          ENDIF
          IF(ISTRAN(5).GE.1)THEN
            DO NT=1,NTOX
              FNTOX(MLTM,NT)='TOX' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
              OPEN(51,FILE=FNTOX(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(51,STATUS='DELETE')
              OPEN(51,FILE=FNTOX(MLTM,NT),STATUS='UNKNOWN')
              WRITE (51,100) TITLE5
              WRITE (51,101) CLTMSR(MLTM)
              WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
              WRITE (51,102) CTUNIT
              CLOSE(51)
              FNTXWT(MLTM,NT)='TXWT' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
              OPEN(51,FILE=FNTXWT(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(51,STATUS='DELETE')
              OPEN(51,FILE=FNTXWT(MLTM,NT),STATUS='UNKNOWN')
              WRITE (51,100) TFNTXWT
              WRITE (51,101) CLTMSR(MLTM)
              WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
              WRITE (51,102) CTUNIT
              CLOSE(51)
              FNTXBT(MLTM,NT)='TXBT' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
              OPEN(51,FILE=FNTXBT(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(51,STATUS='DELETE')
              OPEN(51,FILE=FNTXBT(MLTM,NT),STATUS='UNKNOWN')
              WRITE (51,100) TFNTXBT
              WRITE (51,101) CLTMSR(MLTM)
              WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
              WRITE (51,102) CTUNIT
              CLOSE(51)
              IF(MTMSRC(MLTM).EQ.2)THEN
                FNTXWF(MLTM,NT)='TXWF' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
                OPEN(51,FILE=FNTXWF(MLTM,NT),STATUS='UNKNOWN')
                CLOSE(51,STATUS='DELETE')
                OPEN(51,FILE=FNTXWF(MLTM,NT),STATUS='UNKNOWN')
                WRITE (51,100) TFNTXWF
                WRITE (51,101) CLTMSR(MLTM)
                WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
                WRITE (51,102) CTUNIT
                CLOSE(51)
                FNTXWC(MLTM,NT)='TXWC' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
                OPEN(51,FILE=FNTXWC(MLTM,NT),STATUS='UNKNOWN')
                CLOSE(51,STATUS='DELETE')
                OPEN(51,FILE=FNTXWC(MLTM,NT),STATUS='UNKNOWN')
                WRITE (51,100) TFNTXWC
                WRITE (51,101) CLTMSR(MLTM)
                WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
                WRITE (51,102) CTUNIT
                CLOSE(51)
                FNTXWP(MLTM,NT)='TXWP' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
                OPEN(51,FILE=FNTXWP(MLTM,NT),STATUS='UNKNOWN')
                CLOSE(51,STATUS='DELETE')
                OPEN(51,FILE=FNTXWP(MLTM,NT),STATUS='UNKNOWN')
                WRITE (51,100) TFNTXWP
                WRITE (51,101) CLTMSR(MLTM)
                WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
                WRITE (51,102) CTUNIT
                CLOSE(51)
                FNTXBF(MLTM,NT)='TXBF' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
                OPEN(51,FILE=FNTXBF(MLTM,NT),STATUS='UNKNOWN')
                CLOSE(51,STATUS='DELETE')
                OPEN(51,FILE=FNTXBF(MLTM,NT),STATUS='UNKNOWN')
                WRITE (51,100) TFNTXBF
                WRITE (51,101) CLTMSR(MLTM)
                WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
                WRITE (51,102) CTUNIT
                CLOSE(51)
                FNTXBC(MLTM,NT)='TXBC' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
                OPEN(51,FILE=FNTXBC(MLTM,NT),STATUS='UNKNOWN')
                CLOSE(51,STATUS='DELETE')
                OPEN(51,FILE=FNTXBC(MLTM,NT),STATUS='UNKNOWN')
                WRITE (51,100) TFNTXBC
                WRITE (51,101) CLTMSR(MLTM)
                WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
                WRITE (51,102) CTUNIT
                CLOSE(51)
                FNTXBP(MLTM,NT)='TXBP' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
                OPEN(51,FILE=FNTXBP(MLTM,NT),STATUS='UNKNOWN')
                CLOSE(51,STATUS='DELETE')
                OPEN(51,FILE=FNTXBP(MLTM,NT),STATUS='UNKNOWN')
                WRITE (51,100) TFNTXBP
                WRITE (51,101) CLTMSR(MLTM)
                WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
                WRITE (51,102) CTUNIT
                CLOSE(51)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
        IF(MTMSRA(MLTM).EQ.1)THEN
          FNAVV(MLTM)='AVVTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(61,FILE=FNAVV(MLTM),STATUS='UNKNOWN')
          CLOSE(61,STATUS='DELETE')
          OPEN(61,FILE=FNAVV(MLTM),STATUS='UNKNOWN')
          WRITE (61,100) TITLE6
          WRITE (61,101) CLTMSR(MLTM)
          WRITE (61,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (61,102) CTUNIT
          CLOSE(61)
          FNAVB(MLTM)='AVBTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(71,FILE=FNAVB(MLTM),STATUS='UNKNOWN')
          CLOSE(71,STATUS='DELETE')
          OPEN(71,FILE=FNAVB(MLTM),STATUS='UNKNOWN')
          WRITE (71,100) TITLE7
          WRITE (71,101) CLTMSR(MLTM)
          WRITE (71,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (71,102) CTUNIT
          CLOSE(71)
          FNDIS(MLTM)='DISPR' // CNTMSR(MLTM) // '.OUT'   !DRBC, LZ, 12/14/2020, ADD TKE DISSIPATION RATE OUTPUT    
          OPEN(711,FILE=FNDIS(MLTM),STATUS='UNKNOWN')
          CLOSE(711,STATUS='DELETE')
          OPEN(711,FILE=FNDIS(MLTM),STATUS='UNKNOWN')
          WRITE (711,100) TITLE20
          WRITE (711,101) CLTMSR(MLTM)
          WRITE (711,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (711,102) CTUNIT
          CLOSE(711)          
        ENDIF
        IF(MTMSRP(MLTM).EQ.1)THEN
          FNSEL(MLTM)='SELTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(11,FILE=FNSEL(MLTM),STATUS='UNKNOWN')
          CLOSE(11,STATUS='DELETE')
          OPEN(11,FILE=FNSEL(MLTM),STATUS='UNKNOWN')
          WRITE (11,100) TITLE11
          WRITE (11,101) CLTMSR(MLTM)
          WRITE (11,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (11,102) CTUNIT
          CLOSE(11)
        ENDIF
        IF(MTMSRUE(MLTM).EQ.1)THEN
          FNUVE(MLTM)='UVETS' // CNTMSR(MLTM) // '.OUT'
          OPEN(21,FILE=FNUVE(MLTM),STATUS='UNKNOWN')
          CLOSE(21,STATUS='DELETE')
          OPEN(21,FILE=FNUVE(MLTM),STATUS='UNKNOWN')
          WRITE (21,100) TITLE12
          WRITE (21,101) CLTMSR(MLTM)
          WRITE (21,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (21,102) CTUNIT
          CLOSE(21)
        ENDIF
        IF(MTMSRUT(MLTM).EQ.1)THEN
          FNUVT(MLTM)='UVTTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(31,FILE=FNUVT(MLTM),STATUS='UNKNOWN')
          CLOSE(31,STATUS='DELETE')
          OPEN(31,FILE=FNUVT(MLTM),STATUS='UNKNOWN')
          WRITE (31,100) TITLE13
          WRITE (31,101) CLTMSR(MLTM)
          WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (31,102) CTUNIT
          CLOSE(31)
        ENDIF
        IF(MTMSRU(MLTM).GE.1)THEN
          FNU3D(MLTM)='U3DTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(41,FILE=FNU3D(MLTM),STATUS='UNKNOWN')
          CLOSE(41,STATUS='DELETE')
          OPEN(41,FILE=FNU3D(MLTM),STATUS='UNKNOWN')
          WRITE (41,100) TITLE14
          WRITE (41,101) CLTMSR(MLTM)
          WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (41,102) CTUNIT
          CLOSE(41)
          FNV3D(MLTM)='V3DTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(51,FILE=FNV3D(MLTM),STATUS='UNKNOWN')
          CLOSE(51,STATUS='DELETE')
          OPEN(51,FILE=FNV3D(MLTM),STATUS='UNKNOWN')
          WRITE (51,100) TITLE15
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          CLOSE(51)
        ENDIF
        IF(MTMSRQE(MLTM).EQ.1)THEN
          FNQQE(MLTM)='QQETS' // CNTMSR(MLTM) // '.OUT'
          OPEN(61,FILE=FNQQE(MLTM),STATUS='UNKNOWN')
          CLOSE(61,STATUS='DELETE')
          OPEN(61,FILE=FNQQE(MLTM),STATUS='UNKNOWN')
          WRITE (61,100) TITLE16
          WRITE (61,101) CLTMSR(MLTM)
          WRITE (61,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (61,102) CTUNIT
          CLOSE(61)
        ENDIF
        IF(MTMSRQ(MLTM).EQ.1)THEN
          FNQ3D(MLTM)='Q3DTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(71,FILE=FNQ3D(MLTM),STATUS='UNKNOWN')
          CLOSE(71,STATUS='DELETE')
          OPEN(71,FILE=FNQ3D(MLTM),STATUS='UNKNOWN')
          WRITE (71,100) TITLE17
          WRITE (71,101) CLTMSR(MLTM)
          WRITE (71,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (71,102) CTUNIT
          CLOSE(71)
        ENDIF
      ENDDO
!
      JSTMSR=0
!
!----------------------------------------------------------------------C
!
  300 CONTINUE
!
!----------------------------------------------------------------------C
!
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
      ELSE
        TIME=TIMESEC/TCTMSR
      ENDIF
!
      FOURDPI=4./PI
!
! **  STEP CURRENT TIME INTERVALS FOR WRITE SCENARIOS
!
      DO NTSSS=1,NTSSTSP
       DO MTSSS=1,MTSSTSP(NTSSS)
        IF(TIME.GE.TSSTRT(MTSSS,NTSSS))THEN
        IF(TIME.LE.TSSTOP(MTSSS,NTSSS))THEN
          MTSCUR(NTSSS)=MTSSS
        ENDIF
        ENDIF
       ENDDO
      ENDDO
!
      IF(MDCHH.GE.1)THEN
	  DO L=2,LA
	    QCHANUIJ(L)=0.0
	    QCHANVIJ(L)=0.0
	  ENDDO
        DO NMD=1,MDCHH
          LMDCHHT=LMDCHH(NMD)
          LMDCHUT=LMDCHU(NMD)
          LMDCHVT=LMDCHV(NMD)
          IF(MDCHTYP(NMD).EQ.1)THEN
            QCHANUIJ(LMDCHHT)=QCHANUIJ(LMDCHHT)+QCHANU(NMD)
            QCHANUIJ(LMDCHUT)=QCHANUIJ(LMDCHUT)-QCHANU(NMD)
          ENDIF
          IF(MDCHTYP(NMD).EQ.2)THEN
            QCHANVIJ(LMDCHHT)=QCHANVIJ(LMDCHHT)+QCHANV(NMD)
            QCHANVIJ(LMDCHVT)=QCHANVIJ(LMDCHVT)-QCHANV(NMD)
          ENDIF
          IF(MDCHTYP(NMD).EQ.3)THEN
            QCHANUIJ(LMDCHHT)=QCHANUIJ(LMDCHHT)+QCHANU(NMD)
            QCHANUIJ(LMDCHUT)=QCHANUIJ(LMDCHUT)-QCHANU(NMD)
            QCHANVIJ(LMDCHHT)=QCHANVIJ(LMDCHHT)+QCHANV(NMD)
            QCHANVIJ(LMDCHVT)=QCHANVIJ(LMDCHVT)-QCHANV(NMD)
          ENDIF
	  ENDDO
	ENDIF
!
      DO MLTM=1,MLTMSR
       NTSSS=NTSSSS(MLTM)
       MTSCC=MTSCUR(NTSSS)
       IF(TIME.GE.TSSTRT(MTSCC,NTSSS))THEN
       IF(TIME.LE.TSSTOP(MTSCC,NTSSS))THEN
        I=ILTMSR(MLTM)
        J=JLTMSR(MLTM)
        L=LIJ(I,J)
        LN=LNC(L)
        IF(MTMSRC(MLTM).GE.1)THEN
          IF(ISTRAN(1).GE.1)THEN
            OPEN(11,FILE=FNSAL(MLTM),POSITION='APPEND')
            WRITE(11,201)TIME,(SAL(L,K),K=1,KC)
            CLOSE(11)
          ENDIF
          IF(ISTRAN(2).GE.1)THEN
            OPEN(21,FILE=FNTEM(MLTM),POSITION='APPEND')
            WRITE(21,201)TIME,(TEM(L,K),K=1,KC)
            CLOSE(21)
          ENDIF
          IF(ISTRAN(3).GE.1)THEN
            OPEN(31,FILE=FNDYE(MLTM),POSITION='APPEND')
            WRITE(31,201)TIME,(DYE(L,K),K=1,KC)
            CLOSE(31)
          ENDIF
          IF(ISTRAN(4).GE.1)THEN
            OPEN(31,FILE=FNSFL(MLTM),POSITION='APPEND')
            WRITE(31,201)TIME,(SFL(L,K),K=1,KC)
            CLOSE(31)
          ENDIF
          IF(ISTRAN(6).GE.1)THEN
            OPEN(41,FILE=FNSED(MLTM),POSITION='APPEND')
            IF(ISNDAL.EQ.2)THEN
		    SEDBTMP=SEDBT(L,KBT(L))+SEDBT(L,KBT(L)-1)
	      ELSE
		    SEDBTMP=SEDBT(L,KBT(L))
	      ENDIF
            IF(ISEDVW.GE.98)THEN
		    WRITE(41,201)TIME,SEDBTMP,(SEDT(L,K),K=1,KC),FLOCDIA(L,1),WSETFLOC(L,1),SEDFLOCDIA(L,1)
            ELSE
		    WRITE(41,201)TIME,SEDBTMP,(SEDT(L,K),K=1,KC)
	      ENDIF
            CLOSE(41)
          ENDIF
          IF(ISTRAN(7).GE.1)THEN
	      DO NX=1,NSND
            OPEN(41,FILE=FNSND(MLTM,NX),POSITION='APPEND')
            IF(ISNDAL.EQ.2)THEN
		    SNDBTMP=SNDB(L,KBT(L),NX)+SNDB(L,KBT(L)-1,NX)
	      ELSE
		    SNDBTMP=SNDB(L,KBT(L),NX)
	      ENDIF
            WRITE(41,201)TIME,SNDBTMP,(SND(L,K,NX),K=1,KC),SNDEQSAV(L,NX)
            CLOSE(41)
	      ENDDO
	      DO NX=1,NSND
            OPEN(41,FILE=FNSBL(MLTM,NX),POSITION='APPEND')
!            CQBEDLOADX=0.
!            CQBEDLOADY=0.
            IF(UHDYE(L).NE.0.0)CQBEDLOADX(L,NX)=QSBDLDX(L,NX)/UHDYE(L)
            IF(VHDXE(L).NE.0.0)CQBEDLOADY(L,NX)=QSBDLDY(L,NX)/VHDXE(L)
            WRITE(41,201)TIME,QSBDLDX(L,NX),QSBDLDY(L,NX),CQBEDLOADX(L,NX),CQBEDLOADY(L,NX),SNDFBL(L,NX)
            CLOSE(41)
	      ENDDO
!            OPEN(41,FILE=FNSND(MLTM,NX),POSITION='APPEND')
!            WRITE(41,201)TIME,SNDBT(L,KBT(L)),(SNDT(L,K),K=1,KC)
!            CLOSE(41)
          ENDIF
          IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
            OPEN(41,FILE=FNBED(MLTM),POSITION='APPEND')
	      KTMP=KBT(L)
	      KTMP1=KBT(L)-1
	      KTMP1=MAX(KTMP1,1)
	      NSXD=NSED+NSND
!            WRITE(41,221)TIME,KTMP,HBED(L,KTMP),HBED(L,KTMP1),
            WRITE(41,221)TIME,KTMP,HBED(L,KTMP),HBED(L,KTMP1),HBEDA(L),VDRBED(L,KTMP),(VFRBED(L,KTMP,NX),NX=1,NSXD)
            CLOSE(41)
          ENDIF
          IF(ISTRAN(8).GE.1)THEN
            OPEN(41,FILE=FNDOX(MLTM),POSITION='APPEND')
            WRITE(41,201)TIME,(WQV(L,K,19),K=1,KC)
            CLOSE(41)
            OPEN(42,FILE=FNTOC(MLTM),POSITION='APPEND')
            WRITE(42,201)TIME,(WQV(L,K,6),K=1,KC)
            CLOSE(42)
            OPEN(43,FILE=FNNHX(MLTM),POSITION='APPEND')
            WRITE(43,201)TIME,(WQV(L,K,14),K=1,KC)
            CLOSE(43)
          ENDIF
          IF(ISTRAN(5).GE.1)THEN
            DO NT=1,NTOX
!
              OPEN(51,FILE=FNTOX(MLTM,NT),POSITION='APPEND')
              OPEN(52,FILE=FNTXWT(MLTM,NT),POSITION='APPEND')
              OPEN(56,FILE=FNTXBT(MLTM,NT),POSITION='APPEND')
              IF(MTMSRC(MLTM).EQ.2)THEN
			  OPEN(53,FILE=FNTXWF(MLTM,NT),POSITION='APPEND')
                OPEN(54,FILE=FNTXWC(MLTM,NT),POSITION='APPEND')
                OPEN(55,FILE=FNTXWP(MLTM,NT),POSITION='APPEND')
                OPEN(57,FILE=FNTXBF(MLTM,NT),POSITION='APPEND')
                OPEN(58,FILE=FNTXBC(MLTM,NT),POSITION='APPEND')
                OPEN(59,FILE=FNTXBP(MLTM,NT),POSITION='APPEND')
              ENDIF
!
	        NDOC=NSED+NSND+1
              DO K=1,KC
	          SEDSND(K)=SEDT(L,K)+SNDT(L,K)
	          TXWF(K)=TOXFDFW(L,K,NT)*TOX(L,K,NT)
	          TXWC(K)=TOXCDFW(L,K,NT)*TOX(L,K,NT)
	          TXWP(K)=TOXPFTW(L,K,NT)*TOX(L,K,NT)
	        ENDDO
              DO K=1,KB
	          SEDSNDB(K)=SEDBT(L,K)+SNDBT(L,K)
	          TXBF(K)=0.
	          TXBC(K)=0.
	          TXBP(K)=0.
                TXBT(K)=0.
              ENDDO
              DO K=1,KBT(L)
	          PORH(K)=1.0/PORBED(L,K)
	          TXBF(K)=TOXFDFB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
	          TXBC(K)=TOXCDFB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
	          TXBP(K)=TOXPFTB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
	          TXBT(K)=TOXB(L,K,NT)/HBED(L,K)
	        ENDDO
!
	        K=KBT(L)
              WRITE(51,201)TIME,TXBT(K),TXBF(K),TXBC(K),TXBP(K),TOX(L,1,NT),TXWF(1),TXWC(1),TXWP(1)
!                WRITE(51,201)TIME,TOXFDFB(L,K,NT),TOXCDFB(L,K,NT),
!     &          TOXPFTB(L,K,NT),TOXFDFW(L,1,NT),TOXCDFW(L,1,NT),
!     &          TOXPFTW(L,1,NT)
!
              DO K=1,KC
	          IF(SEDSND(K).GT.0.0)TXWP(K)=1000.*TXWP(K)/SEDSND(K)
	        ENDDO
!
              DO K=1,KBT(L)
	          TXBP(K)=TXBP(K)*HBED(L,K)
	          TXBF(K)=TXBF(K)*PORH(K)
	          TXBC(K)=TXBC(K)*PORH(K)
	          IF(SEDSNDB(K).GT.0.0)TXBP(K)=1000.*TXBP(K)/SEDSNDB(K)
	        ENDDO
              WRITE(52,201)TIME,(TOX(L,K,NT),K=1,KC)
              WRITE(56,201)TIME,(TXBT(K),K=1,KB)
              IF(MTMSRC(MLTM).EQ.2)THEN
                WRITE(53,201)TIME,(TXWF(K),K=1,KC)
                WRITE(54,201)TIME,(TXWC(K),K=1,KC)
                WRITE(55,201)TIME,(TXWP(K),K=1,KC)
                WRITE(57,201)TIME,(TXBF(K),K=1,KB)
                WRITE(58,201)TIME,(TXBC(K),K=1,KB)
                WRITE(59,201)TIME,(TXBP(K),K=1,KB)
              ENDIF
!
              CLOSE(51)
              CLOSE(52)
              CLOSE(53)
              CLOSE(54)
              CLOSE(55)
              CLOSE(56)
              CLOSE(57)
              CLOSE(58)
              CLOSE(59)
!
            ENDDO
!            DO NT=1,NTOX
!            OPEN(51,FILE=FNTOXPF(MLTM,NT),POSITION='APPEND')
!            WRITE(51,201)TIME,TOXPFTB(L,KBT(L),NT),
!     &                  (TOXPFTW(L,K,NT),K=1,KC)
!            CLOSE(51)
!            ENDDO
!            DO NT=1,NTOX
!            OPEN(51,FILE=FNTOXFD(MLTM,NT),POSITION='APPEND')
!            WRITE(51,201)TIME,TOXPFTB(L,KBT(L),NT),
!     &                  (TOXPFTW(L,K,NT),K=1,KC)
!            CLOSE(51)
!            ENDDO
          ENDIF
        ENDIF
        IF(MTMSRA(MLTM).EQ.1)THEN
          OPEN(61,FILE=FNAVV(MLTM),POSITION='APPEND')
          OPEN(71,FILE=FNAVB(MLTM),POSITION='APPEND')
          OPEN(711,FILE=FNDIS(MLTM),POSITION='APPEND')          !DRBC, LZ, 12/14/2020, ADD TKE DISSIPATION RATE OUTPUT 
           DO K=1,KS
           ATMP(K)=10000.*AV(L,K)*HP(L)
           ENDDO
          WRITE(61,201)TIME,(ATMP(K),K=1,KS)
           DO K=1,KS
           ATMP(K)=10000.*AB(L,K)*HP(L)
           ENDDO
          WRITE(71,201)TIME,(ATMP(K),K=1,KS)
          WRITE(711,201)TIME,(DISPRATE(L,K),K=1,KC)             !DRBC, LZ, 12/14/2020, ADD TKE DISSIPATION RATE OUTPUT     
          CLOSE(61)
          CLOSE(71)
          CLOSE(711)          
        ENDIF
        IF(MTMSRP(MLTM).EQ.1)THEN
          OPEN(11,FILE=FNSEL(MLTM),POSITION='APPEND')
          PPTMP=HP(L)+BELV(L)
	    TMPVAL=VDWASTE(L)/DXYP(L)
!          HHTMP=PPTMP-BELV(L)
!          GWELTMP=AGWELV(L)-BELAGW(L)
          QVAL=0.0
	    DO K=1,KC
              QVAL=QVAL+QWSEDA(L,K)
          ENDDO
          WRITE(11,211)TIME,PPTMP,HP(L),BELV(L),HBEDA(L),ZELBEDA(L), TMPVAL,VDWASTE(L),QVAL
          CLOSE(11)
        ENDIF
        IF(MTMSRUE(MLTM).EQ.1)THEN
          OPEN(21,FILE=FNUVE(MLTM),POSITION='APPEND')
          UTMP1=50.*(UHDYE(L+1)+UHDYE(L))/(DYP(L)*HP(L))
          VTMP1=50.*(VHDXE(LN)+VHDXE(L))/(DXP(L)*HP(L))
          IF(SPB(L).EQ.0)THEN
            UTMP1=2.*UTMP1
            VTMP1=2.*VTMP1
          ENDIF
          UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
          VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
          UTMP1=5000.*(TBX(L+1)+TBX(L))
          VTMP1=5000.*(TBY(LN)+TBY(L))
          TBEAST=CUE(L)*UTMP1+CVE(L)*VTMP1
          TBNORT=CUN(L)*UTMP1+CVN(L)*VTMP1
!          TAUBDYN=10000.*QQ(L,0)/CTURB2
          TAUBDYN=10000.*TAUB(L)
!          UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
!          VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
          TAUB2=TAUBDYN*TAUBDYN
          IF(ISWAVE.GT.0)THEN
            CURANG=ATAN2(VTMP,UTMP)
            TAUW1=10000.*QQWV1(L)
            TAUW2=10000.*QQWV2(L)
            TAUB2=TAUB2+0.5*(TAUW2*TAUW2)+FOURDPI*TAUBDYN*TAUW2*COS(CURANG-WACCWE(L))
          END IF
          TAUB2=MAX(TAUB2,0.)
          TAUTOT=SQRT(TAUB2)
	    VELMAG=UTMP*UTMP+VTMP*VTMP
	    TMPDRAG=0.0
          TAUBSEDDYN=10000.*TAUBSED(L)
          TAUBSNDDYN=10000.*TAUBSND(L)
          TAUSURFX=10000.*(CUE(L)*TSX(L)+CVE(L)*TSY(L))
          TAUSURFY=10000.*(CUN(L)*TSX(L)+CVN(L)*TSY(L))
	    TAUSURF=SQRT(TAUSURFX*TAUSURFX+TAUSURFY*TAUSURFY)
          BTAUSURFX=10000.*(CUE(L)*FBODYFX(L,KC)+CVE(L)*FBODYFY(L,KC))
          BTAUSURFY=10000.*(CUN(L)*FBODYFX(L,KC)+CVN(L)*FBODYFY(L,KC))
	    IF(VELMAG.GT.0.0)TMPDRAG=TAUTOT/VELMAG
!          WRITE(21,201)TIME,UTMP,VTMP,TBEAST,TBNORT,TAUB,TAUW1,
!     &                 TAUW2,TAUTOT,TMPDRAG
          WRITE(21,201)TIME,UTMP,VTMP,TBEAST,TBNORT,TAUBDYN,
     &    TAUBSEDDYN,TAUBSNDDYN,WINDST(L),WNDVELE(L),WNDVELN(L),
     &    TAUSURF,TAUSURFX,TAUSURFY,BTAUSURFX,BTAUSURFY
!    &        TAUBSEDDYN,TAUBSNDDYN,WINDST(L),TAUSURF,TAUSURFX,TAUSURFY
          CLOSE(21)
        ENDIF
        IF(MTMSRUT(MLTM).EQ.1)THEN
          OPEN(31,FILE=FNUVT(MLTM),POSITION='APPEND')
	    QBEDLOADX=0.
	    QBEDLOADY=0.
          CQBEDLOADXT=0.
          CQBEDLOADYT=0.
          DO NX=1,NSND
            QBEDLOADX=QBEDLOADX+QSBDLDX(L,NX)
            QBEDLOADY=QBEDLOADY+QSBDLDY(L,NX)
            CQBEDLOADXT=CQBEDLOADXT+CQBEDLOADX(L,NX)
            CQBEDLOADYT=CQBEDLOADYT+CQBEDLOADY(L,NX)
	    ENDDO
!          IF(UHDYE(L).NE.0.0)CQBEDLOADX=QBEDLOADX/UHDYE(L)
!          IF(VHDXE(L).NE.0.0)CQBEDLOADY=QBEDLOADY/VHDXE(L)
! jmh modified to write transport by layer in addition to depth average transport
! orginal output with bed load commented out
!          WRITE(31,201)TIME,UHDYE(L),VHDXE(L),QBEDLOADX,QBEDLOADY,
!     &                 CQBEDLOADXT,CQBEDLOADYT
!
      DO K=1,KC
        ATMP(K)=DZC(K)*UHDY(L,K)
        BTMP(K)=DZC(K)*VHDX(L,K)
	ENDDO
      WRITE(31,201)TIME,UHDYE(L),VHDXE(L),(ATMP(K),K=1,KC),(BTMP(K),K=1,KC)
!
! jmh end modified to write transport by layer in addition to depth average transport
          CLOSE(31)
        ENDIF
        IF(MTMSRU(MLTM).GE.1)THEN
          OPEN(41,FILE=FNU3D(MLTM),POSITION='APPEND')
          OPEN(51,FILE=FNV3D(MLTM),POSITION='APPEND')
          RUVTMP=50.
          IF(SPB(L).EQ.0) RUVTMP=100.
          DO K=1,KC
           UTMP1=RUVTMP*(U(L+1,K)+U(L,K))
           VTMP1=RUVTMP*(V(LN,K)+V(L,K))
           ATMP(K)=CUE(L)*UTMP1+CVE(L)*VTMP1
           BTMP(K)=CUN(L)*UTMP1+CVN(L)*VTMP1
          ENDDO
	    IF(MTMSRU(MLTM).EQ.1)THEN
            WRITE(41,201)TIME,(ATMP(K),K=1,KC)
            WRITE(51,201)TIME,(BTMP(K),K=1,KC)
          ENDIF
	    IF(MTMSRU(MLTM).EQ.2)THEN
            WRITE(41,201)TIME,(UHDY(L,K),K=1,KC)
            WRITE(51,201)TIME,(W(L,K),K=0,KC)
          ENDIF
!          WRITE(51,201)TIME,(W(L,K),K=0,KC)
          CLOSE(41)
          CLOSE(51)
        ENDIF
        IF(MTMSRQE(MLTM).EQ.1)THEN
          OPEN(61,FILE=FNQQE(MLTM),POSITION='APPEND')
          QRNT=DXYP(L)*RAINT(L)
	    QEVT=DXYP(L)*EVAPT(L)
          WRITE(61,201)TIME,QSUME(L),QRNT,QEVT,EVAPSW(L),EVAPGW(L),RIFTR(L)
!     &                 ,QCHANUIJ(L),QCHANVIJ(L),QDWASTE(L),VDWASTE(L)
          CLOSE(61)
        ENDIF
        IF(MTMSRQ(MLTM).EQ.1)THEN
          OPEN(71,FILE=FNQ3D(MLTM),POSITION='APPEND')
          WRITE (71,201)TIME,(QSUM(L,K),K=1,KC)
          CLOSE(71)
        ENDIF
       ENDIF
       ENDIF
      ENDDO
!
!**********************************************************************C
!
  100 FORMAT(A80)
  101 FORMAT('  AT LOCATION  ',A20)
  102 FORMAT('  TIME IN FIRST COLUMN HAS UNITS OF ',A10)
  103 FORMAT('  CELL I,J = ',2I5)
  201 FORMAT(F12.5,28E15.7)                                    !hnr
  221 FORMAT(F12.5,I5,28E15.7)                                 !hnr
  211 FORMAT(F12.5,28E15.7)                                    !hnr
!
!**********************************************************************C
!
      RETURN
      END