!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE WQ3DINP
!
!**********************************************************************C
!
!  READ WATER QUALITY SUBMODEL INPUT FILES
!
!  ORGINALLY CODED BY K.-Y. PARK
!  OPTIMIZED AND MODIFIED BY J. M. HAMRICK
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
! hugo rodriguez  5/2009    add bmd water quality file
!
!----------------------------------------------------------------------C
!
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      CHARACTER*3 CWQHDR(NWQVM)
!
      CHARACTER*11  HHMMSS
!
      CHARACTER*12 FNWQSS(22)
!
      DATA IWQTICI,IWQTAGR,IWQTSTL,IWQTSUN,IWQTBEN,IWQTPSL,IWQTNPL/7*0/
      DATA ISMTICI/0/
      IWQTICI=IWQTICI
      IWQTAGR=IWQTAGR
      IWQTSTL=IWQTSTL
      IWQTSUN=IWQTSUN
      IWQTBEN=IWQTBEN
      IWQTPSL=IWQTPSL
      IWQTNPL=IWQTNPL
      ISMTICI=ISMTICI
!
      FNWQSS( 1)='wqCHLaC1.OUT'
      FNWQSS( 2)='wqCHLaG2.OUT'
      FNWQSS( 3)='wqCHLaD3.OUT'
      FNWQSS( 4)='wqRPOCar.OUT'
      FNWQSS( 5)='wqLPOCar.OUT'
      FNWQSS( 6)='wqTDOCar.OUT'
      FNWQSS( 7)='wqRPOPhr.OUT'
      FNWQSS( 8)='wqLPOPhr.OUT'
      FNWQSS( 9)='wqTDOPhr.OUT'
      FNWQSS(10)='wqP04Tot.OUT'
      FNWQSS(11)='wqRPONit.OUT'
      FNWQSS(12)='wqLPONit.OUT'
      FNWQSS(13)='wqTDONit.OUT'
      FNWQSS(14)='wqNH4Nit.OUT'
      FNWQSS(15)='wqNO3Nit.OUT'
      FNWQSS(16)='wqSilicA.OUT'
      FNWQSS(17)='wqSilicU.OUT'
      FNWQSS(18)='wqChemOD.OUT'
      FNWQSS(19)='wqDisOxy.OUT'
      FNWQSS(20)='wqTAMeta.OUT'
      FNWQSS(21)='wqFColBa.OUT'
      FNWQSS(22)='wqBenAlg.OUT'
!
      DO NF=1,22
	  OPEN(99,FILE=FNWQSS(NF),STATUS='UNKNOWN')
	  CLOSE(99,STATUS='DELETE')
	ENDDO
!
!TT      OPEN(1,FILE='WQWCTS.OUT',STATUS='UNKNOWN')
!TT      CLOSE(1,STATUS='DELETE')
!
      OPEN(1,FILE='WQ3D.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
! **  HARDWIRE BY PASS OF RATE COEFFICIENT MAPS
!
! MRM      ISWQCMAP=0
! MRM      ISWQSMAP=0
!
      NWQKCNT=0
      NWQKDPT=1
!
      UHEQ(1)=0.0
      UHEQ(LC)=0.0
      DO ND=1,NDMWQ
       LF=2+(ND-1)*LDMWQ
       LL=LF+LDM-1
       DO L=LF,LL
        UHEQ(L)=1.0
       ENDDO
      ENDDO
!H-
      TINDAY = 0.0
      TINDAY=TINDAY
      ITNWQ = 0
!
      RKCWQ = 1.0/REAL(KC)
      DO K=1,KC
        WQHT(K)=REAL(KC-K)*RKCWQ
      ENDDO
!
!XH      INWQRST=21
!XH      INWQICI=21
!XH      INWQAGR=22
!XH      INWQSTL=23
!XH      INWQSUN=24
!XH      INWQBEN=25
!XH      INWQPSL=26
!XH      INWQNPL=27
!XH      IWQONC=30
!XH      IWQORST=31
!
!      WQTSNAME(1)  = 'CHL'
!      WQTSNAME(2)  = 'TPC'
!      WQTSNAME(3)  = 'DOC'
!      WQTSNAME(4)  = 'TPP'
!      WQTSNAME(5)  = 'DOP'
!      WQTSNAME(6)  = 'P4T'
!      WQTSNAME(7)  = 'P4D'
!      WQTSNAME(8)  = 'APC'
!      WQTSNAME(9)  = 'TNN'
!      WQTSNAME(10) = 'DON'
!      WQTSNAME(11) = 'NH4'
!      WQTSNAME(12) = 'NO3'
!      WQTSNAME(13) = 'TSI'
!      WQTSNAME(14) = 'SUU'
!      WQTSNAME(15) = 'SAA'
!      WQTSNAME(16) = 'SAD'
!      WQTSNAME(17) = 'COD'
!      WQTSNAME(18) = 'DOO'
!      WQTSNAME(19) = 'TAM'
!      WQTSNAME(20) = 'TMP'
!      WQTSNAME(21) = 'FCB'
!
      WQTSNAME(1)  = 'CHC'
      WQTSNAME(2)  = 'CHG'
      WQTSNAME(3)  = 'CHD'
      WQTSNAME(4)  = 'ROC'
      WQTSNAME(5)  = 'LOC'
      WQTSNAME(6)  = 'DOC'
!     WQTSNAME(7)  = 'P4D'
!     WQTSNAME(8)  = 'ROP'
!     WQTSNAME(9)  = 'LOP'
!     WQTSNAME(10) = 'DOP'
      WQTSNAME(7)  = 'ROP'
      WQTSNAME(8)  = 'LOP'
      WQTSNAME(9)  = 'DOP'
      WQTSNAME(10) = 'PO4'
      WQTSNAME(11) = 'RON'
      WQTSNAME(12) = 'LON'
      WQTSNAME(13) = 'DON'
      WQTSNAME(14) = 'NHX'
      WQTSNAME(15) = 'NOX'
      WQTSNAME(16) = 'SUU'
      WQTSNAME(17) = 'SAA'
      WQTSNAME(18) = 'COD'
      WQTSNAME(19) = 'DOX'
      WQTSNAME(20) = 'TAM'
      WQTSNAME(21) = 'FCB'
      WQTSNAME(22) = 'MAC'

      DO M=0,NWQPS
       WQPSQ(M)=0.0
       WQPSQC(M)=0.0
       DO J=1,NWQV
        WQWPSLC(M,J)=0.0
       ENDDO
      ENDDO
!
      DO K=1,KC
       IWQPSC(1,K)=0
       WQDSQ(1,K)=0.0
       IWQPSC(LC,K)=0
       WQDSQ(LC,K)=0.0
      ENDDO
!
      DO ND=1,NDMWQ
       LF=2+(ND-1)*LDMWQ
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         IWQPSC(L,K)=0
         IWQPSV(L,K)=0
         WQDSQ(L,K)=0.0
        ENDDO
       ENDDO
      ENDDO
!
      DO J=1,NWQV
       DO K=1,KC
        DO L=1,LC
         WQWDSL(L,K,J)=0.0
         WQWPSL(L,K,J)=0.0
        ENDDO
       ENDDO
      ENDDO
!
      CALL RWQC1(IWQDT)
!      CALL RWQC2
!      CALL RWQMAP
!
!
      OPEN(1,FILE='WQWCTS.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WQWCTS.OUT',STATUS='UNKNOWN')
!
      OPEN(99,FILE='WQWCSNAPSHOT.OUT',STATUS='UNKNOWN')
      CLOSE(99,STATUS='DELETE')
      OPEN(99,FILE='WQWCSNAPSHOT.OUT',STATUS='UNKNOWN')
!
      NWQVOUT=0
      DO NW=1,NWQV
       IF(ISTRWQ(NW).EQ.1)THEN
         NWQVOUT=NWQVOUT+1
         CWQHDR(NWQVOUT)=WQTSNAME(NW)
       ENDIF
      ENDDO
!
! Sen Bai
! Add macroalgae output
	NWQVOUT=NWQVOUT+1
      CWQHDR(NWQVOUT)=WQTSNAME(NWQV+1)
! End of Sen Bai
      WRITE(1,1969)(CWQHDR(NW),NW=1,NWQVOUT)
      WRITE(99,1970)NWQVOUT,(CWQHDR(NW),NW=1,NWQVOUT)
!
 1969 FORMAT('C   I    J    K    TIME',7X,A3,8X,A3,8X,A3,                                         &
       8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,                                                       &
       8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,                                                       &
       8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3)
 1970 FORMAT(I5,7X,A3,8X,A3,8X,A3,                                                                &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,                                                      &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,                                                      &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3)
!
! 1969 FORMAT('C   I    J    K    TIME       CHL        TPC',
!     &       '        DOC        TPP        DOP        P4T',
!     &       '        P4D        APC        TNN        DON',
!     &       '        NH4        NO3        TSI        SUU',
!     &       '        SAA        SAD        COD        DOO',
!     &       '        TAM        TMP        FCB        MALG')
!
      CLOSE(1)
      CLOSE(99)
!
! **  INITIALIZE DIURNAL DO ANALYSIS
!
      IF(NDDOAVG.GE.1)THEN
        OPEN(1,FILE='DIURNDO.OUT')
        CLOSE(1,STATUS='DELETE')
        DO K=1,KC
         DO L=2,LA
          DDOMAX(L,K)=-1.E6
          DDOMIN(L,K)=1.E6
         ENDDO
        ENDDO
      ENDIF
!
! **  INITIALIZE LIGHT EXTINCTION ANALYSIS
!
      IF(NDLTAVG.GE.1)THEN
        OPEN(1,FILE='LIGHT.OUT')
        CLOSE(1,STATUS='DELETE')
        NDLTCNT=0
        DO K=1,KC
         DO L=2,LA
          RLIGHTT(L,K)=0.
          RLIGHTC(L,K)=0.
         ENDDO
        ENDDO
      ENDIF
!
!
      IF(IWQICI.EQ.2) CALL RWQRST
!
      IF(IWQBEN.EQ.1)THEN
        CALL SMINIT
        CALL SMRIN1(IWQDT)
!        CALL SMRIN2
!        CALL RSMMAP
        IF(ISMICI.EQ.2) CALL RSMRST
      ENDIF
!
!      OPEN(1,FILE='TFILE',STATUS='UNKNOWN')
!      CALL TIME(HHMMSS)
!      CLOSE(1,STATUS='DELETE')
!      OPEN(1,FILE='TFILE',STATUS='UNKNOWN')
!      WRITE(1,100) HHMMSS
!      CLOSE(1)
!
  100 FORMAT('  TIME = ',A11,' HH.MM.SS.HH')
!MRM   100 FORMAT('  TIME = ',A11,' HH:MM:SS.SS')
!
      RETURN
      END