!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
!TT   SUBROUTINE RWQSUN
!
!**********************************************************************C
!
! READ IN TEMPORALLY VARYING PARAMETERS FOR DAILY SOLAR RADIATION (WQI0)
! AND FRACTIONAL DAYLENGTH (WQFD) (UNIT INWQSUN).
!
!**********************************************************************C
!
!TT     INCLUDE 'EFDC.PAR'
!TT      INCLUDE 'EFDC.CMN'
!
!TT      CHARACTER TITLE(3)*79, SUNCONT*3
!
!TT      OPEN(1,FILE=SUNFN,STATUS='UNKNOWN')
!TT      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
!
!TT      IF(IWQTSUN.EQ.0)THEN
!TT        READ(1,50) (TITLE(M),M=1,3)
!TT        WRITE(2,999)
!TT        WRITE(2,50) (TITLE(M),M=1,3)
!TT      ENDIF
!
!TT      WRITE(2,60)'* IO & FD AT            ', IWQTSUN,
!TT     *  ' TH DAY FROM MODEL START'
!
!TT      READ(1,999)
!TT      READ(1,50) TITLE(1)
!TT      WRITE(2,50) TITLE(1)
!TT      READ(1,53) WQI0,WQFD
!TT      WRITE(2,53) WQI0,WQFD
!
!TT      IF(IWQTSUN.EQ.0)THEN
!TT        WQI1 = WQI0
!TT        WQI2 = WQI0
!TT      ENDIF
!
!TT      READ(1,52) IWQTSUN, SUNCONT
!TT      WRITE(2,52) IWQTSUN, SUNCONT
!TT      IF(SUNCONT.EQ.'END')THEN
!TT        CLOSE(1)
!TT        IWQSUN = 0
!TT      ENDIF
!
!TT      CLOSE(2)
!
!TT  999 FORMAT(1X)
!TT   50 FORMAT(A79)
!TT   52 FORMAT(I7, 1X, A3)
!TT   53 FORMAT(10F8.3)
!TT   60 FORMAT(/, A24, I5, A24)
!
!TT      RETURN
!TT      END
!
!**********************************************************************C
!
      SUBROUTINE RWQSUN
!
!**********************************************************************C
!
! **  NEW VERSION BY J. M. HAMRICK  7 APRIL 1997
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
! **  READS AND INTERPOLATES DAILY AVERAGE SOLAR RADIATION AND
! **  DAYLIGHT FRACTION
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      IF(ITNWQ.GT.0) GOTO 1000
!
!**********************************************************************C
!
! **  READ IN DAILY AVERAGE SOLAR SW RAD SERIES FROM FILE 'SUNDAY.INP'
!
!----------------------------------------------------------------------C
!
      OPEN(1,FILE='SUNDAY.INP',STATUS='UNKNOWN')
!
! **  SKIP OVER TITLE AND AND HEADER LINES
!
      DO IS=1,7
      READ(1,1)
      ENDDO
!
      M=0
      ISPAR=1
!      MCSUNDAY=1 TEMP USE ISPAR FOR MCSUNDAY
      READ(1,*,IOSTAT=ISO)NSUNDAY,TCSUNDAY,TASUNDAY,RMULADJ,ADDADJ
      IF(ISO.GT.0) GOTO 900
      DO M=1,NSUNDAY
        READ(1,*,IOSTAT=ISO)TSSRD(M),SOLSRD(M),SOLFRD(M)
        IF(ISO.GT.0) GOTO 900
        TSSRD(M)=TCSUNDAY*( TSSRD(M)+TASUNDAY )
        SOLSRD(M)=RMULADJ*(SOLSRD(M)+ADDADJ) * PARADJ
      ENDDO
!
      CLOSE(1)
!
      GOTO 901
!
  900 CONTINUE
      WRITE(6,601)M
      STOP
!
  901 CONTINUE
!
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR FILE SUNDAY.INP ')
!
!**********************************************************************C
!
 1000 CONTINUE
!
!**********************************************************************C
!
! **  DAILY AVERAGE SOLAR SW RADIATION INTERPOLTATION FOR WATER QUALITY
!
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      ELSE
        TIME=TIMESEC/86400.
      ENDIF
!
       M1=ISPAR
!      TEMP USE ISPAR FOR MCSUNDAY
!
  100  CONTINUE
       M2=M1+1
       IF(TIME.GT.TSSRD(M2))THEN
         M1=M2
         GOTO 100
        ELSE
         ISPAR=M1
!      TEMP USE ISPAR FOR MCSUNDAY
       ENDIF
!
       TDIFF=TSSRD(M2)-TSSRD(M1)
       WTM1=(TSSRD(M2)-TIME)/TDIFF
       WTM2=(TIME-TSSRD(M1))/TDIFF
       SOLSRDT=WTM1*SOLSRD(M1)+WTM2*SOLSRD(M2)
       SOLFRDT=WTM1*SOLFRD(M1)+WTM2*SOLFRD(M2)
!
!**********************************************************************C
!
      RETURN
      END