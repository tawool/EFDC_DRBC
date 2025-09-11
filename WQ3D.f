C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
C **  WATER QUALITY MODEL SUBROUTINE START HERE
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQ3D
C
C**********************************************************************C
C
C  CONTROL SUBROUTINE FOR WATER QUALITY MODEL
C
C  ORGINALLY CODED BY K.-Y. PARK
C  OPTIMIZED AND MODIFIED BY J. M. HAMRICK
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C  hugo rodriguez     added bmd file for wq
c  hugo rodriguez   10/2010 change format of algaegro and settling time varying file to morton's bentic format
C
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C     CHARACTER*11  HHMMSS
C
      DATA IWQTICI,IWQTAGR,IWQTSTL,IWQTSUN,IWQTBEN,IWQTPSL,IWQTNPL/7*0/
      DATA ISMTICI/0/
C
C**********************************************************************C
C
      IWQTSUN=IWQTSUN
      IWQTBEN=IWQTBEN
      IWQTPSL=IWQTPSL
      IWQTNPL=IWQTNPL
C
C **  READ INITIAL CONDITIONS
C
      IF(IWQICI.EQ.1 .AND. ITNWQ.EQ.IWQTICI) CALL RWQICI(IWQTICI)
C
C **  UPDATE TIME SERIES BOUNDARY CONDITIONS
C
      CALL RWQCSR
C
C **  READ TIME/SPACE VARYING ALGAE PARAMETERS
C
C  IF SIMULATION TIME IS >= THE NEXT TIME IN THE AGR FILE.
      IF(IWQAGR .EQ. 1)THEN                                                !HNR
        IF(ISDYNSTP.EQ.0)THEN
          TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        ELSE
          TIMTMP=TIMESEC/86400.
        ENDIF
      ENDIF                                                                !HNR
C
C **  READ TIME/SPACE VARYING SETTLING VELOCITIES
C
c      IF(IWQSTL.EQ.1 .AND. ITNWQ.EQ.IWQTSTL) CALL RWQSTL(IWQTSTL)
C  CALL SPATIALLY AND TIME VARYING BENTHIC FLUX HERE.  ONLY CALL RWQstl1
C  IF SIMULATION TIME IS >= THE NEXT TIME IN THE AGR FILE.
      IF(IWQSTL .EQ. 1)THEN                                                !HNR
        IF(ISDYNSTP.EQ.0)THEN
          TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        ELSE
          TIMTMP=TIMESEC/86400.
        ENDIF
        IF(TIMTMP .GE. stlDAY)THEN                                         !HNR
          CALL RWQSTL1(TIMTMP)                                             !HNR
        ENDIF                                                              !HNR
      ENDIF                                                                !HNR

C ADDED BY M. MORTON 03/29/98
C  IF SIMULATION TIME IS >= THE NEXT TIME IN THE BENTHIC FILE.
C
C **  UPDATE SOLAR RADIATION INTENSITY
C UPDATE OCCURS ONLY WHEN THE SIMULATION DAY CHANGES.
C   WQI1 = SOLAR RADIATION ON PREVIOUS DAY
C   WQI2 = SOLAR RADIATION TWO DAYS AGO
C   WQI3 = SOLAR RADIATION THREE DAYS AGO
C
      ISTPDAY = INT((86400.0/TIDALP)*REAL(NTSPTC))
      IF(MOD(N, ISTPDAY) .EQ. 0)THEN
        WQI3 = WQI2
        WQI2 = WQI1
        WQI1 = WQI0OPT
        WQI0OPT = 0.0
      ENDIF
C
C **  READ SOLAR RADIATION INTENSITY AND DAYLIGHT LENGTH
C
C NOTE:  IWQSUN=1 CALLS SUBROUTINE RWQSUN WHICH READS THE DAILY
C                 SOLAR RADIATION DATA FROM FILE SUNDAY.INP WHICH
C                 ARE IN UNITS OF LANGLEYS/DAY.
C        IWQSUN=2 USES THE HOURLY SOLAR RADIATION DATA FROM ASER.INP
C                 AND CONVERT WATTS/M**2 TO LANGLEYS/DAY USING 2.065.
C                 AND ADJUST FOR PHOTOSYNTHETIC ACTIVE RADIATION BY 0.43
C
      IF(IWQSUN.EQ.1)THEN
        CALL RWQSUN
        WQI0=SOLSRDT
        WQFD=SOLFRDT
      ENDIF
C
      IF(IWQSUN.EQ.2)THEN
        SOLARAVG=0.
        DO L=2,LA
          SOLARAVG=SOLARAVG+SOLSWRT(L)
        ENDDO
        SOLARAVG=SOLARAVG/FLOAT(LA-1)
        WQI0 = PARADJ*2.065*SOLARAVG
        WQI0OPT = MAX(WQI0OPT, WQI0)
        WQFD=1.
      ENDIF
C
C **  READ BENTHIC FLUX IF REQUIRED
C
C      IF(IWQBEN.EQ.2 .AND. ITNWQ.EQ.IWQTBEN) CALL RWQBEN(IWQTBEN)          !MRM
C ADDED BY M. MORTON 03/29/98
C  IF SIMULATION TIME IS >= THE NEXT TIME IN THE BENTHIC FILE.
      IF(IWQBEN .EQ. 2)THEN                                                !MRM
      IF(ISDYNSTP.EQ.0)THEN
        TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      ELSE
        TIMTMP=TIMESEC/86400.
      ENDIF
      ENDIF                                                                 !MRM
C
C **  UPDATE POINT SOURCE LOADINGS
C
      IF(IWQPSL.GE.1) CALL RWQPSL
C
C
C **  READ NON POINT SOURCE LOADING
C
*     IF(IWQNPL.EQ.1 .AND. ITNWQ.EQ.IWQTNPL) CALL RWQNPL(IWQTNPL)
C
C **  READ SEDIMENT MODEL INITIAL CONDITION
C
      IF(IWQBEN.EQ.1)THEN
        IF(ISMICI.EQ.1 .AND. ITNWQ.EQ.ISMTICI) CALL RSMICI(ISMTICI)
      ENDIF
C
C **  UPDATE TIME STEP IN DAYS FOR DYNAMIC TIME STEPPING
C
      IF(IS2TIM.GE.1) THEN
        IF(ISDYNSTP.NE.0)THEN
          DTWQ=DTDYN/86400.
      DTWQO2=0.5*DTWQ                   !added by John for dyn time step
        END IF
	ENDIF
C
C **  UPDATE TIME IN DAYS
C
      ITNWQ = ITNWQ + 2
      TINDAY = TINDAY + DTWQ
C
C **  UPDATE OLD CONCENTRATIONS
C
CTT       DO NW=1,NWQV
CTT         DO K=1,KC
CTT           DO L=2,LA
CTT              WQVO(L,K,NW)=WQV(L,K,NW)
CTT           ENDDO
CTT         ENDDO
CTT       ENDDO
C
C MRM KEEP TRACK OF D.O. TRANSPORT COMPONENT (I.E., THE NEW D.O.
C   FOLLOWING THE CALL TO CALWQC MINUS OLD D.O. BEFORE THE CALL).
C   FIRST SUBTRACT THE OLD D.O. HERE:
C
      DO K=1,KC
        DO L=2,LA
          XMRM = WQV(L,K,19)*DTWQ*DZC(K)*HP(L)
          XDOTRN(L,K) = XDOTRN(L,K) - XMRM
          XDOALL(L,K) = XDOALL(L,K) - XMRM
        ENDDO
      ENDDO
C
C
C MRM KEEP TRACK OF D.O. TRANSPORT COMPONENT (I.E., THE NEW D.O.
C   FOLLOWING THE CALL TO CALWQC MINUS OLD D.O. BEFORE THE CALL).
C   NOW ADD THE NEW D.O. HERE:
C
      DO K=1,KC
        DO L=2,LA
          XMRM = WQV(L,K,19)*DTWQ*DZC(K)*HP(L)
          XDOTRN(L,K) = XDOTRN(L,K) + XMRM
          XDOALL(L,K) = XDOALL(L,K) + XMRM
        ENDDO
      ENDDO
C
C **  LOAD WQV INTO WQVO FOR REACTION CALCULATION
C
C J.S. ADD MACLGAE
C
       NMALG=0
       IF(IDNOTRVA.GT.0) NMALG=1
       DO NW=1,NWQV+NMALG
         DO K=1,KC
           DO L=2,LA
              WQVO(L,K,NW)=WQV(L,K,NW)
           ENDDO
         ENDDO
       ENDDO
C
C **  UPDATE WATER COLUMN KINETICS AND SEDIMENT MODEL
C **  OVER LONGER TIME INTERVALS THAN PHYSICAL TRANSPORT
C **  IF NWQKDPT .GT. 1
C
      NWQKCNT=NWQKCNT+1
      IF(NWQKCNT.EQ.NWQKDPT)THEN
        NWQKCNT=0
      ENDIF
C
C **    CALCULATE KINETIC SOURCES AND SINKS
C
        IF(ISCRAY.EQ.0)THEN
          TTMP=SECNDS(0.0)
         ELSE
          T1TMP=SECOND( )
          CALL TIMEF(WT1TMP)
        ENDIF
C
C
        IF(ISCRAY.EQ.0)THEN
          TWQKIN=TWQKIN+SECNDS(TTMP)
         ELSE
          T2TMP=SECOND( )
          CALL TIMEF(WT2TMP)
          TWQKIN=TWQKIN+T2TMP-T1TMP
          WTWQKIN=WTWQKIN+(WT2TMP-WT1TMP)*0.001
        ENDIF
C
C **    DIAGNOSE NEGATIVE CONCENTRATIONS
C
        CALL WWQNC
C
C **    WRITE TIME SERIES
C
        IF(ITNWQ.GE.IWQTSB .AND. ITNWQ.LE.IWQTSE)THEN
C           IF(MOD(ITNWQ,IWQTSDT).EQ.0) CALL WWQTS(TINDAY)
c           IF(MOD(ITNWQ,IWQTSDT).EQ.0) CALL WWQTS                        
           IF(MOD(ITNWQ,IWQTSDT).EQ.0) then                               
c             CALL WWQTS                                                  
C          CALL WWQTSBIN
        ENDIF
C
C **    WRITE SNAP SHOT  IWQSSB,IWQSSE,IWQSSDT
C
        IF(ITNWQ.GE.IWQSSB .AND. ITNWQ.LE.IWQSSE)THEN
           IF(MOD(ITNWQ,IWQSSDT).EQ.0) CALL WWQSNAPSHOT
        ENDIF
C
C **    CALL SEDIMENT MODEL
C
        IF(IWQBEN.EQ.1)THEN
C
          IF(ISCRAY.EQ.0)THEN
            TTMP=SECNDS(0.0)
           ELSE
            T1TMP=SECOND( )
            CALL TIMEF(WT1TMP)
          ENDIF
C
          CALL SMMBE
C
          IF(ISCRAY.EQ.0)THEN
            TWQSED=TWQSED+SECNDS(TTMP)
           ELSE
            T2TMP=SECOND( )
            CALL TIMEF(WT2TMP)
            TWQSED=TWQSED+T2TMP-T1TMP
            WTWQSED=WTWQSED+(WT2TMP-WT1TMP)*0.001
          ENDIF
C
          IF(ISMTS.GE.1)THEN
C
C **      WRITE SEDIMENT MODEL TIME SERIES
C
            IF(ITNWQ.GE.ISMTSB .AND. ITNWQ.LE.ISMTSE)THEN
C              IF(MOD(ITNWQ,ISMTSDT).EQ.0) CALL WSMTS(TINDAY)
              IF(MOD(ITNWQ,ISMTSDT).EQ.0) CALL WSMTS
            ENDIF
C
          ENDIF
C
C **      WRITE SEDIMENT MODEL FLUXES TO BINARY FILE:
C
          IF(ITNWQ.GE.ISMTSB .AND. ITNWQ.LE.ISMTSE)THEN
            CALL WSMTSBIN
          ENDIF
        ENDIF
C
      ENDIF
C
C **  ENDIF ON KINETIC AND SEDIMENT UPDATE
C
C **  INSERT TIME CALL
C
C     CALL TIME(HHMMSS)
C     OPEN(1,FILE='TFILE',STATUS='UNKNOWN',POSITION='APPEND')
C     WRITE(1,100) HHMMSS
C     CLOSE(1)
C
C **  WRITE RESTART FILES
C
C      IF(IWQRST.EQ.1) CALL WWQRST
C      IF(IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST
C
CXH      CLOSE(IWQONC)
C
  100 FORMAT('  TIME = ',A11,' HH.MM.SS.HH')
  600 FORMAT(' ITNWQ,IWQTSDT,ITMP,IWQTSB,TWQTSE = ',5I5)
C
      RETURN
      END
