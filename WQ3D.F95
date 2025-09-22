!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
! **  WATER QUALITY MODEL SUBROUTINE START HERE
!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE WQ3D
!
!**********************************************************************C
!
!  CONTROL SUBROUTINE FOR WATER QUALITY MODEL
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
!  hugo rodriguez     added bmd file for wq
!  hugo rodriguez   10/2010 change format of algaegro and settling time varying file to morton's bentic format
!
!----------------------------------------------------------------------C
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!     CHARACTER*11  HHMMSS
!
      DATA IWQTICI,IWQTAGR,IWQTSTL,IWQTSUN,IWQTBEN,IWQTPSL,IWQTNPL/7*0/
      DATA ISMTICI/0/
!
!**********************************************************************C
!
      IWQTSUN=IWQTSUN
      IWQTBEN=IWQTBEN
      IWQTPSL=IWQTPSL
      IWQTNPL=IWQTNPL
!
! **  READ INITIAL CONDITIONS
!
      IF(IWQICI.EQ.1 .AND. ITNWQ.EQ.IWQTICI) CALL RWQICI(IWQTICI)
!
! **  UPDATE TIME SERIES BOUNDARY CONDITIONS
!
      CALL RWQCSR
!
! **  READ TIME/SPACE VARYING ALGAE PARAMETERS
!
!  IF SIMULATION TIME IS >= THE NEXT TIME IN THE AGR FILE.
      IF(IWQAGR .EQ. 1)THEN                                                !HNR
        IF(ISDYNSTP.EQ.0)THEN
          TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        ELSE
          TIMTMP=TIMESEC/86400.
        ENDIF
      ENDIF                                                                !HNR
!
! **  READ TIME/SPACE VARYING SETTLING VELOCITIES
!
!      IF(IWQSTL.EQ.1 .AND. ITNWQ.EQ.IWQTSTL) CALL RWQSTL(IWQTSTL)
!  CALL SPATIALLY AND TIME VARYING BENTHIC FLUX HERE.  ONLY CALL RWQstl1
!  IF SIMULATION TIME IS >= THE NEXT TIME IN THE AGR FILE.
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

! ADDED BY M. MORTON 03/29/98
!  IF SIMULATION TIME IS >= THE NEXT TIME IN THE BENTHIC FILE.
!
! **  UPDATE SOLAR RADIATION INTENSITY
! UPDATE OCCURS ONLY WHEN THE SIMULATION DAY CHANGES.
!   WQI1 = SOLAR RADIATION ON PREVIOUS DAY
!   WQI2 = SOLAR RADIATION TWO DAYS AGO
!   WQI3 = SOLAR RADIATION THREE DAYS AGO
!
      ISTPDAY = INT((86400.0/TIDALP)*REAL(NTSPTC))
      IF(MOD(N, ISTPDAY) .EQ. 0)THEN
        WQI3 = WQI2
        WQI2 = WQI1
        WQI1 = WQI0OPT
        WQI0OPT = 0.0
      ENDIF
!
! **  READ SOLAR RADIATION INTENSITY AND DAYLIGHT LENGTH
!
! NOTE:  IWQSUN=1 CALLS SUBROUTINE RWQSUN WHICH READS THE DAILY
!                 SOLAR RADIATION DATA FROM FILE SUNDAY.INP WHICH
!                 ARE IN UNITS OF LANGLEYS/DAY.
!        IWQSUN=2 USES THE HOURLY SOLAR RADIATION DATA FROM ASER.INP
!                 AND CONVERT WATTS/M**2 TO LANGLEYS/DAY USING 2.065.
!                 AND ADJUST FOR PHOTOSYNTHETIC ACTIVE RADIATION BY 0.43
!
      IF(IWQSUN.EQ.1)THEN
        CALL RWQSUN
        WQI0=SOLSRDT
        WQFD=SOLFRDT
      ENDIF
!
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
!
! **  READ BENTHIC FLUX IF REQUIRED
!
!      IF(IWQBEN.EQ.2 .AND. ITNWQ.EQ.IWQTBEN) CALL RWQBEN(IWQTBEN)          !MRM
! ADDED BY M. MORTON 03/29/98
!  IF SIMULATION TIME IS >= THE NEXT TIME IN THE BENTHIC FILE.
      IF(IWQBEN .EQ. 2)THEN                                                !MRM
      IF(ISDYNSTP.EQ.0)THEN
        TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      ELSE
        TIMTMP=TIMESEC/86400.
      ENDIF
      ENDIF                                                                 !MRM
!
! **  UPDATE POINT SOURCE LOADINGS
!
      IF(IWQPSL.GE.1) CALL RWQPSL
!
!
! **  READ NON POINT SOURCE LOADING
!
!     IF(IWQNPL.EQ.1 .AND. ITNWQ.EQ.IWQTNPL) CALL RWQNPL(IWQTNPL)
!
! **  READ SEDIMENT MODEL INITIAL CONDITION
!
      IF(IWQBEN.EQ.1)THEN
        IF(ISMICI.EQ.1 .AND. ITNWQ.EQ.ISMTICI) CALL RSMICI(ISMTICI)
      ENDIF
!
! **  UPDATE TIME STEP IN DAYS FOR DYNAMIC TIME STEPPING
!
      IF(IS2TIM.GE.1) THEN
        IF(ISDYNSTP.NE.0)THEN
          DTWQ=DTDYN/86400.
      DTWQO2=0.5*DTWQ                   !added by John for dyn time step
        END IF
	ENDIF
!
! **  UPDATE TIME IN DAYS
!
      ITNWQ = ITNWQ + 2
      TINDAY = TINDAY + DTWQ
!
! **  UPDATE OLD CONCENTRATIONS
!
!TT       DO NW=1,NWQV
!TT         DO K=1,KC
!TT           DO L=2,LA
!TT              WQVO(L,K,NW)=WQV(L,K,NW)
!TT           ENDDO
!TT         ENDDO
!TT       ENDDO
!
! MRM KEEP TRACK OF D.O. TRANSPORT COMPONENT (I.E., THE NEW D.O.
!   FOLLOWING THE CALL TO CALWQC MINUS OLD D.O. BEFORE THE CALL).
!   FIRST SUBTRACT THE OLD D.O. HERE:
!
      DO K=1,KC
        DO L=2,LA
          XMRM = WQV(L,K,19)*DTWQ*DZC(K)*HP(L)
          XDOTRN(L,K) = XDOTRN(L,K) - XMRM
          XDOALL(L,K) = XDOALL(L,K) - XMRM
        ENDDO
      ENDDO
!
!
! MRM KEEP TRACK OF D.O. TRANSPORT COMPONENT (I.E., THE NEW D.O.
!   FOLLOWING THE CALL TO CALWQC MINUS OLD D.O. BEFORE THE CALL).
!   NOW ADD THE NEW D.O. HERE:
!
      DO K=1,KC
        DO L=2,LA
          XMRM = WQV(L,K,19)*DTWQ*DZC(K)*HP(L)
          XDOTRN(L,K) = XDOTRN(L,K) + XMRM
          XDOALL(L,K) = XDOALL(L,K) + XMRM
        ENDDO
      ENDDO
!
! **  LOAD WQV INTO WQVO FOR REACTION CALCULATION
!
! J.S. ADD MACLGAE
!
       NMALG=0
       IF(IDNOTRVA.GT.0) NMALG=1
       DO NW=1,NWQV+NMALG
         DO K=1,KC
           DO L=2,LA
              WQVO(L,K,NW)=WQV(L,K,NW)
           ENDDO
         ENDDO
       ENDDO
!
! **  UPDATE WATER COLUMN KINETICS AND SEDIMENT MODEL
! **  OVER LONGER TIME INTERVALS THAN PHYSICAL TRANSPORT
! **  IF NWQKDPT .GT. 1
!
      NWQKCNT=NWQKCNT+1
      IF(NWQKCNT.EQ.NWQKDPT)THEN
        NWQKCNT=0
      ENDIF
!
! **    CALCULATE KINETIC SOURCES AND SINKS
!
        IF(ISCRAY.EQ.0)THEN
          TTMP=SECNDS(0.0)
         ELSE
          T1TMP=SECOND( )
          CALL TIMEF(WT1TMP)
        ENDIF
!
!
        IF(ISCRAY.EQ.0)THEN
          TWQKIN=TWQKIN+SECNDS(TTMP)
         ELSE
          T2TMP=SECOND( )
          CALL TIMEF(WT2TMP)
          TWQKIN=TWQKIN+T2TMP-T1TMP
          WTWQKIN=WTWQKIN+(WT2TMP-WT1TMP)*0.001
        ENDIF
!
! **    DIAGNOSE NEGATIVE CONCENTRATIONS
!
        CALL WWQNC
!
! **    WRITE TIME SERIES
!
        IF(ITNWQ.GE.IWQTSB .AND. ITNWQ.LE.IWQTSE)THEN
!           IF(MOD(ITNWQ,IWQTSDT).EQ.0) CALL WWQTS(TINDAY)
!           IF(MOD(ITNWQ,IWQTSDT).EQ.0) CALL WWQTS                        
           IF(MOD(ITNWQ,IWQTSDT).EQ.0) then                               
!             CALL WWQTS                                                  
!          CALL WWQTSBIN
        ENDIF
!
! **    WRITE SNAP SHOT  IWQSSB,IWQSSE,IWQSSDT
!
        IF(ITNWQ.GE.IWQSSB .AND. ITNWQ.LE.IWQSSE)THEN
           IF(MOD(ITNWQ,IWQSSDT).EQ.0) CALL WWQSNAPSHOT
        ENDIF
!
! **    CALL SEDIMENT MODEL
!
        IF(IWQBEN.EQ.1)THEN
!
          IF(ISCRAY.EQ.0)THEN
            TTMP=SECNDS(0.0)
           ELSE
            T1TMP=SECOND( )
            CALL TIMEF(WT1TMP)
          ENDIF
!
          CALL SMMBE
!
          IF(ISCRAY.EQ.0)THEN
            TWQSED=TWQSED+SECNDS(TTMP)
           ELSE
            T2TMP=SECOND( )
            CALL TIMEF(WT2TMP)
            TWQSED=TWQSED+T2TMP-T1TMP
            WTWQSED=WTWQSED+(WT2TMP-WT1TMP)*0.001
          ENDIF
!
          IF(ISMTS.GE.1)THEN
!
! **      WRITE SEDIMENT MODEL TIME SERIES
!
            IF(ITNWQ.GE.ISMTSB .AND. ITNWQ.LE.ISMTSE)THEN
!              IF(MOD(ITNWQ,ISMTSDT).EQ.0) CALL WSMTS(TINDAY)
              IF(MOD(ITNWQ,ISMTSDT).EQ.0) CALL WSMTS
            ENDIF
!
          ENDIF
!
! **      WRITE SEDIMENT MODEL FLUXES TO BINARY FILE:
!
          IF(ITNWQ.GE.ISMTSB .AND. ITNWQ.LE.ISMTSE)THEN
            CALL WSMTSBIN
          ENDIF
        ENDIF
!
      ENDIF
!
! **  ENDIF ON KINETIC AND SEDIMENT UPDATE
!
! **  INSERT TIME CALL
!
!     CALL TIME(HHMMSS)
!     OPEN(1,FILE='TFILE',STATUS='UNKNOWN',POSITION='APPEND')
!     WRITE(1,100) HHMMSS
!     CLOSE(1)
!
! **  WRITE RESTART FILES
!
!      IF(IWQRST.EQ.1) CALL WWQRST
!      IF(IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST
!
!XH      CLOSE(IWQONC)
!
  100 FORMAT('  TIME = ',A11,' HH.MM.SS.HH')
  600 FORMAT(' ITNWQ,IWQTSDT,ITMP,IWQTSB,TWQTSE = ',5I5)
!
      RETURN
      END