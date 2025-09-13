!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
!TT      OLD SUBROUTINE RWQPSL
!
!**********************************************************************C
!
! READ IN TEMPORALLY VARYING POINT SOURCE INPUT (UNIT INWQPSL).
!: IN (KG/D) EXCEPT XPSQ(M^3/S), XPO2(G/M^3),
!                  XPTAM(KMOL/D), XPFCB(MPN/100ML).
!: IN GES, LOAD IS IN (G/D) EXCEPT TAM IN (MOL/D),
!                  FCB IN (MPN/100ML)*(M^3/S).
!: TO CONVERT KG/D TO G/D, (XW KG/D)*(10^3 G/KG) = (CONV1*XW G/D) WITH
!   CONV1=1.0E3.
!: FOR O2, (XPSQ M^3/S)*(XPO2 G/M^3)*(86400 S/D)=(CONV2*XPSQ*XPO2 G/D)
!   WITH CONV2=8.64E4.
!: FOR TAM, (XPTAM KMOL/D)*(10^3 MOL/KMOL) = (CONV1*XPTAM MOL/D).
!: FOR FCB, (XPFCB MPN/100ML)*(XPSQ M^3/S)*(86400 S/D) =
!   (CONV2*XPSQ*XPFCB (MPN/100ML)*M^3/D).
!
!**********************************************************************C
!
!TT     INCLUDE 'EFDC.PAR'
!TT      INCLUDE 'EFDC.CMN'
!
!TT      PARAMETER (CONV1=1.0E3,CONV2=8.64E4)
!TT      CHARACTER TITLE(3)*79, PSLCONT*3
!
!TT      OPEN(1,FILE=PSLFN,STATUS='UNKNOWN')
!TT      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
!
!TT      IF(IWQTPSL.EQ.0)THEN
!TT        READ(1,50) (TITLE(M),M=1,3)
!TT        WRITE(2,999)
!TT        WRITE(2,50) (TITLE(M),M=1,3)
!TT        READ(1,999)
!TT        READ(1,94) IWQPS
!TT        WRITE(2,83)'* NUMBER OF CELLS FOR POINT SOURCE INPUT      = ',
!TT     *    IWQPS
!TT      ENDIF
!TT      WRITE(2,60)'* POINT SOURCE INPUT AT ', IWQTPSL,
!TT     *  ' TH TC FROM MODEL START'
!
!TT      READ(1,999)
!TT      READ(1,50) (TITLE(M),M=1,3)
!TT      WRITE(2,50) (TITLE(M),M=1,3)
!TT      DO M=1,IWQPS
!TT        READ(1,94) I,J,K,WQPSQ(M),(WQWPSL(M,NW),NW=1,NWQV)
!TT        WRITE(2,94) I,J,K,WQPSQ(M),(WQWPSL(M,NW),NW=1,NWQV)
!TT        IF(IJCT(I,J).LT.1 .OR. IJCT(I,J).GT.8)THEN
!TT          PRINT*, 'I, J, K, M, WQPSQ = ', I,J,K,M,WQPSQ(M)
!TT          STOP 'ERROR!! INVALID (I,J) IN FILE 1'
!TT        ENDIF
!TT        L=LIJ(I,J)
!TT        IWQPSC(L,K)=M
!TT        DO NW=1,18
!TT          WQWPSL(M,NW) = WQWPSL(M,NW) * CONV1
!TT        ENDDO
!TT        WQTT = WQPSQ(M)*CONV2
!TT        WQWPSL(M,19) = WQWPSL(M,19) * WQTT
!TT        WQWPSL(M,20) = WQWPSL(M,20) * CONV1
!TT        WQWPSL(M,NWQV) = WQWPSL(M,NWQV) * WQTT
!TT      ENDDO
!
!TT      READ(1,52) IWQTPSL, PSLCONT
!TT      WRITE(2,52) IWQTPSL, PSLCONT
!TT      IF(PSLCONT.EQ.'END')THEN
!TT        CLOSE(1)
!TT        IWQPSL = 0
!TT      ENDIF
!
!TT      CLOSE(2)
!
!TT  999 FORMAT(1X)
!TT   50 FORMAT(A79)
!TT   52 FORMAT(I7, 1X, A3)
!TT   60 FORMAT(/, A24, I5, A24)
!TT   83 FORMAT(A48, I5)
!TT   94 FORMAT(3I5,1X,7F8.3, /, 16X,7F8.3, /, 16X, 8F8.3)
!
!TT      RETURN
!TT      END
!
!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE RWQPSL
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
! **  ALL LOADS ARE IN KG/DAY EXPECT COLIFORM IN MPN/DAY
! **  INTERNAL CONVERSION TO GM/DAY FOR FIRST 20 STATE VARIABLES
! **  FIX FOR TOTAL ACTIVE METAL AND COLIFORM IN FUTURE VERSIONS
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
       DIMENSION RLDTMP(21)
!
!**********************************************************************C
!
      IF(ITNWQ.GT.0) GOTO 1000
!
!**********************************************************************C
!
! **  READ IN LOADING SERIES FROM FILE 'WQPSL.INP'
!
!----------------------------------------------------------------------C
!
      IF( NPSTMSR.GE.1)THEN
        OPEN(1,FILE='WQPSL.INP',STATUS='UNKNOWN')
!
! **  SKIP OVER TITLE AND AND HEADER LINES
!
        DO IS=1,13
          READ(1,1)
        ENDDO
!
        DO NS=1,NPSTMSR
         MWQPTLT(NS)=1
         READ(1,*,IOSTAT=ISO)MWQPSR(NS),TCWQPSR(NS),TAWQPSR(NS),RMULADJ,ADDADJ
         IF(ISO.GT.0) GOTO 900
         RMULADJ=1000.*RMULADJ
         ADDADJ=ADDADJ
         DO M=1,MWQPSR(NS)
          READ(1,*,IOSTAT=ISO)TWQPSER(M,NS),(RLDTMP(NW),NW=1,7)
          IF(ISO.GT.0) GOTO 900
          READ(1,*,IOSTAT=ISO)(RLDTMP(NW),NW=8,14)
          IF(ISO.GT.0) GOTO 900
          READ(1,*,IOSTAT=ISO)(RLDTMP(NW),NW=15,21)
          IF(ISO.GT.0) GOTO 900
          TWQPSER(M,NS)=TWQPSER(M,NS)+TAWQPSR(NS)
!  FOR FECAL COLIFORM BACTERIA:
          WQPSSER(M,21,NS)=RMULADJ*RLDTMP(21)
          DO NW=1,20
            WQPSSER(M,NW,NS)=RMULADJ*RLDTMP(NW)
          ENDDO
          DO NW=1,21
            WQPSSER(M,NW,NS)=MAX(WQPSSER(M,NW,NS),0.0)
          ENDDO
         ENDDO
        ENDDO
!
        CLOSE(1)
      ENDIF
!
!     WRITE(6,602)
!
      GOTO 901
!
  900 CONTINUE
      WRITE(6,601)NS,M
      STOP
!
  901 CONTINUE
!
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQPS TIME SERIES, NSER,MDATA = ',2I5)
  602 FORMAT(' READ OF FILE WQPSL.INP SUCCESSFUL'/)
!
!**********************************************************************C
!
 1000 CONTINUE
!
!**********************************************************************C
!
! **  INITIALIZE NULL SERIES LOADING TO ZERO
!
      DO NW=1,NWQV
       WQPSSRT(NW,0)=0.
      ENDDO
!
!**********************************************************************C
!
! **  LOADING SERIES INTERPOLTATION
!
      DO NS=1,NPSTMSR
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCWQPSR(NS)
      ELSE
        TIME=TIMESEC/TCWQPSR(NS)
      ENDIF
!
       M1=MWQPTLT(NS)
  100  CONTINUE
       M2=M1+1
       IF(TIME.GT.TWQPSER(M2,NS))THEN
         M1=M2
         GOTO 100
        ELSE
         MWQPTLT(NS)=M1
       ENDIF
!
       TDIFF=TWQPSER(M2,NS)-TWQPSER(M1,NS)
       WTM1=(TWQPSER(M2,NS)-TIME)/TDIFF
       WTM2=(TIME-TWQPSER(M1,NS))/TDIFF
        DO NW=1,NWQV
         WQPSSRT(NW,NS)=WTM1*WQPSSER(M1,NW,NS)+WTM2*WQPSSER(M2,NW,NS)
        ENDDO
!      WRITE(6,6000)N,CSERTWQ(1,NS,NW),CSERTWQ(KC,NS,NW)
!
      ENDDO
!
!
      IF(ITNWQ.EQ.0)THEN
!
       OPEN(1,FILE='WQPSLT.DIA',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='WQPSLT.DIA',STATUS='UNKNOWN')
!
       WRITE(1,112)N,TIME
!
       DO NS=1,NPSTMSR
            WRITE(1,111)NS,(WQPSSRT(NW,NS),NW=1,NWQV)
       ENDDO
!
      CLOSE(1)
!
      ENDIF
!
!**********************************************************************C
!
! **  COMBINE CONSTANT AND TIME VARIABLE PS LOADS
!
! M.R. MORTON 02/20/1999
! MODIFIED SO MULTIPLE POINT SOURCES CAN BE ADDED TO ANY GRID CELL
! AND ANY LAYER (HAD TO CHANGE WQWPSL ARRAY FROM 2D TO 3D).
!MRM      DO NW=1,NWQV
!MRM       DO K=1,KC
!MRM        DO L=2,LA
!MRM         WQWPSL(IWQPSC(L,K),NW)=WQWPSLC(IWQPSC(L,K),NW)
!MRM     &                        +WQPSSRT(NW,IWQPSV(L,K))
!MRM        ENDDO
!MRM       ENDDO
!MRM      ENDDO
!
!MRM  INITIALIZE COMBINED PSL ARRAY TO ZERO:
!
      DO NW=1,NWQV
       DO K=1,KC
        DO L=2,LA
         WQWPSL(L,K,NW) = 0.0
        ENDDO
       ENDDO
      ENDDO
!
!MRM  ADD CONSTANT AND VARIABLE PSLS TO APPROPRIATE GRID CELLS:
!
      IF(ITNWQ.EQ.0)THEN
!
       OPEN(1,FILE='WQPSL.DIA',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='WQPSL.DIA',STATUS='UNKNOWN')
       WRITE(1,112)N,TIME
! 
      ENDIF
!
      DO NS=1,IWQPS
        L = LIJ(ICPSL(NS), JCPSL(NS))
        K = KCPSL(NS)
        ITMP = MVPSL(NS)
        IF(ITNWQ.EQ.0) WRITE(1,121)NS,L,ICPSL(NS),JCPSL(NS),K,ITMP
        IF(K.GE.1)THEN
          DO NW=1,NWQV
           WQWPSL(L,K,NW) = WQWPSL(L,K,NW) + WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP)
          ENDDO
         ELSE

! original even distribution of load was for sigma KC is the total layers
! Sen Bai modified to add even distribution of load when k=0 for GVC code
        if (IGRIDV .eq. 0) then  
          TMPVAL=1./FLOAT(KC)
          DO KK=1,KC
          DO NW=1,NWQV
           WQWPSL(L,KK,NW) = WQWPSL(L,KK,NW) + TMPVAL*( WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP) )
          ENDDO
          ENDDO
! when    IGRIDV=1 
        else 
          TMPVAL=1.0/ (FLOAT(KC)/GVCSCLP(L))
		TMPVAL=1./ (GVCSCLPI(L)*FLOAT(KC))
		DO KK=KC,KC-int(GVCSCLPI(L)*KC)+1,-1
          DO NW=1,NWQV
           WQWPSL(L,KK,NW) = WQWPSL(L,KK,NW) + TMPVAL*( WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP) )
          ENDDO
          ENDDO
        end if 


        ENDIF
      ENDDO
!
      IF(ITNWQ.EQ.0)THEN
!
!       OPEN(1,FILE='WQPSL.DIA',STATUS='UNKNOWN')
!       CLOSE(1,STATUS='DELETE')
!       OPEN(1,FILE='WQPSL.DIA',STATUS='UNKNOWN')
!
!       WRITE(1,112)N,TIME
!
       DO L=2,LA
        ITMP=IWQPSC(L,1)
        IF(ITMP.GT.0)THEN
          DO K=1,KC
            WRITE(1,110)ITMP,IL(L),JL(L),K,(WQWPSL(L,K,NW),NW=1,NWQV)
          ENDDO
        ENDIF
       ENDDO
!
!       DO K=1,KC
!        DO L=2,LA
!         ITMP=IWQPSC(L,K)
!         IF(ITMP.GT.0)THEN
!            WRITE(1,110)IL(L),JL(L),K,(WQWPSL(L,K,NW),NW=1,NWQV)
!         ENDIF
!        ENDDO
!       ENDDO
!
      CLOSE(1)
!
      ENDIF
!
  110 FORMAT(1X,4I4,2X,7E12.4,/,19X,7E12.4,/,19X,7E12.4)
  111 FORMAT(1X,I4,2X,7E12.4,/,7X,7E12.4,/,7X,7E12.4)
  112 FORMAT(' N, TIME = ', I10, F12.5/)
  121 FORMAT(' NS,L,I,J,K,ITMP = ', 6I5/)
!
!**********************************************************************C
!
      RETURN
      END