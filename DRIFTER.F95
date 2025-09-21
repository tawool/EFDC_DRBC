!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE DRIFTER
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
!**********************************************************************C
!
! **  SUBROUTINE DRIFTER CALCULATES THREE DIMENSIONAL TRAJECTORIES
! **  OF NEUTRALLY BUOYANT PARTICLES
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      CHARACTER*80 TITLE
      integer iterac
      TITLE='DRIFTER DATA'
      IF(JSPD.NE.1) GOTO 5
!
      OPEN(96,FILE='DRIFTER.OUT',STATUS='UNKNOWN')
      CLOSE(96,STATUS='DELETE')
      OPEN(96,FILE='DRIFTER.OUT',STATUS='UNKNOWN')
      OPEN(1,FILE='efdcdrifter.log',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
      JSPD=0
      WRITE(96,99) TITLE
      WRITE(96,98) NPD,KC
   98 FORMAT(2I10)
   99 FORMAT(A80)
  200 FORMAT(I10)
  201 FORMAT(2I4,1X,12F10.6)
  202 FORMAT(I5,I20,2I4,1X,12F10.6)   
!
!----------------------------------------------------------------------C
!
      DO NP=1,NPD
      I=NINT(RI(NP))
      J=NINT(RJ(NP))
      L=LIJ(I,J)
      ZS=(RK(NP)-0.5)/FLOAT(KC)
      DLOND=CDLON1+(CDLON2*RI(NP)+CDLON3)/60.
      DLATD=CDLAT1+(CDLAT2*RJ(NP)+CDLAT3)/60.
      WRITE(96,202)NP,N,I,J,ZS,HP(L)    
      ENDDO
!
      CLOSE(96)
!
!----------------------------------------------------------------------C
!
    5 CONTINUE
!
!----------------------------------------------------------------------C
!
! **  LOOP OVER THE NUMBER OF PARTICLES
!
      DO NP=1,NPD
!
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
      ELSE
        DELT=DTDYN
      ENDIF
!
      TOLD=0.
      iterac=1
!
  100 CONTINUE
!
      DTMIN=DELT
      TX=DELT
      TY=DELT
      TZ=DELT
      TNEW=TOLD+DTMIN
!
!----------------------------------------------------------------------C
!
      I=NINT(RI(NP))
      J=NINT(RJ(NP))
      K=NINT(RK(NP))
      L=LIJ(I,J)
!
      IF (L.EQ.0) GOTO 1234     
!
      LN=LNC(L)
      U1L=DTI*DXIU(L)*U1(L,K)
      U1LE=DTI*DXIU(L+1)*U1(L+1,K)
      V1L=DTI*DYIV(L)*V1(L,K)
      V1LN=DTI*DYIV(LN)*V1(LN,K)
      W1K=DTI*DZI*W1(L,K)/H1P(L)
      W1KM=DTI*DZI*W1(L,K-1)/H1P(L)
!
      UL=DTI*DXIU(L)*U(L,K)
      ULE=DTI*DXIU(L+1)*U(L+1,K)
      VL=DTI*DYIV(L)*V(L,K)
      VLN=DTI*DYIV(LN)*V(LN,K)
      WK=DTI*DZI*W(L,K)*HPI(L)
      WKM=DTI*DZI*W(L,K-1)*HPI(L)
!
      RIM=FLOAT(I)-0.5
      RJM=FLOAT(J)-0.5
      RKM=FLOAT(K)-0.5
      RIP=FLOAT(I)+0.5
      RJP=FLOAT(J)+0.5
      RKP=FLOAT(K)+0.5
!
      UOLD=U1L*(RI(NP)-RIP)*(TOLD-DT)-U1LE*(RI(NP)-RIM)*(TOLD-DT)-UL*(RI(NP)-RIP)*TOLD+ULE        &
           *(RI(NP)-RIM)*TOLD
!
      VOLD=V1L*(RJ(NP)-RJP)*(TOLD-DT)-V1LN*(RJ(NP)-RJM)*(TOLD-DT)-VL*(RJ(NP)-RJP)*TOLD+VLN        &
           *(RJ(NP)-RJM)*TOLD
!
      WOLD=W1KM*(RK(NP)-RKP)*(TOLD-DT)-W1K*(RK(NP)-RKM)*(TOLD-DT)-WKM*(RK(NP)-RKP)*TOLD+WK        &
           *(RK(NP)-RKM)*TOLD
!
!----------------------------------------------------------------------C
!
      TOP=0.5*DTMIN*((U1LE*RIM-U1L*RIP)*(TNEW-DT)+(UL*RIP-ULE*RIM)*TNEW+UOLD)+RI(NP)
!
      BOT=(U1LE-U1L)*(TNEW-DT)+(UL-ULE)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
!
      RINEW=TOP/BOT
!
!----------------------------------------------------------------------C
!
      TOP=0.5*DTMIN*((V1LN*RJM-V1L*RJP)*(TNEW-DT)+(VL*RJP-VLN*RJM)*TNEW+VOLD)+RJ(NP)
!
      BOT=(V1LN-V1L)*(TNEW-DT)+(VL-VLN)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
!
      RJNEW=TOP/BOT
!
!----------------------------------------------------------------------C
!
      TOP=0.5*DTMIN*((W1K*RKM-W1KM*RKP)*(TNEW-DT)+(WKM*RKP-WK*RKM)*TNEW+WOLD)+RK(NP)
!
      BOT=(W1K-W1KM)*(TNEW-DT)+(WKM-WK)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
!
      RKNEW=TOP/BOT
!
!----------------------------------------------------------------------C
!
! **  DETERMINE TIME INCREMENT IN X DIRECTION
!
      IF(RINEW.GT.RIP)THEN
!
      AAQ=ULE-U1LE
      BQ=(ULE*TOLD-U1LE*(TOLD-DT)+UOLD)
      CQ=-2.*(RIP-RI(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TX=TQ
      ELSE
       BQ=(ULE*TOLD-U1LE*(TOLD-DT)+UOLD)/AAQ
       CQ=-2.*(RIP-RI(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 10
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TX=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TX=TQ
       ENDIF
      ENDIF
!
      ENDIF
!
   10 CONTINUE
!
      IF(RINEW.LT.RIM)THEN
!
      AAQ=UL-U1L
      BQ=(UL*TOLD-U1L*(TOLD-DT)+UOLD)
      CQ=-2.*(RIM-RI(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TX=TQ
      ELSE
       BQ=(UL*TOLD-U1L*(TOLD-DT)+UOLD)/AAQ
       CQ=-2.*(RIM-RI(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 20
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TX=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TX=TQ
       ENDIF
      ENDIF
!
      ENDIF
!
   20 CONTINUE
!
!----------------------------------------------------------------------C
!
! **  DETERMINE TIME INCREMENT IN Y DIRECTION
!
      IF(RJNEW.GT.RJP)THEN
!
      AAQ=VLN-V1LN
      BQ=(VLN*TOLD-V1LN*(TOLD-DT)+VOLD)
      CQ=-2.*(RJP-RJ(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TY=TQ
      ELSE
       BQ=(VLN*TOLD-V1LN*(TOLD-DT)+VOLD)/AAQ
       CQ=-2.*(RJP-RJ(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 30
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TY=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TY=TQ
       ENDIF
      ENDIF
!
      ENDIF
!
   30 CONTINUE
!
      IF(RJNEW.LT.RJM)THEN
!
      AAQ=VL-V1L
      BQ=(VL*TOLD-V1L*(TOLD-DT)+VOLD)
      CQ=-2.*(RJM-RJ(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TY=TQ
      ELSE
       BQ=(VL*TOLD-V1L*(TOLD-DT)+VOLD)/AAQ
       CQ=-2.*(RJM-RJ(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 40
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TY=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TY=TQ
       ENDIF
      ENDIF
!
      ENDIF
!
   40 CONTINUE
!
!----------------------------------------------------------------------C
!
! **  DETERMINE TIME INCREMENT IN Z DIRECTION
!
      IF(RKNEW.GT.RKP)THEN
      AAQ=WK-W1K
      BQ=(WK*TOLD-W1K*(TOLD-DT)+WOLD)
      CQ=-2.*(RKP-RK(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TZ=TQ
      ELSE
       BQ=(WK*TOLD-W1K*(TOLD-DT)+WOLD)/AAQ
       CQ=-2.*(RKP-RK(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 50
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TZ=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TZ=TQ
       ENDIF
      ENDIF
!
      ENDIF
!
   50 CONTINUE
!
      IF(RKNEW.LT.RKM)THEN
!
      AAQ=WKM-W1KM
      BQ=(WKM*TOLD-W1KM*(TOLD-DT)+WOLD)
      CQ=-2.*(RKM-RK(NP))
!
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TZ=TQ
      ELSE
       BQ=(WKM*TOLD-W1KM*(TOLD-DT)+WOLD)/AAQ
       CQ=-2.*(RKM-RK(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 60
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TZ=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TZ=TQ
       ENDIF
      ENDIF
!
      ENDIF
!
   60 CONTINUE
!
!----------------------------------------------------------------------C
!
      IF(TX.LT.DTMIN) DTMIN=TX
      IF(TY.LT.DTMIN) DTMIN=TY
      IF(TZ.LT.DTMIN) DTMIN=TZ

      if(iterac.gt.10000) then
      open(1,file='efdcdrifter.log',POSITION='APPEND',STATUS='UNKNOWN')
	write(1,11)N,np,tnew,told,dtmin
	close(1)
	goto 1234
	end if
11        format('at time step N=',I8,'  drifter ',I3,' has more than 10000 iterations',3f10.4)


!
      TNEW=TOLD+DTMIN
!
!----------------------------------------------------------------------C
!
      TOP=0.5*DTMIN*((U1LE*RIM-U1L*RIP)*(TNEW-DT)+(UL*RIP-ULE*RIM)*TNEW+UOLD)+RI(NP)
!
      BOT=(U1LE-U1L)*(TNEW-DT)+(UL-ULE)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
!
      RI(NP)=TOP/BOT
!
!----------------------------------------------------------------------C
!
      TOP=0.5*DTMIN*((V1LN*RJM-V1L*RJP)*(TNEW-DT)+(VL*RJP-VLN*RJM)*TNEW+VOLD)+RJ(NP)
!
      BOT=(V1LN-V1L)*(TNEW-DT)+(VL-VLN)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
!
      RJ(NP)=TOP/BOT
!
!----------------------------------------------------------------------C
!
      TOP=0.5*DTMIN*((W1K*RKM-W1KM*RKP)*(TNEW-DT)+(WKM*RKP-WK*RKM)*TNEW+WOLD)+RK(NP)
!
      BOT=(W1K-W1KM)*(TNEW-DT)+(WKM-WK)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
!
      RK(NP)=TOP/BOT
!
      IF(ISPGNS.GE.1)THEN
        itest=nint(ri(np))
        jtest=nint(rj(np))
        DO NPNS=1,NPNSBP
          if(itest.eq.ISPNS(NPNS).and.jtest.eq.JSPNS(NPNS)-1) then
            ri(np)=ri(np)-float(ispns(NPNS))+float(inpns(NPNS))
            rj(np)=rj(np)-float(jspns(NPNS)-1)+float(jnpns(NPNS))
          end if
          if(itest.eq.INPNS(Npns).and.jtest.eq.JNPNS(Npns)+1) then
            ri(np)=ri(np)-float(inpns(NPNS))+float(ispns(NPNS))
            rj(np)=rj(np)-float(jnpns(NPNS)+1)+float(jspns(NPNS))
          end if
        END DO
      END IF

!
!----------------------------------------------------------------------C
!
      DELT=DELT-DTMIN
      TOLD=TNEW
!
      IF(DELT.GT.0.00001) then
      iterac=iterac+1
      GOTO 100
      end if
!
!----------------------------------------------------------------------C
!
1234  ENDDO
!
!**********************************************************************C
!
      IF(NCPD.EQ.NWPD)THEN
!
      OPEN(96,FILE='DRIFTER.OUT',POSITION='APPEND',STATUS='UNKNOWN')
!
      NCPD=1
      DO NP=1,NPD
      I=NINT(RI(NP))
      J=NINT(RJ(NP))
      L=LIJ(I,J)
          IF (L.EQ.0) THEN
            write(96,201)0,0,0.,0.,0.,0.     
            GOTO 5342
          END IF
      ZS=(RK(NP)-0.5)/FLOAT(KC)
      DLOND=CDLON1+(CDLON2*RI(NP)+CDLON3)/60.
      DLATD=CDLAT1+(CDLAT2*RJ(NP)+CDLAT3)/60.
      WRITE(96,202)NP,N,I,J,ZS,HP(L)        
5342      ENDDO
!
      ELSE
      NCPD=NCPD+1
!
      CLOSE(96)
!
      ENDIF
!
!**********************************************************************C
!
      RETURN
      END