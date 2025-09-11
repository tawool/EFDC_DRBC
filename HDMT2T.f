C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE HDMT2T
C
C **  SUBROUTINE HDMT2T EXECUTES THE FULL HYDRODYNAMIC AND MASS
C **   TRANSPORT TIME INTERGATION USING A TWO TIME LEVEL SCHEME
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C
C 05/01/2002        John Hamrick       05/01/2002       John Hamrick
C  modified calls to calbal and budget subroutines
C  added calls to bal2t1, bal2t4, bal2t5
C 05/02/2002        John Hamrick       05/01/2002       John Hamrick
C  modified calculation of cell center bed stress (stored as QQ(l,0))
C  for cells have source/sinks
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION ISSBCP(LCM),WCOREW(LCM),WCORNS(LCM),LCORNER(LCM),
     &          LCORNWE(LCM),LCORNSN(LCM)
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        TTMP=SECOND( )
        CALL TIMEF(WTTMP)
      ENDIF
C
      ISTL=2
      FOURDPI=4./PI
      ISTL=2
      IS2TL=1
      ICALLTP=0
C
C**********************************************************************C
C
C **  SET FLAGS FOR CORNER CELL BED STRESS CORRECTIONS
C
      IF(ISCORTBC.GE.1) THEN
C
C **  SET FLAG FOR CELLS HAVING VOLUME SOURCE OR SINKS
C
      DO L=1,LC
	  ISSBCP(L)=0
	ENDDO
C
	DO L=2,LA
	  IF(RSSBCE(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCW(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCN(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCS(L).GT.1.5)ISSBCP(L)=1
	ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  REINITIALIZE VARIABLES
C
      DO L=2,LA
       H1P(L)=HP(L)
       H1U(L)=HU(L)
       H1UI(L)=HUI(L)
       H1V(L)=HV(L)
       H1VI(L)=HVI(L)
       UHDY1E(L)=UHDYE(L)
       VHDX1E(L)=VHDXE(L)
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
       U1(L,K)=U(L,K)
       V1(L,K)=V(L,K)
       UHDY1(L,K)=UHDY(L,K)
       VHDX1(L,K)=VHDX(L,K)
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  INITIALIZE COURANT NUMBER DIAGNOSTICS
C
      DO K=1,KC
      DO L=2,LA
       CFLUUU(L,K)=0.
       CFLVVV(L,K)=0.
       CFLWWW(L,K)=0.
       CFLCAC(L,K)=0.
      ENDDO
      ENDDO
C
C**********************************************************************C
C
      ILOGC=0
C
C**********************************************************************C
C
C **  CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING
C **  CALCULATE VELOCITY GRADIENTS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
C      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
       TMPVAL=1.+SUBO(L)+SVBO(L)+SUBO(L)*SVBO(L)
      H1C(L)=(HP(L)+SUBO(L)*HP(L-1)+SVBO(L)*HP(LS)
     &         +SUBO(L)*SVBO(L)*HP(LSW))/TMPVAL
      HMC(L)=0.25*(HMP(L)+HMP(L-1)+HMP(LS)+HMP(LSW))
      UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))
     &         +HP(L)*(U(L+1,1)+U(L,1)))*HVI(L)
      U1V(L)=0.25*(H1P(LS)*(U1(LSE,1)+U1(LS,1))
     &          +H1P(L)*(U1(L+1,1)+U1(L,1)))*H1VI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,1)+V(L-1,1))
     &         +HP(L)*(V(LN,1)+V(L,1)))*HUI(L)
      V1U(L)=0.25*(H1P(L-1)*(V1(LNW,1)+V1(L-1,1))
     &          +H1P(L)*(V1(LN,1)+V1(L,1)))*H1UI(L)
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
C
      IF(ISWAVE.EQ.1) CALL WAVEBL
      IF(ISWAVE.EQ.2) CALL WAVESXY
C
C**********************************************************************C
C
C **  FIRST CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
C
      CALL CALTBXY(ISTL,IS2TL)
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES
C
      IF(ISHDMF.GE.1) CALL CALHDMF
C
C**********************************************************************C
C
C **  CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N
C
C----------------------------------------------------------------------C
C
      N=-1
      CALL CALTSXY
C
       DO L=2,LA
       TBX1(L)=(AVCON1*H1UI(L)+STBX(L)*SQRT(V1U(L)*V1U(L)
     &        +U1(L,1)*U1(L,1)))*U1(L,1)
       TBY1(L)=(AVCON1*H1VI(L)+STBY(L)*SQRT(U1V(L)*U1V(L)
     &        +V1(L,1)*V1(L,1)))*V1(L,1)
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
C
C**********************************************************************C
C
C **  SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
C
      CALL CALTBXY(ISTL,IS2TL)
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE STRESSES
C
C----------------------------------------------------------------------C
C
       DO L=2,LA
       TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)
     &       +U(L,1)*U(L,1)))*U(L,1)
       TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)
     &       +V(L,1)*V(L,1)))*V(L,1)
       ENDDO
C
      N=0
      CALL CALTSXY
C
C----------------------------------------------------------------------C
C
C **  SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
C
      DO L=2,LA
        HDFUFX(L)=1.
        HDFUFY(L)=1.
        HDFUF(L)=1.
      ENDDO
C
      IF(ISBSDFUF.GE.1)THEN
cjh060305      HDFUFM=1.E-12
      HDFUFM=0.0
C
      DO L=2,LA
        LS=LSC(L)
        HDFUFX(L)=HDFUFM+G*SUB(L)*HU(L)*(BELV(L-1)-BELV(L))*DXIU(L)
        HDFUFY(L)=HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )-BELV(L))*DYIV(L)
      ENDDO
C
cjh060305      DO L=2,LA
cjh060305        HDFUFX(L)=TBX(L)/HDFUFX(L)
cjh060305        HDFUFY(L)=TBY(L)/HDFUFY(L)
cjh060305      ENDDO
C
cjh060305      DO L=2,LA
cjh060305        HDFUFX(L)=MAX(HDFUFX(L),-1.0)
cjh060305        HDFUFY(L)=MAX(HDFUFY(L),-1.0)
cjh060305      ENDDO
C
cjh060305      DO L=2,LA
cjh060305        HDFUFX(L)=MIN(HDFUFX(L),1.0)
cjh060305        HDFUFY(L)=MIN(HDFUFY(L),1.0)
cjh060305      ENDDO
C
      DO L=2,LA
	  IF(HDFUFX(L).GT.0.0)THEN
          HDFUFX(L)=TBX(L)/HDFUFX(L)
        ELSE
          HDFUFX(L)=1.0
	  ENDIF
	  IF(HDFUFY(L).GT.0.0)THEN
          HDFUFY(L)=TBY(L)/HDFUFY(L)
        ELSE
          HDFUFY(L)=1.0
	  ENDIF
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
C
C----------------------------------------------------------------------C
C
C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
C
      IF(ISWAVE.EQ.0)THEN
C
C----------------------------------------------------------------------c
C
       IF(ISCORTBC.EQ.0) THEN
C
C
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO
C

       DO L=2,LA
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       ENDDO
C
       DO L=2,LA
       QQ1(L,0 )=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX1(L))**2
     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY1(L))**2)
       QQ1(L,KC)=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX1(L))**2
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY1(L))**2)
       ENDDO
C
       DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
C
       DO L=2,LA
       QQ(L,0 )=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2
     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)
       ENDDO
C
       DO L=2,LA
	   TAUBSED(L)=QQ(L,0 )/CTURB2
	   TAUBSND(L)=QQ(L,0 )/CTURB2
       ENDDO
C
       ENDIF
C
C----------------------------------------------------------------------c
C
       IF(ISCORTBC.GE.1) THEN
C
       IF(ISCORTBCD.GE.1)THEN
         OPEN(1,FILE='ADJSTRESSE.OUT')
         CLOSE(1,STATUS='DELETE')
 	 ENDIF
C
       OPEN(1,FILE='TBCORINIT.OUT')
C
       DO L=2,LA
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       ENDDO
C
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO
C
       DO L=2,LA
         IF(ISSBCP(L).EQ.0)THEN
	     IF(SUB(L+1).LT.0.5) WCOREST(L)=FSCORTBCV(L)
	     IF(SUB(L).LT.0.5) WCORWST(L)=FSCORTBCV(L)
	     IF(SVB(LNC(L)).LT.0.5) WCORNTH(L)=FSCORTBCV(L)
	     IF(SVB(L).LT.0.5) WCORSTH(L)=FSCORTBCV(L)
	   ENDIF
	 ENDDO
C
       DO L=2,LA
	   WCOREW(L)=1./(WCOREST(L)+WCORWST(L))
	   WCORNS(L)=1./(WCORNTH(L)+WCORSTH(L))
	 ENDDO
C
       DO L=2,LA
         WCOREST(L)=WCOREST(L)*WCOREW(L)
	   WCORWST(L)=WCORWST(L)*WCOREW(L)
	   WCORNTH(L)=WCORNTH(L)*WCORNS(L)
	   WCORSTH(L)=WCORSTH(L)*WCORNS(L)
       ENDDO
C
       DO L=2,LA
C      QQ1(L,0 )=0.5*CTURB2*SQRT(
C     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX1(L))**2
C     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY1(L))**2)
       QQ1(L,0 )=CTURB2*SQRT(
     &       (RSSBCE(L)*WCOREST(L)*TVAR3E(L)
     &       +RSSBCW(L)*WCORWST(L)*TBX1(L))**2
     &      +(RSSBCN(L)*WCORNTH(L)*TVAR3N(L)
     &        +RSSBCS(L)*WCORSTH(L)*TBY1(L))**2)
       QQ1(L,KC)=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX1(L))**2
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY1(L))**2)
       ENDDO
C
       DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
C
       DO L=2,LA
C      QQ(L,0 )=0.5*CTURB2*SQRT(
C     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2
C     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2)
         QQ(L,0 )=CTURB2*SQRT(
     &       (RSSBCE(L)*WCOREST(L)*TVAR3E(L)
     &       +RSSBCW(L)*WCORWST(L)*TBX(L))**2
     &      +(RSSBCN(L)*WCORNTH(L)*TVAR3N(L)
     &       +RSSBCS(L)*WCORSTH(L)*TBY(L))**2)
         QQ(L,KC)=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)
       ENDDO
C
       DO L=2,LA
	   TAUBSED(L)=QQ(L,0 )/CTURB2
	   TAUBSND(L)=QQ(L,0 )/CTURB2
       ENDDO
C
       DO L=2,LA
         IF(WCORSTH(L).LT.0.49.OR.WCORSTH(L).GT.0.51)THEN
           IF(WCORWST(L).LT.0.49.OR.WCORWST(L).GT.0.51)THEN
             WRITE(1,3678)IL(L),JL(L),WCORWST(L),WCOREST(L),
     &       WCORSTH(L),WCORNTH(L)
	     ENDIF
	   ENDIF
       ENDDO
C
       CLOSE(1)
C
       ENDIF
C
C----------------------------------------------------------------------c
C
      ENDIF
C
C     ENDIF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
C
C----------------------------------------------------------------------C
C
C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
C
      IF(ISWAVE.GE.1)THEN
C
       DO L=2,LA
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       ENDDO
C
       DO L=2,LA
         TAUBC2=0.25*( (TVAR3E(L)+TBX1(L))**2
     &                        +(TVAR3N(L)+TBY1(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U1(L+1,1)+U1(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V1(LN,1)+V1(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     &        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ1(L,0 )=CTURB2*SQRT(TAUB2)
         QQ1(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX1(L))**2
     &                         +(TVAR3S(L)+TSY1(L))**2)
       ENDDO
C
       DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
C
       DO L=2,LA
         TAUBC2=0.25*( (TVAR3E(L)+TBX(L))**2
     &                        +(TVAR3N(L)+TBY(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     &        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ(L,0 )=CTURB2*SQRT(TAUB2)
         QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     &                         +(TVAR3S(L)+TSY(L))**2)
       ENDDO
C
      ENDIF
C
C     ENDIF
C
C**********************************************************************C
C
C **   SET SWITCHES FOR TWO TIME LEVEL INTEGRATION
C
       ISTL=2
       IS2TL=1
       DELT=DT
       DELTD2=DT/2.
       DZDDELT=DZ/DELT
       ROLD=0.
       RNEW=1.
C
C**********************************************************************C
C**********************************************************************C
C
C **  BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT
C **  CALCULATION
C
C **  SET CYCLE COUNTER AND CALL TIMER
C
      NTIMER=0
      N=0
      TIMEDAY=TCON*TBEGIN/86400.
C
C **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
C **  SUPPORTED ON OTHER SYSTEMS.
C
      CALL TIMELOG(N,TIMEDAY)
C     CALL DTIME (TARRAY)
C     WRITE(9,200)N, TARRAY(1),TARRAY(2)
C     CALL FLUSH(9)
  200 FORMAT('  N=',I10,5X,'USER TIME=',E12.4,5X,'SYSTEM TIME=',E12.4)
      NTIMER=1
C
      NINCRMT=1
C
C----------------------------------------------------------------------C
C
CXX DYNSTEP      DO 1000 N=1,NTS
C
 1001 CONTINUE
      IF(N.GE.NTS)GO TO 1000
C
      IF(ISDYNSTP.EQ.0)THEN
        N=N+1
        ETIMESEC=DT*FLOAT(N)
        ETIMEDAY=DT*FLOAT(N)/86400.
        TIMESEC=(DT*FLOAT(N)+TCON*TBEGIN)
        TIMEDAY=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      ELSE
        IF(IDRYTBP.EQ.0)THEN
          CALL CALSTEP
	  ELSE
	    CALL CALSTEPD
	  ENDIF
        DELT=DTDYN
        DELTD2=DTDYN/2.
        DZDDELT=DZ/DTDYN
        N=N+NINCRMT
        ETIMESEC=DT*FLOAT(N)
        ETIMEDAY=(DT*FLOAT(N))/86400.
        TIMESEC=(DT*FLOAT(N)+TCON*TBEGIN)
        TIMEDAY=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      ENDIF
C
      IF(ILOGC.EQ.NTSMMT)THEN
!        CLOSE(8,STATUS='DELETE')                             !hnr 7/27/2009
!        OPEN(8,FILE='EFDCLOG.OUT',STATUS='UNKNOWN')          !hnr 7/27/2009
        IF(ISDRY.GT.0)THEN
          OPEN(1,FILE='DRYWET.LOG',STATUS='UNKNOWN')
          CLOSE(1,STATUS='DELETE')
        ENDIF
        IF(ISCFL.EQ.1)THEN
          OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN')
          CLOSE(1,STATUS='DELETE')
        ENDIF
        ILOGC=0
      ENDIF
C
      IF(ISDYNSTP.EQ.0)THEN
        ILOGC=ILOGC+1
      ELSE
        ILOGC=ILOGC+NINCRMT
      ENDIF
C
      IF(N.LE.NLTS) SNLT=0.
      IF(N.GT.NLTS.AND.N.LE.NTTS)THEN
         NTMP1=N-NLTS
         NTMP2=NTTS-NLTS+1
         SNLT=FLOAT(NTMP1)/FLOAT(NTMP2)
      ENDIF
      IF(N.GT.NTTS) SNLT=1.
C
      IF(N.LE.NTSVB)THEN
       GP=GPO*(FLOAT(N)/FLOAT(NTSVB))
      ELSE
       GP=GPO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  INITIALIZE VOLUME, MASS, MOMENTUM, AND ENERGY BALANCE
C
C      IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
C         CALL CALBAL1
C         NTMP=MOD(N,2)
C         IF(NTMP.EQ.0)THEN
C           CALL CBALEV1
C          ELSE
C           CALL CBALOD1
C         ENDIF
C       ENDIF
C
C  ** INITIALIZE SEDIMENT BUDGET CALCULATION   (DLK 10/15)
C
C      IF(NCTBC.NE.NTSTBC.AND.ISSBAL.GE.1)THEN
C         CALL BUDGET1
C         NTMP=MOD(N,2)
C         IF(NTMP.EQ.0)THEN
C           CALL BUDGEV1
C          ELSE
C           CALL BUDGOD1
C         ENDIF
C       ENDIF
C
C **  INITIALIZE TWO-TIME LEVEL BALANCES
C
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T1
	  ENDIF
	ENDIF
C
C----------------------------------------------------------------------C
C
C **  REENTER HERE FOR TWO TIME LEVEL CORRECTION
C
  500 CONTINUE
C
C**********************************************************************C
C
C **  CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
      IF(KC.GT.1)THEN
        IF(ISQQ.EQ.1)THEN
	    IF(ISTOPT(0).EQ.0)CALL CALAVBOLD (ISTL)
	    IF(ISTOPT(0).GE.1)CALL CALAVB (ISTL)
        ENDIF
        IF(ISQQ.EQ.2) CALL CALAVB2 (ISTL)
      ENDIF
      IF(ISCRAY.EQ.0)THEN
        TAVB=TAVB+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TAVB=TAVB+T2TMP-T1TMP
        WTAVB=WTAVB+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
C
        IF(ISWAVE.EQ.1) CALL WAVEBL
        IF(ISWAVE.EQ.2) CALL WAVESXY
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT MOMENTUM EQUATION TERMS
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
      IF(IS2TIM.EQ.1) CALL CALEXP2T
      IF(IS2TIM.EQ.2) CALL CALIMP2T
C
      IF(ISCRAY.EQ.0)THEN
        TCEXP=TCEXP+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCEXP=TCEXP+T2TMP-T1TMP
        WTCEXP=WTCEXP+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS, CONCENTRATIONS,
C **  VEGETATION CHARACTERISTICS AND SURFACE ELEVATIONS
C
      CALL CALCSER (ISTL)
      CALL CALVEGSER (ISTL)
      CALL CALQVS (ISTL)
      PSERT(0)=0.
      IF(NPSER.GE.1) CALL CALPSER (ISTL)
C
C**********************************************************************C
C
C **  WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
C     SETUP CALL
C
C----------------------------------------------------------------------C
C
      IF(ISWSEDA.GT.0.OR.ISUVDA.GT.0) CALL PUVDASM(ISTL,1)
C
C**********************************************************************C
C
C **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
C **  MODE FORCING (S&M SOLUTION ONLY)
C
C----------------------------------------------------------------------C
C
      IF(ISCDMA.GE.3.AND.ISCDMA.LE.8)THEN
c      IF(ISTL.EQ.3)THEN
C
       DO L=2,LA
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
C
      CALL CALTSXY
C
       DO L=2,LA
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
C
      DO L=2,LA
        FXE(L)=FXE(L)+DT*SUB(L)*DYU(L)*(TSX(L)-TSX1(L))
        FYE(L)=FYE(L)+DT*SVB(L)*DXV(L)*(TSY(L)-TSY1(L))
      ENDDO
C
c      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
      IF(ISCHAN.EQ.0.AND.ISDRY.EQ.0) CALL CALPUV2T
      IF(ISCHAN.GE.1.OR.ISDRY.GE.1) CALL CALPUV2C
C
      IF(ISCRAY.EQ.0)THEN
        TPUV=TPUV+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TPUV=TPUV+T2TMP-T1TMP
        WTPUV=WTPUV+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  WRITE DIAGNOSTICS
C
C----------------------------------------------------------------------C
C
C **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
C **  SUPPORTED ON OTHER SYSTEMS.
C
      IF(ISLOG.GE.1)THEN
!      WRITE(8,17)N,ITER,RSQ,CFMAX,AVMAX,ABMIN,ABMAX,ABMIN     !hnr 7/27/2009
C     CALL FLUSH(8)
      ENDIF
C
C  17 FORMAT('  N,ITER,AVMA,AVMI,ABMA,ABMI',2I5,4(1X,F8.4))
   17 FORMAT('  N,ITER,RSQ,CFMAX,AVMA,AVMI,ABMA,ABMI',
     &        I10,I5,2E12.4,4(1X,F8.4))
C
      ERRMAX=MAX(ERRMAX,ERR)
      ERRMIN=MIN(ERRMIN,ERR)
      ITRMAX=MAX(ITRMAX,ITER)
      IRRMIN=MIN(ITRMIN,ITER)
C
C**********************************************************************C
C
C **  ADVANCE INTERNAL VARIABLES FOR THREE TIME LEVEL STEP
C
C----------------------------------------------------------------------C
C
C2T      IF(ISTL.EQ.3)THEN
C
      DO K=1,KC
      DO L=2,LA
      UHDY2(L,K)=UHDY1(L,K)
      UHDY1(L,K)=UHDY(L,K)
      VHDX2(L,K)=VHDX1(L,K)
      VHDX1(L,K)=VHDX(L,K)
      U2(L,K)=U1(L,K)
      V2(L,K)=V1(L,K)
      U1(L,K)=U(L,K)
      V1(L,K)=V(L,K)
      W2(L,K)=W1(L,K)
      W1(L,K)=W(L,K)
      ENDDO
      ENDDO
C
C2T      ENDIF
C
C**********************************************************************C
C
C **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
C **  MODE FORCING
C
C----------------------------------------------------------------------C
C
      IF(ISCDMA.LE.2.OR.ISCDMA.GE.9)THEN
c     IF(ISTL.EQ.3)THEN
C
       DO L=2,LA
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
C
      CALL CALTSXY
C
       DO L=2,LA
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
C
c     ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W
C
C----------------------------------------------------------------------C
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
      IF(KC.GT.1)THEN
        CALL CALUVW (ISTL,IS2TL)
      ELSE
         DO L=2,LA
          UHDY(L,1)=UHDYE(L)
          U(L,1)=UHDYE(L)*HUI(L)*DYIU(L)
          VHDX(L,1)=VHDXE(L)
          V(L,1)=VHDXE(L)*HVI(L)*DXIV(L)
          W(L,1)=0.
         ENDDO
        CALL CALUVW (ISTL,IS2TL)
      ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TUVW=TUVW+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TUVW=TUVW+T2TMP-T1TMP
        WTUVW=WTUVW+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
C     DIAGNOSTIC CALL
C
C----------------------------------------------------------------------C
C
      IF(ISWSEDA.GT.0.OR.ISUVDA.GT.0) CALL PUVDASM(ISTL,2)
C
C**********************************************************************C
C
C **  CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS
C **  AT TIME LEVEL (N+1)
C
C----------------------------------------------------------------------C
C
      CALL CALCONC (ISTL,IS2TL)
C
C----------------------------------------------------------------------C
C
      DO K=1,KB
      DO L=1,LC
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      ENDDO
      ENDDO
C
      DO NS=1,NSED
       DO K=1,KB
       DO L=1,LC
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
       ENDDO
       ENDDO
      ENDDO
C
      DO NS=1,NSND
       DO K=1,KB
       DO L=1,LC
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
       ENDDO
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=1,LC
        SEDT(L,K)=0.
        SNDT(L,K)=0.
       ENDDO
      ENDDO
C
      DO NS=1,NSED
       DO K=1,KC
        DO L=1,LC
         SEDT(L,K)=SEDT(L,K)+SED(L,K,NS)
        ENDDO
       ENDDO
      ENDDO
C
      DO NS=1,NSND
       DO K=1,KC
        DO L=1,LC
         SNDT(L,K)=SNDT(L,K)+SND(L,K,NS)
        ENDDO
       ENDDO
      ENDDO
C
CJH5/13/97      DO NT=1,NTOX
CJH5/13/97       DO K=1,KC
CJH5/13/97        DO L=1,LC
CJH5/13/97         TOXPFT(L,K,NT)=0.
CJH5/13/97        ENDDO
CJH5/13/97       ENDDO
CJH5/13/97      ENDDO
C
CJH5/13/97      DO NT=1,NTOX
CJH5/13/97       DO NS=1,NSED+NSND
CJH5/13/97        DO K=1,KC
CJH5/13/97         DO L=1,LC
CJH5/13/97          TOXPFT(L,K,NT)=TOXPFT(L,K,NT)+TOXPF(L,K,NS,NT)
CJH5/13/97         ENDDO
CJH5/13/97        ENDDO
CJH5/13/97       ENDDO
CJH5/13/97      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CHECK RANGE OF SALINITY AND DYE CONCENTRATION
C
      IF(ISMMC.EQ.1)THEN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(SAL(L,K).GT.SALMAX)THEN
       SALMAX=SAL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(SAL(L,K).LT.SALMIN)THEN
       SALMIN=SAL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,6001)N
      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(DYE(L,K).GT.SALMAX)THEN
       SALMAX=DYE(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(DYE(L,K).LT.SALMIN)THEN
       SALMIN=DYE(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN
C
!      WRITE(8,6004)SALMAX,IMAX,JMAX,KMAX     !hnr 7/27/2009
!      WRITE(8,6005)SALMIN,IMIN,JMIN,KMIN      !hnr 7/27/2009
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(SFL(L,K).GT.SALMAX)THEN
       SALMAX=SFL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(SFL(L,K).LT.SALMIN)THEN
       SALMIN=SFL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,6006)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6007)SALMIN,IMIN,JMIN,KMIN
C
      ENDIF
C
C
      IF(ISMMC.EQ.2)THEN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(TEM(L,K).GT.SALMAX)THEN
       SALMAX=TEM(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(TEM(L,K).LT.SALMIN)THEN
       SALMIN=TEM(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,6001)N
      WRITE(6,6008)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6009)SALMIN,IMIN,JMIN,KMIN
C
      ENDIF
C
 6001 FORMAT('  N=',I10)
 6002 FORMAT('  SALMAX=',F14.4,5X,'I,J,K=',(3I10))
 6003 FORMAT('  SALMIN=',F14.4,5X,'I,J,K=',(3I10))
 6004 FORMAT('  DYEMAX=',F14.4,5X,'I,J,K=',(3I10))
 6005 FORMAT('  DYEMIN=',F14.4,5X,'I,J,K=',(3I10))
 6006 FORMAT('  SFLMAX=',F14.4,5X,'I,J,K=',(3I10))
 6007 FORMAT('  SFLMIN=',F14.4,5X,'I,J,K=',(3I10))
 6008 FORMAT('  TEMMAX=',F14.4,5X,'I,J,K=',(3I10))
 6009 FORMAT('  TEMMIN=',F14.4,5X,'I,J,K=',(3I10))
C
C**********************************************************************C
C
C **  CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT
C **  CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOULBE TIME
C **  STEP TRANSPORT FIELD
C
C----------------------------------------------------------------------C
C
      IF(ISTRAN(4).GE.1.OR.ISTRAN(8).GE.1)THEN
C2T      NTMP=MOD(N,2)
C2T      IF(NTMP.EQ.0.AND.ISTL.EQ.3)THEN
C
C **  CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION
C
       DO L=2,LA
        HWQ(L)=HP(L)
        WWQ(L,0)=0.
       ENDDO
C
      DO K=1,KC
       DO L=2,LA
        UHDYWQ(L,K)=UHDY2(L,K)
        VHDXWQ(L,K)=VHDX2(L,K)
        UWQ(L,K)=U2(L,K)
        VWQ(L,K)=V2(L,K)
        WWQ(L,K)=W2(L,K)
       ENDDO
      ENDDO
C
C     ADD CHANNEL INTERACTIONS
C

      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        IF(MDCHTYP(NMD).EQ.1)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     &    +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
          HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))
     &    -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
        ENDIF
        IF(MDCHTYP(NMD).EQ.2)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     &    +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
          HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))
     &    -DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
        ENDIF
        IF(MDCHTYP(NMD).EQ.3)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     &    +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
     &    +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
          HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))
     &    -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
          HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))
     &    -DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
        ENDIF
        ENDDO
      ENDIF
C
C     END ADD CHANNEL INTERACTIONS
C
C      IF(ISTRAN(6).GE.1) CALL CALWQC(2)
      IF(ISTRAN(8).GE.1) CALL WQ3D
      IF(ISTRAN(4).GE.1) CALL CALSFT(2)
C
       DO L=2,LA
        H2WQ(L)=HWQ(L)
       ENDDO
C
C2T      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING
C **  AN EQUATION OF STATE
C
C      IF(ISTL.EQ.3)THEN
       DO K=1,KC
       DO L=2,LA
       B1(L,K)=B(L,K)
       ENDDO
       ENDDO
C      ENDIF
C
      IF(BSC.GT.1.E-6)THEN
        CALL CALBUOY
	ELSE
        DO K=1,KC
          DO L=2,LA
            B(L,K)=0.
          ENDDO
        ENDDO
	ENDIF
C
C      IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
C         CALL CALBAL4
C         NTMP=MOD(N,2)
C         IF(NTMP.EQ.0)THEN
C           CALL CBALEV4
C          ELSE
C           CALL CBALOD4
C         ENDIF
C      ENDIF
C
C **  CALL TWO-TIME LEVEL BALANCES
C
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T4
	  ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
C      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
       TMPVAL=1.+SUBO(L)+SVBO(L)+SUBO(L)*SVBO(L)
      H1C(L)=(HP(L)+SUBO(L)*HP(L-1)+SVBO(L)*HP(LS)
     &         +SUBO(L)*SVBO(L)*HP(LSW))/TMPVAL
      UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))
     &         +HP(L)*(U(L+1,1)+U(L,1)))*HVI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,1)+V(L-1,1))
     &         +HP(L)*(V(LN,1)+V(L,1)))*HUI(L)
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES
C **  AT TIME LEVEL (N)
C
C      IF(ISTL.NE.2.AND.ISHDMF.GE.1) CALL CALHDMF
      IF(ISHDMF.GE.1) CALL CALHDMF
C
C**********************************************************************C
C
C **  UPDATE BOTTOM STRESSES AND SURFACE AND BOTTOM TURBULENT
C **  INTENSITIES
C
C----------------------------------------------------------------------C
C
C      IF(ISTL.EQ.2)THEN
C
C      DO K=1,KC
C       DO L=2,LA
C        QQ(L,K)=SQRT(QQ(L,K)*QQ1(L,K))
C       ENDDO
C      ENDDO
C
C      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
      IF(ISCDMA.EQ.2)THEN
C
       DO L=2,LA
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,0)=QQ(L,0)+QQ1(L,0)
        QQ2(L,KC)=QQ(L,KC)+QQ1(L,KC)
        QQ1(L,0)=QQ(L,0)
        QQ1(L,KC)=QQ(L,KC)
       ENDDO
C
      ELSE
C
       DO L=2,LA
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,0)=QQ1(L,0)+QQ1(L,0)
        QQ2(L,KC)=QQ1(L,KC)+QQ1(L,KC)
        QQ1(L,0)=QQ(L,0)
        QQ1(L,KC)=QQ(L,KC)
       ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE BOTTOM STRESS AT LEVEL (N+1)
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
      CALL CALTBXY(ISTL,IS2TL)
C
       DO L=2,LA
        TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)
     &       +U(L,1)*U(L,1)))*U(L,1)
        TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)
     &       +V(L,1)*V(L,1)))*V(L,1)
       ENDDO
C
C**********************************************************************C
C
C **  SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
C
      IF(ISBSDFUF.GE.1)THEN
cjh060305      HDFUFM=1.E-12
      HDFUFM=0.0
C
      DO L=2,LA
        LS=LSC(L)
        HDFUFX(L)=HDFUFM+G*SUB(L)*HU(L)*(BELV(L-1)-BELV(L))*DXIU(L)
        HDFUFY(L)=HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )-BELV(L))*DYIV(L)
      ENDDO
C
cjh060305      DO L=2,LA
cjh060305        HDFUFX(L)=TBX(L)/HDFUFX(L)
cjh060305        HDFUFY(L)=TBY(L)/HDFUFY(L)
cjh060305      ENDDO
C
cjh060305      DO L=2,LA
cjh060305        HDFUFX(L)=MAX(HDFUFX(L),-1.0)
cjh060305        HDFUFY(L)=MAX(HDFUFY(L),-1.0)
cjh060305      ENDDO
C
cjh060305      DO L=2,LA
cjh060305        HDFUFX(L)=MIN(HDFUFX(L),1.0)
cjh060305        HDFUFY(L)=MIN(HDFUFY(L),1.0)
cjh060305      ENDDO
C
      DO L=2,LA
	  IF(HDFUFX(L).GT.0.0)THEN
          HDFUFX(L)=TBX(L)/HDFUFX(L)
        ELSE
          HDFUFX(L)=1.0
	  ENDIF
	  IF(HDFUFY(L).GT.0.0)THEN
          HDFUFY(L)=TBY(L)/HDFUFY(L)
        ELSE
          HDFUFY(L)=1.0
	  ENDIF
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
C
C----------------------------------------------------------------------C
C
C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
C
      IF(ISWAVE.EQ.0)THEN
C
C----------------------------------------------------------------------c
C
       IF(ISCORTBC.EQ.0) THEN
C
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO
C
       DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
C
       DO L=2,LA
       QQ(L,0 )=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2
     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)
       ENDDO
C
       ENDIF
C
C----------------------------------------------------------------------c
C
       IF(ISCORTBC.GE.1) THEN
C
       IF(ISCORTBCD.GE.1)THEN
         NTMPVAL=MOD(N,NTSPTC)
	   IF(NTMPVAL.EQ.0)THEN
           OPEN(1,FILE='ADJSTRESSE.OUT',ACCESS='APPEND')
 	   ENDIF
 	 ENDIF
C
       DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
C
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO
C
       DO L=2,LA
         IF(ISSBCP(L).EQ.0)THEN
	     IF(SUB(L+1).LT.0.5)WCOREST(L)=FSCORTBCV(L)
	     IF(SUB(L).LT.0.5)WCORWST(L)=FSCORTBCV(L)
	     IF(SVB(LNC(L)).LT.0.5)WCORNTH(L)=FSCORTBCV(L)
	     IF(SVB(L).LT.0.5)WCORSTH(L)=FSCORTBCV(L)
	   ENDIF
	 ENDDO
C
       DO L=2,LA
	   WCOREW(L)=1./(WCOREST(L)+WCORWST(L))
	   WCORNS(L)=1./(WCORNTH(L)+WCORSTH(L))
	 ENDDO
C
       DO L=2,LA
         WCOREST(L)=WCOREST(L)*WCOREW(L)
	   WCORWST(L)=WCORWST(L)*WCOREW(L)
	   WCORNTH(L)=WCORNTH(L)*WCORNS(L)
	   WCORSTH(L)=WCORSTH(L)*WCORNS(L)
       ENDDO
C
       DO L=2,LA
C
C       QQ(L,0 )=0.5*CTURB2*SQRT(
C     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2
C     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2)
       QQ(L,0 )=CTURB2*SQRT(
     &       (RSSBCE(L)*WCOREST(L)*TVAR3E(L)
     &       +RSSBCW(L)*WCORWST(L)*TBX(L))**2
     &      +(RSSBCN(L)*WCORNTH(L)*TVAR3N(L)
     &       +RSSBCS(L)*WCORSTH(L)*TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT(
     &           (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)
       ENDDO
C
      IF(ISCORTBCD.GE.1.AND.NTMPVAL.EQ.0)THEN
C
      DO L=2,LA
        LCORNER(L)=0
	  KCORNER=0
	  IF(WCORWST(L).GT.0.505)THEN
	    KCORNER=KCORNER+1
	    LCORNWE(L)=L-1
	  ENDIF
	  IF(WCOREST(L).GT.0.505)THEN
	    KCORNER=KCORNER+1
	    LCORNWE(L)=L+1
	  ENDIF
	  IF(WCORNTH(L).GT.0.505)THEN
	    KCORNER=KCORNER+1
	    LCORNSN(L)=LNC(L)
	  ENDIF
	  IF(WCORSTH(L).GT.0.505)THEN
	    KCORNER=KCORNER+1
	    LCORNSN(L)=LSC(L)
	  ENDIF
	  IF(KCORNER.EQ.2)LCORNER(L)=1
      ENDDO
C
      NCORCELLS=0
	DO L=2,LA
	  NCORCELLS=NCORCELLS+LCORNER(L)
	ENDDO
C
      WRITE(1,3675)TIMEDAY,NCORCELLS
C
      DO L=2,LA
        IF(LMASKDRY(L))THEN
	  IF(LCORNER(L).EQ.1)THEN
	    LWE=LCORNWE(L)
	    LSN=LCORNSN(L)
          TAUTMP=QQ(L,0)/CTURB2
          TAUTMPWE=QQ(LWE,0)/CTURB2
          TAUTMPSN=QQ(LSN,0)/CTURB2
 	    WRITE(1,3677)IL(L),JL(L),TAUTMP,TAUBSND(L),
     &                 TAUBSED(L)
 	    WRITE(1,3676)IL(LWE),JL(LWE),TAUTMPWE,TAUBSND(LWE),
     &                 TAUBSED(LWE)
 	    WRITE(1,3676)IL(LSN),JL(LSN),TAUTMPSN,TAUBSND(LSN),
     &                 TAUBSED(LSN)
C          WRITE(1,3678)IL(L),JL(L),SUB(L),SUB(L+1),SVB(L),
C     &                                   SVB(LNC(L))
C          WRITE(1,3679)RSSBCW(L),RSSBCE(L),RSSBCS(L),RSSBCN(L)
C          WRITE(1,3679)WTMPWST,WTMPEST,WTMPSTH,WTMPNTH
C          WRITE(1,3680)TBX(L),TVAR3E(L),TBY(L),TVAR3N(L),
C     &             TAUTMP,TAUBSND(L)
C	    DO NX=1,NSND
C	      BLTMPVAL=QSBDLDX(L+1,NX)-QSBDLDX(L,NX)
C     &            +QSBDLDY(LNC(L),NX)-QSBDLDY(L,NX)
C            WRITE(1,3681)QSBDLDX(L,NX),QSBDLDX(L+1,NX),QSBDLDY(L,NX),
C     &             QSBDLDY(LNC(L),NX),BLTMPVAL,TAUR(NX+1)
C	    ENDDO
	  ENDIF
	  ENDIF
      ENDDO
C
      ENDIF

      CLOSE(1)
C
      ENDIF
C
C----------------------------------------------------------------------c
C
      ENDIF
C
 3678 FORMAT(2I6,4F13.3)
 3679 FORMAT(12x,4F13.3)
 3680 FORMAT(12x,6F13.5)
 3681 FORMAT(12X,5E13.4,F13.5)
 3677 FORMAT('CORNER',2I5,5E14.5)
 3676 FORMAT(6X,2I5,5E14.5)
 3675 FORMAT(F11.3,I6,' TIME IN DAYS AND NUMBER OF CORNERS')
C
C     ENDIF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
C
C----------------------------------------------------------------------C
C
C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
C
      IF(ISWAVE.GE.1)THEN
C
       DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
C
       DO L=2,LA
         TAUBC2=0.25*( (TVAR3E(L)+TBX(L))**2
     &                        +(TVAR3N(L)+TBY(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     &        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ(L,0 )=CTURB2*SQRT(TAUB2)
         QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     &                         +(TVAR3S(L)+TSY(L))**2)
       ENDDO
C
      ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TTBXY=TTBXY+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TTBXY=TTBXY+T2TMP-T1TMP
        WTTBXY=WTTBXY+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C     ENDIF
C
C**********************************************************************C
C
C **  CALCULATE TURBULENT INTENSITY SQUARED
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
      IF(KC.GT.1)THEN
        IF(ISQQ.EQ.1)THEN
	    IF(ISTOPT(0).EQ.0)CALL CALQQ2TOLD (ISTL)
	    IF(ISTOPT(0).GE.1)CALL CALQQ2T (ISTL)
        ENDIF
        IF(ISQQ.EQ.2) CALL CALQQ2 (ISTL)
      ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TQQQ=TQQQ+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TQQQ=TQQQ+T2TMP-T1TMP
        WTQQQ=WTQQQ+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE MEAN MASS TRANSPORT FIELD
C
      IF(ISSSMMT.NE.2)THEN
        IF(ISICM.GE.1)THEN
          NTMP=MOD(N,2)
          IF(ISTL.EQ.3.AND.NTMP.EQ.0) CALL CALMMT
        ENDIF
      ENDIF
C
C      IF(ISSSMMT.NE.2) CALL CALMMT
C
C**********************************************************************C
C
C **  HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
C **  IF NCTBC EQ NTSTBC APPLY TRAPEZOIDAL CORRECTION
C
C----------------------------------------------------------------------C
C
C2T      IF(NCTBC.EQ.NTSTBC)THEN
C2T       NCTBC=0
C2T       ISTL=2
C2T       DELT=DT
C2T       DELTD2=0.5*DT
C2T       DZDDELT=DZ/DELT
C2T       ROLD=0.5
C2T       RNEW=0.5
C2T       GOTO 500
C2T      ELSE
C2T       NCTBC=NCTBC+1
C2T       ISTL=3
C2T       DELT=DT2
C2T       DELTD2=DT
C2T       DZDDELT=DZ/DELT
C2T       ROLD=0.
C2T       RNEW=1.
C2T      ENDIF
C
C **  WRITE TO TIME SERIES FILES
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C
CDYN      IF(ISTMSR.GE.1)THEN
CDYN        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
CDYN          IF(NCTMSR.EQ.NWTMSR)THEN
CDYN            CALL TMSR
CDYN            ICALLTP=1
CDYN            NCTMSR=1
CDYN           ELSE
CDYN            NCTMSR=NCTMSR+1
CDYN          ENDIF
CDYN        ENDIF
CDYN      ENDIF
C
C
      IF(ISTMSR.GE.1)THEN
c        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
          IF(NCTMSR.GE.NWTMSR)THEN
            CALL TMSR
            NDIFF=NWTMSR-NCTMSR
            ICALLTP=1
            NCTMSR=NINCRMT+NDIFF
           ELSE
            NCTMSR=NCTMSR+NINCRMT
          ENDIF
c        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  WRITE TO DUMP FILES
C
      IF(ISDUMP.GE.1)THEN                                                           !HNR_GHD 7/2019 BMD2                     
        IF(TIME.GE.TSDUMP.AND.TIME.LE.TEDUMP)THEN            !hnr
c          IF(NCDUMP.EQ.NSDUMP)THEN                           !hnr
          IF(NCDUMP.GE.NSDUMP)THEN                           !hnr
            CALL DUMP                                        !hnr
c            NDIFF=NSDUMP-NCDUMP                              !hnr
            ICALLTP=1                                        !hnr
            NCDUMP=1                                         !hnr
c            NCDUMP=NINCRMT+NDIFF                             !hnr
           ELSE                                              !hnr
            NCDUMP=NCDUMP+1                                  !hnr
c            NCDUMP=NCDUMP+NINCRMT                            !hnr
          ENDIF                                              !hnr
        ENDIF                                                !hnr
      ENDIF                                                  !hnr
C
C**********************************************************************C
C
C **  WRITE TO bmd2 OUTPUT FILE                                                 !HNR_GHD 7/2019 BMD2    !HNR_GHD 6/2022 remove bmd file keep only bmd2                 
C                                                                    
      IF (bmdflag.GT.0) THEN                                          
        IF (TIME.GE.TIBMD.AND.TIME.LE.TEBMD) THEN                  
          IF (NCBMD.EQ.NSBMD) THEN                                 
              CALL BMD2_DUMP
            NCBMD=1                                                 
          ELSE                                                      
            NCBMD=NCBMD+1                                          
          END IF                                                     
        END IF                                                       
      END IF                                                         
C
C**********************************************************************C
C
C **  WRITE TOXIC FLUX ACCUMULATION FILES
C
      IF(ISTOXALL.EQ.1)THEN
        MSTOXALL=MSTOXALL+1
        IF(N.GE.MSTOXALL)THEN
          CALL TOXALL
          MSTOXALL=MSTOXALL+(NTSPTC/NSTOXALL)
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C ** WRITE BOTTOM VELOCITIES TO FILE FOR AJ MINE SEDIMENT BED MODEL
C
C     (DLK - ADDED 3/25/97)
      IF(IHYDOUT.GE.1)THEN
        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
          IF(NCHYDOUT.EQ.NWHYDOUT)THEN
            CALL HYDOUT
            NCHYDOUT=1
           ELSE
            NCHYDOUT=NCHYDOUT+1
          ENDIF
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  OUTPUT ZERO DIMENSION VOLUME BALANCE
C
C----------------------------------------------------------------------C
C
C      IF(ISDRY.GE.1.AND.ISDRY.LT.98)THEN
C      IF(ICALLTP.EQ.1)THEN
C        OPEN(1,FILE='ZVOLBAL.OUT',POSITION='APPEND',STATUS='UNKNOWN')
C        DO LS=1,LORMAX
C        IF(VOLZERD.GE.VOLSEL(LS).AND.VOLZERD.LT.VOLSEL(LS+1))THEN
C           WTM=VOLSEL(LS+1)-VOLZERD
C           WTMP=VOLZERD-VOLSEL(LS)
C           DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
C           WTM=WTM/DELVOL
C           WTMP=WTMP/DELVOL
C           SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
C           ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
C        ENDIF
C        ENDDO
C        IF(ISDYNSTP.EQ.0)THEN
C          TIME=DT*FLOAT(N)+TCON*TBEGIN
C          TIME=TIME/TCTMSR
C        ELSE
C          TIME=TIMESEC/TCTMSR
C        ENDIF
C        WRITE(1,5304) TIME,SELZERD,ASFZERD,VOLZERD,VETZERD
C        CLOSE(1)
C      ENDIF
C      ENDIF
C      ICALLTP=0
C
 5304 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
C
C**********************************************************************C
C
C **  WRITE VERTICAL SCALAR FIELD PROFILES
C
      IF(ISVSFP.EQ.1)THEN
        IF(N.GE.NBVSFP.AND.N.LE.NSVSFP)THEN
          CALL VSFP
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE MEAN MASS TRANSPORT FIELD
C
      IF(ISSSMMT.NE.2)THEN
        IF(ISICM.EQ.0) CALL CALMMT
      ENDIF
C
C      IF(ISSSMMT.NE.2) CALL CALMMT
C
C**********************************************************************C
C
C **  ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES
C
      IF(ISPD.EQ.1)THEN
        IF(N.GE.NPDRT) CALL DRIFTER
      ENDIF
      IF(ISLRPD.GE.1)THEN
        IF(ISCRAY.EQ.0)THEN
          T1TMP=SECNDS(0.0)
         ELSE
          T1TMP=SECOND( )
          CALL TIMEF(WT1TMP)
        ENDIF
        IF(ISLRPD.LE.2)THEN
          IF(N.GE.NLRPDRT(1)) CALL LAGRES
        ENDIF
        IF(ISLRPD.GE.3)THEN
          IF(N.GE.NLRPDRT(1)) CALL GLMRES
        ENDIF
        IF(ISCRAY.EQ.0)THEN
          TLRPD=TLRPD+SECNDS(T1TMP)
         ELSE
          T2TMP=SECOND( )
          CALL TIMEF(WT2TMP)
          TLRPD=TLRPD+T2TMP-T1TMP
          WTLRPD=WTLRPD+(WT2TMP-WT1TMP)*0.001
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE VOLUME MASS, MOMENTUM AND ENERGY BALANCES
C
C      IF(ISBAL.GE.1)THEN
C         CALL CALBAL5
C         NTMP=MOD(N,2)
C         IF(NTMP.EQ.0)THEN
C           CALL CBALEV5
C          ELSE
C           CALL CBALOD5
C         ENDIF
C       ENDIF
C
C   SEDIMENT BUDGET CALCULATION     (DLK 10/15)
C
C       IF(ISSBAL.GE.1)THEN
C       CALL BUDGET5
C       ENDIF
C       NTMP=MOD(N,2)
C       IF(NTMP.EQ.0)THEN
C         CALL BUDGEV5
C        ELSE
C         CALL BUDGOD5
C       ENDIF
C
C **  CALL TWO-TIME LEVEL BALANCES
C
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T5
	  ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  PERFORM AN M2 TIDE HARMONIC ANALYSIS EVERY 2 M2 PERIODS
C
      IF(ISHTA.EQ.1) CALL CALHTA
C
C**********************************************************************C
C
C **  CALCULATE DISPERSION COEFFICIENTS
C
C     IF(N.GE.NDISP)THEN
      IF(N.GE.NDISP.AND.NCTBC.EQ.1)THEN
       IF(ISDISP.EQ.2) CALL CALDISP2
       IF(ISDISP.EQ.3) CALL CALDISP3
      ENDIF
C
C**********************************************************************C
C
C **  PERFORM LEAST SQUARES HARMONIC ANALYSIS AT SELECTED LOCATIONS
C
      IF(ISLSHA.EQ.1.AND.N.EQ.NCLSHA)THEN
       CALL LSQHARM
       NCLSHA=NCLSHA+(NTSPTC/24)
      ENDIF
C
C**********************************************************************C
C
C **  PRINT INTERMEDIATE RESULTS
C
C----------------------------------------------------------------------C
C
      IF(NPRINT .EQ. NTSPP)THEN
       NPRINT=1
       CALL OUTPUT1
      ELSE
       NPRINT=NPRINT+1
      ENDIF
C
C**********************************************************************C
C
C **  WRITE TO TIME VARYING GRAPHICS FILES
C
C----------------------------------------------------------------------C
C
CDYN      IF(N.EQ.NCPPH.AND.ISPPH.EQ.1)THEN
C
      IF(N.GE.NCPPH.AND.ISPPH.GE.1)THEN
       CALL SURFPLT
       NCPPH=NCPPH+(NTSPTC/NPPPH)
      ENDIF
C
C
C----------------------------------------------------------------------C
C
CDYN      IF(N.EQ.NCBPH.AND.ISBPH.EQ.1)THEN
C
      IF(N.GE.NCBPH.AND.ISBPH.GE.1)THEN
c	IF(ISBEXP.EQ.0)THEN
       CALL BEDPLTH
       NCBPH=NCBPH+(NTSPTC/NPBPH)
c      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
CDYN      IF(N.EQ.NCVPH.AND.ISVPH.GE.1)THEN
C
      IPLTTMP=0
      IF(ISVPH.EQ.1.OR.ISVPH.EQ.2)IPLTTMP=1
      IF(N.GE.NCVPH.AND.IPLTTMP.EQ.1)THEN
       CALL VELPLTH
       NCVPH=NCVPH+(NTSPTC/NPVPH)
      ENDIF
C
C----------------------------------------------------------------------C
C
CDYN      IF(N.EQ.NCVPV.AND.ISVPV.GE.1)THEN
C
      IF(N.GE.NCVPV.AND.ISVPV.GE.1)THEN
       CALL VELPLTV
       NCVPV=NCVPV+(NTSPTC/NPVPV)
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC
        TVAR1S(L,K)=TOX(L,K,1)
       ENDDO
      ENDDO
C
CDYN      IF(N.EQ.NCSPH(1).AND.ISSPH(1).EQ.1)THEN
C
C      IF(N.GE.NCSPH(1).AND.ISSPH(1).EQ.1)THEN
C       IF(ISTRAN(1).GE.1) CALL SALPLTH (1,SAL)
C       NCSPH(1)=NCSPH(1)+(NTSPTC/NPSPH(1))
C      ENDIF
C
C      IF(N.GE.NCSPH(2).AND.ISSPH(2).EQ.1)THEN
C       IF(ISTRAN(2).GE.1) CALL SALPLTH (2,TEM)
C       NCSPH(2)=NCSPH(2)+(NTSPTC/NPSPH(2))
C      ENDIF
C
C      IF(N.GE.NCSPH(3).AND.ISSPH(3).EQ.1)THEN
C       IF(ISTRAN(3).GE.1) CALL SALPLTH (3,DYE)
C       NCSPH(3)=NCSPH(3)+(NTSPTC/NPSPH(3))
C      ENDIF
C
C      IF(N.GE.NCSPH(4).AND.ISSPH(4).EQ.1)THEN
C       IF(ISTRAN(4).GE.1) CALL SALPLTH (4,SFL)
C       NCSPH(4)=NCSPH(4)+(NTSPTC/NPSPH(4))
C      ENDIF
C
C      IF(N.GE.NCSPH(5).AND.ISSPH(5).EQ.1)THEN
C       IF(ISTRAN(5).GE.1) CALL SALPLTH (5,TVAR1S)
C       NCSPH(5)=NCSPH(5)+(NTSPTC/NPSPH(5))
C      ENDIF
C
C      IF(N.GE.NCSPH(6).AND.ISSPH(6).EQ.1)THEN
C       IF(ISTRAN(6).GE.1) CALL SALPLTH (6,SEDT)
C       NCSPH(6)=NCSPH(6)+(NTSPTC/NPSPH(6))
C      ENDIF
C
C      IF(N.GE.NCSPH(7).AND.ISSPH(7).EQ.1)THEN
C       IF(ISTRAN(7).GE.1) CALL SALPLTH (7,SNDT)
C       NCSPH(7)=NCSPH(7)+(NTSPTC/NPSPH(7))
C      ENDIF
C
C
      IPLTTMP=0
      IF(ISSPH(1).EQ.1.OR.ISSPH(1).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(1).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(1).GE.1) CALL SALPLTH (1,SAL)
       NCSPH(1)=NCSPH(1)+(NTSPTC/NPSPH(1))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(2).EQ.1.OR.ISSPH(2).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(2).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(2).GE.1) CALL SALPLTH (2,TEM)
       NCSPH(2)=NCSPH(2)+(NTSPTC/NPSPH(2))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(3).EQ.1.OR.ISSPH(3).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(3).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(3).GE.1) CALL SALPLTH (3,DYE)
       NCSPH(3)=NCSPH(3)+(NTSPTC/NPSPH(3))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(4).EQ.1.OR.ISSPH(4).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(4).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(4).GE.1) CALL SALPLTH (4,SFL)
       NCSPH(4)=NCSPH(4)+(NTSPTC/NPSPH(4))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(5).EQ.1.OR.ISSPH(5).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(5).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(5).GE.1) CALL SALPLTH (5,TVAR1S)
       NCSPH(5)=NCSPH(5)+(NTSPTC/NPSPH(5))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(6).EQ.1.OR.ISSPH(6).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(6).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(6).GE.1) CALL SALPLTH (6,SEDT)
       NCSPH(6)=NCSPH(6)+(NTSPTC/NPSPH(6))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(7).EQ.1.OR.ISSPH(7).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(7).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(7).GE.1) CALL SALPLTH (7,SNDT)
       NCSPH(7)=NCSPH(7)+(NTSPTC/NPSPH(7))
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO ITMP=1,7
      IF(N.GE.NCSPV(ITMP).AND.ISSPV(ITMP).GE.1)THEN
       CALL SALPLTV(ITMP)
       NCSPV(ITMP)=NCSPV(ITMP)+(NTSPTC/NPSPV(ITMP))
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  WRITE EFDC EXPLORER FORMAT OUTPUT
C
      IF(ISSPH(8).EQ.1.OR.ISBEXP.EQ.1)THEN
      IF(N.GE.NCSPH(8))THEN
       CALL EEXPOUT(0)
       NCSPH(8)=NCSPH(8)+(NTSPTC/NPSPH(8))
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  WRITE TO TIME VARYING 3D HDF GRAPHICS FILES
C
C----------------------------------------------------------------------C
C
      IF(N.EQ.NC3DO.AND.IS3DO.EQ.1)THEN
       CALL OUT3D
       NC3DO=NC3DO+(NTSPTC/NP3DO)
      ENDIF
C
C**********************************************************************C
C
C **  WRITE RESTART FILE EVERY ISRESTO M2 TIDAL CYCLES
C
      IF(ISRESTO.GE.1)THEN
        NRESTO=ISRESTO*NTSPTC
        ISSREST=MOD(N,NRESTO)
        IF(ISSREST.EQ.0)THEN
          CALL RESTOUT(0)
          IF(ISTRAN(8).GE.1)THEN
            IF(IWQRST.EQ.1) CALL WWQRST
            IF(IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST
          ENDIF
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  RECORD TIME
C
C **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
C **  SUPPORTED ON OTHER SYSTEMS.
C
      IF(NTIMER.EQ.NTSPTC)THEN
      CALL TIMELOG(N,TIMEDAY)
C     CALL DTIME (TARRAY)
C     WRITE(9,200)N, TARRAY(1),TARRAY(2)
C     CALL FLUSH(9)
      NTIMER=1
      ELSE
      NTIMER=NTIMER+1
      ENDIF
C
C**********************************************************************C
C
      IF(ISHOW.EQ.1) CALL SHOWVAL1
      IF(ISHOW.EQ.2) CALL SHOWVAL2
      IF(ISHOW.EQ.3) CALL SHOWVAL3
C
C**********************************************************************C
C
      GOTO 1001
 1000 CONTINUE
C
CXX DYNSTEP 1000 CONTINUE
C
C**********************************************************************C
C
C **  TIME LOOP COMPLETED
C
      IF(ISCRAY.EQ.0)THEN
        THDMT=THDMT+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        THDMT=T2TMP-TTMP
        WTHDMT=(WT2TMP-WTTMP)*0.001
      ENDIF
C
C**********************************************************************C
C**********************************************************************C
C
C **  CALCULATE VECTOR POTENTIAL AND VECTOR POTENTIAL TRANSPORTS
C **  USING RESULTS OF THE HARMONIC ANALYSIS
C
C     IF(ISVPTHA.NE.1) GOTO 2000
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LS=LSC(L)
C     VPZ(L,K)=TCVP*SUB(L)*SUB(LS)*SVB(L)*SVB(L-1)*HMC(L)*
C    &         ((AMSU(L,K)+AMSU(LS,K))*(AMCV(L,K)+AMCV(L-1,K))
C    &         -(AMCU(L,K)+AMCU(LS,K))*(AMSV(L,K)+AMSV(L-1,K)))
C     ENDDO
C     ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)
C     VPX(L,K)=TCVP*((AMSV(L,K)+AMSV(L,K+1))*(AMCW(L,K)+AMCW(LS,K))
C    &              -(AMCV(L,K)+AMCV(L,K+1))*(AMSW(L,K)+AMSW(LS,K)))
C     VPY(L,K)=TCVP*((AMSW(L,K)+AMSW(L-1,K))*(AMCU(L,K)+AMCU(L,K+1))
C    &              -(AMCW(L,K)+AMCW(L-1,K))*(AMSU(L,K)+AMSU(L,K+1)))
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LS=LSC(L)
C     LN=LNC(L)
C     UVPT(L,K)=(VPZ(LN,K)-VPZ(L,K))/DYU(L)-DZI*(VPY(L,K)-VPY(L,K-1))
C     VVPT(L,K)=DZI*(VPX(L,K)-VPX(L,K-1))-(VPZ(L+1,K)-VPZ(L,K))/DXV(L)
C     ENDDO
C     ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)
C     LN=LNC(L)
C     WVPT(L,K)=(VPY(L+1,K)-VPY(L,K))/DXP(L)
C    &         -(VPX(LN,K)-VPX(L,K))/DYP(L)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
 2000 CONTINUE
C
C**********************************************************************C
C
C **  PRINT FINAL RESULTS
C
      CALL OUTPUT2
C
C**********************************************************************C
C
C **  WRITE RESTART FILE
C
      IF(ISRESTO.EQ.-1.OR.ISRESTO.EQ.-11)THEN
        CALL RESTOUT(0)
        IF(ISTRAN(8).GE.1)THEN
          IF(IWQRST.EQ.1) CALL WWQRST
          IF(IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST
        ENDIF
      ENDIF
      IF(ISRESTO.EQ.-2)THEN
        CALL RESTMOD
      ENDIF
C
C**********************************************************************C
C
C **  COMPLETE LEAST SQUARES HARMONIC ANALYSIS
C
      LSLSHA=1
      IF(ISLSHA.EQ.1) CALL LSQHARM
C
C**********************************************************************C
C
C **  OUTPUT COURANT NUMBER DIAGNOSTICS
C
      OPEN(1,FILE='CFLMAX.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='CFLMAX.OUT')
C
      DO L=2,LA
       WRITE(1,1991)IL(L),JL(L),(CFLUUU(L,K),K=1,KC)
       WRITE(1,1992)(CFLVVV(L,K),K=1,KC)
       WRITE(1,1992)(CFLWWW(L,K),K=1,KC)
       WRITE(1,1992)(CFLCAC(L,K),K=1,KC)
      ENDDO
C
      CLOSE(1)
C
 1991 FORMAT(2I5,52F8.3)
 1992 FORMAT(10X,52F8.3)
 1993 FORMAT(2I5,E13.5)
C
C**********************************************************************C
C
C **  OUTPUT COSMETIC VOLUME LOSSES FORM DRY CELLS
C
      IF(NDRYSTP.LT.0) THEN
C
      OPEN(1,FILE='DRYLOSS.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='DRYLOSS.OUT')
C
      DO L=2,LA
       WRITE(1,1993)IL(L),JL(L),VDWASTE(L)
      ENDDO
C
      CLOSE(1)
C
	ENDIF
C
C
C **  OUTPUT FINAL MASS AND VOLUME BALANCES
C
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T5
	  ENDIF
	ENDIF
C
C**********************************************************************C
C
      RETURN
      END
