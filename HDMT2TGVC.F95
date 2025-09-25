!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE HDMT2TGVC
!
! **  SUBROUTINE HDMT EXECUTES THE FULL HYDRODYNAMIC AND MASS TRANSPORT
! **  TIME INTERGATION
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
! 05/60/2009        John Hamrick        05/60/2009        John Hamrick
!
!----------------------------------------------------------------------C
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      DIMENSION ISSBCP(LCM),WCOREW(LCM),WCORNS(LCM),LCORNER(LCM),LCORNWE(LCM),LCORNSN(LCM)
!
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        TTMP=SECOND( )
        CALL TIMEF(WTTMP)
      ENDIF
!
      FOURDPI=4./PI
      ISTL=2
      IS2TL=1
      ICALLTP=0
!
!**********************************************************************C
!
! **  SET FLAGS FOR CORNER CELL BED STRESS CORRECTIONS
!
      IF(ISCORTBC.GE.1) THEN
!
! **  SET FLAG FOR CELLS HAVING VOLUME SOURCE OR SINKS
!
      DO L=1,LC
	  ISSBCP(L)=0
	ENDDO
!
	DO L=2,LA
	  IF(RSSBCE(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCW(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCN(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCS(L).GT.1.5)ISSBCP(L)=1
	ENDDO
!
      ENDIF
!
!**********************************************************************C
!
! **  REINITIALIZE VARIABLES
!
      DO L=2,LA
       H1P(L)=HP(L)
       H1U(L)=HU(L)
       H1UI(L)=HUI(L)
       H1V(L)=HV(L)
       H1VI(L)=HVI(L)
       UHDY1E(L)=UHDYE(L)
       VHDX1E(L)=VHDXE(L)
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
       U1(L,K)=U(L,K)
       V1(L,K)=V(L,K)
       UHDY1(L,K)=UHDY(L,K)
       VHDX1(L,K)=VHDX(L,K)
      ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  INITIALIZE COURNT NUMBER DIAGNOSTICS
!
      DO K=1,KC
      DO L=2,LA
       CFLUUU(L,K)=0.
       CFLVVV(L,K)=0.
       CFLWWW(L,K)=0.
       CFLCAC(L,K)=0.
      ENDDO
      ENDDO
!
!**********************************************************************C
!
      ILOGC=0
!
!**********************************************************************C
!
! **  CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING
! **  CALCULATE VELOCITY GRADIENTS
!
!----------------------------------------------------------------------C
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      DO L=2,LA
	KBU=KGVCU(L)
	KBV=KGVCV(L)
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
      HMC(L)=0.25*(HMP(L)+HMP(L-1)+HMP(LS)+HMP(LSW))
      UV(L)=0.25*(HP(LS)*(U(LSE,KBV)+U(LS,KBV))+HP(L)*(U(L+1,KBV)+U(L,KBV)))*HVI(L)
      U1V(L)=0.25*(H1P(LS)*(U1(LSE,KBV)+U1(LS,KBV))+H1P(L)*(U1(L+1,KBV)+U1(L,KBV)))*H1VI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,KBU)+V(L-1,KBU))+HP(L)*(V(LN,KBU)+V(L,KBU)))*HUI(L)
      V1U(L)=0.25*(H1P(L-1)*(V1(LNW,KBU)+V1(L-1,KBU))+H1P(L)*(V1(LN,KBU)+V1(L,KBU)))*H1UI(L)
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
!
      IF(ISWAVE.EQ.1) CALL WAVEBL
      IF(ISWAVE.EQ.2) CALL WAVESXY
!
!**********************************************************************C
!
! **  FIRST CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
!
      CALL CALTBXY(ISTL,IS2TL)
!
!**********************************************************************C
!
! **  CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES
!
      IF(ISHDMF.GE.1) CALL CALHDMF
!
!**********************************************************************C
!
! **  CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N
!
!----------------------------------------------------------------------C
!
      N=-1
      CALL CALTSXY
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      DO L=2,LA
	  KBU=KGVCU(L)
	  KBV=KGVCV(L)
        TBX1(L)=(AVCON1*H1UI(L)+STBX(L)*SQRT(V1U(L)*V1U(L)+U1(L,KBU)*U1(L,KBU)))*U1(L,KBU)
        TBY1(L)=(AVCON1*H1VI(L)+STBY(L)*SQRT(U1V(L)*U1V(L)+V1(L,KBV)*V1(L,KBV)))*V1(L,KBV)
        TSX1(L)=TSX(L)
        TSY1(L)=TSY(L)
      ENDDO
!
!**********************************************************************C
!
! **  SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      CALL CALTBXY(ISTL,IS2TL)
!
!**********************************************************************C
!
! **  SET BOTTOM AND SURFACE STRESSES
!
!----------------------------------------------------------------------C
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      DO L=2,LA
	  KBU=KGVCU(L)
	  KBV=KGVCV(L)
        TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)+U(L,KBU)*U(L,KBU)))*U(L,KBU)
        TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)+V(L,KBV)*V(L,KBV)))*V(L,KBV)
      ENDDO
!
      N=0
      CALL CALTSXY
!
!----------------------------------------------------------------------C
!
! **  SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
!
      DO L=2,LA
        HDFUFX(L)=1.
        HDFUFY(L)=1.
        HDFUF(L)=1.
      ENDDO
!
      IF(ISBSDFUF.GE.1)THEN
!jh060305      HDFUFM=1.E-12
      HDFUFM=0.0

      DO L=2,LA
        LS=LSC(L)
        HDFUFX(L)=HDFUFM+G*SUB(L)*HU(L)*(BELV(L-1)-BELV(L))*DXIU(L)
        HDFUFY(L)=HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )-BELV(L))*DYIV(L)
      ENDDO
!
!jh060305      DO L=2,LA
!jh060305        HDFUFX(L)=TBX(L)/HDFUFX(L)
!jh060305        HDFUFY(L)=TBY(L)/HDFUFY(L)
!jh060305      ENDDO
!
!jh060305      DO L=2,LA
!jh060305        HDFUFX(L)=MAX(HDFUFX(L),-1.0)
!jh060305        HDFUFY(L)=MAX(HDFUFY(L),-1.0)
!jh060305      ENDDO
!
!jh060305      DO L=2,LA
!jh060305        HDFUFX(L)=MIN(HDFUFX(L),1.0)
!jh060305        HDFUFY(L)=MIN(HDFUFY(L),1.0)
!jh060305      ENDDO
!
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
!
      ENDIF
!
!**********************************************************************C
!
! **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
!
!----------------------------------------------------------------------C
!
!     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
!
      IF(ISWAVE.EQ.0)THEN
!
!----------------------------------------------------------------------c
!
       IF(ISCORTBC.EQ.0) THEN
!
!
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO
!

      DO L=2,LA
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
      ENDDO
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      DO L=2,LA
	 KBPM=KGVCP(L)-1
       QQ1(L,KBPM)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX1(L))**2+(RSSBCN(L)          &
                   *TVAR3N(L)+RSSBCS(L)*TBY1(L))**2)
       QQ1(L,KC)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX1(L))**2+(RSSBCN(L)            &
                 *TVAR3S(L)+RSSBCS(L)*TSY1(L))**2)
      ENDDO
!
      DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      DO L=2,LA
	 KBPM=KGVCP(L)-1
       QQ(L,KBPM)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2+(RSSBCN(L)*TVAR3N(L)  &
                  +RSSBCS(L)*TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2+(RSSBCN(L)*TVAR3S(L)    &
                +RSSBCS(L)*TSY(L))**2)
      ENDDO
!
       DO L=2,LA
	   KBPM=KGVCP(L)-1
	   TAUBSED(L)=QQ(L,KBPM)/CTURB2
	   TAUBSND(L)=QQ(L,KBPM)/CTURB2
       ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------c
!
       IF(ISCORTBC.GE.1) THEN

       IF(ISCORTBCD.GE.1)THEN
         OPEN(1,FILE='ADJSTRESSE.OUT')
         CLOSE(1,STATUS='DELETE')
 	 ENDIF
!
       OPEN(1,FILE='TBCORINIT.OUT')
!
       DO L=2,LA
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       ENDDO
!
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO
!
       DO L=2,LA
         IF(ISSBCP(L).EQ.0)THEN
	     IF(SUB(L+1).LT.0.5) WCOREST(L)=FSCORTBCV(L)
	     IF(SUB(L).LT.0.5) WCORWST(L)=FSCORTBCV(L)
	     IF(SVB(LNC(L)).LT.0.5) WCORNTH(L)=FSCORTBCV(L)
	     IF(SVB(L).LT.0.5) WCORSTH(L)=FSCORTBCV(L)
	   ENDIF
	 ENDDO
!
       DO L=2,LA
	   WCOREW(L)=1./(WCOREST(L)+WCORWST(L))
	   WCORNS(L)=1./(WCORNTH(L)+WCORSTH(L))
	 ENDDO
!
       DO L=2,LA
         WCOREST(L)=WCOREST(L)*WCOREW(L)
	   WCORWST(L)=WCORWST(L)*WCOREW(L)
	   WCORNTH(L)=WCORNTH(L)*WCORNS(L)
	   WCORSTH(L)=WCORSTH(L)*WCORNS(L)
       ENDDO
!
       DO L=2,LA
!      QQ1(L,0 )=0.5*CTURB2*SQRT(
!     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX1(L))**2
!     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY1(L))**2)
	 KBPM=KGVCP(L)-1
       QQ1(L,KBPM)=CTURB2*SQRT((RSSBCE(L)*WCOREST(L)*TVAR3E(L)+RSSBCW(L)*WCORWST(L)*TBX1(L))**2   &
                   +(RSSBCN(L)*WCORNTH(L)*TVAR3N(L)+RSSBCS(L)*WCORSTH(L)*TBY1(L))**2)
       QQ1(L,KC)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX1(L))**2+(RSSBCN(L)*TVAR3S(L)  &
                 +RSSBCS(L)*TSY1(L))**2)
       ENDDO

       DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
!
       DO L=2,LA
!      QQ(L,0 )=0.5*CTURB2*SQRT(
!     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2
!     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2)
	   KBPM=KGVCP(L)-1
         QQ(L,KBPM)=CTURB2*SQRT((RSSBCE(L)*WCOREST(L)*TVAR3E(L)+RSSBCW(L)*WCORWST(L)*TBX(L))**2   &
                    +(RSSBCN(L)*WCORNTH(L)*TVAR3N(L)+RSSBCS(L)*WCORSTH(L)*TBY(L))**2)
         QQ(L,KC)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2+(RSSBCN(L)*TVAR3S(L)+ &
                  RSSBCS(L)*TSY(L))**2)
       ENDDO
!
       DO L=2,LA
	   KBPM=KGVCP(L)-1
	   TAUBSED(L)=QQ(L,KBPM)/CTURB2
	   TAUBSND(L)=QQ(L,KBPM)/CTURB2
       ENDDO
!
       DO L=2,LA
         IF(WCORSTH(L).LT.0.49.OR.WCORSTH(L).GT.0.51)THEN
           IF(WCORWST(L).LT.0.49.OR.WCORWST(L).GT.0.51)THEN
             WRITE(1,3678)IL(L),JL(L),WCORWST(L),WCOREST(L),WCORSTH(L),WCORNTH(L)
	     ENDIF
	   ENDIF
       ENDDO

       CLOSE(1)

       ENDIF
!C
!C----------------------------------------------------------------------c
!C
      ENDIF
!C
!C**********************************************************************C
!C
!C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
!C
!C----------------------------------------------------------------------C
!C
!C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
!C
      IF(ISWAVE.GE.1)THEN
!C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       ENDDO
      ENDDO
!C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
	   KBPM=KGVCP(L)-1
         TAUBC2=0.25*( (TVAR3E(L)+TBX1(L))**2+(TVAR3N(L)+TBY1(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U1(L+1,KBPM)+U1(L,KBPM))+1.E-12
         VTMP=0.5*STCUV(L)*(V1(LN,KBPM)+V1(L,KBPM))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))+FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ1(L,0 )=CTURB2*SQRT(TAUB2)
         QQ1(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX1(L))**2+(TVAR3S(L)+TSY1(L))**2)
       ENDDO
      ENDDO
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
      ENDDO
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
	   KBPM=KGVCP(L)-1
         TAUBC2=0.25*( (TVAR3E(L)+TBX(L))**2+(TVAR3N(L)+TBY(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U(L+1,KBPM)+U(L,KBPM))+1.E-12
         VTMP=0.5*STCUV(L)*(V(LN,KBPM)+V(L,KBPM))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))+FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ(L,0 )=CTURB2*SQRT(TAUB2)
         QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2+(TVAR3S(L)+TSY(L))**2)
       ENDDO
      ENDDO
!
      ENDIF
!
!     ENDIF
!
!**********************************************************************C
!
! **   SET SWITCHES FOR TWO TIME LEVEL INTEGRATION
!
       ISTL=2
       IS2TL=1
       DELT=DT
       DELTD2=DT/2.
       DZDDELT=DZ/DELT
       ROLD=0.
       RNEW=1.
!
!**********************************************************************C
!**********************************************************************C
!
! **  BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT
! **  CALCULATION
!
! **  SET CYCLE COUNTER AND CALL TIMER
!
      NTIMER=0
      N=0
      TIMEDAY=TCON*TBEGIN/86400.
!
! **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
! **  SUPPORTED ON OTHER SYSTEMS.
!
      CALL TIMELOG(N,TIMEDAY)
!     CALL DTIME (TARRAY)
!     WRITE(9,200)N, TARRAY(1),TARRAY(2)
!     CALL FLUSH(9)
  200 FORMAT('  N=',I10,5X,'USER TIME=',E12.4,5X,'SYSTEM TIME=',E12.4)
      NTIMER=1
!
      NINCRMT=1
!
!----------------------------------------------------------------------C
!
!      DO 1000 N=1,NTS
!
 1001 CONTINUE
      IF(N.GE.NTS)GO TO 1000
!
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
!
      IF(ILOGC.EQ.NTSMMT)THEN
!        CLOSE(8,STATUS='DELETE')                              !hnr 7/27/2009
!        OPEN(8,FILE='EFDCLOG.OUT',STATUS='UNKNOWN')           !hnr 7/27/2009
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
!
      IF(ISDYNSTP.EQ.0)THEN
        ILOGC=ILOGC+1
      ELSE
        ILOGC=ILOGC+NINCRMT
      ENDIF
!
      IF(N.LE.NLTS) SNLT=0.
      IF(N.GT.NLTS.AND.N.LE.NTTS)THEN
         NTMP1=N-NLTS
         NTMP2=NTTS-NLTS+1
         SNLT=FLOAT(NTMP1)/FLOAT(NTMP2)
      ENDIF
      IF(N.GT.NTTS) SNLT=1.
!
      IF(N.LE.NTSVB)THEN
       GP=GPO*(FLOAT(N)/FLOAT(NTSVB))
      ELSE
       GP=GPO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  INITIALIZE VOLUME, MASS, MOMENTUM, AND ENERGY BALANCE
!
!      IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
!         CALL CALBAL1
!         NTMP=MOD(N,2)
!         IF(NTMP.EQ.0)THEN
!           CALL CBALEV1
!          ELSE
!           CALL CBALOD1
!         ENDIF
!       ENDIF
!
!  ** INITIALIZE SEDIMENT BUDGET CALCULATION   (DLK 10/15)
!
!      IF(NCTBC.NE.NTSTBC.AND.ISSBAL.GE.1)THEN
!         CALL BUDGET1
!         NTMP=MOD(N,2)
!         IF(NTMP.EQ.0)THEN
!           CALL BUDGEV1
!          ELSE
!           CALL BUDGOD1
!         ENDIF
!       ENDIF
!
! **  INITIALIZE TWO-TIME LEVEL BALANCES
!
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T1
	  ENDIF
	ENDIF
!
!----------------------------------------------------------------------C
!
! **  REENTER HERE FOR TWO TIME LEVEL CORRECTION
!
  500 CONTINUE
!
!**********************************************************************C
!
! **  CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)
!
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
      IF(KC.GT.1)THEN
        IF(ISQQ.EQ.1)THEN
	    IF(ISTOPT(0).EQ.0)CALL CALAVBOLD (ISTL)
!
! MODIFIED FOR GVC - JMH 08/03/04
!
	    IF(ISTOPT(0).GE.1)CALL CALAVBGVC (ISTL)
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
!
!**********************************************************************C
!
! **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
!
        IF(ISWAVE.EQ.1) CALL WAVEBL
        IF(ISWAVE.EQ.2) CALL WAVESXY
!
!**********************************************************************C
!
! **  CALCULATE EXPLICIT MOMENTUM EQUATION TERMS
!
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
!
      IF(IS2TIM.EQ.1) CALL CALEXP2TGVC
!      IF(IS2TIM.EQ.2) CALL CALIMP2T
!
      IF(ISCRAY.EQ.0)THEN
        TCEXP=TCEXP+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCEXP=TCEXP+T2TMP-T1TMP
        WTCEXP=WTCEXP+(WT2TMP-WT1TMP)*0.001
      ENDIF
!
!**********************************************************************C
!
! **  UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS, CONCENTRATIONS,
! **  VEGETATION CHARACTERISTICS AND SURFACE ELEVATIONS
!
      CALL CALCSER (ISTL)
      CALL CALVEGSER (ISTL)
      CALL CALQVS (ISTL)
      PSERT(0)=0.
      IF(NPSER.GE.1) CALL CALPSER (ISTL)
!
!**********************************************************************C
!
! **  WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
!     SETUP CALL
!
!----------------------------------------------------------------------C
!
      IF(ISWSEDA.GT.0.OR.ISUVDA.GT.0) CALL PUVDASM(ISTL,1)
!
!**********************************************************************C
!
! **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
! **  MODE FORCING (S&M SOLUTION ONLY)
!
!----------------------------------------------------------------------C
!
      IF(ISCDMA.GE.3.AND.ISCDMA.LE.8)THEN
!      IF(ISTL.EQ.3)THEN
!
      DO L=2,LA
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
      ENDDO
!
      CALL CALTSXY
!
      DO L=2,LA
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
!
      DO L=2,LA
        FXE(L)=FXE(L)+DT*SUB(L)*DYU(L)*(TSX(L)-TSX1(L))
        FYE(L)=FYE(L)+DT*SVB(L)*DXV(L)*(TSY(L)-TSY1(L))
      ENDDO
!
      ENDIF
!      ENDIF
!
!**********************************************************************C
!
! **  SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE
!
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
!
      IF(ISCHAN.EQ.0.AND.ISDRY.EQ.0) CALL CALPUV2TGVC
      IF(ISCHAN.GE.1.OR.ISDRY.GE.1) CALL CALPUV2C
!
      IF(ISCRAY.EQ.0)THEN
        TPUV=TPUV+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TPUV=TPUV+T2TMP-T1TMP
        WTPUV=WTPUV+(WT2TMP-WT1TMP)*0.001
      ENDIF
!
!**********************************************************************C
!
! **  WRITE DIAGNOSTICS
!
!----------------------------------------------------------------------C
!
! **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
! **  SUPPORTED ON OTHER SYSTEMS.
!
      IF(ISLOG.GE.1)THEN
!      WRITE(8,17)N,ITER,RSQ,CFMAX,AVMAX,ABMIN,ABMAX,ABMIN     !hnr 7/27/2009
!     CALL FLUSH(8)
      ENDIF
!
!  17 FORMAT('  N,ITER,AVMA,AVMI,ABMA,ABMI',2I5,4(1X,F8.4))
   17 FORMAT('  N,ITER,RSQ,CFMAX,AVMA,AVMI,ABMA,ABMI',I7,I5,2E12.4,4(1X,E13.4))
!
      ERRMAX=MAX(ERRMAX,ERR)
      ERRMIN=MIN(ERRMIN,ERR)
      ITRMAX=MAX(ITRMAX,ITER)
      IRRMIN=MIN(ITRMIN,ITER)
!
!**********************************************************************C
!
! **  ADVANCE INTERNAL VARIABLES FOR THREE TIME LEVEL STEP
!
!----------------------------------------------------------------------C
!
!      IF(ISTL.EQ.3)THEN
!
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
!
!      ENDIF
!
!**********************************************************************C
!
! **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
! **  MODE FORCING
!
!----------------------------------------------------------------------C
!
      IF(ISCDMA.LE.2.OR.ISCDMA.GE.9)THEN
!     IF(ISTL.EQ.3)THEN
!
      DO L=2,LA
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
      ENDDO
!
      CALL CALTSXY
!
      DO L=2,LA
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
      ENDDO
!
!     ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W
!
!----------------------------------------------------------------------C
!
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
!
      IF(KC.GT.1)THEN
        CALL CALUVWGVC (ISTL,IS2TL)
      ELSE
        DO L=2,LA
          UHDY(L,1)=UHDYE(L)
          U(L,1)=UHDYE(L)*HUI(L)*DYIU(L)
          VHDX(L,1)=VHDXE(L)
          V(L,1)=VHDXE(L)*HVI(L)*DXIV(L)
          W(L,1)=0.
        ENDDO
        CALL CALUVWGVC (ISTL,IS2TL)
      ENDIF
!
      IF(ISCRAY.EQ.0)THEN
        TUVW=TUVW+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TUVW=TUVW+T2TMP-T1TMP
        WTUVW=WTUVW+(WT2TMP-WT1TMP)*0.001
      ENDIF
!
!**********************************************************************C
!
! **  WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
!     DIAGNOSTIC CALL
!
!----------------------------------------------------------------------C
!
      IF(ISWSEDA.GT.0.OR.ISUVDA.GT.0) CALL PUVDASM(ISTL,2)
!
!**********************************************************************C
!
! **  CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS
! **  AT TIME LEVEL (N+1)
!
!----------------------------------------------------------------------C
!
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      CALL CALCONCGVC (ISTL,IS2TL)
!
!----------------------------------------------------------------------C
!
      DO K=1,KB
      DO L=1,LC
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      ENDDO
      ENDDO
!
      DO NS=1,NSED
       DO K=1,KB
       DO L=1,LC
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
       ENDDO
       ENDDO
      ENDDO
!
      DO NS=1,NSND
       DO K=1,KB
       DO L=1,LC
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
       ENDDO
       ENDDO
      ENDDO
!
      DO K=1,KC
       DO L=1,LC
        SEDT(L,K)=0.
        SNDT(L,K)=0.
       ENDDO
      ENDDO
!
      DO NS=1,NSED
       DO K=1,KC
        DO L=1,LC
         SEDT(L,K)=SEDT(L,K)+SED(L,K,NS)
        ENDDO
       ENDDO
      ENDDO
!
      DO NS=1,NSND
       DO K=1,KC
        DO L=1,LC
         SNDT(L,K)=SNDT(L,K)+SND(L,K,NS)
        ENDDO
       ENDDO
      ENDDO
!
!JH5/13/97      DO NT=1,NTOX
!JH5/13/97       DO K=1,KC
!JH5/13/97        DO L=1,LC
!JH5/13/97         TOXPFT(L,K,NT)=0.
!JH5/13/97        ENDDO
!JH5/13/97       ENDDO
!JH5/13/97      ENDDO
!
!JH5/13/97      DO NT=1,NTOX
!JH5/13/97       DO NS=1,NSED+NSND
!JH5/13/97        DO K=1,KC
!JH5/13/97         DO L=1,LC
!JH5/13/97          TOXPFT(L,K,NT)=TOXPFT(L,K,NT)+TOXPF(L,K,NS,NT)
!JH5/13/97         ENDDO
!JH5/13/97        ENDDO
!JH5/13/97       ENDDO
!JH5/13/97      ENDDO
!
!----------------------------------------------------------------------C
!
! **  CHECK RANGE OF SALINITY AND DYE CONCENTRATION
!
      IF(ISMMC.EQ.1)THEN
!
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
!
      WRITE(6,6001)N
      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN
!
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
!
      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN
!
!      WRITE(8,6004)SALMAX,IMAX,JMAX,KMAX          !hnr 7/27/2009
!      WRITE(8,6005)SALMIN,IMIN,JMIN,KMIN          !hnr 7/27/2009
!
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
!
      WRITE(6,6006)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6007)SALMIN,IMIN,JMIN,KMIN
!
      ENDIF
!
!
      IF(ISMMC.EQ.2)THEN
!
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
!
      WRITE(6,6001)N
      WRITE(6,6008)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6009)SALMIN,IMIN,JMIN,KMIN
!
      ENDIF
!
 6001 FORMAT('  N=',I10)
 6002 FORMAT('  SALMAX=',F14.4,5X,'I,J,K=',(3I10))
 6003 FORMAT('  SALMIN=',F14.4,5X,'I,J,K=',(3I10))
 6004 FORMAT('  DYEMAX=',F14.4,5X,'I,J,K=',(3I10))
 6005 FORMAT('  DYEMIN=',F14.4,5X,'I,J,K=',(3I10))
 6006 FORMAT('  SFLMAX=',F14.4,5X,'I,J,K=',(3I10))
 6007 FORMAT('  SFLMIN=',F14.4,5X,'I,J,K=',(3I10))
 6008 FORMAT('  TEMMAX=',F14.4,5X,'I,J,K=',(3I10))
 6009 FORMAT('  TEMMIN=',F14.4,5X,'I,J,K=',(3I10))
!
!**********************************************************************C
!
! **  CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT
! **  CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOULBE TIME
! **  STEP TRANSPORT FIELD
!
!----------------------------------------------------------------------C
!
      IF(ISTRAN(4).GE.1.OR.ISTRAN(8).GE.1)THEN
!2T      NTMP=MOD(N,2)
!2T      IF(NTMP.EQ.0.AND.ISTL.EQ.3)THEN
!
! **  CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION
!
       DO L=2,LA
        HWQ(L)=HP(L)
        WWQ(L,0)=0.
       ENDDO
!
      DO K=1,KC
       DO L=2,LA
        UHDYWQ(L,K)=UHDY2(L,K)
        VHDXWQ(L,K)=VHDX2(L,K)
        UWQ(L,K)=U2(L,K)
        VWQ(L,K)=V2(L,K)
        WWQ(L,K)=W2(L,K)
       ENDDO
      ENDDO
!
!
      IF(ISDIVEX.EQ.3)THEN
!
      DO K=1,KC
       DO L=2,LA
        UHDY2E(L)=UHDY2E(L)+GVCSCLU(L)*UHDYWQ(L,K)*DZC(K)
        VHDX2E(L)=VHDX2E(L)+GVCSCLV(L)*VHDXWQ(L,K)*DZC(K)
       ENDDO
      ENDDO
!
      DO L=2,LA
       TVAR3E(L)=GVCSCLP(L)*UHDY2E(L+1)
       TVAR3N(L)=GVCSCLP(L)*VHDX2E(LNC(L))
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
        TVAR2E(L,K)=GVCSCLU(L+1   )*UHDYWQ(L+1   ,K)
        TVAR2N(L,K)=GVCSCLV(LNC(L))*VHDXWQ(LNC(L),K)
      ENDDO
      ENDDO
!
	  OPEN(1,FILE='DIV3GVC.OUT')
	  CLOSE(1,STATUS='DELETE')
	  OPEN(1,FILE='DIV3GVC.OUT')
        K=KC
        DO L=2,LA
          TVAR3S(L)=(WWQ(L,K-1)-DZC(K)*(TVAR2E(L,K)-GVCSCLU(L)*UHDYWQ(L,K)-TVAR3E(L)+GVCSCLP(L)*  &
                    UHDY2E(L)+TVAR2N(L,K)-GVCSCLV(L)*VHDXWQ(L,K)-TVAR3N(L)+GVCSCLP(L)*VHDX2E(L))  &
                    *DXYIP(L))+( QSUM(L,K)-DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
        ENDDO
	  DO L=2,LA
	    WRITE(1,1492)IL(L),JL(L),TVAR3S(L)
	  ENDDO
	  CLOSE(1)
!
	ENDIF
!
!     ADD CHANNEL INTERACTIONS
!

      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        IF(MDCHTYP(NMD).EQ.1)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))+DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
          HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))-DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
        ENDIF
        IF(MDCHTYP(NMD).EQ.2)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))+DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
          HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))-DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
        ENDIF
        IF(MDCHTYP(NMD).EQ.3)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))+DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))+DT2*             &
          DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
          HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))-DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
          HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))-DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
        ENDIF
        ENDDO
      ENDIF
!
!     END ADD CHANNEL INTERACTIONS
!
!      IF(ISTRAN(6).GE.1) CALL CALWQC(2)
      IF(ISTRAN(8).GE.1) CALL WQ3D
      IF(ISTRAN(4).GE.1) CALL CALSFT(2)
!
      DO L=2,LA
        H2WQ(L)=HWQ(L)
      ENDDO
!
      ENDIF
!
 1492 FORMAT(2I6,E15.7)
!
!**********************************************************************C
!
! **  UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING
! **  AN EQUATION OF STATE
!
!      IF(ISTL.EQ.3)THEN
       DO K=1,KC
       DO L=2,LA
       B1(L,K)=B(L,K)
       ENDDO
       ENDDO
!      ENDIF
!
      IF(BSC.GT.1.E-6)THEN
        CALL CALBUOY
	ELSE
        DO K=1,KC
          DO L=2,LA
            B(L,K)=0.
          ENDDO
        ENDDO
	ENDIF
!
!      IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
!         CALL CALBAL4
!         NTMP=MOD(N,2)
!         IF(NTMP.EQ.0)THEN
!           CALL CBALEV4
!          ELSE
!           CALL CBALOD4
!         ENDIF
!      ENDIF
!
! **  CALL TWO-TIME LEVEL BALANCES
!
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T4
	  ENDIF
	ENDIF
!
!**********************************************************************C
!
! **  CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)
!
!----------------------------------------------------------------------C
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      DO L=2,LA
	KBU=KGVCU(L)
	KBV=KGVCV(L)
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
!      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
       TMPVAL=1.+SUBO(L)+SVBO(L)+SUBO(L)*SVBO(L)
      H1C(L)=(HP(L)+SUBO(L)*HP(L-1)+SVBO(L)*HP(LS)+SUBO(L)*SVBO(L)*HP(LSW))/TMPVAL
      UV(L)=0.25*(HP(LS)*(U(LSE,KBV)+U(LS,KBV))+HP(L)*(U(L+1,KBV)+U(L,KBV)))*HVI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,KBU)+V(L-1,KBU))+HP(L)*(V(LN,KBU)+V(L,KBU)))*HUI(L)
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES
! **  AT TIME LEVEL (N)
!
      IF(ISTL.NE.2.AND.ISHDMF.GE.1) CALL CALHDMF
!
!**********************************************************************C
!
! **  UPDATE BOTTOM STRESSES AND SURFACE AND BOTTOM TURBULENT
! **  INTENSITIES
!
!----------------------------------------------------------------------C
!
!      IF(ISTL.EQ.2)THEN
!
!      DO K=1,KC
!       DO L=2,LA
!        QQ(L,K)=SQRT(QQ(L,K)*QQ1(L,K))
!       ENDDO
!      ENDDO
!
!      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(ISTL.EQ.3)THEN
      IF(ISCDMA.EQ.2)THEN
!
! MODIFIED FOR GVC - JMH 08/03/04
!

      DO L=2,LA
	  KBPM=KGVCP(L)-1
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,KBPM)=QQ(L,KBPM)+QQ1(L,KBPM)
        QQ2(L,KC)=QQ(L,KC)+QQ1(L,KC)
        QQ1(L,KBPM)=QQ(L,KBPM)
        QQ1(L,KC)=QQ(L,KC)
      ENDDO
!
      ELSE
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      DO L=2,LA
	  KBPM=KGVCP(L)-1
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,KBPM)=QQ1(L,KBPM)+QQ1(L,KBPM)
        QQ2(L,KC)=QQ1(L,KC)+QQ1(L,KC)
        QQ1(L,KBPM)=QQ(L,KBPM)
        QQ1(L,KC)=QQ(L,KC)
      ENDDO
!
      ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE BOTTOM STRESS AT LEVEL (N+1)
!
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      CALL CALTBXY(ISTL,IS2TL)
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      DO L=2,LA
	  KBU=KGVCU(L)
	  KBV=KGVCV(L)
        TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)+U(L,KBU)*U(L,KBU)))*U(L,KBU)
        TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)+V(L,KBV)*V(L,KBV)))*V(L,KBV)
      ENDDO
!
!**********************************************************************C
!
! **  SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
!
      IF(ISBSDFUF.GE.1)THEN
!jh060305      HDFUFM=1.E-12
      HDFUFM=0.0
!
      DO L=2,LA
        LS=LSC(L)
        HDFUFX(L)=HDFUFM+G*SUB(L)*HU(L)*(BELV(L-1)-BELV(L))*DXIU(L)
        HDFUFY(L)=HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )-BELV(L))*DYIV(L)
      ENDDO
!
!jh060305      DO L=2,LA
!jh060305        HDFUFX(L)=TBX(L)/HDFUFX(L)
!jh060305        HDFUFY(L)=TBY(L)/HDFUFY(L)
!jh060305      ENDDO
!
!jh060305      DO L=2,LA
!jh060305        HDFUFX(L)=MAX(HDFUFX(L),-1.0)
!jh060305        HDFUFY(L)=MAX(HDFUFY(L),-1.0)
!jh060305      ENDDO
!
!jh060305      DO L=2,LA
!jh060305        HDFUFX(L)=MIN(HDFUFX(L),1.0)
!jh060305        HDFUFY(L)=MIN(HDFUFY(L),1.0)
!jh060305      ENDDO
!
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
!
      ENDIF
!
!**********************************************************************C
!
! **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
!
!----------------------------------------------------------------------C
!
!     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
!
      IF(ISWAVE.EQ.0)THEN
!
!----------------------------------------------------------------------c
!
       IF(ISCORTBC.EQ.0) THEN
!
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO
!
      DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
      ENDDO
!
! MODIFIED FOR GVC - JMH 08/03/04
!

      DO L=2,LA
	 KBPM=KGVCP(L)-1
       QQ(L,KBPM)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2+(RSSBCN(L)*TVAR3N(L)  &
                  +RSSBCS(L)*TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2+(RSSBCN(L)*TVAR3S(L)    &
                +RSSBCS(L)*TSY(L))**2)
      ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------c
!
       IF(ISCORTBC.GE.1) THEN
!
       IF(ISCORTBCD.GE.1)THEN
         NTMPVAL=MOD(N,NTSPTC)
	   IF(NTMPVAL.EQ.0)THEN
           OPEN(1,FILE='ADJSTRESSE.OUT',ACCESS='APPEND')
 	   ENDIF
 	 ENDIF
!
       DO L=2,LA
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
!
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO
!
       DO L=2,LA
         IF(ISSBCP(L).EQ.0)THEN
	     IF(SUB(L+1).LT.0.5)WCOREST(L)=FSCORTBCV(L)
	     IF(SUB(L).LT.0.5)WCORWST(L)=FSCORTBCV(L)
	     IF(SVB(LNC(L)).LT.0.5)WCORNTH(L)=FSCORTBCV(L)
	     IF(SVB(L).LT.0.5)WCORSTH(L)=FSCORTBCV(L)
	   ENDIF
	 ENDDO
!
       DO L=2,LA
	   WCOREW(L)=1./(WCOREST(L)+WCORWST(L))
	   WCORNS(L)=1./(WCORNTH(L)+WCORSTH(L))
	 ENDDO
!
       DO L=2,LA
         WCOREST(L)=WCOREST(L)*WCOREW(L)
	   WCORWST(L)=WCORWST(L)*WCOREW(L)
	   WCORNTH(L)=WCORNTH(L)*WCORNS(L)
	   WCORSTH(L)=WCORSTH(L)*WCORNS(L)
       ENDDO
!
       DO L=2,LA
!
!       QQ(L,0 )=0.5*CTURB2*SQRT(
!     &           (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2
!     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2)
	 KBPM=KGVCP(L)-1
       QQ(L,KBPM)=CTURB2*SQRT((RSSBCE(L)*WCOREST(L)*TVAR3E(L)+RSSBCW(L)*WCORWST(L)*TBX(L))**2     &
                  +(RSSBCN(L)*WCORNTH(L)*TVAR3N(L)+RSSBCS(L)*WCORSTH(L)*TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2+(RSSBCN(L)*TVAR3S(L)    &
                +RSSBCS(L)*TSY(L))**2)
       ENDDO
!
      IF(ISCORTBCD.GE.1.AND.NTMPVAL.EQ.0)THEN
!
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
!
      NCORCELLS=0
	DO L=2,LA
	  NCORCELLS=NCORCELLS+LCORNER(L)
	ENDDO
!
      WRITE(1,3675)TIMEDAY,NCORCELLS
!
      DO L=2,LA
        IF(LMASKDRY(L))THEN
	  IF(LCORNER(L).EQ.1)THEN
	    LWE=LCORNWE(L)
	    LSN=LCORNSN(L)
          TAUTMP=QQ(L,0)/CTURB2
          TAUTMPWE=QQ(LWE,0)/CTURB2
          TAUTMPSN=QQ(LSN,0)/CTURB2
 	    WRITE(1,3677)IL(L),JL(L),TAUTMP,TAUBSND(L),TAUBSED(L)
 	    WRITE(1,3676)IL(LWE),JL(LWE),TAUTMPWE,TAUBSND(LWE),TAUBSED(LWE)
 	    WRITE(1,3676)IL(LSN),JL(LSN),TAUTMPSN,TAUBSND(LSN),TAUBSED(LSN)
!          WRITE(1,3678)IL(L),JL(L),SUB(L),SUB(L+1),SVB(L),
!     &                                   SVB(LNC(L))
!          WRITE(1,3679)RSSBCW(L),RSSBCE(L),RSSBCS(L),RSSBCN(L)
!          WRITE(1,3679)WTMPWST,WTMPEST,WTMPSTH,WTMPNTH
!          WRITE(1,3680)TBX(L),TVAR3E(L),TBY(L),TVAR3N(L),
!     &             TAUTMP,TAUBSND(L)
!	    DO NX=1,NSND
!	      BLTMPVAL=QSBDLDX(L+1,NX)-QSBDLDX(L,NX)
!     &            +QSBDLDY(LNC(L),NX)-QSBDLDY(L,NX)
!            WRITE(1,3681)QSBDLDX(L,NX),QSBDLDX(L+1,NX),QSBDLDY(L,NX),
!     &             QSBDLDY(LNC(L),NX),BLTMPVAL,TAUR(NX+1)
!	    ENDDO
	  ENDIF
	  ENDIF
      ENDDO
!
      ENDIF

      CLOSE(1)
!
      ENDIF
!
!----------------------------------------------------------------------c
!
      ENDIF
!
 3678 FORMAT(2I6,4F13.3)
 3679 FORMAT(12x,4F13.3)
 3680 FORMAT(12x,6F13.5)
 3681 FORMAT(12X,5E13.4,F13.5)
 3677 FORMAT('CORNER',2I5,5E14.5)
 3676 FORMAT(6X,2I5,5E14.5)
 3675 FORMAT(F11.3,I6,' TIME IN DAYS AND NUMBER OF CORNERS')
!
!**********************************************************************C
!
! **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
!
!----------------------------------------------------------------------C
!
!     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
!
      IF(ISWAVE.GE.1)THEN
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
      ENDDO
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         TAUBC2=0.25*( (TVAR3E(L)+TBX(L))**2+(TVAR3N(L)+TBY(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))+FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ(L,0 )=CTURB2*SQRT(TAUB2)
         QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2+(TVAR3S(L)+TSY(L))**2)
       ENDDO
      ENDDO
!
      ENDIF
!
!     ENDIF
!
!**********************************************************************C
!
! **  CALCULATE TURBULENT INTENSITY SQUARED
!
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
      IF(KC.GT.1)THEN
        IF(ISQQ.EQ.1)THEN
	    IF(ISTOPT(0).EQ.0)CALL CALQQ2TOLD (ISTL)
	    IF(ISTOPT(0).GE.1)CALL CALQQ2TGVC (ISTL)
        ENDIF
        IF(ISQQ.EQ.2) CALL CALQQ2 (ISTL)
      ENDIF
!
      IF(ISCRAY.EQ.0)THEN
        TQQQ=TQQQ+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TQQQ=TQQQ+T2TMP-T1TMP
        WTQQQ=WTQQQ+(WT2TMP-WT1TMP)*0.001
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE MEAN MASS TRANSPORT FIELD
!
      IF(ISSSMMT.NE.2)THEN
        IF(ISICM.GE.1)THEN
          NTMP=MOD(N,2)
          IF(ISTL.EQ.3.AND.NTMP.EQ.0) CALL CALMMT
        ENDIF
      ENDIF
!
!      IF(ISSSMMT.NE.2) CALL CALMMT
!
!**********************************************************************C
!
! **  HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
! **  IF NCTBC EQ NTSTBC APPLY TRAPEZOIDAL CORRECTION
!
!----------------------------------------------------------------------C
!
!2T      IF(NCTBC.EQ.NTSTBC)THEN
!2T       NCTBC=0
!2T       ISTL=2
!2T       DELT=DT
!2T       DELTD2=0.5*DT
!2T       DZDDELT=DZ/DELT
!2T       ROLD=0.5
!2T       RNEW=0.5
!2T       GOTO 500
!2T      ELSE
!2T       NCTBC=NCTBC+1
!2T       ISTL=3
!2T       DELT=DT2
!2T       DELTD2=DT
!2T       DZDDELT=DZ/DELT
!2T       ROLD=0.
!2T       RNEW=1.
!2T      ENDIF
!
! **  WRITE TO TIME SERIES FILES
!
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
!
!DYN      IF(ISTMSR.GE.1)THEN
!DYN        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
!DYN          IF(NCTMSR.EQ.NWTMSR)THEN
!DYN            CALL TMSR
!DYN            ICALLTP=1
!DYN            NCTMSR=1
!DYN           ELSE
!DYN            NCTMSR=NCTMSR+1
!DYN          ENDIF
!DYN        ENDIF
!DYN      ENDIF
!
!
      IF(ISTMSR.GE.1)THEN
!        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
          IF(NCTMSR.GE.NWTMSR)THEN
            CALL TMSR
            NDIFF=NWTMSR-NCTMSR
            ICALLTP=1
            NCTMSR=NINCRMT+NDIFF
           ELSE
            NCTMSR=NCTMSR+NINCRMT
          ENDIF
!        ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  WRITE TO DUMP FILES
!
      IF(ISDUMP.GE.1)THEN                                                      !HNR_GHD 7/2019 BMD2                     
        IF(TIME.GE.TSDUMP.AND.TIME.LE.TEDUMP)THEN            !hnr
          IF(NCDUMP.GE.NSDUMP)THEN                           !hnr
            CALL DUMP                                        !hnr
!            NDIFF=NSDUMP-NCDUMP                              !hnr
            ICALLTP=1                                        !hnr
            NCDUMP=1                                         !hnr
!            NCDUMP=NINCRMT+NDIFF                             !hnr
           ELSE                                              !hnr
            NCDUMP=NCDUMP+1                                  !hnr
!            NCDUMP=NCDUMP+NINCRMT                            !hnr
          ENDIF                                              !hnr
        ENDIF                                                !hnr
      ENDIF                                                  !hnr
!
!**********************************************************************C
!
! **  WRITE TO bmd2 OUTPUT FILE                                                    !HNR_GHD 7/2019 BMD2     !HNR_GHD 6/2022 remove bmd file keep only bmd2                
!                                                                 
      IF (bmdflag.GE.1) THEN                                          
        IF (TIME.GE.TIBMD.AND.TIME.LE.TEBMD) THEN                  
          IF (NCBMD.EQ.NSBMD) THEN                                 
               CALL BMD2_DUMP
            NCBMD=1                                                 
          ELSE                                                      
            NCBMD=NCBMD+1                                          
          END IF                                                     
        END IF                                                       
      END IF                                                         
!
!**********************************************************************C
!
!**********************************************************************C
!
! **  WRITE TOXIC FLUX ACCUMULATION FILES
!
      IF(ISTOXALL.EQ.1)THEN
        MSTOXALL=MSTOXALL+1
        IF(N.GE.MSTOXALL)THEN
          CALL TOXALL
          MSTOXALL=MSTOXALL+(NTSPTC/NSTOXALL)
        ENDIF
      ENDIF
!
!**********************************************************************C
!
! ** WRITE BOTTOM VELOCITIES TO FILE FOR AJ MINE SEDIMENT BED MODEL
!     (DLK - ADDED 3/25/97)
!
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
!
!**********************************************************************C
!
! **  OUTPUT ZERO DIMENSION VOLUME BALANCE
!
!----------------------------------------------------------------------C
!
!      IF(ISDRY.GE.1)THEN
!        OPEN(1,FILE='ZVOLBAL.OUT',POSITION='APPEND',STATUS='UNKNOWN')
!        DO LS=1,LORMAX
!        IF(VOLZERD.GE.VOLSEL(LS).AND.VOLZERD.LT.VOLSEL(LS+1))THEN
!           WTM=VOLSEL(LS+1)-VOLZERD
!           WTMP=VOLZERD-VOLSEL(LS)
!           DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
!           WTM=WTM/DELVOL
!           WTMP=WTMP/DELVOL
!           SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
!           ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
!        ENDIF
!        ENDDO
!        IF(ISDYNSTP.EQ.0)THEN
!          TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
!        ELSE
!          TIME=TIMESEC/TCTMSR
!        ENDIF
!        WRITE(1,5304)TIME,SELZERD,ASFZERD,VOLZERD,VETZERD
!        CLOSE(1)
!      ENDIF
      ICALLTP=0
!
 5304 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
!
!**********************************************************************C
!
! **  WRITE VERTICAL SCALAR FIELD PROFILES
!
      IF(ISVSFP.EQ.1)THEN
        IF(N.GE.NBVSFP.AND.N.LE.NSVSFP)THEN
          CALL VSFP
        ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE MEAN MASS TRANSPORT FIELD
!
      IF(ISSSMMT.NE.2)THEN
        IF(ISICM.EQ.0) CALL CALMMT
      ENDIF
!
!      IF(ISSSMMT.NE.2) CALL CALMMT
!
!**********************************************************************C
!
! **  ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES
!
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
!
!**********************************************************************C
!
! **  CALCULATE VOLUME MASS, MOMENTUM AND ENERGY BALANCES
!
!      IF(ISBAL.GE.1)THEN
!         CALL CALBAL5
!         NTMP=MOD(N,2)
!         IF(NTMP.EQ.0)THEN
!           CALL CBALEV5
!          ELSE
!           CALL CBALOD5
!         ENDIF
!       ENDIF
!
!   SEDIMENT BUDGET CALCULATION     (DLK 10/15)
!
!       IF(ISSBAL.GE.1)THEN
!       CALL BUDGET5
!       ENDIF
!       NTMP=MOD(N,2)
!       IF(NTMP.EQ.0)THEN
!         CALL BUDGEV5
!        ELSE
!         CALL BUDGOD5
!       ENDIF
!
! **  CALL TWO-TIME LEVEL BALANCES
!
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T5
	  ENDIF
	ENDIF
!
!**********************************************************************C
!
! **  PERFORM AN M2 TIDE HARMONIC ANALYSIS EVERY 2 M2 PERIODS
!
!HARM      IF(ISHTA.EQ.1) CALL CALHTA
!
!**********************************************************************C
!
! **  CALCULATE DISPERSION COEFFICIENTS
!
!     IF(N.GE.NDISP)THEN
      IF(N.GE.NDISP.AND.NCTBC.EQ.1)THEN
       IF(ISDISP.EQ.2) CALL CALDISP2
       IF(ISDISP.EQ.3) CALL CALDISP3
      ENDIF
!
!**********************************************************************C
!
! **  PERFORM LEAST SQUARES HARMONIC ANALYSIS AT SELECTED LOCATIONS
!
      IF(ISLSHA.EQ.1.AND.N.EQ.NCLSHA)THEN
       CALL LSQHARM
       NCLSHA=NCLSHA+(NTSPTC/24)
      ENDIF
!
!**********************************************************************C
!
! **  PRINT INTERMEDIATE RESULTS
!
!----------------------------------------------------------------------C
!
      IF(NPRINT .EQ. NTSPP)THEN
       NPRINT=1
       CALL OUTPUT1
      ELSE
       NPRINT=NPRINT+1
      ENDIF
!
!**********************************************************************C
!
! **  WRITE TO TIME VARYING GRAPHICS FILES
!
!----------------------------------------------------------------------C
!
!DYN      IF(N.EQ.NCPPH.AND.ISPPH.EQ.1)THEN
!
      IF(N.GE.NCPPH.AND.ISPPH.GE.1)THEN
       CALL SURFPLT
       NCPPH=NCPPH+(NTSPTC/NPPPH)
      ENDIF
!
!
!----------------------------------------------------------------------C
!
!DYN      IF(N.EQ.NCBPH.AND.ISBPH.EQ.1)THEN
!
      IF(N.GE.NCBPH.AND.ISBPH.GE.1)THEN
!	IF(ISBEXP.EQ.0)THEN
       CALL BEDPLTH
       NCBPH=NCBPH+(NTSPTC/NPBPH)
!      ENDIF
      ENDIF
!
!----------------------------------------------------------------------C
!
!DYN      IF(N.EQ.NCVPH.AND.ISVPH.GE.1)THEN
!
      IPLTTMP=0
      IF(ISVPH.EQ.1.OR.ISVPH.EQ.2)IPLTTMP=1
      IF(N.GE.NCVPH.AND.IPLTTMP.EQ.1)THEN
       CALL VELPLTH
       NCVPH=NCVPH+(NTSPTC/NPVPH)
      ENDIF
!
!----------------------------------------------------------------------C
!
!DYN      IF(N.EQ.NCVPV.AND.ISVPV.GE.1)THEN
!
      IF(N.GE.NCVPV.AND.ISVPV.GE.1)THEN
       CALL VELPLTV
       NCVPV=NCVPV+(NTSPTC/NPVPV)
      ENDIF
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=1,LC
        TVAR1S(L,K)=TOX(L,K,1)
       ENDDO
      ENDDO
!
      IPLTTMP=0
      IF(ISSPH(1).EQ.1.OR.ISSPH(1).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(1).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(1).GE.1) CALL SALPLTH (1,SAL)
       NCSPH(1)=NCSPH(1)+(NTSPTC/NPSPH(1))
      ENDIF
!
      IPLTTMP=0
      IF(ISSPH(2).EQ.1.OR.ISSPH(2).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(2).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(2).GE.1) CALL SALPLTH (2,TEM)
       NCSPH(2)=NCSPH(2)+(NTSPTC/NPSPH(2))
      ENDIF
!
      IPLTTMP=0
      IF(ISSPH(3).EQ.1.OR.ISSPH(3).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(3).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(3).GE.1) CALL SALPLTH (3,DYE)
       NCSPH(3)=NCSPH(3)+(NTSPTC/NPSPH(3))
      ENDIF
!
      IPLTTMP=0
      IF(ISSPH(4).EQ.1.OR.ISSPH(4).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(4).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(4).GE.1) CALL SALPLTH (4,SFL)
       NCSPH(4)=NCSPH(4)+(NTSPTC/NPSPH(4))
      ENDIF
!
      IPLTTMP=0
      IF(ISSPH(5).EQ.1.OR.ISSPH(5).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(5).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(5).GE.1) CALL SALPLTH (5,TVAR1S)
       NCSPH(5)=NCSPH(5)+(NTSPTC/NPSPH(5))
      ENDIF
!
      IPLTTMP=0
      IF(ISSPH(6).EQ.1.OR.ISSPH(6).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(6).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(6).GE.1) CALL SALPLTH (6,SEDT)
       NCSPH(6)=NCSPH(6)+(NTSPTC/NPSPH(6))
      ENDIF
!
      IPLTTMP=0
      IF(ISSPH(7).EQ.1.OR.ISSPH(7).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(7).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(7).GE.1) CALL SALPLTH (7,SNDT)
       NCSPH(7)=NCSPH(7)+(NTSPTC/NPSPH(7))
      ENDIF
!
!----------------------------------------------------------------------C
!
      DO ITMP=1,7
      IF(N.GE.NCSPV(ITMP).AND.ISSPV(ITMP).GE.1)THEN
       CALL SALPLTV(ITMP)
       NCSPV(ITMP)=NCSPV(ITMP)+(NTSPTC/NPSPV(ITMP))
      ENDIF
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  WRITE EFDC EXPLORER FORMAT OUTPUT
!
      IF(ISSPH(8).EQ.1.OR.ISBEXP.EQ.1)THEN
      IF(N.GE.NCSPH(8))THEN
       CALL EEXPOUT(0)
       NCSPH(8)=NCSPH(8)+(NTSPTC/NPSPH(8))
      ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  WRITE TO TIME VARYING 3D HDF GRAPHICS FILES
!
!----------------------------------------------------------------------C
!
      IF(N.EQ.NC3DO.AND.IS3DO.EQ.1)THEN
       CALL OUT3D
       NC3DO=NC3DO+(NTSPTC/NP3DO)
      ENDIF
!
!**********************************************************************C
!
! **  WRITE RESTART FILE EVERY ISRESTO M2 TIDAL CYCLES
!
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
!
!**********************************************************************C
!
! **  RECORD TIME
!
! **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
! **  SUPPORTED ON OTHER SYSTEMS.
!
      IF(NTIMER.EQ.NTSPTC)THEN
      CALL TIMELOG(N,TIMEDAY)
!     CALL DTIME (TARRAY)
!     WRITE(9,200)N, TARRAY(1),TARRAY(2)
!     CALL FLUSH(9)
      NTIMER=1
      ELSE
      NTIMER=NTIMER+1
      ENDIF
!
!**********************************************************************C
!
! MODIFIED FOR GVC - JMH 08/03/04
!
      IF(ISHOW.EQ.1) CALL SHOWVAL1
      IF(ISHOW.EQ.2) CALL SHOWVAL2
      IF(ISHOW.EQ.3) CALL SHOWVAL3
!
!**********************************************************************C
!
      GOTO 1001
 1000 CONTINUE
!
!**********************************************************************C
!
! **  TIME LOOP COMPLETED
!
      IF(ISCRAY.EQ.0)THEN
        THDMT=THDMT+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        THDMT=T2TMP-TTMP
        WTHDMT=(WT2TMP-WTTMP)*0.001
      ENDIF
!
!**********************************************************************C
!**********************************************************************C
!
! **  CALCULATE VECTOR POTENTIAL AND VECTOR POTENTIAL TRANSPORTS
! **  USING RESULTS OF THE HARMONIC ANALYSIS
!
!     IF(ISVPTHA.NE.1) GOTO 2000
!
!----------------------------------------------------------------------C
!
!     DO K=1,KC
!     DO L=2,LA
!     LS=LSC(L)
!     VPZ(L,K)=TCVP*SUB(L)*SUB(LS)*SVB(L)*SVB(L-1)*HMC(L)*
!    &         ((AMSU(L,K)+AMSU(LS,K))*(AMCV(L,K)+AMCV(L-1,K))
!    &         -(AMCU(L,K)+AMCU(LS,K))*(AMSV(L,K)+AMSV(L-1,K)))
!     ENDDO
!     ENDDO
!
!     DO K=1,KS
!     DO L=2,LA
!     LS=LSC(L)
!     VPX(L,K)=TCVP*((AMSV(L,K)+AMSV(L,K+1))*(AMCW(L,K)+AMCW(LS,K))
!    &              -(AMCV(L,K)+AMCV(L,K+1))*(AMSW(L,K)+AMSW(LS,K)))
!     VPY(L,K)=TCVP*((AMSW(L,K)+AMSW(L-1,K))*(AMCU(L,K)+AMCU(L,K+1))
!    &              -(AMCW(L,K)+AMCW(L-1,K))*(AMSU(L,K)+AMSU(L,K+1)))
!     ENDDO
!     ENDDO
!
!----------------------------------------------------------------------C
!
!     DO K=1,KC
!     DO L=2,LA
!     LS=LSC(L)
!     LN=LNC(L)
!     UVPT(L,K)=(VPZ(LN,K)-VPZ(L,K))/DYU(L)-DZI*(VPY(L,K)-VPY(L,K-1))
!     VVPT(L,K)=DZI*(VPX(L,K)-VPX(L,K-1))-(VPZ(L+1,K)-VPZ(L,K))/DXV(L)
!     ENDDO
!     ENDDO
!
!     DO K=1,KS
!     DO L=2,LA
!     LS=LSC(L)
!     LN=LNC(L)
!     WVPT(L,K)=(VPY(L+1,K)-VPY(L,K))/DXP(L)
!    &         -(VPX(LN,K)-VPX(L,K))/DYP(L)
!     ENDDO
!     ENDDO
!
!----------------------------------------------------------------------C
!
 2000 CONTINUE
!
!**********************************************************************C
!
! **  PRINT FINAL RESULTS
!
      CALL OUTPUT2
!
!**********************************************************************C
!
! **  WRITE RESTART FILE
!
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
!
!**********************************************************************C
!
! **  COMPLETE LEAST SQUARES HARMONIC ANALYSIS
!
      LSLSHA=1
      IF(ISLSHA.EQ.1) CALL LSQHARM
!
!**********************************************************************C
!
! **  OUTPUT COURANT NUMBER DIAGNOSTICS
!
      OPEN(1,FILE='CFLMAX.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='CFLMAX.OUT')
!
      DO L=2,LA
       WRITE(1,1991)IL(L),JL(L),(CFLUUU(L,K),K=1,KC)
       WRITE(1,1992)(CFLVVV(L,K),K=1,KC)
       WRITE(1,1992)(CFLWWW(L,K),K=1,KC)
       WRITE(1,1992)(CFLCAC(L,K),K=1,KC)
      ENDDO
!
      CLOSE(1)
!
 1991 FORMAT(2I5,52F7.2)
 1992 FORMAT(10X,52F7.2)
 1993 FORMAT(2I5,E13.5)
!
!**********************************************************************C
!
! **  OUTPUT COSMETIC VOLUME LOSSES FORM DRY CELLS
!
      IF(NDRYSTP.LT.0) THEN
!
      OPEN(1,FILE='DRYLOSS.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='DRYLOSS.OUT')
!
      DO L=2,LA
       WRITE(1,1993)IL(L),JL(L),VDWASTE(L)
      ENDDO
!
      CLOSE(1)
!
	ENDIF
!
!
! **  OUTPUT FINAL MASS AND VOLUME BALANCES
!
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T5
	  ENDIF
	ENDIF
!
!**********************************************************************C
!
      RETURN
      END
