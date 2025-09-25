!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE HDMTGVC
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
! 05/01/2002        John Hamrick       05/01/2002       John Hamrick
!  modified calls to calbal and budget subroutines
!
! 12/18/2020        DRBC LZ - ADD CALCULATION OF TURBULENCE DISSIPATION RATE "DISPRATE"
! 01/14/2021        DRBC LZ - PERFORM LOG-LINEAR EXTRAPOLATION FOR "DISPRATE"
! 03/24/2022        DRBC LZ - MODIFY "DISPRATE" FOR 1-D SCENARIO
!
!----------------------------------------------------------------------C
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      REAL :: X_INTERP(3), Y_INTERP(3)  ! DRBC, LZ, 12/12/2020, INTERPOLATION ARRAYS      
!
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        TTMP=SECOND( )
        CALL TIMEF(WTTMP)
      ENDIF
!
      FOURDPI=4./PI
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
      ISTL=3
	IS2TL=0
!
! MODIFIED FOR GVC - JMH 08/03/04
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
!      IF(ISVEG.GE.1.AND.KC.GT.1)
!      DO ND=1,NDM
!       LF=2+(ND-1)*LDM
!       LL=LF+LDM-1
!       DO L=LF,LL
!       TBX1(L)=TBX1(L)+0.5*FXVEG(L,1)*SQRT(V1U(L)*V1U(L)
!     &        +U1(L,1)*U1(L,1))*U1(L,1)
!       TBY1(L)=TBY1(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*U1V(L)
!     &        +V1(L,1)*V1(L,1))*V1(L,1)
!       ENDDO
!      ENDDO
!      ENDIF
!
!**********************************************************************C
!
! **  SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
!
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
!      IF(ISVEG.GE.1.AND.KC.GT.1)
!      DO ND=1,NDM
!       LF=2+(ND-1)*LDM
!       LL=LF+LDM-1
!       DO L=LF,LL
!       TBX(L)=TBX(L)+0.5*FXVEG(L,1)*SQRT(VU(L)*VU(L)
!     &        +U(L,1)*U(L,1))*U(L,1)
!       TBY(L)=TBY(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*UV(L)
!     &        +V(L,1)*V(L,1))*V(L,1)
!       ENDDO
!      ENDDO
!      ENDIF
!
      N=0
      CALL CALTSXY
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
       QQ1(L,KBPM)=0.5*CTURB2*SQRT((TVAR3E(L)+TBX1(L))**2+(TVAR3N(L)+TBY1(L))**2)
       QQ1(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX1(L))**2+(TVAR3S(L)+TSY1(L))**2)
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
       QQ(L,KBPM)=0.5*CTURB2*SQRT((TVAR3E(L)+TBX(L))**2+(TVAR3N(L)+TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2+(TVAR3S(L)+TSY(L))**2)
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
!     ENDIF
!
!**********************************************************************C
!
! **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
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
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       ENDDO
      ENDDO
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         TAUBC2=0.25*( (TVAR3E(L)+TBX1(L))**2+(TVAR3N(L)+TBY1(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U1(L+1,1)+U1(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V1(LN,1)+V1(L,1))
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
! **   SET SWITCHES FOR THREE TIME LEVEL STEP
!
       ISTL=3
       IS2TL=0
       DELT=DT2
       DELTD2=DT
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
!----------------------------------------------------------------------C
!
      DO 1000 N=1,NTS
!
        TIMESEC=(DT*FLOAT(N)+TCON*TBEGIN)
        TIMEDAY=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
!
!      NLOGTMP=2*NLOGTMP
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
      ILOGC=ILOGC+1
!     WRITE(8,6626)ILOGC
 6626 FORMAT(' ILOGC = ',I10)
!
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
!
      IF(IS2TIM.EQ.0) THEN
!
      IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
         CALL CALBAL1
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0)THEN
           CALL CBALEV1
          ELSE
           CALL CBALOD1
         ENDIF
       ENDIF
!
!  ** INITIALIZE SEDIMENT BUDGET CALCULATION   (DLK 10/15)
!
      IF(NCTBC.NE.NTSTBC.AND.ISSBAL.GE.1)THEN
         CALL BUDGET1
!         NTMP=MOD(N,2)
!         IF(NTMP.EQ.0)THEN
!           CALL BUDGEV1
!          ELSE
!           CALL BUDGOD1
!         ENDIF
       ENDIF
!
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
      IF(ISTL.EQ.3)THEN
        IF(ISWAVE.EQ.1) CALL WAVEBL
        IF(ISWAVE.EQ.2) CALL WAVESXY
      ENDIF
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
! NOTES ON VARIOUS VERSIONS OF CALEXP
!
!  CALEXP   -  PRODUCTION VERSION WITH HORIZONTAL MOMENTUM SOURCE
!              AND 3D IMPLICIT VEGETATION DRAG
!                OPERATES WITH CALPUV2
!                OPERATES WITH CALPUV5
!                OPERATES WITH CALPUV9
!                OPERATES WITH CALPUVA
!                OPERATES WITH CALPUVB
!                OPERATES WITH CALPUVC
!                OPERATES WITH CALPUVD
!                OPERATES WITH CALPUVE
!  CALEXP1  -  OLD VERSION OF CALEXP WITH HORIZONTAL DOMAIN
!              DECOMPOSITION FOR SINGLE LAYER, NOT NEEDED FOR LINKING
!  CALEXP2  -  ????????????????????????????
!                OPERATES WITH CALPUV6
!                OPERATES WITH CALPUV7
!                OPERATES WITH CALPUV8
!  CALEXP3  -  SMOLARKIEWCZ-MARGOLIN SCHEME
!                OPERATES WITH CALPUV3
!  CALEXP9  -  PROTOTYPE OF COSMIC ADVECTION SCHEME BUILT FROM CALEXP
!              TO BE REPLACE OF MODIFIED BY COSEXP
!                OPERATES WITH CALPUV2
!                OPERATES WITH CALPUV5
!                OPERATES WITH CALPUV9
!                OPERATES WITH CALPUVA
!                OPERATES WITH CALPUVB
!                OPERATES WITH CALPUVC
!                OPERATES WITH CALPUVD
!                OPERATES WITH CALPUVE
!  CALEXP1D -  TWO TIME LEVEL 1D VERSION CALLED FROM HDMT1D
!  CALEXP2T -  TWO TIME LEVEL 3D VERSION CALLED FROM HDMT2T
!
!X     IF(KC.EQ.1.AND.NDM.GE.2)THEN
!X        CALL CALEXP1 (ISTL)
!X       ELSE
!
! MODIFIED FOR GVC - JMH 08/03/04
!
        IF(ISCDMA.LE.2) CALL CALEXPGVC (ISTL)
!        IF(ISCDMA.EQ.3) CALL CALEXP3 (ISTL)
!        IF(ISCDMA.EQ.4) CALL CALEXP3 (ISTL)
        IF(ISCDMA.EQ.5) CALL CALEXP2 (ISTL)
        IF(ISCDMA.EQ.6) CALL CALEXP2 (ISTL)
        IF(ISCDMA.EQ.9) CALL CALEXP9 (ISTL)
!X      ENDIF
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
      IF(ISTL.EQ.3)THEN
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
      ENDIF
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
! NOTES ON VARIOUS VERSIONS OF CALPUV
!
!  CALPUV    -  RIGID LID SOLVER
!  CALPUVA   -  EXPERIMENTAL
!  CALPUVB   -  EXPERIMENTAL
!  CALPUVC   -  EXPERIMENTAL
!  CALPUVD   -  EXPERIMENTAL
!  CALPUVE   -  EVERGLADES DRYING AND WETTING VERSION
!  CALPUVF   -  EXPERIMENTAL
!  CALPUV2   -  OLD NONDRYING VERSION
!  CALPUV3   -  SMOLARKIEWCZ-MARGOLIN SCHEME
!  CALPUV5   -  PRODUCTION VERSION FOR DRYING AND WETTING
!  CALPUV6   -  EXPERIMENTAL DRYING AND WETTING
!  CALPUV7   -  EXPERIMENTAL DRYING AND WETTING
!  CALPUV8   -  EXPERIMENTAL DRYING AND WETTING
!  CALPUV9   -  PRODUCTION VERSION
!  CALPUV1D  -  TWO TIME LEVEL 1D VERSION CALLED FROM HDMT1D
!  CALPUV2T  -  TWO TIME LEVEL 3D VERSION CALLED FROM HDMT2T
!
!JH      IF(ISDRY.EQ.-1)THEN
!JH         CALL CALPUV(ISTL)
!JH         GOTO 5555
!JH      ENDIF
!
!JH      IF(ISEVER.GE.1)THEN
!JH        CALL CALPUVE(ISTL)
!JH        GOTO 5555
!JH      ENDIF
!
! NOTE THE FOLLOWING BLOCK OF OPTIONS OPERATE WITH CALEXP OR CALEXP9
!
!JH      IF(ISCDMA.LE.2.OR.ISCDMA.GE.9)THEN
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.0) CALL CALPUV2(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.2) CALL CALPUV5(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.3) CALL CALPUV2(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.4) CALL CALPUV5(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.9)THEN
!
! MODIFIED FOR GVC - JMH 08/03/04
!
           IF(ISCHAN.EQ.0.AND.ISDRY.EQ.0) CALL CALPUV9GVC(ISTL)
           IF(ISCHAN.GE.1.OR.ISDRY.GE.1) CALL CALPUV9C(ISTL)
!JH        ENDIF
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.99) CALL CALPUV9(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.10) CALL CALPUVA(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.11) CALL CALPUVB(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.12) CALL CALPUVC(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.13) CALL CALPUVD(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.14) CALL CALPUVF(ISTL)
!JH        IF(ISDRY.EQ.1.OR.ISDRY.EQ.2) CALL CALPUV5(ISTL)
!JH        IF(ISDRY.EQ.11.OR.ISDRY.EQ.12) CALL CALPUV5(ISTL)
!JH        IF(ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV5(ISTL)
!JH        IF(ISDRY.EQ.99) CALL CALPUV5(ISTL)
!JH        GOTO 5555
!JH      ENDIF
!
! NOTE THE FOLLOWING BLOCK OF OPTIONS OPERATES WITH CALEXP3
!
!JH      IF(ISCDMA.GT.2.AND.ISCDMA.LT.5)THEN
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.0) CALL CALPUV3(ISTL)
!JH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.3) CALL CALPUV3(ISTL)
!JH        GOTO 5555
!JH      ENDIF
!
! NOTE THE FOLLOWING BLOCK OF OPTIONS OPERATES WITH CALEXP2
!
!JH      IF(ISCDMA.GE.5.AND.ISCDMA.LE.8)THEN
!JH        IF(ISDRY.EQ.1.OR.ISDRY.EQ.2) CALL CALPUV6(ISTL)
!JH        IF(ISDRY.EQ.11.OR.ISDRY.EQ.12) CALL CALPUV6(ISTL)
!JH        IF(ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV6(ISTL)
!JH        IF(ISDRY.EQ.99) CALL CALPUV7(ISTL)
!X      IF(ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV6(ISTL) !7 MOVED TO 8
!       IF(ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV8(ISTL)
!JH      ENDIF
!
 5555 CONTINUE
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
!      WRITE(8,17)N,ITER,RSQ,CFMAX,AVMAX,ABMIN,ABMAX,ABMIN          !hnr 7/27/2009
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
      IF(ISTL.EQ.3)THEN
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
      ENDIF
!
!**********************************************************************C
!
! **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
! **  MODE FORCING
!
!----------------------------------------------------------------------C
!
      IF(ISCDMA.LE.2.OR.ISCDMA.GE.9)THEN
      IF(ISTL.EQ.3)THEN
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
      ENDIF
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
      IF(KC.GT.1)THEN
!
! MODIFIED FOR GVC - JMH 08/03/04
!
        CALL CALUVWGVC (ISTL,IS2TL)
       ELSE
        DO L=2,LA
          UHDY(L,1)=UHDYE(L)
          U(L,1)=UHDYE(L)*HUI(L)*DYIU(L)
          VHDX(L,1)=VHDXE(L)
          V(L,1)=VHDXE(L)*HVI(L)*DXIV(L)
          W(L,1)=0.
        ENDDO
!
! MODIFIED FOR GVC - JMH 08/03/04
!
        CALL CALUVWGVC (ISTL,IS2TL)
      ENDIF
      IF(ISCRAY.EQ.0)THEN
        TUVW=TUVW+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TUVW=TUVW+T2TMP-T1TMP
        WTUVW=WTUVW+(WT2TMP-WT1TMP)*0.001
      ENDIF
!
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
      ITMP=0
      IF(ISTRAN(4).GE.1) ITMP=1
      IF(ISTRAN(8).GE.1) ITMP=1
      IF(ISICM.GE.1) ITMP=1
!
      IF(ITMP.EQ.1)THEN
      NTMP=MOD(N,2)
      IF(NTMP.EQ.0.AND.ISTL.EQ.3)THEN
!
! **  CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION
!
!jh11/15/06      DO L=2,LA
!jh11/15/06        UHDY2E(L)=0.
!jh11/15/06        VHDX2E(L)=0.
!jh11/15/06        HWQ(L)=0.25*(H2P(L)+2.*H1P(L)+HP(L))
!jh11/15/06      ENDDO
!jh11/15/06C
!jh11/15/06      DO K=1,KC
!jh11/15/06       DO L=2,LA
!jh11/15/06        UHDYWQ(L,K)=0.5*(UHDY1(L,K)+UHDY2(L,K))
!jh11/15/06        VHDXWQ(L,K)=0.5*(VHDX1(L,K)+VHDX2(L,K))
!jh11/15/06        UHDY2E(L)=UHDY2E(L)+GVCSCLU(L)*UHDYWQ(L,K)*DZC(K)
!jh11/15/06        VHDX2E(L)=VHDX2E(L)+GVCSCLV(L)*VHDXWQ(L,K)*DZC(K)
!jh11/15/06        UWQ(L,K)=0.5*(U1(L,K)+U2(L,K))
!jh11/15/06        VWQ(L,K)=0.5*(V1(L,K)+V2(L,K))
!jh11/15/06       ENDDO
!jh11/15/06      ENDDO
!jh11/15/06C
!jh11/15/06      DO L=2,LA
!jh11/15/06       TVAR3E(L)=GVCSCLP(L)*UHDY2E(L+1)
!jh11/15/06       TVAR3N(L)=GVCSCLP(L)*VHDX2E(LNC(L))
!jh11/15/06      ENDDO
!jh11/15/06C
!jh11/15/06      DO K=1,KC
!jh11/15/06      DO L=2,LA
!jh11/15/06        TVAR2E(L,K)=GVCSCLU(L+1   )*UHDYWQ(L+1   ,K)
!jh11/15/06        TVAR2N(L,K)=GVCSCLV(LNC(L))*VHDXWQ(LNC(L),K)
!jh11/15/06      ENDDO
!jh11/15/06      ENDDO
!jh11/15/06C
!jh11/15/06      DO K=1,KS
!jh11/15/06        DO L=2,LA
!jh11/15/06          WWQ(L,K)=SWB3D(L,K)*(WWQ(L,K-1)
!jh11/15/06     &       -DZC(K)*(TVAR2E(L,K)-GVCSCLU(L)*UHDYWQ(L,K)
!jh11/15/06     &               -TVAR3E(L)+GVCSCLP(L)*UHDY2E(L)
!jh11/15/06     &               +TVAR2N(L,K)-GVCSCLV(L)*VHDXWQ(L,K)
!jh11/15/06     &               -TVAR3N(L)+GVCSCLP(L)*VHDX2E(L))*DXYIP(L))
!jh11/15/06     &    +SWB3D(L,K)*( QSUM(L,K)-DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
!jh11/15/06        ENDDO
!jh11/15/06      ENDDO
!jh11/15/06C
!jh11/15/06      IF(ISDIVEX.EQ.3)THEN
!jh11/15/06	  OPEN(1,FILE='DIV3GVC.OUT')
!jh11/15/06	  CLOSE(1,STATUS='DELETE')
!jh11/15/06	  OPEN(1,FILE='DIV3GVC.OUT')
!jh11/15/06        K=KC
!jh11/15/06        DO L=2,LA
!jh11/15/06          TVAR3S(L)=(WWQ(L,K-1)
!jh11/15/06     &       -DZC(K)*(TVAR2E(L,K)-GVCSCLU(L)*UHDYWQ(L,K)
!jh11/15/06     &               -TVAR3E(L)+GVCSCLP(L)*UHDY2E(L)
!jh11/15/06     &               +TVAR2N(L,K)-GVCSCLV(L)*VHDXWQ(L,K)
!jh11/15/06     &               -TVAR3N(L)+GVCSCLP(L)*VHDX2E(L))*DXYIP(L))
!jh11/15/06     &    +( QSUM(L,K)-DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
!jh11/15/06        ENDDO
!jh11/15/06	  DO L=2,LA
!jh11/15/06	    WRITE(1,1492)IL(L),JL(L),TVAR3S(L)
!jh11/15/06	  ENDDO
!jh11/15/06	  CLOSE(1)
!jh11/15/06	ENDIF
!jh11/15/06c
!jh11/15/06      DO L=2,LA
!jh11/15/06       TVAR3E(L)=UHDY2E(L+1)
!jh11/15/06       TVAR3N(L)=VHDX2E(LNC(L))
!jh11/15/06      ENDDO
!jh11/15/06C
!jh11/15/06      DO L=2,LA
!jh11/15/06        HPPTMP=H2WQ(L)+DT2*DXYIP(L)*(QSUME(L)
!jh11/15/06     &      -(TVAR3E(L)-UHDY2E(L)+TVAR3N(L)-VHDX2E(L)))
!jh11/15/06        HWQ(L)=SPB(L)*HPPTMP+(1.-SPB(L))*HWQ(L)
!jh11/15/06      ENDDO
!
!jh11/15/06 - temporarily modified to use 3t hydro transport
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
          TVAR3S(L)=(WWQ(L,K-1)-DZC(K)*(TVAR2E(L,K)-GVCSCLU(L)*UHDYWQ(L,K)-TVAR3E(L)+GVCSCLP(L)   &
                    *UHDY2E(L)+TVAR2N(L,K)-GVCSCLV(L)*VHDXWQ(L,K)-TVAR3N(L)+GVCSCLP(L)*VHDX2E(L)) &
                    *DXYIP(L))+( QSUM(L,K)-DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
        ENDDO
	  DO L=2,LA
	    WRITE(1,1492)IL(L),JL(L),TVAR3S(L)
	  ENDDO
	  CLOSE(1)
	ENDIF
!
!jh11/15/06
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
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))+DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))                  &
          +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
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
      ENDIF
!
 1492 FORMAT(2I6,E15.7)
!
!**********************************************************************C
!
! **  UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING
! **  AN EQUATION OF STATE
!
      IF(ISTL.EQ.3)THEN
       DO K=1,KC
       DO L=2,LA
       B1(L,K)=B(L,K)
       ENDDO
       ENDDO
      ENDIF
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
      IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
         CALL CALBAL4
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0)THEN
           CALL CBALEV4
          ELSE
           CALL CBALOD4
         ENDIF
      ENDIF
!
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
      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
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
!
!GVCDIAG
!      IF(N.LE.6)THEN
!      DO L=2,LA
!	  WRITE(8,881)N,IL(L),JL(L),KGVCU(L),TBX(L),STBX(L),
!     &             U(L,KGVCU(L)),VU(L)
!      ENDDO
!      ENDIF
!  881 FORMAT('TBX GVC ',4I5,10E14.6)
!GVCDIAG
!
!
!      IF(ISVEG.GE.1.AND.KC.GT.1)
!      DO ND=1,NDM
!       LF=2+(ND-1)*LDM
!       LL=LF+LDM-1
!       DO L=LF,LL
!       TBX(L)=TBX(L)+0.5*FXVEG(L,1)*SQRT(VU(L)*VU(L)
!     &        +U(L,1)*U(L,1))*U(L,1)
!       TBY(L)=TBY(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*UV(L)
!     &        +V(L,1)*V(L,1))*V(L,1)
!       ENDDO
!      ENDDO
!      ENDIF
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
       QQ(L,KBPM)=0.5*CTURB2*SQRT((TVAR3E(L)+TBX(L))**2+(TVAR3N(L)+TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2+(TVAR3S(L)+TSY(L))**2)
      ENDDO
!
      ENDIF
!
!     ENDIF
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
	    IF(ISTOPT(0).EQ.0)CALL CALQQ1OLD (ISTL)
!
! MODIFIED FOR GVC - JMH 08/03/04
!
	    IF(ISTOPT(0).GE.1)CALL CALQQ1GVC (ISTL)
        ENDIF
        IF(ISQQ.EQ.2) CALL CALQQ2 (ISTL)
      ENDIF
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
! **  CALCULATE TURBULENCE DISSIPATION RATE, 12/16/2020, 3/24/2022, DRBC LZ
!
      DO L=2, LA
      IF(LMASKDRY(L))THEN  !WATER CELL
      
        IF (KGVCP(L).EQ.KC.AND.LGVCP(L,KC)) THEN  !ONE ACTIVE LAYER, USING LOG-LAW
          LN=LNC(L)
          Ubmd=0.5*(TBX(L+1)+TBX(L))
	  Vbmd=0.5*(TBY(LN)+TBY(L))
          xx=(ubmd*ubmd+vbmd*vbmd)**0.5  
          DISPRATE(L,KC)=xx**1.5/(VKC*HP(L))
        ELSE
          DO K=1,KS        !TWO OR MORE ACTIVE LAYERS, USING TURBULENT MODEL OUTPUT
            IF(LGVCP(L,K))THEN
              DISPRATE(L,K)=QQ(L,K)**1.5/(CTURB*DML(L,K))
            ENDIF
          ENDDO
          
          !EXTRAPOLATE TKE DISSIPATION RATE AT KC, I.E., AIR-WATER INTERFACE
          IF ((KC-KGVCP(L)).EQ.1) THEN         !TWO ACTIVE LAYERS
            DISPRATE(L,KC)=DISPRATE(L,KS)      !APPROXIMATE DISSIPATION RATE IN WATER COLUMN TO AIR-WATER INTERFACE
 
          ELSE                                 !THREE OR MORE ACTIVE LAYERS   
            IORDER=1                           !LINEAR EXTRAPOLATION
            XSUM=0.0
            X_INTERP=0.0
            Y_INTERP=0.0
            
            DO KK=1,IORDER+1
              IF(LGVCP(L,KC-KK))THEN
                XSUM=XSUM+DZC(KC-KK+1)*GVCSCLP(L)  !LAYER VOLUME MODIFIED FOR GVC
                X_INTERP(KK)=XSUM                  !X&Y ARRAYS ARE POSITIVE DOWNWARDS FROM AIR-WATER INTERFACE
                Y_INTERP(KK)=LOG10(DISPRATE(L,KC-KK))  !LOG-LINEAR EXTRAPOLATION TO AVOID NEGATIVE VALUE
              ENDIF
            ENDDO
            
            DISPRTMP=alagrange(X_INTERP,Y_INTERP,IORDER) !EXTRAPOLATION AT X=0 USING FUNCTION LAGRANGE
            DISPRATE(L,KC)=10.0**(DISPRTMP)
          ENDIF
            
        ENDIF
      
      ENDIF
      ENDDO
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
      IF(NCTBC.EQ.NTSTBC)THEN
       NCTBC=0
       ISTL=2
       DELT=DT
       DELTD2=0.5*DT
       DZDDELT=DZ/DELT
       ROLD=0.5
       RNEW=0.5
       GOTO 500
      ELSE
       NCTBC=NCTBC+1
       ISTL=3
       DELT=DT2
       DELTD2=DT
       DZDDELT=DZ/DELT
       ROLD=0.
       RNEW=1.
      ENDIF
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
      IF(ISTMSR.GE.1)THEN
!        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
          IF(NCTMSR.GE.NWTMSR)THEN
            CALL TMSR
            ICALLTP=1
            NCTMSR=1
           ELSE
            NCTMSR=NCTMSR+1
          ENDIF
!        ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  WRITE TO DUMP FILES
!
      IF(ISDUMP.GE.1)THEN                                                !HNR_GHD 7/2019 BMD2                     
        IF(TIME.GE.TSDUMP.AND.TIME.LE.TEDUMP)THEN         !hnr
          IF(NCDUMP.GE.NSDUMP)THEN                        !hnr
            CALL DUMP                                     !hnr
            ICALLTP=1                                     !hnr
            NCDUMP=1                                      !hnr
           ELSE                                           !hnr
            NCDUMP=NCDUMP+1                               !hnr
          ENDIF                                           !hnr
        ENDIF                                             !hnr
      ENDIF                                               !hnr
!                                                          !hnr
!**********************************************************************C
!
! **  WRITE TO bmd2 OUTPUT FILE                                                 !HNR_GHD 7/2019 BMD2   !HNR_GHD 6/2022 remove bmd file keep only bmd2                  
!                                                                 
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
!                                                          !hnr
!**********************************************************************C
!
!**********************************************************************C
!
! ** WRITE BOTTOM VELOCITIES TO FILE FOR AJ MINE SEDIMENT BED MODEL
!     (DLK - ADDED 3/25/97)
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
      IF(ISBAL.GE.1)THEN
         CALL CALBAL5
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0)THEN
           CALL CBALEV5
          ELSE
           CALL CBALOD5
         ENDIF
       ENDIF
!
!   SEDIMENT BUDGET CALCULATION     (DLK 10/15)
!
       IF(ISSBAL.GE.1)THEN
       CALL BUDGET5
       ENDIF
!       NTMP=MOD(N,2)
!       IF(NTMP.EQ.0)THEN
!         CALL BUDGEV5
!        ELSE
!         CALL BUDGOD5
!       ENDIF
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
      IF(N.GE.NCPPH.AND.ISPPH.GE.1)THEN
       CALL SURFPLT
       NCPPH=NCPPH+(NTSPTC/NPPPH)
      ENDIF
!
!----------------------------------------------------------------------C
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
!----------------------------------------------------------------------C
!
      DO ITMP=1,7
      IF(N.GE.NCSPV(ITMP).AND.ISSPV(ITMP).GE.1)THEN
       CALL SALPLTV(ITMP)
       NCSPV(ITMP)=NCSPV(ITMP)+(NTSPTC/NPSPV(ITMP))
      ENDIF
      ENDDO
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
!
!**********************************************************************C
!
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
!
!**********************************************************************C
!
      RETURN
      END