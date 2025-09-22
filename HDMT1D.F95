!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE HDMT1D
!
! **  SUBROUTINE HDMT1D EXECUTES THE FULL HYDRODYNAMIC AND MASS
! **  TRANSPORT TIME INTERGATION FOR 1D CHANNEL MODE USING A 2TL
! **  INTEGRATION
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
!  added calls to bal2t1, bal2t4, bal2t5
!----------------------------------------------------------------------C
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        TTMP=SECOND( )
        CALL TIMEF(WTTMP)
      ENDIF
!!
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
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
      HMC(L)=0.25*(HMP(L)+HMP(L-1)+HMP(LS)+HMP(LSW))
      UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))+HP(L)*(U(L+1,1)+U(L,1)))*HVI(L)
      U1V(L)=0.25*(H1P(LS)*(U1(LSE,1)+U1(LS,1))+H1P(L)*(U1(L+1,1)+U1(L,1)))*H1VI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,1)+V(L-1,1))+HP(L)*(V(LN,1)+V(L,1)))*HUI(L)
      V1U(L)=0.25*(H1P(L-1)*(V1(LNW,1)+V1(L-1,1))+H1P(L)*(V1(LN,1)+V1(L,1)))*H1UI(L)
      ENDDO
!
!**********************************************************************C
!
! **  INITIALIZE DEPTH, AREA, WETTED PERIMETER AND SURFACE WIDTH
!
      CALL CALAREA(1)
!DIAG      CALL SURFPLT
      CALL CALAREA(2)
!DIAG      CALL SURFPLT
!
      DO L=2,LA
       IF(LCT(L).EQ.6) H2WQ(L)=FADYP(L)/SRFYP(L)
       IF(LCT(L).EQ.7) H2WQ(L)=FADXP(L)/SRFXP(L)
      ENDDO
!
      DO L=2,LA
        FADXP2(L)=FADXP1(L)
        FADYP2(L)=FADYP1(L)
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
!
!     IF(ISWAVE.EQ.1) CALL WAVEBL
!     IF(ISWAVE.EQ.2) CALL WAVESXY
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
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TBX1(L)=(AVCON1*H1UI(L)+STBX(L)*SQRT(V1U(L)*V1U(L)+U1(L,1)*U1(L,1)))*U1(L,1)
       TBY1(L)=(AVCON1*H1VI(L)+STBY(L)*SQRT(U1V(L)*U1V(L)+V1(L,1)*V1(L,1)))*V1(L,1)
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
!
      CALL CALTBXY(ISTL,IS2TL)
!
!**********************************************************************C
!
! **  SET BOTTOM AND SURFACE STRESSES
!
!----------------------------------------------------------------------C
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1)))*U(L,1)
       TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1)))*V(L,1)
       ENDDO
      ENDDO
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
       QQ1(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX1(L))**2+(TVAR3N(L)+TBY1(L))**2)
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
       QQ(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX(L))**2+(TVAR3N(L)+TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2+(TVAR3S(L)+TSY(L))**2)
       ENDDO
      ENDDO
!
!     ENDIF
!
!**********************************************************************C
!
! **   SET SWITCHES FOR TWO TIME LEVEL STEP
!
       ISTL=2
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
!----------------------------------------------------------------------C
!
      DO 1000 N=1,NTS
!
        TIMESEC=(DT*FLOAT(N)+TCON*TBEGIN)
        TIMEDAY=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
!
!DIAG      IF(N.EQ.1.AND.ISDSOLV.GE.1)THEN
!DIAG        OPEN(1,FILE='CHANHDMT.DIA',STATUS='UNKNOWN')
!DIAG        CLOSE(1,STATUS='DELETE')
!DIAG      ENDIF
!
!
      NLOGTMP=2*NLOGTMP
      IF(ILOGC.EQ.NTSMMT)THEN
!        CLOSE(8,STATUS='DELETE')                      !hnr 7/27/2009
!        OPEN(8,FILE='EFDCLOG.OUT',STATUS='UNKNOWN')   !hnr 7/27/2009
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
         IF(N.EQ.1.OR.ISTOPT(0).EQ.0) CALL CALAVB (ISTL)
         IF(N.GT.1.AND.ISTOPT(0).GT.1) CALL CALAVB2 (ISTL)
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
!     IF(ISTL.EQ.3)THEN
!       IF(ISWAVE.EQ.1) CALL WAVEBL
!       IF(ISWAVE.EQ.2) CALL WAVESXY
!     ENDIF
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
      CALL CALEXP1D (ISTL)
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
! **  UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS, CONCENTRATIONS
! **  AND SURFACE ELEVATIONS
!
      CALL CALCSER (ISTL)
      CALL CALQVS (ISTL)
      PSERT(0)=0.
      IF(NPSER.GE.1) CALL CALPSER (ISTL)
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
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
      ENDDO
!
      CALL CALTSXY
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
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
      CALL CALPUV1D(ISTL)
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
!      WRITE(8,17)N,ITER,RSQ,CFMAX,AVMAX,ABMIN,ABMAX,ABMIN       !hnr 7/27/2009
!     CALL FLUSH(8)
      ENDIF
!
!  17 FORMAT('  N,ITER,AVMA,AVMI,ABMA,ABMI',2I5,4(1X,F8.4))
   17 FORMAT('  N,ITER,RSQ,CFMAX,AVMA,AVMI,ABMA,ABMI',I7,I5,2E12.4,4(1X,F8.4))
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
!1D      IF(ISTL.EQ.3)THEN
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
!1D      ENDIF
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
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
      ENDDO
!
      CALL CALTSXY
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
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
!DIAG      OPEN(1,FILE='CHANHDMT.DIA',POSITION='APPEND')
!
!DIAG      WRITE(1,1111)N,ISTL
!DIAG      WRITE(1,1112)
!
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
!
      IF(KC.GT.1)THEN
        CALL CALUVW (ISTL,1)
       ELSE
        DO ND=1,NDM
         LF=2+(ND-1)*LDM
         LL=LF+LDM-1
         DO L=LF,LL
          UHDY(L,1)=UHDYE(L)
!1D          U(L,1)=UHDYE(L)*HUI(L)*DYIU(L)
          U(L,1)=UHDYE(L)/FADYU(L)
          VHDX(L,1)=VHDXE(L)
!1D          V(L,1)=VHDXE(L)*HVI(L)*DXIV(L)
          V(L,1)=VHDXE(L)/FADXV(L)
!DIAG          WRITE(1,1113)IL(L),JL(L),UHDYE(L),UHDY(L,1),U(L,1),
!DIAG     &                 UHDY1E(L),UHDY1(L,1),U1(L,1)
          W(L,1)=0.
         ENDDO
        ENDDO
        CALL CALUVW (ISTL,1)
      ENDIF
!
!
!DIAG      WRITE(1,1111)N,ISTL
!DIAG      WRITE(1,1112)
!DIAG      DO L=2,LA
!DIAG        WRITE(1,1113)IL(L),JL(L),UHDYE(L),UHDY(L,1),U(L,1)
!DIAG      ENDDO
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
!DIAG      CLOSE(1)
!
 1111 FORMAT(2I10)
 1112 FORMAT(' I,J,QX,QX,U,AU')
 1113 FORMAT(2I5,6F10.3)
!
!**********************************************************************C
!
! **  CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS
! **  AT TIME LEVEL (N+1)
!
!----------------------------------------------------------------------C
!
      CALL CALCONC (ISTL,1)
!
!----------------------------------------------------------------------C
!
      IF(KC.EQ.1)THEN
!
      DO K=1,KB
      DO L=1,LC
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
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
      ENDIF
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

      WRITE(6,6001)N
      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN

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

      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN

!      WRITE(8,6004)SALMAX,IMAX,JMAX,KMAX    !hnr 7/27/2009
!      WRITE(8,6005)SALMIN,IMIN,JMIN,KMIN    !hnr 7/27/2009

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

      WRITE(6,6006)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6007)SALMIN,IMIN,JMIN,KMIN

      ENDIF


      IF(ISMMC.EQ.2)THEN

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

      WRITE(6,6001)N
      WRITE(6,6008)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6009)SALMIN,IMIN,JMIN,KMIN

      ENDIF

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
!1D      NTMP=MOD(N,2)
!1D      IF(NTMP.EQ.0.AND.ISTL.EQ.3)THEN
!
! **  CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION
!
!1D      DO ND=1,NDM
!1D       LF=2+(ND-1)*LDM
!1D       LL=LF+LDM-1
!1D       DO L=LF,LL
!1D        UHDY2E(L)=0.
!1D        VHDX2E(L)=0.
!1D        HWQ(L)=0.25*(H2P(L)+2.*H1P(L)+HP(L))
!1D       ENDDO
!1D      ENDDO
!
       K=1
       DO L=2,LA
        UHDYWQ(L,K)=0.5*(UHDY(L,K)+UHDY1(L,K))
        VHDXWQ(L,K)=0.5*(VHDX(L,K)+VHDX1(L,K))
       ENDDO

      DO L=2,LA
       IF(LCT(L).EQ.6) H2WQ(L)=FADYP(L)/SRFYP(L)
       IF(LCT(L).EQ.7) H2WQ(L)=FADXP(L)/SRFXP(L)
      ENDDO

      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3E(L)=UHDY2E(L+1)
       TVAR3N(L)=VHDX2E(LNC(L))
       ENDDO
      ENDDO

      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
        TVAR2E(L,K)=UHDYWQ(L+1   ,K)
        TVAR2N(L,K)=VHDXWQ(LNC(L),K)
        ENDDO
       ENDDO
      ENDDO

      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
        WWQ(L,K)=SWB(L)*(WWQ(L,K-1)-DZC(K)*(TVAR2E(L,K)-UHDYWQ(L,K)-TVAR3E(L)+UHDY2E(L)+TVAR2N(L,K) &
                 -VHDXWQ(L,K)-TVAR3N(L)+VHDX2E(L))*DXYIP(L)+SWB(L)*( QSUM(L,K)-DZC(K)*QSUME(L))     &
                 *DXYIP(L)
        ENDDO
       ENDDO
      ENDDO

      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        HPPTMP=H2WQ(L)+DT2*DXYIP(L)*(QSUME(L)-(TVAR3E(L)-UHDY2E(L)+TVAR3N(L)-VHDX2E(L)))
        HWQ(L)=SPB(L)*HPPTMP+(1.-SPB(L))*HWQ(L)
       ENDDO
      ENDDO
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
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))+DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))+DT2              &
                           *DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
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
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        H2WQ(L)=HWQ(L)
       ENDDO
      ENDDO
!
!1D      ENDIF
      ENDIF
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
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
      UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))+HP(L)*(U(L+1,1)+U(L,1)))*HVI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,1)+V(L-1,1))+HP(L)*(V(LN,1)+V(L,1)))*HUI(L)
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
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,0)=QQ(L,0)+QQ1(L,0)
        QQ2(L,KC)=QQ(L,KC)+QQ1(L,KC)
        QQ1(L,0)=QQ(L,0)
        QQ1(L,KC)=QQ(L,KC)
       ENDDO
      ENDDO
!
      ELSE
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,0)=QQ1(L,0)+QQ1(L,0)
        QQ2(L,KC)=QQ1(L,KC)+QQ1(L,KC)
        QQ1(L,0)=QQ(L,0)
        QQ1(L,KC)=QQ(L,KC)
       ENDDO
      ENDDO
!
      ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE BOTTOM STRESS AT LEVEL (N+1)
!
      CALL CALTBXY(ISTL,IS2TL)
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1)))*U(L,1)
        TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1)))*V(L,1)
       ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
!
!----------------------------------------------------------------------C
!
!     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
!
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
       QQ(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX(L))**2+(TVAR3N(L)+TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2+(TVAR3S(L)+TSY(L))**2)
       ENDDO
      ENDDO
!
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
       IF(ISQQ.EQ.1) CALL CALQQ1 (ISTL)
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
! **  HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
! **  IF NCTBC EQ NTSTBC APPLY TRAPEZOIDAL CORRECTION
!
!----------------------------------------------------------------------C
!
!      IF(NCTBC.EQ.NTSTBC)THEN
!       NCTBC=0
!       ISTL=2
!       DELT=DT
!       DELTD2=0.5*DT
!       DZDDELT=DZ/DELT
!       ROLD=0.5
!       RNEW=0.5
!       GOTO 500
!      ELSE
!       NCTBC=NCTBC+1
!       ISTL=3
!       DELT=DT2
!       DELTD2=DT
!       DZDDELT=DZ/DELT
!       ROLD=0.
!       RNEW=1.
!      ENDIF
!
!**********************************************************************C
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
        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
          IF(NCTMSR.EQ.NWTMSR)THEN
            CALL TMSR
            ICALLTP=1
            NCTMSR=1
           ELSE
            NCTMSR=NCTMSR+1
          ENDIF
        ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  WRITE TO DUMP FILES
!
      IF(ISDUMP.GE.1)THEN                                                                  !HNR_GHD 7/2019 BMD2                     
        IF(TIME.GE.TSDUMP.AND.TIME.LE.TEDUMP)THEN                     !hnr
          IF(NCDUMP.EQ.NSDUMP)THEN                                    !hnr
            CALL DUMP                                                 !hnr
            ICALLTP=1                                                 !hnr
            NCDUMP=1                                                  !hnr
           ELSE                                                       !hnr
            NCDUMP=NCDUMP+1                                           !hnr
          ENDIF                                                       !hnr
        ENDIF                                                         !hnr
      ENDIF                                                           !hnr
!
!**********************************************************************C
!
! **  WRITE TO bmd2 OUTPUT FILE                                         !HNR_GHD 6/2022 BMD2                     
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
      IF(ISDRY.GE.1.AND.ICALLTP.EQ.1)THEN
        OPEN(1,FILE='ZVOLBAL.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        DO LS=1,LORMAX
        IF(VOLZERD.GE.VOLSEL(LS).AND.VOLZERD.LT.VOLSEL(LS+1))THEN
           WTM=VOLSEL(LS+1)-VOLZERD
           WTMP=VOLZERD-VOLSEL(LS)
           DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
           WTM=WTM/DELVOL
           WTMP=WTMP/DELVOL
           SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
           ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
        ENDIF
        ENDDO
        IF(ISDYNSTP.EQ.0)THEN
          TIME=DT*FLOAT(N)+TCON*TBEGIN
          TIME=TIME/TCTMSR
        ELSE
          TIME=TIMESEC/TCTMSR
        ENDIF
        WRITE(1,5304) TIME,SELZERD,ASFZERD,VOLZERD,VETZERD
        CLOSE(1)
      ENDIF
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
      IF(ISSSMMT.NE.2) CALL CALMMT
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
      IF(N.EQ.NCPPH.AND.ISPPH.EQ.1)THEN
       CALL SURFPLT
       NCPPH=NCPPH+(NTSPTC/NPPPH)
      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(N.EQ.NCVPH.AND.ISVPH.GE.1)THEN
       CALL VELPLTH
       NCVPH=NCVPH+(NTSPTC/NPVPH)
      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(N.EQ.NCVPV.AND.ISVPV.GE.1)THEN
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
      IF(N.EQ.NCSPH(1).AND.ISSPH(1).EQ.1)THEN
       IF(ISTRAN(1).GE.1) CALL SALPLTH (1,SAL)
       NCSPH(1)=NCSPH(1)+(NTSPTC/NPSPH(1))
      ENDIF
!
      IF(N.EQ.NCSPH(2).AND.ISSPH(2).EQ.1)THEN
       IF(ISTRAN(2).GE.1) CALL SALPLTH (2,TEM)
       NCSPH(2)=NCSPH(2)+(NTSPTC/NPSPH(2))
      ENDIF
!
      IF(N.EQ.NCSPH(3).AND.ISSPH(3).EQ.1)THEN
       IF(ISTRAN(3).GE.1) CALL SALPLTH (3,DYE)
       NCSPH(3)=NCSPH(3)+(NTSPTC/NPSPH(3))
      ENDIF
!
      IF(N.EQ.NCSPH(4).AND.ISSPH(4).EQ.1)THEN
       IF(ISTRAN(4).GE.1) CALL SALPLTH (4,SFL)
       NCSPH(4)=NCSPH(4)+(NTSPTC/NPSPH(4))
      ENDIF
!
      IF(N.EQ.NCSPH(5).AND.ISSPH(5).EQ.1)THEN
       IF(ISTRAN(5).GE.1) CALL SALPLTH (5,TVAR1S)
       NCSPH(5)=NCSPH(5)+(NTSPTC/NPSPH(5))
      ENDIF
!
      IF(N.EQ.NCSPH(6).AND.ISSPH(6).EQ.1)THEN
       IF(ISTRAN(6).GE.1) CALL SALPLTH (6,SEDT)
       NCSPH(6)=NCSPH(6)+(NTSPTC/NPSPH(6))
      ENDIF
!
      IF(N.EQ.NCSPH(7).AND.ISSPH(7).EQ.1)THEN
       IF(ISTRAN(7).GE.1) CALL SALPLTH (7,SNDT)
       NCSPH(7)=NCSPH(7)+(NTSPTC/NPSPH(7))
      ENDIF
!
!----------------------------------------------------------------------C
!
      DO ITMP=1,7
      IF(N.EQ.NCSPV(ITMP).AND.ISSPV(ITMP).GE.1)THEN
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
      RETURN
      END