!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALEXP (ISTL)
!
! **  SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
! **  THIS SUBROUTINE IS CURRENT PRODUCTION VERSION
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      DIMENSION TMPVEC1(KCM),TMPVEC2(KCM),FXTMP(KCM),FYTMP(KCM)
      DIMENSION QMCSINKX(LCM,KCM),QMCSOURX(LCM,KCM),QMCSINKY(LCM,KCM),QMCSOURY(LCM,KCM)
      DIMENSION FUHJ(LCM,KCM),FVHJ(LCM,KCM),DZPC(LCM,KCM)
!
!**********************************************************************C
!
       DELT=DT2
       DELTD2=DT
      IF(ISTL.EQ.2)THEN
       DELT=DT
       DELTD2=0.5*DT
      ENDIF
!
      DELTI=1./DELT
!
      IF(N.EQ.1)THEN
        OPEN(1,FILE='MFLUX.DIA')
        CLOSE(1,STATUS='DELETE')
      ENDIF
!
      IF(N.LE.4)THEN
        OPEN(1,FILE='MFLUX.DIA',POSITION='APPEND')
      ENDIF
!
!**********************************************************************C
!
! **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
!
!----------------------------------------------------------------------C
!
      FCAXE(1)=0.
      FCAYE(1)=0.
      FCAX1E(1)=0.
      FCAY1E(1)=0.
      FXE(1)=0.
      FYE(1)=0.
!
      FCAXE(LC)=0.
      FCAYE(LC)=0.
      FCAX1E(LC)=0.
      FCAY1E(LC)=0.
      FXE(LC)=0.
      FYE(LC)=0.
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        FCAXE(L)=0.
        FCAYE(L)=0.
        FCAX1E(L)=0.
        FCAY1E(L)=0.
        FXE(L)=0.
        FYE(L)=0.
       ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  SELECT ADVECTIVE FLUX FORM
!
!----------------------------------------------------------------------C
!
      IF(ISTL.EQ.2) GOTO 200
      IF(ISTL.EQ.3)THEN
       IF(ISCDMA.EQ.0) GOTO 300
       IF(ISCDMA.EQ.1) GOTO 400
       IF(ISCDMA.EQ.2) GOTO 350
      ENDIF
!
!**********************************************************************C
!
! **  TWO TIME LEVEL STEP
! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
! **  AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N 
!
  200 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=2,LA
        UHDY2(L,K)=0.5*(UHDY1(L,K)+UHDY(L,K))
        VHDX2(L,K)=0.5*(VHDX1(L,K)+VHDX(L,K))
        W2(L,K)=W1(L,K)+W(L,K)
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.5*(UHDY2(L,K)+UHDY2(LS,K))
       UHB=0.5*(UHDY2(L,K)+UHDY2(L+1,K))
       VHC=0.5*(VHDX2(L,K)+VHDX2(L-1,K))
       VHB=0.5*(VHDX2(L,K)+VHDX2(LN,K))
!
       FUHU(L,K)=MAX(UHB,0.)*U1(L,K)+MIN(UHB,0.)*U1(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)+MIN(VHC,0.)*U1(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)+MIN(UHC,0.)*V1(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V1(L,K)+MIN(VHB,0.)*V1(LN,K)
       FUHJ(L,K)=0.
       FVHJ(L,K)=0.
       ENDDO
      ENDDO
!
! ADD RETURN FLOW MOMENTUM FLUX
!
      DO NWR=1,NQWR
       IF(NQWRMFU(NWR).GT.0)THEN
         IU=IQWRU(NWR)
         JU=JQWRU(NWR)
         KU=KQWRU(NWR)
         LU=LIJ(IU,JU)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=QMF*QMF/(H1P(LU)*DZC(KU)*DZC(KU)*BQWRMFU(NWR))
         IF(NQWRMFU(NWR).EQ.1)  FUHJ(LU     ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.2)  FVHJ(LU     ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.3)  FUHJ(LU+1   ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.4)  FVHJ(LNC(LU),KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.-1) FUHJ(LU     ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-2) FVHJ(LU     ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-3) FUHJ(LU+1   ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-4) FVHJ(LNC(LU),KU)=QUMF 
       ENDIF
       IF(NQWRMFD(NWR).GT.0)THEN
         ID=IQWRD(NWR)
         JD=JQWRD(NWR)
         KD=KQWRD(NWR)
         LD=LIJ(ID,JD)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=QMF*QMF/(H1P(LD)*DZC(KD)*DZC(KD)*BQWRMFD(NWR))
         IF(NQWRMFD(NWR).EQ.1)  FUHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.2)  FVHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.3)  FUHJ(LD+1   ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.4)  FVHJ(LNC(LD),KD)=QUMF
         IF(NQWRMFD(NWR).EQ.-1) FUHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-2) FVHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-3) FUHJ(LD+1   ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-4) FVHJ(LNC(LD),KD)=-QUMF
!         IF(N.LE.4)THEN
!           WRITE(1,1112)N,NWR,NS,ID,JD,KD,NQWRMFD(NWR),H1P(LD),QMF,
!     &                  QUMF,FUHJ(LD,KD),FVHJ(LD,KD)
!         ENDIF
       ENDIF
      ENDDO
!
! ** HARDWIRE FOR PEACH BOTTOM
!
!      DO K=1,KC
!       FVHV(535,K)=700./H1P(535)
!      ENDDO
!
! ** END HARDWIRE FOR PEACH BOTTOM
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.25*DXYU(L)*(W2(L,K)+W2(L-1,K))
       WV=0.25*DXYV(L)*(W2(L,K)+W2(LS,K))
       FWU(L,K)=MAX(WU,0.)*U1(L,K)+MIN(WU,0.)*U1(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V1(L,K)+MIN(WV,0.)*V1(L,K+1)
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      GOTO 500
!
!**********************************************************************C
!
! **  THREE TIME LEVEL (LEAP-FROG) STEP
! **  WITH TRANSPORT AT (N) AND TRANSPORTED FIELD AT (N-1)
!
  300 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.5*(UHDY(L,K)+UHDY(LS,K))
       UHB=0.5*(UHDY(L,K)+UHDY(L+1,K))
       VHC=0.5*(VHDX(L,K)+VHDX(L-1,K))
       VHB=0.5*(VHDX(L,K)+VHDX(LN,K))
       FUHU(L,K)=MAX(UHB,0.)*U1(L,K)+MIN(UHB,0.)*U1(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)+MIN(VHC,0.)*U1(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)+MIN(UHC,0.)*V1(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V1(L,K)+MIN(VHB,0.)*V1(LN,K)
       FUHJ(L,K)=0.
       FVHJ(L,K)=0.
       ENDDO
      ENDDO
!
! ADD RETURN FLOW MOMENTUM FLUX
!
      DO NWR=1,NQWR
       IF(NQWRMFU(NWR).GT.0)THEN
         IU=IQWRU(NWR)
         JU=JQWRU(NWR)
         KU=KQWRU(NWR)
         LU=LIJ(IU,JU)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=QMF*QMF/(H1P(LU)*DZC(KU)*DZC(KU)*BQWRMFU(NWR))
         IF(NQWRMFU(NWR).EQ.1)  FUHJ(LU     ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.2)  FVHJ(LU     ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.3)  FUHJ(LU+1   ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.4)  FVHJ(LNC(LU),KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.-1) FUHJ(LU     ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-2) FVHJ(LU     ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-3) FUHJ(LU+1   ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-4) FVHJ(LNC(LU),KU)=QUMF 
       ENDIF
       IF(NQWRMFD(NWR).GT.0)THEN
         ID=IQWRD(NWR)
         JD=JQWRD(NWR)
         KD=KQWRD(NWR)
         LD=LIJ(ID,JD)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=QMF*QMF/(H1P(LD)*DZC(KD)*DZC(KD)*BQWRMFD(NWR))
         IF(NQWRMFD(NWR).EQ.1)  FUHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.2)  FVHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.3)  FUHJ(LD+1   ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.4)  FVHJ(LNC(LD),KD)=QUMF
         IF(NQWRMFD(NWR).EQ.-1) FUHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-2) FVHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-3) FUHJ(LD+1   ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-4) FVHJ(LNC(LD),KD)=-QUMF
!         IF(N.LE.4)THEN
!           WRITE(1,1112)N,NWR,NS,ID,JD,KD,NQWRMFD(NWR),H1P(LD),QMF,
!     &                  QUMF,FUHJ(LD,KD),FVHJ(LD,KD)
!         ENDIF
       ENDIF
      ENDDO
!
! ** HARDWIRE FOR PEACH BOTTOM
!
!      DO K=1,KC
!       FVHV(535,K)=700./H1P(535)
!      ENDDO
!
! ** END HARDWIRE FOR PEACH BOTTOM
!
!----------------------------------------------------------------------C
C
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.5*DXYU(L)*(W(L,K)+W(L-1,K))
       WV=0.5*DXYV(L)*(W(L,K)+W(LS,K))
       FWU(L,K)=MAX(WU,0.)*U1(L,K)+MIN(WU,0.)*U1(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V1(L,K)+MIN(WV,0.)*V1(L,K+1)
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      GOTO 500
!
!**********************************************************************C
!
! **  THREE TIME LEVEL (LEAP-FROG) STEP
! **  FIRST HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE
! **  WITH TRANSPORT AT (N-1/2) AND TRANSPORTED FIELD AT (N-1)
! **  SECOND HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE  
! **  WITH TRANSPORT AT (N+1/2) AND TRANSPORTED FIELD AT (N)
!
  350 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=2,LA
       U2(L,K)=U1(L,K)+U(L,K)
       V2(L,K)=V1(L,K)+V(L,K)
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.25*(UHDY(L,K)+UHDY(LS,K))
       UHB=0.25*(UHDY(L,K)+UHDY(L+1,K))
       VHC=0.25*(VHDX(L,K)+VHDX(L-1,K))
       VHB=0.25*(VHDX(L,K)+VHDX(LN,K))
       FUHU(L,K)=MAX(UHB,0.)*U2(L,K)+MIN(UHB,0.)*U2(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U2(LS,K)+MIN(VHC,0.)*U2(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V2(L-1,K)+MIN(UHC,0.)*V2(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V2(L,K)+MIN(VHB,0.)*V2(LN,K)
       FUHJ(L,K)=0.
       FVHJ(L,K)=0.
       ENDDO
      ENDDO
!
!     DO K=1,KC
!     DO L=2,LA
!     LN=LNC(L)
!     LS=LSC(L)               
!     UHC=0.125*(UHDY(L,K)+UHDY(LS,K)+UHDY1(L,K)+UHDY1(LS,K))
!     UHB=0.125*(UHDY(L,K)+UHDY(L+1,K)+UHDY1(L,K)+UHDY1(L+1,K))
!     VHC=0.125*(VHDX(L,K)+VHDX(L-1,K)+VHDX1(L,K)+VHDX1(L-1,K))
!     VHB=0.125*(VHDX(L,K)+VHDX(LN,K)+VHDX1(L,K)+VHDX1(LN,K))
!     FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
!    &         +MIN(UHB,0.)*U1(L+1,K)
!     FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)
!    &         +MIN(VHC,0.)*U1(L,K)
!     FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)
!    &         +MIN(UHC,0.)*V1(L,K)
!     FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
!    &         +MIN(VHB,0.)*V1(LN,K)
!     ENDDO
!     ENDDO
!
!     DO K=1,KC
!     DO L=2,LA
!     LN=LNC(L)
!     LS=LSC(L)               
!     UHC=0.125*(3.*UHDY(L,K)+3.*UHDY(LS,K)-UHDY1(L,K)-UHDY1(LS,K))
!     UHB=0.125*(3.*UHDY(L,K)+3.*UHDY(L+1,K)-UHDY1(L,K)-UHDY1(L+1,K))
!     VHC=0.125*(3.*VHDX(L,K)+3.*VHDX(L-1,K)-VHDX1(L,K)-VHDX1(L-1,K))
!     VHB=0.125*(3.*VHDX(L,K)+3.*VHDX(LN,K)-VHDX1(L,K)-VHDX1(LN,K))
!     FUHU(L,K)=FUHU(L,K)+MAX(UHB,0.)*U(L,K)
!    &                   +MIN(UHB,0.)*U(L+1,K)
!     FVHU(L,K)=FVHU(L,K)+MAX(VHC,0.)*U(LS,K)
!    &                   +MIN(VHC,0.)*U(L,K)
!     FUHV(L,K)=FUHV(L,K)+MAX(UHC,0.)*V(L-1,K)
!    &                   +MIN(UHC,0.)*V(L,K)
!     FVHV(L,K)=FVHV(L,K)+MAX(VHB,0.)*V(L,K)
!    &                   +MIN(VHB,0.)*V(LN,K)
!     ENDDO
!     ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.25*DXYU(L)*(W(L,K)+W(L-1,K))
       WV=0.25*DXYV(L)*(W(L,K)+W(LS,K))
       FWU(L,K)=MAX(WU,0.)*U2(L,K)+MIN(WU,0.)*U2(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V2(L,K)+MIN(WV,0.)*V2(L,K+1)
       ENDDO
      ENDDO
!
!     DO K=1,KS
!     DO L=2,LA
!     LS=LSC(L)
!     WU=0.125*DXYU(L)*(W(L,K)+W(L-1,K)+W1(L,K)+W1(L-1,K))
!     WV=0.125*DXYV(L)*(W(L,K)+W(LS,K)+W1(L,K)+W1(LS,K))
!     FWU(L,K)=MAX(WU,0.)*U1(L,K)
!    &        +MIN(WU,0.)*U1(L,K+1)
!     FWV(L,K)=MAX(WV,0.)*V1(L,K)
!    &        +MIN(WV,0.)*V1(L,K+1)
!     ENDDO
!     ENDDO
!
!     DO K=1,KS
!     DO L=2,LA
!     LS=LSC(L)
!     WU=0.125*DXYU(L)*(3.*W(L,K)+3.*W(L-1,K)-W1(L,K)-W1(L-1,K))
!     WV=0.125*DXYV(L)*(3.*W(L,K)+3.*W(LS,K)-W1(L,K)-W1(LS,K))
!     FWU(L,K)=FWU(L,K)+MAX(WU,0.)*U(L,K)
!    &                 +MIN(WU,0.)*U(L,K+1)
!     FWV(L,K)=FWV(L,K)+MAX(WV,0.)*V(L,K)
!    &                 +MIN(WV,0.)*V(L,K+1)
!     ENDDO
!     ENDDO
!
!----------------------------------------------------------------------C
!
      GOTO 500
!
!**********************************************************************C
!
! **  THREE TIME LEVEL (LEAP-FROG) STEP
! **  CALCULATE ADVECTIVE FLUXES BY CENTRAL DIFFERENCE WITH TRANSPORT 
! **  AT (N) AND TRANSPORTED FIELD AT (N)
!
  400 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)  
      FUHU(L,K)=0.25*(UHDY(L+1,K)+UHDY(L,K))*(U(L+1,K)+U(L,K))
      FVHU(L,K)=0.25*(VHDX(L,K)+VHDX(L-1,K))*(U(L,K)+U(LS,K))
      FUHV(L,K)=0.25*(UHDY(L,K)+UHDY(LS,K))*(V(L,K)+V(L-1,K))
      FVHV(L,K)=0.25*(VHDX(L,K)+VHDX(LN,K))*(V(LN,K)+V(L,K))
      FUHJ(L,K)=0.0
      FVHJ(L,K)=0.0
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)    
      FWU(L,K)=0.25*DXYU(L)*(W(L,K)+W(L-1,K))*(U(L,K+1)+U(L,K))
      FWV(L,K)=0.25*DXYV(L)*(W(L,K)+W(LS,K))*(V(L,K+1)+V(L,K))
      ENDDO
      ENDDO
!
!**********************************************************************C
!
  500 CONTINUE
!
      DO K=1,KC
       DO L=1,LA
       FUHU(L,K)=STCUV(L)*FUHU(L,K)
       FVHV(L,K)=STCUV(L)*FVHV(L,K)
       ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS
!
!----------------------------------------------------------------------C
!
      IF(ISDCCA.EQ.0)THEN
!
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)          
!       CAC(L,K)=( FCORC(L)*DXYP(L)
!     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
!     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       CAC(L,K)=( FCORC(L)*DXYP(L)+0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L)-0.5*SNLT*(U(L+1,K)+U(L,K))*  &
                DXDJ(L) )*HP(L)
       ENDDO
      ENDDO
!
      ELSE
!
      CFMAX=CF
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)          
!       CAC(L,K)=( FCORC(L)*DXYP(L)
!     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
!     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       CAC(L,K)=(FCORC(L)*DXYP(L)+0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L)-0.5*SNLT*(U(L+1,K)+U(L,K))*  &
                DXDJ(L))*HP(L)
       CFEFF=ABS(CAC(L,K))*DXYIP(L)*HPI(L)
       CFMAX=MAX(CFMAX,CFEFF)
       ENDDO
      ENDDO
!
      ENDIF
!
!      IF(N.EQ.2)THEN
!       OPEN(1,FILE='CORC.DIA')
!       CLOSE(2,STATUS='DELETE')
!       OPEN(1,FILE='CORC.DIA')
!       K=1
!       DO L=2,LA
!       LN=LNC(L)          
!       WRITE(1,1111)IL(L),JL(L),LN,V(LN,K),V(L,K),DYU(L+1),DYU(L),
!     &        U(L+1,K),U(L,K),DXV(LN),DXV(L),HP(L),CAC(L,K)
!       ENDDO
!       CLOSE(1)
!      ENDIF
 1111 FORMAT(3I5,10E13.4)
!
!**********************************************************************C
!
! **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS 
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=1,LC
        FCAX1(L,K)=0.
        FCAY1(L,K)=0.     
       ENDDO
      ENDDO
!
      IF(ISTL.EQ.3)THEN
!      IF(IRVEC.EQ.10.OR.IRVEC.EQ.11)THEN
      IF(IRVEC.EQ.14)THEN
!
      DO L=1,LC
       TVAR3E(L)=0.
       TVAR3N(L)=0.
       TVAR3W(L)=0.
       TVAR3S(L)=0.
       TVAR1S(L,1)=0.
       TVAR2S(L,1)=0.
      ENDDO
!
      DO K=1,KC
      DZTMP=DZC(K) 
      DO L=2,LA
       TVAR3E(L)=TVAR3E(L)+DZTMP*U(L,K)
       TVAR3N(L)=TVAR3N(L)+DZTMP*V(L,K)
       TVAR3W(L)=TVAR3W(L)+DZTMP*U1(L,K)
       TVAR3S(L)=TVAR3S(L)+DZTMP*V1(L,K)
      ENDDO
      ENDDO
!
!       CAC(L,K)=( FCORC(L)*DXYP(L)
!     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
!     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
!
!       TVAR1S=CAC @ N
!       TVAR2S=CAC @ N-1
!
!       TVAR3E=U @ N
!       TVAR3N=V @ N
!       TVAR3W=U @ N-1
!       TVAR3S=V @ N-1
!
      DO L=2,LA
       LN=LNC(L)          
       TVAR1S(L,1)=(FCORC(L)*DXYP(L)+0.5*SNLT*(TVAR3N(LN)+TVAR3N(L))*(DYU(L+1)-DYU(L))-0.5*SNLT* &
                   (TVAR3E(L+1)+TVAR3E(L))*(DXV(LN)-DXV(L)) )*HP(L)
       TVAR2S(L,1)=(FCORC(L)*DXYP(L)+0.5*SNLT*(TVAR3S(LN)+TVAR3S(L))*(DYU(L+1)-DYU(L))-0.5*SNLT* &
                   (TVAR3W(L+1)+TVAR3W(L))*(DXV(LN)-DXV(L)) )*H1P(L)
      ENDDO
!
!       FCAX(L,K)=ROLD*FCAX(L,K)
!     &          +0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
!     &                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
!       FCAY(L,K)=ROLD*FCAY(L,K)
!     &          +0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
!     &                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
!
!
!       TVAR1S=CAC @ N
!       TVAR2S=CAC @ N-1
!
!       TVAR3E=U @ N
!       TVAR3N=V @ N
!       TVAR3W=U @ N-1
!       TVAR3S=V @ N-1
!
      DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX1E(L)=0.0625*SCAX(L)*(TVAR2S(L,1)*(TVAR3S(LN)+TVAR3S(L))+TVAR2S(L-1,1)*(TVAR3S(LNW)+   &
                 TVAR3S(L-1)) )
       FCAY1E(L)=0.0625*SCAY(L)*(TVAR2S(L,1)*(TVAR3W(L+1)+TVAR3W(L))+TVAR2S(LS,1)*(TVAR3W(LSE)+   &
                 TVAR3W(LS)) )
      ENDDO
!
      DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX1E(L)=FCAX1E(L)-0.125*SCAX(L)*(TVAR1S(L,1)*(TVAR3N(LN)+TVAR3N(L))+TVAR1S(L-1,1)*       &
                (TVAR3N(LNW)+TVAR3N(L-1)))
       FCAY1E(L)=FCAY1E(L)-0.125*SCAY(L)*(TVAR1S(L,1)*(TVAR3E(L+1)+TVAR3E(L))+TVAR1S(LS,1)*       &
                 (TVAR3E(LSE)+TVAR3E(LS)) )
      ENDDO
!
      ENDIF
      ENDIF
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX(L,K)=ROLD*FCAX(L,K)+0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))+CAC(L-1,K)*(V(LNW,K)+&
                 V(L-1,K)))
       FCAY(L,K)=ROLD*FCAY(L,K)+0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))+CAC(LS,K)*(U(LSE,K)+&
                 U(LS,K)))
       FX(L,K)=SAAX(L)*(FUHU(L,K)-FUHU(L-1,K)+FVHU(LN,K)-FVHU(L,K)+FUHJ(L,K) )
       FY(L,K)=SAAY(L)*(FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K)+FVHJ(L,K) )
       ENDDO
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
       FCAX(L,K)=SUB(L)*FCAX(L,K)
       FCAY(L,K)=SVB(L)*FCAY(L,K)
       FX(L,K)=SUB(L)*FX(L,K)
       FY(L,K)=SVB(L)*FY(L,K)
       ENDDO
      ENDDO
!
!      IF(N.EQ.2)THEN
!       OPEN(1,FILE='CORC.DIA',POSITION='APPEND')
!        L=2
!        K=1
!        LN=LNC(L)
!        LNW=LNWC(L)
!        WRITE(1,1111)IL(L),JL(L),LNW,SCAX(L),CAC(L,K),V(LN,K),V(L,K),
!     &                              CAC(L-1,K),V(LNW,K),V(L-1,K)
!       DO L=2,LA
!       WRITE(1,1111)L,IL(L),JL(L),(FCAX(L,K),K=1,KC)
!       ENDDO
!       DO L=2,LA
!       WRITE(1,1111)L,IL(L),JL(L),(FCAY(L,K),K=1,KC)
!       ENDDO
!      ENDIF
!
!**********************************************************************C
!
! **  ADD VEGETATION DRAG TO HORIZONTAL ADVECTIVE ACCELERATIONS
!
!----------------------------------------------------------------------C
! 
      IF(ISVEG.GE.1)THEN
!
      DO L=2,LA
        FXVEGE(L)=0.
        FYVEGE(L)=0.
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
        LW=L-1
        LE=L+1
        LS=LSC(L)
        LN=LNC(L)
        LNW=LNWC(L)
        LSE=LSEC(L)
        VTMPATU=0.25*(V1(L,K)+V1(LW,K)+V1(LN,K)+V1(LNW,K))
        UTMPATV=0.25*(U1(L,K)+U1(LE,K)+U1(LS,K)+U1(LSE,K))
        UMAGTMP=SQRT( U1(L,K)*U1(L,K)+VTMPATU*VTMPATU )
        VMAGTMP=SQRT( UTMPATV*UTMPATV+V1(L,K)*V1(L,K) )
        FXVEG(L,K)=UMAGTMP*SUB(L)*DXYU(L)*FXVEG(L,K)
        FYVEG(L,K)=VMAGTMP*SVB(L)*DXYV(L)*FYVEG(L,K)
        FXVEGE(L)=FXVEGE(L)+FXVEG(L,K)*DZC(K)
        FYVEGE(L)=FYVEGE(L)+FYVEG(L,K)*DZC(K)
       ENDDO
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
        FXVEG(L,K)=FXVEG(L,K)*U1(L,K)
        FYVEG(L,K)=FYVEG(L,K)*V1(L,K)
        FX(L,K)=FX(L,K)+FXVEG(L,K)-FXVEGE(L)*U1(L,K)
        FY(L,K)=FY(L,K)+FYVEG(L,K)-FYVEGE(L)*V1(L,K)
       ENDDO
      ENDDO
!
      DO L=2,LA
        FXVEGE(L)=DXYIU(L)*FXVEGE(L)/H1U(L)
        FYVEGE(L)=DXYIV(L)*FYVEGE(L)/H1V(L)
      ENDDO
!
      ENDIF
!
!
!**********************************************************************C
!
! **  ADD HORIZONTAL MOMENTUN DIFFUSION TO ADVECTIVE ACCELERATIONS
!
!----------------------------------------------------------------------C
! 
      IF(ISHDMF.GE.1)THEN
!
      DO K=1,KC
       DO L=2,LA
       LS=LSC(L)
       LN=LNC(L)
       FX(L,K)=FX(L,K)-SDX(L)*(FMDUX(L,K)-FMDUX(L-1,K)+FMDUY(LN,K)-FMDUY(L,K))
       FY(L,K)=FY(L,K)-SDY(L)*(FMDVX(L+1,K)-FMDVX(L,K)+FMDVY(L,K)-FMDVY(LS,K))
       ENDDO
      ENDDO
!
      ENDIF
!
!**********************************************************************C
!
! **  ADD BODY FORCE TO ADVECTIVE ACCELERATIONS
! **  DISTRIBUTE UNIFORMLY OVER ALL LAYERS IF ISBODYF=1
! **  DISTRIBUTE OVER SURFACE LAYER IF ISBODYF=2
!
!----------------------------------------------------------------------C
!
      IF(ISBODYF.EQ.1.OR.ISUVDA.GE.1)THEN
!      
        DO K=1,KC
          DZICK=1./DZC(K)
          DO L=2,LA
            FX(L,K)=FX(L,K)-DXYU(L)*FBODYFX(L,K)
            FY(L,K)=FY(L,K)-DXYV(L)*FBODYFY(L,K)
          ENDDO
        ENDDO
!
      ENDIF
!
      IF(ISBODYF.EQ.2)THEN
!
        DZICKC=1./DZC(KC)
        DO L=2,LA
          FX(L,KC)=FX(L,KC)-DZICKC*DXYU(L)*FBODYFX(L,KC)
          FY(L,KC)=FY(L,KC)-DZICKC*DXYV(L)*FBODYFY(L,KC)
        ENDDO
!
      ENDIF
!
!**********************************************************************C
!
! ** ADD EXPLICIT NONHYDROSTATIC PRESSURE
!
      IF(KC.GT.1.AND.ISPNHYDS.GE.1) THEN
!
      TMPVAL=2./(DZC(1)+DZC(2))  
      DO L=2,LA
        DZPC(L,1)=TMPVAL*(PNHYDS(L,2)-PNHYDS(L,1))
      ENDDO
!
      TMPVAL=2./(DZC(KC)+DZC(KC-1))  
      DO L=2,LA
        DZPC(L,KC)=TMPVAL*(PNHYDS(L,KC)-PNHYDS(L,KC-1))
      ENDDO

      IF(KC.GE.3)THEN
	DO K=2,KS
	  TMPVAL=2./(DZC(K+1)+2.*DZC(K)+DZC(K-1))
        DO L=2,LA
          DZPC(L,K)=TMPVAL*(PNHYDS(L,K+1)-PNHYDS(L,K-1)) 
        ENDDO
      ENDDO
	ENDIF
!
      DO K=1,KC
        DO L=2,LA
	    LS=LSC(L)
	    DZPU=0.5*(DZPC(L,K)+DZPC(L-1,K))
	    DZPV=0.5*(DZPC(L,K)+DZPC(LS ,K))
          FX(L,K)=FX(L,K)+SUB(L)*DYU(L)*(HU(L)*(PNHYDS(L,K)-PNHYDS(L-1,K))-( BELV(L)-BELV(L-1)+  &
          ZZ(K)*(HP(L)-HP(L-1)) )*DZPU )           
          FY(L,K)=FY(L,K)+SVB(L)*DXV(L)*(HV(L)*(PNHYDS(L,K)-PNHYDS(LS,K))-( BELV(L)-BELV(LS )+   &
                  ZZ(K)*(HP(L)-HP(LS )) )*DZPV )           
        ENDDO
      ENDDO
!
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE EXTERNAL ACCELERATIONS
!
!----------------------------------------------------------------------C
! 
      DO K=1,KC
       DO L=2,LA
        FCAXE(L)=FCAXE(L)+FCAX(L,K)*DZC(K)
        FCAYE(L)=FCAYE(L)+FCAY(L,K)*DZC(K)
        FXE(L)=FXE(L)+FX(L,K)*DZC(K)
        FYE(L)=FYE(L)+FY(L,K)*DZC(K)
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
! 
! **  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.
!
      IF(ISWAVE.GE.1)THEN
!
      IF(N.LT.NTSWV)THEN
         TMPVAL=FLOAT(N)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVAL)
       ELSE
        WVFACT=1.0
      ENDIF        
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FXE(L)=FXE(L)+WVFACT*SAAX(L)*FXWAVE(L,K)*DZC(K)
         FYE(L)=FYE(L)+WVFACT*SAAY(L)*FYWAVE(L,K)*DZC(K)
        ENDDO
       ENDDO
      ENDDO
!
      ENDIF
!
!**********************************************************************C
!
! **  COMPLETE CALCULATION OF INTERNAL ADVECTIVE ACCELERATIONS
!
!----------------------------------------------------------------------C
! 
      DO K=1,KC
       DO L=2,LA
         FX(L,K)=FX(L,K)+SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
         FY(L,K)=FY(L,K)+SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)
        ENDDO
      ENDDO
!
!
!GVCDIAG
!      IF(N.LE.6)THEN
!      DO L=2,LA
!	  DO K=1,KC
!	    FXTMP(K)=DZC(K)*FX(L,K)
!	  ENDDO
!	  WRITE(8,891)N,IL(L),JL(L),FXE(L),(FXTMP(K),K=1,KC)
!      ENDDO
!      DO L=2,LA
!	  DO K=1,KC
!	    FXTMP(K)=DZC(K)*FCAX(L,K)
!	  ENDDO
!	  WRITE(8,892)N,IL(L),JL(L),FCAXE(L),(FXTMP(K),K=1,KC)
!      ENDDO
!      ENDIF
!  891 FORMAT('FINAL FX.FXE SIG ',3I5,10E14.6)
!  892 FORMAT('FINAL FC,FCE SIG ',3I5,10E14.6)
!GVCDIAG
!----------------------------------------------------------------------C
! 
! **  ADD NET WAVE REYNOLDS STRESSES TO INTERNAL ADVECTIVE ACCEL.
!
      IF(ISWAVE.GE.1)THEN
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FX(L,K)=FX(L,K)+WVFACT*SAAX(L)*FXWAVE(L,K)
         FY(L,K)=FY(L,K)+WVFACT*SAAY(L)*FYWAVE(L,K)
        ENDDO
       ENDDO
      ENDDO
!
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR
! **  THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP
! **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
!
!----------------------------------------------------------------------C
!
!      IINTPG=0
!
!     ORIGINAL
!
      IF(IINTPG.EQ.0)THEN
!
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(L-1,K+1))*DZC(K+1)      &
                +(B(L,K)-B(L-1,K))*DZC(K))-(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*                 &
                (BELV(L)-BELV(L-1)+Z(K)*(HP(L)-HP(L-1))) )
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*(HV(L)*((B(L,K+1)-B(LS,K+1))*DZC(K+1)         &
                +(B(L,K)-B(LS,K))*DZC(K))-(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*                    &
                (BELV(L)-BELV(LS)+Z(K)*(HP(L)-HP(LS))) )
      ENDDO
      ENDDO
!
	ENDIF
!
!     JACOBIAN 2 LAYERS
!
      IF(IINTPG.EQ.1.AND.KC.LE.2)THEN
!
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(L-1,K+1))*DZC(K+1)      &
                +(B(L,K)-B(L-1,K))*DZC(K))-(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*(BELV(L)-        &
                BELV(L-1)+Z(K)*(HP(L)-HP(L-1))))
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)       &
                +(B(L,K)-B(LS,K))*DZC(K))-(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*(BELV(L)-BELV(LS)+  &
                Z(K)*(HP(L)-HP(LS))) )
      ENDDO
      ENDDO
!
	ENDIF
!
!     JACOBIAN 
!
      IF(IINTPG.EQ.1.AND.KC.GT.2)THEN
!
      K=1
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*(HU(L)*( 0.25*(B(L,K+2)-B(L-1,K+2))*DZC(K+2)  &
                +0.75*(B(L,K+1)-B(L-1,K+1))*DZC(K+1)+(B(L,K)-B(L-1,K))*DZC(K))-0.25*(B(L,K+2)     &
                -B(L,K+1)+B(L-1,K+2)-B(L-1,K+1))*(BELV(L)-BELV(L-1)+Z(K+1)*(HP(L)-HP(L-1)))       &
                -0.50*(B(L,K+2)-B(L,K+1)+B(L-1,K+2)-B(L-1,K+1))*(BELV(L)-BELV(L-1)+Z(K+1)*(HP(L)  &
                -HP(L-1))) )
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*(HV(L)*(0.25*(B(L,K+2)-B(LS ,K+2))*DZC(K+2)   &
                +0.75*(B(L,K+1)-B(LS ,K+1))*DZC(K+1)+(B(L,K)-B(LS,K))*DZC(K))-0.25*(B(L,K+2)      &
                -B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*(BELV(L)-BELV(LS)+Z(K+1)*(HP(L)-HP(LS)))-0.50    &
                *(B(L,K+1)-B(L,K)+B(LS ,K+1)-B(LS,K))*(BELV(L)-BELV(LS)+Z(K)*(HP(L)-HP(LS))))
      ENDDO
!
      K=KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*(HU(L)*((B(L,K+1)-B(L-1,K+1))*DZC(K+1)        &
                +0.75*(B(L,K)-B(L-1,K))*DZC(K)+0.25*(B(L,K-1)-B(L-1,K-1))*DZC(K-1))               &
                -0.50*(B(L,K+1)-B(L,K+1)+B(L-1,K+1)-B(L-1,K+1))*(BELV(L)-BELV(L-1)+Z(K+1)*(HP(L)  &
                -HP(L-1)))-0.25*(B(L,K)-B(L,K-1)+B(L-1,K)-B(L-1,K-1))*(BELV(L)-BELV(L-1)+Z(K-1)*  &
                (HP(L)-HP(L-1))) ) 
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*(HV(L)*((B(L,K+1)-B(LS ,K+1))*DZC(K+1)        &
                +0.75*(B(L,K)-B(LS,K))*DZC(K)+0.25*(B(L,K-1)-B(LS ,K-1))*DZC(K-1))-0.50*(B(L,K+1) &
                -B(L,K)+B(LS ,K+1)-B(LS,K))*(BELV(L)-BELV(LS )+Z(K  )*(HP(L)-HP(LS)))-0.25*(B(L,K)&
                -B(L,K-1)+B(LS,K)-B(LS,K-1))*(BELV(L)-BELV(LS )+Z(K-1)*(HP(L)-HP(LS))))
      ENDDO
!
      IF(KC.GT.3)THEN
      DO K=2,KS-1
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*(HU(L)*(0.25*(B(L,K+2)-B(L-1,K+2))*DZC(K+2)   &
                +0.75*(B(L,K+1)-B(L-1,K+1))*DZC(K+1)+0.75*(B(L,K)-B(L-1,K))*DZC(K)+0.25*(B(L,K-1) &
                -B(L-1,K-1))*DZC(K-1))-0.25*(B(L,K+2)-B(L,K+1)+B(L-1,K+2)-B(L-1,K+1))*(BELV(L)    &
                -BELV(L-1)+Z(K+1)*(HP(L)-HP(L-1)))-0.50*(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*    &
                (BELV(L)-BELV(L-1)+Z(K  )*(HP(L)-HP(L-1)))-0.25*(B(L,K)-B(L,K-1)+B(L-1,K)         &
                -B(L-1,K-1))*(BELV(L)-BELV(L-1)+Z(K-1)*(HP(L)-HP(L-1))))
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*(HV(L)*(0.25*(B(L,K+2)-B(LS,K+2))*DZC(K+2)    &
                +0.75*(B(L,K+1)-B(LS ,K+1))*DZC(K+1)+0.75*(B(L,K)-B(LS,K))*DZC(K)+0.25*(B(L,K-1)  &
                -B(LS,K-1))*DZC(K-1))-0.25*(B(L,K+2)-B(L,K+1)+B(LS,K+2)-B(LS,K+1))*(BELV(L)       &
                -BELV(LS)+Z(K+1)*(HP(L)-HP(LS)))-0.50*(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*        &
                (BELV(L)-BELV(LS)+Z(K)*(HP(L)-HP(LS)))-0.25*(B(L,K)-B(L,K-1)+B(LS,K)-B(LS ,K-1))* &
                (BELV(L)-BELV(LS )+Z(K-1)*(HP(L)-HP(LS ))) )
      ENDDO
      ENDDO
      ENDIF
!
      ENDIF
!
!     FINITE VOLUME
!
      IF(IINTPG.EQ.2)THEN
!
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*((HP(L)*B(L,K+1)-HP(L-1)*B(L-1,K+1))*DZC(K+1) &
                +(HP(L)*B(L,K)-HP(L-1)*B(L-1,K))*DZC(K))-RNEW*SBX(L)*GP*(BELV(L)-BELV(L-1))*      &
                (HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(L-1)*B(L-1,K+1)-HP(L-1)*B(L-1,K))-RNEW*SBX(L)*GP  &
                *(HP(L)-HP(L-1))*(HP(L)*ZZ(K+1)*B(L,K+1)-HP(L)*ZZ(K)*B(L,K)+HP(L-1)*ZZ(K+1)       &
                *B(L-1,K+1)-HP(L-1)*ZZ(K)*B(L-1,K))
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*((HP(L)*B(L,K+1)-HP(LS)*B(LS,K+1))*DZC(K+1)   &
                +( HP(L)*B(L,K)-HP(LS)*B(LS,K))*DZC(K))-RNEW*SBY(L)*GP*(BELV(L)-BELV(LS))*        &
                ( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LS)*B(LS,K+1)-HP(LS)*B(LS,K))-RNEW*SBY(L)*GP     &
                *(HP(L)-HP(LS))*(HP(L)*ZZ(K+1)*B(L,K+1)-HP(L)*ZZ(K)*B(L,K)+HP(LS)*ZZ(K+1)*        &
                B(LS,K+1)-HP(LS)*ZZ(K)*B(LS,K) )
      ENDDO
      ENDDO
!
      ENDIF
!
!     IF(N.EQ.1)THEN
!       OPEN(1,FILE='BUOY.DIA',STATUS='UNKNOWN')
!       DO L=2,LA
!        DO K=1,KS
!        TMP3D(K)=SUBO(L)*FBBX(L,K)
!        ENDDO
!       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
!        DO K=1,KS
!        TMP3D(K)=SVBO(L)*FBBY(L,K)
!        ENDDO
!       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
!       ENDDO
!       CLOSE(1)
!     ENDIF
!
! 1111 FORMAT(2I5,2X,8E12.4)       
!
!**********************************************************************C
!
! **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
       RCDZF=CDZF(K)
       DO L=2,LA
        DU(L,K)=RCDZF*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI+DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)          &
                +FBBX(L,K)+SNLT*(FX(L,K)-FX(L,K+1))) )
        DV(L,K)=RCDZF*(H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI+DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K) &
                +SNLT*(FY(L,K)-FY(L,K+1))) )
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        DU(L,K)=SUB(L)*DU(L,K)
        DV(L,K)=SVB(L)*DV(L,K)
       ENDDO
      ENDDO
C
      IF(ISTL.EQ.2)THEN
C 
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
        DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      IF(N.LE.4)THEN
        CLOSE(1)
      ENDIF
C
 1112 FORMAT('N,NW,NS,I,J,K,NF,H,Q,QU,FUU,FVV=',/,2X,7I5,5E12.4)
C
C**********************************************************************C
C
      RETURN
      END