C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALEXP2TGVC
C
C **  SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
C **  THIS SUBROUTINE IS CURRENT PRODUCTION VERSION
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
C 05/60/2009        John Hamrick        05/60/2009        John Hamrick
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      DIMENSION TMPVEC1(KCM),TMPVEC2(KCM),FXTMP(KCM)
      DIMENSION QMCSINKX(LCM,KCM),QMCSOURX(LCM,KCM),
     &          QMCSINKY(LCM,KCM),QMCSOURY(LCM,KCM)
      DIMENSION FUHJ(LCM,KCM),FVHJ(LCM,KCM),DZPC(LCM,KCM)
C
C**********************************************************************C
C
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
        DELTD2=0.5*DT
        DELTI=1./DELT
      ELSE
        DELT=DTDYN
        DELTD2=0.5*DTDYN
        DELTI=1./DELT
      END IF 
c
      ISTL=2  
C
C**********************************************************************C
C
C      IF(N.EQ.1)THEN
C        OPEN(1,FILE='MFLUX.DIA')
C        CLOSE(1,STATUS='DELETE')
C      ENDIF
C
C      IF(N.LE.4)THEN
C        OPEN(1,FILE='MFLUX.DIA',POSITION='APPEND')
C      ENDIF
C
C**********************************************************************C
C
C **  INITIALIZE MOMENTUM FLUXES AND CORIOLIS TERMS
C **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
	DO L=1,LC
	  FUHU(L,K)=0.0
	  FUHV(L,K)=0.0
	  FVHU(L,K)=0.0
	  FVHV(L,K)=0.0
	  FUHJ(L,K)=0.0
	  FVHJ(L,K)=0.0
	  FWU(L,K)=0.0
	  FWV(L,K)=0.0
	  CAC(L,K)=0.0
	  FCAX(L,K)=0.0
	  FCAY(L,K)=0.0
	  FCAX1(L,K)=0.0
	  FCAY1(L,K)=0.0
	  FX(L,K)=0.0
	  FY(L,K)=0.0
C	  FBBX(L,K)=0.0
C	  FBBY(L,K)=0.0
	  DU(L,K)=0.0
	  DV(L,K)=0.0
      ENDDO
	ENDDO
C
      IF(IS2TLPG.EQ.0)THEN
      DO K=1,KC
	DO L=1,LC
	  FBBX(L,K)=0.0
	  FBBY(L,K)=0.0
      ENDDO
	ENDDO
      ENDIF
C
      DO L=1,LC
        FCAXE(L)=0.
        FCAYE(L)=0.
        FCAX1E(L)=0.
        FCAY1E(L)=0.
        FXE(L)=0.
        FYE(L)=0.
      ENDDO
C
C**********************************************************************C
C
C **  TWO TIME LEVEL STEP
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N 
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
        DO L=2,LA
          UHDY2(L,K)=UHDY(L,K)
          VHDX2(L,K)=VHDX(L,K)
          W2(L,K)=2.*W(L,K)
        ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C MODIFIED FOR GVC - JMH 08/04/04
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.5*(GVCSCLU(L)*UHDY2(L,K)+GVCSCLU(LS )*UHDY2(LS,K))
       UHB=0.5*(GVCSCLU(L)*UHDY2(L,K)+GVCSCLU(L+1)*UHDY2(L+1,K))
       VHC=0.5*(GVCSCLV(L)*VHDX2(L,K)+GVCSCLV(L-1)*VHDX2(L-1,K))
       VHB=0.5*(GVCSCLV(L)*VHDX2(L,K)+GVCSCLV(LN )*VHDX2(LN,K))
C
       FUHU(L,K)=MAX(UHB,0.)*U(L,K)
     &         +MIN(UHB,0.)*U(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U(LS,K)
     &         +MIN(VHC,0.)*U(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V(L-1,K)
     &         +MIN(UHC,0.)*V(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V(L,K)
     &         +MIN(VHB,0.)*V(LN,K)
       FUHJ(L,K)=0.
       FVHJ(L,K)=0.
       ENDDO
      ENDDO
C
C ADD RETURN FLOW MOMENTUM FLUX
C
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
         TMPANG=0.017453*ANGWRMFD(NWR)
         TMPANG=COS(TMPANG)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=TMPANG*QMF*QMF/(H1P(LD)*DZC(KD)*DZC(KD)*BQWRMFD(NWR))
         IF(NQWRMFD(NWR).EQ.1)  FUHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.2)  FVHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.3)  FUHJ(LD+1   ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.4)  FVHJ(LNC(LD),KD)=QUMF
         IF(NQWRMFD(NWR).EQ.-1) FUHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-2) FVHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-3) FUHJ(LD+1   ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-4) FVHJ(LNC(LD),KD)=-QUMF
C         IF(N.LE.4)THEN
C           WRITE(1,1112)N,NWR,NS,ID,JD,KD,NQWRMFD(NWR),H1P(LD),QMF,
C     &                  QUMF,FUHJ(LD,KD),FVHJ(LD,KD)
C         ENDIF
       ENDIF
      ENDDO
C
C ** HARDWIRE FOR PEACH BOTTOM
C
C      DO K=1,KC
C       FVHV(535,K)=700./H1P(535)
C      ENDDO
C
C ** END HARDWIRE FOR PEACH BOTTOM
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.25*DXYU(L)*(W2(L,K)+W2(L-1,K))
       WV=0.25*DXYV(L)*(W2(L,K)+W2(LS,K))
       FWU(L,K)=MAX(WU,0.)*U(L,K)
     &        +MIN(WU,0.)*U(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V(L,K)
     &        +MIN(WV,0.)*V(L,K+1)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C ** BLOCK MOMENTUM FLUX ON LAND SIDE OF TRIANGULAR CELLS
C
      DO K=1,KC
        DO L=1,LA
          FUHU(L,K)=STCUV(L)*FUHU(L,K)
          FVHV(L,K)=STCUV(L)*FVHV(L,K)
        ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS
C
C----------------------------------------------------------------------C
C
C MODIFIED FOR GVC - JMH 08/04/04
C
      DO K=1,KC
        DO L=1,LC
          CAC(L,K)=0.0
        ENDDO
      ENDDO
C
      IF(ISDCCA.EQ.0)THEN
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)          
C       CAC(L,K)=( FCORC(L)*DXYP(L)
C     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
C     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       CAC(L,K)=( FCORC(L)*DXYP(L)
     &        +0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L)
     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*DXDJ(L) )*HP(L)*GVCSCLP(L)
       ENDDO
      ENDDO
C
      ELSE
C
      CFMAX=CF
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)          
C       CAC(L,K)=( FCORC(L)*DXYP(L)
C     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
C     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       CAC(L,K)=( FCORC(L)*DXYP(L)
     &        +0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L)
     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*DXDJ(L) )*HP(L)*GVCSCLP(L)
       CFEFF=ABS(CAC(L,K))*DXYIP(L)*HPI(L)
       CFMAX=MAX(CFMAX,CFEFF)
       ENDDO
      ENDDO
C
      ENDIF
C
      IF(ISDCCA.EQ.1)THEN
      IF(N.EQ.NTS)THEN
       OPEN(1,FILE='CORC1.DIA')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='CORC1.DIA')
       K=1
       DO L=2,LA
       LN=LNC(L)          
       WRITE(1,1111)IL(L),JL(L),LN,V(LN,K),V(L,K),DYU(L+1),DYU(L),
     &        U(L+1,K),U(L,K),DXV(LN),DXV(L),HP(L),CAC(L,K)
       ENDDO
       CLOSE(1)
      ENDIF
      ENDIF
C
 1111 FORMAT(3I5,10E13.4)
 1113 FORMAT(2I5,10E13.4)
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS 
C
C----------------------------------------------------------------------C
C
C **  STANDARD CALCULATION
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX(L,K)=0.25*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
     &                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
       FCAY(L,K)=0.25*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
     &                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
       FX(L,K)=SAAX(L)*(FUHU(L,K)-FUHU(L-1,K)+FVHU(LN,K)-FVHU(L,K)
     &                 +FUHJ(L,K) )
       FY(L,K)=SAAY(L)*(FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K)
     &                 +FVHJ(L,K) )
       ENDDO
      ENDDO
C
C MODIFIED FOR GVC - JMH 08/04/04
C
      DO K=1,KC
       DO L=1,LA
       FX(L,K)=SUB3D(L,K)*FX(L,K)
       FY(L,K)=SVB3D(L,K)*FY(L,K)
       FCAX(L,K)=SUB3D(L,K)*FCAX(L,K)
       FCAY(L,K)=SVB3D(L,K)*FCAY(L,K)
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CORIOLIS-CURVATURE DIAGNOSTICS
C
      IF(ISDCCA.EQ.1)THEN
      IF(N.EQ.NTS)THEN
       OPEN(1,FILE='CORC2.DIA')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='CORC2.DIA')
       K=1
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       WRITE(1,1113)IL(L),JL(L),CAC(L,K),V(LN,K),V(L,K),
     &                              CAC(L-1,K),V(LNW,K),V(L-1,K)
       ENDDO
       CLOSE(1)
      ENDIF
      ENDIF
C
      IF(ISDCCA.EQ.1)THEN
      IF(N.EQ.NTS)THEN
       OPEN(1,FILE='CORC3.DIA')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='CORC3.DIA')
       K=1
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       WRITE(1,1113)IL(L),JL(L),CAC(L,K),U(L+1,K),U(L,K),
     &                              CAC(LS,K),U(LSE,K),U(LS,K)
       ENDDO
       CLOSE(1)
      ENDIF
      ENDIF
C
      IF(ISDCCA.EQ.1)THEN
      IF(N.EQ.NTS)THEN
       OPEN(1,FILE='CORC4.DIA')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='CORC4.DIA')
       DO L=2,LA
       WRITE(1,1113)IL(L),JL(L),(FCAX(L,K),K=1,KC)
       ENDDO
       DO L=2,LA
       WRITE(1,1113)IL(L),JL(L),(FCAY(L,K),K=1,KC)
       ENDDO
       CLOSE(1)
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  ADD VEGETATION DRAG TO HORIZONTAL ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      IF(ISVEG.GE.1)THEN
C
      DO L=2,LA
        FXVEGE(L)=0.
        FYVEGE(L)=0.
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        LW=L-1
        LE=L+1
        LS=LSC(L)
        LN=LNC(L)
        LNW=LNWC(L)
        LSE=LSEC(L)
        VTMPATU=0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))
        UTMPATV=0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))
        UMAGTMP=SQRT( U(L,K)*U(L,K)+VTMPATU*VTMPATU )
        VMAGTMP=SQRT( UTMPATV*UTMPATV+V1(L,K)*V(L,K) )
        FXVEG(L,K)=UMAGTMP*SUB3D(L,K)*DXYU(L)*FXVEG(L,K)
        FYVEG(L,K)=VMAGTMP*SVB3D(L,K)*DXYV(L)*FYVEG(L,K)
        FXVEGE(L)=FXVEGE(L)+FXVEG(L,K)*DZC(K)
        FYVEGE(L)=FYVEGE(L)+FYVEG(L,K)*DZC(K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        FXVEG(L,K)=FXVEG(L,K)*U(L,K)
        FYVEG(L,K)=FYVEG(L,K)*V(L,K)
        FX(L,K)=FX(L,K)+FXVEG(L,K)-FXVEGE(L)*U(L,K)
        FY(L,K)=FY(L,K)+FYVEG(L,K)-FYVEGE(L)*V(L,K)
       ENDDO
      ENDDO
C
      DO L=2,LA
        FXVEGE(L)=DXYIU(L)*FXVEGE(L)/HU(L)
        FYVEGE(L)=DXYIV(L)*FYVEGE(L)/HV(L)
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  ADD HORIZONTAL MOMENTUN DIFFUSION TO ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      IF(ISHDMF.GE.1)THEN
C
c  below loop from calexpgvc.for
c
C      DO K=1,KC
C       DO L=2,LA
C       LS=LSC(L)
C       LN=LNC(L)
C       FX(L,K)=FX(L,K)-SUB3D(L,K)*SDX(L)*
C     &         (FMDUX(L,K)-FMDUX(L-1,K)+FMDUY(LN,K)-FMDUY(L,K))
C       FY(L,K)=FY(L,K)-SVB3D(L,K)*SDY(L)*
C     &         (FMDVX(L+1,K)-FMDVX(L,K)+FMDVY(L,K)-FMDVY(LS,K))
C       ENDDO
C      ENDDO
C
C  below loop from calexp2t.for
C
C      DO K=1,KC
C       DO L=2,LA
C         IF(LMASKDRY(L))THEN
C          FX(L,K)=FX(L,K)
C     &        -SDX(L)*(FMDUX(L,K)+FMDUY(L,K))
C          FY(L,K)=FY(L,K)
C     &        -SDY(L)*(FMDVX(L,K)+FMDVY(L,K))
C         ENDIF
C       ENDDO
C      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  ADD BODY FORCE TO ADVECTIVE ACCELERATIONS
C **  DISTRIBUTE UNIFORMLY OVER ALL LAYERS IF ISBODYF=1
C **  DISTRIBUTE OVER SURFACE LAYER IF ISBODYF=2
C
C----------------------------------------------------------------------C
C
      IF(ISBODYF.EQ.1.OR.ISUVDA.GE.1)THEN
C      
        DO K=1,KC
          DZICK=1./DZC(K)
          DO L=2,LA
            FX(L,K)=FX(L,K)-GVCSCLU(L)*DXYU(L)*FBODYFX(L,K)
            FY(L,K)=FY(L,K)-GVCSCLV(L)*DXYV(L)*FBODYFY(L,K)
          ENDDO
        ENDDO
C
      ENDIF
C
      IF(ISBODYF.EQ.2)THEN
C
        DZICKC=1./DZC(KC)
        DO L=2,LA
          FX(L,KC)=FX(L,KC)-DZICKC*GVCSCLU(L)*DXYU(L)*FBODYFX(L,KC)
          FY(L,KC)=FY(L,KC)-DZICKC*GVCSCLV(L)*DXYV(L)*FBODYFY(L,KC)
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C ** ADD EXPLICIT NONHYDROSTATIC PRESSURE
C
      IF(KC.GT.1.AND.ISPNHYDS.GE.1) THEN
C
      TMPVAL=2./(DZC(1)+DZC(2))  
      DO L=2,LA
        DZPC(L,1)=TMPVAL*(PNHYDS(L,2)-PNHYDS(L,1))
      ENDDO
C
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
C
      DO K=1,KC
        DO L=2,LA
	    LS=LSC(L)
	    DZPU=0.5*(DZPC(L,K)+DZPC(L-1,K))
	    DZPV=0.5*(DZPC(L,K)+DZPC(LS ,K))
          FX(L,K)=FX(L,K)+SUB3D(L,K)*DYU(L)*
     &           ( HU(L)*(PNHYDS(L,K)-PNHYDS(L-1,K))
     &           -( BELV(L)-BELV(L-1)+ZZ(K)*(HP(L)-HP(L-1)) )*DZPU )           
          FY(L,K)=FY(L,K)+SVB3D(L,K)*DXV(L)*
     &           ( HV(L)*(PNHYDS(L,K)-PNHYDS(LS ,K))
     &           -( BELV(L)-BELV(LS )+ZZ(K)*(HP(L)-HP(LS )) )*DZPV )           
        ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      DO K=1,KC
       DO L=2,LA
        FCAXE(L)=FCAXE(L)+FCAX(L,K)*DZC(K)
        FCAYE(L)=FCAYE(L)+FCAY(L,K)*DZC(K)
        FXE(L)=FXE(L)+FX(L,K)*DZC(K)
        FYE(L)=FYE(L)+FY(L,K)*DZC(K)
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.
C
      IF(ISWAVE.GE.1)THEN
C
      IF(N.LT.NTSWV)THEN
         TMPVAL=FLOAT(N)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVAL)
       ELSE
        WVFACT=1.0
      ENDIF        
C
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
C
      ENDIF
C
C**********************************************************************C
C
C **  COMPLETE CALCULATION OF INTERNAL ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C
C MODIFIED FOR GVC - JMH 08/04/04
C
      DO K=1,KC
       DO L=2,LA
         FX(L,K)=FX(L,K)
     &          +SAAX(L)*SUB3D(L,K)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
         FY(L,K)=FY(L,K)
     &          +SAAY(L)*SVB3D(L,K)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)
        ENDDO
      ENDDO
C
C
CGVCDIAG
c      IF(N.LE.6)THEN
c      DO L=2,LA
c	  DO K=1,KC
c	    FXTMP(K)=DZC(K)*FX(L,K)
c	  ENDDO
c	  WRITE(8,891)N,IL(L),JL(L),FXE(L),(FXTMP(K),K=1,KC)
c      ENDDO
c      DO L=2,LA
c	  DO K=1,KC
c	    FXTMP(K)=DZC(K)*FCAX(L,K)
c	  ENDDO
c	  WRITE(8,892)N,IL(L),JL(L),FCAXE(L),(FXTMP(K),K=1,KC)
c      ENDDO
c      ENDIF
c  891 FORMAT('FINAL FX.FXE GVC ',3I5,10E14.6)
c  892 FORMAT('FINAL FC,FCE GVC ',3I5,10E14.6)
CGVCDIAG
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO INTERNAL ADVECTIVE ACCEL.
C
      IF(ISWAVE.GE.1)THEN
C
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
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR
C **  THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP
C **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
C
C----------------------------------------------------------------------C
C
c      IINTPG=0
C
C     ORIGINAL
C
      IF(IS2TLPG.EQ.1)THEN
	  ROLDPG=-0.5
	  RNEWPG=1.5
	ELSE
	  ROLDPG=0.0
	  RNEWPG=1.0
	ENDIF
C
      IF(IINTPG.EQ.0)THEN
C
C MODIFIED FOR GVC - JMH 08/04/04
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLDPG*FBBX(L,K)+RNEWPG*SBX(L)*GP*HU(L)*GVCSCLU(L)
     &           *( HU(L)*GVCSCLU(L)*( 
     &                    SUB3D(L,K+1)*(B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     &                   +SUB3D(L,K  )*(B(L,K  )-B(L-1,K  ))*DZC(K  ) )
     &           -(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*
     &            (BELV(L)-BELV(L-1)+HP(L)-HP(L-1)-(1.-Z(K))
     &           *(GVCSCLP(L)*HP(L)-GVCSCLP(L-1)*HP(L-1))) )
      FBBY(L,K)=ROLDPG*FBBY(L,K)+RNEWPG*SBY(L)*GP*HV(L)*GVCSCLV(L)
     &           *( HV(L)*GVCSCLV(L)*( 
     &                    SVB3D(L,K+1)*(B(L,K+1)-B(LS,K+1))*DZC(K+1)
     &                   +SVB3D(L,K  )*(B(L,K  )-B(LS,K  ))*DZC(K  ) )
     &           -(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*
     &            (BELV(L)-BELV(LS)+HP(L)-HP(LS)-(1.-Z(K))
     &           *(GVCSCLP(L)*HP(L)-GVCSCLP(LS)*HP(LS))) )
      ENDDO
      ENDDO
C
	ENDIF
C
C     JACOBIAN 2 LAYERS
C
      IF(IINTPG.EQ.1.AND.KC.LE.2)THEN
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=SBX(L)*GP*HU(L)*
     &            ( HU(L)*( (B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     &                     +(B(L,K)-B(L-1,K))*DZC(K) )
     &           -(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*
     &            (BELV(L)-BELV(L-1)+Z(K)*(HP(L)-HP(L-1))) )
      FBBY(L,K)=SBY(L)*GP*HV(L)*
     &            ( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)
     &                     +(B(L,K)-B(LS,K))*DZC(K) )
     &           -(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*
     &            (BELV(L)-BELV(LS)+Z(K)*(HP(L)-HP(LS))) )
      ENDDO
      ENDDO
C
	ENDIF
C
C     JACOBIAN 
C
      IF(IINTPG.EQ.1.AND.KC.GT.2)THEN
C
      K=1
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=SBX(L)*GP*HU(L)*
     &            ( HU(L)*( 0.25*(B(L,K+2)-B(L-1,K+2))*DZC(K+2)
     &                     +0.75*(B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     &                     +     (B(L,K  )-B(L-1,K  ))*DZC(K  ) )
     &           -0.25*(B(L,K+2)-B(L,K+1)+B(L-1,K+2)-B(L-1,K+1))*
     &            (BELV(L)-BELV(L-1)+Z(K+1)*(HP(L)-HP(L-1)))
     &           -0.50*(B(L,K+2)-B(L,K+1)+B(L-1,K+2)-B(L-1,K+1))*
     &            (BELV(L)-BELV(L-1)+Z(K+1)*(HP(L)-HP(L-1))) )
      FBBY(L,K)=SBY(L)*GP*HV(L)*
     &            ( HV(L)*( 0.25*(B(L,K+2)-B(LS ,K+2))*DZC(K+2)
     &                     +0.75*(B(L,K+1)-B(LS ,K+1))*DZC(K+1) 
     &                     +     (B(L,K  )-B(LS ,K  ))*DZC(K  ) )
     &           -0.25*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*
     &            (BELV(L)-BELV(LS)+Z(K+1)*(HP(L)-HP(LS)))
     &           -0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))*
     &            (BELV(L)-BELV(LS)+Z(K  )*(HP(L)-HP(LS))) )
      ENDDO
C
      K=KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=SBX(L)*GP*HU(L)*
     &            ( HU(L)*(      (B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     &                     +0.75*(B(L,K  )-B(L-1,K  ))*DZC(K  )
     &                     +0.25*(B(L,K-1)-B(L-1,K-1))*DZC(K-1) )
     &            -0.50*(B(L,K+1)-B(L,K+1)+B(L-1,K+1)-B(L-1,K+1))*
     &            (BELV(L)-BELV(L-1)+Z(K+1)*(HP(L)-HP(L-1)))
     &            -0.25*(B(L,K  )-B(L,K-1)+B(L-1,K  )-B(L-1,K-1))*
     &            (BELV(L)-BELV(L-1)+Z(K-1)*(HP(L)-HP(L-1))) ) 
      FBBY(L,K)=SBY(L)*GP*HV(L)*
     &            ( HV(L)*(      (B(L,K+1)-B(LS ,K+1))*DZC(K+1)
     &                     +0.75*(B(L,K  )-B(LS ,K  ))*DZC(K  ) 
     &                     +0.25*(B(L,K-1)-B(LS ,K-1))*DZC(K-1) )
     &            -0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))*
     &            (BELV(L)-BELV(LS )+Z(K  )*(HP(L)-HP(LS)))
     &            -0.25*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))*
     &            (BELV(L)-BELV(LS )+Z(K-1)*(HP(L)-HP(LS))) )
      ENDDO
C
      IF(KC.GT.3)THEN
      DO K=2,KS-1
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=SBX(L)*GP*HU(L)*
     &            ( HU(L)*( 0.25*(B(L,K+2)-B(L-1,K+2))*DZC(K+2)
     &                     +0.75*(B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     &                     +0.75*(B(L,K  )-B(L-1,K  ))*DZC(K  )
     &                     +0.25*(B(L,K-1)-B(L-1,K-1))*DZC(K-1) )
     &           -0.25*(B(L,K+2)-B(L,K+1)+B(L-1,K+2)-B(L-1,K+1))*
     &            (BELV(L)-BELV(L-1)+Z(K+1)*(HP(L)-HP(L-1)))
     &           -0.50*(B(L,K+1)-B(L,K  )+B(L-1,K+1)-B(L-1,K  ))*
     &            (BELV(L)-BELV(L-1)+Z(K  )*(HP(L)-HP(L-1)))
     &           -0.25*(B(L,K  )-B(L,K-1)+B(L-1,K  )-B(L-1,K-1))*
     &            (BELV(L)-BELV(L-1)+Z(K-1)*(HP(L)-HP(L-1))) )
      FBBY(L,K)=SBY(L)*GP*HV(L)*
     &            ( HV(L)*( 0.25*(B(L,K+2)-B(LS ,K+2))*DZC(K+2)
     &                     +0.75*(B(L,K+1)-B(LS ,K+1))*DZC(K+1) 
     &                     +0.75*(B(L,K  )-B(LS ,K  ))*DZC(K  ) 
     &                     +0.25*(B(L,K-1)-B(LS ,K-1))*DZC(K-1) )
     &           -0.25*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*
     &            (BELV(L)-BELV(LS)+Z(K+1)*(HP(L)-HP(LS)))
     &           -0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))*
     &            (BELV(L)-BELV(LS)+Z(K  )*(HP(L)-HP(LS)))
     &           -0.25*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))*
     &            (BELV(L)-BELV(LS )+Z(K-1)*(HP(L)-HP(LS ))) )
      ENDDO
      ENDDO
      ENDIF
C
      ENDIF
C
C     FINITE VOLUME
C
      IF(IINTPG.EQ.2)THEN
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=SBX(L)*GP*HU(L)*
     &            ( ( HP(L)*B(L,K+1)-HP(L-1)*B(L-1,K+1) )*DZC(K+1)
     &             +( HP(L)*B(L,K  )-HP(L-1)*B(L-1,K  ) )*DZC(K  ) )
     &         -RNEW*SBX(L)*GP*(BELV(L)-BELV(L-1))*
     &          ( HP(L)*B(L,K+1)-HP(L)*B(L,K)
     &           +HP(L-1)*B(L-1,K+1)-HP(L-1)*B(L-1,K) )
     &         -RNEW*SBX(L)*GP*(HP(L)-HP(L-1))*
     &          ( HP(L)*ZZ(K+1)*B(L,K+1)-HP(L)*ZZ(K)*B(L,K)
     &           +HP(L-1)*ZZ(K+1)*B(L-1,K+1)-HP(L-1)*ZZ(K)*B(L-1,K) )
      FBBY(L,K)=SBY(L)*GP*HV(L)*
     &            ( ( HP(L)*B(L,K+1)-HP(LS )*B(LS ,K+1) )*DZC(K+1)
     &             +( HP(L)*B(L,K  )-HP(LS )*B(LS ,K  ) )*DZC(K  ) )
     &         -RNEW*SBY(L)*GP*(BELV(L)-BELV(LS ))*
     &          ( HP(L)*B(L,K+1)-HP(L)*B(L,K)
     &           +HP(LS)*B(LS ,K+1)-HP(LS)*B(LS ,K) )
     &         -RNEW*SBY(L)*GP*(HP(L)-HP(LS ))*
     &          ( HP(L)*ZZ(K+1)*B(L,K+1)-HP(L)*ZZ(K)*B(L,K) 
     &           +HP(LS)*ZZ(K+1)*B(LS ,K+1)-HP(LS)*ZZ(K)*B(LS ,K) )
      ENDDO
      ENDDO
C
      ENDIF
C
C     IF(N.EQ.1)THEN
C       OPEN(1,FILE='BUOY.DIA',STATUS='UNKNOWN')
C       DO L=2,LA
C        DO K=1,KS
C        TMP3D(K)=SUBO(L)*FBBX(L,K)
C        ENDDO
C       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
C        DO K=1,KS
C        TMP3D(K)=SVBO(L)*FBBY(L,K)
C        ENDDO
C       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
C       ENDDO
C       CLOSE(1)
C     ENDIF
C
C 1111 FORMAT(2I5,2X,8E12.4)       
C
C
C MODIFIED FOR GVC - JMH 08/04/04   COMMENTED OUT JMH 02/07/06
C
C      DO K=1,KS
C       DO L=1,LA
C       FBBX(L,K)=SUB3D(L,K)*FBBX(L,K)
C       FBBY(L,K)=SVB3D(L,K)*FBBY(L,K)
C       ENDDO
C      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
C
C----------------------------------------------------------------------C
C
C MODIFIED FOR GVC - JMH 08/04/04
C
C      CDZF(K)=DZC(K)*DZC(K+1)/(DZC(K)+DZC(K+1))
C
      DO K=1,KS
       RCDZF=CDZF(K)
       DO L=2,LA
        DU(L,K)=RCDZF*( GVCSCLU(L)*HU(L)*(U(L,K+1)-U(L,K))*DELTI
     &           +DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)+FBBX(L,K)
     &           +SNLT*(FX(L,K)-FX(L,K+1))) )
        DV(L,K)=RCDZF*( GVCSCLV(L)*HV(L)*(V(L,K+1)-V(L,K))*DELTI
     &           +DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K)
     &           +SNLT*(FY(L,K)-FY(L,K+1))) )
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
C        DU(L,K)=SUB3D(L,K)*DU(L,K)
C        DV(L,K)=SVB3D(L,K)*DV(L,K)
        DU(L,K)=SUB(L)*DU(L,K)
        DV(L,K)=SVB(L)*DV(L,K)
       ENDDO
      ENDDO
C
c      IF(ISTL.EQ.2)THEN
C
C **  CHECK THE PURPOSE OF THIS LOOP WHICH ADD SURFACE STRESS TO 
C **  TO KS INTERFACE DU AND DV EQUATIONS FOR TWO TIME LEVEL 
C **  CORRECTION STEP - JMH 08/04/04
C 
C      DO L=2,LA
C        DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
C        DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
C      ENDDO
C
c      ENDIF
C
C**********************************************************************C
C
C      IF(N.LE.4)THEN
C        CLOSE(1)
C      ENDIF
C
 1112 FORMAT('N,NW,NS,I,J,K,NF,H,Q,QU,FUU,FVV=',/,2X,7I5,5E12.4)
C
C**********************************************************************C
C
      RETURN
      END
