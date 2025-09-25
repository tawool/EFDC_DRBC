!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALEXP1D (ISTL)
!
! **  SUBROUTINE CALEXP1D CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
! **  FOR 1D CHANNEL MODE
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
!**********************************************************************C
!
! **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
!
!----------------------------------------------------------------------C
!
      FXE(1)=0.
      FYE(1)=0.
!
      FXE(LC)=0.
      FYE(LC)=0.
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
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
        GOTO 300
!       IF(ISCDMA.EQ.0) GOTO 300
!       IF(ISCDMA.EQ.1) GOTO 400
!       IF(ISCDMA.EQ.2) GOTO 350
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
!      DO K=1,KC
!       DO L=2,LA
!        UHDY2(L,K)=UHDY(L,K)
!        VHDX2(L,K)=VHDX(L,K)
!        W2(L,K)=W(L,K)
!       ENDDO
!      ENDDO
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
!
       FUHU(L,K)=MAX(UHB,0.)*U(L,K)+MIN(UHB,0.)*U(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U(LS,K)+MIN(VHC,0.)*U(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V(L-1,K)+MIN(UHC,0.)*V(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V(L,K)+MIN(VHB,0.)*V(LN,K)
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
!      DO K=1,KS
!       DO L=2,LA
!       LS=LSC(L)
!       WU=0.25*DXYU(L)*(W2(L,K)+W2(L-1,K))
!       WV=0.25*DXYV(L)*(W2(L,K)+W2(LS,K))
!       FWU(L,K)=MAX(WU,0.)*U1(L,K)
!     &        +MIN(WU,0.)*U1(L,K+1)
!       FWV(L,K)=MAX(WV,0.)*V1(L,K)
!     &        +MIN(WV,0.)*V1(L,K+1)
!       ENDDO
!      ENDDO
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
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
!      DO K=1,KS
!       DO L=2,LA
!       LS=LSC(L)
!       WU=0.5*DXYU(L)*(W(L,K)+W(L-1,K))
!       WV=0.5*DXYV(L)*(W(L,K)+W(LS,K))
!       FWU(L,K)=MAX(WU,0.)*U1(L,K)
!     &        +MIN(WU,0.)*U1(L,K+1)
!       FWV(L,K)=MAX(WV,0.)*V1(L,K)
!     &        +MIN(WV,0.)*V1(L,K+1)
!       ENDDO
!      ENDDO
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
!      DO K=1,KS
!       DO L=2,LA
!       LS=LSC(L)
!       WU=0.25*DXYU(L)*(W(L,K)+W(L-1,K))
!       WV=0.25*DXYV(L)*(W(L,K)+W(LS,K))
!       FWU(L,K)=MAX(WU,0.)*U2(L,K)
!     &        +MIN(WU,0.)*U2(L,K+1)
!       FWV(L,K)=MAX(WV,0.)*V2(L,K)
!     &        +MIN(WV,0.)*V2(L,K+1)
!       ENDDO
!      ENDDO
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
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
!      DO K=1,KS
!      DO L=2,LA
!      LS=LSC(L)    
!      FWU(L,K)=0.25*DXYU(L)*(W(L,K)+W(L-1,K))*(U(L,K+1)+U(L,K))
!      FWV(L,K)=0.25*DXYV(L)*(W(L,K)+W(LS,K))*(V(L,K+1)+V(L,K))
!      ENDDO
!      ENDDO
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
! **  CALCULATE ADVECTIVE ACCELERATIONS 
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       FX(L,K)=SAAX(L)*(FUHU(L,K)-FUHU(L-1,K)+FVHU(LN,K)-FVHU(L,K))
       FY(L,K)=SAAY(L)*(FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K))
       ENDDO
      ENDDO
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
      IF(ISWAVE.GE.2)THEN
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
!      DO K=1,KC
!       DO L=2,LA
!         FX(L,K)=FX(L,K)+SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
!         FY(L,K)=FY(L,K)+SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)
!        ENDDO
!      ENDDO
!
!----------------------------------------------------------------------C
! 
! **  ADD NET WAVE REYNOLDS STRESSES TO INTERNAL ADVECTIVE ACCEL.
!
      IF(ISWAVE.GE.2)THEN
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
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
 FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*(HU(L)*((B(L,K+1)-B(L-1,K+1))*DZC(K+1)             &
           +(B(L,K)-B(L-1,K))*DZC(K))+(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*                      &
           (HMP(L)-HMP(L-1)-Z(K)*(HP(L)-HP(L-1))) )
 FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)            &
           +(B(L,K)-B(LS,K))*DZC(K) )+(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*                        &
           (HMP(L)-HMP(LS)-Z(K)*(HP(L)-HP(LS))))  
      ENDDO
      ENDDO
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
 1111 FORMAT(2I5,2X,8E12.4)       
!
!**********************************************************************C
!
! **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
!
!----------------------------------------------------------------------C
!
!      DO K=1,KS
!       RCDZF=CDZF(K)
!       DO L=2,LA
!        DU(L,K)=RCDZF*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI
!     &           +DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)+FBBX(L,K)
!     &           +SNLT*(FX(L,K)-FX(L,K+1))) )
!        DV(L,K)=RCDZF*( H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI
!     &           +DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K)
!     &           +SNLT*(FY(L,K)-FY(L,K+1))) )
!       ENDDO
!      ENDDO
!
!      IF(ISTL.EQ.2)THEN
! 
!      DO ND=1,NDM
!       LF=2+(ND-1)*LDM
!       LL=LF+LDM-1
!       DO L=LF,LL
!        DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
!        DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
!       ENDDO
!      ENDDO
!
!      ENDIF
!
!**********************************************************************C
!
      RETURN
      END
