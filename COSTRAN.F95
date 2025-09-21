!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE COSTRAN (ISTL,IS2TL,MVAR,M,CON,CON1)
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 03/05/2002        john hamrick       03/05/2002       john hamrick
!  added dynamic time stepping
! 02/07/2020        DRBC LZ - CORRECT A TYPO IN A IF-BLOCK OF KC=1. LINE 277
!
!----------------------------------------------------------------------C
!
! **  SUBROUTINE COSTRAN CALCULATES THE ADVECTIVE
! **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO
! **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES
! **  THE NUMBER OF TIME LEVELS IN THE STEP
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION CON(LCM,KCM), CON1(LCM,KCM)
      DIMENSION CONTMX(LCM,KCM),CONTMN(LCM,KCM)
      DIMENSION DELCX(LCM,KCM),DELCY(LCM,KCM),DELCZ(LCM,KCM),CONCX(LCM,KCM),CONCY(LCM,KCM),       &
                CONCZ(LCM,KCM),CONCXY(LCM,KCM),CONCXZ(LCM,KCM),CONCYZ(LCM,KCM),CONCYX(LCM,KCM),   &
                CONCZX(LCM,KCM),CONCZY(LCM,KCM),CONCXYZ(LCM,KCM),CONCYZX(LCM,KCM),CONCZXY(LCM,KCM)
      DIMENSION FQCPAD(LCM,KCM),QSUMPAD(LCM,KCM),QSUMNAD(LCM,KCM)
!     DIMENSION  CH(LCM,KCM), CMAX(LCM,KCM), CMIN(LCM,KCM),
!    &           CONT(LCM,KCM), CON2(LCM,KCM),
!
!**********************************************************************C
!
      BSMALL=1.0E-6
!BUG    WFQC=0.0
!
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT2
        DELTA=DT2
        IF(ISCDCA(MVAR).EQ.2) DELTA=DT
          DELTD2=DT
          S3TL=1.0
          S2TL=0.0
          ISUD=1
        IF(ISTL.NE.3)THEN
          DELT=DT
          DELTA=DT
          DELTD2=0.5*DT
          S3TL=0.0
          S2TL=1.0
          ISUD=0
        ENDIF
      ELSE
        DELT=DTDYN
        DELTA=DTDYN
        DELTD2=0.5*DTDYN
        S3TL=0.0
        S2TL=1.0
        ISUD=0
      END IF
!
      DELTA4=0.25*DELTA
!
      IF(ISLTMT.GE.1) ISUD=1
!
      DO K=1,KC
       DO L=1,LC
        FUHU(L,K)=0.
        FUHV(L,K)=0.
        FVHU(L,K)=0.
        FVHV(L,K)=0.
        UUU(L,K)=0.
        VVV(L,K)=0.
        DU(L,K)=0.
        DV(L,K)=0.
       ENDDO
      ENDDO
!
      IF(IS2TL.EQ.1)THEN
        ISUD=1
        DO K=1,KC
        DO L=1,LC
          CON1(L,K)=CON(L,K)
        ENDDO
        ENDDO
      ENDIF
!
!**********************************************************************C
!
        DO K=1,KC
         DO L=1,LC
         CONT(L,K)=0.
         CMAX(L,K)=0.
         CMIN(L,K)=0.
         ENDDO
        ENDDO
!
!**********************************************************************C
!
! **  CALCULATED EXTERNAL SOURCES AND SINKS
!
!----------------------------------------------------------------------C
!
      CALL CALFQC (ISTL,IS2TL,MVAR,M,CON,CON1,FQCPAD,QSUMPAD,QSUMNAD)
!     IF(ISTRAN(M).EQ.1) CALL CALFQC (ISTL,M,CON,CON1)
!     IF(ISTRAN(M).EQ.3) CALL CALFQC (ISTL,M,CON,CON1)
!     IF(M.EQ.4)THEN
!      IF(ISTOPT(4).EQ.2) CALL CALSED2(ISTL,0.5,CON1)
!      NTSWVD=ISWVSD
!      IF(ISTOPT(4).EQ.13) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.14.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.15) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.16.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
!     ENDIF
!
!**********************************************************************C
!
! **  BEGIN COMBINED ADVECTION SCHEME
!
!**********************************************************************C
!
! **  INTERMEDIATE ADVECTION CALCULATIONS
! **  CX,CY,CZ
!
      DO K=1,KC
      DO L=2,LA
       DELCX(L,K)=CON1(L,K)-CON1(L-1   ,K)
       DELCY(L,K)=CON1(L,K)-CON1(LSC(L),K)
      ENDDO
      ENDDO
!
      IF(KC.GE.2)THEN
        DO K=1,KS
        DO L=2,LA
         DELCZ(L,K)=CON1(L,K+1)-CON1(L,K)
        ENDDO
        ENDDO
      ENDIF
!
      DO K=1,KC
      DO L=2,LA
      CONCX(L,K)=CON1(L,K)+RCOSMICX(L,K)*(COSMICXP(L,K)*DELCX(L,K)+COSMICXN(L+1,K)*DELCX(L+1,K))
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      CONCY(L,K)=CON1(L,K)+RCOSMICY(L,K)*(COSMICYP(L,K)*DELCY(L,K)+COSMICYN(LNC(L),K)             &
                 *DELCY(LNC(L),K))
      ENDDO
      ENDDO
!
      IF(KC.EQ.1)THEN
        DO L=2,LA
         CONCZ(L,1)=CON1(L,1)
        ENDDO
      ENDIF
!
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCZ(L,1)=CON1(L,1)+RCOSMICZ(L,1)*(COSMICZN(L,1)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCZ(L,KC)=CON1(L,KC)+RCOSMICZ(L,KC)*(COSMICZP(L,KS)*DELCZ(L,KS) )
        ENDDO
      ENDIF
!
      IF(KC.GE.3)THEN
        DO K=2,KS
        DO L=2,LA
         CONCZ(L,K)=CON1(L,K)+RCOSMICZ(L,K)*(COSMICZP(L,K-1)*DELCZ(L,K-1)+COSMICZN(L,K)*DELCZ(L,K))
        ENDDO
        ENDDO
      ENDIF
!
!**********************************************************************C
!
! **  INTERMEDIATE ADVECTION CALCULATIONS
! **  CXY,CXZ
!
      DO K=1,KC
      DO L=2,LA
       DELCY(L,K)=CONCX(L,K)-CONCX(LSC(L),K)
      ENDDO
      ENDDO
!
      IF(KC.GE.2)THEN
        DO K=1,KS
        DO L=2,LA
         DELCZ(L,K)=CONCX(L,K+1)-CONCX(L,K)
        ENDDO
        ENDDO
      ENDIF
!
      DO K=1,KC
      DO L=2,LA
      CONCXY(L,K)=CONCX(L,K)+RCOSMICY(L,K)*(COSMICYP(L,K)*DELCY(L,K)+COSMICYN(LNC(L),K)           &
                  *DELCY(LNC(L),K) )
      ENDDO
      ENDDO
!
      IF(KC.EQ.1)THEN
        DO L=2,LA
         CONCXZ(L,1)=CONCX(L,1)
        ENDDO
      ENDIF
!
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCXZ(L,1)=CONCX(L,1)+RCOSMICZ(L,1)*(COSMICZN(L,1)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCXZ(L,KC)=CONCX(L,KC)+RCOSMICZ(L,KC)*(COSMICZP(L,KS)*DELCZ(L,KS) )
        ENDDO
      ENDIF
!
      IF(KC.GE.3)THEN
        DO K=2,KS
        DO L=2,LA
         CONCXZ(L,K)=CONCX(L,K)+RCOSMICZ(L,K)*(COSMICZP(L,K-1)*DELCZ(L,K-1)+COSMICZN(L,K)*DELCZ(L,K))
        ENDDO
        ENDDO
      ENDIF
!
!
!**********************************************************************C
!
! **  INTERMEDIATE ADVECTION CALCULATIONS
! **  CYZ,CYX
!
      DO K=1,KC
      DO L=2,LA
       DELCX(L,K)=CONCY(L,K)-CONCY(L-1,K)
      ENDDO
      ENDDO
!
      IF(KC.GE.2)THEN
        DO K=1,KS
        DO L=2,LA
         DELCZ(L,K)=CONCY(L,K+1)-CONCY(L,K)
        ENDDO
        ENDDO
      ENDIF
!
      DO K=1,KC
      DO L=2,LA
      CONCYX(L,K)=CONCY(L,K)+RCOSMICX(L,K)*(COSMICXP(L,K)*DELCX(L,K)+COSMICXN(L+1,K)*DELCX(L+1,K))
      ENDDO
      ENDDO
!
      IF(KC.EQ.1)THEN
        DO L=2,LA
!LZ      CONCZ(L,1)=CONCY(L,1)
         CONCYZ(L,1)=CONCY(L,1)                      !DRBC, LZ, 2/7/2020
        ENDDO
      ENDIF
!
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCYZ(L,1)=CONCY(L,1)+RCOSMICZ(L,1)*(COSMICZN(L,1)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCYZ(L,KC)=CONCY(L,KC)+RCOSMICZ(L,KC)*(COSMICZP(L,KS)*DELCZ(L,KS) )
        ENDDO
      ENDIF
!
      IF(KC.GE.3)THEN
        DO K=2,KS
        DO L=2,LA
         CONCYZ(L,K)=CONCY(L,K)+RCOSMICZ(L,K)*(COSMICZP(L,K-1)*DELCZ(L,K-1)+COSMICZN(L,K)*DELCZ(L,K))
        ENDDO
        ENDDO
      ENDIF
!
!**********************************************************************C
!
! **  INTERMEDIATE ADVECTION CALCULATIONS
! **  CZX,CZY
!
      DO K=1,KC
      DO L=2,LA
       DELCX(L,K)=CONCZ(L,K)-CONCZ(L-1   ,K)
       DELCY(L,K)=CONCZ(L,K)-CONCZ(LSC(L),K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      CONCZX(L,K)=CONCZ(L,K)+RCOSMICX(L,K)*(COSMICXP(L,K)*DELCX(L,K)+COSMICXN(L+1,K)*DELCX(L+1,K))
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      CONCZY(L,K)=CONCZ(L,K)+RCOSMICY(L,K)*(COSMICYP(L,K)*DELCY(L,K)+COSMICYN(LNC(L),K)           &
                  *DELCY(LNC(L),K) )
      ENDDO
      ENDDO
!
!**********************************************************************C
!
      DO K=1,KC
      DO L=2,LA
        CONCXYZ(L,K)=( 2.*CON1(L,K)+CONCY(L,K)+CONCYZ(L,K)+CONCZ(L,K)+CONCZY(L,K) )/6.
        CONCYZX(L,K)=( 2.*CON1(L,K)+CONCZ(L,K)+CONCZX(L,K)+CONCX(L,K)+CONCXZ(L,K) )/6.
        CONCZXY(L,K)=( 2.*CON1(L,K)+CONCY(L,K)+CONCYX(L,K)+CONCX(L,K)+CONCXY(L,K) )/6.
      ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  ADVECTIVE FLUX CALCULATION
! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
! **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
! **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CONCXYZ(L-1,K)+MIN(UHDY2(L,K),0.)*CONCXYZ(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CONCYZX(LSC(L),K)+MIN(VHDX2(L,K),0.)*CONCYZX(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(W2(L,K),0.)*CONCZXY(L,K)+MIN(W2(L,K),0.)*CONCZXY(L,K+1)
      ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  STANDARD ADVECTION CALCULATION
!
!----------------------------------------------------------------------C
!
! **  IF ISACAC EQ 0 INCLUDE FQC MASS SOURCES IN UPDATE
!
! BEGIN IF ON TRANSPORT OPTION CHOICE
!
      IF(ISCDCA(MVAR).EQ.0.)THEN
!
! BEGIN IF ON TIME LEVEL CHOICE FOR ISCDAC=0
!
      IF(ISTL.EQ.2)THEN
!
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
!X      CH(L,K)=CONT(L,K)*H1P(L)
!X     &       +DELT*( (WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*H1P(L)+DELT*((RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)            &
              -FVHU(LNC(L),K))*DXYIP(L)+(FWU(L,K-1)-FWU(L,K))*RDZIC )
      ENDDO
      ENDDO
!
      IF(ISFCT(MVAR).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON1(L,K)
       ENDDO
       ENDDO
      ENDIF
!
! ELSE ON TIME LEVEL CHOICE FOR ISCDAC=0
!
      ELSE
!
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
!X      CH(L,K)=CONT(L,K)*H2P(L)
!X     &       +DELT*( (WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*H2P(L)+DELT*( ( RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)            &
              -FVHU(LNC(L),K))*DXYIP(L)+(FWU(L,K-1)-FWU(L,K))*RDZIC )
      ENDDO
      ENDDO
!
      IF(ISFCT(MVAR).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON(L,K)
       ENDDO
       ENDDO
      ENDIF
!
      ENDIF
!
! ENDIF ON TIME LEVEL CHOICE FOR ISCDAC=0
!
      IF(ISUD.EQ.1)THEN
      DO K=1,KC
      DO L=2,LA
      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
      ENDDO
      ENDDO
      ENDIF
!
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      CONT(L,K)=0.0
      ENDDO
      ENDDO
!
! **  ADD REMAINING SEDIMENT SETTLING AND FLUX
!
!     IF(M.EQ.4)THEN
!      IF(ISTOPT(4).EQ.2) CALL CALSED2(ISTL,0.5,CON )
!      NTSWVD=ISWVSD
!      IF(ISTOPT(4).EQ.13) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.14.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.15) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.16.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
!     ENDIF
!
! **  IF ISACAC NE 0 DO NOT INCLUDE FQC MASS SOURCES IN UPDATE
!
! ELSE ON TRANSPORT OPTION CHOICE
!
      ELSE
!
! BEGIN IF ON TIME LEVEL CHOICE FOR ISCDAC.NE.0
!
      IF(ISTL.EQ.2)THEN
!
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
!X      CH(L,K)=CONT(L,K)*H1P(L)
!X     &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*H1P(L)+DELT*((RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K) +FVHU(L,K)             &
              -FVHU(LNC(L),K))*DXYIP(L)+(FWU(L,K-1)-FWU(L,K))*RDZIC)
      ENDDO
      ENDDO
!
      IF(ISFCT(MVAR).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON1(L,K)
       ENDDO
       ENDDO
      ENDIF
!
! ELSE ON TIME LEVEL CHOICE FOR ISCDAC.NE.0
!
      ELSE
!
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
!X      CH(L,K)=CONT(L,K)*H2P(L)
!X     &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*H2P(L)+DELT*((RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)              &
              -FVHU(LNC(L),K))*DXYIP(L)+(FWU(L,K-1)-FWU(L,K))*RDZIC)
      ENDDO
      ENDDO
!
      IF(ISFCT(MVAR).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON(L,K)
       ENDDO
       ENDDO
      ENDIF
!
      ENDIF
!
! ENDIF ON TIME LEVEL CHOICE FOR ISCDAC.NE.0
!
      IF(ISUD.EQ.1)THEN
      DO K=1,KC
      DO L=2,LA
      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
      ENDDO
      ENDDO
      ENDIF
!
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      CONT(L,K)=0.0
      ENDDO
      ENDDO
!
      ENDIF
!
! ENDIF ON TRANSPORT OPTION CHOICE
!
!**********************************************************************C
!
! **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW
! **  CONCENTRATION AT OPEN BOUNDARIES
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBS
      NSID=NCSERS(LL,M)
      L=LCBS(LL)
      LN=LNC(L)
!
      IF(VHDX2(LN,K).LT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)-FVHU(LN,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)-FVHU(LN,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L) &
                              +DELT*(0.5*VHDX2(LN,K)*(CON1(L,K)+CON(L,K))-FVHU(LN,K))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBSTMP=CBS(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBSTMP) CON(L,K)=CBSTMP
       CLOS(LL,K,M)=CON(L,K)
       NLOS(LL,K,M)=N
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
!      WRITE(6,6001)N,K,CBT
       NMNLO=N-NLOS(LL,K,M)
       IF(NMNLO.GE.NTSCRS(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOS(LL,K,M)+(CBT-CLOS(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
       ENDIF
      ENDIF
!
      ENDDO
      ENDDO
!
 6001 FORMAT('N,K,CBTS = ',2I10,F12.3)
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBW
      NSID=NCSERW(LL,M)
      L=LCBW(LL)
!
      IF(UHDY2(L+1,K).LT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(UHDY2(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2)CTMP=CON1(L,K)+DELT*(UHDY2(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)   &
                                 *HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)      &
                 *HPI(L)+DELT*(0.5*UHDY2(L+1,K)*(CON1(L,K)+CON(L,K))-FUHU(L+1,K))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBWTMP=CBW(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBWTMP) CON(L,K)=CBWTMP
       CLOW(LL,K,M)=CON(L,K)
       NLOW(LL,K,M)=N
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOW(LL,K,M)
       IF(NMNLO.GE.NTSCRW(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOW(LL,K,M)+(CBT-CLOW(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
       ENDIF
      ENDIF
!!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBE
      NSID=NCSERE(LL,M)
      L=LCBE(LL)
!
      IF(UHDY2(L,K).GT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(FUHU(L,K)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FUHU(L,K)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L) &
                              +DELT*(FUHU(L,K)-0.5*UHDY2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBETMP=CBE(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBETMP) CON(L,K)=CBETMP
       CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOE(LL,K,M)
       IF(NMNLO.GE.NTSCRE(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOE(LL,K,M)+(CBT-CLOE(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
       ENDIF
      ENDIF
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBN
      NSID=NCSERN(LL,M)
      L=LCBN(LL)
      LS=LSC(L)
!
      IF(VHDX2(L,K).GT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(FVHU(L,K)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FVHU(L,K)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L) &
                              +DELT*(FVHU(L,K)-0.5*VHDX2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBNTMP=CBN(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBNTMP) CON(L,K)=CBNTMP
       CLON(LL,K,M)=CON(L,K)
       NLON(LL,K,M)=N
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
!      WRITE(6,6002)N,K,CBT
       NMNLO=N-NLON(LL,K,M)
       IF(NMNLO.GE.NTSCRN(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLON(LL,K,M)+(CBT-CLON(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
       ENDIF
      ENDIF
!
      ENDDO
      ENDDO
!
 6002 FORMAT('N,K,CBTN = ',2I10,F12.3)
!
!**********************************************************************C
!
! **  MODIFIY VERTICAL MASS DIFFUSION IF ANTI-DIFFUSIVE ADVECTIVE
! **  IS TURNED OFF
!
!      IF(ISADAC(1).EQ.0.AND.ISCDCA(1).NE.1)THEN
!      IF(KC.GE.2)THEN
!
!      DO K=1,KS
!      DO L=2,LA
!       AB(L,K)=AB(L,K)-W2(L,K)*DZG(K)
!       AB(L,K)=MAX(AB(L,K),ABO)
!      ENDDO
!      ENDDO
!
!      ENDIF
!      ENDIF
!
!**********************************************************************C
!
! **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION
!
      IF(ISADAC(MVAR).EQ.0) RETURN
      IF(ISCDCA(MVAR).EQ.1) RETURN
!
!**********************************************************************C
!
! **  STANDARD ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      AUHU=ABS(UHDY2(L,K))
      AVHV=ABS(VHDX2(L,K))
      UTERM=AUHU*(1.-DELTA*AUHU*DXYIU(L)*HUI(L))*(CON(L,K)-CON(L-1,K))
      VTERM=AVHV*(1.-DELTA*AVHV*DXYIV(L)*HVI(L))*(CON(L,K)-CON(LS,K))
      IF(ISADAC(MVAR).GE.2)THEN
        SSCORUE=DELTA*RDZIC*DXYIP(L)*HPI(L)*(FQCPAD(L,K)-QSUMPAD(L,K)*CON(L,K))
        SSCORUW=DELTA*RDZIC*DXYIP(L-1)*HPI(L-1)*(FQCPAD(L-1,K)-QSUMPAD(L-1,K)*CON(L-1,K))
        SSCORVN=DELTA*RDZIC*DXYIP(L)*HPI(L)*(FQCPAD(L,K)-QSUMPAD(L,K)*CON(L,K))
        SSCORVS=DELTA*RDZIC*DXYIP(LS)*HPI(LS)*(FQCPAD(LS,K)-QSUMPAD(LS ,K)*CON(LS ,K))
        SSCORU=MAX(UHDY2(L,K),0.0)*SSCORUW+MIN(UHDY2(L,K),0.0)*SSCORUE
        SSCORV=MAX(VHDX2(L,K),0.0)*SSCORVS+MIN(VHDX2(L,K),0.0)*SSCORVN
	  UTERM=UTERM+SSCORU
	  VTERM=VTERM+SSCORV
	ENDIF
      IF(ISFCT(MVAR).GE.2)THEN
        FUHU(L,K)=0.5*UTERM
        FVHU(L,K)=0.5*VTERM
        IF(ISFCT(MVAR).EQ.3)THEN
          FUHU(L,K)=UTERM
          FVHU(L,K)=VTERM
        ENDIF
      ELSE
      UHU=UTERM/(CON(L,K)+CON(L-1,K)+BSMALL)
      VHV=VTERM/(CON(L,K)+CON(LS,K)+BSMALL)
      FUHU(L,K)=MAX(UHU,0.)*CON(L-1,K)+MIN(UHU,0.)*CON(L,K)
      FVHU(L,K)=MAX(VHV,0.)*CON(LS,K)+MIN(VHV,0.)*CON(L,K)
      ENDIF
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(W2(L,K))
      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)*HPI(L))*(CON(L,K+1)-CON(L,K))
      IF(ISADAC(MVAR).GE.2)THEN
        SSCORWA=DELTA*DZIG(K+1)*HPI(L)*DXYIP(L)*(FQCPAD(L,K+1)-QSUMPAD(L,K+1)*CON(L,K+1))
        SSCORWB=DELTA*DZIG(K)*HPI(L)*DXYIP(L)*(FQCPAD(L,K  )-QSUMPAD(L,K  )*CON(L,K  ))
        SSCORW=MAX(W2(L,K),0.0)*SSCORWB+MIN(W2(L,K),0.0)*SSCORWA
	  WTERM=WTERM+SSCORW
      ENDIF
      IF(ISFCT(MVAR).GE.2)THEN
        FWU(L,K)=0.5*WTERM
        IF(ISFCT(MVAR).EQ.3)THEN
          FWU(L,K)=WTERM
        ENDIF
      ELSE
        WW=WTERM/(CON(L,K+1)+CON(L,K)+BSMALL)
        FWU(L,K)=MAX(WW,0.)*CON(L,K)+MIN(WW,0.)*CON(L,K+1)
      ENDIF
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! ** SET ANTIDIFFUSIVE FLUXES TO ZERO FOR SOURCE CELLS JMH20MARCH97
!
!LH
!LH   DO K=1,KC
!LH   DO L=2,LA
!LH     CONT(L,K)=1.
!LH     ABSQSUM=ABS(QSUM(L,K))
!LH     IF(ABSQSUM.GT.1.E-12) CONT(L,K)=0.
!LH   ENDDO
!LH   ENDDO
!
!LH   DO K=1,KC
!LH   DO L=2,LA
!LH       TMPVAL=CONT(L,K)
!LH       LN=LNC(L)
!LH       FUHU(L  ,K)=TMPVAL*FUHU(L  ,K)
!LH       FUHU(L+1,K)=TMPVAL*FUHU(L+1,K)
!LH       FVHU(L  ,K)=TMPVAL*FVHU(L  ,K)
!LH       FVHU(LN ,K)=TMPVAL*FVHU(LN ,K)
!LH       FWU(L,K  )=TMPVAL*FWU(L,K  )
!LH       FWU(L,K-1)=TMPVAL*FWU(L,K-1)
!LH   ENDDO
!LH   ENDDO
!
!  ** ANTIDIFFUSION TURNED OFF FOR SOURCE CELLS
!
      IF(ISADAC(MVAR).EQ.1)THEN
      DO K=1,KC
      DO L=2,LA
!        IF(ABS(QSUM(L,K)).GT.1.E-12)THEN
        IF(QSUMPAD(L,K).GT.0.0)THEN
          LN=LNC(L)
          FUHU(L  ,K)=0.
          FUHU(L+1,K)=0.
          FVHU(L  ,K)=0.
          FVHU(LN ,K)=0.
          FWU(L,K  )=0.
          FWU(L,K-1)=0.
          CONT(L,K)=0.
        ELSE
          CONT(L,K)=1.
        ENDIF
      ENDDO
      ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! ** SET ANTIDIFFUSIVE FLUXES TO ZERO FOR OPEN BOUNDARY CELLS
!
      DO K=1,KC
!
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
!     DO K=1,KC
      FVHU(LN,K)=0.0
!     ENDDO
      ENDDO
!
      DO LL=1,NCBW
      L=LCBW(LL)
!     DO K=1,KC
      FUHU(L+1,K)=0.0
!     ENDDO
      ENDDO
!
      DO LL=1,NCBE
      L=LCBE(LL)
!     DO K=1,KC
      FUHU(L,K)=0.0
!     ENDDO
      ENDDO
!
      DO LL=1,NCBN
      L=LCBN(LL)
!     DO K=1,KC
      FVHU(L,K)=0.0
!     ENDDO
      ENDDO
!
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE AND APPLY FLUX CORRECTED TRANSPORT LIMITERS
!
      IF(ISFCT(MVAR).EQ.0) GOTO 1100
!
!----------------------------------------------------------------------C
!
! **  DETERMINE MAX AND MIN CONCENTRATIONS
!
      DO K=1,KC
      DO L=1,LC
      CONTMX(L,K)=0.0
      CONTMN(L,K)=0.0
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      CONTMX(L,K)=MAX(CON(L,K),CON2(L,K))
      CONTMN(L,K)=MIN(CON(L,K),CON2(L,K))
      ENDDO
      ENDDO
!
      DO L=2,LA
      CMAX(L,1)=MAX(CONTMX(L,1),CONTMX(L,2))
      CMAX(L,KC)=MAX(CONTMX(L,KS),CONTMX(L,KC))
      CMIN(L,1)=MIN(CONTMN(L,1),CONTMN(L,2))
      CMIN(L,KC)=MIN(CONTMN(L,KS),CONTMN(L,KC))
      ENDDO
!
      DO K=2,KS
      DO L=2,LA
      CMAXT=MAX(CONTMX(L,K-1),CONTMX(L,K+1))
      CMAX(L,K)=MAX(CONTMX(L,K),CMAXT)
      CMINT=MIN(CONTMN(L,K-1),CONTMN(L,K+1))
      CMIN(L,K)=MIN(CONTMN(L,K),CMINT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      CWMAX=SUB(L)*CONTMX(L-1,K)
      CEMAX=SUB(L+1)*CONTMX(L+1,K)
      CSMAX=SVB(L)*CONTMX(LS,K)
      CNMAX=SVB(LN)*CONTMX(LN,K)
      CMAXT=MAX(CNMAX,CEMAX)
      CMAXT=MAX(CMAXT,CSMAX)
      CMAXT=MAX(CMAXT,CWMAX)
      CMAX(L,K)=MAX(CMAX(L,K),CMAXT)
      CWMIN=SUB(L)*CONTMN(L-1,K)+1.E+6*(1.-SUB(L))
      CEMIN=SUB(L+1)*CONTMN(L+1,K)+1.E+6*(1.-SUB(L+1))
      CSMIN=SVB(L)*CONTMN(LS,K)+1.E+6*(1.-SVB(L))
      CNMIN=SVB(LN)*CONTMN(LN,K)+1.E+6*(1.-SVB(LN))
      CMINT=MIN(CNMIN,CEMIN)
      CMINT=MIN(CMINT,CSMIN)
      CMINT=MIN(CMINT,CWMIN)
      CMIN(L,K)=MIN(CMIN(L,K),CMINT)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  SEPARATE POSITIVE AND NEGATIVE FLUXES PUTTING NEGATIVE FLUXES
! **  INTO FUHV, FVHV, AND FWV
!
      DO K=1,KC
      DO L=2,LA
      FUHV(L,K)=MIN(FUHU(L,K),0.)
      FUHU(L,K)=MAX(FUHU(L,K),0.)
      FVHV(L,K)=MIN(FVHU(L,K),0.)
      FVHU(L,K)=MAX(FVHU(L,K),0.)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      FWV(L,K)=MIN(FWU(L,K),0.)
      FWU(L,K)=MAX(FWU(L,K),0.)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  CALCULATE INFLUX AND OUTFLUX IN CONCENTRATION UNITS AND LOAD
! **  INTO DU AND DV, THEN ADJUCT VALUES AT BOUNDARIES
!
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
      LN=LNC(L)
      DU(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L,K)-FUHV(L+1,K)+FVHU(L,K)-FVHV(LN,K))+RDZIC           &
              *(FWU(L,K-1)-FWV(L,K)) )*HPI(L)
      DV(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L+1,K)-FUHV(L,K)+FVHU(LN,K)-FVHV(L,K))+RDZIC           &
              *(FWU(L,K)-FWV(L,K-1)) )*HPI(L)
      ENDDO
      ENDDO
!
      DO K=1,KC
!
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
!     DO K=1,KC
      DU(LN,K)=0.
      DV(LN,K)=0.
!     ENDDO
      ENDDO
!
      DO LL=1,NCBW
      L=LCBW(LL)
!     DO K=1,KC
      DU(L+1,K)=0.
      DV(L+1,K)=0.
!     ENDDO
      ENDDO
!
      DO LL=1,NCBE
      L=LCBE(LL)
!     DO K=1,KC
      DU(L-1,K)=0.
      DV(L-1,K)=0.
!     ENDDO
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
!     DO K=1,KC
      DU(LS,K)=0.
      DV(LS,K)=0.
!     ENDDO
      ENDDO
!
      ENDDO
!
!X      DO NS=1,NQSIJ
!X       L=LQS(NS)
!X       NQSTMP=NQSERQ(NS)
!X       DO K=1,KC
!X        QQQTMP=ABS(QSS(K,NS)+QSERT(K,NQSTMP))
!X        IF(QQQTMP.GE.1.E-12)THEN
!X          DU(L,K)=0.
!X          DV(L,K)=0.
!         DU(L+1,K)=0.
!         DV(LNC(L),K)=0.
!X        ENDIF
!X       ENDDO
!X      ENDDO
!
!X      DO NCTL=1,NQCTL
!X       IU=IQCTLU(NCTL)
!X       JU=JQCTLU(NCTL)
!X       LU=LIJ(IU,JU)
!X       ID=IQCTLD(NCTL)
!X       JD=JQCTLD(NCTL)
!X       IF(ID.EQ.0.AND.JD.EQ.0)THEN
!X         LD=LC
!X        ELSE
!X         LD=LIJ(ID,JD)
!X       ENDIF
!X       DO K=1,KC
!X        QQQTMP=ABS(QCTLT(K,NCTL))
!X        IF(QQQTMP.GE.1.E-12)THEN
!X          DU(LU,K)=0.
!X          DV(LU,K)=0.
!         DU(LU+1,K)=0.
!         DV(LNC(LU),K)=0.
!X          DU(LD,K)=0.
!X          DV(LD,K)=0.
!         DU(LD+1,K)=0.
!         DV(LNC(LD),K)=0.
!X        ENDIF
!X       ENDDO
!X      ENDDO
!
!X      DO NWR=1,NQWR
!X       IU=IQWRU(NWR)
!X       JU=JQWRU(NWR)
!X       KU=KQWRU(NWR)
!X       ID=IQWRD(NWR)
!X       JD=JQWRD(NWR)
!X       KD=KQWRD(NWR)
!X       LU=LIJ(IU,JU)
!X       LD=LIJ(ID,JD)
!X        NQSTMP=NQWRSERQ(NWR)
!X        QQQTMP=ABS(QWR(NWR)+QWRSERT(NQSTMP))
!X        IF(QQQTMP.GE.1.E-12)THEN
!X          DU(LU,KU)=0.
!X          DV(LU,KU)=0.
!         DU(LU+1,K)=0.
!         DV(LNC(LU),K)=0.
!X          DU(LD,KD)=0.
!X          DV(LD,KD)=0.
!         DU(LD+1,K)=0.
!         DV(LNC(LD),K)=0.
!X        ENDIF
!X      ENDDO
!
!----------------------------------------------------------------------C
!
! **  CALCULATE BETA COEFFICIENTS WITH BETAUP AND BETADOWN IN DU AND DV
!
      DO K=1,KC
      DO L=2,LA
      IF(DU(L,K).GT.0.) DU(L,K)=(CMAX(L,K)-CON(L,K))/(DU(L,K)+BSMALL)
      DU(L,K)=MIN(DU(L,K),1.)
      IF(DV(L,K).GT.0.) DV(L,K)=(CON(L,K)-CMIN(L,K))/(DV(L,K)+BSMALL)
      DV(L,K)=MIN(DV(L,K),1.)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  LIMIT FLUXES
!
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,K)=MIN(DV(L-1,K),DU(L,K))*FUHU(L,K)+MIN(DU(L-1,K),DV(L,K))*FUHV(L,K)
      FVHU(L,K)=MIN(DV(LS,K),DU(L,K))*FVHU(L,K)+MIN(DU(LS,K),DV(L,K))*FVHV(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MIN(DV(L,K),DU(L,K+1))*FWU(L,K)+MIN(DU(L,K),DV(L,K+1))*FWV(L,K)
      ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  ANTI-DIFFUSIVE ADVECTION CALCULATION
!
 1100 CONTINUE
!
!----------------------------------------------------------------------C
!
! **  ADD SECOND CALL TO EXPERIMENTAL SED CALCULATION FOR
! **  ANTIDIFFUSIVE CORRECTION
! !!  THIS IS A TEMPORARY TRIAL IF DOESN'T WORK MAKE CALL
! !!  FOR ALL CASES AT PREVIOUS LOCATION
! !!  NOT ALSO THAT FQC HAS BEEN MOVED TO HERE FOR ISADAC=1
! !!  IF PROBLEMS ARISE PUT BACK IN AS BEFORE AT OLD LOCATION
! !!  RIGHT AFTER SECOND CALL TO CALSED2
!
!     IF(M.EQ.4) CALL CALSED2(ISTL,0.5,CON)
!
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
!X      CH(L,K)=CON(L,K)*HP(L)
!X     &       +DELT*( (WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON(L,K)*HP(L)+DELT*((FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)      &
              +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      ENDDO
      ENDDO
!
!     CH(L,K)=CONT(L,K)*HP(L)
!    &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
!    &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
!    &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
!
!RAY      DO K=1,KC
!RAY      DO L=2,LA
!RAY      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
!RAY      ENDDO
!RAY      ENDDO
!
! **  ADD REMAINING SEDIMENT SETTLING AND FLUX
!
!     IF(M.EQ.4)THEN
!      IF(ISTOPT(6).EQ.2) CALL CALSED2(ISTL,0.5)
!      NTSWVD=ISWVSD
!      IF(ISTOPT(4).EQ.13) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.14.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.15) CALL CALSED3(ISTL,0.5,SED)
!      IF(ISTOPT(4).EQ.16.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
!     ENDIF
!
!**********************************************************************C
!
! **  DIAGNOSE FCT SCHEME
!
!----------------------------------------------------------------------C
!
      IF(ISFCT(MVAR).EQ.99)THEN
      WRITE(6,6110)N
!
      DO K=1,KC
      DO L=2,LA
      CCMAX=SCB(L)*(CON(L,K)-CMAX(L,K))
      IF(CCMAX.GT.0.)THEN
       WRITE(6,6111)CON(L,K),CMAX(L,K),IL(L),JL(L),K
      ENDIF
      CCMIN=SCB(L)*(CMIN(L,K)-CON(L,K))
      IF(CCMIN.GT.0.)THEN
       WRITE(6,6112)CMIN(L,K),CON(L,K),IL(L),JL(L),K
      ENDIF
      ENDDO
      ENDDO
!
      ENDIF
!
 6110 FORMAT('  FCT DIAGNOSTICS AT N = ',I5)
 6111 FORMAT('  CON = ',E12.4,3X,'CMAX = ',E12.4,3X,'I,J,K=',(3I10))
 6112 FORMAT('  CMIN = ',E12.4,3X,'CON = ',E12.4,3X,'I,J,K=',(3I10))
!
!----------------------------------------------------------------------C
!
      RETURN
      END