!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALTRAN1D (ISTL,MVAR,M,CON,CON1)
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
! **  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE
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
      DIMENSION FQCPAD(LCM,KCM),QSUMPAD(LCM,KCM),QSUMNAD(LCM,KCM)
!     DIMENSION  CH(LCM,KCM), CMAX(LCM,KCM), CMIN(LCM,KCM),
!    &         CONT(LCM,KCM), CON2(LCM,KCM),
!
!**********************************************************************C
!
!DIA      IF(N.EQ.1)THEN
!DIA        OPEN(1,FILE='CALTRAN1D.DIA')
!DIA        CLOSE(1,STATUS='DELETE')
!DIA        OPEN(1,FILE='CALTRAN1D.DIA')
!DIA        WRITE(1,1110)N,ISTL
!DIA       ELSE
!DIA        OPEN(1,FILE='CALTRAN1D.DIA',POSITION='APPEND')
!DIA        WRITE(1,1110)N,ISTL
!DIA      ENDIF
!
 1110 FORMAT(' N,ISTL = ',2I5)
!
      BSMALL=1.0E-6
!BUG    WFQC=0.0
!
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
       ISUD=1
      ENDIF
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
        FUHV(L,K)=0.
        UUU(L,K)=0.
        VVV(L,K)=0.
        DU(L,K)=0.
        DV(L,K)=0.
       ENDDO
      ENDDO
!
       DO L=1,LC
        IF(LCT(L).EQ.6) AREANEW(L)=FADYP(L)
        IF(LCT(L).EQ.6) AREAOLD(L)=FADYP1(L)
        IF(LCT(L).EQ.7) AREANEW(L)=FADXP(L)
        IF(LCT(L).EQ.7) AREAOLD(L)=FADXP1(L)
       ENDDO
!
!**********************************************************************C
! 
        DO K=1,KC
         DO L=1,LC
         CONT(L,K)=0.
         CMAX(L,K)=0.
         CMIN(L,K)=0.
         CON1(L,K)=CON(L,K)
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
! **  SELECT TRANSPORT OPTION, ISPLIT=1 FOR HORIZONTAL-VERTICAL
! **  OPERATOR SPLITTING
!
!     IF(ISPLIT(M).EQ.1) GOTO 1000
!
!**********************************************************************C
!**********************************************************************C
!
! **  BEGIN COMBINED ADVECTION SCHEME
!
!**********************************************************************C
!
! **  ADVECTIVE FLUX CALCULATION
!
      IF(ISTL.EQ.2) GOTO 300
      IF(ISCDCA(MVAR).EQ.0) GOTO 300
      IF(ISCDCA(MVAR).EQ.1) GOTO 400
      IF(ISCDCA(MVAR).EQ.2) GOTO 350
!
!**********************************************************************C
!
! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
! **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
! **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
!
  300 CONTINUE
!
!----------------------------------------------------------------------C
!
!X      IF(ISTL.EQ.2)THEN
!
!X      DO K=1,KC
!X      DO L=2,LA
!X      CONT(L,K)=CON1(L,K)+DELT*0.5*FQC(L,K)*DXYIP(L)/H1P(L)
!X      ENDDO
!X      ENDDO
!
!X      ELSE
!
!X      DO K=1,KC
!X      DO L=2,LA
!X      CONT(L,K)=CON1(L,K)+DELT*0.5*FQC(L,K)*DXYIP(L)/H2P(L)
!X      ENDDO
!X      ENDDO
!
!X      ENDIF
!
!X      WFQC=0.5
!
!----------------------------------------------------------------------C
!
!X      DO K=1,KC    
!X      DO L=2,LA
!X      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CONT(L-1,K)
!X     &         +MIN(UHDY2(L,K),0.)*CONT(L,K)
!X      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CONT(LSC(L),K)
!X     &         +MIN(VHDX2(L,K),0.)*CONT(L,K)
!X      ENDDO
!X      ENDDO
!X      DO K=1,KS
!X      DO L=2,LA
!X      FWU(L,K)=MAX(W2(L,K),0.)*CONT(L,K)
!X     &        +MIN(W2(L,K),0.)*CONT(L,K+1)
!X      ENDDO
!X      ENDDO
!
!1D      DO K=1,KC    
!1D      DO L=2,LA
!1D      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CON1(L-1,K)
!1D     &         +MIN(UHDY2(L,K),0.)*CON1(L,K)
!1D      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CON1(LSC(L),K)
!1D     &         +MIN(VHDX2(L,K),0.)*CON1(L,K)
!1D      ENDDO
!1D      ENDDO
!1D      DO K=1,KS
!1D      DO L=2,LA
!1D      FWU(L,K)=MAX(W2(L,K),0.)*CON1(L,K)
!1D     &        +MIN(W2(L,K),0.)*CON1(L,K+1)
!1D      ENDDO
!1D      ENDDO
!
      DO K=1,KC    
      DO L=2,LA
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CON1(L-1,K)+MIN(UHDY2(L,K),0.)*CON1(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CON1(LSC(L),K)+MIN(VHDX2(L,K),0.)*CON1(L,K)
      ENDDO
      ENDDO
!
!DIA      WRITE(1,1111) IL(2),JL(2),K,AREAOLD(2),AREANEW(2),
!DIA     &             CON1(2,1),UHDY2(2,1),UHDY2(3,1),FQC(2,1)
!DIA     &             CON1(2,1),FUHU(2,1),FUHU(3,1),FQC(2,1)
!
      GOTO 500
!
 1111 FORMAT(3I5,6E12.4)
!
!**********************************************************************C
!
! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
! **  AVERAGED BETWEEN  (N-1) AND (N+1) AND ADVECTED FIELD AVERAGED 
! **  BETWEEN AT (N-1) AND (N) IF ISTL 3 ONLY
!
  350 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
      CONT(L,K)=0.5*(CON(L,K)+CON1(L,K))+DELT*0.5*FQC(L,K)*DXYIP(L)/H2P(L)
      ENDDO
      ENDDO
!
!X     WFQC=0.5
!
!----------------------------------------------------------------------C
!
      DO K=1,KC    
      DO L=2,LA
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CONT(L-1,K)+MIN(UHDY2(L,K),0.)*CONT(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CONT(LSC(L),K)+MIN(VHDX2(L,K),0.)*CONT(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(W2(L,K),0.)*CONT(L,K)
     &        +MIN(W2(L,K),0.)*CONT(L,K+1)
      ENDDO
      ENDDO
!
      GOTO 500
!
!**********************************************************************C
!
! **  CALCULATE ADVECTIVE FLUXES BY CENTRAL DIFFERENCE WITH TRANSPORT 
! **  AVERAGED BETWEEN (N+1) AND (N-1) AND TRANSPORTED FIELD AT (N)
!
  400 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
      CONT(L,K)=CON1(L,K)
      ENDDO
      ENDDO
!
!X     WFQC=1.0
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)    
      FUHU(L,K)=0.5*UHDY2(L,K)*(CON(L,K)+CON(L-1,K))
      FVHU(L,K)=0.5*VHDX2(L,K)*(CON(L,K)+CON(LS,K))
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
!
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      IF(VHDX2(LN,K).LT.0.) FVHU(LN,K)=VHDX2(LN,K)*CON1(LN,K)
      ENDDO
!
      DO LL=1,NCBW
      L=LCBW(LL)
      IF(UHDY2(L+1,K).LT.0.) FUHU(L+1,K)=UHDY2(L+1,K)*CON1(L+1,K)
      ENDDO
!
      DO LL=1,NCBE
      L=LCBE(LL)
      IF(UHDY2(L,K).GT.0.) FUHU(L,K)=UHDY2(L,K)*CON1(L-1,K)
      ENDDO
!
      DO LL=1,NCBN
      L=LCBN(LL)
      LS =LSC(L)
      IF(VHDX2(L,K).GT.0.) FVHU(L,K)=VHDX2(L,K)*CON1(LS,K)
      ENDDO
!
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=0.5*W2(L,K)*(CON(L,K+1)+CON(L,K))
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX FOR CENTRAL DIFF
!
      IF(ISHDMF.EQ.2) CALL CALDIFF (ISTL,M,CON1)
!     IF(ISTRAN(M).GE.1) CALL CALDIFF (ISTL,M,CON1)
!
!----------------------------------------------------------------------C
!
! **  ADD SECOND CALL TO EXPERIMENTAL SED CALCULATION FOR NO 
! **  ANTIDIFFUSIVE CORRECTION
!
!     IF(ISCDCA(MVAR).EQ.0)THEN
!     IF(M.EQ.4)THEN
!       IF(ISTRAN(4).LE.2) CALL CALSED2(ISTL,0.5,CON1)
!       IF(ISTRAN(4).EQ.3) CALL CALSED3(ISTL,0.5,CON1)
!     ENDIF
!     ENDIF
!
!**********************************************************************C
!
! **  STANDARD ADVECTION CALCULATION
!
  500 CONTINUE
!
!----------------------------------------------------------------------C
!
! **  IF ISACAC EQ 0 INCLUDE FQC MASS SOURCES IN UPDATE
!
! BEGIN IF ON TRANSPORT OPTION CHOICE
!
      IF(ISCDCA(MVAR).EQ.0)THEN
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
      CH(L,K)=CON1(L,K)*AREAOLD(L)+DELT*( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LNC(L),K)) &
              *DXYIP(L)
!     &       +DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
!     &                        +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
!     &                        +(FWU(L,K-1)-FWU(L,K))*RDZIC )
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
      CH(L,K)=CON1(L,K)*AREAOLD(L)+DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)              &
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
!1D      IF(ISUD.EQ.1)THEN
!1D      DO K=1,KC
!1D      DO L=2,LA
!1D      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
!1D      ENDDO
!1D      ENDDO
!1D      ENDIF
!
!1D      DO K=1,KC
!1D      DO L=2,LA
!1D      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
!1D      ENDDO
!1D      ENDDO
!
      DO K=1,KC
      DO L=2,LA
!      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      CON(L,K)=SCB(L)*(CH(L,K)/AREANEW(L))+(1.-SCB(L))*CON(L,K)
!DIA      IF(L.EQ.2) WRITE(1,1111) IL(L),JL(L),K,CON(L,K),
!DIA     &             CON1(L,K)
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
      CH(L,K)=CON1(L,K)*AREAOLD(L)+DELT*( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)                &
              -FVHU(LNC(L),K))*DXYIP(L)
!1D     &       +DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
!1D     &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
!1D     &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
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
      CH(L,K)=CON1(L,K)*AREAOLD(L)+DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)              &
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
! ENDIF ON TIME LEVEL CHOICE FOR ISCDAC.NE.0
!
!1D      IF(ISUD.EQ.1)THEN
!1D      DO K=1,KC
!1D      DO L=2,LA
!1D      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
!1D      ENDDO
!1D      ENDDO
!1D      ENDIF
!
      DO K=1,KC
      DO L=2,LA
!      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      CON(L,K)=SCB(L)*(CH(L,K)/AREANEW(L))+(1.-SCB(L))*CON(L,K)
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
  700 CONTINUE
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
        CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)-FVHU(LN,K))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)-FVHU(LN,K))*DXYIP(L)     &
                                   /AREANEW(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*AREAOLD(L)   &
             /AREANEW(L)+DELT*(0.5*VHDX2(LN,K)*(CON1(L,K)+CON(L,K))-FVHU(LN,K))*DXYIP(L)/AREANEW(L)
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
        CTMP=CON1(L,K)+DELT*(UHDY2(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(UHDY2(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)   &
                                   /AREANEW(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*AREAOLD(L)   &
           /AREANEW(L)+DELT*(0.5*UHDY2(L+1,K)*(CON1(L,K)+CON(L,K))-FUHU(L+1,K))*DXYIP(L)/AREANEW(L)
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
!
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
        CTMP=CON1(L,K)+DELT*(FUHU(L,K)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)/AREANEW(L)
!DIA        WRITE(1,3331)N,UHDY2(L,K),CTMP,CON1(L,K),CON(L,K)
        CON1(L,K)=CON(L,K)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FUHU(L,K)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)       &
                                   /AREANEW(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*AREAOLD(L)   &
               /AREANEW(L)+DELT*(FUHU(L,K)-0.5*UHDY2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBETMP=CBE(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBETMP) CON(L,K)=CBETMP
       CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
!DIA       WRITE(1,3332)N,NLOE(LL,K,M),CBETMP,CLOE(LL,K,M)
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOE(LL,K,M)
       IF(NMNLO.GE.NTSCRE(LL))THEN
        CON(L,K)=CBT
!DIA        WRITE(1,3333)N,NMNLO,NLOE(LL,K,M),CBT
       ELSE
        CON(L,K)=CLOE(LL,K,M)+(CBT-CLOE(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
!DIA        WRITE(1,3334)N,NMNLO,NLOE(LL,K,M),CBT,CLOE(LL,K,M)
       ENDIF
      ENDIF
!
      ENDDO
      ENDDO
!
 3331 FORMAT(' W1 ',I5,4F12.4)
 3332 FORMAT(' W2 ',2I5,4F12.4)
 3333 FORMAT(' W3 ',2I5,4F12.4)
 3334 FORMAT(' W4 ',3I5,4F12.4)
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
        CTMP=CON1(L,K)+DELT*(FVHU(L,K)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FVHU(L,K)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)       &
                                   /AREANEW(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*AREAOLD(L)   &
           /AREANEW(L)+DELT*(FVHU(L,K)-0.5*VHDX2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)/AREANEW(L)
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
!DIA      CLOSE(1)
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
      LS=LSC(L)
      UUU(L,K)=U2(L,K)*(CON(L,K)-CON(L-1,K))*DXIU(L)
      VVV(L,K)=V2(L,K)*(CON(L,K)-CON(LS,K))*DYIV(L)
      ENDDO
      ENDDO
!
      DO K=1,KS
      RDZIG=DZIG(K)
      DO L=2,LA
      WWW(L,K)=W2(L,K)*(CON(L,K+1)-CON(L,K))*HPI(L)*RDZIG
      ENDDO
      ENDDO
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
!1D      UTERM=AUHU*(1.-DELTA*AUHU*DXYIU(L)*HUI(L))*(CON(L,K)-CON(L-1,K))
!1D      VTERM=AVHV*(1.-DELTA*AVHV*DXYIV(L)*HVI(L))*(CON(L,K)-CON(LS,K))
      UTERM=AUHU*(1.-DELTA*AUHU*DXYIU(L)*HUI(L))*(CON(L,K)-CON(L-1,K))
      VTERM=AVHV*(1.-DELTA*AVHV*DXYIV(L)*HVI(L))*(CON(L,K)-CON(LS,K))
!     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
!     AHU(L,K)=MAX(AHTMP,0.)
!     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
!     AHV(L,K)=MAX(AHTMP,0.)
!     UTERM=UTERM*(CON(L,K)-CON(L-1,K))
!     VTERM=VTERM*(CON(L,K)-CON(LS,K))
      UTERM=UTERM-DELTA4*UHDY2(L,K)*(VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K)+WWW(L,K)+WWW(L-1,K) &
             +WWW(L-1,K-1)+WWW(L,K-1))
      VTERM=VTERM-DELTA4*VHDX2(L,K)*(UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K)+WWW(L,K)+WWW(LS,K)  &
            +WWW(LS,K-1)+WWW(L,K-1))
      UHU=UTERM/(CON(L,K)+CON(L-1,K)+BSMALL)
      VHV=VTERM/(CON(L,K)+CON(LS,K)+BSMALL)
!     AUHU=ABS(UHU)
!     AVHV=ABS(VHV)
!     UTERM=AUHU*(1.-DELTD2*AUHU/(DXYU(L)*HU(L)))
!     VTERM=AVHV*(1.-DELTD2*AVHV/(DXYV(L)*HV(L)))
!     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
!     AHU(L,K)=MAX(AHTMP,0.)
!     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
!     AHV(L,K)=MAX(AHTMP,0.)
      FUHU(L,K)=MAX(UHU,0.)*CON(L-1,K)+MIN(UHU,0.)*CON(L,K)
      FVHU(L,K)=MAX(VHV,0.)*CON(LS,K)+MIN(VHV,0.)*CON(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(W2(L,K))
      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)*HPI(L))*(CON(L,K+1)-CON(L,K))
      WTERM=WTERM-DELTA4*W2(L,K)*(UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)+VVV(L,K)+VVV(LN,K)  &
            +VVV(LN,K+1)+VVV(L,K+1))
      WW=WTERM/(CON(L,K+1)+CON(L,K)+BSMALL)
      FWU(L,K)=MAX(WW,0.)*CON(L,K)+MIN(WW,0.)*CON(L,K+1)
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
      DO K=1,KC
!DIR& IVDEP
      DO L=2,LA
        IF(ABS(QSUM(L,K)).GT.1.E-12)THEN
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
      DV(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L+1,K)-FUHV(L,K)+FVHU(LN,K)-FVHV(L,K))+RDZIC*(FWU(L,K) &
              -FWV(L,K-1)) )*HPI(L)
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
!
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
      CH(L,K)=CON(L,K)*HP(L)+DELT*( (FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)     &
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
      IF(ISFCT(MVAR).EQ.2)THEN
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