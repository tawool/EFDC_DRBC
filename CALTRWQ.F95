!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALTRWQ (M,NW,CON,CON1)
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
!    &         CONT(LCM,KCM), CON2(LCM,KCM)
!
!**********************************************************************C
!
      MVAR=M
!
      BSMALL=1.0E-6
!
      DELT=DT2
      DELTA=DT2
      DELTD2=DT
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      ISTL=3
!
      IF(IS2TIM.GE.1) THEN
	  ISTL=2
        IF(ISDYNSTP.EQ.0)THEN
          DELT=DT
          DELTA=DT
          DELTD2=0.5*DT
          S3TL=0.0
          S2TL=1.0
          ISUD=0
         ELSE
          DELT=DTDYN
          DELTA=DTDYN
          DELTD2=0.5*DTDYN
          S3TL=0.0
          S2TL=1.0
          ISUD=0
        END IF
	ENDIF
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
!      CALL CALFQC (ISTL,MVAR,M,CON,CON1)
      CALL CALFQC (ISTL,0,MVAR,M,CON,CON1,FQCPAD,QSUMPAD,QSUMNAD) 
!
!WQJH     IF(ISTRAN(M).GE.1)THEN
!WQJH       IF(M.EQ.7) CALL CALFQC (ISTL,M,CON,CON1)
!WQJH       IF(M.EQ.8)THEN
!WQJH         DO K=1,KC 
!WQJH          DO L=2,LA
!WQJH           FQC(L,K)=0.
!WQJH          ENDDO
!WQJH         ENDDO
!WQJH       ENDIF
!WQJH     ENDIF
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
! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
! **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
! **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
!WQJH     CONT(L,K)=CON1(L,K)+DELT*0.5*FQC(L,K)*DXYIP(L)/H2WQ(L)
      CONT(L,K)=CON1(L,K)
      ENDDO
      ENDDO
!
!WQJH     WFQC=0.5
!
!----------------------------------------------------------------------C
!
      DO K=1,KC    
      DO L=2,LA
      LS=LSC(L)      
      FUHU(L,K)=MAX(UHDYWQ(L,K),0.)*CONT(L-1,K)+MIN(UHDYWQ(L,K),0.)*CONT(L,K)
      FVHU(L,K)=MAX(VHDXWQ(L,K),0.)*CONT(LS,K)+MIN(VHDXWQ(L,K),0.)*CONT(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(WWQ(L,K),0.)*CONT(L,K)+MIN(WWQ(L,K),0.)*CONT(L,K+1)
      ENDDO
      ENDDO
!
!     DO K=1,KS
!     DO L=2,LA
!     FWU(L,K)=0.5*WWQ(L,K)*(CON(L,K+1)+CON(L,K))
!     ENDDO
!     ENDDO
!
!**********************************************************************C
!
! **  CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX
!
!----------------------------------------------------------------------C
!
!     IF(ISTRAN(M).GE.1) CALL CALDIFF (ISTL,M,CON1)
!
!**********************************************************************C
!
! **  ADVECTION CALCULATION
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
      LN=LNC(L)      
!WQJH      CH(L,K)=CONT(L,K)*H2WQ(L)
!WQJH     &       +DELT*((WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CONT(L,K)*H2WQ(L)+DELT*((RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)) &
              *DXYIP(L)+(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      ENDDO
      ENDDO
!
      IF(ISFCT(M).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON1(L,K)
       ENDDO
       ENDDO
      ENDIF
!
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)/HWQ(L)+(1.-SCB(L))*CON1(L,K)
      CONT(L,K)=0.0
      ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW 
! **  CONCENTRATION AT OPEN BOUNDARIES FORM M=4
!
      IF(M.EQ.4)THEN
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBS
      NSID=NCSERS(LL,M)
      L=LCBS(LL)
      LN=LNC(L)
!
      IF(VHDXWQ(LN,K).LT.0.)THEN
       CTMP=CON1(L,K)+DELT*(VHDXWQ(LN,K)*CON1(L,K)-FVHU(LN,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
!      IF(CON(L,K).GT.CBS(LL,1,M)) CON(L,K)=CBS(LL,1,M)
       CLOS(LL,K,M)=CON(L,K)        
       NLOS(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
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
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBW
      NSID=NCSERW(LL,M)
      L=LCBW(LL)      
!
      IF(UHDYWQ(L+1,K).LT.0.)THEN
       CTMP=CON1(L,K)+DELT*(UHDYWQ(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
!      IF(CON(L,K).GT.CBW(LL,1,M)) CON(L,K)=CBW(LL,1,M)
       CLOW(LL,K,M)=CON(L,K)
       NLOW(LL,K,M)=N
      ELSE
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
      IF(UHDYWQ(L,K).GT.0.)THEN
       CTMP=CON1(L,K)+DELT*(FUHU(L,K)-UHDYWQ(L,K)*CON1(L,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
!      IF(CON(L,K).GT.CBE(LL,1,M)) CON(L,K)=CBE(LL,1,M)
       CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
      ELSE
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
      IF(VHDXWQ(L,K).GT.0.)THEN
       CTMP=CON1(L,K)+DELT*(FVHU(L,K)-VHDXWQ(L,K)*CON1(L,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
!      IF(CON(L,K).GT.CBN(LL,1,M)) CON(L,K)=CBN(LL,1,M)
       CLON(LL,K,M)=CON(L,K)
       NLON(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
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
!----------------------------------------------------------------------C
!
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW 
! **  CONCENTRATION AT OPEN BOUNDARIES FOR WQVARS, M=8
!
      IF(M.EQ.8)THEN
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NWQOBS
      NSID=IWQOBS(LL,NW)
      L=LIJ( IWQCBS(LL),JWQCBS(LL) )
      LN=LNC(L)
!
      IF(VHDXWQ(LN,K).LT.0.)THEN
       CTMP=CON1(L,K)+DELT*(VHDXWQ(LN,K)*CON1(L,K)-FVHU(LN,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
!      IF(CON(L,K).GT.CBS(LL,1,M)) CON(L,K)=CBS(LL,1,M)
       CWQLOS(LL,K,NW)=CON(L,K)        
       NWQLOS(LL,K,NW)=N
      ELSE
       CBT=WTCI(K,1)*WQOBCS(LL,1,NW)+WTCI(K,2)*WQOBCS(LL,2,NW)+CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOS(LL,K,NW)
       IF(NMNLO.GE.NTSCRS(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOS(LL,K,NW)+(CBT-CWQLOS(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
       ENDIF
      ENDIF
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NWQOBW
      NSID=IWQOBW(LL,NW)
      L=LIJ( IWQCBW(LL),JWQCBW(LL) )
!
      IF(UHDYWQ(L+1,K).LT.0.)THEN
       CTMP=CON1(L,K)+DELT*(UHDYWQ(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
!      IF(CON(L,K).GT.CBW(LL,1,M)) CON(L,K)=CBW(LL,1,M)
       CWQLOW(LL,K,NW)=CON(L,K)
       NWQLOW(LL,K,NW)=N
      ELSE
       CBT=WTCI(K,1)*WQOBCW(LL,1,NW)+WTCI(K,2)*WQOBCW(LL,2,NW)+CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOW(LL,K,NW)
       IF(NMNLO.GE.NTSCRW(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOW(LL,K,NW)+(CBT-CWQLOW(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
       ENDIF
      ENDIF
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NWQOBE
      NSID=IWQOBE(LL,NW)
      L=LIJ( IWQCBE(LL),JWQCBE(LL) )
!
      IF(UHDYWQ(L,K).GT.0.)THEN
       CTMP=CON1(L,K)+DELT*(FUHU(L,K)-UHDYWQ(L,K)*CON1(L,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
!      IF(CON(L,K).GT.CBE(LL,1,M)) CON(L,K)=CBE(LL,1,M)
       CWQLOE(LL,K,NW)=CON(L,K)
       NWQLOE(LL,K,NW)=N
      ELSE
       CBT=WTCI(K,1)*WQOBCE(LL,1,NW)+WTCI(K,2)*WQOBCE(LL,2,NW)+CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOE(LL,K,NW)
       IF(NMNLO.GE.NTSCRE(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOE(LL,K,NW)+(CBT-CWQLOE(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
       ENDIF
      ENDIF
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NWQOBN
      NSID=IWQOBN(LL,NW)
      L=LIJ( IWQCBN(LL),JWQCBN(LL) )
      LS=LSC(L)
!
      IF(VHDXWQ(L,K).GT.0.)THEN
       CTMP=CON1(L,K)+DELT*(FVHU(L,K)-VHDXWQ(L,K)*CON1(L,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
!      IF(CON(L,K).GT.CBN(LL,1,M)) CON(L,K)=CBN(LL,1,M)
       CWQLON(LL,K,NW)=CON(L,K)
       NWQLON(LL,K,NW)=N
      ELSE
       CBT=WTCI(K,1)*WQOBCN(LL,1,NW)+WTCI(K,2)*WQOBCN(LL,2,NW)+CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLON(LL,K,NW)
       IF(NMNLO.GE.NTSCRN(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLON(LL,K,NW)+(CBT-CWQLON(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
       ENDIF
      ENDIF
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      ENDIF
!
!**********************************************************************C
!
! **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION 
!
      IF(ISADAC(MVAR).EQ.0) RETURN
      IF(ISCDCA(MVAR).EQ.1) RETURN
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=1,LC
      UUU(L,K)=0.0
      VVV(L,K)=0.0
	WWW(L,K)=0.0
      ENDDO
      ENDDO
!
      DO L=1,LC
	WWW(L,0)=0.0
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      UUU(L,K)=UWQ(L,K)*(CON(L,K)-CON(L-1,K))*DXIU(L)
      VVV(L,K)=VWQ(L,K)*(CON(L,K)-CON(LS,K))*DYIV(L)
      ENDDO
      ENDDO
!
      DO K=1,KS
      RDZIG=DZIG(K)
      DO L=2,LA
      LS=LSC(L)
      WWW(L,K)=WWQ(L,K)*(CON(L,K+1)-CON(L,K))*RDZIG/HWQ(L)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)          
      LNW=LNWC(L)
      LSE=LSEC(L)      
      AUHU=ABS(UHDYWQ(L,K))
      AVHV=ABS(VHDXWQ(L,K))
!      UTERM=AUHU*(1.-DELTA*AUHU/(DXYU(L)*HU(L)))*(CON(L,K)-CON(L-1,K))
!      VTERM=AVHV*(1.-DELTA*AVHV/(DXYV(L)*HV(L)))*(CON(L,K)-CON(LS,K))
      UTERM=AUHU*(CON(L,K)-CON(L-1,K))
      VTERM=AVHV*(CON(L,K)-CON(LS,K))
      IF(ISADAC(MVAR).GE.2)THEN
        SSCORUE=DELTA*RDZIC*DXYIP(L)*HWQI(L)*(FQCPAD(L,K)-QSUMPAD(L,K)*CON(L,K))
        SSCORUW=DELTA*RDZIC*DXYIP(L-1)*HWQI(L-1)*(FQCPAD(L-1,K)-QSUMPAD(L-1,K)*CON(L-1,K))
        SSCORVN=DELTA*RDZIC*DXYIP(L)*HWQI(L)*(FQCPAD(L,K)-QSUMPAD(L,K)*CON(L,K))
        SSCORVS=DELTA*RDZIC*DXYIP(LS)*HWQI(LS)*(FQCPAD(LS,K)-QSUMPAD(LS,K)*CON(LS,K))
        SSCORU=MAX(UHDYWQ(L,K),0.0)*SSCORUW+MIN(UHDYWQ(L,K),0.0)*SSCORUE
        SSCORV=MAX(VHDXWQ(L,K),0.0)*SSCORVS+MIN(VHDXWQ(L,K),0.0)*SSCORVN
	  UTERM=UTERM+SSCORU
	  VTERM=VTERM+SSCORV
	ENDIF
!     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
!     AHU(L,K)=MAX(AHTMP,0.)
!     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
!     AHV(L,K)=MAX(AHTMP,0.)
!     UTERM=UTERM*(CON(L,K)-CON(L-1,K))
!     VTERM=VTERM*(CON(L,K)-CON(LS,K))
!
!ADD1031      UTERM=UTERM-0.25*DELTA*UHDYWQ(L,K)*
!ADD1031     &      (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K)
!ADD1031     &      +WWW(L,K)+WWW(L-1,K)+WWW(L-1,K-1)+WWW(L,K-1))
!ADD1031     VTERM=VTERM-0.25*DELTA*VHDXWQ(L,K)*
!ADD1031     &      (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K)
!ADD1031     &      +WWW(L,K)+WWW(LS,K)+WWW(LS,K-1)+WWW(L,K-1))
!
!
      IF(UHDYWQ(L,K).GE.0.0)THEN
        UTERM=UTERM-0.5*DELTA*UHDYWQ(L,K)*(VVV(LNW,K)+VVV(L-1,K)+WWW(L-1,K)+WWW(L-1,K-1)          &
              +UUU(L,K)+UUU(L-1,K))
      ELSE
        UTERM=UTERM-0.5*DELTA*UHDYWQ(L,K)*(VVV(LN,K)+VVV(L,K)+WWW(L,K)+WWW(L,K-1)+UUU(L,K)+UUU(L+1,K))
	ENDIF
!
      IF(VHDXWQ(L,K).GE.0.0)THEN
        VTERM=VTERM-0.5*DELTA*VHDXWQ(L,K)*(UUU(LS,K)+UUU(LSE,K)+WWW(LS,K)+WWW(LS,K-1)+VVV(LS,K)   &
              +VVV(L,K))
      ELSE
        VTERM=VTERM-0.5*DELTA*VHDXWQ(L,K)*(UUU(L,K)+UUU(L+1,K)+WWW(L,K)+WWW(L,K-1)+VVV(LN,K)+VVV(L,K))
	ENDIF
!
!
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
      ENDIF
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(WWQ(L,K))
!      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)/HWQ(L))*(CON(L,K+1)-CON(L,K))
      WTERM=AWW*(CON(L,K+1)-CON(L,K))
!
      IF(ISADAC(MVAR).GE.2)THEN
        SSCORWA=DELTA*DZIG(K+1)*HWQI(L)*DXYIP(L)*(FQCPAD(L,K+1)-QSUMPAD(L,K+1)*CON(L,K+1))
        SSCORWB=DELTA*DZIG(K)*HWQI(L)*DXYIP(L)*(FQCPAD(L,K)-QSUMPAD(L,K)*CON(L,K))
        !SSCORWB=DELTA*DZIG(K)*HWQI(L)*DXYIP(L)*(FQCPAD(L,K)-QSUMPAD(L,K)*CON(L,K)
        SSCORW=MAX(WWQ(L,K),0.0)*SSCORWB+MIN(WWQ(L,K),0.0)*SSCORWA
	  WTERM=WTERM+SSCORW
      ENDIF
!
!CC1031      WTERM=WTERM-0.25*DELTA*WWQ(L,K)*
!CC1031     &      (UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)
!CC1031     &      +VVV(L,K)+VVV(LN,K)+VVV(LN,K+1)+VVV(L,K+1))
!
      IF(WWQ(L,K).GE.0.0)THEN
        WTERM=WTERM-0.5*DELTA*WWQ(L,K)*(UUU(L,K)+UUU(L+1,K)+VVV(L,K)+VVV(LN,K)+WWW(L,K)+WWW(L,K-1))
	ELSE
        WTERM=WTERM-0.5*DELTA*WWQ(L,K)*(UUU(L+1,K+1)+UUU(L,K+1)+VVV(LN,K+1)+VVV(L,K+1)+WWW(L,K)   &
              +WWW(L,K+1))
	ENDIF
!
!
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
      FVHU(LN,K)=0.0
      ENDDO
!
      DO LL=1,NCBW
      L=LCBW(LL)
      FUHU(L+1,K)=0.0
      ENDDO
!
      DO LL=1,NCBE
      L=LCBE(LL)
      FUHU(L,K)=0.0
      ENDDO
!
      DO LL=1,NCBN
      L=LCBN(LL)
      FVHU(L,K)=0.0
      ENDDO
!
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE AND APPLY FLUX CORRECTED TRANSPORT LIMITERS
!
      IF(ISFCT(MVAR).EQ.0) GOTO 600
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
!
!1031      DO L=2,LA
!1031      CMAX(L,1)=MAX(CONT(L,1),CONT(L,2))
!1031      CMAX(L,KC)=MAX(CONT(L,KS),CONT(L,KC))
!1031      ENDDO
!
!1031      DO K=2,KS
!1031      DO L=2,LA
!1031      CMAX(L,K)=MAX(CONT(L,K-1),CONT(L,K),CONT(L,K+1))
!1031      ENDDO
!1031      ENDDO
!
!1031      DO K=1,KC
!1031      DO L=2,LA
!1031      LS=LSC(L)
!1031      LN=LNC(L)
!1031      CWMAX=SUB(L)*CONT(L-1,K) 
!1031      CEMAX=SUB(L+1)*CONT(L+1,K)
!1031      CSMAX=SVB(L)*CONT(LS,K)
!1031      CNMAX=SVB(LN)*CONT(LN,K)
!1031      CMAX(L,K)=MAX(CMAX(L,K),CNMAX,CEMAX,CSMAX,CWMAX)
!1031      ENDDO
!1031      ENDDO
!
!1031      DO K=1,KC
!1031      DO L=2,LA
!1031      CONT(L,K)=MIN(CON(L,K),CON2(L,K))
!1031      ENDDO
!1031      ENDDO
!
!1031      DO L=2,LA
!1031      CMIN(L,1)=MIN(CONT(L,1),CONT(L,2))
!1031      CMIN(L,KC)=MIN(CONT(L,KS),CONT(L,KC))
!1031      ENDDO
!
!1031      DO K=2,KS
!1031      DO L=2,LA
!1031      CMIN(L,K)=MIN(CONT(L,K-1),CONT(L,K),CONT(L,K+1))
!1031      ENDDO
!1031      ENDDO
!
!1031      DO K=1,KC
!1031      DO L=2,LA
!1031      LS=LSC(L)
!1031      LN=LNC(L)
!1031      CWMIN=SUB(L)*CONT(L-1,K)+1.E+6*(1.-SUB(L))
!1031      CEMIN=SUB(L+1)*CONT(L+1,K)+1.E+6*(1.-SUB(L+1))
!1031      CSMIN=SVB(L)*CONT(LS,K)+1.E+6*(1.-SVB(L))
!1031      CNMIN=SVB(LN)*CONT(LN,K)+1.E+6*(1.-SVB(LN))
!1031      CMIN(L,K)=MIN(CMIN(L,K),CNMIN,CEMIN,CSMIN,CWMIN)
!1031      ENDDO
!1031      ENDDO
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
      DU(L,K)=DELT*SCB(L)*(DXYIP(L)*(FUHU(L,K)-FUHV(L+1,K)+FVHU(L,K)-FVHV(LN,K))+RDZIC*(FWU(L,K-1) &
              -FWV(L,K)))*HWQI(L)
      DV(L,K)=DELT*SCB(L)*(DXYIP(L)*(FUHU(L+1,K)-FUHV(L,K)+FVHU(LN,K)-FVHV(L,K))+RDZIC*(FWU(L,K)   &
              -FWV(L,K-1)))*HWQI(L)
      ENDDO
      ENDDO
!
      DO K=1,KC
!
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      DU(LN,K)=0.
      DV(LN,K)=0.
      ENDDO
!
      DO LL=1,NCBW
      L=LCBW(LL)
      DU(L+1,K)=0.
      DV(L+1,K)=0.
      ENDDO
!
      DO LL=1,NCBE
      L=LCBE(LL)
      DU(L-1,K)=0.
      DV(L-1,K)=0.
      ENDDO
!
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      DU(LS,K)=0.
      DV(LS,K)=0.
      ENDDO

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
  600 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      CH(L,K)=CON(L,K)*HWQ(L)+DELT*((FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K))*DXYIP(L)        &
              +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)/HWQ(L)+(1.-SCB(L))*CON(L,K)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  DIAGNOSE FCT SCHEME
!
      IF(ISFCT(M).EQ.99)THEN
      WRITE(6,6010)N
!
      DO K=1,KC
      DO L=2,LA
      CCMAX=SCB(L)*(CON(L,K)-CMAX(L,K))
      IF(CCMAX.GT.0.)THEN
       WRITE(6,6011)CON(L,K),CMAX(L,K),IL(L),JL(L),K
      ENDIF
      CCMIN=SCB(L)*(CMIN(L,K)-CON(L,K))
      IF(CCMIN.GT.0.)THEN
       WRITE(6,6012)CMIN(L,K),CON(L,K),IL(L),JL(L),K
      ENDIF
      ENDDO
      ENDDO
!
      ENDIF
!
 6010 FORMAT('  FCT DIAGNOSTICS AT N = ',I5)
 6011 FORMAT('  CON = ',E12.4,3X,'CMAX = ',E12.4,3X,'I,J,K=',(3I10))
 6012 FORMAT('  CMIN = ',E12.4,3X,'CON = ',E12.4,3X,'I,J,K=',(3I10))
!
!----------------------------------------------------------------------C
!
      RETURN
      END
