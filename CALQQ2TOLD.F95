!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALQQ2TOLD (ISTL)
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 01/12/2002        John Hamrick       01/12/2002       John Hamrick
!  fixed dynamic time stepping
!----------------------------------------------------------------------C
!
! **  SUBROUTINE CALQQ CALCULATES THE TURBULENT INTENSITY SQUARED AT 
! **  TIME LEVEL (N+1).  THE VALUE OF ISTL INDICATES THE NUMBER OF 
! **  TIME LEVELS INVOLVED
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!DIAG      DIMENSION UUUTMP(LCM),EQTMP(LCM)
!
!**********************************************************************C
!
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
        DELTD2=0.5*DT
        DELTI=1./DELT
      ELSE
        DELT=DTDYN
        DELTD2=0.5*DTDYN
        DELTI=1./DELT
      END IF   
!
      S3TL=1.0
      S2TL=0.0
      BSMALL=1.E-12
!
      DO K=1,KS
        DO L=2,LA
          QQ2(L,K)=QQ(L,K)+QQ(L,K)
          QQL2(L,K)=QQL(L,K)+QQL(L,K)
        ENDDO
      ENDDO
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
        FWQQ(L,K)=0.
        FWQQL(L,K)=0.
       ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH TRANSPORT
! **  AVERAGED BETWEEN (N) AND (N+1) AND TRANSPORTED FIELD AT (N) OR
! **  TRANSPORT BETWEEN (N-1) AND (N+1) AND TRANSPORTED FIELD AT (N-1)
! **  FOR ISTL EQUAL TO 2 AND 3 RESPECTIVELY 
!
! **  FUHQQ=FUHU, FVHQQ=FVHU, FUHQQL=FUHV, FVHQQL=FVHV
!
!----------------------------------------------------------------------C
!
!      IF(ISTL.EQ.2)THEN
!
      IF(IDRYTBP.EQ.0)THEN
!
      DO K=1,KC
      DO L=2,LA
        WB=0.5*DXYP(L)*(W2(L,K-1)+W2(L,K))
        FWQQ(L,K)=MAX(WB,0.)*QQ(L,K-1)+MIN(WB,0.)*QQ(L,K)
        FWQQL(L,K)=MAX(WB,0.)*QQL(L,K-1)+MIN(WB,0.)*QQL(L,K)
      ENDDO
      ENDDO
!
      ELSE
!
      DO K=1,KC
      DO L=2,LA
      IF(LMASKDRY(L))THEN
        WB=0.5*DXYP(L)*(W2(L,K-1)+W2(L,K))
        FWQQ(L,K)=MAX(WB,0.)*QQ(L,K-1)+MIN(WB,0.)*QQ(L,K)
        FWQQL(L,K)=MAX(WB,0.)*QQL(L,K-1)+MIN(WB,0.)*QQL(L,K)
      ENDIF
      ENDDO
      ENDDO
!
      ENDIF
!
!
!      ELSE
!
!      IF(ISCDCA(0).EQ.1)THEN
!
!      DO K=1,KC
!      DO L=2,LA
!        WB=0.25*DXYP(L)*(W2(L,K-1)+W2(L,K))
!        FWQQ(L,K)=WB*(QQ(L,K-1)+QQ(L,K))
!        FWQQL(L,K)=WB*(QQL(L,K-1)+QQL(L,K))
!      ENDDO
!      ENDDO
!
!      ELSE
!
!      DO K=1,KC
!      DO L=2,LA
!        WB=0.25*DXYP(L)*(W2(L,K-1)+W2(L,K))
!        FWQQ(L,K)=MAX(WB,0.)*QQ2(L,K-1)
!     &         +MIN(WB,0.)*QQ2(L,K)
!        FWQQL(L,K)=MAX(WB,0.)*QQL2(L,K-1)
!     &          +MIN(WB,0.)*QQL2(L,K)
!      ENDDO
!      ENDDO
!
!      ENDIF
!
!      ENDIF
!
!----------------------------------------------------------------------C
!
!      IF(ISTL.EQ.2)THEN
!
!
      IF(IDRYTBP.EQ.0)THEN
!
      DO K=1,KS
      DO L=2,LA
        LS=LSC(L)
        UHUW=0.5*(UHDY2(L,K)+UHDY2(L,K+1))
        VHVW=0.5*(VHDX2(L,K)+VHDX2(L,K+1))
        FUHU(L,K)=MAX(UHUW,0.)*QQ(L-1,K)+MIN(UHUW,0.)*QQ(L,K)
        FVHU(L,K)=MAX(VHVW,0.)*QQ(LS,K)+MIN(VHVW,0.)*QQ(L,K)
        FUHV(L,K)=MAX(UHUW,0.)*QQL(L-1,K)+MIN(UHUW,0.)*QQL(L,K)
        FVHV(L,K)=MAX(VHVW,0.)*QQL(LS,K)+MIN(VHVW,0.)*QQL(L,K)
      ENDDO
      ENDDO
!
      ELSE
!
      DO K=1,KS
      DO L=2,LA
      IF(LMASKDRY(L))THEN
        LS=LSC(L)
        UHUW=0.5*(UHDY2(L,K)+UHDY2(L,K+1))
        VHVW=0.5*(VHDX2(L,K)+VHDX2(L,K+1))
        FUHU(L,K)=MAX(UHUW,0.)*QQ(L-1,K)+MIN(UHUW,0.)*QQ(L,K)
        FVHU(L,K)=MAX(VHVW,0.)*QQ(LS,K)+MIN(VHVW,0.)*QQ(L,K)
        FUHV(L,K)=MAX(UHUW,0.)*QQL(L-1,K)+MIN(UHUW,0.)*QQL(L,K)
        FVHV(L,K)=MAX(VHVW,0.)*QQL(LS,K)+MIN(VHVW,0.)*QQL(L,K)
      ENDIF
      ENDDO
      ENDDO
!
	ENDIF
!
!      ELSE
!
!      IF(ISCDCA(0).EQ.1)THEN
!
!      DO K=1,KS
!      DO L=2,LA
!        LS=LSC(L)
!        UHUW=0.25*(UHDY2(L,K)+UHDY2(L,K+1))
!        VHVW=0.25*(VHDX2(L,K)+VHDX2(L,K+1))
!        FUHU(L,K)=UHUW*(QQ(L-1,K)+QQ(L,K))
!        FVHU(L,K)=VHVW*(QQ(LS,K)+QQ(L,K))
!        FUHV(L,K)=UHUW*(QQL(L-1,K)+QQL(L,K))
!        FVHV(L,K)=VHVW*(QQL(LS,K)+QQL(L,K))
!      ENDDO
!      ENDDO
!
!      ELSE
!
!      DO K=1,KS
!      DO L=2,LA
!        LS=LSC(L)
!        UHUW=0.25*(UHDY2(L,K)+UHDY2(L,K+1))
!        VHVW=0.25*(VHDX2(L,K)+VHDX2(L,K+1))
!        FUHU(L,K)=MAX(UHUW,0.)*QQ2(L-1,K)
!     &         +MIN(UHUW,0.)*QQ2(L,K)
!        FVHU(L,K)=MAX(VHVW,0.)*QQ2(LS,K)
!     &         +MIN(VHVW,0.)*QQ2(L,K)
!        FUHV(L,K)=MAX(UHUW,0.)*QQL2(L-1,K)
!     &         +MIN(UHUW,0.)*QQL2(L,K)
!        FVHV(L,K)=MAX(VHVW,0.)*QQL2(LS,K)
!     &         +MIN(VHVW,0.)*QQL2(L,K)
!      ENDDO
!      ENDDO
!
!      ENDIF
!
!      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE PRODUCTION, LOAD BOUNDARY CONDITIONS AND SOLVE 
! **  TRANSPORT EQUATIONS
!
! **  FUHQQ=FUHU, FVHQQ=FVHU, FUHQQL=FUHV, FVHQQL=FVHV
!
! **  CU1=CUQ, CU2=CUQL, UUU=QQH, VVV=QQLH
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
!
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      IF(FVHU(LN,K).GT.0)THEN
      FVHU(LN,K)=0.0
      FVHV(LN,K)=0.0
      ENDIF
      ENDDO
!
      DO LL=1,NCBW
      L=LCBW(LL)
      IF(FUHU(L+1,K).GT.0)THEN
      FUHU(L+1,K)=0.0
      FUHV(L+1,K)=0.0
      ENDIF
      ENDDO
!
      DO LL=1,NCBE
      L=LCBE(LL)
      IF(FUHU(L,K).LT.0.)THEN
      FUHU(L,K)=0.0
      FUHV(L,K)=0.0
      ENDIF
      ENDDO
!
      DO LL=1,NCBN
      L=LCBN(LL)
      IF(FVHU(L,K).LT.0.)THEN
      FVHU(L,K)=0.0
      FVHV(L,K)=0.0
      ENDIF
      ENDDO
!
      ENDDO
!
!----------------------------------------------------------------------C
!
      IF(ISWAVE.LE.1)THEN
!
!      IF(ISTL.EQ.2)THEN
!
!
      IF(IDRYTBP.EQ.0)THEN
!
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      UUU(L,K)=QQ(L,K)*H1P(L)+DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)+(FWQQ(L,K)         &
               -FWQQ(L,K+1))*DZIG(K))*DXYIP(L)
!      VVV(L,K)=QQL1(L,K)*H1P(L)
      VVV(L,K)=QQL(L,K)*H1P(L)+DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)+(FWQQL(L,K)       &
               -FWQQL(L,K+1))*DZIG(K))*DXYIP(L)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      UUU(L,K)=MAX(UUU(L,K),0.)
      VVV(L,K)=MAX(VVV(L,K),0.)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PQQB=AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
      PQQU=AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2+AV(L,K)*DZIGSD4(K)         &
           *(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 
!JMH      RFNUM=-PQQB/PQQU
      PQQ=DELT*(PQQB+PQQU)
      UUU(L,K)=UUU(L,K)+2.*PQQ
      VVV(L,K)=VVV(L,K)+CTE1*DML(L,K)*PQQ
      ENDDO
      ENDDO
!
      ELSE
!
!
      DO K=1,KS
      DO L=2,LA
      IF(LMASKDRY(L))THEN
      LN=LNC(L)
      UUU(L,K)=QQ(L,K)*H1P(L)+DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)+(FWQQ(L,K)         &
               -FWQQ(L,K+1))*DZIG(K))*DXYIP(L)
!      VVV(L,K)=QQL1(L,K)*H1P(L)
      VVV(L,K)=QQL(L,K)*H1P(L)+DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)+(FWQQL(L,K)       &
               -FWQQL(L,K+1))*DZIG(K))*DXYIP(L)
      ENDIF
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      IF(LMASKDRY(L))THEN
      UUU(L,K)=MAX(UUU(L,K),0.)
      VVV(L,K)=MAX(VVV(L,K),0.)
      ENDIF
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      IF(LMASKDRY(L))THEN
      LN=LNC(L)
      LS=LSC(L)
      PQQB=AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
      PQQU=AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2+AV(L,K)*DZIGSD4(K)         &
           *(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 
!JMH      RFNUM=-PQQB/PQQU
      PQQ=DELT*(PQQB+PQQU)
      UUU(L,K)=UUU(L,K)+2.*PQQ
      VVV(L,K)=VVV(L,K)+CTE1*DML(L,K)*PQQ
      ENDIF
      ENDDO
      ENDDO
!
      ENDIF
!
!      ELSE
!
!      DO K=1,KS
!      DO L=2,LA
!      LN=LNC(L)
!      UUU(L,K)=QQ1(L,K)*H2P(L)
!     &      +DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)
!     &      +(FWQQ(L,K)-FWQQ(L,K+1))*DZIG(K))*DXYIP(L)
!      VVV(L,K)=QQL1(L,K)*H2P(L)
!     &      +DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)
!     &      +(FWQQL(L,K)-FWQQL(L,K+1))*DZIG(K))*DXYIP(L)
!      ENDDO
!      ENDDO
!
!      DO K=1,KS
!      DO L=2,LA
!      UUU(L,K)=MAX(UUU(L,K),0.)
!      VVV(L,K)=MAX(VVV(L,K),0.)
!      ENDDO
!      ENDDO
!
!      DO K=1,KS
!      DO L=2,LA
!      LN=LNC(L)
!      PQQB=AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
!      PQQU=AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
!     &    +AV(L,K)*DZIGSD4(K)*(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2
!CJMH      RFNUM=-PQQB/PQQU
!      PQQ=DELT*(PQQB+PQQU)
!      UUU(L,K)=UUU(L,K)+2.*PQQ
!      VVV(L,K)=VVV(L,K)+CTE1*DML(L,K)*PQQ
!      ENDDO
!      ENDDO
!
!      ENDIF
      ENDIF
!
!----------------------------------------------------------------------C
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
!      IF(ISTL.EQ.2)THEN
!
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         TVAR1W(L,K)=(WVDTKEM(K)*WVDISP(L,K)+WVDTKEP(K)*WVDISP(L,K+1))
        ENDDO
       ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PQQ=DELT*( AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))+AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)        &
          -U(L+1,K)+U(L,K+1)-U(L,K))**2+AV(L,K)*DZIGSD4(K)*(V(LN,K+1)-V(LN,K)+V(L,K+1)            &
          -V(L,K))**2+WVFACT*TVAR1W(L,K) )
      UUU(L,K)=QQ1(L,K)*H1P(L)+DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)+(FWQQ(L,K)        &
               -FWQQ(L,K+1))*DZIG(K))*DXYIP(L)+2.*PQQ
      VVV(L,K)=QQL1(L,K)*H1P(L)+DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)+(FWQQL(L,K)      &
               -FWQQL(L,K+1))*DZIG(K))*DXYIP(L)+CTE1*DML(L,K)*PQQ
      ENDDO
      ENDDO
!
!      ELSE
!
!      DO ND=1,NDM
!       LF=2+(ND-1)*LDM
!       LL=LF+LDM-1
!       DO K=1,KS
!        DO L=LF,LL
!         TVAR1W(L,K)=(WVDTKEM(K)*WVDISP(L,K)
!     &               +WVDTKEP(K)*WVDISP(L,K+1))
!        ENDDO
!       ENDDO
!      ENDDO
!
!      DO K=1,KS
!      DO L=2,LA
!      LN=LNC(L)
!      LS=LSC(L)
!      PQQ=DELT*(AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
!     &    +AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
!     &    +AV(L,K)*DZIGSD4(K)*(V(LN,K+1)-V(LN,K)+V(L,K+1)-V(L,K))**2
!     &    +WVFACT*TVAR1W(L,K) )
!      UUU(L,K)=QQ1(L,K)*H2P(L)
!     &      +DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)
!     &      +(FWQQ(L,K)-FWQQ(L,K+1))*DZIG(K))*DXYIP(L)+2.*PQQ
!      VVV(L,K)=QQL1(L,K)*H2P(L)
!     &      +DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)
!     &      +(FWQQL(L,K)-FWQQL(L,K+1))*DZIG(K))*DXYIP(L)
!     &      +CTE1*DML(L,K)*PQQ
!      ENDDO
!      ENDDO
!
!      ENDIF
!
      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(KC.LE.2)THEN
!
!
      IF(IDRYTBP.EQ.0)THEN
!
      DO L=2,LA
        CLQTMP=-DELT*CDZKK(1)*AQ(L,1)*HPI(L)
        CUQTMP=-DELT*CDZKKP(1)*AQ(L,2)*HPI(L)
        CMQTMP=1.-CLQTMP-CUQTMP+2.*DELT*SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L))
        CMQLTMP=1.-CLQTMP-CUQTMP+DELT*(SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L)))*(1.+CTE4      &
                *DML(L,1)*DML(L,1)*FPROX(1))
        EQ=1./CMQTMP
        EQL=1./CMQLTMP
        CU1(L,1)=CUQTMP*EQ
        CU2(L,1)=CUQTMP*EQL
        UUU(L,1)=(UUU(L,1)-CLQTMP*HP(L)*QQ(L,0)-CUQTMP*HP(L)*QQ(L,KC))*EQ
        VVV(L,1)=VVV(L,1)*EQL
      ENDDO
!
      ELSE
!
      DO L=2,LA
	IF(LMASKDRY(L))THEN
        CLQTMP=-DELT*CDZKK(1)*AQ(L,1)*HPI(L)
        CUQTMP=-DELT*CDZKKP(1)*AQ(L,2)*HPI(L)
        CMQTMP=1.-CLQTMP-CUQTMP+2.*DELT*SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L))
        CMQLTMP=1.-CLQTMP-CUQTMP+DELT*(SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L)))*(1.+CTE4      &
                *DML(L,1)*DML(L,1)*FPROX(1))
        EQ=1./CMQTMP
        EQL=1./CMQLTMP
        CU1(L,1)=CUQTMP*EQ
        CU2(L,1)=CUQTMP*EQL
        UUU(L,1)=(UUU(L,1)-CLQTMP*HP(L)*QQ(L,0)-CUQTMP*HP(L)*QQ(L,KC))*EQ
        VVV(L,1)=VVV(L,1)*EQL
	ENDIF
      ENDDO

!
      ENDIF

!
!DIAG      IF(N.EQ.300)THEN
!DIAG        OPEN(1,FILE='TURB.DIA',STATUS='UNKNOWN')
!DIAG        CLOSE(1,STATUS='DELETE')
!DIAG        OPEN(1,FILE='TURB.DIA',STATUS='UNKNOWN')
!DIAG        WRITE(1,110)
!DIAG        DO L=2,LA
!DIAG         QQTMP=UUU(L,1)*HPI(L)
!DIAG         WRITE(1,111) IL(L),JL(L),QQ(L,0),QQTMP,
!DIAG     &                QQ(L,KC),UUUTMP(L),EQTMP(L)
!DIAG        ENDDO
!DIAG        CLOSE(1)
!DIAG      ENDIF
!
      ENDIF
!
      IF(KC.GT.2)THEN
!
!
      IF(IDRYTBP.EQ.0)THEN
!

      DO L=2,LA
	  CLQTMP=-DELT*CDZKK(1)*AQ(L,1)*HPI(L)
        CUQTMP=-DELT*CDZKKP(1)*AQ(L,2)*HPI(L)
        CMQTMP=1.-CLQTMP-CUQTMP+2.*DELT*SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L))
        CMQLTMP=1.-CLQTMP-CUQTMP+DELT*(SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L)))*(1.+CTE4      &
                *DML(L,1)*DML(L,1)*FPROX(1))
        EQ=1./CMQTMP
        EQL=1./CMQLTMP
        CU1(L,1)=CUQTMP*EQ
        CU2(L,1)=CUQTMP*EQL
        UUU(L,1)=(UUU(L,1)-CLQTMP*HP(L)*QQ(L,0))*EQ
        VVV(L,1)=VVV(L,1)*EQL
        CUQTMP=-DELT*CDZKKP(KS)*AQ(L,KC)*HPI(L)
        UUU(L,KS)=UUU(L,KS)-CUQTMP*HP(L)*QQ(L,KC)
      ENDDO
!
      DO K=2,KS
        DO L=2,LA
         CLQTMP=-DELT*CDZKK(K)*AQ(L,K)*HPI(L)
         CUQTMP=-DELT*CDZKKP(K)*AQ(L,K+1)*HPI(L)
         CMQTMP=1.-CLQTMP-CUQTMP+2.*DELT*SQRT(QQ(L,K))/(CTURBB1(L,K)*DML(L,K)*HP(L))
         CMQLTMP=1.-CLQTMP-CUQTMP+DELT*(SQRT(QQ(L,K))/(CTURBB1(L,K)*DML(L,K)*HP(L)))*(1.+CTE4     &
                *DML(L,K)*DML(L,K)*FPROX(K))
         EQ=1./(CMQTMP-CLQTMP*CU1(L,K-1))
         EQL=1./(CMQLTMP-CLQTMP*CU2(L,K-1))
         CU1(L,K)=CUQTMP*EQ
         CU2(L,K)=CUQTMP*EQL     
!        IF(EQ.0) PAUSE
         UUU(L,K)=(UUU(L,K)-CLQTMP*UUU(L,K-1))*EQ
         VVV(L,K)=(VVV(L,K)-CLQTMP*VVV(L,K-1))*EQL
        ENDDO
      ENDDO
!
      DO K=KS-1,1,-1
        DO L=2,LA
         UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)
         VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)
        ENDDO
      ENDDO
!
      ELSE
!
      DO L=2,LA
	IF(LMASKDRY(L))THEN
	  CLQTMP=-DELT*CDZKK(1)*AQ(L,1)*HPI(L)
        CUQTMP=-DELT*CDZKKP(1)*AQ(L,2)*HPI(L)
        CMQTMP=1.-CLQTMP-CUQTMP+2.*DELT*SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L))
        CMQLTMP=1.-CLQTMP-CUQTMP+DELT*(SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L)))*(1.+CTE4      &
                *DML(L,1)*DML(L,1)*FPROX(1))
        EQ=1./CMQTMP
        EQL=1./CMQLTMP
        CU1(L,1)=CUQTMP*EQ
        CU2(L,1)=CUQTMP*EQL
        UUU(L,1)=(UUU(L,1)-CLQTMP*HP(L)*QQ(L,0))*EQ
        VVV(L,1)=VVV(L,1)*EQL
        CUQTMP=-DELT*CDZKKP(KS)*AQ(L,KC)*HPI(L)
        UUU(L,KS)=UUU(L,KS)-CUQTMP*HP(L)*QQ(L,KC)
        ENDIF
      ENDDO
!
      DO K=2,KS
        DO L=2,LA
	   IF(LMASKDRY(L))THEN
         CLQTMP=-DELT*CDZKK(K)*AQ(L,K)*HPI(L)
         CUQTMP=-DELT*CDZKKP(K)*AQ(L,K+1)*HPI(L)
         CMQTMP=1.-CLQTMP-CUQTMP+2.*DELT*SQRT(QQ(L,K))/(CTURBB1(L,K)*DML(L,K)*HP(L))
         CMQLTMP=1.-CLQTMP-CUQTMP+DELT*(SQRT(QQ(L,K))/(CTURBB1(L,K)*DML(L,K)*HP(L)))*(1.+CTE4     &
                 *DML(L,K)*DML(L,K)*FPROX(K))
         EQ=1./(CMQTMP-CLQTMP*CU1(L,K-1))
         EQL=1./(CMQLTMP-CLQTMP*CU2(L,K-1))
         CU1(L,K)=CUQTMP*EQ
         CU2(L,K)=CUQTMP*EQL     
!        IF(EQ.0) PAUSE
         UUU(L,K)=(UUU(L,K)-CLQTMP*UUU(L,K-1))*EQ
         VVV(L,K)=(VVV(L,K)-CLQTMP*VVV(L,K-1))*EQL
         ENDIF
        ENDDO
      ENDDO
!
      DO K=KS-1,1,-1
        DO L=2,LA
        IF(LMASKDRY(L))THEN
         UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)
         VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)
      ENDIF
        ENDDO
      ENDDO
!
!
      ENDIF
!
      ENDIF
!
!----------------------------------------------------------------------C
!
!
      IF(IDRYTBP.EQ.0)THEN
!

      DO K=1,KS
      DO L=2,LA
      QQ1(L,K)=QQ(L,K)
      QQHDH=UUU(L,K)*HPI(L)
      QQ(L,K)=MAX(QQHDH,QQMIN)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      QQL1(L,K)=QQL(L,K)
      QQHDH=VVV(L,K)*HPI(L)
      QQL(L,K)=MAX(QQHDH,QQLMIN)
      DMLTMP=QQL(L,K)/QQ(L,K)
      DMLTMP=MAX(DMLTMP,DMLMIN)
      DELB=B(L,K)-B(L,K+1)
      IF(DELB.GT.0.0.AND.ISLLIM.EQ.2)THEN
       DMLMAX=SQRT(RIQMAX)*SQRT(QQ(L,K)/(G*HP(L)*DZIG(K)*DELB))
       DML(L,K)=MIN(DMLMAX,DMLTMP)
	 QQL(L,K)=QQ(L,K)*DML(L,K)
      ELSE
       DML(L,K)=DMLTMP
      ENDIF
      ENDDO
      ENDDO
!
      ELSE
!

      DO K=1,KS
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      QQ1(L,K)=QQ(L,K)
      QQHDH=UUU(L,K)*HPI(L)
      QQ(L,K)=MAX(QQHDH,QQMIN)
      ENDIF
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      QQL1(L,K)=QQL(L,K)
      QQHDH=VVV(L,K)*HPI(L)
      QQL(L,K)=MAX(QQHDH,QQLMIN)
      DMLTMP=QQL(L,K)/QQ(L,K)
      DMLTMP=MAX(DMLTMP,DMLMIN)
      DELB=B(L,K)-B(L,K+1)
      IF(DELB.GT.0.0.AND.ISLLIM.EQ.2)THEN
       DMLMAX=SQRT(RIQMAX)*SQRT(QQ(L,K)/(G*HP(L)*DZIG(K)*DELB))
       DML(L,K)=MIN(DMLMAX,DMLTMP)
	 QQL(L,K)=QQ(L,K)*DML(L,K)
      ELSE
       DML(L,K)=DMLTMP
      ENDIF
      ENDIF
      ENDDO
      ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------C
!----------------------------------------------------------------------C
!
!      QQMXSV=-1.E+12
!      QQMNSV=1.E+12
!      QQLMXSV=-1.E+12
!      QQLMNSV=1.E+12
!      
!      DO K=1,KS
!       DO L=2,LA
!        QQ1(L,K)=S2TL*QQ1(L,K)+S3TL*QQ(L,K)
!        QQL1(L,K)=S2TL*QQL1(L,K)+S3TL*QQL(L,K)
!        QQTMP=UUU(L,K)*HPI(L)
!        QQLTMP=VVV(L,K)*HPI(L)
!      QQMXSV=MAX(QQTMP,QQMXSV)
!      QQMNSV=MIN(QQTMP,QQMNSV)
!      IF(QQMNSV.LT.0.)THEN
!        WRITE(8,601)IL(L),JL(L),K,QQMNSV
!      ENDIF
!      QQLMXSV=MAX(QQLTMP,QQLMXSV)
!      QQLMNSV=MIN(QQLTMP,QQLMNSV)
!        QQ(L,K)=MAX(QQTMP,QQMIN)
!        QQL(L,K)=MAX(QQLTMP,QQLMIN)
!        DMLTMP=QQL(L,K)/QQ(L,K)
!        DML(L,K)=MAX(DMLTMP,DMLMIN)
!       ENDDO
!      ENDDO
!
!      IF(ISTOPT(0).EQ.1)THEN
!        DO K=1,KS
!        DO L=2,LA
!          DELB=B(L,K)-B(L,K+1)
!          IF(DELB.GT.0.0)THEN
!            DMLMAX=CTE3*SQRT(QQ(L,K)/(G*HP(L)*DZIG(K)*DELB))
!            DML(L,K)=MIN(DMLMAX,DML(L,K))
!          ENDIF
!        ENDDO
!        ENDDO
!      ENDIF
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      QQ(L,K)=QQ(LN,K)
      QQL(L,K)=QQL(LN,K)
      DML(L,K)=DML(LN,K)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
      DO LL=1,NCBW
      L=LCBW(LL)
      QQ(L,K)=QQ(L+1,K)
      QQL(L,K)=QQL(L+1,K)
      DML(L,K)=DML(L+1,K)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
      DO LL=1,NCBE
      L=LCBE(LL)
      QQ(L,K)=QQ(L-1,K)
      QQL(L,K)=QQL(L-1,K)
      DML(L,K)=DML(L-1,K)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KS
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      QQ(L,K)=QQ(LS,K)
      QQL(L,K)=QQL(LS,K)
      DML(L,K)=DML(LS,K)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
  110 FORMAT('    I    J   QQ BOT        QQ MID        QQ SURF','       PROD+ADV      1./DIAGON')
  111 FORMAT(2I5,5E14.5)
!
  600 FORMAT('N,QX,QN,QLX,QLN,CX,CN=',I5,6E12.4)
  601 FORMAT('NEG QQ I,J,K,QQ=',3I5,E13.5)
!
!----------------------------------------------------------------------C
!
      RETURN
      END