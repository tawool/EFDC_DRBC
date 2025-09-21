!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CONGRADC (ISTL)
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
! **  SUBROUTINE CONGRAD SOLVES THE EXTERNAL MODE BY A CONJUGATE
! **  GRADIENT SCHEME
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      DIMENSION PNORTH(LCM),PSOUTH(LCM),TMPCG(LCM)
!
!     DIMENSION CGS(LCM),CGW(LCM),CGE(LCM),CGN(LCM),
!    &          PCG(LCM),RCG(LCM),APCG(LCM)
!
!**********************************************************************C
!
!     IF(ISTL.EQ.2)THEN
!      DO L=2,LA
!      CIT=1./CCC(L)
!      CGS(L)=CIT*CCS(L)
!      CGW(L)=CIT*CCW(L)
!      CGE(L)=CIT*CCE(L)
!      CGN(L)=CIT*CCN(L)
!      ENDDO
!     ELSE
!      DO L=2,LA
!      CIT=1./CC(L)
!      CGS(L)=CIT*CS(L)
!      CGW(L)=CIT*CW(L)
!      CGE(L)=CIT*CE(L)
!      CGN(L)=CIT*CN(L)
!      ENDDO
!     ENDIF
!
!
!     IF(ISTL.EQ.2)THEN
!      DO L=2,LA
!      CG(L)=CCC(L)
!      CGS(L)=CCS(L)
!      CGW(L)=CCW(L)
!      CGE(L)=CCE(L)
!      CGN(L)=CCN(L)
!      ENDDO
!     ELSE
!      DO L=2,LA
!      CG(L)=CC(L)
!      CGS(L)=CS(L)
!      CGW(L)=CW(L)
!      CGE(L)=CE(L)
!      CGN(L)=CN(L)
!      ENDDO
!     ENDIF
!
!     DO L=2,LA
!     LN=LNC(L)
!     LS=LSC(L)
!     FPTMP(L)=FP(L)-CG(L)*PAM(L)-CGS(L)*PAM(LS)-CGW(L)*PAM(L-1)
!    &        -CGE(L)*PAM(L+1)-CGN(L)*PAM(LN)
!     ENDDO
!
!**********************************************************************C
!
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
!
      DO L=2,LA
	  PNORTH(L)=P(LNC(L))
	  PSOUTH(L)=P(LSC(L))
	ENDDO
!
      DO L=2,LA
!       LN=LNC(L)
!       LS=LSC(L)
!       RSD=CCC(L)*P(L)+CCS(L)*PSOUTH(L)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
!    &        +CCN(L)*PNORTH(L)-FPTMP(L)
!       RCG(L)=-RSD
        RCG(L)=FPTMP(L)-CCC(L)*P(L)-CCN(L)*PNORTH(L)-CCS(L)*PSOUTH(L)-CCW(L)*P(L-1)-CCE(L)*P(L+1)
      ENDDO
!
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
!         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            RCG(LCHNU)=RCG(LCHNU)+CCCCHH(NMD)*P(LHOST)
            RCG(LHOST)=RCG(LHOST)+CCCCHH(NMD)*P(LCHNU)
          ENDIF
!         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            RCG(LCHNV)=RCG(LCHNV)+CCCCHH(NMD)*P(LHOST)
            RCG(LHOST)=RCG(LHOST)+CCCCHH(NMD)*P(LCHNV)
          ENDIF
        ENDDO
      ENDIF
!
      DO L=2,LA
        PCG(L)=RCG(L)*CCCI(L)
      ENDDO
!
      RPCG=0.
      DO L=2,LA
      RPCG=RPCG+RCG(L)*PCG(L)
      ENDDO
!
      ITER=0
!
!**********************************************************************C
!
  100 CONTINUE
!
      ITER=ITER+1
!
      DO L=2,LA
	  PNORTH(L)=PCG(LNC(L))
	  PSOUTH(L)=PCG(LSC(L))
	ENDDO
!
      DO L=2,LA
        APCG(L)=CCC(L)*PCG(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L)+CCW(L)*PCG(L-1)+CCE(L)*PCG(L+1)
      ENDDO
!
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
!         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            APCG(LCHNU)=APCG(LCHNU)+CCCCHH(NMD)*PCG(LHOST)
            APCG(LHOST)=APCG(LHOST)+CCCCHH(NMD)*PCG(LCHNU)
          ENDIF
!         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            APCG(LCHNV)=APCG(LCHNV)+CCCCHH(NMD)*PCG(LHOST)
            APCG(LHOST)=APCG(LHOST)+CCCCHH(NMD)*PCG(LCHNV)
          ENDIF
        ENDDO
      ENDIF
!
      PAPCG=0.
!OLD     RPCG=0.
!
      DO L=2,LA
      PAPCG=PAPCG+APCG(L)*PCG(L)
!OLD     RPCG=RPCG+RCG(L)*PCG(L)
      ENDDO
!
      ALPHA=RPCG/PAPCG
!
      DO L=2,LA
      P(L)=P(L)+ALPHA*PCG(L)
      ENDDO
!
!**********************************************************************C
!
      DO L=2,LA
        RCG(L)=RCG(L)-ALPHA*APCG(L)
      ENDDO
!
      DO L=2,LA
        TMPCG(L)=CCCI(L)*RCG(L)
      ENDDO
!
      RPCGN=0.
      RSQ=0.
      DO L=2,LA
        RPCGN=RPCGN+RCG(L)*TMPCG(L)
        RSQ=RSQ+RCG(L)*RCG(L)
      ENDDO
!
      IF(RSQ .LE. RSQM) GOTO 200
!
      IF(ITER .GE. ITERM)THEN
       WRITE(6,600)
	 DO L=2,LA
!	   WRITE(8,800)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),       !hnr 7/27/2009
!     &                 FPTMP(L)             !hnr 7/27/2009
       ENDDO
       STOP
      ENDIF
!
!OLD     BETA=0.
!
!OLD     DO L=2,LA
!OLD     BETA=BETA+RCG(L)*APCG(L)
!OLD     ENDDO
!
!OLD     BETA=-BETA/PAPCG
      BETA=RPCGN/RPCG
      RPCG=RPCGN
!
      DO L=2,LA
!OLD     PCG(L)=RCG(L)+BETA*PCG(L)
      PCG(L)=TMPCG(L)+BETA*PCG(L)
      ENDDO
!
      GOTO 100
!
  600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
!
!**********************************************************************C
!
! ** CALCULATE FINAL RESIDUAL
!
  200 CONTINUE
!
      DO L=2,LA
	  PNORTH(L)=P(LNC(L))
	  PSOUTH(L)=P(LSC(L))
	ENDDO
!
      RSQ=0.
!
      DO L=2,LA
      RCG(L)=CCC(L)*P(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L)+CCW(L)*P(L-1)+CCE(L)*P(L+1)-FPTMP(L)
      ENDDO
!
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
!         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            RCG(LCHNU)=RCG(LCHNU)-CCCCHH(NMD)*P(LHOST)
            RCG(LHOST)=RCG(LHOST)-CCCCHH(NMD)*P(LCHNU)
          ENDIF
!         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            RCG(LCHNV)=RCG(LCHNV)-CCCCHH(NMD)*P(LHOST)
            RCG(LHOST)=RCG(LHOST)-CCCCHH(NMD)*P(LCHNV)
          ENDIF
        ENDDO
      ENDIF
!
      DO L=2,LA
        RCG(L)=RCG(L)*CCCI(L)
      ENDDO
      DO L=2,LA
        RSQ=RSQ+RCG(L)*RCG(L)
      ENDDO
!
      IF(ISCRAY.EQ.0)THEN
        TCONG=TCONG+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCONG=TCONG+T2TMP-T1TMP
        WTCONG=WTCONG+(WT2TMP-WT1TMP)*0.001
      ENDIF
!
  800 FORMAT(2I6,6E13.4)
!
      RETURN
      END