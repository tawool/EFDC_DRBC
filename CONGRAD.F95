!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CONGRAD (ISTL)
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
!**********************************************************************C
!
      IF(LA.EQ.2)THEN
	  P(2)=FPTMP(2)/CCC(2)
	  RETURN
	ENDIF
!
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
!
!**********************************************************************C
!
      DO L=2,LA
	  PNORTH(L)=P(LNC(L))
	  PSOUTH(L)=P(LSC(L))
	ENDDO
!
      DO L=2,LA
        RCG(L)=FPTMP(L)-CCC(L)*P(L)-CCN(L)*PNORTH(L)-CCS(L)*PSOUTH(L)-CCW(L)*P(L-1)-CCE(L)*P(L+1)
      ENDDO
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
      PAPCG=0.
!
      DO L=2,LA
        PAPCG=PAPCG+APCG(L)*PCG(L)
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
	    CDIADOM=CCC(L)+CCE(L)+CCN(L)+CCS(L)+CCW(L)
!          WRITE(8,808)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),   !hnr 7/27/2009
!     &                CDIADOM,FPTMP(L),HU(L),HV(L)   !hnr 7/27/2009
        END DO
!  	  CLOSE(8)   !hnr 7/27/2009
	  CALL RESTOUT(1)
        STOP
      ENDIF
!
      BETA=RPCGN/RPCG
      RPCG=RPCGN
!
      DO L=2,LA
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
!
      DO L=2,LA
      RCG(L)=CCC(L)*P(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L)+CCW(L)*P(L-1)+CCE(L)*P(L+1)-FPTMP(L)
      ENDDO
!
      DO L=2,LA
        RCG(L)=RCG(L)*CCCI(L)
      ENDDO
      DO L=2,LA
        RSQ=RSQ+RCG(L)*RCG(L)
      ENDDO
!
!**********************************************************************C
!
      IF(ISDSOLV.GE.1.AND.N.EQ.1)THEN
	  OPEN(1,FILE='EQCOEFCONG.OUT')
        DO L=2,LA
	    CDIADOM=CCC(L)+CCE(L)+CCN(L)+CCS(L)+CCW(L)
      WRITE(1,808)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),CDIADOM,FPTMP(L),HP(L),HU(L),HV(L)
        END DO
	  CLOSE(1)
	ENDIF
!
!**********************************************************************C
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
  800 FORMAT(I5,8E13.4)
  808 FORMAT(2I5,10E13.4)
!
      RETURN
      END