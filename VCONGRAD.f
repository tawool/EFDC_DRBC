!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE VCONGRAD(ISTL)
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
! **  SUBROUTINE VCONGRAD SOLVES THE EXTERNAL MODE LINEAR SYSTEM
! **  BY A PRECONDITIONED CONJUGATE GRADIENT SCHEME USING
! **  ALTIVEC VECTOR PROCESSING
!
! **  THE LINEAR SYSTEM P IS SPECIFIED ON A 5 POINT STENCIL BY
!
!     CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)+CCN(L)*P(LN)
!     =FPTMP(L)
!
!     WHERE LN=LNC(L) AND LS=LSC(L) ARE LOOK UP TABLES FOR THE NORTH
!     AND SOUTH STENCIL POINTS.
!
!     THE ACTUAL SIZE OF THE SYSTEM IS LC-2, WITH VARIABLES OCCUPING
!     POSITIONS (2:LC-1) IN THE ARRAYS WHICH MUST BE DIMENSIONED
!     AT LCM = OR > THAN LC.  THE LOCATIONS 1 AND LC ARE DUMMIES SUCH
!     THAT L-1 AND L+1 ARE DEFINED AND ARE DEFINED AND CAN BE SET TO
!     ZERO.  VECTOR ALTIVEC OPERATIONS MUST THE CARRIED OUT ON THE
!     FULL STORAGE DIMENSION, LCM WHICH IS AN INTEGER MULTIPE OF 4, OF
!     THE ARRAYS FOR PROPER DIMENSIONING OF THE C VECTOR SUBROUTINE
!     ARGUMENTS AS VECTOR FLOATS.
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
! ** HAMRICK'S C VECTOR SUBROUTINES MUST BE LINKED WITH EXECUTABLE
! ** FOR VECTOR EXECUTION ON MAC G4.  FOR PC'S, A DUMMY SUBROUTINE
! ** CVMULAD2 IS INCLUDED IN THE SOURCE CODE
!
!     CVMULAD2(VAL(N),A,B,C,D) => [A]=[B]+[C]*[D]
!
!     [A],[B],[C] AND [D] ARE 1D REAL*4 ARRAYS OF DIMENSION (1:LCM)
!
!**********************************************************************C
!
! **  START TIMER
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
! **  INITIALIZATION
!
      IVECDIG=0
      LCM4=LCM/4
!
      P(1)=0.0
      PW(1)=0.0
      PE(1)=0.0
      PS(1)=0.0
      PN(1)=0.0
      RCG(1)=0.0
      PCG(1)=0.0
      APCG(1)=0.0
      RA4(1)=0.0
      RTMP(1)=0.0
      RTMP1(1)=0.0
      RNULL(1)=0.0
      FPTMP(1)=0.0
      PCGW(1)=0.0
      PCGE(1)=0.0
      PCGS(1)=0.0
      PCGN(1)=0.0
      CCC(1)=0.0
      CCW(1)=0.0
      CCE(1)=0.0
      CCS(1)=0.0
      CCN(1)=0.0
      CCCI(1)=0.0
!
      DO L=LC,LCM
        P(L)=0.0
        PW(L)=0.0
        PE(L)=0.0
        PS(L)=0.0
        PN(L)=0.0
        RCG(L)=0.0
        PCG(L)=0.0
        APCG(L)=0.0
        RA4(L)=0.0
        RTMP(L)=0.0
        RTMP1(L)=0.0
        RNULL(L)=0.0
        FPTMP(L)=0.0
        PCGW(L)=0.0
        PCGE(L)=0.0
        PCGS(L)=0.0
        PCGN(L)=0.0
        CCC(L)=0.0
        CCW(L)=0.0
        CCE(L)=0.0
        CCS(L)=0.0
        CCN(L)=0.0
        CCCI(L)=0.0
      ENDDO
!
      DO L=2,LA
        LN=LNC(L)
        LS=LSC(L)
        PW(L)=P(L-1)
        PE(L)=P(L+1)
        PS(L)=P(LS)
        PN(L)=P(LN)
      ENDDO
!
!**********************************************************************C
!
!     NOTE: LA=LC-1
!
!     DO L=2,LA
!       LN=LNC(L)
!       LS=LSC(L)
!       RCG(L)=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
!    &     +CCN(L)*P(LN)-FPTMP(L)
!     ENDDO
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 0'      !hnr 7/27/2009
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
!          WRITE(8,800)L,CCW(L),CCC(L),CCE(L),CCS(L),CCCI(L),CCN(L),P(L)     !hnr 7/27/2009
        END DO
      ENDIF
!
      DO L=2,LA
        RCG(L)=-FPTMP(L)
        RNULL(L)=0.0
        RTMP(L)=0.0
        RTMP1(L)=0.0
      END DO
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 1'      !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)THEN
!        DO L=1,LC
!          WRITE(8,800)L,RCG(L)     !hnr 7/27/2009
!        END DO
!      ENDIF
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,RCG,CCC,P)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 2'     !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)THEN
!        DO L=1,LC
!          WRITE(8,800)L,RTMP(L)       !hnr 7/27/2009
!        END DO
!      ENDIF
!
!val      CALL CVMULAD2(VAL(LCM),RCG,RTMP,CCS,PS)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 3'      !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)THEN
!        DO L=1,LC
!          WRITE(8,800)L,RCG(L)       !hnr 7/27/2009
!        END DO
!      ENDIF
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,RCG,CCW,PW)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 4'     !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)THEN
!        DO L=1,LC
!          WRITE(8,800)L,RTMP(L)       !hnr 7/27/2009
!        END DO
!      ENDIF
!
!val      CALL CVMULAD2(VAL(LCM),RCG,RTMP,CCE,PE)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 5'  !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)THEN
!        DO L=1,LC
!          WRITE(8,800)L,RCG(L)  !hnr 7/27/2009
!        END DO
!      ENDIF
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,RCG,CCN,PN)
!
!     DO L=2,LA
!       RCG(L)=-RCG(L)
!       PCG(L)=RCG(L)*CCCI(L)
!     ENDDO
!
      DO L=2,LA
        RCG(L)=-RTMP(L)
      END DO
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 6'   !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)THEN
!        DO L=1,LC
!          WRITE(8,800)L,RCG(L)      !hnr 7/27/2009
!        END DO
!      ENDIF
!
!val      CALL CVMULAD2(VAL(LCM),PCG,RNULL,RCG,CCCI)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 7'   !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)THEN
!        DO L=1,LC
!          WRITE(8,800)L,PCG(L)     !hnr 7/27/2009
!        END DO
!      ENDIF
!
!     RPCG=0.
!     DO L=2,LA
!       RPCG=RPCG+RCG(L)*RCG(L)*CCCI(L)
!     ENDDO
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,RNULL,RCG,CCCI)
!val      CALL CVMULAD2(VAL(LCM),RTMP1,RNULL,RCG,RTMP)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 8'    !hnr 7/27/2009
!
      RPCG=0.
      DO L=2,LA
        RPCG=RPCG+RTMP1(L)
      ENDDO
!
!      IF(IVECDIG.EQ.1)WRITE(8,800)LA,RPCG      !hnr 7/27/2009
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
        LN=LNC(L)
        LS=LSC(L)
        PCGS(L)=PCG(LS)
        PCGW(L)=PCG(L-1)
        PCGE(L)=PCG(L+1)
        PCGN(L)=PCG(LN)
      ENDDO
!
!     DO L=2,LA
!       LN=LNC(L)
!       LS=LSC(L)
!       APCG(L)=CCC(L)*PCG(L)+CCS(L)*PCG(LS)+CCW(L)*PCG(L-1)
!    &         +CCE(L)*PCG(L+1)+CCN(L)*PCG(LN)
!     ENDDO
!
!val      CALL CVMULAD2(VAL(LCM),APCG,RNULL,CCC,PCG)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 9'  !hnr 7/27/2009
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,APCG,CCS,PCGS)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 10'  !hnr 7/27/2009
!
!val      CALL CVMULAD2(VAL(LCM),APCG,RTMP,CCW,PCGW)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 11'  !hnr 7/27/2009
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,APCG,CCE,PCGE)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 12' !hnr 7/27/2009
!
!val      CALL CVMULAD2(VAL(LCM),APCG,RTMP,CCN,PCGN)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 13'  !hnr 7/27/2009
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
!          WRITE(8,800)L,APCG(L)   !hnr 7/27/2009
        END DO
      END IF
!
!     PAPCG=0.
!     DO L=2,LA
!       PAPCG=PAPCG+APCG(L)*PCG(L)
!     END DO
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,RNULL,APCG,PCG)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 14'     !hnr 7/27/2009
!
      PAPCG=0.
      DO L=2,LA
        PAPCG=PAPCG+RTMP(L)
        RTMP1(L)=P(L)
      ENDDO
!
!      IF(IVECDIG.EQ.1)WRITE(8,800)LA,PAPCG     !hnr 7/27/2009
!
      ALPHA=RPCG/PAPCG
!
!     DO L=2,LA
!       P(L)=P(L)+ALPHA*PCG(L)
!     ENDDO
!
      DO L=2,LA
        RA4(L)=ALPHA
      END DO
!
!val      CALL CVMULAD2(VAL(LCM),P,RTMP1,RA4,PCG)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 15'    !hnr 7/27/2009
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
!          WRITE(8,800)L,P(L)     !hnr 7/27/2009
        END DO
      END IF
!
!**********************************************************************C
!
!     DO L=2,LA
!       RCG(L)=RCG(L)-ALPHA*APCG(L)
!     ENDDO
!
      DO L=2,LA
        RA4(L)=-ALPHA
        RTMP(L)=RCG(L)
      END DO
!
!val      CALL CVMULAD2(VAL(LCM),RCG,RTMP,RA4,APCG)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 16'    !hnr 7/27/2009
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
!          WRITE(8,800)L,RCG(L)  !hnr 7/27/2009
        END DO
      END IF
!
!     RPCGN=0.
!     DO L=2,LA
!       RPCGN=RPCGN+RCG(L)*RCG(L)*CCCI(L)
!     ENDDO
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,RNULL,RCG,CCCI)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 17'    !hnr 7/27/2009
!
!val      CALL CVMULAD2(VAL(LCM),RTMP1,RNULL,RTMP,RCG)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 18'     !hnr 7/27/2009
!
      RPCGN=0.
      DO L=2,LA
        RPCGN=RPCGN+RTMP1(L)
        RTMP(L)=RCG(L)
      ENDDO
!
!      IF(IVECDIG.EQ.1)WRITE(8,800)LA,RPCGN   !hnr 7/27/2009
!
!     RSQ=0.
!     DO L=2,LA
!       RSQ=RSQ+RCG(L)*RCG(L)
!     ENDDO
!
!val      CALL CVMULAD2(VAL(LCM),RTMP1,RNULL,RTMP,RCG)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 19'  !hnr 7/27/2009
!
      RSQ=0.
      DO L=2,LA
        RSQ=RSQ+RTMP1(L)
      ENDDO
!
!      IF(IVECDIG.EQ.1)WRITE(8,800)LA,RSQ  !hnr 7/27/2009
!
      IF(RSQ.LE.RSQM)GOTO 200
!
      IF(ITER.GE.ITERM)THEN
        WRITE(6,600)
        STOP
      ENDIF
!
      BETA=RPCGN/RPCG
      RPCG=RPCGN
!
!     DO L=2,LA
!       PCG(L)=CCCI(L)*RCG(L)+BETA*PCG(L)
!     ENDDO
!
      DO L=2,LA
        RA4(L)=BETA
      ENDDO
!
!val      CALL CVMULAD2(VAL(LCM),RTMP,RNULL,PCG,RA4)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 20'    !hnr 7/27/2009
!
!val      CALL CVMULAD2(VAL(LCM),PCG,RTMP,CCCI,RCG)
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 21'    !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)THEN
!        DO L=1,LC
!          WRITE(8,800)L,PCG(L)      !hnr 7/27/2009
!        ENDDO
!      ENDIF
!
      GOTO 100
!
  600 FORMAT(1X,'MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
!
!**********************************************************************C
!
! ** CALCULATE FINAL RESIDUAL
!
  200 CONTINUE
!
      RSQ=0.
!
      DO L=2,LA
        LN=LNC(L)
        LS=LSC(L)
        RSD=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
     $        +CCN(L)*P(LN)-FPTMP(L)
        RSD=RSD*CCCI(L)
        RSQ=RSQ+RSD*RSD
      ENDDO
!
!      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 22'  !hnr 7/27/2009
!      IF(IVECDIG.EQ.1)WRITE(8,800)LA,RSQ     !hnr 7/27/2009
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
!
!**********************************************************************C
!
      RETURN
      END