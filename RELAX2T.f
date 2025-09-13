!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE RELAX2T
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 02/15/2002        John Hamrick       02/15/2002       John Hamrick
!  added this subroutine relax2t
!----------------------------------------------------------------------C
!
! **  SUBROUTINE RELAX SOLVES THE FINITE DIFFERENCE FORM
! **  OF A PSEUDO HEMHOLTZ EQUATION
! **
! **              CS(L)*P(LS)+CW(L)*P(L-1)                
! **              +CC(L)*P(L)+CE(L)*P(L+1)                
! **              +CN(L)*P(LN) = FP(L)                    
! **                                                      
! **  BY SUCCESSIVE OVER RELAXATION USING A RED-BLACK ORDERING 
! **  WITH CONVERGENCE MEASURED BY A GLOBAL SQUARE ERROR RSQ.  
! **  NON-CONVERGENCE IS SIGNALED WHEN THE ITERATIONS EXCEED A 
! **  MAXIMUM.                                                 
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      RJ2=RP
!
!      PAVG=0.0
!      DO L=2,LA
!      PAVG=PAVG+P(L)
!      ENDDO
!	PAVG=PAVG/FLOAT(LA-1)
!
!      DO L=2,LA
!        FPTMP(L)=FPTMP(L)-PAVG*(CCC(L)+CCS(L)+CCW(L)+CCE(L)+CCN(L))
!        P(L)=P(L)-PAVG
!      ENDDO
!
      FPSQ=0.
      DO L=2,LA
        FPSQ=FPSQ+FPTMP(L)*FPTMP(L)
      ENDDO
!
      ITER=1
!
  200 CONTINUE
      RSQ=0.
!
!**********************************************************************C
!
! **  RED CELL LOOP
!
      IF(ITER.EQ.1) RPT=1.0
	IF(ITER.GT.1) RPT=1.0/(1.0-0.25*RJ2*RPT)
!
      DO L=2,LA
      K=IL(L)+JL(L)
      IVAL=MOD(K,2)
      IF(IVAL.EQ.0)THEN
        LN=LNC(L)
        LS=LSC(L)
        RSD=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)+CCN(L)*P(LN)-FPTMP(L)
        P(L)=P(L)-RPT*RSD/CCC(L)
        RSQ=RSQ+RSD*RSD
      ENDIF
      ENDDO
!
!**********************************************************************C
!
! **  BLACK CELL LOOP
!
!
      IF(ITER.EQ.1) RPT=1.0/(1.0-0.5*RJ2)
	IF(ITER.GT.1) RPT=1.0/(1.0-0.25*RJ2*RPT)
!
      DO L=2,LA
      K=IL(L)+JL(L)
      IVAL=MOD(K,2)
      IF(IVAL.NE.0)THEN
        LN=LNC(L)
        LS=LSC(L)
        RSD=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
     &        +CCN(L)*P(LN)-FPTMP(L)
        P(L)=P(L)-RPT*RSD/CCC(L)
        RSQ=RSQ+RSD*RSD
      ENDIF
      ENDDO
!
!**********************************************************************C
!
! **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
!
      RSQ=SQRT(RSQ)/SQRT(FPSQ)
      IF(RSQ .LE. RSQM) GOTO 400
!
! **  CHECK MAXIMUM ITERATION CRITERIA
!
      IF(ITER .GE. ITERM)THEN
       WRITE(6,600)
       WRITE(6,601)RSQ
!       WRITE(8,600)         !hnr 7/27/2009
!       WRITE(8,601)RSQ    !hnr 7/27/2009
       DO L=2,LA
!         WRITE(8,800)L,CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),P(L),FPTMP(L)     !hnr 7/27/2009
       ENDDO

       STOP
      ENDIF
!  
      ITER=ITER+1
      GOTO 200
!
  400 CONTINUE
!
!      DO L=2,LA
!	 P(L)=P(L)+PAVG
!      END DO
!
!**********************************************************************C
!
  600 FORMAT(' MAX ITERATIONS EXCEEDED IN EXTERNAL SOLUTION, relax2t')
  601 FORMAT(' RSQ = ',E14.5)
  800 FORMAT(I6,7E13.5)
!
!**********************************************************************C
!
      RETURN
      END