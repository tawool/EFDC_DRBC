!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE ADJMMT
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
! **  SUBROUTINE ADJMMT ADJUST THE MEAN MASS TRANSPORT FIELD
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
! **  SEPARATE INTERNAL AND EXTERNAL MODES AND AVERAGE TO MIDDLE OF
! **  TIME INTERVAL
!
!----------------------------------------------------------------------C
!
      DO L=1,LC
      UHDY2E(L)=0.5*(UHDYE(L)+UHDY1E(L))
      VHDX2E(L)=0.5*(VHDXE(L)+VHDX1E(L))
      ENDDO
!
      DO K=1,KC
      DO L=1,LC
      UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY1(L,K))-UHDY2E(L)
      VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX1(L,K))-VHDX2E(L)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      W2(L,K)=W2(L,K-1)-DZC(K)*DXYIP(L)*(UHDY2(L+1,K)-UHDY2(L,K)+VHDX2(LN,K)-VHDX2(L,K))
      ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE AND OUTPUT THE DIVERGENCE OF THE EXTERNAL
! **  EULERIAN RESIDUAL TRANSPORT FIELD
!
!----------------------------------------------------------------------C
!
  100 CONTINUE
!
      RNTCMMT=FLOAT(NTSMMT)/FLOAT(NTSPTC)
      DO L=2,LA
      LN=LNC(L)
      DIVERTE(L)=SPB(L)*( ((HLPF(L)-H1P(L))/(RNTCMMT*TIDALP))*DXYP(L)   &
           -QSUMELPF(L)+UHDY2E(L+1)-UHDY2E(L)+VHDX2E(LN)-VHDX2E(L) )
      ENDDO
!
      DIVERTEG=0.
      DO L=2,LA
      DIVERTEG=DIVERTEG+DIVERTE(L)
      FP(L)=-CC(L)*DIVERTE(L)
      ENDDO
!
!----------------------------------------------------------------------C

      DIVMAX=-1.E-20
      DIVMIN=1.E+20
      DO L=2,LA
      IF(DIVERTE(L).GT.DIVMAX)THEN
       DIVMAX=DIVERTE(L)
       IMAX=IL(L)
       JMAX=JL(L)
      ENDIF
      IF(DIVERTE(L).LT.DIVMIN)THEN
       DIVMIN=DIVERTE(L)
       IMIN=IL(L)
       JMIN=JL(L)
      ENDIF
      ENDDO
!
      ITER=0
      WRITE(6,6001)ITER,DIVERTEG
      WRITE(6,6002)DIVMAX,IMAX,JMAX
      WRITE(6,6003)DIVMIN,IMIN,JMIN
!
 6001 FORMAT('  ITER=',I10,3X,'DIVERTEG=',E14.4)
 6002 FORMAT('  DIVMAX=',E14.4,5X,'I,J=',(2I10))
 6003 FORMAT('  DIVMIN=',E14.4,5X,'I,J=',(2I10))
!
!**********************************************************************C
!
! **  INSERT BOUNDARY CONDITIONS AND SAVE OUTFLOWS IN UHDY1E AND VHDX1E
!
!----------------------------------------------------------------------C
!
! **  OPEN FILE TO WRITE DIAGNOSTICS
!
      OPEN(1,FILE='ADJMMT.DIA',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='ADJMMT.DIA',STATUS='UNKNOWN')
!
      QE=0.
      QW=0.
      QN=0.
      QS=0.
      QT=0.
      WTS=0.
      WTW=0.
      WTE=0.
      WTN=0.
      WTT=0.
!
      DO LL=1,NPBE
      L=LPBE(LL)
      P(L)=0.
      QE=QE+UHDY2E(L)
      WTE=WTE+ABS(UHDY2E(L))
      ENDDO
!
      DO LL=1,NPBW
      L=LPBW(LL)
      P(L)=0.
      QW=QW+UHDY2E(L+1)
      WTW=WTW+ABS(UHDY2E(L+1))
      ENDDO
!
      DO LL=1,NPBN
      L=LPBN(LL)
      P(L)=0.
      QN=QN+VHDX2E(L)
      WTN=WTN+ABS(VHDX2E(L))
      ENDDO
!
      DO LL=1,NPBS
      L=LPBS(LL)
      P(L)=0.
      LN=LNC(L)
      QS=QS+VHDX2E(LN)
      WTS=WTS+ABS(VHDX2E(LN))
      ENDDO
!
      QE=ABS(QE)
      QW=ABS(QW)
      QN=ABS(QN)
      QS=ABS(QS)
      QT=1./(QE+QW+QS+QN)
      QE=QE*QT*DIVERTEG/WTE
      QW=QW*QT*DIVERTEG/WTW
      QN=QN*QT*DIVERTEG/WTN
      QS=QS*QT*DIVERTEG/WTS
!
      DO LL=1,NPBE
      L=LPBE(LL)
      WRITE(1,1120)IL(L),JL(L),UHDY2E(L)
      FP(L-1)=FP(L-1)+QE*ABS(UHDY2E(L))*CC(L-1)
      UHDY1E(L)=UHDY2E(L)-QE*ABS(UHDY2E(L))
      WRITE(1,1120)IL(L),JL(L),UHDY2E(L)
      WRITE(1,1112)
      ENDDO
!
      DO LL=1,NPBW
      L=LPBW(LL)
      WRITE(1,1130)IL(L),JL(L),UHDY2E(L+1)
      FP(L+1)=FP(L+1)+QW*ABS(UHDY2E(L+1))*CC(L+1)
      UHDY1E(L+1)=UHDY2E(L+1)+QW*ABS(UHDY2E(L+1))
      WRITE(1,1130)IL(L),JL(L),UHDY2E(L+1)
      WRITE(1,1112)
      ENDDO
!
      DO LL=1,NPBN
      L=LPBN(LL)
      LS=LSC(L)
      WRITE(1,1140)IL(L),JL(L),VHDX2E(L)
      FP(LS)=FP(LS)+QN*ABS(VHDX2E(L))*CC(LS)
      VHDX1E(L)=VHDX2E(L)-QN*ABS(VHDX2E(L))
      WRITE(1,1140)IL(L),JL(L),VHDX2E(L)
      WRITE(1,1112)
      ENDDO
!
      DO LL=1,NPBS
      L=LPBS(LL)
      LN=LNC(L)
      WRITE(1,1150)IL(L),JL(L),VHDX2E(LN)
      FP(LN)=FP(LN)+QS*ABS(VHDX2E(LN))*CC(LN)
      VHDX1E(LN)=VHDX2E(LN)+QS*ABS(VHDX2E(LN))
      WRITE(1,1150)IL(L),JL(L),VHDX2E(LN)
      WRITE(1,1112)
      ENDDO
!
! **  WRITE DIAGNOSTICS
!
      FPSUM=0.
      WRITE(1,1112)
      WRITE(1,1110)
      DO L=2,LA
      FPSUM=FPSUM+FP(L)
      WRITE(1,1111)IL(L),JL(L),CC(L),CS(L),CW(L),CE(L),CN(L),FP(L),QSUMELPF(L)
      ENDDO
      WRITE(1,1112)
      WRITE(1,1113)FPSUM
!
      CLOSE(1)
!
 1110 FORMAT('    I    J          CC          CS          CW','          CE          CN          FP    QSUMELPF',/)
 1111 FORMAT(2I5,7(1X,E12.4))
 1112 FORMAT(//)
 1113 FORMAT('FPSUM =  ',E12.4)
 1120 FORMAT('EAST , I,J,UHDY2E = ',2I5,2X,E12.4)
 1130 FORMAT('WEST , I,J,UHDY2E = ',2I5,2X,E12.4)
 1140 FORMAT('NORTH, I,J,VHDX2E = ',2I5,2X,E12.4)
 1150 FORMAT('SOUTH, I,J,VHDX2E = ',2I5,2X,E12.4)
!
!**********************************************************************C
!
! **  SOLVE THE FINITE DIFFERENCE FORM OF THE POTENTIAL EQUATION
! **
! **              CS(L)*P(LS)+CW(L)*P(L-1)
! **              +P(L)+CE(L)*P(L+1)
! **              +CN(L)*P(LN) = FP(L)
! **
! **  BY SUCCESSIVE OVER RELAXATION USING A RED-BLACK ORDERING
! **  WITH CONVERGENCE MEASURED BY A GLOBAL SQUARE ERROR RSQADJ
! **  NON-CONVERGENCE IS SIGNALED WHEN THE ITERATIONS EXCEED ITMADJ
!
!----------------------------------------------------------------------C
!
      ITER=1
!
  200 CONTINUE
      RSQE=0.
!
!----------------------------------------------------------------------C
!
! **  RED CELL LOOP
!
      DO LL=1,NRC
      L=LRC(LL)
      LN=LNC(L)
      LS=LSC(L)
      RSDE=CS(L)*P(LS)+CW(L)*P(L-1)+P(L)+CE(L)*P(L+1)+CN(L)*P(LN)-FP(L)
      P(L)=P(L)-RPADJ*RSDE
      RSQE=RSQE+RSDE*RSDE
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  BLACK CELL LOOP
!
      DO LL=1,NBC
      L=LBC(LL)
      LN=LNC(L)
      LS=LSC(L)
      RSDE=CS(L)*P(LS)+CW(L)*P(L-1)+P(L)+CE(L)*P(L+1)+CN(L)*P(LN)-FP(L)
      P(L)=P(L)-RPADJ*RSDE
      RSQE=RSQE+RSDE*RSDE
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
!
      IF(RSQE.LE.RSQMADJ) GOTO 800
!
!----------------------------------------------------------------------C
!
! **  CHECK MAXIMUM ITERATION CRITERIA
!
      IF(ITER .GE. ITRMADJ) GOTO 800
!     IF(ITER .GE. ITRMADJ)THEN
!
!     DO L=2,LA
!     LS=LSC(L)
!     UHDY2E(L)=UHDY2E(L)+HRU(L)*(P(L)-P(L-1))
!     VHDX2E(L)=VHDX2E(L)+HRV(L)*(P(L)-P(LS))
!     ENDDO
!
!     DO LL=1,NPBE
!     L=LPBE(LL)
!     UHDY2E(L)=UHDY1E(L)
!     ENDDO
!
!     DO LL=1,NPBW
!     L=LPBW(LL)
!     UHDY2E(L+1)=UHDY1E(L+1)
!     ENDDO
!
!     DO LL=1,NPBN
!     L=LPBN(LL)
!     VHDX2E(L)=VHDX1E(L)
!     ENDDO
!
!     DO LL=1,NPBS
!     L=LPBS(LL)
!     LN=LNC(L)
!     VHDX2E(LN)=VHDX1E(LN)
!     ENDDO
!
!     GOTO 100
!     ENDIF
!
      ITER=ITER+1
      GOTO 200
!
!**********************************************************************C
!
! **  CORRECT EULERIAN RESIDUAL TRANSPORT VELOCITY FIELD
!
!----------------------------------------------------------------------C
!
  800 CONTINUE
!
      DO L=2,LA
      LS=LSC(L)
      UHDY2E(L)=UHDY2E(L)+HRU(L)*(P(L)-P(L-1))
      VHDX2E(L)=VHDX2E(L)+HRV(L)*(P(L)-P(LS))
      ENDDO
!
      DO LL=1,NPBE
      L=LPBE(LL)
      UHDY2E(L)=UHDY1E(L)
      ENDDO
!
      DO LL=1,NPBW
      L=LPBW(LL)
      UHDY2E(L+1)=UHDY1E(L+1)
      ENDDO
!
      DO LL=1,NPBN
      L=LPBN(LL)
      VHDX2E(L)=VHDX1E(L)
      ENDDO
!
      DO LL=1,NPBS
      L=LPBS(LL)
      LN=LNC(L)
      VHDX2E(LN)=VHDX1E(LN)
      ENDDO
!
      DO K=1,KC
      DO L=1,LC
      UHDY2(L,K)=UHDY2(L,K)+UHDY2E(L)
      VHDX2(L,K)=VHDX2(L,K)+VHDX2E(L)
      U2(L,K)=UHDY2(L,K)/(HMU(L)*DYU(L))
      V2(L,K)=VHDX2(L,K)/(HMV(L)*DXV(L))
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      RNTCMMT=FLOAT(NTSMMT)/FLOAT(NTSPTC)
      DO L=2,LA
      LN=LNC(L)
      DIVERTE(L)=SPB(L)*( ((HLPF(L)-H1P(L))/(RNTCMMT*TIDALP))*DXYP(L)-QSUMELPF(L)+UHDY2E(L+1)     &
                 -UHDY2E(L)+VHDX2E(LN)-VHDX2E(L) )
      ENDDO
!
      DIVERTEG=0.
      DO L=2,LA
      DIVERTEG=DIVERTEG+DIVERTE(L)
      ENDDO
!
      DIVMAX=-1.E-20
      DIVMIN=1.E+20
      DO L=2,LA
      IF(DIVERTE(L).GT.DIVMAX)THEN
       DIVMAX=DIVERTE(L)
       IMAX=IL(L)
       JMAX=JL(L)
      ENDIF
      IF(DIVERTE(L).LT.DIVMIN)THEN
       DIVMIN=DIVERTE(L)
       IMIN=IL(L)
       JMIN=JL(L)
      ENDIF
      ENDDO
!
      WRITE(6,6001)ITER,DIVERTEG
      WRITE(6,6002)DIVMAX,IMAX,JMAX
      WRITE(6,6003)DIVMIN,IMIN,JMIN
!
!**********************************************************************C
!
! **  CALCULATE VECTOR POTENTIAL TRANSPORT AND LAGRANGIAN RESIDUAL
! **  TRANSPORT
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
!     UVPT(L,K)=(HMC(LN)*VPZ(LN,K)-HMC(L)*VPZ(L,K))/DYU(L)
!    &          -DZIC(K)*(VPY(L,K)-VPY(L,K-1))
!     VVPT(L,K)=DZIC(K)*(VPX(L,K)-VPX(L,K-1))
!    &          -(HMC(L+1)*VPZ(L+1,K)-HMC(L)*VPZ(L,K))/DXV(L)
      UVPT(L,K)=(VPZ(LN,K)-VPZ(L,K))/DYU(L)-DZIC(K)*(VPY(L,K)-VPY(L,K-1))
      VVPT(L,K)=DZIC(K)*(VPX(L,K)-VPX(L,K-1))-(VPZ(L+1,K)-VPZ(L,K))/DXV(L)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      WVPT(L,K)=(VPY(L+1,K)-VPY(L,K))/DXP(L)-(VPX(LN,K)-VPX(L,K))/DYP(L)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=1,LC
      UHLPF(L,K)=UHDY2(L,K)/DYU(L)
      VHLPF(L,K)=VHDX2(L,K)/DXV(L)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=1,LC
      UHDY2(L,K)=UHDY2(L,K)+UVPT(L,K)*DYU(L)
      VHDX2(L,K)=VHDX2(L,K)+VVPT(L,K)*DXV(L)
      ENDDO
      ENDDO
!
      DO K=1,KS
      DO L=1,LC
      W2(L,K)=W2(L,K)+WVPT(L,K)
      W(L,K)=W2(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=1,LC
      U2(L,K)=UHDY2(L,K)/(HMU(L)*DYU(L))
      V2(L,K)=VHDX2(L,K)/(HMU(L)*DXV(L))
      U(L,K)=U2(L,K)
      V(L,K)=V2(L,K)
      ENDDO
      ENDDO
!
!**********************************************************************C
!
      RETURN
      END