!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      REAL FUNCTION VALKH(HFFDG)
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      IF(HFFDG.LE.0.02)THEN
        VALKH=HFFDG*HFFDG
        RETURN
      ENDIF
!
      IF(HFFDG.GE.10.)THEN
        VALKH=HFFDG
        RETURN
      ENDIF
!
      DO NTAB=2,1001
      FTMPM1=FUNKH(NTAB-1)
      FTMP  =FUNKH(NTAB  )
      IF(FTMPM1.LE.HFFDG.AND.HFFDG.LT.FTMP)THEN
        VALKH=RKHTAB(NTAB)-(RKHTAB(NTAB)-RKHTAB(NTAB-1))*(FTMP-HFFDG)/(FTMP-FTMPM1)
        RETURN
      ENDIF
      ENDDO
!
      IF(NTAB.EQ.1001)THEN
        WRITE(6,600) RKHTAB(1001)
!        WRITE(8,600) RKHTAB(1001)     !hnr 7/27/2009
        STOP
      ENDIF
!
! **  INITIALIZE WAVE DISPERSION RELATION TABLE
!
!      DO N=1,501
!       RKH(N)=0.02*FLOAT(N-1)
!       FRKH(N)=RKH(N)*TANH(RKH(N))
!      ENDDO
!
!      WFRKH=WVFRQ*WVFRQ*DEP(N)/9.8
!      IF(WFRKH.LE.0.02) WRKH=WFRKH*WFRKH
!      IF(WFRKH.GE.10.) WRKH=WFRKH
!      IF(WFRKH.GT.0.02.AND.WFRKH.LT.10.)THEN
!        DO M=1,500
!         IF(WFRKH.GT.FRKH(M).AND.WFRKH.LT.FRKH(M+1))THEN
!           DRDF=(RKH(M+1)-RKH(M))/(FRKH(M+1)-FRKH(M))
!           WRKH=DRDF*(WFRKH-FRKH(M))+RKH(M)
!           GOTO 200
!         ENDIF
!        ENDDO
!      ENDIF
!  200 CONTINUE
!
!
!
  600 FORMAT(' WAVE DISPERSION TABLE OUT OF BOUNDS KH = ',E12.4)
!
      RETURN
      END