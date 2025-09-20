!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALPSER (ISTL)
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
!**********************************************************************C
!
! ** SUBROUTINE CALPSER UPDATES TIME VARIABLE SURFACE ELEVATION 
! ** BOUNDARY CONDITIONS
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      PSERT(0)=0.
      PSERST(0)=0.
      PSERZDS(0)=0.
	PSERZDF(0)=0.
!
      DO NS=1,NPSER
!
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)/TCPSER(NS)+TBEGIN*(TCON/TCPSER(NS))
      ELSE
        TIME=TIMESEC/TCPSER(NS)
      ENDIF
      M1=MPTLAST(NS)
  100 CONTINUE
      M2=M1+1
      IF(TIME.GT.TPSER(M2,NS))THEN
       M1=M2
       GOTO 100
      ELSE
       MPTLAST(NS)=M1
      ENDIF      
!
      TDIFF=TPSER(M2,NS)-TPSER(M1,NS)
      WTM1=(TPSER(M2,NS)-TIME)/TDIFF
      WTM2=(TIME-TPSER(M1,NS))/TDIFF
      PSERT(NS)=WTM1*PSER(M1,NS)+WTM2*PSER(M2,NS)
      PSERST(NS)=WTM1*PSERS(M1,NS)+WTM2*PSERS(M2,NS)
!     WRITE(6,6000)N,PSERT(NS)
!
      ENDDO
!
 6000 FORMAT('N, PSERT = ',I6,4X,F12.4)
!
!**********************************************************************C
!
      RETURN
      END