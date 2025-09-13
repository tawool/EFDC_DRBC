!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE WWQRST
!
!**********************************************************************C
!
! **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
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
!
!**********************************************************************C
!
! WRITE SPATIAL DISTRIBUTIONS AT THE END OF SIMULATION TO UNIT IWQORST.
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!      REAL XTEM
!      LOGICAL FEXIST
!
! WRITE ASCII RESTART FILE:
!
      OPEN(1,FILE='WQWCRST.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WQWCRST.OUT',STATUS='UNKNOWN')
!
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
      WRITE(1,101) N,TIME
      WRITE(1,102)
!
! J.S.
      NWQV0=NWQV
      IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
      DO L=2,LA
        DO K=1,KC
          WRITE(1,90) L,K,(WQV(L,K,NW),NW=1,NWQV0)
        ENDDO
      ENDDO
!
      CLOSE(1)
!
! ALSO WRITE BINARY RESTART FILE:
!
!      INQUIRE(FILE='WQWCRST.BIN', EXIST=FEXIST)
!      IF(FEXIST)THEN
!        OPEN(UNIT=1, FILE='WQWCRST.BIN', ACCESS='TRANSPARENT',
!     +    FORM='UNFORMATTED', STATUS='UNKNOWN')
!        CLOSE(UNIT=1, STATUS='DELETE')
!      ENDIF
!      OPEN(UNIT=1, FILE='WQWCRST.BIN', ACCESS='TRANSPARENT',
!     +   FORM='UNFORMATTED', STATUS='UNKNOWN')
!      WRITE(1) N, TIME
!      NWQV0=NWQV
!      IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
!      DO L=2,LA
!        DO K=1,KC
!          WRITE(1) L, K
!          DO NW=1,NWQV0
!            WRITE(1) WQV(L,K,NW)
!          ENDDO
!        ENDDO
!      ENDDO
!      CLOSE(1)
!
   90 FORMAT(2I5, 1P, 22E12.4)
  101 FORMAT('CC  WQ RESTART FILE TIME STEP, TIME = ',I10,F12.5)
  102 FORMAT('C   L    K  BC          BD          BG          ',
     &       'RPOC        LPOC        DOC         ',
     &       'RPOP        LPOP        DOP         PTO4        ',
     &       'RPON        LPON        DON         AMN         ',
     &       'NIT         SU          SA          COD         ',
     &       'DO          TAM         FCB        MALG')
!
      RETURN
      END