!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE DEPPLT
!
! **  SUBROUTINE DEPPLT WRITES A FILE TO CONTOUR PLOT DEPTH
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      CHARACTER*80 TITLE
!
!**********************************************************************C
!
      OPEN(1,FILE='BELVCON.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BELVCON.OUT',STATUS='UNKNOWN')
!
      TITLE='BOTTOM ELEVATION CONTOURS'
!
      NTIME=1
      LINES=LA-1
      LEVELS=1
      DBS=0.
!
      WRITE (1,99) TITLE
      WRITE (1,101)LINES,LEVELS
      WRITE (1,250)DBS
!
      WRITE (1,100)NTIME
!
      DO L=2,LA
      WRITE(1,200)IL(L),JL(L),DLON(L),DLAT(L),BELV(L)
      ENDDO
!
      CLOSE(1)
!
!**********************************************************************C
!
   99 FORMAT(A80)
  100 FORMAT(I10)
  101 FORMAT(2I10)
  200 FORMAT(2I4,1X,2F15.3,F10.3)
!  200 FORMAT(2I4,1X,10F12.6)
  250 FORMAT(12F10.6)
!
!**********************************************************************C
!
      RETURN
      END