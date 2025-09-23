!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE SKIPCOMM(IUNIT, CC)
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
! SKIPS OVER COMMENT LINES IN INPUT FILES
!
      INTEGER IUNIT
      CHARACTER CC*1, LINE*80
100   READ(IUNIT, *, END=999) LINE
      IF(LINE(1:1) .EQ. CC) GOTO 100
      IF(LINE(1:1) .EQ. 'C') GOTO 100
      IF(LINE(1:1) .EQ. 'C') GOTO 100
      BACKSPACE(IUNIT)
999   RETURN
      END