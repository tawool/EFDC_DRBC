!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY PAUL M CRAIG ON 21 JULY 2003
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
      SUBROUTINE TIMELOG(N,TIMEDAY)
!
      CHARACTER*8 MRMDATE,MRMTIME*10

! WRITE OUT MODEL TIME STEP AND SUN/PC SYSTEM CLOCK TIME TO TIME.LOG FILE:
!JH      CALL DATE(MRMDATE)
      CALL TIME(MRMTIME)
      ! *** WRITE OUT MODEL TIME STEP AND SYSTEM CLOCK TIME TO TIME.LOG
!      CALL DATE_AND_TIME(MRMDATE,MRMTIME)
      WRITE(9,100)N,TIMEDAY,MRMDATE,MRMTIME

  100 FORMAT(' ','N =',I10,5X,'TIMEDAY =',F12.4,5X,'DATE = ',A8,5X,'TIME = ',A10)
      RETURN
      END