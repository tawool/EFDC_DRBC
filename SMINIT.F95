!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE SMINIT
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!XH      INSMICI=40
!XH      INSMRST=40
!XH      ISMORST=45
!XH      ISMOZB=46
!
      SMTSNAME(1) = 'SOM'
      SMTSNAME(2) = 'SIM'
      SMTSNAME(3) = 'SBF'
!
      DO L=2,LA
        SMHYST(L)=.FALSE.
      ENDDO
!
      RETURN
      END