!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      REAL FUNCTION CSNDSET(SND,SDEN,IOPT)
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
! **  CALCULATES HINDERED SETTLING CORRECTION FOR CLASS NS NONCOHESIVE
! **  SEDIMENT
!
      INCLUDE 'EFDC.PAR'
!
      ROPT=FLOAT(IOPT)
      CSNDSET=(1.-SDEN*SND)**ROPT
!
      RETURN
      END