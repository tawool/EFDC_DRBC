!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      FUNCTION FSTRSE(VOID,BMECH1,BMECH2,BMECH3)
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
! 11/15/2001        john hamrick       11/15/2001       john hamrick
!  added standard exponential form constitutive relationship
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
!
!----------------------------------------------------------------------C
!
      INCLUDE 'EFDC.PAR'
!
! **  FSTRSE IS WATER SPECIFIC WEIGHT NORMALIZED EFFECTIVE STRESS
!
!     UNITS ARE METERS
!
      IF(BMECH1.GT.0.0)THEN
!
!      st=sto*exp(-(void-voidso)/voidss)
!      BMECH1=sto
!      BMECH2=voidso
!      BMECH3=voidss
!
        TMP=-(VOID-BMECH2)/BMECH3
        FSTRSE=BMECH1*EXP(TMP)
!
      ELSE
!
!        FSTRSELOG=-0.0147351*(VOID**3)+0.333957*(VOID**2)-3.28661*VOID
!     &          +8.90864
        FSTRSELOG=-0.0147351*(VOID**3)+0.311854*(VOID**2)-2.96371*VOID+7.34698
        FSTRSE=EXP(FSTRSELOG)
!
      END IF
!
      RETURN
      END