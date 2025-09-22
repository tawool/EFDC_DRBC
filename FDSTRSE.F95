!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      FUNCTION FDSTRSE(VOID,BMECH1,BMECH2,BMECH3)
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 11/15/2001        john hamrick       11/15/2001       john hamrick
!  added standard exponential form constitutive relationship
!
!----------------------------------------------------------------------C
!
      INCLUDE 'EFDC.PAR'
!
! **  FDSTRSE IS COMPRESSION LENGTH SCALE
!
!     =-D(EFFECTIVE STRESS/WATER SPECIFIC WEIGHT)/D(VOID RATIO)
!
!     IE DERIVATIVE OF SPECIFIC WEIGHT NORMALIZED EFFECTIVE
!        STRESS WITH RESPECT TO VOID RATIO
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
        TMP=-(VOID-BMECH2)/BMECH3
        FDSTRSE=(BMECH1/BMECH3)*EXP(TMP)
!
      ELSE
!
!        FDSTRSEL=-0.00835084*(VOID**3)+0.207061*(VOID**2)
!     &          -2.57907*VOID+ 8.79215
        FSTRSEL=-0.0147351*(VOID**3)+0.311854*(VOID**2)-2.96371*VOID+7.34698
        DFSTRSEL=-0.0442053*(VOID**2)+0.623708*VOID-2.96371
        FDSTRSE=DFSTRSEL*EXP(FSTRSEL)
!
      END IF
!
      RETURN
      END