!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      FUNCTION FHYDCN(VOID,BMECH4,BMECH5,BMECH6,IBMECHK)
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
!
      INCLUDE 'EFDC.PAR'
!
! **  FHYDCN IS HYDRAULIC CONDUCTIVITY DIVIDED BY 1+VOID RATIO
!
!     UNITS ARE METERS/SECOND
!
      IF(BMECH4.GT.0.0)THEN
!
!      k=ko*exp((void-voidko)/voidkk)
!      BMECH4=ko
!      BMECH5=voidko
!      BMECH6=voidkk
!
        TMP=(VOID-BMECH5)/BMECH6
        IF(IBMECHK.EQ.0) THEN
          FHYDCN=BMECH4*EXP(TMP)/(1.+VOID)
        ELSE
          FHYDCN=BMECH4*EXP(TMP)
        END IF
!
      ELSE
!
        FHYDCNLOG=0.00816448*(VOID**3)-0.232453*(VOID**2)+2.5759*VOID-28.581
        FHYDCN=EXP(FHYDCNLOG)
!
      END IF
!
      RETURN
      END