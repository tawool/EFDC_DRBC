!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      FUNCTION SETSTVEL(D,SSG)
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
      INCLUDE 'EFDC.PAR'
!
! **  NONCOHEASIVE SEDIMENT SETTLING AND SHIELDS CRITERIA
! **  USING VAN RIJN'S EQUATIONS
!
      VISC=1.E-6
      GP=(SSG-1.)*9.82
      GPD=GP*D
      SQGPD=SQRT(GPD)
      RD=SQGPD*D/VISC
!
! **  SETTLING VELOCITY  
!
      IF(D.LT.1.0E-4)THEN
        WSET=SQGPD*RD/18.
      ENDIF
!
      IF(D.GE.1.0E-4.AND.D.LT.1.E-3)THEN
        TMP=SQRT(1.+0.01*RD*RD)-1.
        WSET=10.0*SQGPD*TMP/RD
      ENDIF
!
      IF(D.GE.1.E-3)THEN
        WSET=1.1*SQGPD
      ENDIF
!
      SETSTVEL=WSET
!
      RETURN
      END