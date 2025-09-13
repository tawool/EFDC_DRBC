!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE SETSHLD(TSC,THETA,D,SSG,DSR,USC)
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
      TMP=GP/(VISC*VISC)
      DSR=D*(TMP**0.333333)
      GPD=GP*D
!
! **  SHIELDS 
!
      IF(DSR.LE.4.0)THEN
        THETA=0.24/DSR
      ENDIF
!
      IF(DSR.GT.4.0.AND.DSR.LE.10.0)THEN
        THETA=0.14/(DSR**0.64)
      ENDIF
!
      IF(DSR.GT.10.0.AND.DSR.LE.20.0)THEN
        THETA=0.04/(DSR**0.1)
      ENDIF
!
      IF(DSR.GT.20.0.AND.DSR.LE.150.0)THEN
        THETA=0.013*(DSR**0.29)
      ENDIF
!
      IF(DSR.GT.150.0)THEN
        THETA=0.055
      ENDIF
!
      TSC=GPD*THETA
      USC=SQRT(TSC)
!
      RETURN
      END