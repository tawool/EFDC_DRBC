!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      REAL FUNCTION CSEDRESS(DENBULK,WRSPO,VDRO,VDR,VDRC,IOPT)
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
! **  CALCULATES SURFACE EROSION RATE OF COHESIVE
! **  SEDIMENT AS A FUNCTION OF BED BULK DENSITY
!
! **  IOPT=1  BASED ON
! **
! **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT ERODIBILITY
! **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
! **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
!
! **  IOPT=2 & 3 BASED ON J. M. HAMRICK'S MODIFICATION OF
! **
! **  SANFORD, L.P., AND J. P. Y. MAA, 2001: A UNIFIED EROSION FORMULATION
! **  FOR FINE SEDIMENT, MARINE GEOLOGY, 179, 9-23.
!
! **  IOPT=4 & 5 BASED ON J. M. HAMRICK'S PARAMETERIZATION OF SEDFLUME
!     TEST DATA
! **

      IF(IOPT.EQ.1)THEN
        DENBULK=0.001*DENBULK
        IF(DENBULK.LE.1.065)THEN
          CSEDRESS=0.62
        ELSE
          TMP=0.198/(DENBULK-1.0023)
          TMP=EXP(TMP)
          CSEDRESS=6.4E-4*(10.**TMP)
        ENDIF
      ENDIF
!
      IF(IOPT.EQ.2)THEN
        CSEDRESS=WRSPO*(1.+VDRO)/(1.+VDR)
      ENDIF
!
      IF(IOPT.EQ.3)THEN
        CSEDRESS=WRSPO*(1.+VDRO)/(1.+VDRC)
      ENDIF
!
!      IF(IOPT.EQ.3)THEN
!        CSEDRESS=WRSPO*(1.+VDRO)/(1.+VDR)
!      ENDIF
!
      IF(IOPT.EQ.4)THEN
	  TMPVAL=(1.+VDRO)/(1.+VDR)
	  FACTOR=EXP(-TMPVAL)
        CSEDRESS=FACTOR*WRSPO*(1.+VDRO)/(1.+VDR)
      ENDIF
!
      IF(IOPT.EQ.5)THEN
	  TMPVAL=(1.+VDRO)/(1.+VDRC)
	  FACTOR=EXP(-TMPVAL)
        CSEDRESS=FACTOR*WRSPO*(1.+VDRO)/(1.+VDRC)
      ENDIF
!
      IF(IOPT.GE.99)THEN
        CSEDRESS=WRSPO
      ENDIF
!
      RETURN
      END