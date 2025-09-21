!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      REAL FUNCTION CSEDRESB(DENBULK,WRSPO,VDRO,VDR,VDRC,IOPT)
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
! **  CALCULATES BULK EROSION RATE OF COHESIVE
! **  SEDIMENT AS A FUNCTION OF BED BULK DENSITY AND OTHER VARIABLES
! **  CURRENT OPTIONS SHOULD NOT BE USED
!
! **  IOPT=1  BASED ON
! **
! **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT ERODIBILITY
! **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
! **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
!
! **  IOPT=2  BASED ON J. M. HAMRICK'S MODIFICATION OF
! **
! **  SANFORD, L.P., AND J. P. Y. MAA, 2001: A UNIFIED EROSION FORMULATION
! **  FOR FINE SEDIMENT, MARINE GEOLOGY, 179, 9-23.
!
      IF(IOPT.GE.1)CSEDRESB=0.0
!
!      IF(IOPT.EQ.1)THEN
!        DENBULK=0.001*DENBULK
!        IF(DENBULK.LE.1.065)THEN
!          CSEDRESB=0.62
!        ELSE
!          TMP=0.198/(DENBULK-1.0023)
!          TMP=EXP(TMP)
!          CSEDRESB=6.4E-4*(10.**TMP)
!        ENDIF
!      ENDIF
!
!      IF(IOPT.EQ.2)THEN
!        CSEDRESB=WRSPO*(1.+VDRO)/(1.+VDR)
!      ENDIF
!
!      IF(IOPT.EQ.3)THEN
!        CSEDRESB=WRSPO*(1.+VDRO)/(1.+VDRC)
!      ENDIF
!
!      IF(IOPT.EQ.99)THEN
!        CSEDRESB=WRSPO
!      ENDIF
!
      RETURN
      END