!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      REAL FUNCTION CSEDTAUB(DENBULK,TDUM,V1,V2,V3,IOPT)
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
! **  CALCULATES CRITIAL STRESS FOR BULK OR MASS EROSION OF COHESIVE
! **  SEDIMENT AS A FUNCTION OF BED BULK DENSITY
!
! **  IOPT=1  BASED ON
! **
! **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT DRODIBILITY
! **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
! **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
!
! **  IOPT=2  BASED ON
! **
! **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT DRODIBILITY
! **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
! **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
!
      IF(IOPT.EQ.1)THEN
        DENBULK=0.001*DENBULK
        IF(DENBULK.LE.1.013)THEN
          CSEDTAUB=0.0
        ELSE
          CSEDTAUB=0.001*(9.808*DENBULK-9.934)
        ENDIF
      ENDIF
!
      IF(IOPT.EQ.2)THEN
        DENBULK=0.001*DENBULK
        IF(DENBULK.LE.1.013)THEN
          CSEDTAUB=0.0
        ELSE
          CSEDTAUB=0.001*(9.808*DENBULK-9.934)
        ENDIF
      ENDIF
!
      RETURN
      END