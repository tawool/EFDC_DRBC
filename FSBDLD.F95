!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      REAL FUNCTION FSBDLD(DIASED,GPDIASED,D50,DEP,PEXP,PHID,CSHIELDS,SBDLDP,ISOPT)
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 12/03/2001        john hamrick       12/03/2001       john hamrick
!  modifed to include alternate formulas for bed load phi
!
!----------------------------------------------------------------------C
!
      INCLUDE 'EFDC.PAR'
!
! **  CALCULATES DIMNSIONLESS BED LOAD TRANSPORT COEFFICIENT
!
!     DIASED = SAND GRAIN DIAMETER
!     GPDIASED    = SAND GRAIN SPECIFIC GRAVITY
!     CSHIELDS     = SAND GRAIN SETTLING VELOCITY
!
! **  ISOPT=0  USE CONSTANT VALUE
!
      IF(ISOPT.EQ.0) FSBDLD=SBDLDP
!
! **  ISOPT=1  BASED ON
! **
! **  VAN RIJN, L. C., 1984: SEDIMENT TRANSPORT, PART I: BED
! **  LOAD TRANSPORT, J. HYDRAULIC ENGINEERING, 110, 1431-1455.
!
      IF(ISOPT.EQ.1)THEN
        RD=DIASED*SQRT(GPDIASED)*1.E6
        RD=(1./RD)**0.2
        TMP=CSHIELDS**2.1
        FSBDLD=0.053*RD/TMP
      ENDIF
!
! **  ISOPT=2  BASED ON MODIFIED ENGULAND-HANSEN FORMULA
! **  REFERENCE TO BE ADDED
!
      IF(ISOPT.EQ.2)THEN
        TMP1=(DEP/D50)**0.33333
        TMP2=(PEXP/PHID)**1.125
        FSBDLD=2.0367*TMP1*TMP2
      ENDIF
!
! **  ISOPT=3  BASED ON WU, WANG AND JIA, J. HYDR. RES. V38, 3000
!
      IF(ISOPT.EQ.3)THEN
        TMP1=0.03*((PHID/PEXP)**0.6)
        TMP2=TMP1**2.2
        FSBDLD=0.0053/TMP2
      ENDIF
!

      RETURN
      END