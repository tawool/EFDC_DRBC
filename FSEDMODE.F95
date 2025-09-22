!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      FUNCTION FSEDMODE(WS,USTOT,USGRN,RSNDM,ISNDM1,ISNDM2,IMODE)
!
! **  FSEDMODE SETS BEDLOAD (IMODE=1) AND SUSPENDED LOAD (IMODE=2)
! **  TRANSPORT FRACTIONS
!
      INCLUDE 'EFDC.PAR'
!
!     CHOOSE BETWEEEN TOTAL STRESS SHEAR VELOCITY AND GRAIN STRESS
!     SHEAR VELOCITY
!
      IF(ISNDM2.EQ.0)THEN
        US=USTOT
      ELSE
        US=USGRN
      ENDIF
!
      USDWS=US/WS
!
!     CHOOSE BETWEEN MODE OPTIONS
!
!     ISNDM1=0 SET BOTH BEDLOAD AND SUSPENDED LOAD FRACTIONS TO 1.0
!     this option allows maximum bed load and suspened load transport
!     independent of each other
!
      IF(ISNDM1.EQ.0) FSEDMODE=1.0
!
!     ISNDM1=1 SET BEDLOAD FRACTION TO 1. USE BINARY RELATIONSHIP FOR SUSPENDED
!     this option allows maximum bedload transport. suspended load is zero
!     for usdws.lt.rsndm or maximum for usdws.ge.rsndm
!     ie:  fsedmodeBL = 1 and fsedmodeSL = 0 or 1
!
      IF(ISNDM1.EQ.1)THEN
        FSEDMODE=0.
        IF(IMODE.EQ.1)THEN
          FSEDMODE=1.0
        ELSE
          IF(USDWS.GE.RSNDM)FSEDMODE=1.
        END IF
      END IF
!
!     ISNDM1=2 SET BEDLOAD FRACTION TO 1, USE LINEAR RELATIONSHIP FOR SUSPENDED
!     this option allows maximum bedload transport. suspended load ranges from zero
!     to maximum in a linear manner
!     ie:  fsedmodeBL = 1 and 0 < fsedmodeSL < 1
!
      IF(ISNDM1.EQ.2)THEN
        IF(IMODE.EQ.1)THEN
          FSEDMODE=1.0
        ELSE
          TMPVAL=((USDWS)-0.4)/9.6
          TMPVAL=MIN(TMPVAL,1.0)
          TMPVAL=MAX(TMPVAL,0.0)
          FSEDMODE=TMPVAL
        END IF
      END IF
!
!     ISNDM1=3 USE BINARY RELATIONSHIP FOR BEDLOAD AND SUSPENDED LOAD
!     this option does not allow simulataneous transport of bedload
!     and suspended load.  for usdws.lt.rsndm maximum bedload transport
!     occurs.  otherwise maximum suspended load occurs
!     ie.  either: fsedmodeBL =1 fsedmodeSL = 0
!              or: fsedmodeBL =0 fsedmodeSL = 1
!
      IF(ISNDM1.EQ.3)THEN
        FSEDMODE=0.
        IF(IMODE.EQ.1)THEN
          IF(USDWS.LT.RSNDM)FSEDMODE=1.
        ELSE
          IF(USDWS.GE.RSNDM)FSEDMODE=1.
        END IF
      END IF
!
!     ISNDM1=4 USE LINEAR RELATIONSHIP FOR BEDLOAD AND SUSPENDED LOAD
!     this options partitions total load between suspended load and
!     bedload using the linear suspended load relationships
!     ie fsedmodeBL + fsedmodeSL = 1
!
      IF(ISNDM1.EQ.4)THEN
        TMPVAL=((USDWS)-0.4)/9.6
        TMPVAL=MIN(TMPVAL,1.0)
        TMPVAL=MAX(TMPVAL,0.0)
        IF(IMODE.EQ.1)THEN
          FSEDMODE=1.-TMPVAL
        ELSE
          FSEDMODE=TMPVAL
        END IF
      END IF
!
      RETURN
      END