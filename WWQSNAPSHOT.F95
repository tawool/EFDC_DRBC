!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE WWQSNAPSHOT
!
!**********************************************************************C
!
! **  SHEN'S MODIFICATION TO OUTPUT MACROALGAE
!
! **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
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
! hugo rodriguez  10/2010  algae chla/C time variable in algae kinetic file
!
!**********************************************************************C
!
! WRITE TIME-SERIES OUTPUT: WQCHLX=1/WQCHLX
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      CHARACTER*12 FNWQSS(22)
      DIMENSION WQVOUT(NWQVM+6)
!
      FNWQSS( 1)='wqCHLaC1.OUT'
      FNWQSS( 2)='wqCHLaG2.OUT'
      FNWQSS( 3)='wqCHLaD3.OUT'
      FNWQSS( 4)='wqRPOCar.OUT'
      FNWQSS( 5)='wqLPOCar.OUT'
      FNWQSS( 6)='wqTDOCar.OUT'
      FNWQSS( 7)='wqRPOPhr.OUT'
      FNWQSS( 8)='wqLPOPhr.OUT'
      FNWQSS( 9)='wqTDOPhr.OUT'
      FNWQSS(10)='wqP04Tot.OUT'
      FNWQSS(11)='wqRPONit.OUT'
      FNWQSS(12)='wqLPONit.OUT'
      FNWQSS(13)='wqTDONit.OUT'
      FNWQSS(14)='wqNH4Nit.OUT'
      FNWQSS(15)='wqNO3Nit.OUT'
      FNWQSS(16)='wqSilicA.OUT'
      FNWQSS(17)='wqSilicU.OUT'
      FNWQSS(18)='wqChemOD.OUT'
      FNWQSS(19)='wqDisOxy.OUT'
      FNWQSS(20)='wqTAMeta.OUT'
      FNWQSS(21)='wqFColBa.OUT'
      FNWQSS(22)='wqBenAlg.OUT'
!
!
!      TINDAY=TINDAY
!
      OPEN(1,FILE='WQWCSNAPSHOT.OUT',STATUS='UNKNOWN',POSITION='APPEND')
!
      IF(ISDYNSTP.EQ.0)THEN
        TIMTMP=DT*FLOAT(N)+TCON*TBEGIN
        TIMTMP=TIMTMP/TCTMSR  
      ELSE
        TIMTMP=TIMESEC/TCTMSR  
      ENDIF
!
      DO K=1,KC
        WRITE(1,101)TIMTMP,K
	  DO L=2,LA
	  iz=IWQZMAP(L,K)                                               !hnr  10/2010
          DO NW=1,NWQV+1
            WQVO(L,K,NW)=WQVO(L,K,NW)*0.5
          ENDDO
          NWQOUT=0
          IF(ISTRWQ(1).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(L,K,1)*WQCHLC(iz)                      !hnr 10/2010
          ENDIF
          IF(ISTRWQ(2).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(L,K,2)*WQCHLD(iz)                      !hnr  10/2010
          ENDIF
          IF(ISTRWQ(3).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(L,K,3)*WQCHLG(iz)                      !hnr  10/2010
          ENDIF
          DO NW=4,NWQV 
            IF(ISTRWQ(NW).EQ.1)THEN
              NWQOUT=NWQOUT+1
              WQVOUT(NWQOUT)=WQVO(L,K,NW)
            ENDIF
          ENDDO
          IF(IDNOTRVA.GT.0)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(L,K,IDNOTRVA)
          ENDIF
          DO NW=1,NWQV+1
            WQVO(L,K,NW)=WQVO(L,K,NW)*2.0
          ENDDO
	    WRITE(1,102)(WQVOUT(NW),NW=1,NWQOUT)
	  ENDDO
	ENDDO
!
      CLOSE(1)
!
  101 FORMAT(F12.5,I6)
  102 FORMAT(25E12.4)
!
      RETURN
      END