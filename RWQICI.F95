!
!**********************************************************************C
!
      SUBROUTINE RWQICI(IWQTICI)
!
!**********************************************************************C
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
! READ IN SPATIALLY AND/OR TEMPORALLY VARYING ICS (UNIT INWQICI).
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      DIMENSION XWQV(NWQVM)
      CHARACTER TITLE(3)*79, ICICONT*3
!
      OPEN(1,FILE=ICIFN,STATUS='UNKNOWN')
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
!
      IF(IWQTICI.EQ.0)THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      ENDIF
!
      WRITE(2,60)'* INITIAL CONDITIONS AT ', IWQTICI,' TH DAY FROM MODEL START'
!
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
      DO M=2,LA
        READ(1,84) I,J,K,(XWQV(NW),NW=1,NWQV)
        IF(IJCT(I,J).LT.1 .OR. IJCT(I,J).GT.8)THEN
          PRINT*, 'I, J, K, LINE# = ', I,J,K,M-1
          STOP 'ERROR!! INVALID (I,J) IN FILE 1'
        ENDIF
        L=LIJ(I,J)
        DO NW=1,NWQV
          WQV(L,K,NW)=XWQV(NW)
        ENDDO
        WRITE(2,84) I,J,K,(WQV(L,K,NW),NW=1,NWQV)
      ENDDO
!
!: WQCHLX=1/WQCHLX
!
      DO L=2,LA
        DO K=1,KC
        IZ=IWQZMAP(L,K)                                                               !hnr 10/2010
          WQCHL(L,K) = WQV(L,K,1)*WQCHLC(iz) + WQV(L,K,2)*WQCHLD(iz) + WQV(L,K,3)*WQCHLG(iz)                                                   !hnr 10/2010
          IF(IWQSRP.EQ.1)THEN
            O2WQ = MAX(WQV(L,K,19), 0.0)
            WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQV(L,K,20) )
            WQTAMP(L,K) = WQV(L,K,20) - WQTAMD
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*WQTAMP(L,K))
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*WQTAMP(L,K))
           ELSE IF(IWQSRP.EQ.2)THEN
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*SEDT(L,K))
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*SEDT(L,K))
           ELSE
            WQPO4D(L,K) = WQV(L,K,10)
            WQSAD(L,K)  = WQV(L,K,17)
          ENDIF
        ENDDO
      ENDDO
!
      READ(1,52) IWQTICI, ICICONT
      WRITE(2,52) IWQTICI, ICICONT
      IF(ICICONT.EQ.'END')THEN
        CLOSE(1)
        IWQICI = 0
      ENDIF
!
      CLOSE(2)
!
  999 FORMAT(1X)
   50 FORMAT(A79)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
   84 FORMAT(3I5, 21E12.4)
!
      RETURN
      END