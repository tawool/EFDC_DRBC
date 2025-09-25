!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE RSMICI(ISMTICI)
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
!
!**********************************************************************C
!
! READ IN SPATIALLY AND/OR TEMPORALLY VARYING ICS (UNIT INSMICI).
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      DIMENSION XSMPON(NSMGM),XSMPOP(NSMGM),XSMPOC(NSMGM)
      CHARACTER TITLE(3)*79, ICICONT*3
!
      OPEN(1,FILE='WQSDICI.INP',STATUS='OLD')
!
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
!
      IF(ISMTICI.EQ.0)THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      ENDIF
!
      WRITE(2,60)'* INITIAL CONDITIONS AT ', ISMTICI,' TH DAY FROM MODEL START'
!
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
!
      DO M=2,LA
!QUESTION        READ(INSMRST,90) I,J,(XSMPON(NW),NW=1,NSMG),
        READ(1,90) I,J,(XSMPON(NW),NW=1,NSMG),(XSMPOP(NW),NW=1,NSMG),(XSMPOC(NW),NW=1,NSMG),      &
                   XSM1NH4,XSM2NH4,XSM2NO3,XSM2PO4,XSM2H2S,XSMPSI,XSM2SI,XSMBST,XSMT
        WRITE(2,90) I,J,(XSMPON(NW),NW=1,NSMG),(XSMPOP(NW),NW=1,NSMG),(XSMPOC(NW),NW=1,NSMG),     &
                    XSM1NH4,XSM2NH4,XSM2NO3,XSM2PO4,XSM2H2S,XSMPSI,XSM2SI,XSMBST,XSMT
        IF(IJCT(I,J).LT.1 .OR. IJCT(I,J).GT.8)THEN
          PRINT*, 'I, J, LINE# = ', I,J,M-1
          STOP 'ERROR!! INVALID (I,J) IN FILE WQSDICI.INP'
        ENDIF
        L=LIJ(I,J)
        DO MM=1,NSMG
          SMPON(L,MM)=XSMPON(MM)
          SMPOP(L,MM)=XSMPOP(MM)
          SMPOC(L,MM)=XSMPOC(MM)
        ENDDO
        SM1NH4(L)=XSM1NH4
        SM2NH4(L)=XSM2NH4
        SM2NO3(L)=XSM2NO3
        SM2PO4(L)=XSM2PO4
        SM2H2S(L)=XSM2H2S
        SMPSI(L) =XSMPSI
        SM2SI(L) =XSM2SI
        SMBST(L) =XSMBST
        SMT(L)   =XSMT
      ENDDO
!
      READ(1,52) ISMTICI, ICICONT
      WRITE(2,52) ISMTICI, ICICONT
!
      IF(ICICONT.EQ.'END')THEN
        CLOSE(1)
        ISMICI = 0
      ENDIF
!
      CLOSE(2)
!
  999 FORMAT(1X)
   50 FORMAT(A79)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
   84 FORMAT(3I5, 20F8.4, F8.2)
   90 FORMAT(2I5, 25E14.4)
!
      RETURN
      END
