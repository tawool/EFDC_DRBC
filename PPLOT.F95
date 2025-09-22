!
!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!  
      SUBROUTINE PPLOT (IPT)
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
!**********************************************************************C
!  
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!  
      DIMENSION BNDU(51),BNDL(51)
!
      CHARACTER BLANK,ASTER,LET1(51),LET2(51),CHARY(ICM,JCM)
!
      DATA BLANK/' '/
      DATA ASTER/'*'/
!
      DATA LET1/'A',' ','B',' ','C',' ','D',' ','E',' ','F',' ',
     &         'G',' ','H',' ','I',' ','J',' ','K',' ','L',' ',
     &         'M',' ','N',' ','O',' ','P',' ','Q',' ','R',' ',
     &         'S',' ','T',' ','U',' ','V',' ','W',' ','X',' ',
     &         'Y',' ','Z'/                                    
!
      DATA LET2/'A','A','B','B','C','C','D','D','E','E','F','F',
     &         'G','G','H','H','I','I','J','J','K','K','L','L',
     &         'M','M','N','N','O','O','P','P','Q','Q','R','R',
     &         'S','S','T','T','U','U','V','V','W','W','X','X',
     &         'Y','Y','Z'/                                    
!
      PMAX=-99999.
      PMIN= 99999.
!
      DO L=2,LA
      IF(PAM(L) .GT. PMAX) PMAX=PAM(L)
      IF(PAM(L) .LT. PMIN) PMIN=PAM(L)
      ENDDO
!
      RNBAN=FLOAT(NBAN)
      PINV=(PMAX-PMIN)/RNBAN
!
! **  SET BOUND ARRAYS FOR REFERENCE TABLE
!
      BNDU(1)=PMAX
      BNDL(1)=PMAX-PINV
!
      DO M=2,NBAN
      MM=M-1
      BNDU(M)=BNDL(MM)
      BNDL(M)=BNDU(M)-PINV
      ENDDO
!
      IF(IPT.EQ.1)THEN
       DO M=1,NBAN
       WRITE (7,10) BNDU(M),LET1(M),BNDL(M)
       ENDDO
      ELSE
       DO M=1,NBAN
       WRITE (7,10) BNDU(M),LET2(M),BNDL(M)
       ENDDO
      ENDIF
!
   10 FORMAT (5X,E12.4,5X,A1,5X,E12.4)
!
!     WRITE (7,11)
   11 FORMAT (////)
      WRITE(7,12)
   12 FORMAT(1H1)
!
! **  LOAD CHARACTER ARRAY
!
      DO J=1,JC
      DO I=1,IC
      IF(IJCT(I,J) .NE. 9) CHARY(I,J)=BLANK
      IF(IJCT(I,J) .EQ. 9) CHARY(I,J)=ASTER
      ENDDO
      ENDDO
!
      BNDU(1)=BNDU(1)+1.
      BNDL(NBAN)=BNDL(NBAN)-1.
!
      DO L=2,LA
      I=IL(L)
      J=JL(L)
       IF(IPT.EQ.1)THEN
        DO M=1,NBAN
        IF(PAM(L).LT.BNDU(M).AND.PAM(L).GE.BNDL(M)) CHARY(I,J)=LET1(M)
        ENDDO
       ELSE
        DO M=1,NBAN
        IF(PAM(L).LT.BNDU(M).AND.PAM(L).GE.BNDL(M)) CHARY(I,J)=LET2(M)
        ENDDO
       ENDIF
      ENDDO
!
      DO JJ=1,JC,120
      JS=JJ
      JE=JJ+119
      IF(JE.GT.JC) JE=JC
      WRITE(7,22)JS,JE
      DO I=1,IC
      WRITE (7,20) I,(CHARY(I,J),J=JS,JE)
      ENDDO 
      ENDDO
!
!  20 FORMAT (10X,80A1)
   20 FORMAT (1X,I3,2X,120A1)
   22 FORMAT(1H1,'VALUES FOR J=',I5,2X,'TO J=',I5,//)
!
      RETURN
      END