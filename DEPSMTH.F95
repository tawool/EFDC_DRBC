!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE DEPSMTH
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
      IF(WSMH.GT.0.) GOTO 1000
!
! **  SMOOTH BOTTOM ELEVATION
!
      WSMH=-WSMH
!
      DO NSM=1,NSHMAX
!
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9)THEN
        I=IL(L)
        J=JL(L)
        HTN=BELV(LNC(L))
        HTS=BELV(LSC(L))
        IF(IJCT(I,J+1).EQ.9) HTN=BELV(L)
        IF(IJCT(I,J-1).EQ.9) HTS=BELV(L)
        HTMP(L)=(1.-WSMH)*BELV(L)+0.5*WSMH*(HTN+HTS)
       ENDIF
      ENDDO
!
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9)THEN
        I=IL(L)
        J=JL(L)
        HTE=HTMP(L+1)
        HTW=HTMP(L-1)
        IF(IJCT(I+1,J).EQ.9) HTE=HTMP(L)
        IF(IJCT(I-1,J).EQ.9) HTW=HTMP(L)
        BELV(L)=(1.-WSMH)*HTMP(L)+0.5*WSMH*(HTE+HTW)
       ENDIF
      ENDDO
!
      ENDDO
!
! **  SMOOTH DEPTH
!
      DO NSM=1,NSHMAX
!
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9)THEN
        I=IL(L)
        J=JL(L)
        HTN=HMP(LNC(L))
        HTS=HMP(LSC(L))
        IF(IJCT(I,J+1).EQ.9) HTN=HMP(L)
        IF(IJCT(I,J-1).EQ.9) HTS=HMP(L)
        HTMP(L)=(1.-WSMH)*HMP(L)+0.5*WSMH*(HTN+HTS)
       ENDIF
      ENDDO
!
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9)THEN
        I=IL(L)
        J=JL(L)
        HTE=HTMP(L+1)
        HTW=HTMP(L-1)
        IF(IJCT(I+1,J).EQ.9) HTE=HTMP(L)
        IF(IJCT(I-1,J).EQ.9) HTW=HTMP(L)
        HMP(L)=(1.-WSMH)*HTMP(L)+0.5*WSMH*(HTE+HTW)
       ENDIF
      ENDDO
!
      ENDDO
!
      GOTO 2000
!
!**********************************************************************C
!
 1000 CONTINUE
!
! **  SMOOTH BOTTOM ELEVATION
!
      DO L=2,LA
       HTMP(L)=BELV(L)
      ENDDO
!
      DO NSM=1,NSHMAX
!
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9)THEN
        I=IL(L)
        J=JL(L)
        HTN=HTMP(LNC(L))
        HTS=HTMP(LSC(L))
        HTE=HTMP(L+1)
        HTW=HTMP(L-1)
        IF(IJCT(I  ,J+1).EQ.9) HTN=HTMP(L)
        IF(IJCT(I  ,J-1).EQ.9) HTS=HTMP(L)
        IF(IJCT(I+1,J  ).EQ.9) HTE=HTMP(L)
        IF(IJCT(I-1,J  ).EQ.9) HTW=HTMP(L)
        BELV(L)=(1.-WSMH)*HTMP(L)+0.25*WSMH*(HTN+HTS+HTE+HTW)
       ENDIF
      ENDDO
!
      DO L=2,LA
       HTMP(L)=BELV(L)
      ENDDO
!
      ENDDO
!
! **  SMOOTH DEPTH
!
      DO L=2,LA
       HTMP(L)=HMP(L)
      ENDDO
!
      DO NSM=1,NSHMAX
!
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9)THEN
        I=IL(L)
        J=JL(L)
        HTN=HTMP(LNC(L))
        HTS=HTMP(LSC(L))
        HTE=HTMP(L+1)
        HTW=HTMP(L-1)
        IF(IJCT(I  ,J+1).EQ.9) HTN=HTMP(L)
        IF(IJCT(I  ,J-1).EQ.9) HTS=HTMP(L)
        IF(IJCT(I+1,J  ).EQ.9) HTE=HTMP(L)
        IF(IJCT(I-1,J  ).EQ.9) HTW=HTMP(L)
        HMP(L)=(1.-WSMH)*HTMP(L)+0.25*WSMH*(HTN+HTS+HTE+HTW)
       ENDIF
      ENDDO
!
      DO L=2,LA
       HTMP(L)=HMP(L)
      ENDDO
!
      ENDDO
!
!**********************************************************************C
!
 2000 CONTINUE
!
       OPEN(1,FILE='NEWDXDY.INP',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='NEWDXDY.INP',STATUS='UNKNOWN')
       DO L=2,LA
        WRITE(1,339)IL(L),JL(L),DXP(L),DYP(L),HMP(L),BELV(L),ZBR(L)
       ENDDO
       CLOSE(1)
!
  339 FORMAT(1X,I5,2X,I5,2X,E12.4,2X,E12.4,2X,E12.4,2X,E12.4,2X,E12.4)
!
      RETURN
      END