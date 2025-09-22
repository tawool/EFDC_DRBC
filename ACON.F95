!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE ACON(ITVAL)

!**********************************************************************C
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
      PARAMETER (NJELM=2,NATDM=1)
!
!**********************************************************************C
!
!      DIMENSION ZAD(KCM,NATDM),  TAD(KCM,NATDM), UAGD(KCM,NATDM),
!     &         VAGD(KCM,NATDM), WAGD(KCM,NATDM),SALAD(KCM,NATDM),
!     &        TEMAD(KCM,NATDM),DYEAD(KCM,NATDM),SFLAD(KCM,NATDM),
!     &        TOXAD(KCM,NTXM,NATDM),
!     &        SEDAD(KCM,NSCM,NATDM),SNDAD(KCM,NSNM,NATDM)
!
!      DIMENSION TOXA(NTXM),SEDA(NSCM),SNDA(NSNM)
!
!      COMMON/ADATA/ ZAD,  TAD, UAGD,VAGD, WAGD,SALAD,
!     &        TEMAD,DYEAD,SFLAD,TOXAD,SEDAD,SNDAD
!
!      COMMON/AVAL/ NAZD,NATD,ZVAL,TVAL,UAG, VAG, WAG,SALA,
!     &            TEMA,DYEA,SFLA,TOXA,SEDA,SNDA
!
!**********************************************************************C
!
      IF(ZVAL.LT.ZAD(1,ITVAL))THEN
        UAG=UAGD(1,ITVAL)
        VAG=VAGD(1,ITVAL)
        WAG=WAGD(1,ITVAL)
        SALA=SALAD(1,ITVAL)
        TEMA=TEMAD(1,ITVAL)
!JHFIX-REMOVE        DYEA=SALAD(1,ITVAL)
!JHFIX-ADD THE FOLLOWING 1 LINE(S)
        DYEA=DYEAD(1,ITVAL)
        SFLA=SFLAD(1,ITVAL)
        DO NT=1,NTOX
         TOXA(NT)=TOXAD(1,NT,ITVAL)
        ENDDO
        DO NS=1,NSED
         SEDA(NS)=SEDAD(1,NS,ITVAL)
        ENDDO
        DO NX=1,NSND
         SNDA(NX)=SNDAD(1,NX,ITVAL)
        ENDDO
        RETURN
      ENDIF
!
      IF(ZVAL.GE.ZAD(NAZD,ITVAL))THEN
        UAG=UAGD(NAZD,ITVAL)
        VAG=VAGD(NAZD,ITVAL)
        WAG=WAGD(NAZD,ITVAL)
        SALA=SALAD(NAZD,ITVAL)
        TEMA=TEMAD(NAZD,ITVAL)
        DYEA=DYEAD(NAZD,ITVAL)
        SFLA=SFLAD(NAZD,ITVAL)
        DO NT=1,NTOX
         TOXA(NT)=TOXAD(NAZD,NT,ITVAL)
        ENDDO
        DO NS=1,NSED
         SEDA(NS)=SEDAD(NAZD,NS,ITVAL)
        ENDDO
        DO NX=1,NSND
         SNDA(NX)=SNDAD(NAZD,NX,ITVAL)
        ENDDO
        RETURN
      ENDIF
!
      NZ=1
 1000 CONTINUE
      NZP=NZ+1
      IF(ZVAL.GE.ZAD(NZ,ITVAL).AND.ZVAL.LT.ZAD(NZP,ITVAL))THEN
        DZI=1./(ZAD(NZP,ITVAL)-ZAD(NZ,ITVAL))
        WTNZ=DZI*(ZAD(NZP,ITVAL)-ZVAL)
        WTNZP=DZI*(ZVAL-ZAD(NZ,ITVAL))
        UAG=WTNZ*UAGD(NZ,ITVAL)+WTNZP*UAGD(NZP,ITVAL)
        VAG=WTNZ*VAGD(NZ,ITVAL)+WTNZP*VAGD(NZP,ITVAL)
        WAG=WTNZ*WAGD(NZ,ITVAL)+WTNZP*WAGD(NZP,ITVAL)
        SALA=WTNZ*SALAD(NZ,ITVAL)+WTNZP*SALAD(NZP,ITVAL)
        TEMA=WTNZ*TEMAD(NZ,ITVAL)+WTNZP*TEMAD(NZP,ITVAL)
        DYEA=WTNZ*DYEAD(NZ,ITVAL)+WTNZP*DYEAD(NZP,ITVAL)
        SFLA=WTNZ*SFLAD(NZ,ITVAL)+WTNZP*SFLAD(NZP,ITVAL)
        DO NT=1,NTOX
         TOXA(NT)=WTNZ*TOXAD(NZ,NT,ITVAL)+WTNZP*TOXAD(NZP,NT,ITVAL)
        ENDDO
        DO NS=1,NSED
         SEDA(NS)=WTNZ*SEDAD(NZ,NS,ITVAL)+WTNZP*SEDAD(NZP,NS,ITVAL)
        ENDDO
        DO NX=1,NSND
         SNDA(NX)=WTNZ*SNDAD(NZ,NX,ITVAL)+WTNZP*SNDAD(NZP,NX,ITVAL)
        ENDDO
        RETURN
       ELSE
        NZ=NZ+1
        GOTO 1000
      ENDIF
!
!**********************************************************************C
!
      RETURN
      END