!
!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE BAL2T3B(IBALSTDT)
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
! 05/01/2002        john hamrick       05/01/2002       john hamrick
!  subroutine added for 2 time-level balances including sed,snd,tox
! 06/06/2002        john hamrick       06/06/2002       john hamrick
!  modified snd mass balance with respect to bed load outflow
! 06/07/2002        john hamrick       06/07/2002       john hamrick
!  added qdwaste to water mass balance
!----------------------------------------------------------------------C
!
! **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM,
! **  AND ENERGY BALANCES
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
      ELSE
        DELT=DTDYN
      END IF
!
!**********************************************************************C
!
!     DIMENSION CONT(LCM,KCM)
!
!**********************************************************************C
!
! **  ACCUMULATE INTERNAL SOURCES AND SINKS
!
!----------------------------------------------------------------------C
!
      IF(IBALSTDT.EQ.1)THEN
        DO L=2,LA
          WVOLOUT=WVOLOUT-DTSED*QMORPH(L)
          BVOLOUT=BVOLOUT+DTSED*QMORPH(L)
          VOLMORPH2T=VOLMORPH2T+DTSED*QMORPH(L)
!          WVOLOUT=WVOLOUT-DTSED*(QWTRBEDA(L)+QWTRBEDA1(L)
!    &                          +QWTRBED(L,KBT(L)))
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(ISTRAN(5).GE.1)THEN
!
      DO NT=1,NTOX
      M=MSVTOX(NT)
!SO   WRITE(8,*)'nt m ',NT,M
!
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=TOX(L,K,NT)
       ENDDO
      ENDDO
!
!
!  TOXBLB2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE
!  DUE TO BED LOAD TRANSPORT OUT OF DOMAIN
!
      IF(IBALSTDT.EQ.1)THEN
        IF(NSBDLDBC.GT.0) THEN
          TOXBLB2T(NT)=TOXBLB2T(NT)+DTSED*TOXBLB(NT)
        ENDIF
!
!  TOXFLUXW2T(NT) IS WATER COLUMN SIDE TOXIC FLUX DUE TO SUSPENDED LOAD
!    (POSITIVE INTO WATER COLUMN)
!  TOXFLUXB2T(NT) IS BED SIDE TOXIC FLUX DUE TO SUSPENDED LOAD (POSITIVE INTO WATER COLUMN)
!  TADFLUX2T(NT) IS PORE WATER ADVECTION+DIFFUSION FLUX (POSITIVE INTO WATER COLUMN)
!  TOXFBL2T(NT) IS NET TOXIC FLUX FROM BED ASSOCIATED WITH BED LOAD TRANSPORT
!    (SHOULD EQUAL TOXBLB2T(NT)
!
        DO L=2,LA
!        TOXFTMP=TOXF(L,0,NT)+TOXFB(L,KBT(L),NT)
!        TOXFLUXW2T(NT)=TOXFLUXW2T(NT)+DELT*DXYP(L)*TOXFTMP
!        TOXFLUXB2T(NT)=TOXFLUXB2T(NT)+DELT*DXYP(L)*(TOXFTMP
!     &     +TOXFBL(L,NT))
!        TADFLUX2T(NT)=TADFLUX2T(NT)+DELT*DXYP(L)*TADFLUX(L,NT)
          TOXFLUXW2T(NT)=TOXFLUXW2T(NT)+DTSED*DXYP(L)*TOXF(L,0,NT)
          TOXFLUXB2T(NT)=TOXFLUXB2T(NT)+DTSED*DXYP(L)*TOXFB(L,KBT(L),NT)
          TADFLUX2T(NT)=TADFLUX2T(NT)+DTSED*DXYP(L)*TADFLUX(L,NT)
!        TOXFBL2T(NT)=TOXFBL2T(NT)+DELT*DXYP(L)*TOXFBL(L,NT)
        ENDDO
!
        TOXFBL2T(NT)=TOXFBL2T(NT)+DTSED*TOXFBLT(NT)
!
        IF(ISBKERO.GE.1)THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            TOXFLUXB2T(NT)=TOXFLUXB2T(NT)+DTSED*DXYP(LBANK)*TOXFBEBKB(LBANK,NT)
            TOXFLUXB2T(NT)=TOXFLUXB2T(NT)+DTSED*DXYP(LCHAN)*TOXFBECHB(LCHAN,NT)
          ENDDO
        ENDIF
!
      ENDIF
!
      ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(ISTRAN(6).GE.1)THEN
!
      DO NSX=1,NSED
      M=MSVSED(NSX)
!
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=SED(L,K,NSX)
       ENDDO
      ENDDO
!
! SEDFLUX2T(NSX) IS IS NET COHESIVE MASS FLUX POSITIVE FROM BED
!   TO WATER COLUMN
!
      IF(IBALSTDT.EQ.1)THEN
!
        DO L=2,LA
          SEDFLUX2T(NSX)=SEDFLUX2T(NSX)+DTSED*DXYP(L)*SEDF(L,0,NSX)
        ENDDO
!
        IF(ISBKERO.GE.1)THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
             SEDFLUX2T(NSX)=SEDFLUX2T(NSX)+DTSED*DXYP(LBANK)*SEDFBEBKB(LBANK,NSX)
            SEDFLUX2T(NSX)=SEDFLUX2T(NSX)+DTSED*DXYP(LCHAN)*SEDFBECHB(LCHAN,NSX)
          ENDDO
        ENDIF
!
      ENDIF
!
      ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(ISTRAN(7).GE.1)THEN
!
      DO NSX=1,NSND
      M=MSVSND(NSX)
!
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=SND(L,K,NSX)
       ENDDO
      ENDDO
!
!  SBLOUT2T(NSX) IS NET NONCOHESIVE SEDIMENT MASS GOING OUT OF DOMAIN DUE
!  DUE TO BED LOAD TRANSPORT OUT OF DOMAIN
!
      IF(IBALSTDT.EQ.1)THEN
      IF(NSBDLDBC.GT.0) THEN
        DO NSB=1,NSBDLDBC
          LUTMP=LSBLBCU(NSB)
          LDTMP=LSBLBCD(NSB)
          IF(LDTMP.EQ.0) THEN
              SBLOUT2T(NSX)=SBLOUT2T(NSX)+DTSED*QSBDLDOT(LUTMP,NSX)
          ENDIF
        ENDDO
      ENDIF
      ENDIF
!
!  SNDFLUX2T(NSX) IS NET NONCOHESIVE SEDIMENT FLUX DUE TO SUSPENDED LOAD
!    (POSITIVE INTO WATER COLUMN)
!  SNDFBL2T(NSX) IS NET NONCOHESIVE SEDIMENT FLUX FROM BED ASSOCIATED WITH
!    BED LOAD TRANSPORT (SHOULD EQUAL SBLOUT2T(NSX))
!
!      TMPVAL=SNDFBL2T(NSX)
      IF(IBALSTDT.EQ.1)THEN
!
        DO L=2,LA
          SNDFLUX2T(NSX)=SNDFLUX2T(NSX)+DTSED*DXYP(L)*(SNDF(L,0,NSX)-SNDFBL(L,NSX))
          SNDFBL2T(NSX)=SNDFBL2T(NSX)+DTSED*DXYP(L)*SNDFBL(L,NSX)
        ENDDO
!
        IF(ISBKERO.GE.1)THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
             SNDFLUX2T(NSX)=SNDFLUX2T(NSX)+DTSED*DXYP(LBANK)*SNDFBEBKB(LBANK,NSX)
            SNDFLUX2T(NSX)=SNDFLUX2T(NSX)+DTSED*DXYP(LCHAN)*SNDFBECHB(LCHAN,NSX)
          ENDDO
        ENDIF
!
      ENDIF
!      TMPVAL=SNDFBL2T(NSX)-TMPVAL
!       WRITE(8,800)N,NSX,SNDFBL2T(NSX),TMPVAL
!
      ENDDO
!
      ENDIF
!
  800 FORMAT('N,NS,SNDFBL2T,DEL',2I5,2E14.5)
!
!----------------------------------------------------------------------C
!
      RETURN
      END