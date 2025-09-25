!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
     SUBROUTINE CALBED
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
!----------------------------------------------------------------------C
!
!**********************************************************************C
!
! **  SUBROUTINE CALSND CALCULATES NONCOHESIVER SEDIMENT SETTLING,
! **  DEPOSITION AND RESUSPENSION AND IS CALLED FOR SSEDTOX
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!zzdiff      COMMON/PMC/CCSHEAR(LCM),DSTAR(NSTM),BDLDFACTOR(NSTM)  
!zzdiff      REAL*4 TMPVAL,TMPVAL1,TMPVAL2  
!zzdiff      REAL*4 DSTAR  
!zzdiff      REAL*4 CCUSTAR(LCM)  
!
      COMMON/SSEDTOX1/ CTMPDRY(LCM),CSHIELDS50(LCM),                                              & 
                     USTAR(LCM),UCELLCTR(LCM),VCELLCTR(LCM),                                      & 
                     QWBDTOP(LCM),QSBDTOP(LCM),ZETATOP(LCM),                                      & 
                     ZBEDGT(LCM),QSBDLDP(LCM),RBPSBL(LCM),                                        & 
                     TOXBBALO(LCM),TOXBBALN(LCM),TOXWBALO(LCM),                                   & 
                     ROUSE(LCM),ZEQ(LCM),ZEQI(LCM),ZEQD(LCM),                                     & 
                     ZEQDI(LCM),SNDEQ(LCM),SNDEQB(LCM),                                           & 
                     TOXWBALN(LCM),FACSUSL(LCM),FACBEDL(LCM),                                     & 
                     PEXP(LCM,NSNM),PHID(LCM,NSNM),                                               & 
                     SEDFPA(LCM,NSCM),SNDFPA(LCM,NSNM)
!
      COMMON/SSEDTOX1A/ CBEDTOTAL(LCM),QCELLCTR(LCM),HGDH(LCM),                                   &
                       FRACCOH(LCM,KBM),FRACNON(LCM,KBM),USTARSED(LCM),                           &
                       USTARSND(LCM),QWATPA(LCM),QSSDPA(LCM)
!
      COMMON/SSEDTOX2/ CDECAYW(LCM,KCM),CDECAYB(LCM,KBM),STRSE(LCM,KBM),                          &
                      HYDCN(LCM,KBM),COEFK(LCM,KBM),COEFSK(LCM,KBM),                              &
                      PRESE(LCM,KBM),PRESH(LCM,KBM),PREST(LCM,KBM),                               &
                      STRST(LCM,KBM),DZBTR1(LCM,KBM),SEDDIA50(LCM,KBM),                           &
                      SEDBALL(LCM,KBM),ZBEDC(LCM,KBM),SGSM1(LCM,KBM),                             &
                      DSTRSE(LCM,KBM),DZBTR(LCM,KBM),STRSEM(LCM,KBM),                             &
                      SEDDIA90(LCM,KBM),SEDGEOSTD(LCM,KBM),                                       &
                      SEDDIAGS(LCM,KBM),ZOTOTAL(LCM),ZOGRAIN(LCM)
!
      COMMON/SSEDTOX3/ ALOW(LCM,KBM+1),BMNN(LCM,KBM+1),CUPP(LCM,KBM+1),                           &
                      RRHS(LCM,KBM+1),TOXTMP(LCM,KBM+1),                                          &
                      GAMTMP(LCM,KBM+1),ACOEF(LCM,0:KBM),                                         &
                      QCOEF(LCM,0:KBM),ZBEDG(LCM,0:KBM)
!
      COMMON/SSEDTOX4/ WSETA(LCM,0:KSM,NSTM),SEDS(LCM,KCM,NSCM),                                  &
                      SNDS(LCM,KCM,NSNM),TOXS(LCM,KCM,NTXM),                                      &
                      SEDBS(LCM,KBM,NSCM),SNDBS(LCM,KBM,NSNM),                                    &
                      TOXBS(LCM,KBM,NTXM)
!
      COMMON/SSEDTOX5/ NSP2(NTXM),DERRB(KBM),STRESSS(0:KSM)
!
      COMMON/SSEDTOX6/ DELT,DELTI,DSEDGMM,FOURDPI,SEDMDGM,S2TL,S3TL,CORDT,DIASED,GPDIASED,BEDEX
!
      COMMON/SSEDTOX7/ ISTL,IS2TL,ISUD
!
      DIMENSION SNDHYDCN(LCM,KBM)
!
!**********************************************************************C
!
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
        HBEDMIN=1.E-4
        IF(ISTRAN(7).GE.1)THEN
          HBEDMIN=MAX(HBEDMIN,SNDDMX)
        END IF
!
!----------------------------------------------------------------------C
!
! ** CONSTANT POROSITY BED
!
        IF(IBMECH.EQ.0)THEN
!
! ** UPDATE TOP LAYER THICKNESS TO MAINTIAN CONSTANT POROSITY
!
          VOIDCON1=BEDPORC/(1.-BEDPORC)

          DO L=2,LA
            K=KBT(L)
            HBEDTMP=(1.+VOIDCON1)*HBED(L,K)/(1.+VDRBED(L,K))
            TMPVALO=VDRBED(L,K)*HBED(L,K)/(1.+VDRBED(L,K))
            TMPVALN=VOIDCON1*HBEDTMP/(1.+VOIDCON1)
            QWBDTOP(L)=DELTI*(TMPVALO-TMPVALN)
            HBED(L,K)=HBEDTMP
            QWTRBED(L,K)=QWBDTOP(L)+QGW(L)/DXYP(L)
          ENDDO
!
          DO K=0,KBT(L)-1
            DO L=2,LA
              QWTRBED(L,K)=QGW(L)/DXYP(L)
            ENDDO
          END DO
!
! ** ADD OR REMOVE LAYERS
!
!          CALL CALBLAY
!
        ENDIF
!
!----------------------------------------------------------------------C
!
! ** SIMPLE CONSOLIDATING BED
!
        IF(IBMECH.EQ.1)THEN
!
! ** DETERMINE TIME DIFFERENCE AND UPDATE VOID RATIO
!
! BEGIN JMH FIXED IBMECH.EQ.1  OPTION 12/30/02
!
! **  IF SEDVRDT.GT.0.0001 CONSOLIDATE TO SEDVRM (THE MINIMUM VOID RATIO)
!
          IF(SEDVRDT.GT.0.0001)THEN
              TMPEXP=EXP(-DELT/SEDVRDT)
            DO K=1,KB
              DO L=2,LA
                IF(K.LE.KBT(L))THEN
                  VDRBED1(L,K)=VDRBED(L,K)
                  HBED1(L,K)=HBED(L,K)
                  VDRBED(L,K)=SEDVDRM+(VDRBED1(L,K)-SEDVDRM)*TMPEXP
                  TMPTOP=1.+VDRBED(L,K)
                  TMPBOT=1.+VDRBED1(L,K)
                  HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
                ENDIF
              ENDDO
            ENDDO
          ENDIF
!
! **  IF SEDVRDT.GT.0.0001 CONSOLIDATE TO SEDVRM INSTANTANEOUSLY
!
          IF(SEDVRDT.GE.0.0.AND.SEDVRDT.LE.0.0001)THEN
              TMPEXP=0.0
            DO K=1,KB
              DO L=2,LA
                IF(K.LE.KBT(L))THEN
                  VDRBED1(L,K)=VDRBED(L,K)
                  HBED1(L,K)=HBED(L,K)
                  VDRBED(L,K)=SEDVDRM
                  TMPTOP=1.+VDRBED(L,K)
                  TMPBOT=1.+VDRBED1(L,K)
                  HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
                ENDIF
              ENDDO
            ENDDO
          ENDIF
!
!
! **  IF SEDVRDT.LT.0.0 MAINTAIN INITIAL VOID RATIO (SAVED IN VDRBED2)
!
!          IF(SEDVRDT.LT.0.0) THEN
!              TMPEXP=1.0
!           DO K=1,KB
!              DO L=2,LA
!                IF(K.LE.KBT(L))THEN
!                  VDRBED1(L,K)=VDRBED(L,K)
!                  HBED1(L,K)=HBED(L,K)
!                  VDRBED(L,K)=VDRBED2(L,K)
!                  TMPTOP=1.+VDRBED(L,K)
!                  TMPBOT=1.+VDRBED1(L,K)
!                  HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
!                ENDIF
!              ENDDO
!            ENDDO
!          ENDIF
!
          IF(SEDVRDT.LT.0.0) THEN
              TMPEXP=1.0
!            DO K=1,KB
              DO L=2,LA
                K=KBT(L)
                VDRBED1(L,K)=VDRBED(L,K)
                HBED1(L,K)=HBED(L,K)
                VDRBED(L,K)=VDRBED2(L,K)
!                IF(VDRBED(L,K).NE.VDRBED1(L,K))THEN
                  TMPTOP=1.+VDRBED(L,K)
                  TMPBOT=1.+VDRBED1(L,K)
                  HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
!               ENDIF
              ENDDO
              DO L=2,LA
                K=KBT(L)-1
                IF(K.GT.0)THEN
                  VDRBED1(L,K)=VDRBED(L,K)
                  HBED1(L,K)=HBED(L,K)
                  VDRBED(L,K)=VDRBED2(L,K)
!                  IF(VDRBED(L,K).NE.VDRBED1(L,K))THEN
                    TMPTOP=1.+VDRBED(L,K)
                    TMPBOT=1.+VDRBED1(L,K)
                    HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
!                 ENDIF
                ENDIF
              ENDDO
!            ENDDO
          ENDIF
!
! END JMH FIXED IBMECH.EQ.1  OPTION 12/30/02
!
! ** UPDATE POROSITY
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
                PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
              ENDIF
            ENDDO
          ENDDO
!
! ** UPDATE PORE WATER FLOWS
!
          DO L=2,LA
            QWTRBED(L,0)=QGW(L)/DXYP(L)
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                TMPVAL=HBED(L,K)/(1.+VDRBED(L,K))
                QWTRBED(L,K)=QWTRBED(L,K-1)-DELTI*TMPVAL*(VDRBED(L,K)-VDRBED1(L,K))
              ENDIF
            ENDDO
          ENDDO
!
! ** ADD OR REMOVE LAYERS
!
!          CALL CALBLAY
!
        ENDIF
!
!----------------------------------------------------------------------C
!
! ** CONSOLIDATE TO VOID RATIO/POROSITY DEFINED AS FUNCTION OF
!    COHESIVE SEDIMENT FRACTION, FRACCOH(L,K)
!
        IF(IBMECH.EQ.99)THEN
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                VDRBED1(L,K)=VDRBED(L,K)
                HBED1(L,K)=HBED(L,K)
                POROTMP=BMECH1*(FRACCOH(L,K)**BMECH2)+BMECH3
                VDRBED(L,K)=POROTMP/(1.-POROTMP)
                TMPTOP=1.+VDRBED(L,K)
                TMPBOT=1.+VDRBED1(L,K)
                HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
              ENDIF
            ENDDO
          ENDDO
!
! ** UPDATE POROSITY
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
                PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
              ENDIF
            ENDDO
          ENDDO
!
! ** UPDATE PORE WATER FLOWS
!
          DO L=2,LA
            QWTRBED(L,0)=QGW(L)/DXYP(L)
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                TMPVAL=HBED(L,K)/(1.+VDRBED(L,K))
                QWTRBED(L,K)=QWTRBED(L,K-1)-DELTI*TMPVAL*(VDRBED(L,K)-VDRBED1(L,K))
              ENDIF
            ENDDO
          ENDDO
!
        ENDIF
!
!----------------------------------------------------------------------C
!
! ** SET FLAG FOR FINITE STRAIN CONSOLIDATION
!
        IFLAG=0
        IF(ISTRAN(6).GE.1)IFLAG=IFLAG+1
        IF(ISTRAN(7).GE.1)IFLAG=IFLAG+1
!
!----------------------------------------------------------------------C
!
! ** FINITE STRAIN CONSOLIDATING HOMOGENEOUS BED
!
        IF(IBMECH.GE.2.AND.IBMECH.LT.98)THEN
        IF(IFLAG.EQ.1)THEN
!
          WDENKGM3=1.E3
          WDENGMM3=1.E6
!
! ++  SET PHYSICAL VERTICAL COORDINATES OF THE BED
!
!     ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
!     ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
!
          DO L=2,LA
            ZBEDG(L,0)=ZELBEDA(L)
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
              ELSE  
                ZBEDG(L,K)=HBED(L,K)  
              ENDIF
            ENDDO
          ENDDO
!
          DO L=2,LA
            ZBEDGT(L)=ZBEDG(L,KBT(L))
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
              ELSE  
                ZBEDC(L,K)=0.5*ZBEDG(L,K)  
              ENDIF
            ENDDO
          ENDDO
!
! ++  CALCULATE TRANSFORMED THICKNESS OF BED LAYERS
!     DZBTR = TRANSFORMED THICKNESS OF BED LAYER
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
                DZBTR1(L,K)=HBED1(L,K)/(1.+VDRBED1(L,K))
              ENDIF
            ENDDO
          ENDDO
!
! ++  CALCULATE WATER SPECIFIC WEIGHT NORMALIZED
!       EFFECTIVE STRESS USING FSTRSE
!     CALCULATE DERIVATIVE OF EFFECTIVE STRESS 
!       WITH RESPECT TO VOID RATIO, DSTRSE USING
!       FUNCTION FDSTRSE
!     CALCULATE HYDRAULIC CONDUCTIVITY DIVIED BY (1+VOID),
!       HYDCN =USING FUNCTION FHYDCN 
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                STRSE(L,K)=FSTRSE(VDRBED(L,K),BMECH1,BMECH2,BMECH3)
                DSTRSE(L,K)=FDSTRSE(VDRBED(L,K),BMECH1,BMECH2,BMECH3)
                HYDCN(L,K)=FHYDCN(VDRBED(L,K),BMECH4,BMECH5,BMECH6,IBMECHK)
              ENDIF
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              KBTM1=KBT(L)-1
              IF(K.LE.KBTM1)THEN
                TMPVAL=( DZBTR(L,K)/HYDCN(L,K) )+( DZBTR(L,K+1)/HYDCN(L,K+1) )
                COEFK(L,K)=(DZBTR(L,K)+DZBTR(L,K+1))/TMPVAL
! fix 0602                STTMP=DSTRSE(L,K)*HYDCN(L,K)
! fix 0602                STTMPP=DSTRSE(L,K+1)*HYDCN(L,K+1)
! fix 0602                TMPVAL=( DZBTR(L,K)/STTMP )
! fix 0602     &          +( DZBTR(L,K+1)/STTMPP )
! fix 0602                COEFSK(L,K)=(DZBTR(L,K)+DZBTR(L,K+1))/TMPVAL
                DSTRESET=(DZBTR(L,K)*DSTRSE(L,K+1)+DZBTR(L,K+1)*DSTRSE(L,K))/(DZBTR(L,K)+DZBTR(L,K+1))
                COEFSK(L,K)=DSTRESET*COEFK(L,K)
              ENDIF
              IF(K.EQ.KBT(L))THEN  
                COEFK(L,K)=HYDCN(L,K)
                COEFSK(L,K)=DSTRSE(L,K)*HYDCN(L,K)
              ENDIF
            ENDDO
          ENDDO
!
! ++  CALCULATE PRESSURE COMPONENTS
!
          DO K=1,KB
            DO L=2,LA
              SGSM1(L,K)=0.
            ENDDO
          ENDDO
!
          IF(ISTRAN(6).GT.0)THEN
            DO NS=1,NSED
              DO K=1,KB
                DO L=2,LA
                  IF(K.LE.KBT(L))THEN
                    SGSM1(L,K)=SGSM1(L,K)+SSG(NS)*VFRBED(L,K,NS)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(ISTRAN(7).GT.0)THEN
            DO NX=1,NSND
              NS=NSED+NX
              DO K=1,KB
                DO L=2,LA
                  IF(K.LE.KBT(L))THEN
                    SGSM1(L,K)=SGSM1(L,K)+SSG(NS)*VFRBED(L,K,NS)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                SGSM1(L,K)=SGSM1(L,K)-1.
              ENDIF
            ENDDO
          ENDDO
!
! ++  NEW IMPLICIT CONSOLIDATION SOLUTION BEGINS HERE
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LT.KBT(L))THEN
                ACOEF(L,K)=2.*COEFSK(L,K)/(DZBTR(L,K)+DZBTR(L,K+1))
                TMPVALK=DZBTR(L,K)*SGSM1(L,K)
                TMPVALKP=DZBTR(L,K+1)*SGSM1(L,K+1)
                QCOEF(L,K)=(TMPVALK+TMPVALKP)*COEFK(L,K)/(DZBTR(L,K)+DZBTR(L,K+1))
              ELSE
                ACOEF(L,K)=0.0
                QCOEF(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
          DO L=2,LA
            K=KBT(L)
            ACOEF(L,K)=2.*COEFSK(L,K)/DZBTR(L,K)
            QCOEF(L,K)=SGSM1(L,K)*COEFK(L,K)
            QCOEF(L,0)=QGW(L)*DXYIP(L)
            QWTRBED(L,0)=QGW(L)*DXYIP(L)
          ENDDO
!
          DO L=2,LA
            ALOW(L,1)=0.
            CUPP(L,KBT(L))=0.
            DO K=1,KBT(L)-1
              CUPP(L,K)=-DELT*ACOEF(L,K)/DZBTR(L,K)
            ENDDO
            DO K=2,KBT(L)
              ALOW(L,K)=-DELT*ACOEF(L,K-1)/DZBTR(L,K)
            ENDDO
            DO K=1,KBT(L)
              BMNN(L,K)=1.0-ALOW(L,K)-CUPP(L,K)
            ENDDO
            K=KBT(L)
            BMNN(L,K)=BMNN(L,K)+DELT*ACOEF(L,K)/DZBTR(L,K)
            DO K=1,KBT(L)
              RRHS(L,K)=VDRBED(L,K)+DELT*(QCOEF(L,K-1)-QCOEF(L,K))/DZBTR(L,K)
              VDRBED1(L,K)=VDRBED(L,K)
              HBED1(L,K)=HBED(L,K)
            ENDDO
          ENDDO
!
          IF(IBMECH.EQ.2) THEN
            DO L=2,LA
              K=KBT(L)
!              RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*SEDVDRD/DZBTR(L,K)
              RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*VDRDEPO(1)/DZBTR(L,K)
            ENDDO
          ENDIF
!
          IF(IBMECH.EQ.3) THEN
            DO L=2,LA
              K=KBT(L)
              RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*(VDRBED(L,K)+(STRSE(L,K)/DSTRSE(L,K)))/DZBTR(L,K)
            ENDDO
          ENDIF
!
!          DO L=2,LA
!           DO K=1,KBT(L)
!            WRITE(8,*)K,ALOW(L,K),BMNN(L,K),CUPP(L,K),RRHS(L,K)
!           END DO
!          END DO
!
          DO L=2,LA
            BETTMP=BMNN(L,1)
            TOXTMP(L,1)=RRHS(L,1)/BETTMP
            DO KK=2,KBT(L)
              GAMTMP(L,KK)=CUPP(L,KK-1)/BETTMP
              BETTMP=BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
              TOXTMP(L,KK)=(RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))/BETTMP
            ENDDO
            DO KK=KBT(L)-1,1,-1
               TOXTMP(L,KK)=TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LT.KBT(L))THEN
                QWTRBED(L,K)=-ACOEF(L,K)*(TOXTMP(L,K+1)-TOXTMP(L,K))+QCOEF(L,K)
              ELSE
                QWTRBED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
          IF(IBMECH.EQ.2) THEN
            DO L=2,LA
              K=KBT(L)
              QWTRBED(L,K)=-ACOEF(L,K)*(SEDVDRD-TOXTMP(L,K))+QCOEF(L,K)
            ENDDO
          ELSE
            DO L=2,LA
              K=KBT(L)
              QWTRBED(L,K)=-ACOEF(L,K)*(VDRBED(L,K)+(STRSE(L,K)/DSTRSE(L,K))-TOXTMP(L,K))+QCOEF(L,K)
            ENDDO
          END IF
!
! ++  CALCULATE VOID RATIOS  
!     VDRBED =  VOID RATIO OF BED LAYER 
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                VDRBED(L,K)=VDRBED1(L,K)-DELT*(QWTRBED(L,K)-QWTRBED(L,K-1))/DZBTR1(L,K)
              ELSE
                VDRBED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
! ++  UPDATE LAYER THICKNESS
!     HBED = BED LAYER THICKNESS
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                HBED(L,K)=HBED1(L,K)*(1.+VDRBED(L,K))/(1.+VDRBED1(L,K))
              ELSE
                HBED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
          DO L=2,LA
            QWTRBED(L,0)=QGW(L)/DXYP(L)
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                QWTRBED(L,K)=QWTRBED(L,K-1)-DELTI*(HBED(L,K)-HBED1(L,K))
              ENDIF
            ENDDO
          ENDDO
!
!     ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
!     ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
              ELSE  
                ZBEDG(L,K)=HBED(L,K)  
              ENDIF
            ENDDO
          ENDDO
!
          DO L=2,LA
            ZBEDGT(L)=ZBEDG(L,KBT(L))
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
              ELSE  
                ZBEDC(L,K)=0.5*ZBEDG(L,K)  
              ENDIF
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
              ENDIF
            ENDDO
          ENDDO
!
! ** UPDATE POROSITY
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
                PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
              ELSE
                PORBED(L,K)=0.0
                PORBED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
! ** ADD OR REMOVE LAYERS
!
!          CALL CALBLAY
!
        ENDIF
        ENDIF
!
!----------------------------------------------------------------------C
!
! ** FINITE STRAIN CONSOLIDATING NON-HOMOGENEOUS BED
!
        IF(IBMECH.GE.2.AND.IBMECH.LT.98)THEN
        IF(IFLAG.EQ.2)THEN
!
          WDENKGM3=1.E3
          WDENGMM3=1.E6
!
! ++  SET PHYSICAL VERTICAL COORDINATES OF THE BED
!
!     ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
!     ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
!
          DO L=2,LA
            ZBEDG(L,0)=ZELBEDA(L)
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
              ELSE  
                ZBEDG(L,K)=HBED(L,K)  
              ENDIF
            ENDDO
          ENDDO
!
          DO L=2,LA
            ZBEDGT(L)=ZBEDG(L,KBT(L))
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
              ELSE  
                ZBEDC(L,K)=0.5*ZBEDG(L,K)  
              ENDIF
            ENDDO
          ENDDO
!
! ++  CALCULATE TRANSFORMED THICKNESS OF BED LAYERS
!     DZBTR = TRANSFORMED THICKNESS OF BED LAYER
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
                DZBTR1(L,K)=HBED1(L,K)/(1.+VDRBED1(L,K))
              ENDIF
            ENDDO
          ENDDO
!
! ++  CALCULATE WATER SPECIFIC WEIGHT NORMALIZED
!       EFFECTIVE STRESS USING FSTRSE
!     CALCULATE DERIVATIVE OF EFFECTIVE STRESS 
!       WITH RESPECT TO VOID RATIO, DSTRSE USING
!       FUNCTION FDSTRSE
!     CALCULATE HYDRAULIC CONDUCTIVITY DIVIED BY (1+VOID),
!       HYDCN =USING FUNCTION FHYDCN 
!
! ++  NONCOHESIVE HYDRAULIC CONDUCTIVITY BASED ON R. R. RUMER, CHAP 3,
! ++  EQ 13 AND 14, IN 'FLOW THROUGH POROUS MEDIA' ED. R. J. M. DE WIEST, 
! ++  ACADEMIC PRESS, 1969
!
!     KH=2854*(n**2)*(d**2)/(1-n) = 2854*(e**2)*(d**2)/(1+e)
!     for KH in m/s and d in meters
!
 
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                STRSE(L,K)=FSTRSE(VDRBEDSED(L,K),BMECH1,BMECH2,BMECH3)
                DSTRSE(L,K)=FDSTRSE(VDRBEDSED(L,K),BMECH1,BMECH2,BMECH3)
                HYDCN(L,K)=FHYDCN(VDRBEDSED(L,K),BMECH4,BMECH5,BMECH6,IBMECHK)
                TMPVAL=VDRBEDSND(L,K)/(1.+VDRBEDSND(L,K))
                SNDHYDCN(L,K)=2854.0*TMPVAL*TMPVAL*(SEDDIA50(L,K)**2)
              ENDIF
            ENDDO
          ENDDO
! 
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                COEFK(L,K)=(FRACCOH(L,K)/HYDCN(L,K))+(FRACNON(L,K)/SNDHYDCN(L,K))
               ENDIF
            ENDDO
          ENDDO
! 
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                HYDCN(L,K)=1./COEFK(L,K)
              ENDIF
            ENDDO
          ENDDO
! 
          DO K=1,KB
            DO L=2,LA
              COEFK(L,K)=0.
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              KBTM1=KBT(L)-1
              IF(K.LE.KBTM1)THEN
                TMPVAL=( DZBTR(L,K)/HYDCN(L,K))+( DZBTR(L,K+1)/HYDCN(L,K+1) )
                COEFK(L,K)=(DZBTR(L,K)+DZBTR(L,K+1))/TMPVAL
! fix 0602                STTMP=DSTRSE(L,K)*HYDCN(L,K)
! fix 0602                STTMPP=DSTRSE(L,K+1)*HYDCN(L,K+1)
! fix 0602                TMPVAL=( DZBTR(L,K)/STTMP )
! fix 0602     &          +( DZBTR(L,K+1)/STTMPP )
! fix 0602                COEFSK(L,K)=(DZBTR(L,K)+DZBTR(L,K+1))/TMPVAL
                DSTRESET=(DZBTR(L,K)*DSTRSE(L,K+1)+DZBTR(L,K+1)*DSTRSE(L,K))/(DZBTR(L,K)+DZBTR(L,K+1))
                COEFSK(L,K)=DSTRESET*COEFK(L,K)
              ENDIF
              IF(K.EQ.KBT(L))THEN  
                COEFK(L,K)=HYDCN(L,K)
                COEFSK(L,K)=DSTRSE(L,K)*HYDCN(L,K)
              ENDIF
            ENDDO
          ENDDO
!
! ++  CALCULATE PRESSURE COMPONENTS
!
          DO K=1,KB
            DO L=2,LA
              SGSM1(L,K)=0.
            ENDDO
          ENDDO
!
          IF(ISTRAN(6).GT.0)THEN
            DO NS=1,NSED
              DO K=1,KB
                DO L=2,LA
                  IF(K.LE.KBT(L))THEN
                    SGSM1(L,K)=SGSM1(L,K)+SSG(NS)*VFRBED(L,K,NS)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(ISTRAN(7).GT.0)THEN
            DO NX=1,NSND
              NS=NSED+NX
              DO K=1,KB
                DO L=2,LA
                  IF(K.LE.KBT(L))THEN
                    SGSM1(L,K)=SGSM1(L,K)+SSG(NS)*VFRBED(L,K,NS)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                SGSM1(L,K)=SGSM1(L,K)-1.
              ENDIF
            ENDDO
          ENDDO
!
! ++  NEW IMPLICIT CONSOLIDATION SOLUTION BEGINS HERE
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LT.KBT(L))THEN
                ACOEF(L,K)=2.*COEFSK(L,K)/(DZBTR(L,K)+DZBTR(L,K+1))
                TMPVALK=DZBTR(L,K)*SGSM1(L,K)
                TMPVALKP=DZBTR(L,K+1)*SGSM1(L,K+1)
                QCOEF(L,K)=(TMPVALK+TMPVALKP)*COEFK(L,K)/(DZBTR(L,K)+DZBTR(L,K+1))
              ELSE
                ACOEF(L,K)=0.0
                QCOEF(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
          DO L=2,LA
            K=KBT(L)
            ACOEF(L,K)=2.*COEFSK(L,K)/DZBTR(L,K)
            QCOEF(L,K)=SGSM1(L,K)*COEFK(L,K)
            QCOEF(L,0)=QGW(L)*DXYIP(L)
            QWTRBED(L,0)=QGW(L)*DXYIP(L)
          ENDDO
!
          DO L=2,LA
            ALOW(L,1)=0.
            CUPP(L,KBT(L))=0.
            DO K=1,KBT(L)-1
              CUPP(L,K)=-DELT*ACOEF(L,K)/DZBTR(L,K)
            ENDDO
            DO K=2,KBT(L)
              ALOW(L,K)=-DELT*ACOEF(L,K-1)/DZBTR(L,K)
            ENDDO
            DO K=1,KBT(L)
              BMNN(L,K)=FRACCOH(L,K)-ALOW(L,K)-CUPP(L,K)
            ENDDO
            K=KBT(L)
            BMNN(L,K)=BMNN(L,K)+DELT*ACOEF(L,K)/DZBTR(L,K)
            DO K=1,KBT(L)
              RRHS(L,K)=FRACCOH(L,K)*VDRBEDSED(L,K)+DELT*(QCOEF(L,K-1)-QCOEF(L,K))/DZBTR(L,K)
              VDRBED1(L,K)=VDRBED(L,K)
              HBED1(L,K)=HBED(L,K)
            ENDDO
          ENDDO
!
          IF(IBMECH.EQ.2) THEN
            DO L=2,LA
              K=KBT(L)
              RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*VDRDEPO(1)/DZBTR(L,K)
            ENDDO
          ENDIF
!
          IF(IBMECH.EQ.3) THEN
            DO L=2,LA
              K=KBT(L)
              RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*(VDRBEDSED(L,K)+(STRSE(L,K)/DSTRSE(L,K)))/DZBTR(L,K)
            ENDDO
          ENDIF
!
!          DO L=2,LA
!           DO K=1,KBT(L)
!            WRITE(8,*)K,ALOW(L,K),BMNN(L,K),CUPP(L,K),RRHS(L,K)
!           END DO
!          END DO
!
          DO L=2,LA
            BETTMP=BMNN(L,1)
            TOXTMP(L,1)=RRHS(L,1)/BETTMP
            DO KK=2,KBT(L)
              GAMTMP(L,KK)=CUPP(L,KK-1)/BETTMP
              BETTMP=BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
              TOXTMP(L,KK)=(RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))/BETTMP
            ENDDO
            DO KK=KBT(L)-1,1,-1
               TOXTMP(L,KK)=TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LT.KBT(L))THEN
                QWTRBED(L,K)=-ACOEF(L,K)*(TOXTMP(L,K+1)-TOXTMP(L,K))+QCOEF(L,K)
              ELSE
                QWTRBED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
          IF(IBMECH.EQ.2) THEN
            DO L=2,LA
              K=KBT(L)
              QWTRBED(L,K)=-ACOEF(L,K)*(VDRDEPO(1)-TOXTMP(L,K))+QCOEF(L,K)
            ENDDO
          ELSE
            DO L=2,LA
              K=KBT(L)
              QWTRBED(L,K)=-ACOEF(L,K)*(VDRBEDSED(L,K)+(STRSE(L,K)/DSTRSE(L,K))-TOXTMP(L,K))+QCOEF(L,K)
            ENDDO
          END IF
!
! ++  CALCULATE VOID RATIOS  
!     VDRBED =  VOID RATIO OF BED LAYER 
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
	          IF(FRACCOH(L,K).GT.0.0)THEN
                  VDRBEDSED(L,K)=VDRBEDSED(L,K)-DELT*(QWTRBED(L,K)-QWTRBED(L,K-1))/(FRACCOH(L,K)*DZBTR1(L,K))
                  VDRBEDSED(L,K)=MAX(VDRBEDSED(L,K),VDRBEDSND(L,K))
                ELSE
                  VDRBEDSED(L,K)=0.0
                ENDIF
              ELSE
                VDRBEDSED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
! ++  UPDATE VOID RATIO
! 
          DO K=1,KB
           DO L=2,LA
              IF(K.LE.KBT(L))THEN
                VDRBED(L,K)=FRACCOH(L,K)*VDRBEDSED(L,K)+FRACNON(L,K)*VDRBEDSND(L,K)
              ELSE
                VDRBED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
! 
! ++  UPDATE LAYER THICKNESS
!     HBED = BED LAYER THICKNESS
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                HBED(L,K)=HBED1(L,K)*(1.+VDRBED(L,K))/(1.+VDRBED1(L,K))
              ELSE
                HBED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
          DO L=2,LA
            QWTRBED(L,0)=QGW(L)/DXYP(L)
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                QWTRBED(L,K)=QWTRBED(L,K-1)-DELTI*(HBED(L,K)-HBED1(L,K))
              ENDIF
            ENDDO
          ENDDO
!
!     ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
!     ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
              ELSE  
                ZBEDG(L,K)=HBED(L,K)  
              ENDIF
            ENDDO
          ENDDO
!
          DO L=2,LA
            ZBEDGT(L)=ZBEDG(L,KBT(L))
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
              ELSE  
                ZBEDC(L,K)=0.5*ZBEDG(L,K)  
              ENDIF
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
              ENDIF
            ENDDO
          ENDDO
!
! ** UPDATE POROSITY
!
          DO K=1,KB
            DO L=2,LA
              IF(K.LE.KBT(L))THEN
                PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
                PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
              ELSE
                PORBED(L,K)=0.0
                PORBED(L,K)=0.0
              ENDIF
            ENDDO
          ENDDO
!
! ** ADD OR REMOVE LAYERS
!
!          CALL CALBLAY
!
        ENDIF
        ENDIF
!
!----------------------------------------------------------------------C
!
! ++ UPDATE BULK DENSITY
!
        DO K=1,KB
          DO L=2,LA
            SEDBT(L,K)=0.0
            SNDBT(L,K)=0.0
          ENDDO
        ENDDO
!
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(ISTRAN(7).GE.1)THEN
          DO NS=1,NSND
            DO K=1,KB
              DO L=2,LA
                SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
          
        DO K=1,KB
          DO L=2,LA
            IF(HBED(L,K).GT.0.)THEN
              BDENBED(L,K)=1000.*PORBED(L,K)+0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)
            ELSE
              BDENBED(L,K)=0.
            ENDIF
          ENDDO
        ENDDO
!
!----------------------------------------------------------------------C
!
      ENDIF
!
 8011 FORMAT('CBED',I12,6E13.5)    

!**********************************************************************C
!
      RETURN
      END
