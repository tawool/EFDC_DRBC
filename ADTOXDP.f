!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE ADTOXDP(NT,DELTI,PARTDIF)
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
! **  SUBROUTINE ADTOXDP SOLVES TOXIC PORE WATER ADVECTION AND
! **  DIFFUSION IN DOUBLE PRECISION
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
!      REAL*4
!       DELTI                   PASSED AS ARGUMENT    REQUIRES LOCAL DP VERSION
!      BETTMP                  LOCAL
!       DIFBWFAC                LOCAL
!
!      REAL*4
!      DZ(KCM)                 IN GLOBAL COMMON      REQUIRES LOCAL DP VERSION
!
!      REAL*4
!      DIFTOX(NTXM)            IN GLOBAL COMMON      REQUIRES LOCAL DP VERSION
!      DIFTOXS(NTXM)           IN GLOBAL COMMON      REQUIRES LOCAL DP VERSION
!
!      REAL*4
!      HP(LCM)                 IN GLOBAL COMMON      REQUIRES LOCAL DP VERSION
!      TOXBBALO(LCM)           LOCAL
!      DIFTOXBW(LCM)           LOCAL
!      TOXBBALN(LCM)           LOCAL
!      TOXWBALO(LCM)           LOCAL
!      TOXWBALN(LCM)           LOCAL
!
!       REAL*4
!      HBED(LCM,KBM)           IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
!      QWTRBED(LCM,0:KBM)      IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
!       PARTDIF(LCM,KBM)        PASSED AS ARGUMENT   REQUIRES LOCAL DP VERSION
!       PORBED(LCM,KBM)         IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
!
!      REAL*4
!      ALOW(LCM,KBM+1)         LOCAL
!      BMNN(LCM,KBM+1)         LOCAL
!      CUPP(LCM,KBM+1)         LOCAL
!      RRHS(LCM,KBM+1)         LOCAL
!      GAMTMP(LCM,KBM+1)       LOCAL
!      TOXTMP(LCM,KBM+1)       LOCAL
!
!      REAL*4
!      TADFLUX(LCM,NTXM)       IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION  RETURN SP
!      CONGW(LCM,NSTVM)        IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
!                                                FOR CONGW(L,NT+4)
!
!      REAL*4
!      TOX(LCM,KCM,NTXM)       IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION  RETURN SP
!      TOXPFTW(LCM,KCM,NTXM)   IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
!
!      REAL*4
!      TOXB(LCM,KBM,NT)        IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION  RETURN SP
!      TOXPFTB(LCM,KKBM,NTXM)  IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
!
!**********************************************************************C
!
!  **  DECLARE PASSED ARGUMENTS
!
      REAL*4 DELTI
!
      REAL*4 PARTDIF
      DIMENSION PARTDIF(LCM,KBM)
!
!  **  DECLARE LOCAL DOUBLE PRECISION VARIABLES
!
      REAL*8 DELTI_DP,BETTMP,DIFBWFAC
!
      REAL*8 DZC_DP
      DIMENSION DZC_DP(KCM)
!
      REAL*8 DIFTOX_DP,DIFTOXS_DP
      DIMENSION DIFTOX_DP(NTXM),DIFTOXS_DP(NTXM)
!
      REAL*8 DIFTOXBW,TOXBBALO,TOXBBALN,TOXWBALO,TOXWBALN,HP_DP
      DIMENSION DIFTOXBW(LCM),TOXBBALO(LCM),TOXBBALN(LCM),TOXWBALO(LCM),TOXWBALN(LCM),HP_DP(LCM)
!
      REAL*8 PARTDIF_DP,HBED_DP,PORBED_DP
      DIMENSION PARTDIF_DP(LCM,KBM),HBED_DP(LCM,KBM),PORBED_DP(LCM,KBM)
!
      REAL*8 QWTRBED_DP
      DIMENSION QWTRBED_DP(LCM,0:KBM)
!
      REAL*8 ALOW,BMNN,CUPP,RRHS,GAMTMP,TOXTMP
      DIMENSION ALOW(LCM,KBM+1),BMNN(LCM,KBM+1),CUPP(LCM,KBM+1),RRHS(LCM,KBM+1),GAMTMP(LCM,KBM+1),TOXTMP(LCM,KBM+1)
!
      REAL*8 TADFLUX_DP
      DIMENSION TADFLUX_DP(LCM,NTXM)
!
      REAL*8 CONGW_DP
      DIMENSION CONGW_DP(LCM,NSTVM)
!
      REAL*8 TOX_DP,TOXPFTW_DP
      DIMENSION TOX_DP(LCM,KCM,NTXM),TOXPFTW_DP(LCM,KCM,NTXM)
!
      REAL*8 TOXB_DP,TOXPFTB_DP
      DIMENSION TOXB_DP(LCM,KBM,NTXM),TOXPFTB_DP(LCM,KBM,NTXM)
!
!**********************************************************************C
!
!  **  LOAD DOUBLE PRECISION VARIABLS
!
       DELTI_DP=DBLE(DELTI)
       DZC_DP(1)=DBLE(DZC(1))
       DIFTOX_DP(NT)=DBLE(DIFTOX(NT))
       DIFTOXS_DP(NT)=DBLE(DIFTOXS(NT))
!
       DO L=2,LA
         HP_DP(L)=DBLE(HP(L))
         QWTRBED_DP(L,0)=DBLE(QWTRBED(L,0))
         TOX_DP(L,1,NT)=DBLE(TOX(L,1,NT))
         TOXPFTW_DP(L,1,NT)=DBLE(TOXPFTW(L,1,NT))
!         TADFLUX_DP(L,NT)=DBLE(TADFLUX(L,NT))
         CONGW_DP(L,NT+4)=DBLE(CONGW(L,NT+4))
       ENDDO
!
       DO K=1,KB
       DO L=2,LA
         HBED_DP(L,K)=DBLE(HBED(L,K))
         QWTRBED_DP(L,K)=DBLE(QWTRBED(L,K))
         PORBED_DP(L,K)=DBLE(PORBED(L,K))
         PARTDIF_DP(L,K)=DBLE(PARTDIF(L,K))
         TOXB_DP(L,K,NT)=DBLE(TOXB(L,K,NT))
         TOXPFTB_DP(L,K,NT)=DBLE(TOXPFTB(L,K,NT))
       ENDDO
       ENDDO
!
!**********************************************************************C
!
          DO L=2,LA
            DIFTOXBW(L)=0.0
          ENDDO
!
          DO L=2,LA
            IF(LMASKDRY(L)) DIFTOXBW(L)=DIFTOXS_DP(NT)
          ENDDO
!
            DO L=2,LA
!
              DIFBWFAC=2./HBED_DP(L,KBT(L))
              IF(ISDIFBW(NT).EQ.1)DIFBWFAC=1.0
              TOXBBALO(L)=0.
              KBTP1=KBT(L)+1
              KBTM1=KBT(L)-1
              ALOW(L,1)=0.
              CUPP(L,KBTP1)=0.
!
              DO K=1,KBTM1
!###############################################################################
! HQI Change to implement particle mixing as a rate in units of length/time
! instead of a mixing coefficient in units of length**2/time
! RM, 09/28/05
!                CUPP(L,K)=MIN(QWTRBED_DP(L,K),0.)
!     &          -(DIFTOX_DP(NT)+PARTDIF_DP(L,K))*(PORBED_DP(L,K)
!     &                  +PORBED_DP(L,K+1))/
!     &                  (HBED_DP(L,K)+HBED_DP(L,K+1))
	CUPP(L,K)=DMIN1(QWTRBED_DP(L,K),0.)-((DIFTOX_DP(NT)/(HBED_DP(L,K)+HBED_DP(L,K+1)))+PARTDIF_DP(L,K))*(PORBED_DP(L,K)+PORBED_DP(L,K+1))
!###############################################################################
              ENDDO
  CUPP(L,KBT(L))=AMIN1(QWTRBED_DP(L,KBT(L)),0.)-DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
!
              DO K=2,KBT(L)
!###############################################################################
! HQI Change to implement particle mixing as a rate in units of length/time
! instead of a mixing coefficient in units of length**2/time
! RM, 09/28/05
!                ALOW(L,K)=-MAX(QWTRBED_DP(L,K-1),0.)
!     &          -(DIFTOX_DP(NT)+PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)
!     &                  +PORBED_DP(L,K))
!     &                  /(HBED_DP(L,K-1)+HBED_DP(L,K))
	ALOW(L,K)=-AMAX1(QWTRBED_DP(L,K-1),0.)-((DIFTOX_DP(NT)/(HBED_DP(L,K-1)+HBED_DP(L,K)))+ARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)+PORBED_DP(L,K))
!###############################################################################
              ENDDO
	ALOW(L,KBTP1)=-AMAX1(QWTRBED_DP(L,KBT(L)),0.)-DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
!
              DO K=1,KBT(L)
                BMNN(L,K)=DELTI_DP*HBED_DP(L,K)*PORBED_DP(L,K)/(1.-TOXPFTB_DP(L,K,NT))
              ENDDO
              BMNN(L,KBTP1)=DELTI_DP*DZC_DP(1)*HP_DP(L)/(1.-TOXPFTW_DP(L,1,NT))
!
!###############################################################################
! HQI Change to implement particle mixing as a rate in units of length/time
! instead of a mixing coefficient in units of length**2/time
! RM, 09/28/05
!              BMNN(L,1)=BMNN(L,1)+MAX(QWTRBED_DP(L,1),0.)
!     &        +(DIFTOX_DP(NT)+PARTDIF_DP(L,1))*(PORBED_DP(L,2)
!     &                  +PORBED_DP(L,1))/
!     &                  (HBED_DP(L,2)+HBED_DP(L,1))
	BMNN(L,1)=BMNN(L,1)+AMAX1(QWTRBED_DP(L,1),0.) +((DIFTOX_DP(NT)/(HBED_DP(L,2)+HBED_DP(L,1)))+ PARTDIF_DP(L,1))*(PORBED_DP(L,2)+PORBED_DP(L,1))
!###############################################################################
              DO K=2,KBTM1
!###############################################################################
! HQI Change to implement particle mixing as a rate in units of length/time
! instead of a mixing coefficient in units of length**2/time
! RM, 09/28/05
!                BMNN(L,K)=BMNN(L,K)+MAX(QWTRBED_DP(L,K),0.)
!     &         +(DIFTOX_DP(NT)+PARTDIF_DP(L,K))*(PORBED_DP(L,K+1)
!     &                  +PORBED_DP(L,K))/
!     &                  (HBED_DP(L,K+1)+HBED_DP(L,K))
!     &                         -MIN(QWTRBED_DP(L,K-1),0.)
!     &         +(DIFTOX_DP(NT)+PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)
!     &                  +PORBED_DP(L,K))/
!     &                  (HBED_DP(L,K-1)+HBED_DP(L,K))
	BMNN(L,K)=BMNN(L,K)+AMAX1(QWTRBED_DP(L,K),0.) +((DIFTOX_DP(NT)/(HBED_DP(L,K+1)+HBED_DP(L,K)))+ PARTDIF_DP(L,K))*(PORBED_DP(L,K+1)+PORBED_DP(L,K))-AMIN1(QWTRBED_DP(L,K-1),0.)+((DIFTOX_DP(NT)/(HBED_DP(L,K-1)+HBED_DP(L,K)))+PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)+PORBED_DP(L,K))
!###############################################################################
              ENDDO
              K=KBT(L)
!###############################################################################
! HQI Change to implement particle mixing as a rate in units of length/time
! instead of a mixing coefficient in units of length**2/time
! RM, 09/28/05
!              BMNN(L,K)=BMNN(L,K)+MAX(QWTRBED_DP(L,K),0.)
!     &           +DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
!     &                         -MIN(QWTRBED_DP(L,K-1),0.)
!     &           +(DIFTOX_DP(NT)+PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)
!     &                  +PORBED_DP(L,K))/(HBED_DP(L,K-1)+HBED_DP(L,K))
	MNN(L,K)=BMNN(L,K)+AMAX1(QWTRBED_DP(L,K),0.) +DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L)) -AMIN1(QWTRBED_DP(L,K-1),0.) +((DIFTOX_DP(NT)/(HBED_DP(L,K-1)+HBED_DP(L,K))) +PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)+PORBED_DP(L,K))
!###############################################################################
              K=KBTP1
	MNN(L,K)=BMNN(L,K)-AMIN1(QWTRBED_DP(L,K-1),0.) +DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
!
              DO K=1,KBT(L)
                RRHS(L,K)=DELTI_DP*TOXB_DP(L,K,NT)
                TOXBBALO(L)=TOXBBALO(L)+TOXB_DP(L,K,NT)
              ENDDO
	RHS(L,1)=RRHS(L,1)+AMAX1(QWTRBED_DP(L,0),0.) *CONGW_DP(L,NT+4)
              RRHS(L,KBTP1)=DELTI_DP*DZC_DP(1)*HP_DP(L)*TOX_DP(L,1,NT)
              TOXWBALO(L)=DZC_DP(1)*HP_DP(L)*TOX_DP(L,1,NT)
!
            ENDDO
!
! **  TRI-DIAGONAL SOLVER
!
            DO L=2,LA
              KBTP1=KBT(L)+1
              BETTMP=BMNN(L,1)
              TOXTMP(L,1)=RRHS(L,1)/BETTMP
              DO KK=2,KBTP1
                GAMTMP(L,KK)=CUPP(L,KK-1)/BETTMP
                BETTMP=BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
                TOXTMP(L,KK)=(RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))/ BETTMP
              ENDDO
              DO KK=KBT(L),1,-1
                TOXTMP(L,KK)=TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)
              ENDDO
            ENDDO
!
! **  CONVERT SCALED SOLUTION VARIABLES AND CALCULATE FINAL MASS
!
            DO L=2,LA
              TOXBBALN(L)=0.0
              KBTP1=KBT(L)+1
              DO K=1,KBT(L)
!                TOXBSMB(L,K)=TOXB_DP(L,K,NT)
                TOXB_DP(L,K,NT)=HBED_DP(L,K)*PORBED_DP(L,K)*TOXTMP(L,K)/(1.-TOXPFTB_DP(L,K,NT))
                TOXBBALN(L)=TOXBBALN(L)+TOXB_DP(L,K,NT)
              ENDDO
!              TOXSMB(L)=TOX_DP(L,1,NT)
              TOX_DP(L,1,NT)=TOXTMP(L,KBTP1)/(1.-TOXPFTW_DP(L,1,NT))
              TOXWBALN(L)=DZC_DP(1)*HP_DP(L)*TOX_DP(L,1,NT)
            ENDDO
!
! **  CALCULATED ADVECTION-DIFFUSION FLUX FOR MASS BALANCE
!
            DO L=2,LA
              TADFLUX_DP(L,NT)=DELTI_DP*(DZC_DP(1)*HP_DP(L) *TOX_DP(L,1,NT)-TOXWBALO(L))
            ENDDO
!
!**********************************************************************C
!
! **  RETURN UPDATED VARIABLES IN SINGLE PRECISION BACK TO GLOBAL
! **  COMMON
!
!
       DO L=2,LA
         TOX(L,1,NT)=SNGL(TOX_DP(L,1,NT))
         TADFLUX(L,NT)=SNGL(TADFLUX_DP(L,NT))
       ENDDO
!
       DO K=1,KB
       DO L=2,LA
         TOXB(L,K,NT)=SNGL(TOXB_DP(L,K,NT))
       ENDDO
       ENDDO
!
!**********************************************************************C
!
      RETURN
      END
