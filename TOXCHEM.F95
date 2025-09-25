!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE TOXCHEM
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
                      USTAR(LCM),UCELLCTR(LCM),VCELLCTR(LCM),                                     &
                      QWBDTOP(LCM),QSBDTOP(LCM),ZETATOP(LCM),                                     &
                      ZBEDGT(LCM),QSBDLDP(LCM),RBPSBL(LCM),                                       &
                      TOXBBALO(LCM),TOXBBALN(LCM),TOXWBALO(LCM),                                  &
                      ROUSE(LCM),ZEQ(LCM),ZEQI(LCM),ZEQD(LCM),                                    &
                      ZEQDI(LCM),SNDEQ(LCM),SNDEQB(LCM),                                          &
                      TOXWBALN(LCM),FACSUSL(LCM),FACBEDL(LCM),                                    &
                      PEXP(LCM,NSNM),PHID(LCM,NSNM),                                              &
                      SEDFPA(LCM,NSCM),SNDFPA(LCM,NSNM)
!
      COMMON/SSEDTOX1A/ CBEDTOTAL(LCM),QCELLCTR(LCM),HGDH(LCM),FRACCOH(LCM,KBM),                  &
                        FRACNON(LCM,KBM),USTARSED(LCM),USTARSND(LCM),QWATPA(LCM),QSSDPA(LCM)
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
      COMMON/SSEDTOX6/ DELT,DELTI,DSEDGMM,FOURDPI,SEDMDGM,S2TL,S3TL,                              &
                      CORDT,DIASED,GPDIASED,BEDEX
!
      COMMON/SSEDTOX7/ ISTL,IS2TL,ISUD
!
! ** BEGIN ADDITIONS REQUIRED BY HOUSATONIC OPTION
!
      COMMON/TOXCHEM_TEMP/ W_TEMP, SWIND, DIFFA, DIFFW,A_TEMP, HENRY(NTXM), ATMOS, VOLA_TXW(LCM)
! ===>  From Wasp Volatilization Code (1 February 1990)
      REAL KLT, KLA, KAW, KL_RIV, KL_LAKE_O, KL_LAKE_M, KG
      REAL KL, VEL(LCM)
! ===>  End Addition   JAH
!
! ** END ADDITIONS REQUIRED BY HOUSATONIC OPTION
!
!**********************************************************************C
!
      IF(ISHOUSATONIC.EQ.0)THEN
      IF(ISTRAN(5).GE.1)THEN
        DO NT=1,NTOX
!
!----------------------------------------------------------------------C
!
! **     NOTES:
!
!        BULK DECAY COEFFICIENT
!
!          RKTOXWT=RKTOXW(NT)    !*(TEM(L,K)-TKTOXW(NT))**ALPTOX
!          RKTOXBT=RKTOXB(NT)    !*(TEM(L,1)-TRTOXB(NT))**ALPTOX
!
!        VOLITIZATION
!
!          VOLTOX(NT)
!          RMOLTX(NT)=MOLECULAR WEIGHT
!
!        PHOTOLOSIS
!
!          RKTOXP(NT)=BASE RATE
!          SKTOXP(NT)=SOLAR RADIATION AT BASE RATE
!
          DO K=1,KC
            DO L=2,LA
              CDECAYW(L,K)=1./(1.+DELT*RKTOXW(NT))
            ENDDO
          ENDDO
!
          DO K=1,KC
            DO L=2,LA
              TOX(L,K,NT)=CDECAYW(L,K)*TOX(L,K,NT)
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              CDECAYB(L,K)=1./(1.+DELT*RKTOXB(NT))
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              TOXB(L,K,NT)=CDECAYB(L,K)*TOXB(L,K,NT)
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------C
!
        ENDDO
      ENDIF
	ENDIF
!
!**********************************************************************C
!
      IF(ISHOUSATONIC.EQ.1)THEN
      IF(ISTRAN(5).GE.1)THEN
        DO NT=1,NTOX
!
!----------------------------------------------------------------------C
!
! **     NOTES:
!
!        BULK DECAY COEFFICIENT
!
!          RKTOXWT=RKTOXW(NT)    !*(TEM(L,K)-TKTOXW(NT))**ALPTOX
!          RKTOXBT=RKTOXB(NT)    !*(TEM(L,1)-TRTOXB(NT))**ALPTOX
!
!
!     VOLATILIZATION FROM WASP 5 VALATILIZATION CODE (1 FEBRUARY 1990,
!     VOLAT.FOR)
!     SET VISCOUS SUBLAYER AND DRAG COEFFICIENTS FOR VOLATILIZATION
          XLAM2 = 4.
          CDRAG = 0.0011
          HENRY(NT)=3.85E-4
          RMOLTX(NT)=320.39
          ATMOS=0.
          W_TEMP=10.0 ! LT mean temp.
          A_TEMP=10.0 ! LT mean temp.
          SWIND=3.5756
          KLT=1.024
          DIFFA=1.9E-4*(RMOLTX(NT)**(-2./3.))
          DIFFW=2.2E-8*(RMOLTX(NT)**(-2./3.))
          DO L=2,LA
             UTMP1=50.*(UHDYE(L+1)+UHDYE(L))/(DYP(L)*HP(L))
             VTMP1=50.*(VHDXE(L)+VHDXE(L))/(DXP(L)*HP(L))
             IF(SPB(L).EQ.0)THEN
                UTMP1=2.*UTMP1
                VTMP1=2.*VTMP1
             ENDIF
             UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
             VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
             VEL(L)=SQRT(UTMP*UTMP+VTMP*VTMP)/100.
          ENDDO
!      VOLUME OPTIONS COMPUTED FOR WASP OPTIONS 3, 4, AND 5
!      COMPUTE VIA ALL METHODS FOR STREAMS/RIVERS AND LAKES/RESERVOIRS.  PICK
!      MAXIMUM FOR KL
!
!      STREAM, RIVER OR ESTUARY
          KG=100.  ! (Gas transfer Coefficient for flowing systems, M/D
          DO L=2,LA
            STP20 = W_TEMP - 20.
            KAW = HENRY(NT)/(8.206E-05*(W_TEMP + 273.15))
!      CHOOSE KL EQUATION BASED ON DEPTH & VELOCITY:
            XX = SQRT (32./RMOLTX(NT))
!      OWENS:
            IF (HP(L) .LT. 0.61) THEN
              XKL = XX*5.349*(VEL(L)**.67)/(HP(L)**.85)
            ELSE
!      O'CONNOR-DOBBINS:
              IF (VEL(L) .LT. 0.518 .OR. HP(L) .GT. 13.584*VEL(L)**2.9135) THEN
                XKL = SQRT (DIFFW*VEL(L)/HP(L))*86400.
              ELSE
!      CHURCHILL:
                XKL = XX*5.049*VEL(L)**.969/HP(L)**.673
              ENDIF
            ENDIF
            XKGH = KG * KAW
!jah
            KL_RIV = (1./(1./(XKL+1.e-12) + 1./XKGH))*KLT**STP20
!jah        KL_RIV = (1./(1./XKL + 1./XKGH))*KLT**STP20
!            if((L.EQ.103).AND.(TIMESEC-2.8918080E+08.GT.0)) then
!               write(*,*)VEL(L), HP(L), HENRY(NT), W_TEMP
!               write(*,*)RMOLTX(NT), DIFFW, DIFFA, KLT, STP20
!               write(*,*)"Stream equation"
!               write(*,*)KL_RIV, XKL, XKGH
!            endif
!
!              RESERVOIR, LAKE OR POND
!
!     COMPUTE DENSITY (G/ML) AND VISCOSITY (M**2/S) OF AIR AND WATER
!
            DENA = 0.001293/(1. + 0.00367*A_TEMP)
            DENW = 1. - 8.8E-05*W_TEMP
            XNUA = (1.32 + 0.009*A_TEMP)*10.0
            XNUW = (10.**(1301./(998.333 + 8.1855*STP20 + 0.00585*STP20**2) -3.30233)/DENW)*1.E-04
!
!     COMPUTE SCHMIDT NUMBERS FOR AIR AND WATER
!
            SCA = XNUA/DIFFA
            SCW = XNUW/DIFFW
!
!     COMPUTE AIR LAYER AND WATER LAYER TRANSFER RATES
!          USE O'CONNOR OR MACKAY AS PER USER OPTION
!
!
!     O'CONNOR:
!        NOTE: 0.905 = (VON KARMON CONSTANT)**0.333
!
            USTAR_VOL = SQRT (CDRAG)*SWIND
            XKL = USTAR_VOL*SQRT (DENA/DENW)*(.905/XLAM2)*(1./SCW)**.666 + 1.0E-9
            XKG = USTAR_VOL*(.905/XLAM2)*((1./SCA)**.666) + 1.0E-9
            XKGH = XKG*KAW
            KL_LAKE_O = (1./(1./XKL + 1./XKGH))*KLT**STP20*86400.

!            if((L.EQ.103).AND.(TIMESEC-2.8918080E+08.GT.0)) then
!               write(*,*)"OConnor Lake equation"
!               write(*,*)KL_LAKE_O, XKL, XKG, USTAR_VOL,XNUA,XNUW
!               write(*,*)SCA,SCW,USTAR_VOL
!            endif

!     MACKAY:
!
            USTAR_VOL = .01*SWIND*SQRT (6.1 + 0.63*SWIND)
            IF (USTAR_VOL .GT. 0.3) XKL = USTAR_VOL*0.00341*(1./SCW)**.5 + 1.E-06
            IF (USTAR_VOL .LE. 0.3) XKL = USTAR_VOL**2.2*0.0144*(1./SCW)**.5 + 1.E-06
            XKG = USTAR_VOL*.0462*(1./SCA)**.666 + 1.E-03
            XKGH = XKG*KAW
            KL_LAKE_M = 0.
!            KL_LAKE_M = (1./(1./XKL + 1./XKGH))*KLT**STP20*86400.

!            if((L.EQ.103).AND.(TIMESEC-2.8918080E+08.GT.0)) then
!               write(*,*)"Mackay Lake equation"
!               write(*,*)KL_LAKE_M, XKL, XKGH, XKG, USTAR_VOL
!               stop
!            endif
!
!     Take max volatilization rate of all 3 equations
            KL=MAX(KL_LAKE_O,KL_LAKE_M,KL_RIV)
            KLA = KL/HP(L)/86400.     ! KLA IN SECONDS**-1
!     Volatilization flux in ug/l/s
            VOLA_TXW(L) = KLA*((TOXFDFW(L,KC,NT)*TOX(L,KC,NT)) - ATMOS/KAW)
!            CVOLA_TXW = 1./(1.+DELT*KLA)
!           UPDATE SURFACE CELL CONCENTRATIONS ONLY
!            TOX(L,KC,NT)=CVOLA_TXW*TOX(L,K,NT)
            TOX(L,KC,NT)=TOX(L,KC,NT)-(VOLA_TXW(L)*DELT)
!
            IF(ISTOXALL.EQ.1.AND.NT.EQ.1)THEN
            ATOXVOL(L,1)=ATOXVOL(L,1)-DXYP(L)*HP(L)*(VOLA_TXW(L)*DELT)
            ENDIF
!
          ENDDO
!
!        PHOTOLOSIS
!
!          RKTOXP(NT)=BASE RATE
!          SKTOXP(NT)=SOLAR RADIATION AT BASE RATE
!
          DO K=1,KC
            DO L=2,LA
              CDECAYW(L,K)=1./(1.+DELT*RKTOXW(NT))
            ENDDO
          ENDDO
!
          DO K=1,KC
            DO L=2,LA
              TOX(L,K,NT)=CDECAYW(L,K)*TOX(L,K,NT)
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              CDECAYB(L,K)=1./(1.+DELT*RKTOXB(NT))
            ENDDO
          ENDDO
!
          DO K=1,KB
            DO L=2,LA
              TOXB(L,K,NT)=CDECAYB(L,K)*TOXB(L,K,NT)
            ENDDO
          ENDDO
!
!----------------------------------------------------------------------C
!
        ENDDO
      ENDIF
      ENDIF
!
!**********************************************************************C
!
      RETURN
      END
