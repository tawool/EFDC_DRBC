!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALBUOY
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
! **  CALBUOY CALCULATES THE BUOYANCY USING MELLOR'S APPROXIMATION
! **  TO THE UNESCO EQUATION OF STATE (MELLOR, G.L., J. ATM AND OCEAN
! **  TECH, VOL 8, P 609)
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      IF(IBSC.EQ.1) GOTO 1000
!
      ISPCOR=0
!
! **  DENSITY RHOO AT P=0, S=0, AND T=TEMO
!
      RHOO=999.842594+6.793952E-2*TEMO-9.095290E-3*TEMO*TEMO+1.001685E-4*TEMO*TEMO*               &
           TEMO-1.120083E-6*TEMO*TEMO*TEMO*TEMO+6.536332E-9*TEMO*TEMO*TEMO*TEMO*TEMO
!
! **  DENSITY RHO=B(L,K) AT P=0, S=SAL(L,K), AND T=TEM(L,K)
!
      IF(ISTRAN(1).EQ.0.AND.ISTRAN(2).EQ.0)THEN
!
      DO K=1,KC
      DO L=2,LA
      B(L,K)=RHOO
      ENDDO
      ENDDO
!
      ENDIF
C
      IF(ISTRAN(1).GE.1.AND.ISTRAN(2).EQ.0)THEN
!
      DO K=1,KC
      DO L=2,LA
      SAL(L,K)=MAX(SAL(L,K),0.)
      SSTMP=SAL(L,K)
      TTMP=TEMO
      B(L,K)=RHOO+SSTMP*(0.824493-4.0899E-3*TTMP+7.6438E-5*TTMP*TTMP-8.2467E-7*TTMP*TTMP*TTMP+    &
             5.3875E-9*TTMP*TTMP*TTMP*TTMP)+SQRT(SSTMP)*SSTMP*(-5.72466E-3+1.0227E-4*TTMP         &
             -1.6546E-6*TTMP*TTMP)+4.8314E-4*SSTMP*SSTMP
      ENDDO
      ENDDO
!
      ENDIF
!
      IF(ISTRAN(1).EQ.0.AND.ISTRAN(2).GE.1)THEN
!
      DO K=1,KC
      DO L=2,LA
      TTMP=TEM(L,K)
      B(L,K)=999.842594+6.793952E-2*TTMP-9.095290E-3*TTMP*TTMP+1.001685E-4*TTMP*TTMP*TTMP         &
             -1.120083E-6*TTMP*TTMP*TTMP*TTMP+6.536332E-9*TTMP*TTMP*TTMP*TTMP*TTMP
      ENDDO
      ENDDO
!
      ENDIF
!
      IF(ISTRAN(1).GE.1.AND.ISTRAN(2).GE.1)THEN
!
      DO K=1,KC
      DO L=2,LA
      SAL(L,K)=MAX(SAL(L,K),0.)
      SSTMP=SAL(L,K)
      TTMP=TEM(L,K)
      RHTMP=999.842594+6.793952E-2*TTMP-9.095290E-3*TTMP*TTMP+1.001685E-4*TTMP*TTMP*TTMP          &
            -1.120083E-6*TTMP*TTMP*TTMP*TTMP+6.536332E-9*TTMP*TTMP*TTMP*TTMP*TTMP
      B(L,K)=RHTMP+SSTMP*(0.824493-4.0899E-3*TTMP+7.6438E-5*TTMP*TTMP-8.2467E-7*TTMP*TTMP*TTMP    &
             +5.3875E-9*TTMP*TTMP*TTMP*TTMP)+SQRT(SSTMP)*SSTMP*(-5.72466E-3+1.0227E-4*TTMP        &
             -1.6546E-6*TTMP*TTMP)+4.8314E-4*SSTMP*SSTMP
      ENDDO
      ENDDO
!
      ENDIF
!
! **  APPLY MELLOR'S PRESSURE CORRECTION
!
      IF(ISPCOR.EQ.1)THEN
!
      DO K=1,KC
      DO L=2,LA
      PRES=RHOO*G*HP(L)*(1.-ZZ(K))*1.E-6
      CCON=1449.2+1.34*(SAL(L,K)-35.)+4.55*TEM(L,K)-0.045*TEM(L,K)*TEM(L,K)+0.00821*PRES+15.E-9   &
           *PRES*PRES
      TMP=PRES/(CCON*CCON)
      B(L,K)=B(L,K)+1.E+4*TMP*(1.-0.2*TMP)
      ENDDO
      ENDDO
!
      ENDIF 
!
! **  REPLACE DENSITY B(L,K) WITH BUOYANCY B(L,K)
!
      DO K=1,KC
      DO L=2,LA
      B(L,K)=(B(L,K)/RHOO)-1.
      ENDDO
      ENDDO
!
! **  APPLY LOW SEDIMENT CONCENTRATION CORRECTION TO BUOYANCY
!
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
!
      DO K=1,KC
       DO L=2,LA
        TVAR1S(L,K)=0.
        TVAR1W(L,K)=0.
       ENDDO
      ENDDO
!
      ENDIF
!
      IF(ISTRAN(6).GE.1)THEN
!
      DO NS=1,NSED
      DO K=1,KC
      DO L=2,LA
        TVAR1S(L,K)=TVAR1S(L,K)+SDEN(NS)*SED(L,K,NS)
        TVAR1W(L,K)=TVAR1W(L,K)+(SSG(NS)-1.)*SDEN(NS)*SED(L,K,NS)
      ENDDO
      ENDDO
      ENDDO
!
      ENDIF
!
!      IF(ISTRAN(6).GE.1)THEN
!
!      DO NS=1,NSED
!      DO K=1,KC
!      DO L=2,LA
!        B(L,K)=B(L,K)*(1.-SDEN(NS)*SED(L,K,NS))
!     &        +(SSG(NS)-1.)*SDEN(NS)*SED(L,K,NS)
!      ENDDO
!      ENDDO
!      ENDDO
!
!      ENDIF
!
      IF(ISTRAN(7).GE.1)THEN
!
      DO NN=1,NSND
      NS=NN+NSED
      DO K=1,KC
      DO L=2,LA
        TVAR1S(L,K)=TVAR1S(L,K)+SDEN(NS)*SND(L,K,NN)
        TVAR1W(L,K)=TVAR1W(L,K)+(SSG(NS)-1.)*SDEN(NS)*SND(L,K,NN)
      ENDDO
      ENDDO
      ENDDO
!
      ENDIF
!
!      IF(ISTRAN(7).GE.1)THEN
!
!      DO NN=1,NSND
!      NS=NN+NSED
!      DO K=1,KC
!      DO L=2,LA
!        B(L,K)=B(L,K)*(1.-SDEN(NS)*SND(L,K,NN))
!     &        +(SSG(NS)-1.)*SDEN(NS)*SND(L,K,NN)
!      ENDDO
!      ENDDO
!      ENDDO
!
!      ENDIF
!
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
!
      DO K=1,KC
       DO L=2,LA
        B(L,K)=B(L,K)*(1.-TVAR1S(L,K))+TVAR1W(L,K)
       ENDDO
      ENDDO
!
      ENDIF
!
      GOTO 2000
!
!**********************************************************************C
!
!     DENSITY AS A LINEAR FUNCTION OF SALINITY ONLY.  FOR DIAGNOSTIC
!     PURPOSES ONLY
!
 1000 CONTINUE
!
      DO K=1,KC
      DO L=2,LA
      B(L,K)=0.00075*SAL(L,K)
      ENDDO
      ENDDO
!      
!**********************************************************************C
!
 2000 CONTINUE
!
      RETURN
      END