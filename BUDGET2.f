!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE BUDGET2
!
! **  ADDED BY DON KINGERY, CH2M-HILL ON 15 OCTOBER 1996
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
! **  SUBROUTINES BUDGETN CALCULATE SEDIMENT BUDGET (TOTAL SEDIMENTS)
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)

      VOLMOUT=VOLMOUT-VHDX2(LN,K)*DZC(K)
      SMASSOUT=SMASSOUT-MIN(VHDX2(LN,K),0.)*SAL1(LN,K)*DZC(K)-MAX(VHDX2(LN,K),0.)*SAL1(L,K)*DZC(K)

!
!
!  ADDED SEDIMENT FLUXES         DLK 9/27
!
       DO NT=1,NSED
          SEDOUT=SEDOUT-MIN(VHDX2(LN,K),0.)*SED1(LN,K,NT)*DZC(K)-MAX(VHDX2(LN,K),0.)*SED1(L,K,NT)*DZC(K)
       ENDDO
       DO NT=1,NSND
        SEDOUT=SEDOUT-MIN(VHDX2(LN,K),0.)*SND1(LN,K,NT)*DZC(K)-MAX(VHDX2(LN,K),0.)*SND1(L,K,NT)*DZC(K)
      ENDDO
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)

      VOLMOUT=VOLMOUT-UHDY2(L+1,K)*DZC(K)
      SMASSOUT=SMASSOUT-MIN(UHDY2(L+1,K),0.)*SAL1(L+1,K)*DZC(K)-MAX(UHDY2(L+1,K),0.)*SAL1(L,K)*DZC(K)

!
!  ADDED SEDIMENT FLUXES         DLK 10/15
!
      DO NT=1,NSED
        SEDOUT=SEDOUT-MIN(UHDY2(L+1,K),0.)*SED1(L+1,K,NT)*DZC(K)-MAX(UHDY2(L+1,K),0.)*SED1(L,K,NT)*DZC(K)
      ENDDO
      DO NT=1,NSND
        SEDOUT=SEDOUT-MIN(UHDY2(L+1,K),0.)*SND1(L+1,K,NT)*DZC(K)-MAX(UHDY2(L+1,K),0.)*SND1(L,K,NT)*DZC(K)
      ENDDO
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)

      VOLMOUT=VOLMOUT+UHDY2(L,K)*DZC(K)
      SMASSOUT=SMASSOUT+MIN(UHDY2(L,K),0.)*SAL1(L,K)*DZC(K)+MAX(UHDY2(L,K),0.)*SAL1(L-1,K)*DZC(K)
!
!  ADDED SEDIMENT FLUXES         DLK 10/15
!
      DO NT=1,NSED
        SEDOUT=SEDOUT+MIN(UHDY2(L,K),0.)*SED1(L,K,NT)*DZC(K)+MAX(UHDY2(L,K),0.)*SED1(L-1,K,NT)*DZC(K)
      ENDDO
      DO NT=1,NSND
        SEDOUT=SEDOUT+MIN(UHDY2(L,K),0.)*SND1(L,K,NT)*DZC(K)+MAX(UHDY2(L,K),0.)*SND1(L-1,K,NT)*DZC(K)
      ENDDO
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)

      VOLMOUT=VOLMOUT+VHDX2(L,K)*DZC(K)
      SMASSOUT=SMASSOUT+MIN(VHDX2(L,K),0.)*SAL1(L,K)*DZC(K)+MAX(VHDX2(L,K),0.)*SAL1(LS,K)*DZC(K)
!
!  ADDED SEDIMENT FLUXES         DLK 10/15
!
      DO NT=1,NSED
        SEDOUT=SEDOUT+MIN(VHDX2(L,K),0.)*SED1(L,K,NT)*DZC(K)+MAX(VHDX2(L,K),0.)*SED1(LS,K,NT)*DZC(K)
      ENDDO
      DO NT=1,NSND
        SEDOUT=SEDOUT+MIN(VHDX2(L,K),0.)*SND1(L,K,NT)*DZC(K)+MAX(VHDX2(L,K),0.)*SND1(LS,K,NT)*DZC(K)
      ENDDO
!
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  ACCUMULATE FLUX OF SED POSITIVE FROM BED TO WATER COLUMN
!
!       DO NS=1,NSED
!        DO K=1,KB 
!        DO L=2,LA 
!         VOLBW3(L,K)=VOLBW3(L,K)+SCB(L)*DXYP(L)*SEDF(L,0,NS)
!        ENDDO
!       ENDDO 
!       ENDDO 
!
!       DO NS=1,NSND
!        DO K=1,KB 
!        DO L=2,LA 
!         VOLBW3(L,K)=VOLBW3(L,K)+SCB(L)*DXYP(L)*SNDF(L,0,NS)
!        ENDDO
!       ENDDO
!       ENDDO 
!
!**********************************************************************C
!
      RETURN
      END