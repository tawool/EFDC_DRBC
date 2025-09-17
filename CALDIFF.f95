!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALDIFF (ISTL,M,CON1)
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
! **  SUBROUTINE CALDIFF CALCULATES THE HORIZONTAL DIFFUSIVE 
! **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO
! **  A REVISEDED VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL  
! **  INDICATES THE NUMBER OF TIME LEVELS IN THE STEP
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION CON1(LCM,KCM)
!
!**********************************************************************C
!
! **  HORIZONTAL DIFFUSIVE FLUX CALCULATION
!
!----------------------------------------------------------------------C
!
      DO K=1,KC    
      DO L=2,LA
      LS=LSC(L)      
      FUHU(L,K)=FUHU(L,K)+0.5*SUB(L)*DYU(L)*HU(L)*(AH(L,K)+AH(L-1,K))*(CON1(L-1,K)-CON1(L,K))*     &
                DXIU(L)
      FVHU(L,K)=FVHU(L,K)+0.5*SVB(L)*DXV(L)*HV(L)*(AH(L,K)+AH(LS,K))*(CON1(LS,K)-CON1(L,K))*DYIV(L)
      ENDDO
      ENDDO
!
!**********************************************************************C
!
      RETURN
      END