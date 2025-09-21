!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CBALEV4
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
! **  SUBROUTINES CBALEV CALCULATE GLOBAL VOLUME, MASS, MOMENTUM,
! **  AND ENERGY BALANCES
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
! **  CALCULATE MOMENTUM AND ENERGY DISSIPATION
!
!----------------------------------------------------------------------C
!
      DO L=2,LA
      LN=LNC(L)
!     UUEOUTE=UUEOUTE+SPB(L)*SPB(L+1)*DXYU(L)
!    &      *(U(L,1)*TBX(L)-U(L,KC)*TSX(L))
!     VVEOUTE=VVEOUTE+SPB(L)*SPB(LN)*DXYV(L)
!    &      *(V(L,1)*TBY(L)-V(L,KC)*TSY(L))
      UUEOUTE=UUEOUTE+0.5*SPB(L)*DXYP(L)*(U(L,1)*TBX(L)+U(L+1,1)*TBX(L+1)-U(L,KC)*TSX(L)-U(L+1,KC) &
              *TSX(L+1))
      VVEOUTE=VVEOUTE+0.5*SPB(L)*DXYP(L)*(V(L,1)*TBY(L)+V(LN,1)*TBX(LN)-V(L,KC)*TSY(L)-V(LN,KC)    &
              *TSX(LN))
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      DUTMP=0.5*( U(L,K+1)+U(L+1,K+1)-U(L,K)-U(L+1,K) )
      DVTMP=0.5*( V(L,K+1)+V(LN,K+1)-V(L,K)-V(LN,K) )
      UUEOUTE=UUEOUTE+SPB(L)*2.0*DXYP(L)*AV(L,K)*( DUTMP*DUTMP )/(DZC(K+1)+DZC(K))
      VVEOUTE=VVEOUTE+SPB(L)*2.0*DXYP(L)*AV(L,K)*( DVTMP*DVTMP )/(DZC(K+1)+DZC(K))
      BBEOUTE=BBEOUTE+SCB(L)*DXYP(L)*HP(L)*GP*AB(L,K)*(B(L,K+1)-B(L,K))
      ENDDO
      ENDDO
!
!**********************************************************************C
!
      RETURN
      END