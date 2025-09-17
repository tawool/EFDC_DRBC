!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALEBIGVC
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
! **  CALEBI CALCULATES THE EXTERNAL BUOYANCY INTEGRALS
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!     DIMENSION CH(LCM,KCM)
!
!**********************************************************************C
!
!      DO L=2,LA
!      BI1(L)=0.
!      BI2(L)=0.
!      CH(L,KC)=DZC(KC)*B(L,KC)
!      BE(L)=GP*DZC(KC)*B(L,KC)
!      ENDDO
!
!      DO K=KS,1,-1
!      DO L=2,LA
!	 IF(K.GE.KGVCP(L))THEN
!        CH(L,K)=CH(L,K+1)+DZC(K)*B(L,K)
!        BE(L)=BE(L)+GP*DZC(K)*B(L,K)
!	 ENDIF
!      ENDDO
!      ENDDO
!
!      DO K=1,KC
!      DO L=2,LA
!	 IF(K.GE.KGVCP(L))THEN
!        BI1(L)=BI1(L)+GP*DZC(K)*(CH(L,K)-0.5*DZC(K)*B(L,K))
!        BI2(L)=BI2(L)+GP*DZC(K)*(CH(L,K)+Z(K-1)*B(L,K))
!	 ENDIF
!      ENDDO
!      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      BI1GVC(L,K)=0.
      BI2GVC(L,K)=0.
      BEGVC(L,K)=0.0
      ENDDO
      ENDDO
!
      DO L=2,LA
      CH(L,KC)=DZC(KC)*B(L,KC)
      BEGVC(L,KC)=GP*DZC(KC)*B(L,KC)
      ENDDO
!
      DO K=KS,1,-1
      DO L=2,LA
	IF(K.GE.KGVCP(L))THEN
        CH(L,K)=CH(L,K+1)+DZC(K)*B(L,K)
        BEGVC(L,K)=BEGVC(L,K+1)+GP*DZC(K)*B(L,K)
	ENDIF
      ENDDO
      ENDDO
!
      DO KBOT=1,KC
      DO K=KBOT,KC
      DO L=2,LA
	IF(K.GE.KGVCP(L))THEN
        BI1GVC(L,KBOT)=BI1GVC(L,KBOT)+GP*DZC(K)*(CH(L,K)-0.5*DZC(K)*B(L,K))
        BI2GVC(L,KBOT)=BI2GVC(L,KBOT)+GP*DZC(K)*(CH(L,K)+Z(K-1)*B(L,K))
	ENDIF
      ENDDO
      ENDDO
      ENDDO
!
!**********************************************************************C
!
      RETURN
      END