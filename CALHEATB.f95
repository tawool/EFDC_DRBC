!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALHEATB(ISTL)
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
! **  SUBROUTINE CALHEATB CALCULATES TEMPERATURE DISTRIBUTION IN
! **  BED FOR THERMAL SIMULATION
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!
      DIMENSION BETABED(LCM),BETBED(LCM),DIFFTIMI(LCM)
      DIMENSION ABEDTEM(LCM,KBHM),BBEDTEM(LCM,KBHM),CBEDTEM(LCM,KBHM),RBEDTEM(LCM,KBHM),          &
                UBEDTEM(LCM,KBHM),GBEDTEM(LCM,KBHM),DISTINT(LCM,KBHM-1)
!
!**********************************************************************C
!
!      DELT=DT2
      DELT=FLOAT(NTSTBC)*DT
      S3TL=1.0
      S2TL=0.0
	NVAL=NTSTBC-1
      IF(IS2TIM.EQ.1)THEN
	  NVAL=1
        IF(ISDYNSTP.EQ.0)THEN
          DELT=DT
        ELSE
          DELT=DTDYN
        END IF
        S3TL=0.0
        S2TL=1.0
      ENDIF
!
      BELTMPMAX=-1.E6
      BELTMPMIN=1.E6
	DO L=1,LA
	  BELTMPMAX=MIN(BELTMPMAX,BELV(L))
	  BELTMPMIN=MIN(BELTMPMIN,BELV(L))
      ENDDO
!
      RBADJ=1.0
	BELVDIF=ABS(BELTMPMAX-BELTMPMIN)
	IF(BELVDIF.GT.DABEDT) RBADJ=0.0
!
!**********************************************************************C
!
! **  INITIAILIZE TEMPERATURE DISTRIBUTION IN BED
!
      IF(ISBEDTEMI.EQ.1.AND.N.EQ.NVAL)THEN
!
      DO L=1,LA
        BETABED(L)=(DABEDT+RBADJ*(BELV(L)-BELTMPMIN))/SQRT(10.03823E+6*HTBED2)
      ENDDO
!
      DZBED=1./FLOAT(KBH-1)
	ZVAL=-0.5*DZBED
	DO K=1,KBH-1
	  ZVAL=ZVAL+DZBED
	  DO L=2,LA
	    DISTINT(L,K)=EXP(BETABED(L)*(ZVAL-1.))
	  ENDDO
	ENDDO
!
      DO K=1,KBH-1
	  DO L=2,LA
	    TEMB(L,K)=TBEDIT(L)+(TEM(L,1)-TBEDIT(L))*DISTINT(L,K)
        ENDDO
	ENDDO
!
      OPEN(1,FILE='TEMBINIT.OUT')
	WRITE(1,101)N,NVAL,DELT
	DO L=2,LA
	  TMPTHICK=DABEDT+RBADJ*(BELV(L)-BELTMPMIN)
	  WRITE(1,111)L,IL(L),JL(L),TEM(L,1),TBEDIT(L),(TEMB(L,K),K=1,KBH-1),TMPTHICK
	ENDDO
	CLOSE(1)
!
      ENDIF
!
  101 FORMAT(2I5,E13.4)
  111 FORMAT(3I5,32F10.3)
!
!**********************************************************************C
!
! **  SET TRIDIAGONAL EQUATION COEFFICIENTS
!
      DIFTW=1.4e-7
      DZBED=1./FLOAT(KBH-1)
!
	DO L=2,LA
        DIFFTIMI(L)=HTBED2/((DABEDT+RBADJ*(BELV(L)-BELTMPMIN))*DZBED)**2
	ENDDO
!
!              HTBED2=THERMAL DIFFUSIVITY OF BED M**2/SEC
!
      DO L=2,LA
      ABEDTEM(L,1)=-DIFFTIMI(L)
	ENDDO
!
	DO K=2,KBH-1
	DO L=2,LA
	  ABEDTEM(L,K)=-DIFFTIMI(L)
	ENDDO
	ENDDO
!
      DO L=2,LA
	  THICKBED=DZBED*(DABEDT+RBADJ*(BELV(L)-BELTMPMIN))
	  THICKWAT=DZC(1)*HP(L)
	  RATIO=THICKBED/THICKWAT
	  TMPCOND=2.*DIFTW/(THICKBED*THICKWAT)
        UBED=0.5*( U(L,1)+U(L+1,1) )
        VBED=0.5*( V(L,1)+V(LNC(L),1) )
        USPD=SQRT( UBED*UBED+VBED*VBED )
	  TMPCONV=HTBED1*USPD/THICKBED
	  ABEDTEM(L,KBH)=-RATIO*(TMPCOND+TMPCONV)
	  CBEDTEM(L,KBH-1)=-(TMPCOND+TMPCONV)
	ENDDO
!               HTBED1=DIMENSIONLESS CONVECTION COEFFICIENT
!
	DO K=1,KBH-2
	DO L=2,LA
	  CBEDTEM(L,K)=-DIFFTIMI(L)
	ENDDO
	ENDDO
!
	DO L=2,LA
        CBEDTEM(L,KBH)=0.
	ENDDO
!
      DO K=1,KBH
	DO L=2,LA
	  BBEDTEM(L,K)=(1./DELT)-ABEDTEM(L,K)-CBEDTEM(L,K)
	ENDDO
	ENDDO
!
	DO L=2,LA
        RBEDTEM(L,1)=(1./DELT)*TEMB(L,1)+DIFFTIMI(L)*TBEDIT(L)
      ENDDO
!
      DO K=2,KBH-1
	DO L=2,LA
	  RBEDTEM(L,K)=(1./DELT)*TEMB(L,K)
	ENDDO
	ENDDO
!
	DO L=2,LA
	  RBEDTEM(L,KBH)=(1./DELT)*TEM(L,1)
	ENDDO
!
!**********************************************************************C
!
! **  SOLVE TRIDIAGONAL EQUATIONS
!
      DO L=2,LA
        BETBED(L)=BBEDTEM(L,1)
	  UBEDTEM(L,1)=RBEDTEM(L,1)/BETBED(L)
	ENDDO
!
      DO K=2,KBH
	  DO L=2,LA
	    GBEDTEM(L,K)=CBEDTEM(L,K-1)/BETBED(L)
	    BETBED(L)=BBEDTEM(L,K)-ABEDTEM(L,K)*GBEDTEM(L,K)
	    UBEDTEM(L,K)=(RBEDTEM(L,K)-ABEDTEM(L,K)*UBEDTEM(L,K-1))
     &                /BETBED(L)
	  ENDDO
	ENDDO
!
      DO K=KBH-1,1,-1
	  DO L=2,LA
	    UBEDTEM(L,K)=UBEDTEM(L,K)-GBEDTEM(L,K+1)*UBEDTEM(L,K+1)
	  ENDDO
	ENDDO
!
!**********************************************************************C
!
! **  SET NEW BED AND BOTTOM WATER TEMPERATURES
!
      DO K=1,KBH
	  DO L=2,LA
	    TEMB(L,K)=UBEDTEM(L,K)
	  ENDDO
	ENDDO
!
      DO L=2,LA
        TEM(L,1)=UBEDTEM(L,KBH)
      ENDDO
!
!**********************************************************************C
!
      RETURN
      END