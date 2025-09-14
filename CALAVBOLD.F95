!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALAVBOLD (ISTL)
!
! **  SUBROUTINE CALAV CALCULATES VERTICAL VISCOSITY AND DIFFUSIVITY
! **  USING GLAPERIN ET AL'S MODIFICATION OF THE MELLOR-YAMADA MODEL
! **  (NOTE AV, AB, AND AQ ARE ACTUALLY DIVIDED BY H)
! **  IF ISGA=1 VALUES ARE GEOMETRIC AVERAGES WITH THE PREVIOUS VALUES
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 03/19/2002        John Hamrick       03/19/2002       John Hamrick
!  added drycell bypass and consistent initialization of dry values
!----------------------------------------------------------------------C
!
!**********************************************************************C 
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!     DIMENSION QQI(LCM)
!
!**********************************************************************C 
!
!   SHTOP    =      0.4939
!   SHBOT    =     34.6764
!   SMTOP1   =      0.3933
!   SMTOP2   =      7.8464
!   SMBOT1   =     34.6764
!   SMBOT2   =      6.1272
!   RLIMIT   =      0.0233
!   SHMIN    =      0.0934
!   SMMIN    =      0.1099
!   SHMAX    =      5.2073
!   SMMAX    =      4.9639
!
!  GALPERIN stability functions
!
!
      SFAV0= 0.392010
	SFAV1= 7.760050
	SFAV2=34.676440
	SFAV3= 6.127200
	SFAB0= 0.493928
	SFAB1=34.676440
	RIQMIN=-0.999/SFAB1
!
!
      QQIMAX=1./QQMIN
      AVMAX=AVO
      ABMAX=ABO
      AVMIN=10.
      ABMIN=10.
!     RIQMIN=-1./44.
!030705      RIQMIN=-0.023
!      RIQMAX=0.28
      RAVBTMP=1.
      IF(ISAVBMN.GE.1) RAVBTMP=0.
!
      DO K=1,KC
	DO L=1,LC
	IF(IMASKDRY(L).EQ.1)THEN
      AV(L,K)=AVO*HPI(L)
      AB(L,K)=ABO*HPI(L)
      ENDIF
	ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
!      IF(ISTL.EQ.3)THEN
      IF(ISFAVB.EQ.0)THEN
!
      DO K=1,KS
       DO L=2,LA
	 IF(LMASKDRY(L))THEN
!       QQI(L)=1./(QQMIN+QQ(L,K))
       QQI(L)=1./QQ(L,K)
       QQI(L)=MIN(QQI(L),QQIMAX)
	 ENDIF
       ENDDO
      DO L=2,LA
	 IF(LMASKDRY(L))THEN
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)*(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      IF(ISLLIM.GE.1) RIQ=MIN(RIQ,RIQMAX)
!      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))
!      SFAB=0.5/(1.+36.*RIQ)
!      SFAV=0.3933*(1.+7.8464*RIQ)/((1.+34.6764*RIQ)*(1.+6.1272*RIQ))
!      SFAB=0.4939/(1.+34.6764*RIQ)
!
            SFAV=SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB=SFAB0/(1.+SFAB1*RIQ)
!
      AB(L,K)=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AV(L,K)=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AV(L,K))
       ABMAX=MAX(ABMAX,AB(L,K))
       AVMIN=MIN(AVMIN,AV(L,K))
       ABMIN=MIN(ABMIN,AB(L,K))
      AV(L,K)=AV(L,K)*HPI(L)
      AB(L,K)=SCB(L)*AB(L,K)*HPI(L)
      ENDIF
      ENDDO
      ENDDO
!
      ENDIF
      IF(ISFAVB.EQ.1)THEN
!
      DO K=1,KS
       DO L=2,LA
	 IF(LMASKDRY(L))THEN
!       QQI(L)=1./(QQMIN+QQ(L,K))
       QQI(L)=1./QQ(L,K)
       QQI(L)=MIN(QQI(L),QQIMAX)
	ENDIF
       ENDDO
      DO L=2,LA
	 IF(LMASKDRY(L))THEN
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)*(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      IF(ISLLIM.GE.1) RIQ=MIN(RIQ,RIQMAX)
!      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))
!      SFAB=0.5/(1.+36.*RIQ)
!     SFAV=0.3933*(1.+7.8464*RIQ)/((1.+34.6764*RIQ)*(1.+6.1272*RIQ))
!     SFAB=0.4939/(1.+34.6764*RIQ)
!
            SFAV=SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB=SFAB0/(1.+SFAB1*RIQ)
!
      ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AVTMP)
       ABMAX=MAX(ABMAX,ABTMP)
       AVMIN=MIN(AVMIN,AVTMP)
       ABMIN=MIN(ABMIN,ABTMP)
      AV(L,K)=0.5*(AV(L,K)+AVTMP*HPI(L))
      AB(L,K)=SCB(L)*0.5*(AB(L,K)+ABTMP*HPI(L))
      ENDIF
      ENDDO
      ENDDO
!
      ENDIF
      IF(ISFAVB.EQ.2)THEN
!
      DO K=1,KS
       DO L=2,LA
	 IF(LMASKDRY(L))THEN
!       QQI(L)=1./(QQMIN+QQ(L,K))
       QQI(L)=1./QQ(L,K)
       QQI(L)=MIN(QQI(L),QQIMAX)
	ENDIF
       ENDDO
      DO L=2,LA
	 IF(LMASKDRY(L))THEN
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)*(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      IF(ISLLIM.GE.1) RIQ=MIN(RIQ,RIQMAX)
!      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))
!      SFAB=0.5/(1.+36.*RIQ)
!      SFAV=0.3933*(1.+7.8464*RIQ)/((1.+34.6764*RIQ)*(1.+6.1272*RIQ))
!      SFAB=0.4939/(1.+34.6764*RIQ)
!
            SFAV=SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB=SFAB0/(1.+SFAB1*RIQ)
!
      ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AVTMP)
       ABMAX=MAX(ABMAX,ABTMP)
       AVMIN=MIN(AVMIN,AVTMP)
       ABMIN=MIN(ABMIN,ABTMP)
      AV(L,K)=SQRT(AV(L,K)*AVTMP*HPI(L))
      AB(L,K)=SCB(L)*SQRT(AB(L,K)*ABTMP*HPI(L))
      ENDIF
      ENDDO
      ENDDO
!
      ENDIF
!      ENDIF
!
!      IF(ISTL.EQ.2)THEN!
!      LL=LF+LDM-1!
!      DO K=1,KS
!       DO L=LF,LL
!       QQI(L)=1./(QQMIN+QQ(L,K))
!       ENDDO
!      DO L=LF,LL
!      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)
!     &    *(B(L,K+1)-B(L,K))*QQI(L)
!      RIQ=MAX(RIQ,RIQMIN)
!      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))
!      SFAB=0.5/(1.+36.*RIQ)
!      ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
!      AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
!       AVMAX=MAX(AVMAX,AVTMP)
!       ABMAX=MAX(ABMAX,ABTMP)
!       AVMIN=MIN(AVMIN,AVTMP)
!       ABMIN=MIN(ABMIN,ABTMP)
!      AV(L,K)=SQRT(AV(L,K)*AVTMP*HPI(L))
!      AB(L,K)=SCB(L)*SQRT(AB(L,K)*ABTMP*HPI(L))
!      ENDDO
!      ENDDO
!      ENDDO
!
!      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(ISAVBMN.GE.1)THEN
        DO K=1,KS
        DO L=2,LA
         AVTMP=AVMN*HPI(L)
         ABTMP=ABMN*HPI(L)
         AV(L,K)=MAX(AV(L,K),AVTMP)
         AB(L,K)=MAX(AB(L,K),ABTMP)
        ENDDO
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
!

      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)      
      AVUI(L,K)=2./(AV(L,K)+AV(L-1,K))
      AVVI(L,K)=2./(AV(L,K)+AV(LS,K))
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
!      IF(ISTL.EQ.3)THEN
!
      DO K=2,KS
      DO L=2,LA
!      AQ(L,K)=0.255*(AV(L,K-1)+AV(L,K))
      AQ(L,K)=0.205*(AV(L,K-1)+AV(L,K))
      ENDDO
      ENDDO
!
      DO L=2,LA
!      AQ(L,1)=0.255*AV(L,1)
!      AQ(L,KC)=0.255*AV(L,KS)
      AQ(L,1)=0.205*AV(L,1)
      AQ(L,KC)=0.205*AV(L,KS)
      ENDDO
!
!      ELSE
!
!      DO K=2,KS
!      DO L=2,LA
!      AQTMP=0.255*(AV(L,K-1)+AV(L,K))
!      AQTMP=0.205*(AV(L,K-1)+AV(L,K))
!      AQ(L,K)=SQRT(AQ(L,K)*AQTMP)
!      AQ(L,K)=AQTMP
!      ENDDO
!      ENDDO
!
!      DO L=2,LA
!      AQTMP=0.255*AV(L,1)
!      AQTMP=0.205*AV(L,1)
!      AQ(L,1)=SQRT(AQ(L,1)*AQTMP)
!      AQ(L,1)=AQTMP
!      AQTMP=0.255*AV(L,KS)
!      AQTMP=0.205*AV(L,KS)
!      AQ(L,KC)=SQRT(AQ(L,KC)*AQTMP)
!      AQ(L,KC)=AQTMP
!      ENDDO
!
!      ENDIF
!
!----------------------------------------------------------------------C
!
      RETURN
      END