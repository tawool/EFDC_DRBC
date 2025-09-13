!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALAVB2 (ISTL)
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
!
!----------------------------------------------------------------------C
!
!**********************************************************************C 
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!     DIMENSION QQI(LCM)
!
      DATA ATURB1,ATURB2,TURBC1/0.92,0.74,0.08/
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
      AVMAX=AVO
      ABMAX=ABO
      AVMIN=10.
      ABMIN=10.
!     RIQMIN=-1./44.
!030705      RIQMIN=-0.023
      RAVBTMP=1.
      IF(ISAVBMN.GE.1) RAVBTMP=0.
!
!----------------------------------------------------------------------C
!
!
      DO K=1,KS
       DO L=2,LA
         CTURBB1(L,K)=CTURB
         CTURBB2(L,K)=CTURB2B
       ENDDO
      ENDDO
!
      DO K=1,KS
       DO L=2,LA
         DELBTMP=(B(L,K+1)-B(L,K))*DZIG(K)
         RITMP=-GP*HP(L)*DELBTMP/QQ(L,K)
         IF(RITMP.GT.0.)THEN
           RIQ=DML(L,K)*DML(L,K)*RITMP
           BFUN=EXP(-3.11*RIQ)
           CTURBB1(L,K)=CTURB/(BFUN+1.E-16)
           IF(BBT(L,K).GT.0.)THEN
             TMPVAL=DELBTMP*DELBTMP/(RITMP*BBT(L,K))
             CTURBB2(L,K)=CTURB2B/(1.+0.61*(1.-BFUN)*TMPVAL)
           ENDIF
         ENDIF
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
!      IF(ISTL.EQ.3)THEN
      IF(ISFAVB.EQ.0)THEN
!
      DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=LF+LDM-1
      DO K=1,KS
       DO L=LF,LL
       QQI(L)=1./QQ(L,K)
       ENDDO
      DO L=LF,LL
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)*(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      TMPVAL1=1.-( 6.*ATURB1/CTURBB1(L,K) )
      SBTOP=ATURB2*TMPVAL1
      SBBOT=3.*ATURB2*( 6.*ATURB1+CTURBB2(L,K) )
      SVTOP=ATURB1*( TMPVAL1-3.*TURBC1)
      TMPVAL2=TMPVAL1*( CTURBB2(L,K)-3.*ATURB2 )
      TMPVAL3=-3.*TURBC1*( CTURBB2(L,K)+6.*ATURB1)
      SVTOP2=3.*ATURB2*(TMPVAL2+TMPVAL3)/SVTOP
      SVBOT=9.*ATURB1*ATURB2
      SFAV=SVTOP*(1.+SVTOP2*RIQ)/((1.+SVBOT*RIQ)*(1.+SBBOT*RIQ))
      SFAB=SBTOP/(1.+SBBOT*RIQ)
      AB(L,K)=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AV(L,K)=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AV(L,K))
       ABMAX=MAX(ABMAX,AB(L,K))
       AVMIN=MIN(AVMIN,AV(L,K))
       ABMIN=MIN(ABMIN,AB(L,K))
      AV(L,K)=AV(L,K)*HPI(L)
      AB(L,K)=SCB(L)*AB(L,K)*HPI(L)
      ENDDO
      ENDDO
      ENDDO
!
      ENDIF
      IF(ISFAVB.EQ.1)THEN
!
      DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=LF+LDM-1
      DO K=1,KS
       DO L=LF,LL
       QQI(L)=1./QQ(L,K)
       ENDDO
      DO L=LF,LL
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)*(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      TMPVAL1=1.-( 6.*ATURB1/CTURBB1(L,K) )
      SBTOP=ATURB2*TMPVAL1
      SBBOT=3.*ATURB2*( 6.*ATURB1+CTURBB2(L,K) )
      SVTOP=ATURB1*( TMPVAL1-3.*TURBC1)
      TMPVAL2=TMPVAL1*( CTURBB2(L,K)-3.*ATURB2 )
      TMPVAL3=-3.*TURBC1*( CTURBB2(L,K)+6.*ATURB1)
      SVTOP2=3.*ATURB2*(TMPVAL2+TMPVAL3)/SVTOP
      SVBOT=9.*ATURB1*ATURB2
      SFAV=SVTOP*(1.+SVTOP2*RIQ)/((1.+SVBOT*RIQ)*(1.+SBBOT*RIQ))
      SFAB=SBTOP/(1.+SBBOT*RIQ)
      ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AVTMP)
       ABMAX=MAX(ABMAX,ABTMP)
       AVMIN=MIN(AVMIN,AVTMP)
       ABMIN=MIN(ABMIN,ABTMP)
      AV(L,K)=0.5*(AV(L,K)+AVTMP*HPI(L))
      AB(L,K)=SCB(L)*0.5*(AB(L,K)+ABTMP*HPI(L))
      ENDDO
      ENDDO
      ENDDO
!
      ENDIF
      IF(ISFAVB.EQ.2)THEN
!
      DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=LF+LDM-1
      DO K=1,KS
       DO L=LF,LL
       QQI(L)=1./QQ(L,K)
       ENDDO
      DO L=LF,LL
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)*(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      TMPVAL1=1.-( 6.*ATURB1/CTURBB1(L,K) )
      SBTOP=ATURB2*TMPVAL1
      SBBOT=3.*ATURB2*( 6.*ATURB1+CTURBB2(L,K) )
      SVTOP=ATURB1*( TMPVAL1-3.*TURBC1)
      TMPVAL2=TMPVAL1*( CTURBB2(L,K)-3.*ATURB2 )
      TMPVAL3=-3.*TURBC1*( CTURBB2(L,K)+6.*ATURB1)
      SVTOP2=3.*ATURB2*(TMPVAL2+TMPVAL3)/SVTOP
      SVBOT=9.*ATURB1*ATURB2
      SFAV=SVTOP*(1.+SVTOP2*RIQ)/((1.+SVBOT*RIQ)*(1.+SBBOT*RIQ))
      SFAB=SBTOP/(1.+SBBOT*RIQ)
      ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AVTMP)
       ABMAX=MAX(ABMAX,ABTMP)
       AVMIN=MIN(AVMIN,AVTMP)
       ABMIN=MIN(ABMIN,ABTMP)
      AV(L,K)=SQRT(AV(L,K)*AVTMP*HPI(L))
      AB(L,K)=SCB(L)*SQRT(AB(L,K)*ABTMP*HPI(L))
      ENDDO
      ENDDO
      ENDDO
!
      ENDIF
!      ENDIF
!
!      IF(ISTL.EQ.2)THEN
!
!      DO ND=1,NDM
!      LF=2+(ND-1)*LDM
!      LL=LF+LDM-1
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
      IF(ISTL.EQ.3)THEN
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
      ELSE
!
      DO K=2,KS
      DO L=2,LA
!      AQTMP=0.255*(AV(L,K-1)+AV(L,K))
      AQTMP=0.205*(AV(L,K-1)+AV(L,K))
!      AQ(L,K)=SQRT(AQ(L,K)*AQTMP)
      AQ(L,K)=AQTMP
      ENDDO
      ENDDO
!
      DO L=2,LA
!      AQTMP=0.255*AV(L,1)
      AQTMP=0.205*AV(L,1)
!      AQ(L,1)=SQRT(AQ(L,1)*AQTMP)
      AQ(L,1)=AQTMP
!      AQTMP=0.255*AV(L,KS)
      AQTMP=0.205*AV(L,KS)
!      AQ(L,KC)=SQRT(AQ(L,KC)*AQTMP)
      AQ(L,KC)=AQTMP
      ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------C
!
      RETURN
      END