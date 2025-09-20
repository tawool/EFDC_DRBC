!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALHEAT(ISTL)
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
!  HNR_GHD, 7/2023   INCREASE TEMPERATURE BASED ON HEAT LOAD     
!
!----------------------------------------------------------------------C
!
!**********************************************************************C
!
! **  SUBROUTINE CALHEAT CALCULATES SURFACE AND INTERNAL HEAT SOURCES 
! **  AND SINKS IN THE HEAT (TEM) TRANSPORT EQUATION
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      DIMENSION TEMOLD(LCM,KCM)
!
!**********************************************************************C
!
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      IF(ISTL.EQ.2)THEN
       DELT=DT
       S3TL=0.0
       S2TL=1.0
      ENDIF
!
!      IF(ISTL.EQ.2)THEN
!       DELT=DT
!       S3TL=0.0
!       S2TL=1.0
!      ENDIF
!
      IF(IS2TIM.EQ.1)THEN
        IF(ISDYNSTP.EQ.0)THEN
          DELT=DT
        ELSE
          DELT=DTDYN
        END IF
        S3TL=0.0
        S2TL=1.0
      ENDIF
!
!**********************************************************************C
!
      IF(ISTOPT(2).EQ.1)THEN
!
      DO K=1,KCM
      DO L=2,LA
        TEMOLD(L,K)=TEM(L,K)
      ENDDO
      ENDDO
!
      DO L=2,LA
        TVAR3S(L)=0.
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  ADSORB SW SOLR RAD TO ALL LAYERS AND BED
!
      IF(IASWRAD.EQ.0)THEN
!
       DO L=2,LA
        SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/(1.+0.00412*TEM(L,KC))))
!    &    *(1+1.E-6*PATMT(L)*(4.5+0.0006*TEM(L,KC)*TEM(L,KC)))
        HDEP=MAX(HP(L),0.)
        WNDTMP=WINDST(L)
!X        IF(HP(L).LT.HWET) WNDTMP=0.
        TVAR1S(L,KC)=HDEP*TEM(L,KC)-(DELT*DZIC(KC))*( 1.312E-14*((TEM(L,KC)+273.)**4)             &
                    *(0.39-0.05*SQRT(VPA(L)))*(1.-.8*CLOUDT(L))+5.248E-14*((TEM(L,KC)+273.)**3)   &
                    *(TEM(L,KC)-TATMT(L))+CCNHTT(L)*0.288E-3*WNDTMP*(TEM(L,KC)-TATMT(L))+CLEVAP(L)&
                    *0.445*WNDTMP*(SVPW-VPA(L))/PATMT(L))+(DELT*DZIC(KC))*0.2385E-6*SOLSWRT(L)    &
                    *(FSWRATF*EXP(SWRATNF*HDEP*(Z(KC)-1.))+(1.-FSWRATF)*EXP(SWRATNS*HDEP          &
                    *(Z(KC)-1.))-FSWRATF*EXP(SWRATNF*HDEP*(Z(KC-1)-1.))-(1.-FSWRATF)              &
                    *EXP(SWRATNS*HDEP*(Z(KC-1)-1.)))
       ENDDO

!     &   +1.5*0.288E-6*WNDTMP*(TEM(L,KC)-TATMT(L))
!     &   +1.5*0.445E-3*WNDTMP*(SVPW-VPA(L))/PATMT(L) )

!
       IF(KC.GT.1)THEN
        DO K=1,KS
         DO L=2,LA
          HDEP=MAX(HP(L),0.)
          TVAR1S(L,K)=HDEP*TEM(L,K)+(DELT*DZIC(K))*0.2385E-6*SOLSWRT(L)*(FSWRATF*EXP(SWRATNF*HDEP &
                      *(Z(K)-1.))+(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(K)-1.))-FSWRATF*EXP(SWRATNF    &
                      *HDEP*(Z(K-1)-1.))-(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(K-1)-1.)))
         ENDDO
        ENDDO
       ENDIF
!
! adsorb remaining solar radiation into the bottom layer
!
       IF(DABEDT.GT.0.0)THEN
       DO L=2,LA
	  K=1
        HDEP=MAX(HP(L),0.)
        TVAR1S(L,K)=TVAR1S(L,K)+(DELT*DZIC(K))*0.2385E-6*SOLSWRT(L)*(FSWRATF*EXP(SWRATNF*HDEP     &
                    *(Z(K-1)-1.))+(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(K-1)-1.)))
       ENDDO
       ENDIF

!
!       DO L=2,LA
!        HDEP=MAX(HP(L),0.)
!        UBED=0.5*( U(L,1)+U(L+1,1) )
!        VBED=0.5*( V(L,1)+V(LNC(L),1) )
!        USPD=SQRT( UBED*UBED+VBED*VBED )
!        TMPVAL=(HTBED1*USPD+HTBED2)*(TEM(L,1)-TEMB(L,KBHM))
!        TVAR1S(L,1)=TVAR1S(L,1)-DELT*DZIC(1)*TMPVAL
!        TEMB(L,KBHM)=TEMB(L,KBHM)+DELT*TMPVAL/DABEDT
!     &        +(DELT/DABEDT)*0.2385E-6*SOLSWRT(L)*(
!     &        +FSWRATF*EXP(SWRATNF*HDEP*(Z(0)-1.))
!     &        +(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(0)-1.)) )
!       ENDDO
!
       DO K=1,KC
        DO L=2,LA
         IF(HP(L).GT.0.) TEM(L,K)=HPI(L)*TVAR1S(L,K)
	   TEM(L,K)=MAX(TEM(L,K),0.0)
        ENDDO
       ENDDO
!
       IF(ISICE.EQ.1)THEN
       DO K=1,KC
        DO L=2,LA
         TEM(L,K)=RICECOVL(L)*TEMPICE+(1.-RICECOVL(L))*TEM(L,K)
        ENDDO
       ENDDO
	 ENDIF
!
       IF(ISDRY.GT.0)THEN
         DO K=1,KC
           DO L=2,LA
             IF(IMASKDRY(L).EQ.1.) TEM(L,K)=TATMT(L)
           ENDDO
         ENDDO
!         DO L=2,LA
!           IF(IMASKDRY(L).EQ.1.) TEMB(L,KBHM)=TATMT(L)
!         ENDDO
       ENDIF
!
       NTSTBCM1=NTSTBC-1
	 IF(DABEDT.GT.0.0) THEN
         IF(IS2TIM.EQ.1) CALL CALHEATB(ISTL)	
	   IF(IS2TIM.EQ.0) THEN
	     IF(NCTBC.EQ.NTSTBCM1) CALL CALHEATB(ISTL)	
         ENDIF
       ENDIF
!
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  ADSORB SW SOLR RAD TO SURFACE LAYER
!
      IF(IASWRAD.EQ.1)THEN
!
       DO L=2,LA
        SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/(1.+0.00412*TEM(L,KC))))
!    &    *(1+1.E-6*PATMT(L)*(4.5+0.0006*TEM(L,KC)*TEM(L,KC)))
        HDEP=MAX(HP(L),0.)
        WNDTMP=WINDST(L)
!X        IF(HP(L).LT.HWET) WNDTMP=0.
        TVAR1S(L,KC)=HDEP*TEM(L,KC)-(DELT*DZIC(KC))*( 1.312E-14*((TEM(L,KC)+273.)**4)*(0.39-0.05  &
                     *SQRT(VPA(L)))*(1.-.8*CLOUDT(L))+5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)  &
                     -TATMT(L))+CCNHTT(L)*0.288E-3*WNDTMP*(TEM(L,KC)-TATMT(L))+CLEVAP(L)*0.445    &
                     *WNDTMP*(SVPW-VPA(L))/PATMT(L))+(DELT*DZIC(KC))*0.2385E-6*SOLSWRT(L)
       ENDDO
!
       DO L=2,LA
         IF(HP(L).GT.0.) TEM(L,KC)=HPI(L)*TVAR1S(L,KC)
	   TEM(L,KC)=MAX(TEM(L,KC),0.0)
       ENDDO
!
       IF(ISICE.EQ.1)THEN
       DO K=1,KC
        DO L=2,LA
         TEM(L,K)=RICECOVL(L)*TEMPICE+(1.-RICECOVL(L))*TEM(L,K)
        ENDDO
       ENDDO
	 ENDIF
!
       IF(ISDRY.GT.0)THEN
         DO K=1,KC
           DO L=2,LA
             IF(IMASKDRY(L).EQ.1.) TEM(L,K)=TATMT(L)
           ENDDO
         ENDDO
!         DO L=2,LA
!           IF(IMASKDRY(L).EQ.1.) TEMB(L,KBHM)=TATMT(L)
!         ENDDO
       ENDIF
C      
       NTSTBCM1=NTSTBC-1
	 IF(DABEDT.GT.0.0) THEN
         IF(IS2TIM.EQ.1) CALL CALHEATB(ISTL)	
	   IF(IS2TIM.EQ.0) THEN
	     IF(NCTBC.EQ.NTSTBCM1) CALL CALHEATB(ISTL)	
         ENDIF
       ENDIF
!
      ENDIF
!
! **  ENDIF ISOPT(2) EQ 1
!
      ENDIF
!
  600 FORMAT(4I5,2E12.4)
!
!**********************************************************************C
!
      IF(ISTOPT(2).EQ.2)THEN
!
! ** IMPLEMENT EXTERNALLY SPECIFIED EQUILIBRIUM TEMPERATURE FROMULATION
!
      TMPKC=DELT/DZC(KC)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TEM(L,KC)=TEM(L,KC)-TMPKC*CLOUDT(L)*HPI(L)*(TEM(L,KC)-TATMT(L))
	  TEM(L,KC)=MAX(TEM(L,KC),0.0)
       ENDDO
      ENDDO
!
       IF(ISICE.EQ.1)THEN
       DO K=1,KC
        DO L=2,LA
         TEM(L,K)=RICECOVL(L)*TEMPICE+(1.-RICECOVL(L))*TEM(L,K)
        ENDDO
       ENDDO
	 ENDIF
!
      ENDIF
!
!**********************************************************************C
!
      IF(ISTOPT(2).EQ.3)THEN
!
! ** IMPLEMENT CONSTANT COEFFICIENT EQUILIBRIUM TEMPERATURE FROMULATION
!
      DTHEQT=DELT*HEQT*FLOAT(KC)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TEM(L,KC)=TEM(L,KC)-DTHEQT*HPI(L)*(TEM(L,KC)-TEMO)
	  TEM(L,KC)=MAX(TEM(L,KC),0.0)
       ENDDO
      ENDDO
!
       IF(ISICE.EQ.1)THEN
       DO K=1,KC
        DO L=2,LA
         TEM(L,K)=RICECOVL(L)*TEMPICE+(1.-RICECOVL(L))*TEM(L,K)
        ENDDO
       ENDDO
	 ENDIF
!
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  START EXTERNAL TIME SERIES LOAD FOR TEMPERATURE    HNR_GHD 7/2023

      DO NL=1,NLIJ(2)
        L=LLDS(NL,2)
        NLSTMP=NSERL(NL,2)
        DO K=1,KC
          WSH=CALWSPHEAT(TEM(L,K),SAL(L,K))
          RHOL=FUNDEN(SAL(L,K),0.0,TEM(L,K))
          DENOM=WSH*RHOL*DXYP(L)*H1P(L)*DZC(K)
          DELTEMP(L,K)=SERLT(K,NLSTMP,2)*dt/DENOM                  !LOAD SERLT IN J/S  
          TEM(L,K)=TEM(L,K)+DELTEMP(L,K)
        END DO
      END DO

! **  END EXTERNAL TIME SERIES LOAD FOR TEMPERATURE    HNR_GHD 7/2023
!
!----------------------------------------------------------------------C
!
!**********************************************************************C
!
      RETURN
      END