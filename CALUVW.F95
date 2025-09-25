!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALUVW (ISTL,IS2TL)
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
! **  CALCULATE THE INTERNAL SOLUTION AT TIME LEVEL (N+1)
! **  THE VALUE OF ISTL INDICATES THE NUMBER OF TIME LEVELS IN THE STEP
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
!
!      IF(ISTL.EQ.3)THEN
!       DELT=DT2
!       DELTD2=DT
!      ELSE
!       DELT=DT
!       DELTD2=0.5*DT
!      ENDIF
!
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT2
        DELTD2=DT
        IF(ISTL.EQ.2)THEN
          DELT=DT
          DELTD2=0.5*DT
        ENDIF
        DELTI=1./DELT
      ELSE
        DELT=DTDYN
        DELTD2=0.5*DTDYN
        DELTI=1./DELT
      ENDIF
!
      IF(ISCDMA.EQ.3)THEN
        DELT=0.5*DELT
        DELTD2=0.5*DELTD2
        DELTI=1./DELT
      ENDIF
!
      IF(ISCDMA.EQ.4)THEN
        DELT=0.5*DELT
        DELTD2=0.5*DELTD2
        DELTI=1./DELT
      ENDIF
!
      IF(KC.EQ.1) GOTO 30
!
!**********************************************************************C
!
! **  CALCULATE BOTTOM FRICTION COEFFICIENT
!
!----------------------------------------------------------------------C
!
      IF(ISTL.EQ.3)THEN
!
       DO L=2,LA
        RCX(L)=AVCON1/H1U(L)+STBX(L)*SQRT(U1(L,1)*U1(L,1)+V1U(L)*V1U(L))
        RCY(L)=AVCON1/H1V(L)+STBY(L)*SQRT(U1V(L)*U1V(L)+V1(L,1)*V1(L,1))
       ENDDO
!
!GVCDIAG
!      IF(N.LE.4)THEN
!      DO L=2,LA
!	  KBTMPU=1
!	  WRITE(8,8881)N,IL(L),JL(L),KBTMPU,RCX(L),H1U(L),STBX(L),
!     &               U1(L,1),V1U(L)
!      ENDDO
!      ENDIF
! 8881 FORMAT('KBU,RCX,H1U,ST,U,V SIG ',4I5,10E14.6)
!GVCDIAG
!
!
!      IF(ISVEG.GE.1.AND.KC.GT.1)
!      DO ND=1,NDM
!       LF=2+(ND-1)*LDM
!       LL=LF+LDM-1
!       DO L=LF,LL
!        RCX(L)=RCX(L)+0.5*FXVEG(L,1)*SQRT(U1(L,1)*U1(L,1)
!     &        +V1U(L)*V1U(L))
!        RCY(L)=RCY(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*U1V(L)
!     &        +V1(L,1)*V1(L,1))
!       ENDDO
!      ENDDO
!      ENDIF
!
      ELSE
!
       IF(AVCON1.LT.0.00001)THEN
       DO L=2,LA
        Q1=SQRT(U1(L,1)*U1(L,1)+V1U(L)*V1U(L))
        Q2=SQRT(U(L,1)*U(L,1)+VU(L)*VU(L))
        RCX(L)=STBX(L)*SQRT(Q1*Q2)
        Q1=SQRT(U1V(L)*U1V(L)+V1(L,1)*V1(L,1))
        Q2=SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
        RCY(L)=STBY(L)*SQRT(Q1*Q2)
       ENDDO
	 ELSE
       DO L=2,LA
        Q1=SQRT(U1(L,1)*U1(L,1)+V1U(L)*V1U(L))
        Q2=SQRT(U(L,1)*U(L,1)+VU(L)*VU(L))
        RCX(L)=AVCON1/SQRT(H1U(L)*HU(L))+STBX(L)*SQRT(Q1*Q2)
        Q1=SQRT(U1V(L)*U1V(L)+V1(L,1)*V1(L,1))
        Q2=SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
        RCY(L)=AVCON1/SQRT(H1V(L)*HV(L))+STBY(L)*SQRT(Q1*Q2)
       ENDDO
	 ENDIF
!
!      IF(ISVEG.GE.1.AND.KC.GT.1)
!      DO ND=1,NDM
!       LF=2+(ND-1)*LDM
!       LL=LF+LDM-1
!       DO L=LF,LL
!        Q1=SQRT(U1(L,1)*U1(L,1)+V1U(L)*V1U(L))
!        Q2=SQRT(U(L,1)*U(L,1)+VU(L)*VU(L))
!        RCX(L)=RCX(L)+0.5*FXVEG(L,1)*SQRT(Q1*Q2)
!        Q1=SQRT(U1V(L)*U1V(L)+V1(L,1)*V1(L,1))
!        Q2=SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
!        RCY(L)=RCY(L)+0.5*FYVEG(L,1)*SQRT(Q1*Q2)
!       ENDDO
!      ENDDO
!      ENDIF
!
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE THE U AND V SHEARS
!
!----------------------------------------------------------------------C
!
!
!GVCDIAG
!      IF(N.LE.4)THEN
!      DO L=2,LA
!	  WRITE(8,881)N,IL(L),JL(L),RCX(L),UHE(L),(DU(L,K),K=1,KS)
!      ENDDO
!      ENDIF
!  881 FORMAT('R,UHE,DU SIG ',3I5,10E14.6)
!GVCDIAG

      RCDZM=CDZM(1)*DELTI
      RCDZU=CDZU(1)
      RCDZL=CDZL(1)
       DO L=2,LA
        CMU=1.+RCDZM*HU(L)*AVUI(L,1)
        CMV=1.+RCDZM*HV(L)*AVVI(L,1)
        EU=1./CMU
        EV=1./CMV
        CU1(L,1)=RCDZU*EU
        CU2(L,1)=RCDZU*EV
        DU(L,1)=(DU(L,1)-RCDZL*RCX(L)*UHE(L)*HUI(L))*EU
        DV(L,1)=(DV(L,1)-RCDZL*RCY(L)*VHE(L)*HVI(L))*EV
        UUU(L,1)=EU
        VVV(L,1)=EV
       ENDDO
!
       DO K=2,KS
        RCDZM=CDZM(K)*DELTI
        RCDZU=CDZU(K)
        RCDZL=CDZL(K)
        DO L=2,LA
         CMU=1.+RCDZM*HU(L)*AVUI(L,K)
         CMV=1.+RCDZM*HV(L)*AVVI(L,K)
         EU=1./(CMU-RCDZL*CU1(L,K-1))
         EV=1./(CMV-RCDZL*CU2(L,K-1))
         CU1(L,K)=RCDZU*EU
         CU2(L,K)=RCDZU*EV
         DU(L,K)=(DU(L,K)-RCDZL*DU(L,K-1))*EU
         DV(L,K)=(DV(L,K)-RCDZL*DV(L,K-1))*EV
         UUU(L,K)=-RCDZL*UUU(L,K-1)*EU
         VVV(L,K)=-RCDZL*VVV(L,K-1)*EV
        ENDDO
       ENDDO
!
       DO K=KS-1,1,-1
        DO L=2,LA
         DU(L,K)=DU(L,K)-CU1(L,K)*DU(L,K+1)
         DV(L,K)=DV(L,K)-CU2(L,K)*DV(L,K+1)
         UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)
         VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)
        ENDDO
       ENDDO
!
!GVCDIAG
!      IF(N.LE.4)THEN
!      DO L=2,LA
!	  WRITE(8,884)N,IL(L),JL(L),(DU(L,K),K=1,KS)
!      ENDDO
!      DO L=2,LA
!	  WRITE(8,885)N,IL(L),JL(L),(UUU(L,K),K=1,KS)
!      ENDDO
!      ENDIF
!  884 FORMAT('BKSUB DU  SIG ',3I5,10E14.6)
!  885 FORMAT('BKSUB UUU SIG ',3I5,10E14.6)
!GVCDIAG
!

       DO L=2,LA
        AAU(L)=0.
        AAV(L)=0.
        BBU(L)=1.
        BBV(L)=1.
       ENDDO
!
       DO K=1,KS
        RCDZR=CDZR(K)
        DO L=2,LA
         CRU=RCDZR*RCX(L)*AVUI(L,K)
         CRV=RCDZR*RCY(L)*AVVI(L,K)
         AAU(L)=AAU(L)+CRU*DU(L,K)
         AAV(L)=AAV(L)+CRV*DV(L,K)
         BBU(L)=BBU(L)+CRU*UUU(L,K)
         BBV(L)=BBV(L)+CRV*VVV(L,K)
        ENDDO
       ENDDO
!
       DO L=2,LA
        AAU(L)=AAU(L)/BBU(L)
        AAV(L)=AAV(L)/BBV(L)
       ENDDO
!
      DO K=1,KS
       RDZG=DZG(K)
       DO L=2,LA
        DU(L,K)=RDZG*HU(L)*AVUI(L,K)*(DU(L,K)-AAU(L)*UUU(L,K))
        DV(L,K)=RDZG*HV(L)*AVVI(L,K)*(DV(L,K)-AAV(L)*VVV(L,K))
       ENDDO
      ENDDO
!
!GVCDIAG
!      IF(N.LE.4)THEN
!      DO L=2,LA
!	  WRITE(8,882)N,IL(L),JL(L),(DU(L,K),K=1,KS)
!      ENDDO
!      ENDIF
!  882 FORMAT('SMW DU SIG ',3I5,10E14.6)
!GVCDIAG
!**********************************************************************C
!
! **  CALCULATED U AND V
!
! **  DUSUM+UHE=UHE, DVSUM+VHE=VHE
!
!----------------------------------------------------------------------C
!
       DO K=1,KS
        RCDZD=CDZD(K)
        DO L=2,LA
         UHE(L)=UHE(L)+RCDZD*DU(L,K)
         VHE(L)=VHE(L)+RCDZD*DV(L,K)
        ENDDO
       ENDDO
!
       DO L=2,LA
        UHDY(L,KC)=UHE(L)*SUB(L)
        VHDX(L,KC)=VHE(L)*SVB(L)
       ENDDO
!
       DO K=KS,1,-1
        DO L=2,LA
         UHDY(L,K)=UHDY(L,K+1)-DU(L,K)*SUB(L)
         VHDX(L,K)=VHDX(L,K+1)-DV(L,K)*SVB(L)
        ENDDO
       ENDDO
!
      DO K=1,KC
       DO L=2,LA
        U(L,K)=UHDY(L,K)*HUI(L)
        V(L,K)=VHDX(L,K)*HVI(L)
        UHDY(L,K)=UHDY(L,K)*DYU(L)
        VHDX(L,K)=VHDX(L,K)*DXV(L)
       ENDDO
      ENDDO
!
!GVCDIAG
!      IF(N.LE.4)THEN
!      DO L=2,LA
!	  WRITE(8,883)N,IL(L),JL(L),(UHDY(L,K),K=1,KC)
!      ENDDO
!      ENDIF
!  883 FORMAT('UHDY SIG ',3I5,10E14.6)
!GVCDIAG
!
! **  ADD ADJUSTMENT TO 3D HORIZONTAL TRANSPORT
!
      DO L=2,LA
        TVAR3E(L)=0.
        TVAR3N(L)=0.
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
        TVAR3E(L)=TVAR3E(L)+UHDY(L,K)*DZC(K)
        TVAR3N(L)=TVAR3N(L)+VHDX(L,K)*DZC(K)
       ENDDO
      ENDDO
!
      UERMX=-1.E+12
      UERMN=1.E+12
      VERMX=-1.E+12
      VERMN=1.E+12
      DO L=2,LA
        TVAR3E(L)=TVAR3E(L)-UHDYE(L)
        TVAR3N(L)=TVAR3N(L)-VHDXE(L)
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
        UHDY(L,K)=UHDY(L,K)-TVAR3E(L)*DZIC(K)
        VHDX(L,K)=VHDX(L,K)-TVAR3N(L)*DZIC(K)
       ENDDO
      ENDDO
!
! **  RESET VELOCITIES
!
      DO L=2,LA
       UHE(L)=0.
       VHE(L)=0.
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
        UHE(L)=UHE(L)+UHDY(L,K)*DZC(K)
        VHE(L)=VHE(L)+VHDX(L,K)*DZC(K)
        U(L,K)=UHDY(L,K)*HUI(L)
        V(L,K)=VHDX(L,K)*HVI(L)
       ENDDO
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
        U(L,K)=U(L,K)*DYIU(L)
        V(L,K)=V(L,K)*DXIV(L)
       ENDDO
      ENDDO
!
      DO L=2,LA
       UHE(L)=UHE(L)*DYIU(L)
       VHE(L)=VHE(L)*DXIV(L)
      ENDDO
!
! **  UNCOMMENT BELOW TO WRITE CONTINUITY DIAGNOSITCS
!
!CDIG      DO L=2,LA
!CDIG        TVAR3E(L)=0.
!CDIG        TVAR3N(L)=0.
!CDIG      ENDDO
!
!CDIG      DO K=1,KC
!CDIG       DO L=2,LA
!CDIG        TVAR3E(L)=TVAR3E(L)+UHDY(L,K)*DZC(K)
!CDIG        TVAR3N(L)=TVAR3N(L)+VHDX(L,K)*DZC(K)
!CDIG       ENDDO
!CDIG      ENDDO
!
!CDIG      UERMX=-1.E+12
!CDIG      UERMN=1.E+12
!CDIG      VERMX=-1.E+12
!CDIG      VERMN=1.E+12
!CDIG      DO L=2,LA
!CDIG        TVAR3E(L)=TVAR3E(L)-UHDYE(L)
!CDIG        TVAR3N(L)=TVAR3N(L)-VHDXE(L)
!CDIG        IF(TVAR3E(L).GT.UERMX)THEN
!CDIG          UERMX=TVAR3E(L)
!CDIG          LUMX=L
!CDIG        ENDIF
!CDIG        IF(TVAR3E(L).LT.UERMN)THEN
!CDIG          UERMN=TVAR3E(L)
!CDIG          LUMN=L
!CDIG        ENDIF
!CDIG        IF(TVAR3N(L).GT.VERMX)THEN
!CDIG          VERMX=TVAR3N(L)
!CDIG          LVMX=L
!CDIG        ENDIF
!CDIG        IF(TVAR3N(L).LT.VERMN)THEN
!CDIG          VERMN=TVAR3N(L)
!CDIG          LVMN=L
!CDIG        ENDIF
!CDIG      ENDDO
!
!CDIG      WRITE(6,6661)IL(LUMX),JL(LUMX),UERMX
!CDIG      WRITE(6,6662)IL(LUMN),JL(LUMN),UERMN
!CDIG      WRITE(6,6663)IL(LVMX),JL(LVMX),VERMX
!CDIG      WRITE(6,6664)IL(LVMN),JL(LVMN),VERMN
!
 6661 FORMAT(' I,J,UHDYERMX = ',2I5,E14.5)
 6662 FORMAT(' I,J,UHDYERMN = ',2I5,E14.5)
 6663 FORMAT(' I,J,VHDYERMX = ',2I5,E14.5)
 6664 FORMAT(' I,J,VHDYERMX = ',2I5,E14.5)
!
!**********************************************************************C
!
! **  CALCULATE W
!
!----------------------------------------------------------------------C
!
      IF(ISTL.EQ.3)THEN
!
       DO L=2,LA
        TVAR3E(L)=UHDYE(L+1   )
        TVAR3N(L)=VHDXE(LNC(L))
        TVAR3W(L)=UHDY2E(L+1   )
        TVAR3S(L)=VHDX2E(LNC(L))
       ENDDO
!
       DO K=1,KC
        DO L=2,LA
         TVAR1E(L,K)=UHDY(L+1   ,K)
         TVAR1N(L,K)=VHDX(LNC(L),K)
         TVAR1W(L,K)=UHDY2(L+1   ,K)
         TVAR1S(L,K)=VHDX2(LNC(L),K)
        ENDDO
       ENDDO
!
       DO K=1,KS
        DO L=2,LA
         LN=LNC(L)
               
         W(L,K)=SWB(L)*( W(L,K-1)-W2(L,K)+W2(L,K-1)                                               &
            -DZC(K)*(TVAR1E(L,K)-UHDY(L,K)-TVAR3E(L)+UHDYE(L)                                     &
            +TVAR1W(L,K)-UHDY2(L,K)-TVAR3W(L)+UHDY2E(L)                                          &
            +TVAR1N(L,K)-VHDX(L,K)-TVAR3N(L)+VHDXE(L)                                            &
            +TVAR1S(L,K)-VHDX2(L,K)-TVAR3S(L)+VHDX2E(L) )*DXYIP(L) )                             &
            +2.*SWB(L)*( QSUM(L,K)-DZC(K)*QSUME(L) )*DXYIP(L)
        ENDDO
       ENDDO
!
! **  UNCOMMENT BELOW FOR CONTINUITY DIAGNOSTICS
!
!CDIG      WSFMAX=-1.E+12
!CDIG      WSFMIN=1.E+12
!CDIG      SURFOUT=0.
!CDIG      K=KC
!CDIG      DO L=2,LA
!CDIG      IF(SWB(L).GT.0.5)THEN
!CDIG         WSURF=SWB(L)*( W(L,K-1)-W2(L,K)+W2(L,K-1)
!CDIG     &  -DZC(K)*(TVAR1E(L,K)-UHDY(L,K)-TVAR3E(L)+UHDYE(L)
!CDIG     &       +TVAR1W(L,K)-UHDY2(L,K)-TVAR3W(L)+UHDY2E(L)
!CDIG     &       +TVAR1N(L,K)-VHDX(L,K)-TVAR3N(L)+VHDXE(L)
!CDIG     &       +TVAR1S(L,K)-VHDX2(L,K)-TVAR3S(L)+VHDX2E(L) )*DXYIP(L) )
!CDIG     &       +2.*SWB(L)*( QSUM(L,K)-DZC(K)*QSUME(L) )*DXYIP(L)
!CDIG      SURFOUT=SURFOUT+WSURF*DXYP(L)
!CDIG      IF(WSURF.GT.WSFMAX)THEN
!CDIG        WSFMAX=WSURF
!CDIG        IMAX=IL(L)
!CDIG        JMAX=JL(L)
!CDIG      ENDIF
!CDIG      IF(WSURF.LT.WSFMIN)THEN
!CDIG        WSFMIN=WSURF
!CDIG        IMIN=IL(L)
!CDIG        JMIN=JL(L)
!CDIG      ENDIF
!CDIG      ENDIF
!CDIG      ENDDO
!
!CDIG      L=LIJ(IMAX,JMAX)
!CDIG      WSFMAX=WSFMAX*DXYP(L)
!CDIG      L=LIJ(IMIN,JMIN)
!CDIG      WSFMIN=WSFMIN*DXYP(L)
!CDIG      WRITE(6,601)IMAX,JMAX,WSFMAX
!CDIG      WRITE(6,602)IMIN,IMAX,WSFMIN
!CDIG      WRITE(6,603)SURFOUT
!
      ENDIF
!
      IF(ISTL.EQ.2)THEN
!
       DO L=2,LA
        TVAR3E(L)=UHDYE(L+1   )
        TVAR3N(L)=VHDXE(LNC(L))
        TVAR3W(L)=UHDY1E(L+1   )
        TVAR3S(L)=VHDX1E(LNC(L))
       ENDDO
!
       DO K=1,KC
        DO L=2,LA
         TVAR1E(L,K)=UHDY(L+1   ,K)
         TVAR1N(L,K)=VHDX(LNC(L),K)
         TVAR1W(L,K)=UHDY1(L+1   ,K)
         TVAR1S(L,K)=VHDX1(LNC(L),K)
        ENDDO
       ENDDO
!
       DO K=1,KS
        DO L=2,LA
         LN=LNC(L)
         W(L,K)=SWB(L)*( W(L,K-1)-W1(L,K)+W1(L,K-1)                                               &
       -DZC(K)*(TVAR1E(L,K)-UHDY(L,K)-TVAR3E(L)+UHDYE(L)                                          &
            +TVAR1W(L,K)-UHDY1(L,K)-TVAR3W(L)+UHDY1E(L)                                           &
            +TVAR1N(L,K)-VHDX(L,K)-TVAR3N(L)+VHDXE(L)                                             &
            +TVAR1S(L,K)-VHDX1(L,K)-TVAR3S(L)+VHDX1E(L) )*DXYIP(L) )                              &
            +2.*SWB(L)*( QSUM(L,K)-DZC(K)*QSUME(L) )*DXYIP(L)
        ENDDO
       ENDDO
!
!CDIG      WSFMAX=-1.E+12
!CDIG      WSFMIN=1.E+12
!CDIG      SURFOUT=0.
!CDIG      K=KC
!CDIG      DO L=2,LA
!CDIG      IF(SWB(L).GT.0.5)THEN
!CDIG         WSURF=SWB(L)*( W(L,K-1)-W1(L,K)+W1(L,K-1)
!CDIG     &  -DZC(K)*(TVAR1E(L,K)-UHDY(L,K)-TVAR3E(L)+UHDYE(L)
!CDIG     &       +TVAR1W(L,K)-UHDY1(L,K)-TVAR3W(L)+UHDY1E(L)
!CDIG     &       +TVAR1N(L,K)-VHDX(L,K)-TVAR3N(L)+VHDXE(L)
!CDIG     &       +TVAR1S(L,K)-VHDX1(L,K)-TVAR3S(L)+VHDX1E(L) )*DXYIP(L) )
!CDIG     &       +2.*SWB(L)*( QSUM(L,K)-DZC(K)*QSUME(L) )*DXYIP(L)
!CDIG      SURFOUT=SURFOUT+WSURF*DXYP(L)
!CDIG      IF(WSURF.GT.WSFMAX)THEN
!CDIG        WSFMAX=WSURF
!CDIG        IMAX=IL(L)
!CDIG        JMAX=JL(L)
!CDIG      ENDIF
!CDIG      IF(WSURF.LT.WSFMIN)THEN
!CDIG        WSFMIN=WSURF
!CDIG        IMIN=IL(L)
!CDIG        JMIN=JL(L)
!CDIG      ENDIF
!CDIG      ENDIF
!CDIG      ENDDO
!
!CDIG      L=LIJ(IMAX,JMAX)
!CDIG      WSFMAX=WSFMAX*DXYP(L)
!CDIG      L=LIJ(IMIN,JMIN)
!CDIG      WSFMIN=WSFMIN*DXYP(L)
!CDIG      WRITE(6,601)IMAX,JMAX,WSFMAX
!CDIG      WRITE(6,602)IMIN,IMAX,WSFMIN
!CDIG      WRITE(6,603)SURFOUT
!
      ENDIF
!
  601 FORMAT(' IMAX,JMAX,QWSFMAX = ',2I5,E14.5)
  602 FORMAT(' IMIN,JMIN,QWSFMIN = ',2I5,E14.5)
  603 FORMAT(' TOTAL SURF Q ERR = ',E14.5)
!
!**********************************************************************C
!
! **  CALCULATE U AND V ON OPEN BOUNDARIES
!
   30 CONTINUE
!
!----------------------------------------------------------------------C
!
      IF(IS1DCHAN.EQ.0)THEN
!
      DO K=1,KC
      DO LL=1,NCBS
	IF(ISPBS(LL).LE.1)THEN
      L=LCBS(LL)
      LN=LNC(L)
      LNN=LNC(LN)
      IF(LN.NE.LC)THEN
        VHDX(LN,K)=VHDX(LNN,K)-VHDXE(LNN)+VHDXE(LN)
        V(LN,K)=VHDX(LN,K)/(HV(LN)*DXV(LN))
       ELSE
        VHDX(LN,K)=0.
        V(LN,K)=0.
      ENDIF
	ENDIF
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBW
	IF(ISPBW(LL).LE.1)THEN
      L=LCBW(LL)
      LP=L+1
      LPP=L+2
      IF(LP.NE.LC)THEN
        UHDY(LP,K)=UHDY(LPP,K)-UHDYE(LPP)+UHDYE(LP)
        U(LP,K)=UHDY(LP,K)/(HU(LP)*DYU(LP))
       ELSE
        UHDY(LP,K)=0.
        U(LP,K)=0.
      ENDIF
	ENDIF
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBE
	IF(ISPBE(LL).LE.1)THEN
      L=LCBE(LL)
      UHDY(L,K)=UHDY(L-1,K)-UHDYE(L-1)+UHDYE(L)
      U(L,K)=UHDY(L,K)/(HU(L)*DYU(L))
	ENDIF
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBN
	IF(ISPBN(LL).LE.1)THEN
      L=LCBN(LL)
      LS=LSC(L)
      VHDX(L,K)=VHDX(LS,K)-VHDXE(LS)+VHDXE(L)
      V(L,K)=VHDX(L,K)/(HV(L)*DXV(L))
	ENDIF
      ENDDO
      ENDDO
!
      ENDIF
!
      IF(IS1DCHAN.GE.1)THEN
!
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      LNN=LNC(LN)
      VHDX(LN,K)=VHDX(LNN,K)-VHDXE(LNN)+VHDXE(LN)
      V(LN,K)=VHDX(LN,K)/FADXV(LN)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)
      LN=LNC(L)
      UHDY(L+1,K)=UHDY(L+2,K)-UHDYE(L+2)+UHDYE(L+1)
      U(L+1,K)=UHDY(L+1,K)/FADYU(L+1)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      UHDY(L,K)=UHDY(L-1,K)-UHDYE(L-1)+UHDYE(L)
      U(L,K)=UHDY(L,K)/FADYU(L)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      VHDX(L,K)=VHDX(LS,K)-VHDXE(LS)+VHDXE(L)
      V(L,K)=VHDX(L,K)/FADXV(L)
      ENDDO
      ENDDO
!
      ENDIF
!
!**********************************************************************C
!
! **  CALCULATE AVERAGE CELL FACE TRANSPORTS FOR SALT, TEMPERATURE AND
! **  SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2
!
!----------------------------------------------------------------------C
!
      IF(ISCDMA.LE.4)THEN
      IF(ISTL.EQ.2)THEN
!
        DO K=1,KC
         DO L=2,LA
          UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY1(L,K))
          VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX1(L,K))
          U2(L,K)=0.5*(U(L,K)+U1(L,K))
          V2(L,K)=0.5*(V(L,K)+V1(L,K))
          W2(L,K)=0.5*(W(L,K)+W1(L,K))
         ENDDO
        ENDDO
!
       ELSE
!
        DO K=1,KC
         DO L=2,LA
!         DU(L,K)=0.25*(UHDY(L,K)-UHDY2(L,K))
!         DV(L,K)=0.25*(VHDX(L,K)-VHDX2(L,K))
          UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY2(L,K))
          VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX2(L,K))
          U2(L,K)=0.5*(U(L,K)+U2(L,K))
          V2(L,K)=0.5*(V(L,K)+V2(L,K))
!         DW(L,K)=0.25*(W(L,K)-W2(L,K))
          W2(L,K)=0.5*(W(L,K)+W2(L,K))
         ENDDO
        ENDDO
!
      ENDIF
      ENDIF
!
!
      IF(ISCDMA.GE.5)THEN
      IF(ISTL.EQ.2)THEN
!
        DO K=1,KC
         DO L=2,LA
          UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY1(L,K))
          VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX1(L,K))
          U2(L,K)=0.5*(U(L,K)+U1(L,K))
          V2(L,K)=0.5*(V(L,K)+V1(L,K))
          W2(L,K)=0.5*(W(L,K)+W1(L,K))
         ENDDO
        ENDDO
!
       ELSE
!
        DO K=1,KC
         DO L=2,LA
!         DU(L,K)=0.25*(UHDY(L,K)-UHDY2(L,K))
!         DV(L,K)=0.25*(VHDX(L,K)-VHDX2(L,K))
          UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY2(L,K))
          VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX2(L,K))
          U2(L,K)=0.5*(U(L,K)+U2(L,K))
          V2(L,K)=0.5*(V(L,K)+V2(L,K))
!         DW(L,K)=0.25*(W(L,K)-W2(L,K))
          W2(L,K)=0.5*(W(L,K)+W2(L,K))
         ENDDO
        ENDDO
!
      ENDIF
      ENDIF
!
!
      IF(ISWVSD.GE.1)THEN
!
        DO K=1,KC
         DO L=2,LA
          UHDY2(L,K)=UHDY2(L,K)+DYU(L)*UVPT(L,K)
          VHDX2(L,K)=VHDX2(L,K)+DXV(L)*VVPT(L,K)
          U2(L,K)=U2(L,K)+UVPT(L,K)/HMU(L)
          V2(L,K)=V2(L,K)+VVPT(L,K)/HMV(L)
          W2(L,K)=W2(L,K)+WVPT(L,K)
         ENDDO
        ENDDO
!
      ENDIF
!
! **  ADDITIONAL 3D CONTINUITY ADJUSTED ADDED BELOW
!
      IF(KC.GT.1)THEN
!
      DO L=2,LA
        TVAR3E(L)=0.
        TVAR3N(L)=0.
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
        TVAR3E(L)=TVAR3E(L)+UHDY2(L,K)*DZC(K)
        TVAR3N(L)=TVAR3N(L)+VHDX2(L,K)*DZC(K)
       ENDDO
      ENDDO
!
      IF(ISTL.EQ.3)THEN
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H2P(L)+DELT*DXYIP(L)*(QSUME(L)-TVAR3E(L+1)+TVAR3E(L)-TVAR3N(LN) +TVAR3N(L))
        IF(ISGWIE.GE.1) HPPTMP=HPPTMP-DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        HPI(L)=1./HP(L)
        ENDDO
       ELSE
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H1P(L)+DELT*DXYIP(L)*(QSUME(L)-TVAR3E(L+1)+TVAR3E(L)-TVAR3N(LN) +TVAR3N(L))
        IF(ISGWIE.GE.1) HPPTMP=HPPTMP-DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        HPI(L)=1./HP(L)
        ENDDO
      ENDIF
!
      IF(MDCHH.GE.1)THEN
        RLAMN=QCHERR
        RLAMO=1.-RLAMN
        DO NMD=1,MDCHH
        LHOST=LMDCHH(NMD)
        LCHNU=LMDCHU(NMD)
        LCHNV=LMDCHV(NMD)
        IF(MDCHTYP(NMD).EQ.1)THEN
          TMPVAL=DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUN(NMD))
          HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)
          HP(LCHNU)=HP(LCHNU)-TMPVAL*DXYIP(LCHNU)
          HPI(LHOST)=1./HP(LHOST)
          HPI(LCHNU)=1./HP(LCHNU)
        ENDIF
        IF(MDCHTYP(NMD).EQ.2)THEN
          TMPVAL=DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVN(NMD))
          HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)
          HP(LCHNV)=HP(LCHNV)-TMPVAL*DXYIP(LCHNV)
          HPI(LHOST)=1./HP(LHOST)
          HPI(LCHNV)=1./HP(LCHNV)
        ENDIF
        ENDDO
      ENDIF
!
      ENDIF
!
!GVCDIAG
!      IF(N.LE.NTSPTC)THEN
!      DO L=2,LA
!	  WRITE(8,891)N,IL(L),JL(L),UHDYE(L),(U(L,K),K=1,KC)
!      ENDDO
!      DO L=2,LA
!	  WRITE(8,892)N,IL(L),JL(L),HP(L),(W(L,K),K=1,KC)
!      ENDDO
!      ENDIF
!  891 FORMAT('FINAL QX,U SIG ',3I5,10E14.6)
!  892 FORMAT('FINAL HP,W SIG ',3I5,10E14.6)
!GVCDIAG
!**********************************************************************C
!
! **  ACCUMULTATE MAX COURANT NUMBERS
!
!hnr      DO K=1,KC
!hnr      DO L=2,LA
!hnr       CFLUUUT=DELT*ABS(DXIU(L)*U(L,K))
!hnr       CFLUUU(L,K)=MAX(CFLUUUT,CFLUUU(L,K))
!hnr       CFLVVVT=DELT*ABS(DYIV(L)*V(L,K))
!hnr       CFLVVV(L,K)=MAX(CFLVVVT,CFLVVV(L,K))
!hnr       CFLWWWT=DELT*ABS(HPI(L)*DZIG(K)*W(L,K))
!hnr       CFLWWW(L,K)=MAX(CFLWWWT,CFLWWW(L,K))
!hnr       CFLCACT=DELT*ABS(CAC(L,K)*DXYIP(L)*HPI(L))
!hnr       CFLCAC(L,K)=MAX(CFLCACT,CFLCAC(L,K))
!hnr      ENDDO
!hnr      ENDDO
      DO K=1,KC                                               !hnr
      DO L=2,LA                                               !hnr
       CFLUUU(L,K)=DELT*ABS(DXIU(L)*U(L,K))                   !hnr
       CFLVVV(L,K)=DELT*ABS(DYIV(L)*V(L,K))                   !hnr
       CFLWWW(L,K)=DELT*ABS(HPI(L)*DZIG(K)*W(L,K))            !hnr
       CFLCAC(L,K)=DELT*ABS(CAC(L,K)*DXYIP(L)*HPI(L))         !hnr
      ENDDO                                                   !hnr
      ENDDO                                                   !hnr
!
!**********************************************************************C
!
! ** CALCULATE NONHYDROSTATIC PRESSURE
!
      IF(KC.GT.1.AND.ISPNHYDS.GE.1) CALL CALPNHS
!
!**********************************************************************C
!
! **  WRITE TO DIAGNOSTIC FILE CFL.OUT WITH DIAGNOSTICS OF MAXIMUM
! **  TIME STEP
! **  SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2
!
!----------------------------------------------------------------------C
!
      IF(ISCFL.GE.1.AND.ISTL.EQ.3)THEN
!
        OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN',POSITION='APPEND')
        IF(ISCFLM.GE.1.AND.N.EQ.1)THEN
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')
          CLOSE(2,STATUS='DELETE')
          DO L=1,LC
          ICFLMP(L)=0
          ENDDO
        ENDIF
!
        DTCFL=1.E+18
!
        K=1
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))
        UETMP=ABS(DXIU(L+1)*U2(L+1,K))
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))
        WBTMP=0.
        WTTMP=ABS(HPI(L)*DZIC(K)*W2(L,K))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)+1.0E-12
        DTMAX=0.5/DTMAXI
        IF(DTMAX.LT.DTCFL)THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
          KCFL=K
        ENDIF
        ENDDO
!
        IF(KC.GT.1)THEN
        K=KC
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))
        UETMP=ABS(DXIU(L+1)*U2(L+1,K))
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))
        WTTMP=0.
        WBTMP=ABS(HPI(L)*DZIC(K)*W2(L,K-1))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)+1.0E-12
        DTMAX=0.5/DTMAXI
        IF(DTMAX.LT.DTCFL)THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
          KCFL=K
        ENDIF
        ENDDO
        ENDIF
!
        IF(KC.GT.2)THEN
        DO K=2,KS
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))
        UETMP=ABS(DXIU(L+1)*U2(L+1,K))
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))
        WBTMP=ABS(HPI(L)*DZIC(K)*W2(L,K-1))
        WTTMP=ABS(HPI(L)*DZIC(K)*W2(L,K  ))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)+1.0E-12
        DTMAX=0.5/DTMAXI
        IF(DTMAX.LT.DTCFL)THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
          KCFL=K
        ENDIF
        ENDDO
        ENDDO
        ENDIF
!
        IVAL=MOD(N,ISCFL)
        IDTCFL=NINT(DTCFL)
        IF(ISCFL.EQ.1) WRITE(1,1212)DTCFL,N,ICFL,JCFL,KCFL
        IF(ISCFL.GE.2.AND.IVAL.EQ.0 )WRITE(1,1213)IDTCFL
        IF(ISCFLM.GE.1 )THEN
          LTMP=LIJ(ICFL,JCFL)
          ICFLMP(LTMP)=ICFLMP(LTMP)+1
        ENDIF
!
        IF(ISCFLM.GE.1.AND.N.EQ.NTS)THEN
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')
          TMPVALN=1./FLOAT(NTS)
          DO L=2,LA
           TMPVAL=TMPVALN*FLOAT(ICFLMP(L))
           WRITE(2,1214)IL(L),JL(L),ICFLMP(L),TMPVAL
          ENDDO
          CLOSE(2)
        ENDIF
!
        CLOSE(1)
      ENDIF
!
!----------------------------------------------------------------------C
!
      IF(ISCFL.GE.1.AND.IS1DCHAN.GE.1)THEN
!
        OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN',POSITION='APPEND')
        IF(ISCFLM.GE.1.AND.N.EQ.1)THEN
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')
          CLOSE(2,STATUS='DELETE')
          DO L=1,LC
          ICFLMP(L)=0
          ENDDO
        ENDIF
!
        DTCFL=1.E+18
!
        K=1
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))
        UETMP=ABS(DXIU(L+1)*U2(L+1,K))
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+1.0E-12
        DTMAX=0.5/DTMAXI
        IF(DTMAX.LT.DTCFL)THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
        ENDIF
        ENDDO
!
        IVAL=MOD(N,ISCFL)
        IDTCFL=NINT(DTCFL)
        IF(ISCFL.EQ.1) WRITE(1,1212)DTCFL,N,ICFL,JCFL,KCFL
        IF(ISCFL.GE.2.AND.IVAL.EQ.0 )WRITE(1,1213)IDTCFL
        IF(ISCFLM.GE.1 )THEN
          LTMP=LIJ(ICFL,JCFL)
          ICFLMP(LTMP)=ICFLMP(LTMP)+1
        ENDIF
!
        IF(ISCFLM.GE.1.AND.N.EQ.NTS)THEN
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')
          TMPVALN=1./FLOAT(NTS)
          DO L=2,LA
           TMPVAL=TMPVALN*FLOAT(ICFLMP(L))
           WRITE(2,1214)IL(L),JL(L),ICFLMP(L),TMPVAL
          ENDDO
          CLOSE(2)
        ENDIF
!
        CLOSE(1)
      ENDIF
!
 1212 FORMAT(' MAX TIME STEP =',F10.2,' SEC FOR N,I,J,K =',I8,3I5)
 1213 FORMAT(I4)
 1214 FORMAT(2I5,I12,F10.2)
!
!**********************************************************************C
!
      RETURN
      END
