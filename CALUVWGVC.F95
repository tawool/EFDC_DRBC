!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALUVWGVC (ISTL,IS2TL)
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
	  KBTMPU=KGVCU(L)
	  KBTMPV=KGVCV(L)
        RCX(L)=AVCON1/H1U(L)+STBX(L)*SQRT(U1(L,KBTMPU)*U1(L,KBTMPU)+V1U(L)*V1U(L))
        RCY(L)=AVCON1/H1V(L)+STBY(L)*SQRT(U1V(L)*U1V(L)+V1(L,KBTMPV)*V1(L,KBTMPV))
       ENDDO
!
!GVCDIAG
!      IF(N.LE.4)THEN
!      DO L=2,LA
!	  KBTMPU=KGVCU(L)
!	  WRITE(8,8881)N,IL(L),JL(L),KBTMPU,RCX(L),H1U(L),STBX(L),
!     &               U1(L,KBTMPU),V1U(L)
!      ENDDO
!      ENDIF
! 8881 FORMAT('KBU,RCX,H1U,ST,U,V GVC ',4I5,10E14.6)
!GVCDIAG
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
	  KBTMPU=KGVCU(L)
	  KBTMPV=KGVCV(L)
        Q1=SQRT(U1(L,KBTMPU)*U1(L,KBTMPU)+V1U(L)*V1U(L))
        Q2=SQRT(U(L,KBTMPU)*U(L,KBTMPU)+VU(L)*VU(L))
        RCX(L)=STBX(L)*SQRT(Q1*Q2)
        Q1=SQRT(U1V(L)*U1V(L)+V1(L,KBTMPV)*V1(L,KBTMPV))
        Q2=SQRT(UV(L)*UV(L)+V(L,KBTMPV)*V(L,KBTMPV))
        RCY(L)=STBY(L)*SQRT(Q1*Q2)
       ENDDO
	 ELSE
       DO L=2,LA
	  KBTMPU=KGVCU(L)
	  KBTMPV=KGVCV(L)
        Q1=SQRT(U1(L,KBTMPU)*U1(L,KBTMPU)+V1U(L)*V1U(L))
        Q2=SQRT(U(L,KBTMPU)*U(L,KBTMPU)+VU(L)*VU(L))
        RCX(L)=AVCON1/SQRT(H1U(L)*HU(L))+STBX(L)*SQRT(Q1*Q2)
        Q1=SQRT(U1V(L)*U1V(L)+V1(L,KBTMPV)*V1(L,KBTMPV))
        Q2=SQRT(UV(L)*UV(L)+V(L,KBTMPV)*V(L,KBTMPV))
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
! **  working definitions
!
!           CDZM(K)=0.5*DZC(K)*DZC(K+1)
!           CDZU(K)=-DZC(K)/(DZC(K)+DZC(K+1))
!           CDZL(K)=-DZC(K+1)/(DZC(K)+DZC(K+1))
!
!
!           CDZR(1)=DZC(1)-1.
!           CDZD(1)=DZC(1)
!           DO K=2,KS
!           CDZR(K)=DZC(K)+CDZR(K-1)
!           CDZD(K)=DZC(K)+CDZD(K-1)
!           ENDDO
!
!           DO K=1,KS
!           CDZR(K)=CDZR(K)*DZG(K)*CDZL(1)
!           ENDDO
!
!
! **  first interface at k=1, cases for null and kbu/kbv =1
!
!GVCDIAG
!      IF(N.LE.4)THEN
!      DO L=2,LA
!	  WRITE(8,881)N,IL(L),JL(L),RCX(L),UHE(L),(DU(L,K),K=1,KS)
!      ENDDO
!      ENDIF
!  881 FORMAT('R,UHE,DU GVC ',3I5,10E14.6)
!GVCDIAG

!
      RCDZM=CDZM(1)*DELTI
      RCDZU=CDZU(1)
      RCDZL=CDZL(1)
       DO L=2,LA
	   IF(LGVCU(L,1))THEN
           CMU=1.+RCDZM*GVCSCLU(L)*HU(L)*AVUI(L,1)
           EU=1./CMU
           CU1(L,1)=RCDZU*EU
           DU(L,1)=(DU(L,1)-RCDZL*RCX(L)*UHE(L)*HUI(L))*EU
           UUU(L,1)=EU
         ELSE
           CMU=1.
!           EU=0.
           CU1(L,1)=0.
           DU(L,1)=0.
           UUU(L,1)=0.
         ENDIF
       ENDDO
       DO L=2,LA
	   IF(LGVCV(L,1))THEN
           CMV=1.+RCDZM*GVCSCLV(L)*HV(L)*AVVI(L,1)
           EV=1./CMV
           CU2(L,1)=RCDZU*EV
           DV(L,1)=(DV(L,1)-RCDZL*RCY(L)*VHE(L)*HVI(L))*EV
           VVV(L,1)=EV
         ELSE
           CMV=1.
!           EV=0.
           CU2(L,1)=0.
           DV(L,1)=0.
           VVV(L,1)=0.
         ENDIF
       ENDDO
!
       DO K=2,KS
        RCDZM=CDZM(K)*DELTI
        RCDZU=CDZU(K)
        RCDZL=CDZL(K)
        DO L=2,LA
	    IF(LGVCU(L,K))THEN
	      IF(K.EQ.KGVCU(L))THEN
              CMU=1.+RCDZM*GVCSCLU(L)*HU(L)*AVUI(L,K)
              EU=1./(CMU-RCDZL*CU1(L,K-1))
              CU1(L,K)=RCDZU*EU
              DU(L,K)=(DU(L,K)-RCDZL*RCX(L)*UHE(L)*HUI(L))*EU
              UUU(L,K)=EU
		  ELSE
              CMU=1.+RCDZM*GVCSCLU(L)*HU(L)*AVUI(L,K)
              EU=1./(CMU-RCDZL*CU1(L,K-1))
              CU1(L,K)=RCDZU*EU
              DU(L,K)=(DU(L,K)-RCDZL*DU(L,K-1))*EU
              UUU(L,K)=-RCDZL*UUU(L,K-1)*EU
		  ENDIF
          ELSE
            CMU=1.
!            EU=0.
            CU1(L,K)=0.
            DU(L,K)=0.
            UUU(L,K)=0.
          ENDIF
        ENDDO
        DO L=2,LA
	    IF(LGVCV(L,K))THEN
	      IF(K.EQ.KGVCV(L))THEN
              CMV=1.+RCDZM*GVCSCLV(L)*HV(L)*AVVI(L,K)
              EV=1./(CMV-RCDZL*CU2(L,K-1))
              CU2(L,K)=RCDZU*EV
              DV(L,K)=(DV(L,K)-RCDZL*RCY(L)*VHE(L)*HVI(L))*EV
              VVV(L,K)=EV
		  ELSE
              CMV=1.+RCDZM*GVCSCLV(L)*HV(L)*AVVI(L,K)
              EV=1./(CMV-RCDZL*CU2(L,K-1))
              CU2(L,K)=RCDZU*EV
              DV(L,K)=(DV(L,K)-RCDZL*DV(L,K-1))*EV
              VVV(L,K)=-RCDZL*VVV(L,K-1)*EV
		  ENDIF
          ELSE
            CMV=1.
!            EV=0.
            CU2(L,K)=0.
            DV(L,K)=0.
            VVV(L,K)=0.
          ENDIF
        ENDDO
       ENDDO
!
! **  back substitution
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
!  884 FORMAT('BKSUB DU  GVC ',3I5,10E14.6)
!  885 FORMAT('BKSUB UUU GVC ',3I5,10E14.6)
!GVCDIAG
!
! **  set SHERMAN-MORRISON BACK SUBSTITUTION
!
       DO L=2,LA
        AAU(L)=0.
        AAV(L)=0.
        BBU(L)=1.
        BBV(L)=1.
       ENDDO
!
!           CDZR(1)=DZC(1)-1.
!
!           DO K=2,KS
!           CDZR(K)=DZC(K)+CDZR(K-1)
!           ENDDO
!
!           DO K=1,KS
!           CDZR(K)=CDZR(K)*DZG(K)*CDZL(1)
!           ENDDO
!
!           CDZL(K)=-DZC(K+1)/(DZC(K)+DZC(K+1))
!
       DO K=1,KS
        DO L=2,LA
	   IF(LGVCU(L,K))THEN
           CRU=CDZRGVCU(L,K)*RCX(L)*AVUI(L,K)
           AAU(L)=AAU(L)+CRU*DU(L,K)
           BBU(L)=BBU(L)+CRU*UUU(L,K)
	   ENDIF
        ENDDO
       ENDDO
!
       DO K=1,KS
        DO L=2,LA
	   IF(LGVCV(L,K))THEN
           CRV=CDZRGVCV(L,K)*RCY(L)*AVVI(L,K)
           AAV(L)=AAV(L)+CRV*DV(L,K)
           BBV(L)=BBV(L)+CRV*VVV(L,K)
	   ENDIF
        ENDDO
       ENDDO
!
       DO L=2,LA
        AAU(L)=AAU(L)/BBU(L)
        AAV(L)=AAV(L)/BBV(L)
       ENDDO
C
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
!  882 FORMAT('SMW DU GVC ',3I5,10E14.6)
!GVCDIAG
!
!     note that du and dv are actually the depth time the velocity
!     difference
!
!**********************************************************************C
!
! **  CALCULATED U AND V
!
! **  DUSUM+UHE=UHE, DVSUM+VHE=VHE
!
!----------------------------------------------------------------------C
!
!
!        CDZD(1)=DZC(1)
!         DO K=2,KS
!         CDZD(K)=DZC(K)+CDZD(K-1)
!         ENDDO
!
!	   CDZDGVCU(L,K)=GVCSCLU(L)*(CDZD(K) FIXED FOR GVC)
!	   CDZDGVCU(L,K)=GVCSCLV(L)*(CDZD(K) FIXED FOR GVC)
!
!     solve for the surface layer velocity * depth and temporarily
!     store results in uhe and vhe
!
       DO K=1,KS
!        RCDZD=CDZD(K)
        DO L=2,LA
!	   IF(LGVCU(L,K)) UHE(L)=UHE(L)+RCDZD*GVCSCLU(L)*DU(L,K)
	   IF(LGVCU(L,K)) UHE(L)=UHE(L)+CDZDGVCU(L,K)*DU(L,K)
        ENDDO
       ENDDO
!
       DO K=1,KS
!        RCDZD=CDZD(K)
        DO L=2,LA
!         IF(LGVCV(L,K)) VHE(L)=VHE(L)+RCDZD*GVCSCLV(L)*DV(L,K)
         IF(LGVCV(L,K)) VHE(L)=VHE(L)+CDZDGVCV(L,K)*DV(L,K)
        ENDDO
       ENDDO
!
!      move surface layer velocity times depth to uhdy and vhdx
!
       DO L=2,LA
        UHDY(L,KC)=UHE(L)*SUB3D(L,KC)
        VHDX(L,KC)=VHE(L)*SVB3D(L,KC)
       ENDDO
!
       DO K=KS,1,-1
        DO L=2,LA
         IF(LGVCU(L,K))THEN
	     UHDY(L,K)=SUB3D(L,K)*(UHDY(L,K+1)-DU(L,K))
         ELSE
	     UHDY(L,K)=0.0
         ENDIF
        ENDDO
       ENDDO
!
       DO K=KS,1,-1
        DO L=2,LA
         IF(LGVCV(L,K))THEN
	     VHDX(L,K)=SVB3D(L,K)*(VHDX(L,K+1)-DV(L,K))
         ELSE
	     VHDX(L,K)=0.0
         ENDIF
        ENDDO
       ENDDO
!
      DO K=1,KC
       DO L=2,LA
        U(L,K)=SUB3D(L,K)*UHDY(L,K)*HUI(L)
        V(L,K)=SVB3D(L,K)*VHDX(L,K)*HVI(L)
        UHDY(L,K)=SUB3D(L,K)*UHDY(L,K)*DYU(L)
        VHDX(L,K)=SVB3D(L,K)*VHDX(L,K)*DXV(L)
       ENDDO
      ENDDO
!
!GVCDIAG
!      IF(N.LE.4)THEN
!      DO L=2,LA
!	  WRITE(8,883)N,IL(L),JL(L),(UHDY(L,K),K=1,KC)
!      ENDDO
!      ENDIF
!  883 FORMAT('UHDY GVC ',3I5,10E14.6)
!GVCDIAG
!
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
        TVAR3E(L)=TVAR3E(L)+GVCSCLU(L)*UHDY(L,K)*DZC(K)
        TVAR3N(L)=TVAR3N(L)+GVCSCLV(L)*VHDX(L,K)*DZC(K)
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
        IF(LGVCU(L,K)) UHDY(L,K)=UHDY(L,K)-GVCSCLUI(L)*TVAR3E(L)*DZIC(K)
       ENDDO
      ENDDO
!
      DO K=1,KC
       DO L=2,LA
        IF(LGVCV(L,K)) VHDX(L,K)=VHDX(L,K)-GVCSCLVI(L)*TVAR3N(L)*DZIC(K)
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
        UHE(L)=UHE(L)+GVCSCLU(L)*UHDY(L,K)*DZC(K)
        VHE(L)=VHE(L)+GVCSCLV(L)*VHDX(L,K)*DZC(K)
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
!      DO L=2,LA
!        TVAR3E(L)=0.
!        TVAR3N(L)=0.
!      ENDDO
!
!      DO K=1,KC
!       DO L=2,LA
!        TVAR3E(L)=TVAR3E(L)+GVCSCLU(L)*UHDY(L,K)*DZC(K)
!        TVAR3N(L)=TVAR3N(L)+GVCSCLV(L)*VHDX(L,K)*DZC(K)
!       ENDDO
!      ENDDO
!
!      UERMX=-1.E+12
!      UERMN=1.E+12
!      VERMX=-1.E+12
!      VERMN=1.E+12
!      DO L=2,LA
!        TVAR3E(L)=TVAR3E(L)-UHDYE(L)
!        TVAR3N(L)=TVAR3N(L)-VHDXE(L)
!        IF(TVAR3E(L).GT.UERMX)THEN
!          UERMX=TVAR3E(L)
!          LUMX=L
!        ENDIF
!        IF(TVAR3E(L).LT.UERMN)THEN
!          UERMN=TVAR3E(L)
!          LUMN=L
!        ENDIF
!        IF(TVAR3N(L).GT.VERMX)THEN
!          VERMX=TVAR3N(L)
!          LVMX=L
!        ENDIF
!        IF(TVAR3N(L).LT.VERMN)THEN
!          VERMN=TVAR3N(L)
!          LVMN=L
!        ENDIF
!      ENDDO
!
!      WRITE(8,6661)IL(LUMX),JL(LUMX),UERMX
!      WRITE(8,6662)IL(LUMN),JL(LUMN),UERMN
!      WRITE(8,6663)IL(LVMX),JL(LVMX),VERMX
!      WRITE(8,6664)IL(LVMN),JL(LVMN),VERMN
!
 6660 FORMAT(' INT MOD ERR AT = ',I8)
 6661 FORMAT(' I,J,UHDYERMX = ',2I5,10E12.4)
 6662 FORMAT(' I,J,UHDYERMN = ',2I5,10E12.4)
 6663 FORMAT(' I,J,VHDXERMX = ',2I5,10E12.4)
 6664 FORMAT(' I,J,VHDXERMX = ',2I5,10E12.4)
 6665 FORMAT(' I,J,UMXSTUFF = ',2I5,10E12.4)
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
	  LN=LNC(L)
	  LE=L+1
        TVAR3E(L)=GVCSCLP(L)*UHDYE(LE)
        TVAR3N(L)=GVCSCLP(L)*VHDXE(LN)
        TVAR3W(L)=GVCSCLP(L)*UHDY2E(LE)
        TVAR3S(L)=GVCSCLP(L)*VHDX2E(LN)
       ENDDO
!
!      mode here

       DO K=1,KC
        DO L=2,LA
	   LN=LNC(L)
	   LE=L+1
         TVAR1E(L,K)=GVCSCLU(LE)*UHDY(LE,K)
         TVAR1N(L,K)=GVCSCLV(LN)*VHDX(LN,K)
         TVAR1W(L,K)=GVCSCLU(LE)*UHDY2(LE,K)
         TVAR1S(L,K)=GVCSCLV(LN)*VHDX2(LN,K)
        ENDDO
       ENDDO
!
!020906       DO K=1,KS
!        DO L=2,LA
!	   IF(K.GE.KGVCP(L))THEN
!         LN=LNC(L)
!         W(L,K)=SWB3D(L,K)*( W(L,K-1)-W2(L,K)+W2(L,K-1)
!     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
!     &          -TVAR3E(L)+GVCSCLU(L)*UHDYE(L)
!     &          +TVAR1W(L,K)-GVCSCLU(L)*UHDY2(L,K)
!     &          -TVAR3W(L)+GVCSCLU(L)*UHDY2E(L)
!     &          +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
!     &          -TVAR3N(L)+GVCSCLV(L)*VHDXE(L)
!     &          +TVAR1S(L,K)-GVCSCLV(L)*VHDX2(L,K)
!     &          -TVAR3S(L)+GVCSCLV(L)*VHDX2E(L)
!     &          )*DXYIP(L) )
!     &          +2.*SWB3D(L,K)*( QSUM(L,K)
!     &          -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
!        ENDIF
!        ENDDO
!       ENDDO
!
       DO K=1,KS
        DO L=2,LA
	   IF(K.GE.KGVCP(L))THEN
         LN=LNC(L)
         W(L,K)=SWB3D(L,K)*( W(L,K-1)-W2(L,K)+W2(L,K-1)                                          &
       -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)                                                 &
               -TVAR3E(L)+GVCSCLP(L)*UHDYE(L)                                                    &
               +TVAR1W(L,K)-GVCSCLU(L)*UHDY2(L,K)                                                &
               -TVAR3W(L)+GVCSCLP(L)*UHDY2E(L)                                                   &
               +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)                                                 &
               -TVAR3N(L)+GVCSCLP(L)*VHDXE(L)                                                    &
               +TVAR1S(L,K)-GVCSCLV(L)*VHDX2(L,K)                                                &
               -TVAR3S(L)+GVCSCLP(L)*VHDX2E(L)                                                   &
               )*DXYIP(L) )                                                                      &
               +2.*SWB3D(L,K)*( QSUM(L,K)                                                        &
               -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
        ENDIF
        ENDDO
       ENDDO
!
! **  UNCOMMENT BELOW FOR CONTINUITY DIAGNOSTICS
!
!      WSFMAX=-1.E+12
!      WSFMIN=1.E+12
!      SURFOUT=0.
!      K=KC
!      DO L=2,LA
!      IF(SWB(L).GT.0.5)THEN
!         WSURF=SWB(L)*( W(L,K-1)-W2(L,K)+W2(L,K-1)
!     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
!     &   -TVAR3E(L)+GVCSCLP(L)*UHDYE(L)
!     &       +TVAR1W(L,K)-GVCSCLU(L)*UHDY2(L,K)
!     &   -TVAR3W(L)+GVCSCLP(L)*UHDY2E(L)
!     &       +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
!     &   -TVAR3N(L)+GVCSCLP(L)*VHDXE(L)
!     &       +TVAR1S(L,K)-GVCSCLV(L)*VHDX2(L,K)
!     &   -TVAR3S(L)+GVCSCLP(L)*VHDX2E(L) )*DXYIP(L) )
!     &       +2.*SWB(L)*( QSUM(L,K)-
!     &   DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
!      SURFOUT=SURFOUT+WSURF*DXYP(L)
!      IF(WSURF.GT.WSFMAX)THEN
!        WSFMAX=WSURF
!        IMAX=IL(L)
!        JMAX=JL(L)
!      ENDIF
!      IF(WSURF.LT.WSFMIN)THEN
!        WSFMIN=WSURF
!        IMIN=IL(L)
!        JMIN=JL(L)
!      ENDIF
!      ENDIF
!      ENDDO
!
!      L=LIJ(IMAX,JMAX)
!      WSFMAX=WSFMAX*DXYP(L)
!      L=LIJ(IMIN,JMIN)
!      WSFMIN=WSFMIN*DXYP(L)
!      WRITE(8,601)IMAX,JMAX,WSFMAX
!      WRITE(8,602)IMIN,JMIN,WSFMIN
!      WRITE(8,603)SURFOUT
!
      ENDIF
!
      IF(ISTL.EQ.2)THEN
!
       DO L=2,LA
	  LN=LNC(L)
	  LE=L+1
        TVAR3E(L)=GVCSCLP(L)*UHDYE(LE)
        TVAR3N(L)=GVCSCLP(L)*VHDXE(LN)
        TVAR3W(L)=GVCSCLP(L)*UHDY1E(LE)
        TVAR3S(L)=GVCSCLP(L)*VHDX1E(LN)
       ENDDO
!
       DO K=1,KC
        DO L=2,LA
	   LN=LNC(L)
	   LE=L+1
         TVAR1E(L,K)=GVCSCLU(LE)*UHDY(LE,K)
         TVAR1N(L,K)=GVCSCLV(LN)*VHDX(LN,K)
         TVAR1W(L,K)=GVCSCLU(LE)*UHDY1(LE,K)
         TVAR1S(L,K)=GVCSCLV(LN)*VHDX1(LN,K)
        ENDDO
       ENDDO
!
!020906       DO K=1,KS
!        DO L=2,LA
!	   IF(K.GE.KGVCP(L))THEN
!         LN=LNC(L)
!         W(L,K)=SWB3D(L,K)*( W(L,K-1)-W1(L,K)+W1(L,K-1)
!     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
!     &          -TVAR3E(L)+GVCSCLU(L)*UHDYE(L)
!     &          +TVAR1W(L,K)-GVCSCLU(L)*UHDY1(L,K)
!     &          -TVAR3W(L)+GVCSCLU(L)*UHDY1E(L)
!     &          +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
!     &          -TVAR3N(L)+GVCSCLV(L)*VHDXE(L)
!     &          +TVAR1S(L,K)-GVCSCLV(L)*VHDX1(L,K)
!     &          -TVAR3S(L)+GVCSCLV(L)*VHDX1E(L)
!     &          )*DXYIP(L) )
!     &          +2.*SWB3D(L,K)*( QSUM(L,K)
!     &          -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
!         ENDIF
!        ENDDO
!       ENDDO
!
       DO K=1,KS
        DO L=2,LA
	   IF(K.GE.KGVCP(L))THEN
         LN=LNC(L)
         W(L,K)=SWB3D(L,K)*( W(L,K-1)-W1(L,K)+W1(L,K-1)-DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K) &
                -TVAR3E(L)+GVCSCLP(L)*UHDYE(L)+TVAR1W(L,K)-GVCSCLU(L)*UHDY1(L,K)-TVAR3W(L)       &
                +GVCSCLP(L)*UHDY1E(L)+TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)-TVAR3N(L)+GVCSCLP(L)      &
                *VHDXE(L)+TVAR1S(L,K)-GVCSCLV(L)*VHDX1(L,K)-TVAR3S(L)+GVCSCLP(L)*VHDX1E(L))      &
                *DXYIP(L))+2.*SWB3D(L,K)*(QSUM(L,K)-DZC(K)*GVCSCLP(L)*QSUME(L))*DXYIP(L)
         ENDIF
        ENDDO
       ENDDO
!
!  UNCOMMENT BELOW FOR SURFACE LAYER CONTINUITY DIAGNOSTIC
!
!      WSFMAX=-1.E+12
!      WSFMIN=1.E+12
!      SURFOUT=0.
!      K=KC
!      DO L=2,LA
!      IF(SWB(L).GT.0.5)THEN
!         WSURF=SWB(L)*( W(L,K-1)-W1(L,K)+W1(L,K-1)
!     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
!     &   -TVAR3E(L)+GVCSCLP(L)*UHDYE(L)
!     &       +TVAR1W(L,K)-GVCSCLU(L)*UHDY1(L,K)
!     &   -TVAR3W(L)+GVCSCLP(L)*UHDY1E(L)
!     &       +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
!     &   -TVAR3N(L)+GVCSCLP(L)*VHDXE(L)
!     &       +TVAR1S(L,K)-GVCSCLV(L)*VHDX1(L,K)
!     &   -TVAR3S(L)+GVCSCLP(L)*VHDX1E(L) )*DXYIP(L) )
!     &       +2.*SWB(L)*( QSUM(L,K)
!     &   -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
!      SURFOUT=SURFOUT+WSURF*DXYP(L)
!      IF(WSURF.GT.WSFMAX)THEN
!        WSFMAX=WSURF
!        IMAX=IL(L)
!        JMAX=JL(L)
!      ENDIF
!      IF(WSURF.LT.WSFMIN)THEN
!        WSFMIN=WSURF
!        IMIN=IL(L)
!        JMIN=JL(L)
!      ENDIF
!      ENDIF
!      ENDDO
!
!      L=LIJ(IMAX,JMAX)
!      WSFMAX=WSFMAX*DXYP(L)
!      L=LIJ(IMIN,JMIN)
!      WSFMIN=WSFMIN*DXYP(L)
!      WRITE(8,601)IMAX,JMAX,WSFMAX
!      WRITE(8,602)IMIN,JMIN,WSFMIN
!      WRITE(8,603)SURFOUT
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
        VHDX(LN,K)=SVB3D(LN,K)*(VHDX(LNN,K)-VHDXE(LNN)+VHDXE(LN))
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
        UHDY(LP,K)=SUB3D(LP,K)*(UHDY(LPP,K)-UHDYE(LPP)+UHDYE(LP))
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
      UHDY(L,K)=SUB3D(L,K)*(UHDY(L-1,K)-UHDYE(L-1)+UHDYE(L))
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
      VHDX(L,K)=SVB3D(L,K)*(VHDX(LS,K)-VHDXE(LS)+VHDXE(L))
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
          V2(L,K)=U2(L,K)+VVPT(L,K)/HMV(L)
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
        TVAR3E(L)=TVAR3E(L)+GVCSCLU(L)*UHDY2(L,K)*DZC(K)
        TVAR3N(L)=TVAR3N(L)+GVCSCLV(L)*VHDX2(L,K)*DZC(K)
       ENDDO
      ENDDO
!
      IF(ISTL.EQ.3)THEN
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H2P(L)+DELT*DXYIP(L)*(QSUME(L)-TVAR3E(L+1)+TVAR3E(L)-TVAR3N(LN) +TVAR3N(L))
        IF(ISGWIE.GE.1) HPPTMP=HPPTMP-DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
	  ERRHPP=ABS(HP(L)-HPPTMP)
!	  IF(ERRHPP.GT.0.00001.AND.SPB(L).GT.0.5)
!     &    WRITE(8,6999)N,IL(L),JL(L),HPPTMP,HP(L)
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        HPI(L)=1./HP(L)
        ENDDO
       ELSE
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H1P(L)+DELT*DXYIP(L)*(QSUME(L)-TVAR3E(L+1)+TVAR3E(L)-TVAR3N(LN) +TVAR3N(L))
        IF(ISGWIE.GE.1) HPPTMP=HPPTMP-DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
	  ERRHPP=ABS(HP(L)-HPPTMP)
!	  IF(ERRHPP.GT.0.00001.AND.SPB(L).GT.0.5)
!     &    WRITE(8,6999)N,IL(L),JL(L),HPPTMP,HP(L)
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
!  891 FORMAT('FINAL QX,U GVC ',3I5,10E14.6)
!  892 FORMAT('FINAL HP,W GVC ',3I5,10E14.6)
!GVCDIAG
!
 6999 FORMAT('  DEP ADJ ERR ',3I5,2E14.6)
!
!**********************************************************************C
!
! **  ACCUMULTATE MAX COURANT NUMBERS
!
!hnr  DO K=1,KC
!hnr  DO L=2,LA
!hnr       CFLUUUT=DELT*ABS(DXIU(L)*U(L,K))
!hnr       CFLUUU(L,K)=MAX(CFLUUUT,CFLUUU(L,K))
!hnr       CFLVVVT=DELT*ABS(DYIV(L)*V(L,K))
!hnr       CFLVVV(L,K)=MAX(CFLVVVT,CFLVVV(L,K))
!hnr      CFLWWWT=DELT*ABS(HPI(L)*DZIG(K)*W(L,K))
!hnr       CFLWWW(L,K)=MAX(CFLWWWT,CFLWWW(L,K))
!hnr       CFLCACT=DELT*ABS(CAC(L,K)*DXYIP(L)*HPI(L))
!hnr       CFLCAC(L,K)=MAX(CFLCACT,CFLCAC(L,K))
!hnr      ENDDO
!hnr      ENDDO
      DO K=1,KC
      DO L=2,LA
       CFLUUU(L,K)=DELT*ABS(DXIU(L)*U(L,K))
       CFLVVV(L,K)=DELT*ABS(DYIV(L)*V(L,K))
       CFLWWW(L,K)=DELT*ABS(HPI(L)*DZIG(K)*W(L,K))
       CFLCAC(L,K)=DELT*ABS(CAC(L,K)*DXYIP(L)*HPI(L))
      ENDDO
      ENDDO
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