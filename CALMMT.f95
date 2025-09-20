!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALMMT
!------------------------------------------------------------------------------
! PURPOSE:
!
!   Subroutine CALMMTF calculates the mean mass transport field
!
! VARIABLE LIST:
!
! MODIFICATION HISTORY:
!
!   Date       Author             Comments
!   ---------- ------------------ ---------------------------------------------
!   02/01/1990 John M. Hamrick    Orignial author
!   02/12/2001 Mike Morton        Cleaned code, reformatted code
!   06/27/2002 John M. Hamrick    Make changes for ICM and WASP interfaces
!   04/17/2006 Hugo N Rodriguez   Added ulpf and vlpf to normal wasp output to be used with hybrid
!   06/24/2019 Hugo N Rodriguez   ADJUSTED VALUES FOR DRY CELLS WHEN W/D IS USED
!   12/20/2020 DRBC LZ            AVERAGED TURBULENCE DISSIPATION RATE "DISPRATE"
!
!------------------------------------------------------------------------------
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
! **  INITIALIZE CE-QUAL-ICM INTERFACE
!
      IF(ISICM.GE.1.AND.JSWASP.EQ.1) CALL CEQICM
!
!**********************************************************************C
!
      IF (NMMT.GT.1) GO TO 100
!
!----------------------------------------------------------------------C
!
!     INITIALIZE LOW PASS FILTERED VARIABLES AND DISPLACEMENTS
!
!----------------------------------------------------------------------C
!
      IF(NTSMMT.LT.NTSPTC) THEN
        DO L=1,LC
          HLPF(L)=0.
          QSUMELPF(L)=0.
          UELPF(L)=0.
          VELPF(L)=0.
          RAINLPF(L)=0.
          EVPSLPF(L)=0.
          EVPGLPF(L)=0.
          RINFLPF(L)=0.
          GWLPF(L)=0.
        END DO

        DO NSC=1,NSED
          DO K=1,KB
            DO L=1,LC
              SEDBLPF(L,K,NSC)=0.
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KB
            DO L=1,LC
              SNDBLPF(L,K,NSN)=0.
            END DO
          END DO
        END DO
        DO NT=1,NTOX
          DO K=1,KB
            DO L=1,LC
              TOXBLPF(L,K,NT)=0.
            END DO
          END DO
        END DO

        HLPF(1)=HMIN
        HLPF(LC)=HMIN

        DO K=1,KS
          DO L=1,LC
            ABLPF(L,K)=0.
            ABEFF(L,K)=0.
            WLPF(L,K)=0.
          END DO
        END DO

        DO K=1,KC                    !DRBC, LZ, 12/20/2020
          DO L=1,LC
            DISPRLPF(L,K)=0.
          END DO
        END DO

        DO K=1,KC
          DO L=1,LC
            AHULPF(L,K)=0.
            AHVLPF(L,K)=0.
            SALLPF(L,K)=0.
            TEMLPF(L,K)=0.
            SFLLPF(L,K)=0.
            DYELPF(L,K)=0.
            UHLPF(L,K)=0.
            VHLPF(L,K)=0.
            ULPF(L,K)=0.  
            VLPF(L,K)=0.  
            QSUMLPF(L,K)=0.
          END DO
        END DO

        DO NSC=1,NSED
          DO K=1,KC
            DO L=1,LC
              SEDLPF(L,K,NSC)=0.
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KC
            DO L=1,LC
              SNDLPF(L,K,NSN)=0.
            END DO
          END DO
        END DO
        DO NT=1,NTOX
          DO K=1,KC
            DO L=1,LC
              TOXLPF(L,K,NT)=0.
            END DO
          END DO
        END DO

        DO NT=1,NTOX
          DO NS=1,NSED+NSND
            DO K=1,KC
              DO L=1,LC
                TXPFLPF(L,K,NS,NT)=0.
              END DO
            END DO
          END DO
        END DO

        DO NS=1,NQSER
          DO K=1,KC
            QSRTLPP(K,NS)=0.
            QSRTLPN(K,NS)=0.
          END DO
        END DO

        DO NS=1,NQCTL
          DO K=1,KC
            QCTLTLP(K,NS)=0.
          END DO
        END DO

        DO NMD=1,MDCHH
          QCHNULP(NMD)=0.
          QCHNVLP(NMD)=0.
        END DO

        DO NWR=1,NQWR
          QWRSERTLP(NWR)=0.
        END DO

      ELSE

        DO L=1,LC
          HLPF(L)=0.
          QSUMELPF(L)=0.
          UELPF(L)=0.
          VELPF(L)=0.
          RAINLPF(L)=0.
          EVPSLPF(L)=0.
          EVPGLPF(L)=0.
          RINFLPF(L)=0.
          GWLPF(L)=0.
        END DO

        DO NSC=1,NSED
          DO K=1,KB
            DO L=1,LC
              SEDBLPF(L,K,NSC)=0.
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KB
            DO L=1,LC
              SNDBLPF(L,K,NSN)=0.
            END DO
          END DO
        END DO
        DO NT=1,NTOX
          DO K=1,KB
            DO L=1,LC
              TOXBLPF(L,K,NT)=0.
            END DO
          END DO
        END DO

        HLPF(1)=HMIN
        HLPF(LC)=HMIN

        DO K=1,KS
          DO L=1,LC
            ABLPF(L,K)=0.
            WIRT(L,K)=0.
            WLPF(L,K)=0.
            WTLPF(L,K)=0.
          END DO
        END DO

        DO K=1,KC                    !DRBC, LZ, 12/20/2020
          DO L=1,LC
            DISPRLPF(L,K)=0.
          END DO
        END DO

        DO K=1,KC
          DO L=1,LC
            AHULPF(L,K)=0.
            AHVLPF(L,K)=0.
            SALLPF(L,K)=0.
            TEMLPF(L,K)=0.
            SFLLPF(L,K)=0.
            DYELPF(L,K)=0.
            UHLPF(L,K)=0.
            UIRT(L,K)=0.
            ULPF(L,K)=0.
            UTLPF(L,K)=0.
            VHLPF(L,K)=0.
            QSUMLPF(L,K)=0.
            VIRT(L,K)=0.
            VLPF(L,K)=0.
            VTLPF(L,K)=0.
          END DO
        END DO

        DO NSC=1,NSED
          DO K=1,KC
            DO L=1,LC
              SEDLPF(L,K,NSC)=0.
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KC
            DO L=1,LC
              SNDLPF(L,K,NSN)=0.
            END DO
          END DO
        END DO
        DO NT=1,NTOX
          DO K=1,KC
            DO L=1,LC
              TOXLPF(L,K,NT)=0.
            END DO
          END DO
        END DO

        DO NT=1,NTOX
          DO NS=1,NSED+NSND
            DO K=1,KC
              DO L=1,LC
                TXPFLPF(L,K,NS,NT)=0.
              END DO
            END DO
          END DO
        END DO

        DO NS=1,NQSER
          DO K=1,KC
            QSRTLPP(K,NS)=0.
            QSRTLPN(K,NS)=0.
          END DO
        END DO

        DO NS=1,NQCTL
          DO K=1,KC
            QCTLTLP(K,NS)=0.
          END DO
        END DO

        DO NMD=1,MDCHH
          QCHNULP(NMD)=0.
          QCHNVLP(NMD)=0.
        END DO

        DO NWR=1,NQWR
          QWRSERTLP(NWR)=0.
        END DO

      END IF

!**********************************************************************C
!
! **  ACCUMULATE FILTERED VARIABLES AND DISPLACEMENTS
!
!----------------------------------------------------------------------C
!
  100 CONTINUE
!
      IF(NTSMMT.LT.NTSPTC) THEN

        DO L=2,LA
          LN=LNC(L)
          
          IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
              HLPF(L)=HLPF(L)+HP(L)
          ELSE
              HLPF(L)=HLPF(L)+HDRY
          END IF
          
          QSUMELPF(L)=QSUMELPF(L)+QSUME(L)
          UTMP1=0.5*(UHDYE(L+1)+UHDYE(L))/(DYP(L)*HP(L))
          VTMP1=0.5*(VHDXE(LN)+VHDXE(L))/(DXP(L)*HP(L))
          UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
          VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
          UELPF(L)=UELPF(L)+UTMP
          VELPF(L)=VELPF(L)+VTMP
          RAINLPF(L)=RAINLPF(L)+DXYP(L)*RAINT(L)
        END DO

        IF (ISGWIE.EQ.0) THEN
          DO L=2,LA
            EVPSLPF(L)=EVPSLPF(L)+DXYP(L)*EVAPT(L)
            EVPGLPF(L)=0.
            RINFLPF(L)=0.
            GWLPF(L)=0.
          END DO
        ELSE
          DO L=2,LA
            EVPSLPF(L)=EVPSLPF(L)+EVAPSW(L)
            EVPGLPF(L)=EVPGLPF(L)+EVAPGW(L)
            RINFLPF(L)=RINFLPF(L)+RIFTR(L)
            GWLPF(L)=GWLPF(L)+AGWELV(L)
          END DO
        END IF

        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
              TOXBLPF(L,K,NT)=TOXBLPF(L,K,NT)+TOXB(L,K,NT)
            END DO
          END DO
        END DO
        DO NSC=1,NSED
          DO K=1,KB
            DO L=2,LA
              SEDBLPF(L,K,NSC)=SEDBLPF(L,K,NSC)+SEDB(L,K,NSC)
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KB
            DO L=2,LA
              SNDBLPF(L,K,NSN)=SNDBLPF(L,K,NSN)+SNDB(L,K,NSN)
            END DO
          END DO
        END DO

        IF(ISWASP.EQ.99.OR.ISICM.GE.1) THEN
          DO K=1,KS
            DO L=2,LA
              ABLPF(L,K)=ABLPF(L,K)+(AB(L,K)*HP(L))
              ABEFF(L,K)=ABEFF(L,K)+AB(L,K)*(SAL(L,K+1)-SAL(L,K))
              WLPF(L,K)=WLPF(L,K)+W(L,K)
            END DO
          END DO
        ELSE
          DO K=1,KS
            DO L=2,LA
              ABEFF(L,K)=ABEFF(L,K)+AB(L,K)*(SAL(L,K+1)-SAL(L,K))
              
              IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                  WLPF(L,K)=WLPF(L,K)+W(L,K)
                  ABLPF(L,K)=ABLPF(L,K)+AB(L,K)                   
              ELSE
                  WLPF(L,K)=WLPF(L,K)+0.0
                  ABLPF(L,K)=ABLPF(L,K)+0.0                   
              END IF
              
            END DO
          END DO
        END IF

        DO K=1,KC                                     ! DRBC, LZ, 12/20/2020
          DO L=2,LA 
            IF(ISCDRY(L).EQ.0) THEN 
              DISPRLPF(L,K)=DISPRLPF(L,K)+DISPRATE(L,K)
            ELSE
              DISPRLPF(L,K)=DISPRLPF(L,K)+0.0
            END IF
          END DO
        END DO

        DO K=1,KC
          DO L=2,LA
            LS=LSC(L)
            SFLLPF(L,K)=SFLLPF(L,K)+SFL(L,K)
            DYELPF(L,K)=DYELPF(L,K)+DYE(L,K)
            SALLPF(L,K)=SALLPF(L,K)+SAL(L,K)

            IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
              TEMLPF(L,K)=TEMLPF(L,K)+TEM(L,K)
            ELSE
              TEMLPF(L,K)=TEMLPF(L,K)+TATMT(L)
            END IF
            IF((ISCDRY(L).EQ.1).OR.(ISCDRY(L-1).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
              AHULPF(L,K)=AHULPF(L,K)+0.0
              UHLPF(L,K)=UHLPF(L,K)+0.0
            ELSE
              AHULPF(L,K)=AHULPF(L,K)+0.5*(AH(L,K)+AH(L-1,K))
              UHLPF(L,K)=UHLPF(L,K)+UHDY(L,K)/DYU(L)
            END IF
            IF((ISCDRY(L).EQ.1).OR.(ISCDRY(LS).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
              AHVLPF(L,K)=AHVLPF(L,K)+0.0
              VHLPF(L,K)=VHLPF(L,K)+0.0
            ELSE
              AHVLPF(L,K)=AHVLPF(L,K)+0.5*(AH(L,K)+AH(LS,K))
              VHLPF(L,K)=VHLPF(L,K)+VHDX(L,K)/DXV(L)
            END IF
            
            
            ULPF(L,K)=ULPF(L,K)+U(L,K)       
            VLPF(L,K)=VLPF(L,K)+V(L,K)       
            QSUMLPF(L,K)=QSUMLPF(L,K)+QSUM(L,K)
          END DO
        END DO

        DO NT=1,NTOX
          DO K=1,KC
            DO L=2,LA
              TOXLPF(L,K,NT)=TOXLPF(L,K,NT)+TOX(L,K,NT)
            END DO
          END DO
        END DO
        DO NSC=1,NSED
          DO K=1,KC
            DO L=2,LA
              SEDLPF(L,K,NSC)=SEDLPF(L,K,NSC)+SED(L,K,NSC)
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KC
            DO L=2,LA
              SNDLPF(L,K,NSN)=SNDLPF(L,K,NSN)+SND(L,K,NSN)
            END DO
          END DO
        END DO

        DO NT=1,NTOX
          DO NS=1,NSED+NSND
            DO K=1,KC
              DO L=1,LC
                TXPFLPF(L,K,NS,NT)=TXPFLPF(L,K,NS,NT)+TOXPFW(L,K,NS,NT)
              END DO
            END DO
          END DO
        END DO

        DO NS=1,NQSER
          DO K=1,KC
            QSRTLPP(K,NS)=QSRTLPP(K,NS)+MAX(QSERT(K,NS),0.)
            QSRTLPN(K,NS)=QSRTLPN(K,NS)+MIN(QSERT(K,NS),0.)
          END DO
        END DO

        DO NS=1,NQCTL
          DO K=1,KC
            QCTLTLP(K,NS)=QCTLTLP(K,NS)+QCTLT(K,NS)
          END DO
        END DO

        DO NMD=1,MDCHH
          QCHNULP(NMD)=QCHNULP(NMD)+QCHANU(NMD)
          QCHNVLP(NMD)=QCHNVLP(NMD)+QCHANV(NMD)
        END DO

        DO NWR=1,NQWR
          QWRSERTLP(NWR)=QWRSERTLP(NWR)+QWRSERT(NWR)
        END DO

      ELSE

        DO L=2,LA
          LN=LNC(L)
          
          IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
              HLPF(L)=HLPF(L)+HP(L)
          ELSE
              HLPF(L)=HLPF(L)+HDRY
          END IF
          
          QSUMELPF(L)=QSUMELPF(L)+QSUME(L)
          UTMP1=0.5*(UHDYE(L+1)+UHDYE(L))/(DYP(L)*HP(L))
          VTMP1=0.5*(VHDXE(LN)+VHDXE(L))/(DXP(L)*HP(L))
          UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
          VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
          UELPF(L)=UELPF(L)+UTMP
          VELPF(L)=VELPF(L)+VTMP
          RAINLPF(L)=RAINLPF(L)+DXYP(L)*RAINT(L)
        END DO

        IF (ISGWIE.EQ.0) THEN
          DO L=2,LA
            EVPSLPF(L)=EVPSLPF(L)+DXYP(L)*EVAPT(L)
            EVPGLPF(L)=0.
            RINFLPF(L)=0.
            GWLPF(L)=0.
          END DO
        ELSE
          DO L=2,LA
            EVPSLPF(L)=EVPSLPF(L)+EVAPSW(L)
            EVPGLPF(L)=EVPGLPF(L)+EVAPGW(L)
            RINFLPF(L)=RINFLPF(L)+RIFTR(L)
            GWLPF(L)=GWLPF(L)+AGWELV(L)
          END DO
        END IF

        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
              TOXBLPF(L,K,NT)=TOXBLPF(L,K,NT)+TOXB(L,K,NT)
            END DO
          END DO
        END DO
        DO NSC=1,NSED
          DO K=1,KB
            DO L=2,LA
              SEDBLPF(L,K,NSC)=SEDBLPF(L,K,NSC)+SEDB(L,K,NSC)
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KB
            DO L=2,LA
              SNDBLPF(L,K,NSN)=SNDBLPF(L,K,NSN)+SNDB(L,K,NSN)
            END DO
          END DO
        END DO

        IF(ISWASP.EQ.99.OR.ISICM.GE.1) THEN
          DO K=1,KS
            DO L=2,LA
              ABLPF(L,K)=ABLPF(L,K)+(AB(L,K)*HP(L))
              ABEFF(L,K)=ABEFF(L,K)+AB(L,K)*(SAL(L,K+1)-SAL(L,K))
              WIRT(L,K)=WIRT(L,K)+DT*W(L,K)
              WLPF(L,K)=WLPF(L,K)+W(L,K)
              WTLPF(L,K)=WTLPF(L,K)+DT*(FLOAT(NMMT)-0.5)*W(L,K)
            END DO
          END DO
        ELSE
          DO K=1,KS
            DO L=2,LA
              
              IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                  WLPF(L,K)=WLPF(L,K)+W(L,K)
                  ABLPF(L,K)=ABLPF(L,K)+AB(L,K)                   
              ELSE
                  WLPF(L,K)=WLPF(L,K)+0.0
                  ABLPF(L,K)=ABLPF(L,K)+0.0
              END IF
              
              ABEFF(L,K)=ABEFF(L,K)+AB(L,K)*(SAL(L,K+1)-SAL(L,K))
              WIRT(L,K)=WIRT(L,K)+DT*W(L,K)
              WTLPF(L,K)=WTLPF(L,K)+DT*(FLOAT(NMMT)-0.5)*W(L,K)
            END DO
          END DO
        END IF

        DO K=1,KC                                     ! DRBC, LZ, 12/20/2020
          DO L=2,LA 
            IF(ISCDRY(L).EQ.0) THEN 
              DISPRLPF(L,K)=DISPRLPF(L,K)+DISPRATE(L,K)
            ELSE
              DISPRLPF(L,K)=DISPRLPF(L,K)+0.0
            END IF
          END DO
        END DO

        DO K=1,KC
          DO L=2,LA
            LS=LSC(L)
            SFLLPF(L,K)=SFLLPF(L,K)+SFL(L,K)
            DYELPF(L,K)=DYELPF(L,K)+DYE(L,K)
            SALLPF(L,K)=SALLPF(L,K)+SAL(L,K)
            
            IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
              TEMLPF(L,K)=TEMLPF(L,K)+TEM(L,K)
            ELSE
              TEMLPF(L,K)=TEMLPF(L,K)+TATMT(L)
            END IF
            IF((ISCDRY(L).EQ.1).OR.(ISCDRY(L-1).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
              AHULPF(L,K)=AHULPF(L,K)+0.0
              UHLPF(L,K)=UHLPF(L,K)+0.0
            ELSE
              AHULPF(L,K)=AHULPF(L,K)+0.5*(AH(L,K)+AH(L-1,K))
              UHLPF(L,K)=UHLPF(L,K)+UHDY(L,K)/DYU(L)
            END IF
            IF((ISCDRY(L).EQ.1).OR.(ISCDRY(LS).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
              AHVLPF(L,K)=AHVLPF(L,K)+0.0
              VHLPF(L,K)=VHLPF(L,K)+0.0
            ELSE
              AHVLPF(L,K)=AHVLPF(L,K)+0.5*(AH(L,K)+AH(LS,K))
              VHLPF(L,K)=VHLPF(L,K)+VHDX(L,K)/DXV(L)
            END IF
            
            UIRT(L,K)=UIRT(L,K)+DT*U(L,K)
            ULPF(L,K)=ULPF(L,K)+U(L,K)
            UTLPF(L,K)=UTLPF(L,K)+DT*(FLOAT(NMMT)-0.5)*U(L,K)
            QSUMLPF(L,K)=QSUMLPF(L,K)+QSUM(L,K)
            VIRT(L,K)=VIRT(L,K)+DT*V(L,K)
            VLPF(L,K)=VLPF(L,K)+V(L,K)
            VTLPF(L,K)=VTLPF(L,K)+DT*(FLOAT(NMMT)-0.5)*V(L,K)
          END DO
        END DO

        DO NT=1,NTOX
          DO K=1,KC
            DO L=2,LA
              TOXLPF(L,K,NT)=TOXLPF(L,K,NT)+TOX(L,K,NT)
            END DO
          END DO
        END DO
        DO NSC=1,NSED
          DO K=1,KC
            DO L=2,LA
              SEDLPF(L,K,NSC)=SEDLPF(L,K,NSC)+SED(L,K,NSC)
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KC
            DO L=2,LA
              SNDLPF(L,K,NSN)=SNDLPF(L,K,NSN)+SND(L,K,NSN)
            END DO
          END DO
        END DO

        DO NT=1,NTOX
          DO NS=1,NSED+NSND
            DO K=1,KC
              DO L=1,LC
                TXPFLPF(L,K,NS,NT)=TXPFLPF(L,K,NS,NT)+TOXPFW(L,K,NS,NT)
              END DO
            END DO
          END DO
        END DO

        DO NS=1,NQSER
          DO K=1,KC
            QSRTLPP(K,NS)=QSRTLPP(K,NS)+MAX(QSERT(K,NS),0.)
            QSRTLPN(K,NS)=QSRTLPN(K,NS)+MIN(QSERT(K,NS),0.)
          END DO
        END DO

        DO NS=1,NQCTL
          DO K=1,KC
            QCTLTLP(K,NS)=QCTLTLP(K,NS)+QCTLT(K,NS)
          END DO
        END DO

        DO NMD=1,MDCHH
          QCHNULP(NMD)=QCHNULP(NMD)+QCHANU(NMD)
          QCHNVLP(NMD)=QCHNVLP(NMD)+QCHANV(NMD)
        END DO

        DO NWR=1,NQWR
          QWRSERTLP(NWR)=QWRSERTLP(NWR)+QWRSERT(NWR)
        END DO

        DO K=1,KS
          DO L=2,LA
            LS=LSC(L)
            VPX(L,K)=VPX(L,K)+0.25*(V(L,K+1)+V(L,K))*(WIRT(L,K)+WIRT(LS,K))
            VPY(L,K)=VPY(L,K)+0.25*(W(L,K)+W(L-1,K))*(UIRT(L,K+1)+UIRT(L,K))
          END DO
        END DO

        DO K=1,KC
          DO L=2,LA
            LS=LSC(L)
            VPZ(L,K)=VPZ(L,K)+0.25*(U(L,K)+U(LS,K))*(VIRT(L,K)+VIRT(L-1,K))
          END DO
        END DO

      END IF
!
!**********************************************************************C
!
! **  CHECK FOR END OF FILTER
!
      IF (NMMT.LT.NTSMMT) GO TO 200
!
!**********************************************************************C
!
! **  COMPLETE THE FILTERING
!
!----------------------------------------------------------------------C
!
      FLTWT=1./FLOAT(NTSMMT)

      IF(ISICM.GE.1) FLTWT=2.*FLTWT

      IF(NTSMMT.LT.NTSPTC) THEN

        DO L=2,LA
          HLPF(L)=FLTWT*HLPF(L)
          QSUMELPF(L)=FLTWT*QSUMELPF(L)
          UELPF(L)=FLTWT*UELPF(L)
          VELPF(L)=FLTWT*VELPF(L)
          RAINLPF(L)=FLTWT*RAINLPF(L)
          EVPSLPF(L)=FLTWT*EVPSLPF(L)
          EVPGLPF(L)=FLTWT*EVPGLPF(L)
          RINFLPF(L)=FLTWT*RINFLPF(L)
          GWLPF(L)=FLTWT*GWLPF(L)
        END DO

        DO NSC=1,NSED
          DO K=1,KB
            DO L=2,LA
              SEDBLPF(L,K,NSC)=SEDBLPF(L,K,NSC)*FLTWT
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KB
            DO L=2,LA
              SNDBLPF(L,K,NSN)=SNDBLPF(L,K,NSN)*FLTWT
            END DO
          END DO
        END DO
        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
              TOXBLPF(L,K,NT)=TOXBLPF(L,K,NT)*FLTWT
            END DO
          END DO
        END DO

        DO K=1,KS
          DO L=2,LA
            ABLPF(L,K)=FLTWT*ABLPF(L,K)
            ABEFF(L,K)=FLTWT*ABEFF(L,K)
            WLPF(L,K)=FLTWT*WLPF(L,K)
          END DO
        END DO

        DO K=1,KC                              !DRBC, LZ, 12/20/2020
          DO L=2,LA
            DISPRLPF(L,K)=FLTWT*DISPRLPF(L,K)
          END DO
        END DO

        DO K=1,KC
          DO L=2,LA
            AHULPF(L,K)=AHULPF(L,K)*FLTWT
            AHVLPF(L,K)=AHVLPF(L,K)*FLTWT
            SALLPF(L,K)=SALLPF(L,K)*FLTWT
            TEMLPF(L,K)=TEMLPF(L,K)*FLTWT
            SFLLPF(L,K)=SFLLPF(L,K)*FLTWT
            DYELPF(L,K)=DYELPF(L,K)*FLTWT
            UHLPF(L,K)=FLTWT*UHLPF(L,K)
            VHLPF(L,K)=FLTWT*VHLPF(L,K)
            ULPF(L,K)=FLTWT*ULPF(L,K) 
            VLPF(L,K)=FLTWT*VLPF(L,K) 
            QSUMLPF(L,K)=FLTWT*QSUMLPF(L,K)
          END DO
        END DO

        DO NSC=1,NSED
          DO K=1,KC
            DO L=2,LA
              SEDLPF(L,K,NSC)=SEDLPF(L,K,NSC)*FLTWT
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KC
            DO L=2,LA
              SNDLPF(L,K,NSN)=SNDLPF(L,K,NSN)*FLTWT
            END DO
          END DO
        END DO
        DO NT=1,NTOX
          DO K=1,KC
            DO L=2,LA
              TOXLPF(L,K,NT)=TOXLPF(L,K,NT)*FLTWT
            END DO
          END DO
        END DO

        DO NT=1,NTOX
          DO NS=1,NSED+NSND
            DO K=1,KC
              DO L=1,LC
                TXPFLPF(L,K,NS,NT)=TXPFLPF(L,K,NS,NT)*FLTWT
              END DO
            END DO
          END DO
        END DO

        DO NS=1,NQSER
          DO K=1,KC
            QSRTLPP(K,NS)=FLTWT*QSRTLPP(K,NS)
            QSRTLPN(K,NS)=FLTWT*QSRTLPN(K,NS)
          END DO
        END DO

        DO NS=1,NQCTL
          DO K=1,KC
            QCTLTLP(K,NS)=FLTWT*QCTLTLP(K,NS)
          END DO
        END DO

        DO NMD=1,MDCHH
          QCHNULP(NMD)=FLTWT*QCHNULP(NMD)
          QCHNVLP(NMD)=FLTWT*QCHNVLP(NMD)
        END DO

        DO NWR=1,NQWR
          QWRSERTLP(NWR)=FLTWT*QWRSERTLP(NWR)
        END DO

      ELSE

        DO L=2,LA
          HLPF(L)=FLTWT*HLPF(L)
          QSUMELPF(L)=FLTWT*QSUMELPF(L)
          UELPF(L)=FLTWT*UELPF(L)
          VELPF(L)=FLTWT*VELPF(L)
          RAINLPF(L)=FLTWT*RAINLPF(L)
          EVPSLPF(L)=FLTWT*EVPSLPF(L)
          EVPGLPF(L)=FLTWT*EVPGLPF(L)
          RINFLPF(L)=FLTWT*RINFLPF(L)
          GWLPF(L)=FLTWT*GWLPF(L)
        END DO

        DO NSC=1,NSED
          DO K=1,KB
            DO L=2,LA
              SEDBLPF(L,K,NSC)=SEDBLPF(L,K,NSC)*FLTWT
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KB
            DO L=2,LA
              SNDBLPF(L,K,NSN)=SNDBLPF(L,K,NSN)*FLTWT
            END DO
          END DO
        END DO
        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
              TOXBLPF(L,K,NT)=TOXBLPF(L,K,NT)*FLTWT
            END DO
          END DO
        END DO

        DO K=1,KS
          DO L=2,LA
            ABLPF(L,K)=FLTWT*ABLPF(L,K)
            ABEFF(L,K)=FLTWT*ABEFF(L,K)
            VPX(L,K)=FLTWT*VPX(L,K)
            VPY(L,K)=FLTWT*VPY(L,K)
            WLPF(L,K)=FLTWT*WLPF(L,K)
            WTLPF(L,K)=FLTWT*WTLPF(L,K)
          END DO
        END DO

        DO K=1,KC                              !DRBC, LZ, 12/20/2020
          DO L=2,LA
            DISPRLPF(L,K)=FLTWT*DISPRLPF(L,K)
          END DO
        END DO

        DO K=1,KC
          DO L=2,LA
            AHULPF(L,K)=AHULPF(L,K)*FLTWT
            AHVLPF(L,K)=AHVLPF(L,K)*FLTWT
            SALLPF(L,K)=FLTWT*SALLPF(L,K)
            TEMLPF(L,K)=FLTWT*TEMLPF(L,K)
            SFLLPF(L,K)=FLTWT*SFLLPF(L,K)
            DYELPF(L,K)=FLTWT*DYELPF(L,K)
            UHLPF(L,K)=FLTWT*UHLPF(L,K)
            ULPF(L,K)=FLTWT*ULPF(L,K)
            UTLPF(L,K)=FLTWT*UTLPF(L,K)
            VHLPF(L,K)=FLTWT*VHLPF(L,K)
            QSUMLPF(L,K)=FLTWT*QSUMLPF(L,K)
            VLPF(L,K)=FLTWT*VLPF(L,K)
            VTLPF(L,K)=FLTWT*VTLPF(L,K)
            VPZ(L,K)=FLTWT*VPZ(L,K)
          END DO
        END DO

        DO NSC=1,NSED
          DO K=1,KC
            DO L=2,LA
              SEDLPF(L,K,NSC)=SEDLPF(L,K,NSC)*FLTWT
            END DO
          END DO
        END DO
        DO NSN=1,NSND
          DO K=1,KC
            DO L=2,LA
              SNDLPF(L,K,NSN)=SNDLPF(L,K,NSN)*FLTWT
            END DO
          END DO
        END DO
        DO NT=1,NTOX
          DO K=1,KC
            DO L=2,LA
              TOXLPF(L,K,NT)=TOXLPF(L,K,NT)*FLTWT
            END DO
          END DO
        END DO

        DO NT=1,NTOX
          DO NS=1,NSED+NSND
            DO K=1,KC
              DO L=1,LC
                TXPFLPF(L,K,NS,NT)=TXPFLPF(L,K,NS,NT)*FLTWT
              END DO
            END DO
          END DO
        END DO


        DO NS=1,NQSER
          DO K=1,KC
            QSRTLPP(K,NS)=FLTWT*QSRTLPP(K,NS)
            QSRTLPN(K,NS)=FLTWT*QSRTLPN(K,NS)
          END DO
        END DO

        DO NS=1,NQCTL
          DO K=1,KC
            QCTLTLP(K,NS)=FLTWT*QCTLTLP(K,NS)
          END DO
        END DO

        DO NMD=1,MDCHH
          QCHNULP(NMD)=FLTWT*QCHNULP(NMD)
          QCHNVLP(NMD)=FLTWT*QCHNVLP(NMD)
        END DO

        DO NWR=1,NQWR
          QWRSERTLP(NWR)=FLTWT*QWRSERTLP(NWR)
        END DO
!
!     CALCULATE THE VECTOR POTENTIAL COMPONENTS
!
        DO K=1,KS
          DO L=2,LA
            LS=LSC(L)
            VPX(L,K)=VPX(L,K) -0.25*(VTLPF(L,K+1)+VTLPF(L,K))*(WLPF(L,K)+WLPF(LS,K))
            VPY(L,K)=VPY(L,K) -0.25*(WTLPF(L,K)+WTLPF(L-1,K))*(ULPF(L,K+1)+ULPF(L,K))
          END DO
        END DO

        DO K=1,KC
          DO L=2,LA
            LS=LSC(L)
            VPZ(L,K)=VPZ(L,K) -0.25*(UTLPF(L,K)+UTLPF(LS,K))*(VLPF(L,K)+VLPF(L-1,K))
            VPZ(L,K)=VPZ(L,K)*HMC(L)*SUB(L)*SUB(LS)*SVB(L)*SVB(L-1)
          END DO
        END DO
!
!     CALCULATE VECTOR POTENTIAL TRANSPORT VELOCITY
!
        DO K=1,KC
          DO L=2,LA
            LS=LSC(L)
            LN=LNC(L)
            UVPT(L,K)=(VPZ(LN,K)-VPZ(L,K))/DYU(L) -DZIC(K)*(VPY(L,K)-VPY(L,K-1))
            VVPT(L,K)=DZIC(K)*(VPX(L,K)-VPX(L,K-1)) -(VPZ(L+1,K)-VPZ(L,K))/DXV(L)
          END DO
        END DO
!
        DO K=1,KS
          DO L=2,LA
            LS=LSC(L)
            LN=LNC(L)
            WVPT(L,K)=(VPY(L+1,K)-VPY(L,K))/DXP(L)-(VPX(LN,K)-VPX(L,K))/DYP(L)
          END DO
        END DO

      END IF
!
!     ADJUST TRANSPORTS AT TIDAL ELEVATION BOUNDARY CELLS
!
      QXW=0.
      QXWVP=0.

      DO K=1,KC
        DO LL=1,NPBW
          L=LPBW(LL)
          QXW=QXW+UHLPF(L+1,K)*DZC(K)*DYU(L+1)
          QXWVP=QXWVP+UVPT(L+1,K)*DZC(K)*DYU(L+1)
        END DO
      END DO

      QXE=0.
      QXEVP=0.
      DO K=1,KC
        DO LL=1,NPBE
          L=LPBE(LL)
          QXE=QXE+UHLPF(L,K)*DZC(K)*DYU(L)
          QXEVP=QXEVP+UVPT(L,K)*DZC(K)*DYU(L)
        END DO
      END DO

      QYS=0.
      QYSVP=0.
      DO K=1,KC
        DO LL=1,NPBS
          L=LPBS(LL)
          LN=LNC(L)
          QYS=QYS+VHLPF(LN,K)*DZC(K)*DXV(LN)
          QYSVP=QYSVP+VVPT(LN,K)*DZC(K)*DXV(LN)
        END DO
      END DO

      QYN=0.
      QYNVP=0.
      DO K=1,KC
        DO LL=1,NPBN
          L=LPBN(LL)
          LN=LNC(L)
          QYN=QYN+VHLPF(L,K)*DZC(K)*DXV(L)
          QYNVP=QYNVP+VVPT(L,K)*DZC(K)*DXV(L)
        END DO
      END DO
!
!**********************************************************************C
!
! **  OUTPUT RESIDUAL TRANSPORT TO FILE restran.out
!
      IF (ISSSMMT.EQ.1.AND.N.LT.NTS) GO TO 198
!
!----------------------------------------------------------------------C
!
      IF (ISRESTR.EQ.1) THEN

        IF (JSRESTR.EQ.1) THEN
          OPEN(98,FILE='RESTRAN.OUT',STATUS='UNKNOWN')
          CLOSE(98,STATUS='DELETE')
          OPEN(98,FILE='RESTRAN.OUT',STATUS='UNKNOWN')
          JSRESTR=0
        ELSE
          OPEN(98,FILE='RESTRAN.OUT',ACCESS='APPEND',STATUS='UNKNOWN')
        END IF

        IF(NTSMMT.LT.NTSPTC) THEN
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            WRITE(98,907)HMP(L),HLPF(L),QSUMELPF(L)
            WRITE(98,907)(UHLPF(L,K),K=1,KC)
            WRITE(98,907)(VHLPF(L,K),K=1,KC)
            WRITE(98,907)(AHULPF(L,K),K=1,KC)
            WRITE(98,907)(AHVLPF(L,K),K=1,KC)
            WRITE(98,907)(SALLPF(L,K),K=1,KC)
            WRITE(98,907)(ABLPF(L,K),K=1,KS)
            WRITE(98,907)(ABEFF(L,K),K=1,KS)
          END DO
        ELSE
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            WRITE(98,907)HMP(L),HLPF(L),QSUMELPF(L)
            WRITE(98,907)(UHLPF(L,K),K=1,KC)
            WRITE(98,907)(VHLPF(L,K),K=1,KC)
            WRITE(98,907)(VPZ(L,K),K=1,KC)
            WRITE(98,907)(AHULPF(L,K),K=1,KC)
            WRITE(98,907)(AHVLPF(L,K),K=1,KC)
            WRITE(98,907)(SALLPF(L,K),K=1,KC)
            WRITE(98,907)(VPX(L,K),K=1,KS)
            WRITE(98,907)(VPY(L,K),K=1,KS)
            WRITE(98,907)(ABLPF(L,K),K=1,KS)
!      WRITE(98,907)(ABEFF(L,K),K=1,KS)
          END DO
        END IF

        CLOSE(98)

      END IF

  907 FORMAT(12E12.4)
!  
!**********************************************************************C
!
! **  OUTPUT TO WASP COMPATIABLE FILES
!
      IF(ISWASP.gt.4.and.igridv.eq.0) CALL WASPHYDROLINK 
      IF(ISWASP.gt.4.and.igridv.eq.1) CALL WASPHYDROLINKgvc
!
      IF(ISRCA.GE.1) CALL RCAHQ
      IF(ISICM.GE.1) CALL CEQICM
!
!**********************************************************************C
!
  198 CONTINUE
!
!**********************************************************************C
!
! **  WRITE GRAPHICS FILES FOR RESIDUAL VARIABLES
!
      IF (ISSSMMT.EQ.1.AND.N.LT.NTS) GO TO 199
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SALINITY CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
!
      IF (ISRSPH(1).EQ.1.AND.ISTRAN(1).GE.1) THEN
           CALL RSALPLTH(1,SALLPF)
      END IF

      IF (ISRSPH(2).EQ.1.AND.ISTRAN(2).GE.1) THEN
           CALL RSALPLTH(2,TEMLPF)
      END IF

      IF (ISRSPH(3).EQ.1.AND.ISTRAN(3).GE.1) THEN
           CALL RSALPLTH(3,DYELPF)
      END IF

      IF (ISRSPH(4).EQ.1.AND.ISTRAN(4).GE.1) THEN
           CALL RSALPLTH(4,SFLLPF)
      END IF

      DO K=2,KB
        DO L=2,LA
          SEDBTLPF(L,K)=0.
          SNDBTLPF(L,K)=0.
        END DO
      END DO

      DO K=1,KC
        DO L=2,LA
          TVAR1S(L,K)=TOXLPF(L,K,1)
          SEDTLPF(L,K)=0.
          SNDTLPF(L,K)=0.
        END DO
      END DO

      IF (ISRSPH(5).EQ.1.AND.ISTRAN(5).GE.1) THEN
        DO NT=1,NTOX
           CALL RSALPLTH(5,TVAR1S)
        END DO
      END IF

      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            SEDBTLPF(L,K)=SEDBTLPF(L,K)+SEDBLPF(L,K,NS)
          END DO
        END DO
      END DO

      DO NS=1,NSED
        DO K=1,KC
          DO L=2,LA
            SEDTLPF(L,K)=SEDTLPF(L,K)+SEDLPF(L,K,NS)
          END DO
        END DO
      END DO

      IF (ISRSPH(6).EQ.1.AND.ISTRAN(6).GE.1) THEN
        DO NSC=1,NSED
            CALL RSALPLTH(6,SEDTLPF)
        END DO
      END IF

      DO NS=1,NSND
        DO K=1,KB
          DO L=2,LA
            SNDBTLPF(L,K)=SNDBTLPF(L,K)+SNDBLPF(L,K,NS)
          END DO
        END DO
      END DO

      DO NS=1,NSND
        DO K=1,KC
          DO L=2,LA
            SNDTLPF(L,K)=SNDTLPF(L,K)+SNDLPF(L,K,NS)
          END DO
        END DO
      END DO

      IF (ISRSPH(7).EQ.1.AND.ISTRAN(7).GE.1) THEN
        DO NSN=1,NSND
           CALL RSALPLTH(7,SNDTLPF)
        END DO
      END IF
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES:
! **  SUBROUTINE RVELPLTH
!
      IF(ISRVPH.GE.1) CALL RVELPLTH
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SURFACE ELEVATION PLOTTING IN HORIZONTAL PLANES:
! **  SUBROUTINE RVELPLTH
!
      IF(ISRPPH.EQ.1) CALL RSURFPLT
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SALINITY AND VERTICAL MASS DIFFUSIVITY CONTOURING IN
! **  3 VERTICAL PLANES:  SUBROUTINE RSALPLTV
!
      DO ITMP=1,7
      IF(ISRSPV(ITMP).GE.1) CALL RSALPLTV(ITMP)
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL NORMAL AND TANGENTIAL VELOCITY CONTOURING AND AND
! **  TANGENTIAL VELOCITY VECTOR PLOTTING IN VERTICAL PLANES:
! **  SUBROUTINE RVELPLTV
!
      IF(ISRVPV.GE.1) CALL RVELPLTV
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL 3D SCALAR AND VECTOR OUTPUT FILES
!
      IF(ISR3DO.GE.1) CALL ROUT3D
!
!----------------------------------------------------------------------C
!
  199 CONTINUE
!
!**********************************************************************C
!
!     RESET COUNTER
!
      NMMT=0
   200 CONTINUE

      IF(ISICM.GE.1) THEN
        NMMT=NMMT+2
      ELSE
        NMMT=NMMT+1
      END IF
!
! write debug file
!
!      OPEN(98,FILE='W&D_DEBUG.OUT',POSITION='APPEND',STATUS='UNKNOWN')                  !HNR DEBUG
!      L=LIJ(7,2)
!      TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
!      WRITE(98,201)7,2,5,L,N,TIMTMP,HP(L),HDRY,HWET,HMIN,UHDY(L,5),
!     &VHDX(L,5),ISCDRY(L),ISDRY,ITIMSOL,IS2TIM,IGRIDV
!201   FORMAT(4I3,I9,F12.6,F7.4,3F7.2,2F10.6,I10,4I8)
!      CLOSE(98)
!      
!----------------------------------------------------------------------c
! ** end SUBROUTINE CALMMT
!----------------------------------------------------------------------c
      RETURN
      END