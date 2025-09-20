!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CALFQC (ISTL,IS2TL,MVAR,M,CON,CON1,FQCPAD,QSUMPAD,QSUMNAD) 
!
!                                    MVAR=VARIALBE ID 1-8
!                                    M=ARRAY STORAGE ID FOR TIME SERIES
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 7/2023 HUGO RODRIGUEZ  ADDED TIME SERIES OF LOAD FOR SALT,TEMP(HEAT),DYE AND SHELLFISH
!----------------------------------------------------------------------C
!
! **  SUBROUTINE CALFQC CALCULATES MASS SOURCES AND SINKS ASSOCIATED
! **  WITH CONSTANT AND TIME SERIES INFLOWS AND OUTFLOWS; CONTROL 
! **  STRUCTURE INFLOWS AND OUTLOWS; WITHDRAWAL AND RETURN STRUCTURE 
! **  OUTFLOWS; AND  EMBEDED CHANNEL INFLOWS AND OUTFLOWS
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION CON(LCM,KCM),CON1(LCM,KCM),CONQ(LCM,KCM)
      DIMENSION FQCPAD(LCM,KCM),QSUMPAD(LCM,KCM),QSUMNAD(LCM,KCM)
!
!**********************************************************************C
!
! **  INITIALIZE VOLUMETRIC SOURCE-SINK FLUXES AND AUXILLARY VARIABLES
!
!----------------------------------------------------------------------C
!
      IF(ISTL.EQ.3)THEN
      DO K=1,KC
       DO L=1,LC
         CONQ(L,K)=CON(L,K)
       ENDDO 
      ENDDO
	ENDIF
!
      IF(ISTL.EQ.2.AND.IS2TL.EQ.0)THEN
      DO K=1,KC
       DO L=1,LC
         CONQ(L,K)=0.5*(CON(L,K)+CON1(L,K))
       ENDDO 
      ENDDO
	ENDIF
!
      IF(ISTL.EQ.2.AND.IS2TL.EQ.1)THEN
      DO K=1,KC
       DO L=1,LC
         CONQ(L,K)=0.5*(3.*CON(L,K)-CON1(L,K))
       ENDDO 
      ENDDO
	ENDIF
!
      DO K=1,KC
       DO L=1,LC
       FQC(L,K)=0.
       FQCPAD(L,K)=0
	 QSUMPAD(L,K)=0.
	 QSUMNAD(L,K)=0.
       ENDDO 
      ENDDO
!
      IF(MVAR.EQ.4) GOTO 1000     !SHELLFISH LARVAE
      IF(MVAR.EQ.8) GOTO 1500
!
!**********************************************************************C
!
! **  CORRECT CONCENTRATION FOR WATER SURFACE ELEVATION DATA ASSIMILATION
!
      IF(ISWSEDA.GT.0)THEN
        DO K=1,KC
          DO L=1,LC
            FQC(L,K)=FQC(L,K)+QWSEDA(L,K)*CONQ(L,K)
	      QSUMPAD(L,K)=QSUMPAD(L,K)+MAX(QWSEDA(L,K),0.)
	      QSUMNAD(L,K)=QSUMNAD(L,K)+MIN(QWSEDA(L,K),0.)
          ENDDO
        ENDDO
      ENDIF
!
!**********************************************************************C
!
      IF(ISTL.EQ.2.AND.IS2TL.EQ.1)THEN
!
!     CALCULATE FOR TWO TIME LEVEL INTEGRATION
!
!----------------------------------------------------------------------C
!
! **  STANDARD VOLUMETRICS SOURCE SINK LOCATIONS (2TL)
!
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       TMPVAL=QFACTOR(NS)*QSERT(K,NQSTMP)
       FQC(L,K)=FQC(L,K)+MAX(QSS(K,NS),0.)*CQS(K,NS,M)+MIN(QSS(K,NS),0.)*CONQ(L,K)+MAX(TMPVAL,0.) &
                *CSERT(K,NCSTMP,M)+MIN(TMPVAL,0.)*CONQ(L,K)
!    &         +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
!    &         +MIN(QSS(K,NS),0.)*0.5*(CON(L,K)+CON1(L,K)) 
!    &         +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
!    &         +MIN(QSERT(K,NQSTMP),0.)*0.5*(CON(L,K)+CON1(L,K)) 
        FQCPAD(L,K)=FQCPAD(L,K)+MAX(QSS(K,NS),0.)*CQS(K,NS,M)+MAX(TMPVAL,0.)*CSERT(K,NCSTMP,M)
	  QSUMPAD(L,K)=QSUMPAD(L,K)+MAX(QSS(K,NS),0.)+MAX(TMPVAL,0.)
	  QSUMNAD(L,K)=QSUMNAD(L,K)+MIN(QSS(K,NS),0.)+MIN(TMPVAL,0.)
	  ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  START EXTERNAL TIME SERIES LOAD LOCATIONS (2TL)    HNR_GHD 7/2023
!
!LOAD ASSIGNMENT IS DONE SEPARATELY BY VARIABLE 
        IF(MVAR.EQ.1) THEN  !SALT LOAD
          DO NL=1,NLIJ(1)
            L=LLDS(NL,1)
            NLSTMP=NSERL(NL,1)
            DO K=1,KC
              FQC(L,K)=FQC(L,K)+SERLT(K,NLSTMP,1)    !LOAD IN KG/S OF SALT  
              FQCPAD(L,K)=FQCPAD(L,K)+SERLT(K,NLSTMP,1)
            END DO
          END DO
        END IF
        
!      IF(MVAR.EQ.2) !TEMPERATURE IS DONE DIFFERENTLY BY INCREASING T BY DELTATEMP DIRECTLY
      
        IF(MVAR.EQ.3) THEN
          DO NL=1,NLIJ(3)
            L=LLDS(NL,3)
            NLSTMP=NSERL(NL,3)
            DO K=1,KC
              FQC(L,K)=FQC(L,K)+SERLT(K,NLSTMP,3)    !LOAD IN KG/S OF DYE  
              FQCPAD(L,K)=FQCPAD(L,K)+SERLT(K,NLSTMP,3) 
            END DO
        END DO
      END IF
!
! **  END EXTERNAL LOADS LOCATIONS (2TL)    HNR_GHD 7/2023
!
!----------------------------------------------------------------------C
!
! **  JET-PLUME VOLUMETRICS SOURCE SINK LOCATIONS (2TL)
!
      IF(NQJPIJ.GT.0)THEN
      DO NJP=1,NQJPIJ
!
      IF(ICALJP(NJP).EQ.1)THEN
	 RPORTS=FLOAT(NPORTJP(NJP))
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
! QVJPTMP=TIME SERIES DISCHARGE FROM JET-PLUME
       QVJPTMP=0.
       DO K=1,KC
        QVJPTMP=QVJPTMP+QSERT(K,NQSERJP(NJP))
       ENDDO
! QCJPTMP=ENTRAINMENT FLUX
	 QCJPTMP=0.
	 QVJPENT=0.
! REMOVE ENTRAINMENT FLUX AND CALCULATE TOTAL ENTRAIMENT FLUX
	 DO K=1,KC
	   FQC(LJP,K)=FQC(LJP,K)-RPORTS*QJPENT(K,NJP)*CONQ(LJP,K)
	   QCJPTMP=QCJPTMP+QJPENT(K,NJP)*CONQ(LJP,K)
	   QVJPENT=QVJPENT+QJPENT(K,NJP)
	   QSUMNAD(LJP,K)=QSUMNAD(LJP,K)-RPORTS*QJPENT(K,NJP)
	 ENDDO
! PLACE JET FLUX AND ENTRAINMENT FLUX IS EFFECTIVE LAYER
       FQC(LJP,KTMP)=FQC(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QQCJP(NJP)*CQCJP(1,NJP,M)+RPORTS*QVJPTMP &
                     *CSERT(1,NCSERJP(NJP,M),M)
       FQCPAD(LJP,KTMP)=FQCPAD(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QQCJP(NJP)*CQCJP(1,NJP,M)+RPORTS   &
                        *QVJPTMP*CSERT(1,NCSERJP(NJP,M),M)
	 QSUMPAD(LJP,KTMP)=QSUMPAD(LJP,KTMP)+RPORTS*QVJPENT+RPORTS*QQCJP(NJP)+RPORTS*QVJPTMP
      ENDIF
!           
      IF(ICALJP(NJP).EQ.2)THEN
	 RPORTS=FLOAT(NPORTJP(NJP))
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
	 NS=NQWRSERJP(NJP)
       LU=LIJ(IUPCJP(NJP),JUPCJP(NJP))
       KU=KUPCJP(NJP)
       CONUP=CONQ(LU,KU)
	 QCJPTMP=0.
	 QVJPENT=0.
! REMOVE ENTRAIMENT FLUX AND CALCULATE TOTAL ENTRAINMENT
	 DO K=1,KC
	   FQC(LJP,K)=FQC(LJP,K)-RPORTS*QJPENT(K,NJP)*CONQ(LJP,K)
	   QCJPTMP=QCJPTMP+QJPENT(K,NJP)*CONQ(LJP,K)
	   QVJPENT=QVJPENT+QJPENT(K,NJP)
	   QSUMNAD(LJP,K)=QSUMNAD(LJP,K)-RPORTS*QJPENT(K,NJP)
	 ENDDO
!  PLACE ENTRAINMENT, CONSTANT AND TIME SERIES FLUXES IN EFFECTIVE CELL
       FQC(LJP,KTMP)=FQC(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QWRCJP(NJP)*(CWRCJP(NJP,M)+CONUP)        &
                     +RPORTS*QWRSERT(NS)*(CQWRSERT(NS,M)+CONUP)
       FQCPAD(LJP,KTMP)=FQCPAD(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QWRCJP(NJP)*(CWRCJP(NJP,M)+CONUP)  &
                        +RPORTS*QWRSERT(NS)*(CQWRSERT(NS,M)+CONUP)
	 QSUMPAD(LJP,KTMP)=QSUMPAD(LJP,KTMP)+RPORTS*QVJPENT+RPORTS*QWRCJP(NJP)+RPORTS*QWRSERT(NS)
! REMOVAL WITHDRAWAL FROM UPSTREAM CELL
       FQC(LU,KU)=FQC(LU,KU)-RPORTS*QWRCJP(NJP)*CONUP-RPORTS*QWRSERT(NS)*CONUP
	 QSUMNAD(LU,KU)=QSUMNAD(LU,KU)-RPORTS*QWRCJP(NJP)-RPORTS*QWRSERT(NS)
      ENDIF
!
      ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  CONTROL STRUCTURES (2TL)
!
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF(ID.EQ.0.AND.JD.EQ.0)THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      ENDIF
       DO K=1,KC
       FQC(LU,K)=FQC(LU,K)-QCTLT(K,NCTL)*CONQ(LU,K)
       FQC(LD,K)=FQC(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K) 
       FQCPAD(LD,K)=FQCPAD(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K) 
       QSUMPAD(L,K)=QSUMPAD(L,K)+RQWD*QCTLT(K,NCTL)
       QSUMNAD(L,K)=QSUMNAD(L,K)-QCTLT(K,NCTL)
!      FQC(LU,K)=FQC(LU,K)
!    &          -QCTLT(K,NCTL)*0.5*(CON(LU,K)+CON1(LU,K)) 
!      FQC(LD,K)=FQC(LD,K)
!    &          +RQWD*QCTLT(K,NCTL)*0.5*(CON(LU,K)+CON1(LU,K)) 
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  WITHDRAWAL CONCENTRATION RISE RETURN (2TL)
!
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
       NQSTMP=NQWRSERQ(NWR)
       NCSTMP=NQWRSERQ(NWR)
       FQC(LU,KU)=FQC(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))*CONQ(LU,KU)
       FQC(LD,KD)=FQC(LD,KD)+QWR(NWR)*(CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)*(CONQ(LU,KU)      &
                  +CQWRSERT(NCSTMP,M))
       FQCPAD(LD,KD)=FQCPAD(LD,KD)+QWR(NWR)*(CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)*(CONQ(LU,KU)&
                     +CQWRSERT(NCSTMP,M))
       QSUMPAD(LD,KD)=QSUMPAD(LD,KD)+QWR(NWR)+QWRSERT(NQSTMP)
       QSUMNAD(LU,KU)=QSUMNAD(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))
!      FQC(LU,KU)=FQC(LU,K)
!    &          -(QWR(K,NWR)+QSERT(K,NQSTMP))*0.5*(CON(LU,K)
!    &                                            +CON1(LU,K)) 
!      FQC(LD,KD)=FQC(LD,K)
!    &          +QWR(K,NWR)*(0.5*(CON(LU,K)+CON1(LU,K))+CQWR(K,NWR,M))
!    &          +QSERT(K,NQSTMP)*(0.5*(CON(LU,K)+CON1(LU,K))
!    &                                +CSERT(K,NCSTMP,M))
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  SUBGRID SCALE CHANNEL EXCHANGE (2TL)
!
      IF(MDCHH.GE.1)THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF(MDCHTYP(NMD).EQ.1)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.2)THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.3)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        ENDIF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)+MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)+MIN(QUKTMP,0.)               &
                       *CONQ(LMDCHHT,K)+MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)+MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)-MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)-MIN(QUKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)-MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)-MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
!    &             +0.5*MAX(QUKTMP,0.)*(CON1(LMDCHUT,K)+CON(LMDCHUT,K))
!    &             +0.5*MIN(QUKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
!    &             +0.5*MAX(QVKTMP,0.)*(CON1(LMDCHVT,K)+CON(LMDCHVT,K))
!    &             +0.5*MIN(QVKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
!       FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
!    &             -0.5*MAX(QUKTMP,0.)*(CON1(LMDCHUT,K)+CON(LMDCHUT,K))
!    &             -0.5*MIN(QUKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
!       FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
!    &             -0.5*MAX(QVKTMP,0.)*(CON1(LMDCHVT,K)+CON(LMDCHVT,K))
!    &             -0.5*MIN(QVKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
        ENDDO
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  GROUNDWATER, EVAP, RAINFALL (2TL)
!
      IF(ISHOUSATONIC.EQ.0)THEN
      IF(ISGWIE.NE.0)THEN
        DO L=2,LA
         FQC(L,1)=FQC(L,1)-RIFTR(L)*CONQ(L,1)
        ENDDO
      ENDIF
	ENDIF
!
      IF(M.EQ.2)THEN
        IF(ISTOPT(2).EQ.0.OR.ISTOPT(2).EQ.3)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TEMO*DXYP(L)
          ENDDO
        ENDIF
        IF(ISTOPT(2).EQ.1.OR.ISTOPT(2).EQ.2)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TATMT(L)*DXYP(L)
           FQCPAD(L,KC)=FQCPAD(L,KC)+RAINT(L)*TATMT(L)*DXYP(L)
           QSUMPAD(L,KC)=QSUMPAD(L,KC)+RAINT(L)*DXYP(L)
          ENDDO
        ENDIF
      ENDIF
!
      IF(M.EQ.2)THEN
        IF(ISTOPT(2).EQ.0)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)-EVAPSW(L)*CONQ(L,KC)
          ENDDO
        ENDIF
      ENDIF
!
!----------------------------------------------------------------------C
!
      ENDIF
!
!**********************************************************************C
!
      IF(ISTL.EQ.2.AND.IS2TL.EQ.0)THEN
!
!     CALCULATE FOR TWO TIME LEVEL CORRECTION TO 3 TL INTEGRATION
!
!----------------------------------------------------------------------C
!
! **  STANDARD VOLUMETRICS SOURCE SINK LOCATIONS (2TL)
!
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       TMPVAL=QFACTOR(NS)*QSERT(K,NQSTMP)
       FQC(L,K)=FQC(L,K)+MAX(QSS(K,NS),0.)*CQS(K,NS,M)+MIN(QSS(K,NS),0.)*CONQ(L,K)+MAX(TMPVAL,0.) &
                *CSERT(K,NCSTMP,M)+MIN(TMPVAL,0.)*CONQ(L,K)
!    &         +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
!    &         +MIN(QSS(K,NS),0.)*0.5*(CON(L,K)+CON1(L,K)) 
!    &         +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
!    &         +MIN(QSERT(K,NQSTMP),0.)*0.5*(CON(L,K)+CON1(L,K)) 
        FQCPAD(L,K)=FQCPAD(L,K)+MAX(QSS(K,NS),0.)*CQS(K,NS,M)+MAX(TMPVAL,0.)*CSERT(K,NCSTMP,M)
	  QSUMPAD(L,K)=QSUMPAD(L,K)+MAX(QSS(K,NS),0.)+MAX(TMPVAL,0.)
	  QSUMNAD(L,K)=QSUMNAD(L,K)+MIN(QSS(K,NS),0.)+MIN(TMPVAL,0.)
	  ENDDO
      ENDDO
!----------------------------------------------------------------------C
!
! **  START EXTERNAL TIME SERIES LOAD LOCATIONS (2TL)    HNR_GHD 7/2023
!
!LOAD ASSIGNMENT IS DONE SEPARATELY BY VARIABLE 
        IF(MVAR.EQ.1) THEN  !SALT LOAD
          DO NL=1,NLIJ(1)
            L=LLDS(NL,1)
            NLSTMP=NSERL(NL,1)
            DO K=1,KC
              FQC(L,K)=FQC(L,K)+SERLT(K,NLSTMP,1)    !LOAD IN KG/S OF SALT  
              FQCPAD(L,K)=FQCPAD(L,K)+SERLT(K,NLSTMP,1)
            END DO
          END DO
        END IF
        
!      IF(MVAR.EQ.2) !TEMPERATURE IS DONE DIFFERENTLY BY INCREASING T BY DT DIRECTLY
      
        IF(MVAR.EQ.3) THEN
          DO NL=1,NLIJ(3)
            L=LLDS(NL,3)
            NLSTMP=NSERL(NL,3)
            DO K=1,KC
              FQC(L,K)=FQC(L,K)+SERLT(K,NLSTMP,3)    !LOAD IN KG/S OF DYE  
              FQCPAD(L,K)=FQCPAD(L,K)+SERLT(K,NLSTMP,3) 
            END DO
        END DO
      END IF
!
! **  END EXTERNAL LOADS LOCATIONS (2TL)    HNR_GHD 7/2023
!
!----------------------------------------------------------------------C
!
!----------------------------------------------------------------------C
!
! **  JET-PLUME VOLUMETRICS SOURCE SINK LOCATIONS (2TL)
!
      IF(NQJPIJ.GT.0)THEN
      DO NJP=1,NQJPIJ
!
      IF(ICALJP(NJP).EQ.1)THEN
	 RPORTS=FLOAT(NPORTJP(NJP))
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
! QVJPTMP=TIME SERIES DISCHARGE FROM JET-PLUME
       QVJPTMP=0.
       DO K=1,KC
        QVJPTMP=QVJPTMP+QSERT(K,NQSERJP(NJP))
       ENDDO
! QCJPTMP=ENTRAINMENT FLUX
	 QCJPTMP=0.
	 QVJPENT=0.
! REMOVE ENTRAINMENT FLUX AND CALCULATE TOTAL ENTRAIMENT FLUX
	 DO K=1,KC
	   FQC(LJP,K)=FQC(LJP,K)-RPORTS*QJPENT(K,NJP)*CONQ(LJP,K)
	   QCJPTMP=QCJPTMP+QJPENT(K,NJP)*CONQ(LJP,K)
	   QVJPENT=QVJPENT+QJPENT(K,NJP)
	   QSUMNAD(LJP,K)=QSUMNAD(LJP,K)-RPORTS*QJPENT(K,NJP)
	 ENDDO
! PLACE JET FLUX AND ENTRAINMENT FLUX IS EFFECTIVE LAYER
       FQC(LJP,KTMP)=FQC(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QQCJP(NJP)*CQCJP(1,NJP,M)+RPORTS*QVJPTMP &
                     *CSERT(1,NCSERJP(NJP,M),M)
       FQCPAD(LJP,KTMP)=FQCPAD(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QQCJP(NJP)*CQCJP(1,NJP,M)+RPORTS   &
                        *QVJPTMP*CSERT(1,NCSERJP(NJP,M),M)
	 QSUMPAD(LJP,KTMP)=QSUMPAD(LJP,KTMP)+RPORTS*QVJPENT+RPORTS*QQCJP(NJP)+RPORTS*QVJPTMP
      ENDIF
!           
      IF(ICALJP(NJP).EQ.2)THEN
	 RPORTS=FLOAT(NPORTJP(NJP))
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
	 NS=NQWRSERJP(NJP)
       LU=LIJ(IUPCJP(NJP),JUPCJP(NJP))
       KU=KUPCJP(NJP)
       CONUP=CONQ(LU,KU)
	 QCJPTMP=0.
	 QVJPENT=0.
! REMOVE ENTRAIMENT FLUX AND CALCULATE TOTAL ENTRAINMENT
	 DO K=1,KC
	   FQC(LJP,K)=FQC(LJP,K)-RPORTS*QJPENT(K,NJP)*CONQ(LJP,K)
	   QCJPTMP=QCJPTMP+QJPENT(K,NJP)*CONQ(LJP,K)
	   QVJPENT=QVJPENT+QJPENT(K,NJP)
	   QSUMNAD(LJP,K)=QSUMNAD(LJP,K)-RPORTS*QJPENT(K,NJP)
	 ENDDO
!  PLACE ENTRAINMENT, CONSTANT AND TIME SERIES FLUXES IN EFFECTIVE CELL
       FQC(LJP,KTMP)=FQC(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QWRCJP(NJP)*(CWRCJP(NJP,M)+CONUP)+RPORTS &
                     *QWRSERT(NS)*(CQWRSERT(NS,M)+CONUP)
       FQCPAD(LJP,KTMP)=FQCPAD(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QWRCJP(NJP)*(CWRCJP(NJP,M)+CONUP)  &
                        +RPORTS*QWRSERT(NS)*(CQWRSERT(NS,M)+CONUP)
	 QSUMPAD(LJP,KTMP)=QSUMPAD(LJP,KTMP)+RPORTS*QVJPENT+RPORTS*QWRCJP(NJP)+RPORTS*QWRSERT(NS)
! REMOVAL WITHDRAWAL FROM UPSTREAM CELL
       FQC(LU,KU)=FQC(LU,KU)-RPORTS*QWRCJP(NJP)*CONUP-RPORTS*QWRSERT(NS)*CONUP
	 QSUMNAD(LU,KU)=QSUMNAD(LU,KU)-RPORTS*QWRCJP(NJP)-RPORTS*QWRSERT(NS)
      ENDIF
!
      ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  CONTROL STRUCTURES (2TL)
!
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF(ID.EQ.0.AND.JD.EQ.0)THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      ENDIF
       DO K=1,KC
       FQC(LU,K)=FQC(LU,K)-QCTLT(K,NCTL)*CONQ(LU,K)
       FQC(LD,K)=FQC(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K) 
       FQCPAD(LD,K)=FQCPAD(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K) 
       QSUMPAD(L,K)=QSUMPAD(L,K)+RQWD*QCTLT(K,NCTL)
       QSUMNAD(L,K)=QSUMNAD(L,K)-QCTLT(K,NCTL)
!      FQC(LU,K)=FQC(LU,K)
!    &          -QCTLT(K,NCTL)*0.5*(CON(LU,K)+CON1(LU,K)) 
!      FQC(LD,K)=FQC(LD,K)
!    &          +RQWD*QCTLT(K,NCTL)*0.5*(CON(LU,K)+CON1(LU,K)) 
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  WITHDRAWAL CONCENTRATION RISE RETURN (2TL)
!
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
       NQSTMP=NQWRSERQ(NWR)
       NCSTMP=NQWRSERQ(NWR)
       FQC(LU,KU)=FQC(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))*CONQ(LU,KU)
       FQC(LD,KD)=FQC(LD,KD)+QWR(NWR)*(CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)*(CONQ(LU,KU)      &
                  +CQWRSERT(NCSTMP,M))
       FQCPAD(LD,KD)=FQCPAD(LD,KD)+QWR(NWR)*(CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)*(CONQ(LU,KU)&
                     +CQWRSERT(NCSTMP,M))
       QSUMPAD(LD,KD)=QSUMPAD(LD,KD)+QWR(NWR)+QWRSERT(NQSTMP)
       QSUMNAD(LU,KU)=QSUMNAD(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))
!      FQC(LU,KU)=FQC(LU,K)
!    &          -(QWR(K,NWR)+QSERT(K,NQSTMP))*0.5*(CON(LU,K)
!    &                                            +CON1(LU,K)) 
!      FQC(LD,KD)=FQC(LD,K)
!    &          +QWR(K,NWR)*(0.5*(CON(LU,K)+CON1(LU,K))+CQWR(K,NWR,M))
!    &          +QSERT(K,NQSTMP)*(0.5*(CON(LU,K)+CON1(LU,K))
!    &                                +CSERT(K,NCSTMP,M))
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  SUBGRID SCALE CHANNEL EXCHANGE (2TL)
!
      IF(MDCHH.GE.1)THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF(MDCHTYP(NMD).EQ.1)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.2)THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.3)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        ENDIF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)+MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)+MIN(QUKTMP,0.)               &
                       *CONQ(LMDCHHT,K)+MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)+MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)-MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)-MIN(QUKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)-MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)-MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
!    &             +0.5*MAX(QUKTMP,0.)*(CON1(LMDCHUT,K)+CON(LMDCHUT,K))
!    &             +0.5*MIN(QUKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
!    &             +0.5*MAX(QVKTMP,0.)*(CON1(LMDCHVT,K)+CON(LMDCHVT,K))
!    &             +0.5*MIN(QVKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
!       FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
!    &             -0.5*MAX(QUKTMP,0.)*(CON1(LMDCHUT,K)+CON(LMDCHUT,K))
!    &             -0.5*MIN(QUKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
!       FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
!    &             -0.5*MAX(QVKTMP,0.)*(CON1(LMDCHVT,K)+CON(LMDCHVT,K))
!    &             -0.5*MIN(QVKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
        ENDDO
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  GROUNDWATER, EVAP, RAINFALL (2TL)
!
      IF(ISHOUSATONIC.EQ.0)THEN
      IF(ISGWIE.NE.0)THEN
        DO L=2,LA
         FQC(L,1)=FQC(L,1)-RIFTR(L)*CONQ(L,1)
        ENDDO
      ENDIF
      ENDIF
!
      IF(M.EQ.2)THEN
        IF(ISTOPT(2).EQ.0.OR.ISTOPT(2).EQ.3)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TEMO*DXYP(L)
          ENDDO
        ENDIF
        IF(ISTOPT(2).EQ.1.OR.ISTOPT(2).EQ.2)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TATMT(L)*DXYP(L)
           FQCPAD(L,KC)=FQCPAD(L,KC)+RAINT(L)*TATMT(L)*DXYP(L)
           QSUMPAD(L,KC)=QSUMPAD(L,KC)+RAINT(L)*DXYP(L)
          ENDDO
        ENDIF
      ENDIF
!
      IF(M.EQ.2)THEN
        IF(ISTOPT(2).EQ.0)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)-EVAPSW(L)*CONQ(L,KC)
          ENDDO
        ENDIF
      ENDIF
!
!----------------------------------------------------------------------C
!
      ENDIF
!
!**********************************************************************C
!
      IF(ISTL.EQ.3)THEN
!
!     CALCULATE FOR THREE TIME LEVELS
!
!----------------------------------------------------------------------C
!
! **  STANDARD VOLUMETRICS SOURCE SINK LOCATIONS (3TL)
!
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
	 TMPVAL=QFACTOR(NS)*QSERT(K,NQSTMP)
       FQC(L,K)=FQC(L,K)+MAX(QSS(K,NS),0.)*CQS(K,NS,M)+MIN(QSS(K,NS),0.)*CONQ(L,K)+MAX(TMPVAL,0.) &
                *CSERT(K,NCSTMP,M)+MIN(TMPVAL,0.)*CONQ(L,K)
!    &         +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
!    &         +MIN(QSS(K,NS),0.)*CON(L,K) 
!    &         +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
!    &         +MIN(QSERT(K,NQSTMP),0.)*CON(L,K)
        FQCPAD(L,K)=FQCPAD(L,K)+MAX(QSS(K,NS),0.)*CQS(K,NS,M)+MAX(TMPVAL,0.)*CSERT(K,NCSTMP,M)
	  QSUMPAD(L,K)=QSUMPAD(L,K)+MAX(QSS(K,NS),0.)+MAX(TMPVAL,0.)
	  QSUMNAD(L,K)=QSUMNAD(L,K)+MIN(QSS(K,NS),0.)+MIN(TMPVAL,0.)
       ENDDO
      ENDDO
!----------------------------------------------------------------------C
!
! **  START EXTERNAL TIME SERIES LOAD LOCATIONS (2TL)    HNR_GHD 7/2023
!
!LOAD ASSIGNMENT IS DONE SEPARATELY BY VARIABLE 
        IF(MVAR.EQ.1) THEN  !SALT LOAD
          DO NL=1,NLIJ(1)
            L=LLDS(NL,1)
            NLSTMP=NSERL(NL,1)
            DO K=1,KC
              FQC(L,K)=FQC(L,K)+SERLT(K,NLSTMP,1)    !LOAD IN KG/S OF SALT  
              FQCPAD(L,K)=FQCPAD(L,K)+SERLT(K,NLSTMP,1)
            END DO
          END DO
        END IF
        
!      IF(MVAR.EQ.2) !TEMPERATURE IS DONE DIFFERENTLY BY INCREASING T BY DT DIRECTLY
      
        IF(MVAR.EQ.3) THEN
          DO NL=1,NLIJ(3)
            L=LLDS(NL,3)
            NLSTMP=NSERL(NL,3)
            DO K=1,KC
              FQC(L,K)=FQC(L,K)+SERLT(K,NLSTMP,3)    !LOAD IN KG/S OF DYE  
              FQCPAD(L,K)=FQCPAD(L,K)+SERLT(K,NLSTMP,3) 
            END DO
        END DO
      END IF
!
! **  END EXTERNAL LOADS LOCATIONS (2TL)    HNR_GHD 7/2023
!
!----------------------------------------------------------------------C
!
!----------------------------------------------------------------------C
!
! **  JET-PLUME VOLUMETRICS SOURCE SINK LOCATIONS (3TL)
!
! **  JET-PLUME VOLUMETRICS SOURCE SINK LOCATIONS (2TL)
!
      IF(NQJPIJ.GT.0)THEN
      DO NJP=1,NQJPIJ
!
      IF(ICALJP(NJP).EQ.1)THEN
	 RPORTS=FLOAT(NPORTJP(NJP))
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
! QVJPTMP=TIME SERIES DISCHARGE FROM JET-PLUME
       QVJPTMP=0.
       DO K=1,KC
        QVJPTMP=QVJPTMP+QSERT(K,NQSERJP(NJP))
       ENDDO
! QCJPTMP=ENTRAINMENT FLUX
	 QCJPTMP=0.
	 QVJPENT=0.
! REMOVE ENTRAINMENT FLUX AND CALCULATE TOTAL ENTRAIMENT FLUX
	 DO K=1,KC
	   FQC(LJP,K)=FQC(LJP,K)-RPORTS*QJPENT(K,NJP)*CONQ(LJP,K)
	   QCJPTMP=QCJPTMP+QJPENT(K,NJP)*CONQ(LJP,K)
	   QVJPENT=QVJPENT+QJPENT(K,NJP)
	   QSUMNAD(LJP,K)=QSUMNAD(LJP,K)-RPORTS*QJPENT(K,NJP)
	 ENDDO
! PLACE JET FLUX AND ENTRAINMENT FLUX IS EFFECTIVE LAYER
       FQC(LJP,KTMP)=FQC(LJP,KTMP)+RPORTS*QCJPTMP                                                 &! ENTRAINMENT
                     +RPORTS*QQCJP(NJP)*CQCJP(1,NJP,M)                                            &! CONSTANT DISCHARGE
                     +RPORTS*QVJPTMP*CSERT(1,NCSERJP(NJP,M),M)                                     !TIME SERIES DISCHARGE
       FQCPAD(LJP,KTMP)=FQCPAD(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QQCJP(NJP)*CQCJP(1,NJP,M)+RPORTS   &
                        *QVJPTMP*CSERT(1,NCSERJP(NJP,M),M)
	 QSUMPAD(LJP,KTMP)=QSUMPAD(LJP,KTMP)+RPORTS*QVJPENT+RPORTS*QQCJP(NJP)+RPORTS*QVJPTMP
      ENDIF
!           
!           
      IF(ICALJP(NJP).EQ.2)THEN
	 RPORTS=FLOAT(NPORTJP(NJP))
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
	 NS=NQWRSERJP(NJP)
       LU=LIJ(IUPCJP(NJP),JUPCJP(NJP))
       KU=KUPCJP(NJP)
       CONUP=CONQ(LU,KU)
	 QCJPTMP=0.
	 QVJPENT=0.
! REMOVE ENTRAIMENT FLUX AND CALCULATE TOTAL ENTRAINMENT
	 DO K=1,KC
	   FQC(LJP,K)=FQC(LJP,K)-RPORTS*QJPENT(K,NJP)*CONQ(LJP,K)
	   QCJPTMP=QCJPTMP+QJPENT(K,NJP)*CONQ(LJP,K)
	   QVJPENT=QVJPENT+QJPENT(K,NJP)
	   QSUMNAD(LJP,K)=QSUMNAD(LJP,K)-RPORTS*QJPENT(K,NJP)
	 ENDDO
!  PLACE ENTRAINMENT, CONSTANT AND TIME SERIES FLUXES IN EFFECTIVE CELL
       FQC(LJP,KTMP)=FQC(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QWRCJP(NJP)*(CWRCJP(NJP,M)+CONUP)+RPORTS &
                     *QWRSERT(NS)*(CQWRSERT(NS,M)+CONUP)
       FQCPAD(LJP,KTMP)=FQCPAD(LJP,KTMP)+RPORTS*QCJPTMP+RPORTS*QWRCJP(NJP)*(CWRCJP(NJP,M)+CONUP)  &
                        +RPORTS*QWRSERT(NS)*(CQWRSERT(NS,M)+CONUP)
	 QSUMPAD(LJP,KTMP)=QSUMPAD(LJP,KTMP)+RPORTS*QVJPENT+RPORTS*QWRCJP(NJP)+RPORTS*QWRSERT(NS)
! REMOVAL WITHDRAWAL FROM UPSTREAM CELL
       FQC(LU,KU)=FQC(LU,KU)-RPORTS*QWRCJP(NJP)*CONUP-RPORTS*QWRSERT(NS)*CONUP
	 QSUMNAD(LU,KU)=QSUMNAD(LU,KU)-RPORTS*QWRCJP(NJP)-RPORTS*QWRSERT(NS)
      ENDIF
!
      ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  CONTROL STRUCTURES (3TL)
!
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF(ID.EQ.0.AND.JD.EQ.0)THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      ENDIF
       DO K=1,KC
       FQC(LU,K)=FQC(LU,K)-QCTLT(K,NCTL)*CONQ(LU,K)
       FQC(LD,K)=FQC(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K)
       FQCPAD(LD,K)=FQCPAD(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K) 
       QSUMPAD(L,K)=QSUMPAD(L,K)+RQWD*QCTLT(K,NCTL)
       QSUMNAD(L,K)=QSUMNAD(L,K)-QCTLT(K,NCTL)
!      FQC(LU,K)=FQC(LU,K)
!    &          -QCTLT(K,NCTL)*CON(LU,K)
!      FQC(LD,K)=FQC(LD,K)
!    &          +RQWD*QCTLT(K,NCTL)*CON(LU,K)
       ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  WITHDRAWAL CONCENTRATION RISE RETURN (3TL)
!
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
       NQSTMP=NQWRSERQ(NWR)
       NCSTMP=NQWRSERQ(NWR)
       FQC(LU,KU)=FQC(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))*CONQ(LU,KU)
       FQC(LD,KD)=FQC(LD,KD)+QWR(NWR)*(CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)*(CONQ(LU,KU)+     &
                  CQWRSERT(NCSTMP,M))
       FQCPAD(LD,KD)=FQCPAD(LD,KD)+QWR(NWR)*(CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)*(CONQ(LU,KU)&
                     +CQWRSERT(NCSTMP,M))
       QSUMPAD(LD,KD)=QSUMPAD(LD,KD)+QWR(NWR)+QWRSERT(NQSTMP)
       QSUMNAD(LU,KU)=QSUMNAD(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))
!      FQC(LU,K)=FQC(LU,K)
!    &          -(QWR(K,NWR)+QSERT(K,NQSTMP))*CON(LU,K)
!      FQC(LD,K)=FQC(LD,K)
!    &          +QWR(K,NWR)*(CON(LU,K)+CQWR(K,NWR,M))
!    &          +QSERT(K,NQSTMP)*(CON(LU,K)+CSERT(K,NCSTMP,M))
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  SUBGRID SCALE CHANNEL EXCHANGE (3TL)
!
      IF(MDCHH.GE.1)THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF(MDCHTYP(NMD).EQ.1)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.2)THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.3)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        ENDIF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)                      &
                    +MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)        &          
                    +MIN(QUKTMP,0.)*CONQ(LMDCHHT,K)        &
                    +MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)        &
                    +MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)-MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)-MIN(QUKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)-MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)-MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
!       FQC(LMDCHHT,K)=FQC(LMDCHHT,K)
!    &               +MAX(QUKTMP,0.)*CON(LMDCHUT,K)
!    &               +MIN(QUKTMP,0.)*CON(LMDCHHT,K)
!    &               +MAX(QVKTMP,0.)*CON(LMDCHVT,K)
!    &               +MIN(QVKTMP,0.)*CON(LMDCHHT,K)
!       FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
!    &               -MAX(QUKTMP,0.)*CON(LMDCHUT,K)
!    &               -MIN(QUKTMP,0.)*CON(LMDCHHT,K)
!       FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
!    &               -MAX(QVKTMP,0.)*CON(LMDCHVT,K)
!    &               -MIN(QVKTMP,0.)*CON(LMDCHHT,K)
        ENDDO
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  GROUNDWATER, EVAP, RAINFALL (3TL)
!
      IF(ISHOUSATONIC.EQ.0)THEN
      IF(ISGWIE.NE.0)THEN
        DO L=2,LA
         FQC(L,1)=FQC(L,1)-RIFTR(L)*CONQ(L,1)
        ENDDO
      ENDIF
      ENDIF
!
      IF(M.EQ.2)THEN
        IF(ISTOPT(2).EQ.0.OR.ISTOPT(2).EQ.3)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TEMO*DXYP(L)
          ENDDO
        ENDIF
        IF(ISTOPT(2).EQ.1.OR.ISTOPT(2).EQ.2)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TATMT(L)*DXYP(L)
           FQCPAD(L,KC)=FQCPAD(L,KC)+RAINT(L)*TATMT(L)*DXYP(L)
           QSUMPAD(L,KC)=QSUMPAD(L,KC)+RAINT(L)*DXYP(L)
          ENDDO
        ENDIF
      ENDIF
!
      IF(M.EQ.2)THEN
        IF(ISTOPT(2).EQ.0)THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)-EVAPSW(L)*CONQ(L,KC)
          ENDDO
        ENDIF
      ENDIF
!
!----------------------------------------------------------------------C
!
      ENDIF
!
      GOTO 2000
!
!**********************************************************************C
!
!     TWO TIME LEVEL SFL SOURCE-SINK INTERPOLATION
!
 1000 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       FQC(L,K)=FQC(L,K)+MAX(QSS(K,NS),0.)*CQS(K,NS,M)+MIN(QSS(K,NS),0.)*CONQ(L,K)                &
                +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)+MIN(QSERT(K,NQSTMP),0.)*CONQ(L,K)
       ENDDO
      ENDDO
! **  START EXTERNAL TIME SERIES LOAD SFL (2TL)    HNR_GHD 7/2023
!
        IF(MVAR.EQ.4) THEN
          DO NL=1,NLIJ(4)
            L=LLDS(NL,4)
            NLSTMP=NSERL(NL,4)
            DO K=1,KC
              FQC(L,K)=FQC(L,K)+SERLT(K,NLSTMP,4)    !LOAD IN KG/S OF SFL
              FQCPAD(L,K)=FQCPAD(L,K)+SERLT(K,NLSTMP,4) 
            END DO
        END DO
      END IF
! **  END EXTERNAL TIME SERIES LOAD SFL (2TL)    HNR_GHD 7/2023
      
      
!
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF(ID.EQ.0.AND.JD.EQ.0)THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      ENDIF
       DO K=1,KC
       FQC(LU,K)=FQC(LU,K)-QCTLT(K,NCTL)*CONQ(LU,K)
       FQC(LD,K)=FQC(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K)
       ENDDO
      ENDDO
!
      DO NWR=1,NQWR
      IU=IQWRU(NWR)
      JU=JQWRU(NWR)
      KU=KQWRU(NWR)
      ID=IQWRD(NWR)
      JD=JQWRD(NWR)
      KD=KQWRD(NWR)
      LU=LIJ(IU,JU)
      LD=LIJ(ID,JD)
      NQSTMP=NQWRSERQ(NWR)
      NCSTMP=NQWRSERQ(NWR)
       FQC(LU,KU)=FQC(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))*CONQ(LU,KU)
       FQC(LD,KD)=FQC(LD,KD)+QWR(NWR)*(SFLKILL*CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)*(SFLKILL  &
                  *CONQ(LU,KU)+CQWRSERT(NCSTMP,M))
!      FQC(LU,K)=FQC(LU,K)
!    &          -(QWR(K,NWR)+QSERT(K,NQSTMP))*CON(LU,K)
!      FQC(LD,K)=FQC(LD,K)
!    &          +QWR(K,NWR)*(CON(LU,K)+CQWR(K,NWR,M))
!    &          +QSERT(K,NQSTMP)*(CON(LU,K)+CSERT(K,NCSTMP,M))
      ENDDO
!
      IF(MDCHH.GE.1)THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF(MDCHTYP(NMD).EQ.1)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.2)THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.3)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        ENDIF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)+MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)+MIN(QUKTMP,0.)               &
                       *CONQ(LMDCHHT,K)+MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)+MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)-MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)-MIN(QUKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)-MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)-MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
        ENDDO
        ENDDO
      ENDIF
!
      GOTO 2000
!
!**********************************************************************C
!
!     TWO TIME LEVEL WATER QUALITY SOURCE-SINK INTERPOLATION
!
 1500 CONTINUE
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
       DO L=1,LC
         CONQ(L,K)=CON1(L,K)
       ENDDO 
      ENDDO
!
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
       DO K=1,KC
	 TMPVAL=QFACTOR(NS)*QSERT(K,NQSTMP)
       FQC(L,K)=FQC(L,K)+MIN(QSS(K,NS),0.)*CONQ(L,K)+MIN(TMPVAL,0.)*CONQ(L,K)
	 QSUMPAD(L,K)=QSUMPAD(L,K)+MAX(QSS(K,NS),0.)+MAX(TMPVAL,0.)
	 QSUMNAD(L,K)=QSUMNAD(L,K)+MIN(QSS(K,NS),0.)+MIN(TMPVAL,0.)
       ENDDO
      ENDDO
!
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF(ID.EQ.0.AND.JD.EQ.0)THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      ENDIF
       DO K=1,KC
       FQC(LU,K)=FQC(LU,K)-QCTLT(K,NCTL)*CONQ(LU,K)
       FQC(LD,K)=FQC(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K)
       FQCPAD(LD,K)=FQCPAD(LD,K)+RQWD*QCTLT(K,NCTL)*CONQ(LU,K) 
       QSUMPAD(L,K)=QSUMPAD(L,K)+RQWD*QCTLT(K,NCTL)
       QSUMNAD(L,K)=QSUMNAD(L,K)-QCTLT(K,NCTL)
       ENDDO
      ENDDO
!
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
       NQSTMP=NQWRSERQ(NWR)
       NCSTMP=NQWRSERQ(NWR)
        FQC(LU,KU)=FQC(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))*CONQ(LU,KU)
        FQC(LD,KD)=FQC(LD,KD)+QWR(NWR)*(CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)*(CONQ(LU,KU)     &
                   +CQWRSERT(NCSTMP,M))
        FQCPAD(LD,KD)=FQCPAD(LD,KD)+QWR(NWR)*(CONQ(LU,KU)+CQWR(NWR,M))+QWRSERT(NQSTMP)            &
                      *(CONQ(LU,KU)+CQWRSERT(NCSTMP,M))
        QSUMPAD(LD,KD)=QSUMPAD(LD,KD)+QWR(NWR)+QWRSERT(NQSTMP)
        QSUMNAD(LU,KU)=QSUMNAD(LU,KU)-(QWR(NWR)+QWRSERT(NQSTMP))
      ENDDO
!
      IF(MDCHH.GE.1)THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF(MDCHTYP(NMD).EQ.1)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.2)THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        ENDIF
        IF(MDCHTYP(NMD).EQ.3)THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        ENDIF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)                    &
                    +MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)      &
                    +MIN(QUKTMP,0.)*CONQ(LMDCHHT,K)      &
                    +MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)      &
                    +MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)-MAX(QUKTMP,0.)*CONQ(LMDCHUT,K)-MIN(QUKTMP,0.)*CONQ(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)-MAX(QVKTMP,0.)*CONQ(LMDCHVT,K)-MIN(QVKTMP,0.)*CONQ(LMDCHHT,K)
        ENDDO
        ENDDO
      ENDIF
!
!**********************************************************************C
!
 2000 CONTINUE
!
!      DO K=1,KC
!      TMPVAL=DZIC(K)
!       DO L=1,LC
!       FQC(L,K)=TMPVAL*FQC(L,K)
!       FQCPAD(L,K)=TMPVAL*FQCPAD(L,K)
!       QSUMPAD(L,K)=TMPVAL*QSUMPAD(L,K)
!       QSUMNAD(L,K)=TMPVAL*QSUMNAD(L,K)
!       ENDDO 
!      ENDDO
!
!**********************************************************************C
!
      RETURN
      END