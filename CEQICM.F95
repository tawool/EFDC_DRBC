!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE CEQICM
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
! **  SUBROUTINE FOR INTERFACING CE-QUAL-ICM MODEL
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!      PARAMETER ( NQINFLM=96 )
!
      DIMENSION DZRCA(KCM+1),DZZRCA(KCM+1),QINRCA(NQSIJM),QINTFL(NQINFLM,KCM),TMPARA(ICM,JCM)
      DIMENSION TMPICMF(2*LCM*KCM)
!
!**********************************************************************C
!
! **  READ CONTROL FILES AND WRITE TIME INVARIANT FILES ON FIRST ENTRY
!
      IF(JSWASP.EQ.0) GOTO 1000
!
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      ELSE
        TIME=TIMESEC/86400.
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  READ MAIN CONTROL FILE
!
      OPEN(1,FILE='EFDC.ICM',STATUS='UNKNOWN')
!
      DO NSKIP=1,9
       READ(1,100)
      ENDDO
      READ(1,*)ISDICM,IAUXICM,ISTICM,IDMPCL,JDMPCL,NCICM,NFICM
!
      CLOSE(1)
!
!----------------------------------------------------------------------C
!
! **  READ CELL MAPPING FILE
!
      OPEN(1,FILE='EFDC_C_ICM.INP',STATUS='UNKNOWN')
!
      DO NSKIP=1,7
       READ(1,100)
      ENDDO
!
      DO M=1,NCICM
       READ(1,*)IDUM,JDUM,LCEFDC(M),KCEFDC(M),MDUM
      ENDDO
!
      CLOSE(1)
!
!----------------------------------------------------------------------C
!
! **  READ FLOW MAPPING FILE
!
      OPEN(1,FILE='EFDC_F_ICM.INP',STATUS='UNKNOWN')
!
      DO NSKIP=1,8
       READ(1,100)
      ENDDO
!
      DO M=1,NFICM
       READ(1,*)IDUM,JDUM,LFEFDC(M),KFEFDC(M),MDUM,IDRICM(M)
      ENDDO
!
      CLOSE(1)
!
!----------------------------------------------------------------------C
!
! **  WRITE I,J INDICES DEFINING FLOWS BETWEEN ARBITARY CELLS
! **  (POSTIVE FLOW DIRECTION DEFINED FROM FIRST TO SECOND I,J PAIR)
!
      IF(IAUXICM.GE.1)THEN
!
      OPEN(1,FILE='FLWMAP.INP',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='FLWMAP.INP',STATUS='UNKNOWN')
!
      NINTFL=0
!
      IF(NQSIJ.GT.0)THEN
        DO NS=1,NQSIJ
         NINTFL=NINTFL+1
         WRITE(1,101)IDMPCL,JDMPCL,IQS(NS),JQS(NS)
        ENDDO
      ENDIF
!
      IF(NQCTL.GT.0)THEN
        DO NCTL=1,NQCTL
         NINTFL=NINTFL+1
         IF(IQCTLD(NCTL).GT.0)THEN
           WRITE(1,101) IQCTLU(NCTL),JQCTLU(NCTL),IQCTLD(NCTL),JQCTLD(NCTL)
          ELSE
           WRITE(1,101)IQCTLU(NCTL),JQCTLU(NCTL),IDMPCL,JDMPCL
         ENDIF
        ENDDO
      ENDIF
!
      IF(NQWR.GT.0)THEN
        DO NWR=1,NQWR
         NINTFL=NINTFL+1
         WRITE(1,101)IQWRU(NWR),JQWRU(NWR),IQWRD(NWR),JQWRD(NWR)
        ENDDO
      ENDIF
!
      IF(MDCHH.GT.0)THEN
        DO NMD=1,MDCHH
         IF(IMDCHU(NMD).GT.1)THEN
           NINTFL=NINTFL+1
           WRITE(1,101)IMDCHU(NMD),JMDCHU(NMD),IMDCHH(NMD),JMDCHH(NMD)
         ENDIF
        ENDDO
        DO NMD=1,MDCHH
         IF(IMDCHV(NMD).GT.1)THEN
           NINTFL=NINTFL+1
           WRITE(1,101)IMDCHV(NMD),JMDCHV(NMD),IMDCHH(NMD),JMDCHH(NMD)
         ENDIF
        ENDDO
      ENDIF
!
      CLOSE(1)
!
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE EXTERNAL INFLOW LOCATIONS
!
      OPEN(1,FILE='INFLOWIJ.DAT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='INFLOWIJ.DAT',STATUS='UNKNOWN')
!
      WRITE(1,110)
      WRITE(1,111)
!
      IF(NQSIJ.GE.1)THEN
        DO NS=1,NQSIJ
         WRITE(1,112)NS,IQS(NS),JQS(NS)
        ENDDO
      ENDIF
!
      CLOSE(1)
!
!----------------------------------------------------------------------C
!
! **  INITIALIZE INTERFACE FILES
!
      OPEN(1,FILE='EFDCHYD.INP',FORM='UNFORMATTED',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='EFDCHYD.INP',FORM='UNFORMATTED',STATUS='UNKNOWN')
!
      IF(ISDICM.EQ.1)THEN
        OPEN(2,FILE='EFDCHYD.ASC',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        OPEN(2,FILE='EFDCHYD.ASC',STATUS='UNKNOWN')
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE SIGMA GRID FRACTIONAL LAYER THICKNESS (SURFACE DOWN FOR ICM)
!
      WRITE(1) (DZC(K),K=KC,1,-1)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2001)TIME
        WRITE(2,2002)
        DO K=KC,1,-1
         KICM=KC-K+1
         WRITE(2,293)KICM,DZC(K),K
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE CELL HORIZONTAL SURFACE AREAS (LOOP OVER SURFACE LAYER CELLS)
!
      WRITE(1) (DXYP(L),L=2,LA)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2003)
        DO L=2,LA
         LM=L-1
         WRITE(2,293)LM,DXYP(L),L
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE CELL HORIZONTAL DEPTH AVERAGE FLOW FACE AREAS
!     (LOOP OVER SURFACE LAYER CELLS)
!
      NFICMS=NFICM/KC
!
      DO M=1,NFICMS
       L=LFEFDC(M)
       IF(IDRICM(M).EQ.1) TMPICMF(M)=DYU(L)*HMU(L)
       IF(IDRICM(M).EQ.2) TMPICMF(M)=DXV(L)*HMV(L)
      ENDDO
!
      WRITE(1) (TMPICMF(M),M=1,NFICMS)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2004)
        DO M=1,NFICMS
         WRITE(2,293)M,TMPICMF(M),IDRICM(M),LFEFDC(M)
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE CELL HORIZONTAL X DIRECTION LENGTHS
!     (LOOP OVER SURFACE LAYER CELLS)
!
      WRITE(1) (DXP(L),L=2,LA)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2005)
        DO L=2,LA
         LM=L-1
         WRITE(2,293)LM,DXP(L),L
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE CELL HORIZONTAL Y DIRECTION LENGTHS
!     (LOOP OVER SURFACE LAYER CELLS)
!
      WRITE(1) (DYP(L),L=2,LA)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2006)
        DO L=2,LA
         LM=L-1
         WRITE(2,293)LM,DYP(L),L
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE CELL HORIZONTAL DEPTH AVERAGE INITIAL VOLUMES
!     (LOOP OVER SURFACE LAYER CELLS)
!
      DO L=2,LA
       TVAR3C(L-1)=DXYP(L)*H2WQ(L)
       HTMP(L)=H2WQ(L)
      ENDDO
!
      LTMP=LA-1
      WRITE(1) (TVAR3C(L),L=1,LTMP)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2007)
        DO L=2,LA
         LM=L-1
         WRITE(2,293)LM,TVAR3C(LM),L
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
      CLOSE(1)
      IF(ISDICM.EQ.1) CLOSE(2)
!
!----------------------------------------------------------------------C
!
! **  INITIALIZE OTHER FILES TO RECEIVE TIME VARYING DATA
!
      IF(IAUXICM.EQ.1)THEN
!
      OPEN(1,FILE='INFLOW.DAT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
      OPEN(1,FILE='HYDRLGY.INP',FORM='UNFORMATTED',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
      IF(ISDICM.EQ.1)THEN
      OPEN(1,FILE='HYDRLGY.ASC',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      ENDIF
!
      OPEN(1,FILE='INTFLW.INP',FORM='UNFORMATTED',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
      IF(ISDICM.EQ.1)THEN
      OPEN(1,FILE='INTFLW.ASC',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      ENDIF
!
      ENDIF
!
      IF(ISTICM.EQ.1)THEN
        OPEN(1,FILE='SALICM.INP',FORM='UNFORMATTED',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='TEMICM.INP',FORM='UNFORMATTED',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
!
      IF(ISTICM.EQ.1.AND.ISDICM.EQ.1)THEN
        OPEN(1,FILE='SALICM.ASC',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='TEMICM.ASC',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
!
      OPEN(1,FILE='EFDCFLW.INP',FORM='UNFORMATTED',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
      OPEN(1,FILE='EFDCRME.INP',FORM='UNFORMATTED',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
      IF(ISDICM.EQ.1)THEN
        OPEN(1,FILE='EFDCFLW.ASC',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EFDCRME.ASC',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE SUMMARY DATA TO LOG FILE
!
      OPEN(1,FILE='EFDCICM.LOG',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='EFDCICM.LOG',STATUS='UNKNOWN')
!
      ICICM1=IDMPCL
      ICICM2=JDMPCL
      ICICM3=ISDICM
      ICICM4=IAUXICM
      NCICM1=NINTFL
!
      WRITE(1,2001) TIME
      WRITE(1,102)IC,JC,KC
      WRITE(1,103)NINTFL
      WRITE(1,104)IDMPCL,JDMPCL
      WRITE(1,2002)
      WRITE(1,2003)
      WRITE(1,2004)
      WRITE(1,2005)
      WRITE(1,2006)
      WRITE(1,2007)
!
      CLOSE(1)
!
!----------------------------------------------------------------------C
!
      JSWASP=0
      RETURN
!
!**********************************************************************C
!
 1000 CONTINUE
!
! **  SET VECTOR POTENTIAL TRANSPORT SWITCH TO INCLUDE VECTOR POTENTIAL
!     TRANSPORT IF AVERGING PERIOD NTSMMT IS GE REFERENCE
!     TIME PERIOD NTSPTC
!
      SVPT=0.
      IF(NTSMMT.GE.NTSPTC) SVPT=1.
!
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      ELSE
        TIME=TIMESEC/86400.
      ENDIF
      NMID=N-(NTSMMT/2)
      TIMMID=(DT*FLOAT(NMID)+TCON*TBEGIN)/86400.
!
!----------------------------------------------------------------------C
!
! **  WRITE TIME AT END OF AVERAGING PERIOD TO EFDCICM.LOG
!
      OPEN(1,FILE='EFDCICM.LOG',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(1,106)TIME
      WRITE(1,2008)
      WRITE(1,2009)
      WRITE(1,2010)
      WRITE(1,2011)
      CLOSE(1)
!
!----------------------------------------------------------------------C
!
! **  WRITE INFLOWS AT END OF AVERAGING PERIOD TO INFLOW.DAT
!
      IF(IAUXICM.GE.1)THEN
!
      OPEN(1,FILE='INFLOW.DAT',STATUS='UNKNOWN',POSITION='APPEND')
!
      DO NS=1,NQSIJ
       QINRCA(NS)=0.
      ENDDO
!
      DO NS=1,NQSIJ
       NQSTMP=NQSERQ(NS)
       IF(NQSTMP.GT.0)THEN
         DO K=1,KC
          QINRCA(NS)=QINRCA(NS)+QSRTLPP(K,NQSTMP)+MAX(QSS(K,NS),0.)
         ENDDO
        ELSE
         DO K=1,KC
          QINRCA(NS)=QINRCA(NS)+MAX(QSS(K,NS),0.)
         ENDDO
       ENDIF
      ENDDO
!
      WRITE(1,120)TIME,(QINRCA(NS),NS=1,NQSIJ)
!
      CLOSE(1)
!
      ENDIF
!
!C----------------------------------------------------------------------C
!
! **  WRITE INTERNAL FLOWS TO INTFLW.INP
!
      IF(IAUXICM.GE.1)THEN
!
      OPEN(1,FILE='INTFLW.INP',FORM='UNFORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
!
      IF(ISDICM.EQ.1)THEN
       OPEN(2,FILE='INTFLW.ASC',STATUS='UNKNOWN',POSITION='APPEND')
       WRITE(2,106)TIME
      ENDIF
!
      DO K=1,KC
       DO NS=1,NCRCA1
        QINTFL(NS,K)=0.
       ENDDO
      ENDDO
!
      NINTFL=0
!
      IF(NQSIJ.GT.0)THEN
      DO NS=1,NQSIJ
       NINTFL=NINTFL+1
       NQSTMP=NQSERQ(NS)
       IF(NQSTMP.GT.0)THEN
         DO K=1,KC
          QINTFL(NINTFL,K)=QINTFL(NINTFL,K)+QSRTLPN(K,NQSTMP)+MIN(QSS(K,NS),0.)
         ENDDO
        ELSE
         DO K=1,KC
          QINTFL(NINTFL,K)=QINTFL(NINTFL,K)+MIN(QSS(K,NS),0.)
         ENDDO
       ENDIF
      ENDDO
      ENDIF
!
      IF(NQCTL.GT.0)THEN
        DO NCTL=1,NQCTL
         NINTFL=NINTFL+1
         DO K=1,KC
          QINTFL(NINTFL,K)=QINTFL(NINTFL,K)+QCTLTLP(K,NCTL)
         ENDDO
        ENDDO
      ENDIF
!
      IF(NQWR.GT.0)THEN
        DO NWR=1,NQWR
         NQSTMP=NQWRSERQ(NWR)
         NINTFL=NINTFL+1
         DO K=1,KC
         IF(K.EQ.KQWRU(NWR))THEN
           QINTFL(NINTFL,K)=QINTFL(NINTFL,K)+QWRSERTLP(NQSTMP)+QWR(NWR)
         ENDIF
         ENDDO
        ENDDO
      ENDIF
!
      IF(MDCHH.GT.0)THEN
        DO NMD=1,MDCHH
         IF(IMDCHU(NMD).GT.1)THEN
           NINTFL=NINTFL+1
           DO K=1,KC
            QINTFL(NINTFL,K)=QINTFL(NINTFL,K)+QCHNULP(NMD)*DZC(K)
           ENDDO
         ENDIF
        ENDDO
        DO NMD=1,MDCHH
         IF(IMDCHV(NMD).GT.1)THEN
           NINTFL=NINTFL+1
           DO K=1,KC
            QINTFL(NINTFL,K)=QINTFL(NINTFL,K)+QCHNULP(NMD)*DZC(K)
           ENDDO
         ENDIF
        ENDDO
      ENDIF
!
      WRITE(1) QINTFL
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,213)
        DO NS=1,NINTFL
        WRITE(2,211)NS,(QINTFL(NS,K),K=1,KC)
        ENDDO
      ENDIF
!
      CLOSE(1)
      IF(ISDICM.EQ.1) CLOSE(2)
!
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE HYDRODYNAMIC DATA TO EFDCHDY.INP
!
      OPEN(1,FILE='EFDCHYD.INP',FORM='UNFORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(1) TIME
!
      OPEN(3,FILE='EFDCFLW.INP',FORM='UNFORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(3) TIME
!
      OPEN(13,FILE='EFDCRME.INP',FORM='UNFORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(13) TIME
!
      IF(ISDICM.EQ.1)THEN
        OPEN(2,FILE='EFDCHYD.ASC',STATUS='UNKNOWN',POSITION='APPEND')
        WRITE(2,106) TIME
        OPEN(4,FILE='EFDCFLW.ASC',STATUS='UNKNOWN',POSITION='APPEND')
        WRITE(4,106) TIME
        OPEN(14,FILE='EFDCRME.ASC',STATUS='UNKNOWN',POSITION='APPEND')
        WRITE(14,106) TIME
      ENDIF
!
      TAVGTMP=FLOAT(NTSMMT)*DT
!
!----------------------------------------------------------------------C
!
! **  WRITE CELL HORIZONTAL DEPTH AVERAGE VOLUMES
!     (LOOP OVER SURFACE LAYER CELLS)
!
      DO L=1,LC
       TVAR3E(L)=0.
       TVAR3W(L)=0.
       TVAR3S(L)=0.
       TVAR3N(L)=0.
       TVAR3C(L)=0.
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
       TVAR3W(L)=TVAR3W(L)+DYU(L)*UHLPF(L,K)
       TVAR2W(L,K)=DYU(L)*UHLPF(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
       LE=L+1
       TVAR3E(L)=TVAR3E(L)+DYU(LE)*UHLPF(LE,K)
       TVAR2E(L,K)=DYU(LE)*UHLPF(LE,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
       TVAR3S(L)=TVAR3S(L)+DXV(L)*VHLPF(L,K)
       TVAR2S(L,K)=DXV(L)*VHLPF(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
       LN=LNC(L)
       TVAR3N(L)=TVAR3N(L)+DXV(LN)*VHLPF(LN,K)
       TVAR2N(L,K)=DXV(LN)*VHLPF(LN,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
       TVAR3N(L)=TVAR3N(L)*DZC(K)
       TVAR3S(L)=TVAR3S(L)*DZC(K)
       TVAR3E(L)=TVAR3E(L)*DZC(K)
       TVAR3W(L)=TVAR3W(L)*DZC(K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
       TVAR2N(L,K)=TVAR2N(L,K)*DZC(K)
       TVAR2S(L,K)=TVAR2S(L,K)*DZC(K)
       TVAR2E(L,K)=TVAR2E(L,K)*DZC(K)
       TVAR2W(L,K)=TVAR2W(L,K)*DZC(K)
      ENDDO
      ENDDO
!
!      DO L=2,LA
!       TVAR3C(L)=TVAR3C(L)+RAINLPF(L)-EVPSLPF(L)
!      ENDDO
!
      TMPVAL=DT*FLOAT(NTSMMT)
!
      DO L=2,LA
       WLPF(L,0)=0.
       WLPF(L,KC)=0.
       TVAR2C(L,0)=0.
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
       WLPF(L,K)=WLPF(L,K-1)+TVAR2W(L,K)-TVAR2E(L,K)+TVAR2S(L,K)-TVAR2N(L,K)+QSUMLPF(L,K)
      ENDDO
      ENDDO
!
      DO L=2,LA
       TVAR3C(L)=WLPF(L,KS)+TVAR2W(L,KC)-TVAR2E(L,KC)+TVAR2S(L,KC)-TVAR2N(L,KC)+QSUMLPF(L,KC)
      ENDDO
!
      DO K=1,KS
      DO L=2,LA
       TVAR2C(L,K)=TVAR2C(L,K-1)+DZC(K)*TVAR3C(L)
      ENDDO
      ENDDO
!
!
      DO K=1,KS
      DO L=2,LA
       WLPF(L,K)=WLPF(L,K)-TVAR2C(L,K)
      ENDDO
      ENDDO
!
      DO L=2,LA
       TVAR3C(L)=TMPVAL*TVAR3C(L)
      ENDDO
!
      DO L=2,LA
       TVAR3C(L)=DXYP(L)*HTMP(L)+TVAR3C(L)
      ENDDO
!
      DO L=2,LA
       HTMP(L)=TVAR3C(L)*DXYIP(L)
      ENDDO
!
      WRITE(1) (TVAR3C(L),L=2,LA)
!
      DO L=2,LA
       QSUMLPF(L,KC)=QSUMLPF(L,KC)-RAINLPF(L)+EVPSLPF(L)
      ENDDO
!
      DO K=KC,1,-1
       WRITE(3) (QSUMLPF(L,K),L=2,LA)
      ENDDO
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2008)
        WRITE(4,401)
        DO L=2,LA
         LM=L-1
         TMPVAL=DXYP(L)*HWQ(L)
!C         WRITE(2,294)LM,TVAR3C(L),TMPVAL,L
         WRITE(2,2294)LM,L,HP(L),HWQ(L),HTMP(L)
         WRITE(4,402)LM,L,(QSUMLPF(L,K),K=1,KC)
!         IF(TVAR3C(L).LT.0.0)THEN
!           LP=L+1
!           LN=LNC(L)
!           DO K=1,KC
!             RVAL1= DYU(L)*UHLPF(L,K)*DZC(K)
!             RVAL2=-DYU(LP)*UHLPF(LP,K)*DZC(K)
!             RVAL3= DXV(L)*VHLPF(L,K)*DZC(K)
!             RVAL4=-DXV(LN)*VHLPF(LN,K)*DZC(K)
!             RVAL5= QSUMLPF(L,K)
!             RVAL6=RVAL1+RVAL2+RVAL3+RVAL4+RVAL5
!             WRITE(2,222)L,K,RVAL1,RVAL2,RVAL3,RVAL4,RVAL5,RVAL6
!           ENDDO
!         ENDIF
        ENDDO
      ENDIF
!
      DO L=2,LA
       TVAR3C(L)=RAINLPF(L)-EVPSLPF(L)
      ENDDO
!
      WRITE(13) (TVAR3C(L),L=2,LA)
!
      DO L=2,LA
       LM=L-1
       WRITE(14,402)LM,L,TVAR3C(L)
      ENDDO
!
  222 FORMAT(' ERROR ',2I5,6F12.2)
  401 FORMAT(/,'  LICM     L    QSUMLPF(L,K) K=1,KC',/)
  402 FORMAT(2I6,12E13.5)
 2294 FORMAT(2I6,4F12.5)
!
!----------------------------------------------------------------------C
!
! **  WRITE HORIZONTAL FLOWS
!
      DO M=1,NFICM
       L=LFEFDC(M)
       K=KFEFDC(M)
       IF(IDRICM(M).EQ.1)THEN
         TMPICMF(M)=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
       ENDIF
       IF(IDRICM(M).EQ.2)THEN
         TMPICMF(M)=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
       ENDIF
      ENDDO
!
      WRITE(1) (TMPICMF(M),M=1,NFICM)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2009)
        DO M=1,NFICM
         L=LFEFDC(M)
         K=KFEFDC(M)
         WRITE(2,293)M,TMPICMF(M),IDRICM(M),L,K
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE VERTICAL DIFFUSIVITY
!
      MM=0
      DO K=KS,1,-1
       DO L=2,LA
        MM=MM+1
        TMPICMF(MM)=ABLPF(L,K)
       ENDDO
      ENDDO
!
      WRITE(1) (TMPICMF(M),M=1,MM)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2010)
        MM=0
        DO K=KS,1,-1
         DO L=2,LA
          MM=MM+1
          WRITE(2,293)MM,ABLPF(L,K),L,K
         ENDDO
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  DETERMINATION OF VERTICAL VELOCITY AND VOLUME FLUX
!
! **  LOAD NET DEPTH INTEGRATED INFLOWS INTO TVAR3E AND
! **  OUT FLOWS INTO TVAR3N
!
!JH      DO L=2,LA
!JH       LN=LNC(L)
!JH       TVAR3E(L)=0.
!JH       TVAR3N(L)=0.
!JH       DO K=1,KC
!JH        TVAR3E(L)=TVAR3E(L)+DZC(K)*
!JH     &          ( DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))
!JH     &           +DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K)) )
!JH        TVAR3N(L)=TVAR3N(L)+DZC(K)*
!JH     &          ( DYU(L+1)*(UHLPF(L+1,K)+SVPT*UVPT(L+1,K))
!JH     &           +DXV(LN )*(VHLPF(LN ,K)+SVPT*VVPT(LN ,K)) )
!JH       ENDDO
!JH      ENDDO
!
!      IF(MDCHH.GE.1)THEN
!        DO NMD=1,MDCHH
!        IF(MDCHTYP(NMD).EQ.1)THEN
!          TVAR3E(LMDCHH(NMD))=TVAR3E(LMDCHH(NMD))
!     &                   +QCHNULP(NMD)
!          TVAR3N(LMDCHU(NMD))=TVAR3N(LMDCHU(NMD))
!     &                   +QCHNULP(NMD)
!        ENDIF
!        IF(MDCHTYP(NMD).EQ.2)THEN
!          TVAR3E(LMDCHH(NMD))=TVAR3E(LMDCHH(NMD))
!     &                   +QCHNVLP(NMD)
!          TVAR3N(LMDCHV(NMD))=TVAR3N(LMDCHV(NMD))
!     &                   +QCHNVLP(NMD)
!        ENDIF
!        IF(MDCHTYP(NMD).EQ.3)THEN
!          TVAR3E(LMDCHH(NMD))=TVAR3E(LMDCHH(NMD))
!     &                   +QCHNULP(NMD)+QCHNVLP(NMD)
!          TVAR3N(LMDCHU(NMD))=TVAR3N(LMDCHU(NMD))
!     &                   +QCHNULP(NMD)+QCHNVLP(NMD)
!        ENDIF
!        ENDDO
!      ENDIF
!
!      DO L=2,LA
!       WLPF(L,0)=0.
!      ENDDO
!
!      DO K=1,KC
!       DO L=2,LA
!        TVAR1N(L,K)=0.
!        TVAR1S(L,K)=0.
!        TVAR1E(L,K)=0.
!        TVAR1W(L,K)=0.
!        TVAR2C(L,K)=DZC(K)*QSUMELPF(L)
!       ENDDO
!      ENDDO
!
!      DO K=1,KC
!       DO L=2,LA
!        LN=LNC(L)
!        TVAR1N(L,K)=DXV(LN)*VHLPF(LN,K)
!        TVAR1S(L,K)=DXV(L)*VHLPF(L,K)
!        TVAR1E(L,K)=DYU(L+1)*UHLPF(L+1,K)
!        TVAR1W(L,K)=DYU(L)*UHLPF(L,K)
!       ENDDO
!      ENDDO
!
!      DO K=1,KC
!       DO L=2,LA
!        LN=LNC(L)
!        TVAR1N(L,K)=DZC(K)*TVAR1N(L,K)
!        TVAR1S(L,K)=DZC(K)*TVAR1S(L,K)
!        TVAR1E(L,K)=DZC(K)*TVAR1E(L,K)
!        TVAR1W(L,K)=DZC(K)*TVAR1W(L,K)
!       ENDDO
!      ENDDO
!
!
! **  CALCULATE QZ (STORED IN WLPF(L,K)
!
!      IF(KC.GT.1)THEN
!      DO K=1,KC
!       DO L=2,LA
!        LN=LNC(L)
!         WLPF(L,K)=SWB(L)*( WLPF(L,K-1)
!     &      +TVAR1W(L,K)-TVAR3W(L)
!     &      -TVAR1E(L,K)+TVAR3E(L)
!     &      +TVAR1S(L,K)-TVAR3S(L)
!     &      -TVAR1N(L,K)+TVAR3N(L)
!     &      +QSUMLPF(L,K)-TVAR2C(L,K) )
!       ENDDO
!      ENDDO
!      ENDIF
!
! **  WRITE VERTICAL VOLUME FLUX
!
      MM=0
      DO K=KS,1,-1
       DO L=2,LA
        MM=MM+1
        TMPICMF(MM)=WLPF(L,K)
       ENDDO
      ENDDO
!
      WRITE(1) (TMPICMF(M),M=1,MM)
!
      IF(ISDICM.EQ.1)THEN
        WRITE(2,2011)
        MM=0
        DO K=KS,1,-1
         DO L=2,LA
          MM=MM+1
          WRITE(2,293)MM,WLPF(L,K),L,K
         ENDDO
        ENDDO
        MM=0
        DO L=2,LA
         WRITE(2,293)MM,WLPF(L,KC),L,KC
        ENDDO
      ENDIF
!
!----------------------------------------------------------------------C
!
      CLOSE(1)
      CLOSE(3)
      CLOSE(13)
      IF(ISDICM.EQ.1) CLOSE(2)
      IF(ISDICM.EQ.1) CLOSE(4)
      IF(ISDICM.EQ.1) CLOSE(14)
!
!----------------------------------------------------------------------C
!
! **  WRITE SALINITY
!
      IF(ISTICM.EQ.1)THEN
      OPEN(1,FILE='SALICM.INP',FORM='UNFORMATTED',POSITION='APPEND')
!
      MM=0
      DO K=KS,1,-1
       DO L=2,LA
        MM=MM+1
        TMPICMF(MM)=SALLPF(L,K)
       ENDDO
      ENDDO
!
      WRITE(1) (TMPICMF(M),M=1,MM)
!
      CLOSE(1)
      ENDIF
!
      IF(ISTICM.EQ.1.AND.ISDICM.EQ.1)THEN
      OPEN(1,FILE='SALICM.ASC',POSITION='APPEND')
      WRITE(1,106) TIME
      WRITE(1,2012)
!
      MM=0
      DO K=KS,1,-1
       DO L=2,LA
        MM=MM+1
        WRITE(2,293)MM,SALLPF(L,K),L,K
       ENDDO
      ENDDO
!
      CLOSE(1)
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE TEMPERATURE
!
      IF(ISTICM.EQ.1)THEN
      OPEN(1,FILE='TEMICM.INP',FORM='UNFORMATTED',POSITION='APPEND')
!
      MM=0
      DO K=KS,1,-1
       DO L=2,LA
        MM=MM+1
        TMPICMF(MM)=TEMLPF(L,K)
       ENDDO
      ENDDO
!
      WRITE(1) (TMPICMF(M),M=1,MM)
!!
      CLOSE(1)
      ENDIF
!
      IF(ISTICM.EQ.1.AND.ISDICM.EQ.1)THEN
      OPEN(1,FILE='TEMICM.ASC',POSITION='APPEND')
      WRITE(1,106) TIME
      WRITE(1,2013)
!
      MM=0
      DO K=KS,1,-1
       DO L=2,LA
        MM=MM+1
        WRITE(2,293)MM,TEMLPF(L,K),L,K
       ENDDO
      ENDDO
!
      CLOSE(1)
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  WRITE HYDROLOGY
!
      IF(IAUXICM.EQ.1)THEN
!
      OPEN(1,FILE='HYDRLGY.INP',FORM='UNFORMATTED',POSITION='APPEND')
!
      WRITE(1)TIME
!
      MM=0
      DO L=2,LA
       MM=MM+1
       TVAR3C(MM)=RAINLPF(L)
      ENDDO
!
      WRITE(1) (TVAR3C(M),M=1,MM)
!
      MM=0
      DO L=2,LA
       MM=MM+1
       TVAR3C(MM)=EVPSLPF(L)
      ENDDO
!
      WRITE(1) (TVAR3C(M),M=1,MM)
!
      MM=0
      DO L=2,LA
       MM=MM+1
       TVAR3C(MM)=RINFLPF(L)
      ENDDO
!
      WRITE(1) (TVAR3C(M),M=1,MM)
!
      CLOSE(1)
!
      IF(ISDICM.EQ.1)THEN
        OPEN(2,FILE='HYDRLGY.ASC',POSITION='APPEND')
        WRITE(2,106)TIME
        WRITE(2,212)
!
        DO L=2,LA
         WRITE(2,200)L,IL(L),JL(L),RAINLPF(L),EVPSLPF(L),EVPGLPF(L),RINFLPF(L),GWLPF(L)
        ENDDO
!
        CLOSE(2)
      ENDIF
!
      ENDIF
!
!----------------------------------------------------------------------C
!
  100 FORMAT(120X)
  101 FORMAT(4I10)
  102 FORMAT(/,' NROW,NCOL,NLAYR = ',3I10/)
  103 FORMAT(/,' NO INTERNAL FLOWS, NINTFL (LINES) IN FLWMAP.INP = ',I10/)
  104 FORMAT(/,' ROW, COLUMN INDICES OF DUMP CELL = ',2I10/)
  105 FORMAT(/,' SIMULATION STARTING TIME IN DAYS = ',F12.6/)
  106 FORMAT(/,' TIME IN DAYS AT END OF AVERAGING PERIOD = ',F12.6/)
  110 FORMAT(' LOCATION OF INFLOWS ',/)
  111 FORMAT(' INFLOW #   ROW INDEX   COLUMN INDEX ',/)
  112 FORMAT(2X,I5,7X,I5,7X,I5)
  120 FORMAT(F12.6,13F12.4)
  200 FORMAT(3I5,6E14.6)
  201 FORMAT(' L,I(ROW),J(COL),QX(I,J,K),K=1,KC ',/)
  202 FORMAT(' L,I(ROW),J(COL),QY(I,J,K),K=1,KC ',/)
  203 FORMAT(' L,I(ROW),J(COL),QZ(I,J,K),K=1,KS ',/)
  204 FORMAT(' L,I(ROW),J(COL),AX(I,J,K),K=1,KC ',/)
  205 FORMAT(' L,I(ROW),J(COL),AY(I,J,K),K=1,KC ',/)
  206 FORMAT(' L,I(ROW),J(COL),AZ(I,J,K),K=1,KS ',/)
  207 FORMAT(' L,I(ROW),J(COL),SELS(I,J),SELE(I,J),DSEL(I,J) ',/)
  208 FORMAT(' L,I(ROW),J(COL),SAL(I,J,K),K=1,KC ',/)
  209 FORMAT(' L,I(ROW),J(COL),TEM(I,J,K),K=1,KC ',/)
  210 FORMAT(//)
  211 FORMAT(I5,2X,6E15.6)
  212 FORMAT(' L,I(ROW),J(ROW),RAINLPF(I,J),EVPSLPF(I,J),EVPGLPF(I,J),RINFLPF(I,J),GWLPF(I,J) ',/)
  213 FORMAT(' NQINTFL,QINTFL ',/)
  215 FORMAT(' L,I(ROW),J(COL),SURFELV START AVG INTERVAL',/)
  216 FORMAT(' L,I(ROW),J(COL),DEL SURFELV OVER INTERVAL',/)
  291 FORMAT(I8,F8.4)
  292 FORMAT(I8,E13.5)
  293 FORMAT(I8,E13.5,3I8)
  294 FORMAT(I8,E13.5,E13.5,3I8)
 2001 FORMAT(/,' TIME AT ICM INTERFACE INITIALIZATION = ',F12.4,/)
 2002 FORMAT(/,' SIGMA LAYER FRACTIONAL THICKNESS: KICM, DZ, KEFDC',/)
 2003 FORMAT(/,' HORIZONTAL CELL SURFACE AREAS, TOP LAYER : LICM, AREA',' LEFDC',/)
 2004 FORMAT(/,' DEPTH INTEGRATED HORIZONTAL FLOW FACE AREAS: NF, AREA',' IDIR, LEFDC',/)
 2005 FORMAT(/,' X DIRECTION CELL LENTHS, TOP LAYER: LICM, DX, LEFDC',/)
 2006 FORMAT(/,' Y DIRECTION CELL LENTHS, TOP LAYER: LICM, DY, LEFDC',/)
 2007 FORMAT(/,' INITIAL DEPTH AVERGE CELL VOLUMES: LICM, VOL, LEFDC',/)
 2008 FORMAT(/,' DEPTH AVERGE CELL VOLUMES: LICM, VOL, LEFDC',/)
 2009 FORMAT(/,' HORIZONTAL FLOWS: NF, FLOW, IDIR, LEFDC, KEFDC'/)
 2010 FORMAT(/,' VERTICAL DIFFUSIVITY: ND, DIFF, LEFDC, KEFDC',/)
 2011 FORMAT(/,' VERTICAL SIGMA FLOW: ND, FLOWZ, LEFDC, KEFDC',/)
 2012 FORMAT(/,' HYDRO SALINITY: NICM, SAL, LEFDC, KEFDC',/)
 2013 FORMAT(/,' HYDRO SALINITY: NICM, TEM, LEFDC, KEFDC',/)
!
!**********************************************************************C
!
      RETURN
      END
