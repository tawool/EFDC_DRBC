C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALQVS (ISTL)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C
C 4/10/06           Hugo N Rodriguez   Add a lake correcting flow
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
C ** SUBROUTINE CALQVS UPDATES TIME VARIABLE VOLUME SOURCES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      DELT=DT2
      IF(ISTL.EQ.2) DELT=DT
      IF(ISDYNSTP.GT.0) DELT=DTDYN

C
C**********************************************************************C
C
C **  INITIALIZE NULL (0) FLOW SERIES
C
      GWSERT(0)=0.
      QWRSERT(0)=0.
      DO K=1,KC
        QSERT(K,0)=0.
        QCTLT(K,0)=0.
        QCTLTO(K,0)=0.
      ENDDO
C
      NCTMP=4+NSED+NSND+NTOX
      DO NC=1,NCTMP
        GWCSERT(0,NC)=0.
      ENDDO
C
      DO L=2,LA
        QGW(L)=0.0
      END DO
C
      DO NC=1,NCTMP
        DO L=2,LA
          CONGW(L,NC)=0.0
        END DO
      END DO
C
C **  INITIALIZE TOTAL FLOW SERIES
C
      DO L=1,LC
      QSUME(L)=0.
      ENDDO
C
      DO K=1,KC
       DO L=1,LC
       QSUM(L,K)=0.
       ENDDO
      ENDDO
C
C **  ADD FLOW CORRECT FOR WATER SURFACE ELEVATION DATA ASSIMILATION
C
      IF(ISWSEDA.GT.0)THEN
        DO K=1,KC
          DO L=1,LC
            QSUM(L,K)=QWSEDA(L,K)
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  VOLUME SOURCE/SINK INTERPOLATION
C
      DO NS=1,NQSER
C
        IF(ISTL.EQ.2)THEN
          IF(ISDYNSTP.EQ.0)THEN
            TIME=DT*(FLOAT(N)-0.5)/TCQSER(NS)
     &            +TBEGIN*(TCON/TCQSER(NS))
          ELSE
            TIME=TIMESEC/TCQSER(NS)
          ENDIF
        ELSE
          IF(ISDYNSTP.EQ.0)THEN
            TIME=DT*FLOAT(N-1)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
          ELSE
            TIME=TIMESEC/TCQSER(NS)
          ENDIF
        ENDIF
C
        M1=MQTLAST(NS)
  100   CONTINUE
        M2=M1+1
        IF(TIME.GT.TQSER(M2,NS))THEN
          M1=M2
          GOTO 100
        ELSE
          MQTLAST(NS)=M1
        ENDIF
C
        TDIFF=TQSER(M2,NS)-TQSER(M1,NS)
        WTM1=(TQSER(M2,NS)-TIME)/TDIFF
        WTM2=(TIME-TQSER(M1,NS))/TDIFF
        DO K=1,KC
          QSERT(K,NS)=WTM1*QSER(M1,K,NS)+WTM2*QSER(M2,K,NS)
        ENDDO
C
      ENDDO
C
      IF(N.EQ.1)THEN
        DO LL=1,NQSIJ
          L=LQS(LL)
          ITYP=LCT(L)
          IF(ITYP.LE.0.OR.ITYP.GE.8)THEN
            WRITE(6,6111)LL,IQS(LL),JQS(LL),L,ITYP
!            WRITE(8,6111)LL,IQS(LL),JQS(LL),L,ITYP   !hnr 7/27/2009
          ENDIF
        ENDDO
      ENDIF
C
      DO LL=1,NQSIJ
        NS=NQSERQ(LL)
        L=LQS(LL)
        DO K=1,KC
	    IF(LGVCP(L,K))
     &    QSUM(L,K)=QSUM(L,K)+RQSMUL(LL)*(QSS(K,LL)
     &                                   +QFACTOR(LL)*QSERT(K,NS))
        ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  GROUNDWATER SOURCE/SINK INTERPOLATION
C
      IF(NGWSER.GE.1)THEN
c
        NCTMP=4+NSED+NSND+NTOX
C
        DO NS=1,NGWSER
C
          IF(ISTL.EQ.2)THEN
            IF(ISDYNSTP.EQ.0)THEN
              TIME=DT*(FLOAT(N)-0.5)/TCGWSER(NS)
     &            +TBEGIN*(TCON/TCGWSER(NS))
            ELSE
              TIME=TIMESEC/TCGWSER(NS)
            ENDIF
          ELSE
            IF(ISDYNSTP.EQ.0)THEN
              TIME=DT*FLOAT(N-1)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
            ELSE
              TIME=TIMESEC/TCGWSER(NS)
            ENDIF
          ENDIF
C
          M1=MGWTLAST(NS)
  700     CONTINUE
          M2=M1+1
          IF(TIME.GT.TGWSER(M2,NS))THEN
            M1=M2
            GOTO 700
          ELSE
            MGWTLAST(NS)=M1
          ENDIF
C
          TDIFF=TGWSER(M2,NS)-TGWSER(M1,NS)
          WTM1=(TGWSER(M2,NS)-TIME)/TDIFF
          WTM2=(TIME-TGWSER(M1,NS))/TDIFF
          GWSERT(NS)=WTM1*GWSER(M1,NS)+WTM2*GWSER(M2,NS)
          DO NC=1,NCTMP
            GWCSERT(NC,NS)=WTM1*GWCSER(M1,NC,NS)+WTM2*GWCSER(M2,NC,NS)
          END DO
C
        ENDDO
C
        DO L=2,LA
          QGW(L)=GWFAC(L)*GWSERT(NGWSL(L))
        END DO
C
        DO NC=1,NCTMP
          DO L=2,LA
            CONGW(L,NC)=GWCSERT(NC,NGWSL(L))
          END DO
        END DO
C
      ENDIF
C
C**********************************************************************C
C
C **  CONTROL STRUCTURES AND TIDAL INLETS
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
      DO NCTL=1,NQCTL
        IF(NQCTYP(NCTL).LE.1)THEN
          NCTLT=NQCTLQ(NCTL)
          RQDW=1.
          IU=IQCTLU(NCTL)
          JU=JQCTLU(NCTL)
          LU=LIJ(IU,JU)
          HUP=HP(LU)+BELV(LU)+HCTLUA(NCTLT)
	    IF(NQCTYP(NCTL).EQ.-1)HUP=HP(LU)
          ID=IQCTLD(NCTL)
          JD=JQCTLD(NCTL)
          IF(ID.EQ.0.AND.JD.EQ.0)THEN
            LD=LC
            HDW=0.
            RQDW=0.
          ELSE
            LD=LIJ(ID,JD)
            HDW=HP(LD)+BELV(LD)+HCTLDA(NCTLT)
	      IF(NQCTYP(NCTL).EQ.-1)THEN
	        LD=LC
              HDW=0.
              RQDW=0.
	      ENDIF
          ENDIF
      DELH=HCTLUM(NCTLT)*HUP-HCTLDM(NCTLT)*HDW
CJMH INSERT SET POINT
      IF(NQCTYP(NCTL).EQ.0.AND.AQCTL(NCTLT).GT.0.0)THEN
        IF(HUP.LT.AQCTL(NCTLT)) DELH=-100.
      ENDIF
CJMH INSERT SET POINT
CJMH LIMIT FLOW FROM DRY CELL
      IF(DELH.LE.0.OR.HP(LU).LT.HWET)THEN
        DO K=1,KC
        QCTLT(K,NCTL)=0.
        ENDDO
       ELSE
        IF(NQCTYP(NCTL).EQ.1)DELH=SQRT(DELH)
        M1=0
        M2=1
 500    M1=M1+1
        M2=M2+1
        IF(M2.GT.MQCTL(NCTLT))THEN
           WRITE(6,6666)
           WRITE(6,6667)NCTL,NCTLT,IU,JU,ID,JD
           WRITE(6,6668)HUP,HP(LU),HDW,HP(LD)
!           WRITE(8,6666)   !hnr 7/27/2009
!           WRITE(8,6667)NCTL,NCTLT,IU,JU,ID,JD   !hnr 7/27/2009
!           WRITE(8,6668)HUP,HP(LU),HDW,HP(LD)   !hnr 7/27/2009
           STOP
        ENDIF
        IF(DELH.GE.HDIFCTL(M1,NCTLT).AND.DELH.LE.HDIFCTL(M2,NCTLT))THEN
          TDIFF=HDIFCTL(M2,NCTLT)-HDIFCTL(M1,NCTLT)
          WTM1=(HDIFCTL(M2,NCTLT)-DELH)/TDIFF
          WTM2=(DELH-HDIFCTL(M1,NCTLT))/TDIFF
           DO K=1,KC
           QCTLT(K,NCTL)=WTM1*QCTL(M1,1,K,NCTLT)
     &                  +WTM2*QCTL(M2,1,K,NCTLT)
           ENDDO
         ELSE
          GOTO 500
        ENDIF
      ENDIF
      IF(NQCTYP(NCTL).EQ.1)THEN
      IF(ISTL.EQ.3)THEN
        DO K=1,KC
         QCTLST(K,NCTL)=QCTLT(K,NCTL)
         TMPVAL=QCTLTO(K,NCTL)
     &       +DT*AQCTL(NCTLT)*QCTLST(K,NCTL)*QCTLST(K,NCTL)
         QCTLT(K,NCTL)=TMPVAL/(1.+DT*AQCTL(NCTLT)*QCTLTO(K,NCTL))
         QCTLTO(K,NCTL)=QCTLT(K,NCTL)
         QCTLSTO(K,NCTL)=QCTLST(K,NCTL)
        ENDDO
       ELSE
        DO K=1,KC
         QCTLST(K,NCTL)=QCTLT(K,NCTL)
         TMPVAL=QCTLTO(K,NCTL)
     &       +DT*AQCTL(NCTLT)*QCTLST(K,NCTL)*QCTLST(K,NCTL)
         QCTLT(K,NCTL)=TMPVAL/(1.+DT*AQCTL(NCTLT)*QCTLTO(K,NCTL))
         QCTLT(K,NCTL)=0.5*(QCTLT(K,NCTL)+QCTLTO(K,NCTL))
        ENDDO
      ENDIF
      ENDIF
	QCTLMAX=(HP(LU)-HDRY)*DXYP(L)/(DELT*FLOAT(KC))
      DO K=1,KC
       QCTLT(K,NCTL)=MIN(QCTLT(K,NCTL),QCTLMAX)
      ENDDO
      DO K=1,KC
	  IF(LGVCP(LU,K))
     &    QSUM(LU,K)=QSUM(LU,K)-RQCMUL(NCTL)*QCTLT(K,NCTL)
	  IF(LGVCP(LD,K))
     &    QSUM(LD,K)=QSUM(LD,K)+RQCMUL(NCTL)*RQDW*QCTLT(K,NCTL)
      ENDDO
      ENDIF
      ENDDO
C
      DO NCTL=1,NQCTL
      IF(NQCTYP(NCTL).EQ.2)THEN
      NCTLT=NQCTLQ(NCTL)
      RQDW=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      HUP=HP(LU)+BELV(LU)+HCTLUA(NCTLT)
CJMH LIMIT FLOW FROM DRY CELL
      IF(HUP.LT.HDIFCTL(1,NCTLT).OR.HP(LU).LT.HWET)THEN
        DO K=1,KC
         QCTLT(K,NCTL)=0.
        ENDDO
        GOTO 560
      ENDIF
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      LD=LIJ(ID,JD)
      HDW=HP(LD)+BELV(LD)+HCTLDA(NCTLT)
      HTMPD=HDIFCTD(1,NCTLT)+0.001
      HDW=MAX(HDW,HTMPD)
      MU1=0
      MU2=1
      MD1=0
      MD2=1
 555  MU1=MU1+1
      MU2=MU1+1
      IF(MU2.GT.MQCTL(NCTLT))THEN
        WRITE(6,6676)
        WRITE(6,6677)NCTL,NCTLT,IU,JU,ID,JD
        WRITE(6,6678)HUP,HP(LU),HDW,HP(LD)
        WRITE(6,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),
     &               HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
!        WRITE(8,6676)   !hnr 7/27/2009
!        WRITE(8,6677)NCTL,NCTLT,IU,JU,ID,JD   !hnr 7/27/2009
!        WRITE(8,6678)HUP,HP(LU),HDW,HP(LD)   !hnr 7/27/2009
!        WRITE(8,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),   !hnr 7/27/2009
!     &               HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)   !hnr 7/27/2009
        STOP
      ENDIF
      IF(HUP.GE.HDIFCTL(MU1,NCTLT).AND.HUP.LE.HDIFCTL(MU2,NCTLT))THEN
        TDIFFU=HDIFCTL(MU2,NCTLT)-HDIFCTL(MU1,NCTLT)
        WTM1U=(HDIFCTL(MU2,NCTLT)-HUP)/TDIFFU
        WTM2U=(HUP-HDIFCTL(MU1,NCTLT))/TDIFFU
       ELSE
        GOTO 555
      ENDIF
 556  MD1=MD1+1
      MD2=MD1+1
      IF(MD2.GT.MQCTL(NCTLT))THEN
        WRITE(6,6686)
        WRITE(6,6687)NCTL,NCTLT,IU,JU,ID,JD
        WRITE(6,6688)HUP,HP(LU),HDW,HP(LD)
        WRITE(6,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),
     &               HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
!        WRITE(8,6686)   !hnr 7/27/2009
!        WRITE(8,6687)NCTL,NCTLT,IU,JU,ID,JD   !hnr 7/27/2009
!        WRITE(8,6688)HUP,HP(LU),HDW,HP(LD)   !hnr 7/27/2009
!        WRITE(8,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),   !hnr 7/27/2009
!     &               HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)   !hnr 7/27/2009
        STOP
      ENDIF
      IF(HDW.GE.HDIFCTD(MD1,NCTLT).AND.HDW.LE.HDIFCTD(MD2,NCTLT))THEN
        TDIFFD=HDIFCTD(MD2,NCTLT)-HDIFCTD(MD1,NCTLT)
        WTM1D=(HDIFCTD(MD2,NCTLT)-HDW)/TDIFFD
        WTM2D=(HDW-HDIFCTD(MD1,NCTLT))/TDIFFD
       ELSE
        GOTO 556
      ENDIF
      DO K=1,KC
      QCTLT(K,NCTL)=WTM1U*( WTM1D*QCTL(MU1,MD1,K,NCTLT)
     &                    +WTM2D*QCTL(MU1,MD2,K,NCTLT) )
     &             +WTM2U*( WTM1D*QCTL(MU2,MD1,K,NCTLT)
     &                    +WTM2D*QCTL(MU2,MD2,K,NCTLT) )
      ENDDO
  560 CONTINUE
  	QCTLMAX=(HP(LU)-HDRY)*DXYP(L)/(DELT*FLOAT(KC))
      DO K=1,KC
       QCTLT(K,NCTL)=MIN(QCTLT(K,NCTL),QCTLMAX)
      ENDDO
      DO K=1,KC
	  IF(LGVCP(LU,K))
     &    QSUM(LU,K)=QSUM(LU,K)-RQCMUL(NCTL)*QCTLT(K,NCTL)
	  IF(LGVCP(LD,K))
     &    QSUM(LD,K)=QSUM(LD,K)+RQCMUL(NCTL)*RQDW*QCTLT(K,NCTL)
      ENDDO
      ENDIF
      ENDDO
C
      IF(ISCRAY.EQ.0)THEN
        TQCTL=TQCTL+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TQCTL=TQCTL+T2TMP-T1TMP
        WTQCLT=WTQCLT+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  FLOW WITHDRAWAL AND RETURN
C
      NTMP=4+NSED+NSND+NTOX
C
      DO NC=1,NTMP
        CQWRSERT(0,NC)=0.
      ENDDO
C
      DO NS=1,NQWRSR
C
       IF(ISTL.EQ.2)THEN
         IF(ISDYNSTP.EQ.0)THEN
           TIME=DT*(FLOAT(N)-0.5)/TCQWRSR(NS)
     &            +TBEGIN*(TCON/TCQWRSR(NS))
         ELSE
           TIME=TIMESEC/TCQWRSR(NS)
         ENDIF
       ELSE
         IF(ISDYNSTP.EQ.0)THEN
           TIME=DT*FLOAT(N-1)/TCQWRSR(NS)+TBEGIN*(TCON/TCQWRSR(NS))
         ELSE
           TIME=TIMESEC/TCQWRSR(NS)
         ENDIF
       ENDIF
C
       M1=MQWRTLST(NS)
  200  CONTINUE
       M2=M1+1
       IF(TIME.GT.TQWRSER(M2,NS))THEN
         M1=M2
         GOTO 200
        ELSE
         MQWRTLST(NS)=M1
       ENDIF
C
       TDIFF=TQWRSER(M2,NS)-TQWRSER(M1,NS)
       WTM1=(TQWRSER(M2,NS)-TIME)/TDIFF
       WTM2=(TIME-TQWRSER(M1,NS))/TDIFF
       QWRSERT(NS)=WTM1*QWRSER(M1,NS)+WTM2*QWRSER(M2,NS)
       DO NC=1,NTMP
        CQWRSERT(NS,NC)=WTM1*CQWRSER(M1,NS,NC)+WTM2*CQWRSER(M2,NS,NC)
       ENDDO
C
      ENDDO
C
      IF(NQWR.GT.0)THEN
      DO NWR=1,NQWR
      IU=IQWRU(NWR)
      JU=JQWRU(NWR)
      KU=KQWRU(NWR)
      ID=IQWRD(NWR)
      JD=JQWRD(NWR)
      KD=KQWRD(NWR)
      LU=LIJ(IU,JU)
      LD=LIJ(ID,JD)
      NS=NQWRSERQ(NWR)
	IF(LGVCP(LU,KU))
     &  QSUM(LU,KU)=QSUM(LU,KU)-QWR(NWR)-QWRSERT(NS)
	IF(LGVCP(LD,KD))
     &  QSUM(LD,KD)=QSUM(LD,KD)+QWR(NWR)+QWRSERT(NS)
      ENDDO
	ENDIF
C
C**********************************************************************C
C
C **  CALL JPEFDC AND PLACE JET-PLUME VOLUMES SOURCES
C
      IF(NQJPIJ.GT.0.AND.N.EQ.1) CALL JPEFDC
C
      IF(NQJPIJ.GT.0.AND.ISTL.EQ.3)THEN
        IF(NUDJPC(1).EQ.NUDJP(1))THEN
          CALL JPEFDC
          NUDJPC(1)=1
         ELSE
          NUDJPC(1)=NUDJPC(1)+1
        ENDIF
      ENDIF
C
      IF(NQJPIJ.GT.0.AND.IS2TIM.GE.1)THEN
        IF(NUDJPC(1).EQ.NUDJP(1))THEN
          CALL JPEFDC
          NUDJPC(1)=1
         ELSE
          NUDJPC(1)=NUDJPC(1)+1
        ENDIF
      ENDIF
C
C **  PLACE JET-PLUME VOLUMES SOURCES
C
      IF(NQJPIJ.GT.0)THEN
      DO NJP=1,NQJPIJ
C
      IF(ICALJP(NJP).EQ.1)THEN
	 RPORTS=FLOAT(NPORTJP(NJP))
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
C QVJPTMP = JETPLUME DISCHARGE PER PORT
       QVJPTMP=QQCJP(NJP)
       DO K=1,KC
        QVJPTMP=QVJPTMP+QSERT(K,NQSERJP(NJP))
       ENDDO
C SUBTRACT THE ENTRAINMENT FROM EACH LAYER
	 DO K=1,KC
	   IF(LGVCP(LJP,K))
     &     QSUM(LJP,K)=QSUM(LJP,K)-RPORTS*QJPENT(K,NJP)
	 ENDDO
C PLACE DISCHARGE AND TOTAL ENTRAINMENT AT EFFECTIVE LOCATION
	 IF(LGVCP(LJP,KTMP))
     &   QSUM(LJP,KTMP)=QSUM(LJP,KTMP)+RPORTS*(QVJPTMP+QJPENTT(NJP))
      ENDIF
C
      IF(ICALJP(NJP).EQ.2)THEN
	 RPORTS=FLOAT(NPORTJP(NJP))
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
C QVJPTMP = JETPLUME DISCHARGE PER PORT
       QVJPTMP=QWRCJP(NJP)+QWRSERT(NQWRSERJP(NJP))
C SUBTRACT ENTRAIMENT FROM EACH LAYER
       DO K=1,KC
	   IF(LGVCP(LJP,K))
     &     QSUM(LJP,K)=QSUM(LJP,K)-RPORTS*QJPENT(K,NJP)
	 ENDDO
C PLACE DISCHARGE AND TOTAL ENTRAINMENT AT EFFECTIVE LOCATION
	 IF(LGVCP(LJP,KTMP))
     &   QSUM(LJP,KTMP)=QSUM(LJP,KTMP)+RPORTS*(QVJPTMP+QJPENTT(NJP))
C REMOVE DISCHARGE FROM UPSTREAM INTAKE CELL
       LU=LIJ(IUPCJP(NJP),JUPCJP(NJP))
       KU=KUPCJP(NJP)
	 IF(LGVCP(LU,KU))
     &   QSUM(LU,KU)=QSUM(L,KU)-RPORTS*QVJPTMP
      ENDIF
C
      ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  GROUND WATER INTERACTION, EVAPORATION AND RAINFALL
C
      IF(ISGWIE.EQ.0)THEN
        IF(EVAPCVT.LT.0.)THEN
          DO L=2,LA
           SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/
     &              (1.+0.00412*TEM(L,KC))))
          EVAPT(L)=CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW-VPA(L))/PATMT(L)
		IF(HP(L).LT.HWET) EVAPT(L)=0.
          IF(LGVCP(L,KC))
     &      QSUM(L,KC)=QSUM(L,KC)+DXYP(L)*(RAINT(L)-EVAPT(L))
          ENDDO
         ELSE
          DO L=2,LA
		 IF(HP(L).LT.HWET) EVAPT(L)=0.
           IF(LGVCP(L,KC))
     &       QSUM(L,KC)=QSUM(L,KC)+DXYP(L)*(RAINT(L)-EVAPT(L))
          ENDDO
        ENDIF
       ELSE
        DO L=2,LA
        	 IF(LGVCP(L,KC))
     &   QSUM(L,KC)=QSUM(L,KC)+DXYP(L)*RAINT(L)
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  LAKE CORRECTING FLOW ADDED                 !hnr
C                                                !hnr
c      IF(ISLAKE.GE.1)THEN                        !hnr
c        LK=LIJ(ILK,JLK)                          !hnr
C        DO K=1,KC                                !hnr
C          QSUM(LK,K)=QSUM(LK,K)+QLK(K)           !hnr
C        end do                                   !hnr
C      ENDIF                                      !hnr
C                                                !hnr
C**********************************************************************C
C
C **  DETERMINE NET EXTERNAL VOLUME SOURCE/SINK
C
      DO K=1,KC
       DO L=1,LC
       QSUME(L)=QSUME(L)+QSUM(L,K)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  UPDATE ZERO DIMENSION VOLUME BALANCE
C
C----------------------------------------------------------------------C
C
C     IF(ISDRY.GE.1.AND.ISTL.EQ.3)THEN
C       VOLADD=0.
C       DO L=2,LA
C       IF(SPB(L).NE.0)THEN
C         VOLADD=VOLADD+QSUME(L)
C       ENDIF
C       ENDDO
C       VOLADD=VOLADD*DT
C       VOLZERD=VOLZERD+VOLADD
C       VETZERD=VETZERD+VOLADD+DT*QETTMP
C     ENDIF
C
C5303 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
C
C**********************************************************************C
C
C **  WRITE DIAGNOSTIC FILE FOR VOLUME SOURCES,SINKS, ETC
C
      ITMPD=0
      IF(ISDIQ.EQ.2.AND.ISTL.EQ.2) ITMPD=1
      IF(ISDIQ.EQ.1) ITMPD=1
C
      NTT=4+NTOX+NSED+NSND
C
      IF(ITMPD.EQ.1)THEN
       IF(N.EQ.NTSPTC.OR.N.EQ.1)THEN
         OPEN(1,FILE='QDIAG.OUT',STATUS='UNKNOWN')
         CLOSE(1,STATUS='DELETE')
         OPEN(1,FILE='QDIAG1.OUT',STATUS='UNKNOWN')
         CLOSE(1,STATUS='DELETE')
         OPEN(1,FILE='QDIAG1.OUT',STATUS='UNKNOWN')
        ELSE
         OPEN(1,FILE='QDIAG.OUT',POSITION='APPEND',STATUS='UNKNOWN')
       ENDIF
       WRITE(1,101)N
C
        DO LL=1,NQSIJ
        NQSTMP=NQSERQ(LL)
        NCSTMP=NCSERQ(LL,1)
        L=LQS(LL)
        I=IL(L)
        J=JL(L)
        WRITE(1,102)I,J
        WRITE(1,216)LL,L,(QSS(K,LL),K=1,KC)
        DO NT=1,NTT
          WRITE(1,217)LL,NT,(CQS(K,LL,NT),K=1,KC)
        ENDDO
        WRITE(1,104)
        WRITE(1,105)I,J
        WRITE(1,206)LL,L,(QSERT(K,NQSTMP),K=1,KC)
        DO NT=1,NTT
         NCSTMP=NCSERQ(LL,NT)
         WRITE(1,207)LL,NT,NCSTMP,(CSERT(K,NCSTMP,NT),K=1,KC)
         ENDDO
        WRITE(1,104)
        ENDDO
C
        DO NCTL=1,NQCTL
        IU=IQCTLU(NCTL)
        JU=JQCTLU(NCTL)
        ID=IQCTLD(NCTL)
        JD=JQCTLD(NCTL)
        IF(IU.EQ.0.AND.JU.EQ.0)THEN
          LU=0
          HUP=0.
         ELSE
          LU=LIJ(IU,JU)
          HUP=HP(LU)+BELV(LU)+HCTLUA(NCTLT)
          IF(NQCTYP(NCTL).EQ.-1)HUP=HP(LU)
        ENDIF
        IF(ID.EQ.0.AND.JD.EQ.0)THEN
          LD=0
          HDW=0.
         ELSE
          LD=LIJ(ID,JD)
          HDW=HP(LD)+BELV(LD)+HCTLDA(NCTLT)
        ENDIF
        WRITE(1,107)IU,JU,LU,NCTLT,HUP
         DO K=1,KC
         WRITE(1,108)K,QCTLT(K,NCTL)
         ENDDO
        WRITE(1,104)
        WRITE(1,109)ID,JD,LD,NCTLT,HDW
         DO K=1,KC
         WRITE(1,108)K,QCTLT(K,NCTL)
         ENDDO
        WRITE(1,104)
        ENDDO
C
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
        WRITE(1,110)IU,JU
        WRITE(1,111)KU,QWR(NWR),CQWR(NWR,1),CQWR(NWR,2)
        WRITE(1,104)
        WRITE(1,112)ID,JD
        WRITE(1,111)KD,QWR(NWR),CQWR(NWR,1),CQWR(NWR,2)
        WRITE(1,104)
        WRITE(1,113)IU,JU
        WRITE(1,114)KU,QWRSERT(NQSTMP),CQWRSERT(NQSTMP,1),
     &                 CQWRSERT(NQSTMP,2)
        WRITE(1,104)
        WRITE(1,115)ID,JD
        WRITE(1,114)KD,QWRSERT(NQSTMP),CQWRSERT(NQSTMP,1),
     &                 CQWRSERT(NQSTMP,2)
        WRITE(1,104)
        ENDDO
C
        CLOSE(1)
      ENDIF
C
  101 FORMAT('  SOURCE/SINK DIAGNOSTICS AT TIME STEP =',I8,//)
  102 FORMAT(3X,'CONST NQSIJ SOURCE/SINK FLOW AT I =',I5,' J =',I5,/)
  103 FORMAT(5X,'K =',I5,5X,'QSS(K) = ',E12.4,5X,'CQS(K,1) = ',E12.4,
     &       5X,'CQS(K,5) = ',E12.4)
  203 FORMAT(5X,'K =',I5,5X,'QSS(K) = ',E12.4,5X,'CQS(K, ) = ',
     &       5X, 12E12.4)
  104 FORMAT(/)
  105 FORMAT(3X,'TIME VAR NQSIJ SOURCE/SINK FLOW AT I =',I5,' J=',I5,/)
  106 FORMAT(5X,'K =',I5,5X,'QSERT(K) = ',E12.4,
     &       5X,'CSERT(K,1) = ',E12.4,5X,'CSERT(K,5) = ',E12.4)
  206 FORMAT(5X,'NQ,LQ     =',2I4,7X,'QSERT() = ',12E12.4)
  207 FORMAT(5X,'NQ,NT,NCQ =',3I4,3X,'CSERT() = ',12E12.4)
  216 FORMAT(5X,'NQ,LQ =',2I4,3X,'QSS() = ',12E12.4)
  217 FORMAT(5X,'NQ,NT =',2I4,3X,'CQS() = ',12E12.4)
  107 FORMAT(3X,'UPSTRM CONTROLED SINK FLOW AT I =',I5,' J =',I5,
     &          ' L =',I5,'  NQCTLT =',I5,'  HUP = ',E12.4/)
  108 FORMAT(5X,'K =',I5,5X,'QCTL(K) = ',2E12.4)
  109 FORMAT(3X,'DWNSTRM CONTROLED SOURCE FLOW AT I =',I5,' J =',I5,
     &          ' L =',I5,'  NQCTLT =',I5,'  HDW = ',E12.4/)
  110 FORMAT(3X,'UPSTRM CONST WITHDRW SINK FLOW AT I =',I5,' J =',I5,/)
  111 FORMAT(5X,'K =',I5,5X,'QWR(K) = ',E12.4,
     &       5X,'CQWR(1) = ',E12.4,5X,'CQWR(2) = ',E12.4)
  112 FORMAT(3X,'DWNSTRM CONST RETN SOURCE FLOW AT I =',I5,' J =',I5,/)
  113 FORMAT(3X,'UPSTRM VAR WITHDRW SINK FLOW AT I =',I5,' J =',I5,/)
  114 FORMAT(5X,'K =',I5,5X,'QSERT(K) = ',E12.4,
     &       5X,'CSERT(K,1) = ',E12.4,5X,'CSERT(K,2) = ',E12.4)
  115 FORMAT(3X,'DWNSTRM VAR RETN SOURCE FLOW AT I =',I5,' J =',I5,/)
 6666 FORMAT(' SINGLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS ')
 6667 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
 6668 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
 6676 FORMAT(' DOUBLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS, UP ')
 6677 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
 6678 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
 6679 FORMAT(' HUF,HUL,HDF,HDL = ',4(2X,E12.4))
 6686 FORMAT(' DOUBLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS, DW ')
 6687 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
 6688 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
 6111 FORMAT(' INVALID NQSIJ LOCATION, NQSIJ,I,J = ',7I5)
C
C**********************************************************************C
C
      RETURN
      END
