C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SSEDTOXGVC(ISTLX,IS2TLX,CORDTX)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
C **  SUBROUTINE SSEDTOXGVC SHIFT VARIABLES INTO POSITION TO USE
C **  SIGMA VERSION OF SSEDTOX AND DAUGHTER SUBROUTINES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      COMMON/SSEDTOX1A/ CBEDTOTAL(LCM),QCELLCTR(LCM),HGDH(LCM),
     &                  FRACCOH(LCM,KBM),FRACNON(LCM,KBM),USTARSED(LCM),
     &                  USTARSND(LCM),QWATPA(LCM),QSSDPA(LCM)
C
C**********************************************************************C
C
C     DIAGNOSTICS
C
      NTSTBCM=NTSTBC-1
      IF(N.EQ.NTSTBCM)THEN
	  OPEN(1,FILE='SSEDTOXGVC1.OUT')
	  CLOSE(1,STATUS='DELETE')
	  OPEN(1,FILE='SSEDTOXGVC1.OUT')
	  DO L=2,LA
	    WRITE(1,101)IL(L),JL(L),KGVCP(L),(SED(L,K,1),K=1,KC)
	  ENDDO
	  CLOSE(1)
	ENDIF
C
C**********************************************************************C
C
C     SHIFT DOWN
C
      DO L=2,LA
 	  KTMP=KC+1-KGVCP(L)
	  KTMPP=KTMP+1
	  KO=KGVCP(L)-1
	  IF(KGVCP(L).GT.1)THEN
c
	  KBPM=KGVCP(L)-1
        QQ(L,0)=QQ(L,KBPM)
C
	  DO K=1,KTMP
          U(L,K)=U(L,K+KO)
          U(L,K)=V(L,K+KO)
          SAL(L,K)=SAL(L,K+KO)
          TEM(L,K)=TEM(L,K+KO)
          DYE(L,K)=DYE(L,K+KO)
          SFL(L,K)=SFL(L,K+KO)
          SEDFLOCDIA(L,K)=SEDFLOCDIA(L,K+KO)
		FLOCDIA(L,K)=FLOCDIA(L,K+KO)
          WSETFLOC(L,K)=WSETFLOC(L,K+KO)
          STPOCW(L,K)=STPOCW(L,K+KO)
          STDOCW(L,K)=STDOCW(L,K+KO)
          DO NS=1,NSED
            SED(L,K,NS)=SED(L,K+KO,NS)
          ENDDO
          DO NS=1,NSND
            SND(L,K,NS)=SND(L,K+KO,NS)
          ENDDO
          DO NS=1,NSED+NSND
            STFPOCW(L,K,NS)=STFPOCW(L,K+KO,NS)
          ENDDO
          DO NT=1,NTOX
            TOX(L,K,NT)=TOX(L,K+KO,NT)
          ENDDO
	  ENDDO
c
        IF(KTMPP.LE.KC)THEN
	    DO K=KTMPP,KC
            U(L,K)=0.0
            V(L,K)=0.0
            SAL(L,K)=0.0
            TEM(L,K)=0.0
            DYE(L,K)=0.0
            SFL(L,K)=0.0
            SEDFLOCDIA(L,K)=0.0
            FLOCDIA(L,K)=0.0
            WSETFLOC(L,K)=0.0
            STPOCW(L,K)=0.0
            STDOCW(L,K)=0.0
            DO NS=1,NSED
              SED(L,K,NS)=0.0
            ENDDO
            DO NS=1,NSND
              SND(L,K,NS)=0.0
            ENDDO
            DO NS=1,NSED+NSND
              STFPOCW(L,K,NS)=0
            ENDDO
            DO NT=1,NTOX
              TOX(L,K,NT)=0.0
            ENDDO
          ENDDO
	  ENDIF
c
        ENDIF
      ENDDO
C
      IF(N.EQ.NTSTBCM)THEN
	  OPEN(1,FILE='SSEDTOXGVC2.OUT')
	  CLOSE(1,STATUS='DELETE')
	  OPEN(1,FILE='SSEDTOXGVC2.OUT')
	  DO L=2,LA
	    WRITE(1,101)IL(L),JL(L),KGVCP(L),(SED(L,K,1),K=1,KC)
	  ENDDO
	  CLOSE(1)
	ENDIF
C
      IF(N.EQ.NTSTBCM)THEN
	  DO L=2,LA
          TAUB(L)=QQ(L,0)/CTURB2
	    TAUBSED(L)=TAUB(L)
	    TAUBSND(L)=TAUB(L)
	  ENDDO
      ENDIF

C**********************************************************************C
C
C     CALL SSEDTOX
C
	IF(ISHOUSATONIC.EQ.0)THEN
        CALL SSEDTOX(ISTLX,IS2TLX,CORDTX)
      ELSE
        CALL SSEDTOXHOUS(ISTLX,IS2TLX,CORDTX)
      ENDIF
C
C
      IF(N.EQ.NTSTBCM)THEN
	  OPEN(1,FILE='SSEDTOXGVC3.OUT')
	  CLOSE(1,STATUS='DELETE')
	  OPEN(1,FILE='SSEDTOXGVC3.OUT')
	  DO L=2,LA
	    WRITE(1,101)IL(L),JL(L),KGVCP(L),(SED(L,K,1),K=1,KC)
	  ENDDO
	  CLOSE(1)
	ENDIF
C
      IF(N.EQ.NTSTBCM)THEN
	  OPEN(1,FILE='SSEDTOXGVC0.OUT')
	  CLOSE(1,STATUS='DELETE')
	  OPEN(1,FILE='SSEDTOXGVC0.OUT')
	  DO L=2,LA
	    WRITE(1,102)IL(L),JL(L),KGVCP(L),QCELLCTR(L),CBEDTOTAL(L),
     &                                TAUB(L),TAUBSED(L),TAUBSND(L)
c	    WRITE(1,103)IL(L),JL(L),KGVCP(L),RSSBCW(L),WCORWST(L),
c     &                                     RSSBCE(L),WCOREST(L),
c     &                                     RSSBCS(L),WCORSTH(L),
c     &                                     RSSBCN(L),WCORNTH(L)
	  ENDDO
	  CLOSE(1)
	ENDIF
C
C**********************************************************************C
C
C     SHIFT UP
C
      NSORB=NSED+NSND+2
C
      DO L=2,LA
 	  KTMP=KC+1-KGVCP(L)
	  KTMPP=KTMP+1
	  KO=KGVCP(L)-1
	  IF(KGVCP(L).GT.1)THEN
C
        DO K=KC,KGVCP(L),-1
	    U(L,K)=U(L,K-KO)
	    V(L,K)=V(L,K-KO)
	    SAL(L,K)=SAL(L,K-KO)
	    TEM(L,K)=TEM(L,K-KO)
	    DYE(L,K)=DYE(L,K-KO)
	    SFL(L,K)=SFL(L,K-KO)
          SEDFLOCDIA(L,K)=SEDFLOCDIA(L,K-KO)
		FLOCDIA(L,K)=FLOCDIA(L,K-KO)
          WSETFLOC(L,K)=WSETFLOC(L,K-KO)
          STPOCW(L,K)=STPOCW(L,K-KO)
          STDOCW(L,K)=STDOCW(L,K-KO)
          DO NS=1,NSED
	      SED(L,K,NS)=SED(L,K-KO,NS)
          ENDDO
          DO NS=1,NSND
	      SND(L,K,NS)=SND(L,K-KO,NS)
          ENDDO
          DO NS=1,NSED+NSND
	      STFPOCW(L,K,NS)=STFPOCW(L,K-KO,NS)
          ENDDO
          DO NT=1,NTOX
	      TOX(L,K,NT)=TOX(L,K-KO,NT)
	      TOXFDFW(L,K,NT)=TOXFDFW(L,K-KO,NT)
	      TOXCDFW(L,K,NT)=TOXCDFW(L,K-KO,NT)
	      TOXPFTW(L,K,NT)=TOXPFTW(L,K-KO,NT)
          ENDDO
          DO NT=1,NTOX
	      DO NS=1,NSORB
	        TOXPFW(L,K,NS,NT)=TOXPFW(L,K-KO,NS,NT)
	      ENDDO
          ENDDO
        ENDDO
C
        DO K=1,KGVCP(L)-1
          U(L,K)=0.0
          V(L,K)=0.0
          SAL(L,K)=0.0
          TEM(L,K)=0.0
          DYE(L,K)=0.0
          SFL(L,K)=0.0
          SEDFLOCDIA(L,K)=0.0
		FLOCDIA(L,K)=0.0
          WSETFLOC(L,K)=0.0
          STPOCW(L,K)=0.0
          STDOCW(L,K)=0.0
          DO NS=1,NSED
	      SED(L,K,NS)=0.0
          ENDDO
          DO NS=1,NSND
	      SND(L,K,NS)=0.0
          ENDDO
          DO NS=1,NSED+NSND
	      STFPOCW(L,K,NS)=0.0
          ENDDO
          DO NT=1,NTOX
	      TOX(L,K,NT)=0.0
	      TOXFDFW(L,K,NT)=0.0
	      TOXCDFW(L,K,NT)=0.0
	      TOXPFTW(L,K,NT)=0.0
          ENDDO
          DO NT=1,NTOX
	      DO NS=1,NSORB
	        TOXPFW(L,K,NS,NT)=0.0
	      ENDDO
          ENDDO
        ENDDO
C
        ENDIF
      ENDDO
C
      IF(N.EQ.NTSTBCM)THEN
	  OPEN(1,FILE='SSEDTOXGVC4.OUT')
	  CLOSE(1,STATUS='DELETE')
	  OPEN(1,FILE='SSEDTOXGVC4.OUT')
	  DO L=2,LA
	    WRITE(1,101)IL(L),JL(L),KGVCP(L),(SED(L,K,1),K=1,KC)
	  ENDDO
	  CLOSE(1)
	ENDIF
C
C**********************************************************************C
C
  101 FORMAT(3I5,25E10.2)
  102 FORMAT(3I5,8E14.5)
  103 FORMAT(3I5,8F8.2)
C
C**********************************************************************C
C
      RETURN
	END
