       SUBROUTINE BMD2_DUMP
C
C **  CREATED BY HUGO N RODRIGUEZ ON july 2, 2019
C
C **  SUBROUTINE BMD_DUMP2 WRITES FULL FIELD OF MODEL VARIABLES
C **  AT SPECIFIED TIME INTERVALS IN A BMD2 FILE FORMAT
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      INCLUDE 'ALLSET.INT'
C
C**********************************************************************C
	REAL*8 BMD_TIME
	REAL*4 xx
      REAL*4, ALLOCATABLE :: DATA_OUT(:,:)
      
C
      BMD_WRITE=BMD_WRITE+1
      IF(BMD_WRITE.GT.NumBMD2Writes) RETURN

      BMD_TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
C
      ALLOCATE(DATA_OUT(NUMSEG,NUMVAR))
      
      ISEG=0
	do kk=1, ibmdseg
        l=isegin(kk)
        DO K=1,KC
          IF(.NOT.(IGRIDV.EQ.1.AND.K.LT.KGVCP(L))) THEN
            ISEG=ISEG+1
            IVAR=0

            IF (jjf(1).eq.1) THEN
              IVAR=IVAR+1
              DATA_OUT(ISEG,IVAR)=P(l) * GI
            ENDIF

            IF (jjf(2).eq.1) THEN
              IVAR=IVAR+1
              DATA_OUT(ISEG,IVAR)=HP(L)
            ENDIF

            IF (jjf(3).eq.1) THEN
             IVAR=IVAR+1
             xx=DXP(L)*DYP(L)*HP(L)/(ABS(UHDYE(L))+ABS(VHDXE(L)))/3600.0
             IF(UHDYE(L).EQ.0.0.AND.VHDXE(L).EQ.0.0) xx=1.0E10
             DATA_OUT(ISEG,IVAR)=XX
            ENDIF

            IF (jjf(4).eq.1) THEN
              IVAR=IVAR+1
              cflmax=0.0
              do kki=1,kc
                IF(IGRIDV.EQ.1.AND.K.LT.KGVCP(L)) GOTO 10
                xx=max(cflmax,cfluuu(l,kki))
                xx=max(xx,cflvvv(l,kki))
                xx=max(xx,cflcac(l,kki))
                if (k.ne.kc) then
                  xx=max(xx,cflwww(l,kki))
                end if
10              CONTINUE                    
              end do
              DATA_OUT(ISEG,IVAR)=XX
            ENDIF

            IF (jjf(5).eq.1)  THEN
              IVAR=IVAR+1
              LN=LNC(l)
              Ubmd=500.*(TBX(L+1)+TBX(L))
              Vbmd=500.*(TBY(LN)+TBY(L))
              xx=(ubmd*ubmd+vbmd*vbmd)**0.5    
              DATA_OUT(ISEG,IVAR)=XX
            ENDIF
         
            IF (jjf(6).eq.1)  THEN
              IVAR=IVAR+1
		    IF(IGRIDV.EQ.0) then
                XX=hp(l)*DZC(k)                !LAYER DEPTH FOR SIGMA GRID
              ELSE
                XX=hp(l)*DZC(k)*gvcsclp(l)                !LAYER DEPTH FOR GVC
              END IF
              DATA_OUT(ISEG,IVAR)=XX
            ENDIF

            IF (jjf(7).eq.1) THEN
              IVAR=IVAR+1
              DATA_OUT(ISEG,IVAR)=SAL(L,K)
            ENDIF

            IF (jjf(8).eq.1) THEN
              IVAR=IVAR+1
              DATA_OUT(ISEG,IVAR)=TEM(L,K)
            ENDIF
              
            IF (jjf(9).eq.1) THEN
              IVAR=IVAR+1
              DATA_OUT(ISEG,IVAR)=DYE(L,K)
            ENDIF

            IF (jjf(10).eq.1) THEN
              IVAR=IVAR+1
              xx=sedt(l,k)+sndt(l,k)
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
c
c             velocities U and V at EAST AND SOUTH FACE
            IF(IGRIDV.EQ.0) then                  !for sigma grid
              RUVTMP=50.
              IF(SPB(L).EQ.0) RUVTMP=100.
	        LNBMD=LNC(l)
 	        ubmd= (u(l,k)+u(l+1,k))*RUVTMP    !VELOCITY U AT CENTER TO CALCULATE EAST AND NORTH VELOCITIES AT CENTER
	        vbmd= (v(l,k)+v(LNBMD,k))*RUVTMP   !VELOCITY v AT CENTER TO CALCULATE EAST AND NORTH VELOCITIES AT CENTER
            ELSE
c                    velocities U and V at cell center for GVC grid
              IF (LGVCU(L,K)) THEN
                IF (LGVCU(L+1,K)) THEN
                  ubmd= (u(l,k)+u(l+1,k))/2.*100.
                ELSE
                  ubmd= u(l,k)*100.
                END IF
              ELSE
                IF (LGVCU(L+1,K)) THEN
                  ubmd= u(l+1,k)*100.
                ELSE
                  ubmd= 0.0
                END IF
              END IF
              IF (LGVCV(L,K)) THEN
	          LNBMD=LNC(l)
                IF (LGVCV(LNBMD,K)) THEN
                  vbmd= (v(l,k)+v(LNBMD,k))/2.*100.
                ELSE
                  vbmd= v(l,k)*100.
                END IF
              ELSE
	          LNBMD=LNC(l)
                IF (LGVCV(LNBMD,K)) THEN
                  vbmd= v(LNBMD,k)*100.
                ELSE
                  vbmd= 0.0
                END IF
              END IF
            END IF
C
            IF (jjf(11).eq.1) THEN
              IVAR=IVAR+1
              xx=U(L,K)*100.0
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
               
            IF (jjf(12).eq.1) THEN
              IVAR=IVAR+1
              xx=V(L,K)*100.0
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF

c             velocities East and North at cell center
            IF (jjf(13).eq.1) THEN
              IVAR=IVAR+1
              xx=CUE(L)*ubmd+CVE(L)*vbmd
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF

            IF (jjf(14).eq.1) THEN
              IVAR=IVAR+1
              xx=CUN(L)*ubmd+CVN(L)*vbmd
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF

            IF (jjf(15).eq.1) THEN
              IVAR=IVAR+1
              xx=(ubmd*ubmd+vbmd*vbmd)**0.5
              DATA_OUT(ISEG,IVAR)=xx     !SPEED AT CELL CENTER
            ENDIF
C
c             Flow U (bewtween nodes I,J and I,J+1) and Flow V (bewtween nodes I,J and I+1,J)

            IF (jjf(16).eq.1) THEN
              IVAR=IVAR+1
		    IF(IGRIDV.EQ.0) THEN
 	          xx = u(l,k)*dyu(l)*hu(l)*dzc(k)
 	        ELSE
                XX=U(L,K)*DYU(L)*HU(L)*DZC(K)*GVCSCLU(L)
              END IF
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
C
            IF (jjf(17).eq.1) THEN
              IVAR=IVAR+1
		    IF(IGRIDV.EQ.0) THEN
	          xx= v(l,k)*dxv(l)*hv(l)*dzc(k)
	        ELSE
                xx=v(l,k)*DXV(L)*HV(L)*dzc(k)*gvcsclV(l)
                if (K.LT.kgvcp(l)) xx=-999.0
              END IF
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
C	        
            IF (jjf(18).eq.1) THEN
              IVAR=IVAR+1
              DATA_OUT(ISEG,IVAR)=CFLUUU(L,K)
            ENDIF

            IF (jjf(19).eq.1) THEN
              IVAR=IVAR+1
              DATA_OUT(ISEG,IVAR)=CFLVVV(L,K)
            ENDIF

            IF (jjf(20).eq.1) THEN
              IVAR=IVAR+1
              DATA_OUT(ISEG,IVAR)=CFLCAC(L,K)
            ENDIF
C
            IF (jjf(21).eq.1) THEN                       !volume of the cell at the layer level
              IVAR=IVAR+1
		    IF(IGRIDV.EQ.0) THEN
                xx=(DZC(K)*HP(L))*DXYP(L)
              ELSE
                XX=hp(l)*DZC(k)*gvcsclp(l)*DXYP(L)                !LAYER VOLUME MODIFIED FOR GVC
              END IF
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
c
            IF (jjf(22).eq.1) THEN
              IVAR=IVAR+1
              xx=w(l,k)*100.
    	        IF(K.EQ.KC) xx=-999.0
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
c
            IF (jjf(23).eq.1) THEN
              IVAR=IVAR+1
              xx=dxyp(l)*w(l,k)
    	        IF(K.EQ.KC) xx=-999.0
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
C
            IF (jjf(24).eq.1) THEN
              IVAR=IVAR+1
              xx=cflwww(L,k)
    	        IF(K.EQ.KC) xx=-999.0
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
C
            IF (jjf(25).eq.1) THEN
              IVAR=IVAR+1
 		    IF(IGRIDV.EQ.0) THEN
                xx=ab(l,k)*hp(l)*DZC(k)
              ELSE
                xx=ab(l,k)*hp(l)*DZC(k)*gvcsclp(l)
              END IF  
    	        IF(K.EQ.KC) xx=-999.0
              DATA_OUT(ISEG,IVAR)=xx
            ENDIF
          ENDIF
        end do
      end do

      Call bmd2setnextframe(HANDLEBMD2,BMD_TIME,DATA_OUT,iError)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
 
      
      return
	end
