      SUBROUTINE WASPHYDROLINKgvc
!
! DATE MODIFIED     BY                 DATE APPROVED    BY
!
! 05/08/2019        DRBC LZ - SET inumsegconsts VALUE FOR SALINITY
! 06/13/2019        DRBC LZ - ADDED TEMPORAL SALINITY OUTPUT TO HYDRO-LINKAGE FILE
! 06/24/2019        HNR_GHD - ADJUSTED VALUES FOR DRY CELLS WHEN W/D IS USED
! 06/28/2019        DRBC LZ - ADDED INITIAL SALINITY OUTPUT TO HYDRO-LINKAGE FILE DURING A HOT-START
! 07/23/2019        DRBC LZ - ADJUSTED THE LOGIC SEQUENCE IF A SOUTH CELL IS WITH ID 8
! 08/26/2020        DRBC LZ - MODIFIED MAX DISPERSION CHECK TO AVOID FEEDBACK TO OTHER EFDC SUBROUTINES
! 12/20/2020        DRBC LZ - ADDED TKE DISSIPATION RATE
! 04/21/2022        HNR_GHD - UNLINKED DISPRATE FROM SALINITY - NOT IT CAN BE SETUP WHEN SALINITY IS NOT SIMULATED
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      INCLUDE 'ALLset.INT'
      parameter(nf=44000,nsg=20000)
!
      DIMENSION LAUX(ICM,JCM,KCM)
      INTEGER LU,LD,IU,ID,JU,JD,NAUX,nlay,SUR
      INTEGER istartyear,istartmonth,istartday,istarthour,istartminute
     *,istartsecond,Ihl_debug,Ihl_mode,inumsegconsts,j,
     *FLAGWASPBC(NQSERM,KCM)
      CHARACTER*6  sn(lcm),sn1
      CHARACTER*20 HYDFIL
      CHARACTER*23 segname
      CHARACTER*17 SN2
      character*80 DES,MOD
      character*3 Itext, Jtext, Ktext

	REAL*8  AUX,tmin,thour,tsec
        REAL  rinterval,crnu(nf),brintt(nf),flow(nf)
        REAL  SegVolume(nsg),SegDepth(nsg),SegVel(nsg)
        REAL  SegSalt(nsg),SegTemp(nsg),vol1,vol2,ad,adcoeff,abwmax
        real  abwmx(nsg),SegDispr(nsg)
        REAL  ABTMP      !DRBC LZ, 8/26/2020
!
!
!**********************************************************************C
!
! **  READ CONTROL DATA FILE EFDC.WSP
!
!----------------------------------------------------------------------C
!
      SVPT=1.
      IF(NTSMMT.LT.NTSPTC)SVPT=0.
!
      IF(JSWASP.EQ.1) THEN

        OPEN(1,FILE='EFDC.WSP',STATUS='UNKNOWN')
        WRITE(6,*)'EFDC.WSP opened'

!
!1**  READ CELL VOLUME PARAMETERS
!
        READ(1,1)
        READ(1,1)
        READ(1,*) IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP
!
!2**  READ DIFFUSION PARAMETERS
!
        READ(1,1)
        READ(1,1)
        READ(1,*) NRFLD,SCALR,CONVR,ISNKH,adcoeff,abwmax

        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          abwmx(l)=abwmax
        END DO
!
!3**  READ ADVECTION PARAMETERS
!
        READ(1,1)
        READ(1,1)
        READ(1,*) IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD,IDAYS
!
!4**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)
!
        READ(1,1)
        READ(1,1)
        READ(1,*) DEPSED,TDINTS,SEDIFF, WSS1, WSS2, WSS3
!
        do i=1,5
          read(1,*,err=11)
        end do
11    continue
        DO LT=2,LALT
          read(1,*,err=12)i,j,abwmax
          l=lij(i,j)
          abwmx(l)=abwmax
        END DO
12     continue

        CLOSE(1)

        OPEN(1,FILE='ABmax.txt',STATUS='UNKNOWN')
        write(1,*)'    I    J     ABmax'
        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          write(1,21)illt(lt),jllt(lt),abwmx(l)
        END DO
        close(1)
21      format(2I5,f10.6)
        WRITE(6,*)'EFDC.WSP read succesfully and ABmax.txt written'
!
    1   FORMAT (80X)
!

        IF(NQSER.GE.1)THEN
          OPEN(1,FILE='QSER.INP',STATUS='UNKNOWN')
!
          DO IS=1,14
           READ(1,1)
          ENDDO
!
          DO NS=1,NQSER
            READ(1,*,IOSTAT=ISO)ISTYP, MQSER(NS),TCQSER(NS),TAQSER(NS),RMULADJ,ADDADJ,ICHGQS
            IF(ISO.GT.0) GOTO 860
            IF(ISTYP.EQ.1)THEN
              READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
              do k=1,kc
                if(wkq(k).lt.1.e-6) flagwaspbC(ns,k)=1
              end do
              IF(ISO.GT.0) GOTO 860
              DO M=1,MQSER(NS)
                READ(1,*,IOSTAT=ISO)TQSER(M,NS),QSERTMP
                IF(ISO.GT.0) GOTO 860
              END DO
            ELSE
              DO M=1,MQSER(NS)
                READ(1,*,IOSTAT=ISO)TQSER(M,NS),(QSER(M,K,NS), K=1,KC)
                IF(ISO.GT.0) GOTO 860
              END DO
            END IF
          END DO
        END IF
        CLOSE(1)
        goto 862
860     WRITE(6,861)
861     FORMAT('  READ ERROR FOR FILE QSER.INP ')
862     continue
!
!=======================================================================
!      hlopen parameters
!=======================================================================
        Ihl_handle = 0
        Ihl_debug = 0
        Ihl_mode = 1	!Ihl_mode=0 to READ from dll hyd file, =1 to WRITE to dll hyd file
        call hlsetdebug(Ihl_debug)
        call hlopen(HYDFIL,Ihl_mode,Ihl_handle,ierror)
        if(ierror.gt.0) then
          call hlgetlasterror(errstring)
          write(6,6000)ierror,errstring,1
          stop
        end if
        call hlsetlanguage(Ihl_handle,1,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,2
          stop
        end if
        DES=' 666  '//char(0)
        call hladddescription(Ihl_handle,0,DES,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,3
          stop
        end if
        MOD='Created by: JL666  '//char(0)
        call hladddescription(Ihl_handle,1,MOD,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,4
          stop
        end if
        call hlsetcreator(Ihl_handle, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,5
          stop
        end if
        j=1461*(iyear+4800+int((imon-14)/12))/4+367*(imon-2-12*
     &    int((imon-14)/12))/12-(3*((iyear+4800+int((imon-14)/12)
     &    +100)/100))/4+iday-32075+tbegin+idays
	istartyear=100*(int(4*(j+68569)/146097)-49)+int(4000*
     &   (int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/
     &   1461001)+int(int(80*int(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*
     &   int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
	istartmonth=int(80*int(int(j+68569-(146097*int(4*int(j+68569)/
     &   146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+
     &   68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)+2-12*int(int(80*
     &   int(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)
     &   -1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
	istartday=int(int(j+68569-(146097*int(4*int(j+68569)/146097)
     &   +3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)+1)/1461001)/4+31)-2447*int(80*int(
     &   int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)
     &   -1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/
     &   146097)+3)/4)+1)/1461001)/4+31)/2447)/80
	istarthour=0
	istartminute=0
	istartsecond=0

	    IBEGIN=IDAYS*NTSPTC
	    AUX=FLOAT(IBEGIN)/FLOAT(NTSMMT)
	    NAUX=INT(AUX)
	      IF(AUX.GT.NAUX) then
	        AUX=(NAUX+1.-AUX)
                tsec=AUX*NTSMMT*TCON/NTSPTC
                istartsecond=INT(tsec)
                IF(tsec.GT.60) THEN
                  TMIN=tsec/60.0
                  istartminute=INT(TMIN)
                  istartsecond=INT((TMIN-istartminute)*60.)
                  IF(istartminute.GT.60) THEN
                    THOUR=FLOAT(istartminute)/60.0
                    istarthour=INT(THOUR)
                    istartminute=INT((THOUR-istarthour)*60.)
                  ENDIF
                ENDIF
              ENDIF

        call hlsetseedmoment(Ihl_handle,istartmonth,istartday,
     +     istartyear,istarthour,istartminute,istartsecond,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,6
          stop
        end if
	call hlsetnumlayers(Ihl_handle,kc,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,7
          stop
        end if

        NJUN=0
        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          nlay=kc-kgvcp(l)+1
          NJUN=NJUN+nlay
        END DO
        if(njun.gt.nsg) then
        write(6,500)NJUN,NSG
        STOP
        END IF
500     FORMAT('THE NUMBER OF WASP SEGMENTS IN YOUR APPLICATION',I6,1x,'IS GREATER THAN THE ARRAY DIMENSION:',I7)
        call hlsetnumsegments(Ihl_handle,njun,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,8
          stop
        end if
        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          sn(l)='      '//char(0)
        END DO
        OPEN(94,FILE='segname.inp',iostat=ios,STATUS='old')
        IF(ios.eq.0) then
          do i=1,4
            read(94,*)
          end do
          do kk=1,la
            READ(94,*,err=111)I,J,SN1
            L=LIJ(I,J)
            sn(l)=sn1
          end do
        end if
111     CLOSE(94)
        I=0

        DO K=KC,1,-1
          DO LT=2,LALT
            l=lij(illt(lt),jllt(lt))
            if (LGVCP(L,K)) THEN
    	      I=I+1
	      LAUX(ILLT(LT),JLLT(LT),K)=I
              if(il(l).ge.100) then
                WRITE(itext,"(I3)")IL(L)
              elseif(il(l).ge.10) then
                WRITE(itext,"(I2)")IL(L)
              else
                WRITE(itext,"(I1)")IL(L)
              end if
              if(jl(l).ge.100) then
                WRITE(jtext,"(I3)")jL(L)
              elseif(jl(l).ge.10) then
                WRITE(jtext,"(I2)")jL(L)
              else
                WRITE(jtext,"(I1)")jL(L)
              end if
              if(k.ge.100) then
                WRITE(ktext,"(I3)")k
              elseif(k.ge.10) then
                WRITE(ktext,"(I2)")k
              else
                WRITE(ktext,"(I1)")k
              end if
	      sn2=' I='//itext//' J='//jtext//' K='//ktext
    	      segname=sn2//sn(l)//char(0)
	      call hlsetsegname(ihl_handle,i,segname,ierror)
            END IF
          END DO
        END DO

        OPEN(1,FILE='seginfo.wsp',STATUS='UNKNOWN')
        DO K=KC,2,-1
          DO LT=2,LALT
            l=lij(illt(lt),jllt(lt))
            IF (LGVCP(L,K)) THEN
              LU=LAUX(ILLT(LT),JLLT(LT),K)
              IF (LGVCP(L,K-1)) THEN
	        LD=LAUX(ILLT(LT),JLLT(LT),K-1)
              ELSE
	        LD=0
              ENDIF
              SUR=2
              IF(K.EQ.KC) SUR=1
              Write(1,*)LU,SUR,LD
            END IF
          END DO
        END DO

        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          IF (LGVCP(L,1)) THEN
            LU=LAUX(ILLT(LT),JLLT(LT),1)
            LD=0
            SUR=2
            Write(1,*)LU,SUR,LD
          END IF
        END DO
        close(1)
        NCHNH=0
        NCHNV=0
        DO K=KC,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (LGVCU(L,K)) THEN
              NCHNH=NCHNH+1
            END IF
            IF (IJCTLT(I+1,J).EQ.8) THEN
              IF (LGVCU(LIJ(I+1,J),K)) THEN
                NCHNH=NCHNH+1
              END IF
            END IF
            IF (LGVCV(L,K)) THEN
              NCHNH=NCHNH+1
            END IF
            IF (IJCTLT(I,J+1).EQ.8) THEN
              IF (LGVCV(LIJ(I,J+1),K)) THEN
                NCHNH=NCHNH+1
              END IF
            END IF
          END DO
        END DO
        DO LT=2,LALT
          DO k=1,KS
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (LGVCW(L,K)) THEN
              NCHNV=NCHNV+1
            END IF
          END DO
        END DO
        NCHN=NCHNH+NCHNV


	DO LQ=1,NQSIJ
          I=IQS(LQ)
          J=JQS(LQ)
          L=LIJ(I,J)
	  IF(LIJLT(I,J).GT.0) THEN
            NS=NQSERQ(LQ)
	    DO K=KGVCP(L),KC
	      IF(flagwaspbC(ns,k).NE.1) THEN
	        NCHN=NCHN+1
	      END IF
	    END DO
	  END IF
	END DO
	NQ=NQCTL
	DO LQ=1,NQCTL
	  I=IQCTLU(LQ)
	  J=JQCTLU(LQ)
          L=LIJ(I,J)
	  IF(LIJLT(I,J).GT.0) THEN
	    DO K=KGVCP(L),KC
              NCHN=NCHN+1
	    END DO
	  END IF
	END DO

        if(NCHN.gt.NF) then
        write(6,600)NCHN,NF
        STOP
        END IF
600     FORMAT('THE NUMBER OF WASP FLOWS IN YOUR APPLICATION',I6,1X,'IS GREATER THAN THE ARRAY DIMENSION:',I7)
        call hlsetnumflowpaths(Ihl_handle, NCHN, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,9
          stop
        end if
        LWASP=0
        DO K=KC,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (LGVCU(L,K)) THEN
              LWASP=LWASP+1
              LDTM=LAUX(I,J,K)
              LUTM=LDTM-1
              IF (IJCTLT(I-1,J).EQ.8) LUTM=0
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,10
                  stop
               end if
            END IF
            IF (IJCTLT(I+1,J).EQ.8) THEN
              IF (LGVCU(L+1,K)) THEN
	        LWASP=LWASP+1
                LDTM=0
                LUTM=LAUX(I,J,K)
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,11
                  stop
               end if
              END IF
            END IF
          END DO

          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (LGVCV(L,K)) THEN
              LWASP=LWASP+1
              LSLT=LSCLT(LT)
              LDTM=LAUX(I,J,K)
!
!LZ              LUTM=LAUX(ILLT(LSLT),JLLT(LSLT),K)
!LZ              IF(IJCTLT(I,J-1).EQ.8) LUTM=0
!
              IF(IJCTLT(I,J-1).EQ.8) THEN          !DRBC LZ, 7/23/2019 ADJUST THE LOGIC SEQUENCE IF A SOUTH CELL IS WITH ID 8
                LUTM=0
              ELSE
                LUTM=LAUX(ILLT(LSLT),JLLT(LSLT),K)
              END IF
!
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,2,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,12
                  stop
               end if
            END IF
            IF (IJCTLT(I,J+1).EQ.8) THEN
              LN=LNC(L)
              IF (LGVCV(LN,K)) THEN
                LWASP=LWASP+1
                LDTM=0
                LUTM=LAUX(I,J,K)
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,2,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,13
                  stop
               end if
              END IF
            END IF
          END DO
        END DO

        DO K=KC,1,-1
          DO LT=1,NQSIJ
            I=IQS(LT)
            J=JQS(LT)
            L=LIJ(I,J)
            IF(LIJLT(I,J).EQ.0) GOTO 100
            NS=NQSERQ(Lt)
            IF(flagwaspbC(ns,k).EQ.1) GOTO 100
            IF (LGVCP(L,K)) THEN
              LWASP=LWASP+1
              LDTM=LAUX(I,J,K)
              LUTM=0
              call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
              if(ierror .gt. 0)then
                call hlgetlasterror(errstring)
                write(6,6000) ierror, errstring,14
                stop
              end if
            END IF
100       END DO
        END DO

        DO K=KC,1,-1
          DO LT=1,NQCTL
            I=IQCTLU(LT)
            J=JQCTLU(LT)
            LUTM=LAUX(I,J,K)
            I=IQCTLD(LT)
            J=JQCTLD(LT)
            LDTM=LAUX(I,J,K)
	      IF(LUTM.EQ.0.AND.LDTM.EQ.0) GOTO 200
              LWASP=LWASP+1
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,15
                  stop
               end if
200       END DO
        END DO

        IF(KC.GT.1)THEN
          DO K=KS,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF (LGVCW(L,K)) THEN
                LWASP=LWASP+1
                LUTM=LAUX(I,J,K+1)
                LDTM=LAUX(I,J,K)
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,3,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,16
                  stop
               end if
              END IF
            END DO
          END DO
        ENDIF

	inumsegconsts=3  !volume,depth,velocity
        if(istran(2).GE.1) inumsegconsts=4   !water temperature modeled in EFDC and transfered
	if(istran(1).GE.1) inumsegconsts=5   !salinity modeled in EFDC and transfered, 05/08/2019, DRBC LZ
	inumsegconsts=6                      !TKE DISSIPATION RATE, 12/20/2012, DRBC LZ       !HNR_GHD - 4/2022 HOTWIRED TO 6 
        call hlsetnumsegconsts(Ihl_handle, inumsegconsts, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,17
          stop
        end if
        call hlsetnumfpconsts(Ihl_handle, 3, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,18
          stop
        end if
        call hlsetsegconsttype(Ihl_handle, 1, 0, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,19
          stop
        end if
        call hlsetsegconsttype(Ihl_handle, 2, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,20
          stop
        end if
        call hlsetsegconsttype(Ihl_handle, 3, 2, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,21
          stop
        end if
        if(istran(2).GE.1) then
        call hlsetsegconsttype(Ihl_handle, 4, 3, ierror)
        if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,22
          stop
        end if
        end if
        if(istran(1).GE.1) then
           call hlsetsegconsttype(Ihl_handle, 5, 4, ierror)
           if(ierror .gt. 0)then
             call hlgetlasterror(errstring)
             write(6,6000) ierror, errstring,23
             stop
           end if
		 end if
!========================================================================
!              Li Need to Add this Block
!========================================================================
!	  if(istran(1).GE.1) then                                !HNR_GHD 4/2022  UNLINKED DISPRATE FROM SALINITY 
           call hlsetsegconsttype(Ihl_handle, 6, 5, ierror)
           if(ierror .gt. 0)then
             call hlgetlasterror(errstring)
             write(6,6000) ierror, errstring,235
             stop
           end if
!        end if

        call hlsetfpconsttype(Ihl_handle, 1, 0, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,24
          stop
        end if
        call hlsetfpconsttype(Ihl_handle, 1, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,25
          stop
        end if
        call hlsetfpconsttype(Ihl_handle, 1, 2, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,26
          stop
        end if
	call hlsetvartimestep(Ihl_handle,0,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,27
          stop
        end if
	call hlsethydtimestep(Ihl_handle,dt,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,28
          stop
        end if
	rinterval=(dt*NTSMMT)/86400.
	call hlsetupdateint(Ihl_handle,rinterval,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,29
          stop
        end if
	call hlsethydtowaspratio(Ihl_handle,NTSMMT,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,30
          stop
        end if
        IF(IDAYS.EQ.0) THEN
          IF(ISRESTI.EQ.0) THEN
          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF (LGVCP(L,K)) THEN
                LWASP=LWASP+1
                SegVel(LWASP)=0.0

              IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                SegDepth(LWASP)=HP(L)*DZC(K)*GVCSCLP(L)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=SAL(L,K)
                SegTemp(LWASP)=TEM(L,K)
                SegDispr(LWASP)=DISPRATE(L,K)        !DRBC LZ 12/20/2020
              ELSE
                SegDepth(LWASP)=HDRY*DZC(K)*GVCSCLP(L)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=0.0
                SegTemp(LWASP)=TATMT(L)
                SegDispr(LWASP)=0.0                  !DRBC LZ 12/20/2020
              END IF

                IF(NTSMMT.LT.NTSPTC) THEN

                IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                  SegDepth(LWASP)=HMP(L)*DZC(K)*GVCSCLP(L)
                  SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                  SegSalt(LWASP)=SALINIT(L,K)
                  SegTemp(LWASP)=TEMINIT(L,K)
                  SegDispr(LWASP)=DISPRATE(L,K)        !DRBC LZ 12/20/2020
                ELSE
                  SegDepth(LWASP)=HDRY*DZC(K)*GVCSCLP(L)
                  SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                  SegSalt(LWASP)=0.0
                  SegTemp(LWASP)=TATMT(L)
                  SegDispr(LWASP)=0.0                  !DRBC LZ 12/20/2020
                END IF

                END IF
              END IF
            END DO
          END DO
          DO I=1,NCHN
	    FLOW(I)=0.0
	    CRNU(I)=0.0
	    BRINTT(I)=0.0
          END DO
	  call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,31
            stop
          end if
	  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,32
            stop
          end if
	  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,33
            stop
          end if
          if(istran(2).GE.1) then
	    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,34
              stop
            end if
          end if
          if(istran(1).GE.1) then
	    call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,35
              stop
            end if
          end if
	  call hlsetseginfo(Ihl_handle,6,SegDispr,ierror)   !DRBC LZ, 12/20/2020, OUTPUT TKE DISSIPATION RATE
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,351
            stop
          end if                    
	  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,36
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,37
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,38
            stop
          end if
	  call hlmomentcomplete(Ihl_Handle,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,39
            stop
          end if
          ELSE
          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF (LGVCP(L,K)) THEN
                LWASP=LWASP+1
                LN=LNC(L)
                VELX=0.5*(U(L,K)+U(L+1,K))
                VELY=0.5*(V(L,K)+V(LN,K))
                VELZ=0.5*(W(L,K-1)+W(L,K))
                VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
              
              IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                SegVel(LWASP)=VELMAG
                SegDepth(LWASP)=HP(L)*DZC(K)*GVCSCLP(L)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=SAL(L,K)
                SegTemp(LWASP)=TEM(L,K)
                SegDispr(LWASP)=DISPRATE(L,K)        !DRBC LZ 12/20/2020
              ELSE
                SegVel(LWASP)=0.0
                SegDepth(LWASP)=HDRY*DZC(K)*GVCSCLP(L)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=0.0
                SegTemp(LWASP)=TATMT(L)
                SegDispr(LWASP)=0.0                  !DRBC LZ 12/20/2020
              END IF
              
              END IF
            END DO
          END DO
! Advection and dispersion in the X-direction:
          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              ADDLW=0.0
              IF (LGVCU(L,K)) THEN
                vol1=DYP(L)*AH(L,K)*DZC(K)*GVCSCLP(L)*HP(L)/DXP(L)
                vol2=DYP(L-1)*AH(L-1,K)*DZC(K)*GVCSCLP(L-1)*HP(L-1)/DXP(L-1)
                ADDLW=0.5*(vol1+vol2)
                  vol1=DXYP(L)*HP(L)*DZC(K)*GVCSCLP(L)         
                  vol2=DXYP(L-1)*HP(L-1)*DZC(K)*GVCSCLP(L-1)       
                  if (vol1.LT.vol2) vol2=vol1                    
                  ad=vol2/dt*adcoeff                             
                  if (addlw.GT.ad) then                          
                    addlw=ad                                     
                  end if                                         
                  if(ishdmf.lt.2) then                           
                    addlw=0.                                     
                  end if                                         
                LWASP=LWASP+1
                
                IF((ISCDRY(L).EQ.1).OR.(ISCDRY(L-1).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
                  FLOW(LWASP)=0.0
	            CRNU(LWASP)=0.0
  	            BRINTT(LWASP)=0.0
                ELSE
                  vol1=DYP(L)*HP(L)*DZC(K)*GVCSCLP(L)
                  vol2=DYP(L-1)*HP(L-1)*DZC(K)*GVCSCLP(L-1)
                  FLOW(LWASP)=U(L,K)*0.5*(vol1+vol2)
	            vol1=1.0/DXP(L)
	            vol2=1.0/DXP(L-1)
	            CRNU(LWASP)=U(L,K)*0.5*(vol1+vol2)
  	            BRINTT(LWASP)=ADDLW
                END IF
                
              END IF
              IF (IJCTLT(I+1,J).EQ.8) THEN
                IF (LGVCU(L+1,K)) THEN
                  LWASP=LWASP+1
                  
                IF((ISCDRY(L).EQ.1).OR.(ISCDRY(L+1).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
	            FLOW(LWASP)=0.0
	            CRNU(LWASP)=0.0
                  BRINTT(LWASP)=0.0
                ELSE
                  vol1=DYP(L)*HP(L)*DZC(K)*GVCSCLP(L)
                  vol2=DYP(L+1)*HP(L+1)*DZC(K)*GVCSCLP(L+1)
                  FLOW(LWASP)=U(L+1,K)*0.5*(vol1+vol2)
                  vol1=1.0/DXP(L)
	            vol2=1.0/DXP(L+1)
	            CRNU(LWASP)=U(L+1,K)*0.5*(vol1+vol2)
                  BRINTT(LWASP)=0.0             
                END IF
                
                END IF
              END IF
            END DO
!
! Advection and dispersion in the Y-direction:
!
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              ADDLS=0.0
              LS=LSC(L)
              IF (LGVCV(L,K)) THEN
                vol1=DXP(L)*AH(L,K)*DZC(K)*GVCSCLP(L)*HP(L)/DYP(L)
                vol2=DXP(LS)*AH(LS,K)*DZC(K)*GVCSCLP(LS)*HP(LS)/DYP(LS)
                ADDLS=0.5*(vol1+vol2)
                  vol1=DXYP(L)*HP(L)*DZC(K)*GVCSCLP(L)         
                  vol2=DXYP(Ls)*HP(Ls)*DZC(K)*GVCSCLP(Ls)      
                  if (vol1.LT.vol2) vol2=vol1                    
                  ad=vol2/dt*adcoeff                             
                  if (addls.GT.ad) then                          
                    addls=ad                                     
                  end if                                         
                  if(ishdmf.lt.2) then                           
                    addls=0.                                     
                  end if                                         
	        LWASP=LWASP+1
              
                IF((ISCDRY(L).EQ.1).OR.(ISCDRY(LS).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
                  FLOW(LWASP)=0.0
	            CRNU(LWASP)=0.0
                  BRINTT(LWASP)=0.0
                ELSE
                  vol1=DXP(L)*HP(L)*DZC(K)*GVCSCLP(L)
                  vol2=DXP(LS)*HP(LS)*DZC(K)*GVCSCLP(LS)
                  FLOW(LWASP)=V(L,K)*0.5*(vol1+vol2)
                  vol1=1.0/DYP(L)
	            vol2=1.0/DYP(LS)
	            CRNU(LWASP)=V(L,K)*0.5*(vol1+vol2)
                  BRINTT(LWASP)=ADDLS
                END IF
                 
              END IF
              IF (IJCTLT(I,J+1).EQ.8) THEN
                LN=LNC(L)
                IF (LGVCV(LN,K)) THEN
  	          LWASP=LWASP+1
              
                IF((ISCDRY(L).EQ.1).OR.(ISCDRY(LN).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
                  FLOW(LWASP)=0.0
	            CRNU(LWASP)=0.0
                  BRINTT(LWASP)=0.0
                ELSE
                  vol1=DXP(L)*HP(L)*DZC(K)*GVCSCLP(L)
                  vol2=DXP(LN)*HP(LN)*DZC(K)*GVCSCLP(LN)
                  FLOW(LWASP)=V(LN,K)*0.5*(vol1+vol2)
                  vol1=1.0/DYP(L)
	            vol2=1.0/DYP(LN)
	            CRNU(LWASP)=V(LN,K)*0.5*(vol1+vol2)
	            BRINTT(LWASP)=0.0                                
                END IF
                
                END IF
              END IF
            END DO
          END DO

! Advection in input flows
          DO K=KC,1,-1
            DO LT=1,nqsij
              I=Iqs(LT)
              J=Jqs(LT)
              IF(LIJLT(I,J).EQ.0) GOTO 310
              NS=NQSERQ(Lt)
              L=LQS(Lt)
              IF(flagwaspbC(ns,k).EQ.1) GOTO 310                       
              IF (LGVCP(L,K)) THEN
                LWASP=LWASP+1
                FLOW(LWASP)=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
                CRNU(LWASP)=flow(LWASP)/DXp(L)/dyp(l)/(dzc(k)*HP(L)*GVCSCLP(L))
                BRINTT(LWASP)=0.0
              END IF
310         END DO
          END DO
!  Advection in structure flows
          DO K=KC,1,-1
            DO LT=1,nqctl
              Iu=Iqctlu(LT)
              Ju=Jqctlu(LT)
              Lu=Lij(iu,ju)
              Id=Iqctld(LT)
              Jd=Jqctld(LT)
              Ld=Lij(id,jd)
              IF(LU.EQ.0.AND.LD.EQ.0) GOTO 410
              flowx=RQCMUL(Lt)*QCTLT(K,Lt)
              UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*HP(L)*GVCSCLP(L))
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWX
	      CRNU(LWASP)=UDDXTMP
              BRINTT(LWASP)=0.0
410         END DO
          END DO
!
! Advection and dispersion in the Z-direction:
!
          IF (KC.GT.1) THEN
            DO K=KS,1,-1
              DO LT=2,LALT
                I=ILLT(LT)
                J=JLLT(LT)
                L=LIJ(I,J)
                ADDL=ab(l,k)*hp(l)                         
                if (addl.gt.abwmx(l)) then
! ORIGINAL CODE                  ab(l,k)=abwmx(l)/hp(l)  
                  ABTMP=abwmx(l)/hp(l)                      !8/26/2020, DRBC LZ, AVOID FEEDBACK TO OTHER EFDC SUBROUTINES
                ELSE                                        !8/26/2020, DRBC LZ 
                  ABTMP=AB(L,K)                             !8/26/2020, DRBC LZ 
                end if
                addl=0.0
                IF (LGVCW(L,K)) THEN
! ORIGINAL CODE                  ADDL=DXYP(L)*AB(L,K)*DZIG(K)*GVCSCLPI(L) 
                  ADDL=DXYP(L)*ABTMP*DZIG(K)*GVCSCLPI(L)    !8/26/2020, DRBC LZ              
                  vol1=DXYP(L)*HP(L)*DZC(K)*GVCSCLP(L)               
                  vol2=DXYP(L)*HP(L)*DZC(K+1)*GVCSCLP(L)             
                  if (vol1.LT.vol2) vol2=vol1                        
                  ad=adcoeff*vol2/dt                                 
                  if (addl.GT.ad) then                               
                    addl=ad                                          
                  end if                                             
                  LWASP=LWASP+1

                IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                  FLOW(LWASP)=-DXYP(L)*W(L,K)
                  CRNU(LWASP)=W(L,K)*DZIG(K)*GVCSCLPI(L)/HP(L)
	            BRINTT(LWASP)=ADDL
                ELSE
  	            FLOW(LWASP)=0.0
                  CRNU(LWASP)=0.0
	            BRINTT(LWASP)=0.0
                END IF
                
                END IF
              END DO
            END DO
          END IF
	  call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,40
            stop
          end if
	  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,41
            stop
          end if
	  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,42
            stop
          end if
          if(istran(2).GE.1) then
	    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,43
              stop
            end if
          end if
          if(istran(1).GE.1) then                  ! 06/28/2019 DRBC LZ, OUTPUT INITIAL SALINITY DURING A HOT-START
            call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
            if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,44
              stop
            end if
          end if           
          call hlsetseginfo(Ihl_handle,6,SegDispr,ierror)   !DRBC LZ, 12/20/2020, OUTPUT TKE DISSIPATION RATE
	  if(ierror .gt. 0)then
	    call hlgetlasterror(errstring)
	    write(6,6000) ierror, errstring,441
	    stop
          end if 
	  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,45
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,46
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,47
            stop
          end if
	  call hlmomentcomplete(Ihl_Handle,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,48
            stop
          end if
        END IF
      END IF
      GOTO 3000
      END IF
      IF(N.LT.IBEGIN) GOTO 3000
!
!----------------------------------------------------------------------C
!
! **  WRITE TIME STEP DATA
!
! Advection and dispersion in the X-direction:
!
      LWASP=0
      DO K=KC,1,-1
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          ADDLW=0.0
          IF (LGVCU(L,K)) THEN
            LW=L-1
            vol1=DYP(L)*DZC(K)*HLPF(L)*GVCSCLP(L)
            vol2=DYP(LW)*DZC(K)*HLPF(LW)*GVCSCLP(LW)
            ADDLW=0.5*(vol1+vol2)*AHULPF(L,K)*DXIU(L)
            vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                
            vol2=DXYP(LW)*HLPF(LW)*DZC(K)*GVCSCLP(LW)             
            if (vol1.LT.vol2) vol2=vol1                           
            ad=adcoeff*vol2/dt                                    
            if (addlw.GT.ad) then                                 
              addlw=ad                                            
            end if                                                
            if(ishdmf.lt.2) then                                  
              addlw=0.                                            
            end if                                                
            vol1=DYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
            vol2=DYP(L-1)*HLPF(L-1)*DZC(K)*GVCSCLP(L-1)
            FLOWX=UHLPF(L,K)*DYU(L)*DZC(K)*GVCSCLU(L)
            UDDXTMP=UHLPF(L,K)*DXIU(L)/HU(L)
            LWASP=LWASP+1
            FLOW(LWASP)=FLOWX
            CRNU(LWASP)=UDDXTMP
            BRINTT(LWASP)=ADDLW
          END IF
          IF (IJCTLT(I+1,J).EQ.8) THEN
            IF (LGVCU(L+1,K)) THEN
              vol1=DYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
              vol2=DYP(L+1)*HLPF(L+1)*DZC(K)*GVCSCLP(L+1)
              FLOWX=UHLPF(L+1,K)*DYU(L+1)*DZC(K)*GVCSCLU(L+1)
              UDDXTMP=UHLPF(L+1,K)*DXIU(L+1)/HU(L+1)
              vol1=DYP(L)*DZC(K)*HLPF(L)*GVCSCLP(L)
              vol2=DYP(L+1)*DZC(K)*HLPF(L+1)*GVCSCLP(L+1)
              ADDLW=0.5*(vol1+vol2)*AHULPF(L+1,K)*DXIU(L)
              vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                
              vol2=DXYP(L+1)*HLPF(L+1)*DZC(K)*GVCSCLP(L+1)             
              if (vol1.LT.vol2) vol2=vol1                           
              ad=adcoeff*vol2/dt                                    
              if (addlw.GT.ad) then                                 
                addlw=ad                                            
              end if                                                
              if(ishdmf.lt.2) then                                  
                addlw=0.                                            
              end if                                                
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWX
	      CRNU(LWASP)=UDDXTMP
	      BRINTT(LWASP)=ADDLW
            END IF
          END IF
        END DO
!
! Advection and dispersion in the Y-direction:
!
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          ADDLS=0.0
          IF (LGVCV(L,K)) THEN
            LS=LSC(L)
            vol1=DXP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                       
            vol2=DXP(LS)*HLPF(LS)*DZC(K)*GVCSCLP(LS)                     
            ADDLS=AHVLPF(L,K)*0.5*(vol1+vol2)*DYIV(L)
            vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                       
            vol2=DXYP(LS)*HLPF(LS)*DZC(K)*GVCSCLP(LS)                    
            if (vol1.LT.vol2) vol2=vol1                                  
            ad=adcoeff*vol2/dt                                           
            if (addls.GT.ad) then                                        
              addls=ad                                                   
            end if                                                       
            if(ishdmf.lt.2) then                                         
              addls=0.                                                   
            end if                                                       
            vol1=DXP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
            vol2=DXP(LS)*HLPF(LS)*DZC(K)*GVCSCLP(LS)
            FLOWY=VHLPF(L,K)*DXV(L)*DZC(K)*GVCSCLV(L)
            VDDYTMP=VHLPF(L,K)*DYIV(L)/HV(L)
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWY
	    CRNU(LWASP)=VDDYTMP
	    BRINTT(LWASP)=ADDLS
          END IF
          IF (IJCTLT(I,J+1).EQ.8) THEN
            LN=LNC(L)
            IF (LGVCV(LN,K)) THEN
              vol1=DXP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
              vol2=DXP(LN)*HLPF(LN)*DZC(K)*GVCSCLP(LN)
              FLOWY=VHLPF(LN,K)*DXV(LN)*DZC(K)*GVCSCLV(LN)
              VDDYTMP=VHLPF(LN,K)*DYIV(LN)/HV(LN)
              vol1=DXP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                       
              vol2=DXP(LN)*HLPF(LN)*DZC(K)*GVCSCLP(LN)                     
              ADDLS=AHVLPF(LN,K)*0.5*(vol1+vol2)*DYIV(LN)
              vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                       
              vol2=DXYP(LN)*HLPF(LN)*DZC(K)*GVCSCLP(LN)                    
              if (vol1.LT.vol2) vol2=vol1                                  
              ad=adcoeff*vol2/dt                                           
              if (addls.GT.ad) then                                        
                addls=ad                                                   
              end if                                                       
              if(ishdmf.lt.2) then                                         
                addls=0.                                                   
              end if                                                       
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWY
              CRNU(LWASP)=VDDYTMP
	      BRINTT(LWASP)=ADDLS
            END IF
          END IF
        END DO
      END DO

! Advection in input flows
      DO K=KC,1,-1
        DO LT=1,nqsij
          I=Iqs(LT)
          J=Jqs(LT)
          IF(LIJLT(I,J).EQ.0) GOTO 300
          NS=NQSERQ(Lt)
          L=LQS(Lt)
          IF(flagwaspbC(ns,k).EQ.1) GOTO 300                       
          IF (LGVCP(L,K)) THEN
            flowx=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
            UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*GVCSCLP(L)*HLPF(L))
            LWASP=LWASP+1
	    FLOW(LWASP)=FLOWX
	    CRNU(LWASP)=UDDXTMP
            BRINTT(LWASP)=0.0
          END IF
300     END DO
      END DO
!  Advection and dispersion in structure flows
      DO K=KC,1,-1
        DO LT=1,nqctl
          Iu=Iqctlu(LT)
          Ju=Jqctlu(LT)
          Lu=Lij(iu,ju)
          Id=Iqctld(LT)
          Jd=Jqctld(LT)
          Ld=Lij(id,jd)
          IF(LU.EQ.0.AND.LD.EQ.0) GOTO 400
          flowx=RQCMUL(Lt)*QCTLT(K,Lt)
          UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*GVCSCLP(L)*HmP(L))
          if(ishdmf.lt.2) then                                  
            addls=0.                                            
          end if                                                
          LWASP=LWASP+1
	  FLOW(LWASP)=FLOWX
	  CRNU(LWASP)=UDDXTMP
	  BRINTT(LWASP)=0.0
400     END DO
      END DO
!
! Advection and dispersion in the Z-direction:
!
      IF (KC.GT.1) THEN
        DO K=KS,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            ADDL=ablpf(l,k)*hlpf(l)                                
            if(addl.gt.abwmx(l)) then
! ORIGINAL CODE               ablpf(l,k)=abwmx(l)/hlpf(l) 
              ABTMP=abwmx(l)/hlpf(l)                      !8/26/2020, DRBC LZ, AVOID FEEDBACK TO OTHER EFDC SUBROUTINES
            ELSE                                          !8/26/2020, DRBC LZ
              ABTMP=ABLPF(L,K)                            !8/26/2020, DRBC LZ
            end if
            ADDL=0.0
            IF (LGVCW(L,K)) THEN
! ORIGINAL CODE              ADDL=DXYP(L)*ABLPF(L,K)*DZIG(K)*GVCSCLPI(L)    
              ADDL=DXYP(L)*ABTMP*DZIG(K)*GVCSCLPI(L)      !8/26/2020, DRBC LZ        
              vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)               
              vol2=DXYP(L)*HLPF(L)*DZC(K+1)*GVCSCLP(L)             
              if (vol1.LT.vol2) vol2=vol1                          
              ad=adcoeff*vol2/dt                                   
              if (addl.GT.ad) then                                 
                addl=ad                                            
              end if                                               
              FLOWZ=-DXYP(L)*wlpf(l,k)
              WDDZTMP=wlpf(l,k)*DZIG(K)*GVCSCLPI(L)/HLPF(L)
              KPTMP=K+1
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWZ
	      CRNU(LWASP)=WDDZTMP
	      BRINTT(LWASP)=ADDL
            END IF
          END DO
        END DO
      END IF
!
! Segment Properties:
!
      LCELTMP=0
      DO K=KC,1,-1
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          LN=LNC(L)
          IF (LGVCP(L,K)) THEN
            LCELTMP=LCELTMP+1
            VOLUM=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)

            IF(NTSMMT.LT.NTSPTC) THEN
              
              IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                VOLUM=DXYP(L)*HP(L)*DZC(K)*GVCSCLP(L)
              ELSE
                VOLUM=DXYP(L)*HDRY*DZC(K)*GVCSCLP(L)
              END IF
            END IF
          
            DEPTH=HLPF(L)*DZC(K)*GVCSCLP(L)
            VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K) +UHLPF(L+1,K)+SVPT*UVPT(L+1,K))/HLPF(L)
            VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K) +VHLPF(LN,K)+SVPT*VVPT(LN,K) )/HLPF(L)
            VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1) +WLPF(L,K)+SVPT*WVPT(L,K) )
            VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
            SegSalt(LCELTMP)=SALLPF(L,K)
            SegTemp(LCELTMP)=TEMLPF(L,K)
	    SegDispr(LCELTMP)=DISPRLPF(L,K)        !DRBC LZ 12/20/2020            
	    SegVolume(LCELTMP)=VOLUM
	    SegDepth(LCELTMP)=DEPTH
	    SegVel(LCELTMP)=VELMAG
          END IF
        END DO
      END DO
!
	   call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,49
            stop
         end if
	   call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,50
            stop
         end if
	   call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,51
            stop
         end if
       if(istran(2).GE.1) then
	   call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,52
            stop
         end if
       end if
       if(istran(1).GE.1) then                  ! 06/13/2019 DRBC LZ, OUTPUT TEMPORAL SALINITY
           call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
           if(ierror .gt. 0)then
             call hlgetlasterror(errstring)
             write(6,6000) ierror, errstring,53
             stop
           end if
       end if 
       call hlsetseginfo(Ihl_handle,6,SegDispr,ierror)   !DRBC LZ, 12/20/2020, OUTPUT TKE DISSIPATION RATE
       if(ierror .gt. 0)then
         call hlgetlasterror(errstring)
         write(6,6000) ierror, errstring,531
         stop
       end if       
	   call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,54
            stop
         end if
	   call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,55
            stop
         end if
	   call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,56
            stop
         end if
	   call hlmomentcomplete(Ihl_Handle,ierror)
         if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,57
            stop
         end if

!
3000  JSWASP=0
6000  format('Error ',I10, ' : ', A30,I3)
      RETURN
      END