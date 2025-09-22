      SUBROUTINE WASPHYDROLINK
!
! DATE MODIFIED     BY                 DATE APPROVED    BY  
!
! 05/08/2019        DRBC LZ - ADD INITIAL  SALINITY OUTPUT TO HYDRO-LINKAGE FILE
! 06/13/2019        DRBC LZ - ADD TEMPORAL SALINITY OUTPUT TO HYDRO-LINKAGE FILE
! 06/24/2019        HNR_GHD - ADJUSTED VALUES FOR DRY CELLS WHEN W/D IS USED
! 06/28/2019        DRBC LZ - ADD INITIAL SALINITY OUTPUT TO HYDRO-LINKAGE FILE DURING A HOT-START
! 12/20/2021        DRBC LZ - ADDED TKE DISSIPATION RATE 
! 04/21/2022        HNR_GHD - ADDED THE SETTIN OF DISSIPATION RATE IN THE HYD FILE
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      INCLUDE 'ALLset.INT'

      parameter(nf=44000,nsg=20000)
!
      DIMENSION LAUX(ICM,JCM,KCM)
      INTEGER LU,LD,IU,ID,JU,JD,NAUX
      INTEGER istartyear,istartmonth,istartday,istarthour,istartminute,istartsecond,              &
              Ihl_debug,Ihl_mode,inumsegconsts,j,FLAGWASPBC(NQSERM,KCM)
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
11      continue
        DO LT=2,LALT                               
          read(1,*,err=12)i,j,abwmax               
          l=lij(i,j)                               
          abwmx(l)=abwmax                          
        END DO                                     
12      continue                                   

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
                READ(1,*,IOSTAT=ISO)                   
                IF(ISO.GT.0) GOTO 860                                     
              END DO                                                      
            ELSE                                                          
              DO M=1,MQSER(NS)                                            
                READ(1,*,IOSTAT=ISO)    
                IF(ISO.GT.0) GOTO 860                                     
              END DO                                                      
            END IF                                                        
          END DO                                                          
        END IF                                                            
!                                                                       
        CLOSE(1)                                                          
        goto 862                                                          
860     WRITE(6,861)                                                      
861     FORMAT('  READ ERROR FOR FILE QSER.INP ')                         
862     continue                                                          
!
!**********************************************************************C
!
! **  DEFINE EFDC-WASP CORRESPONDENCE AND INITIALIZE FILES
!
!----------------------------------------------------------------------C
!
        Ihl_handle = 0
        Ihl_debug = 0
        Ihl_mode = 1
        call hlsetdebug(Ihl_debug)
        call hlopen(HYDFIL,Ihl_mode,Ihl_handle,ierror)
        if(ierror.gt.0) then
          call hlgetlasterror(errstring)
          write(6,6000)ierror,errstring
          stop
        end if
        call hlsetlanguage(Ihl_handle,1,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        if(DTSSFAC.gt.0.) then
          call hlsetvartimestep(Ihl_handle,1,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
         else
           call hlsetvartimestep(Ihl_handle,0,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
        end if
        DES=' 666 '//char(0)
        call hladddescription(Ihl_handle,0,DES,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
!=======================================================================
!     Store the modeler name
!=======================================================================
        MOD='Created by: JL666 '//char(0)
        call hladddescription(Ihl_handle,1,MOD,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
!=======================================================================
!     Set the creator
!=======================================================================
        call hlsetcreator(Ihl_handle, 1, ierror)   !1 for fortran program
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        j=1461*(iyear+4800+int((imon-14)/12))/4+367*(imon-2-12*                                   &
         int((imon-14)/12))/12-(3*((iyear+4800+int((imon-14)/12)                                  &
         +100)/100))/4+iday-32075+tbegin+idays
	istartyear=100*(int(4*(j+68569)/146097)-49)+int(4000*                                           &
        (int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/                                 &
        1461001)+int(int(80*int(int(j+68569-(146097*int(4*int(j+68569)                            &
        /146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*                                        &
        int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
	istartmonth=int(80*int(int(j+68569-(146097*int(4*int(j+68569)/                                  &
        146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+                             &
        68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)+2-12*int(int(80*                             &
        int(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)                                  &
        -1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)                                    &
        /146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
	istartday=int(int(j+68569-(146097*int(4*int(j+68569)/146097)                                    &
        +3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)                              &
        /146097)+3)/4)+1)/1461001)/4+31)-2447*int(80*int(                                         &
        int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)                                      &
        -1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/                                   &
        146097)+3)/4)+1)/1461001)/4+31)/2447)/80
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

        call hlsetseedmoment(Ihl_handle,istartmonth,istartday,istartyear,istarthour,              &
                             istartminute,istartsecond,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
	call hlsetnumlayers(Ihl_handle,kc,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        NJUN=KC*(LCLT-2)
        NJUN=NJUN+1            
        if(njun.gt.nsg) then
        write(6,500)NJUN,NSG
        STOP
        END IF
500     FORMAT('THE NUMBER OF WASP SEGMENTS IN YOUR APPLICATION',I6,1x,                           &
               'IS GREATER THAN THE ARRAY DIMENSION:',I7)
        call hlsetnumsegments(Ihl_handle,njun,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if

        OPEN(94,FILE='settling.wsp',iostat=ios,STATUS='unknown')             
        IF(ios.eq.0) then                                                    
          write(94,*)'   Flow Function'                                       
          write(94,221)njun-1                                                
221       format(I12)                                                        
          do k=kc,2,-1                                                      
            do l=1,lclt-2                                                  
              l1=(lclt-2)*(kc-k)+l                                               
              l2=l1+(lclt-2)                                                
              write(94,223,err=222)l1,l2,1.0                         
            end do                                                  
          end do                                                 
          do l=1,lclt-2                                                
            l1=(lclt-2)*(kc-1)+l                                                   
            write(94,223,err=222)l1,njun,1.0                                
          end do                                                    
        end if                                                      
223     format(2I7,f10.2)                                            
222     CLOSE(94)                                                   
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
    	    I=I+1
	    LAUX(ILLT(LT),JLLT(LT),K)=I
	    l=lij(illt(lt),jllt(lt))
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
          END DO
        END DO
        i=i+1                                                     
        segname='bottom settling segment'                         
        call hlsetsegname(ihl_handle,i,segname,ierror)            
        NCHNH=0
        NCHNV=0
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          NCHNH=NCHNH+INT(SUBO(L))
          IF (IJCTLT(I+1,J).EQ.8) THEN
            IF (SUBO(L+1).EQ.1.) NCHNH=NCHNH+1
          END IF
          NCHNH=NCHNH+INT(SVBO(L))
          IF (IJCTLT(I,J+1).EQ.8) THEN
            IF (SVBO(LNC(L)).EQ.1.) NCHNH=NCHNH+1
          END IF
          NCHNV=NCHNV+INT(SWB(L))
        END DO
        NCHN=KC*NCHNH+(KC-1)*NCHNV

	NQ=NQSIJ
	DO L=1,NQSIJ
	  IF(LIJLT(IQS(L),JQS(L)).EQ.0) NQ=NQ-1
	END DO
        NCHN=NCHN+KC*NQ
!
        NQ=0                                        
	DO L=1,NQSIJ                                
            I=IQS(L)                                
            J=JQS(L)                                
	      IF(LIJLT(I,J).EQ.0) GOTO 1001         
              NS=NQSERQ(L)                          
	  DO K=1,KC                                 
	    IF(flagwaspbC(ns,k).EQ.1) NQ=NQ+1       
	  END DO                                    
1001    END DO                                      
        NCHN=NCHN-NQ
!
	NQ=NQCTL
	DO L=1,NQCTL
	  IF(LIJLT(IQCTLU(L),JQCTLU(L)).EQ.0) THEN
	    IF(LIJLT(IQCTLD(L),JQCTLD(L)).EQ.0) NQ=NQ-1
	  END IF
	END DO
        NCHN=NCHN+KC*NQ

	NQ=NQwr                                                   
	DO L=1,NQwr
	  IF(LIJLT(IQWRU(L),JQwrU(L)).EQ.0) THEN
	    IF(LIJLT(IQwrD(L),JQwrD(L)).EQ.0) NQ=NQ-1
	  END IF
	END DO
        NCHN=NCHN+NQ

      if(iswseda.gt.0) then      
        nchn=nchn+nlwseda*KC
      end if

        if(NCHN.gt.NF) then
        write(6,600)NCHN,NF
        STOP
        END IF
600     FORMAT('THE NUMBER OF WASP FLOWS IN YOUR APPLICATION',I6,1X,                              &
               'IS GREATER THAN THE ARRAY DIMENSION:',I7)
        call hlsetnumflowpaths(Ihl_handle, NCHN, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        OPEN(94,FILE='wasp_bound_seg.wsp',iostat=ios,STATUS='unknown') 
        write(94,*)'Number   wasp segment      I      J      K'        
        LCLTM2=LCLT-2
        LWASP=0
        lbound=0
        DO K=KC,1,-1
          KMUL=KC-K
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (SUBO(L).EQ.1.) THEN
              LWASP=LWASP+1
              LDTM=LT-1+KMUL*LCLTM2
              LUTM=LDTM-1
601         format(I7,I15,3I7)                                      
              IF (IJCTLT(I-1,J).EQ.8) then                          
                LUTM=0                                              
                lbound=lbound+1                                     
                write(94,601)lbound,ldtm,i,j,k                      
              end if                                                
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
            END IF

            IF (IJCTLT(I+1,J).EQ.8) THEN
              IF (SUBO(L+1).EQ.1.) THEN
	        LWASP=LWASP+1
                LDTM=0
                LUTM=LT-1+KMUL*LCLTM2
                lbound=lbound+1                                      
                write(94,601)lbound,lutm,i,j,k                       

               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
              END IF
            END IF
          END DO

          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (SVBO(L).EQ.1.) THEN
              LWASP=LWASP+1
              LSLT=LSCLT(LT)
              LDTM=LT-1+KMUL*LCLTM2
              LUTM=LSLT-1+KMUL*LCLTM2
              IF(IJCTLT(I,J-1).EQ.8) then
                LUTM=0
                lbound=lbound+1                                    
                write(94,601)lbound,ldtm,i,j,k                     
               end if
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,2,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
            END IF
            IF (IJCTLT(I,J+1).EQ.8) THEN
              LN=LNC(L)
              IF (SVBO(LN).EQ.1.) THEN
                LWASP=LWASP+1
                LSLT=LSCLT(LT)
                LDTM=0
                LUTM=LT-1+KMUL*LCLTM2
                lbound=lbound+1                                       
                write(94,601)lbound,lutm,i,j,k                        
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,2,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
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
	      IF(LIJLT(I,J).EQ.0) GOTO 100
              NS=NQSERQ(Lt)                                          
                IF(flagwaspbC(ns,k).EQ.1) GOTO 100                   
            LWASP=LWASP+1
            LDTM=LAUX(I,J,K)
            LUTM=0
                lbound=lbound+1                                   
                write(94,601)lbound,ldtm,i,j,k                    
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
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
              if(lutm.eq.0) then
                lbound=lbound+1                                   
                write(94,601)lbound,ldtm,IQCTLD(LT),JQCTLD(LT),k  
              end if
              if(ldtm.eq.0) then
                lbound=lbound+1                                   
                write(94,601)lbound,lutm,IQCTLu(LT),JQCTLu(LT),k  
              end if
              LWASP=LWASP+1
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
200       END DO
        END DO

        DO LT=1,NQwr                                  
          I=IQwrU(LT)
          J=JQwrU(LT)
          k=KQwrU(LT)
          LUTM=LAUX(I,J,K)
          I=IQwrD(LT)
          J=JQwrD(LT)
          K=KQwrD(LT)
          LDTM=LAUX(I,J,K)
          IF(LUTM.EQ.0.AND.LDTM.EQ.0) GOTO 201
          LWASP=LWASP+1
          call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
201     END DO

        IF(ISWSEDA.GT.0) THEN      
         DO K=KC,1,-1
          DO LT=1,NLWSEDA
            LUTM=0
            I=ICWSEDA(LT)
            J=ICWSEDA(LT)
            LDTM=LAUX(I,J,K)

            lbound=lbound+1                                  
            write(94,601)lbound,ldtm,I,J,k                   

              LWASP=LWASP+1
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
          END DO
         END DO
        END IF

        IF(KC.GT.1)THEN
          DO K=KS,1,-1
            KMUL1=KS-K
            KMUL2=KMUL1+1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF (SWB(L).EQ.1.) THEN
                LWASP=LWASP+1
                LUTM=LT-1+KMUL1*LCLTM2
                LDTM=LT-1+KMUL2*LCLTM2
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,3,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
              END IF
            END DO
          END DO
        ENDIF
	inumsegconsts=3  !volume,depth,velocity
	if(istran(2).GE.1) inumsegconsts=4   !water temperature modeled in EFDC and transfered   
	if(istran(1).GE.1) inumsegconsts=5   !salinity modeled in EFDC and transfered
	inumsegconsts=6                      !TKE DISSIPATION RATE, 12/22/2012, DRBC LZ      !HNR_GHD -  hotwired to 6

        call hlsetnumsegconsts(Ihl_handle, inumsegconsts, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        call hlsetnumfpconsts(Ihl_handle, 3, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        call hlsetsegconsttype(Ihl_handle, 1, 0, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        call hlsetsegconsttype(Ihl_handle, 2, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        call hlsetsegconsttype(Ihl_handle, 3, 2, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        if(istran(2).GE.1) then
        call hlsetsegconsttype(Ihl_handle, 4, 3, ierror)
        if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        end if
        if(istran(1).GE.1) then                  ! 05/08/2019 DRBC LZ, OUTPUT SALINITY
         call hlsetsegconsttype(Ihl_handle, 5, 4, ierror)
         if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,23         
          stop
         end if
        end if   
      
        call hlsetsegconsttype(Ihl_handle, 6, 5, ierror)    !HNR_GHD 4/2022   added setting dipsrate as a constituent
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,235
            stop
          end if
      
        call hlsetfpconsttype(Ihl_handle, 1, 0, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        call hlsetfpconsttype(Ihl_handle, 1, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        call hlsetfpconsttype(Ihl_handle, 1, 2, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
      if(DTSSFAC.gt.0.) then
      call hlsettimestep(Ihl_handle,dtdyn,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        else
	call hlsethydtimestep(Ihl_handle,dt,ierror)
      call hlsettimestep(Ihl_handle,dt,ierror)   
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        endif
      if(DTSSFAC.gt.0.) then                      
      call hlsettimestep(Ihl_handle,dtdyn,ierror) 
        if(ierror .gt. 0)then                     
          call hlgetlasterror(errstring)          
          write(6,6000) ierror, errstring         
          stop                                    
        end if                                    
        else                                      
      call hlsethydtimestep(Ihl_handle,dt,ierror) 
        if(ierror .gt. 0)then                     
          call hlgetlasterror(errstring)          
          write(6,6000) ierror, errstring         
          stop                                    
        end if                                    
      endif                                       
      if(DTSSFAC.gt.0.) then                                    
      else                                                      
        rinterval=(dt*NTSMMT)/86400.                            
      call hlsetupdateint(Ihl_handle,rinterval,ierror)          
        if(ierror .gt. 0)then                                   
          call hlgetlasterror(errstring)                        
          write(6,6000) ierror, errstring                       
          stop                                                  
        end if                                                  
      end if                                                    
	call hlsethydtowaspratio(Ihl_handle,NTSMMT,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
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
              LWASP=LWASP+1
              SegVel(LWASP)=0.0
              
              IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                SegDepth(LWASP)=HP(L)*DZC(K)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=SAL(L,K)
                SegTemp(LWASP)=TEM(L,K)
                SegDispr(LWASP)=DISPRATE(L,K)        
              ELSE
                SegDepth(LWASP)=HDRY*DZC(K)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=0.0              
                SegTemp(LWASP)=TATMT(L)
                SegDispr(LWASP)=0.0                  !DRBC LZ 12/22/2020
              END IF
              
              IF(NTSMMT.LT.NTSPTC) THEN
                  
                IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                  SegDepth(LWASP)=HMP(L)*DZC(K)
                  SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                  SegSalt(LWASP)=SALINIT(L,K)
                  SegTemp(LWASP)=TEMINIT(L,K)
                  SegDispr(LWASP)=DISPRATE(L,K)        !DRBC LZ 12/22/2020
                ELSE
                  SegDepth(LWASP)=HDRY*DZC(K)
                  SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                  SegSalt(LWASP)=0.0
                  SegTemp(LWASP)=TATMT(L)
                  SegDispr(LWASP)=0.0                  !DRBC LZ 12/22/2020
                END IF
                
              END IF
            END DO
          END DO

          SegVel(LWASP+1)=0.0                 
          SegDepth(LWASP+1)=10.0
          SegVolume(LWASP+1)=1.0e12
          SegSalt(LWASP+1)=0.00
          SegTemp(LWASP+1)=SegTemp(LWASP)
          SegDispr(LWASP+1)=0.0

          DO I=1,NCHN
	    FLOW(I)=0.0
	    CRNU(I)=0.0
	    BRINTT(I)=0.0
          END DO
	  call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
          if(istran(2).GE.1) then
	    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              stop
            end if
          end if
          if(istran(1).GE.1) then                  ! 05/08/2019 DRBC LZ, OUTPUT INITIAL SALINITY
	    call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,35
              stop
            end if
          end if
	  call hlsetseginfo(Ihl_handle,6,SegDispr,ierror)   !DRBC LZ, 12/22/2020, OUTPUT TKE DISSIPATION RATE
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,351
            stop
          end if           
	  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlmomentcomplete(Ihl_Handle,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
          ELSE
          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              LN=LNC(L)
              LWASP=LWASP+1
              VELX=0.5*(U(L,K)+U(L+1,K))
              VELY=0.5*(V(L,K)+V(LN,K))
              VELZ=0.5*(W(L,K-1)+W(L,K))
              VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
              
              IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
                SegVel(LWASP)=VELMAG
                SegDepth(LWASP)=HP(L)*DZC(K)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=SAL(L,K)
                SegTemp(LWASP)=TEM(L,K)
                SegDispr(LWASP)=DISPRATE(L,K)        !DRBC LZ 12/22/2020
              ELSE
                SegVel(LWASP)=0.0
                SegDepth(LWASP)=HDRY*DZC(K)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=0.0
                SegTemp(LWASP)=TATMT(L)
                SegDispr(LWASP)=0.0                  !DRBC LZ 12/22/2020
              END IF
              
            END DO
          END DO
          SegVel(LWASP+1)=0.0       
          SegDepth(LWASP+1)=10.0
          SegVolume(LWASP+1)=1.0e12
          SegSalt(LWASP+1)=0.00
          SegTemp(LWASP+1)=SegTemp(LWASP)
          SegDispr(LWASP+1)=0.0

          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              ADDLW=0.0
              IF (SUB(L).EQ.1.) THEN
                LW=L-1
                ADDLW=DYU(L)*0.5*(AH(L,K)+AH(L-1,K))*DZC(K)*0.5*(HP(L)+HP(LW))*DXIU (L)
                  vol1=DXYP(L)*HP(L)*DZC(K)         
                  vol2=DXYP(LW)*HP(LW)*DZC(K)       
                  if (vol1.LT.vol2) vol2=vol1                    
                  ad=vol2/dt*adcoeff                                     
                  if (addlw.GT.ad) then                          
                    addlw=ad                             
                  end if                                         
                  if(ishdmf.lt.2) then                           
                    addlw=0.                                     
                  end if                                         

              END IF
              IF (SUBO(L).EQ.1.) THEN
                LWASP=LWASP+1
                
                IF((ISCDRY(L).EQ.1).OR.(ISCDRY(L-1).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
                  FLOW(LWASP)=0.0
	            CRNU(LWASP)=0.0
  	            BRINTT(LWASP)=0.0
                ELSE
                  FLOW(LWASP)=UHDY(L,K)*DZC(K)
	            CRNU(LWASP)=2.*UHDY(L,K)*DYIU(L)*DXIU(L)/(HP(L)+HP(L-1))
  	            BRINTT(LWASP)=ADDLW
                END IF
                
              END IF
              IF (IJCTLT(I+1,J).EQ.8) THEN
                IF (SUBO(L+1).EQ.1.) THEN
                  LWASP=LWASP+1
                  
                IF((ISCDRY(L).EQ.1).OR.(ISCDRY(L+1).EQ.1)) THEN    !HNR_GHD 6/2019 W/D
	            FLOW(LWASP)=0.0
	            CRNU(LWASP)=0.0
                  BRINTT(LWASP)=0.0
                ELSE
	            FLOW(LWASP)=UHDY(L+1,K)*DZC(K)
	            CRNU(LWASP)=2.*UHDY(L+1,K)*DYIU(L+1)*DXIU(L+1)/(HP(L)+HP(L+1))
                  BRINTT(LWASP)=0.0
                END IF
                
                END IF
              END IF
            END DO
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              ADDLS=0.0
              LS=LSC(L)
              IF (SVB(L).EQ.1.) THEN
                ADDLS=DXV(L)*0.5*(AH(L,K)+AH(LS,K))*DZC(K)*0.5*(HP(L) +HP(LS))*DYIV (L)
                  vol1=DXYP(L)*HP(L)*DZC(K)                   
                  vol2=DXYP(Ls)*HP(Ls)*DZC(K)                 
                  if (vol1.LT.vol2) vol2=vol1                   
                  ad=vol2/dt*adcoeff                                    
                  if (addls.GT.ad) then                         
                    addls=ad                                
                  end if                                        
                  if(ishdmf.lt.2) then                           
                    addls=0.                                     
                  end if                                         
              END IF
              IF (SVBO(L).EQ.1.) THEN
	        LWASP=LWASP+1
              
                IF((ISCDRY(L).EQ.1).OR.(ISCDRY(LS).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
                  FLOW(LWASP)=0.0
	            CRNU(LWASP)=0.0
                  BRINTT(LWASP)=0.0
                ELSE
                  FLOW(LWASP)=VHDX(L,K)*DZC(K)
	            CRNU(LWASP)=2.*VHDX(L,K)*DYIV(L)*DXIV(L)/(HP(L)+HP(LS))
                  BRINTT(LWASP)=ADDLS
                END IF
                
              END IF
              IF (IJCTLT(I,J+1).EQ.8) THEN
                LN=LNC(L)
                IF (SVBO(LN).EQ.1.) THEN
  	          LWASP=LWASP+1
                
                IF((ISCDRY(L).EQ.1).OR.(ISCDRY(LN).EQ.1)) THEN                !HNR_GHD 6/2019 W/D
	            FLOW(LWASP)=0.0
                  CRNU(LWASP)=0.0
	            BRINTT(LWASP)=0.0                                
                ELSE
	            FLOW(LWASP)=VHDX(LN,K)*DZC(K)
                  CRNU(LWASP)=2.*VHDX(LN,K)*DYIV(LN)*DXIV(LN)/(HP(L)+HP(LN))
	            BRINTT(LWASP)=addls                                
                END IF
                
                END IF
              END IF
            END DO
          END DO

          DO K=KC,1,-1
            DO LT=1,nqsij
              I=Iqs(LT)
              J=Jqs(LT)
              IF(LIJLT(I,J).EQ.0) GOTO 310
              NS=NQSERQ(Lt)
              L=LQS(Lt)
                IF(flagwaspbC(ns,k).EQ.1) GOTO 310                       
              LWASP=LWASP+1
              FLOW(LWASP)=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
              CRNU(LWASP)=flow(LWASP)/DXp(L)/dyp(l)/(dzc(k)*HP(L))
              BRINTT(LWASP)=0.0
310         END DO
          END DO
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
              UDDXTMP=flowx/DXp(Lu)/dyp(lu)/(dzc(k)*HP(Lu))
              IF(iu.eq.id) THEN
                addls=dxv(lu)*ah(lu,k)*dzc(k)*0.5*(hp(lu)+hp(ld))*dyiv(lu)
              ELSE
                addls=dyu(lu)*ah(lu,k)*dzc(k)*0.5*(hp(lu)+hp(ld))*dxiu(lu)
              END IF
                  if(ishdmf.lt.2) then                           
                    addls=0.                                     
                  end if                                         
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWX
	      CRNU(LWASP)=UDDXTMP
              BRINTT(LWASP)=ADDLS
410         END DO
          END DO

            DO LT=1,nqwr
              Iu=Iqwru(LT)
              Ju=Jqwru(LT)
              Ku=Kqwru(LT)
              Lu=Lij(iu,ju)
              Id=Iqwrd(LT)
              Jd=Jqwrd(LT)
              Kd=Kqwrd(LT)
              Ld=Lij(id,jd)
              NS=NQWRSERQ(lt)
              IF(LU.EQ.0.AND.LD.EQ.0) GOTO 411
              flowx=QWR(lt)+QWRSERT(NS)
              UDDXTMP=flowx/DXp(Lu)/dyp(lu)/(dzc(ku)*HP(Lu))
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWX
	      CRNU(LWASP)=UDDXTMP
            BRINTT(LWASP)=0.0
411         END DO

!  Advection and dispersion in assimilation flows
        IF(ISWSEDA.GT.0) THEN     
         DO K=KC,1,-1
          DO LT=1,NLWSEDA
            I=ICWSEDA(LT)
            J=ICWSEDA(LT)
            L=LIJ(I,J)
            LWASP=LWASP+1
            FLOW(LWASP)=QWSEDA(L,K)
            CRNU(LWASP)=flow(LWASP)/DXp(L)/dyp(l)/(dzc(k)*HP(L))
            BRINTT(LWASP)=0.0
          END DO
         END DO
        END IF

!
! Advection and dispersion in the Z-direction:
!
          IF (KC.GT.1) THEN
            DO K=KS,1,-1
              DO LT=2,LALT
                I=ILLT(LT)
                J=JLLT(LT)
                L=LIJ(I,J)
                addl=0.0
                ADDL1=ab(l,k)*hp(l) 
                if (addl1.gt.abwmx(l)) then
                  addl1=abwmx(l)         
                end if
                IF (SPB(L).EQ.1.) THEN
                  ADDL=DXYP(L)*addl1/hp(l)*DZIG(K)
                  vol1=DXYP(L)*HP(L)*DZC(K)                 
                  vol2=DXYP(L)*HP(L)*DZC(K+1)               
                  if (vol1.LT.vol2) vol2=vol1                 
                  ad=adcoeff*vol2/dt                                  
                  if (addl.GT.ad) then                        
                    addl=ad                               
                  end if                                      
                END IF
                IF (SWB(L).EQ.1) THEN
                  LWASP=LWASP+1

                IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
  	            FLOW(LWASP)=-DXYP(L)*W(L,K)
                  CRNU(LWASP)=W(L,K)*DZIG(K)/HP(L)
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
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
! next only if temperature is transfered
          if(istran(2).GE.1) then
	    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              stop
            end if
          end if
          if(istran(1).GE.1) then                  ! 06/28/2019 DRBC LZ, OUTPUT INITIAL SALINITY DURING A HOT-START
            call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
            if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              stop
            end if
          end if           
          call hlsetseginfo(Ihl_handle,6,SegDispr,ierror)   !DRBC LZ, 12/22/2020, OUTPUT TKE DISSIPATION RATE
	  if(ierror .gt. 0)then
	    call hlgetlasterror(errstring)
	    write(6,6000) ierror, errstring,441
	    stop
          end if           
	  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlmomentcomplete(Ihl_Handle,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
        END IF
!      END INITIAL CONDITIONS WHEN IDAYS=0
       END IF
      GOTO 3000
      END IF
      IF(N.LT.IBEGIN) GOTO 3000
!
!----------------------------------------------------------------------C
!
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
          IF (SUB(L).EQ.1.) THEN
            LW=L-1
            ADDLW=DYU(L)*AHULPF(L,K)*DZC(K)*0.5*(HLPF(L) +HLPF(LW))*DXIU (L)
            vol1=DXYP(L)*HLPF(L)*DZC(K)                           
            vol2=DXYP(LW)*HLPF(LW)*DZC(K)                         
            if (vol1.LT.vol2) vol2=vol1                           
            ad=adcoeff*vol2/dt                                            
            if (addlw.GT.ad) then                                 
              addlw=ad                                        
            end if                                                
                  if(ishdmf.lt.2) then                           
                    addlw=0.                                     
                  end if                                         
          END IF
          IF (SUBO(L).EQ.1.) THEN
            TMPVAL=UHLPF(L,K)+SVPT*UVPT(L,K)
            FLOWX=DYU(L)*TMPVAL*DZC(K)
            UDDXTMP=2.*TMPVAL*DXIU(L)/(HLPF(L)+HLPF(L-1))
            LWASP=LWASP+1
            FLOW(LWASP)=FLOWX
            CRNU(LWASP)=UDDXTMP
            BRINTT(LWASP)=ADDLW
          END IF
          IF (IJCTLT(I+1,J).EQ.8) THEN
            IF (SUBO(L+1).EQ.1.) THEN
              TMPVAL=UHLPF(L+1,K)+SVPT*UVPT(L+1,K)
              FLOWX=DYU(L+1)*TMPVAL*DZC(K)
              UDDXTMP=2.*TMPVAL*DXIU(L+1)/(HLPF(L+1)+HLPF(L))
              IPTMP=I+1
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
          IF (SVB(L).EQ.1.) THEN
            LS=LSC(L)
            ADDLS=DXV(L)*AHVLPF(L,K)*DZC(K)*0.5*(HLPF(L) +HLPF(LS))*DYIV (L)
            vol1=DXYP(L)*HLPF(L)*DZC(K)                          
            vol2=DXYP(LS)*HLPF(LS)*DZC(K)                        
            if (vol1.LT.vol2) vol2=vol1                          
            ad=adcoeff*vol2/dt                                           
            if (addls.GT.ad) then                                
              addls=ad                                       
            end if                                               
                  if(ishdmf.lt.2) then                           
                    addls=0.                                     
                  end if                                         
          END IF
          IF (SVBO(L).EQ.1.) THEN
            TMPVAL=VHLPF(L,K)+SVPT*VVPT(L,K)
            FLOWY=DXV(L)*TMPVAL*DZC(K)
            VDDYTMP=2.*TMPVAL*DYIV(L)/(HLPF(L)+HLPF(LSC(L)))
            JMTMP=J-1
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWY
	    CRNU(LWASP)=VDDYTMP
	    BRINTT(LWASP)=ADDLS
          END IF
          IF (IJCTLT(I,J+1).EQ.8) THEN
            LN=LNC(L)
            IF (SVBO(LN).EQ.1.) THEN
              TMPVAL=VHLPF(LN,K)+SVPT*VVPT(LN,K)
              FLOWY=DXV(LN)*TMPVAL*DZC(K)
              VDDYTMP=2.*TMPVAL*DYIV(LN)/(HLPF(LN)+HLPF(L))
              JPTMP=J+1
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWY
	    CRNU(LWASP)=VDDYTMP
	    BRINTT(LWASP)=ADDLS
            END IF
          END IF
        END DO
      END DO

! Advection and dispersion in input flows
      DO K=KC,1,-1
        DO LT=1,nqsij
          I=Iqs(LT)
          J=Jqs(LT)
          IF(LIJLT(I,J).EQ.0) GOTO 300
          NS=NQSERQ(Lt)
          L=LQS(Lt)
            IF(flagwaspbC(ns,k).EQ.1) GOTO 300                       
          flowx=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
          UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*HmP(L))
          addlw=0.0                                        
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWX
	    CRNU(LWASP)=UDDXTMP
	    BRINTT(LWASP)=ADDLW
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
          UDDXTMP=flowx/DXp(Lu)/dyp(lu)/(dzc(k)*HmP(Lu))
          IF(iu.eq.id) THEN
            addls=dxv(lu)*ahvlpf(lu,k)*dzc(k)*0.5*(hlpf(lu)+hlpf(ld))*dyiv(lu)
          ELSE
            addls=dyu(lu)*ahulpf(lu,k)*dzc(k)*0.5*(hlpf(lu)+hlpf(ld))*dxiu(lu)
          END IF
                  if(ishdmf.lt.2) then                           
                    addls=0.                                     
                  end if                                         
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWX
	    CRNU(LWASP)=UDDXTMP
	    BRINTT(LWASP)=ADDLs
400     END DO
      END DO

!  Advection and dispersion in w/r flows       
        DO LT=1,nqwr
          Iu=Iqwru(LT)
          Ju=Jqwru(LT)
          Ku=Kqwrd(LT)
          Lu=Lij(iu,ju)
          Id=Iqwrd(LT)
          Jd=Jqwrd(LT)
          Kd=Kqwrd(LT)
          Ld=Lij(id,jd)
          NS=NQWRSERQ(lt)
          IF(LU.EQ.0.AND.LD.EQ.0) GOTO 405
          flowx=QWR(lt)+QWRSERT(NS)
          UDDXTMP=flowx/DXp(Lu)/dyp(lu)/(dzc(ku)*HmP(Lu))
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWX
	    CRNU(LWASP)=UDDXTMP
	    BRINTT(LWASP)=0.0
405     END DO

        IF(ISWSEDA.GT.0) THEN      
         DO K=KC,1,-1
          DO LT=1,NLWSEDA
            I=ICWSEDA(LT)
            J=ICWSEDA(LT)
            L=LIJ(I,J)
            LWASP=LWASP+1
            FLOW(LWASP)=QWSEDA(L,K)
            CRNU(LWASP)=flow(LWASP)/DXp(L)/dyp(l)/(dzc(k)*HP(L))
            BRINTT(LWASP)=0.0
          END DO
         END DO
        END IF
!
! Advection and dispersion in the Z-direction:
!
      IF (KC.GT.1) THEN
        DO K=KS,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            ADDL=0.0
            addl1=ablpf(l,k)*hlpf(l)  
            if(addl1.gt.abwmx(l)) then
              addl1=abwmx(l)   
            end if
            IF (SPB(L).EQ.1.) THEN
              ADDL=DXYP(L)*Addl1/HLPF(L)*DZIG(K)
              vol1=DXYP(L)*HLPF(L)*DZC(K)           
              vol2=DXYP(L)*HLPF(L)*DZC(K+1)         
              if (vol1.LT.vol2) vol2=vol1           
              ad=adcoeff*vol2/dt                            
              if (addl.GT.ad) then                  
                addl=ad                         
              end if                                
            END IF
            IF (SWB(L).EQ.1) THEN
              TMPVAL=WLPF(L,K)+SVPT*WVPT(L,K)
              FLOWZ=-DXYP(L)*TMPVAL
              WDDZTMP=TMPVAL*DZIG(K)/HLPF(L)
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
          LCELTMP=LCELTMP+1
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          LN=LNC(L)
          VOLUM=DXYP(L)*HLPF(L)*DZC(K)

          IF(NTSMMT.LT.NTSPTC) THEN
              
            IF(ISCDRY(L).EQ.0) THEN                !HNR_GHD 6/2019 W/D
              VOLUM=DXYP(L)*HP(L)*DZC(K)
            ELSE
              VOLUM=DXYP(L)*HDRY*DZC(K)
            END IF
          END IF
          
          DEPTH=HLPF(L)*DZC(K)
          VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K) +UHLPF(L+1,K)+SVPT*UVPT(L+1,K))/HLPF(L)
          VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K) +VHLPF(LN,K)+SVPT*VVPT(LN,K) )/HLPF(L)
          VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1) +WLPF(L,K)+SVPT*WVPT(L,K) )
          VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
            SegSalt(LCELTMP)=SALLPF(L,K)
            SegTemp(LCELTMP)=TEMLPF(L,K)
            SegDispr(LCELTMP)=DISPRLPF(L,K)        !DRBC LZ 12/22/2020 
	    SegVolume(LCELTMP)=VOLUM
	    SegDepth(LCELTMP)=DEPTH
	    SegVel(LCELTMP)=VELMAG
        END DO
      END DO

          SegVel(LCELTMP+1)=0.0               
          SegDepth(LCELTMP+1)=10.0
          SegVolume(LCELTMP+1)=1.0e12
          SegSalt(LCELTMP+1)=0.00
          SegTemp(Lceltmp+1)=SegTemp(Lceltmp)
          SegDispr(LCELTMP+1)=0.0

!
         if(DTSSFAC.gt.0.) then
             call hlsettimestep(Ihl_handle,dtdyn,ierror)    
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
         endif

	   call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
       if(istran(2).GE.1) then
	   call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
       end if
       if(istran(1).GE.1) then                  ! 06/13/2019 DRBC LZ, OUTPUT TEMPORAL SALINITY
           call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
           if(ierror .gt. 0)then
             call hlgetlasterror(errstring)
             write(6,6000) ierror, errstring,35
             stop
           end if
       end if        
       call hlsetseginfo(Ihl_handle,6,SegDispr,ierror)   !DRBC LZ, 12/22/2020, OUTPUT TKE DISSIPATION RATE
       if(ierror .gt. 0)then
         call hlgetlasterror(errstring)
         write(6,6000) ierror, errstring,531
         stop
       end if        
	   call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlmomentcomplete(Ihl_Handle,ierror)
         if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if

!
3000  JSWASP=0
6000  format('Error ',I10, ' : ', A)
      RETURN
      END