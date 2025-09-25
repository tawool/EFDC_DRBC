      SUBROUTINE INIT_BMD2
!**********************************************************************C
!
! **  CREATED BY HUGO N RODRIGUEZ ON june 2019 FOR BMD2 OUTPUT
!
!
! **  SUBROUTINE INIT_BMD2 INITIALIZES THE BMD2 FILE
!**********************************************************************C
!
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      INCLUDE 'allset.int'
      PARAMETER(IKMAX=25)
!
      CHARACTER*99 ITEXT, JTEXT, KTEXT
      CHARACTER*30 A 
      CHARACTER*99 DMPFIL,DMPFIL1 
      CHARACTER*99 D(IKMAX),Uni(IKMAX),SEGNAME 
      INTEGER OPTFLAG
    
      OPEN(94,FILE='BMDOUTPUT.CTL',IOSTAT=IOS,STATUS='OLD')
      IF(IOS.NE.0) RETURN
      DO I=1,4
        READ(94,*)
      END DO
      READ(94,*)BMDFLAG
      
      IF(BMDFLAG.EQ.0) THEN
        CLOSE(94)
        RETURN
      ENDIF
       DO I=1,6
        READ(94,*)
      END DO
      
      READ(94,*)DTBMD,TIBMD,TEBMD,DMPFIL
      NSBMD=NINT(DTBMD*3600.0/TCON*NTSPTC)
      NCBMD=1
        D(01) = 'Surface Elev'//char(0)
        D(02) = 'Water Depth'//char(0)
        D(03) = 'Flush Time'//char(0)
        D(04) = 'Courant Max'//char(0)
        D(05) = 'Shear Stress'//char(0)
        D(06) = 'Layer Depth'//char(0)
        D(07) = 'Salt'//char(0)
        D(08) = 'Temp'//char(0)
        D(09) = 'Dye'//char(0)
        D(10) = 'TSS'//char(0)
        D(11) = 'Vel U'//char(0)
        D(12) = 'Vel V'//char(0)
        D(13) = 'Vel East'//char(0)
        D(14) = 'Vel North'//char(0)
        D(15) = 'Speed'//char(0)
        D(16) = 'Flow U'//char(0)
        D(17) = 'Flow V'//char(0)
        D(18) = 'Courant U'//char(0)
        D(19) = 'Courant V'//char(0)
        D(20) = 'Courant Coriolis'//char(0)
        D(21) = 'Volume'//char(0)
        D(22) = 'Vel W'//char(0)
        D(23) = 'Flow W'//char(0)
        D(24) = 'Courant W'//char(0)
        D(25) = 'Vert Diffusion Coef'//char(0)

        UNI(01) = 'meters'//char(0)
        UNI(02) = 'meters'//char(0)
        UNI(03) = 'hours'//char(0)
        UNI(04) = ''//char(0)
        UNI(05) = 'Pa'//char(0)
        UNI(06) = 'meters'//char(0)
        UNI(07) = 'PSU'//char(0)
        UNI(08) = 'C'//char(0)
        UNI(09) = 'Conc.'//char(0)
        UNI(10) = 'mg/L'//char(0)
        UNI(11) = 'cm/s'//char(0)
        UNI(12) = 'cm/s'//char(0)
        UNI(13) = 'cm/s'//char(0)
        UNI(14) = 'cm/s'//char(0)
        UNI(15) = 'cm/s'//char(0)
        UNI(16) = 'm3/s'//char(0)
        UNI(17) = 'm3/s'//char(0)
        UNI(18) = ''//char(0)
        UNI(19) = ''//char(0)
        UNI(20) = ''//char(0)
        UNI(21) = 'm3'//char(0)
        UNI(22) = 'cm/s'//char(0)
        UNI(23) = 'm3/s'//char(0)
        UNI(24) = ''//char(0)
        UNI(25) = 'm2/s'//char(0)

       DO I=1,10
         READ(94,*)
       END DO
     
       READ(94,*)IK
      
      if(ik.ne.IKMAX) then
        write(*,*)'Error in file bmdoutput.ctl. Number of BMD2 output variables:',ik,             &
        ' different than maximum for this compilation: ', ikmax, 'Program will stop'
        WRITE(*,*)'The variables in this compilation are:'
        write(*,*)'           K   VARIABLE'
        DO K=1,IKMAX
          WRITE(*,*)K,D(K)
        END DO
        stop
      endif
      
       DO KK=1,IKMAX
         READ(94,*)I,A,JJF(I)
       END DO

	if(kc.eq.1) then
	  do i=22,25
	    jjf(i)=0
	  end do
      endif

      IF(ISTRAN(1).eq.0) jjf(7)=0                            ! dont output salinity if it is not modeled even if the user inadvertently tries to output it
      IF(ISTRAN(2).eq.0) jjf(8)=0                            ! dont output temperature if it is not modeled even if the user inadvertently tries to output it
      IF(ISTRAN(3).eq.0) jjf(9)=0                            ! dont output dye if it is not modeled even if the user inadvertently tries to output it
      IF(ISTRAN(6).eq.0.and.ISTRAN(7).eq.0) jjf(10)=0        ! dont output TSS if sediment is not modeled even if the user inadvertently tries to output it
     
       
      DO I=1,3
        READ(94,*)
      END DO

      READ(94,*)IBMDSEG  
      
      IF(IBMDSEG.GT.0) THEN  
        DO KK=1,IBMDSEG
          READ(94,*)I,J
          ISEGIN(KK)=LIJ(I,J)
          IF(ISEGIN(KK).EQ.0) THEN
            WRITE(0,10)I,J
10          format('Segment ',i5,',',i5,' invalid, program will stop')
            WRITE(0,*)
            STOP
          ENDIF
        END DO
      ELSE
        IBMDSEG= (LA-1)
        DO KK=1,IBMDSEG
          ISEGIN(KK)=KK+1
        END DO
      END IF
          
      NUMVAR=0
      DO I=1,IKMAX
	  IF(JJF(I).EQ.1)NUMVAR=NUMVAR+1
      ENDDO
      IF(NUMVAR.EQ.0) THEN
        JJF(1)=1
        NUMVAR=1
      ENDIF
	
      NUMSEG=0
      DO I=1,IBMDSEG 
        IF(IGRIDV.EQ.0) THEN
          NUMSEG=NUMSEG+KC 
        ELSE
          NUMSEG=NUMSEG+KC-KGVCP(ISEGIN(I))+1 
        ENDIF
      ENDDO
    !------------------------------------------------------------------------------------------------------------
    ! Set BMD2 Performance Options
    !------------------------------------------------------------------------------------------------------------
      
      CALL BMD2SETLOGLEVELTHRESHOLD(0,IERROR)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
      
      Call bmd2setproperty("MAXVIEWBUFFERSIZE"//char(0),'8GB',IERROR)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
      Call bmd2setproperty("MAXTHREADS"//char(0),'4',IERROR)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
      Call bmd2setproperty("MAXTVSVIEWFRAMES"//char(0),'700',IERROR)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
      Call bmd2setproperty("MINTVSVIEWFRAMES"//char(0),'600',IERROR)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
      Call bmd2setproperty("MAXVSTVIEWFRAMES"//char(0),'1',IERROR)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF

      CALL BMD2SETINDEXBASE(1,IERROR)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
      
      DMPFIL1 = TRIM(DMPFIL)//'.BMD2'//CHAR(0)

      NumBMD2Writes = (NTC*TIDALP)/3600./DTBMD
      BMD_WRITE=0
      IHOUR=0
      IMIN=0
      ISEC=0
               
      OPTFLAG=2

      CALL BMD2CREATE(DMPFIL1,NUMSEG,NUMVAR,NumBMD2Writes,                   &
                      IYEAR,IMON,IDAY,IHOUR,IMIN,ISEC,OPTFLAG,HANDLEBMD2,IERROR) 
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
                
      WRITE(*,*)'BMD2 FILE SUCCESSFULLY CREATED.'
      WRITE(*,*)
      
      IVAR=0
      DO I=1,IKMAX
        IF(JJF(I).EQ.1) THEN
          IVAR=IVAR+1
          CALL BMD2SETVARIABLEMETADATA(HANDLEBMD2,IVAR,"NAME"//char(0),D(I),IERROR)
          IF(iError .gt. 0)Then
            Call bmd2getlasterror(80,errMessage)
            Write(*,*)errMessage
            Stop
          End IF
          CALL BMD2SETVARIABLEMETADATA(HANDLEBMD2,IVAR,"UNITS"//char(0),UNI(I),IERROR)
          IF(iError .gt. 0)Then
            Call bmd2getlasterror(80,errMessage)
            Write(*,*)errMessage
            Stop
          End IF
        ENDIF
      END DO 

      ISEG=0
      DO KK=1,ibmdseg
        l=isegin(kk)
        DO K=1,KC
          ! segment name I=0001,J=0001,K=001 I=9999,J=9999,K=999. KEEP FORMAT USED BY LCW SO WE DONT DISRUPT WITH GRAPH - COVERS MOST CASES ALTHOUGH SOME APPLICANTION HAVE J LARGER THAN 9999 
          IF (.NOT. (IGRIDV.EQ.1 .AND. K.LT.KGVCP(L))) THEN
            ISEG=ISEG+1
            WRITE(SEGNAME,40) IL(L),JL(L),K
40          FORMAT('I=',I4.4,',J=',I4.4,',K=',I3.3)
                          
            SEGNAME=TRIM(SEGNAME)//CHAR(0)
            CALL BMD2SETSEGMENTMETADATA(HANDLEBMD2,ISEG,"NAME"//char(0),SEGNAME,IERROR)
          IF(iError .gt. 0)Then
            Call bmd2getlasterror(80,errMessage)
            Write(*,*)errMessage
            Stop
          End IF
             
           ! ADDED METADATA FOR SEGMENT: I,J,K 
            WRITE(ITEXT,50)IL(L)
            ITEXT=TRIM(ITEXT)//CHAR(0)
            WRITE(JTEXT,50)JL(L)
            JTEXT=TRIM(JTEXT)//CHAR(0)
            WRITE(KTEXT,50)K
            KTEXT=TRIM(KTEXT)//CHAR(0)
50          FORMAT(I9)
            CALL BMD2SETSEGMENTMETADATA(HANDLEBMD2,ISEG,"I"//CHAR(0),ITEXT,IERROR)
          IF(iError .gt. 0)Then
            Call bmd2getlasterror(80,errMessage)
            Write(*,*)errMessage
            Stop
          End IF
            CALL BMD2SETSEGMENTMETADATA(HANDLEBMD2,ISEG,"J"//CHAR(0),JTEXT,IERROR)
          IF(iError .gt. 0)Then
            Call bmd2getlasterror(80,errMessage)
            Write(*,*)errMessage
            Stop
          End IF
            CALL BMD2SETSEGMENTMETADATA(HANDLEBMD2,ISEG,"K"//CHAR(0),KTEXT,IERROR)
          IF(iError .gt. 0)Then
            Call bmd2getlasterror(80,errMessage)
            Write(*,*)errMessage
            Stop
          End IF
                           
          ENDIF
        ENDDO 
      ENDDO 

      CALL BMD2BEGINFRAMES(HANDLEBMD2,IERROR)
      IF(iError .gt. 0)Then
        Call bmd2getlasterror(80,errMessage)
        Write(*,*)errMessage
        Stop
      End IF
!      if(tbegin.ge.tibmd) CALL BMD2_DUMP    !WRITE INITIAL CONDITIONS at tbegin
      
      CLOSE(94)
      
      RETURN
      END
