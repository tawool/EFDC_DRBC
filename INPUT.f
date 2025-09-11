C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE INPUT(TITLE)
C
C**********************************************************************C
C
C **  SUBROUTINE INPUT READS ALL INPUT DATA EXCEPT DATA IN LXLY.INP,
C **  MASK.INP AND RESTART.INP
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 11/07/2001        john hamrick       11/07/2001       john hamrick
C  added body force switch isbodyf (c14) and read of input file
c  fbody.inp containing the body force fbodyfx and fbodyfy
C 11/09/2001        john hamrick       11/09/2001       john hamrick
C  modified bed mechanics coefficient set on (C38)
C  modified bedload functional relationship to more general form (C42A)
C  added bedload outflow/recirculation boundary condtion switch (C42A)
C   and bed load boundary condition file SEDBLBC.INP
C 11/15/2001        john hamrick       11/15/2001       john hamrick
C  added switch IBMECHK on C38 to control form of function FHYDCN
C 11/19/2001        John Hamrick       11/19/2001       John Hamrick
C  replaced LSBDLDBC with LSBLBCU and LSBLBCD and modified format of
C  SEDBLBC.INP
C 11/29/2001        John Hamrick       11/29/2001       John Hamrick
C  fixed the specification of LSBLBCD for no recirculation case
C 12/03/2001        John Hamrick       12/03/2001       John Hamrick
C  initialized NSBDLDBC
C 01/11/2002        John Hamrick       01/11/2002       John Hamrick
C  added ick2cor,ck2uum,ck2vvm,ck2uvm,ck2uuc,ck2vvc,ck2uvc,ck2fcx,
C  ck2fcy,zbrwall,belv1
C 01/31/2002        John Hamrick       01/31/2002       John Hamrick
C  modified toxic-organic carbon sorption options on C44 and added
C  additional input for poc fractions on c45b and c45c
C 02/11/2002        John Hamrick       02/11/2002       John Hamrick
C  moved subgrid scale channel mapping to follow read of modchan.inp
C 03/05/2002        John Hamrick       03/05/2002       John Hamrick
C  modified isdry option to include transport bypass
C 05/22/2002        John Hamrick       05/22/2002       John Hamrick
C  added sed-tox debug flag ISDTXBUG
C
C 7/2023  HUGO RODRIGUEZ    INITIALIZE ARRAYS FOR LOAD TIME SERIES (here since ainit.f is after this call       
C 7/2023  HUGO RODRIGUEZ    ADDED READ OF TIME SERIES OF LOADS      
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION RMULADS(NSTVM),ADDADS(NSTVM)
      DIMENSION QSERSM(NDQSER,KCM)
      character*132 message_out
C
C     DIMENSION CQSE(5),WKQ(KCM)
C
C**********************************************************************C
C
      CHARACTER*80 TEXT,TITLE
      CHARACTER*10 CDUM
	CHARACTER*1 CHARSPACE,CHARSKP1,CHARSKP2
      INTEGER LOADTSFLAG
      REAL    SERLTMP
C
C**********************************************************************C
C
      CHARSKP1='C'
	CHARSKP2='c'
      G=9.81
      PI=3.1415926535898
      PI2=2.*PI
    1 FORMAT (120X)
    2 FORMAT (A80)
C
C**********************************************************************C
C
C----------------------------------------------------------------------C
C
C **  START LOAD TIME SERIES INITIALIZATION    HNR_GHD 7/2023
C       located here since AINIT.F is after this subroutine
C
      DO NC=1,4
        NLSER(NC)=0
        NLIJ(NC)=0
        DO NL=1,NLIJM
          ILDS(NLIJM,NC)=0
          JLDS(NLIJM,NC)=0
          LLDS(NLIJM,NC)=0
          NSERL(NLIJM,NC)=0
        END DO
        DO NL=1,NLSERM
          MLSER(NL,NC)=0
          TCLSER(NL,NC)=0.0
          TALSER(NL,NC)=0.0
	    MLTLAST(NL,NC)=1
          DO K=1,KCM
           SERLT(K,0:NL,NC)=0.0
          END DO
          DO ND=1,NDLSER
            TLSER(ND,NL,NC)=0.0
            DO K=1,KCM
              SERL(ND,K,NL,NC)=0.0
            END DO
          END DO
        END DO
      END DO
      
      DO K=1,KCM
        DO NL=1,NLIJM
          DELTEMP(NL,K)=0.0
        END DO
      END DO
C
C **  END LOAD TIME SERIES INITIALIZATION    HNR_GHD 7/2023
C----------------------------------------------------------------------C
C**********************************************************************C
C
C **  READ MAIN INPUT FILE EFDC.INP
C
      OPEN(1,FILE='EFDC.INP',STATUS='UNKNOWN')
C
C**********************************************************************C
C
C1**  READ TITLE CARD
C
      NCARD=1
      DO NSKIP=1,18
      READ(1,1)
      ENDDO
      READ(1,2) TITLE
      WRITE(7,1002)NCARD
      WRITE(7,2) TITLE
C
C1A**  READ MODE OPTIONS
C
      NCARD=1
      DO NSKIP=1,18
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)IGRIDH,INESTH,IGRIDV,ITIMSOL,ISHOUSATONIC  !IS2TLPG
      WRITE(7,1002)NCARD
      WRITE(7,*) IGRIDH,INESTH,IGRIDV,ITIMSOL,ISHOUSATONIC  !IS2TLPG
      IF(ISO.GT.0) GOTO 100
	IS2TIM=ITIMSOL
C ** TEST, ACTIVATE AND MOVE LOCATIONS FOR IS2TLPG
      IS2TLPG=0
C
C2**  READ RESTART AND DIAGNOSTIC SWITCHES
C
      NCARD=2
      DO NSKIP=1,25
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISRESTI,ISRESTO,ISRESTR,ISPAR,ISLOG,ISDIVEX,
     &          ISNEGH,ISMMC,ISBAL,IDUM,ISHOW
      WRITE(7,1002)NCARD
      WRITE(7,*) ISRESTI,ISRESTO,ISRESTR,ISPAR,ISLOG,ISDIVEX,
     &          ISNEGH,ISMMC,ISBAL,IDUM,ISHOW
      IF(ISO.GT.0) GOTO 100
C
C3**  READ RELAXATION PARAMETERS
C
      NCARD=3
      DO NSKIP=1,24
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) RP,RSQM,ITERM,IRVEC,RPADJ,
     &                 RSQMADJ,ITRMADJ,ITERHPM,IDRYCK,ISDSOLV,FILT3TL
      WRITE(7,1002)NCARD
      WRITE(7,*) RP,RSQM,ITERM,IRVEC,RPADJ,
     &                 RSQMADJ,ITRMADJ,ITERHPM,IDRYCK,ISDSOLV,FILT3TL
      IF(ISO.GT.0) GOTO 100
C
C4**  READ LONGTERM MASS TRANSPORT INTEGRATION ONLY SWITCHES
C
      NCARD=4
      DO NSKIP=1,20
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISLTMT,ISSSMMT,ISLTMTS,ISIA,RPIA,RSQMIA,
     &                     ITRMIA,ISAVEC
      WRITE(7,1002)NCARD
      WRITE(7,*) ISLTMT,ISSSMMT,ISLTMTS,ISIA,RPIA,RSQMIA,
     &                     ITRMIA,ISAVEC
      IF(ISO.GT.0) GOTO 100
C
C5**  READ MOMENTUM ADVECTION AND DIFFUSION SWITCHES AND MISC
C     SWITCHES
C
      NCARD=5
      DO NSKIP=1,32
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISCDMA,ISHDMF,ISDISP,ISWASP,ISDRY,
     *            ISQQ,ISRLID,ISVEG,ISVEGL,ISITB,ISEVER,IINTPG
      WRITE(7,1002)NCARD
      WRITE(7,*) ISCDMA,ISHDMF,ISDISP,ISWASP,ISDRY,
     *            ISQQ,ISRLID,ISVEG,ISVEGL,ISITB,ISEVER,IINTPG
      IF(ISO.GT.0) GOTO 100
      IDRYTBP=0
      IF(ISDRY.LT.0)THEN
	  ISDRY=ABS(ISDRY)
	  IDRYTBP=1
      ENDIF
      IF(ISWASP.EQ.99)ISICM=1
      IF(ISRLID.EQ.1) ISDRY=-1
C      IF(ISDRY.GE.1) IRVEC=2
      IF(ISWASP.EQ.10)ISRCA=1
      JSWAVE=0
      IS1DCHAN=0
      IF(ISCDMA.EQ.10) IS1DCHAN=1
C
      ISCOSMIC=0
      NCARD=6
      DO NSKIP=1,22
      READ(1,1)
      ENDDO
      DO N=0,8
      READ(1,*,IOSTAT=ISO) ISTRAN(N),ISTOPT(N),ISCDCA(N),ISADAC(N),
     &          ISFCT(N),ISPLIT(N),ISADAH(N),ISADAV(N),ISCI(N),ISCO(N)
      IF(ISCDCA(N).GE.4) ISCOSMIC=1
      WRITE(7,1002)NCARD
      WRITE(7,*) ISTRAN(N),ISTOPT(N),ISCDCA(N),ISADAC(N),
     &          ISFCT(N),ISPLIT(N),ISADAH(N),ISADAV(N),ISCI(N),ISCO(N)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C7**  READ TIME-RELATED INTEGER PARAMETERS
C
      NCARD=7
      DO NSKIP=1,19
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) NTC,NTSPTC,NLTC,NTTC,NTCPP,NTSTBC,NTCNB,
     &          NTCVB,NTSMMT,NFLTMT,NDRYSTP,iyear,imon,iday            !hnr
      WRITE(7,1002)NCARD
      WRITE(7,*) NTC,NTSPTC,NLTC,NTTC,NTCPP,NTSTBC,NTCNB,
     &          NTCVB,NTSMMT,NFLTMT,NDRYSTP,iyear,imon,iday              !hnr
      IF(ISO.GT.0) GOTO 100
      if(ntsptc.gt.ntsm) then                                          !hnr
        write(*,109)ntsptc,ntsm                                        !hnr
109     FORMAT('The number of time steps per reference time period'         !hnr
     $'in your application',I8,'  is larger than the maximum allowed of'    !hnr
     $,I8,'  for this compilation')                                         !hnr
        write(*,*)'The program will stop card 7'                                   !hnr
        stop                                                                !hnr
      end if                                                                !hnr
C
C8**  READ TIME-RELATED REAL PARAMETERS
C
      NCARD=8
      DO NSKIP=1,16
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) TCON,TBEGIN,TIDALP,CF,ISCORV,ISDCCA,
     &                     ISCFL,ISCFLM,DTSSFAC,DTSSDHDT
      WRITE(7,1002)NCARD
      WRITE(7,*) TCON,TBEGIN,TIDALP,CF,ISCORV,ISDCCA,
     &                     ISCFL,ISCFLM,DTSSFAC,DTSSDHDT
      IF(ISO.GT.0) GOTO 100
      IF(DTSSFAC.GT.0.0)THEN
        ISDYNSTP=1
        ntsmmt=1          !1 time step output if variable time step  hnr
      ELSE
        ISDYNSTP=0
      ENDIF
      IF(IS2TIM.EQ.0)ISDYNSTP=0
C
C9**  READ HORIZONTAL SPACE RELATED AND SMOOTHING PARAMETERS
C
      NCARD=9
      DO NSKIP=1,27
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) IC,JC,LC,LVC,ISCLO,NDM,LDM,ISMASK,
     &                ISPGNS,NSHMAX,NSBMAX,WSMH,WSMB
      WRITE(7,1002)NCARD
      WRITE(7,*) IC,JC,LC,LVC,ISCLO,NDM,LDM,ISMASK,
     &                ISPGNS,NSHMAX,NSBMAX,WSMH,WSMB
      IF(ISO.GT.0) GOTO 100
      IS2LMC=0
      IF(KC.LT.0) THEN
         KC=-KC
         IS2LMC=1
      ENDIF
C
      NWETCELLSNEW=LVC
	NWETCELLSOLD=LVC
      if(ic.gt.icm) then                                               !hnr
        write(*,110)ic,icm                                             !hnr
110     FORMAT('The number of Columns (I index) in your application'   !hnr
     $  ,I5,'  is larger than the maximum allowed of',I5,'  for this ' !hnr
     $   'compilation')                                                !hnr
        write(*,*)'The program will stop card 9_1'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(jc.gt.jcm) then                                               !hnr
        write(*,120)jc,jcm                                             !hnr
120     FORMAT('The number of Rows (J index) in your application',I5,  !hnr
     $  '  is larger than the maximum allowed of',I5,'  for this '     !hnr
     &   'compilation')                                                !hnr
        write(*,*)'The program will stop card 9_2'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(lc.gt.lcm) then                                               !hnr
        write(*,130)lvc,lcm-2                                          !hnr
130     FORMAT('The number of water cells in your application',I6,     !hnr
     $  '  is larger than the maximum allowed of',I6,                  !hnr
     $  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 9_3'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
C
C9A**  READ VERTICAL SPACE RELATED  PARAMETERS
C
      NCARD=901
      DO NSKIP=1,14
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)KC,KSIG,ISETGVC,SELVREF,BELVREF,ISGVCCK
      WRITE(7,1002)NCARD
      WRITE(7,*) KC,KSIG,ISETGVC,SELVREF,BELVREF,ISGVCCK
      IF(ISO.GT.0) GOTO 100
      if(kc.gt.kcm) then                                               !hnr
        write(*,140)kc,kcm                                             !hnr
140     FORMAT('The number of layers in your application',I3,          !hnr
     $  '  is larger than the maximum allowed of',I3,                  !hnr
     $  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 9A'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr

C
C **  DOMAIN DECOMPOSITION CHECKS FOR HORIZONTAL LOOPS
C
      LCM2T=LC-2
      IF(NDM.EQ.1) LDM=LCM2T
      IF(NDM.GE.2) NCHECK=NDM*LDM
      IF(NDM.GE.2)THEN
        IF(NCHECK.NE.LCM2T)THEN
          WRITE(6,6774)
          STOP
        ENDIF
      ENDIF
C
 6774 FORMAT(' INCONSISTENT DOMAIN DECOMPOSITION NDM,LDW ON CARD 9')
C
C **  ENDDOMAIN DECOMPOSITION CHECKS FOR HORIZONTAL LOOPS
C
C      IF(KC.GE.2.AND.ISVEG.EQ.0) ISITB=0
C     IF(KC.GE.2) ISITB=0
C
C10*  READ LAYER THICKNESS IN VERTICAL
C
      NCARD=10
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      DO K=1,KC
        READ(1,*,IOSTAT=ISO)KDUM,DZC(K)
        WRITE(7,1002)NCARD
        WRITE(7,*)KDUM,DZC(K)
        IF(ISO.GT.0) GOTO 100
      ENDDO
C
C11*  READ GRID, ROUGHNESS, MASKING AND DEPTH PARAMETERS
C
      NCARD=11
      DO NSKIP=1,18
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) DX,DY,DXYCVT,IMDXDY,ZBRADJ,ZBRCVRT,HMIN,
     &                  HADJ,HCVRT,HDRY,HWET,BELADJ,BELCVRT
      WRITE(7,1002)NCARD
      WRITE(7,*) DX,DY,DXYCVT,IMDXDY,ZBRADJ,ZBRCVRT,HMIN,
     &                  HADJ,HCVRT,HDRY,HWET,BELADJ,BELCVRT
      IF(ISO.GT.0) GOTO 100
C     READ(1,*,IOSTAT=ISO) DX,DY,DXYCVT,ZBRADJ,ZBRCVRT,HMIN,HADJ,HCVRT,
C    &                  HDRY,HWET,BELADJ,BELCVRT,RMWET
C
C11A* READ TWO-LAYER MOMENTUM FLUX AND CURVATURE ACCELERATION
C     CORRECTION FACTORS
C
      NCARD=11
      DO NSKIP=1,16
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ICK2COR,CK2UUM,CK2VVM,CK2UVM,CK2UUC,
     &                  CK2VVC,CK2UVC,CK2FCX,CK2FCY
      WRITE(7,1002)NCARD
      WRITE(7,*) ICK2COR,CK2UUM,CK2VVM,CK2UVM,CK2UUC,
     &                  CK2VVC,CK2UVC,CK2FCX,CK2FCY
      IF(ISO.GT.0) GOTO 100
      IF(ICK2COR.GE.1) THEN
        IS2LMC=ICK2COR
      END IF
C
C11B* READ CORNER CELL BOTTOM STRESS CORRECTION OPTIONS
C
      NCARD=11
      DO NSKIP=1,10
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)ISCORTBC,ISCORTBCD,FSCORTBC
      WRITE(7,1002)NCARD
      WRITE(7,*) ISCORTBC,ISCORTBCD,FSCORTBC
      IF(ISO.GT.0) GOTO 100
C
C12*  READ TURBULENT DIFFUSION PARAMETERS
C
      NCARD=12
      DO NSKIP=1,16
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) AHO,AHD,AVO,ABO,AVMN,ABMN,VISMUD,AVCON,
     &  ZBRWALL
        WRITE(7,1002)NCARD
        WRITE(7,*) AHO,AHD,AVO,ABO,AVMN,ABMN,VISMUD,AVCON,ZBRWALL
      IF(ISO.GT.0) GOTO 100
C
C12A*  READ TURBULENCE CLOSURE OPTIONS
C
      NCARD=121
      DO NSKIP=1,25
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)ISTOPT(0),ISSQL,ISAVBMN,ISFAVB,ISINWV,
     &        ISLLIM,IFPROX,ISVTURB,VTURBEFF
        WRITE(7,1002)NCARD
        WRITE(7,*)ISTOPT(0),ISSQL,ISAVBMN,ISFAVB,ISINWV,
     &        IFPROX,ISVTURB,VTURBEFF
      IF(ISO.GT.0) GOTO 100
C
C13*  READ TURBULENCE CLOSURE PARAMETERS
C
      NCARD=13
      DO NSKIP=1,17
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) VKC,CTURB,CTURB2B,CTE1,CTE2,CTE3,CTE4,
     &                      CTE5,RIQMAX,QQMIN,QQLMIN,DMLMIN
      WRITE(7,1002)NCARD
      WRITE(7,*) VKC,CTURB,CTURB2B,CTE1,CTE2,CTE3,CTE4,
     &                      CTE5,RIQMAX,QQMIN,QQLMIN,DMLMIN
      IF(ISO.GT.0) GOTO 100
C
C14*  READ TIDAL & ATMOSPHERIC FORCING, GROUND WATER
C     AND SUBGRID CHANNEL PARAMETERS
C
      NCARD=14
      DO NSKIP=1,18
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) MTIDE,NWSER,NASER,ISGWIT,ISCHAN,ISWAVE,
     &     ITIDASM,ISPERC,ISBODYF,ISPNHYDS
      WRITE(7,1002)NCARD
      WRITE(7,*) MTIDE,NWSER,NASER,ISGWIT,ISCHAN,ISWAVE,ITIDASM,
     &     ISPERC,ISBODYF,ISPNHYDS
      ISWCBL=0
      ISWVSD=0
C      IF(ISWAVE.GT.0) ISWCBL=1
      IF(ISO.GT.0) GOTO 100
      if(mtide.gt.npform) then                                         !hnr
        write(*,150)mtide,npform                                       !hnr
150     FORMAT('The number of periodic forcing functions (MTIDE) '     !hnr
     $  'in your application',I6,                                      !hnr
     $  '  is larger than the maximum allowed of',I6,                  !hnr
     $  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card14_1'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(nwser.gt.nwserm) then                                         !hnr
        write(*,160)nwser,nwserm                                       !hnr
160     FORMAT('The number of Wind Time Series in your application'    !hnr
     $  ,I3,'  is larger than the maximum allowed of',I3,              !hnr
     $  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 14_2'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(naser.gt.naserm) then                                         !hnr
        write(*,170)ic,icm                                             !hnr
170     FORMAT('The number of Atmospheric Condition Time Series '      !hnr
     $  'in your application',I3,                                      !hnr
     $  '  is larger than the maximum allowed of',I3,                  !hnr
     $  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 14_3'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
C
C15*  READ PERIODIC FORCING (TIDAL) CONSTITUENT SYMBOLS AND PERIODS
C
      NCARD=15
      DO NSKIP=1,7
      READ(1,1)
      ENDDO
      DO M=1,MTIDE
      READ(1,*,IOSTAT=ISO) SYMBOL(M),TCP(M)
      WRITE(7,1002)NCARD
      WRITE(7,*) SYMBOL(M),TCP(M)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C16*  READ SURFACE ELEVATION OR PRESSURE BOUNDARY CONDITION PARAMETERS
C
      NCARD=16
      DO NSKIP=1,17
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) NPBS,NPBW,NPBE,NPBN,NPFOR,NPFORT,
     &                     NPSER,PDGINIT
      WRITE(7,1002)NCARD
      WRITE(7,*) NPBS,NPBW,NPBE,NPBN,NPFOR,NPFORT,NPSER,PDGINIT
      IF(ISO.GT.0) GOTO 100
      if(npbs.gt.npbsm) then                                           !hnr
        write(*,180)npbs,npbsm                                         !hnr
180     FORMAT('The number of WSE BC South Boundary '                  !hnr
     $  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 16_1'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(npbw.gt.npbwm) then                                           !hnr
        write(*,190)npbw,npbwm                                         !hnr
190     FORMAT('The number of WSE BC West Boundary '                   !hnr
     $  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     $  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 16_2'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(npbn.gt.npbnm) then                                           !hnr
        write(*,200)npbn,npbnm                                         !hnr
200     FORMAT('The number of WSE BC North Boundary'                   !hnr
     $  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 16_3'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(npbe.gt.npbem) then                                           !hnr
        write(*,210)npbe,npbem                                         !hnr
210     FORMAT('The number of WSE BC East Boundary '                   !hnr
     &  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     &   '  for this compilation')                                     !hnr
        write(*,*)'The program will stop card 16_4'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(npser.gt.npserm) then                                         !hnr
        write(*,220)npser,npserm                                       !hnr
220     FORMAT('The number of PSER Time Series'                        !hnr
     &  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     $  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 16_5'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
C
C17*  READ PERIODIC FORCING (TIDAL) SURFACE ELEVATION OR
C     PRESSURE BOUNDARY CONDITION FORCINGS
C
      IF(NPFORT.GE.1)THEN
	  OPEN(2,FILE='TIDALBC.OUT')
	  CLOSE(2,STATUS='DELETE')
	  OPEN(2,FILE='TIDALBC.OUT')
      ENDIF
C
      NCARD=17
      DO NSKIP=1,15
      READ(1,1)
      ENDDO
      DO NP=1,NPFOR
      DO M=1,MTIDE
	IF(NPFORT.EQ.0)THEN
        READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M)
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M)
        IF(ISO.GT.0) GOTO 100
      ENDIF
	IF(NPFORT.GE.1)THEN
        READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),PFX0(NP,M)
        RAD=PI2*PFPH(NP,M)/TCP(M)
        CPFAM0(NP,M)=PFAM(NP,M)*COS(RAD)
        SPFAM0(NP,M)=PFAM(NP,M)*SIN(RAD)
	  WRITE(2,2068)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),
     &                         CPFAM0(NP,M),SPFAM0(NP,M)
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),
     &                      CPFAM0(NP,M),SPFAM0(NP,M)
        IF(ISO.GT.0) GOTO 100
      ENDIF
	IF(NPFORT.GE.1)THEN
        READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),PFX1(NP,M)
        RAD=PI2*PFPH(NP,M)/TCP(M)
        CPFAM1(NP,M)=PFAM(NP,M)*COS(RAD)-CPFAM0(NP,M)
        SPFAM1(NP,M)=PFAM(NP,M)*SIN(RAD)-SPFAM0(NP,M)
	  WRITE(2,2068)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),
     &                         CPFAM1(NP,M),SPFAM1(NP,M)
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),
     &                         CPFAM1(NP,M),SPFAM1(NP,M)
        CPFAM2(NP,M)=0.0
        SPFAM2(NP,M)=0.0
      ENDIF
	IF(NPFORT.EQ.2)THEN
        READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),PFX2(NP,M)
        RAD=PI2*PFPH(NP,M)/TCP(M)
	  IF(PFX2(NP,M).GT.0.0)THEN
          CPFAM2(NP,M)=PFAM(NP,M)*COS(RAD)-CPFAM0(NP,M)
          SPFAM2(NP,M)=PFAM(NP,M)*SIN(RAD)-SPFAM0(NP,M)
	  ELSE
          CPFAM2(NP,M)=0.
          SPFAM2(NP,M)=0.
        ENDIF
	  WRITE(2,2068)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),
     &                         CPFAM2(NP,M),SPFAM2(NP,M)
      ENDIF
      ENDDO
      ENDDO
c
C18*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
C     ON SOUTH OPEN BOUNDARIES
C
      NCARD=18
      DO NSKIP=1,13
      READ(1,1)
      ENDDO
	IF(NPFORT.EQ.0)THEN
        DO L=1,NPBS
          READ(1,*,IOSTAT=ISO) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
            RAD=PI2*PFPH(NPFORS,M)/TCP(M)
            AMP=G*PFAM(NPFORS,M)
            PCBS(L,M)=AMP*COS(RAD)
            PSBS(L,M)=AMP*SIN(RAD)
          ENDDO
        ENDDO
      ENDIF
	IF(NPFORT.EQ.1)THEN
        DO L=1,NPBS
          READ(1,*,IOSTAT=ISO) IPBS(L),JPBS(L),ISPBS(L),NPFORS,
     &                         NPSERS(L),NPSERS1(L),TPCOORDS(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L),
     &               NPSERS1(L),TPCOORDS(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
            PCBS(L,M)=CPFAM0(NPFORS,M)+TPCOORDS(L)*CPFAM1(NPFORS,M)
     &               +TPCOORDS(L)*TPCOORDS(L)*CPFAM2(NPFORS,M)
            PSBS(L,M)=SPFAM0(NPFORS,M)+TPCOORDS(L)*SPFAM1(NPFORS,M)
     &               +TPCOORDS(L)*TPCOORDS(L)*SPFAM2(NPFORS,M)
            TMPAMP=SQRT(PCBS(L,M)*PCBS(L,M)+PSBS(L,M)*PSBS(L,M))
		  TMPPHS=ATAN2(PSBS(L,M),PCBS(L,M))
		  TMPPHS=TMPPHS*TCP(M)/PI2
		  IF(TMPPHS.LT.0.0)TMPPHS=TMPPHS+TCP(M)
            WRITE(2,2069)TPCOORDS(L),TMPAMP,TMPPHS,PCBS(L,M),PSBS(L,M),
     &                   IPBS(L),JPBS(L),SYMBOL(M)
            PCBS(L,M)=G*PCBS(L,M)
            PSBS(L,M)=G*PSBS(L,M)
          ENDDO
        ENDDO
      ENDIF
	IF(NPFORT.EQ.2)THEN
        DO L=1,NPBS
          READ(1,*,IOSTAT=ISO) IPBS(L),JPBS(L),ISPBS(L),NPFORS,
     &                         NPSERS(L),NPSERS1(L),TPCOORDS(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L),
     &               NPSERS1(L),TPCOORDS(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
	      BOTTOM=PFX2(NPFORS,M)*(1.0-PFX2(NPFORS,M))
	      TOP1=TPCOORDS(L)*PFX2(NPFORS,M)*(TPCOORDS(L)-PFX2(NPFORS,M))
	      TOP2=TPCOORDS(L)*(1.0-TPCOORDS(L))
	      IF(BOTTOM.EQ.0.0)THEN
	        TOP1=TPCOORDS(L)
	        TOP2=TPCOORDS(L)*TPCOORDS(L)
            ELSE
	        TOP1=TOP1/BOTTOM
	        TOP2=TOP2/BOTTOM
            ENDIF
            PCBS(L,M)=CPFAM0(NPFORS,M)+TOP1*CPFAM1(NPFORS,M)
     &               +TOP2*CPFAM2(NPFORS,M)
            PSBS(L,M)=SPFAM0(NPFORS,M)+TOP1*SPFAM1(NPFORS,M)
     &               +TOP2*SPFAM2(NPFORS,M)
            TMPAMP=SQRT(PCBS(L,M)*PCBS(L,M)+PSBS(L,M)*PSBS(L,M))
		  TMPPHS=ATAN2(PSBS(L,M),PCBS(L,M))
		  TMPPHS=TMPPHS*TCP(M)/PI2
		  IF(TMPPHS.LT.0.0)TMPPHS=TMPPHS+TCP(M)
            WRITE(2,2069)TPCOORDS(L),TMPAMP,TMPPHS,PCBS(L,M),PSBS(L,M),
     &                   IPBS(L),JPBS(L),SYMBOL(M)
            PCBS(L,M)=G*PCBS(L,M)
            PSBS(L,M)=G*PSBS(L,M)
          ENDDO
        ENDDO
      ENDIF
C
 2068 FORMAT(I4,3X,A2,5X,E14.4,3E14.5,5X,2I5)
 2069 FORMAT(5E14.5,2I8,3X,A2)
C
C19*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
C     ON WEST OPEN BOUNDARIES
C
      NCARD=19
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
	IF(NPFORT.EQ.0)THEN
        DO L=1,NPBW
          READ(1,*,IOSTAT=ISO) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
            RAD=PI2*PFPH(NPFORW,M)/TCP(M)
            AMP=G*PFAM(NPFORW,M)
            PCBW(L,M)=AMP*COS(RAD)
            PSBW(L,M)=AMP*SIN(RAD)
          ENDDO
        ENDDO
      ENDIF
	IF(NPFORT.EQ.1)THEN
        DO L=1,NPBW
          READ(1,*,IOSTAT=ISO) IPBW(L),JPBW(L),ISPBW(L),NPFORW,
     &                         NPSERW(L),NPSERW1(L),TPCOORDW(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L),
     &               NPSERW1(L),TPCOORDW(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
            PCBW(L,M)=CPFAM0(NPFORW,M)+TPCOORDW(L)*CPFAM1(NPFORW,M)
     &               +TPCOORDW(L)*TPCOORDW(L)*CPFAM2(NPFORW,M)
            PSBW(L,M)=SPFAM0(NPFORW,M)+TPCOORDW(L)*SPFAM1(NPFORW,M)
     &               +TPCOORDW(L)*TPCOORDW(L)*SPFAM2(NPFORW,M)
            TMPAMP=SQRT(PCBW(L,M)*PCBW(L,M)+PSBW(L,M)*PSBW(L,M))
            TMPPHS=0.0
            IF(TMPAMP.GT.0.0) TMPPHS=ATAN2(PSBW(L,M),PCBW(L,M))
		  TMPPHS=TMPPHS*TCP(M)/PI2
		  IF(TMPPHS.LT.0.0)TMPPHS=TMPPHS+TCP(M)
            WRITE(2,2069)TPCOORDW(L),TMPAMP,TMPPHS,PCBW(L,M),PSBW(L,M),
     &                   IPBW(L),JPBW(L),SYMBOL(M)
            PCBW(L,M)=G*PCBW(L,M)
            PSBW(L,M)=G*PSBW(L,M)
          ENDDO
        ENDDO
      ENDIF
	IF(NPFORT.EQ.2)THEN
        DO L=1,NPBW
          READ(1,*,IOSTAT=ISO) IPBW(L),JPBW(L),ISPBW(L),NPFORW,
     &                         NPSERW(L),NPSERW1(L),TPCOORDW(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L),
     &               NPSERW1(L),TPCOORDW(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
	      BOTTOM=PFX2(NPFORW,M)*(1.0-PFX2(NPFORW,M))
	      TOP1=TPCOORDW(L)*PFX2(NPFORW,M)*(TPCOORDW(L)-PFX2(NPFORW,M))
	      TOP2=TPCOORDW(L)*(1.0-TPCOORDW(L))
	      IF(BOTTOM.EQ.0.0)THEN
	        TOP1=TPCOORDW(L)
	        TOP2=TPCOORDW(L)*TPCOORDW(L)
            ELSE
	        TOP1=TOP1/BOTTOM
	        TOP2=TOP2/BOTTOM
            ENDIF
            PCBW(L,M)=CPFAM0(NPFORW,M)+TOP1*CPFAM1(NPFORW,M)
     &               +TOP2*CPFAM2(NPFORW,M)
            PSBW(L,M)=SPFAM0(NPFORW,M)+TOP1*SPFAM1(NPFORW,M)
     &               +TOP2*SPFAM2(NPFORW,M)
            TMPAMP=SQRT(PCBW(L,M)*PCBW(L,M)+PSBW(L,M)*PSBW(L,M))
            TMPPHS=0.0
            IF(TMPAMP.GT.0.0) TMPPHS=ATAN2(PSBW(L,M),PCBW(L,M))
		  TMPPHS=TMPPHS*TCP(M)/PI2
		  IF(TMPPHS.LT.0.0)TMPPHS=TMPPHS+TCP(M)
            WRITE(2,2069)TPCOORDW(L),TMPAMP,TMPPHS,PCBW(L,M),PSBW(L,M),
     &                   IPBW(L),JPBW(L),SYMBOL(M)
            PCBW(L,M)=G*PCBW(L,M)
            PSBW(L,M)=G*PSBW(L,M)
          ENDDO
        ENDDO
      ENDIF
C
C20*  READ PERIODIC FORCING (TIDAL)ELEVATION BOUNDARY CONDTIONS
C     ON EAST OPEN BOUNDARIES
C
      NCARD=20
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
	IF(NPFORT.EQ.0)THEN
        DO L=1,NPBE
          READ(1,*,IOSTAT=ISO) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
            RAD=PI2*PFPH(NPFORE,M)/TCP(M)
            AMP=G*PFAM(NPFORE,M)
            PCBE(L,M)=AMP*COS(RAD)
            PSBE(L,M)=AMP*SIN(RAD)
          ENDDO
        ENDDO
      ENDIF
	IF(NPFORT.EQ.1)THEN
        DO L=1,NPBE
          READ(1,*,IOSTAT=ISO) IPBE(L),JPBE(L),ISPBE(L),NPFORE,
     &                         NPSERE(L),NPSERE1(L),TPCOORDE(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L),
     &               NPSERE1(L),TPCOORDE(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
            PCBE(L,M)=CPFAM0(NPFORE,M)+TPCOORDE(L)*CPFAM1(NPFORE,M)
     &               +TPCOORDE(L)*TPCOORDE(L)*CPFAM2(NPFORE,M)
            PSBE(L,M)=SPFAM0(NPFORE,M)+TPCOORDE(L)*SPFAM1(NPFORE,M)
     &               +TPCOORDE(L)*TPCOORDE(L)*SPFAM2(NPFORE,M)
            TMPAMP=SQRT(PCBE(L,M)*PCBE(L,M)+PSBE(L,M)*PSBE(L,M))
            TMPPHS=0.0
            IF(TMPAMP.GT.0.0) TMPPHS=ATAN2(PSBE(L,M),PCBE(L,M))
		  TMPPHS=TMPPHS*TCP(M)/PI2
		  IF(TMPPHS.LT.0.0)TMPPHS=TMPPHS+TCP(M)
            WRITE(2,2069)TPCOORDE(L),TMPAMP,TMPPHS,PCBE(L,M),PSBE(L,M),
     &                   IPBE(L),JPBE(L),SYMBOL(M)
            PCBE(L,M)=G*PCBE(L,M)
            PSBE(L,M)=G*PSBE(L,M)
          ENDDO
        ENDDO
      ENDIF
	IF(NPFORT.EQ.2)THEN
        DO L=1,NPBE
          READ(1,*,IOSTAT=ISO) IPBE(L),JPBE(L),ISPBE(L),NPFORE,
     &                         NPSERE(L),NPSERE1(L),TPCOORDE(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L),
     &               NPSERE1(L),TPCOORDE(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
	      BOTTOM=PFX2(NPFORE,M)*(1.0-PFX2(NPFORE,M))
	      TOP1=TPCOORDE(L)*PFX2(NPFORE,M)*(TPCOORDE(L)-PFX2(NPFORE,M))
	      TOP2=TPCOORDE(L)*(1.0-TPCOORDE(L))
	      IF(BOTTOM.EQ.0.0)THEN
	        TOP1=TPCOORDE(L)
	        TOP2=TPCOORDE(L)*TPCOORDE(L)
            ELSE
	        TOP1=TOP1/BOTTOM
	        TOP2=TOP2/BOTTOM
            ENDIF
            PCBE(L,M)=CPFAM0(NPFORE,M)+TOP1*CPFAM1(NPFORE,M)
     &               +TOP2*CPFAM2(NPFORE,M)
            PSBE(L,M)=SPFAM0(NPFORE,M)+TOP1*SPFAM1(NPFORE,M)
     &               +TOP2*SPFAM2(NPFORE,M)
            TMPAMP=SQRT(PCBE(L,M)*PCBE(L,M)+PSBE(L,M)*PSBE(L,M))
            TMPPHS=0.0
            IF(TMPAMP.GT.0.0) TMPPHS=ATAN2(PSBE(L,M),PCBE(L,M))
		  TMPPHS=TMPPHS*TCP(M)/PI2
		  IF(TMPPHS.LT.0.0)TMPPHS=TMPPHS+TCP(M)
            WRITE(2,2069)TPCOORDE(L),TMPAMP,TMPPHS,PCBE(L,M),PSBE(L,M),
     &                   IPBE(L),JPBE(L),SYMBOL(M)
            PCBE(L,M)=G*PCBE(L,M)
            PSBE(L,M)=G*PSBE(L,M)
          ENDDO
        ENDDO
      ENDIF
C
C21*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
C     ON NORTH OPEN BOUNDARIES
C
      NCARD=21
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
	IF(NPFORT.EQ.0)THEN
        DO L=1,NPBN
          READ(1,*,IOSTAT=ISO) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
            RAD=PI2*PFPH(NPFORN,M)/TCP(M)
            AMP=G*PFAM(NPFORN,M)
            PCBN(L,M)=AMP*COS(RAD)
            PSBN(L,M)=AMP*SIN(RAD)
          ENDDO
        ENDDO
      ENDIF
	IF(NPFORT.GE.1)THEN
        DO L=1,NPBN
          READ(1,*,IOSTAT=ISO) IPBN(L),JPBN(L),ISPBN(L),NPFORN,
     &                         NPSERN(L),NPSERN1(L),TPCOORDN(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L),
     &               NPSERN1(L),TPCOORDN(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
            PCBN(L,M)=CPFAM0(NPFORN,M)+TPCOORDN(L)*CPFAM1(NPFORN,M)
     &               +TPCOORDN(L)*TPCOORDN(L)*CPFAM2(NPFORN,M)
            PSBN(L,M)=SPFAM0(NPFORN,M)+TPCOORDN(L)*SPFAM1(NPFORN,M)
     &               +TPCOORDN(L)*TPCOORDN(L)*SPFAM2(NPFORN,M)
            TMPAMP=SQRT(PCBN(L,M)*PCBN(L,M)+PSBN(L,M)*PSBN(L,M))
            TMPPHS=0.0
            IF(TMPAMP.GT.0.0) TMPPHS=ATAN2(PSBN(L,M),PCBN(L,M))
		  TMPPHS=TMPPHS*TCP(M)/PI2
		  IF(TMPPHS.LT.0.0)TMPPHS=TMPPHS+TCP(M)
            WRITE(2,2069)TPCOORDN(L),TMPAMP,TMPPHS,PCBN(L,M),PSBN(L,M),
     &                   IPBN(L),JPBN(L),SYMBOL(M)
            PCBN(L,M)=G*PCBN(L,M)
            PSBN(L,M)=G*PSBN(L,M)
         ENDDO
        ENDDO
      ENDIF
	IF(NPFORT.EQ.2)THEN
        DO L=1,NPBN
          READ(1,*,IOSTAT=ISO) IPBN(L),JPBN(L),ISPBN(L),NPFORN,
     &                         NPSERN(L),NPSERN1(L),TPCOORDN(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L),
     &               NPSERN1(L),TPCOORDN(L)
          IF(ISO.GT.0) GOTO 100
          DO M=1,MTIDE
	      BOTTOM=PFX2(NPFORN,M)*(1.0-PFX2(NPFORN,M))
	      TOP1=TPCOORDN(L)*PFX2(NPFORN,M)*(TPCOORDN(L)-PFX2(NPFORN,M))
	      TOP2=TPCOORDN(L)*(1.0-TPCOORDN(L))
	      IF(BOTTOM.EQ.0.0)THEN
	        TOP1=TPCOORDN(L)
	        TOP2=TPCOORDN(L)*TPCOORDN(L)
            ELSE
	        TOP1=TOP1/BOTTOM
	        TOP2=TOP2/BOTTOM
            ENDIF
            PCBN(L,M)=CPFAM0(NPFORN,M)+TOP1*CPFAM1(NPFORN,M)
     &               +TOP2*CPFAM2(NPFORN,M)
            PSBN(L,M)=SPFAM0(NPFORN,M)+TOP1*SPFAM1(NPFORN,M)
     &               +TOP2*SPFAM2(NPFORN,M)
            TMPAMP=SQRT(PCBN(L,M)*PCBN(L,M)+PSBN(L,M)*PSBN(L,M))
            TMPPHS=0.0
            IF(TMPAMP.GT.0.0) TMPPHS=ATAN2(PSBN(L,M),PCBN(L,M))
		  TMPPHS=TMPPHS*TCP(M)/PI2
		  IF(TMPPHS.LT.0.0)TMPPHS=TMPPHS+TCP(M)
            WRITE(2,2069)TPCOORDN(L),TMPAMP,TMPPHS,PCBN(L,M),PSBN(L,M),
     &                   IPBN(L),JPBN(L),SYMBOL(M)
            PCBN(L,M)=G*PCBN(L,M)
            PSBN(L,M)=G*PSBN(L,M)
          ENDDO
        ENDDO
      ENDIF
C
      IF(NPFORT.GE.1)THEN
	  CLOSE(2)
      ENDIF
C
C21A* READ WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
C
      NCARD=211
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
C
      READ(1,*,IOSTAT=ISO)ISWSEDA,NLWSEDA,ISUVDA,NLUVDA,NUVSER
      WRITE(7,1002)NCARD
      WRITE(7,*) ISWSEDA,NLWSEDA,ISUVDA,NLUVDA,NUVSER
      IF(ISO.GT.0) GOTO 100
C
      DO N=1,NUVSER
	  MUVTLAST(N)=1
	ENDDO
C
C21B* WATER SURFACE ELEVATION DATA ASSIMILATION
C
      NCARD=212
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
C
      IF(ISWSEDA.GT.0.0)THEN
	  DO N=1,NLWSEDA
          READ(1,*,IOSTAT=ISO)ICWSEDA(N),JCWSEDA(N),NWSESERA(N),
     &                        TSWSEDA(N)
        ENDDO
        WRITE(7,1002)NCARD
	  DO N=1,NLWSEDA
          WRITE(7,*)ICWSEDA(N),JCWSEDA(N),NWSESERA(N),
     &                        TSWSEDA(N)
        ENDDO
        IF(ISO.GT.0) GOTO 100
	ENDIF
C
C21C* VELOCITY DATA ASSIMILATION
C
      NCARD=213
      DO NSKIP=1,13
      READ(1,1)
      ENDDO
C
      IF(ISUVDA.GT.0.0)THEN
	  IF(NLUVDA.GT.0.AND.NUVSER.GT.0)THEN
	    DO N=1,NLUVDA
            READ(1,*,IOSTAT=ISO)ICUVDA(N),JCUVDA(N),NUVSERA(N),
     &        TSUVDA(N),FSUVDA(N),IWUVDA(N),IRVUDA(N),RRUVDA(N)
          ENDDO
          WRITE(7,1002)NCARD
	    DO N=1,NLUVDA
            WRITE(7,*)ICUVDA(N),JCUVDA(N),NUVSERA(N),
     &        TSUVDA(N),FSUVDA(N),IWUVDA(N),IRVUDA(N),RRUVDA(N)
          ENDDO
          IF(ISO.GT.0) GOTO 100
	  ENDIF
	  IF(NLUVDA.GT.0.AND.NUVSER.EQ.0)THEN
	    N=1
          READ(1,*,IOSTAT=ISO)IDUM,JDUM,NUVSERA(N),
     &        TSUVDA(N),FSUVDA(N),IWUVDA(N),IRVUDA(N),RRUVDA(N)
          WRITE(7,1002)NCARD
          WRITE(7,*)IDUM,JDUM,NUVSERA(N),
     &        TSUVDA(N),FSUVDA(N),IWUVDA(N),IRVUDA(N),RRUVDA(N)
          IF(ISO.GT.0) GOTO 100
	  ENDIF
	ENDIF
C
C22*  READ NUM OF SEDIMENT AMD TOXICS AND NUM OF CONCENTRATION TIME SERIES
C
      NCARD=22
      DO NSKIP=1,19
      READ(1,1)
      ENDDO
C                                          SAL      TEM      DYE
      READ(1,*,IOSTAT=ISO) NTOX,NSED,NSND,NCSER(1),NCSER(2),NCSER(3),
     &           NCSER(4),NTOXSER,NSEDSER,NSNDSER,ISSBAL
      WRITE(7,1002)NCARD
      WRITE(7,*) NTOX,NSED,NSND,NCSER(1),NCSER(2),NCSER(3),
     &           NCSER(4),NTOXSER,NSEDSER,NSNDSER,ISSBAL
      IF(ISO.GT.0) GOTO 100
      MTMP=4
      DO N=1,NTOX
       MTMP=MTMP+1
       MSVTOX(N)=MTMP
      ENDDO
      DO N=1,NSED
       MTMP=MTMP+1
       MSVSED(N)=MTMP
      ENDDO
      DO N=1,NSND
       MTMP=MTMP+1
       MSVSND(N)=MTMP
      ENDDO
      DO N=1,NTOX
       M=MSVTOX(N)
       NCSER(M)=NTOXSER
      ENDDO
      DO N=1,NSED
       M=MSVSED(N)
       NCSER(M)=NSEDSER
      ENDDO
      DO N=1,NSND
       M=MSVSND(N)
       NCSER(M)=NSNDSER
      ENDDO
      if(ncser(1).gt.ncserm) then                                      !hnr
        write(*,230)ncser(1),ncserm                                    !hnr
230     FORMAT('The number of Salinity Time Series '                   !hnr
     &  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 22_1'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(ncser(2).gt.ncserm) then                                      !hnr
        write(*,240)ncser(2),ncserm                                    !hnr
240     FORMAT('The number of Temperature Time Series '                !hnr
     &  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     *  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 22_2'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(ncser(3).gt.ncserm) then                                      !hnr
        write(*,250)ncser(3),ncserm                                    !hnr
250     FORMAT('The number of Dye Time Series in your application'     !hnr
     &  ,I5,'  is larger than the maximum allowed of',I5,              !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 22_3'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
C
C23*  READ VELOCITY, VOL SOUR/SINK, FLOW CONTROL, & WITHDRAW/RETURN DATA
C
      NCARD=23
      DO NSKIP=1,17
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) NQSIJ,NQJPIJ,NQSER,NQCTL,
     &                  NQCTLT,NQWR,NQWRSR,ISDIQ
      WRITE(7,1002)NCARD
      WRITE(7,*) NVBS,NUBW,NUBE,NVBN,NQSIJ,NQJPIJ,NQSER,NQCTL,
     &                  NQCTLT,NQWR,NQWRSR,ISDIQ
      IF(ISO.GT.0) GOTO 100
      if(nqsij.gt.nqsijm) then                                         !hnr
        write(*,260)nqsij,nqsijm                                       !hnr
260     FORMAT('The number of Inflow Locations in your application'    !hnr
     &  ,I5,'  is larger than the maximum allowed of',I5,              !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 23_1'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(nqser.gt.nqserm) then                                         !hnr
        write(*,270)nqser,nqserm                                       !hnr
270     FORMAT('The number of Inflow Time Series in your application'  !hnr
     &  ,I5,'  is larger than the maximum allowed of',I5,              !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 23_2'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(nqctl.gt.nqctlm) then                                         !hnr
        write(*,280)nqctl,nqctlm                                       !hnr
280     FORMAT('The number of Pressure Controlled Pairs '              !hnr
     &  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     *  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 23_3'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(nqctlt.gt.nqcttm) then                                        !hnr
        write(*,290)nqctlt,nqcttm                                      !hnr
290     FORMAT('The number of Pressure Controlled Tables '             !hnr
     &  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     *  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 23_4'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(nqwr.gt.nqwrm) then                                           !hnr
        write(*,300)nqwr,nqwrm                                         !hnr
300     FORMAT('The number of Specified Withdrawal/Return Pairs '      !hnr
     &  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     *  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 23_5'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(nqwrsr.gt.nqwrsrm) then                                       !hnr
        write(*,310)nqwrsr,nqwrsrm                                     !hnr
310     FORMAT('The number of Time Series Specifying '                 !hnr
     &  'Withdrawal/Return Pairs in your application',I5,              !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     *  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 23_6'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
C
C24*  READ VOLUMN SOURCE/SINK LOCATIONS, MAGNITUDES, & VOL & CONC SERIES
C
      NCARD=24
      DO NSKIP=1,27
      READ(1,1)
      ENDDO
      DO L=1,NQSIJ
       READ(1,*,IOSTAT=ISO)IQS(L),JQS(L),QSSE,NQSMUL(L),NQSMF(L),
     &          NQSERQ(L),NCSERQ(L,1),NCSERQ(L,2),NCSERQ(L,3),
     &          NCSERQ(L,4),NTOXSRQ,NSEDSRQ,NSNDSRQ,QFACTOR(L)
       WRITE(7,1002)NCARD
       WRITE(7,*)IQS(L),JQS(L),QSSE,NQSMUL(L),NQSMF(L),
     &          NQSERQ(L),NCSERQ(L,1),NCSERQ(L,2),NCSERQ(L,3),
     &          NCSERQ(L,4),NTOXSRQ,NSEDSRQ,NSNDSRQ,QFACTOR(L)
       IF(ISO.GT.0) GOTO 100
       DO K=1,KC
       QSS(K,L)=QSSE*DZC(K)
       ENDDO
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERQ(L,M)=NTOXSRQ
       ENDDO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERQ(L,M)=NSEDSRQ
       ENDDO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERQ(L,M)=NSNDSRQ
       ENDDO
      ENDDO
C
C25*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRATIONS
C     SAL,TEM,DYE,SFL,TOX(1 TO NOTX)
C
      NCARD=25
      DO NSKIP=1,12
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NQSIJ
      READ(1,*,IOSTAT=ISO) (CQSE(M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQSE(M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
       DO MS=1,MMAX
        DO K=1,KC
        CQS(K,L,MS)=CQSE(MS)
        ENDDO
       ENDDO
      ENDDO
C
C26*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRATIONS
C     SED(1 TO NSED),SND(1 TO NSND)
C
      NCARD=26
      DO NSKIP=1,13
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NQSIJ
      READ(1,*,IOSTAT=ISO) (CQSE(M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQSE(M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
       DO MS=MMIN,MMAX
        DO K=1,KC
        CQS(K,L,MS)=CQSE(MS)
        ENDDO
       ENDDO
      ENDDO
C
C27*  READ JET/PLUME SOURCE LOCATIONS AND PARAMETERS
C
      NCARD=27
      DO NSKIP=1,19
      READ(1,1)
      ENDDO
      DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) IDUM,ICALJP(L),IQJP(L),JQJP(L),KQJP(L),
     &   NPORTJP(L),XJETL(L),YJETL(L),ZJET(L),PHJET(L),THJET(L),
     &                      DJET(L),CFRD(L),DJPER(L)
      WRITE(7,1002)NCARD
      WRITE(7,*) IDUM,ICALJP(L),IQJP(L),JQJP(L),KQJP(L),
     &  NPORTJP(L),XJETL(L),YJETL(L),ZJET(L),PHJET(L),THJET(L),
     &                   DJET(L),CFRD(L),DJPER(L)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C28*  READ JET/PLUME SOURCE LOCATIONS AND PARAMETERS
C
      NCARD=28
      DO NSKIP=1,21
      READ(1,1)
      ENDDO
      DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) IDUM,NJEL(L),NJPMX(L),ISENT(L),ISTJP(L),
     &   NUDJP(L),IOUTJP(L),NZPRJP(L),ISDJP(L),IUPCJP(L),
     &   JUPCJP(L),KUPCJP(L)
      WRITE(7,1002)NCARD
      WRITE(7,*) IDUM,NJEL(L),NJPMX(L),ISENT(L),ISTJP(L),
     &   NUDJP(L),IOUTJP(L),NZPRJP(L),ISDJP(L),IUPCJP(L),
     &   JUPCJP(L),KUPCJP(L)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C29*  READ ADDITIONAL JET/PLUME PARAMETERS
C
      NCARD=29
      DO NSKIP=1,17
      READ(1,1)
      ENDDO
      DO L=1,NQJPIJ
       READ(1,*,IOSTAT=ISO) IDUM,QQCJP(L),NQSERJP(L),NQWRSERJP(L),
     &   NCSERJP(L,1),NCSERJP(L,2),NCSERJP(L,3),
     &   NCSERJP(L,4),NTXSRJP,NSDSRJP,NSNSRJP
       WRITE(7,1002)NCARD
       WRITE(7,*) IDUM,QQCJP(L),NQSERJP(L),NQWRSERJP(L),
     &   NCSERJP(L,1),NCSERJP(L,2),NCSERJP(L,3),
     &   NCSERJP(L,4),NTXSRJP,NSDSRJP,NSNSRJP
       NUDJPC(L)=1
       IF(ISO.GT.0) GOTO 100
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERJP(L,M)=NTXSRJP
       ENDDO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERJP(L,M)=NSDSRJP
       ENDDO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERJP(L,M)=NSNSRJP
       ENDDO
	 IF(ICALJP(L).EQ.2)THEN
	   QWRCJP(L)=QQCJP(L)
	   QQCJP(L)=0.
	 ELSE
	   QWRCJP(L)=0.
	 ENDIF
      ENDDO
      IF(NQJPIJ.GT.1)THEN
        DO L=2,NQJPIJ
         NUDJP(L)=NUDJP(1)
        ENDDO
      ENDIF
C
C30*  READ TIME CONSTANT INFLOW CONCENTRATIONS FOR TIME CONSTANT
C     JET/PLUME SOURCES SAL,TEM,DYE,SFL,TOX(1 TO NOTX)
C
      NCARD=30
      DO NSKIP=1,12
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) (CQSE(M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQSE(M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
	IF(ICALJP(L).EQ.1)THEN
       DO MS=1,MMAX
	  CWRCJP(L,MS)=0.0
        DO K=1,KC
        CQCJP(K,L,MS)=CQSE(MS)
        ENDDO
       ENDDO
	ELSE
       DO MS=1,MMAX
	  CWRCJP(L,MS)=CQSE(MS)
        DO K=1,KC
        CQCJP(K,L,MS)=0.0
        ENDDO
       ENDDO
	ENDIF
      ENDDO
C
C31*  READ TIME CONSTANT INFLOW CONCENTRATIONS FOR TIME CONSTANT
C     JET/PLUME SOURCES SED(1 TO NSED),SND(1 TO NSND)
C
      NCARD=31
      DO NSKIP=1,13
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) (CQSE(M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQSE(M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
	IF(ICALJP(L).EQ.1)THEN
       DO MS=MMIN,MMAX
	  CWRCJP(L,MS)=0.
        DO K=1,KC
        CQCJP(K,L,MS)=CQSE(MS)
        ENDDO
       ENDDO
	ELSE
       DO MS=MMIN,MMAX
	  CWRCJP(L,MS)=CQSE(MS)
        DO K=1,KC
        CQCJP(K,L,MS)=0.
        ENDDO
       ENDDO
	ENDIF
      ENDDO
C
C32*  READ SURF ELEV OR PRESS DEPENDENT FLOW CONTROL STRUCTURE INFO
C
      NCARD=32
      DO NSKIP=1,32
      READ(1,1)
      ENDDO
      DO L=1,NQCTL
      READ(1,*,IOSTAT=ISO)IQCTLU(L),JQCTLU(L),IQCTLD(L),JQCTLD(L),
     &                   NQCTYP(L),NQCTLQ(L),NQCMUL(L),NQCMFU(L),
     &                   NQCMFD(L),BQCMFU(L),BQCMFD(L)
      WRITE(7,1002)NCARD
      WRITE(7,*)IQCTLU(L),JQCTLU(L),IQCTLD(L),JQCTLD(L),
     &                   NQCTYP(L),NQCTLQ(L),NQCMUL(L),NQCMFU(L),
     &                   NQCMFD(L),BQCMFU(L),BQCMFD(L)
C     &                   NQCMFD(L),IQCAX(L),JQCAX(L)
      IF(ISO.GT.0) GOTO 100
       DO K=1,KC
         QCTLTO(K,L)=0.
         QCTLT(K,L)=0.
       ENDDO
      ENDDO
C
C33*  READ FLOW WITHDRAWAL, HEAT OR MATERIAL ADDITION, FLOW RETURN DATA
C
      NCARD=33
      DO NSKIP=1,27
      READ(1,1)
      ENDDO
      DO L=1,NQWR
      READ(1,*,IOSTAT=ISO)IQWRU(L),JQWRU(L),KQWRU(L),
     &                    IQWRD(L),JQWRD(L),KQWRD(L),QWR(L),
     &    NQWRSERQ(L),NQWRMFU(L),NQWRMFD(L),BQWRMFU(L),BQWRMFD(L),
     &    ANGWRMFD(L)
      WRITE(7,1002)NCARD
      WRITE(7,*)IQWRU(L),JQWRU(L),KQWRU(L),
     &                     IQWRD(L),JQWRD(L),KQWRD(L),QWR(L),
     &    NQWRSERQ(L),NQWRMFU(L),NQWRMFD(L),BQWRMFU(L),BQWRMFD(L),
     &    ANGWRMFD(L)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C34*  READ TIME CONSTANT WITHDRAWAL,ADD,RETURN, CONCENTRATION INCREASES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=34
      DO NSKIP=1,10
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NQWR
      READ(1,*,IOSTAT=ISO) (CQWR(L,MS),MS=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQWR(L,MS),MS=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C35*  READ TIME CONSTANT WITHDRAWAL,ADD,RETURN, CONCENTRATION INCREASES
C     SED(1 TO NSED),SND(1 TO NSND)
C
      NCARD=35
      DO NSKIP=1,7
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NQWR
      READ(1,*,IOSTAT=ISO) (CQWR(L,MS),MS=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQWR(L,MS),MS=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C36*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
C
      NCARD=36
      DO NSKIP=1,42
      READ(1,1)
      ENDDO
      IF(NSED.GT.0.OR.NSND.GT.0)THEN
      READ(1,*,IOSTAT=ISO)ISEDINT,ISEDBINT,ISEDWC,ISMUD,ISNDWC,ISEDVW,
     &                     ISNDVW,KB,ISDTXBUG
      WRITE(7,1002)NCARD
      WRITE(7,*)ISEDINT,ISEDBINT,ISEDWC,ISMUD,ISNDWC,ISEDVW,
     &                     ISNDVW,KB,ISDTXBUG
      IF(ISO.GT.0) GOTO 100
      ENDIF
C
C36A*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
C
      COEFTSBL=4.0
C
      NCARD=361
      DO NSKIP=1,29
      READ(1,1)
      ENDDO
      IF(NSED.GT.0.OR.NSND.GT.0)THEN
      READ(1,*,IOSTAT=ISO)ISBEDSTR,ISBSDIAM,ISBSDFUF,COEFTSBL,VISMUDST,
     &                    ISBKERO
      WRITE(7,1002)NCARD
      WRITE(7,*)ISBEDSTR,ISBSDIAM,ISBSDFUF,COEFTSBL,VISMUDST,ISBKERO
      IF(ISO.GT.0) GOTO 100
	ENDIF
C
C36B*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
C
      NCARD=362
      DO NSKIP=1,19
      READ(1,1)
      ENDDO
      IF(NSED.GT.0.OR.NSND.GT.0)THEN
      READ(1,*,IOSTAT=ISO)ISEDAL,ISNDAL,IALTYP,IALSTUP,
     &                     ISEDEFF,HBEDAL,COEHEFF,COEHEFF2
      WRITE(7,1002)NCARD
      WRITE(7,*)ISEDAL,ISNDAL,IALTYP,IALSTUP,
     &                     HBEDAL,COEHEFF,COEHEFF2
      IF(ISO.GT.0) GOTO 100
      ENDIF
C
C37*  BED MECHANICAL PROPERTIES PARAMETER SET 1
C
      NCARD=37
      DO NSKIP=1,32
        READ(1,1)
      ENDDO
      IF(NSED.GT.0.OR.NSND.GT.0)THEN
        READ(1,*,IOSTAT=ISO)ISEDDT,IBMECH,IMORPH,HBEDMAX,BEDPORC,
     &                      SEDMDMX,SEDMDMN,SEDVDRD,SEDVDRM,SEDVRDT
        WRITE(7,1002)NCARD
        WRITE(7,*)ISEDDT,IBMECH,IMORPH,HBEDMAX,BEDPORC,
     &                      SEDMDMX,SEDMDMN,SEDVDRD,SEDVDRM,SEDVRDT
        IF(ISO.GT.0) GOTO 100
        ISEDDTC=0
        IF(IBMECH.EQ.0) THEN
          SNDVDRD=BEDPORC/(1.-BEDPORC)
          SEDVDRD=BEDPORC/(1.-BEDPORC)
          SEDVDRM=SEDVDRD
        END IF
        IF(IBMECH.GE.1) THEN
          ISEDBINT=1
          IF(ISEDINT.EQ.1) ISEDINT=3
          ISEDINT=MAX(ISEDINT,2)
        END IF
      ENDIF
C
      SNDVDRD=BEDPORC/(1.-BEDPORC)
C
	DO NS=1,NSED
	  VDRDEPO(NS)=SEDVDRD
	ENDDO
C
	DO NS=1,NSND
	  NX=NS+NSED
	  VDRDEPO(NX)=SNDVDRD
	ENDDO
C
C38*  BED MECHANICAL PROPERTIES PARAMETER SET 2
C
      NCARD=38
      DO NSKIP=1,16
      READ(1,1)
      ENDDO
      IF(NSED.GT.0.OR.NSND.GT.0)THEN
      READ(1,*,IOSTAT=ISO)IBMECHK,BMECH1,BMECH2,BMECH3,BMECH4,BMECH5,
     &                     BMECH6
      WRITE(7,1002)NCARD
      WRITE(7,*)IBMECHK,BMECH1,BMECH2,BMECH3,BMECH4,BMECH5,BMECH6
      IF(ISO.GT.0) GOTO 100
      ENDIF
C
C39*  READ COHESIVE SEDIMENT PARAMETER SET 1 REPEAT DATA LINE NSED TIMES
C
      NCARD=39
      DO NSKIP=1,24
      READ(1,1)
      ENDDO
      IF(NSED.GT.0)THEN
      DO N=1,NSED
       READ(1,*,IOSTAT=ISO)SEDO(N),SEDBO(N),SDEN(N),SSG(N),
     &        WSEDO(N),SEDN(N),SEXP(N),TAUD(N),ISEDSCOR(N),ISPROBDEP(N)
       WRITE(7,1002)NCARD
      WRITE(7,*)SEDO(N),SEDBO(N),SDEN(N),SSG(N),
     &        WSEDO(N),SEDN(N),SEXP(N),TAUD(N),ISEDSCOR(N),ISPROBDEP(N)
       IF(ISO.GT.0) GOTO 100
       SEDDIA(N)=0.
      ENDDO
      ENDIF
C
      IF(NSED.GT.0)THEN
        DO L=1,LC
	    TAUDS(L)=TAUD(1)
        ENDDO
	ENDIF
C
C40*  READ COHESIVE SEDIMENT PARAMETER SET 2 REPEAT DATA LINE NSED TIMES
C
      NCARD=40
      DO NSKIP=1,29
      READ(1,1)
      ENDDO
      IF(NSED.GT.0)THEN
      DO N=1,NSED
       READ(1,*,IOSTAT=ISO)IWRSP(N),IWRSPB(N),WRSPO(N),TAUR(N),TAUN(N),
     &      TEXP(N),VDRRSPO(N),COSEDHID(N)
       WRITE(7,1002)NCARD
       WRITE(7,*)IWRSP(N),IWRSPB(N),WRSPO(N),TAUR(N),TAUN(N),TEXP(N),
     &           VDRRSPO(N),COSEDHID(N)
       IF(ISO.GT.0) GOTO 100
C
       ISNDEQ(N)=0
      ENDDO
      ENDIF
C
C     HOUSATONIC VERSION FOR CRITICAL STRESS
C     HQI Change, 08/25/03, SO and RM
C     Change to implement critical shear stress options
C
      IF(ISHOUSATONIC.EQ.1)THEN
      IF(NSED.GT.0)THEN
       if(IWRSP(N).EQ.999) then
          open(1001,file='TAU_CRIT_COH.INP',status='OLD')
          read(1001,*,IOSTAT=ISO)KBINPUT
          do l = 2, LC-1
             read(1001,*,IOSTAT=ISO) (TAUCRCOH(l,k),k=1,KBINPUT)
             do k=KBINPUT+1,KB
                TAUCRCOH(l,k)=TAUCRCOH(l,KBINPUT)
             enddo
c             write(7,7777)l,(TAUCRCOH(l,k),k=1,KB)
          enddo
          close(1001)
          IF(ISO.GT.0) GOTO 100
       endif
	ENDIF
	ENDIF
c
 7777 format(i6,20e12.4)
C
C     SITE SPECIFIC RESUSPENSION INFORMATION BASED ON HAMRICK'S
C     ANALYSIS OF SEDFLUME CORES JANUARY 2008
C
      IF(NSED.GT.0)THEN
      IF(IWRSP(N).GE.99)THEN
C
        OPEN(99,FILE='SSCOHSEDPROP.INP')
        DO NSKIP=1,12
	    READ(99,1)
        ENDDO
	  READ(99,*,IOSTAT=ISO)NCOHSEDS,NCOHSEDL,RMULADJ
	  DO NSITES=1,NCOHSEDS
	    READ(99,*,IOSTAT=ISO)IDUM,RUMLADJC,TAUDSS(NSITES)
	    READ(99,*,IOSTAT=ISO)IDUM,(WRSPOSS(NL,NSITES),NL=1,NCOHSEDL)
          DO NL=1,NCOHSEDL
	      WRSPOSS(NL,NSITES)=RUMLADJC*RMULADJ*WRSPOSS(NL,NSITES)
          ENDDO
	    READ(99,*,IOSTAT=ISO)IDUM,(TAURSS(NL,NSITES),NL=1,NCOHSEDL)
	    READ(99,*,IOSTAT=ISO)IDUM,(TAUNSS(NL,NSITES),NL=1,NCOHSEDL)
          READ(99,*,IOSTAT=ISO)IDUM,(TEXPSS(NL,NSITES),NL=1,NCOHSEDL)
	    READ(99,*,IOSTAT=ISO)IDUM,(DEPBBSS(NL,NSITES),NL=1,NCOHSEDL)
        ENDDO
        CLOSE(99)
        IF(ISO.GT.0) THEN
	    NCARD=4099
          GOTO 100
        ENDIF

      ENDIF
	ENDIF
C
C41*  READ NONCOHESIVE SEDIMENT PARAMETER SET 1 REPEAT DATA LINE NSND TIMES
C
      NCARD=41
      DO NSKIP=1,20
      READ(1,1)
      ENDDO
      IF(NSND.GT.0)THEN
      DO NX=1,NSND
       N=NX+NSED
       READ(1,*,IOSTAT=ISO)SEDO(N),SEDBO(N),SDEN(N),SSG(N),SEDDIA(N),
     &                  WSEDO(N),SEDN(N),SEXP(N),TAUD(N),ISEDSCOR(N)
C
C IF SETTLING VELOCITY IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
C
       IF(WSEDO(N).LT.0.0)THEN
         WSEDO(N)=SETSTVEL(SEDDIA(N),SSG(N))
       ENDIF
C
       WRITE(7,1002)NCARD
       WRITE(7,*)SEDO(N),SEDBO(N),SDEN(N),SSG(N),SEDDIA(N),
     &                  WSEDO(N),SEDN(N),SEXP(N),TAUD(N),ISEDSCOR(N)
       IF(ISO.GT.0) GOTO 100
      ENDDO
      ENDIF
C
C42*  READ NONCOHESIVE SEDIMENT PARAMETER SET 2 REPEAT DATA LINE NSND TIMES
C
      NCARD=42
      DO NSKIP=1,44
      READ(1,1)
      ENDDO
      IF(NSND.GT.0)THEN
      DO NX=1,NSND
       N=NX+NSED
       READ(1,*,IOSTAT=ISO)ISNDEQ(N),ISBDLD(N),TAUR(N),TAUN(N),
     &                     TCSHIELDS(N),ISLTAUC(N),IBLTAUC(N),
     &                     IROUSE(NX),ISNDM1(NX),ISNDM2(NX),RSNDM(NX)
C     &                     TEXP(N),ISLTAUC(N),IBLTAUC(N)
C
C IF TAUR(N) IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
C
C       TAUR:     CRITICAL SHIELDS STRESS IN (M/S)**2   (ISNDEQ=2)
C       TAUN:     EQUAL TO TAUR FOR NONCHOESIVE SED TRANS  (ISNDEQ=2)
C       TEXP:     CRITICAL SHIELDS PARAMETER  (ISNDEQ=2)
C
       DSTAR=0.0
	 USTAR=0.0
       IF(TAUR(N).LT.0.0)THEN
C         CALL SETSHLD(TAUR(N),TEXP(N),SEDDIA(N),SSG(N),DSTAR,USTAR)
         CALL SETSHLD(TAUR(N),TCSHIELDS(N),SEDDIA(N),SSG(N),DSTAR,USTAR)
         TAUN(N)=TAUR(N)
       ENDIF
C
C IF TAUR(N) IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
C
       WRITE(7,1002)NCARD
C       WRITE(7,*)ISNDEQ(N),TAUR(N),TAUN(N),TEXP(N),SEDDIA(N),SSG(N),
       WRITE(7,*)ISNDEQ(N),TAUR(N),TAUN(N),TCSHIELDS(N),SEDDIA(N),
     &           SSG(N),DSTAR,USTAR
       IF(ISO.GT.0) GOTO 100
       IWRSP(N)=0
       WRSPO(N)=0
      ENDDO
      ENDIF
C
C42A*  READ NONCOHESIVE SEDIMENT BED LOAD PARAMETERS
C
      NCARD=42
      DO NSKIP=1,19
      READ(1,1)
      ENDDO
      IF(NSND.GT.0)THEN
      DO NS=1,NSND
       READ(1,*,IOSTAT=ISO)ISBDLDBC,SBDLDA(NS),SBDLDB(NS),
     &           SBDLDG1(NS),SBDLDG2(NS),SBDLDG3(NS),SBDLDG4(NS),
     &             SBDLDP(NS),ISBLFUC,BLBSNT
       WRITE(7,1002)NCARD
       WRITE(7,*)ISBDLDBC,SBDLDA(NS),SBDLDB(NS),
     &           SBDLDG1(NS),SBDLDG2(NS),SBDLDG3(NS),SBDLDG4(NS),
     &             SBDLDP(NS),ISBLFUC,BLBSNT
       IF(ISO.GT.0) GOTO 100
      ENDDO
C      DO NX=1,NSND
C       N=NX+NSED
c       READ(1,*,IOSTAT=ISO)ISBDLDBC,SBDLDA,SBDLDB,SBDLDG1,
c     &           SBDLDG2,SBDLDG3,SBDLDG4,SBDLDP,ISBLFUC,BLBSNT
c       WRITE(7,1002)NCARD
c       WRITE(7,*)ISBDLDBC,SBDLDA,SBDLDB,SBDLDG1,SBDLDG2,SBDLDG3,
c     &           SBDLDG4,SBDLDP,ISBLFUC,BLBSNT
c       IF(ISO.GT.0) GOTO 100
C      ENDDO
      ENDIF
C
C43*  READ TOXIC CONTAMINANT INITIAL CONDITIONS AND PARAMETERS
C
      NCARD=43
      DO NSKIP=1,23
      READ(1,1)
      ENDDO
      IF(NTOX.GT.0)THEN
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM,ITXINT(NT),ITXBDUT(NT),TOXINTW(NT),
     &      TOXINTB(NT),RKTOXW(NT),TKTOXW(NT),RKTOXB(NT),TRTOXB(NT)
        WRITE(7,1002)NCARD
      WRITE(7,*)NDUM,ITXINT(NT),ITXBDUT(NT),TOXINTW(NT),
     &      TOXINTB(NT),RKTOXW(NT),TKTOXW(NT),RKTOXB(NT),TRTOXB(NT)
        IF(ISO.GT.0) GOTO 100
      ENDDO
      ENDIF
C
C44*  READ TOXIC CONTAMINANT PARAMETERS
C
      NCARD=44
      DO NSKIP=1,24
      READ(1,1)
      ENDDO
      IF(NTOX.GT.0)THEN
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM,ISTOC(NT),VOLTOX(NT),RMOLTX(NT),
     &                       RKTOXP(NT),SKTOXP(NT),DIFTOX(NT),
     &                       DIFTOXS(NT),PDIFTOX (NT),DPDIFTOX(NT)
        WRITE(7,1002)NCARD
      WRITE(7,*)NDUM,ISTOC(NT),VOLTOX(NT),RMOLTX(NT),RKTOXP(NT),
     &                       SKTOXP(NT),DIFTOX(NT),
     &                       DIFTOXS(NT),PDIFTOX (NT),DPDIFTOX(NT)
        IF(ISO.GT.0) GOTO 100
        ISPMXZ(NT)=0
        IF(PDIFTOX(NT).LT.0) ISPMXZ(NT)=1
	  ISDIFBW(NT)=0
	  IF(DIFTOXS(NT).LT.0.0)THEN
	    DIFTOXS(NT)=ABS(DIFTOXS(NT))
	    ISDIFBW(NT)=1
	  ENDIF
      ENDDO
      ENDIF
C
C44A*  READ TOXIC CONTAMINANT PARAMETERS
C
      NCARD=441
      DO NSKIP=1,13
      READ(1,1)
      ENDDO
      IF(NTOX.GT.0)THEN
        READ(1,*,IOSTAT=ISO)IADTOXDP,IADTOXCOR,ISTOXALL,NSTOXALL
        WRITE(7,1002)NCARD
        WRITE(7,*)IADTOXDP,IADTOXCOR,ISTOXALL,NSTOXALL
        IF(ISO.GT.0) GOTO 100
      ENDIF
C
      JSTOXALL=0
      IF(NSTOXALL.GT.0) MSTOXALL=NTSPTC/NSTOXALL
C
C
C45*  READ TOXIC CONTAMINANT-SEDIMENT INTERACTION PARAMETERS
C
      NCARD=45
      DO NSKIP=1,19
      READ(1,1)
      ENDDO
      IF(NTOX.GT.0)THEN
      DO NT=1,NTOX
       IF(NSED.GT.0)THEN
         DO N=1,NSED
         READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,
     &         ITXPARW(N,NT),TOXPARW(N,NT),CONPARW(N,NT),
     &         ITXPARB(N,NT),TOXPARB(N,NT),CONPARB(N,NT)
         WRITE(7,1002)NCARD
         WRITE(7,*)NDUM1,NDUM2,
     &         ITXPARW(N,NT),TOXPARW(N,NT),CONPARW(N,NT),
     &         ITXPARB(N,NT),TOXPARB(N,NT),CONPARB(N,NT)
         IF(ISO.GT.0) GOTO 100
         ENDDO
       ENDIF
       IF(NSND.GT.0)THEN
         DO NX=1,NSND
         N=NX+NSED
         READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,
     &         ITXPARW(N,NT),TOXPARW(N,NT),CONPARW(N,NT),
     &         ITXPARB(N,NT),TOXPARB(N,NT),CONPARB(N,NT)
         WRITE(7,1002)NCARD
         WRITE(7,*)NDUM1,NDUM2,
     &         ITXPARW(N,NT),TOXPARW(N,NT),CONPARW(N,NT),
     &         ITXPARB(N,NT),TOXPARB(N,NT),CONPARB(N,NT)
         IF(ISO.GT.0) GOTO 100
         ENDDO
       ENDIF
      ENDDO
      ENDIF
C
C45A*  READ TOXIC CONTAMINANT ORGANIC CARBON PARAMETERS
C
      NCARD=451
      DO NSKIP=1,25
      READ(1,1)
      ENDDO
      IF(NTOX.GT.0)THEN
c	  NTMP=0
c	  DO NT=1,NTOX
c	    IF(ISTOC(NT).GT.0)NTMP=1
c	  ENDDO
c	  IF(NTMP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,
     & 	  STDOCWC,STPOCWC,STDOCBC,STPOCBC
          WRITE(7,1002)NCARD
          WRITE(7,*)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,
     & 	  STDOCWC,STPOCWC,STDOCBC,STPOCBC
          IF(ISO.GT.0) GOTO 100
c        ENDIF
      ENDIF
C

C45B* READ TOXIC CONTAMINANT-ORGANIC CARBON INTERACTION PARAMETERS
C
      NCARD=452
      DO NSKIP=1,19
      READ(1,1)

      ENDDO
      IF(NTOX.GT.0)THEN
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,
     &         ITXPARWC(1,NT),TOXPARWC(1,NT),CONPARWC(1,NT),
     &         ITXPARBC(1,NT),TOXPARBC(1,NT),CONPARBC(1,NT)
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM1,NDUM2,
     &         ITXPARWC(1,NT),TOXPARWC(1,NT),CONPARWC(1,NT),
     &         ITXPARBC(1,NT),TOXPARBC(1,NT),CONPARBC(1,NT)
        IF(ISO.GT.0) GOTO 100
        READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,
     &         ITXPARWC(2,NT),TOXPARWC(2,NT),CONPARWC(2,NT),
     &         ITXPARBC(2,NT),TOXPARBC(2,NT),CONPARBC(2,NT)
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM1,NDUM2,
     &         ITXPARWC(2,NT),TOXPARWC(2,NT),CONPARWC(2,NT),
     &         ITXPARBC(2,NT),TOXPARBC(2,NT),CONPARBC(2,NT)
        IF(ISO.GT.0) GOTO 100
      ENDDO
      ENDIF
C
C45C* READ TOXIC CONTAMINANT-ORGANIC CARBON WATER COLUMN POC FRACTIONS
C
      NCARD=453
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      IF(NTOX.GT.0)THEN
      NTMP=NSED+NSND
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM1,(FPOCWST(NS,NT),NS=1,NTMP)
        IF(ISO.GT.0) GOTO 100
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM1,(FPOCWST(NS,NT),NS=1,NTMP)
      ENDDO
      ENDIF
C
C reset inorganic sediment partition coefficients based on
C fraction of POC assigned to each sediment class in water columnn
C

      DO NT=1,NTOX
cjmh040203        IF(ISTOC(NT).GE.1.AND.ISTOC(NT).LE.2)THEN
        IF(ISTOC(NT).GE.2)THEN
          IF(NSED.GT.0)THEN
            DO NS=1,NSED
              ITXPARW(NS,NT)=0
C              TOXPARW(NS,NT)=FPOCWST(NS,NT)*TOXPARWC(2,NT)
              TOXPARW(NS,NT)=TOXPARWC(2,NT)
              CONPARW(NS,NT)=0.
            ENDDO
          ENDIF
          IF(NSND.GT.0)THEN
            DO NX=1,NSND
              NS=NSED+NX
              ITXPARW(NS,NT)=0
C              TOXPARW(NS,NT)=FPOCWST(NS,NT)*TOXPARWC(2,NT)
              TOXPARW(NS,NT)=TOXPARWC(2,NT)
              CONPARW(NS,NT)=0.
            ENDDO
          ENDIF
        ENDIF
      ENDDO
C
C45D* READ TOXIC CONTAMINANT-ORGANIC CARBON SED BED POC FRACTIONS
C
      NCARD=454
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      IF(NTOX.GT.0)THEN
      NTMP=NSED+NSND
      WRITE(7,1002)NCARD
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM1,(FPOCBST(NS,NT),NS=1,NTMP)
        IF(ISO.GT.0) GOTO 100
        WRITE(7,*)NDUM1,(FPOCBST(NS,NT),NS=1,NTMP)
      ENDDO
      ENDIF
C
C reset inorganic sediment partition coefficients based on
C fraction of POC assigned to each sediment class in sediment bed
C
      DO NT=1,NTOX
cjmn040203        IF(ISTOC(NT).GE.1.AND.ISTOC(NT).LE.2)THEN
        IF(ISTOC(NT).GE.2)THEN
          IF(NSED.GT.0)THEN
            DO NS=1,NSED
              ITXPARB(NS,NT)=0
C              TOXPARB(NS,NT)=FPOCBST(NS,NT)*TOXPARBC(2,NT)
              TOXPARB(NS,NT)=TOXPARBC(2,NT)
              CONPARB(NS,NT)=0
            ENDDO
          ENDIF
          IF(NSND.GT.0)THEN
            DO NX=1,NSND
              NS=NSED+NX
              ITXPARB(NS,NT)=0
C              TOXPARB(NS,NT)=FPOCBST(NS,NT)*TOXPARBC(2,NT)
              TOXPARB(NS,NT)=TOXPARBC(2,NT)
              CONPARB(NS,NT)=0
            ENDDO
          ENDIF
        ENDIF
      ENDDO
C
C46*  READ BUOYANCY, TEMPERATURE, DYE DATA AND CONCENTRATION BC DATA
C
      NCARD=46
      DO NSKIP=1,20
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)BSC,TEMO,HEQT,ISBEDTEMI,KBH,RKDYE,
     &                    NCBS,NCBW,NCBE,NCBN
      WRITE(7,1002)NCARD
      WRITE(7,*)BSC,TEMO,HEQT,ISBEDTEMI,KBH,RKDYE,NCBS,NCBW,NCBE,NCBN
      IF(ISO.GT.0) GOTO 100
      IF(BSC.EQ.2)THEN
        BSC=1.
        IBSC=1
       ELSE
        IBSC=0
      ENDIF
C
C46a*  READ ICE EFFECTS
C
      NCARD=461
      DO NSKIP=1,24
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)ISICE,ISICECOV,ISICETHK,NISER,ICETHKFUN,
     &     DYICEBEG,DYICEEND,DYICEM1,DYICEM2,RICETHKMAX,TEMPICE
      WRITE(7,1002)NCARD
      WRITE(7,*)ISICE,ISICECOV,ISICETHK,NISER,ICETHKFUN,DYICEBEG,
     &     DYICEEND,DYICEM1,DYICEM2,RICETHKMAX,TEMPICE
      IF(ISO.GT.0) GOTO 100
C
C47*  READ LOCATIONS OF CONC BC'S ON SOUTH BOUNDARIES
C
      NCARD=47
      DO NSKIP=1,16
      READ(1,1)
      ENDDO
      DO L=1,NCBS
       READ(1,*,IOSTAT=ISO) ICBS(L),JCBS(L),NTSCRS(L),
     &   NCSERS(L,1),NCSERS(L,2),NCSERS(L,3),NCSERS(L,4),
     &   NTOXSRC,NSEDSRC,NSNDSRC
       WRITE(7,1002)NCARD
      WRITE(7,*) ICBS(L),JCBS(L),NTSCRS(L),
     &   NCSERS(L,1),NCSERS(L,2),NCSERS(L,3),NCSERS(L,4),
     &   NTOXSRC,NSEDSRC,NSNDSRC
       IF(ISO.GT.0) GOTO 100
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERS(L,M)=NTOXSRC
       ENDDO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERS(L,M)=NSEDSRC
       ENDDO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERS(L,M)=NSNDSRC
       ENDDO
      ENDDO
C
C48*  READ CONSTANT BOTTOM CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=48
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) (CBS(L,1,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBS(L,1,M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C49*  READ CONSTANT BOTTOM CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      NCARD=49
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) (CBS(L,1,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBS(L,1,M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C50*  READ CONSTANT SURFACE CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=50
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) (CBS(L,2,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBS(L,2,M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C51*  READ CONSTANT SURFACE CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      NCARD=51
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) (CBS(L,2,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBS(L,2,M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C52*  READ LOCATIONS OF CONC BC'S ON WEST BOUNDARIES
C
      NCARD=52
      DO NSKIP=1,16
      READ(1,1)
      ENDDO
      DO L=1,NCBW
       READ(1,*,IOSTAT=ISO) ICBW(L),JCBW(L),NTSCRW(L),
     &   NCSERW(L,1),NCSERW(L,2),NCSERW(L,3),NCSERW(L,4),
     &   NTOXSRC,NSEDSRC,NSNDSRC
       WRITE(7,1002)NCARD
      WRITE(7,*) ICBW(L),JCBW(L),NTSCRW(L),
     &   NCSERW(L,1),NCSERW(L,2),NCSERW(L,3),NCSERW(L,4),
     &   NTOXSRC,NSEDSRC,NSNDSRC
       IF(ISO.GT.0) GOTO 100
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERW(L,M)=NTOXSRC
       ENDDO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERW(L,M)=NSEDSRC
       ENDDO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERW(L,M)=NSNDSRC
       ENDDO
      ENDDO
C
C53*  READ CONSTANT BOTTOM CONCENTRATION ON WEST CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=53
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) (CBW(L,1,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBW(L,1,M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C54*  READ CONSTANT BOTTOM CONCENTRATION ON WEST CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      NCARD=54
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) (CBW(L,1,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBW(L,1,M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C55*  READ CONSTANT SURFACE CONCENTRATION ON WEST CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=55
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) (CBW(L,2,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBW(L,2,M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C56*  READ CONSTANT SURFACE CONCENTRATION ON WEST CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      NCARD=56
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) (CBW(L,2,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBW(L,2,M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C57*  READ LOCATIONS OF CONC BC'S ON EAST BOUNDARIES
C
      NCARD=57
      DO NSKIP=1,16
      READ(1,1)
      ENDDO
      DO L=1,NCBE
       READ(1,*,IOSTAT=ISO) ICBE(L),JCBE(L),NTSCRE(L),
     &   NCSERE(L,1),NCSERE(L,2),NCSERE(L,3),NCSERE(L,4),
     &   NTOXSRC,NSEDSRC,NSNDSRC
       WRITE(7,1002)NCARD
      WRITE(7,*) ICBE(L),JCBE(L),NTSCRE(L),
     &   NCSERE(L,1),NCSERE(L,2),NCSERE(L,3),NCSERE(L,4),
     &   NTOXSRC,NSEDSRC,NSNDSRC
       IF(ISO.GT.0) GOTO 100
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERE(L,M)=NTOXSRC
       ENDDO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERE(L,M)=NSEDSRC
       ENDDO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERE(L,M)=NSNDSRC
       ENDDO
      ENDDO
C
C58*  READ CONSTANT BOTTOM CONCENTRATION ON EAST CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=58
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) (CBE(L,1,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBE(L,1,M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C59*  READ CONSTANT BOTTOM CONCENTRATION ON EAST CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      NCARD=59
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) (CBE(L,1,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBE(L,1,M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C60*  READ CONSTANT SURFACE CONCENTRATION ON EAST CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=60
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) (CBE(L,2,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBE(L,2,M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C61*  READ CONSTANT SURFACE CONCENTRATION ON EAST CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      NCARD=61
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) (CBE(L,2,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBE(L,2,M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C62*  READ LOCATIONS OF CONC BC'S ON NORTH BOUNDARIES
C
      NCARD=62
      DO NSKIP=1,16
      READ(1,1)
      ENDDO
      DO L=1,NCBN
       READ(1,*,IOSTAT=ISO) ICBN(L),JCBN(L),NTSCRN(L),
     &   NCSERN(L,1),NCSERN(L,2),NCSERN(L,3),NCSERN(L,4),
     &   NTOXSRC,NSEDSRC,NSNDSRC
       WRITE(7,1002)NCARD
      WRITE(7,*) ICBN(L),JCBN(L),NTSCRN(L),
     &   NCSERN(L,1),NCSERN(L,2),NCSERN(L,3),NCSERN(L,4),
     &   NTOXSRC,NSEDSRC,NSNDSRC
       IF(ISO.GT.0) GOTO 100
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERN(L,M)=NTOXSRC
       ENDDO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERN(L,M)=NSEDSRC
       ENDDO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERN(L,M)=NSNDSRC
       ENDDO
      ENDDO
C
C63*  READ CONSTANT BOTTOM CONCENTRATION ON NORTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=63
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) (CBN(L,1,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBN(L,1,M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C64*  READ CONSTANT BOTTOM CONCENTRATION ON NORTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      NCARD=64
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) (CBN(L,1,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBN(L,1,M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C65*  READ CONSTANT SURFACE CONCENTRATION ON NORTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      NCARD=65
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      MMAX=4+NTOX
      DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) (CBN(L,2,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBN(L,2,M),M=1,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C66*  READ CONSTANT SURFACE CONCENTRATION ON NORTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      NCARD=66
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) (CBN(L,2,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBN(L,2,M),M=MMIN,MMAX)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C66A*  READ CONCENTRATION DATA ASSIMILATION PARAMETERS
C
      NCARD=66
      DO NSKIP=1,8
       READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) NLCDA,TSCDA,(ISCDA(K),K=1,7)
      WRITE(7,1002)NCARD
      WRITE(7,*) NLCDA,TSCDA,(ISCDA(K),K=1,7)
      IF(ISO.GT.0) GOTO 100
C
C66B*  READ CONCENTRATION DATA ASSIMILATION LOCATIONS AND
C      SERIES IDENTIFIERS
C
      NCARD=66
      DO NSKIP=1,15
       READ(1,1)
      ENDDO
      IF(NLCDA.GT.0)THEN
        WRITE(7,1002)NCARD
        DO L=1,NLCDA
         READ(1,*,IOSTAT=ISO) ITPCDA(L),ICDA(L),JCDA(L),
     &               ICCDA(L),JCCDA(L),IDIRCDA(L),(NCSERA(L,K),K=1,7)
         WRITE(7,*)  ITPCDA(L),ICDA(L),JCDA(L),
     &               ICCDA(L),JCCDA(L),IDIRCDA(L),(NCSERA(L,K),K=1,7)
         IF(ISO.GT.0) GOTO 100
        ENDDO
      ENDIF
C
C67*  READ NEUTRALLY BUOYANT PARTICLE DRIFTER DATA
C
      NCARD=67
      DO NSKIP=1,22
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISPD,NPD,NPDRT,NWPD,ISLRPD,ILRPD1,ILRPD2,
     &                   JLRPD1, JLRPD2, MLRPDRT,IPLRPD
      WRITE(7,1002)NCARD
      WRITE(7,*) ISPD,NPD,NPDRT,NWPD,ISLRPD,ILRPD1,ILRPD2,
     &                   JLRPD1, JLRPD2, MLRPDRT,IPLRPD
      IF(ISO.GT.0) GOTO 100
C
C68*  READ NEUTRALLY BUOYANT PARTICLE INITIAL POSITIONS
C
      NCARD=68
      DO NSKIP=1,8
      READ(1,1)
      ENDDO
      DO NP=1,NPD
      READ(1,*,IOSTAT=ISO) RI(NP),RJ(NP),RK(NP)
      WRITE(7,1002)NCARD
      WRITE(7,*) RI(NP),RJ(NP),RK(NP)
      ENDDO
C
C69*  CONSTANTS FOR LONGITUDE AND LATITUDE OF CELL CENTERS
C
      NCARD=69
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) CDLON1,CDLON2,CDLON3,CDLAT1,CDLAT2,CDLAT3
      WRITE(7,1002)NCARD
      WRITE(7,*) CDLON1,CDLON2,CDLON3,CDLAT1,CDLAT2,CDLAT3
      IF(ISO.GT.0) GOTO 100
C
C70*  CONTROLS FOR WRITING ASCII OR BINARY DUMP FILES
C
      NCARD=70
      DO NSKIP=1,20
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)ISDUMP,ISADMP,NSDUMP,TSDUMP,TEDUMP,ISDMPP,
     &          ISDMPU,ISDMPW,ISDMPT,IADJDMP
      WRITE(7,1002)NCARD
      WRITE(7,*)ISDUMP,ISADMP,NSDUMP,TSDUMP,TEDUMP,ISDMPP,
     &          ISDMPU,ISDMPW,ISDMPT,IADJDMP
      IF(ISO.GT.0) GOTO 100
      JSDUMP=1
      NCDUMP=1
C
C71*  CONTROLS FOR HORIZONTAL PLANE SCALAR FIELD CONTOURING
C
      NCARD=71
      DO NSKIP=1,15
      READ(1,1)
      ENDDO
      DO N=1,7
       READ(1,*,IOSTAT=ISO) ISSPH(N),NPSPH(N),ISRSPH(N),ISPHXY(N)
       WRITE(7,1002)NCARD
      WRITE(7,*) ISSPH(N),NPSPH(N),ISRSPH(N),ISPHXY(N)
      ENDDO
      IF(ISO.GT.0) GOTO 100
C
      ISSPH(8)=0
	NPSPH(8)=0
      DO N=1,7
	  IF(ISSPH(N).GE.1.AND.ISPHXY(N).EQ.3)ISSPH(8)=1
	ENDDO
	IF(ISSPH(8).EQ.1)THEN
	  DO N=1,7
	    NPSPH(8)=MAX(NPSPH(N),NPSPH(8))
	  ENDDO
	ENDIF
C
C71A*  CONTROLS FOR HORIZONTAL PLANE SEDIMENT BED PROPERTIES
C
      NCARD=71
      DO NSKIP=1,22
      READ(1,1)
      ENDDO
       READ(1,*,IOSTAT=ISO) ISBPH,ISBEXP,NPBPH,ISRBPH,ISBBDN,ISBLAY,
     &                ISBPOR,ISBSED,ISBSND,ISBVDR,ISBARD
       WRITE(7,1002)NCARD
      WRITE(7,*) ISBPH,ISBEXP,NPBPH,ISRBPH,ISBBDN,ISBLAY,
     &                      ISBPOR,ISBSED,ISBSND,ISBVDR,ISBARD
      IF(ISO.GT.0) GOTO 100
	IF(ISBEXP.GE.1) NPSPH(8)=MAX(NPSPH(8),NPBPH)
	JSBPH=1
	JSBPHA=1
C
C71B*  CONTROLS FOR FOOD CHAIN MODEL OUTPUT
C
      NCARD=71
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
       READ(1,*,IOSTAT=ISO) ISFDCH,NFDCHZ,HBFDCH,TFCAVG
       WRITE(7,1002)NCARD
      WRITE(7,*) ISFDCH,NFDCHZ,HBFDCH,TFCAVG
      IF(ISO.GT.0) GOTO 100
	JSFDCH=1
C
C72*  CONTROLS FOR HORIZONTAL PLANE SURFACE ELEVATION OR PRESSURE
C     PLOTTING
C
      NCARD=72
      DO NSKIP=1,14
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISPPH,NPPPH,ISRPPH,IPPHXY
      WRITE(7,1002)NCARD
      WRITE(7,*) ISPPH,NPPPH,ISRPPH,IPPHXY
      IF(ISO.GT.0) GOTO 100
C
C73*  CONTROLS FOR HORIZONTAL PLANE VELOCITY PLOTTING
C
      NCARD=73
      DO NSKIP=1,14
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISVPH,NPVPH,ISRVPH,IVPHXY
      WRITE(7,1002)NCARD
      WRITE(7,*) ISVPH,NPVPH,ISRVPH,IVPHXY
      IF(ISO.GT.0) GOTO 100
C
C74*  CONTROLS FOR VERTICAL PLANE SCALAR FIELD CONTOURING
C
      NCARD=74
      DO NSKIP=1,14
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISECSPV,NPSPV(1),ISSPV(1),ISRSPV(1),
     &                     ISHPLTV(1)
      WRITE(7,1002)NCARD
      WRITE(7,*) ISECSPV,NPSPV(1),ISSPV(1),ISRSPV(1),
     &                     ISHPLTV(1)
      SHPLTV(1)=FLOAT(ISHPLTV(1))
      SBPLTV(1)=1.0-SHPLTV(1)
      DO N=2,7
      READ(1,*,IOSTAT=ISO) IDUMMY,NPSPV(N),ISSPV(N),ISRSPV(N),
     &                     ISHPLTV(N)
      WRITE(7,1002)NCARD
      WRITE(7,*) IDUMMY,NPSPV(N),ISSPV(N),ISRSPV(N),
     &                     ISHPLTV(N)
      SHPLTV(N)=FLOAT(ISHPLTV(N))
      SBPLTV(N)=1.0-SHPLTV(N)
      ENDDO
      IF(ISO.GT.0) GOTO 100
C
C75*  MORE CONTROLS FOR VERTICAL PLANE SCALAR FIELD CONTOURING
C
      NCARD=75
      DO NSKIP=1,8
      READ(1,1)
      ENDDO
      DO IS=1,ISECSPV
      READ(1,*,IOSTAT=ISO) DUM,NIJSPV(IS),CCTITLE(10+IS)
      WRITE(7,1002)NCARD
      WRITE(7,*) DUM,NIJSPV(IS),CCTITLE(10+IS)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C76*  I,J LOCATIONS DEFINING VERTICAL PLANE FOR CONTOURING
C
      NCARD=76
      DO NSKIP=1,8
      READ(1,1)
      ENDDO
      DO IS=1,ISECSPV
      DO NPP=1,NIJSPV(IS)
      READ(1,*,IOSTAT=ISO) DUM,ISPV(NPP,IS),JSPV(NPP,IS)
      WRITE(7,1002)NCARD
      WRITE(7,*) DUM,ISPV(NPP,IS),JSPV(NPP,IS)
      IF(ISO.GT.0) GOTO 100
      ENDDO
      ENDDO
C
C77*  CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
C
      NCARD=77
      DO NSKIP=1,11
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISECVPV,NPVPV,ISVPV,ISRVPV
      WRITE(7,1002)NCARD
      WRITE(7,*) ISECVPV,NPVPV,ISVPV,ISRVPV
      IF(ISO.GT.0) GOTO 100
C
C78*  MORE CONTROLS FOR VERTICAL PLANE VELOCITY PLOTTING
C
      NCARD=78
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      DO IS=1,ISECVPV
      READ(1,*,IOSTAT=ISO) DUM,NIJVPV(IS),ANGVPV(IS),CVTITLE(10+IS)
      WRITE(7,1002)NCARD
      WRITE(7,*) DUM,NIJVPV(IS),ANGVPV(IS),CVTITLE(10+IS)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C79*  I,J LOCATIONS FOR VERTICAL PLANE VELOCITY PLOTTING
C
      NCARD=79
      DO NSKIP=1,8
      READ(1,1)
      ENDDO
      DO IS=1,ISECVPV
      DO NPP=1,NIJVPV(IS)
      READ(1,*,IOSTAT=ISO) DUM,IVPV(NPP,IS),JVPV(NPP,IS)
      WRITE(7,1002)NCARD
      WRITE(7,*) DUM,IVPV(NPP,IS),JVPV(NPP,IS)
      IF(ISO.GT.0) GOTO 100
      ENDDO
      ENDDO
C
C80*  CONTROLS FOR 3-D GRAPHICS ARRAYS
C
      NCARD=80
      DO NSKIP=1,39
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)IS3DO,ISR3DO,NP3DO,KPC,NWGG,I3DMIN,I3DMAX,
     &                 J3DMIN,J3DMAX,I3DRW,SELVMAX,BELVMIN
      WRITE(7,1002)NCARD
      WRITE(7,*)IS3DO,ISR3DO,NP3DO,KPC,NWGG,I3DMIN,I3DMAX,
     &                 J3DMIN,J3DMAX,I3DRW,SELVMAX,BELVMIN
      IF(ISO.GT.0) GOTO 100
      NCALL3D=0
      NRCAL3D=0
C
C81*  CONTROLS FOR 3-D GRAPHICS ARRAYS
C
      NCARD=81
      DO NSKIP=1,15
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)CDUM,IS3DUUU,JS3DUUU,UUU3DMA,UUU3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DUUU,JS3DUUU,UUU3DMA,UUU3DMI
      IF(ISO.GT.0) GOTO 100
      READ(1,*,IOSTAT=ISO)CDUM,IS3DVVV,JS3DVVV,VVV3DMA,VVV3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DVVV,JS3DVVV,VVV3DMA,VVV3DMI
      IF(ISO.GT.0) GOTO 100
      READ(1,*,IOSTAT=ISO)CDUM,IS3DWWW,JS3DWWW,WWW3DMA,WWW3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DWWW,JS3DWWW,WWW3DMA,WWW3DMI
      IF(ISO.GT.0) GOTO 100
      READ(1,*,IOSTAT=ISO)CDUM,IS3DSAL,JS3DSAL,SAL3DMA,SAL3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DSAL,JS3DSAL,SAL3DMA,SAL3DMI
      IF(ISO.GT.0) GOTO 100
      READ(1,*,IOSTAT=ISO)CDUM,IS3DTEM,JS3DTEM,TEM3DMA,TEM3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DTEM,JS3DTEM,TEM3DMA,TEM3DMI
      IF(ISO.GT.0) GOTO 100
      READ(1,*,IOSTAT=ISO)CDUM,IS3DDYE,JS3DDYE,DYE3DMA,DYE3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DDYE,JS3DDYE,DYE3DMA,DYE3DMI
      IF(ISO.GT.0) GOTO 100
      READ(1,*,IOSTAT=ISO)CDUM,IS3DSED,JS3DSED,SED3DMA,SED3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DSED,JS3DSED,SED3DMA,SED3DMI
      IF(ISO.GT.0) GOTO 100
      READ(1,*,IOSTAT=ISO)CDUM,IS3DSND,JS3DSND,SND3DMA,SND3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DSND,JS3DSND,SND3DMA,SND3DMI
      IF(ISO.GT.0) GOTO 100
      READ(1,*,IOSTAT=ISO)CDUM,IS3DTOX,JS3DTOX,TOX3DMA,TOX3DMI
      WRITE(7,1002)NCARD
      WRITE(7,*)CDUM,IS3DTOX,JS3DTOX,TOX3DMA,TOX3DMI
      IF(ISO.GT.0) GOTO 100
C
C82*  READ NUMBER OF PERIODIC FORCING (TIDAL) CONSTITUENT AND HARMONIC
C     ANALYSIS PARAMETERS
C
      NCARD=82
      DO NSKIP=1,10
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) ISLSHA,MLLSHA,NTCLSHA,ISLSTR,ISHTA
      WRITE(7,1002)NCARD
      WRITE(7,*) ISLSHA,MLLSHA,NTCLSHA,ISLSTR,ISHTA
      IF(ISO.GT.0) GOTO 100
      if(mllsha.gt.mlm) then                                           !hnr
        write(*,320)mllsha,mlm                                         !hnr
320     FORMAT('The number of Harmonic Analysis Locations '            !hnr
     *  'in your application',I5,                                      !hnr
     $  '  is larger than the maximum allowed of',I5,                  !hnr
     *  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop car 82'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
C
C83*  READ HARMONIC ANALYSIS AND TIME SERIES LOCATIONS AND SWITCHES
C
      NCARD=83
      DO NSKIP=1,12
      READ(1,1)
      ENDDO
      DO M=1,MLLSHA
      READ(1,*,IOSTAT=ISO) ILLSHA(M),JLLSHA(M),LSHAP(M),LSHAB(M),
     &          LSHAUE(M),LSHAU(M),CLSL(M)
      WRITE(7,1002)NCARD
      WRITE(7,*) ILLSHA(M),JLLSHA(M),LSHAP(M),LSHAB(M),
     &          LSHAUE(M),LSHAU(M),CLSL(M)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C84*  CONTROLS FOR SAVING TIME SERIES
C
      NCARD=84
      DO NSKIP=1,17
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,
     &               TCTMSR
      WRITE(7,1002)NCARD
      WRITE(7,*)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,
     &               TCTMSR
C NEXT LINE ADDED BY DON KINERY
C    &              ,IHYDOUT,NWHYDOUT
      IF(ISO.GT.0) GOTO 100
C     MLTMSR=MIN(MLTMSR,18)
      JSTMSR=1
      NCTMSR=1
C NEXT TWO LINES ADDED BY DON KINERY
      JSHYDOUT=1
      NCHYDOUT=1
C
      if(mltmsr.gt.mltmsrm) then                                       !hnr
        write(*,330)mltmsr,mltmsrm                                      !hnr
330     FORMAT('The number of Output Time Series Locations '           !hnr
     &  'in your application',I6,                                      !hnr
     $  '  is larger than the maximum allowed of',I6,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 84_1'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      if(ntsstsp.gt.ntsstspm) then                                     !hnr
        write(*,340)ntsstsp,ntsstspm                                   !hnr
340     FORMAT('The number of Time Series start/stop scenarios '       !hnr
     &  'in your application',I3,                                      !hnr
     $  '  is larger than the maximum allowed of',I3,                  !hnr
     *  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop card 84_2'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
C85*  CONTROLS FOR SAVING TIME SERIES
C
      NCARD=85
      DO NSKIP=1,7
      READ(1,1)
      ENDDO
      DO ITSSS=1,NTSSTSP
       READ(1,*,IOSTAT=ISO)IDUM,MTSSTSP(ITSSS)
       WRITE(7,1002)NCARD
      WRITE(7,*)IDUM,MTSSTSP(ITSSS)
       IF(ISO.GT.0) GOTO 100
      ENDDO
      if(mtsstsp(itsss).gt.mtsstspm) then                              !hnr
        write(*,350)mtsstsp(itsss),itsss,mtsstspm                      !hnr
350     FORMAT('The number of stop/start pairs Locations '             !hnr
     &  'in your application',I6,                                      !hnr
     $  '  for scenario',I6,'  is larger than the maximum allowed of', !hnr
     *  I6,'  for this compilation')                                   !hnr
        write(*,*)'The program will stop card 85'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
C
C86*  CONTROLS FOR SAVING TIME SERIES
C
      NCARD=86
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      DO ITSSS=1,NTSSTSP
       DO MTSSS=1,MTSSTSP(ITSSS)
       READ(1,*,IOSTAT=ISO)IDUM,IDUM,TSSTRT(MTSSS,ITSSS),
     &                               TSSTOP(MTSSS,ITSSS)
       WRITE(7,1002)NCARD
      WRITE(7,*)IDUM,IDUM,TSSTRT(MTSSS,ITSSS),
     &                               TSSTOP(MTSSS,ITSSS)
       IF(ISO.GT.0) GOTO 100
       ENDDO
      ENDDO
C
C87*  CONTROLS FOR SAVING TIME SERIES
C
      NCARD=87
      DO NSKIP=1,17
      READ(1,1)
      ENDDO
      DO M=1,MLTMSR
      READ(1,*,IOSTAT=ISO)ILTMSR(M),JLTMSR(M),NTSSSS(M),MTMSRP(M),
     &       MTMSRC(M),MTMSRA(M),MTMSRUE(M),MTMSRUT(M),MTMSRU(M),
     &                 MTMSRQE(M),MTMSRQ(M),CLTMSR(M)
      WRITE(7,1002)NCARD
      WRITE(7,*)ILTMSR(M),JLTMSR(M),NTSSSS(M),MTMSRP(M),
     &       MTMSRC(M),MTMSRA(M),MTMSRUE(M),MTMSRUT(M),MTMSRU(M),
     &                 MTMSRQE(M),MTMSRQ(M),CLTMSR(M)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C88*  CONTROLS FOR EXTRACTING INSTANTANEOUS VERTICAL SCALAR FIELD
C     PROFILES
C
      NCARD=88
      DO NSKIP=1,10
      READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO)ISVSFP,MDVSFP,MLVSFP,TMVSFP,TAVSFP
      WRITE(7,1002)NCARD
      WRITE(7,*)ISVSFP,MDVSFP,MLVSFP,TMVSFP,TAVSFP
      IF(ISO.GT.0) GOTO 100
      JSVSFP=1
C
C89*  SAMPLING DEPTHS FOR EXTRACTING INST VERTICAL SCALAR FIELD
C     PROFILES
C
      NCARD=89
      DO NSKIP=1,7
      READ(1,1)
      ENDDO
      DO M=1,MDVSFP
      READ(1,*,IOSTAT=ISO)IDUM,DMVSFP(M)
      WRITE(7,1002)NCARD
      WRITE(7,*)IDUM,DMVSFP(M)
      IF(ISO.GT.0) GOTO 100
      ENDDO
C
C90*  HORIZONTAL SPACE-TIME LOCATIONS FOR SAMPLING
C
      NCARD=90
      DO NSKIP=1,9
      READ(1,1)
      ENDDO
      DO M=1,MLVSFP
      READ(1,*,IOSTAT=ISO)IDUM,TIMVSFP(M),IVSFP(M),JVSFP(M)
      WRITE(7,1002)NCARD
      WRITE(7,*)IDUM,TIMVSFP(M),IVSFP(M),JVSFP(M)
      IF(ISO.GT.0) GOTO 100
      ENDDO
      NTS=NTC*NTSPTC
      NBVSFP=NTC*NTSPTC
      NSVSFP=0
      DO M=1,MLVSFP
      TIMVSFP(M)=TMVSFP*(TIMVSFP(M)+TAVSFP)
      ENDDO
      DT=TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)
      DO M=1,MLVSFP
      NTMP=NINT( (TIMVSFP(M)-TCON*TBEGIN)/DT )
      NTMP=MIN(NTMP,NTS)
      NBVSFP=MIN(NBVSFP,NTMP)-1
      NSVSFP=MAX(NSVSFP,NTMP)+1
      NTVSFP(M)=NTMP
      ENDDO
      DO M=1,MLVSFP
      TIMVSFP(M)=(TIMVSFP(M)/TMVSFP)-TAVSFP
      ENDDO
C
C**********************************************************************C
C
      IF(ISHOUSATONIC.EQ.1)THEN
C
C91
C
C#############################################################################
C     HQI change to input settling velocity function parameters
C     RM 06/24/04
C   Reference settling velocity and concentration above which concentration
C   enhanced occurs and below which constant settling velocity
C
      DO NSKIP=1,10
         READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) SED_CRIT, CONSTWS1, CONSTWS2
      IF(ISO.GT.0) THEN
         WRITE(*,*)'Reference settling velocity not input'
         STOP
      ENDIF
C#############################################################################
C
C92
C
C#############################################################################
C     HQI change to input parameters to maintain a constant TSS in the reach
C     5D backwaters
C     RM 07/30/04

      DO NSKIP=1,8
         READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) LABKWTR, SED_MIN
      IF(ISO.GT.0) THEN
         WRITE(*,*)'Min TSS in backwaters not input'
         STOP
      ENDIF
C
C     READ CELL INDEXES OF CELLS IN BACKWATERS
C
      IF(LABKWTR.GT.0)THEN
      OPEN(101,FILE='BKWTR.INP',STATUS='UNKNOWN')
      READ(101,*)
      DO IS=1,LABKWTR
         READ(101,*)IJBKWTR(IS)
      ENDDO
      CLOSE(101)
      ENDIF
C
C#############################################################################
C
C93
C
C#############################################################################
C     HQI change to input spatially varying cohesive erosion function parameters
C     RM 10/03/05

      DO NSKIP=1,10
         READ(1,1)
      ENDDO
      READ(1,*,IOSTAT=ISO) COEFF_US, EXPO_US, COEFF_DS, EXPO_DS
      IF(ISO.GT.0) THEN
         WRITE(*,*)'Cohesive resuspension parameters not input'
         STOP
      ENDIF
C#############################################################################
C
C94
C
C#############################################################################
C     HQI change to input spatially foodchain averaging depths
C     RM 10/05/05

      DO NSKIP=1,9
         READ(1,1)
      ENDDO
C
      IF(ISFDCH.GT.0)THEN
      READ(1,*,IOSTAT=ISO) FCMHBDUS, FCMHBDDS, FCMHBDWP
      IF(ISO.GT.0) THEN
         WRITE(*,*)'Foodchain exposure depths not input'
         STOP
      ENDIF
      ENDIF
C
C      DO NSKIP=1,5
C         READ(1,1)
C      ENDDO
C      READ(1,*,IOSTAT=ISO)ISTOXALL,NSTOXALL
C      JSTOXALL=0
C      MSTOXALL=NTSPTC/NSTOXALL
C#############################################################################
C
      ENDIF
C
C**********************************************************************C
C
      GOTO 2000
C
C**********************************************************************C
C
C **  WRITE INPUT ERROR MESSAGES AND TERMINATE RUN
C

  100 WRITE(*,1001)NCARD
!      WRITE(8,1001)NCARD   !hnr 7/27/2009
      WRITE(7,1001)NCARD
      STOP
 2000 CONTINUE
C
C**********************************************************************C
C
C **  NOW REWIND UNIT 1 & READ IN AS CHARACTER TO WRITE TO UNIT 7
C
C----------------------------------------------------------------------C
C
      REWIND (1)
   21 READ(1,22,END=24) TEXT
      WRITE (7,23) TEXT
      GOTO 21
   24 CONTINUE
      CLOSE(1)
   22 FORMAT (A80)
   23 FORMAT (1X,A80)
C
C**********************************************************************C
C
      IF(ISDRY.GT.0)THEN
C
      OPEN(1,FILE='WETDRYCHG.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WETDRYCHG.OUT',POSITION='APPEND')
      TIME=TCON*TBEGIN/86400.
      WRITE(1,1234)TIME,NWETCELLSNEW,NWETCELLSOLD
C
	ENDIF
C
 1234 FORMAT(F12.4,2I8)
C
C**********************************************************************C
C
C **  READ CELL TYPES FROM FILES CELL.INP AND CELLLT.INP
C
C----------------------------------------------------------------------C
C
      OPEN(1,FILE='CELL.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES AND DETERMINE FILE FORMAT
C
      DO IS=1,4
      READ(1,1)
      ENDDO
      READ(1,*)JCTMP
      CLOSE(1)
      OPEN(1,FILE='CELL.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,4
      READ(1,1)
      ENDDO
C
      IF(JCTMP.NE.JC)THEN
C **    READ OLD FILE FORMAT
        JACROSS=JC
        IF(JC.GT.640)JACROSS=640
        DO JT=1,JC,JACROSS
        JF=JT
        JLAST=JT+JACROSS-1
        IF(JLAST.GT.JC) JLAST=JC
        WRITE (7,8)JF,JLAST
        DO I=1,IC
        READ(1,6,IOSTAT=ISO) (IJCT(I,J),J=JF,JLAST)
        IF(ISO.GT.0) GOTO 800
        WRITE (7,16) (IJCT(I,J),J=JF,JLAST)
        ENDDO
        WRITE(7,15)
        ENDDO
       ELSE
C **    READ NEW FILE FORMAT
C       READ(1,66)JCTMP
        IF(IC.GT.640)THEN
          IACROSS=640
          DO IT=1,IC,IACROSS
           IFIRST=IT
           ILAST=IT+IACROSS-1
           IF(ILAST.GT.IC) ILAST=IC
           WRITE (7,88)IFIRST,ILAST
           IF(JC.LE.999)THEN
             DO J=JC,1,-1
             READ(1,66,IOSTAT=ISO)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)
             IF(ISO.GT.0) GOTO 800
             WRITE (7,166)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)
             ENDDO
           ELSE
             DO J=JC,1,-1
             READ(1,660,IOSTAT=ISO)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)
             IF(ISO.GT.0) GOTO 800
             WRITE (7,1660)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)
             ENDDO
           ENDIF
           WRITE(7,15)
          ENDDO
         ELSE
          IFIRST=1
          ILAST=IC
          WRITE (7,88)IFIRST,ILAST
           IF(JC.LE.999)THEN
             DO J=JC,1,-1
             READ(1,66,IOSTAT=ISO)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)
             IF(ISO.GT.0) GOTO 800
             WRITE (7,166)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)
             ENDDO
           ELSE
             DO J=JC,1,-1
             READ(1,660,IOSTAT=ISO)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)
             IF(ISO.GT.0) GOTO 800
             WRITE (7,1660)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)
             ENDDO
           ENDIF
          WRITE(7,15)
        ENDIF
      ENDIF
C
      CLOSE(1)
C
C **  IF JCTMP NE JC WRITE NEWCELL.INP FILE
C
      IF(JCTMP.NE.JC)THEN
        OPEN(1,FILE='NEWCELL.INP',STATUS='UNKNOWN')
        IACROSS=IC
        IF(IC.GT.640)IACROSS=640
        DO IT=1,IC,IACROSS
        IFIRST=IT
        ILAST=IT+IACROSS-1
        IF(ILAST.GT.IC) ILAST=IC
        DO J=JC,1,-1
        WRITE(1,66)J,(IJCT(I,J),I=IFIRST,ILAST)
        ENDDO
        ENDDO
        CLOSE(1)
      ENDIF
C
C------------------------------------------------------------------------C
C
      IF(ISWASP.GT.0)THEN
      OPEN(1,FILE='CELLLT.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES AND DETERMINE FILE FORMAT
C
      DO IS=1,4
      READ(1,1)
      ENDDO
!      READ(1,660)JCTMP        !hnr
      READ(1,*)JCTMP           !hnr
      CLOSE(1)
      OPEN(1,FILE='CELLLT.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,4
      READ(1,1)
      ENDDO
C
      IF(JCTMP.NE.JC)THEN
C **    READ OLD FILE FORMAT
        JACROSS=JC
        IF(JC.GT.640)JACROSS=640
        DO JT=1,JC,JACROSS
        JF=JT
        JLAST=JT+JACROSS-1
        IF(JLAST.GT.JC) JLAST=JC
        WRITE (7,8)JF,JLAST
        DO I=1,IC
        READ(1,6,IOSTAT=ISO) (IJCTLT(I,J),J=JF,JLAST)
        IF(ISO.GT.0) GOTO 800
        WRITE (7,16) (IJCTLT(I,J),J=JF,JLAST)
        ENDDO
        WRITE(7,15)
        ENDDO
       ELSE
C **    READ NEW FILE FORMAT
C       READ(1,66)JCTMP
        IF(IC.GT.640)THEN
          IACROSS=640
          DO IT=1,IC,IACROSS
           IFIRST=IT
           ILAST=IT+IACROSS-1
           IF(ILAST.GT.IC) ILAST=IC
           WRITE (7,88)IFIRST,ILAST
           IF(JC.LE.999)THEN
             DO J=JC,1,-1
!             READ(1,66,IOSTAT=ISO)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)               !hnr
             READ(1,66,IOSTAT=ISO)JDUMY,(IJCTlt(I,J),I=IFIRST,ILAST)              !hnr
             IF(ISO.GT.0) GOTO 800
!             WRITE (7,166)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)                        !hnr
             WRITE (7,166)JDUMY,(IJCT(ltI,J),I=IFIRST,ILAST)                         !hnr
             ENDDO
           ELSE
             DO J=JC,1,-1
!             READ(1,660,IOSTAT=ISO)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)              !hnr
             READ(1,660,IOSTAT=ISO)JDUMY,(IJCTlt(I,J),I=IFIRST,ILAST)               !hnr
             IF(ISO.GT.0) GOTO 800
!             WRITE (7,1660)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)                      !hnr
             WRITE (7,1660)JDUMY,(IJCT(ltI,J),I=IFIRST,ILAST)                       !hnr
             ENDDO
           ENDIF
           WRITE(7,15)
          ENDDO
         ELSE
          IFIRST=1
          ILAST=IC
          WRITE (7,88)IFIRST,ILAST
           IF(JC.LE.999)THEN
             DO J=JC,1,-1
!             READ(1,66,IOSTAT=ISO)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)             !hnr
             READ(1,66,IOSTAT=ISO)JDUMY,(IJCTlt(I,J),I=IFIRST,ILAST)              !hnr
             IF(ISO.GT.0) GOTO 800
!             WRITE (7,166)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)              !hnr
             WRITE (7,166)JDUMY,(IJCTlt(I,J),I=IFIRST,ILAST)              !hnr
             ENDDO
           ELSE
             DO J=JC,1,-1
!             READ(1,660,IOSTAT=ISO)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)              !hnr
             READ(1,660,IOSTAT=ISO)JDUMY,(IJCTlt(I,J),I=IFIRST,ILAST)              !hnr
             IF(ISO.GT.0) GOTO 800
!             WRITE (7,1660)JDUMY,(IJCT(I,J),I=IFIRST,ILAST)              !hnr
             WRITE (7,1660)JDUMY,(IJCTlt(I,J),I=IFIRST,ILAST)              !hnr
             ENDDO
           ENDIF
          WRITE(7,15)
        ENDIF
      ENDIF
C
      CLOSE(1)
	ENDIF
C
      IF(ISLTMT.GE.1)THEN
      DO J=1,JC
      DO I=1,IC
      IJCT(I,J)=IJCTLT(I,J)
      ENDDO
      ENDDO
      ENDIF
C
    8 FORMAT ('   CELL TYPE ARRAY,J=',I5,2X,'TO J=',I5,//)
   88 FORMAT ('   CELLLT TYPE ARRAY,I=',I5,2X,'TO I=',I5,//)
C
C**********************************************************************C
C
C **  IF ISPGNS GE 1, READ IN NORTH-SOUTH BOUNDARY CELLS FROM
C **  FILE MAPPGNS.INP TO SPECIFY A PERIODIC DOMAIN IN THE NORTH-SOUTH
C **  DIRECTION
C
      IF(ISPGNS.GE.1)THEN
C
      OPEN(1,FILE='MAPPGNS.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,8
      READ(1,1)
      ENDDO
C
      READ(1,*,IOSTAT=ISO) NPNSBP
      IF(ISO.GT.0) GOTO 950
C
      DO NPNS=1,NPNSBP
      READ(1,*,IOSTAT=ISO) ISPNS(NPNS),JSPNS(NPNS),
     &                     INPNS(NPNS),JNPNS(NPNS)
      IF(ISO.GT.0) GOTO 950
      ENDDO
C
      CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  GENERATE CELL MAPPINGS
C
      CALL CELLMAP
C
C**********************************************************************C
C
C **  READ IN CELL CENTER DEPTHS AND BOTTOM BED ELEVATION FOR CARTESIAN
C **  OR MIXED CATERSIAN CURVILINEAR GRID FORM FILE DEPTH.INP.  IF THE
C **  NUMBER OF VARIALBLE CELLS EQUALS NUMBER OF WATER CELLS OR
C **  ISCLO=1 DEPTHS ARE READ FROM THE FILE DXDY.INP
C
C----------------------------------------------------------------------C
C
C     LCM2=LC-2
C     IF(LCM2.GT.LVC)THEN
C
C     OPEN(1,FILE='DEPTH.INP',STATUS='UNKNOWN')
C     GOTO 999
C
C     WRITE (7,9)
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
C     DO IS=1,5
C     READ(1,1)
C     ENDDO
C
C     IF(JC.LE.30)THEN
C       DO I=1,IC
C       READ(1,7,IOSTAT=ISO) (RTMP1IJ(I,J),J=1,JC)
C       IF(ISO.GT.0) GOTO 820
C       WRITE (7,17) (RTMP1IJ(I,J),J=1,JC)
C       ENDDO
C      ELSE
C       DO JJ=1,JC,30
C       JF=JJ
C       JLAST=JF+29
C       IF(JLAST.GT.JC) JLAST=JC
C       DO I=1,IC
C       READ(1,7,IOSTAT=ISO) (RTMP1IJ(I,J),J=JF,JLAST)
C       IF(ISO.GT.0) GOTO 820
C       WRITE (7,17) (RTMP1IJ(I,J),J=JF,JLAST)
C       ENDDO
C       WRITE (7,15)
C       ENDDO
C     ENDIF
C
C     DO L=2,LA
C     HMP(L)=RTMP1IJ(IL(L),JL(L)
C     ENDDO
C
C **  SKIP OVER BLANK, TITLE AND AND HEADER LINES
C
C     DO IS=1,5
C     READ(1,1)
C     ENDDO
C
C     IF(JC.LE.30)THEN
C       DO I=1,IC
C       READ(1,7,IOSTAT=ISO) (RTMP1IJ(I,J),J=1,JC)
C       IF(ISO.GT.0) GOTO 820
C       WRITE (7,17) (RTMP1IJ(I,J),J=1,JC)
C       ENDDO
C      ELSE
C       DO JJ=1,JC,30
C       JF=JJ
C       JLAST=JF+29
C       IF(JLAST.GT.JC) JLAST=JC
C       DO I=1,IC
C       READ(1,7,IOSTAT=ISO) (RTMP1IJ(I,J),J=JF,JLAST)
C       IF(ISO.GT.0) GOTO 820
C       WRITE (7,17) (RTMP1IJ(I,J),J=JF,JLAST)
C       ENDDO
C       WRITE (7,15)
C       ENDDO
C     ENDIF
C
C     DO L=2,LA
C     BELV(L)=RTMP1IJ(IL(L),JL(L))
C     ENDDO
C
C **  SKIP OVER BLANK, TITLE AND AND HEADER LINES
C
C     DO IS=1,5
C     READ(1,1)
C     ENDDO
C
C     IF(JC.LE.30)THEN
C       DO I=1,IC
C       READ(1,7,IOSTAT=ISO) (RTMP1IJ(I,J),J=1,JC)
C       IF(ISO.GT.0) GOTO 820
C       WRITE (7,17) (RTMP1IJ(I,J),J=1,JC)
C       ENDDO
C      ELSE
C       DO JJ=1,JC,30
C       JF=JJ
C       JLAST=JF+29
C       IF(JLAST.GT.JC) JLAST=JC
C       DO I=1,IC
C       READ(1,7,IOSTAT=ISO) (RTMP1IJ(I,J),J=JF,JLAST)
C       IF(ISO.GT.0) GOTO 820
C       WRITE (7,17) (RTMP1IJ(I,J),J=JF,JLAST)
C       ENDDO
C       WRITE (7,15)
C       ENDDO
C     ENDIF
C
C     DO L=2,LA
C     ZBR(L)=RTMP1IJ(IL(L),JL(L))
C     ENDDO
C
C     CLOSE(1)
C
C     ENDIF
C
   15 FORMAT (/)
C   6 FORMAT (60I2)
    6 FORMAT (640I1)
   66 FORMAT (I3,2X,640I1)
  660 FORMAT (I4,2X,640I1)
    9 FORMAT (/,' DEPTH ARRAY:',//)
   16 FORMAT (1X,640I1)
  166 FORMAT (1X,I3,2X,640I1)
 1660 FORMAT (1X,I4,2X,640I1)
    7 FORMAT (30F4.1)
   17 FORMAT(1X,30F4.1)
C
C**********************************************************************C
C
C **  READ CURVILINEAR-ORTHOGONAL OR VARIABLE CELL DATA FROM FILE
C **  DXDY.INP
C
C----------------------------------------------------------------------C
C
C **  INITIALIZE CELL DIMENSIONS TO CONSTANT CARTESIAN OR DUMMY VALUES
C
      DO L=1,LC
      DXP(L)=DX*DXYCVT
      DYP(L)=DY*DXYCVT
      ZBR(L)=ZBRADJ
      ENDDO
C
C **  READ IN DX, DY, DEPTH AND BOTTOM ELEVATION AT CELL CENTERS OF
C **  VARIABLE CELLS
C
      IF(LVC.GT.0)THEN
      OPEN(1,FILE='DXDY.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,4
      READ(1,1)
      ENDDO
C
      IF(ISVEG.EQ.0)THEN
        DO LT=1,LVC
C       READ(1,*,IOSTAT=ISO)I,J,DXIJ,DYIJ,HIJ,BELVIJ
        READ(1,*,IOSTAT=ISO)I,J,DXIJ,DYIJ,HIJ,BELVIJ,ZBRIJ
        IF(ISO.GT.0) GOTO 830
        L=LIJ(I,J)
        DXP(L)=DXYCVT*DXIJ
        DYP(L)=DXYCVT*DYIJ
        HMP(L)=HADJ + HCVRT*HIJ
        HMP(L)=MAX(HMP(L),HMIN)
        BELV(L)=BELADJ + BELCVRT*BELVIJ
        BELV1(L)=BELADJ + BELCVRT*BELVIJ
        ZBR(L)=ZBRADJ + ZBRCVRT*ZBRIJ
        ENDDO
       ELSE
        DO LT=1,LVC
C       READ(1,*,IOSTAT=ISO)I,J,DXIJ,DYIJ,HIJ,BELVIJ
        READ(1,*,IOSTAT=ISO)I,J,DXIJ,DYIJ,HIJ,BELVIJ,ZBRIJ,MVEGIJT
        IF(ISO.GT.0) GOTO 830
        L=LIJ(I,J)
        DXP(L)=DXYCVT*DXIJ
        DYP(L)=DXYCVT*DYIJ
        HMP(L)=HADJ + HCVRT*HIJ
        HMP(L)=MAX(HMP(L),HMIN)
        BELV(L)=BELADJ + BELCVRT*BELVIJ
        BELV1(L)=BELADJ + BELCVRT*BELVIJ
        ZBR(L)=ZBRADJ + ZBRCVRT*ZBRIJ
        MVEGL(L)=MVEGIJT
        ENDDO
      ENDIF
C
C ** MODIFY DX AND DY FOR 1D CHANNELS
C
      IF(IS1DCHAN.GT.0)THEN
       DO L=2,LA
         IF(LCT(L).EQ.6) DYP(L)=1.
         IF(LCT(L).EQ.7) DXP(L)=1.
       ENDDO
      ENDIF
C
      CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE MODDXDY.INP TO MODIFY INPUT VALUES OF DX AND DY
C
      IF(IMDXDY.GT.0)THEN
      OPEN(1,FILE='MODDXDY.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,4
      READ(1,1)
      ENDDO
C
      READ(1,*) NMDXDY
      IF(NMDXDY.GE.1)THEN
        DO NMD=1,NMDXDY
        READ(1,*)ITMP,JTMP,RMDX,RMDY
        LTMP=LIJ(ITMP,JTMP)
        DXP(LTMP)=RMDX*DXP(LTMP)
        DYP(LTMP)=RMDY*DYP(LTMP)
        ENDDO
      ENDIF
C
      CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE MODCHAN.INP TO INSERT SUBGRID CHANNELS INTO
C **  HOST CELLS
C
      MDCHH=0
      IF(ISCHAN.GT.0)THEN
      OPEN(1,FILE='MODCHAN.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,8
      READ(1,1)
      ENDDO
C
      IF(ISCHAN.EQ.1)THEN
      READ(1,*) MDCHH,MDCHHD,MDCHHD2
      READ(1,*) MDCHITM,MDCHHQ,QCHERR
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        READ(1,*)MDCHTYP(NMD),IMDCHH(NMD),JMDCHH(NMD),
     &                        IMDCHU(NMD),JMDCHU(NMD),
     &                        IMDCHV(NMD),JMDCHV(NMD)
        QCHANU(NMD)=0.
        QCHANUN(NMD)=0.
        QCHANV(NMD)=0.
        QCHANVN(NMD)=0.
        ENDDO
      ENDIF
      ENDIF
C
      IF(ISCHAN.EQ.2)THEN
      READ(1,*) MDCHH,MDCHHD,MDCHHD2
      READ(1,*) MDCHITM,MDCHHQ,QCHERR
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        READ(1,*)MDCHTYP(NMD),IMDCHH(NMD),JMDCHH(NMD),
     &                        IMDCHU(NMD),JMDCHU(NMD),
     &                        IMDCHV(NMD),JMDCHV(NMD),
     &                        CHANLEN(NMD),PMDCH(NMD)
        QCHANU(NMD)=0.
        QCHANUN(NMD)=0.
        QCHANV(NMD)=0.
        QCHANVN(NMD)=0.
        ENDDO
      ENDIF
      ENDIF
C
      CLOSE(1)
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        LMDCHH(NMD)=LIJ(IMDCHH(NMD),JMDCHH(NMD))
        IF(IMDCHU(NMD).EQ.1.AND.JMDCHU(NMD).EQ.1)THEN
          LMDCHU(NMD)=1
         ELSE
          LMDCHU(NMD)=LIJ(IMDCHU(NMD),JMDCHU(NMD))
        ENDIF
        IF(IMDCHV(NMD).EQ.1.AND.JMDCHV(NMD).EQ.1)THEN
          LMDCHV(NMD)=1
         ELSE
          LMDCHV(NMD)=LIJ(IMDCHV(NMD),JMDCHV(NMD))
        ENDIF
        ENDDO
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE CHANSEC.INP FOR 1-D CHANNEL CROSS SECTION DATA
C
      IF(IS1DCHAN.GT.0)THEN
      OPEN(1,FILE='CHANSEC.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,21
      READ(1,1)
      ENDDO
C
      DO LL=1,LC-2
C
      READ(1,*) ISEC,JSEC,ISECDAT,NXYSDATT,BELVBT,RMULADJ
      L=LIJ(ISEC,JSEC)
      NXYSDAT(L)=NXYSDATT
      BELVB(L)=BELVBT
      IF(LCT(L).EQ.6)THEN
        DO ND=1,NXYSDAT(L)
         READ(1,*)EHXYS(ND,L),AREADY(ND,L),WPERDY(ND,L),SURFDY(ND,L)
         AREADX(ND,L)=1.
         WPERDX(ND,L)=1.
         SURFDX(ND,L)=0.
        ENDDO
      ENDIF
      IF(LCT(L).EQ.7)THEN
        DO ND=1,NXYSDAT(L)
         READ(1,*)EHXYS(ND,L),AREADX(ND,L),WPERDX(ND,L),SURFDX(ND,L)
         AREADY(ND,L)=1.
         WPERDY(ND,L)=1.
         SURFDY(ND,L)=0.
        ENDDO
      ENDIF
      BELVB(L)=RMULADJ*BELVB(L)
      DO ND=1,NXYSDAT(L)
       EHXYS(ND,L)=RMULADJ*EHXYS(ND,L)
       AREADX(ND,L)=RMULADJ*RMULADJ*AREADX(ND,L)
       WPERDX(ND,L)=RMULADJ*WPERDX(ND,L)
       SURFDX(ND,L)=RMULADJ*SURFDX(ND,L)
       AREADY(ND,L)=RMULADJ*RMULADJ*AREADY(ND,L)
       WPERDY(ND,L)=RMULADJ*WPERDY(ND,L)
       SURFDY(ND,L)=RMULADJ*SURFDY(ND,L)
      ENDDO
      IF(ISECDAT.EQ.0)THEN
        DO ND=1,NXYSDAT(L)
         EHXYS(ND,L)=EHXYS(ND,L)-BELVB(L)
        ENDDO
      ENDIF
C
      ENDDO
C
      CLOSE(1)
      ENDIF
C
C
C**********************************************************************C
C
C **  OPEN FILE CHANJUN.INP FOR 1-D CHANNEL JUNCTION DATA
C
      IF(IS1DCHAN.GT.0)THEN
C
      OPEN(1,FILE='CHANJUN.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,17
      READ(1,1)
      ENDDO
C
      NJUNX=0
      NJUNY=0
      READ(1,*) NJUNTMP
CDIA      WRITE(6,6543) NJUNTMP
C
      IF(NJUNTMP.GT.0)THEN
C
      DO NJ=1,NJUNTMP
       READ(1,*) IJUNTMP,JJUNTMP,JDIR,JUNTYP
       IF(JDIR.EQ.1)THEN
         NJUNX=NJUNX+1
         LJUNX(NJUNX)=LIJ(IJUNTMP,JJUNTMP)
         JUNTPX(NJUNX)=JUNTYP
       ENDIF
       IF(JDIR.EQ.2)THEN
         NJUNY=NJUNY+1
         LJUNY(NJUNY)=LIJ(IJUNTMP,JJUNTMP)
         JUNTPY(NJUNY)=JUNTYP
       ENDIF
      ENDDO
C
      ENDIF
C
      CLOSE(1)
C
      ENDIF
C
 6543 FORMAT(' NJUNTMP = ',I5)
C
C**********************************************************************C
C
C **  OPEN FILE GWATER.INP TO SPECIFY GROUNDWATER INTERACTION
C **  BY INFILTRATION AND EVAPOTRANSPIRATION
C
      ISGWIE=0
      IF(ISGWIT.EQ.1)THEN
      OPEN(1,FILE='GWATER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,6
      READ(1,1)
      ENDDO
C
      READ(1,*) ISGWIE
C     WRITE(6,339)ISGWIE
      IF(ISGWIE.GE.1)THEN
        READ(1,*) DAGWZ,RNPOR,RIFTRM
C       WRITE(6,907) DAGWZ,RNPOR,RIFTRM
       ELSE
        DAGWZ=0.0
        RNPOR=1.E-12
        RIFTRM=0.0
      ENDIF
C
      CLOSE(1)
      ENDIF
C
  339 FORMAT(2I5,6F14.5)
C
C**********************************************************************C
C
C **  OPEN FILE FBODY.INP TO READ IN SPATIALLY VARYING BODY FORCES
C
      IF(ISBODYF.GE.1)THEN
        OPEN(1,FILE='FBODY.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,7
          READ(1,1)
        ENDDO
C
        READ(1,*)CVTFACX,CVTFACY
        DO LL=2,LA
          READ(1,*)ITMP,JTMP,FBODY1,FBODY2
          L=LIJ(ITMP,JTMP)
          IF(ISBODYF.EQ.1)THEN
	      DO K=1,KC
              FBODYFX(L,K)=CVTFACX*FBODY1
              FBODYFY(L,K)=CVTFACY*FBODY2
            ENDDO
          ENDIF
          IF(ISBODYF.EQ.2)THEN
	      DO K=1,KC-1
              FBODYFX(L,K)=0.0
              FBODYFY(L,K)=0.0
            ENDDO
            FBODYFX(L,KC)=CVTFACX*FBODY1
            FBODYFY(L,KC)=CVTFACY*FBODY2
          ENDIF
        END DO
	  DO K=1,KC
          FBODYFX(1,K)=0.0
          FBODYFY(1,K)=0.0
          FBODYFX(LC,K)=0.0
          FBODYFY(LC,K)=0.0
        END DO
        DO L=1,LC
          FBODYFXI(L)=0.0
          FBODYFYI(L)=0.0
        END DO
C
        CLOSE(1)
      ELSE
	  DO K=1,KC
          DO L=1,LC
            FBODYFX(L,K)=0.0
            FBODYFY(L,K)=0.0
          END DO
        END DO
        DO L=1,LC
          FBODYFXI(L)=0.0
          FBODYFYI(L)=0.0
        END DO
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE SEDBLBC.INP TO READ IN SEDIMENT BEDLOAD OUTFLOW
C **  OR RECIRCULATION BOUNDARY CONDITIONS
C
      NSBDLDBC=0
C
      IF(ISBDLDBC.GE.1)THEN
        OPEN(1,FILE='SEDBLBC.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,12
          READ(1,1)
        ENDDO
C
        READ(1,*)NSBDLDBC
        DO N=1,NSBDLDBC
          READ(1,*)ITMPU,JTMPU,ITMPD,JTMPD,ISDBLDIR(N)
          LSBLBCU(N)=LIJ(ITMPU,JTMPU)
          IF(ITMPD.GT.0.AND.JTMPD.GT.0) THEN
            LSBLBCD(N)=LIJ(ITMPD,JTMPD)
          ELSE
            LSBLBCD(N)=0
          ENDIF
        ENDDO
C
        CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE GWMAP.INP TO SPECIFY GROUNDWATER INTERACTION BY
C **  AMBIENT GROUNDWATER FLOW
C
      IF(ISGWIT.EQ.2)THEN
        OPEN(1,FILE='GWMAP.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,10
          READ(1,1)
        ENDDO
C
        DO LL=2,LA
          READ(1,*)ITMP,JTMP,IVALUE,RVALUE
          L=LIJ(ITMP,JTMP)
          NGWSL(L)=IVALUE
          GWFAC(L)=RVALUE
        END DO
C
        CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  READ IN SPATIALLY VARYING SEDIMENT ROUGHNESS HEIGHT FOR
C **  DETERMINING GRAIN STRESS
C
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
      IF(ISBEDSTR.EQ.3)THEN
C
         OPEN(1,FILE='SEDROUGH.INP')
	   DO IS=1,2
            READ(1,1)
         ENDDO
         DO L=2,LC-1
           READ(1,*) LDUM,IDUM,JDUM,ZBRSED(L)
         ENDDO
	   CLOSE(1)
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE DOCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
C **  DISSOLVED ORGANIC CARBON IN WATER COLUMN
C
      IVAL=0
	DO NT=1,NTOX
	  IF(ISTOC(NT).EQ.1.OR.ISTOC(NT).EQ.2)IVAL=1
	ENDDO
C
      IF(IVAL.EQ.1)THEN
C
      IF(ISTDOCW.EQ.1)THEN
        OPEN(1,FILE='DOCW.INP',STATUS='UNKNOWN')
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (STDOCW(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 854
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STDOCW(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 854
          ENDDO
        ENDIF
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE POCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
C **  PARTICULATE ORGANIC CARBON IN WATER COLUMN
C
      IVAL=0
	DO NT=1,NTOX
	  IF(ISTOC(NT).EQ.1)IVAL=1
	ENDDO
C
      IF(IVAL.EQ.1)THEN
C
      IF(ISTPOCW.EQ.1)THEN
        OPEN(1,FILE='POCW.INP',STATUS='UNKNOWN')
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (STPOCW(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 854
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STPOCW(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 854
          ENDDO
        ENDIF
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE FPOCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
C **  PARTICULATE ORGANIC CARBON FRACTION FOR EACH SEDIMENT CLASS
C **  IN WATER COLUMN
C
      IVAL=0
	DO NT=1,NTOX
	  IF(ISTOC(NT).EQ.2.OR.ISTOC(NT).EQ.3)IVAL=1
	ENDDO
C
      IF(IVAL.EQ.1)THEN
C
      IF(ISTPOCW.EQ.3)THEN
        OPEN(1,FILE='FPOCW.INP',STATUS='UNKNOWN')
	  DO NS=1,NSED+NSND
          DO IS=1,8
            READ(1,1)
          ENDDO
          READ(1,*)ISALTYP
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STFPOCW(L,K,NS),K=1,KC)
              IF(ISO.GT.0) GOTO 854
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                            (STFPOCW(L,K,NS),K=1,KC)
              IF(ISO.GT.0) GOTO 854
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE DOCB.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
C **  DISSOLVED ORGANIC CARBON IN SEDIMENT BED
C
      IVAL=0
	DO NT=1,NTOX
	  IF(ISTOC(NT).EQ.1.OR.ISTOC(NT).EQ.2)IVAL=1
	ENDDO
C
      IF(IVAL.EQ.1)THEN
C
      IF(ISTDOCB.EQ.1)THEN
        OPEN(1,FILE='DOCB.INP',STATUS='UNKNOWN')
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)ISALTYP,IREAD,KBINPUT
        DO K=1,KB
	    DO L=2,LC-1
	      STDOCB(L,K)=0.0
	    ENDDO
	  ENDDO
	  IF(IREAD.EQ.0)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) STDOCB(L,1)
              IF(ISO.GT.0) GOTO 856
	        DO K=2,KB
	          STDOCB(L,K)=STDOCB(L,1)
	        ENDDO
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,STDOCB(L,1)
              IF(ISO.GT.0) GOTO 856
	        DO K=2,KB
	          STDOCB(L,K)=STDOCB(L,1)
	        ENDDO
            ENDDO
          ENDIF
        ENDIF
	  IF(IREAD.EQ.1)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STDOCB(L,K),K=1,KB)
              IF(ISO.GT.0) GOTO 856
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                          (STDOCB(L,K),K=1,KB)
              IF(ISO.GT.0) GOTO 856
            ENDDO
          ENDIF
        ENDIF
	  IF(IREAD.EQ.2)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STDOCB(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 856
	        DO K=KBINPUT,KB
	          STDOCB(L,K)=STDOCB(L,KBINPUT)
	        ENDDO
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                          (STDOCB(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 856
	        DO K=KBINPUT,KB
	          STDOCB(L,K)=STDOCB(L,KBINPUT)
	        ENDDO
            ENDDO
          ENDIF
        ENDIF
        CLOSE(1)
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE POCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
C **  PARTICULATE ORGANIC CARBON IN BED
C
      IVAL=0
	DO NT=1,NTOX
	  IF(ISTOC(NT).EQ.1)IVAL=1
	ENDDO
C
      IF(IVAL.EQ.1)THEN
C
      IF(ISTPOCB.EQ.1)THEN
        OPEN(1,FILE='POCB.INP',STATUS='UNKNOWN')
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)ISALTYP,IREAD,KBINPUT
        DO K=1,KB
	    DO L=2,LC-1
	      STPOCB(L,K)=0.0
	    ENDDO
	  ENDDO
	  IF(IREAD.EQ.0)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) STPOCB(L,1)
              IF(ISO.GT.0) GOTO 856
	        DO K=2,KB
	          STPOCB(L,K)=STPOCB(L,1)
	        ENDDO
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,STPOCB(L,1)
              IF(ISO.GT.0) GOTO 856
	        DO K=2,KB
	          STPOCB(L,K)=STPOCB(L,1)
	        ENDDO
            ENDDO
          ENDIF
        ENDIF
	  IF(IREAD.EQ.1)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STPOCB(L,K),K=1,KB)
              IF(ISO.GT.0) GOTO 856
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                          (STPOCB(L,K),K=1,KB)
              IF(ISO.GT.0) GOTO 856
            ENDDO
          ENDIF
        ENDIF
	  IF(IREAD.EQ.2)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STPOCB(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 856
	        DO K=KBINPUT,KB
	          STPOCB(L,K)=STPOCB(L,KBINPUT)
	        ENDDO
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                          (STPOCB(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 856
	        DO K=KBINPUT,KB
	          STPOCB(L,K)=STPOCB(L,KBINPUT)
	        ENDDO
            ENDDO
          ENDIF
        ENDIF
        CLOSE(1)
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  OPEN FILE FPOCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
C **  PARTICULATE ORGANIC CARBON FRACTION FOR EACH SEDIMENT CLASS
C **  IN BED
C
      IVAL=0
	DO NT=1,NTOX
	  IF(ISTOC(NT).EQ.2.OR.ISTOC(NT).EQ.3)IVAL=1
	ENDDO
C
      IF(IVAL.EQ.1)THEN
C
      IF(ISTPOCB.EQ.3)THEN
        OPEN(1,FILE='FPOCB.INP',STATUS='UNKNOWN')
	  DO NS=1,NSED+NSND
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)ISALTYP,IREAD,KBINPUT
        DO K=1,KB
	    DO L=2,LC-1
	      STFPOCB(L,K,NS)=0.0
	    ENDDO
	  ENDDO
	  IF(IREAD.EQ.0)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) STFPOCB(L,1,NS)
              IF(ISO.GT.0) GOTO 856
	        DO K=2,KB
	          STFPOCB(L,K,NS)=STFPOCB(L,1,NS)
	        ENDDO
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,STFPOCB(L,1,NS)
              IF(ISO.GT.0) GOTO 856
	        DO K=2,KB
	          STFPOCB(L,K,NS)=STFPOCB(L,1,NS)
	        ENDDO
            ENDDO
          ENDIF
        ENDIF
	  IF(IREAD.EQ.1)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STFPOCB(L,K,NS),K=1,KB)
              IF(ISO.GT.0) GOTO 856
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                          (STFPOCB(L,K,NS),K=1,KB)
              IF(ISO.GT.0) GOTO 856
            ENDDO
          ENDIF
        ENDIF
	  IF(IREAD.EQ.2)THEN
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STFPOCB(L,K,NS),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 856
	        DO K=KBINPUT,KB
	          STFPOCB(L,K,NS)=STFPOCB(L,KBINPUT,NS)
	        ENDDO
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                          (STFPOCB(L,K,NS),K=1,KBINPUT,NS)
              IF(ISO.GT.0) GOTO 856
	        DO K=KBINPUT,KB
	          STFPOCB(L,K,NS)=STFPOCB(L,KBINPUT,NS)
	        ENDDO
            ENDDO
          ENDIF
        ENDIF
	  ENDDO
        CLOSE(1)
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C###########################################################################
C HQI Change to include sptially varying, but time constant bulk foc
C FPOCB  - Bulk foc from data
C PFPOCB - Pseudo foc from data, to be used for all partitioning calculations
C RM, 02/29/04
C**********************************************************************C
C
C **  OPEN FILE FOCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
C **  PARTICULATE ORGANIC CARBON IN BED AND PSEUDO-POC IN BED
C
      IF(ISTPOCB.EQ.4)THEN
C
        OPEN(1,FILE='FOCB.INP',STATUS='UNKNOWN')
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)ISALTYP,IREAD,KBINPUT
        DO K=1,KB
          DO L=2,LC-1
            FPOCB(L,K)=0.0
          ENDDO
        ENDDO
        DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(FPOCB(L,K),K=1,KBINPUT)
           DO K=1,KBINPUT
              FPOCB(L,K) = FPOCB(L,K)/1000000.
           ENDDO
           DO K=KBINPUT+1,KB
              FPOCB(L,K)=FPOCB(L,KBINPUT)
           ENDDO
           IF(ISO.GT.0) GOTO 856
        ENDDO
        CLOSE(1)
C
        OPEN(1,FILE='PSEUDO_FOCB.INP',STATUS='UNKNOWN')
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)ISALTYP,IREAD,KBINPUT
        DO K=1,KB
          DO L=2,LC-1
            PFPOCB(L,K)=0.0
          ENDDO
        ENDDO
        DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(PFPOCB(L,K),K=1,KBINPUT)
           DO K=1,KBINPUT
              PFPOCB(L,K) = PFPOCB(L,K)/1000000.
           ENDDO
           DO K=KBINPUT+1,KB
              PFPOCB(L,K)=PFPOCB(L,KBINPUT)
           ENDDO
           IF(ISO.GT.0) GOTO 856
        ENDDO
        CLOSE(1)
      ENDIF
C
C###########################################################################
C**********************************************************************C
C
C **  WRITE NEW DXDY FILES TO UPDATE OLD VERSIONS OF THE MODEL
C **  OR INCORPORATE MODIFIED DX'S AND DY'S DIRECTLY INTO
C **  DXDY.INP FILE
C
C     DXTMP=20.
C     DYTMP=20.
      OPEN(1,FILE='NEWDXDY.INP',STATUS='UNKNOWN')
      DO J=1,JC
      DO I=1,IC
      L=LIJ(I,J)
      IF(IJCT(I,J).GE.1.AND.IJCT(I,J).LT.9)THEN
        WRITE(1,339)IL(L),JL(L),DXP(L),DYP(L),HMP(L),BELV(L),
     &              ZBR(L)
C       WRITE(1,339)I,J,DXTMP,DYTMP,H(I,J),BELVIJ(I,J),
C    &              ZBRIJ(I,J)
      ENDIF
      ENDDO
      ENDDO
      CLOSE(1)
C
C**********************************************************************C
C
C **  READ IN INITIAL SALINITY, TEMPERATURE, DYE, SED, SND, TOX
C **  FOR COLD STARTS FORM FILE XXXX.INP
C
C----------------------------------------------------------------------C
C
C **  SALINITY
C
      DO K=1,KC
      DO L=2,LA
       SALINIT(L,K)=0.
      ENDDO
      ENDDO
C
      IF(ISTRAN(1).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(1).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(1).GE.1)THEN
        OPEN(1,FILE='SALT.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,4
         READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (SALINIT(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 840
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SALINIT(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 840
          ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  TEMPERATURE
C
      DO K=1,KC
      DO L=2,LA
       TEMINIT(L,K)=TEMO
      ENDDO
      ENDDO
C
      IF(ISTRAN(2).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(2).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(2).EQ.1)THEN
C
        OPEN(1,FILE='TEMP.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,4
         READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
          READ(1,*,IOSTAT=ISO) (TEMINIT(L,K),K=1,KC)
          IF(ISO.GT.0) GOTO 842
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(TEMINIT(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 842
          ENDDO
        ENDIF
        CLOSE(1)
C
        IF(ISRESTI.EQ.0.AND.ISBEDTEMI.EQ.0)THEN
          OPEN(1,FILE='TEMPB.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
          DO IS=1,4
            READ(1,1)
          ENDDO
          READ(1,*)ISALTYP
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (TEMB(L,K),K=1,KBH)
              IF(ISO.GT.0) GOTO 842
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(TEMB(L,K),K=1,KBH)
              IF(ISO.GT.0) GOTO 842
            ENDDO
          ENDIF
          CLOSE(1)
        ENDIF
C
      ENDIF
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  DYE
C
      DO K=1,KC
      DO L=2,LA
       DYEINIT(L,K)=0.
      ENDDO
      ENDDO
C
C
      IF(ISTRAN(3).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(3).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(3).GE.1)THEN
        OPEN(1,FILE='DYE.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,4
        READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (DYEINIT(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 844
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(DYEINIT(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 844
          ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  SFL
C
      DO K=1,KC
      DO L=2,LA
       SFLINIT(L,K)=0.
      ENDDO
      ENDDO
C
      IF(ISTRAN(4).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(4).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(4).GE.1)THEN
        OPEN(1,FILE='SFL.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,4
        READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (SFLINIT(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 846
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SFLINIT(L,K),K=1,KC)
           IF(ISO.GT.0) GOTO 846
          ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  TOXICS
C
      IF(ISTRAN(5).EQ.0) THEN
C
      DO NT=1,NTOX
      DO K=1,KC
      DO L=2,LA
       TOXINIT(L,K,NT)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      DO NT=1,NTOX
      DO K=1,KB
      DO L=2,LA
       TOXBINIT(L,K,NT)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C
      IF(ISTRAN(5).GE.1) THEN
C
      DO NT=1,NTOX
      DO K=1,KC
      DO L=2,LA
       TOXINIT(L,K,NT)=TOXINTW(NT)
      ENDDO
      ENDDO
      ENDDO
C
      DO NT=1,NTOX
      DO K=1,KB
      DO L=2,LA
       TOXBINIT(L,K,NT)=TOXINTB(NT)
      ENDDO
      ENDDO
      ENDDO
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(5).EQ.0) IISTMP=0
C
      IF(IISTMP.EQ.0.AND.ISTRAN(5).GE.1)THEN
      IF(ISLTMT.EQ.0)THEN
        OPEN(1,FILE='TOXW.INP',STATUS='UNKNOWN')
        IF(ITXINT(1).EQ.1.OR.ITXINT(1).EQ.3)THEN
C        IF(ITXINT(NT).EQ.1.OR.ITXINT(NT).EQ.3)THEN
        DO NT=1,NTOX
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
        READ(1,1)
        ENDDO
        READ(1,*)ISALTYP,ITOXWU(NT)
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (TOXINIT(L,K,NT),K=1,KC)
           IF(ISO.GT.0) GOTO 848
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(TOXINIT(L,K,NT),K=1,KC)
           IF(ISO.GT.0) GOTO 848
          ENDDO
        ENDIF
        ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
      ENDIF
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(5).EQ.0) IISTMP=0
C
      IF(IISTMP.EQ.0.AND.ISTRAN(5).GE.1)THEN
      IF(ISLTMT.EQ.0.)THEN
        OPEN(1,FILE='TOXB.INP',STATUS='UNKNOWN')
        IF(ITXINT(1).EQ.2.OR.ITXINT(1).EQ.3)THEN
C        IF(ITXINT(NT).EQ.2.OR.ITXINT(NT).EQ.3)THEN
        DO NT=1,NTOX
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
        READ(1,1)
        ENDDO
        READ(1,*)ISALTYP,ITOXBU(NT),KBINPUT
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (TOXBINIT(L,K,NT),K=1,KBINPUT)
           IF(ISO.GT.0) GOTO 852
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                        (TOXBINIT(L,K,NT),K=1,KBINPUT)
           IF(ISO.GT.0) GOTO 852
          ENDDO
        ENDIF
        ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT
C
      IF(ISTRAN(6).EQ.0)THEN
C
      DO NS=1,NSED
      DO K=1,KC
      DO L=2,LA
       SEDINIT(L,K,NS)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      DO NS=1,NSED
      DO K=1,KB
      DO L=2,LA
       SEDBINIT(L,K,NS)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      END IF
C
      IF(ISTRAN(6).GE.1)THEN
C
      DO NS=1,NSED
      DO K=1,KC
      DO L=2,LA
       SEDINIT(L,K,NS)=SEDO(NS)
      ENDDO
      ENDDO
      ENDDO
C
      DO NS=1,NSED
      DO K=1,KB
      DO L=2,LA
       SEDBINIT(L,K,NS)=SEDBO(NS)
      ENDDO
      ENDDO
      ENDDO
C
      ITXINTT=0
      IF(ISEDINT.EQ.1) ITXINTT=1
      IF(ISEDINT.EQ.3) ITXINTT=1
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
C
      IF(IISTMP.EQ.0.AND.ISTRAN(6).GE.1)THEN
      IF(ISLTMT.EQ.0.AND.ITXINTT.GE.1)THEN
        OPEN(1,FILE='SEDW.INP',STATUS='UNKNOWN')
        DO NS=1,NSED
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
        READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (SEDINIT(L,K,NS),K=1,KC)
           IF(ISO.GT.0) GOTO 854
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SEDINIT(L,K,NS),K=1,KC)
           IF(ISO.GT.0) GOTO 854
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
      ENDIF
C
      ITXINTT=0
      IF(ISEDINT.EQ.2) ITXINTT=1
      IF(ISEDINT.EQ.3) ITXINTT=1
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
C
      IF(IISTMP.EQ.0.AND.ISTRAN(6).GE.1)THEN
      IF(ISLTMT.EQ.0.AND.ITXINTT.GE.1)THEN
        OPEN(1,FILE='SEDB.INP',STATUS='UNKNOWN')
        DO NS=1,NSED
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
        READ(1,1)
        ENDDO
        READ(1,*)ISALTYP,ISEDBU(NS),KBINPUT
        DO K=1,KB
	    DO L=2,LC-1
	      SEDBINIT(L,K,NS)=0.0
	    ENDDO
	  ENDDO
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (SEDBINIT(L,K,NS),K=1,KBINPUT)
           IF(ISO.GT.0) GOTO 856
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &                          (SEDBINIT(L,K,NS),K=1,KBINPUT)
           IF(ISO.GT.0) GOTO 856
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  NON-COHESIVE SEDIMENT
C
      IF(ISTRAN(7).EQ.0)THEN
C
      DO NX=1,NSND
      DO K=1,KC
      DO L=2,LA
       SNDINIT(L,K,NX)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      DO NX=1,NSND
      DO K=1,KB
      DO L=2,LA
       SNDBINIT(L,K,NX)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      END IF
C
      IF(ISTRAN(7).GE.1)THEN
C
      DO NX=1,NSND
      NS=NX+NSED
      DO K=1,KC
      DO L=2,LA
       SNDINIT(L,K,NX)=SEDO(NS)
      ENDDO
      ENDDO
      ENDDO
C
      DO NX=1,NSND
      NS=NX+NSED
      DO K=1,KB
      DO L=2,LA
       SNDBINIT(L,K,NX)=SEDBO(NS)
      ENDDO
      ENDDO
      ENDDO
C
      ITXINTT=0
      IF(ISEDINT.EQ.1) ITXINTT=1
      IF(ISEDINT.EQ.3) ITXINTT=1
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
C
      IF(IISTMP.EQ.0.AND.ISTRAN(7).GE.1)THEN
      IF(ISLTMT.EQ.0.AND.ITXINTT.GE.1)THEN
        OPEN(1,FILE='SNDW.INP',STATUS='UNKNOWN')
        DO NX=1,NSND
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
        READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (SNDINIT(L,K,NX),K=1,KC)
           IF(ISO.GT.0) GOTO 858
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SNDINIT(L,K,NX),K=1,KC)
           IF(ISO.GT.0) GOTO 858
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
      ENDIF
C
      ITXINTT=0
      IF(ISEDINT.EQ.2) ITXINTT=1
      IF(ISEDINT.EQ.3) ITXINTT=1
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
C
      IF(IISTMP.EQ.0.AND.ISTRAN(7).GE.1)THEN
      IF(ISLTMT.EQ.0.AND.ITXINTT.GE.1)THEN
        OPEN(1,FILE='SNDB.INP',STATUS='UNKNOWN')
        DO NX=1,NSND
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
        READ(1,1)
        ENDDO
        DO K=1,KB
	    DO L=2,LC-1
	      SNDBINIT(L,K,NX)=0.0
	    ENDDO
	  ENDDO
        READ(1,*)ISALTYP,ISNDBU(NX),KBINPUT
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (SNDBINIT(L,K,NX),K=1,KBINPUT)
           IF(ISO.GT.0) GOTO 862
          ENDDO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &     (SNDBINIT(L,K,NX),K=1,KBINPUT)
           IF(ISO.GT.0) GOTO 862
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C  ** SEDIMENT BED MECHANICAL INITIAL CONDITIONS
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
C
      IF(ISTRAN(6).GT.0.OR.ISTRAN(7).GT.0) THEN
C
C  ** BED LAYER THICKNESS
C
      IF(IISTMP.EQ.0.AND.IBMECH.GE.1)THEN
        OPEN(1,FILE='BEDLAY.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)IBEDLAYU,ISALTYP,KBINPUT
        IF(IBEDLAYU.GT.0) THEN
	    DO K=1,KB
	      DO L=2,LC-1
	        BEDLINIT(L,K)=0.0
	      ENDDO
	    ENDDO
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (BEDLINIT(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 862
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &        (BEDLINIT(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 862
            ENDDO
          ENDIF
        ENDIF
        CLOSE(1)
      ENDIF
C
C  ** BED LAYER BULK DENSITY
C
      IF(IISTMP.EQ.0.AND.IBMECH.GE.1)THEN
        OPEN(1,FILE='BEDBDN.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
          READ(1,1)
        ENDDO
        READ(1,*)IBEDBDNU,ISALTYP,KBINPUT
        IF(IBEDBDNU.GT.0)THEN
	    DO K=1,KB
	      DO L=2,LC-1
	        BEDBINIT(L,K)=0.0
	      ENDDO
	    ENDDO
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (BEDBINIT(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 862
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &        (BEDBINIT(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 862
            ENDDO
          ENDIF
        ENDIF
        CLOSE(1)
      ENDIF
C
C  ** BED LAYER DRY DENSITY, POROSITY OR VOID RATIO
C
      IF(IISTMP.EQ.0.AND.IBMECH.GE.1)THEN
        OPEN(1,FILE='BEDDDN.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,8
        READ(1,1)
        ENDDO
        READ(1,*)IBEDDDNU,ISALTYP,KBINPUT
        IF(IBEDDDNU.GT.0)THEN
	    DO K=1,KB
	      DO L=2,LC-1
	        BEDDINIT(L,K)=0.0
	      ENDDO
	    ENDDO
          IF(ISALTYP.EQ.0)THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (BEDDINIT(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 862
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &        (BEDDINIT(L,K),K=1,KBINPUT)
              IF(ISO.GT.0) GOTO 862
             ENDDO
          ENDIF
        ENDIF
        CLOSE(1)
      ENDIF
C
C  ** CONSOLIDATION MAP
C
      IF(IBMECH.EQ.9)THEN
        OPEN(1,FILE='CONSOLMAP.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,6
        READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) LCONSOL(L)
            IF(ISO.GT.0) GOTO 862
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,LCONSOL(L)
            IF(ISO.GT.0) GOTO 862
          ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
   19 FORMAT (/,' INITIAL BUOYANCY ARRAY:',//)
  907 FORMAT(12F6.2)
C
C**********************************************************************C
C
C **  READ IN HORIZONTAL VELOCITY TIME SERIES FOR DATA ASSIMILATION
C **  FILE UVSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NUVSER.GE.1)THEN
      OPEN(1,FILE='UVSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,19
      READ(1,1)
      ENDDO
C
      DO NS=1,NUVSER
        READ(1,*,IOSTAT=ISO) ITYPE,MUVSER(NS),TCUVSER(NS),TAUVSER(NS),
     &                   RMULADJ,ADDADJ
        IF(ISO.GT.0) GOTO 850
        IF(ITYPE.EQ.-1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 860
          DO M=1,MUVSER(NS)
            READ(1,*,IOSTAT=ISO)TUVSER(M,NS),USERTMP,VSERTMP
            IF(ISO.GT.0) GOTO 860
            TUVSER(M,NS)=TUVSER(M,NS)+TAUVSER(NS)
            USERTMP=(RMULADJ*(USERTMP+ADDADJ))
            VSERTMP=(RMULADJ*(VSERTMP+ADDADJ))
            DO K=1,KC
			USER(M,K,NS)=USERTMP*WKQ(K)
			VSER(M,K,NS)=VSERTMP*WKQ(K)
            ENDDO
          ENDDO
        ENDIF
        IF(ITYPE.EQ.0)THEN
          DO M=1,MUVSER(NS)
            READ(1,*,IOSTAT=ISO)TUVSER(M,NS),(USER(M,K,NS),K=1,KC),
     &                                      (VSER(M,K,NS),K=1,KC)
            IF(ISO.GT.0) GOTO 860
            TUVSER(M,NS)=TUVSER(M,NS)+TAUVSER(NS)
            DO K=1,KC
              USER(M,K,NS)=RMULADJ*(USER(M,K,NS)+ADDADJ)
              VSER(M,K,NS)=RMULADJ*(VSER(M,K,NS)+ADDADJ)
            ENDDO
          ENDDO
        ENDIF
        IF(ITYPE.GE.1)THEN
          DO M=1,MUVSER(NS)
            READ(1,*,IOSTAT=ISO)TUVSER(M,NS),USERTMP,VSERTMP
            IF(ISO.GT.0) GOTO 860
            TUVSER(M,NS)=TUVSER(M,NS)+TAUVSER(NS)
            USERTMP=(RMULADJ*(USERTMP+ADDADJ))
            VSERTMP=(RMULADJ*(VSERTMP+ADDADJ))
            DO K=1,KC
              USER(M,K,NS)=-99999.9
              VSER(M,K,NS)=-99999.9
            ENDDO
	      USER(M,ITYPE,NS)=USERTMP
            VSER(M,ITYPE,NS)=VSERTMP
          ENDDO
        ENDIF
      ENDDO
C
      CLOSE(1)
      ENDIF
C
 7001 FORMAT(2I6,F12.5,50F10.5)
C
C**********************************************************************C
C
C **  READ IN HORIZONTAL VELOCITY SPATIAL FIELD FOR DATA ASSIMILATION
C **  FILE UVDAMAP.INP
C
C----------------------------------------------------------------------C
C
      IF(NLUVDA.GT.0.AND.NUVSER.eq.0)THEN
      OPEN(1,FILE='UVDAMAP.INP',STATUS='UNKNOWN')
      OPEN(99,FILE='UVDAMAP.OUT',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,7
      READ(1,1)
      ENDDO
C
      UVMAX=0.
	IMAX=0
	JMAX=0
      NL=0
      DO NLT=1,NLUVDA
	  READ(1,*)ITMP,JTMP,USTMP,VSTMP
	  IF(IJCT(ITMP,JTMP).NE.0.AND.IJCT(ITMP,JTMP).NE.9)THEN
	    NL=NL+1
	    ICUVDA(NL)=ITMP
		JCUVDA(NL)=JTMP
		USERT(1,NL)=USTMP
		VSERT(1,NL)=VSTMP
	    UVTMP=SQRT(USTMP*USTMP+VSTMP*VSTMP)
	    WRITE(99,1999)ITMP,JTMP,USTMP,VSTMP,UVTMP
	    IF(UVTMP.GT.UVMAX)THEN
	      UVMAX=UVTMP
	      IMAX=ITMP
	      JMAX=JTMP
	    ENDIF
	    DO K=1,KC
	      USERT(K,NL)=USERT(1,NL)
	      VSERT(K,NL)=VSERT(1,NL)
          ENDDO
          NUVSERA(NL)=NL
          TSUVDA(NL)=TSUVDA(1)
          FSUVDA(NL)=FSUVDA(1)
          IWUVDA(NL)=IWUVDA(1)
          IRVUDA(NL)=IRVUDA(1)
          RRUVDA(NL)=RRUVDA(1)
        ENDIF
      ENDDO
	NLUVDA=NL
!	WRITE(8,*)'I,J,UVMAX=',IMAX,JMAX,UVMAX  !hnr 7/27/2009
C
      CLOSE(1)
	CLOSE(99)
      ENDIF
C
 1999 FORMAT(2I6,3F10.3)
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY SURFACE ELEVATION TIME SERIES FROM THE
C **  FILE PSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NPSER.GE.1)THEN
      OPEN(1,FILE='PSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,10
      READ(1,1)
      ENDDO
C
      DO NS=1,NPSER
        READ(1,*,IOSTAT=ISO) ITYPE,MPSER(NS),TCPSER(NS),TAPSER(NS),
     &                   RMULADJ,ADDADJ,PSERZDF(NS)
        IF(ISO.GT.0) GOTO 850
        PSERZDF(NS)=G*PSERZDF(NS)
	  IF(ITYPE.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) RMULADJS,ADDADJS,PSERZDS(NS)
	    IF(ISO.GT.0) GOTO 850
          PSERZDS(NS)=G*PSERZDS(NS)
        ELSE
          RMULADJS=0
		ADDADJS=0.
		PSERZDS(NS)=0.
        ENDIF
        IF(ITYPE.EQ.0)THEN
          DO M=1,MPSER(NS)
            READ(1,*,IOSTAT=ISO)TPSER(M,NS),PSERTMP
            IF(ISO.GT.0) GOTO 850
            TPSER(M,NS)=TPSER(M,NS)+TAPSER(NS)
            PSER(M,NS)=G*(PSERTMP+ADDADJ)*RMULADJ
            PSERS(M,NS)=0.0
	    ENDDO
        ENDIF
        IF(ITYPE.EQ.1)THEN
          DO M=1,MPSER(NS)
            READ(1,*,IOSTAT=ISO)TPSER(M,NS),PSERTMP,PSERTMPS
            IF(ISO.GT.0) GOTO 850
            TPSER(M,NS)=TPSER(M,NS)+TAPSER(NS)
            PSER(M,NS)=G*(PSERTMP+ADDADJ)*RMULADJ
            PSERS(M,NS)=G*(PSERTMPS+ADDADJS)*RMULADJS
	    ENDDO
        ENDIF
      ENDDO
C
      CLOSE(1)
      ENDIF
C
CDIA      WRITE(6,6776)'READ PSER.INP'
C
 6776 FORMAT(A20)
C
C**********************************************************************C
C
C **  READ IN VOLUMETRIC SOURCE OR RIVER INFLOW TIME SERIES FROM THE
C **  FILE QSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NQSER.GE.1)THEN
      OPEN(1,FILE='QSER.INP',STATUS='UNKNOWN')
C      OPEN(2,FILE='QSER.CK',STATUS='UNKNOWN')
C      CLOSE(2,STATUS='DELETE')
C      OPEN(2,FILE='QSER.CK',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,14
        READ(1,1)
        ENDDO
C
        DO NS=1,NQSER
         ISMOOTH=0
C        NTMP=0
C **  FIX FOR EXTRA HEADERS IN SOME QSER.INP FILES
         READ(1,6789)CHARSPACE
	   IF(CHARSPACE.NE.CHARSKP1)THEN
	     IF(CHARSPACE.NE.CHARSKP2) BACKSPACE(1)
	   ENDIF
         READ(1,6789)CHARSPACE
	   IF(CHARSPACE.NE.CHARSKP1)THEN
	     IF(CHARSPACE.NE.CHARSKP2) BACKSPACE(1)
	   ENDIF
C
        READ(1,*,IOSTAT=ISO)ISTYP, MQSER(NS),TCQSER(NS),TAQSER(NS),
     &                   RMULADJ,ADDADJ,ICHGQS
        IF(MQSER(NS).LT.0)THEN
           ISMOOTH=1
           MQSER(NS)=-MQSER(NS)
        ENDIF
        IF(ISO.GT.0) GOTO 860
      if(mqser(ns).gt.ndqser) then                                     !hnr
        write(*,370)mqser(ns),ns,ndqser                                !hnr
370     FORMAT('The number of data pairs in your application',I8,      !hnr
     $  '  for Qser time series',I6,                                   !hnr
     &  '  is larger than the maximum allowed of',I8,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop qser file'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 860
           DO M=1,MQSER(NS)
C           NTMP=NTMP+1
           READ(1,*,IOSTAT=ISO)TQSER(M,NS),QSERTMP
C          WRITE(2,2222)NS,NTMP,TQSER(M,NS),QSERTMP
           IF(ISO.GT.0) GOTO 860
           TQSER(M,NS)=TQSER(M,NS)+TAQSER(NS)
           QSERTMP=(RMULADJ*(QSERTMP+ADDADJ))
           IF(ICHGQS.EQ.1) QSERTMP=MAX(QSERTMP,0.0)
           IF(ICHGQS.EQ.-1) QSERTMP=MIN(QSERTMP,0.0)
            DO K=1,KC
            QSER(M,K,NS)=QSERTMP*WKQ(K)
            ENDDO
           ENDDO
         ELSE
          DO M=1,MQSER(NS)
          READ(1,*,IOSTAT=ISO)TQSER(M,NS),(QSER(M,K,NS), K=1,KC)
          IF(ISO.GT.0) GOTO 860
          TQSER(M,NS)=TQSER(M,NS)+TAQSER(NS)
           DO K=1,KC
           QSER(M,K,NS)=RMULADJ*(QSER(M,K,NS)+ADDADJ)
           IF(ICHGQS.EQ.1) QSER(M,K,NS)=MAX(QSER(M,K,NS),0.0)
           IF(ICHGQS.EQ.-1) QSER(M,K,NS)=MIN(QSER(M,K,NS),0.0)
           ENDDO
          ENDDO
        ENDIF
        IF(ISMOOTH.EQ.1)THEN
          DO K=1,KC
          DO M=1,MQSER(NS)
            QSERSM(N,K)=QSER(M,K,NS)
          ENDDO
          ENDDO
          DO K=1,KC
          DO M=2,MQSER(NS)-1
            QSER(N,K,NS)=0.25*(QSERSM(M-1,K)+QSERSM(M+1,K))
     &                  +0.5*QSERSM(M,K)
          ENDDO
          ENDDO
        ENDIF
        ENDDO
C
      CLOSE(1)
C      CLOSE(2)
      ENDIF
C
 2222 FORMAT(2I5,F12.7,F12.4)
 6789 FORMAT(A1)
C
CDIA      WRITE(6,6776)'READ QSER.INP'
C
C**********************************************************************C
C
C **  READ IN FLOW WITHDRAWL-RETURN FLOW AND CONCENTRATION RISE
C **  TIME SERIES FROM THE FILE QWRS.INP
C
C----------------------------------------------------------------------C
C
      IF(NQWRSR.GE.1)THEN
      OPEN(1,FILE='QWRS.INP',STATUS='UNKNOWN')
C      OPEN(2,FILE='QSER.CK',STATUS='UNKNOWN')
C      CLOSE(2,STATUS='DELETE')
C      OPEN(2,FILE='QSER.CK',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        NCTMP=4+NSED+NSND+NTOX
C
        DO IS=1,16
        READ(1,1)
        ENDDO
C
        DO NS=1,NQWRSR
C        NTMP=0
        READ(1,*,IOSTAT=ISO)ISTYP,MQWRSR(NS),TCQWRSR(NS),TAQWRSR(NS),
     &                   RMULADJ,ADDADJ
        IF(ISO.GT.0) GOTO 865
        IF(ISTYP.EQ.0)THEN
          DO NC=1,NCTMP
           DO M=1,MQWRSR(NS)
            CQWRSER(M,NS,NC)=0.
           ENDDO
          ENDDO
          DO M=1,MQWRSR(NS)
           READ(1,*,IOSTAT=ISO)TQWRSER(M,NS),QWRSER(M,NS)
           IF(ISO.GT.0) GOTO 865
           TQWRSER(M,NS)=TQWRSER(M,NS)+TAQWRSR(NS)
           QWRSER(M,NS)=(RMULADJ*(QWRSER(M,NS)+ADDADJ))
          ENDDO
         ELSE
          DO M=1,MQWRSR(NS)
           READ(1,*,IOSTAT=ISO)TQWRSER(M,NS),QWRSER(M,NS),
     &                  (CQWRSER(M,NS,NC),NC=1,NCTMP)
           IF(ISO.GT.0) GOTO 865
           TQWRSER(M,NS)=TQWRSER(M,NS)+TAQWRSR(NS)
           QWRSER(M,NS)=(RMULADJ*(QWRSER(M,NS)+ADDADJ))
          ENDDO
        ENDIF
        ENDDO
C
      CLOSE(1)
C
      ENDIF
C
CDIA      WRITE(6,6776)'READ QWRS.INP'
C
C**********************************************************************C
C
C **  READ IN GROUNDWATER INFLOW/OUTFLOW AND CONCENTRATION TIME
C **  SERIES FROM THE FILE GWSER.INP
C
C----------------------------------------------------------------------C
C
      IF(ISGWIT.EQ.2)THEN
        OPEN(1,FILE='GWSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        NCTMP=4+NSED+NSND+NTOX
C
        DO IS=1,14
        READ(1,1)
        ENDDO
C
        READ(1,*)NGWSER
        IF(NGWSER.GT.0)THEN
          DO NS=1,NGWSER
            READ(1,*,IOSTAT=ISO)MGWSER(NS),TCGWSER(NS),TAGWSER(NS),
     &                   RMULADJ,ADDADJ
            IF(ISO.GT.0) GOTO 865
            DO M=1,MGWSER(NS)
              READ(1,*,IOSTAT=ISO)TGWSER(M,NS),GWSER(M,NS),
     &                  (GWCSER(M,NS,NC),NC=1,NCTMP)
              IF(ISO.GT.0) GOTO 865
              TGWSER(M,NS)=TGWSER(M,NS)+TAGWSER(NS)
              GWSER(M,NS)=(RMULADJ*(GWSER(M,NS)+ADDADJ))
            ENDDO
          ENDDO
        ENDIF
C
        CLOSE(1)
C
      ENDIF
C
C
C**********************************************************************C
C
C **  READ IN SPATIAL MAPS AND TIME SERIES FOR EXTERNAL SPECIFICATION OF DISSOLVED AND
C **  PARTICULATE ORGANIC CARBON FOR USE IN TOXIC CONTAMINANT SORPTION
C
C----------------------------------------------------------------------C
C
C **  DISSOLVED ORGANIC CARBON
c      IF(ISTOC.GE.1)THEN
c        OPEN(1,FILE='OCSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
C
c        DO IS=1,14
c        READ(1,1)
c        ENDDO
C
c        READ(1,*)NOCSER
c        IF(NOCSER.GT.0)THEN
c          DO NS=1,NOCSER
c            READ(1,*,IOSTAT=ISO)MOCSER(NS),TCOCSER(NS),TAOCSER(NS),
c     &                   RMULADJ,ADDADJ
c            IF(ISO.GT.0) GOTO 865
c            DO M=1,MCOSER(NS)
c              READ(1,*,IOSTAT=ISO)TOCSER(M,NS),DOCWSER(M,NS),
c     &          DOCBSER(M,NS),POCWSER(M,NS),POCBSER(M,NS)
c              IF(ISO.GT.0) GOTO 865
c              TOCSER(M,NS)=TOCSER(M,NS)+TACOSER(NS)
c              DOCWSER(M,NS)=(RMULADJ*(DOCWSER(M,NS)+ADDADJ))
c              POCWSER(M,NS)=(RMULADJ*(POCWSER(M,NS)+ADDADJ))
c              DOCBSER(M,NS)=(RMULADJ*(DOCBSER(M,NS)+ADDADJ))
c              POCWSER(M,NS)=(RMULADJ*(POCBSER(M,NS)+ADDADJ))
c            ENDDO
c          ENDDO
c        ENDIF
C
c        CLOSE(1)
C
c      ENDIF
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE SALINITY TIME SERIES
C **  FROM THE FILE SSER.INP
C
C----------------------------------------------------------------------C
C
 8888 format(3i5,2f10.2)
c
      IF(NCSER(1).GE.1)THEN
        OPEN(1,FILE='SSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      ENDDO
C
        NC=1
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     &                   TACSER(NS,NC),RMULADJ,ADDADJ
        IF(ISO.GT.0) GOTO 870
      if(mcser(ns,nc).gt.ndcser) then                                  !hnr
        write(*,390)mcser(ns,nc),ns,ndcser                             !hnr
390     FORMAT('The number of data pairs in your application',I8,      !hnr
     $  '  for Sser time series',I6,                                   !hnr
     &  '  is larger than the maximum allowed of',I8,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop sser file'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 870
           DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP

!	     if(m.eq.1)write(8,8888)nc,ns,m,TCSER(M,NS,NC),CSERTMP    !hnr 7/27/2009
!	     if(m.eq.MCSER(NS,NC))write(8,8888)nc,ns,m,              !hnr 7/27/2009
!     &                         TCSER(M,NS,NC),CSERTMP             !hnr 7/27/2009

           IF(ISO.GT.0) GOTO 870
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
            DO K=1,KC
            CSER(M,K,NS,NC)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            ENDDO
           ENDDO
         ELSE
          DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
          IF(ISO.GT.0) GOTO 870
          TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
           CSER(M,K,NS,NC)=RMULADJ*(CSER(M,K,NS,NC)+ADDADJ)
           ENDDO
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE TEMPERATURE TIME
C **  SERIES FROM THE FILE TSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NCSER(2).GE.1)THEN
        OPEN(1,FILE='TSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      ENDDO
C
        NC=2
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     &                   TACSER(NS,NC),RMULADJ,ADDADJ
        IF(ISO.GT.0) GOTO 880
      if(mcser(ns,nc).gt.ndcser) then                                  !hnr
        write(*,400)mcser(ns,nc),ns,ndcser                             !hnr
400     FORMAT('The number of data pairs in your application',I8,      !hnr
     $  '  for Tser time series',I6,                                   !hnr
     &  '  is larger than the maximum allowed of',I8,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop tser file'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 880
           DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP

!	     if(m.eq.1)write(8,8888)nc,ns,m,TCSER(M,NS,NC),CSERTMP   !hnr 7/27/2009
!	     if(m.eq.MCSER(NS,NC))write(8,8888)nc,ns,m,       !hnr 7/27/2009
!     &                         TCSER(M,NS,NC),CSERTMP              !hnr 7/27/2009

           IF(ISO.GT.0) GOTO 880
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
            DO K=1,KC
            CSER(M,K,NS,NC)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            ENDDO
           ENDDO
         ELSE
          DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
          IF(ISO.GT.0) GOTO 880
          TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
           CSER(M,K,NS,NC)=RMULADJ*(CSER(M,K,NS,NC)+ADDADJ)
           ENDDO
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE DYE CONCENTRATION
C **  TIME SERIES FROM THE FILE DSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NCSER(3).GE.1)THEN
        OPEN(1,FILE='DSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      ENDDO
C
        NC=3
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     &                   TACSER(NS,NC),RMULADJ,ADDADJ
        IF(ISO.GT.0) GOTO 890
      if(mcser(ns,nc).gt.ndcser) then                                  !hnr
        write(*,410)mcser(ns,nc),ns,ndcser                             !hnr
410     FORMAT('The number of data pairs in your application',I8,      !hnr
     $  '  for Dser time series',I6,                                   !hnr
     &  '  is larger than the maximum allowed of',I8,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop dser file'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 890
           DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP
           IF(ISO.GT.0) GOTO 890
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
            DO K=1,KC
            CSER(M,K,NS,NC)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            ENDDO
           ENDDO
         ELSE
          DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
          IF(ISO.GT.0) GOTO 890
          TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
           CSER(M,K,NS,NC)=RMULADJ*(CSER(M,K,NS,NC)+ADDADJ)
           ENDDO
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE COHESIVE SEDIMENT
C **  CONCENTRATION TIME SERIES FROM THE FILE SDSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NSED.GT.0)THEN
C
      NFSED=MSVSED(1)
      IF(NCSER(NFSED).GE.1)THEN
        OPEN(1,FILE='SDSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      ENDDO
C
        NC=MSVSED(1)
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     &                   TACSER(NS,NC),RMULADS(1),ADDADS(1)
        IF(ISO.GT.0) GOTO 900
        IF(NSED.GT.1)THEN
          DO NT=2,NSED
           READ(1,*,IOSTAT=ISO)RMULADS(NT),ADDADS(NT)
           IF(ISO.GT.0) GOTO 900
           NTT=NT-1
           MCSER(NS,NC+NTT)=MCSER(NS,NC)
           TCCSER(NS,NC+NTT)=TCCSER(NS,NC)
           TACSER(NS,NC+NTT)=TACSER(NS,NC)
          ENDDO
        ENDIF
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 900
          DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP
           IF(ISO.GT.0) GOTO 900
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
             CSER(M,K,NS,NC)=(RMULADS(1)*(CSERTMP+ADDADS(1)))*WKQ(K)
           ENDDO
           DO NT=2,NSED
            NTT=NT-1
            TCSER(M,NS,NC+NTT)=TCSER(M,NS,NC)
            READ(1,*,IOSTAT=ISO)CSERTMP
            IF(ISO.GT.0) GOTO 900
            DO K=1,KC
              CSER(M,K,NS,NC+NTT)
     &                      =(RMULADS(NT)*(CSERTMP+ADDADS(NT)))*WKQ(K)
            ENDDO
           ENDDO
          ENDDO
         ELSE
          DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
           IF(ISO.GT.0) GOTO 900
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
            CSER(M,K,NS,NC)=RMULADS(1)*(CSER(M,K,NS,NC)+ADDADS(1))
           ENDDO
           DO NT=2,NSED
            NTT=NT-1
            TCSER(M,NS,NC+NTT)=TCSER(M,NS,NC)
            READ(1,*,IOSTAT=ISO)(CSER(M,K,NS,NC+NTT), K=1,KC)
            IF(ISO.GT.0) GOTO 900
            DO K=1,KC
             CSER(M,K,NS,NC+NTT)
     &                   =RMULADS(NT)*(CSER(M,K,NS,NC+NTT)+ADDADS(NT))
            ENDDO
           ENDDO
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
      ENDIF
C
C **  CHECK SEDIMENT SERIES
C
C        OPEN(2,FILE='SEDCK.OUT')
C        CLOSE(2,STATUS='DELETE')
C        OPEN(2,FILE='SEDCK.OUT')
C        DO NC=MSVSED(1),MSVSED(NSED)
C         DO NS=1,NCSER(NC)
C          DO M=1,MCSER(NS,NC)
C           WRITE(2,2001)NC,NS,M,TCSER(M,NS,NC),CSER(M,1,NS,NC)
C          ENDDO
C         ENDDO
C        ENDDO
C        CLOSE(2)
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE NONCOHESIVE SEDIMENT
C **  CONCENTRATION TIME SERIES FROM THE FILE SNSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NSND.GT.0)THEN
C
      NFSND=MSVSND(1)
      IF(NCSER(NFSND).GE.1)THEN
        OPEN(1,FILE='SNSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      ENDDO
C
        NC=MSVSND(1)
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     &                   TACSER(NS,NC),RMULADS(1),ADDADS(1)
        IF(ISO.GT.0) GOTO 902
        IF(NSND.GT.1)THEN
          DO NT=2,NSND
           READ(1,*,IOSTAT=ISO)RMULADS(NT),ADDADS(NT)
           IF(ISO.GT.0) GOTO 902
           NTT=NT-1
           MCSER(NS,NC+NTT)=MCSER(NS,NC)
           TCCSER(NS,NC+NTT)=TCCSER(NS,NC)
           TACSER(NS,NC+NTT)=TACSER(NS,NC)
          ENDDO
        ENDIF
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 902
          DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP
           IF(ISO.GT.0) GOTO 902
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
             CSER(M,K,NS,NC)=(RMULADS(1)*(CSERTMP+ADDADS(1)))*WKQ(K)
           ENDDO
           DO NT=2,NSND
            NTT=NT-1
            TCSER(M,NS,NC+NTT)=TCSER(M,NS,NC)
            READ(1,*,IOSTAT=ISO)CSERTMP
            IF(ISO.GT.0) GOTO 902
            DO K=1,KC
              CSER(M,K,NS,NC+NTT)
     &                     =(RMULADS(NT)*(CSERTMP+ADDADS(NT)))*WKQ(K)
            ENDDO
           ENDDO
          ENDDO
         ELSE
          DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
           IF(ISO.GT.0) GOTO 902
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
            CSER(M,K,NS,NC)=RMULADS(1)*(CSER(M,K,NS,NC)+ADDADS(1))
           ENDDO
           DO NT=2,NSND
            NTT=NT-1
            TCSER(M,NS,NC+NTT)=TCSER(M,NS,NC)
            READ(1,*,IOSTAT=ISO)(CSER(M,K,NS,NC+NTT), K=1,KC)
            IF(ISO.GT.0) GOTO 902
            DO K=1,KC
             CSER(M,K,NS,NC+NTT)
     &                  =RMULADS(NT)*(CSER(M,K,NS,NC+NTT)+ADDADS(NT))
            ENDDO
           ENDDO
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
      ENDIF
C
C **  CHECK SEDIMENT SERIES
C
C        OPEN(2,FILE='SNDCK.OUT')
C        CLOSE(2,STATUS='DELETE')
C        OPEN(2,FILE='SNDCK.OUT')
C        DO NC=MSVSND(1),MSVSED(NSND)
C         DO NS=1,NCSER(NC)
C          DO M=1,MCSER(NS,NC)
C           WRITE(2,2001)NC,NS,M,TCSER(M,NS,NC),CSER(M,1,NS,NC)
C          ENDDO
C         ENDDO
C        ENDDO
C        CLOSE(2)
C
 2001 FORMAT(3I5,2F12.5)
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE TOXIC CONTAMINANT
C **  CONCENTRATION TIME SERIES FROM THE FILE TXSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NTOX.GT.0)THEN
C
      NFTOX=MSVTOX(1)
      IF(NCSER(NFTOX).GE.1)THEN
        OPEN(1,FILE='TXSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      ENDDO
C
        NC=MSVTOX(1)
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     &                   TACSER(NS,NC),RMULADS(1),ADDADS(1)
        IF(ISO.GT.0) GOTO 904
        IF(NTOX.GT.1)THEN
          DO NT=2,NTOX
           READ(1,*,IOSTAT=ISO)RMULADS(NT),ADDADS(NT)
           IF(ISO.GT.0) GOTO 904
           NTT=NT-1
           MCSER(NS,NC+NTT)=MCSER(NS,NC)
           TCCSER(NS,NC+NTT)=TCCSER(NS,NC)
           TACSER(NS,NC+NTT)=TACSER(NS,NC)
          ENDDO
        ENDIF
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 904
          DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP
           IF(ISO.GT.0) GOTO 904
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
             CSER(M,K,NS,NC)=(RMULADS(1)*(CSERTMP+ADDADJ))*WKQ(K)
           ENDDO
           DO NT=2,NTOX
            NTT=NT-1
            TCSER(M,NS,NC+NTT)=TCSER(M,NS,NC)
            READ(1,*,IOSTAT=ISO)CSERTMP
            IF(ISO.GT.0) GOTO 904
            DO K=1,KC
              CSER(M,K,NS,NC+NTT)
     &                     =(RMULADS(NT)*(CSERTMP+ADDADS(NT)))*WKQ(K)
            ENDDO
           ENDDO
          ENDDO
         ELSE
          DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
           IF(ISO.GT.0) GOTO 904
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
            CSER(M,K,NS,NC)=RMULADS(1)*(CSER(M,K,NS,NC)+ADDADJ)
           ENDDO
           DO NT=2,NTOX
            NTT=NT-1
            TCSER(M,NS,NC+NTT)=TCSER(M,NS,NC)
            READ(1,*,IOSTAT=ISO)(CSER(M,K,NS,NC+NTT), K=1,KC)
            IF(ISO.GT.0) GOTO 904
            DO K=1,KC
             CSER(M,K,NS,NC+NTT)=RMULADS(NT)*(CSER(M,K,NS,NC+NTT)
     &                                       +ADDADJ)
            ENDDO
           ENDDO
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE SHELL FISH LARVAE
C **  TIME SERIES FROM THE FILE SFSER.INP
C
C----------------------------------------------------------------------C
C
      IF(NCSER(4).GE.1)THEN
        OPEN(1,FILE='SFSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      ENDDO
C
        NC=7
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     &                   TACSER(NS,NC),RMULADJ,ADDADJ
        IF(ISO.GT.0) GOTO 910
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 910
           DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP
           IF(ISO.GT.0) GOTO 910
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
            DO K=1,KC
            CSER(M,K,NS,NC)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            ENDDO
           ENDDO
         ELSE
          DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
          IF(ISO.GT.0) GOTO 910
          TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
           DO K=1,KC
            CSER(M,K,NS,NC)=RMULADJ*(CSER(M,K,NS,NC)+ADDADJ)
           ENDDO
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
C*************************************************************************************************C
C
C ** START READ LOAD INPUT CONTROL FILE AND LOAD TIME SERIES FOR SALT, TEMP(HEAT), DYE, SHELLFISH LARVE 
C
C**********************************************************************C    START HNR_GHD ADDED ON 7/2023
C
      OPEN(94,FILE='LOADTIMESERIES.CTL',IOSTAT=IOS,STATUS='OLD')
      IF(IOS.NE.0) GOTO 500
      DO I=1,4
        READ(94,*)
      END DO
      READ(94,*)LOADTSFLAG
      IF(LOADTSFLAG.EQ.0) THEN
        CLOSE(94)
        GOTO 500
      ENDIF

      DO I=1,9
        READ(94,*)
      ENDDO
C                           SAL      TEM      DYE      SF
      READ(94,*,IOSTAT=ISO) NLSER(1),NLSER(2),NLSER(3),NLSER(4)
      IF(NLSER(1).GT.NLSERM) THEN
        WRITE(*,510)NLSER(1),NLSERM
510     FORMAT('The number of Load Salinity Time Series '             
     &  'in your application',I5,                                     
     $  '  is larger than the maximum allowed of',I5,                 
     &  '  for this compilation')                                     
        WRITE(*,*)'The program will stop check file LOADTIMESERIES.CTL'
        STOP                                                          
      END IF   
      IF(NLSER(2).GT.NLSERM) THEN
        WRITE(*,511)NLSER(2),NLSERM
511     FORMAT('The number of Load Temperature Time Series '            
     &  'in your application',I5,                                     
     $  '  is larger than the maximum allowed of',I5,                 
     &  '  for this compilation')                                     
        WRITE(*,*)'The program will stop check file LOADTIMESERIES.CTL'
        STOP                                                          
      END IF                                                          
      IF(NLSER(3).GT.NLSERM) THEN
        WRITE(*,512)NLSER(3),NLSERM
512     FORMAT('The number of Load Dye Time Series '             
     &  'in your application',I5,                                     
     $  '  is larger than the maximum allowed of',I5,                 
     &  '  for this compilation')                                     
        WRITE(*,*)'The program will stop check file LOADTIMESERIES.CTL'
        STOP                                                          
      END IF                                                          
      IF(NLSER(4).GT.NLSERM) THEN
        WRITE(*,513)NLSER(4),NLSERM
513     FORMAT('The number of Load Shellfish Larvae Time Series ' 
     &  'in your application',I5,                                     
     $  '  is larger than the maximum allowed of',I5,                 
     &  '  for this compilation')                                     
        WRITE(*,*)'The program will stop check file LOADTIMESERIES.CTL'
        STOP                                                          
      END IF                                                          
      
      DO I=1,9
        READ(94,*)
      ENDDO
      READ(94,*,IOSTAT=ISO)(NLIJ(I),I=1,4)
      DO I=1,4
        IF(NLIJ(I).GT.NLIJM) THEN                        
          WRITE(*,514)NLIJ(I),NLIJM                                    
514       FORMAT('The number of LOAD DISCHARGE Locations NLIJ(',I5,')',
     & ' is larger than the maximum allowed of',I5,
     & '  for this compilation')                                      
        WRITE(*,*)'The program will stop check file LOADTIMESERIES.CTL'
        STOP                                                           
        END IF  
      END DO

      DO NL=1,4
        DO J=1,8
          READ(94,*)
        ENDDO
        DO L=1,NLIJ(NL)
         READ(94,*,IOSTAT=ISO)ILDS(L,NL),JLDS(L,NL),NSERL(L,NL)
         IF(ISO.GT.0) THEN
          WRITE(*,5104)NLIJ(I)                                   
5104       FORMAT('The number of I,J LOAD Locations for NLIJ(',I3,')',
     & ' is different than the number indicated')                                      
         WRITE(*,*)'The program will stop check file LOADTIMESERIES.CTL'
          STOP                                                           
         END IF  
          II=ILDS(L,NL)
          JJ=JLDS(L,NL)
          LLDS(L,NL)=LIJ(II,JJ)
        END DO 
      END DO
      CLOSE(94)
      WRITE(*,*)'FILE LOADTIMESERIES.CTL WAS ACTIVE AND READ CORRECTLY'
       
C    READ LOAD TIME SERIES FILES
       
      DO NC=1,4
        IF(NC.EQ.1.AND.NLSER(NC).GE.1) THEN
          OPEN(1,FILE='LSSER.INP',STATUS='UNKNOWN')
        ELSE IF(NC.EQ.2.AND.NLSER(NC).GE.1) THEN
          OPEN(1,FILE='LTSER.INP',STATUS='UNKNOWN')
        ELSE IF(NC.EQ.3.AND.NLSER(NC).GE.1) THEN
          OPEN(1,FILE='LDSER.INP',STATUS='UNKNOWN')
        ELSE IF(NC.EQ.4.AND.NLSER(NC).GE.1) THEN
          OPEN(1,FILE='LSFSER.INP',STATUS='UNKNOWN')
        ELSE
          CYCLE                                       !skips to the next iteration of a DO loop
        END IF
        
        DO I=1,15  
         READ(1,*)
        ENDDO
C
        DO NS=1,NLSER(NC)
          READ(1,*,IOSTAT=ISO)ISTYP,MLSER(NS,NC),TCLSER(NS,NC),
     &                   TALSER(NS,NC),RMULADJ,ADDADJ
           IF(ISO.GT.0) GOTO 8700  
          IF(MLSER(NS,NC).gt.NDLSER) then                             
            WRITE(*,515)MLSER(NS,NC),NS,NC,NDLSER                        
515         FORMAT('The number of data pairs in your application',I8, 
     $      '  for LOAD time series',I6,'for variable NC=',I2,                         
     &      '  is larger than the maximum allowed of',I8,             
     &      '  for this compilation')                                 
            WRITE(*,*)'The program will stop'   
            WRITE(*,*)'(NC=1 SALT, NC=2 HEAT, NC=3 DYE, NC=4 SFL)' 
            STOP                                                      
          END IF                                                      
          IF(ISTYP.EQ.1)THEN
            READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
            IF(ISO.GT.0) GOTO 8700  
            DO M=1,MLSER(NS,NC)
              READ(1,*,IOSTAT=ISO)TLSER(M,NS,NC),SERLTMP
              IF(ISO.GT.0) GOTO 8700  
              TLSER(M,NS,NC)=TLSER(M,NS,NC)+TALSER(NS,NC)
               DO K=1,KC
                SERL(M,K,NS,NC)=(RMULADJ*(SERLTMP+ADDADJ))*WKQ(K)
              ENDDO
            ENDDO
          ELSE
            DO M=1,MLSER(NS,NC)
             READ(1,*,IOSTAT=ISO)TLSER(M,NS,NC),(SERL(M,K,NS,NC),K=1,KC)
             IF(ISO.GT.0) GOTO 8700  
             TLSER(M,NS,NC)=TLSER(M,NS,NC)+TALSER(NS,NC)
             DO K=1,KC
              SERL(M,K,NS,NC)=RMULADJ*(SERL(M,K,NS,NC)+ADDADJ)
             ENDDO
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      END DO
        
500   CONTINUE     !SKIP TO HERE IF THERE IS NO LOAD TIME SERIES CONTROL FILE OR THE USER DON'T WANT TO READ THE CTL FILE
C**************************************************************************************************C
C
C ** END READ LOAD INPUT CONTROL FILE AND LOAD TIME SERIES FOR SALT, TEMP(HEAT), DYE, SHELLFISH LARVE 
C
C*************************************************************^^^^^^*********C    END HNR_GHD ADDED ON 7/2023
C
C      
C****************************************************************C
C
C **  READ IN FREE SURFACE ELEVATION OR PRESSURE CONTROLLED FLOW
C **  SPECIFICATION FROM THE FILE QCTL.INP
C
C **  THE FLOW IS GIVEN BY:
C
C     HUP=HP(L)+BELV(L)+HCTLUA(NS)=ADJUSTED ELEVATION OF UPSTREAM CELL
C         FREE SURFACE
C
C     HDW=HP(L)+BELV(L)+HCTLDA(NS)=ADJUSTED ELEVATION OF DOWNSTREAM CELL
C         FREE SURFACE
C
C     DELH=HCTLUM(NS)*HUP-HCTLDM(NS)*HDW
C
C     IF(DELH.LE.0)THEN
C        FLOW=0
C      ELSE
C        ENTER QCTL(M,K,NS) VS HDIFCTL(M,NS) TABLE WITH DELH TO GIVE
C        FLOW=QCTL(M,K,NS)
C     ENDIF
C
C----------------------------------------------------------------------C
C
      IF(NQCTL.GE.1)THEN
        OPEN(1,FILE='QCTL.INP',STATUS='UNKNOWN')
        OPEN(99,FILE='QCTLCK.INP',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE='QCTLCK.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,14
        READ(1,1)
        ENDDO
C
        DO NS=1,NQCTLT
        READ(1,*, IOSTAT=ISO)ISTYP,MQCTL(NS),HCTLUA(NS),HCTLUM(NS),
     &              HCTLDA(NS),HCTLDM(NS),RMULADJ,ADDADJ,AQCTL(NS)
        WRITE(99,991)NS
        WRITE(99,992)ISTYP,MQCTL(NS),HCTLUA(NS),HCTLUM(NS),
     &              HCTLDA(NS),HCTLDM(NS),RMULADJ,ADDADJ,AQCTL(NS)
        IF(ISO.GT.0) GOTO 920
C
        IF(ISTYP.EQ.0)THEN
          DO M=1,MQCTL(NS)
           READ(1,*,IOSTAT=ISO) HDIFCTL(M,NS),(QCTL(M,1,K,NS),K=1,KC)
           IF(ISO.GT.0) GOTO 920
           DO K=1,KC
            QCTL(M,1,K,NS)=RMULADJ*(QCTL(M,1,K,NS)+ADDADJ)
           ENDDO
          ENDDO
        ENDIF
C
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 920
          DO M=1,MQCTL(NS)
           READ(1,*,IOSTAT=ISO) HDIFCTL(M,NS),QCTLTMP
           IF(ISO.GT.0) GOTO 920
           DO K=1,KC
            QCTL(M,1,K,NS)=RMULADJ*(QCTLTMP+ADDADJ)*WKQ(K)
           ENDDO
          ENDDO
        ENDIF
C
        IF(ISTYP.EQ.2)THEN
          DO MD=1,MQCTL(NS)
          DO MU=1,MQCTL(NS)
           READ(1,*,IOSTAT=ISO) HDIFCTL(MU,NS),HDIFCTD(MD,NS),
     &                          (QCTL(MU,MD,K,NS),K=1,KC)
           IF(ISO.GT.0) GOTO 920
           DO K=1,KC
            QCTL(MU,MD,K,NS)=RMULADJ*(QCTL(MU,MD,K,NS)+ADDADJ)
           ENDDO
          ENDDO
          ENDDO
        ENDIF
C
        IF(ISTYP.EQ.3)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 920
          DO MD=1,MQCTL(NS)
          DO MU=1,MQCTL(NS)
           READ(1,*,IOSTAT=ISO) HDIFCTL(MU,NS),HDIFCTD(MD,NS),QCTLTMP
           IF(ISO.GT.0) GOTO 920
           DO K=1,KC
            QCTL(MU,MD,K,NS)=RMULADJ*(QCTLTMP+ADDADJ)*WKQ(K)
           ENDDO
          ENDDO
          ENDDO
        ENDIF
C
        IF(ISTYP.LE.1)THEN
         DO M=1,MQCTL(NS)
          WRITE(99,993)M,HDIFCTL(M,NS),(QCTL(M,1,K,NS),K=1,KC)
         ENDDO
        ENDIF
        IF(ISTYP.GE.2)THEN
         DO MD=1,MQCTL(NS)
         DO MU=1,MQCTL(NS)
          WRITE(99,994)MU,MD,HDIFCTL(MU,NS),HDIFCTD(MD,NS),
     &                   (QCTL(MU,MD,K,NS),K=1,KC)
         ENDDO
         ENDDO
        ENDIF
C
        ENDDO
C
        CLOSE(1)
        CLOSE(99)
      ENDIF
C
  991 FORMAT(/,' CONTROL TABLE NS =',I5,/)
  992 FORMAT(2I5,7F10.4)
  993 FORMAT(I5,11F10.4)
  994 FORMAT(2I5,11F10.4)
 1001 FORMAT(/,' READ ERROR FROM FILE EFDC.INP ON CARD ',I5/)
 1002 FORMAT(/,' INPUT ECHO NCARD = ',I5/)
C
C**********************************************************************C
C
C **  READ IN WITHDRAWAL, ADD HEAT OR MATERIAL, RETURN TIME SERIES
C **  FROM THE FILE QWRSER.INP
C
C----------------------------------------------------------------------C
C
C     IF(NQCOOL.GE.1)THEN
C       OPEN(1,FILE='QCOOL.INP',STATUS='UNKNOWN')
C       DO NS=1,NQCOOL
C       READ(1,*,ERR=930)ISTYP,MQLSER(NS),TCQLSER(NS),TAQLSER(NS),
C    &                   RMULADJ,ADDADJ
C       IF(ISTYP.EQ.1)THEN
C         READ(1,*,ERR=930) (WKQ(K),K=1,KC)
C          DO M=1,MQLSER(NS)
C          READ(1,*,ERR=930)TQLSER(M,NS),QLSERTMP,TRISETMP
C          TQLSER(M,NS)=TQLSER(M,NS)+TAQLSER(NS)
C           DO K=1,KC
C           QLSER(M,K,NS)=RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
C            DO NC=1,5
C            CQLSER(M,K,NS,NC)=0.
C            ENDDO
C           CQLSER(M,K,NS,2)=TRISETMP
C           ENDDO
C          ENDDO
C        ELSE
C         DO M=1,MCSER(NS,NC)
C         READ(1,*,ERR=910)TQLSER(M,NS),TRISETMP,(QLSER(M,K,NS),K=1,KC)
C         TQLSER(M,NS)=TQLSER(M,NS)+TAQLSER(NS)
C          DO K=1,KC
C          QLSER(M,K,NS)=RMULADJ*(QLSER(M,K,NS)+ADDADJ)
C            DO NC=1,5
C            CQLSER(M,K,NS,NC)=0.
C            ENDDO
C          CQLSER(M,K,NS,2)=TRISETMP
C          ENDDO
C         ENDDO
C       ENDIF
C       ENDDO
C       CLOSE(1)
C     ENDIF
C
C**********************************************************************C
C
C **  READ IN ATMOSPHERIC AND WEATHER CONDITION TIME SERIES FROM THE
C **  FILE ASER.INP
C
C----------------------------------------------------------------------C
C
C     INITIALIZE ATMOSPHERIC VARIABLES
C
      DO L=2,LA
       PATMT(L)=1000.
       TATMT(L)=0.
       RAINT(L)=0.
       EVAPT(L)=0.
       SOLSWRT(L)=0.
       CLOUDT(L)=0.
       SVPA(L)=0.
       RHA(L)=0.
       VPA(L)=0.
       CLEVAP(L)=0.
       CCNHTT(L)=0.
      ENDDO
C
      TBEDITTMP=TEMO
C
      IF(NASER.GT.0)THEN
      OPEN(1,FILE='ASER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,36
       READ(1,1)
      ENDDO
C
      DO N=1,NASER
C
      READ(1,*,IOSTAT=ISO) MASER(N),TCASER(N),TAASER(N),IRELH(N),
     &                RAINCVT,EVAPCVT,SOLRCVT,CLDCVT
c      write(*,*) MASER(N),TCASER(N),TAASER(N),IRELH(N),
c     &                RAINCVT,EVAPCVT,SOLRCVT,CLDCVT      
      IF(ISO.GT.0) GOTO 940
      if(maser(n).gt.ndaser) then                                      !hnr
        write(*,380)maser(n),n,ndaser                                  !hnr
380     FORMAT('The number of data pairs in your application',I7,      !hnr
     $  '  for Aser time series',I3,                                   !hnr
     &  '  is larger than the maximum allowed of',I7,                  !hnr
     &  '  for this compilation')                                      !hnr
        write(*,*)'The program will stop aser file'                    !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      READ(1,*,IOSTAT=ISO) IASWRAD,REVC,RCHC,SWRATNF,
     &  SWRATNS,FSWRATF,DABEDT,TBEDITTMP,HTBED1,HTBED2
c      write(*,*) IASWRAD,REVC,RCHC,SWRATNF,
c     &  SWRATNS,FSWRATF,DABEDT,TBEDITTMP,HTBED1,HTBED2      
      IF(ISO.GT.0) GOTO 940
C
      DO M=1,MASER(N)
      READ(1,*,IOSTAT=ISO)TASER(M,N),PATM(M,N),TDRY(M,N),
     &       TWET(M,N),RAIN(M,N),EVAP(M,N),SOLSWR(M,N),CLOUD(M,N)
c      WRITE(*,*)TASER(M,N),PATM(M,N),TDRY(M,N),
c     &       TWET(M,N),RAIN(M,N),EVAP(M,N),SOLSWR(M,N),CLOUD(M,N)      
      IF(ISO.GT.0) GOTO 940
      ENDDO
C
      DO M=1,MASER(N)
      TASER(M,N)=TASER(M,N)+TAASER(N)
      RAIN(M,N)=RAINCVT*RAIN(M,N)
      EVAP(M,N)=EVAPCVT*EVAP(M,N)
      SOLSWR(M,N)=SOLRCVT*SOLSWR(M,N)
      CLOUD(M,N)=CLDCVT*CLOUD(M,N)
      ENDDO
C
      ENDDO
C
      CLOSE(1)
      ENDIF
C
      IF(NASER.LE.1)THEN
        DO L=2,LA
         ATMWHT(L,1)=1.
        ENDDO
      ENDIF
C
      IF(NASER.GT.1)THEN
        OPEN(1,FILE='ATMMAP.INP',STATUS='UNKNOWN')
        DO IS=1,4
          READ(1,1)
        ENDDO
        DO L=2,LA
          READ(1,*)LD,ID,JD,(ATMWHT(L,N),N=1,NASER)
        ENDDO
        CLOSE(1)
      ENDIF
C
      IF(TBEDITTMP.GE.0.0)THEN
	DO L=1,LC
	  TBEDIT(L)=TBEDITTMP
	ENDDO
	ENDIF
C
      IF(TBEDITTMP.LT.0.0)THEN
        OPEN(1,FILE='TEMPBBOT.INP',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,4
          READ(1,1)
        ENDDO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0)THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) TBEDIT(L)
            IF(ISO.GT.0) GOTO 842
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,TBEDIT(L)
            IF(ISO.GT.0) GOTO 842
          ENDDO
        ENDIF
        CLOSE(1)
        TBEDIT(1)=TBEDIT(1)
        TBEDIT(LC)=TBEDIT(LC-1)
	ENDIF
C
      DO K=1,KBHM
      TEMB(1,K)=TBEDIT(1)
      TEMB(LC,K)=TBEDIT(LC)
      TEMB1(1,K)=TBEDIT(1)
      TEMB1(LC,K)=TBEDIT(LC)
	ENDDO
C
C      IF(ISRESTI.EQ.0)THEN
      IF(ISRESTI.EQ.0.AND.ISBEDTEMI.EQ.1)THEN
      DO K=1,KBHM
        DO L=2,LA
         TEMB(L,K)=TBEDIT(L)
         TEMB1(L,K)=TBEDIT(L)
        ENDDO
        ENDDO
      ENDIF
C
CDIA      WRITE(6,6776)'READ ASER.INP'
C
C**********************************************************************C
C
C **  READ IN ABOVE WATER SURFACE WIND TIME SERIES FROM THE
C **  FILE WSER.INP OF TIME VARYING WIND FORCING FUNCTION FROM
C **  FILE WFSER.INP
C
C----------------------------------------------------------------------C
C
C     INITIALIZE WIND STRESS AND SPEED
C
      DO L=2,LA
        WINDST(L)=0.
        TSX(L)=0.
        TSY(L)=0.
      ENDDO
C
      IF(NWSER.GT.0)THEN
      OPEN(1,FILE='WSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
C      DO IS=1,19
      DO IS=1,17
      READ(1,1)
      ENDDO
C
      DO N=1,NWSER
C
      READ(1,*,IOSTAT=ISO) MWSER(N),TCWSER(N),TAWSER(N),WINDSCT,
     &                      ISWDINT(N)
C     &                      ISWDINT(N),WSHELE,WSHELN
      IF(ISO.GT.0) GOTO 940
C
      if(mwser(n).gt.ndwser) then                                      !hnr
        write(*,490)mwser(n),n,ndwser                                  !hnr
490     FORMAT('The number of data pairs in your application',I7,      !hnr
     $  ' for Wser time series',I3,                                    !hnr
     &  ' is larger than the maximum allowed of',I7,                   !hnr
     &  ' for this compilation')                                       !hnr
        write(*,*)'The program will stop wser file'                              !hnr
        stop                                                           !hnr
      end if                                                           !hnr
      DO M=1,MWSER(N)
      READ(1,*,IOSTAT=ISO)TWSER(M,N),WINDS(M,N),WINDD(M,N)
      IF(ISO.GT.0) GOTO 940
      ENDDO
C
      DO M=1,MWSER(N)
       TWSER(M,N)=TWSER(M,N)+TAWSER(N)
      ENDDO
C
      IF(ISWDINT(N).LE.1)THEN
      DO M=1,MWSER(N)
       WINDS(M,N)=WINDSCT*WINDS(M,N)
      ENDDO
      ENDIF
C
      IF(ISWDINT(N).EQ.1)THEN
      DO M=1,MWSER(N)
        IF(WINDD(M,N).LE.180.)THEN
          WINDD(M,N)=WINDD(M,N)+180.
          IF(WINDD(M,N).EQ.360.) WINDD(M,N)=0.
         ELSE
          WINDD(M,N)=WINDD(M,N)-180.
          IF(WINDD(M,N).EQ.360.) WINDD(M,N)=0.
        ENDIF
      ENDDO
      ENDIF
C
      IF(ISWDINT(N).EQ.2)THEN
      DO M=1,MWSER(N)
       WINDS(M,N)=WINDSCT*WINDS(M,N)
       WINDD(M,N)=WINDSCT*WINDD(M,N)
      ENDDO
      ENDIF
C
      ENDDO
C
      CLOSE(1)
      ENDIF
C
      IF(NWSER.EQ.1)THEN
        DO L=2,LA
         WNDWHT(L,1,1)=1.
        ENDDO
      ENDIF
C
      IF(NWSER.GT.1)THEN
        OPEN(1,FILE='WNDMAP.INP',STATUS='UNKNOWN')
        DO IS=1,6
          READ(1,1)
        ENDDO
        READ(1,*)NWNDMAP
	  DO NW=1,NWNDMAP
	    READ(1,*)TWNDMAPBEG(NW),TWNDMAPEND(NW)
          DO L=2,LA
              READ(1,*)LD,ID,JD,(WNDWHT(L,N,NW),N=1,NWSER)
          ENDDO
        ENDDO
        CLOSE(1)
      ENDIF
C
C ** READ WIND FIELD FORCING FUNCTION
C
      IF(NWSER.LE.-1) THEN
C
        OPEN(1,FILE='WFSER.INP')
C
C **    SIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,29
          READ(1,1)
        ENDDO
C
        READ(1,*,IOSTAT=ISO) MWSERTMP,MWSER(1),TCWSER(1),TAWSER(1),
     &          RMULADJ,MCOEFWF,XWFMUL,XWFADD,YWFMUL,YWFADD
        IF(ISO.GT.0) GOTO 940
C
        DO M=1,MWSERTMP
	    READ(1,*)XWFKER(M),YWFKER(M)
        ENDDO
        DO M=1,MWSERTMP
	    XWFKER(M)=XWFMUL*(XWFKER(M)+XWFADD)
		YWFKER(M)=YWFMUL*(YWFKER(M)+YWFADD)
        ENDDO
C
        DO M=1,MWSER(1)
          READ(1,*,IOSTAT=ISO)TWSER(M,1),(WFCOEF(M,K),K=1,MCOEFWF)
	    DO K=1,MCOEFWF
	      WFCOEF(M,K)=RMULADJ*WFCOEF(M,K)
          ENDDO
          IF(ISO.GT.0) GOTO 940
        ENDDO
C
        DO M=1,MWSER(1)
          TWSER(M,1)=TWSER(M,1)+TAWSER(1)
        ENDDO
C
        CLOSE(1)
C
      ENDIF
C
C      IF(NWSER.EQ.-2) THEN
C
C        OPEN(1,FILE='WFSER.OUT')
C        DO M=1,MWSER(1)
C          WRITE(1,1981)TWSER(M,1),(WFCOEF(M,K),K=1,MCOEFWF)
C        ENDDO
C        CLOSE(1)
C
C      ENDIF
C
 1981 FORMAT(F12.4,12E13.5)
C
CDIA      WRITE(6,6776)'READ WSER.INP'
C
C**********************************************************************C
C
C **  READ IN EXTERNALLY SPECIFIED ICE COVER INFORMATION
C **  FROM THE FILE ICECOVER.INP
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
        RICECOVL(L)=0.0
	  RICETHKL(L)=0.0
	  RICEWHT(L,1)=1.0
      ENDDO
C
      IF(ISICE.EQ.1.AND.NISER.GE.1)THEN
      OPEN(1,FILE='ICECOVER.INP')
      DO IS=1,14
	  READ(1,1)
	ENDDO
C
      DO N=1,NISER
	  RICECOVT(N)=0.
	  RICETHKT(N)=0.
        MITLAST(N)=1
C
	  READ(1,*,IOSTAT=ISO)MISER(N),TCISER(N),TAISER(N),
     &                      RMULADJCOV,RMULADJTHK
        IF(ISO.GT.0) GOTO 950
C
	  DO M=1,MISER(N)
	    READ(1,*,IOSTAT=ISO)TISER(M,N),RICECOVS(M,N),RICETHKS(M,N)
          IF(ISO.GT.0) GOTO 950
        ENDDO
C
	  DO M=1,MISER(N)
          TISER(M,N)=TISER(M,N)+TAISER(N)
		RICECOVS(M,N)=RMULADJCOV*RICECOVS(M,N)
		RICETHKS(M,N)=RMULADJTHK*RICETHKS(M,N)
        ENDDO
C
      ENDDO
C
      CLOSE(1)
	ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISICE.EQ.1.AND.NISER.GT.1)THEN
        OPEN(1,FILE='ICEMAP.INP')
        DO IS=1,4
          READ(1,1)
        ENDDO
        DO L=2,LA
          READ(1,*)LD,ID,JD,(RICEWHT(L,N),N=1,NISER)
        ENDDO
        CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  READ IN SHELL FISH LARAVE BEHAVIOR DATA
C **  FROM THE FILE SFBSER.INP
C
C----------------------------------------------------------------------C
C
      IF(ISTRAN(4).GE.1)THEN
        OPEN(1,FILE='SFBSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,18
      READ(1,1)
      ENDDO
C
      READ(1,*,IOSTAT=ISO) MSFSER,TCSFSER,TASFSER,TSRSF,TSSSF,
     &                     ISSFLDN,ISSFLFE,SFLKILL
      IF(ISO.GT.0) GOTO 960
C
      DO M=1,MSFSER
      READ(1,*,IOSTAT=ISO) TSFSER(M),RKDSFL(M),WSFLST(M),WSFLSM(M),
     &                     DSFLMN(M),DSFLMX(M),SFNTBE(M),SFATBT(M)
      IF(ISO.GT.0) GOTO 960
      ENDDO
C
      CLOSE(1)
C
      ENDIF
C
C**********************************************************************C
C
C **  READ VEGETATION DATA FROM VEGE.INP AND VEGSER.INP
C
      IF(ISVEG.GE.1)THEN
C
      OPEN(1,FILE='VEGE.INP',STATUS='UNKNOWN')
C
      DO NS=1,12
      READ(1,1)
      ENDDO
C
      READ(1,*)MVEGTYP,MVEGOW,NVEGSER,UVEGSCL
      DO M=1,MVEGTYP
      READ(1,*,ERR=3120)IDUM,NVEGSERV(M),RDLPSQ(M),BPVEG(M),HPVEG(M),
     &         ALPVEG(M),BETVEG(M),GAMVEG(M),SCVEG(M)
      BDLTMP=BPVEG(M)*BPVEG(M)*RDLPSQ(M)
      PVEGX(M)=1.-BETVEG(M)*BDLTMP
      PVEGY(M)=1.-BETVEG(M)*BDLTMP
      PVEGZ(M)=1.-ALPVEG(M)*BDLTMP
      BDLPSQ(M)=BPVEG(M)*RDLPSQ(M)
      ENDDO
C
      CLOSE(1)
C
      IF(NVEGSER.GT.0)THEN
	  DO M=1,NVEGSER
          MVEGTLAST(M)=1
	  ENDDO
	ENDIF
C
      ENDIF
C
      GOTO 3122
 3120 WRITE(6,3121)
 3121 FORMAT('  READ ERROR FOR FILE VEGE.INP ')
      STOP
 3122 CONTINUE

C
      IF(NVEGSER.GE.1)THEN
        OPEN(1,FILE='VEGSER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,8
        READ(1,1)
        ENDDO
C
        DO NS=1,NVEGSER
          READ(1,*,IOSTAT=ISO) MVEGSER(NS),TCVEGSER(NS),TAVEGSER(NS)
          IF(ISO.GT.0) GOTO 7120
          DO M=1,MVEGSER(NS)
            READ(1,*,IOSTAT=ISO)TVEGSER(M,NS),VEGSERR(M,NS),
     &      VEGSERB(M,NS),VEGSERH(M,NS)
            IF(ISO.GT.0) GOTO 7120
            TVEGSER(M,NS)=TVEGSER(M,NS)+TAVEGSER(NS)
          ENDDO
        ENDDO
        CLOSE(1)
C
C **  REINITIALIZE CLASSES HAVING TIME SERIES INFORMATION
C
      DO M=1,MVEGTYP
        IF(NVEGSERV(M).GT.0)THEN
	    NS=NVEGSERV(M)
	    RDLPSQ(M)=VEGSERR(1,NS)
		BPVEG(M)=VEGSERB(1,NS)
		HPVEG(M)=VEGSERH(1,NS)
          BDLTMP=BPVEG(M)*BPVEG(M)*RDLPSQ(M)
          PVEGX(M)=1.-BETVEG(M)*BDLTMP
          PVEGY(M)=1.-BETVEG(M)*BDLTMP
          PVEGZ(M)=1.-ALPVEG(M)*BDLTMP
          BDLPSQ(M)=BPVEG(M)*RDLPSQ(M)
	  ENDIF
      ENDDO
C
      ENDIF
C
      GOTO 7122
 7120 WRITE(6,7121)
 7121 FORMAT('  READ ERROR FOR FILE VEGSER.INP ')
      STOP
 7122 CONTINUE
C
C**********************************************************************C
C
C  READ MASS BALANCE ADJUSTMENT FOR A LAKE PARAMETERS       HNR
C
      islake=0                                                         !hnr
      OPEN(1,FILE='efdclakecorrection.inp',iostat=ios,STATUS='old')    !hnr
      IF(ios.eq.0) then                                                !hnr
C                                                                      !hnr
        DO NSKIP=1,13                                                  !hnr
          READ(1,1)                                                    !hnr
        END DO                                                         !hnr
        READ(1,*,IOSTAT=ISO)ISLAKE,ILK,JLK,PSERLK,ILKTS,NWLKTS,        !hnr
     &           SMTHF,NCLK                                            !hnr
        IF(ISO.GT.0) then                                              !hnr
          write(6,*)'error in efdclakecorrection.inp'                  !hnr
          stop                                                         !hnr
        end if                                                         !hnr
      end if                                                           !hnr
      JSLKTS=1                                                         !hnr
      NCLKTS=1                                                         !hnr
      NCL=1                                                            !hnr
C
C
C**********************************************************************C
C
C **  READ BANK EROSION MAP AND TIME SERIES FILE
C
      IF(ISBKERO.GE.1)THEN
C
        OPEN(1,FILE='BEMAP.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,10
        READ(1,1)
        ENDDO
C
        READ(1,*)NBEPAIR,NBESER
C
        DO NS=1,NBEPAIR
          READ(1,*)IBANKBE(NS),JBANKBE(NS),ICHANBE(NS),JCHANBE(NS),
     &             NBESERN(NS),FBESER(NS)
        ENDDO
C
        CLOSE(1)
C
      ENDIF
C
C
      IF(ISBKERO.GE.1.AND.NBESER.GT.0)THEN
C
        OPEN(1,FILE='BESER.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,10
        READ(1,1)
        ENDDO
C
        DO NS=1,NBESER
C
          READ(1,*)MBESER(NS),TCBESER(NS),TABESER(NS),RMULADJ,ADDADJ
          MBETLAST(NS)=1
C
          DO M=1,MBESER(NS)
            READ(1,*)TBESER(M,NS),BESER(M,NS),FWCBESER(M,NS)
            TBESER(M,NS)=TBESER(M,NS)+TABESER(NS)
            BESER(M,NS)=RMULADJ*(BESER(M,NS)+ADDADJ)
          ENDDO
C
        ENDDO
C
        CLOSE(1)
C
      ENDIF
C
C**********************************************************************C
C
C **  READ ZONALLY VARYING SEDIMENT BED PARTICLE MIXING
C
C----------------------------------------------------------------------C
C
      ITMPPMX=0
      DO NT=1,NTOX
        IF(ISPMXZ(NT).EQ.1)ITMPPMX=1
      ENDDO
C
      IF(ITMPPMX.EQ.1)THEN
C
        OPEN(1,FILE='PARTMIX.INP')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,7
        READ(1,1)
        ENDDO
C
C#######################################################################
C     HQI change to input multiplication scale factor for particle mixing rate
C     RM 10/06/05
c        READ(1,*)NPMXZ,NPMXPTS
         READ(1,*)NPMXZ,NPMXPTS,PMIXSF
C#######################################################################
        DO NZ=1,NPMXZ
          DO NP=1,NPMXPTS
              READ(1,*)PMXDEPTH(NP,NZ),PMXCOEF(NP,NZ)
C#######################################################################
C     HQI change to input multiplication scale factor for particle mixing rate
C     RM 10/06/05
              PMXCOEF(NP,NZ) = PMIXSF*PMXCOEF(NP,NZ)
C#######################################################################
          ENDDO
        ENDDO
C
        CLOSE(1)
C
C
        OPEN(1,FILE='PMXMAP.INP')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,4
        READ(1,1)
        ENDDO
C
        DO L=2,LC-1
           READ(1,*)LDUM,IDUM,JDUM,LPMXZ(L)
        ENDDO
C
        CLOSE(1)
C
      ENDIF
C
C**********************************************************************C
C
C **  READ IN HARMONIC TIDE GAUGE DATA FOR TIDE ELEVAVTION
C **  ASSIMULATION FROM FILE TIDASM.INP
C
C----------------------------------------------------------------------C
C
C      IF(ITIDASM.GE.1)THEN
C        OPEN(1,FILE='TIDASM.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
C      DO IS=1,25
C      READ(1,1)
C      ENDDO
C
C      READ(1,*,IOSTAT=ISO) MLTIDAMS, MTTIDASM
C      IF(ISO.GT.0) GOTO 960
C
C      DO ML=1,MLTIDAMS
C      READ(1,*,IOSTAT=ISO) IITIDASM(ML),JJTIDASM(ML),TTIDASM(ML)
C      IF(ISO.GT.0) GOTO 970
C        DO MT=1,MTIDE
C          READ(1,*,IOSTAT=ISO)PFAM(ML,MT), PFPH(ML,MT)
C          IF(ISO.GT.0) GOTO 970
C        ENDDO
C      ENDDO
C
C      CLOSE(1)
C
C      ENDIF
C
C**********************************************************************C
C
      GOTO 3000
C
C**********************************************************************C
C
C **  WRITE READ ERROR FOR OTHER INPUT FILES AND TERMINATE RUN
C
  800 WRITE(6,801)
!      WRITE(8,801)  !hnr 7/27/2009
  801 FORMAT('  READ ERROR FOR FILE CELL.INP ')
      STOP
  810 WRITE(6,811)
!      WRITE(8,811)    !hnr 7/27/2009
  811 FORMAT('  READ ERROR FOR FILE CELLLT.INP ')
      STOP
  820 WRITE(6,821)
!      WRITE(8,821)   !hnr 7/27/2009
  821 FORMAT('  READ ERROR FOR FILE DEPTH.INP ')
      STOP
  830 WRITE(6,831)
!      WRITE(8,831)   !hnr 7/27/2009
  831 FORMAT('  READ ERROR FOR FILE DXDY.INP ')
      STOP
  840 WRITE(6,841)
!      WRITE(8,841)    !hnr 7/27/2009
  841 FORMAT('  READ ERROR FOR FILE SALT.INP ')
      STOP
  842 WRITE(6,843)
!      WRITE(8,843)  !hnr 7/27/2009
  843 FORMAT('  READ ERROR FOR FILE TEMP.INP ')
      STOP
  844 WRITE(6,845)
!      WRITE(8,845)   !hnr 7/27/2009
  845 FORMAT('  READ ERROR FOR FILE DYE.INP ')
      STOP
  846 WRITE(6,847)
!      WRITE(8,847)    !hnr 7/27/2009
  847 FORMAT('  READ ERROR FOR FILE SFL.INP ')
      STOP
  848 WRITE(6,849)
!      WRITE(8,849)    !hnr 7/27/2009
  849 FORMAT('  READ ERROR FOR FILE TOXW.INP ')
      STOP
  852 WRITE(6,853)
!      WRITE(8,853)    !hnr 7/27/2009
  853 FORMAT('  READ ERROR FOR FILE TOXB.INP ')
      STOP
  850 WRITE(6,851)
!      WRITE(8,851)   !hnr 7/27/2009
  851 FORMAT('  READ ERROR FOR FILE PSER.INP ')
      STOP
  854 WRITE(6,855)
!      WRITE(8,855)    !hnr 7/27/2009
  855 FORMAT('  READ ERROR FOR FILE SEDW.INP ')
      STOP
  856 WRITE(6,857)
!      WRITE(8,857)    !hnr 7/27/2009
  857 FORMAT('  READ ERROR FOR FILE SEDB.INP ')
      STOP
  858 WRITE(6,859)
!      WRITE(8,859)    !hnr 7/27/2009
  859 FORMAT('  READ ERROR FOR FILE SNDW.INP ')
      STOP
  862 WRITE(6,863)
!      WRITE(8,863)   !hnr 7/27/2009
  863 FORMAT('  READ ERROR FOR FILE SNDB.INP ')
      STOP
  860 WRITE(6,861)
!      WRITE(8,861)   !hnr 7/27/2009
  861 FORMAT('  READ ERROR FOR FILE QSER.INP ')
      STOP
  865 WRITE(6,866)
!      WRITE(8,866)   !hnr 7/27/2009
  866 FORMAT('  READ ERROR FOR FILE QWRS.INP ')
      STOP
  870 WRITE(6,871)
!      WRITE(8,871)   !hnr 7/27/2009
  871 FORMAT('  READ ERROR FOR FILE SSER.INP ')
      STOP
  880 WRITE(6,881)
!      WRITE(8,881)   !hnr 7/27/2009
  881 FORMAT('  READ ERROR FOR FILE TSER.INP ')
      STOP
  890 WRITE(6,891)
!      WRITE(8,891)   !hnr 7/27/2009
  891 FORMAT('  READ ERROR FOR FILE DSER.INP ')
      STOP
  900 WRITE(6,901)
!      WRITE(8,901)   !hnr 7/27/2009
  901 FORMAT('  READ ERROR FOR FILE SDSER.INP ')
      STOP
  902 WRITE(6,903)
!      WRITE(8,903)   !hnr 7/27/2009
  903 FORMAT('  READ ERROR FOR FILE SNSER.INP ')
      STOP
  904 WRITE(6,905)
!      WRITE(8,905)   !hnr 7/27/2009
  905 FORMAT('  READ ERROR FOR FILE TXSER.INP ')
      STOP
  910 WRITE(6,911)
!      WRITE(8,911)   !hnr 7/27/2009
  911 FORMAT('  READ ERROR FOR FILE SFSER.INP ')
      STOP
  920 WRITE(6,921)
!      WRITE(8,921)     !hnr 7/27/2009
  921 FORMAT('  READ ERROR FOR FILE QCTL.INP ')
      STOP
C 930 WRITE(6,931)
C 931 FORMAT('  READ ERROR FOR FILE QCOOL.INP ')
C     STOP
  940 WRITE(6,941)
!      WRITE(8,941)   !hnr 7/27/2009
  941 FORMAT('  READ ERROR FOR FILE ASER.INP ')
      STOP
  950 WRITE(6,951)
!      WRITE(8,951)  !hnr 7/27/2009
  951 FORMAT('  READ ERROR FOR FILE MAPPGNS.INP ')
      STOP
  960 WRITE(6,961)
!      WRITE(8,961)  !hnr 7/27/2009
  961 FORMAT('  READ ERROR FOR FILE SFBSER.INP ')
      STOP
  970 WRITE(6,971)
!      WRITE(8,971)   !hnr 7/27/2009
  971 FORMAT('  READ ERROR FOR FILE TIDASM.INP ')
      STOP

C*************************************************************^^^^^^*********C    START HNR_GHD ADDED ON 7/2023
8700  WRITE(*,*)'READ ERROR FOR LOAD FILE LxSER.INP FOR VARIABLE NC=',NC
      WRITE(*,*)'(NC=1 SALT, NC=2 HEAT, NC=3 DYE, NC=4 SFL)' 
      STOP
C*************************************************************^^^^^^*********C    END HNR_GHD ADDED ON 7/2023

      
C
 3000 CONTINUE
C
C**********************************************************************C
C
      RETURN
      END
