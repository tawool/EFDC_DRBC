!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
! **  FILE FOR EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!     NON-CRAY VERSION
!
!**********************************************************************C
!
      PROGRAM AAEFDC
!
!**********************************************************************C
!
! **  WELCOME TO THE ENVIRONMENTAL FLUID DYNAMICS COMPUTER CODE SERIES
! **  DEVELOPED BY JOHN M. HAMRICK.  THE EFDC CODE WAS ORGINALLY
! **  DEVELOPED AT VIRGINIA INSTITUTE OF MARINE SCIENCE
! **  /SCHOOL OF MARINE SCIENCE, THE COLLEGE OF
! **  WILLIAM AND MARY, GLOUCESTER POINT, VA 23062
!
! **  THIS SOURCE FILE IS A DIRECT RELEASE BY THE DEVELOPER
! **  AND DIFFERS SIGNIFICANTLY FROM PRE 1 MARCH 96 VIMS RELEASES OF
! **  EFDC AND POST 1 MARCH 96 VIMS RELEASES OF HEM3D (THE NAMED
! **  CURRENTLY USED BY VIMS FOR THE VERSION OF EFDC EXISTING AT
! **  THE TIME OF MY DEPARTURE) WITH RESPECT TO ERROR FIXES AND
! **  APPLICATION CAPABILITIES
!
! **  EFDC IS CURRENTLY MAINTAINED BY TETRA TECH, INC., WITH PRIMARY
! **  SUPPORT FORM THE US ENVIRONMENTAL PROTECTION AGENCY
!
! **  ENVIRONMENTAL FLUID DYNAMICS CODE AND EFDC ARE
! **  TRADEMARKS OF JOHN M. HAMRICK, PH.D., P.E.
!
! **  EFDC SOLVES THE 3D REYNOLDS AVERAGED NAVIER-STOKES
! **  EQUATIONS (WITH HYDROSTATIC AND BOUSINESSQ APPROXIMATIONS) AND
! **  TRANSPORT EQUATIONS FOR TURBULENT INTENSITY, TURBULENT
! **  INTENSITYXLENGHT SCALE, SALINITY (OR WATER VAPOR CONTENT),
! **  TEMPERATURE, AN INERT TRACER (CALLED DYE), A DYNAMICALLY ACTIVE
! **  SUSPENDED SETTLING PARTICLE FIELD (CALLED SEDIMENT).  A FREE
! **  SURFACE OR RIGID LID IS PRESENT ON THE VERTICAL BOUNDARY Z=1
! **  IN THE SIGMA STRETCHED VERTICAL COORDINATE.  THE HORIZONTAL
! **  COORDINATE SYSTEM IS CURVILINEAR AND ORTHOGONAL.
!
! **  THE NUMERICAL SOLUTION SCHEME IS ON A SPATIALLY STAGGERED MAC
! **  OR C GRID AND THE TIME INTEGRATION USES A THREE TIME LEVEL
! **  LEAPFROG INTEGRATION WITH PERIODIC TRAPEZOIDAL CORRECTIONS TO
! **  SUPPRESS THE COMPUTATIONAL MODE AND REDUCE NOISE.
! **  SPATIAL SOLUTION OF THE EXTERNAL MODE FOR THE FREE SURFACE
! **  ELEVATION OR KINEMATIC PRESSURE UNDER THE RIGID LID IS BY
! **  RED-BLACK SUCCESSIVE OVER RELAXATION (RB SOR) OR CONJUGATE
! **  GRADIENT SOLUTION OF A PSEUDO-HEMHOLTZ EQUATION.  THE INTERNAL
! **  SOLUTION IS IMPLICIT FOR THE VERTICAL SHEAR OR VELOCITY STRUCTURE.
! **  A NUMBER OF OPTIONS ARE AVAILABLE FOR REPRESENTING THE ADVECTIVE
! **  TRANSPORT TERMS IN THE MOMENTUM AND SCALAR TRANSPORT EQUATIONS.
!
! **  PRIMARY DOCUMENTATION INCLUDES:
!
!     HAMRICK, J. M., 1992A:  A THREE-DIMENSIONAL ENVIRONMENTAL
!     FLUID DYNAMICS COMPUTER CODE: THEORETICAL AND COMPUTATIONAL
!     ASPECTS. THE COLLEGE OF WILLIAM AND MARY, VIRGINIA INSTITUTE
!     OF MARINE SCIENCE, SPECIAL REPORT 317, 63 PP.
!
!     HAMRICK, J. M., 1996A:  USERS MANUAL FOR THE ENVIRONMENTAL
!     FLUID DYNAMIC COMPUTER CODE. THE COLLEGE OF WILLIAM AND MARY,
!     VIRGINIA INSTITUTE OF MARINE SCIENCE, SPECIAL REPORT 328, 224 PP.
!
!     PARK, K., A. Y. KUO, J. SHEN, AND J. M. HAMRICK, 1995:
!     A THREE-DIMENSIONAL HYDRODYNAMIC-EUTROPHICATION MODEL (HEM3D):
!     DESCRIPTION OF WATER QUALITY AND SEDIMENT PROCESSES SUBMODELS.
!     THE COLLEGE OF WILLIAM AND MARY, VIRGINIA INSTITUTE OF MARINE
!     SCIENCE. SPECIAL REPORT 327, 113 PP.
!
!     TETRA TECH, INC., 1999B: THEORETICAL AND COMPUTATIONAL ASPECTS
!     OF SEDIMENT AND CONTAMINANT TRANSPORT IN THE EFDC MODEL.
!     A REPORT TO THE U.S. ENVIRONMENTAL PROTECTION AGENCY,
!     TETRA TECH, INC., FAIRFAX, VA.
!
! **  ADDITIONAL REFERENCES TO MODEL APPLICATIONS ARE AVAILABLE
! **  FROM THE DEVELOPER
!
! **  CHANGES MADE TO THIS CODE BY UNAUTHORIZED PERSONS WILL BE
! **  SUPPORTED ON A COST REIMBURSED BASIS ONLY.  SUPPORT IS AVAILABLE
! **  FROM JOHN M. HAMRICK, 2520 WEST WHITTAKER CLOSE
! **  WILLIAMSBURG, VA, TEL. 804-258-0608, FAX 804-258-9698
! **  EMAIL: HAM@VISI.NET
!
! **  THE AUTHOR AND TETRA TECH, INC.  ASSUME NO LIABILITY FOR USE
! **  OF THIS CODE FOR ENVIRONMENTAL AND ENGINEERING STUDIES.
!
! **  THE FOLLOWING FILES ARE NECESSARY TO COMPILE THIS CODE:
!
! **      EFDC.CMN
! **      EFDC.PAR
!
! **  THIS CODE HAS BEEN COMPLIED ON SUN SPARC, HP 9000/700, AND SGI
! **  WORKSTATIONS, VAX VMS AND DEC ALPHA SYSTEMS,
! **  CRAY Y/MP AND C90 SYSTEMS, MACINTOSH SYSTEMS USING LSI AND
! **  ABSOFT FORTRAN, AND 486 AND PENTIUM PC'S USING LAHEY
! **  PROFESSIONAL FORTRAN AND MICROSOFT POWERSTATION
! **  LINES IN THE CODE BEGINNING WITH CDHP IMPLEMENT THE
! **  HP 9000/700 SERIES VECTOR LIBRARY AND MAY BE ACTIVATES BY
! **  REPLACING CDHP WITH 4 BLANK SPACES.  VAX EXTENSION TIMING
! **  UTILITIES USING THE CALL SECNDS MAY BE ACTIVATED BY UNCOMMENTING
! **  ALL LINES CONTAINING SECNDS.
! **  TO RUN ON CRAY SYSTEMS, REPLACE ACCESS='APPEND' WITH
! **  POSITION='APPEND' IN APPROPRIATE FILE OPEN STATEMENTS
!
! **  THE FOLLOWING FILES MAY BE NECESSARY TO RUN THIS CODE:
!
! **      EFDC.INP
! **      CELL.INP
! **      CELLLT.INP
! **      DXDY.INP
! **      LXLY.INP
! **      SALT.INP OR RESTART.INP OR RESTRAN.INP
! **      TEMP.INP OR RESTART.INP OR RESTRAN.INP
! **      ASER.INP
! **      WSER.INP
! **      QSER.INP
! **      PSER.INP
! **      SSER.INP
! **      TSER.INP
! **      DSER.INP
! **      SDSER.INP
! **      SNSER.INP
! **      TXSER.INP
! **      SFSER.INP
! **      TXSER.INP
! **      QCTL.INP
! **      MASK.INP
! **      SHOW.INP
! **      VEGE.INP
! **      MODDXDY.INP
! **      MODCHAN.INP
! **      GWATER.INP
!
!**********************************************************************C
!
! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
!
! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
!
!----------------------------------------------------------------------C
!
! CHANGE RECORD
! DATE MODIFIED     BY                 DATE APPROVED    BY
! 02/15/2002        John Hamrick       02/15/2002       John Hamrick
!  modified calculation of dxdj and dydi
! 03/19/2002        John Hamrick       03/19/2002       John Hamrick
!  eliminated code section which adjust restart in information for
!  wet and dry see 'fixed moved to restin1'
!
!----------------------------------------------------------------------C
!
!**********************************************************************C
!

	INCLUDE 'EFDC.PAR'
	INCLUDE 'EFDC.CMN'
	INCLUDE 'ALLSET.INT'
	CHARACTER*80 TITLE
!
	DIMENSION CPUTIME(2)
!
!**********************************************************************C
!
! **  SET IEEE EXCEPTION TRAPS (SUN SYSTEMS ONLY)
!
!     CALL SETUP_IEEE
!
!**********************************************************************C
!
      TIMEEND=SECNDS(0.0)
!
      CALL WELCOME
!
!**********************************************************************C
!
! **  OPEN OUTPUT FILES
!
!----------------------------------------------------------------------C
!
      OPEN(7,FILE='EFDC.OUT',STATUS='UNKNOWN')
      OPEN(8,FILE='EFDCLOG.OUT',STATUS='UNKNOWN')
      OPEN(9,FILE='TIME.LOG',STATUS='UNKNOWN')
      OPEN(10,FILE='DRYWET.LOG',STATUS='UNKNOWN')
      OPEN(1,FILE='VSFP.OUT',STATUS='UNKNOWN')
!
      CLOSE(7,STATUS='DELETE')
      CLOSE(8,STATUS='DELETE')
      CLOSE(9,STATUS='DELETE')
      CLOSE(10,STATUS='DELETE')
      CLOSE(1,STATUS='DELETE')
!
      OPEN(7,FILE='EFDC.OUT',STATUS='UNKNOWN')
!      OPEN(8,FILE='EFDCLOG.OUT',STATUS='UNKNOWN')   !hnr 7/27/2009
      OPEN(9,FILE='TIME.LOG',STATUS='UNKNOWN')
!
      OPEN(1,FILE='SEDIAG.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
      OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
!
      OPEN(1,FILE='NEGSEDSND.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')

!
!**********************************************************************C
!
! **  SET HARDWIRES TO BE REMOVED IN SUBSEQUENT VERSIONS
!
      ISHOUSATONIC=0   ! SET TO 1 TO ACTIVATE HOUSATONIC OPTIONS
!
!**********************************************************************C
!
! **  CALL INPUT SUBROUTINE
!
      CALL INPUT(TITLE)
!
!**********************************************************************C
!
! **  CALL SUBROUTINE TO ADJUST, CONVERT AND SMOOTH DEPTH FIELD
!
      IF(NSHMAX.GE.1) CALL DEPSMTH
!
!**********************************************************************C
!
! **  SET TIME RELATED PARAMETERS
!
! **  THE PARAMETER NTC=NUMBER OF TIME CYCLES, CONTROLS
! **  THE LENGTH OF RUN (NUMBER OF TIME STEPS)
!
!----------------------------------------------------------------------C
!
      ISCRAY=0
!
      TCYCLE=0.0
      TLRPD=0.0
      THDMT=0.0
      TVDIF=0.0
      TCGRS=0.0
      TCONG=0.0
      TSADV=0.0
      TRELAXV=0.0
      TPUV=0.0
      TCEXP=0.0
      TAVB=0.0
      TUVW=0.0
      TQQQ=0.0
      TQCLT=0.0
      TWQADV=0.0
      TWQDIF=0.0
      TWQKIN=0.0
      TWQSED=0.0
      TSSED=0.0
	TMISC=0.0
	TTBXY=0.0
	TOUTPUT=0.0
      WTLRPD=0.0
      WTHDMT=0.0
      WTVDIF=0.0
      WTCGRS=0.0
      WTCONG=0.0
      WTSADV=0.0
      WTRELAXV=0.0
      WTPUV=0.0
      WTCEXP=0.0
      WTAVB=0.0
      WTUVW=0.0
      WTQQQ=0.0
      WTQCLT=0.0
      WTWQADV=0.0
      WTWQDIF=0.0
      WTWQKIN=0.0
      WTWQSED=0.0
      WTSSED=0.0
	WTMISC=0.0
	WTTBXY=0.0
	WTOUTPUT=0.0
      CFMAX=CF
      PI=3.1415926535898
      NBAN=49
      TPN=FLOAT(NTSPTC)
      NTS=NTC*NTSPTC/NFLTMT
      NLTS=NTSPTC*NLTC
      NTTS=NTSPTC*NTTC
      SNLT=0.
      NCTBC=1
      NPRINT=1
      NTSPP=NTCPP*NTSPTC/NFLTMT
      NTSNB=NTCNB*NTSPTC
      NTSVB=NTCVB*NTSPTC
      ITRMAX=0
      ITRMIN=1000
      ERRMAX=1E-9
      ERRMIN=1000.
      NMMT=1
      NBAL=1
      NBALE=1
      NBALO=1
      NBUD=1
      NHAR=1
      NTSPTC2=2*NTSPTC/NFLTMT
!     NDISP=NTS-NTSPTC+1
      NDISP=NTS-NTSPTC+2
      NSHOWR=0
      NSHOWC=0
!
      DO NS=1,NASER
        MATLAST(NS)=1
      ENDDO
!
      DO NS=1,NWSER
        MWTLAST(NS)=1
      ENDDO
!
      DO NS=1,NPSER
        MPTLAST(NS)=1
      ENDDO
!
      DO NS=1,NQSER
        MQTLAST(NS)=1
      ENDDO
!
      DO NS=1,NQWRSR
        MQWRTLST(NS)=1
      ENDDO
!
      NTMP=4+NSED+NSND+NTOX
      DO NC=1,NTMP
      DO NN=1,NCSER(NC)
        MCTLAST(NN,NC)=1
      ENDDO
      ENDDO
!
      MSFTLST=1
!
!**********************************************************************C
!
! **  SET CONTROLS FOR WRITING TO INSTANTANEOUS 2D SCALAR CONTOURING
! **  AND 2D VELOCITY VECTOR PLOTTING FILES
!
!----------------------------------------------------------------------C
!
! **  SCALAR FIELD CONTOURING IN HORIZONTAL PLANES: SUBROUTINE SALPLTH
!
      DO N=1,7
      IF(ISSPH(N).EQ.0)THEN
        NCSPH(N)=0
        JSSPH(N)=0
      ENDIF
      IF(ISSPH(N).EQ.2)THEN
        NCSPH(N)=NTS-(NTSPTC-(NTSPTC/NPSPH(N)))/NFLTMT
        JSSPH(N)=1
      ENDIF
      IF(ISSPH(N).EQ.1)THEN
        NCSPH(N)=NTSPTC/NPSPH(N)/NFLTMT
        JSSPH(N)=1
      ENDIF
      ENDDO
!
      IF(ISSPH(8).EQ.1)NCSPH(8)=NTSPTC/NPSPH(8)/NFLTMT
!
!----------------------------------------------------------------------C
!
! **  SEDIMENT BED PROPERTIES CONTOURING IN HORIZONTAL
! **  PLANES: SUBROUTINE SBEDPLT
!
      IF(ISBPH.EQ.0)THEN
        NCBPH=0
        JSBPH=0
      ENDIF
      IF(ISBPH.EQ.2)THEN
        NCBPH=NTS-(NTSPTC-(NTSPTC/NPBPH))/NFLTMT
        JSBPH=1
      ENDIF
      IF(ISBPH.EQ.1)THEN
        NCBPH=NTSPTC/NPBPH/NFLTMT
        JSBPH=1
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  FREE SURFACE ELEVATION OR PRESSURE CONTOURING IN HORIZONTAL
! **  PLANES: SUBROUTINE SURFPLT
!
      IF(ISPPH.EQ.0)THEN
        NCPPH=0
        JSPPH=0
      ENDIF
      IF(ISPPH.EQ.2)THEN
        NCPPH=NTS-(NTSPTC-(NTSPTC/NPPPH))/NFLTMT
        JSPPH=1
      ENDIF
      IF(ISPPH.EQ.1)THEN
        NCPPH=NTSPTC/NPPPH/NFLTMT
        JSPPH=1
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES: SUBROUTINE VELPLTH
!
      IF(ISVPH.EQ.0)THEN
        NCVPH=0
        JSVPH=0
      ENDIF
      IF(ISVPH.EQ.2)THEN
        NCVPH=NTS-(NTSPTC-(NTSPTC/NPVPH))/NFLTMT
        JSVPH=1
      ENDIF
      IF(ISVPH.EQ.1)THEN
        NCVPH=NTSPTC/NPVPH/NFLTMT
        JSVPH=1
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  SCALAR FIELD CONTOURING IN VERTICAL PLANES: SUBROUTINE SALPLTV
!
      DO N=1,7
      IF(ISSPV(N).EQ.0)THEN
        NCSPV(N)=0
        JSSPV(N)=0
      ENDIF
      IF(ISSPV(N).EQ.2)THEN
        NCSPV(N)=NTS-(NTSPTC-(NTSPTC/NPSPV(N)))/NFLTMT
        JSSPV(N)=1
         DO IS=1,ISECSPV
         CCTITLE(20+IS)=CCTITLE(10+IS)
         CCTITLE(30+IS)=CCTITLE(10+IS)
         CCTITLE(40+IS)=CCTITLE(10+IS)
         CCTITLE(50+IS)=CCTITLE(10+IS)
         ENDDO
      ENDIF
      IF(ISSPV(N).EQ.1)THEN
        NCSPV(N)=NTSPTC/NPSPV(N)/NFLTMT
        JSSPV(N)=1
         DO IS=1,ISECSPV
         CCTITLE(20+IS)=CCTITLE(10+IS)
         CCTITLE(30+IS)=CCTITLE(10+IS)
         CCTITLE(40+IS)=CCTITLE(10+IS)
         CCTITLE(50+IS)=CCTITLE(10+IS)
         ENDDO
      ENDIF
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  NORMAL VELOCITY CONTOURING AND TANGENTIAL VELOCITY VECTOR
! **  PLOTTING IN VERTICAL PALNES: SUBROUTINE VELPLTV
!
      IF(ISVPV.EQ.0)THEN
        NCVPV=0
        JSVPV=0
      ENDIF
      IF(ISVPV.EQ.2)THEN
        NCVPV=NTS-(NTSPTC-(NTSPTC/NPVPV))/NFLTMT
        JSVPV=1
         DO IS=1,ISECVPV
         CVTITLE(20+IS)=CVTITLE(10+IS)
         CVTITLE(30+IS)=CVTITLE(10+IS)
         CVTITLE(40+IS)=CVTITLE(10+IS)
         CVTITLE(50+IS)=CVTITLE(10+IS)
         CVTITLE(60+IS)=CVTITLE(10+IS)
         CVTITLE(70+IS)=CVTITLE(10+IS)
         CVTITLE(80+IS)=CVTITLE(10+IS)
         CVTITLE(90+IS)=CVTITLE(10+IS)
         ENDDO
      ENDIF
      IF(ISVPV.EQ.1)THEN
        NCVPV=NTSPTC/NPVPV/NFLTMT
        JSVPV=1
         DO IS=1,ISECVPV
         CVTITLE(20+IS)=CVTITLE(10+IS)
         CVTITLE(30+IS)=CVTITLE(10+IS)
         CVTITLE(40+IS)=CVTITLE(10+IS)
         CVTITLE(50+IS)=CVTITLE(10+IS)
         CVTITLE(60+IS)=CVTITLE(10+IS)
         CVTITLE(70+IS)=CVTITLE(10+IS)
         CVTITLE(80+IS)=CVTITLE(10+IS)
         CVTITLE(90+IS)=CVTITLE(10+IS)
         ENDDO
       ELSE
        NCVPV=0
        JSVPV=0
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  THREE-DIMENSIONAL HDF FORMAT GRAPHICS FILES: SUBROUTINE OUT3D
!
      IF(IS3DO.EQ.1)THEN
       NC3DO=NTS-(NTSPTC-(NTSPTC/NP3DO))/NFLTMT
      ENDIF
!
!**********************************************************************C
!
! **  SET CONTROLS FOR WRITING TO FILTERED, AVERAGED OR RESIDUAL
! **  2D SCALAR CONTOURING AND 2D VELOCITY VECTOR PLOTTING FILES
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATION
! **  CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
!
      DO N=1,7
      IF(ISRSPH(N).GE.1) JSRSPH(N)=1
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES:
! **  SUBROUTINE RVELPLTH
!
      IF(ISRVPH.GE.1) JSRVPH=1
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SURFACE ELEVATION PLOTTING IN HORIZONTAL PLANES:
! **  SUBROUTINE RVELPLTH
!
      IF(ISRPPH.GE.1) JSRPPH=1
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SCALAR FIELD CONTOURING IN VERTICAL
! **  PLANES: SUBROUTINE RSALPLTV
!
      DO N=1,7
      IF(ISRSPV(N).GE.1) JSRSPV(N)=1
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL NORMAL AND TANGENTIAL VELOCITY CONTOURING AND AND
! **  TANGENTIAL VELOCITY VECTOR PLOTTING IN VERTICAL PLANES:
! **  SUBROUTINE RVELPLTV
!
      IF(ISRVPV.GE.1) JSRVPV=1
!
!**********************************************************************C
!
! **  SET CONTROLS FOR WRITING TO DRIFTER, HARMONIC ANALYSIS,
! **  RESIDUAL TRANSPORT, AND BLANCE OUTPUT FILES
!
!----------------------------------------------------------------------C
!
      JSPD=1
      NCPD=1
!
      JSLSHA=1
      IF(ISLSHA.EQ.1)THEN
       LSLSHA=0
       NCLSHA=NTS-NTCLSHA*NTSPTC
      ENDIF
!
      IF(ISRESTR.EQ.1) JSRESTR=1
      JSWASP=0
      IF(ISWASP.GE.1) JSWASP=1
!
      IF(ISBAL.GE.1)THEN
        JSBAL=1
        JSBALO=1
        JSBALE=1
      ENDIF
      JSSBAL=1
!
!**********************************************************************C
!
! **  SET CONTROL FOR CALCULATION OF LAGRANGIAN MEAN VELOCITIY FIELDS
! **  BY PARTICLE TRACKING
!
!----------------------------------------------------------------------C
!
      IF(ISLRPD.GE.1)THEN
        NLRPDRT(1)=NTS-NTSPTC-(NTSPTC-(NTSPTC/MLRPDRT))/NFLTMT
        DO M=2,MLRPDRT
        NLRPDRT(M)=NLRPDRT(M-1)+(NTSPTC/MLRPDRT)
        ENDDO
        JSLRPD=1
       ELSE
        NLRPDRT(1)=NTS+2
        JSLRPD=0
      ENDIF
!
!**********************************************************************C
!
! **  SET SOME CONSTANTS
!
!----------------------------------------------------------------------C
!
      JSTBXY=0
!
      CTURB2=CTURB**0.667
      CTURB3=CTURB**0.333
!
      KS=KC-1
      IF(KS.EQ.0) KS=1
      DZI=FLOAT(KC)
      DZ=1./DZI
      DZS=DZ*DZ
!
      DT=TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)
      DTI=1./DT
      DT2=2.*DT
      DTMIN=DT
!
      AVCON1=2.*(1.-AVCON)*DZI*AVO
      G=9.81
      GPO=G*BSC
      GI=1./G
      GID2=.5*GI
      PI=3.1415926535898
      PI2=2.*PI
!
      TCVP=0.0625*TIDALP/PI
!
!**********************************************************************C
!
! **  SET CONSTANTS FOR M2 TIDAL CYCLE HARMONIC ANALYSIS
!
!----------------------------------------------------------------------C
!
      IF(ISHTA.GT.0)THEN
!
      AC=0.
      AS=0.
      ACS=0.
	TSHIFT=(TBEGIN*TCON/DT)+FLOAT(NTC-2)*NTSPTC
!
      DO N=1,NTSPTC
      TNT=FLOAT(N)+TSHIFT
      NP=NTSPTC+N
      WC(N)=COS(2.*PI*TNT/TPN)
      WS(N)=SIN(2.*PI*TNT/TPN)
      WC(NP)=WC(N)
      WS(NP)=WS(N)
      AC=AC + 2.*WC(N)*WC(N)
      AS=AS + 2.*WS(N)*WS(N)
!     ACS=ACS +2.*WC(N)*WS(N)
      ACS=0.
      WC2(N)=COS(4.*PI*TNT/TPN)
      WS2(N)=SIN(4.*PI*TNT/TPN)
      WC2(NP)=WC2(N)
      WS2(NP)=WC2(N)
      AC2=AC2 + 2.*WC2(N)*WC2(N)
      AS2=AS2 + 2.*WS2(N)*WS2(N)
!     ACS2=ACS2 +2.*WC2(N)*WS2(N)
      ACS2=0.
      ENDDO
!
      DET=AC*AS-ACS*ACS
      AS=AS/DET
      AC=AC/DET
      ACS=ACS/DET
      DET=AC2*AS2-ACS2*ACS2
      AS2=AS2/DET
      AC2=AC2/DET
      ACS2=ACS2/DET
!
      ENDIF
!
!**********************************************************************C
!
! **  SET WEIGHTS FOR SALINITY AND TEMPERATURE BOUNDARY INTERPOLATION
!
!----------------------------------------------------------------------C
!
      IF(KC.GT.1)THEN
      DO K=1,KC
      WTCI(K,1)=FLOAT(K-KC)/FLOAT(1-KC)
      WTCI(K,2)=FLOAT(K-1)/FLOAT(KC-1)
      ENDDO
      ELSE
      WTCI(1,1)=0.5
      WTCI(1,2)=0.5
      ENDIF
!
!**********************************************************************C
!
! **  INITIALIZE ARRAYS
!
      CALL AINIT
!
!**********************************************************************C
!
! **  READ IN XLON AND YLAT OR UTME AND UTMN OF CELL CENTERS OF
! **  CURVILINEAR PORTION OF THE  GRID
!
      IF(ISCLO.EQ.1)THEN
!
      OPEN(1,FILE='LXLY.INP',STATUS='UNKNOWN')
!
      DO NS=1,4
      READ(1,1111)
      ENDDO
 1111 FORMAT(80X)
!
      IF(ISCORV.EQ.1)THEN
       DO LL=1,LVC
       READ(1,*,ERR=3000)I,J,XLNUTME,YLTUTMN,CCUE,CCVE,CCUN,CCVN,TMPVAL,TMPCOR
       L=LIJ(I,J)
       DLON(L)=XLNUTME
       DLAT(L)=YLTUTMN
       ANG1=ATAN2(CCUN,CCUE)
       ANG2=ATAN2(-CCVE,CCVN)
       ANG=0.5*(ANG1+ANG2)
       CUE(L)=COS(ANG)
       CVE(L)=-SIN(ANG)
       CUN(L)=SIN(ANG)
       CVN(L)=COS(ANG)
       WINDSTKA(L)=TMPVAL
       FCORC(L)=TMPCOR
       DETTMP=1./( CUE(L)*CVN(L)-CUN(L)*CVE(L) )
       IF(DETTMP.EQ.0.0)THEN
         WRITE(6,6262)
         WRITE(6,6263)IL(L),JL(L)
         STOP
       ENDIF
       ENDDO
      ELSE
       DO LL=1,LVC
       READ(1,*,ERR=3000)I,J,XLNUTME,YLTUTMN,CCUE,CCVE,CCUN,CCVN,TMPVAL
       L=LIJ(I,J)
       DLON(L)=XLNUTME
       DLAT(L)=YLTUTMN
       ANG1=ATAN2(CCUN,CCUE)
       ANG2=ATAN2(-CCVE,CCVN)
       ANG=0.5*(ANG1+ANG2)
       CUE(L)=COS(ANG)
       CVE(L)=-SIN(ANG)
       CUN(L)=SIN(ANG)
       CVN(L)=COS(ANG)
       WINDSTKA(L)=TMPVAL
       FCORC(L)=CF
       DETTMP=1./( CUE(L)*CVN(L)-CUN(L)*CVE(L) )
       IF(DETTMP.EQ.0.0)THEN
         WRITE(6,6262)
         WRITE(6,6263)IL(L),JL(L)
         STOP
       ENDIF
       ENDDO
      ENDIF
!
 6262 FORMAT('  SINGULAR INVERSE TRANSFORM FROM E,N TO CURV X,Y')
 6263 FORMAT('  I,J =',2I10/)
!
      CLOSE(1)
!
      ENDIF
!
      FCORC(1)=FCORC(2)
      FCORC(LC)=FCORC(LA)
!
      GOTO 3002
 3000 WRITE(6,3001)
 3001 FORMAT('  READ ERROR FOR FILE LXLY.INP ')
      STOP
 3002 CONTINUE
!
!     COSTMP=COSD(15.)
!     SINTMP=SIND(15.)
!     SINNEG=-SINTMP
      OPEN(1,FILE='NEWLXLY.INP',STATUS='UNKNOWN')
      IF(ISCORV.EQ.1)THEN
       DO L=2,LA
        WRITE(1,1112)IL(L),JL(L),DLON(L),DLAT(L),CUE(L),CVE(L),CUN(L),CVN(L),WINDSTKA(L),FCORC(L)
       ENDDO
      ELSE
       DO L=2,LA
        WRITE(1,1112)IL(L),JL(L),DLON(L),DLAT(L),CUE(L),CVE(L),CUN(L), CVN(L),WINDSTKA(L)
       ENDDO
      ENDIF
      CLOSE(1)
!
      ZERO=0.
      OPEN(1,FILE='LIJMAP.OUT',STATUS='UNKNOWN')
      DO L=2,LA
      WRITE(1,1113)L,IL(L),JL(L),ZERO
      ENDDO
      CLOSE(1)
!
 1112 FORMAT (2I5,2F12.4,6F12.7)
 1113 FORMAT (3I5,F10.2)
!
!**********************************************************************C
!
! **  SET CORNER CELL STRESS CORRECTION
!
      DO L=2,LA
	  FSCORTBCV(L)=0.0
	ENDDO
!
      IF(ISCORTBC.GE.1)THEN
      DO L=2,LA
	  FSCORTBCV(L)=FSCORTBC
	ENDDO
	ENDIF
!
      IF(ISCORTBC.EQ.2)THEN
!
      OPEN(1,FILE='CORNERC.INP')
!
      DO NS=1,4
      READ(1,1111)
      ENDDO
!
      READ(1,*)NTMP
	DO NT=1,NTMP
	  READ(1,*)I,J,TMPVAL
	  L=LIJ(I,J)
	  FSCORTBCV(L)=TMPVAL
	ENDDO
!
      CLOSE(1)
!
      ENDIF
!
!**********************************************************************C!
!
! **  READ SPATIAL AVERAGING MAP FOR FOOD CHAIN MODEL OUTPUT
!
      IF(ISFDCH.EQ.1)THEN
!
	DO L=1,LC
	  MFDCHZ(L)=0
	ENDDO
!
      OPEN(1,FILE='FOODCHAIN.INP')
!
      DO NS=1,4
      READ(1,1111)
      ENDDO
!
      READ(1,*)NFDCHIJ
	DO LT=1,NFDCHIJ
	  READ(1,*)I,J,ITMPVAL
	  L=LIJ(I,J)
	  MFDCHZ(L)=ITMPVAL
	ENDDO
!
      CLOSE(1)
!
      ENDIF
!
!**********************************************************************C
!
! **  READ IN COUNTER CLOCKWISE ANGLE FROM EAST SPECIFYING
! **  PRINCIPAL FLOOD FLOW DIRECTION
!
      IF(ISTRAN(4).GE.1.AND.ISSFLFE.GE.1)THEN
!
      OPEN(1,FILE='FLDANG.INP',STATUS='UNKNOWN')
!
      DO LL=2,LA
      READ(1,*,ERR=3130)I,J,ANGTMP1,ANGTMP2
      L=LIJ(I,J)
      ACCWFLD(L,1)=0.0174533*ANGTMP1
      ACCWFLD(L,2)=0.0174533*ANGTMP2
      ENDDO
!
      CLOSE(1)
!
      ENDIF
!
      GOTO 3132
 3130 WRITE(6,3131)
 3131 FORMAT('  READ ERROR FOR FILE FLDANG.INP ')
      STOP
 3132 CONTINUE
!
!**********************************************************************C
!
! **  SET BOUNDARY CONDITION SWITCHES
!
      CALL SETBCS
!
!**********************************************************************C
!
! **  SET VERTICAL GRID DEPENDENT ARRAYS AND HARDWIRE DIMENSIONLESS
! **  MIXING LENGTH
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
!     DZC(K)=DZ
         DZIC(K)=1./DZC(K)
      ENDDO
!
      DZIG(0)=0.
      !DZIG(KC)=0.
      DZIG=0.
      DO K=1,KS
      DZG(K)=0.5*(DZC(K)+DZC(K+1))
      DZIG(K)=1./DZG(K)
      DZIGSD4(K)=0.25*DZIG(K)*DZIG(K)
      CDZU(K)=-DZC(K)/(DZC(K)+DZC(K+1))
      CDZL(K)=-DZC(K+1)/(DZC(K)+DZC(K+1))
      CDZF(K)=DZC(K)*DZC(K+1)/(DZC(K)+DZC(K+1))
      CDZM(K)=0.5*DZC(K)*DZC(K+1)
      ENDDO
!
      CDZR(1)=DZC(1)-1.
      CDZD(1)=DZC(1)
      DO K=2,KS
      CDZR(K)=DZC(K)+CDZR(K-1)
      CDZD(K)=DZC(K)+CDZD(K-1)
      ENDDO
!
      DO K=1,KS
      CDZR(K)=CDZR(K)*DZG(K)*CDZL(1)
      ENDDO
!
      CDZKMK(1)=0.
      DO K=2,KC
      CDZKMK(K)=DZIG(K-1)*DZIC(K)
      ENDDO
!
      DO K=1,KS
      CDZKK(K)=DZIC(K)*DZIG(K)
      CDZKKP(K)=DZIG(K)*DZIC(K+1)
      ENDDO
      CDZKK(KC)=0.
!
      Z(0)=0.
      IF(KC.GT.1)THEN
      DO K=1,KS
      Z(K)=Z(K-1)+DZC(K)
      ZZ(K)=Z(K)-0.5*DZC(K)
	IF(IFPROX.EQ.0)FPROX(K)=0.
      IF(IFPROX.EQ.1)FPROX(K)=1./(VKC*Z(K)*(1.-Z(K)))**2
      IF(IFPROX.EQ.2)FPROX(K)=(1./(VKC*Z(K))**2)+CTE5*(1./(VKC*(1.-Z(K)))**2)/(CTE4+0.00001)
      ENDDO
      ENDIF
!
      Z(KC)=Z(KS)+DZC(KC)
      ZZ(KC)=Z(KC)-0.5*DZC(KC)
!
      IF(ISRESTI.EQ.0)THEN
       DO K=0,KC
       DO L=1,LC
       DML(L,K)=VKC*Z(K)*(1.-Z(K))
       ENDDO
       ENDDO
      ENDIF
!
!**********************************************************************C
!
! **  READ GVC SPECIFIC FILES AND SET REAL AND LOGICAL MASK FOR GVC
!
      IF(IGRIDV.EQ.1) CALL SETGVC
!
!**********************************************************************C
!
! **  CALCUATE CURVATURE METRICS (NEW ADDITION)
!
!----------------------------------------------------------------------C
!
      DO L=1,LC
       DYDI(L)=0.
       DXDJ(L)=0.
      ENDDO
!
! ** DYDI
!
!      DO L=2,LA
!      IF(SUBO(L+1).GT.0.5) DYDI(L)=DYU(L+1)-DYU(L)
!      IF(SUBO(L+1).LT.0.5)THEN
!        IF(SUBO(L).LT.0.5) DYDI(L)=0.0
!        IF(SUBO(L).GT.0.5)THEN
!          DDYDDDX=2.*(DYP(L)-DYP(L-1))/(DXP(L)+DXP(L-1))
!          DYUP1=DYP(L)+0.5*DDYDDDX*DXP(L)
!          DYDI(L)=DYUP1-DYU(L)
!        ENDIF
!      ENDIF
!      ENDDO
!
      DO L=2,LA
	  I=IL(L)
	  J=JL(L)
	  IF(IJCT(I-1,J).GE.1.AND.IJCT(I-1,J).LE.5)THEN
	    IF(IJCT(I+1,J).GE.1.AND.IJCT(I+1,J).LE.5)THEN
	      DYDI(L)=DYU(L+1)-DYU(L)
	    ELSE
            DDYDDDX=2.*(DYP(L)-DYP(L-1))/(DXP(L)+DXP(L-1))
            DYUP1=DYP(L)+0.5*DDYDDDX*DXP(L)
            DYDI(L)=DYUP1-DYU(L)
	    END IF
	  ELSE
	    IF(IJCT(I+1,J).GE.1.AND.IJCT(I+1,J).LE.5)THEN
            DDYDDDX=2.*(DYP(L+1)-DYP(L))/(DXP(L+1)+DXP(L))
            DYUM1=DYP(L)-0.5*DDYDDDX*DXP(L)
            DYDI(L)=DYU(L)-DYUM1
	    ELSE
	      DYDI(L)=0.0
	    END IF
	  END IF
      ENDDO
!
! ** DXDJ
!
!      DO L=2,LA
!      LN=LNC(L)
!      LS=LSC(L)
!      IF(SVBO(LN).GT.0.5) DXDJ(L)=DXV(LN)-DXV(L)
!      IF(SVBO(LN).LT.0.5)THEN
!        IF(SVBO(L).LT.0.5) DXDJ(L)=0.0
!        IF(SVBO(L).GT.0.5)THEN
!          DDXDDDY=2.*(DXP(L)-DXP(LS))/(DYP(L)+DYP(LS))
!          DXVLN=DXP(L)+0.5*DDXDDDY*DYP(L)
!          DXDJ(L)=DXVLN-DXV(L)
!        ENDIF
!      ENDIF
!      ENDDO
!
      DO L=2,LA
	  LN=LNC(L)
	  LS=LSC(L)
	  I=IL(L)
	  J=JL(L)
	  IF(IJCT(I,J-1).GE.1.AND.IJCT(I,J-1).LE.5)THEN
	    IF(IJCT(I,J+1).GE.1.AND.IJCT(I,J+1).LE.5)THEN
	      DXDJ(L)=DXV(LN)-DXV(L)
	    ELSE
            DDXDDDY=2.*(DXP(L)-DXP(LS))/(DYP(L)+DYP(LS))
            DXVLN=DXP(L)+0.5*DDXDDDY*DYP(L)
            DXDJ(L)=DXVLN-DXV(L)
	    END IF
	  ELSE
	    IF(IJCT(I,J+1).GE.1.AND.IJCT(I,J+1).LE.5)THEN
            DDXDDDY=2.*(DXP(LN)-DXP(L))/(DYP(LN)+DYP(L))
            DXVLS=DXP(L)-0.5*DDXDDDY*DYP(L)
            DXDJ(L)=DXV(L)-DXVLS
	    ELSE
	      DXDJ(L)=0.0
	    END IF
	  END IF
      ENDDO
!
!**********************************************************************C
!
! **  READ RESTART CONDITIONS OR INITIALIZE SCALAR FIELDS
!
!     ISRESTI.EQ.10 READS AND OLD RESTART FILE GENERATED BY
!     PRE SEPTEMBER 8, 1992 VERSIONS OF EFDC.FOR
!
!----------------------------------------------------------------------C
!
      IF(ISLTMT.EQ.0)THEN
      IF(ISRESTI.GE.1)THEN
       IF(ISRESTI.EQ.1) CALL RESTIN1
       IF(ISRESTI.EQ.2) CALL RESTIN2
       IF(ISRESTI.EQ.10) CALL RESTIN10
      ENDIF
      IF(ISRESTI.EQ.-1) CALL RESTIN1
      ENDIF
!
      IF(ISRESTI.GE.1)THEN
	 DO K=1,KC
	 DO L=2,LA
	 FLOCDIA(L,K)=0.001
	 SEDFLOCDIA(L,K)=FLOCDIA(L,K)*SED(L,K,1)
	 ENDDO
	 ENDDO
	ENDIF
!----------------------------------------------------------------------C
!
! **  INTIALIZE SALINITY FIELD IF NOT READ IN FROM RESTART FILE
!
      IF(ISTRAN(1).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(1).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(1).GE.1)THEN
!
       NREST=0
!
       DO K=1,KC
       DO L=2,LA
       SAL(L,K)=SALINIT(L,K)
       SAL1(L,K)=SALINIT(L,K)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,1)=SALINIT(L,K)
       NLOS(LL,K,1)=0
       IF(NCSERS(LL,1).EQ.0)THEN
        SAL(L,K)=WTCI(K,1)*CBS(LL,1,1)+WTCI(K,2)*CBS(LL,2,1)
        SAL1(L,K)=SAL(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,1)=SALINIT(L,K)
       NLOW(LL,K,1)=0
       IF(NCSERW(LL,1).EQ.0)THEN
        SAL(L,K)=WTCI(K,1)*CBW(LL,1,1)+WTCI(K,2)*CBW(LL,2,1)
        SAL1(L,K)=SAL(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,1)=SALINIT(L,K)
       NLOE(LL,K,1)=0
       IF(NCSERE(LL,1).EQ.0)THEN
        SAL(L,K)=WTCI(K,1)*CBE(LL,1,1)+WTCI(K,2)*CBE(LL,2,1)
        SAL1(L,K)=SAL(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,1)=SALINIT(L,K)
       NLON(LL,K,1)=0
       IF(NCSERN(LL,1).EQ.0)THEN
        SAL(L,K)=WTCI(K,1)*CBN(LL,1,1)+WTCI(K,2)*CBN(LL,2,1)
        SAL1(L,K)=SAL(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       OPEN(1,FILE='NEWSALT.INP',STATUS='UNKNOWN')
       IONE=1
       WRITE(1,9101)IONE
       DO L=2,LC-1
       WRITE(1,9102)L,IL(L),JL(L),(SAL(L,K),K=1,KC)
       ENDDO
       CLOSE(1)
!
      ENDIF
      ENDIF
	ENDIF
!
 9101 FORMAT(I5)
 9102 FORMAT(3I5,28F8.2)
!
! **  INTIALIZE TEMP FIELD IF NOT READ IN FROM RESTART FILE
!
      IF(ISTRAN(2).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(2).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(2).EQ.1)THEN
!
       NREST=0
!
       DO K=1,KC
       DO L=2,LA
       TEM(L,K)=TEMINIT(L,K)
       TEM1(L,K)=TEM(L,K)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,2)=TEMINIT(L,K)
       NLOS(LL,K,2)=0
       IF(NCSERS(LL,2).EQ.0)THEN
        TEM(L,K)=WTCI(K,1)*CBS(LL,1,2)+WTCI(K,2)*CBS(LL,2,2)
        TEM1(L,K)=TEM(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,2)=TEMINIT(L,K)
       NLOW(LL,K,2)=0
       IF(NCSERW(LL,2).EQ.0)THEN
        TEM(L,K)=WTCI(K,1)*CBW(LL,1,2)+WTCI(K,2)*CBW(LL,2,2)
        TEM1(L,K)=TEM(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,2)=TEMINIT(L,K)
       NLOE(LL,K,2)=0
       IF(NCSERE(LL,2).EQ.0)THEN
        TEM(L,K)=WTCI(K,1)*CBE(LL,1,2)+WTCI(K,2)*CBE(LL,2,2)
        TEM1(L,K)=TEM(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,2)=TEMINIT(L,K)
       NLON(LL,K,2)=0
       IF(NCSERN(LL,2).EQ.0)THEN
        TEM(L,K)=WTCI(K,1)*CBN(LL,1,2)+WTCI(K,2)*CBN(LL,2,2)
        TEM1(L,K)=TEM(L,K)
       ENDIF
       ENDDO
       ENDDO
!
      ENDIF
      ENDIF
	ENDIF
!
! **  INTIALIZE TEMPERATURE BC IF NOT READ IN FROM RESTART FILE
!     AND CONSTANT INTIAL CONDITION IS USED
!
      IF(ISTRAN(2).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(2).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(2).GE.2)THEN
       M=2
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=TEMO
       NLOS(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=TEMO
       NLOW(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=TEMO
       NLOE(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=TEMO
       NLON(LL,K,M)=0
       ENDDO
       ENDDO
!
      ENDIF
      ENDIF
	ENDIF
!
! **  INTIALIZE DYE FIELD IF NOT READ IN FROM RESTART FILE
!
      IF(ISTRAN(3).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(3).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(3).GE.1)THEN
!
       NREST=0
!
       DO K=1,KC
       DO L=2,LA
       DYE(L,K)=DYEINIT(L,K)
       DYE1(L,K)=DYE(L,K)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,3)=DYEINIT(L,K)
       NLOS(LL,K,3)=0
       IF(NCSERS(LL,3).EQ.0)THEN
        DYE(L,K)=WTCI(K,1)*CBS(LL,1,3)+WTCI(K,2)*CBS(LL,2,3)
        DYE1(L,K)=DYE(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,3)=DYEINIT(L,K)
       NLOW(LL,K,3)=0
       IF(NCSERW(LL,3).EQ.0)THEN
        DYE(L,K)=WTCI(K,1)*CBW(LL,1,3)+WTCI(K,2)*CBW(LL,2,3)
        DYE1(L,K)=DYE(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,3)=DYEINIT(L,K)
       NLOE(LL,K,3)=0
       IF(NCSERE(LL,3).EQ.0)THEN
        DYE(L,K)=WTCI(K,1)*CBE(LL,1,3)+WTCI(K,2)*CBE(LL,2,3)
        DYE1(L,K)=DYE(L,K)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,3)=DYEINIT(L,K)
       NLON(LL,K,3)=0
       IF(NCSERN(LL,3).EQ.0)THEN
        DYE(L,K)=WTCI(K,1)*CBN(LL,1,3)+WTCI(K,2)*CBN(LL,2,3)
        DYE1(L,K)=DYE(L,K)
       ENDIF
       ENDDO
       ENDDO
!
      ENDIF
      ENDIF
	ENDIF
!
! **  INTIALIZE DYE BC IF NOT READ IN FROM RESTART FILE
! **  AND CONSTANT INITIAL CONDITIONS ARE USED
!
      IF(ISTRAN(3).GE.1)THEN
	IFLAG=0
      IF(ISRESTI.EQ.0)IFLAG=1
	IF(ISRESTI.EQ.1.AND.ISCI(3).EQ.0)IFLAG=1
      IF(IFLAG.EQ.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(3).EQ.0)THEN
!
       M=3
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=0.
       NLOS(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=0.
       NLOW(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=0.
       NLOE(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=0.
       NLON(LL,K,M)=0
       ENDDO
       ENDDO
!
      ENDIF
      ENDIF
	ENDIF
!
! **  INTIALIZE TOX AND BC IF NOT READ IN FROM RESTART FILE
! **  AND VARIABLE INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(5).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(5).EQ.1)THEN
       DO NT=1,NTOX
       IF(ITXINT(NT).EQ.1.OR.ITXINT(NT).EQ.3)THEN
       M=4+NT

       DO K=1,KC
       DO L=2,LA
       TOX(L,K,NT)=TOXINIT(L,K,NT)
       TOX1(L,K,NT)=TOX(L,K,NT)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,M)=TOXINIT(L,K,NT)
       NLOS(LL,K,M)=0
       IF(NCSERS(LL,M).EQ.0)THEN
        TOX(L,K,NT)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
        TOX1(L,K,NT)=TOX(L,K,NT)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,M)=TOXINIT(L,K,NT)
       NLOW(LL,K,M)=0
       IF(NCSERW(LL,M).EQ.0)THEN
        TOX(L,K,NT)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
        TOX1(L,K,NT)=TOX(L,K,NT)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,M)=TOXINIT(L,K,NT)
       NLOE(LL,K,M)=0
       IF(NCSERE(LL,3).EQ.0)THEN
        TOX(L,K,NT)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
        TOX1(L,K,NT)=TOX(L,K,NT)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,M)=TOXINIT(L,K,NT)
       NLON(LL,K,M)=0
       IF(NCSERN(LL,M).EQ.0)THEN
        TOX(L,K,NT)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
        TOX1(L,K,NT)=TOX(L,K,NT)
       ENDIF
       ENDDO
       ENDDO
!
      ENDIF
      ENDDO
      ENDIF
!
! **  INTIALIZE TOX BC IF NOT READ IN FROM RESTART FILE
! **  AND CONSTANT INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(5).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(5).EQ.1)THEN
       DO NT=1,NTOX
       IF(ITXINT(NT).EQ.0.OR.ITXINT(NT).EQ.2)THEN
       M=4+NT
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=TOXINTW(NT)
       NLOS(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=TOXINTW(NT)
       NLOW(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=TOXINTW(NT)
       NLOE(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=TOXINTW(NT)
       NLON(LL,K,M)=0
       ENDDO
       ENDDO
!
      ENDIF
      ENDDO
      ENDIF
!
! **  INTIALIZE TOX BED IF NOT READ IN FROM RESTART FILE
! **  AND VARIABLE INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(5).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(5).EQ.1)THEN
       DO NT=1,NTOX
       IF(ITXINT(NT).EQ.2.OR.ITXINT(NT).EQ.3)THEN
!
       DO K=1,KB
       DO L=2,LA
       TOXB(L,K,NT)=TOXBINIT(L,K,NT)
       TOXB1(L,K,NT)=TOXB(L,K,NT)
       ENDDO
       ENDDO
!
       ENDIF
       ENDDO
      ENDIF
!
! **  INTIALIZE SED AND BC IF NOT READ IN FROM RESTART FILE
! **  AND VARIABLE INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(6).EQ.1)THEN
       IF(ISEDINT.EQ.1.OR.ISEDINT.EQ.3)THEN
       DO NS=1,NSED
       M=4+NTOX+NS

       DO K=1,KC
       DO L=2,LA
       SED(L,K,NS)=SEDINIT(L,K,NS)
       SED1(L,K,NS)=SED(L,K,NS)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,M)=SEDINIT(L,K,NS)
       NLOS(LL,K,M)=0
       IF(NCSERS(LL,M).EQ.0)THEN
        SED(L,K,NS)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
        SED1(L,K,NS)=SED(L,K,NS)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,M)=SEDINIT(L,K,NS)
       NLOW(LL,K,M)=0
       IF(NCSERW(LL,M).EQ.0)THEN
        SED(L,K,NS)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
        SED1(L,K,NS)=SED(L,K,NS)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,M)=SEDINIT(L,K,NS)
       NLOE(LL,K,M)=0
       IF(NCSERE(LL,3).EQ.0)THEN
        SED(L,K,NS)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
        SED1(L,K,NS)=SED(L,K,NS)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,M)=SEDINIT(L,K,NS)
       NLON(LL,K,M)=0
       IF(NCSERN(LL,M).EQ.0)THEN
        SED(L,K,NS)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
        SED1(L,K,NS)=SED(L,K,NS)
       ENDIF
       ENDDO
       ENDDO
!
      ENDDO
!
	 DO K=1,KC
	 DO L=2,LA
	 FLOCDIA(L,K)=0.001
	 SEDFLOCDIA(L,K)=FLOCDIA(L,K)*SED(L,K,1)
	 ENDDO
	 ENDDO
!
      ENDIF
      ENDIF
!
! **  INTIALIZE SED BC IF NOT READ IN FROM RESTART FILE AND
! **  CONSTANT INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(6).EQ.1)THEN
      IF(ISEDINT.EQ.0.OR.ISEDINT.EQ.2)THEN
      DO NS=1,NSED
       M=4+NTOX+NS
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=SEDO(NS)
       NLOS(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=SEDO(NS)
       NLOW(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=SEDO(NS)
       NLOE(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=SEDO(NS)
       NLON(LL,K,M)=0
       ENDDO
       ENDDO
!
      ENDDO
      ENDIF
      ENDIF
!
! **  INTIALIZE SED BED IF NOT READ IN FROM RESTART FILE
! **  AND VARIABLE INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(6).EQ.1)THEN
       IF(ISEDINT.EQ.2.OR.ISEDINT.EQ.3)THEN
!
       DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
       SEDB(L,K,NS)=SEDBINIT(L,K,NS)
       SEDB1(L,K,NS)=SEDB(L,K,NS)
       ENDDO
       ENDDO
       ENDDO
!
       ENDIF
      ENDIF
!
! **  INTIALIZE SND AND BC IF NOT READ IN FROM RESTART FILE
! **  AND VARIABLE INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(7).EQ.1)THEN
       IF(ISEDINT.EQ.1.OR.ISEDINT.EQ.3)THEN
       DO NS=1,NSND
       M=4+NTOX+NSED+NS

       DO K=1,KC
       DO L=2,LA
       SND(L,K,NS)=SNDINIT(L,K,NS)
       SND1(L,K,NS)=SND(L,K,NS)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,M)=SNDINIT(L,K,NS)
       NLOS(LL,K,M)=0
       IF(NCSERS(LL,M).EQ.0)THEN
        SND(L,K,NS)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
        SND1(L,K,NS)=SND(L,K,NS)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,M)=SNDINIT(L,K,NS)
       NLOW(LL,K,M)=0
       IF(NCSERW(LL,M).EQ.0)THEN
        SND(L,K,NS)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
        SND1(L,K,NS)=SND(L,K,NS)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,M)=SNDINIT(L,K,NS)
       NLOE(LL,K,M)=0
       IF(NCSERE(LL,3).EQ.0)THEN
        SND(L,K,NS)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
        SND1(L,K,NS)=SND(L,K,NS)
       ENDIF
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,M)=SNDINIT(L,K,NS)
       NLON(LL,K,M)=0
       IF(NCSERN(LL,M).EQ.0)THEN
        SND(L,K,NS)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
        SND1(L,K,NS)=SED(L,K,NS)
       ENDIF
       ENDDO
       ENDDO
!
      ENDDO
      ENDIF
      ENDIF
!
! **  INTIALIZE SND BC IF NOT READ IN FROM RESTART FILE AND
! **  CONSTANT INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(7).EQ.1)THEN
      IF(ISEDINT.EQ.0.OR.ISEDINT.EQ.2)THEN
       DO NX=1,NSND
       NS=NSED+NX
       M=4+NTOX+NSED+NX
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=SEDO(NS)
       NLOS(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=SEDO(NS)
       NLOW(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=SEDO(NS)
       NLOE(LL,K,M)=0
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=SEDO(NS)
       NLON(LL,K,M)=0
       ENDDO
       ENDDO
!
      ENDDO
      ENDIF
      ENDIF
!
! **  INTIALIZE SND BED IF NOT READ IN FROM RESTART FILE
! **  AND VARIABLE INITIAL CONDITIONS ARE USED
!
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
!
      IF(IISTMP.EQ.0.AND.ISTRAN(7).EQ.1)THEN
       IF(ISEDINT.EQ.2.OR.ISEDINT.EQ.3)THEN
!
       DO NX=1,NSND
       DO K=1,KB
       DO L=2,LA
       SNDB(L,K,NX)=SNDBINIT(L,K,NX)
       SNDB1(L,K,NX)=SNDB(L,K,NX)
       ENDDO
       ENDDO
       ENDDO
!
       ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  CORRECT CONCENTRATIONS FOR GVC
!
      IF(IGRIDV.EQ.1)THEN
!
       IF(ISTRAN(1).GE.1)THEN
         DO K=1,KC
         DO L=1,LC
          IF(K.LT.KGVCP(L)) SAL(L,K)=0.0
          IF(K.LT.KGVCP(L)) SAL1(L,K)=0.0
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(2).GE.1)THEN
         DO K=1,KC
         DO L=1,LC
          IF(K.LT.KGVCP(L)) TEM(L,K)=0.0
          IF(K.LT.KGVCP(L)) TEM1(L,K)=0.0
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(3).GE.1)THEN
         DO K=1,KC
         DO L=1,LC
          IF(K.LT.KGVCP(L)) DYE(L,K)=0.0
          IF(K.LT.KGVCP(L)) DYE1(L,K)=0.0
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(5).GE.1)THEN
         DO NT=1,NTOX
          DO K=1,KC
          DO L=1,LC
           IF(K.LT.KGVCP(L)) TOX(L,K,NT)=0.0
          ENDDO
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(6).GE.1)THEN
         DO NS=1,NSED
          DO K=1,KC
          DO L=1,LC
           IF(K.LT.KGVCP(L)) SED(L,K,NS)=0.0
          ENDDO
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(7).GE.1)THEN
        DO NS=1,NSND
         DO K=1,KC
         DO L=1,LC
          IF(K.LT.KGVCP(L)) SND(L,K,NS)=0.0
         ENDDO
        ENDDO
         ENDDO
       ENDIF
!
      ENDIF
!
!**********************************************************************C
!
! **  INITIALIZE SEDIMENT BED
!
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1) CALL BEDINIT
!     IF(NSED.GE.1.OR.NSND.GE.1) CALL BEDINIT
!
!**********************************************************************C
!
! **  INITIALIZE BUOYANCY AND EQUATION OF STATE
!
      CALL CALBUOY
!
!**********************************************************************C
!
! **  INITIALIZE SFL IF(ISRESTI.EQ.0.AND ISTRAN(4).GE.1)
!
      IF(ISRESTI.EQ.0.AND.ISTRAN(4).GE.1)THEN
      IF(ISLTMT.EQ.0.AND.ISTOPT(4).EQ.11)THEN
!
       DO K=1,KC
       DO L=1,LC
       SFL(L,K)=SAL(L,K)
       SFL2(L,K)=SAL(L,K)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,5)=SALINIT(L,K)
       NLOS(LL,K,5)=0
       SFL(L,K)=WTCI(K,1)*CBS(LL,1,5)+WTCI(K,2)*CBS(LL,2,5)
       SFL2(L,K)=SFL(L,K)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,5)=SALINIT(L,K)
       NLOW(LL,K,5)=0
       SFL(L,K)=WTCI(K,1)*CBW(LL,1,5)+WTCI(K,2)*CBW(LL,2,5)
       SFL2(L,K)=SFL(L,K)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,5)=SALINIT(L,K)
       NLOE(LL,K,5)=0
       SFL(L,K)=WTCI(K,1)*CBE(LL,1,5)+WTCI(K,2)*CBE(LL,2,5)
       SFL2(L,K)=SFL(L,K)
       ENDDO
       ENDDO
!
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,5)=SALINIT(L,K)
       NLON(LL,K,5)=0
       SFL(L,K)=WTCI(K,1)*CBN(LL,1,5)+WTCI(K,2)*CBN(LL,2,5)
       SFL2(L,K)=SFL(L,K)
       ENDDO
       ENDDO
!
      ENDIF
      ENDIF
!
!**********************************************************************C
!
! **  ACTIVATE DYE TRACER CONTINUITY CHECK
!
!----------------------------------------------------------------------C
!
      IF(ISMMC.EQ.1)THEN
!
      DO K=1,KC
      DO L=1,LC
      DYE(L,K)=1.
      DYE1(L,K)=1.
      ENDDO
      ENDDO
!
      DO K=1,KC
!
      DO LL=1,NCBS
      CLOS(LL,K,3)=1.
      NLOS(LL,K,3)=0
      ENDDO
!
      DO LL=1,NCBW
      CLOW(LL,K,3)=1.
      NLOW(LL,K,3)=0
      ENDDO
!
      DO LL=1,NCBE
      CLOE(LL,K,3)=1.
      NLOE(LL,K,3)=0
      ENDDO
!
      DO LL=1,NCBN
      CLON(LL,K,3)=1.
      NLON(LL,K,3)=0
      ENDDO
!
      ENDDO
!
      ENDIF
!
!**********************************************************************C
!
! **  MASK CELLS TO BE CONVERTED FROM WATER TO LAND
!
! **  CALL MOVED TO SETBCS ON 23 JAN 2004
!
!      IF(ISMASK.EQ.1) CALL CELLMASK
!
!**********************************************************************C
!
! **  SET VERTICAL GRID DEPENDENT ARRAYS AND HARDWIRE DIMENSIONLESS
! **  MIXING LENGTH
!
!     MOVED TO JUST BEFORE CALL SETGVC
!
!**********************************************************************C
!
! **  INITIALIZE UNSTRETCHING PROCEDURE
!
!----------------------------------------------------------------------C
!
      DZPC=(SELVMAX-BELVMIN)/FLOAT(KPC)
!
      ZP(0)=BELVMIN
      DO KP=1,KPC
      ZP(KP)=ZP(KP-1)+DZPC
      ENDDO
!
      DO KP=1,KPC
      ZZP(KP)=0.5*(ZP(KP)+ZP(KP-1))
      ENDDO
!
      DO L=2,LA
      TMP=(BELV(L)-BELVMIN)/DZPC
      KPB(L)=NINT(0.5+TMP)
      ENDDO
!
!**********************************************************************C
!
! **  CALCULATE CONSTANT HORIZONTAL SPATIAL ARRAYS
!
!----------------------------------------------------------------------C
!
      DO L=2,LA
      DXYU(L)=DXU(L)*DYU(L)
      DXYV(L)=DXV(L)*DYV(L)
      DXYP(L)=STCAP(L)*DXP(L)*DYP(L)
      DXIU(L)=1./DXU(L)
      DYIU(L)=1./DYU(L)
      DXIV(L)=1./DXV(L)
      DYIV(L)=1./DYV(L)
      DXYIP(L)=1./(STCAP(L)*DXP(L)*DYP(L))
      DXYIU(L)=1./(DXU(L)*DYU(L))
      DXYIV(L)=1./(DXV(L)*DYV(L))
      HRU(L)=SUB(L)*HMU(L)*DYU(L)*DXIU(L)
      HRV(L)=SVB(L)*HMV(L)*DXV(L)*DYIV(L)
!      HRUO(L)=SUB(L)*DYU(L)*DXIU(L)
!      HRVO(L)=SVB(L)*DXV(L)*DYIV(L)
      HRUO(L)=SUBO(L)*DYU(L)*DXIU(L)
      HRVO(L)=SVBO(L)*DXV(L)*DYIV(L)
      SBX(L)=0.5*SBX(L)*DYU(L)
      SBY(L)=0.5*SBY(L)*DXV(L)
      SBXO(L)=SBX(L)
      SBYO(L)=SBY(L)
      SNLPX(L)=GID2*SNLPX(L)*DYU(L)
      SNLPY(L)=GID2*SNLPY(L)*DXV(L)
      ENDDO
!
!**********************************************************************C
!
! **  DEACTIVATE DRY CELLS
!
!     IF(ISDRY.GE.1.AND.ISDRY.LE.3)THEN
!       OPEN(1,FILE='DRYWET.LOG',POSITION='APPEND',STATUS='UNKNOWN')
!       DO L=2,LA
!       IF(HP(L).LE.HDRY)THEN
!         LN=LNC(L)
!         NTMP=0
!         WRITE(1,6902)NTMP,IL(L),JL(L),HP(L),H1P(L),H2P(L)
!         WRITE(6,6902)NTMP,IL(L),JL(L),HP(L),H1P(L),H2P(L)
!         WRITE(8,6902)NTMP,IL(L),JL(L),HP(L),H1P(L),H2P(L)
!         NATDRY(L)=0
!         ISCDRY(L)=2
!         SUB(L)=0.
!         SUB(L+1)=0.
!         SVB(L)=0.
!         SVB(LN)=0.
!         SBX(L)=0.
!         SBX(L+1)=0.
!         SBY(L)=0.
!         SBY(LN)=0.
!       ENDIF
!       ENDDO
!       CLOSE(1)
!     ENDIF
!
 6902 FORMAT('  DRYING AT N,I,J =',I10,2I6,'  H,H1,H2 =',3(2X,E12.4))
!
!**********************************************************************C
!
! **  CHECK FOR DRYING AND SET SWITCHES ON RESTART
!    'fixed moved to restin1'
!
!----------------------------------------------------------------------C
!!c      IF(ISDRY.GE.1.AND.ISRESTI.EQ.1)THEN
!      OPEN(1,FILE='DRYWET.LOG',POSITION='APPEND',STATUS='UNKNOWN')
!
!      DO L=2,LA
!      LS=LSC(L)
!      LN=LNC(L)
!      HUTMPP=0.5*(HP(L)+HP(L-1))
!      IF(HUTMPP.LE.HUWET(L))THEN
!        SUB(L)=0.
!        SBX(L)=0.
!      ENDIF
!      HUTMPP=0.5*(HP(L)+HP(L+1))
!      IF(HUTMPP.LE.HUWET(L+1))THEN
!        SUB(L+1)=0.
!        SBX(L+1)=0.
!      ENDIF
!      HVTMPP=0.5*(HP(L)+HP(LS))
!      IF(HVTMPP.LE.HVWET(L))THEN
!        SVB(L)=0.
!        SBY(L)=0.
!      ENDIF
!      HVTMPP=0.5*(HP(L)+HP(LN))
!      IF(HVTMPP.LE.HVWET(LN))THEN
!        SVB(LN)=0.
!        SBY(LN)=0.
!      ENDIF
!      IF(HP(L).LE.HDRY)THEN
!        SUB(L)=0.
!        SVB(L)=0.
!        SUB(L+1)=0.
!        SVB(LN)=0.
!        SBX(L)=0.
!        SBY(L)=0.
!        SBX(L+1)=0.
!        SBY(LN)=0.
!      ENDIF
!      ENDDO
!
!      CLOSE(1)
!      ENDIF
!
!**********************************************************************C
!
! **  INITIALIZE ZERO DIMENSION VOLUME BALANCE
!
!----------------------------------------------------------------------C
!
      IF(ISDRY.GE.1.AND.ISDRY.LE.98)THEN
        OPEN(1,FILE='ZVOLBAL.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='AVSEL.OUT',STATUS='UNKNOWN')
        LPBTMP=0
        DO L=2,LA
        ISLUSED(L)=0
        IF(SPB(L).EQ.0)THEN
          LPBTMP=LPBTMP+1
          ISLUSED(L)=1
        ENDIF
        LORDER(L)=0
        ENDDO
        ISLUSED(1)=1
        ISLUSED(LC)=1
        LORMAX=LC-2-LPBTMP
        DO LS=1,LORMAX
        BELMIN=100000.
          DO L=2,LA
           IF(SPB(L).NE.0.AND.ISLUSED(L).NE.1)THEN
             IF(BELV(L).LT.BELMIN)THEN
               LBELMIN=L
               BELMIN=BELV(L)
             ENDIF
           ENDIF
           ENDDO
         LORDER(LS)=LBELMIN
         ISLUSED(LBELMIN)=1
        ENDDO
        WRITE(1,5300)
        LS=1
        L=LORDER(LS)
        BELSURF(LS)=BELV(L)
        ASURFEL(LS)=DXYP(L)
        VOLSEL(LS)=0.
        WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
        DO LS=2,LORMAX
          L=LORDER(LS)
          BELSURF(LS)=BELV(L)
          ASURFEL(LS)=ASURFEL(LS-1)+DXYP(L)
          VOLSEL(LS)=VOLSEL(LS-1)+0.5*(BELSURF(LS)-BELSURF(LS-1))*(ASURFEL(LS)+ASURFEL(LS-1))
        WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
        ENDDO
        LS=LORMAX+1
        BELSURF(LS)=BELV(L)+10.0
        ASURFEL(LS)=ASURFEL(LS-1)
        VOLSEL(LS)=VOLSEL(LS-1)+0.5*(BELSURF(LS)-BELSURF(LS-1))*(ASURFEL(LS)+ASURFEL(LS-1))
        WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
        VOLZERD=0.
        VOLLDRY=0.
        DO L=2,LA
        IF(SPB(L).NE.0)THEN
          VOLZERD=VOLZERD+DXYP(L)*HP(L)
          IF(HP(L).GT.HDRY) VOLLDRY=VOLLDRY+DXYP(L)*HP(L)
        ENDIF
        ENDDO
        DO LS=1,LORMAX
        IF(VOLZERD.GE.VOLSEL(LS).AND.VOLZERD.LT.VOLSEL(LS+1))THEN
           WTM=VOLSEL(LS+1)-VOLZERD
           WTMP=VOLZERD-VOLSEL(LS)
           DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
           WTM=WTM/DELVOL
           WTMP=WTMP/DELVOL
           SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
           ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
        ENDIF
        ENDDO
        VETZERD=VOLZERD
        WRITE(1,5302)
        WRITE(1,5303) SELZERD,ASFZERD,VOLZERD,VOLLDRY
        CLOSE(1)
      ENDIF
!
 5300 FORMAT('   M    BELSURF     ASURFEL     ','   VOLSEL',/)
 5301 FORMAT(1X,I3,2X,F10.5,2X,E12.4,2X,E12.4)
 5302 FORMAT(/)
 5303 FORMAT(2X,F10.5,3(2X,E12.4))
!
!**********************************************************************C
!
! **  INITIALIZE ELEVATION OF ACTIVE GROUNDWATER ZONE FOR COLD START
!
!----------------------------------------------------------------------C
!
!       WRITE(6,5301) ISGWIE,DAGWZ,RNPOR,RIFTRM
!
        IF(ISGWIE.GE.1.AND.ISRESTI.EQ.0)THEN
!
        DO L=2,LA
          IF(HP(L).GT.HDRY)THEN
            AGWELV(L)=BELV(L)
           ELSE
            IF(BELAGW(L).LT.SELZERD)THEN
              AGWELV(L)=SELZERD
              AGWELV(L)=MIN(AGWELV(L),BELV(L))
             ELSE
              AGWELV(L)=BELAGW(L)
            ENDIF
          ENDIF
        ENDDO
!
        DO L=2,LA
        AGWELV1(L)=AGWELV(L)
        AGWELV2(L)=AGWELV(L)
        ENDDO
!
        OPEN(1,FILE='GWELV.OUT',STATUS='UNKNOWN')
        WRITE(1,5400)
        WRITE(1,5402)
        DO L=2,LA
        WRITE(1,5401)IL(L),JL(L),BELV(L),BELAGW(L),AGWELV(L)
        ENDDO
        CLOSE(1)
!
        ENDIF
!
!
 5400 FORMAT('   I   J    BELELV      BELAGW     ','  AGWELV',/)
 5401 FORMAT(1X,2I5,2X,F10.5,2X,F10.5,2X,F10.5)
 5402 FORMAT(/)
!
!**********************************************************************C
!
! **  CALCULATE CONSTANT C ARRAYS FOR EXTERNAL P SOLUTION
! **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV
! **  DXYIP=1/(DXP*DYP)
!
!----------------------------------------------------------------------C
!
      IF(ISLTMT.EQ.0)THEN
!
!----------------------------------------------------------------------C
!
      IF(IRVEC.NE.9)THEN
!
      DO L=2,LA
      CC(L)=1.
      CCC(L)=1.
      ENDDO
!
      IF(ISRLID.EQ.1)THEN
       DO L=2,LA
       CC(L)=0.
       CCC(L)=0.
       IF(SPB(L).EQ.0.) CC(L)=1.
       IF(SPB(L).EQ.0.) CCC(L)=1.
       ENDDO
      ENDIF
!
      DO L=2,LA
      LN=LNC(L)
      C1=-G*DT*DT*SPB(L)*DXYIP(L)
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(L+1)
      CN(L)=C1*HRV(LN)
      CC(L)=CC(L)-CS(L)-CW(L)-CE(L)-CN(L)
      CCI(L)=1./CC(L)
      CCS(L)=0.25*CS(L)
      CCW(L)=0.25*CW(L)
      CCE(L)=0.25*CE(L)
      CCN(L)=0.25*CN(L)
      CCC(L)=CCC(L)-CCS(L)-CCW(L)-CCE(L)-CCN(L)
      CCCI(L)=1./CCC(L)
      ENDDO
!
      DO LR=1,NRC
      L=LRC(LR)
      CCSR(LR)=CCS(L)*CCCI(L)
      CCWR(LR)=CCW(L)*CCCI(L)
      CCER(LR)=CCE(L)*CCCI(L)
      CCNR(LR)=CCN(L)*CCCI(L)
      CSR(LR)=CS(L)*CCI(L)
      CWR(LR)=CW(L)*CCI(L)
      CER(LR)=CE(L)*CCI(L)
      CNR(LR)=CN(L)*CCI(L)
      ENDDO
!
      DO LB=1,NBC
      L=LBC(LB)
      CCSB(LB)=CCS(L)*CCCI(L)
      CCWB(LB)=CCW(L)*CCCI(L)
      CCEB(LB)=CCE(L)*CCCI(L)
      CCNB(LB)=CCN(L)*CCCI(L)
      CSB(LB)=CS(L)*CCI(L)
      CWB(LB)=CW(L)*CCI(L)
      CEB(LB)=CE(L)*CCI(L)
      CNB(LB)=CN(L)*CCI(L)
      ENDDO
!
      DO LL=1,NPBW
      L=LPBW(LL)
      IF(ISPBW(LL).EQ.0)THEN
       IF(ISRED(L).EQ.1)THEN
        LR=LLRC(L)
        CER(LR)=0.
        CCER(LR)=0.
       ELSE
        LB=LLBC(L)
        CEB(LB)=0.
        CCEB(LB)=0.
       ENDIF
      ENDIF
      ENDDO
!
      DO LL=1,NPBE
      L=LPBE(LL)
      IF(ISPBE(LL).EQ.0)THEN
       IF(ISRED(L).EQ.1)THEN
        LR=LLRC(L)
        CWR(LR)=0.
        CCWR(LR)=0.
       ELSE
        LB=LLBC(L)
        CWB(LB)=0.
        CCWB(LB)=0.
       ENDIF
      ENDIF
      ENDDO
!
      DO LL=1,NPBS
      L=LPBS(LL)
      IF(ISPBS(LL).EQ.0)THEN
       IF(ISRED(L).EQ.1)THEN
        LR=LLRC(L)
        CNR(LR)=0.
        CCNR(LR)=0.
       ELSE
        LB=LLBC(L)
        CNB(LB)=0.
        CCNB(LB)=0.
       ENDIF
      ENDIF
      ENDDO
!
      DO LL=1,NPBN
      L=LPBN(LL)
      IF(ISPBN(LL).EQ.0)THEN
       IF(ISRED(L).EQ.1)THEN
        LR=LLRC(L)
        CSR(LR)=0.
        CCSR(LR)=0.
       ELSE
        LB=LLBC(L)
        CSB(LB)=0.
        CCSB(LB)=0.
       ENDIF
      ENDIF
      ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------C
!
      ELSE
!
!----------------------------------------------------------------------C
!
      DO L=2,LA
      LE=L+1
      LN=LNC(L)
      C1=-SPB(L)/(HRV(L)+HRU(L)+HRU(LE)+HRV(LN))
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(LE)
      CN(L)=C1*HRV(LN)
      CC(L)=C1
      ENDDO
!
      DO LL=1,NPBE
      L=LPBE(LL)-1
      LE=L+1
      LN=LNC(L)
      C1=(HRV(L)+HRU(L)+HRV(LN))
      IF(C1.EQ.0.) C1=1.
      C1=-SPB(L)/C1
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=0.
      CN(L)=C1*HRV(LN)
      CC(L)=C1
      ENDDO
!
      DO LL=1,NPBW
      L=LPBW(LL)+1
      LE=L+1
      LN=LNC(L)
      C1=(HRV(L)+HRU(LE)+HRV(LN))
      IF(C1.EQ.0.) C1=1.
      C1=-SPB(L)/C1
      CS(L)=C1*HRV(L)
      CW(L)=0.
      CE(L)=C1*HRU(LE)
      CN(L)=C1*HRV(LN)
      CC(L)=C1
      ENDDO
!
      DO LL=1,NPBN
      L=LPBN(LL)
      L=LSC(L)
      LE=L+1
      LN=LNC(L)
      C1=(HRV(L)+HRU(L)+HRU(LE))
      IF(C1.EQ.0.) C1=1.
      C1=-SPB(L)/C1
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(LE)
      CN(L)=0.
      CC(L)=C1
      ENDDO
!
      DO LL=1,NPBS
      L=LPBS(LL)
      L=LNC(L)
      LE=L+1
      LN=LNC(L)
      C1=(HRU(L)+HRU(LE)+HRV(LN))
      IF(C1.EQ.0.) C1=1.
      C1=-SPB(L)/C1
      CS(L)=0.
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(LE)
      CN(L)=C1*HRV(LN)
      CC(L)=C1
      ENDDO
!
!----------------------------------------------------------------------C
!
      ENDIF
!
!**********************************************************************C
!
! **  SMOOTH INITIAL SALINITY
!
!----------------------------------------------------------------------C
!
      IF(NSBMAX.GE.1)THEN
      CALL SALTSMTH
      ENDIF
!
!**********************************************************************C
!
! **  OUTPUT INITIAL DEPTH AND SALINITY FIELDS
!
!----------------------------------------------------------------------C
!
! **  PLOT SMOOTHED CELL CENTER STATIC DEPTHS
!
      DO L=2,LA
      PAM(L)=HMP(L)
      ENDDO
      WRITE (7,16)
      CALL PPLOT (2)
!
!----------------------------------------------------------------------C
!
      CALL DEPPLT
!
!----------------------------------------------------------------------C
!
! **  PLOT INITIAL SALINITY IN SURFACE AND BOTTOM LAYERS
!
      DO KK=1,KC,KS
      DO L=2,LA
      PAM(L)=SAL(L,KK)
      ENDDO
      WRITE (7,316) KK
      CALL PPLOT (1)
      ENDDO
!
!----------------------------------------------------------------------C
!
   16 FORMAT (1H1,' CELL CENTER STATIC DEPTHS',//)
  316 FORMAT (1H1,'INITIAL SALINITY IN LAYER',I5,//)
!
!----------------------------------------------------------------------C
!
!
!**********************************************************************C
!
! **  INITIALIZE SALINITY AND TEMPATURE DATA ASSIMILATION
!
      DO J=1,NLDAM
	DO I=1,NDDAM
        FSALASM(I,J)=0.0
	  FVOLASM(I,J)=0.0
        FTEMASM(I,J)=0.0
	ENDDO
	ENDDO
!
!**********************************************************************C
!
! **  INITIALIZE WATER QUALITY MODEL AND READ INPUT
!
      IF(ISTRAN(8).GE.1) CALL WQ3DINP
!
!**********************************************************************C
!
! **  INITIALIZE EFDC EXPLORER OUTPUT
!
      IF(ISSPH(8).EQ.1.OR.ISBEXP.EQ.1) CALL EEXPOUT(1)
!
!**********************************************************************C
!
! **  INITIALIZE BMD2 FILE                            !HNR_GHD 6/2022 BMD2
!
        CALL INIT_BMD2                                       !HNR_GHD 6/2022 BMD2
!
!**********************************************************************C
!
! **  SELECT FULL HYDRODYNAMIC AND MASS TRANSPORT CALCULATION OR
! **  LONG-TERM MASS TRANSPORT CALCULATION
!
      IF(ISLTMT.EQ.0)THEN
         IF(IS1DCHAN.EQ.0)THEN
           IF(IS2TIM.EQ.0) THEN
		   IF(IGRIDV.EQ.0) CALL HDMT
		   IF(IGRIDV.EQ.1) CALL HDMTGVC
           ENDIF
           IF(IS2TIM.GE.1) THEN
		   IF(IGRIDV.EQ.0) CALL HDMT2T
		   IF(IGRIDV.EQ.1) CALL HDMT2TGVC
           ENDIF
         ENDIF
         IF(IS1DCHAN.GE.1) CALL HDMT1D
      ENDIF
      IF(ISLTMT.GE.1) CALL LTMT
!
!**********************************************************************C
!
! **  OUTPUT TIMING INFORMATION
!
! **  EE TIMING INFO
!
      TIMEEND=SECNDS(TIMEEND)
!  TCPU=DTIME(CPUTIME)       !hnr
      CALL TIMELOG(N,TIMEDAY)
!
! **  END EE TIMING INFO
!
! **  NORMAL SCREEN AND EFDCLOG.OUT FILE OUTPUT
!
      WRITE(6,5989)
      WRITE(6,5995)THDMT,TCONG
      WRITE(6,5996)TPUV,TSSED
      WRITE(6,5997)TCEXP,TAVB
      WRITE(6,5998)TUVW,TQQQ
      WRITE(6,5999)TVDIF,TSADV
      WRITE(6,5994)TLRPD,TTBXY
      WRITE(6,5993)TWQDIF,TWQADV
      WRITE(6,5992)TWQKIN,TWQSED
      WRITE(6,5991)CPUTIME(1),CPUTIME(2)
      WRITE(6,5990)TIMEEND,TCPU
!
!      WRITE(8,5989)   !hnr 7/27/2009
!      WRITE(8,5995)THDMT,TCONG   !hnr 7/27/2009
!      WRITE(8,5996)TPUV,TSSED   !hnr 7/27/2009
!      WRITE(8,5997)TCEXP,TAVB   !hnr 7/27/2009
!      WRITE(8,5998)TUVW,TQQQ   !hnr 7/27/2009
!      WRITE(8,5999)TVDIF,TSADV   !hnr 7/27/2009
!      WRITE(8,5994)TLRPD,TTBXY   !hnr 7/27/2009
!      WRITE(8,5993)TWQDIF,TWQADV   !hnr 7/27/2009
!      WRITE(8,5992)TWQKIN,TWQSED   !hnr 7/27/2009
!      WRITE(8,5991)CPUTIME(1),CPUTIME(2)   !hnr 7/27/2009
!      WRITE(8,5990)TIMEEND,TCPU   !hnr 7/27/2009
!
!      WRITE(8,5989)   !hnr 7/27/2009
!      WRITE(8,*)THDMT,TCONG   !hnr 7/27/2009
!      WRITE(8,*)TPUV,TSSED   !hnr 7/27/2009
!      WRITE(8,*)TCEXP,TAVB   !hnr 7/27/2009
!      WRITE(8,*)TUVW,TQQQ   !hnr 7/27/2009
!      WRITE(8,*)TVDIF,TSADV   !hnr 7/27/2009
!      WRITE(8,*)TLRPD,TTBXY   !hnr 7/27/2009
!      WRITE(8,*)TWQDIF,TWQADV   !hnr 7/27/2009
!      WRITE(8,*)TWQKIN,TWQSED   !hnr 7/27/2009
!
      IF(ISCRAY.NE.0)THEN
        WRITE(6,5995)WTHDMT,WTCONG
        WRITE(6,5996)WTPUV,WTSSED
        WRITE(6,5997)WTCEXP,WTAVB
        WRITE(6,5998)WTUVW,WTQQQ
        WRITE(6,5999)WTVDIF,WTSADV
        WRITE(6,5994)WTLRPD,WTTBXY
        WRITE(6,5993)WTWQDIF,WTWQADV
        WRITE(6,5992)WTWQKIN,WTWQSED
!        WRITE(8,5995)WTHDMT,WTCONG   !hnr 7/27/2009
!        WRITE(8,5996)WTPUV,WTSSED   !hnr 7/27/2009
!        WRITE(8,5997)WTCEXP,WTAVB   !hnr 7/27/2009
!        WRITE(8,5998)WTUVW,WTQQQ   !hnr 7/27/2009
!        WRITE(8,5999)WTVDIF,WTSADV   !hnr 7/27/2009
!        WRITE(8,5994)WTLRPD,WTTBXY   !hnr 7/27/2009
!        WRITE(8,5993)WTWQDIF,WTWQADV   !hnr 7/27/2009
!        WRITE(8,5992)WTWQKIN,WTWQSED   !hnr 7/27/2009
	ENDIF
!
 5989 FORMAT(//' TIMING INFORMATION IN SECONDS'/)
 5990 FORMAT(' ELLAPSED TIME = ',F14.4,'  CPU TIME    = ',F14.4/)
 5991 FORMAT(' USER CPU      = ',F14.4,'  SYSTEM CPU  = ',F14.4/)
 5992 FORMAT(' TWQKIN = ',F14.4,'  TWQSED = ',F14.4/)
 5993 FORMAT(' TWQDIF = ',F14.4,'  TWQADV = ',F14.4/)
 5994 FORMAT(' TLRPD  = ',F14.4,'  TTBXY  = ',F14.4/)
 5995 FORMAT(' THDMT  = ',F14.4,'  TCONG  = ',F14.4/)
 5996 FORMAT(' TPUV   = ',F14.4,'  TSSED  = ',F14.4/)
 5997 FORMAT(' TCEXP  = ',F14.4,'  TAVB   = ',F14.4/)
 5998 FORMAT(' TUVW   = ',F14.4,'  TQQQ   = ',F14.4/)
 5999 FORMAT(' TVDIF  = ',F14.4,'  TSADV  = ',F14.4/)
!
! **  END NORMAL SCREEN AND EFDCLOG.OUT FILE OUTPUT
!
! **  TIMING OUTPUT REQUIRED BY EE
! **  PLEASE CONFINE EE RELATED CHANGES TO FOLLOWING BLOCK OF CODE
!
      TIMEEND=TIMEEND/3600.
	TCPU=TCPU/3600.
!
      WRITE(9,7995)WTHDMT,WTCONG
      WRITE(9,7996)WTPUV,WTCGRS
      WRITE(9,7997)WTCEXP,WTAVB
      WRITE(9,7998)WTUVW,WTQQQ
      WRITE(9,7999)WTVDIF,WTSADV
      WRITE(9,7994)WTLRPD
      WRITE(9,7993)WTWQDIF,WTWQADV
      WRITE(9,7992)WTWQKIN,WTWQSED
!
      WRITE(9,6995)THDMT/3600.,TCONG/3600.
      WRITE(9,6996)TPUV/3600.,TSSED/3600.
      WRITE(9,6997)TCEXP/3600.,TAVB/3600.
      WRITE(9,6998)TUVW/3600.,TQQQ/3600.
      WRITE(9,6999)TVDIF/3600.,TSADV/3600.
      WRITE(9,7000)TLRPD/3600.,TMISC/3600.
      WRITE(9,7001)CPUTIME(1)/3600.,CPUTIME(2)/3600.
      WRITE(9,7002) REAL(TIMEEND),TCPU
!
 6995 FORMAT(//'***TIMING (Hours)',/,'T HDMT       = ',F14.4,'  T CONG GRAD  = ',F14.4)
 6996 FORMAT('T P&UV VELS  = ',F14.4,'  T SSEDTOX    = ',F14.4)
 6997 FORMAT('T EXPLICIT   = ',F14.4,'  T C VERT V&D = ',F14.4)
 6998 FORMAT('T CALC UVW   = ',F14.4,'  T TURB QQQ   = ',F14.4)
 6999 FORMAT('T VERT DFUSN = ',F14.4,'  T ADV TRANSP = ',F14.4)
 7000 FORMAT('TLRPD        = ',F14.4,'  T MISC TIME  = ',F14.4)
 7001 FORMAT('CPU USER     = ',F14.4,'  CPU SYSTEM   = ',F14.4)
 7002 FORMAT('ELAPSED TIME = ',F14.4,'  CPU TIME     = ',F14.4)
!
 7992 FORMAT(' WTWQKIN = ',F14.4,'  WTWQSED = ',F14.4/)
 7993 FORMAT(' WTWQDIF = ',F14.4,'  WTWQADV = ',F14.4/)
 7994 FORMAT(' WTLRPD = ',F14.4/)
 7995 FORMAT(' WTHDMT = ',F14.4,'  WTCONG = ',F14.4/)
 7996 FORMAT(' WTPUV  = ',F14.4,'  WTCGRS = ',F14.4/)
 7997 FORMAT(' WTCEXP = ',F14.4,'  WTAVB  = ',F14.4/)
 7998 FORMAT(' WTUVW  = ',F14.4,'  WTQQQ  = ',F14.4/)
 7999 FORMAT(' WTVDIF = ',F14.4,'  WTSADV = ',F14.4/)
!
! **  PLEASE CONFINE EE RELATED CHANGES TO PRECEDING BLOCK OF CODE
! **  END TIMING OUTPUT REQUIRED BY EE
!
!**********************************************************************C
!
! **  CLOSE OUTPUT  FILES
!
!----------------------------------------------------------------------C
!

      IF(BMDFLAG.GT.0) THEN         !HNR_GHD 6/2022 BMD2 closing
        WRITE(*,*)'OPTIMIZING AND CLOSING BMD2 FILE'
        CALL BMD2ENDFRAMES(HANDLEBMD2,IERROR)
        IF(iError .gt. 0)Then
          Call bmd2getlasterror(80,errMessage)
          Write(*,*)errMessage
          Stop
        End IF
        CALL BMD2CLOSE(HANDLEBMD2,IERROR)
        IF(iError .gt. 0)Then
          Call bmd2getlasterror(80,errMessage)
          Write(*,*)errMessage
          Stop
        End IF
        WRITE(*,*)'BMD2 FILE SUCCESSFULLY CLOSED.'
      ENDIF                !HNR_GHD 7/2019 BMD2

      IF (ISWASP.EQ.8) THEN                                   !hnr
        write(6,*)'  '                                        !hnr
        CALL hlclose(Ihl_handle,ierror)                       !hnr
        if(ierror .gt. 0)then                                 !hnr
          call hlgetlasterror(errstring)                      !hnr
          write(6,6000) ierror, errstring,'CLOSING'           !hnr
6000      format('Error ',I10, ' : ', A30,A10)                !hnr
          stop                                                !hnr
        end if                                                !hnr
      END IF                                                  !hnr

      CLOSE(7)
!      CLOSE(8)   !hnr 7/27/2009
      CLOSE(9)
!
!**********************************************************************C
!
      STOP
      END
