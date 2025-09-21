!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      REAL FUNCTION CSEDTAUS(DENBULK,TAUCO,VDRO,VDR,VDRC,IOPT,L)
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
      INCLUDE 'EFDC.PAR'
!
!#######################################################################
!  HQI Added, 11/18/2003, HAMRICK COMMENTED OUT SINCE NOT NEEDED FOR
!  FOR 11/24 VERSION OF IOPT= 99 THAT IS ACTIVE AS OF 01/08/2004
!
!JH      INCLUDE 'EFDC.CMN'
!JH      REAL rSEDPHI(4),rSIGPHI,D50SIG,D90SIG,rSNDBT,cohcon,rfracoh
!JH      REAL rTVAR3W,rTVAR3E,Z90
!JH      real rfrac1,rfrac2,rfrac3,ubnd,lbnd,usiz,lsiz
!
!#######################################################################
!
! **  CALCULATES CRITIAL STRESS FOR SURFACE EROSION OF COHESIVE
! **  SEDIMENT AS A FUNCTION OF BED BULK DENSITY
!
! **  IOPT=1  BASED ON
! **
! **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT ERODIBILITY
! **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
! **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
!
! **  IOPT=2 & 3  BASED ON J. M. HAMRICK'S MODIFICATION OF
! **
! **  SANFORD, L.P., AND J. P. Y. MAA, 2001: A UNIFIED EROSION FORMULATION
! **  FOR FINE SEDIMENT, MARINE GEOLOGY, 179, 9-23.
!
! **  IOPT=4  BASED ON J. M. HAMRICK'S PARAMETERIZATION OF SEDFLUME
!     TEST DATA
! **
!
      IF(IOPT.EQ.1)THEN
        DENBULK=0.001*DENBULK
        IF(DENBULK.LE.1.065)THEN
          CSEDTAUS=0.0
        ELSE
          TMP=(DENBULK-1.065)**0.2
          CSEDTAUS=0.001*(0.883*TMP+0.05)
        ENDIF
      ENDIF
!
      IF(IOPT.EQ.2)THEN
        CSEDTAUS=TAUCO*(1.+VDRO)/(1.+VDR)
      ENDIF
!
      IF(IOPT.EQ.3)THEN
        CSEDTAUS=TAUCO*(1.+VDRO)/(1.+VDRC)
      ENDIF
!
!      IF(IOPT.EQ.3)THEN
!        CSEDTAUS=TAUCO*EXP((1.+VDRO)/(1.+VDRC))
!      ENDIF
!
      IF(IOPT.EQ.4)THEN
!        CSEDTAUS=TAUCO*EXP((1.+VDRO)/(1.+VDRC))
        CSEDTAUS=TAUCO
      ENDIF
!
      IF(IOPT.EQ.5)THEN
!        CSEDTAUS=TAUCO*EXP((1.+VDRO)/(1.+VDRC))
        CSEDTAUS=TAUCO
      ENDIF
!

!#######################################################################
!  HQI Change, 08/25/03, and 11/24/03  SO and RM
!  Change to implement critical shear stress option
!  CSEDTAUS is tau/rho with tau in dyne/cm^2
!  IWRSP(1) = 99 for shear stress as a function of bulk density from
!                sed-flume experiments
      IF(IOPT.GE.99)THEN
!     Bulk density in g/cc and Tau in N/m**2
         IF(L.LE.265) then
!     Woods Pond cells
!     CSEDTAUS=0.015*((DENBULK/1000.)**13.1)/1000
!     CSEDTAUS=1.2/1000.
            CSEDTAUS=0.2/1000.
         ELSE
!     North of Woods Pond cells
!     CSEDTAUS=8.0E-6*((DENBULK/1000.)**17.4)/1000.
!     CSEDTAUS=0.6/1000.
            CSEDTAUS=0.4/1000.
         ENDIF
      ENDIF
!
!#######################################################################
!
!#######################################################################
!  HQI Added, 11/18/2003
!  Compute the D90 of a cell bed from the fractions of each grain size
!  and the computed D50's (i.e. the Deff for the four classes.
!
!      IF(IOPT.EQ.9999)THEN
!      KTOP=KBT(L)
!      SEDDIA(1)=22.0*1.e-6 !convert micron to meter
!
!      rfracoh= SEDB(L,KTOP,1)/( SEDB(L,KTOP,1)+SNDB(L,KTOP,1)
!     &                          +SNDB(L,KTOP,2)+SNDB(L,KTOP,3) )
!      rfrac1 = rfracoh + SNDB(L,KTOP,1)
!     &  /( SEDB(L,KTOP,1)+SNDB(L,KTOP,1)+SNDB(L,KTOP,2)+SNDB(L,KTOP,3) )
!      rfrac2 = rfrac1 + SNDB(L,KTOP,2)
!     &  /( SEDB(L,KTOP,1)+SNDB(L,KTOP,1)+SNDB(L,KTOP,2)+SNDB(L,KTOP,3) )
!      rfrac3 = rfrac2 + SNDB(L,KTOP,3)
!     &  /( SEDB(L,KTOP,1)+SNDB(L,KTOP,1)+SNDB(L,KTOP,2)+SNDB(L,KTOP,3) )
!
!      ubnd = 999.
!      lbnd = -999.
!
! 6400, 570, 160, 63 in microns
!
!      if ( rfrac3 .GT. 0.9 .AND. rfrac3 .LT. ubnd ) THEN
!        ubnd = rfrac3
!        usiz = SEDDIA(4)
!      END if
!      if ( rfrac2 .GT. 0.9 .AND. rfrac2 .LT. ubnd ) THEN
!        ubnd = rfrac2
!        usiz = SEDDIA(3)
!      END if
!      if ( rfrac1 .GT. 0.9 .AND. rfrac1 .LT. ubnd ) THEN
!        ubnd = rfrac1
!        usiz = SEDDIA(2)
!      END if
!      if ( rfracoh .GT. 0.9 .AND. rfracoh .LT. ubnd ) THEN
!        ubnd = rfracoh
!        usiz = SEDDIA(1)
!      END if
!
!      if ( rfracoh .LT. 0.9 .AND. rfracoh .GT. lbnd ) THEN
!        lbnd = rfracoh
!        lsiz = SEDDIA(1)
!      END if
!      if ( rfrac1 .LT. 0.9 .AND. rfrac1 .GT. lbnd ) THEN
!        lbnd = rfrac1
!        lsiz = SEDDIA(2)
!      END if
!      if ( rfrac2 .LT. 0.9 .AND. rfrac2 .GT. lbnd ) THEN
!        lbnd = rfrac2
!        lsiz = SEDDIA(3)
!      END if
!      if ( rfrac3 .LT. 0.9 .AND. rfrac3 .GT. lbnd ) THEN
!        lbnd = rfrac3
!        lsiz = SEDDIA(4)
!      END if
!
!      IF ( lbnd .LT. 0. ) THEN
!        lbnd = 0.
!        lsize = 0.
!      END IF
!      D90SIG = lsiz + (usiz-lsiz)*(0.9-lbnd)/(ubnd-lbnd)
!
!
!
!  SEDPHIC ** GP-CONVERT SEDIMENT DIAMETERS IN M TO MM AND SET PHI SIZE
!
!      SEDDIA(1)=22.0*1.e-6 !convert micron to meter
!      DO NS=1,NSED+NSND
!        rSEDPHI(NS)=-LOG(1000.*SEDDIA(NS))/LOG(2.)
!      END DO
!
! **  GP - SET MEAN PHI FOR TOP LAYER OF BED
!
!      rTVAR3W=0.
!      rTVAR3E=0.
!      rSIGPHI=0.
!
!      DO NX=1,NSND
!        NS=NX+NSED
!        rTVAR3W=rTVAR3W+rSEDPHI(NS)*SNDB(L,KTOP,NX)
!        rTVAR3E=rTVAR3E+SNDB(L,KTOP,NX)
!      END DO
!      rTVAR3W=rTVAR3W+rSEDPHI(1)*SEDB(L,KTOP,1)
!      rTVAR3E=rTVAR3E+SEDB(L,KTOP,1)
!
!      IF (rTVAR3E.LE.0.) THEN
!        rTVAR3E=1.
!      ELSE
!        rTVAR3W=rTVAR3W/rTVAR3E
!      END IF
!
!      DO NX=1,NSND
!        NS=NX+NSED
!        rSIGPHI = rSIGPHI + ((rSEDPHI(NS)-rTVAR3W)**2)
!     &                     *SNDB(L,KTOP,NX)
!      END DO
!      rSIGPHI = rSIGPHI + ((rSEDPHI(1)-rTVAR3W)**2)
!     &                     *SEDB(L,KTOP,1)
!      rSIGPHI = rSIGPHI/rTVAR3E
!     write(6,*) rSIGPHI, rTVAR3E
!     stop
!
!      rSIGPHI=SQRT(rSIGPHI)
!      rSIGPHI=2.**(rSIGPHI)
!
!
! **  SET MEAN D50
!
!      D50SIG=0.
!      rSNDBT=0.
!
!      DO NX=1,NSND
!       NS=NX+NSED
!       D50SIG=D50SIG+SNDB(L,KTOP,NX)*LOG(SEDDIA(NS))
!       !D50SIG=D50SIG+SNDB(L,KTOP,NX)*(SEDDIA(NS))
!       rSNDBT=rSNDBT+SNDB(L,KTOP,NX)
!      ENDDO
!      D50SIG=D50SIG+SEDB(L,KTOP,1)*LOG(SEDDIA(1))
!      !D50SIG=D50SIG+SEDB(L,KTOP,1)*(SEDDIA(1))
!      rSNDBT=rSNDBT+SEDB(L,KTOP,1)
!
!      D50SIG=EXP(D50SIG/rSNDBT)
!      !D50SIG=D50SIG/rSNDBT
!
!
! ** Compute the D90 from the standard deviation of grain size
!    distributionin bed
!
!
!      Z90 = 1.281551  !(Z-score for the 90th percentile)
!     D90SIG= (((D50SIG*1000.) + (Z90*rSIGPHI)))/1000.
!
! ** Cohesive concentration
!      cohcon= (SEDB(L,KTOP,1)*1e-6)/HBED(L,KTOP)
!
!      rfracoh = rfracoh*100.
!
!      CSEDTAUS=(0.36*((D90SIG/D50SIG)**0.948803))
!     &         *(0.613976*(rfracoh**(-0.109768)))
!     &         *(0.846585*(cohcon**0.217732))
!
!      CSEDTAUS = CSEDTAUS/1000.
!     Limit tau critical to 2 to 50 dynes/cm^2 based on Earl  Hayter
!     recommendation  09/12/03
!      if (CSEDTAUS.LT.(2./10000.)) then
!        CSEDTAUS = 2./10000.
!      end if
!
!     !write(6,'(f12.5)') D50SIG, D90SIG, rfracoh, COHCON, CSEDTAUS
!     write(6,'(f12.5)')CSEDTAUS
!     IF ( L .EQ. 4393 ) THEN
!       stop
!     END IF
!
!      ENDIF
!
!#######################################################################

      RETURN
      END