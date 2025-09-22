!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
!      SUBROUTINE WWQTS(TINDAY)
      SUBROUTINE WWQTS
!
!**********************************************************************C
!
! **  SHEN'S MODIFICATION TO OUTPUT MACROALGAE
!
! **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
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
! hugo rodriguez  10/2010  algae chla/C time variable in algae kinetic file
!
!**********************************************************************C
!
! WRITE TIME-SERIES OUTPUT: WQCHLX=1/WQCHLX
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
      DIMENSION WQVOUT(NWQVM+6)
!
!      TINDAY=TINDAY
!
      OPEN(1,FILE='WQWCTS.OUT',STATUS='UNKNOWN',POSITION='APPEND')
!
      IF(ISDYNSTP.EQ.0)THEN
        TIMTMP=DT*FLOAT(N)+TCON*TBEGIN
        TIMTMP=TIMTMP/TCTMSR  
      ELSE
        TIMTMP=TIMESEC/TCTMSR  
      ENDIF
!
      DO M=1,IWQTS
        DO K=1,KC
          LL=LWQTS(M)
          IZ=IWQZMAP(LL,K)
!          DO NW=1,NWQV
          DO NW=1,NWQV+1    ! Sen Bai
            WQVO(LL,K,NW) = WQVO(LL,K,NW)*0.5
          ENDDO
          NWQOUT=0
          IF(ISTRWQ(1).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(LL,K,1)*WQCHLC(iz)   !hnr 10/2010
          ENDIF
          IF(ISTRWQ(2).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(LL,K,2)*WQCHLD(iz)   !hnr 10/2010
          ENDIF
          IF(ISTRWQ(3).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(LL,K,3)*WQCHLG(iz)   !hnr 10/2010
          ENDIF
!          DO NW=4,NWQV
           DO NW=4,NWQV    ! HAMRICK CORRECTED   +1  ! Sen Bai 
            IF(ISTRWQ(NW).EQ.1)THEN
              NWQOUT=NWQOUT+1
              WQVOUT(NWQOUT)=WQVO(LL,K,NW)
            ENDIF
          ENDDO
! Sen Bai
          IF(IDNOTRVA.GT.0)THEN
	      NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(LL,K,IDNOTRVA)
          ENDIF
! End of Sen Bai

! ** ADD TEMPERATURE AND SALINITY SATUARATION
          NWQOUT=NWQOUT+1
          WQVOUT(NWQOUT)=TEM(LL,K)
          NWQOUT=NWQOUT+1
          WQVOUT(NWQOUT)=SAL(LL,K)
! ** ADD DO SATUARATION
          NWQOUT=NWQOUT+1
          TVAL1=1./(TEM(LL,K)+273.15)
          TVAL2=TVAL1*TVAL1
          TVAL3=TVAL1*TVAL2
          TVAL4=TVAL2*TVAL2
          RLNSAT1=-139.3441+(1.575701E+5*TVAL1)-(6.642308E+7*TVAL2)+(1.2438E+10*TVAL3)            &
                  -(8.621949E+11*TVAL4)
          RLNSAT2=RLNSAT1-SAL(LL,K)*( 1.7674E-2-(1.0754E+1*TVAL1)+(2.1407E+3*TVAL2) )
          WQVOUT(NWQOUT) = EXP(RLNSAT2)
! ** ADD FLOW INDUCED REAERATION
          NWQOUT=NWQOUT+1
          WQVOUT(NWQOUT) = 0.
          IF(IWQKA(IZ) .EQ. 2)THEN
            UMRM = 0.5*( U(LL,K) + U(LL+1   ,K) )
            VMRM = 0.5*( V(LL,K) + V(LNC(LL),K) )
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
            WQVOUT(NWQOUT) = WQKRO(IZ) * XMRM**0.5 / HP(LL)**0.5
          ENDIF
! ** ADD WIND INDUCED REAERATION
          NWQOUT=NWQOUT+1
          WINDREA = WINDST(LL)
          WQVOUT(NWQOUT)=0.728*SQRT(WINDREA)+(0.0372*WINDREA-0.317)*WINDREA
          WRITE(1,71) IL(LL),JL(LL),K,TIMTMP,(WQVOUT(NWOUT),NWOUT=1,NWQOUT)
! PUT BACK IN BEGINNING FORM
          DO NW=1,NWQV+1    ! Sen Bai
            WQVO(LL,K,NW) = WQVO(LL,K,NW)*2.0
          ENDDO
!
        ENDDO
      ENDDO
!
   71 FORMAT(3I5,F11.5,1X, 25E11.3)
!
!     K=1
!OLD      DO M=1,IWQTS
!OLD         DO K=1,KC
!OLD          LL=LWQTS(M)
!OLD          DO NW=1,NWQV
!OLD            WQVO(LL,K,NW) = WQVO(LL,K,NW)*0.5
!OLD          ENDDO
!OLD          CHLWQ = WQVO(LL,K,1)*WQCHLC + WQVO(LL,K,2)*WQCHLD
!OLD     *      + WQVO(LL,K,3)*WQCHLG
!OLD          CHLC = WQVO(LL,K,1)*WQCHLC
!OLD          CHLD = WQVO(LL,K,2)*WQCHLD
!OLDC          CHLCD = CHLC + CHLD
!OLD          CHLG = WQVO(LL,K,3)*WQCHLG
!OLD          TBWQ = WQVO(LL,K,1)+WQVO(LL,K,2)+WQVO(LL,K,3)
!OLD          TOCWQ = WQVO(LL,K,4)+WQVO(LL,K,5)+WQVO(LL,K,6) + TBWQ
!OLD          IF(IWQSRP.EQ.1)THEN
!OLD            O2WQ = MAX(WQVO(LL,K,19), 0.0)
!OLD            TAMDWQ = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQVO(LL,K,20) )
!OLD            TAMPWQ = WQVO(LL,K,20) - TAMDWQ
!OLD            PO4DWQ = WQVO(LL,K,10) / (1.0 + WQKPO4P*TAMPWQ)
!OLD            SADWQ = WQVO(LL,K,17) / (1.0 + WQKSAP*TAMPWQ)
!OLD           ELSE IF(IWQSRP.EQ.2)THEN
!OLD            PO4DWQ = WQVO(LL,K,10) / (1.0 + WQKPO4P*SEDT(LL,K))
!OLD            SADWQ = WQVO(LL,K,17) / (1.0 + WQKSAP*SEDT(LL,K))
!OLD           ELSE
!OLD            PO4DWQ = WQVO(LL,K,10)
!OLD            SADWQ  = WQVO(LL,K,17)
!OLD          ENDIF
!OLD          XPO4DWQ = MAX(PO4DWQ,0.0)
!OLD          APCWQ = 1.0 / (WQCP1PRM + WQCP2PRM*EXP(-WQCP3PRM*XPO4DWQ))
!OLD          TPWQ = WQVO(LL,K,7)+WQVO(LL,K,8)+WQVO(LL,K,9)+WQVO(LL,K,10)
!OLD     *      + APCWQ*TBWQ
!OLD          TNWQ = WQVO(LL,K,11)+WQVO(LL,K,12)+WQVO(LL,K,13)+WQVO(LL,K,14)
!OLD     *      +WQVO(LL,K,15) + WQANCC*WQVO(LL,K,1)+WQANCD*WQVO(LL,K,2)
!OLD     *      +WQANCG*WQVO(LL,K,3)
!OLD          IF(IWQSI.EQ.1)THEN
!OLD            TSIWQ = WQVO(LL,K,16)+WQVO(LL,K,17) + WQASCD*WQVO(LL,K,2)
!OLD           ELSE
!OLD            TSIWQ = 0.0
!OLD          ENDIF
!
! M. MORTON 03/11/96 ADDED BOD5 TERM FOR OUTPUT:
!OLD          IZ=IWQZMAP(LL,K)
!OLD          BOD5 = 2.67 * ( WQVO(LL,K,5)*(1.0-EXP(-5.0*WQKLC))
!OLD     +           + WQVO(LL,K,6)*(1.0-EXP(-5.0*WQKDC(IZ)))
!OLD     +           + WQVO(LL,K,18)*(1.0-EXP(-5.0*WQKCD(IZ)))
!OLD     +           + WQVO(LL,K,1)*(1.0-EXP(-5.0*WQBMRC(IZ)))
!OLD     +           + WQVO(LL,K,2)*(1.0-EXP(-5.0*WQBMRD(IZ)))
!OLD     +           + WQVO(LL,K,3)*(1.0-EXP(-5.0*WQBMRG(IZ))) )
!OLD     +           + 4.57 * WQVO(LL,K,14)*(1.0-EXP(-5.0*WQNITM))
!OLD
!MRM          BOD5=0.0       ! JI, 10/2/97
! MRM MODIFIED OUTPUT FORMAT: REMOVED TAM AND TMP AND ADDED
! DIATOMS (CHLD) AND GREEN ALGAE (CHLG); REMOVED APCWQ AND
! ADDED CYANOBACTERIA (CHLC):
!OLD          IF(IDNOTRVA.EQ.0)THEN
!OLD          WRITE(1,71) IL(LL),JL(LL),K,TIMTMP, CHLWQ,TOCWQ,WQVO(LL,K,6),
!OLD     +      TPWQ, WQVO(LL,K,9), WQVO(LL,K,10), PO4DWQ, CHLC, TNWQ,
!OLD     +      (WQVO(LL,K,NW),NW=13,15), TSIWQ, (WQVO(LL,K,NW),NW=16,17),
!OLD     +      SADWQ, BOD5, WQVO(LL,K,19), CHLD, CHLG, WQVO(LL,K,NWQV)
!OLD         ELSE
! HHTMP = WATER DEPTH (METERS)
!          HHTMP = GI*P(LL) - BELV(LL)
! CHLM  = MACROALGAE BIOMASS IN MICROGRAMS/SQUARE METER:
!          CHLM = WQVO(LL,K,NWQV+1) * WQCHLM * HHTMP
! CHLM IN UG/L AS FOLLOWS:
!OLD          CHLM = WQVO(LL,K,NWQV+1) * WQCHLM
!OLD          WRITE(1,71) IL(LL),JL(LL),K,TIMTMP, CHLWQ,TOCWQ,WQVO(LL,K,6),
!OLD     +      TPWQ, WQVO(LL,K,9), WQVO(LL,K,10), PO4DWQ, CHLC, TNWQ,
!OLD     +      (WQVO(LL,K,NW),NW=13,15), TSIWQ, (WQVO(LL,K,NW),NW=16,17),
!OLD     +      SADWQ, BOD5, WQVO(LL,K,19), CHLD, CHLG,WQVO(LL,K,NWQV),
!OLD     +      CHLM
! MRM     +      WQVO(LL,K,NWQV+1)
!OLD          ENDIF
!          WRITE(1,71) IL(LL),JL(LL),K, TIMTMP,CHLWQ,TOCWQ,WQVO(LL,K,6),
!     *      TPWQ,WQVO(LL,K,9),WQVO(LL,K,10),PO4DWQ,APCWQ,TNWQ,
!     *      (WQVO(LL,K,NW),NW=13,15),TSIWQ,(WQVO(LL,K,NW),NW=16,17),
!     *      SADWQ,(WQVO(LL,K,NW),NW=18,20),TAMPWQ,WQVO(LL,K,NWQV)
!OLD         ENDDO
!OLD      ENDDO
!
      CLOSE(1)
!
!   71 FORMAT(3I5,F11.5, 1P, 23E11.3)
!
      RETURN
      END