!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE OUTPUT2
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
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
! **  OUTPUT RESULTS OF RELAXATION SOLUTION
!
      WRITE (7,40) RP
   40 FORMAT (1H1,' RESULTS OF RELAX SOLUTION - RP=',F5.2,//)
!
      WRITE (7,41)
   41 FORMAT (' GLOBAL SQUARED ERROR',//)
!
      WRITE(7,43)ERRMAX,ERRMIN
   43 FORMAT('ERRMAX =',3X,E12.4,5X,'ERRMIN =',3X,E12.4)
!
!     DO N=1,NTS,10
!     WRITE (7,20) N, (ERR(NN),NN=N,N+9)
!     ENDDO
   20 FORMAT (1X,I5,3X,10E12.4)
!
      WRITE(7,40)RP
      WRITE (7,42)
   42 FORMAT (' ITERATIONS TO CONVERGENCE',//)
!
      WRITE(7,44)ITRMAX,ITRMIN
   44 FORMAT('ITRMAX =',I5,5X,'ITRMIN =',I5)
!
!     DO N=1,NTS,10
!     WRITE (7,21) N, (ITR(NN),NN=N,N+9)
!     ENDDO
   21 FORMAT (1X,I5,5X,10I10)
!
   30 FORMAT (10E12.4)
!
!**********************************************************************C
!
! **  OUTPUT HARMONIC ANALYSIS
!
!----------------------------------------------------------------------C
!
      DO L=2,LA
      PAM(L)=(AMCP(L)*AMCP(L)+AMSP(L)*AMSP(L))**.5
      IF(AMSP(L).EQ.0.0.AND.AMCP(L).EQ.0.0)THEN
       PPH(L)=999999.
      ELSE
       PPH(L)=ATAN2(AMSP(L),AMCP(L))
      ENDIF
      ENDDO
!
!----------------------------------------------------------------------C
!
      DO L=2,LA
      PAM(L)=PAM(L)*GI
      ENDDO
      WRITE (7,55)
      CALL PPLOT (1)
!
   55 FORMAT (1H1,'TIDAL SURFACE DISPLACEMENT AMPLITUDE IN METERS',//)
!
!----------------------------------------------------------------------C
!
      DO L=2,LA
      PAM(L)=0.5*TIDALP*PPH(L)/PI
      ENDDO
      WRITE(7,588)
      CALL PPLOT (1)
!
  588 FORMAT (1H1,'TIDAL SURFACE DISPLACEMENT PHASE IN SEC',//)
!
!----------------------------------------------------------------------C
!
!     DO KK=1,KC,KS
!
!     DO L=2,LA
!     PAM(L)=SQRT(AMCU(L,KK)*AMCU(L,KK)
!    &               +AMSU(L,KK)*AMSU(L,KK))
!     ENDDO
!     WRITE (7,56) KK
!     CALL PPLOT (2)
!
!     ENDDO
!
!  56 FORMAT (1H1,' X TIDAL VELOCITY AMPLITUDE, M/S, LAYER',I5,//)
!
!----------------------------------------------------------------------C
!
!     DO KK=1,KC,KS
!
!     DO L=2,LA
!     LN=LNC(L)
!     PAM(L)=SQRT(AMCV(L,KK)*AMCV(L,KK)
!    &               +AMSV(L,KK)*AMSV(L,KK))
!     ENDDO
!     WRITE (7,57) KK
!     CALL PPLOT (2)
!
!     ENDDO
!
!  57 FORMAT (1H1,' Y TIDAL VELOCITY AMPLITUDE, M/S, LAYER',I5,//)
!
!**********************************************************************C
!
! **  PRINTED OUTPUT OF P,U,AND V AMPLITUDES
!
!     WRITE(7,71)
!  71 FORMAT(1H1,'   L   IL   JL',9X,'PA',12X,'PP',13X,'U',14X,'V',//)
!
!     DO L=2,LA
!     WRITE(7,72)L,IL(L),JL(L),PAM(L),PPH(L),U(L),V(L)
!     ENDDO
   72 FORMAT(3I5,4(3X,E12.4))
!
!**********************************************************************C
!
! **  OUTPUT VECTOR POTENTIAL TRANSPORT VELOCITY
!
!     IF(ISVPTHA.NE.1) GOTO 100
!
!----------------------------------------------------------------------C
!
!     DO KK=1,KC,KS
!
!     DO L=2,LA
!     PAM(L)=0.5*(UVPT(L,KK)+UVPT(L+1,KK))/HMU(L)
!     ENDDO
!     WRITE (7,1458) KK
!     CALL PPLOT (2)
!
!     ENDDO
!
 1458 FORMAT(1H1,' X VECTOR POTENTIAL TRANSPORT VEL, M/S, LAYER',I5,//)
!
!----------------------------------------------------------------------C
!
!     DO KK=1,KC,KS
!
!     DO L=2,LA
!     LN=LNC(L)
!     PAM(L)=0.5*(VVPT(L,KK)+VVPT(LN,KK))/HMV(L)
!     ENDDO
!     WRITE (7,1459) KK
!     CALL PPLOT (2)
!
!     ENDDO
!
 1459 FORMAT(1H1,' Y VECTOR POTENTIAL TRANSPORT VEL, M/S, LAYER',I5,//)
!
!**********************************************************************C

  100 CONTINUE
!
      RETURN
      END
