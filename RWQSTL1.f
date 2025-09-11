C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQSTL1(timtmp)
C
C**********************************************************************C
C
C **  subroutine to introduce Morton's bentic time varying format 
c          for the settling velocities
C
C **  created on October 2010 by Hugo Rodriguez 
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C
C----------------------------------------------------------------------C
C
C
C**********************************************************************C
C
C READ IN SPATIALLY AND/OR TEMPORALLY VARYING PARAMETERS FOR SETTLING
C VELOCITIES OF ALGAE, RPOM, LPOM & PARTICULATE METAL (UNIT INWQSTL).
C ALSO SPATIALLY/TEMPORALLY VARYING REAERATION ADJUSTMENT FACTOR.
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      CHARACTER TITLE(3)*79
C
      OPEN(1,FILE=STLFN,STATUS='UNKNOWN')
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
 
     
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
        read(1,*)
        read(1,*)
      WRITE(2,60)'* SETTLING VELOCITY AT  ', timtmp,
     *  ' of model run'
C     SEQUENTIALLY READ THROUGH settling velocity FILE UNTIL THE APPROPRIATE
C     TIME IS FOUND:
C
C      sDAY   = CURRENT DAY AT WHICH settling IS IN EFFECT
C      stlDAY = NEXT DAY AT WHICH settling CHANGES (PASSED TO MAIN PROGRAM)
C

10    READ(1, *, END=15) stlDAY
      IF(stlDAY .GT. TIMTMP) GOTO 20
      sDAY = stlDAY

      DO I=1,IWQZ
        READ(1,*,end=15) MM,WQWSC(mm),WQWSD(mm),WQWSG(mm),WQWSRP(mm),
     *    WQWSLP(mm),WQWSS(mm)
      ENDDO
      goto 10

C     UNEXPECTED END-OF-FILE ENCOUNTERED:
15    WRITE(2,16) stlFN
16    FORMAT(//,' ************* WARNING *************',/,
     +          ' END-OF-FILE ENCOUNTERED IN FILE: ', A20,/,/
     +          ' settling velocities SET TO VALUES CORRESPONDING ',
     +          ' TO LAST DAY IN FILE.',/)

20    CONTINUE

      WRITE(2, 48) sDAY
48    FORMAT(/,' DAY IN settling velocity FILE: ',F10.5,/,
     + '    IZ    WSc      WSd     WSg    WSrp    WSlp'
     + '    WSs    WSM REAERF')
C
      DO I=1,IWQZ
        write(2,51) i,WQWSC(i),WQWSD(i),WQWSG(i),WQWSRP(i),
     *    WQWSLP(i),WQWSS(i)
      ENDDO

      CLOSE(2)
      CLOSE(1)
C
   50 FORMAT(A79)
   51 FORMAT(I8, 10F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, f10.5, A13)
  999 FORMAT(1X)
C
      RETURN
      END
