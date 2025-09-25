!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE SURFPLTE
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
! **  SUBROUTINE SURFPLT WRITES FILES TO CONTOUR FREE SURFACE
! **  ELEVATION
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      integer IVER
      CHARACTER*80 TITLE
!
!**********************************************************************C
!
! **  OUTPUT EFDC EXPLORER FORMAT
!
!      IF(IPPHXY.EQ.3)THEN
!
      LINES=LA-1
      IF(JSPPH.EQ.1)THEN
        OPEN(92,FILE='SURFCON.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
!######################################################################
!     HQI Change for EFDC Explorer compatibility
!    .        FORM='UNFORMATTED')
 
!######################################################################
        CLOSE(92,STATUS='DELETE')
        OPEN(92,FILE='SURFCON.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
!    .        FORM='UNFORMATTED')
        
!######################################################################
!     HQI Change for EFDC Explorer compatibility     
        IVER=101
!       IVER=102
!######################################################################

        WRITE(92)IVER,IC,JC,LINES
        JSPPH=0
!        CLOSE(92)
      ENDIF
!
!      OPEN(92,FILE='SURFCON.OUT',POSITION='APPEND'
!     &             ,FORM='UNFORMATTED')
!
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
!
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
      ELSE
        DELT=DTDYN
      ENDIF
      WRITE (92)N,TIME,DELT
!
      DO L=2,LA
        WRITE(92)HP(L)
      ENDDO
!
!      CLOSE(92)
!
!      ENDIF
!
!**********************************************************************C
!
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,9E14.5)
  201 FORMAT(9E14.5)
  250 FORMAT(12E12.4)
!MRM  200 FORMAT(2I5,1X,1P,8E13.5)
!MRM  250 FORMAT(1P,12E11.3)
!
!**********************************************************************C
!
      RETURN
      END
