!
!**********************************************************************C
!**********************************************************************C
!**********************************************************************C
!
      SUBROUTINE LTMT
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
! **  SUBROUTINE LTMT EXECUTES A LONG-TERM MASS TRANSPORT
! **  TIME INTEGRATION 
!
!**********************************************************************C
!
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
!
!**********************************************************************C
!
      OPEN(99,FILE='RESTRAN.INP',STATUS='UNKNOWN')
!
!**********************************************************************C
!
! **  READ FIRST MEAN MASS TRANSPORT FIELD AND INITIALIZED NEXT
!
!----------------------------------------------------------------------C
!
      CALL RESTRAN
!
      DO L=2,LA
      H1P(L)=HLPF(L)
      HP(L)=HLPF(L)
      UHDY1E(L)=0.
      UHDYE(L)=0.
      VHDX1E(L)=0.
      VHDXE(L)=0.
      ENDDO
!
!     IF(ISLTMT.EQ.1)THEN
        DO K=1,KS
        DO L=2,LA
        AB(L,K)=ABLPF(L,K)
        ENDDO
        ENDDO
!      ELSE
!       DO K=1,KS
!       DO L=2,LA
!       AB(L,K)=ABEFF(L,K)
!       ENDDO
!       ENDDO
!     ENDIF
!
      DO K=1,KC
      DO L=1,LC
      UHDY1(L,K)=UHLPF(L,K)*DYU(L)
      VHDX1(L,K)=VHLPF(L,K)*DXV(L)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=1,LC
      UHDY(L,K)=UHDY1(L,K)
      UHDY1E(L)=UHDY1E(L)+UHDY1(L,K)*DZC(K)
      VHDX(L,K)=VHDX1(L,K)
      VHDX1E(L)=VHDX1E(L)+VHDX1(L,K)*DZC(K)
      ENDDO
      ENDDO
!
      DO L=1,LC
      UHDYE(L)=UHDY1E(L)
      VHDXE(L)=VHDX1E(L)
      ENDDO
!
!----------------------------------------------------------------------C
!
      IF(ISLTMTS.EQ.0)THEN
!
      DO K=1,KC
      DO L=1,LC
      SAL(L,K)=0.
      SAL1(L,K)=0.
      ENDDO
      ENDDO
!
      ELSE
!
      DO K=1,KC
      DO L=1,LC
      SAL(L,K)=SALLPF(L,K)
      SAL1(L,K)=SALLPF(L,K)
      ENDDO
      ENDDO
!
      ENDIF
!
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      SAL1(L,K)=SALLPF(L,K)
      SAL(L,K)=SALLPF(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      SAL1(L,K)=SALLPF(L,K)
      SAL(L,K)=SALLPF(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      SAL1(L,K)=SALLPF(L,K)
      SAL(L,K)=SALLPF(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      SAL1(L,K)=SALLPF(L,K)
      SAL(L,K)=SALLPF(L,K)
      ENDDO
      ENDDO
!
!**********************************************************************C
!
! **  READ SECOND MEAN MASS TRANSPORT FIELD FOR UNSTEADY CASE 
!
!----------------------------------------------------------------------C
!
      IF(ISSSMMT.EQ.0)THEN
!
      CALL RESTRAN
!
      RNTCMMT=FLOAT(NTSMMT)/FLOAT(NTSPTC)
      DO L=2,LA
      HP(L)=H1P(L)+(HLPF(L)-H1P(L))/(RNTCMMT*FLOAT(NTSPTC))
      UHDYE(L)=0.
      VHDXE(L)=0.
      ENDDO
!
!     IF(ISLTMT.EQ.1)THEN
        DO K=1,KS
        DO L=2,LA
        AB(L,K)=SQRT(AB(L,K)*ABLPF(L,K))
        ENDDO
        ENDDO
!      ELSE
!       DO K=1,KS
!       DO L=2,LA
!       AB(L,K)=SQRT(AB(L,K)*ABEFF(L,K))
!       ENDDO
!       ENDDO
!     ENDIF
!
      DO K=1,KC
      DO L=1,LC
      UHDY(L,K)=UHLPF(L,K)*DYU(L)
      VHDX(L,K)=VHLPF(L,K)*DXV(L)
      SAL(L,K)=SALLPF(L,K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=1,LC
      UHDYE(L)=UHDYE(L)+UHDY(L,K)*DZC(K)
      VHDXE(L)=VHDXE(L)+VHDX(L,K)*DZC(K)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      ENDIF
!
!**********************************************************************C
!     
! **  ADJUST THE STEADY OR INITIAL MEAN MASS TRANSPORT FIELD  
!
      CALL ADJMMT
!
!**********************************************************************C
!
! **  INITIALIZE COUNTERS
!
      NCTS=0
!
!**********************************************************************C
!**********************************************************************C
!
! **  BEGIN THE TIME INTEGRATION LOOP
!
!----------------------------------------------------------------------C
!
      DO 1000 N=1,NTS
      NCTS=NCTS+1
!
!**********************************************************************C
!
! **  ADVANCE CONCENTRATION FIELDS
!
!----------------------------------------------------------------------C
!
      IF(ISCDCA(1).NE.1)THEN
        CALL CALCONC (2,0)
       ELSE
        CALL CALCONC (3,0)
      ENDIF
      IF(ISHOW.GE.1) CALL SHOWVAL1
!
!----------------------------------------------------------------------C
!
! **  CHECK RANGE OF SALINITY AND DYE CONCENTRATION
!
      IF(ISMMC.EQ.1)THEN
!
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(SAL(L,K).GT.SALMAX)THEN
       SALMAX=SAL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(SAL(L,K).LT.SALMIN)THEN
       SALMIN=SAL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
!
      WRITE(6,6001)N
      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN
!
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(DYE(L,K).GT.SALMAX)THEN
       SALMAX=DYE(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(DYE(L,K).LT.SALMIN)THEN
       SALMIN=DYE(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
!
      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN
!
      ENDIF
!
 6001 FORMAT('  N=',I10)
 6002 FORMAT('  SALMAX=',F14.4,5X,'I,J,K=',(3I10))
 6003 FORMAT('  SALMIN=',F14.4,5X,'I,J,K=',(3I10))
 6004 FORMAT('  DYEMAX=',F14.4,5X,'I,J,K=',(3I10))
 6005 FORMAT('  DYEMIN=',F14.4,5X,'I,J,K=',(3I10))
!
!**********************************************************************C
!
! **  ADVANCE MMT FIELDS
!
!----------------------------------------------------------------------C
!
      IF(ISSSMMT.EQ.0)THEN
!
      IF(NCTS.EQ.NTSMMT)THEN
      NCTS=0
!
      DO L=2,LA
      H1P(L)=HLPF(L)
      UHDY1E(L)=UHDYE(L)
      VHDX1E(L)=VHDXE(L)
      ENDDO
!
      DO K=1,KC
      DO L=2,LA
      SAL1(L,K)=SALLPF(L,K)
      UHDY1(L,K)=UHDY(L,K)
      VHDX1(L,K)=VHDX(L,K)
      ENDDO
      ENDDO
!
      CALL RESTRAN
!
      RNTCMMT=FLOAT(NTSMMT)/FLOAT(NTSPTC)
      DO L=2,LA
      HP(L)=H1P(L)+(HLPF(L)-H1P(L))/(RNTCMMT*FLOAT(NTSPTC))
      UHDYE(L)=0.
      VHDXE(L)=0.
      ENDDO
!
!     IF(ISLTMT.EQ.1)THEN
        DO K=1,KS
        DO L=2,LA
        AB(L,K)=SQRT(AB(L,K)*ABLPF(L,K))
        ENDDO
        ENDDO
!      ELSE
!       DO K=1,KS
!       DO L=2,LA
!       AB(L,K)=SQRT(AB(L,K)*ABEFF(L,K))
!       ENDDO
!       ENDDO
!     ENDIF
!
      DO K=1,KC
      DO L=1,LC
      UHDY(L,K)=UHLPF(L,K)*DYU(L)
      VHDX(L,K)=VHLPF(L,K)*DXV(L)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO L=1,LC
      UHDYE(L)=UHDYE(L)+UHDY(L,K)*DZC(K)
      VHDXE(L)=VHDXE(L)+VHDX(L,K)*DZC(K)
      ENDDO
      ENDDO
!
      CALL ADJMMT
!
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/(RNTCMMT*FLOAT(NTSPTC))
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      ELSE
!
      DO L=2,LA
      HP(L)=H1P(L)+(HLPF(L)-H1P(L))*FLOAT(NCTS+1)/FLOAT(NTSMMT)
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))*FLOAT(NCTS+1)/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))*FLOAT(NCTS+1)/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))*FLOAT(NCTS+1)/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))*FLOAT(NCTS+1)/FLOAT(NTSMMT)
      ENDDO
      ENDDO
!
      ENDIF
!
      ENDIF
!
!**********************************************************************C
!
 1000 CONTINUE
!
!**********************************************************************C
!**********************************************************************C
!
! **  WRITE GRAPHICS FILES FOR RESIDUAL VARIABLES
!
!----------------------------------------------------------------------C
!
      DO K=1,KC
      DO L=1,LC
      SALLPF(L,K)=SAL(L,K)
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SALINITY CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
!
      IF(ISRSPH(1).EQ.1.AND.ISTRAN(1).GE.1)THEN
       CALL RSALPLTH(1,SALLPF)
      ENDIF
!
      IF(ISRSPH(2).EQ.1.AND.ISTRAN(2).GE.1)THEN
       CALL RSALPLTH(2,TEMLPF)
      ENDIF
!
      IF(ISRSPH(3).EQ.1.AND.ISTRAN(3).GE.1)THEN
       CALL RSALPLTH(3,DYELPF)
      ENDIF
!
      IF(ISRSPH(4).EQ.1.AND.ISTRAN(4).GE.1)THEN
       CALL RSALPLTH(4,SEDTLPF)
      ENDIF
!
      IF(ISRSPH(4).EQ.1.AND.ISTRAN(5).GE.1)THEN
       CALL RSALPLTH(5,SNDTLPF)
      ENDIF
!
!FIX      IF(ISRSPH(5).EQ.1.AND.ISTRAN(6).GE.1)THEN
!FIX       CALL RSALPLTH(6,TOXLPF)
!FIX      ENDIF
!
      IF(ISRSPH(7).EQ.1.AND.ISTRAN(7).GE.1)THEN
       CALL RSALPLTH(7,SFLLPF)
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES:
! **  SUBROUTINE RVELPLTH
!
      IF(ISRVPH.EQ.1)THEN
       CALL RVELPLTH
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SURFACE ELEVATION PLOTTING IN HORIZONTAL PLANES:
! **  SUBROUTINE RVELPLTH
!
      IF(ISRPPH.EQ.1)THEN
       CALL RSURFPLT
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL SALINITY AND VERTICAL MASS DIFFUSIVITY CONTOURING IN 
! **  3 VERTICAL PLANES:  SUBROUTINE RSALPLTV
!
      IF(ISRSPV(1).GE.1)THEN
       CALL RSALPLTV(1)
      ENDIF
!
!----------------------------------------------------------------------C
!
! **  RESIDUAL NORMAL AND TANGENTIAL VELOCITY CONTOURING AND AND 
! **  TANGENTIAL VELOCITY VECTOR PLOTTING IN VERTICAL PLANES:
! **  SUBROUTINE RVELPLTV
!
      IF(ISRVPV.GE.1)THEN
       CALL RVELPLTV
      ENDIF
!
!**********************************************************************C
!
      CLOSE(99) 
!
!**********************************************************************C
!
      RETURN
      END