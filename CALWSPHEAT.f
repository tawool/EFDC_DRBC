C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION CALWSPHEAT(TEM,SAL)
C
C **  CALWSPHEAT CALCULATES THE WATER SPECIFIC HEAT IN SI UNITS (J/KgC)
C **  FOR SEAWATER AT SEA LEVEL AS A FUNCTIONS OF SALT AND TEMP
C **  USING THE EQUATION FROM      
C **  D.T. Jamieson, J.S. Tudhope, R. Morris and G. Cartwright,
C **  Physical properties of sea water solutions: heat capacity,
C **  Desalination, 7(1) (1969) 23–30. HNR_GHD 7/2023
C **  as given by:
C **  Sharqawy, Mostafa H., Lienhard, John H., and Zubair, Syed M. (2010). 
C **  Thermophysical properties of seawater: a review of existing correlations and data. 
C **  Desalination and Water Treatment. 16 (2010) 354–380
C
C **  CREATED BY HUGO RODRIGUEZ ON JULY 2023
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
      INCLUDE 'EFDC.PAR'
C
C**********************************************************************C
C
      SS=SAL
      TT=TEM
      IF(TT.LT.0.0) TT=0.0
      IF(TT.GT.90.0) TT=90.0
      IF(SS.LT.0.0) SS=0.0
      IF(SS.GT.150.0) SS=150.0
       
       
      TT=1.00024*(TT+273.15)
C
      AS = 5.328 - 9.76E-2*SS + 4.04E-4*SS*SS
      BS = -6.913E-3 + 7.351E-4*SS - 3.15E-6*SS*SS
      CS = 9.6E-6 - 1.927E-6*SS + 8.23E-9*SS*SS
      DS = 2.5E-9 + 1.666E-9*SS - 7.125E-12*SS*SS
       
      CALWSPHEAT = 1000.0*(AS + BS*TT + CS*TT*TT + DS*TT*TT*TT)
C
      RETURN
      END
