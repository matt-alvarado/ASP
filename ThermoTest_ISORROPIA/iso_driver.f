C     Driver for the ASP comparison: reads case totals (mol/m3), T, metastable flag;
C     loops RH 0.90 -> 0.20 by 0.02 and prints gas, liquid and solid results.
      PROGRAM DRV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WI(8), WT(8), GAS(3), AERLIQ(15), AERSLD(19),
     &          CNTRL(2), OTHER(9)
      CHARACTER*15 SCASI
      READ (*,*) (WI(I),I=1,8), TEMPI, METS
      CNTRL(1) = 0.D0
      CNTRL(2) = DBLE(METS)
      RHI = 0.92D0
      DO 10 J = 1, 36
         RHI = RHI - 0.02D0
         CALL ISOROPIA (WI, RHI, TEMPI, CNTRL, WT, GAS, AERLIQ,
     &                  AERSLD, SCASI, OTHER)
         WRITE (*,'(F6.3,3ES14.6,15ES14.6,19ES14.6,2X,A15)') RHI,
     &      (GAS(I),I=1,3), (AERLIQ(I),I=1,15), (AERSLD(I),I=1,19),
     &      SCASI
 10   CONTINUE
      END
