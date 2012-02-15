!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! binsolu.h
!! This file calculates the amount of water associated with the aqueous organic
!! species. The data for each compound comes from a UNIFAC calculation of the 
!! water activity in a binary prganic-water solution. Values are linearly 
!! interpolated.
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. FUNCTION HydrophilicOrganicWaterContent(InParticle, Temp, RH)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!This function calculates the water content (in mol/particle) caused by the
!hydrophilic oprganics. Based on previous work by Betty Pun and Rob Griffin:
! see comments
REAL*8 FUNCTION HydrophilicOrganicWaterContent(InParticle, Temp, RH)
      
!/**************************************************************************
! Include file: binsolu.h  
!
!Purpose: concentration of solute (micromol solute / microgram water)
!    in a binary solution at RH = Aw = 0, 0.1, 0.2, 0.3 etc
!
!
!Revision history: Developed by Betty Pun, December 1999, under CARB funding
!                  Reconfigured for Fortran and only 5 Type A species by
!              Rob Griffin January, 2000
!                  Reconfigured for combined philic-phobic module, April, 2001
!              Matt Alvarado April, 2006   
!               Reordered organic species (SOA1 - SO10)
!			   Matt Alvarado June, 2006   
!               Added organic species (POA1-POA8, LEVO, CBIO, CPD1-CPD3)
!				Revised data using corrected UNIFAC values
!***********************************************************************/
!        graduation of RH scale in molalbin definition, 
!      starting with RH = 0 (first index of first dim) 
!                     ending with 1 (last index of first dim)
! 
		
	USE Chemistry, ONLY : AqMolecularMass
	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle
	REAL*8				   :: Temp, RH

	!! Internal Variables
	INTEGER        :: I, LowerIndex
	REAL*8		   :: MOLALBIN(0:10,1:23), RHGrad
	real*8		   :: MOLALRH(1:23), ORGCONC(1:23), LowerMolal, UpperMolal
	REAL*8         :: a, b, SUM 
	  
	  RHGRAD = 0.1

!       binary solution molality (umol/ug water) corresponding to RH 
!       at 10% graduations --> 11 rows
!       NOTE: This must have same order as HydrophilicOrgChems!
!
      !POA1, C29 alkane, MW = 408
	  MOLALBIN(0,1) = 555.6
      MOLALBIN(1,1) = 241.3
      MOLALBIN(2,1) = 120.0
      MOLALBIN(3,1) = 79.61
      MOLALBIN(4,1) = 59.44
      MOLALBIN(5,1) = 47.34
      MOLALBIN(6,1) = 39.25
      MOLALBIN(7,1) = 33.48
      MOLALBIN(8,1) = 29.16
      MOLALBIN(9,1) = 25.79
      MOLALBIN(10,1)= 0.0
	        
	 !POA2, C4 diacid, MW = 118.0
	  MOLALBIN(0,2) = 555.6
      MOLALBIN(1,2) = 0.5467
      MOLALBIN(2,2) = 0.2518
      MOLALBIN(3,2) = 0.1529
      MOLALBIN(4,2) = 0.1028
      MOLALBIN(5,2) = 7.213E-2
      MOLALBIN(6,2) = 5.104E-2
      MOLALBIN(7,2) = 3.517E-2
      MOLALBIN(8,2) = 2.212E-2
      MOLALBIN(9,2) = 9.968E-3
      MOLALBIN(10,2)= 0.0

	 !POA3, naphthalene diacid, MW = 216.0
	  MOLALBIN(0,3) = 555.6
      MOLALBIN(1,3) = 1.348
      MOLALBIN(2,3) = 0.6469
      MOLALBIN(3,3) = 0.4126
	  MOLALBIN(4,3) = 0.2950
      MOLALBIN(5,3) = 0.2241
      MOLALBIN(6,3) = 0.1764
      MOLALBIN(7,3) = 0.1420
      MOLALBIN(8,3) = 0.1057
      MOLALBIN(9,3) = 9.483E-2
      MOLALBIN(10,3)= 0.0

	 !POA4, 5 ring PAH, MW = 276.0
	  MOLALBIN(0,4) = 555.6
      MOLALBIN(1,4) = 230.8
      MOLALBIN(2,4) = 115.2
      MOLALBIN(3,4) = 76.57
	  MOLALBIN(4,4) = 57.22
      MOLALBIN(5,4) = 45.63
      MOLALBIN(6,4) = 37.89
      MOLALBIN(7,4) = 32.36
      MOLALBIN(8,4) = 28.22
      MOLALBIN(9,4) = 25.00
      MOLALBIN(10,4)= 0.0

	 !POA5, 5 ring branched alkane, MW = 412.0
	  MOLALBIN(0,5) = 555.6
      MOLALBIN(1,5) = 321.2
      MOLALBIN(2,5) = 159.4
      MOLALBIN(3,5) = 106.0
	  MOLALBIN(4,5) = 79.07
      MOLALBIN(5,5) = 63.05
      MOLALBIN(6,5) = 52.32
      MOLALBIN(7,5) = 44.65
      MOLALBIN(8,5) = 38.91
      MOLALBIN(9,5) = 34.43
      MOLALBIN(10,5)= 0.0

	 !POA6, carboxybenzoic acid, MW = 166.0
	  MOLALBIN(0,6) = 555.6
      MOLALBIN(1,6) = 0.9824
      MOLALBIN(2,6) = 0.4646
      MOLALBIN(3,6) = 0.2914
	  MOLALBIN(4,6) = 0.2042
      MOLALBIN(5,6) = 0.1514
      MOLALBIN(6,6) = 0.1156
      MOLALBIN(7,6) = 8.940E-2
      MOLALBIN(8,6) = 6.894E-2
      MOLALBIN(9,6) = 5.189E-2
      MOLALBIN(10,6)= 0.0

	 !POA7, C18 acid, MW = 284.0
	  MOLALBIN(0,7) = 555.6
      MOLALBIN(1,7) = 4.429
      MOLALBIN(2,7) = 2.159
      MOLALBIN(3,7) = 1.401
	  MOLALBIN(4,7) = 1.022
      MOLALBIN(5,7) = 0.7948
      MOLALBIN(6,7) = 0.6428
      MOLALBIN(7,7) = 0.5339
      MOLALBIN(8,7) = 0.4521
      MOLALBIN(9,7) = 0.3882
      MOLALBIN(10,7)= 0.0

	 !POA8, 2 ring branched alkane, MW = 390.0
	  MOLALBIN(0,8) = 555.6
      MOLALBIN(1,8) = 273.0
      MOLALBIN(2,8) = 135.6
      MOLALBIN(3,8) = 90.09
	  MOLALBIN(4,8) = 67.20
      MOLALBIN(5,8) = 53.53
      MOLALBIN(6,8) = 40.41
      MOLALBIN(7,8) = 37.87
      MOLALBIN(8,8) = 32.99
      MOLALBIN(9,8) = 29.19
      MOLALBIN(10,8)= 0.0

	  !SOA1, C2 diacid, MW = 90
	  MOLALBIN(0,9) = 555.6
      MOLALBIN(1,9) = 0.49396
      MOLALBIN(2,9) = 0.22574
      MOLALBIN(3,9) = 0.13569
      MOLALBIN(4,9) = 0.09018
      MOLALBIN(5,9) = 0.06232
      MOLALBIN(6,9) = 0.04319
      MOLALBIN(7,9) = 0.0289
      MOLALBIN(8,9) = 0.0174
      MOLALBIN(9,9) = 0.00755
      MOLALBIN(10,9)= 0.0
      
	  !SOA2, C8 diacid, MW = 184.
	  MOLALBIN(0,10) = 555.6
      MOLALBIN(1,10) = 0.45742
      MOLALBIN(2,10) = 0.21651
      MOLALBIN(3,10) = 0.13549
      MOLALBIN(4,10) = 0.09439
      MOLALBIN(5,10) = 0.06918
      MOLALBIN(6,10) = 0.05184
      MOLALBIN(7,10) = 0.03891
      MOLALBIN(8,10) = 0.02848
      MOLALBIN(9,10) = 0.01921
      MOLALBIN(10,10)= 0.0
      
	  !SOA3, dialdehyde, MW = 154
	  MOLALBIN(0,11) = 555.6
      MOLALBIN(1,11) = 0.86272
      MOLALBIN(2,11) = 0.40057
      MOLALBIN(3,11) = 0.24605
      MOLALBIN(4,11) = 0.16855
      MOLALBIN(5,11) = 0.1216
      MOLALBIN(6,11) = 0.08988
      MOLALBIN(7,11) = 0.06671
      MOLALBIN(8,11) = 0.0486
      MOLALBIN(9,11) = 0.03329
      MOLALBIN(10,11)= 0.0

	  !SOA4, monoacid, MW= 186
      MOLALBIN(0,12) = 555.6
      MOLALBIN(1,12) = 0.55832
      MOLALBIN(2,12) = 0.263
      MOLALBIN(3,12) = 0.16412
      MOLALBIN(4,12) = 0.11418
      MOLALBIN(5,12) = 0.08375
      MOLALBIN(6,12) = 0.06308
      MOLALBIN(7,12) = 0.0478
      MOLALBIN(8,12) = 0.03574
      MOLALBIN(9,12) = 0.02543
      MOLALBIN(10,12)= 0.0

	  !SOA5, monoaldehyde, MW = 186
      MOLALBIN(0,13) = 555.6
      MOLALBIN(1,13) = 0.92947
      MOLALBIN(2,13) = 0.42461
      MOLALBIN(3,13) = 0.25726
      MOLALBIN(4,13) = 0.17392
      MOLALBIN(5,13) = 0.12435
      MOLALBIN(6,13) = 0.09149
      MOLALBIN(7,13) = 0.06809
      MOLALBIN(8,13) = 0.05035
      MOLALBIN(9,13) = 0.03601
      MOLALBIN(10,13)= 0.0
	  
	  !SOA6, like nitrobenzoic acid, MW = 211 NOTE: Not Griffin values: fixed!
	  MOLALBIN(0,14) = 555.6
      MOLALBIN(1,14) = 3.454
      MOLALBIN(2,14) = 1.671
      MOLALBIN(3,14) = 1.077
      MOLALBIN(4,14) = 0.7786
      MOLALBIN(5,14) = 0.5994
      MOLALBIN(6,14) = 0.4794
      MOLALBIN(7,14) = 0.3933
      MOLALBIN(8,14) = 0.3283
      MOLALBIN(9,14) = 0.2773
      MOLALBIN(10,14) = 0.0

	  !SOA7, like benzoic axid, MW = 178
      MOLALBIN(0,15) = 555.6
      MOLALBIN(1,15) = 1.98908
      MOLALBIN(2,15) = 0.95252
      MOLALBIN(3,15) = 0.60644
      MOLALBIN(4,15) = 0.43292
      MOLALBIN(5,15) = 0.32835
      MOLALBIN(6,15) = 0.25820
      MOLALBIN(7,15) = 0.20764
      MOLALBIN(8,15) = 0.16923
      MOLALBIN(9,15) = 0.13830
      MOLALBIN(10,15) = 0.0

	  !SOA8, PAH, MW = 217
      MOLALBIN(0,16) = 555.6
      MOLALBIN(1,16) = 63.80219
      MOLALBIN(2,16) = 31.65389
      MOLALBIN(3,16) = 20.93034
      MOLALBIN(4,16) = 15.57215
      MOLALBIN(5,16) = 12.35612
      MOLALBIN(6,16) = 10.21267
      MOLALBIN(7,16) = 8.68083
      MOLALBIN(8,16) = 7.53173
      MOLALBIN(9,16) = 6.63831
      MOLALBIN(10,16) = 0.0

	  !SOA9, C16, MW = 303
      MOLALBIN(0,17) = 555.6
      MOLALBIN(1,17) = 2.81667
      MOLALBIN(2,17) = 1.37478
      MOLALBIN(3,17) = 0.89336
      MOLALBIN(4,17) = 0.65205
      MOLALBIN(5,17) = 0.50677
      MOLALBIN(6,17) = 0.40951
      MOLALBIN(7,17) = 0.33965
      MOLALBIN(8,17) = 0.28692
      MOLALBIN(9,17) = 0.24558
      MOLALBIN(10,17) = 0.0
      
	  !SO10, biological, MW=217 NOTE: Not Griffin values: fixed!
	  MOLALBIN(0,18) = 555.6
      MOLALBIN(1,18) = 2.084
      MOLALBIN(2,18) = 1.017
      MOLALBIN(3,18) = 0.6592
      MOLALBIN(4,18) = 0.4793
      MOLALBIN(5,18) = 0.3705
      MOLALBIN(6,18) = 0.2971
      MOLALBIN(7,18) = 0.2441
      MOLALBIN(8,18) = 0.2036
      MOLALBIN(9,18) = 0.1715
      MOLALBIN(10,18) = 0.0

	  !LEVO, levoglucosan, MW=162.1
	  MOLALBIN(0,19) = 555.6
      MOLALBIN(1,19) = 0.2010
      MOLALBIN(2,19) = 0.1055
      MOLALBIN(3,19) = 7.079E-2
      MOLALBIN(4,19) = 5.189E-2
      MOLALBIN(5,19) = 3.952E-2
      MOLALBIN(6,19) = 3.046E-2
      MOLALBIN(7,19) = 2.324E-2
      MOLALBIN(8,19) = 1.689E-2
      MOLALBIN(9,19) = 1.082E-2
      MOLALBIN(10,19) = 0.0

	  !CBIO, cellobiosan, MW=324.3
	  MOLALBIN(0,20) = 555.6
      MOLALBIN(1,20) = 0.1071
      MOLALBIN(2,20) = 5.622E-2
      MOLALBIN(3,20) = 3.782E-2
      MOLALBIN(4,20) = 2.785E-2
      MOLALBIN(5,20) = 2.136E-2
      MOLALBIN(6,20) = 1.663E-2
      MOLALBIN(7,20) = 1.289E-2
      MOLALBIN(8,20) = 9.701E-3
      MOLALBIN(9,20) = 6.667E-3
      MOLALBIN(10,20) = 0.0

	  !CPD1, biomass burning neutral, MW=330.3
	  MOLALBIN(0,21) = 555.6
      MOLALBIN(1,21) = 32.02
      MOLALBIN(2,21) = 15.82
      MOLALBIN(3,21) = 10.42
      MOLALBIN(4,21) = 7.718
      MOLALBIN(5,21) = 6.097
      MOLALBIN(6,21) = 5.015
      MOLALBIN(7,21) = 4.242
      MOLALBIN(8,21) = 3.663
      MOLALBIN(9,21) = 3.211
      MOLALBIN(10,21) = 0.0
      
	  !CPD2, biomass burning monoacid, MW=332.4
	  MOLALBIN(0,22) = 555.6
      MOLALBIN(1,22) = 0.7109
      MOLALBIN(2,22) = 0.3391
      MOLALBIN(3,22) = 0.2149
      MOLALBIN(4,22) = 0.1525
      MOLALBIN(5,22) = 0.1148
      MOLALBIN(6,22) = 8.946E-2
      MOLALBIN(7,22) = 7.108E-2
      MOLALBIN(8,22) = 5.702E-2
      MOLALBIN(9,22) = 4.574E-2
      MOLALBIN(10,22) = 0.0

	  !CPD3, biomass burning humic acid, MW=730.7
	  MOLALBIN(0,23) = 555.6
      MOLALBIN(1,23) = 0.2802
      MOLALBIN(2,23) = 0.1336
      MOLALBIN(3,23) = 8.455E-2
      MOLALBIN(4,23) = 5.989E-2
      MOLALBIN(5,23) = 4.497E-2
      MOLALBIN(6,23) = 3.489E-2
      MOLALBIN(7,23) = 2.756E-2
      MOLALBIN(8,23) = 2.191E-2
      MOLALBIN(9,23) = 1.732E-2
      MOLALBIN(10,23) = 0.0
      

	LowerIndex = INT(10.0*RH) 
	b = MOD(10*RH, 1.0)
	!WRITE(*,*) "Lower Index", LowerIndex,  "b: ", b
	a = 1-b
	
	Sum = 0.0
	DO I = 1, 23
		LowerMolal = MOLALBIN(LowerIndex,I)
		UpperMolal = MOLALBIN(LowerIndex+1,I)
		MOLALRH(I) = (a*LowerMolal + b*UpperMolal) !umol/ug = mol org/g H2O
		!IF(I .EQ. 9) WRITE(*,*) "Molality: ", MOLALRH(I)*1000, "RH: ",  RH
	
		ORGCONC(I) = InParticle%AqOrgChems(I)*InParticle%NumberofParticles !Mol org/cm3
		
		!ZSR formula
		Sum = Sum + (ORGCONC(I)/MOLALRH(I)) !g H2O /cm3 air
	
	END DO

	
	Sum = Sum/AqMolecularMass(1) !Water content in mol/cm3 air
	!WRITE(*,*) "Water mass: ", AqMolecularMass(1)
	HydrophilicOrganicWaterContent = Sum/InParticle%NumberofParticles !water content in mol/particle
	!WRITE(*,*) "OrgWater: ", HydrophilicOrganicWaterContent*InParticle%NumberofParticles, "RH: ", RH

END FUNCTION HydrophilicOrganicWaterContent



