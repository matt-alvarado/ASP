!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! StepASP.f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!! 06/27  2007   Matt Alvarado     Added InputConcentrations
!!								   Added OutputConcentrations
!!                                 Added InitializeASP
!! 07/16  2007   Matt Alvarado     Set any aerosol concentration below
!!									1.0e-30 mol/particle to 0.
!! 07/18  2007   Matt Alvarado     Corrected output concentrations for
!!									solid electrolytes
!! 08/31  2007   Matt Alvarado     Added aerosol size to interface
!! 10/09  2007   Matt Alvarado     Changed zero criteria for aerosol concentrations in 
!!								   InputConcentrations to 1.0e-40 mol/particle
!! 10/12  2007   Matt Alvarado     Removed Scale from photorates: now in chemtrop.F
!! 02/19  2009   Matt Alvarado     Changed DFPORT to IFPORT (to use ifort compiler)
!! 08/26  2010   Matt Alvarado     Comment out IFPORT (to use pgf90 compiler)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES						        !!
!! 1. ModelParameters						!!
!! 2. GridPointFields						!!
!! 3. Aerosols								!!
!! 4. Condensation							!!
!! 5. Time									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!	
!! 1. SUBROUTINE StepASPOnce ()										!!
!! 2. SUBROUTINE InputConcentrations
!! 3. SUBROUTINE OutputConcentrations
!! 4. SUBROUTINE InitializeASP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE StepASP

		IMPLICIT NONE

		PRIVATE

		PUBLIC ::	ASPInterface, InitializeASP

	CONTAINS


	SUBROUTINE StepASPOnce (TimeStep)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine recalculates the gas-phase reactions rates,  !!
	!! then steps gas chemistry, condensation, coagulation, and     !!
	!! updates the optical properties of the aerosol.				!!
	!! All calculations are done in gridbox 1,1,1 MJA 052307	    !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


		!USE IFPORT,				 ONLY : RTC, ETIME
		USE	ModelParameters
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: TimeStep
		
		!! Internal Variables
		INTEGER :: I, NumBins, XP, YP, ZP
		REAL*8 :: X(3), TEMP
	
		TYPE(Particle), POINTER :: Current
		
		
		!Set grid point
		XP = 1
		YP = 1
		ZP = 1

		!Calculate Temperature
		X    = LocateGridPoint(XP,YP,ZP,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
					
		!Step Gas Chemistry
		WRITE(*,*) "Before StepGasChem"
		CALL StepGasChemistry(TimeStep,1,XP,YP,ZP)
		WRITE(*,*) "StepGasChem Okay"
			
		!Step Condensation, if desired
		current => particles(XP,YP,ZP)%first
		IF(.TRUE.) THEN
			!CALL FindAqElectrolyteEquilibriumForGridPoint (XP, YP, ZP, UpdateThermo = .TRUE.) 
			CALL StepCondensationAll (TimeStep, 1, XP, YP, ZP)
			!Note above calls RegridAerosol automatically
			!with recalculates radius and density
			!WRITE(*,*) "Cond Okay"
		END IF

		!Step Coagulation, if desired
		IF (.TRUE.) THEN	
			CALL StepSectionalCoagulationJacobson(XP,YP,ZP,TimeStep)
			!Note that above already calls RegridAerosol
			!which recalculates radius and density
			!WRITE(*,*) "Coag Okay"
		END IF

		!Step Optical Properties Particles
		CALL SortAerosolAtGridPointForCoagulation (XP, YP, ZP) 
		current => particles(XP,YP,ZP)%first
		I = 1
		DO WHILE(associated(current))
    			
			!Force particle temperature to be env. Temperature
			!WARNING: This is a kludge, since the program is calculating
			!a huge temperature for the largest particles, and I'm not sure why
			current%Temperature = TEMP
					
			!Recalculate optical parameters
			WRITE(*,*) "Call Optical, Particle #", I
			CALL ShellRefIndAndRad(Current)
			WRITE(*,*) "Optical Okay"
				
			I = I + 1
		
			current => current%next
				
		END DO !END STEP OPTICAL PROPERTIES
	
		RETURN
	END SUBROUTINE StepASPOnce

	SUBROUTINE InputConcentrations(Temp, Press, Dens, GasConc, &
									AeroNumConc, AeroMassConc, PhotoRates)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine takes in Temperature, pressure,				!!
	!! gas and aerosol concentrations, and photolysis rates from	!!
	!! a 3D model and updates the appropriate concentrations in ASP !!
	!! grid box 1,1,1. MJA 052307									!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
		!USE IFPORT,				ONLY : RTC, ETIME
		USE	ModelParameters,	ONLY : HowManyBins
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8, INTENT(IN) :: Temp  !(in K)
		REAL*8, INTENT(IN) :: Press !(in mbar)
		REAL*8, INTENT(IN) :: Dens  !(in kg/m3 dry air)
		REAL*8 :: GasConc(HowManyEvolveGasChems)
		REAL*8, INTENT(IN) :: AeroNumConc(HowManyBins)
		REAL*8 :: AeroMassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
		REAL*8, INTENT(IN) :: PhotoRates(12)

		!! Internal Variables
		INTEGER :: I, J, L, RxnIndex, NumAq, GG
		REAL*8 :: HeteroNumberConc(HowManyBins), HeteroRadii(HowManyBins)
		REAL*8 :: Scale

		TYPE(Particle), POINTER :: Current
		
		!Set temperature
		CALL SetTempField (1, 1, 1, Temp)

		!Set Pressure
		CALL SetPressField (1, 1, 1, Press)
		
		!Update chemical concentrations
		DO I = 1, HowManyEvolveGasChems
			IF(GasConc(I) .LT. 0.) GasConc(I) = 0.
		END DO 
		CALL UpdateChemicalConcentrations(1, 1, 1, GasConc(1:HowManyEvolveGasChems))
		
		!WRITE(*,*) "Before Aero Chems Input", GetRelativeHumidityFromGrid (1, 1, 1)
		!Update aerosol concentrations
		CALL SortAerosolAtGridPointForCoagulation (1, 1, 1) 
		current => particles(1,1,1)%first
		I = 1
		
		NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
		DO WHILE(associated(current))
    		current%numberofparticles = AeroNumConc(I)
			IF(current%Numberofparticles .GT. 0.) THEN
				IF(current%dry) THEN
					current%dry = .FALSE.
				END IF
				DO J = 1, NumAq
					IF(AeroMassConc(I,J) .LT. 0.) AeroMassConc(I,J) = 0.
					Current%AqChems(J) = AeroMassConc(I,J)*Dens*1.0e-12 &
					                      /(current%Numberofparticles*AqMolecularMass(J))
					IF(Current%AqChems(J) .LT. 1.0e-40) Current%AqChems(J) = 0.
					
					!WRITE(*,*) I, J, AeroMassConc(I,J), Current%AqChems(J), Dens, current%Numberofparticles, AqMolecularMass(J)
					!PAUSE
				END DO
				!WRITE(*,*) current%AqChems(8)+current%AqChems(12)
				
				
				DO J = 1, HowManyOrgChems
						IF(AeroMassConc(I,J+NumAq) .LT. 0.) AeroMassConc(I,J+NumAq) = 0.
						Current%OrgChems(J) = AeroMassConc(I,J+NumAq)*Dens*1.0e-12 &
					                      /(current%Numberofparticles*OrgMolecularMass(J))
						IF(Current%OrgChems(J) .LT. 1.0e-40) Current%OrgChems(J) = 0.
				
				
				END DO
			
				DO J = 1, HowManyAqOrgChems
						IF(AeroMassConc(I,J+NumAq+HowManyOrgChems) .LT. 0.) AeroMassConc(I,J+NumAq+HowManyOrgChems) = 0.
						Current%AqOrgChems(J) = AeroMassConc(I,J+NumAq+HowManyOrgChems)*Dens*1.0e-12 &
					                      /(current%Numberofparticles*AqOrgMolecularMass(J))					
						IF(Current%AqOrgChems(J) .LT. 1.0e-40) Current%AqOrgChems(J) = 0.
				
				END DO
			ELSE !No particles
				current%numberofparticles = 0.
				DO J = 1, NumAq
					Current%AqChems(J) = 0.
				END DO
				
				DO J = 1, HowManyOrgChems
					Current%OrgChems(J) = 0.
				END DO
				
				DO J = 1, HowManyAqOrgChems
					Current%AqOrgChems(J) = 0.
				END DO				
			END IF
			!Force particle temperature to be env. Temperature
			current%Temperature = Temp
				
			!Equilibrate, recalculate radius and density
			!WRITE(*,*) "Before Recauculate Radius", current%effectiveradius
			CALL FindElectrolyteEquilibrium (current, UpdateThermo=.TRUE., FirstEquilibration=.TRUE., ReturnType=GG)
			CALL RecalculateRadius(current)		
			!WRITE(*,*) "After Recauculate Radius", current%effectiveradius
			I = I + 1
			current => current%next
				
		END DO
		!WRITE(*,*) "After Aero Chems Input"
		
		!Recalculate all reaction rates
		CALL RecalculateReactionRates(IFGAS,IFAQ)
		
		!WRITE(*,*) "Before Photorates"
		!Reset Photolysis Rates
		DO L = 1, GasPhaseReactionRateArraySize(1)
						
			!! Only do photolysis reactions
			IF (INT(GasPhaseReactionRates(L,1)) .EQ. 1 .AND. INT(GasPhaseReactionRates(L,3)) .EQ. 0.) THEN	
								
				!Find reaction index number
				RxnIndex = INT(GasPhaseReactionRates(L,4))

				SELECT CASE (RxnIndex)
					CASE(4) !NO2 => NO + O
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(1)
						WRITE(*,*) "ASP Photo", GasPhaseChemicalRates(1,1,1,L)
					CASE(5) !NO3 => NO + O2
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(2)
					CASE(6) !NO3 => NO2 + O
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(3)
					
					CASE(3) !O3 => O + O2
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(4)
					CASE(2) !O3 => O1D + O2
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(5)
					CASE(12) !HONO => Products
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(6)

					CASE(11) !H2O2 => 2OH
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(7)
					CASE(15) !HCHO => 2HO2 + CO
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(8)
					CASE(16) !HCHO => H2 + CO
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(9)

					CASE(17) !ALD2 => Products
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(10)
					CASE(23) !KETL => Products
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(11)
					CASE(22) !MGLY => Products
						GasPhaseChemicalRates(1,1,1,L) = PhotoRates(12)

				END SELECT
			END IF
		END DO
		!WRITE(*,*) "After Photorates"
		!Get particle number conc. and radii for hetero. rate calculation
		current => particles(1,1,1)%first
		I = 1
		DO WHILE(associated(current))
			HeteroNumberConc(I) = current%NumberofParticles
			HeteroRadii(I) = current%EffectiveRadius/100.0 !convert from cm to m
			!ResetHeteroRates expects radii in units of m
			current => current%next
			I = I+1
		END DO
				
		!Calculate Heterogeneous rates from size dist. info
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HowManyBins, .FALSE.)
		!WRITE(*,*) "Hetero Okay"

	END SUBROUTINE InputConcentrations

	SUBROUTINE OutputConcentrations(Temp, Press, Dens, GasConc, &
									AeroNumConc, AeroMassConc, &
									ExtCoeff, SingScat, Assym, Radius)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine outputs Temperature, pressure,				!!
	!! gas and aerosol concentrations, and optical properties from	!!
	!! a 3D model and updates the appropriate concentrations in ASP !!
	!! grid box 1,1,1. MJA 052307									!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
		!USE IFPORT,				ONLY : RTC, ETIME
		USE	ModelParameters,	ONLY : HowManyBins
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: Temp, Press, Dens, X(3)
		REAL*8, DIMENSION(18) :: ExtCoeff, SingScat, Assym
		REAL*8 :: GasConc(HowManyEvolveGasChems)
		REAL*8 :: AeroNumConc(HowManyBins)
		REAL*8 :: AeroMassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
		REAL*8 :: Radius(HowManyBins)
		
		!! Internal Variables
		INTEGER :: I, J, NumAq
		REAL*8 :: SaltConc(HowManyAqChems)
	
		TYPE(Particle), POINTER :: Current
		
		
		!Output Temperature and Pressure
		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		Temp = GetTemp(X(1),X(2),X(3))
		Press = GetPress(X(1),X(2),X(3))
	
		!Output Gas chemical concentrations
		DO I = 1, HowManyEvolveGasChems
			GasConc(I) = GridGasChem(1,1,1,I)
		END DO

		!Update aerosol concentrations
		CALL SortAerosolAtGridPointForCoagulation (1, 1, 1) 
		current => particles(1,1,1)%first
		I = 1
		NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
		DO WHILE(associated(current))
    		AeroNumConc(I) = current%numberofparticles
			Radius(I) = current%effectiveradius
			
			!Water
			AeroMassConc(I,1) = current%AqChems(1)*(current%Numberofparticles*AqMolecularMass(1)) &
					                    /(Dens*1.0e-12)
												
			!Do this step so that we only pass neutral electrolytes
			!to dynamics code - helps prevent odd charge balance problems
			!WRITE(*,*) "Before output HypotheticalElectrolyteConcentrations"
			SaltConc = HypotheticalElectrolyteConcentrations (current)
			!WRITE(*,*) SaltConc(5)+SaltConc(8)+SaltConc(12), current%AqChems(5)+current%AqChems(8)+current%Aqchems(12)
		
			IF(current%Numberofparticles .GT. 0.) THEN
				DO J = 2, HowManyAqChems
					AeroMassConc(I,J) = (SaltConc(J)+Current%AqChems(J))*(current%Numberofparticles*AqMolecularMass(J)) &
					                    /(Dens*1.0e-12)
				END DO
				
				DO J = HowManyAqChems+1, NumAq
					AeroMassConc(I,J) = 0.
				END DO

				DO J = 1, HowManyOrgChems
						AeroMassConc(I,J+NumAq) = Current%OrgChems(J)*(current%Numberofparticles*OrgMolecularMass(J)) &
					                    /(Dens*1.0e-12)
				END DO
				
				DO J = 1, HowManyAqOrgChems
						AeroMassConc(I,J+NumAq+HowManyOrgChems) = Current%AqOrgChems(J)*(current%Numberofparticles*AqOrgMolecularMass(J)) &
					                    /(Dens*1.0e-12)
				END DO
			ELSE !No particles
				AeroNumConc(I) = 0.
				DO J = 1, HowManyAqChems
					AeroMassConc(I,J) = 0.
				END DO
				
				DO J = 1, HowManyOrgChems
					AeroMassConc(I,J+HowManyAqChems) = 0.
				END DO
				
				DO J = 1, HowManyAqOrgChems
					AeroMassConc(I,J+HowManyAqChems+HowManyOrgChems) = 0.
				END DO				
			END IF
		
			I = I + 1
			current => current%next
				
		END DO

		!Output Average Aerosol optical properties
		CALL AerosolOptProp (1, 1, 1, ExtCoeff, SingScat, Assym)
		!WRITE(*,*) ExtCoeff(1), SingScat(1), Assym(1)
		!PAUSE

	END SUBROUTINE OutputConcentrations

	SUBROUTINE InitializeASP()

	!! Include modules
	!USE IFPORT,				 ONLY : RTC, ETIME
	USE	ModelParameters
	USE InfrastructuralCode	
	USE GridPointFields
	USE Chemistry
	USE Aerosols		
	USE Condensation
	USE OutputRoutines
	USE Coagulation
	
	implicit none

	!! Define the local variables
	integer :: i, j, k, l, seedarray(2),q, r, NumBins
        INTEGER :: time(8)
	REAL*4  :: etimearray(2)
	TYPE(Particle), POINTER :: Current, Env
	CHARACTER (len = 8)     :: ErrorDate
	CHARACTER (len = 10)    :: ErrorTime

	!! Initialize the random number generator 
        call DATE_AND_TIME(values=time)     ! Get the current time (for pgi compiler)
        seedarray(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
        !seedarray(1) = rtc() !for ifort compiler
        seedarray(2) = FLOOR(SQRT(REAL(seedarray(1))))
	CALL Random_Seed(put=seedarray)

		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize (DBLE(1.0*cm), DBLE(1.0*cm), DBLE(1.0*cm))
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()
	
	!! Send a Welcome Message
	CALL Transcript ("**********************************************************")
	CALL Transcript ("** Welcome to the ASP Model developed by Matt Alvarado  **")
	CALL Transcript ("** (mjalvara@mit.edu).                                  **")
	CALL Transcript ("** ASP is an updated, expanded version of MELAM         **")
	CALL Transcript ("** by H. D. Steele (donnan@mit.edu).                    **")
	CALL Transcript ("**********************************************************")
	CALL Date_and_Time(ErrorDate, ErrorTime)
	CALL Transcript ("*********************************************")
	CALL Transcript ("** Run Started at                          **")
	CALL Transcript ("** "//ErrorTime(1:2)//":"//ErrorTime(3:4)//":"//ErrorTime(5:6)//" on "//ErrorDate(5:6)//		&
					  "/"//ErrorDate(7:8)//"/"//ErrorDate(1:4)//" "//"                 "//"**")
	CALL Transcript ("*********************************************")
				
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams

		!! INITIALIZE the DISSOLUTION routine
		CALL SetAllDissolution

		!! INITIALIZE the PARTICLES	
		CALL InitializeParticlesSectional
				
		!!WARNING Force initial water equilibrium between gas and aerosol
		!!(Uses open system for water eq.)
		!CALL EquilibrateInternallyAtGridPoint(1,1,1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
		CALL RegridAerosol (1, 1, 1)

	END SUBROUTINE InitializeASP

	SUBROUTINE ASPInterface(Timestep, Temp, Press, Dens, GasConc, &
									NumConc, MassConc, PhotoRates, &
									ExtCoeff, SingScat, Assym, Radius)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine connects ASP to an external 3D Eulerian model 
	!! MJA 062707									!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
		!USE IFPORT,				ONLY : RTC, ETIME
		USE	ModelParameters,	ONLY : HowManyBins
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: Timestep !(in s)
		REAL*8 :: Temp  !(in K)
		REAL*8 :: Press !(in mbar)
		REAL*8 :: Dens  !(in kg/m3 dry air)
		REAL*8 :: GasConc(HowManyEvolveGasChems)
		REAL*8 :: NumConc(HowManyBins)
		REAL*8 :: MassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
		REAL*8 :: Radius(HowManyBins)
		REAL*8 :: PhotoRates(12)
		REAL*8, DIMENSION(18) :: ExtCoeff, SingScat, Assym
		INTEGER :: I
	
		!WRITE(*,*) "Before Conc Input"
		CALL InputConcentrations(Temp, Press, Dens, GasConc, &
									NumConc, MassConc, PhotoRates)
		!WRITE(*,*) "After Conc Input"
		CALL StepASPOnce (TimeStep)
		!WRITE(*,*) "Before Conc Output"
		CALL OutputConcentrations(Temp, Press, Dens, GasConc, &
									NumConc, MassConc, &
									ExtCoeff, SingScat, Assym, Radius)
		!WRITE(*,*) "After Conc Output"

	END SUBROUTINE ASPInterface
END MODULE StepASP
