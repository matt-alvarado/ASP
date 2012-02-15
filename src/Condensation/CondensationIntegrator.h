!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! CondensationIntegrator.h
!! This is the numerical integration routine for the flux-limited
!! kinetic dissolution of inorganic compounds into the aerosol aqueous phase.
!! As the main routine calls LSODES, the required ODE and Jacobian 
!! functions are given here as well
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY					      			!!
!!									      	!!
!! Month  Year   Name              Description			        	!!
!! 07     2006   Matt Alvarado     Began Update History				!!
!! 08/28  2006   Matt Alvarado     Commented out some write statements		!!
!! 08/30  2006   Matt Alvarado     1. Edited StepCondensationAll to call	!!
!!						RegridAerosol at end.		!!
!!				   2. Edited StepCondensation to exit		!!
!!					       water eq. loop after 5 iterations!!
!!				  3. Removed commented code from		!!
!!					       StepCondensation			!!
!!				  4. Fixed Jacobian for NO3 and Cl		!!
!! 09/08  2006   Matt Alvarado     Changed water eq. loop in StepCondensation	!!
!! 09/13  2006   Matt Alvarado     Changed Ammonia eq loop in StepCondensation  !!
!! 09/18  2006   Matt Alvarado     Added eq check to StepCondensationAll		!!
!! 10/03  2006   Matt Alvarado     Updated GasSatConc and NH3SatConc			!!
!! 10/13  2006   Matt Alvarado     Created BigCondensationOdeEvaluator and		!!
!!				     	BigCondensationJacobianEvaluator		!!
!! 10/17  2006   Matt Alvarado     Moved EmptyGridCell check and aerosol sorting to 
!!				      	StepCondensationAll()
!! 07/12  2007   Matt Alvarado     Have inorganic condensation other than
!!								   NH3 and H2SO4 shut down when RH is less than 30%
!!									To remove oscillations, exit water eq loop
!!									in StepCondensation after 50 iterations if 
!!									change is less than 1%.
!! 07/13  2007   Matt Alvarado     Remove coupling to solid salt in AcidGasSatConc 
!!									and NH3SatConc when metastable
!!									or when above DRH for salt.
!! 07/16  2007   Matt Alvarado     Exit water equilibrium in StepCondensation 
!!									after 1000 iterations to prevent oscilating
!!									nonconvergence at low RH.
!!									Try forcing NH3 condensation to be flag 10
!!									in NH3SatConc
!!									Try forcing HCl, HNO3 condensation to be flag 10
!!									in AcidGasSatConc
!!									Add flag to catch when Cl- and NO3- are near 
!!									zero and evaporating to BigCondensationODEEvaluator
!!									and BigCondensationJacobianEvaluator
!!									Set HClCoeff(I) to 0.0 for all particles
!!									when Cl is not present (skips this condensation!)
!! 07/20  2007   Matt Alvarado     Added flag to skip calculation of AcidGasSatConc
!!									in BigCondODE and BigConJac if RH is less than 30%
!! 07/30  2007   Matt Alvarado     Put in all flags for condensation, except skip NH3 flag 22
!!									(NH42SO4 formation) in NH3SatConc
!!									Increased atol to reduce charge balance errors
!! 07/31  2007   Matt Alvarado     Put code in BigCondensationODEEvaluator to catch
!!									NH4+ below 0 and recover by dissociating NH4NO3
!!									or NH4Cl 
!! 08/01  2007   Matt Alvarado     Added code to just exit gas-solid eq at end of 
!!									SUBROUTINE BigStepInorgCondensation after max iterations
!!									Added NH42SO4 and NH4HSO4 to NH4 below zero flag
!! 08/06  2007   Matt Alvarado     Forced NH3SatConc to use gas-solid equilibrium when
!!									coupled to NO3- or Cl-.
!! 08/07  2007   Matt Alvarado     Added line to NH3SatConc to force Flag 10 if
!!									RH is below 30% (and thus we are doing gas-solid eq
!!									for Cl- and NO3-)
!!									Fixed gas-solid eq loop at end of BigStepInorgCondensation
!! 08/08  2007   Matt Alvarado      Removed forcing of gas-solid in NH3SatConc, because this caused 
!!									the 2D simulations to crash
!! 08/13  2007   Matt Alvarado      Fixed sign error in BigCondensationJacobianEvaluator (Case 7, HClFlag=10, d/dCl-(dH+/dt) 
!!									Fixed bug in Gas-Solid loop in BigStepInorgCondensation
!! 08/17  2007   Matt Alvarado      Try Jacobson gas-solid idea in GasSatConcs (See 17.14, p. 594)
!!									Removed all 20 flags from NH3SatConc and AcidGasSatConc
!!									Created GasSatConcAll
!!									Fixed below 0 checks
!! 08/20  2007   Matt Alvarado      Fixed GasSatConcAll
!!									Added Cl- below 0 check to BigStepInorgCondensation
!! 08/21  2007   Matt Alvarado      Redid below 0 checks in BigStepInorgCondensation
!!									Set flag to limit HClSat and HNO3Sat to total HCl
!!									and HNO3 available.
!! 08/22  2007   Matt Alvarado      Added a second H+ below 0 check to BigCondensationODEEvaluator
!!									to handle the case where both H+ and NH4+ are below 0
!! 09/01  2007   Matt Alvarado      Forced below 0 H+ recovery to use bigger of NH4+ and H2O eq estimates
!! 09/03  2007   Matt Alvarado      Removed NH4+ from Cl- and NO3- below zero loops in LSODES
!!									Hopefully will fix charge balance errors.
!! 09/04  2007   Matt Alvarado      Made flag for low RH set to 20%
!!									Added in new Jacobian for gas-solid rxns
!!									Changes to GasSatConcAll, BigCondensationJacobianEvaluator, 
!!									and Ia and Ja in BigStepInorgCondensation
!!                                  Added Na, Mg and Ca Cl and NO3 salts to Clsum and NO3sum
!!									in BigODE and BigJac.
!!									Added Na, Mg and Ca salts to below 0 check at end of
!!									BigStepInorgCondensation
!! 09/05  2007   Matt Alvarado      Removed forcing of gas-solid eq when NH4Cl or NH4NO3 are present
!! 09/18  2007   Matt Alvarado      Forced NO3NearZero to be true when SkipNO3 is true,
!!									as was already done for Cl
!! 09/21  2007   Matt Alvarado      Cleaned up Jacobian and ODE, removed AcidGasSatConc and
!!									NH3SatGasConc, and DissociateAmmoniaSalt
!!									Added Tiny to Jacobain and GasSatConcAll calcs
!! 09/27  2007   Matt Alvarado      Removed 1.0e-4 molec/particle and evap check for Cl- and NO3-
!! 10/01  2007   Matt Alvarado      Added EquilibrateInorgGridPoint option to StepCondensation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE StepCondensationAll ()								!!
!! 2. SUBROUTINE BigStepInorgCondensation ()						!!
!! 3. SUBROUTINE BigCondensationODEEvaluator()						!!
!! 4. SUBROUTINE BigCondensationJacobianEvaluator()					!!
!! 5. SUBROUTINE DissolutionFactors									!!
!! 6. SUBROUTINE CondensationODEEvaluator()							!!
!! 7. SUBROUTINE CondensationJacobianEvaluator()					!!
!! 8. SUBROUTINE StepCondensation ()
!! 9. SUBROUTINE GasSatConcAll										!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!! This allows the main program to call all the condensation integrators together
!! Also regrid aerosol automatically
SUBROUTINE StepCondensationAll (TimestepSize, NumTimeSteps, XP, YP, ZP)

	USE InfrastructuralCode
	USE ModelParameters
	USE Aerosols, ONLY :  Particles, Particle,			&
						  ParticleArray,				&
						  RegridAerosol,				&
						  SortAerosolAtGridPointForCondensation
						  							
	IMPLICIT NONE

	!! External Variables
	REAL*8  :: TimeStepSize		!! Should be in Seconds
	INTEGER :: NumTimeSteps, XP, YP, ZP

	!! Internal Variables
	INTEGER :: I
	REAL*8 :: S1MassFracs(2), S2MassFracs(6), S3MassFracs(2), &
				S4MassFracs(3), S5MassFracs(2), S6MassFracs(3), &
				S7MassFracs(8), S8MassFracs(3), S9MassFracs(4), &
				S10MassFracs(4), InitialBurden(10)
	REAL*4  :: etimearray(2)
	REAL*8  :: CondTime
	LOGICAL :: EmptyGridCell
	TYPE(Particle), POINTER :: FirstParticle, CurrentParticle
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! If there are no aerosol at the gridcell (or if they are all hydrophobic),  !!
	!! then there is no need to condense                                          !!
	EmptyGridCell = .TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CurrentParticle => Particles(XP,YP,ZP)%First

5	IF (ASSOCIATED(CurrentParticle)) THEN
		IF ((CurrentParticle%Sectional .AND. CurrentParticle%NumberOfParticles .EQ. 0) .OR. CurrentParticle%Dry) THEN
			IF(ASSOCIATED(CurrentParticle%Next)) THEN
				CurrentParticle => CurrentParticle%Next			
				GOTO 5
			END IF
		ELSE
			EmptyGridCell = .FALSE.
		END IF
	END IF
	
	
	IF (.NOT.EmptyGridCell) THEN
	
		!If doing equilibrium (09/18/06 MJA)
		IF(.NOT. DissolutionEquilibriumOrFlux) THEN
			CALL StepCondensation (TimestepSize, 1, XP, YP, ZP)
	
		!If doing flux
		ELSE
			IF (ThermoBulkMode) &
			CALL ERROR("Can't try to use the kinetic based condensation and dissolution routine when running in bulk mode "// &
				   "(specified in AerosolModes.in).  Please remedy.")

			!! Sort the gridcell, pushing all empty sections to the end
			CALL SortAerosolAtGridPointForCondensation (XP, YP, ZP)
		
			!Calculate surrogate organic concentrations
			IF(DoHydrophobicOrgDissolutionFlag .OR. DoHydrophilicOrgDissolutionFlag) &
				CALL CalculateSurrogateGasConc(XP, YP, ZP, S1MassFracs, S2MassFracs, &
					S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
					S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)
	
			DO I = 1, NumTimeSteps
				!WRITE(*,*) "StepDissolutionAll Check 0"
				IF(DoDissolutionFlag .OR. DoCondensationFlag) &
					CALL StepCondensation (TimestepSize, 1, XP, YP, ZP)
		
				!WRITE(*,*) "StepDissolutionAll Check 1"
				IF(DoHydrophobicOrgDissolutionFlag) &
					CALL StepOrgCondensation (TimestepSize, 1, XP, YP, ZP)
		
				!WRITE(*,*) "StepDissolutionAll Check 2"
				IF(DoHydrophilicOrgDissolutionFlag) &
					CALL StepAqOrgCondensation (TimestepSize, 1, XP, YP, ZP)
				!WRITE(*,*) "StepDissolutionAll Check 3"
			END DO

			!Update individual gas-phase concentrations
			IF(DoHydrophobicOrgDissolutionFlag .OR. DoHydrophilicOrgDissolutionFlag) &
				CALL UpdateGasPhaseOrganics(XP, YP, ZP, S1MassFracs, S2MassFracs, &
					S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
					S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)
	
		END IF

		!WRITE(*,*) "StepDissolutionAll Check 4"
		!Move particles across bin boundaries, if necessary
		CALL RegridAerosol (XP, YP, ZP)		
	
	END IF !Empty Grid Cell
			
	!WRITE(*,*) "StepDissolutionAll Check 5"
END SUBROUTINE StepCondensationAll

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Call LSODES or the equilibrium routine to integrate forward	 !!
!! the concentration of aerosol / gas-phase water at the grid-	 !!
!! cell level.  It must construct a joint array of all of the	 !!
!! aerosol plus the grid cell for each of the condensing species !!
!! XP, YP, and ZP are the gridpoint numbers in gas grid space.	 !!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BigStepInorgCondensation (TimestepSize, NumTimeSteps, XP, YP, ZP)

	USE InfrastructuralCode, ONLY : INT2STR, ERROR,					&
									Warn, Transcript,				&
									REAL2STR

	USE GridPointFields,	 ONLY : LocateGridPoint,				&
									getM,							&
									GasPhaseChemGridIndex,			&
									GetGridCellChemBurden,			&
									ReplaceGridCellChemBurden,		&
									GetTemp,						&
									GetRelativeHumidity,			&
									AddToGridCellChemBurden,		&
									GridGasChem,					&
									GetPress,						&
									GetSatVapPress,					&
									GetRelativeHumidityfromGrid,    &
									GetAirDensity

	USE Aerosols,            ONLY : Particles, Particle,			&
									ParticleArray,					&
									RecalculateRadius,				&
									EqCoeff,						&
									Molality,						&
									FindElectrolyteEquilibrium,		&
									CalculateIonicStrength,			&
									KusikMeissner,					&
									SortAerosolAtGridPointForCondensation,			&
									EquilibriumWaterContent, &
									FindAqElectrolyteEquilibriumForGridPoint, &
									EquilibrateAqueousDissociationReaction, &
									DUMPPARTICLECONTENTSATERROR

	USE Chemistry,			 ONLY : HowManyGasChems,				&
									HowManyEvolveGasChems,			&
									HowManyEvolveAqChems,			&
									GasPhaseChemicalNames,			&
									HowManyCondensingChems,			&
									AqEquilibriaList,				&
									HowManyAqChems,					&
									HowManyAqCations,				&
									HowManyAqAnions, &
									FindChem,						&
									AqPhaseChemicalNames, &
									HowManyAqEqReactions

	USE ModelParameters,     ONLY : Avogadro,								&
									AerosolWaterEquilibriumRHThreshold,		&
									WaterEquilibriumIndex,					&
									ThermoBulkMode,							&
									DissolutionEquilibriumOrFlux,			&
									AqThermoEquilibriumError,				&
									WaterContentPrecision,					&
									micron,									&
									AmmoniaFlux,							&
									DoHydrophobicOrgDissolutionFlag,        &
									DoHydrophilicOrgDissolutionFlag,		&
									ProtonIndex, HydroxyIndex, &
									Metastable, AQTHERMONUMBALLOWEDITERATIONS


	IMPLICIT NONE

	!! External Variables
	REAL*8  :: TimeStepSize		!! Should be in Seconds
	INTEGER :: NumTimeSteps, XP, YP, ZP

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!! LSODES Related Internal Variables !!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8,  ALLOCATABLE	:: y(:), rwork(:), Dum1(:), Dum2(:), Dum3(:), Dum4(:)
	REAL*8					:: t,outtime,rtol,atol
	INTEGER, ALLOCATABLE	:: iwork(:)
	INTEGER					:: itol,itask,istate,iopt,mf,neq,LIW,LRW
		
	!! Non-LSODES Local Variables
	INTEGER				:: I, J, K, L, LLLL, C, Q, allocation_error, Rxn, HowManyLinks, YLen, GG,					&
						   HowManyNonZeroJacobianTerms, GasChemIndex1, TimeStep,	&
						   ReturnType, WhichDominantIon, WhichSecondIon , Flag, Dum, RxnType, &
						   InternalReturnType, NH4Index, HIndex, OHIndex, GasSolidRxn, NH3i, NH4i, &
						   SO4Index, NO3Index, ClIndex, NH3Index, JaIndex, CaseFlag
	REAL*8				:: II, JJ, KK, MM, ChemTransferArray(HowManyGasChems), Temperature,	&
						   RH, XX(3), OrgWaterContent, Store, Store3,Store4, GasPhaseBurden, &
						   Electrolyte, Change, Sum
	REAL*8, ALLOCATABLE :: InitialConcsArray(:), ChemStorageArray(:)
	LOGICAL :: EmptyGridCell
	REAL*8  :: BeginTime,EndTime
	REAL*4  :: etimearray(2)
	REAL*8  :: WaterTime, CondTime, AmmoniaTime, GasPhaseBurden2

	LOGICAL, PARAMETER :: Scaffolding = .FALSE., OpenSystem = .FALSE., ForceEquilibrium = .TRUE.
	LOGICAL :: Equilibrated

	TYPE(Particle), POINTER :: FirstParticle, CurrentParticle

	!! SCSCSCSC
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript(">>>Entering BigStepInorgCondensation<<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Find GRID POINT Values for NON-AEROSOL !!
	!! DEPENDENT Parameters                   !!
	XX          = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
	Temperature = GetTemp(XX(1),XX(2),XX(3))	
	RH		    = GetRelativeHumidity(XX(1),XX(2),XX(3))
	WRITE(*,*) "Temperature: ", Temperature
	WRITE(*,*) "Press: ", GetPress(XX(1),XX(2),XX(3))
	!WRITE(*,*) "Sat Vap Press: ", GetSatVapPress (XX(1), XX(2), XX(3))
	WRITE(*,*) "RH: ", RH
	WRITE(*,*) "Other RH: ", GetRelativeHumidityFromGrid(1,1,1)
	WRITE(*,*) "Air Density: ", GetAirDensity(XX(1),XX(2),XX(3))
	!WRITE(*,*) "Water (molec/cm3): ", GridGasChem(XP,YP,ZP,1)


	
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! IF EQUILIBRIUM RATHER THAN FLUX, DO THAT HERE AND RETURN !!
	IF (.NOT. DissolutionEquilibriumOrFlux .AND. RH .LE. AerosolWaterEquilibriumRHThreshold) THEN
		WRITE(*,*) "Calling Equilibrium Routine"
		CALL EquilibrateGridPointAll(XP, YP, ZP, .TRUE., ForceWaterEquilibrium=.FALSE., WaterOpenSystem=.FALSE.)
		RETURN
	END IF

	IF (Scaffolding) CALL TRANSCRIPT("Using LSODES to integrate dissolution.")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Everything should now be sorted           !!
	!! Back up to the head of the aerosol list   !!
	FirstParticle   => Particles(XP,YP,ZP)%First
	CurrentParticle => FirstParticle

	!! COUNT the number of LINKS in the gridcell that we should condense upon
	HowManyLinks  =  1
15	IF(ASSOCIATED(CurrentParticle%Next)) THEN
		CurrentParticle => CurrentParticle%Next
		IF (CurrentParticle%NumberOfParticles .GT. 0 .AND. .NOT.CurrentParticle%Dry) THEN
			HowManyLinks  = HowManyLinks + 1
			GOTO 15
		END IF
	ENDIF
	CurrentParticle => FirstParticle


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set the lengths of the various !!
	!! arrays and allocate them !!!!!!!!
	YLen	= 12 + HowManyLinks*16	  !! As per the definition below
	NEQ		= 4  + 6*HowManyLinks	  !! Total number of equations

	!! Temporarily set this to the maximum value (which is also for reaction type 4)
	HowManyNonZeroJacobianTerms = 36*HowManyLinks+10 !44*HowManyLinks + 10
	LRW  = 60 +  25*NEQ + 3*HowManyNonZeroJacobianTerms
	LIW  = 31 + neq + HowManyNonZeroJacobianTerms

	!! Allocate the LSODES working arrays
	ALLOCATE (Y(ylen), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of Y could not proceed in StepCondensation()")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Y-VECTOR Structure										!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!    INTEGRATED TERMS (1 to 4 + 6*HowManyLinks)										!!
	!!    (1)								            : H2SO4 Gas Phase Concentration		!!
	!!    (2)								            : HNO3 Gas Phase Concentration		!!
	!!    (3)								            : HCl Gas Phase Concentration		!!
	!!    (4)								            : NH3 Gas Phase Concentration		!!	
	!!    (5 to 4 + HowManyLinks)				        : SO4-- Aerosol Concentrations		!!
	!!    (5 + HowManyLinks to 4 + 2*HowManyLinks)      : NO3- Aerosol Concentrations		!!
	!!    (5 + 2*HowManyLinks to 4 + 3*HowManyLinks)    : Cl- Aerosol Concentrations		!!
	!!    (5 + 3*HowManyLinks to 4 + 4*HowManyLinks)    : NH4+ Aerosol Concentrations		!!
	!!    (5 + 4*HowManyLinks to 4 + 5*HowManyLinks)    : NH3(aq) Aerosol Concentrations	!!
	!!    (5 + 5*HowManyLinks to 4 + 6*HowManyLinks)    : H+ Aerosol Concentrations			!!
	!!																						!!
	!!	  PARAMETER TERMS																	!!
	!!	  (5 + 6*HowManyLinks to ??)														!!
	!!																						!!
	!!    (5 + 6*HowManyLinks to 7 + 6*HowManyLinks)    : Location (XP, YP, ZP)				!!
	!!    (8 + 7*HowManyLinks to 7 + 8*HowManyLinks)    : (S/H) for HNO3(g) <=> H+ & NO3-	!!
	!!    (8 + 8*HowManyLinks to 7 + 9*HowManyLinks)    : (S/H) for HCl(g) <=> H+ & Cl-		!!
	!!    (8 + 9*HowManyLinks to 7 + 10*HowManyLinks)   : (S/H) for NH3(g) & H+(aq) <=> NH4+!!
	!!    (8 + 10*HowManyLinks to 7 + 11*HowManyLinks)  : (S/H) for NH3(g) <=> NH3 (aq)		!!
	!!	  (8 + 11*HowManyLinks to 7 + 12*HowManyLinks)  : NumberofParticles per bin			!!
	!!	  (8 + 12*HowManyLinks to 7 + 13*HowManyLinks)  : H2SO4 Rate Coefficient			!!
	!!	  (8 + 13*HowManyLinks to 7 + 14*HowManyLinks)  : HNO3 Rate Coefficient				!!
	!!	  (8 + 14*HowManyLinks to 7 + 15*HowManyLinks)  : HCl Rate Coefficient				!!
	!!	  (8 + 15*HowManyLinks to 7 + 16*HowManyLinks)  : NH3 Rate Coefficient				!!
	!!    (8 + 16*HowManyLinks)							: H2SO4 Rxn Index					!!
	!!    (9 + 16*HowManyLinks)							: HNO3 Rxn Index					!!
	!!    (10 + 16*HowManyLinks)						: HCl Rxn Index						!!
	!!    (11 + 16*HowManyLinks)						: NH3(g) & H+ <=> NH4+ Rxn Index	!!
	!!    (12 + 16*HowManyLinks)						: NH3(g) <=> NH3(aq) Rxn Index		!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! HowManyLinks is the number of (non-zero) sections in the size dist					!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	ALLOCATE (rwork(LRW), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of RWORK could not proceed in StepCondensation()")

	ALLOCATE (iwork(LIW), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of IWORK could not proceed in StepCondensation()")

	ALLOCATE (InitialConcsArray(NEQ), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of InitialConcsArray could not proceed in StepCondensation()")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! ChemStorgaeArray is indexed the same as Y (see above)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ALLOCATE (ChemStorageArray(NEQ), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of ChemStorageArray could not proceed in StepCondensation()")

	!Dummy arrays (included so we don't have to adjust DissolutionFactors())
	ALLOCATE (Dum1(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of Dum1 could not proceed in StepCondensation()")
	ALLOCATE (Dum2(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of Dum2 could not proceed in StepCondensation()")
	ALLOCATE (Dum3(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of Dum3 could not proceed in StepCondensation()")
	ALLOCATE (Dum4(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of Dum4 could not proceed in StepCondensation()")


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Initialize the LSODES Work Vectors !!
	DO I = 1, LRW !!!!!!!!!!!!!!!!!!!!!!!!!!
		RWORK(I) = 0.
	END DO
	DO I = 1, LIW
		IWORK(I) = 0
	END DO

	!! Initialize the Y vector
	DO J = 1, YLen
		Y(J) = 0.
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Take the # of Particles counts from each of the aerosol particles !!
	CurrentParticle => Particles(XP,YP,ZP)%First !!!!!!!!!!!!!!!!!!!!!!!!!!
	DO J = 1, HowManyLinks
		Y(11*HowManyLinks+7+J) = CurrentParticle%NumberOfParticles
		CurrentParticle		  => CurrentParticle%Next
	END DO

	IF (Scaffolding) CALL TRANSCRIPT("Cond. Check 1")
		
	!!!!!!!!!!!!!!!!!!!!!!!!!
	!! LOOP OVER TIMESTEPS !!
	DO TimeStep = 1, NumTimeSteps

		!! Integrate to this time
		OutTime =  TimeStep * TimeStepSize

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Initialize Chem Storage Array !!
		DO I = 1, NEQ
			ChemStorageArray(I) = 0.
		END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Set LSODES-Related Variables (Notes Excerpted !!
		!! from LSODES file)							 !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! itask  = an index specifying the task to be performed.				!!
		!!          input only.  itask has the following values and meanings.	!!
		!!          1  means normal computation of output values of y(t) at		!!
		!!             t = tout (by overshooting and interpolating).			!!
		!!          2  means take one step only and return.						!!
		!!          3  means stop at the first internal mesh point at or		!!
		!!             beyond t = tout and return.								!!
		!!          4  means normal computation of output values of y(t) at		!!
		!!             t = tout but without overshooting t = tcrit.				!!
		!!             tcrit must be input as rwork(1).  tcrit may be equal to	!!
		!!             or beyond tout, but not behind it in the direction of	!!
		!!             integration.  this option is useful if the problem		!!
		!!             has a singularity at or beyond t = tcrit.				!!
		!!          5  means take one step, without passing tcrit, and return.	!!
		!!             tcrit must be input as rwork(1).							!!
		!!																		!!
		!!          note..  if itask = 4 or 5 and the solver reaches tcrit		!!
		!!          (within roundoff), it will return t = tcrit (exactly) to	!!
		!!          indicate this (unless itask = 4 and tout comes before tcrit,!!
		!!          in which case answers at t = tout are returned first).		!!
		itask  = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! istate = an index used for input and output to specify the
		!!         the state of the calculation.
		!!
		!!         on input, the values of istate are as follows.
		!!          1  means this is the first call for the problem
		!!             (initializations will be done).  see note below.
		!!          2  means this is not the first call, and the calculation
		!!             is to continue normally, with no change in any input
		!!             parameters except possibly tout and itask.
		!!             (if itol, rtol, and/or atol are changed between calls
		!!             with istate = 2, the new values will be used but not
		!!             tested for legality.)
		!!          3  means this is not the first call, and the
		!!             calculation is to continue normally, but with
		!!             a change in input parameters other than
		!!             tout and itask.  changes are allowed in
		!!             neq, itol, rtol, atol, iopt, lrw, liw, mf,
		!!             the conditional inputs ia and ja,
		!!             and any of the optional inputs except h0.
		!!             in particular, if miter = 1 or 2, a call with istate = 3
		!!             will cause the sparsity structure of the problem to be
		!!             recomputed (or reread from ia and ja if moss = 0).
		!!          note..  a preliminary call with tout = t is not counted
		!!          as a first call here, as no initialization or checking of
		!!          input is done.  (such a call is sometimes useful for the
		!!          purpose of outputting the initial conditions.)
		!!          thus the first call for which tout .ne. t requires
		!!          istate = 1 on input.
		!!
		!!          on output, istate has the following values and meanings.
		!!           1  means nothing was done, as tout was equal to t with
		!!              istate = 1 on input.  (however, an internal counter was
		!!              set to detect and prevent repeated calls of this type.)
		!!           2  means the integration was performed successfully.
		!!          -1  means an excessive amount of work (more than mxstep
		!!              steps) was done on this call, before completing the
		!!              requested task, but the integration was otherwise
		!!              successful as far as t.  (mxstep is an optional input
		!!              and is normally 500.)  to continue, the user may
		!!              simply reset istate to a value .gt. 1 and call again
		!!              (the excess work step counter will be reset to 0).
		!!              in addition, the user may increase mxstep to avoid
		!!              this error return (see below on optional inputs).
		!!          -2  means too much accuracy was requested for the precision
		!!              of the machine being used.  this was detected before
		!!              completing the requested task, but the integration
		!!              was successful as far as t.  to continue, the tolerance
		!!              parameters must be reset, and istate must be set
		!!              to 3.  the optional output tolsf may be used for this
		!!              purpose.  (note.. if this condition is detected before
		!!              taking any steps, then an illegal input return
		!!              (istate = -3) occurs instead.)
		!!          -3  means illegal input was detected, before taking any
		!!              integration steps.  see written message for details.
		!!              note..  if the solver detects an infinite loop of calls
		!!              to the solver with illegal input, it will cause
		!!              the run to stop.
		!!          -4  means there were repeated error test failures on
		!!              one attempted step, before completing the requested
		!!              task, but the integration was successful as far as t.
		!!              the problem may have a singularity, or the input
		!!              may be inappropriate.
		!!          -5  means there were repeated convergence test failures on
		!!              one attempted step, before completing the requested
		!!              task, but the integration was successful as far as t.
		!!              this may be caused by an inaccurate jacobian matrix,
		!!              if one is being used.
		!!          -6  means ewt(i) became zero for some i during the
		!!              integration.  pure relative error control (atol(i)=0.0)
		!!              was requested on a variable which has now vanished.
		!!              the integration was successful as far as t.
		!!          -7  means a fatal error return flag came from the sparse
		!!              solver cdrv by way of prjs or slss (numerical
		!!              factorization or backsolve).  this should never happen.
		!!              the integration was successful as far as t.
		!!
		!!          note.. an error return with istate = -1, -4, or -5 and with
		!!          miter = 1 or 2 may mean that the sparsity structure of the
		!!          problem has changed significantly since it was last
		!!          determined (or input).  in that case, one can attempt to
		!!          complete the integration by setting istate = 3 on the next
		!!          call, so that a new structure determination is done.
		!!
		!!          note..  since the normal output value of istate is 2,
		!!          it does not need to be reset for normal continuation.
		!!          also, since a negative input value of istate will be
		!!          regarded as illegal, a negative output value requires the
		!!          user to change it, and possibly other inputs, before
		!!          calling the solver again.
		istate = 1	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! iopt   = an integer flag to specify whether or not any optional		!!
		!!          inputs are being used on this call.  input only.			!!
		!!          the optional inputs are listed separately below.			!!
		!!          iopt = 0 means no optional inputs are being used.			!!
		!!                   default values will be used in all cases.			!!
		!!          iopt = 1 means one or more optional inputs are being used.	!!
		iopt   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! mf     = the method flag.  used only for input.					!!
		!!          the standard choices for mf are..						!!
		!!            mf = 10  for a nonstiff problem,						!!
		!!            mf = 21 or 22 for a stiff problem with ia/ja supplied	!!
		!!                     (21 if jac is supplied, 22 if not),			!!
		!!            mf = 121 for a stiff problem with jac supplied,		!!
		!!                     but not ia/ja,								!!
		!!            mf = 222 for a stiff problem with neither ia/ja nor	!!
		!!                     jac supplied.								!!
		!!          the sparseness structure can be changed during the		!!
		!!          problem by making a call to lsodes with istate = 3.		!!
		mf = 21  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!h0      rwork(5)  the step size to be attempted on the first step.
        !!            the default value is determined by the solver.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		rwork(5) = 1.e-10

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! hmax    rwork(6)  the maximum absolute step size allowed.
		!!                   the default value is infinite.
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		rwork(6) = 1.0

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! hmin     = the minimum allowed timestep !!
!		rwork(7) = 1.e-10 !!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! atol   = an absolute error tolerance parameter, either a scalar or	!!
		!!          an array of length neq.  input only.						!!
		!!																		!!
		!!             the input parameters itol, rtol, and atol determine		!!
		!!          the error control performed by the solver.  the solver will	!!
		!!          control the vector e = (e(i)) of estimated local errors		!!
		!!          in y, according to an inequality of the form				!!
		!!                      rms-norm of ( e(i)/ewt(i) )   .le.   1,			!!
		!!          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),			!!
		!!          and the rms-norm (root-mean-square norm) here is			!!
		!!          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))	!!
		!!          is a vector of weights which must always be positive, and	!!
		!!          the values of rtol and atol should all be non-negative.		!!
		!!          the following table gives the types (scalar/array) of		!!
		!!          rtol and atol, and the corresponding form of ewt(i).		!!
		!!																		!!
		!!             itol    rtol       atol          ewt(i)					!!
		!!              1     scalar     scalar     rtol*abs(y(i)) + atol		!!
		!!              2     scalar     array      rtol*abs(y(i)) + atol(i)	!!
		!!              3     array      scalar     rtol(i)*abs(y(i)) + atol	!!
		!!              4     array      array      rtol(i)*abs(y(i)) + atol(i)	!!
		!!																		!!
		!!          when either of these parameters is a scalar, it need not	!!
		!!          be dimensioned in the user-s calling program.				!!
		!!																		!!
		!!          if none of the above choices (with itol, rtol, and atol		!!
		!!          fixed throughout the problem) is suitable, more general		!!
		!!          error controls can be obtained by substituting				!!
		!!          user-supplied routines for the setting of ewt and/or for	!!
		!!          the norm calculation.  see part iv below.					!!
		!!																		!!
		!!          if global errors are to be estimated by making a repeated	!!
		!!          run on the same problem with smaller tolerances, then all	!!
		!!          components of rtol and atol (i.e. of ewt) should be scaled	!!
		!!          down uniformly.												!!
		itol   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! This is a very stiff problem, as aerosol and the gas phase should 
		!! not be the same scale at all.  Screw absolute tolerance.
		!! Could also explore using a vector of atol's
		!!
		!! However, if the initial concentration is zero, some absolute tolerance
		!! is necessary (an error of 0 is unacceptable).
		!!
		!! The units of this are molecules
		atol   = 1.0e-9

		!Relative Error Tolerance
		rtol   = (AqThermoEquilibriumError/10000.) **(1./NumTimeSteps) ! AqThermoEquilibriumError
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! This is the maximum number of internal steps allowed. !!
		!! The default is 500 and that's not enough.			 !!
		iwork(6) = 300000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Tell the program which gridpoint it's operating at: !!
		!WRITE(*,*) "Location: ", XP, YP, ZP
		Y(6*HowManyLinks + 5) = XP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Y(6*HowManyLinks + 6) = YP
		Y(6*HowManyLinks + 7) = ZP
	
		!! Reset timestep to the correct position
		T = (TimeStep - 1) * TimeStepSize

		!! IF the flag to equilibrate ammonia has been selected, 
		!! call error
		IF (AmmoniaFlux .EQ. .FALSE. ) THEN
			CALL ERROR("Ammonia must be done in flux mode in BigStepInorgCondensation()")
		END IF
		
		!! Find the indices of the reactions
		DO I = 1, HowManyDissolutionReactions
			IF (DissolutionData(I,1) .EQ. FindChem("H2SO4",0)) THEN
				Y(8 + 16*HowManyLinks) = I
			ELSE IF (DissolutionData(I,1) .EQ. FindChem("HNO3",0)) THEN
				Y(9 + 16*HowManyLinks) = I
			ELSE IF (DissolutionData(I,1) .EQ. FindChem("HCl",0)) THEN
				Y(10 + 16*HowManyLinks) = I
			ELSE IF (DissolutionData(I,1) .EQ. FindChem("NH3",0)) THEN
				IF (DissolutionData(I,3) .EQ. FindChem("NH3",1)) THEN
					Y(12 + 16*HowManyLinks) = I
				ELSEIF (DissolutionData(I,3) .EQ. FindChem("NH4+",1)) THEN
					Y(11 + 16*HowManyLinks) = I
				END IF
			END IF
		END DO		
	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICALS										!!
		!! The Y-Vector we will integrate with and push through LSODES looks like this:			!!
	 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!    INTEGRATED TERMS (1 to 4 + 6*HowManyLinks)										!!
		!!    (1)								            : H2SO4 Gas Phase Concentration		!!
		!!    (2)								            : HNO3 Gas Phase Concentration		!!
		!!    (3)								            : HCl Gas Phase Concentration		!!
		!!    (4)								            : NH3 Gas Phase Concentration		!!	
		!!    (5 to 4 + HowManyLinks)				        : SO4-- Aerosol Concentrations		!!
		!!    (5 + HowManyLinks to 4 + 2*HowManyLinks)      : NO3- Aerosol Concentrations		!!
		!!    (5 + 2*HowManyLinks to 4 + 3*HowManyLinks)    : Cl- Aerosol Concentrations		!!
		!!    (5 + 3*HowManyLinks to 4 + 4*HowManyLinks)    : NH4+ Aerosol Concentrations		!!
		!!    (5 + 4*HowManyLinks to 4 + 5*HowManyLinks)    : NH3(aq) Aerosol Concentrations	!!
		!!    (5 + 5*HowManyLinks to 4 + 6*HowManyLinks)    : H+ Aerosol Concentrations			!!
		!!																						!!
		!!	  PARAMETER TERMS																	!!
		!!	  (5 + 6*HowManyLinks to ??)														!!
		!!																						!!
		!!    (5 + 6*HowManyLinks to 7 + 6*HowManyLinks)    : Location (XP, YP, ZP)				!!
		!!    (8 + 7*HowManyLinks to 7 + 8*HowManyLinks)    : (S/H) for HNO3(g) <=> H+ & NO3-	!!
		!!    (8 + 8*HowManyLinks to 7 + 9*HowManyLinks)    : (S/H) for HCl(g) <=> H+ & Cl-		!!
		!!    (8 + 9*HowManyLinks to 7 + 10*HowManyLinks)   : (S/H) for NH3(g) & H+(aq) <=> NH4+!!
		!!    (8 + 10*HowManyLinks to 7 + 11*HowManyLinks)  : (S/H) for NH3(g) <=> NH3 (aq)		!!
		!!	  (8 + 11*HowManyLinks to 7 + 12*HowManyLinks)  : NumberofParticles per bin			!!
		!!	  (8 + 12*HowManyLinks to 7 + 13*HowManyLinks)  : H2SO4 Rate Coefficient			!!
		!!	  (8 + 13*HowManyLinks to 7 + 14*HowManyLinks)  : HNO3 Rate Coefficient				!!
		!!	  (8 + 14*HowManyLinks to 7 + 15*HowManyLinks)  : HCl Rate Coefficient				!!
		!!	  (8 + 15*HowManyLinks to 7 + 16*HowManyLinks)  : NH3 Rate Coefficient				!!
		!!    (8 + 16*HowManyLinks)							: H2SO4 Rxn Index					!!
		!!    (9 + 16*HowManyLinks)							: HNO3 Rxn Index					!!
		!!    (10 + 16*HowManyLinks)						: HCl Rxn Index						!!
		!!    (11 + 16*HowManyLinks)						: NH3(g) & H+ <=> NH4+ Rxn Index	!!
		!!    (12 + 16*HowManyLinks)						: NH3(g) <=> NH3(aq) Rxn Index		!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! HowManyLinks is the number of (non-zero) sections in the size dist					!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Initialize the Gas Concentrations from the local chemical concentrations: !! 
		Y(1) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(8 + 16*HowManyLinks),1)))
		Y(2) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(9 + 16*HowManyLinks),1)))
		Y(3) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(10 + 16*HowManyLinks),1)))
		Y(4) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(11 + 16*HowManyLinks),1)))

		WRITE(*,*) "In H2SO4", Y(1), Y(1)*1e9/GetM (XX(1),XX(2),XX(3)), GetM (XX(1),XX(2),XX(3))
		WRITE(*,*) "In HNO3", Y(2), Y(2)*1e9/GetM (XX(1),XX(2),XX(3)), GetM (XX(1),XX(2),XX(3))
		WRITE(*,*) "In HCl", Y(3), Y(3)*1e9/GetM (XX(1),XX(2),XX(3)), GetM (XX(1),XX(2),XX(3))	
		WRITE(*,*) "In NH3", Y(4), Y(4)*1e9/GetM (XX(1),XX(2),XX(3)), GetM (XX(1),XX(2),XX(3))
		!! The condensation rates
		Y(12*HowManyLinks+8:13*HowManyLinks+7) = CondensationRateCoefficients (XP, YP, ZP, INT(Y(8 + 16*HowManyLinks)), HowManyLinks)
		Y(13*HowManyLinks+8:14*HowManyLinks+7) = CondensationRateCoefficients (XP, YP, ZP, INT(Y(9 + 16*HowManyLinks)), HowManyLinks)
		Y(14*HowManyLinks+8:15*HowManyLinks+7) = CondensationRateCoefficients (XP, YP, ZP, INT(Y(10 + 16*HowManyLinks)), HowManyLinks)
		Y(15*HowManyLinks+8:16*HowManyLinks+7) = CondensationRateCoefficients (XP, YP, ZP, INT(Y(11 + 16*HowManyLinks)), HowManyLinks)

		IF (Scaffolding) CALL TRANSCRIPT("Cond. Check 2")
		!! Specify Ia and Ja
		DO I = 1, NEQ
			IF (I .LE. 4) THEN !Gases
				CaseFlag = I
			ELSEIF (I .GT. 4 .AND. I .LE. 4 + HowManyLinks) THEN !SO4--
				CaseFlag = 5
			ELSE IF (I .GT. 4 + HowManyLinks .AND. I .LE. 4 + 2*HowManyLinks) THEN
				CaseFlag = 6	 
			ELSE IF (I .GT. 4 + 2*HowManyLinks .AND. I .LE. 4 + 3*HowManyLinks) THEN
				CaseFlag = 7	 
			ELSE IF (I .GT. 4 + 3*HowManyLinks .AND. I .LE. 4 + 4*HowManyLinks) THEN
				CaseFlag = 8	 
			ELSE IF (I .GT. 4 + 4*HowManyLinks .AND. I .LE. 4 + 5*HowManyLinks) THEN
				CaseFlag = 9	 
			ELSE IF (I .GT. 4 + 5*HowManyLinks .AND. I .LE. 4 + 6*HowManyLinks) THEN
				CaseFlag = 10	 
			END IF

			SELECT CASE (CaseFlag)

			!d/d(H2SO4(g))
			Case(1)
				
				!Ia for this column
				IWORK(30+I) = 1

				!Ja's for this column
				JaIndex = 31 + NEQ !This keeps track of the index for 
								   !Ja in IWORK
				!d(H2SO4(g))/dt
				IWORK(JaIndex + 1) = 1
				JaIndex = JaIndex + 1

				!d(SO4--)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + J
					JaIndex = JaIndex + 1
				END DO	

				!d(H+)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
					JaIndex = JaIndex + 1
				END DO
			
			!d/d(HNO3(g))
			CASE(2)
				
				!Ia for this column
				IWORK(30+I) = 2

				!Ja's for this column
				!d(HNO3(g))/dt
				IWORK(JaIndex + 1) = 2
				JaIndex = JaIndex + 1
				
				!d(HCl(g))/dt
				IWORK(JaIndex + 1) = 3
				JaIndex = JaIndex + 1
				
				!d(NH3(g))/dt
				IWORK(JaIndex + 1) = 4
				JaIndex = JaIndex + 1

				!d(NO3-)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + HowManyLinks +J
					JaIndex = JaIndex + 1
				END DO	

				!d(Cl-)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 2*HowManyLinks +J
					JaIndex = JaIndex + 1
				END DO	

				!d(NH4+)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 3*HowManyLinks + J
					JaIndex = JaIndex + 1
				END DO

				!d(H+)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
					JaIndex = JaIndex + 1
				END DO
			
			!d/d(HCl(g))
			CASE(3)
				
				!Ia for this column
				IWORK(30+I) = 2

				!Ja's for this column
				!d(HNO3(g))/dt
				IWORK(JaIndex + 1) = 2
				JaIndex = JaIndex + 1
				
				!d(HCl(g))/dt
				IWORK(JaIndex + 1) = 3
				JaIndex = JaIndex + 1
				
				!d(NH3(g))/dt
				IWORK(JaIndex + 1) = 4
				JaIndex = JaIndex + 1

				!d(NO3-)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + HowManyLinks +J
					JaIndex = JaIndex + 1
				END DO	

				
				!d(Cl-)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 2*HowManyLinks +J
					JaIndex = JaIndex + 1
				END DO	

				!d(NH4+)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 3*HowManyLinks + J
					JaIndex = JaIndex + 1
				END DO

				!d(H+)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
					JaIndex = JaIndex + 1
				END DO

			!d/d(NH3(g))
			CASE(4)
				
				!Ia for this column
				IWORK(30+I) = 2

				!Ja's for this column
				!d(HNO3(g))/dt
				IWORK(JaIndex + 1) = 2
				JaIndex = JaIndex + 1
				
				!d(HCl(g))/dt
				IWORK(JaIndex + 1) = 3
				JaIndex = JaIndex + 1
				
				!d(NH3(g))/dt
				IWORK(JaIndex + 1) = 4
				JaIndex = JaIndex + 1

				!d(NO3-)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + HowManyLinks +J
					JaIndex = JaIndex + 1
				END DO	

				!d(Cl-)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 2*HowManyLinks +J
					JaIndex = JaIndex + 1
				END DO	

				!d(NH4+)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 3*HowManyLinks + J
					JaIndex = JaIndex + 1
				END DO

				!d(NH3(aq))/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 4*HowManyLinks + J
					JaIndex = JaIndex + 1
				END DO
				
				!d(H+)/dt
				DO J = 1, HowManyLinks
					IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
					JaIndex = JaIndex + 1
				END DO
			
			!d/d(SO4--)
			CASE(5)
				J = I - 4

				!Ia for this column
				IWORK(30+I) = 4

				!Ja's for this column
				!d(NH3(g))/dt
				IWORK(JaIndex + 1) = 4
				JaIndex = JaIndex + 1
			
				!d(NH4+)/dt
				IWORK(JaIndex + 1) = 4 + 3*HowManyLinks + J
				JaIndex = JaIndex + 1
				
				!d(H+)/dt
				IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
				JaIndex = JaIndex + 1
				
			!d/d(NO3-)
			CASE(6)
				J = I - (4 + HowManyLinks)

				!Ia for this column
				IWORK(30+I) = 2

				!Ja's for this column
				!d(HNO3(g))/dt
				IWORK(JaIndex + 1) = 2
				JaIndex = JaIndex + 1
			
				!d(NO3-)/dt
				IWORK(JaIndex + 1) = 4 + HowManyLinks +J
				JaIndex = JaIndex + 1
	
				!d(H+)/dt
				IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
				JaIndex = JaIndex + 1

			!d/d(Cl-)
			CASE(7)
				J = I - (4 +2*HowManyLinks)

				!Ia for this column
				IWORK(30+I) = 3

				!Ja's for this column			
				!d(HCl(g))/dt
				IWORK(JaIndex + 1) = 3
				JaIndex = JaIndex + 1

				!d(Cl-)/dt
				IWORK(JaIndex + 1) = 4 + 2*HowManyLinks + J
				JaIndex = JaIndex + 1
				
				!d(H+)/dt
				IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
				JaIndex = JaIndex + 1
			
			!d/d(NH4+)
			CASE(8)
				J = I - (4 +3*HowManyLinks)

				!Ia for this column
				IWORK(30+I) = 4
			
				!d(NH3(g))/dt
				IWORK(JaIndex + 1) = 4
				JaIndex = JaIndex + 1
	
				!d(NH4+)/dt
				IWORK(JaIndex + 1) = 4 + 3*HowManyLinks + J
				JaIndex = JaIndex + 1
				
				!d(H+)/dt
				IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
				JaIndex = JaIndex + 1
				
			!d/d(NH3)
			CASE(9)
				J = I - (4 +4*HowManyLinks)

				!Ia for this column
				IWORK(30+I) = 4

				!Ja's for this column
				!d(NH3(g))/dt
				IWORK(JaIndex + 1) = 4
				JaIndex = JaIndex + 1

				!d(NH3(aq))/dt
				IWORK(JaIndex + 1) = 4 + 4*HowManyLinks + J
				JaIndex = JaIndex + 1
			
			!d/d(H+)
			CASE(10)
				J = I - (4 +5*HowManyLinks)

				!Ia for this column
				IWORK(30+I) = 2

				!Ja's for this column
				!d(HNO3(g))/dt
				IWORK(JaIndex + 1) = 2
				JaIndex = JaIndex + 1

				!d(HCl(g))/dt
				IWORK(JaIndex + 1) = 3
				JaIndex = JaIndex + 1
			
				!d(NH3(g))/dt
				IWORK(JaIndex + 1) = 4
				JaIndex = JaIndex + 1

				!d(NO3-)/dt
				IWORK(JaIndex + 1) = 4 + HowManyLinks +J
				JaIndex = JaIndex + 1

				!d(Cl-)/dt
				IWORK(JaIndex + 1) = 4 + 2*HowManyLinks + J
				JaIndex = JaIndex + 1
	
				!d(NH4+)/dt
				IWORK(JaIndex + 1) = 4 + 3*HowManyLinks + J
				JaIndex = JaIndex + 1
				
				!d(H+)/dt
				IWORK(JaIndex + 1) = 4 + 5*HowManyLinks + J
				JaIndex = JaIndex + 1

			END SELECT
		END DO

		!Set Y vector for Aerosol Concentration (molecules/particle)		
		SO4Index = FindChem("SO4--", 1)
		NO3Index = FindChem("NO3-", 1)
		ClIndex = FindChem("Cl-",1)
		NH4Index = FindChem("NH4+",1)
		NH3Index = FindChem("NH3",1)
		HIndex = FindChem("H+",1)
		
		CurrentParticle => FirstParticle
		DO I = 1, HowManyLinks  !! Loop over the aerosol
					
			Y(4 + I) = CurrentParticle%AqChems(SO4Index)*Avogadro
			Y(4 + HowManyLinks +I) = CurrentParticle%AqChems(NO3Index)*Avogadro
			Y(4 + 2*HowManyLinks +I) = CurrentParticle%AqChems(ClIndex)*Avogadro
			Y(4 + 3*HowManyLinks +I) = CurrentParticle%AqChems(NH4Index)*Avogadro
			Y(4 + 4*HowManyLinks +I) = CurrentParticle%AqChems(NH3Index)*Avogadro
			Y(4 + 5*HowManyLinks +I) = CurrentParticle%AqChems(HIndex)*Avogadro
			CurrentParticle => CurrentParticle%Next	
		END DO 
		
		IF (Scaffolding) CALL TRANSCRIPT("Cond. Check 3")
		!! Add in SoverH values to Y vector
		!HNO3 Rxn
		CALL DissolutionFactors (XP, YP, ZP, INT(Y(9 + 16*HowManyLinks)), &
								HowManyLinks, Temperature, & 
								INT(DissolutionData(Y(9 + 16*HowManyLinks),11)),	&
							   Dum1, Dum2, Dum3, &
							   Y(8 + 7*HowManyLinks:7 + 8*HowManyLinks), Dum4, ReturnType)
		!HCl Rxn
		CALL DissolutionFactors (XP, YP, ZP, INT(Y(10 + 16*HowManyLinks)), &
								HowManyLinks, Temperature, & 
								INT(DissolutionData(Y(10 + 16*HowManyLinks),11)),	&
							   Dum1, Dum2, Dum3, &
							   Y(8 + 8*HowManyLinks:7 + 9*HowManyLinks), Dum4, ReturnType)

		!NH3(g) & H+ <=> NH4+ Rxn
		CALL DissolutionFactors (XP, YP, ZP, INT(Y(11 + 16*HowManyLinks)), &
								HowManyLinks, Temperature, & 
								INT(DissolutionData(Y(11 + 16*HowManyLinks),11)),	&
							   Dum1, Dum2, Dum3, &
							   Y(8 + 9*HowManyLinks:7 + 10*HowManyLinks), Dum4, ReturnType)

		!NH3(g) <=> NH3(aq) Rxn
		CALL DissolutionFactors (XP, YP, ZP, INT(Y(12 + 16*HowManyLinks)), &
								HowManyLinks, Temperature, & 
								INT(DissolutionData(Y(12 + 16*HowManyLinks),11)),	&
							   Dum1, Dum2, Dum3, &
							   Y(8 + 10*HowManyLinks:7 + 11*HowManyLinks), Dum4, ReturnType)

		IF (Scaffolding) CALL TRANSCRIPT("Cond. Check 4")
		!! Store the initial concentrations for later
		DO J = 1, NEQ
			InitialConcsArray(J) = Y(J)
			WRITE(*,*) J, Y(J)
		END DO
		DO J = NEQ+1, YLEN
			WRITE(*,*) J, Y(J)
		END DO

		!! Each step starts as if a new run
		istate = 1

		!IF (Scaffolding) 
		WRITE(*,*) "Before BigCondensation LSODES"
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
30		CALL LSODES (BigCondensationODEEvaluator,neq,y,t,outtime,itol,rtol,atol,itask, &
					 istate,iopt,rwork,lrw,iwork,liw,BigCondensationJacobianEvaluator,mf)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
		!IF (SCAFFOLDING) 
		WRITE(*,*) "After BigCondensation LSODES"
		!! If the maximum number of steps occurred, then warn the user (in Scaffolding
		!! mode) and push the integration back into LSODES until the appropriate end time 
		!! is reached.
		IF (ISTATE .EQ. -1 ) THEN !.AND. T < OUTTIME) THEN
			!IF (Scaffolding) &
			!CALL WARN("Warning! Nominal maximun number of steps in LSODES was exceeded in Condensation Integration for "	&
			!	     //TRIM(GasPhaseChemicalNames(GasChemIndex1))//" ("//trim(int2str(iwork(6)))//" steps)... Recycling")
			CALL TRANSCRIPT("LSODES past max number of steps. Sent back.")
			ISTATE = 2
			GOTO 30
		END IF

		!! An Error will occur if the array lengths were too small
		IF (IWORK(17) > LRW)  &
			CALL ERROR("The Given Length for LSODES' Working Array RWORK (accessing in StepGasChemistry())"//		&
					  " Should have been of size "//TRIM(INT2STR(iwork(17)))//" but was only "//					&
					  TRIM(INT2STR(LRW))//" instead.  Fix the assignment in EvolveGasChemistry's SetEvolveGasChemConstants().")

		IF (ISTATE .LT. -1) THEN
			WRITE(*,*) "Time: ", T
			WRITE(*,*) "OutTime: ", OUTTIME 
			DO J = 1, NEQ
				WRITE(*,*) J, Y(J)
			END DO
			DO J = NEQ+1, YLEN
				WRITE(*,*) J, Y(J)
			END DO
			CALL ERROR("In BigStepInorgCondensation(), LSODES encountered a problem during the integration and returned ISTATE = "//	&
					TRIM(INT2STR(ISTATE))//".  Please investigate.")	
		END IF
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Replace the grid point gas chemical concentrations !!
		!Skip replace if gas is to be held constant
		IF(INT(DissolutionData(Y(8 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
			CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(8 + 16*HowManyLinks),1)), Y(1))
		END IF
		IF(INT(DissolutionData(Y(9 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
			CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(9 + 16*HowManyLinks),1)), Y(2))
		END IF
		IF(INT(DissolutionData(Y(10 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
			CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(10 + 16*HowManyLinks),1)), Y(3))
		END IF
		IF(INT(DissolutionData(Y(11 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
			CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(11 + 16*HowManyLinks),1)), Y(4))
		END IF
		
		IF (Scaffolding) WRITE(*,*) "After Replace Gas Chems" 
		WRITE(*,*) "Y(1) :", Y(1)
		WRITE(*,*) "Y(2) :", Y(2)
		WRITE(*,*) "Y(3) :", Y(3)
		WRITE(*,*) "Y(4) :", Y(4)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Walk the linked list and replace the chemicals in each particle !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CurrentParticle => FirstParticle
		DO I = 1, HowManyLinks  !! Loop over the aerosol
					
			
			CurrentParticle%AqChems(SO4Index) = Y(4 + I)/Avogadro
			CurrentParticle%AqChems(NO3Index) = Y(4 + HowManyLinks +I)/Avogadro
			CurrentParticle%AqChems(ClIndex) = Y(4 + 2*HowManyLinks +I)/Avogadro
			CurrentParticle%AqChems(NH4Index) = Y(4 + 3*HowManyLinks +I)/Avogadro
			CurrentParticle%AqChems(NH3Index) = Y(4 + 4*HowManyLinks +I)/Avogadro
			CurrentParticle%AqChems(HIndex) = Y(4 + 5*HowManyLinks +I)/Avogadro 
	
			WRITE(*,*) "Before Replacement"
			DO J = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
				WRITE(*,*) trim(AqPhaseChemicalNames(J)), " : ", CurrentParticle%AqChems(J)
			END DO
				
			!Check for below 0 concentrations, recover conserving charge and mass
			
			!NO3- below 0
			IF(CurrentParticle%AqChems(NO3Index) .LT. 0.0) THEN
				Store = CurrentParticle%AqChems(NO3Index)
				IF(CurrentParticle%AqChems(FindChem("NH4NO3",1)) .GT. ABS(Store)) THEN
					CurrentParticle%AqChems(NO3Index) = 0.0
					CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) - Store
					CurrentParticle%AqChems(FindChem("NH4NO3",1))  = &
					      CurrentParticle%AqChems(FindChem("NH4NO3",1)) + Store
				ELSE
					CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) &
								+ CurrentParticle%AqChems(FindChem("NH4NO3",1))
					CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) &
								+ CurrentParticle%AqChems(FindChem("NH4NO3",1))
					CurrentParticle%AqChems(FindChem("NH4NO3",1))  = 0.0
					Store = CurrentParticle%AqChems(NO3Index)
					IF(CurrentParticle%AqChems(FindChem("KNO3",1)) .GT. ABS(Store)) THEN
						CurrentParticle%AqChems(NO3Index) = 0.0
						CurrentParticle%AqChems(FindChem("K+",1)) = CurrentParticle%AqChems(FindChem("K+",1)) - Store
						CurrentParticle%AqChems(FindChem("KNO3",1))  = &
					      CurrentParticle%AqChems(FindChem("KNO3",1)) + Store
					ELSE
						CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) &
								+ CurrentParticle%AqChems(FindChem("KNO3",1))
						CurrentParticle%AqChems(FindChem("K+",1)) = CurrentParticle%AqChems(FindChem("K+",1)) &
								+ CurrentParticle%AqChems(FindChem("KNO3",1))
						CurrentParticle%AqChems(FindChem("KNO3",1))  = 0.0
						Store = CurrentParticle%AqChems(NO3Index)
						
						IF(FindChem("NaNO3",1, SuppressError=.TRUE.) .NE. 0) THEN
						  IF(CurrentParticle%AqChems(FindChem("NaNO3",1)) .GT. ABS(Store)) THEN
							CurrentParticle%AqChems(NO3Index) = 0.0
							CurrentParticle%AqChems(FindChem("Na+",1)) = CurrentParticle%AqChems(FindChem("Na+",1)) - Store
							CurrentParticle%AqChems(FindChem("NaNO3",1))  = &
							CurrentParticle%AqChems(FindChem("NaNO3",1)) + Store
							Store = 0.0
						  ELSE
						    CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) &
								+ CurrentParticle%AqChems(FindChem("NaNO3",1))
						    CurrentParticle%AqChems(FindChem("Na+",1)) = CurrentParticle%AqChems(FindChem("Na+",1)) &
								+ CurrentParticle%AqChems(FindChem("NaNO3",1))
						    CurrentParticle%AqChems(FindChem("NaNO3",1))  = 0.0
						    Store = CurrentParticle%AqChems(NO3Index)
						  END IF
						END IF
						
						IF(FindChem("Mg(NO3)2",1, SuppressError=.TRUE.) .NE. 0) THEN
						  IF(CurrentParticle%AqChems(FindChem("Mg(NO3)2",1)) .GT. 0.5*ABS(Store)) THEN
							CurrentParticle%AqChems(NO3Index) = 0.0
							CurrentParticle%AqChems(FindChem("Mg++",1)) = CurrentParticle%AqChems(FindChem("Mg++",1)) - 0.5*Store
							CurrentParticle%AqChems(FindChem("Mg(NO3)2",1))  = &
							CurrentParticle%AqChems(FindChem("Mg(NO3)2",1)) + 0.5*Store
						   	Store = 0.0
						  ELSE
						    CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) &
								+ 2*CurrentParticle%AqChems(FindChem("Mg(NO3)2",1))
						    CurrentParticle%AqChems(FindChem("Mg++",1)) = CurrentParticle%AqChems(FindChem("Mg++",1)) &
								+ CurrentParticle%AqChems(FindChem("Mg(NO3)2",1))
						    CurrentParticle%AqChems(FindChem("Mg(NO3)2",1))  = 0.0
						    Store = CurrentParticle%AqChems(NO3Index)
						  END IF
						END IF

						IF(FindChem("Ca(NO3)2",1, SuppressError=.TRUE.) .NE. 0) THEN
						  IF(CurrentParticle%AqChems(FindChem("Ca(NO3)2",1)) .GT. 0.5*ABS(Store)) THEN
							CurrentParticle%AqChems(NO3Index) = 0.0
							CurrentParticle%AqChems(FindChem("Ca++",1)) = CurrentParticle%AqChems(FindChem("Ca++",1)) - 0.5*Store
							CurrentParticle%AqChems(FindChem("Ca(NO3)2",1))  = &
							CurrentParticle%AqChems(FindChem("Ca(NO3)2",1)) + 0.5*Store
						    Store = 0.0
						  ELSE
						    CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) &
								+ 2*CurrentParticle%AqChems(FindChem("Ca(NO3)2",1))
						    CurrentParticle%AqChems(FindChem("Ca++",1)) = CurrentParticle%AqChems(FindChem("Ca++",1)) &
								+ CurrentParticle%AqChems(FindChem("Ca(NO3)2",1))
						    CurrentParticle%AqChems(FindChem("Ca(NO3)2",1))  = 0.0
						    Store = CurrentParticle%AqChems(NO3Index)
						  END IF
						END IF

						IF(Store .LT. 0.0) THEN
							CurrentParticle%AqChems(NO3Index) = 0.0
							CurrentParticle%AqChems(HIndex) = CurrentParticle%AqChems(HIndex) - Store 
							Y(2) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(9 + 16*HowManyLinks),1)))
							Y(2) = Y(2) + Store*Avogadro*Y(7 + 11*HowManyLinks + I)
							IF(Y(2) .LT. 0.0) CALL ERROR("HNO3(g) went below 0 in BigStepInorgCondensation.")
							IF(INT(DissolutionData(Y(9 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
								CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(9 + 16*HowManyLinks),1)), Y(2))
							END IF
						END IF
					END IF
				END IF
			END IF
			
			!Cl- below 0
			IF(CurrentParticle%AqChems(ClIndex) .LT. 0.0) THEN
				Store = CurrentParticle%AqChems(ClIndex)
				IF(CurrentParticle%AqChems(FindChem("NH4Cl",1)) .GT. ABS(Store)) THEN
					CurrentParticle%AqChems(ClIndex) = 0.0
					CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) - Store
					CurrentParticle%AqChems(FindChem("NH4Cl",1))  = &
					      CurrentParticle%AqChems(FindChem("NH4Cl",1)) + Store
				ELSE
					CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) &
								+ CurrentParticle%AqChems(FindChem("NH4Cl",1))
					CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) &
								+ CurrentParticle%AqChems(FindChem("NH4Cl",1))
					CurrentParticle%AqChems(FindChem("NH4Cl",1))  = 0.0
					Store = CurrentParticle%AqChems(ClIndex)
					IF(CurrentParticle%AqChems(FindChem("KCl",1)) .GT. ABS(Store)) THEN
						CurrentParticle%AqChems(ClIndex) = 0.0
						CurrentParticle%AqChems(FindChem("K+",1)) = CurrentParticle%AqChems(FindChem("K+",1)) - Store
						CurrentParticle%AqChems(FindChem("KCl",1))  = &
					      CurrentParticle%AqChems(FindChem("KCl",1)) + Store
					ELSE
						CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) &
								+ CurrentParticle%AqChems(FindChem("KCl",1))
						CurrentParticle%AqChems(FindChem("K+",1)) = CurrentParticle%AqChems(FindChem("K+",1)) &
								+ CurrentParticle%AqChems(FindChem("KCl",1))
						CurrentParticle%AqChems(FindChem("KCl",1))  = 0.0
						Store = CurrentParticle%AqChems(ClIndex)

						IF(FindChem("NaCl",1, SuppressError=.TRUE.) .NE. 0) THEN
						  IF(CurrentParticle%AqChems(FindChem("NaCl",1)) .GT. ABS(Store)) THEN
							CurrentParticle%AqChems(ClIndex) = 0.0
							CurrentParticle%AqChems(FindChem("Na+",1)) = CurrentParticle%AqChems(FindChem("Na+",1)) - Store
							CurrentParticle%AqChems(FindChem("NaCl",1))  = &
							CurrentParticle%AqChems(FindChem("NaCl",1)) + Store
							Store = 0.0
						  ELSE
						    CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) &
								+ CurrentParticle%AqChems(FindChem("NaCl",1))
						    CurrentParticle%AqChems(FindChem("Na+",1)) = CurrentParticle%AqChems(FindChem("Na+",1)) &
								+ CurrentParticle%AqChems(FindChem("NaCl",1))
						    CurrentParticle%AqChems(FindChem("NaCl",1))  = 0.0
						    Store = CurrentParticle%AqChems(ClIndex)
						  END IF
						END IF
						
						IF(FindChem("MgCl2",1, SuppressError=.TRUE.) .NE. 0) THEN
						  IF(CurrentParticle%AqChems(FindChem("MgCl2",1)) .GT. 0.5*ABS(Store)) THEN
							CurrentParticle%AqChems(ClIndex) = 0.0
							CurrentParticle%AqChems(FindChem("Mg++",1)) = CurrentParticle%AqChems(FindChem("Mg++",1)) - 0.5*Store
							CurrentParticle%AqChems(FindChem("MgCl2",1))  = &
							CurrentParticle%AqChems(FindChem("MgCl2",1)) + 0.5*Store
						   	Store = 0.0
						  ELSE
						    CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) &
								+ 2*CurrentParticle%AqChems(FindChem("MgCl2",1))
						    CurrentParticle%AqChems(FindChem("Mg++",1)) = CurrentParticle%AqChems(FindChem("Mg++",1)) &
								+ CurrentParticle%AqChems(FindChem("MgCl2",1))
						    CurrentParticle%AqChems(FindChem("MgCl2",1))  = 0.0
						    Store = CurrentParticle%AqChems(ClIndex)
						  END IF
						END IF

						IF(FindChem("CaCl2",1, SuppressError=.TRUE.) .NE. 0) THEN
						  IF(CurrentParticle%AqChems(FindChem("CaCl2",1)) .GT. 0.5*ABS(Store)) THEN
							CurrentParticle%AqChems(ClIndex) = 0.0
							CurrentParticle%AqChems(FindChem("Ca++",1)) = CurrentParticle%AqChems(FindChem("Ca++",1)) - 0.5*Store
							CurrentParticle%AqChems(FindChem("CaCl2",1))  = &
							CurrentParticle%AqChems(FindChem("CaCl2",1)) + 0.5*Store
						    Store = 0.0
						  ELSE
						    CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) &
								+ 2*CurrentParticle%AqChems(FindChem("CaCl2",1))
						    CurrentParticle%AqChems(FindChem("Ca++",1)) = CurrentParticle%AqChems(FindChem("Ca++",1)) &
								+ CurrentParticle%AqChems(FindChem("CaCl2",1))
						    CurrentParticle%AqChems(FindChem("CaCl2",1))  = 0.0
						    Store = CurrentParticle%AqChems(ClIndex)
						  END IF
						END IF

						IF(Store .LT. 0.0) THEN
							CurrentParticle%AqChems(ClIndex) = 0.0
						    CurrentParticle%AqChems(HIndex) = CurrentParticle%AqChems(HIndex) - Store 
							Y(3) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(10 + 16*HowManyLinks),1)))
							Y(3) = Y(3) + Store*Avogadro*Y(7 + 11*HowManyLinks + I)
							IF(Y(3) .LT. 0.0) CALL ERROR("HCl(g) went below 0 in BigStepInorgCondensation.")
							IF(INT(DissolutionData(Y(10 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
								CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(10 + 16*HowManyLinks),1)), Y(3))
							END IF
						END IF
					END IF
				END IF
			END IF
			
			IF(CurrentParticle%AqChems(NH3Index) .LT. 0.0) CurrentParticle%AqChems(NH3Index) = 0.0
			
			!NH4+ below 0
			IF(CurrentParticle%AqChems(NH4Index) .LT. 0.0) THEN
				Store = CurrentParticle%AqChems(NH4Index)
				IF(CurrentParticle%AqChems(FindChem("NH4Cl",1)) .GT. ABS(Store)) THEN
					CurrentParticle%AqChems(NH4Index) = 0.0
					CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) - Store
					CurrentParticle%AqChems(FindChem("NH4Cl",1))  = &
					      CurrentParticle%AqChems(FindChem("NH4Cl",1)) + Store
				ELSE
					CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) &
								+ CurrentParticle%AqChems(FindChem("NH4Cl",1))
					CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) &
								+ CurrentParticle%AqChems(FindChem("NH4Cl",1))
					CurrentParticle%AqChems(FindChem("NH4Cl",1))  = 0.0
					Store = CurrentParticle%AqChems(NH4Index)
					IF(CurrentParticle%AqChems(FindChem("NH4NO3",1)) .GT. ABS(Store)) THEN
						CurrentParticle%AqChems(NH4Index) = 0.0
						CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) - Store
						CurrentParticle%AqChems(FindChem("NH4NO3",1))  = &
					      CurrentParticle%AqChems(FindChem("NH4NO3",1)) + Store
					ELSE
						CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) &
								+ CurrentParticle%AqChems(FindChem("NH4NO3",1))
						CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) &
								+ CurrentParticle%AqChems(FindChem("NH4NO3",1))
						CurrentParticle%AqChems(FindChem("NH4NO3",1))  = 0.0
						Store = CurrentParticle%AqChems(NH4Index)
						IF(CurrentParticle%AqChems(FindChem("(NH4)2SO4",1)) .GT. 0.5*ABS(Store)) THEN
							CurrentParticle%AqChems(NH4Index) = 0.0
							CurrentParticle%AqChems(SO4Index) = CurrentParticle%AqChems(SO4Index) - 0.5*Store
							CurrentParticle%AqChems(FindChem("(NH4)2SO4",1))  = &
								CurrentParticle%AqChems(FindChem("(NH4)2SO4",1)) + 0.5*Store
						ELSE
							CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) &
									+ 2*CurrentParticle%AqChems(FindChem("(NH4)2SO4",1))
							CurrentParticle%AqChems(SO4Index) = CurrentParticle%AqChems(SO4Index) &
									+ CurrentParticle%AqChems(FindChem("(NH4)2SO4",1))
							CurrentParticle%AqChems(FindChem("(NH4)2SO4",1))  = 0.0
							Store = CurrentParticle%AqChems(NH4Index)
							IF(CurrentParticle%AqChems(FindChem("NH4HSO4",1)) .GT. ABS(Store)) THEN
								CurrentParticle%AqChems(NH4Index) = 0.0
								CurrentParticle%AqChems(FindChem("HSO4-",1)) = CurrentParticle%AqChems(FindChem("HSO4-",1)) - Store
								CurrentParticle%AqChems(FindChem("NH4HSO4",1))  = &
									CurrentParticle%AqChems(FindChem("NH4HSO4",1)) + Store
							ELSE
								CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) &
									+ CurrentParticle%AqChems(FindChem("NH4HSO4",1))
								CurrentParticle%AqChems(FindChem("HSO4-",1)) = CurrentParticle%AqChems(FindChem("HSO4-",1)) &
									+ CurrentParticle%AqChems(FindChem("NH4HSO4",1))
								CurrentParticle%AqChems(FindChem("NH4HSO4",1))  = 0.0
								Store = CurrentParticle%AqChems(NH4Index)
								CurrentParticle%AqChems(NH4Index) = 0.0
								CurrentParticle%AqChems(HIndex) = CurrentParticle%AqChems(HIndex) + Store 
								Y(4) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(11 + 16*HowManyLinks),1)))
								Y(4) = Y(4) + Store*Avogadro*Y(7 + 11*HowManyLinks + I)
								IF(Y(4) .LT. 0.0) CALL ERROR("NH3(g) went below 0 in BigStepInorgCondensation.")
								IF(INT(DissolutionData(Y(11 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
									CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(11 + 16*HowManyLinks),1)), Y(4))
								END IF
							END IF
						END IF
					END IF
				END IF
			END IF
			
			!H+ Below 0
			!Use NH4+, then HNO3(g), then HCl(g)
			WRITE(*,*) "H+ Conc: ", CurrentParticle%AqChems(HIndex)
			IF(CurrentParticle%AqChems(HIndex) .LT. 0.0) THEN
				CALL TRANSCRIPT("Use NH4+")
				Store = CurrentParticle%AqChems(HIndex)
				IF(CurrentParticle%AqChems(NH4Index) .GT. ABS(Store)) THEN
					CurrentParticle%AqChems(HIndex) = 0.0
					CurrentParticle%AqChems(NH3Index) = CurrentParticle%AqChems(NH3Index) - Store
					CurrentParticle%AqChems(NH4Index) = CurrentParticle%AqChems(NH4Index) + Store
				ELSE
					CALL TRANSCRIPT("Use HNO3(g)")
					CurrentParticle%AqChems(HIndex) = CurrentParticle%AqChems(HIndex) + CurrentParticle%AqChems(NH4Index)
					CurrentParticle%AqChems(NH3Index) = CurrentParticle%AqChems(NH3Index) + CurrentParticle%AqChems(NH4Index)
					CurrentParticle%AqChems(NH4Index) = 0.0
					Store = CurrentParticle%AqChems(HIndex)
					Y(2) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(9 + 16*HowManyLinks),1)))
					IF(Y(2) .GT. ABS(Store*Avogadro*Y(7 + 11*HowManyLinks + I))) THEN
						Y(2) = Y(2) + Store*Avogadro*Y(7 + 11*HowManyLinks + I)
						CurrentParticle%AqChems(HIndex) = CurrentParticle%AqChems(HIndex) - Store
						CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) - Store
						IF(Y(2) .LT. 0.0) CALL ERROR("HNO3(g) went below 0 in BigStepInorgCondensation.")
						IF(INT(DissolutionData(Y(9 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
							CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(9 + 16*HowManyLinks),1)), Y(2))
						END IF
					ELSE	
						CALL TRANSCRIPT("Use HCl(g)")
						CurrentParticle%AqChems(HIndex) = CurrentParticle%AqChems(HIndex) + Y(2)/Avogadro/Y(7 + 11*HowManyLinks + I)
						CurrentParticle%AqChems(NO3Index) = CurrentParticle%AqChems(NO3Index) + Y(2)/Avogadro/Y(7 + 11*HowManyLinks + I)
						Y(2) = 0.0
						IF(INT(DissolutionData(Y(9 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
							CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(9 + 16*HowManyLinks),1)), Y(2))
						END IF
						Store = CurrentParticle%AqChems(HIndex)
					
						Y(3) = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(10 + 16*HowManyLinks),1)))
						IF(Y(3) .GT. ABS(Store*Avogadro*Y(7 + 11*HowManyLinks + I))) THEN
							Y(3) = Y(3) + Store*Avogadro*Y(7 + 11*HowManyLinks + I)
							CurrentParticle%AqChems(HIndex) = CurrentParticle%AqChems(HIndex) - Store
							CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) - Store
							IF(Y(3) .LT. 0.0) CALL ERROR("HCl(g) went below 0 in BigStepInorgCondensation.")
							IF(INT(DissolutionData(Y(10 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
								CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(10 + 16*HowManyLinks),1)), Y(3))
							END IF
						ELSE	
							CurrentParticle%AqChems(HIndex) = CurrentParticle%AqChems(HIndex) + Y(3)/Avogadro/Y(7 + 11*HowManyLinks + I)
							CurrentParticle%AqChems(ClIndex) = CurrentParticle%AqChems(ClIndex) + Y(3)/Avogadro/Y(7 + 11*HowManyLinks + I)
							Y(3) = 0.0
							IF(INT(DissolutionData(Y(10 + 16*HowManyLinks),1)) .LE. HowManyEvolveGasChems) THEN
								CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Y(10 + 16*HowManyLinks),1)), Y(3))
							END IF
						END IF
					END IF
				END IF
			END IF



			!! Re-Equilibrate
			!IF (Scaffolding) 
			WRITE(*,*) "Before Electrolyte Eq"				
			WRITE(*,*) "H+: ", CurrentParticle%AqChems(ProtonIndex)
			CALL FindElectrolyteEquilibrium(CurrentParticle, UpdateThermo = .TRUE., &
                                                        FirstEquilibration=.FALSE.,ReturnType=GG)
			!IF (Scaffolding) WRITE(*,*) "Electrolyte Eq Time: ", AmmoniaTime

			!WRITE(*,*) "Cl- :   ", CurrentParticle%AqChems(ClIndex)*CurrentParticle%NumberofParticles*Avogadro, &
			!			"     HCl(g):   ", Y(3)
			CurrentParticle => CurrentParticle%Next

		END DO !! Quit loop over particles for replacement

	END DO !! Stop Looping Over the Time Steps
	!WRITE(*,*) "Outside LSODES loop"
	!! Get rid of everything that was allocated in here.
	DEALLOCATE (Y,RWORK,IWORK,InitialConcsArray,ChemStorageArray, Dum1, Dum2, Dum3, Dum4)
		
	!! If relative humidity is below 20%, and aerosol is not metastable,
	!! equilibrate gas-solid reactions for particles since we didn't
	!! condense for them (MJA, 07/12/07)
	IF (RH .LT. 0.20 .AND. NOT(METASTABLE)) THEN
		Equilibrated = .FALSE.
		DO WHILE (.NOT.Equilibrated)
		K = K + 1
		Equilibrated = .TRUE.

		!! Loop over the reactions
		DO I = 1, HowManyDissolutionReactions
			!Only gas-solid rxns				
			IF (DissolutionData(I,11) .EQ. 6.) THEN
				!WRITE(*,*) I
				!Get the gas phase concentration from the grid
				GasPhaseBurden = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(I,1)))
				GasPhaseBurden2 = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(I,2)))
			
				!Equilibrate gas-solid reactions (Type 6) using surface area to remove non-uniqueness problem
				CALL EquilibrateGasSolidReactions(XP, YP, ZP, GasPhaseBurden, I, Temperature, ReturnType, &
							GasPhaseBurden2 = GasPhaseBurden2, OptErrorTolerance=AqThermoEquilibriumError/2.)
		     
                                IF (ReturnType .NE. 2) THEN
					Equilibrated    = .FALSE.
				END IF
				
				!! Replace the gas phase chemical concentrations
				IF (INT(DissolutionData(I,1)) .LE. HowManyEvolveGasChems) THEN
					CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(I,1)), GasPhaseBurden)
				END IF

				IF (INT(DissolutionData(I,2)) .LE. HowManyEvolveGasChems .AND. DissolutionData(I,11) .EQ. 6.)	&
					CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(I,2)), GasPhaseBurden2)
			
			END IF
		END DO



		IF (K .GT. AqThermoNumbAllowedIterations) THEN
			!CALL DumpParticleContentsAtError(Particles(XP,YP,ZP)%First, "UnEquilableParticle.txt") 
			!CALL WARN ("BigStepInorgCondensation in CondensationIntegrator.h tried to equilibrate all of the gas-solid "// &
			!		    "reactions and failed. You may try to fix this by increasing AqThermoNumbAllowedIterations,"// &
			!			" the maximum  allowed iterations of this type of routine.")
			EXIT
		END IF
		
		END DO !Do while
	END IF

	!! SCSCSCSC
	IF (Scaffolding) CALL TRANSCRIPT (">>>Exiting BigStepInorgCondensation<<<")
	!WRITE(*,*) "StepCondensation took "//trim(real2str(CondTime))//" seconds."
		
	RETURN
END SUBROUTINE BigStepInorgCondensation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This program provides the ODE ydot vector to LSODES. !!
!! Rewritten for coupled inorganic case by Matt Alvarado!!
!! 10/04/06												!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BigCondensationODEEvaluator (neq, t, y, ydot)

		USE Aerosols,            ONLY : Particles, Particle, ParticleArray, &
										FindElectrolyteEquilibrium, &
										EquilibrateAqueousDissociationReaction, &
										KusikMeissner, AerosolPH, &
										DumpParticleContentsAtError, &
										EqCoeff, Molality
		USE Chemistry,			 ONLY : HowManyGasChems, FindChem, &
										HowManyAqEqReactions, &
										AqPhaseChemicalNames, &
										AqEquilibriaList
		USE InfrastructuralCode, ONLY : Transcript, ERROR, REAL2STR, INT2STR
		USE ModelParameters,     ONLY : Avogadro, AqThermoEquilibriumError, &
										Metastable, HydroxyIndex, WaterMolecMass
		USE GridPointFields,	 ONLY : LocateGridPoint,				&
										getM,							&
										GasPhaseChemGridIndex,			&
										GetTemp, GetRelativeHumidity

		IMPLICIT NONE

		!! External Variable (prescribed by LSODES)
		REAL*8		:: t, y, ydot
		INTEGER		:: neq
		DIMENSION	:: y(1), ydot(1)

		!! Internal Variables
		INTEGER :: I, J, C, HowManyLinks, XP, YP, ZP, allocation_error
		INTEGER :: H2SO4RxnIndex, HNO3RxnIndex, HClRxnIndex
		INTEGER ::  NH3RxnIndex1, NH3RxnIndex2 
		Integer :: 	KNO3i, KCli, KHSO4i, &
				Nai, NaNO3i, NaCli, NaHSO4i, Na2SO4i, &
				NH3i, NH4i, NH4NO3i, NH4Cli, NH4HSO4i, NH42SO4i, &
				Cli, NO3i, SO4i, HSO4i, K2SO4i, Ki, Hi, OHi, Q, &
				MgNO32i, MgCl2i, CaCl2i, CaNO32i
		
		REAL*8				:: Temperature, RH, XX(3), Change, Store
		REAL*8, ALLOCATABLE :: H2SO4Coeff(:) ! 1/s
		REAL*8, ALLOCATABLE :: HNO3Coeff(:) ! 1/s
		REAL*8, ALLOCATABLE :: HClCoeff(:) ! 1/s
		REAL*8, ALLOCATABLE :: NH3Coeff(:) ! 1/s
		REAL*8, ALLOCATABLE :: NumParticlesPerBin(:) ! particles/cm3
		REAL*8, ALLOCATABLE :: HNO3Sat(:), HClSat(:), NH3Sat(:) !molecules/cm3
		INTEGER, ALLOCATABLE :: HNO3Flag(:), HClFlag(:), NH3Flag(:)
		REAL*8, ALLOCATABLE :: HNO3SoverH(:), HClSoverH(:)
		REAL*8, ALLOCATABLE :: NH3SoverH(:), NH3SoverH2(:) 
		REAL*8, ALLOCATABLE :: HClJac(:), HNO3Jac(:), NH3Jac(:)
		LOGICAL, ALLOCATABLE :: ClNearZero(:), NO3NearZero(:)

		LOGICAL				:: EqFlag, WriteFlag, SkipCl, SkipNO3
		LOGICAL, PARAMETER :: Scaffolding =  .FALSE.
		REAL*8 :: ClSum, NO3Sum, Dum1, Dum2

		TYPE(Particle), POINTER :: CurrentParticle
		
		IF (SCAFFOLDING) print *, ">>> In BigCondensationODEEvaluator <<<"

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICALS										!!
		!! The Y-Vector we will integrate with and push through LSODES looks like this:			!!
	 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!    INTEGRATED TERMS (1 to 4 + 6*HowManyLinks)										!!
		!!    (1)								            : H2SO4 Gas Phase Concentration		!!
		!!    (2)								            : HNO3 Gas Phase Concentration		!!
		!!    (3)								            : HCl Gas Phase Concentration		!!
		!!    (4)								            : NH3 Gas Phase Concentration		!!	
		!!    (5 to 4 + HowManyLinks)				        : SO4-- Aerosol Concentrations		!!
		!!    (5 + HowManyLinks to 4 + 2*HowManyLinks)      : NO3- Aerosol Concentrations		!!
		!!    (5 + 2*HowManyLinks to 4 + 3*HowManyLinks)    : Cl- Aerosol Concentrations		!!
		!!    (5 + 3*HowManyLinks to 4 + 4*HowManyLinks)    : NH4+ Aerosol Concentrations		!!
		!!    (5 + 4*HowManyLinks to 4 + 5*HowManyLinks)    : NH3(aq) Aerosol Concentrations	!!
		!!    (5 + 5*HowManyLinks to 4 + 6*HowManyLinks)    : H+ Aerosol Concentrations			!!
		!!																						!!
		!!	  PARAMETER TERMS																	!!
		!!	  (5 + 6*HowManyLinks to ??)														!!
		!!																						!!
		!!    (5 + 6*HowManyLinks to 7 + 6*HowManyLinks)    : Location (XP, YP, ZP)				!!
		!!    (8 + 7*HowManyLinks to 7 + 8*HowManyLinks)    : (S/H) for HNO3(g) <=> H+ & NO3-	!!
		!!    (8 + 8*HowManyLinks to 7 + 9*HowManyLinks)    : (S/H) for HCl(g) <=> H+ & Cl-		!!
		!!    (8 + 9*HowManyLinks to 7 + 10*HowManyLinks)   : (S/H) for NH3(g) & H+(aq) <=> NH4+!!
		!!    (8 + 10*HowManyLinks to 7 + 11*HowManyLinks)  : (S/H) for NH3(g) <=> NH3 (aq)		!!
		!!	  (8 + 11*HowManyLinks to 7 + 12*HowManyLinks)  : NumberofParticles per bin			!!
		!!	  (8 + 12*HowManyLinks to 7 + 13*HowManyLinks)  : H2SO4 Rate Coefficient			!!
		!!	  (8 + 13*HowManyLinks to 7 + 14*HowManyLinks)  : HNO3 Rate Coefficient				!!
		!!	  (8 + 14*HowManyLinks to 7 + 15*HowManyLinks)  : HCl Rate Coefficient				!!
		!!	  (8 + 15*HowManyLinks to 7 + 16*HowManyLinks)  : NH3 Rate Coefficient				!!
		!!    (8 + 16*HowManyLinks)							: H2SO4 Rxn Index					!!
		!!    (9 + 16*HowManyLinks)							: HNO3 Rxn Index					!!
		!!    (10 + 16*HowManyLinks)						: HCl Rxn Index						!!
		!!    (11 + 16*HowManyLinks)						: NH3(g) & H+ <=> NH4+ Rxn Index	!!
		!!    (12 + 16*HowManyLinks)						: NH3(g) <=> NH3(aq) Rxn Index		!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! HowManyLinks is the number of (non-zero) sections in the size dist					!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		HowManyLinks = (neq - 4)/6
				
		!Allocate parameter arrays
		ALLOCATE (H2SO4Coeff(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of H2SO4RateCoeff could not proceed in CondensationODEEvaluator")
		ALLOCATE (HNO3Coeff(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3RateCoeff could not proceed in CondensationODEEvaluator")
		ALLOCATE (HClCoeff(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HClRateCoeff could not proceed in CondensationODEEvaluator")
		ALLOCATE (NH3Coeff(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of NH3RateCoeff could not proceed in CondensationODEEvaluator")
		
		ALLOCATE (NumParticlesPerBin(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of NumParticlesPerBin could not proceed in CondensationODEEvaluator")
		
		ALLOCATE (HNO3Flag(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3Flag could not proceed in CondensationODEEvaluator")
		ALLOCATE (HClFlag(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HClFlag could not proceed in CondensationODEEvaluator")
		ALLOCATE (NH3Flag(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of NH3Flag could not proceed in CondensationODEEvaluator")

		ALLOCATE (HNO3Sat(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3Sat could not proceed in CondensationODEEvaluator")
		ALLOCATE (HClSat(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HClSat could not proceed in CondensationODEEvaluator")
		ALLOCATE (NH3Sat(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of NH3Sat could not proceed in CondensationODEEvaluator")
			
		ALLOCATE (HNO3SoverH(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3SoverH could not proceed in CondensationODEEvaluator")
		ALLOCATE (HClSoverH(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HClSoverH could not proceed in CondensationODEEvaluator")
		ALLOCATE (NH3SoverH(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of NH3SoverH could not proceed in CondensationODEEvaluator")
		ALLOCATE (NH3SoverH2(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of NH3SoverH2 could not proceed in CondensationODEEvaluator")
		
		ALLOCATE (ClNearZero(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of ClNearZero could not proceed in CondensationODEEvaluator")
		ALLOCATE (NO3NearZero(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of NO3NearZero could not proceed in CondensationODEEvaluator")
	
		ALLOCATE (HNO3Jac(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3Sat could not proceed in CondensationODEEvaluator")
		ALLOCATE (HClJac(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of HClSat could not proceed in CondensationODEEvaluator")
		ALLOCATE (NH3Jac(HowManyLinks), stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of NH3Sat could not proceed in CondensationODEEvaluator")

		
		XP = INT(Y(5 + 6*HowManyLinks))
		YP = INT(Y(6 + 6*HowManyLinks))
		ZP = INT(Y(7 + 6*HowManyLinks))

		XX          = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
		Temperature = GetTemp(XX(1),XX(2),XX(3))
		RH =GetRelativeHumidity(XX(1),XX(2),XX(3))
		WriteFlag = .FALSE. 
	
		H2SO4RxnIndex = Y(8 + 16*HowManyLinks)
		HNO3RxnIndex = Y(9 + 16*HowManyLinks)
		HClRxnIndex = Y(10 + 16*HowManyLinks)
		NH3RxnIndex1 = Y(11 + 16*HowManyLinks)
		NH3RxnIndex2 = Y(12 + 16*HowManyLinks)
			
		!Find chemical indices for later
		Ki = FindChem("K+", 1)
		KNO3i = FindChem("KNO3", 1)
		KCli = FindChem("KCl", 1)
		KHSO4i = FindChem("KHSO4", 1)
		K2SO4i = FindChem("K2SO4", 1)
		
		Nai = FindChem("Na+", 1, SuppressError=.TRUE.)
		NaNO3i = FindChem("NaNO3", 1, SuppressError=.TRUE.)
		NaCli = FindChem("NaCl", 1, SuppressError=.TRUE.)
		NaHSO4i = FindChem("NaHSO4", 1, SuppressError=.TRUE.)
		Na2SO4i = FindChem("Na2SO4", 1, SuppressError=.TRUE.)
		
		MgNO32i = FindChem("Mg(NO3)2", 1, SuppressError=.TRUE.)
		MgCl2i = FindChem("MgCl2", 1, SuppressError=.TRUE.)
		CaNO32i = FindChem("Ca(NO3)2", 1, SuppressError=.TRUE.)
		CaCl2i = FindChem("CaCl2", 1, SuppressError=.TRUE.)

		NH3i = FindChem("NH3", 1)
		NH4i = FindChem("NH4+", 1)
		NH4NO3i = FindChem("NH4NO3", 1)
		NH4Cli = FindChem("NH4Cl", 1)
		NH4HSO4i = FindChem("NH4HSO4", 1)
		NH42SO4i = FindChem("(NH4)2SO4", 1)

		Cli = FindChem("Cl-", 1)
		NO3i = FindChem("NO3-", 1)
		SO4i = FindChem("SO4--", 1)
		HSO4i = FindChem("HSO4-", 1)
		Hi = FindChem("H+", 1)
		OHi = FindChem("OH-", 1)

		
		!Check if any Cl- present (MJA, 071607)
		Clsum = Y(3) !molecules/cm3
		CurrentParticle => Particles(XP,YP,ZP)%First
		DO I = 1, HowManyLinks
			Clsum = Clsum + (Y(4 + 2*HowManyLinks +I) &
							 + CurrentParticle%AqChems(NH4Cli)*Avogadro &
							 + CurrentParticle%AqChems(KCli)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + I)
			IF (NaCli .NE. 0) THEN
				Clsum = Clsum + (CurrentParticle%AqChems(NaCli)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + I)
			END IF
			IF (MgCl2i .NE. 0) THEN
				Clsum = Clsum + (2.*CurrentParticle%AqChems(MgCl2i)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + I)
			END IF
			IF (CaCl2i .NE. 0) THEN
				Clsum = Clsum + (2.*CurrentParticle%AqChems(CaCl2i)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + I)
			END IF		
		
			CurrentParticle		  => CurrentParticle%Next
		END DO
		IF(Clsum .LE. 1.0) THEN !less than 1 molecule/cm3 of Cl
								  !protects against numerical errors
								  !in advection
			SkipCl = .TRUE.	
			!WRITE(*,*) "SkipCl True.", Clsum
		ELSE
			SkipCl = .FALSE.
			!WRITE(*,*) "SkipCl False.", Clsum, Y(3), Y(4 + 2*HowManyLinks +1)*Y(7 + 11*HowManyLinks + 1)
		END IF

		!Check if any NO3- present (MJA, 073007)
		NO3sum = Y(2) !molecules/cm3
		CurrentParticle => Particles(XP,YP,ZP)%First
		DO I = 1, HowManyLinks
			NO3sum = NO3sum + (Y(4 + HowManyLinks +I) &
							 + CurrentParticle%AqChems(NH4NO3i)*Avogadro &
							 + CurrentParticle%AqChems(KNO3i)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + I)
			IF (NaNO3i .NE. 0) THEN
				NO3sum = NO3sum + (CurrentParticle%AqChems(NaNO3i)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + I)
			END IF
			IF (MgNO32i .NE. 0) THEN
				NO3sum = NO3sum + (2.*CurrentParticle%AqChems(MgNO32i)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + I)
			END IF
			IF (CaNO32i .NE. 0) THEN
				NO3sum = NO3sum + (2.*CurrentParticle%AqChems(CaNO32i)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + I)
			END IF		
	
			CurrentParticle		  => CurrentParticle%Next
		END DO
		
		IF(NO3sum .LE. 1.0) THEN !less than 1 molecule/cm3 of Cl
								  !protects against numerical errors
								  !in advection
			SkipNO3 = .TRUE.	
			!WRITE(*,*) "SkipNO3 True.", NO3sum
		ELSE
			SkipNO3 = .FALSE.
			!WRITE(*,*) "SkipNO3 False.", NO3sum, Y(2), Y(4 + HowManyLinks +1)*Y(7 + 11*HowManyLinks + 1)
		END IF
		
		EqFlag = .FALSE.
		CurrentParticle => Particles(XP,YP,ZP)%First
		DO I = 1, HowManyLinks
						
			!Corrected Henry's Law Constants
			HNO3SoverH(I) = Y(7 + 7*HowManyLinks + I)
			IF(SkipNO3) HNO3SoverH(I) = 1.0
			
			HClSoverH(I) = Y(7 + 8*HowManyLinks + I)
			IF(SkipCl) HClSoverH(I) = 1.0
			
			NH3SoverH(I) = Y(7 + 9*HowManyLinks + I)
			NH3SoverH2(I) = Y(7 + 10*HowManyLinks + I)
			
			NumParticlesPerBin(I) = Y(7 + 11*HowManyLinks + I)

			H2SO4Coeff(I) = Y(7 + 12*HowManyLinks + I)
			HNO3Coeff(I) = Y(7 + 13*HowManyLinks + I)
			IF(SkipNO3) HNO3Coeff(I) = 0.0

			HClCoeff(I) = Y(7 + 14*HowManyLinks + I)
			IF(SkipCl) HClCoeff(I) = 0.0
			
			NH3Coeff(I) = Y(7 + 15*HowManyLinks + I)
						
			!WRITE(*,*) "Try 1"
			DO J = 1, 4+ 6*HowManyLinks
				!WRITE(*,*) J, Y(J)
			END DO 
			
			!Maybe eq will fix problem? (Matt, 10/23/06)
			CurrentParticle%AqChems(SO4i) = Y(4 + I)/Avogadro
			CurrentParticle%AqChems(NO3i) = Y(4 + HowManyLinks +I)/Avogadro
			CurrentParticle%AqChems(Cli) = Y(4 + 2*HowManyLinks +I)/Avogadro
			CurrentParticle%AqChems(NH4i) = Y(4 + 3*HowManyLinks +I)/Avogadro
			CurrentParticle%AqChems(NH3i) = Y(4 + 4*HowManyLinks +I)/Avogadro
			CurrentParticle%AqChems(Hi) = Y(4 + 5*HowManyLinks +I)/Avogadro
							
			IF (Y(4 + 5*HowManyLinks +I) .LT. 0.) THEN
				CALL TRANSCRIPT("H+ below zero in BigCondODE")
				!CALL TRANSCRIPT("Temperature:  "//trim(real2str(Temperature))//"  ")
				!CALL TRANSCRIPT("RH:  "//trim(real2str(RH))//"  ")
				!CALL TRANSCRIPT("Y(H+):  "//trim(real2str(Y(4 + 5*HowManyLinks +I)))//"  ")
				
				Store = Y(4 + 5*HowManyLinks +I)
				!Use NH3(g) + H+(aq) => NH4+(aq) equilibrium to estimate H+ conc.
				Dum1 = Y(4 + 3*HowManyLinks +I)*Y(7 + 9*HowManyLinks + I)/Y(4)
				!use H2O(aq) <=> H+ & OH-
				Dum2 = Avogadro*EqCoeff (1, CurrentParticle, 1) &
												*CurrentParticle%WaterActivity*CurrentParticle%AqChems(1) &
					                           *WaterMolecMass/(CurrentParticle%GammaMixed(1)**2) &
											   /Molality(FindChem("OH-",1), CurrentParticle) &
											   /1000.
				Y(4 + 5*HowManyLinks +I) = MAX(Dum1,Dum2)
				IF(Dum1 .GT. Dum2) THEN
					CALL TRANSCRIPT("Use NH4+ Eq.")
				ELSE
					CALL TRANSCRIPT("Use H2O Eq")
				END IF
			
				CALL TRANSCRIPT(" New Y(H+):  "//trim(real2str(Y(4 + 5*HowManyLinks +I)))//"  ")
				Change = Store - Y(4 + 5*HowManyLinks +I) !Change is negative
				!Writeflag = .TRUE.
				
				CALL TRANSCRIPT("Just fix.")
				CurrentParticle%AqChems(Hi) = CurrentParticle%AqChems(Hi) - Change/Avogadro
			END IF

			IF (Y(4 + HowManyLinks +I) .LT. 0. .AND. NOT(SkipNO3)) THEN
				CALL TRANSCRIPT("NO3- below zero in BigCondODE")
				!WRITE(*,*) Y(4 + HowManyLinks +I), I, Y(7 + 11*HowManyLinks +I)
				!WRITE(*,*) "pH: ", AerosolPH(CurrentParticle)
				!WRITE(*,*) "RH: ", RH 
				Writeflag = .TRUE.
				Store = Y(4 + HowManyLinks +I)
				Y(4 + HowManyLinks +I) = Y(2)/Y(4 + 5*HowManyLinks +I)/Y(7 + 7*HowManyLinks + I)
				Change = Store - Y(4 + HowManyLinks +I) !Change is negative!
				CALL TRANSCRIPT(" New Y(NO3-):  "//trim(real2str(Y(4 + HowManyLinks +I)))//"  ")
				!WRITE(*,*) "Change: ", Change/Avogadro
				
				!CALL TRANSCRIPT("Just fix.")
				CurrentParticle%AqChems(NO3i) = CurrentParticle%AqChems(NO3i) - Change/Avogadro
				
			END IF
			
			IF (Y(4 + 2*HowManyLinks +I) .LT. 0. .AND. NOT(SkipCl)) THEN
				CALL TRANSCRIPT("Cl- below zero in BigCondODE")
				!WRITE(*,*) Y(4 + 2*HowManyLinks +I), I 
				!Writeflag = .TRUE.
				Store = Y(4 + 2*HowManyLinks +I)
				Y(4 + 2*HowManyLinks +I) = Y(3)/Y(4 + 5*HowManyLinks +I)/Y(7 + 8*HowManyLinks + I)
				Change = Store - Y(4 + 2*HowManyLinks +I) !Change is negative!
				!WRITE(*,*) "Change: ", Change
				CALL TRANSCRIPT(" New Y(Cl-):  "//trim(real2str(Y(4 + 2*HowManyLinks +I)))//"  ")
				
				CALL TRANSCRIPT("Just fix.")
				CurrentParticle%AqChems(Cli) = CurrentParticle%AqChems(Cli) - Change/Avogadro
				
			END IF

			IF (Y(4 + 3*HowManyLinks +I) .LT. 0.) THEN
				CALL TRANSCRIPT("NH4+ below zero in BigCondODE")
				WRITE(*,*) Y(4 + 3*HowManyLinks +I), I 
				WRITE(*,*) "RH: ", RH 
				Writeflag = .TRUE.
				Store = Y(4 + 3*HowManyLinks +I)
				Y(4 + 3*HowManyLinks +I) = Y(4)*Y(4 + 5*HowManyLinks +I)/Y(7 + 9*HowManyLinks + I)
				Change = Store - Y(4 + 3*HowManyLinks +I) !Change is negative!
				CALL TRANSCRIPT(" New Y(NH4+):  "//trim(real2str(Y(4 + 2*HowManyLinks +I)))//"  ")

				CALL TRANSCRIPT("Just fix.")
				CurrentParticle%AqChems(NH4i) = CurrentParticle%AqChems(NH4i) - Change/Avogadro
			END IF
		
			Y(4 + I) = CurrentParticle%AqChems(SO4i)*Avogadro
			Y(4 + HowManyLinks +I) = CurrentParticle%AqChems(NO3i)*Avogadro
			Y(4 + 2*HowManyLinks +I) = CurrentParticle%AqChems(Cli)*Avogadro
			Y(4 + 3*HowManyLinks +I) = CurrentParticle%AqChems(NH4i)*Avogadro
			Y(4 + 4*HowManyLinks +I) = CurrentParticle%AqChems(NH3i)*Avogadro
			Y(4 + 5*HowManyLinks +I) = CurrentParticle%AqChems(Hi)*Avogadro
				
			!WRITE(*,*) "Try 2"
			DO J = 1, 4+ 6*HowManyLinks
				!WRITE(*,*) J, Y(J)
			END DO 
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! This is where we actually calculate the   !!
			!! new gas sat concs.                        !!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			 CALL GasSatConcAll(Y, I, HowManyLinks, CurrentParticle, HClRxnIndex, HNO3RxnIndex, &
			                    NH3RxnIndex1, Temperature, &
						        SkipCl, SkipNO3, &
                                HClSoverH(I), HNO3SoverH(I), NH3SoverH(I), HClSat(I), HClFlag(I), &
								HNO3Sat(I), HNO3Flag(I), NH3Sat(I), NH3Flag(I), &
								HClJac(I), HNO3Jac(I), NH3Jac(I))
			  !Check to make sure saturation concentration
			  !is not bigger than total amount
			  IF(HClSat(I) .GT. ClSum) THEN
				!WRITE(*,*) "HClSat is total HCl!"	
				HClSat(I) = ClSum
			  END IF
			  IF(HNO3Sat(I) .GT. NO3Sum) THEN
				!WRITE(*,*) "HNO3Sat is total HNO3!"	
				HNO3Sat(I) = NO3Sum
			  END IF
			  IF(WriteFlag) CALL TRANSCRIPT("HNO3Flag:  "//trim(int2str(HNO3Flag(I)))//"  ")
			  IF(WriteFlag) CALL TRANSCRIPT("HClFlag:  "//trim(int2str(HClFlag(I)))//"  ")
			  IF(WriteFlag) CALL TRANSCRIPT("NH3Flag:  "//trim(int2str(NH3Flag(I)))//"  ")
			
			!WRITE(*,*) "Cl-: ", Y(4 + 2*HowManyLinks +I)
			IF(SkipCl) THEN
				ClNearZero(I) = .TRUE.
				!WRITE(*,*) "Cl- near 0!"
			ELSE
				ClNearZero(I) = .FALSE.
			END IF
			
			!WRITE(*,*) "NO3-: ", Y(4 + HowManyLinks +I)
			IF(SkipNO3) THEN
				NO3NearZero(I) = .TRUE.
                WRITE(*,*) "NO3- near 0!"
			ELSE
				NO3NearZero(I) = .FALSE.
			END IF
						
			CurrentParticle		  => CurrentParticle%Next
		
		END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! The gas phase concentration changes as !!
		!! the negative of the sum of the other	  !!
		!! rates. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		YDOT(1) = 0.
		YDOT(2) = 0.
		YDOT(3) = 0.
		YDOT(4) = 0.

		!Initialize dH+/dt to 0 for summation
		DO I = 1, HowManyLinks
			YDOT(4 + 5*HowManyLinks + I) = 0.
		END DO

		!H2SO4 Condensation
		DO I = 1, HowManyLinks
			
			!SO4-- Concentrations
			YDOT(4 + I) = H2SO4Coeff(I)*Y(1) ! The Gas Phase Conc
						   
			!H+ Concentrations
			YDOT(4 + 5*HowManyLinks + I) = YDOT(4 + 5*HowManyLinks + I) &
										+ 2.*YDOT(4 + I)

			!H2SO4(g) Concentrations
			YDOT(1) = YDOT(1) - YDOT(4 + I) * NumParticlesPerBin(I)
			
		END DO

		!HNO3 Condensation
		!Shut down HNO3 condensation if RH is less than 20% (MJA, 07/12/07)
		IF(RH .GE. 0.20 .OR. METASTABLE) THEN
		DO I = 1, HowManyLinks

			!Skip for a given particle if NO3 is 
			!close to 0 and we are evaporating
			IF(NOT(NO3NearZero(I))) THEN

			!WRITE(*,*) "HNO3 ", Y(2), HNO3Sat(I), HNO3Flag(I)
			!NO3- Concentrations
			YDOT(4 + HowManyLinks + I) = HNO3Coeff(I)*(Y(2) & ! The Gas Phase Conc
											- HNO3Sat(I) )		   
			!H+ Concentrations
			YDOT(4 + 5*HowManyLinks + I) = YDOT(4 + 5*HowManyLinks + I) &
										+ YDOT(4 + HowManyLinks + I)
			!HNO3(g) Concentrations
			YDOT(2) = YDOT(2) - YDOT(4 + HowManyLinks + I) &
								* NumParticlesPerBin(I)
		
			END IF
		END DO
		END IF
				
		!HCl Condensation
		!Shut down HCl condensation if RH is less than 20% (MJA, 07/12/07)
		IF(RH .GE. 0.20 .OR. METASTABLE) THEN
		DO I = 1, HowManyLinks
			!Skip for a given particle if Cl is 
			!close to 0 and we are evaporating
			IF(NOT(ClNearZero(I))) THEN
			
			!Cl- Concentrations
			YDOT(4 + 2*HowManyLinks + I) = HClCoeff(I)*(Y(3) & ! The Gas Phase Conc
											- HClSat(I) )
			!WRITE(*,*) "dCl-/dt :", YDOT(4 + 2*HowManyLinks + I), HClCoeff
			
			!H+ Concentrations
			YDOT(4 + 5*HowManyLinks + I) = YDOT(4 + 5*HowManyLinks + I) &
										+ YDOT(4 + 2*HowManyLinks + I)
			!HCl(g) Concentrations
			YDOT(3) = YDOT(3) - YDOT(4 + 2*HowManyLinks + I) &
								* NumParticlesPerBin(I)
			
			END IF
		END DO
		END IF
		
		!NH3 Condensation
		DO I = 1, HowManyLinks

			!WRITE(*,*) "NH3 ", Y(4), NH3Sat(I), NH3Flag(I)
			!NH4+ Concentrations
			YDOT(4 + 3*HowManyLinks + I) = NH3Coeff(I)*(Y(4) & ! The Gas Phase Conc
											- NH3Sat(I) )
													   
			!H+ Concentrations (Note minus sign!)
			YDOT(4 + 5*HowManyLinks + I) = YDOT(4 + 5*HowManyLinks + I) &
										- YDOT(4 + 3*HowManyLinks + I)
			
			!NH3(aq) Concentrations
			YDOT(4 + 4*HowManyLinks + I) = NH3Coeff(I)*(Y(4) & ! The Gas Phase Conc
											- NH3SoverH2(I)*Y(4 + 4*HowManyLinks + I))
			
			!NH3(g) Concentrations
			YDOT(4) = YDOT(4) -  (YDOT(4 + 3*HowManyLinks + I) &
								+ YDOT(4 + 4*HowManyLinks + I)) &
								* NumParticlesPerBin(I)
		END DO

		DEALLOCATE(H2SO4Coeff, HNO3Coeff, HClCoeff, NH3Coeff)
		DEALLOCATE(HNO3Sat, HClSat, NH3Sat)
		DEALLOCATE(HNO3Flag, HClFlag, NH3Flag)
		DEALLOCATE(HNO3SoverH, HClSoverH, NH3SoverH, NH3SoverH2)		

	RETURN
END SUBROUTINE BigCondensationODEEvaluator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate a JACOBIAN for Condensation !!
!!                                       !!
!! As per LSODES's specification, only   !!
!! calculate at the stated row I         !!
!! Rewritten to do all inroagncis together,
!! Matt Alvarado 10/13/06				 !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BigCondensationJacobianEvaluator (neq, t, y, I, ia, ja, pdj)

	USE GridPointFields,     ONLY : GasPhaseChemicalRates
	USE InfrastructuralCode, ONLY : ERROR, Transcript
	USE Aerosols,            ONLY : Particles, Particle, ParticleArray
	USE GridPointFields,	 ONLY : LocateGridPoint,				&
										getM,							&
										GasPhaseChemGridIndex,			&
										GetTemp, GetRelativeHumidity
	USE ModelParameters
	USE Chemistry,			 ONLY : FindChem, HowManyAqEqReactions
	IMPLICIT NONE

	!! External Variables
	REAL*8		:: t, y, pdj
	INTEGER		:: I, neq, ia, ja
	DIMENSION	:: y(1), ia(1), ja(1), pdj(1)

	!! Internal Variables
	INTEGER :: J, Q, C, HowManyLinks, XP, YP, ZP, allocation_error
	INTEGER :: H2SO4RxnIndex, HNO3RxnIndex, HClRxnIndex
	INTEGER ::  NH3RxnIndex1, NH3RxnIndex2, CaseFlag
	Integer :: 	KNO3i, KCli, KHSO4i, &
				Nai, NaNO3i, NaCli, NaHSO4i, Na2SO4i, &
				NH3i, NH4i, NH4NO3i, NH4Cli, NH4HSO4i, NH42SO4i, &
				Cli, NO3i, SO4i, HSO4i, K2SO4i, Ki, Hi, &
				MgNO32i, MgCl2i, CaCl2i, CaNO32i
	

	REAL*8				:: Temperature, XX(3), Change, Store, RH
	REAL*8, ALLOCATABLE :: H2SO4Coeff(:) ! 1/s
	REAL*8, ALLOCATABLE :: HNO3Coeff(:) ! 1/s
	REAL*8, ALLOCATABLE :: HClCoeff(:) ! 1/s
	REAL*8, ALLOCATABLE :: NH3Coeff(:) ! 1/s
	REAL*8, ALLOCATABLE :: NumParticlesPerBin(:) ! particles/cm3
	REAL*8, ALLOCATABLE :: HNO3Sat(:), HClSat(:), NH3Sat(:) !molecules/cm3
	INTEGER, ALLOCATABLE :: HNO3Flag(:), HClFlag(:), NH3Flag(:) 
	REAL*8, ALLOCATABLE :: HNO3SoverH(:), HClSoverH(:)
	REAL*8, ALLOCATABLE :: NH3SoverH(:), NH3SoverH2(:) 
	REAL*8, ALLOCATABLE :: HClJac(:), HNO3Jac(:), NH3Jac(:)
	LOGICAL, ALLOCATABLE :: ClNearZero(:), NO3NearZero(:)

	TYPE(Particle), POINTER :: CurrentParticle

	LOGICAL				:: EqFlag, SkipCl, SkipNO3
	LOGICAL, PARAMETER  :: Scaffolding =  .FALSE.
	REAL*8 :: ClSum, NO3Sum
	REAL*8, PARAMETER :: Tiny = 1.0e-10 !molecules/particle 

	IF (SCAFFOLDING) print *, ">>> In BigCondensationJacobianEvaluator <<<"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICALS										!!
	!! The Y-Vector we will integrate with and push through LSODES looks like this:			!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!    INTEGRATED TERMS (1 to 4 + 6*HowManyLinks)										!!
	!!    (1)								            : H2SO4 Gas Phase Concentration		!!
	!!    (2)								            : HNO3 Gas Phase Concentration		!!
	!!    (3)								            : HCl Gas Phase Concentration		!!
	!!    (4)								            : NH3 Gas Phase Concentration		!!	
	!!    (5 to 4 + HowManyLinks)				        : SO4-- Aerosol Concentrations		!!
	!!    (5 + HowManyLinks to 4 + 2*HowManyLinks)      : NO3- Aerosol Concentrations		!!
	!!    (5 + 2*HowManyLinks to 4 + 3*HowManyLinks)    : Cl- Aerosol Concentrations		!!
	!!    (5 + 3*HowManyLinks to 4 + 4*HowManyLinks)    : NH4+ Aerosol Concentrations		!!
	!!    (5 + 4*HowManyLinks to 4 + 5*HowManyLinks)    : NH3(aq) Aerosol Concentrations	!!
	!!    (5 + 5*HowManyLinks to 4 + 6*HowManyLinks)    : H+ Aerosol Concentrations			!!
	!!																						!!
	!!	  PARAMETER TERMS																	!!
	!!	  (5 + 6*HowManyLinks to ??)														!!
	!!																						!!
	!!    (5 + 6*HowManyLinks to 7 + 6*HowManyLinks)    : Location (XP, YP, ZP)				!!
	!!    (8 + 7*HowManyLinks to 7 + 8*HowManyLinks)    : (S/H) for HNO3(g) <=> H+ & NO3-	!!
	!!    (8 + 8*HowManyLinks to 7 + 9*HowManyLinks)    : (S/H) for HCl(g) <=> H+ & Cl-		!!
	!!    (8 + 9*HowManyLinks to 7 + 10*HowManyLinks)   : (S/H) for NH3(g) & H+(aq) <=> NH4+!!
	!!    (8 + 10*HowManyLinks to 7 + 11*HowManyLinks)  : (S/H) for NH3(g) <=> NH3 (aq)		!!
	!!	  (8 + 11*HowManyLinks to 7 + 12*HowManyLinks)  : NumberofParticles per bin			!!
	!!	  (8 + 12*HowManyLinks to 7 + 13*HowManyLinks)  : H2SO4 Rate Coefficient			!!
	!!	  (8 + 13*HowManyLinks to 7 + 14*HowManyLinks)  : HNO3 Rate Coefficient				!!
	!!	  (8 + 14*HowManyLinks to 7 + 15*HowManyLinks)  : HCl Rate Coefficient				!!
	!!	  (8 + 15*HowManyLinks to 7 + 16*HowManyLinks)  : NH3 Rate Coefficient				!!
	!!    (8 + 16*HowManyLinks)							: H2SO4 Rxn Index					!!
	!!    (9 + 16*HowManyLinks)							: HNO3 Rxn Index					!!
	!!    (10 + 16*HowManyLinks)						: HCl Rxn Index						!!
	!!    (11 + 16*HowManyLinks)						: NH3(g) & H+ <=> NH4+ Rxn Index	!!
	!!    (12 + 16*HowManyLinks)						: NH3(g) <=> NH3(aq) Rxn Index		!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! HowManyLinks is the number of (non-zero) sections in the size dist					!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	HowManyLinks = (neq - 4)/6
				
	!Allocate parameter arrays
	ALLOCATE (H2SO4Coeff(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of H2SO4RateCoeff could not proceed in CondensationODEEvaluator")
	ALLOCATE (HNO3Coeff(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3RateCoeff could not proceed in CondensationODEEvaluator")
	ALLOCATE (HClCoeff(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HClRateCoeff could not proceed in CondensationODEEvaluator")
	ALLOCATE (NH3Coeff(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of NH3RateCoeff could not proceed in CondensationODEEvaluator")
	
	ALLOCATE (NumParticlesPerBin(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of NumParticlesPerBin could not proceed in CondensationODEEvaluator")
	
	ALLOCATE (HNO3Sat(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3Sat could not proceed in CondensationODEEvaluator")
	ALLOCATE (HClSat(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HClSat could not proceed in CondensationODEEvaluator")
	ALLOCATE (NH3Sat(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of NH3Sat could not proceed in CondensationODEEvaluator")
		
	ALLOCATE (HNO3Flag(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3Flag could not proceed in CondensationODEEvaluator")
	ALLOCATE (HClFlag(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HClFlag could not proceed in CondensationODEEvaluator")
	ALLOCATE (NH3Flag(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of NH3Flag could not proceed in CondensationODEEvaluator")
		
	ALLOCATE (HNO3SoverH(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3SoverH could not proceed in CondensationODEEvaluator")
	ALLOCATE (HClSoverH(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HClSoverH could not proceed in CondensationODEEvaluator")
	ALLOCATE (NH3SoverH(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of NH3SoverH could not proceed in CondensationODEEvaluator")
	ALLOCATE (NH3SoverH2(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of NH3SoverH2 could not proceed in CondensationODEEvaluator")
		
	ALLOCATE (ClNearZero(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of ClNearZero could not proceed in CondensationODEEvaluator")
	ALLOCATE (NO3NearZero(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of NO3NearZero could not proceed in CondensationODEEvaluator")
	
	ALLOCATE (HNO3Jac(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HNO3Sat could not proceed in CondensationODEEvaluator")
	ALLOCATE (HClJac(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of HClSat could not proceed in CondensationODEEvaluator")
	ALLOCATE (NH3Jac(HowManyLinks), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of NH3Sat could not proceed in CondensationODEEvaluator")

	
	XP = INT(Y(5 + 6*HowManyLinks))
	YP = INT(Y(6 + 6*HowManyLinks))
	ZP = INT(Y(7 + 6*HowManyLinks))

	XX          = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
	Temperature = GetTemp(XX(1),XX(2),XX(3))
	RH =GetRelativeHumidity(XX(1),XX(2),XX(3)) 
	
	H2SO4RxnIndex = Y(8 + 16*HowManyLinks)
	HNO3RxnIndex = Y(9 + 16*HowManyLinks)
	HClRxnIndex = Y(10 + 16*HowManyLinks)
	NH3RxnIndex1 = Y(11 + 16*HowManyLinks)
	NH3RxnIndex2 = Y(12 + 16*HowManyLinks)

	!Find chemical indices for later
	Ki = FindChem("K+", 1)
	KNO3i = FindChem("KNO3", 1)
	KCli = FindChem("KCl", 1)
	KHSO4i = FindChem("KHSO4", 1)
	K2SO4i = FindChem("K2SO4", 1)
		
	Nai = FindChem("Na+", 1, SuppressError=.TRUE.)
	NaNO3i = FindChem("NaNO3", 1, SuppressError=.TRUE.)
	NaCli = FindChem("NaCl", 1, SuppressError=.TRUE.)
	NaHSO4i = FindChem("NaHSO4", 1, SuppressError=.TRUE.)
	Na2SO4i = FindChem("Na2SO4", 1, SuppressError=.TRUE.)
	
	MgNO32i = FindChem("Mg(NO3)2", 1, SuppressError=.TRUE.)
	MgCl2i = FindChem("MgCl2", 1, SuppressError=.TRUE.)
	CaNO32i = FindChem("Ca(NO3)2", 1, SuppressError=.TRUE.)
	CaCl2i = FindChem("CaCl2", 1, SuppressError=.TRUE.)
		
	NH3i = FindChem("NH3", 1)
	NH4i = FindChem("NH4+", 1)
	NH4NO3i = FindChem("NH4NO3", 1)
	NH4Cli = FindChem("NH4Cl", 1)
	NH4HSO4i = FindChem("NH4HSO4", 1)
	NH42SO4i = FindChem("(NH4)2SO4", 1)

	Cli = FindChem("Cl-", 1)
	NO3i = FindChem("NO3-", 1)
	SO4i = FindChem("SO4--", 1)
	HSO4i = FindChem("HSO4-", 1)
	Hi = FindChem("H+", 1)

	!Check if any Cl- present (MJA, 071607)
	Clsum = Y(3) !molecules/cm3
	CurrentParticle => Particles(XP,YP,ZP)%First
	DO Q = 1, HowManyLinks
		Clsum = Clsum + (Y(4 + 2*HowManyLinks +Q) &
						 + CurrentParticle%AqChems(NH4Cli)*Avogadro &
						 + CurrentParticle%AqChems(KCli)*Avogadro &
						 ) *Y(7 + 11*HowManyLinks + Q)
	
		IF (NaCli .NE. 0) THEN
			Clsum = Clsum + (CurrentParticle%AqChems(NaCli)*Avogadro &
							 ) *Y(7 + 11*HowManyLinks + Q)
		END IF
		IF (MgCl2i .NE. 0) THEN
			Clsum = Clsum + (2.*CurrentParticle%AqChems(MgCl2i)*Avogadro &
						 ) *Y(7 + 11*HowManyLinks + Q)
		END IF
		IF (CaCl2i .NE. 0) THEN
			Clsum = Clsum + (2.*CurrentParticle%AqChems(CaCl2i)*Avogadro &
						 ) *Y(7 + 11*HowManyLinks + Q)
		END IF		
	
		CurrentParticle		  => CurrentParticle%Next
	END DO
	IF(Clsum .LE. 1.0) THEN !less than 1 molecule/cm3 of Cl
							  !protects against numerical errors
							  !in advection
		SkipCl = .TRUE.	
		!WRITE(*,*) "SkipCl True.", Clsum
	ELSE
		SkipCl = .FALSE.
		!WRITE(*,*) "SkipCl False.", Clsum, Y(3), Y(4 + 2*HowManyLinks +1)*Y(7 + 11*HowManyLinks + 1)
	END IF

	!Check if any NO3- present (MJA, 073007)
	NO3sum = Y(2) !molecules/cm3
	CurrentParticle => Particles(XP,YP,ZP)%First
	DO Q = 1, HowManyLinks
		NO3sum = NO3sum + (Y(4 + HowManyLinks +Q) &
						 + CurrentParticle%AqChems(NH4NO3i)*Avogadro &
						 + CurrentParticle%AqChems(KNO3i)*Avogadro &
						 ) *Y(7 + 11*HowManyLinks + Q)
		IF (NaNO3i .NE. 0) THEN
			NO3sum = NO3sum + (CurrentParticle%AqChems(NaNO3i)*Avogadro &
						 ) *Y(7 + 11*HowManyLinks + Q)
		END IF
		IF (MgNO32i .NE. 0) THEN
			NO3sum = NO3sum + (CurrentParticle%AqChems(MgNO32i)*Avogadro &
						 ) *Y(7 + 11*HowManyLinks + Q)
		END IF
		IF (CaNO32i .NE. 0) THEN
			NO3sum = NO3sum + (CurrentParticle%AqChems(CaNO32i)*Avogadro &
						 ) *Y(7 + 11*HowManyLinks + Q)
		END IF	
		
		
		CurrentParticle		  => CurrentParticle%Next
	END DO
	
	IF(NO3sum .LE. 1.0) THEN !less than 1 molecule/cm3 of Cl
							  !protects against numerical errors
							  !in advection
		SkipNO3 = .TRUE.	
		!WRITE(*,*) "SkipNO3 True.", NO3sum
	ELSE
		SkipNO3 = .FALSE.
		!WRITE(*,*) "SkipNO3 False.", NO3sum, Y(2), Y(4 + HowManyLinks +1)*Y(7 + 11*HowManyLinks + 1)
	END IF
	
	CurrentParticle => Particles(XP,YP,ZP)%First
	DO Q = 1, HowManyLinks
			
		!Corrected Henry's Law Constants
		HNO3SoverH(Q) = Y(7 + 7*HowManyLinks + Q)
		IF(SkipNO3) HNO3SoverH(Q) = 1.0

		HClSoverH(Q) = Y(7 + 8*HowManyLinks + Q)
		IF(SkipCl) HClSoverH(Q) = 1.0
		
		NH3SoverH(Q) = Y(7 + 9*HowManyLinks + Q)
		NH3SoverH2(Q) = Y(7 + 10*HowManyLinks + Q)
			
		NumParticlesPerBin(Q) = Y(7 + 11*HowManyLinks + Q)

		H2SO4Coeff(Q) = Y(7 + 12*HowManyLinks + Q)
		
		HNO3Coeff(Q) = Y(7 + 13*HowManyLinks + Q)
		IF(SkipNO3) HNO3Coeff(Q) = 0.0
		
		HClCoeff(Q) = Y(7 + 14*HowManyLinks + Q)
		IF(SkipCl) HClCoeff(Q) = 0.0
		
		NH3Coeff(Q) = Y(7 + 15*HowManyLinks + Q)

		CALL GasSatConcAll(Y, Q, HowManyLinks, CurrentParticle, HClRxnIndex, HNO3RxnIndex, &
			               NH3RxnIndex1, Temperature, &
						    SkipCl, SkipNO3, &
                            HClSoverH(Q), HNO3SoverH(Q), NH3SoverH(Q), HClSat(Q), HClFlag(Q), &
						    HNO3Sat(Q), HNO3Flag(Q), NH3Sat(Q), NH3Flag(Q), &
								HClJac(Q), HNO3Jac(Q), NH3Jac(Q))
		!Check to make sure saturation concentration
		!is not bigger than total amount
		IF(HClSat(Q) .GT. ClSum) THEN
			!WRITE(*,*) "HClSat is total HCl!"	
			HClSat(Q) = ClSum
		END IF
		IF(HNO3Sat(Q) .GT. NO3Sum) THEN
			!WRITE(*,*) "HNO3Sat is total HNO3!"	
			HNO3Sat(Q) = NO3Sum
		END IF

		IF(SkipCl) THEN
			ClNearZero(Q) = .TRUE.
		ELSE
			ClNearZero(Q) = .FALSE.
		END IF
		IF(SkipNO3) THEN
			NO3NearZero(Q) = .TRUE.
		ELSE
			NO3NearZero(Q) = .FALSE.
		END IF
		
		CurrentParticle		  => CurrentParticle%Next
	END DO


	IF (I .LE. 4) THEN !Gases
		CaseFlag = I
	ELSEIF (I .GT. 4 .AND. I .LE. 4 + HowManyLinks) THEN !SO4--
		CaseFlag = 5
	ELSE IF (I .GT. 4 + HowManyLinks .AND. I .LE. 4 + 2*HowManyLinks) THEN
		CaseFlag = 6	 
	ELSE IF (I .GT. 4 + 2*HowManyLinks .AND. I .LE. 4 + 3*HowManyLinks) THEN
		CaseFlag = 7	 
	ELSE IF (I .GT. 4 + 3*HowManyLinks .AND. I .LE. 4 + 4*HowManyLinks) THEN
		CaseFlag = 8	 
	ELSE IF (I .GT. 4 + 4*HowManyLinks .AND. I .LE. 4 + 5*HowManyLinks) THEN
		CaseFlag = 9	 
	ELSE IF (I .GT. 4 + 5*HowManyLinks .AND. I .LE. 4 + 6*HowManyLinks) THEN
		CaseFlag = 10	 
	END IF
	
	!! REMEMBER :: I is the index of the d/d(conc) (column), J the index of d(conc)/dt (row)
	SELECT CASE (CaseFlag)	
	
	
	!!I = 1 is for d/d(H2SO4(g))
	CASE(1) 
		
		PDJ(1) = 0.
		DO J = 1, HowManyLinks
			
			!d(H2SO4(g))/dt
			PDJ(1) = PDJ(1) - H2SO4Coeff(J)*NumParticlesPerBin(J)
			
			!d(SO4--)/dt
			PDJ(4 + J) = H2SO4Coeff(J)

			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = 2*H2SO4Coeff(J)
	
		END DO

	!!I = 2 is for d/d(HNO3(g))
	CASE(2) 
		
		PDJ(2) = 0.
		PDJ(3) = 0.
		PDJ(4) = 0.

		DO J = 1, HowManyLinks
			
			!Cut out HNO3 and HCl condensation if RH is .LT. 0.20
			!MJA, 07/12/07
			IF((RH .GE. 0.20 .OR. METASTABLE) .AND. NOT(NO3NearZero(J))) THEN
			!d(HNO3(g))/dt
			PDJ(2) = PDJ(2) - HNO3Coeff(J)*NumParticlesPerBin(J)
			
			!d(NO3-)/dt
			PDJ(4 + HowManyLinks + J) = HNO3Coeff(J)

			!d(Cl-)/dt
			PDJ(4 + 2*HowManyLinks + J) = 0.
			
			!d(NH4+)/dt
			PDJ(4 + 3*HowManyLinks + J) = 0.
			
			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = HNO3Coeff(J)
			
			
			IF (HNO3Flag(J) .EQ. 33 .OR. HNO3Flag(J) .EQ. 34) THEN !NH3-HCl-HNO3 gas-solid
				!In these cases, HNO3Sat is a function of HNO3(g)
				!d(HNO3(g))/dt
				PDJ(2) = PDJ(2) + HNO3Coeff(J)*NumParticlesPerBin(J)*HNO3Jac(J)
			
				!d(NO3-)/dt
				PDJ(4 + HowManyLinks + J) = PDJ(4 + HowManyLinks + J)-HNO3Coeff(J)*HNO3Jac(J)

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)-HNO3Coeff(J)*HNO3Jac(J)

			END IF

			IF (HClFlag(J) .EQ. 33) THEN !NH3-HCl-HNO3 gas-solid
				!In this case, HClSat is a function of HNO3(g)
				!d(HCl(g))/dt
				PDJ(3) = PDJ(3) + HClCoeff(J)*NumParticlesPerBin(J)*HClJac(J)
			
				!d(Cl-)/dt
				PDJ(4 + 2*HowManyLinks + J) = PDJ(4 + 2*HowManyLinks + J)-HClCoeff(J)*HClJac(J)

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)-HClCoeff(J)*HClJac(J)
			END IF

			IF (NH3Flag(J) .EQ. 33. .OR. NH3Flag(J) .EQ. 34.) THEN !NH3-HCl-HNO3 gas-solid
				!In these cases, NH3Sat is a function of HNO3(g)
				!d(NH3(g))/dt
				PDJ(4) = PDJ(4) + NH3Coeff(J)*NumParticlesPerBin(J)*(-1.0*NH3Jac(J))
			
				!d(NH4+)/dt
				PDJ(4 + 3*HowManyLinks + J) = PDJ(4 + 3*HowManyLinks + J)-NH3Coeff(J)*(-1.0*NH3Jac(J))

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)+NH3Coeff(J)*(-1.0*NH3Jac(J))
			END IF


			END IF
		END DO	

	!!I = 3 is for d/d(HCl(g))
	CASE(3) 
		
		PDJ(2) = 0.
		PDJ(3) = 0.
		PDJ(4) = 0.

		DO J = 1, HowManyLinks
			
			!Cut out HNO3 and HCl condensation if RH is .LT. 0.20
			!MJA, 07/12/07
			IF((RH .GE. 0.20 .OR. METASTABLE) .AND. NOT(ClNearZero(J))) THEN
			!d(HCl(g))/dt
			PDJ(3) = PDJ(3) - HClCoeff(J)*NumParticlesPerBin(J)
			
			!d(Cl-)/dt
			PDJ(4 + 2*HowManyLinks + J) = HClCoeff(J)

			!d(NH4+)/dt
			PDJ(4 + 3*HowManyLinks + J) = 0.

			!d(NO3-)/dt
			PDJ(4 + HowManyLinks + J) = 0.
			
			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = HClCoeff(J)
						
			IF (HNO3Flag(J) .EQ. 33) THEN !NH3-HCl-HNO3 gas-solid
				!In these cases, HNO3Sat is a function of HCl(g)
				!d(HNO3(g))/dt
				PDJ(2) = PDJ(2) + HNO3Coeff(J)*NumParticlesPerBin(J)*HNO3Jac(J)
			
				!d(NO3-)/dt
				PDJ(4 + HowManyLinks + J) = PDJ(4 + HowManyLinks + J)-HNO3Coeff(J)*HNO3Jac(J)

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)-HNO3Coeff(J)*HNO3Jac(J)

			END IF

			IF (HClFlag(J) .EQ. 33 .OR. HClFlag(J) .EQ. 35) THEN !NH3-HCl-HNO3 gas-solid
				!In these cases, HClSat is a function of HCl(g)
				!d(HCl(g))/dt
				PDJ(3) = PDJ(3) + HClCoeff(J)*NumParticlesPerBin(J)*HClJac(J)
			
				!d(Cl-)/dt
				PDJ(4 + 2*HowManyLinks + J) = PDJ(4 + 2*HowManyLinks + J)-HClCoeff(J)*HClJac(J)

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)-HClCoeff(J)*HClJac(J)
			END IF

			IF (NH3Flag(J) .EQ. 33. .OR. NH3Flag(J) .EQ. 35.) THEN !NH3-HCl-HNO3 gas-solid
				!In these cases, NH3Sat is a function of HCl(g)
				!d(NH3(g))/dt
				PDJ(4) = PDJ(4) + NH3Coeff(J)*NumParticlesPerBin(J)*(-1.0*NH3Jac(J))
			
				!d(NH4+)/dt
				PDJ(4 + 3*HowManyLinks + J) = PDJ(4 + 3*HowManyLinks + J)-NH3Coeff(J)*(-1.0*NH3Jac(J))

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)+NH3Coeff(J)*(-1.0*NH3Jac(J))
			END IF
	
			END IF
		END DO	

	!!I = 4 is for d/d(NH3(g))
	CASE(4) !!!!!!!!!!!!!!
		
		PDJ(2) = 0.
		PDJ(3) = 0.
		PDJ(4) = 0.
	
		DO J = 1, HowManyLinks
			
			!d(NH3(g))/dt
			PDJ(4) = PDJ(4) - NH3Coeff(J)*NumParticlesPerBin(J)

			!d(NO3-)/dt
			PDJ(4 + HowManyLinks + J) = 0.

			!d(Cl-)/dt
			PDJ(4 + 2*HowManyLinks + J) = 0.

			!d(NH4+)/dt
			PDJ(4 + 3*HowManyLinks + J) = NH3Coeff(J)

			!d(NH3(aq))/dt
			PDJ(4 + 4*HowManyLinks + J) = NH3Coeff(J)

			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = -1.0*NH3Coeff(J)

			
			!Cut out HNO3 and HCl condensation if RH is .LT. 0.20
			!MJA, 07/12/07
			IF(RH .GE. 0.20 .OR. METASTABLE) THEN			
			IF ((HNO3Flag(J) .EQ. 33 .OR. HNO3Flag(J) .EQ. 34) .AND. NOT(NO3NearZero(J))) THEN !NH3-HCl-HNO3 gas-solid
				!In these cases, HNO3Sat is a function of NH3(g)
				!d(HNO3(g))/dt
				PDJ(2) = PDJ(2) + HNO3Coeff(J)*NumParticlesPerBin(J)*(-1.0*HNO3Jac(J))
			
				!d(NO3-)/dt
				PDJ(4 + HowManyLinks + J) = PDJ(4 + HowManyLinks + J)-HNO3Coeff(J)*(-1.0*HNO3Jac(J))

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)-HNO3Coeff(J)*(-1.0*HNO3Jac(J))

			END IF

			IF ((HClFlag(J) .EQ. 33 .OR. HClFlag(J) .EQ. 35) .AND. NOT(ClNearZero(J))) THEN !NH3-HCl-HNO3 gas-solid
				!In these cases, HClSat is a function of NH3(g)
				!d(HCl(g))/dt
				PDJ(3) = PDJ(3) + HClCoeff(J)*NumParticlesPerBin(J)*(-1.0*HClJac(J))
			
				!d(Cl-)/dt
				PDJ(4 + 2*HowManyLinks + J) = PDJ(4 + 2*HowManyLinks + J)-HClCoeff(J)*(1.0*HClJac(J))

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)-HClCoeff(J)*(-1.0*HClJac(J))
			END IF

			IF (NH3Flag(J) .EQ. 33. .OR. NH3Flag(J) .EQ. 34. .OR. NH3Flag(J) .EQ. 35.) THEN !NH3-HCl-HNO3 gas-solid
				!In these cases, NH3Sat is a function of NH3(g)
				!d(NH3(g))/dt
				PDJ(4) = PDJ(4) + NH3Coeff(J)*NumParticlesPerBin(J)*(NH3Jac(J))
			
				!d(NH4+)/dt
				PDJ(4 + 3*HowManyLinks + J) = PDJ(4 + 3*HowManyLinks + J)-NH3Coeff(J)*(NH3Jac(J))

				!d(H+)/dt
				PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J)+NH3Coeff(J)*(NH3Jac(J))
			END IF
			
			END IF
		END DO

	!!For d/d(SO4--)
	CASE(5)

		J = I - 4
		
		PDJ(4) = 0.
		
		!d(NH4+)/dt
		PDJ(4 + 3*HowManyLinks + J) = 0.

		!d(H+)/dt
		PDJ(4 + 5*HowManyLinks + J) = 0.
	
	!!For d/d(NO3-)
	CASE(6) 

		J = I - (4 + HowManyLinks)
		
		PDJ(2) = 0.

		!Cut out HNO3 and HCl condensation if RH is .LT. 0.20
		!MJA, 07/12/07
		IF(RH .GE. 0.20 .OR. METASTABLE) THEN
		!d(NO3-)/dt
		PDJ(4 + HowManyLinks + J) = 0.
			
		!d(H+)/dt
		PDJ(4 + 5*HowManyLinks + J) = 0.

		IF (HNO3Flag(J) .EQ. 10 .AND. NOT(NO3NearZero(J))) THEN

			!d(HNO3(g))/dt
			PDJ(2) = PDJ(2) + NumParticlesPerBin(J) &
						*HNO3Coeff(J)*HNO3Sat(J)/(Y(I)+Tiny)
				
			!d(NO3-)/dt
			PDJ(4 + HowManyLinks + J) = PDJ(4 + HowManyLinks + J) &
										- HNO3Coeff(J)*HNO3Sat(J)/(Y(I)+Tiny)
			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J) &
											- HNO3Coeff(J)*HNO3Sat(J)/(Y(I)+Tiny)
		END IF

		END IF
	!!For d/d(Cl-)
	CASE(7) 

		J = I - (4 + 2*HowManyLinks)
		
		PDJ(3) = 0.

		!Cut out HNO3 and HCl condensation if RH is .LT. 0.20
		!MJA, 07/12/07
		IF(RH .GE. 0.20 .OR. METASTABLE) THEN
		!d(Cl-)/dt
		PDJ(4 + 2*HowManyLinks + J) = 0.
			
		!d(H+)/dt
		PDJ(4 + 5*HowManyLinks + J) = 0.

		IF (HClFlag(J) .EQ. 10 .AND. NOT(ClNearZero(J))) THEN

			!d(HCl(g))/dt
			PDJ(3) = PDJ(3) + NumParticlesPerBin(J) &
						*HClCoeff(J)*HClSat(J)/(Y(I)+Tiny)

			!d(Cl-)/dt
			PDJ(4 + 2*HowManyLinks + J) = PDJ(4 + 2*HowManyLinks + J) &
										- HClCoeff(J)*HClSat(J)/(Y(I)+Tiny)
			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J) &
										- HClCoeff(J)*HClSat(J)/(Y(I)+Tiny)	
	
		
		END IF

		END IF

	!!For d/d(NH4+)
	CASE(8) 

		J = I - (4 + 3*HowManyLinks)

		PDJ(4) = 0.

		!d(NH4+)/dt
		PDJ(4 + 3*HowManyLinks + J) = 0.

		!d(H+)/dt
		PDJ(4 + 5*HowManyLinks + J) = 0.

			
		IF (NH3Flag(J) .EQ. 10) THEN
			
			!d(NH3(g))/dt
			PDJ(4) = PDJ(4) + NumParticlesPerBin(J) &
						*NH3Coeff(J)*NH3Sat(J)/(Y(I)+Tiny)
		
			!d(NH4+)/dt
			PDJ(4 + 3*HowManyLinks + J) = PDJ(4 + 3*HowManyLinks + J) &
										 - NH3Coeff(J)*NH3Sat(J)/(Y(I)+Tiny)

			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J) &
										  + NH3Coeff(J)*NH3Sat(J)/(Y(I)+Tiny)
		
		END IF

	!!For d/d(NH3(aq))
	CASE(9) 

		J = I - (4 + 4*HowManyLinks)
		
		PDJ(4) = NH3Coeff(J)*NH3SoverH2(J)*NumParticlesperBin(J)

		PDJ(4 + 4*HowManyLinks + J) = -1.0*NH3Coeff(J)*NH3SoverH2(J)

	!!For d/d(H+)
	CASE(10) 

		J = I - (4 + 5*HowManyLinks)

		PDJ(2) = 0.
		PDJ(3) = 0.
		PDJ(4) = 0.

		!d(NO3-)/dt
		PDJ(4 + HowManyLinks + J) = 0.
		
		!d(Cl-)/dt
		PDJ(4 + 2*HowManyLinks + J) = 0.
			
		!d(NH4+)/dt
		PDJ(4 + 3*HowManyLinks + J) = 0.

		!d(H+)/dt
		PDJ(4 + 5*HowManyLinks + J) = 0.

		!Cut out HNO3 and HCl condensation if RH is .LT. 0.20
		!MJA, 07/12/07
		IF(RH .GE. 0.20 .OR. METASTABLE) THEN
		IF (HNO3Flag(J) .LT. 30 .AND. NOT(NO3NearZero(J))) THEN
		!Not gas-solid
			
			!d(HNO3(g))/dt
			PDJ(2) = PDJ(2) + NumParticlesPerBin(J) &
						*HNO3Coeff(J)*HNO3Sat(J)/(Y(I)+Tiny)

			!d(NO3-)/dt
			PDJ(4 + HowManyLinks + J) = PDJ(4 + HowManyLinks + J) &
										- HNO3Coeff(J)*HNO3Sat(J)/(Y(I)+Tiny)
		
			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J) &
										  - HNO3Coeff(J)*HNO3Sat(J)/(Y(I)+Tiny)

		END IF


		IF (HClFlag(J) .LT. 30 .AND. NOT(ClNearZero(J))) THEN
			
			!d(HCl(g))/dt
			PDJ(3) = PDJ(3) + NumParticlesPerBin(J) &
						*HClCoeff(J)*HClSat(J)/(Y(I)+Tiny)

			!d(Cl-)/dt
			PDJ(4 + 2*HowManyLinks + J) = PDJ(4 + 2*HowManyLinks + J) &
										- HClCoeff(J)*HClSat(J)/(Y(I)+Tiny)
		
			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J) &
										  - HClCoeff(J)*HClSat(J)/(Y(I)+Tiny)

		END IF
		END IF

		IF (NH3Flag(J) .LT. 30) THEN

			!d(NH3(g))/dt
			PDJ(4) = PDJ(4) - NumParticlesPerBin(J) &
						*NH3Coeff(J)*NH3Sat(J)/(Y(I)+Tiny)

			!d(NH4+)/dt
			PDJ(4 + 3*HowManyLinks + J) = PDJ(4 + 3*HowManyLinks + J) &
										  + NH3Coeff(J)*NH3Sat(J)/(Y(I)+Tiny)

			!d(H+)/dt
			PDJ(4 + 5*HowManyLinks + J) = PDJ(4 + 5*HowManyLinks + J) &
										  - NH3Coeff(J)*NH3Sat(J)/(Y(I)+Tiny)

		END IF

	END SELECT

	DEALLOCATE(H2SO4Coeff, HNO3Coeff, HClCoeff, NH3Coeff)
	DEALLOCATE(HNO3Sat, HClSat, NH3Sat)
	DEALLOCATE(HNO3Flag, HClFlag, NH3Flag)
	DEALLOCATE(HNO3SoverH, HClSoverH, NH3SoverH, NH3SoverH2)		

	RETURN
END SUBROUTINE BigCondensationJacobianEvaluator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DissolutionFactors --									!!
!!															!!
!! Give the dissolution routines everything they need to	!!
!! calculate one iteration of one chemical, including:		!!
!!															!!
!!  1. S' / Effective Henry's Law Ratios					!!
!!  2. Initial Appropriate Aerosol Concentrations			!!
!!															!!
!! ReturnType = 1 if fine, -1 if there is none of			!!
!! the chemical in the gridpoint							!!
SUBROUTINE DissolutionFactors (XP, YP, ZP, WhichRxn, HowManyLinks, Temperature, ReactionType,	&
							   InitialConcs1, InitialConcs2, InitialConcs3, SoverH, SoverH2, ReturnType)


	USE GridPointFields, ONLY :		GasPhaseChemicalRates,			&
									GetSatVapBurden,				&
									LocateGridPoint,				&
									GasPhaseChemGridIndex,			&
									GetGridCellChemBurden

	USE Aerosols,        ONLY :		Particles, Particle,			&
									ParticleArray, Molality

	USE ModelParameters, ONLY :		RstarMB, moles, grams, watermolecmass,	&
									ProtonIndex, HydroxyIndex,				&
									Avogadro

	USE Chemistry,       ONLY :		HowManyAqChems,HowManyAqCations,		&
									HowManyAqAnions,						&
									AqEquilibriaList,						&
									GasPhaseChemicalNames,					&
									AqPhaseChemicalNames,					&
									AqCationNames,AqAnionNames, & 
									FindChem, HowManyAqOrgChems

	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	!! External Variables, in
	INTEGER :: XP, YP, ZP, WhichRxn, HowManyLinks
	REAL*8  :: Temperature

	!! External Variables, out
	REAL*8  :: SoverH(HowManyLinks), SoverH2(HowManyLinks), InitialConcs1(HowManyLinks), &
			   InitialConcs2(HowManyLinks), InitialConcs3(HowManyLinks), RealFlag(HowManyLinks), GasSaturationConc(HowManyLinks)
	INTEGER :: ReturnType, ReactionType, Flag, IonA, IonB, ExpA, ExpB

	!! Internal Variables
	INTEGER :: InternalHowManyLinks, I, J, Allocation_Error, GasChemIndex
	REAL*8  :: II, XX(3), H2SO4AqIndex, AnionSatConc,  NH4SaturationConc

	TYPE(Particle), POINTER :: CurrentParticle

	!! Identify the gas index of the condensing chemical
	GasChemIndex = DissolutionData(WhichRxn,1)

	!! An initial guess
	ReturnType = -1

	!! If Water Condensation, we need a couple of extra items...
	IF (ReactionType .EQ. 0) THEN
		XX = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
		II = GetSatVapBurden(XX(1),XX(2),XX(3))
	END IF

	!! If there is a positive gas phase concentration, then return positive return flag.
	IF (ReactionType .EQ. 0.) ReturnType = 1

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Take the burdens and # of Particles counts from each of the aerosol particles !!
	CurrentParticle => Particles(XP,YP,ZP)%First !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO J = 1, HowManyLinks

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Pre-Compute the Henry's Law Coefficient !!
		!! Fill these placeholders (for non-water  !!
		!! condensation) with S' / H', which will  !!
		!! be held constant through the itnegration!!
		!!										   !!
		!!										   !!
		!! Possible reaction types: !!
		!!										   !!!!!!!!!!!!!!!
		!! 0. Water condensation (a special case for which		!!
		!!    we don't specify H's)								!!
		!!														!!
		!! 1. Direct dissolution of a non-dissociating species  !!
		!!														!!
		!! 2. Dissolution of a dissociating	species				!!
		!!														!!
		!! 3. Coupled reactions of types 1 & 2					!!
		!!														!!
		!! 4. Dissolution of a species that reacts with an ion  !!
		!!    upon entry.										!!
		!!														!!
		!! 5. Coupled reactions of types 1 & 4					!!
		!! 7. Sulfate condensation
		!!NOTE: 6 (gas-solid reactions) is skipped
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		SELECT CASE (ReactionType)

		!! Water condensing
		CASE (0)

			SoverH(J) = II * SatVapPressRatio (CurrentParticle, 1)

			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) * Avogadro

			!! Put the total solute concentration here...
			InitialConcs2(J) = 0.
			DO I = HowManyAqChems+1, HowManyAqChems + HowManyAqCations + HowManyAqAnions
				InitialConcs2(J) = InitialConcs2(J) + CurrentParticle%AqChems(I)
			END DO
			DO I = 1, HowManyAqOrgChems 
				InitialConcs2(J) = InitialConcs2(J) + CurrentParticle%AqOrgChems(I)
			END DO
			InitialConcs2(J) = InitialConcs2(J) * Avogadro

			InitialConcs3(J) = 0.

		!! Direct dissolution
		CASE (1)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) .GT. 0.)	&
			ReturnType = 1

			!! NOTE: Factor of 1000 has units of g H2O/L soln (approximate!)
			SoverH(J) =																	 &
					SatVapPressRatio (CurrentParticle, INT(DissolutionData(WhichRxn,3))) &
					* 1000. / HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature	 &
					/ CurrentParticle%AqChems(1) / moles * grams / WaterMolecMass

			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) * Avogadro
			InitialConcs2(J) = 0.
			InitialConcs3(J) = 0.

		!! Dissociating dissolution
		CASE (2)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,4)) .GT. 0.)						&
			ReturnType = 1

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! The reaction doesn't use  the Henry's Law formulation		!!
			!! because we assume it condenses directly into the dissociated	!!
			!! form, which we must now formulate...	!!!!!!!!!!!!!!!!!!!!!!!!!!
			IonA = INT(DissolutionData(WhichRxn,3))
			IonB = INT(DissolutionData(WhichRxn,4))
			ExpA = INT(DissolutionData(WhichRxn,13))
			ExpB = INT(DissolutionData(WhichRxn,14))

			!! This one is to be multiplied by Y(IonA)**ExpA * Y(IonB)**ExpB
			!! NOTE: Factor of 1000 has units of g H2O/L soln (approximate!)
			SoverH(J) =																						 &
						SatVapPressRatio (CurrentParticle , INT(DissolutionData(WhichRxn,3))) 				 &
						* (1000. * grams / moles / WaterMolecMass / CurrentParticle%AqChems(1))**(ExpA+ExpB) &
						/ HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature							 & 
						*(CurrentParticle%GammaMixed(DissolutionData(WhichRxn,12))/Avogadro)**				 & ! 3 lines of
						(AqEquilibriaList(DissolutionData(WhichRxn,12),4)+									 & ! mean act
						 AqEquilibriaList(DissolutionData(WhichRxn,12),5))									 &
						 * Avogadro 
			
			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) * Avogadro
			InitialConcs2(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,4)) * Avogadro
			InitialConcs3(J) = 0.
			
		!! Dissociating and direct dissolution (NOT USED! Matt Alvarado 10/04/06)
		CASE (3)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,4)) .GT. 0.)						&
			ReturnType = 1

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! The reaction doesn't use  the Henry's Law formulation		!!
			!! because we assume it condenses directly into the dissociated	!!
			!! form, which we must now formulate...	!!!!!!!!!!!!!!!!!!!!!!!!!!
			IonA = INT(DissolutionData(WhichRxn,3))
			IonB = INT(DissolutionData(WhichRxn,4))
			ExpA = INT(DissolutionData(WhichRxn,13))
			ExpB = INT(DissolutionData(WhichRxn,14))

			!! This one is to be multiplied by Y(IonA)**ExpA * Y(IonB)**ExpB
			SoverH(J) =																						 &
						SatVapPressRatio (CurrentParticle , INT(DissolutionData(WhichRxn,3))) 				 &
						* (1000. * grams / moles / WaterMolecMass / CurrentParticle%AqChems(1))**(ExpA+ExpB) &
						/ HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature							 & 
						*(CurrentParticle%GammaMixed(DissolutionData(WhichRxn,12))/Avogadro)**				 & ! 3 lines of
						(AqEquilibriaList(DissolutionData(WhichRxn,12),4)+									 & ! mean act
						 AqEquilibriaList(DissolutionData(WhichRxn,12),5))									 &
						 * Avogadro 
	
			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) * Avogadro
			InitialConcs2(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,4)) * Avogadro


			SoverH2(J) = SatVapPressRatio (CurrentParticle,								 &
					INT(DissolutionData(DissolutionData(WhichRxn,15),3)))				 &
					* 1000. / HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature	 &
					/ CurrentParticle%AqChems(1) / moles * grams / WaterMolecMass		 

			!! Load the Y vector with the appropriate aerosol !!
			!! concentration.								  !!
			InitialConcs3(J) = CurrentParticle%AqChems(DissolutionData(DissolutionData(WhichRxn,15),3)) * Avogadro


		!! Dissolution that incorporates an aqueous ion (NH3(g) & H+(aq) => NH4+(aq))
		CASE (4)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,2)) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) .GT. 0.)						&
			ReturnType = 1

			IonA = FLOOR(DissolutionData(WhichRxn,8))
			IonB = FLOOR(DissolutionData(WhichRxn,9))
			ExpA = ANINT(10000.*(DissolutionData(WhichRxn,8)-FLOOR(DissolutionData(WhichRxn,8))))
			ExpB = ANINT(10000.*(DissolutionData(WhichRxn,9)-FLOOR(DissolutionData(WhichRxn,9))))

			IF (CurrentParticle%AqChems(DissolutionData(WhichRxn,2)) .LE. 0) THEN
				CALL ERROR ("A dissolution reaction (specified in Dissolution.in) specifies that the gas "//		&
							TRIM(GasPhaseChemicalNames(DissolutionData(WhichRxn,1)))//								&
							" dissolves and incorporates an aqueous ion on the way in.  However, that ion "//		&
							"has a zero concentration, which we think shouldn't happen and is a sign of worse things happening.")
			ELSE

				!! Should be multiplied by IonA and Divided by IonB
				SoverH(J) =																	 &
						SatVapPressRatio (CurrentParticle, INT(DissolutionData(WhichRxn,3))) &
						*CurrentParticle%GammaMixed(IonA)**ExpA								 & 
						/CurrentParticle%GammaMixed(IonB)**ExpB								 &
						/ HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature			 &
						* Avogadro 
				!WRITE(*,*) IonA, CurrentParticle%GammaMixed(IonA)
				!WRITE(*,*) IonB, CurrentParticle%GammaMixed(IonB), CurrentParticle%IonicStr
				!PAUSE

				!! Load the Y vector with the appropriate aerosol    !!
				!! concentration.									 !!
				InitialConcs1(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) * Avogadro
				InitialConcs2(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,2)) * Avogadro
				InitialConcs3(J) = 0.
			
			END IF

		!! Dissolution that incorporates an aqueous ion (NH3(g) & H+(aq) => NH4+(aq)) coupled with direct diss.
		CASE (5)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,2)) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) .GT. 0.)						&
			ReturnType = 1

			IonA = FLOOR(DissolutionData(WhichRxn,8))
			IonB = FLOOR(DissolutionData(WhichRxn,9))
			ExpA = ANINT(10000.*(DissolutionData(WhichRxn,8)-FLOOR(DissolutionData(WhichRxn,8))))
			ExpB = ANINT(10000.*(DissolutionData(WhichRxn,9)-FLOOR(DissolutionData(WhichRxn,9))))

			IF (CurrentParticle%AqChems(DissolutionData(WhichRxn,2)) .LE. 0) THEN
				CALL ERROR ("A dissolution reaction (specified in Dissolution.in) specifies that the gas "//	&
							TRIM(GasPhaseChemicalNames(DissolutionData(WhichRxn,1)))//							&
							" dissolves and incorporates an aqueous ion on the way in.  However, that ion "//   &
							"has a zero concentration, which we think shouldn't happen and is a sign of worse things happening.")
			ELSE

				!! Should be multiplied by IonA and Divided by IonB
				SoverH(J) =																	 &
						SatVapPressRatio (CurrentParticle, INT(DissolutionData(WhichRxn,3))) &
						*CurrentParticle%GammaMixed(IonA)**ExpA								 & 
						/CurrentParticle%GammaMixed(IonB)**ExpB								 &
						/ HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature			 &
						* Avogadro 
	
				!! Load the Y vector with the appropriate aerosol    !!
				!! concentration.									 !!
				InitialConcs1(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) * Avogadro
				InitialConcs2(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,2)) * Avogadro


				SoverH2(J) = SatVapPressRatio (CurrentParticle,								 &
						INT(DissolutionData(DissolutionData(WhichRxn,15),3)))				 &
						* 1000. / HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature	 &
						/ CurrentParticle%AqChems(1) / moles * grams / WaterMolecMass

				!! Load the Y vector with the appropriate aerosol !!
				!! concentration.								  !!
				InitialConcs3(J) = CurrentParticle%AqChems(DissolutionData(DissolutionData(WhichRxn,15),3)) * Avogadro

			END IF
		
		!! Sulfate condensing
		CASE (7)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(DissolutionData(WhichRxn,4)) .GT. 0.)						&
			ReturnType = 1

			H2SO4AqIndex = FindChem("H2SO4", 1)
			SoverH(J) = 1.0

			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,3)) * Avogadro
			InitialConcs2(J) = CurrentParticle%AqChems(DissolutionData(WhichRxn,4)) * Avogadro
			InitialConcs3(J) = 0.
		
		END SELECT

		CurrentParticle => CurrentParticle%Next
	END DO

	RETURN
END SUBROUTINE DissolutionFactors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate a JACOBIAN for Condensation !!
!!                                       !!
!! As per LSODES's specification, only   !!
!! calculate at the stated row I         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CondensationJacobianEvaluator (neq, t, y, I, ia, ja, pdj)

	USE GridPointFields,     ONLY : GasPhaseChemicalRates
	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	!! External Variables
	REAL*8		:: t, y, pdj
	INTEGER		:: I, neq, ia, ja
	DIMENSION	:: y(1), ia(1), ja(1), pdj(1)

	!! Internal Variables
	INTEGER :: J, HowManyLinks

	LOGICAL, PARAMETER :: SCAFFOLDING = .FALSE.

	!WRITE(*,*) "Entering Condensation Jacobian."
	!! Determine HowManyLinks by backing out of the -1 reference point
	IF (Y(2*neq+9) .EQ. -1) THEN
		HowManyLinks = (neq - 1)/3
	ELSE IF (Y(3*neq+8) .EQ. -1) THEN
		HowManyLinks = (neq - 1)/2
	ELSE IF (Y(6*neq+5) .EQ. -1) THEN
		HowManyLinks = neq - 1
	ELSE
		CALL ERROR ("Couldn't back howmanylinks out of neq in CondensationODEEvaluator()...")
	END IF

	!!!NOTE: I = 1 is the gas-phase equation, I .GT. 1 is an aerosol index (which repeats
	!!!if more than 1 aerocol phase conpound is considered.)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! If this is for the gas phase !!
	!! chemical concentration then  !!
	!! it is different              !!
	IF (I .EQ. 1) THEN  !!!!!!!!!!!!!!

		!! Every case is the same here
		DO J = 2, HowManyLinks+1
			PDJ(J)                = Y(3*HowManyLinks+6+J) !! k_i 

			IF (Y(6*HowManyLinks+10) .GT. 1) &
				PDJ(J+HowManyLinks) = Y(3*HowManyLinks+6+J) !! k_i 

			IF (Y(6*HowManyLinks+10) .EQ. 3 .OR. Y(6*HowManyLinks+10) .EQ. 5) &
				PDJ(J+2*HowManyLinks) = Y(3*HowManyLinks+6+J) !! k_i 
		END DO

		DO J = 2, HowManyLinks+1
			PDJ(1) = PDJ(1) - PDJ(J) * Y(4*HowManyLinks+6+J)  !! - k_i * (# of Particles)
		END DO


	!! There are several structures for the columns beyond the first
	!! that depend on reaction type
	ELSE

		SELECT CASE (INT(Y(6*HowManyLinks+10)))

		!! WATER CONDENSATION
		CASE (0)
				     !!		       k_j                Henry's Law  
			PDJ(I) = -1 * Y(3*HowManyLinks+6+I) * Y(5*HowManyLinks+6+I) *		&
					 (1.-Y(I)/(Y(I)+Y(I+HowManyLinks)))/(Y(I)+Y(I+HowManyLinks))  !! - d2 c_j / dt d c_j 

					 !!			   (# of Particles)
			PDJ(1) = -1. * PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

		CASE (1) !! Dissolution of type one 

			!Normal Dissolution
					!!		       k_j                Henry's Law  
				PDJ(I) = -1 * Y(3*HowManyLinks+6+I) * Y(5*HowManyLinks+6+I)

					!!			       (# of Particles)
				PDJ(1) = -1. * PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 
		
		CASE (2) !! Dissolution of type two

			!derivative wrt H+
			IF (I .LE. HowManyLinks+1) THEN

				
				IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 0 &
				   .OR. INT(Y(7*HowManyLinks+11+I)) .EQ. 1) THEN
						 !!		       k_j                Henry's Law  
					PDJ(I) = -1 * Y(3*HowManyLinks+6+I) * Y(5*HowManyLinks+6+I) *			&
						 DissolutionData(Y(6*HowManyLinks+12),13)*Y(I)					&
						 **(DissolutionData(Y(6*HowManyLinks+12),13)-1.) *				&
						 Y(I+HowManyLinks)**DissolutionData(Y(6*HowManyLinks+12),14)

					PDJ(I+HowManyLinks) = PDJ(I)

						  !!	              (# of Particles)
					PDJ(1) = -1. * PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

				!Use calculated sat anion conc (Y(I+2*HowManyLinks))
				ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 2) THEN
						 !!		       k_j                Henry's Law  
					PDJ(I) = -1 * Y(3*HowManyLinks+6+I) * Y(5*HowManyLinks+6+I) *			&
						 DissolutionData(Y(6*HowManyLinks+12),13)*Y(I)					&
						 **(DissolutionData(Y(6*HowManyLinks+12),13)-1.) *				&
						 Y(I+2*HowManyLinks)**DissolutionData(Y(6*HowManyLinks+12),14)

					PDJ(I+HowManyLinks) = PDJ(I)

						  !!	              (# of Particles)
					PDJ(1) = -1. * PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

				!Gas Solid reaction, all terms = 0
				ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 3) THEN
					
					PDJ(I) = 0.0
					PDJ(I+HowManyLinks) = PDJ(I)

						  !!	              (# of Particles)
					PDJ(1) = -1. * PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

				END IF 
			
			!Derivative wrt anion (NO3- or Cl-)
			ELSE
				!WRITE(*,*) "1608 of CondInteg",Y(6*HowManyLinks+11+I), INT(Y(6*HowManyLinks+11+I))
				IF(INT(Y(6*HowManyLinks+11+I)) .EQ. 0 &
				   .OR. INT(Y(6*HowManyLinks+11+I)) .EQ. 1) THEN

						 ! !		       k_j                 Henry's Law  
					PDJ(I) = -1. * Y(2*HowManyLinks+6+I) * Y(4*HowManyLinks+6+I) *			&
						 Y(I-HowManyLinks)**DissolutionData(Y(6*HowManyLinks+12),13) *	&
						 DissolutionData(Y(6*HowManyLinks+12),14)*Y(I)					&
						 **(DissolutionData(Y(6*HowManyLinks+12),14)-1.)

					PDJ(I-HowManyLinks) = PDJ(I)
						 !!			            (# of Particles)
					PDJ(1) = -1. * PDJ(I) * Y(3*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

				!Use calculated sat anion conc (Y(I+2*HowManyLinks)), so all 0
				ELSE IF(INT(Y(6*HowManyLinks+11+I)) .EQ. 2 &
				       .OR. INT(Y(6*HowManyLinks+11+I)) .EQ. 3) THEN
				
					PDJ(I) = 0.0
					PDJ(I+HowManyLinks) = PDJ(I)

						  !!	              (# of Particles)
					PDJ(1) = -1. * PDJ(I) * Y(3*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

				END IF 
			
			END IF

		CASE (3) !! Dissolution of type two coupled to type one

			IF (I .LE. HowManyLinks+1) THEN

						 !!		       k_j                Henry's Law  
				PDJ(I) = -1 * Y(3*HowManyLinks+6+I) * Y(5*HowManyLinks+6+I) *			&
						 DissolutionData(Y(6*HowManyLinks+12),13)*Y(I)					&
						 **(DissolutionData(Y(6*HowManyLinks+12),13)-1.) *				&
						 Y(I+HowManyLinks)**DissolutionData(Y(6*HowManyLinks+12),14)

				PDJ(I+HowManyLinks) = PDJ(I)

						 !!	              (# of Particles)
				PDJ(1) = -1. * PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

			ELSE IF (I .LE. 2*HowManyLinks+1) THEN

						 !!		       k_j                 Henry's Law  
				PDJ(I) = -1. * Y(2*HowManyLinks+6+I) * Y(4*HowManyLinks+6+I) *			&
						 Y(I-HowManyLinks)**DissolutionData(Y(6*HowManyLinks+12),13) *	&
						 DissolutionData(Y(6*HowManyLinks+12),14)*Y(I)					&
						 **(DissolutionData(Y(6*HowManyLinks+12),14)-1.)

				PDJ(I-HowManyLinks) = PDJ(I)
						 !!			                              (# of Particles)
				PDJ(1) = -1. * PDJ(I) * Y(3*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

			ELSE

						 !!		       k_j                Henry's Law  
				PDJ(I) = -1 * Y(HowManyLinks+6+I) * Y(3*HowManyLinks+6+I)

						 !!			       (# of Particles)
				PDJ(1) = -1. * PDJ(I) * Y(2*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

			END IF

		CASE (4)

			IF (I .LE. HowManyLinks+1) THEN !Derivative w.r.t NH4+

				IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 0 &
				   .OR. INT(Y(7*HowManyLinks+11+I)) .EQ. 1) THEN
						 
							!!		       k_j                Henry's Law  
					PDJ(I) = -1. * Y(3*HowManyLinks+6+I) * Y(5*HowManyLinks+6+I) / Y(I+HowManyLinks)

					PDJ(I+HowManyLinks) = -1 * PDJ(I)

						 !!			                                 (# of Particles)
					PDJ(1) =  PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

				ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 2) THEN
					!Using coupled NH4 saturation conc
					PDJ(I) = 0.0
					PDJ(I+HowManyLinks) = -1 * PDJ(I)

						 !!			                                 (# of Particles)
					PDJ(1) =  PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 
				
				ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 3) THEN
					!Using gas-solid reaction
					PDJ(I) = 0.0
					PDJ(I+HowManyLinks) = -1 * PDJ(I)

						 !!			                                 (# of Particles)
					PDJ(1) =  PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 
			
				
				END IF
			
			ELSE !Derivative w.r.t. H+
				
				IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 0 &
				   .OR. INT(Y(7*HowManyLinks+11+I)) .EQ. 1) THEN
							!!		       k_j                Henry's Law  
					PDJ(I) = Y(2*HowManyLinks+6+I) * Y(4*HowManyLinks+6+I) * Y(I-HowManyLinks) / Y(I) / Y(I)

					PDJ(I-HowManyLinks) = -1. * PDJ(I)

						 !!			              (# of Particles)
					PDJ(1) =  PDJ(I-HowManyLinks) * Y(3*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 
				ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 2) THEN
						!!		       k_j                Henry's Law         Calculated Saturation NH4 conc
					PDJ(I) = Y(2*HowManyLinks+6+I) * Y(4*HowManyLinks+6+I) * Y(I+HowManyLinks) / Y(I) / Y(I)

					PDJ(I-HowManyLinks) = -1. * PDJ(I)

						 !!			              (# of Particles)
					PDJ(1) =  PDJ(I-HowManyLinks) * Y(3*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 
				ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 3) THEN
					!Using gas-solid reaction
					PDJ(I) = 0.0
					PDJ(I+HowManyLinks) = -1 * PDJ(I)

						 !!			                                 (# of Particles)
					PDJ(1) =  PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 		
				END IF
			END IF


		CASE (5)

			!dF/dc(NH4+)
			IF (I .LE. HowManyLinks+1) THEN

						 !!		       k_j                Henry's Law  
				PDJ(I) = -1. * Y(3*HowManyLinks+6+I) * Y(5*HowManyLinks+6+I) / Y(I+HowManyLinks)
				PDJ(I+HowManyLinks) = -1 * PDJ(I)

						 !!			                                 (# of Particles)
				PDJ(1) =  PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

			!dF/dc(H+)
			ELSE IF (I .GT. HowManyLinks+1 .AND. I .LE. 2*HowManyLinks+1) THEN

						 !!		 k_j                Henry's Law  
				PDJ(I) = Y(2*HowManyLinks+6+I) * Y(4*HowManyLinks+6+I) * Y(I-HowManyLinks) / Y(I) / Y(I)
				PDJ(I-HowManyLinks) = -1. * PDJ(I)

						 !!			              (# of Particles)
				PDJ(1) =  PDJ(I-HowManyLinks) * Y(3*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

			!dF/dc(NH3(aq))
			ELSE

						 !!		       k_j                Henry's Law  
				PDJ(I) = -1 * Y(HowManyLinks+6+I) * Y(3*HowManyLinks+6+I)

						 !!			       (# of Particles)
				PDJ(1) = -1. * PDJ(I) * Y(2*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

			END IF

		!! SULFATE CONDENSATION
		!! Assumes equilibrium gas-phase concentration is 0.0
		CASE (7)
				      
			IF (I .LE. HowManyLinks+1) THEN

						 
				PDJ(I) = 0.0

				PDJ(I+HowManyLinks) = PDJ(I)

						 !!	              (# of Particles)
				PDJ(1) = -1. * PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

			ELSE

						 
				PDJ(I) = 0.0

				PDJ(I-HowManyLinks) = PDJ(I)
						 !!			                              (# of Particles)
				PDJ(1) = -1. * PDJ(I) * Y(3*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 

			END IF

		END SELECT
		END IF

	RETURN
END SUBROUTINE CondensationJacobianEvaluator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This program provides the ODE ydot vector to LSODES. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CondensationODEEvaluator (neq, t, y, ydot)

		USE Aerosols,            ONLY : Particles, Particle, ParticleArray
		USE Chemistry,			 ONLY : HowManyGasChems
		USE InfrastructuralCode, ONLY : Transcript, ERROR, REAL2STR, INT2STR
		USE ModelParameters,     ONLY : Avogadro

		IMPLICIT NONE

		!! External Variable (prescribed by LSODES)
		REAL*8		:: t, y, ydot
		INTEGER		:: neq
		DIMENSION	:: y(1), ydot(1)

		!! Internal Variables
		INTEGER :: I, HowManyLinks

		LOGICAL, PARAMETER :: Scaffolding =  .FALSE.

		IF (SCAFFOLDING) print *, ">>> In CondensationODEEvaluator <<<"

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICAL																!!
		!! The Y-Vector we will integrate with and push through LSODES looks like this:								!!
 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!    (1)								    : Gas Phase Concenatration of Chemical DissolutionData(Rxn,1)	!!
		!!    (2 - HowManyLinks+1)				    : Aerosol Concentrations of Chemical DissolutionData(Rxn,3)		!!
		!!	  (HowManyLinks+2   - 2*HowManyLinks+1)	: Aerosol Concentrations of Chemical DissolutionData(Rxn,2) or	!!
		!!											  DissolutionData(Rxn,4), as appropriate						!!
		!!    (2*HowManyLinks+2 - 3*HowManyLinks+1) : Concentration of undissociated species for coupled rxn		!!
		!!    (3*HowManyLinks+2)					: Intentionally Empty											!!
		!!    (3*HowManyLinks+3 - 3*HowManyLinks+5) : Grid Point Specification (x,y,z)								!!
		!!    (3*HowManyLinks+6 - 3*HowManyLinks+7)	: GasChemIndex 1 & 2											!!
		!!    (3*HowManyLinks+8 - 4*HowManyLinks+7) : Condensation Rates for Chemical One							!!
		!!	  (4*HowManyLinks+8 - 5*HowManyLinks+7) : How Many Particles in Each									!!
		!!	  (5*HowManyLinks+8 - 6*HowManyLinks+9) : Pre-computed effective Henry's Law Coeffs						!!
		!!    (6*HowManyLinks+10)					: Flag for what type of reaction								!!
		!!    (6*HowManyLinks+11)					: "-1" (as a reference for finding HowManyLinks)				!!
		!!    (6*HowManyLinks+12)					: Dissolution Reaction Number									!!
		!!    (6*HowManyLinks+13 - 7*HowManyLinks+12) : Pre-computed effective Henry's Law Coeffs for coupled rxn   !!
		!!    (7*HowManyLinks+13 - 8*HowManyLink+12)  : Flags for HNO3 and HCl dissolution	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! HowManyLinks which the number of aerosol / sections														!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF (Y(2*neq+9) .EQ. -1) THEN
			HowManyLinks = (neq - 1)/3
		ELSE IF (Y(3*neq+8) .EQ. -1) THEN
			HowManyLinks = (neq - 1)/2
		ELSE IF (Y(6*neq+5) .EQ. -1) THEN
			HowManyLinks = neq - 1
		ELSE
			CALL ERROR ("Couldn't back howmanylinks out of neq in CondensationODEEvaluator()...")
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! The gas phase concentration changes as !!
		!! the negative of the sum of the other	  !!
		!! rates. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		YDOT(1) = 0.

		!WRITE(*,*) "Case: ", Y(6*HowManyLinks+10), 
		SELECT CASE (INT(Y(6*HowManyLinks+10)))

		!! WATER CONDENSATION
		CASE (0)

			DO I = 2, neq
				YDOT(I) = Y(3*HowManyLinks+6+I) *			& ! The Condensation Rate
						  (Y(1) -							& ! The Gas Phase Conc
						   Y(I)/(Y(I)+Y(I+HowManyLinks))*	& ! Raout's Law
						   Y(5*HowManyLinks+6+I))			  ! S_w' * C_s,w

				YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I)

			END DO

		!! DISSOLUTION of a single, non-dissociating species
		CASE (1)

			DO I = 2, neq
				
				YDOT(I) = Y(3*HowManyLinks+6+I) *	&		! The Condensation Rate
						  (Y(1) -					&		! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) * Y(I))	! S_i' * c_s,i / H_i

				YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I)
				
			END DO


		!! DISSOLUTION of a single, dissociating species
		CASE (2)



			DO I = 2, neq-HowManyLinks
				
				IF( INT(Y(7*HowManyLinks+11+I)) .EQ. 0 &
			        .OR. INT(Y(7*HowManyLinks+11+I)) .EQ. 1) THEN

				     YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1) -					& ! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) *	& ! S/H'
   						   Y(I)**DissolutionData(Y(6*HowManyLinks+12),13) *				&
						   Y(I+HowManyLinks)**DissolutionData(Y(6*HowManyLinks+12),14))

				     YDOT(I+HowManyLinks) = YDOT(I)

				     YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I) / DissolutionData(Y(6*HowManyLinks+12),13)
				
				!Uses calculated saturation anion conc rather than real anion conc
				ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 2) THEN
					  
					 YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1) -					& ! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) *	& ! S/H'
   						   Y(I)**DissolutionData(Y(6*HowManyLinks+12),13) *				&
						   Y(I+2*HowManyLinks)**DissolutionData(Y(6*HowManyLinks+12),14))

				     YDOT(I+HowManyLinks) = YDOT(I)

				     YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I) / DissolutionData(Y(6*HowManyLinks+12),13)
				
				ELSE IF (INT(Y(7*HowManyLinks+11+I)) .EQ. 3) THEN
					YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1) -					& ! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) ) !Saturation Gas Phase Conc (for gas-solid rxn)

							!IF(SCAFFOLDING) 
							!WRITE(*,*) "Gas Driving Force: ", Y(1)-Y(5*HowManyLinks+6+I)
							!WRITE(*,*) "Gas Conc  ", TRIM(REAL2STR(Y(1)))
							!WRITE(*,*) "Gas-solid Sat Conc.  ", TRIM(REAL2STR(Y(5*HowManyLinks+6+I)))
					 YDOT(I+HowManyLinks) = YDOT(I)

				     YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I) / DissolutionData(Y(6*HowManyLinks+12),13)

				END IF
			
			END DO

		!! DISSOLUTION of a coupled reactions of types 1 & 2
		CASE (3)

			DO I = 2, 1+HowManyLinks
				YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1) -					& ! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) *	&
   						   Y(I)**DissolutionData(Y(6*HowManyLinks+12),13) *				&
						   Y(I+HowManyLinks)**DissolutionData(Y(6*HowManyLinks+12),14))

				YDOT(I+HowManyLinks) = YDOT(I)

				YDOT(2*HowManyLinks+I)	=			&
						  Y(3*HowManyLinks+6+I) *	&		! The Condensation Rate
						  (Y(1) -					&		! The Gas Phase Conc
						   Y(6*HowManyLinks+11+I)	&		! S_i' * c_s,i / H_i
						   * Y(2*HowManyLinks+I))

				YDOT(1) = YDOT(1) - (YDOT(I) / DissolutionData(Y(6*HowManyLinks+12),13) + &
						  YDOT(2*HowManyLinks+I)) * Y(4*HowManyLinks+6+I) 
			END DO


		!! DISSOLUTION of a single species that bonds with an ion on the way in
		CASE (4)

			IF( INT(Y(7*HowManyLinks+11+I)) .EQ. 0 &
			        .OR. INT(Y(7*HowManyLinks+11+I)) .EQ. 1) THEN

				DO I = 2, 1+HowManyLinks
					YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1) -					& ! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) *	&
   						   Y(I) / Y(I+HowManyLinks)) 

					YDOT(I+HowManyLinks) = -1 * YDOT(I)

					YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I)
				END DO
			
			ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 2) THEN
			!Uses saturation NH4 concentration
				DO I = 2, 1+HowManyLinks
					YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1) -					& ! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) *	&
   						   Y(I+2*HowManyLinks) / Y(I+HowManyLinks)) 

					YDOT(I+HowManyLinks) = -1 * YDOT(I)

					YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I)
				END DO
			ELSE IF(INT(Y(7*HowManyLinks+11+I)) .EQ. 3) THEN
			!Uses Gas-solid reaction
				DO I = 2, 1+HowManyLinks
					YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1) -					& ! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I)) 

					YDOT(I+HowManyLinks) = -1 * YDOT(I)

					YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I)
				END DO
			
			END IF
			
		!! Coupled reactions of types 1 & 4
		CASE (5)

			DO I = 2, 1+HowManyLinks
				!NH4+ formed by ion binding reaction
				YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1) -					& ! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) *	& ! S/H'
   						   Y(I) / Y(I+HowManyLinks))  ! c(NH4+)/c(H+)

				IF (SCAFFOLDING) WRITE(*,*) "ki", Y(3*HowManyLinks+6+I)
				!IF (SCAFFOLDING) WRITE(*,*) "Gas Conc (molec/cm3)", Y(1)
				!IF (SCAFFOLDING) WRITE(*,*) "S/H", Y(5*HowManyLinks+6+I)
				IF (SCAFFOLDING) WRITE(*,*) "Conc ratio", Y(I) / Y(I+HowManyLinks)
				IF (SCAFFOLDING) WRITE(*,*) "Grad", Y(1) - Y(5*HowManyLinks+6+I) *	& ! S/H'
   													Y(I) / Y(I+HowManyLinks)
				IF (SCAFFOLDING) WRITE(*,*) "dH+/dt", -1.*YDOT(I)
				IF (SCAFFOLDING) WRITE(*,*) "H+ conc", Y(I+HowManyLinks)
				IF (SCAFFOLDING) WRITE(*,*) "NH+ conc", Y(I)
				
				
				!H+ destroyed by ion binding reaction
				YDOT(I+HowManyLinks) = -1 * YDOT(I)

				!NH3(aq) formed by direct dissolution
				YDOT(2*HowManyLinks+I)	=			&
						  Y(3*HowManyLinks+6+I) *	&		! The Condensation Rate
						  (Y(1) -					&		! The Gas Phase Conc
						   Y(6*HowManyLinks+11+I)	&		! S_i' / H_i
						   * Y(2*HowManyLinks+I))           ! c(NH3(aq))

				YDOT(1) = YDOT(1) - (YDOT(I) + YDOT(2*HowManyLinks+I)) * Y(4*HowManyLinks+6+I) 

			END DO

		!! SULFATE CONDENSATION
		!! Assumes that H2SO4 vapor pressure = 0 
		CASE (7)

			DO I = 2, neq-HowManyLinks
				!WRITE(*,*) Y(3*HowManyLinks+6+I)
				YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1)) * &					! The Gas Phase Conc
						   DissolutionData(Y(6*HowManyLinks+12),13) !Cation Stoichiometry
				YDOT(I+HowManyLinks) = YDOT(I)*(DissolutionData(Y(6*HowManyLinks+12),14)/DissolutionData(Y(6*HowManyLinks+12),13))

				YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I) / DissolutionData(Y(6*HowManyLinks+12),13)
			END DO

		END SELECT

	RETURN
END SUBROUTINE CondensationODEEvaluator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Call LSODES or the equilibrium routine to integrate forward	 !!
!! the concentration of aerosol / gas-phase water at the grid-	 !!
!! cell level.  It must construct a joint array of all of the	 !!
!! aerosol plus the grid cell for each of the condensing species !!
!! XP, YP, and ZP are the gridpoint numbers in gas grid space.	 !!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE StepCondensation (TimestepSize, NumTimeSteps, XP, YP, ZP)

	USE InfrastructuralCode, ONLY : INT2STR, ERROR,					&
									Warn, Transcript,				&
									REAL2STR

	USE GridPointFields,	 ONLY : LocateGridPoint,				&
									getM,							&
									GasPhaseChemGridIndex,			&
									GetGridCellChemBurden,			&
									ReplaceGridCellChemBurden,		&
									GetTemp,						&
									GetRelativeHumidity,			&
									    AddToGridCellChemBurden, &
                                                                        GetSatVapBurden

	USE Aerosols,            ONLY : Particles, Particle,			&
									ParticleArray,					&
									RecalculateRadius,				&
									EqCoeff,						&
									Molality,						&
									FindElectrolyteEquilibrium,		&
									CalculateIonicStrength,			&
									KusikMeissner,					&
									SortAerosolAtGridPointForCondensation,			&
									EquilibriumWaterContent, &
									FindAqElectrolyteEquilibriumForGridPoint, &
									EquilibrateAqueousDissociationReaction

	USE Chemistry,			 ONLY : HowManyGasChems,				&
									HowManyEvolveGasChems,			&
									HowManyEvolveAqChems,			&
									GasPhaseChemicalNames,			&
									HowManyCondensingChems,			&
									AqEquilibriaList,				&
									HowManyAqChems,					&
									HowManyAqCations,				&
									FindChem,						&
									AqPhaseChemicalNames, &
									HowManyAqEqReactions

	USE ModelParameters,     ONLY : Avogadro,								&
									AerosolWaterEquilibriumRHThreshold,		&
									WaterEquilibriumIndex,					&
									ThermoBulkMode,							&
									DissolutionEquilibriumOrFlux,			&
									AqThermoEquilibriumError,				&
									WaterContentPrecision,					&
									micron,									&
									AmmoniaFlux,							&
									DoHydrophobicOrgDissolutionFlag,        &
									DoHydrophilicOrgDissolutionFlag,		&
									ProtonIndex, &
									CombinedInorgDiss, &
									ZSROrgWater


	IMPLICIT NONE

	!! External Variables
	REAL*8  :: TimeStepSize		!! Should be in Seconds
	INTEGER :: NumTimeSteps, XP, YP, ZP

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!! LSODES Related Internal Variables !!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8,  ALLOCATABLE	:: y(:), rwork(:)
	REAL*8					:: t,outtime,rtol,atol
	INTEGER, ALLOCATABLE	:: iwork(:)
	INTEGER					:: itol,itask,istate,iopt,mf,neq,LIW,LRW
		
	!! Non-LSODES Local Variables
	INTEGER				:: I, J, K, L, LLLL, C, allocation_error, Rxn, HowManyLinks, YLen, GG,					&
						   HowManyNonZeroJacobianTerms, GasChemIndex1, TimeStep,	&
						   ReturnType, WhichDominantIon, WhichSecondIon , Flag, Dum, RxnType, &
						   InternalReturnType, NH4Index, HIndex, OHIndex, GasSolidRxn
	REAL*8				:: II, JJ, KK, MM, ChemTransferArray(HowManyGasChems), Temperature,	&
									    RH, XX(3), NH3Index, OrgWaterContent, Store, Store3,Store4, GasPhaseBurden, Electrolyte, WaterSatBurden
	REAL*8, ALLOCATABLE :: InitialConcsArray(:), ChemStorageArray(:,:)
	LOGICAL :: EmptyGridCell
	REAL*8  :: BeginTime,EndTime, ABCD
	REAL*4  :: etimearray(2)
	REAL*8  :: WaterTime, CondTime, AmmoniaTime

	LOGICAL, PARAMETER :: Scaffolding = .FALSE., OpenSystem = .FALSE., ForceEquilibrium = .TRUE.
	LOGICAL :: Equilibrated, UseEq

	TYPE(Particle), POINTER :: FirstParticle, CurrentParticle

	!! SCSCSCSC
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript(">>>Entering StepCondensationRxn<<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Find GRID POINT Values for NON-AEROSOL !!
	!! DEPENDENT Parameters                   !!
	XX          = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
	Temperature = GetTemp(XX(1),XX(2),XX(3))
	RH		    = GetRelativeHumidity(XX(1),XX(2),XX(3))
		
	IF (Scaffolding) CALL TRANSCRIPT("Using LSODES to integrate dissolution.")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Everything should now be sorted           !!
	!! Back up to the head of the aerosol list   !!
	FirstParticle   => Particles(XP,YP,ZP)%First
	CurrentParticle => FirstParticle

	!! COUNT the number of LINKS in the gridcell that we should condense upon
	HowManyLinks  =  1
15	IF(ASSOCIATED(CurrentParticle%Next)) THEN
		CurrentParticle => CurrentParticle%Next
		IF (CurrentParticle%NumberOfParticles .GT. 0 .AND. .NOT.CurrentParticle%Dry) THEN
			HowManyLinks  = HowManyLinks + 1
			GOTO 15
		END IF
	ENDIF
	CurrentParticle => FirstParticle


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set the lengths of the various !!
	!! arrays and allocate them !!!!!!!!
	YLen	= 12 + HowManyLinks*8	  !! As per the definition below
	NEQ		= 1  + 3*HowManyLinks	  !! Maximum value is for reaction type 5.  Allocate based on this.

	!! Temporarily set this to the maximum value (which is also for reaction type 4)
	HowManyNonZeroJacobianTerms = 11*HowManyLinks + 1
	LRW  = 60 +  25*NEQ + 3*HowManyNonZeroJacobianTerms
	LIW  = 31 + neq + HowManyNonZeroJacobianTerms

	!! Allocate the LSODES working arrays
	ALLOCATE (Y(ylen), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of Y could not proceed in StepCondensation()")

	ALLOCATE (rwork(LRW), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of RWORK could not proceed in StepCondensation()")

	ALLOCATE (iwork(LIW), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of IWORK could not proceed in StepCondensation()")

	ALLOCATE (InitialConcsArray(NEQ), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of InitialConcsArray could not proceed in StepCondensation()")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This is indexed as follows:												!!
	!! ChemStorageArray(Rxn, 1)         = Change in Gas Phase Conc 1			!!
	!! ChemStorageArray(Rxn, 2)         = Change in Aerosol Conc 1 (1)			!!
	!! ...																		!!
	!! ChemStorageArray(Rxn, HML+1)     = Change in Aerosol Conc 1 (HowManyLinks)!
	!! ChemStorageArray(Rxn, HML+2)     = Change in Aerosol Conc 2 (1)			!!
	!! ...																		!!
	!! ChemStorageArray(Rxn, 2*HML+1)   = Change in Aerosol Conc 2 (HowManyLinks)!
	!! ChemStorageArray(Rxn, 2*HML+2)   = Change in Aerosol Conc 3 (1)			!!
	!! ...																		!!
	!! ChemStorageArray(Rxn, 3*HML+1)   = Change in Aerosol Conc 3 (HowManyLinks)!
	ALLOCATE (ChemStorageArray(HowManyDissolutionReactions,NEQ), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of ChemStorageArray could not proceed in StepCondensation()")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Initialize the LSODES Work Vectors !!
	DO I = 1, LRW !!!!!!!!!!!!!!!!!!!!!!!!!!
		RWORK(I) = 0.
	END DO
	DO I = 1, LIW
		IWORK(I) = 0
	END DO

	!! Initialize the Y vector
	DO J = 1, YLen
		Y(J) = 0.
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Take the # of Particles counts from each of the aerosol particles !!
	CurrentParticle => Particles(XP,YP,ZP)%First !!!!!!!!!!!!!!!!!!!!!!!!!!
	DO J = 1, HowManyLinks
		Y(4*HowManyLinks+7+J) = CurrentParticle%NumberOfParticles
		CurrentParticle		  => CurrentParticle%Next
	END DO


	IF (Scaffolding) CALL TRANSCRIPT("Cond. Check 1")

		
	!!!!!!!!!!!!!!!!!!!!!!!!!
	!! LOOP OVER TIMESTEPS !!
	DO TimeStep = 1, NumTimeSteps

		!! Integrate to this time
		OutTime =  TimeStep * TimeStepSize

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Initialize Chem Storage Array !!
		DO I = 1, NEQ
		DO J = 1, HowManyDissolutionReactions
			ChemStorageArray(J,I) = 0.
		END DO  ;  END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Set LSODES-Related Variables (Notes Excerpted !!
		!! from LSODES file)							 !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! itask  = an index specifying the task to be performed.				!!
		!!          input only.  itask has the following values and meanings.	!!
		!!          1  means normal computation of output values of y(t) at		!!
		!!             t = tout (by overshooting and interpolating).			!!
		!!          2  means take one step only and return.						!!
		!!          3  means stop at the first internal mesh point at or		!!
		!!             beyond t = tout and return.								!!
		!!          4  means normal computation of output values of y(t) at		!!
		!!             t = tout but without overshooting t = tcrit.				!!
		!!             tcrit must be input as rwork(1).  tcrit may be equal to	!!
		!!             or beyond tout, but not behind it in the direction of	!!
		!!             integration.  this option is useful if the problem		!!
		!!             has a singularity at or beyond t = tcrit.				!!
		!!          5  means take one step, without passing tcrit, and return.	!!
		!!             tcrit must be input as rwork(1).							!!
		!!																		!!
		!!          note..  if itask = 4 or 5 and the solver reaches tcrit		!!
		!!          (within roundoff), it will return t = tcrit (exactly) to	!!
		!!          indicate this (unless itask = 4 and tout comes before tcrit,!!
		!!          in which case answers at t = tout are returned first).		!!
		itask  = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! istate = an index used for input and output to specify the
		!!         the state of the calculation.
		!!
		!!         on input, the values of istate are as follows.
		!!          1  means this is the first call for the problem
		!!             (initializations will be done).  see note below.
		!!          2  means this is not the first call, and the calculation
		!!             is to continue normally, with no change in any input
		!!             parameters except possibly tout and itask.
		!!             (if itol, rtol, and/or atol are changed between calls
		!!             with istate = 2, the new values will be used but not
		!!             tested for legality.)
		!!          3  means this is not the first call, and the
		!!             calculation is to continue normally, but with
		!!             a change in input parameters other than
		!!             tout and itask.  changes are allowed in
		!!             neq, itol, rtol, atol, iopt, lrw, liw, mf,
		!!             the conditional inputs ia and ja,
		!!             and any of the optional inputs except h0.
		!!             in particular, if miter = 1 or 2, a call with istate = 3
		!!             will cause the sparsity structure of the problem to be
		!!             recomputed (or reread from ia and ja if moss = 0).
		!!          note..  a preliminary call with tout = t is not counted
		!!          as a first call here, as no initialization or checking of
		!!          input is done.  (such a call is sometimes useful for the
		!!          purpose of outputting the initial conditions.)
		!!          thus the first call for which tout .ne. t requires
		!!          istate = 1 on input.
		!!
		!!          on output, istate has the following values and meanings.
		!!           1  means nothing was done, as tout was equal to t with
		!!              istate = 1 on input.  (however, an internal counter was
		!!              set to detect and prevent repeated calls of this type.)
		!!           2  means the integration was performed successfully.
		!!          -1  means an excessive amount of work (more than mxstep
		!!              steps) was done on this call, before completing the
		!!              requested task, but the integration was otherwise
		!!              successful as far as t.  (mxstep is an optional input
		!!              and is normally 500.)  to continue, the user may
		!!              simply reset istate to a value .gt. 1 and call again
		!!              (the excess work step counter will be reset to 0).
		!!              in addition, the user may increase mxstep to avoid
		!!              this error return (see below on optional inputs).
		!!          -2  means too much accuracy was requested for the precision
		!!              of the machine being used.  this was detected before
		!!              completing the requested task, but the integration
		!!              was successful as far as t.  to continue, the tolerance
		!!              parameters must be reset, and istate must be set
		!!              to 3.  the optional output tolsf may be used for this
		!!              purpose.  (note.. if this condition is detected before
		!!              taking any steps, then an illegal input return
		!!              (istate = -3) occurs instead.)
		!!          -3  means illegal input was detected, before taking any
		!!              integration steps.  see written message for details.
		!!              note..  if the solver detects an infinite loop of calls
		!!              to the solver with illegal input, it will cause
		!!              the run to stop.
		!!          -4  means there were repeated error test failures on
		!!              one attempted step, before completing the requested
		!!              task, but the integration was successful as far as t.
		!!              the problem may have a singularity, or the input
		!!              may be inappropriate.
		!!          -5  means there were repeated convergence test failures on
		!!              one attempted step, before completing the requested
		!!              task, but the integration was successful as far as t.
		!!              this may be caused by an inaccurate jacobian matrix,
		!!              if one is being used.
		!!          -6  means ewt(i) became zero for some i during the
		!!              integration.  pure relative error control (atol(i)=0.0)
		!!              was requested on a variable which has now vanished.
		!!              the integration was successful as far as t.
		!!          -7  means a fatal error return flag came from the sparse
		!!              solver cdrv by way of prjs or slss (numerical
		!!              factorization or backsolve).  this should never happen.
		!!              the integration was successful as far as t.
		!!
		!!          note.. an error return with istate = -1, -4, or -5 and with
		!!          miter = 1 or 2 may mean that the sparsity structure of the
		!!          problem has changed significantly since it was last
		!!          determined (or input).  in that case, one can attempt to
		!!          complete the integration by setting istate = 3 on the next
		!!          call, so that a new structure determination is done.
		!!
		!!          note..  since the normal output value of istate is 2,
		!!          it does not need to be reset for normal continuation.
		!!          also, since a negative input value of istate will be
		!!          regarded as illegal, a negative output value requires the
		!!          user to change it, and possibly other inputs, before
		!!          calling the solver again.
		istate = 1	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! iopt   = an integer flag to specify whether or not any optional		!!
		!!          inputs are being used on this call.  input only.			!!
		!!          the optional inputs are listed separately below.			!!
		!!          iopt = 0 means no optional inputs are being used.			!!
		!!                   default values will be used in all cases.			!!
		!!          iopt = 1 means one or more optional inputs are being used.	!!
		iopt   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! mf     = the method flag.  used only for input.					!!
		!!          the standard choices for mf are..						!!
		!!            mf = 10  for a nonstiff problem,						!!
		!!            mf = 21 or 22 for a stiff problem with ia/ja supplied	!!
		!!                     (21 if jac is supplied, 22 if not),			!!
		!!            mf = 121 for a stiff problem with jac supplied,		!!
		!!                     but not ia/ja,								!!
		!!            mf = 222 for a stiff problem with neither ia/ja nor	!!
		!!                     jac supplied.								!!
		!!          the sparseness structure can be changed during the		!!
		!!          problem by making a call to lsodes with istate = 3.		!!
		mf = 21  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!h0      rwork(5)  the step size to be attempted on the first step.
        !!            the default value is determined by the solver.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		rwork(5) = 1.e-10

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! hmax    rwork(6)  the maximum absolute step size allowed.
		!!                   the default value is infinite.
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		rwork(6) = 1.0

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! hmin     = the minimum allowed timestep !!
!		rwork(7) = 1.e-10 !!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! atol   = an absolute error tolerance parameter, either a scalar or	!!
		!!          an array of length neq.  input only.						!!
		!!																		!!
		!!             the input parameters itol, rtol, and atol determine		!!
		!!          the error control performed by the solver.  the solver will	!!
		!!          control the vector e = (e(i)) of estimated local errors		!!
		!!          in y, according to an inequality of the form				!!
		!!                      rms-norm of ( e(i)/ewt(i) )   .le.   1,			!!
		!!          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),			!!
		!!          and the rms-norm (root-mean-square norm) here is			!!
		!!          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))	!!
		!!          is a vector of weights which must always be positive, and	!!
		!!          the values of rtol and atol should all be non-negative.		!!
		!!          the following table gives the types (scalar/array) of		!!
		!!          rtol and atol, and the corresponding form of ewt(i).		!!
		!!																		!!
		!!             itol    rtol       atol          ewt(i)					!!
		!!              1     scalar     scalar     rtol*abs(y(i)) + atol		!!
		!!              2     scalar     array      rtol*abs(y(i)) + atol(i)	!!
		!!              3     array      scalar     rtol(i)*abs(y(i)) + atol	!!
		!!              4     array      array      rtol(i)*abs(y(i)) + atol(i)	!!
		!!																		!!
		!!          when either of these parameters is a scalar, it need not	!!
		!!          be dimensioned in the user-s calling program.				!!
		!!																		!!
		!!          if none of the above choices (with itol, rtol, and atol		!!
		!!          fixed throughout the problem) is suitable, more general		!!
		!!          error controls can be obtained by substituting				!!
		!!          user-supplied routines for the setting of ewt and/or for	!!
		!!          the norm calculation.  see part iv below.					!!
		!!																		!!
		!!          if global errors are to be estimated by making a repeated	!!
		!!          run on the same problem with smaller tolerances, then all	!!
		!!          components of rtol and atol (i.e. of ewt) should be scaled	!!
		!!          down uniformly.												!!
		itol   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! This is a very stiff problem, as aerosol and the gas phase should 
		!! not be the same scale at all.  Screw absolute tolerance.
		!! Could also explore using a vector of atol's
		!!
		!! However, if the initial concentration is zero, some absolute tolerance
		!! is necessary (an error of 0 is unacceptable).
		!!
		!! The units of this are molecules
		atol   = 1.0

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Tell the program which gridpoint it's operating at: !!
		Y(3*HowManyLinks+3) = XP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Y(3*HowManyLinks+4) = YP
		Y(3*HowManyLinks+5) = ZP


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! IF EQUILIBRIUM RATHER THAN FLUX, DO THAT HERE AND RETURN !!
		IF (.NOT. DissolutionEquilibriumOrFlux .AND. RH .LE. AerosolWaterEquilibriumRHThreshold) THEN
			WRITE(*,*) "Calling Equilibrium Routine"
			CALL EquilibrateGridPointAll(XP, YP, ZP, .TRUE., ForceWaterEquilibrium=.FALSE., WaterOpenSystem=.FALSE.)
			RETURN
		END IF
		
		
		IF (Scaffolding) WRITE(*,*) "Before Dissolution Loop"
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! LOOP OVER EACH DISSOLUTION REACTION !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		DO Rxn = 1, HowManyDissolutionReactions

		!If we are doing the combined inorganic dissolution
		IF(CombinedInorgDiss .AND. (Rxn .GE. 2 .AND. Rxn .LE. 8)) THEN

			WRITE(*,*) "Rxn: ", Rxn

			!!When this is true, do eq for NH3, HNO3 and HCl
			!!Else, do flux mode
			UseEq = .TRUE.
	
			!Call the combined integrator only once
			IF(Rxn .EQ. 2) THEN
				IF(NOT(UseEq)) THEN
					CALL BigStepInorgCondensation (TimestepSize, NumTimeSteps, XP, YP, ZP)
				END IF
			ELSE IF(Rxn .EQ. 3 .AND. UseEq) THEN
				CALL EquilibrateInorgGridPoint(XP, YP, ZP)
				CYCLE
			ELSE	
				!Skip all the H2SO4, HNO3, HCl, and NH3 rxns
				!since they are done above
				WRITE(*,*) "Cycling"
				CYCLE
			END IF

		END IF
		
		!Send an error if this tries to do ammonia in flux form
		IF(AmmoniaFlux .EQ. .TRUE. .AND. CombinedInorgDiss .EQ. .FALSE.) &
			CALL ERROR("Operator split method only works if Ammonia is done as eq.")


		!! The allowed error for water content may be different than for
		!! the trace species.
		IF (Rxn .EQ. 1) THEN
			rtol   = WaterContentPrecision **(1./ NumTimeSteps) ! AqThermoEquilibriumError
		ELSE
			rtol   = (AqThermoEquilibriumError) **(1./NumTimeSteps) ! AqThermoEquilibriumError
		END IF

		!! If the reaction is coupled to another one, the direct dissolution
		!! is done at the same time as the other one, meaning we skip this
		!! direct one, UNLESS the reaction is NH3 & H+ => NH4+ when NH4+ is low
		!! (that's the Dum flag). Cycle...
		IF (SCAFFOLDING) WRITE(*,*) "Dum", Dum
		!IF (DissolutionData(Rxn, 15) .LT. 0) CYCLE

		!! Reset timestep to the correct position
		T = (TimeStep - 1) * TimeStepSize
		
		!! If Water is condensing and the RH Flag is below the critical value,
		!! then simply equilibrate it. (Uses total water content: org + inorg)
		IF (Rxn .EQ. 1 .AND. RH .LE. AerosolWaterEquilibriumRHThreshold) THEN
			Equilibrated = .FALSE.
			
			DO WHILE (NOT(Equilibrated))
				CurrentParticle => Particles(XP,YP,ZP)%First
				Equilibrated = .TRUE.
			
				DO J = 1, HowManyLinks
					IF (SCAFFOLDING) WRITE(*,*) "Before Equilibrate, particle #", J
				
					!Skip empty particles
					IF(CurrentParticle%NumberofParticles .LE. 0.) CYCLE

					LLLL = 1
					ReturnType = 1
					DO WHILE (ReturnType .NE. 2)
						IF (LLLL .GT. 900) WRITE(*,*) LLLL, "Equilibrating Water! Water Content = ", CurrentParticle%AqChems(1), RH
						
						ABCD = CurrentParticle%AqChems(1)
						IF(ZSROrgWater) THEN
							CALL EqAllWater(CurrentParticle, ReturnType, OpenSystem)
						ELSE
							!Ross Eqn
							CALL EquilibriumWaterContentAll (CurrentParticle, OpenSystem, ReturnType, ForceEquilibrium=.FALSE.)
						END IF

						IF(ReturnType .NE. 2) Equilibrated = .FALSE.
						IF(LLLL .GT. 50 .AND. &
						  ABS((ABCD-CurrentParticle%AqChems(1))/CurrentParticle%AqChems(1)) .LE. 0.01) THEN
						  ReturnType  = 2
						  Equilibrated = .TRUE.
						END IF
						
						!This is here to prevent the aqeuous equilibrium
						!From hanging up the loop
						!IF(ReturnType .EQ. 2 .AND. LLLL .GT. 50) THEN
						!	CALL FindElectrolyteEquilibrium(CurrentParticle, UpdateThermo = .TRUE.)
						!	EXIT
						!END IF

						!Equilibrate organics across phases	
						!IF(DoHydrophobicOrgDissolutionFlag .EQ. .TRUE. .AND. DoHydrophilicOrgDissolutionFlag .EQ. .TRUE.) THEN
						!	CALL EquilibrateOrganicParticle (CurrentParticle, Temperature, InternalReturnType)	 
						!	IF (SCAFFOLDING) WRITE(*,*) "Org Return", InternalReturnType 
						!	IF(InternalReturnType .NE. 2) THEN
						!		ReturnType = 1
						!	END IF
						!END IF
					
						!We loop over individual inorganic reactions to improve the speed
						!of the water equilibrium convergence
						!CALL FindElectrolyteEquilibrium(CurrentParticle, UpdateThermo = .TRUE., FirstEquilibration = .TRUE., &
						!			ReturnType = InternalReturnType)
						!IF(InternalReturnType .NE. 2) THEN
						!		ReturnType = 1
						!END IF
						
						!CALL KusikMeissner (CurrentParticle)
						!DO C = 1, HowManyAqEqReactions
						!WRITE(*,*) "Eq Rxn #", C
						!	CALL EquilibrateAqueousDissociationReaction (C, CurrentParticle, ReturnType = InternalReturnType)
						!	!CALL KusikMeissner (CurrentParticle)
						!	IF(InternalReturnType .NE. 2) THEN
						!		ReturnType = 1
						!	END IF
						!END DO


						LLLL = LLLL+1
						IF(LLLL .GT. 1000) THEN
							CALL WARN("Over 1000 water iterations in StepCondensation().")
							CALL WARN("RH equals :"//trim(real2str(RH)))
						END IF
					END DO
						
					!IF (SCAFFOLDING) 
					!WRITE(*,*) "Equilibrate water for particle #", J, LLLL
					CurrentParticle => CurrentParticle%Next
				END DO
			END DO
			CYCLE
		END IF
		

		!! IF the flag to equilibrate ammonia has been selected, just equilibrate it.
		!! This is for NH3(g) & H+(aq) <=> NH4+(aq) and NH3(g) <=> NH3(aq)
		IF (GasPhaseChemicalNames(DissolutionData(Rxn,1)) .EQ. "NH3" &
			.AND. (INT(DissolutionData(Rxn,11)) .EQ.  4 .OR. INT(DissolutionData(Rxn,11)) .EQ.  1) &
			.AND. AmmoniaFlux .EQ. .FALSE. ) THEN
			Equilibrated = .FALSE.
			LLLL = 1

			DO WHILE (NOT(Equilibrated))
				CurrentParticle => Particles(XP,YP,ZP)%First
				Equilibrated = .TRUE.
				DO J = 1, HowManyLinks
					IF (SCAFFOLDING) WRITE(*,*) "Before Ammonia Equilibrate, particle #", J
				
					IF (SCAFFOLDING) WRITE(*,*) "Equilibrating Ammonia!"
					!WRITE(*,*) "Particle # ", J, "Iteration: ", LLLL
					!CALL KusikMeissner (CurrentParticle)
					XX = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
					Temperature = GetTemp(XX(1),XX(2),XX(3))
					GasPhaseBurden = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Rxn,1)))
					WaterSatBurden = GetSatVapBurden (XX(1), XX(2), XX(3))
                                        CALL EquilibrateDissolutionReaction (GasPhaseBurden, Rxn, CurrentParticle, &
                                        Temperature, ReturnType, WaterSatBurden = WaterSatBurden, GasPhaseBurden2=0.0, OptErrorTolerance=AqThermoEquilibriumError/2.) 
					IF(ReturnType .NE. 2) Equilibrated = .FALSE.
					!WRITE(*,*) "ReturnType = ", ReturnType
					CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(Rxn,1)), GasPhaseBurden)		
					IF(ReturnType .NE. 2) THEN
						!CALL FindElectrolyteEquilibrium(CurrentParticle, .TRUE.)
					END IF
					IF (SCAFFOLDING) WRITE(*,*) "Equilibrate Ammonia for particle #", J
					CurrentParticle => CurrentParticle%Next
				END DO
				LLLL = LLLL + 1
				!After 20 total iterations, just stop
				IF(LLLL .GT. 20) EXIT
			END DO
			WRITE(*,*) "NH3 eq. took "//trim(real2str(AmmoniaTime))//" seconds."
			WRITE(*,*) "Number of iterations: ", LLLL
			CYCLE	
		END IF
		
		!! Identify the gas index of the condensing chemical
		GasChemIndex1 = DissolutionData(Rxn,1)
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! This is the maximum number of internal steps allowed. !!
		!! The default is 500 and that's not enough.			 !!
		iwork(6) = 300000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICAL																!!
			!! The Y-Vector we will integrate with and push through LSODES looks like this:								!!
	 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!    (1)								    : Gas Phase Concenatration of Chemical DissolutionData(Rxn,1)	!!
			!!    (2 - HowManyLinks+1)				    : Aerosol Concentrations of Chemical DissolutionData(Rxn,3)		!!
			!!	  (HowManyLinks+2   - 2*HowManyLinks+1)	: Aerosol Concentrations of Chemical DissolutionData(Rxn,2) or	!!
			!!											  DissolutionData(Rxn,4), as appropriate						!!
			!!    (2*HowManyLinks+2 - 3*HowManyLinks+1) : Concentration of undissociated species for coupled rxn		!!
			!!    (3*HowManyLinks+2)					: Intentionally Empty											!!
			!!    (3*HowManyLinks+3 - 3*HowManyLinks+5) : Grid Point Specification (x,y,z)								!!
			!!    (3*HowManyLinks+6 - 3*HowManyLinks+7)	: GasChemIndex 1 & 2											!!
			!!    (3*HowManyLinks+8 - 4*HowManyLinks+7) : Condensation Rates for Chemical One							!!
			!!	  (4*HowManyLinks+8 - 5*HowManyLinks+7) : How Many Particles in Each									!!
			!!	  (5*HowManyLinks+8 - 6*HowManyLinks+9) : Pre-computed effective Henry's Law Coeffs						!!
			!!    (6*HowManyLinks+10)					: Flag for what type of reaction								!!
			!!    (6*HowManyLinks+11)					: "-1" (as a reference for finding HowManyLinks)				!!
			!!    (6*HowManyLinks+12)					: Dissolution Reaction Number									!!
			!!    (6*HowManyLinks+13 - 7*HowManyLinks+12) : Pre-computed effective Henry's Law Coeffs for coupled rxn   !!
			!!    (7*HowManyLinks+13 - 8*HowManyLink+12)  : Flags for HNO3 and HCl dissolution
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! HowManyLinks which the number of aerosol / sections														!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Initialize the Concentrations from the local chemical concentrations: !! 
			Y(1) = GetGridCellChemBurden (XP, YP, ZP, GasChemIndex1) !!!!!!!!!!!!!!!!!!

			IF (DissolutionData(Rxn,11) .EQ. 5) &
			Y(3*HowManyLinks+2) = GetGridCellChemBurden (XP, YP, ZP, GasChemIndex1)

			!! Tell LSODES which chemical we are to consider
			Y(3*HowManyLinks+6) = GasChemIndex1

			!! The condensation rates
			Y(3*HowManyLinks+8:4*HowManyLinks+7) = CondensationRateCoefficients (XP, YP, ZP, Rxn, HowManyLinks)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Reaction Type !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! There are several reaction types in this context, !!
			!! which are different than elsewhere:				 !!!!
			!!  0. Water Condensation								!!
			!!  1. Direct Dissolution of a single chemical			!!
			!!  2. Dissolution of a single species that dissociates	!!
			!!  3. Coupled reactions of types 1 & 2					!!
			!!  4. Dissolution of a species that binds with an ion	!!
			!!  5. Coupled reactions of types 1 & 4	
			!!  7. Sulfate condensation (eq. gas conc. is 0)				!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			Y(6*HowManyLinks+10) = DissolutionData(Rxn,11) 

			IF (SCAFFOLDING) WRITE(*,*) "Dissolution Reaction #", TRIM(INT2STR(Rxn)), "Type: ", TRIM(Real2STR(Y(6*HowManyLinks+10)))
			
			IF(Y(6*HowManyLinks+10) .EQ. 6.) THEN
				IF (Scaffolding) CALL WARN("StepCondensation skips gas-solid reactions.")
				CYCLE
			END IF

			!! A "-1" as reference for later subroutines
			Y(6*HowManyLinks+11) = -1

			!! The reaction number
			Y(6*HowManyLinks+12) = Rxn


			!! Respecify NEQ based on reaction type
			!! And then describe the Jacobian's structure (Ia and Ja, goes into IWORK)
			SELECT CASE(INT(Y(6*HowManyLinks+10)))

			CASE (0:1,8,9)
				NEQ = HowManyLinks+1

				!! Specify Ia and Ja
				IWORK(30+1)     = 1       !! first Ia
				IWORK(30+NEQ+1) = 3*neq-1 !! Final Ia Term

				IWORK(32+NEQ)   = 1		  !! first Ja
						
				DO J = 1, NEQ-1
					IWORK(31+J) = NEQ + (J*2-1)  !! Ia
					IWORK(32+NEQ+J)       = J+1  !! Ja (fill the 1st column)
					IWORK(31+(NEQ+J)*2-1) = 1    !! Ja (fill the 1st row)
					IWORK(31+(NEQ+J)*2)   = J+1  !! Ja (fill the diagonal)
				END DO

			CASE (2,4,7)

				NEQ           = 2*HowManyLinks+1
				IWORK(31+NEQ) = 4*neq-2 !! Final Ia Term

				!! Specify Ia and Ja
				IWORK(31)     = 1       !! first Ia
				IWORK(32+NEQ) = 1	    !! first Ja

				DO J = 2, 2*HowManyLinks+1

					IWORK(30+J)					 = -4 + 2 * HowManyLinks + 3*J      !! Ia
					IWORK(32+2*HowManyLinks+J)   = J    !! Ja (fill the 1st column)
					IWORK(28+4*HowManyLinks+3*J) = 1    !! Ja (fill the 1st row)

					IF (J .LE. HowManyLinks+1) THEN
						IWORK(29+4*HowManyLinks+3*J) = J			   !! Ja (fill the 1st 1st chem's diagonal)
						IWORK(30+4*HowManyLinks+3*J) = J+HowManyLinks  !! Ja (fill the 1st 2nd chem's diagonal)
					ELSE
						IWORK(29+4*HowManyLinks+3*J) = J-HowManyLinks !! Ja (fill the 1st 1st chem's diagonal)
						IWORK(30+4*HowManyLinks+3*J) = J              !! Ja (fill the 1st 2nd chem's diagonal)
					END IF
				END DO

			CASE (3,5)

				!! If this reaction is coupled to a type one reaction, deal with that here
				NEQ           = 3*HowManyLinks+1
				IWORK(31+NEQ) = 11*HowManyLinks + 2 !! Final Ia Term

				!! Specify Ia and Ja
				IWORK(31)     = 1       !! first Ia
				IWORK(32+NEQ) = 1	    !! first Ja

				DO J = 2, 2*HowManyLinks+1
					IWORK(30+J)					 = -4 + 3 * HowManyLinks + 3*J      !! Ia
					IWORK(32+3*HowManyLinks+J)   = J    !! Ja (fill the 1st column)
					IWORK(28+6*HowManyLinks+3*J) = 1    !! Ja (fill the 1st row)

					IF (J .LE. HowManyLinks+1) THEN
						IWORK(29+6*HowManyLinks+3*J) = J			   !! Ja (fill the 1st 1st chem's diagonal)
						IWORK(30+6*HowManyLinks+3*J) = J+HowManyLinks  !! Ja (fill the 1st 2nd chem's diagonal)
					ELSE
						IWORK(29+6*HowManyLinks+3*J) = J-HowManyLinks !! Ja (fill the 1st 1st chem's diagonal)
						IWORK(30+6*HowManyLinks+3*J) = J              !! Ja (fill the 1st 2nd chem's diagonal)
					END IF
				END DO

				DO J = 1, HowManyLinks
					IWORK(31+2*HowManyLinks+J) = 9 * HowManyLinks + J*2  !! Ia
					IWORK(32+NEQ+2*HowManyLinks+J) = J+1+2*HowManyLinks  !! Ja (fill the 1st column)
					IWORK(30 + 12*HowManyLinks+2+2*J) = 1                !! Ja (fill the 1st row)
					IWORK(30 + 12*HowManyLinks+3+2*J)   = 2*HowManyLinks+J+1  !! Ja (fill the diagonal)
				END DO

			END SELECT


			!! Fill the remainder of the Y vector off-site
			CALL DissolutionFactors (XP, YP, ZP, Rxn, HowManyLinks, Temperature, INT(Y(6*HowManyLinks+10)), 	 &
			      Y(2:1+HowManyLinks), Y(HowManyLinks+2:2*HowManyLinks+1), Y(2*HowManyLinks+2:3*HowManyLinks+1), &  ! initial concs
				  Y(HowManyLinks*5+8:HowManyLinks*6+7), Y(HowManyLinks*6+13:HowManyLinks*7+12), ReturnType)
				  			
			!For NH3(g) & H+ => NH4+
			Dum = 0
			IF(INT(Y(6*HowManyLinks+10)) .EQ. 4 .OR. &
			   INT(Y(6*HowManyLinks+10)) .EQ. 5) THEN
					
				DO J = 1, HowManyLinks
					IF(Y(1+J) .LE. 10.0) THEN !If there is less than 10 molecules per particle of NH4+
						Dum = 1
					END IF
				END DO

				IF (Dum == 1) THEN
					IF (SCAFFOLDING) WRITE(*,*) "Skipping the reaction NH3 & H+ => NH4+, NH4+ concentration too low."
					CYCLE
				END IF
			END IF	



			!! Store the concentrations for later
			DO J = 1, NEQ
				InitialConcsArray(J) = Y(J)
			END DO

			!! Check to see if there is any of the chemical in the system
			II = 0. ; JJ = 0.
			DO K = 1, HowManyLinks
				II = II + Y(1+K)
				JJ = JJ + Y(1+HowManyLinks+K)
				KK = KK + Y(1+2*HowManyLinks+K)
			END DO
			
			!! Cycle if not
			SELECT CASE (INT(DissolutionData(Rxn,11)))
			CASE (:1)
				IF (Y(1) .EQ. 0 .AND. II .EQ. 0) CYCLE
			CASE (2)
				IF (Y(1) .EQ. 0 .AND. (II .EQ. 0 .OR. JJ .EQ. 0)) CYCLE
			CASE (3)
				IF (Y(1) .EQ. 0 .AND. (II .EQ. 0 .OR. JJ .EQ. 0) .AND. KK .EQ. 0) CYCLE
			CASE (4)
				IF ((Y(1) .EQ. 0 .OR. JJ .EQ. 0) .AND. II .EQ. 0) CYCLE
			CASE (5)
				IF ((Y(1) .EQ. 0 .OR. JJ .EQ. 0) .AND. II .EQ. 0 .AND. KK .EQ. 0) CYCLE
			END SELECT

		
			!! Each chemical starts as if a new run
			istate = 1

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
30			CALL LSODES (CondensationODEEvaluator,neq,y,t,outtime,itol,rtol,atol,itask, &
						 istate,iopt,rwork,lrw,iwork,liw,CondensationJacobianEvaluator,mf)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF (SCAFFOLDING) WRITE(*,*) "After Condensation LSODES"
			!! If the maximum number of steps occurred, then warn the user (in Scaffolding
			!! mode) and push the integration back into LSODES until the appropriate end time 
			!! is reached.
			IF (ISTATE .EQ. -1 ) THEN !.AND. T < OUTTIME) THEN
				IF (Scaffolding) &
				CALL WARN("Warning! Nominal maximun number of steps in LSODES was exceeded in Condensation Integration for "	&
				 	     //TRIM(GasPhaseChemicalNames(GasChemIndex1))//" ("//trim(int2str(iwork(6)))//" steps)... Recycling")
				ISTATE = 2
				GOTO 30
			END IF

			!! An Error will occur if the array lengths were too small
			IF (IWORK(17) > LRW)  &
				CALL ERROR("The Given Length for LSODES' Working Array RWORK (accessing in StepGasChemistry())"//		&
						  " Should have been of size "//TRIM(INT2STR(iwork(17)))//" but was only "//					&
						  TRIM(INT2STR(LRW))//" instead.  Fix the assignment in EvolveGasChemistry's SetEvolveGasChemConstants().")

			IF (ISTATE .LT. -1) &
			CALL ERROR("In StepCondensation(), LSODES encountered a problem during the integration and returned ISTATE = "//	&
						TRIM(INT2STR(ISTATE))//".  Please investigate.")
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Store the output concentrations until we !!
			!! have considered all of the chemicals		!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! This is indexed as follows:												!!
			!! ChemStorageArray(Rxn, 1)         = Change in Gas Phase Conc 1			!!
			!! ChemStorageArray(Rxn, 2)         = Change in Aerosol Conc 1 (1)			!!
			!!																			!!
			!! ChemStorageArray(Rxn, HML+1)     = Change in Aerosol Conc 1 (HowManyLinks)!
			!! ChemStorageArray(Rxn, HML+2)     = Change in Aerosol Conc 2 (1)			!!
			!!																			!!
			!! ChemStorageArray(Rxn, 2*HML+1)   = Change in Aerosol Conc 2 (HowManyLinks)!
			!! ChemStorageArray(Rxn, 2*HML+2)   = Change in Gas Phase Conc 2			!!
			DO K = 1, NEQ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				ChemStorageArray(Rxn,K) = Y(K) - InitialConcsArray(K)
			END DO

		END DO  !! Ends the loop over all the dissolution reactions
		IF (Scaffolding) WRITE(*,*) "After Dissolution Loop"

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Replace the background chemical concentrations !!

		!! Load an appropriate vector to send to the gas phase gridpoint
		DO I = 1, HowManyDissolutionReactions
			IF (SCAFFOLDING) WRITE(*,*) I, ChemStorageArray(I,1)
			!Skip replace if gas is to be held constant
			IF(INT(DissolutionData(I,1)) .LE. HowManyEvolveGasChems) THEN
				CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(I,1)), ChemStorageArray(I,1))
			END IF
		END DO
		IF (Scaffolding) WRITE(*,*) "After Replace Gas Chems"


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Walk the linked list and replace the chemicals in each particle !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CurrentParticle => FirstParticle

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Replace the chemicals now at the end of the !!
		!! integrations in one pass through the list.  !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO I = 1, HowManyLinks  !! Loop over the aerosol
		DO J = 1, HowManyDissolutionReactions

			!! How we replace depends on the 
			!! reaction type

			IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(DissolutionData(J,3)) + &
															ChemStorageArray(J,I+1) / Avogadro

			CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) + &
															ChemStorageArray(J,I+1) / Avogadro

			IF(CurrentParticle%AqChems(DissolutionData(J,3)) .LT. 0.0) THEN
				!WRITE(*,*) "Negative Cation error, reaction #", J
				!WRITE(*,*) AqPhaseChemicalNames(DissolutionData(J,3)), CurrentParticle%AqChems(DissolutionData(J,3))
			END IF

			RxnType = INT(DissolutionData(J,11)) 
				
			SELECT CASE (RxnType)
				
				CASE (2)
				
				IF(CurrentParticle%AqChems(DissolutionData(J,4)) .LT. 0.0) THEN
					!WRITE(*,*) "Negative Anion error, reaction #", J
					!WRITE(*,*) AqPhaseChemicalNames(DissolutionData(J,4)), CurrentParticle%AqChems(DissolutionData(J,4))
				END IF

				IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(DissolutionData(J,4)) + &
															ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
				CurrentParticle%AqChems(DissolutionData(J,4)) = CurrentParticle%AqChems(DissolutionData(J,4)) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro

			20	IF(CurrentParticle%AqChems(DissolutionData(J,3)) .LT. 0.0 .OR. &
				   CurrentParticle%AqChems(DissolutionData(J,4)) .LT. 0.0) THEN
					!CALL Transcript("In CondensationIntegrator, an acid reaction went below zero.")
					
					!Reset concs to previous value
					!Cation
					CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) - &
															ChemStorageArray(J,I+1) / Avogadro
					!Anion
					CurrentParticle%AqChems(DissolutionData(J,4)) = CurrentParticle%AqChems(DissolutionData(J,4)) - &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
					!Gas
					CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)
					
					Electrolyte = 0.0
					DO C = 1, HowManyDissolutionReactions
						!Find the corresponding gas-solid reaction, if it exists
						IF(DissolutionData(J,1) .EQ. DissolutionData(C,2)) THEN
							Electrolyte = DissolutionData(C,3)
							GasSolidRxn = C
						END IF
					END DO
					
					IF(INT(Electrolyte) .NE. 0. .AND. & !Gas-solid Reaction exists
					   CurrentParticle%AqChems(Electrolyte) .GT. ABS(ChemStorageArray(J,I+1) / Avogadro)) THEN !There is enough salt
						!WRITE(*,*) "Salt Evaporation!"
		
						!Evaporate Electrolyte (remember chemstoragearray is negative here!)
						CurrentParticle%AqChems(Electrolyte) = CurrentParticle%AqChems(Electrolyte) + ChemStorageArray(J,I+1) / Avogadro
						
						!Add acid gas to gas phase
						CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), -1.0*ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)
						
						!Add ammonia to gas phase
						CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(GasSolidRxn,1)), -1.0*ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)
						
						!Equilibrate gas-phase ammonia (to prevent oscillations)
						DO C = 1, HowManyDissolutionReactions
							IF (GasPhaseChemicalNames(DissolutionData(C,1)) .EQ. "NH3" &
								.AND. INT(DissolutionData(C,11)) .EQ.  4 &
								.AND. AmmoniaFlux .EQ. .FALSE. ) THEN
								IF (SCAFFOLDING) WRITE(*,*) "Before Ammonia Equilibrate, particle #", J
				
								IF (SCAFFOLDING) WRITE(*,*) "Equilibrating Ammonia!"
								CALL KusikMeissner (CurrentParticle)
								XX = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
								Temperature = GetTemp(XX(1),XX(2),XX(3))
								GasPhaseBurden = GetGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(C,1)))
                                                                WaterSatBurden = GetSatVapBurden (XX(1), XX(2), XX(3))
								CALL EquilibrateDissolutionReaction (GasPhaseBurden, C, CurrentParticle, Temperature, ReturnType, WaterSatBurden = WaterSatBurden, GasPhaseBurden2=0.0, OptErrorTolerance=AqThermoEquilibriumError/2.)
								CALL ReplaceGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(C,1)), GasPhaseBurden)		
							END IF
						END DO

						!Cycle (should just skip rest of check)
						GOTO 20
					END IF
						
					!If ammonia salt doesn't fix it,
					!decide if anion or cation is limiting
					Store3 = CurrentParticle%AqChems(DissolutionData(J,3))
					Store4 = CurrentParticle%AqChems(DissolutionData(J,4))

					IF(Store3 .LE. Store4) THEN !Cation limiting
										
						!WRITE(*,*) "Cation limiting! This is always H+"
						IF (Store3 .NE. 0.0) THEN
							!WRITE(*,*) "Store3 not zero!"
							!Only allow ions to reach 0, not go under
							!Cation
							!WRITE(*,*) "Cation: ", AqPhaseChemicalNames(DissolutionData(J,3))
							CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) - &
															0.9*Store3
							!WRITE(*,*) "Corrected Cation conc: ", CurrentParticle%AqChems(DissolutionData(J,3))
							!This keeps H+ above zero to prevent errors with ammonia
							!Anion
							!WRITE(*,*) "Anion: ", AqPhaseChemicalNames(DissolutionData(J,4))
							CurrentParticle%AqChems(DissolutionData(J,4)) = CurrentParticle%AqChems(DissolutionData(J,4)) - &
																0.9*Store3
							!Acid(Gas)
							!WRITE(*,*) "Gas: ", GasPhaseChemicalNames(INT(DissolutionData(J,1)))
							CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), 0.9*Store3*Avogadro*CurrentParticle%NumberofParticles)
						
						ELSE IF(Store3 .EQ. 0.0 .AND. Store4 .NE. 0.0) THEN
							!WRITE(*,*) "Store3 zero!"
							!Push cation (H+) slightly above 0
							!Cation
							!WRITE(*,*) "Skippy Cation: ", AqPhaseChemicalNames(DissolutionData(J,3))
							CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) + &
															0.01*Store4
							!WRITE(*,*) "Corrected Cation conc: ", CurrentParticle%AqChems(DissolutionData(J,3))
							!This keeps H+ above zero to prevent errors with ammonia
							!Anion
							!WRITE(*,*) "Anion: ", AqPhaseChemicalNames(DissolutionData(J,4))
							CurrentParticle%AqChems(DissolutionData(J,4)) = CurrentParticle%AqChems(DissolutionData(J,4)) + &
																0.01*Store4
							!Acid(Gas)
							!WRITE(*,*) "Gas: ", GasPhaseChemicalNames(INT(DissolutionData(J,1)))
							CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), -0.01*Store4*Avogadro*CurrentParticle%NumberofParticles)
							
						END IF

					ELSE !Anion Limiting

							WRITE(*,*) "Anion Limiting!"
							WRITE(*,*) "Crude Correction!"
							!Cation
							CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) - &
															Store4
							!Anion
							CurrentParticle%AqChems(DissolutionData(J,4)) = CurrentParticle%AqChems(DissolutionData(J,4)) - &
																Store4
							!Acid(Gas)
							CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), Store4*Avogadro*CurrentParticle%NumberofParticles)
						END IF
					END IF
					
				
				CASE (3) !No zero check!
				CurrentParticle%AqChems(DissolutionData(J,4)) = CurrentParticle%AqChems(DissolutionData(J,4)) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro


				CurrentParticle%AqChems(DissolutionData(DissolutionData(J,15),3)) =								&
															CurrentParticle%AqChems(DissolutionData(DissolutionData(J,15),3)) + &
															ChemStorageArray(J,I+2*HowManyLinks+1) / Avogadro

				CASE (4)
				IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(DissolutionData(J,2)) + ChemStorageArray(J,I+HowManyLinks+1) / Avogadro 
				
				CurrentParticle%AqChems(DissolutionData(J,2)) = CurrentParticle%AqChems(DissolutionData(J,2)) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro

				IF (CurrentParticle%AqChems(DissolutionData(J,2)) .LT. 0) THEN
					
					!CALL Transcript("In CondensationIntegrator, there was a cation concentration that went below zero.")
					!Reset concs to previous value
					!NH4+
					CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) - &
															ChemStorageArray(J,I+1) / Avogadro
					!H+
					CurrentParticle%AqChems(DissolutionData(J,2)) = CurrentParticle%AqChems(DissolutionData(J,2)) - &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
					!NH3(Gas)
					CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)

					!Only allow H+ to reach 0, not go under
					!Store previous H+
					Store = CurrentParticle%AqChems(DissolutionData(J,2))
					!NH4+
					CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) + &
															Store
					!H+
					CurrentParticle%AqChems(DissolutionData(J,2)) = CurrentParticle%AqChems(DissolutionData(J,2)) - &
																Store
					!NH3(Gas)
					CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), -1.*Store*Avogadro*CurrentParticle%NumberofParticles)

				END IF
				
				CASE (5)
				IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(DissolutionData(J,2)) +  ChemStorageArray(J,I+HowManyLinks+1) / Avogadro 
				CurrentParticle%AqChems(DissolutionData(J,2)) = CurrentParticle%AqChems(DissolutionData(J,2)) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
		
				IF (CurrentParticle%AqChems(DissolutionData(J,2)) .LT. 0) THEN
					!Reset concs to previous value
					!NH4+
					CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) - &
															ChemStorageArray(J,I+1) / Avogadro
					!H+
					CurrentParticle%AqChems(DissolutionData(J,2)) = CurrentParticle%AqChems(DissolutionData(J,2)) - &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
					!NH3(Gas)
					CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)

					!Only allow H+ to reach 0, not go under
					!NH4+
					CurrentParticle%AqChems(DissolutionData(J,3)) = CurrentParticle%AqChems(DissolutionData(J,3)) + &
															CurrentParticle%AqChems(DissolutionData(J,2))
					!H+
					CurrentParticle%AqChems(DissolutionData(J,2)) = CurrentParticle%AqChems(DissolutionData(J,2)) - &
																CurrentParticle%AqChems(DissolutionData(J,2))
					!NH3(Gas)
					CALL AddToGridCellChemBurden (XP, YP, ZP, INT(DissolutionData(J,1)), -1.*CurrentParticle%AqChems(DissolutionData(J,2))*Avogadro*CurrentParticle%NumberofParticles)

				END IF
				
				IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(DissolutionData(DissolutionData(J,15),3)) + ChemStorageArray(J,I+2*HowManyLinks+1) / Avogadro 
				CurrentParticle%AqChems(DissolutionData(DissolutionData(J,15),3)) =								&
															CurrentParticle%AqChems(DissolutionData(DissolutionData(J,15),3)) + &
															ChemStorageArray(J,I+2*HowManyLinks+1) / Avogadro
				
				CASE (7)
				CurrentParticle%AqChems(DissolutionData(J,4)) = CurrentParticle%AqChems(DissolutionData(J,4)) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
			
			END SELECT
		END DO 
		
		!! Re-Equilibrate
		!IF (Scaffolding) 
		!WRITE(*,*) "Before Electrolyte Eq"				
		!WRITE(*,*) "H+: ", CurrentParticle%AqChems(ProtonIndex)
		CALL FindElectrolyteEquilibrium(CurrentParticle, UpdateThermo = .TRUE., FirstEquilibration=.FALSE., ReturnType=GG)
		!IF (Scaffolding) 
		!WRITE(*,*) "After Electrolyte Eq"
		!PAUSE

		CurrentParticle => CurrentParticle%Next

		END DO !! Quit loop over particles for replacement

		END DO !! Stop Looping Over the Time Steps

		!! Get rid of everything that was allocated in here.
		DEALLOCATE (Y,RWORK,IWORK,InitialConcsArray,ChemStorageArray)

		!! SCSCSCSC
		IF (Scaffolding) CALL TRANSCRIPT (">>>Exiting StepCondensation<<<")
		!WRITE(*,*) "StepCondensation took "//trim(real2str(CondTime))//" seconds."
		
	RETURN
END SUBROUTINE StepCondensation

!! Matt wrote this subroutine to calculate the correct saturation concentration
!! of HNO3, HCl, and NH3 when the dissolution of these species is coupled to gas-solid equilibrium reactions ( (c) 2006)
SUBROUTINE GasSatConcAll(Y, ParticleIndex, NumBins, InParticle, HClRxn, HNO3Rxn, NH3Rxn, Temperature, &
						 SkipCl, SkipNO3, &
                         HClSoverH, HNO3SoverH, NH3SoverH, HClSat, HClFlag, HNO3Sat, HNO3Flag, NH3Sat, NH3Flag, &
						 HClJacTerm, HNO3JacTerm, NH3JacTerm)

	USE GridPointFields, ONLY :		GasPhaseChemicalRates,			&
									GetSatVapBurden,				&
									LocateGridPoint,				&
									GasPhaseChemGridIndex,			&
									GetGridCellChemBurden, &
									GridGasChem, ClosestGridPoint, &
									GetRelativeHumidity

	USE Aerosols,        ONLY :		Particles, Particle,			&
									ParticleArray, Molality, EqCoeff, &
									DeliquesenceRH, EquilibriumConstantsRatio, &
									AerosolPh

	USE ModelParameters, ONLY :		RstarMB, moles, grams, watermolecmass,	&
									ProtonIndex, HydroxyIndex,				&
									Avogadro, Metastable, MinimumWater

	USE Chemistry,       ONLY :		HowManyAqChems,HowManyAqCations,		&
									HowManyAqAnions,						&
									AqEquilibriaList,                       &
									HowManyAqEqReactions,						&
									GasPhaseChemicalNames,					&
									AqPhaseChemicalNames,					&
									AqCationNames,AqAnionNames, & 
									FindChem

	USE InfrastructuralCode, ONLY : ERROR, WARN, TRANSCRIPT, Int2Str, Real2Str

	IMPLICIT NONE

	!! External Variables, in
	TYPE(Particle), POINTER :: InParticle
	INTEGER :: HClRxn, HNO3Rxn, NH3Rxn, ParticleIndex, NumBins
	REAL*8  :: Temperature, HClSoverH, HNO3SoverH, NH3SoverH, Y
	DIMENSION :: Y(1) !gases in molecules/cm3, aerosol in molecules/particle
	LOGICAL :: SkipCl, SkipNO3

	!! External Variables, out
	REAL*8, INTENT(OUT) ::  HClSat, HNO3Sat, NH3Sat
	INTEGER, INTENT(OUT) :: HClFlag, HNO3Flag, NH3Flag
    REAL*8, INTENT(OUT) ::  HClJacTerm, HNO3JacTerm, NH3JacTerm

	!! Internal Variables
	INTEGER :: I, J, NH4NO3Index, NH4ClIndex, NH4NO3GasIndex, NH4ClGasIndex
	REAL*8 :: TempSatConc, RxnSatConc, Keq_NH4NO3, Keq_NH4Cl, CatConc, WaterFactor, & 
	           P_NH3, TempGasSatConc, AnionConc, ProtonConc, NH4Conc, NH3GasConc, &
			   DRH_NH4NO3, DRH_NH4Cl, RH, XP, YP, ZP, P_o, P_Acid, P_NH3Sat, P_AcidSat, &
			   HNO3GasConc, HClGasConc, NO3Conc, ClConc, TempHClSat, TempHNO3Sat, &
			   TempNH3Sat, P_HNO3, P_HCl, P_HNO3Sat, P_HClSat, Dummy
	LOGICAL :: UseNH310
	REAL*8, PARAMETER :: Tiny = 1.0e-10 

	!Check to make sure reaction data is correct
	IF (GasPhaseChemicalNames(DissolutionData(HNO3Rxn,1)) .NE. "HNO3") THEN
		CALL ERROR("In GasSatConcAll (in CondensationIntegrator.h) you have the incorrect index &
		            for the dissolution of HNO3.")
	ELSE IF(GasPhaseChemicalNames(DissolutionData(HClRxn,1)) .NE. "HCl") THEN
		CALL ERROR("In GasSatConcAll (in CondensationIntegrator.h) you have the incorrect index &
		            for the dissolution of HCl.")
	ELSE IF(GasPhaseChemicalNames(DissolutionData(NH3Rxn,1)) .NE. "NH3") THEN
		CALL ERROR("In GasSatConcAll (in CondensationIntegrator.h) you have the incorrect index &
		            for the dissolution of NH3.")
	END IF

	ProtonConc = Y(4 + 5*NumBins + ParticleIndex) !molecules/particle
	IF(ProtonConc .LT. 0.) THEN
		WRITE(*,*) "ProtonConc (molec/part): ", ProtonConc
		CALL Error("Proton Conc below 0 in GasSatConcAll.")
	END IF

	
	HNO3GasConc = Y(2) !molecules/cm3
	HClGasConc = Y(3) !molecules/cm3
	NH3GasConc = Y(4) !molecules/cm3
	
	ProtonConc = Y(4 + 5*NumBins + ParticleIndex) !molecules/particle
	NH4Conc	   = Y(4 + 3*NumBins + ParticleIndex) !molecules/particle
	NO3Conc    = Y(4 + NumBins + ParticleIndex) !molecules/particle
	ClConc     = Y(4 + 2*NumBins + ParticleIndex) !molecules/particle

	XP = Y(5+6*NumBins)
	YP = Y(6+6*NumBins)
	ZP = Y(7+7*NumBins)
	RH = GetRelativeHumidity (XP, YP, ZP)

	!! NOTE: Flag = 10 is normal dissolution
	!!       Flag = 30 is when dissolution is coupled to a gas-solid reaction, with NH3(g) term
	!!       Flag = 31 is when dissolution is coupled to a gas-solid reaction, with HNO3(g) term
	!!       Flag = 32 is when dissolution is coupled to a gas-solid reaction, with HCl(g) term

	!Assume Flag = 10 to start
	HNO3Flag = 10
	HClFlag = 10
	NH3Flag = 10
	HClJacTerm = 0.
	HNO3JacTerm = 0.
	NH3JacTerm = 0.

	!Can't couple to solids if metastable!
	!RH check guarantees flag 10 for NH3 below 20% RH
	IF(Metastable .OR. RH .LT. 0.20) THEN
		HNO3Sat = HNO3SoverH*NO3Conc*ProtonConc
		IF (HNO3Sat .LT. 0.) CALL ERROR("Negative HNO3 Sat. Conc. in GasSatConcAll.")
		HClSat = HClSoverH*ClConc*ProtonConc
		IF (HClSat .LT. 0.) CALL ERROR("Negative HCl Sat. Conc. in GasSatConcAll.")
		NH3Sat = NH3SoverH*NH4Conc/(ProtonConc+Tiny)
		IF (NH3Sat .LT. 0.) CALL ERROR("Negative NH3 Sat. Conc. in GasSatConcAll.")		
		RETURN
	END IF
	
	!Calculate gas-sat conc for coupled rxn
	!WRITE(*,*) "Gas Conc Calc: ", SoverH, TempSatConc, ProtonConc 
	TempHNO3Sat = HNO3SoverH*NO3Conc*ProtonConc
	TempHClSat = HClSoverH*ClConc*ProtonConc
	TempNH3Sat = NH3SoverH*NH4Conc/(ProtonConc+Tiny)
	
	P_NH3 = NH3GasConc*RstarMB*Temperature/Avogadro
	P_HNO3 = HNO3GasConc*RstarMB*Temperature/Avogadro
	P_HCl = HClGasConc*RstarMB*Temperature/Avogadro

	!Check for coupling to gas-solid reactions
	DO I = 1, HowManyDissolutionReactions
	
		!Check for gas-solid reactions 
		IF(DissolutionData(I,11) .EQ. 6) THEN
			
			IF(GasPhaseChemicalNames(DissolutionData(I,2)) .EQ. "HNO3") THEN
				NH4NO3GasIndex = I
				Keq_NH4NO3 = HenrysLaw (Temperature, I)
				!Check to make sure RH is below DRH for electrolyte
				DO J = 1, HowManyAqEqReactions
					IF(DissolutionData(I,3) .EQ. AqEquilibriaList(J,1)) &
							DRH_NH4NO3 = DeliquesenceRH (J, InParticle)
							
				END DO
			ELSE IF(GasPhaseChemicalNames(DissolutionData(I,2)) .EQ. "HCl") THEN
				NH4ClGasIndex = I
				Keq_NH4Cl = HenrysLaw (Temperature, I)
				!Check to make sure RH is below DRH for electrolyte
				DO J = 1, HowManyAqEqReactions
					IF(DissolutionData(I,3) .EQ. AqEquilibriaList(J,1)) &
							DRH_NH4Cl = DeliquesenceRH (J, InParticle)
				END DO
			ELSE
				CALL ERROR("In GasSatConcAll, found a bad gas-solid reaction.")
			END IF
				
		END IF

	END DO

	!Find NH4NO3 and NH4Cl aqueous reactions
	DO I = 1, HowManyAqEqReactions
		IF(AqEquilibriaList(I,1) .GT. 0) THEN 
		IF(AqPhaseChemicalNames(AqEquilibriaList(I,1)) .EQ. "NH4NO3") THEN
			NH4NO3Index = I
		ELSE IF (AqPhaseChemicalNames(AqEquilibriaList(I,1)) .EQ. "NH4Cl") THEN
			NH4ClIndex = I
		END IF
		END IF
	END DO		
		
	!If above DRH for both, just do metastable
	IF(DRH_NH4NO3 .LE. RH .AND. DRH_NH4Cl .LE. RH) THEN
		
		HNO3Sat = HNO3SoverH*NO3Conc*ProtonConc
		IF (HNO3Sat .LT. 0.) CALL ERROR("Negative HNO3 Sat. Conc. in GasSatConcAll.")
		HClSat = HClSoverH*ClConc*ProtonConc
		IF (HClSat .LT. 0.) CALL ERROR("Negative HCl Sat. Conc. in GasSatConcAll.")
		NH3Sat = NH3SoverH*NH4Conc/(ProtonConc+Tiny)
		IF (NH3Sat .LT. 0.) CALL ERROR("Negative NH3 Sat. Conc. in GasSatConcAll.")
		RETURN
	
	!If skipping both, do metastable for NH3, dummy for others
	ELSE IF(SkipCl .AND. SkipNO3) THEN
		HClSat = 1.0
		HNO3Sat = 1.0
		NH3Sat = NH3SoverH*NH4Conc/(ProtonConc+Tiny)
		IF (NH3Sat .LT. 0.) CALL ERROR("Negative NH3 Sat. Conc. in GasSatConcAll.")
		RETURN

	!First case, do all gases
	ELSE IF(DRH_NH4NO3 .GT. RH .AND. DRH_NH4Cl .GT. RH .AND. NOT(SkipCl) .AND. NOT(SkipNO3)) THEN
	
		!Idea from Jacobson, Sec 17.14, p. 594
		P_o = P_NH3 - P_HCl - P_HNO3
		P_NH3Sat = 0.5*P_o +0.5*SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4NO3 + 1.0/Keq_NH4Cl))) 
		NH3Sat = P_NH3Sat*Avogadro/RstarMB/Temperature
		
		!Ad hoc fix for low NH3 conditions (MJA, 09/21/07)
		IF(NH3Sat .GT. Y(4)) THEN
			!WRITE(*,*) "Hi Matt!"
			UseNH310 = .TRUE.
			P_NH3Sat = Y(4)*RstarMB*Temperature/Avogadro
		END IF
		
		P_HNO3Sat = 1.0/(Keq_NH4NO3*P_NH3Sat+Tiny)
		P_HClSat = 1.0/(Keq_NH4Cl*P_NH3Sat+Tiny)
		!NH3Sat = P_NH3Sat*Avogadro/RstarMB/Temperature
		HNO3Sat = P_HNO3Sat*Avogadro/RstarMB/Temperature
		HClSat = P_HClSat*Avogadro/RstarMB/Temperature		
		
		!If Gas-Solid reaction give lower gas sat conc, use that instead
		!Only use if RH is below DRH of electrolyte
		!WRITE(*,*) "HNO3Sat (ppbv): ", P_HNO3Sat*1.0e9/1000.0
		IF(HNO3Sat .LT. TempHNO3Sat ) THEN !&
			!.OR. InParticle%AqChems(FindChem("NH4NO3", 1)) .GT. 1.0e-24 ) THEN
			HNO3Flag = 33
			HNO3JacTerm = (0.5/(Keq_NH4NO3*P_NH3Sat*P_NH3Sat+Tiny)) &
			              *(1.0 + P_o/(Tiny+SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4NO3 + 1.0/Keq_NH4Cl)))))
		ELSE
			HNO3Sat = TempHNO3Sat
		END IF
	
		!If Gas-Solid reaction give lower gas sat conc, use that instead
		!Only use if RH is below DRH of electrolyte
		!WRITE(*,*) "Test: ", HClSat, TempHClSat
		IF(HClSat .LT. TempHClSat ) THEN !.OR. InParticle%AqChems(FindChem("NH4Cl", 1)) .GT. 1.0e-24) THEN
			HClFlag = 33
			HClJacTerm = (0.5/(Keq_NH4Cl*P_NH3Sat*P_NH3Sat+Tiny)) &
			              *(1.0 + P_o/(Tiny+SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4NO3 + 1.0/Keq_NH4Cl)))))
		ELSE
			HClSat = TempHClSat
		END IF
	
		!If Gas-Solid reaction give lower gas sat conc, use that instead
		!Only use if RH is below DRH of electrolyte
		IF(NH3Sat .LT. TempNH3Sat ) THEN !.AND. NOT(UseNH310) ) THEN
			NH3Flag = 33
			NH3JacTerm = (0.5) &
			              *(1.0 + P_o/(Tiny+SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4NO3 + 1.0/Keq_NH4Cl)))))
		ELSE
			NH3Sat = TempNH3Sat
		END IF
	
		!WRITE(*,*) "HCl: ", P_HClSat*1.0e6, HClFlag
		!WRITE(*,*) "HNO3: ", P_HNO3Sat*1.0e6, HNO3Flag
		!WRITE(*,*) "NH3: ", P_NH3Sat*1.0e6, NH3Flag

		IF (HNO3Sat .LT. 0. ) THEN
			CALL WARN("Negative HNO3 Sat. Conc. in GasSatConcAll.")
		END IF
		IF (HClSat .LT. 0. ) THEN
			CALL WARN("Negative HCl Sat. Conc. in GasSatConcAll.")
		END IF
		IF (NH3Sat .LT. 0. ) THEN
			CALL WARN("Negative NH3 Sat. Conc. in GasSatConcAll.")
		END IF
		RETURN
	
	!Second case, just do NH4NO3
	ELSE IF(SkipCl .OR. DRH_NH4Cl .LE. RH) THEN
		
		!Idea from Jacobson, Sec 17.14, p. 594
		P_o = P_NH3 - P_HNO3
		P_NH3Sat = 0.5*P_o +0.5*SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4NO3))) 
		P_HNO3Sat = 1.0/(Keq_NH4NO3*P_NH3Sat+Tiny)
		NH3Sat = P_NH3Sat*Avogadro/RstarMB/Temperature
		HNO3Sat = P_HNO3Sat*Avogadro/RstarMB/Temperature


		!If Gas-Solid reaction give lower gas sat conc, use that instead
		!Only use if RH is below DRH of electrolyte
		IF(HNO3Sat .LT. TempHNO3Sat) THEN
			HNO3Flag = 34
		    HNO3JacTerm = (0.5/(Keq_NH4NO3*P_NH3Sat*P_NH3Sat+Tiny)) &
			              *(1.0 + P_o/(Tiny+SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4NO3)))))
	
		ELSE
			HNO3Sat = TempHNO3Sat
		END IF
	
		IF(SkipCl) THEN
			HClFlag = 10
			HClSat = 1.0
		ELSE
			HClFlag = 10
			HClSat = TempHClSat
		END IF
	
		!If Gas-Solid reaction give lower gas sat conc, use that instead
		!Only use if RH is below DRH of electrolyte
		IF(NH3Sat .LT. TempNH3Sat) THEN
			NH3Flag = 34
		    NH3JacTerm = (0.5) &
			              *(1.0 + P_o/(Tiny+SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4NO3)))))
	
		ELSE
			NH3Sat = TempNH3Sat
		END IF
	
		IF (HNO3Sat .LT. 0. ) THEN
			CALL WARN("Negative HNO3 Sat. Conc. in GasSatConcAll.")
		END IF
		IF (HClSat .LT. 0. ) THEN
			CALL WARN("Negative HCl Sat. Conc. in GasSatConcAll.")
		END IF
		IF (NH3Sat .LT. 0. ) THEN
			CALL WARN("Negative NH3 Sat. Conc. in GasSatConcAll.")
		END IF
		RETURN
	
	!Third case, just do NH4Cl
	ELSE IF(SkipNO3 .OR. DRH_NH4NO3 .LE. RH) THEN
		
		!Idea from Jacobson, Sec 17.14, p. 594
		P_o = P_NH3 - P_HCl
		P_NH3Sat = 0.5*P_o +0.5*SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4Cl))) 
		P_HClSat = 1.0/(Keq_NH4Cl*P_NH3Sat+Tiny)
		NH3Sat = P_NH3Sat*Avogadro/RstarMB/Temperature
		HClSat = P_HClSat*Avogadro/RstarMB/Temperature
				
		IF(SkipNO3) THEN
			HNO3Flag = 10
			HNO3Sat = 1.0
		ELSE
			HNO3Flag = 10
			HNO3Sat = TempHNO3Sat
		END IF
		
		!If Gas-Solid reaction give lower gas sat conc, use that instead
		!Only use if RH is below DRH of electrolyte
		IF(HClSat .LT. TempHClSat) THEN
			HClFlag = 35
		    HNO3JacTerm = (0.5/(Keq_NH4Cl*P_NH3Sat*P_NH3Sat+Tiny)) &
			              *(1.0 + P_o/(Tiny+SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4Cl)))))

		ELSE
			HClSat = TempHClSat
		END IF
	
		!If Gas-Solid reaction give lower gas sat conc, use that instead
		!Only use if RH is below DRH of electrolyte
		IF(NH3Sat .LT. TempNH3Sat) THEN
			NH3Flag = 35
		    NH3JacTerm = (0.5) &
			              *(1.0 + P_o/(Tiny+SQRT(ABS((P_o**2.0) + 4.0*(1.0/Keq_NH4Cl)))))

		
		ELSE
			NH3Sat = TempNH3Sat
		END IF

		IF (HNO3Sat .LT. 0. ) THEN
			CALL WARN("Negative HNO3 Sat. Conc. in GasSatConcAll.")
		END IF
		IF (HClSat .LT. 0. ) THEN
			CALL WARN("Negative HCl Sat. Conc. in GasSatConcAll.")
		END IF
		IF (NH3Sat .LT. 0. ) THEN
			CALL WARN("Negative NH3 Sat. Conc. in GasSatConcAll.")
		END IF
		RETURN

	END IF

END SUBROUTINE GasSatConcAll
