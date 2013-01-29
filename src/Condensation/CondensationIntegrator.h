!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! CondensationIntegrator.h
!! This is the numerical integration routine for the flux-limited
!! kinetic dissolution of inorganic compounds into the aerosol aqueous phase.
!! As the main routine calls LSODES, the required ODE and Jacobian 
!! functions are given here as well
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY					      	             !!
!!									     !!
!! Month  Year   Name              Description			             !!
!! 07     2006   Matt Alvarado     Began Update History		             !!
!! 08/28  2006   Matt Alvarado     Commented out some write statements       !!
!! 08/30  2006   Matt Alvarado     1. Edited StepCondensationAll to call     !!
!!						RegridAerosol at end.	     !!
!!				   2. Edited StepCondensation to exit	     !!
!!					    water eq. loop after 5 iterations!!
!!				  3. Removed commented code from	     !!
!!					       StepCondensation		     !!
!!				  4. Fixed Jacobian for NO3 and Cl	     !!
!! 09/08  2006   Matt Alvarado     Changed water eq. loop in StepCondensation!!
!! 09/13  2006   Matt Alvarado     Changed Ammoniaeq loop in StepCondensation!!
!! 09/18  2006   Matt Alvarado     Added eq check to StepCondensationAll     !!
!! 10/03  2006   Matt Alvarado     Updated GasSatConc and NH3SatConc	     !!
!! 10/13  2006   Matt Alvarado     Created BigCondensationOdeEvaluator and   !!
!!				     	BigCondensationJacobianEvaluator     !!
!! 10/17  2006   Matt Alvarado     Moved EmptyGridCell check and aerosol 
!!                                      sorting to 
!!				      	StepCondensationAll()
!! 07/12  2007   Matt Alvarado     Have inorganic condensation other than
!!			NH3 and H2SO4 shut down when RH is less than 30%
!!			To remove oscillations, exit water eq loop
!!			in StepCondensation after 50 iterations if 
!!			change is less than 1%.
!! 07/13  2007   Matt Alvarado Remove coupling to solid salt in AcidGasSatConc 
!!				and NH3SatConc when metastable
!!				or when above DRH for salt.
!! 07/16  2007   Matt Alvarado  Exit water equilibrium in StepCondensation 
!!				after 1000 iterations to prevent oscilating
!!				nonconvergence at low RH.
!!				Try forcing NH3 condensation to be flag 10
!!				in NH3SatConc
!!			Try forcing HCl, HNO3 condensation to be flag 10
!!				in AcidGasSatConc
!!				Add flag to catch when Cl- and NO3- are near 
!!			zero and evaporating to BigCondensationODEEvaluator
!!			and BigCondensationJacobianEvaluator
!!			Set HClCoeff(I) to 0.0 for all particles
!!			when Cl is not present (skips this condensation!)
!! 07/20  2007   Matt Alvarado Added flag to skip calculation of AcidGasSatConc
!!			in BigCondODE and BigConJac if RH is less than 30%
!! 07/30  2007   Matt Alvarado    Put in all flags for condensation, 
!!                               except skip NH3 flag 22
!!				(NH42SO4 formation) in NH3SatConc
!!				Increased atol to reduce charge balance errors
!! 07/31  2007   Matt Alvarado     Put code in BigCondensationODEEvaluator 
!!                 to catch NH4+ below 0 and recover by dissociating NH4NO3
!!		   or NH4Cl 
!! 08/01  2007   Matt Alvarado  Added code to just exit gas-solid eq at end of 
!!		SUBROUTINE BigStepInorgCondensation after max iterations
!!			Added NH42SO4 and NH4HSO4 to NH4 below zero flag
!! 08/06  2007   Matt Alvarado     Forced NH3SatConc to use gas-solid 
!!   equilibrium when coupled to NO3- or Cl-.
!! 08/07  2007   Matt Alvarado     Added line to NH3SatConc to force Flag 10 if
!!			RH is below 30% (and thus we are doing gas-solid eq
!!			for Cl- and NO3-)
!!		Fixed gas-solid eq loop at end of BigStepInorgCondensation
!! 08/08  2007   Matt Alvarado      Removed forcing of gas-solid in NH3SatConc,
!!     because this caused the 2D simulations to crash
!! 08/13  2007   Matt Alvarado      Fixed sign error in 
!!                                  BigCondensationJacobianEvaluator 
!!                                  (Case 7, HClFlag=10, d/dCl-(dH+/dt) 
!!			Fixed bug in Gas-Solid loop in BigStepInorgCondensation
!! 08/17  2007   Matt Alvarado      Try Jacobson gas-solid idea in GasSatConcs 
!!                                  (See 17.14, p. 594)
!!		Removed all 20 flags from NH3SatConc and AcidGasSatConc
!!		Created GasSatConcAll
!!		Fixed below 0 checks
!! 08/20  2007   Matt Alvarado      Fixed GasSatConcAll
!!		Added Cl- below 0 check to BigStepInorgCondensation
!! 08/21  2007   Matt Alvarado Redid below 0 checks in BigStepInorgCondensation
!!			Set flag to limit HClSat and HNO3Sat to total HCl
!!			and HNO3 available.
!! 08/22  2007   Matt Alvarado    Added a second H+ below 0 check to 
!!                                BigCondensationODEEvaluator
!!		to handle the case where both H+ and NH4+ are below 0
!! 09/01  2007   Matt Alvarado    Forced below 0 H+ recovery to use bigger 
!!                                of NH4+ and H2O eq estimates
!! 09/03  2007   Matt Alvarado    Removed NH4+ from Cl- and NO3- below 
!!                                zero loops in LSODES
!!				Hopefully will fix charge balance errors.
!! 09/04  2007   Matt Alvarado      Made flag for low RH set to 20%
!!			Added in new Jacobian for gas-solid rxns
!!		Changes to GasSatConcAll, BigCondensationJacobianEvaluator, 
!!			and Ia and Ja in BigStepInorgCondensation
!!                 Added Na, Mg and Ca Cl and NO3 salts to Clsum and NO3sum
!!		in BigODE and BigJac.
!!		Added Na, Mg and Ca salts to below 0 check at end of
!!		BigStepInorgCondensation
!! 09/05  2007   Matt Alvarado      Removed forcing of gas-solid eq when 
!!                                   NH4Cl or NH4NO3 are present
!! 09/18  2007   Matt Alvarado      Forced NO3NearZero to be true when SkipNO3 
!!                                  is true,as was already done for Cl
!! 09/21  2007   Matt Alvarado      Cleaned up Jacobian and ODE, removed 
!!                                  AcidGasSatConc and NH3SatGasConc, and 
!!                                  DissociateAmmoniaSalt
!!				Added Tiny to Jacobain and GasSatConcAll calcs
!! 09/27  2007   Matt Alvarado      Removed 1.0e-4 molec/particle and
!!                                  evap check for Cl- and NO3-
!! 10/01  2007   Matt Alvarado      Added EquilibrateInorgGridPoint option 
!!                                  to StepCondensation
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/17  2012   Matt Alvarado     Removed out-of-date and unused subroutines!!
!!                                 GasSatConcAll,                            !!
!!                                  BigCondensationJacobianEvaluator,        !!
!!                                 BigCondensationODEEvaluator, and          !!
!!                                 BigStepInorgCondensation.                 !!
!! 01/24  2013   Matt Alvarado     Removed CalculateSurrogateGasConc and     !!
!!                                  UpdateGasPhaseOrganics from              !!
!!                                  StepCondensationAll                      !!
!! 01/25  2013   Matt Alvarado    Added UNIFAC cals to StepCondensationAll   !!
!! 01/28  2013   Matt Alvarado    Removed call to EquilibriumWaterContentAll !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE StepCondensationAll ()					!!
!! 2. SUBROUTINE DissolutionFactors					!!
!! 3. SUBROUTINE StepCondensation ()                                    !!
!! 4. SUBROUTINE CondensationJacobianEvaluator                          !!
!! 5. SUBROUTINE CondensationODEEvaluator                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This allows the main program to call all the condensation integrators 
!! together. Also regrids aerosol automatically.
SUBROUTINE StepCondensationAll (TimestepSize, NumTimeSteps)

	USE InfrastructuralCode
	USE ModelParameters
	USE Aerosols, ONLY :  Particles, Particle,			&
			      ParticleArray,				&
			      RegridAerosol,				&
			      SortAerosolAtGridPointForCondensation
						  		
	IMPLICIT NONE

	!! External Variables
	REAL*8  :: TimeStepSize		!! Should be in Seconds
	INTEGER :: NumTimeSteps

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
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! If there are no aerosol at the gridcell (or if they are all hydrophobic),  !
!! then there is no need to condense                                          !
	EmptyGridCell = .TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CurrentParticle => Particles%First

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
			CALL StepCondensation (TimestepSize, 1)
	
		!If doing flux
		ELSE
			IF (ThermoBulkMode) &
			CALL ERROR("Can't try to use the kinetic based condensation and dissolution routine when running in bulk mode "// &
				   "(specified in AerosolModes.in).  Please remedy.")

			!! Sort the gridcell, pushing all empty sections to the end
			CALL SortAerosolAtGridPointForCondensation ()
		
			!Calculate surrogate organic concentrations
			!IF(DoHydrophobicOrgDissolutionFlag .OR. DoHydrophilicOrgDissolutionFlag) &
				!CALL CalculateSurrogateGasConc(S1MassFracs, S2MassFracs, &
				!	S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
				!	S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)
	
			DO I = 1, NumTimeSteps
				!WRITE(*,*) "StepDissolutionAll Check 0"
				IF(DoDissolutionFlag .OR. DoCondensationFlag) &
					CALL StepCondensation (TimestepSize, 1)
		
				!WRITE(*,*) "StepDissolutionAll Check 1"
				IF(DoHydrophobicOrgDissolutionFlag) THEN
                                    CALL UpdateHydrophobicUNIFAC_all
				    CALL StepOrgCondensation (TimestepSize, 1)
                                ENDIF
		
				!WRITE(*,*) "StepDissolutionAll Check 2"
				IF(DoHydrophilicOrgDissolutionFlag) THEN
                                   CALL UpdateHydrophilicUNIFAC_all
				   CALL StepAqOrgCondensation (TimestepSize, 1)
                                ENDIF
				!WRITE(*,*) "StepDissolutionAll Check 3"
			END DO

			!Update individual gas-phase concentrations
			!IF(DoHydrophobicOrgDissolutionFlag .OR. &
                        !     DoHydrophilicOrgDissolutionFlag) &
			!	CALL UpdateGasPhaseOrganics(S1MassFracs, & 
                        !                S2MassFracs, S3MassFracs, &
                        !                 S4MassFracs, S5MassFracs, &
                        !                 S6MassFracs, S7MassFracs, &
                        !                 S8MassFracs, S9MassFracs, &
                        !                 S10MassFracs, InitialBurden)
		END IF

		!WRITE(*,*) "StepDissolutionAll Check 4"
		!Move particles across bin boundaries, if necessary
		CALL RegridAerosol ()		
	
	END IF !Empty Grid Cell
			
	!WRITE(*,*) "StepDissolutionAll Check 5"
END SUBROUTINE StepCondensationAll

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DissolutionFactors --				    !!
!!							    !!
!! Give the dissolution routines everything they need to    !!
!! calculate one iteration of one chemical, including:	    !!
!!							    !!
!!  1. S' / Effective Henry's Law Ratios		    !!
!!  2. Initial Appropriate Aerosol Concentrations	    !!
!!							    !!
!! ReturnType = 1 if fine, -1 if there is none of	    !!
!! the chemical in the gridpoint			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DissolutionFactors (WhichRxn, HowManyLinks, Temperature, & 
                               ReactionType, InitialConcs1, InitialConcs2, & 
                               InitialConcs3, SoverH, SoverH2, ReturnType)


	USE GridPointFields, ONLY :	GasPhaseChemicalRates,		&
					GetSatVapBurden,	        &
                                        GetGridCellChemBurden

	USE Aerosols,        ONLY :	Particles, Particle,		&
					ParticleArray, Molality

	USE ModelParameters, ONLY :	RstarMB, moles, grams, watermolecmass,&
					ProtonIndex, HydroxyIndex,	&
					Avogadro

	USE Chemistry,       ONLY :	HowManyAqChems,HowManyAqCations,&
					HowManyAqAnions,		&
					AqEquilibriaList,		&
					GasPhaseChemicalNames,	        &
					AqPhaseChemicalNames,		&
					AqCationNames,AqAnionNames,     & 
					FindChem, HowManyAqOrgChems

	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	!! External Variables, in
	INTEGER :: WhichRxn, HowManyLinks
	REAL*8  :: Temperature

	!! External Variables, out
	REAL*8  :: SoverH(HowManyLinks), SoverH2(HowManyLinks), &
                   InitialConcs1(HowManyLinks), &
	           InitialConcs2(HowManyLinks), InitialConcs3(HowManyLinks),&
                   RealFlag(HowManyLinks), GasSaturationConc(HowManyLinks)
	INTEGER :: ReturnType, ReactionType, Flag, IonA, IonB, ExpA, ExpB

	!! Internal Variables
	INTEGER :: InternalHowManyLinks, I, J, Allocation_Error, GasChemIndex
	REAL*8  :: II, H2SO4AqIndex, AnionSatConc,  NH4SaturationConc

	TYPE(Particle), POINTER :: CurrentParticle

	!! Identify the gas index of the condensing chemical
	GasChemIndex = DissolutionData(WhichRxn,1)

	!! An initial guess
	ReturnType = -1

	!! If Water Condensation, we need a couple of extra items...
	IF (ReactionType .EQ. 0) THEN
		II = GetSatVapBurden()
	END IF

	!! If there is a positive gas phase concentration, 
        !! then return positive return flag.
	IF (ReactionType .EQ. 0.) ReturnType = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Take the burdens and # of Particles from each of the aerosol particles !!
	CurrentParticle => Particles%First !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO J = 1, HowManyLinks

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Pre-Compute the Henry's Law Coefficient          !!
		!! Fill these placeholders (for non-water           !!
		!! condensation) with S' / H', which will           !!
		!! be held constant through the itnegration         !!
		!!					            !!
		!!					            !!
		!! Possible reaction types:                         !!
		!!                                                  !!
		!! 0. Water condensation (a special case for which  !!
		!!    we don't specify H's)			    !!
		!!						    !!
		!! 1. Dissolution of a non-dissociating species     !!
		!!						    !!
		!! 2. Dissolution of a dissociating species	    !!
		!!						    !!
		!! 3. Coupled reactions of types 1 & 2		    !!
		!!						    !!
		!! 4. Dissolution of a species that reacts with an  !!
		!!    ion upon entry.				    !!
		!!						    !!
		!! 5. Coupled reactions of types 1 & 4		    !!
		!! 7. Sulfate condensation                          !!
		!!NOTE: 6 (gas-solid reactions) is skipped          !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		SELECT CASE (ReactionType)

		!! Water condensing
		CASE (0)

			SoverH(J) = II * SatVapPressRatio (CurrentParticle, 1)

			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.				     !!
			InitialConcs1(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) * Avogadro

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

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.)	&
			ReturnType = 1

			!! NOTE: Factor of 1000 has units of g H2O/L soln (approximate!)
			SoverH(J) =																	 &
					SatVapPressRatio (CurrentParticle, INT(DissolutionData(WhichRxn,3))) &
					* 1000. / HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature	 &
					/ CurrentParticle%AqChems(1) / moles * grams / WaterMolecMass

			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) * Avogadro
			InitialConcs2(J) = 0.
			InitialConcs3(J) = 0.

		!! Dissociating dissolution
		CASE (2)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,4))) .GT. 0.)						&
			ReturnType = 1

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! The reaction doesn't use the Henry's Law formulation      !!
		!! because we assume it condenses directly into the dissociated
		!! form, which we must now formulate...	!!!!!!!!!!!!!!!!!!!!!!!
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
						*(CurrentParticle%GammaMixed(INT(DissolutionData(WhichRxn,12)))/Avogadro)**				 & ! 3 lines of
						(AqEquilibriaList(INT(DissolutionData(WhichRxn,12)),4)+									 & ! mean act
						 AqEquilibriaList(INT(DissolutionData(WhichRxn,12)),5))									 &
						 * Avogadro 
			
			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) * Avogadro
			InitialConcs2(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,4))) * Avogadro
			InitialConcs3(J) = 0.
			
		!! Dissociating and direct dissolution (NOT USED! Matt Alvarado 10/04/06)
		CASE (3)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,4))) .GT. 0.)						&
			ReturnType = 1

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! The reaction doesn't use  the Henry's Law formulation     !!
		!! because we assume it condenses directly into the dissociated
		!! form, which we must now formulate...	!!!!!!!!!!!!!!!!!!!!!!!
			IonA = INT(DissolutionData(WhichRxn,3))
			IonB = INT(DissolutionData(WhichRxn,4))
			ExpA = INT(DissolutionData(WhichRxn,13))
			ExpB = INT(DissolutionData(WhichRxn,14))

			!! This one is to be multiplied by Y(IonA)**ExpA * Y(IonB)**ExpB
			SoverH(J) =																						 &
						SatVapPressRatio (CurrentParticle , INT(DissolutionData(WhichRxn,3))) 				 &
						* (1000. * grams / moles / WaterMolecMass / CurrentParticle%AqChems(1))**(ExpA+ExpB) &
						/ HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature							 & 
						*(CurrentParticle%GammaMixed(INT(DissolutionData(WhichRxn,12)))/Avogadro)**				 & ! 3 lines of
						(AqEquilibriaList(INT(DissolutionData(WhichRxn,12)),4)+									 & ! mean act
						 AqEquilibriaList(INT(DissolutionData(WhichRxn,12)),5))									 &
						 * Avogadro 
	
			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) * Avogadro
			InitialConcs2(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,4))) * Avogadro


			SoverH2(J) = SatVapPressRatio (CurrentParticle,								 &
					INT(DissolutionData(INT(DissolutionData(WhichRxn,15)),3)))				 &
					* 1000. / HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature	 &
					/ CurrentParticle%AqChems(1) / moles * grams / WaterMolecMass		 

			!! Load the Y vector with the appropriate aerosol !!
			!! concentration.								  !!
			InitialConcs3(J) = CurrentParticle%AqChems(INT(DissolutionData(INT(DissolutionData(WhichRxn,15)),3))) * Avogadro


		!! Dissolution that incorporates an aqueous ion (NH3(g) & H+(aq) => NH4+(aq))
		CASE (4)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,2))) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.)						&
			ReturnType = 1

			IonA = FLOOR(DissolutionData(WhichRxn,8))
			IonB = FLOOR(DissolutionData(WhichRxn,9))
			ExpA = ANINT(10000.*(DissolutionData(WhichRxn,8)-FLOOR(DissolutionData(WhichRxn,8))))
			ExpB = ANINT(10000.*(DissolutionData(WhichRxn,9)-FLOOR(DissolutionData(WhichRxn,9))))

			IF (CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,2))) .LE. 0) THEN
				CALL ERROR ("A dissolution reaction (specified in Dissolution.in) specifies that the gas "//		&
							TRIM(GasPhaseChemicalNames(INT(DissolutionData(WhichRxn,1))))//								&
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
				InitialConcs1(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) * Avogadro
				InitialConcs2(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,2))) * Avogadro
				InitialConcs3(J) = 0.
			
			END IF

		!! Dissolution that incorporates an aqueous ion (NH3(g) & H+(aq) => NH4+(aq)) coupled with direct diss.
		CASE (5)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,2))) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.)						&
			ReturnType = 1

			IonA = FLOOR(DissolutionData(WhichRxn,8))
			IonB = FLOOR(DissolutionData(WhichRxn,9))
			ExpA = ANINT(10000.*(DissolutionData(WhichRxn,8)-FLOOR(DissolutionData(WhichRxn,8))))
			ExpB = ANINT(10000.*(DissolutionData(WhichRxn,9)-FLOOR(DissolutionData(WhichRxn,9))))

			IF (CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,2))) .LE. 0) THEN
				CALL ERROR ("A dissolution reaction (specified in Dissolution.in) specifies that the gas "//	&
							TRIM(GasPhaseChemicalNames(INT(DissolutionData(WhichRxn,1))))//							&
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
				InitialConcs1(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) * Avogadro
				InitialConcs2(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,2))) * Avogadro


				SoverH2(J) = SatVapPressRatio (CurrentParticle,								 &
						INT(DissolutionData(INT(DissolutionData(WhichRxn,15)),3)))				 &
						* 1000. / HenrysLaw (Temperature, WhichRxn) / RstarMB / Temperature	 &
						/ CurrentParticle%AqChems(1) / moles * grams / WaterMolecMass

				!! Load the Y vector with the appropriate aerosol !!
				!! concentration.								  !!
				InitialConcs3(J) = CurrentParticle%AqChems(INT(DissolutionData(INT(DissolutionData(WhichRxn,15)),3))) * Avogadro

			END IF
		
		!! Sulfate condensing
		CASE (7)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.	&
			    .AND. CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,4))) .GT. 0.)						&
			ReturnType = 1

			H2SO4AqIndex = FindChem("H2SO4", 1)
			SoverH(J) = 1.0

			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) * Avogadro
			InitialConcs2(J) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,4))) * Avogadro
			InitialConcs3(J) = 0.
		
		END SELECT

		CurrentParticle => CurrentParticle%Next
	END DO

	RETURN
END SUBROUTINE DissolutionFactors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Call LSODES or the equilibrium routine to integrate forward	 !!
!! the concentration of aerosol / gas-phase water.      	 !!
!! It must construct a joint array of all of the	         !!
!! aerosol plus each of the condensing species                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE StepCondensation (TimestepSize, NumTimeSteps)

	USE InfrastructuralCode, ONLY : INT2STR, ERROR,		&
					Warn, Transcript,	&
					REAL2STR

	USE GridPointFields,	 ONLY : getM,			        &
					GetGridCellChemBurden,		&
					ReplaceGridCellChemBurden,	&
					GetTemp,			&
					GetRelativeHumidity,	        &
					AddToGridCellChemBurden,        &
                                        GetSatVapBurden

	USE Aerosols,            ONLY : Particles, Particle,		&
					ParticleArray,			&
					RecalculateRadius,		&
					EqCoeff,			&
					Molality,			&
					FindElectrolyteEquilibrium,     &
					CalculateIonicStrength,	        &
				        KusikMeissner,			&
					SortAerosolAtGridPointForCondensation,&
				FindAqElectrolyteEquilibriumForGridPoint, &
				EquilibrateAqueousDissociationReaction

	USE Chemistry,		 ONLY : HowManyGasChems,		&
					HowManyEvolveGasChems,		&
					HowManyEvolveAqChems,		&
					GasPhaseChemicalNames,		&
					HowManyCondensingChems,		&
					AqEquilibriaList,		&
					HowManyAqChems,			&
					HowManyAqCations,		&
					FindChem,			&
					AqPhaseChemicalNames,           &
					HowManyAqEqReactions

	USE ModelParameters,     ONLY : Avogadro,			&
					AerosolWaterEquilibriumRHThreshold,&
					WaterEquilibriumIndex,		&
					ThermoBulkMode,			&
					DissolutionEquilibriumOrFlux,   &
					AqThermoEquilibriumError,	&
					WaterContentPrecision,		&
					micron,				&
					AmmoniaFlux,			&
					DoHydrophobicOrgDissolutionFlag,     &
					DoHydrophilicOrgDissolutionFlag,&
					ProtonIndex, &
					CombinedInorgDiss, &
					ZSROrgWater


	IMPLICIT NONE

	!! External Variables
	REAL*8  :: TimeStepSize		!! Should be in Seconds
	INTEGER :: NumTimeSteps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! LSODES Related Internal Variables !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8,  ALLOCATABLE	:: y(:), rwork(:)
	REAL*8			:: t,outtime,rtol,atol
	INTEGER, ALLOCATABLE	:: iwork(:)
	INTEGER			:: itol,itask,istate,iopt,mf,neq,LIW,LRW
		
	!! Non-LSODES Local Variables
	INTEGER	                :: I, J, K, L, LLLL, C, allocation_error, & 
                                   Rxn, HowManyLinks, YLen, GG,		  &
				   HowManyNonZeroJacobianTerms,           &
                                   GasChemIndex1, TimeStep,	          &
				   ReturnType, WhichDominantIon,          &
                                   WhichSecondIon , Flag, Dum, RxnType,   &
				   InternalReturnType, NH4Index, HIndex,  &
                                   OHIndex, GasSolidRxn
	REAL*8			:: II, JJ, KK, MM,                        &
                                   ChemTransferArray(HowManyGasChems),    &
                                   Temperature,	RH, NH3Index,             &
                                   OrgWaterContent, Store, Store3,Store4, &
                                   GasPhaseBurden, Electrolyte, WaterSatBurden
	REAL*8, ALLOCATABLE :: InitialConcsArray(:), ChemStorageArray(:,:)
	LOGICAL :: EmptyGridCell
	REAL*8  :: BeginTime,EndTime, ABCD
	REAL*4  :: etimearray(2)
	REAL*8  :: WaterTime, CondTime, AmmoniaTime

	LOGICAL, PARAMETER :: Scaffolding = .FALSE., OpenSystem = .FALSE., &
                               ForceEquilibrium = .TRUE.
	LOGICAL :: Equilibrated, UseEq

	TYPE(Particle), POINTER :: FirstParticle, CurrentParticle

	!! SCSCSCSC
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript(">>>Entering StepCondensationRxn<<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Find GRID POINT Values for NON-AEROSOL !!
	!! DEPENDENT Parameters                   !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Temperature = GetTemp()
	RH = GetRelativeHumidity()
		
	IF (Scaffolding) CALL TRANSCRIPT("Using LSODES to integrate dissolution.")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Everything should now be sorted           !!
	!! Back up to the head of the aerosol list   !!
	FirstParticle   => Particles%First
	CurrentParticle => FirstParticle

	!! COUNT the number of LINKS that we should condense upon
	HowManyLinks  =  1
15	IF(ASSOCIATED(CurrentParticle%Next)) THEN
		CurrentParticle => CurrentParticle%Next
		IF (CurrentParticle%NumberOfParticles .GT. 0 .AND. &
                   .NOT.CurrentParticle%Dry) THEN
			HowManyLinks  = HowManyLinks + 1
			GOTO 15
		END IF
	ENDIF
	CurrentParticle => FirstParticle


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set the lengths of the various !!
	!! arrays and allocate them !!!!!!!!
	YLen	= 12 + HowManyLinks*8	  !! As per the definition below
	NEQ		= 1  + 3*HowManyLinks	  
        !! Maximum value is for reaction type 5.  Allocate based on this.

	!! Temporarily set this to the maximum value 
        !! (which is also for reaction type 4)
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
!! This is indexed as follows:						    !!
!! ChemStorageArray(Rxn, 1)         = Change in Gas Phase Conc 1	    !!
!! ChemStorageArray(Rxn, 2)         = Change in Aerosol Conc 1 (1)	    !!
!! ...									    !!
!! ChemStorageArray(Rxn,HML+1)     = Change in Aerosol Conc 1 (HowManyLinks)!!
!! ChemStorageArray(Rxn,HML+2)     = Change in Aerosol Conc 2 (1)	    !!
!! ...									    !!
!! ChemStorageArray(Rxn,2*HML+1)   = Change in Aerosol Conc 2 (HowManyLinks)!!
!! ChemStorageArray(Rxn,2*HML+2)   = Change in Aerosol Conc 3 (1)	    !!
!! ...									    !!
!! ChemStorageArray(Rxn,3*HML+1)   = Change in Aerosol Conc 3 (HowManyLinks)!!
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
	CurrentParticle => Particles%First !!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set LSODES-Related Variable (Notes Excerpted !!
!! from LSODES file)			        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! itask  = an index specifying the task to be performed.		!!
!!          input only.  itask has the following values and meanings.	!!
!!          1  means normal computation of output values of y(t) at	!!
!!             t = tout (by overshooting and interpolating).		!!
!!          2  means take one step only and return.			!!
!!          3  means stop at the first internal mesh point at or	!!
!!             beyond t = tout and return.				!!
!!          4  means normal computation of output values of y(t) at	!!
!!             t = tout but without overshooting t = tcrit.		!!
!!             tcrit must be input as rwork(1).  tcrit may be equal to	!!
!!             or beyond tout, but not behind it in the direction of	!!
!!             integration.  this option is useful if the problem	!!
!!             has a singularity at or beyond t = tcrit.		!!
!!          5  means take one step, without passing tcrit, and return.	!!
!!             tcrit must be input as rwork(1).				!!
!!								     	!!
!!          note..  if itask = 4 or 5 and the solver reaches tcrit	!!
!!          (within roundoff), it will return t = tcrit (exactly) to	!!
!!          indicate this (unless itask = 4 and tout comes before tcrit,!!
!!          in which case answers at t = tout are returned first).	!!
		itask  = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
		istate = 1	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! iopt   = an integer flag to specify whether or not any optional      !!
!!          inputs are being used on this call.  input only.		!!
!!          the optional inputs are listed separately below.		!!
!!          iopt = 0 means no optional inputs are being used.		!!
!!                   default values will be used in all cases.		!!
!!          iopt = 1 means one or more optional inputs are being used.	!!
		iopt   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mf     = the method flag.  used only for input.		    !!
!!          the standard choices for mf are..			    !!
!!            mf = 10  for a nonstiff problem,			    !!
!!            mf = 21 or 22 for a stiff problem with ia/ja supplied !!
!!                     (21 if jac is supplied, 22 if not),	    !!
!!            mf = 121 for a stiff problem with jac supplied,	    !!
!!                     but not ia/ja,				    !!
!!            mf = 222 for a stiff problem with neither ia/ja nor   !!
!!                     jac supplied.				    !!
!!          the sparseness structure can be changed during the	    !!
!!          problem by making a call to lsodes with istate = 3.	    !!
		mf = 21  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
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
!!          an array of length neq.  input only.			!!
!!							       		!!
!!             the input parameters itol, rtol, and atol determine	!!
!!          the error control performed by the solver.  the solver will	!!
!!          control the vector e = (e(i)) of estimated local errors	!!
!!          in y, according to an inequality of the form		!!
!!                      rms-norm of ( e(i)/ewt(i) )   .le.   1,		!!
!!          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),		!!
!!          and the rms-norm (root-mean-square norm) here is		!!
!!          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))	!!
!!          is a vector of weights which must always be positive, and	!!
!!          the values of rtol and atol should all be non-negative.	!!
!!          the following table gives the types (scalar/array) of	!!
!!          rtol and atol, and the corresponding form of ewt(i).	!!
!!									!!
!!             itol    rtol       atol          ewt(i)			!!
!!              1     scalar     scalar     rtol*abs(y(i)) + atol	!!
!!              2     scalar     array      rtol*abs(y(i)) + atol(i)	!!
!!              3     array      scalar     rtol(i)*abs(y(i)) + atol	!!
!!              4     array      array      rtol(i)*abs(y(i)) + atol(i)	!!
!!									!!
!!          when either of these parameters is a scalar, it need not	!!
!!          be dimensioned in the user-s calling program.		!!
!!									!!
!!          if none of the above choices (with itol, rtol, and atol	!!
!!          fixed throughout the problem) is suitable, more general	!!
!!          error controls can be obtained by substituting		!!
!!          user-supplied routines for the setting of ewt and/or for	!!
!!          the norm calculation.  see part iv below.			!!
!!									
!!          if global errors are to be estimated by making a repeated	!!
!!         run on the same problem with smaller tolerances, then all	!!
!!          components of rtol and atol (i.e. of ewt) should be scaled	!!
!!          down uniformly.					        !!
		itol   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! This is a very stiff problem, as aerosol and the gas phase should 
	!! not be the same scale at all. 
	!!
	!! If the initial concentration is zero, some absolute tolerance
	!! is necessary (an error of 0 is unacceptable).
	!!
	!! The units of this are molecules
		atol   = 1.0

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Tell the program which gridpoint it's operating at: !!
                !! 02-17-2012 MJA Now useless, set to 1,1,1 to keep    !!
                !! array structure                                     !!
		Y(3*HowManyLinks+3) = 1. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Y(3*HowManyLinks+4) = 1.
		Y(3*HowManyLinks+5) = 1.


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! IF EQUILIBRIUM RATHER THAN FLUX, DO THAT HERE AND RETURN !!
		IF (.NOT. DissolutionEquilibriumOrFlux .AND. RH .LE. AerosolWaterEquilibriumRHThreshold) THEN
			WRITE(*,*) "Calling Equilibrium Routine"
			CALL EquilibrateGridPointAll(.TRUE., ForceWaterEquilibrium=.FALSE., WaterOpenSystem=.FALSE.)
			RETURN
		END IF
		
		
		IF (Scaffolding) WRITE(*,*) "Before Dissolution Loop"
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! LOOP OVER EACH DISSOLUTION REACTION !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! 02-17-2012 MJA 
                !! A. For H2O (Rxn = 1) and H2SO4 (Rxn = 2), the condensation 
                !!    is done here in StepCondensation (see below).
                !! B. For all other condensation/dissolution reactions
                !!       (i.e., Rxn = 3 or more) 
                !!        then these reactions are equilibrated using 
                !!        EquilibrateInorgGridPoint in 
                !!        CondensationRelatedFunctions.h. 
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO Rxn = 1, HowManyDissolutionReactions
   
		!If we are doing the combined inorganic dissolution
		IF(CombinedInorgDiss .AND. (Rxn .GE. 3)) THEN
			!Call the equilibrator only once
			IF(Rxn .EQ. 3) THEN
				CALL EquilibrateInorgGridPoint()
				CYCLE
			ELSE	
				!Skip all the other rxns
				!since they are equilibrated above
				CYCLE
			END IF
		END IF
		
		!Send an error if this tries to do ammonia in flux form
	     IF(AmmoniaFlux .EQ. .TRUE. .AND. CombinedInorgDiss .EQ. .FALSE.) CALL ERROR("Operator split method only works if Ammonia is done as eq.")


		!! The allowed error for water content may be different than
		!! for H2SO4
		IF (Rxn .EQ. 1) THEN
			rtol   = WaterContentPrecision **(1./ NumTimeSteps)
		ELSE
			rtol   = (AqThermoEquilibriumError) **(1./NumTimeSteps)
		END IF

	!! If the reaction is coupled to another one, the direct dissolution
	!! is done at the same time as the other one, meaning we skip this
	!! direct one, UNLESS the reaction is NH3 & H+ => NH4+ when NH4+ is low
	!! (that's the Dum flag). Cycle...

		!! Reset timestep to the correct position
		T = (TimeStep - 1) * TimeStepSize
		
	!! If Water is condensing and the RH Flag is below the critical value,
	!! then simply equilibrate it. (Uses total water content: org + inorg)
		IF (Rxn .EQ. 1 .AND.  &
                    RH .LE. AerosolWaterEquilibriumRHThreshold) THEN
			Equilibrated = .FALSE.
			
			DO WHILE (NOT(Equilibrated))
				CurrentParticle => Particles%First
				Equilibrated = .TRUE.
			
				DO J = 1, HowManyLinks
				
				        !Skip empty particles
				        IF(CurrentParticle%NumberofParticles .LE. 0.) CYCLE

					LLLL = 1
					ReturnType = 1
					DO WHILE (ReturnType .NE. 2)
					
						ABCD = CurrentParticle%AqChems(1)

						CALL EqAllWater(CurrentParticle, ReturnType, OpenSystem)

						IF(ReturnType .NE. 2) Equilibrated = .FALSE.
						IF(LLLL .GT. 50 .AND. &
						  ABS((ABCD-CurrentParticle%AqChems(1))/CurrentParticle%AqChems(1)) .LE. 0.01) THEN
						  ReturnType  = 2
						  Equilibrated = .TRUE.
						END IF

						LLLL = LLLL+1
						IF(LLLL .GT. 1000) THEN
							CALL WARN("Over 1000 water iterations in StepCondensation().")
							CALL WARN("RH equals :"//trim(real2str(RH)))
						END IF
					END DO
						
					CurrentParticle => CurrentParticle%Next
				END DO
			END DO
			CYCLE
		END IF
				
		!! Identify the gas index of the condensing chemical
		GasChemIndex1 = DissolutionData(Rxn,1)
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! This is the maximum number of internal steps allowed. !!
		!! The default is 500 and that's not enough.		 !!
		iwork(6) = 300000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICAL				     !!
!! The Y-Vector we will integrate in LSODES looks like this:		     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    (1): Gas Phase Concenatration of Chemical DissolutionData(Rxn,1)	     !!
!!    (2 - HowManyLinks+1)                                                   !!
!!       : Aerosol Concentrations of Chemical DissolutionData(Rxn,3)	     !!
!!    (HowManyLinks+2   - 2*HowManyLinks+1)	                             !!
!!       : Aerosol Concentrations of Chemical DissolutionData(Rxn,2) or	     !!
!!	    DissolutionData(Rxn,4), as appropriate			     !!
!!    (2*HowManyLinks+2 - 3*HowManyLinks+1)                                  !!
!!       : Concentration of undissociated species for coupled rxn	     !!
!!    (3*HowManyLinks+2) : Intentionally Empty				     !!
!!    (3*HowManyLinks+3 - 3*HowManyLinks+5) : Grid Point Specification       !!
!!        02-17-2012 MJA Now useless, set to 1,1,1			     !!
!!    (3*HowManyLinks+6 - 3*HowManyLinks+7)	: GasChemIndex 1 & 2	     !!
!!    (3*HowManyLinks+8 - 4*HowManyLinks+7)                                  !!
!!       : Condensation Rates for Chemical One	   			     !!
!!    (4*HowManyLinks+8 - 5*HowManyLinks+7)                                  !!
!!       : How Many Particles in Each					     !!
!!    (5*HowManyLinks+8 - 6*HowManyLinks+9)                                  !!
!!       : Pre-computed effective Henry's Law Coeffs		             !!
!!    (6*HowManyLinks+10) : Flag for what type of reaction		     !!
!!    (6*HowManyLinks+11) : "-1" (as a reference for finding HowManyLinks)   !!
!!    (6*HowManyLinks+12) : Dissolution Reaction Number			     !!
!!    (6*HowManyLinks+13 - 7*HowManyLinks+12)                                !!
!!       : Pre-computed effective Henry's Law Coeffs for coupled rxn         !!
!!    (7*HowManyLinks+13 - 8*HowManyLink+12)                                 !!
!!       : Flags for HNO3 and HCl dissolution                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! HowManyLinks which the number of aerosol / sections			     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialize the Concentrations from the local chemical concentrations: !! 
			Y(1) = GetGridCellChemBurden (GasChemIndex1)

			IF (DissolutionData(Rxn,11) .EQ. 5) &
			Y(3*HowManyLinks+2) = &
                            GetGridCellChemBurden (GasChemIndex1)

			!! Tell LSODES which chemical we are to consider
			Y(3*HowManyLinks+6) = GasChemIndex1

			!! The condensation rates
			Y(3*HowManyLinks+8:4*HowManyLinks+7) = &
                            CondensationRateCoefficients (Rxn, HowManyLinks)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Reaction Type !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! There are several reaction types in this context, !!
			!! which are different than elsewhere:		     !!
			!!  0. Water Condensation			     !!
			!!  1. Direct Dissolution of a single chemical	     !!
			!!  2. Dissolution of a single species that dissociates
			!!  3. Coupled reactions of types 1 & 2		     !!
			!!  4. Dissolution of a species that binds with      !!
                        !!      an ion                                       !!
			!!  5. Coupled reactions of types 1 & 4	             !!
			!!  7. Sulfate condensation (eq. gas conc. is 0)     !!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
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
			CALL DissolutionFactors (Rxn, HowManyLinks, Temperature, INT(Y(6*HowManyLinks+10)), 	 &
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
30			CALL LSODES (CondensationODEEvaluator,neq,y,t, &
                                     outtime,itol,rtol,atol,itask, &
		                     istate,iopt,rwork,lrw,iwork,liw, &
                                     CondensationJacobianEvaluator,mf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF (SCAFFOLDING) WRITE(*,*) "After Condensation LSODES"
!! If the maximum number of steps occurred, then warn the user 
!! and push the integration back into LSODES until the appropriate end time 
!! is reached.
			IF (ISTATE .EQ. -1 ) THEN 
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
!! have considered all of the chemicals     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This is indexed as follows:						     !!
!! ChemStorageArray(Rxn, 1)         = Change in Gas Phase Conc 1	     !!
!! ChemStorageArray(Rxn, 2)         = Change in Aerosol Conc 1 (1)           !!
!!									     !!
!! ChemStorageArray(Rxn, HML+1)     = Change in Aerosol Conc 1 (HowManyLinks)!!
!! ChemStorageArray(Rxn, HML+2)     = Change in Aerosol Conc 2 (1)	     !!
!!									     !!
!! ChemStorageArray(Rxn, 2*HML+1)   = Change in Aerosol Conc 2 (HowManyLinks)!!
!! ChemStorageArray(Rxn, 2*HML+2)   = Change in Gas Phase Conc 2	     !!
			DO K = 1, NEQ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		          ChemStorageArray(Rxn,K) = Y(K) - InitialConcsArray(K)
			END DO

		END DO  !! Ends the loop over all the dissolution reactions
		IF (Scaffolding) WRITE(*,*) "After Dissolution Loop"

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Replace the background chemical concentrations !!
		DO I = 1, HowManyDissolutionReactions
			IF (SCAFFOLDING) WRITE(*,*) I, ChemStorageArray(I,1)
			!Skip replace if gas is to be held constant
			IF(INT(DissolutionData(I,1)) .LE. HowManyEvolveGasChems) THEN
				CALL AddToGridCellChemBurden (INT(DissolutionData(I,1)), ChemStorageArray(I,1))
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

			IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(INT(DissolutionData(J,3))) + &
															ChemStorageArray(J,I+1) / Avogadro

			CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) + &
															ChemStorageArray(J,I+1) / Avogadro

			IF(CurrentParticle%AqChems(INT(DissolutionData(J,3))) .LT. 0.0) THEN
				!WRITE(*,*) "Negative Cation error, reaction #", J
				!WRITE(*,*) AqPhaseChemicalNames(DissolutionData(J,3)), CurrentParticle%AqChems(DissolutionData(J,3))
			END IF

			RxnType = INT(DissolutionData(J,11)) 
				
			SELECT CASE (RxnType)
				
				CASE (2)
				
				IF(CurrentParticle%AqChems(INT(DissolutionData(J,4))) .LT. 0.0) THEN
					!WRITE(*,*) "Negative Anion error, reaction #", J
					!WRITE(*,*) AqPhaseChemicalNames(DissolutionData(J,4)), CurrentParticle%AqChems(DissolutionData(J,4))
				END IF

				IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(INT(DissolutionData(J,4))) + &
															ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
				CurrentParticle%AqChems(INT(DissolutionData(J,4))) = CurrentParticle%AqChems(INT(DissolutionData(J,4))) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro

			20	IF(CurrentParticle%AqChems(INT(DissolutionData(J,3))) .LT. 0.0 .OR. &
				   CurrentParticle%AqChems(INT(DissolutionData(J,4))) .LT. 0.0) THEN
					!CALL Transcript("In CondensationIntegrator, an acid reaction went below zero.")
					
					!Reset concs to previous value
					!Cation
					CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) - &
															ChemStorageArray(J,I+1) / Avogadro
					!Anion
					CurrentParticle%AqChems(INT(DissolutionData(J,4))) = CurrentParticle%AqChems(INT(DissolutionData(J,4))) - &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
					!Gas
					CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)
					
					Electrolyte = 0.0
					DO C = 1, HowManyDissolutionReactions
						!Find the corresponding gas-solid reaction, if it exists
						IF(DissolutionData(J,1) .EQ. DissolutionData(C,2)) THEN
							Electrolyte = DissolutionData(C,3)
							GasSolidRxn = C
						END IF
					END DO
					
					IF(INT(Electrolyte) .NE. 0. .AND. & !Gas-solid Reaction exists
					   CurrentParticle%AqChems(INT(Electrolyte)) .GT. ABS(ChemStorageArray(J,I+1) / Avogadro)) THEN !There is enough salt
						!WRITE(*,*) "Salt Evaporation!"
		
						!Evaporate Electrolyte (remember chemstoragearray is negative here!)
						CurrentParticle%AqChems(INT(Electrolyte)) = CurrentParticle%AqChems(INT(Electrolyte)) + ChemStorageArray(J,I+1) / Avogadro
						
						!Add acid gas to gas phase
						CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), -1.0*ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)
						
						!Add ammonia to gas phase
						CALL AddToGridCellChemBurden (INT(DissolutionData(GasSolidRxn,1)), -1.0*ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)
						
						!Equilibrate gas-phase ammonia (to prevent oscillations)
						DO C = 1, HowManyDissolutionReactions
							IF (GasPhaseChemicalNames(INT(DissolutionData(C,1))) .EQ. "NH3" &
								.AND. INT(DissolutionData(C,11)) .EQ.  4 &
								.AND. AmmoniaFlux .EQ. .FALSE. ) THEN
								IF (SCAFFOLDING) WRITE(*,*) "Before Ammonia Equilibrate, particle #", J
				
								IF (SCAFFOLDING) WRITE(*,*) "Equilibrating Ammonia!"
								CALL KusikMeissner (CurrentParticle)
								Temperature = GetTemp()
								GasPhaseBurden = GetGridCellChemBurden (INT(DissolutionData(C,1)))
                                                                WaterSatBurden = GetSatVapBurden ()
								CALL EquilibrateDissolutionReaction (GasPhaseBurden, C, CurrentParticle, Temperature, ReturnType, WaterSatBurden = WaterSatBurden, GasPhaseBurden2=0.0, OptErrorTolerance=AqThermoEquilibriumError/2.)
								CALL ReplaceGridCellChemBurden (INT(DissolutionData(C,1)), GasPhaseBurden)		
							END IF
						END DO

						!Cycle (should just skip rest of check)
						GOTO 20
					END IF
						
					!If ammonia salt doesn't fix it,
					!decide if anion or cation is limiting
					Store3 = CurrentParticle%AqChems(INT(DissolutionData(J,3)))
					Store4 = CurrentParticle%AqChems(INT(DissolutionData(J,4)))

					IF(Store3 .LE. Store4) THEN !Cation limiting
										
						!WRITE(*,*) "Cation limiting! This is always H+"
						IF (Store3 .NE. 0.0) THEN
							!WRITE(*,*) "Store3 not zero!"
							!Only allow ions to reach 0, not go under
							!Cation
							!WRITE(*,*) "Cation: ", AqPhaseChemicalNames(DissolutionData(J,3))
							CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) - &
															0.9*Store3
							!WRITE(*,*) "Corrected Cation conc: ", CurrentParticle%AqChems(DissolutionData(J,3))
							!This keeps H+ above zero to prevent errors with ammonia
							!Anion
							!WRITE(*,*) "Anion: ", AqPhaseChemicalNames(DissolutionData(J,4))
							CurrentParticle%AqChems(INT(DissolutionData(J,4))) = CurrentParticle%AqChems(INT(DissolutionData(J,4))) - &
																0.9*Store3
							!Acid(Gas)
							!WRITE(*,*) "Gas: ", GasPhaseChemicalNames(INT(DissolutionData(J,1)))
							CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), 0.9*Store3*Avogadro*CurrentParticle%NumberofParticles)
						
						ELSE IF(Store3 .EQ. 0.0 .AND. Store4 .NE. 0.0) THEN
							!WRITE(*,*) "Store3 zero!"
							!Push cation (H+) slightly above 0
							!Cation
							!WRITE(*,*) "Skippy Cation: ", AqPhaseChemicalNames(DissolutionData(J,3))
							CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) + &
															0.01*Store4
							!WRITE(*,*) "Corrected Cation conc: ", CurrentParticle%AqChems(DissolutionData(J,3))
							!This keeps H+ above zero to prevent errors with ammonia
							!Anion
							!WRITE(*,*) "Anion: ", AqPhaseChemicalNames(DissolutionData(J,4))
							CurrentParticle%AqChems(INT(DissolutionData(J,4))) = CurrentParticle%AqChems(INT(DissolutionData(J,4))) + &
																0.01*Store4
							!Acid(Gas)
							!WRITE(*,*) "Gas: ", GasPhaseChemicalNames(INT(DissolutionData(J,1)))
							CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), -0.01*Store4*Avogadro*CurrentParticle%NumberofParticles)
							
						END IF

					ELSE !Anion Limiting

							WRITE(*,*) "Anion Limiting!"
							WRITE(*,*) "Crude Correction!"
							!Cation
							CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) - &
															Store4
							!Anion
							CurrentParticle%AqChems(INT(DissolutionData(J,4))) = CurrentParticle%AqChems(INT(DissolutionData(J,4))) - &
																Store4
							!Acid(Gas)
							CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), Store4*Avogadro*CurrentParticle%NumberofParticles)
						END IF
					END IF
					
				
				CASE (3) !No zero check!
				CurrentParticle%AqChems(INT(DissolutionData(J,4))) = CurrentParticle%AqChems(INT(DissolutionData(J,4))) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro


				CurrentParticle%AqChems(INT(DissolutionData(INT(DissolutionData(J,15)),3))) =								&
															CurrentParticle%AqChems(INT(DissolutionData(INT(DissolutionData(J,15)),3))) + &
															ChemStorageArray(J,I+2*HowManyLinks+1) / Avogadro

				CASE (4)
				IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(INT(DissolutionData(J,2))) + ChemStorageArray(J,I+HowManyLinks+1) / Avogadro 
				
				CurrentParticle%AqChems(INT(DissolutionData(J,2))) = CurrentParticle%AqChems(INT(DissolutionData(J,2))) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro

				IF (CurrentParticle%AqChems(INT(DissolutionData(J,2))) .LT. 0) THEN
					
					!CALL Transcript("In CondensationIntegrator, there was a cation concentration that went below zero.")
					!Reset concs to previous value
					!NH4+
					CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) - &
															ChemStorageArray(J,I+1) / Avogadro
					!H+
					CurrentParticle%AqChems(INT(DissolutionData(J,2))) = CurrentParticle%AqChems(INT(DissolutionData(J,2))) - &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
					!NH3(Gas)
					CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)

					!Only allow H+ to reach 0, not go under
					!Store previous H+
					Store = CurrentParticle%AqChems(INT(DissolutionData(J,2)))
					!NH4+
					CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) + &
															Store
					!H+
					CurrentParticle%AqChems(INT(DissolutionData(J,2))) = CurrentParticle%AqChems(INT(DissolutionData(J,2))) - &
																Store
					!NH3(Gas)
					CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), -1.*Store*Avogadro*CurrentParticle%NumberofParticles)

				END IF
				
				CASE (5)
				IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(INT(DissolutionData(J,2))) +  ChemStorageArray(J,I+HowManyLinks+1) / Avogadro 
				CurrentParticle%AqChems(INT(DissolutionData(J,2))) = CurrentParticle%AqChems(INT(DissolutionData(J,2))) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
		
				IF (CurrentParticle%AqChems(INT(DissolutionData(J,2))) .LT. 0) THEN
					!Reset concs to previous value
					!NH4+
					CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) - &
															ChemStorageArray(J,I+1) / Avogadro
					!H+
					CurrentParticle%AqChems(INT(DissolutionData(J,2))) = CurrentParticle%AqChems(INT(DissolutionData(J,2))) - &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
					!NH3(Gas)
					CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), ChemStorageArray(J,I+1)*CurrentParticle%NumberofParticles)

					!Only allow H+ to reach 0, not go under
					!NH4+
					CurrentParticle%AqChems(INT(DissolutionData(J,3))) = CurrentParticle%AqChems(INT(DissolutionData(J,3))) + &
															CurrentParticle%AqChems(INT(DissolutionData(J,2)))
					!H+
					CurrentParticle%AqChems(INT(DissolutionData(J,2))) = CurrentParticle%AqChems(INT(DissolutionData(J,2))) - &
																CurrentParticle%AqChems(INT(DissolutionData(J,2)))
					!NH3(Gas)
					CALL AddToGridCellChemBurden (INT(DissolutionData(J,1)), -1.*CurrentParticle%AqChems(INT(DissolutionData(J,2)))*Avogadro*CurrentParticle%NumberofParticles)

				END IF
				
				IF (Scaffolding) WRITE(*,*) J, CurrentParticle%AqChems(INT(DissolutionData(INT(DissolutionData(J,15)),3))) + ChemStorageArray(J,I+2*HowManyLinks+1) / Avogadro 
				CurrentParticle%AqChems(INT(DissolutionData(INT(DissolutionData(J,15)),3))) =								&
															CurrentParticle%AqChems(INT(DissolutionData(INT(DissolutionData(J,15)),3))) + &
															ChemStorageArray(J,I+2*HowManyLinks+1) / Avogadro
				
				CASE (7)
				CurrentParticle%AqChems(INT(DissolutionData(J,4))) = CurrentParticle%AqChems(INT(DissolutionData(J,4))) + &
																ChemStorageArray(J,I+HowManyLinks+1) / Avogadro
			
			END SELECT
		END DO 
		
		!! Re-Equilibrate
		CALL FindElectrolyteEquilibrium(CurrentParticle, UpdateThermo = .TRUE., FirstEquilibration=.FALSE., ReturnType=GG)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate a JACOBIAN for Condensation !!
!!                                       !!
!! As per LSODES's specification, only   !!
!! calculate at the stated row I         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 02-17-2012 MJA Since this is only used for water and sulfate, all other!!
!! reaction types besides 0, 1, and 7 will now give errors.               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
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
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Reaction Type !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!  0. Water Condensation			     !!
	!!  1. Direct Dissolution of a single chemical	     !!
	!!  2. Dissolution of a single species that dissociates
	!!  3. Coupled reactions of types 1 & 2		     !!
	!!  4. Dissolution of a species that binds with      !!
        !!      an ion                                       !!
	!!  5. Coupled reactions of types 1 & 4	             !!
	!!  7. Sulfate condensation (eq. gas conc. is 0)     !!
        !!  Only 0, 1, and 7 types should be done here.      !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
		
		CASE (2:6) !! These types not included
                   CALL ERROR("Cannot deal with reactions of type 2-6 in CondensationJacobianEvaluator")

        
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 02-17-2012 MJA Since this is only used for water and sulfate, all other!!
!! reaction types besides 0, 1, and 7 will now give errors.               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE CondensationODEEvaluator (neq, t, y, ydot)

		USE Aerosols,            ONLY : Particles, Particle, &
                                                 ParticleArray
		USE Chemistry,		 ONLY : HowManyGasChems
		USE InfrastructuralCode, ONLY : Transcript, ERROR, & 
                                                 REAL2STR, INT2STR
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICAL				     !!
!! The Y-Vector we will integrate in LSODES looks like this:		     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    (1): Gas Phase Concenatration of Chemical DissolutionData(Rxn,1)	     !!
!!    (2 - HowManyLinks+1)                                                   !!
!!       : Aerosol Concentrations of Chemical DissolutionData(Rxn,3)	     !!
!!    (HowManyLinks+2   - 2*HowManyLinks+1)	                             !!
!!       : Aerosol Concentrations of Chemical DissolutionData(Rxn,2) or	     !!
!!	    DissolutionData(Rxn,4), as appropriate			     !!
!!    (2*HowManyLinks+2 - 3*HowManyLinks+1)                                  !!
!!       : Concentration of undissociated species for coupled rxn	     !!
!!    (3*HowManyLinks+2) : Intentionally Empty				     !!
!!    (3*HowManyLinks+3 - 3*HowManyLinks+5) : Grid Point Specification       !!
!!        02-17-2012 MJA Now useless, set to 1,1,1			     !!
!!    (3*HowManyLinks+6 - 3*HowManyLinks+7)	: GasChemIndex 1 & 2	     !!
!!    (3*HowManyLinks+8 - 4*HowManyLinks+7)                                  !!
!!       : Condensation Rates for Chemical One	   			     !!
!!    (4*HowManyLinks+8 - 5*HowManyLinks+7)                                  !!
!!       : How Many Particles in Each					     !!
!!    (5*HowManyLinks+8 - 6*HowManyLinks+9)                                  !!
!!       : Pre-computed effective Henry's Law Coeffs		             !!
!!    (6*HowManyLinks+10) : Flag for what type of reaction		     !!
!!    (6*HowManyLinks+11) : "-1" (as a reference for finding HowManyLinks)   !!
!!    (6*HowManyLinks+12) : Dissolution Reaction Number			     !!
!!    (6*HowManyLinks+13 - 7*HowManyLinks+12)                                !!
!!       : Pre-computed effective Henry's Law Coeffs for coupled rxn         !!
!!    (7*HowManyLinks+13 - 8*HowManyLink+12)                                 !!
!!       : Flags for HNO3 and HCl dissolution                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! HowManyLinks which the number of aerosol / sections			     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        
        !! There are several structures for the columns beyond the first
	!! that depend on reaction type
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Reaction Type !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!  0. Water Condensation			     !!
	!!  1. Direct Dissolution of a single chemical	     !!
	!!  2. Dissolution of a single species that dissociates
	!!  3. Coupled reactions of types 1 & 2		     !!
	!!  4. Dissolution of a species that binds with      !!
        !!      an ion                                       !!
	!!  5. Coupled reactions of types 1 & 4	             !!
	!!  7. Sulfate condensation (eq. gas conc. is 0)     !!
        !!  Only 0, 1, and 7 types should be done here.      !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! WATER CONDENSATION
		CASE (0)

			DO I = 2, neq
				YDOT(I) = Y(3*HowManyLinks+6+I) *			& ! The Condensation Rate
						  (Y(1) -				& ! The Gas Phase Conc
						   Y(I)/(Y(I)+Y(I+HowManyLinks))*	& ! Raout's Law
						   Y(5*HowManyLinks+6+I))		! S_w' * C_s,w

				YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I)

			END DO

		!! DISSOLUTION of a single, non-dissociating species
		CASE (1)

			DO I = 2, neq
				
				YDOT(I) = Y(3*HowManyLinks+6+I) *	&		! The Condensation Rate
						  (Y(1) -		&		! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) * Y(I))	! S_i' * c_s,i / H_i

				YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I)
				
			END DO


		!! DISSOLUTION of a single, dissociating species
		CASE (2:6)
                     CALL ERROR("Cannot deal with reactions of type 2-6 in CondensationODEEvaluator")

		!! SULFATE CONDENSATION
		!! Assumes that H2SO4 vapor pressure = 0 
		CASE (7)

			DO I = 2, neq-HowManyLinks
				!WRITE(*,*) Y(3*HowManyLinks+6+I)
				YDOT(I) = Y(3*HowManyLinks+6+I) *	& ! The Condensation Rate
						  (Y(1)) * &		 ! The Gas Phase Conc
					   DissolutionData(INT(Y(6*HowManyLinks+12)),13) !Cation Stoichiometry
				YDOT(I+HowManyLinks) = YDOT(I)*(DissolutionData(INT(Y(6*HowManyLinks+12)),14)/DissolutionData(INT(Y(6*HowManyLinks+12)),13))

				YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I) / DissolutionData(INT(Y(6*HowManyLinks+12)),13)
			END DO

		END SELECT

	RETURN
END SUBROUTINE CondensationODEEvaluator
