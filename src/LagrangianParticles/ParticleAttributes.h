!2345678901234567890123456789012345678901234567890123456789012345678901234567890
!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ParticleAttributes.h
!! Contains functions for calculating particle properties (density, radius, mass
!! curvature correction, Knudsen number, Reynolds number, etc.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							      	!!
!!									       	!!
!! Month  Year   Name              Description					!!
!! 07     2006   Matt Alvarado     Began Update History				!!
!! 07/24  2006   Matt Alvarado     Removed "ISNAN" to fit pgf90			!!
!! 09/07  2006   Matt Alvarado     Changed RecalculateRadius to give lower edge	!!
!!				      when number conc. = 0.			!!
!! 11/01  2006   Matt Alvarado     Created ShellRefIndAndRad
!! 05/18  2007   Matt Alvarado     Revised ShellRefIndAndRad
!!				   Created AerosolOptProp
!! 07/13  2007   Matt Alvarado     Added line to  ShellRefIndAndRad
!!	 			      to keep core rad below shell rad
!! 07/17  2007   Matt Alvarado     Added line to SurfaceTension
!!				      to keep value above that of pure Ethanol
!! 07/18  2007   Matt Alvarado     Removed calls to ParticleDensity from
!!				      CurvatureCorrection and OrgCurvatureCorrection
!!				      Just use values in particle data structure instead
!!				   In ParticleDensity, if ion loop fails, force
!!				      density of ionic solution to 1.0 g/cm3
!! 07/26  2007   Matt Alvarado     Fixed AbsCoreRad calculation in ShellRefIndAndRad
!! 09/07  2007   Matt Alvarado     Fixed shell refractive index calcultion
!! 09/17  2007   Matt Alvarado     Removed write statement from ParticleDensity
!! 10/01  2007   Matt Alvarado     Set ParticleDensity to skip bins with 0 particles
!! 10/03  2007   Matt Alvarado     Set Denomsum limit in SurfaceTension to 1.0e40
!! 10/04  2007   Matt Alvarado     Protected ShellRefIndAndRad from near zero number concentrations
!! 10/11  2007   Matt Alvarado     Set HypotheticalElectrolyteConcentrations to warn of charge
!!                                 imbalance, but not stop.
!! 10/15  2007   Matt Alvarado     Fixed Assymetry parameter in ShellRefIndAndRad
!! 08/30  2010   Matt Alvarado     Changed ShellRefIndandRad, AerosolOpticalDepth,
!!                                   and AerosolOptProp to allow 62 photolysis wavelengths
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE ResetAllOriginPressureLevels (OriginPressureLevel, X, Y, Z)
!! 2. FUNCTION ParticleKnudsen (InParticle)
!! 3. FUNCTION ParticleKnudsenForEnergy (InParticle)
!! 4. FUNCTION CunninghamSlipCorrectionFactor (InParticle)
!! 5. FUNCTION AerosolPH (InParticle)
!! 6. FUNCTION ParticleDensity (InParticle) !Sets radius and effective radius
!! 7. FUNCTION ParticleMass (InParticle)
!! 8. FUNCTION GetPatwardhandAndKumarWeightings (InParticle)
!! 9. FUNCTION HypotheticalElectrolyteConcentrations (InParticle)
!!10. FUNCTION TerminalVelocity (InParticle)
!!11. SUBROUTINE RecalculateRadius (InParticle) !Calls particle density
!!12. FUNCTION CurvatureCorrection (InParticle)
!!13. FUNCTION SurfaceTension (InParticle)
!!14. FUNCTION ReynoldsNumber (InParticle)
!!15. REAL*8 FUNCTION OrgCurvatureCorrection (InParticle)
!!16. REAL*8 FUNCTION MolalityAqCarbon(InParticle)
!!17. REAL*8 FUNCTION AqueousSolutionMass
!!18. SUBROUTINE ShellRefIndAndRad(InParticle)
!!19. REAL*8 FUNCTION AerosolOpticalDepth(XP, YP, ZP)
!!20. SUBROUTINE AerosolOptProp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
SUBROUTINE ResetAllOriginPressureLevels (OriginPressureLevel, X, Y, Z)

	IMPLICIT NONE

	Type(Particle),POINTER :: Particle
	INTEGER :: X, Y, Z
	REAL    :: OriginPressureLevel

	Particle => Particles(X,Y,Z)%First

10	IF (ASSOCIATED(Particle)) Particle%OriginPressureLevel = OriginPressureLevel
		
	IF (ASSOCIATED(Particle%Next)) THEN
		Particle => Particle%Next
		GOTO 10
	END IF

	RETURN
END SUBROUTINE ResetAllOriginPressureLevels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Knudsen Number of the Particle	!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleKnudsen (InParticle)

	USE GridPointFields, ONLY : MeanFreePathOfAir

	IMPLICIT NONE

	Type(Particle),POINTER :: InParticle

	ParticleKnudsen = MeanFreePathOfAir (InParticle%Location(1),InParticle%Location(2),InParticle%Location(3)) &
					  / InParticle%EffectiveRadius

END FUNCTION ParticleKnudsen


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Knudsen Number Of Energy of the Particle !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleKnudsenForEnergy (InParticle)

	USE GridPointFields, ONLY : MolecularThermalDiffusivity, GetThermalVelocityOfAir

    IMPLICIT NONE

	Type(Particle),POINTER :: InParticle

	ParticleKnudsenForEnergy = 3 *																			  &
			MolecularThermalDiffusivity(InParticle%Location(1),InParticle%Location(2),InParticle%Location(3)) &
			/ GetThermalVelocityOfAir(InParticle%Location(1),InParticle%Location(2),InParticle%Location(3))   &
			/ InParticle%EffectiveRadius

	RETURN
END FUNCTION ParticleKnudsenForEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the CUNNINGHAM SLIP CORRECTION FACTOR, which !!
!! has a number of possible forms.  This function selects !!
!! the appropriate one and returns the value.			  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CunninghamSlipCorrectionFactor (InParticle)


    IMPLICIT NONE

	Type(Particle),POINTER :: InParticle
	REAL*8		   :: Knudsen

	Knudsen = ParticleKnudsen (InParticle)

	!! From Kasten 1968 (via Jacobson, 1997)
	CunninghamSlipCorrectionFactor = 1 + Knudsen * (1.249 + 0.42 * dexp(-0.87 / Knudsen))
	
	RETURN
END FUNCTION CunninghamSlipCorrectionFactor



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns the pH of an aerosol. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION AerosolPH (InParticle)

	USE ModelParameters, ONLY : ProtonIndex

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle

	!!Internal Variables
	REAL*8 :: ProtonMolarity

	IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) THEN
		AerosolPH = 14.
		RETURN
	END IF

	IF (InParticle%AqChems(ProtonIndex) .GT. 0) THEN
		!Molarity in mol/L solution
		ProtonMolarity = InParticle%AqChems(ProtonIndex)*InParticle%SolutionDensity*1000/AqueousSolutionMass(InParticle)
		
		AerosolPH = -1.*LOG10(ProtonMolarity)
	ELSE
		AerosolPH = 14.
	END IF

END FUNCTION AerosolPH



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the PARTICLE DENSITY for a lagrangian particle !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleDensity (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro,SmallestAerosolPossible
	USE Chemistry,			 ONLY : IFSOLID, AqMolecularMass,HowManyAqChems,HowManyAqCations,		&
									HowManyAqAnions, HowManyAqEqReactions, AqEquilibriaList,		&
									AqAnionCharge, AqCationCharge , &
									HowManyOrgChems, OrgMolecularMass, OrgDensity, SolidSaltDensity, &
									HowManyAqOrgChems, AqOrgMolecularMass, AqOrgDensity
	USE InfrastructuralCode, ONLY : ERROR,REAL2STR, WARN

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I,K, SkipFlag, C
	REAL*8  :: Volume, Mass, ionMass, PSIij, TempElectrolytes(HowManyAqChems),	&
			   BinaryDensity, NumerSum, DenomSum, TempMolality,					&
  			   PatKumWeights(HowManyAqEqReactions,3), XX, YY, Residual, dR, &
			    SolutionDensity, InsolMass, InsolVolume, InsolDensity, InorgSolnDensity
 
	REAL*8, PARAMETER :: ResidualThreshold = 0.0025 !**3.

	IF (InParticle%NumberOfParticles .LE. 0. ) THEN
		ParticleDensity = 1.0
		InParticle%EmbryoRadius    = 0.
		InParticle%EffectiveRadius = 0.
		InParticle%InsolubleRadius = 0.
		InParticle%ParticleDensity = 1.
		InParticle%SolutionDensity = 1.
		InParticle%InsolubleDensity = 1.
		InParticle%InorgSolnDensity = 1.
		RETURN
	END IF


	Mass = 0 ; Volume = 0; ionMass = 0; NumerSum = 0; DenomSum = 0

	!! Get electrolyte equivalents of solution containing dissociated ions,
	!! which assumes a two-parameter fit proscribed in Resch, 1995 
	!! (assuming that the density of water does not change with Temperature).
	!! This is a SOLUTION density, which includes the incorporated water.
	!! So we scale up by the mass of the water as well!
	TempElectrolytes = HypotheticalElectrolyteConcentrations (InParticle)
	!For some reason, This routine returns NaN if there isn't a write
	!statement here. MJA, 070507
	WRITE(1234,*) TempElectrolytes(1)

	!! Get the weightings for the Patwardhan and Kumar mixing rule
	PatKumWeights = GetPatwardhandAndKumarWeightings (InParticle)
	
	!! Loop over each reaction and consider the appropriate equation
	DO I = 1, HowManyAqEqReactions
		
		!! Some quasi-electrolytes (such as HSO4-) are not quantifyable in this
		!! density parameterization and should have been pushed to upper level
		!! electrolytes anyways (such as H2SO4).  So we ignore those electrolytes!
		IF (AqEquilibriaList(I,1) .LE. HowManyAqChems .AND. AqEquilibriaList(I,1) .GT. 0.) THEN ! GT 0 disallows water

			!! This is the equivalent molality of the associated electrolyte
			!Write(*,*) "Temp Electrolytes: ", TempElectrolytes(AqEquilibriaList(I,1))
                        TempMolality = TempElectrolytes(AqEquilibriaList(I,1)) / InParticle%AqChems(1) / &
						   (moles / grams * AqMolecularMass(1)) * 1000.

			Mass    = Mass    + AqMolecularMass(AqEquilibriaList(I,1))*TempElectrolytes(AqEquilibriaList(I,1))
			ionMass = ionMass + AqMolecularMass(AqEquilibriaList(I,1))*TempElectrolytes(AqEquilibriaList(I,1))
	
			IF(IsNaN (TempMolality)) THEN
				WRITE(*,*) "Problem in Temp Molality"
				WRITE(*,*) InParticle%NumberOfParticles, I, AqEquilibriaList(I,13)
				WRITE(*,*) TempElectrolytes(AqEquilibriaList(I,1)), TempMolality, InParticle%AqChems(1)
			END IF

			!WRITE(*,*) TempMolality, I
			!! Binary density is an exponential function of the two parameters input into the system
			IF(TempMolality .LT. 100.0) THEN
				
				!WRITE(*,*) "Before Exp, Rxn: ", I
				!WRITE(*,*) "TempMolality: ", TempMolality
				BinaryDensity = AqEquilibriaList(I,12) + EXP(-1.*AqEquilibriaList(I,13)*TempMolality)*		&
						   (1.*grams/cm/cm/cm-AqEquilibriaList(I,12))
				
				!This prevents the binary density of NH4OH from going below 0.9 g/cm3
				IF(AqEquilibriaList(I,13) .LT. 0. .AND. BinaryDensity .LT. 0.9) BinaryDensity = 0.9
						
			ELSE IF(TempMolality .GE. 100.0 .AND. AqEquilibriaList(I,13) .GT. 0.) THEN
				!If the molality gets very large, and exponent goes to zero
				!(This keeps EXP from giving an overflow error)
				BinaryDensity = AqEquilibriaList(I,12)
			ELSE !For NH4OH at high molality, keep it above 0.9 g/cm3
				BinaryDensity = 0.9	
			END IF

			!! Employ Mixing Rule of Patwardhan and Kumar (AIChE J. 39, 711-714) to Mix these ions appropriately
			PSIij =  1000. * PatKumWeights(I,1) + AqMolecularMass(AqEquilibriaList(I,1)) * TempMolality
			IF(IsNaN(PatKumWeights(I,1))) THEN
				CALL WARN("Problem in PatKumWeights! Get NaN!")
			END IF
			!! Composite Density is the Numerator Sum divided by the Denomenator Sum, although we will track it as Volume
			NumerSum = NumerSum + PSIij
			DenomSum = DenomSum + PSIij / BinaryDensity
		END IF
	END DO

	!WRITE(*,*) "After Ion loop"
	!PAUSE
	!! Now use the composite density to add volume to the system
	Volume = Volume + (ionMass + AqMolecularMass(1)*InParticle%AqChems(1)) * DenomSum / NumerSum
	Mass   = Mass   + AqMolecularMass(1)*InParticle%AqChems(1)

	!WRITE(*,*) "Ion solution Volume: ", Volume, DenomSum, NumerSum
	
	!If for some reason the above doesn't work, set inorganic solution density to 1.0
	IF(ISNaN(Volume)) THEN
		CALL WARN("Setting ion solution density to 1.0 g/cm3")
		Volume = Mass 
	END IF

	!This only includes water and inorganic ions (needed for optical property calculation)
	InorgSolnDensity = Mass/Volume
	!WRITE(*,*) "Inorganic Solution Density: ", InorgSolnDensity

	!Now add in contributions from Hydrophilic Organics (assumed constant density)
	!! This assumes a fixed density for these aqueous phase chemicals
	DO I = 1, HowManyAqOrgChems
		Mass   = Mass   + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
		Volume = Volume + InParticle%AqOrgChems(I)*AqOrgMolecularMass(I)/AqOrgDensity(I) 
	END DO
	
	!WRITE(*,*) "After Aq Org loop", Volume
	IF(ISNaN(Volume)) STOP

	!! This is purely the solution density (water, ions, non-solid electrolytes (e.g. NH3), aq organics)
	SolutionDensity =  Mass/Volume

	!!Calculate radius of the solution embryo
	InParticle%EmbryoRadius = (0.75 * Volume / pi)**0.33333333333
	!WRITE(*,*) "Embryo Radius: ", InParticle%EmbryoRadius
	
	!! Consider the contributions of Solid salts (except ions and water)
	InsolVolume = 0.
	InsolMass = 0.
	DO I = 2, HowManyAqChems
		SkipFlag = 0
		!Skip non-solid electrolytes (NH3, CO2, etc).
		DO C = 1, HowManyAqEqReactions
			IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				SkipFlag = 1
			END IF
		END DO
			
		IF(SkipFlag .EQ. 0) THEN
			InsolMass   = InsolMass   + AqMolecularMass(I)*InParticle%AqChems(I)
			InsolVolume = InsolVolume + InParticle%AqChems(I)*AqMolecularMass(I)/SolidSaltDensity(I)

			Mass   = Mass   + AqMolecularMass(I)*InParticle%AqChems(I)
			Volume = Volume + InParticle%AqChems(I)*AqMolecularMass(I)/SolidSaltDensity(I)
		END IF
	END DO
	
	!WRITE(*,*) "After Solid Salts loop", Volume
	IF(ISNaN(Volume)) STOP
	
	!Now add in controbutions from hydrophobic organics	
	DO I = 1, HowManyOrgChems		
		InsolMass   = InsolMass   + OrgMolecularMass(I)*InParticle%OrgChems(I)
		InsolVolume = InsolVolume + InParticle%OrgChems(I)*OrgMolecularMass(I)/OrgDensity(I) 
	
		!WRITE(*,*) I, InsolVolume
		Mass   = Mass   + OrgMolecularMass(I)*InParticle%OrgChems(I)
		Volume = Volume + InParticle%OrgChems(I)*OrgMolecularMass(I)/OrgDensity(I) 
	END DO
		
	!WRITE(*,*) "After Hydrophobic Org loop", Volume
	IF(ISNaN(Volume)) STOP
	
	!!This is the average density of the non-aqeous part (solid salts, BC, and hydrophobic organics)
	InsolDensity = InsolMass/InsolVolume
	InParticle%InsolubleRadius = (0.75 * InsolVolume / pi)**0.33333333333

	!! This is the average density of the particle
	ParticleDensity = Mass / Volume
	IF(Mass .LT. 0.0 ) WRITE(*,*) Mass, Volume

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! RESET EFFECTIVE RADIUS, too. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	InParticle%EffectiveRadius = (0.75 * Volume / pi)**0.33333333333
	!WRITE(*,*) "Volumes: ", Volume, InsolVolume
	!! If this returns nonsense, then probably it wasn't ever pushed into thermodynamic equilibrium
	!! and there are no ions, which will lead to infinite numbers!
	!IF (IsNaN(InParticle%EffectiveRadius)) THEN
	!	!CALL DumpParticleContentsAtError(InParticle, "NoIonsInParticleDensityParticle.txt",InDensity=.FALSE.,		&
	!	!								 InDoWaterActivity=.FALSE.,InDoSurfaceTension=.FALSE.)
	!	IF (InParticle%Sectional .AND. InParticle%NumberOfParticles .LE. 0) &
	!	CALL ERROR ("Called ParticleDensity() for a sectional distribution with negative particle concentration.")
	!	IF (InParticle%AqChems(1) .EQ. 0) THEN
	!		CALL ERROR ("Called ParticleDensity() for an aerosol with no water which is really not ok.  Please look into it!")
	!	ELSE
	!		WRITE(*,*) InParticle%AqChems(1), InParticle%NumberOfParticles, Mass
	!		CALL ERROR ("Called ParticleDensity() for an aerosol with no ions, leading me to believe that it is not in "// &
	!					"thermodynamic equilibrium (there should either be ions or no water), or that no equilibrium "//   &
	!					"reactions were defined for any of the constituent species.  Please look into it!")
	!	END IF
	!END IF

	InParticle%InorgSolnDensity = InorgSolnDensity 
	InParticle%SolutionDensity = SolutionDensity
	InParticle%InsolubleDensity = InsolDensity
	InParticle%ParticleDensity = ParticleDensity

	!WRITE(*,*) "End Particle Density", Volume, Mass
	
	
	RETURN
END FUNCTION ParticleDensity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Particle Mass !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleMass (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro
	USE Chemistry,			 ONLY : AqMolecularMass, HowManyAqChems,	&
									HowManyAqAnions, HowManyAqCations,	&
									AqCationMass, AqAnionMass, &
									HowManyOrgChems, OrgMolecularMass, &
									HowManyAqOrgChems, AqOrgMolecularMass 

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I

	ParticleMass = 0. 

	!! Consider the contributions of Aqueous Chemicals (except ions)
	!! This assumes a fixed density for these aqueous phase chemicals
	DO I = 1, HowManyAqChems
		ParticleMass   = ParticleMass + AqMolecularMass(I)*InParticle%AqChems(I)
	END DO

	DO I = 1, HowManyAqCations
		ParticleMass   = ParticleMass + AqCationMass(I)*InParticle%AqChems(HowManyAqChems+I)
	END DO

	DO I = 1, HowManyAqAnions
		ParticleMass   = ParticleMass + AqAnionMass(I)*InParticle%AqChems(HowManyAqChems+HowManyAqCations+I)
	END DO

	DO I = 1, HowManyOrgChems
		ParticleMass   = ParticleMass + OrgMolecularMass(I)*InParticle%OrgChems(I)
	END DO

	DO I = 1, HowManyAqOrgChems
		ParticleMass   = ParticleMass + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
	END DO

	RETURN
END FUNCTION ParticleMass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get a vector of weightings for the Patwardhan and Kumar	   !!
!! y_ij factors corrected for partial dissociation, incomplete !!
!! reaction sets, etc.										   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns 1). y_ij  2). Cation Index 3). Anion Index		   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Patwardhan and Kumar (AIChE J. 39, 711-714)				   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GetPatwardhandAndKumarWeightings (InParticle)

	USE Chemistry, ONLY : HowManyAqEqReactions,			&
						  HowManyAqChems,				&
						  HowManyAqCations,				&
						  HowManyAqAnions,				&
						  AqEquilibriaList,				&
						  AqCationCharge,				&
						  AqAnionCharge, &
							AqPhaseChemicalNames
	USE InfrastructuralCode, ONLY : WARN

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: GetPatwardhandAndKumarWeightings(HowManyAqEqReactions,3)
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: K, CatIndex, AnIndex
	REAL*8  :: II, JJ,														  &
			   CationicStrFr(HowManyAqCations),AnionicStrFr(HowManyAqAnions), &
			   CationChgFrac(HowManyAqCations),AnionChgFrac(HowManyAqAnions)

	!! Update the ionic strength
	CALL CalculateIonicStrength(InParticle)
	!WRITE(*,*) "Ionic Strength: ", InParticle%IonicStr
	IF(IsNAN(InParticle%IonicStr)) CALL WARN("Ionic Strength Error in PKweightings!")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make Charge Fraction Vectors !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Calculate Cationic First:
	II = 0.
	DO K = 1, HowManyAqCations
		CationChgFrac(K) = CationMolality(K, InParticle) * AqCationCharge(K)
		II               = II + CationChgFrac(K)
	END DO
	DO K = 1, HowManyAqCations
		CationChgFrac(K) = CationChgFrac(K) /II
		IF(IsNAN(CationChgFrac(K))) CALL WARN("CationChgFrac error in PKweightings!")
	END DO

	!! And Anionic Next:
	II = 0.
	DO K = 1, HowManyAqAnions
		AnionChgFrac(K) = AnionMolality(K, InParticle) * AqAnionCharge(K)
		II              = II + AnionChgFrac(K)
	END DO
	DO K = 1, HowManyAqAnions
		AnionChgFrac(K) = AnionChgFrac(K) / II
		IF(IsNAN(AnionChgFrac(K))) CALL WARN("AnionChgFrac error in PKweightings!")
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make Ionic Strength Fraction Vectors !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (Inparticle%IonicStr .EQ. 0.) THEN
		DO K = 1, HowManyAqCations
			CationicStrFr(K) = 0.
		END DO
		DO K = 1, HowManyAqAnions
			AnionicStrFr(K)  = 0.
		END DO
	ELSE
		DO K = 1, HowManyAqCations
			CationicStrFr(K) = CationMolality(K, InParticle) * AqCationCharge(K) * AqCationCharge(K) / 2. / InParticle%IonicStr
			IF(IsNAN(CationicStrFr(K))) CALL WARN("CationicStrFr error in PKweightings!")

		END DO
		DO K = 1, HowManyAqAnions
			AnionicStrFr(K)  = AnionMolality(K, InParticle)  * AqAnionCharge(K) * AqAnionCharge(K) /2. / InParticle%IonicStr
			IF(IsNAN(AnionicStrFr(K))) CALL WARN("AnionicStrFr error in PKweightings!")
		END DO
	END IF

	
	!! This will be the sum of all the y_ij Patwardhan and Kumar mixing rule terms
	JJ = 0.

	!! Loop over all reactions
	DO K = 1, HowManyAqEqReactions
		

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! For partial dissociation, consider the reactions as if the partially !!
		!! dissociated form exchanges directly with the fully dissociated form, !!
		!! which is not how it would have been written.							!!
		IF(AqEquilibriaList(K,1) .LE. HowManyAqChems) THEN ! electrolyte uncharged
			CatIndex = AqEquilibriaList(K,2)
			AnIndex  = AqEquilibriaList(K,3)
		ELSE											   ! electrolyte charged

			!! Is the partially dissociated "electrolyte" a cation or an anion?
			IF (AqEquilibriaList(K,1) .LE. HowManyAqChems+HowManyAqCations) THEN  ! It's a Cation
				CatIndex = AqEquilibriaList(K,1)-HowManyAqChems
				AnIndex  = AqEquilibriaList(K,3)
			ELSE																  ! It's an Anion
				CatIndex = AqEquilibriaList(K,2)
				AnIndex  = AqEquilibriaList(K,1)-HowManyAqChems-HowManyAqCations
			END IF
		END IF
			
		!! Store the indices used to calculate this	
		GetPatwardhandAndKumarWeightings(K,2) = REAL(CatIndex)
		GetPatwardhandAndKumarWeightings(K,3) = REAL(AnIndex)

		
		!Exclude Levitocite reaction from this calculation
		IF(AqEquilibriaList(K,22) .NE. 0) THEN
			IF (AqPhaseChemicalNames(AqEquilibriaList(K,1)) .EQ. "(NH4)3H(SO4)2") THEN
				GetPatwardhandAndKumarWeightings(K,1) = 0.
			END IF
		
		!Normal Calculation
		ELSE
			!! Calculate the Appropriate Weightings for the Patwardhan and Kumar (AIChE J. 39, 711-714)
			!! Mixed Activity Coefficient:
			!!
			!!   ln(Aw_mix) = Sum[(i,j), (Isf[i]*Cf[j]+Isf[j]*Cf[i])*ln(Aw[ij])]
			GetPatwardhandAndKumarWeightings(K,1) = AnionicStrFr(AnIndex)*CationChgFrac(CatIndex) + &
												CationicStrFr(CatIndex)*AnionChgFrac(AnIndex)
		
		END IF
		
		JJ = JJ + GetPatwardhandAndKumarWeightings(K,1)
	END DO

	!! Normalize the y_ij estimates to equal one
	DO K = 1, HowManyAqEqReactions
		GetPatwardhandAndKumarWeightings(K,1) = GetPatwardhandAndKumarWeightings(K,1) / JJ 
	END DO

	RETURN

END FUNCTION GetPatwardhandAndKumarWeightings


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Estimate Dissociated Electrolyte Concentrations !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For some routines, it is necessary to find the hypothetical
!! electrolyte concentrations that a set of dissociated ions 
!! would make up.  This function provides a vector of such electrolyte
!! abundances, in the form of a vector matching the AqChems vector
!! in the particle structure.
!!
!! The algorithm is taken from Jacobson's Fundamentals of Atmospheric
!! Modeling, 1999 page 494
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION HypotheticalElectrolyteConcentrations (InParticle)

	USE Chemistry, ONLY : HowManyAqChems, HowManyAqCations, HowManyAqAnions, HowManyAqEqReactions, AqEquilibriaList, & 
AqPhaseChemicalNames, AqAnionNames, AqCationNames
	USE ModelParameters, ONLY : AqThermoNumbAllowedIterationsForHypLoop,moles,WaterEquilibriumIndex
	USE InfrastructuralCode, ONLY : ERROR, REAL2STR, TRANSCRIPT

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: HypotheticalElectrolyteConcentrations(HowManyAqChems)
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	REAL*8  :: TempChem(HowManyAqChems+HowManyAqCations+HowManyAqAnions), DeltaMoles, Sum, ContentSum, ExcludeRxnInit
	INTEGER :: I, J, K, C, ExcludeRxn 

	!! counters and settings
	K = 0
	ExcludeRxn = 0

	!! First set all electrolytes to zero, unless they are non-solid
10	ContentSum=0
	DO I = 1, HowManyAqChems
		TempChem(I) = 0.
		
		!Check for non-solid electrolyte (NH3, CO2, etc), and add in.
		DO C = 1, HowManyAqEqReactions
			IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				TempChem(I) = InParticle%AqChems(I)
			END IF
		END DO			
		 
		ContentSum = ContentSum + InParticle%AqChems(I)
	END DO

	!WRITE(*,*) "Check 1 ContentSum: ", ContentSum

	!! Then take the Ion Concentrations as is.
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		TempChem(I) = InParticle%AqChems(I)
		ContentSum  = ContentSum + InParticle%AqChems(I)
	END DO

	!WRITE(*,*) "Check 2 ContentSum: ", ContentSum

	!WRITE(*,*) "Check 3"

	!! Loop over the reactions and adjust the equations
	DO J = 1, 2
	DO I = 1, HowManyAqEqReactions
		  !WRITE(*,*) I, ExcludeRxn, AqEquilibriaList(I,1), HowManyAqChems
		IF (I .EQ. ExcludeRxn) THEN
			CYCLE

		!! Some reactions may have ions serving as electrolytes
		!! We must identify these, and then push them the other direction
		ELSE IF (AqEquilibriaList(I,1) .GT. HowManyAqChems) THEN

			!WRITE(*,*) TempChem (AqEquilibriaList(I,1))
			!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), TempChem (HowManyAqChems+AqEquilibriaList(I,2)), TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,4)
			!WRITE(*,*) AqAnionNames(AqEquilibriaList(I,3)), TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)), TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,5)
		
			TempChem (HowManyAqChems+AqEquilibriaList(I,2))					   =		&
					TempChem (HowManyAqChems+AqEquilibriaList(I,2)) + TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,4)
			TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)) =		&
					TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)) + &
					TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,5)
			TempChem (AqEquilibriaList(I,1))                                   = TempChem(AqEquilibriaList(I,1)) - TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,5)
			!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), TempChem (HowManyAqChems+AqEquilibriaList(I,2))
			!WRITE(*,*) AqAnionNames(AqEquilibriaList(I,3)), TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3))
			!WRITE(*,*) TempChem (AqEquilibriaList(I,1))


		!! Other Reactions are only water dissociations
		ELSE IF (AqEquilibriaList(I,1) .LE. 0) THEN

			!! The adjustment parameter is limited by the concentration of each ion divided
			!! by the appropriate stoicheometric coefficient.
			DeltaMoles  = MIN(TempChem (HowManyAqChems+AqEquilibriaList(I,2)) / AqEquilibriaList(I,4),					&
						      TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)) / AqEquilibriaList(I,5))

			!! Now adjust all of the concentrations appropriately
			IF (DeltaMoles .GT. 0) THEN
				!! Don't adjust back towards water at all, just get rid of it
				TempChem (HowManyAqChems+AqEquilibriaList(I,2))					   = &
					TempChem (HowManyAqChems+AqEquilibriaList(I,2)) - DeltaMoles * AqEquilibriaList(I,4)
				TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)) = &
					TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)) - DeltaMoles * AqEquilibriaList(I,5)
			END IF

		!! And still others are canonical
		ELSE

			IF (AqEquilibriaList(I,22) .EQ. 0) THEN
			
				!! The adjustment parameter is limited by the concentration of each ion divided
				!! by the appropriate stoicheometric coefficient.
				DeltaMoles  = MIN(TempChem (HowManyAqChems+AqEquilibriaList(I,2)) / AqEquilibriaList(I,4),					&
						      TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)) / AqEquilibriaList(I,5))
                                !WRITE(*,*) "Delta Moles: ", DeltaMoles
				!! Now adjust all of the concentrations appropriately
				IF (DeltaMoles .GT. 0) THEN
					TempChem (AqEquilibriaList(I,1))                           = TempChem(AqEquilibriaList(I,1)) + DeltaMoles
					TempChem (HowManyAqChems+AqEquilibriaList(I,2))			   = TempChem (HowManyAqChems+AqEquilibriaList(I,2)) - &
																			 DeltaMoles * AqEquilibriaList(I,4)
					TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)) = TempChem (HowManyAqChems+HowManyAqCations + &
																	AqEquilibriaList(I,3)) - DeltaMoles * AqEquilibriaList(I,5)
				END IF
		

			ELSE
				!Levitocite case: Skip this
				CYCLE						
			END IF
		
		END IF
	END DO
	END DO

	!WRITE(*,*) "Check 4"

	!! Check to make sure have accounted for all of the ions
	SUM = 0.
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		!IF(IsNan(TempChem(I))) THEN
		!	WRITE(*,*) "Error in HypElec!"
		!	STOP
		!END IF
		SUM = SUM + TempChem(I)
	END DO

	!! This is an error
	!! I changed zero limit to 1.0e-3 from 1.0e-12 to eliminate some problems (Matt Alvarado, 5/10/2006)
	!! Note that this does not error so long as the amount of ions left is less than 1/1000 of the total aqueous moles!
	!! For dry particles, just runs!
	IF (SUM .GT. ContentSum*1.e-3) THEN
	!! It should not be caught in an infinite loop!

		WRITE(*,*) SUM, ContentSum*1.0e-3

		!! Show us what was wrong if the calculation doesn't work right.  (Probably the ContentSum*1.e-12 is not right below.)
		CALL TRANSCRIPT("")
		CALL TRANSCRIPT("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		CALL TRANSCRIPT("+++ There are problems in Hypothetical Electrolyte Concentrations +++")
		CALL TRANSCRIPT("")

		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions

		IF (TempChem(I) .gt. 0) then
		IF (I .LE. HowManyAqChems) CALL TRANSCRIPT(TRIM(AqPhaseChemicalNames(I))//TRIM(REAL2STR(TempChem(I)/moles)))
		IF (I .GT. HowManyAqChems .AND. I .LE. HowManyAqChems+HowManyAqCations) &
				CALL TRANSCRIPT(TRIM(AqCationNames(I-HowManyAqChems))//TRIM(REAL2STR(TempChem(I)/moles)))
		IF (I .GT. HowManyAqChems+HowManyAqCations .AND. I .LE. HowManyAqChems+HowManyAqCations+HowManyAqAnions) &
				CALL TRANSCRIPT(TRIM(AqAnionNames(I-HowManyAqChems-HowManyAqCations))//TRIM(REAL2STR(TempChem(I)/moles)))
		END IF

		END DO

		CALL DumpParticleContentsAtError(InParticle, "FailedHypotheticalElectrolyteConcentrationsParticle.txt",	&
						 InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., &
                                                 InDoWaterActivity=.FALSE.,InDensity=.FALSE.,InDoSurfaceTension=.FALSE., &
                                                 InIonicStrength=.TRUE.)
		!Removed this so that routines warns of charge imbalance 
		!in transcript, but does not stop
		!CALL ERROR ("HypotheticalElectrolyteConcentrations() is stuck in what appears to be an infinite loop.")

	END IF

	!! Transfer the good stuff into a shorter vector for output
	DO I = 1, HowManyAqChems
		!WRITE(*,*) I, TempChem(I)
		HypotheticalElectrolyteConcentrations(I) = TempChem(I)
	END DO

	RETURN
END FUNCTION HypotheticalElectrolyteConcentrations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the TERMINAL FALL VELOCITY for particles !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION TerminalVelocity (InParticle)

	USE GridPointFields, ONLY : GetAirDensity,GetDynamicViscosityOfAir

	IMPLICIT NONE

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Internal Vairables
	REAL*8 :: Dum

	TerminalVelocity = ReynoldsNumber (InParticle) * 0.5  &
					   * GetDynamicViscosityOfAir(InParticle%Location(1),InParticle%Location(2),InParticle%Location(3))  &
					   / GetAirDensity(InParticle%Location(1),InParticle%Location(2),InParticle%Location(3))			 &
					   / InParticle%EffectiveRadius

	RETURN
END FUNCTION TerminalVelocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Recalculate the Radius of the Particle !!
SUBROUTINE RecalculateRadius (InParticle)

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, IFAQ

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Type local variables
	REAL*8  :: II

	IF (InParticle%NumberOfParticles .EQ. 0. ) THEN
		IF (InParticle%Edges(1) .GT. 0.) THEN
			InParticle%EmbryoRadius    = InParticle%Edges(1)
			InParticle%InsolubleRadius = 0.
			InParticle%EffectiveRadius = InParticle%Edges(1)
		ELSE
			InParticle%EmbryoRadius    = 0.01*InParticle%Edges(2)
			InParticle%InsolubleRadius = 0.
			InParticle%EffectiveRadius = 0.01*InParticle%Edges(2)
		END IF

		InParticle%ParticleDensity = 1.
		InParticle%SolutionDensity = 1.
		InParticle%InsolubleDensity = 1.
		InParticle%InorgSolnDensity = 1.
		RETURN
	END IF

	!! Calling the Density Routine Resets the Radius, so just do that
	II = ParticleDensity (InParticle)	

	RETURN
END SUBROUTINE RecalculateRadius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The correction to pressure from curvature, used for water activity and in condensation. !!
!! It is also known as the "Kelvin Effect"												   !!
!!																						   !!
!! The optional argument PerturbWater is a multiple by which to change the water content   !!
!! so we can get a sense of the derivative of the water residual.						   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CurvatureCorrection (InParticle) 

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : Rstar, micron, ThermoBulkMode, IgnoreCorrectionstoDiffusivity
	USE GridPointFields, ONLY : SurfaceTensionofWater

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER  :: InParticle
	REAL*8 :: TotalMoles
	INTEGER :: I

	!!Type Internal variables
	REAL*8 :: TotalDens
	!WRITE(*,*) "Curvature Started", SurfaceTension(InParticle), ParticleDensity(InParticle)
	
	!Ignore this for bulk mode, no particles, or if ordered to.
	IF (ThermoBulkMode .OR. InParticle%NumberofParticles .LE. 0.0 .OR. IgnoreCorrectionstoDiffusivity) THEN
		CurvatureCorrection = 1.
		
	ELSE
		!! Do this in stages because ParticleDensity forces a radius + solution density update..."
		!WRITE(*,*) "Before Dens"
		TotalDens  = InParticle%SolutionDensity
		!WRITE(*,*) "After Dens" 
		
		!WRITE(*,*) "Curvature Correction: ", InParticle%NumberofParticles

		CurvatureCorrection = 2.*AqMolecularMass(1)/(InParticle%Temperature*Rstar*InParticle%SolutionDensity)
		!WRITE(*,*) CurvatureCorrection, InParticle%EmbryoRadius, SurfaceTension(InParticle)
		
		
		CurvatureCorrection = EXP(CurvatureCorrection*SurfaceTension(InParticle)/InParticle%EmbryoRadius)

	END IF

	!WRITE(*,*) "Curvature Okay", CurvatureCorrection 
	RETURN
	
END FUNCTION CurvatureCorrection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Surface Tension depends on Temperature and Electrolytic Properties  !!
!! of the Aerosol's constituent species.							   !!
!!																	   !!
!! Follow Li and Lu (2001): Chem. Eng. Sci., 56: 2879-2888, who use	   !!
!! a Gibbs dividing surface, define surface excess, etc.  It is a 2-   !!
!! parameter-per-electrolyte system and the data is taken in when the  !!
!! electrolytes are defined.										   !!
!!																	   !!
!! Many other systems seem to use the sig = sig_w + Beta I formulation !!
!! from Pruppacher and Klett (1997), but that is significantly less    !!
!! effective than this technique.									   !!
!!
!! Units are dyn / cm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION SurfaceTension (InParticle)

	USE GridPointFields, ONLY : SurfaceTensionOfWater, GetRelativeHumidity
	USE Chemistry, ONLY : AqEquilibriaList, HowManyAqEqReactions
	USE ModelParameters 
	USE InfrastructuralCode, ONLY : Error, REAL2STR
	IMPLICIT NONE

	!! Input Variables
	TYPE (Particle),POINTER :: InParticle

	INTEGER :: I
	REAL*8  :: II, Sum1, DenomSum, TotalIonMolality, Test

	IF (InParticle%NumberOfParticles .EQ. 0.0) RETURN

	!! CALL FindElectrolyteEquilibrium (InParticle)
	CALL CalculateIonicStrength(InParticle)
    CALL KusikMeissner(InParticle)

	IF(UseDonnanSurfaceTension) THEN !If true, Donnan; If False, Jacobson formulas
		!!Donnan's Parameterization
		!! 1. + SUM (K_i * Activity_i)
		DenomSum = 1.
		DO I = 1, HowManyAqEqReactions
			DenomSum = DenomSum + MeanActivity(I, InParticle) * AqEquilibriaList(I,15)
		END DO

		!! Probably this only happens if the RH is actually 0.
		IF (Denomsum .GT. 1.0e40) THEN
			WRITE(*,*) DenomSum, InParticle%numberofParticles
			CALL DumpParticleContentsAtError(InParticle, "BadActivitiesInSurfaceTensionParticle.txt", &
                                                 InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., &
                                                 InDoWaterActivity=.TRUE.,InDensity=.TRUE.,InDoSurfaceTension=.FALSE., &
                                                 InIonicStrength=.TRUE.)
			CALL ERROR ("In SurfaceTension(), the mean activity coefficients appear to have grown without bounds.  "// &
					"This is likely because the water is leaching out of the system towards zero.  Perhaps the "// &
					"relative humidity for the particle is too small?  I think it is "//trim(real2str(100*         &
					GetRelativeHumidity(InParticle%location(1),InParticle%location(2),InParticle%location(3))))//  &
					"%.  Another possibility is that you are in the water equilibration routine and it is dumping "// &
					"water from the particle to try to reduce the water activity, and it got too dry.")
		END IF

		!! Loop over each reaction and consider the appropriate equation
		Sum1 = 0.
		DO I = 1, HowManyAqEqReactions

			!! Ignore dissociation of 
			!IF ((AqEquilibriaList(I,1) .GT. HowManyAqChems) .OR. (AqEquilibriaList(I,1) .LT. 0)) CYCLE

			!! AqIonGrid(I,14) is Surface Excess and AqIonGrid(I,15) is Equilibrium Adsorption Coefficient
			!WRITE(*,*) "Mean Act.: ", MeanActivity(I, InParticle)
			!WRITE(*,*) "Ka: ", AqEquilibriaList(I,15)
			!WRITE(*,*) "DenomSum: ", DenomSum
			Sum1 = Sum1 + AqEquilibriaList(I,14) * LOG (1. - AqEquilibriaList(I,15)	&  
				* MeanActivity(I, InParticle) / DenomSum) 

		END DO

		!! Surface Tesion of the Solution is the Surface Tension of Water + R * T * Sum1
		SurfaceTension = SurfaceTensionOfWater(InParticle%Location(1),InParticle%Location(2),InParticle%Location(3)) &
	                 + Rstar*InParticle%Temperature*Sum1
	
		!Now, add in effect of aqueous organics using formula in
		! Jacobson, "Fund. of Atm. Modeling", 2nd ed., p. 534
 		SurfaceTension = SurfaceTension - 0.0187*InParticle%Temperature*LOG(1+628.14*MolalityAqCarbon(InParticle))
	
	ELSE !Use formulas from Jacobson, "Fund. of Atm. Modeling", 2nd ed., p. 534		
		SurfaceTension = SurfaceTensionOfWater(InParticle%Location(1),InParticle%Location(2),InParticle%Location(3)) 
		!Ions
		TotalIonMolality = 0.
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			TotalIonMolality = TotalIonMolality + Molality(I, InParticle)
		END DO
		SurfaceTension = SurfaceTension + 1.7*TotalIonMolality
		!Organics
		SurfaceTension = SurfaceTension - 0.0187*InParticle%Temperature*LOG(1+628.14*MolalityAqCarbon(InParticle))

	END IF

	!!Force a minimum Surface tension
	!!Use the surface tension of Benzene at 298 K 
	!!(Also used in OrgCurvatureCorrection as the surface tension		
	!!See Seinfeld and Pandis (1998), Table 10.2, p. 558
	IF(SurfaceTension .LT. 28.21) SurfaceTension = 28.21

	!! Units are dyn / cm
	InParticle%SurfaceTension = SurfaceTension

	RETURN
END FUNCTION SurfaceTension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the PARTICLE REYNOLDS NUMBER !!
!!										  !!!!!!!!!!!!!!!!!
!! Uses a three-regime (slip, continuum, and transition) !!
!! from Beard's parameterization. (see Jacobson, 1997)   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tested against Fig. 16.3 in Jacobson, 1999 and it     !!
!! matches well.										 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ReynoldsNumber (InParticle)

	USE GridPointFields, ONLY : GetAirDensity,GetDynamicViscosityOfAir
	USE ModelParameters, ONLY : g, grams, cm
	IMPLICIT NONE

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Type local variables
	REAL*8 :: Dv, Cs, AirDensity, X, Y

	Dv		   = GetDynamicViscosityOfAir(InParticle%Location(1),InParticle%Location(2),InParticle%Location(3))
	Cs         = CunninghamSlipCorrectionFactor(InParticle)
	AirDensity = GetAirDensity(InParticle%Location(1),InParticle%Location(2),InParticle%Location(3))

	!! Estimate the Reynolds Number and then Calculate the Real Value, 
	!! this is the estimated terminal velocity
	ReynoldsNumber = 2*InParticle%EffectiveRadius*InParticle%EffectiveRadius*(InParticle%ParticleDensity-AirDensity)*g*Cs/(9*DV)

	!! The Estimated Reynolds Number
	ReynoldsNumber = 2.*InParticle%EffectiveRadius*AirDensity*ReynoldsNumber/DV 

	!! The temporary is the Reynolds Number unless Re > 0.01
	IF (ReynoldsNumber > 0.01 .and. ReynoldsNumber < 300) THEN

		X = dLog((InParticle%ParticleDensity-AirDensity)*AirDensity*InParticle%EffectiveRadius**3.*32.*g/(3.*DV*DV))

		!! The Values here taken from Prupp and Klett, p.417 and originate from Beard's paper (1976)
		Y = -.318657d1+(0.992696d0)*X-(0.153193d-2)*X**2.-(0.987059d-3)*X**3.	&
			-(0.578878d-3)*X**4.+(0.855176d-4)*X**5.-(0.327815d-5)*X**6.

		ReynoldsNumber = Cs*dExp(Y)

	ENDIF

	IF (ReynoldsNumber .GE. 300.) THEN

		!! Use X to signify Np**1/6
		X = (InParticle%SurfaceTension**3.*AirDensity*AirDensity/(DV**4. * (InParticle%ParticleDensity-AirDensity) * g))**(1./6.)

		Y = dLog(16./3.*X*InParticle%EffectiveRadius*InParticle%EffectiveRadius*	&
				 (InParticle%ParticleDensity-AirDensity)*g/InParticle%SurfaceTension)

		ReynoldsNumber = Cs*X*dExp(-5.00015+5.23778*Y-2.04914*Y**2.+0.475294*Y**3.-0.0542819*Y**4.+0.00238449*Y**5.)

	END IF

	RETURN
END FUNCTION ReynoldsNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The correction to pressure from curvature of the organic phase.						   !!
!! It is also known as the "Kelvin Effect"												   !!
!! This uses approximate values for molecular mass and surface tension
!! of the organic phase. Since the curvature correction is mainly for
!! water above RH = 98%, this should be of minor importance.																					   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION OrgCurvatureCorrection (InParticle) 

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : Rstar, micron, ThermoBulkMode, IgnoreCorrectionstoDiffusivity
	USE GridPointFields, ONLY : SurfaceTensionofWater

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER  :: InParticle
	REAL*8 :: TotalMoles
	INTEGER :: I

	!!Type Internal variables
	REAL*8 :: TotalDens, SurfTens, Radius
	
	!Ignore this for bulk mode or no insoluble species
	IF (ThermoBulkMode .OR. IgnoreCorrectionstoDiffusivity) THEN
		OrgCurvatureCorrection = 1.

	ELSE
		TotalDens = InParticle%InsolubleDensity	
		!!Assume a molecular mass of 200 g/mol
		OrgCurvatureCorrection = 2.*200.0/(InParticle%Temperature*Rstar*InParticle%InsolubleDensity)

		!Keep radius above 1 nm for calculation
		IF(InParticle%InsolubleRadius .LE. 0.001*micron) THEN
			Radius = 0.001*micron
		ELSE
			Radius = InParticle%InsolubleRadius
		END IF
		
		!!Assume surface Tension of Benzene, Seinfeld and Pandis, p. 522
		SurfTens = 28.21
		OrgCurvatureCorrection = EXP(OrgCurvatureCorrection*SurfTens/Radius)

	END IF

	!WRITE(*,*) "OrgCurvature Okay"
	RETURN
	
END FUNCTION OrgCurvatureCorrection

!Molality of carbon (mol/kg H2O) in aqueous phase for surface tension calculation (Matt Alvarado, 6/2006)
REAL*8 FUNCTION MolalityAqCarbon(InParticle)

	USE Chemistry,				ONLY :	AqOrgNumCarbon
	USE ModelParameters
	
	IMPLICIT NONE

	REAL*8  :: Total
	INTEGER :: I
	TYPE(PARTICLE),POINTER :: InParticle

	
	Total = 0.0
	
	DO I = 1, HowManyAqOrgChems
		Total = Total + InParticle%AqOrgChems(I)*AqOrgNumCarbon(I)
	END DO
	
	MolalityAqCarbon = Total / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

	RETURN

END FUNCTION MolalityAqCarbon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Mass of aqueous phase (g/particle) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION AqueousSolutionMass (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro
	USE Chemistry,			 ONLY : AqMolecularMass, HowManyAqChems,	&
									HowManyAqAnions, HowManyAqCations,	&
									AqCationMass, AqAnionMass, &
									HowManyAqOrgChems, AqOrgMolecularMass, &
									HowManyAqEqReactions, AqEquilibriaList

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I, C

	AqueousSolutionMass = 0. 

	!Water
	AqueousSolutionMass   = AqueousSolutionMass + AqMolecularMass(1)*InParticle%AqChems(1)

	!Non-solid Electrolytes
	DO I = 2, HowManyAqChems
		!Check for non-solid electrolyte (NH3, CO2, etc), and add in.
		DO C = 1, HowManyAqEqReactions
			IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				AqueousSolutionMass = AqueousSolutionMass + AqMolecularMass(I)*InParticle%AqChems(I)
			END IF
		END DO					 
	END DO

	DO I = 1, HowManyAqCations
		AqueousSolutionMass   = AqueousSolutionMass + AqCationMass(I)*InParticle%AqChems(HowManyAqChems+I)
	END DO

	DO I = 1, HowManyAqAnions
		AqueousSolutionMass   = AqueousSolutionMass + AqAnionMass(I)*InParticle%AqChems(HowManyAqChems+HowManyAqCations+I)
	END DO

	DO I = 1, HowManyAqOrgChems
		AqueousSolutionMass   = AqueousSolutionMass + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
	END DO

	RETURN
END FUNCTION AqueousSolutionMass

!! This subroutine calculates the volume-average refractive indices for
!! the reflective shell and then calls DMILAY for the efficienies and 
!! assymetry parameter
SUBROUTINE ShellRefIndAndRad(InParticle)
! 08/30/2010 Changed for the Photolyis Wavebands MJA
	
	USE ModelParameters
	USE Chemistry, ONLY : AqRefractiveIndex, OrgRefractiveIndex, &
						AqOrgRefractiveIndex, AqMolecularMass, &
						SolidSaltDensity, AqPhaseChemicalNames, &
						OrgMolecularMass, OrgDensity, &
						AqOrgMolecularMass, AqOrgDensity

	USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, ERROR
	
	IMPLICIT NONE
	
	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!!External Subroutines ..
    EXTERNAL  DMiLay

	!! Internal Variables
	INTEGER, PARAMETER  :: MaxAng = 100
 	INTEGER				:: I, JJ, NumAng, ib
	REAL*8				:: SingleParticleVolume, VolumeConc, Factor,xxx
	REAL*8				:: Er, Ei, SumEr, SumEi, Wavelength, Sum
	REAL*8				:: SolnMolWeight, MolalRefraction, MolVol, TotalVolConc, Rshell, &
						   ShellVolumeConc
	REAL*8, ALLOCATABLE	:: AqVolConc(:), &
							   OrgVolConc(:), AqOrgVolConc(:), MolFrac(:)
	COMPLEX				:: SolarShellInd, SolarCoreInd, IRShellInd, IRCoreInd
	REAL*8				:: Angles(MaxAng), CosAng(MaxAng), M1( MaxAng, 2), &
							M2( MaxAng, 2 ), S21( MaxAng, 2 ), &
							D21( MaxAng, 2 ), SolutionRefracIndex(2), SolutionVolConc

	!! BUILD N ARRAY for Species Volume Concentrations
	ALLOCATE (AqVolConc(HowManyAqChems+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate AqVolConc at size "	&
						  //TRIM(INT2STR(HowManyAqChems+HowManyAqCations+HowManyAqAnions)))

	!! BUILD N ARRAY for Species Volume Concentrations
	ALLOCATE (OrgVolConc(HowManyOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate OrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyOrgChems)))

	!! BUILD N ARRAY for Species Volume Concentrations
	ALLOCATE (AqOrgVolConc(HowManyAqOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate AqOrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))

	!! BUILD N ARRAY for Species Volume Concentrations
	ALLOCATE (MolFrac(1+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate MolFrac at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))
	
	!WRITE(*,*) "Number of Particles: ", InParticle%numberofparticles 
	!WRITE(*,*) "Size (um): ", InParticle%effectiveradius*1.0e4
	IF (InParticle%numberofparticles .GT. 1.0e-6) THEN
		SingleParticleVolume = (4.0*Pi/3.0)*(InParticle%effectiveradius)**3 !cm3
		
		VolumeConc = SingleParticleVolume*InParticle%numberofparticles !cm3/cm3
		
		!Convert aqueous concentrations to cm3 solution/cm3 air
		
		!Water
		AqVolConc(1)  = InParticle%AqChems(1) * InParticle%NumberOfParticles*AqMolecularMass(1)/InParticle%SolutionDensity
		SolutionVolConc = AqVolConc(1)
		TotalVolConc = AqVolConc(1)
		
				
		!Solid Salts
		DO I = 2, HowManyAqChems
			AqVolConc(I)  = InParticle%AqChems(I) * InParticle%NumberOfParticles*AqMolecularMass(I)/SolidSaltDensity(I)
			TotalVolConc = TotalVolConc + AqVolConc(I)
			!WRITE(*,*) I, "TotalVolConc", TotalVolConc
			IF(AqVolConc(I) .LT. 0.0) WRITE(*,*) "Negative AqVolume Conc!"	
			IF(AqVolConc(I) .GT. 1.0e30) WRITE(*,*) "Error", AqPhaseChemicalNames(I), AqVolConc(I)
		END DO

		!Inorganic Ions
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(I)  = InParticle%AqChems(I) * InParticle%NumberOfParticles*AqMolecularMass(I)/InParticle%SolutionDensity
			TotalVolConc = TotalVolConc + AqVolConc(I)
			!WRITE(*,*) I, "TotalVolConc", TotalVolConc
			SolutionVolConc = SolutionVolConc + AqVolConc(I)
			IF(AqVolConc(I) .LT. 0.0) WRITE(*,*) "Negative AqVolume Conc!"	
		END DO

		!Convert org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyOrgChems
			OrgVolConc(I)  = InParticle%OrgChems(I) * InParticle%NumberOfParticles*OrgMolecularMass(I)/OrgDensity(I)
			TotalVolConc = TotalVolConc + OrgVolConc(I)
			!WRITE(*,*) I, "TotalVolConc", TotalVolConc
			IF(OrgVolConc(I) .LT. 0.0) WRITE(*,*) "Negative OrgVolume Conc!"	
		END DO
	
		!Convert hydrophilic org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(I)  = InParticle%AqOrgChems(I) * InParticle%NumberOfParticles*AqOrgMolecularMass(I)/AqOrgDensity(I)
			TotalVolConc = TotalVolConc + AqOrgVolConc(I)
			!WRITE(*,*) I, "TotalVolConc", TotalVolConc
			IF(AqOrgVolConc(I) .LT. 0.0) WRITE(*,*) "Negative AqOrgVolume Conc!"	
		END DO

		!Calculate shell volume conc (Total - BC)
		ShellVolumeConc = VolumeConc - OrgVolConc(1)
			
	ELSE !No particles, use dummy numbers
		InParticle%ExtEff = 1.0
		InParticle%ScaEff = 1.0
		InParticle%BackScaEff = 1.0
		InParticle%AssymParam = 1.0
        RETURN          
	END IF

	!Calculate solar refractive index
	!Use the volume average dielectric constant mixing rule
	!See Jacobson, "Fund. of Atm. Modeling", 2nd ed, p. 310
	SumEr = 0.
	SumEi = 0.
	
	!Calculate average refractive index for ionic solution
	!Following Tang et. al., JGR V.102, D19, pp. 23269-23275, 1997. 
        !and Tang, JGR V.102, D2, pp. 1883-1893, 1997.)
		!Step 1: Calculate mole fractions in solution
		Sum = InParticle%AqChems(1)
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			Sum = Sum + InParticle%AqChems(I)
		END DO
		MolFrac(1) = InParticle%AqChems(1)/Sum
		DO I = 2, 1+HowManyAqCations+HowManyAqAnions
			MolFrac(I) = InParticle%AqChems(HowManyAqChems+I-1)/Sum
		END DO
		
		!Step 2: Calculate average molecular weight	
		SolnMolWeight = MolFrac(1)*AqMolecularMass(1)
		DO I = 2, 1+HowManyAqCations+HowManyAqAnions
			SolnMolWeight = SolnMolWeight + MolFrac(I)*AqMolecularMass(HowManyAqChems+I-1)
		END DO
		
		!Step 3: Calculate average molal refraction
		MolalRefraction = MolFrac(1)*3.717 !See Tang et. al, 1997
		DO I = 2, 1+HowManyAqCations+HowManyAqAnions
			MolalRefraction = MolalRefraction + MolFrac(I)*AqRefractiveIndex(HowManyAqChems+I-1,1)
		END DO
		
		!Step 4: Calculate molal volume
		MolVol = SolnMolWeight/InParticle%InorgSolnDensity

		!Step 5: Calculate Refractive Index
		xxx = MolalRefraction/MolVol
		SolutionRefracIndex(1) = SQRT((2.*xxx+1.)/(1.-xxx))
		SolutionRefracIndex(2) = 1.0e-9 !Imaginary index for water
	
		!Step 6: Calculate Er and Ei
		Er = SolutionRefracIndex(1)**2 - SolutionRefracIndex(2)**2
		Ei = 2*SolutionRefracIndex(1)*SolutionRefracIndex(2)
		SumEr = SumEr + (SolutionVolConc/ShellVolumeConc)*Er
		SumEi = SumEi + (SolutionVolConc/ShellVolumeConc)*Ei
	
	!Solid salts
	DO I = 2, HowManyAqChems
		Er = AqRefractiveIndex(I,1)**2 - AqRefractiveIndex(I,2)**2
		Ei = 2*AqRefractiveIndex(I,1)*AqRefractiveIndex(I,2)
		SumEr = SumEr + (AqVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqVolConc(I)/ShellVolumeConc)*Ei
	END DO
	
	!Organics, skipping BC
	DO I = 2, HowManyOrgChems
		Er = OrgRefractiveIndex(I,1)**2 - OrgRefractiveIndex(I,2)**2
		Ei = 2*OrgRefractiveIndex(I,1)*OrgRefractiveIndex(I,2)
		SumEr = SumEr + (OrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (OrgVolConc(I)/ShellVolumeConc)*Ei
	END DO

	!Hydrophilic Organics
	DO I = 1, HowManyAqOrgChems
		Er = AqOrgRefractiveIndex(I,1)**2 - AqOrgRefractiveIndex(I,2)**2
		Ei = 2*AqOrgRefractiveIndex(I,1)*AqOrgRefractiveIndex(I,2)
		SumEr = SumEr + (AqOrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqOrgVolConc(I)/ShellVolumeConc)*Ei
	END DO

	Factor = SQRT(SumEr**2 + SumEi**2)
	InParticle%ShellRealRefracInd(1) = SQRT((Factor + SumEr)/2.0)
	InParticle%ShellImagRefracInd(1) = SQRT((Factor - SumEr)/2.0)
	
	
	!Calculate Core radius
	InParticle%AbsCoreRad = OrgVolConc(1)/InParticle%NumberofParticles
	InParticle%AbsCoreRad = (0.75*InParticle%AbsCoreRad/Pi)**(1.0/3.0) !In cm
	IF(InParticle%AbsCoreRad .GT. InParticle%EffectiveRadius) InParticle%AbsCoreRad = 0.9999*InParticle%EffectiveRadius
	!WRITE(*,*) "Rshell by TotalVolConc: ", (0.75*TotalVolConc/InParticle%NumberofParticles/Pi)**(1.0/3.0)
	!Rshell = (0.75*TotalVolConc/InParticle%NumberofParticles/Pi)**(1.0/3.0)
	!WRITE(*,*) "Rcore: ", InParticle%AbsCoreRad
	!WRITE(*,*) "Rshell: ", InParticle%EffectiveRadius
									

	SolarShellInd = CMPLX(InParticle%ShellRealRefracInd(1), -1.0*InParticle%ShellImagRefracInd(1))
	SolarCoreInd = CMPLX(OrgRefractiveIndex(1,1), -1.0*OrgRefractiveIndex(1,2))
	!WRITE(*,*) "Shell Ind: ", SolarShellInd
	!WRITE(*,*) "Core Ind: ",  SolarCoreInd
	NumAng = 0.

	!Calculate IR refractive index
	!Use the volume average dielectric constant mixing rule
	!See Jacobson, "Fund. of Atm. Modeling", 2nd ed, p. 310
	SumEr = 0.
	SumEi = 0.
	
	!Water and solid salts
	DO I = 1, HowManyAqChems
		Er = AqRefractiveIndex(I,3)**2 - AqRefractiveIndex(I,4)**2
		Ei = 2*AqRefractiveIndex(I,3)*AqRefractiveIndex(I,4)
		SumEr = SumEr + (AqVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqVolConc(I)/ShellVolumeConc)*Ei
	END DO

	!Assume ions the same as H2O
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		Er = AqRefractiveIndex(1,3)**2 - AqRefractiveIndex(1,4)**2
		Ei = 2*AqRefractiveIndex(1,3)*AqRefractiveIndex(1,4)
		SumEr = SumEr + (AqVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqVolConc(I)/ShellVolumeConc)*Ei
	END DO
	
	!Organics, skipping BC
	DO I = 2, HowManyOrgChems
		Er = OrgRefractiveIndex(I,3)**2 - OrgRefractiveIndex(I,4)**2
		Ei = 2*OrgRefractiveIndex(I,3)*OrgRefractiveIndex(I,4)
		SumEr = SumEr + (OrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (OrgVolConc(I)/ShellVolumeConc)*Ei
	END DO

	!Hydrophilic Organics
	DO I = 1, HowManyAqOrgChems
		Er = AqOrgRefractiveIndex(I,3)**2 - AqOrgRefractiveIndex(I,4)**2
		Ei = 2*AqOrgRefractiveIndex(I,3)*AqOrgRefractiveIndex(I,4)
		SumEr = SumEr + (AqOrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqOrgVolConc(I)/ShellVolumeConc)*Ei
	END DO

	Factor = SQRT(SumEr**2 + SumEi**2)
	InParticle%ShellRealRefracInd(2) = SQRT((Factor + SumEr)/2.0)
	InParticle%ShellImagRefracInd(2) = SQRT((Factor - SumEr)/2.0)

	IRShellInd = CMPLX(InParticle%ShellRealRefracInd(2), -1.0*InParticle%ShellImagRefracInd(2))
	IRCoreInd = CMPLX(OrgRefractiveIndex(1,3), -1.0*OrgRefractiveIndex(1,4))
	NumAng = 0.
		
	!Call DMiLay for all 18 wavelength bins used in CRM6, using proxy wavelength
	!(usually midpoint of band range)
	!In this radiation scheme,  six  and  12 bands are selected for solar 
        ! and thermal IR regions, respectively. The spectral division is below: 
        ! Band Range              Refractive Index  Proxy Wavelength
	! 0.2 - 0.7 um              Solar (550 nm)      550 nm
	! 0.7 - 1.3 um              Solar (550 nm)      1.0 um
	! 1.3 - 1.9 um	            Solar (550 nm)      1.6 um
	! 1.9 - 2.5 um              Solar (550 nm)      2.2 um
	! 2.5 -3.5 um               Solar (550 nm)      3.0 um
	! 3.5 - 4.0 um              Solar (550 nm)      3.75 um  
	! 2200 - 1900 cm**-1,       IR (10 um)		4.79 um
	! 1900 - 1700 cm**-1,	    IR (10 um)		5.45 um 
	! 1700 -1400  cm**-1,	    IR (10 um)		6.51 um
	!  1400 - 1250 cm**-1,	    IR (10 um)		7.57 um
	!  1250 - 1100 cm**-1,	    IR (10 um)		8.55 um
	! 1100 - 980 cm**-1,	    IR (10 um)		9.65 um
        ! 980 - 800 cm**-1,	    IR (10 um)	       11.35 um
	!  800 - 670 cm**-1,	    IR (10 um)	       13.71 um
	!  670 - 540 cm**-1,	    IR (10 um)	       16.72 um
	! 540 - 400 cm**-1,	    IR (10 um)	       21.76 um
	!  400 - 280 cm**-1,	    IR (10 um) 	       30.36 um
	!  280 - 0 cm**-1,	    IR (10 um)	       71.43 um
	
	!DO ib = 1, 18
	DO ib = 1, 62 !For photolysis between 290 nm and 410 nm, plus 550 and 700 nm
		
                IF (ib .LE. 60) THEN 
                       Wavelength = 2*PI/((290+(i-1)*2)*nm)
                ELSE IF (ib .EQ. 61) THEN
                       Wavelength = 2*PI/(550*nm)
                ELSE 
                       Wavelength = 2*PI/(700*nm)
                END IF

                !Choose correct proxy wavelength for CRM6
		!SELECT CASE(ib)
		!	CASE(1)
		!		Wavelength = 2*PI/(550*nm)
		!	CASE(2)
		!		Wavelength = 2*PI/(1.0*micron)
		!	CASE(3)
		!		Wavelength = 2*PI/(1.6*micron)		
		!	CASE(4)
		!		Wavelength = 2*PI/(2.2*micron)
		!	CASE(5)
		!		Wavelength = 2*PI/(3.0*micron)				
		!	CASE(6)
		!		Wavelength = 2*PI/(3.75*micron)
		!	CASE(7)
		!		Wavelength = 2*PI/(4.79*micron)
		!	CASE(8)
		!		Wavelength = 2*PI/(5.45*micron)		
		!	CASE(9)
		!		Wavelength = 2*PI/(6.51*micron)
		!	CASE(10)
		!		Wavelength = 2*PI/(7.57*micron)				
		!	CASE(11)
		!		Wavelength = 2*PI/(8.55*micron)
		!	CASE(12)
		!		Wavelength = 2*PI/(9.65*micron)				
		!	CASE(13)
	       	!               Wavelength = 2*PI/(11.35*micron)
		!	CASE(14)
		!		Wavelength = 2*PI/(13.71*micron)				
		!	CASE(15)
		!		Wavelength = 2*PI/(16.72*micron)				
		!	CASE(16)
		!		Wavelength = 2*PI/(21.76*micron)
		!	CASE(17)
		!		Wavelength = 2*PI/(30.36*micron)				
		!	CASE(18)
		!		Wavelength = 2*PI/(71.43*micron)
		!	END SELECT

		!IF (ib .LT. 7) THEN !Solar bands
				
			CALL DMiLay( InParticle%AbsCoreRad, &
				 InParticle%EffectiveRadius, &
				 Wavelength, &
				 SolarShellInd, & 
				 SolarCoreInd, &
				  CosAng, &
				  NumAng, &
				  InParticle%ExtEff(ib), &
				  InParticle%ScaEff(ib), &
				  InParticle%BackScaEff(ib), &
				  InParticle%AssymParam(ib), & 
                        M1, M2, S21, D21, MAXANG )

		!ELSE !IR bands
		!	CALL DMiLay( InParticle%AbsCoreRad, &
		!		 InParticle%EffectiveRadius, &
		!		 Wavelength, &
		!		 IRShellInd, & 
		!		 IRCoreInd, &
		!		  CosAng, &
		!		  NumAng, &
		!		  InParticle%ExtEff(ib), &
		!		  InParticle%ScaEff(ib), &
		!		  InParticle%BackScaEff(ib), &
		!		  InParticle%AssymParam(ib), & 
                !        M1, M2, S21, D21, MAXANG )
		!END IF
			
		!DMiLay gives the assymetry parameter times the scattering efficiency
		!So divide by scattering efficiency here (MJA, 10/15/07)
		InParticle%AssymParam(ib) = InParticle%AssymParam(ib)/InParticle%ScaEff(ib)

	END DO

	RETURN
END SUBROUTINE ShellRefIndAndRad

!Calculates aerosol optical depth for a grid point at 550 nm.
!Can only be used in inversion mode! 
REAL*8 FUNCTION AerosolOpticalDepth (X, Y, Z)
! 08/30/2010 Changed for photolysis radiation bins MJA
	
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE

	!External Variables
	INTEGER :: X, Y, Z
	
	Type(Particle),POINTER :: Particle

	AerosolOpticalDepth = 0.
	
	Particle => Particles(X,Y,Z)%First

10	IF (ASSOCIATED(Particle)) THEN

		AerosolOpticalDepth = AerosolOpticalDepth + Particle%NumberofParticles*(PI*Particle%EffectiveRadius**2)*Particle%ExtEff(61)
		
		IF (ASSOCIATED(Particle%Next)) THEN
			Particle => Particle%Next
			GOTO 10
		END IF
	
	END IF
	
	AerosolOpticalDepth = AerosolOpticalDepth*InversionHeight
	
	RETURN
END FUNCTION AerosolOpticalDepth

!Calculates total aerosol extinction coefficient (m**-1) 
!average single scattering albedo, and 
!average assymetry factor 
!for a grid point at 18 radiative bands.
SUBROUTINE AerosolOptProp (X, Y, Z, ExtCoeff, SingScat, Assym)
! 08/30/2010 Changed for the photolysis wavebands MJA
	
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE

	!Output Variables
	REAL*8, DIMENSION(18) :: ExtCoeff, SingScat, Assym
	
	!External Variables
	INTEGER :: X, Y, Z
	
	!Internal variables
	INTEGER :: ib
	REAL*8, DIMENSION(18) :: ScatCoeff
	REAL*8 :: AssymSum, xxx
	Type(Particle),POINTER :: Particle

	!DO ib = 1, 18
        DO ib = 1, 62
		ExtCoeff(ib) = 0.0
		ScatCoeff(ib) = 0.0
		AssymSum = 0.0
		Particle => Particles(X,Y,Z)%First

	10	IF (ASSOCIATED(Particle)) THEN
			xxx = Particle%NumberofParticles*(PI*Particle%EffectiveRadius**2)
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*Particle%ExtEff(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*Particle%ScaEff(ib)
			AssymSum = AssymSum + xxx*Particle%ScaEff(ib)*Particle%AssymParam(ib)
					
			IF (ASSOCIATED(Particle%Next)) THEN
				Particle => Particle%Next
				GOTO 10
			END IF
	
		END IF

		!Formula from p.327, eqn 9.103 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		SingScat(ib) = ScatCoeff(ib)/ExtCoeff(ib)
		IF(ib .EQ. 61) WRITE(*,*) "SingScat: ", SingScat(1)
		
		!Formula from p.329, eqn 9.117 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		Assym(ib) = AssymSum/ScatCoeff(ib)
		IF(ib .EQ. 61) WRITE(*,*) "Assym: ", Assym(1)

		!Convert ExtCoeff from 1/cm to 1/m
		ExtCoeff(ib) = 100*ExtCoeff(ib)
		IF(ib .EQ. 61) WRITE(*,*) "ExtCoeff: ", ExtCoeff(1)

	END DO
	
	RETURN
END SUBROUTINE AerosolOptProp
