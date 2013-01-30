!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! CondensationRelatedFunctions.h
!! Contains many functions needed by the numerical integration routines,
!! and includes the main equilibration routine, as well as many equilibration 
!! subroutines (water, sulfate, gas-solid reactions, internal equilibrium).

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 08/30  2006   Matt Alvarado     EquilibrateInternallyatGridPoint now stops!!
!!				checking for inorg eq after 5 iterations     !!
!! 08/31  2006	 Matt Alvarado     Edited EquilibrateDissolutionReaction     !!
!!				   to return if NH3(aq) approaches 0.	     !!
!! 09/07  2006   Matt Alvarado     1. Added RecalculateRadius to	     !!
!!					EquilibrateInternallyatGridPoint     !!
!!				   2. Created RidderMethod root finder	     !!
!! 09/08  2006   Matt Alvarado    Created  EquilibriumWaterContentAll	     !!
!!				UpdateWaterActivityAll and WaterResidualAll  !!
!! 08/01  2007   Matt Alvarado    Removed low water checks from              !!
!!				  EquilibrateGridPointAll, since they        !!
!!				  sometimes lead to unphysical answers.      !!
!! 09/28  2007   Matt Alvarado    Created EquilibrateInorgGridPoint          !!
!! 10/01  2007   Matt Alvarado    Increased "already in eq" error tolerance  !!
!!				  in EquilibrateGasSolidReactions
!!				  to 1%
!!		Protected EquilibrateGasSolidReactions against 0 particles
!!
!! 10/03  2007   Matt Alvarado     Added line to EquilibrateInorgGridPoint 
!!                                 to skip mostly empty bins
!!				Added same to EquilibrateGridPointAll
!!				Reduced gas-solid error tolerance to 0.5%
!! 10/04  2007   Matt Alvarado     In EquilibrateInorgGridPoint, go to lower !!
!!                          error tolerance after 10 iterations instead of 25!!
!! 10/05  2007   Matt Alvarado     Changed EquilibrateDissolutionReaction to !!
!!                                  just exit after trying to equilibrate    !!
!!                                 for 2000 iterations                       !!
!!		Set EquilibrateGasSolid to just exit after 1000 iterations   !!
!! 10/11  2007   Matt Alvarado     Changed EquilibrateDissolutionReaction    !!
!!                                 to just exit *with returntype=2* after    !!
!!                                  trying to equilibrate for 2000 iterations!!
!! 10/14  2010   Matt Alvarado     Removed "OPTIONAL" statements             !!
!! 02/17  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/29  2012   Matt Alvarado     Changed EquilibrateInternallyatGridPoint  !!
!!                                   to skip dry (BC only) particles.        !!
!! 01/25  2013   Matt Alvarado     Changed EquilibrateGridPointAll to not use!!
!!                                  CACM SOA functions UpdateGasPhaseOrganics!!
!!                                  and CalculateSurrogateGasConc            !!
!! 01/28  2013   Matt Alvarado    Removed  EquilibriumWaterContentAll	     !!
!!				UpdateWaterActivityAll and WaterResidualAll  !!
!!                              Also removed call to EquilibriumWaterContent !!
!!                                 in EquilibrateDissolutionReaction.        !!
!! 01/29  2013   Matt Alvarado    Removed Subroutine DryParticleQuickly
!! 01/30  2013   Matt Alvarado    Added checks to EquilibrateDissolutionReaction
!!                                to keep ions above 0.
!!                                Added ElectrolyteEq step to 
!!                                 EquilibrateInorgGridPoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		     !!
!! 1. FUNCTION GetThermalVelocity (ChemIndex, ChemName)
!! 2. FUNCTION CondensationRateCoefficients (ReactionNumb, NumbLinks)
!! 3. FUNCTION HenrysLaw (Temperature, DissolutionRxn)
!! 4. FUNCTION SatVapPressRatio (InParticle , ChemIndex )
!! 5. SUBROUTINE EquilibrateDissolutionReaction (GasPhaseBurden, WhichRxn, & !!
!!                 InParticle, Temperature, ReturnType, OptErrorTolerance, & !!
!!		    WaterSatBurden, GasPhaseBurden2)                         !!
!! 6. FUNCTION DissolutionEquilibriumConstantsRatio (GasPhaseBurden,       & !!
!!                                    WhichRxn, Temperature, InParticle,   & !!
!!					GasPhaseBurden2, WaterSatBurden)     !!
!! 7. SUBROUTINE EquilibrateSulfate()
!! 8. SUBROUTINE EquilibrateGasSolidReactions()
!! 9. SUBROUTINE EquilibrateGridPointAll
!!10. SUBROUTINE EqAllWater
!!11. SUBROUTINE EquilibrateInternallyatGridPoint
!!12. FUNCTION RidderMethod
!!13. SUBROUTINE EquilibrateInorgGridPoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The THERMAL VELOCITY of SOME NON-WATER CHEMICAL !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetThermalVelocity (ChemIndex)

	USE ModelParameters,	 ONLY : ThmVelFac, cm, Avogadro
	USE Chemistry,	ONLY : GasMolecularMass, FindChem, HowManyAqEqReactions
	USE InfrastructuralCode, ONLY : WARN, ERROR, INT2STR
	USE GridPointFields,     ONLY : GetTemp

	IMPLICIT NONE

	!! External Variables
	INTEGER		:: ChemIndex

	!! Internal Variables
	INTEGER :: I
 
        I = ChemIndex

	GetThermalVelocity = ThmVelFac * SQRT(Avogadro * GetTemp() / &
                              GasMolecularMass(I))

	RETURN
END FUNCTION GetThermalVelocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate MASS FLUX CONDENSATION RATE COEFFICIENT condensation !!
!! and dissolution.  This is separated from CondensationRates()   !!
!! to make it easier to calculate both the Jacobian and the ODEs  !!
!! without having to recalculate all of this at each step.        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION CondensationRateCoefficients (ReactionNumb, NumbLinks)

	USE Chemistry,		 ONLY : EnthalpyOfVaporization,		&
					GasMolecularmass,		&
					MolecularDiffusionCoefficient

	USE GridPointFields,     ONLY : GetAirDensity,		        &
					GetCpm,				&
					GetDynamicViscosityOfAir,       &
					GetTemp,			&
					MolecularThermalConductivity,	&
					GetSatVapConcentration

	USE InfrastructuralCode, ONLY : ERROR, INT2STR

	USE Aerosols,		 ONLY : Particle, Particles,		&
					ParticleKnudsenForEnergy,	&
					ReynoldsNumber

	USE ModelParameters,     ONLY : AirMolecMass,		        &
				        Avogadro,		        &
					Pi,				&
					ThermalAccommodationCoeff,	&
					Rstar,				&
					ChemScale,			&
					ThermoBulkMode,			&
					IgnoreCorrectionstoDiffusivity

    IMPLICIT NONE

	!! Define input and output variables
	INTEGER	:: ReactionNumb, NumbLinks
	REAL*8	:: CondensationRateCoefficients(NumbLinks)

	TYPE(Particle),POINTER :: Current
	INTEGER ::  I, ChemIndex
	REAL*8  ::      AirDensity,		&
			AirTemp,		&
			Cpm,			&	! Moist specific heat
			Dv,			&	! Diffusion Coeffienct
			DvPrime,		&
			DynVisc,		&
			EnergyKnud,		&	
			GasKnud,		&
			GasThermalVel,	        &
			KappaD,			&	! Thermal conductivity 
			KappaDPrime,	        &
			ReynoldsNumb,	        &
			SchmidtNumber,	        &
			X,			&
			GasMeanFreePath,        &
			SatWaterVapConc
	
	LOGICAL :: Scaffolding = .FALSE.

	!! Get the primary condensing chemical
	ChemIndex = DissolutionData(ReactionNumb,1)

	!! We don't allow yet the combination of two gas phase
	!! species dissolving into and out of aerosol, but 
	!! the input file can parse them.
	!! barf here if one is specified.
	IF (DissolutionData(ReactionNumb,11) .EQ. 5) CALL ERROR ("Haven't yet implemented two-species dissolution")

	!! Retrieve or Calculate the Variables that are the 
	!! same for all aerosol
	AirTemp		= GetTemp()
	AirDensity	= GetAirDensity()
	Dv		= MolecularDiffusionCoefficient(ChemIndex)
	DynVisc		= GetDynamicViscosityOfAir()
	GasThermalVel   = GetThermalVelocity(ChemIndex)
	SatWaterVapConc	= GetSatVapConcentration()


	!! Non water condensation
	IF (ReactionNumb .GT. 0) THEN
		Cpm    = GetCpm()
		KappaD = MolecularThermalConductivity()
	END IF


	!! Constant is 32. / 3. / Pi
	GasMeanFreePath = 3.395305453 * Dv * AirMolecMass / GasThermalVel /&
			(AirMolecMass + GasMolecularMass(ChemIndex) / Avogadro)


	!! Get the first particle
	Current => Particles%First

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Loop over all of the particles !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO I = 1, NumbLinks

	!! Walk the chain and err if can't find the next link
	IF (I .GT. 1) THEN
		IF (.NOT. ASSOCIATED(Current%Next)) &
			CALL ERROR("Miscounted the number of particles in the argument to CondensationRateCoefficients().  This just shouldn't happen.")
		Current => Current%Next
	END IF

	!! If it's a section, make sure there's something in it
	!! Or if the particle can't attract water
	IF ((Current%Sectional .AND. Current%NumberOfParticles .EQ. 0) .OR. Current%Dry) THEN
		CondensationRateCoefficients(I) = 0.
		CYCLE
	END IF

	!! The ones that are not the same for each
	ReynoldsNumb = ReynoldsNumber(Current)    ! defined for a particle by its radius
	EnergyKnud	 = ParticleKnudsenForEnergy(Current)


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Water is not the same other Chems. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! We are calculating k_i for a Dissolution ODE that looks like this: !!
!!                                                                    !!
!!  d c_i          /               c_i   \                            !!
!!  -----  =  k_i |  C_i - S'_i * -----   |                           !!
!!    dt           \                H'   /                            !!
!!								      !!
!! Condensation is similar, as noted below.			      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The Dissolution Rate Coefficient   !!
	!! is this (condensation is corrected !!
	!! further, see below):		      !!
	!!				      !!
	!! k_i = 4 pi r D_i'                  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The MOLECULAR DIFFUSION COEFFICIENT must be corrected !!
	!! for collisional geometry and ventilation.		 !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! The Gas Knudsen Number 
	GasKnud	= GasMeanFreePath / Current%EffectiveRadius 
	SchmidtNumber = DynVisc / AirDensity / Dv

	X = ReynoldsNumb**0.5 * SchmidtNumber**0.3333333333333

	!! Following Pruppacher and Klett, use this X-Factor to calculate the 
	!! VENTILLATION CORRECTION
	IF (X .LE. 1.4) THEN
		DvPrime = Dv*(1.+0.108*X*X)
	ELSE
		DvPrime = Dv*(0.78+0.308*X)
	END IF

	!! Correct Dv away from the continuum solution for 
        !! COLLISIONAL GEOMETRY and STICKING PROBABILITY
	!! Following Jacobson, 1999, (textbook) who follows 
        !! Fuchs and Sutugin 1971 and Pruppacher and Klett 1997
	DvPrime = DvPrime/(1. + ((1.33+0.71/GasKnud)/(1.+1./GasKnud)+ &
			  4*(1.-DissolutionData(ReactionNumb,10))/    &
			  (3*DissolutionData(ReactionNumb,10)))*GasKnud)


	!Ignore corrections if ordered to
	IF(IgnoreCorrectionstoDiffusivity) DvPrime = Dv
	
	!! Assemble the k_i for that section / particle
	CondensationRateCoefficients(I) = Current%EffectiveRadius*4*Pi*DvPrime


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! It is water!  And so we have more work to do... !!
	IF (ChemIndex .EQ. 1) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! For water, the condensation rate is: !!
	!!                                      !!
	!!  d c_w                               !!
	!!  -----  =  k_w C_s,w ( RH - S'_w )   !!
	!!    dt                                !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The THERMAL CONDUCTIVITY OF AIR must also be corrected !!
	!! for collisional geometry and ventilation.			  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Again play this X trick, this tipe the 
	!! parenthetical is the Prandtl Number
	X = ReynoldsNumb**0.5 * (DynVisc * Cpm / KappaD)**0.3333333333333

	!! Following Pruppacher and Klett, use this X-Factor to calculate the 
	!! VENTILLATION CORRECTION

	IF (X .LE. 1.4) THEN
		KappaDPrime = KappaD*(1.+0.108*X*X)
	ELSE
		KappaDPrime = KappaD*(0.78+0.308*X)
	END IF

	!! Correct for the Collisional Geometry
	KappaDPrime = KappaDPrime/(1. + ((1.33+0.71/EnergyKnud)/	&
		(1.+1./EnergyKnud) + 4*(1.-ThermalAccommodationCoeff) /	&
		(3*ThermalAccommodationCoeff))*EnergyKnud)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! But we do find the Condensation Rate Coefficient for H2O:	        !!
!!								        !!
!!			              4 pi r D_w'		        !!
!! k_w = ------------------------------------------------------	        !!
!!		  /  D_w' L_w M_w S_w' C_s,w  / L_w M_w	    \	   \	!!
!!		 |  -----------------------  | --------- - 1 | + 1  | V	!!
!!		  \	    kappa_d'  T	      \  R*  T      /	   /    !!
!!									!!
!!						                        !!
!! is a more complex form because we don't let ourselves		!!
!! neglect L_w, as we did before, and so have a big term in		!!
!! the denominator.							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CondensationRateCoefficients(I) = CondensationRateCoefficients(I) &
	/((DvPrime * EnthalpyOfVaporization(1) * GasMolecularMass(1) *	&
	SatVapPressRatio(Current, 1) *	SatWaterVapConc / chemscale /	&
	Avogadro / KappaDPrime / AirTemp) *	(EnthalpyOfVaporization(1) * &
	GasMolecularMass(1) / ChemScale / Rstar / AirTemp - 1.) + 1.)

	END IF	!! End of water if
	END DO  !! Move to another particle

	RETURN
END FUNCTION  CondensationRateCoefficients


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Henry's Law Coefficient !!
!! This is indexed by the aqueous    !!
!! phase chemical index.	     !!
!! Puts out in Moles / MB, so make   !!
!! sure to convert to the weird      !!
!! pressure units.		     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION HenrysLaw (Temperature, DissolutionRxn)

	USE ModelParameters, ONLY : EqCoeffTref 

	IMPLICIT NONE

	REAL*8  :: Temperature, TempRatio
	INTEGER :: DissolutionRxn

	TempRatio = EqCoeffTref / Temperature


	HenrysLaw = DissolutionData(DissolutionRxn,5)			&
	  	* EXP(DissolutionData(DissolutionRxn,6) * (TempRatio - 1.0) + &
	DissolutionData(DissolutionRxn,7) * (LOG(TempRatio) - TempRatio + 1.))

	
	RETURN
END FUNCTION HenrysLaw


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Correction to the Near-Surface Concentrations  !!
!! for use in the condensation routines.  For water vapor, this !!
!! may either be the kelvin or a more complicated effect.  For	!!
!! other species, it is a simple correction based on the kelvin !!
!! effect.							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION SatVapPressRatio (InParticle , ChemIndex )

	USE Chemistry,	ONLY :	HowManyGasChems,	&
				HowManyAqChems,		&
				HowManyAqCations,	&
				HowManyAqAnions,	&
				EnthalpyOfVaporization, &
				AqPhaseChemicalNames,	&	
				AqCationNames,AqAnionNames


	USE Aerosols,   ONLY : CurvatureCorrection,    &
				FindElectrolyteEquilibrium, &
				recalculateradius,ParticleMass

	USE ModelParameters, ONLY : ThermoBulkMode, RStar,  &
                                    Pi, moles, micron, WaterMolecMass,&
				OutputDeckSubDir, CosContactAngle

	USE InfrastructuralCode, ONLY : ERROR,real2str

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle
	INTEGER				   :: ChemIndex

	!! Internal Variables
	INTEGER :: I
	REAL*8  :: Radiative, XX, YY, dRdV
	real*8  :: ii,jj,kk
	
	!! Thermo Mode doesn't use this correction
	IF (ThermoBulkMode) THEN

		SatVapPressRatio = 1.
		RETURN

	ELSE

		
		IF (ChemIndex .EQ. 1) THEN
		

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! WATER CONDENSATION !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The water activity term is not addressed here,	!!
!! it happens either in the condensation calculation	!!
!! itself (where we use Raoult's Law) or in the		!!
!! equilibration routine (where we use a more compelx	!!
!! form).  We then here choose either a Kohler or	!!
!! a Gorbunov and Hamilton type approach.		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			IF (InParticle%InsolubleRadius .EQ. 0.) THEN	

			!! The Kelvin Effect is abstracted to a function 
			!! and housed in "ParticleAttributes.h"
			    SatVapPressRatio  = CurvatureCorrection(InParticle)
			
			ELSE
				
			!! Use Gorbunov and Hamilton 1998 and 1999 approach
			!! Calculate the dimentionless parameters
				XX = InParticle%InsolubleRadius / &
                                      InParticle%EmbryoRadius
				YY = Sqrt(1. + XX*XX - 2.*XX*CosContactAngle)

				!! If this is the Kohler Example, there may 
                                !! be problems if we use Fletcher's Approx.
				IF (InParticle%InsolubleRadius .EQ. 0. .OR. & 
                                 CosContactAngle+1. .LT. 0.0000001) THEN

					!! The core is not involved in this 
                                        !! case, so make the simple calculation
					dRdV = 0.25/Pi/InParticle%EmbryoRadius/InParticle%EmbryoRadius

				ELSE

					!! This is a problem, meaning either no
                                        !! or negative volume for the embryo
					IF ((ABS(CosContactAngle -1.)<0.0000001) .AND.			&
						(InParticle%InsolubleRadius .GE. InParticle%EffectiveRadius))   &
						CALL ERROR("In SUBROUTINE SatVapPressRatio (), a water condensation related function that " //	&
						           "calculates S' using the Gorbunov and Hamilton (1997) formulation, found there was " // &
								   "negative or zero solution volume, which should never occur.")

					!! See Fletcher 1962 or Gorbunov and Hamilton, 1997 for this statement
					dRdV =	1./ ((Pi*InParticle%InsolubleRadius * InParticle%InsolubleRadius*						&
							((2./XX/XX-																			&
							 (1.-CosContactAngle*CosContactAngle)**2.*XX*XX/YY/YY/YY/YY/YY+ &  
							 2*(1.-CosContactAngle*XX)/YY/XX/XX +                        			&
							 (1.-CosContactAngle*CosContactAngle)/YY/YY/YY *				&
							 (1.-(CosContactAngle-XX)*(1.-CosContactAngle*XX)*XX/YY/YY)))))

				END IF


	!! S' from Gorbunov and Hamilton, 1997 using their eq. (21) except
	!! we do not assume that the density of the solution is the density
	!! of water (we use the calculated density)
				SatVapPressRatio = DEXP(										&
				InParticle%SurfaceTension / Rstar / InParticle%Temperature *	&
				2.*Pi * InParticle%InsolubleRadius *							&
				WaterMolecMass / InParticle%SolutionDensity * dRdV *			&
				((2./XX*(1.+1./YY) -											&
				CosContactAngle/YY + (XX-CosContactAngle) &
				/YY/YY/YY*(1.-XX*CosContactAngle))  & !d(As/sd)/dre
				-CosContactAngle * XX*XX/YY * (1.-((XX-CosContactAngle)/YY)**2.))) !d(Aa/s)/dre

			END IF
		ELSE

!! The Kelvin Effect is abstracted to a function 
!! housed in "ParticleAttributes.h"
			!WRITE(*,*) "Before CurveCorr"
			SatVapPressRatio  = CurvatureCorrection(InParticle)
			!WRITE(*,*) "After CurveCorr"
		END IF	
	END IF

	RETURN
END FUNCTION SatVapPressRatio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! --- EQUILIBRATE ONE AQUEOUS DISSOCIATION REACTION ---       !!
!!							       !!
!! Iterate a particular reaction until it reaches equilibrium. !!
!!							       !!
!! ReturnType: 1 is normal, 2 is "already was in eq"	       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateDissolutionReaction (GasPhaseBurden, WhichRxn, & 
               InParticle, Temperature, ReturnType, &
	       WaterSatBurden, GasPhaseBurden2, OptErrorTolerance)

	USE Chemistry,	ONLY :	AqEquilibriaList,		&
				HowManyAqChems,			&
				HowManyAqCations,		&
				AqPhaseChemicalNames,	&
				HowManyEvolveGasChems,	&
				GasPhaseChemicalNames

	USE ModelParameters,	ONLY :	AqThermoEquilibriumError,	&
					WaterContentPrecision,		&
					AqThermoNumbAllowedIterations,	&
					AerosolWaterEquilibriumRHThreshold,&
					Avogadro,			&
					moles, grams, WaterMolecMass, micron 

	USE GridPointFields,	ONLY :  GetRelativeHumidity

	USE InfrastructuralCode,ONLY :	INT2STR, REAL2STR,		&
					Transcript,			&
					ERROR, WARN

	USE Aerosols,		ONLY :	DumpParticleContentsAtError,     &
					EqCoeff,			 &
					MeanActivity,                    &
					EquilibrateAqueousDissociationReaction

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReturnType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature
	REAL*8 :: OptErrorTolerance, WaterSatBurden, GasPhaseBurden2

	!! Local Variables
	INTEGER :: I, GG
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, ActivityCoefficients, &
                   ErrorTolerance, DissociatedScaling, DissociatedExp
	LOGICAL :: WaterRxn

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .TRUE. ! .TRUE.

	IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

        !WRITE(*,*) "Entering EquilibrateDissolutionReaction, Reaction # ", WhichRxn

	ErrorTolerance = OptErrorTolerance

        !!This checks for solid-gas reactions that are already in
        !!equilibrium, as the solid conc = 0.0 and the gas-phase is subsaturated
	IF (INT(DissolutionData(WhichRxn,11)) .EQ. 6.)  THEN
	   EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, & 
                        WhichRxn, Temperature, InParticle, &
                        WaterSatBurden=WaterSatBurden, &
                        GasPhaseBurden2=GasPhaseBurden2)

		IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		   .AND. InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) &
                    .LT. 1e-27) THEN !There is no solid left
			ReturnType = 2
			RETURN
		END IF

	END IF	

        !!This checks for Type 2 reactions that are subsaturated in the gas-phase, 
        !!but there is not enough of an ion to evaporate more.
        IF (INT(DissolutionData(WhichRxn,11)) .EQ. 2.)  THEN
               EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)

		IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		   .AND. (InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) &
                    .LT. 1.0e-28 .OR. InParticle%AqChems(INT(DissolutionData(WhichRxn,4))) &
                    .LT. 1.0e-28)) THEN !There is not enough ions left
                        !WRITE(*,*) "Too Small Caught!"
			ReturnType = 2
			RETURN
		END IF
        END IF 

	!This checks to make sure the NH3(g) <=> NH3(aq)
	!doesn't crash because there isn't enough NH3(aq)
	IF (INT(DissolutionData(WhichRxn,11)) .EQ. 1. .AND. &
      GasPhaseChemicalNames(INT(DissolutionData(WhichRxn,1))) .EQ. "NH3")  THEN
		!WRITE(*,*) "Check Started for Reaction: ", WhichRxn
		EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)

		!WRITE(*,*) "Ratio Okay for Reaction: ", WhichRxn
		IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
			.AND. InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .LT. 1e-26) THEN !There is no NH3 left
			ReturnType = 2
			RETURN
		END IF
		
	END IF	

        !WRITE(*,*) "After initial checks"
	
	!! Presume innocence in the form of a successful return.
	ReturnType = 1
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine generally follows the Mass Flux Iteration method reviewed!! 
!! in Jacobson, 1999 (the textbook) and reviewed earlier in              !!
!! Jacobson et al. 1996, and Villars 1959                                !!
!!                                                                       !!
!! Jacobson, M.Z., A. Tabazadeh, and R.P. Turco, Simulating equilibrium  !!
!!       within aerosols and nonequilibrium between gases and aerosols,  !!
!!       Journal of Geophysical Research, 101 (D4), 9079-9091, 1996.     !!
!!									 !!
!! Villars, D.S., A method of successive approximations for computing    !!
!!    combustion equilibria on a high speed digital computer,            !!
!!    Journal of Physical Chemistry, 63, 521-5, 1959.			 !!
!!									 !!
!! STEP 1: Find the most aberrant ratio of concentration to              !!
!!         stoicheometric coefficient  (cf., Jacobson 1999, p.497)	 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Check for already equilibrated (either by none of 
!! the species being present or explicit equilibrium)
!!
!! And then pre-calculate the Qnumer and Qdenomer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There are several reactions types:				             !!
!!									     !!
!! 0. Water Condensation.						     !!
!! 1. Direct dissolution of a single species not incorporating dissociation. !!
!! 2. Dissolution of a single species that then dissociates directly.	     !!
!! 4. Dissolution of a species that binds with an ion on the way in.	     !!
!! 5. Dissolution of two gas phase species to form 1 or 2 aqueous species.   !!
!!									     !!
!! Identify which we're dealing with and go from there.(Type 3 doesn't apply.)!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SELECT CASE (INT(DissolutionData(WhichRxn,11)))

	CASE(0)

		EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)

		IF  (ABS(EqRatio-1.) .LT. AqThermoEquilibriumError) THEN
			ReturnType = 2
			RETURN
		END IF

                !! If we want to include water equilibration in this, 
                !! then equilibrate it and return
		IF (GetRelativeHumidity() .GT. &
		    AerosolWaterEquilibriumRHThreshold) &
		    CALL ERROR("Passed EquilibrateDissolutionReaction() in CondensationRelatedFunctions.h the water equilibration " //  &
				   "reaction even though the RH is above the threshold you specify for allowable equilibration ("//			&
				   TRIM(REAL2STR(AerosolWaterEquilibriumRHThreshold))//").  You must filter for that outside of this routine.")

		RETURN

	CASE(1)

		!WRITE(*,*) "Reaction Check 1", GasPhaseBurden, WhichRxn, Temperature
		EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)
		!WRITE(*,*) "Reaction Check 2"
		IF  (ABS(EqRatio-1.) .LT. AqThermoEquilibriumError) THEN
			ReturnType = 2
			RETURN
		END IF

		Qnumer = InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))
		Qdenom = GasPhaseBurden / Avogadro / InParticle%NumberOfParticles
		!WRITE(*,*) "Reaction Check 3"
		
	CASE(2)
		!WRITE(*,*) "Reaction Check 1"
		EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)
		!WRITE(*,*) "Reaction Check 2"

		!! If all of the concentrations are zero, then equilibrium has effectively been established.
		IF  ((ABS(EqRatio-1.) .LT. AqThermoEquilibriumError)) THEN
			ReturnType = 2
			RETURN
		END IF
		Qnumer = MIN( InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) / DissolutionData(WhichRxn,13),	&
					  InParticle%AqChems(INT(DissolutionData(WhichRxn,4))) / DissolutionData(WhichRxn,14) )

		Qdenom = GasPhaseBurden / Avogadro 	/ InParticle%NumberOfParticles
		!WRITE(*,*) "Reaction Check 3"

	CASE(4)
		!WRITE(*,*) "Reaction Check 1"
		EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)
		!WRITE(*,*) "Reaction Check 2"

		!! If all of the concentrations are zero, then equilibrium has effectively been established.
		IF  ((ABS(EqRatio-1) .LT. AqThermoEquilibriumError)) THEN
			ReturnType = 2
			RETURN
		END IF

		Qnumer = InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))
		Qdenom = MIN(GasPhaseBurden / Avogadro  / InParticle%NumberOfParticles,	&
					 InParticle%AqChems(INT(DissolutionData(WhichRxn,2))))
		IF(Qnumer .LT. 0. .OR. Qdenom .LT. 0.) CALL ERROR("Conc below 0 in EqDissRxn!")
		!WRITE(*,*) "Reaction Check 3"
	
	CASE(6)

		!WRITE(*,*) "Before Eq Ratio"
		EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)
		!WRITE(*,*) "Initial EqRatio: ", EqRatio
		
		IF  ((ABS(EqRatio-1) .LT. AqThermoEquilibriumError)) THEN
			ReturnType = 2
			RETURN
		END IF

		Qnumer = InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) 
		Qdenom = MIN(GasPhaseBurden,GasPhaseBurden2)
		Qdenom = Qdenom / Avogadro / InParticle%NumberOfParticles
		!WRITE(*,*) "Type 5", Qnumer, Qdenom
	END SELECT

	!WRITE(*,*) "Check 4"
	!! I counts the number of iterations
	I = 1

	!! STEP 2: Calculate the Mass Step Sizes by which to Correct Concentrations
	Z  = (Qdenom + Qnumer) / 2.
	dX = Qdenom - Z				!! The Mass Flux Factor
	!WRITE(*,*) "dx: ", dX, " Z: ", Z, "Qdenom: ", Qdenom, "Qnumer: ", Qnumer

	!! Loop over Corrections until convergence if there is anything to do
	IF (.NOT. (Qnumer .EQ. 0 .AND. Qdenom .EQ. 0)) THEN
	DO WHILE ((ABS(EqRatio-1) .GT. ErrorTolerance))

	!WRITE(*,*) "Check 5"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 3: Adjust Molalities: !!
	!!							  !!
	!! -- Ions First			  !!
	SELECT CASE (INT(DissolutionData(WhichRxn,11)))

	CASE (1)

		InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) = InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) + dX
		IF (DissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems) &
			GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles

	CASE (2)

		InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))	= InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))	&
											+ dX * DissolutionData(WhichRxn,13)

		InParticle%AqChems(INT(DissolutionData(WhichRxn,4)))	= InParticle%AqChems(INT(DissolutionData(WhichRxn,4)))	&
											+ dX * DissolutionData(WhichRxn,14)

		IF (DissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems)	&
			GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles
	
	CASE (4)

		InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))	= InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) + dX
		InParticle%AqChems(INT(DissolutionData(WhichRxn,2)))	= InParticle%AqChems(INT(DissolutionData(WhichRxn,2))) - dX
		IF (DissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems)		&
			GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles

	CASE (6)

		InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))	= InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) + dX
		IF (DissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems)		&
		  GasPhaseBurden  = GasPhaseBurden  - dX * Avogadro * InParticle%NumberOfParticles
		IF (DissolutionData(WhichRxn,2) .LE. HowManyEvolveGasChems)		&
		  GasPhaseBurden2 = GasPhaseBurden2 - dX * Avogadro * InParticle%NumberOfParticles

		!WRITE(*,*) "Molalities Okay"
	END SELECT


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) "Check 7"
	Z = 0.5*Z

	!! The recalculation is based on the equilibrium ratio
	IF (INT(DissolutionData(WhichRxn,11)) .EQ. 6) THEN
		EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)
		!WRITE(*,*) "Ratio Okay"
	ELSE
		!WRITE(*,*) "Check 8", GasPhaseBurden, WhichRxn
		EqRatio = DissolutionEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle, &
				  WaterSatBurden=WaterSatBurden, GasPhaseBurden2=GasPhaseBurden2)
		!WRITE(*,*) "Check 9"
	END IF

	!WRITE(*,*) "Check 10"
	IF (EqRatio .GE. 1.+ErrorTolerance) dX = -1.*Z
	IF (EqRatio .LE. 1.-ErrorTolerance) dX = Z
	!WRITE(*,*) EqRatio, Z, dX

	!! Count the number of attempts before convergence
	I = I+1

	IF (INT(DissolutionData(WhichRxn,11)) .EQ. 6 .AND. InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .LT. 1.0e-27) THEN 
                !Much Less than a molecule per particle
		!WRITE(*,*) "Solid Electrolyte ", trim(AqPhaseChemicalNames(DissolutionData(WhichRxn,3))), " is negligible.", I		
		RETURN
	END IF

	!This checks to make sure the NH3(g) <=> NH3(aq)
	!doesn't crash because there isn't enough NH3(aq)
	IF (DissolutionData(WhichRxn,11) .EQ. 1. .AND. &
		GasPhaseChemicalNames(INT(DissolutionData(WhichRxn,1))) .EQ. "NH3")  THEN		
		!WRITE(*,*) "Ratio Okay for Reaction: ", WhichRxn
		IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
			.AND. InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .LT. 1e-26) THEN !There is no NH3 left
			RETURN
		END IF	
	END IF	

        !Check to make sure that there is some of the ions present for Type 2 rxns
        IF (DissolutionData(WhichRxn,11) .EQ. 2.)  THEN
		IF ((EqRatio .GT. 1.0) & !The gas-phase is sub-saturated
		   .AND. (InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))+ dX * DissolutionData(WhichRxn,13) &
                    .LT. 1.0e-28 .OR. InParticle%AqChems(INT(DissolutionData(WhichRxn,4)))+ dX * DissolutionData(WhichRxn,13) &
                    .LT. 1.0e-28)) THEN !dX will force AqChems very small or negative on next step
                        !WRITE(*,*) "Not enough!", EqRatio, InParticle%AqChems(INT(DissolutionData(WhichRxn,3))), InParticle%AqChems(INT(DissolutionData(WhichRxn,4))), dX

                        !If not negative on next step, do one step and return
                        IF(InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))+ dX * DissolutionData(WhichRxn,13) &
                         .GT. 0.0 .AND. InParticle%AqChems(INT(DissolutionData(WhichRxn,4)))+ dX * DissolutionData(WhichRxn,13) &
                         .GT. 0.0) THEN
                            !Write(*,*) "Postive branch"
			    InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) = InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))	&
											+ dX * DissolutionData(WhichRxn,13)

			    InParticle%AqChems(INT(DissolutionData(WhichRxn,4))) = InParticle%AqChems(INT(DissolutionData(WhichRxn,4)))	&
											+ dX * DissolutionData(WhichRxn,14)

		            IF (DissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems)	&
			         GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles
                            ReturnType = 2
                            RETURN

                        !!If next step is negative, set dX so that minimum is 1.0e-28 and return
                        ELSE
                            !Write(*,*) "Negative branch"
                            IF (InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))+ dX * DissolutionData(WhichRxn,13) &
                                .LE. 0.0) THEN
                                 dX = (1.0e-28-InParticle%AqChems(INT(DissolutionData(WhichRxn,3))))/ DissolutionData(WhichRxn,13)
                            ELSE
                                 dX = (1.0e-28-InParticle%AqChems(INT(DissolutionData(WhichRxn,4))))/ DissolutionData(WhichRxn,14)
                            ENDIF

			    InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) = InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))	&
											+ dX * DissolutionData(WhichRxn,13)

			    InParticle%AqChems(INT(DissolutionData(WhichRxn,4))) = InParticle%AqChems(INT(DissolutionData(WhichRxn,4)))	&
											+ dX * DissolutionData(WhichRxn,14)

		            IF (DissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems)	&
			         GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles
                            ReturnType = 2
                            RETURN
                       ENDIF
		END IF
        END IF

	!Just exit if no convergence after 2000 iterations(MJA, 10/05/07)
        ! Change to 200 (MJA, 01/30/2013)
	IF (I .GT. 100) THEN
		!WRITE(*,*) "Over 100 iterations in EquilibrateDissolutionReaction for Reaction # ", WhichRxn
		ReturnType = 2
		EXIT
		
	END IF

	!WRITE(*,*) "Check 11"
	END DO ; END IF

	RETURN
END SUBROUTINE EquilibrateDissolutionReaction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Check the Ratio between the Existing and Ideal Equilibrium   !!
!! Constants.  This is used by several other routines.		!!
!!								!!
!! When the equations are read, the program tries to figure out !!
!! whether an uncharged electrolyte is dissociating (in which   !!
!! case the denominator activity coefficient is simply 1), and  !!
!! if not it stores an equivalent activity coefficient ratio.   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This should only be used in the context of EquilibrateDissolutionReaction!  
!! Otherwise there is no check for lack of the species.
REAL*8 FUNCTION DissolutionEquilibriumConstantsRatio (GasPhaseBurden, &
                       WhichRxn, Temperature, InParticle,	&
			GasPhaseBurden2, WaterSatBurden)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyAqChems, &
                                    HowManyAqCations, AqPhaseChemicalNames
	USE Aerosols,        ONLY : MeanActivity

	USE ModelParameters, ONLY : RStarMB, moles, grams, & 
                                    WaterMolecMass, Avogadro

	USE InfrastructuralCode, ONLY : Error

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReactionType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature
	REAL*8 :: GasPhaseBurden2, WaterSatBurden 

	!! Internal Variables
	REAL*8	:: Eq, ExpA, ExpB, P1, P2, MeanAct ! AwPow,
        INTEGER :: IonA, IonB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There are several reactions types:				             !!
!!									     !!
!! 0. Water Condensation.						     !!
!! 1. Direct dissolution of a single species not inclorporating dissociation.!!
!! 2. Dissolution of a single species that then dissociates directly.	     !!
!! 4. Dissolution of a species that binds with an ion on the way in.	     !!
!! 6. Two gas phase species forming a solid.				     !!
!! Identify which we're dealing with and go from there.(Type 3 doesn't apply.)!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! See if the reaction is empty !!
	SELECT CASE (INT(DissolutionData(WhichRxn,11)))

	CASE (1)
		!WRITE(*,*) "Flag 1"
		!! If there is none of the chemical to deal with
		IF (InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .EQ. 0. .AND.	&
			GasPhaseBurden .EQ. 0) THEN
			DissolutionEquilibriumConstantsRatio = 1.
			RETURN
		END IF
		!WRITE(*,*) "Flag 2"


	CASE (2)
		!WRITE(*,*) "Flag 1"
		!! If there is none of the chemical to deal with
		IF ((InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .EQ. 0.  .OR.	  &
			 InParticle%AqChems(INT(DissolutionData(WhichRxn,4))) .EQ. 0.) .AND.  &
			 GasPhaseBurden .EQ. 0.) THEN
			DissolutionEquilibriumConstantsRatio = 1.
			RETURN
		END IF
		!WRITE(*,*) "Flag 2"

	CASE (4)

		!! If there is none of the chemical to deal with
		!WRITE(*,*) "Flag 1", WhichRxn, AqPhaseChemicalNames(DissolutionData(WhichRxn,2)), GasPhaseBurden !InParticle%AqChems(DissolutionData(WhichRxn,2))
		IF ( (InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .EQ. 0.  .OR.	  &
			 InParticle%AqChems(INT(DissolutionData(WhichRxn,2))) .EQ. 0.) .AND.   &
			 GasPhaseBurden .EQ. 0.) THEN
			DissolutionEquilibriumConstantsRatio = 1.
			!WRITE(*,*) "Flag 2"
			RETURN
		END IF

	CASE (6)

		!! If there is none of the chemical to deal with
		IF ((InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .EQ. 0.) .AND.   &
			 (GasPhaseBurden .EQ. 0. .OR. GasPhaseBurden2 .EQ. 0)) THEN
			DissolutionEquilibriumConstantsRatio = 1.
			RETURN
		END IF

	END SELECT


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Prepare the needed quantities !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Flag 3"
	!! Get the molality of the electrolyte 
        !! and the equilibrium at this temperature
	Eq = HenrysLaw (Temperature, WhichRxn)

	!! If the solution is set to blow up, then establish a large artificial value.
	IF (GasPhaseBurden .EQ. 0. .OR. (Eq .EQ. 0. .AND. WhichRxn .NE. 1)) THEN
		DissolutionEquilibriumConstantsRatio = 1.e12
		RETURN
	END IF
	!WRITE(*,*) "Flag 4"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make the primary calculation !!
	SELECT CASE (INT(DissolutionData(WhichRxn,11)))

	CASE (0)

		DissolutionEquilibriumConstantsRatio = WaterSatBurden * SatVapPressRatio (InParticle, 1) / GasPhaseBurden

	CASE (1)		
		!WRITE(*,*) "Before Ratio", INT(DissolutionData(WhichRxn,3))

		DissolutionEquilibriumConstantsRatio =			&
	SatVapPressRatio (InParticle, INT(DissolutionData(WhichRxn,3))) &
	* 1000. * InParticle%AqChems(INT(DissolutionData(WhichRxn,3)))	&
	/ Eq / RstarMB / Temperature					&
	/ InParticle%AqChems(1) / moles * grams / WaterMolecMass	&
		* Avogadro / GasPhaseBurden 
				
		!WRITE(*,*) "Ratio", DissolutionEquilibriumConstantsRatio
	CASE (2)
		!WRITE(*,*) "Flag 5"
		IonA = INT(DissolutionData(WhichRxn,3))
		IonB = INT(DissolutionData(WhichRxn,4))
		ExpA = INT(DissolutionData(WhichRxn,13))
		ExpB = INT(DissolutionData(WhichRxn,14))
		!WRITE(*,*) "Flag 6", INT(DissolutionData(WhichRxn,3)) 
		
		DissolutionEquilibriumConstantsRatio =										&
					SatVapPressRatio (InParticle , INT(DissolutionData(WhichRxn,3))) 			&
					* InParticle%AqChems(IonA)**ExpA * InParticle%AqChems(IonB)**ExpB			&
					* (1000. * grams / moles / WaterMolecMass / InParticle%AqChems(1))**(ExpA+ExpB) &
					/ Eq / RstarMB / Temperature								& 
					*InParticle%GammaMixed(INT(DissolutionData(WhichRxn,12)))**				& ! 3 lines of
					(AqEquilibriaList(INT(DissolutionData(WhichRxn,12)),4)+					& ! mean act
					 AqEquilibriaList(INT(DissolutionData(WhichRxn,12)),5))					&
					 * Avogadro / GasPhaseBurden
			
	CASE (4)

		IonA = FLOOR(DissolutionData(WhichRxn,8))
		IonB = FLOOR(DissolutionData(WhichRxn,9))
		ExpA = ANINT(10000.*(DissolutionData(WhichRxn,8)-FLOOR(DissolutionData(WhichRxn,8))))
		ExpB = ANINT(10000.*(DissolutionData(WhichRxn,9)-FLOOR(DissolutionData(WhichRxn,9))))

		IF (InParticle%AqChems(INT(DissolutionData(WhichRxn,2))) .LE. 0) THEN
			DissolutionEquilibriumConstantsRatio = 1.e12
		ELSE  IF (InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .LE. 0) THEN 
			DissolutionEquilibriumConstantsRatio = 0.0000000001
		ELSE
			DissolutionEquilibriumConstantsRatio =										&
					SatVapPressRatio (InParticle, INT(DissolutionData(WhichRxn,3))) *	&
					InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) /					&
					InParticle%AqChems(INT(DissolutionData(WhichRxn,2)))						&
					*InParticle%GammaMixed(IonA)**ExpA									& 
					/InParticle%GammaMixed(IonB)**ExpB									&
					/ Eq / RstarMB / Temperature										&
					* Avogadro / GasPhaseBurden
			IF(DissolutionEquilibriumConstantsRatio .GT. 1.0e12) THEN
				DissolutionEquilibriumConstantsRatio = 1.e12
			ELSE IF (DissolutionEquilibriumConstantsRatio .LT. 1.0e-10) THEN
				DissolutionEquilibriumConstantsRatio = 0.0000000001
			END IF
		END IF

	CASE (6)

		!WRITE(*,*) "HenrysLaw", Eq
		IF (InParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .EQ. 0.) THEN
			DissolutionEquilibriumConstantsRatio = 1.e12
		ELSE

		!WRITE(*,*) "Eq Const", Eq
		!STOP
		
		!Partial Pressures of gas-phase species in mbar (Burden = molecules/cm3)
		P1 = GasPhaseBurden*RstarMB*Temperature/Avogadro &
		     * SatVapPressRatio (InParticle, INT(DissolutionData(WhichRxn,3)))
		
		P2 = GasPhaseBurden2*RstarMB*Temperature/Avogadro &
		     * SatVapPressRatio (InParticle, INT(DissolutionData(WhichRxn,3)))

		!Note that we ignore curvature correction here!
		DissolutionEquilibriumConstantsRatio =	 1.0 / Eq / P1 / P2

		END IF

	END SELECT
!WRITE(*,*) "Flag 8"
	RETURN
END FUNCTION DissolutionEquilibriumConstantsRatio

SUBROUTINE EquilibrateSulfate()

!! This subroutine pushes all of the H2SO4 in the gas-phase into
!! the aerosol phase. The fraction of gas-phase sulfate added to each
!! size bin is the ratio of the bin mass transfer coefficient
!! to the sum of the coefficients of all bins.

!! NOTE that this does NOT model the kinetic transfer of sulfate!
!! All gas-phase H2SO4 is forced to condense instantaneously by this
!! routine. This routine should NOT be used to do kinetic 
!! mass-transfer. 

	!! Include modules
	USE InfrastructuralCode, ONLY : SetFileHandleCounter, &
					Transcript, &
					ERROR, &
					REAL2STR, &
					INT2STR, &
					GetFileHandle

	USE Chemistry,		 ONLY : HowManyGasChems, &
				        HowManyEvolveGasChems, &
					GasPhaseChemicalNames, &
					HowManyAqChems, &
					HowManyAqCations, &
					HowManyAqAnions, &
					AqPhaseChemicalNames, &	
					AqCationNames, &
					AqAnionNames, &
					FindChem, &
					MolecularDiffusionCoefficient

	USE GridPointFields,	 ONLY : GridGasChem, &
                                        UpdateChemicalConcentrations, &
					GetM, &
					GetTemp, &
					MeanFreePathofAir, &
					GetDynamicViscosityofAir, &
					GetAirDensity

	USE Aerosols,		 ONLY : Particle, &
					Particles,  &
				FindAqElectrolyteEquilibriumForGridPoint, &
					FindElectrolyteEquilibrium, &
					KusikMeissner, &
					UpdateWaterActivity, &
					AerosolPH, &
					SortAerosolAtGridPointForCoagulation, &
					ParticleDensity, &
					AerosolModes, &
					ReynoldsNumber
	USE	ModelParameters
		
	IMPLICIT NONE

	!! Internal Variables
	integer :: i, j, k, l, seedarray(2),q, r, OutFH(4), NumBins, A, S, &
                   HeteroBins, H2SO4GasIndex, SO4AqIndex,         &
                   HNO3GasIndex, NumAqChems, &
	           NH3GasIndex, HAqIndex
	real*8  :: ii, jj, kk, Dum1, Dum2
	real*8  ::  Mgas, Temp, H2SO4Conc, TransferSum, &
	            FractionofSulfate, Knudsen, InvKn, AccomCoeff, &
				omega, Re, GasMeanFreePath, DynamVisc, Dens, &
				Dv, Sc, DumNum, VentFac, Total

	REAL*8, ALLOCATABLE	::	Chem(:), NumConc(:), SulfTrans(:)

	TYPE(Particle), POINTER :: CurrentParticle


		!Count # of bins, # of Aq. Chems
		CurrentParticle => particles%first
		NumBins = 0
 8		if (associated(currentparticle)) then
			NumBins = NumBins + 1
			NumAqChems = size(currentparticle%AqChems)
			if (associated(currentparticle%next)) then
				currentparticle => currentparticle%next
				goto 8
			end if
		end if
		ALLOCATE(NumConc(NumBins))
		ALLOCATE(SulfTrans(NumBins))
		
		!Find H2SO4 in gas and aqueous phases
		H2SO4GasIndex = FindChem("H2SO4", 0)
		SO4AqIndex = FindChem("SO4--", 1)
		HAqIndex = FindChem("H+", 1)

		!Calculate gas conc of H2SO4 in mol/cm3
		H2SO4Conc = GridGasChem(H2SO4GasIndex)/ Avogadro !mol/cm3
			
		!Calculate sum of the mass transfer coefficients for all particles
		CurrentParticle => particles%first
		TransferSum = 0				
		I = 1
 33		if (associated(currentparticle)) then
    			Write(*,*) "Number of Particles", CurrentParticle%NumberofParticles					
		   IF(CurrentParticle%NumberofParticles .GT. 1e-9) THEN
                                        
		        !Calculate correction for collision geometry 
		        !and sticking probability (omega)
		        GasMeanFreePath = MeanFreePathOfAir()
		        Knudsen = GasMeanFreePath / CurrentParticle%EffectiveRadius 
		        InvKn = 1.0/Knudsen
	                AccomCoeff = 0.65 !From Poschl et al.
	                omega = (1 + ((1.33 + 0.71*InvKn)/(1.0 + InvKn) + 4*(1-AccomCoeff)/3.0/AccomCoeff)*Knudsen)**(-1)
                                       
		        !Calculate the Ventilation Factor
			Re = ReynoldsNumber (CurrentParticle)
			DynamVisc = GetDynamicViscosityOfAir()
			Dens = GetAirDensity()
			Dv = MolecularDiffusionCoefficient (H2SO4GasIndex)						
			Sc = DynamVisc/Dens/Dv
			DumNum = (Re**(0.5))*(Sc**(0.33333333333))
			IF (DumNum .LE. 1.4) THEN
				VentFac = 1 + 0.108*DumNum*DumNum
			ELSE
				VentFac = 0.78 + 0.308*DumNum
			END IF

			!Get number concentration of particles
			NumConc(I) = currentparticle%numberofparticles
						
			!Calculate mass transfer coefficient
			SulfTrans(I) = 4*Pi*Dv*omega*VentFac*currentparticle%numberofparticles*currentparticle%effectiveradius
			TransferSum = TransferSum + SulfTrans(I)

		!Ignore bins with no particles
		ELSE
			SulfTrans(I) = 0.0
			TransferSum = TransferSum + SulfTrans(I)
		END IF
				
		if (associated(currentparticle%next)) then
			currentparticle => currentparticle%next
			I = I+ 1
			goto 33
		end if
	end if 
				
		
	!Check to make sure gas-phase fractions sum to 1
	Total = 0.0
	DO I = 1, NumBins
		FractionofSulfate = SulfTrans(I)/TransferSum
		Total = Total + FractionofSulfate
	END DO
	
	IF(Total .LT. 0.999999999 .OR. Total .GT. 1.000000001) THEN
                WRITE(*,*) "Fraction Sum: ", Total
		CALL ERROR("Error in EquilibrateSulfate(). Fractions do not add to 1.")
	END IF

	!Allocate gas-phase sulfate to bins
	!The fraction of gas-phase sulfate added to each size bin
	!Is the ratio of the bin mass transfer coefficient
	!to the sum of the coefficients of all bins.
	!Here, all the sulfate is added as SO4-- so that
	!it is dissociated when EquilibrateGridPoint is called
	 
	CurrentParticle => particles%first
	I = 1
20      if (associated(currentparticle)) then
    								
		FractionofSulfate = SulfTrans(I)/TransferSum

		!Add sulfate
		currentparticle%AqChems(SO4AqIndex) = currentparticle%AqChems(SO4AqIndex) + H2SO4Conc*FractionofSulfate/NumConc(I)

		!AddProtons
		currentparticle%AqChems(HAqIndex) = currentparticle%AqChems(HAqIndex) + 2*H2SO4Conc*FractionofSulfate/NumConc(I)
			
		I = I + 1

		if (associated(currentparticle%next)) then
			currentparticle => currentparticle%next
			goto 20
		end if
	end if 
		
	!Update gas H2SO4 concentration by setting it equal to 0.
	ALLOCATE(CHEM(HowManyGasChems))
	DO q = 1, HowManyGasChems
		Chem(q) = GridGasChem(q) !molecules/cm3
	END DO
	
	Chem(H2SO4GasIndex) = 0.0
	CALL UpdateChemicalConcentrations(Chem(1:HowManyEvolveGasChems))
	
	RETURN
END SUBROUTINE EquilibrateSulfate

SUBROUTINE EquilibrateGasSolidReactions(GasPhaseBurden, WhichRxn, & 
                   Temperature, ReturnType, &
		   GasPhaseBurden2, OptErrorTolerance)

!! This subroutine calculates the bulk equilibrium of NH4Cl and NH4NO3 
!! and then weights the mass transfer of those species by the surface area
!! of the particle bins.

	!! Include modules
	USE InfrastructuralCode, ONLY : SetFileHandleCounter, &
					Transcript, &
					ERROR, &
					REAL2STR, &
					INT2STR, &
					GetFileHandle

	USE Chemistry,		 ONLY : HowManyGasChems, &
					HowManyEvolveGasChems, &
					GasPhaseChemicalNames, &
					HowManyAqChems, &
					HowManyAqCations, &
					HowManyAqAnions, &
					AqPhaseChemicalNames, &	
					AqCationNames, &
					AqAnionNames, &
					FindChem, &
					MolecularDiffusionCoefficient

	USE GridPointFields,	 ONLY : GridGasChem, &
					UpdateChemicalConcentrations, &
					GetM, &
					GetTemp, &
					MeanFreePathofAir, &
					GetDynamicViscosityofAir, &
					GetAirDensity

	USE Aerosols,		 ONLY : Particle, &
				        Particles,  &
				FindAqElectrolyteEquilibriumForGridPoint, &
					FindElectrolyteEquilibrium, &
					KusikMeissner, &
					UpdateWaterActivity, &
					AerosolPH, &
					SortAerosolAtGridPointForCoagulation, &
					ParticleDensity, &
					AerosolModes, &
					ReynoldsNumber
	USE	ModelParameters
		
	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReturnType
	REAL*8 :: GasPhaseBurden, Temperature
	REAL*8 :: OptErrorTolerance, GasPhaseBurden2

	!! Local Variables
	INTEGER :: I
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, ActivityCoefficients, &
                   ErrorTolerance, DissociatedScaling, DissociatedExp
	LOGICAL :: WaterRxn

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .TRUE. ! .TRUE.

	integer :: j, k, l, seedarray(2),q, r, OutFH(4), NumBins,  &
	           Gas, Aq
	real*8  :: jj, kk, P1, P2, Eq, Total
	real*8  ::  Mgas, Temp, AreaSum, SolidSum, NewSolidSum, &
	            DeltaSolid, FractionofMass, ZeroCheck, &
				RemainingSolidtoEvaporate

	REAL*8, ALLOCATABLE	::	Chem(:), NumConc(:), SurfArea(:)

	TYPE(Particle), POINTER :: CurrentParticle


		!Count # of bins
		CurrentParticle => particles%first
		NumBins = 0
		8 if (associated(currentparticle)) then
			NumBins = NumBins + 1
			if (associated(currentparticle%next)) then
				currentparticle => currentparticle%next
				goto 8
			end if
		end if
		ALLOCATE(NumConc(NumBins))
		ALLOCATE(SurfArea(NumBins))
		
		!Calculate total surface area and surface area of each bin		
		CurrentParticle => particles%first
		AreaSum = 0				
		I = 1
		30 if (associated(currentparticle)) then
    								
				IF(CurrentParticle%NumberofParticles .GT. 0.) THEN
					
					SurfArea(I) = 4*Pi*currentparticle%numberofparticles*(currentparticle%effectiveradius**2.0)
					AreaSum = AreaSum + SurfArea(I)

				!Ignore bins with no particles
				ELSE
					SurfArea(I) = 0.0
					AreaSum = AreaSum + SurfArea(I)
				END IF
				
				if (associated(currentparticle%next)) then
					currentparticle => currentparticle%next
					I = I+ 1
					goto 30
				end if
			end if 

	!! Calculate the total ammount of the solid species initially
	!! present in the aerosol.
		CurrentParticle => particles%first
		SolidSum = 0				
		40 if (associated(currentparticle)) then
    								
				!Units of mol/cm3
				IF (CurrentParticle%NumberofParticles .GT. 0.) THEN
					SolidSum = SolidSum + CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3)))*CurrentParticle%NumberofParticles
					!WRITE(*,*) "SolidSum", WhichRxn, SolidSum
				END IF
				if (associated(currentparticle%next)) then
					currentparticle => currentparticle%next
					I = I+ 1
					goto 40
				end if
			end if 
		NewSolidSum = SolidSum
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!           CALCULATE BULK EQUILIBRIUM               !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	!! Set the error tolerance according to the 
	!! subroutine's optinoal input.
	!! THIS IS NOT THE TOLERANCE FOR THE EQUILIBRATION
	!! CHECK (WHICH USES THE GLOBAL VALUE), RATHER IT 
	!! IS THE VALUE WE EQUILIBRATE TO IF THAT CHECK FAILS	
	ErrorTolerance = OptErrorTolerance

	IF (DissolutionData(WhichRxn,11) .NE. 6.)			   &
		CALL ERROR("You can only call EquilibrateGasSolidReactions() for a reaction that involves two gas phase species (reaction "// &
				   "type = 6)")
		
!!Matt Alvarado added this to check for solid-gas reactions that are already in
!!Equilibrium, as the solid conc = 0.0 and the gas-phase is subsaturated	
	IF (DissolutionData(WhichRxn,11) .EQ. 6.)  THEN
		!WRITE(*,*) "Check Started for Reaction: ", WhichRxn
		Eq = HenrysLaw (Temperature, WhichRxn)
		P1 = GasPhaseBurden*RstarMB*Temperature/Avogadro
		P2 = GasPhaseBurden2*RstarMB*Temperature/Avogadro

		EqRatio = 1.0 / Eq / P1 / P2
		!WRITE(*,*) "Ratio Okay for Reaction: ", WhichRxn
		IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
			.AND. NewSolidSum .LT. 1e-27) THEN !There is no solid left
			ReturnType = 2
			RETURN
		END IF
		

	END IF	
	
	!! Presume innocence in the form of a successful return.
	ReturnType = 1
	!WRITE(*,*) "Start EquilibrateDissolutionReaction for Reaction: ", WhichRxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine generally follows the Mass Flux Iteration method reviewed!! 
!! in Jacobson, 1999 (the textbook) and reviewed earlier in              !!
!! Jacobson et al. 1996, and Villars 1959                                !!
!!                                                                       !!
!! Jacobson, M.Z., A. Tabazadeh, and R.P. Turco, Simulating equilibrium  !!
!!       within aerosols and nonequilibrium between gases and aerosols,  !!
!!       Journal of Geophysical Research, 101 (D4), 9079-9091, 1996.     !!
!!									 !!
!! Villars, D.S., A method of successive approximations for computing    !!
!!    combustion equilibria on a high speed digital computer,            !!
!!    Journal of Physical Chemistry, 63, 521-5, 1959.			 !!
!!									 !!
!! STEP 1: Find the most aberrant ratio of concentration to              !!
!!         stoicheometric coefficient  (cf., Jacobson 1999, p.497)	 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Check for already equilibrated (either by none of 
!! the species being present or explicit equilibrium)
!!
!! And then pre-calculate the Qnumer and Qdenomer

	!WRITE(*,*) "Before Eq Ratio"
	!Partial Pressures of gas-phase species in mbar
	Eq = HenrysLaw (Temperature, WhichRxn)
	P1 = GasPhaseBurden*RstarMB*Temperature/Avogadro
	P2 = GasPhaseBurden2*RstarMB*Temperature/Avogadro

	EqRatio = 1.0 / Eq / P1 / P2
	
	!WRITE(*,*) "Initial EqRatio: ", EqRatio
		
	IF  ((ABS(EqRatio-1) .LT. 0.005)) THEN
		ReturnType = 2
		RETURN
	END IF

	Qnumer = NewSolidSum !mol/cm3
	Qdenom = MIN(GasPhaseBurden,GasPhaseBurden2) !molecules/cm3
	Qdenom = Qdenom / Avogadro  !mol/cm3
	!WRITE(*,*) "Type 6", Qnumer, Qdenom

	!WRITE(*,*) "Check 4"
	!! I counts the number of iterations
	I = 1

	!! STEP 2: Calculate the Mass Step Sizes by which to Correct Concentrations
	Z  = (Qdenom + Qnumer) / 2.
	dX = Qdenom - Z				!! The Mass Flux Factor
	!WRITE(*,*) "dx: ", dX, " Z: ", Z
	!! Loop over Corrections until convergence if there is anything to do
	IF (.NOT. (Qnumer .EQ. 0 .AND. Qdenom .EQ. 0)) THEN
	DO WHILE ((ABS(EqRatio-1) .GT. ErrorTolerance))

	!WRITE(*,*) "Check 5"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 3: Adjust Molalities: !!
	!!			      !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	NewSolidSum	= NewSolidSum + dX
	IF (DissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems)		&
		 GasPhaseBurden  = GasPhaseBurden  - dX * Avogadro 
	IF (DissolutionData(WhichRxn,2) .LE. HowManyEvolveGasChems)		&
		 GasPhaseBurden2 = GasPhaseBurden2 - dX * Avogadro 

	!WRITE(*,*) GasPhaseBurden, GasPhaseBurden2
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) "Check 7"
	Z = 0.5*Z

	!! The recalculation is based on the equilibrium ratio
	Eq = HenrysLaw (Temperature, WhichRxn)
	P1 = GasPhaseBurden*RstarMB*Temperature/Avogadro
	P2 = GasPhaseBurden2*RstarMB*Temperature/Avogadro

	EqRatio = 1.0 / Eq / P1 / P2
	
	IF (EqRatio .GE. 1.+ErrorTolerance) dX = -1.*Z
	IF (EqRatio .LE. 1.-ErrorTolerance) dX = Z
	!WRITE(*,*) EqRatio, Z, dX

	!! Count the number of attempts before convergence
	I = I+1

	!Exit if all solid evaporates
	IF (INT(DissolutionData(WhichRxn,11)) .EQ. 6 .AND. NewSolidSum .LT. 1.0e-27) THEN !Less than a molecule per cm3
		EXIT
	END IF

	!! Just exit if no convergence (MJA, 100507)
	IF (I .GT. 1000) THEN
		EXIT
	END IF

	!WRITE(*,*) "Check 11"
	END DO ; END IF
	
	! DeltaSolid is the ammount of the solid species that
	! must change phases to reach equilibrium. It is positive if
	! the gas-phase species are condensing, and negative if they are evaporating.
	DeltaSolid = NewSolidSum - SolidSum

	IF (DeltaSolid .GE. 0.0) THEN !Condensing

		!Check to make sure fractions sum to 1
		Total = 0.0
		DO I = 1, NumBins
			Total = Total + SurfArea(I)/AreaSum
		END DO
		IF (Total .LT. 0.999999 .OR. Total .GT. 1.000001) THEN
			CALL ERROR("Error in EquilibrateGasSolidReactions(). Area fractions do not sum to 1.")
		END IF
		
		
		CurrentParticle => particles%first
		I = 1
		20 if (associated(currentparticle)) then
    								
			FractionofMass = SurfArea(I)/AreaSum

			!Add solid
			IF(CurrentParticle%NumberofParticles .GT. 0.0) THEN
			CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) + DeltaSolid*FractionofMass/CurrentParticle%NumberofParticles
			END IF

			I = I + 1

			if (associated(currentparticle%next)) then
				currentparticle => currentparticle%next
				goto 20
			end if
		end if 
		
	ELSE IF (DeltaSolid .LT. 0.0) THEN !Evaporating

		!Check if all solid is to evaporate
		IF (NewSolidSum .LT. 1e-27) THEN
			CurrentParticle => particles%first
			I = 1
			2 if (associated(currentparticle)) then

				!Remove all solid
				CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) = 0.0
				I = I + 1

				if (associated(currentparticle%next)) then
					currentparticle => currentparticle%next
					goto 2
				end if
			end if 
			RETURN
		END IF !All Evaporating
				
		
		!If some of the particles evaporated all their solid before equilibrating,
		!we need to look at the particles that still have solid to supply the rest
		!needed to reach equilibrium.
		
		RemainingSolidtoEvaporate = -1.0*DeltaSolid
		J = 0
		
		DO WHILE (RemainingSolidtoEvaporate .GT. 1e-30)

			!Calculate the total surface area for the particles that still
			!have the solid
			CurrentParticle => particles%first
			AreaSum = 0				
			I = 1
			3 if (associated(currentparticle)) then
    								
				!Particle has some solid left
				IF(CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.) THEN
					
					AreaSum = AreaSum + SurfArea(I)

				END IF
				
				if (associated(currentparticle%next)) then
					currentparticle => currentparticle%next
					I = I+ 1
					goto 3
				end if
			end if 

			!Check to make sure fractions sum to 1
			Total = 0.0
			
			CurrentParticle => particles%first
			I = 1
			111 if (associated(currentparticle)) then
				
				IF(CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.) THEN
					Total = Total + SurfArea(I)/AreaSum
				END IF

				I = I + 1

				if (associated(currentparticle%next)) then
					currentparticle => currentparticle%next
					goto 111
				end if
			end if 
			
			IF (Total .LT. 0.999999 .OR. Total .GT. 1.000001) THEN
				CALL ERROR("Error in EquilibrateGasSolidReactions(). Area fractions do not sum to 1.")		
			END IF
			
			CurrentParticle => particles%first
			I = 1
			1 if (associated(currentparticle)) then

				IF(CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) .GT. 0.) THEN
					FractionofMass = SurfArea(I)/AreaSum
				ELSE
					FractionofMass = 0.0 !For particles with no solid
				END IF
				
				!Evaporate solid, checking to make sure it does not go below 0
				ZeroCheck = CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) &
																   - RemainingSolidtoEvaporate*FractionofMass/CurrentParticle%NumberofParticles
				IF (ZeroCheck .GE. 0.0) THEN 
					CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) = ZeroCheck
					RemainingSolidtoEvaporate = RemainingSolidtoEvaporate - RemainingSolidtoEvaporate*FractionofMass
				ELSE
					RemainingSolidtoEvaporate = RemainingSolidtoEvaporate - CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3)))*CurrentParticle%NumberofParticles 
					CurrentParticle%AqChems(INT(DissolutionData(WhichRxn,3))) = 0.0
				END IF	 
				
				I = I + 1

				if (associated(currentparticle%next)) then
					currentparticle => currentparticle%next
					goto 1
				end if
			end if 
			
			J = J+1
			IF (J .GT. 1000) THEN
				WRITE(*,*) RemainingSolidtoEvaporate
				CALL ERROR("Infinite loop in EquilibrateGasSolidReactions()")
			END IF
		
		END DO
 	
	END IF
	
	RETURN		
		
END SUBROUTINE EquilibrateGasSolidReactions


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! EQUILIBRATE ALL OF THE PARTICLES AT A GRIDPOINT WITH THE	   !!
!! GAS PHASE, INCLUDING ORGANICS.				   !!
!!								   !!
!! The default is to not include water equilibration in this       !!
!! routine                                                         !!
!!                                                                 !!
!! ForceWaterEquilibrium will make the routine use the equilibrium !!
!! approach no matter whether the RH is above the threshold that   !!
!! would tell the routine to use the integration approach.         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateGridPointAll(EquilibrateWater, &
                ForceWaterEquilibrium, WaterOpenSystem)

	USE Aerosols, ONLY :    Particles,		&
			        Particle,		&
				ParticleArray,		&
				KusikMeissner,		&
				FindElectrolyteEquilibrium,		&
				ReformElectrolyte,	&
				RecalculateRadius,	&
				SurfaceTension,		&
				DumpParticleContentsAtError,		&
				SortAerosolAtGridPointForCondensation,  &
				PrintAqEqReaction


	USE GridPointFields, ONLY : GetTemp,		        &
				    GetGridCellChemBurden,	&
				    ReplaceGridCellChemBurden,	&
				    GetRelativeHumidity,        &
                                    GetSatVapBurden

	USE Chemistry,	 ONLY : HowManyEvolveGasChems, &
				HowManyAqEqReactions, &
				AqEquilibriaList, &
				HowManyAqChems, &
				HowManyAqCations, &
				HowManyAqAnions

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterations,	&
					AqThermoEquilibriumError, &
					ThermoBulkMode, &
					WaterContentPrecision, &
					DoDissolutionFlag, &
					DoHydrophobicOrgDissolutionFlag, &
					DoHydrophilicOrgDissolutionFlag, &
					ZSROrgWater, Metastable

	USE Time, ONLY : CurrentTime, BeginningTime

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR

	IMPLICIT NONE

	!! External Variables
	LOGICAL :: EquilibrateWater, ForceWaterEquilibrium, WaterOpenSystem

	!! Internal Variables
	LOGICAL :: Equilibrated, RxnEquilibrated, LocalEquilibrateWater
	LOGICAL :: LocalForceWaterEquilibrium, LocalWaterOpenSystem
	REAL*8  :: GasPhaseBurden, GasPhaseBurden2, Temperature, &
                   II, TotalMoles, WaterSatBurden
	REAL*8  :: S1MassFracs(2), S2MassFracs(6), S3MassFracs(2), &
		S4MassFracs(3), S5MassFracs(2), S6MassFracs(3), &
		S7MassFracs(8), S8MassFracs(4), S9MassFracs(4), &
		S10MassFracs(4), InitialBurden(10), OrgWaterContent, InitWater
	INTEGER :: I, J, K, L, M, Q, HowManyLinks, ReturnType, OrgRxnIndex, &
                  AqRxnIndex, StoreReturn
	TYPE(Particle),POINTER :: Current

	LocalEquilibrateWater = EquilibrateWater

	LocalForceWaterEquilibrium = ForceWaterEquilibrium

	LocalWaterOpenSystem = WaterOpenSystem

	!! Sort the aerosol, pushing all empty sections to the end
	CALL SortAerosolAtGridPointForCondensation ()

	!! Count the number of links
	HowManyLinks = 0
	Current => Particles%First
	IF (ASSOCIATED(Current)) THEN
10		IF (Current%NumberOfParticles .GT. 1e-6 ) &
			HowManyLinks = HowManyLinks + 1
			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 10
			END IF
	END IF

	!Equilibrate gas-phase sulfate (using mass-transfer rate to remove non-uniqueness)
	IF(DoDissolutionFlag) CALL EquilibrateSulfate()
		
	!Calculate surrogate concentrations and 
	!mass fractions of organics contributing to surrogates, and
	!Store initial surrogate gas phase burdens in a dummy vector
	!IF(DoHydrophobicOrgDissolutionFlag.OR.DoHydrophilicOrgDissolutionFlag)&
	!	CALL CalculateSurrogateGasConc(S1MassFracs, S2MassFracs, &
        !        S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	!	S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, &
        !        InitialBurden)

	Temperature = GetTemp()

	Equilibrated = .FALSE.
	K = 0

	!! The main loop calculates whether it's equilibrated
	DO WHILE (.NOT.Equilibrated)

		K = K + 1
		!WRITE(*,*) "Number of Interations of EquilibrateGridPointAll", K
		Equilibrated = .TRUE.

		!! Loop over the reactions
		DO I = 1, HowManyDissolutionReactions+1+HowManyOrganicDissolutionReactions+HowManyAqOrganicDissolutionReactions			
		
			L = 0
			RxnEquilibrated = .FALSE.
			
			!First, equilibrate the inorganic reactions
			IF (I .LE. HowManyDissolutionReactions + 1 .AND. DoDissolutionFlag) THEN
				
				!Get the gas phase concentration from the grid
				IF (I .LE. HowManyDissolutionReactions) THEN
					GasPhaseBurden = GetGridCellChemBurden ( INT(DissolutionData(I,1)))
					IF (DissolutionData(I,11) .EQ. 6.) GasPhaseBurden2 = GetGridCellChemBurden ( INT(DissolutionData(I,2)))
				
					!Equilibrate gas-solid reactions (Type 6) using surface area to remove non-uniqueness problem
					IF (DissolutionData(I,11) .EQ. 6.) THEN
						CALL EquilibrateGasSolidReactions( GasPhaseBurden, I, Temperature, ReturnType, &
										   GasPhaseBurden2 = GasPhaseBurden2, OptErrorTolerance=AqThermoEquilibriumError/2.)
						IF (ReturnType .NE. 2) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF
				
					END IF
				
				END IF
			
				
			
				!! Loop over all of the particles until 
				!! this particular reaction is happy
				DO WHILE (.NOT.RxnEquilibrated)

					RxnEquilibrated = .TRUE.

					Current => Particles%First
					DO J = 1, HowManyLinks
									
						!Skip bins with too few particles
						DO WHILE(Current%numberofparticles .LT. 1.0e-6)
							Current => Current%Next
						END DO

						!Count total moles of water soluble material in particle
						TotalMoles = 0.0
						DO Q = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
							TotalMoles = TotalMoles + Current%AqChems(Q)
						END DO
					
						IF (I .EQ. 1) THEN

							!! If water should be equilibrated, equilibrate it
							!! equilibrate the first few times through but don't let it stagnate
							!! (This step is ripe for oscillating non-convergence)
							IF (LocalEquilibrateWater) THEN
								!WRITE(*,*), "Before EqAllWater"
								CALL EqAllWater(Current, ReturnType = ReturnType,  &
												   OpenSystem = LocalWaterOpenSystem)
                                                                !WRITE(*,*), "After EqAllWater"
							ELSE
								ReturnType = 2
							END IF

						ELSE IF (I .EQ. 2) THEN

							!!This is the hard-coded sulfate condensation reaction
							!!Ignore here, as this is handled by the subroutine EquilibrateSulfate
							ReturnType = 2

						ELSE IF (I .GT. 2 .AND. I .LE. HowManyDissolutionReactions) THEN
 
							
							!WRITE(*,*) "Rxn: ", I, "Type", DissolutionData(I,11) 
							!! If the electrolyte is attached to an electrolyte dissociation reaction,
							!! but the dissolution is based on its electrolyte, then reform the electrolyte
							!! at step one to get it out of the aerosol phase.
							!! (only if it is the initial time and step)
							IF (CurrentTime .EQ. BeginningTime		&		!! It is the initial time
								.AND. K .EQ. 1 .AND. L .EQ. 1		&		!! It is the first step
								.AND. DissolutionData(I,11) .EQ. 1	&		!! It is condensation of an electrolyte w/out dissociation
								.AND. DissolutionData(I,12) .GT. 1) THEN    !! There is a related dissociation reaction identified
								
								CALL ReformElectrolyte(INT(DissolutionData(I,12)), Current)
								
							END IF

							!If water content is too low, skip dissolution reactions
							!IF (((ThermoBulkMode .AND. Current%AqChems(1) .LT. 1.0e-17) & !Dry Bulk
							!	.OR. ( NOT(ThermoBulkMode) .AND. Current%AqChems(1) .LT. 0.0001*TotalMoles)) .AND. & !Dry Particle (mole fraction of H2O < 0.01%)
							!	DissolutionData(I,11) .NE. 6.) CYCLE    !Not a gas-solid rxn

							IF (DissolutionData(I,11) .EQ. 6.) THEN
								CYCLE !Already did gas-solid (Type 6) reactions
							
							ELSE
								!Equilibrate the reaction for one particle
								!WRITE(*,*) I, GasphaseBurden
								WaterSatBurden = GetSatVapBurden ()
                                                                !WRITE(*,*) "Before EquilibrateDissolutionReaction"
                                                                CALL EquilibrateDissolutionReaction (GasPhaseBurden, I, &
                                                                      Current, Temperature, ReturnType, &
                                                                      WaterSatBurden=WaterSatBurden, GasPhaseBurden2=0.0, OptErrorTolerance=AqThermoEquilibriumError/2.)				
			                                        !WRITE(*,*) "After EquilibrateDissolutionReaction"

							END IF

							!! Make a more refined guess if still equilibrating late in the game
							!! (That is, force it to equilibrate much closer to 1 each time through)
							IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
								WaterSatBurden = GetSatVapBurden ()
                                                                IF (DissolutionData(I,11) .EQ. 6.) THEN
									CYCLE !Already equilibrated gas-solid reactions
								ELSE
									CALL EquilibrateDissolutionReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, WaterSatBurden=WaterSatBurden, GasPhaseBurden2=0.0, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
								END IF
							END IF

						!! HMDR+1 is a dummy saying "do the aqueous rxns."
						ELSE IF (I .EQ. HowManyDissolutionReactions+1) THEN

							!! Recalculation of KM can make it much more difficult to 
							!! equilibrate (it is not formally necessarily possible).
							!! But it is worth trying to do.  

							!If water content is too low, skip dissociation reactions
							!IF (((ThermoBulkMode .AND. Current%AqChems(1) .LT. 1.0e-17) & !Dry Bulk
							!	.OR. ( NOT(ThermoBulkMode) .AND. Current%AqChems(1) .LT. 0.0001*TotalMoles))) CYCLE !Dry Particle (mole fraction of H2O < 0.01%) 

							
							!WRITE(*,*) "Calling Electrolyte Eq"
							CALL FindElectrolyteEquilibrium(Current, .TRUE., &
                                                            FirstEquilibration=.FALSE., ReturnType=ReturnType)
							
							!If Electrolyte equilibrium hasn't been reached
							!after 50 times through particle list or 5 times through a single reaction
							!, just give up on it.
							IF (L .GT. 5 .OR. K .GT. 50) THEN							
								ReturnType = 2
							END IF

						END IF

						IF (ReturnType .NE. 2) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF

						!! If not the aqueous electrolyte equilibration, then update the thermo
						IF (I .LE. HowManyDissolutionReactions) CALL KusikMeissner (Current)

						L = L+1
						Current => Current%Next
					END DO !Particles

				

					IF (L .GT. AqThermoNumbAllowedIterations) THEN
						CALL DumpParticleContentsAtError(Particles%First, &
                                                     "UnEquilableParticle.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.)

						CALL ERROR ("EquilibrateGridPoint() in CondensationRelatedFunctions.h tried to equilibrate all "     // &
					            "of the particles (stored as "//TRIM(INT2STR(HowManyLinks))//" links) at a grid point "  // &
								"to dissolution reaction "//TRIM(INT2STR(I))//" (first is water, the others are indexed "// &
								"as input in Dissolution.in) and failed.  The contents of the first particle in that "   // &
								"section are printed as UnEquilableParticle.txt")
					END IF
				END DO !RxnEq

				!Matt Changed this!
				IF (I .LE. HowManyDissolutionReactions) THEN
					!! Replace the gas phase chemical concentrations
					IF (INT(DissolutionData(I,1)) .LE. HowManyEvolveGasChems) THEN
						CALL ReplaceGridCellChemBurden ( INT(DissolutionData(I,1)), GasPhaseBurden)
					END IF

					IF (INT(DissolutionData(I,2)) .LE. HowManyEvolveGasChems .AND. DissolutionData(I,11) .EQ. 6.)	&
						CALL ReplaceGridCellChemBurden ( INT(DissolutionData(I,2)), GasPhaseBurden2)
				END IF
						
			!Second, equilibrate the hydrophobic organic reactions
			ELSE IF (I .GT. HowManyDissolutionReactions + 1 .AND. I .LE. HowManyDissolutionReactions+1+HowManyOrganicDissolutionReactions &
			         .AND. DoHydrophobicOrgDissolutionFlag) THEN
				
				!Get the proper index for the organic dissolution data
				OrgRxnIndex = I - (HowManyDissolutionReactions + 1)
				!WRITE(*,*) "Org Rxn #: ", OrgRxnIndex

				!Get the gas-phase concentrations from the grid
				GasPhaseBurden = GetGridCellChemBurden ( INT(OrganicDissolutionData(OrgRxnIndex,1)))
				
				!! Loop over all of the particles until 
				!! this particular reaction is happy
				DO WHILE (.NOT.RxnEquilibrated)

					RxnEquilibrated = .TRUE.

					Current => Particles%First
					DO J = 1, HowManyLinks
						
						!Skip bins with too few particles
						DO WHILE(Current%numberofparticles .LT. 1.0e-6)
							Current => Current%Next
						END DO			
						
						!Equilibrate one hydrophobic organic reaction for one particle
                                                !WRITE(*,*) "Before Org Eq"
						CALL EquilibrateHydrophobicReaction (GasPhaseBurden, OrgRxnIndex, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/2.)							
					         !WRITE(*,*) "After Org Eq"
						!! Make a more refined guess if still equilibrating late in the game
						!! (That is, force it to equilibrate much closer to 1 each time through)
						IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
							CALL EquilibrateHydrophobicReaction (GasPhaseBurden, OrgRxnIndex, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
						END IF

						!WRITE(*,*) "Past Equil"
						IF (ReturnType .NE. 2 ) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF

						!Recalculate activity coefficients
						!WRITE(*,*) "Before UNIFAC"
                                                CALL UpdateHydrophobicUNIFAC(Current, Temperature)
                                                !WRITE(*,*) "After UNIFAC"
							L = L+1
					Current => Current%Next
					END DO

					IF (L .GT. AqThermoNumbAllowedIterations) THEN
						CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
						CALL ERROR ("EquilibrateGridPointAll() in CondensationRelatedFunctions.h tried to equilibrate all "     // &
					            "of the particles (stored as "//TRIM(INT2STR(HowManyLinks))//" links) at a grid point "  // &
								"to hydrophobic organic dissolution reaction "//TRIM(INT2STR(OrgRxnIndex))//" (As indexed "// &
								"as input in OrganicDissolution.in) and failed.  The contents of the first particle in that "   // &
								"section are printed as UnEquilableParticle-Hydrophobic.txt")
					END IF
				END DO
			
				!! Replace the surrogate gas phase chemical concentrations
				IF (INT(OrganicDissolutionData(OrgRxnIndex,1)) .LE. HowManyEvolveGasChems)										&
					CALL ReplaceGridCellChemBurden ( INT(OrganicDissolutionData(OrgRxnIndex,1)), GasPhaseBurden)			
		
			ELSE IF (I .GT. HowManyDissolutionReactions+1+HowManyOrganicDissolutionReactions &
					.AND. DoHydrophilicOrgDissolutionFlag) THEN

				AqRxnIndex = I - (HowManyDissolutionReactions+1+HowManyOrganicDissolutionReactions)
				GasPhaseBurden = GetGridCellChemBurden ( INT(AqOrganicDissolutionData(AqRxnIndex,1)))

				!WRITE(*,*) "Aq Org Rxn #: ", AqRxnIndex

				!! Loop over all of the particles until 
				!! this particular reaction is happy
				DO WHILE (.NOT.RxnEquilibrated)

					RxnEquilibrated = .TRUE.

					Current => Particles%First
					DO J = 1, HowManyLinks
						
						!Skip bins with too few particles
						DO WHILE(Current%numberofparticles .LT. 1.0e-6)
							Current => Current%Next
						END DO			
						
						!Equilibrate one reaction
						CALL EquilibrateHydrophilicReaction (GasPhaseBurden, AqRxnIndex, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/2.)							

						!! Make a more refined guess if still equilibrating late in the game
						!! (That is, force it to equilibrate much closer to 1 each time through)
						IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
							CALL EquilibrateHydrophilicReaction (GasPhaseBurden, AqRxnIndex, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
						END IF

						IF (ReturnType .NE. 2 ) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF

						!Recalculate activity coefficients
						!WRITE(*,*) "Before Update"
						CALL UpdateHydrophilicUNIFAC(Current, Temperature)
						!WRITE(*,*) "After Update"

							L = L+1
					Current => Current%Next
					END DO

					IF (L .GT. AqThermoNumbAllowedIterations) THEN
						CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
						CALL ERROR ("EquilibrateGridPointAll() in CondensationRelatedFunctions.h tried to equilibrate all "     // &
					            "of the particles (stored as "//TRIM(INT2STR(HowManyLinks))//" links) at a grid point "  // &
								"to hydrophilic organic dissolution reaction "//TRIM(INT2STR(AqRxnIndex))//" (As indexed "// &
								"as input in HydrophilicOrganicDissolution.in) and failed.  The contents of the first particle in that "   // &
								"section are printed as UnEquilableParticle-Hydrophilic.txt")
					END IF
				END DO
			
				!! Replace the surrogate gas phase chemical concentrations
				IF (INT(AqOrganicDissolutionData(AqRxnIndex,1)) .LE. HowManyEvolveGasChems)										&
					CALL ReplaceGridCellChemBurden ( INT(AqOrganicDissolutionData(AqRxnIndex,1)), GasPhaseBurden)			
					
			END IF	 
			
		END DO !Rxn #
			
		IF (K .GT. AqThermoNumbAllowedIterations) THEN
			CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
			CALL ERROR ("EquilibrateGridPoint() in CondensationRelatedFunctions.h tried to equilibrate all of the dissolution "// &
					    "reactions, aerosol water content, and aerosol ionic content and failed.  Usually this implies an "// &
						"unrealistic starting condition.  You may try to fix this by increasing AqThermoNumbAllowedIterations,"// &
						" the maximum  allowed iterations of this type of routine.  Another problem may be an oscillation of "// &
						"water content and a particular dissolving species.")
		END IF

	END DO !NotinEq

	!Update Gas Phase Concentrations of individual organics
	!IF(DoHydrophobicOrgDissolutionFlag .OR. DoHydrophilicOrgDissolutionFlag) &
	!	CALL UpdateGasPhaseOrganics( S1MassFracs, S2MassFracs, &
        !    S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	!		S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)	
	
	RETURN
END SUBROUTINE EquilibrateGridPointAll

!!This subroutine finds the total water equilibrium for a particle.
!!Note that it forces each particle to contain at least 600 molecules of H2O.
!! (c) Matt Alvarado, 2006
SUBROUTINE EqAllWater(Current, ReturnType, OpenSystem)

	USE Aerosols, ONLY : Particles,					&
				Particle,				&
				ParticleArray,				&
				EquilibriumWaterContentAmount,     &
				DumpParticleContentsAtError, &
				FindElectrolyteEquilibrium

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterations,	&
					AqThermoEquilibriumError, &
					WaterContentPrecision, &
					Avogadro, &
					MinimumWater

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR

	USE GridPointFields, ONLY : GetTemp,				&
					GetGridCellChemBurden,		&
					ReplaceGridCellChemBurden, &
					GetRelativeHumidity, &
					GetGasBurden	
	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: Current
	INTEGER :: ReturnType
	LOGICAL:: OpenSystem

	!! Internal Variables
	LOGICAL :: Equilibrated
	REAL*8  :: OrgWaterContent, InorgWaterContent, CalcWater, EqRatio
	REAL*8  :: GasWaterContent, Temperature, RH, DeltaWater
	REAL*8  :: Qnumer, Qdenom, dX, Z_MFI, NewBurden, ErrorTolerance, &
				Upper, Lower, Initial, Final
	INTEGER :: I
	LOGICAL :: ForceEquilibrium, Scaffolding

	Scaffolding = .FALSE.

	IF(Scaffolding) WRITE(*,*) "Check 1"
	IF(Current%NumberofParticles .LE. 0.) THEN
		Current%AqChems(1) = 0.
		RETURN
	END IF
	
	Temperature = GetTemp()
	RH = GetRelativeHumidity()

	IF(Scaffolding) WRITE(*,*) "Check 2"
	!Calculate equilibrium orgabnic water content at this RH
	OrgWaterContent = HydrophilicOrganicWaterContent(Current, Temperature, RH)
	
	  IF(Scaffolding) WRITE(*,*) "Check 3", OrgWaterContent
	!Calculate equilibrium inorganic water content at this RH
	ForceEquilibrium = .TRUE.
	CALL EquilibriumWaterContentAmount (Current, RH, ReturnType, ForceEquilibrium, InorgWaterContent)
	IF(Scaffolding) WRITE(*,*) "Check 3a", InorgWaterContent

	!Use ZSR approx to calculate total eq water content
	CalcWater = OrgWaterContent + InorgWaterContent

	!Calculate equilibrium ratio
	EqRatio = CalcWater/Current%AqChems(1)
	!WRITE(*,*) "Input", EqRatio

	IF(Scaffolding) WRITE(*,*) "Check 4"
	!Check if we are already at equilibrium
	IF  (ABS(EqRatio-1.) .LT. 0.001) THEN
		ReturnType = 2 !Says was already in eq.
		RETURN
	END IF

	!WRITE(*,*) "Water Check", RH, Current%AqChems(1), EqRatio
	!Return if water content is approaching zero (fix at 5.0e-21 moles/particle).
        !Minimum embryo radius 2.8 nanometers.
	IF(Current%AqChems(1) .LT. MinimumWater .AND. EqRatio .LT. 1.0) THEN
		!WRITE(*,*) "MinimumWater! RH = ", RH 
		ReturnType = 2 !Says was already in eq.
		RETURN
	END IF

	!If not already in equilibrium, calculate new water content
	IF(OpenSystem) THEN 
		!! OPEN SYSTEM
		!! Calculate and return the total water content 
		!! for the particle. This is CalcWater
		
		!WRITE(*,*) "Open System! Hi Matt!"
		Current%AqChems(1) = CalcWater
		ReturnType = 1 !Now in equilibrium
		RETURN
	
	ELSE 	
		!! CLOSED SYSTEM
		!Just make current water equal calc water and return
		GasWaterContent = GetGasBurden (1)/Avogadro/Current%NumberofParticles		
		!Scaled to be mol H2O in gas phase per particle
				
		DeltaWater = CalcWater - Current%AqChems(1)
		!If CalcWater is below minimum, return minimum
		IF(CalcWater .LT. 0.99*MinimumWater) DeltaWater = 0.99*MinimumWater - Current%AqChems(1)
		
		Current%AqChems(1) = Current%AqChems(1) + DeltaWater		
		GasWaterContent = GasWaterContent - DeltaWater
		NewBurden = GasWaterContent*Avogadro*Current%NumberofParticles
		CALL ReplaceGridCellChemBurden ( 1, NewBurden)
			
		!Get new RH
		RH = GetRelativeHumidity()

		!Calculate equilibrium organic water content at this RH
		OrgWaterContent = HydrophilicOrganicWaterContent(Current, Temperature, RH)
	
		!Calculate equilibrium inorganic water content at this RH
		CALL EquilibriumWaterContentAmount (Current, RH, ReturnType, ForceEquilibrium, InorgWaterContent)	
			
		!Add total eq water content
		CalcWater = OrgWaterContent + InorgWaterContent

		!Calculate equilibrium ratio
		EqRatio = Current%AqChems(1)/CalcWater

		I = I+1
		
		GasWaterContent = GetGasBurden (1)/Avogadro/Current%NumberofParticles		

		ReturnType = 1	
		RETURN
	
	END IF

END SUBROUTINE EqAllWater

!! This subroutine only equilibrates all the particles internally.
!! In this routine, only water is allowed to go between 
!! the gas and particle phases.
!! The default is not to equilibrate water, but this can be changed
!! by setting "EquilibrateWater" to .TRUE.
!! Matt Alvarado, 08/28/06
SUBROUTINE EquilibrateInternallyatGridPoint(EquilibrateWater, WaterOpenSystem)

	USE Aerosols,    ONLY : Particles,		&
				Particle,		&
				ParticleArray,		&
				KusikMeissner,		&
				FindElectrolyteEquilibrium,	&
				ReformElectrolyte,	&
				RecalculateRadius,	&
				SurfaceTension,		&
				DumpParticleContentsAtError,	&
				SortAerosolAtGridPointForCondensation, &
				PrintAqEqReaction

	USE GridPointFields, ONLY : GetTemp,		&
				    GetGridCellChemBurden,	&
				    ReplaceGridCellChemBurden,	&
				    GetRelativeHumidity

	USE Chemistry,		 ONLY : HowManyEvolveGasChems, &
					HowManyAqEqReactions, &
					AqEquilibriaList, &
					HowManyAqChems, &
					HowManyAqCations, &
					HowManyAqAnions

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterations,	&
					AqThermoEquilibriumError, &
					ThermoBulkMode, &
					WaterContentPrecision, &
					DoDissolutionFlag, &
					DoHydrophobicOrgDissolutionFlag, &
					DoHydrophilicOrgDissolutionFlag, &
					ZSROrgWater

	USE Time, ONLY : CurrentTime, BeginningTime

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR,Transcript

	IMPLICIT NONE

	!! External Variables
	LOGICAL :: EquilibrateWater, WaterOpenSystem

	!! Internal Variables
	LOGICAL :: Equilibrated, RxnEquilibrated, LocalEquilibrateWater
	LOGICAL :: LocalWaterOpenSystem, ParticleEquilibrated
	REAL*8  :: Temperature, II
	INTEGER :: I, K, HowManyLinks, ReturnType
	REAL*4  :: etimearray(2), etime
	TYPE(Particle),POINTER :: Current

	LocalEquilibrateWater = EquilibrateWater
	
	LocalWaterOpenSystem = WaterOpenSystem
	
	!! Sort the gridcell, pushing all empty sections to the end
	CALL SortAerosolAtGridPointForCondensation ()
	!! Count the number of links
	HowManyLinks = 0
	Current => Particles%First
	IF (ASSOCIATED(Current)) THEN
10		IF (Current%NumberOfParticles .GT. 0. ) &
			HowManyLinks = HowManyLinks + 1
			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 10
			END IF
	END IF

	Temperature = GetTemp()

	Equilibrated = .FALSE.
	K = 0

	!! The main loop calculates whether it's equilibrated
	DO WHILE (.NOT.Equilibrated)

		K = K + 1
		Equilibrated = .TRUE.

		Current => particles%first
		ParticleEquilibrated = .TRUE.
		DO I = 1, HowManyLinks

                        !Skip dry partciles (BC only)
                        IF (Current%Dry) THEN
                            if (associated(current%next)) then
				current => current%next
			    end if
                            CYCLE
                        ENDIF
  
			!Water Eq.
			CALL EqAllWater(Current, ReturnType, LocalWaterOpenSystem)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ii = etime(etimearray)
!	     CALL TRANSCRIPT("")
!	     CALL TRANSCRIPT("CPU time after Water Eq "//trim(real2str(ii)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			IF(ReturnType .NE. 2) ParticleEquilibrated = .FALSE.
			!WRITE(*,*) "WaterReturnType", ReturnType, Current%AqChems(1)	
				
			CALL FindElectrolyteEquilibrium (Current, UpdateThermo=.TRUE., FirstEquilibration=.TRUE., ReturnType = ReturnType)
			!Give up on inorganic eq after 10 iterations
			IF(ReturnType .NE. 2 .AND. K .LE. 5) ParticleEquilibrated = .FALSE.
			!WRITE(*,*) "InOrganicReturnType", ReturnType	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ii = etime(etimearray)
!	     CALL TRANSCRIPT("")
!	     CALL TRANSCRIPT("CPU time after Electrolyte Eq "//trim(real2str(ii)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!				
			CALL EquilibrateOrganicParticle (Current, Temperature, ReturnType, OptErrorTolerance=AqThermoEquilibriumError/2.)
			IF(ReturnType .NE. 2) ParticleEquilibrated = .FALSE.
			!WRITE(*,*) "OrganicReturnType", ReturnType	
			
			IF(NOT(ParticleEquilibrated)) Equilibrated = .FALSE.
			if (associated(current%next)) then
				current => current%next
			end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ii = etime(etimearray)
!	     CALL TRANSCRIPT("")
!	     CALL TRANSCRIPT("CPU time after Org Eq "//trim(real2str(ii)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		END DO
		
		IF (K .GT. 200) THEN
			CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
			CALL ERROR ("EquilibrateInternallyatGridPoint() in CondensationRelatedFunctions.h tried to equilibrate "// &
					    "aerosol water content and aerosol internal equilibrium and failed.  Usually this implies an "// &
						"unrealistic starting condition.  You may try to fix this by increasing"// &
						" the maximum  allowed iterations of this type of routine.  Another problem may be an oscillation of "// &
						"water content.")
		END IF

	END DO !NotinEq

	!Recalculate Radii of particles (including empty ones)
	Current => particles%first
	DO WHILE (associated(current))
		CALL RecalculateRadius(Current)
		current => current%next
	END DO

	RETURN
END SUBROUTINE EquilibrateInternallyatGridPoint


REAL*8 FUNCTION RidderMethod(func,x1,x2,xacc)
!!This is a root-finder based on Ridder's Method for root finding
!!from Sec. 9.2, pp. 351-352 of 
!!"Numerical Recipies in FORTRAN: The Art of Scientific Computing"
!! 2nd ed. by W.H. Press et al. (ISBN:0 521 43064 X)
!! x1 and x2 are boundaries around the root, while xacc is the required accuracy
!! func is the function we are trying to zero
	USE InfrastructuralCode, ONLY : ERROR
	
	IMPLICIT NONE

	!External variables
	REAL*8 :: x1, x2, xacc, func
	EXTERNAL func
	 
	!Internal variables
	INTEGER :: J
	REAL fh,fl,fm,fnew,s,xh,xl,xm,xnew

	fl = func(x1)
	fh = func(x2)
	IF ((fl .GT. 0. .AND. fh .LT. 0.) .OR. (fh .GT. 0. .AND. fl .LT. 0.)) THEN
		xl = x1
		xh = x2
		RidderMethod = -1.0e30 !An unlikely dummy value
		DO J = 1, 100 !Only 100 iterations allowed!
			xm = 0.5*(xl+xh)
			fm = func(xm)
			s = sqrt(fm**2-fl*fh)
			IF (s .eq. 0.) RETURN
			xnew = xm + (xm-xl)*(sign(1., fl-fh)*fm/s)
			IF (abs(xnew-RidderMethod) .le. xacc) RETURN
			RidderMethod = xnew
			fnew = func(RidderMethod)
			IF (fnew .eq. 0.) RETURN
			IF (sign(fm,fnew) .ne. fm) THEN
				xl = xm
				fl = fm
				xh = RidderMethod
				fh = fnew
			ELSE IF (sign(fl,fnew) .ne. fl) THEN
				xh = RidderMethod
				fh = fnew
			ELSE IF (sign(fh,fnew) .ne. fh) THEN
				xl = RidderMethod
				fl = fnew
			ELSE
				CALL ERROR("RidderMethod - this shouldn't happen.")
			END IF
			IF (abs(xh-xl) .le. xacc) RETURN
		END DO
		CALL ERROR("Ridder Method exceeded the maximum number of iterations")
	
	ELSE IF (fl .EQ. 0.) THEN
		RidderMethod = x1
	ELSE IF (fh .EQ. 0.) THEN
		RidderMethod = x2
	ELSE
		CALL ERROR("No Root found in brackets in RidderMethod.")
	END IF
	RETURN

END FUNCTION RidderMethod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! EQUILIBRATE ALL OF THE PARTICLES AT A GRIDPOINT WITH THE    !!
!! GAS PHASE, NH3, HNO3 and HCl only, no water or H2SO4.       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateInorgGridPoint()

	USE Aerosols,    ONLY : Particles,		&
			        Particle,		&
				ParticleArray,		&
				KusikMeissner,		&
				FindElectrolyteEquilibrium,	&
				ReformElectrolyte,	&
				RecalculateRadius,	&
				SurfaceTension,			&
				DumpParticleContentsAtError,	&
				SortAerosolAtGridPointForCondensation, &
				PrintAqEqReaction


	USE GridPointFields, ONLY : GetTemp,		&
				    GetGridCellChemBurden,		&
				    ReplaceGridCellChemBurden,	&
				    GetRelativeHumidity, &
                                    GetSatVapBurden

	USE Chemistry,		 ONLY : HowManyEvolveGasChems, &
					HowManyAqEqReactions, &
					AqEquilibriaList, &
					HowManyAqChems, &
					HowManyAqCations, &
					HowManyAqAnions

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterations,	&
					AqThermoEquilibriumError, &
					ThermoBulkMode, &
					WaterContentPrecision, &
					DoDissolutionFlag, &
					DoHydrophobicOrgDissolutionFlag, &
					DoHydrophilicOrgDissolutionFlag, &
					ZSROrgWater, Metastable

	USE Time, ONLY : CurrentTime, BeginningTime

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR

	IMPLICIT NONE

	!! Internal Variables
	LOGICAL :: Equilibrated, RxnEquilibrated, LocalEquilibrateWater
	LOGICAL :: LocalForceWaterEquilibrium, LocalWaterOpenSystem
	REAL*8  :: GasPhaseBurden, GasPhaseBurden2, Temperature, XX(3), II, TotalMoles, WaterSatBurden
	INTEGER :: I, J, K, L, M, Q, HowManyLinks, ReturnType, OrgRxnIndex, AqRxnIndex, StoreReturn
	TYPE(Particle),POINTER :: Current

	!! Sort the gridcell, pushing all empty sections to the end
	CALL SortAerosolAtGridPointForCondensation ()

	!! Count the number of links
	HowManyLinks = 0
	Current => Particles%First
	IF (ASSOCIATED(Current)) THEN
10		IF (Current%NumberOfParticles .GT. 1e-6 ) &
			HowManyLinks = HowManyLinks + 1
			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 10
			END IF
	END IF
		
	Temperature = GetTemp()

	Equilibrated = .FALSE.
	K = 0

	!! The main loop calculates whether it's equilibrated
	DO WHILE (.NOT.Equilibrated)

		K = K + 1
		
		Equilibrated = .TRUE.

		!! Loop over the reactions
		DO I = 3, HowManyDissolutionReactions+1			
			L = 0
			RxnEquilibrated = .FALSE.
			!WRITE(*,*) "Reaction: ", I
			!First, equilibrate the inorganic reactions
			IF (I .LE. HowManyDissolutionReactions + 1 .AND. DoDissolutionFlag) THEN
				
				!Get the gas phase concentration from the grid
				IF (I .LE. HowManyDissolutionReactions) THEN
					GasPhaseBurden = GetGridCellChemBurden ( INT(DissolutionData(I,1)))
					IF (DissolutionData(I,11) .EQ. 6.) GasPhaseBurden2 = GetGridCellChemBurden ( INT(DissolutionData(I,2)))
				
					!Equilibrate gas-solid reactions (Type 6) using surface area to remove non-uniqueness problem
					IF (DissolutionData(I,11) .EQ. 6.) THEN
						WaterSatBurden = GetSatVapBurden ()
                                                CALL EquilibrateGasSolidReactions( GasPhaseBurden, I, Temperature, ReturnType, &
										   OptErrorTolerance=AqThermoEquilibriumError/10., &
										   GasPhaseBurden2 = GasPhaseBurden2)
						IF (ReturnType .NE. 2) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF
				
					END IF
				
				END IF
			
				!! Loop over all of the particles until 
				!! this particular reaction is happy
				DO WHILE (.NOT.RxnEquilibrated)

					RxnEquilibrated = .TRUE.

					Current => Particles%First
					DO J = 1, HowManyLinks
									
						!Skip bins with too few particles
						DO WHILE(Current%numberofparticles .LT. 1.0e-6)
							Current => Current%Next
						END DO

						!WRITE(*,*) "Particle: ", J
						!Count total moles of water soluble material in particle
						TotalMoles = 0.0
						DO Q = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
							TotalMoles = TotalMoles + Current%AqChems(Q)
						END DO
					
						IF (I .GT. 2 .AND. I .LE. HowManyDissolutionReactions) THEN
 
							
							!WRITE(*,*) "Rxn: ", I, "Type", DissolutionData(I,11) 
							!! If the electrolyte is attached to an electrolyte dissociation reaction,
							!! but the dissolution is based on its electrolyte, then reform the electrolyte
							!! at step one to get it out of the aerosol phase.
							!! (only if it is the initial time and step)
							IF (CurrentTime .EQ. BeginningTime		&		!! It is the initial time
								.AND. K .EQ. 1 .AND. L .EQ. 1		&		!! It is the first step
								.AND. DissolutionData(I,11) .EQ. 1	&		!! It is condensation of an electrolyte w/out dissociation
								.AND. DissolutionData(I,12) .GT. 1) THEN    !! There is a related dissociation reaction identified
								
								CALL ReformElectrolyte(INT(DissolutionData(I,12)), Current)
								
							END IF

							!If water content is too low, skip dissolution reactions
							!IF (((ThermoBulkMode .AND. Current%AqChems(1) .LT. 1.0e-17) & !Dry Bulk
							!	.OR. ( NOT(ThermoBulkMode) .AND. Current%AqChems(1) .LT. 0.0001*TotalMoles)) .AND. & !Dry Particle (mole fraction of H2O < 0.01%)
							!	DissolutionData(I,11) .NE. 6.) CYCLE    !Not a gas-solid rxn

							IF (DissolutionData(I,11) .EQ. 6.) THEN
								CYCLE !Already did gas-solid (Type 6) reactions
							
							ELSE
								!Equilibrate the reaction for one particle
							        !WRITE(*,*) "Before EquilibrateDissolutionReaction"
								WaterSatBurden = GetSatVapBurden () 
                                                                CALL EquilibrateDissolutionReaction (GasPhaseBurden, I, &
                                                                      Current, Temperature, ReturnType, &
                                                                      WaterSatBurden=WaterSatBurden, GasPhaseBurden2=0.0, OptErrorTolerance=AqThermoEquilibriumError/2.)							
								!WRITE(*,*) "After EquilibrateDissolutionReaction"
							END IF

							!! Make a more refined guess if still equilibrating late in the game
							!! (That is, force it to equilibrate much closer to 1 each time through)
							IF ((K .GT. 10 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
								IF (DissolutionData(I,11) .EQ. 6.) THEN
									CYCLE !Already equilibrated gas-solid reactions
								ELSE
									WaterSatBurden = GetSatVapBurden () 
                                                                        !WRITE(*,*) "Before EquilibrateDissolutionReaction2"
                                                                        CALL EquilibrateDissolutionReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, &
									WaterSatBurden=WaterSatBurden, GasPhaseBurden2=0.0, OptErrorTolerance=AqThermoEquilibriumError/10.)
                                                                        !WRITE(*,*) "After EquilibrateDissolutionReaction2"
								END IF
							END IF

						!! HMDR+1 is a dummy saying "do the aqueous rxns."
						ELSE IF (I .EQ. HowManyDissolutionReactions+1) THEN
                                                        !Don't need to update activities - that was done 
                                                        !after previous dissolution reactions
							CALL FindElectrolyteEquilibrium(Current, .FALSE., &
                                                           FirstEquilibration=.FALSE., ReturnType=ReturnType)

							!If Electrolyte equilibrium hasn't been reached
							!after 50 times through particle list or 5 times through a single reaction
							!, just give up on it.
							IF (L .GT. 5 .OR. K .GT. 50) THEN							
								ReturnType = 2
							END IF

						END IF

						IF (ReturnType .NE. 2) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF
                                                !WRITE(*,*) "Before KusikMeissner"
						!! If not the aqueous electrolyte equilibration, then update the thermo
						IF (I .LE. HowManyDissolutionReactions) CALL KusikMeissner (Current)
                                                !WRITE(*,*) "After KusikMeissner"

						L = L+1
						Current => Current%Next
					END DO !Particles

				        !WRITE(*,*) "Number Eq steps for Reaction # ", I, " = ", L

					IF (L .GT. AqThermoNumbAllowedIterations) THEN
						!Just exit after max iterations
						RxnEquilibrated = .TRUE.
					END IF
				END DO !RxnEq

				!Matt Changed this!
				IF (I .LE. HowManyDissolutionReactions) THEN

					!! Replace the gas phase chemical concentrations
					IF (INT(DissolutionData(I,1)) .LE. HowManyEvolveGasChems) THEN
						CALL ReplaceGridCellChemBurden ( INT(DissolutionData(I,1)), GasPhaseBurden)
					END IF

					IF (INT(DissolutionData(I,2)) .LE. HowManyEvolveGasChems .AND. DissolutionData(I,11) .EQ. 6.)	&
						CALL ReplaceGridCellChemBurden ( INT(DissolutionData(I,2)), GasPhaseBurden2)
				END IF	 
			END IF
		
		END DO !Rxn #
			
		!Just exit after 500 iterations
		IF (K .GT. 500) THEN
			Equilibrated = .TRUE.
		END IF

	END DO !NotinEq
	
	RETURN
END SUBROUTINE EquilibrateInorgGridPoint

logical function isnan(a)
   IMPLICIT NONE
   real ::  a
      if (a.ne.a) then
         isnan = .true.
      else
         isnan = .false.
      end if
   return
end 
