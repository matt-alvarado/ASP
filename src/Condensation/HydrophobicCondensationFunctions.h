!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! HydrophobicCondensationFunctions.h
!! Contains many functions needed by the numerical integration routines
!! for organic dissolution in the aqueous and hydrophobic phases.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						             !!
!!								   	     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History		             !!
!! 08/28  2006   Matt Alvarado     1. Changed zero flag for	             !!
!!				   EquilibrateOrganicReaction		     !!
!!				   2. Set a minimum value for		     !!
!!				      TotalMolesOM (10^-30)		     !!
!! 08/29  2006   Matt Alvarado	   Reset minimum value for		     !!
!!				      TotalMolesOM (10^-25)		
!! 09/08  2006   Matt Alvarado     Created BinaryWaterActivityOrganic
!! 07/10  2007   Matt Alvarado     Rewrote CalculateSurrogateGasConc
!!				   to set all negative values to 0.
!! 07/17  2007   Matt Alvarado     Edited TotalMolesOM to apply
!!				  minimum flag correctly.
!! 09/10  2007   Matt Alvarado     Edited EquilibrateOrganicParticle
!!				  to just exit after 500 iterations
!! 05/26  2010   Matt Alvarado    Just use unity activity coefficients
!! 10/15  2010   Matt Alvarado    Remove Optional Arguments
!! 02/17  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 01/02  2013   Matt Alvarado    Replaced UR21 with CH3COCO2H               !!
!! 01/23  2013   Matt Alvarado    Added flag to SetOrganicDissolution to     !!
!!                                  decide between two vaopr pressure methods!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		     !!
!! 1. FUNCTION SaturationVaporPressure (Temperature, OrganicDissolutionRxn)  !!
!! 2. FUNCTION TotalMolesOM(InParticle)                                      !!
!! 3. SUBROUTINE UpdateHydrophobicUNIFAC(InParticle, Temperature)
!! 4. SUBROUTINE CalculateSurrogateGasConc(lots)
!! 5. SUBROUTINE UpdateGasPhaseOrganics(lots)
!! 6. SUBROUTINE EquilibrateHydrophobicPhaseAtGridPoint()
!! 7. SUBROUTINE EquilibrateHydrophobicReaction (GasPhaseBurden, WhichRxn, &
!!               InParticle, Temperature, ReturnType, OptErrorTolerance)
!! 8. FUNCTION HydrophobicEquilibriumConstantsRatio (GasPhaseBurden, & 
!!             WhichRxn, Temperature, InParticle)
!! 9. SUBROUTINE EquilibrateHydrophilicReaction
!! 10. FUNCTION HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, & 
!!              WhichRxn, Temperature, InParticle)		
!! 11. FUNCTION EffectiveOrganicHenrysLaw (Temperature, ProtonMolality, & 
!!              AqOrganicDissolutionRxn)
!! 12. SUBROUTINE UpdateHydrophilicUNIFAC(InParticle, Temperature)
!! 13. FUNCTION TotalAqMoles (InParticle)
!! 14. SUBROUTINE CalcInfiniteDilutionHydrophilicUNIFAC
!! 15. SUBROUTINE EquilibrateAllOrganicsAtGridPoint
!! 16. SUBROUTINE EquilibrateOrganicReaction
!! 17. FUNCTION TwoPhaseOrgEqRatio(WhichAqRxn, WhichOrgRxn, InParticle, & 
!!                 Temperature)
!! 18. SUBROUTINE EquilibrateOrganicParticle
!! 19. FUNCTION BinaryWaterActivityOrganic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate SaturationVaporPressure !!
!! This is indexed by the organic    !!
!! dissolution reaction. 	     !!
!! Puts out in MBar, so make         !!
!! sure to convert the weird         !!
!! pressure units.		     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION SaturationVaporPressure  (Temperature, OrganicDissolutionRxn)

	USE ModelParameters, ONLY : EqCoeffTref 

	IMPLICIT NONE

	REAL*8  :: Temperature, TempRatio
	INTEGER :: OrganicDissolutionRxn


        IF(OrganicDissolutionData(OrganicDissolutionRxn,2) .LT. 0.005) THEN
	    !Boiling point divided by Temperature
	    TempRatio = OrganicDissolutionData(OrganicDissolutionRxn,5)/Temperature

	    !Units of mbar, see Myrdal and Yalkowski and Schwarzenbach et al.,2003
	    !"Environmental Organic Chemistry"
            SaturationVaporPressure  = 1000*EXP ( -1.0*(21.2 + 0.3*OrganicDissolutionData(OrganicDissolutionRxn,6) &
	                                 + 177.0*OrganicDissolutionData(OrganicDissolutionRxn,7) ) &
					 *(TempRatio - 1.0) + (10.8 + 0.25*OrganicDissolutionData(OrganicDissolutionRxn,6)) &
					 *LOG(TempRatio) ) / OrganicDissolutionData(OrganicDissolutionRxn,8)
        ELSE
        
            SaturationVaporPressure = OrganicDissolutionData(OrganicDissolutionRxn,14)*(300.0/Temperature) & 
                                  *exp(OrganicDissolutionData(OrganicDissolutionRxn,15)*(1000.0/8.314)*(1.0/300.0-1.0/Temperature))
        END IF   

	RETURN
END FUNCTION SaturationVaporPressure 

!Moles per particle of hydrophobic organics
!Note that this value is forced to be above 1.0e-25 to
!keep hydrophobic equilibrium constants from dividing by ~0, and
!giving an error in LSODES
REAL*8 FUNCTION TotalMolesOM(InParticle) 

	USE Chemistry,				ONLY :	HowManyOrgChems
	
	IMPLICIT NONE

	REAL*8  :: Total
	INTEGER :: I
	TYPE(PARTICLE),POINTER :: InParticle

	Total = 0.0
	!Skip I = 1 since that is hard-wired for BC
	DO I = 2, HowManyOrgChems
		Total = Total + InParticle%OrgChems(I)
	END DO

	!Force total moles of OM to be above 1.0e-25 moles/particle
	!(prevents zeros in equilibrium coefficients)
	IF(Total .GE. 1.0E-25) THEN
		TotalMolesOM = Total
	ELSE
		TotalMolesOM = 1.0E-25
	END IF
	

	RETURN

END FUNCTION TotalMolesOM

!! This subroutine provides the link between the main model and the 
!! UNIFAC hydrophobic driver file soa_unidriv.f
!! Note that the order of the hydrophobic organic compounds in 
!! unifacparam.h and in OrgPhaseChems.in must be the same! 
SUBROUTINE UpdateHydrophobicUNIFAC(InParticle, Temperature)

	USE Chemistry,				ONLY :	HowManyOrgChems

	IMPLICIT NONE

	INTEGER :: I, N
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8, ALLOCATABLE :: XPASS(:), GAMMA(:)
	REAL*8 :: TotalOM, Temperature
	

	IF (1 .EQ. 0) THEN !UNIFAC
		N = 23 !Number of compounds
		ALLOCATE(XPASS(0:N-1), GAMMA(0:N-1))
	
		!Calculate mole fractions, putting them in an array ordered from 0:N-1
		TotalOM = TotalMolesOM(InParticle)
		
		!If no organic is present, set to 1.0
		IF(TotalOM .LE. 0.0) THEN
			DO I = 0, N-1
				InParticle%HydrophobicActivityCoeffs(I+2) = 1.0
			END DO
			RETURN
		END IF

		DO I = 0, N-1
			!Note the I+2 since OrgChems is ordered from 1 to 24 with BC as the first compound
			XPASS(I) = InParticle%OrgChems(I+2)/TotalOM
		END DO

		!Call UNIFAC, placing activity coefficients in GAMMA(0:N-1)
		CALL UNIDRIV(XPASS,GAMMA,N,Temperature)
	
		!Place activity coefficients into particle array
		DO I = 0, N-1
			InParticle%HydrophobicActivityCoeffs(I+2) = GAMMA(I)
		END DO
	
	ELSE !Assume activity coefficients equal unity

		N = HowManyOrgChems-1

		DO I = 0, N-1
			InParticle%HydrophobicActivityCoeffs(I+2) = 1.0
		END DO
		
	END IF

	RETURN
END SUBROUTINE UpdateHydrophobicUNIFAC


!! This subroutine takes the gas-phase concentrations of condensible
!! organic coumpounds and orders them into 10 surrogate compounds
!! for gas-particle partitioning calculations. The mixing ratio (ppbv)
!! of the surrogate compounds is calculated to conserve system mass.
!! That is, the mixing ratio of each compound is converted to a mass 
!! concentration, which is summed for each surrogate compound. 
!! The ppbv of the surrogate is then calculated using the MW of the surrogate. 
SUBROUTINE CalculateSurrogateGasConc(SOA1MassFracs, SOA2MassFracs, &
           SOA3MassFracs, SOA4MassFracs, SOA5MassFracs, SOA6MassFracs, &
           SOA7MassFracs, SOA8MassFracs, SOA9MassFracs, &
           SO10MassFracs, InitialBurden)

	USE Chemistry,		ONLY : FindChem, GasMolecularMass

	USE GridPointFields,	ONLY : GridGasChem, &
					GetGasChem, &
					ReplaceGridCellChemBurden

	USE ModelParameters,     ONLY : Avogadro 

	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	INTEGER :: UR21i, UR28i, SOA1i, &
			   RP13i, RP17i, RP18i, UR26i, UR29i, UR30i, SOA2i, &
			   RPR9i, RP12i, SOA3i, &
			   UR3i, UR8i, UR23i, SOA4i, &
			   UR7i, UR17i, SOA5i, &
			   AP1i, AP6i, UR31i, SOA6i, &
			   ADACi, ARACi, RPR7i, RP14i, RP19i, UR2i, &
                           UR14i, UR27i, SOA7i, &
			   AP10i, UR11i, UR15i, SOA8i, &
			   AP11i, AP12i, UR20i, UR34i, SOA9i, &
			   AP7i, AP8i, UR5i, UR6i, SO10i, I

	REAL*8  :: SOA1Temp, SOA2Temp, SOA3Temp, SOA4Temp, SOA5Temp, &
			SOA6Temp, SOA7Temp, SOA8Temp, SOA9Temp, SO10Temp, &
			SOA1MassFracs(2), SOA2MassFracs(6), SOA3MassFracs(2), &
			SOA4MassFracs(3), SOA5MassFracs(2), SOA6MassFracs(3), &
			SOA7MassFracs(8), SOA8MassFracs(3), SOA9MassFracs(4), &
			SO10MassFracs(4), InitialBurden(10),  Sum
	
	UR21i = FindChem("CH3COCO2H", 0)
	UR28i = FindChem("UR28", 0)
	SOA1i = FindChem("SOA1", 0)

	!Mass Conc of surrogate, ug/m3
	SOA1Temp = (1e12/Avogadro)*(GridGasChem(UR21i)*GasMolecularMass(UR21i) + &
				GridGasChem(UR28i)*GasMolecularMass(UR28i))
		
	!If the surrogate concentration = 0, tell model to evenly
	!divide up mass fractions, else calculate based on gas-phase concentrations
	IF (SOA1Temp .LE. 0.0) THEN
		SOA1Temp = 0.
		SOA1MassFracs(1) = 0.5
		SOA1MassFracs(2) = 0.5
	ELSE
		SOA1MassFracs(1) = (1e12/Avogadro)*(GridGasChem(UR21i)*GasMolecularMass(UR21i)) &
						/SOA1Temp
		SOA1MassFracs(2) = (1e12/Avogadro)*(GridGasChem(UR28i)*GasMolecularMass(UR28i)) &
						/SOA1Temp
	END IF

	RP13i = FindChem("RP13", 0)
	RP17i = FindChem("RP17", 0)
	RP18i = FindChem("RP18", 0)
	UR26i = FindChem("UR26", 0)
	UR29i = FindChem("UR29", 0)
	UR30i = FindChem("UR30", 0)
	SOA2i = FindChem("SOA2", 0)

	!Mass Conc of surrogate, ug/m3
	SOA2Temp = (1e12/Avogadro)*(GridGasChem(RP13i)*GasMolecularMass(RP13i) + &
			   GridGasChem(RP17i)*GasMolecularMass(RP17i) + &
			   GridGasChem(RP18i)*GasMolecularMass(RP18i) + &
			   GridGasChem(UR26i)*GasMolecularMass(UR26i) + &
			   GridGasChem(UR29i)*GasMolecularMass(UR29i) + &
				GridGasChem(UR30i)*GasMolecularMass(UR30i))

	IF (SOA2Temp .LE. 0.0) THEN
		SOA2Temp = 0.
		SOA2MassFracs(1) = 1.0/6.0
		SOA2MassFracs(2) = 1.0/6.0
		SOA2MassFracs(3) = 1.0/6.0
		SOA2MassFracs(4) = 1.0/6.0
		SOA2MassFracs(5) = 1.0/6.0
		SOA2MassFracs(6) = 1.0/6.0
	ELSE
		SOA2MassFracs(1) = (1e12/Avogadro)*(GridGasChem(RP13i)*GasMolecularMass(RP13i)) &
						/SOA2Temp
		SOA2MassFracs(2) = (1e12/Avogadro)*(GridGasChem(RP17i)*GasMolecularMass(RP17i)) &
						/SOA2Temp
		SOA2MassFracs(3) = (1e12/Avogadro)*(GridGasChem(RP18i)*GasMolecularMass(RP18i)) &
						/SOA2Temp
		SOA2MassFracs(4) = (1e12/Avogadro)*(GridGasChem(UR26i)*GasMolecularMass(UR26i)) &
						/SOA2Temp
		SOA2MassFracs(5) = (1e12/Avogadro)*(GridGasChem(UR29i)*GasMolecularMass(UR29i)) &
						/SOA2Temp
		SOA2MassFracs(6) = (1e12/Avogadro)*(GridGasChem(UR30i)*GasMolecularMass(UR30i)) &
						/SOA2Temp
	END IF
		
	RPR9i = FindChem("RPR9", 0)
	RP12i = FindChem("RP12", 0)
	SOA3i = FindChem("SOA3", 0)

	!Mass Conc of surrogate, ug/m3
	SOA3Temp = (1e12/Avogadro)*(GridGasChem(RPR9i)*GasMolecularMass(RPR9i) + &
			   GridGasChem(RP12i)*GasMolecularMass(RP12i))
	
	IF (SOA3Temp .LE. 0.0) THEN
		SOA3Temp = 0.
		SOA3MassFracs(1) = 0.5
		SOA3MassFracs(2) = 0.5
	ELSE
		SOA3MassFracs(1) = (1e12/Avogadro)*(GridGasChem(RPR9i)*GasMolecularMass(RPR9i)) &
						/SOA3Temp
		SOA3MassFracs(2) = (1e12/Avogadro)*(GridGasChem(RP12i)*GasMolecularMass(RP12i)) &
						/SOA3Temp
	END IF

	UR3i = FindChem("UR3", 0)
	UR8i = FindChem("UR8", 0)
	UR23i = FindChem("UR23", 0)
	SOA4i = FindChem("SOA4", 0)

	!Mass Conc of surrogate, ug/m3
	SOA4Temp = (1e12/Avogadro)*(GridGasChem(UR3i)*GasMolecularMass(UR3i) + &
			   GridGasChem(UR8i)*GasMolecularMass(UR8i) + &
			   GridGasChem(UR23i)*GasMolecularMass(UR23i))

	IF (SOA4Temp .LE. 0.0) THEN
		SOA4Temp = 0.
		SOA4MassFracs(1) = 1.0/3.0
		SOA4MassFracs(2) = 1.0/3.0
		SOA4MassFracs(3) = 1.0/3.0
	ELSE
		SOA4MassFracs(1) = (1e12/Avogadro)*(GridGasChem(UR3i)*GasMolecularMass(UR3i)) &
						/SOA4Temp
		SOA4MassFracs(2) = (1e12/Avogadro)*(GridGasChem(UR8i)*GasMolecularMass(UR8i)) &
						/SOA4Temp
		SOA4MassFracs(3) = (1e12/Avogadro)*(GridGasChem(UR23i)*GasMolecularMass(UR23i)) &
						/SOA4Temp
	END IF

	UR7i = FindChem("UR7", 0)
	UR17i = FindChem("UR17", 0)
	SOA5i = FindChem("SOA5", 0)

	!Mass Conc of surrogate, ug/m3
	SOA5Temp = (1e12/Avogadro)*(GridGasChem(UR7i)*GasMolecularMass(UR7i) + &
			   GridGasChem(UR17i)*GasMolecularMass(UR17i))

	IF (SOA5Temp .LE. 0.0) THEN
		SOA5Temp = 0.
		SOA5MassFracs(1) = 0.5
		SOA5MassFracs(2) = 0.5
	ELSE
		SOA5MassFracs(1) = (1e12/Avogadro)*(GridGasChem(UR7i)*GasMolecularMass(UR7i)) &
						/SOA5Temp
		SOA5MassFracs(2) = (1e12/Avogadro)*(GridGasChem(UR17i)*GasMolecularMass(UR17i)) &
						/SOA5Temp
	END IF

	AP1i = FindChem("AP1", 0)
	AP6i = FindChem("AP6", 0)
	UR31i = FindChem("UR31", 0)
	SOA6i = FindChem("SOA6", 0)

	!Mass Conc of surrogate, ug/m3
	SOA6Temp = (1e12/Avogadro)*(GridGasChem(AP1i)*GasMolecularMass(AP1i) + &
			   GridGasChem(AP6i)*GasMolecularMass(AP6i) + &
			   GridGasChem(UR31i)*GasMolecularMass(UR31i))

	IF (SOA6Temp .LE. 0.0) THEN
		SOA6Temp = 0.
		SOA6MassFracs(1) = 1.0/3.0
		SOA6MassFracs(2) = 1.0/3.0
		SOA6MassFracs(3) = 1.0/3.0
	ELSE
		SOA6MassFracs(1) = (1e12/Avogadro)*(GridGasChem(AP1i)*GasMolecularMass(AP1i)) &
						/SOA6Temp
		SOA6MassFracs(2) = (1e12/Avogadro)*(GridGasChem(AP6i)*GasMolecularMass(AP6i)) &
						/SOA6Temp
		SOA6MassFracs(3) = (1e12/Avogadro)*(GridGasChem(UR31i)*GasMolecularMass(UR31i)) &
						/SOA6Temp
	END IF

	ADACi = FindChem("ADAC", 0)
	ARACi = FindChem("ARAC", 0)
	RPR7i = FindChem("RPR7", 0)
	RP14i = FindChem("RP14", 0)
	RP19i = FindChem("RP19", 0)
	UR2i  = FindChem("UR2", 0)
	UR14i = FindChem("UR14", 0)
	UR27i = FindChem("UR27", 0)
	SOA7i = FindChem("SOA7", 0)

	!Mass Conc of surrogate, ug/m3
	SOA7Temp = (1e12/Avogadro)*(GridGasChem(ADACi)*GasMolecularMass(ADACi) + &
			   GridGasChem(ARACi)*GasMolecularMass(ARACi) + &
			   GridGasChem(RPR7i)*GasMolecularMass(RPR7i) + &
			   GridGasChem(RP14i)*GasMolecularMass(RP14i) + &
			   GridGasChem(RP19i)*GasMolecularMass(RP19i) + &
				GridGasChem(UR2i)*GasMolecularMass(UR2i) + &
				GridGasChem(UR14i)*GasMolecularMass(UR14i) + &
			GridGasChem(UR27i)*GasMolecularMass(UR27i))

	IF (SOA7Temp .LE. 0.0) THEN
		SOA7Temp = 0.
		SOA7MassFracs(1) = 0.125
		SOA7MassFracs(2) = 0.125
		SOA7MassFracs(3) = 0.125
		SOA7MassFracs(4) = 0.125
		SOA7MassFracs(5) = 0.125
		SOA7MassFracs(6) = 0.125
		SOA7MassFracs(7) = 0.125
		SOA7MassFracs(8) = 0.125
	ELSE
		SOA7MassFracs(1) = (1e12/Avogadro)*(GridGasChem(ADACi)*GasMolecularMass(ADACi)) &
						/SOA7Temp
		SOA7MassFracs(2) = (1e12/Avogadro)*(GridGasChem(ARACi)*GasMolecularMass(ARACi)) &
						/SOA7Temp
		SOA7MassFracs(3) = (1e12/Avogadro)*(GridGasChem(RPR7i)*GasMolecularMass(RPR7i)) &
						/SOA7Temp
		SOA7MassFracs(4) = (1e12/Avogadro)*(GridGasChem(RP14i)*GasMolecularMass(RP14i)) &
						/SOA7Temp
		SOA7MassFracs(5) = (1e12/Avogadro)*(GridGasChem(RP19i)*GasMolecularMass(RP19i)) &
						/SOA7Temp
		SOA7MassFracs(6) = (1e12/Avogadro)*(GridGasChem(UR2i)*GasMolecularMass(UR2i)) &
						/SOA7Temp
		SOA7MassFracs(7) = (1e12/Avogadro)*(GridGasChem(UR14i)*GasMolecularMass(UR14i)) &
						/SOA7Temp
		SOA7MassFracs(8) = (1e12/Avogadro)*(GridGasChem(UR27i)*GasMolecularMass(UR27i)) &
						/SOA7Temp
	
		Sum = 0.0
		DO I = 1,8
			Sum = SOA7MassFracs(I) + Sum
		END DO

		IF (Sum .GT. 1.000001 .OR. Sum .LT. 0.9999999) THEN
			WRITE(*,*) Sum
			CALL ERROR("SOA7 Mass fractions do not sum to 1!")
		END IF

	END IF
	
	AP10i = FindChem("AP10", 0)
	UR11i = FindChem("UR11", 0)
	UR15i = FindChem("UR15", 0)
	SOA8i = FindChem("SOA8", 0)

	!Mass Conc of surrogate, ug/m3
	SOA8Temp = (1e12/Avogadro)*(GridGasChem(AP10i)*GasMolecularMass(AP10i) + &
			   GridGasChem(UR11i)*GasMolecularMass(UR11i) + &
			   GridGasChem(UR15i)*GasMolecularMass(UR15i))

	IF (SOA8Temp .LE. 0.0) THEN
		SOA8Temp = 0.
		SOA8MassFracs(1) = 1.0/3.0
		SOA8MassFracs(2) = 1.0/3.0
		SOA8MassFracs(3) = 1.0/3.0
	ELSE
		SOA8MassFracs(1) = (1e12/Avogadro)*(GridGasChem(AP10i)*GasMolecularMass(AP10i)) &
						/SOA8Temp
		SOA8MassFracs(2) = (1e12/Avogadro)*(GridGasChem(UR11i)*GasMolecularMass(UR11i)) &
						/SOA8Temp
		SOA8MassFracs(3) = (1e12/Avogadro)*(GridGasChem(UR15i)*GasMolecularMass(UR15i)) &
						/SOA8Temp
	END IF

	AP11i = FindChem("AP11", 0)
	AP12i = FindChem("AP12", 0)
	UR20i = FindChem("UR20", 0)
	UR34i = FindChem("UR34", 0)
	SOA9i = FindChem("SOA9", 0)

	!Mass Conc of surrogate, ug/m3
	SOA9Temp = (1e12/Avogadro)*(GridGasChem(AP11i)*GasMolecularMass(AP11i) + &
			   GridGasChem(AP12i)*GasMolecularMass(AP12i) + &
			   GridGasChem(UR20i)*GasMolecularMass(UR20i) + &
			   GridGasChem(UR34i)*GasMolecularMass(UR34i))

	IF (SOA9Temp .LE. 0.0) THEN
		SOA9Temp = 0.
		SOA9MassFracs(1) = 0.25
		SOA9MassFracs(2) = 0.25
		SOA9MassFracs(3) = 0.25
		SOA9MassFracs(4) = 0.25
	ELSE
		SOA9MassFracs(1) = (1e12/Avogadro)*(GridGasChem(AP11i)*GasMolecularMass(AP11i)) &
						/SOA9Temp
		SOA9MassFracs(2) = (1e12/Avogadro)*(GridGasChem(AP12i)*GasMolecularMass(AP12i)) &
						/SOA9Temp
		SOA9MassFracs(3) = (1e12/Avogadro)*(GridGasChem(UR20i)*GasMolecularMass(UR20i)) &
						/SOA9Temp
		SOA9MassFracs(4) = (1e12/Avogadro)*(GridGasChem(UR34i)*GasMolecularMass(UR34i)) &
						/SOA9Temp
	END IF
		
	AP7i = FindChem("AP7", 0)
	AP8i = FindChem("AP8", 0)
	UR5i = FindChem("UR5", 0)
	UR6i = FindChem("UR6", 0)
	SO10i = FindChem("SO10", 0)

	!Mass Conc of surrogate, ug/m3
	SO10Temp = (1e12/Avogadro)*(GridGasChem(AP7i)*GasMolecularMass(AP7i) + &
			   GridGasChem(AP8i)*GasMolecularMass(AP8i) + &
			   GridGasChem(UR5i)*GasMolecularMass(UR5i) + &
			   GridGasChem(UR6i)*GasMolecularMass(UR6i))

	IF (SO10Temp .LE. 0.0) THEN
		SO10Temp = 0.
		SO10MassFracs(1) = 0.25
		SO10MassFracs(2) = 0.25
		SO10MassFracs(3) = 0.25
		SO10MassFracs(4) = 0.25
	ELSE
		SO10MassFracs(1) = (1e12/Avogadro)*(GridGasChem(AP7i)*GasMolecularMass(AP7i)) &
						/SO10Temp
		SO10MassFracs(2) = (1e12/Avogadro)*(GridGasChem(AP8i)*GasMolecularMass(AP8i)) &
						/SO10Temp
		SO10MassFracs(3) = (1e12/Avogadro)*(GridGasChem(UR5i)*GasMolecularMass(UR5i)) &
						/SO10Temp
		SO10MassFracs(4) = (1e12/Avogadro)*(GridGasChem(UR6i)*GasMolecularMass(UR6i)) &
						/SO10Temp
	END IF

	!Convert concentrations to molecules/cm3
	SOA1Temp = (Avogadro/1e12)*(SOA1Temp/GasMolecularMass(SOA1i))
	SOA2Temp = (Avogadro/1e12)*(SOA2Temp/GasMolecularMass(SOA2i))
	SOA3Temp = (Avogadro/1e12)*(SOA3Temp/GasMolecularMass(SOA3i))
	SOA4Temp = (Avogadro/1e12)*(SOA4Temp/GasMolecularMass(SOA4i))
	SOA5Temp = (Avogadro/1e12)*(SOA5Temp/GasMolecularMass(SOA5i))
	SOA6Temp = (Avogadro/1e12)*(SOA6Temp/GasMolecularMass(SOA6i))
	SOA7Temp = (Avogadro/1e12)*(SOA7Temp/GasMolecularMass(SOA7i))
	SOA8Temp = (Avogadro/1e12)*(SOA8Temp/GasMolecularMass(SOA8i))
	SOA9Temp = (Avogadro/1e12)*(SOA9Temp/GasMolecularMass(SOA9i))
	SO10Temp = (Avogadro/1e12)*(SO10Temp/GasMolecularMass(SO10i))

	!Return initial gas-phase burden (in molecules/cm3)
	InitialBurden(1) = SOA1Temp
	InitialBurden(2) = SOA2Temp
	InitialBurden(3) = SOA3Temp
	InitialBurden(4) = SOA4Temp
	InitialBurden(5) = SOA5Temp
	InitialBurden(6) = SOA6Temp
	InitialBurden(7) = SOA7Temp
	InitialBurden(8) = SOA8Temp
	InitialBurden(9) = SOA9Temp
	InitialBurden(10) = SO10Temp

	!Replace burdens (updates gas-phase concentrations of surogates)
	CALL ReplaceGridCellChemBurden ( SOA1i, SOA1Temp)
	CALL ReplaceGridCellChemBurden ( SOA2i, SOA2Temp)
	CALL ReplaceGridCellChemBurden ( SOA3i, SOA3Temp)
	CALL ReplaceGridCellChemBurden ( SOA4i, SOA4Temp)
	CALL ReplaceGridCellChemBurden ( SOA5i, SOA5Temp)
	CALL ReplaceGridCellChemBurden ( SOA6i, SOA6Temp)
	CALL ReplaceGridCellChemBurden ( SOA7i, SOA7Temp)
	CALL ReplaceGridCellChemBurden ( SOA8i, SOA8Temp)
	CALL ReplaceGridCellChemBurden ( SOA9i, SOA9Temp)
	CALL ReplaceGridCellChemBurden ( SO10i, SO10Temp)

END SUBROUTINE CalculateSurrogateGasConc 

SUBROUTINE UpdateGasPhaseOrganics(S1MassFracs, S2MassFracs, &
                S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
		S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs,& 
                InitialBurden)

	USE Chemistry,		 ONLY : FindChem, GasMolecularMass

	USE GridPointFields,	 ONLY : GridGasChem, &
					ReplaceGridCellChemBurden, &
					GetGridCellChemBurden

	USE ModelParameters,     ONLY : Avogadro 

	IMPLICIT NONE

	INTEGER ::   I, UR21i, UR28i, SOA1i, &
	             RP13i, RP17i, RP18i, UR26i, UR29i, UR30i, SOA2i, &
		     RPR9i, RP12i, SOA3i, &
		     UR3i, UR8i, UR23i, SOA4i, &
			   UR7i, UR17i, SOA5i, &
			   AP1i, AP6i, UR31i, SOA6i, &
			   ADACi, ARACi, RPR7i, RP14i, &
                           RP19i, UR2i, UR14i, UR27i, SOA7i, &
			   AP10i, UR11i, UR15i, SOA8i, &
			   AP11i, AP12i, UR20i, UR34i, SOA9i, &
			   AP7i, AP8i, UR5i, UR6i, SO10i, HowManyNonLumpedRxns

	REAL*8  :: 	S1MassFracs(2), S2MassFracs(6), S3MassFracs(2), &
			S4MassFracs(3), S5MassFracs(2), S6MassFracs(3), &
			S7MassFracs(8), S8MassFracs(3), S9MassFracs(4), &
			S10MassFracs(4), InitialBurden(10), DeltaBurden(10), &
		DeltaMassConc(10), S1DeltaMassConc(2), S2DeltaMassConc(6), &
		S3DeltaMassConc(2), S4DeltaMassConc(3), S5DeltaMassConc(2), &
		S6DeltaMassConc(3), S7DeltaMassConc(8), S8DeltaMassConc(3), &
				S9DeltaMassConc(4), S10DeltaMassConc(4), &
				UR21DeltaBurden, UR28DeltaBurden, &
			   RP13DeltaBurden, RP17DeltaBurden, RP18DeltaBurden, &
                          UR26DeltaBurden, UR29DeltaBurden, UR30DeltaBurden, &
			   RPR9DeltaBurden, RP12DeltaBurden,  &
			   UR3DeltaBurden, UR8DeltaBurden, UR23DeltaBurden,  &
			   UR7DeltaBurden, UR17DeltaBurden,  &
			   AP1DeltaBurden, AP6DeltaBurden, UR31DeltaBurden,  &
			   ADACDeltaBurden, ARACDeltaBurden, RPR7DeltaBurden, &
                           RP14DeltaBurden, RP19DeltaBurden, UR2DeltaBurden, &
                           UR14DeltaBurden, UR27DeltaBurden, &
	AP10DeltaBurden, UR11DeltaBurden, UR15DeltaBurden, UR19DeltaBurden,  &
	AP11DeltaBurden, AP12DeltaBurden, UR20DeltaBurden, UR34DeltaBurden,  &
	AP7DeltaBurden, AP8DeltaBurden, UR5DeltaBurden, UR6DeltaBurden, &
			   UR21NewBurden, UR28NewBurden, &
	RP13NewBurden, RP17NewBurden, RP18NewBurden, UR26NewBurden, & 
                          UR29NewBurden, UR30NewBurden, &
			   RPR9NewBurden, RP12NewBurden,  &
			   UR3NewBurden, UR8NewBurden, UR23NewBurden,  &
			   UR7NewBurden, UR17NewBurden,  &
			   AP1NewBurden, AP6NewBurden, UR31NewBurden,  &
			   ADACNewBurden, ARACNewBurden, RPR7NewBurden, &
                           RP14NewBurden, RP19NewBurden, UR2NewBurden, &
                           UR14NewBurden, UR27NewBurden, &
			   AP10NewBurden, UR11NewBurden, UR15NewBurden,  &
		AP11NewBurden, AP12NewBurden, UR20NewBurden, UR34NewBurden,  &
			AP7NewBurden, AP8NewBurden, UR5NewBurden, UR6NewBurden


	!This tells the lumping to skip the first 8 non-lumped reactions
	!You must adjust this if you change the reaction number or order!
	HowManyNonLumpedRxns = 8
	
	!Calculate change in the gas-phase burdens of the surrogate species
	DO I = 1, 10
		DeltaBurden(I) = InitialBurden(I) - GetGridCellChemBurden ( INT(OrganicDissolutionData(I+HowManyNonLumpedRxns,1)))
	END DO

	!Convert Burden to mass concentration (ug/m3)
	DO I = 1, 10
		DeltaMassConc(I) = (1.0e12/Avogadro)*(DeltaBurden(I)*GasMolecularMass(INT(OrganicDissolutionData(I+HowManyNonLumpedRxns,1))))
	END DO

	!Calculate the change in mass conc for the gas-phase species
	!based on their initial mass fraction contribution to the surrogate
	S1DeltaMassConc(1) = S1MassFracs(1)*DeltaMassConc(1)
	S1DeltaMassConc(2) = S1MassFracs(2)*DeltaMassConc(1)

	S2DeltaMassConc(1) = S2MassFracs(1)*DeltaMassConc(2)
	S2DeltaMassConc(2) = S2MassFracs(2)*DeltaMassConc(2)
	S2DeltaMassConc(3) = S2MassFracs(3)*DeltaMassConc(2)
	S2DeltaMassConc(4) = S2MassFracs(4)*DeltaMassConc(2)
	S2DeltaMassConc(5) = S2MassFracs(5)*DeltaMassConc(2)
	S2DeltaMassConc(6) = S2MassFracs(6)*DeltaMassConc(2)

	S3DeltaMassConc(1) = S3MassFracs(1)*DeltaMassConc(3)
	S3DeltaMassConc(2) = S3MassFracs(2)*DeltaMassConc(3)
	
	S4DeltaMassConc(1) = S4MassFracs(1)*DeltaMassConc(4)
	S4DeltaMassConc(2) = S4MassFracs(2)*DeltaMassConc(4)
	S4DeltaMassConc(3) = S4MassFracs(3)*DeltaMassConc(4)
	
	S5DeltaMassConc(1) = S5MassFracs(1)*DeltaMassConc(5)
	S5DeltaMassConc(2) = S5MassFracs(2)*DeltaMassConc(5)

	S6DeltaMassConc(1) = S6MassFracs(1)*DeltaMassConc(6)
	S6DeltaMassConc(2) = S6MassFracs(2)*DeltaMassConc(6)
	S6DeltaMassConc(3) = S6MassFracs(3)*DeltaMassConc(6)
	
	S7DeltaMassConc(1) = S7MassFracs(1)*DeltaMassConc(7)
	S7DeltaMassConc(2) = S7MassFracs(2)*DeltaMassConc(7)
	S7DeltaMassConc(3) = S7MassFracs(3)*DeltaMassConc(7)
	S7DeltaMassConc(4) = S7MassFracs(4)*DeltaMassConc(7)
	S7DeltaMassConc(5) = S7MassFracs(5)*DeltaMassConc(7)
	S7DeltaMassConc(6) = S7MassFracs(6)*DeltaMassConc(7)
	S7DeltaMassConc(7) = S7MassFracs(7)*DeltaMassConc(7)
	S7DeltaMassConc(8) = S7MassFracs(8)*DeltaMassConc(7)

	S8DeltaMassConc(1) = S8MassFracs(1)*DeltaMassConc(8)
	S8DeltaMassConc(2) = S8MassFracs(2)*DeltaMassConc(8)
	S8DeltaMassConc(3) = S8MassFracs(3)*DeltaMassConc(8)

	S9DeltaMassConc(1) = S9MassFracs(1)*DeltaMassConc(9)
	S9DeltaMassConc(2) = S9MassFracs(2)*DeltaMassConc(9)
	S9DeltaMassConc(3) = S9MassFracs(3)*DeltaMassConc(9)
	S9DeltaMassConc(4) = S9MassFracs(4)*DeltaMassConc(9)

	S10DeltaMassConc(1) = S10MassFracs(1)*DeltaMassConc(10)
	S10DeltaMassConc(2) = S10MassFracs(2)*DeltaMassConc(10)
	S10DeltaMassConc(3) = S10MassFracs(3)*DeltaMassConc(10)
	S10DeltaMassConc(4) = S10MassFracs(4)*DeltaMassConc(10)

	!Update gas-phase concentrations of primary organics
	
	!Surrogate 1
	UR21i = FindChem("CH3COCO2H", 0)
	UR28i = FindChem("UR28", 0)

	UR21DeltaBurden = (Avogadro/1.0e12)*(S1DeltaMassConc(1)/GasMolecularMass(UR21i))
	UR21NewBurden =  GetGridCellChemBurden( UR21i)-UR21DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR21i, UR21NewBurden)
	
	UR28DeltaBurden = (Avogadro/1.0e12)*(S1DeltaMassConc(2)/GasMolecularMass(UR28i))
	UR28NewBurden =  GetGridCellChemBurden( UR28i)-UR28DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR28i, UR28NewBurden)

	!Surrogate 2
	RP13i = FindChem("RP13", 0)
	RP17i = FindChem("RP17", 0)
	RP18i = FindChem("RP18", 0)
	UR26i = FindChem("UR26", 0)
	UR29i = FindChem("UR29", 0)
	UR30i = FindChem("UR30", 0)

	RP13DeltaBurden = (Avogadro/1.0e12)*(S2DeltaMassConc(1)/GasMolecularMass(RP13i))
	RP13NewBurden =  GetGridCellChemBurden( RP13i)-RP13DeltaBurden
	CALL ReplaceGridCellChemBurden ( RP13i, RP13NewBurden)

	RP17DeltaBurden = (Avogadro/1.0e12)*(S2DeltaMassConc(2)/GasMolecularMass(RP17i))
	RP17NewBurden =  GetGridCellChemBurden( RP17i)-RP17DeltaBurden
	CALL ReplaceGridCellChemBurden ( RP17i, RP17NewBurden)

	RP18DeltaBurden = (Avogadro/1.0e12)*(S2DeltaMassConc(3)/GasMolecularMass(RP18i))
	RP18NewBurden =  GetGridCellChemBurden( RP18i)-RP18DeltaBurden
	CALL ReplaceGridCellChemBurden ( RP18i, RP18NewBurden)
	
	UR26DeltaBurden = (Avogadro/1.0e12)*(S2DeltaMassConc(4)/GasMolecularMass(UR26i))
	UR26NewBurden =  GetGridCellChemBurden( UR26i)-UR26DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR26i, UR26NewBurden)

	UR29DeltaBurden = (Avogadro/1.0e12)*(S2DeltaMassConc(5)/GasMolecularMass(UR29i))
	UR29NewBurden =  GetGridCellChemBurden( UR29i)-UR29DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR29i, UR29NewBurden)

	UR30DeltaBurden = (Avogadro/1.0e12)*(S2DeltaMassConc(6)/GasMolecularMass(UR30i))
	UR30NewBurden =  GetGridCellChemBurden( UR30i)-UR30DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR30i, UR30NewBurden)


	!Surrogate 3
	RPR9i = FindChem("RPR9", 0)
	RP12i = FindChem("RP12", 0)

	RPR9DeltaBurden = (Avogadro/1.0e12)*(S3DeltaMassConc(1)/GasMolecularMass(RPR9i))
	RPR9NewBurden =  GetGridCellChemBurden( RPR9i)-RPR9DeltaBurden
	CALL ReplaceGridCellChemBurden ( RPR9i, RPR9NewBurden)

	RP12DeltaBurden = (Avogadro/1.0e12)*(S3DeltaMassConc(2)/GasMolecularMass(RP12i))
	RP12NewBurden =  GetGridCellChemBurden( RP12i)-RP12DeltaBurden
	CALL ReplaceGridCellChemBurden ( RP12i, RP12NewBurden)

	!Surrogate 4
	UR3i = FindChem("UR3", 0)
	UR8i = FindChem("UR8", 0)
	UR23i = FindChem("UR23", 0)

	UR3DeltaBurden = (Avogadro/1.0e12)*(S4DeltaMassConc(1)/GasMolecularMass(UR3i))
	UR3NewBurden =  GetGridCellChemBurden( UR3i)-UR3DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR3i, UR3NewBurden)

	UR8DeltaBurden = (Avogadro/1.0e12)*(S4DeltaMassConc(2)/GasMolecularMass(UR8i))
	UR8NewBurden =  GetGridCellChemBurden( UR8i)-UR8DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR8i, UR8NewBurden)

	UR23DeltaBurden = (Avogadro/1.0e12)*(S4DeltaMassConc(3)/GasMolecularMass(UR23i))
	UR23NewBurden =  GetGridCellChemBurden( UR23i)-UR23DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR23i, UR23NewBurden)

	!Surrogate 5
	UR7i = FindChem("UR7", 0)
	UR17i = FindChem("UR17", 0)

	UR7DeltaBurden = (Avogadro/1.0e12)*(S5DeltaMassConc(1)/GasMolecularMass(UR7i))
	UR7NewBurden =  GetGridCellChemBurden( UR7i)-UR7DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR7i, UR7NewBurden)
	
	UR17DeltaBurden = (Avogadro/1.0e12)*(S5DeltaMassConc(2)/GasMolecularMass(UR17i))
	UR17NewBurden =  GetGridCellChemBurden( UR17i)-UR17DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR17i, UR17NewBurden)
	
	!Surrogate 6
	AP1i = FindChem("AP1", 0)
	AP6i = FindChem("AP6", 0)
	UR31i = FindChem("UR31", 0)

	AP1DeltaBurden = (Avogadro/1.0e12)*(S6DeltaMassConc(1)/GasMolecularMass(AP1i))
	AP1NewBurden =  GetGridCellChemBurden( AP1i)-AP1DeltaBurden
	CALL ReplaceGridCellChemBurden ( AP1i, AP1NewBurden)
	
	AP6DeltaBurden = (Avogadro/1.0e12)*(S6DeltaMassConc(2)/GasMolecularMass(AP6i))
	AP6NewBurden =  GetGridCellChemBurden( AP6i)-AP6DeltaBurden
	CALL ReplaceGridCellChemBurden ( AP6i, AP6NewBurden)

	UR31DeltaBurden = (Avogadro/1.0e12)*(S6DeltaMassConc(3)/GasMolecularMass(UR31i))
	UR31NewBurden =  GetGridCellChemBurden( UR31i)-UR31DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR31i, UR31NewBurden)

	!Surrogate 7
	ADACi = FindChem("ADAC", 0)
	ARACi = FindChem("ARAC", 0)
	RPR7i = FindChem("RPR7", 0)
	RP14i = FindChem("RP14", 0)
	RP19i = FindChem("RP19", 0)
	UR2i  = FindChem("UR2", 0)
	UR14i = FindChem("UR14", 0)
	UR27i = FindChem("UR27", 0)

	ADACDeltaBurden = (Avogadro/1.0e12)*(S7DeltaMassConc(1)/GasMolecularMass(ADACi))
	ADACNewBurden =  GetGridCellChemBurden( ADACi)-ADACDeltaBurden
	CALL ReplaceGridCellChemBurden ( ADACi, ADACNewBurden)

	ARACDeltaBurden = (Avogadro/1.0e12)*(S7DeltaMassConc(2)/GasMolecularMass(ARACi))
	ARACNewBurden =  GetGridCellChemBurden( ARACi)-ARACDeltaBurden
	CALL ReplaceGridCellChemBurden ( ARACi, ARACNewBurden)

	RPR7DeltaBurden = (Avogadro/1.0e12)*(S7DeltaMassConc(3)/GasMolecularMass(RPR7i))
	RPR7NewBurden =  GetGridCellChemBurden( RPR7i)-RPR7DeltaBurden
	CALL ReplaceGridCellChemBurden ( RPR7i, RPR7NewBurden)

	RP14DeltaBurden = (Avogadro/1.0e12)*(S7DeltaMassConc(4)/GasMolecularMass(RP14i))
	RP14NewBurden =  GetGridCellChemBurden( RP14i)-RP14DeltaBurden
	CALL ReplaceGridCellChemBurden ( RP14i, RP14NewBurden)

	RP19DeltaBurden = (Avogadro/1.0e12)*(S7DeltaMassConc(5)/GasMolecularMass(RP19i))
	RP19NewBurden =  GetGridCellChemBurden( RP19i)-RP19DeltaBurden
	CALL ReplaceGridCellChemBurden ( RP19i, RP19NewBurden)

	UR2DeltaBurden = (Avogadro/1.0e12)*(S7DeltaMassConc(6)/GasMolecularMass(UR2i))
	UR2NewBurden =  GetGridCellChemBurden( UR2i)-UR2DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR2i, UR2NewBurden)

	UR14DeltaBurden = (Avogadro/1.0e12)*(S7DeltaMassConc(7)/GasMolecularMass(UR14i))
	UR14NewBurden =  GetGridCellChemBurden( UR14i)-UR14DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR14i, UR14NewBurden)

	UR27DeltaBurden = (Avogadro/1.0e12)*(S7DeltaMassConc(8)/GasMolecularMass(UR27i))
	UR27NewBurden =  GetGridCellChemBurden( UR27i)-UR27DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR27i, UR27NewBurden)

	!Surrogate 8
	AP10i = FindChem("AP10", 0)
	UR11i = FindChem("UR11", 0)
	UR15i = FindChem("UR15", 0)

	AP10DeltaBurden = (Avogadro/1.0e12)*(S8DeltaMassConc(1)/GasMolecularMass(AP10i))
	AP10NewBurden =  GetGridCellChemBurden( AP10i)-AP10DeltaBurden
	CALL ReplaceGridCellChemBurden ( AP10i, AP10NewBurden)

	UR11DeltaBurden = (Avogadro/1.0e12)*(S8DeltaMassConc(2)/GasMolecularMass(UR11i))
	UR11NewBurden =  GetGridCellChemBurden( UR11i)-UR11DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR11i, UR11NewBurden)

	UR15DeltaBurden = (Avogadro/1.0e12)*(S8DeltaMassConc(3)/GasMolecularMass(UR15i))
	UR15NewBurden =  GetGridCellChemBurden( UR15i)-UR15DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR15i, UR15NewBurden)

	!Surrogate 9
	AP11i = FindChem("AP11", 0)
	AP12i = FindChem("AP12", 0)
	UR20i = FindChem("UR20", 0)
	UR34i = FindChem("UR34", 0)

	AP11DeltaBurden = (Avogadro/1.0e12)*(S9DeltaMassConc(1)/GasMolecularMass(AP11i))
	AP11NewBurden =  GetGridCellChemBurden( AP11i)-AP11DeltaBurden
	CALL ReplaceGridCellChemBurden ( AP11i, AP11NewBurden)

	AP12DeltaBurden = (Avogadro/1.0e12)*(S9DeltaMassConc(2)/GasMolecularMass(AP12i))
	AP12NewBurden =  GetGridCellChemBurden( AP12i)-AP12DeltaBurden
	CALL ReplaceGridCellChemBurden ( AP12i, AP12NewBurden)

	UR20DeltaBurden = (Avogadro/1.0e12)*(S9DeltaMassConc(3)/GasMolecularMass(UR20i))
	UR20NewBurden =  GetGridCellChemBurden( UR20i)-UR20DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR20i, UR20NewBurden)

	UR34DeltaBurden = (Avogadro/1.0e12)*(S9DeltaMassConc(4)/GasMolecularMass(UR34i))
	UR34NewBurden =  GetGridCellChemBurden( UR34i)-UR34DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR34i, UR34NewBurden)

	!Surrogate 10
	AP7i = FindChem("AP7", 0)
	AP8i = FindChem("AP8", 0)
	UR5i = FindChem("UR5", 0)
	UR6i = FindChem("UR6", 0)

	AP7DeltaBurden = (Avogadro/1.0e12)*(S10DeltaMassConc(1)/GasMolecularMass(AP7i))
	AP7NewBurden =  GetGridCellChemBurden( AP7i)-AP7DeltaBurden
	CALL ReplaceGridCellChemBurden ( AP7i, AP7NewBurden)

	AP8DeltaBurden = (Avogadro/1.0e12)*(S10DeltaMassConc(2)/GasMolecularMass(AP8i))
	AP8NewBurden =  GetGridCellChemBurden( AP8i)-AP8DeltaBurden
	CALL ReplaceGridCellChemBurden ( AP8i, AP8NewBurden)
	
	UR5DeltaBurden = (Avogadro/1.0e12)*(S10DeltaMassConc(3)/GasMolecularMass(UR5i))
	UR5NewBurden =  GetGridCellChemBurden( UR5i)-UR5DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR5i, UR5NewBurden)
	
	UR6DeltaBurden = (Avogadro/1.0e12)*(S10DeltaMassConc(4)/GasMolecularMass(UR6i))
	UR6NewBurden =  GetGridCellChemBurden( UR6i)-UR6DeltaBurden
	CALL ReplaceGridCellChemBurden ( UR6i, UR6NewBurden)

	RETURN
END SUBROUTINE UpdateGasPhaseOrganics

SUBROUTINE EquilibrateHydrophobicPhaseAtGridPoint()

	USE Aerosols,    ONLY : SortAerosolAtGridPointForCondensation, &
				Particles,			&
				Particle, &
				DumpParticleContentsAtError
	
	USE Chemistry,	 ONLY : HowManyEvolveGasChems, &
				GasMolecularMass
	
	USE GridPointFields, ONLY : GetTemp,	&
			            GetGridCellChemBurden,	&
				    ReplaceGridCellChemBurden

	USE ModelParameters,     ONLY : AqThermoNumbAllowedIterations,	&
					AqThermoEquilibriumError, &
					ThermoBulkMode, &
					Avogadro

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR
	
	IMPLICIT NONE

	INTEGER ::  HowManyLinks
	REAL*8  :: S1MassFracs(2), S2MassFracs(6), S3MassFracs(2), &
			S4MassFracs(3), S5MassFracs(2), S6MassFracs(3), &
			S7MassFracs(8), S8MassFracs(4), S9MassFracs(4), &
			S10MassFracs(4), InitialBurden(10) 
	TYPE(Particle),POINTER :: Current

	LOGICAL :: Equilibrated, RxnEquilibrated
	REAL*8  :: GasPhaseBurden, Temperature,  II
	INTEGER :: I, J, K, L, M, Q, ReturnType

	
	!Calculate surrogate concentrations and 
	!mass fractions of organics contributing to surrogates, and
	!Store initial surrogate gas phase burdens in a dummy vector
	CALL CalculateSurrogateGasConc(S1MassFracs, S2MassFracs, &
                S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)

	!! Sort the gridcell, pushing all empty sections to the end
	CALL SortAerosolAtGridPointForCondensation ()

	!! Count the number of links
	HowManyLinks = 0
	Current => Particles%First
	IF (ASSOCIATED(Current)) THEN
10		IF (Current%NumberOfParticles .GT. 0 ) &
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
		DO I = 1, HowManyOrganicDissolutionReactions						 
			
			!WRITE (*,*) "Reaction: ", I
			L = 0
			RxnEquilibrated = .FALSE.
				
			GasPhaseBurden = GetGridCellChemBurden ( INT(OrganicDissolutionData(I,1)))
		!WRITE(*,*) GasPhaseBurden
			!! Loop over all of the particles until 
			!! this particular reaction is happy
			DO WHILE (.NOT.RxnEquilibrated)

				RxnEquilibrated = .TRUE.

				Current => Particles%First
				DO J = 1, HowManyLinks
									
					!Equilibrate one reaction
					!WRITE(*,*) "Rxn: ", I
					CALL EquilibrateHydrophobicReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, OptErrorTolerance=AqThermoEquilibriumError/2.)							
					
					!! Make a more refined guess if still equilibrating late in the game
					!! (That is, force it to equilibrate much closer to 1 each time through)
					IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
						CALL EquilibrateHydrophobicReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
					END IF

					!WRITE(*,*) "Past Equil"
					IF (ReturnType .NE. 2 ) THEN
						RxnEquilibrated = .FALSE.
						Equilibrated    = .FALSE.
					END IF

					!WRITE(*,*) "Call UNIFAC"
					!Recalculate activity coefficients
					CALL UpdateHydrophobicUNIFAC(Current, Temperature)
					!WRITE(*,*) "Past UNIFAC"

					L = L+1
					Current => Current%Next
				END DO

				IF (L .GT. AqThermoNumbAllowedIterations) THEN
					CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
					CALL ERROR ("EquilibrateHydrophobicPhaseAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all "     // &
					            "of the particles (stored as "//TRIM(INT2STR(HowManyLinks))//" links) at a grid point "  // &
								"to organic dissolution reaction "//TRIM(INT2STR(I))//" (As indexed "// &
								"as input in OrganicDissolution.in) and failed.  The contents of the first particle in that "   // &
								"section are printed as UnEquilableParticle-Hydrophobic.txt")
				END IF
			END DO
			
			!! Replace the surrogate gas phase chemical concentrations
			IF (INT(OrganicDissolutionData(I,1)) .LE. HowManyEvolveGasChems)										&
				CALL ReplaceGridCellChemBurden ( INT(OrganicDissolutionData(I,1)), GasPhaseBurden)			
		END DO



		IF (K .GT. AqThermoNumbAllowedIterations) THEN
			CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
			CALL ERROR ("EquilibrateHydrophobicPhaseAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all of the hydrophobic dissolution "// &
					    "reactions and failed.  Usually this implies an "// &
						"unrealistic starting condition.  You may try to fix this by increasing AqThermoNumbAllowedIterations,"// &
						" the maximum  allowed iterations of this type of routine.")
		END IF

	END DO

	!WRITE(*,*) "Before org update"
	!Update Concentrations of individual organics
	CALL UpdateGasPhaseOrganics( S1MassFracs, S2MassFracs, &
                S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)
	!WRITE(*,*) "aFTER org update"

	RETURN

END SUBROUTINE EquilibrateHydrophobicPhaseAtGridPoint




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! --- EQUILIBRATE ONE HYDROPHOBIC DISSOLUTION REACTION ---    !!
!!							       !!
!! Iterate a particular reaction until it reaches equilibrium. !!
!!							       !!
!! ReturnType: 1 is normal, 2 is "already was in eq"	       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateHydrophobicReaction (GasPhaseBurden, WhichRxn, &
            InParticle, Temperature, ReturnType, OptErrorTolerance)

	USE Chemistry,		ONLY :	AqEquilibriaList,		&
					HowManyOrgChems,		&
					OrgPhaseChemicalNames,	&
					HowManyEvolveGasChems

	USE ModelParameters,	ONLY :	AqThermoEquilibriumError,	&
					WaterContentPrecision,		&
					AqThermoNumbAllowedIterations,	&
					AerosolWaterEquilibriumRHThreshold,&
					Avogadro,&
					moles, grams, WaterMolecMass, micron 

	USE GridPointFields,	ONLY :  GetRelativeHumidity

	USE InfrastructuralCode,ONLY :	INT2STR, REAL2STR,		&
					Transcript,			&
					ERROR, WARN

	USE Aerosols,		ONLY :	EquilibriumWaterContent,	 &
					DumpParticleContentsAtError, &
					EqCoeff,			&
					MeanActivity,&
					EquilibrateAqueousDissociationReaction

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReturnType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature
	REAL*8 :: OptErrorTolerance

	!! Local Variables
	INTEGER :: I
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, ActivityCoefficients,&
                   ErrorTolerance, DissociatedScaling, DissociatedExp
	LOGICAL :: WaterRxn

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .TRUE. ! .TRUE.

	IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN


	ErrorTolerance = OptErrorTolerance
	
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
!! There is only one reaction type:					     !!
!!									     !!
!! 8. Direct dissolution of a single species not inclorporating dissociation.!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Reaction Check 1", GasPhaseBurden, WhichRxn, Temperature
	EqRatio = HydrophobicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!WRITE(*,*) "Reaction Check 2"
	IF  (ABS(EqRatio-1.) .LT. AqThermoEquilibriumError) THEN
		ReturnType = 2
		RETURN
	END IF

	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) .LT. 1e-30) THEN !There is no organic left
		ReturnType = 2
		RETURN
	END IF

	!Check if all organic has condensed, but org-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The gas-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-30) THEN !There is no organic left in gas
		ReturnType = 2
		RETURN
	END IF

	Qnumer = InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3)))
	Qdenom = GasPhaseBurden / Avogadro  / InParticle%NumberOfParticles
	!WRITE(*,*) "Reaction Check 3", Qnumer, Qdenom

	!WRITE(*,*) "Check 4"
	!! I counts the number of iterations
	I = 1

	!! STEP 2: Calculate the Mass Step Sizes by which to Correct Concentrations
	Z  = (Qdenom + Qnumer) / 2.
	dX = Qdenom - Z				!! The Mass Flux Factor
	
	!! Loop over Corrections until convergence if there is anything to do
	IF (.NOT. (Qnumer .EQ. 0 .AND. Qdenom .EQ. 0)) THEN
	DO WHILE ((ABS(EqRatio-1) .GT. ErrorTolerance))

	!WRITE(*,*) "Check 5"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 3: Adjust Molalities: !!
	!!							  !!
	!! -- Ions First			  !!

	InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) = InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) + dX
	IF (OrganicDissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems) &
			GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles

	!WRITE(*,*) InParticle%OrgChems(OrganicDissolutionData(WhichRxn,3)), GasPhaseBurden
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) "Check 7"
	Z = 0.5*Z

	!! The recalculation is based on the equilibrium ratio

	!WRITE(*,*) "Check 8", GasPhaseBurden, WhichRxn
	EqRatio = HydrophobicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!WRITE(*,*) WhichRxn, EqRatio
	!WRITE(*,*) "Check 9"

	!WRITE(*,*) "Check 10"
	IF (EqRatio .GE. 1.+ErrorTolerance) dX = -1.*Z
	IF (EqRatio .LE. 1.-ErrorTolerance) dX = Z
	!IF (WhichRxn .EQ. 9) WRITE(*,*) EqRatio, InParticle%OrgChems(OrganicDissolutionData(WhichRxn,3)), GasPhaseBurden

	!! Count the number of attempts before convergence
	I = I+1

	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) .LT. 1e-30) THEN !There is no organic left
		RETURN
	END IF

	!Check if all organic has condensed, but org-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The gas-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-30) THEN !There is no organic left in gas
		RETURN
	END IF
	
	!Error if no convergence
	IF (I .GT. 500) THEN
		CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle-HydrophobicDissolutionError.txt", &
											 InDoSurfaceTension=.FALSE.,InDoWaterActivity=.FALSE., InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., &
                                                     InDensity=.TRUE.,  InIonicStrength=.TRUE.)
		WRITE(*,*) "Org Reaction #: ", WhichRxn, "Eq. Ratio: ", EqRatio
		!WRITE(*,*) GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles
		CALL ERROR("Apparent Non Convergence in EquilibrateHydrophobicReaction().  Trying to converge "// &
				   "towards equilibrium in hydrophobic dissolution reaction number "//TRIM(INT2STR(WhichRxn))// ".")

	END IF

	!WRITE(*,*) "Check 11"
	END DO ; END IF

	RETURN
END SUBROUTINE EquilibrateHydrophobicReaction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Check the Ratio between the Existing and Ideal Equilibrium   !!
!! Constants for hydrophobic reactions.                         !!
!! This is used by several other routines.		        !!
!!								!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This should only be used in the context of EquilibrateHydrophobicReaction!  
!! Otherwise there is no check for lack of the species.
REAL*8 FUNCTION HydrophobicEquilibriumConstantsRatio (GasPhaseBurden, & 
                 WhichRxn, Temperature, InParticle)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyOrgChems

	USE ModelParameters, ONLY : RStarMB, moles, grams, &
                                    WaterMolecMass, Avogadro

	USE InfrastructuralCode, ONLY : Error

	USE Aerosols, ONLY : OrgCurvatureCorrection

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReactionType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature, ActCoeff, OrgMolFrac

	!! Internal Variables
	REAL*8	:: Pvap, P1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There is only one reaction type:					     !!
!!									     !!
!! 8. Direct dissolution of a single species not incorporating dissociation. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! See if the reaction is empty !!
	IF (InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) .EQ. 0. & 
            .AND. GasPhaseBurden .EQ. 0) THEN
		
		HydrophobicEquilibriumConstantsRatio = 1.
		RETURN
	END IF
	!WRITE(*,*) "Flag 2"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Prepare the needed quantities !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Flag 3"
	!! Get the mole fraction of the condensed organic,
	!! the activity coefficient of the condensed organic,
	!! the saturation vapor pressure,
	!! and the partial pressure and the equilibrium at this temperature
	Pvap = SaturationVaporPressure (Temperature, WhichRxn)
	
	!WRITE(*,*) "Flag 3a"
	ActCoeff = InParticle%HydrophobicActivityCoeffs(INT(OrganicDissolutionData(WhichRxn,3)))
	!WRITE(*,*) WhichRxn, ActCoeff
	!WRITE(*,*) "Flag 3b"
	OrgMolFrac = InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3)))/TotalMolesOM(InParticle)
	!WRITE(*,*) "Flag 3c"
	P1 = GasPhaseBurden*RstarMB*Temperature/Avogadro
	!WRITE(*,*) "Flag 3d"	
	!! If the solution is set to blow up, then establish a large artificial value.
	IF (GasPhaseBurden .LE. 0.) THEN
		HydrophobicEquilibriumConstantsRatio = 1.e12
		RETURN
	END IF
	!WRITE(*,*) "Flag 4"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make the primary calculation !!	
		
	HydrophobicEquilibriumConstantsRatio = OrgCurvatureCorrection(InParticle)*OrgMolFrac*ActCoeff*Pvap/P1
				
	!IF(WhichRxn .EQ. 23) WRITE(*,*) "Phobic Ratio", HydrophobicEquilibriumConstantsRatio, GasPhaseBurden
	!WRITE(*,*) "Ratio", HydrophobicEquilibriumConstantsRatio
	!WRITE(*,*) "Flag 8"
	RETURN
END FUNCTION HydrophobicEquilibriumConstantsRatio


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! --- EQUILIBRATE ONE HYDROPHIlIC DISSOLUTION REACTION ---    !!
!!							       !!
!! Iterate a particular reaction until it reaches equilibrium. !!
!!							       !!
!! ReturnType: 1 is normal, 2 is "already was in eq"	       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateHydrophilicReaction (GasPhaseBurden, WhichRxn, & 
            InParticle, Temperature, ReturnType, OptErrorTolerance)

	USE Chemistry,	ONLY :	AqEquilibriaList,		&
				HowManyOrgChems,		&
				OrgPhaseChemicalNames,	&
				HowManyEvolveGasChems

	USE ModelParameters,	ONLY :	AqThermoEquilibriumError,	&
					WaterContentPrecision,		&
					AqThermoNumbAllowedIterations,	&
					AerosolWaterEquilibriumRHThreshold,&
					Avogadro,			&
					moles, grams, WaterMolecMass, micron 

	USE GridPointFields,	ONLY :  GetRelativeHumidity

	USE InfrastructuralCode,	ONLY :	INT2STR, REAL2STR,	&
						Transcript,		&
						ERROR, WARN

	USE Aerosols,		ONLY :	EquilibriumWaterContent,	 &
					DumpParticleContentsAtError, &
					EqCoeff,			&
					MeanActivity,&
					EquilibrateAqueousDissociationReaction

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReturnType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature
	REAL*8 :: OptErrorTolerance

	!! Local Variables
	INTEGER :: I
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, ActivityCoefficients, ErrorTolerance, DissociatedScaling, DissociatedExp
	LOGICAL :: WaterRxn

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .TRUE. ! .TRUE.

	IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	!! Set the error tolerance according to the 
	!! subroutine's optinoal input.
	!! THIS IS NOT THE TOLERANCE FOR THE EQUILIBRATION
	!! CHECK (WHICH USES THE GLOBAL VALUE), RATHER IT 
	!! IS THE VALUE WE EQUILIBRATE TO IF THAT CHECK FAILS	
	!IF (PRESENT(OptErrorTolerance)) THEN
	!	ErrorTolerance = OptErrorTolerance
	!ELSE
		ErrorTolerance = AqThermoEquilibriumError / 2.
	!END IF
	
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
!! There is only one reaction type:					     !!
!!									     !!
!! 9. Direct dissolution of a single species not inclorporating dissociation.!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Reaction Check 1", GasPhaseBurden, WhichRxn, Temperature
	EqRatio = HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!WRITE(*,*) "Reaction Check 2"
	IF  (ABS(EqRatio-1.) .LT. AqThermoEquilibriumError) THEN
		ReturnType = 2
		RETURN
	END IF

	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) .LT. 1e-30) THEN !There is no organic left
		ReturnType = 2
		RETURN
	END IF

	!Check if all organic has condensed, but aq-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The aq-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-30) THEN !There is no organic left in gas
		ReturnType = 2
		RETURN
	END IF


	Qnumer = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3)))
	Qdenom = GasPhaseBurden / Avogadro / InParticle%NumberOfParticles
	!WRITE(*,*) "Reaction Check 3"

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
	!!							  !!
	!! -- Ions First			  !!

	InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) + dX
	IF (AqOrganicDissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems) &
			GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) "Check 7"
	Z = 0.5*Z

	!! The recalculation is based on the equilibrium ratio

	!WRITE(*,*) "Check 8", GasPhaseBurden, WhichRxn
	EqRatio = HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!WRITE(*,*) "Check 9"

	!WRITE(*,*) "Check 10"
	IF (EqRatio .GE. 1.+ErrorTolerance) dX = -1.*Z
	IF (EqRatio .LE. 1.-ErrorTolerance) dX = Z
	!WRITE(*,*) EqRatio, Z, dX

	!! Count the number of attempts before convergence
	I = I+1

	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) .LT. 1e-30) THEN !There is no organic left
		RETURN
	END IF

	!Check if all organic has condensed, but aq-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The aq-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-30) THEN !There is no organic left in gas
		RETURN
	END IF
	
	!Error if no convergence
	IF (I .GT. 500) THEN
		CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle-HydrophilicDissolutionError.txt", &
											 InDoSurfaceTension=.FALSE.,InDoWaterActivity=.FALSE., InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., &
                                                     InDensity=.TRUE.,  InIonicStrength=.TRUE.)
		WRITE(*,*) "Aq Org Reaction #: ", WhichRxn, "Eq. Ratio: ", EqRatio, InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3)))
		WRITE(*,*) GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles
		CALL ERROR("Apparent Non Convergence in EquilibrateHydrophilicReaction().  Trying to converge "// &
				   "towards equilibrium in hydrophilic dissolution reaction number "//TRIM(INT2STR(WhichRxn))// ".")

	END IF

	!WRITE(*,*) "Check 11"
	END DO ; END IF

	!WRITE(*,*) "Check 12"
	RETURN
END SUBROUTINE EquilibrateHydrophilicReaction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Check the Ratio between the Existing and Ideal Equilibrium   !!
!! Constants for hydrophilic reactions.                         !!
!! This is used by several other routines.		        !!
!!								!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This should only be used in the context of EquilibrateHydrophilicReaction!  
!! Otherwise there is no check for lack of the species.
REAL*8 FUNCTION HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, & 
                  WhichRxn, Temperature, InParticle)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyAqOrgChems

	USE ModelParameters, ONLY : RStarMB, moles, grams, & 
                                    WaterMolecMass, Avogadro, ProtonIndex

	USE InfrastructuralCode, ONLY : Error
	USE Aerosols, ONLY : Molality, CurvatureCorrection

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReactionType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature, ActCoeff, OrgMolality, & 
                  ProtonMolality

	!! Internal Variables
	REAL*8	:: Henry,  P1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There is only one reaction type:					     !!
!!									     !!
!! 9. Direct dissolution of a single species not incorporating dissociation. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! See if the reaction is empty !!

	!WRITE(*,*) "Flag 1"
	!! If there is none of the chemical to deal with
	!WRITE(*,*) WhichRxn
	!WRITE(*,*) InParticle%AqOrgChems(AqOrganicDissolutionData(WhichRxn,3))
	!WRITE(*,*) GasPhaseBurden
	IF (InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) .EQ. 0. .AND.	&
		GasPhaseBurden .EQ. 0) THEN
		
		HydrophilicEquilibriumConstantsRatio = 1.
		RETURN
	END IF
	!WRITE(*,*) "Flag 2"

	!! If the solution is set to blow up, then establish a large artificial value.
	IF (GasPhaseBurden .LE. 0.) THEN
		HydrophilicEquilibriumConstantsRatio = 1.e12
		RETURN
	END IF
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Prepare the needed quantities !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Flag 3"
	
	!Calculate the proton molality
	!WRITE(*,*) ProtonIndex
	ProtonMolality = Molality (ProtonIndex, InParticle)
	
	!Calculate the effective Henry's Law Constant (mol/kg H2O/mbar)
	Henry = EffectiveOrganicHenrysLaw  (Temperature, ProtonMolality, WhichRxn)
	
	!Retrieve the hydrophilic activity coefficient (corrected for infinite dilution convention)
	ActCoeff = InParticle%HydrophilicActivityCoeffs(INT(AqOrganicDissolutionData(WhichRxn,3)))

	!Partial pressure in mbar
	P1 = GasPhaseBurden*RstarMB*Temperature/Avogadro

	!Calculate molality of species in aqueous phase
	OrgMolality = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make the primary calculation !!
		
	HydrophilicEquilibriumConstantsRatio = CurvatureCorrection(InParticle)*ActCoeff*OrgMolality/(P1*Henry)
				
	!IF(WhichRxn .EQ. 23) WRITE(*,*) "Ratio", HydrophilicEquilibriumConstantsRatio, GasPhaseBurden
	!WRITE(*,*) "Flag 8"
	RETURN
END FUNCTION HydrophilicEquilibriumConstantsRatio


!Calculates effective henry's law for hydrophilic organics 
!in units of mol/kg/mbar 
REAL*8 FUNCTION EffectiveOrganicHenrysLaw  (Temperature, ProtonMolality, AqOrganicDissolutionRxn)

	USE ModelParameters, ONLY : EqCoeffTref
	USE InfrastructuralCode, ONLY : Error 

	IMPLICIT NONE

	REAL*8  :: Temperature, ProtonMolality, KA1, KA2, Henry, AcidFactor
	INTEGER :: AqOrganicDissolutionRxn, I
	
	I = AqOrganicDissolutionRxn
	
	!in mol/kg
	KA1 = AqOrganicDissolutionData(i,7)*EXP((AqOrganicDissolutionData(i,8)/0.001097) &
											*((1.0/Temperature)-(1.0/298.15)))
	!in mol/kg
	KA2 = AqOrganicDissolutionData(i,9)*EXP((AqOrganicDissolutionData(i,12)/0.001097) &
											*((1.0/Temperature)-(1.0/298.15)))
	
	!Correction to Henry's law due to acid dissociation
	AcidFactor = 1 + KA1/ProtonMolality + KA1*KA2/(ProtonMolality**2)
		
	IF(AcidFactor .LT. 1.0) CALL ERROR("Error in EffectiveOrganicHenrysLaw - corrected value is less than uncorrected.")

	!in (ug i/ ug H2O) / (ug i / m3 air)
	Henry = AqOrganicDissolutionData(i,5)
	
	!Conversion to units of mol i/ kg H20 / mbar
	Henry = (Henry*1e6)/8.314e-5/298.15
	
	!Temperature Dependence
	Henry = Henry*EXP((AqOrganicDissolutionData(i,6)/0.001097) &
						*((1.0/Temperature)-(1.0/298.15)))
		
	EffectiveOrganicHenrysLaw  = Henry*AcidFactor
	
	RETURN
END FUNCTION EffectiveOrganicHenrysLaw

!! This subroutine provides the link between the main model and the 
!! UNIFAC hydrophilic driver file soa_unidriva.f
!! Note that the order of the hydrophilic organic compounds in 
!! unifacparama.h and in HydrophilicOrgChems.in must be the same! 
SUBROUTINE UpdateHydrophilicUNIFAC(InParticle, Temperature)

	USE Chemistry,				ONLY :	HowManyAqOrgChems

	IMPLICIT NONE

	INTEGER :: I, N
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8, ALLOCATABLE :: XPASS(:), GAMMA(:)
	REAL*8 :: TotalMoles, Temperature, InfDiluActCoeff(1:23)
	

		!Calculate infinite dilution activity coefficients
		CALL CalcInfiniteDilutionHydrophilicUNIFAC(Temperature, InfDiluActCoeff)
		
		N = 24 !Number of compounds
		ALLOCATE(XPASS(0:N-1), GAMMA(0:N-1))
	
		!Calculate mole fractions, putting them in an array ordered from 0:N-1
		TotalMoles = TotalAqMoles(InParticle)
		
		DO I = 0, N-2
			XPASS(I) = InParticle%AqOrgChems(I+1)/TotalMoles
			!WRITE(*,*) I+1, XPASS(I)
		END DO

		!Water is the last conpound (index 23)
			XPASS(N-1) = InParticle%AqChems(1)/TotalMoles

		!Call UNIFAC, placing activity coefficients in GAMMA(0:N-1)
		CALL UNIDRIVA(XPASS,GAMMA,N,Temperature)
	
		!Place activity coefficients into particle array
		DO I = 0, N-2
			InParticle%HydrophilicActivityCoeffs(I+1) = GAMMA(I)/InfDiluActCoeff(I+1)
		END DO


END SUBROUTINE UpdateHydrophilicUNIFAC

!Moles per particle of hydrophilic organics and water
REAL*8 FUNCTION TotalAqMoles(InParticle) 

	USE Chemistry,				ONLY :	HowManyAqOrgChems
	
	IMPLICIT NONE

	REAL*8  :: Total
	INTEGER :: I
	TYPE(PARTICLE),POINTER :: InParticle

	
	Total = InParticle%AqChems(1) !Water Content
	
	DO I = 1, HowManyAqOrgChems
		Total = Total + InParticle%AqOrgChems(I)
	END DO
	
	TotalAqMoles = Total

	RETURN

END FUNCTION TotalAqMoles

!! This subroutine provides the link between the main model and the 
!! UNIFAC hydrophilic driver file soa_unidriva.f
!! Note that the order of the hydrophilic organic compounds in 
!! unifacparama.h and in HydrophilicOrgChems.in must be the same! 
SUBROUTINE CalcInfiniteDilutionHydrophilicUNIFAC(Temperature, InfDiluActCoeff)

	USE Chemistry,				ONLY :	HowManyAqOrgChems

	IMPLICIT NONE

	INTEGER :: I, N, J
	REAL*8, ALLOCATABLE :: XPASS(:), GAMMA(:)
	REAL*8 :: TotalMoles, Temperature, InfDiluActCoeff(1:23)
	

		N = 24 !Number of compounds
		ALLOCATE(XPASS(0:N-1), GAMMA(0:N-1))
	
		DO I = 1, N-1
			
			!Assign mole fractions
			DO J = 1,N
				IF (J == I) THEN !Dilute Compound
					XPASS(J-1) = 1.0e-6
				ELSE IF (J == N) THEN !Water
					XPASS(J-1) = 1.0 - 1.0e-6
				ELSE
					XPASS(J-1) = 0.0
				END IF
			END DO

			!Call UNIFAC, placing activity coefficients in GAMMA(0:N-1)
			CALL UNIDRIVA(XPASS,GAMMA,N,Temperature)
		
			InfDiluActCoeff(I) = GAMMA(I-1)
		END DO

		RETURN

END SUBROUTINE CalcInfiniteDilutionHydrophilicUNIFAC

SUBROUTINE EquilibrateAllOrganicsAtGridPoint()

	USE Aerosols,    ONLY : SortAerosolAtGridPointForCondensation, &
				Particles,				&
				Particle, &
				DumpParticleContentsAtError
	
	USE Chemistry,	 ONLY : HowManyEvolveGasChems, &
				GasMolecularMass
	
	USE GridPointFields, ONLY : GetTemp,	&
				    GetGridCellChemBurden,		&
				    ReplaceGridCellChemBurden

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterations,	&
					AqThermoEquilibriumError, &
					ThermoBulkMode, &
					Avogadro

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR
	
	IMPLICIT NONE

	INTEGER ::  HowManyLinks
	REAL*8  :: S1MassFracs(2), S2MassFracs(6), S3MassFracs(2), &
			S4MassFracs(3), S5MassFracs(2), S6MassFracs(3), &
			S7MassFracs(8), S8MassFracs(4), S9MassFracs(4), &
			S10MassFracs(4), InitialBurden(10) 
	TYPE(Particle),POINTER :: Current

	LOGICAL :: Equilibrated, RxnEquilibrated
	REAL*8  :: GasPhaseBurden, Temperature, II
	INTEGER :: I, J, K, L, M, Q, ReturnType, AqRxnIndex

	
	!Calculate surrogate concentrations and 
	!mass fractions of organics contributing to surrogates, and
	!Store initial surrogate gas phase burdens in a dummy vector
	CALL CalculateSurrogateGasConc( S1MassFracs, S2MassFracs, &
                S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)

	!! Sort the gridcell, pushing all empty sections to the end
	CALL SortAerosolAtGridPointForCondensation ()

	!! Count the number of links
	HowManyLinks = 0
	Current => Particles%First
	IF (ASSOCIATED(Current)) THEN
10		IF (Current%NumberOfParticles .GT. 0 ) &
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
		DO I = 1, HowManyOrganicDissolutionReactions+HowManyAqOrganicDissolutionReactions						 
			
			!WRITE (*,*) "Reaction: ", I
			L = 0
			RxnEquilibrated = .FALSE.
				
			!Hydrophobic Reactions
			IF (I .LE. HowManyOrganicDissolutionReactions) THEN
			
				GasPhaseBurden = GetGridCellChemBurden ( INT(OrganicDissolutionData(I,1)))
				!WRITE(*,*) GasPhaseBurden
				!! Loop over all of the particles until 
				!! this particular reaction is happy
				DO WHILE (.NOT.RxnEquilibrated)

					RxnEquilibrated = .TRUE.

					Current => Particles%First
					DO J = 1, HowManyLinks
									
						!Equilibrate one reaction
						!WRITE(*,*) "Rxn: ", I
						CALL EquilibrateHydrophobicReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, OptErrorTolerance=AqThermoEquilibriumError/2.)							
					
						!! Make a more refined guess if still equilibrating late in the game
						!! (That is, force it to equilibrate much closer to 1 each time through)
						IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
							CALL EquilibrateHydrophobicReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
						END IF

						!WRITE(*,*) "Past Equil"
						IF (ReturnType .NE. 2 ) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF

						!WRITE(*,*) "Call UNIFAC"
						!Recalculate activity coefficients
						CALL UpdateHydrophobicUNIFAC(Current, Temperature)
						!WRITE(*,*) "Past UNIFAC"

							L = L+1
					Current => Current%Next
					END DO

					IF (L .GT. AqThermoNumbAllowedIterations) THEN
						CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
						CALL ERROR ("EquilibrateAllOrganicsAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all "     // &
					            "of the particles (stored as "//TRIM(INT2STR(HowManyLinks))//" links) at a grid point "  // &
								"to hydrophobic organic dissolution reaction "//TRIM(INT2STR(I))//" (As indexed "// &
								"as input in OrganicDissolution.in) and failed.  The contents of the first particle in that "   // &
								"section are printed as UnEquilableParticle-Hydrophobic.txt")
					END IF
				END DO
			
				!! Replace the surrogate gas phase chemical concentrations
				IF (INT(OrganicDissolutionData(I,1)) .LE. HowManyEvolveGasChems)										&
					CALL ReplaceGridCellChemBurden ( INT(OrganicDissolutionData(I,1)), GasPhaseBurden)			
		
			!Hydrophilic Reactions
			ELSE
		
				AqRxnIndex = I - HowManyOrganicDissolutionReactions
				GasPhaseBurden = GetGridCellChemBurden ( INT(AqOrganicDissolutionData(AqRxnIndex,1)))
				!WRITE(*,*) AqRxnIndex
				!IF (AqRxnIndex .EQ. 1) WRITE(*,*) GasPhaseBurden
				!! Loop over all of the particles until 
				!! this particular reaction is happy
				DO WHILE (.NOT.RxnEquilibrated)

					RxnEquilibrated = .TRUE.

					Current => Particles%First
					DO J = 1, HowManyLinks
									
						!Equilibrate one reaction
						!WRITE(*,*) "Rxn: ", I, AqRxnIndex
						CALL EquilibrateHydrophilicReaction (GasPhaseBurden, AqRxnIndex, Current, Temperature, ReturnType, OptErrorTolerance=AqThermoEquilibriumError/2.)							
						!WRITE(*,*) "Equilibrate Okay"

						!! Make a more refined guess if still equilibrating late in the game
						!! (That is, force it to equilibrate much closer to 1 each time through)
						IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
							CALL EquilibrateHydrophilicReaction (GasPhaseBurden, AqRxnIndex, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
						END IF

						!WRITE(*,*) "Past Equil"
						IF (ReturnType .NE. 2 ) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF

						!WRITE(*,*) "Call UNIFAC"
						!Recalculate activity coefficients
						CALL UpdateHydrophilicUNIFAC(Current, Temperature)
						!WRITE(*,*) "Past UNIFAC"

							L = L+1
					Current => Current%Next
					END DO

					IF (L .GT. AqThermoNumbAllowedIterations) THEN
						CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
						CALL ERROR ("EquilibrateAllOrganicsAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all "     // &
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
		END DO



		IF (K .GT. AqThermoNumbAllowedIterations) THEN
			CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
			CALL ERROR ("EquilibrateHydrophobicPhaseAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all of the hydrophobic dissolution "// &
					    "reactions and failed.  Usually this implies an "// &
						"unrealistic starting condition.  You may try to fix this by increasing AqThermoNumbAllowedIterations,"// &
						" the maximum  allowed iterations of this type of routine.")
		END IF

	END DO

	!WRITE(*,*) "Before org update"
	!Update Gas Phase Concentrations of individual organics
	CALL UpdateGasPhaseOrganics(S1MassFracs, S2MassFracs, &
                S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)
	!WRITE(*,*) "After org update"

	RETURN

END SUBROUTINE EquilibrateAllOrganicsAtGridPoint

REAL*8 function dev(xin, Q, RH, Temp)
   
   IMPLICIT NONE
   
   !External Variables
   real*8, INTENT(IN) :: xin, RH, Temp
   integer, INTENT(IN) :: Q

   !Internal Variables
   REAL*8 :: GAMMA(0:23), XPASS(0:23)
	 
	 XPASS(Q) = xin
	 XPASS(23) = 1.0 - XPASS(Q)
	 
	 CALL UNIDRIVA(XPASS,GAMMA,24,Temp)
	 
	 dev=GAMMA(23)*XPASS(23)-RH/100.
	 return

end function dev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate until a single organic in two phases (org and aqueous)!! 
!! is in equilibrium.                                            !!
!!								 !!
!! ReturnType: 1 is normal, 2 is "already was in eq"		 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateOrganicReaction (WhichAqRxn, WhichOrgRxn, &
               InParticle, Temperature, ReturnType, OptErrorTolerance)

	USE Chemistry,	ONLY :	AqEquilibriaList,		&
				HowManyOrgChems,		&
				OrgPhaseChemicalNames,	&
				HowManyEvolveGasChems

	USE ModelParameters
			
	USE GridPointFields,		ONLY : 	GetRelativeHumidity

	USE InfrastructuralCode,	ONLY :	INT2STR, REAL2STR,	&
						Transcript,	&
						ERROR, WARN

	USE Aerosols,		ONLY :	EquilibriumWaterContent, &
					DumpParticleContentsAtError, &
					EqCoeff,		&
					MeanActivity, &
					EquilibrateAqueousDissociationReaction

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichAqRxn, WhichOrgRxn, ReturnType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: Temperature
	REAL*8 :: OptErrorTolerance

	!! Local Variables
	INTEGER :: I
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, &
                    ActivityCoefficients, ErrorTolerance, &
                    DissociatedScaling, DissociatedExp
	LOGICAL :: WaterRxn

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .TRUE. ! .TRUE.

	!If there are no particles, dissolution isn't allowed into one phase
	!or org phase is empty, skip this step
	IF (InParticle%NumberOfParticles .EQ. 0.) RETURN
	IF (NOT(DoHydrophobicOrgDissolutionFlag) .OR. NOT(DoHydrophilicOrgDissolutionFlag)) RETURN
	IF (TotalMolesOM(InParticle) .LE. 0.0) RETURN 

	!! Set the error tolerance according to the 
	!! subroutine's input.
	ErrorTolerance = OptErrorTolerance

	
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

	!If none of the organic is present, return
	IF(InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn, 3))) .EQ. 0.0 &
	   .AND. InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn, 3))) .EQ. 0.0) THEN
		ReturnType = 2
		RETURN
	END IF

	!Check if reaction is already in equilibrium
	EqRatio = TwoPhaseOrgEqRatio(WhichAqRxn, WhichOrgRxn, InParticle, Temperature)

	IF  (ABS(EqRatio-1.) .LT. AqThermoEquilibriumError) THEN
		ReturnType = 2
		RETURN
	END IF

	!If the concentration in either phase is going to 0, and the
	!other phase is subsaturated, return.
	IF((InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))) .LT. 1.0e-30 &
		.AND. EqRatio .GT. 1.0) .OR. &
		(InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) .LT. 1.0e-30 &
		.AND. EqRatio .LT. 1.0) ) THEN
		ReturnType = 2
		RETURN
	END IF

	Qnumer = InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3)))
	Qdenom = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3)))
	!WRITE(*,*) "Reaction Check 3", Qnumer, Qdenom

	!WRITE(*,*) "Check 4"
	!! I counts the number of iterations
	I = 1

	!! STEP 2: Calculate the Mass Step Sizes by which to Correct Concentrations
	Z  = (Qdenom + Qnumer) / 2.
	dX = Qdenom - Z				!! The Mass Flux Factor
	
	!! Loop over Corrections until convergence if there is anything to do
	IF (.NOT. (Qnumer .EQ. 0 .AND. Qdenom .EQ. 0)) THEN
	DO WHILE ((ABS(EqRatio-1) .GT. ErrorTolerance))

	!WRITE(*,*) "Check 5"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 3: Adjust Molalities: !!
	!!							  !!
	!! -- Ions First			  !!

	InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))) = InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))) + dX
	InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) - dX
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) "Check 7"
	Z = 0.5*Z

	!! The recalculation is based on the equilibrium ratio

	!WRITE(*,*) "Check 8", GasPhaseBurden, WhichRxn
	EqRatio = TwoPhaseOrgEqRatio(WhichAqRxn, WhichOrgRxn, InParticle, Temperature)
	!WRITE(*,*) WhichAqRxn, EqRatio
	!WRITE(*,*) "Check 9"

	!WRITE(*,*) "Check 10"
	IF (EqRatio .GE. 1.+ErrorTolerance) dX = -1.*Z
	IF (EqRatio .LE. 1.-ErrorTolerance) dX = Z
	
	!If the concentration in either phase is going to 0, and the
	!other phase is subsaturated, return.
	IF((InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))) .LT. 1.0e-30 &
		.AND. EqRatio .GT. 1.0) .OR. &
		(InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) .LT. 1.0e-30 &
		.AND. EqRatio .LT. 1.0) ) THEN
		RETURN
	END IF
	
	!! Count the number of attempts before convergence
	I = I+1

	!Error if no convergence
	IF (I .GT. 500) THEN
		WRITE(*,*) EqRatio, InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))), &
		            InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3)))
		CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle-HydrophobicDissolutionError.txt", &
											 InDoSurfaceTension=.FALSE.,InDoWaterActivity=.FALSE., InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE.,  &
                                                     InDensity=.TRUE., InIonicStrength=.TRUE.)
		CALL ERROR("Apparent Non Convergence in EquilibrateOrganicReaction.  Trying to converge "// &
				   "towards equilibrium in hydrophilic internal reaction number "//TRIM(INT2STR(WhichAqRxn))// ".")

	END IF

	!WRITE(*,*) "Check 11"
	END DO ; END IF

	RETURN
END SUBROUTINE EquilibrateOrganicReaction

REAL*8 FUNCTION TwoPhaseOrgEqRatio(WhichAqRxn, WhichOrgRxn, & 
                                   InParticle, Temperature)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyOrgChems

	USE ModelParameters, ONLY : RStarMB, moles, grams, WaterMolecMass, &
                                    Avogadro, &
				    ProtonIndex

	USE InfrastructuralCode, ONLY : Error

	USE Aerosols, ONLY : Molality, CurvatureCorrection, &
                             OrgCurvatureCorrection

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichAqRxn, WhichOrgRxn
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: Temperature

	!! Internal Variables
	REAL*8	:: Pvap, OrgActCoeff, OrgMolFrac, ProtonMolality, &
				Henry, AqOrgMolality, AqActCoeff
	

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Prepare the needed quantities !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Flag 3"
	!! Get the mole fraction of the condensed organic,
	!! the activity coefficient of the condensed organic,
	!! the saturation vapor pressure,
	!! and the partial pressure and the equilibrium at this temperature
	Pvap = SaturationVaporPressure (Temperature, WhichOrgRxn)
	
	OrgActCoeff = InParticle%HydrophobicActivityCoeffs(INT(OrganicDissolutionData(WhichOrgRxn,3)))

	OrgMolFrac = InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3)))/TotalMolesOM(InParticle)

	!Calculate the proton molality
	!WRITE(*,*) "ProtonIndex", ProtonIndex

	ProtonMolality = Molality (ProtonIndex, InParticle)
	
	!Calculate the effective Henry's Law Constant (mol/kg H2O/mbar)
	Henry = EffectiveOrganicHenrysLaw  (Temperature, ProtonMolality, WhichAqRxn)
	
	!Retrieve the hydrophilic activity coefficient (corrected for infinite dilution convention)
	AqActCoeff = InParticle%HydrophilicActivityCoeffs(INT(AqOrganicDissolutionData(WhichAqRxn,3)))

	!Calculate molality of species in aqueous phase
	AqOrgMolality = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

	!Prevents a return of infinity if there is no organic in aqueous phase
	IF(AqOrgMolality .EQ. 0.0) THEN
		TwoPhaseOrgEqRatio = 1.0e10
		RETURN
	END IF
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make the primary calculation !!	
		
	TwoPhaseOrgEqRatio = OrgCurvatureCorrection(InParticle)*OrgMolFrac*OrgActCoeff*Pvap*Henry/(CurvatureCorrection(InParticle)*AqActCoeff*AqOrgMolality)
				
	!WRITE(*,*) "Ratio", TwoPhaseOrgEqRatio
	!WRITE(*,*) "Flag 8"
	RETURN


END FUNCTION TwoPhaseOrgEqRatio

SUBROUTINE EquilibrateOrganicParticle (InParticle, Temperature, &
                                      ReturnType, OptErrorTolerance)

	USE Chemistry,	ONLY :	AqEquilibriaList,		&
				HowManyOrgChems,		&
				OrgPhaseChemicalNames,	&
				HowManyEvolveGasChems

	USE ModelParameters
			
	USE GridPointFields,	ONLY :  GetRelativeHumidity

	USE InfrastructuralCode,ONLY :	INT2STR, REAL2STR,	&
					Transcript,	&
					ERROR, WARN

	USE Aerosols,		ONLY :	EquilibriumWaterContent,	 &
					DumpParticleContentsAtError, &
					EqCoeff,		&
					MeanActivity,&
					EquilibrateAqueousDissociationReaction

	IMPLICIT NONE

	!! External Variables
	INTEGER :: ReturnType, allocation_error
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: Temperature
	REAL*8 :: OptErrorTolerance

	!! Local Variables
	INTEGER :: I, K, C, InternalReturnType
	REAL*8  :: II, ErrorTolerance
	INTEGER, ALLOCATABLE :: OrgIndex(:)

	!! SCSCSCSCSC
	LOGICAL :: Eq, Scaffolding = .TRUE. ! .TRUE.

	ReturnType = 2
	!If there are no particles, dissolution isn't allowed into one phase
	!or org phase is empty, skip this step
	IF (InParticle%NumberOfParticles .EQ. 0.) RETURN
	IF (NOT(DoHydrophobicOrgDissolutionFlag) .OR. NOT(DoHydrophilicOrgDissolutionFlag)) RETURN
	!WRITE(*,*) "Total OM in EquilibrateOrganicParticle ", TotalMolesOM(InParticle)
	IF (TotalMolesOM(InParticle) .LE. 0.0) RETURN 

	!! Set the error tolerance according to the input
	ErrorTolerance = OptErrorTolerance
	
	!! Presume innocence in the form of a successful return.
	ReturnType = 1

	!! Allocate the LSODES working arrays
	ALLOCATE (OrgIndex(HowManyAqOrganicDissolutionReactions), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of OrgIndex could not proceed in EquilibrateOrganicParticle")

	
	!Find matching hydrophobic reaction for hydrophilic reaction
	DO I = 1, HowManyAqOrganicDissolutionReactions
		OrgIndex(I) = 0
		DO C = 1, HowManyOrganicDissolutionReactions
			IF(OrganicDissolutionData(C,1) .EQ. AqOrganicDissolutionData(I,1)) THEN
				OrgIndex(I) = C
			END IF
		END DO
	END DO

	
	K = 0
	Eq = .FALSE.

	DO WHILE(NOT(Eq))

		ReturnType = 2
		Eq = .TRUE.
		DO I = 1, HowManyAqOrganicDissolutionReactions
			CALL EquilibrateOrganicReaction (I, OrgIndex(I), InParticle, Temperature, InternalReturnType, ErrorTolerance)
			IF(InternalReturnType .NE. 2) Eq = .FALSE.
			IF(K .EQ. 0 .AND. InternalReturnType .NE. 2) ReturnType = 1
		END DO
	
	
	!! Count the number of attempts before convergence
	K = K+1

	!Error if no convergence
	!Actually, just exit (MJA, 091107)
	IF (K .GT. 500) THEN
		Eq = .TRUE.
	END IF

	END DO

	RETURN
END SUBROUTINE EquilibrateOrganicParticle

REAL*8 FUNCTION BinaryWaterActivityOrganic(InParticle, WhichRxn, Temp)
!! This function calculates the water activity for a binary solution 
!! of aqueous organic and water (uding UNIFAC) when the organic 
!! is at the same 
!! molality as in the complex solution. This is later used to 
!! calculate the total water activity of the 
!! organic-inorganic solution.
!! Created by Matt Alvarado, 09/08/2006

	IMPLICIT NONE

	!Input Variables
	TYPE(PARTICLE),POINTER :: InParticle
	INTEGER :: WhichRxn
	REAL*8 :: Temp !Temperature

	!Internal Variables
	INTEGER :: I, J, N
	REAL*8, ALLOCATABLE :: XPASS(:), GAMMA(:)
	REAL*8 :: xorg, Ratio

	!I is the index of the organic in AqOrgChems
	I = AqOrganicDissolutionData(WhichRxn,3)
	Ratio = InParticle%AqOrgChems(I) &
	           / InParticle%AqChems(1)

	!Mole fraction of organic in binary solution at same molality
	!as complex solution
	xorg = Ratio/(1 + Ratio)

	N = 24 !Number of compounds
	ALLOCATE(XPASS(0:N-1), GAMMA(0:N-1))
	!Assign mole fractions
	DO J = 1,N
		IF (J == I) THEN !Organic Compound
			XPASS(J-1) = xorg
		ELSE IF (J == N) THEN !Water
			XPASS(J-1) = 1.0 - xorg
		ELSE
			XPASS(J-1) = 0.0
		END IF
	END DO

	!Call UNIFAC, placing activity coefficients in GAMMA(0:N-1)
	CALL UNIDRIVA(XPASS,GAMMA,N,Temp)
		
	BinaryWaterActivityOrganic = GAMMA(N-1)*XPASS(N-1)

	RETURN
	
END FUNCTION BinaryWaterActivityOrganic
