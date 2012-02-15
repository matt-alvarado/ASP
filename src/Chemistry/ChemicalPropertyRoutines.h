!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ChemicalPropertyRoutines.h
!! This calculates the molecular diffusion coefficient of all gases.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History	
!! 10/15  2010   Matt Alvarado     Removed Optional Arguments						!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:					!!
!! 1. FUNCTION MolecularDiffusionCoefficient (XP,YP,ZP, ChemIndex, ChemName)	!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! DIFFUSIVITY OF TRACE CHEMICAL SPECIES as a function of temperature and pressure. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION MolecularDiffusionCoefficient (XP,YP,ZP, ChemIndex)

		USE ModelParameters,     ONLY : cm,AirMolecMass,Angstrom,Rstar,Avogadro,Pi,ChemScale
		USE InfrastructuralCode, ONLY : Warn, Error, INT2STR
		USE GridPointFields,	 ONLY : GetAirDensity, GetTemp
		IMPLICIT NONE

		!! External Variables
		REAL*8 :: XP,YP,ZP
		INTEGER		:: ChemIndex
		!CHARACTER*(*), OPTIONAL :: ChemName

		!! Internal Variables
		INTEGER :: i

		!IF (PRESENT(ChemIndex)) 
                 i = ChemIndex
		!IF (PRESENT(ChemName))  i = FindChem(ChemName,GasPhase)

		!IF (PRESENT(ChemIndex).AND.PRESENT(ChemName)) THEN
		!	IF (i .EQ. ChemIndex) THEN
		!		CALL WARN("Both a Chem Index ("//TRIM(INT2STR(ChemIndex))//") and a Chemical Name ("//  &
		!				  TRIM(ChemName)//") were passed to GetThermalVelocityOfAir.  They signify the "// &
		!				  "same chemical, which is acceptable, but this is bad form nonetheless.")
		!	ELSE
		!		CALL ERROR("Both a Chem Index ("//TRIM(INT2STR(ChemIndex))//") and a Chemical Name ("//TRIM(ChemName)// &
		!				   ") were passed to GetThermalVelocityOfAir.  They specify different chemicals.  Please fix.")
		!	END IF
		!END IF

		MolecularDiffusionCoefficient = 0.375/(4.5*Angstrom)**2./GetAirDensity(XP,YP,ZP)*		&
										SQRT(Rstar*GetTemp(XP,YP,ZP)*AirMolecMass*ChemScale/	&
										Avogadro/GasMolecularMass(i)*(GasMolecularMass(i)+		&
										Avogadro*AirMolecMass)/2/Pi)

		RETURN
	END FUNCTION MolecularDiffusionCoefficient
