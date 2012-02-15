!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Temperature.h
!! This file contains the calculation and retrieval     
!! functions for various thermodynamics fields of the gridded
!! environment.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. SUBROUTINE SetTempField (XP, YP, ZP, Temp)
!! 2. FUNCTION GetThermalVelocityOfAir (XP, YP, ZP)
!! 3. FUNCTION MolecularThermalDiffusivity (XP, YP, ZP) 
!! 4. FUNCTION MolecularThermalConductivity (XP, YP, ZP)
!! 5. FUNCTION GetDynamicViscosityOfAir (XP, YP, ZP) RESULT (DynVisc)
!! 6. FUNCTION GetTemp (XP, YP, ZP) RESULT (Temp)
!! 7. FUNCTION GetTempFromGrid (XP, YP, ZP) RESULT (Temp)
!! 8. FUNCTION MeanFreePathOfAir (XP, YP, ZP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set the initial TEMPERATURE field.  This will be ammended or	  !!
	!! forced eventually to reflect an actual environment.			  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetTempField (XP, YP, ZP, Temp)

		implicit none

		!! Inputs Variables:
		integer :: XP, YP, ZP
		real*8  :: Temp

		GridTemp(XP, YP, ZP) = Temp

		RETURN
	END SUBROUTINE SetTempField

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The THERMAL VELOCITY of AIR is simply proportional to the sq root	!!
	!! of the temperature.  Derivation is from gas kinetic theory, setting	!!
	!! average kinetic energy of air molec to 4/pi kB T.					!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetThermalVelocityOfAir (XP, YP, ZP)

		use ModelParameters, ONLY : ThmVelFac, cm, AirMolecMass

		implicit none

		real*8, intent(in) :: XP, YP, ZP ! inputs

		GetThermalVelocityOfAir = ThmVelFac * SQRT(GetTemp(XP, YP, ZP) / AirMolecMass )

		RETURN
	END FUNCTION GetThermalVelocityOfAir 

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The MOLECULAR THERMAL DIFFUSIVITY, calculated from List (1984)'s !!
	!! Smithsonian Meteorological Tables fit.							!!
	!!																	!!
	!! Dh = Kappa_d / (rho cpd)											!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION MolecularThermalDiffusivity (XP, YP, ZP) 

		USE ModelParameters, ONLY : m, joules, cpd

		IMPLICIT NONE

		REAL*8 :: XP,YP,ZP

		!! The molecular thermal diffusivity (J cm-1 s-1 K-1)
		MolecularThermalDiffusivity = MolecularThermalConductivity (XP, YP, ZP) / GetAirDensity(XP,YP,ZP) / cpd

	END FUNCTION MolecularThermalDiffusivity 

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The MOLECULAR THERMAL CONDUCTIVITY, calculated from List (1984)'s!!
	!! Smithsonian Meteorological Tables fit.							!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION MolecularThermalConductivity (XP, YP, ZP) 

		USE ModelParameters, ONLY : joules, cm

		IMPLICIT NONE

		REAL*8 :: XP,YP,ZP

		!! List 1984's empirical fit to Kd (J cm-1 s-1 K-1):
		MolecularThermalConductivity = ( 2.3807d-4 + 7.1128d-7 * (GetTemp(XP,YP,ZP)-273.16) )* joules / cm

	END FUNCTION MolecularThermalConductivity 


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The DYNAMIC VISCOSITY of AIR is calculated using Sutherland's Equation,!! 
	!! taken from List, 1984 (Smithsonian Met Tables) via Jacobson's Modeling !!
	!! book.  The quantity itself is the ratio of shearing stress to shear.	  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetDynamicViscosityOfAir (XP, YP, ZP) RESULT (DynVisc)

		use ModelParameters
		implicit none

		real*8, intent(in) :: XP, YP, ZP ! inputs
		real*8			   :: T			 ! internal

		T = GetTemp (XP, YP, ZP)

		DynVisc = 1.4962864d-5 * T**1.5 / (T + 120.) * grams / cm

	END FUNCTION GetDynamicViscosityOfAir 


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the TEMPERATURE										  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetTemp (XP, YP, ZP) RESULT (Temp)

		implicit none

		real*8, intent(in) :: XP, YP, ZP

		!! Call the interpolating function to land the correct
		!! temperature from the gridded array.
		Temp = InterpolateGrid (XP, YP, ZP, GridTemp, TempGridIndex)

		return
	END FUNCTION GetTemp
 
	!! Same as GetTemp but not via interpolation
 	REAL*8 FUNCTION GetTempFromGrid (XP, YP, ZP) RESULT (Temp)
		
		USE InfrastructuralCode, ONLY : ERROR

		IMPLICIT NONE

		INTEGER :: XP, YP, ZP

		IF (.NOT.AllSameGrid) &
		CALL ERROR("The subroutine GetTempFromGrid() was run with a grid structure that does not space grid-cells the same "// &
				   "for temperature, chemistry, pressure, and aerosol.  You can't do this with this function, use GetTemp() "// &
				   "instead.")
	
		Temp = GridTemp(XP, YP, ZP)

		RETURN
	END FUNCTION GetTempFromGrid



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Find the Mean Free Path at a Particular Location !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION MeanFreePathOfAir (XP, YP, ZP)

	USE ModelParameters, ONLY : grams, cm, lengthscale, massscale

	IMPLICIT NONE
	
	!! External
	REAL*8 :: XP,YP,ZP

	!! Internal
	REAL*8 :: T 

	T = GetTemp(XP,YP,ZP)

	!! The mean free path is equal to 2 * the dynamic viscosity of air
	!! divided by (density of air * the thermal velocity of air)
	MeanFreePathOfAir = 2*GetDynamicViscosityOfAir(XP,YP,ZP)/ &
						(GetThermalVelocityOfAir(XP,YP,ZP) *	&
						GetAirDensity(XP,YP,ZP)*lengthscale)

	!! This result does include an implicit factor of lengthscale!

	RETURN
	END FUNCTION MeanFreePathOfAir

