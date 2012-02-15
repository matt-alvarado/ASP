!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Water.h
!! This file contains the calculation and retrieval     
!! functions for contains the calculation and retrieval 
!! functions for the Saturation Vapor Pressure of Water in Air

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!! 07/13  2007   Matt Alvarado     Fixed Relative Humidity functions
!! 07/30  2007   Matt Alvarado     Fixed more relative humidity functions
!! 08/29  2007   Matt Alvarado     Fixed mixing ratio and air density functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. FUNCTION GetSatVapPress (XP, YP, ZP) RESULT (SatVapPress)
!! 2. FUNCTION GetSatVapBurden(XP, YP, ZP)
!! 3. FUNCTION GetSatVapConcentration (XP, YP, ZP)
!! 4. FUNCTION GetSatMixingRatio (XP, YP, ZP) RESULT (SatMixRatio)
!! 5. FUNCTION GetMixingRatio (XP, YP, ZP, RHin) RESULT (MixRatio)
!! 6. FUNCTION GetAirDensity (XP, YP, ZP) RESULT (Density)
!! 7. FUNCTION GetRelativeHumidity (XP, YP, ZP)
!! 8. FUNCTION GetRelativeHumidityFromGrid (X, Y, Z)
!! 9. FUNCTION GetRelativeHumidityFromBurden (XP, YP, ZP, WaterBurden)
!!10. SUBROUTINE SetAllRelativeHumidities (RH)
!!11. SUBROUTINE SetRelativeHumidity (XP,YP,ZP,RH)
!!12. FUNCTION SurfaceTensionOfWater (XP,YP,ZP)
!!13. FUNCTION DiffusionCoefficientOfWater (XP,YP,ZP)
!!14. FUNCTION GetCpm (XP,YP,ZP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Saturation Vapor Pressure Values.  At some point,
	!! this may be updated to take interpolated values.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!! This function checks out (07/13/07, MJA)
	REAL*8 FUNCTION GetSatVapPress (XP, YP, ZP) RESULT (SatVapPress)

		USE ModelParameters, ONLY : mbar

		implicit none

		real*8, intent(in) :: XP, YP, ZP

		real*8  :: Temp

		Temp  = GetTemp (XP,YP,ZP)
		SatVapPress = 6.112 * dexp(17.67*(Temp-273.15) /(Temp-29.65)) * mbar

	RETURN
	END FUNCTION GetSatVapPress

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Saturation Vapor Concentration Values.  At some point,
	!! this may be updated to take interpolated values.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetSatVapBurden(XP, YP, ZP)

		IMPLICIT NONE
	
		REAL*8, INTENT(IN) :: XP, YP, ZP

		GetSatVapBurden = GetSatVapConcentration (XP, YP, ZP)

		RETURN
	END FUNCTION GetSatVapBurden

	REAL*8 FUNCTION GetSatVapConcentration (XP, YP, ZP) 

		USE ModelParameters, ONLY : eps, ChemScale

		IMPLICIT NONE

		real*8, intent(in) :: XP, YP, ZP
		real*8 :: SatVapPress

		SatVapPress			   = GetSatVapPress (XP, YP, ZP)

		GetSatVapConcentration = ChemScale * SatVapPress/					&
								(GetPress (XP, YP, ZP)) * GetM(XP,YP,ZP)

	RETURN
	END FUNCTION GetSatVapConcentration

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Saturation Mass Mixing Ratio of water wv_sat
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetSatMixingRatio (XP, YP, ZP) RESULT (SatMixRatio)

		use ModelParameters, ONLY : eps, ppb
		implicit none

		real*8, intent(in) :: XP, YP, ZP
		real*8			   :: SatVapPress	! internal

		SatVapPress = GetSatVapPress (XP, YP, ZP)

		SatMixRatio = SatVapPress*eps/					&			! ppb; See Emanuel, 1994
					  (GetPress (XP, YP, ZP)- SatVapPress) * 1.0d9 * ppb

	RETURN
	END FUNCTION GetSatMixingRatio

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Mass Mixing Ratio  of Water Vapor, wv
	!! See 2.31 of Jacobson, 2005 (2nd ed.)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetMixingRatio (XP, YP, ZP, RHin) RESULT (MixRatio)

		use ModelParameters, ONLY : eps, ppb
		implicit none

		real*8, intent(in) :: XP, YP, ZP
		real*8, optional   :: RHin
		real*8			   :: SatVapPress, RH

		IF (Present(RHin)) THEN
			RH = RHin
		ELSE
			RH = GetRelativeHumidity (XP, YP, ZP)
		END IF

		SatVapPress = GetSatVapPress (XP, YP, ZP)
		MixRatio = RH * SatVapPress *eps / (GetPress (XP, YP, ZP) - RH * SatVapPress) ! Mixing ratio (-)
		!! removed a factor of 1e9 from this, 11/12/2002

	RETURN
	END FUNCTION GetMixingRatio

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Air Density								!!
	!!														!!
	!! Inputs include location and Relative Humidity (RH)	!!
	!! Output is Air Density in (g / lengthscale^3)			!!
	!! Tested against Example 2.3 in Jacobson				!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetAirDensity (XP, YP, ZP) RESULT (Density)

		use ModelParameters
		implicit none

		real*8, intent(in) :: XP, YP, ZP
		real*8			   :: MixRatio, RH

		MixRatio = GetMixingRatio (XP, YP, ZP)

		Density = GetPress (XP, YP, ZP) * (1. + MixRatio) /		&
				  (GetTemp (XP, YP, ZP) * Rdry * (1. + MixRatio/eps))

	RETURN
	END FUNCTION GetAirDensity


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Invert the Concentration of Water to get the local relative humidity !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetRelativeHumidity (XP, YP, ZP)

		USE ModelParameters, ONLY : eps, ChemScale

		IMPLICIT NONE

		REAL*8 :: XP, YP, ZP

		!! Local Variables
		REAL*8 :: WaterConc

		WaterConc = GetGasChem(XP,YP,ZP,1)

		GetRelativeHumidity = GetPress(XP,YP,ZP) * WaterConc / &
							  (GetSatVapPress(XP,YP,ZP) * (GetM(XP,YP,ZP)))


		RETURN
	END FUNCTION GetRelativeHumidity

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Same as GetRelativeHumidity, but at a particular grid point instead  !!
	!! of location.															!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetRelativeHumidityFromGrid (X, Y, Z)

		USE InfrastructuralCode, ONLY : ERROR
		USE ModelParameters, ONLY : eps, ChemScale

		IMPLICIT NONE 

		INTEGER :: X,Y,Z

		REAL*8 :: XX(3)

		!! Local Variables
		REAL*8 :: T, WaterConc

		IF (.NOT.AllSameGrid) &
		CALL ERROR("The subroutine GetRelativeHumidityFromGrid() was run with a grid structure that does not space "// &
				   "grid-cells the same for temperature, chemistry, pressure, and aerosol.  You can't do this with "// &
				   "this function, use GetTemp() instead.")

		XX = LocateGridPoint(X, Y, Z, 1) 

		WaterConc = GetGasChem(XX(1),XX(2),XX(3),1)

		GetRelativeHumidityFromGrid = GetPress(XX(1),XX(2),XX(3)) * WaterConc / &
							  (GetSatVapPress(XX(1),XX(2),XX(3)) * GetM(XX(1),XX(2),XX(3)))


		RETURN
	END FUNCTION GetRelativeHumidityFromGrid


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! There are times in the condensation routine that we are holding a		!!
	!! water burden for a particular gridcell in a vector that is running		!!
	!! through LSODES.  We want to be able to get a relative humidity value		!!
	!! from this without actually altering the Eulerian Array.  So this takes	!!
	!! that burden, diagnoses the other environmental parameters, and then		!!
	!! delivers a relative humidity value.										!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetRelativeHumidityFromBurden (XP, YP, ZP, WaterBurden)

		USE ModelParameters, ONLY : eps, ChemScale
		REAL*8 :: XP, YP, ZP

		!! Local Variables
		REAL*8 :: T, WaterBurden
		
		GetRelativeHumidityFromBurden =																							&
				GetPress(XP,YP,ZP) * WaterBurden  /							&
				(GetSatVapPress(XP,YP,ZP) * (GetM(XP,YP,ZP)))

		RETURN
	END FUNCTION GetRelativeHumidityFromBurden

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Reset all of the relative humidities at once !!
!! primarily for debugging but also for simple  !!
!! experiments.									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetAllRelativeHumidities (RH)

	IMPLICIT NONE

	REAL*8  :: RH
	INTEGER :: I, J, K

	DO I = 1, XGridPoints(GasPhaseChemGridIndex)
	DO J = 1, YGridPoints(GasPhaseChemGridIndex)
	DO K = 1, ZGridPoints(GasPhaseChemGridIndex)
		CALL SetRelativeHumidity (I,J,K,RH)
	END DO; END DO; END DO

END SUBROUTINE SetAllRelativeHumidities 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Reset the relative humidity at a given  !!
!! gridcell location to a particular value !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetRelativeHumidity (XP,YP,ZP,RH)

	USE ModelParameters, ONLY : eps, ChemScale, ppbscale

	IMPLICIT NONE

	INTEGER :: XP, YP, ZP
	REAL*8  :: RH, XX(3)

	!! Local Variables
	REAL*8 :: T, Press, WaterMR, SPress

	!! translate the grid cell number into a 3-space location
	XX     = LocateGridPoint(XP, YP, ZP, GasPhaseChemGridIndex)

	SPress = GetSatVapPress(XX(1),XX(2),XX(3))

	!! Maybe put a MM on the top?
	GridGasChem(XP,YP,ZP,1) = ChemScale * SPress * RH * GetM(XX(1),XX(2),XX(3)) / (GetPress(XX(1),XX(2),XX(3)))

	RETURN

END SUBROUTINE SetRelativeHumidity

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This approximation for SURFACE TENSION of Water (dyn / cm) is	!!
	!! from Jacobson (15.2)  from Prupp & Klett fit of a polynomian		!!
	!! expression to Dorsch and Hacker (1951) see books for ref.		!!
	!! It is only really good for [-40,40] but use it everywhere here.	!!
	!! Further dependence of surface tension on electrolytic properties !!
	!! is treated by SurfaceTension() in ParticleAttributes.h		    !! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION SurfaceTensionOfWater (XP,YP,ZP)

		USE InfrastructuralCode, ONLY : Warn, REAL2STR
		USE ModelParameters,     ONLY : dyncm

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: XP,YP,ZP

		!! Local Variables
		REAL*8 :: Tc

		Tc = GetTemp(XP,YP,ZP) - 273.15

		!! Supercooled or not Supercooled
		IF (Tc < 0) THEN
			SurfaceTensionOfWater = 75.93+0.115*Tc+0.06818*(Tc**2)+6.511d-3*(Tc**3)+		&
					   2.933d-4*(Tc**4)+6.283d-6*(Tc**5)+5.285d-8*(Tc**6)
		ELSE
			SurfaceTensionOfWater = 76.1 - 0.155*Tc
		END IF

		SurfaceTensionOfWater = SurfaceTensionOfWater * dyncm

		RETURN
	END FUNCTION SurfaceTensionOfWater

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! DIFFUSIVITY OF WATER VAPOR as a functino of temperature and pressure. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION DiffusionCoefficientOfWater (XP,YP,ZP)

		USE ModelParameters, ONLY : Atmospheres, cm

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: XP,YP,ZP

		!! Calculate the Diffusion Coefficient of Water (cm^2/s)
		!! See Seinfeld and Pandis p.801 or Prupp & Klett
		DiffusionCoefficientOfWater = (0.211*Atmospheres/GetPress(XP,YP,ZP))*(GetTemp(XP,YP,ZP)/273.)**1.94 * cm * cm

		RETURN
	END FUNCTION DiffusionCoefficientOfWater

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SPECIFIC HEAT OF MOIST AIR AT CONSTANT PRESSURE !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetCpm (XP,YP,ZP)

		USE Modelparameters, ONLY : Cpd, CpV

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: XP,YP,ZP

		!! Internal Variables
		REAL*8 :: MR

		MR = GetMixingRatio(XP,YP,ZP)

		GetCpm = (Cpd + MR * CpV) / (1. + MR)

		RETURN
	END FUNCTION GetCpm