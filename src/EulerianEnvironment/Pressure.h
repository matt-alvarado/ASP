!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Pressure.h
!! Contains subroutines for inputting and outputting pressure data,
!! changing temperature and pressure, calcualing M (total gas) 
!! concentration, and expanding grid point volumes.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. SUBROUTINE SetPressField (XP, YP, ZP, Press)
!! 2. FUNCTION   GetPress (XP, YP, ZP) RESULT (Press)
!! 3. FUNCTION   GetPressFromGrid (XP, YP, ZP) RESULT (Press)
!! 4. FUNCTION   GetM (XP, YP, ZP) RESULT (M)
!! 5. SUBROUTINE ChangeTempPress (NewTemperature, NewPressure)
!! 6. SUBROUTINE ExpandGPVolume (VolumeScale)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set the initial PRESSURE field.  This will be ammended or	  !!
	!! forced eventually to reflect an actual environment.			  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetPressField (XP, YP, ZP, Press)

		USE ModelParameters, ONLY : PressureScale

		implicit none

		!! Inputs Variables:
		integer :: XP, YP, ZP
		real*8  :: Press

		GridPress(XP, YP, ZP) = Press * PressureScale

		RETURN
	END SUBROUTINE SetPressField


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the PRESSURE	!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetPress (XP, YP, ZP) RESULT (Press)

		implicit none

		real*8, intent(in)	:: XP, YP, ZP	! Input

		!! Call the interpolating function to land the correct
		!! temperature from the gridded array.
		Press = InterpolateGrid (XP, YP, ZP, GridPress, PressGridIndex)

		return
	END FUNCTION GetPress

	REAL*8 FUNCTION GetPressFromGrid (XP, YP, ZP) RESULT (Press)

		IMPLICIT NONE

		INTEGER	:: XP, YP, ZP	! Input

		Press = GridPress(XP, YP, ZP)

		return
	END FUNCTION GetPressFromGrid


	!!!!!!!!!!!!!!!!
	!! Retrieve M !!
	!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetM (XP, YP, ZP) RESULT (M)
		
		USE ModelParameters, ONLY : PressureScale

		implicit none

		real*8, intent(in)	:: XP, YP, ZP	! Input

		!! Use the Ideal Gas Law to scale a reference M to local values
		M =	2.688e19 * (GetPress(XP,YP,ZP)/1013./PressureScale)/(GetTemp (XP,YP,ZP)/273.)

		return
	END FUNCTION GetM

