!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Initialize.h
!! This file sets the temperature and pressure grids
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:			!!
!! 1. SUBROUTINE InitializeGridPointFields (T, Press)					!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE InitializeGridPointFields (T, Press)

		USE InfrastructuralCode, ONLY : Error
		implicit none

		!! Define the input variables
		real*8  :: T, Press

		!! Define variables for internal use
		integer   :: i, j, k, allocation_error

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Initialize Primary Variables First!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		CALL SetGrids()		! This sets the extent of each of the grid meshes

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Initialize TEMPERATURE !!
		allocation_error = 0
		allocate (GridTemp(XGridPoints(TempGridIndex),YGridPoints(TempGridIndex),ZGridPoints(TempGridIndex)), &
				  stat = allocation_error)

				! Allocation Error !
				if (allocation_error > 0) CALL ERROR ("Allocation of GridTemp could not proceed in InitializeGridPointFields()")

		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Initialize PRESSURE !!
		allocation_error = 0
		allocate (GridPress(XGridPoints(PressGridIndex),YGridPoints(PressGridIndex),ZGridPoints(PressGridIndex)), &
				  stat = allocation_error)

				! Allocation Error !
				if (allocation_error > 0) CALL ERROR ("Allocation of GridPress could not proceed in InitializeGridPointFields()")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! POPULATE the fields allocated above with initial values	!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (1 == 1) then

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Populate the TEMPERATURE grid !!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			DO i = 1,XGridPoints(TempGridIndex)
				DO j = 1,YGridPoints(TempGridIndex)
					DO k = 1,ZGridPoints(TempGridIndex)
						CALL SetTempField (i, j, k, T)
					END DO
				END DO
			END DO


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Populate the PRESSURE grid !!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			DO i = 1,XGridPoints(PressGridIndex)
				DO j = 1,YGridPoints(PressGridIndex)
					DO k = 1,ZGridPoints(PressGridIndex)
						CALL SetPressField (i, j, k, Press)
					END DO
				END DO
			END DO

			AllSamePressure = .TRUE.


		END IF

		RETURN
	END SUBROUTINE InitializeGridPointFields