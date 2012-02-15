!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! SetGrid.h
!! This file sets grids and their gridp point volumes.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							       	!!
!!							      			!!
!! Month  Year   Name              Description				        !!
!! 07     2006   Matt Alvarado     Began Update History				!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	 !!
!! 1. SUBROUTINE SetGrids ()					 !!
!! 2. SUBROUTINE SetGridPointVolumes ()			         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine sets the number of gridcells for each	!!
	!! of the grids in the Eulerian Section of the Array.  It	!!
	!! should eventually be updated to read from a deck.		!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetGrids ()

		IMPLICIT NONE

		!!! TEMPERATURE !!!
		XGridPoints(TempGridIndex) = 1
		YGridPoints(TempGridIndex) = 1
		ZGridPoints(TempGridIndex) = 1

		!!! PRESSURE !!!
		XGridPoints(PressGridIndex)	= 1
		YGridPoints(PressGridIndex)	= 1
		ZGridPoints(PressGridIndex)	= 1

		!!! GAS PHASE CHEMISTRY !!!
		XGridPoints(GasPhaseChemGridIndex) = 1
		YGridPoints(GasPhaseChemGridIndex) = 1
		ZGridPoints(GasPhaseChemGridIndex) = 1

		CALL SetGridPointVolumes

		!! Make a Logical Fn to see if everything
		!! is working on the same grid scale
		IF ((XGridPoints(TempGridIndex) .EQ. XGridPoints(PressGridIndex))        .AND. &
			(XGridPoints(TempGridIndex) .EQ. XGridPoints(GasPhaseChemGridIndex)) .AND. &
			(YGridPoints(TempGridIndex) .EQ. YGridPoints(PressGridIndex))        .AND. &
			(YGridPoints(TempGridIndex) .EQ. YGridPoints(GasPhaseChemGridIndex)) .AND. &
			(ZGridPoints(TempGridIndex) .EQ. ZGridPoints(PressGridIndex))        .AND. &
			(ZGridPoints(TempGridIndex) .EQ. ZGridPoints(GasPhaseChemGridIndex))) THEN

			AllSameGrid = .TRUE.
		ELSE
			AllSameGrid = .FALSE.
		END IF

	END SUBROUTINE SetGrids

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set Grid Point Volumes !!
	SUBROUTINE SetGridPointVolumes ()

		USE ModelParameters, ONLY : DomainX, DomainY, DomainZ
	
		IMPLICIT NONE

		GridPointVolume(TempGridIndex) = &
			DomainX*DomainY*DomainZ/XGridPoints(TempGridIndex)/YGridPoints(TempGridIndex)/ZGridPoints(TempGridIndex)
		GridPointVolume(PressGridIndex) = &
			DomainX*DomainY*DomainZ/XGridPoints(PressGridIndex)/YGridPoints(PressGridIndex)/ZGridPoints(PressGridIndex)
		GridPointVolume(GasPhaseChemGridIndex) = &
			DomainX*DomainY*DomainZ/XGridPoints(GasPhaseChemGridIndex)/YGridPoints(GasPhaseChemGridIndex)/ &
			ZGridPoints(GasPhaseChemGridIndex)

		RETURN
	END SUBROUTINE SetGridPointVolumes
