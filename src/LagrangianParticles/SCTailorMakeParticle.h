!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Donnan Steele -- MELAM -- (c) 2001   !!
!! Update Information:					!!  
!! 08     2005   Matt Alvarado	   Added comments
!! 07/24  2006   Matt Alvarado     Change "Pointer => NULL()" to				!!
!!									NULLIFY(POINTER) to fit pgf90							!!
!!				
!!										!!			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!	
!! 1. FUNCTION SCTailorMakeParticle(WhichPopulation, Radius)        !!
!! 2. FUNCTION SCTailorMakeWaterParticle(Radius)					!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Create a particle with given dry rad !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION SCTailorMakeParticle(WhichPopulation, Radius)

	USE GridPointFields,	 ONLY : ClosestGridPoint, GasPhaseChemGridIndex, GetTemp
	USE ModelParameters,	 ONLY : LagrSectSizeCutoff,			&
									lengthscale,				&
									Pi, EulerMass,				&
									micron,	moles,				&
									AverageAerosolDensity,		&
									DomainX, DomainY, DomainZ
	USE Chemistry,			 ONLY : HowManyAqChems,				&
									HowManyAqCations,			&
									HowManyAqAnions,			&
									AqMolecularMass,			&
									HowManyAqEqReactions
	USE InfrastructuralCode, ONLY : REAL2STR, Transcript, Warn, Error

	IMPLICIT NONE

	!! Input Variables
	INTEGER :: WhichPopulation
	REAL*8  :: Radius

	!! Local Variables
	REAL*8 :: II, JJ, KK
	TYPE(Particle), POINTER :: SCTailorMakeParticle, CurrentParticle
	INTEGER				    :: allocation_error, I,J,X,Y,Z
	LOGICAL, PARAMETER		:: Scaffolding = .FALSE. !.TRUE.

	IF (SCAFFOLDING) CALL Transcript(">>>Entering MakeParticle()<<<")
	IF (SCAFFOLDING) CALL Transcript("")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! For sampling, we want to use a random number to find a number of standard	!!
	!! deviations away from the mean value.  To do this, we must invert the error	!!
	!! function.  To do this, we use the following numerical trick, which inverts	!!
	!! an error function for a series with a mean at the Euler-Masscheroni Constant	!!
	!! and a standard deviation of pi/sqrt(6):										!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	II = Radius

	!! If particle is large enough to be tracked, 
	!! create the New Particle and allocate its arrays
	ALLOCATE (SCTailorMakeParticle, stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of SCTailorMakeParticle Failed in MakeParticle()")
	ALLOCATE (SCTailorMakeParticle%AqChems(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of SCTailorMakeParticle%AqChems Failed in MakeParticle()")

	!! Allocate the Appropriate Electrolytic Arrays
	ALLOCATE (SCTailorMakeParticle%GammaMixed(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of SCTailorMakeParticle%GammaMixed Failed in MakeParticle()")

	!! Initialize
	DO I = 1, HowManyAqEqReactions
		SCTailorMakeParticle%GammaMixed(I) = 0
	END DO

	SCTailorMakeParticle%IonicStr		= 1.
	SCTailorMakeParticle%EmbryoRadius	= II
	SCTailorMakeParticle%InsolubleRadius	= 0.
	SCTailorMakeParticle%EffectiveRadius	= II

	!! Create a Location for the Particle
	CALL Random_Number(II)
	SCTailorMakeParticle%location(1) = II*DomainX
	CALL Random_Number(II)
	SCTailorMakeParticle%location(2) = II*DomainY
	CALL Random_Number(II)
	SCTailorMakeParticle%location(3) = II*DomainZ

	!! Place the particle at the local gridpoint
	CALL ClosestGridPoint(SCTailorMakeParticle%location(1), SCTailorMakeParticle%location(2), SCTailorMakeParticle%location(3), &
						  GasPhaseChemGridIndex, X, Y, Z) 

	!! For Now, the Initial Temperature of the Particle is that of the 
	!! Local Environment (1/2002)
	SCTailorMakeParticle%Temperature = GetTemp(SCTailorMakeParticle%Location(1),SCTailorMakeParticle%Location(2),	&
											   SCTailorMakeParticle%Location(3))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Create a Velocity for the Particle.  It is some random		!!
	!! percentage of the terminal velocity, randomly apportioned	!!
	!! between the three directions.								!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL Random_Number(II)
	CALL Random_Number(JJ)
	CALL Random_Number(KK)

	!! Set initial velocity at half or more of terminal velocity
	KK = (1.+KK)/2.  ! * TerminalVelocity(SCTailorMakeParticle%Radius)

	!! Set initial angles
	II = 2.*Pi*II	 ! Flat Angle is out of 360 degrees
	JJ = Pi*JJ-Pi/2. ! Azimuthal angle is out of 180 degrees

	!! Change this into two random angles and then move to cosines and sines from there:
	SCTailorMakeParticle%velocity(1) = KK*DCOS(II)*DCOS(JJ)
	SCTailorMakeParticle%velocity(2) = KK*DSIN(II)*DCOS(JJ)
	SCTailorMakeParticle%velocity(3) = KK*DSIN(JJ)

	IF (SCAFFOLDING) CALL Transcript("Effective Radius          : ("//TRIM(REAL2STR(SCTailorMakeParticle%EffectiveRadius)))
	IF (SCAFFOLDING) CALL Transcript("Velocity (X,Y,Z): ("//TRIM(REAL2STR(SCTailorMakeParticle%velocity(1)))//","//		&
									 TRIM(REAL2STR(SCTailorMakeParticle%velocity(2)))//","//							&
									 TRIM(REAL2STR(SCTailorMakeParticle%velocity(3)))//")")
	
	!! Load in the chemistry assignments, stored as Moles.
	!! Calculate a Mass of each chemical type (using AerosolModes, which is simply a mass fraction value)
	!! Assume, for this calculation, that the mean desity is AverageAerosolDensity (=1).

	!! assign a token amount of water so it doesn't freak out!
	SCTailorMakeParticle%AqChems(1) = 100 *	&
							 4./3.*Pi*(SCTailorMakeParticle%EffectiveRadius)**3.		&
							 *AverageAerosolDensity/				&
							 AqMolecularMass(1) !*moles  !! commenting this out because get factor squared when varying
	DO I = 2, HowManyAqChems
		SCTailorMakeParticle%AqChems(I) = AerosolModes(WhichPopulation,3+I) *	&
								 4./3.*Pi*(SCTailorMakeParticle%EffectiveRadius)**3.		&
								 *AverageAerosolDensity/				&
								 AqMolecularMass(I) !*moles  !! commenting this out because get factor squared when varying
	END DO

	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		SCTailorMakeParticle%AqChems(I) = 0.
	END DO


	IF (SCAFFOLDING) CALL Transcript(">>>Exiting MakeParticle()<<<")
	
	RETURN

END FUNCTION SCTailorMakeParticle


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Create a particle with given wet rad !!
!! which will then barf if do therm stf !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION SCTailorMakeWaterParticle(Radius)

	USE GridPointFields,	 ONLY : ClosestGridPoint, GasPhaseChemGridIndex, GetTemp
	USE ModelParameters,	 ONLY : LagrSectSizeCutoff,			&
									lengthscale,				&
									Pi, EulerMass,				&
									micron,	moles,				&
									AverageAerosolDensity,		&
									DomainX, DomainY, DomainZ
	USE Chemistry,			 ONLY : HowManyAqChems,				&
									HowManyAqCations,			&
									HowManyAqAnions,			&
									AqMolecularMass,			&
									HowManyAqEqReactions
	USE InfrastructuralCode, ONLY : REAL2STR, Transcript, Warn, Error

	IMPLICIT NONE

	!! Input Variables
	INTEGER :: WhichPopulation
	REAL*8  :: Radius

	!! Local Variables
	REAL*8 :: II, JJ, KK
	TYPE(Particle), POINTER :: SCTailorMakeWaterParticle, CurrentParticle
	INTEGER				    :: allocation_error, I,J,X,Y,Z
	LOGICAL, PARAMETER		:: Scaffolding = .FALSE. !.TRUE.

	IF (SCAFFOLDING) CALL Transcript(">>>Entering MakeParticle()<<<")
	IF (SCAFFOLDING) CALL Transcript("")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! For sampling, we want to use a random number to find a number of standard	!!
	!! deviations away from the mean value.  To do this, we must invert the error	!!
	!! function.  To do this, we use the following numerical trick, which inverts	!!
	!! an error function for a series with a mean at the Euler-Masscheroni Constant	!!
	!! and a standard deviation of pi/sqrt(6):										!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	II = Radius

	!! If particle is large enough to be tracked, 
	!! create the New Particle and allocate its arrays
	ALLOCATE (SCTailorMakeWaterParticle, stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of SCTailorMakeParticle Failed in MakeParticle()")
	ALLOCATE (SCTailorMakeWaterParticle%AqChems(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of SCTailorMakeParticle%AqChems Failed in MakeParticle()")

	!! Allocate the Appropriate Electrolytic Arrays
	ALLOCATE (SCTailorMakeWaterParticle%GammaMixed(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of SCTailorMakeParticle%GammaMixed Failed in MakeParticle()")

	!! Initialize
	DO I = 1, HowManyAqEqReactions
		SCTailorMakeWaterParticle%GammaMixed(I) = 0
	END DO

	SCTailorMakeWaterParticle%IonicStr		= 0.
	SCTailorMakeWaterParticle%EmbryoRadius		= II
	SCTailorMakeWaterParticle%InsolubleRadius	= 0.
	SCTailorMakeWaterParticle%EffectiveRadius	= II


	SCTailorMakeWaterParticle%SolutionDensity		= 1.
	SCTailorMakeWaterParticle%ParticleDensity		= 1.


	SCTailorMakeWaterParticle%Sectional = .FALSE.
	SCTailorMakeWaterParticle%NumberOfParticles	= 1

	!! Create a Location for the Particle
	SCTailorMakeWaterParticle%location(1) = 0.00001
	SCTailorMakeWaterParticle%location(2) = 0.00001
	SCTailorMakeWaterParticle%location(3) = 0.00001

	!! Place the particle at the local gridpoint
	CALL ClosestGridPoint(SCTailorMakeWaterParticle%location(1), SCTailorMakeWaterParticle%location(2), &
						  SCTailorMakeWaterParticle%location(3), GasPhaseChemGridIndex, X, Y, Z) 

	!! For Now, the Initial Temperature of the Particle is that of the 
	!! Local Environment (1/2002)
	SCTailorMakeWaterParticle%Temperature = GetTemp(SCTailorMakeWaterParticle%Location(1),SCTailorMakeWaterParticle%Location(2), &
													SCTailorMakeWaterParticle%Location(3))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Create a Velocity for the Particle.  It is some random		!!
	!! percentage of the terminal velocity, randomly apportioned	!!
	!! between the three directions.								!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL Random_Number(II)
	CALL Random_Number(JJ)
	CALL Random_Number(KK)

	!! Set initial velocity at half or more of terminal velocity
	KK = (1.+KK)/2.  ! * TerminalVelocity(SCTailorMakeParticle%Radius)

	!! Set initial angles
	II = 2.*Pi*II	 ! Flat Angle is out of 360 degrees
	JJ = Pi*JJ-Pi/2. ! Azimuthal angle is out of 180 degrees

	!! Change this into two random angles and then move to cosines and sines from there:
	SCTailorMakeWaterParticle%velocity(1) = KK*DCOS(II)*DCOS(JJ)
	SCTailorMakeWaterParticle%velocity(2) = KK*DSIN(II)*DCOS(JJ)
	SCTailorMakeWaterParticle%velocity(3) = KK*DSIN(JJ)

	IF (SCAFFOLDING) CALL Transcript("Effective Radius          : ("//TRIM(REAL2STR(SCTailorMakeWaterParticle%EffectiveRadius)))
	IF (SCAFFOLDING) CALL Transcript("Velocity (X,Y,Z): ("//TRIM(REAL2STR(SCTailorMakeWaterParticle%velocity(1)))//","//   &
									  TRIM(REAL2STR(SCTailorMakeWaterParticle%velocity(2)))//","//						   &
									  TRIM(REAL2STR(SCTailorMakeWaterParticle%velocity(3)))//")")
	

	!! assign a token amount of water so it doesn't freak out!
	SCTailorMakeWaterParticle%AqChems(1) = 4./3.*Pi*(SCTailorMakeWaterParticle%EffectiveRadius)**3./AqMolecularMass(1)


	SCTailorMakeWaterParticle%AqChems(2) = 2./moles

	DO I = 2, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		SCTailorMakeWaterParticle%AqChems(I) = 0.
	END DO

	SCTailorMakeWaterParticle%EffectiveRadius  = SCTailorMakeWaterParticle%EffectiveRadius

	NULLIFY(SCTailorMakeWaterParticle%Next)
	!SCTailorMakeWaterParticle%Next => Null()

	RETURN
END FUNCTION SCTailorMakeWaterParticle


