!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! AqOrgCondensationIntegrator.h
!! This is the numerical integration routine for the flux-limited
!! kinetic dissolution of organic compounds into the aerosol aqueous phase.
!! As the main routine calls LSODES, the required ODE and Jacobian 
!! functions are given here as well
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!! 08/30  2006   Matt Alvarado     1. Removed water eq from						!!
!!									StepAqOrgCondensation, added threshold check!!
!!									before doing internal org eq				!!
!!								   2. Created APDNonDiss							!!									!!
!! 10/17  2006   Matt Alvarado     Moved EmpyGridCell check and aerosol sorting to 
!!									StepCondensationAll()
!! 07/17  2007   Matt Alvarado     Changed zero criteria to cycle if 
!!									less than 1 molecule/cm3 of compound
!!									Added between-phase org equilibrium to end of 
!!									aqorgintegration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE StepAqOrgCondensation ()
!! 2. SUBROUTINE AqOrgCondensationODEEvaluator()
!! 3. SUBROUTINE AqOrgCondensationJacobianEvaluator()
!! 4. SUBROUTINE AqOrgDissolutionFactors
!! 5. FUNCTION AqOrgCondensationRateCoefficients (XP, YP, ZP, ReactionNumb, NumbLinks)					!!
!! 6. SUBROUTINE APDNonDiss(Y, Timestep, neq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Call LSODES or the equilibrium routine to integrate forward	 !!
!! the concentration of aerosol / gas-phase water at the grid-	 !!
!! cell level.  It must construct a joint array of all of the	 !!
!! aerosol plus the grid cell for each of the condensing species !!
!! XP, YP, and ZP are the gridpoint numbers in gas grid space.	 !!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE StepAqOrgCondensation (TimestepSize, NumTimeSteps, XP, YP, ZP)

	USE InfrastructuralCode, ONLY : INT2STR, ERROR,					&
									Warn, Transcript,				&
									REAL2STR

	USE GridPointFields,	 ONLY : LocateGridPoint,				&
									getM,							&
									GasPhaseChemGridIndex,			&
									GetGridCellChemBurden,			&
									ReplaceGridCellChemBurden,		&
									GetTemp,						&
									GetRelativeHumidity,			&
									AddToGridCellChemBurden

	USE Aerosols,            ONLY : Particles, Particle,			&
									ParticleArray,					&
									RecalculateRadius,				&
									EqCoeff,						&
									Molality,						&
									CalculateIonicStrength,			&
									KusikMeissner,					&
									SortAerosolAtGridPointForCondensation,			&
									EquilibriumWaterContent, &
									FindAqElectrolyteEquilibriumForgridPoint

	USE Chemistry,			 ONLY : HowManyGasChems,				&
									HowManyEvolveGasChems,			&
									HowManyEvolveAqChems,			&
									GasPhaseChemicalNames,			&
									HowManyCondensingChems,			&
									AqEquilibriaList,				&
									HowManyAqChems,					&
									HowManyAqCations, &
									FindChem

	USE ModelParameters

	IMPLICIT NONE

	!! External Variables
	REAL*8  :: TimeStepSize		!! Should be in Seconds
	INTEGER :: NumTimeSteps, XP, YP, ZP

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!! LSODES Related Internal Variables !!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8,  ALLOCATABLE	:: y(:), rwork(:)
	REAL*8					:: t,outtime,rtol,atol
	INTEGER, ALLOCATABLE	:: iwork(:)
	INTEGER					:: itol,itask,istate,iopt,mf,neq,LIW,LRW
		
	!! Non-LSODES Local Variables
	INTEGER				:: I, J, K,allocation_error, Rxn, HowManyLinks, YLen,					&
						   HowManyNonZeroJacobianTerms, GasChemIndex1, TimeStep,	&
						   ReturnType, WhichDominantIon, WhichSecondIon , Flag
	REAL*8				:: II, JJ, KK, ChemTransferArray(HowManyGasChems), Temperature,	&
						   RH, XX(3), NH3Index, Store
	REAL*8, ALLOCATABLE :: InitialConcsArray(:), ChemStorageArray(:,:)
	LOGICAL :: EmptyGridCell
	REAL*8  :: BeginTime,EndTime

	LOGICAL, PARAMETER :: Scaffolding = .FALSE., OpenSystem = .FALSE., ForceEquilibrium = .TRUE.


	TYPE(Particle), POINTER :: FirstParticle, CurrentParticle

	!! SCSCSCSC
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript(">>>Entering StepCondensation<<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Find GRID POINT Values for NON-AEROSOL !!
	!! DEPENDENT Parameters                   !!
	XX          = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
	Temperature = GetTemp(XX(1),XX(2),XX(3))
	RH		    = GetRelativeHumidity(XX(1),XX(2),XX(3))
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! IF EQUILIBRIUM RATHER THAN FLUX, JUST RETURN !!
	!! (Equilibrium of all components is handled in StepCondensation)
	IF (.NOT. DissolutionEquilibriumOrFlux .AND. RH .LE. AerosolWaterEquilibriumRHThreshold) THEN			
		RETURN
	END IF

	IF (Scaffolding) CALL TRANSCRIPT("Using LSODES to integrate organic dissolution into aqueous phase.")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Everything should now be sorted           !!
	!! Back up to the head of the aerosol list   !!
	FirstParticle   => Particles(XP,YP,ZP)%First
	CurrentParticle => FirstParticle
		
	!! COUNT the number of LINKS in the gridcell that we should condense upon
	HowManyLinks  =  1
15	IF(ASSOCIATED(CurrentParticle%Next)) THEN
		CurrentParticle => CurrentParticle%Next
		IF (CurrentParticle%NumberOfParticles .GT. 0 .AND. .NOT.CurrentParticle%Dry) THEN
			HowManyLinks  = HowManyLinks + 1
			GOTO 15
		END IF
	ENDIF
	CurrentParticle => FirstParticle


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set the lengths of the various !!
	!! arrays and allocate them !!!!!!!!
	YLen	= 12 + HowManyLinks*8	  !! As per the definition below
	NEQ		= 1  + 3*HowManyLinks	  !! Maximum value is for reaction type 5.  Allocate based on this.

	!! Temporarily set this to the maximum value (which is also for reaction type 4)
	HowManyNonZeroJacobianTerms = 11*HowManyLinks + 1
	LRW  = 60 +  25*NEQ + 3*HowManyNonZeroJacobianTerms
	LIW  = 31 + neq + HowManyNonZeroJacobianTerms

	!! Allocate the LSODES working arrays
	ALLOCATE (Y(ylen), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of Y could not proceed in StepCondensation()")

	ALLOCATE (rwork(LRW), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of RWORK could not proceed in StepCondensation()")

	ALLOCATE (iwork(LIW), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of IWORK could not proceed in StepCondensation()")

	ALLOCATE (InitialConcsArray(NEQ), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of InitialConcsArray could not proceed in StepCondensation()")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This is indexed as follows:												!!
	!! ChemStorageArray(Rxn, 1)         = Change in Gas Phase Conc 1			!!
	!! ChemStorageArray(Rxn, 2)         = Change in Aerosol Conc 1 (1)			!!
	!! ...																		!!
	!! ChemStorageArray(Rxn, HML+1)     = Change in Aerosol Conc 1 (HowManyLinks)!
	!! ChemStorageArray(Rxn, HML+2)     = Change in Aerosol Conc 2 (1)			!!
	!! ...																		!!
	!! ChemStorageArray(Rxn, 2*HML+1)   = Change in Aerosol Conc 2 (HowManyLinks)!
	!! ChemStorageArray(Rxn, 2*HML+2)   = Change in Aerosol Conc 3 (1)			!!
	!! ...																		!!
	!! ChemStorageArray(Rxn, 3*HML+1)   = Change in Aerosol Conc 3 (HowManyLinks)!
	ALLOCATE (ChemStorageArray(HowManyAqOrganicDissolutionReactions,NEQ), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of ChemStorageArray could not proceed in StepCondensation()")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Initialize the LSODES Work Vectors !!
	DO I = 1, LRW !!!!!!!!!!!!!!!!!!!!!!!!!!
		RWORK(I) = 0.
	END DO
	DO I = 1, LIW
		IWORK(I) = 0
	END DO

	!! Initialize the Y vector
	DO J = 1, YLen
		Y(J) = 0.
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Take the # of Particles counts from each of the aerosol particles !!
	CurrentParticle => Particles(XP,YP,ZP)%First !!!!!!!!!!!!!!!!!!!!!!!!!!
	DO J = 1, HowManyLinks
		Y(4*HowManyLinks+7+J) = CurrentParticle%NumberOfParticles
		CurrentParticle		  => CurrentParticle%Next
	END DO


	IF (Scaffolding) CALL TRANSCRIPT("Cond. Check 1")

		
	!!!!!!!!!!!!!!!!!!!!!!!!!
	!! LOOP OVER TIMESTEPS !!
	DO TimeStep = 1, NumTimeSteps

		!! Integrate to this time
		OutTime =  TimeStep * TimeStepSize

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Initialize Chem Storage Array !!
		DO I = 1, NEQ
		DO J = 1, HowManyAqOrganicDissolutionReactions
			ChemStorageArray(J,I) = 0.
		END DO  ;  END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Set LSODES-Related Variables (Notes Excerpted !!
		!! from LSODES file)							 !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! itask  = an index specifying the task to be performed.				!!
		!!          input only.  itask has the following values and meanings.	!!
		!!          1  means normal computation of output values of y(t) at		!!
		!!             t = tout (by overshooting and interpolating).			!!
		!!          2  means take one step only and return.						!!
		!!          3  means stop at the first internal mesh point at or		!!
		!!             beyond t = tout and return.								!!
		!!          4  means normal computation of output values of y(t) at		!!
		!!             t = tout but without overshooting t = tcrit.				!!
		!!             tcrit must be input as rwork(1).  tcrit may be equal to	!!
		!!             or beyond tout, but not behind it in the direction of	!!
		!!             integration.  this option is useful if the problem		!!
		!!             has a singularity at or beyond t = tcrit.				!!
		!!          5  means take one step, without passing tcrit, and return.	!!
		!!             tcrit must be input as rwork(1).							!!
		!!																		!!
		!!          note..  if itask = 4 or 5 and the solver reaches tcrit		!!
		!!          (within roundoff), it will return t = tcrit (exactly) to	!!
		!!          indicate this (unless itask = 4 and tout comes before tcrit,!!
		!!          in which case answers at t = tout are returned first).		!!
		itask  = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! istate = an index used for input and output to specify the
		!!         the state of the calculation.
		!!
		!!         on input, the values of istate are as follows.
		!!          1  means this is the first call for the problem
		!!             (initializations will be done).  see note below.
		!!          2  means this is not the first call, and the calculation
		!!             is to continue normally, with no change in any input
		!!             parameters except possibly tout and itask.
		!!             (if itol, rtol, and/or atol are changed between calls
		!!             with istate = 2, the new values will be used but not
		!!             tested for legality.)
		!!          3  means this is not the first call, and the
		!!             calculation is to continue normally, but with
		!!             a change in input parameters other than
		!!             tout and itask.  changes are allowed in
		!!             neq, itol, rtol, atol, iopt, lrw, liw, mf,
		!!             the conditional inputs ia and ja,
		!!             and any of the optional inputs except h0.
		!!             in particular, if miter = 1 or 2, a call with istate = 3
		!!             will cause the sparsity structure of the problem to be
		!!             recomputed (or reread from ia and ja if moss = 0).
		!!          note..  a preliminary call with tout = t is not counted
		!!          as a first call here, as no initialization or checking of
		!!          input is done.  (such a call is sometimes useful for the
		!!          purpose of outputting the initial conditions.)
		!!          thus the first call for which tout .ne. t requires
		!!          istate = 1 on input.
		!!
		!!          on output, istate has the following values and meanings.
		!!           1  means nothing was done, as tout was equal to t with
		!!              istate = 1 on input.  (however, an internal counter was
		!!              set to detect and prevent repeated calls of this type.)
		!!           2  means the integration was performed successfully.
		!!          -1  means an excessive amount of work (more than mxstep
		!!              steps) was done on this call, before completing the
		!!              requested task, but the integration was otherwise
		!!              successful as far as t.  (mxstep is an optional input
		!!              and is normally 500.)  to continue, the user may
		!!              simply reset istate to a value .gt. 1 and call again
		!!              (the excess work step counter will be reset to 0).
		!!              in addition, the user may increase mxstep to avoid
		!!              this error return (see below on optional inputs).
		!!          -2  means too much accuracy was requested for the precision
		!!              of the machine being used.  this was detected before
		!!              completing the requested task, but the integration
		!!              was successful as far as t.  to continue, the tolerance
		!!              parameters must be reset, and istate must be set
		!!              to 3.  the optional output tolsf may be used for this
		!!              purpose.  (note.. if this condition is detected before
		!!              taking any steps, then an illegal input return
		!!              (istate = -3) occurs instead.)
		!!          -3  means illegal input was detected, before taking any
		!!              integration steps.  see written message for details.
		!!              note..  if the solver detects an infinite loop of calls
		!!              to the solver with illegal input, it will cause
		!!              the run to stop.
		!!          -4  means there were repeated error test failures on
		!!              one attempted step, before completing the requested
		!!              task, but the integration was successful as far as t.
		!!              the problem may have a singularity, or the input
		!!              may be inappropriate.
		!!          -5  means there were repeated convergence test failures on
		!!              one attempted step, before completing the requested
		!!              task, but the integration was successful as far as t.
		!!              this may be caused by an inaccurate jacobian matrix,
		!!              if one is being used.
		!!          -6  means ewt(i) became zero for some i during the
		!!              integration.  pure relative error control (atol(i)=0.0)
		!!              was requested on a variable which has now vanished.
		!!              the integration was successful as far as t.
		!!          -7  means a fatal error return flag came from the sparse
		!!              solver cdrv by way of prjs or slss (numerical
		!!              factorization or backsolve).  this should never happen.
		!!              the integration was successful as far as t.
		!!
		!!          note.. an error return with istate = -1, -4, or -5 and with
		!!          miter = 1 or 2 may mean that the sparsity structure of the
		!!          problem has changed significantly since it was last
		!!          determined (or input).  in that case, one can attempt to
		!!          complete the integration by setting istate = 3 on the next
		!!          call, so that a new structure determination is done.
		!!
		!!          note..  since the normal output value of istate is 2,
		!!          it does not need to be reset for normal continuation.
		!!          also, since a negative input value of istate will be
		!!          regarded as illegal, a negative output value requires the
		!!          user to change it, and possibly other inputs, before
		!!          calling the solver again.
		istate = 1	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! iopt   = an integer flag to specify whether or not any optional		!!
		!!          inputs are being used on this call.  input only.			!!
		!!          the optional inputs are listed separately below.			!!
		!!          iopt = 0 means no optional inputs are being used.			!!
		!!                   default values will be used in all cases.			!!
		!!          iopt = 1 means one or more optional inputs are being used.	!!
		iopt   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! mf     = the method flag.  used only for input.					!!
		!!          the standard choices for mf are..						!!
		!!            mf = 10  for a nonstiff problem,						!!
		!!            mf = 21 or 22 for a stiff problem with ia/ja supplied	!!
		!!                     (21 if jac is supplied, 22 if not),			!!
		!!            mf = 121 for a stiff problem with jac supplied,		!!
		!!                     but not ia/ja,								!!
		!!            mf = 222 for a stiff problem with neither ia/ja nor	!!
		!!                     jac supplied.								!!
		!!          the sparseness structure can be changed during the		!!
		!!          problem by making a call to lsodes with istate = 3.		!!
		mf = 21  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! hmin     = the minimum allowed timestep !!
!		rwork(7) = 1.e-10 !!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! atol   = an absolute error tolerance parameter, either a scalar or	!!
		!!          an array of length neq.  input only.						!!
		!!																		!!
		!!             the input parameters itol, rtol, and atol determine		!!
		!!          the error control performed by the solver.  the solver will	!!
		!!          control the vector e = (e(i)) of estimated local errors		!!
		!!          in y, according to an inequality of the form				!!
		!!                      rms-norm of ( e(i)/ewt(i) )   .le.   1,			!!
		!!          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),			!!
		!!          and the rms-norm (root-mean-square norm) here is			!!
		!!          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))	!!
		!!          is a vector of weights which must always be positive, and	!!
		!!          the values of rtol and atol should all be non-negative.		!!
		!!          the following table gives the types (scalar/array) of		!!
		!!          rtol and atol, and the corresponding form of ewt(i).		!!
		!!																		!!
		!!             itol    rtol       atol          ewt(i)					!!
		!!              1     scalar     scalar     rtol*abs(y(i)) + atol		!!
		!!              2     scalar     array      rtol*abs(y(i)) + atol(i)	!!
		!!              3     array      scalar     rtol(i)*abs(y(i)) + atol	!!
		!!              4     array      array      rtol(i)*abs(y(i)) + atol(i)	!!
		!!																		!!
		!!          when either of these parameters is a scalar, it need not	!!
		!!          be dimensioned in the user-s calling program.				!!
		!!																		!!
		!!          if none of the above choices (with itol, rtol, and atol		!!
		!!          fixed throughout the problem) is suitable, more general		!!
		!!          error controls can be obtained by substituting				!!
		!!          user-supplied routines for the setting of ewt and/or for	!!
		!!          the norm calculation.  see part iv below.					!!
		!!																		!!
		!!          if global errors are to be estimated by making a repeated	!!
		!!          run on the same problem with smaller tolerances, then all	!!
		!!          components of rtol and atol (i.e. of ewt) should be scaled	!!
		!!          down uniformly.												!!
		itol   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! This is a very stiff problem, as aerosol and the gas phase should 
		!! not be the same scale at all.  Screw absolute tolerance.
		!! Could also explore using a vector of atol's
		!!
		!! However, if the initial concentration is zero, some absolute tolerance
		!! is necessary (an error of 0 is unacceptable).
		!!
		!! The units of this are molecules
		atol   = 0.1

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Tell the program which gridpoint it's operating at: !!
		Y(3*HowManyLinks+3) = XP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Y(3*HowManyLinks+4) = YP
		Y(3*HowManyLinks+5) = ZP


		IF (Scaffolding) WRITE(*,*) "Before Hydrophilic Organic Dissolution Loop"
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! LOOP OVER EACH DISSOLUTION REACTION !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		DO Rxn = 1, HowManyAqOrganicDissolutionReactions
		IF (Scaffolding) WRITE(*,*) "Aq Org Reaction #", Rxn

		rtol   = (AqThermoEquilibriumError/10.) **(1./NumTimeSteps) ! AqThermoEquilibriumError
	

		!! If the reaction is coupled to another one, the direct dissolution
		!! is done at the same time as the other one, meaning we skip this
		!! direct one.  Cycle...
		IF (AqOrganicDissolutionData(Rxn, 15) .LT. 0) CYCLE

		!! Reset timestep to the correct position
		T = (TimeStep - 1) * TimeStepSize

		!! Identify the gas index of the condensing chemical
		GasChemIndex1 = AqOrganicDissolutionData(Rxn,1)
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! This is the maximum number of internal steps allowed. !!
		!! The default is 500 and that's not enough.			 !!
		iwork(6) = 6000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICAL																!!
			!! The Y-Vector we will integrate with and push through LSODES looks like this:								!!
	 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!    (1)								    : Gas Phase Concenatration of Chemical AqOrganicDissolutionData(Rxn,1)	!!
			!!    (2 - HowManyLinks+1)				    : Aerosol Concentrations of Chemical AqOrganicDissolutionData(Rxn,3)		!!
			!!	  (HowManyLinks+2   - 2*HowManyLinks+1)	: Aerosol Concentrations of Chemical AqOrganicDissolutionData(Rxn,2) or	!!
			!!											  AqOrganicDissolutionData(Rxn,4), as appropriate						!!
			!!    (2*HowManyLinks+2 - 3*HowManyLinks+1) : Concentration of undissociated species for coupled rxn		!!
			!!    (3*HowManyLinks+2)					: Intentionally Empty											!!
			!!    (3*HowManyLinks+3 - 3*HowManyLinks+5) : Grid Point Specification (x,y,z)								!!
			!!    (3*HowManyLinks+6 - 3*HowManyLinks+7)	: GasChemIndex 1 & 2											!!
			!!    (3*HowManyLinks+8 - 4*HowManyLinks+7) : Condensation Rates for Chemical One							!!
			!!	  (4*HowManyLinks+8 - 5*HowManyLinks+7) : How Many Particles in Each									!!
			!!	  (5*HowManyLinks+8 - 6*HowManyLinks+9) : Pre-computed effective partitioning coefficients						!!
			!!    (6*HowManyLinks+10)					: Flag for what type of reaction								!!
			!!    (6*HowManyLinks+11)					: "-1" (as a reference for finding HowManyLinks)				!!
			!!    (6*HowManyLinks+12)					: AqOrgDissolution Reaction Number									!!
			!!    (6*HowManyLinks+13 - 7*HowManyLinks+12) : Blank													   !!
			!!    (7*HowManyLinks+13 - 8*HowManyLink+12)  : Blank
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! HowManyLinks which the number of aerosol / sections														!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Initialize the Concentrations from the local chemical concentrations: !! 
			Y(1) = GetGridCellChemBurden (XP, YP, ZP, GasChemIndex1) !!!!!!!!!!!!!!!!!!

			!! Tell LSODES which chemical we are to consider
			Y(3*HowManyLinks+6) = GasChemIndex1

			!! The condensation rates
			Y(3*HowManyLinks+8:4*HowManyLinks+7) =AqOrgCondensationRateCoefficients (XP, YP, ZP, Rxn, HowManyLinks)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Reaction Type !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! There are several reaction types in this context, !!
			!! which are different than elsewhere:				 !!!!
			!!  0. Water Condensation								!!
			!!  1. Direct Dissolution of a single chemical			!!
			!!  2. Dissolution of a single species that dissociates	!!
			!!  3. Coupled reactions of types 1 & 2					!!
			!!  4. Dissolution of a species that binds with an ion	!!
			!!  5. Coupled reactions of types 1 & 4	
			!!  7. Sulfate condensation (eq. gas conc. is 0)				!!
			Y(6*HowManyLinks+10) = AqOrganicDissolutionData(Rxn,11) !!!!!!!!!!!

			!WRITE(*,*) "AqOrgDissolution Reaction #", TRIM(INT2STR(Rxn)), "Type: ", TRIM(Real2STR(Y(6*HowManyLinks+10)))
			
			!! A "-1" as reference for later subroutines
			Y(6*HowManyLinks+11) = -1

			!! The reaction number
			Y(6*HowManyLinks+12) = Rxn


			!! Respecify NEQ based on reaction type
			!! And then describe the Jacobian's structure (Ia and Ja, goes into IWORK)
			IF (INT(Y(6*HowManyLinks+10)) .EQ. 9) THEN

				NEQ = HowManyLinks+1

				!! Specify Ia and Ja
				IWORK(30+1)     = 1       !! first Ia
				IWORK(30+NEQ+1) = 3*neq-1 !! Final Ia Term

				IWORK(32+NEQ)   = 1		  !! first Ja
						
				DO J = 1, NEQ-1
					IWORK(31+J) = NEQ + (J*2-1)  !! Ia
					IWORK(32+NEQ+J)       = J+1  !! Ja (fill the 1st column)
					IWORK(31+(NEQ+J)*2-1) = 1    !! Ja (fill the 1st row)
					IWORK(31+(NEQ+J)*2)   = J+1  !! Ja (fill the diagonal)
				END DO

			ELSE
				CALL ERROR("All hydrophilic organic dissolution reactions should be of type 9.")
			END IF


			!! Fill the remainder of the Y vector off-site
			CALL AqOrgDissolutionFactors (XP, YP, ZP, Rxn, HowManyLinks, Temperature, INT(Y(6*HowManyLinks+10)), 	 &
			      Y(2:1+HowManyLinks), Y(HowManyLinks+2:2*HowManyLinks+1), Y(2*HowManyLinks+2:3*HowManyLinks+1), &  ! initial concs
				  Y(HowManyLinks*5+8:HowManyLinks*6+7), Y(HowManyLinks*6+13:HowManyLinks*7+12), ReturnType, Y(HowManyLinks*7+13:HowManyLinks*8+12))


			!! Store the concentrations for later
			DO J = 1, NEQ
				InitialConcsArray(J) = Y(J)
			END DO

			!! Check to see if there is any of the chemical in the system
			II = 0. ; JJ = 0.
			DO K = 1, HowManyLinks
				II = II + Y(1+K)*Y(4*HowManyLinks+7+K)
				JJ = JJ + Y(1+HowManyLinks+K)*Y(4*HowManyLinks+7+K)
				KK = KK + Y(1+2*HowManyLinks+K)*Y(4*HowManyLinks+7+K)
			END DO
	
			!! Cycle if conc. too low
			SELECT CASE (INT(AqOrganicDissolutionData(Rxn,11)))
			CASE (:1, 8, 9)
				IF (Y(1)+II .LT. 1.0) CYCLE
			CASE (2)
				IF (Y(1) .EQ. 0 .AND. (II .EQ. 0 .OR. JJ .EQ. 0)) CYCLE
			CASE (3)
				IF (Y(1) .EQ. 0 .AND. (II .EQ. 0 .OR. JJ .EQ. 0) .AND. KK .EQ. 0) CYCLE
			CASE (4)
				IF ((Y(1) .EQ. 0 .OR. JJ .EQ. 0) .AND. II .EQ. 0) CYCLE
			CASE (5)
				IF ((Y(1) .EQ. 0 .OR. JJ .EQ. 0) .AND. II .EQ. 0 .AND. KK .EQ. 0) CYCLE
			END SELECT

		
			!! Each chemical starts as if a new run
			istate = 1

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
30			CALL LSODES (AqOrgCondensationODEEvaluator,neq,y,t,outtime,itol,rtol,atol,itask, &
						 istate,iopt,rwork,lrw,iwork,liw,AqOrgCondensationJacobianEvaluator,mf)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
			IF (Scaffolding) WRITE(*,*) "After AqOrgCondensation LSODES", Rxn
			!! If the maximum number of steps occurred, then warn the user (in Scaffolding
			!! mode) and push the integration back into LSODES until the appropriate end time 
			!! is reached.
			IF (ISTATE .EQ. -1 ) THEN !.AND. T < OUTTIME) THEN
				IF (Scaffolding) &
				CALL WARN("Warning! Nominal maximun number of steps in LSODES was exceeded in AqOrgCondensation Integration for "	&
				 	     //TRIM(GasPhaseChemicalNames(GasChemIndex1))//" ("//trim(int2str(iwork(6)))//" steps)... Recycling")
				ISTATE = 2
				GOTO 30
			END IF

			!! An Error will occur if the array lengths were too small
			IF (IWORK(17) > LRW)  &
				CALL ERROR("The Given Length for LSODES' Working Array RWORK (accessing in StepAqOrgCondensation())"//		&
						  " Should have been of size "//TRIM(INT2STR(iwork(17)))//" but was only "//					&
						  TRIM(INT2STR(LRW))//" instead.  Fix the assignment in EvolveGasChemistry's SetEvolveGasChemConstants().")

			IF (ISTATE .LT. -1) &
			CALL ERROR("In StepAqOrgCondensation(), LSODES encountered a problem during the integration and returned ISTATE = "//	&
						TRIM(INT2STR(ISTATE))//".  Please investigate.")
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Store the output concentrations until we !!
			!! have considered all of the chemicals		!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! This is indexed as follows:												!!
			!! ChemStorageArray(Rxn, 1)         = Change in Gas Phase Conc 1			!!
			!! ChemStorageArray(Rxn, 2)         = Change in Aerosol Conc 1 (1)			!!
			!!																			!!
			!! ChemStorageArray(Rxn, HML+1)     = Change in Aerosol Conc 1 (HowManyLinks)!
			!! ChemStorageArray(Rxn, HML+2)     = Change in Aerosol Conc 2 (1)			!!
			!!																			!!
			!! ChemStorageArray(Rxn, 2*HML+1)   = Change in Aerosol Conc 2 (HowManyLinks)!
			!! ChemStorageArray(Rxn, 2*HML+2)   = Change in Gas Phase Conc 2			!!
			DO K = 1, NEQ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				ChemStorageArray(Rxn,K) = Y(K) - InitialConcsArray(K)
			END DO

		END DO  !! Ends the loop over all the dissolution reactions
		IF (Scaffolding) WRITE(*,*) "After Dissolution Loop"

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Replace the background chemical concentrations !!

		!! Load an appropriate vector to send to the gas phase gridpoint
		DO I = 1, HowManyAqOrganicDissolutionReactions
			!Skip replace if gas is to be held constant
			IF(INT(AqOrganicDissolutionData(I,1)) .LE. HowManyEvolveGasChems) THEN
				CALL AddToGridCellChemBurden (XP, YP, ZP, INT(AqOrganicDissolutionData(I,1)), ChemStorageArray(I,1))
			END IF
		END DO
		IF (Scaffolding) WRITE(*,*) "After Replace Gas Chems"


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Walk the linked list and replace the chemicals in each particle !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CurrentParticle => FirstParticle

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Replace the chemicals now at the end of the !!
		!! integrations in one pass through the list.  !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO I = 1, HowManyLinks  !! Loop over the aerosol
		DO J = 1, HowManyAqOrganicDissolutionReactions

			CurrentParticle%AqOrgChems(AqOrganicDissolutionData(J,3)) = CurrentParticle%AqOrgChems(AqOrganicDissolutionData(J,3)) + &
															ChemStorageArray(J,I+1) / Avogadro
			!Check for zero
			IF(CurrentParticle%AqOrgChems(AqOrganicDissolutionData(J,3)) .LT. 0.) THEN
				Store = CurrentParticle%AqOrgChems(AqOrganicDissolutionData(J,3))
				CurrentParticle%AqOrgChems(AqOrganicDissolutionData(J,3)) = 0.
				CALL AddToGridCellChemBurden (XP, YP, ZP, INT(AqOrganicDissolutionData(I,1)), Store*Avogadro*CurrentParticle%NumberofParticles)
			END IF

		END DO 

		!Recalculate activity coefficients
		CALL UpdateHydrophilicUNIFAC(CurrentParticle, Temperature)
		IF (Scaffolding) WRITE(*,*) "After UNIFAC"	

		!Equilibrate organics across phases
		!NOTE: This should only be done ifBoth phases are allowed 
		IF(DoHydrophobicOrgDissolutionFlag .EQ. .TRUE. .AND. DoHydrophilicOrgDissolutionFlag .EQ. .TRUE.) THEN
				CALL EquilibrateOrganicParticle (CurrentParticle, Temperature, ReturnType, OptErrorTolerance=AqThermoEquilibriumError/2.)	 
		END IF
		IF (Scaffolding) WRITE(*,*) "After Org Eq"

		CurrentParticle => CurrentParticle%Next

		END DO !! Quit loop over particles for replacement

	END DO !! Stop Looping Over the Time Steps

	!! Get rid of everything that was allocated in here.
	DEALLOCATE (Y,RWORK,IWORK,InitialConcsArray,ChemStorageArray)
		
	!! SCSCSCSC
	IF (Scaffolding) CALL TRANSCRIPT (">>>Exiting StepAqOrgCondensation<<<")

	RETURN
END SUBROUTINE StepAqOrgCondensation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This program provides the ODE ydot vector to LSODES. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AqOrgCondensationODEEvaluator (neq, t, y, ydot)

		USE Aerosols,            ONLY : Particles, Particle, ParticleArray
		USE Chemistry,			 ONLY : HowManyGasChems
		USE InfrastructuralCode, ONLY : Transcript, ERROR
		USE ModelParameters,     ONLY : Avogadro

		IMPLICIT NONE

		!! External Variable (prescribed by LSODES)
		REAL*8		:: t, y, ydot
		INTEGER		:: neq
		DIMENSION	:: y(1), ydot(1)

		!! Internal Variables
		INTEGER :: I, HowManyLinks

		LOGICAL, PARAMETER :: Scaffolding =  .FALSE.

IF (SCAFFOLDING) print *, ">>> In AqOrgCondensationODEEvaluator <<<"

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET THE Y-VECTOR FOR THE CONDENSING CHEMICAL																!!
		!! The Y-Vector we will integrate with and push through LSODES looks like this:								!!
 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!    (1)								    : Gas Phase Concenatration of Chemical AqOrganicDissolutionData(Rxn,1)	!!
		!!    (2 - HowManyLinks+1)				    : Aerosol Concentrations of Chemical AqOrganicDissolutionData(Rxn,3)		!!
		!!	  (HowManyLinks+2   - 2*HowManyLinks+1)	: Aerosol Concentrations of Chemical AqOrganicDissolutionData(Rxn,2) or	!!
		!!											  AqOrganicDissolutionData(Rxn,4), as appropriate						!!
		!!    (2*HowManyLinks+2 - 3*HowManyLinks+1) : Concentration of undissociated species for coupled rxn		!!
		!!    (3*HowManyLinks+2)					: Intentionally Empty											!!
		!!    (3*HowManyLinks+3 - 3*HowManyLinks+5) : Grid Point Specification (x,y,z)								!!
		!!    (3*HowManyLinks+6 - 3*HowManyLinks+7)	: GasChemIndex 1 & 2											!!
		!!    (3*HowManyLinks+8 - 4*HowManyLinks+7) : Condensation Rates for Chemical One							!!
		!!	  (4*HowManyLinks+8 - 5*HowManyLinks+7) : How Many Particles in Each									!!
		!!	  (5*HowManyLinks+8 - 6*HowManyLinks+9) : Pre-computed effective Henry's Law Coeffs						!!
		!!    (6*HowManyLinks+10)					: Flag for what type of reaction								!!
		!!    (6*HowManyLinks+11)					: "-1" (as a reference for finding HowManyLinks)				!!
		!!    (6*HowManyLinks+12)					: Dissolution Reaction Number									!!
		!!    (6*HowManyLinks+13 - 7*HowManyLinks+12) : Pre-computed effective Henry's Law Coeffs for coupled rxn   !!
		!!    (7*HowManyLinks+13 - 8*HowManyLink+12)  : Flags for HNO3 and HCl dissolution	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! HowManyLinks which the number of aerosol / sections														!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF (Y(2*neq+9) .EQ. -1) THEN
			HowManyLinks = (neq - 1)/3
		ELSE IF (Y(3*neq+8) .EQ. -1) THEN
			HowManyLinks = (neq - 1)/2
		ELSE IF (Y(6*neq+5) .EQ. -1) THEN
			HowManyLinks = neq - 1
		ELSE
			CALL ERROR ("Couldn't back howmanylinks out of neq in CondensationODEEvaluator()...")
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! The gas phase concentration changes as !!
		!! the negative of the sum of the other	  !!
		!! rates. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		YDOT(1) = 0.

		!WRITE(*,*) "Case: ", Y(6*HowManyLinks+10), 
		SELECT CASE (INT(Y(6*HowManyLinks+10)))

		!! DISSOLUTION of a single, non-dissociating organic species into aqueous phase
		CASE (9)

			DO I = 2, neq
				YDOT(I) = Y(3*HowManyLinks+6+I) *	&		! The Condensation Rate
						  (Y(1) -					&		! The Gas Phase Conc
						   Y(5*HowManyLinks+6+I) * Y(I))	! S_i' * c_s,i / H_i

				YDOT(1) = YDOT(1) - YDOT(I) * Y(4*HowManyLinks+6+I)
			END DO


		END SELECT

	RETURN
END SUBROUTINE AqOrgCondensationODEEvaluator



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate a JACOBIAN for Condensation !!
!!                                       !!
!! As per LSODES's specification, only   !!
!! calculate at the stated row I         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AqOrgCondensationJacobianEvaluator (neq, t, y, I, ia, ja, pdj)

	USE GridPointFields,     ONLY : GasPhaseChemicalRates
	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	!! External Variables
	REAL*8		:: t, y, pdj
	INTEGER		:: I, neq, ia, ja
	DIMENSION	:: y(1), ia(1), ja(1), pdj(1)

	!! Internal Variables
	INTEGER :: J, HowManyLinks

	LOGICAL, PARAMETER :: SCAFFOLDING = .FALSE.

	IF(SCAFFOLDING) WRITE(*,*) "Entering AqOrgCondensation Jacobian."
	!! Determine HowManyLinks by backing out of the -1 reference point
	IF (Y(2*neq+9) .EQ. -1) THEN
		HowManyLinks = (neq - 1)/3
	ELSE IF (Y(3*neq+8) .EQ. -1) THEN
		HowManyLinks = (neq - 1)/2
	ELSE IF (Y(6*neq+5) .EQ. -1) THEN
		HowManyLinks = neq - 1
	ELSE
		CALL ERROR ("Couldn't back howmanylinks out of neq in CondensationODEEvaluator()...")
	END IF

	!!!NOTE: I = 1 is the gas-phase equation, I .GT. 1 is an aerosol index (which repeats
	!!!if more than 1 aerocol phase conpound is considered.)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! If this is for the gas phase !!
	!! chemical concentration then  !!
	!! it is different              !!
	IF (I .EQ. 1) THEN  !!!!!!!!!!!!!!

		!! Every case is the same here
		DO J = 2, HowManyLinks+1
			PDJ(J)                = Y(3*HowManyLinks+6+J) !! k_i 

			IF (Y(6*HowManyLinks+10) .GT. 1) &
				PDJ(J+HowManyLinks) = Y(3*HowManyLinks+6+J) !! k_i 

			IF (Y(6*HowManyLinks+10) .EQ. 3 .OR. Y(6*HowManyLinks+10) .EQ. 5) &
				PDJ(J+2*HowManyLinks) = Y(3*HowManyLinks+6+J) !! k_i 
		END DO

		DO J = 2, HowManyLinks+1
			PDJ(1) = PDJ(1) - PDJ(J) * Y(4*HowManyLinks+6+J)  !! - k_i * (# of Particles)
		END DO


	!! There are several structures for the columns beyond the first
	!! that depend on reaction type
	ELSE

		SELECT CASE (INT(Y(6*HowManyLinks+10)))

		CASE (9) !! Dissolution of organic into aqueous phase 

				     !!		       k_j                Henry's Law  
			PDJ(I) = -1 * Y(3*HowManyLinks+6+I) * Y(5*HowManyLinks+6+I)

					 !!			       (# of Particles)
			PDJ(1) = -1. * PDJ(I) * Y(4*HowManyLinks+6+I) !! - d2 C_j / dt d c_j 


		END SELECT
		END IF

	RETURN
END SUBROUTINE AqOrgCondensationJacobianEvaluator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! AqOrgDissolutionFactors --									!!
!!															!!
!! Give the dissolution routines everything they need to	!!
!! calculate one iteration of one chemical, including:		!!
!!															!!
!!  1. Effective Eq Constant								!!
!!  2. Initial Appropriate Aerosol Concentrations			!!
!!															!!
!! ReturnType = 1 if fine, -1 if there is none of			!!
!! the chemical in the gridpoint							!!
SUBROUTINE AqOrgDissolutionFactors (XP, YP, ZP, WhichRxn, HowManyLinks, Temperature, ReactionType,	&
							   InitialConcs1, InitialConcs2, InitialConcs3, SoverH, SoverH2, ReturnType, RealFlag)


	USE GridPointFields, ONLY :		GasPhaseChemicalRates,			&
									GetSatVapBurden,				&
									LocateGridPoint,				&
									GasPhaseChemGridIndex,			&
									GetGridCellChemBurden

	USE Aerosols,        ONLY :		Particles, Particle,			&
									ParticleArray, Molality

	USE ModelParameters, ONLY :		RstarMB, moles, grams, watermolecmass,	&
									ProtonIndex, HydroxyIndex,				&
									Avogadro

	USE Chemistry,       ONLY :		HowManyAqChems,HowManyAqCations,		&
									HowManyAqAnions,						&
									AqEquilibriaList,						&
									GasPhaseChemicalNames,					&
									AqPhaseChemicalNames,					&
									AqCationNames,AqAnionNames, & 
									FindChem

	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	!! External Variables, in
	INTEGER :: XP, YP, ZP, WhichRxn, HowManyLinks
	REAL*8  :: Temperature

	!! External Variables, out
	REAL*8  :: SoverH(HowManyLinks), SoverH2(HowManyLinks), InitialConcs1(HowManyLinks), &
			   InitialConcs2(HowManyLinks), InitialConcs3(HowManyLinks), RealFlag(HowManyLinks)
	INTEGER :: ReturnType, ReactionType, Flag, IonA, IonB, ExpA, ExpB

	!! Internal Variables
	INTEGER :: InternalHowManyLinks, I, J, Allocation_Error, GasChemIndex
	REAL*8  :: II, XX(3), H2SO4AqIndex, AnionSatConc, GasSaturationConc, ProtonMolality

	TYPE(Particle), POINTER :: CurrentParticle

	!! Identify the gas index of the condensing chemical
	GasChemIndex = AqOrganicDissolutionData(WhichRxn,1)

	!! An initial guess
	ReturnType = -1

	!! If Water Condensation, we need a couple of extra items...
	IF (ReactionType .EQ. 0) THEN
		XX = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
		II = GetSatVapBurden(XX(1),XX(2),XX(3))
	END IF

	!! If there is a positive gas phase concentration, then return positive return flag.
	IF (ReactionType .EQ. 0.) ReturnType = 1

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Take the burdens and # of Particles counts from each of the aerosol particles !!
	CurrentParticle => Particles(XP,YP,ZP)%First !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO J = 1, HowManyLinks

		SELECT CASE (ReactionType)

		
		!! Direct dissolution into AqOrganic phase
		CASE (9)

			IF (ReturnType .EQ. -1 .AND. CurrentParticle%AqOrgChems(AqOrganicDissolutionData(WhichRxn,3)) .GT. 0.)	&
			ReturnType = 1

			ProtonMolality = Molality (ProtonIndex, CurrentParticle)
	
			SoverH(J) =																	 &
					SatVapPressRatio (CurrentParticle, INT(AqOrganicDissolutionData(WhichRxn,3))) &
					* CurrentParticle%HydrophilicActivityCoeffs(AqOrganicDissolutionData(WhichRxn,3)) &
					* 1000 /  EffectiveOrganicHenrysLaw  (Temperature, ProtonMolality, WhichRxn) &
					/ RstarMB / Temperature / CurrentParticle%AqChems(1) / watermolecmass

			!! Load the Y vector with the appropriate aerosol    !!
			!! concentration.									 !!
			InitialConcs1(J) = CurrentParticle%AqOrgChems(AqOrganicDissolutionData(WhichRxn,3)) * Avogadro
			InitialConcs2(J) = 0.
			InitialConcs3(J) = 0.
		
		END SELECT

		CurrentParticle => CurrentParticle%Next
	END DO

	RETURN
END SUBROUTINE AqOrgDissolutionFactors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate MASS FLUX CONDENSATION RATE COEFFICIENT for org disolution
!! into hydrophobic organic phase. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION AqOrgCondensationRateCoefficients (XP, YP, ZP, ReactionNumb, NumbLinks)

	USE Chemistry,		 ONLY : EnthalpyOfVaporization,			&
								GasMolecularmass,				&
								MolecularDiffusionCoefficient

	USE GridPointFields, ONLY : LocateGridPoint,				&
								GasPhaseChemGridIndex,			&
								GetAirDensity,					&
								GetCpm,							&
								GetDynamicViscosityOfAir,		&
								GetTemp,						&
								MolecularThermalConductivity,	&
								GetSatVapConcentration

	USE InfrastructuralCode, ONLY : ERROR, INT2STR

	USE Aerosols,		 ONLY : Particle, Particles,		&
								ParticleKnudsenForEnergy,	&
								ReynoldsNumber

	USE ModelParameters, ONLY : AirMolecMass,				&
								Avogadro,					&
								Pi,							&
								ThermalAccommodationCoeff,	&
								Rstar,						&
								ChemScale,					&
								ThermoBulkMode,				&
								IgnoreCorrectionstoDiffusivity

    IMPLICIT NONE

	!! Define input and output variables
	INTEGER	:: ReactionNumb, NumbLinks, XP, YP, ZP
	REAL*8	:: AqOrgCondensationRateCoefficients(NumbLinks)

	TYPE(Particle),POINTER :: Current
	INTEGER ::  I, ChemIndex
	REAL*8  ::	AirDensity,		&
				AirTemp,		&
				Cpm,			&	! Moist specific heat
				Dv,				&	! Diffusion Coeffienct of Condensing Gas
				DvPrime,		&
				DynVisc,		&
				EnergyKnud,		&	
				GasKnud,		&
				GasThermalVel,	&
				KappaD,			&	! Thermal conductivity of condensing gas
				KappaDPrime,	&
				ReynoldsNumb,	&
				SchmidtNumber,	&
				X,				&
				GasMeanFreePath,&
				SatWaterVapConc, &
				XX(3)
	
	LOGICAL :: Scaffolding = .FALSE.

	!! Locate the gridpoint in 3-space
	XX = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)

	!! Get the primary condensing chemical
	ChemIndex = AqOrganicDissolutionData(ReactionNumb,1)

	!! Retrieve or Calculate the Variables that are the 
	!! same for all aerosol
	AirTemp			= GetTemp(XX(1),XX(2),XX(3))
	AirDensity		= GetAirDensity(XX(1),XX(2),XX(3))
	Dv				= MolecularDiffusionCoefficient(XX(1),XX(2),XX(3),ChemIndex)
	DynVisc			= GetDynamicViscosityOfAir(XX(1),XX(2),XX(3))
	GasThermalVel   = GetThermalVelocity(XX(1),XX(2),XX(3),ChemIndex)
	SatWaterVapConc	= GetSatVapConcentration(XX(1),XX(2),XX(3))


	!! Non water condensation
	IF (ReactionNumb .GT. 0) THEN
		Cpm	   = GetCpm(XX(1),XX(2),XX(3))
		KappaD = MolecularThermalConductivity(XX(1),XX(2),XX(3))
	END IF


	!! Constant is 32. / 3. / Pi
	GasMeanFreePath = 3.395305453 * Dv * AirMolecMass / GasThermalVel /			&
					  (AirMolecMass + GasMolecularMass(ChemIndex) / Avogadro)


	!! Get the first particle
	Current => Particles(XP,YP,ZP)%First

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Loop over all of the particles !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO I = 1, NumbLinks

	!! Walk the chain and err if can't find the next link
	IF (I .GT. 1) THEN
		IF (.NOT. ASSOCIATED(Current%Next)) &
			CALL ERROR("Miscounted the number of particles at point "//TRIM(INT2STR(XP))//", "//TRIM(INT2STR(YP))//", "//	&
					   TRIM(INT2STR(ZP))//" in the argument to CondensationRateCoefficients().  This just shouldn't happen.")
		Current => Current%Next
	END IF

	!! If it's a section, make sure there's something in it
	!! Or if the particle can't attract water
	IF ((Current%Sectional .AND. Current%NumberOfParticles .EQ. 0) .OR. Current%Dry) THEN
		AqOrgCondensationRateCoefficients(I) = 0.
		CYCLE
	END IF

	!! The ones that are not the same for each
	ReynoldsNumb = ReynoldsNumber(Current)    ! defined for a particle by its radius
	EnergyKnud	 = ParticleKnudsenForEnergy(Current)


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Water is not the same other Chems. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Wa are calculating k_i for a Dissolution ODE that looks like this: !!
	!!                                                                    !!
	!!  d c_i          /               c_i   \                            !!
	!!  -----  =  k_i |  C_i - S'_i * -----   |                           !!
	!!    dt           \                H'   /                            !!
	!!																	  !!
	!! Condensation is similar, as noted below.							  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The Dissolution Rate Coefficient   !!
	!! is this (condensation is corrected !!
	!! further, see below):				  !!
	!!									  !!
	!! k_i = 4 pi r D_i'                  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The MOLECULAR DIFFUSION COEFFICIENT must be corrected !!
	!! for collisional geometry and ventilation.			 !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! The Gas Knudsen Number, (place this after grabbing the uncorrected Dv)
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

	!! Correct Dv away from the continuum solution for COLLISIONAL GEOMETRY and STICKING PROBABILITY
	!! Following Jacobson, 1999, who follows Fuchs and Sutugin 1971 and Pruppacher and Klett 1997
	DvPrime = DvPrime/(1. + ((1.33+0.71/GasKnud)/(1.+1./GasKnud)+ &
			  4*(1.-AqOrganicDissolutionData(ReactionNumb,10))/			  &
			  (3*AqOrganicDissolutionData(ReactionNumb,10)))*GasKnud)

	!Ignore corrections if ordered to
	IF(IgnoreCorrectionstoDiffusivity) DvPrime = Dv

	!! Assemble the k_i for that section / particle
	AQOrgCondensationRateCoefficients(I) = Current%EffectiveRadius*4*Pi*DvPrime

	END DO  !! Move to another particle

	RETURN
END FUNCTION  AqOrgCondensationRateCoefficients

SUBROUTINE APDNonDiss(Y, Timestep, neq)
!!Copyright Matt Alvarado, 8/30/2006
!!Based on Jacobson's APD scheme for non-dissociating compounds
!!(See Jacobson, M.Z., "Fundamentals of Atm. Modeling", 2nd ed., 2005
!! pp. 583-585)
!!Integrates condensation for a single aqueous organic forward one time step

	IMPLICIT NONE

	!Input Variables
	REAL*8		:: Y, Timestep
	INTEGER		:: neq
	DIMENSION	:: Y(1)

	!Internal Variables
	REAL*8		:: Sum1, Sum2
	REAL*8		:: Factor(neq-1), Recip(neq-1)
	INTEGER		:: I, HowManyLinks

	HowManyLinks = neq - 1
	
	Sum1 = 0.0
	Sum2 = 0.0
	
	DO I = 2,neq
				    !S/H		
	Factor(I-1) = EXP(-Timestep*Y(5*HowManyLinks+6+I) &
				      *Y(3*HowManyLinks+6+I)) !k
		Sum1 = Sum1 + Y(I)*(1.0-Factor(I-1))
		Recip(I-1) = (1.0/Y(5*HowManyLinks+6+I)) !H/S
		Sum2 = Sum2 + Recip(I-1)*(1.0-Factor(I-1))
	END DO 
	
	!Gas-phase concentration
	Y(1) = (Y(1) + Sum1)/(1+Sum2)
	
	!Aerosol phase concentrations
	DO I = 2, neq
		Y(I) = (Y(I) - Recip(I-1)*Y(1))*Factor(I-1)
		Y(I) = Y(I) + Recip(I-1)*Y(1)
	END DO
	RETURN
END SUBROUTINE APDNonDiss
