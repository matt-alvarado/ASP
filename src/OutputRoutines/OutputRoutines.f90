!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! OutputRoutines.f90
!! This file contains many useful ways for summarizing aerosol
!! data for output to ASCII text files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!! 08/29  2006   Matt Alvarado     1. Changed TotalPOA and TotalSOA				!!
!!									calculation in SizeDistOutput				!!
!!								   2. Changed currentparticle to current		!!
!!									 in SizeDistOutput							!!
!! 09/01  2006   Matt Alvarado     1. Added TotalOrgCar calculation to			!!
!!										SizeDistOutput							!!
!! 09/07  2006   Matt Alvarado     Updated SpillBeans to spill BoundaryParticles!!
!! 11/01  2006   Matt Alvarado     Added Optical properties to SpillBeans		!!
!! 10/03  2007   Matt Alvarado     Added lines to SpillBeans to help with 
!!									initializing CRM6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES						        !!
!! 1. ModelParameters						!!
!! 2. InfrastructuralCode					!!
!! 3. GridPointFields						!!
!! 4. Chemistry								!!
!! 5. Aerosols								!!
!! 6. Condensation							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:			!!	
!! 1. SUBROUTINE EachAerosolOutput (OutputFileName)
!! 2. SUBROUTINE BinEntireDomainForPrinting (OutputFileName, inHowManyBins)
!! 3. SUBROUTINE SummarizeDomainMassForOutput (FileName, FirstCall)
!! 4. SUBROUTINE BinLagrangianParticlesForOutput (XP,YP,ZP,NumberOfBins,SmallEdge,LargeEdge, FH, FirstComment)
!! 5. SUBROUTINE SpillBeans(XP, YP, ZP)
!! 6. SUBROUTINE SizeDistOutput(XP, YP, ZP, t, Filename)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE OutputRoutines

PRIVATE
PUBLIC :: BinEntireDomainForPrinting,		&
		  SummarizeDomainMassForOutput,		&
		  BinLagrangianParticlesForOutput,	&
		  EachAerosolOutput, &
		  SpillBeans, &
		  SizeDistOutput

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Produces a data file that is appropriate		 !!
!! for printing the aerosol distribution of the  !!
!! entire domain ensemble aerosol.               !!
SUBROUTINE EachAerosolOutput (OutputFileName)

	USE GridPointFields,	 ONLY : XGridPoints, YGridPoints, ZGridPoints, GasPhaseChemGridIndex

	USE ModelParameters,	 ONLY : OutputDeckSubDir, micron, grams, moles

	USE Aerosols,            ONLY : Particles, Particle, AerosolPH
	
	USE Chemistry,           ONLY : HowManyAqChems,					&
									HowManyAqCations,				&
									HowManyAqAnions,				&
									AqPhaseChemicalNames,			&
									AqCationNames,AqAnionNames

	USE Condensation,        ONLY : IsActivated

	USE InfrastructuralCode, ONLY : WARN, ERROR, INT2STR, REAL2STR,	 &
									GetFileHandle, ReturnFileHandle

	IMPLICIT NONE

	CHARACTER *(*) :: OutputFileName
	TYPE(Particle), POINTER :: Current
	LOGICAL :: ChemicalPresent(HowManyAqChems+HowManyAqCations+HowManyAqAnions)
	INTEGER :: I, X, Y, Z, FH

	!! Fill the Logical Array "ChemicalPresent()" with trues if the
	!! chemical is found in at least one aerosol particle.  This will
	!! be used in a second pass through the particles for outputting 
	!! their contents.
	DO I = 2, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		ChemicalPresent(I) = .FALSE.
	END DO
	ChemicalPresent(1) = .TRUE. !! Assume water is present

	DO X = 1, XGridPoints(GasPhaseChemGridIndex)
	DO Y = 1, YGridPoints(GasPhaseChemGridIndex)
	DO Z = 1, ZGridPoints(GasPhaseChemGridIndex)

		Current => Particles(X,Y,Z)%First

10		IF (ASSOCIATED(Current)) THEN
			DO I = 2, HowManyAqChems+HowManyAqCations+HowManyAqAnions
				IF (Current%AqChems(I) .GT. 0.) ChemicalPresent(I) = .TRUE.
			END DO
		END IF

		IF (ASSOCIATED(Current%Next)) THEN
			Current => Current%Next
			GOTO 10
		END IF
	END DO ; END DO ; END DO

!DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
!	Print *, I,ChemicalPresent(I) 
!END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Output the traits of each particle !!
	!! to a file                          !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FH = GetFileHandle()
    OPEN(UNIT=FH, FILE=TRIM(OutputDeckSubDir)//TRIM(OutputFileName)//".csv")
	WRITE (FH, '(a)', advance='no') "ParticleID, "//			 &
									"Particle Distribution #,"// &
									"Number Of Particles,"//	 &
									"Activated (1 is yes 0 is no),"// &
									"Effective Radius,"//		 &
									"Solution Embryo Radius,"//	 &
									"InsolCoreRadius,"//		 &
									"CosContactAngle,"//		 &
									"InsolCoreMass,"//			 &
									"Solution Density,"//		 &
									"Particle Density,"//		 &
									"Surface Tension,"//		 &
									"PH,"//						 &
									"Origin Pressure Level"

	DO I = 1, HowManyAqChems
		IF (ChemicalPresent(I)) &
		WRITE (FH, '(a)', advance='no') ","//TRIM(AqPhaseChemicalNames(I))
	END DO
	DO I = 1, HowManyAqCations
		IF (ChemicalPresent(I+HowManyAqChems)) &
		WRITE (FH, '(a)', advance='no') ","//TRIM(AqCationNames(I))
	END DO
	DO I = 1, HowManyAqAnions
		IF (ChemicalPresent(I+HowManyAqChems+HowManyAqCations)) &
		WRITE (FH, '(a)', advance='no') ","//TRIM(AqAnionNames(I))
	END DO

	WRITE (FH, '(a)') ""

	DO X = 1, XGridPoints(GasPhaseChemGridIndex)
	DO Y = 1, YGridPoints(GasPhaseChemGridIndex)
	DO Z = 1, ZGridPoints(GasPhaseChemGridIndex)

		Current => Particles(X,Y,Z)%First

20		IF (ASSOCIATED(Current)) THEN

			WRITE (FH, '(a)', advance='no') TRIM(INT2STR(Current%ParticleID))//", "//			&
										    TRIM(INT2STR(Current%ParticleDistribution))//", "// &
											TRIM(REAL2STR(Current%NumberOfParticles,8))//", "

			
			WRITE (FH, '(a)', advance='no') TRIM(REAL2STR(Current%EffectiveRadius/micron,8))//", "//	&
											TRIM(REAL2STR(Current%EmbryoRadius/micron,8))//", "//			&
											TRIM(REAL2STR(Current%SolutionDensity,8))//", "//			&
											TRIM(REAL2STR(Current%ParticleDensity,8))//", "//			&
											TRIM(REAL2STR(Current%SurfaceTension,8))//", "//			&
											TRIM(REAL2STR(AerosolPH(Current)))//", "//					&
											TRIM(REAL2STR(Current%OriginPressureLevel))

			DO I = 1, HowManyAqChems
				IF (ChemicalPresent(I)) &
				WRITE (FH, '(a)', advance='no') ","//TRIM(REAL2STR(Current%AqChems(I)/moles))
			END DO
			DO I = 1, HowManyAqCations
				IF (ChemicalPresent(I+HowManyAqChems)) &
				WRITE (FH, '(a)', advance='no') ","//TRIM(REAL2STR(Current%AqChems(I+HowManyAqChems)/moles))
			END DO
			DO I = 1, HowManyAqAnions
				IF (ChemicalPresent(I+HowManyAqChems+HowManyAqCations)) &
				WRITE (FH, '(a)', advance='no') ","//TRIM(REAL2STR(Current%AqChems(I+HowManyAqChems+HowManyAqCations)/moles))
			END DO

			WRITE (FH, '(a)') ""

		END IF


		IF (ASSOCIATED(Current%Next)) THEN
			Current => Current%Next
			GOTO 20
		END IF
	END DO ; END DO ; END DO

	CLOSE(FH)
	CALL ReturnFileHandle(FH)

	RETURN
END SUBROUTINE EachAerosolOutput 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Produces a data file that is appropriate		 !!
!! for printing the aerosol distribution of the  !!
!! entire domain ensemble aerosol.               !!
SUBROUTINE BinEntireDomainForPrinting (OutputFileName, inHowManyBins)

	USE GridPointFields,	 ONLY : XGridPoints, YGridPoints, ZGridPoints, GasPhaseChemGridIndex

	USE Aerosols,            ONLY : Particles, Particle,					&
									ParticleArray,							&
									RecalculateRadius,						&
									HypotheticalElectrolyteConcentrations 

	USE Chemistry,			 ONLY : HowManyAqChems,					&
									HowManyAqCations,				&
									HowManyAqAnions,				&
									AqMolecularMass,				&
									AqCationMass, AqAnionMass,		&
									AqPhaseChemicalNames,			&
									AqCationNames,AqAnionNames


	USE InfrastructuralCode, ONLY : WARN, ERROR, INT2STR, REAL2STR,	 &
									GetFileHandle, ReturnFileHandle
	
	USE ModelParameters,	 ONLY : OutputDeckSubDir, micron

	IMPLICIT NONE

	!! External Variables
	INTEGER	       :: inHowManyBins
	CHARACTER *(*) :: OutputFileName

	!! Internal Variables
	INTEGER ::  I, Bin, X, Y, Z,					&
				HowManyBins,						&
				HowManyAerosol,						&
				HowManyElectroEquiv,				&
				HowManyContents

	REAL*8  ::	SmallestRadius, LargestRadius,		&
				BinBoundaries(inHowManyBins+1),		&
				RadiusRatio,						&
				NumbHistogram(inHowManyBins),		&
				ElectrolyteEquiv(inHowManyBins,HowManyAqChems), &
				MassHistogram(inHowManyBins,HowManyAqChems+HowManyAqCations+HowManyAqAnions), &
				HypotheticalElectroConcs(HowManyAqChems),	&
				MaxNum, MinNum,						&
				MaxEquiv, MinEquiv,					&
				MaxContent, MinContent

	LOGICAL  :: HypotheticalContent(HowManyAqChems), ActualContent(HowManyAqChems+HowManyAqCations+HowManyAqAnions)

	TYPE(Particle), POINTER :: Current


	HowManyBins = inHowManyBins

	!! Make a first pass to determine the largest and smallest radius
	Current => Particles(1,1,1)%First
	CALL RecalculateRadius(Current)
	SmallestRadius = Current%EffectiveRadius
	LargestRadius = Current%EffectiveRadius

	DO X = 1, XGridPoints(GasPhaseChemGridIndex)
	DO Y = 1, YGridPoints(GasPhaseChemGridIndex)
	DO Z = 1, ZGridPoints(GasPhaseChemGridIndex)

		Current => Particles(X,Y,Z)%First
		
10		IF (ASSOCIATED(Current)) THEN

			HowManyAerosol = HowManyAerosol + Current%NumberOfParticles

			IF (Current%EffectiveRadius .LT. SmallestRadius) SmallestRadius = Current%EffectiveRadius
			IF (Current%EffectiveRadius .GT. LargestRadius)  LargestRadius  = Current%EffectiveRadius

			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 10
			END IF
		END IF
	END DO ; END DO ; END DO


	!! Check to see if the specified number of bins makes sense with the number of particles
	IF (HowManyAerosol/10. .LT. HowManyBins) THEN
		X = MAX (1, FLOOR(HowManyAerosol / 10.))
		CALL WARN("In BinEntireDomainForPrinting(), you asked for "//TRIM(INT2STR(HowManyBins))//		&
				  " size bins.  However, there are only "//TRIM(INT2STR(HowManyAerosol))//				&
				  " aerosol.  Resetting the number of bins to "//TRIM(INT2STR(X))//".")
		HowManyBins = X
	END IF

	!! Get the volume ratio for the bins
	RadiusRatio = (LargestRadius / SmallestRadius) ** (1. / HowManyBins)

	!! And set the boundaries of the radius bins
	BinBoundaries(1) = SmallestRadius
	DO I = 1, HowManyBins
		BinBoundaries(I+1) = BinBoundaries(I) * RadiusRatio 
	END DO
	DO I = 1, HowManyBins+1
		BinBoundaries(I)   = BinBoundaries(I) / micron
	END DO

	!! Initialize histogram outputs
	DO I = 1, HowManyBins

		NumbHistogram(I) = 0
		
		DO X = 1, HowManyAqChems
			ElectrolyteEquiv(I,X) = 0.
		END DO

		DO X = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			MassHistogram(I,X) = 0.
		END DO
	END DO

	!! Tally everything
	DO X = 1, XGridPoints(GasPhaseChemGridIndex)
	DO Y = 1, YGridPoints(GasPhaseChemGridIndex)
	DO Z = 1, ZGridPoints(GasPhaseChemGridIndex)

		Current => Particles(X,Y,Z)%First
		
20		IF (ASSOCIATED(Current)) THEN

			!! The largest radius particle get places out by the inversion, so correct
			Bin = FLOOR (DLOG(Current%EffectiveRadius / SmallestRadius) / DLOG(RadiusRatio)) + 1
			IF (Bin .EQ. HowManyBins + 1) Bin = HowManyBins

			NumbHistogram(Bin) = NumbHistogram(Bin)+ANINT(Current%NumberOfParticles)
			
			!! Hypothetical electrolyte concentrations
			!! (translate to mg / m3)
			HypotheticalElectroConcs = HypotheticalElectrolyteConcentrations (Current)
			DO I = 1, HowManyAqChems
				ElectrolyteEquiv(Bin,I) = ElectrolyteEquiv(Bin,I) + Current%NumberOfParticles*(Current%AqChems(I)	&
										  + HypotheticalElectroConcs(I)) * &
										  AqMolecularMass(I) * 1.e12
			END DO

			!! Actual concentrations including ions
			!! (translate to mg / m3)
			DO I = 1, HowManyAqChems
				MassHistogram(Bin,I) = MassHistogram(Bin,I) + Current%AqChems(I)	&
									   * AqMolecularMass(I) * 1.e12 * Current%NumberOfParticles
			END DO
			DO I = 1, HowManyAqCations
				MassHistogram(Bin,HowManyAqChems+I) =								&
										MassHistogram(Bin,HowManyAqChems+I) +		&
										Current%AqChems(HowManyAqChems+I) *			&
										AqCationMass(I) * 1.e12 * Current%NumberOfParticles
			END DO
			DO I = 1, HowManyAqAnions
				MassHistogram(Bin,HowManyAqChems+HowManyAqCations+I) =							&
										MassHistogram(Bin,HowManyAqChems+HowManyAqCations+I) +	&
										Current%AqChems(HowManyAqChems+HowManyAqCations+I) *	&
										AqAnionMass(I) * 1.e12 * Current%NumberOfParticles
			END DO

			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 20
			END IF
		END IF
	END DO ; END DO ; END DO

	MaxNum		= NumbHistogram(1)
	MinNum		= NumbHistogram(1)
	MaxEquiv	= ElectrolyteEquiv(1,1)
	MinEquiv	= ElectrolyteEquiv(1,1)
	MaxContent	= MassHistogram(1,1)
	MinContent	= MassHistogram(1,1)

	!! Sum everything to stack the plots
	DO X = 1, HowManyBins

		DO Y = 1, HowManyAqChems-1
			DO Z = Y+1, HowManyAqChems
				IF(ElectrolyteEquiv(X,Y) .GT. 0.) ElectrolyteEquiv(X,Y) = ElectrolyteEquiv(X,Y) + ElectrolyteEquiv(X,Z)
		END DO ; END DO 

		!! And general rescaleing
		DO Y = 1, HowManyAqChems
			ElectrolyteEquiv(X,Y) = ElectrolyteEquiv(X,Y) / (DLOG10(BinBoundaries(X+1)) - DLOG10(BinBoundaries(X)))
		END DO

		DO Y = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions - 1
			DO Z = Y+1, HowManyAqChems + HowManyAqCations + HowManyAqAnions 
				IF(MassHistogram(X,Y) .GT. 0.) MassHistogram(X,Y) = MassHistogram(X,Y) + MassHistogram(X,Z)
		END DO ; END DO 

		!! And general rescaleing
		DO Y = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions - 1
			IF (MassHistogram(X,Y) .GT. 0) &
				MassHistogram(X,Y) = MassHistogram(X,Y) / (DLOG10(BinBoundaries(X+1)) - DLOG10(BinBoundaries(X)))
		END DO

		!! Find High and Lows for each type
		IF (MaxNum	   .LT. NumbHistogram(X))		MaxNum		= NumbHistogram(X)
		IF (MinNum	   .GT. NumbHistogram(X)											&
				.AND. NumbHistogram(X) .GT. 0.)		MinNum		= NumbHistogram(X)
		IF (MaxEquiv   .LT. ElectrolyteEquiv(X,1))	MaxEquiv	= ElectrolyteEquiv(X,1)
		IF (MaxContent .LT. MassHistogram(X,1))		MaxContent	= MassHistogram(X,1)

		DO Y = 1, HowManyAqChems
			IF (MinEquiv .GT. ElectrolyteEquiv(X,Y)										&
				.AND. ElectrolyteEquiv(X,Y) .GT. 0.) MinEquiv	= ElectrolyteEquiv(X,Y)
		END DO

		DO Y = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions 
			IF (MinContent .GT. MassHistogram(X,Y)										&
				.AND. MassHistogram(X,Y) .GT. 0.) MinContent	= MassHistogram(X,Y)
		END DO
	END DO

	!! Double the Radii into Diameters
	DO X = 1, HowManyBins+1
		BinBoundaries(X) = BinBoundaries(X)*2.
	END DO

	MaxNum			= 10.**CEILING(DLOG10(MaxNum))
	MinNum			= 10.**FLOOR  (DLOG10(MinNum))
	MaxEquiv		= 10.**CEILING(DLOG10(MaxEquiv))
	MinEquiv		= 10.**FLOOR  (DLOG10(MinEquiv))
	MaxContent		= 10.**CEILING(DLOG10(MaxContent))
	MinContent		= 10.**FLOOR  (DLOG10(MinContent))
	LargestRadius	= 10.**CEILING(DLOG10(LargestRadius*2./micron))
	SmallestRadius	= 10.**FLOOR  (DLOG10(SmallestRadius*2./micron))

	!! Figure out which content is there (don't assume one  
	!! bin is representative, although probably could)
	DO X = 1, HowManyAqChems
		HypotheticalContent(X) = .FALSE.
	END DO

	!! Figure out which species are present
	DO X = 1, HowManyBins
		DO Y = 1, HowManyAqChems
			IF (ElectrolyteEquiv(X,Y) .GT. 0.) HypotheticalContent(Y) = .TRUE.
	END DO ; END DO

	HowManyElectroEquiv = 0
	DO Y = 1, HowManyAqChems
		IF (HypotheticalContent(Y)) HowManyElectroEquiv = HowManyElectroEquiv + 1
	END DO

	DO X = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		ActualContent(X) = .FALSE.
	END DO

	DO X = 1, HowManyBins
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (MassHistogram(X,Y) .GT. 0.) ActualContent(Y) = .TRUE.
	END DO ; END DO

	HowManyContents     = 0
	DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		IF (ActualContent(Y)) HowManyContents = HowManyContents+1
	END DO

	!! Output the tallied data as a MATLAB file
	I = GetFileHandle()
    OPEN(UNIT=I, FILE=TRIM(OutputDeckSubDir)//TRIM(OutputFileName)//".m")

	!! The Number Histogram
	WRITE (I,'(a)') "% Output from MELAM produced by SUBROUTINE BinEntireDomainForPrinting()"
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') "clear MassHistogram ElectrolyteEquiv Radius NumHistogram ;"
	WRITE (I,'(a)') "clear C D Content ;"
	WRITE (I,'(a)') "clear ActualContentNames HypotheticalContentNames;"
	WRITE (I,'(a)') "clf;"
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') ""

	WRITE (I,'(a)') "% Number Histogram Information"
	WRITE (I,'(a)') "% bin number, radius, number conc"

	WRITE (I,'(a)') "NumHistogram = ["

	DO X = 1, HowManyBins

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !//","
		WRITE (I,'(a)')              TRIM(REAL2STR(DBLE(1.e-50)))//"; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !//","
		WRITE (I,'(a)')              TRIM(REAL2STR(MAX(1.e-50,NumbHistogram(X))))//"; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))//" " !", "
		WRITE (I,'(a)')				 TRIM(REAL2STR(MAX(1.e-50,NumbHistogram(X))))//";"

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))//" " !//","
		WRITE (I,'(a)')              TRIM(REAL2STR(DBLE(1.e-50)))//"; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !//","
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(DBLE(1.e-50)))

		IF (X .LT. HowManyBins) THEN
			WRITE (I,'(a)') "; "
		ELSE
			WRITE (I,'(a)') "]; "
		END IF
	END DO

	WRITE (I,'(a)') "subplot(3,1,1)"
	WRITE (I,'(a)') "loglog(NumHistogram(:,2),NumHistogram(:,3),'k');"
	WRITE (I,'(a)') "hold on;"
	WRITE (I,'(a)') "fill(NumHistogram(:,2),NumHistogram(:,3),1);"
	WRITE (I,'(a)') "AXIS(["//TRIM(REAL2STR(SmallestRadius))//" "//TRIM(REAL2STR(LargestRadius))//" "//	&
					TRIM(REAL2STR(MinNum))//" "//TRIM(REAL2STR(MaxNum))//" ])"
	WRITE (I,'(a)') "xlabel('Particle Diameter (microns)');"
	WRITE (I,'(a)') "ylabel('d N / d log D');"
	WRITE (I,'(a)') "title('Number Concentration');"

	WRITE (I,'(a)') "hold off;"

	!! The Hypothetical Electrolyte Content
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') ""

	WRITE (I,'(a)') "% Hypothetical Electrolyte Concentration Histogram Information"
	WRITE (I,'(a)') "% First cite the names:"
	WRITE (I,'(a)') ""


	WRITE (I,'(a)') ""
	WRITE (I,'(a)') "% bin number, radius, Mass Conc of Each (Hypothetical for Electrolytes) Content..."

	WRITE (I,'(a)') "ElectrolyteEquiv = ["

	DO X = 1, HowManyBins

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " ! ", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50"
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " ! ", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " "//TRIM(REAL2STR(MAX(1.e-50,ElectrolyteEquiv(X,Y)),8)) 
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " "//TRIM(REAL2STR(MAX(1.e-50,ElectrolyteEquiv(X,Y)),8)) 
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " ! ", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))//" " !", "
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50"
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " ! ", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50"
		END DO

		IF (X .LT. HowManyBins) THEN
			WRITE (I,'(a)') "; "
		ELSE
			WRITE (I,'(a)') "]; "
		END IF
	END DO

	WRITE (I,'(a)') "subplot(3,1,2)"
	WRITE(I,'(a)',advance='no') "HypotheticalContentNames = {"
	DO X = 1, HowManyAqChems
		IF (HypotheticalContent(X)) THEN
			WRITE(I,'(a)',advance='no') "'"//TRIM(AqPhaseChemicalNames(X))//"' ; "
		END IF
	END DO
	WRITE(I,'(a)') "};"
	WRITE (I,'(a)') "C = [1:"//TRIM(INT2STR(HowManyElectroEquiv))//"];"
	WRITE (I,'(a)') "loglog(ElectrolyteEquiv(:,2),ElectrolyteEquiv(:,3:"//TRIM(INT2STR(2+HowManyElectroEquiv))//"),'k');"
	WRITE (I,'(a)') "AXIS(["//TRIM(REAL2STR(SmallestRadius))//" "//TRIM(REAL2STR(LargestRadius))//" "//TRIM(REAL2STR(MinEquiv))	&
					//" "//TRIM(REAL2STR(MaxEquiv))//" ])"
	WRITE (I,'(a)') "xlabel('Particle Diameter (microns)');"
	WRITE (I,'(a)') "ylabel('d Mass / d log D');"
	WRITE (I,'(a)') "title('Particle Mass Content as Reconstituted Electrolytes');"
	WRITE (I,'(a)') "hold on;"
	WRITE (I,'(a)') "H=fill(ElectrolyteEquiv(:,2),ElectrolyteEquiv(:,3:"//TRIM(INT2STR(2+HowManyElectroEquiv))//"),C);"
	WRITE (I,'(a)') "legend(H,HypotheticalContentNames);"
	WRITE (I,'(a)') "hold off;"

	!! The Actual Content
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') ""

	WRITE (I,'(a)') "% Electrolyte & Other Concentration Histogram Information"
	WRITE (I,'(a)') "% First cite the names:"
	WRITE (I,'(a)') ""

	!! Then the content
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') "% bin number, radius, Mass Conc of Each (Hypothetical for Electrolytes) Content..."

	WRITE (I,'(a)') "Content = ["

	DO X = 1, HowManyBins
		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50 "
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " "//TRIM(REAL2STR(MAX(1.e-50,MassHistogram(X,Y)),8)) !", "
		END DO
		WRITE (I,'(a)') "; "


		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " "//TRIM(REAL2STR(MAX(1.e-50,MassHistogram(X,Y)),8)) !", "
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))//" " !", "
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50 "
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50 "
		END DO

		IF (X .LT. HowManyBins) THEN
			WRITE (I,'(a)') "; "
		ELSE
			WRITE (I,'(a)') ""
			WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
			WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))
			DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
				IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50 "
			END DO
			WRITE (I,'(a)') "]; "
		END IF
	END DO

	WRITE (I,'(a)') "subplot(3,1,3)"
	!! Put out names
	WRITE(I,'(a)',advance='no') "ActualContentNames = {"
	DO X = 1, HowManyAqChems
		IF (ActualContent(X)) WRITE(I,'(a)',advance='no') "'"//TRIM(AqPhaseChemicalNames(X))//"'; "
	END DO
	DO X = 1, HowManyAqCations
		IF (ActualContent(HowManyAqChems+X)) WRITE(I,'(a)',advance='no') "'"//TRIM(AqCationNames(X))//"'; "
	END DO
	DO X = 1, HowManyAqAnions
		IF (ActualContent(HowManyAqChems+HowManyAqCations+X)) WRITE(I,'(a)',advance='no') "'"//TRIM(AqAnionNames(X))//"'; "
	END DO
	WRITE (I,'(a)') "};"
	WRITE (I,'(a)') "D = [1:"//TRIM(INT2STR(HowManyContents))//"];"
	WRITE (I,'(a)') "loglog(Content(:,2),Content(:,3:"//TRIM(INT2STR(2+HowManyContents))//"),'k');"
	WRITE (I,'(a)') "AXIS(["//TRIM(REAL2STR(SmallestRadius))//" "//TRIM(REAL2STR(LargestRadius))//" "//		&
					TRIM(REAL2STR(MinContent))//" "//TRIM(REAL2STR(MaxContent))//" ])"
	WRITE (I,'(a)') "xlabel('Particle Diameter (microns)');"
	WRITE (I,'(a)') "ylabel('d Mass / d log D');"
	WRITE (I,'(a)') "title('Particle Mass Content');"
	WRITE (I,'(a)') "hold on;"
	WRITE (I,'(a)') "H = fill(Content(:,2),Content(:,3:"//TRIM(INT2STR(2+HowManyContents))//"),D);"
	WRITE (I,'(a)') "legend(H,ActualContentNames);"
	WRITE (I,'(a)') "hold off;"

	!! Close the file
	CLOSE(I)
	CALL ReturnFileHandle(I)
	
	RETURN
END SUBROUTINE BinEntireDomainForPrinting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Splash everything into one bin by constituent !!
!! and print in text file, vaguely mirroring     !!
!! ISORROPIA's output type.                      !!
!!												 !!!!!!
!! FOR REPEATED CALLS, ONLY CALL AFTER EQUILIBRATION !!
!! OR IT WILL INCREASE THE NUMBER OF COLUMNS FROM ONE!!
!! CALL TO ANOTHER.									 !!
SUBROUTINE SummarizeDomainMassForOutput (FileName, FirstCall)

	USE GridPointFields,	ONLY :	XGridPoints,			&
									YGridPoints,			&
									ZGridPoints,			&
									GasPhaseChemGridIndex,	&
									GetGasChemVector,		&
									GetRelativeHumidity,	&
									GetTemp,				&
									LocateGridPoint

	USE Aerosols,           ONLY :	Particles, Particle, &
									ParticleArray

	USE Chemistry,			ONLY :  HowManyAqChems,					&
									HowManyAqCations,				&
									HowManyAqAnions,				&
									AqMolecularMass,				&
									GasMolecularMass,				&
									AqCationMass, AqAnionMass,		&
									GasPhaseChemicalNames,			&
									AqPhaseChemicalNames,			&
									AqCationNames,AqAnionNames,		&
									HowManyGasChems

	USE ModelParameters,    ONLY :  Avogadro, OutputDeckSubDir

	USE InfrastructuralCode, ONLY : GetFileHandle,					&
									ReturnFileHandle,				&
									ERROR,							&
									REAL2STR

	IMPLICIT NONE

	!! External Variables
	!INTEGER,OPTIONAL         :: FileHandle
	CHARACTER *(*), OPTIONAL :: FileName
	LOGICAL, OPTIONAL :: FirstCall

	!! Internal Variables
	INTEGER :: I, X, Y, Z, FH
	REAL*8  :: AerosolMass(HowManyAqChems+HowManyAqCations+HowManyAqAnions), &
			   GasMass(HowManyGasChems), GridPointGasVec(HowManyGasChems),	 &
			   RH, Temperature, NumbAerosol, XX(3)

	TYPE(Particle), POINTER :: Current


	!! Error check the inputs
	!IF (.NOT.PRESENT(FileName) .AND. .NOT.PRESENT(FileHandle) .OR.  &
	!         PRESENT(FileName) .AND.      PRESENT(FileHandle))      &
	!	CALL ERROR("Error in call to SummarizeDomainMassForOutput(), you must supply either a file handle or a file name "//	&
	!			   "and supplied either neither or both.  Please check the syntax.")

	!IF (PRESENT(FileHandle)) FH = FileHandle
	IF (PRESENT(FileName)) THEN
		FH = GetFileHandle()
		OPEN(UNIT=FH, FILE=TRIM(OutputDeckSubDir)//TRIM(FileName)//".txt")
	END IF

	!! Initialize
	DO I = 1, HowManyGasChems
		GasMass(I) = 0.
	END DO

	DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		AerosolMass(I) = 0.
	END DO

	NumbAerosol = 0.
	RH          = 0.
	Temperature = 0.

	!! Tally everything
	DO X = 1, XGridPoints(GasPhaseChemGridIndex)
	DO Y = 1, YGridPoints(GasPhaseChemGridIndex)
	DO Z = 1, ZGridPoints(GasPhaseChemGridIndex)

		GridPointGasVec = GetGasChemVector (X, Y, Z)

		XX          = LocateGridPoint(X,Y,Z,GasPhaseChemGridIndex)
		RH          = RH + GetRelativeHumidity(XX(1),XX(2),XX(3))
		Temperature = Temperature + GetTemp(XX(1),XX(2),XX(3))

		DO I = 1, HowManyGasChems
			GasMass(I) = GasMass(I) + GridPointGasVec(I) &
					     * GasMolecularMass(I) &
						 * 1.e12 / Avogadro
		END DO

		Current => Particles(X,Y,Z)%First
10		IF (ASSOCIATED(Current)) THEN

			NumbAerosol = NumbAerosol + Current%NumberOfParticles
			DO I = 1, HowManyAqChems
				AerosolMass(I) =	AerosolMass(I) + Current%AqChems(I)	&
									* AqMolecularMass(I) * 1.e12 * Current%NumberOfParticles
			END DO
			DO I = 1, HowManyAqCations
				AerosolMass(HowManyAqChems+I) =							&
									AerosolMass(HowManyAqChems+I) +		&
									Current%AqChems(HowManyAqChems+I) *	&
									AqCationMass(I) * 1.e12 * Current%NumberOfParticles
			END DO
			DO I = 1, HowManyAqAnions
				AerosolMass(HowManyAqChems+HowManyAqCations+I) =						 &
									AerosolMass(HowManyAqChems+HowManyAqCations+I) +	 &
									Current%AqChems(HowManyAqChems+HowManyAqCations+I) * &
									AqAnionMass(I) * 1.e12 * Current%NumberOfParticles
			END DO

			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 10
			END IF
		END IF
	END DO ; END DO ; END DO

	!! Scale RH, Temp, and Number Conc
	RH          = RH / XGridPoints(GasPhaseChemGridIndex) / YGridPoints(GasPhaseChemGridIndex) / ZGridPoints(GasPhaseChemGridIndex)
	Temperature = Temperature / XGridPoints(GasPhaseChemGridIndex) / YGridPoints(GasPhaseChemGridIndex) / &
				  ZGridPoints(GasPhaseChemGridIndex)
	NumbAerosol = NumbAerosol / XGridPoints(GasPhaseChemGridIndex) / YGridPoints(GasPhaseChemGridIndex) / &
					ZGridPoints(GasPhaseChemGridIndex) 

	!! Now You Have all of the data.  
	!! Output the tallied data as an EXCEL type comma delimited file.
	IF (PRESENT(FirstCall) .AND. FirstCall .EQ. .TRUE.) THEN

		!! RH 
		WRITE (FH,'(a)',advance='no') "Mean RH, "

		!! Temperature
		WRITE (FH,'(a)',advance='no') "Mean Temperature, "

		!! Number content
		WRITE (FH,'(a)',advance='no') "Aerosol Numb (n/cm3), "
		
		!! Print names of content
		DO I = 1, HowManyGasChems
			IF (GasMass(I) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(GasPhaseChemicalNames(I))//" (g) (microg / m3), "
		END DO

		IF (AerosolMass(1) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(AqPhaseChemicalNames(1))//" (l) (microg / m3), "

		DO I = 2, HowManyAqChems
			IF (AerosolMass(I) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(AqPhaseChemicalNames(I))//" (aq) (microg / m3), "
		END DO
		DO I = 1, HowManyAqCations
			IF (AerosolMass(HowManyAqChems+I) .GT. 0.) &
				WRITE (FH,'(a)',advance='no') TRIM(AqCationNames(I))//" (aq) (microg / m3), "
		END DO
		DO I = 1, HowManyAqAnions
			IF (AerosolMass(HowManyAqChems+HowManyAqCations+I) .GT. 0.) &
				WRITE (FH,'(a)',advance='no') TRIM(AqAnionNames(I))//" (aq) (microg / m3), "
		END DO

		WRITE (FH,'(a)') ""
	END IF

	!! Report the environmental variables
	WRITE (FH,'(a)',advance='no') &
		TRIM(REAL2STR(RH,8))//", "//TRIM(REAL2STR(Temperature,8))//", "//TRIM(REAL2STR(NumbAerosol,8))//", "

	!! Report the current concentrations / mass loadings
	DO I = 1, HowManyGasChems
		IF (GasMass(I) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(REAL2STR(GasMass(I),8))//", "
	END DO
	DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		IF (AerosolMass(I) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(REAL2STR(AerosolMass(I),8))//", "
	END DO
	WRITE (FH,'(a)') ""

	!! Close the file if we're controlling it locally
	IF (PRESENT(FileName)) THEN
		CLOSE(FH)
		CALL ReturnFileHandle(FH)
	END IF

END SUBROUTINE SummarizeDomainMassForOutput


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Take Lagrangian Particles and put them in set of !!
!! bins for outputting                              !!
!!                                                  !!
!! Endpoints are treated just as they are for the   !!
!! sectional definitions: the lowest bin is from    !!
!! zero to an endpoint and the largest is from an   !!
!! endpoint to infinity.                            !!
SUBROUTINE BinLagrangianParticlesForOutput (XP,YP,ZP,NumberOfBins,SmallEdge,LargeEdge, FH, FirstComment)

	USE InfrastructuralCode, ONLY : ERROR, REAL2STR
	USE Aerosols,            ONLY : Particles, Particle, RecalculateRadius
	USE ModelParameters, ONLY : LagrangianOrNot,km

	IMPLICIT NONE

	INTEGER :: XP, YP, ZP, NumberOfBins,FH
	REAL*8  :: SmallEdge, LargeEdge, BinRadiusRatio
	REAL*8  :: BinEdges(NumberOfBins,2),BinVol(NumberOfBins),BinNumb(NumberOfBins)
	INTEGER :: I
	TYPE(Particle), POINTER :: Runner
	CHARACTER *(*) :: FirstComment

	!! Exception Handle on Settings
	IF (.NOT.LagrangianOrNot) &
		CALL ERROR ("Called SCBinLagrangianParticlesForOutput() even though we are not tracking lagrangian " // &
					"particles on this run.  This is not allowed.")

	IF (NumberOfBins .LE. 1) &
		CALL ERROR ("SCBinLagrangianParticlesForOutput() is only to be used for number distributions of aerosol " // &
				    "(at least two bins).  You've asked for one or fewer bins to be calculated.  We don't accommodate that.")

	IF (NumberOfBins .EQ. 2 .AND. SmallEdge .NE. LargeEdge) &
		CALL ERROR ("If you want to use only two bins in SCBinLagrangianParticlesForOutput(), you have to set the small " // &
					"bin boundary and the large bin boundary equal.  The structure will go from zero to this dividor and " // &
					"then from that dividor to infinity (radius).")

	DO I = 1, NumberOfBins
		BinVol(I) = 0.
		BinNumb(I) = 0.
	END DO


	BinEdges(1,1) = 0.
	BinEdges(1,2) = SmallEdge

	IF (NumberOfBins .EQ. 2) THEN
		!! Calculate the Bin Edges
		BinEdges(2,1) = SmallEdge
		BinEdges(2,2) = 1.*km
	ELSE
		!! Calculate the radius ratio
		BinRadiusRatio = (LargeEdge/SmallEdge)**(1./(NumberOfBins-2.))
		DO I = 2, NumberOfBins-1
			BinEdges(I,1) = BinEdges(I-1,2)
			BinEdges(I,2) = BinEdges(I,1) * BinRadiusRatio
		END DO

		BinEdges(NumberOfBins-1,2) = LargeEdge !! just so this is exactly as asked with no rounding error
		BinEdges(NumberOfBins,1)   = LargeEdge
		BinEdges(NumberOfBins,2)   = 1.*km
	END IF


	!! Count up the bins
	Runner => Particles(XP,YP,ZP)%First

100	IF (ASSOCIATED(Runner)) THEN

		CALL RecalculateRadius(Runner)

		DO I = 1, NumberOfBins
			IF (Runner%EffectiveRadius .GE. BinEdges(I,1) .AND. &
				Runner%EffectiveRadius .LT. BinEdges(I,2)) EXIT
		END DO

		BinVol(I)  = BinVol(I)  + Runner%EffectiveRadius**3.
		BinNumb(I) = BinNumb(I) + Runner%NumberOfParticles

	END IF
	
	IF (ASSOCIATED(Runner%Next)) THEN
		Runner => Runner%Next
		GOTO 100
	END IF

	WRITE (FH,'(a)',advance='no') TRIM(FirstComment)

	DO I = 1, NumberOfBins
		WRITE (FH,'(a)',advance='no') ","//TRIM(REAL2STR(BinNumb(I)))
	END DO

	DO I = 1, NumberOfBins
		IF (BinNumb(I) .GT. 0.) THEN
			BinVol(I) = (BinVol(I)/BinNumb(I))**(1./3.)
		ELSE
			BinVol(I) = (BinEdges(I,1)+BinEdges(I,2))/2.
		END IF

		WRITE (FH,'(a)',advance='no') ","//TRIM(REAL2STR(BinVol(I)))

	END DO

	WRITE (FH,'(a)') ""


	RETURN
END SUBROUTINE BinLagrangianParticlesForOutput

SUBROUTINE SpillBeans(XP, YP, ZP, Env)
	
	USE InfrastructuralCode, ONLY : ERROR, REAL2STR, INT2STR
	USE Aerosols,            ONLY : Particles, &
									Particle, &
									RecalculateRadius, &
									AerosolPh, &
									CurvatureCorrection, &
									OrgCurvatureCorrection, &
									BoundaryParticles, &
									HypotheticalElectrolyteConcentrations
	USE Chemistry,			 ONLY : AqPhaseChemicalNames, &
									OrgPhaseChemicalNames, &
									HowManyOrgChems, &
									HowManyAqOrgChems, &
									AqOrgPhaseChemicalNames, &
									AqMolecularMass, &
									OrgMolecularMass, &
									AqOrgMolecularMass, &
									HowManyAqChems, FindChem, &
									GasPhaseBackground
	USE ModelParameters,	 ONLY : micron 
	USE GridPointFields, ONLY : GasPhaseChemGridIndex, LocateGridPoint, &
								GetAirDensity, GetM, &
								GridGasChem

	Implicit None

	!External
	INTEGER		::	XP, YP, ZP
	LOGICAL, OPTIONAL :: Env !True if we want the boundary distribution

	!Internal
	INTEGER		:: I, NumAqChems, NumOrgChems, A
	REAL :: XX(3), Dens, ppbmCO, Mgas
	REAL*8 :: SaltConc(HowManyAqChems)
	
	TYPE(Particle), POINTER :: CurrentParticle
			
	! spill the beans...
	IF (NOT(PRESENT(Env))) THEN
		CurrentParticle => particles(XP,YP,ZP)%first
	ELSE
		IF (Env) THEN
			CurrentParticle => BoundaryParticles%first
		ELSE
			CurrentParticle => particles(XP,YP,ZP)%first
		END IF
	END IF

	XX = LocateGridPoint (XP, YP, ZP, GasPhaseChemGridIndex)
	Dens = GetAirDensity(XX(1),XX(2),XX(3)) !g/cm3
	Mgas = GetM   (XX(1),XX(2),XX(3))
	WRITE(*,*) "Dens: ", Dens
	ppbmCO = (GridGasChem(1,1,1,FindChem("CO",0)) - &
	          GasPhaseBackground(FindChem("CO",0)))*(1e9/Mgas)*(28/28.973)
	WRITE(*,*) "ppbm CO: ", ppbmCO
			
	NumAqChems = size(currentparticle%AqChems)
	NumOrgChems = HowManyOrgChems
	i = 1
	99 if (associated(currentparticle)) then
		print *, ""
		print *, "Number ",trim(INT2STR(I)),",  ",currentparticle%numberofparticles,currentparticle%sectional
		print *, "Effective Radius: ",trim(real2str(CurrentParticle%EffectiveRadius/micron)),",  (",trim(real2str(CurrentParticle%Edges(1)/micron))," to ",trim(real2str(CurrentParticle%Edges(2)/micron)),")"
		print *, "Embryo Radius: ",trim(real2str(CurrentParticle%EmbryoRadius/micron))
		print *, "Insoluble Radius: ",trim(real2str(CurrentParticle%InsolubleRadius/micron))
		print *, "Absorbing Core Radius: ",trim(real2str(CurrentParticle%AbsCoreRad/micron))
		print *, "Shell Real Ref. Ind. 550 nm",trim(real2str(currentparticle%ShellRealRefracInd(1)))
		print *, "Shell Imag. Ref. Ind. 550 nm",trim(real2str(currentparticle%ShellImagRefracInd(1)))
		print *, "Extinction Efficiency: 550 nm",trim(real2str(CurrentParticle%ExtEff(1)))
		print *, "Scattering Efficiency: 550 nm",trim(real2str(currentparticle%ScaEff(1)))
		print *, "BackScattering Efficiency: 550 nm",trim(real2str(currentparticle%BackScaEff(1)))
		print *, "Assymetry Parameter: 550 nm",trim(real2str(currentparticle%AssymParam(1)))

	
		print *, "Temperature ",trim(real2str(currentparticle%Temperature))
		print *, "Ionic Strength ",trim(real2str(currentparticle%IonicStr))
		print *, "Water Activity ",trim(real2str(currentparticle%WaterActivity))
		print *, "SolutionDensity ",trim(real2str(currentparticle%SolutionDensity))
		!print *, "InsolubleDensity ",trim(real2str(currentparticle%InsolubleDensity))
		print *, "ParticleDensity ",trim(real2str(currentparticle%ParticleDensity))
		print *, "Ph ", trim(real2str(AerosolPh(currentparticle)))
		print *, "Curvature Correction ", trim(real2str(CurvatureCorrection(currentparticle)))
		print *, "OrgCurvature Correction ", trim(real2str(OrgCurvatureCorrection(currentparticle)))	
		print *, "Solution Surface Tension", trim(real2str(currentparticle%SurfaceTension))
		
		print *, "Temperature ", trim(real2str(currentparticle%temperature))
		
		SaltConc = HypotheticalElectrolyteConcentrations (currentparticle)
		DO A = 1, NumAqChems		
			IF(currentparticle%aqchems(A) .GT. 1e-40) THEN
				write(*,*) trim(int2str(A))," : ",AqPhaseChemicalNames(A) ," (mol/cm3 air): ", trim(real2str(currentparticle%aqchems(A)*currentparticle%numberofparticles))
				!IF(A .LE. HowManyAqChems) THEN
				!write(*,*) trim(int2str(A))," : ",AqPhaseChemicalNames(A) ," (ppbm): ", trim(real2str((SaltConc(A)+currentparticle%aqchems(A))*currentparticle%numberofparticles* &
				!																						AqMolecularMass(A)*1.0e9/Dens))
				!END IF
				!IF(A .LE. HowManyAqChems) THEN
				!write(*,*) trim(int2str(A))," : ",AqPhaseChemicalNames(A) ," (ppbm/ppbm CO): ", trim(real2str((SaltConc(A)+currentparticle%aqchems(A))*currentparticle%numberofparticles* &
				!																						AqMolecularMass(A)*1.0e9/Dens/ppbmCO))
				!END IF
			END IF
		END DO
		
		DO A = 1, NumOrgChems		
			IF(currentparticle%Orgchems(A) .GT. 1e-40) THEN
				write(*,*) trim(int2str(A))," : ",OrgPhaseChemicalNames(A) ," (mol/cm3 air): ", trim(real2str(currentparticle%Orgchems(A)*currentparticle%numberofparticles))
				!write(*,*) trim(int2str(A))," : ",OrgPhaseChemicalNames(A) ," (ppbm): ", trim(real2str(currentparticle%Orgchems(A)*currentparticle%numberofparticles* &
				!																							OrgMolecularMass(A)*1.0e9/Dens))
				!write(*,*) trim(int2str(A))," : ",OrgPhaseChemicalNames(A) ," (ppbm/ppbm CO): ", trim(real2str(currentparticle%Orgchems(A)*currentparticle%numberofparticles* &
				!																							OrgMolecularMass(A)*1.0e9/Dens/ppbmCO))

			END IF
		END DO
		DO A = 1, HowManyAqOrgChems		
			IF(currentparticle%AqOrgchems(A) .GT. 1e-40) THEN
				write(*,*) trim(int2str(A))," : ",AqOrgPhaseChemicalNames(A) ," (mol/cm3 air): ", trim(real2str(currentparticle%AqOrgchems(A)*currentparticle%numberofparticles))
				!write(*,*) trim(int2str(A))," : ",AqOrgPhaseChemicalNames(A) ," (ppbm): ", trim(real2str(currentparticle%AqOrgchems(A)*currentparticle%numberofparticles* &
				!																							AqOrgMolecularMass(A)*1.0e9/Dens))
				!write(*,*) trim(int2str(A))," : ",AqOrgPhaseChemicalNames(A) ," (ppbm/ppbm CO): ", trim(real2str(currentparticle%AqOrgchems(A)*currentparticle%numberofparticles* &
				!																							AqOrgMolecularMass(A)*1.0e9/Dens/ppbmCO))
			END IF
		END DO
		I = I + 1

		if (associated(currentparticle%next)) then
			currentparticle => currentparticle%next
			goto 99
		end if
	end if !Spill Beans

END SUBROUTINE SpillBeans

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine outputs a file at a given time step that                     !!
!! Contains the information necessary to plot the mass and number distributions. !!
!! First written by Matt Alvarado, 10/19/05.                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SizeDistOutput(XP, YP, ZP, t, Filename)
	
	USE InfrastructuralCode, ONLY : ERROR, REAL2STR, INT2STR, GetFileHandle
	USE Aerosols,            ONLY : Particles, &
									Particle, &
									RecalculateRadius, &
									AerosolPh
	USE Chemistry,			 ONLY : AqPhaseChemicalNames, &
									OrgPhaseChemicalNames, &
									GasPhaseChemicalNames, &
									HowManyOrgChems, &
									FindChem, &
									OrgMolecularMass, &
									AqMolecularMass, &
									HowManyAqOrgChems, &
									AqOrgMolecularMass, &
									AqOrgNumCarbon
	USE ModelParameters,	 ONLY : micron, OutputDeckSubDir, &
									Monodisperse, ThermoBulkMode

	USE GridPointFields,	 ONLY : GridGasChem, &
									GetM, &
									LocateGridPoint, &
									GasPhaseChemGridIndex
	
	Implicit None

	!External
	INTEGER		::	XP, YP, ZP
	REAL*8      :: t
	CHARACTER *(*), OPTIONAL :: FileName

	!Internal
	INTEGER		:: I, NumAqChems, NumOrgChems, A, J

	TYPE(Particle), POINTER :: Current
	Integer :: 	KNO3i, KCli, KHSO4i, &
				Nai, NaNO3i, NaCli, NaHSO4i, Na2SO4i, &
				NH3i, NH4i, NH4NO3i, NH4Cli, NH4HSO4i, NH42SO4i, &
				Cli, NO3i, SO4i, HSO4i, K2SO4i, PotassiumIndex, &
				K2SO4Index, FH, NumBins	, COi, Traceri
	REAL*8 :: TotalPOA, TotalSOA, TotalK, TotalNa, TotalNH4, TotalCl, &
			  TotalNO3, TotalSO4, DeltaLogRad, X(3), Mgas, TotalOrgCar		
					
	
	!Find chemical indices for later
	PotassiumIndex = FindChem("K+", 1)
	KNO3i = FindChem("KNO3", 1)
	KCli = FindChem("KCl", 1)
	KHSO4i = FindChem("KHSO4", 1)
	K2SO4Index = FindChem("K2SO4", 1)
	K2SO4i = K2SO4Index
		
	Nai = FindChem("Na+", 1)
	NaNO3i = FindChem("NaNO3", 1)
	NaCli = FindChem("NaCl", 1)
	NaHSO4i = FindChem("NaHSO4", 1)
	Na2SO4i = FindChem("Na2SO4", 1)
		
	NH3i = FindChem("NH3", 1)
	NH4i = FindChem("NH4+", 1)
	NH4NO3i = FindChem("NH4NO3", 1)
	NH4Cli = FindChem("NH4Cl", 1)
	NH4HSO4i = FindChem("NH4HSO4", 1)
	NH42SO4i = FindChem("(NH4)2SO4", 1)

	Cli = FindChem("Cl-", 1)
	NO3i = FindChem("NO3-", 1)
	SO4i = FindChem("SO4--", 1)
	HSO4i = FindChem("HSO4-", 1)

	Traceri = FindChem("Trace", 0)
	COi = FindChem("CO", 0)

	!Count # of bins, # of Aq. Chems
	Current => particles(1,1,1)%first
	NumBins = 0
	I = 1
	7 if (associated(Current)) then
		NumBins = NumBins + 1
		NumAqChems = size(Current%AqChems)
		NumOrgChems = size(Current%OrgChems)
		IF (I .EQ. 2) THEN
			DeltaLogRad = log10(Current%Edges(2)*1.0e4) - log10(Current%Edges(1)*1.0e4)
		END IF
			
		if (associated(Current%next)) then
			Current => Current%next
			I = I + 1
			goto 7
		end if
	end if

	FH = GetFileHandle()
	
	OPEN(UNIT = FH, FILE=TRIM(OutputDeckSubDir)//TRIM(Filename), STATUS="REPLACE")
	WRITE(FH, FMT='(A4,A11)', ADVANCE='NO') "Time", " "
	WRITE(FH, FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(Traceri)
	WRITE(FH, FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(COi)

	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "Edge 1", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "Edge 2", " "
	
	WRITE(FH,FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "Number", " "

	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "OrgCar", " "
			
	WRITE(FH,FMT='(A5, A8)', ADVANCE='NO') "TotBC", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotPOA", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotSOA", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotH2O", " "
	WRITE(FH,FMT='(A4, A9)', ADVANCE='NO') "TotK", " "
	WRITE(FH,FMT='(A5, A8)', ADVANCE='NO') "TotNa", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotNH4", " "
	WRITE(FH,FMT='(A5, A8)', ADVANCE='NO') "TotCl", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotNO3", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotSO4", " "
			
	WRITE(FH, FMT='(A1)') ""
	
	X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
	Mgas = GetM   (X(1),X(2),X(3))
	
	Current => particles(1,1,1)%first
		DO I = 1, NumBins
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') t/60
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') GridGasChem(1,1,1,Traceri)*1e9/Mgas	
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') GridGasChem(1,1,1,COi)*1e9/Mgas
			
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%Edges(1)*1.0e4
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%Edges(2)*1.0e4
			
			IF(Current%effectiveradius*1.0e4 .GT. 0.0) THEN
				WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%effectiveradius*1.0e4
			ELSE
				IF (I .EQ. NumBins) THEN			
					IF(Monodisperse .EQ. .FALSE. .AND. ThermoBulkMode .EQ. .FALSE.) THEN
						WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') 10**(log10(Current%Edges(1)*1.0e4)+0.5*DeltaLogRad)
					ELSE
						WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') 0.0
					END IF
				ELSE
					IF(Monodisperse .EQ. .FALSE..AND. ThermoBulkMode .EQ. .FALSE.) THEN
						WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') 10**(log10(Current%Edges(2)*1.0e4)-0.5*DeltaLogRad)
					ELSE
						WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') 0.0
					END IF
				END IF
			END IF
			
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%numberofparticles
						
			TotalOrgCar = 0.0
			DO J = 2, NumOrgChems
				TotalOrgCar = TotalOrgCar + Current%AqOrgChems(J-1)*AqOrgNumCarbon(J-1)*Current%numberofparticles*1.0e12*12.0
				TotalOrgCar = TotalOrgCar + Current%OrgChems(J)*AqOrgNumCarbon(J-1)*Current%numberofparticles*1.0e12*12.0
			END DO
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalOrgCar

			!BC Mass Conc (ug/m3)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%Orgchems(1)*Current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
			!POA Mass Conc (ug/m3)
			TotalPOA = 0.
			DO J = 2, 9
				TotalPOA = TotalPOA + Current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalPOA = TotalPOA + Current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			DO J = 20, NumOrgChems
				TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

			!SOA Mass Conc (ug/m3)
			TotalSOA = 0.
			DO J = 10, 19
				TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalSOA = TotalSOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
			!Water Mass Conc. (ug/m3)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%Aqchems(1)*Current%numberofparticles*1.0e12*AqMolecularMass(1)
					
			!K+ Mass Conc (ug/m3)
			TotalK = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
					TotalK = TotalK + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				ELSE IF (J .EQ. K2SO4Index) THEN
					TotalK = TotalK + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalK = TotalK*AqMolecularMass(PotassiumIndex)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
			!Na+ Mass Conc (ug/m3)
			TotalNa = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
					TotalNa = TotalNa + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				ELSE IF (J .EQ. Na2SO4i) THEN
					TotalNa = TotalNa + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNa = TotalNa*AqMolecularMass(Nai)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
			!NH4+ Mass Conc (ug/m3)
			TotalNH4 = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
					TotalNH4 = TotalNH4 + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				ELSE IF (J .EQ. NH42SO4i) THEN
					TotalNH4 = TotalNH4 + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

			!Cl- Mass Conc (ug/m3)
			TotalCl = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
						TotalCl = TotalCl + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalCl = TotalCl*AqMolecularMass(Cli)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

			!NO3- Mass Conc (ug/m3)
			TotalNO3 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
					TotalNO3 = TotalNO3 + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

			!SO4- Mass Conc (ug/m3)
			TotalSO4 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
				    .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
					TotalSO4 = TotalSO4 + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

			
			WRITE(FH, FMT='(A1)') ""
					
			IF (associated(Current%next)) THEN
				Current => Current%next
			ELSE
				CYCLE
			END IF

		END DO

		CLOSE(UNIT = FH)


END SUBROUTINE SizeDistOutput

END MODULE OutputRoutines

