!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ReadInputDeck.f90
!! This file holds parameters and constants germaine to the 
!! aerosol definitions and internal processes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES						        !!
!! 1. Chemistry								!!
!! 2. InfrastructuralCode					!!
!! 3. GridPointFields						!!				
!! 4. ModelParameters						!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!	
!! 1. SUBROUTINE ReadMainInputDeck()								!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! READ INPUT DECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Read the main input deck and set the appropriate !!
	!! parameters throughout the model.					!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ReadMainInputDeck ()

		USE Chemistry,			 ONLY : SetPhaseIfs
		USE InfrastructuralCode, ONLY : GetFileHandle,		&
										ReturnFileHandle,	&
										GetLine,			&
										GetToken,			&
										IsInteger,			&
										IsReal,				&
										STR2INT,			&
										STR2REAL,			&
										StripToken,			&
										EOF,				&
										ERROR,WARN,			&
										Transcript
		USE GridPointFields,     ONLY : InitializeGridPointFields,	&
										SetWaterVaporField 
		USE ModelParameters

		IMPLICIT NONE

		INTEGER   :: I,FH,FHT
		LOGICAL   :: IFGAS,IFAQ,IFSOLID
		CHARACTER(len=512) :: Token
		CHARACTER(len=512)::CurrLine,StoreLine
		REAL*8    :: T,P,II, RH, BRH 

		FHT = GetFileHandle()
		FH  = GetFileHandle()
	        OPEN(UNIT=FH, FILE='InorganicInputDecks/ASPInputDeck.in', STATUS='OLD')

		!!! -- Specify Input Deck Location
		CALL SetInputDeckSubDir("InorganicInputDecks/")
		
		!!! -- Specify Output Deck Location
		CALL SetOutputDeckSubDir("OutputFiles/")
                CALL SetTranscriptFH(FHT)
		OPEN(UNIT=FHT, FILE=TRIM("OutputFiles/Transcript.TXT"))
		

		!! -- Get TEMPERATURE
		Token = GetLine(FH)

		!! Assume it is an input deck if not a real
		IF (IsReal(Token) .EQ. 0) THEN
			CALL ERROR("In ASPInputDeck.in, You tried to load Temperature from an input "//&
					   "deck.  As much as I'd like to oblige, that isn't implemented yet.")
		ELSE
			T = STR2REAL(Token)
		END IF

		!! -- Get PRESSURE
		Token = GetLine(FH)

		!! Assume it is an input deck if not a real
		IF (IsReal(Token) .EQ. 0) THEN
			CALL ERROR("In ASPInputDeck.in, You tried to load Pressure from an input "//&
					   "deck.  As much as I'd like to oblige, that isn't implemented yet.")
		ELSE
			P = STR2REAL(Token)
		END IF

		!! -- Get INITIAL RELATIVE HUMIDITY
		Token = GetLine(FH)

		!! Assume it is an input deck if not a real
		IF (IsReal(Token) .EQ. 0) THEN
			CALL ERROR("In ASPInputDeck.in, You tried to load RH from an input "//&
					   "deck.  As much as I'd like to oblige, that isn't implemented yet.")
		ELSE
			II = STR2REAL(Token)
			RH = II
		END IF

		!! -- Get BACKGROUND RELATIVE HUMIDITY
		Token = GetLine(FH)

		!! Assume it is an input deck if not a real
		IF (IsReal(Token) .EQ. 0) THEN
			CALL ERROR("In ASPInputDeck.in, You tried to load RH from an input "//&
					   "deck.  As much as I'd like to oblige, that isn't implemented yet.")
		ELSE
			II = STR2REAL(Token)
			BRH = II
		END IF

		!! Set the grid structure, t and p
		CALL InitializeGridPointFields (T,P)
		CALL SetWaterVaporField (RH)
		CALL SetBackgroundRH(BRH)
		
		!!! -- Specify Which Chemical Phases are to be included
		CurrLine = GetLine(FH)

			!!! Gas Phase
			CALL GetToken(CurrLine,";",Token)

			IF (IsInteger(Token) .EQ. 0) &
				CALL ERROR("Must Supply an Integer in ASPInputDeck.in to indicate "//&
						   "whether Gas Phase Chemistry Should be Included")
			I = STR2INT(Token)
			IF (I < 0 .OR. I > 1) &
				CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
						   "whether Gas Phase Chemistry Should be Included.  Received "//StripToken(Token))
			IF (I .EQ. 0) THEN
				IFGAS = .FALSE.
				CALL SetGasChemistryFlag (.FALSE.)
			ELSE
				IFGAS = .TRUE. 
				CALL SetGasChemistryFlag (.TRUE.)
			END IF

			!!! Aqueous Phase
			CALL GetToken(CurrLine,";",Token)

			IF (IsInteger(Token) .EQ. 0) &
				CALL ERROR("Must Supply an Integer in ASPInputDeck.in to indicate"//&
						   "whether Aqueous Phase Chemistry Should be Included")
			I = STR2INT(Token)
			IF (I < 0 .OR. I > 1) &
				CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
						   "whether Aqueous Phase Chemistry Should be Included.  Received "//StripToken(Token))
			IF (I .EQ. 0) THEN
				IFAQ = .FALSE.
				CALL SetAqChemistryFlag (.FALSE.)
			ELSE
				IFAQ = .TRUE.
				CALL SetAqChemistryFlag (.TRUE.)
			END IF

		CALL SetPhaseIfs(IFGAS,IFAQ)
		
		!!! Poly or Mono Disperse?
		Token = GetLine(FH)

		IF (IsInteger(Token) .EQ. 0) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Aeresol are Monodisperse.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Aeresol are Monodisperse.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetMonodisperse(.FALSE.)
		ELSE
			CALL SetMonodisperse(.TRUE.)
		END IF

		!!! Condensation?
		Token = GetLine(FH)

		IF (IsInteger(Token) .EQ. 0) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Condensation is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Condensation is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetCondensationFlag(.FALSE.)
		ELSE
			CALL SetCondensationFlag(.TRUE.)
		END IF

		!!! Inorganic Dissolution?
		Token = GetLine(FH)

		IF (IsInteger(Token) .EQ. 0) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Dissolution is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Dissolution is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetDissolutionFlag(.FALSE.)
		ELSE
			CALL SetDissolutionFlag(.TRUE.)
		END IF

		!!! Hydrophobic Organic Dissolution?
		Token = GetLine(FH)

		IF (IsInteger(Token) .EQ. 0) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Hydrophobic Organic Dissolution is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Hydrophobic Organic Dissolution is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetHydrophobicOrgDissolutionFlag(.FALSE.)
		ELSE
			CALL SetHydrophobicOrgDissolutionFlag(.TRUE.)
		END IF

		!!! Hydrophilic Organic Dissolution?
		Token = GetLine(FH)

		IF (IsInteger(Token) .EQ. 0) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Hydrophilic Organic Dissolution is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Hydrophilic Organic Dissolution is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetHydrophilicOrgDissolutionFlag(.FALSE.)
		ELSE
			CALL SetHydrophilicOrgDissolutionFlag(.TRUE.)
		END IF

		!! Is dissolution equilibrium or flux?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether disolution is an equilibrium 0 or flux 1 process. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether disolution is an equilibrium (0) or flux (1) process. Received "//StripToken(Token))		
		CALL SetDissolutionEquilibriumOrFlux (I)
		
		!! Is Ammonia dissolution equilibrium or flux?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether disolution is an equilibrium 0 or flux 1 process. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether ammonia disolution is an equilibrium (0) or flux (1) process. Received "//StripToken(Token))	
		IF (I .EQ. 0) THEN
			CALL SetAmmoniaFlux(.FALSE.)
		ELSE
			CALL SetAmmoniaFlux(.TRUE.)
		END IF
		
		
		!!Metastable Inorganic Aerosol?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether inorganic aerosol is metastable (1) or solid salts are allowed (0). Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether inorganic aerosol is metastable (1) or solid salts are allowed (0). Received "//StripToken(Token))	
		IF (I .EQ. 0) THEN
			CALL SetMetastable(.FALSE.)
		ELSE
			CALL SetMetastable(.TRUE.)
		END IF
		
		!!! Coagulation?
		Token = GetLine(FH)

		IF (IsInteger(Token) .EQ. 0) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Coagulation is considered.  Received "//StripToken(Token))
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether Coagulation is considered.  Received "//StripToken(Token))
		IF (I .EQ. 0) THEN
			CALL SetCoagulationFlag(.FALSE.)
		ELSE
			CALL SetCoagulationFlag(.TRUE.)
		END IF

		!!Constant Coagulation Kernel?
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether coagulationkernel is constant (1) or calculated (0) process. Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether coagulation kernel is constant (1) or calculated (0) process. Received "//StripToken(Token))	
		IF (I .EQ. 0) THEN
			CALL SetConstantKernel(.FALSE.)
		ELSE
			CALL SetConstantKernel(.TRUE.)
		END IF

		!!Ignore corrections to diffusivity? (Also sets curvature corrections to 1.0!)
		Token = GetLine(FH)
		IF (.NOT.IsInteger(Token)) &
			CALL ERROR ("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether diffusivity corrections are ignored (1) or included (0). Received "//StripToken(Token))	
		I = STR2INT(Token)
		IF (I < 0 .OR. I > 1) &
			CALL ERROR("Must Supply either a 0 or 1 ASPInputDeck.in to indicate"//&
					   "whether diffusivity corrections are ignored (1) or included (0).. Received "//StripToken(Token))	
		IF (I .EQ. 0) THEN
			CALL SetIgnoreDiffuseCorrect(.FALSE.)
		ELSE
			CALL SetIgnoreDiffuseCorrect(.TRUE.)
		END IF
			
	END SUBROUTINE ReadMainInputDeck
