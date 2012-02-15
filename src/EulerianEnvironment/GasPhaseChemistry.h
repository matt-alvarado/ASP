!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! GasPhaseChemistry.h
!! This file contains routines for inputing and outputing gas phase
!! chemical concnetrations to the grid points, allocating the 
!! gas phase chemical grid, conveting units, and calculating 
!! gas phase reaction rates.
!!
!! WARNING!: MELAM considers concentrations different from burdens,
!! with Burden = Concentration*GridPointVolume. ASP ignores the 
!! distinction, so burden and concentration are equivalent, and 
!! in units of molecules/cm^{3}. 
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History	
!! 07/30  2007   Matt Alvarado     Removed spurious eps factor from RH calc
!!                                 in SetWaterVaporField						!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. SUBROUTINE ResetGasPhaseReactionRates (ReactionRates,ArraySize)
!!		NOTE: Any changes made to rate constant functions here must be reflected
!!		in subroutine MakeODESet in Chemistry/ChemistryParametersInitRoutines.h
!!
!! 2. SUBROUTINE AllocateGasChemistryGrid (NumbChems)
!! 3. SUBROUTINE SetHomogeneousChemFieldPPB (Index, Value)
!! 4. SUBROUTINE SetHomogeneousChemFieldMGperM3 (Index, Value, MolecMass)
!! 5. SUBROUTINE SetChemFieldFromDeck (Index, FileName)
!! 6. FUNCTION	 GetGasChemVector (XP, YP, ZP) RESULT (Y)
!! 7. SUBROUTINE UpdateChemicalBurdens (XP, YP, ZP, Y)
!! 8. SUBROUTINE UpdateChemicalConcentrations (XP, YP, ZP, Y)
!! 9. FUNCTION   GetGasChems (XP, YP, ZP) RESULT (Y)
!!10. FUNCTION   GetGasChem (XP, YP, ZP, SingleChemIndex) RESULT (Y)
!!11. FUNCTION   GetGasBurden (XP, YP, ZP, SingleChemIndex) RESULT (Y)
!!12. FUNCTION   GetGasChemFromBurden (Burden)
!!13. FUNCTION GetGridCellChemBurden (X, Y, Z, SingleChemIndex)
!!14. SUBROUTINE AddToGridCellChemBurden (X, Y, Z, SingleChemIndex, DeltaChemBurden)
!!15. SUBROUTINE ReplaceGridCellChemBurden (X, Y, Z, SingleChemIndex, ChemBurden)
!!16. FUNCTION   GetPartialPressure (XP, YP, ZP, SingleChemIndex)
!!17. SUBROUTINE SetWaterVaporField (RH)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Each timestep the reaction rates are recalculated !!
	!! This subroutine performs that service.			 !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ResetGasPhaseReactionRates (ReactionRates,ArraySize)

		Use ModelParameters
		USE InfrastructuralCode, ONLY : Transcript, ERROR, Int2Str, Real2Str
		IMPLICIT NONE

		REAL*8,  INTENT(in) :: ReactionRates(:,:)
		INTEGER, INTENT(in) :: ArraySize(2)

		INTEGER :: I,J,K,L
		REAL*8	:: Temp,Mgas,X(3),k0,ki, Y_o, Y_o_term, Y_inf_term, Z, F, &
					alpha, beta, m_o, m_inf, Ratio, Scale, P, k1, k2, k3, kf, Keq, &
					KO, KINF, MO

		!! Loop over all gridcells and reactions
		DO I = 1, XGridPoints(GasPhaseChemGridIndex)
			DO J = 1, YGridPoints(GasPhaseChemGridIndex)
				DO K = 1, ZGridPoints(GasPhaseChemGridIndex)
					DO L = 1, ArraySize(1)
						
						!! Use Select to Differentiate Reaction Types
						SELECT CASE (INT(ReactionRates(L,1)))	! Number of Bodies in the reaction

						!! One body reactions (units are per-second)
						CASE(1)
							!PAN Degradation case (NOT USED! Use type 4!)
							IF (INT(ReactionRates(L,3)) .EQ. 1) THEN
								X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
								TEMP = GetTemp(X(1),X(2),X(3))
								Mgas = GetM   (X(1),X(2),X(3))
								k0   = ReactionRates(L,4)*EXP(ReactionRates(L,5)/Temp)* Mgas
								ki	 = ReactionRates(L,6)*EXP(ReactionRates(L,7)/Temp)
								GasPhaseChemicalRates(I,J,K,L) =							&
									(k0*ki)/(k0+ki) * ReactionRates(L,8)**				&
									(1./(1.+(DLOG10(k0/ki)**2.)))
								
							!Isomerization Case
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 2) THEN
								GasPhaseChemicalRates(I,J,K,L) = ReactionRates(L,4)          *		&
														 Temp**ReactionRates(L,5)    *		&
														 EXP(ReactionRates(L,6)/Temp)**	&
														 ReactionRates(L,7)
							
							!Photolysis Case
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 0) THEN
								GasPhaseChemicalRates(I,J,K,L) = ReactionRates(L,3)
						
							!Heterogeneous Case
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 3) THEN
								!CALL TRANSCRIPT("Hetergeneous Reaction detected!")
								!CALL TRANSCRIPT("Uptake coefficient = "//TRIM(REAL2STR(ReactionRates(L,4)))//"")
								!CALL TRANSCRIPT("Molecular Weight = "//TRIM(REAL2STR(ReactionRates(L,5)))//" g/mol")
								!CALL TRANSCRIPT("Gas Diffusivity = "//TRIM(REAL2STR(ReactionRates(L,6)))//" m2/s")
								!Temporary - call ResetHeteroRates to set rate!
								GasPhaseChemicalRates(I,J,K,L) = 0.0 

							!Equilibrium Degradation case
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 4) THEN
								X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
								TEMP = GetTemp(X(1),X(2),X(3))
								Mgas = GetM   (X(1),X(2),X(3))
								
								k0   = ReactionRates(L,4)*(Temp/300.)**(-1.*ReactionRates(L,5)) * Mgas
								ki	 = ReactionRates(L,6)*(Temp/300.)**(-1.*ReactionRates(L,7))
								kf   =	(k0*ki)/(k0+ki) * ReactionRates(L,8)**				&
									(1./(1.+(DLOG10(k0/ki)**2.)))
								Keq = ReactionRates(L,9)*EXP(ReactionRates(L,10)/TEMP)
								
								GasPhaseChemicalRates(I,J,K,L) = kf/Keq
						
							!Fixed Rate Case
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 5) THEN
								GasPhaseChemicalRates(I,J,K,L) = ReactionRates(L,4)
							
						
							!HNO4 Thermal degradation
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 6) THEN
								 X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
								TEMP = GetTemp(X(1),X(2),X(3))
								MO = GetM   (X(1),X(2),X(3))

								 KO = 4.1E-5*EXP(-10649.2/TEMP)
								 KINF = 5.7E+15*EXP(-11172.6/TEMP)
								 GasPhaseChemicalRates(I,J,K,L) = KO*MO/(1.0+(KO*MO/KINF)) &
													             *0.5**(1.0+(LOG10(KO*MO/KINF))**2)**(-1)
							
							!N2O5 Thermal degradation
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 7) THEN
								 X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
								TEMP = GetTemp(X(1),X(2),X(3))
								MO = GetM   (X(1),X(2),X(3))
								 
								 KO = 1.0E-3*(300.0/TEMP)**3.5*EXP(-11001.5/TEMP)
								 KINF = 9.7E+14*(TEMP/300.0)**0.1*EXP(-11082.0/TEMP)
								 GasPhaseChemicalRates(I,J,K,L) = KO*MO/(1.0+(KO*MO/KINF)) &
													             *0.45**(1.0+(LOG10(KO*MO/KINF))**2)**(-1)
							
							!O1D Thermal deactivation
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 8) THEN
								 X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
								TEMP = GetTemp(X(1),X(2),X(3))
								 
								 GasPhaseChemicalRates(I,J,K,L) = (1.53e11/Temp)*exp(95.6/Temp)
														
							!O + O2 + M => O3 + M 
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 9) THEN
								 X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
								TEMP = GetTemp(X(1),X(2),X(3))
								 
								 GasPhaseChemicalRates(I,J,K,L) = (5.53e16/(Temp**4.8))
							
							!RAD* + O2  
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 10) THEN
								 X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
								TEMP = GetTemp(X(1),X(2),X(3))
								 
								 GasPhaseChemicalRates(I,J,K,L) = 0.9*4.62e7/Temp/60
							
							ELSE
								CALL ERROR("Error in ResetGasPhaseReactionRates - your 1st order rate constant doesn't make sense")								
							END IF
						
						
						!! Two body reactions (units are 1 / ppb / s)
						CASE(2)
						
							X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
							TEMP = GetTemp(X(1),X(2),X(3))
							Mgas = GetM   (X(1),X(2),X(3))
						
							!Normal 2-Body Case
							IF (INT(ReactionRates(L,3)) .EQ. 0) THEN

								GasPhaseChemicalRates(I,J,K,L) = ReactionRates(L,3)          *		&
														 Temp**ReactionRates(L,4)    *		&
														 EXP(ReactionRates(L,5)/Temp)**	&
														 ReactionRates(L,6)
							!Organic Nitrate Case
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 1 .OR. INT(ReactionRates(L,3)) .EQ. 2) THEN
								!!Formula taken from Carter and Atkinson (1989)
								
								alpha = 1.94E-22
								beta = 0.97
								Y_o = alpha*EXP(beta*ReactionRates(L,8))
								m_o = 0.
								m_inf = 8.1

								Y_o_term = Y_o*Mgas*(Temp/300.0)**(-m_o)
								Y_inf_term = 0.826*(Temp/300)**(-m_inf)

								Z = (1+(log10(Y_o_term/Y_inf_term))**(2.0))**(-1.0)
								F = 0.411
								Ratio = ReactionRates(L,9)*(Y_o_term/(1+(Y_o_term/Y_inf_term)))*F**(Z)

								!Special scaling for a few reactions
								Ratio = ReactionRates(L,10)*Ratio + ReactionRates(L,11)
								
								IF(INT(ReactionRates(L,3)) .EQ. 1) THEN
									Scale = Ratio/(1.0 + Ratio)
								ELSE
									Scale = 1.0 - Ratio/(1.0 + Ratio)
								END IF
								
								!IF (ReactionRates(L,10) .EQ. 0.5 .AND. ReactionRates(L,8) .EQ. 10.0) THEN
								!	WRITE(*,*) Scale
									
								!	STOP
								!END IF

								!CALL TRANSCRIPT("Scale for following rxn: "//TRIM(REAL2STR(Scale)))
								GasPhaseChemicalRates(I,J,K,L) = Scale*ReactionRates(L,4)          *		&
														 Temp**ReactionRates(L,5)    *		&
														 EXP(ReactionRates(L,6)/Temp)**	&
														 ReactionRates(L,7)
							
							!This is a CO + OH reaction
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 3) THEN	
								k1 = 1.3e-13
								k2 = 3.19e-33
								GasPhaseChemicalRates(I,J,K,L) = k1 + k2*Mgas
							
							!This is a k = k1 + k2*M case 
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 4) THEN
								k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
								k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
								GasPhaseChemicalRates(I,J,K,L) = k1 + k2*Mgas
							
							!This is a k = k1 + k3*M*(1+k3*M/k2) case 
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 5) THEN
								k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
								k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
								k3 = ReactionRates(L,8)*exp(ReactionRates(L,9)/TEMP)
								GasPhaseChemicalRates(I,J,K,L) = k1 + k3*Mgas*(1 + k3*Mgas/k2)
							
							!This is a heterogeneous reaction 
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 6) THEN
								!Temporary - call ResetHeteroRates to set rate!
								GasPhaseChemicalRates(I,J,K,L) = 0.0
								WRITE(*,*) "2nd order heter rxn detected!"
								 
				
							END IF

						!! Three body reactions, using the troe pressure-broadening 
						!! From from p.1303 of Seinfeld & Pandis
						!! This yields an effective second-order rate constant.
						!! Rates are set in 1 / ppb / ppb / sec 
						CASE(3)
						
							X    = LocateGridPoint(I,J,K,GasPhaseChemGridIndex)	
							TEMP = GetTemp(X(1),X(2),X(3))
							Mgas = GetM   (X(1),X(2),X(3))
						
							!Case for third-order rate constant
							!(like for O + O2 + M => O3 + M) (NOT USED - see first order class 9)
							IF (INT(ReactionRates(L,3)) .EQ. 1) THEN
								GasPhaseChemicalRates(I,J,K,L) = ReactionRates(L,4)*(Temp/300.)**(-1.*ReactionRates(L,5)) * Mgas
						
							!Hardwired OH + NO2 + M Reaction
							ELSE IF (INT(ReactionRates(L,3)) .EQ. 2) THEN
								KO = 2.85E-30*(300.0/TEMP)**2.67
								KINF = 3.13E-11
								MO = Mgas
								GasPhaseChemicalRates(I,J,K,L) = 0.937*KO*MO/(1.0+(0.937*KO*MO/KINF)) &
     							*(EXP(-1.*TEMP/363.0))**(1.0+ &
     							((LOG10(0.937*KO*MO/KINF)-0.12)/(0.85-0.67*LOG10( &
     							EXP(-1.*TEMP/363.0))))**2)**(-1)

							
							!Calculate pseudo-2nd order rate constant
							ELSE
							
								k0   = ReactionRates(L,3)*(Temp/300.)**(-1.*ReactionRates(L,4)) * Mgas
								ki	 = ReactionRates(L,5)*(Temp/300.)**(-1.*ReactionRates(L,6))
								GasPhaseChemicalRates(I,J,K,L) =							&
									(k0*ki)/(k0+ki) * ReactionRates(L,7)**				&
									(1./(1.+(DLOG10(k0/ki)**2.)))
							END IF
						END SELECT

						!CALL TRANSCRIPT("Reaction #"//TRIM(INT2STR(L))//" Rate Constant:"//TRIM(REAL2STR(GasPhaseChemicalRates(I,J,K,L)))//"")
					END DO
				END DO
			END DO
		END DO

	END SUBROUTINE ResetGasPhaseReactionRates


	!!! --- Allocate the gas phase chemical grids --- !!!
	!!! --- Called from  Chemistry Parameters.f90 --- !!!
	SUBROUTINE AllocateGasChemistryGrid (NumbChems)

		USE InfrastructuralCode, ONLY : ERROR
		IMPLICIT NONE

		INTEGER :: NumbChems, Allocation_Error, I, J, K
		
		!! Bring the monster to life
		ALLOCATE (GridGasChem(XGridPoints(GasPhaseChemGridIndex),YGridPoints(GasPhaseChemGridIndex),	&
				  ZGridPoints(GasPhaseChemGridIndex), NumbChems),stat = allocation_error)

				! Allocation Error !
				IF (allocation_error > 0) &
				CALL ERROR ("Allocation of GridGasChem could not proceed in AllocateGasChemistryGrid()")

		!! Water Mixing Ratios have been temporarily stored and are now put
		!! in their final holding place
		IF(.NOT.ALLOCATED(TemporaryWaterArray)) &
		CALL ERROR ("Somehow, you called for an initialization of gas phase chemistry before the relative humidity "// &
					"grid was set, which should have happened while the main input deck was being read in.")
		DO I = 1, XGridPoints(GasPhaseChemGridIndex)
			DO J = 1, YGridPoints(GasPhaseChemGridIndex)
				DO K = 1, ZGridPoints(GasPhaseChemGridIndex)
					GridGasChem(I,J,K,1) = TemporaryWaterArray(I,J,K)
				END DO
			END DO
		END DO
		DEALLOCATE(TemporaryWaterArray)

	RETURN
	END SUBROUTINE AllocateGasChemistryGrid

	SUBROUTINE AllocateGasReactionGrid (NumbRates)

		USE InfrastructuralCode, ONLY : ERROR

		IMPLICIT NONE

		INTEGER :: NumbRates, Allocation_Error

		!! Bring the monster to life
		ALLOCATE (GasPhaseChemicalRates(XGridPoints(GasPhaseChemGridIndex),YGridPoints(GasPhaseChemGridIndex),	&
				  ZGridPoints(GasPhaseChemGridIndex), NumbRates),stat = allocation_error)

				! Allocation Error !
				IF (allocation_error > 0) CALL ERROR ("Allocation of GridGasChem could not proceed in AllocateGasReactionGrid()")

	RETURN
	END SUBROUTINE AllocateGasReactionGrid
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Initialize all of the grid points with the same initial value !!
	!! input value of the chemical concentration is in ppbv			 !!
	!! Stored value is in total abundance.                           !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetHomogeneousChemFieldPPB (Index, Value)

		USE InfrastructuralCode, ONLY : ERROR, Transcript
		USE ModelParameters,	 ONLY : ppbscale,chemscale
		implicit none

		integer :: X, Y, Z, Index
		real*8  :: Value, XX(3), Mgas

		logical :: Scaffolding = .FALSE. !.TRUE.

		!! SCSCSCSC
		IF (SCAFFOLDING) CALL TRANSCRIPT(">>Entering SetHomogeneousChemField<<")

		IF (XGridPoints(GasPhaseChemGridIndex) == 0 .or.	&
			YGridPoints(GasPhaseChemGridIndex) == 0 .or.	&
			ZGridPoints(GasPhaseChemGridIndex) == 0)	    &
			CALL ERROR ("You must call InitializeGridPointFields() before SetGasPhaseChemistry ()"//   &
						"or the program will crash. ")

		!! Loop over the entire domain and set each point to the specified value.
		DO X = 1, XGridPoints(GasPhaseChemGridIndex)
			DO Y = 1, YGridPoints(GasPhaseChemGridIndex)
				DO Z = 1, ZGridPoints(GasPhaseChemGridIndex)
					XX   = LocateGridPoint(X,Y,Z,GasPhaseChemGridIndex)	
					Mgas = GetM (XX(1),XX(2),XX(3))
					GridGasChem(X, Y, Z, Index) = Value*ppbscale*Mgas/ChemScale
				END DO
			END DO
		END DO

		!! SCSCSCSC
		IF (SCAFFOLDING) CALL TRANSCRIPT(">>Exiting SetHomogeneousChemField<<")

	RETURN
	END SUBROUTINE SetHomogeneousChemFieldPPB

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Initialize all of the grid points with the same initial value !!
	!! input value of the chemical concentration is in ppbv			 !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetHomogeneousChemFieldMGperM3 (Index, Value, MolecMass)

		USE InfrastructuralCode, ONLY : ERROR, Transcript
		USE ModelParameters,	 ONLY : chemscale,grams,avogadro
		implicit none

		integer :: X, Y, Z, Index
		real*8  :: Value, XX(3), Mgas
		real*8  :: MolecMass			!! Molecular mass of gas chem # Index -- circular ref to access chemistry module here.

		logical :: Scaffolding = .FALSE. !.TRUE.

		!! SCSCSCSC
		IF (SCAFFOLDING) CALL TRANSCRIPT(">>Entering SetHomogeneousChemField<<")

		IF (XGridPoints(GasPhaseChemGridIndex) == 0 .or.	&
			YGridPoints(GasPhaseChemGridIndex) == 0 .or.	&
			ZGridPoints(GasPhaseChemGridIndex) == 0)	    &
			CALL ERROR ("You must call InitializeGridPointFields() before SetGasPhaseChemistry ()"//   &
						"or the program will crash. ")

		!! Loop over the entire domain and set each point to the specified value.
		DO X = 1, XGridPoints(GasPhaseChemGridIndex)
			DO Y = 1, YGridPoints(GasPhaseChemGridIndex)
				DO Z = 1, ZGridPoints(GasPhaseChemGridIndex)
					GridGasChem(X, Y, Z, Index) = Value/1.e12/	&
												  MolecMass*grams/chemscale*Avogadro
				END DO
			END DO
		END DO

		!! SCSCSCSC
		IF (SCAFFOLDING) CALL TRANSCRIPT(">>Exiting SetHomogeneousChemField<<")

	RETURN
	END SUBROUTINE SetHomogeneousChemFieldMGperM3


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Initialize a Gas Phase Chemical from an Input Deck !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetChemFieldFromDeck (Index, FileName)
		
		USE ModelParameters
		USE InfrastructuralCode

		implicit none

		!! In variables
		integer				 :: Index
		character*(*)		 :: FileName

		!! Local variables
		integer				 :: X, Y, Z, FH
		character (len=50)	 :: Token
		character (len=1024) :: CurrLine
		real*8				 :: XX(3), Mgas

		!! -- This is one potential error -- !!
		IF (XGridPoints(GasPhaseChemGridIndex) == 0 .or.	&
			YGridPoints(GasPhaseChemGridIndex) == 0 .or.	&
			ZGridPoints(GasPhaseChemGridIndex) == 0)		&
				CALL ERROR ("You must call InitializeGridPointFields() before SetGasPhaseChemistry ()"//  &
				"or the program will crash.")

		!! Prepare File and Open It
		FH = GetFileHandle()
	    OPEN(UNIT=FH, FILE=TRIM(InputDeckSubDir)//TRIM(FileName), STATUS='OLD')

		!! How do I error protect on this?
		DO Z = 1, ZGridPoints(GasPhaseChemGridIndex)
			DO Y = 1, YGridPoints(GasPhaseChemGridIndex)

				!! Get the next data line from the input file (this function 
				!! is defined in StringCommands.f90)
				CurrLine = GetLine(FH)

				!! End of Files cause stoppage
				IF (CurrLine .eq. EOF) &
						CALL ERROR ("End Of File of Input Deck "//TRIM(StripToken(FileName))//"reached when seeking "// &
						"range of X values at Grid Points (Y,Z)= ("//TRIM(INT2STR(Y))//","//			 &
						TRIM(INT2STR(Z))//"). Expected dimensions for this input deck are (#X,#Y,#Z) = ("&
						//TRIM(INT2STR(XGridPoints(GasPhaseChemGridIndex)))//","//						 &
						TRIM(INT2STR(YGridPoints(GasPhaseChemGridIndex)))//","//  &
						TRIM(INT2STR(ZGridPoints(GasPhaseChemGridIndex)))//")")

				!! Each line contains all of the X input values for the particular Y and Z grid point
				DO X = 1, XGridPoints(GasPhaseChemGridIndex)

					CALL GetToken(CurrLine, ";", Token)

					!! If we want more tokens than are given, put an error
					IF (LEN_TRIM(StripToken(Token)) == 0) &
						CALL ERROR("End of line reached in Input Deck "//TRIM(StripToken(FileName))//"when seeking "&
							       //"the value at X grid point #"//TRIM(INT2STR(X))//" at Grid Points (Y,Z)= ("//  &
								   TRIM(INT2STR(Y))//","//TRIM(INT2STR(Z))//  &
								   ").  Expected dimensions for this input deck are (#X,#Y,#Z) = ("//&
								   TRIM(INT2STR(XGridPoints(GasPhaseChemGridIndex)))//","//  &
								   TRIM(INT2STR(YGridPoints(GasPhaseChemGridIndex)))//&
								   ","//TRIM(INT2STR(ZGridPoints(GasPhaseChemGridIndex)))//")")

					XX   = LocateGridPoint(X,Y,Z,GasPhaseChemGridIndex)	
					Mgas = GetM   (XX(1),XX(2),XX(3))

					!! If there are no errors, set the value
					GridGasChem(X,Y,Z,index) = STR2REAL(Token)*ppbscale*Mgas

				END DO
			END DO
		END DO

		CLOSE(FH)
		CALL ReturnFileHandle(FH)

	RETURN
	END SUBROUTINE SetChemFieldFromDeck

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Give a vector of the current chemical concentrations to the outside world. !!
	!!  XP,YP,and ZP in this case are grid points.								  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION GetGasChemVector (XP, YP, ZP) RESULT (Y)

		IMPLICIT NONE

		!! External Variables
		INTEGER	:: XP, YP, ZP
		REAL*8  :: Y(ylen)

		!! Internal Variables
		INTEGER   :: I

		DO I = 1, HowManyGasChems
			Y(I) = GridGasChem(XP, YP, ZP, I)
		END DO

		IF (ylen > HowManyGasChems) THEN
			DO I = HowManyGasChems+1, ylen
				Y(I) = 0.
			END DO
		END IF
	END FUNCTION GetGasChemVector

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Once a routine calculates new values for the chemicals at a grid point !!
	!! (which it does locally), it places those values back in the grid.      !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Call it like this (assuming the vector you're using is longer than the !!
	!! evolving gas chemical length:										  !!
	!!    CALL UpdateChemicalConcentrations(XP,YP,ZP,Y(1:HowManyEvolveGasChems))
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE UpdateChemicalBurdens (XP, YP, ZP, Y)

		IMPLICIT NONE

		!! External Variables
		INTEGER :: XP, YP, ZP
		REAL*8  :: Y(HowManyEvolveGasChems)

		!! Internal Variables
		INTEGER :: I

		!! We only replace those chemicals that are allowed to evolve 
		!! (which are listed first, by construction)
		DO I = 1, HowManyEvolveGasChems
			IF (Y(I) .EQ. -1) CYCLE
			CALL ReplaceGridCellChemBurden (XP, YP, ZP, I, Y(I))
		END DO

	END SUBROUTINE

	SUBROUTINE UpdateChemicalConcentrations (XP, YP, ZP, Y)

		IMPLICIT NONE

		!! External Variables
		INTEGER :: XP, YP, ZP
		REAL*8  :: Y(HowManyEvolveGasChems)

		!! Internal Variables
		INTEGER :: I

		!! We only replace those chemicals that are allowed to evolve 
		!! (which are listed first, by construction)
		DO I = 1, HowManyEvolveGasChems
			GridGasChem(XP, YP, ZP, I) = Y(I)
		END DO

	END SUBROUTINE

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Chemistry Vector Y() Interpolated to Any Point !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION GetGasChems (XP, YP, ZP) RESULT (Y)

		IMPLICIT NONE

		REAL*8				:: Y(HowManyGasChems)	! Output
		REAL*8, INTENT(in)	:: XP, YP, ZP			! Input
		INTEGER				:: i

		!! Interpolate each chemical using the interpolation function in InterpolateGrid.f90
		DO i = 1,HowManyGasChems
			Y(i) = InterpolateGrid (XP, YP, ZP, GridGasChem(:,:,:,i), GasPhaseChemGridIndex)
		END DO

		RETURN
	END FUNCTION GetGasChems

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve a single Gas Phase Chemical Interpolated to Any Point. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetGasChem (XP, YP, ZP, SingleChemIndex) RESULT (Y)

		implicit none

		REAL*8, intent(in)	:: XP, YP, ZP			! Input
		INTEGER, intent(in) :: SingleChemIndex

		!! Interpolate the chemical using the interpolation function in InterpolateGrid.f90
		Y = InterpolateGrid (XP, YP, ZP, GridGasChem(:,:,:,SingleChemIndex), GasPhaseChemGridIndex)

		RETURN
	END FUNCTION GetGasChem

	REAL*8 FUNCTION GetGasBurden (XP, YP, ZP, SingleChemIndex) RESULT (Y)

		implicit none

		REAL*8, intent(in)	:: XP, YP, ZP			! Input
		INTEGER, intent(in) :: SingleChemIndex

		!! Interpolate the chemical using the interpolation function in InterpolateGrid.f90
		Y = InterpolateGrid (XP, YP, ZP, GridGasChem(:,:,:,SingleChemIndex), GasPhaseChemGridIndex)

		RETURN
	END FUNCTION GetGasBurden

	REAL*8 FUNCTION GetGasChemFromBurden (Burden)

		IMPLICIT NONE

		REAL*8, intent (in) :: Burden

		GetGasChemFromBurden = Burden

		RETURN
	END FUNCTION GetGasChemFromBurden

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve a single Gas Phase Chemical Total Burden within a Grid Cell !!
	!! (this is just the concentration - I've removed the scaling by volume (Matt Alvarado, 4/12/06).									!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Input location is A GRID CELL INDEX rather than a location.			!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetGridCellChemBurden (X, Y, Z, SingleChemIndex)

		USE InfrastructuralCode, ONLY : ERROR

		IMPLICIT NONE

		INTEGER, intent(in) :: SingleChemIndex, X, Y, Z
		
		IF (X .GT. XGridPoints(GasPhaseChemGridIndex) .OR. &
            Y .GT. YGridPoints(GasPhaseChemGridIndex) .OR. &
		    Z .GT. ZGridPoints(GasPhaseChemGridIndex)) &
		CALL ERROR ("In GetGridCellChemBurden(), called for a grid cell that is beyond the domain of the model.")

		GetGridCellChemBurden = GridGasChem(X,Y,Z,SingleChemIndex)

		RETURN
	END FUNCTION GetGridCellChemBurden

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This one replaces that same burden number when condensation is over !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE AddToGridCellChemBurden (X, Y, Z, SingleChemIndex, DeltaChemBurden)

		implicit none

		REAL*8,  intent(in) :: DeltaChemBurden
		INTEGER, intent(in) :: SingleChemIndex, X, Y, Z

		GridGasChem(X,Y,Z,SingleChemIndex) = GridGasChem(X,Y,Z,SingleChemIndex) + &
											 DeltaChemBurden

		RETURN
	END SUBROUTINE AddToGridCellChemBurden

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This one replaces that same burden number when condensation is over !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ReplaceGridCellChemBurden (X, Y, Z, SingleChemIndex, ChemBurden)

		implicit none

		REAL*8,  intent(in) :: ChemBurden
		INTEGER, intent(in) :: SingleChemIndex, X, Y, Z

		GridGasChem(X,Y,Z,SingleChemIndex) = ChemBurden

		RETURN
	END SUBROUTINE ReplaceGridCellChemBurden


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve a single Gas Phase Chemical Partial Pressure Interpolated to Any Point. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetPartialPressure (XP, YP, ZP, SingleChemIndex)

		USE ModelParameters, ONLY : kB

		implicit none

		REAL*8, intent(in)	:: XP, YP, ZP			! Input
		INTEGER, intent(in) :: SingleChemIndex

		!! Interpolate the chemical using the interpolation function in InterpolateGrid.f90
		!! And calculate partial pressure using Dalton's Law
		GetPartialPressure = kB * GetTemp(XP,YP,ZP) * GetM(XP,YP,ZP) * &
							 InterpolateGrid (XP, YP, ZP, GridGasChem(:,:,:,SingleChemIndex), GasPhaseChemGridIndex)

		RETURN
	END FUNCTION GetPartialPressure

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Populate the Gas Phase Water Vapor Fields !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetWaterVaporField (RH)

		USE InfrastructuralCode, ONLY : ERROR
		USE ModelParameters,     ONLY : ppbscale, chemscale

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: RH

		!! Internal Variables
		INTEGER :: I,J,K,allocation_error
		REAL*8  :: XX(3),SPress !, MM

		!! Bring the monster to life
		ALLOCATE (TemporaryWaterArray(XGridPoints(GasPhaseChemGridIndex),YGridPoints(GasPhaseChemGridIndex),	&
				  ZGridPoints(GasPhaseChemGridIndex)),stat = allocation_error)
		IF (allocation_error > 0) CALL ERROR ("Allocation of TemporaryWaterArray could not proceed in SetWaterVaporField()")

		DO I = 1, XGridPoints(GasPhaseChemGridIndex)
			DO J = 1, YGridPoints(GasPhaseChemGridIndex)
				DO K = 1, ZGridPoints(GasPhaseChemGridIndex)
					XX = LocateGridPoint(I, J, K, GasPhaseChemGridIndex)
					SPress = GetSatVapPress(XX(1),XX(2),XX(3))
					TemporaryWaterArray(I,J,K)  = ChemScale * SPress * RH * GetM(XX(1),XX(2),XX(3)) / &
												  (GetPress(XX(1),XX(2),XX(3)))  
				END DO
			END DO
		END DO

	END SUBROUTINE SetWaterVaporField
