!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ASP Model								        !!
!!  Developed by Matt Alvarado (mjalvara@mit.edu) for his thesis		!!
!!  with Prof. R.G. Prinn (rprinn@mit.edu) at MIT. (2004-2008)			!!
!!										!!
!!  ASP is an updated, expanded version of the					!!
!!  Mixed Eulerian-Lagrangian Aerosol Model (MELAM) developed by 		!!
!!	Donnan Steele (donnan@mit.edu) for his thesis project			!!
!!  with Prof. R.G.Prinn at MIT. (2000-2004)					!!
!!										!!
!!	+ Always compile promoting all reals to *8				!!
!!										!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY								!!
!!										!!
!! Month  Year   Name              Description					!!
!! 07     2006   Matt Alvarado     Began Update History				!!
!! 08/25  2006   Matt Alvarado     Changed GOTO loops to DO WHILE loops         !!
!! 09/06  2006   Matt Alvarado     Added particle deposition to size dist. test !!
!! 09/07  2006   Matt Alvarado     Added gas deposition to size dist. test	!!
!! 10/03  2007   Matt Alvarado     Minor Changes                                !!
!! 08/27  2010   Matt Alvarado     Minor Changes                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM ASP

	!! Include modules
	!USE IFPORT,				 ONLY : RTC, ETIME
	USE ModelParameters
	USE InfrastructuralCode	
	USE GridPointFields
	USE Chemistry
	USE Aerosols		
	USE Condensation
	USE OutputRoutines
	USE Coagulation
	
	implicit none

	!! Define the local variables
	integer :: i, j, k, l, seedarray(2),q, r, NumBins, A, S, HeteroBins, &
	           Gas, Aq, HNO3GasIndex, NumAqChems, NH3GasIndex, PotassiumIndex, &
			   HClGasIndex, SulfateIndex, BisulfateIndex, K2SO4Index, SOA2GasIndex, NumOrgChems
	Integer :: 	KNO3i, KCli, KHSO4i, UR21i, &
				Nai, NaNO3i, NaCli, NaHSO4i, Na2SO4i, &
				NH3i, NH4i, NH4NO3i, NH4Cli, NH4HSO4i, NH42SO4i, &
				Cli, NO3i, SO4i, HSO4i, K2SO4i,LevAeroIndex, AmSulfAeroIndex, H2SO4AeroIndex, N, HNO3Index, &
				NH3Gas, HNO3Gas, HClGas, H2SO4Gas, NumAqOrgChems, ReturnType
	real*8  :: ii, jj, kk, y(200), Dum1
	real*8  :: T_mix, Tau_Dep, slope, b, ChemTimestep, MixTimestep, &
				CoagTimeStep, CondTimeStep, t, Mgas, Temp, &
				X(3), NewRH, HNO3Conc, NH3Conc, ut, &
				Potassium, HClConc, TotalNumber, TotalVolume, TotalSulfateMass, TotalPotassiumMass, SOA2Conc, &
				TotalPOA, TotalSOA, TotalK, TotalNa, TotalNH4, TotalCl, TotalNO3, TotalSO4, &
				NucRate, MolsinNucParticle, GridTotalSOA, RH, tiny, x1, x2, x3, f1, f2, f3, chg, xs, molalbin, &
				K2SO4Mass, Kmass, TotalMass, TimeStepSize, PressEq, P1, EnvMass, InitialMass, NewMass

	REAL*4  :: etimearray(2)
	REAL*8, ALLOCATABLE	::	XPASS(:), GAMMA(:), Chem(:), NumConc(:), &
				   HeteroNumberConc(:), HeteroRadii(:), &
							InitNumConc(:), EnvNumConc(:)
	REAL*8, DIMENSION(18) :: ExtCoeff, SingScat, Assym
	REAL*8	:: Term(250), Total, Yield
	INTEGER :: GP(3), clock, time(8)
        INTEGER, ALLOCATABLE :: outFH(:), OFH(:)

	TYPE (DiffEqTerm), POINTER :: CurrentTerm
	TYPE (ChemTerm),   POINTER :: CurrentChem

							
	
	LOGICAL :: DoTestCoagulation, DoHumidityTest
	
	TYPE(Particle), POINTER :: Current, Env
	

	
	IF(.FALSE.) CALL EulerianBox()

	!! Initialize the random number generator 
        call DATE_AND_TIME(values=time)     ! Get the current time (for pgi compiler)
        seedarray(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
        !seedarray(1) = RTC() !for ifort compiler
	seedarray(2) = FLOOR(SQRT(REAL(seedarray(1))))
	CALL Random_Seed(put=seedarray)

	!! Send a Welcome Message
	CALL Transcript ("**********************************************************")
	CALL Transcript ("** Welcome to the ASP Model developed by Matt Alvarado  **")
	CALL Transcript ("** (mjalvara@mit.edu).                                  **")
	CALL Transcript ("** ASP is an updated, expanded version of MELAM         **")
	CALL Transcript ("** by H. D. Steele (donnan@mit.edu).                    **")
	CALL Transcript ("**********************************************************")
	
	IF(.FALSE.) THEN !Gas Chemistry Only
		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		!Small domain compensates for large particle #
		CALL SetDomainSize (1*cm, 1*cm, 1*cm)
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()
	
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL ReadPhotolysisRates
	
		!Calculate reaction rates for non-photolysis and non-hetero gas rxns
		CALL RecalculateReactionRates(IFGAS,IFAQ)
		
		!Calculate initial hetero. Rxn rates
		HeteroBins = 1
		ALLOCATE(HeteroNumberConc(HeteroBins))
		ALLOCATE(HeteroRadii(HeteroBins))

		HeteroNumberConc(1) = 1.5e5 !Guess for Williams fire
                !1.5e5 !Guess for Williams fire
                !3.0e4 !Guess for Alaska
		!3.1e4 !Otavi (Haywood et al, 2003)  
		!1.433e5 !Timbavati (Trentmann et al, 2005)
		HeteroRadii(1) = 1.0e-7 ! in m
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		
		!Calculate Initial Photolysis Rates
		ut =  19.0 !Yokelson Williams fire, 11 AM local time
                !19.0 !Yokelson Lake McKay fire, 11 AM local time
                !2.51 !ARCTAS Lake McKay Late Pass time in hours 
                ! 0.87 !Alaska Universal time in hours
		!11.00 !Otavi Initial Universal Time in hours
		!8.75 !Timbavati Initial Universal Time in hours
	
                CALL RecalculatePhotoRates (ut, .TRUE.)
      
		ALLOCATE(CHEM(HowManyGasChems))

		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
		Mgas = GetM   (X(1),X(2),X(3))
         

		ALLOCATE(OutFH(3))
		OutFH(1) = GetFileHandle()
		OutFH(2) = GetFileHandle()
		OutFH(3) = GetFileHandle()

		OPEN(UNIT=OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
		WRITE(OutFH(1), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT=OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
		WRITE(OutFH(2), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT=OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")
		WRITE(OutFH(3), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
	
		DO q = 1, HowManyGasChems
			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE 
				WRITE(OutFH(3),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			END IF
		END DO
	
		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""
	
		t=0
		ChemTimestep = 1.0
		MixTimestep = 1*ChemTimestep !Should always be an integer multiple!
		DO r = 1, INT(390*60/MixTimestep)+1
		
			!Write Output to file in units of ppb
			IF (MOD(t,60.0)<=0.05) THEN
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') t
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') t
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') t
		
				DO q = 1, HowManyGasChems
		
					Chem(q) = GridGasChem(1,1,1,q) !molecules/cm3

					IF (q .LE. 70) THEN
						WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
						WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE 
						WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					END IF
				END DO
		
				WRITE(OutFH(1), FMT='(A1)') ""
				WRITE(OutFH(2), FMT='(A1)') ""
				WRITE(OutFH(3), FMT='(A1)') ""
			END IF !Write Output

			
			!NOTE: Right now, pressure and temperature are assumed to be 
			!constant, so we don't recalculate the rest of the reaction rates
			!in order to save computer time (including this doubles computer time)
			!CALL RecalculateReactionRates(IFGAS,IFAQ)
			!WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
			!after this step, or those rates will be set to 0!
			
			!Recalculate Photolysis Rates
			CALL RecalculatePhotoRates(ut, .FALSE.)

			!Recalculate Heterogeneous Rates
			CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .FALSE.)
			
			!Step Gas Chemistry
			CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep),1,1,1)

			Slope = 2.0 
			b = 212.22 !Yokelson Williams Fire (California, Nov. 17 2009) estimated from yo = 1km, Ky = 1178 m2/s
                        !7884 !ARCTAS Lake McKay Late Pass from Jason St. Clair peak CO data  
                                 !(Ky = 8632 m2/s, yo = 16.5 km)
                        !212.22 !Yokelson Williams Fire (California, Nov. 17 2009) estimated from yo = 1km, Ky = 1178 m2/s
                        !1387.1  !Alaska Base
                        !1731.0 !Arctas Early pass
                        !1387.1 !Alaska Base
			!7000.0 !Alaska Slow
			!350.0 !Alaska Fast
	
			!156.25 !Otavi
			
			!142 + (6.0/7.0) !Timbavati

			!400 !Timbavati slow
			!66 + (2.0/3.0) !Timbavati Fast  
			!1.0e50 !No Dilution 	
			!12000	 
			T_mix = Slope*t+b
		
			!Step Dilution of Chems
			DO q = 1, HowManyGasChems
				Chem(q) = GridGasChem(1,1,1,q)
				Chem(q) = (Chem(q) + (GasPhaseBackground(q)*MixTimestep/T_mix))/(1+ MixTimestep/T_mix)
			END DO		
		
			!Update chemical concentrations
			CALL UpdateChemicalConcentrations(1, 1, 1, Chem(1:HowManyEvolveGasChems))
			
			!Step Dilution of "particles" for heterogeneous chemistry
			DO q = 1, HeteroBins
				HeteroNumberConc(q) = HeteroNumberConc(q)/(1+ MixTimestep/T_mix)
			END DO
			
			t = t + MixTimestep
			ut = ut + (MixTimestep/3600.0)
	        IF(ut .GE. 24.0) ut = 0.0
				
		END DO !Time Loop

	END IF !Gas Chemistry Only

	IF (.TRUE.) THEN !BULK PARTICLE TESTS
		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize (1*cm, 1*cm, 1*cm)
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()
	
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL ReadPhotolysisRates
		
		CALL RecalculateReactionRates(IFGAS,IFAQ)
		
		!Set heterogeneous Reactions
		HeteroBins = 1
		ALLOCATE(HeteroNumberConc(HeteroBins))
		ALLOCATE(HeteroRadii(HeteroBins))

		HeteroNumberConc(1) = 1.5e5 !Guess for Williams fire
                !1.433e5 !Timbavati (Trentmann et al, 2005)
		!3.1e4 !Otavi (Haywood et al, 2003)
		HeteroRadii(1) = 1.0e-7
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		
		!Calculate Initial Photolysis Rates
		ut =  2.5 !ARCTAS Lake McKay Late Pass time in hours 
                !8.75 !Timbavati Initial Universal Time in hours
		!11.00 !Otavi Initial Universal Time in hours
		CALL RecalculatePhotoRates (ut, .TRUE.)

		ALLOCATE(CHEM(HowManyGasChems))

		!! INITIALIZE the DISSOLUTION routine
		CALL SetAllDissolution
		!WRITE(*,*) "Dissolution Okay"
		!STOP

		!! INITIALIZE the PARTICLES
		CALL InitializeParticlesSectional	
		current => particles(1,1,1)%first
		Dum1 = HydrophilicOrganicWaterContent(current, 298.15, 0.5555)
		
		
		WRITE(*,*) "Check 1", current%EffectiveRadius
		!!WARNING Force initial equilibrium between gas and aerosol
		!!(Uses open system for water eq.)
		CALL EquilibrateInternallyAtGridPoint(1,1,1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
                 !CALL EquilibrateGridPointAll(1, 1, 1, .TRUE., .TRUE., .TRUE.)
		WRITE(*,*) "Check 2"
                CALL EquilibrateHydrophobicPhaseAtGridPoint(1,1,1)
		WRITE(*,*) "Check 3"
		
		CALL SpillBeans(1,1,1)

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
		
		!Count # of bins, # of Aq. Chems
		current => particles(1,1,1)%first
		NumBins = 0
		DO WHILE(associated(current))
			NumBins = NumBins + 1
			NumAqChems = size(current%AqChems)
			NumOrgChems = size(current%OrgChems)
                        NumAqOrgChems = size(current%AqOrgChems)
		    current => current%next
		END DO

		ALLOCATE(NumConc(NumBins))
		
		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
		Mgas = GetM   (X(1),X(2),X(3))

		ALLOCATE(OutFH(4))
		OutFH(1) = GetFileHandle()
		OutFH(2) = GetFileHandle()
		OutFH(3) = GetFileHandle()
		OutFH(4) = GetFileHandle()

		OPEN(UNIT = OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
		WRITE(OutFH(1), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
		WRITE(OutFH(2), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")
		WRITE(OutFH(3), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(4), FILE=TRIM(OutputDeckSubDir)//"BulkAerosolOut.txt", STATUS="REPLACE")
		WRITE(OutFH(4), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		
		DO q = 1, HowManyGasChems
			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE 
				WRITE(OutFH(3),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			END IF
		END DO

		IF (NumBins .GT. 1) CALL ERROR("The current output routine needs a bulk aerosol to work.") 
	
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "Number", " "
		WRITE(OutFH(4),FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
		WRITE(OutFH(4),FMT='(A10, A3)', ADVANCE='NO') "Ionic Str.", " "
		WRITE(OutFH(4),FMT='(A8, A5)', ADVANCE='NO') "H2O Act.", " "
		WRITE(OutFH(4),FMT='(A10, A3)', ADVANCE='NO') "Sol. Dens.", " "
		WRITE(OutFH(4),FMT='(A11, A2)', ADVANCE='NO') "Part. Dens.", " "
		WRITE(OutFH(4),FMT='(A2, A11)', ADVANCE='NO') "Ph", " "
		WRITE(OutFH(4),FMT='(A4, A9)', ADVANCE='NO') "Temp", " "
		DO I = 1, NumAqChems
			WRITE(OutFH(4),FMT='(A13)', ADVANCE='NO') AqPhaseChemicalNames(I)
                        WRITE(*,*) I, AqPhaseChemicalNames(I)
		END DO
		DO I = 1, NumOrgChems		
					WRITE(OutFH(4),FMT='(A13)', ADVANCE='NO') OrgPhaseChemicalNames(I)  
                                        WRITE(*,*) I, OrgPhaseChemicalNames(I)
		END DO
		DO I = 1, NumAqOrgChems		
					WRITE(OutFH(4),FMT='(A13)', ADVANCE='NO') AqOrgPhaseChemicalNames(I)
                                        WRITE(*,*) I, AqOrgPhaseChemicalNames(I)
		END DO
		WRITE(OutFH(4),FMT='(A5, A8)', ADVANCE='NO') "TotBC", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotPOA", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotSOA", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotH2O", " "
		WRITE(OutFH(4),FMT='(A4, A9)', ADVANCE='NO') "TotK", " "
		WRITE(OutFH(4),FMT='(A5, A8)', ADVANCE='NO') "TotNa", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotNH4", " "
		WRITE(OutFH(4),FMT='(A5, A8)', ADVANCE='NO') "TotCl", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotNO3", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotSO4", " "
					
		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""
		WRITE(OutFH(4), FMT='(A1)') ""

		t=0
		Slope = 2.0 
		b = 7884 !ARCTAS Lake McKay Late Pass from Jason St. Clair peak CO data  
                                 !(Ky = 8632 m2/s, yo = 16.5 km)
                !142 + (6.0/7.0) !Timbavati
		!400.0 !Timbavati Slow
		!1.0e50 !No Dilution
		!66 + (2.0/3.0) !Timbavati Fast 
		!156.25 !Otavi
		Tau_Dep = 1.0e50 !No Deposition
				!1600.0	
		ChemTimestep = 1.0
		MixTimestep = 1*ChemTimestep !Should always be an integer multiple!
		CondTimeStep = 60*MixTimestep
		
		!Run the main loop with above timesteps
		DO r = 1, INT(119*60/MixTimestep)

			!Write Output to file every minute in units of ppb
			IF (MOD(t,60.0)<=0.05) THEN
				
				!Write current time in minutes
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				
				!Write gas phase concentrations in ppb
				DO q = 1, HowManyGasChems
		
					Chem(q) = GridGasChem(1,1,1,q) !molecules/cm3

					IF (q .LE. 70) THEN
						WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
						WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE 
						WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					END IF
				END DO

				
				IF (NumBins .GT. 1) CALL ERROR("The current output routine needs a bulk aerosol to work.") 
	
				current => particles(1,1,1)%first
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature
				!ug/m3 air
				DO I = 1, NumAqChems
					WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(I)*current%numberofparticles*1.0e12*AqMolecularMass(I)
				END DO
				DO I = 1, NumOrgChems		
					WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(I)*current%numberofparticles*1.0e12*OrgMolecularMass(I)
				END DO
	                        DO I = 1, NumAqOrgChems		
					WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqOrgchems(I)*current%numberofparticles*1.0e12*AqOrgMolecularMass(I)
				END DO
				!BC Mass Conc (ug/m3)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(1)*current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
			!POA Mass Conc (ug/m3)
			TotalPOA = 0.
			DO J = 2, 9
				TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			DO J = 20, NumOrgChems
				TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
			END DO

                        DO J = 20, NumAqOrgChems
                                TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
                        END DO

			WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

			!SOA Mass Conc (ug/m3)
			TotalSOA = 0.
			DO J = 10, 19
				TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalSOA = TotalSOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
				!Water Mass Conc. (ug/m3)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Aqchems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1)
					
				!K+ Mass Conc (ug/m3)
				TotalK = 0.		
				DO J = 1, NumAqChems
					IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
						TotalK = TotalK + current%Aqchems(J)*current%numberofparticles*1.0e12
					ELSE IF (J .EQ. K2SO4Index) THEN
						TotalK = TotalK + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalK = TotalK*AqMolecularMass(PotassiumIndex)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
				!Na+ Mass Conc (ug/m3)
				TotalNa = 0.		
				DO J = 1, NumAqChems
					IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
						TotalNa = TotalNa + current%Aqchems(J)*current%numberofparticles*1.0e12
					ELSE IF (J .EQ. Na2SO4i) THEN
						TotalNa = TotalNa + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalNa = TotalNa*AqMolecularMass(Nai)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
				!NH4+ Mass Conc (ug/m3)
				TotalNH4 = 0.		
				DO J = 1, NumAqChems
					IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
						TotalNH4 = TotalNH4 + current%Aqchems(J)*current%numberofparticles*1.0e12
					ELSE IF (J .EQ. NH42SO4i) THEN
						TotalNH4 = TotalNH4 + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

				!Cl- Mass Conc (ug/m3)
				TotalCl = 0.		
				DO J = 1, NumAqChems
					IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
						TotalCl = TotalCl + current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalCl = TotalCl*AqMolecularMass(Cli)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

				!NO3- Mass Conc (ug/m3)
				TotalNO3 = 0.
				DO J = 1, NumAqChems
					IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
						TotalNO3 = TotalNO3 + current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

				!SO4- Mass Conc (ug/m3)
				TotalSO4 = 0.
				DO J = 1, NumAqChems
					IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
						.OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
						TotalSO4 = TotalSO4 + current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

				!Advance to next line
				WRITE(OutFH(1), FMT='(A1)') ""
				WRITE(OutFH(2), FMT='(A1)') ""
				WRITE(OutFH(3), FMT='(A1)') ""
				WRITE(OutFH(4), FMT='(A1)') ""
			END IF !Write Output
		
			!Calculate Temperature
			X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
			TEMP = GetTemp(X(1),X(2),X(3))

			!NOTE: Right now, pressure and temperature are assumed to be 
			!constant, so we don't recalculate the rest of the reaction rates
			!CALL RecalculateReactionRates(IFGAS,IFAQ)
			!WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
			!after this step, or those rates will be set to 0!
			
			!Recalculate Photolysis Rates
			CALL RecalculatePhotoRates(ut, .FALSE.)

			!Recalculate Heterogeneous Rates
			CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .FALSE.)
						
			!Step Gas Chemistry
			CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep),1,1,1)
			
			!Only Equilibrate every 60 s (one second before the minute)
			IF (MOD(t,CondTimeStep) .LT. 59.01 .AND. MOD(t,CondTimeStep) .GT. 58.99) THEN

				!Equilibrate Gases and particles
				CALL EquilibrateSulfate(1, 1, 1)
				CALL EquilibrateGridPointAll(1,1,1,EquilibrateWater = .TRUE., ForceWaterEquilibrium=.FALSE., WaterOpenSystem=.FALSE.)
				!CALL EquilibrateGridPointAll(1,1,1,EquilibrateWater = .FALSE.)
				CALL EquilibrateHydrophobicPhaseAtGridPoint(1,1,1)
				!WRITE(*,*) "Eq. Okay"

			END IF
			
			!Calculate T_mix
			T_mix = Slope*t+b
		
			!Step Dilution of Chems
			DO q = 1, HowManyGasChems
				Chem(q) = GridGasChem(1,1,1,q)
				Chem(q) = (Chem(q) + (GasPhaseBackground(q)*MixTimestep/T_mix))/(1+ MixTimestep/T_mix)
			END DO	

			!Update chemical concentrations
			CALL UpdateChemicalConcentrations(1, 1, 1, Chem(1:HowManyEvolveGasChems))
			
			!Step Dilution of Particles
			current => particles(1,1,1)%first
			I = 1
			DO WHILE(associated(current))
									
				!Get Number Concentration of particles in particles/cm3
				NumConc(I) = current%numberofparticles; 
	
				!Dilute NumConc (assume no background for now)
				NumConc(I) = NumConc(I)/(1+ MixTimestep/T_mix + MixTimestep/Tau_Dep)

				!Update Particle list
				current%numberofparticles = NumConc(I); 

				!Force particle temperature to be env. Temperature
				!WARNING: This is a kludge, since the program is calculating
				!a huge temperature for the largest particles, and I'm not sure why
				current%Temperature = TEMP
				
				I = I + 1

				
				current => current%next
			END DO !END STEP PARTICLE DILUTION			
			
			!Spill Beans
			!CALL SpillBeans(1,1,1)
									
			!Step Dilution of "particles" for heterogeneous chemistry
			!and recalculate heterogeneous reaction rates
			DO q = 1, HeteroBins
				HeteroNumberConc(q) = HeteroNumberConc(q)/(1+ MixTimestep/T_mix + MixTimestep/Tau_Dep)
			END DO
						
			!Update time
			t = t + MixTimestep

			IF (MOD(t,60.0)<=0.05) THEN
				!WRITE(*,*) trim(real2str(t/60)), " minutes done"
			END IF

		END DO

		!Spill Beans
		CALL SpillBeans(1,1,1)
		    
	END IF !Bulk Particle Tests

	IF (.FALSE.) THEN !Size Dist PARTICLE TESTS
		
		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize (1*cm, 1*cm, 1*cm)
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()
	
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL ReadPhotolysisRates
		
		CALL RecalculateReactionRates(IFGAS,IFAQ)
			
		!Calculate Initial Photolysis Rates
		ut = 8.75 !Timbavati Initial Universal Time in hours
		!11.00 !Otavi Initial Universal Time in hours
		CALL RecalculatePhotoRates (ut, .TRUE.)

		ALLOCATE(CHEM(HowManyGasChems))
	
		!! INITIALIZE the DISSOLUTION routine
		CALL SetAllDissolution

		!! INITIALIZE the PARTICLES	
		CALL InitializeParticlesSectional
		CALL RegridAerosol (1, 1, 1)
		CALL SpillBeans(1,1,1)
		!WRITE(*,*) "AOD: ", AerosolOpticalDepth (1, 1, 1)

		!! Call Boundary distributions
		CALL ReadBoundaryDistributionsOrganic
		CALL PopBoundaryPart
		CALL SortBoundaryPart
		CALL SpillBeans(1,1,1,Env=.TRUE.)
		
		!!WARNING Force initial water equilibrium between gas and aerosol
		!!(Uses open system for water eq.)
		CALL EquilibrateInternallyAtGridPoint(1,1,1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
		!CALL EquilibrateGridPointAll(1, 1, 1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
		CALL RegridAerosol (1, 1, 1)
		CALL SpillBeans(1,1,1)
	
	
		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
		Mgas = GetM   (X(1),X(2),X(3))
		
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

		!Count # of bins, # of Aq. Chems
		current => particles(1,1,1)%first
		Env => BoundaryParticles%first
		NumBins = 0
		Potassium = 0
		DO WHILE(associated(current))
			NumBins = NumBins + 1
			NumAqChems = size(current%AqChems)
			NumOrgChems = size(current%OrgChems)
			NumAqOrgChems = size(current%AqOrgChems)
			Potassium = Potassium + current%AqChems(PotassiumIndex)*current%numberofparticles
			Potassium = Potassium + current%AqChems(KNO3i)*current%numberofparticles
			Potassium = Potassium + 2*current%AqChems(K2SO4i)*current%numberofparticles
			Potassium = Potassium + current%AqChems(KCli)*current%numberofparticles
			Potassium = Potassium + current%AqChems(KHSO4i)*current%numberofparticles		
			Env => Env%next
			current => current%next
		END DO
		WRITE(*,*) "Potassium (ug/m3): ", Potassium*AqMolecularMass(PotassiumIndex)*1.0E12

		ALLOCATE(NumConc(NumBins))
		ALLOCATE(InitNumConc(NumBins))
		ALLOCATE(EnvNumConc(NumBins))
		
		!Set heterogeneous Reactions
		HeteroBins = NumBins
		ALLOCATE(HeteroNumberConc(HeteroBins))
		ALLOCATE(HeteroRadii(HeteroBins))

		!Get particle number conc. and radii for hetero. rate calculation
		current => particles(1,1,1)%first
		I = 1
		DO WHILE(associated(current))
			HeteroNumberConc(I) = current%NumberofParticles
			HeteroRadii(I) = current%EffectiveRadius/100.0 !convert from cm to m
			!ResetHeteroRates expects radii in units of m
			current => current%next
			I = I+1
		END DO
				
		!Calculate Heterogeneous rates from size dist. info
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		CALL RecalculatePhotoRates(ut, .TRUE.)
					
		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
		Mgas = GetM   (X(1),X(2),X(3))

		ALLOCATE(OutFH(3+NumBins))
		DO I = 1, 3+Numbins
			OutFH(I) = GetFileHandle()
		END DO
	
		OPEN(UNIT = OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
		WRITE(OutFH(1), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
		WRITE(OutFH(2), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")
		WRITE(OutFH(3), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
	
		DO q = 1, HowManyGasChems
			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE 
				WRITE(OutFH(3),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			END IF
		END DO
		WRITE(OutFH(3),FMT='(A8,A5)', ADVANCE='NO') "ExtCoeff", " "
		WRITE(OutFH(3),FMT='(A8,A5)', ADVANCE='NO') "SingScat", " "
		WRITE(OutFH(3),FMT='(A5,A8)', ADVANCE='NO') "Assym", " "

		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""

		DO I = 4, 3 + NumBins
			OPEN(UNIT = OutFH(I), FILE=TRIM(OutputDeckSubDir)//"Aerosol"//TRIM(INT2STR(I-3))//".txt", STATUS="REPLACE")
			WRITE(OutFH(I), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "Number", " "
			WRITE(OutFH(I),FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
			WRITE(OutFH(I),FMT='(A10, A3)', ADVANCE='NO') "Ionic Str.", " "
			WRITE(OutFH(I),FMT='(A8, A5)', ADVANCE='NO') "H2O Act.", " "
			WRITE(OutFH(I),FMT='(A10, A3)', ADVANCE='NO') "Sol. Dens.", " "
			WRITE(OutFH(I),FMT='(A11, A2)', ADVANCE='NO') "Part. Dens.", " "
			WRITE(OutFH(I),FMT='(A2, A11)', ADVANCE='NO') "Ph", " "
			WRITE(OutFH(I),FMT='(A4, A9)', ADVANCE='NO') "Temp", " "
			DO J = 1, NumAqChems
				WRITE(OutFH(I),FMT='(A13)', ADVANCE='NO') AqPhaseChemicalNames(J)
			END DO
			DO J = 1, NumOrgChems		
				!WRITE(OutFH(I),FMT='(A13)', ADVANCE='NO') OrgPhaseChemicalNames(J)
			END DO
			
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotBC", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotPOA", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotSOA", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotH2O", " "
			WRITE(OutFH(I),FMT='(A4, A9)', ADVANCE='NO') "TotK", " "
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotNa", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotNH4", " "
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotCl", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotNO3", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotSO4", " "
			
			WRITE(OutFH(I), FMT='(A1)') ""
		END DO

		t=0	
		!Write current time in minutes
		WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
		WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
					
		!Write gas phase concentrations in ppb
		DO q = 1, HowManyGasChems
		
			Chem(q) = GridGasChem(1,1,1,q) !molecules/cm3

			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			ELSE 
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			END IF
		END DO
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') 0.0
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') 0.0
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') 0.0

		!Advance to next line
		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""
			
		current => particles(1,1,1)%first
		DO I = 4, 3+NumBins
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature
			!mol/particle
			DO J = 1, NumAqChems
				WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(J)
			END DO
			DO J = 1, NumOrgChems		
				!WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(J)
			END DO
								
			!BC Mass Conc (ug/m3)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(1)*current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
			!POA Mass Conc (ug/m3)
			TotalPOA = 0.
			DO J = 2, 9
				TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			DO J = 20, NumOrgChems
				TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

			!SOA Mass Conc (ug/m3)
			TotalSOA = 0.
			DO J = 10, 19
				TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalSOA = TotalSOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
			!Water Mass Conc. (ug/m3)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Aqchems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1)
					
			!K+ Mass Conc (ug/m3)
			TotalK = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
					TotalK = TotalK + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. K2SO4Index) THEN
					TotalK = TotalK + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalK = TotalK*AqMolecularMass(PotassiumIndex)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
			!Na+ Mass Conc (ug/m3)
			TotalNa = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
					TotalNa = TotalNa + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. Na2SO4i) THEN
					TotalNa = TotalNa + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNa = TotalNa*AqMolecularMass(Nai)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
			!NH4+ Mass Conc (ug/m3)
			TotalNH4 = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
					TotalNH4 = TotalNH4 + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. NH42SO4i) THEN
					TotalNH4 = TotalNH4 + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

			!Cl- Mass Conc (ug/m3)
			TotalCl = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
						TotalCl = TotalCl + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalCl = TotalCl*AqMolecularMass(Cli)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

			!NO3- Mass Conc (ug/m3)
			TotalNO3 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
					TotalNO3 = TotalNO3 + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

			!SO4- Mass Conc (ug/m3)
			TotalSO4 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
				    .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
					TotalSO4 = TotalSO4 + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

			
			WRITE(OutFH(I), FMT='(A1)') ""
					
			IF (associated(current%next)) THEN
				current => current%next
			ELSE
				CYCLE
			END IF

		END DO
		
		!Output size distribution
		CALL SizeDistOutput(1,1,1,t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
	
		Slope = 2.0 
		b = 142 + (6.0/7.0) !Timbavati
			!1.0e50 !No Dilution 
			!142 + (6.0/7.0) !Timbavati
			!66 + (2.0/3.0) !Timbavati Fast
			!400.0 !Timbavati Slow 
			!2500.0 !Timba very slow (Kyy = 100 m2/s, yo = 1km)
			!156.25 !Otavi
			
		ChemTimestep = 1.
		MixTimestep = 1*ChemTimestep !Should always be an integer multiple!
		CoagTimeStep = 60.
		CondTimeStep = 20*ChemTimestep
		
		!Run the main loop with above timesteps
		DO r = 1, INT(61*60/MixTimestep)
			!WRITE(*,*) "In Loop"

			!Change CondTimeStep after first minute
			IF (r .GT. INT(60/Mixtimestep)) CondTimeStep = 20. 
			
			!Calculate Temperature
			X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
			TEMP = GetTemp(X(1),X(2),X(3))

			!NOTE: Right now, pressure and temperature are assumed to be 
			!constant, so we don't have to recalculate the rest of the reaction rates
			!CALL RecalculateReactionRates(IFGAS,IFAQ)
			!WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
			!after this step, or those rates will be set to 0!
			
			!Recalculate Photolysis Rates
			CALL RecalculatePhotoRates(ut, .FALSE.)
			!WRITE(*,*) "Photo Okay"

			!Get particle number conc. and radii for hetero. rate calculation
			current => particles(1,1,1)%first
			I = 1
			DO WHILE(associated(current))
				HeteroNumberConc(I) = current%NumberofParticles
				HeteroRadii(I) = current%EffectiveRadius/100.0 !convert from cm to m
				!ResetHeteroRates expects radii in units of m
				current => current%next
				I = I+1
			END DO
				
			!Calculate Heterogeneous rates from size dist. info
			CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .FALSE.)
			!WRITE(*,*) "Hetero Okay"
						
			!Step Gas Chemistry
			CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep),1,1,1)
			!WRITE(*,*) "StepGasChem Okay"
			
			!Step Condensation
			IF (MOD(t,CondTimeStep) .LT. 0.01 .AND. t .NE. 0.0) THEN
				!Step Condensation, if desired
				IF(.TRUE.) THEN
					CALL RegridAerosol (1, 1, 1)
					CALL FindAqElectrolyteEquilibriumForGridPoint (1, 1, 1, UpdateThermo = .TRUE.) 
					CALL StepCondensationAll (CondTimeStep, 1, 1, 1, 1)
					!Note above calls RegridAerosol automatically
					!WRITE(*,*) t, "Cond Okay"
				END IF
			END IF

			!Step Coagulation
			IF (MOD(t,CoagTimeStep) .LT. 0.01 .AND. t .NE. 0.0) THEN

				!Step Coagulation, if desired
				DoTestCoagulation = .TRUE.
				IF (DoTestCoagulation) THEN	
					CALL StepSectionalCoagulationJacobson(1,1,1,CoagTimeStep)
					!Note that above already calls RegridAerosol
					!WRITE(*,*) t, "Coag Okay"
				END IF
					
			END IF

			!Calculate T_mix
			T_mix = Slope*t+b
		
			!Step Dilution of Chems
			DO q = 1, HowManyGasChems
				Chem(q) = GridGasChem(1,1,1,q)
				Chem(q) = (Chem(q) + (GasPhaseBackground(q)*MixTimestep/T_mix))/(1+ MixTimestep/T_mix + MixTimestep*GasPhaseDepVel(q)/InversionHeight)
			END DO	

			!Update chemical concentrations
			CALL UpdateChemicalConcentrations(1, 1, 1, Chem(1:HowManyEvolveGasChems))

			!Step Dilution and Deposition of Particles
			CALL SortAerosolAtGridPointForCoagulation (1, 1, 1) 
			current => particles(1,1,1)%first
			Env => BoundaryParticles%first
			I = 1
			DO WHILE(associated(current))
    			
				IF(Env%Edges(1) .NE. Current%Edges(1)) THEN
					!WRITE(*,*) Env%Edges(1), Current%Edges(1)
					!CALL ERROR("Env. Aerosol not lined up!")					
				END IF
				CALL RecalculateRadius(current)
				!Get Number Concentration of particles in particles/cm3
				InitNumConc(I) = current%numberofparticles 
				EnvNumConc(I) = Env%numberofparticles
				!Dilute and Dep number of particles (assume no background for now)
				current%numberofparticles = (InitNumConc(I) + (EnvNumConc(I)*MixTimestep/T_mix))/ &
							(1.0 + MixTimestep*(1.0/T_mix + TerminalVelocity(Current)/InversionHeight))
				!WRITE(*,*) I, InitNumConc(I), current%numberofparticles, TerminalVelocity(Current)
				IF(current%Numberofparticles .GT. 0.) THEN
					DO J = 1, NumAqChems
						InitialMass = InitNumConc(I)*current%AqChems(J)
						EnvMass = EnvNumConc(I)*Env%AqChems(J)
						NewMass = (InitialMass + EnvMass*MixTimestep/T_mix)/ &
										  (1+ MixTimestep*(1.0/T_mix + TerminalVelocity(Current)/InversionHeight))
						Current%AqChems(J) = NewMass/current%numberofparticles
						!IF(current%Numberofparticles .LT. 1.0) WRITE(*,*) I,J, Current%AqChems(J)
						!IF(IsNAN(Current%AqChems(J))) WRITE(*,*) "Error in Particle Diluton! ", I, J
					END DO
				
					DO J = 1, NumOrgChems
						InitialMass = InitNumConc(I)*current%OrgChems(J)
						EnvMass = EnvNumConc(I)*Env%OrgChems(J)
						NewMass = (InitialMass + EnvMass*MixTimestep/T_mix)/ &
										  (1+ MixTimestep*(1.0/T_mix + TerminalVelocity(Current)/InversionHeight))
						Current%OrgChems(J) = NewMass/current%numberofparticles
					END DO
				
					DO J = 1, NumAqOrgChems
						InitialMass = InitNumConc(I)*current%AqOrgChems(J)
						EnvMass = EnvNumConc(I)*Env%AqOrgChems(J)
						NewMass = (InitialMass + EnvMass*MixTimestep/T_mix)/ &
										  (1+ MixTimestep*(1.0/T_mix + TerminalVelocity(Current)/InversionHeight))
						Current%AqOrgChems(J) = NewMass/current%numberofparticles
					END DO
				ELSE !No particles
					DO J = 1, NumAqChems
						Current%AqChems(J) = 0.
						!IF(IsNAN(Current%AqChems(J))) WRITE(*,*) "Error in Particle Diluton! ", I, J
					END DO
				
					DO J = 1, NumOrgChems
						Current%OrgChems(J) = 0.
					END DO
				
					DO J = 1, NumAqOrgChems
						Current%AqOrgChems(J) = 0.
					END DO				
				END IF
				!Force particle temperature to be env. Temperature
				!WARNING: This is a kludge, since the program is calculating
				!a huge temperature for the largest particles, and I'm not sure why
				current%Temperature = TEMP
				
				!Recalculate optical parameters
				!WRITE(*,*) "Call Optical"
				CALL ShellRefIndAndRad(Current)
				!WRITE(*,*) "Optical Okay"
				
				I = I + 1
		
				current => current%next
				Env => Env%next
				
			END DO !END STEP PARTICLE DILUTION
			!WRITE(*,*) "Before AOD"
			!WRITE(*,*) "AOD: ", AerosolOpticalDepth (1, 1, 1)	
			!WRITE(*,*) "Dilution Okay"	
			!CALL SpillBeans(1,1,1)
			CALL AerosolOptProp (1, 1, 1, ExtCoeff, SingScat, Assym)
			
			!Write Output to file every minute 
			IF (MOD(t,10.0) .LT. 0.01) THEN
				
				!Write current time in minutes
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
					
				!Write gas phase concentrations in ppb
				DO q = 1, HowManyGasChems
		
					Chem(q) = GridGasChem(1,1,1,q) !molecules/cm3

					IF (q .LE. 70) THEN
						WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
						WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE 
						WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					END IF
				END DO
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') ExtCoeff(1)
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') SingScat(1)
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Assym(1)

				!Advance to next line
				WRITE(OutFH(1), FMT='(A1)') ""
				WRITE(OutFH(2), FMT='(A1)') ""
				WRITE(OutFH(3), FMT='(A1)') ""
			
				GridTotalSOA = 0.0
				current => particles(1,1,1)%first
				DO I = 4, 3+NumBins
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature	
					!mol/particle
					DO J = 1, NumAqChems
						WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(J)
					END DO
					DO J = 1, NumOrgChems		
						!WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(J)
					END DO
					
					!BC Mass Conc (ug/m3)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(1)*current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
					!POA Mass Conc (ug/m3)
					TotalPOA = 0.
					DO J = 2, 9
						TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
						TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
					END DO
					DO J = 20, NumOrgChems
						TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
						TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
					END DO
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

					!SOA Mass Conc (ug/m3)
					TotalSOA = 0.
					DO J = 10, 19
						TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
						TotalSOA = TotalSOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
					END DO
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
					!Water Mass Conc. (ug/m3)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Aqchems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1)
			
					!K+ Mass Conc (ug/m3)
					TotalK = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
							TotalK = TotalK + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. K2SO4Index) THEN
							TotalK = TotalK + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalK = TotalK*AqMolecularMass(PotassiumIndex)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
					!Na+ Mass Conc (ug/m3)
					TotalNa = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
							TotalNa = TotalNa + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. Na2SO4i) THEN
							TotalNa = TotalNa + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNa = TotalNa*AqMolecularMass(Nai)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
					!NH4+ Mass Conc (ug/m3)
					TotalNH4 = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
							TotalNH4 = TotalNH4 + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. NH42SO4i) THEN
							TotalNH4 = TotalNH4 + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

					!Cl- Mass Conc (ug/m3)
					TotalCl = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
							TotalCl = TotalCl + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalCl = TotalCl*AqMolecularMass(Cli)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

					!NO3- Mass Conc (ug/m3)
					TotalNO3 = 0.
					DO J = 1, NumAqChems
						IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
							TotalNO3 = TotalNO3 + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

					!SO4- Mass Conc (ug/m3)
					TotalSO4 = 0.
					DO J = 1, NumAqChems
						IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
						    .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
							TotalSO4 = TotalSO4 + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

					WRITE(OutFH(I), FMT='(A1)') ""

					IF (associated(current%next)) THEN
						current => current%next
					ELSE
						CYCLE
					END IF

				END DO
			END IF !Write Output

			IF (MOD(t,600.0)<=0.05 .OR. ABS(t-540) <=0.05 .OR. ABS(t-2820) <=0.05) THEN
				IF(INT(t) .NE. 0.) THEN
					CALL SizeDistOutput(1,1,1,t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
				END IF
			END IF
			
			!Update time
			t = t + MixTimestep
			ut = ut + MixTimeStep/3600.

			IF (MOD(t,CoagTimeStep)<=0.01) THEN
				WRITE(*,*) trim(real2str(t/60)), " minutes done"
			END IF

		END DO

		
		Dohumiditytest=.TRUE.
		    
	END IF !Size Dist. Particle Tests
	
	IF (.FALSE.) THEN !THERMO TESTS (Compare to ISORROPIA)
		
		IF(NOT(DohumidityTest)) THEN
		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize (1*cm, 1*cm, 1*cm)
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()
	
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL RecalculateReactionRates(IFGAS,IFAQ)
		
		!Set heterogeneous Reactions
		HeteroBins = 1
		ALLOCATE(HeteroNumberConc(HeteroBins))
		ALLOCATE(HeteroRadii(HeteroBins))

		HeteroNumberConc(1) = 1.433e5 !Set for Timbavati Conditions
		HeteroRadii(1) = 1.0e-7
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		
		ALLOCATE(CHEM(HowManyGasChems))

		!! INITIALIZE the DISSOLUTION routine
		CALL SetAllDissolution

		WRITE(*,*) "Before Initialization"
		
		!! INITIALIZE the PARTICLES
		CALL InitializeParticlesSectional	
		
		WRITE(*,*) "After Initialization"
		END IF
		!Spill Beans
		CALL SpillBeans(1,1,1)
                 
		!Count # of bins, # of Aq. Chems
		current => particles(1,1,1)%first
		NumBins = 0
		DO WHILE(associated(current))
			NumBins = NumBins + 1
			NumAqChems = size(current%AqChems)
			current => current%next
		END DO	
		IF(NOT(DohumidityTest)) ALLOCATE(NumConc(NumBins))
		
		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
		Mgas = GetM   (X(1),X(2),X(3))
                

		ALLOCATE(OFH(4))
		OFH(1) = GetFileHandle()
		OFH(4) = GetFileHandle()

		OPEN(UNIT = OFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOutBulk.txt", STATUS="REPLACE")
		WRITE(OFH(1), FMT='(A2,A13)', ADVANCE='NO') "RH", " "
		OPEN(UNIT = OFH(4), FILE=TRIM(OutputDeckSubDir)//"BulkAerosolOut.txt", STATUS="REPLACE")
		WRITE(OFH(4), FMT='(A2,A13)', ADVANCE='NO') "RH", " "
		
                WRITE(*,*) "After File Handles"

		DO q = 1, HowManyGasChems
			IF (q .LE. 70) THEN
				WRITE(OFH(1),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			END IF
		END DO
                WRITE(*,*) "After first write"
                
		
                IF (NumBins .GT. 1) CALL ERROR("The current output routine needs a bulk aerosol to work.")
                
                WRITE(OFH(4),FMT='(A10, A3)', ADVANCE='NO') "Ionic Str.", " "
		WRITE(OFH(4),FMT='(A8, A5)', ADVANCE='NO') "H2O Act.", " "
		WRITE(OFH(4),FMT='(A10, A3)', ADVANCE='NO') "Sol. Dens.", " "
		WRITE(OFH(4),FMT='(A11, A2)', ADVANCE='NO') "Part. Dens.", " "
		WRITE(OFH(4),FMT='(A2, A11)', ADVANCE='NO') "Ph", " "
		WRITE(OFH(4),FMT='(A4, A9)', ADVANCE='NO') "Temp", " "
		WRITE(OFH(4),FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
		!WRITE(OFH(4),FMT='(A10, A3)', ADVANCE='NO') "Scat. Eff.", " "
		!WRITE(OFH(4),FMT='(A12, A1)', ADVANCE='NO') "Scat. Coeff.", " "
		WRITE(OFH(4),FMT='(A11, A2)', ADVANCE='NO') "Ext. Coeff.", " "
		WRITE(OFH(4),FMT='(A10, A3)', ADVANCE='NO') "Scat. Alb.", " "
		WRITE(OFH(4),FMT='(A6, A7)', ADVANCE='NO') "Assym.", " "
		WRITE(OFH(4),FMT='(A12, A1)', ADVANCE='NO') "Ext. Coeff.2", " "
		WRITE(OFH(4),FMT='(A12, A1)', ADVANCE='NO') "Ext. Coeff.3", " "
	        WRITE(OFH(4),FMT='(A12, A1)', ADVANCE='NO') "Ext. Coeff.4", " "
		WRITE(OFH(4),FMT='(A12, A1)', ADVANCE='NO') "Ext. Coeff.5", " "
		WRITE(OFH(4),FMT='(A12, A1)', ADVANCE='NO') "Ext. Coeff.6", " "

		DO I = 1, NumAqChems
			WRITE(OFH(4),FMT='(A13)', ADVANCE='NO') AqPhaseChemicalNames(I)
		END DO
                WRITE(*,*) "After aero write"
                

		WRITE(OFH(1), FMT='(A1)') ""
		WRITE(OFH(4), FMT='(A1)') ""
		Write(*,*) Mgas
		!Spill Beans
		CALL SpillBeans(1,1,1)
                WRITE(*,*) "After Spill Beans"
                

		NewRH = 0.82
		DO J = 1, 26
			NewRH = NewRH - 0.02		
			CALL SetRelativeHumidity (1,1,1,NewRH)
			
			WRITE(*,*) "Before Equilibrate"								
			!CALL EquilibrateInternallyAtGridPoint(1,1,1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
					!Spill Beans
                        CALL EquilibrateGridPointAll(1,1,1, .TRUE., .TRUE., .TRUE.)
			WRITE(*,*) "Equilibrate Okay"
			
			

			!Recalculate optical parameters
			current => particles(1,1,1)%first
			DO I = 1, NumBins
				CALL ShellRefIndAndRad(Current)
				current => current%next
			END DO

			current => particles(1,1,1)%first
				!Write current RH
				WRITE(OFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') NewRH
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') NewRH
				
				!Write gas phase concentrations in ppb
				DO q = 1, HowManyGasChems
		
					Chem(q) = GridGasChem(1,1,1,q) !molecules/cm3

					IF (q .LE. 70) THEN
						WRITE(OFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas !ppb
					END IF
				END DO

				!IF (NumBins .GT. 1) CALL ERROR("The current output routine needs a bulk aerosol to work.") 
	
				CALL AerosolOptProp (1, 1, 1, ExtCoeff, SingScat, Assym)
		
				current => particles(1,1,1)%first
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%EffectiveRadius
				!WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%ScaEff
				!WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%ScaEff*(PI*current%EffectiveRadius**2)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') ExtCoeff(1)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') SingScat(1)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') Assym(1)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') ExtCoeff(2)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') ExtCoeff(3)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') ExtCoeff(4)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') ExtCoeff(5)
				WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') ExtCoeff(6)

				!ug/m3 air
				DO I = 1, NumAqChems
					WRITE(OFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(I)*current%numberofparticles*1.0e12*AqMolecularMass(I)
				END DO
		
				!Advance to next line
				WRITE(OFH(1), FMT='(A1)') ""
				WRITE(OFH(4), FMT='(A1)') ""
			
			!Spill Beans
			CALL SpillBeans(1,1,1)

			Gas = 0
			Aq = 1
			HNO3GasIndex = FindChem("HNO3", Gas)
			!Calculate gas conc of HNO3 in mol/cm3
			HNO3Conc = GridGasChem(1,1,1,HNO3GasIndex) !molec/cm3
			!WRITE(*,*) "HNO3 gas (molecule/cm3): ", trim(real2str(HNO3Conc))
			
			NH3GasIndex = FindChem("NH3", Gas)
			!Calculate gas conc of NH3 in mol/cm3
			NH3Conc = GridGasChem(1,1,1,NH3GasIndex) !molec/cm3
			!WRITE(*,*) "NH3 gas (molecule/cm3): ", trim(real2str(NH3Conc))

			HClGasIndex = FindChem("HCl", Gas)
			!Calculate gas conc of HCl in mol/cm3
			HClConc = GridGasChem(1,1,1,HClGasIndex) !molec/cm3
			!WRITE(*,*) "HCl gas (molecule/cm3): ", trim(real2str(HClConc))

			SOA2GasIndex = FindChem("SOA2", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			!WRITE(*,*) "SOA2 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

			SOA2GasIndex = FindChem("RP13", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			!WRITE(*,*) "RP13 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

			SOA2GasIndex = FindChem("RP17", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			!WRITE(*,*) "RP17 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

			SOA2GasIndex = FindChem("RP18", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			!WRITE(*,*) "RP18 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

			SOA2GasIndex = FindChem("UR26", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			!WRITE(*,*) "UR26 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

			SOA2GasIndex = FindChem("UR29", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			!WRITE(*,*) "UR29 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

			SOA2GasIndex = FindChem("UR30", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			!WRITE(*,*) "UR30 gas (molecule/cm3): ", trim(real2str(SOA2Conc))
			
                        SOA2GasIndex = FindChem("IVOC1", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC1 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

                        SOA2GasIndex = FindChem("IVOC2", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC2 gas (molecule/cm3): ", trim(real2str(SOA2Conc))
                        
                         SOA2GasIndex = FindChem("IVOC3", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC3 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

			 SOA2GasIndex = FindChem("IVOC4", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC4 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

 SOA2GasIndex = FindChem("IVOC5", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC5 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

 SOA2GasIndex = FindChem("IVOC6", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC6 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

 SOA2GasIndex = FindChem("IVOC7", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC7 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

 SOA2GasIndex = FindChem("IVOC8", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC8 gas (molecule/cm3): ", trim(real2str(SOA2Conc))

 SOA2GasIndex = FindChem("IVOC9", Gas)
			!Calculate gas conc of SOA2 in mol/cm3
			SOA2Conc = GridGasChem(1,1,1,SOA2GasIndex) !molec/cm3
			WRITE(*,*) "IVOC9 gas (molecule/cm3): ", trim(real2str(SOA2Conc))


LevAeroIndex = FindChem("(NH4)3H(SO4)2", Aq)
			AmSulfAeroIndex = FindChem("(NH4)2SO4", Aq)
			H2SO4AeroIndex = FindChem("H2SO4", Aq)
			
		END DO
	
	END IF !THERMO TESTS (Compare to ISORROPIA)

	IF (.FALSE.) THEN !COAG TESTS 
		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize (1*cm, 1*cm, 1*cm)
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()
	
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		
		!! INITIALIZE the PARTICLES
		CALL InitializeParticlesSectional

		
		!! Sort the gridcell, sections in order and first
		CALL SortAerosolAtGridPointForCoagulation (1, 1, 1)
		CALL SpillBeans(1,1,1)

		!Count # of bins
		current => particles(1,1,1)%first
		NumBins = 0
		DO WHILE(associated(current))
			NumBins = NumBins + 1
			current => current%next
		END DO	
		ALLOCATE(NumConc(NumBins))

		!Calculate Temperature
		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
			
		ALLOCATE(OFH(4))
		OutFH(4) = GetFileHandle()

		OPEN(UNIT = OutFH(4), FILE=TRIM(OutputDeckSubDir)//"CoagTestOut.txt", STATUS="REPLACE")
		WRITE(OutFH(4), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		
		DO Q = 1, NumBins		
			WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "Number", " "
			WRITE(OutFH(4),FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
		END DO

		WRITE(OutFH(4), FMT='(A1)') ""
			
		!Spill Beans
		CALL RegridAerosol(1, 1, 1)
		CALL SpillBeans(1,1,1)

		t=0
		CALL SizeDistOutput(1,1,1,t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
		MixTimestep = 60
		CoagTimeStep = 60 
		
		current => particles(1,1,1)%first
		I = 1
		TotalNumber = 0.0
		TotalVolume = 0.0
		TotalMass = 0.0
		TotalPotassiumMass = 0.0
		PotassiumIndex = FindChem("K+", 1)
		K2SO4Index = FindChem("K2SO4", 1)
		DO WHILE(associated(current))
			!current%Temperature = TEMP
			TotalNumber = TotalNumber + current%numberofparticles
			TotalVolume = TotalVolume + (4.0/3.0)*Pi*(current%effectiveradius**3)*current%numberofparticles
			IF(current%numberofparticles .GT. 0.0) THEN
				TotalMass = TotalMass + ParticleMass(current)*current%numberofparticles
				TotalPotassiumMass = TotalPotassiumMass + (current%aqchems(PotassiumIndex) + 2*current%aqchems(K2SO4Index)) &
			                      *AqMolecularMass(PotassiumIndex)*current%numberofparticles
			END IF
			I = I + 1
			current => current%next
		END DO	
	
		!WRITE(*,*) "Total Initial Number Conc: ", TotalNumber
		!WRITE(*,*) "Total Initial Volume Conc: ", TotalVolume
		!WRITE(*,*) "Total Initial Particle Mass: ", TotalMass
		!WRITE(*,*) "Total Initial Potassium Mass: ", TotalPotassiumMass
		
		!Run the main loop with above timesteps
		DO r = 1, INT(10*60*60/MixTimestep)
			
			
			
			IF (MOD(t,60.0)<=0.05) THEN !Write every 10 minutes
				!Write current time in minutes
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				
			    CALL SortAerosolAtGridPointForCoagulation (1, 1, 1)
				
				current => particles(1,1,1)%first
				DO WHILE(associated(current))
					WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
					WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius/micron
					current => current%next
				END DO
				
				!Advance to next line
				WRITE(OutFH(4), FMT='(A1)') ""
			END IF
			
			current => particles(1,1,1)%first
			I = 1
			TotalNumber = 0.0
			TotalVolume = 0.0
			TotalPotassiumMass = 0.0
			K2SO4mass = 0.0 
			Kmass = 0.0
			PotassiumIndex = FindChem("K+", 1)
			K2SO4Index = FindChem("K2SO4", 1)
			DO WHILE(associated(current))
				IF (current%numberofparticles .GT. 0.0) THEN
					TotalNumber = TotalNumber + current%numberofparticles
					TotalVolume = TotalVolume + (4.0/3.0)*Pi*(current%effectiveradius**3)*current%numberofparticles
					!WRITE(*,*) I, current%radius/micron, current%aqchems(K2SO4Index)
					K2SO4mass = K2SO4mass + 2*current%aqchems(K2SO4Index)*AqMolecularMass(PotassiumIndex)*current%numberofparticles
					Kmass = Kmass + current%aqchems(PotassiumIndex)*AqMolecularMass(PotassiumIndex)*current%numberofparticles
					TotalPotassiumMass = TotalPotassiumMass + (current%aqchems(PotassiumIndex) + 2*current%aqchems(K2SO4Index)) &
			                      *AqMolecularMass(PotassiumIndex)*current%numberofparticles
				END IF
				I = I + 1
				current => current%next
			END DO 

			!IF(r .EQ. 2) THEN
				!WRITE(*,*) "Total Number Conc: ", TotalNumber
				!WRITE(*,*) "Total Volume Conc: ", TotalVolume
				!WRITE(*,*) "K2SO4 potassium mass", K2SO4Mass
				!WRITE(*,*) "K+ mass", Kmass
				!WRITE(*,*) "Total Potassium Mass: ", TotalPotassiumMass
			!END IF

			!Coagulate
			CALL StepSectionalCoagulationJacobson(1,1,1,CoagTimeStep)
			
			!CALL SpillBeans(1,1,1)
			!Update time
			t = t + MixTimestep

			IF (MOD(t,3600.0)<=0.05) THEN
					CALL SizeDistOutput(1,1,1,t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
			END IF

			IF (MOD(t,60.0)<=0.05) THEN
				!WRITE(*,*) trim(real2str(t/60)), " minutes done"
			END IF

		END DO

		current => particles(1,1,1)%first
		I = 1
		TotalNumber = 0.0
		TotalVolume = 0.0
		TotalPotassiumMass = 0.0
		TotalMass = 0.0
		PotassiumIndex = FindChem("K+", 1)
		K2SO4Index = FindChem("K2SO4", 1)
		DO WHILE(associated(current))
			IF(current%numberofparticles .GT. 0.0) THEN
				TotalNumber = TotalNumber + current%numberofparticles
				TotalVolume = TotalVolume + (4.0/3.0)*Pi*(current%effectiveradius**3)*current%numberofparticles
				TotalPotassiumMass = TotalPotassiumMass + (current%aqchems(PotassiumIndex) + 2*current%aqchems(K2SO4Index)) &
			                      *AqMolecularMass(PotassiumIndex)*current%numberofparticles
				TotalMass = TotalMass + ParticleMass(current)*current%numberofparticles
			END IF

			I = I + 1
				current => current%next
		END DO 

		!Spill Beans
		!CALL SpillBeans(1,1,1)
		!WRITE(*,*) "Total Final Number Conc: ", TotalNumber
		!WRITE(*,*) "Total Final Volume Conc: ", TotalVolume
		!WRITE(*,*) "Total Final Particle Mass: ", TotalMass
		!WRITE(*,*) "Total Final Potassium Mass: ", TotalPotassiumMass
			
	END IF !COAG TESTS

	IF (.FALSE.) THEN !Condensation (flux) TESTS
		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize (1*cm, 1*cm, 1*cm)
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()
	
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL ReadPhotolysisRates
		
		CALL RecalculateReactionRates(IFGAS,IFAQ)
			
		!Calculate Initial Photolysis Rates
		ut = 8.75 !Timbavati Initial Universal Time in hours
		!11.00 !Otavi Initial Universal Time in hours
		CALL RecalculatePhotoRates (ut, .TRUE.)

		ALLOCATE(CHEM(HowManyGasChems))

		!! INITIALIZE the DISSOLUTION routine
		CALL SetAllDissolution

		!! INITIALIZE the PARTICLES	
		CALL InitializeParticlesSectional
		
		!!WARNING Force initial water equilibrium between gas and aerosol
		!!(Uses open system for water eq.)
		CALL EquilibrateInternallyAtGridPoint(1,1,1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
		CALL RegridAerosol (1, 1, 1)
		CALL SpillBeans(1,1,1)

		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
		Mgas = GetM   (X(1),X(2),X(3))
		
		!PressEq = SaturationVaporPressure (TEMP, 1)
		!WRITE(*,*) "Press Eq: ", PressEq
		!WRITE(*,*) "Diff: " , MolecularDiffusionCoefficient(X(1),X(2),X(3),INT(OrganicDissolutionData(1,1)))
		!STOP

		!Find chemical indices for later
		NH3Gas = FindChem("NH3", 0)
		HNO3Gas = FindChem("HNO3", 0)
		HClGas = FindChem("HCl", 0)
		H2SO4Gas = FindChem("H2SO4", 0)

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

		!Count # of bins, # of Aq. Chems, and equilibrate water
		current => particles(1,1,1)%first
		NumBins = 0
		DO WHILE(associated(current))
			NumBins = NumBins + 1
			NumAqChems = size(current%AqChems)
			NumOrgChems = size(current%OrgChems)
			NumAqOrgChems = size(current%AqOrgChems)
			CALL EqAllWater(Current, ReturnType, .TRUE.) !Open system!
			current => current%next
		END DO
		ALLOCATE(NumConc(NumBins))	

		!Set heterogeneous Reactions
		HeteroBins = NumBins
		ALLOCATE(HeteroNumberConc(HeteroBins))
		ALLOCATE(HeteroRadii(HeteroBins))

		!Get particle number conc. and radii for hetero. rate calculation
		current => particles(1,1,1)%first
		I = 1
		DO WHILE(associated(current))
			HeteroNumberConc(I) = current%NumberofParticles
			HeteroRadii(I) = current%EffectiveRadius/100.0 !convert from cm to m
			!ResetHeteroRates expects radii in units of m
			current => current%next
			I = I+1
		END DO
				
		!Calculate Heterogeneous rates from size dist. info
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		CALL RecalculatePhotoRates(ut, .TRUE.)
					
		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
		Mgas = GetM   (X(1),X(2),X(3))

		ALLOCATE(OutFH(3+NumBins))
		DO I = 1, 3+Numbins
			OutFH(I) = GetFileHandle()
		END DO
	
		OPEN(UNIT = OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
		WRITE(OutFH(1), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
		WRITE(OutFH(2), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")
		WRITE(OutFH(3), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
	
		DO q = 1, HowManyGasChems
			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE 
				WRITE(OutFH(3),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			END IF
		END DO

		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""

		DO I = 4, 3 + NumBins
			OPEN(UNIT = OutFH(I), FILE=TRIM(OutputDeckSubDir)//"Aerosol"//TRIM(INT2STR(I-3))//".txt", STATUS="REPLACE")
			WRITE(OutFH(I), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "Number", " "
			WRITE(OutFH(I),FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
			WRITE(OutFH(I),FMT='(A10, A3)', ADVANCE='NO') "Ionic Str.", " "
			WRITE(OutFH(I),FMT='(A8, A5)', ADVANCE='NO') "H2O Act.", " "
			WRITE(OutFH(I),FMT='(A10, A3)', ADVANCE='NO') "Sol. Dens.", " "
			WRITE(OutFH(I),FMT='(A11, A2)', ADVANCE='NO') "Part. Dens.", " "
			WRITE(OutFH(I),FMT='(A2, A11)', ADVANCE='NO') "Ph", " "
			WRITE(OutFH(I),FMT='(A4, A9)', ADVANCE='NO') "Temp", " "
			DO J = 1, NumAqChems
				WRITE(OutFH(I),FMT='(A13)', ADVANCE='NO') AqPhaseChemicalNames(J)
			END DO
			DO J = 1, NumOrgChems		
				!WRITE(OutFH(I),FMT='(A13)', ADVANCE='NO') OrgPhaseChemicalNames(J)
			END DO
			DO J = 1, NumAqOrgChems		
				!WRITE(OutFH(I),FMT='(A13)', ADVANCE='NO') AqOrgPhaseChemicalNames(J)
			END DO
					
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotBC", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotPOA", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotSOA", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotH2O", " "
			WRITE(OutFH(I),FMT='(A4, A9)', ADVANCE='NO') "TotK", " "
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotNa", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotNH4", " "
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotCl", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotNO3", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotSO4", " "
			
			WRITE(OutFH(I), FMT='(A1)') ""
		END DO

		t=0	
		!Write current time in minutes
		WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
		WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
					
		!Write gas phase concentrations in ppb
		DO q = 1, HowManyGasChems
		
			Chem(q) = GridGasChem(1,1,1,q) !molecules/cm3

			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			ELSE 
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			END IF
		END DO
		
		!Advance to next line
		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""
			
		current => particles(1,1,1)%first
		DO I = 4, 3+NumBins
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius/micron
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature
			!mol/particle
			DO J = 1, NumAqChems
				WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(J)
			END DO
			DO J = 1, NumOrgChems		
				!WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(J)
			END DO
			DO J = 1, NumAqOrgChems		
				!WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqOrgchems(J)
			END DO
								
			!BC Mass Conc (ug/m3)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(1)*current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
			!POA Mass Conc (ug/m3)
			TotalPOA = 0.
			DO J = 2, 9
				TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
			END DO
			DO J = 1, 8
				TotalPOA = TotalPOA + current%AqOrgchems(J)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J)
			END DO
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

			!SOA Mass Conc (ug/m3)
			TotalSOA = 0.
			DO J = 10, 19
				TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
			END DO
			DO J = 9, 18
				TotalSOA = TotalSOA + current%AqOrgchems(J)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J)
			END DO
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
			!Water Mass Conc. (ug/m3)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Aqchems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1)
					
			!K+ Mass Conc (ug/m3)
			TotalK = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
					TotalK = TotalK + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. K2SO4Index) THEN
					TotalK = TotalK + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalK = TotalK*AqMolecularMass(PotassiumIndex)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
			!Na+ Mass Conc (ug/m3)
			TotalNa = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
					TotalNa = TotalNa + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. Na2SO4i) THEN
					TotalNa = TotalNa + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNa = TotalNa*AqMolecularMass(Nai)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
			!NH4+ Mass Conc (ug/m3)
			TotalNH4 = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
					TotalNH4 = TotalNH4 + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. NH42SO4i) THEN
					TotalNH4 = TotalNH4 + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

			!Cl- Mass Conc (ug/m3)
			TotalCl = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
						TotalCl = TotalCl + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalCl = TotalCl*AqMolecularMass(Cli)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

			!NO3- Mass Conc (ug/m3)
			TotalNO3 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
					TotalNO3 = TotalNO3 + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

			!SO4- Mass Conc (ug/m3)
			TotalSO4 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
				    .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
					TotalSO4 = TotalSO4 + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

			
			WRITE(OutFH(I), FMT='(A1)') ""
					
			IF (associated(current%next)) THEN
				current => current%next
			ELSE
				CYCLE
			END IF

		END DO
		
		!Output size distribution
		CALL SizeDistOutput(1,1,1,t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
			
		WRITE(*,*) "Initial NH3 : ", trim(real2str(GridGasChem(1,1,1,NH3Gas)*1e9/Mgas)), " ppb"
		WRITE(*,*) "Initial HNO3 : ", trim(real2str(GridGasChem(1,1,1,HNO3Gas)*1e9/Mgas)), " ppb"
		WRITE(*,*) "Initial HCl : ", trim(real2str(GridGasChem(1,1,1,HClGas)*1e9/Mgas)), " ppb"
		WRITE(*,*) "Initial H2SO4 : ", trim(real2str(GridGasChem(1,1,1,H2SO4Gas)*1e9/Mgas)), " ppb"
				
		!Run the main loop with above timesteps
		TimeStepSize = 60.0
		DO r = 1, INT(30*60/TimeStepSize)
			
			!Call Condensation
			!WRITE(*,*) "Before StepCondensationAll"
			!CALL FindAqElectrolyteEquilibriumForGridPoint (1, 1, 1, UpdateThermo = .TRUE.) 
			CALL StepCondensationAll (TimeStepSize, 1, 1, 1, 1)
	
			!PressEq = SaturationVaporPressure (TEMP, 1)
			!P1 = GetGridCellChemBurden (1, 1, 1, INT(OrganicDissolutionData(1,1)))*RstarMB*Temp/Avogadro
			!WRITE(*,*) "Press Diff (mbar): ", P1-PressEq
			!WRITE(*,*) "Molecular Weight: ", GasMolecularMass(INT(OrganicDissolutionData(1,1)))
			!WRITE(*,*) "Diff: " , MolecularDiffusionCoefficient(X(1),X(2),X(3),INT(OrganicDissolutionData(1,1)))
			!STOP
			
			
			!Decrease Temperature
			!X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
			!TEMP = GetTemp(X(1),X(2),X(3))
			!Temp = Temp -0.01*TimeStepSize
			!CALL SetTempField (1, 1, 1, Temp)
			!WRITE(*,*) "Temperature: ", GetTemp(X(1),X(2),X(3))
			!WRITE(*,*) "RH :", GetRelativeHumidity (X(1), X(2), X(3))
			!Mgas = GetM (X(1),X(2),X(3))

			!WRITE(*,*) "Before regrid"
			
			!Move particles across bin boundaries, if necessary!
			!(Also updates temperature)
			CALL RegridAerosol (1, 1, 1)

			!Update time
			t = t + TimeStepSize

			!WRITE(*,*) "Before output"
					
			!Write Output to file every minute 
			IF (MOD(t,60.0) .LT. 0.01 .OR. DissolutionEquilibriumorFlux .EQ. .FALSE.) THEN
				
				!Write current time in minutes
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
					
				!Write gas phase concentrations in ppb
				DO q = 1, HowManyGasChems
		
					Chem(q) = GridGasChem(1,1,1,q) !molecules/cm3

					IF (q .LE. 70) THEN
						WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
						WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE 
						WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					END IF
				END DO
		
				!Advance to next line
				WRITE(OutFH(1), FMT='(A1)') ""
				WRITE(OutFH(2), FMT='(A1)') ""
				WRITE(OutFH(3), FMT='(A1)') ""
			
				GridTotalSOA = 0.0
				current => particles(1,1,1)%first
				DO I = 4, 3+NumBins
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius/micron
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature
					!mol/particle
					DO J = 1, NumAqChems
						WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(J)
					END DO
					DO J = 1, NumOrgChems		
						!WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(J)
					END DO
					DO J = 1, NumAqOrgChems		
						!WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqOrgchems(J)
					END DO
					
					!BC Mass Conc (ug/m3)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(1)*current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
					!POA Mass Conc (ug/m3)
					TotalPOA = 0.
					DO J = 2, 9
						TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
					END DO
					DO J = 1, 8
						TotalPOA = TotalPOA + current%AqOrgchems(J)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J)
					END DO
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

					!SOA Mass Conc (ug/m3)
					TotalSOA = 0.
					DO J = 10, 19
						TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
					END DO
					DO J = 9, 18
						TotalSOA = TotalSOA + current%AqOrgchems(J)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J)
					END DO
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
					!Water Mass Conc. (ug/m3)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Aqchems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1)
			
					!K+ Mass Conc (ug/m3)
					TotalK = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
							TotalK = TotalK + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. K2SO4Index) THEN
							TotalK = TotalK + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalK = TotalK*AqMolecularMass(PotassiumIndex)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
					!Na+ Mass Conc (ug/m3)
					TotalNa = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
							TotalNa = TotalNa + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. Na2SO4i) THEN
							TotalNa = TotalNa + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNa = TotalNa*AqMolecularMass(Nai)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
					!NH4+ Mass Conc (ug/m3)
					TotalNH4 = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
							TotalNH4 = TotalNH4 + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. NH42SO4i) THEN
							TotalNH4 = TotalNH4 + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

					!Cl- Mass Conc (ug/m3)
					TotalCl = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
							TotalCl = TotalCl + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalCl = TotalCl*AqMolecularMass(Cli)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

					!NO3- Mass Conc (ug/m3)
					TotalNO3 = 0.
					DO J = 1, NumAqChems
						IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
							TotalNO3 = TotalNO3 + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

					!SO4- Mass Conc (ug/m3)
					TotalSO4 = 0.
					DO J = 1, NumAqChems
						IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
						    .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
							TotalSO4 = TotalSO4 + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

					WRITE(OutFH(I), FMT='(A1)') ""

					IF (associated(current%next)) THEN
						current => current%next
					ELSE
						CYCLE
					END IF

				END DO
			END IF !Write Output
			
			!WRITE(*,*) "After output"
			
			IF (MOD(t,60.0)<=0.05 .OR. DissolutionEquilibriumorFlux .EQ. .FALSE.) THEN
				IF(DissolutionEquilibriumorFlux .EQ. .TRUE.) THEN
					CALL SizeDistOutput(1,1,1,t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
				ELSE
					CALL SizeDistOutput(1,1,1,t,"SizeDist2.txt")
				END IF
			END IF

			!WRITE(*,*) "After size dist"
			
			!IF (MOD(t,60.0)<=0.01) THEN
				WRITE(*,*) trim(real2str(t/60)), " minutes done"
			!END IF
			!If just doing equilibrium, exit here
			IF(DissolutionEquilibriumOrFlux .EQ. .False.) EXIT
		
		END DO

		CALL SpillBeans(1,1,1)
		WRITE(*,*) "Final NH3 : ", trim(real2str(GridGasChem(1,1,1,NH3Gas)*1e9/Mgas)), " ppb"
		WRITE(*,*) "Final HNO3 : ", trim(real2str(GridGasChem(1,1,1,HNO3Gas)*1e9/Mgas)), " ppb"
		WRITE(*,*) "Final HCl : ", trim(real2str(GridGasChem(1,1,1,HClGas)*1e9/Mgas)), " ppb"
		WRITE(*,*) "Final H2SO4 : ", trim(real2str(GridGasChem(1,1,1,H2SO4Gas)*1e9/Mgas)), " ppb"
		WRITE(*,*) "Final POA1 : ", trim(real2str(GridGasChem(1,1,1,INT(OrganicDissolutionData(1,1)))*1e9/Mgas)), " ppb"
		    
	END IF !Condensation (flux) Tests
	
	IF (.FALSE.) THEN !MOVING CENTER TESTS 
		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize (1*cm, 1*cm, 1*cm)
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()

		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL RecalculateReactionRates(IFGAS,IFAQ)
		
		!Set heterogeneous Reactions
		HeteroBins = 1
		ALLOCATE(HeteroNumberConc(HeteroBins))
		ALLOCATE(HeteroRadii(HeteroBins))

		HeteroNumberConc(1) = 1.433e5 !Set for Timbavati Conditions
		HeteroRadii(1) = 1.0e-7
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		
		ALLOCATE(CHEM(HowManyGasChems))

		!! INITIALIZE the DISSOLUTION routine
		CALL SetDissolution
		
		!! INITIALIZE the PARTICLES
		CALL InitializeParticlesSectional
			
		!Spill Beans
		CALL SpillBeans(1,1,1)

		DO r = 1, 100
			current => particles(1,1,1)%first
			NumAqChems = size(current%AqChems)
			i = 1
			DO WHILE(associated(current))
				IF(I .EQ. 1) THEN
					I = I+ 1
					current => current%next
					CYCLE !don't shrink smallest
				END IF
				!WRITE(*,*) NumAqChems
				!Shrink Particles
				DO J = 2, NumAqChems
					current%AqChems(J) = 0.5*current%AqChems(J)
				END DO
				DO J = 1, 24
						current%OrgChems(J) = 0.5*current%OrgChems(J)
				END DO
				DO J = 1, 23
					current%AqOrgChems(J) = 0.5*current%AqOrgChems(J)
				END DO
				CALL FindElectrolyteEquilibrium(Current, UpdateThermo=.TRUE., FirstEquilibration=.FALSE., ReturnType=ReturnType)
				CALL EqAllWater(Current, ReturnType, OpenSystem=.TRUE.)
			
				
				current => current%next
			END DO 
		
			
			!Equilibrate Aq. Phase	
			!CALL FindAqElectrolyteEquilibriumForGridPoint (1, 1, 1, .TRUE.)

			!Resize Particles
			CALL RegridAerosol (1, 1, 1)
		END DO

		CALL RegridAerosol (1, 1, 1)

		!Spill Beans
		CALL SpillBeans(1,1,1)
			
	END IF !MOVING CENTER TESTS


	IF (.FALSE.) THEN !Chamber TESTS
		!!!!!!!!!!!!!!!!!!!!!!!!!
		!! SET the DOMAIN SIZE !!
		CALL SetDomainSize (1*cm, 1*cm, 1*cm)
	
		!! Call subroutines to set module data
		CALL SetFileHandleCounter(12)
		CALL ReadMainInputDeck()

		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL ReadPhotolysisRates
		
		CALL RecalculateReactionRates(IFGAS,IFAQ)
		
		ALLOCATE(CHEM(HowManyGasChems))

		!! INITIALIZE the DISSOLUTION routine
		CALL SetDissolution
		CALL SetOrganicDissolution
		!WRITE(*,*) "Dissolution Okay"
		!STOP

		!! INITIALIZE the PARTICLES
		CALL InitializeParticlesSectional	
		
		!WRITE(*,*) "Check 1"
		!!WARNING Force initial equilibrium between gas and aerosol
		!!(Uses open system for water eq.)
		!CALL EquilibrateGridPointAll(1,1,1, .TRUE., .TRUE., .TRUE.)
		CALL EquilibrateHydrophobicPhaseAtGridPoint(1,1,1)
		!WRITE(*,*) "Check 2"
		
		CALL SpillBeans(1,1,1)

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
		
		!Count # of bins, # of Aq. Chems
		current => particles(1,1,1)%first
		NumBins = 0
		DO WHILE(associated(current))
			NumBins = NumBins + 1
			NumAqChems = size(current%AqChems)
			NumOrgChems = size(current%OrgChems)
			current => current%next
		END DO
		ALLOCATE(NumConc(NumBins))
		
		X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
		TEMP = GetTemp(X(1),X(2),X(3))
		Mgas = GetM   (X(1),X(2),X(3))

		ALLOCATE(OutFH(4))
		OutFH(1) = GetFileHandle()
		OutFH(2) = GetFileHandle()
		OutFH(3) = GetFileHandle()
		OutFH(4) = GetFileHandle()

		OPEN(UNIT = OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
		WRITE(OutFH(1), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
		WRITE(OutFH(2), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")
		WRITE(OutFH(3), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT = OutFH(4), FILE=TRIM(OutputDeckSubDir)//"BulkAerosolOut.txt", STATUS="REPLACE")
		WRITE(OutFH(4), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		
		DO q = 1, HowManyGasChems
			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE 
				WRITE(OutFH(3),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			END IF
		END DO

		IF (NumBins .GT. 1) CALL ERROR("The current output routine needs a bulk aerosol to work.") 
	
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "Number", " "
		WRITE(OutFH(4),FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
		WRITE(OutFH(4),FMT='(A10, A3)', ADVANCE='NO') "Ionic Str.", " "
		WRITE(OutFH(4),FMT='(A8, A5)', ADVANCE='NO') "H2O Act.", " "
		WRITE(OutFH(4),FMT='(A10, A3)', ADVANCE='NO') "Sol. Dens.", " "
		WRITE(OutFH(4),FMT='(A11, A2)', ADVANCE='NO') "Part. Dens.", " "
		WRITE(OutFH(4),FMT='(A2, A11)', ADVANCE='NO') "Ph", " "
		WRITE(OutFH(4),FMT='(A4, A9)', ADVANCE='NO') "Temp", " "
		DO I = 1, NumAqChems
			!WRITE(OutFH(4),FMT='(A13)', ADVANCE='NO') AqPhaseChemicalNames(I)
		END DO
		DO I = 1, NumOrgChems		
					WRITE(OutFH(4),FMT='(A13)', ADVANCE='NO') OrgPhaseChemicalNames(I)
		END DO
		
		WRITE(OutFH(4),FMT='(A5, A8)', ADVANCE='NO') "TotBC", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotPOA", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotSOA", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "Yield ", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotH2O", " "
		WRITE(OutFH(4),FMT='(A4, A9)', ADVANCE='NO') "TotK", " "
		WRITE(OutFH(4),FMT='(A5, A8)', ADVANCE='NO') "TotNa", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotNH4", " "
		WRITE(OutFH(4),FMT='(A5, A8)', ADVANCE='NO') "TotCl", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotNO3", " "
		WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "TotSO4", " "
					
		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""
		WRITE(OutFH(4), FMT='(A1)') ""

		ChemTimestep = 0.5
		MixTimestep = 120.0*ChemTimestep !Should always be an integer multiple!
		!Run the main loop with above timesteps
		DO r = 1, INT(72*60*60/MixTimestep+1)

			!Write Output to file every minute in units of ppb
			IF (MOD(t,3600.0)<=0.05) THEN
				
				!Write current time in minutes
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') t/3600
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') t/3600
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') t/3600
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') t/3600
				
				!Write gas phase concentrations in ppb
				DO q = 1, HowManyGasChems
		
					Chem(q) = GridGasChem(1,1,1,q) !molecules/cm3

					IF (q .LE. 70) THEN
						WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
						WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE 
						WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					END IF
				END DO

				
				IF (NumBins .GT. 1) CALL ERROR("The current output routine needs a bulk aerosol to work.") 
	
				current => particles(1,1,1)%first
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature
				!ug/m3 air
				DO I = 1, NumAqChems
					!WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(I)*current%numberofparticles*1.0e12*AqMolecularMass(I)
				END DO
				DO I = 1, NumOrgChems		
					WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(I)*current%numberofparticles*1.0e12*OrgMolecularMass(I)
				END DO

				!BC Mass Conc (ug/m3)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(1)*current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
				!POA Mass Conc (ug/m3)
				TotalPOA = 0.
				DO J = 2, 9
					TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				END DO
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

				!SOA Mass Conc (ug/m3)
				TotalSOA = 0.
				DO J = 10, NumOrgChems
					TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				END DO
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA

                                Yield = Yield + TotalSOA
                                !!Extra code for Jimenez test, MJA 11-10-2010
                                !!Holds OA constant by setting SOA to 0                               
                                DO J = 10, NumOrgChems
					current%Orgchems(J)=0
				END DO
				!! End ExtraCode
                                WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') Yield
                                

	
				!Water Mass Conc. (ug/m3)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%Aqchems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1)
					
				!K+ Mass Conc (ug/m3)
				TotalK = 0.		
				DO J = 1, NumAqChems
					IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
						TotalK = TotalK + current%Aqchems(J)*current%numberofparticles*1.0e12
					ELSE IF (J .EQ. K2SO4Index) THEN
						TotalK = TotalK + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalK = TotalK*AqMolecularMass(PotassiumIndex)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
				!Na+ Mass Conc (ug/m3)
				TotalNa = 0.		
				DO J = 1, NumAqChems
					IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
						TotalNa = TotalNa + current%Aqchems(J)*current%numberofparticles*1.0e12
					ELSE IF (J .EQ. Na2SO4i) THEN
						TotalNa = TotalNa + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalNa = TotalNa*AqMolecularMass(Nai)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
				!NH4+ Mass Conc (ug/m3)
				TotalNH4 = 0.		
				DO J = 1, NumAqChems
					IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
						TotalNH4 = TotalNH4 + current%Aqchems(J)*current%numberofparticles*1.0e12
					ELSE IF (J .EQ. NH42SO4i) THEN
						TotalNH4 = TotalNH4 + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

				!Cl- Mass Conc (ug/m3)
				TotalCl = 0.		
				DO J = 1, NumAqChems
					IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
						TotalCl = TotalCl + current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalCl = TotalCl*AqMolecularMass(Cli)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

				!NO3- Mass Conc (ug/m3)
				TotalNO3 = 0.
				DO J = 1, NumAqChems
					IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
						TotalNO3 = TotalNO3 + current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

				!SO4- Mass Conc (ug/m3)
				TotalSO4 = 0.
				DO J = 1, NumAqChems
					IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
						.OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
						TotalSO4 = TotalSO4 + current%Aqchems(J)*current%numberofparticles*1.0e12
					END IF
				END DO
				!Convert from umol/m3 to ug/m3
				TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
				WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

				!Advance to next line
				WRITE(OutFH(1), FMT='(A1)') ""
				WRITE(OutFH(2), FMT='(A1)') ""
				WRITE(OutFH(3), FMT='(A1)') ""
				WRITE(OutFH(4), FMT='(A1)') ""
			END IF !Write Output
		
			!Calculate Temperature
			X    = LocateGridPoint(1,1,1,GasPhaseChemGridIndex)	
			TEMP = GetTemp(X(1),X(2),X(3))

			!NOTE: Right now, pressure and temperature are assumed to be 
			!constant, so we don't recalculate the rest of the reaction rates
			!CALL RecalculateReactionRates(IFGAS,IFAQ)
			!WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
			!after this step, or those rates will be set to 0!
			
			!Step Gas Chemistry
			
			
			UR21i = FindChem("UR21", 0)
			!WRITE(*,*) "Before Chem", GridGasChem(1,1,1,UR21i)
			CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep),1,1,1)			
			!WRITE(*,*) "After Chem", GridGasChem(1,1,1,UR21i)
			!Only Equilibrate every 60 s (one second before the minute)
			!IF (MOD(t,CondTimeStep) .LT. 59.01 .AND. MOD(t,CondTimeStep) .GT. 58.99) THEN

				!Equilibrate Gases and particles
				!CALL EquilibrateGridPointAll(1,1,1, .TRUE., .TRUE., .TRUE.)
				!WRITE(*,*) "Before Eq", GridGasChem(1,1,1,UR21i)
				CALL EquilibrateHydrophobicPhaseAtGridPoint(1,1,1)
				!WRITE(*,*) "After Eq", GridGasChem(1,1,1,UR21i)
			!END IF
								
			!Update time
			t = t + MixTimestep

			IF (MOD(t,60.0)<=0.05) THEN
				!WRITE(*,*) trim(real2str(t/60)), " minutes done"
			END IF


                        !!Extra code for Jimenez test, MJA 11-10-2010
                        


		END DO

		!Spill Beans
		CALL SpillBeans(1,1,1)
		    
	END IF !Chamber Tests

	IF(.FALSE.) THEN !ZSR Values
		
	  DO I = 1, 9
		RH = I * 10.0	
		Q = 7 + 6
		Temp = 298.15
		 
		tiny = 1.e-06
		x1 = 0.0001
		x3 = 0.9999

		f1 = dev(x1, Q, RH, Temp)
		f3 = dev(x3, Q, RH, Temp)
		
10		chg = x3-x1

		x2 = (x1+x3)/2.
		f2 = dev(x2, Q, RH, Temp)

		if (chg .lt. tiny) then
			xs = x2
			goto 333
		endif

		if (abs(f2) .lt. tiny) then
         xs = x2
         goto 333
		endif

		if (f1*f2 .lt. 0.) then
         x3 = x2
         f3 = f2
		else
         if (f2*f3 .lt. 0.) then
             x1 = x2
             f1 = f2
         else
            print*,'no root'
         endif
		endif

		goto 10

  333	continue

		print*,xs
		molalbin = xs/18./(1.-xs)
		print*,RH, molalbin

  222	continue
	  END DO		

		
		
		
			!Call UNIFACA, placing activity coefficients in GAMMA(0:N-1) 
			   !OrgMassFrac	 !Water Activity		
			!WRITE(*,*) XPASS(Q)/(1-XPASS(Q))/18.0, Gamma(N-1)*XPASS(N-1) 
			!WRITE(*,*) "Mole frac: ", XPASS(Q), XPASS(N-1)
			!WRITE(*,*) "Gamma: ", Gamma(Q), Gamma(N-1)
			!WRITE(*,*) XPASS(Q), Gamma(N-1), Gamma(Q)
			!XPASS(Q) = XPASS(Q) + 0.01
			!XPASS(N-1) = XPASS(N-1) - 0.01
	
	END IF
	
	!ii = etime(etimearray)

	CALL TRANSCRIPT("")
	CALL TRANSCRIPT("Total Execution Time was "//trim(real2str(ii)))

END PROGRAM ASP
