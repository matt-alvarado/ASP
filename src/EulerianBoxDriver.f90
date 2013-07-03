!! 11/06  2012  Matt Alvarado Updated intitalization to make parallel to others in
!!                            DevelopmentDriver.f90
!! 07/03  2013  Matt Alvarado Major Rewrite to make it work.
SUBROUTINE EulerianBox()
	!! Include modules
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
				K2SO4Mass, Kmass, TotalMass, TimeStepSize, PressEq, P1, EnvMass, InitialMass, NewMass, &
				Xlength, Ylength, Height, FireCO, FluxCO, DeltaCO, spinup

	REAL*4  :: etimearray(2)
	REAL*8, ALLOCATABLE	::	XPASS(:), GAMMA(:), Chem(:), NumConc(:), &
							HeteroNumberConc(:), HeteroRadii(:), &
							InitNumConc(:), EnvNumConc(:), Flux(:), SourceNumConc(:), SourceMassConc(:,:), &
							AqMassFlux(:,:), OrgMassFlux(:,:), AqOrgMassFlux(:,:), NumFlux(:)
        REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym,BackScatCoeff,SubExtCoeff, SubSSA, &
                                          ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                                          SubExtCoeffSep, SubSSASep
        INTEGER, ALLOCATABLE :: OutFH(:), OFH(:)

							
	
	LOGICAL :: DoTestCoagulation, DoHumidityTest
	
	TYPE(Particle), POINTER :: Current, Env

	!Molar emision ratios to CO for Timbavati Fire
	!See Trentman, 2005 
	real*8, ALLOCATABLE :: Gas_ER(:) 
        real*8, ALLOCATABLE :: Aero_Num_ER(:), Aero_AqMass_ER(:,:)
        real*8, ALLOCATABLE :: Aero_OrgMass_ER(:,:), Aero_AqOrgMass_ER(:,:)


        spinup = 600 !spinup time in s
	Xlength = boxlength !in km
	Ylength = boxwidth !in km
	Height = inversionheight !in km
		
	FireCO = fire_co_emis !kg/s
	!Want molecules/m2/s
	FluxCO = FireCO/XLength/Ylength !kg/km2/s
	FluxCO = FluxCO*6.022e23*1000/28.0 !molecules/km2/s
	
	! INITIALIZE THE CHEMISTRY 
	CALL SetChemistryParams

        ! INITIALIZE THE REACTION RATES
	CALL RecalculateReactionRates(IFGAS,IFAQ)			
 	CALL ReadPhotolysisRates
	ut = UTCStartTime
	CALL RecalculatePhotoRates (ut, .TRUE.)

	ALLOCATE(CHEM(HowManyGasChems))
	ALLOCATE(Gas_ER(HowManyGasChems))

        !!This subroutine (a) calculates the ERs for all gasses 
        !!and places them in Gas_ER(:) and (b) sets the initial
        !!gas concentrations in the box to the background values
        CALL InitializeGasEulerianBox(HowManyGasChems, Gas_ER) 
	
	!! INITIALIZE the DISSOLUTION routine
	CALL SetAllDissolution

	!! INITIALIZE the PARTICLES	
	CALL InitializeParticlesSectional
	CALL RegridAerosol ()
	CALL SpillBeans()

	!! Call Boundary distributions
	CALL ReadBoundaryDistributionsOrganic
	CALL PopulateParticlesSectionsRightAway(.TRUE.)
	CALL SortBoundaryPart
	CALL SpillBeans(Env=.TRUE.)
		
	!!WARNING Force initial water equilibrium between gas and aerosol
	!!(Uses open system for water eq.)
	CALL EquilibrateInternallyAtGridPoint(EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
	
	CALL RegridAerosol ()
	CALL SpillBeans()
		
	TEMP = GetTemp()
	Mgas = GetM   ()

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

		!Count # of bins, # of Aq., Org, and AqOrg Chems
		current => particles%first
		NumBins = 0
		DO WHILE(associated(current))
			NumBins = NumBins + 1
			NumAqChems = size(current%AqChems)
			NumOrgChems = size(current%OrgChems)
			NumAqOrgChems = size(current%AqOrgChems)
			current => current%next
		END DO
	




		ALLOCATE(NumConc(NumBins))
		ALLOCATE(InitNumConc(NumBins))
		ALLOCATE(EnvNumConc(NumBins))
		
		!Set heterogeneous Reactions
		HeteroBins = NumBins
		ALLOCATE(HeteroNumberConc(HeteroBins))
		ALLOCATE(HeteroRadii(HeteroBins))

		!Calculate gas fluxes (molecules/km2/s)
		ALLOCATE(Flux(HowManyGasChems))
		DO q = 1, HowManyGasChems
			Flux(q) = 0
		END DO
		
		q = FindChem("CO", 0)
		Flux(q) = FluxCO
		WRITE(*,*) Flux(q), q
		



		ALLOCATE(sourcenumconc(NumBins))
		ALLOCATE(sourcemassconc(NumBins, 44))
		sourcenumconc(1) = 14.555 !number/cm3/ppbm CO, 
		sourcemassconc(1,1) = 5.259e-6 !ppbm H2O/ppbm CO
		sourcemassconc(1,3) = 8.747e-15 !ppbm NH3/ppbm CO
		sourcemassconc(1,6) = 2.655e-5 !ppbm NH4Cl/ppbm CO
		sourcemassconc(1,7) = 2.065e-6 !ppbm NH4NO3/ppbm CO
		sourcemassconc(1,8) = 8.849e-7 !NH42SO4 ppbm/ppbm CO
		sourcemassconc(1,10) = 1.239e-5 !ppbm KCl/ppbm CO
		sourcemassconc(1,12) = 4.140e-6 !ppbm K2SO4/ppbm CO
		sourcemassconc(1,24) = 4.107e-6 !ppbm BC/ppbm CO
		sourcemassconc(1,25) = 1.133e-5 !ppbm POA1/ppbm CO
		sourcemassconc(1,44) = 3.265e-5 !ppbm CBIO/ppbm CO
	
		sourcenumconc(2) = 6.571 !number/cm3/ppbm CO, 
		sourcemassconc(2,1) = 7.217e-4 !ppbm H2O/ppbm CO
		sourcemassconc(2,3) = 1.578e-12 !ppbm NH3/ppbm CO
		sourcemassconc(2,6) = 3.211e-3 !ppbm NH4Cl/ppbm CO
		sourcemassconc(2,7) = 2.498e-4 !ppbm NH4NO3/ppbm CO
		sourcemassconc(2,8) = 1.070e-4 !NH42SO4 ppbm/ppbm CO
		sourcemassconc(2,10) = 1.499e-3 !ppbm KCl/ppbm CO
		sourcemassconc(2,12) = 5.008e-4 !ppbm K2SO4/ppbm CO
		sourcemassconc(2,24) = 4.968e-4 !ppbm BC/ppbm CO
		sourcemassconc(2,25) = 1.371e-3 !ppbm POA1/ppbm CO
		sourcemassconc(2,44) = 3.950e-3 !ppbm CBIO/ppbm CO

		sourcenumconc(3) = 0.8053 !number/cm3/ppbm CO, 
		sourcemassconc(3,1) = 2.302e-3 !ppbm H2O/ppbm CO
		sourcemassconc(3,3) = 4.951e-12 !ppbm NH3/ppbm CO
		sourcemassconc(3,6) = 0.010 !ppbm NH4Cl/ppbm CO
		sourcemassconc(3,7) = 7.836e-4 !ppbm NH4NO3/ppbm CO
		sourcemassconc(3,8) = 3.358e-4 !NH42SO4 ppbm/ppbm CO
		sourcemassconc(3,10) = 4.702e-3 !ppbm KCl/ppbm CO
		sourcemassconc(3,12) = 1.571e-3 !ppbm K2SO4/ppbm CO
		sourcemassconc(3,24) = 1.558e-3 !ppbm BC/ppbm CO
		sourcemassconc(3,25) = 4.301e-3 !ppbm POA1/ppbm CO
		sourcemassconc(3,44) = 0.012 !ppbm CBIO/ppbm CO

		
		!Set Initial Aerosol Concentrations to background
		CALL SortAerosolAtGridPointForCoagulation () 
		current => particles%first
		Env => BoundaryParticles%first
		I = 1
		DO WHILE(associated(current))
    			
			IF(Env%Edges(1) .NE. Current%Edges(1)) THEN
				!WRITE(*,*) Env%Edges(1), Current%Edges(1)
				!CALL ERROR("Env. Aerosol not lined up!")					
			END IF

			!Get Number Concentration of particles in particles/cm3
			current%numberofparticles = Env%numberofparticles
			DO J = 1, NumAqChems	
				Current%AqChems(J) = Env%AqChems(J)
			END DO
			DO J = 1, NumOrgChems	
				Current%OrgChems(J) = Env%OrgChems(J)
			END DO	
			DO J = 1, NumAqOrgChems	
				Current%AqOrgChems(J) = Env%AqOrgChems(J)
			END DO	
			CALL RecalculateRadius(current)
			I = I + 1
		
			current => current%next
			Env => Env%next	
		END DO
	
		TEMP = GetTemp()
		Mgas = GetM   ()
		deltaCO = FireCO/Xlength/Ylength
	    deltaCO = 1.0e9*deltaCO/(GetAirDensity()*1000) !ppbm CO*m/s

		ALLOCATE(NumFlux(NumBins))
		DO q = 1, NumBins
			NumFlux(q) = 0.
		END DO
		!particles/m2/s
		NumFlux(1) = sourcenumconc(1)*deltaCO*1.0e6
		NumFlux(2) = sourcenumconc(2)*deltaCO*1.0e6
		NumFlux(3) = sourcenumconc(3)*deltaCO*1.0e6

		ALLOCATE(AqMassFlux(NumBins,NumAqChems))
		DO q = 1, NumBins
			DO j = 1, NumAqChems
				AqMassFlux(q,j) = 0.
			END DO
		END DO
		!mol/m2/s
		q = FindChem("H2O",1)
		AqMassFlux(1,q) = sourcemassconc(1,1)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(2,q) = sourcemassconc(2,1)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(3,q) = sourcemassconc(3,1)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)

		q = FindChem("NH3",1)
		AqMassFlux(1,q) = sourcemassconc(1,3)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(2,q) = sourcemassconc(2,3)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(3,q) = sourcemassconc(3,3)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)

		q = FindChem("NH4Cl",1)
		AqMassFlux(1,q) = sourcemassconc(1,6)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(2,q) = sourcemassconc(2,6)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(3,q) = sourcemassconc(3,6)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)

		q = FindChem("NH4NO3",1)
		AqMassFlux(1,q) = sourcemassconc(1,7)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(2,q) = sourcemassconc(2,7)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(3,q) = sourcemassconc(3,7)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)

		q = FindChem("(NH4)2SO4",1)
		AqMassFlux(1,q) = sourcemassconc(1,8)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(2,q) = sourcemassconc(2,8)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(3,q) = sourcemassconc(3,8)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)

		q = FindChem("KCl",1)
		AqMassFlux(1,q) = sourcemassconc(1,10)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(2,q) = sourcemassconc(2,10)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(3,q) = sourcemassconc(3,10)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)

		q = FindChem("K2SO4",1)
		AqMassFlux(1,q) = sourcemassconc(1,12)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(2,q) = sourcemassconc(2,12)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)
		AqMassFlux(3,q) = sourcemassconc(3,12)*deltaCO*1.0e-3*GetAirDensity( )/AqMolecularMass(q)

		ALLOCATE(OrgMassFlux(NumBins,NumOrgChems))
		DO q = 1, NumBins
			DO j = 1, NumOrgChems
				OrgMassFlux(q,j) = 0.
			END DO
		END DO

		q = FindChem("BC",2)
		OrgMassFlux(1,q) = sourcemassconc(1,24)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)
		OrgMassFlux(2,q) = sourcemassconc(2,24)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)
		OrgMassFlux(3,q) = sourcemassconc(3,24)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)

		q = FindChem("POA1",2)
		OrgMassFlux(1,q) = sourcemassconc(1,25)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)
		OrgMassFlux(2,q) = sourcemassconc(2,25)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)
		OrgMassFlux(3,q) = sourcemassconc(3,25)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)

		q = FindChem("CBIO",2)
		OrgMassFlux(1,q) = sourcemassconc(1,44)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)
		OrgMassFlux(2,q) = sourcemassconc(2,44)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)
		OrgMassFlux(3,q) = sourcemassconc(3,44)*deltaCO*1.0e-3*GetAirDensity( )/OrgMolecularMass(q)

		ALLOCATE(AqOrgMassFlux(NumBins,NumAqOrgChems))
		DO q = 1, NumBins
			DO j = 1, NumAqOrgChems
				AqOrgMassFlux(q,j) = 0.
			END DO
		END DO

		!Get particle number conc. and radii for hetero. rate calculation
		current => particles%first
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
					
		ALLOCATE(OutFH(3+NumBins))
		DO I = 1, 3+Numbins
			OutFH(I) = GetFileHandle()
		END DO
	
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
		
			Chem(q) = GridGasChem(q) !molecules/cm3

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
			
		current => particles%first
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
		CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
			
		ChemTimestep = 1.
		MixTimestep = 1*ChemTimestep !Should always be an integer multiple!
		CoagTimeStep = 60.
		CondTimeStep = 20*ChemTimestep
		
		!Run the main loop with above timesteps
		DO r = 1, INT(71*60/MixTimestep)
			!WRITE(*,*) "In Loop"

			!Change CondTimeStep after first minute
			IF (r .GT. INT(60/Mixtimestep)) CondTimeStep = 20. 
			
			!Calculate Temperature
			TEMP = GetTemp()

			!NOTE: Right now, pressure and temperature are assumed to be 
			!constant, so we don't have to recalculate the rest of the reaction rates
			!CALL RecalculateReactionRates(IFGAS,IFAQ)
			!WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
			!after this step, or those rates will be set to 0!
			
			!Recalculate Photolysis Rates
			CALL RecalculatePhotoRates(ut, .FALSE.)
			!WRITE(*,*) "Photo Okay"

			!Get particle number conc. and radii for hetero. rate calculation
			current => particles%first
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
			CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep))
			!WRITE(*,*) "StepGasChem Okay"
			
			!Step Condensation
			IF (MOD(t,CondTimeStep) .LT. 0.01 .AND. t .NE. 0.0) THEN
				!Step Condensation, if desired
				IF(.TRUE.) THEN
					CALL RegridAerosol ()
					CALL FindAqElectrolyteEquilibriumForGridPoint (UpdateThermo = .TRUE.) 
					CALL StepCondensationAll (CondTimeStep, 1)
					!Note above calls RegridAerosol automatically
					!WRITE(*,*) t, "Cond Okay"
				END IF
			END IF

			!Step Coagulation
			IF (MOD(t,CoagTimeStep) .LT. 0.01 .AND. t .NE. 0.0) THEN

				!Step Coagulation, if desired
				DoTestCoagulation = .TRUE.
				IF (DoTestCoagulation) THEN	
					CALL StepSectionalCoagulationJacobson(CoagTimeStep)
					!Note that above already calls RegridAerosol
					!WRITE(*,*) t, "Coag Okay"
				END IF
					
			END IF

		
			!Step Dilution of Chems
			!Flux should have units of molecules/m2/s
			DO q = 1, HowManyGasChems
				Chem(q) = GridGasChem(q)
				IF(t .LT. spinup) THEN
					Chem(q) = (Chem(q) + (WindSpeed*GasPhaseBackground(q)*MixTimestep/Xlength)) &
				          /(1+ (MixTimestep*GasPhaseDepVel(q)/(Height*100)) + (MixTimestep*Windspeed/Xlength))
		
				ELSE
					Chem(q) = (Chem(q) + (WindSpeed*GasPhaseBackground(q)*MixTimestep/Xlength) + (Flux(q)*MixTimestep/Height/1.0e6)) &
				          /(1+ (MixTimestep*GasPhaseDepVel(q)/(Height*100)) + (MixTimestep*Windspeed/Xlength))
				END IF
			END DO	
			!Update chemical concentrations
			CALL UpdateChemicalConcentrations(Chem(1:HowManyEvolveGasChems))

			!Step Dilution and Deposition of Particles
			CALL SortAerosolAtGridPointForCoagulation () 
			current => particles%first
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
				!Flux should have units of particles/m2/s
				IF(t .LT. spinup) THEN
					current%numberofparticles = (InitNumConc(I) + (WindSpeed*EnvNumConc(I)*MixTimestep/Xlength)) &
							/(1.0 + MixTimestep*(TerminalVelocity(Current)/(Height*100) + Windspeed/Xlength))
				ELSE
					current%numberofparticles = (InitNumConc(I) + (WindSpeed*EnvNumConc(I)*MixTimestep/Xlength) &
				                             + (NumFlux(I)*MixTimestep/Height/1.0e6)) &
							/(1.0 + MixTimestep*(TerminalVelocity(Current)/(Height*100) + Windspeed/Xlength))
				END IF
				!WRITE(*,*) I, InitNumConc(I), current%numberofparticles, TerminalVelocity(Current)
				
				!AqMassFlux in units of mol/m2/s
				IF(current%Numberofparticles .GT. 0.) THEN
					DO J = 1, NumAqChems
						InitialMass = InitNumConc(I)*current%AqChems(J)
						EnvMass = EnvNumConc(I)*Env%AqChems(J)
						IF(t .LT. spinup) THEN
							NewMass = (InitialMass + Windspeed*EnvMass*MixTimestep/Xlength)/ &
										  (1+ MixTimestep*(TerminalVelocity(Current)/(Height*100) + Windspeed/Xlength))
						ELSE
							NewMass = (InitialMass + Windspeed*EnvMass*MixTimestep/Xlength &
						            + AqMassFlux(I,J)*MixTimestep/Height/1.0e6)/ &
										  (1+ MixTimestep*(TerminalVelocity(Current)/(Height*100) + Windspeed/Xlength))
						END IF
						Current%AqChems(J) = NewMass/current%numberofparticles
						!IF(current%Numberofparticles .LT. 1.0) WRITE(*,*) I,J, Current%AqChems(J)
						!IF(IsNAN(Current%AqChems(J))) WRITE(*,*) "Error in Particle Diluton! ", I, J
					END DO

				!OrgMassFlux in units of mol/m2/s
					DO J = 1, NumOrgChems
						InitialMass = InitNumConc(I)*current%OrgChems(J)
						EnvMass = EnvNumConc(I)*Env%OrgChems(J)
						IF (t .LT. spinup) THEN
							NewMass = (InitialMass + Windspeed*EnvMass*MixTimestep/Xlength)/ &
										  (1+ MixTimestep*(TerminalVelocity(Current)/(Height*100) + Windspeed/Xlength))
			
						ELSE
							NewMass = (InitialMass + Windspeed*EnvMass*MixTimestep/Xlength &
						            + OrgMassFlux(I,J)*MixTimestep/Height/1.0e6)/ &
										  (1+ MixTimestep*(TerminalVelocity(Current)/(Height*100) + Windspeed/Xlength))
						END IF
						Current%OrgChems(J) = NewMass/current%numberofparticles
					END DO
				
				!AqOrgMassFlux in units of mol/m2/s
					DO J = 1, NumAqOrgChems
						InitialMass = InitNumConc(I)*current%AqOrgChems(J)
						EnvMass = EnvNumConc(I)*Env%AqOrgChems(J)
						IF (t .LT. spinup) THEN
							NewMass = (InitialMass + Windspeed*EnvMass*MixTimestep/Xlength)/ &
										  (1+ MixTimestep*(TerminalVelocity(Current)/(Height*100) + Windspeed/Xlength))
						ELSE
							NewMass = (InitialMass + Windspeed*EnvMass*MixTimestep/Xlength &
						            + AqOrgMassFlux(I,J)*MixTimestep/Height/1.0e6)/ &
										  (1+ MixTimestep*(TerminalVelocity(Current)/(Height*100) + Windspeed/Xlength))
						END IF
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
			!WRITE(*,*) "AOD: ", AerosolOpticalDepth ()	
			!WRITE(*,*) "Dilution Okay"	
			!CALL SpillBeans()
			CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSSA, 0)
			
			!Write Output to file every minute 
			IF (MOD(t,10.0) .LT. 0.01) THEN
				
				!Write current time in minutes
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
					
				!Write gas phase concentrations in ppb
				DO q = 1, HowManyGasChems
		
					Chem(q) = GridGasChem(q) !molecules/cm3

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
				current => particles%first
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
					CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
				END IF
			END IF
			
			!Update time
			t = t + MixTimestep
			ut = ut + MixTimeStep/3600.

			IF (MOD(t,CoagTimeStep)<=0.01) THEN
				WRITE(*,*) trim(real2str(t/60)), " minutes done"
			END IF

		END DO


END SUBROUTINE EulerianBox

SUBROUTINE InitializeGasEulerianBox(NumGasChems, Gas_ER) 

	!! Include modules
	USE ModelParameters
	USE InfrastructuralCode	
	USE GridPointFields
	USE Chemistry
	USE Aerosols		
	USE Condensation
	USE OutputRoutines
	USE Coagulation
	
	implicit none

        REAL*8, INTENT(IN) :: NumGasChems
        REAL*8, INTENT(INOUT) ::  Gas_ER(:)

  

	!Set Inital Gas Concentrations to Background
!	DO q = 1, NumGasChems
!		Chem(q) = GasPhaseBackground(q)
!	END DO

	!Update chemical concentrations
!	CALL UpdateChemicalConcentrations(Chem(1:HowManyEvolveGasChems))

END SUBROUTINE
