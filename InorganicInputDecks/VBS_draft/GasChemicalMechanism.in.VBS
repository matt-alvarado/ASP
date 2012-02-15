!! GasChemicalMechanism.in
!! ASP
!! Matt Alvarado, 2004-2006 (c)
!! Based on MELAM by Donnan Steele, 2000-2004 (c)
!! This is a gas phase mechanism input deck for ASP.
!! Each chemical must be typed the same way each time, and 
!! must match one of the chemicals in GasPhaseChems.in
!! Use + to separate compounds, an unreacting M to indicate pressure 
!! dependence, and => to signify the reactant / product boundary
!! You can use stoichiometric factors (you must put whitespace
!! between the factor and the chemical), but be careful with negatives; 
!! the must follow the form:
!!
!!	RPR5   + OH      => O3 + -1.0 HO2 + UR14 + H2O                
!!
!! 
!! Water is referred to at "H2O" - you may use it even though it
!! is not specified in the chemical input decks.
!!
!! After specifying the reaction, put a semi-colon and specify the
!! appropriate rate coefficients
!!
!!  For 1ST ORDER RXNS, the reaction rate is proportional to concentration.
!!  So:
!!    1. Reaction Rate (sec^-1) (should be less than 1)
!!	
!!	*IF the first term is 0., this is a photolysis reaction, and the
!!	the second term should be the index of the appropriate column in PhotoRates.in
!!	(ignoring the first time column). These should correspond to the reaction indices
!!	of TUV v4.1		
!!
!!	*OR if the first term is 2., that is a flag for an isomerization, and next 4 terms follow
!!	the two-body formula k = A * T^B exp(C/T)^D but rate is in 1/s
!!
!!    *OR if the first term is 3., that is a flag for a heterogeneous reaction, and the following
!!    term is the accomodation coefficient, and the third term is the molecular weight in g/mol
!!    of the reactant, and the fourth is an unused dummy variable (used to be 
!!    the gas-phase diffusivity in m2/s - now this is calculated in the code). 
!!
!!	*OR if the first term is 4., that is a flag for a thermal degradation reaction. This reaction is
!!	assumed to be in equilibrium with a third-order formation reaction with a pseudo-2nd order
!!	reaction rate of 
!!     k_f = (k_0*M* k_inf) / (k_0*M + k_inf) * F^((1+log10(k_0*M/k_inf)^2)^-1) (cm3/molecules/s)
!!  	Each of k_0 and k_inf are temperature dependent, as follows:
!!     k(T) = k(300)*(T/300)^(-n)
!!  	With appropriate definitions of F.  So, specify with appropriate units:
!!	          1. k_0(300)
!!		  2. n for k_0
!!		  3. k_inf(300)
!!		  4. n for k_inf
!!		  5. F_c
!!  	The last two terms are for the equilibrium constant Keq (cm3/molecules)
!!	Keq = A*exp(C/T)
!!	So the rate constant is k = k_f/Keq 
!!   
!!    *OR if the first term is 5., the following term is the fixed rate in 1/sec. 
!!	
!!	*OR if the first term is 6. or 7., that's the flag for the HNO4 and N2O5 thermal
!!	degradation reactions, respectively, hardwired at the CACM value. 
!!
!!	*OR if the first term is 8., that's the flag for the O1D reaction, 
!!	hardwired to the CACM value.
!!	*OR if the first term is 9., that's the flag for the O+O2+M => O3 + M reaction, 
!!	hardwired to the (first order) CACM value.
!!    *OR if the first terms is 10 that's a hardwired RAD* + O2 reaction.
!!	
!!  For Second Order Reactions, assume the form:
!!	k = A * T^B exp(C/T)^D
!!  All units are appropriate combinations of molecules, cm, seconds, and K.
!!  So, specify:  1. "A": Pre-Exponential Factor
!!		  2. "B": Direct Temperature Dependence
!!		  3. "C": Arrhenius Temperature Dependence
!!		  4. "D": Arrhenius Off-Factor for Temperature Depedence 
!!	
!!	*BUT if the reaction can form an organic nitrate, use flag 1. for nitrate formation 
!!	and flag 2. for normal case. Terms 2 through 5 are as for all other 
!!	2-body reactions. The sixth is the number of carbons in the RO2 radical.
!!	The seventh is the scaling factor to be used (if primary or tertiary RO2 radical).
!!    The eigth and ninth terms are scaling factors for the ratio that are used in CACM
!!    for RO2_9, RO2_17, RO2_24 and RO2_27 - otherwise they are 1.0 and 0.0 and have no effect. 
!!	See Carter and Atkinson, J. Atm. Chem 8 (1989) 165-173 for details 
!!	of scaling for these reactions and the CACM rates file for the special scaling factors.
!!
!!	*OR if the first term is 3.0, that is the flag for the CO + OH reaction.
!!	That rate constant is hard-wired as 1.5e-13*(1+0.6*P) where P is in bar.
!!
!!    *OR if the first term is 4., that is the flag for a reaction rate constant of the form
!!     k = k1 + k2*M
!!	Where k1 and k2 follow Arrhenius dependence
!!	k(T) = A*exp(C/T)
!!     So the next four rate terms are:
!!		  1. A for k1 (units cm3/molecule/s_
!!		  2. C for k1 (units K)
!!		  3. A for k2 (units cm6/molecules2/s)
!!		  4. C for k2 (units K)
!!
!!    *OR if the first term is 5., that is the flag for a reaction rate constant of the form
!!     k = k1 + k3*M(1 + k3*M/k2)
!!	Where k1 and k2 follow Arrhenius dependence
!!	k(T) = A*exp(C/T)
!!     So the next four rate terms are:
!!		  1. A for k1 (units cm3/molecule/s)
!!		  2. C for k1 (units K)
!!		  3. A for k2 (units cm3/molecules/s)
!!		  4. C for k2 (units K)
!!		  3. A for k3 (units cm6/molecules2/s)
!!		  4. C for k3 (units K)
!!     *OR if the first term is 6., that is the flag for a 2nd order heterogeneous rxn.
!!         The next term is the reaction rate in cm3*cm/molecule/s, and must be multiplied by
!!		the particle surface area concentration (in cm2/cm3).
!!
!!  For Third order reactions, assume the form:
!!     k = (k_0*M* k_inf) / (k_0*M + k_inf) * F^((1+log10(k_0*M/k_inf)^2)^-1)
!!  Each of k_0 and k_inf are temperature dependent, as follows:
!!     k(T) = k(300)*(T/300)^(-n)
!!  With appropriate definitions of F.  So, specify with appropriate units:
!!	          1. k_0
!!		  2. n for k_0
!!		  3. k_inf
!!		  4. n for k_inf
!!		  5. F_c
!!  *BUT if the first term is 1, just calculate using the thirds order rate constant 
!!  of the form k(T) = k(300)*(T/300)^(-n) where k has units of cm^6/(molecules)^2/s.
!!  (Usually only for O + O2 + M => O3 + M)
!!  *OR if the first term is 2, that's the hard-wired OH+NO2+M reaction rate from CACM
!!	See GasPhaseChemistry.h for more details on reaction rate calculations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PHOTOLYSIS REACTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    These rates must be adjusted for each different model run. 
!!    Here, photolysis rate indices are from TUV v. 4.6 
!!    There must be 24 hours of rate constants in PhotoRates.in corresponding
!!    to these reactions with increments every 15 minutes.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

NO2              => NO    + O	      				; 0. ; 4. 
NO3		     => NO 	  + O2            		; 0. ; 5.
NO3		     => NO2   + O             			; 0. ; 6.
O3		     => O     + O2           			; 0. ; 2.
O3		     => O1D   + O2           			; 0. ; 1. 
HONO             => 0.9 OH  + 0.9 NO + 0.1 HO2 + 0.1 NO2        ; 0. ; 7.
H2O2		     => OH    + OH            			; 0. ; 3.
HCHO             => HO2   + HO2 + CO  				; 0. ; 8.
HCHO             => H2    + CO	      			        ; 0. ; 8.
!Used Acetaldehyde for photolysis rate (propanal higher - index 20)
ALD2             => CO + HO2 + RO2_5 + RO2_T                    ; 0. ; 10.
!Assumed the same as ALD2
RPR1             => CO + HO2 + RO2_20 + RO2_T                   ; 0. ; 10.
!Assumed the same as ALD2
RPR3             => CO + 2.0 HO2 + UR4                          ; 0. ; 10.
!Assumed the same as ALD2
RPR8             => CO + HO2 + RO2_9 + RO2_T                    ; 0. ; 10.
!Assumed the same as ALD2
RPR8             => HO2 + RO2_57 + RO2_T                        ; 0. ; 10.
!Assumed the same as ALD2
RP10             => UR25                                        ; 0. ; 10.
!Assumed photolysis same as acetone
KETL		     => RO2_5 + RO2_8 + 2.0 RO2_T		; 0. ; 11.
!Assumed photolysis same as acetone
KETH		     => RO2_5 + RO2_8 + 2.0 RO2_T		; 0. ; 11. 
MGLY             => CO + HO2 + RO2_8 + RO2_T                    ; 0. ; 12. 
!Assumed the same as MGLY
RP16             => 2.0 CO + OH + HO2                           ; 0. ; 12. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!HETERO. RXNS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!NO2               => 0.5 HONO + 0.5 HNO3                    ; 3. ; 1.0e-3; 46.0;  1.0e-5
!NO2               => 1.0 HONO 					 ; 3. ; 1.0e-3; 46.0;  1.0e-5
!HO2               => 0.5 H2O2                               ; 3. ; 0.1;    33.0;  1.0e-5
!NO3               => HNO3                                   ; 3. ; 1.0e-3; 62.0;  1.0e-5
!N2O5              => 2.0 HNO3                               ; 3. ; 0.1;    108.0; 1.0e-5    
!HO2                => ACID                                  ; 3. ; 0.1;    33.0; 1.0e-5
!O3                => ACID                                   ; 3. ; 2.0e-4;  48.0; 1.0e-5
!MEOH + NO2        => HO2 + HCHO + HNO3 + NO + -1.0 NO2		; 6. ; 4.33e-12
!MEOH + OH         => HO2 + HCHO + H2O					; 6. ; 3.68e-8
!SO2			=> H2SO4						;3.   ; 2.5e-4; 64.0; 1.0e-5

!Furan Reactions  - stoichiometry from MCM of University of Leeds (http://mcm.leeds.ac.uk/MCM/)
FURAN  + OH        =>  HO2 + RP20			    ; 6.2e-11; 	0.; 	0.; 		1. !2-methylfuran from Bierbach
RP20   + OH        => 0.17 RO2_59 + 0.83 RO2_60 + RO2_T     ; 5.6e-11;  0.;   0.;         1. !Rate for 4-oxo-2-pentenal (from methyl furan_
RP20               => 0.4 UR25 + 0.6 RO2_60 + 0.6 HO2	+ 0.6 RO2_T      ; 0. ; 10. 
RP20 + NO3         => RO2_60 + HNO3 + RO2_T                 ; 5.76e-12 ; 0.; -1862; 1. !Rate from Leeds
RP20 + O3 => 1.0675 GLY + 0.82 HO2 + 0.57 OH + 1.265 CO + 0.0675 H2O2 + 0.125 HCHO + 0.1125 RP16 ; 2.0e-18 ; 0.; 0.; 1. 
RO2_59 + NO        => NO2 + HO2 + 2.0 GLY                  ; 2.46e-12; 0.;   180.0;      1. !Assumes 4 carbons, Jenkins
RO2_59 + HO2       => OH  + HO2 + 2.0 GLY                  ; 3.14e-13; 0.;   800.2;      1. !CACM default rate
RO2_59 + RO2_T     => HO2 + 2.0 GLY + RO2_T                ; 1.0e-15; 0.;   	0.;      	1. !CACM default rate
RO2_60 + NO        => NO2 + HO2 + 0.4 GLY + 0.4 CO + 0.6 UR24 ; 2.46e-12; 0.;   180.0;      1. !Assumes 4 carbons, Jenkins
RO2_60 + HO2       => 0.29 O3 + 0.29 ALD2 + 0.71 OH + 0.284 GLY + 0.284 HO2 + 0.284 CO + 0.426 HO2 + 0.426 UR24; 3.14e-13; 0.;   800.2;      1. !CACM default rate
RO2_60 + RO2_T     => HO2 + 0.4 GLY + 0.4 CO + 0.6 UR24  + RO2_T ; 1.0e-15; 0.;   	0.;      	1. !CACM default rate
RO2_60 + NO2 + M   => PN11 + M					; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5;0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PN11  		 => NO2 + RO2_60 + RO2_T	            ; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
PN11 + OH          => GLY + CO + CO + NO2                  ; 3.7e-11; 0. ; 0.; 1. !Rate from Leeds


!!!!!!!!!!!!!!Other Reactions - These are adjusted automatically!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Non-typical reactions
O => O3 ; 9. !O + O2 + M expressed as a first order reaction

!Reaction Rate caluclated using JPL, not CACM
CO    + OH       => HO2   + CO2           			; 3.
!See Table 5 of CACM (Rxn 23)
HNO3  + OH       => NO3   + H2O           			; 5.; 7.2e-15; 785.1; 4.1e-16; 1439.4; 1.9e-33; 724.7 
!See Table 5 of CACM (Rxn 31)
HO2   + HO2      => H2O2  + O2	      			; 4.; 2.2e-13; 598.9; 1.85e-33; 981.4
!Using equilibrium const. from JPL (K = 2.7E-27*EXP(11000/T))
N2O5             => NO2   + NO3					;7. !; 4. ;2.8e-30; 3.5; 2.0e-12; 0.2; 0.45; 2.7e-27; 11000.0 
!Using equilibrium const. from JPL (K = 2.1E-27*EXP(10900/T))
HNO4             => NO2   + HO2					;6. !; 4.; 1.8e-31; 3.2; 4.7e-12; 0.0; 0.60; 2.1e-27; 10900.0

O	+ NO2      => NO    + O2					; 6.5e-12; 	0.; 	119.8;	1.
O     + NO2 + M  => NO3   + M	      				; 9.0e-32; 	2.; 	2.2e-11; 	0. ; 0.8
O3    + NO       => NO2   + O2	      			; 1.8e-12; 	0.;	-1368.9; 	1.
NO2   + O3	     => NO3   + O2					; 1.4e-13;  0.;   -2471.1;	1.
NO    + NO3      => 2.0 NO2 						; 1.8e-11;  0.;   110.7;      1.
NO    + NO       => NO2   + NO2					; 5.09e-18; -1.0;	528.4;	1.

NO2   + NO3 + M  => N2O5  + M						; 2.8e-30;	3.5;	2.0e-12;	0.2;	0.45

NO2   + NO3      => NO    + NO2 + O2      			; 4.5e-14;  0.;	-1258.2;	1.
O1D   + H2O      => OH    + OH	      			; 2.2e-10; 	0.;	0.;		1.

O1D               => O							; 8.

NO    + OH + M   => HONO  + M             			; 7.0e-31;	2.6;	3.6e-11;	0.1;	0.6

OH    + NO2 + M  => HNO3  + M	      				; 2. 

O3    + OH       => HO2   + O2            			; 1.9e-12;	0.;	-1001.5;	1.
HO2   + NO       => NO2   + OH	      			; 3.41e-12; 0.;	271.8;	1.
NO2   + HO2 + M  => HNO4  + M             			; 1.8e-31;	3.2;	4.7e-12;	0.;	0.6

HNO4  + OH       => NO2   + O2 + H2O      			; 1.5e-12;	0.;	362.4;	1.
HO2   + O3       => OH    + O2  + O2  				; 1.4e-14; 	0.;	-598.9;	1.

NO3   + HO2      => 0.8 NO2   + 0.2 HNO3 + 0.8 OH  + O2	; 4.0e-12;	0.;	0.;		1.
O     + O3       => O2    + O2            			; 8.0e-12;	0.;	-2058.4;	1.
OH    + SO2 + M  => H2SO4 + HO2 + M    				; 4.1e-31 ; 3.3; 2.0e-12 ; 0.; 0.45

H2O2  + OH       => HO2   + H2O					; 2.91e-12; 0.;	-161.;	1.
O     + NO       => NO2                  				; 6.75e-6; -2.6;	0.;		1.
HONO  + OH       => NO2   + H2O           			; 2.7e-12;	0.;	261.7;	1.
NO3   + OH       => NO2   + HO2           			; 2.0e-11;	0.;	0.;		1.
NO3   + NO3      => NO2   + NO2 + O2      			; 8.5e-13;	0.;	-2450.9;	1.
OH    + HO2      => H2O   + O2            			; 4.8e-11;	0.;	251.6;	1.

CH4   + OH       => RO2_1 + RO2_T + H2O				; 2.66e-12;	0.;	-1800.2;	1.

HCHO  + OH       => HO2   + CO  + H2O 				; 1.2e-14; 	1.;	286.9;	1.
HCHO  + NO3      => HNO3  + CO  + HO2				; 2.0e-12;	0.;	-2430.8;	1.

MEOH  + OH       => HO2   + HCHO + H2O    			; 6.0e-18;	2.0;	170.1;	1.

ETHE  + OH       => RO2_2 + RO2_T					; 1.96e-12; 0.; 	437.8;	1.
ETHE  + NO3      => RO2_3 + RO2_T					; 4.89e-18; 2.;	-2282.3;	1.
ETHE  + O3	     => 0.63 CO + 0.12 HO2 + 0.12 OH + 0.37 ACID + 1.0 HCHO + 0.14 H2O ; 9.14e-15; 0.; -2580.3; 1.
ETHE  + O        => 0.6 CO + HO2 + 0.6 RO2_1 + 0.4 RO2_4 + RO2_T ; 7.3e-13; 0.; 0.; 1.

!The product ratios are weak functions of temperature - the values at 300 K are used
ETOH  + OH       => 0.95 HO2 + 0.95 ALD2 + 0.05 RO2_2 + 0.05 RO2_T + H2O; 6.18e-18; 2.; 532.; 1. 

OLEL  + OH       => RO2_2 + RO2_T					; 5.86e-12; 0.; 	500.3; 	1.
OLEL  + NO3      => RO2_3 + RO2_T                           ; 1.0e-13;	0.;	-800.2;	1.
OLEL + O3 => 0.56 CO + 0.2 CO2 + 0.36 OH + 0.28 HO2 + 0.5 HCHO + 0.5 ALD2 + 0.24 ACID + 0.1 ALKL + 0.28 RO2_5 + 0.28 RO2_T ; 1.0e-17; 0.; 0.; 1.	
OLEL + O  => 0.5 ALKL + 0.4 ALD2 + 0.1 RO2_4 + 0.1 RO2_5 + 0.2 RO2_T ; 4.66e-12; 0.; 0.; 1.

ALKL  + OH       => RO2_5 + RO2_T + H2O    			; 3.91e-12; 0.; 	0.; 		1.

ALD2  + OH       => RO2_6 + RO2_T + H2O				; 6.91e-12; 0.; 	250.;		1.
ALD2  + NO3	     => HNO3 + RO2_6 + RO2_T				; 3.0e-13;	0.;	-1427.;	1. 

KETL  + OH       => RO2_7 + RO2_T + H2O                     ; 4.91e-12; 0.; 	0.; 		1. 

ISOP + OH 	     => 0.66 RO2_9 + 0.34 RO2_10 + RO2_T		; 2.55e-11;	0.;	410.2;	1.
ISOP + NO3	     => 0.66 RO2_11 + 0.34 RO2_12 + RO2_T		; 3.02e-12;	0.;	-445.9;	1.

!I had to split this into 2 reactions, but stoichiometry still works
ISOP + O3 => 0.068 CO2 + 0.461 CO + 0.5 HCHO + 0.664 OH + 0.366 HO2 + 0.054 OLEL + 0.121 ACID + 0.389 MVK + 0.17 MCR ; 7.86e-15;	0.;	-1912.9;	1.
ISOP + O3 => ISOP + O3 + 0.271 RO2_13 + 0.095 RO2_14 + 0.366 RO2_T ; 7.86e-15;	0.;	-1912.9;	1.

ISOP + O         => 0.925 OLEL + 0.075 ALD2			; 3.5e-11;	0.;	0.;		1.

!From structure-activity relationship (Kwok and Atkinson, 1995)
ALCH + OH        => RO2_2 + RO2_T + H2O                     ; 6.87e-18; 2.; 	909.16;	1.

!From structure-activity relationship (Kwok and Atkinson, 1995)
KETH + OH	     => RO2_16 + RO2_T + H2O				; 2.03e-17; 2.;	311.1;	1.

AROO  + NO3      => HNO3 + RAD1					; 3.77e-12; 0.;	0.;		1.

!From structure-activity relationship (Kwok and Atkinson, 1995)
AROO  + OH       => 0.16 HO2 + 0.16 AROO + 0.1 RO2_17 + 0.1 RO2_T + 0.74 RAD2 + 0.1 H2O ; 2.264e-10; 0.; 0.; 1.

!From structure-activity relationship (Kwok and Atkinson, 1995)
OLEH  + OH       => RO2_18 + RO2_T                          ; 1.77e-17; 2.;	924.46;	1. 
OLEH  + NO3      => RO2_19 + RO2_T                          ; 3.02e-19; 2.;   -376.04;    1.
!I had to split this into 2 reactions, but stoichiometry still works
OLEH  + O3       => 0.56 CO + 0.2 CO2 + 0.36 OH + 0.28 HO2 + 0.5 HCHO + 0.5 RPR1 + 0.12 ACID ; 3.02e-23; 2.; 424.16; 1.
OLEH  + O3       => OLEH  + O3 + 0.12 UR1 + 0.1 ALKM + 0.28 RO2_20 + 0.28 RO2_T ; 3.02e-23; 2.; 424.16; 1.
OLEH  + O        => 0.5 ALKM + 0.4 RPR1 + 0.1 RO2_4 + 0.1 RO2_20 + 0.2 RO2_T ; 1.41e-17; 2.; 424.16; 1.

!From structure-activity relationship (Kwok and Atkinson, 1995)
ALKM  + OH       => RO2_20 + RO2_T + H2O                    ; 1.46e-17; 2.;	616.95;	1.

AROL  + OH       => 0.16 HO2 + 0.16 AROO + 0.1 RO2_21 + 0.74 RAD3 + 0.1 RO2_T + 0.1 H2O ; 3.27e-11; 0.; 0.; 1.

!From structure-activity relationship (Kwok and Atkinson, 1995)
AROH  + OH       => 0.16 HO2 + 0.16 AROO + 0.84 RAD4        ; 3.36e-17; 2.; 	484.94; 	1.

ARAL  + NO3      => HNO3 + O3 -1.0 HO2 + ARAC               ; 1.4e-12; 0.; 	-1872.2;	1.

!Stoichiometry calculated at 298K from expressions in Rob Griffin's CACM Rates file 
ARAL  + OH       => -0.58 HO2 + 0.16 RPR2 + 0.74 O3 + 0.74 ARAC + 0.00454 RO2_22 + 0.00454 RO2_T + 0.09546 RAD5 + 0.74454 H2O ; 1.29e-11; 0.; 0.; 1.  

!From structure-activity relationship (Kwok and Atkinson, 1995)
ARAC   + OH      => 0.16 HO2 + 0.16 UR2 + 0.1 RO2_23 + 0.74 RAD6 + 0.1 RO2_T + 0.1 H2O ;1.055e-12 ; 0.; 0.; 1.

BIOL   + OH      => RO2_24 + RO2_T                          ; 1.7e-10; 0.; 	0.; 		1.
BIOL   + NO3     => RO2_25 + RO2_T                          ; 1.46e-11; 0.;	0.;		1.
BIOL   + O3      => 0.445 CO + 0.055 H2O2 + 0.445 HO2 + 0.89 OH + 0.055 UR3 + 0.445 UR4 + 0.055 RPR3 + 0.445 RO2_26 + 0.445 RO2_T; 2.5e-16; 0.; 0.; 1.
BIOL   + O       => 0.75 UR5 + 0.25 UR6                     ; 1.35e-10; 0.; 	-500.3; 	1.

BIOH   + OH      => RO2_27 + RO2_T                          ; 1.77e-10; 0.; 	0.; 		1.
BIOH   + NO3     => RO2_28 + RO2_T                          ; 2.91e-11; 0.;	0.;		1.
BIOH   + O3      => 0.445 CO + 0.055 H2O2 + 0.89 OH + 0.055 UR7 + 0.055 UR8 + 0.445 RO2_29 + 0.445 RO2_30 + 0.89 RO2_T ; 1.4e-16; 0.; 0.; 1.
BIOH   + O       => 0.75 UR9 + 0.25 UR10                     ; 8.59e-11 ; 0.; 0.; 	1.

PAH    + OH      => 0.16 HO2 + 0.16 UR11 + 0.1 RO2_31 + 0.74 RAD7 + 0.1 RO2_T + 0.1 H2O ; 7.7e-11 ; 0.; 0.; 1.

!From structure-activity relationship (Kwok and Atkinson, 1995)
ALKH   + OH      => RO2_32 + RO2_T + H2O                    ; 6.57e-17; 2.; 361.; 1.
 
RO2_T + HO2      => HO2							; 3.41e-13; 0.;	800.2;	1.
RO2_T + NO       => NO							; 4.2e-12; 0.;	181.2;	1.
RO2_T + RO2_T    => RO2_T						; 1.0e-15; 0.;	0.;		1.

RAD2  	     => RO2_33 + RO2_T					;10.
RAD3  	     => RO2_34 + RO2_T					;10.
RAD4  	     => RO2_35 + RO2_T					;10.
RAD5  	     => RO2_36 + RO2_T					;10.
RAD6  	     => RO2_37 + RO2_T					;10.
RAD7  	     => RO2_38 + RO2_T					;10.

RAD1  + NO2	     => RPR4						; 5.1e-11;	0.;	0.;		1.
RAD2  + NO2	     => RPR4 + H2O					; 5.1e-11;	0.;	0.;		1.
RAD3  + NO2	     => UR12 + H2O					; 5.1e-11;	0.;	0.;		1.
RAD4  + NO2	     => UR13 + H2O					; 5.1e-11;	0.;	0.;		1.
RAD5  + NO2	     => RPR5 + H2O					; 5.1e-11;	0.;	0.;		1.
RAD6  + NO2	     => UR14 + H2O					; 5.1e-11;	0.;	0.;		1.
RAD7  + NO2	     => UR15 + H2O					; 5.1e-11;	0.;	0.;		1.
RO2_1  + NO      => NO2   + HO2  + HCHO   			; 4.09e-12;	0.;	180.2;	1.
RO2_1  + RO2_T   => HO2 + HCHO + RO2_T + O2 			; 1.0e-15; 0.; 	0.; 		1.
RO2_1	 + HO2     => HO2   + OH   + HCHO   			; 3.41e-13;	0.;	800.2;	1.
RO2_2  + NO      => NO2 + HO2 + HCHO + ALD2 			; 2.455e-12; 0.;	180.2;	1.
RO2_2  + RO2_T   => HO2 + HCHO + ALD2 + RO2_T + O2		; 1.0e-15; 0.; 0.; 	1.
RO2_2  + HO2     => OH + HO2 + HCHO + ALD2 			; 3.41e-13; 0.;	800.2;	1.
RO2_3  + NO      => 2.0 NO2 + HCHO + ALD2 			; 2.45e-12; 0.;	180.2;	1.
RO2_3  + RO2_T   => NO2 + HO2 + HCHO + ALD2 + O2 + RO2_T 	; 1.0e-15; 0.;	800.2;	1.
RO2_3  + HO2     => NO2 + HO2 + OH + HCHO + ALD2		; 3.41e-13;	0.; 800.2; 1.
RO2_4  + NO      => NO2 + CO + HO2 + HCHO 			; 3.45e-12;	0.;	180.2;	1.
RO2_4  + RO2_T   => CO + HO2 + HCHO + RO2_T + O2 		;  1.0e-15; 0.; 	0.; 		1.
RO2_4  + HO2     => CO + HO2 + OH + HCHO  			;  3.41e-13; 0.;	800.2;	1.
RO2_5  + NO      => ALKL 						; 1.; 2.91e-12; 0.;	180.2;	1.; 3.; 0.4 ; 1.0; 0.0 !Primary!; 1.; 2.91e-12; 0.;	180.2;	1.; 3.; 0.4 ; 1.0; 0.0 !Primary
RO2_5  + NO      => NO2 + HO2 + ALD2 				; 2.; 2.91e-12; 0.;	180.2;	1.; 3.; 0.4 ; 1.0; 0.0 !Primary
RO2_5  + RO2_T   => HO2 + ALD2 + RO2_T 				; 1.0e-15; 0.;	0.;		1.
RO2_5  + HO2     => HO2 + OH + ALD2					; 3.41e-13; 0.;	800.2;	1.
RO2_6  + NO      => NO2 + CO2 + RO2_5 + RO2_T               ; 1.11e-11; 0.; 	180.2;	1.
RO2_6  + NO2 + M => PAN1 + M						; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN1             => RO2_6 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_6  + HO2     => O3 + ACID 					; 9.743e-14; 0.;	800.2;	1.
RO2_6  + RO2_T   => CO2 + RO2_5 + 2.0 RO2_T + O2		; 1.0e-15; 0.;	0.;		1.
RO2_7  + NO      => NO2 + ALD2 + RO2_8 + RO2_T			; 2.45e-12; 0.;	180.2;	1.
RO2_7  + RO2_T   => ALD2 + RO2_8 + 2.0 RO2_T + O2           ; 1.0e-15; 0.;	0.;		1.
RO2_7  + HO2     => OH + ALD2 + RO2_8 + RO2_T			; 3.41e-13; 0.;	800.2;	1.
RO2_8  + NO      => NO2 + CO2 + RO2_1 + RO2_T			; 1.11e-11; 0.;	180.2;	1.
RO2_8  + NO2 + M => PAN2 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5 ; 0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN2             => RO2_8 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_8  + HO2     => O3 + ACID						; 9.743e-14; 0.;	800.2;	1.
RO2_8  + RO2_T   => CO2 + RO2_1 + 2.0 RO2_T + O2		; 1.0e-15;  0.;	0.;		1. 
RO2_9  + NO      => OLEL 						; 1.; 2.08e-12; 0.;	180.2;	1.; 5.; 0.3 ; 0.5; 0.5 !Tertiary
RO2_9  + NO      => NO2 + HO2 + HCHO + MVK			; 2.; 2.08e-12; 0.;	180.2;	1.; 5.; 0.3 ; 0.5; 0.5 !Tertiary
RO2_9  + RO2_T   => HO2 + MVK + HCHO + RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_9  + HO2     => HO2 + OH + MVK + HCHO                   ; 3.41e-13; 0.;	800.2;	1.
RO2_10 + NO      => NO2 + HO2 + HCHO + MCR                  ; 2.08e-12; 0.;	180.2;	1.
RO2_10 + RO2_T   => HO2 + HCHO + MCR + RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_10 + HO2     => HO2 + OH + HCHO + MCR                   ; 3.41e-13; 0.;	800.2;	1.
RO2_11 + NO      => 2.0 NO2 + HCHO + MVK				; 2.08e-12; 0.;	180.2;	1.
RO2_11 + RO2_T   => NO2 + HCHO + MVK + RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_11 + HO2     => NO2 + OH + HCHO + MVK                   ; 3.41e-13; 0.;	800.2;	1.
RO2_12 + NO      => 2.0 NO2 + HCHO + MCR				; 2.08e-12; 0.;	180.2;	1.
RO2_12 + RO2_T   => NO2 + HCHO + MCR + RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_12 + HO2     => NO2 + OH + HCHO + MCR				; 3.41e-13; 0.;	800.2;	1.
RO2_13 + NO      => NO2 + HCHO + RO2_39 + RO2_T             ; 2.45e-12; 0.;	180.2;	1.
RO2_13 + RO2_T   => HCHO + RO2_39 + 2.0 RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_13 + HO2     => HCHO + OH + RO2_39 + RO2_T              ; 3.41e-13; 0.;	800.2;	1.
RO2_39 + NO      => NO2 + CO2 + RO2_14 + RO2_T              ; 1.11e-11; 0.;	180.2;	1.
RO2_39 + NO2 + M => PAN3 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5; 0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN3             => RO2_39 + NO2 + RO2_T 	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_39 + HO2     => O3 + 0.5 OLEL + 0.5 ACID                ; 9.743e-14; 0.;	800.2;	1.
RO2_39 + RO2_T   => CO2 + RO2_14 + 2.0 RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_14 + NO      => OLEL 						; 1.; 2.91e-12; 0.;	180.2;	1.; 3.; 0.3 ; 1.0; 0.0 !Tertiary
RO2_14 + NO      => NO2 + RO2_7 + RO2_T				; 2.; 2.91e-12; 0.;	180.2;	1.; 3.; 0.3 ; 1.0; 0.0 !Tertiary
RO2_14 + RO2_T   => RO2_7 + 2.0 RO2_T + O2                  ; 1.0e-15;  0.;	0.;		1. 
RO2_14 + HO2     => OH + RO2_7 + RO2_T				; 3.41e-13; 0.;	800.2;	1.

RO2_16 + NO      => NO2 + ALD2 + RO2_8 + RO2_T              ; 1.48e-12; 0.; 	180.2; 	1.
RO2_16 + RO2_T   => ALD2 + RO2_8 + 2.0 RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_16 + HO2     => OH + ALD2 + RO2_8 + RO2_T               ; 3.41e-13; 0.;	800.2;	1.
RO2_17 + NO      => AP1 						; 1.; 1.248e-12; 0.;	180.2;	1.; 8.; 0.4 ; 0.4 ; 0.6 !Primary
RO2_17 + NO      => NO2 + HO2 + RPR2				; 2.; 1.248e-12; 0.;	180.2;      1.; 8.; 0.4 ; 0.4 ; 0.6 !Primary
RO2_17 + RO2_T   => HO2 + RPR2 + RO2_T + O2                 ;  1.0e-15;  0.;	0.;		1. 
RO2_17 + HO2     => HO2 + OH + RPR2                         ; 3.41e-13; 0.;	800.2;	1.
RO2_18 + NO      => AP2 						; 1.; 1.25e-12; 0.;	180.2;	1.; 8.; 1.0 ; 1.0 ; 0.0 !Secondary
RO2_18 + NO      => NO2 + HO2 + UR16				; 2.; 1.25e-12; 0.;	180.2;      1.; 8.; 1.0 ; 1.0 ; 0.0 !Secondary
RO2_18 + RO2_T   => HO2 + UR16 + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_18 + HO2     => HO2 + OH + UR16                         ; 3.41e-13; 0.;	800.2;	1.
RO2_19 + NO      => 2.0 NO2 + HCHO + RPR1				; 1.05e-12; 0.;	180.2;      1.
RO2_19 + RO2_T   => NO2 + HCHO + RPR1 + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_19 + HO2     => NO2 + OH + HCHO + RPR1                  ; 3.41e-13; 0.;	800.2;	1.
RO2_20 + NO      => AP3 						; 1.; 1.25e-12; 0.;	180.2;	1.; 8.; 1.0 ; 1.0 ; 0.0 !Secondary
RO2_20 + NO      => NO2 + RO2_18 + RO2_T				; 2.; 1.25e-12; 0.;	180.2;      1.; 8.; 1.0 ; 1.0 ; 0.0 !Secondary
RO2_20 + RO2_T   => RO2_18 + 2.0 RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_20 + HO2     => OH + RO2_18 + RO2_T                     ; 3.41e-13; 0.;	800.2;	1.
RO2_21 + NO      => AP4 						; 1.; 1.05e-12; 0.;	180.2;	1.; 9.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_21 + NO      => NO2 + HO2 + ARAL				; 2.; 1.05e-12; 0.;	180.2;      1.; 9.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_21 + RO2_T   => HO2 + ARAL + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_21 + HO2     => HO2 + OH + ARAL					; 3.41e-13; 0.;	800.2;	1.
RO2_22 + NO      => AP5 						; 1.; 1.25e-12; 0.;	180.2;	1.; 8.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_22 + NO      => NO2 + HO2 + RPR6				; 2.; 1.25e-12; 0.;	180.2;      1.; 8.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_22 + RO2_T   => HO2 + RPR6 + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_22 + HO2     => HO2 + OH + RPR6 	                  ; 3.41e-13; 0.;	800.2;	1.
RO2_23 + NO      => AP6 						; 1.; 1.25e-12; 0.;	180.2;	1.; 8.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_23 + NO      => NO2 + HO2 + RPR7				; 2.; 1.25e-12; 0.;	180.2;      1.; 8.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_23 + RO2_T   => HO2 + RPR7 + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_23 + HO2     => HO2 + OH + RPR7 				; 3.41e-13; 0.;	800.2;	1.
RO2_24 + NO      => AP7 						; 1.; 8.89e-13; 0.;	180.2;	1.; 10.; 0.3 ; 0.5 ; 0.5 !Tertiary
RO2_24 + NO      => NO2 + HO2 + RPR3				; 2.; 8.89e-13; 0.;	180.2;      1.; 10.; 0.3 ; 0.5 ; 0.5 !Tertiary
RO2_24 + RO2_T   => HO2 + RPR3 + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_24 + HO2     => HO2 + OH + RPR3 				; 3.41e-13; 0.;	800.2;	1.
RO2_25 + NO      => 2.0 NO2 + RPR3                          ; 8.89e-13; 0.; 	180.2;	1.
RO2_25 + RO2_T   => NO2 + RPR3 + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_25 + HO2     => NO2 + OH + RPR3 				; 3.41e-13; 0.;	800.2;	1.
RO2_26 + NO      => NO2 + UR17 + RO2_8 + RO2_T			; 8.89e-13; 0.; 	180.2;	1.
RO2_26 + RO2_T   => UR17 + RO2_8 + 2.0 RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_26 + HO2     => UR17 + OH + RO2_8 + RO2_T			; 3.41e-13; 0.;	800.2;	1.
RO2_27 + NO      => AP8 						; 1.; 8.89e-13; 0.;	180.2;	1.; 10.; 0.3 ; 0.5; 0.5 !Tertiary
RO2_27 + NO      => NO2 + HO2 + 0.3 UR7 + 0.7 UR10		; 2.; 8.89e-13; 0.;	180.2;      1.; 10.; 0.3 ; 0.5; 0.5 !Tertiary
RO2_27 + RO2_T   => HO2 + 0.3 UR7 + 0.7 UR10 + RO2_T + O2   ; 1.0e-15;  0.;	0.;		1. 
RO2_27 + HO2     => HO2 + OH + 0.3 UR7 + 0.7 UR10	 	; 3.41e-13; 0.;	800.2;	1.
RO2_28 + NO      => 2.0 NO2 + 0.3 UR7 + 0.7 UR10            ; 8.89e-13; 0.; 	180.2;	1.
RO2_28 + RO2_T   => NO2 + 0.3 UR7 + 0.7 UR10 + RO2_T + O2   ; 1.0e-15;  0.;	0.;		1. 
RO2_28 + HO2     => NO2 + OH + 0.3 UR7 + 0.7 UR10 	      ; 3.41e-13; 0.;	800.2;	1.
RO2_29 + NO      => AP9 						; 1.; 1.05e-12; 0.;	180.2;	1.; 9.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_29 + NO      => NO2 + RO2_40 + RO2_T				; 2.; 1.05e-12; 0.;	180.2;      1.; 9.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_29 + RO2_T   => RO2_40 + 2.0 RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_29 + HO2     => OH + RO2_40 + RO2_T				; 3.41e-13; 0.;	800.2;	1.
RO2_40 + NO      => NO2 + RPR8 + RO2_8 + RO2_T			; 1.05e-12; 0.;	180.2;      1.
RO2_40 + RO2_T   => RPR8 + RO2_8 + 2.0 RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_40 + HO2     => OH + RPR8 + RO2_8 + RO2_T			; 3.41e-13; 0.;	800.2;	1.
RO2_30 + NO      => NO2 + UR18 + RO2_8 + RO2_T              ; 8.89e-13; 0.; 	180.2;	1.
RO2_30 + RO2_T   => UR18 + RO2_8 + 2.0 RO2_T + O2           ; 1.0e-15;  0.;	0.;		1. 
RO2_30 + HO2     => OH + UR18 + RO2_8 + RO2_T              	; 3.41e-13; 0.;	800.2;	1.
RO2_31 + NO      => AP10 						; 1.; 6.32e-13; 0.;	180.2;	1.; 12.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_31 + NO      => NO2 + HO2 + UR19				; 2.; 6.32e-13; 0.;	180.2;      1.; 12.; 0.4 ; 1.0 ; 0.0 !Primary
RO2_31 + RO2_T   => HO2 + UR19 + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_31 + HO2     => HO2 + OH + UR19					; 3.41e-13; 0.;	800.2;	1.
RO2_32 + NO      => AP11 						; 1.; 3.2e-13; 0.;	180.2;	1.; 16.; 1.0 ; 1.0 ; 0.0 !Secondary
RO2_32 + NO      => NO2 + RO2_41 + RO2_T              	; 2.; 3.2e-13; 0.;	180.2;      1.; 16.; 1.0 ; 1.0 ; 0.0 !Secondary
RO2_32 + RO2_T   => RO2_41 + 2.0 RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_32 + HO2     => OH + RO2_41 + RO2_T        	      	; 3.41e-13; 0.;	800.2;	1.
RO2_41 + NO      => AP12 						; 1.; 3.2e-13; 0.;	180.2;	1.; 16.; 1.0 ; 1.0 ; 0.0 !Secondary
RO2_41 + NO      => NO2 + HO2 + UR20				; 2.; 3.2e-13; 0.;	180.2;      1.; 16.; 1.0 ; 1.0 ; 0.0 !Secondary
RO2_41 + RO2_T   => HO2 + UR20 + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_41 + HO2     => HO2 + OH + UR20			      	; 3.41e-13; 0.;	800.2;	1.
RO2_33	     => RO2_42 + RO2_T                          ; 2.; 1.97e42; -11.4; -9460.0; 1.  
RO2_33 + NO      => NO2 + 0.3 HO2 + 0.3 RPR9 + 0.7 UR35 + 0.7 RO2_1 + 0.7 RO2_T		; 1.25e-12; 0.;	180.2;      1.
RO2_33 + RO2_T   => 0.3 HO2 + 0.3 RPR9 + 0.7 UR35 + 1.7 RO2_T + 0.7 RO2_1 + O2  ; 1.0e-15;  0.;	0.;		1. 
RO2_33 + HO2     => 0.3 HO2 + 0.3 OH + 0.3 RPR9 + 0.7 UR35 + 0.7 RO2_1 + 0.7 RO2_T          ; 3.41e-13; 0.;	800.2;	1.
RO2_42 + NO      => NO2 + HO2 + RP10 + MGLY			; 1.25e-12; 0.;	180.2;      1.
RO2_42 + RO2_T   => HO2 + RP10 + MGLY + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_42 + HO2     => HO2 + OH + RP10 + MGLY                  ; 3.41e-13; 0.;	800.2;	1.
RO2_34	     => RO2_43 + RO2_T                          ; 2.; 1.97e42; -11.4; -9460.0; 1.  
RO2_34 + NO      => NO2 + 0.3 HO2 + 0.3 RP11 + 0.7 UR35 + 0.7 RO2_1 + 0.7 RO2_T		; 1.05e-12; 0.;	180.2;      1.
RO2_34 + RO2_T   => 0.3 HO2 + 0.3 RP11 + 0.7 UR35 + 0.7 RO2_1 + 1.7 RO2_T + O2  ; 1.0e-15;  0.;	0.;		1. 
RO2_34 + HO2     => 0.3 HO2 + 0.3 OH + 0.3 RP11 + 0.7 UR35 + 0.7 RO2_1 + 0.7 RO2_T          ; 3.41e-13; 0.;	800.2;	1.
RO2_43 + NO      => NO2 + HO2 + RP10 + MGLY			; 1.05e-12; 0.;	180.2;      1.
RO2_43 + RO2_T   => HO2 + RP10 + MGLY + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_43 + HO2     => HO2 + OH + RP10 + MGLY                  ; 3.41e-13; 0.;	800.2;	1.
RO2_35	     => RO2_44 + RO2_T                          ; 2.; 1.804e42; -11.4; -9460.0; 1.  
RO2_35 + NO      => NO2 + 0.5 HO2 + 0.5 RP11 + 0.5 RO2_5 + 0.5 RO2_T + 0.5 UR35		; 8.89e-13; 0.;	180.2;      1.
RO2_35 + RO2_T   => 0.5 HO2 + 0.5 RP11 + 0.5 RO2_5 + 0.5 UR35 + 1.5 RO2_T + O2  ; 1.0e-15;  0.;	0.;		1. 
RO2_35 + HO2     => 0.5 HO2 + 0.5 OH + 0.5 RP11 + 0.5 UR35 + RO2_25 + RO2_T         ; 3.41e-13; 0.;	800.2;	1.
RO2_44 + NO      => NO2 + HO2 + RP10 + MGLY			; 8.89e-13; 0.;	180.2;      1.
RO2_44 + RO2_T   => HO2 + RP10 + MGLY + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_44 + HO2     => HO2 + OH + RP10 + MGLY                  ; 3.41e-13; 0.;	800.2;	1.
RO2_36	     => RO2_45 + RO2_T                          ; 2.; 1.97e42; -11.4; -9460.0; 1.  
RO2_36 + NO      => NO2 + HO2 + 0.3 RP12 + 0.7 UR35		; 1.25e-12; 0.;	180.2;      1.
RO2_36 + RO2_T   => HO2 + 0.3 RP12 + 0.7 UR35 + RO2_T + O2  ; 1.0e-15;  0.;	0.;		1. 
RO2_36 + HO2     => HO2 + OH + 0.3 RP12 + 0.7 UR35          ; 3.41e-13; 0.;	800.2;	1.
RO2_45 + NO      => NO2 + HO2 + RP10 + MGLY			; 1.25e-12; 0.;	180.2;      1.
RO2_45 + RO2_T   => HO2 + RP10 + MGLY + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_45 + HO2     => HO2 + OH + RP10 + MGLY                  ; 3.41e-13; 0.;	800.2;	1.
RO2_37	     => RO2_46 + RO2_T                          ; 2.; 1.97e42; -11.4; -9460.0; 1.  
RO2_37 + NO      => NO2 + HO2 + 0.3 RP13 + 0.7 UR35	      ; 1.25e-12; 0.;	180.2;      1.
RO2_37 + RO2_T   => HO2 + 0.3 RP13 + 0.7 UR35 + RO2_T + O2  ; 1.0e-15;  0.;	0.;		1. 
RO2_37 + HO2     => HO2 + OH + 0.3 RP13 + 0.7 UR35          ; 3.41e-13; 0.;	800.2;	1.
RO2_46 + NO      => NO2 + HO2 + RP10 + MGLY			; 1.25e-12; 0.;	180.2;      1.
RO2_46 + RO2_T   => HO2 + RP10 + MGLY + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_46 + HO2     => HO2 + OH + RP10 + MGLY                  ; 3.41e-13; 0.;	800.2;	1.
RO2_38	     => RO2_47 + RO2_T                          ; 2.; 1.97e42; -11.4; -9460.0; 1.  
RO2_38 + NO      => NO2 + HO2 + RP14				; 6.32e-13; 0.;	180.2;      1.
RO2_38 + RO2_T   => HO2 + RP14 + RO2_T + O2                 ; 1.0e-15;  0.;	0.;		1. 
RO2_38 + HO2     => HO2 + OH + RP14                         ; 3.41e-13; 0.;	800.2;	1.
RO2_47 + NO      => NO2 + HO2 + RP15 + MGLY			; 6.32e-13; 0.;	180.2;      1.
RO2_47 + RO2_T   => HO2 + RP15 + MGLY + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_47 + HO2     => HO2 + OH + RP15 + MGLY                  ; 3.41e-13; 0.;	800.2;	1.

MGLY   + OH      => RO2_48 + RO2_T + H2O                    ; 1.72e-11; 0.;   0.; 		1.
MGLY   + NO3     => HNO3 + RO2_48 + RO2_T                   ; 1.4e-12; 	0.; 	-1897.3;	1.

RO2_48 + NO      => NO2 + CO2 + RO2_8 + RO2_T               ; 1.11e-11; 0.; 	180.2;	1.
RO2_48 + NO2 + M => PAN4 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5;0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN4             => RO2_48 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_48 + HO2     => O3 + UR21 					; 9.743e-14; 0.;	800.2;	1.
RO2_48 + RO2_T   => CO2 + RO2_8 + 2.0 RO2_T + O2		; 1.0e-15; 0.;	0.;		1.

MVK    + OH      => RO2_49 + RO2_T                          ; 4.14e-12;	0.;	452.9;	1.
!Had to split into 2 reactions; but stoichiometry works
MVK    + O3      => 0.56 CO + 0.2 CO2 + 0.28 HO2 + 0.36 OH + 0.5 MGLY + 0.5 HCHO + 0.12 ACID ; 7.5e-16; 0.; -1519.9; 1.
MVK    + O3      => MVK + O3 + 0.1 ALD2 + 0.12 UR21 + 0.28 RO2_8 + 0.28 RO2_T + 0.2 H2O ; 7.5e-16; 0.; -1519.9; 1.
MVK    + O       => 0.85 KETL + 0.15 RO2_4 + 0.15 RO2_8 + 0.3 RO2_T ; 4.32e-12; 0.; 0.; 	1.

RO2_49 + NO      => NO2 + HO2 + MGLY + HCHO			; 2.45e-12; 0.;	180.2;      1.
RO2_49 + RO2_T   => HO2 + MGLY + HCHO + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_49 + HO2     => HO2 + OH + MGLY + HCHO                  ; 3.41e-13; 0.;	800.2;	1.

MCR    + OH      => 0.3 RO2_50 + 0.3 H2O + 0.7 RO2_51 + RO2_T ; 1.86e-11; 0.;	176.1;	1.
MCR    + NO3     => 0.3 HNO3 + 0.3 RO2_50 + 0.7 RO2_52 + RO2_T ; 1.5e-12; 0.; 	-1726.2;	1.
MCR    + O3      => 0.41 CO + 0.41 HO2 + 0.82 OH + 0.5 HCHO + 0.59 MGLY + 0.09 ACID + 0.41 RO2_53 + 0.41 RO2_T ; 1.36e-15; 0.; -2113.7; 1.
MCR    + O       => 0.15 CO + 0.15 HO2 + 0.85 ALD2 + 0.15 RO2_7 + 0.15 RO2_T ;6.34e-12; 0.; 0.; 1.
 
RO2_50 + NO      => NO2 + CO2 + RO2_14 + RO2_T              ; 1.11e-11; 0.; 	180.2;	1.
RO2_50 + NO2 + M => PAN5 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5;0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN5             => RO2_50 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_50 + HO2     => O3 + 0.5 ACID + 0.5 OLEL 		      ; 9.743e-14; 0.;	800.2;	1.
RO2_50 + RO2_T   => CO2 + RO2_14 + 2.0 RO2_T + O2		; 1.0e-15; 0.;	0.;		1.
RO2_51 + NO      => NO2 + HO2 + MGLY + HCHO			; 2.45e-12; 0.;	180.2;      1.
RO2_51 + RO2_T   => HO2 + MGLY + HCHO + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_51 + HO2     => HO2 + OH + MGLY + HCHO                  ; 3.41e-13; 0.;	800.2;	1.
RO2_52 + NO      => 2.0 NO2 + MGLY + HCHO				; 2.45e-12; 0.;	180.2;      1.
RO2_52 + RO2_T   => NO2 + MGLY + HCHO + RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_52 + HO2     => NO2 + OH + MGLY + HCHO                  ; 3.41e-13; 0.;	800.2;	1.
RO2_53 + NO      => NO2 + HCHO + RO2_54 + RO2_T			; 2.91e-12; 0.;	180.2;      1.
RO2_53 + RO2_T   => HCHO + RO2_54 + 2.0 RO2_T + O2          ; 1.0e-15;  0.;	0.;		1. 
RO2_53 + HO2     => OH + HCHO + RO2_54 + RO2_T              ; 3.41e-13; 0.;	800.2;	1.
RO2_54 + NO      => CO2 + CO + NO2 + HO2				; 1.11e-11; 0.; 	180.2;	1.
RO2_54 + NO2 + M => PAN6 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5;0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN6             => RO2_54 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_54 + HO2     => O3 + RP16 					; 9.743e-14; 0.;	800.2;	1.
RO2_54 + RO2_T   => CO2 + CO + HO2 + RO2_T + O2			; 1.0e-15; 0.;	0.;		1.
 
RPR1   + OH      => RO2_55 + RO2_T + H2O                    ; 5.71e-17; 2.;	578.8;	1.
RPR1   + NO3     => HNO3 + RO2_55 + RO2_T                   ; 2.48e-18; 2.;	-1099.;	1.
RO2_55 + NO      => CO2 + NO2 + RO2_20 + RO2_T			; 1.11e-11; 0.; 	180.2;	1.
RO2_55 + NO2 + M => PAN7 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5;0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN7             => RO2_55 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_55 + HO2     => O3 + UR1	 					; 9.743e-14; 0.;	800.2;	1.
RO2_55 + RO2_T   => CO2 + RO2_20 +  2.0 RO2_T + O2		; 1.0e-15; 0.;	0.;		1.
RPR2   + OH      => O3 + -1.0 HO2 + UR2 + H2O                 ; 1.291e-11; 0.; 	0.;		1.
RPR3   + OH      => RO2_56 + RO2_T + H2O                    ; 4.85e-17; 2.;	682.75;	1.
RPR3   + NO3     => HNO3 + RO2_56 + RO2_T 			; 2.105e-18; 2.; 	-994.25;	1.
RO2_56 + NO      => CO2 + NO2 + HO2 + UR4				; 1.11e-11; 0.; 	180.2;	1.
RO2_56 + NO2 + M => PAN8 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5;0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN8             => RO2_56 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_56 + HO2     => O3 + UR3	 					; 9.743e-14; 0.;	800.2;	1.
RO2_56 + RO2_T   => CO2 + HO2 + UR4 +  RO2_T + O2		; 1.0e-15; 	0.;	0.;		1.
RPR4   + NO3     => HNO3 + RAD8                             ; 3.773e-12; 0.;	0.;		1.
RAD8   + NO2     => UR22 + H2O                              ; 2.2955e-11; 	0.;	151.0;	1.
RPR5   + OH      => O3 + -1.0 HO2 + UR14 + H2O                ; 1.291e-11;	0.;	0.;		1.
RPR6   + OH      => O3 + -1.0 HO2 + RPR7 + H2O                ; 1.291e-11;	0.;	0.;		1.
RPR7   + OH      => O3 + -1.0 HO2 + ADAC + H2O                ; 1.291e-11;	0.;	0.;		1.
RPR8   + OH      => RO2_57 + RO2_T + H2O                    ; 1.145e-10; 0.;	0.;		1.
RPR8   + NO3     => HNO3 + RO2_57 + RO2_T 			; 4.97e-12;	0.;	-1677.;	1.
RO2_57 + NO      => CO2 + NO2 + RO2_9 + RO2_T			; 1.11e-11; 0.; 	180.2;	1.
RO2_57 + NO2 + M => PAN9 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5;0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PAN9             => RO2_57 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_57 + HO2     => O3 + UR23	 					; 9.743e-14; 0.;	800.2;	1.
RO2_57 + RO2_T   => CO2 + RO2_9 + 2.0 RO2_T + O2		; 1.0e-15; 	0.;	0.;		1.
RPR9   + OH      => O3 + -1.0 HO2 + RP17 + H2O 			; 2.414e-10; 0.;	0.;		1.
RP10   + OH      => HO2 + UR24 + H2O                        ; 1.098e-10; 0.;	0.;		1.
RP11   + OH      => 0.67 O3 + -0.67 HO2 + 0.67 UR26 + 0.33 UR24 + 0.67 H2O + 0.33 RO2_13 ; 1.964e-10; 0.;	0.;		1.
RP12   + OH      => O3 + -1.0 HO2 + RP13 + H2O			; 2.407e-10; 0.;  0.;         1.
RP13   + OH      => O3 + -1.0 HO2 + RP18 + H2O			; 2.189e-10; 0.;  0.;         1.
RP14   + OH      => O3 + -1.0 HO2 + RP19 + H2O			; 1.452e-10; 0.;  0.;         1.
RP15   + OH      => O3 + -1.0 HO2 + UR27 + H2O			; 1.30e-11;  0.;  0.;         1.
RP16   + OH      => RO2_58 + RO2_T + H2O                    ; 1.72e-11;  0.;	0.;		1.
RP16   + NO3     => HNO3 + RO2_58 + RO2_T                   ; 1.4e-12; 	0.; 	-1897.3;	1.
RO2_58 + NO      => CO + CO2 + NO2 + OH				; 1.11e-11; 0.; 	180.2;	1.
RO2_58 + NO2 + M => PN10 + M						; 9.7e-29 ; 5.6; 9.3e-12 ; 1.5;0.6
!Keq from Rob Griffin, Keq = 8.62e-29*exp(13954.0)
PN10             => RO2_58 + NO2 + RO2_T	      		; 4. ; 9.7e-29 ; 5.6 ; 9.3e-12 ; 1.5 ; 0.6; 8.62e-29; 13954.0
RO2_58 + HO2     => O3 + UR28	 					; 9.743e-14; 0.;	800.2;	1.
RO2_58 + RO2_T   => CO + CO2 + OH +  RO2_T + O2		      ; 1.0e-15; 	0.;	0.;		1.
RP17   + OH      => O3 + -1.0 HO2 + UR29 + H2O			; 2.195e-10; 0.;  0.;         1.
RP18   + OH      => O3 + -1.0 HO2 + UR30 + H2O			; 1.970e-10; 0.;  0.;         1.
RP19   + OH      => O3 + -1.0 HO2 + UR31 + H2O			; 1.336e-10; 0.;  0.;         1.

IVOC2 + OH => IVOC1 ; 4.0e-11; 0.;  0.;         1.
IVOC3 + OH => IVOC2 ; 4.0e-11; 0.;  0.;         1.
IVOC4 + OH => IVOC3 ; 4.0e-11; 0.;  0.;         1. 
IVOC5 + OH => IVOC4 ; 4.0e-11; 0.;  0.;         1.
IVOC6 + OH => IVOC5 ; 4.0e-11; 0.;  0.;         1.
IVOC7 + OH => IVOC6 ; 4.0e-11; 0.;  0.;         1. 
IVOC8 + OH => IVOC7 ; 4.0e-11; 0.;  0.;         1.
IVOC9 + OH => IVOC8 ; 4.0e-11; 0.;  0.;         1.

END-OF-FILE-012343210-ELIF-FO-DNE


