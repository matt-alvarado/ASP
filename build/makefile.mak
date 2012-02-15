##############################################################
# Makefile for MELAM for Donnan's dual xeon desktop system.  #
# (c) Donnan Steele -- 2003				     #
##############################################################

###############
##### Options #

EXENAM = melam   		# Executable Name

SRCDIR = .\SOURCECODE  		# Source Directory
OBJDIR = .\OBJECTS     		# Object Directory
RUNDIR = .             		# Run Directory
DRIVER = Driver			# Name of the Driver to be used

# Give a default Configuration, assuming none specified
!IF "$(CFG)" == ""
CFG    = Standard
!MESSAGE No configuration specified. Defaulting to Standard.
!ENDIF 

# Give a default Configuration, assuming none specified
!IF "$(FORTRANCOMPILER)" == ""
FORTRANCOMPILER    = Intel
!MESSAGE No compiler specified. Defaulting to Intel...
!MESSAGE ...you may also specify 'Compaq' from the command line
!MESSAGE ...by setting the FORTRANCOMPILER flag
!MESSAGE 
!ENDIF 

################################
# Give Comments on the options #
!MESSAGE Source Code Found in Relative Directory Specification $(SRCDIR)
!MESSAGE Objects  Stored in Relative  Directory  Specification $(SRCDIR)
!MESSAGE 

###
### There is only a stantard (release mode) configuration for now.  
### Still wrap it in if's.
###

!IF "$(CFG)" == "Standard"

ALL : "$(RUNDIR)\$(EXENAM).exe"

CLEAN :
	-@erase "$(OBJDIR)\AEROSOLS.mod"
	-@erase "$(OBJDIR)\Chemistry.mod"
	-@erase "$(OBJDIR)\Chemistry.obj"
	-@erase "$(OBJDIR)\Condensation.mod"
	-@erase "$(OBJDIR)\Condensation.obj"
	-@erase "$(OBJDIR)\dgaus8.obj"
	-@erase "$(OBJDIR)\$(DRIVER).obj"
	-@erase "$(OBJDIR)\Gibbs-DuhemEvaluationFunction.obj"
	-@erase "$(OBJDIR)\GridPointFields.mod"
	-@erase "$(OBJDIR)\GridPointFields.obj"
	-@erase "$(OBJDIR)\InfrastructuralCode.mod"
	-@erase "$(OBJDIR)\InfrastructuralCode.obj"
	-@erase "$(OBJDIR)\LagrangianAerosols.obj"
	-@erase "$(OBJDIR)\lsodes.obj"
	-@erase "$(OBJDIR)\indexx.obj"
	-@erase "$(OBJDIR)\PPND16.obj"
	-@erase "$(OBJDIR)\CDFNor.obj"
	-@erase "$(OBJDIR)\ModelParameters.mod"
	-@erase "$(OBJDIR)\ModelParameters.obj"
	-@erase "$(OBJDIR)\ParcelModels.mod"
	-@erase "$(OBJDIR)\ParcelModels.obj"
	-@erase "$(OBJDIR)\ReadInputDeck.obj"
	-@erase "$(OBJDIR)\StepMELAM.mod"
	-@erase "$(OBJDIR)\StepMELAM.obj"
	-@erase "$(OBJDIR)\Time.mod"
	-@erase "$(OBJDIR)\Time.obj"
	-@erase "$(OBJDIR)\Coagulation.mod"
	-@erase "$(OBJDIR)\Coagulation.obj"
	-@erase "$(OBJDIR)\OutputRoutines.obj"
	-@erase "$(RUNDIR)\$(EXENAM).exe"

"$(RUNDIR)" :
    if not exist "$(RUNDIR)/$(NULL)" mkdir "$(RUNDIR)"

##################################################
## Choose the Compiler and set flags accordingly #
!IF "$(FORTRANCOMPILER)" == "Intel"

F90=df.exe
F90_PROJ=/compile_only /nologo /architecture:generic /fltconsistency /fpconstant /fpe:3 /math_library:fast /optimize:5 /real_size:64 /integer_size:64 /tune:p6p /warn:nofileopt /module:$(OBJDIR) /object:$(OBJDIR)
F90_OBJS=$(OBJDIR)

!ELSEIF "$(FORTRANCOMPILER)" == "Compaq"

F90=df.exe
### F90_PROJ=/compile_only /nologo /architecture:generic /fltconsistency /fpconstant /fpe:3 /math_library:fast /optimize:5 /real_size:64 /tune:p6p /warn:nofileopt /module:$(OBJDIR) /object:$(OBJDIR)

F90_PROJ=/check:bounds /compile_only /debug:full /nologo /real_size:64 /traceback /warn:argument_checking /warn:nofileopt /module:$(OBJDIR) /object:$(OBJDIR)

F90_OBJS=$(OBJDIR)



!ENDIF

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

#####################################################
## Define everything to do with the linking system ##
!MESSAGE
LINK32=link.exe
LINK32_FLAGS= kernel32.lib /nologo /stack:5000000000 /subsystem:console /incremental:no /pdb:"$(RUNDIR)\$(EXENAM).pdb" /machine:I386 /out:"$(RUNDIR)\$(EXENAM).exe"
LINK32_OBJS= \
	"$(OBJDIR)\dgaus8.obj" \
	"$(OBJDIR)\lsodes.obj" \
	"$(OBJDIR)\PPND16.obj" \
	"$(OBJDIR)\CDFNor.obj" \
	"$(OBJDIR)\indexx.obj" \
	"$(OBJDIR)\GridPointFields.obj" \
	"$(OBJDIR)\InfrastructuralCode.obj" \
	"$(OBJDIR)\Chemistry.obj" \
	"$(OBJDIR)\Gibbs-DuhemEvaluationFunction.obj" \
	"$(OBJDIR)\LagrangianAerosols.obj" \
	"$(OBJDIR)\Condensation.obj" \
	"$(OBJDIR)\OutputRoutines.obj" \
	"$(OBJDIR)\Driver.obj" \
	"$(OBJDIR)\ModelParameters.obj" \
	"$(OBJDIR)\ParcelModels.obj" \
	"$(OBJDIR)\ReadInputDeck.obj" \
	"$(OBJDIR)\Coagulation.obj" \
	"$(OBJDIR)\StepMELAM.obj" \
	"$(OBJDIR)\Time.obj"

"$(RUNDIR)\$(EXENAM).exe" : "$(RUNDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF

################################################
# Now define the locations of the source files #

!IF "$(CFG)" == "Standard" 

SOURCE=$(SRCDIR)\ComputationalSubroutines\dgaus8.f
"$(OBJDIR)\dgaus8.obj" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\AEROSOLS.mod"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\ComputationalSubroutines\lsodes.f
"$(OBJDIR)\lsodes.obj" : $(SOURCE) "$(OBJDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\ComputationalSubroutines\indexx.f
"$(OBJDIR)\indexx.obj" : $(SOURCE) "$(OBJDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\ComputationalSubroutines\PPND16.f
"$(OBJDIR)\PPND16.obj" : $(SOURCE) "$(OBJDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\ComputationalSubroutines\CDFNor.f90
"$(OBJDIR)\CDFNor.obj" : $(SOURCE) "$(OBJDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\ModelParameters.f90
F90_MODOUT=\
	"ModelParameters"
"$(OBJDIR)\ModelParameters.obj"		"$(OBJDIR)\ModelParameters.mod" : $(SOURCE) "$(OBJDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\EulerianEnvironment\GridPointFields.f90
F90_MODOUT=\
	"GridPointFields"
"$(OBJDIR)\GridPointFields.obj"		"$(OBJDIR)\GridPointFields.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\ModelParameters.mod" "$(SRCDIR)\EulerianEnvironment\GasPhaseChemistry.h" "$(SRCDIR)\EulerianEnvironment\Initialize.h" "$(SRCDIR)\EulerianEnvironment\InterpolateGrid.h" "$(SRCDIR)\EulerianEnvironment\Pressure.h" "$(SRCDIR)\EulerianEnvironment\SetGrid.h" "$(SRCDIR)\EulerianEnvironment\Temperature.h" "$(SRCDIR)\EulerianEnvironment\Water.h"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\InfrastructuralCode\InfrastructuralCode.f90
F90_MODOUT=\
	"InfrastructuralCode"
"$(OBJDIR)\InfrastructuralCode.obj"	"$(OBJDIR)\InfrastructuralCode.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\ModelParameters.mod" "$(SRCDIR)\InfrastructuralCode\Counters.h" "$(SRCDIR)\InfrastructuralCode\InputOutputCommands.h" "$(SRCDIR)\InfrastructuralCode\StringCommands.h"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\Chemistry\Chemistry.f90
F90_MODOUT=\
	"Chemistry"
"$(OBJDIR)\Chemistry.obj"	"$(OBJDIR)\Chemistry.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\ModelParameters.mod" "$(OBJDIR)\Time.mod" "$(SRCDIR)\Chemistry\AqPhaseChemistryInits.h" "$(SRCDIR)\Chemistry\AqPhaseReactionIntegration.h" "$(SRCDIR)\Chemistry\ChemicalPropertyRoutines.h" "$(SRCDIR)\Chemistry\ChemistryDerivativeStructures.h" "$(SRCDIR)\Chemistry\ChemistryParametersInitRoutines.h" "$(SRCDIR)\Chemistry\EvolveGasChemistry.h" "$(SRCDIR)\Chemistry\GasPhaseChemistryInits.h" 
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE="$(SRCDIR)\ComputationalSubroutines\Gibbs-DuhemEvaluationFunction.f90"
"$(OBJDIR)\Gibbs-DuhemEvaluationFunction.obj" : $(SOURCE) "$(OBJDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\LagrangianParticles\LagrangianAerosols.f90
F90_MODOUT=\
	"Aerosols"
"$(OBJDIR)\LagrangianAerosols.obj"	"$(OBJDIR)\AEROSOLS.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\Chemistry.mod" "$(OBJDIR)\ModelParameters.mod" "$(OBJDIR)\Time.mod" "$(SRCDIR)\LagrangianParticles\InitializationAndSamplingRoutines.h" "$(SRCDIR)\LagrangianParticles\ParticleAttributes.h" "$(SRCDIR)\LagrangianParticles\ParticleStructure.h" "$(SRCDIR)\LagrangianParticles\SCTailorMakeParticle.h" "$(SRCDIR)\LagrangianParticles\Thermodynamics.h" "$(SRCDIR)\LagrangianParticles\SortAerosol.h"
 	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\Condensation\Condensation.f90
F90_MODOUT=\
	"Condensation"
"$(OBJDIR)\Condensation.obj"	"$(OBJDIR)\Condensation.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\AEROSOLS.mod" "$(OBJDIR)\Chemistry.mod" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\ModelParameters.mod" "$(SRCDIR)\Condensation\CondensationInitialization.h" "$(SRCDIR)\Condensation\CondensationIntegrator.h" "$(SRCDIR)\Condensation\CondensationRelatedFunctions.h" 
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\Coagulation\Coagulation.f90
F90_MODOUT=\
	"Coagulation"
"$(OBJDIR)\Coagulation.obj" 	"$(OBJDIR)\Coagulation.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\AEROSOLS.mod" "$(OBJDIR)\OutputRoutines.mod" "$(OBJDIR)\ModelParameters.mod" "$(OBJDIR)\Time.mod" "$(SRCDIR)\Coagulation\CoagulationKernels.h" "$(SRCDIR)\Coagulation\CoagulationSC.h"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\OutputRoutines\OutputRoutines.f90
F90_MODOUT=\
	"OutputRoutines"
"$(OBJDIR)\OutputRoutines.obj"	"$(OBJDIR)\OutputRoutines.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\AEROSOLS.mod" "$(OBJDIR)\Chemistry.mod" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\ModelParameters.mod" "$(SRCDIR)\OutputRoutines\OutputRoutines.f90" 
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\StepMELAM.f90
F90_MODOUT=\
	"StepMELAM"
"$(OBJDIR)\StepMELAM.obj"	"$(OBJDIR)\StepMELAM.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\AEROSOLS.mod"  "$(OBJDIR)" "$(OBJDIR)\Time.mod" "$(OBJDIR)\Chemistry.mod" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\ModelParameters.mod" "$(SRCDIR)\OutputRoutines\OutputRoutines.f90" 
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\ParcelModels\ParcelModels.f90
F90_MODOUT=\
	"ParcelModels"
"$(OBJDIR)\ParcelModels.obj" 	"$(OBJDIR)\ParcelModels.mod": $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\Coagulation.mod" "$(OBJDIR)\Chemistry.mod" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\AEROSOLS.mod" "$(OBJDIR)\OutputRoutines.mod" "$(OBJDIR)\Condensation.mod" "$(OBJDIR)\ModelParameters.mod" "$(OBJDIR)\Time.mod" "$(SRCDIR)\ParcelModels\ConstantUpdraft.h" "$(SRCDIR)\ParcelModels\EmanuelUpdraft.h" "$(SRCDIR)\ParcelModels\ParcelSupportingFunctions.h"   "$(OBJDIR)\StepMELAM.mod"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\$(DRIVER).f90
"$(OBJDIR)\$(DRIVER).obj" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\Coagulation.mod" "$(OBJDIR)\Chemistry.mod" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\AEROSOLS.mod" "$(OBJDIR)\OutputRoutines.mod" "$(OBJDIR)\Condensation.mod" "$(OBJDIR)\ModelParameters.mod" "$(OBJDIR)\Time.mod" "$(OBJDIR)\ParcelModels.mod" "$(OBJDIR)\StepMELAM.mod" 
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\ReadInputDeck.f90
"$(OBJDIR)\ReadInputDeck.obj" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\Chemistry.mod" "$(OBJDIR)\InfrastructuralCode.mod" "$(OBJDIR)\GridPointFields.mod" "$(OBJDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=$(SRCDIR)\Time.f90
F90_MODOUT=\
	"Time"
"$(OBJDIR)\Time.obj"	"$(OBJDIR)\Time.mod" : $(SOURCE) "$(OBJDIR)" "$(OBJDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)

!ENDIF 