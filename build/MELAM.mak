# Microsoft Developer Studio Generated NMAKE File, Based on MELAM.dsp
SHELL = /bin/sh

# F90 compiler options
FFLAGS = f90 -cpp -64 -O3 -OPT:reorg_common=off

# F90 Options -- multiprocessor 
#F90 = $(FFLAGS) -mp -Dmultitask

# F90 Options -- single processor
F90 = $(FFLAGS)

OBJS = \
	CDFNor.o \
	Chemistry.o \
	Coagulation.o \
	Condensation.o \
	dgaus8.o \
	DMiLay.o \
	MIEV0.o \
	MieSupport.o \
	EulerianBoxDriver.o \
	Gibbs-DuhemEvaluationFunction.o \
	GridPointFields.o \
	indexx.o \
	InfrastructuralCode.o \
	LagrangianAerosols.o \
	lsodes.o \
	ModelParameters.o \
	OutputRoutines.o \
	PPND16.o \
	ReadInputDeck.o \
	soa_unidriv.o \
	soa_unidriva.o \
	Time.o \
	unifac.o

#=============================================================================
#  Executable with FAST-J (default)
#=============================================================================
asp 	      : $(OBJS) 
		$(F90)  $(OBJS) -o asp

#==============================================================================
#  Dependencies Listing
#==============================================================================
CDFNor.o 		: .\SourceCode\ComputationalSubroutines\CDFNor.f90
dgaus8.o                : .\SourceCode\ComputationalSubroutines\dgaus8.f
indexx.o		: .\SourceCode\ComputationalSubroutines\indexx.f
lsodes.o		: .\SourceCode\ComputationalSubroutines\lsodes.f
PPND16.o		: .\SourceCode\ComputationalSubroutines\PPND16.f
GridPointFields.o       : .\SourceCode\EulerianEnvironment\GridPointFields.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"GridPointFields"


"$(INTDIR)\GridPointFields.obj"	"$(INTDIR)\GridPointFields.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\GridPointFields.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\SourceCode\InfrastructuralCode\InfrastructuralCode.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"InfrastructuralCode"


"$(INTDIR)\InfrastructuralCode.obj"	"$(INTDIR)\InfrastructuralCode.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\InfrastructuralCode.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\SourceCode\Chemistry\Chemistry.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"Chemistry"


"$(INTDIR)\Chemistry.obj"	"$(INTDIR)\Chemistry.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\GridPointFields.mod" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\ModelParameters.mod" "$(INTDIR)\Time.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\Chemistry.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\SourceCode\LagrangianParticles\DMiLay.f

"$(INTDIR)\DMiLay.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=.\SourceCode\LagrangianParticles\MIEV0.f

"$(INTDIR)\MIEV0.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)

SOURCE=.\SourceCode\LagrangianParticles\MieSupport.f

"$(INTDIR)\MieSupport.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=".\SourceCode\ComputationalSubroutines\Gibbs-DuhemEvaluationFunction.f90"

"$(INTDIR)\Gibbs-DuhemEvaluationFunction.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\SourceCode\LagrangianParticles\LagrangianAerosols.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"Aerosols"


"$(INTDIR)\LagrangianAerosols.obj"	"$(INTDIR)\Aerosols.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\GridPointFields.mod" "$(INTDIR)\Chemistry.mod" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\ModelParameters.mod" "$(INTDIR)\Time.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\LagrangianAerosols.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\SourceCode\Condensation\Condensation.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"Condensation"


"$(INTDIR)\Condensation.obj"	"$(INTDIR)\Condensation.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\Aerosols.mod" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\Chemistry.mod" "$(INTDIR)\GridPointFields.mod" "$(INTDIR)\ModelParameters.mod" "$(INTDIR)\Time.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\Condensation.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\SourceCode\UNIFAC\soa_unidriv.f

"$(INTDIR)\soa_unidriv.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\SourceCode\UNIFAC\soa_unidriva.f

"$(INTDIR)\soa_unidriva.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\SourceCode\UNIFAC\unifac.f

"$(INTDIR)\unifac.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\SourceCode\OutputRoutines\OutputRoutines.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"OutputRoutines"


"$(INTDIR)\OutputRoutines.obj"	"$(INTDIR)\OutputRoutines.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\GridPointFields.mod" "$(INTDIR)\Aerosols.mod" "$(INTDIR)\Chemistry.mod" "$(INTDIR)\Condensation.mod" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\OutputRoutines.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\SourceCode\Coagulation\Coagulation.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"Coagulation"


"$(INTDIR)\Coagulation.obj"	"$(INTDIR)\Coagulation.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\Aerosols.mod" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\GridPointFields.mod" "$(INTDIR)\Chemistry.mod" "$(INTDIR)\OutputRoutines.mod" "$(INTDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\Coagulation.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\SourceCode\DevelopmentDriver.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"


"$(INTDIR)\DevelopmentDriver.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\GridPointFields.mod" "$(INTDIR)\Chemistry.mod" "$(INTDIR)\Aerosols.mod" "$(INTDIR)\Condensation.mod" "$(INTDIR)\OutputRoutines.mod" "$(INTDIR)\Coagulation.mod" "$(INTDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"

!ENDIF 

SOURCE=.\EulerianBoxDriver.f90

"$(INTDIR)\EulerianBoxDriver.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\GridPointFields.mod" "$(INTDIR)\Chemistry.mod" "$(INTDIR)\Aerosols.mod" "$(INTDIR)\Condensation.mod" "$(INTDIR)\OutputRoutines.mod" "$(INTDIR)\Coagulation.mod" "$(INTDIR)\ModelParameters.mod"


SOURCE=.\SourceCode\ModelParameters.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"ModelParameters"


"$(INTDIR)\ModelParameters.obj"	"$(INTDIR)\ModelParameters.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\ModelParameters.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\SourceCode\ReadInputDeck.f90

"$(INTDIR)\ReadInputDeck.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\Chemistry.mod" "$(INTDIR)\InfrastructuralCode.mod" "$(INTDIR)\GridPointFields.mod" "$(INTDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\SourceCode\Time.f90

!IF  "$(CFG)" == "MELAM - Win32 Release"

F90_MODOUT=\
	"Time"


"$(INTDIR)\Time.obj"	"$(INTDIR)\Time.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\ModelParameters.mod"
	$(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "MELAM - Win32 Debug"


"$(INTDIR)\Time.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 




"$(INTDIR)\PPND16.obj" : $(SOURCE) "$(INTDIR)"


Chemistry.o 		: Chemistry.f90			
	Coagulation.o \
	Condensation.o \
	dgaus8.o \
	DMiLay.o \
	MIEV0.o \
	MieSupport.o \
	EulerianBoxDriver.o \
	Gibbs-DuhemEvaluationFunction.o \
	GridPointFields.o \
	indexx.o \
	InfrastructuralCode.o \
	LagrangianAerosols.o \
	lsodes.o \
	ModelParameters.o \
	OutputRoutines.o \
	PPND16.o \
	ReadInputDeck.o \
	soa_unidriv.o \
	soa_unidriva.o \
	Time.o \
	unifac.o
BLKSLV.o                   : BLKSLV.f                   jv_mie.h 
CLDSRF.o                   : CLDSRF.f                   cmn_fj.h CMN_SIZE define.h jv_cmn.h 
CO_strat_pl.o              : CO_strat_pl.f              CMN_SIZE define.h 
EFOLD.o                    : EFOLD.f                    
FLINT.o          


#=============================================================================
#  Other Makefile Commands
#=============================================================================
clean:
	rm -rf *.o *.mod geos rii_files

.SUFFIXES: .f .F .f90 .F90
.f.o:			; $(F90) -c $*.f
.F.o:			; $(F90) -c $*.F
.f90.o:                 ; $(F90) -c -freeform $*.f90 
.F90.o:                 ; $(F90) -c -freeform $*.F90


"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=df.exe
F90_PROJ=/check:bounds /compile_only /debug:full /nologo /real_size:64 /traceback /warn:argument_checking /warn:nofileopt /module:"Debug/" /object:"Debug/" /pdbfile:"Debug/DF60.PDB" 
F90_OBJS=.\Debug/

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\MELAM.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\MELAM.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\MELAM.pdb" /debug /machine:I386 /out:"$(OUTDIR)\MELAM.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\CDFNor.obj" \
	"$(INTDIR)\dgaus8.obj" \
	"$(INTDIR)\indexx.obj" \
	"$(INTDIR)\lsodes.obj" \
	"$(INTDIR)\PPND16.obj" \
	"$(INTDIR)\GridPointFields.obj" \
	"$(INTDIR)\InfrastructuralCode.obj" \
	"$(INTDIR)\Chemistry.obj" \
	"$(INTDIR)\DMiLay.obj" \
	"$(INTDIR)\MIEV0.obj" \
	"$(INTDIR)\MieSupport.obj" \
	"$(INTDIR)\Gibbs-DuhemEvaluationFunction.obj" \
	"$(INTDIR)\LagrangianAerosols.obj" \
	"$(INTDIR)\Condensation.obj" \
	"$(INTDIR)\soa_unidriv.obj" \
	"$(INTDIR)\soa_unidriva.obj" \
	"$(INTDIR)\unifac.obj" \
	"$(INTDIR)\OutputRoutines.obj" \
	"$(INTDIR)\Coagulation.obj" \
	"$(INTDIR)\EulerianBoxDriver.obj" \
	"$(INTDIR)\ModelParameters.obj" \
	"$(INTDIR)\ReadInputDeck.obj" \
	"$(INTDIR)\Time.obj"

"$(OUTDIR)\MELAM.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("MELAM.dep")
!INCLUDE "MELAM.dep"
!ELSE 
!MESSAGE Warning: cannot find "MELAM.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "MELAM - Win32 Release" || "$(CFG)" == "MELAM - Win32 Debug"

	$(F90) $(F90_PROJ) $(SOURCE)




