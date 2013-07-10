!! ASP (c), 2004-2013, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! InfrastructuralCode.f90
!! This is the main source file for functions useful
!! in programming, like counters, input and output commands,
!! and commands to convert numbers to strings and vice versa
!! This is the source file linking many header files together.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						    !!
!!								    !!
!! Month  Year   Name              Description			    !!
!! 07     2006   Matt Alvarado     Began Update History		    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES				    !!
!! 1. ModelParameters			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES						    !!
!! 1. InputOutputCommands.h					    !!
!! 2. StringCommands.h						    !!
!! 3. Counters.h					            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE InfrastructuralCode

PUBLIC
PRIVATE FH, AID

!! This is the current file handle
INTEGER :: FH

!! This is the current aerosol id
INTEGER :: AID

!! This is the generic end-of-file statement
CHARACTER (len=33) :: EOF = "END-OF-FILE-012343210-ELIF-FO-DNE"

CONTAINS

INCLUDE "InputOutputCommands.h"
INCLUDE "StringCommands.h"
INCLUDE "Counters.h"

END MODULE InfrastructuralCode
