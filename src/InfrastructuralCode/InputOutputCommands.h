!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! InputOutputCommands.h
!! Errors, Warnings, and Transcripts are handled here, as well as
!! grabbing an input line from an ASCII file 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!! 01/04  2013   Matt Alvarado     Expanded length of allowed lines from 
!!                                   512 to 1024 characters. Needed to 
!!                                   allow 50 photo rxns.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	!!
!! 1. FUNCTION GetLine(FileHandle) RESULT (OutLine)				!!
!! 2. SUBROUTINE Error(ErrorMessage)							!!
!! 3. SUBROUTINE Warn(WarningMessage)							!!
!! 4. SUBROUTINE Transcript(Message)							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This Function Gets the next real line stripped of leading	!!
!! whitespace and trailing comments								!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER (len=1024) FUNCTION GetLine(FileHandle) RESULT (OutLine)

	USE ModelParameters

	IMPLICIT NONE
	INTEGER :: i, FileHandle, Status

	!! Take a new line from the input deck that is not a comment and has some data
5	READ(FileHandle,"(a)",iostat=Status),OutLine
        !WRITE(*,*) "Line Read"

	IF (Status < 0) THEN
		OutLine = EOF					    ! This is an END OF FILE, defined in ModelParameters
	ELSE

	!! Get rid of trailing and leading spaces / tabs
	OutLine = TRIM(StripToken(OutLine))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Strip leading & trailing whitespace and comments	    !
	IF (LEN_TRIM(OutLine) == 0) GOTO 5			    ! A Blank Line
	!!							    !
	IF (OutLine(1:1) == "!") &				    !
		GOTO 5						    ! A Comment Line
	!!							    !
	if (INDEX(OutLine,"!") > 1) &				    ! Remove trailing comments
		OutLine = OutLine(1:INDEX(OutLine,"!")-1)	    !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ENDIF

END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Return an error message that stops the code !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Error(ErrorMessage) 

	IMPLICIT NONE

	!! Input Variables
	CHARACTER*(*) :: ErrorMessage
	CHARACTER (len = 8)     :: ErrorDate
	CHARACTER (len = 10)    :: ErrorTime

	CALL Warn ("ERROR! "//ErrorMessage//" Stopping...")
	CALL Date_and_Time(ErrorDate, ErrorTime)
	CALL Transcript ("*********************************************")
	CALL Transcript ("** Error occured at                         **")
	CALL Transcript ("** "//ErrorTime(1:2)//":"//ErrorTime(3:4)//":"//ErrorTime(5:6)//" on "//ErrorDate(5:6)//		&
					  "/"//ErrorDate(7:8)//"/"//ErrorDate(1:4)//" "//"                 "//"**")
	CALL Transcript ("*********************************************")

	STOP

END SUBROUTINE Error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Print a Warning to the user !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Warn(WarningMessage) 

	IMPLICIT NONE

	CHARACTER*(*) :: WarningMessage

	!! Internal Variables
	INTEGER :: WarnLen, I

	WarnLen = LEN_TRIM(StripToken(WarningMessage))

	CALL Transcript ("")
	CALL Transcript ("*******************************************************")

	!! Break Line if Necessary
10	IF (WarnLen > 45) THEN
		DO I = 40,WarnLen
			IF (WarningMessage(i:i) .EQ. " " .OR. WarningMessage(i:i) .EQ. "	") EXIT
		END DO

		CALL Transcript ("*** "//TRIM(WarningMessage(1:I)))
		WarningMessage = StripToken(WarningMessage(I+1:WarnLen))
		WarnLen = LEN_TRIM(WarningMessage)
		IF (WarnLen > 45) GOTO 10
	END IF

	IF (WarnLen > 0) CALL Transcript ("*** "//TRIM(WarningMessage))

	CALL Transcript ("*******************************************************")
	CALL Transcript ("")


END SUBROUTINE Warn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Record a string to the screen and to the transcript !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Transcript(Message)

	USE ModelParameters, ONLY : TranscriptFH
	IMPLICIT NONE

	CHARACTER*(*) :: Message

	IF (TranscriptFH > 0) WRITE (TranscriptFH,*) TRIM(Message)
	PRINT *, TRIM(Message) !Keep Standard output for DEMMUCOM (MJA, 040207)

END SUBROUTINE Transcript
