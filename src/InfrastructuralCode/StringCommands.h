!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! StringCommands.h
!! This module contains subroutines that deal with	
!! string parsing, tokenizing, stripping, etcetera.	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY								!!
!!										!!
!! Month  Year   Name              Description					!!
!! 07     2006   Matt Alvarado     Began Update History			        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		    !!
!! 1. FUNCTION StripToken (InLine) RESULT (OutLine)			    !!
!! 2. FUNCTION INT2STR ( InInt ) RESULT ( OutString )			    !!
!! 3. FUNCTION CountDigits ( InInt ) RESULT ( Digits )			    !!
!! 4. FUNCTION REAL2STR ( InReal, inDecimalPlaces ) RESULT ( OutString )    !!
!! 5. SUBROUTINE GetToken(InLine, Delimitor, OutLine)			    !!
!! 6. FUNCTION STR2INT( InString ) RESULT ( OutInt )			    !!
!! 7. FUNCTION STR2REAL(InString) RESULT (OutReal)			    !!
!! 8. FUNCTION IsInteger( InString ) RESULT ( Flag )			    !!
!! 9. FUNCTION IsReal( InString ) RESULT ( Flag )			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!! This Strips leading and trailing whitespace from a string !!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CHARACTER (len = 512) FUNCTION StripToken (InLine) RESULT (OutLine)
		CHARACTER*(*) :: InLine
		INTEGER       :: i, j

		!! String might be an empty string
		IF (LEN_TRIM(Inline) .EQ. 0) THEN
			OutLine = ""
			RETURN
		END IF

                

		!! Strip Leading Whitespace
		i=1
		DO WHILE (InLine(i:i) == " " .or. InLine(i:i) == "	") ! The second option has a tab character in it
			IF (i == len(InLine)) EXIT			   ! This is a blank line with tabs in it
			i = i+1						   ! Look for first nonblank character
			IF (i==len(InLine)) THEN								
				OutLine = ""										
			END IF													
		END DO														
																	
		!! Strip Trailing Whitespace								
		j=LEN(InLine)												
		DO WHILE (InLine(j:j) == " " .or. InLine(j:j) == "	") ! The second option has a tab character in it
			IF (j .le. 1) EXIT				   ! Trip trailing spaces / tabs off
			j = j-1						   ! 
		END DO							   
																	
		OutLine = InLine(i:j)					   ! Trim Whitespace

		RETURN
	END FUNCTION StripToken

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Create a String from an Integer !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CHARACTER (len = 512) FUNCTION INT2STR ( InInt ) RESULT ( OutString )

		IMPLICIT NONE

		!! External Variables
		INTEGER :: InInt

		!! Internal Variables
		INTEGER :: Digits, I, J, K, L, WorkInt

		!! Initialize 
		OutString = ""
		WorkInt = InInt

		!! Zero is a special case, given that we're playing with LOG's.
		IF (InInt .EQ. 0) THEN
			OutString = "0"
			RETURN
		END IF

		!! Calculate Number of Digits in String
		IF (InInt < 0) WorkInt = -1*WorkInt
		Digits = CountDigits(WorkInt)

		!! Loop Over the Digits to FILL OutString
		DO I = DIGITS-1, 0, -1
			K = Digits - I
			J = FLOOR(WorkInt/10.**I)
			OutString(K:K) = CHAR(J+48)
			WorkInt = WorkInt - J*10.**I
		END DO

		IF (InInt < 0.) OutString = "-"//TRIM(OutString)

		RETURN
	END FUNCTION INT2STR

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Count the number of digits in an integer !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER FUNCTION CountDigits ( InInt ) RESULT ( Digits )

		IMPLICIT NONE

		INTEGER :: InInt
		
		IF (InInt .EQ. 0 .OR. InInt .EQ. 1) THEN
			Digits = 1
			GOTO 20
		END IF

		Digits = CEILING(DLog10(ABS(DFLOAT(InInt))))
		IF (MODULO(InInt,10) .EQ. 0 .AND. Digits .NE. CEILING(dlog10(DFLOAT(InInt+1)))) Digits = Digits + 1

20		RETURN
	END FUNCTION CountDigits

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Create a String from an Real !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CHARACTER (len = 512) FUNCTION REAL2STR ( InReal, inDecimalPlaces ) RESULT ( OutString )

		IMPLICIT NONE

		!! External Variables
		REAL*8			  :: InReal
		INTEGER, OPTIONAL :: inDecimalPlaces

		!! Internal Variables
		INTEGER :: I,LenDec,DecimalPlaces, e
		REAL*8  :: WorkReal
		CHARACTER (len=64) :: Decimals, IntDecimalPlaces

	    !! the minimum dlog10(InReal) before correcting using scientific notation
		REAL*8, PARAMETER :: MinUnSci = -2.

		!! Initialize 
		OutString = ""

		!! See if Number of Digits is Passed
		!! If Not, Install Default Value
		IF (PRESENT(inDecimalPlaces)) THEN
			DecimalPlaces = inDecimalPlaces
		ELSE
			DecimalPlaces = 3
		END IF

		IF (DecimalPlaces < 0) CALL ERROR ("Negative Number of Decimal Places Passed "//&
								"to Real2Str, Which Makes Me Very Unhappy!!! Stopping.")

		IF (DecimalPlaces > 12) CALL ERROR ("Ridiculously Large Number of Decimal Places Suggested "//&
								"to Real2Str, Which Makes Me Very Unhappy!!! Stopping.")

		IF (InReal .LT. 0.0) THEN
			WorkReal = -1.*InReal
		ELSE
			WorkReal = InReal
		END IF
		
		!! Deal with zero or very near zero
		IF (WorkReal .LT. 1.e-40) WorkReal = 0.
		IF (WorkReal .EQ. 0.) THEN
			OutString = "0."
			RETURN
		END IF

		!! Determine if should print in scientific notation on the big end
		e = CountDigits(Floor(WorkReal)) - 1
		IF (e .GT. 1) THEN
			WorkReal = WorkReal / 10.**e

!print *, "DING2",e,WorkReal

		!! Now check to see if it's really small
		ELSE

			IF (WorkReal .LT. 1 ) THEN
			e = FLOOR(dlog10(WorkReal))
			IF (e .LT. MinUnSci) THEN
				WorkReal = WorkReal / 10.**e
			END IF
			END IF

		END IF

		!! See if the decimal representation is too big
		LenDec = LEN_TRIM(INT2STR(INT(ANINT((WorkReal-FLOOR(WorkReal))*10.**DecimalPlaces))))

		!! If rounds up to 1
		IF (LenDec .GT. DecimalPlaces) WorkReal = ANINT(WorkReal)

		!! Print the Integral Portion
		OutString = TRIM(INT2STR(FLOOR(WorkReal)))
		WorkReal  = WorkReal - FLOOR(WorkReal)

		!! If rounds up to one, pass zero beyond the decimal place
		IF (LenDec .GT. DecimalPlaces) THEN
			WorkReal = 0.
			LenDec   = 1
		END IF

		IF (DecimalPlaces > 0) THEN

			!! Install a Decimal Place
			I = LEN_TRIM(OutString) + 1
			OutString(I:I) = "."

			!! Print the After Decimal Places
			Decimals  = TRIM(INT2STR(INT(ANINT(WorkReal*10.**DecimalPlaces))))

			IntDecimalPlaces=""
			DO I = 1, DecimalPlaces-LenDec
				IntDecimalPlaces(I:I) = '0'
			END DO

			!! Create final string
			OutString = TRIM(OutString) // TRIM(IntDecimalPlaces) // TRIM(StripToken(Decimals))
		END IF

		!! Add trimmings onto string
		IF (InReal < 0.) OutString = "-"//TRIM(OutString)
		IF (e .GT. 1 .OR. e .LT. MinUnSci)    OutString = TRIM(OutString)//"e"//INT2STR(e)

		RETURN
	END FUNCTION REAL2STR

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This Grabs a token from a string up to some delimiter			!!
	!! This Strips leading and trailing whitespace from a string		!!
	!!																	!!
	!! The token is considered to begin at the beginning of the string	!!
	!! and ends either at the end of the string or at the delimtor,		!!
	!! whichever comes first											!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 SUBROUTINE GetToken(InLine, Delimitor, OutLine)

		CHARACTER*(*) :: InLine, Delimitor, OutLine
		INTEGER		  :: i, ld, lt

		!! Measure the Delimitor
		ld = LEN(Delimitor)

		!! Measure the Input Line
		lt = LEN_TRIM(InLine)

		!! Chercher the delimitor in the string
		DO i = 1,lt+1-ld
			IF (InLine(i:i+ld-1) .eq. Delimitor) EXIT
		END DO

		!! If the string does not include the delimitor, then it is either
		!! an empty string or a single-token incomplete string.  Either way
		!! we may set the return values as follows
		IF (i .ge. lt+1-ld) THEN
			!! If the string lengths are wrong, errors will result
			IF (LEN(OutLine) .lt. LEN_TRIM(StripToken(InLine))) &
				CALL ERROR("Output string in call to GetToken is of insufficient length")
			OutLine   = TRIM(StripToken(InLine))! The whole line is the token
			InLine  = ""						! The instring is updated to be blank
		ELSE
		!! If the token exists, then return that token and strip it from the 
		!! InLine.
			!! If the string lengths are wrong, errors will result
			IF (LEN(OutLine) .lt. LEN_TRIM(StripToken(InLine(1:i-1)))) &
			CALL ERROR("Output string in call to GetToken is of insufficient length")
			OutLine  = StripToken(InLine(1:i-1))
			InLine   = StripToken(InLine(i+ld:lt)) ! This has to be shorter than the imput line, 
		END IF					       ! which is why this formulation is copacetic

		RETURN

	END SUBROUTINE GetToken
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The following two functions STR2INT and STR2REAL are adaptations of		!!
	!! (very useful) C native routines that are used in analyzing input decks	!!
	!! at various points in this model code.									!!
	!! The are each adapted from GNU licencable routines on www.envpro.ncsc.org !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER FUNCTION STR2INT( InString ) RESULT ( OutInt )

		IMPLICIT NONE

		!! Input typing
		CHARACTER*(*) :: InString

		!! Local variable typing
		CHARACTER (len = 100)   :: CurrString
		INTEGER					:: Sign,		&
								   I, J, K,		&
								   IC, I0

		!! Dump leading and trailing spaces and tabs
		CurrString = TRIM(StripToken(InString))

		!! Make sure the input string is not blank
		IF (LEN_TRIM(CurrString) == 0) CALL ERROR("Empty input string into STR2INT")

		!! I is the position variable with which we nibble
		i = 1

		!! Adjust for sign
		Sign = 1
		IF(CurrString(i:i) .eq. '-' ) THEN		! Sign would be at i = 1 since we used
			Sign = -1							! StripToken()
			i    = i + 1
		END IF
		IF(CurrString(i:i) .eq. '+' ) i = i + 1

		OutInt	= 0         !  accumulate as long as there are digits.
		K		= 0
		I0		= ICHAR('0')

		!! Now total the integer
		DO  J = I, LEN_TRIM(CurrString)
			IC = ICHAR(CurrString(J:J) ) - I0
			IF (IC .LT. 0  .OR.  IC .GT. 9 ) EXIT	! This is a non-integral value.  Bugger out of the loop.
			OutInt = 10 * OutInt +  IC
			k = k+1
		END DO

		!! Either prepare the result or fail.
		IF ( K .GT. 0 ) THEN
			OutInt = SIGN * OutInt
		ELSE
			CALL ERROR("Bad Format in  >>>"//TRIM(CurrString)//"<<<"//&
					   "Call to STR2INT Failed.  Stopping")
		END IF
        
	END FUNCTION

	!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!! Now STR2REAL !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION STR2REAL(InString) RESULT (OutReal)

		IMPLICIT NONE

		!! InPut Typing
		CHARACTER*(*) :: InString
            
		!! Local Variables
		CHARACTER (len = 100)  :: CurrString
		INTEGER				   :: I, N, P, IOS, J, K
		CHARACTER*8			   :: FMT
		CHARACTER (len=15)	   :: Numbers = "0123456789.-+de"

		!! Strip Leading and Trailing Whitespace / Tabs
		CurrString = TRIM(StripToken(InString))
            
		!! Make sure the input string is not blank
		IF (LEN_TRIM(CurrString) == 0) CALL ERROR ("Empty input string into STR2REAL")

		!! Look for Bad Characters and Set a Flag if You See One
		K = 0
		DO J = 1, LEN_TRIM(CurrString) 
			IF (INDEX (Numbers, CurrString(J:J)) == 0) K = 1
		END DO
		IF (K == 1) &								! Non-Numerical Numbers
			CALL ERROR ("STR2REAL Received the Following Confusing String:"//">>"//TRIM(CurrString)//"<<")

		!! Measure the String to get an appropriate format statement for
		!! the string to real read
		P = INDEX(CurrString,'.')
		N = LEN_TRIM(CurrString)
		IF ( P .GT. 0 ) THEN
		WRITE( FMT, 94010 ) N, N - P
		ELSE
		WRITE( FMT, 94010 ) N, 0
		END IF

		!! Push the string value into real form, using the constructed format statement
		READ( CurrString, FMT, IOSTAT = IOS ) OutReal

		!! If there is some error, admit it and stop the program
		IF( IOS .NE. 0 )   &
			CALL ERROR("Error reading REAL from >>>"//TRIM(CurrString)//"<<<"//"In STR2REAL.  Stopping.")
        
		RETURN

94010   FORMAT( '(G', I2.2, '.', I2.2, ')' )
94020   FORMAT( 3A, I7 )

	END FUNCTION

	!! IsInteger and IsReal test to see whether a given character 
	!! string contains a viable integer or real value.
	!!		0 if NO 
	!!		1 if YES
	LOGICAL FUNCTION IsInteger( InString ) RESULT ( Flag )

		IMPLICIT NONE

		!! Input typing
		CHARACTER*(*) :: InString

		!! Local variable typing
		INTEGER				:: I, J, LT
		CHARACTER (len=128) :: CurrString

		!! Dump leading and trailing spaces and tabs
		CurrString = StripToken(InString)
		LT		   = len_trim(CurrString)
		Flag	   = .TRUE.
		J		   = 1

		!! If the string is too long, then it is denied integral status
		IF (LT > 128 .OR. LT ==0 ) THEN
			Flag = .FALSE.
			RETURN
		END IF

		!! Leading Signs are acceptable
		IF(CurrString(J:J) .eq. '-' .OR. CurrString(J:J) .eq. '+') J = J + 1

		!! Check to make sure that each character is a number
		DO I = J, LT
			IF (INDEX ('1234567890', CurrString(I:I)) == 0) Flag = .FALSE.
		END DO

	END FUNCTION IsInteger

	LOGICAL FUNCTION IsReal( InString ) RESULT ( Flag )

		IMPLICIT NONE

		!! Input typing
		CHARACTER*(*) :: InString

		!! Local variable typing
		LOGICAL				:: PassBy
		INTEGER				:: I, J, LT, DOT, EXP
		CHARACTER (len=128) :: CurrString

		!! Dump leading and trailing spaces and tabs
		CurrString = StripToken(InString)
		LT		   = len_trim(CurrString)
		Flag	   = .TRUE.
		J		   = 1
		Dot		   = 0
		Exp		   = 0
		PassBy	   = .FALSE.

		!! If the string is too long, then it is denied real number status
		IF (LT > 128 .OR. LT ==0 ) THEN
			Flag = .FALSE.
			RETURN
		END IF

		!! Leading Signs are acceptable
		IF(CurrString(J:J) .eq. '-' .OR. CurrString(J:J) .eq. '+') J = J + 1

		!! Check to make sure that each character is a number
		DO I = J, LT

 		  IF (PassBy) THEN
			PassBy = .FALSE.
			Cycle
		  END IF

		  !! A real may only have a particular set of constants in it, and some
		  !! may not be repeated (I allow leading zeros here, although maybe that's
		  !! not what people really want).
		  IF (DOT == 0) THEN

			IF (EXP == 0) THEN
				IF (INDEX ('1234567890.de', CurrString(I:I)) == 0) THEN
					Flag = .FALSE.
					EXIT
				END IF
				IF (INDEX ('de', CurrString(I:I)) >  0) THEN
					Exp  = 1
					IF (CurrString(I+1:I+1) .EQ. '+' .OR.	&! Acceptable to have a sign for the scientific notation
						CurrString(I+1:I+1) .EQ. '-') PassBy = .TRUE.
				END IF
				IF (INDEX ('.' , CurrString(I:I)) >  0) Dot  = 1
			ELSE !! EXP == 1
				IF (INDEX ('1234567890.', CurrString(I:I)) == 0) THEN
					Flag = .FALSE.
					EXIT
				END IF
				IF (INDEX ('.', CurrString(I:I)) >  0) Dot  = 1
			END IF
		  ELSE !! Dot == 1
			IF (EXP == 0) THEN
				IF (INDEX ('1234567890de', CurrString(I:I)) == 0) THEN
					Flag = .FALSE.
					EXIT
				END IF
				IF (INDEX ('de', CurrString(I:I)) >  0) THEN
					Exp  = 1
					IF (CurrString(I+1:I+1) .EQ. '+' .OR.	&! Acceptable to have a sign for the scientific notation
						CurrString(I+1:I+1) .EQ. '-') PassBy = .TRUE.
				END IF
			ELSE !! EXP == 1
				IF (INDEX ('1234567890.', CurrString(I:I)) == 0) THEN
					Flag = .FALSE.
					EXIT
				END IF
			END IF
		  END IF
		END DO

	END FUNCTION
