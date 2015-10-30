!==================================================================================================
MODULE SetPrecision  ! Sets Precision of REAL and INTEGER Variables
							! Written by John Mahaffy
                		! Modified by Gino Banco (2008)
!==================================================================================================
IMPLICIT NONE
        
! Define two parameters that set precision of floating point arithmetic and number of digits in integers respectively
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(13,307)	! floating point precision (double)
!INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(6)			! floating point precision (single)
INTEGER, PARAMETER :: lng = KIND(10000000)					! maximum integer value ("long")

!================================================    
END MODULE SetPrecision
!================================================
