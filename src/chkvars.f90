! --------------------------------------------------------------------------
! chkvars.f90: An auxiliary function for variable check.
! --------------------------------------------------------------------------
! 
! USAGE:
! 
! call chkvars (nobs, nvars, x, ju)
! 
! INPUT ARGUMENTS:
! 
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    y(no) = response variable. This argument should be a two-level factor {-1, 1} 
!            for classification.
!    
! OUTPUT:
!
!    ju(nvars) = flag of predictor variables
!                ju(j) = 0 => this predictor has zero variance
!                ju(j) = 1 => this predictor does not have zero variance
! 
! LICENSE: GNU GPL (version 2 or later)
! 
! AUTHORS:
!    Yi Yang (yi@stat.umn.edu) and Hui Zou (hzou@stat.umn.edu), 
!    School of Statistics, University of Minnesota.
! 
! REFERENCES:
!    Yang, Y. and Zou, H (2012). An Efficient Algorithm for Computing The HHSVM and Its Generalizations.
!    Journal of Computational and Graphical Statistics. To be accepted after minor revision. 

! --------------------------------------------------
SUBROUTINE chkvars (nobs, nvars, x, ju)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: ju (nvars)
      DOUBLE PRECISION :: x (nobs, nvars)
    ! - - - local declarations - - -
      INTEGER :: i
      INTEGER :: j
      DOUBLE PRECISION :: t
! - - - begin - - -
      DO j = 1, nvars
         ju (j) = 0
         t = x (1, j)
         DO i = 2, nobs
            IF (x(i, j) /= t) THEN
               ju (j) = 1
               EXIT
            END IF
         END DO
      END DO
END SUBROUTINE chkvars
