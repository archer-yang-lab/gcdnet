! DESCRIPTION:
!
!    These functions are minor modifications from the `glmnet` package:
!
!    Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
!    Regularization Paths for Generalized Linear Models via Coordinate Descent.
!    Journal of Statistical Software, 33(1), 1-22.
!    URL: https://www.jstatsoft.org/v33/i01/.
!
! --------------------------------------------------------------------------
! standard: An auxiliary function for standardize x matrix.
! --------------------------------------------------------------------------
!
! USAGE:
!
! CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
!
! INPUT ARGUMENTS:
!
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an
!                     observation vector.
!    ju(nvars) = flag of predictor variables
!                ju(j) = 0 => this predictor has zero variance
!                ju(j) = 1 => this predictor does not have zero variance
!    isd = standarization flag:
!          isd = 0 => do not standardize predictor variables
!          isd = 1 => standardize predictor variables
!          NOTE: no matter isd is 1 or 0, matrix x is always centered by column.
!                That is, col.mean(x) = 0.
!
! OUTPUT:
!
!    x(nobs, nvars) = standarized matrix x
!    xmean(nvars) = column mean of x matrix
!    xnorm(nvars) = column standard deviation of x matrix
!    maj(nvars) = column variance of x matrix
!
! --------------------------------------------------------------------------
! chkvars: An auxiliary function for variable check.
! --------------------------------------------------------------------------
!
! USAGE:
!
! CALL chkvars(nobs, nvars, x, ju)
!
! INPUT ARGUMENTS:
!
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an
!                     observation vector.
!    y(nobs) = response variable. This argument should be a two-level factor
!              {-1, 1} for classification.
!
! OUTPUT:
!
!    ju(nvars) = flag of predictor variables
!                ju(j) = 0 => this predictor has zero variance
!                ju(j) = 1 => this predictor does not have zero variance
!

! --------------------------------------------------------------------------
SUBROUTINE standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
  ! ------------------------------------------------------------------------
  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: nobs
  INTEGER :: nvars
  INTEGER :: isd
  INTEGER :: intr
  INTEGER :: ju(nvars)
  DOUBLE PRECISION :: xmsq
  DOUBLE PRECISION :: xvar
  DOUBLE PRECISION :: x(nobs, nvars)
  DOUBLE PRECISION :: xmean(nvars)
  DOUBLE PRECISION :: xnorm(nvars)
  DOUBLE PRECISION :: maj(nvars)
  ! -------- LOCAL DECLARATIONS -------- !
  INTEGER :: j
  ! -------- STANDARDIZATION -------- !
  IF (intr == 0) THEN
     DO j = 1, nvars
        IF (ju(j) == 1) THEN
           xmean(j) = 0.0D0
           maj(j) = DOT_PRODUCT(x(:,j), x(:,j)) / nobs
           IF (isd == 1) THEN
              xmsq = (SUM(x(:,j)) / nobs)**2
              xvar = maj(j) - xmsq
              xnorm(j) = SQRT(xvar)
              x(:,j) = x(:,j) / xnorm(j)
              maj(j) = 1.0D0 + xmsq / xvar
           END IF
        END IF
     END DO
  ELSE
     DO j = 1, nvars
        IF (ju(j) == 1) THEN
           xmean(j) = SUM(x(:,j)) / nobs ! MEAN
           x(:,j) = x(:,j) - xmean(j)
           maj(j) = DOT_PRODUCT(x(:,j), x(:,j)) / nobs
           IF (isd == 1) THEN
              xnorm(j) = SQRT(maj(j)) ! STANDARD DEVIATION
              x(:,j) = x(:,j) / xnorm(j)
              maj(j) = 1.0D0
           END IF
        END IF
     END DO
  END IF
END SUBROUTINE standard


! -------------------------------------------------------------------------
SUBROUTINE chkvars(nobs, nvars, x, ju)
  ! -----------------------------------------------------------------------
  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: nobs
  INTEGER :: nvars
  INTEGER :: ju(nvars)
  DOUBLE PRECISION :: x(nobs, nvars)
  ! -------- LOCAL DECLARATIONS -------- !
  INTEGER :: i
  INTEGER :: j
  DOUBLE PRECISION :: t
  ! -------- VARIABLE CHECKING -------- !
  DO j = 1, nvars
     ju(j) = 0
     t = x(1, j)
     DO i = 2, nobs
        IF (x(i, j) /= t) THEN
           ju(j) = 1
           EXIT
        END IF
     END DO
  END DO
END SUBROUTINE chkvars
