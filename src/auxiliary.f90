! DESCRIPTION: 
!
!    These functions are minor modifications from the glmnet package:
!
!    Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). 
!    Regularization Paths for Generalized Linear Models via Coordinate Descent. 
!    Journal of Statistical Software, 33(1), 1-22. 
!    URL http://www.jstatsoft.org/v33/i01/.
!
! --------------------------------------------------------------------------
! standard: An auxiliary function for standardize x matrix.
! --------------------------------------------------------------------------
!
! USAGE:
! 
! call standard (nobs,nvars,x,ju,isd,xmean,xnorm,maj)   
! 
! INPUT ARGUMENTS:
! 
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    ju(nvars) = flag of predictor variables
!                ju(j) = 0 => this predictor has zero variance
!                ju(j) = 1 => this predictor does not have zero variance
!    isd = standarization flag:
!          isd = 0 => do not standardize predictor variables
!          isd = 1 => standardize predictor variables
!          NOTE: no matter isd is 1 or 0, matrix x is always centered by column. That is, col.mean(x) = 0.
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

! --------------------------------------------------
SUBROUTINE standard(nobs,nvars,x,ju,isd,xmean,xnorm,maj)     
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    INTEGER::  nobs
    INTEGER::nvars
    INTEGER::isd
    INTEGER::ju(nvars)
    DOUBLE PRECISION::  x(nobs,nvars)
    DOUBLE PRECISION::xmean(nvars)
    DOUBLE PRECISION::xnorm(nvars)
    DOUBLE PRECISION::maj(nvars)
    ! - - - local declarations - - -
    INTEGER:: j
! - - - begin - - -                                
    DO j=1,nvars                                  
        IF(ju(j)==1) THEN                         
            xmean(j)=sum(x(:,j))/nobs     !mean                        
            x(:,j)=x(:,j)-xmean(j)    
            maj(j)=dot_product(x(:,j),x(:,j))/nobs                                              
              IF(isd==1) THEN
                xnorm(j)=sqrt(maj(j))    !standard deviation               
                x(:,j)=x(:,j)/xnorm(j)
                maj(j)=1.0D0
            ENDIF                                                        
        ENDIF                                     
    ENDDO                             
END SUBROUTINE standard


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
