! --------------------------------------------------------------------------
! standard.f90: An auxiliary function for standardize x matrix.
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
