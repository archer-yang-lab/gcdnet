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
