! --------------------------------------------------------------------------
! hsvmlassoNET.f90: the GCD algorithm for HHSVM.
! --------------------------------------------------------------------------
! 
! USAGE:
! 
! call hsvmlassoNET (delta, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
! & pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
! & nbeta, alam, npass, jerr)
! 
! INPUT ARGUMENTS:
! 
!    delta = the parameter in HHSVM model.
!    lam2 = regularization parameter for the quadratic penalty of the coefficients
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    y(nobs) = response variable. This argument should be a two-level factor {-1, 1} 
!            for classification.
!    jd(jd(1)+1) = predictor variable deletion flag
!                  jd(1) = 0  => use all variables
!                  jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
!    pf(nvars) = relative L1 penalties for each predictor variable
!                pf(j) = 0 => jth variable unpenalized
!    pf2(nvars) = relative L2 penalties for each predictor variable
!                pf2(j) = 0 => jth variable unpenalized
!    dfmax = limit the maximum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the maximum number of variables ever to be nonzero. 
!           For example once beta enters the model, no matter how many 
!           times it exits or re-enters model through the path, it will 
!           be counted only once. 
!    nlam = the number of lambda values
!    flmin = user control of lambda values (>=0)
!            flmin < 1.0 => minimum lambda = flmin*(largest lambda value)
!            flmin >= 1.0 => use supplied lambda values (see below)
!    ulam(nlam) = user supplied lambda values (ignored if flmin < 1.0)
!    eps = convergence threshold for coordinate majorization descent. 
!          Each inner coordinate majorization descent loop continues 
!          until the relative change in any coefficient is less than eps.
!    isd = standarization flag:
!          isd = 0 => regression on original predictor variables
!          isd = 1 => regression on standardized predictor variables
!          Note: output solutions always reference original
!                variables locations and scales.
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
! 
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    b0(nalam) = intercept values for each solution
!    beta(pmax, nalam) = compressed coefficient values for each solution
!    ibeta(pmax) = pointers to compressed coefficients
!    nbeta(nalam) = number of compressed coefficients for each solution
!    alam(nalam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr  = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                    jerr < 7777 => memory allocation error
!                    jerr = 7777 => all used predictors have zero variance
!                    jerr = 10000 => maxval(vp) <= 0.0
!           jerr < 0 => non fatal error - partial output:
!                    Solutions for larger lambdas (1:(k-1)) returned.
!                    jerr = -k => convergence for kth lambda value not reached
!                           after maxit (see above) iterations.
!                    jerr = -10000-k => number of non zero coefficients along path
!                           exceeds pmax (see above) at kth lambda value.
! 
! LICENSE: GNU GPL (version 2 or later)
! 
! AUTHORS:
!    Yi Yang (yiyang@umn.edu), Yuwen Gu (guxxx192@umn.edu) 
!    and Hui Zou (hzou@stat.umn.edu).
!    School of Statistics, University of Minnesota.
! 
! REFERENCES:
!    Yang, Y. and Zou, H. (2012). An Efficient Algorithm for Computing 
!    The HHSVM and Its Generalizations.
!    Journal of Computational and Graphical Statistics, 22, 396-415.


SUBROUTINE hsvmlassoNET (delta, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju(jd(2:(jd(1)+1))) = 0
      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
      CALL hsvmlassoNETpath (delta, lam2, maj, nobs, nvars, x, y, ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr)
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
         IF (isd == 1) THEN
            DO j = 1, nk
               beta (j, l) = beta (j, l) / xnorm (ibeta(j))
            END DO
         END IF
         b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
        & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE hsvmlassoNET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ---------------- hsvmlassoNETpath ---------------- !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE hsvmlassoNETpath (delta, lam2, maj, nobs, nvars, x, y, ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
        ! - - - arg types - - -
      DOUBLE PRECISION, PARAMETER :: big = 9.9D30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0D-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
      INTEGER :: i
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
      ! - - - begin main program - - - !
      ! - - - allocate variables - - - !
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! ---------- Initial setup ---------- !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      r = 0.0D0
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = maj / delta
      IF (flmin < 1.0D0) THEN
         flmin = Max(mfl, flmin)
         alf = flmin ** (1.0D0/(DBLE(nlam)-1.0D0))
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ----------------- lambda loop ------------------- !
      ! ------------ This is the outest loop ------------ !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO l = 1, nlam !!! begin lambda loop (the outest loop)
        IF (flmin >= 1.0D0) THEN
          al = ulam (l)
        ELSE
          IF (l > 2) THEN
            al = al * alf
          ELSE IF (l == 1) THEN
            al = big
          ELSE IF (l == 2) THEN
            al = 0.0D0
            DO i = 1, nobs
              IF (r(i) > 1.0D0) THEN
                dl (i) = 0.0D0
              ELSE IF (r(i) <= (1-delta)) THEN
                dl (i) = - 1.0D0
              ELSE
                dl (i) = (r(i)-1.0D0) / delta
              END IF
            END DO
            DO j = 1, nvars
              IF (ju(j) /= 0) THEN
                IF (pf(j) > 0.0D0) THEN
                  u = dot_product(dl * y, x(:, j))
                  al = Max(al, Abs(u)/pf(j))
                END IF
              END IF
            END DO
            al = al * alf / nobs
          END IF
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ------------------ outer loop -------------------- !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO !!! begin outer loop
          oldbeta(0) = b(0)
          IF (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ----------------- middle loop -------------------- !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          DO
            npass = npass + 1
            dif = 0.0D0
            DO k = 1, nvars !!! begin updating beta's
              IF (ju(k) /= 0) THEN
                oldb = b (k)
                u = 0.0D0
                DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                    dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                    dl (i) = - 1.0D0
                  ELSE
                    dl (i) = (r(i)-1.0D0) / delta
                  END IF
                  u = u + dl(i) * y (i) * x (i, k)
                END DO
                u = maj(k) * b(k) - u / nobs
                v = al * pf (k)
                v = Abs(u) - v
                IF (v > 0.0D0) THEN
                  b(k) = SIGN(v, u) / (maj(k) + pf2(k) * lam2)
                ELSE
                  b(k) = 0.0D0
                END IF

                d = b(k) - oldb
                IF (ABS(d) > 0.0D0) THEN
                  dif = Max (dif, d**2)
                  r = r + y * x (:, k) * d
                  IF (mm(k) == 0) THEN
                    ni = ni + 1 !!! number of variables ever entered the model
                    IF (ni > pmax) EXIT
                    mm(k) = ni
                    m(ni) = k !!! indicate which coefficient is non-zero
                  END IF
                END IF
              END IF
            END DO !!! end updating beta's
            d = 0.0D0
            DO i = 1, nobs !!! begin updating beta0
              IF (r(i) > 1.0D0) THEN
                dl (i) = 0.0D0
              ELSE IF (r(i) <= (1-delta)) THEN
                dl (i) = - 1.0D0
              ELSE
                dl (i) = (r(i)-1.0D0) / delta
              END IF
              d = d + dl (i) * y (i)
            END DO !!! end updating beta0
            d = - 0.5D0 * delta * d / nobs
            IF (ABS(d) > 0.0D0) THEN
              b (0) = b (0) +  d
              r = r + y * d
              dif = Max(dif, d**2)
            END IF
            IF (ni > pmax) EXIT
            IF (dif < eps * delta) EXIT
            IF(npass > maxit) THEN
                 jerr=-l
                 RETURN
            ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! -----------------inner loop ------------------- !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO !!! begin inner loop
              npass = npass + 1
              dif = 0.0D0
              DO j = 1, ni
                k = m(j)
                oldb = b(k)
                u = 0.0D0
                DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                    dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                    dl(i) = - 1.0D0
                  ELSE
                    dl(i) = (r(i)-1.0D0) / delta
                  END IF
                  u = u + dl(i) * y(i) * x(i, k)
                END DO
                u = maj(k) * b(k) - u / nobs
                v = al * pf(k)
                v = Abs(u) - v
                IF (v > 0.0D0) THEN
                  b(k) = sign(v, u) / (maj(k) + pf2(k) * lam2)
                ELSE
                  b(k) = 0.0D0
                END IF
                d = b(k) - oldb
                IF (Abs(d) > 0.0D0) THEN
                  dif = Max (dif, d**2)
                  r = r + y * x(:, k) * d
                END IF
              END DO
              d = 0.0D0
              DO i = 1, nobs
                IF (r(i) > 1.0D0) THEN
                  dl (i) = 0.0D0
                ELSE IF (r(i) <= (1-delta)) THEN
                  dl (i) = - 1.0D0
                ELSE
                  dl (i) = (r(i)-1.0D0) / delta
                END IF
                d = d + dl (i) * y (i)
              END DO
              d = - 0.5D0 * delta * d / nobs
              IF (ABS(d) > 0.0D0) THEN
                b (0) = b (0) + d
                r = r + y * d
                dif = Max (dif, d**2)
              END IF
              IF (dif < eps * delta) EXIT
              IF(npass > maxit) THEN
                   jerr=-l
                   RETURN
              ENDIF
            END DO !!! end inner loop
          END DO !!! end middle loop
          IF (ni > pmax) EXIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! -------------- final check ---------------- !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          vrg = 1
          IF ((b(0)-oldbeta(0))**2 >= eps) vrg = 0
          DO j = 1, ni
            IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
              vrg = 0
              EXIT
            END IF
          END DO
          IF (vrg == 1) EXIT
        END DO !!! end outer loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ----------- final update & save results ------------ !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (ni > pmax) THEN
          jerr = - 10000 - l
          EXIT
        END IF
        IF (ni > 0) beta (1:ni, l) = b (m(1:ni))
        nbeta (l) = ni
        b0 (l) = b (0)
        alam (l) = al
        nalam = l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
        me = count (ABS(beta(1:ni, l)) > 0.0D0)
        IF (me > dfmax) EXIT
      END DO !!! end lambda loop (the outest loop)
      DEALLOCATE (b, oldbeta, r, mm)
      RETURN
END SUBROUTINE hsvmlassoNETpath
