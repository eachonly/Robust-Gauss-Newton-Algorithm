!******************************************************************
!
! Purpose: Optimize SLS Objective Function with Robust Gauss-Newton Algorithm
!
! Programmer: George Kuczera, Youwei Qin, Dmitri Kavetski
! Created: 21 May 2016 at Newcastle, Australia
! Last modified: 15 July 2018 at Nanjing, China
! Copyright, George Kuczera, Youwei Qin, Dmitri Kavetski, 2018-2023. All rights reserved.
!
! References
! * Qin2018a: Youwei Qin, Kavetski Dmitri, George Kuczera (2018),
!            A Robust Gauss-Newton algorithm for the optimization of hydrological models: From standard Gauss-Newton to Robust Gauss-Newton,
!            Water Resources Research, accept

! * Qin2018b: Youwei Qin, Kavetski Dmitri, George Kuczera (2018),
!            A Robust Gauss-Newton algorithm for the optimization of hydrological models: Benchmarking against industry-standard algorithms,
!            Water Resources Research, accept
!
!******************************************************************
! ---
! Input
!   p:          Number of parameters
!   n:          Number of observations in calibration
!   x0:         Initial parameters
!   xLo:        Lower bounds on parameters
!   xHi(:):     Upper bounds on parameters
!   cnv:        Convergence data structure
!   decFile:    Name of dumpfile    
! ---
! Output
!   info:       Run information data structure
!   x(:):       Final parameters
!   error       Error code, 0 = ok
!   message:    Error message
! ---
! Notes
!   This module follows fairly closely to the pseudocode in Qin2018a,
!   Any issues or bugs, please contact the first author(Email:youwei.qin@uon.edu.au)

MODULE rgnMod
   USE constantsMod, ONLY : ik, rk
   IMPLICIT NONE
   PRIVATE 
   PUBLIC :: setDefaultRgnConvergeSettings, rgnConvType, rgn, rgnInfoType
   LOGICAL, PARAMETER :: NO=.false., YES=.true.
   REAL(rk), PARAMETER :: EPS=EPSILON(0.0_rk)
   
   TYPE rgnConvType
      INTEGER(ik) :: iterMax, noReduction, noRelChangeF, noRelChangePar, fail
      REAL(rk) :: noRelChangeFTol, noRelChangeParTol, tolSafe
      INTEGER(ik) :: dumpResults
      CHARACTER(512) :: logFile
   END TYPE rgnConvType
   
   TYPE rgnInfoType
      INTEGER(ik) :: nIter, termFlag, nEval
      REAL(rk) :: f, cpuTime,objTime
   END TYPE rgnInfoType
   
   TYPE rgnSetType
      REAL(ik) :: gtol, stol, ftol, gtolMin, tol, tolFor, alpha, sigma, c1, rho, beta, hLow, hHiFrac, xScale, fScale
      INTEGER(ik) :: nls
   END TYPE rgnSetType
   TYPE (rgnSetType) :: set
     
   INTERFACE Diag
     MODULE PROCEDURE Diag_ele, Diag_mat
   END INTERFACE

CONTAINS
!
! Initialize the RGN converge variables 
SUBROUTINE setDefaultRgnConvergeSettings (cnvSet, iterMax, dump, logFile, fail)
   TYPE (rgnConvType), INTENT(inout)  :: cnvSet
   INTEGER(ik), INTENT(in), OPTIONAL  :: iterMax, fail
   INTEGER(ik), OPTIONAL              :: dump
   CHARACTER(*), INTENT(in), OPTIONAL :: logFile
   !---
   !
      cnvSet%iterMax = 100; IF (PRESENT(iterMax)) cnvSet%iterMax = iterMax
      
      cnvSet%fail = 100000; IF (PRESENT(fail)) cnvSet%fail = fail
      
      cnvSet%noReduction = 4
      
      cnvSet%noRelChangeFTol = 1.0e-5_rk
      cnvSet%noRelChangeF = 5
      
      cnvSet%noRelChangeParTol = 1.0e-5_rk
      cnvSet%noRelChangePar = 5
	  cnvSet%tolSafe = 1.0e-14_rk
      
      cnvSet%dumpResults = 0; IF (PRESENT(dump)) cnvSet%dumpResults = dump
      
      cnvSet%logFile = 'rgnLog.txt'; IF (PRESENT(logFile)) cnvSet%logFile = logFile
            
END SUBROUTINE setDefaultRgnConvergeSettings
!
! Initialize the RGN constants 
SUBROUTINE setRgnConstants (alpha, beta, nls)
   INTEGER(ik), INTENT(in), OPTIONAL :: nls
   REAL(rk), INTENT(in), OPTIONAL :: alpha, beta
   !---
   !
      set%gtol = 1.0e-10_rk
      set%stol = 1.0e-11_rk
      set%ftol = 1.0e-10_rk
      set%gtolMin = 0.1_rk
      set%tol = 0.1_rk
      set%tolFor = 0.1_rk
      set%alpha = 0.001_rk; IF (PRESENT(alpha)) set%alpha = alpha
      set%sigma = 1.0_rk
      set%c1 = 0.0001_rk
      set%rho = 0.6_rk
      set%nls = 4; IF (PRESENT(nls)) set%nls = nls
      set%beta = 10.0_rk; IF (PRESENT(beta)) set%beta = beta
      set%hLow = 1.0e-8_rk
      set%hHiFrac = 0.5_rk
      set%xScale = 10.0_rk
      set%fScale = 1.0_rk
      
END SUBROUTINE setRgnConstants
!
!
! Robust Gauss-Newton code based on Qin2018a   
SUBROUTINE rgn (objFunc, p, n, x0, xLo, xHi, cnv, x, info, error, message, decFile)
   IMPLICIT NONE
   EXTERNAL :: objFunc
   INTEGER(ik), INTENT(in)            :: p        ! Number of parameters
   INTEGER(ik), INTENT(in)            :: n        ! Number of observations in calibration
   REAL(rk), INTENT(in)               :: x0(:)    ! Initial parameters
   REAL(rk), INTENT(in)               :: xLo(:)   ! Lower bounds on parameters
   REAL(rk), INTENT(in)               :: xHi(:)   ! Upper bounds on parameters
   TYPE (rgnConvType), INTENT(in)     :: cnv      ! Convergence data structure
   TYPE (rgnInfoType), INTENT(out)    :: info     ! Run information data structure
   REAL(rk), INTENT(out)              :: x(:)     ! Final parameters
   INTEGER(ik), INTENT(inout)         :: error    ! error code 0= ok
   CHARACTER(*),INTENT(inout)         :: message  ! error message
   CHARACTER(*), INTENT(in), OPTIONAL :: decFile  ! dumpfile name
   INTERFACE
      SUBROUTINE objFunc (nPar, nSim, x, r, f, timeFunc, error, message)
         USE constantsMod, ONLY: ik, rk
         IMPLICIT NONE
         INTEGER(ik), INTENT(in) :: nPar
         INTEGER(ik),INTENT(in):: nSim
         REAL(rk), INTENT(in) :: x(:)
         REAL(rk), INTENT(out) :: r(:)
         REAL(rk), INTENT(out):: f
         REAL(rk),INTENT(out):: timeFunc   
         INTEGER(ik), INTENT(inout):: error
         CHARACTER(100),INTENT(inout) :: message
   INTEGER(ik)::status
      END SUBROUTINE objFunc
   END INTERFACE
   
   INTEGER(ik) :: nIter, i, j, k, m, nrls, nf, iMax, nr, termCode, noReduction, noRelChangeF, noRelChangePar
   INTEGER(ik), ALLOCATABLE :: as(:)
   
   LOGICAL :: forceRelease, flag_ls, xist
   
   REAL(rk) :: f, fl, fh, fBest, gMaxFree, gMaxThawn, minSingFrac, sig, fredExp, ft, fls, cons, scaledG,         &
               scaledS, scaledfExp, scaledfAct, fredAct, fOldBest, maxRelPar, time(2), gradDf, hessDf
   REAL(rk), ALLOCATABLE :: h(:), r(:), rBest(:), rl(:), rh(:), xl(:), xh(:), xBest(:), Ja(:,:), g(:), He(:,:), xScale(:),  &
                            xt(:), xp(:), delX(:), delXRdc(:), HeRdc(:,:), gRdc(:), xls(:), hLo(:), hHi(:),       &
                            x0ldBest(:), tsv(:), fOptSeries(:), delXAct(:)
   INTEGER(ik):: status
   INTEGER(ik), PARAMETER :: BF=0, BL=1, BFL=2, BH=3, BFH=4,                       &
                             NUL_CON=-1, GRAD_CON=0, SEARCH_CON=1, FRED_CON=2
                                !local parameters
   CHARACTER(*),PARAMETER::procnam="rgnMain"
   REAL(rk):: time4fcall,time4fcallAcc
   CHARACTER(20) :: dfm(4)
   !CHARACTER(100) :: mess
   !----
   error=0                              ! Initialize error flag
   message=""                           ! Initialize message
   time4fcall=0.0_rk;time4fcallAcc=0.0_rk
   ! Allocate work arrays
      ALLOCATE (h(p), r(n), rBest(n), rl(n), rh(n), xl(p), xh(p), xBest(p), Ja(n,p), g(p), He(p,p), as(p), xScale(p), &
                xp(p), xt(p), delX(p), xls(p), hLo(p), hHi(p), x0ldbest(p), fOptSeries(cnv%iterMax), delXAct(p),STAT=status)        
      IF (cnv%dumpResults >= 1) THEN
         OPEN (unit=99, file=cnv%logFile, status='unknown')
         WRITE(dfm(1),'(a,i4,a)') '(a,', p,'g15.7)'
         WRITE(dfm(2),'(a,i4,a)') '(a,', p,'i3)'
         WRITE(dfm(3),'(a,i4,a)') '(33x,', p,'g15.7)'
         WRITE(dfm(4),'(a,i4,a)') '(a,', p,'(i4,11x))'
      END IF
   !
   ! Assign constants
      CALL setRgnConstants
      xScale = set%xScale
      hLo = set%hLow; hHi =  set%hHiFrac*(xHi-xlo)
   !
   ! Initialize
      as = BF  ! Assume all initial parameters are free
      nIter = 0
      h = hHi
      forceRelease = NO
      noReduction = 0; noRelChangeF = 0; noRelChangePar = 0
      info%termFlag = 0; info%nEval = 0
      CALL CPU_TIME (time(1))
      x = x0
      CALL objFunc (nPar=p, nSim=n, x=x, r=rBest, f=f, timeFunc=time4fcall, error=error, message=message); info%nEval = info%nEval + 1; IF (error /=0) GO TO 1
      fBest = f; xBest = x
      time4fcallAcc=time4fcallAcc+time4fcall  
      !CALL userRunTimeMessage ('Starting RGN', -1)
   !-------------
   ! RGN iteration loop
      iterLoop : DO
         nIter = nIter + 1
   !
   ! Save best result from previous iteration
         fOldBest = fBest; x0ldbest = xBest
         !CALL userRunTimeMessage (mess, -1)
         IF (cnv%dumpResults >= 1) THEN
            WRITE(99,'(a)')          '-----------------------------------------'
            WRITE(99,'(a,i5)')       'Iteration No.=                 ', nIter
            WRITE(99,'(a,i5)')       'ObjFun Calls=                  ', info%nEval
            WRITE(99,'(a,g15.7)')    'ObjFun Value f=                 ', f
            WRITE(99,dfm(1))         'Parameter set=                  ', x
            WRITE(99,dfm(4))         'Bound Index=                    ', MERGE(-1_ik, MERGE(0_ik, 1_ik, x(1:p) <= xHi(1:p)), x(1:p) < xLo(1:p))
            WRITE(99,dfm(1))         'Sampling Scale h=               ', h
         END IF
   !
   ! Get Jacobian and update best function result
         xh = x; xl = x; r = rBest
         DO k = 1, p
            xh(k) = x(k) + h(k); xh(k) = MIN(xHi(k), xh(k))
            IF (cnv%dumpResults >= 2) WRITE(99,dfm(1)) 'Forward Jacoian sample point:   ', xh
            CALL objFunc (nPar=p, nSim=n, x=xh, r=rh, f=fh, timeFunc=time4fcall, error=error, message=message); info%nEval = info%nEval + 1; IF (error /=0) GO TO 1; CALL updateBest (fh, xh, rh)
            xl(k) = x(k) - h(k); xl(k) = MAX(xLo(k), xl(k))
            time4fcallAcc=time4fcallAcc+time4fcall  
            IF (cnv%dumpResults >= 2) WRITE(99,dfm(1)) 'Backward Jacobian sample Point: ', xl
            CALL objFunc (nPar=p, nSim=n, x=xl, r=rl, f=fl, timeFunc=time4fcall, error=error, message=message); info%nEval = info%nEval + 1; IF (error /=0) GO TO 1; CALL updateBest (fl, xl, rl)
            time4fcallAcc=time4fcallAcc+time4fcall  
            Ja(:,k) = (rh-rl)/(xh(k)-xl(k))
            xh(k) = x(k); xl(k) = x(k)
            IF (cnv%dumpResults >= 2) WRITE(99,'(a,i3,2g15.7)') 'Jacobian matrix column:          ', k, fh, fl
         END DO
   !
   ! Calculate gradient and Hessian
         DO i = 1, p
            g(i) = DOT_PRODUCT(Ja(:,i),r)
            DO j = 1, i
               He(i,j) = DOT_PRODUCT(Ja(:,i),Ja(:,j))
               He(j,i) = He(i,j)
            END DO
         END DO
   !
   ! Perform active set update
         DO k = 1, p
            IF (x(k) <= xLo(k) + 10.0_rk*EPS*MAX(ABS(xLo(k)),xScale(k))) THEN
               as(k) = MERGE(BFL, BL, g(k) < 0.0_rk)
               CALL updateHess (He, k)
            ELSE IF (x(k) >= xHi(k) - 10.0_rk*EPS*MAX(ABS(xHi(k)),xScale(k))) THEN
               as(k) = MERGE(BFH, BH, g(k) > 0.0_rk)
               CALL updateHess (He, k)
            ELSE
               as(k) = BF
            END IF
         END DO
         
         IF (cnv%dumpResults >= 2) THEN
            WRITE(99,'(a,g15.7)')      'Best objective function value f=', fBest
            WRITE(99,dfm(1))           'Best parameter x=               ', xBest
            WRITE(99,dfm(1))           'Gradient g at parameter x=       ', g
            WRITE(99,dfm(1))           'Hessian at x=                    ', He(1,1)
            DO j = 2, p
               WRITE(99,dfm(3))         He(j,1:j)
            END DO
         END IF         
   !
   ! Determine termination code and hence status of forcRelease
         IF (nIter > 1) THEN
            termCode = NUL_CON
            scaledG = 0.0_rk
            DO k = 1, p
               IF (as(k) == BF .or. as(k) == BFL .or. as(k) == BFH) THEN
                  scaledG = MAX (scaledG, ABS(g(k))*MAX(ABS(x(k)),xScale(k))/MAX(f,set%fScale))
               END IF
            END DO
            scaledS = 0.0_rk
            DO k = 1, p
               scaledS = MAX (scaledS, ABS(delXAct(k))/MAX(ABS(x(k)),xScale(k)))
            END DO
            scaledfExp = ABS(fredExp)/MAX(f,set%fScale)
            scaledfAct = ABS(fredAct)/MAX(f,set%fScale)
            
            IF (scaledG <= set%gtol) THEN
               termCode = GRAD_CON
            ELSE IF (scaledS <= set%stol .and. scaledG <= set%gtolMin) THEN
               termCode = SEARCH_CON
            ELSE IF (scaledfExp <= set%ftol .and. scaledfAct <= set%ftol .and. scaledG <= set%gtolMin) THEN
               termCode = FRED_CON
            END IF
            
            nf = SUM(MERGE(1_ik, 0_ik, as == BF))
            IF (nf == 0) THEN
               forceRelease = YES
            ELSE IF (termCode /= NUL_CON) THEN
               forceRelease = YES
            ELSE
               forceRelease = NO            
            END IF
         END IF
   !
   ! Check conditions for releasing parameters
         nrls = SUM(MERGE(1_ik, 0_ik, as == BFL .or. as == BFH))
         nf = SUM(MERGE(1_ik, 0_ik, as == BF))
         IF (nrls > 0) THEN
            IF (nf > 0) THEN
               gMaxFree = MAXVAL(ABS(g)*MAX(ABS(x),xScale), mask = as == BF)
               DO k = 1, p
                  IF (ABS(g(k))*MAX(ABS(x(k)),xScale(k))>= set%tol*gMaxFree .and. (as(k) == BFL .or. as(k) == BFH)) THEN
                     as(k) = BF
                  END IF
               END DO              
            END IF
            IF (forceRelease) THEN
               iMax = 0
               gMaxThawn = 0.0_rk
               DO k = 1, p
                  IF (ABS(g(k)*MAX(ABS(x(k)),xScale(k))) > gMaxThawn .and. (as(k) == BFL .or. as(k) == BFH)) THEN
                     gMaxThawn = ABS(g(k)*MAX(ABS(x(k)),xScale(k))); iMax = k
                  END IF
               END DO
               IF (iMax > 0) as(iMax) = BF
               DO k = 1, p
                  IF (ABS(g(k)*MAX(ABS(x(k)),xScale(k))) > set%tolFor*gMaxThawn .and. (as(k) == BFL .or. as(k) == BFH)) as(k) = BF
               END DO
            END IF           
         END IF
   !
   ! Solve normal equations after removing non-free parameters
         IF (cnv%dumpResults >= 2) WRITE(99,dfm(4)) 'Active set=                     ', as
         nr = SUM(MERGE(1_ik, 0_ik, as == BF))
         ALLOCATE (HeRdc(nr,nr), delXRdc(nr), gRdc(nr), tsv(nr))
         j = 0
         DO k = 1, p
            IF (as(k) == BF) THEN
               j = j + 1; gRdc(j) = g(k)
               m = 0
               DO i = 1, p
                  IF (as(i) == BF) THEN
                     m = m + 1; HeRdc(m,j) = He(i,k)
                  END IF
               END DO
               HeRdc(j,:) = HeRdc(:,j)
            END IF
         END DO
         minSingFrac = set%alpha*SQRT(EPS)
         CALL svdSolve (m=nr, n=nr, A=HeRdc, b=-gRdc, x=delXRdc, TS=tsv, error=error, message=message, minSingFrac=minSingFrac)
         IF (error /=0) GO TO 1
         j = 0
         DO k = 1, p
            IF (as(k) == BF) THEN
               j = j + 1; delX(k) = delXRdc(j)
            ELSE
               delX(k) = 0.0_rk
            END IF
         END DO
         IF (cnv%dumpResults >= 1) WRITE(99,dfm(1)) 'Truncated SV=                   ', tsv
   !
   ! Project search step onto box constraints
         IF (cnv%dumpResults >= 2) WRITE(99,dfm(1)) 'SVD delX=                       ', delX
         xt = x + delX
         xp = MIN(xHi, MAX(xLo, xt))
         delX = xp - x
         IF (cnv%dumpResults >= 2) THEN
            WRITE(99,dfm(1)) 'Projected delX=                 ', delX
            WRITE(99,dfm(1)) 'Projected xp=                   ', xp
         END IF
   !
   ! Update delXRdc for calculation fredExp
         j = 0
         DO k = 1, p
            IF (as(k) == BF) THEN
               j = j + 1
               delXRdc(j) = delX(k)
            END IF
         END DO
   !
   ! Calculate the expected function reduction with constrained delXRdc
         gradDf = DOT_PRODUCT(gRdc,delXRdc)
         hessDf = DOT_PRODUCT(delXRdc,MATMUL(HeRdc,delXRdc)) 
         fredExp = gradDf + 0.5_rk*hessDf
         DEALLOCATE (HeRdc, gRdc, delXRdc, tsv)
   !
   ! Perform inexact line search
         sig = MIN (set%sigma, 1.0_rk)
         cons = set%c1*MIN(0.0_rk, DOT_PRODUCT(delX,g))
         IF (cnv%dumpResults >= 3) WRITE(99,'(a,g15.7)') 'Cons=                            ', cons
         flag_ls = NO
         DO i = 0, set%nls
            xt = x + sig*delX
            CALL objFunc (nPar=p, nSim=n, x=xt, r=rl, f=ft, timeFunc=time4fcall, error=error, message=message); info%nEval = info%nEval + 1; IF (error /=0) GO TO 1
            time4fcallAcc=time4fcallAcc+time4fcall
            IF (cnv%dumpResults >= 3) THEN
                WRITE(99,dfm(1))           'xt=                             ', xt
                WRITE(99,'(a,g15.7)')      'ft=                             ', ft
                WRITE(99,'(a,g15.7)')      'ft+sig=                         ',ft + sig*cons
            END IF
            IF (ft < f + sig*cons) THEN  
               xls = xt; fls = ft; flag_ls = YES
               IF (cnv%dumpResults >= 1) WRITE(99,'(a,i4,a,g15.7)') 'Line search successful at iteration', i ,' with sigma=', sig
               EXIT
            ELSE
               sig = set%rho*sig
            END IF
         END DO
         IF (.not.flag_ls) THEN
            IF (cnv%dumpResults >= 1) WRITE(99,'(a,i4,a,g15.7)') 'Line search failed'
            fls = f; xls = x        
         END IF
         fredAct = fls - f
         delXAct = xls - x ! ! Save actual search step for termination calculation
   !
   ! Update best variables
         IF (fBest < fls) THEN    ! Jacobian evaluation produced better f than line search
            x = xBest; f = fBest
            flag_ls = YES
         ELSE                     ! Line search produced best f
            x = xls; f = fls
            CALL updateBest (fls, xls, rl)
         END IF
   !
   ! Store the value of best objective function for termination check
         fOptSeries(nIter) = fBest      
  !
   ! Update sampling scale
         IF (flag_ls) THEN
!         IF (flag_ls .and. .not.(noReduction > cnv%fail .or. noRelChangeF > cnv%fail .or. noRelChangePar > cnv%fail)) THEN  !GAK enhance
            h = MIN (set%beta*h, hHi)
         ELSE
            h = MAX (h/set%beta, hLo)
         END IF
   !
   ! Check for convergence
         IF (nIter >= cnv%iterMax) THEN
            info%termFlag = 1; EXIT
         END IF     
         IF (nIter > 1) THEN 
            noRelChangeF = 0
            DO k = MAX(1,nIter - cnv%noRelChangeF + 1), nIter 
                IF (ABS((fOptSeries(k)-fOptSeries(nIter))/(fOptSeries(k)+cnv%tolsafe)) <= cnv%noRelChangeFTol) THEN
                  noRelChangeF = noRelChangeF + 1
                END IF
            END DO
            
            IF (noRelChangeF >= cnv%noRelChangeF) THEN
              info%termFlag = 2; EXIT
            ENDIF
            
			   noReduction = MERGE (noReduction+1_ik, 0_ik, f >= fOldBest)
            IF (noReduction >= cnv%noReduction) THEN
               info%termFlag = 3; EXIT
            END IF
 
            maxRelPar = -HUGE(f)
            DO k = 1, p
               IF (as(k) == BF) THEN
                  maxRelPar = MAX (maxRelPar, ABS((x0ldBest(k)-x(k))/(x0ldBest(k)+cnv%tolsafe)))
               END IF
            END DO
            noRelChangePar = MERGE (noRelChangePar+1_ik, 0_ik, maxRelPar >= 0.0_rk .and. maxRelPar < cnv%noRelChangeParTol)
            IF (noRelChangePar >= cnv%noRelChangePar) THEN
               info%termFlag = 4; EXIT
            END IF
         END IF
      END DO iterLoop
   !
   ! Save optional information
      CALL CPU_TIME (time(2)); info%cpuTime = time(2) - time(1);info%objTime=time4fcallAcc
      info%nIter = nIter; info%f = f
      WRITE(message,'(a,i2,a,g15.6)') 'RGN ended with termination code: ', info%termFlag, ' f=', info%f
      IF (cnv%dumpResults >= 1) THEN
         WRITE(99,'(a,i2)')    '>>>>> RGN ended with termination code: ', info%termFlag
         WRITE(99,'(a,i8)')    '      number of function calls:    ', info%nEval
         WRITE(99,'(a,f10.3)') '      cpu time (sec):               ', info%cpuTime
         CLOSE (unit=99)
      END IF
      RETURN
   !
   ! Error states
1    error = 1; message = "f-"//procnam//"RGN objFunc call failed"
CONTAINS

SUBROUTINE updateBest (f, x, r)
   REAL(rk), INTENT(in) :: x(:), r(:), f
   !---
   !
      IF (f < fBest) THEN
         fBest = f
         xBest = x
         rBest = r
      END IF
END SUBROUTINE updateBest

SUBROUTINE updateHess (Hess, k)
   USE constantsMod,ONLY:zero
   REAL(rk), INTENT(inout) :: Hess(:,:)
   INTEGER(ik),INTENT(in)::k
   !---
   REAL(rk)::diagK
   diagK=Hess(k,k)
   Hess(k,:)=zero; Hess(:,k)=zero
   Hess(k,k)=diagK
END SUBROUTINE updateHess


END SUBROUTINE rgn
!
!
SUBROUTINE svdSolve (m, n, A, b, x, Ainv, S, tS, error, message, minSingFrac, minSingVal, cn)
   ! Solves Ax=b using SVD decomposition followed setting singular values to zero and then back substitution
   INTEGER(ik), INTENT(in) :: m, n
   REAL(rk), INTENT(in) :: A(:,:)
   REAL(rk), INTENT(in), OPTIONAL :: b(:)
   REAL(rk), INTENT(out), OPTIONAL :: x(:)
   REAL(rk), INTENT(out), OPTIONAL :: Ainv(:,:)
   REAL(rk), INTENT(out), OPTIONAL :: S(:)
   REAL(rk), INTENT(out), OPTIONAL :: tS(:)
   INTEGER(ik), INTENT(inout) :: error
   CHARACTER(*),INTENT(inout) :: message
   REAL(rk), INTENT(in), OPTIONAL :: minSingFrac
   REAL(rk), INTENT(out), OPTIONAL :: minSingVal
   REAL(rk), INTENT(out), OPTIONAL ::  cn
   INTEGER(ik) :: i, j, k, nite
   REAL(rk), ALLOCATABLE :: U(:,:), SD(:,:), W(:), V(:,:), tmp(:)
   REAL(rk) :: wMin
   !----
   ! Check consistency of dimension
!      error = 0; message = 'ok'
      
      IF (SIZE(A,1) /= m) THEN; error = 1; message = 'm not same as assumed size in A(m,n)' ; RETURN; END IF
      IF (SIZE(A,2) /= n) THEN; error = 1; message = 'n not same as assumed size in A(m,n)' ; RETURN; END IF
      IF (PRESENT(b) .and. PRESENT(x)) THEN
         IF (SIZE(b) /= n)   THEN; error = 1; message = 'n not same as assumed size in b(n)'   ; RETURN; END IF
         IF (SIZE(x) /= n)   THEN; error = 1; message = 'n not same as assumed size in x(n)'   ; RETURN; END IF
      END IF
      IF (PRESENT(S)) THEN
         IF (SIZE(S) /= n)   THEN; error = 1; message = 'n not same as assumed size in S(n)'   ; RETURN; END IF
      END IF
   !
   ! Perform SVD of A
      IF (ALLOCATED(U)) DEALLOCATE(U);     ALLOCATE (U(m,n))
      IF (ALLOCATED(SD)) DEALLOCATE(SD);   ALLOCATE (SD(n,n))
      IF (ALLOCATED(W)) DEALLOCATE(W);     ALLOCATE (W(n))
      IF (ALLOCATED(V)) DEALLOCATE(V);     ALLOCATE (V(n,n))
      IF (ALLOCATED(tmp)) DEALLOCATE(tmp); ALLOCATE (tmp(n))
   ! Singular value decomposition
      CALL svdDecomp (a=A, u=U, s=SD, v=V, nite=nite)
      IF (error /= 0) GO TO 11
   ! Dealing with U and V, in the opposite direction
      U=-U; V=-V   
      FORALL (n=1:SIZE(V,1)) W(n) = SD(n,n)
      IF (PRESENT(S)) S = W
   !
   ! Zero singular values
      IF (PRESENT(minSingFrac)) THEN
         wMin = minSingFrac
      ELSE
         wMin = SQRT(EPS)
      END IF
      wMin = MAXVAL(W)*wMin
      IF (PRESENT(minSingVal)) minSingVal = wMin
      W = MERGE(W, 0.0_rk, W > wMin)
      IF (PRESENT(tS)) tS = W
      IF (PRESENT(cn)) cn = MAXVAL(W)/MINVAL(W)
   !
   ! Get x using back substitution
      IF (PRESENT(b) .and. PRESENT(x)) THEN
         CALL svdBackSub (m=m, n=n, U=U, W=W, V=V, b=b, x=x, error=error, message=message)
      END IF
   !
   ! Get inverse
      IF (PRESENT(Ainv)) THEN
         IF (m == n) THEN
            DO i = 1, n
               DO j = 1, i
                  DO k = 1, n
                     IF (W(k) > 0.0_rk) THEN
                        tmp(k) = V(i,k)/W(k)
                     ELSE
                        tmp(k) = 0.0_rk
                     END IF
                  END DO
                  Ainv(i,j) = DOT_PRODUCT(tmp(1:n),U(j,1:n))
                  Ainv(j,i) = Ainv(i,j)
               END DO
            END DO
         ELSE
            error = 1; message = 'cannot get inverse for non-square matrix A'; GO TO 11
         END IF
      END IF
   !
   ! Clean up
11    DEALLOCATE (U, SD, W, V, tmp)
END SUBROUTINE svdSolve
!
!
  SUBROUTINE svdDecomp(a, u, s, v, nite)
    ! Singular value decomposition
    IMPLICIT NONE
    REAL(rk), DIMENSION(:,:), INTENT(IN) :: a
    REAL(rk), DIMENSION(SIZE(a,1), SIZE(a,2)), INTENT(OUT) :: u
    REAL(rk), DIMENSION(SIZE(a,2), SIZE(a,2)), INTENT(OUT) :: s, v
    INTEGER(ik),INTENT(OUT) :: nite
    REAL(rk), DIMENSION(SIZE(a,1), SIZE(a,2)) :: q1
    REAL(rk), DIMENSION(SIZE(a,1), SIZE(a,1)) :: u1
    REAL(rk), DIMENSION(SIZE(a,2), SIZE(a,2)) :: q, e 
    REAL(rk), DIMENSION(SIZE(a,2)) :: f
    REAL(rk)::err
    INTEGER :: n    
    ! init u,v,u1
    u = 0.0_rk
    v = 0.0_rk
    u1 = 0.0_rk
    FORALL (n=1:SIZE(a,1)) u1(n,n) = 1.0_rk
    FORALL (n=1:SIZE(a,2)) v(n,n) = 1.0_rk
    ! initial state:
    CALL Qr(a, q1, s)
    U = MATMUL(u1, q1)
    CALL Qr(TRANSPOSE(s), q, s)
    v = matmul(v, q)
    ! iterate while converged:
    nite = 1
    DO
       CALL Qr(TRANSPOSE(s), q, s)
       u = MATMUL(u, q)
       CALL Qr(TRANSPOSE(s), q, s)
       v = MATMUL(v, q)
       ! check the error:
       e = Triu(s)
       f = Diag(s)
       err = Norm(RESHAPE(e, [SIZE(e,1)*SIZE(e,2)])) / Norm(f)
       nite = nite + 1
       IF (err < EPS) EXIT
    END DO
  END SUBROUTINE svdDecomp

  FUNCTION Norm(x) RESULT(valu)
    ! L2-norm
    IMPLICIT NONE
    REAL(rk), DIMENSION(:), INTENT(in) :: x
    REAL(rk) :: valu
    valu = SQRT(SUM(x**2))
  END FUNCTION Norm

  FUNCTION Diag_ele(A) RESULT(v)
    ! Diagonal elements
    IMPLICIT NONE
    REAL(rk), DIMENSION(:,:), INTENT(in) :: a
    REAL(rk), DIMENSION(:), ALLOCATABLE :: v
    INTEGER :: i, n
    n = MINVAL([SIZE(a,1), SIZE(a,2)]) 
    IF (ALLOCATED(v)) DEALLOCATE(v);ALLOCATE(v(n)) 
    FORALL(i=1:n) v(i) = a(i,i)
  END FUNCTION Diag_ele

  FUNCTION Diag_mat(v) RESULT(a)
    ! Diagonal matrix
    IMPLICIT NONE
    REAL(rk), DIMENSION(:), INTENT(in) :: v
    REAL(rk), DIMENSION(SIZE(v), SIZE(v)) :: a
    INTEGER(ik) :: i
    a = 0.0_rk
    FORALL(i=1:SIZE(v)) a(i,i) = v(i)
  END FUNCTION Diag_mat

  FUNCTION Triu(a) RESULT(au)
    ! Upper triangular part
    IMPLICIT NONE
    REAL(rk), DIMENSION(:,:), INTENT(in) :: a
    REAL(rk), DIMENSION(SIZE(a,1), SIZE(a,2)) :: au
    INTEGER :: n, m, i, j
    au = 0.0_rk
    m = SIZE(a,1)
    n = SIZE(a,2)
    DO i = 1, m
       DO j = i+1, n
          IF (i+1 <= n) au(i,j) = a(i,j)
       END DO
    END DO
  END FUNCTION Triu
  
  SUBROUTINE Qr(a,q,r)
    ! Modified Gram-Schmidt process
    IMPLICIT NONE
    REAL(rk), DIMENSION(:,:), INTENT(in) :: a
    REAL(rk), DIMENSION(SIZE(a,1), SIZE(a,2)), INTENT(out) :: q
    REAL(rk), DIMENSION(SIZE(a,2), SIZE(a,2)), INTENT(out) :: r
    REAL(rk), DIMENSION(SIZE(a,1), SIZE(a,2)) :: a0
    integer :: k,n
    n = size(a,2)
    q = 0.0_rk
    r = 0.0_rk
    a0 = a
    DO k = 1, n
       r(k,k) = Norm(a0(:,k))
       q(:,k) = a0(:,k) / r(k,k)
       r(k,k+1:n) = MATMUL(q(:,k), a0(:,k+1:n))
       a0(:,k+1:n) = a0(:,k+1:n) - MATMUL(q(:,k:k), r(k:k,k+1:n))
    END DO
  END SUBROUTINE Qr

SUBROUTINE svdBackSub (m, n, U, W, V, b, x, error, message)
   ! Solves Ax=b using SVD back substitution
   ! Singular value decomposition of A(m,n) = U(m,n) * W(n) *Vtranspose (n,n)
   INTEGER(ik), INTENT(in) :: m, n
   REAL(rk), INTENT(in) :: U(:,:), W(:), V(:,:), b(:)
   REAL(rk), INTENT(out) :: x(:)
   INTEGER(ik), INTENT(inout) :: error
   CHARACTER(*),INTENT(inout) :: message
   INTEGER(ik) :: j
   REAL(rk), ALLOCATABLE :: tmp(:)
   !----
   ! Check consistency of dimension
!      error = 0; message = 'ok'

      IF (SIZE(U,1) /= m) THEN; message = 'm not same as assumed size in U(m,n)' ; GO TO 111; END IF
      IF (SIZE(U,2) /= n) THEN; message = 'n not same as assumed size in U(m,n)' ; GO TO 111; END IF
      IF (SIZE(W) /= n)   THEN; message = 'n not same as assumed size in W(n)'   ; GO TO 111; END IF
      IF (SIZE(V,1) /= n) THEN; message = 'n not same as assumed size in V(n,n1)'; GO TO 111; END IF
      IF (SIZE(V,2) /= n) THEN; message = 'n not same as assumed size in V(n1,n)'; GO TO 111; END IF
      IF (SIZE(b) /= n)   THEN; message = 'n not same as assumed size in b(n)'   ; GO TO 111; END IF
      IF (SIZE(x) /= n)   THEN; message = 'n not same as assumed size in x(n)'   ; GO TO 111; END IF
      IF (ALLOCATED(tmp)) DEALLOCATE(tmp); ALLOCATE(tmp(n))
   !
   ! Perform back substitution
      DO j = 1, n
         IF (W(j) /= 0.0_rk) THEN
            tmp(j) = DOT_PRODUCT(U(1:m,j),b(1:m))/W(j)
         ELSE
            tmp(j) = 0.0_rk
         END IF
      END DO
      DO j = 1, n
         x(j) = DOT_PRODUCT(V(j,1:n),tmp(1:n))
      END DO
      DEALLOCATE(tmp)
      RETURN
111     error = 1
END SUBROUTINE svdBackSub

END MODULE rgnMod