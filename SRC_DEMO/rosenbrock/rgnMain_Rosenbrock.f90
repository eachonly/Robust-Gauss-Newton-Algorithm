PROGRAM testRGN
! Purpose: Calibrate 2D Rosenbrock function with Robust Gauss-Newton Algorithm (RGN)
! ---
! Programmer: Youwei Qin, Dmitri Kavetski, Michael Leonard
! Created: July 2018 AD, Hohai University, China.
! Last modified: October 2018 AD, Hohai University, China
! Copyright, Youwei Qin, Dmitri Kavetski, Michael Leonard, 2018-2023. All rights reserved.
! ---
! This is the demo for calibrating Rosenbrock function with RGN
! The core of RGN is recorded in rgn.f90
! The data exchange between RGN and Rosenbrock function is through "objFunc"
! and the sum of least squares objective function value is evaluated and returned to RGN subroutine.
! The public variables were shared through subroutine "constantsMod.f90"
!******************************************************************
   USE rgnMod
   USE constantsMod
   IMPLICIT NONE
   INTEGER(ik), PARAMETER :: p=2, n=2
   REAL(rk) :: x0(p), xLo(p), xHi(p), x(p)
   TYPE (rgnConvType) :: cnv
   INTEGER(ik) :: error
   CHARACTER(256) :: message
   TYPE (rgnInfoType) :: info
   EXTERNAL objFunc
   CHARACTER(20) :: dfm1
   !----
   !Write out the message what is running
   WRITE(*,'(a)') " Calibrating Rosenbrock with RGN, approximate running time 1-2 seconds"
   WRITE(*,*)
   !
!   error=0                                              ! Initialize error flag
!   message=""                                           ! Initialize message
   x0 = [-1.0_rk, 0.0_rk]                               ! Start point of the search, with the optimum at [1.0 1.0]
   xLo = [-1.5_rk, -1.0_rk]                             ! Low bound
   xhi = [ 1.5_rk,  3.0_rk]                             ! Upper bound
   !Give the format
   WRITE(dfm1,'(a,i4,a)')     '(a,', 2,'g15.7)'
   CALL setDefaultRgnConvergeSettings (cnvSet=cnv, dump=10_ik, fail=0_ik)
   !
   ! key input parameters: p is the number of parameters to be optimized
   !                       n is the number of residuals
   CALL rgn (objFunc=objFunc, p=2_ik, n=2_ik, x0=x0, xLo=xlo, xHi=Xhi, cnv=cnv, x=x, info=info, error=error, message=message)
   IF(error /= 0)then
     WRITE(*,*) message
     READ(*,*)
   END IF
   WRITE(*,dfm1)                "Best parameter set:     ", x
   WRITE(*,'(a,g15.7)')         "Best objfunc value:     ", info%f
   WRITE(*,'(a,4x,i0)')         "Number of objfunc calls:", info%nEval
   WRITE(*,'(a,4x,i0)')         "Total iteration:         ", info%nIter
   WRITE(*,'(a,4x,i0)')         "Termination flag:        ", info%termFlag
   WRITE(*,'(a,g15.7)')         "CPU time:                ",info%cpuTime
END PROGRAM testRGN

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
   CHARACTER(*),INTENT(inout) :: message
   INTEGER(ik) :: i
   !---
   !
   !time for evaluating
   REAL(rk)::timeObj(2)
   CALL CPU_TIME (timeObj(1))
   f = 0.0_rk
   r(1) = 1.0_rk-x(1)
   r(2)=10.0_rk*(x(2)-x(1)**2)                      ! Compute residual
   f = f + r(1)**2+r(2)**2                          ! Calculate objective function
   f = f/2.0_rk
   CALL CPU_TIME (timeObj(2))
   timeFunc=timeObj(2)-timeObj(1)
END SUBROUTINE objFunc
