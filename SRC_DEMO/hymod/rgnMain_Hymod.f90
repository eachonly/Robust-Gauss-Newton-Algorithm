PROGRAM testRGN
! Purpose: Calibrate HYMOD parameters with Robust Gauss-Newton Algorithm (RGN)
! ---
! Programmer: Youwei Qin, Dmitri Kavetski, Michael Leonard
! Created: July 2018 AD, Hohai University, China.
! Last modified: October 2018 AD, Hohai University, China
! Copyright, Youwei Qin, Dmitri Kavetski, Michael Leonard, 2018-2023. All rights reserved.
! ---
! This is the demo for calibrating HYMOD with RGN
! The core of RGN is recorded in rgn.f90; the core of HYMOD is recoded in hymod.f90
! The data exchange between RGN and HYMOD is through "objFunc", where the HYMOD is called,
! and the sum of least squares objective function value is evaluated and returned to RGN subroutine.
! The public variables were shared through subroutine "constantsMod.f90"

! Firstly the information about parameters are loaded, such as parameter name, initial value, lower and upper limit;
! Then the input data was loaded, which include the rainfall, ET, and runoff
! Finally, the initial states for storage are loaded.
!******************************************************************
   USE rgnMod
   USE constantsMod, ONLY: ik, rk
   USE hydroDataMod, ONLY: nWarmUp, nData, nPar, nState, rain, pet, obsQ, parName, stateName, stateVal
   IMPLICIT NONE
!  Constants
   CHARACTER(*),PARAMETER::procnam="testMain"
!   INTEGER(ik):: nWarmup,nData, nPar,nState
   INTEGER(ik)::i, status
   REAL(rk),ALLOCATABLE :: x0(:), xLo(:), xHi(:), x(:)
!   REAL(rk),ALLOCATABLE::rain(:),pet(:),obsQ(:)
   TYPE (rgnConvType) :: cnv
   INTEGER(ik) :: error
   CHARACTER(256) :: message
!   CHARACTER(128),ALLOCATABLE::parName(:)
   CHARACTER(256)::dataFileName
   TYPE (rgnInfoType) :: info
   EXTERNAL objFunc
   CHARACTER(20) :: dfm1
   !Write out the message what is running
   WRITE(*,'(a)') "Calibrating Hymod with RGN, approximate running time 10-40 seconds"
   WRITE(*,*)
   !----
!   error=0                                              ! Initialize error flag
!   message=""                                           ! Initialize message
   !Part 1: load files for HYMOD
   !Get the basic information of the model with file 'inputData.txt',
   !which includes the parameter information and initial status of states
   OPEN (UNIT=1,FILE='inputData.txt',STATUS='old',ACTION='READ',IOSTAT=status)
   IF (status /= 0)THEN
       WRITE(*,*) "Error in opening inputData.txt"
       READ(*,*);STOP
   END IF    
   READ(1,*) npar    !number of parameters 
   
   !Allocate the public data
   ALLOCATE(xLo(nPar),xHi(nPar),x0(npar),x(nPar),STAT=status)
   IF (status /= 0) THEN
       WRITE(*,*) "Error in allocating lowPar, highPar, and x0" 
       READ(*,*);STOP
   END IF
   !Initialize the parameters
   xLo=0.0_rk;xHi=0.0_rk
   !Allocate the Parameter-related array. like the name,Lowpar,highPar
   ALLOCATE(parName(nPar),STAT=status)
   IF (status /= 0) THEN
       WRITE(*,*) "Error in allocating name"
       READ(*,*);STOP
   END IF
   !Read the parameters Line by Line
   DO i=1,npar
       READ(1,*) parName(i),xLo(i),xHi(i),x0(i)
   END DO
   !Continue reading the File, the dataFileName is the name of the input file, includes rainfall, runoff, evaporation
   READ(1,*) dataFileName  
   !Read the number of days used for warmup strategy
   READ(1,*) nWarmUp
   CLOSE(UNIT=1)
   
   !Get the basic information of the model states with file 'states.txt',
   OPEN (UNIT=1,FILE='states.txt',STATUS='old',ACTION='READ',IOSTAT=status)
   IF (status /= 0)THEN
       WRITE(*,*) "Error in opening states.txt"
       READ(*,*);STOP
   END IF
   READ(1,*) nState    !number of parameters 
         
   !Allocate the Parameter-related array, such as name
   ALLOCATE(stateName(nState),stateVal(nState),STAT=status)
   IF (status /= 0) THEN
       WRITE(*,*) "Error in allocating stateName"
       READ(*,*);STOP
   END IF
   !Read the parameters Line by Line
   DO i=1,nState
       READ(1,*) stateName(i),stateVal(i)
   END DO
   CLOSE(UNIT=1)
   
   !Get the input and observed response data
   OPEN (UNIT=1,FILE=dataFileName,STATUS='old',ACTION='READ',IOSTAT=status)
   IF (status /= 0) THEN
       WRITE(*,*) "Error in opening datafile"
   END IF
   ! Start to read the file. The first line is the head, need to skip
   READ(1,*)
   nData=0
   !Read the file to get the number of data in the file
   DO
       READ(1,*,IOSTAT=status)
       IF (status /= 0) EXIT
       nData=nData+1
   END DO
!Allocate the array to store the Data
   ALLOCATE(rain(nData),pet(nData),obsQ(nData),STAT=status)
   IF (status /= 0) THEN
       WRITE(*,*) "Error in allocating rainfall-runoff observations"
       READ(*,*);STOP
   END IF
!Read the contents in the File
   REWIND(UNIT=1);READ(1,*)
   DO i=1,nData
       READ(1,*) rain(i),pet(i),obsQ(i)
   END DO    
  CLOSE(UNIT=1)
!Give the format
WRITE(dfm1,'(a,i4,a)')     '(a,', nPar,'g15.7)'
!Part 2: Run RGN
!Initialize the RGN default settings
  CALL setDefaultRgnConvergeSettings (cnvSet=cnv, dump=10_ik, fail=0_ik)
!Call RGN optimization algorithms
! key input parameters: p is the number of parameters to be optimized
!                       n is the number of residuals
  CALL rgn (objFunc=objFunc, p=nPar, n=nData-nWarmUp+1, x0=x0, xLo=xlo, xHi=Xhi, cnv=cnv, x=x, info=info, error=error, message=message)
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
  DEALLOCATE(xLo,xHi,x0,x)
  DEALLOCATE(parName,stateName,stateVal)
  DEALLOCATE(rain,pet,obsQ)
END PROGRAM testRGN

SUBROUTINE objFunc (nPar, nSim, x, r, f, timeFunc, error, message)
   USE constantsMod, ONLY: ik, rk
   USE Hymod_Mod
   USE hydroDataMod, ONLY: nWarmUp,rain,pet,obsQ,stateVal
   IMPLICIT NONE
   INTEGER(ik), INTENT(in) :: nPar                  ! Number of parameters
   INTEGER(ik), INTENT(in) :: nSim                  ! Number of data
   REAL(rk), INTENT(in)    :: x(:)                  ! Parameters
   REAL(rk), INTENT(out)   :: r(:)                  ! Residuals
   REAL(rk), INTENT(out)   :: f                     ! Objective function value = Sum of squared residuals
   REAL(rk),INTENT(out)    :: timeFunc              ! Time for evaluate the objective function
   INTEGER(ik), INTENT(inout)    :: error
   CHARACTER(*),INTENT(inout) :: message
   INTEGER(ik) :: i

   !time for evaluating
   REAL(rk)::timeObj(2)
   !Hymod Parameters and states
   REAL(rk):: S(5)                      ! soil moister states
   REAL(rk):: Smax,b,alpha,Ks,Kq        ! parameters
   REAL(rk):: Qs,Qq,Q                   ! responses
   INTEGER(ik)::status
   LOGICAL:: flexS
   CHARACTER(*),PARAMETER::procnam="ObjFunc"
   !---
   !
   CALL CPU_TIME (timeObj(1))
   flexS=.true.         ! Allow fix of Smax
   !Assign parameters
   Smax=x(1)            ! Maximum storage capacity
   b=x(2)               ! Degree of spatial variability of the soil moisture capacity
   alpha=x(3)           ! Factor distributing the flow between slow and quick release reservoirs
   Ks=x(4)              ! Residence time of the slow release reservoir
   Kq=x(5)              ! Residence time of the quick release reservoir
   
   ! Initialize surfacewater storages and baseflowstorage, and the initial storage of C1,C2,C3
   S(1)=stateVal(1);S(2)=stateVal(2);S(3)=stateVal(3);S(4)=stateVal(4); S(5)=stateVal(5)
   Qs=stateVal(6);Qq=stateVal(7);Q=stateVal(8)
   
   ! check feas of state wrt pars unless flexi-state
   IF(S(1)>Smax.and..not.flexS)THEN
        message="f-"//procnam//"/Soil moisture exceeds"
        error=-10;RETURN
   ENDIF
   ! * Allows convenient initialization and adjustment of states
   if(flexS)then
       if(S(1)>Smax)S(1)=Smax
   endif
    
   !Warnup model so that memory of arbitrary initial conditions is forgotten
   OPEN (UNIT=1,FILE='dischargenew.txt') 
   DO i = 1, nWarmUp-1
       CALL HyMod(precip=rain(i), pet=pet(i), S=S,    &
       Smax=Smax, b=b, alpha=alpha, Ks=Ks, Kq=Kq,   &
       Qs=Qs, Qq=Qq, Q=Q, err=error, message=message)
       IF (error /=0) THEN
           message="f-"//procnam//"/&"//message; RETURN
       END IF
       WRITE(1,*)  Q
   END DO
!
! Compute R^2  for remainder of observed record
   f = 0.0_rk
   DO i = nWarmUp, nWarmUp+nSim-1
       CALL HyMod(precip=rain(i), pet=pet(i), S=S,    &
       Smax=Smax, b=b, alpha=alpha, Ks=Ks, Kq=Kq,   &
       Qs=Qs, Qq=Qq, Q=Q, err=error, message=message)
       WRITE(1,*)  Q
       IF (error /=0) THEN
           message="f-"//procnam//"/&"//message; RETURN
       END IF       
       r(i-nWarmUp+1)=obsQ(i)-Q
       f=f+r(i-nWarmUp+1)**2
   END DO  
  f = f/2.0_rk
  WRITE(1,*) f
  CLOSE(UNIT=1)
  CALL CPU_TIME (timeObj(2))
  timeFunc=timeObj(2)-timeObj(1)
END SUBROUTINE objFunc
