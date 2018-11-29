!******************************************************************
!
! Purpose: Evolves a single step of the full HYMOD model.
! ---
! Programmer: Youwei Qin and Dmitri Kavetski
! Created: Nov 2013 AD, Uni of Adelaide, Australia.
! Last modified: July 2018 AD, Hohai University, China
! Copyright, Youwei Qin, Dmitri Kavetski, 2018-2023. All rights reserved.
! ---
! Performance
!
! IN:   rain, pet, par, f
!
! OUT:  f
! ---
! Notes
!   Any issues or bugs, please contact the first author(Email:youwei.qin@uon.edu.au)
!
MODULE HyMod_Mod
   IMPLICIT NONE
   SAVE
   PRIVATE 
   PUBLIC :: Hymod
   INTEGER, PARAMETER :: mik=8      ! SELECTED_INT_KIND(9)   ! Integer kind
   INTEGER, PARAMETER :: mrk=8      ! Real kind
   INTEGER(mik),PARAMETER::isoil=1,islo=2,iqk1=3,iqk2=4,iqk3=5
      
   CONTAINS
   ! Daily HYMOD model
   SUBROUTINE Hymod (precip, pet, S, &
                     Smax, b, alpha, Ks, Kq, Qs, Qq, Q, err, message)
   IMPLICIT NONE
   !dummy parameters
   REAL(mrk),INTENT(inout)::S(:)                ! states
   REAL(mrk),INTENT(in):: precip, pet           ! forcings
   REAL(mrk),INTENT(in)::Smax,b,alpha,Ks,Kq     ! parameters
   REAL(mrk),INTENT(out)::Qs,Qq,Q               ! responses
   INTEGER(mik),INTENT(inout)::err
   CHARACTER(*),INTENT(inout)::message
   !local parameters
   CHARACTER(*),PARAMETER::procnam="evolveHYMOD"
   REAL(mrk)::Uo,Uq,Us,Ue
   ! Start procedure here
    CALL checkFeasHYMOD(S,Smax,err,message)
    IF(err/=0)THEN
       message="f-"//procnam//"/&"//message; RETURN
    ENDIF
   !The first part is to calculate the sum of direct flow and saturated flow,U0
   CALL PDM(S(isoil),Smax,b,precip,Uo,err,message)
   IF(err/=0)THEN
      message="f-"//procnam//"/&"//message; RETURN
   ENDIF
   Ue=pet                                       !     - ET at max (eg, Wagener HESS)
   Ue=min(Ue,S(isoil))
   S(isoil)=S(isoil)-Ue                         !       final soil storage
   Uq=alpha*Uo; Us=Uo-Uq                        ! *** quick vs slow flow routing
   call nashCascade(S(iqk1:),k=Kq,Vin=Uq,Vout=Qq)      ! quick flow (Cascade Of The Nash)
   call nashCascade(S(islo:islo),k=Ks,Vin=Us,Vout=Qs)  ! slow flow (single linear reservoir)
   Q=Qs+Qq 
   END SUBROUTINE Hymod
   
   PURE SUBROUTINE PDM(S,Smax,b,P,Q,err,message)
   ! Purpose: Single time step of the Probability Distribution Model (PDM).
   !          (finite non-leaky tank - soft threshold)
   ! ---
   ! Programmer: Dmitri Kavetski, Augusta-6 Glenelg
   ! ---
   ! Motivated by Matlab code 'ModPdm01.m' of Hoshin V. Gupta (9/11/2005)
   ! ---
   ! INPUTS
   !   S       = initial storage
   !   Smax    = maximum storage parameter
   !   b       = exponent parameter
   !   P       = precipitation for the time-step
   ! OUTPUTS
   !   S       = final storage
   !   Q       = model outflow (spillage + soft)
   ! ---
   ! Notes
   ! 1. Smax=Cmax/(1+Bexp)=Cmax*Bconst and Bconst=1/(1+Bexp)
   ! 2. Uses S as the primary state variable (faster and cleaner)
   ! 3. HVG used b=log(1-(Bpar/2))/log(0.5) to convert scaled B [0,2] to unscaled b [0,inf)
   ! ---
   USE constantsMod,ONLY:zero,half,one
   IMPLICIT NONE
   ! dummies
   REAL(mrk),INTENT(inout)::S
   REAL(mrk),INTENT(in)::Smax,b,P
   REAL(mrk),INTENT(out)::Q
   INTEGER(mik),INTENT(inout)::err
   CHARACTER(*),INTENT(inout)::message
   ! locals
   REAL(mrk),PARAMETER::Smin=zero,tolMBE=1.e-14_mrk
   REAL(mrk)::Sbeg,frac,OV1,OV2,Pinf,eMB,net
   REAL(mrk)::b1,ob1,Hmax,Hbeg,Hend
   CHARACTER(*),PARAMETER::procnam="PDM"
   ! Start procedure here
   IF(S>Smax)THEN
     message="f-"//procnam//"badIni[S>Smax]"
     err=-10; RETURN
   ELSEIF(S<Smin)THEN
     message="f-"//procnam//"badIni[S<Smin]"
     err=-10; RETURN
   ELSEIF(b<zero)THEN
     message="f-"//procnam//"illegal[b<zero]"
     err=-10; RETURN
   ENDIF
   Sbeg=S; OV2=zero
  ! Moore's description, D'z implementation based on Gupta's Matlab
   b1=b+one; ob1=one/b1; Hmax=Smax*b1  ! max height
   frac=one-(one-S/Smax)**ob1          ! active fraction
   Hbeg=Hmax*frac                      ! initial height
   OV2=MAX(Hbeg+P-Hmax,zero)           ! spillage
   Pinf=P-OV2                          ! P that does not spill
   Hend=MIN(Hbeg+Pinf,Hmax)            ! final height
   S=Smax*(one-(one-Hend/Hmax)**b1)    ! final storage
   OV1=MAX(Sbeg+Pinf-S,zero)           ! soft outflow
   Q=OV1+OV2                           ! total outflow
   IF(S>Smax)THEN
     message="f-"//procnam//"bug?[Send>Smax]"
     err=-10; RETURN
   ELSEIF(S<Smin)THEN
     message="f-"//procnam//"bug?[Send<Smin]"
     err=-10; RETURN
   ENDIF
   net=P-Q
   eMB=Sbeg+net-S                      ! S-balance discrEpency
   IF(ABS(eMB)>tolMBE*MAX(Smax,P,Q))THEN
     message="f-"//procnam//"/MBE[eMB>tol]"
     err=-10; RETURN
   ENDIF
   ! End procedure here
   ENDSUBROUTINE PDM
    
   !----------------------------------------------------
   PURE SUBROUTINE nashCascade(S,k,Vin,Vout)
   ! Purpose: Executes a single time step of the Nash Cascade Model.
   !          (Series of infinite linear leaky tanks)
   ! ---
   ! Programmer: Dmitri Kavetski, Augusta-6 Glenelg
   ! ---
   ! Motivated by Matlab code 'ModNash.m' of Hoshin V. Gupta (9/18/2005)
   ! ---
   ! INPUTS
   !   S     = Initial storage in individual tanks
   !   k     = Leakage rate of the tanks
   !   Vin   = Input to first tank in the series
   ! OUTPUTS
   !   Vout  = Outflow from last tank in series
   !   S     = Final storage in individual tanks
   ! ---
   USE constantsMod,ONLY:zero,half,one
   IMPLICIT NONE
   ! dummies
   REAL(mrk),INTENT(inout)::S(:)
   REAL(mrk),INTENT(in)::k
   REAL(mrk),INTENT(in)::Vin
   REAL(mrk),INTENT(out)::Vout
   ! locals
   INTEGER(mik)::i,n
   REAL(mrk)::q(0:SIZE(S))
   ! Start procedure here
   n=SIZE(S); q(0)=Vin
   DO i=1,n
     q(i)=k*S(i); S(i)=S(i)+q(i-1)-q(i)
   ENDDO
   Vout=q(n)
   ! End procedure here
   ENDSUBROUTINE nashCascade
   
   !----------------------------------------------------
    PURE SUBROUTINE checkFeasHYMOD(S,Smax,err,message)
    ! Purpose: Checks feasibility of HYMOD states and parameters.
    USE constantsMod,ONLY:zero
    IMPLICIT NONE
    ! dummies
    REAL(mrk),INTENT(in)::S(:),Smax
    INTEGER(mik),INTENT(inout)::err
    CHARACTER(*),INTENT(inout)::message
    ! locals
    CHARACTER(*),PARAMETER::procnam="checkFeasHYMOD"
    CHARACTER(*),PARAMETER::fmt1='(a,es15.8,a,es15.8,a)'
    CHARACTER(*),PARAMETER::fmt2='(a,i0,a,es15.8,a)'
    ! Start procedure here
    IF(S(isoil)>Smax)THEN
      WRITE(message,fmt1)"f-"//procnam//"/badIni&
        &[S(soil)>Smax;(",S(isoil),")>(",Smax,")]"
      err=-10; RETURN
    ELSEIF(ANY(S<zero))THEN
      WRITE(message,fmt2)"f-"//procnam//"/badIni[S)<0.0]"
      err=-10; RETURN
    ENDIF
    ! End procedure here
    ENDSUBROUTINE checkFeasHYMOD
   !----------------------------------------------------
END MODULE HyMod_Mod