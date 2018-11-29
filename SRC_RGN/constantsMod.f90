!******************************************************************
!
! Purpose: modules for some most widely used constants

! Author: George Kuczera
! Email: george.kuczera@newcastle.edu.au
! School of Engineering
! University of Newcastle
! Callaghan, NSW, 2308, Australia
!
! Copyright, George Kuczera, 2005-2010. All rights reserved.
!
MODULE constantsMod
   IMPLICIT NONE
   SAVE
   PUBLIC
   INTEGER, PARAMETER :: ik4=4   ! SELECTED_INT_KIND(9)   ! Integer kind
   INTEGER, PARAMETER :: ik=8   ! SELECTED_INT_KIND(9)   ! Integer kind
   INTEGER, PARAMETER :: rk4=4   ! Real*4 kind
   INTEGER, PARAMETER :: rk=8   ! Real kind
   !
   ! Machine constants
   REAL(rk), PARAMETER :: hugeRe = HUGE(1.0_rk)           ! largest real on machine
   REAL(rk), PARAMETER :: tinyRe = TINY(1.0_rk)           ! smallest real on machine
   REAL(rk), PARAMETER :: epsRe = EPSILON(1.0_rk)         ! smallest real epsilon on machine
   INTEGER(ik), PARAMETER :: minExp = MINEXPONENT(1.0_rk) ! smallest exponent
   INTEGER(ik), PARAMETER :: maxExp = MAXEXPONENT(1.0_rk) ! biggest exponent
   INTEGER(ik), PARAMETER :: hugeInt = HUGE(1_ik)            ! largest integer on machine
   !
   ! Other constants
   REAL(rk), PARAMETER :: zero=0.0_rk, half=0.5_rk, one=1.0_rk, two=2.0_rk, four=4.0_rk
   REAL(rk), PARAMETER :: pi = 3.141592653589793238462643383279502884197_rk
   REAL(rk), PARAMETER :: twoPi = 6.283185307179586476925286766559005768394_rk
   !
   ! Script constants and structures
   INTEGER(ik), PARAMETER :: maxScriptLines=700, nCharScriptLine=60, scriptSize=maxScriptLines*nCharScriptLine, codeSize=3*scriptSize, maxnScriptProc=40
   TYPE scriptProcType
      INTEGER(ik) :: nCode
      CHARACTER(20) :: name, nameUC
      CHARACTER(50) :: description
      CHARACTER(scriptSize) :: formula
      INTEGER(ik4) :: code(0:codeSize)
   END TYPE scriptProcType
   
END MODULE constantsMod

