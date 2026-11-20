!===============================================================================!
MODULE MOD_TimeDiscretization
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE TimeDiscretization
  MODULE PROCEDURE TimeDiscretization
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: TimeDiscretization
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretization()
!-------------------------------------------------------------------------------!
USE MOD_Output,             ONLY: WriteMeshToDisk
USE MOD_Output,             ONLY: WriteSolutionToDisk
USE MOD_Equation,           ONLY: TimeStep
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: tEnd
USE MOD_FiniteVolume2D_vars,ONLY: dt_Analyze
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: tGlobal 
USE MOD_FiniteVolume2D_vars,ONLY: timescheme
USE MOD_FiniteVolume2D_vars,ONLY: maxTimeSteps
USE MOD_FiniteVolume2D_vars,ONLY: computationalTime
USE MOD_FiniteVolume2D_vars,ONLY: U, nElemsX, nElemsY, global_min
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: NRelaxedTimesteps
USE MOD_FiniteVolume2D_vars,ONLY: relaxed_dt
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: U
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
#endif

#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
#endif

#ifdef ACTIVEFLUX
!*Postprocessing AF
USE MOD_FiniteVolume2D     ,ONLY: AF_PostProcessing_Subroutine
#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro
#endif
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D     ,ONLY: Overwrite_WC_from_U
#ifdef PATHCONSERVATIVESHOCKDETECTION
USE MOD_FiniteVolume2D     ,ONLY: Shock_Detector_Based_On_Path_Conservative
USE MOD_FiniteVolume2D_vars,ONLY: IsSomeoneFlagged_PCSD
#endif
#endif
USE MOD_FiniteVolume2D     ,ONLY: ComputeTransitionMatricesConservedInput
USE MOD_FiniteVolume2D     ,ONLY: Impose_Symmetric_Update_Y_Axis
USE MOD_Equation           ,ONLY: BoundaryConditions
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL :: t, tAnalyze, dt_min, time_start, time_end, local_min
INTEGER :: iStep, iVar
!-------------------------------------------------------------------------------!

t = 0.0
dt_Analyze = tEnd/REAL(nOutputFiles)
tAnalyze   = t + dt_Analyze
IF (WhichOutput .EQ. 1 .OR. WhichOutput .EQ. 3) THEN
  CALL WriteMeshToDisk()
END IF
CALL WriteSolutionToDisk(t)
IF (t .EQ. tEnd) THEN
  RETURN
END IF

call cpu_time(time_start)
global_min = 100.

iStep = 0
DO WHILE ( iStep .LT. maxTimeSteps)

#ifdef CENTEREDPRIMITIVE
#ifdef PATHCONSERVATIVESHOCKDETECTION
  !*Before overwriting
  CALL Shock_Detector_Based_On_Path_Conservative()
#endif
#ifndef DONOTOVERWRITEPRIMITIVE
  CALL Overwrite_WC_from_U()
#endif
#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
  CALL Compute_min_p_max_ro()
#endif
#endif

#if(1==0)
PRINT*, "I DO NOT WANT TO RUN IN THIS MODE (TIME STEP COMPUTATION AS ALEX)"
PRINT*, "Stopping in timediscretization"
STOP
!*NB: This is needed if I compute dt based on the reconstructed values
CALL BoundaryConditions(t)
CALL ComputeTransitionMatricesConservedInput(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 
#endif

  iStep = iStep + 1
  dt_min = TimeStep(t)

  IF (iStep .LT. NRelaxedTimesteps+1) THEN
    PRINT*, "We are relaxing iteration", iStep
    dt_min=relaxed_dt
    PRINT*, "...at dt", dt_min
  END IF


  dt = MIN(MIN(dt_min,tAnalyze-t),MIN(dt_min,tEnd-t))

  !*=========================================================
  PRINT*, "Iteration", iStep, "with dt", dt, "time", t
  ! IF (iStep .EQ. 3) STOP
  !*=========================================================
!   DO iVar=1,nVar
!     PRINT*, "U  ", iVar, MAXVAL(U(iVar,1:nElemsX,1:nElemsY)), MINVAL(U(iVar,1:nElemsX,1:nElemsY))
! #ifdef ACTIVEFLUX
!     PRINT*, "W_X", iVar, MAXVAL(W_X(iVar,1:nElemsX+1,1:nElemsY)), MINVAL(W_X(iVar,1:nElemsX+1,1:nElemsY))
!     PRINT*, "W_Y", iVar, MAXVAL(W_Y(iVar,1:nElemsX,1:nElemsY+1)), MINVAL(W_Y(iVar,1:nElemsX,1:nElemsY+1))
! #endif
! #ifdef CENTEREDPRIMITIVE
!     PRINT*, "WC ", iVar, MAXVAL(WC(iVar,1:nElemsX,1:nElemsY)), MINVAL(WC(iVar,1:nElemsX,1:nElemsY))
! #endif
!   END DO
!   PRINT*
  !*=========================================================

#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
! PRINT*, "Iteration", iStep, "with dt", dt, "time", t, "Someone flagged?", IsSomeoneFlagged_PCSD
#endif

  SELECT CASE(timescheme)
#if defined(IMEX)
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
    CASE(-1)
      CALL TimeDiscretizationByIMEXEuler(t)
    CASE(-2)
      CALL TimeDiscretizationByIMEXDeC2(t)
    CASE(-3)
      CALL TimeDiscretizationByIMEXDeC3(t)
    CASE(-4)
      CALL TimeDiscretizationByIMEXDeC4(t)
    CASE(-5)
      CALL TimeDiscretizationByIMEXDeC5(t)
    CASE(-11)
      CALL TimeDiscretizationByIMEXEuler_Implicit_Update_Conserved(t)
    CASE(-12)
      CALL TimeDiscretizationByIMEXDeC2_Implicit_Update_Conserved(t)
    CASE(-13)
      CALL TimeDiscretizationByIMEXDeC3_Implicit_Update_Conserved(t)
    CASE(-14)
      CALL TimeDiscretizationByIMEXDeC4_Implicit_Update_Conserved(t)
    CASE(-15)
      CALL TimeDiscretizationByIMEXDeC5_Implicit_Update_Conserved(t)
    CASE(-22)
      CALL TimeDiscretizationByIMEXDeC2_Crank_Nicolson_Conserved(t)
    CASE(-52) !*In this version max_ro and min_p are stage dependent
      CALL IMEXDeC2_Crank_Nicolson_Conserved_rho_p_stage_dependent(t)
#ifdef ACTIVEFLUX
    CASE(-32)
      CALL TimeDiscretizationByIMEXDeC2_Implicit_Euler_Conserved(t)
#endif
#ifdef CENTEREDPRIMITIVE
    CASE(-42)
      CALL TimeDiscretizationByARS222_Crank_Nicolson_Conserved(t)
    CASE(-82)
      CALL IMEXDeC2_rho_p_stage_dependent_overwriting_stages(t)
#endif
#endif

#elif defined(IMEXMOMENTUM)
    CASE(-101)
      CALL TimeDiscretizationByIMEXMOMENTUMEuler(t)
    CASE(-102)
      CALL TimeDiscretizationByIMEXMOMENTUMDeC2(t)
    ! CASE(-103)
    !   CALL TimeDiscretizationByIMEXMOMENTUMDeC3(t)
    ! CASE(-104)
    !   CALL TimeDiscretizationByIMEXMOMENTUMDeC4(t)
    ! CASE(-105)
    !   CALL TimeDiscretizationByIMEXMOMENTUMDeC5(t)
    ! CASE(-122)
    !   CALL TimeDiscretizationByIMEXMOMENTUMDeC2_Crank_Nicolson_Conserved(t)


#else
    CASE(1)
      CALL TimeDiscretizationByForwardEuler(t)
    CASE(2)
      CALL TimeDiscretizationBySSPRK2(t)
    CASE(3)
      CALL TimeDiscretizationBySSPRK3(t)
    CASE(4)
      CALL TimeDiscretizationBySSPRK4(t)
    CASE(5)
      CALL TimeDiscretizationByRK65(t)
    CASE(12)
      CALL TimeDiscretizationByDeC2(t)
    CASE(13)
      CALL TimeDiscretizationByDeC3(t)
    CASE(14)
      CALL TimeDiscretizationByDeC4(t)
    CASE(15)
      CALL TimeDiscretizationByDeC5(t)

#ifdef PATANKAR
    CASE(21)
      CALL TimeDiscretizationByMPEuler(t)
    CASE(22)
      CALL TimeDiscretizationByMPDeC2(t)
    CASE(25)
      CALL TimeDiscretizationByMPDeC5(t)
#endif
#endif
    CASE default
      WRITE(*,*) "Time discretization not implemented"
      STOP
  END SELECT

#ifdef ACTIVEFLUX
#ifndef PRIMITIVEONLY
  CALL AF_PostProcessing_Subroutine(t+dt) !*NB: The input time is needed for the BCs
#endif
#endif

#if defined(CENTEREDPRIMITIVE) && defined(FIRSTORDERPRIMITIVE)
  CALL FirstOrderEvolutionOnlyPrimitive(t) 
#endif

      ! PRINT*, 1, MAXVAL(u(1,1:nElemsX,1:nElemsY)), MINVAL(u(1,1:nElemsX,1:nElemsY))
      ! PRINT*, 2, MAXVAL(u(2,1:nElemsX,1:nElemsY)), MINVAL(u(2,1:nElemsX,1:nElemsY))
      ! PRINT*, 3, MAXVAL(u(3,1:nElemsX,1:nElemsY)), MINVAL(u(3,1:nElemsX,1:nElemsY))
      ! PRINT*, 4, MAXVAL(u(4,1:nElemsX,1:nElemsY)), MINVAL(u(4,1:nElemsX,1:nElemsY))


  t = t + dt


  ! CALL Impose_Symmetric_Update_Y_Axis()

  ! PRINT*, "time ", t, ", time iter ", iStep, ", dt ", dt
  local_min = MINVAL(U(1,1:nElemsX,1:nElemsY))
  global_min = MIN(local_min, global_min)
  ! PRINT*, "min h", local_min, ", min at all time ", global_min
  tGlobal = t
  IF (ABS(tAnalyze-t) .LT. 1.0E-10) THEN
    CALL WriteSolutionToDisk(t)
    tAnalyze = tAnalyze + dt_Analyze
    IF (tAnalyze .GT. tEnd) THEN
      tAnalyze = tEnd
    END IF
  END IF
  IF (ABS(t-tEnd) .LT. 1.0E-10) THEN
    EXIT
  END IF
END DO

call cpu_time(time_end)
computationalTime = time_end - time_start

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretization
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByForwardEuler(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: K0
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
#endif

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: K0_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: K0_Y
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj
!-------------------------------------------------------------------------------!

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif
#ifdef ACTIVEFLUX
K0_X(1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY) !*NB: +1 in X direction
K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1) !*NB: +1 in Y direction
#endif

!--------------------!
! Forward Euler      !
!--------------------!
tStage = t 
CALL FVTimeDerivative(tStage)

U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
WC(1:nVar,1:nElemsX,1:nElemsY) = K0_WC(1:nVar,1:nElemsX,1:nElemsY) + WCt(1:nVar,1:nElemsX,1:nElemsY)*dt
#endif
#ifdef ACTIVEFLUX
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = K0_X(1:nVar,1:nElemsX+1,1:nElemsY) + Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt !*NB:+1 in X direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) + Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt !*NB:+1 in Y direction
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByForwardEuler
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(CENTEREDPRIMITIVE) && defined(FIRSTORDERPRIMITIVE)
SUBROUTINE FirstOrderEvolutionOnlyPrimitive(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj
!-------------------------------------------------------------------------------!

K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Forward Euler      !
!--------------------!
tStage = t 
CALL FVTimeDerivativePrimitiveOnly(tStage,ReconstructionInput=1)

WC(1:nVar,1:nElemsX,1:nElemsY) = K0_WC(1:nVar,1:nElemsX,1:nElemsY) + WCt(1:nVar,1:nElemsX,1:nElemsY)*dt


!-------------------------------------------------------------------------------!
END SUBROUTINE FirstOrderEvolutionOnlyPrimitive
#endif
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationBySSPRK2(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K1_WC
USE MOD_FiniteVolume2D_vars,ONLY: K2_WC
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: K0_X
USE MOD_FiniteVolume2D_vars,ONLY: K1_X
USE MOD_FiniteVolume2D_vars,ONLY: K2_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: K0_Y
USE MOD_FiniteVolume2D_vars,ONLY: K1_Y
USE MOD_FiniteVolume2D_vars,ONLY: K2_Y

USE MOD_FiniteVolume2D,     ONLY: AF_PostProcessing_Subroutine
#endif

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj
!-------------------------------------------------------------------------------!

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
K0_X(1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY) !*NB: +1 in X direction

!*Staggering in Y direction
K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1) !*NB: +1 in Y direction
#endif

!--------------------!
! First Stage        !
!--------------------!
tStage = t + 0.0*dt
CALL FVTimeDerivative(tStage)
K1(1:nVar,1:nElemsX,1:nElemsY) = &
    K0(1:nVar,1:nElemsX,1:nElemsY) &
  + Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K1(1:nVar,1:nElemsX,1:nElemsY)
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
K1_WC(1:nVar,1:nElemsX,1:nElemsY) = &
    K0_WC(1:nVar,1:nElemsX,1:nElemsY) &
  + WCt(1:nVar,1:nElemsX,1:nElemsY)*dt
WC(1:nVar,1:nElemsX,1:nElemsY)  = K1_WC(1:nVar,1:nElemsX,1:nElemsY)
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
K1_X(1:nVar,1:nElemsX+1,1:nElemsY) = &    !*NB: +1 in X direction
    K0_X(1:nVar,1:nElemsX+1,1:nElemsY) &  !*NB: +1 in X direction
  + Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt !*NB: +1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY)  = K1_X(1:nVar,1:nElemsX+1,1:nElemsY) !*NB: +1 in X direction

!*Staggering in Y direction
K1_Y(1:nVar,1:nElemsX,1:nElemsY+1) = &    !*NB: +1 in Y direction
    K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) &  !*NB: +1 in Y direction
  + Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt !*NB: +1 in Y direction 
W_Y(1:nVar,1:nElemsX,1:nElemsY+1)  = K1_Y(1:nVar,1:nElemsX,1:nElemsY+1) !*NB: +1 in Y direction
#endif



!--------------------!
! Second Stage       !
!--------------------!
tStage = t + 1.0*dt
#ifdef POSTPROCESSINGATEACHSTEP
CALL AF_PostProcessing_Subroutine(tStage)
K1(1:nVar,1:nElemsX,1:nElemsY)=U(1:nVar,1:nElemsX,1:nElemsY)
K1_X(1:nVar,1:nElemsX+1,1:nElemsY)=W_X(1:nVar,1:nElemsX+1,1:nElemsY)
K1_Y(1:nVar,1:nElemsX,1:nElemsY+1)=W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = &
    0.5*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.5*K1(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.5*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K2(1:nVar,1:nElemsX,1:nElemsY)
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
K2_WC(1:nVar,1:nElemsX,1:nElemsY) = &
    0.5*K0_WC(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.5*K1_WC(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.5*WCt(1:nVar,1:nElemsX,1:nElemsY)*dt
WC(1:nVar,1:nElemsX,1:nElemsY)  = K2_WC(1:nVar,1:nElemsX,1:nElemsY)
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
K2_X(1:nVar,1:nElemsX+1,1:nElemsY) = &
    0.5*K0_X(1:nVar,1:nElemsX+1,1:nElemsY) &
  + 0.5*K1_X(1:nVar,1:nElemsX+1,1:nElemsY) &
  + 0.5*Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt
W_X(1:nVar,1:nElemsX+1,1:nElemsY)  = K2_X(1:nVar,1:nElemsX+1,1:nElemsY)

!*Staggering in Y direction
K2_Y(1:nVar,1:nElemsX,1:nElemsY+1) = &
    0.5*K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) &
  + 0.5*K1_Y(1:nVar,1:nElemsX,1:nElemsY+1) &
  + 0.5*Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt
W_Y(1:nVar,1:nElemsX,1:nElemsY+1)  = K2_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationBySSPRK2
!===============================================================================!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationBySSPRK3(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K1_WC
USE MOD_FiniteVolume2D_vars,ONLY: K2_WC
USE MOD_FiniteVolume2D_vars,ONLY: K3_WC
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: K0_X
USE MOD_FiniteVolume2D_vars,ONLY: K1_X
USE MOD_FiniteVolume2D_vars,ONLY: K2_X
USE MOD_FiniteVolume2D_vars,ONLY: K3_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: K0_Y
USE MOD_FiniteVolume2D_vars,ONLY: K1_Y
USE MOD_FiniteVolume2D_vars,ONLY: K2_Y
USE MOD_FiniteVolume2D_vars,ONLY: K3_Y

USE MOD_FiniteVolume2D,     ONLY: AF_PostProcessing_Subroutine
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj
!-------------------------------------------------------------------------------!

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
K0_X(1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY) !*NB: +1 in X direction

!*Staggering in Y direction
K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1) !*NB: +1 in Y direction
#endif


!--------------------!
! First Stage        !
!--------------------!
tStage = t + 0.0*dt
CALL FVTimeDerivative(tStage,Optional_Input=1) !*No need for post-processing, solution in input
K1(1:nVar,1:nElemsX,1:nElemsY)= Ut(1:nVar,1:nElemsX,1:nElemsY)
U(1:nVar,1:nElemsX,1:nElemsY) = &
    K0(1:nVar,1:nElemsX,1:nElemsY) &
  + K1(1:nVar,1:nElemsX,1:nElemsY)*dt
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
K1_WC(1:nVar,1:nElemsX,1:nElemsY)= WCt(1:nVar,1:nElemsX,1:nElemsY)
WC(1:nVar,1:nElemsX,1:nElemsY) = &
    K0_WC(1:nVar,1:nElemsX,1:nElemsY) &
  + K1_WC(1:nVar,1:nElemsX,1:nElemsY)*dt
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
K1_X(1:nVar,1:nElemsX+1,1:nElemsY)= Wt_X(1:nVar,1:nElemsX+1,1:nElemsY) !*NB: +1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = &
    K0_X(1:nVar,1:nElemsX+1,1:nElemsY) &
  + K1_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt

!*Staggering in Y direction
K1_Y(1:nVar,1:nElemsX,1:nElemsY+1)= Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1) !*NB: +1 in X direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = &
    K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) &
  + K1_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt
#endif


!--------------------!
! Second Stage       !
!--------------------!
tStage = t + 1.0*dt
#ifdef ACTIVEFLUX
#ifdef POSTPROCESSINGATEACHSTEP
CALL AF_PostProcessing_Subroutine(tStage)
#endif
#endif
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + &
          0.25*K1(1:nVar,1:nElemsX,1:nElemsY)*dt+0.25*K2(1:nVar,1:nElemsX,1:nElemsY)*dt 
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
K2_WC(1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
WC(1:nVar,1:nElemsX,1:nElemsY) = K0_WC(1:nVar,1:nElemsX,1:nElemsY) + &
          0.25*K1_WC(1:nVar,1:nElemsX,1:nElemsY)*dt+0.25*K2_WC(1:nVar,1:nElemsX,1:nElemsY)*dt 
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
K2_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = K0_X(1:nVar,1:nElemsX+1,1:nElemsY) + &
          0.25*K1_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt+0.25*K2_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt 

!*Staggering in Y direction
K2_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) + &
          0.25*K1_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt+0.25*K2_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt 
#endif

!--------------------!
! Third Stage       !
!--------------------!
tStage = t + 0.5*dt
#ifdef ACTIVEFLUX
#ifdef POSTPROCESSINGATEACHSTEP
CALL AF_PostProcessing_Subroutine(tStage)
#endif
#endif
CALL FVTimeDerivative(tStage)
K3(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + &
          1./6.*K1(1:nVar,1:nElemsX,1:nElemsY)*dt+&
          1./6.*K2(1:nVar,1:nElemsX,1:nElemsY)*dt+&
          2./3.*K3(1:nVar,1:nElemsX,1:nElemsY)*dt
#if defined(CENTEREDPRIMITIVE) && !defined(FIRSTORDERPRIMITIVE)
K3_WC(1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
WC(1:nVar,1:nElemsX,1:nElemsY) = K0_WC(1:nVar,1:nElemsX,1:nElemsY) + &
          1./6.*K1_WC(1:nVar,1:nElemsX,1:nElemsY)*dt+&
          1./6.*K2_WC(1:nVar,1:nElemsX,1:nElemsY)*dt+&
          2./3.*K3_WC(1:nVar,1:nElemsX,1:nElemsY)*dt
#endif
#ifdef ACTIVEFLUX
!*Staggering in X direction
K3_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = K0_X(1:nVar,1:nElemsX+1,1:nElemsY) + &
          1./6.*K1_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt+&
          1./6.*K2_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt+&
          2./3.*K3_X(1:nVar,1:nElemsX+1,1:nElemsY)*dt


!*Staggering in Y direction
K3_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) + &
          1./6.*K1_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt+&
          1./6.*K2_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt+&
          2./3.*K3_Y(1:nVar,1:nElemsX,1:nElemsY+1)*dt

#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationBySSPRK3
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationBySSPRK4(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
USE MOD_FiniteVolume2D_vars,ONLY: K4
USE MOD_FiniteVolume2D_vars,ONLY: K5
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj
!-------------------------------------------------------------------------------!

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! First Stage        !
!--------------------!
tStage = t + 0.0*dt
CALL FVTimeDerivative(tStage)
K1(1:nVar,1:nElemsX,1:nElemsY) = &
    1.00000000000000*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.39175222700392*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K1(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Second Stage       !
!--------------------!
tStage = t + 0.39175222700392*dt
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = &
    0.44437049406734*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.55562950593266*K1(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.36841059262959*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K2(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Third Stage        !
!--------------------!
tStage = t + 0.58607968896780*dt
CALL FVTimeDerivative(tStage)
K3(1:nVar,1:nElemsX,1:nElemsY) = &
    0.62010185138540*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.37989814861460*K2(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.25189177424738*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K3(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Fourth Stage       !
!--------------------!
tStage = t + 0.474542364687*dt
CALL FVTimeDerivative(tStage)
K4(1:nVar,1:nElemsX,1:nElemsY) = &
    0.17807995410773*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.82192004589227*K3(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.54497475021237*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K4(1:nVar,1:nElemsX,1:nElemsY)
K5(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Fifth Stage        !
!--------------------!
tStage  = t + 0.93501063100924*dt
CALL FVTimeDerivative(tStage)
U(1:nVar,1:nElemsX,1:nElemsY)  = &
    0.00683325884039*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.51723167208978*K2(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.12759831133288*K3(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.34833675773694*K4(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.08460416338212*K5(1:nVar,1:nElemsX,1:nElemsY)*dt &
  + 0.22600748319395*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationBySSPRK4
!===============================================================================!
!
!!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByRK65(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: UN0!U^(0)
USE MOD_FiniteVolume2D_vars,ONLY: K0 !F(U^(0))
USE MOD_FiniteVolume2D_vars,ONLY: K1 !F(U^(1))
USE MOD_FiniteVolume2D_vars,ONLY: K2 !F(U^(2))
USE MOD_FiniteVolume2D_vars,ONLY: K3 !F(U^(3))
USE MOD_FiniteVolume2D_vars,ONLY: K4 !F(U^(4))
USE MOD_FiniteVolume2D_vars,ONLY: K5 !F(U^(5))
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(6,5) :: ARK
REAL, DIMENSION(6)   :: bRK, cRK
REAL                 :: tStage
INTEGER              :: ii, jj, nStages
!-------------------------------------------------------------------------------!
nStages = 6
ARK = reshape((/ 0., 0.25, 0.125, 0., 0.1875, -0.42857142857142855, &
                 0., 0.,   0.125, 0., -0.375, 1.1428571428571428,   &
                 0., 0.,   0.,   0.5,  0.375, 0.8571428571428571, &
                 0., 0.,   0.,   0. , 0.5625, -1.7142857142857142, &
                 0., 0.,   0.,   0. , 0.    , 1.1428571428571428 /), shape(ARK))
bRK = (/ 0.077777777777777778, 0., 0.35555555555555556, 0.13333333333333333, 0.35555555555555556, 0.077777777777777778/)
cRK = (/ 0.0, 0.25, 0.25, 0.5, 0.75, 1.0 /)

!--------------------!
! Zero  Stage        !
!--------------------!

UN0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
tStage = t + cRK(1)*dt
CALL FVTimeDerivative(tStage)
K0(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! First Stage        !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +ARK(2,1)*K0(1:nVar,1:nElemsX,1:nElemsY)*dt
tStage = t + cRK(2)*dt
CALL FVTimeDerivative(tStage)
K1(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Second Stage       !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(ARK(3,1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + ARK(3,2)*K1(1:nVar,1:nElemsX,1:nElemsY) )*dt

tStage = t + cRK(3)*dt
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = Ut

!--------------------!
! Third Stage        !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(ARK(4,1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + ARK(4,2)*K1(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(4,3)*K2(1:nVar,1:nElemsX,1:nElemsY) )*dt

tStage = t + cRK(4)*dt
CALL FVTimeDerivative(tStage)
K3(1:nVar,1:nElemsX,1:nElemsY) = Ut

!--------------------!
! Fourth Stage       !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(ARK(5,1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + ARK(5,2)*K1(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(5,3)*K2(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(5,4)*K3(1:nVar,1:nElemsX,1:nElemsY) )*dt

tStage = t + cRK(5)*dt
CALL FVTimeDerivative(tStage)
K4(1:nVar,1:nElemsX,1:nElemsY) = Ut

!--------------------!
! Fifth Stage        !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(ARK(6,1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + ARK(6,2)*K1(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(6,3)*K2(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(6,4)*K3(1:nVar,1:nElemsX,1:nElemsY)  &
    + ARK(6,5)*K4(1:nVar,1:nElemsX,1:nElemsY) )*dt

tStage = t + cRK(6)*dt
CALL FVTimeDerivative(tStage)
K5(1:nVar,1:nElemsX,1:nElemsY) = Ut

!--------------------!
! Final Update       !
!--------------------!
U(1:nVar,1:nElemsX,1:nElemsY) = UN0(1:nVar,1:nElemsX,1:nElemsY) &
  +(bRK(1)*K0(1:nVar,1:nElemsX,1:nElemsY) &
    + bRK(2)*K1(1:nVar,1:nElemsX,1:nElemsY)  &
    + bRK(3)*K2(1:nVar,1:nElemsX,1:nElemsY)  &
    + bRK(4)*K3(1:nVar,1:nElemsX,1:nElemsY)  &
    + bRK(5)*K4(1:nVar,1:nElemsX,1:nElemsY)  &
    + bRK(6)*K5(1:nVar,1:nElemsX,1:nElemsY) )*dt
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByRK65
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByDeC2(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 2
INTEGER, PARAMETER   :: Msteps = 2
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!

thetaDeC = reshape((/ 0.0, 0.5, 0.0, 0.5 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0 , 1.0  /)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY)
  Wp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
  Wp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
  Wp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO


CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  FWp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  FWp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
  FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    !*For staggering in X direction
    !*NB:+1 in X direction
    Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

    !*For staggering in X direction
    !*NB:+1 in Y direction
    Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_X(jj,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) + dt*thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
        !*For staggering in X direction
        !*NB:+1 in X direction
        Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_X(jj,1:nVar,1:nElemsX+1,1:nElemsY)

        !*For staggering in X direction
        !*NB:+1 in Y direction
        Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) + dt*thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
        Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
#endif

      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

    END DO
  END IF

  Up=Ua
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)

#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByDeC2
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByDeC3(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(3,3) :: thetaDeC
REAL, DIMENSION(3)   :: betaDeC
REAL, DIMENSION(3)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 3
INTEGER, PARAMETER   :: MSteps = 3
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!

thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.2083333333333333333, &
                      0.1666666666666666666, &
                      0.0000000000000000000, &
                      0.3333333333333333333, &
                      0.6666666666666666666, &
                      0.0000000000000000000, &
                      -0.0416666666666666666, &
                      0.1666666666666666666  /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0, 0.5 , 1.0/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY)
  Wp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
  Wp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
  Wp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO


CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  FWp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  FWp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
  FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    !*For staggering in X direction
    !*NB:+1 in X direction
    Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

    !*For staggering in X direction
    !*NB:+1 in Y direction
    Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_X(jj,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) + dt*thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
        !*For staggering in X direction
        !*NB:+1 in X direction
        Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_X(jj,1:nVar,1:nElemsX+1,1:nElemsY)

        !*For staggering in X direction
        !*NB:+1 in Y direction
        Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) + dt*thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
        Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
#endif

      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

    END DO
  END IF

  Up=Ua
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)

#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByDeC3
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByDeC4(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(3,3) :: thetaDeC
REAL, DIMENSION(3)   :: betaDeC
REAL, DIMENSION(3)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 4
INTEGER, PARAMETER   :: MSteps = 3
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!

thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.2083333333333333333, &
                      0.1666666666666666666, &
                      0.0000000000000000000, &
                      0.3333333333333333333, &
                      0.6666666666666666666, &
                      0.0000000000000000000, &
                      -0.0416666666666666666, &
                      0.1666666666666666666 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0, 0.5 , 1.0/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY)
  Wp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
  Wp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
  Wp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO


CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  FWp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  FWp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
  FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    !*For staggering in X direction
    !*NB:+1 in X direction
    Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

    !*For staggering in X direction
    !*NB:+1 in Y direction
    Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_X(jj,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) + dt*thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
        !*For staggering in X direction
        !*NB:+1 in X direction
        Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_X(jj,1:nVar,1:nElemsX+1,1:nElemsY)

        !*For staggering in X direction
        !*NB:+1 in Y direction
        Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) + dt*thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
        Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
#endif

      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

    END DO
  END IF

  Up=Ua
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)

#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByDeC4
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByDeC5(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(4,4) :: thetaDeC
REAL, DIMENSION(4)   :: betaDeC
REAL, DIMENSION(4)   :: tStage
INTEGER, PARAMETER   :: KCorr = 5, MSteps = 4
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.1103005664791649049, &
                      0.0730327668541684294, &
                      0.0833333333333333287, &
                      0.0000000000000000000, &
                      0.1896994335208351257, &
                      0.4505740308958107176, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      -0.0339073642291439076, &
                      0.2269672331458314485, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      0.0103005664791649201, &
                      -0.0269672331458315692, &
                      0.0833333333333333287 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 0.2763932022500210639 , 0.7236067977499789361 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY)
  Wp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
  Wp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
  Wp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO


CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  FWp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  FWp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
  FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    !*For staggering in X direction
    !*NB:+1 in X direction
    Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

    !*For staggering in X direction
    !*NB:+1 in Y direction
    Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_X(jj,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) + dt*thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
        !*For staggering in X direction
        !*NB:+1 in X direction
        Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_X(jj,1:nVar,1:nElemsX+1,1:nElemsY)

        !*For staggering in X direction
        !*NB:+1 in Y direction
        Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) + dt*thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
        Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(ii,1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
#endif

      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,1:nElemsX+1,1:nElemsY) = Wt_X(1:nVar,1:nElemsX+1,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,1:nElemsX,1:nElemsY+1) = Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

    END DO
  END IF

  Up=Ua
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)

#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByDeC5
!===============================================================================!
!
!
!===============================================================================!
#ifdef PATANKAR
SUBROUTINE TimeDiscretizationByMPEuler(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp !Pij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestUp !Dij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse !Pij sparse
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse  !Dij sparse
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr = 1, MSteps = 2
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, 0.5000000000000000000, 0.0000000000000000000, 0.5000000000000000000 /), shape(thetaDeC))



! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
END DO

CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
  ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
  DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps

    Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
    DO jj = 1,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
    END DO

    CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

  ELSE
    DO ii = 2,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
      DO jj = 1,MSteps
        Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
      END DO

      CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
      ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
      DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
    END DO
  END IF

  Up=Ua

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByMPEuler
!===============================================================================!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByMPDeC2(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp !Pij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestUp !Dij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse !Pij sparse
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse  !Dij sparse
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr = 2, MSteps = 2
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, 0.5000000000000000000, 0.0000000000000000000, 0.5000000000000000000 /), shape(thetaDeC))



! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
END DO

CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
  ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
  DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps

    Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
    DO jj = 1,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
    END DO

    CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

  ELSE
    DO ii = 2,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
      DO jj = 1,MSteps
        Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
      END DO

      CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
      ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
      DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
    END DO
  END IF

  Up=Ua

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByMPDeC2
!===============================================================================!
!
!

!===============================================================================!
SUBROUTINE TimeDiscretizationByMPDeC5(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp !Pij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestUp !Dij(U^(k-1)) sparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse !Pij sparse
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse  !Dij sparse
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(4,4) :: thetaDeC
REAL, DIMENSION(4)   :: betaDeC
REAL, DIMENSION(4)   :: tStage
INTEGER, PARAMETER   :: KCorr = 5, MSteps = 4
INTEGER              :: ii, jj, kk
!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.1103005664791649049, &
                      0.0730327668541684294, &
                      0.0833333333333333287, &
                      0.0000000000000000000, &
                      0.1896994335208351257, &
                      0.4505740308958107176, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      -0.0339073642291439076, &
                      0.2269672331458314485, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      0.0103005664791649201, &
                      -0.0269672331458315692, &
                      0.0833333333333333287 /), shape(thetaDeC))



! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 0.2763932022500210639 , 0.7236067977499789361 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
END DO

CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
  ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
  DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps

    Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
    DO jj = 1,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
    END DO

    CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

  ELSE
    DO ii = 2,MSteps
      Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(1,2:nVar,1:nElemsX,1:nElemsY)
      DO jj = 1,MSteps
        Ua(ii,2:nVar,1:nElemsX,1:nElemsY) = Ua(ii,2:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,2:nVar,1:nElemsX,1:nElemsY)
      END DO

      CALL MODIFIED_PATANKAR_MATRIX_INVERSION(ii, MSteps, thetaDeC, Up(:,1,1:nElemsX,1:nElemsY), Ua(ii,1,1:nElemsX,1:nElemsY))

    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,2:nVar,1:nElemsX,1:nElemsY) = Ut(2:nVar,1:nElemsX,1:nElemsY)
      ProdUp(ii,1:NNZsparse) = ProductionSparse(1:NNZsparse)
      DestUp(ii,1:NNZsparse) = DestructionSparse(1:NNZsparse)
    END DO
  END IF

  Up=Ua

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByMPDeC5
!===============================================================================!
!
!
!===============================================================================!
! JACOBI: mass is split into diagonal (vector) and off diagonal (CRS)
!===============================================================================!
SUBROUTINE MODIFIED_PATANKAR_MATRIX_INVERSION( m_substep, MSteps, thetaDeC, hp, ha)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars
USE MOD_Mesh,               ONLY: GlobalElem
USE MOD_JacobiIteration,    ONLY: jacobi
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: m_substep
INTEGER,INTENT(IN) :: MSteps
REAL, DIMENSION(MSteps,MSteps), INTENT(IN) :: thetaDeC
REAL, DIMENSION(MSteps,nElemsX,nElemsY), INTENT(IN) :: hp
REAL, DIMENSION(nElemsX,nElemsY), INTENT(INOUT) :: ha
REAL, DIMENSION(NNZsparse) :: Mass !OffDiagonal Matrix
REAL, DIMENSION(nElemsX*nElemsY) :: Diagonal
REAL, DIMENSION(nElemsX*nElemsY) :: rhs, haVec
INTEGER                          :: iterations
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER              :: ii, rr, jj,K,L,KKSparseVect, KLSparseVect, iRowStart, iiLPer, jjLPer
!-------------------------------------------------------------------------------!

Mass = 0.
Diagonal = 1.
KLSparseVect=0
iRowStart = 0

DO jj=1,nElemsY
  DO ii=1,nElemsX
    K = GlobalElem(ii,jj)

    KKSparseVect = SparseIndexMat(K,3)
    
    L = GlobalElem(ii,jj-1)
    iiLPer = MODULO(ii-1,nElemsX)+1
    jjLPer = MODULO(jj-2,nElemsY)+1
    
    KLSparseVect = SparseIndexMat(K,1)

    DO rr=1,MSteps
      IF (thetaDeC(m_substep,rr)>0.) THEN
        Mass(KLSparseVect) = Mass(KLSparseVect) - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )
        Diagonal(K)        = Diagonal(K)        + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,ii,    jj      )) )

      ELSE
        Diagonal(K)        = Diagonal(K)        - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,ii,  jj        )) )
        Mass(KLSparseVect) = Mass(KLSparseVect) + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )

      ENDIF
    END DO
    


    L = GlobalElem(ii-1,jj)           
    iiLPer = MODULO(ii-2,nElemsX)+1  
    jjLPer = MODULO(jj-1,nElemsY)+1
    
    KLSparseVect = SparseIndexMat(K,2)

    DO rr=1,MSteps
      IF (thetaDeC(m_substep,rr)>0.) THEN
        Mass(KLSparseVect) = Mass(KLSparseVect) - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )
        Diagonal(K)        = Diagonal(K)        + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,ii,    jj      )) )

      ELSE
        Diagonal(K)        = Diagonal(K)        - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,ii,  jj        )) )
        Mass(KLSparseVect) = Mass(KLSparseVect) + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )

      ENDIF
    END DO


    
    L = GlobalElem(ii+1,jj)
    iiLPer = MODULO(ii,  nElemsX)+1
    jjLPer = MODULO(jj-1,nElemsY)+1
    
    KLSparseVect = SparseIndexMat(K,4)


    DO rr=1,MSteps
      IF (thetaDeC(m_substep,rr)>0.) THEN
        Mass(KLSparseVect) = Mass(KLSparseVect) - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )
        Diagonal(K)        = Diagonal(K)        + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,ii,    jj      )) )

      ELSE
        Diagonal(K)        = Diagonal(K)        - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,ii,  jj        )) )
        Mass(KLSparseVect) = Mass(KLSparseVect) + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )

      ENDIF
    END DO



    L = GlobalElem(ii,jj+1)
    iiLPer = MODULO(ii-1,nElemsX)+1
    jjLPer = MODULO(jj,  nElemsY)+1

    KLSparseVect = SparseIndexMat(K,5)


    DO rr=1,MSteps
      IF (thetaDeC(m_substep,rr)>0.) THEN
        Mass(KLSparseVect) = Mass(KLSparseVect) - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )
        Diagonal(K)        = Diagonal(K)        + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,ii,    jj      )) )

      ELSE
        Diagonal(K)        = Diagonal(K)        - dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(ProdUP(rr,KLSparseVect),hp(m_substep,ii,  jj        )) )
        Mass(KLSparseVect) = Mass(KLSparseVect) + dt*thetaDeC(m_substep,rr)*( PATANKAR_DIV(DestUP(rr,KLSparseVect),hp(m_substep,iiLPer,jjLPer  )) )

      ENDIF
    END DO


    rhs(GlobalElem(ii,jj))=hp(1,ii,jj)
    
  END DO
END DO

haVec=0.


CALL jacobi(RowStart,ColumnsVector,Mass, Diagonal, rhs, haVec, iterations)
JacobiCounter = JacobiCounter +1
JacobiIterations(JacobiCounter) = iterations


DO jj=1,nElemsY
  DO ii=1,nElemsX
    K = GlobalElem(ii,jj)
    ha(ii,jj)=haVec(K)
  END DO
END DO


!
!-------------------------------------------------------------------------------!
END SUBROUTINE MODIFIED_PATANKAR_MATRIX_INVERSION
!===============================================================================!
!
!
!===============================================================================!
REAL FUNCTION  PATANKAR_DIV( prod , h )
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY
IMPLICIT NONE 
REAL, INTENT(IN) :: prod, h

IF ( h .LT. MIN_DENSITY ) THEN
  PATANKAR_DIV = 0.
ELSE
  PATANKAR_DIV = 2. * h * prod / ( h*h + MAX( h*h , MIN_DENSITY ) )
ENDIF

END  FUNCTION  PATANKAR_DIV 
!===============================================================================!
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
SUBROUTINE TimeDiscretizationByIMEXEuler(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: K0_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: K0_Y
#endif

#ifdef CENTEREDPRIMITIVE
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction     ,ONLY: Laplacian_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_Equation           ,ONLY: BoundaryConditions 
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj, indc, indi, indj, its
INTEGER         :: iterations
#ifdef ACTIVEFLUX
REAL            :: initial_guess_X(1:NRows_X)
REAL            :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL            :: initial_guess(1:NRows_WC)
#endif
REAL            :: px,py,dx,dy
!-------------------------------------------------------------------------------!
!*Debug
!-------------------------------------------------------------------------------!
REAL                                :: residual, dxvx, dyvy, partx, party
#ifdef ACTIVEFLUX
REAL                                :: old_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
REAL                                :: old_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                                :: old_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
REAL                                :: tolerance=1e-10
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug

!-------------------------------------------------------------------------------!

#if(1==1)
!*=============================================
!*=============================================
!*=============================================
!*FOR SAFETY CHECK
CALL BoundaryConditions(t)
#ifdef ACTIVEFLUX
old_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)=W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
old_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)=W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
old_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
!*=============================================
!*=============================================
!*=============================================
#endif

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
K0_X(1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY) !*NB: +1 in X direction
K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1) !*NB: +1 in Y direction
#endif
#ifdef CENTEREDPRIMITIVE
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY) 
#endif

tStage = t 
CALL FVTimeDerivative(tStage)

#ifdef ACTIVEFLUX
CALL Impose_BC_on_Wt()
#endif
#ifdef CENTEREDPRIMITIVE
CALL Impose_BC_on_WCt()
#endif
!*----------------------
!*NB: It is super important to have BCs on W_X, W_Y, Wt_X and Wt_Y
!*With this we have them
!*----------------------

!*----------------------
!*At this point
!*Conservative variables are updated normally
!*Primitive variables are updated differently
!*->ro is updated normally
!*->p  is updated solving a linear system
!*->v  is updated normally but with an extra explicit contribution computed from p
!*----------------------

!*----------------------
!*Normal update of conservative variables <=============
!*----------------------
U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + Ut(1:nVar,1:nElemsX,1:nElemsY)*dt


!*----------------------
!*Update of primitive variables <=============
!*----------------------

!*----------------------
!*->Normal update of ro <=============
!*----------------------
#ifdef ACTIVEFLUX
W_X(1,1:nElemsX+1,1:nElemsY) = K0_X(1,1:nElemsX+1,1:nElemsY) + Wt_X(1,1:nElemsX+1,1:nElemsY)*dt !*NB:+1 in X direction
W_Y(1,1:nElemsX,1:nElemsY+1) = K0_Y(1,1:nElemsX,1:nElemsY+1) + Wt_Y(1,1:nElemsX,1:nElemsY+1)*dt !*NB:+1 in Y direction
#endif
#ifdef CENTEREDPRIMITIVE
WC(1,1:nElemsX,1:nElemsY) = K0_WC(1,1:nElemsX,1:nElemsY) + WCt(1,1:nElemsX,1:nElemsY)*dt 
#endif



#ifdef ACTIVEFLUX
!*----------------------
!*IMEX linear system W_X <=============
!*----------------------
!*Assembly IMEX Matrix Euler
CALL IMEX_Matrix_X(dt)
!*Assembly IMEX rhs Euler
CALL IMEX_Euler_rhs_X(Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),t+dt)

!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif



#if(1==0)
ALLOCATE(sol_debug(1:NRows_X))
iterations=0
CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_debug,iterations,initial_guess_X)

#ifdef COUNTJACOBI
JacobiCounter_X = JacobiCounter_X +1
JacobiIterations_X(JacobiCounter_X) = iterations
#endif

DO indi=1,NRows_X
  IF(ABS(sol_X(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_X"
    PRINT*, indi, sol_X(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug)
#endif



!*Copying sol_X into W_X(4,ii,jj)
CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))


!*----------------------
!*IMEX linear system W_Y <=============
!*----------------------
!*Assembly IMEX Matrix Euler
CALL IMEX_Matrix_Y(dt)
!*Assembly IMEX rhs Euler
CALL IMEX_Euler_rhs_Y(Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),t+dt)

!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, 1000, 1e-12, 1e-6 )
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif




#if(1==0)
ALLOCATE(sol_debug(1:NRows_Y))
iterations=0
CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_debug,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
JacobiCounter_Y = JacobiCounter_Y +1
JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

DO indi=1,NRows_Y
  IF(ABS(sol_Y(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_Y"
    PRINT*, indi, sol_Y(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug)
#endif

!*Copying sol_Y into W_Y(4,ii,jj)
CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))

#if(1==0)
!*=============================================
!*=============================================
!*=============================================
dx=MESH_DX(1)
dy=MESH_DX(2)

PRINT*, "SAFETY CHECK W_X"
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    residual=W_X(4,ii,jj)
#ifdef RELAXATION
    residual=residual-(dt/EPS_LM)**2*Gmm*min_p/max_ro*Laplacian_Central_Order2(W_X(4,ii-1:ii+1,jj-1:jj+1),dx,dy)
#else
    residual=residual-(dt       )**2*Gmm*min_p/max_ro*Laplacian_Central_Order2(W_X(4,ii-1:ii+1,jj-1:jj+1),dx,dy)
#endif
    residual=residual-old_W_X(4,ii,jj)
    residual=residual-dt*Wt_X(4,ii,jj)
    dxvx=First_Derivative_Central_Order2(old_W_X(2,ii-1:ii+1,jj),dx)
    dyvy=First_Derivative_Central_Order2(old_W_X(3,ii,jj-1:jj+1),dy)
    residual=residual+dt*Gmm*min_p*(dxvx+dyvy)
    partx=First_Derivative_Central_Order2(Wt_X(2,ii-1:ii+1,jj),dx)
    party=First_Derivative_Central_Order2(Wt_X(3,ii,jj-1:jj+1),dy)
    residual=residual+dt**2*Gmm*min_p*(partx+party)
    IF (ABS(residual) .GT. tolerance) THEN
      PRINT*, "Mismatch in W_X linear system"
      PRINT*, ii, jj, residual, W_X(4,ii,jj)
      STOP
    END IF
    ! STOP
  END DO
END DO


PRINT*, "SAFETY CHECK W_Y"
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    residual=W_Y(4,ii,jj)
#ifdef RELAXATION
    residual=residual-(dt/EPS_LM)**2*Gmm*min_p/max_ro*Laplacian_Central_Order2(W_Y(4,ii-1:ii+1,jj-1:jj+1),dx,dy)
#else
    residual=residual-(dt       )**2*Gmm*min_p/max_ro*Laplacian_Central_Order2(W_Y(4,ii-1:ii+1,jj-1:jj+1),dx,dy)
#endif
    residual=residual-old_W_Y(4,ii,jj)
    residual=residual-dt*Wt_Y(4,ii,jj)
    dxvx=First_Derivative_Central_Order2(old_W_Y(2,ii-1:ii+1,jj),dx)
    dyvy=First_Derivative_Central_Order2(old_W_Y(3,ii,jj-1:jj+1),dy)
    residual=residual+dt*Gmm*min_p*(dxvx+dyvy)
    partx=First_Derivative_Central_Order2(Wt_Y(2,ii-1:ii+1,jj),dx)
    party=First_Derivative_Central_Order2(Wt_Y(3,ii,jj-1:jj+1),dy)
    residual=residual+dt**2*Gmm*min_p*(partx+party)
    IF (ABS(residual) .GT. tolerance) THEN
      PRINT*, "Mismatch in W_Y linear system"
      PRINT*, ii, jj, residual, W_Y(4,ii,jj)
      STOP
    END IF
    ! STOP
  END DO
END DO
!*=============================================
!*=============================================
!*=============================================
#endif

#endif

#ifdef CENTEREDPRIMITIVE
!*----------------------
!*IMEX linear system WC <=============
!*----------------------
!*Assembly IMEX Matrix Euler 
CALL IMEX_Matrix(dt) 
!*Assembly IMEX rhs Euler
CALL IMEX_Euler_rhs(WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),t+dt)

!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess(1:NRows_WC))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif



#if(1==0)
ALLOCATE(sol_debug(1:NRows_WC))
iterations=0
CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_debug,iterations,initial_guess)

#ifdef COUNTJACOBI
JacobiCounter_WC = JacobiCounter_WC +1
JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

DO indi=1,NRows_X
  IF(ABS(sol_WC(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_WC"
    PRINT*, indi, sol_WC(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug)
#endif



!*Copying sol_WC into WC(4,ii,jj)
CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))


#endif




dx=MESH_DX(1)
dy=MESH_DX(2)

#ifdef ACTIVEFLUX
!*----------------------
!*Update velocity with extra explicit contribution <=============
!*----------------------
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    W_X(2:3,ii,jj) = K0_X(2:3,ii,jj) + Wt_X(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(W_X(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(W_X(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    W_X(2,ii,jj)=W_X(2,ii,jj)-(dt/EPS_LM**2)/max_ro*px
    W_X(3,ii,jj)=W_X(3,ii,jj)-(dt/EPS_LM**2)/max_ro*py
#else
    W_X(2,ii,jj)=W_X(2,ii,jj)-(dt)/max_ro*px
    W_X(3,ii,jj)=W_X(3,ii,jj)-(dt)/max_ro*py
#endif
  END DO
END DO

!*----------------------
!*Update velocity with extra explicit contribution <=============
!*----------------------
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    W_Y(2:3,ii,jj) = K0_Y(2:3,ii,jj) + Wt_Y(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(W_Y(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(W_Y(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    W_Y(2,ii,jj)=W_Y(2,ii,jj)-(dt/EPS_LM**2)/max_ro*px
    W_Y(3,ii,jj)=W_Y(3,ii,jj)-(dt/EPS_LM**2)/max_ro*py
#else
    W_Y(2,ii,jj)=W_Y(2,ii,jj)-(dt)/max_ro*px
    W_Y(3,ii,jj)=W_Y(3,ii,jj)-(dt)/max_ro*py
#endif
  END DO
END DO
#endif

#ifdef CENTEREDPRIMITIVE
!*----------------------
!*Update velocity with extra explicit contribution <=============
!*----------------------
DO jj=1,nElemsY
  DO ii=1,nElemsX
    WC(2:3,ii,jj) = K0_WC(2:3,ii,jj) + WCt(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(WC(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(WC(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    WC(2,ii,jj)=WC(2,ii,jj)-(dt/EPS_LM**2)/max_ro*px
    WC(3,ii,jj)=WC(3,ii,jj)-(dt/EPS_LM**2)/max_ro*py
#else
    WC(2,ii,jj)=WC(2,ii,jj)-(dt)/max_ro*px
    WC(3,ii,jj)=WC(3,ii,jj)-(dt)/max_ro*py
#endif
  END DO
END DO
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXEuler
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE IMEX_Matrix_X(delta_t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY

USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_X
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_X
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_X

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: delta_t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER         :: ii, jj
INTEGER         :: inde, indc, ind, indB, indR, indU, indL
REAL            :: dx, dy
REAL            :: Const_dx, Const_dy
!-------------------------------------------------------------------------------!

dx = MESH_DX(1)
dy = MESH_DX(2)


#ifdef NORMALIZEALL
Const_dx =                      Gmm*min_p/max_ro
Const_dy =           (dx/dy)**2*Gmm*min_p/max_ro
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
Const_dx = (delta_t/EPS_LM        )**2*Gmm*min_p/max_ro
Const_dy = (delta_t/EPS_LM*(dx/dy))**2*Gmm*min_p/max_ro
#else
Const_dx = (delta_t        )**2*Gmm*min_p/max_ro
Const_dy = (delta_t*(dx/dy))**2*Gmm*min_p/max_ro
#endif
#else
#ifdef RELAXATION
Const_dx = (delta_t/EPS_LM/dx)**2*Gmm*min_p/max_ro
Const_dy = (delta_t/EPS_LM/dy)**2*Gmm*min_p/max_ro
#else
Const_dx = (delta_t/dx)**2*Gmm*min_p/max_ro
Const_dy = (delta_t/dy)**2*Gmm*min_p/max_ro
#endif
#endif


!*First equations->Inside the domain

indc=0 !*index cell
inde=0 !*index element

DO jj=1,nElemsY
  DO ii=1,nElemsX+1 

    !*C,B,R,U,L
    indc=indc+1
    inde=inde+1

    RowStart_X(indc)=inde

    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF (((From_1Index_To_2Indices_X(indc,1)-ii) .NE. 0) .OR. ((From_1Index_To_2Indices_X(indc,2)-jj) .NE. 0)) THEN
    !   PRINT*, "Problem in numeration for linear system IMEX", indc
    !   PRINT*, ii,jj
    !   PRINT*, From_1Index_To_2Indices_X(indc,:)
    !   STOP
    ! END IF
    ! PRINT*, ii,jj, From_1Index_To_2Indices_X(indc,:)
    !*-----------------------------------------

    !*NB: The 1-index label of ii,jj is indc
    !*We do not need to retrieve this information as we are sweeping the cells in the same order we numerated them 

    !*C
    ! ind=Neighbours_1Index_X(indc,1) !*NB: This is equal to indc
    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF ((indc-ind) .NE. 0) THEN
    !   PRINT*, "Problem"
    !   PRINT*, ind-indc
    !   STOP
    ! END IF
    !*-----------------------------------------
    Columns_X(inde)  = indc
#ifdef NORMALIZEALL
#ifdef RELAXATION
    Values_X(inde)   = ((dx/delta_t)*EPS_LM)**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_X(inde)   = ((dx/delta_t)       )**2+2.0*Const_dx+2.0*Const_dy
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    Values_X(inde)   = dx**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_X(inde)   = 1.0+2.0*Const_dx+2.0*Const_dy
#endif
    Diagonal_X(indc) = Values_X(inde)

    !*B (ii,jj-1)
    inde=inde+1
    indB=Neighbours_1Index_X(indc,2) !*From_2Indices_To_1Index_X(ii,jj-1)
    Columns_X(inde) = indB
    Values_X(inde)  = -Const_dy

    !*R (ii+1,jj)
    inde=inde+1
    indR=Neighbours_1Index_X(indc,3) !*From_2Indices_To_1Index_X(ii+1,jj)
    Columns_X(inde) = indR
    Values_X(inde)  = -Const_dx

    !*U (ii,jj+1)
    inde=inde+1
    indU=Neighbours_1Index_X(indc,4) !*From_2Indices_To_1Index_X(ii,jj+1)
    Columns_X(inde) = indU
    Values_X(inde)  = -Const_dy

    !*L (ii-1,jj)
    inde=inde+1
    indL=Neighbours_1Index_X(indc,5) !*From_2Indices_To_1Index_X(ii-1,jj)
    Columns_X(inde) = indL
    Values_X(inde)  = -Const_dx

    !*-----------------------------------------
    !*SAFETY CHECK
    !*-----------------------------------------
    ! PRINT*
    ! PRINT*, indc 
    ! PRINT*, "C", indc, Neighbours_1Index_X(indc,1)
    ! PRINT*, "B", indB, Neighbours_1Index_X(indc,2)
    ! PRINT*, "R", indR, Neighbours_1Index_X(indc,3)
    ! PRINT*, "U", indU, Neighbours_1Index_X(indc,4)
    ! PRINT*, "L", indL, Neighbours_1Index_X(indc,5)
    !*-----------------------------------------

  END DO
END DO


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1
  inde=inde+1
  RowStart_X(indc) = inde
  Columns_X(inde)  = indc
  Values_X(inde)   = 1.0
  Diagonal_X(indc) = Values_X(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(ii,nElemsY)
      Values_X(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(ii,1)
      Values_X(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  inde=inde+1
  RowStart_X(indc) = inde
  Columns_X(inde)  = indc
  Values_X(inde)   = 1.0
  Diagonal_X(indc) = Values_X(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1) ! Periodic u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(2,jj) !*Tricky 1=Nx+1
      Values_X(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(nElemsX,jj) !*ALSO TRICKY, SIMMETRY W.R.T. Nx+1
      Values_X(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1
  inde=inde+1
  RowStart_X(indc) = inde
  Columns_X(inde)  = indc
  Values_X(inde)   = 1.0
  Diagonal_X(indc) = Values_X(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(ii,1)
      Values_X(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(ii,nElemsY)
      Values_X(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  inde=inde+1
  RowStart_X(indc) = inde
  Columns_X(inde)  = indc
  Values_X(inde)   = 1.0
  Diagonal_X(indc) = Values_X(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(nElemsX,jj) !*Tricky 1=Nx+1
      Values_X(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(2,jj) !*ALSO TRICKY, SIMMETRY W.R.T. 1
      Values_X(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO


RowStart_X(NRows_X+1)=NNZsparse_X+1

!*-----------------------------------------
!*SAFETY CHECK
!*-----------------------------------------
! PRINT*, Columns_X
! PRINT*, Values_X
! PRINT*, RowStart_X
! PRINT*, Diagonal_X
! STOP
!*-----------------------------------------

END SUBROUTINE IMEX_Matrix_X
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE IMEX_Matrix_Y(delta_t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY

USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_Y
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_Y

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: delta_t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER         :: ii, jj
INTEGER         :: inde, indc, ind, indB, indR, indU, indL
REAL            :: dx, dy
REAL            :: Const_dx, Const_dy
!-------------------------------------------------------------------------------!

dx = MESH_DX(1)
dy = MESH_DX(2)


#ifdef NORMALIZEALL
Const_dx =                      Gmm*min_p/max_ro
Const_dy =           (dx/dy)**2*Gmm*min_p/max_ro
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
Const_dx = (delta_t/EPS_LM        )**2*Gmm*min_p/max_ro
Const_dy = (delta_t/EPS_LM*(dx/dy))**2*Gmm*min_p/max_ro
#else
Const_dx = (delta_t        )**2*Gmm*min_p/max_ro
Const_dy = (delta_t*(dx/dy))**2*Gmm*min_p/max_ro
#endif
#else
#ifdef RELAXATION
Const_dx = (delta_t/EPS_LM/dx)**2*Gmm*min_p/max_ro
Const_dy = (delta_t/EPS_LM/dy)**2*Gmm*min_p/max_ro
#else
Const_dx = (delta_t/dx)**2*Gmm*min_p/max_ro
Const_dy = (delta_t/dy)**2*Gmm*min_p/max_ro
#endif
#endif

!*First equations->Inside the domain

indc=0 !*index cell
inde=0 !*index element

DO jj=1,nElemsY+1
  DO ii=1,nElemsX 

    !*C,B,R,U,L
    indc=indc+1
    inde=inde+1

    RowStart_Y(indc)=inde

    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF (((From_1Index_To_2Indices_Y(indc,1)-ii) .NE. 0) .OR. ((From_1Index_To_2Indices_Y(indc,2)-jj) .NE. 0)) THEN
    !   PRINT*, "Problem in numeration for linear system IMEX", indc
    !   PRINT*, ii,jj
    !   PRINT*, From_1Index_To_2Indices_Y(indc,:)
    !   STOP
    ! END IF
    ! PRINT*, ii,jj, From_1Index_To_2Indices_Y(indc,:)
    !*-----------------------------------------

    !*NB: The 1-index label of ii,jj is indc
    !*We do not need to retrieve this information as we are sweeping the cells in the same order we numerated them 

    !*C
    ind=Neighbours_1Index_Y(indc,1) !*NB: This is equal to indc
    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF ((indc-ind) .NE. 0) THEN
    !   PRINT*, "Problem"
    !   PRINT*, ind-indc
    !   STOP
    ! END IF
    !*-----------------------------------------
    Columns_Y(inde)  = indc
#ifdef NORMALIZEALL
#ifdef RELAXATION
    Values_Y(inde)   = ((dx/delta_t)*EPS_LM)**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_Y(inde)   = ((dx/delta_t)       )**2+2.0*Const_dx+2.0*Const_dy
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    Values_Y(inde)   = dx**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_Y(inde)   = 1.0+2.0*Const_dx+2.0*Const_dy
#endif
    Diagonal_Y(indc) = Values_Y(inde)

    !*B (ii,jj-1)
    inde=inde+1
    indB=Neighbours_1Index_Y(indc,2) !*From_2Indices_To_1Index_Y(ii,jj-1)
    Columns_Y(inde) = indB
    Values_Y(inde)  = -Const_dy

    !*R (ii+1,jj)
    inde=inde+1
    indR=Neighbours_1Index_Y(indc,3) !*From_2Indices_To_1Index_Y(ii+1,jj)
    Columns_Y(inde) = indR
    Values_Y(inde)  = -Const_dx

    !*U (ii,jj+1)
    inde=inde+1
    indU=Neighbours_1Index_Y(indc,4) !*From_2Indices_To_1Index_Y(ii,jj+1)
    Columns_Y(inde) = indU
    Values_Y(inde)  = -Const_dy

    !*L (ii-1,jj)
    inde=inde+1
    indL=Neighbours_1Index_Y(indc,5) !*From_2Indices_To_1Index_Y(ii-1,jj)
    Columns_Y(inde) = indL
    Values_Y(inde)  = -Const_dx

    !*-----------------------------------------
    !*SAFETY CHECK
    !*-----------------------------------------
    ! PRINT*
    ! PRINT*, indc 
    ! PRINT*, "C", indc, Neighbours_1Index_Y(indc,1)
    ! PRINT*, "B", indB, Neighbours_1Index_Y(indc,2)
    ! PRINT*, "R", indR, Neighbours_1Index_Y(indc,3)
    ! PRINT*, "U", indU, Neighbours_1Index_Y(indc,4)
    ! PRINT*, "L", indL, Neighbours_1Index_Y(indc,5)
    !*-----------------------------------------

  END DO
END DO


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX
  indc=indc+1
  inde=inde+1
  RowStart_Y(indc) = inde
  Columns_Y(inde)  = indc
  Values_Y(inde)   = 1.0
  Diagonal_Y(indc) = Values_Y(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(ii,nElemsY) !*Tricky 1=Ny+1
      Values_Y(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(ii,2) !*ALSO TRICKY, SIMMETRY W.R.T. 1
      Values_Y(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  inde=inde+1
  RowStart_Y(indc) = inde
  Columns_Y(inde)  = indc
  Values_Y(inde)   = 1.0
  Diagonal_Y(indc) = Values_Y(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1) ! Periodic u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(1,jj) 
      Values_Y(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(nElemsX,jj)
      Values_Y(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1
  inde=inde+1
  RowStart_Y(indc) = inde
  Columns_Y(inde)  = indc
  Values_Y(inde)   = 1.0
  Diagonal_Y(indc) = Values_Y(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(ii,2) !*ALSO TRICKY, SIMMETRY W.R.T. 1
      Values_Y(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(ii,nElemsY) !*Tricky 1=Ny+1
      Values_Y(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1
  inde=inde+1
  RowStart_Y(indc) = inde
  Columns_Y(inde)  = indc
  Values_Y(inde)   = 1.0
  Diagonal_Y(indc) = Values_Y(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(nElemsX,jj)
      Values_Y(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(1,jj)
      Values_Y(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO


RowStart_Y(NRows_Y+1)=NNZsparse_Y+1

!*-----------------------------------------
!*SAFETY CHECK
!*-----------------------------------------
! PRINT*, Columns_Y
! PRINT*, Values_Y
! PRINT*, RowStart_Y
! PRINT*, Diagonal_Y
! STOP
!*-----------------------------------------

END SUBROUTINE IMEX_Matrix_Y
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE IMEX_Matrix(delta_t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY

USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_WC
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_WC
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_WC

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: delta_t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER         :: ii, jj
INTEGER         :: inde, indc, ind, indB, indR, indU, indL
REAL            :: dx, dy
REAL            :: Const_dx, Const_dy
!-------------------------------------------------------------------------------!

dx = MESH_DX(1)
dy = MESH_DX(2)


#ifdef NORMALIZEALL
Const_dx =                      Gmm*min_p/max_ro
Const_dy =           (dx/dy)**2*Gmm*min_p/max_ro
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
Const_dx = (delta_t/EPS_LM        )**2*Gmm*min_p/max_ro
Const_dy = (delta_t/EPS_LM*(dx/dy))**2*Gmm*min_p/max_ro
#else
Const_dx = (delta_t        )**2*Gmm*min_p/max_ro
Const_dy = (delta_t*(dx/dy))**2*Gmm*min_p/max_ro
#endif
#else
#ifdef RELAXATION
Const_dx = (delta_t/EPS_LM/dx)**2*Gmm*min_p/max_ro
Const_dy = (delta_t/EPS_LM/dy)**2*Gmm*min_p/max_ro
#else
Const_dx = (delta_t/dx)**2*Gmm*min_p/max_ro
Const_dy = (delta_t/dy)**2*Gmm*min_p/max_ro
#endif
#endif


!*First equations->Inside the domain

indc=0 !*index cell
inde=0 !*index element

DO jj=1,nElemsY
  DO ii=1,nElemsX

    !*C,B,R,U,L
    indc=indc+1
    inde=inde+1

    RowStart_WC(indc)=inde

    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF (((From_1Index_To_2Indices_WC(indc,1)-ii) .NE. 0) .OR. ((From_1Index_To_2Indices_WC(indc,2)-jj) .NE. 0)) THEN
    !   PRINT*, "Problem in numeration for linear system IMEX", indc
    !   PRINT*, ii,jj
    !   PRINT*, From_1Index_To_2Indices_WC(indc,:)
    !   STOP
    ! END IF
    ! PRINT*, ii,jj, From_1Index_To_2Indices_WC(indc,:)
    !*-----------------------------------------

    !*NB: The 1-index label of ii,jj is indc
    !*We do not need to retrieve this information as we are sweeping the cells in the same order we numerated them 

    !*C
    ! ind=Neighbours_1Index_WC(indc,1) !*NB: This is equal to indc
    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF ((indc-ind) .NE. 0) THEN
    !   PRINT*, "Problem"
    !   PRINT*, ind-indc
    !   STOP
    ! END IF
    !*-----------------------------------------
    Columns_WC(inde)  = indc
#ifdef NORMALIZEALL
#ifdef RELAXATION
    Values_WC(inde)   = ((dx/delta_t)*EPS_LM)**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_WC(inde)   = ((dx/delta_t)       )**2+2.0*Const_dx+2.0*Const_dy
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    Values_WC(inde)   = dx**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_WC(inde)   = 1.0+2.0*Const_dx+2.0*Const_dy
#endif
    Diagonal_WC(indc) = Values_WC(inde)

    !*B (ii,jj-1)
    inde=inde+1
    indB=Neighbours_1Index_WC(indc,2) !*From_2Indices_To_1Index_WC(ii,jj-1)
    Columns_WC(inde) = indB
    Values_WC(inde)  = -Const_dy

    !*R (ii+1,jj)
    inde=inde+1
    indR=Neighbours_1Index_WC(indc,3) !*From_2Indices_To_1Index_WC(ii+1,jj)
    Columns_WC(inde) = indR
    Values_WC(inde)  = -Const_dx

    !*U (ii,jj+1)
    inde=inde+1
    indU=Neighbours_1Index_WC(indc,4) !*From_2Indices_To_1Index_WC(ii,jj+1)
    Columns_WC(inde) = indU
    Values_WC(inde)  = -Const_dy

    !*L (ii-1,jj)
    inde=inde+1
    indL=Neighbours_1Index_WC(indc,5) !*From_2Indices_To_1Index_WC(ii-1,jj)
    Columns_WC(inde) = indL
    Values_WC(inde)  = -Const_dx

    !*-----------------------------------------
    !*SAFETY CHECK
    !*-----------------------------------------
    ! PRINT*
    ! PRINT*, indc 
    ! PRINT*, "C", indc, Neighbours_1Index_WC(indc,1)
    ! PRINT*, "B", indB, Neighbours_1Index_WC(indc,2)
    ! PRINT*, "R", indR, Neighbours_1Index_WC(indc,3)
    ! PRINT*, "U", indU, Neighbours_1Index_WC(indc,4)
    ! PRINT*, "L", indL, Neighbours_1Index_WC(indc,5)
    !*-----------------------------------------

  END DO
END DO


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1
  inde=inde+1
  RowStart_WC(indc) = inde
  Columns_WC(inde)  = indc
  Values_WC(inde)   = 1.0
  Diagonal_WC(indc) = Values_WC(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_WC(inde) = From_2Indices_To_1Index_WC(ii,nElemsY)
      Values_WC(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_WC(inde) = From_2Indices_To_1Index_WC(ii,1)
      Values_WC(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY
  indc=indc+1
  inde=inde+1
  RowStart_WC(indc) = inde
  Columns_WC(inde)  = indc
  Values_WC(inde)   = 1.0
  Diagonal_WC(indc) = Values_WC(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1) ! Periodic u=u
      inde=inde+1
      Columns_WC(inde) = From_2Indices_To_1Index_WC(1,jj) !*ADJUSTED BC
      Values_WC(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_WC(inde) = From_2Indices_To_1Index_WC(nElemsX,jj) !*ADJUSTED BC !*NO CHANGES
      Values_WC(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX
  indc=indc+1
  inde=inde+1
  RowStart_WC(indc) = inde
  Columns_WC(inde)  = indc
  Values_WC(inde)   = 1.0
  Diagonal_WC(indc) = Values_WC(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_WC(inde) = From_2Indices_To_1Index_WC(ii,1)
      Values_WC(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_WC(inde) = From_2Indices_To_1Index_WC(ii,nElemsY)
      Values_WC(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  inde=inde+1
  RowStart_WC(indc) = inde
  Columns_WC(inde)  = indc
  Values_WC(inde)   = 1.0
  Diagonal_WC(indc) = Values_WC(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_WC(inde) = From_2Indices_To_1Index_WC(nElemsX,jj) !*ADJUSTED BC !*NO CHANGES
      Values_WC(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_WC(inde) = From_2Indices_To_1Index_WC(1,jj) !*ADJUSTED BC
      Values_WC(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO


RowStart_WC(NRows_WC+1)=NNZsparse_WC+1

!*-----------------------------------------
!*SAFETY CHECK
!*-----------------------------------------
! PRINT*, Columns_WC
! PRINT*, Values_WC
! PRINT*, RowStart_WC
! PRINT*, Diagonal_WC
! STOP
!*-----------------------------------------

END SUBROUTINE IMEX_Matrix
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE IMEX_Euler_rhs_X(Wt,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), INTENT(IN) :: Wt
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: et, etx, ety !*Extra term=et=(ro*-ro)/ro/ro* divided by (in case relaxation is active) eps_lm**2
REAL                  :: supp(-1:1)
REAL                  :: div_v_grad_v
REAL                  :: ro_p_contribution
REAL                  :: dx_part, dy_part, div_vn
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)

!*RMK: +1 in X direction

rhs_X=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=(W_X(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)*EPS_LM)**2
#else
    rhs_X(indc)=(W_X(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=(W_X(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_X(indc)=W_X(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    !*NB: It is important to have BCs here
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

    ux=First_Derivative_Central_Order2(W_X(2,ii-1:ii+1,jj),dx)
    vy=First_Derivative_Central_Order2(W_X(3,ii,jj-1:jj+1),dy)

    div_vn=ux+vy

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-dt*Gmm*min_p*(div_vn*((dx/dt)*EPS_LM)**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*((dx/dt)*EPS_LM)**2+dy_part*((dx/dt)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_X(indc)=rhs_X(indc)-dt*Gmm*min_p*(div_vn*((dx/dt)       )**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*((dx/dt)       )**2+dy_part*((dx/dt)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=rhs_X(indc)-dt*Gmm*min_p*(div_vn*dx**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_X(indc)=rhs_X(indc)-dt*Gmm*min_p*div_vn                 &
                          &-(dt**2)*Gmm*min_p*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

  END DO !ii
END DO !jj





!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEX_Euler_rhs_X
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE IMEX_Euler_rhs_Y(Wt,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), INTENT(IN) :: Wt
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: et, etx, ety !*Extra term=et=(ro*-ro)/ro/ro* divided by (in case relaxation is active) eps_lm**2
REAL                  :: supp(-1:1)
REAL                  :: div_v_grad_v
REAL                  :: ro_p_contribution
REAL                  :: dx_part, dy_part, div_vn
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)

!*RMK: +1 in Y direction

rhs_Y=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=(W_Y(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)*EPS_LM)**2
#else
    rhs_Y(indc)=(W_Y(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=(W_Y(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_Y(indc)=W_Y(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc


    !*Alex's suggestion with upwinded (already computed operators)
    !*NB: It is important to have BCs here
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

    ux=First_Derivative_Central_Order2(W_Y(2,ii-1:ii+1,jj),dx)
    vy=First_Derivative_Central_Order2(W_Y(3,ii,jj-1:jj+1),dy)

    div_vn=ux+vy

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-dt*Gmm*min_p*(div_vn*((dx/dt)*EPS_LM)**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*((dx/dt)*EPS_LM)**2+dy_part*((dx/dt)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_Y(indc)=rhs_Y(indc)-dt*Gmm*min_p*(div_vn*((dx/dt)       )**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*((dx/dt)       )**2+dy_part*((dx/dt)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=rhs_Y(indc)-dt*Gmm*min_p*(div_vn*dx**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_Y(indc)=rhs_Y(indc)-dt*Gmm*min_p*div_vn                 &
                          &-(dt**2)*Gmm*min_p*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

  END DO !ii
END DO !jj





!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEX_Euler_rhs_Y
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE IMEX_Euler_rhs(Wt,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), INTENT(IN) :: Wt
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: et, etx, ety !*Extra term=et=(ro*-ro)/ro/ro* divided by (in case relaxation is active) eps_lm**2
REAL                  :: supp(-1:1)
REAL                  :: div_v_grad_v
REAL                  :: ro_p_contribution
REAL                  :: dx_part, dy_part, div_vn
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)

!*RMK: +1 in X direction

rhs_WC=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=(WC(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)*EPS_LM)**2
#else
    rhs_WC(indc)=(WC(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=(WC(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_WC(indc)=WC(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    !*NB: It is important to have BCs here
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

    ux=First_Derivative_Central_Order2(WC(2,ii-1:ii+1,jj),dx)
    vy=First_Derivative_Central_Order2(WC(3,ii,jj-1:jj+1),dy)

    div_vn=ux+vy

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)-dt*Gmm*min_p*(div_vn*((dx/dt)*EPS_LM)**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*((dx/dt)*EPS_LM)**2+dy_part*((dx/dt)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_WC(indc)=rhs_WC(indc)-dt*Gmm*min_p*(div_vn*((dx/dt)       )**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*((dx/dt)       )**2+dy_part*((dx/dt)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=rhs_WC(indc)-dt*Gmm*min_p*(div_vn*dx**2)                 &
                          &-(dt**2)*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_WC(indc)=rhs_WC(indc)-dt*Gmm*min_p*div_vn                 &
                          &-(dt**2)*Gmm*min_p*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

  END DO !ii
END DO !jj





!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_WC(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_WC(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_WC(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEX_Euler_rhs
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE IMEX_Euler_rhs_dt_input(Wt,t,delta_t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), INTENT(IN) :: Wt
REAL                                                                         , INTENT(IN) :: t
REAL                                                                         , INTENT(IN) :: delta_t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: et, etx, ety !*Extra term=et=(ro*-ro)/ro/ro* divided by (in case relaxation is active) eps_lm**2
REAL                  :: supp(-1:1)
REAL                  :: div_v_grad_v
REAL                  :: ro_p_contribution
REAL                  :: dx_part, dy_part, div_vn
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)

!*RMK: +1 in X direction

rhs_WC=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=(WC(4,ii,jj) + delta_t*Wt(4,ii,jj))*((dx/delta_t)*EPS_LM)**2
#else
    rhs_WC(indc)=(WC(4,ii,jj) + delta_t*Wt(4,ii,jj))*((dx/delta_t)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=(WC(4,ii,jj) + delta_t*Wt(4,ii,jj))*dx**2
#else
    rhs_WC(indc)=WC(4,ii,jj) + delta_t*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    !*NB: It is important to have BCs here
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

    ux=First_Derivative_Central_Order2(WC(2,ii-1:ii+1,jj),dx)
    vy=First_Derivative_Central_Order2(WC(3,ii,jj-1:jj+1),dy)

    div_vn=ux+vy

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)-delta_t*Gmm*min_p*(div_vn*((dx/delta_t)*EPS_LM)**2)                 &
                          &-(delta_t**2)*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_WC(indc)=rhs_WC(indc)-delta_t*Gmm*min_p*(div_vn*((dx/delta_t)       )**2)                 &
                          &-(delta_t**2)*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=rhs_WC(indc)-delta_t*Gmm*min_p*(div_vn*dx**2)                 &
                          &-(delta_t**2)*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_WC(indc)=rhs_WC(indc)-delta_t*Gmm*min_p*div_vn                 &
                          &-(delta_t**2)*Gmm*min_p*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

  END DO !ii
END DO !jj





!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_WC(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_WC(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_WC(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEX_Euler_rhs_dt_input
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC2(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
#endif

USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro


!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 2
INTEGER, PARAMETER   :: MSteps = 2
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0, 0.5, 0.0, 0.5 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0 , 1.0  /)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

!*Initialiting initial guesses for the linear systems
#ifdef ACTIVEFLUX
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
  Up=Ua

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
    U(1:nVar,1:nElemsX,1:nElemsY)                                       = Up(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif

#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==0)
      !*Stage dependent max_ro and min_p
      PRINT*, "Stage dependent max_ro and min_p"
      PRINT*, "Keep it deactivated"
      STOP
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivative(tStage(ii))

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif


#endif

#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
    FLUXES_U   = 0.0

#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif

#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of conserved variables
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of conserved variables
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Up(1,1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))


#if(1==0)
    !*==========================
    !*IT COULD BE DELETED
    !*DIRICHLET BCs CAUSES ISSUES
    !*==========================

    ! PRINT*, "Check (Dirichlet) BC X"

    indc=(nElemsX+1)*nElemsY

    !*Later equations->Boundary conditions
    !*->Down boundary
    indj=0
    DO indi=1,nElemsX+1 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    !*->Right boundary
    indi=nElemsX+2
    DO indj=1,nElemsY
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*->Top boundary
    indj=nElemsY+1
    DO indi=1,nElemsX+1 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO



    !*->Left boundary
    indi=0
    DO indj=1,nElemsY
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    ! PRINT*, "Check (Dirichlet) BC Y"

    indc=nElemsX*(nElemsY+1)


    !*Later equations->Boundary conditions
    !*->Down boundary
    indj=0
    DO indi=1,nElemsX 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    !*->Right boundary
    indi=nElemsX+1
    DO indj=1,nElemsY+1
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*->Top boundary
    indj=nElemsY+2
    DO indi=1,nElemsX 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO



    !*->Left boundary
    indi=0
    DO indj=1,nElemsY+1
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*==========================
    !*IT COULD BE DELETED
    !*DIRICHLET BCs CAUSES ISSUES
    !*==========================
#endif

#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, MIN(50,NRows-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_WC into Wa_WC(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif



  END DO


END DO


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC2
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(CENTEREDPRIMITIVE)
SUBROUTINE IMEXDeC2_rho_p_stage_dependent_overwriting_stages(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
#endif

USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D     ,ONLY: Overwrite_WC_from_U
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 2
INTEGER, PARAMETER   :: MSteps = 2
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0, 0.5, 0.0, 0.5 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0 , 1.0  /)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

!*Initialiting initial guesses for the linear systems
#ifdef ACTIVEFLUX
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
  Up=Ua

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
    U(1:nVar,1:nElemsX,1:nElemsY)                                       = Up(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif

#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==1)
      !*Stage dependent max_ro and min_p
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivative(tStage(ii))

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif


#endif

#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
    FLUXES_U   = 0.0

#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif

#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of conserved variables
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of conserved variables
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Up(1,1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))


#if(1==0)
    !*==========================
    !*IT COULD BE DELETED
    !*DIRICHLET BCs CAUSES ISSUES
    !*==========================

    ! PRINT*, "Check (Dirichlet) BC X"

    indc=(nElemsX+1)*nElemsY

    !*Later equations->Boundary conditions
    !*->Down boundary
    indj=0
    DO indi=1,nElemsX+1 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    !*->Right boundary
    indi=nElemsX+2
    DO indj=1,nElemsY
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*->Top boundary
    indj=nElemsY+1
    DO indi=1,nElemsX+1 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO



    !*->Left boundary
    indi=0
    DO indj=1,nElemsY
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    ! PRINT*, "Check (Dirichlet) BC Y"

    indc=nElemsX*(nElemsY+1)


    !*Later equations->Boundary conditions
    !*->Down boundary
    indj=0
    DO indi=1,nElemsX 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    !*->Right boundary
    indi=nElemsX+1
    DO indj=1,nElemsY+1
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*->Top boundary
    indj=nElemsY+2
    DO indi=1,nElemsX 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO



    !*->Left boundary
    indi=0
    DO indj=1,nElemsY+1
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*==========================
    !*IT COULD BE DELETED
    !*DIRICHLET BCs CAUSES ISSUES
    !*==========================
#endif

#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, MIN(50,NRows-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_WC into Wa_WC(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifndef DONOTOVERWRITEPRIMITIVE
  !*To prevent double post-processing of last stage.
  !*Remember that there is an extra post-processing of the final solution outside
  IF (kk .NE. KCorr) THEN 
    ! PRINT*, kk
    WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
    U(1:nVar,1:nElemsX,1:nElemsY)  = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
    CALL Overwrite_WC_from_U()
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY)    = U(1:nVar,1:nElemsX,1:nElemsY) 
  END IF
#endif

  END DO


END DO


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE IMEXDeC2_rho_p_stage_dependent_overwriting_stages
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC3(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
#endif

USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro


!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(3,3) :: thetaDeC
REAL, DIMENSION(3)   :: betaDeC
REAL, DIMENSION(3)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 3
INTEGER, PARAMETER   :: MSteps = 3
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.2083333333333333333, &
                      0.1666666666666666666, &
                      0.0000000000000000000, &
                      0.3333333333333333333, &
                      0.6666666666666666666, &
                      0.0000000000000000000, &
                      -0.0416666666666666666, &
                      0.1666666666666666666 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0, 0.5 , 1.0/)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

!*Initialiting initial guesses for the linear systems
#ifdef ACTIVEFLUX
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
  Up=Ua

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
    U(1:nVar,1:nElemsX,1:nElemsY)                                       = Up(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif

#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==0)
      !*Stage dependent max_ro and min_p
      PRINT*, "Stage dependent max_ro and min_p"
      PRINT*, "Keep it deactivated"
      STOP
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivative(tStage(ii))

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif


#endif

#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
    FLUXES_U   = 0.0

#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif

#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of conserved variables
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of conserved variables
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Up(1,1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))

#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, MIN(50,NRows-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_WC into Wa_WC(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif



  END DO


END DO


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC3
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC4(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
#endif

USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro


!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(3,3) :: thetaDeC
REAL, DIMENSION(3)   :: betaDeC
REAL, DIMENSION(3)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 4
INTEGER, PARAMETER   :: MSteps = 3
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.2083333333333333333, &
                      0.1666666666666666666, &
                      0.0000000000000000000, &
                      0.3333333333333333333, &
                      0.6666666666666666666, &
                      0.0000000000000000000, &
                      -0.0416666666666666666, &
                      0.1666666666666666666 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0, 0.5 , 1.0/)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

!*Initialiting initial guesses for the linear systems
#ifdef ACTIVEFLUX
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
  Up=Ua

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
    U(1:nVar,1:nElemsX,1:nElemsY)                                       = Up(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif

#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==0)
      !*Stage dependent max_ro and min_p
      PRINT*, "Stage dependent max_ro and min_p"
      PRINT*, "Keep it deactivated"
      STOP
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivative(tStage(ii))

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif


#endif

#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
    FLUXES_U   = 0.0

#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif

#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of conserved variables
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of conserved variables
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Up(1,1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))

#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, MIN(50,NRows-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_WC into Wa_WC(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif



  END DO


END DO


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC4
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC5(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
#endif

USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro


!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(4,4) :: thetaDeC
REAL, DIMENSION(4)   :: betaDeC
REAL, DIMENSION(4)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 5
INTEGER, PARAMETER   :: MSteps = 4
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.1103005664791649049, &
                      0.0730327668541684294, &
                      0.0833333333333333287, &
                      0.0000000000000000000, &
                      0.1896994335208351257, &
                      0.4505740308958107176, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      -0.0339073642291439076, &
                      0.2269672331458314485, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      0.0103005664791649201, &
                      -0.0269672331458315692, &
                      0.0833333333333333287 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 0.2763932022500210639 , 0.7236067977499789361 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

!*Initialiting initial guesses for the linear systems
#ifdef ACTIVEFLUX
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
  Up=Ua

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
    U(1:nVar,1:nElemsX,1:nElemsY)                                       = Up(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif

#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==0)
      !*Stage dependent max_ro and min_p
      PRINT*, "Stage dependent max_ro and min_p"
      PRINT*, "Keep it deactivated"
      STOP
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivative(tStage(ii))

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif


#endif

#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
    FLUXES_U   = 0.0

#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif

#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of conserved variables
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of conserved variables
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Up(1,1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))

#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, MIN(50,NRows-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_WC into Wa_WC(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif



  END DO


END DO


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif

#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC5
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXEuler_Implicit_Update_Conserved(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly

USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: K0_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: K0_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction     ,ONLY: Laplacian_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_Equation           ,ONLY: BoundaryConditions 
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj, indc, indi, indj, its
INTEGER         :: iterations
#ifdef ACTIVEFLUX
REAL            :: initial_guess_X(1:NRows_X)
REAL            :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL            :: initial_guess_WC(1:NRows_WC)
#endif
REAL            :: px,py,dx,dy
!-------------------------------------------------------------------------------!
!*Debug
!-------------------------------------------------------------------------------!
REAL                                :: residual, dxvx, dyvy, partx, party
#ifdef ACTIVEFLUX
REAL                                :: old_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
REAL                                :: old_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                                :: old_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
REAL                                :: tolerance=1e-10
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug

!-------------------------------------------------------------------------------!

#if(1==1)
!*=============================================
!*=============================================
!*=============================================
!*FOR SAFETY CHECK
CALL BoundaryConditions(t)
#ifdef ACTIVEFLUX
old_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)=W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
old_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)=W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
old_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
!*=============================================
!*=============================================
!*=============================================
#endif

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
K0_X(1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY) !*NB: +1 in X direction
K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1) !*NB: +1 in Y direction
#endif
#ifdef CENTEREDPRIMITIVE
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

tStage = t 
CALL FVTimeDerivativePrimitiveOnly(tStage)

#ifdef ACTIVEFLUX
CALL Impose_BC_on_Wt()
#endif

#ifdef CENTEREDPRIMITIVE
CALL Impose_BC_on_WCt()
#endif


!*----------------------
!*NB: It is super important to have BCs on W_X, W_Y, Wt_X and Wt_Y
!*With this we have them
!*----------------------

!*----------------------
!*At this point
!*Primitive variables are updated as follows
!*->ro is updated normally
!*->p  is updated solving a linear system
!*->v  is updated normally but with an extra explicit contribution computed from p
!*----------------------

!*----------------------
!*Update of primitive variables <=============
!*----------------------

!*----------------------
!*->Normal update of ro <=============
!*----------------------
#ifdef ACTIVEFLUX
W_X(1,1:nElemsX+1,1:nElemsY) = K0_X(1,1:nElemsX+1,1:nElemsY) + Wt_X(1,1:nElemsX+1,1:nElemsY)*dt !*NB:+1 in X direction
W_Y(1,1:nElemsX,1:nElemsY+1) = K0_Y(1,1:nElemsX,1:nElemsY+1) + Wt_Y(1,1:nElemsX,1:nElemsY+1)*dt !*NB:+1 in Y direction
#endif
#ifdef CENTEREDPRIMITIVE
WC(1,1:nElemsX,1:nElemsY) = K0_WC(1,1:nElemsX,1:nElemsY) + WCt(1,1:nElemsX,1:nElemsY)*dt 
#endif


!*=============================
#ifdef ACTIVEFLUX
!*----------------------
!*IMEX linear system W_X <=============
!*----------------------
!*Assembly IMEX Matrix Euler
CALL IMEX_Matrix_X(dt)
!*Assembly IMEX rhs Euler
CALL IMEX_Euler_rhs_X(Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),t+dt)

!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))


!*----------------------
#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif
!*----------------------

!*----------------------
#if(1==0)
ALLOCATE(sol_debug(1:NRows_X))
iterations=0
CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_debug,iterations,initial_guess_X)

#ifdef COUNTJACOBI
JacobiCounter_X = JacobiCounter_X +1
JacobiIterations_X(JacobiCounter_X) = iterations
#endif

DO indi=1,NRows_X
  IF(ABS(sol_X(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_X"
    PRINT*, indi, sol_X(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug)
#endif
!*----------------------


!*Copying sol_X into W_X(4,ii,jj)
CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))


!*----------------------
!*IMEX linear system W_Y <=============
!*----------------------
!*Assembly IMEX Matrix Euler
CALL IMEX_Matrix_Y(dt)
!*Assembly IMEX rhs Euler
CALL IMEX_Euler_rhs_Y(Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),t+dt)

!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))


!*----------------------
#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, 1000, 1e-12, 1e-6 )
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif
!*----------------------


!*----------------------
#if(1==0)
ALLOCATE(sol_debug(1:NRows_Y))
iterations=0
CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_debug,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
JacobiCounter_Y = JacobiCounter_Y +1
JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

DO indi=1,NRows_Y
  IF(ABS(sol_Y(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_Y"
    PRINT*, indi, sol_Y(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug)
#endif
!*----------------------


!*Copying sol_Y into W_Y(4,ii,jj)
CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))

!*----------------------
#if(1==0)
!*=============================================
!*=============================================
!*=============================================
dx=MESH_DX(1)
dy=MESH_DX(2)

PRINT*, "SAFETY CHECK W_X"
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    residual=W_X(4,ii,jj)
#ifdef RELAXATION
    residual=residual-(dt/EPS_LM)**2*Gmm*min_p/max_ro*Laplacian_Central_Order2(W_X(4,ii-1:ii+1,jj-1:jj+1),dx,dy)
#else
    residual=residual-(dt       )**2*Gmm*min_p/max_ro*Laplacian_Central_Order2(W_X(4,ii-1:ii+1,jj-1:jj+1),dx,dy)
#endif
    residual=residual-old_W_X(4,ii,jj)
    residual=residual-dt*Wt_X(4,ii,jj)
    dxvx=First_Derivative_Central_Order2(old_W_X(2,ii-1:ii+1,jj),dx)
    dyvy=First_Derivative_Central_Order2(old_W_X(3,ii,jj-1:jj+1),dy)
    residual=residual+dt*Gmm*min_p*(dxvx+dyvy)
    partx=First_Derivative_Central_Order2(Wt_X(2,ii-1:ii+1,jj),dx)
    party=First_Derivative_Central_Order2(Wt_X(3,ii,jj-1:jj+1),dy)
    residual=residual+dt**2*Gmm*min_p*(partx+party)
    IF (ABS(residual) .GT. tolerance) THEN
      PRINT*, "Mismatch in W_X linear system"
      PRINT*, ii, jj, residual, W_X(4,ii,jj)
      STOP
    END IF
    ! STOP
  END DO
END DO


PRINT*, "SAFETY CHECK W_Y"
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    residual=W_Y(4,ii,jj)
#ifdef RELAXATION
    residual=residual-(dt/EPS_LM)**2*Gmm*min_p/max_ro*Laplacian_Central_Order2(W_Y(4,ii-1:ii+1,jj-1:jj+1),dx,dy)
#else
    residual=residual-(dt       )**2*Gmm*min_p/max_ro*Laplacian_Central_Order2(W_Y(4,ii-1:ii+1,jj-1:jj+1),dx,dy)
#endif
    residual=residual-old_W_Y(4,ii,jj)
    residual=residual-dt*Wt_Y(4,ii,jj)
    dxvx=First_Derivative_Central_Order2(old_W_Y(2,ii-1:ii+1,jj),dx)
    dyvy=First_Derivative_Central_Order2(old_W_Y(3,ii,jj-1:jj+1),dy)
    residual=residual+dt*Gmm*min_p*(dxvx+dyvy)
    partx=First_Derivative_Central_Order2(Wt_Y(2,ii-1:ii+1,jj),dx)
    party=First_Derivative_Central_Order2(Wt_Y(3,ii,jj-1:jj+1),dy)
    residual=residual+dt**2*Gmm*min_p*(partx+party)
    IF (ABS(residual) .GT. tolerance) THEN
      PRINT*, "Mismatch in W_Y linear system"
      PRINT*, ii, jj, residual, W_Y(4,ii,jj)
      STOP
    END IF
    ! STOP
  END DO
END DO
!*=============================================
!*=============================================
!*=============================================
#endif
!*----------------------

#endif

#ifdef CENTEREDPRIMITIVE
!*----------------------
!*IMEX linear system WC <=============
!*----------------------
!*Assembly IMEX Matrix Euler
CALL IMEX_Matrix(dt)
!*Assembly IMEX rhs Euler
CALL IMEX_Euler_rhs(WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),t+dt)


!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))

!*----------------------
#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif
!*----------------------

!*----------------------
#if(1==0)
ALLOCATE(sol_debug(1:NRows_WC))
iterations=0
CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_debug,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
JacobiCounter_WC = JacobiCounter_WC +1
JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

DO indi=1,NRows_WC
  IF(ABS(sol_WC(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in WC"
    PRINT*, indi, sol_WC(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug)
#endif
!*----------------------


!*Copying sol_WC into WC(4,ii,jj)
CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))

#endif



dx=MESH_DX(1)
dy=MESH_DX(2)


!*----------------------
!*Update velocity with extra explicit contribution <=============
!*----------------------
#ifdef ACTIVEFLUX
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    W_X(2:3,ii,jj) = K0_X(2:3,ii,jj) + Wt_X(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(W_X(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(W_X(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    W_X(2,ii,jj)=W_X(2,ii,jj)-(dt/EPS_LM**2)/max_ro*px
    W_X(3,ii,jj)=W_X(3,ii,jj)-(dt/EPS_LM**2)/max_ro*py
#else
    W_X(2,ii,jj)=W_X(2,ii,jj)-(dt)/max_ro*px
    W_X(3,ii,jj)=W_X(3,ii,jj)-(dt)/max_ro*py
#endif
  END DO
END DO

!*----------------------
!*Update velocity with extra explicit contribution <=============
!*----------------------
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    W_Y(2:3,ii,jj) = K0_Y(2:3,ii,jj) + Wt_Y(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(W_Y(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(W_Y(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    W_Y(2,ii,jj)=W_Y(2,ii,jj)-(dt/EPS_LM**2)/max_ro*px
    W_Y(3,ii,jj)=W_Y(3,ii,jj)-(dt/EPS_LM**2)/max_ro*py
#else
    W_Y(2,ii,jj)=W_Y(2,ii,jj)-(dt)/max_ro*px
    W_Y(3,ii,jj)=W_Y(3,ii,jj)-(dt)/max_ro*py
#endif
  END DO
END DO
#endif

#ifdef CENTEREDPRIMITIVE
DO jj=1,nElemsY
  DO ii=1,nElemsX
    WC(2:3,ii,jj) = K0_WC(2:3,ii,jj) + WCt(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(WC(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(WC(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    WC(2,ii,jj)=WC(2,ii,jj)-(dt/EPS_LM**2)/max_ro*px
    WC(3,ii,jj)=WC(3,ii,jj)-(dt/EPS_LM**2)/max_ro*py
#else
    WC(2,ii,jj)=WC(2,ii,jj)-(dt)/max_ro*px
    WC(3,ii,jj)=WC(3,ii,jj)-(dt)/max_ro*py
#endif
  END DO
END DO
#endif

CALL FVTimeDerivativeConservedOnly(t+dt)


U(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXEuler_Implicit_Update_Conserved
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE TimeDiscretizationByARS222_Crank_Nicolson_Conserved(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly

USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K1_WC
USE MOD_FiniteVolume2D_vars,ONLY: K2_WC
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction     ,ONLY: Laplacian_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_Equation           ,ONLY: BoundaryConditions 
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj, indc, indi, indj, its
INTEGER         :: iterations
#ifdef CENTEREDPRIMITIVE
REAL            :: initial_guess_WC(1:NRows_WC)
#endif
REAL            :: px,py,dx,dy
REAL, PARAMETER :: g=1.0-SQRT(2.0)/2.0
REAL, PARAMETER :: d=1.0-1.0/(2.0*g)
!-------------------------------------------------------------------------------!
!*Debug
!-------------------------------------------------------------------------------!
REAL                                :: residual, dxvx, dyvy, partx, party
#ifdef CENTEREDPRIMITIVE
REAL                                :: WC0(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL                                :: WC1(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL                                :: WC2(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL                                :: WCt0(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL                                :: WCt1(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL                                :: WCtcombined(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
REAL                                :: tolerance=1e-10
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug

!-------------------------------------------------------------------------------!

WC0(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)         = 0.0
WC1(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)         = 0.0
WC2(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)         = 0.0
WCt0(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)        = 0.0
WCt1(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)        = 0.0
WCtcombined(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = 0.0


CALL BoundaryConditions(t)
#ifdef CENTEREDPRIMITIVE
WC0(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#ifdef CENTEREDPRIMITIVE
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

!*=============================================
!*STAGE 1 - Implicit Euler with g*dt time increment
!*=============================================

tStage = t 
CALL FVTimeDerivativePrimitiveOnly(tStage)

#ifdef CENTEREDPRIMITIVE
CALL Impose_BC_on_WCt()
WCt0(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
!*----------------------
!*NB: It is super important to have BCs on W_X, W_Y, Wt_X and Wt_Y
!*With this we have them
!*----------------------

!*----------------------
!*At this point
!*Primitive variables are updated as follows
!*->ro is updated normally
!*->p  is updated solving a linear system
!*->v  is updated normally but with an extra explicit contribution computed from p
!*----------------------

!*----------------------
!*Update of primitive variables <=============
!*----------------------

!*----------------------
!*->Normal update of ro <=============
!*----------------------
#ifdef CENTEREDPRIMITIVE
WC(1,1:nElemsX,1:nElemsY) = K0_WC(1,1:nElemsX,1:nElemsY) + WCt(1,1:nElemsX,1:nElemsY)*dt*g
#endif


!*=============================
#ifdef CENTEREDPRIMITIVE
!*----------------------
!*IMEX linear system WC <=============
!*----------------------
!*Assembly IMEX Matrix Euler
CALL IMEX_Matrix(dt*g)
!*Assembly IMEX rhs Euler
CALL IMEX_Euler_rhs_dt_input(WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),t+dt*g,dt*g)


!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))

!*----------------------
#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif
!*----------------------




!*Copying sol_WC into WC(4,ii,jj)
CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))

#endif



dx=MESH_DX(1)
dy=MESH_DX(2)


!*----------------------
!*Update velocity with extra explicit contribution <=============
!*----------------------
#ifdef CENTEREDPRIMITIVE
DO jj=1,nElemsY
  DO ii=1,nElemsX
    WC(2:3,ii,jj) = K0_WC(2:3,ii,jj) + WCt(2:3,ii,jj)*dt*g
    px=First_Derivative_Central_Order2(WC(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(WC(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    WC(2,ii,jj)=WC(2,ii,jj)-(dt*g/EPS_LM**2)/max_ro*px
    WC(3,ii,jj)=WC(3,ii,jj)-(dt*g/EPS_LM**2)/max_ro*py
#else
    WC(2,ii,jj)=WC(2,ii,jj)-(dt*g)/max_ro*px
    WC(3,ii,jj)=WC(3,ii,jj)-(dt*g)/max_ro*py
#endif
  END DO
END DO
#endif


CALL BoundaryConditions(t)
#ifdef CENTEREDPRIMITIVE
WC1(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif


!*=============================================
!*STAGE 2 - More complicated
!*=============================================
tStage = t+dt*g 
CALL FVTimeDerivativePrimitiveOnly(tStage)

#ifdef CENTEREDPRIMITIVE
CALL Impose_BC_on_WCt()
WCt1(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

WCtcombined(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) =       d * WCt0(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) &
                                                                        & + (1.0-d) * WCt1(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)



!*----------------------
!*->Normal update of ro <=============
!*----------------------
#ifdef CENTEREDPRIMITIVE
WC(1,1:nElemsX,1:nElemsY) = K0_WC(1,1:nElemsX,1:nElemsY) + WCtcombined(1,1:nElemsX,1:nElemsY)*dt
#endif



!*=============================
#ifdef CENTEREDPRIMITIVE
!*----------------------
!*IMEX linear system WC <=============
!*----------------------
!*Assembly IMEX rhs
CALL IMEX_ARS222_rhs(WCtcombined,WC0,WC1,t,g,d)

!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector(WC0(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))

!*----------------------
#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif
!*----------------------

!*Copying sol_WC into WC(4,ii,jj)
CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))

#endif



dx=MESH_DX(1)
dy=MESH_DX(2)


!*----------------------
!*Update velocity with extra explicit contribution <=============
!*----------------------
#ifdef CENTEREDPRIMITIVE
DO jj=1,nElemsY
  DO ii=1,nElemsX
    WC(2:3,ii,jj) = K0_WC(2:3,ii,jj) + WCtcombined(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(g*WC(4,ii-1:ii+1,jj)+(1.0-g)*WC1(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(g*WC(4,ii,jj-1:jj+1)+(1.0-g)*WC1(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    WC(2,ii,jj)=WC(2,ii,jj)-(dt/EPS_LM**2)/max_ro*px
    WC(3,ii,jj)=WC(3,ii,jj)-(dt/EPS_LM**2)/max_ro*py
#else
    WC(2,ii,jj)=WC(2,ii,jj)-(dt)/max_ro*px
    WC(3,ii,jj)=WC(3,ii,jj)-(dt)/max_ro*py
#endif
  END DO
END DO
#endif


CALL BoundaryConditions(t)
#ifdef CENTEREDPRIMITIVE
WC2(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

!*---------------------------------------------
!*Crank Nicolson for U
!*---------------------------------------------

!*====================
!*tn
!*====================
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = WC0(1:nVar,1:nElemsX,1:nElemsY)
#endif

CALL FVTimeDerivativeConservedOnly(t)

K0(1:nVar,1:nElemsX,1:nElemsY)=Ut(1:nVar,1:nElemsX,1:nElemsY)


!*====================
!*tn+1
!*====================
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = WC2(1:nVar,1:nElemsX,1:nElemsY)
#endif

CALL FVTimeDerivativeConservedOnly(t+dt)

K1(1:nVar,1:nElemsX,1:nElemsY)=Ut(1:nVar,1:nElemsX,1:nElemsY)


!*====================
!*Actual update
!*====================
U(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + 0.5*(K0(1:nVar,1:nElemsX,1:nElemsY)+K1(1:nVar,1:nElemsX,1:nElemsY))*dt

!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = WC2(1:nVar,1:nElemsX,1:nElemsY)
#endif
!*And U has been computed




!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByARS222_Crank_Nicolson_Conserved
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC2_Implicit_Update_Conserved(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 2
INTEGER, PARAMETER   :: MSteps = 2
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0, 0.5, 0.0, 0.5 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0 , 1.0  /)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

#ifdef ACTIVEFLUX
!*Initialiting initial guesses for the linear systems
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif
    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif


#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

      !*======================================
      !*NB: HERE I NEED TO CALL FV TIME DERIVATIVE OF THE PRIMITIVE VARIABLES ONLY
      !*======================================
      CALL FVTimeDerivativePrimitiveOnly(tStage(ii)) 

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif

#endif


#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif


#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif


    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))
#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

  END DO

END DO


!*====================================
!*NOW THAT PRIMTIIVE VARIABLES HAVE BEEN UPDATED 
!*WE NEED TO UPDATE THE CONSERVED ONES BASED ON THE UPDATED PRIMITIVES
!*====================================
DO ii = 1,MSteps 

#ifdef ACTIVEFLUX
  !*NB: This is in principle not necessary
  CALL Impose_BC_on_W(Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
  CALL Impose_BC_on_WC(Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

  CALL FVTimeDerivativeConservedOnly(tStage(ii))

  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

END DO


!*---------------------------------------------
!*Update last subtimenode
!*---------------------------------------------
ii=MSteps 

!*---------------------------------------------
!*Compute combination of the fluxes
!*---------------------------------------------
FLUXES_U   = 0.0

DO jj = 1,MSteps
  FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
END DO

!*->Normal update of conserved variables
Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC2_Implicit_Update_Conserved
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC3_Implicit_Update_Conserved(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(3,3) :: thetaDeC
REAL, DIMENSION(3)   :: betaDeC
REAL, DIMENSION(3)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 3
INTEGER, PARAMETER   :: MSteps = 3
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.2083333333333333333, &
                      0.1666666666666666666, &
                      0.0000000000000000000, &
                      0.3333333333333333333, &
                      0.6666666666666666666, &
                      0.0000000000000000000, &
                      -0.0416666666666666666, &
                      0.1666666666666666666 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0, 0.5 , 1.0/)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

#ifdef ACTIVEFLUX
!*Initialiting initial guesses for the linear systems
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif
    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif


#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

      !*======================================
      !*NB: HERE I NEED TO CALL FV TIME DERIVATIVE OF THE PRIMITIVE VARIABLES ONLY
      !*======================================
      CALL FVTimeDerivativePrimitiveOnly(tStage(ii)) 

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif

#endif


#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif


#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif


    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))
#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

  END DO

END DO


!*====================================
!*NOW THAT PRIMTIIVE VARIABLES HAVE BEEN UPDATED 
!*WE NEED TO UPDATE THE CONSERVED ONES BASED ON THE UPDATED PRIMITIVES
!*====================================
DO ii = 1,MSteps 

#ifdef ACTIVEFLUX
  !*NB: This is in principle not necessary
  CALL Impose_BC_on_W(Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
  CALL Impose_BC_on_WC(Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

  CALL FVTimeDerivativeConservedOnly(tStage(ii))

  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

END DO


!*---------------------------------------------
!*Update last subtimenode
!*---------------------------------------------
ii=MSteps 

!*---------------------------------------------
!*Compute combination of the fluxes
!*---------------------------------------------
FLUXES_U   = 0.0

DO jj = 1,MSteps
  FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
END DO

!*->Normal update of conserved variables
Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC3_Implicit_Update_Conserved
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC4_Implicit_Update_Conserved(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(3,3) :: thetaDeC
REAL, DIMENSION(3)   :: betaDeC
REAL, DIMENSION(3)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 4
INTEGER, PARAMETER   :: MSteps = 3
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.2083333333333333333, &
                      0.1666666666666666666, &
                      0.0000000000000000000, &
                      0.3333333333333333333, &
                      0.6666666666666666666, &
                      0.0000000000000000000, &
                      -0.0416666666666666666, &
                      0.1666666666666666666 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0, 0.5 , 1.0/)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

#ifdef ACTIVEFLUX
!*Initialiting initial guesses for the linear systems
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif
    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif


#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

      !*======================================
      !*NB: HERE I NEED TO CALL FV TIME DERIVATIVE OF THE PRIMITIVE VARIABLES ONLY
      !*======================================
      CALL FVTimeDerivativePrimitiveOnly(tStage(ii)) 

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif

#endif


#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif


#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif


    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))
#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

  END DO

END DO


!*====================================
!*NOW THAT PRIMTIIVE VARIABLES HAVE BEEN UPDATED 
!*WE NEED TO UPDATE THE CONSERVED ONES BASED ON THE UPDATED PRIMITIVES
!*====================================
DO ii = 1,MSteps 

#ifdef ACTIVEFLUX
  !*NB: This is in principle not necessary
  CALL Impose_BC_on_W(Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
  CALL Impose_BC_on_WC(Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

  CALL FVTimeDerivativeConservedOnly(tStage(ii))

  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

END DO


!*---------------------------------------------
!*Update last subtimenode
!*---------------------------------------------
ii=MSteps 

!*---------------------------------------------
!*Compute combination of the fluxes
!*---------------------------------------------
FLUXES_U   = 0.0

DO jj = 1,MSteps
  FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
END DO

!*->Normal update of conserved variables
Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC4_Implicit_Update_Conserved
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC5_Implicit_Update_Conserved(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(4,4) :: thetaDeC
REAL, DIMENSION(4)   :: betaDeC
REAL, DIMENSION(4)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 5
INTEGER, PARAMETER   :: MSteps = 4
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif

REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0000000000000000000, &
                      0.1103005664791649049, &
                      0.0730327668541684294, &
                      0.0833333333333333287, &
                      0.0000000000000000000, &
                      0.1896994335208351257, &
                      0.4505740308958107176, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      -0.0339073642291439076, &
                      0.2269672331458314485, &
                      0.4166666666666666852, &
                      0.0000000000000000000, &
                      0.0103005664791649201, &
                      -0.0269672331458315692, &
                      0.0833333333333333287 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0000000000000000000 , 0.2763932022500210639 , 0.7236067977499789361 , 1.0000000000000000000/)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

#ifdef ACTIVEFLUX
!*Initialiting initial guesses for the linear systems
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif
    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif

#endif


#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

      !*======================================
      !*NB: HERE I NEED TO CALL FV TIME DERIVATIVE OF THE PRIMITIVE VARIABLES ONLY
      !*======================================
      CALL FVTimeDerivativePrimitiveOnly(tStage(ii)) 

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif

#endif


#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif


#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif


    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------

    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))
#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
#ifdef ACTIVEFLUX
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

  END DO

END DO


!*====================================
!*NOW THAT PRIMTIIVE VARIABLES HAVE BEEN UPDATED 
!*WE NEED TO UPDATE THE CONSERVED ONES BASED ON THE UPDATED PRIMITIVES
!*====================================
DO ii = 1,MSteps 

#ifdef ACTIVEFLUX
  !*NB: This is in principle not necessary
  CALL Impose_BC_on_W(Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
  CALL Impose_BC_on_WC(Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

  CALL FVTimeDerivativeConservedOnly(tStage(ii))

  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

END DO


!*---------------------------------------------
!*Update last subtimenode
!*---------------------------------------------
ii=MSteps 

!*---------------------------------------------
!*Compute combination of the fluxes
!*---------------------------------------------
FLUXES_U   = 0.0

DO jj = 1,MSteps
  FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
END DO

!*->Normal update of conserved variables
Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)

#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC5_Implicit_Update_Conserved
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE TimeDiscretizationByIMEXDeC2_Crank_Nicolson_Conserved(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly

USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
#endif

#ifdef CENTEREDPRIMITIVE
!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
#endif

USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro

USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 2
INTEGER, PARAMETER   :: MSteps = 2
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif
REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0, 0.5, 0.0, 0.5 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0 , 1.0  /)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

#ifdef ACTIVEFLUX
!*Initialiting initial guesses for the linear systems
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif
    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif


#endif

#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==0)
      !*Stage dependent max_ro and min_p
      PRINT*, "Stage dependent max_ro and min_p"
      PRINT*, "Keep it deactivated"
      STOP
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivativePrimitiveOnly(tStage(ii))

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif
#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif


#endif

#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p
#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))
#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, MIN(50,NRows_WC-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_WC into Wa_WC(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))

#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)

#ifdef ACTIVEFLUX
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

#endif

  END DO


END DO


!*---------------------------------------------
!*Crank Nicolson for U
!*---------------------------------------------

!*====================
!*tn
!*====================
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

CALL FVTimeDerivativeConservedOnly(t)

K0(1:nVar,1:nElemsX,1:nElemsY)=Ut(1:nVar,1:nElemsX,1:nElemsY)


!*====================
!*tn+1
!*====================
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif

CALL FVTimeDerivativeConservedOnly(t+dt)

K1(1:nVar,1:nElemsX,1:nElemsY)=Ut(1:nVar,1:nElemsX,1:nElemsY)


!*====================
!*Actual update
!*====================
U(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + 0.5*(K0(1:nVar,1:nElemsX,1:nElemsY)+K1(1:nVar,1:nElemsX,1:nElemsY))*dt

!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif
!*And U has been computed

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC2_Crank_Nicolson_Conserved
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE IMEXDeC2_Crank_Nicolson_Conserved_rho_p_stage_dependent(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly

USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef ACTIVEFLUX
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
#endif

#ifdef CENTEREDPRIMITIVE
!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
#endif

USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro

USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 2
INTEGER, PARAMETER   :: MSteps = 2
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
#ifdef ACTIVEFLUX
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef ACTIVEFLUX
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
#endif
#ifdef CENTEREDPRIMITIVE
REAL                 :: FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
#endif
REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0, 0.5, 0.0, 0.5 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0 , 1.0  /)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

#ifdef ACTIVEFLUX
!*Initialiting initial guesses for the linear systems
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))
#endif
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
#ifdef ACTIVEFLUX
  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y
#endif
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
#ifdef ACTIVEFLUX
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WC(Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
#endif
    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
#ifdef ACTIVEFLUX
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

#ifdef ACTIVEFLUX
    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
        WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_WC(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_WC(1,3,indi,indj-1:indj+1),dy)

        WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO
#endif


#endif

#ifdef ACTIVEFLUX
    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==1)
      !*STAGE DEPENDENT max_ro and min_p ACTIVE
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivativePrimitiveOnly(tStage(ii))

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef ACTIVEFLUX
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif
#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          WCt(2,indi,indj) = WCt(2,indi,indj) - part_x
          WCt(3,indi,indj) = WCt(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_WC(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_WC(ii,3,indi,indj-1:indj+1),dy)

          WCt(4,indi,indj)   = WCt(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO
#endif


#endif

#ifdef ACTIVEFLUX
      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif
#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
#ifdef ACTIVEFLUX
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
#endif
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

#ifdef ACTIVEFLUX
      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
      FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of density
#ifdef ACTIVEFLUX
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
#endif
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p
#ifdef ACTIVEFLUX
    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))
#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*IMEX linear system WC <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs(FLUXES_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES WC"
    its = 10000
    sol_WC = initial_guess_WC
    call pmgmres_ilu_cr(NRows_WC, NNZsparse_WC, RowStart_WC, Columns_WC, Values_WC, sol_WC, rhs_WC, its, MIN(50,NRows_WC-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi WC"
    iterations=0
    CALL jacobi_2(RowStart_WC,Columns_WC,Values_WC,Diagonal_WC,rhs_WC,sol_WC,iterations,initial_guess_WC)

#ifdef COUNTJACOBI
    JacobiCounter_WC = JacobiCounter_WC +1
    JacobiIterations_WC(JacobiCounter_WC) = iterations
#endif

#endif

    !*Copying sol_WC into Wa_WC(ii,4,:,:)
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))

#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)

#ifdef ACTIVEFLUX
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO
#endif

#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi-1:indi+1,indj)-Wp_WC(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,nVar,indi,indj-1:indj+1)-Wp_WC(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_WC(ii,2,indi,indj)=Wa_WC(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_WC(ii,3,indi,indj)=Wa_WC(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

#endif

  END DO


END DO


!*---------------------------------------------
!*Crank Nicolson for U
!*---------------------------------------------

!*====================
!*tn
!*====================
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(1,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(1,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

CALL FVTimeDerivativeConservedOnly(t)

K0(1:nVar,1:nElemsX,1:nElemsY)=Ut(1:nVar,1:nElemsX,1:nElemsY)


!*====================
!*tn+1
!*====================
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif

CALL FVTimeDerivativeConservedOnly(t+dt)

K1(1:nVar,1:nElemsX,1:nElemsY)=Ut(1:nVar,1:nElemsX,1:nElemsY)


!*====================
!*Actual update
!*====================
U(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + 0.5*(K0(1:nVar,1:nElemsX,1:nElemsY)+K1(1:nVar,1:nElemsX,1:nElemsY))*dt

!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
#ifdef ACTIVEFLUX
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)
#endif
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif
!*And U has been computed

!-------------------------------------------------------------------------------!
END SUBROUTINE IMEXDeC2_Crank_Nicolson_Conserved_rho_p_stage_dependent
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE TimeDiscretizationByIMEXDeC2_Implicit_Euler_Conserved(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativePrimitiveOnly
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivativeConservedOnly

USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 2
INTEGER, PARAMETER   :: MSteps = 2
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0, 0.5, 0.0, 0.5 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0 , 1.0  /)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

!*Initialiting initial guesses for the linear systems
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))

DO ii = 1,MSteps

  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------

  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

#ifndef IMEXL2FULLYUPWINDED
    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px/max_ro
        part_y=py/max_ro
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

      END DO
    END DO
#endif

    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated

    DO ii = 1,MSteps

      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

#if(1==0)
      !*Stage dependent max_ro and min_p
      PRINT*, "Stage dependent max_ro and min_p"
      PRINT*, "Keep it deactivated"
      STOP
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivativePrimitiveOnly(tStage(ii))

#ifndef IMEXL2FULLYUPWINDED
      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * min_p  * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px/max_ro
          part_y=py/max_ro
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y


          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * min_p  * (ux + vy)

        END DO
      END DO
#endif


      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated

      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
    DO jj = 1,MSteps

      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of density
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

    !*->Linear system for p

    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_X(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEX_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEX_Matrix_Y(dt*betaDeC(ii)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEX_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))


    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)/max_ro*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)/max_ro*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)/max_ro*py
#endif
      END DO
    END DO


  END DO


END DO


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------

!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)

CALL FVTimeDerivativeConservedOnly(t+dt)
U(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY) + Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXDeC2_Implicit_Euler_Conserved
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE IMEX_DeC_rhs_X(Wt,W0,Wm,betam,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), INTENT(IN) :: Wt, W0, WM
REAL                                                                           , INTENT(IN) :: betam
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: supp(-1:1)
REAL                  :: dx_part, dy_part, delta_t
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)
delta_t=dt*betam !*For normalization !*IMPORTANT: Normalization w.r.t. delta_t, not dt

!*RMK: +1 in X direction

rhs_X=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)*EPS_LM)**2
#else
    rhs_X(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_X(indc)=W0(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_X(indc)=rhs_X(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=rhs_X(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_X(indc)=rhs_X(indc)-(dt**2)*betam*Gmm*min_p*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

    !*We have now two extra contributions: div v and lap p

    !*div v part
    dx_part=First_Derivative_Central_Order2(W0(2,ii-1:ii+1,jj)-Wm(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(W0(3,ii,jj-1:jj+1)-Wm(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-dt*betam*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2)
#else
    rhs_X(indc)=rhs_X(indc)-dt*betam*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2)
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=rhs_X(indc)-dt*betam*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2)
#else
    rhs_X(indc)=rhs_X(indc)-dt*betam*Gmm*min_p*(dx_part+dy_part)
#endif

    !*lap p part
    lapp=Laplacian_Central_Order2(Wm(4,ii-1:ii+1,jj-1:jj+1),dx,dy)

#ifdef NORMALIZEALL
    rhs_X(indc)=rhs_X(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*(lapp*(dx/delta_t)**2)
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-(dt/EPS_LM)**2*betam**2/max_ro*Gmm*min_p*(lapp*dx**2)
#else
    rhs_X(indc)=rhs_X(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*(lapp*dx**2)
#endif    
#else
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-(dt/EPS_LM)**2*betam**2/max_ro*Gmm*min_p*lapp
#else
    rhs_X(indc)=rhs_X(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*lapp
#endif    
#endif    

  END DO !ii
END DO !jj


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_X(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO

!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_X(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_X(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_X(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEX_DeC_rhs_X
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE IMEX_DeC_rhs_Y(Wt,W0,Wm,betam,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), INTENT(IN) :: Wt, W0, Wm
REAL                                                                           , INTENT(IN) :: betam
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: et, etx, ety !*Extra term=et=(ro*-ro)/ro/ro* divided by (in case relaxation is active) eps_lm**2
REAL                  :: supp(-1:1)
REAL                  :: div_v_grad_v
REAL                  :: ro_p_contribution
REAL                  :: dx_part, dy_part, delta_t
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)
delta_t=dt*betam !*For normalization !*IMPORTANT: Normalization w.r.t. delta_t, not dt

!*RMK: +1 in Y direction

rhs_Y=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)*EPS_LM)**2
#else
    rhs_Y(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_Y(indc)=W0(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1


    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_Y(indc)=rhs_Y(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=rhs_Y(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_Y(indc)=rhs_Y(indc)-(dt**2)*betam*Gmm*min_p*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

    !*We have now two extra contributions: div v and lap p

    !*div v part
    dx_part=First_Derivative_Central_Order2(W0(2,ii-1:ii+1,jj)-Wm(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(W0(3,ii,jj-1:jj+1)-Wm(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-dt*betam*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2)
#else
    rhs_Y(indc)=rhs_Y(indc)-dt*betam*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2)
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=rhs_Y(indc)-dt*betam*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2)
#else
    rhs_Y(indc)=rhs_Y(indc)-dt*betam*Gmm*min_p*(dx_part+dy_part)
#endif

    !*lap p part
    lapp=Laplacian_Central_Order2(Wm(4,ii-1:ii+1,jj-1:jj+1),dx,dy)

#ifdef NORMALIZEALL
    rhs_Y(indc)=rhs_Y(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*(lapp*(dx/delta_t)**2)
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-(dt/EPS_LM)**2*betam**2/max_ro*Gmm*min_p*(lapp*dx**2)
#else
    rhs_Y(indc)=rhs_Y(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*(lapp*dx**2)
#endif    
#else
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-(dt/EPS_LM)**2*betam**2/max_ro*Gmm*min_p*lapp
#else
    rhs_Y(indc)=rhs_Y(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*lapp
#endif    
#endif

  END DO !ii
END DO !jj





!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_Y(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_Y(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_Y(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_Y(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEX_DeC_rhs_Y
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE IMEX_DeC_rhs(Wt,W0,Wm,betam,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), INTENT(IN) :: Wt, W0, WM
REAL                                                                           , INTENT(IN) :: betam
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: supp(-1:1)
REAL                  :: dx_part, dy_part, delta_t
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)
delta_t=dt*betam !*For normalization !*IMPORTANT: Normalization w.r.t. delta_t, not dt

!*RMK: +1 in X direction

rhs_WC=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)*EPS_LM)**2
#else
    rhs_WC(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_WC(indc)=W0(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_WC(indc)=rhs_WC(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=rhs_WC(indc)-(dt**2)*betam*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_WC(indc)=rhs_WC(indc)-(dt**2)*betam*Gmm*min_p*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

    !*We have now two extra contributions: div v and lap p

    !*div v part
    dx_part=First_Derivative_Central_Order2(W0(2,ii-1:ii+1,jj)-Wm(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(W0(3,ii,jj-1:jj+1)-Wm(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)-dt*betam*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2)
#else
    rhs_WC(indc)=rhs_WC(indc)-dt*betam*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2)
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=rhs_WC(indc)-dt*betam*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2)
#else
    rhs_WC(indc)=rhs_WC(indc)-dt*betam*Gmm*min_p*(dx_part+dy_part)
#endif

    !*lap p part
    lapp=Laplacian_Central_Order2(Wm(4,ii-1:ii+1,jj-1:jj+1),dx,dy)

#ifdef NORMALIZEALL
    rhs_WC(indc)=rhs_WC(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*(lapp*(dx/delta_t)**2)
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)-(dt/EPS_LM)**2*betam**2/max_ro*Gmm*min_p*(lapp*dx**2)
#else
    rhs_WC(indc)=rhs_WC(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*(lapp*dx**2)
#endif    
#else
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)-(dt/EPS_LM)**2*betam**2/max_ro*Gmm*min_p*lapp
#else
    rhs_WC(indc)=rhs_WC(indc)-(dt       )**2*betam**2/max_ro*Gmm*min_p*lapp
#endif    
#endif    

  END DO !ii
END DO !jj


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEX_DeC_rhs
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE IMEX_ARS222_rhs(Wt,W0,Wm,t,g,d)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), INTENT(IN) :: Wt, W0, WM
REAL                                                                           , INTENT(IN) :: t
REAL                                                                           , INTENT(IN) :: g
REAL                                                                           , INTENT(IN) :: d
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: supp(-1:1)
REAL                  :: dx_part, dy_part, delta_t
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)
delta_t=dt*g !*For normalization !*IMPORTANT: Normalization w.r.t. delta_t, not dt

rhs_WC=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)*EPS_LM)**2
#else
    rhs_WC(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_WC(indc)=W0(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)-(dt**2)*g*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_WC(indc)=rhs_WC(indc)-(dt**2)*g*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=rhs_WC(indc)-(dt**2)*g*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_WC(indc)=rhs_WC(indc)-(dt**2)*g*Gmm*min_p*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

    !*We have now two extra contributions: div v and lap p

    !*div v part
    dx_part=First_Derivative_Central_Order2(g*W0(2,ii-1:ii+1,jj)+(1.0-g)*Wm(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(g*W0(3,ii,jj-1:jj+1)+(1.0-g)*Wm(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)-dt*Gmm*min_p*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2)
#else
    rhs_WC(indc)=rhs_WC(indc)-dt*Gmm*min_p*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2)
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_WC(indc)=rhs_WC(indc)-dt*Gmm*min_p*(dx_part*dx**2+dy_part*dx**2)
#else
    rhs_WC(indc)=rhs_WC(indc)-dt*Gmm*min_p*(dx_part+dy_part)
#endif

    !*lap p part
    lapp=Laplacian_Central_Order2(Wm(4,ii-1:ii+1,jj-1:jj+1),dx,dy)

#ifdef NORMALIZEALL
    rhs_WC(indc)=rhs_WC(indc)+(dt       )**2*(1.0-g)*g/max_ro*Gmm*min_p*(lapp*(dx/delta_t)**2)
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)+(dt/EPS_LM)**2*(1.0-g)*g/max_ro*Gmm*min_p*(lapp*dx**2)
#else
    rhs_WC(indc)=rhs_WC(indc)+(dt       )**2*(1.0-g)*g/max_ro*Gmm*min_p*(lapp*dx**2)
#endif    
#else
#ifdef RELAXATION
    rhs_WC(indc)=rhs_WC(indc)+(dt/EPS_LM)**2*(1.0-g)*g/max_ro*Gmm*min_p*lapp
#else
    rhs_WC(indc)=rhs_WC(indc)+(dt       )**2*(1.0-g)*g/max_ro*Gmm*min_p*lapp
#endif    
#endif    

  END DO !ii
END DO !jj


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_WC(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEX_ARS222_rhs
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#endif
#ifdef IMEXMOMENTUM
SUBROUTINE TimeDiscretizationByIMEXMOMENTUMEuler(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX

!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: K0_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: K0_Y

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y


USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction     ,ONLY: Laplacian_Central_Order2
USE MOD_Equation           ,ONLY: BoundaryConditions 
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y

USE MOD_Equation           ,ONLY: ExactFunction
USE MOD_Equation           ,ONLY: ConsToPrim
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_X
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: tStage, final_time, x(2), Cons(1:4), Prim(1:4)
INTEGER         :: ii, jj, indc, indi, indj, its
INTEGER         :: iterations
REAL            :: initial_guess_X(1:NRows_X)
REAL            :: initial_guess_Y(1:NRows_Y)
REAL            :: px,py,dx,dy
!-------------------------------------------------------------------------------!
!*Debug
!-------------------------------------------------------------------------------!
REAL                                :: residual, dxvx, dyvy, partx, party
REAL                                :: old_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
REAL                                :: old_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
REAL                                :: tolerance=1e-10
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug

!-------------------------------------------------------------------------------!

#if(1==1)
!*=============================================
!*=============================================
!*=============================================
!*FOR SAFETY CHECK
CALL BoundaryConditions(t)
old_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)=W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
old_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)=W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
!*=============================================
!*=============================================
!*=============================================
#endif

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
K0_X(1:nVar,1:nElemsX+1,1:nElemsY) = W_X(1:nVar,1:nElemsX+1,1:nElemsY) !*NB: +1 in X direction
K0_Y(1:nVar,1:nElemsX,1:nElemsY+1) = W_Y(1:nVar,1:nElemsX,1:nElemsY+1) !*NB: +1 in Y direction

tStage = t 
CALL FVTimeDerivative(tStage)

CALL Impose_BC_on_Wt()
!*----------------------
!*NB: It is super important to have BCs on W_X, W_Y, Wt_X and Wt_Y
!*With this we have them
!*----------------------

!*----------------------
!*At this point
!*Conservative variables are updated normally
!*Primitive variables are updated differently
!*->ro is updated normally
!*->p  is updated solving a linear system
!*->rov  is updated normally but with an extra explicit contribution computed from p
!*----------------------

!*----------------------
!*Normal update of conservative variables <=============
!*----------------------
U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + Ut(1:nVar,1:nElemsX,1:nElemsY)*dt


!*----------------------
!*Update of primitive variables <=============
!*----------------------

!*----------------------
!*->Normal update of ro <=============
!*----------------------
W_X(1,1:nElemsX+1,1:nElemsY) = K0_X(1,1:nElemsX+1,1:nElemsY) + Wt_X(1,1:nElemsX+1,1:nElemsY)*dt !*NB:+1 in X direction
W_Y(1,1:nElemsX,1:nElemsY+1) = K0_Y(1,1:nElemsX,1:nElemsY+1) + Wt_Y(1,1:nElemsX,1:nElemsY+1)*dt !*NB:+1 in Y direction


!*----------------------
!*IMEX linear system W_X <=============
!*----------------------
!*Assembly IMEX Matrix Euler
CALL IMEXMOMENTUM_Matrix_X(dt,W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))
!*Assembly IMEX rhs Euler
CALL IMEXMOMENTUM_Euler_rhs_X(Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),t+dt)

!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, 1000, 1e-12, 1e-6 )
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif



#if(1==0)
ALLOCATE(sol_debug(1:NRows_X))
iterations=0
CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_debug,iterations,initial_guess_X)

#ifdef COUNTJACOBI
JacobiCounter_X = JacobiCounter_X +1
JacobiIterations_X(JacobiCounter_X) = iterations
#endif

DO indi=1,NRows_X
  IF(ABS(sol_X(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_X"
    PRINT*, indi, sol_X(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug)
#endif



!*Copying sol_X into W_X(4,ii,jj)
CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))


!*----------------------
!*IMEX linear system W_Y <=============
!*----------------------
!*Assembly IMEX Matrix Euler
CALL IMEXMOMENTUM_Matrix_Y(dt,W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))
!*Assembly IMEX rhs Euler
CALL IMEXMOMENTUM_Euler_rhs_Y(Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),t+dt)

!*Jacobi
iterations=0

!*Initialiting initial guess
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))


#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, 1000, 1e-12, 1e-6 )
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif




#if(1==0)
ALLOCATE(sol_debug(1:NRows_Y))
iterations=0
CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_debug,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
JacobiCounter_Y = JacobiCounter_Y +1
JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

DO indi=1,NRows_Y
  IF(ABS(sol_Y(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_Y"
    PRINT*, indi, sol_Y(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug)
#endif

!*Copying sol_Y into W_Y(4,ii,jj)
CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))


dx=MESH_DX(1)
dy=MESH_DX(2)
!*----------------------
!*Update momentum with extra explicit contribution <=============
!*----------------------
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    W_X(2:3,ii,jj) = K0_X(2:3,ii,jj) + Wt_X(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(W_X(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(W_X(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    W_X(2,ii,jj)=W_X(2,ii,jj)-(dt/EPS_LM**2)*px
    W_X(3,ii,jj)=W_X(3,ii,jj)-(dt/EPS_LM**2)*py
#else
    W_X(2,ii,jj)=W_X(2,ii,jj)-(dt)*px
    W_X(3,ii,jj)=W_X(3,ii,jj)-(dt)*py
#endif
  END DO
END DO

!*----------------------
!*Update momentum with extra explicit contribution <=============
!*----------------------
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    W_Y(2:3,ii,jj) = K0_Y(2:3,ii,jj) + Wt_Y(2:3,ii,jj)*dt
    px=First_Derivative_Central_Order2(W_Y(4,ii-1:ii+1,jj),dx)
    py=First_Derivative_Central_Order2(W_Y(4,ii,jj-1:jj+1),dy)
#ifdef RELAXATION
    W_Y(2,ii,jj)=W_Y(2,ii,jj)-(dt/EPS_LM**2)*px
    W_Y(3,ii,jj)=W_Y(3,ii,jj)-(dt/EPS_LM**2)*py
#else
    W_Y(2,ii,jj)=W_Y(2,ii,jj)-(dt)*px
    W_Y(3,ii,jj)=W_Y(3,ii,jj)-(dt)*py
#endif
  END DO
END DO

!*==========================================================
! final_time=t+dt
! DO jj=1,nElemsY
!   DO ii=1,nElemsX+1
!     x(1:2)=MeshBary_X(1:2,ii,jj)
!     CALL ExactFunction(InitialCondition,final_time,x,Cons)
!     CALL ConsToPrim(Cons,Prim)
!     ! W_X(1,ii,jj)=Prim(1)
!     ! W_X(2,ii,jj)=Prim(2)
!     ! W_X(3,ii,jj)=Prim(3)
!     ! W_X(4,ii,jj)=Prim(4)
!   END DO
! END DO
!*==========================================================

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXMOMENTUMEuler
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE IMEXMOMENTUM_Matrix_X(delta_t,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts

USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_X
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_X
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_X

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: delta_t
REAL,INTENT(IN) :: W( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER         :: ii, jj
INTEGER         :: inde, indc, ind, indB, indR, indU, indL
REAL            :: dx, dy
REAL            :: Const_dx, Const_dy
!-------------------------------------------------------------------------------!

dx = MESH_DX(1)
dy = MESH_DX(2)




!*First equations->Inside the domain

indc=0 !*index cell
inde=0 !*index element

DO jj=1,nElemsY
  DO ii=1,nElemsX+1 

#ifdef NORMALIZEALL
    Const_dx =                      Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy =           (dx/dy)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
    Const_dx = (delta_t/EPS_LM        )**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy = (delta_t/EPS_LM*(dx/dy))**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#else
    Const_dx = (delta_t        )**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy = (delta_t*(dx/dy))**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#endif
#else
#ifdef RELAXATION
    Const_dx = (delta_t/EPS_LM/dx)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy = (delta_t/EPS_LM/dy)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#else
    Const_dx = (delta_t/dx)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy = (delta_t/dy)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#endif
#endif



    !*C,B,R,U,L
    indc=indc+1
    inde=inde+1

    RowStart_X(indc)=inde

    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF (((From_1Index_To_2Indices_X(indc,1)-ii) .NE. 0) .OR. ((From_1Index_To_2Indices_X(indc,2)-jj) .NE. 0)) THEN
    !   PRINT*, "Problem in numeration for linear system IMEX", indc
    !   PRINT*, ii,jj
    !   PRINT*, From_1Index_To_2Indices_X(indc,:)
    !   STOP
    ! END IF
    ! PRINT*, ii,jj, From_1Index_To_2Indices_X(indc,:)
    !*-----------------------------------------

    !*NB: The 1-index label of ii,jj is indc
    !*We do not need to retrieve this information as we are sweeping the cells in the same order we numerated them 

    !*C
    ! ind=Neighbours_1Index_X(indc,1) !*NB: This is equal to indc
    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF ((indc-ind) .NE. 0) THEN
    !   PRINT*, "Problem"
    !   PRINT*, ind-indc
    !   STOP
    ! END IF
    !*-----------------------------------------
    Columns_X(inde)  = indc
#ifdef NORMALIZEALL
#ifdef RELAXATION
    Values_X(inde)   = ((dx/delta_t)*EPS_LM)**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_X(inde)   = ((dx/delta_t)       )**2+2.0*Const_dx+2.0*Const_dy
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    Values_X(inde)   = dx**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_X(inde)   = 1.0+2.0*Const_dx+2.0*Const_dy
#endif
    Diagonal_X(indc) = Values_X(inde)

    !*B (ii,jj-1)
    inde=inde+1
    indB=Neighbours_1Index_X(indc,2) !*From_2Indices_To_1Index_X(ii,jj-1)
    Columns_X(inde) = indB
    Values_X(inde)  = -Const_dy

    !*R (ii+1,jj)
    inde=inde+1
    indR=Neighbours_1Index_X(indc,3) !*From_2Indices_To_1Index_X(ii+1,jj)
    Columns_X(inde) = indR
    Values_X(inde)  = -Const_dx

    !*U (ii,jj+1)
    inde=inde+1
    indU=Neighbours_1Index_X(indc,4) !*From_2Indices_To_1Index_X(ii,jj+1)
    Columns_X(inde) = indU
    Values_X(inde)  = -Const_dy

    !*L (ii-1,jj)
    inde=inde+1
    indL=Neighbours_1Index_X(indc,5) !*From_2Indices_To_1Index_X(ii-1,jj)
    Columns_X(inde) = indL
    Values_X(inde)  = -Const_dx

    !*-----------------------------------------
    !*SAFETY CHECK
    !*-----------------------------------------
    ! PRINT*
    ! PRINT*, indc 
    ! PRINT*, "C", indc, Neighbours_1Index_X(indc,1)
    ! PRINT*, "B", indB, Neighbours_1Index_X(indc,2)
    ! PRINT*, "R", indR, Neighbours_1Index_X(indc,3)
    ! PRINT*, "U", indU, Neighbours_1Index_X(indc,4)
    ! PRINT*, "L", indL, Neighbours_1Index_X(indc,5)
    !*-----------------------------------------

  END DO
END DO


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1
  inde=inde+1
  RowStart_X(indc) = inde
  Columns_X(inde)  = indc
  Values_X(inde)   = 1.0
  Diagonal_X(indc) = Values_X(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(ii,nElemsY)
      Values_X(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(ii,1)
      Values_X(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  inde=inde+1
  RowStart_X(indc) = inde
  Columns_X(inde)  = indc
  Values_X(inde)   = 1.0
  Diagonal_X(indc) = Values_X(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1) ! Periodic u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(2,jj) !*Tricky 1=Nx+1
      Values_X(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(nElemsX,jj) !*ALSO TRICKY, SIMMETRY W.R.T. Nx+1
      Values_X(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1
  inde=inde+1
  RowStart_X(indc) = inde
  Columns_X(inde)  = indc
  Values_X(inde)   = 1.0
  Diagonal_X(indc) = Values_X(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(ii,1)
      Values_X(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(ii,nElemsY)
      Values_X(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  inde=inde+1
  RowStart_X(indc) = inde
  Columns_X(inde)  = indc
  Values_X(inde)   = 1.0
  Diagonal_X(indc) = Values_X(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(nElemsX,jj) !*Tricky 1=Nx+1
      Values_X(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_X(inde) = From_2Indices_To_1Index_X(2,jj) !*ALSO TRICKY, SIMMETRY W.R.T. 1
      Values_X(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO


RowStart_X(NRows_X+1)=NNZsparse_X+1

!*-----------------------------------------
!*SAFETY CHECK
!*-----------------------------------------
! PRINT*, Columns_X
! PRINT*, Values_X
! PRINT*, RowStart_X
! PRINT*, Diagonal_X
! STOP
!*-----------------------------------------

END SUBROUTINE IMEXMOMENTUM_Matrix_X
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE IMEXMOMENTUM_Matrix_Y(delta_t,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts

USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_Y
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_Y

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: delta_t
REAL,INTENT(IN) :: W( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER         :: ii, jj
INTEGER         :: inde, indc, ind, indB, indR, indU, indL
REAL            :: dx, dy
REAL            :: Const_dx, Const_dy
!-------------------------------------------------------------------------------!

dx = MESH_DX(1)
dy = MESH_DX(2)



!*First equations->Inside the domain

indc=0 !*index cell
inde=0 !*index element

DO jj=1,nElemsY+1
  DO ii=1,nElemsX 

#ifdef NORMALIZEALL
    Const_dx =                      Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy =           (dx/dy)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
    Const_dx = (delta_t/EPS_LM        )**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy = (delta_t/EPS_LM*(dx/dy))**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#else
    Const_dx = (delta_t        )**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy = (delta_t*(dx/dy))**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#endif
#else
#ifdef RELAXATION
    Const_dx = (delta_t/EPS_LM/dx)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy = (delta_t/EPS_LM/dy)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#else
    Const_dx = (delta_t/dx)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
    Const_dy = (delta_t/dy)**2*Gmm*W(4,ii,jj)/W(1,ii,jj)
#endif
#endif


    !*C,B,R,U,L
    indc=indc+1
    inde=inde+1

    RowStart_Y(indc)=inde

    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF (((From_1Index_To_2Indices_Y(indc,1)-ii) .NE. 0) .OR. ((From_1Index_To_2Indices_Y(indc,2)-jj) .NE. 0)) THEN
    !   PRINT*, "Problem in numeration for linear system IMEX", indc
    !   PRINT*, ii,jj
    !   PRINT*, From_1Index_To_2Indices_Y(indc,:)
    !   STOP
    ! END IF
    ! PRINT*, ii,jj, From_1Index_To_2Indices_Y(indc,:)
    !*-----------------------------------------

    !*NB: The 1-index label of ii,jj is indc
    !*We do not need to retrieve this information as we are sweeping the cells in the same order we numerated them 

    !*C
    ind=Neighbours_1Index_Y(indc,1) !*NB: This is equal to indc
    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF ((indc-ind) .NE. 0) THEN
    !   PRINT*, "Problem"
    !   PRINT*, ind-indc
    !   STOP
    ! END IF
    !*-----------------------------------------
    Columns_Y(inde)  = indc
#ifdef NORMALIZEALL
#ifdef RELAXATION
    Values_Y(inde)   = ((dx/delta_t)*EPS_LM)**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_Y(inde)   = ((dx/delta_t)       )**2+2.0*Const_dx+2.0*Const_dy
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    Values_Y(inde)   = dx**2+2.0*Const_dx+2.0*Const_dy
#else
    Values_Y(inde)   = 1.0+2.0*Const_dx+2.0*Const_dy
#endif
    Diagonal_Y(indc) = Values_Y(inde)

    !*B (ii,jj-1)
    inde=inde+1
    indB=Neighbours_1Index_Y(indc,2) !*From_2Indices_To_1Index_Y(ii,jj-1)
    Columns_Y(inde) = indB
    Values_Y(inde)  = -Const_dy

    !*R (ii+1,jj)
    inde=inde+1
    indR=Neighbours_1Index_Y(indc,3) !*From_2Indices_To_1Index_Y(ii+1,jj)
    Columns_Y(inde) = indR
    Values_Y(inde)  = -Const_dx

    !*U (ii,jj+1)
    inde=inde+1
    indU=Neighbours_1Index_Y(indc,4) !*From_2Indices_To_1Index_Y(ii,jj+1)
    Columns_Y(inde) = indU
    Values_Y(inde)  = -Const_dy

    !*L (ii-1,jj)
    inde=inde+1
    indL=Neighbours_1Index_Y(indc,5) !*From_2Indices_To_1Index_Y(ii-1,jj)
    Columns_Y(inde) = indL
    Values_Y(inde)  = -Const_dx

    !*-----------------------------------------
    !*SAFETY CHECK
    !*-----------------------------------------
    ! PRINT*
    ! PRINT*, indc 
    ! PRINT*, "C", indc, Neighbours_1Index_Y(indc,1)
    ! PRINT*, "B", indB, Neighbours_1Index_Y(indc,2)
    ! PRINT*, "R", indR, Neighbours_1Index_Y(indc,3)
    ! PRINT*, "U", indU, Neighbours_1Index_Y(indc,4)
    ! PRINT*, "L", indL, Neighbours_1Index_Y(indc,5)
    !*-----------------------------------------

  END DO
END DO


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX
  indc=indc+1
  inde=inde+1
  RowStart_Y(indc) = inde
  Columns_Y(inde)  = indc
  Values_Y(inde)   = 1.0
  Diagonal_Y(indc) = Values_Y(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(ii,nElemsY) !*Tricky 1=Ny+1
      Values_Y(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(ii,2) !*ALSO TRICKY, SIMMETRY W.R.T. 1
      Values_Y(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  inde=inde+1
  RowStart_Y(indc) = inde
  Columns_Y(inde)  = indc
  Values_Y(inde)   = 1.0
  Diagonal_Y(indc) = Values_Y(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1) ! Periodic u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(1,jj) 
      Values_Y(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(nElemsX,jj)
      Values_Y(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1
  inde=inde+1
  RowStart_Y(indc) = inde
  Columns_Y(inde)  = indc
  Values_Y(inde)   = 1.0
  Diagonal_Y(indc) = Values_Y(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(ii,2) !*ALSO TRICKY, SIMMETRY W.R.T. 1
      Values_Y(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(ii,nElemsY) !*Tricky 1=Ny+1
      Values_Y(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1
  inde=inde+1
  RowStart_Y(indc) = inde
  Columns_Y(inde)  = indc
  Values_Y(inde)   = 1.0
  Diagonal_Y(indc) = Values_Y(inde)

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1)   ! Periodic u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(nElemsX,jj)
      Values_Y(inde)  = -1.0
    CASE(2,5) ! Transmissive,Reflecting u=u
      inde=inde+1
      Columns_Y(inde) = From_2Indices_To_1Index_Y(1,jj)
      Values_Y(inde)  = -1.0
    CASE(3,4) ! Inflow,Outflow u=something
      !*Nothing
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO


RowStart_Y(NRows_Y+1)=NNZsparse_Y+1

!*-----------------------------------------
!*SAFETY CHECK
!*-----------------------------------------
! PRINT*, Columns_Y
! PRINT*, Values_Y
! PRINT*, RowStart_Y
! PRINT*, Diagonal_Y
! STOP
!*-----------------------------------------

END SUBROUTINE IMEXMOMENTUM_Matrix_Y
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE IMEXMOMENTUM_Euler_rhs_X(Wt,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), INTENT(IN) :: Wt
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: et, etx, ety !*Extra term=et=(ro*-ro)/ro/ro* divided by (in case relaxation is active) eps_lm**2
REAL                  :: supp(-1:1)
REAL                  :: div_v_grad_v
REAL                  :: ro_p_contribution
REAL                  :: dx_part, dy_part, div_vn
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)

!*RMK: +1 in X direction

rhs_X=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=(W_X(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)*EPS_LM)**2
#else
    rhs_X(indc)=(W_X(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=(W_X(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_X(indc)=W_X(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    !*NB: It is important to have BCs here
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

    ux=First_Derivative_Central_Order2(W_X(2,ii-1:ii+1,jj),dx)
    vy=First_Derivative_Central_Order2(W_X(3,ii,jj-1:jj+1),dy)

    div_vn=ux+vy !*NB: I DO NOT CHANGE THE NAME BUT IT IS THE DIVERGENCE OF THE MOMENTUM BECAUS I ONLY RUN THIS WITH MOMENTUM IN PRIMITIVE VARIABLES

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-dt*Gmm*W_X(4,ii,jj)/W_X(1,ii,jj)*(div_vn*((dx/dt)*EPS_LM)**2)                 &
                          &-(dt**2)*Gmm*W_X(4,ii,jj)/W_X(1,ii,jj)*(dx_part*((dx/dt)*EPS_LM)**2+dy_part*((dx/dt)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_X(indc)=rhs_X(indc)-dt*Gmm*W_X(4,ii,jj)/W_X(1,ii,jj)*(div_vn*((dx/dt)       )**2)                 &
                          &-(dt**2)*Gmm*W_X(4,ii,jj)/W_X(1,ii,jj)*(dx_part*((dx/dt)       )**2+dy_part*((dx/dt)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=rhs_X(indc)-dt*Gmm*W_X(4,ii,jj)/W_X(1,ii,jj)*(div_vn*dx**2)                 &
                          &-(dt**2)*Gmm*W_X(4,ii,jj)/W_X(1,ii,jj)*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_X(indc)=rhs_X(indc)-dt*Gmm*W_X(4,ii,jj)/W_X(1,ii,jj)*div_vn                 &
                          &-(dt**2)*Gmm*W_X(4,ii,jj)/W_X(1,ii,jj)*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

  END DO !ii
END DO !jj





!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEXMOMENTUM_Euler_rhs_X
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE IMEXMOMENTUM_Euler_rhs_Y(Wt,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), INTENT(IN) :: Wt
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: et, etx, ety !*Extra term=et=(ro*-ro)/ro/ro* divided by (in case relaxation is active) eps_lm**2
REAL                  :: supp(-1:1)
REAL                  :: div_v_grad_v
REAL                  :: ro_p_contribution
REAL                  :: dx_part, dy_part, div_vn
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)

!*RMK: +1 in Y direction

rhs_Y=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=(W_Y(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)*EPS_LM)**2
#else
    rhs_Y(indc)=(W_Y(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/dt)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=(W_Y(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_Y(indc)=W_Y(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc


    !*Alex's suggestion with upwinded (already computed operators)
    !*NB: It is important to have BCs here
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

    ux=First_Derivative_Central_Order2(W_Y(2,ii-1:ii+1,jj),dx)
    vy=First_Derivative_Central_Order2(W_Y(3,ii,jj-1:jj+1),dy)

    div_vn=ux+vy !*NB: I DO NOT CHANGE THE NAME BUT IT IS THE DIVERGENCE OF THE MOMENTUM BECAUS I ONLY RUN THIS WITH MOMENTUM IN PRIMITIVE VARIABLES

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-dt*Gmm*W_Y(4,ii,jj)/W_Y(1,ii,jj)*(div_vn*((dx/dt)*EPS_LM)**2)                 &
                          &-(dt**2)*Gmm*W_Y(4,ii,jj)/W_Y(1,ii,jj)*(dx_part*((dx/dt)*EPS_LM)**2+dy_part*((dx/dt)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_Y(indc)=rhs_Y(indc)-dt*Gmm*W_Y(4,ii,jj)/W_Y(1,ii,jj)*(div_vn*((dx/dt)       )**2)                 &
                          &-(dt**2)*Gmm*W_Y(4,ii,jj)/W_Y(1,ii,jj)*(dx_part*((dx/dt)       )**2+dy_part*((dx/dt)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=rhs_Y(indc)-dt*Gmm*W_Y(4,ii,jj)/W_Y(1,ii,jj)*(div_vn*dx**2)                 &
                          &-(dt**2)*Gmm*W_Y(4,ii,jj)/W_Y(1,ii,jj)*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_Y(indc)=rhs_Y(indc)-dt*Gmm*W_Y(4,ii,jj)/W_Y(1,ii,jj)*div_vn                 &
                          &-(dt**2)*Gmm*W_Y(4,ii,jj)/W_Y(1,ii,jj)*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

  END DO !ii
END DO !jj





!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEXMOMENTUM_Euler_rhs_Y
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationByIMEXMOMENTUMDeC2(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: FUp !F(U^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Ua  !U^(k)
USE MOD_FiniteVolume2D_vars,ONLY: Up  !U^(k-1)

!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X !W_X^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X  !F(W_X^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X  !W_X^(k)

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y !W_Y^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y  !F(W_Y^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y  !W_Y^(k)

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions
USE MOD_Equation           ,ONLY: Impose_BC_on_Wt
USE MOD_Equation           ,ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_X
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector_Y
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix_Y

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr

!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y

USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem

USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: ConsToPrim

USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif


!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(2,2) :: thetaDeC
REAL, DIMENSION(2)   :: betaDeC
REAL, DIMENSION(2)   :: tStage
INTEGER, PARAMETER   :: KCorr  = 2
INTEGER, PARAMETER   :: MSteps = 2
INTEGER              :: ii, jj, kk, indc, indi, indj, its
INTEGER              :: iterations
REAL                 :: initial_guess_X(1:NRows_X)
REAL                 :: initial_guess_Y(1:NRows_Y)
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
REAL                 :: FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Defined to include BCs
REAL                 :: FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) !*NB: Defined to include BCs
REAL, DIMENSION(nVar):: Wtemp
REAL, DIMENSION(nVar):: Utemp

!-------------------------------------------------------------------------------!
thetaDeC = reshape((/ 0.0, 0.5, 0.0, 0.5 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0 , 1.0  /)
tStage = t + betaDeC*dt

dx=MESH_DX(1)
dy=MESH_DX(2)

!*---------------------------------------------
!*Initialization
!*---------------------------------------------
!*NB: We start by imposing BCs so that they are imposed on U, W_X and W_Y and we can copy them
!*Essentially we need them on W_X and W_Y
!*---------------------------------------------
CALL BoundaryConditions(t)

!*Initialiting initial guesses for the linear systems
CALL Put_Matrix_In_Vector_X(W_X(4,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_X(1:NRows_X))
CALL Put_Matrix_In_Vector_Y(W_Y(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),initial_guess_Y(1:NRows_Y))

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

  !*For staggering in X direction
  !*NB:+1 in X direction
  Wa_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
  Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

  !*For staggering in X direction
  !*NB:+1 in Y direction
  Wa_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
  Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
END DO

!*---------------------------------------------
!*Iterations
!*---------------------------------------------
DO kk = 1,KCorr

  !*---------------------------------------------
  !*Set p=a so that we'll compute p->a
  !*---------------------------------------------
  Up=Ua

  !*For staggering in X direction
  !*NB:+1 in X direction
  Wp_X=Wa_X

  !*For staggering in Y direction
  !*NB:+1 in Y direction
  Wp_Y=Wa_Y

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
      CALL Impose_BC_on_W(Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),tStage(ii))
    END DO
  END IF

  !*---------------------------------------------
  !*Computation of the fluxes
  !*---------------------------------------------
  IF (kk == 1) THEN !*First iteration, all fluxes are equal to the one in first subtimenode
    !*==============================
    !*IMPORTANT: DO NOT TOUCH U, W_X AND W_Y BEFORE THIS CALL
    !*==============================
    !*SAFER
    U(1:nVar,1:nElemsX,1:nElemsY)                                       = Up(1,1:nVar,1:nElemsX,1:nElemsY)
    W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------

    DO indj=1,nElemsY  
      DO indi=1,nElemsX+1
        px=First_Derivative_Central_Order2(Wp_X(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_X(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px
        part_y=py
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif
        Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
        Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

        !*NB: This is actually the divergence of the momentum because I run with momentum in conserved variables only
        ux=First_Derivative_Central_Order2(Wp_X(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_X(1,3,indi,indj-1:indj+1),dy)

        Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * Wp_X(1,4,indi,indj)/Wp_X(1,1,indi,indj) * (ux + vy)
      END DO
    END DO

    DO indj=1,nElemsY+1  
      DO indi=1,nElemsX

        px=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_Y(1,nVar,indi,indj-1:indj+1),dy)
        part_x=px
        part_y=py
#ifdef RELAXATION
        part_x=part_x/EPS_LM**2
        part_y=part_y/EPS_LM**2
#endif

        Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
        Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y

        !*NB: This is actually the divergence of the momentum because I run with momentum in conserved variables only
        ux=First_Derivative_Central_Order2(Wp_Y(1,2,indi-1:indi+1,indj),dx)
        vy=First_Derivative_Central_Order2(Wp_Y(1,3,indi,indj-1:indj+1),dy)

        Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * Wp_Y(1,4,indi,indj)/Wp_Y(1,1,indi,indj) * (ux + vy)

      END DO
    END DO

    CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in X direction
      W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)


      CALL FVTimeDerivative(tStage(ii))

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
      DO indj=1,nElemsY  
        DO indi=1,nElemsX+1
          px=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px
          part_y=py
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif
          Wt_X(2,indi,indj) = Wt_X(2,indi,indj) - part_x
          Wt_X(3,indi,indj) = Wt_X(3,indi,indj) - part_y

          !*NB: This is actually the divergence of the momentum because I run with momentum in conserved variables only
          ux=First_Derivative_Central_Order2(Wp_X(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_X(ii,3,indi,indj-1:indj+1),dy)

          Wt_X(4,indi,indj)   = Wt_X(4,indi,indj) - Gmm * Wp_X(ii,4,indi,indj)/Wp_X(ii,1,indi,indj) * (ux + vy)
        END DO
      END DO

      DO indj=1,nElemsY+1  
        DO indi=1,nElemsX

          px=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
          part_x=px
          part_y=py
#ifdef RELAXATION
          part_x=part_x/EPS_LM**2
          part_y=part_y/EPS_LM**2
#endif

          Wt_Y(2,indi,indj) = Wt_Y(2,indi,indj) - part_x
          Wt_Y(3,indi,indj) = Wt_Y(3,indi,indj) - part_y

          !*NB: This is actually the divergence of the momentum because I run with momentum in conserved variables only
          ux=First_Derivative_Central_Order2(Wp_Y(ii,2,indi-1:indi+1,indj),dx)
          vy=First_Derivative_Central_Order2(Wp_Y(ii,3,indi,indj-1:indj+1),dy)

          Wt_Y(4,indi,indj)   = Wt_Y(4,indi,indj) - Gmm * Wp_Y(ii,4,indi,indj)/Wp_Y(ii,1,indi,indj) * (ux + vy)

        END DO
      END DO


      CALL Impose_BC_on_Wt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated


      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in X direction
      FWp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FWp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

    END DO
  END IF

  !*---------------------------------------------
  !*Update all subtimenodes but the first one
  !*---------------------------------------------
  DO ii = 2,MSteps !*NB: The first one does not change

    !*---------------------------------------------
    !*If last iteration, do this only for the last subtimenode
    !*---------------------------------------------
    IF ((kk == KCorr) .AND. (ii .NE. MSteps)) THEN
      CYCLE
    END IF

    !*---------------------------------------------
    !*Compute combination of the fluxes
    !*---------------------------------------------
    FLUXES_U   = 0.0
    FLUXES_W_X = 0.0
    FLUXES_W_Y = 0.0
    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)

      !*For staggering in X direction
      !*NB:+1 in X direction
      FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + thetaDeC(ii,jj)*FWp_X(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)

      !*For staggering in Y direction
      !*NB:+1 in Y direction
      FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + thetaDeC(ii,jj)*FWp_Y(jj,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

    END DO

    !*---------------------------------------------
    !*Actual update
    !*---------------------------------------------
    !*->Normal update of conserved variables
    !*->Normal update of density
    !*->Linear system for p
    !*->Modified update for v
    !*---------------------------------------------


    !*->Normal update of conserved variables
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Up(1,1:nVar,1:nElemsX,1:nElemsY) + dt*FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)

    !*->Normal update of density
    Wa_X(ii,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) = Wp_X(1,1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_W_X(1,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
    Wa_Y(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) = Wp_Y(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) + dt*FLUXES_W_Y(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)

    !*->Linear system for p

    !*----------------------
    !*IMEX linear system W_X <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEXMOMENTUM_Matrix_X(dt*betaDeC(ii),Wp_X( ii, 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs DeC
    CALL IMEXMOMENTUM_DeC_rhs_X(FLUXES_W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(1,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),Wp_X(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES X"
    its = 10000
    sol_X = initial_guess_X
    call pmgmres_ilu_cr(NRows_X, NNZsparse_X, RowStart_X, Columns_X, Values_X, sol_X, rhs_X, its, MIN(50,NRows_X-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    !*Jacobi
    ! PRINT*, "Jacobi X"
    iterations=0
    CALL jacobi_2(RowStart_X,Columns_X,Values_X,Diagonal_X,rhs_X,sol_X,iterations,initial_guess_X)

#ifdef COUNTJACOBI
    JacobiCounter_X = JacobiCounter_X +1
    JacobiIterations_X(JacobiCounter_X) = iterations
#endif

#endif

    !*Copying sol_X into Wa_X(ii,4,:,:)
    CALL Put_Vector_In_Matrix_X(sol_X(1:NRows_X),Wa_X(ii,nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))

    !*----------------------
    !*IMEX linear system W_Y <=============
    !*----------------------
    !*Assembly IMEX Matrix DeC
    CALL IMEXMOMENTUM_Matrix_Y(dt*betaDeC(ii),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)) !*NB: The difference with IMEX Euler is the multiplication by betaDeC(ii)
    !*Assembly IMEX rhs Euler
    CALL IMEXMOMENTUM_DeC_rhs_Y(FLUXES_W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),Wp_Y(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),betaDeC(ii),tStage(ii))

#if(1==0)
    PRINT*, "GMRES is slower in this case, switch to Jacobi"
    STOP
    PRINT*, "GMRES Y"
    its = 10000
    sol_Y = initial_guess_Y
    call pmgmres_ilu_cr(NRows_Y, NNZsparse_Y, RowStart_Y, Columns_Y, Values_Y, sol_Y, rhs_Y, its, MIN(50,NRows_Y-1), 1e-13, 1e-6 )
    !*Inner iterations MR must be less than  order of the system N
#else
    ! PRINT*, "Jacobi Y"
    !*Jacobi
    iterations=0
    CALL jacobi_2(RowStart_Y,Columns_Y,Values_Y,Diagonal_Y,rhs_Y,sol_Y,iterations,initial_guess_Y)

#ifdef COUNTJACOBI
    JacobiCounter_Y = JacobiCounter_Y +1
    JacobiIterations_Y(JacobiCounter_Y) = iterations
#endif

#endif

    !*Copying sol_Y into W_Y(ii,4,:,:)
    CALL Put_Vector_In_Matrix_Y(sol_Y(1:NRows_Y),Wa_Y(ii,nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))


#if(1==0)
    !*==========================
    !*IT COULD BE DELETED
    !*DIRICHLET BCs CAUSES ISSUES
    !*==========================

    ! PRINT*, "Check (Dirichlet) BC X"

    indc=(nElemsX+1)*nElemsY

    !*Later equations->Boundary conditions
    !*->Down boundary
    indj=0
    DO indi=1,nElemsX+1 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    !*->Right boundary
    indi=nElemsX+2
    DO indj=1,nElemsY
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*->Top boundary
    indj=nElemsY+1
    DO indi=1,nElemsX+1 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO



    !*->Left boundary
    indi=0
    DO indj=1,nElemsY
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_X(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_X(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_X(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    ! PRINT*, "Check (Dirichlet) BC Y"

    indc=nElemsX*(nElemsY+1)


    !*Later equations->Boundary conditions
    !*->Down boundary
    indj=0
    DO indi=1,nElemsX 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO


    !*->Right boundary
    indi=nElemsX+1
    DO indj=1,nElemsY+1
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*->Top boundary
    indj=nElemsY+2
    DO indi=1,nElemsX 
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO



    !*->Left boundary
    indi=0
    DO indj=1,nElemsY+1
      indc=indc+1

      CALL ExactFunction(&
        InitialCondition,tStage(ii),MeshGP_Y(:,indi,indj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)

      IF (( Wtemp(nVar) .NE. rhs_Y(indc) ) .OR. ( Wtemp(nVar) .NE. Wa_Y(ii,nVar,indi,indj) )) THEN
        PRINT*, "Problem", indi, indj
        STOP
      END IF

    END DO

    !*==========================
    !*IT COULD BE DELETED
    !*DIRICHLET BCs CAUSES ISSUES
    !*==========================
#endif



    !*->Modified update for rov
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY
      DO indi=1,nElemsX+1 !*NB:+1 in X direction
        Wa_X(ii,2:3,indi,indj) = Wp_X(1,2:3,indi,indj) + dt*FLUXES_W_X(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi-1:indi+1,indj)-Wp_X(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_X(ii,nVar,indi,indj-1:indj+1)-Wp_X(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)*py
#else
        Wa_X(ii,2,indi,indj)=Wa_X(ii,2,indi,indj)-(dt          )*betaDeC(ii)*px
        Wa_X(ii,3,indi,indj)=Wa_X(ii,3,indi,indj)-(dt          )*betaDeC(ii)*py
#endif
      END DO
    END DO

    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY+1 !*NB:+1 in Y direction
      DO indi=1,nElemsX
        Wa_Y(ii,2:3,indi,indj) = Wp_Y(1,2:3,indi,indj) + dt*FLUXES_W_Y(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi-1:indi+1,indj)-Wp_Y(ii,nVar,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_Y(ii,nVar,indi,indj-1:indj+1)-Wp_Y(ii,nVar,indi,indj-1:indj+1),dy)
#ifdef RELAXATION
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt/EPS_LM**2)*betaDeC(ii)*py
#else
        Wa_Y(ii,2,indi,indj)=Wa_Y(ii,2,indi,indj)-(dt          )*betaDeC(ii)*px
        Wa_Y(ii,3,indi,indj)=Wa_Y(ii,3,indi,indj)-(dt          )*betaDeC(ii)*py
#endif
      END DO
    END DO


  END DO


END DO


!*---------------------------------------------
!*Final output
!*---------------------------------------------
!*Actually here it is not important to have BCs
!*---------------------------------------------
U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
!*For staggering in X direction
!*NB:+1 in X direction
W_X(1:nVar,1:nElemsX+1,1:nElemsY) = Wa_X(MSteps,1:nVar,1:nElemsX+1,1:nElemsY)

!*For staggering in Y direction
!*NB:+1 in Y direction
W_Y(1:nVar,1:nElemsX,1:nElemsY+1) = Wa_Y(MSteps,1:nVar,1:nElemsX,1:nElemsY+1)


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByIMEXMOMENTUMDeC2
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE IMEXMOMENTUM_DeC_rhs_X(Wt,W0,Wm,betam,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), INTENT(IN) :: Wt, W0, WM
REAL                                                                           , INTENT(IN) :: betam
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: supp(-1:1)
REAL                  :: dx_part, dy_part, delta_t
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)
delta_t=dt*betam !*For normalization !*IMPORTANT: Normalization w.r.t. delta_t, not dt

!*RMK: +1 in X direction

rhs_X=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)*EPS_LM)**2
#else
    rhs_X(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_X(indc)=W0(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-(dt**2)*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_X(indc)=rhs_X(indc)-(dt**2)*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=rhs_X(indc)-(dt**2)*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_X(indc)=rhs_X(indc)-(dt**2)*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

    !*We have now two extra contributions: div v and lap p

    !*div v part
    !*NB: This is the divergence of the momentum actually due to the fact that I run with momentum in primitive variables
    dx_part=First_Derivative_Central_Order2(W0(2,ii-1:ii+1,jj)-Wm(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(W0(3,ii,jj-1:jj+1)-Wm(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-dt*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2)
#else
    rhs_X(indc)=rhs_X(indc)-dt*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2)
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_X(indc)=rhs_X(indc)-dt*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*dx**2+dy_part*dx**2)
#else
    rhs_X(indc)=rhs_X(indc)-dt*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part+dy_part)
#endif

    !*lap p part
    lapp=Laplacian_Central_Order2(Wm(4,ii-1:ii+1,jj-1:jj+1),dx,dy)

#ifdef NORMALIZEALL
    rhs_X(indc)=rhs_X(indc)-(dt       )**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(lapp*(dx/delta_t)**2)
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-(dt/EPS_LM)**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(lapp*dx**2)
#else
    rhs_X(indc)=rhs_X(indc)-(dt       )**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(lapp*dx**2)
#endif    
#else
#ifdef RELAXATION
    rhs_X(indc)=rhs_X(indc)-(dt/EPS_LM)**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*lapp
#else
    rhs_X(indc)=rhs_X(indc)-(dt       )**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*lapp
#endif    
#endif    

  END DO !ii
END DO !jj


!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_X(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO

!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_X(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_X(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_X(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_X(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_X(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_X(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEXMOMENTUM_DeC_rhs_X
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE IMEXMOMENTUM_DeC_rhs_Y(Wt,W0,Wm,betam,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition

USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: First_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Second_Mixed_Derivative_Central_Order2
USE MOD_Reconstruction,     ONLY: Laplacian_Central_Order2

USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction

#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), INTENT(IN) :: Wt, W0, Wm
REAL                                                                           , INTENT(IN) :: betam
REAL                                                                           , INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER               :: ii, jj
INTEGER               :: indc
REAL                  :: dx,dy
REAL                  :: u, v, ux, vx, uy, vy, uxx, vxx, uyy, vyy, uxy, vxy
REAL                  :: p, px, py, lapp
REAL                  :: et, etx, ety !*Extra term=et=(ro*-ro)/ro/ro* divided by (in case relaxation is active) eps_lm**2
REAL                  :: supp(-1:1)
REAL                  :: div_v_grad_v
REAL                  :: ro_p_contribution
REAL                  :: dx_part, dy_part, delta_t
REAL, DIMENSION(nVar) :: Utemp, Wtemp
!-------------------------------------------------------------------------------!

dx=MESH_DX(1)
dy=MESH_DX(2)
delta_t=dt*betam !*For normalization !*IMPORTANT: Normalization w.r.t. delta_t, not dt

!*RMK: +1 in Y direction

rhs_Y=0.0

!*Let us start by the known contribution inside the domain

!*Inside the domain
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)*EPS_LM)**2
#else
    rhs_Y(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*((dx/delta_t)       )**2
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=(W0(4,ii,jj) + dt*Wt(4,ii,jj))*dx**2
#else
    rhs_Y(indc)=W0(4,ii,jj) + dt*Wt(4,ii,jj)
#endif

  END DO !ii
END DO !jj


!*Now let us add the extra contribution to be computed
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1


    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*Alex's suggestion with upwinded (already computed operators)
    dx_part=First_Derivative_Central_Order2(Wt(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(Wt(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-(dt**2)*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_Y(indc)=rhs_Y(indc)-(dt**2)*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2) !*NB: It should be - because in Wt I already have "-"
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=rhs_Y(indc)-(dt**2)*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*dx**2+dy_part*dx**2) !*NB: It should be - because in Wt I already have "-"
#else
    rhs_Y(indc)=rhs_Y(indc)-(dt**2)*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part+dy_part) !*NB: It should be - because in Wt I already have "-"
#endif

    !*We have now two extra contributions: div v and lap p

    !*div v part
    !*NB: This is the divergence of the momentum actually due to the fact that I run with momentum in primitive variables
    dx_part=First_Derivative_Central_Order2(W0(2,ii-1:ii+1,jj)-Wm(2,ii-1:ii+1,jj),dx)
    dy_part=First_Derivative_Central_Order2(W0(3,ii,jj-1:jj+1)-Wm(3,ii,jj-1:jj+1),dy)

#ifdef NORMALIZEALL
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-dt*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*((dx/delta_t)*EPS_LM)**2+dy_part*((dx/delta_t)*EPS_LM)**2)
#else
    rhs_Y(indc)=rhs_Y(indc)-dt*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*((dx/delta_t)       )**2+dy_part*((dx/delta_t)       )**2)
#endif
#elif defined(NORMALIZELINEARSYSTEMS)
    rhs_Y(indc)=rhs_Y(indc)-dt*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part*dx**2+dy_part*dx**2)
#else
    rhs_Y(indc)=rhs_Y(indc)-dt*betam*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(dx_part+dy_part)
#endif

    !*lap p part
    lapp=Laplacian_Central_Order2(Wm(4,ii-1:ii+1,jj-1:jj+1),dx,dy)

#ifdef NORMALIZEALL
    rhs_Y(indc)=rhs_Y(indc)-(dt       )**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(lapp*(dx/delta_t)**2)
#elif defined(NORMALIZELINEARSYSTEMS)
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-(dt/EPS_LM)**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(lapp*dx**2)
#else
    rhs_Y(indc)=rhs_Y(indc)-(dt       )**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*(lapp*dx**2)
#endif    
#else
#ifdef RELAXATION
    rhs_Y(indc)=rhs_Y(indc)-(dt/EPS_LM)**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*lapp
#else
    rhs_Y(indc)=rhs_Y(indc)-(dt       )**2*betam**2*Gmm*Wm(4,ii,jj)/Wm(1,ii,jj)*lapp
#endif    
#endif

  END DO !ii
END DO !jj





!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState1(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_Y(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition B not implemented"
      STOP
  END SELECT
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState2(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_Y(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition R not implemented"
      STOP
  END SELECT
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState3(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_Y(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition T not implemented"
      STOP
  END SELECT
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1

  SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
    CASE(1,2,5) ! Periodic,Transmissive,Reflecting u=u
      !*Do nothing, keep 0
    CASE(3,4) ! Inflow,Outflow u=something
      !*rhs_Y(indc)=PrimRefState4(nVar)
      CALL ExactFunction(&
        InitialCondition,t,MeshGP_Y(:,ii,jj,1,1),Utemp(1:nVar))
      CALL ConsToPrim(Utemp(1:nVar),Wtemp(1:nVar))
      rhs_Y(indc)=Wtemp(nVar)
      ! PRINT*, "TO BE DELETED", rhs_Y(indc)
    CASE DEFAULT
      PRINT*, "Boundary condition L not implemented"
      STOP
  END SELECT
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE IMEXMOMENTUM_DeC_rhs_Y
#endif
#endif
!===============================================================================!
END MODULE MOD_TimeDiscretization
!-------------------------------------------------------------------------------!
